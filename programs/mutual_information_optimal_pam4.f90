! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025  Marco Origlia

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
program mutual_information_optimal_pam4
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Numerically evaluate the Mutual Information of various schemes
    use iso_fortran_env, only: lock_type, dp => real64
    use io_fortran_lib, only: from_file, to_file
    use re2often_noisemapper
    use re2often_utils, only: save_data, make_directory_and_file_name
    use forbear, only: bar_object
    use re2often_mi ! defines a noisemapper_type object
    use quadpack, only: dqags
    implicit none

    ! +---------------------+
    ! | Command line inputs |
    ! +---------------------+
    integer :: argc
    character(len=500), allocatable :: argv(:)

    character(len=250) :: output_root
    character(len=250) :: output_dir
    character(len=250) :: output_name
    character(len=5)   :: header(5)
    integer        :: io
    double precision :: snr(2)         ! Signal to Noise Ratio in dB
    integer :: nsnr           ! Number of SNR points
    integer :: bps            ! Bits per symbol
    logical :: isHard         ! Whether to perform hard reverse reconciliation
    logical :: editConfig     ! Whether to perform the calculation for a specific monotonicity configuration
    integer :: monoConfig     ! Selected monotonicity configuration

    ! +-------------+
    ! | Output data |
    ! +-------------+
    double precision, allocatable, target :: outdata(:,:)[:]

    ! +-------------------+
    ! | Search parameters |
    ! +-------------------+
    double precision :: p_tilde, t_tilde, I_tmp
    integer :: p_i ! index of p_tilde in the grid
    integer :: t_i ! index of t_tilde in the grid
    integer :: p_num, t_num
    double precision :: pres, tres


    ! +-----------------------------------------+
    ! | Numerical integration auxiliary buffers |
    ! +-----------------------------------------+
    real(c_double) :: Abserr
    integer :: Neval, Ier, Limit, Lenw, Last
    integer, allocatable :: Iwork(:)
    real(c_double), allocatable :: Work(:)



    logical, allocatable                  :: snr_done(:)[:]

    ! +---------------------+
    ! | Image-specific data |
    ! +---------------------+
    integer :: me, n_im

    type(bar_object) :: progress_bar

    type(lock_type) :: lck[*]


    ! generic iteration variables
    integer :: i_snr, ii, jj

    if (this_image() == 1) then
        print *, " +----------------------------+"
        print *, " | Mutual Information program |"
        print *, " +----------------------------+"
    end if

    argc = command_argument_count()
    allocate(argv(argc))

    do ii = 1, argc
        call get_command_argument(ii, argv(ii))
    end do


    ! +--------------------+
    ! | Set default values |
    ! +--------------------+
    nsnr = 251
    snr  = [-50, 0]
    bps  = 2
    isHard = .false.
    editConfig = .false.
    pres = 1d-3
    tres = 5d-3
    p_num = floor((0.25d0-pres)/pres, kind(1)) + 1

    ii = 1
    do while(ii <= argc)
        if (argv(ii) == "--nsnr") then
            read(argv(ii+1),*) nsnr
            ii = ii + 2
            ! print *, "nsnr", nsnr
        elseif (argv(ii) == "--snr") then
            read(argv(ii+1), *) snr(1)
            read(argv(ii+2), *) snr(2)
            ii = ii+3
            ! print *, "snr", snr
        ! elseif (argv(ii) == "--bps") then
        !     read(argv(ii+1),*) bps
        !     ii = ii + 2
        !     ! print *, "bps", bps
        elseif (argv(ii) == "--outdir") then
            call get_command_argument(ii+1, output_root)
            ii = ii + 2
        elseif (argv(ii) == "-h") then
            isHard = .true.
            ii = ii + 1
        elseif(argv(ii) == "-c") then
            editConfig = .true.
            read(argv(ii + 1), *) monoConfig
            ii = ii + 2
        else
            print *, "Unrecognized argument: ", argv(ii)
            stop
        end if
    end do

    if (.not. editConfig) then
        ! set default configuration
        monoConfig = 2
        do ii = 0, bps - 1
            monoConfig = monoConfig * (1 + ishft(1, (ishft(1, ii)))) ! useless but nice to know :D
        end do
    end if

    if (.not. isHard) then
        Limit = 100
        Lenw  = 400
        allocate(Iwork(Limit))
        allocate(Work(Lenw))
    end if


    me   = this_image()
    n_im = num_images()

    ! +----------------------------------------+
    ! | Allocation of simulation result arrays |
    ! +----------------------------------------+
    allocate(snr_done(nsnr)[*])
    allocate(outdata(nsnr, 7)[*]) ! SNR [dB], Eb/N0 [dB], I [b/cu], p_opt, t_opt
    if (me==1) then
        snr_done(:) = .false.
        outdata(:, 1) = [(snr(1) + real(ii, dp)*(snr(2)-snr(1))/real(nsnr - 1, dp), ii = 0, nsnr-1)]
        outdata(:, 3) = -1 ! Impossible value for the mutual information
    end if
    ! snrdb        => outdata(:, 1)[1]
    ! snrdb_scaled => outdata(:, 2)[1]
    ! I            => outdata(:, 3)[1]

    nm = noisemapper_create(bps)
    if (.not. isHard) then
        call noisemapper_set_monotonicity(nm)
        ! Allocates the monotonicity configuration and sets the default
        if (editConfig) then
            i_snr = monoConfig ! i_snr is opportunistically used as a temporary variable
            do ii = 0, nm%M-1
                if (iand(i_snr, 1) == 0) then
                    nm%monotonicity_configuration(ii) = .false.
                else
                    nm%monotonicity_configuration(ii) = .true.
                end if
                i_snr = ishft(i_snr, -1)
            end do
        end if
    end if

    if (me == 1) then
        call progress_bar%initialize(&
            filled_char_string='+', prefix_string='SNR points progress |',&
            suffix_string='| ', add_progress_percent=.true.)
        call progress_bar%start
    end if

    sync all

    i_snr = 1
    loop_snr : do while (i_snr .le. nsnr)
        lock(lck[1])
        do while(snr_done(i_snr)[1])
            i_snr = i_snr + 1
            if (i_snr .gt. nsnr) then
                unlock(lck[1])
                exit loop_snr
            end if
        end do
        snr_done(i_snr)[1] = .true.
        unlock(lck[1])
        if (me == 1) then
            call progress_bar%update(current=real(i_snr, 8)/real(nsnr, 8))
        end if

        loop_p_tilde: do p_i = 0, p_num-1
            p_tilde = 0.25d0 + p_i*pres
            call noisemapper_set_symbol_probabilities(nm, &
                [.5d0-p_tilde, p_tilde, p_tilde, .5d0-p_tilde])
            call noisemapper_update_N0_from_snrdb(nm, outdata(i_snr, 1)[1])

            t_num = max(ceiling(nm%N0/tres), 1)

            loop_t_tilde: do t_i = 0, t_num-1
                t_tilde = 2d0 + t_i*tres
                call noisemapper_set_y_thresholds(nm, [-t_tilde, 0d0, t_tilde])

                if (isHard) then
                    call noisemapper_update_hard_reverse_tables(nm)
                    I_tmp = 0

                    do jj = 0, nm%M-1
                        do ii = 0, nm%M-1
                            I_tmp = I_tmp + nm%fwd_probabilities(ii, jj) * nm%probabilities(ii) * &
                                (log0(nm%fwd_probabilities(ii, jj)) - log0(nm%delta_Fy(jj)))
                        end do
                    end do
                    I_tmp = I_tmp/log(2d0)
                else
                    call noisemapper_set_Fy_grids(nm)
                    call dqags(f_soft_reverse, 0d0, 1d0, 1d-12, 1d-6, &
                        I_tmp, Abserr, Neval, Ier, &
                        Limit, Lenw, Last, Iwork, Work)

                    if (Ier /= 0) then
                        print '("Error at ", f10.3, " [dB]: error ", i1)', &
                            outdata(i_snr, 1)[1], Ier
                    end if

                    I_tmp = I_tmp + H_Xhat(nm)
                end if

                if (I_tmp .gt. outdata(i_snr, 3)[1]) then
                    outdata(i_snr, 3)[1] = I_tmp
                    outdata(i_snr, 4)[1] = p_tilde
                    outdata(i_snr, 5)[1] = t_tilde
                    outdata(i_snr, 6)[1] = p_num
                    outdata(i_snr, 7)[1] = t_num
                end if
            end do loop_t_tilde
        end do loop_p_tilde
    end do loop_snr

    ! +-----------------------------------------+
    ! | First Image saves the simulation result |
    ! +-----------------------------------------+
    sync all

    if (me == 1) then
        outdata(:, 2) = outdata(:, 1) - 10*log10(outdata(:,3))

        call make_directory_and_file_name(output_root, bps, .true., isHard, &
            snr, nsnr, 0, 0, 0, 0, output_dir, output_name)
        call execute_command_line("mkdir -p " // trim(output_dir))

        open(newunit=io, file=trim(output_dir) // "/" // trim(output_name) // ".log", &
            status="replace", action="write")

        write(io, '(A, T16, A, T32, A, T48, A, T64, A)') &
            "SNR [dB]", "Eb/N0 [dB]", "I [b/cu]", "p_opt", "t_opt"
        do i_snr = 1, nsnr
            write(io, '(f12.8, T16, f12.9, T32, E12.3E3, T48, E12.3E3, T64, E12.3E3, T80, E12.3E3, T96, E12.3E3)') &
                outdata(i_snr, :)
        end do
        close(io)


        header(1) = trim("SNR")
        header(2) = trim("scSNR")
        header(3) = trim("I")
        header(4) = trim("p_opt")
        header(5) = trim("t_opt")
        call to_file(x=outdata(:,:5), file=trim(output_dir)//"/"//trim(output_name)//".csv", &
            header=header, fmt="e")
    end if
end program mutual_information_optimal_pam4
