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
program gmi
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Numerically evaluate the Generalized Mutual Information
    use iso_fortran_env, only: lock_type, dp => real64
    use io_fortran_lib, only: from_file, to_file
    ! use stdlib_random, only: stdlib_random_seed => random_seed
    ! use stdlib_stats_distribution_normal, only: rvs_normal
    ! use re2often_noise_mapper, only: TNoiseMapper
    use re2often_noisemapper
    use re2often_utils, only: save_data, make_directory_and_file_name
    use forbear, only: bar_object
    use re2often_mi ! defines a noisemapper_type object
    implicit none

    ! +---------------------+
    ! | Command line inputs |
    ! +---------------------+
    integer :: argc
    character(len=500), allocatable :: argv(:)

    character(250) :: output_root
    character(250) :: output_dir
    character(250) :: output_name
    character(5)   :: header(4)
    integer        :: io
    double precision :: snr(2)         ! Signal to Noise Ratio in dB
    integer :: nsnr           ! Number of SNR points
    integer :: bps            ! Bits per symbol
    logical :: isReverse      ! Whether to calculate the M.I. for the reverse reconciliation
    logical :: isHard         ! Whether to perform hard reverse reconciliation
    logical :: uniform_th     ! Whether to set thresholds for uniform output symbols
    logical :: editConfig     ! Whether to perform the calculation for a specific monotonicity configuration
    integer :: monoConfig     ! Selected monotonicity configuration
    logical :: encodingNatural
    logical :: useML

    ! +-------------+
    ! | Output data |
    ! +-------------+
    double precision, allocatable, target :: outdata(:,:)[:]
    logical, allocatable                  :: snr_done(:)[:]

    ! +---------------------+
    ! | Image-specific data |
    ! +---------------------+
    integer :: me, n_im

    type(bar_object) :: progress_bar

    type(lock_type) :: lck[*]


    ! generic iteration variables
    integer :: i_snr, ii

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
    nsnr = 11
    snr  = [3.5d0, 4d0]
    bps  = 2
    isReverse = .false.
    isHard = .false.
    uniform_th = .false.
    editConfig = .false.
    encodingNatural = .false.
    useML = .false.

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
        elseif (argv(ii) == "--bps") then
            read(argv(ii+1),*) bps
            ii = ii + 2
            ! print *, "bps", bps
        elseif (argv(ii) == "--outdir") then
            call get_command_argument(ii+1, output_root)
            ii = ii + 2
        elseif (argv(ii) == "-h") then
            isHard = .true.
            ii = ii + 1
        elseif (argv(ii) == "-r") then
            isReverse = .true.
            ii = ii + 1
        elseif (argv(ii) == "-u") then
            uniform_th = .true.
            ii = ii + 1
        elseif(argv(ii) == "-c") then
            editConfig = .true.
            read(argv(ii + 1), *) monoConfig
            ii = ii + 2
        elseif (argv(ii) == "--natural") then
            encodingNatural = .true.
            ii = ii + 1
        elseif (argv(ii) == "--ml") then
            useML = .true.
            ii = ii + 1
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


    me   = this_image()
    n_im = num_images()

    ! +----------------------------------------+
    ! | Allocation of simulation result arrays |
    ! +----------------------------------------+
    allocate(snr_done(nsnr)[*])
    allocate(outdata(nsnr, 4)[*])
    if (me==1) then
        snr_done(:) = .false.
        outdata(:, 1) = [(snr(1) + real(ii, dp)*(snr(2)-snr(1))/real(nsnr - 1, dp), ii = 0, nsnr-1)]
        outdata(:,2:) = 0
    end if
    ! snrdb        => outdata(:, 1)[1]
    ! snrdb_scaled => outdata(:, 2)[1]
    ! I            => outdata(:, 3)[1]

    nm = noisemapper_create(bps)
    if (encodingNatural) then
        call noisemapper_set_encoding_natural(nm)
    end if
    if (isReverse .and. (.not. isHard)) then
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

        call noisemapper_update_N0_from_snrdb(nm, outdata(i_snr, 1)[1])
        sqrtN0 = sqrt(nm%N0)

        if (isReverse) then
            if (isHard) then
                if (uniform_th) then
                    call noisemapper_set_y_thresholds_uniform(nm)
                    ! outdata(i_snr, 3)[1] = I_hard_reverse_uniform_output_th(outdata(i_snr, 1)[1])
                else
                    call noisemapper_set_y_thresholds(nm)
                    ! outdata(i_snr, 3)[1] = I_hard_reverse_equidistant_th(outdata(i_snr, 1)[1])
                end if
                call noisemapper_update_hard_reverse_tables(nm)
                if (useML) then
                    outdata(i_snr, 3)[1] = I_s_ml_hard_direct(q_ml_hard_direct_prod)
                else
                    outdata(i_snr, 3)[1] = I_s_map_hard_reverse(q_map_hard_product)
                end if
            else
                if (me == 1) then
                    print *, "Not implemented"
                end if
                stop
            end if
        else
            if (useML) then
                if (me == 1) then
                    print *, "Not Implemented"
                end if
                stop
            end if
            outdata(i_snr, 3)[1] = I_s_map_soft_direct(q_map_soft_direct_prod)
        end if
    end do loop_snr

    ! +-----------------------------------------+
    ! | First Image saves the simulation result |
    ! +-----------------------------------------+
    sync all

    if (me == 1) then
        outdata(:, 2) = outdata(:, 1) - 10*log10(outdata(:,3))
        if (uniform_th) then
            output_root = trim(output_root)//"/uniform"
        else
            output_root = trim(output_root)//"/equidis"
        end if
        outdata(:, 4) = 1d0 ! Optimization not implemented

        call make_directory_and_file_name(output_root, bps, isReverse, isHard, &
            snr, nsnr, 0, 0, 0, 0, output_dir, output_name)
        call execute_command_line("mkdir -p " // trim(output_dir))

        open(newunit=io, file=trim(output_dir) // "/" // trim(output_name) // ".log", &
            status="replace", action="write")

        write(io, '(A, T16, A, T32, A, T48, A)') &
            "SNR [dB]", "Eb/N0 [dB]", "I", "s_opt"
        do i_snr = 1, nsnr
            write(io, '(f12.8, T16, f12.9, T32, E12.3E3, T48, E12.3E3)') &
                outdata(i_snr, :)
        end do
        close(io)


        header(1) = trim("SNR")
        header(2) = trim("E_b/N_0")
        header(3) = trim("I")
        header(4) = trim("s_opt")
        call to_file(x=outdata, file=trim(output_dir)//"/"//trim(output_name)//".csv", &
            header=header, fmt="e")
    end if
end program gmi
