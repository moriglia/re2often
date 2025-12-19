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
    use flap ! CLI parser: command_line_interface
    use lincoa_mod
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
    logical :: editOutProbs
    logical :: optimize_s
    integer, allocatable :: encodingVector(:)
    double precision, allocatable :: prob_out(:)
    double precision, allocatable :: thresholds(:)
    double precision :: s(1)
    integer :: M_half
    double precision :: I_neg

    type(command_line_interface) :: cli
    integer :: error

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


    ! +------------+
    ! | CLI parser |
    ! +------------+
    call cli%init(&
        progname = "GMI", &
        version  = "0", &
        authors  = "Marco Origlia", &
        license  = "GPL-3.0-or-later")

    call cli%add(&
        switch='--nsnr', &
        help='NSNR value', &
        required=.true., &
        act='store', &
        error=error)
    call cli%add(switch='--snr', &
        help='SNR range (2)', &
        required=.true., &
        act='store', &
        nargs='2', &
        error=error)
    call cli%add(switch='--bps', &
        help='Bits per symbol', required=.false., &
        def='2', &
        act='store', &
        error=error)
    call cli%add(switch='--outdir', &
        help='Output directory', &
        required=.true., &
        act='store', &
        error=error)
    call cli%add(switch='--hard', &
        help='Hard reconciliation (implies "-r")', &
        required=.false., &
        act='store_true', &
        def='.false.', &
        error=error)
    call cli%add(switch='-r', &
        help='Reverse Reconciliation', &
        required=.false., &
        act='store_true', &
        def='.false.', &
        error=error)
    call cli%add(switch='-u', &
        help='Uniform output probability thresholds', &
        required=.false., &
        act='store_true', &
        def='.false.', &
        error=error)
    call cli%add(switch='-c', &
        help='Edit config and monoConfig value', &
        required=.false., &
        act='store', &
        def='0', &
        error=error)
    call cli%add(switch='--natural', &
        help='Use natural encoding', &
        required=.false., &
        act='store_true', &
        def='.false.', &
        error=error)
    call cli%add(switch='--ml', &
        help='Use maximum likelyhood', &
        required=.false., &
        act='store_true', &
        def='.false.', &
        error=error)
    call cli%add(switch='--encoding', &
        help='Specify the labelling, a permutation of numbers (0,...,2^BPS-1) is expected', &
        required=.false., &
        act='store', &
        nargs='+', &
        def='0', & ! Does not mean anything, just to prevent the cli to complain
        error=error)
    call cli%add(switch='--prob-out', &
         help="Set output probabilities (specify the probabilities of the top M/2-1 but last symbols)",&
         required=.false.,&
         nargs='+',&
         act="store",&
         def='0.25',&
         error=error)
    call cli%add(switch='--opt-s', &
         help="Optimise parameter S",&
         required=.false.,&
         nargs='+',&
         act="store_true",&
         def='.false.',&
         error=error)


    call cli%parse(error=error)
    if (error /= 0) stop

    ! Retrieve values
    call cli%get(switch='--nsnr', val=nsnr)
    call cli%get(switch='--snr', val=snr)
    call cli%get(switch='--bps', val=bps)
    call cli%get(switch='--outdir', val=output_root)
    call cli%get(switch='--hard', val=isHard)
    call cli%get(switch='-r', val=isReverse)
    call cli%get(switch='-u', val=uniform_th)
    call cli%get(switch='--natural', val=encodingNatural)
    call cli%get(switch='--ml', val=useML)
    call cli%get(switch='--opt-s', val=optimize_s)

    if (cli%is_passed(switch='--encoding')) then
        call cli%get_varying(switch='--encoding', val=encodingVector)
    end if
    editOutProbs = cli%is_passed(switch='--prob-out')
    if (editOutProbs) then
        call cli%get_varying(switch='--prob-out', val=prob_out)
    end if

    editConfig = cli%is_passed(switch="-c")

    if (editConfig) then
        call cli%get(switch="-c", val=monoConfig)
    else
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
        outdata(:,2:3) = 0
        outdata(:, 4) = 1d0
    end if
    ! snrdb        => outdata(:, 1)[1]
    ! snrdb_scaled => outdata(:, 2)[1]
    ! I            => outdata(:, 3)[1]

    nm = noisemapper_create(bps)
    M_half = ishft(nm%M, -1)
    if (encodingNatural) then
        call noisemapper_set_encoding_natural(nm)
        if (me == 1 .and. cli%is_passed(switch='--encoding')) then
            print *, "Warning: --natural overrides --encoding (ignored)"
        end if
    elseif( cli%is_passed(switch='--encoding')) then
        call noisemapper_set_encoding_custom(nm, encodingVector)
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
    if (editOutProbs) then
        allocate(thresholds(nm%M-1))
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

        if (editOutProbs) then
            thresholds(M_half) = 0d0
            do ii = 1,M_half-1
                thresholds(M_half+ii) = noisemapper_inverse_Fy_search(nm, 0.5d0 + sum(prob_out(:ii)))
                thresholds(M_half-ii) = -thresholds(M_half + ii)
            end do
            call noisemapper_set_y_thresholds(nm, thresholds)
        else if (uniform_th) then
            call noisemapper_set_y_thresholds_uniform(nm)
        else
            call noisemapper_set_y_thresholds(nm)
        end if

        if (isReverse) then
            if (isHard) then
                call noisemapper_update_hard_reverse_tables(nm)
                if (useML) then
                    outdata(i_snr, 3)[1] = I_s_ml_hard_direct(q_ml_hard_direct_prod)
                else
                    outdata(i_snr, 3)[1] = I_s_map_hard_reverse(q_map_hard_product)
                end if
            else
                call noisemapper_set_Fy_grids(nm)
                if (useML) then
                    if (optimize_s) then
                        s(1) = 1d0
                        call lincoa(gmi_ml_reverse_soft_s, s, f=I_neg, xl=[1d-9])
                        outdata(i_snr, 3)[1] = -I_neg
                        outdata(i_snr, 4)[1] = s(1)
                    else
                        outdata(i_snr, 3)[1] = I_s_ml_soft_reverse(nm)
                    end if
                else
                    if (optimize_s) then
                        s(1) = 1d0
                        call lincoa(gmi_map_reverse_soft_s, s, f=I_neg, xl=[1d-9])
                        outdata(i_snr, 3)[1] = -I_neg
                        outdata(i_snr, 4)[1] = s(1)
                    else
                        outdata(i_snr, 3)[1] = I_s_map_soft_reverse(q_map_soft_reverse_prod)
                    end if
                end if
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
        if (editOutProbs) then
            output_root = trim(output_root)//"/prob_out_custom/"
            write(output_root, '(A, "p", f6.4)') trim(output_root), prob_out(1)
            do ii = 2, M_half-1
                write(output_root, '(A, "_p", f6.4)') trim(output_root), prob_out(ii)
            end do
        else if (uniform_th) then
            output_root = trim(output_root)//"/uniform"
        else
            output_root = trim(output_root)//"/equidis"
        end if
        if (optimize_s) then
            output_root = trim(output_root)//"/opt_s"
        end if

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

contains

    subroutine not_implemented
        print *, "Not implemented"; stop
    end subroutine not_implemented

    subroutine gmi_map_reverse_soft_s(s, I_neg)
        double precision, intent(in)  :: s(:)
        double precision, intent(out) :: I_neg

        I_neg = -I_s_map_soft_reverse(q_map_soft_reverse_prod, s(1))
    end subroutine gmi_map_reverse_soft_s

    subroutine gmi_ml_reverse_soft_s(s, I_neg)
        double precision, intent(in)  :: s(:)
        double precision, intent(out) :: I_neg

        I_neg = -I_s_ml_soft_reverse(nm, s(1))
    end subroutine gmi_ml_reverse_soft_s

end program gmi
