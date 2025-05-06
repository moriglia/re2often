! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2024  Marco Origlia

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
program reverse_reconciliation
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Perform reverse reconciliation Montecarlo Simulation
    use iso_fortran_env, only: dp => real64
    use io_fortran_lib, only: from_file, to_file
    use stdlib_random, only: stdlib_random_seed => random_seed
    use stdlib_stats_distribution_normal, only: rvs_normal
    use stdlib_stats_distribution_uniform, only: rvs_uniform
    ! use re2often_noise_mapper, only: TNoiseMapper
    use re2often_noisemapper
    use re2often_utils, only: save_data, make_directory_and_file_name
    use ldpc_decoder, only: TDecoder
    use forbear, only: bar_object
    implicit none

    ! +---------------------+
    ! | Command line inputs |
    ! +---------------------+
    integer :: argc
    character(len=500), allocatable :: argv(:)

    character(500) :: tanner_file
    character(250) :: output_root
    character(250) :: output_dir
    character(250) :: output_name
    integer        :: io
    double precision :: snr(2)         ! Signal to Noise Ratio in dB
    integer :: nsnr           ! Number of SNR points
    integer :: bps            ! Bits per symbol
    integer :: min_ferr       ! Minimum frame errors
    integer :: max_sim        ! Maximum simulation loops
    integer :: min_sim        ! Minimum simulation loops
    integer :: max_iter       ! Maximum number of LDPC iterations
    logical :: isHard         ! Whether to perform hard reverse reconciliation
    logical :: uniform_th     ! Whether to set thresholds for uniform output symbols
    logical :: tanner_header  ! Whether tanner file has a header
    logical :: onlyinfo       ! Whether to compare only the first N-M bits, instead of whole frame
    logical :: doEve          ! Eve is performing error correction
    logical :: useInterleaver ! Scramble bits

    integer, allocatable :: edge_definition(:,:)

    ! +-------------+
    ! | Output data |
    ! +-------------+
    double precision, allocatable, target :: outdata(:,:)
    double precision, pointer             :: snrdb(:)
    double precision, pointer             :: ber(:)
    double precision, pointer             :: fer(:)

    ! +---------------------+
    ! | Image-specific data |
    ! +---------------------+
    integer :: me, n_im
    integer :: seed(8), time_seed

    ! +----------------------------------+
    ! | Simulation buffers and variables |
    ! +----------------------------------+
    integer, allocatable :: b_err(:)[:] ! Cumulative bit error count per SNR value
    integer, allocatable :: f_err(:)[:] ! Cumulative frame error count per SNR value
    integer, allocatable :: f_cnt(:)[:] ! Total frame count per SNR value

    integer, allocatable :: x_i(:) ! Indexes of generated symbols
    ! double precision, allocatable :: x(:), y(:), lappr(:), lappr_out(:)
    double precision, allocatable :: y(:), lappr(:), lappr_out(:)
    ! generated symbol, Gaussian channel output, LAPPRs before and after decoding
    double precision, allocatable :: nhat(:) ! Soft Metric
    integer, allocatable :: xhat(:) ! Decided symbol
    logical, allocatable :: word(:), synd(:) ! word received by Bob and it's syndrome
    integer :: new_errors ! Number of errors of the current iteration
    integer :: K          ! Number of information bit per frame
    integer :: N_iter     ! Number of iterations at the end of the decoding

    double precision :: sigma
    double precision :: alpha

    type(TDecoder)     :: decoder
    type(noisemapper_type) :: nm

    type(bar_object) :: progress_bar


    ! generic iteration variables
    integer :: i, i_snr, i_frame

    if (this_image() == 1) then
        print *, " +-----------------------------------+"
        print *, " | REVERSE reconciliation simulation |"
        print *, " +-----------------------------------+"
    end if

    argc = command_argument_count()
    allocate(argv(argc))

    do i = 1, argc
        call get_command_argument(i, argv(i))
    end do


    ! +--------------------+
    ! | Set default values |
    ! +--------------------+
    nsnr = 11
    snr  = [3.5d0, 4d0]
    bps  = 2
    min_ferr = 50
    max_sim  = 10000
    min_sim  = 250
    max_iter = 50
    tanner_file = "assets/codes/dvbs2ldpc0.500.csv"
    output_root  = "res/rate1d2"
    isHard = .false.
    uniform_th = .false.
    alpha = 1d0
    onlyinfo = .false.
    tanner_header = .false.
    doEve = .false.
    useInterleaver = .false.

    i = 1
    do while(i <= argc)
        if (argv(i) == "--nsnr") then
            read(argv(i+1),*) nsnr
            i = i + 2
            ! print *, "nsnr", nsnr
        elseif (argv(i) == "--snr") then
            read(argv(i+1), *) snr(1)
            read(argv(i+2), *) snr(2)
            i = i+3
            ! print *, "snr", snr
        elseif (argv(i) == "--bps") then
            read(argv(i+1),*) bps
            i = i + 2
            ! print *, "bps", bps
        elseif (argv(i) == "--minframeerr") then
            read(argv(i+1), *) min_ferr
            i = i + 2
            ! print *, "ferr", min_ferr
        elseif (argv(i) == "--minsimloops") then
            read(argv(i+1), *) min_sim
            i = i + 2
            ! print *, "min_sim", min_sim
        elseif (argv(i) == "--maxsimloops") then
            read(argv(i+1), *) max_sim
            i = i + 2
            ! print *, "max_sim", max_sim
        elseif (argv(i) == "--maxldpciter") then
            read(argv(i+1), *) max_iter
            i = i + 2
            ! print *, "ldpc", max_iter
        elseif (argv(i) == "--tannerfile") then
            call get_command_argument(i+1, tanner_file)
            i = i + 2
            ! print*, trim(tanner_file)
        elseif (argv(i) == "--outdir") then
            call get_command_argument(i+1, output_root)
            i = i + 2
        elseif (argv(i) == "--hard") then
            isHard = .true.
            i = i + 1
        elseif (argv(i) == "--alpha") then
            read(argv(i+1), *) alpha
            i = i + 2
        elseif (argv(i) == "-u") then
            uniform_th = .true.
            i = i + 1
        elseif (argv(i) == "--onlyinfo") then
            onlyinfo = .true.
            i = i + 1
        elseif (argv(i) == "-th") then
            tanner_header = .true.
            i = i + 1
        elseif (argv(i) == "--eve" ) then
            doEve = .true.
            i = i + 1
        elseif (argv(i) == "--interleaver") then
            useInterleaver = .true.
            i = i + 1
        else
            print *, "Unrecognized argument: ", argv(i)
            stop
        end if
    end do

    if (doEve) then
        ! Implies hard reconciliation
        isHard = .true.
    end if


    ! +------------------------------------+
    ! | Initialization of the random seeds |
    ! +------------------------------------+
    me   = this_image()
    n_im = num_images()
    call random_seed(get=seed)
    call sleep(me)
    call system_clock(time_seed)
    seed(1) = mod(seed(1)*sum(seed(me:)), abs(time_seed) + 2) - 12399027
    seed(2) = 467738 * me * (seed(3) - time_seed) + iand(sum(seed(:mod(me, 8))), time_seed)
    seed(3) = (time_seed + me) * (2 + mod(abs(883 + seed(mod(time_seed, 8) + 1)), me)) + &
        mod(sum(seed), max(maxval(seed(:)), 37))
    seed(4) = mod(987654321*time_seed, abs(time_seed - me*me*me) + 3)
    seed(5) = seed(7)/7 - time_seed*me + me ** mod(abs(seed(5)), 29)
    seed(6) = sum(seed(::2) ** abs(time_seed)) + 929812093 / sum(seed(1::2))
    seed(7) = me * seed(6) * 661242 - mod(seed(6) + time_seed**abs(sum(seed)), me)
    seed(8) = time_seed - seed(4)**me + sum(seed(2::2)) ** abs(time_seed) + 329999999/time_seed
    call stdlib_random_seed(sum(seed), time_seed) ! Necessary for the random functions in the stdlib
    call random_seed(put=[me, seed, time_seed])


    ! +----------------------------------------+
    ! | Allocation of simulation result arrays |
    ! +----------------------------------------+
    if (me==1) then
        allocate(outdata(nsnr, 3))
        ber => outdata(:, 2)
        fer => outdata(:, 3)
        outdata(:,2:) = 0
    else
        allocate(outdata(nsnr, 1))
    end if
    snrdb => outdata(:, 1)
    snrdb = [(snr(1) + real(i, dp)*(snr(2)-snr(1))/real(nsnr - 1, dp), i = 0, nsnr-1)]

    allocate(b_err(nsnr)[*])
    allocate(f_err(nsnr)[*])
    allocate(f_cnt(nsnr)[*])

    b_err(:) = 0
    f_err(:) = 0
    f_cnt(:) = 0

    critical
        call from_file(file=tanner_file, into=edge_definition, header=tanner_header)
    end critical


    ! +-----------------------------------------------------------------+
    ! | Allocation of Decoder, Noise Mapper and simulation data buffers |
    ! +-----------------------------------------------------------------+
    decoder = TDecoder(&
        edge_definition(1 , 1), &
        edge_definition(2:, 2), &
        edge_definition(2:, 3))

    deallocate(edge_definition)

    if (onlyinfo) then
        K = decoder%vnum - decoder%cnum
    else
        K = decoder%vnum
    end if

    allocate(x_i(decoder%vnum/bps))
    ! allocate(x(decoder%vnum/bps))
    allocate(y(decoder%vnum/bps))

    allocate(xhat(decoder%vnum/bps))
    allocate(nhat(decoder%vnum/bps))

    allocate(lappr(decoder%vnum))
    allocate(lappr_out(decoder%vnum))

    allocate(word(decoder%cnum))
    allocate(synd(decoder%cnum))

    nm = noisemapper_create(bps)
    if (.not. isHard) then
        call noisemapper_set_monotonicity(nm)
    end if

    ! +-------------------------------+
    ! | Display simulation parameters |
    ! +-------------------------------+
    if (me == 1) then
        call decoder%print
        print '("Each frame carries ", i6, " information bits.")', decoder%vnum - decoder%cnum
        print '("Performing ", i3, " iterations.")', max_iter
        print *, "SNR [dB] range:", snrdb(1), snrdb(size(snrdb))
        print '("Simulation loops: at least", i6, " up to ", i6, " or at least ", i3, " frame errors are found")', &
            min_sim, max_sim, min_ferr
        print '("Constellation has ", i3, " symbols, each carrying ", i3, " bits")', nm%M, nm%bps
        if (isHard) then
            print *, "Alice will use only HARD information"
        else
            print *, "Alice will use SOFT information"
        end if
        if (uniform_th) then
            print *, "Bob is using uniform probability quantization"
        else
            print *, "Bob is using maximum a-posteriori probability quantization"
        end if
    end if

    sync all

    ! +------------+
    ! | Simulation |
    ! +------------+
    if (me == 1) then
        call progress_bar%initialize(&
            filled_char_string='+', prefix_string='SNR points progress |',&
            suffix_string='| ', add_progress_percent=.true.)
        call progress_bar%start
    end if

    loop_snr : do i_snr = 1, nsnr
        call noisemapper_update_N0_from_snrdb(nm, snrdb(i_snr))
        ! if (uniform_th .or. (.not. isHard)) then
        !     call noisemapper_set_Fy_grids(nm)
        ! end if
        if (uniform_th) then
            call noisemapper_set_y_thresholds_uniform(nm)
        else
            call noisemapper_set_y_thresholds(nm)
        end if
        if (isHard) then
            call noisemapper_update_hard_reverse_tables(nm)
        else
            call noisemapper_set_Fy_grids(nm)
        end if

        if ( doEve ) then
            ! Use lappr(bps+1 : 2*bps) to temporarily store the denominator
            ! of the argument of the log, and lappr(1:bps) for the numerator
            lappr(1 : 2*bps) = 0
            do i = 1, bps
                do io = 0, nm%M-1
                    if (nm%s_to_b(io, i-1)) then
                        lappr(bps + i) = lappr(bps + i) + nm%delta_Fy(io)
                    else
                        lappr(i) = lappr(i) + nm%delta_Fy(io)
                    end if
                end do
            end do
            ! Store the lappr in the first bps locations
            lappr(1:bps) = log(lappr(1:bps)) - log(lappr(bps+1:2*bps))


            do i = 1, decoder%vnum/bps-1
                ! copy the lappr for the first symbol to all other symbols
                lappr(i*bps + 1 : (i+1)*bps) = lappr(1:bps)
            end do
        end if

        loop_frame : do i_frame = 1, max_sim
            ! Alice generates random symbols:
            call noisemapper_random_symbol(nm, x_i)
            y = noisemapper_symbol_index_to_value(nm, x_i)

            ! AWGN channel
            y    = rvs_normal(loc=y, scale=nm%sigma)

            ! Bob evaluates the soft metric, takes the decisions
            if (isHard) then
                xhat = noisemapper_decide_symbol(nm, y)
            else
                call noisemapper_generate_soft_metric(nm, y, nhat, xhat)
            end if
            word = noisemapper_symbol_to_word(nm, xhat)
            ! syndrome computation is postponed, so that the same shuffling
            ! applies to the word and the lappr arrays

            ! Alice uses the soft metric to find the output
            if (isHard) then
                if (.not. doEve) then
                    call noisemapper_convert_symbol_to_hard_lappr(nm, x_i, lappr)
                end if
            else
                call noisemapper_soft_reverse_lappr(nm, x_i, nhat, lappr, 1d-12)
            end if
            lappr = alpha*lappr

            if (useInterleaver) then
                call shuffle_word_and_lappr(word, lappr)
            end if
            synd = decoder%word_to_synd(word)
            N_iter = max_iter
            call decoder%decode(lappr, lappr_out, synd, N_iter)

            new_errors = count( (lappr_out(:K) < 0) .neqv. word(:K) )

            critical
                if (new_errors .gt. 0) then
                    b_err(i_snr)[1] = b_err(i_snr)[1] + new_errors
                    f_err(i_snr)[1] = f_err(i_snr)[1] + 1
                end if
                f_cnt(i_snr)[1]     = f_cnt(i_snr)[1] + 1
            end critical

            if (((f_err(i_snr)[1] .ge. min_ferr) .and. &
                (f_cnt(i_snr)[1] .ge. min_sim)) .or. (f_cnt(i_snr)[1] .ge. max_sim) ) then
                if (me==1) then
                    call progress_bar%update(current=real(i_snr, 8)/real(nsnr, 8))
                end if
                exit loop_frame
            end if
        end do loop_frame
        if (i_snr .ge. 3) then
            if (all(b_err(i_snr-2 : i_snr)[1] == 0)) then
                ! Check again after 10 seconds, so that if new errors pop up from other images, we keep helping them
                call sleep(10)
                if (all(b_err(i_snr-2 : i_snr)[1] == 0)) then
                    if (me==1) then
                        call progress_bar%update(current=1d0)
                        call progress_bar%destroy
                    end if
                    exit loop_snr
                end if
            end if
        end if
    end do loop_snr

    ! +-----------------------------------------+
    ! | First Image saves the simulation result |
    ! +-----------------------------------------+
    sync all
    if (me == 1) then
        if (alpha == 1d0) then
            call make_directory_and_file_name(output_root, bps, .true., isHard, &
                snr, nsnr, min_sim, max_sim, max_iter, min_ferr,         &
                output_dir, output_name)
        else
            call make_directory_and_file_name(output_root, bps, .true., isHard, &
                snr, nsnr, min_sim, max_sim, max_iter, min_ferr,         &
                output_dir, output_name, alpha)
            print *, trim(output_dir)
            print *, trim(output_name)
        end if
        call execute_command_line("mkdir -p " // trim(output_dir))
        open(newunit=io, file=trim(output_dir) // "/" // trim(output_name) // ".log", &
            status="replace", action="write")

        write(io, '(A, T16, A, T32, A, T48, A, T64, A, T80, A)') "SNR [dB]", "F_NUM", "ERR", "BER", "FERR", "FER"
        do i_snr = 1, nsnr
            if (b_err(i_snr) == 0) then
                ber(i_snr) = 0
                fer(i_snr) = 0
            else
                ber(i_snr) = real(b_err(i_snr), dp)/real(f_cnt(i_snr), dp)/real(K, dp)
                fer(i_snr) = real(f_err(i_snr), dp)/real(f_cnt(i_snr), dp)
            end if

            write(io, '(f12.3, T16, I6, T32, I10, T48, ES10.3E3, T64, I10, T80, ES10.3E3)') &
                snrdb(i_snr), f_cnt(i_snr), b_err(i_snr), ber(i_snr), f_err(i_snr), fer(i_snr)
        end do
        close(io)

        ! call to_file(x=outdata, file=output_file, header=["SNR", "BER", "FER"], fmt="f")
        if (alpha == 1d0) then
            call save_data(outdata, output_root, bps, .true., isHard, snr, nsnr, min_sim, max_sim, max_iter, min_ferr)
        else
            call save_data(outdata, output_root, bps, .true., isHard, snr, nsnr, min_sim, max_sim, max_iter, min_ferr, alpha)
        end if
    end if

contains

    subroutine shuffle_word_and_lappr(word, lappr)
        logical, intent(inout) :: word(0:)
        double precision, intent(inout) :: lappr(0:size(word)-1)

        double precision :: tmp
        integer :: i, n, j

        n = size(word)

        do i = 0, n-2 ! Number of currently shuffled elements
            ! Number of unshaffled elements is n-i
            j = i + rvs_uniform(n-1-i) ! Pick one of the remaining unshaffled elements
            if (i /= j) then
                ! Swap word bits
                word(i) = word(i) .neqv. word(j)
                word(j) = word(i) .neqv. (.not. word(j))
                word(i) = word(i) .neqv. (.not. word(j))

                ! Swap lappr data
                tmp = lappr(i)
                lappr(i) = lappr(j)
                lappr(j) = tmp
            end if
        end do
    end subroutine shuffle_word_and_lappr
end program reverse_reconciliation
