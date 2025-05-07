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
program direct_reconciliation
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Perform direct reconciliation
    use io_fortran_lib, only: from_file, to_file
    use iso_c_binding, only: dp => c_double
    use ldpc_decoder, only: TDecoder
    use re2often_noisemapper
    use re2often_utils, only: save_data, make_directory_and_file_name
    use forbear, only: bar_object
    use stdlib_random, only: stdlib_random_seed => random_seed
    use stdlib_stats_distribution_normal, only: rvs_normal
    use stdlib_stats_distribution_uniform, only: rvs_uniform
    implicit none

    integer :: argc
    character(len=500), allocatable :: argv(:)

    character(500) :: tanner_file
    character(500) :: output_root
    character(250) :: output_dir
    character(250) :: output_name
    integer :: io
    double precision :: snr(2)         ! Signal to Noise Ratio in dB
    integer :: nsnr           ! Number of SNR points
    integer :: bps            ! Bits per symbol
    integer :: min_ferr       ! Minimum frame errors
    integer :: max_sim        ! Maximum simulation loops
    integer :: min_sim        ! Minimum simulation loops
    integer :: max_iter       ! Maximum number of LDPC iterations
    logical :: tanner_header  ! Whether tanner file has a header
    logical :: onlyinfo       ! Whether to compare only the first N-M bits, instead of whole frame
    logical :: useInterleaver ! Apply random shuffling to bits

    integer, allocatable :: edge_definition(:,:)
    double precision, allocatable, target :: outdata(:,:)

    double precision, pointer :: snrdb(:)
    double precision, pointer :: ber(:)
    double precision, pointer :: fer(:)

    integer, allocatable :: b_err(:)[:]
    integer, allocatable :: f_err(:)[:]
    integer, allocatable :: f_cnt(:)[:]

    integer(c_int), allocatable :: x_i(:)
    ! real(c_double), allocatable :: x(:), y(:), lappr(:), lappr_out(:)
    real(c_double), allocatable :: y(:), lappr(:), lappr_out(:)
    logical, allocatable :: word(:), synd(:)
    integer :: new_errors
    integer :: K, N_iter

    type(TDecoder)     :: decoder
    type(noisemapper_type) :: nm

    integer :: i, i_snr, i_frame
    integer :: me, n_im
    integer :: seed(8), time_seed

    type(bar_object) :: progress_bar

    if (this_image()==1) then
        print *, " +----------------------------------+"
        print *, " | DIRECT reconciliation simulation |"
        print *, " +----------------------------------+"
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
    snr  = [-3, -2]
    bps  = 1
    min_ferr = 50
    max_sim  = 5000
    min_sim  = 250
    max_iter = 50
    tanner_file = "./assets/codes/dvbs2ldpc0.500.csv"
    output_root = "./res/rate1d2"
    onlyinfo = .false.
    tanner_header = .false.
    useInterleaver = .false.


    i = 1
    do while(i <= argc)
        if (argv(i) == "--nsnr") then
            read(argv(i+1),*) nsnr
            i = i + 2
        elseif (argv(i) == "--snr") then
            read(argv(i+1), *) snr(1)
            read(argv(i+2), *) snr(2)
            i = i+3
        elseif (argv(i) == "--bps") then
            read(argv(i+1),*) bps
            i = i + 2
        elseif (argv(i) == "--minframeerr") then
            read(argv(i+1), *) min_ferr
            i = i + 2
        elseif (argv(i) == "--minsimloops") then
            read(argv(i+1), *) min_sim
            i = i + 2
        elseif (argv(i) == "--maxsimloops") then
            read(argv(i+1), *) max_sim
            i = i + 2
        elseif (argv(i) == "--maxldpciter") then
            read(argv(i+1), *) max_iter
            i = i + 2
        elseif (argv(i) == "--tannerfile") then
            call get_command_argument(i+1, tanner_file)
            i = i + 2
        elseif (argv(i) == "--outdir") then
            call get_command_argument(i+1, output_root)
            i = i + 2
        elseif (argv(i) == "--onlyinfo") then
            onlyinfo = .true.
            i = i + 1
        elseif (argv(i) == "-th") then
            tanner_header = .true.
            i = i + 1
        elseif (argv(i) == "--interleaver") then
            useInterleaver = .true.
            i = i + 1
        else
            print *, "Unrecognized argument: ", argv(i)
            stop
        end if
    end do

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
    seed(7) = me * seed(6) * 661242 - (seed(6) + time_seed**abs(sum(seed))) / me
    seed(8) = time_seed - seed(4)**me + sum(seed(2::2)) ** abs(time_seed) + 329999999/time_seed
    call stdlib_random_seed(sum(seed), time_seed) ! Necessary for the random functions in the stdlib
    call random_seed(put=[me, seed, time_seed])
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

    decoder = TDecoder(&
        edge_definition(1 , 1), &
        edge_definition(2:, 2), &
        edge_definition(2:, 3))

    deallocate(edge_definition) ! Save up some memory unused memory

    if (onlyinfo) then
        K = decoder%vnum - decoder%cnum
    else
        K = decoder%vnum
    end if

    allocate(x_i(decoder%vnum/bps))
    ! allocate(x(decoder%vnum/bps))
    allocate(y(decoder%vnum/bps))

    allocate(lappr(decoder%vnum))
    allocate(lappr_out(decoder%vnum))

    allocate(word(decoder%vnum))
    allocate(synd(decoder%cnum))

    nm = noisemapper_create(bps)

    if (me == 1) then
        call progress_bar%initialize(&
            filled_char_string='+', prefix_string='SNR points progress |',&
            suffix_string='| ', add_progress_percent=.true.)
        call progress_bar%start
    end if

    sync all

    loop_snr : do i_snr = 1, nsnr
        call noisemapper_update_N0_from_snrdb(nm, snrdb(i_snr))

        loop_frame : do i_frame = 1, max_sim
            call noisemapper_random_symbol(nm, x_i)
            y    = noisemapper_symbol_index_to_value(nm, x_i)
            word = logical(noisemapper_symbol_to_word(nm, x_i)) ! from logical 1 to logical 4

            y    = rvs_normal(loc=y, scale=nm%sigma)

            call noisemapper_y_to_lappr(nm, y, lappr)

            if (useInterleaver) then
                call shuffle_word_and_lappr(word, lappr)
            end if
            synd = decoder%word_to_synd(word)

            N_iter = max_iter
            call decoder%decode(lappr, lappr_out, synd, N_iter)

            new_errors = count( (lappr_out(:K) < 0) .neqv. word(:K) )

            critical
                if (new_errors > 0) then
                    b_err(i_snr)[1] = b_err(i_snr)[1] + new_errors
                    f_err(i_snr)[1] = f_err(i_snr)[1] + 1
                end if
                f_cnt(i_snr)[1]     = f_cnt(i_snr)[1] + 1
            end critical

            if (((f_err(i_snr)[1] .ge. min_ferr) .and. &
                (f_cnt(i_snr)[1] .ge. min_sim)) .or. &
                (f_cnt(i_snr)[1] .ge. max_sim)) then
                if (me==1) then
                    call progress_bar%update(current=real(i_snr, 8)/real(nsnr, 8))
                end if
                exit loop_frame
            end if
        end do loop_frame
        if ((i_snr .ge. 2)) then
            if (all(b_err(i_snr-1 : i_snr)[1] == 0)) then
                ! Check again after 10 seconds, so that if new errors pop up from other images, we keep helping them
                call sleep(10)
                if (all(b_err(i_snr-1 : i_snr)[1] == 0)) then
                    if (me==1) then
                        call progress_bar%update(current=1d0)
                    end if
                    exit loop_snr
                end if
            end if
        end if
    end do loop_snr

    sync all

    if (me == 1) then
        call progress_bar%destroy
        call make_directory_and_file_name(output_root, bps, .false., .false., &
            snr, nsnr, min_sim, max_sim, max_iter, min_ferr,         &
            output_dir, output_name)
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

            write(io,  '(f12.3, T10, I10, T32, I10, T48, ES10.3E3, T64, I10, T80, ES10.3E3)') &
                snrdb(i_snr), f_cnt(i_snr), b_err(i_snr), ber(i_snr), f_err(i_snr), fer(i_snr)
        end do
        close(io)

        ! call to_file(x=outdata, file=output_file, header=["SNR", "BER", "FER"], fmt="f")
        call save_data(outdata, output_root, bps, .false., .false., snr, nsnr, min_sim, max_sim, max_iter, min_ferr)
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
                word(j) = word(i) .neqv. word(j)
                word(i) = word(i) .neqv. word(j)

                ! Swap lappr data
                tmp = lappr(i)
                lappr(i) = lappr(j)
                lappr(j) = tmp
            end if
        end do
    end subroutine shuffle_word_and_lappr
end program
