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
    use iso_fortran_env, only: dp => real64
    use ldpc_decoder, only: TDecoder
    use re2often_noise_mapper
    use forbear, only: bar_object
    use stdlib_random, only: stdlib_random_seed => random_seed
    use stdlib_stats_distribution_normal, only: rvs_normal
    implicit none

    integer :: argc
    character(len=500), allocatable :: argv(:)

    character(500) :: tanner_file
    character(500) :: output_file
    double precision :: snr(2)         ! Signal to Noise Ratio in dB
    integer :: nsnr           ! Number of SNR points
    integer :: bps            ! Bits per symbol
    integer :: min_ferr       ! Minimum frame errors
    integer :: max_sim        ! Maximum simulation loops
    integer :: min_sim        ! Minimum simulation loops
    integer :: max_iter       ! Maximum number of LDPC iterations

    integer, allocatable :: edge_definition(:,:)
    double precision, allocatable, target :: outdata(:,:)

    double precision, pointer :: snrdb(:)
    double precision, pointer :: ber(:)
    double precision, pointer :: fer(:)

    integer, allocatable :: b_err(:)[:]
    integer, allocatable :: f_err(:)[:]
    integer, allocatable :: f_cnt(:)[:]

    integer, allocatable :: x_i(:)
    double precision, allocatable :: x(:), y(:), lappr(:), lappr_out(:)
    logical, allocatable :: word(:), synd(:)
    double precision :: new_errors
    integer :: K, N_iter

    type(TDecoder)     :: decoder
    type(TNoiseMapper) :: nm

    double precision   :: sigma

    integer :: i, i_snr, i_frame
    integer :: me, n_im
    integer :: seed(8), time_seed

    type(bar_object) :: progress_bar


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
    tanner_file = ""
    output_file = ""

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
        elseif (argv(i) == "--csvout") then
            call get_command_argument(i+1, output_file)
            i = i + 2
            ! print *, trim(output_file)
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
    seed(2) = 467738 * me * (seed(3) - time_seed) + iand(sum(seed(:me)), time_seed)
    seed(3) = (time_seed + me) * (2 + mod(abs(883 + seed(mod(time_seed, 8) + 1)), me)) + &
        mod(sum(seed), max(maxval(seed(:)), 37))
    seed(4) = mod(987654321*time_seed, abs(time_seed - me*me*me) + 3)
    seed(5) = seed(7)/7 - time_seed*me + me ** mod(abs(seed(5)), 29)
    seed(6) = sum(seed(::2) ** abs(time_seed)) + 929812093 / sum(seed(1::2))
    seed(7) = me * seed(6) * 661242 - me / (seed(6) + time_seed**abs(sum(seed)))
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
        call from_file(file=tanner_file, into=edge_definition, header=.true.)
    end critical

    decoder = TDecoder(&
        edge_definition(1 , 1), &
        edge_definition(2:, 2), &
        edge_definition(2:, 3))

    K = decoder%vnum - decoder%cnum

    allocate(x_i(decoder%vnum/bps))
    allocate(x(decoder%vnum/bps))
    allocate(y(decoder%vnum/bps))

    allocate(lappr(decoder%vnum))
    allocate(lappr_out(decoder%vnum))

    allocate(word(decoder%cnum))
    allocate(synd(decoder%cnum))

    nm = TNoiseMapper(bps=bps, N0=1d0)

    if (me == 1) then
       call progress_bar%initialize(&
            filled_char_string='+', prefix_string='SNR points progress |',&
            suffix_string='| ', add_progress_percent=.true.)
       call progress_bar%start
    end if

    loop_snr : do i_snr = 1, nsnr
        call nm%update_N0(nm%E_symbol * 10d0**(-snrdb(i_snr)/10d0))
        sigma = sqrt(nm%N0/2)

        loop_frame : do i_frame = 1, max_sim
            x_i  = nm%random_symbols()
            x    = nm%symbol_index_to_value(x_i)
            word = nm%symbol_to_word(x_i)
            synd = decoder%word_to_synd(word)

            y    = rvs_normal(loc=x, scale=sigma)

            lappr = nm%y_to_lappr(y)

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

            if ((f_err(i_snr)[1] .ge. min_ferr) .and. &
                (f_cnt(i_snr)[1] .ge. min_sim) ) then
                if (me==1) then
                    call progress_bar%update(current=real(i_snr, 8)/real(nsnr, 8))
                end if
                exit loop_frame
            end if
        end do loop_frame
    end do loop_snr

    sync all

    if (me == 1) then
        print '(A, T16, A, T32, A, T48, A, T64, A, T80, A)', "SNR [dB]", "F_NUM", "ERR", "BER", "FERR", "FER"
        do i_snr = 1, nsnr
            if (b_err(i_snr) == 0) then
                ber(i_snr) = 0
                fer(i_snr) = 0
            else
                ber(i_snr) = real(b_err(i_snr), dp)/real(f_cnt(i_snr), dp)/real(K, dp)
                fer(i_snr) = real(f_err(i_snr), dp)/real(f_cnt(i_snr), dp)
            end if

            print '(f6.3, T16, I10, T32, I10, T48, ES10.3E3, T64, I10, T80, ES10.3E3)', &
                snrdb(i_snr), f_cnt(i_snr), b_err(i_snr), ber(i_snr), f_err(i_snr), fer(i_snr)
        end do

        call to_file(x=outdata, file=output_file, header=["SNR", "BER", "FER"], fmt="f")
    end if

end program
