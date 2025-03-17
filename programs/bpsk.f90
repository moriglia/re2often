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
program bpsk
    use io_fortran_lib, only: from_file, to_file
    use iso_c_binding, only: wp => c_double, kl=> c_bool
    use ldpc_decoder, only: TDecoder
    ! use forbear, only: bar_object
    use stdlib_random, only: stdlib_random_seed => random_seed
    use stdlib_stats_distribution_normal, only: rvs_normal
    use re2often_utils
    implicit none

    integer     :: argc
    character(len=500), allocatable :: argv(:)
    character(len=500) :: tanner_file
    character(len=500) :: output_root
    logical     :: tanner_has_header

    real(wp)    :: snr(2)         ! Signal to Noise Ratio in dB
    integer     :: nsnr           ! Number of SNR points
    integer     :: min_ferr       ! Minimum frame errors
    integer     :: max_sim        ! Maximum simulation loops
    integer     :: min_sim        ! Minimum simulation loops
    integer     :: max_iter       ! Maximum number of LDPC iterations

    integer, allocatable :: edge_definition(:,:)
    integer              :: Ne
    real(wp), allocatable, target :: outdata(:,:)

    real(wp), pointer :: snrdb_array(:), ber(:), fer(:)

    integer :: i, me, n_images

    type(TDecoder) :: decoder
    integer :: time_seed
    integer:: seed(8)

    integer :: i_snr, i_frame
    integer :: N_symb, K

    real(wp) :: Es, N0, N0_half, sigma, N0_half_est


    integer, allocatable  :: x_i(:)
    real(wp), allocatable :: x(:)
    real(wp), allocatable :: y(:)
    logical, allocatable  :: word(:)
    logical, allocatable  :: synd(:)
    real(wp), allocatable :: lappr(:)
    real(wp), allocatable :: lappr_updated(:)

    integer :: new_errors, n_ldpc_it

    integer, allocatable :: b_error_count(:)[:]
    integer, allocatable :: f_error_count(:)[:]
    integer, allocatable :: f_count(:)[:]

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
    min_ferr = 50
    max_sim  = 5000
    min_sim  = 250
    max_iter = 50
    tanner_file = "./assets/codes/dvbs2ldpc0.500.csv"
    tanner_has_header = .false.
    output_root = "./res/rate1d2"

    i = 1
    do while(i <= argc)
        if (argv(i) == "--nsnr") then
            read(argv(i+1),*) nsnr
            i = i + 2
        elseif (argv(i) == "--snr") then
            read(argv(i+1), *) snr(1)
            read(argv(i+2), *) snr(2)
            i = i+3
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
        elseif(argv(i) == "-th") then
            tanner_has_header = .true.
            i = i + 1
        else
            print *, "Unrecognized argument: ", argv(i)
            stop
        end if
    end do

    deallocate(argv)
    allocate(outdata(nsnr, 3))


    ! type(bar_object) :: progress_bar

    snrdb_array => outdata(:,1)
    ber         => outdata(:,2)
    fer         => outdata(:,3)

    me       = this_image()
    n_images = num_images()


    call random_seed(get=seed)
    call sleep(me)
    call system_clock(time_seed)
    seed(1) = mod(seed(1)*sum(seed(me:)), abs(time_seed) + 2) - 12399027
    seed(2) = 467738 * me * (seed(3) - time_seed) + iand(sum(seed(:mod(me, 8)+1)), time_seed)
    seed(3) = (time_seed + me) * (2 + mod(abs(883 + seed(mod(time_seed, 8) + 1)), me)) + &
        mod(sum(seed), max(maxval(seed(:)), 37))
    seed(4) = mod(987654321*time_seed, abs(time_seed - me*me*me) + 3)
    seed(5) = seed(7)/7 - time_seed*me + me ** mod(abs(seed(5)), 29)
    seed(6) = sum(seed(::2) ** abs(time_seed)) + 929812093 / sum(seed(1::2))
    seed(7) = me * seed(6) * 661242 - me / (seed(6) + time_seed**abs(sum(seed)))
    seed(8) = time_seed - seed(4)**me + sum(seed(2::2)) ** abs(time_seed) + 329999999/time_seed
    call stdlib_random_seed(sum(seed), time_seed) ! Necessary for the random functions in the stdlib
    call random_seed(put=[me, seed, time_seed])

    critical
        call from_file(&
            file=tanner_file, &
            into=edge_definition, &
            header=tanner_has_header)
    end critical

    decoder = TDecoder(&
        edge_definition(1 , 1), &
        edge_definition(2:, 2), &
        edge_definition(2:, 3))

    if (me==1) then
        print *, edge_definition(1,:)
        call decoder%print
    end if

    deallocate(edge_definition)

    K = decoder%vnum - decoder%cnum

    allocate(x(decoder%vnum))
    allocate(x_i(decoder%vnum))
    allocate(y(decoder%vnum))
    allocate(word(decoder%vnum))
    allocate(synd(decoder%cnum))
    allocate(lappr(decoder%vnum))
    allocate(lappr_updated(decoder%vnum))

    Es = 1

    allocate(b_error_count(nsnr)[*])
    allocate(f_error_count(nsnr)[*])
    allocate(f_count(nsnr)[*])

    b_error_count(:) = 0
    f_error_count(:) = 0
    f_count(:)       = 0

    snrdb_array(:) = [(snr(1) + real(i, wp)*(snr(2) - snr(1))/real(nsnr-1, wp), i=0, nsnr-1)]
    ! if (me == 1) then
    !    call progress_bar%initialize(&
    !         filled_char_string='+', prefix_string='SNR points progress |',&
    !         suffix_string='| ', add_progress_percent=.true.)
    !    call progress_bar%start
    ! end if

    snr_loop : do i_snr = 1, nsnr
        N0 = (10_wp ** (-snrdb_array(i_snr)/10_wp) ) * Es
        N0_half = N0/2_wp
        sigma = sqrt(N0_half)

        print '("[",I2,"] SNR[dB]=",f6.2,3x,"N_0/2=",f9.6)', me, snrdb_array(i_snr), N0_half

        N0_half_est = 0

        frame_loop: do i_frame = 1, max_sim
            call random_number(x)
            do i = 1, decoder%vnum
                if (x(i) < 0.5) then
                    x(i) = -1_wp
                else
                    x(i) = 1_wp
                end if
            end do

            word = x < 0
            synd = decoder%word_to_synd(word)

            y = rvs_normal(loc=x, scale=sigma)
            lappr = 4*y/N0
            N0_half_est = N0_half_est + sum((y-x)**2) / real(decoder%vnum, wp)

            n_ldpc_it = max_iter
            call decoder%decode(lappr, lappr_updated, synd, n_ldpc_it)

            new_errors = count((lappr_updated < 0) .neqv. word)

            critical
                if (new_errors > 0) then
                    f_error_count(i_snr)[1] = f_error_count(i_snr)[1] + 1
                    b_error_count(i_snr)[1] = b_error_count(i_snr)[1] + new_errors
                end if
                f_count(i_snr)[1] = f_count(i_snr)[1] + 1
            end critical

            if ((f_count(i_snr)[1] > min_sim) .and. (f_error_count(i_snr)[1] > min_ferr)) then
                exit frame_loop
            end if
            if (f_count(i_snr)[1] >= max_sim) exit frame_loop
        end do frame_loop

        N0_half_est = N0_half_est / real(i_frame, wp)
        print '("[",I2,"] COMPLETED SNR[dB]=",f6.2,3x,"N_0/2=",f9.6," (estimated)   DELTA=", f9.6)', &
            me, snrdb_array(i_snr), N0_half_est, abs(N0_half - N0_half_est)

        ! if (me==1) then
        !    call progress_bar%update(current=real(i_snr, 8)/real(nsnr, 8))
        ! end if
    end do snr_loop

    sync all

    if (me == 1) then
        print '(A, T16, A, T32, A, T48, A, T64, A, T80, A)', "SNR [dB]", "F_NUM", "ERR", "BER", "FERR", "FER"
        do i_snr = 1, nsnr
            if (b_error_count(i_snr) == 0) then
                ber(i_snr) = 0.0_wp
                fer(i_snr) = 0.0_wp
            else
                ber(i_snr) = real(b_error_count(i_snr), wp)/real(f_count(i_snr), wp)/real(K, wp)
                fer(i_snr) = real(f_error_count(i_snr), wp)/real(f_count(i_snr), wp)
            end if

            print '(f12.3, T14, I10, T32, I10, T48, ES10.3E3, T64, I10, T80, ES10.3E3)', &
                snrdb_array(i_snr), f_count(i_snr), b_error_count(i_snr), ber(i_snr), f_error_count(i_snr), fer(i_snr)
        end do

        call save_data(outdata, trim(output_root)//"/bpsk", 1, &
            .false., .false., snr, nsnr, min_sim, max_sim, max_iter, min_ferr)
        ! call to_file(x=outdata, file="bpsk_50it.csv", header=["SNR", "BER", "FER"], fmt="f")
    end if


end program bpsk
