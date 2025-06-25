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
program mi_opt
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Find optimal parameters for the MI
    use iso_fortran_env, only: &
        lock_type, &
        dp => real64, &
        stdout => output_unit, &
        event_type
    use io_fortran_lib, only: from_file, to_file
    use re2often_noisemapper
    use re2often_utils, only: save_data, make_directory_and_file_name
    ! use forbear, only: bar_object
    use re2often_mi ! defines a noisemapper_type object
    use lincoa_mod
    use quadpack, only: dqags
    use flap ! CLI parser: command_line_interface
    implicit none

    ! +---------------------+
    ! | Command line inputs |
    ! +---------------------+
    integer :: argc
    character(len=500), allocatable :: argv(:)

    character(250) :: output_root
    character(250) :: output_dir
    character(250) :: output_name
    character(20), allocatable   :: header(:)
    integer        :: io_log, io_csv
    double precision :: snr(2)         ! Signal to Noise Ratio in dB
    integer :: nsnr           ! Number of SNR points
    integer :: bps            ! Bits per symbol
    logical :: isReverse      ! Whether to calculate the M.I. for the reverse reconciliation
    logical :: isHard         ! Whether to perform hard reverse reconciliation
    logical :: isGMI          ! Compute the GMI instead of the MI

    integer, allocatable :: monoConfig(:)  ! Monotonicity configuration number list
    integer, allocatable :: encConfig(:,:) ! Encoding configuration list, one per row
    character(250) :: monoFile ! File with monotonicity configurations
    character(250) :: encFile  ! File with symbol permutations

    type(command_line_interface) :: cli
    integer :: error

    ! +-------------------+
    ! | Optimization data |
    ! +-------------------+
    double precision, allocatable, target :: Aineq(:,:)
    double precision, allocatable, target :: bineq(:)
    double precision, allocatable :: opt_array(:)
    double precision, pointer :: A(:,:)
    double precision, pointer :: b(:)
    double precision :: I ! the mutual information !

    ! +-------------+
    ! | Output data |
    ! +-------------+
    double precision, allocatable, target :: outdata(:,:)[:]
    integer :: encoding_index[*]
    integer :: monoconf_index[*]
    integer :: snrdb_index[*]
    integer, parameter :: o_mi = 3
    integer, parameter :: o_th = 4
    integer            :: o_p  ! = o_th + nm%M-1
    integer            :: o_s  ! = o_p + nm%M
    integer            :: o_e
    integer            :: o_c
    type(event_type), allocatable :: snr_done(:)[:]
    integer            :: until_count
    character(500)     :: format_csv
    character(500)     :: format_log
    character(50)      :: tmpstr

    double precision, pointer :: snr_array(:)

    ! +--------------------+
    ! | Computed constants |
    ! +--------------------+
    integer :: M_half


    ! +---------------------+
    ! | Image-specific data |
    ! +---------------------+
    integer :: me, n_im

    type(lock_type) :: lck[*]


    ! generic iteration variables
    integer :: i_snr, ii, i_encoding, i_config

    ! --- Code -----------------------------------------------

#define NOT_IMPLEMENTED() \
    print *, "Not Implemented" ; \
    stop

    me = this_image()
    n_im = num_images()

    if (me == 1) then
        print *, " +--------------------------+"
        print *, " | Optimal MI / GMI program |"
        print *, " +--------------------------+"
    end if

    ! +------------+
    ! | CLI parser |
    ! +------------+
    call cli%init(&
        progname = "mi_opt", &
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
        switch_ab = '-hr', &
        help='Hard reconciliation (implies "--reverse")', &
        required=.false., &
        act='store_true', &
        def='.false.', &
        error=error)
    call cli%add(switch='--reverse', &
        switch_ab='-r', &
        help='Reverse Reconciliation', &
        required=.false., &
        act='store_true', &
        def='.false.', &
        error=error)
    call cli%add(switch='--config-file', &
        help='File containing one configuration number per line', &
        required=.false., &
        act='store', &
        def='0', &
        error=error)
    call cli%add(switch='--encoding-file', &
        help='CSV file containing one permutation of the symbol indexes per line', &
        required=.false., &
        act='store', &
        def='0', &
        error=error)
    call cli%add(switch='--gmi', &
        help='Compute GMI instead of MI', &
        required=.false., &
        act='store_true', &
        def='.false.', &
        error=error)

    call cli%parse(error=error)
    if (error /= 0) stop

    ! Retrieve values
    call cli%get(switch='--nsnr', val=nsnr)
    call cli%get(switch='--snr', val=snr)
    call cli%get(switch='--bps', val=bps)
    call cli%get(switch='--outdir', val=output_root)
    call cli%get(switch='--hard', val=isHard)
    call cli%get(switch='--reverse', val=isReverse)
    call cli%get(switch='--config-file', val=monoFile)
    call cli%get(switch='--encoding-file', val=encFile)
    call cli%get(switch='--gmi', val=isGMI)

    if (isHard) then
        isReverse = .true. ! --hard implies --reverse
    end if

    if (.not. isReverse) then
        NOT_IMPLEMENTED()
    end if

    if (isGMI) then
        if (.not. cli%is_passed(switch="--encoding-file")) then
            print *, "--encoding-file <encoding.csv> is mandatory when --gmi is active"
            stop
        end if
        critical
            call from_file(file=encFile, into=encConfig)
        end critical
    end if
    if (isReverse .and. (.not. isHard)) then
        if (.not. cli%is_passed(switch="--config-file")) then
            print *, "--config-file <config.csv> is mandatory when -r is active and --hard isn't"
            stop
        end if
        critical
            call from_file(file=monoFile, into=monoConfig)
        end critical
    end if

    nm = noisemapper_create(bps)

    M_half = ishft(nm%M, -1)
    o_p = o_th + nm%M-1
    o_s = o_p + nm%M
    o_e = o_s + 1
    o_c = o_e + nm%M
    allocate(outdata(nsnr, o_c)[*])
    allocate(snr_done(nsnr)[*])
    if (me==1) then
        allocate(header(o_c))

        encoding_index = 1
        monoconf_index = 1
        snrdb_index    = 1

        outdata(:, 2 ) = 0
        outdata(:, 4:) = 0
        outdata(:, 3 ) = -1d0


        header(1) = trim("SNR")
        header(2) = trim("E_b/N_0")
        header(3) = trim("I")
        header(o_th) = "\theta"
        call write_header(header(o_th), 1, nm%M-1, header(o_th:))
        header(o_p) = "P"
        call write_header(header(o_p), 0, nm%M-1, header(o_p:))
        header(o_s) = "s"
        header(o_e) = "B"
        call write_header(header(o_e), 0, nm%M-1, header(o_e:))
        header(o_c) = "C"
    end if

    snr_array => outdata(:,1)
    snr_array = [(snr(1) + i_snr*(snr(2) - snr(1))/(nsnr-1), i_snr = 0, nsnr-1)]


    allocate(Aineq(M_half, M_half))
    allocate(bineq(M_half))
    allocate(opt_array(M_half))
    Aineq(:,:) = 0
    do ii = 1, M_half
        Aineq(ii,ii) = -1d0
    end do
    do ii = 2, M_half-1
        ! Note that on the last row we have the bound for "s"
        ! Which is decoupled from the other equations
        Aineq(ii, ii-1) = 1d0
    end do
    bineq(:) = -1d-2

    if (isGMI) then
        A => Aineq
        b => bineq
    else
        ! Exclude the s parameter from the optimization
        A => Aineq(1:M_half-1, 1:M_half-1)
        b => bineq(1:M_half-1)
    end if

    if (me == 1) then
        if (isGMI) then
            output_root = trim(output_root)//"/opt-gmi"
        else
            output_root = trim(output_root)//"/opt-mi"
        end if

        call make_directory_and_file_name(output_root, bps, isReverse, isHard, &
            snr, nsnr, 0, 0, 0, 0, output_dir, output_name)
        call execute_command_line("mkdir -p " // trim(output_dir))

        open(newunit=io_log, file=trim(output_dir) // "/" // trim(output_name) // ".log", &
            status="replace", action="write")
        open(newunit=io_csv, file=trim(output_dir) // "/" // trim(output_name) // ".csv", &
            status="replace", action="write")

        format_log = '(f12.8, 4x, f12.9, 4x, E12,'
        format_csv = '(f0, ",", f0, ",", E,'
        write(io_log, '(A12, 4x, A12, 4x, A12)', advance='no') trim(header(1)), trim(header(2)), trim(header(3))
        write(io_csv, '(A, ",", A, ",", A)', advance='no') trim(header(1)), trim(header(2)), trim(header(3))

        write(format_log, '(A, " ", i0, "(4x, f12.8), ", i0, "(4x, f12.8)" )') trim(format_log), nm%M-1, nm%M
        write(format_csv, '(A, " ", i0, A, i0, A )') &
            trim(format_csv), nm%M-1, '(",", f), ', nm%M, '(",", f)'
        do ii = 1, nm%M-1
            write(tmpstr, '("\theta_", i0)') ii
            write(io_log, '(4x, A12)', advance='no') trim(tmpstr)
            write(io_csv, '(",", A)', advance='no') trim(tmpstr)
        end do
        do ii = 1, nm%M
            write(tmpstr, '("P_", i0)') ii
            write(io_log, '(4x, A12)', advance='no') trim(tmpstr)
            write(io_csv, '(",", A)', advance='no') trim(tmpstr)
        end do

        if (isGMI) then
            write(format_log, '(A)') trim(format_log)//', 4x, E12'
            write(format_csv, '(A)') trim(format_csv)//', ",", E'
            write(format_log, '(A, i0, A)') trim(format_log)//', ', nm%M, '(4x,  i0)'
            write(format_csv, '(A, i0, A)') trim(format_csv)//', ', nm%M, '(",", i0)'

            write(io_log, '(4x, A12)', advance='no') "s"
            write(io_csv, '(",", A)' , advance='no') "s"
            do ii = 0, nm%M-1
                write(tmpstr, '("B_", i0)') ii
                write(io_log, '(4x, A12)', advance='no') trim(tmpstr)
                write(io_csv, '(",", A)'  , advance='no')  trim(tmpstr)
            end do
        end if
        if (isReverse .and. .not. isHard) then
            write(format_log, '(A)') trim(format_log)//', 4x, i0'
            write(format_csv, '(A)') trim(format_csv)//', ",", i0'

            write(io_log, '(4x, A12)', advance='no') "C"
            write(io_csv, '(",",  A)', advance='no') "C"
        end if
        write(format_log, '(A)') trim(format_log)//')'
        write(format_csv, '(A)') trim(format_csv)//')'
        write(io_log, '(A)') ""
        write(io_csv, '(A)') ""
        print *, trim(format_csv)
        print *, trim(format_log)

        flush(io_log)
        flush(io_csv)

        if ( isGMI ) then
            until_count = size(encConfig, 1)
        else
            until_count = 1
        end if
        if ( isReverse .and. .not. isHard ) then
            until_count = until_count * size(monoConfig)
        end if
    end if

    sync all


    if ((me==1) .and. (n_im > 1)) then
        do i_snr = 1, nsnr
            event wait( snr_done(i_snr), until_count=until_count )
            call write_result_to_file(i_snr)
        end do
        goto 100 ! The end :D
    end if

    loop_snr : do while (.true.)
        lock(lck[1])
        ! Exit if we are about to get an SNR index above nsnr
        if (snrdb_index[1] .gt. nsnr) then
            unlock(lck[1])
            exit loop_snr
        end if
        ! Actually get the data
        i_snr = snrdb_index[1]
        if (isGMI) then
            i_encoding = encoding_index[1]
        end if
        if (.not. isHard) then
            i_config = monoconf_index[1]
        end if

        ! Increment monoconf_index
        if (isReverse .and. .not. isHard) then
            monoconf_index[1] = mod(monoconf_index[1], size(monoConfig)) + 1
        end if

        ! Increment encoding_index
        if (isGMI) then
            if (isReverse .and. .not. isHard) then
                if (monoconf_index[1] == 1) then
                    encoding_index[1] = mod(encoding_index[1], size(encConfig, 1)) + 1
                end if
            else ! isHard .or. .not. isReverse
                encoding_index[1] = mod(encoding_index[1], size(encConfig, 1)) + 1
            end if
        end if

        ! Increment snrdb_index
        if (isGMI) then
            if (encoding_index[1] == 1) then
                if ( isReverse .and. .not. isHard ) then
                    if (monoconf_index[1] == 1) then
                        snrdb_index[1] = snrdb_index[1] + 1
                    end if
                else
                    snrdb_index[1] = snrdb_index[1] + 1
                end if
            end if
        else
            if ( isReverse .and. .not. isHard ) then
                if (monoconf_index[1] == 1) then
                    snrdb_index[1] = snrdb_index[1] + 1
                end if
            else
                snrdb_index[1] = snrdb_index[1] + 1
            end if
        end if

        ! Update output
        write(stdout, '("SNR: ", I0, "/", I0, " ")', advance='no') i_snr, nsnr
        if (isGMI) then
            write(stdout, '("ENC: ", I0, "/", I0, " ")', advance='no') i_encoding, size(encConfig, 1)
        end if
        if (.not. isHard) then
            write(stdout, '("CFG: ", I0, "/", I0, " ")', advance='no') i_config, size(monoConfig)
        end if
        write(stdout, "(A)"), ""
        unlock(lck[1])

        if (isGMI) then
            call noisemapper_set_encoding_custom(nm, encConfig(i_encoding, :))
        end if
        if (isReverse .and. .not. isHard) then
            call noisemapper_set_monotonicity(nm, &
                config_int_to_bool(nm, monoConfig(i_config)))
        end if

        ! Common setup
        call noisemapper_update_N0_from_snrdb(nm, snr_array(i_snr))
        call noisemapper_set_y_thresholds_uniform(nm) ! first guess
        opt_array(1:M_half-1) = nm%y_thresholds(M_half+1 : nm%M-1)

        if (isGMI) then
            opt_array(M_half) = 1d0
            if (isHard) then
                call lincoa_optimize_wrapper( calfun_gmi_hard_opt_threshold_s, opt_array )
            elseif (isReverse) then
                call lincoa_optimize_wrapper( calfun_gmi_soft_opt_threshold_s, opt_array )
            else
                NOT_IMPLEMENTED()
            end if
        else
            if (isHard) then
                call lincoa_optimize_wrapper( calfun_mi_hard_opt_threshold, nm%y_thresholds(M_half+1:nm%M-1) )
            elseif (isReverse) then
                call lincoa_optimize_wrapper( calfun_mi_soft_opt_threshold, nm%y_thresholds(M_half+1:nm%M-1) )
            else
                NOT_IMPLEMENTED()
            end if
        end if

        ! update result
        I = -I
        call update_result
        event post(snr_done(i_snr)[1])

        if (me==1 .and. n_im==1) then
            ! Query anyways so that we don't have to go throught all branches
            ! to know whether we completed the computation of the current SNR
            call event_query(snr_done(i_snr), ii)
            if (ii == until_count) then
                call write_result_to_file(i_snr)
            end if
        end if
    end do loop_snr

100 sync all

    ! if (me == 1) then
    !     write(stdout, *) ""
    !     outdata(:, 2) = outdata(:, 1) - 10*log10(outdata(:,3))
    !     if (isGMI) then
    !         output_root = trim(output_root)//"/opt-gmi"
    !     else
    !         output_root = trim(output_root)//"/opt-mi"
    !     end if

    !     call make_directory_and_file_name(output_root, bps, isReverse, isHard, &
    !         snr, nsnr, 0, 0, 0, 0, output_dir, output_name)
    !     call execute_command_line("mkdir -p " // trim(output_dir))

    !     ! open(newunit=io, file=trim(output_dir) // "/" // trim(output_name) // ".log", &
    !     !     status="replace", action="write")

    !     ! ! write(io, '(A, T16, A, T32, A, T48)') &
    !     ! !     "SNR [dB]", "Eb/N0 [dB]", "I"
    !     ! ! do i_snr = 1, nsnr
    !     ! !     write(io, '(f12.8, T16, f12.9, T32, E12.3E3)', advance='no') &
    !     ! !         outdata(i_snr, :3)
    !     ! !     do ii = 1, nm%M-1
    !     ! !         write(io, '(8X, E12.3E3)', advance='no') outdata(i_snr, o_mi+ii)
    !     ! !     end do
    !     ! !     if (isGMI) then
    !     ! !         write(io, '(8X, E12.3E3)') outdata(i_snr, o_s)
    !     ! !     else
    !     ! !         write(io, '(A)') "" ! Just add a newline
    !     ! !     end if
    !     ! ! end do
    !     ! write(io, *) header
    !     ! do i_snr = 1, nsnr
    !     !     write(io, *) outdata(i_snr, :)
    !     ! end do
    !     ! close(io)

    !     ! call to_file(x=outdata, file=trim(output_dir)//"/"//trim(output_name)//".csv", &
    !     !     header=header, fmt="e")
    ! end if


contains

    subroutine lincoa_optimize_wrapper( objective_function, initial_guess )
        interface
            subroutine f(x, y)
                double precision, intent(in) :: x(:)
                double precision, intent(out) :: y
            end subroutine f
        end interface
        procedure(f) :: objective_function
        double precision, intent(inout) :: initial_guess(:)

        call lincoa(&
            objective_function, initial_guess, &
            f=I, &
            Aineq=A, bineq=b, &
            rhobeg=2*nm%sigma, rhoend=1d-6)
    end subroutine lincoa_optimize_wrapper

    subroutine calfun_mi_hard_opt_threshold(theta, I_neg)
        !! Compute \( I(X;\hat{X}) \) with uniform input probabilities
        !! and varying thresholds
        double precision, intent(in) :: theta(:)
        !! Thresholds to be optimized [M/2 + 1, M-1]
        !! The negative thresholds are assumed to be symmetric
        double precision, intent(out) :: I_neg
        !! - I (as the lincoa routine looks for a minimum)

        integer :: i, j

        call noisemapper_set_y_thresholds(nm, [-theta(size(theta):1:-1), 0d0, theta(:)])
        call noisemapper_update_hard_reverse_tables(nm, .true.)

        I_neg = 0
        do j = 0, nm%M-1
            do i = 0, nm%M-1
                I_neg = I_neg + nm%fwd_probabilities(i, j) * nm%probabilities(i) * &
                    (log0(nm%delta_Fy(j)) - log0(nm%fwd_probabilities(i, j)))
            end do
        end do
        I_neg = I_neg/log(2d0)
    end subroutine calfun_mi_hard_opt_threshold


    subroutine calfun_mi_soft_opt_threshold(theta, I_neg)
        !! Compute \( I(X, N;\hat{X}) \) with uniform input probabilities
        !! and varying thresholds
        double precision, intent(in) :: theta(:)
        !! Thresholds to be optimized [M/2 + 1, M-1]
        !! The negative thresholds are assumed to be symmetric
        double precision, intent(out) :: I_neg
        !! - I (as the lincoa routine looks for a minimum)

        real(c_double) :: Abserr
        integer :: Neval, Ier, Limit, Lenw, Last

        integer :: Iwork(100)
        real(c_double) :: Work(400)
        Limit = 100
        Lenw = 400

        call noisemapper_set_y_thresholds(nm, [-theta(M_half-1:1:-1), 0d0, theta(:M_half-1)])
        call noisemapper_set_Fy_grids(nm)

        call dqags(f_soft_reverse, 0d0, 1d0, 1d-12, 1d-6, &
            I_neg, Abserr, Neval, Ier, &
            Limit, Lenw, Last, Iwork, Work)

        if (Ier /= 0) then
            print '("Error at ", f10.3, " [dB]: error ", i1)', snr_array(i_snr), Ier
        end if

        I_neg = - I_neg - H_Xhat(nm)
    end subroutine calfun_mi_soft_opt_threshold



    subroutine calfun_gmi_hard_opt_threshold_s(th_s, I_neg)
        !! Compute \( I_s(X;\hat{\mathbf{B}}) \) with uniform input probabilities
        !! and varying thresholds and s
        double precision, intent(in) :: th_s(:)
        !! Thresholds to be optimized [M/2 + 1, M-1] + s parameters
        !! The negative thresholds are assumed to be symmetric
        double precision, intent(out) :: I_neg
        !! - I (as the lincoa routine looks for a minimum)

        ! double precision, pointer :: ptheta(:)
        ! double precision, pointer :: ps

        associate( ptheta => th_s(1:M_half-1), ps => th_s(M_half) )
            call noisemapper_set_y_thresholds(nm, [-ptheta(size(ptheta):1:-1), 0d0, ptheta(:)])
            call noisemapper_update_hard_reverse_tables(nm, .true.)

            I_neg = -I_s_map_hard_reverse(q_map_hard_product, ps)
        end associate
    end subroutine calfun_gmi_hard_opt_threshold_s


    subroutine calfun_gmi_soft_opt_threshold_s(th_s, I_neg)
        !! Compute \( I_s(X;\hat{\mathbf{B}}) \) with uniform input probabilities
        !! and varying thresholds and s
        double precision, intent(in) :: th_s(:)
        !! Thresholds to be optimized [M/2 + 1, M-1] + s parameters
        !! The negative thresholds are assumed to be symmetric
        double precision, intent(out) :: I_neg
        !! - I (as the lincoa routine looks for a minimum)

        associate( ptheta => th_s(1:M_half-1), ps => th_s(M_half) )
            call noisemapper_set_y_thresholds(nm, [-ptheta(size(ptheta):1:-1), 0d0, ptheta(:)])
            call noisemapper_set_Fy_grids(nm)

            I_neg = -I_s_map_soft_reverse(q_map_soft_reverse_prod, s=ps)
        end associate
    end subroutine calfun_gmi_soft_opt_threshold_s


    function config_int_to_bool(nm, c) result(b)
        type(noisemapper_type), intent(in) :: nm
        integer, intent(in) :: c
        logical(1) :: b(0:nm%bps-1)

        integer :: i

        do i = 0, nm%bps-1
            b(i) = iand(ishft(c, -i), 1) == 1
        end do
    end function config_int_to_bool


    subroutine write_header(tmpstr, start, end, header)
        !! print string to header
        character(20), value, intent(in) :: tmpstr
        integer, intent(in) :: start
        integer, intent(in) :: end
        character(20), intent(inout) :: header(start:end)

        ! character(20) :: tpmcopy
        integer :: i
        do i = end, start, -1
            write(header(i), '(A, "_", I0)') trim(tmpstr), i
        end do

    end subroutine write_header


    subroutine update_result
        critical
            if (I .gt. outdata(i_snr, o_mi)[1]) then
                outdata(i_snr, o_mi)[1] = I
                if (isReverse) then
                    outdata(i_snr, o_th:o_p-1)[1] = nm%y_thresholds
                    outdata(i_snr, o_p:o_s-1)[1]  = nm%delta_Fy
                    if (.not. isHard) then
                        outdata(i_snr, o_c)[1]    = monoConfig(i_config)
                    end if
                end if
                if (isGMI) then
                    outdata(i_snr, o_e:o_c-1)[1]  = encConfig(i_encoding, :)
                    outdata(i_snr, o_s)[1]        = opt_array(M_half)
                end if
            end if
        end critical
    end subroutine update_result


    subroutine write_result_to_file(snr_i)
        integer, intent(in) :: snr_i

        outdata(snr_i, 2) = outdata(snr_i, 1) - 10*log10(outdata(snr_i, 3))
        if (isGMI) then
            if (isReverse) then
                if (isHard) then
                    write(io_log, format_log) outdata(snr_i, :o_s), int(outdata(snr_i, o_e:o_c-1))
                    write(io_csv, format_csv) outdata(snr_i, :o_s), int(outdata(snr_i, o_e:o_c-1))
                else
                    write(io_log, format_log) outdata(snr_i, :o_s), int(outdata(snr_i, o_e:o_c))
                    write(io_csv, format_csv) outdata(snr_i, :o_s), int(outdata(snr_i, o_e:o_c))
                end if
            else
                NOT_IMPLEMENTED()
            end if
        else
            if (isReverse) then
                if (isHard) then
                    write(io_log, format_log) outdata(snr_i, :o_s-1)
                    write(io_csv, format_csv) outdata(snr_i, :o_s-1)
                else
                    write(io_log, format_log) outdata(snr_i, :o_s-1), int(outdata(snr_i, o_c))
                    write(io_csv, format_csv) outdata(snr_i, :o_s-1), int(outdata(snr_i, o_c))
                end if
            else
                NOT_IMPLEMENTED()
            end if
        end if
        flush(io_log)
        flush(io_csv)
    end subroutine write_result_to_file

end program mi_opt
