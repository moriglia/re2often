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
module re2often_noisemapper
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Implementation of alphabet utilities and noise mapping/demapping
    !! Mainteining C interoperability
    use, intrinsic :: iso_c_binding
    use re2often_utils, only: binsearch
    use stdlib_stats_distribution_normal, only: cdf_normal
    implicit none

    type :: noisemapper_type
        !! Descriptor for the alphabet and the noise channel
        integer(c_int) :: bps
        !! Bit per symbol
        integer(c_int) :: M
        !! order of modulation (number of constellation symbols)
        real(c_double), allocatable :: constellation(:)
        !! constellation points
        real(c_double), allocatable :: probabilities(:)
        !! probabilities for each constellation point
        logical(c_bool), allocatable :: s_to_b(:,:)
        !! c_bool pointer to symbol to bit map
        real(c_double) :: E_s
        !! Expected symbol energy (per quadrature, only dealing with PAM)
        real(c_double) :: N0
        !! Expected noise variance (on both quadratures)
        real(c_double) :: sigma
        !! standard deviation of noise (only one quadrature)

        ! +-----------------------------+
        ! | Reverse reconciliation data |
        ! +-----------------------------+
        real(c_double), allocatable :: y_thresholds(:)
        !! Decision thresholds. It ranges from `1` to `M-1`, with
        !! `1` corresponding to the threshold between symbol `0` and
        !! symbol `1`
        real(c_double), allocatable :: Fy_thresholds(:)
        !! Cumulative Density Function of the channel output at the
        !! thresholds. It ranges from `0` to `M`, with `Fy_thresholds(0) = 0`,
        !! `Fy_thresholds(M) = 1`, else `Fy_thresholds(i)` is the CDF
        !! evaluated at `y_thresholds(i)`
        real(c_double), allocatable :: delta_Fy(:)
        !! Probability that the channel output lays in the
        !! decision region of each symbol

        ! +----------------------------------+
        ! | Hard reverse reconciliation data |
        ! +----------------------------------+
        real(c_double), allocatable :: fwd_probabilities(:,:)
        !! Forward transition probabilities (likelihoods):
        !! Location (i, j) contains \(P(\hat{X}=a_j | X=a_i)\)
        real(c_double), allocatable :: reverse_hard_lappr_table(:,:)
        !! table of LAPPRs of the received bits given a transmitted symbol.
        !! Location (i, k) contains the LAPPR(k) given \(X=a_i\).
        !! `i` ranges in `(0, M)`, `k` ranges in `(0, bps)`

        ! real(c_double), allocatable :: bwd_probabilities(:,:)
        ! !! Backward transition probabilities (a posteriori probabilities):
        ! !! Location (i, j) contains \(P(X=a_i | \hat{X}=a_j)\)
        ! !! Mind the inversion of indexes with respect to `fwd_probabilities`

        ! +---------------------------------------+
        ! | Reverse Reconciliation SOFTENING data |
        ! +---------------------------------------+
        logical(c_bool), allocatable :: monotonicity_configuration(:)
        !! Monotonicity configuration `(0:M-1)`: `.false.` means increasing,
        !! `.true.` means decreasing.
        real(c_double), allocatable :: Fy_grid(:)
        !! grid of CDF values taken at equally spaced intervals.
        !! Note that it is 1-based
        real(c_double) :: base_y_grid
        !! First element of the y grid
        real(c_double) :: y_grid_step
        !! step of the y grid
    end type noisemapper_type

    ! +--------------------------------------+
    ! | Interfaces for DIRECT reconciliation |
    ! +--------------------------------------+
    interface noisemapper_y_to_lappr
        module procedure noisemapper_y_to_lappr_single
        module procedure noisemapper_y_to_lappr_array
    end interface noisemapper_y_to_lappr

    interface noisemapper_random_symbol
        module procedure noisemapper_random_symbol_single
        module procedure noisemapper_random_symbol_array
    end interface noisemapper_random_symbol

    ! +---------------------------------------+
    ! | Interfaces for REVERSE reconciliation |
    ! +---------------------------------------+
    interface noisemapper_decide_symbol
        module procedure noisemapper_decide_symbol_single
        module procedure noisemapper_decide_symbol_array
    end interface noisemapper_decide_symbol

    ! +--------------------------------------------+
    ! | Interfaces for SOFT REVERSE reconciliation |
    ! +--------------------------------------------+
    interface noisemapper_generate_soft_metric
        module procedure noisemapper_generate_soft_metric_single
        module procedure noisemapper_generate_soft_metric_array
    end interface noisemapper_generate_soft_metric

    interface noisemapper_soft_reverse_lappr
        module procedure noisemapper_soft_reverse_lappr_single
        module procedure noisemapper_soft_reverse_lappr_array
    end interface noisemapper_soft_reverse_lappr

contains


    module subroutine noisemapper_deallocate(nm)
        !! Destructor for noise mapper
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        if (allocated(nm%constellation)) deallocate(nm%constellation)
        if (allocated(nm%probabilities)) deallocate(nm%probabilities)
        if (allocated(nm%s_to_b       )) deallocate(nm%s_to_b       )
    end subroutine noisemapper_deallocate


    module function noisemapper_create(bps) result(nm)
        !! Create the nm object
        integer(c_int), intent(in) :: bps
        !! bit per symbol
        type(noisemapper_type) :: nm
        !! Noise mapper
        !! This function will not setup N0

        integer :: i, k

        call noisemapper_deallocate(nm)

        nm%bps = bps
        nm%M   = ishft(1, bps)

        allocate(nm%constellation(0:nm%M-1))
        allocate(nm%probabilities(0:nm%M-1))
        nm%constellation = [(real(1-nm%M, c_double) + real(2*i, c_double), &
            i = 0, nm%M-1) ]
        nm%probabilities = 1d0/real(nm%M, c_double)

        nm%E_s = sum(nm%probabilities * abs(nm%constellation)**2)

        allocate(nm%s_to_b(0:nm%M-1 , 0:nm%bps-1))
        do i = 0, nm%M-1
            do k = 0, nm%bps - 1
                nm%s_to_b(i, k) = iand(ishft(ishft(i, -k)+1, -1), 1)==1
            end do
        end do
    end function noisemapper_create


    module subroutine noisemapper_update_N0_from_snrdb(nm, snrdb)
        !! Update N0 based on the value of the SNR
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: snrdb
        !! SNR in dB

        nm%N0 = nm%E_s * (10d0**(-snrdb/10d0))
        nm%sigma = sqrt(nm%N0/2d0)
    end subroutine noisemapper_update_N0_from_snrdb


    module subroutine noisemapper_y_to_lappr_single(nm, y, lappr)
        !! calculate lappr from channel output sample for direct reconciliation
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y
        !! AWGN channel output sample
        real(c_double), intent(out) :: lappr(0:nm%bps - 1)
        !! log-a posteriori-probabilities for the de-mapped bits

        real(c_double) :: den(0:nm%bps-1)
        real(c_double) :: addendum
        integer :: i, k

        den(:)   = 0
        lappr(:) = 0

        do i = 0, nm%M-1
            addendum = nm%probabilities(i) * exp(-(y - nm%constellation(i))**2/nm%N0)
            do k = 0, nm%bps - 1
                if (nm%s_to_b(i, k)) then
                    den(k) = den(k) + addendum
                else
                    lappr(k) = lappr(k) + addendum
                end if
            end do
        end do
        lappr = log(lappr) - log(den)
    end subroutine noisemapper_y_to_lappr_single


    module subroutine noisemapper_y_to_lappr_array(nm, y, lappr)
        !! calculate lappr from set of channel output samples for direct reconciliation
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y(0:)
        !! AWGN channel samples
        real(c_double), intent(out) :: lappr(0:size(y)*nm%bps-1)
        !! LAPPR corresponding to the channel outputs

        integer :: j

        do j = 0, size(y)-1
            call noisemapper_y_to_lappr_single(nm, y(j), lappr(j*nm%bps : (j+1)*nm%bps - 1))
        end do
    end subroutine noisemapper_y_to_lappr_array


    module subroutine noisemapper_random_symbol_single(nm, x_i)
        !! Generate a random symbol
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(out) :: x_i
        !! Random symbol of the constellation (index in 0:M-1)

        double precision :: rnd
        integer :: i

        call random_number(rnd) ! rnd is in [0, 1)

        x_i = nm%M-1
        do i = 0, nm%M-2
            if (rnd .lt. nm%probabilities(i)) then
                x_i = i
                return
            end if
            rnd = rnd - nm%probabilities(i)
        end do
    end subroutine noisemapper_random_symbol_single


    module subroutine noisemapper_random_symbol_array(nm, x_i)
        !! Generate random symbols
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(out) :: x_i(:)
        !! Random symbols of the constellation (index in 0:M-1)

        integer :: j
        do j = 1, size(x_i)
            call noisemapper_random_symbol_single(nm, x_i(j))
        end do
    end subroutine noisemapper_random_symbol_array


    module function noisemapper_symbol_index_to_value(nm, x_i) result (x)
        !! Convert constellation index to point
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(in) :: x_i(:)
        !! Set of constellation indexes
        real(c_double) :: x(size(x_i))
        !! set of constellation points

        integer :: j

        do j = 1, size(x_i)
            x(j) = nm%constellation(x_i(j))
        end do
    end function noisemapper_symbol_index_to_value


    module function noisemapper_symbol_to_word(nm, x_i) result (word)
        !! Convert a set of symbol indexes to a word
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(in) :: x_i(0:)
        !! Set of indexes
        logical(c_bool) :: word(0:size(x_i)*nm%bps - 1)
        !! Word corresponding to the sequence of symbols

        integer :: j

        do j = 0, size(x_i) - 1
            word(j*nm%bps : (j+1)*nm%bps-1) = nm%s_to_b(x_i(j), :)
        end do
    end function noisemapper_symbol_to_word

    ! +-----------------------------------------+
    ! | Common reverse reconciliation functions |
    ! +-----------------------------------------+

    module subroutine noisemapper_deallocate_reverse_common(nm)
        !! Deallocation of common arrays for reverse reconciliation
        type(noisemapper_type), intent(inout) :: nm
        !! noise mapper

        if (allocated(nm%y_thresholds) ) deallocate(nm%y_thresholds)
        if (allocated(nm%Fy_thresholds)) deallocate(nm%Fy_thresholds)
        if (allocated(nm%delta_Fy)     ) deallocate(nm%delta_Fy)
    end subroutine noisemapper_deallocate_reverse_common


    module subroutine noisemapper_allocate_reverse_common(nm)
        !! Allocate the common reverse reconciliation buffers
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        call noisemapper_deallocate_reverse_common(nm)

        allocate(nm%y_thresholds(1:nm%M-1))
        allocate(nm%Fy_thresholds(0:nm%M))
        allocate(nm%delta_Fy(0:nm%M-1))
    end subroutine noisemapper_allocate_reverse_common


    module subroutine noisemapper_set_y_thresholds(nm, thresholds)
        !! Set the decision thresholds
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        real(c_double), intent(in), optional :: thresholds(1:nm%M-1)
        !! y thresholds. If not present, the thresholds are the points
        !! half way between two adjacent constellation points
        !! @warning The order of the array is not checked

        integer :: i


        call noisemapper_allocate_reverse_common(nm)

        ! Set the thresholds
        if (.not. present(thresholds)) then
            nm%y_thresholds = (nm%constellation(0:nm%M-2) + nm%constellation(1:nm%M-1))/2d0
        else
            nm%y_thresholds = thresholds
        end if

        ! Set the CDF at each threshold
        nm%Fy_thresholds(0)    = 0 ! at -\infty
        nm%Fy_thresholds(nm%M) = 1 ! at +\infty
        do i = 1, nm%M-1
            nm%Fy_thresholds(i) = sum( nm%probabilities * cdf_normal(&
                x    = nm%y_thresholds(i), &
                loc  = nm%constellation,   &
                scale= nm%sigma))  ! at y_threshold(i)
        end do

        ! Set the probability of the channel output being in each decision region
        nm%delta_Fy = nm%Fy_thresholds(1:nm%M) - nm%Fy_thresholds(0:nm%M-1)
    end subroutine noisemapper_set_y_thresholds


    module function noisemapper_decide_symbol_single(nm, y) result(x_i)
        !! Take a decision for the received channel output based
        !! on the thresholds
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y
        !! Channel output sample
        integer(c_int) :: x_i
        !! Alphabet index of the decided symbol

        x_i = binsearch(nm%y_thresholds, y)

    end function noisemapper_decide_symbol_single


    module function noisemapper_decide_symbol_array(nm, y) result(x_i)
        !! Take a decision for the set of received channel outputs
        !! based on thresholds
        type(noisemapper_type), intent(in) :: nm
        !! Noisemapper
        real(c_double), intent(in) :: y(:)
        !! Set of input samples
        integer(c_int) :: x_i(size(y))
        !! Decisions

        integer :: i

        do i = 1, size(y)
            x_i(i) = binsearch(nm%y_thresholds, y(i))
        end do
    end function noisemapper_decide_symbol_array

    ! +---------------------------------------+
    ! | Hard reverse reconciliation functions |
    ! +---------------------------------------+

    module subroutine noisemapper_deallocate_reverse_hard(nm)
        !! Deallocate transition probability table and lappr table
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        if (allocated(nm%fwd_probabilities)) deallocate(nm%fwd_probabilities)
        if (allocated(nm%reverse_hard_lappr_table)) deallocate(nm%reverse_hard_lappr_table)
    end subroutine noisemapper_deallocate_reverse_hard


    module subroutine noisemapper_allocate_reverse_hard(nm)
        !! Allocate transition probability table and lappr table
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        call noisemapper_deallocate_reverse_hard(nm)

        allocate(nm%fwd_probabilities(0:nm%M-1, 0:nm%M-1))
        allocate(nm%reverse_hard_lappr_table(0:nm%M-1, 0:nm%bps-1))
    end subroutine noisemapper_allocate_reverse_hard


    module subroutine noisemapper_update_hard_reverse_tables(nm)
        !! Update hard reverse reconciliation tables
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        integer :: i, j, k
        real(c_double) :: denominator(0:nm%M-1, 0:nm%bps-1)

        call noisemapper_allocate_reverse_hard(nm)

        do i = 0, nm%M - 1
            do j = 0, nm%M-2
                nm%fwd_probabilities(i,j) = cdf_normal(&
                    x     = nm%y_thresholds(j+1),      &
                    loc   = nm%constellation(i),       &
                    scale = nm%sigma                   )
            end do
        end do
        nm%fwd_probabilities(:, nm%M-1) = 1
        nm%fwd_probabilities(:, 1:nm%M-1) = nm%fwd_probabilities(:, 1:nm%M-1) &
            - nm%fwd_probabilities(:, 0:nm%M-2)

        denominator(:,:) = 0
        nm%reverse_hard_lappr_table(:,:) = 0
        ! do i = 0, nm%M-1 ! transmitted symbol
        do j = 0, nm%M-1 ! received symbol
            do k = 0, nm%bps-1 ! received bit
                if (nm%s_to_b(j, k)) then
                    denominator(:, k) = denominator(:, k) + nm%fwd_probabilities(:, j)
                else
                    nm%reverse_hard_lappr_table(:, k) = &
                        nm%reverse_hard_lappr_table(:, k) + nm%fwd_probabilities(:, j)
                end if
            end do
        end do
        ! end do

        nm%reverse_hard_lappr_table = log(nm%reverse_hard_lappr_table) - log(denominator)
    end subroutine noisemapper_update_hard_reverse_tables


    module subroutine noisemapper_convert_symbol_to_hard_lappr(nm, x_i, lappr)
        !! Get LAPPR for each symbol from the tables
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(in) :: x_i(0:)
        !! Transmitted symbols
        real(c_double), intent(out) :: lappr(0:nm%bps*size(x_i)-1)
        !! LAPPR array associated with the transmitted sybmols

        integer :: i

        do i = 0, size(x_i) - 1
            lappr(i * nm%bps : (i+1) * nm%bps - 1) = nm%reverse_hard_lappr_table(x_i(i), :)
        end do
    end subroutine noisemapper_convert_symbol_to_hard_lappr


    ! +---------------------------------------------+
    ! | REVERSE RECONCILIATION SOFTENING procedures |
    ! +---------------------------------------------+
    elemental module function noisemapper_Fy(nm, y) result(Fy)
        !! Evaluate the CDF of the output at the given point
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y
        !! Channel output sample
        real(c_double) :: Fy
        !! CDF of the channel output at `y`

        Fy = sum(nm%probabilities * cdf_normal(x=y, loc=nm%constellation, scale=nm%sigma))
    end function noisemapper_Fy


    module subroutine noisemapper_deallocate_reverse_soft(nm)
        !! Deallocate data used for soft reverse reconciliation
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        if (allocated(nm%monotonicity_configuration)) deallocate(nm%monotonicity_configuration)
        if (allocated(nm%Fy_grid)) deallocate(nm%Fy_grid)
    end subroutine noisemapper_deallocate_reverse_soft


    module subroutine noisemapper_set_monotonicity(nm, config)
        !! Set monotonicity configuration
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        logical(c_bool), optional, intent(in) :: config(0:nm%M-1)
        !! Configuration

        if (allocated(nm%monotonicity_configuration)) then
            if (size(nm%monotonicity_configuration) /= nm%M) then
                deallocate(nm%monotonicity_configuration)
            end if
        end if
        if (.not. allocated(nm%monotonicity_configuration)) then
            allocate(nm%monotonicity_configuration(0:nm%M-1))
        end if
        if (present(config)) then
            nm%monotonicity_configuration = config
        else
            nm%monotonicity_configuration(0::2) = .false.
            nm%monotonicity_configuration(1::2) = .true.
        end if
    end subroutine noisemapper_set_monotonicity


    module subroutine noisemapper_generate_soft_metric_single(nm, y, n, xhat)
        !! Generate soft metric from a single channel output sample
        !! and give the decided symbol, too.
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y
        !! Channel output sample
        real(c_double), intent(out) :: n
        !! Soft metric
        integer(c_int), intent(out) :: xhat
        !! decided symbol

        xhat = noisemapper_decide_symbol_single(nm, y)

        n = (noisemapper_Fy(nm, y)- nm%Fy_thresholds(xhat))/nm%delta_Fy(xhat)
        if (nm%monotonicity_configuration(xhat)) then
            n = 1d0 - n
        end if
    end subroutine noisemapper_generate_soft_metric_single


    module subroutine noisemapper_generate_soft_metric_array(nm, y, n, xhat)
        !! Generate soft metric from a set of channel output samples
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y(:)
        !! Channel output sample
        real(c_double), intent(out) :: n(size(y))
        !! Soft metric
        integer(c_int), intent(out) :: xhat(size(y))
        !! decided symbol

        integer :: i

        do i = 1, size(y)
            call noisemapper_generate_soft_metric_single(nm, y(i), n(i), xhat(i))
        end do
    end subroutine noisemapper_generate_soft_metric_array


    module subroutine noisemapper_set_Fy_grids(nm, th)
        !! Setup the grid for the inverse of the CDF
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        real(c_double), intent(in), optional :: th
        !! threshold for the PDF minimum value.
        !! It must be strictly positive, lower than 1

        real(c_double) :: threshold
        real(c_double) :: y_start, y_stop
        integer :: n_points, i

        if (present(th)) then
            threshold = th
        else
            threshold = 1d-9
        end if

        y_stop = nm%constellation(nm%M-1) + sqrt(-2*(nm%sigma**2)*log(threshold))
        y_start = -y_stop

        nm%base_y_grid = y_start
        nm%y_grid_step = 1d-3
        n_points = ceiling(2*y_stop/nm%y_grid_step)

        if (allocated(nm%Fy_grid)) then
            deallocate(nm%Fy_grid)
        end if
        allocate(nm%Fy_grid(n_points))

        nm%Fy_grid = noisemapper_Fy(nm, &
            [(nm%base_y_grid + i*nm%y_grid_step, i=1, n_points)])
    end subroutine noisemapper_set_Fy_grids


    module function noisemapper_invert_soft_metric(nm, n, x_i) result(y)
        !! generate all tentative channel output samples from the received soft metric
        type(noisemapper_type), intent(in) :: nm
        !! noise mapper
        real(c_double), intent(in) :: n
        !! Soft metric
        integer(c_int), intent(in) :: x_i
        !! Hypotetical received symbol alphabet index
        real(c_double) :: y
        !! Tentative channel output samples

        real(c_double) :: Fy
        integer :: idx

        Fy = n
        if (nm%monotonicity_configuration(x_i)) then
            Fy = 1-Fy
        end if
        Fy = nm%Fy_thresholds(x_i) + Fy*nm%delta_Fy(x_i)

        idx = binsearch(nm%Fy_grid, Fy)

        y = nm%base_y_grid + idx*nm%y_grid_step
    end function noisemapper_invert_soft_metric


    module subroutine noisemapper_soft_reverse_lappr_single(nm, x_i, n, lappr)
        !! Calculate the LAPPR from the transmitted symbol and the soft metric
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(in) :: x_i
        !! transmitted symbol alphabet index
        real(c_double), intent(in) :: n
        !! Soft metric
        real(c_double), intent(out) :: lappr(0:nm%bps-1)
        !! LAPPRs

        integer :: i !! alphabet index of tentative received symbol
        integer :: k !! further alphabet/bit index
        real(c_double) :: denominator(0:nm%bps-1)
        real(c_double) :: twoy, x, addendum


        lappr(:) = 0
        denominator(:) = 0

        x = nm%constellation(x_i)

        do i = 0, nm%M-1
            twoy = 2*noisemapper_invert_soft_metric(nm, n, i)

            addendum = 0
            do k = 0, nm%M-1 ! k used as symbol index
                addendum = addendum + nm%probabilities(k) &
                    * exp((twoy - x - nm%constellation(k))*(nm%constellation(k) - x)/nm%N0)
            end do
            addendum = nm%delta_Fy(i) / addendum

            do k = 0, nm%bps - 1 ! k used as bit index
                if (nm%s_to_b(i, k)) then
                    denominator(k) = denominator(k) + addendum
                else
                    lappr(k) = lappr(k) + addendum
                end if
            end do
        end do


        lappr = log(lappr) - log(denominator)
    end subroutine noisemapper_soft_reverse_lappr_single


    module subroutine noisemapper_soft_reverse_lappr_array(nm, x_i, n, lappr)
        !! Calculate the LAPPR from the transmitted symbols and the soft metrics arrays
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(in) :: x_i(0:)
        !! transmitted symbol alphabet index
        real(c_double), intent(in) :: n(0:size(x_i)-1)
        !! Soft metric
        real(c_double), intent(out) :: lappr(0:size(x_i)*nm%bps-1)
        !! LAPPRs

        integer :: i

        do i = 0, size(x_i)-1
            call noisemapper_soft_reverse_lappr_single(nm, x_i(i), n(i), lappr(i*nm%bps : (i+1)*nm%bps-1))
        end do
    end subroutine noisemapper_soft_reverse_lappr_array


    ! +-----------------------------+
    ! | Uniform probabilities at RX |
    ! +-----------------------------+
    module subroutine noisemapper_set_y_thresholds_uniform(nm)
        !! Set thresholds with uniform decision probabilities
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        integer :: i
        real(c_double) :: thresholds(1:nm%M-1)

        do i = 1, nm%M-1
            thresholds(i) = nm%base_y_grid + &
                nm%y_grid_step * binsearch(nm%Fy_grid, real(i, c_double)/real(nm%M, c_double))
        end do
        call noisemapper_set_y_thresholds(nm, thresholds)
    end subroutine noisemapper_set_y_thresholds_uniform
end module re2often_noisemapper
