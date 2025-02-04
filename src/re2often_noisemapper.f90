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
end module re2often_noisemapper
