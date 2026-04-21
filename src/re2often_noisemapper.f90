! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025-2026  Marco Origlia

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
submodule (re2often) re2often_noisemapper
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Implementation of alphabet utilities and noise mapping/demapping
    !! Mainteining C interoperability
    use, intrinsic :: iso_c_binding
    use re2often_utils, only: binsearch
    use stdlib_stats_distribution_normal, only: cdf_normal
    implicit none

contains
    module subroutine noisemapper_finalize(nm)
        !! Wrapper for the deallocate function
        type(noisemapper_type) :: nm
        !! Noisemapper

        call nm%deallocate
    end subroutine noisemapper_finalize

    module subroutine noisemapper_deallocate(nm)
        !! Destructor for noise mapper
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        if (allocated(nm%constellation)) deallocate(nm%constellation)
        if (allocated(nm%probabilities)) deallocate(nm%probabilities)
        if (allocated(nm%s_to_b       )) deallocate(nm%s_to_b       )

        call nm%deallocate_reverse_hard
        call nm%deallocate_reverse_soft
        call nm%deallocate_reverse_common

        if (allocated(nm%alice_bit_priors)) deallocate(nm%alice_bit_priors)
    end subroutine noisemapper_deallocate


    module subroutine noisemapper_set_symbol_probabilities(nm, probabilities)
        !! Allocate and set probability vector for imput constellation symbols
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        real(c_double), intent(in), optional :: probabilities(0:nm%M-1)
        !! Input probabilities

        if (allocated(nm%probabilities) .and. (size(nm%probabilities)/=nm%M)) then
            deallocate(nm%probabilities)
        end if
        if (.not. allocated(nm%probabilities)) then
            allocate(nm%probabilities(0:nm%M-1))
        end if

        if (.not. present(probabilities)) then
            nm%probabilities = 1d0/real(nm%M, c_double)
        else
            nm%probabilities = probabilities
        end if
        nm%E_s = sum(nm%probabilities * abs(nm%constellation)**2)

        call nm%update_alice_priors
    end subroutine noisemapper_set_symbol_probabilities


    module function noisemapper_create(bps) result(nm)
        !! Create the nm object
        integer(c_int), intent(in) :: bps
        !! bit per symbol
        type(noisemapper_type) :: nm
        !! Noise mapper
        !! This function will not setup N0

        integer :: i

        call nm%deallocate()

        nm%bps = bps
        nm%M   = ishft(1, bps)

        allocate(nm%constellation(0:nm%M-1))
        nm%constellation = [(real(1-nm%M, c_double) + real(2*i, c_double), &
            i = 0, nm%M-1) ]

        call nm%set_symbol_probabilities()

        allocate(nm%s_to_b(0:nm%M-1 , 0:nm%bps-1))
        call nm%set_encoding_gray()
    end function noisemapper_create


    module subroutine noisemapper_set_encoding_gray(nm)
        !! Set Gray encoding
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        integer :: i, k
        do i = 0, nm%M-1
            do k = 0, nm%bps - 1
                nm%s_to_b(i, k) = iand(ishft(ishft(i, -k)+1, -1), 1)==1
            end do
        end do

        call nm%update_alice_priors
    end subroutine noisemapper_set_encoding_gray


    module subroutine noisemapper_set_encoding_natural(nm)
        !! Set Gray encoding
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        integer :: i, k
        do i = 0, nm%M-1
            do k = 0, nm%bps - 1
                nm%s_to_b(i, k) = iand(ishft(i, -k), 1)==1
            end do
        end do

        call nm%update_alice_priors
    end subroutine noisemapper_set_encoding_natural


    module subroutine noisemapper_set_encoding_custom(nm, labels)
        !! Set custom encoding labels
        !! @warning: no check on labels being a complete permutation of (0:M-1)
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        integer, intent(in) :: labels(0:nm%M-1)
        !! Encoding Labels

        integer :: i, k
        do i = 0, nm%M-1
            do k = 0, nm%bps-1
                nm%s_to_b(i, k) = iand(ishft(labels(i), -k), 1)==1
            end do
        end do

        call nm%update_alice_priors
    end subroutine noisemapper_set_encoding_custom


    module subroutine noisemapper_update_N0_from_snrdb(nm, snrdb)
        !! Update N0 based on the value of the SNR
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: snrdb
        !! SNR in dB

        nm%N0 = nm%E_s * (10d0**(-snrdb/10d0))
        nm%sigma = sqrt(nm%N0/2d0)
    end subroutine noisemapper_update_N0_from_snrdb


    module subroutine noisemapper_y_to_lappr_single(nm, y, lappr)
        !! calculate lappr from channel output sample for direct reconciliation
        class(noisemapper_type), intent(in) :: nm
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

        do k = 0, nm%bps - 1
            if (den(k) == 0) then
                lappr(k) = 1d100
            elseif (lappr(k) == 0) then
                lappr(k) = -1d100
            else
                lappr(k) = log(lappr(k)) - log(den(k))
            end if
        end do
    end subroutine noisemapper_y_to_lappr_single


    module subroutine noisemapper_y_to_lappr_array(nm, y, lappr)
        !! calculate lappr from set of channel output samples for direct reconciliation
        class(noisemapper_type), intent(in) :: nm
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


    module subroutine noisemapper_update_alice_priors(nm)
        !! Update a priori probabilities for Alice
        !! This function returns with no error nor warning
        !! if either the encoding or the probabilities aren't set
        class(noisemapper_type), intent(inout) :: nm
        !! Noisemapper

        integer :: l
        real(c_double) :: denominator

        if (.not. allocated(nm%s_to_b))        return
        if (.not. allocated(nm%probabilities)) return

        if (allocated(nm%alice_bit_priors)) deallocate(nm%alice_bit_priors)
        allocate(nm%alice_bit_priors(0:nm%bps-1))

        do l = 0, nm%bps-1
            denominator = sum(merge(nm%probabilities, 0d0, nm%s_to_b(:,l)))
            if (denominator .le. 0) then
                nm%alice_bit_priors(l) = 1d100
                cycle
            end if
            nm%alice_bit_priors(l) = 1 - denominator
            if (nm%alice_bit_priors(l) .le. 0) then
                nm%alice_bit_priors(l) = -1d100
                cycle
            end if
            nm%alice_bit_priors(l) = log(nm%alice_bit_priors(l)) - log(denominator)
        end do
    end subroutine noisemapper_update_alice_priors


    module subroutine noisemapper_y_to_llr_single(nm, y, llr)
        !! Calculate LLR from the channel output
        class(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y
        !! AWGN channel output sample
        real(c_double), intent(out) :: llr(0:nm%bps-1)
        !! log-likelyhood ratios of the bits associated to one single symbol transmission

        call noisemapper_y_to_lappr_single(nm, y, llr)
        llr = llr - nm%alice_bit_priors
    end subroutine noisemapper_y_to_llr_single


    module subroutine noisemapper_y_to_llr_array(nm, y, llr)
        !! calculate lappr from set of channel output samples for direct reconciliation
        class(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y(0:)
        !! AWGN channel samples
        real(c_double), intent(out) :: llr(0:size(y)*nm%bps-1)
        !! LAPPR corresponding to the channel outputs

        integer :: j

        do j = 0, size(y)-1
            call noisemapper_y_to_llr_single(nm, y(j), llr(j*nm%bps : (j+1)*nm%bps - 1))
        end do
    end subroutine noisemapper_y_to_llr_array


    module subroutine noisemapper_random_symbol_single(nm, x_i)
        !! Generate a random symbol
        class(noisemapper_type), intent(in) :: nm
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
        class(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(out) :: x_i(:)
        !! Random symbols of the constellation (index in 0:M-1)

        integer :: j
        do j = 1, size(x_i)
            call nm%random_symbol(x_i(j))
        end do
    end subroutine noisemapper_random_symbol_array


    module function noisemapper_symbol_index_to_value(nm, x_i) result (x)
        !! Convert constellation index to point
        class(noisemapper_type), intent(in) :: nm
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
        class(noisemapper_type), intent(in) :: nm
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
        class(noisemapper_type), intent(inout) :: nm
        !! noise mapper

        if (allocated(nm%y_thresholds) ) deallocate(nm%y_thresholds)
        if (allocated(nm%Fy_thresholds)) deallocate(nm%Fy_thresholds)
        if (allocated(nm%delta_Fy)     ) deallocate(nm%delta_Fy)
    end subroutine noisemapper_deallocate_reverse_common


    module subroutine noisemapper_allocate_reverse_common(nm)
        !! Allocate the common reverse reconciliation buffers
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        call nm%deallocate_reverse_common()

        allocate(nm%y_thresholds(1:nm%M-1))
        allocate(nm%Fy_thresholds(0:nm%M))
        allocate(nm%delta_Fy(0:nm%M-1))
    end subroutine noisemapper_allocate_reverse_common


    module subroutine noisemapper_set_y_thresholds(nm, thresholds)
        !! Set the decision thresholds
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        real(c_double), intent(in), optional :: thresholds(1:nm%M-1)
        !! y thresholds. If not present, the thresholds are the points
        !! half way between two adjacent constellation points
        !! @warning The order of the array is not checked

        integer :: i


        call nm%allocate_reverse_common()

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


    impure elemental module function noisemapper_decide_symbol_single(nm, y) result(x_i)
        !! Take a decision for the received channel output based
        !! on the thresholds
        class(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y
        !! Channel output sample
        integer(c_int) :: x_i
        !! Alphabet index of the decided symbol

        x_i = binsearch(nm%y_thresholds, y)
    end function noisemapper_decide_symbol_single


    ! module function noisemapper_decide_symbol_array(nm, y) result(x_i)
    !     !! Take a decision for the set of received channel outputs
    !     !! based on thresholds
    !     class(noisemapper_type), intent(in) :: nm
    !     !! Noisemapper
    !     real(c_double), intent(in) :: y(:)
    !     !! Set of input samples
    !     integer(c_int) :: x_i(size(y))
    !     !! Decisions

    !     integer :: i

    !     do i = 1, size(y)
    !         x_i(i) = binsearch(nm%y_thresholds, y(i))
    !     end do
    ! end function noisemapper_decide_symbol_array

    ! +---------------------------------------+
    ! | Hard reverse reconciliation functions |
    ! +---------------------------------------+

    module subroutine noisemapper_deallocate_reverse_hard(nm)
        !! Deallocate transition probability table and lappr table
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        if (allocated(nm%fwd_probabilities)) deallocate(nm%fwd_probabilities)
        if (allocated(nm%reverse_hard_lappr_table)) deallocate(nm%reverse_hard_lappr_table)
    end subroutine noisemapper_deallocate_reverse_hard


    module subroutine noisemapper_allocate_reverse_hard(nm, skipLapprTable)
        !! Allocate transition probability table and lappr table
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        logical, intent(in), optional :: skipLapprTable
        !! Do not allocate LAPPR table

        call nm%deallocate_reverse_hard()

        allocate(nm%fwd_probabilities(0:nm%M-1, 0:nm%M-1))
        if (present(skipLapprTable)) then
            if (skipLapprTable) then
                return
            end if
        end if
        allocate(nm%reverse_hard_lappr_table(0:nm%M-1, 0:nm%bps-1))
    end subroutine noisemapper_allocate_reverse_hard


    module subroutine noisemapper_update_hard_reverse_tables(nm, skipLapprTable)
        !! Update hard reverse reconciliation tables
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        logical, intent(in), optional :: skipLapprTable
        !! Do not compute LAPPR table

        integer :: i, j, k
        real(c_double) :: denominator(0:nm%M-1, 0:nm%bps-1)

        if (present(skipLapprTable)) then
            call nm%allocate_reverse_hard(skipLapprTable)
        end if

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

        if (present(skipLapprTable)) then
            if (skipLapprTable) then
                return
            end if
        end if

        denominator(:,:) = 0
        nm%reverse_hard_lappr_table(:,:) = 0
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

        do i = 0, nm%M - 1
            do k = 0, nm%bps-1
                if (nm%reverse_hard_lappr_table(i, k) == 0) then
                    nm%reverse_hard_lappr_table(i, k) = -1d100
                elseif (denominator(i, k) == 0) then
                    nm%reverse_hard_lappr_table(i, k) = 1d100
                else
                    nm%reverse_hard_lappr_table(i, k) = &
                        log(nm%reverse_hard_lappr_table(i, k)) - log(denominator(i, k))
                end if
            end do
        end do
    end subroutine noisemapper_update_hard_reverse_tables


    module subroutine noisemapper_convert_symbol_to_hard_lappr(nm, x_i, lappr)
        !! Get LAPPR for each symbol from the tables
        class(noisemapper_type), intent(in) :: nm
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
        class(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y
        !! Channel output sample
        real(c_double) :: Fy
        !! CDF of the channel output at `y`

        Fy = sum(nm%probabilities * cdf_normal(x=y, loc=nm%constellation, scale=nm%sigma))
    end function noisemapper_Fy


    module subroutine noisemapper_deallocate_reverse_soft(nm)
        !! Deallocate data used for soft reverse reconciliation
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        if (allocated(nm%monotonicity_configuration)) deallocate(nm%monotonicity_configuration)
        if (allocated(nm%Fy_grid)) deallocate(nm%Fy_grid)
    end subroutine noisemapper_deallocate_reverse_soft


    module subroutine noisemapper_set_monotonicity_array(nm, config)
        !! Set monotonicity configuration
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        logical(c_bool), intent(in) :: config(0:nm%M-1)
        !! Configuration

        if (allocated(nm%monotonicity_configuration)) then
            if (size(nm%monotonicity_configuration) /= nm%M) then
                deallocate(nm%monotonicity_configuration)
            end if
        end if
        if (.not. allocated(nm%monotonicity_configuration)) then
            allocate(nm%monotonicity_configuration(0:nm%M-1))
        end if
        nm%monotonicity_configuration = config
    end subroutine noisemapper_set_monotonicity_array


    module subroutine noisemapper_set_monotonicity_default(nm)
        !! Set default alternating monotonicity configuration
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        if (allocated(nm%monotonicity_configuration)) then
            if (size(nm%monotonicity_configuration) /= nm%M) then
                deallocate(nm%monotonicity_configuration)
            end if
        end if
        if (.not. allocated(nm%monotonicity_configuration)) then
            allocate(nm%monotonicity_configuration(0:nm%M-1))
        end if
        nm%monotonicity_configuration(0::2) = .false.
        nm%monotonicity_configuration(1::2) = .true.
    end subroutine noisemapper_set_monotonicity_default


    module subroutine noisemapper_set_monotonicity_from_integer(nm, config)
        !! Set monotonicity configuration from index
        class(noisemapper_type), intent(inout) :: nm
        !! Noisemapper
        integer(c_int), intent(in) :: config
        !! Configuration number

        logical(c_bool) :: config_l(0:nm%M-1)
        integer :: k, cfg

        cfg = config
        do k = 0, nm%M-1
            config_l(k) = mod(cfg, 2) == 1
            cfg = ishft(cfg, -1)
        end do
        call noisemapper_set_monotonicity_array(nm, config_l)
    end subroutine noisemapper_set_monotonicity_from_integer



    impure elemental module subroutine noisemapper_generate_soft_metric_single(nm, y, n, xhat)
        !! Generate soft metric from a single channel output sample
        !! and give the decided symbol, too.
        class(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y
        !! Channel output sample
        real(c_double), intent(out) :: n
        !! Soft metric
        integer(c_int), intent(out) :: xhat
        !! decided symbol

        xhat = nm%decide_symbol(y)

        n = (nm%Fy(y)- nm%Fy_thresholds(xhat))/nm%delta_Fy(xhat)
        if (nm%monotonicity_configuration(xhat)) then
            n = 1d0 - n
        end if
    end subroutine noisemapper_generate_soft_metric_single


    module subroutine noisemapper_generate_soft_metric_array(nm, y, n, xhat)
        !! Generate soft metric from a set of channel output samples
        class(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: y(:)
        !! Channel output sample
        real(c_double), intent(out) :: n(size(y))
        !! Soft metric
        integer(c_int), intent(out) :: xhat(size(y))
        !! decided symbol

        integer :: i

        do i = 1, size(y)
            call nm%generate_soft_metric(y(i), n(i), xhat(i))
        end do
    end subroutine noisemapper_generate_soft_metric_array


    module subroutine noisemapper_set_Fy_grids(nm, th)
        !! Setup the grid for the inverse of the CDF
        class(noisemapper_type), intent(inout) :: nm
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
        nm%y_grid_step = 5d-4
        n_points = ceiling(2*y_stop/nm%y_grid_step)

        if (allocated(nm%Fy_grid)) then
            deallocate(nm%Fy_grid)
        end if
        allocate(nm%Fy_grid(n_points))

        nm%Fy_grid = nm%Fy(&
            [(nm%base_y_grid + i*nm%y_grid_step, i=1, n_points)])
    end subroutine noisemapper_set_Fy_grids


    module function noisemapper_invert_soft_metric(nm, n, x_i) result(y)
        !! generate all tentative channel output samples from the received soft metric
        class(noisemapper_type), intent(in) :: nm
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


    module function noisemapper_invert_soft_metric_search(nm, n, x_i, res) result(y)
        !! generate all tentative channel output samples from the received soft metric
        class(noisemapper_type), intent(in) :: nm
        !! noise mapper
        real(c_double), intent(in) :: n
        !! Soft metric
        integer(c_int), intent(in) :: x_i
        !! Hypotetical received symbol alphabet index
        real(c_double), intent(in), optional :: res
        !! Resolution
        real(c_double) :: y
        !! Tentative channel output samples

        real(c_double) :: Fy, resolution
        integer :: idx

        if (present(res)) then
            resolution = res
        else
            resolution = 1d-12
        end if

        Fy = n
        if (nm%monotonicity_configuration(x_i)) then
            Fy = 1-Fy
        end if
        Fy = nm%Fy_thresholds(x_i) + Fy*nm%delta_Fy(x_i)

        idx = binsearch(nm%Fy_grid, Fy)


        if ((idx == 0) .or. (idx == size(nm%Fy_grid))) then
            y = nm%invert_Fy(Fy, res=resolution)
        else
            y = nm%invert_Fy(Fy, res=resolution, &
                ybounds=(nm%base_y_grid + (idx + [0, 1])*nm%y_grid_step))
        end if
    end function noisemapper_invert_soft_metric_search


    module subroutine noisemapper_soft_reverse_lappr_single(nm, x_i, n, lappr, res)
        !! Calculate the LAPPR from the transmitted symbol and the soft metric
        class(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(in) :: x_i
        !! transmitted symbol alphabet index
        real(c_double), intent(in) :: n
        !! Soft metric
        real(c_double), intent(out) :: lappr(0:nm%bps-1)
        !! LAPPRs
        real(c_double), intent(in), optional :: res
        !! If present, it toggles the soft search down to `res` as resolution

        integer :: i !! alphabet index of tentative received symbol
        integer :: k !! further alphabet/bit index
        real(c_double) :: denominator(0:nm%bps-1)
        real(c_double) :: twoy, x, addendum


        lappr(:) = 0
        denominator(:) = 0

        x = nm%constellation(x_i)

        do i = 0, nm%M-1
            if (present(res)) then
                twoy = 2*nm%invert_soft_metric_search(n, i, res)
            else
                twoy = 2*nm%invert_soft_metric(n, i)
            end if

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

        do k = 0, nm%bps - 1
            if (denominator(k) == 0) then
                lappr(k) = 1d100
            elseif (lappr(k) == 0) then
                lappr(k) = -1d100
            else
                lappr(k) = log(lappr(k)) - log(denominator(k))
            end if
        end do
    end subroutine noisemapper_soft_reverse_lappr_single


    module subroutine noisemapper_soft_reverse_lappr_array(nm, x_i, n, lappr, res)
        !! Calculate the LAPPR from the transmitted symbols and the soft metrics arrays
        class(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(in) :: x_i(0:)
        !! transmitted symbol alphabet index
        real(c_double), intent(in) :: n(0:size(x_i)-1)
        !! Soft metric
        real(c_double), intent(out) :: lappr(0:size(x_i)*nm%bps-1)
        !! LAPPRs
        real(c_double), intent(in), optional :: res
        !! If present, it toggles the soft search down to `res` as resolution

        integer :: i

        if (present(res)) then
            do i = 0, size(x_i)-1
                call noisemapper_soft_reverse_lappr(nm, x_i(i), n(i), lappr(i*nm%bps : (i+1)*nm%bps-1), res)
            end do
        else
            do i = 0, size(x_i)-1
                call noisemapper_soft_reverse_lappr(nm, x_i(i), n(i), lappr(i*nm%bps : (i+1)*nm%bps-1))
            end do
        end if
    end subroutine noisemapper_soft_reverse_lappr_array


    ! +-----------------------------+
    ! | Uniform probabilities at RX |
    ! +-----------------------------+
    module function noisemapper_inverse_Fy_search(nm, Fy, res, ybounds) result (y)
        !! Find the `y` value whose CDF is `Fy`, within a certain `res`-olution
        !! on the CDF
        class(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: Fy
        !! CDF value to be inverted
        real(c_double), intent(in), optional :: res
        !! Resolution of the search result
        real(c_double), intent(in), optional :: ybounds(2)
        !! Initial lower and upper bound
        real(c_double) :: y
        !! output channel value whose CDF is `Fy`

        real(c_double) :: y_l, y_h, resolution, Fy_l, Fy_h, Fy_next

        if (Fy == 0) then
           y = -1d300
           return
        else if (Fy==1) then
           y = 1d300
           return
        end if

        ! Setup resolution
        if (present(res)) then
            resolution = res
        else
            resolution = 1d-12
        end if

        if (.not. present(yBounds)) then
            go to 50
        else
            y_l = minval(ybounds)
            y_h = maxval(ybounds)

            Fy_l = nm%Fy(y_l)
            Fy_h = nm%Fy(y_h)
            if ((Fy .lt. Fy_l) .or. (Fy .gt. Fy_h)) then
                go to 50
            end if
            go to 100 ! find_y
        end if

50      Fy_l = nm%Fy(0d0)
        if (Fy_l .gt. Fy) then
            y_h = 0
            y_l = -1
            Fy_h = Fy_l
            Fy_l = nm%Fy(y_l)
            do while (Fy_l .gt. Fy)
                y_h = y_l
                y_l = 2*y_l
                Fy_h = Fy_l
                Fy_l = nm%Fy(y_l)
            end do
        else
            y_l = 0
            y_h = 1
            Fy_h = nm%Fy(y_h)
            do while (Fy_h .lt. Fy)
                y_l = y_h
                y_h = 2*y_h
                Fy_l = Fy_h
                Fy_h = nm%Fy(y_h)
            end do
        end if

100     find_y: do while (.true.)
            y = y_l + (Fy-Fy_l)*(y_h-y_l)/(Fy_h - Fy_l) ! linear interpolation
            Fy_next = nm%Fy(y)
            if (abs(Fy-Fy_next) .lt. resolution) then
                return
            end if
            if (Fy_next .gt. Fy) then
                y_h = y
                Fy_h = Fy_next
            else
                y_l = y
                Fy_l = Fy_next
            end if
        end do find_y
    end function noisemapper_inverse_Fy_search


    module subroutine noisemapper_set_y_thresholds_uniform(nm)
        !! Set thresholds with uniform decision probabilities
        class(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        integer :: i
        real(c_double) :: thresholds(1:nm%M-1)

        do i = 1, nm%M-1
            thresholds(i) = nm%invert_Fy(real(i, c_double)/real(nm%M, c_double))
        end do
        call nm%set_y_thresholds(thresholds)
    end subroutine noisemapper_set_y_thresholds_uniform


    real(c_double) impure elemental module function f_xhat_n_cond_x(nm, n, xhat, x) result(pdf)
        !! PDF of \(N, \hat{X}|X\)
        class(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: n
        !! Soft metric
        integer(c_int), intent(in) :: xhat
        !! Bob's decided symbol (index within constellation)
        integer(c_int), intent(in) :: x
        !! Alice's transmitted symbol (index within constellation)

        integer :: k
        real(c_double) :: a_j, two_y_i

        pdf = 0

        a_j = nm%constellation(x)
        two_y_i = 2*nm%invert_soft_metric_search(n, xhat)
        do k = 0, nm%M-1
            pdf = pdf + nm%probabilities(k) * &
                exp((nm%constellation(k) - a_j)*(two_y_i - nm%constellation(k) - a_j)/nm%N0)
        end do
        pdf = nm%delta_Fy(xhat) / pdf
    end function f_xhat_n_cond_x
end submodule re2often_noisemapper
