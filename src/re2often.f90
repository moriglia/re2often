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
module re2often
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    use, intrinsic :: iso_c_binding
    ! use re2often_utils, only: binsearch
    ! use stdlib_stats_distribution_normal, only: cdf_normal
    implicit none

    type, public :: noisemapper_type
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

        real(c_double), allocatable :: alice_bit_priors(:)
        !! Log a priori probabilities of each bit

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
    contains
        final :: noisemapper_finalize
        procedure, pass :: deallocate => noisemapper_deallocate
        !! Deallocate basic arrays
        procedure, pass :: set_symbol_probabilities => noisemapper_set_symbol_probabilities
        !! Set symbol probabilities
        procedure, pass :: noisemapper_random_symbol_single
        procedure, pass :: noisemapper_random_symbol_array
        generic         :: random_symbol => noisemapper_random_symbol_single, noisemapper_random_symbol_array
        !! Generate random symbols
        procedure, pass :: symbol_index_to_value => noisemapper_symbol_index_to_value
        !! Convert symbol index to constellation value
        procedure, pass :: noisemapper_y_to_lappr_single
        procedure, pass :: noisemapper_y_to_lappr_array
        generic         :: y_to_lappr => noisemapper_y_to_lappr_array, noisemapper_y_to_lappr_single
        !! Compute the LAPPR for direct reconciliation (Bob is to guess Alice's tx sequence)
        procedure, pass :: noisemapper_y_to_llr_single
        procedure, pass :: noisemapper_y_to_llr_array
        generic         :: y_to_llr => noisemapper_y_to_llr_array, noisemapper_y_to_llr_single
        !! Compute the LLR for direct reconciliation
        procedure, pass :: update_alice_priors => noisemapper_update_alice_priors
        !! Update prior log-ratios on Alice's bits for direct reconciliation
        procedure, pass :: decide_symbol => noisemapper_decide_symbol_single
        !! Take a decision on a channel output (Bob's side)
        procedure, pass :: update_N0_from_snrdb => noisemapper_update_N0_from_snrdb
        !! Use SNR (dB) to update the value of the noise spectral density (N_0/2 on each quadrature)
        procedure, pass :: Fy => noisemapper_Fy
        !! Compute the CDF of the output symbol
        procedure, pass :: invert_Fy => noisemapper_inverse_Fy_search
        !! Compute the inverse CDF of the output symbol
        procedure, pass :: set_Fy_grids => noisemapper_set_Fy_grids
        !! Setup a pre-computed array of CDF values for the channel output

        ! +----------------------------+
        ! | Encoding related functions |
        ! +----------------------------+
        procedure, pass :: set_encoding => noisemapper_set_encoding_custom
        procedure, pass :: set_encoding_gray => noisemapper_set_encoding_gray
        procedure, pass :: set_encoding_natural => noisemapper_set_encoding_natural
        procedure, pass :: symbol_to_word => noisemapper_symbol_to_word

        ! +-------------------------------------------+
        ! | Generic reverse reconciliation procedures |
        ! +-------------------------------------------+
        procedure, pass :: deallocate_reverse_common => noisemapper_deallocate_reverse_common
        !! Deallocate arrays used for reverse reconciliation
        procedure, pass :: allocate_reverse_common => noisemapper_allocate_reverse_common
        !! Allocate common arrays for reverse reconciliation
        procedure, pass :: set_y_thresholds => noisemapper_set_y_thresholds
        !! Set the thresholds for hard decisions on channel output (Bob's side)
        procedure, pass :: set_y_thresholds_uniform => noisemapper_set_y_thresholds_uniform
        !! Set adaptive thresholds for uniformly distributed
        !! hard decisions on channel output (Bob's side)

        ! +---------------------------------------------------+
        ! | Reverse reconciliation with hard information only |
        ! +---------------------------------------------------+
        procedure, pass :: deallocate_reverse_hard => noisemapper_deallocate_reverse_hard
        !! Deallocate arrays specifically for hard reverse reconciliation
        procedure, pass :: allocate_reverse_hard => noisemapper_allocate_reverse_hard
        !! Allocate arrays specifically for hard reverse reconciliation
        procedure, pass :: update_hard_reverse_tables => noisemapper_update_hard_reverse_tables
        !! Update lookup tables for hard reverse reconciliation
        !! (fwd_probabilities and reverse_hard_lappr_table)
        procedure, pass :: hard_lappr => noisemapper_convert_symbol_to_hard_lappr
        !! Compute hard lappr from transitted symbol (Alice's side)

        ! +----------------------------------------------+
        ! | Reverse reconciliation with soft information |
        ! +----------------------------------------------+
        procedure, pass :: deallocate_reverse_soft => noisemapper_deallocate_reverse_soft
        !! Deallocate monotonicity configuration and Fy grid
        procedure, pass :: generate_soft_metric => noisemapper_generate_soft_metric_single, &
            noisemapper_generate_soft_metric_array
        !! Use the channel output (Bob's side) to generate the soft metric.
        !! At the same time, this function also provides a discretized version
        !! of the channel output
        procedure, pass :: noisemapper_soft_reverse_lappr_single
        procedure, pass :: noisemapper_soft_reverse_lappr_array
        generic         :: soft_reverse_lappr => noisemapper_soft_reverse_lappr_single, &
            noisemapper_soft_reverse_lappr_array
        !! Compute the LAPPRs (Alice's side) from the channel inputs,
        !! and the soft metric provided by Bob.
        procedure, pass :: set_monotonicity => noisemapper_set_monotonicity_default, &
            noisemapper_set_monotonicity_from_integer, &
            noisemapper_set_monotonicity_array
        !! Set the monotonicity configuration for the transformation functions
        !! to be used by Bob.
        procedure, pass :: invert_soft_metric => noisemapper_invert_soft_metric
        !! Coarse inversion of the soft metric
        procedure, pass :: invert_soft_metric_search => noisemapper_invert_soft_metric_search
        !! Search-based  inversion of the soft metric


        ! +---------------------------------+
        ! | Theoretical entropic quantities |
        ! +---------------------------------+
        procedure, pass :: f_n_xhat_cond_x
        !! Joint distribution of the soft metric and the discretized channel output
        !! given the input symbol
        procedure, pass :: H_Xhat
        !! Entropy of the discretized channel output.
        procedure, pass :: H_Xhat_cond_X
        !! Entropy of the discretized channel output given the input symbol
        procedure, pass :: I_direct
        !! Mutual information of the direct channel
        procedure, pass :: I_hard_reverse
        !! Mutual information of the hard reverse reconciliation scheme
        procedure, pass :: I_hard_reverse_equidistant_th
        !! Mutual information of the hard reverse scheme
        !! with maximum likelyhood thresholds
        procedure, pass :: I_hard_reverse_uniform_output_th
        !! Mutual information of the hard reverse scheme
        !! with thresholds set to have uniform output probability
        procedure, pass :: I_soft_reverse
        !! Mutual information of the soft reverse reconciliation scheme
        procedure, pass :: I_soft_reverse_equidistant_th
        !! Mutual information of the soft reverse reconciliation scheme
        !! with maximum likelyhood thresholds
        procedure, pass :: I_soft_reverse_uniform_output_th
        !! Mutual information of the soft reverse reconciliation scheme
        !! with thresholds yielding uniform output symbol probability

        procedure, pass :: I_s_ml_soft_direct
        procedure, pass :: I_s_ml_hard_direct
        procedure, pass :: I_s_ml_soft_reverse
        procedure, pass :: I_s_map_soft_reverse
        procedure, pass :: I_s_map_soft_direct
        procedure, pass :: I_s_map_hard_reverse
    end type noisemapper_type


    ! +------------------------------------+
    ! | Creation of the noisemapper object |
    ! +------------------------------------+
    interface
        module subroutine noisemapper_finalize(nm)
            !! Wrapper for the deallocate function
            type(noisemapper_type) :: nm
            !! Noisemapper
        end subroutine noisemapper_finalize
        module subroutine noisemapper_deallocate(nm)
            !! Destructor for noise mapper
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
        end subroutine noisemapper_deallocate
        module subroutine noisemapper_set_symbol_probabilities(nm, probabilities)
            !! Allocate and set probability vector for imput constellation symbols
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
            real(c_double), intent(in), optional :: probabilities(0:nm%M-1)
            !! Input probabilities
        end subroutine noisemapper_set_symbol_probabilities
    end interface


    ! +--------------------------------------+
    ! | Interfaces for DIRECT reconciliation |
    ! +--------------------------------------+
    interface noisemapper_y_to_lappr
        module subroutine noisemapper_y_to_lappr_single(nm, y, lappr)
            !! calculate lappr from channel output sample for direct reconciliation
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
            real(c_double), intent(in) :: y
            !! AWGN channel output sample
            real(c_double), intent(out) :: lappr(0:nm%bps - 1)
            !! log-a posteriori-probabilities for the de-mapped bits
        end subroutine noisemapper_y_to_lappr_single
        module subroutine noisemapper_y_to_lappr_array(nm, y, lappr)
            !! calculate lappr from set of channel output samples for direct reconciliation
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
            real(c_double), intent(in) :: y(0:)
            !! AWGN channel samples
            real(c_double), intent(out) :: lappr(0:size(y)*nm%bps-1)
            !! LAPPR corresponding to the channel outputs
        end subroutine noisemapper_y_to_lappr_array
    end interface noisemapper_y_to_lappr

    interface
        module subroutine noisemapper_random_symbol_single(nm, x_i)
            !! Generate a random symbol
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
            integer(c_int), intent(out) :: x_i
            !! Random symbol of the constellation (index in 0:M-1)
        end subroutine noisemapper_random_symbol_single
        module subroutine noisemapper_random_symbol_array(nm, x_i)
            !! Generate random symbols
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
            integer(c_int), intent(out) :: x_i(:)
            !! Random symbols of the constellation (index in 0:M-1)
        end subroutine noisemapper_random_symbol_array
    end interface

    interface
        module subroutine noisemapper_y_to_llr_single(nm, y, llr)
            !! Calculate LLR from the channel output
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
            real(c_double), intent(in) :: y
            !! AWGN channel output sample
            real(c_double), intent(out) :: llr(0:nm%bps-1)
            !! log-likelyhood ratios of the bits associated to one single symbol transmission
        end subroutine noisemapper_y_to_llr_single
        module subroutine noisemapper_y_to_llr_array(nm, y, llr)
            !! Calculate LLR from the channel output for a set of channel output samples
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
            real(c_double), intent(in) :: y(0:)
            !! AWGN channel output sample
            real(c_double), intent(out) :: llr(0:nm%bps*sizeof(y)-1)
            !! log-likelyhood ratios of the bits associated to one single symbol transmission
        end subroutine noisemapper_y_to_llr_array
        module subroutine noisemapper_update_alice_priors(nm)
            !! Update a priori probabilities for Alice
            !! This function returns with no error nor warning
            !! if either the encoding or the probabilities aren't set
            class(noisemapper_type), intent(inout) :: nm
            !! Noisemapper
        end subroutine noisemapper_update_alice_priors
    end interface

    ! +---------------------------------------+
    ! | Interfaces for REVERSE reconciliation |
    ! +---------------------------------------+
    interface
        module subroutine noisemapper_allocate_reverse_common(nm)
            !! Allocate the common reverse reconciliation buffers
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
        end subroutine noisemapper_allocate_reverse_common
        module subroutine noisemapper_deallocate_reverse_common(nm)
            !! Deallocation of common arrays for reverse reconciliation
            class(noisemapper_type), intent(inout) :: nm
            !! noise mapper
        end subroutine noisemapper_deallocate_reverse_common
    end interface
    interface noisemapper_decide_symbol
        impure elemental module function noisemapper_decide_symbol_single(nm, y) result(x_i)
            !! Take a decision for the received channel output based
            !! on the thresholds
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
            real(c_double), intent(in) :: y
            !! Channel output sample
            integer(c_int) :: x_i
            !! Alphabet index of the decided symbol
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
        ! end function noisemapper_decide_symbol_array
    end interface noisemapper_decide_symbol

    ! +--------------------------------------------+
    ! | Interfaces for SOFT REVERSE reconciliation |
    ! +--------------------------------------------+
    interface noisemapper_generate_soft_metric
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
        end subroutine noisemapper_generate_soft_metric_array
    end interface noisemapper_generate_soft_metric

    interface noisemapper_soft_reverse_lappr
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
        end subroutine noisemapper_soft_reverse_lappr_array
    end interface noisemapper_soft_reverse_lappr

    interface
        module subroutine noisemapper_set_monotonicity_default(nm)
            !! Set default alternating monotonicity configuration
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
        end subroutine noisemapper_set_monotonicity_default
        module subroutine noisemapper_set_monotonicity_array(nm, config)
            !! Set monotonicity configuration
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
            logical(c_bool), intent(in) :: config(0:nm%M-1)
            !! Configuration
        end subroutine noisemapper_set_monotonicity_array
        module subroutine noisemapper_set_monotonicity_from_integer(nm, config)
            !! Set monotonicity configuration from index
            class(noisemapper_type), intent(inout) :: nm
            !! Noisemapper
            integer(c_int), intent(in) :: config
            !! Configuration number
        end subroutine noisemapper_set_monotonicity_from_integer
    end interface


    interface
        module function noisemapper_create(bps) result(nm)
            !! Create the nm object
            integer(c_int), intent(in) :: bps
            !! bit per symbol
            type(noisemapper_type) :: nm
            !! Noise mapper
            !! This function will not setup N0
        end function noisemapper_create
        module subroutine noisemapper_update_N0_from_snrdb(nm, snrdb)
            !! Update N0 based on the value of the SNR
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
            real(c_double), intent(in) :: snrdb
            !! SNR in dB
        end subroutine noisemapper_update_N0_from_snrdb
        module subroutine noisemapper_set_y_thresholds_uniform(nm)
            !! Set thresholds with uniform decision probabilities
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
        end subroutine noisemapper_set_y_thresholds_uniform
        module subroutine noisemapper_set_y_thresholds(nm, thresholds)
            !! Set the decision thresholds
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
            real(c_double), intent(in), optional :: thresholds(1:nm%M-1)
            !! y thresholds. If not present, the thresholds are the points
            !! half way between two adjacent constellation points
            !! @warning The order of the array is not checked
        end subroutine noisemapper_set_y_thresholds
        module subroutine noisemapper_update_hard_reverse_tables(nm, skipLapprTable)
            !! Update hard reverse reconciliation tables
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
            logical, intent(in), optional :: skipLapprTable
            !! Do not compute LAPPR table
        end subroutine noisemapper_update_hard_reverse_tables
        module subroutine noisemapper_set_Fy_grids(nm, th)
            !! Setup the grid for the inverse of the CDF
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
            real(c_double), intent(in), optional :: th
            !! threshold for the PDF minimum value.
            !! It must be strictly positive, lower than 1
        end subroutine noisemapper_set_Fy_grids
        elemental module function noisemapper_Fy(nm, y) result(Fy)
            !! Evaluate the CDF of the output at the given point
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
            real(c_double), intent(in) :: y
            !! Channel output sample
            real(c_double) :: Fy
            !! CDF of the channel output at `y`
        end function noisemapper_Fy
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
        end function noisemapper_inverse_Fy_search
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
        end function noisemapper_invert_soft_metric_search
        module subroutine noisemapper_set_encoding_gray(nm)
            !! Set Gray encoding
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
        end subroutine noisemapper_set_encoding_gray
        module subroutine noisemapper_set_encoding_natural(nm)
            !! Set Gray encoding
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
        end subroutine noisemapper_set_encoding_natural
        module subroutine noisemapper_set_encoding_custom(nm, labels)
            !! Set custom encoding labels
            !! @warning: no check on labels being a complete permutation of (0:M-1)
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
            integer, intent(in) :: labels(0:nm%M-1)
            !! Encoding Labels
        end subroutine noisemapper_set_encoding_custom
        module function noisemapper_symbol_index_to_value(nm, x_i) result (x)
            !! Convert constellation index to point
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
            integer(c_int), intent(in) :: x_i(:)
            !! Set of constellation indexes
            real(c_double) :: x(size(x_i))
            !! set of constellation points
        end function noisemapper_symbol_index_to_value
        module function noisemapper_symbol_to_word(nm, x_i) result (word)
            !! Convert a set of symbol indexes to a word
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
            integer(c_int), intent(in) :: x_i(0:)
            !! Set of indexes
            logical(c_bool) :: word(0:size(x_i)*nm%bps - 1)
            !! Word corresponding to the sequence of symbols
        end function noisemapper_symbol_to_word
        module subroutine noisemapper_allocate_reverse_hard(nm, skipLapprTable)
            !! Allocate transition probability table and lappr table
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
            logical, intent(in), optional :: skipLapprTable
            !! Do not allocate LAPPR table
        end subroutine noisemapper_allocate_reverse_hard
        module subroutine noisemapper_deallocate_reverse_hard(nm)
            !! Deallocate transition probability table and lappr table
            class(noisemapper_type), intent(inout) :: nm
            !! Noise mapper
        end subroutine noisemapper_deallocate_reverse_hard
        module subroutine noisemapper_convert_symbol_to_hard_lappr(nm, x_i, lappr)
            !! Get LAPPR for each symbol from the tables
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
            integer(c_int), intent(in) :: x_i(0:)
            !! Transmitted symbols
            real(c_double), intent(out) :: lappr(0:nm%bps*size(x_i)-1)
            !! LAPPR array associated with the transmitted sybmols
        end subroutine noisemapper_convert_symbol_to_hard_lappr
        module subroutine noisemapper_deallocate_reverse_soft(nm)
            !! Deallocate monotonicity configuration and Fy grid
            class(noisemapper_type), intent(inout) :: nm
        end subroutine noisemapper_deallocate_reverse_soft
    end interface


    interface
        ! +------------------------+
        ! | MI submodule interface |
        ! +------------------------+
        real(c_double) impure module elemental function f_n_xhat_cond_x(nm, n, xhat, x) result(pdf)
            !! Initialized noisemapper object
            class(noisemapper_type), intent(in) :: nm
            !! PDF of \(N, \hat{X}|X\)
            real(c_double), intent(in) :: n
            !! Soft metric
            integer(c_int), intent(in) :: xhat
            !! Bob's decided symbol (index within constellation)
            integer(c_int), intent(in) :: x
            !! Alice's transmitted symbol (index within constellation)
        end function f_n_xhat_cond_x
        real(c_double) module function H_Xhat(nm) result(H)
            !! Entropy of the output symbols
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
        end function H_Xhat
        real(c_double) module function H_Xhat_cond_X(nm) result(H)
            !! Entropy of the output symbols
            class(noisemapper_type), intent(in) :: nm
            !! Noise mapper
        end function H_Xhat_cond_X
        real(c_double) module function I_hard_reverse(nm) result (I)
            !! Mutual information of the discrete Input and Output channel
            class(noisemapper_type), intent(inout) :: nm
            !! Initialized noisemapper object
        end function I_hard_reverse
        real(c_double) module function I_hard_reverse_uniform_output_th(nm, snrdb) result (I)
            !! Mutual information of the discrete Input and Output channel
            class(noisemapper_type), intent(inout) :: nm
            !! Initialized noisemapper object
            real(c_double), intent(in) :: snrdb
            !! SNR [dB] at which to evaluate the Mutual information
        end function I_hard_reverse_uniform_output_th
        real(c_double) module function I_hard_reverse_equidistant_th(nm, snrdb) result (I)
            !! Mutual information of the discrete Input and Output channel
            class(noisemapper_type), intent(inout) :: nm
            !! Initialized noisemapper object
            real(c_double), intent(in) :: snrdb
            !! SNR [dB] at which to evaluate the Mutual information
        end function I_hard_reverse_equidistant_th
        real(c_double) module function I_direct(nm, npts) result(I)
            !! Initialized noisemapper object
            class(noisemapper_type), intent(inout) :: nm
            !! Mutual information of the direct reconciliation scheme
            integer(c_int), intent(in), optional :: npts
        end function I_direct
        real(c_double) module function I_soft_reverse(nm) result(I)
            !! Mutual information of the soft reverse reconciliation scheme
            class(noisemapper_type), intent(in) :: nm
            !! Initialized noisemapper object
        end function I_soft_reverse
        real(c_double) module function I_soft_reverse_equidistant_th(nm, snrdb) result(I)
            !! Mutual information of the soft reverse reconciliation scheme
            class(noisemapper_type), intent(inout) :: nm
            !! Initialized noisemapper object
            real(c_double), intent(in) :: snrdb
            !! SNR [dB] at which to calculate the mutual information
        end function I_soft_reverse_equidistant_th
        real(c_double) module function I_soft_reverse_uniform_output_th(nm, snrdb) result(I)
            !! Mutual information of the soft reverse reconciliation scheme
            class(noisemapper_type), intent(inout) :: nm
            !! Initialized noisemapper object
            real(c_double), intent(in) :: snrdb
            !! SNR [dB] at which to calculate the mutual information
        end function I_soft_reverse_uniform_output_th
    end interface

    integer :: dqags_Limit
    interface
        ! +-------------------------+
        ! | GMI submodule interface |
        ! +-------------------------+
        module function I_s_ml_soft_reverse(nm, s) result(I_s)
            !! Compute the GMI in the Maximum-Likelyhood version
            class(noisemapper_type),   intent(in) :: nm
            real(c_double), optional, intent(in) :: s
            real(c_double)                       :: I_s
        end function I_s_ml_soft_reverse
        module function I_s_ml_hard_direct(nm, s) result (I_s)
            class(noisemapper_type), intent(in) :: nm
            real(c_double), intent(in), optional :: s
            real(c_double) :: I_s
        end function I_s_ml_hard_direct
        module function I_s_map_soft_direct(nm, s) result(I_s)
            !! GMI for MAP criterion with hard information only, reverse direction
            class(noisemapper_type), intent(in) :: nm
            !! Initialized noisemapper object
            real(c_double), intent(in), optional :: s
            !! Positive parameter s of the GMI
            real(c_double) :: I_s
            !! GMI(s)
        end function I_s_map_soft_direct
        module function I_s_map_hard_reverse(nm, s) result(I_s)
            !! GMI for MAP criterion with hard information only, reverse direction
            class(noisemapper_type), intent(in) :: nm
            !! Initialized noisemapper object
            real(c_double), intent(in), optional :: s
            !! Positive parameter s of the GMI
            real(c_double) :: I_s
            !! GMI(s)
        end function I_s_map_hard_reverse
        module function I_s_map_soft_reverse(nm, s, useDenominator) result(I_s)
            !! GMI-MAP
            class(noisemapper_type), intent(in)  :: nm
            real(c_double), intent(in), optional :: s
            logical       , intent(in), optional :: useDenominator
            real(c_double)                       :: I_s
        end function I_s_map_soft_reverse
        module function I_s_ml_soft_direct(nm, s) result(I_s)
            !! Compute the GMI-ML for the direct channel in case of generic probability
            class(noisemapper_type),  intent(in) :: nm
            real(c_double), optional, intent(in) :: s
            real(c_double)                       :: I_s
        end function I_s_ml_soft_direct
    end interface
end module re2often
