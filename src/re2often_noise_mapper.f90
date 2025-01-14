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
module re2often_noise_mapper
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! This module contains all Pulse Amplitude Modulation functionalities
    !! for the Soft Reverse Reconciliation
    use iso_fortran_env, only: dp => real64
    use stdlib_stats_distribution_normal, only: cdf_normal
    use re2often_utils, only: binsearch
    implicit none

    private

    public :: TNoiseMapper

    type :: TNoiseMapper
        !! Object with both random symbol generation and noise (de)mapping
        !! for the Reconciliation Softening
        integer, public :: bps
        !! Bit per symbol
        integer, public :: M
        !! Order of modulation (number of symbols in the alphabet)
        double precision, public, allocatable :: constellation(:)
        !! Points of the constellation (0:M-1)
        double precision, public, allocatable :: probabilities(:)
        !! Probability of each constellation point (0:M-1)
        double precision, public :: E_symbol
        !! Symbol Variance
        double precision, public :: N0
        !! Total complex noise variance (Actual on each quadrature is N0/2)
        double precision, public :: sigma
        !! sqrt(N0/2)
        logical, allocatable :: symbol_to_bit_map(:,:)
        !! Bits (rows) associated to each symbol (row index), size is (0:M-1, 0:bps-1)

        double precision, public, allocatable :: y_thresholds(:)
        !! Decision thresholds (1:M-1)
        !! Decision threshold `i` divides the decision
        !! regions of symbols `i-1` and `i`
        double precision, public, allocatable :: Fy_thresholds(:)
        !! Output sample cumulative density function at thresholds (0:M)
        !! with `Fy_thresholds(0)=`\(F_Y(-\infty)\)
        !! and  `Fy_thresholds(M)=`\(F_Y(+\infty)\)
        double precision, public, allocatable :: deltaFy(:)
        !! Probability that the channel output lays in each symbol decision region
        !! calculated as \(P\{y\in D_i\} = F_Y(\sup D_i) - F_Y(\inf D_i)\)
        !! where \(D_i\) denotes the decision interval for symbol \(i\), so
        !! the infimum and sumpremum are the thresholds between symbol \(a_i\) and
        !! the adjacent symbols
        logical, public, allocatable          :: negative_monotonicity(:)
        !! Monotonicity configuration, if true, the softening metric `nhat` decreases
        !! with `y` increasing in the region corresponding to the index of the configuration vector,
        !! otherwise `nhat` increases with `y` increasing
        !! @warning Not implemented: setting a custom configuration

        double precision, private, allocatable :: y_grid(:)
        !! grid of channel output points `(0:npoints-1)`, with `npoints` not saved as part of this type
        double precision, private, allocatable :: F_grid(:)
        !! grid of distribution values corresponding to the `y_grid` points.
        !! Its range starts from 0, i.e. its range is `(0:npoints-1)`, `npoints` is not saved as part of the type

        ! integer, public, allocatable :: Fy_threshold_index(:)
        ! !! Index of the threshold within the grid of the CDF values

        double precision, public, allocatable :: fwd_probability(:,:)
        !! transition probability (TX input to hard RX output), ranges (0:M-1 , 0:M-1)
        !! Rows represent transmitted symbols, columns received symbols
        !! \(P\{\hat{X}=a_j|X=a_i\}=\) `fwd_probabilities(i,j)`
        double precision, public, allocatable :: bwd_probability(:,:)
        !! transition probability (RX hard decided symbol to possible TX symbol),
        !! ranges (0:M-1 , 0:M-1)
        !! Rows represent transmitted symbols, columns received symbols
        !! \(P\{X=a_i|\hat{X}=a_j\}=\) `bwd_probabilities(i,j)`
        double precision, public, allocatable :: lappr_hard(:,:)
        !! Log-a posteriori-probability ratio for hard reverse reconciliation
        !! ranges: (0:M-1 , 0:bps-1)
    contains
        procedure, public, pass :: free_noise_mapper
        procedure, public, pass :: deallocate_thresholds
        procedure, public, pass :: deallocate_transition_probabilities
        procedure, public, pass :: random_symbols
        procedure, public, pass :: symbol_index_to_value
        procedure, public, pass :: update_N0
        procedure, public, pass :: y_to_lappr_single
        procedure, public, pass :: y_to_lappr_array
        generic, public         :: y_to_lappr => y_to_lappr_single, y_to_lappr_array
        procedure, public, pass :: symbol_to_word_single
        procedure, public, pass :: symbol_to_word_array
        generic, public         :: symbol_to_word => symbol_to_word_single, symbol_to_word_array
        procedure, public, pass :: set_y_thresholds
        procedure, public, pass :: set_grid
        procedure, public, pass :: decide_symbol
        procedure, public, pass :: cdf_y
        procedure, public, pass :: generate_soft_metric
        procedure, public, pass :: generate_tentative_channel_sample
        procedure, public, pass :: generate_lappr_single
        procedure, public, pass :: generate_lappr_array
        generic, public         :: generate_lappr => generate_lappr_single, generate_lappr_array
        procedure, public, pass :: update_transition_probabilities
        procedure, public, pass :: convert_symbol_to_hard_lappr_single
        procedure, public, pass :: convert_symbol_to_hard_lappr_array
        generic, public         :: convert_symbol_to_hard_lappr => &
            convert_symbol_to_hard_lappr_single, convert_symbol_to_hard_lappr_array
        final :: TNoiseMapperDestructor
    end type TNoiseMapper


    interface TNoiseMapper
        module procedure TNoiseMapperConstructor
    end interface TNoiseMapper
contains

    subroutine deallocate_thresholds(this)
        !! Deallocation of threshold array
        class(TNoiseMapper), intent(inout) :: this
        !! noise mapper

        if (allocated(this%y_thresholds )) deallocate(this%y_thresholds)
        if (allocated(this%Fy_thresholds)) deallocate(this%Fy_thresholds)
        if (allocated(this%deltaFy      )) deallocate(this%deltaFy)
        if (allocated(this%F_grid       )) deallocate(this%F_grid)
        if (allocated(this%y_grid       )) deallocate(this%y_grid)
        ! if (allocated(this%Fy_threshold_index)) deallocate(this%Fy_threshold_index)
    end subroutine deallocate_thresholds

    subroutine deallocate_transition_probabilities(this)
        !! Deallocate forward and backward transition probability arrays
        class(TNoiseMapper), intent(inout) :: this
        !! noise mapper

        if (allocated(this%fwd_probability)) deallocate(this%fwd_probability)
        if (allocated(this%bwd_probability)) deallocate(this%bwd_probability)
        if (allocated(this%lappr_hard     )) deallocate(this%lappr_hard     )
    end subroutine deallocate_transition_probabilities

    subroutine free_noise_mapper(nm)
        !! Deallocates all allocatable variables within TNoiseMapper type
        class(TNoiseMapper), intent(inout) :: nm
        !! The noise mapper whose allocatable variables are to be deallocated

        if (allocated(nm%constellation))     deallocate(nm%constellation)
        if (allocated(nm%probabilities))     deallocate(nm%probabilities)
        if (allocated(nm%symbol_to_bit_map)) deallocate(nm%symbol_to_bit_map)
        call nm%deallocate_thresholds
        if (allocated(nm%negative_monotonicity)) deallocate(nm%negative_monotonicity)
        call nm%deallocate_transition_probabilities
    end subroutine free_noise_mapper


    subroutine TNoiseMapperDestructor(nm)
        !! Finalizer function for TNoiseMapper (Wrapper for the deallocator function)
        type(TNoiseMapper), intent(inout) :: nm
        !! Noise mapper to be finalized
        call nm%free_noise_mapper
    end subroutine TNoiseMapperDestructor


    function TNoiseMapperConstructor(bps, N0, probabilities) result (nm)
        !! Constructor for the TNoiseMapper class
        integer, intent(in) :: bps
        !! bit per symbol
        double precision, intent(in) :: N0
        !! Total noise variance on both quadratures
        double precision, intent(in), optional :: probabilities(0:ishft(1, bps)-1)
        !! Array of probabilities for each symbol

        type(TNoiseMapper) :: nm
        !! Constructed TNoiseMapper object

        integer :: i, k

        call nm%free_noise_mapper

        nm%bps = bps
        nm%M = ishft(1, bps)

        allocate(nm%constellation(0:nm%M-1))
        allocate(nm%probabilities(0:nm%M-1))
        allocate(nm%symbol_to_bit_map(0:nm%M-1 , 0:bps-1))

        nm%constellation(:) = &
            [(real(1-nm%M, dp) + real(2*i, dp), &
            i=0, nm%M-1 )]
        if (present(probabilities)) then
            nm%probabilities = probabilities
        else
            nm%probabilities(:) = 1/real(nm%M, dp)
        end if

        nm%E_symbol = 0
        do i = 0, nm%M-1
            nm%E_symbol = nm%E_symbol + nm%probabilities(i)*(abs(nm%constellation(i))**2)
        end do

        do i = 0, nm%M-1
            do k = 0, bps-1
                nm%symbol_to_bit_map(i, k) = iand(ishft(ishft(i, -k)+1, -1), 1) == 1
            end do
        end do

        allocate(nm%y_thresholds(1:nm%M-1))
        allocate(nm%Fy_thresholds(0:nm%M))
        ! allocate(nm%Fy_threshold_index(0:nm%M))
        allocate(nm%deltaFy(0:nm%M-1))

        call nm%update_N0(N0)
        ! this will also update the grid and set the thresholds

        allocate(nm%negative_monotonicity(0:nm%M-1))
        nm%negative_monotonicity(0::2) = .false.
        nm%negative_monotonicity(1::2) = .true.
        ! default to alternating configuration
    end function TNoiseMapperConstructor


    impure elemental function random_symbols(this) result(x)
        !! Return a random symbol index of the constellation
        class(TNoiseMapper), intent(in) :: this
        !! TNoiseMapper object to use for probabilities and constellation
        integer :: x
        !! Drawn result as index of constellation symbol: \(i \in \{0,\ldots,M-1\}\)

        double precision :: rnd
        integer :: i
        call random_number(rnd)

        x = this%M-1
        do i = 0, this%M-2
            rnd = rnd - this%probabilities(i)
            if (rnd .lt. 0) then
                x = i
                return
            end if
        end do
    end function random_symbols


    elemental function symbol_index_to_value(this, x_index) result (x_value)
        !! Convert constellation symbol index to symbol constellation value
        class(TNoiseMapper), intent(in) :: this
        !! Noise Mapper
        integer, intent(in) :: x_index
        !! constellation symbol index
        double precision :: x_value
        !! constellation point corresponding to index `x_index` (`x_value` = \(a_\texttt{x_index}\)

        x_value = this%constellation(x_index)
    end function symbol_index_to_value


    function symbol_to_word_single(this, symbol) result(word)
        !! De-map a symbol
        class(TNoiseMapper), intent(in) :: this
        !! Noise Mapper
        integer, intent(in) :: symbol
        !! index of the constellation symbol
        logical :: word(0:this%bps-1)
        !! bit group associated to the symbol

        word = this%symbol_to_bit_map(symbol, :)
    end function symbol_to_word_single


    function symbol_to_word_array(this, symbol) result(word)
        !! De-map a symbol
        class(TNoiseMapper), intent(in) :: this
        !! Noise Mapper
        integer, intent(in) :: symbol(0:)
        !! array of indexes of the constellation symbols
        logical :: word(0:size(symbol)*this%bps-1)
        !! bit word associated to the symbol array

        integer :: i

        do i = 0, size(symbol)-1
            word(i*this%bps : (i+1)*this%bps - 1) = this%symbol_to_bit_map(symbol(i), :)
        end do
    end function symbol_to_word_array


    subroutine update_N0(this, N0)
        !! Update total noise variance
        class(TNoiseMapper), intent(inout) :: this
        !! Noise mapper object
        double precision, intent(in) :: N0
        !! Total noise variance on both quadratures

        this%N0 = max(N0, 0d0)
        this%sigma = sqrt(this%N0/2d0)
        call this%set_grid
        call this%set_y_thresholds(5d-1*(this%constellation(1:this%M-1)+this%constellation(0:this%M-2)))
        ! Always call set_y_tresholds after set_grid

        call this%update_transition_probabilities
    end subroutine update_N0


    function y_to_lappr_single(this, y) result(lappr)
        !! Convert single channel output to lappr
        class(TNoiseMapper), intent(in) :: this
        !! nosie mapper
        double precision, intent(in) :: y
        !! channel sample
        double precision :: lappr(0:this%bps-1)
        !! array of LAPPRs

        double precision :: den(0:this%bps-1)
        double precision :: addendum
        integer :: i, k

        lappr = 0
        den   = 0

        do i = 0, this%M-1
            addendum = this%probabilities(i) * exp(-(y-this%constellation(i))**2/this%N0)
            do k = 0, this%bps-1
                if (this%symbol_to_bit_map(i, k)) then
                    den(k) = den(k) +  addendum
                else
                    lappr(k) = lappr(k) + addendum
                end if
            end do
        end do

        do k = 0, this%bps-1
            if (lappr(k) == 0) then
                lappr(k) = -1d100
            elseif (den(k) == 0) then
                lappr(k) = 1d100
            else
                lappr(k) = log(lappr(k)) - log(den(k))
            end if
        end do
    end function y_to_lappr_single


    function y_to_lappr_array(this, y) result(lappr)
        !! Convert array of channel outputs to array of LAPPRs
        class(TNoiseMapper), intent(in) :: this
        !! nosie mapper
        double precision, intent(in) :: y(0:)
        !! channel samples
        double precision :: lappr(0:size(y)*this%bps-1)
        !! array of LAPPRs

        integer :: i

        do i = 0, size(y)-1
            lappr(i*this%bps : (i+1)*this%bps - 1) = this%y_to_lappr_single(y(i))
        end do
    end function y_to_lappr_array


    subroutine set_y_thresholds(this, y_thresholds)
        !! set decision thresholds
        !! @warning No check on the ordering of the thresholds
        ! !! @warning Call after set_grid, because it needs Fy_threshold_index initialized
        class(TNoiseMapper), intent(inout) :: this
        !! noise mapper
        double precision, intent(in) :: y_thresholds(1:this%M-1)
        !! Set of thresholds

        double precision :: scale
        integer          :: i

        scale = sqrt(this%N0/2d0)

        this%y_thresholds = y_thresholds
        this%Fy_thresholds(1:this%M-1) = 0
        ! do i = 0, this%M-1
        !     this%Fy_thresholds(1:this%M-1) = this%Fy_thresholds(1:this%M-1) + &
        !         this%probabilities(i) * cdf_normal(&
        !         x     = y_thresholds, &
        !         loc   = this%constellation(i), &
        !         scale = scale)
        ! end do
        this%Fy_thresholds(1:this%M-1) = this%cdf_y(y_thresholds)
        this%Fy_thresholds(0) = 0
        this%Fy_thresholds(this%M) = 1
        this%deltaFy = this%Fy_thresholds(1:) - this%Fy_thresholds(:this%M-1)

        ! do i = 1, this%M-1
        !     this%Fy_threshold_index(i) = binsearch(this%F_grid, this%Fy_thresholds(i))
        ! end do
        ! this%Fy_threshold_index(0) = 1
        ! this%Fy_threshold_index(this%M) = size(this%F_grid) - 1
    end subroutine set_y_thresholds


    subroutine set_grid(this, threshold)
        !! Update grid points
        class(TNoiseMapper), intent(inout) :: this
        !! Noise mapper
        double precision, intent(in), optional :: threshold
        !! Threshold to consider the distribution as null

        double precision :: y_h, y_l, th
        integer          :: npoints, i

        if ( present(threshold) .and. (threshold .lt. 1) .and. (threshold .gt. 0)) then
            th = threshold
        else
            th = 1d-6
        end if


        y_h = this%constellation(this%M-1) + sqrt(-this%N0 * (log(th) - log(real(this%M, dp))))
        y_l = -y_h
        npoints = ceiling((y_h - y_l)*1000) + 1 ! 2000 points per step, with step = 2 between consecutive constellation points

        if (allocated(this%F_grid)) deallocate(this%F_grid)
        allocate(this%F_grid(0:npoints-1))
        if (allocated(this%y_grid)) deallocate(this%y_grid)
        allocate(this%y_grid(0:npoints-1))

        this%y_grid = [(y_l + real(i, dp)*(y_h-y_l)/real(npoints-1, dp), i = 0, npoints-1)]
        this%F_grid(1:npoints-2) = this%cdf_y(this%y_grid(1:npoints-2))
        this%F_grid(0) = 0
        this%F_grid(npoints-1) = 1
    end subroutine set_grid


    elemental function decide_symbol(this, y) result (x_i)
        !! Take a hard decision on the given output symbol based on thresholds
        class(TNoiseMapper), intent(in) :: this
        !! Noise Mapper
        double precision, intent(in)    :: y
        !! Channel output sample
        integer :: x_i
        !! Index of decision region where y belongs to

        x_i = binsearch(this%y_thresholds, y)
    end function decide_symbol


    elemental function cdf_y(this, y) result (Fy)
        !! Evaluate the cumulative density function of a channel output symbol
        class(TNoiseMapper), intent(in) :: this
        !! Noise mapper
        double precision,   intent(in) :: y
        !! channel output sample
        double precision               :: Fy
        !! CDF of the channel output sample

        Fy = sum(this%probabilities * &
            cdf_normal(x=y, loc=this%constellation, scale=sqrt(this%N0/2d0)))
    end function cdf_y


    elemental subroutine generate_soft_metric(this, y, nhat, xhat_i)
        !! Generate the softening metric from output sample
        class(TNoiseMapper), intent(in) :: this
        !! Noise mapper
        double precision, intent(in)    :: y
        !! Channel output sample
        double precision, intent(out)   :: nhat
        !! the output soft metric
        integer, intent(out), optional  :: xhat_i
        !! The index of the constellation symbol whose decision region `y` belongs to

        integer :: i

        i = binsearch(this%y_thresholds, y)
        if (this%negative_monotonicity(i)) then
            nhat = (this%Fy_thresholds(i+1) - this%cdf_y(y))/this%deltaFy(i)
        else
            nhat = (this%cdf_y(y) - this%Fy_thresholds(i))/this%deltaFy(i)
        end if

        if (present(xhat_i)) then
            xhat_i = i
        end if
    end subroutine generate_soft_metric


    pure function generate_tentative_channel_sample(this, nhat, xhat_i) result(y)
        !! Demap the soft metric `nhat` to a channel output `y` in the hypothesis that `xhat_i` was received
        class(TNoiseMapper), intent(in) :: this
        !! Noise Mapper
        double precision, intent(in)    :: nhat
        !! soft metric generated by Bob
        integer, intent(in)             :: xhat_i
        !! hypotetically received symbol
        double precision                :: y
        !! channel output whose metric corresponds to `nhat` and whose decision to `xhat_i`

        double precision :: Fy
        integer          :: i

        if (this%negative_monotonicity(xhat_i)) then
            Fy = this%Fy_thresholds(xhat_i + 1) - nhat * this%deltaFy(xhat_i)
        else
            Fy = this%Fy_thresholds(xhat_i)     + nhat * this%deltaFy(xhat_i)
        end if

        ! i = this%Fy_threshold_index(xhat_i) + binsearch(&
        !     this%F_grid(this%Fy_threshold_index(xhat_i) : this%Fy_threshold_index(xhat_i + 1)), &
        !     Fy)

        i = binsearch(this%F_grid(1:), Fy)
        ! slower than the commented code above, but easier to deal with
        ! consider including and testing the indexes of the thresholds
        ! within the grid for a faster search

        y = this%y_grid(i) + (this%y_grid(i+1) - this%y_grid(i)) * (Fy - this%F_grid(i)) / (this%F_grid(i+1) - this%F_grid(i))
    end function generate_tentative_channel_sample


    pure subroutine generate_lappr_single(this, x_i, nhat, lappr)
        !! Use the transmitted symbol and the soft metric to generate the LAPPRs
        class(TNoiseMapper), intent(in) :: this
        !! Noise Mapper
        integer, intent(in) :: x_i
        !! Index of transmitted symbol within the alphabet
        double precision, intent(in) :: nhat
        !! Soft metric generated by Bob
        double precision, intent(out) :: lappr(0:this%bps-1)
        !! LAPPR array of the bits associated with a single transmission

        double precision :: D(0:this%bps-1)
        double precision :: addendum
        double precision :: x ! Transmitted Symbol
        double precision :: yhat ! hypotetical channel output

        integer :: i ! Possible Index of the received symbol within the alphabet
        integer :: l ! Bit index within symbol
        integer :: k ! Generic constellation index

        D(:)     = 0
        lappr(:) = 0

        x        = this%constellation(x_i)

        do i = 0, this%M-1
            ! Reconstruct the channel output in the hypothesis the received symbol index is `i`
            yhat = this%generate_tentative_channel_sample(nhat, i)

            ! Build the term to be added to the numerator or denominator of the lappr
            addendum = 0
            do k = 0, this%M-1
                addendum = addendum + this%probabilities(k) * &
                    exp((2*yhat-this%constellation(k)-x)*(this%constellation(k) - x)/this%N0)
            end do
            addendum = this%deltaFy(i) / addendum

            ! Add the newly calculated lappr component to numerator or denominator
            ! according to the position of the bit within the symbol
            do l = 0, this%bps-1
                if ( this%symbol_to_bit_map(i, l) ) then
                    D(l) = D(l) + addendum
                else
                    lappr(l) = lappr(l) + addendum
                end if
            end do
        end do

        do l = 0, this%bps-1
            if (D(l) == 0) then
                lappr(l) = 1d100
            elseif (lappr(l) == 0) then
                lappr(l) = -1d100
            else
                ! The cases above should never happen, but you never know
                lappr(l) = log(lappr(l)) - log(D(l))
            end if
        end do
    end subroutine generate_lappr_single


    pure subroutine generate_lappr_array(this, x_i, nhat, lappr)
        !! Apply the function above to the entire set of soft metric and transmitted symbol pairs
        class(TNoiseMapper), intent(in) :: this
        !! Noise Mapper
        integer, intent(in) :: x_i(0:)
        !! array indexes of the transmitted symbols, one per channel use
        double precision, intent(in) :: nhat(0:size(x_i)-1)
        !! array of soft metric data, one per channel use
        double precision, intent(out) :: lappr(0:this%bps*size(x_i)-1)
        !! Array of lapprs, `this%bps` elements per channel use

        integer :: j

        do j = 0, size(nhat)-1
            call this%generate_lappr_single(x_i(j), nhat(j), lappr( j*this%bps : (j+1)*this%bps - 1 ) )
        end do
    end subroutine generate_lappr_array


    subroutine update_transition_probabilities(this)
        !! update the transition probability table between channel inputs and outputs
        class(TNoiseMapper), intent(inout) :: this
        !! Noise Mapper

        integer :: i, j

        double precision :: pXhat(0:this%M-1)

        if (.not. allocated(this%fwd_probability)) then
            allocate(this%fwd_probability(0:this%M-1, 0:this%M-1))
        end if
        if (.not. allocated(this%bwd_probability)) then
            allocate(this%bwd_probability(0:this%M-1, 0:this%M-1))
        end if
        if (.not. allocated(this%lappr_hard)) then
            allocate(this%lappr_hard(0:this%M-1, 0:this%bps-1))
        end if

        ! --------------------------------
        ! Forward transition probabilities
        ! --------------------------------
        do i = 0, this%M-1
            do j = 0, this%M-2
                this%fwd_probability(i,j) = cdf_normal(&
                    x=this%y_thresholds(j+1), &
                    loc=this%constellation(i),&
                    scale=this%sigma) ! set to CDF at upper bound of the region
            end do
            this%fwd_probability(i, this%M-1) = 1d0
        end do
        this%fwd_probability(:, 1:this%M-1) = &
            this%fwd_probability(:, 1:this%M-1) - this%fwd_probability(:, 0:this%M-2)
        ! remove lower bound, which is upper bound of previous region

        ! ---------------------------------
        ! Backward transition probabilities
        ! ---------------------------------
        do i = 0, this%M-1
            this%bwd_probability(i, :) = this%fwd_probability(i, :) * this%probabilities(i)
        end do
        pXhat = sum(this%bwd_probability, 1)
        do j = 0, this%M-1
            this%bwd_probability(:, j) = this%bwd_probability(:, j) / pXhat(j)
        end do

        ! -------------------------------------
        ! LAPPR for hard reverse reconciliation
        ! -------------------------------------
        do i = 0, this%M-1
            do j = 0, this%bps-1
                this%lappr_hard(i, j) = &
                    log(sum(this%fwd_probability(i, :), mask=.not.this%symbol_to_bit_map(:,j))) - &
                    log(sum(this%fwd_probability(i, :), mask=this%symbol_to_bit_map(:,j)))
            end do
        end do
    end subroutine update_transition_probabilities


    pure subroutine convert_symbol_to_hard_lappr_single(this, x, lappr)
        !! Generate LAPPR for hard reverse reconciliation from transmitted symbol
        class(TNoiseMapper), intent(in) :: this
        !! Noise mapper
        integer, intent(in) :: x
        !! transmitted symbol
        double precision, intent(out) :: lappr(0:this%bps-1)
        !! LAPPR

        lappr = this%lappr_hard(x, :)
    end subroutine convert_symbol_to_hard_lappr_single


    pure subroutine convert_symbol_to_hard_lappr_array(this, x, lappr)
        !! Generate LAPPR for hard reverse reconciliation
        !! from an array of transmitted symbols
        class(TNoiseMapper), intent(in) :: this
        !! Noise mapper
        integer, intent(in) :: x(0:)
        !! Array of transmitted symbols
        double precision, intent(out) :: lappr(0:size(x)*this%bps-1)
        !! array of LAPPRs associated to the input bits

        integer :: i

        do i = 0, size(x) -1
            lappr(i*this%bps : (i+1)*this%bps-1) = this%lappr_hard(x(i), :)
        end do
    end subroutine convert_symbol_to_hard_lappr_array

end module re2often_noise_mapper
