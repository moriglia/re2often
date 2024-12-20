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
    contains
        procedure, public, pass :: free_noise_mapper
        procedure, public, pass :: deallocate_thresholds
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
        procedure, public, pass :: decide_symbol
        procedure, public, pass :: cdf_y
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
    end subroutine deallocate_thresholds

    subroutine free_noise_mapper(nm)
        !! Deallocates all allocatable variables within TNoiseMapper type
        class(TNoiseMapper), intent(inout) :: nm
        !! The noise mapper whose allocatable variables are to be deallocated

        if (allocated(nm%constellation))     deallocate(nm%constellation)
        if (allocated(nm%probabilities))     deallocate(nm%probabilities)
        if (allocated(nm%symbol_to_bit_map)) deallocate(nm%symbol_to_bit_map)
        call nm%deallocate_thresholds
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
        call nm%update_N0(N0)

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
        call nm%set_y_thresholds(5d-1*(nm%constellation(1:nm%M-1)+nm%constellation(0:nm%M-2)))
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
        class(TNoiseMapper), intent(inout) :: this
        !! noise mapper
        double precision, intent(in) :: y_thresholds(1:this%M-1)
        !! Set of thresholds

        double precision :: scale
        integer          :: i

        scale = sqrt(this%N0/2d0)

        this%y_thresholds = y_thresholds
        this%Fy_thresholds(1:this%M-1) = 0
        do i = 0, this%M-1
            this%Fy_thresholds(1:this%M-1) = this%Fy_thresholds(1:this%M-1) + &
                this%probabilities(i) * cdf_normal(&
                x     = y_thresholds, &
                loc   = this%constellation(i), &
                scale = scale)
        end do
        this%Fy_thresholds(0) = 0
        this%Fy_thresholds(this%M) = 1
    end subroutine set_y_thresholds


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

end module re2often_noise_mapper
