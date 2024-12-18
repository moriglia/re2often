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
        logical(1), allocatable :: symbol_to_bit_map(:,:)
        !! Bits (rows) associated to each symbol (row index), size is (0:M-1, 0:bps-1)
    contains
        procedure, public, pass :: free_noise_mapper
        final :: TNoiseMapperDestructor
    end type TNoiseMapper


    interface TNoiseMapper
        module procedure TNoiseMapperConstructor
    end interface TNoiseMapper
contains

    subroutine free_noise_mapper(nm)
        !! Deallocates all allocatable variables within TNoiseMapper type
        class(TNoiseMapper), intent(inout) :: nm
        !! The noise mapper whose allocatable variables are to be deallocated

        if (allocated(nm%constellation))     deallocate(nm%constellation)
        if (allocated(nm%probabilities))     deallocate(nm%probabilities)
        if (allocated(nm%symbol_to_bit_map)) deallocate(nm%symbol_to_bit_map)
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
        nm%N0 = N0

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
    end function TNoiseMapperConstructor

end module re2often_noise_mapper
