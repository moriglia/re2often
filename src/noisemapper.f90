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
module noisemapper
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Implementation of alphabet utilities and noise mapping/demapping
    !! Mainteining C interoperability
    use, intrinsic :: iso_c_binding
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
    end type noisemapper_type


    interface noisemapper_y_to_lappr
        module procedure noisemapper_y_to_lappr_single
        module procedure noisemapper_y_to_lappr_array
    end interface noisemapper_y_to_lappr

    interface noisemapper_random_symbol
        module procedure noisemapper_random_symbol_single
        module procedure noisemapper_random_symbol_array
    end interface noisemapper_random_symbol

contains


    subroutine noisemapper_deallocate(nm)
        !! Destructor for noise mapper
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper

        if (allocated(nm%constellation)) deallocate(nm%constellation)
        if (allocated(nm%probabilities)) deallocate(nm%probabilities)
        if (allocated(nm%s_to_b       )) deallocate(nm%s_to_b       )
    end subroutine noisemapper_deallocate


    function noisemapper_create(bps) result(nm)
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


    subroutine noisemapper_update_N0_from_snrdb(nm, snrdb)
        !! Update N0 based on the value of the SNR
        type(noisemapper_type), intent(inout) :: nm
        !! Noise mapper
        real(c_double), intent(in) :: snrdb
        !! SNR in dB

        nm%N0 = nm%E_s * (10d0**(-snrdb/10d0))
        nm%sigma = sqrt(nm%N0/2d0)
    end subroutine noisemapper_update_N0_from_snrdb


    subroutine noisemapper_y_to_lappr_single(nm, y, lappr)
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
            addendum = nm%probabilities(i) * exp(-(y - nm%constellation(i))/nm%N0)
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


    subroutine noisemapper_y_to_lappr_array(nm, y, lappr)
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


    subroutine noisemapper_random_symbol_single(nm, x_i)
        !! Generate a random symbol
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(out) :: x_i
        !! Random symbol of the constellation (index in 0:M-1)

        double precision :: rnd

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


    subroutine noisemapper_random_symbol_array(nm, x_i)
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


    function noisemapper_symbol_index_to_value(nm, x_i) result (x)
        !! Convert constellation index to point
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        integer(c_int), intent(in) :: x_i(:)
        !! Set of constellation indexes
        integer(c_double) :: x(size(x_i))
        !! set of constellation points

        integer :: j

        do j = 1, size(x_i)
            x = nm%constellation(x_i(j))
        end do
    end function noisemapper_symbol_index_to_value

end module noisemapper
