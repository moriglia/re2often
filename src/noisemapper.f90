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
end module noisemapper
