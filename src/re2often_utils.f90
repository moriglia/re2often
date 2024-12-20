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
module re2often_utils
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    implicit none

    public :: binsearch

contains
    pure function binsearch(vector, y) result (index)
        !! Retrieve index of datum within a vector
        double precision, intent(in) :: vector(:)
        !! 1-based array \(\mathbf{v}\) to search data from
        double precision, intent(in) :: y
        !! Datum \(y\) to search
        integer :: index
        !! `index` \(i\) such that \(\mathbf{v}_i \leq y \lt \mathbf{v}_{i+1}\).
        !!  `index` is set to 0 if \(y \lt \mathbf{v}_{1}\).
        !!  `index` is set to `size(vector)` if \(y \geq \mathbf{v}_{\texttt{size(vector)}}\)

        integer :: index_u
        integer :: index_l

        if (y .lt. vector(1)) then
            index = 0
            return
        end if
        if (y .ge. vector(size(vector))) then
            index = size(vector)
            return
        end if

        index_u   = size(vector)
        index_l   = 1
        index = ishft(index_u + index_l, -1)
        do while ( (index_u - index_l) .gt. 1)
            if (y .ge. vector(index)) then
                index_l   = index
            else
                index_u   = index
            end if
            index = ishft(index_u + index_l, -1)
        end do
    end function binsearch
end module re2often_utils
