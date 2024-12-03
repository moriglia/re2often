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
module ldpc_edge_list
  implicit none

  public :: TEdgeList
  
  type, public :: TEdgeList
     integer :: N
     integer, allocatable :: data(:)
   contains
     final :: destructor
  end type TEdgeList

  interface TEdgeList
     module procedure TEdgeListConstructor
  end interface TEdgeList
contains
  function TEdgeListConstructor(N) result (buffer)
    integer, intent(in) :: N
    type(TEdgeList) :: buffer

    buffer%N=N
    allocate(buffer%data(N))
  end function TEdgeListConstructor

  subroutine destructor(buffer)
    type(TEdgeList) :: buffer
    
    if (allocated(buffer%data)) deallocate(buffer%data)
  end subroutine destructor
  
end module ldpc_edge_list
