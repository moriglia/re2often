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
program configurations_list
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! List one representative for each configurations class
    implicit none

    integer :: argc, ii, bps, M, other
    character(len=500), allocatable :: argv(:)
    logical :: vFlag

    argc = command_argument_count()
    allocate(argv(argc))

    do ii = 1, argc
        call get_command_argument(ii, argv(ii))
    end do

    bps = 2
    M = 4
    vFlag = .false.

    ii = 1
    do while (ii .le. argc)
        if (argv(ii) == "--bps") then
            ii = ii + 1
            read(argv(ii),*) bps
        elseif (argv(ii) == "--verbose") then
            print *, "Verbose output"
            vFlag = .true.
        else
            print *, "Unrecognized argument"
            stop
        end if
        ii = ii + 1
    end do


    M = ishft(1, bps)
    cfg_loop: do ii = 0, ishft(1, M) - 1

        other = flip(ii, M)
        if (other .lt. ii) then
            call display_verbose(vFlag, "Flipped", other, ii)
            goto 100
        end if

        other = reverse(ii, M)
        if (other .lt. ii) then
            call display_verbose(vFlag, "Reverse", other, ii)
            goto 100
        end if

        other = mirror(ii, M)
        if (other .lt. ii) then
            call display_verbose(vFlag, "Mirrord", other, ii)
            goto 100
        end if

        print '(I0)', ii

100     continue
    end do cfg_loop

contains
    integer function flip(ii, M) result(other)
        integer, intent(in) :: ii
        integer, intent(in) :: M

        other = ishft(1, M) - 1  ! all 1s
        other = ieor(other, ii) ! bitwise xor
    end function flip


    integer function reverse(ii, M) result(other)
        integer, intent(in) :: ii
        integer, intent(in) :: M

        integer :: pos

        other = 0

        do pos = 0, M-1
            other = ior(other, ishft(iand(ishft(ii, -(M-1-pos)), 1), pos))
        end do
    end function reverse


    integer function mirror(ii, M) result(other)
        integer, intent(in) :: ii
        integer, intent(in) :: M

        other = flip(reverse(ii, M), M)
    end function mirror


    subroutine display_verbose(flag, a, b, c)
        logical, intent(in) :: flag
        character(len=*), intent(in) :: a
        integer, intent(in) :: b
        integer, intent(in) :: c

        if (flag) then
            print '(A, T10, I0, T20, I0)', a, b, c
        end if
    end subroutine display_verbose

end program configurations_list
