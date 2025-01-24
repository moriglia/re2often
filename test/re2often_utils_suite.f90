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
module re2often_utils_suite
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Test suite for the utilities
    use iso_fortran_env, only: dp => real64
    use stdlib_stats_distribution_normal, only: cdf_normal
    use re2often_utils
    use testdrive, only : new_unittest, unittest_type, error_type, check
    implicit none

    private

    public :: collect_suite

contains

    subroutine collect_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [&
            new_unittest("Binary search", test_binsearch),&
            new_unittest("Logical to integer conversion", test_logical2integer)]
    end subroutine collect_suite


    subroutine test_binsearch(error)
        type(error_type), intent(out), allocatable :: error

        real(dp) :: arr1(1)
        real(dp) :: arr2(2)
        real(dp) :: arr31(31)
        integer  :: i

        arr1(1) = 3d0
        call check(error, binsearch(arr1, 2d0), 0)
        if (allocated(error)) return
        call check(error, binsearch(arr1, 3d0), 1)
        if (allocated(error)) return
        call check(error, binsearch(arr1, 3.1d0), 1)
        if (allocated(error)) return

        arr2 = real([-2, 4], dp)
        call check(error, binsearch(arr2, -3d0), 0)
        if (allocated(error)) return
        call check(error, binsearch(arr2, real(-2, dp)), 1)
        if (allocated(error)) return
        call check(error, binsearch(arr2, 3.1d0), 1)
        if (allocated(error)) return
        call check(error, binsearch(arr2, real(4, dp)), 2)
        if (allocated(error)) return
        call check(error, binsearch(arr2, 9d0), 2)
        if (allocated(error)) return

        arr31 = [(real(i, dp), i = 1, 31)]
        call check(error, binsearch(arr31, -3d0), 0)
        if (allocated(error)) return
        call check(error, binsearch(arr31, real(1, dp)), 1)
        if (allocated(error)) return
        call check(error, binsearch(arr31, 3.1d0), 3)
        if (allocated(error)) return
        call check(error, binsearch(arr31, real(4, dp)), 4)
        if (allocated(error)) return
        call check(error, binsearch(arr31, 27.5d0), 27)
        if (allocated(error)) return
        call check(error, binsearch(arr31, 32d0), 31)
    end subroutine test_binsearch


    subroutine test_logical2integer(error)
        type(error_type), intent(out), allocatable :: error

        logical :: arr1(0:3)
        logical :: arr2(0:5)

        arr1 = [.false., .false., .false., .true.]
        call check(error, logical2integer(arr1) == 8)
        if (allocated(error)) return

        arr1 = [.false., .true., .true., .false.]
        call check(error, logical2integer(arr1) == 6)
        if (allocated(error)) return

        arr2 = [.true., .true., .false., .false., .false., .false.]
        call check(error, logical2integer(arr2) == 3)
        if (allocated(error)) return

        arr2 = [.false., .true., .false., .false., .true., .false.]
        call check(error, logical2integer(arr2) == 18)
        if (allocated(error)) return

        arr2 = [.false., .true., .false., .false., .true., .true.]
        call check(error, logical2integer(arr2) == 50)
    end subroutine test_logical2integer

end module re2often_utils_suite
