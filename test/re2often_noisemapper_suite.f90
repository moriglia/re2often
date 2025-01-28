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
module re2often_noisemapper_suite
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Test suite for the Noise mapper
    use iso_fortran_env, only: dp => real64
    use stdlib_stats_distribution_normal, only: cdf_normal
    use noisemapper
    use testdrive, only : new_unittest, unittest_type, error_type, check
    implicit none

    private

    public :: collect_suite

contains

    subroutine collect_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [&
            new_unittest("Constructor", test_constructor) &
        ]
    end subroutine collect_suite


    subroutine test_constructor(error)
        type(error_type), allocatable, intent(out) :: error

        type(noisemapper_type) :: nm

        nm = noisemapper_create(3)

        call check(error, nm%bps, 3)
        if (allocated(error)) return

        call check(error, nm%M, 8)
        if (allocated(error)) return

        call check(error, maxval(abs(nm%probabilities - 0.125d0)) < 1d-12)
        if (allocated(error)) return

        call check(error, &
            maxval(abs(nm%constellation - real([-7, -5, -3, -1, 1, 3, 5, 7], dp))) < 1d-12)
        if (allocated(error)) return

        call check(error, nm%E_s, 21d0, thr=1d-12)
        if (allocated(error)) return

        call check(error, all(nm%s_to_b(:, 0) .eqv. &
            [.false., .true., .true., .false., .false., .true., .true., .false.]))
        if (allocated(error)) return
        call check(error, all(nm%s_to_b(:, 1) .eqv. &
            [.false., .false., .true., .true., .true., .true., .false., .false.]))
        if (allocated(error)) return
        call check(error, all(nm%s_to_b(:, 2) .eqv. &
            [.false., .false., .false., .false., .true., .true., .true., .true.]))
        if (allocated(error)) return


        ! Test correct destruction/construction
        call noisemapper_deallocate(nm)
        nm = noisemapper_create(2)

        call check(error, nm%bps, 2)
        if (allocated(error)) return

        call check(error, nm%M, 4)
        if (allocated(error)) return

        call check(error, maxval(abs(nm%probabilities - 0.25d0)) < 1d-12)
        if (allocated(error)) return

        call check(error, &
            maxval(abs(nm%constellation - real([-3, -1, 1, 3], dp))) < 1d-12)
        if (allocated(error)) return

        call check(error, nm%E_s, 5d0, thr=1d-12)
        if (allocated(error)) return

        call check(error, all(nm%s_to_b(:, 0) .eqv. &
            [.false., .true., .true., .false.]))
        if (allocated(error)) return
        call check(error, all(nm%s_to_b(:, 1) .eqv. &
            [.false., .false., .true., .true.]))
    end subroutine test_constructor


    subroutine test_update_N0_from_snrdb(error)
        type(error_type), allocatable, intent(out) :: error

        type(noisemapper_type) :: nm

        nm = noisemapper_create(3)

        call noisemapper_update_N0_from_snrdb(nm, 0d0)

        call check(error, nm%N0, 21d0, thr=1d-12)
        if (allocated(error)) return
        call check(error, nm%sigma, sqrt(10.5d0), thr=1d-12)
        if (allocated(error)) return

        call noisemapper_update_N0_from_snrdb(nm, 10d0*log10(2d0))

        call check(error, nm%N0, 42d0, thr=1d-12)
        if (allocated(error)) return
        call check(error, nm%sigma, sqrt(21d0), thr=1d-12)
    end subroutine test_update_N0_from_snrdb

end module re2often_noisemapper_suite
