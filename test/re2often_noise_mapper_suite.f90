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
module re2often_noise_mapper_suite
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Test suite for the Noise mapper
    use iso_fortran_env, only: dp => real64
    use re2often_noise_mapper
    use testdrive, only : new_unittest, unittest_type, error_type, check
    implicit none

    private

    public :: collect_suite

contains

    subroutine collect_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [new_unittest("Constructor", test_constructor)]
    end subroutine collect_suite


    subroutine test_constructor(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm


        nm = TNoiseMapper(bps=3, N0=0.5d0)

        call check(error, nm%bps, 3)
        if (allocated(error)) return

        call check(error, nm%M, 8)
        if (allocated(error)) return

        call check(error, nm%N0, 0.5d0, thr=1d-12)
        if (allocated(error)) return

        call check(error, maxval(abs(nm%probabilities - 0.125d0)) < 1d-12)
        if (allocated(error)) return

        call check(error, &
            maxval(abs(nm%constellation - real([-7, -5, -3, -1, 1, 3, 5, 7], dp))) < 1d-12)
        if (allocated(error)) return

        call check(error, nm%E_symbol, 21d0, thr=1d-12)
        if (allocated(error)) return

        call check(error, all(nm%symbol_to_bit_map(:, 0) .eqv. &
            [.false., .true., .true., .false., .false., .true., .true., .false.]))
        if (allocated(error)) return
        call check(error, all(nm%symbol_to_bit_map(:, 1) .eqv. &
            [.false., .false., .true., .true., .true., .true., .false., .false.]))
        if (allocated(error)) return
        call check(error, all(nm%symbol_to_bit_map(:, 2) .eqv. &
            [.false., .false., .false., .false., .true., .true., .true., .true.]))
        if (allocated(error)) return


        ! Test correct destruction/construction
        nm = TNoiseMapper(2, 0.75d0, [0.15d0, 0.35d0, 0.3d0, 0.2d0])

        call check(error, nm%bps, 2)
        if (allocated(error)) return

        call check(error, nm%M, 4)
        if (allocated(error)) return

        call check(error, nm%N0, 0.75d0, thr=1d-12)
        if (allocated(error)) return

        call check(error, maxval(abs(nm%probabilities - [0.15d0, 0.35d0, 0.3d0, 0.2d0])) < 1d-12)
        if (allocated(error)) return

        call check(error, &
            maxval(abs(nm%constellation - real([-3, -1, 1, 3], dp))) < 1d-12)
        if (allocated(error)) return

        call check(error, nm%E_symbol, 9d0*0.35d0 + 0.65d0, thr=1d-12)
        if (allocated(error)) return

        call check(error, all(nm%symbol_to_bit_map(:, 0) .eqv. &
            [.false., .true., .true., .false.]))
        if (allocated(error)) return
        call check(error, all(nm%symbol_to_bit_map(:, 1) .eqv. &
            [.false., .false., .true., .true.]))
    end subroutine test_constructor

end module re2often_noise_mapper_suite
