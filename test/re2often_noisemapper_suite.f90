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
    use iso_c_binding
    use stdlib_stats_distribution_normal, only: cdf_normal
    use re2often_noisemapper
    use testdrive, only : new_unittest, unittest_type, error_type, check
    implicit none

    private

    public :: collect_suite

contains

    subroutine collect_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [&
            new_unittest("Constructor", test_constructor), &
            new_unittest("Update N0 from SNR", test_update_N0_from_snrdb), &
            new_unittest("Direct LAPPR", test_y_to_lappr), &
            new_unittest("Symbol index to value", test_symbol_index_to_value), &
            new_unittest("Symbol to word", test_symbol_to_word), &
            new_unittest("Set y thresholds with defaults", test_set_y_thresholds_default), &
            new_unittest("Symbol decision", test_decide_symbol), &
            new_unittest("Uniform thresholds", test_set_y_thresholds_uniform) &
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

        call noisemapper_update_N0_from_snrdb(nm, 10d0*log10(0.5d0))

        call check(error, nm%N0, 42d0, thr=1d-12)
        if (allocated(error)) return
        call check(error, nm%sigma, sqrt(21d0), thr=1d-12)
    end subroutine test_update_N0_from_snrdb


    subroutine test_y_to_lappr(error)
        type(error_type), allocatable, intent(out) :: error

        type(noisemapper_type) :: nm

        real(c_double) :: lappr(0:8)
        real(c_double) :: y(0:2)

        y(0) = 2.97d0

        nm = noisemapper_create(3)

        call noisemapper_update_N0_from_snrdb(nm, 0d0)
        call noisemapper_y_to_lappr(nm, y(0), lappr(0:2))

        call check(error, lappr(0), &
            log(exp(-(y(0)+7)**2/21d0) + exp(-(y(0)+1)**2/21d0) + &
            exp(-(y(0)-7)**2/21d0) + exp(-(y(0)-1)**2/21d0)) - &
            log(exp(-(y(0)+5)**2/21d0) + exp(-(y(0)+3)**2/21d0) + &
            exp(-(y(0)-5)**2/21d0) + exp(-(y(0)-3)**2/21d0)) , &
            thr=1d-12)
        if (allocated(error)) return

        call check(error, lappr(1), &
            log(exp(-(y(0)+7)**2/21d0) + exp(-(y(0)+5)**2/21d0) + &
            exp(-(y(0)-7)**2/21d0) + exp(-(y(0)-5)**2/21d0)) - &
            log(exp(-(y(0)+1)**2/21d0) + exp(-(y(0)+3)**2/21d0) + &
            exp(-(y(0)-1)**2/21d0) + exp(-(y(0)-3)**2/21d0)) , &
            thr=1d-12)
        if (allocated(error)) return

        call check(error, lappr(2), &
            log(exp(-(y(0)+7)**2/21d0) + exp(-(y(0)+5)**2/21d0) + &
            exp(-(y(0)+3)**2/21d0) + exp(-(y(0)+1)**2/21d0)) - &
            log(exp(-(y(0)-1)**2/21d0) + exp(-(y(0)-3)**2/21d0) + &
            exp(-(y(0)-5)**2/21d0) + exp(-(y(0)-7)**2/21d0)) , &
            thr=1d-12)
        if (allocated(error)) return

        y(1:2) = [-3.4d0, 9.1d0]

        ! Re-check y(0)
        call noisemapper_y_to_lappr(nm, y, lappr)
        call check(error, lappr(0), &
            log(exp(-(y(0)+7)**2/21d0) + exp(-(y(0)+1)**2/21d0) + &
            exp(-(y(0)-7)**2/21d0) + exp(-(y(0)-1)**2/21d0)) - &
            log(exp(-(y(0)+5)**2/21d0) + exp(-(y(0)+3)**2/21d0) + &
            exp(-(y(0)-5)**2/21d0) + exp(-(y(0)-3)**2/21d0)) , &
            thr=1d-12)
        if (allocated(error)) return

        call check(error, lappr(1), &
            log(exp(-(y(0)+7)**2/21d0) + exp(-(y(0)+5)**2/21d0) + &
            exp(-(y(0)-7)**2/21d0) + exp(-(y(0)-5)**2/21d0)) - &
            log(exp(-(y(0)+1)**2/21d0) + exp(-(y(0)+3)**2/21d0) + &
            exp(-(y(0)-1)**2/21d0) + exp(-(y(0)-3)**2/21d0)) , &
            thr=1d-12)
        if (allocated(error)) return

        call check(error, lappr(2), &
            log(exp(-(y(0)+7)**2/21d0) + exp(-(y(0)+5)**2/21d0) + &
            exp(-(y(0)+3)**2/21d0) + exp(-(y(0)+1)**2/21d0)) - &
            log(exp(-(y(0)-1)**2/21d0) + exp(-(y(0)-3)**2/21d0) + &
            exp(-(y(0)-5)**2/21d0) + exp(-(y(0)-7)**2/21d0)) , &
            thr=1d-12)
        if (allocated(error)) return


        ! check y(1)
        call check(error, lappr(3), &
            log(exp(-(y(1)+7)**2/21d0) + exp(-(y(1)+1)**2/21d0) + &
            exp(-(y(1)-7)**2/21d0) + exp(-(y(1)-1)**2/21d0)) - &
            log(exp(-(y(1)+5)**2/21d0) + exp(-(y(1)+3)**2/21d0) + &
            exp(-(y(1)-5)**2/21d0) + exp(-(y(1)-3)**2/21d0)) , &
            thr=1d-12)
        if (allocated(error)) return

        call check(error, lappr(4), &
            log(exp(-(y(1)+7)**2/21d0) + exp(-(y(1)+5)**2/21d0) + &
            exp(-(y(1)-7)**2/21d0) + exp(-(y(1)-5)**2/21d0)) - &
            log(exp(-(y(1)+1)**2/21d0) + exp(-(y(1)+3)**2/21d0) + &
            exp(-(y(1)-1)**2/21d0) + exp(-(y(1)-3)**2/21d0)) , &
            thr=1d-12)
        if (allocated(error)) return

        call check(error, lappr(5), &
            log(exp(-(y(1)+7)**2/21d0) + exp(-(y(1)+5)**2/21d0) + &
            exp(-(y(1)+3)**2/21d0) + exp(-(y(1)+1)**2/21d0)) - &
            log(exp(-(y(1)-1)**2/21d0) + exp(-(y(1)-3)**2/21d0) + &
            exp(-(y(1)-5)**2/21d0) + exp(-(y(1)-7)**2/21d0)) , &
            thr=1d-12)
        if (allocated(error)) return

        ! check y(2)
        call check(error, lappr(6), &
            log(exp(-(y(2)+7)**2/21d0) + exp(-(y(2)+1)**2/21d0) + &
            exp(-(y(2)-7)**2/21d0) + exp(-(y(2)-1)**2/21d0)) - &
            log(exp(-(y(2)+5)**2/21d0) + exp(-(y(2)+3)**2/21d0) + &
            exp(-(y(2)-5)**2/21d0) + exp(-(y(2)-3)**2/21d0)) , &
            thr=1d-12)
        if (allocated(error)) return

        call check(error, lappr(7), &
            log(exp(-(y(2)+7)**2/21d0) + exp(-(y(2)+5)**2/21d0) + &
            exp(-(y(2)-7)**2/21d0) + exp(-(y(2)-5)**2/21d0)) - &
            log(exp(-(y(2)+1)**2/21d0) + exp(-(y(2)+3)**2/21d0) + &
            exp(-(y(2)-1)**2/21d0) + exp(-(y(2)-3)**2/21d0)) , &
            thr=1d-12)
        if (allocated(error)) return

        call check(error, lappr(8), &
            log(exp(-(y(2)+7)**2/21d0) + exp(-(y(2)+5)**2/21d0) + &
            exp(-(y(2)+3)**2/21d0) + exp(-(y(2)+1)**2/21d0)) - &
            log(exp(-(y(2)-1)**2/21d0) + exp(-(y(2)-3)**2/21d0) + &
            exp(-(y(2)-5)**2/21d0) + exp(-(y(2)-7)**2/21d0)) , &
            thr=1d-12)
    end subroutine test_y_to_lappr


    subroutine test_symbol_index_to_value(error)
        type(error_type), allocatable, intent(out) :: error

        type(noisemapper_type) :: nm

        integer(c_int) :: x_i(10)
        real(c_double) :: x(10)

        nm = noisemapper_create(2)

        x_i = [0, 2, 3, 3, 1, 3, 1, 0, 2, 2]
        x   = real([-3, 1, 3, 3, -1, 3, -1, -3, 1, 1], c_double)

        call check(error, all(abs(noisemapper_symbol_index_to_value(nm, x_i) - x) .lt. 1d-12))
    end subroutine test_symbol_index_to_value


    subroutine test_symbol_to_word(error)
        type(error_type), allocatable, intent(out) :: error

        type(noisemapper_type) :: nm
        integer(c_int) :: x_i(10)
        logical(c_bool):: word(30)
        logical(c_bool):: word_expected(30)

        nm = noisemapper_create(2)
        x_i = [0, 2, 3, 1, 1, 0, 3, 1, 2, 2]
        word_expected(:20) = [&
            .false., .false.,&
            .true., .true.,  &
            .false., .true., &
            .true., .false., &
            .true., .false., &
            .false., .false.,&
            .false., .true., &
            .true., .false., &
            .true., .true.,  &
            .true., .true.   &
            ]
        word(:20) = noisemapper_symbol_to_word(nm, x_i)
        call check(error, logical(all(word(:20) .eqv. word_expected(:20))))
        if (allocated(error)) then
            print *, word(:20)
            print *, word_expected(:20)
            return
        end if

        nm = noisemapper_create(3) ! should call deallocate within
        x_i = [0, 2, 3, 1, 1, 0, 3, 1, 2, 2]
        word_expected = [&
            .false., .false., .false., &
            .true. , .true. , .false., &
            .false., .true. , .false., &
            .true. , .false., .false., &
            .true. , .false., .false., &
            .false., .false., .false., &
            .false., .true. , .false., &
            .true. , .false., .false., &
            .true. , .true. , .false., &
            .true. , .true. , .false.  &
            ]
        word = noisemapper_symbol_to_word(nm, x_i)
        call check(error, logical(all(word .eqv. word_expected)))
        if (allocated(error)) then
            print *, word
            print *, word_expected
            return
        end if

        x_i = x_i + 4
        word_expected(3::3) = .true.
        word_expected(2::3) = .not. word_expected(2::3)
        word = noisemapper_symbol_to_word(nm, x_i)
        call check(error, logical(all(word .eqv. word_expected)))
        if (allocated(error)) then
            print *, word
            print *, word_expected
            return
        end if
    end subroutine test_symbol_to_word

    ! +------------------------------------------+
    ! | REVERSE reconciliation common procedures |
    ! +------------------------------------------+
    subroutine test_set_y_thresholds_default(error)
        type(error_type), allocatable, intent(out) :: error

        type(noisemapper_type) :: nm

        nm = noisemapper_create(2)

        call noisemapper_update_N0_from_snrdb(nm, 0d0)
        call noisemapper_set_y_thresholds(nm)

        call check(error, all(abs(nm%y_thresholds - real([-2, 0, 2], c_double)) .lt. 1d-12))
        if (allocated(error)) return

        call check(error, all(abs(nm%Fy_thresholds(1:3) &
            - 0.25d0*cdf_normal(loc=-3d0, x=real([-2, 0, 2], c_double), scale=nm%sigma) &
            - 0.25d0*cdf_normal(loc=-1d0, x=real([-2, 0, 2], c_double), scale=nm%sigma) &
            - 0.25d0*cdf_normal(loc=+1d0, x=real([-2, 0, 2], c_double), scale=nm%sigma) &
            - 0.25d0*cdf_normal(loc=+3d0, x=real([-2, 0, 2], c_double), scale=nm%sigma) &
            ) .lt. 1d-12))
        if (allocated(error)) return

        call check(error, nm%Fy_thresholds(0), 0d0, thr=1d-12)
        if (allocated(error)) return

        call check(error, nm%Fy_thresholds(4), 1d0, thr=1d-12)
    end subroutine test_set_y_thresholds_default


    subroutine test_decide_symbol(error)
        type(error_type), allocatable, intent(out) :: error

        type(noisemapper_type) :: nm
        integer :: x_i(10)
        real(c_double) :: y(10)

        nm = noisemapper_create(3)
        call noisemapper_update_N0_from_snrdb(nm, 0d0)
        call noisemapper_set_y_thresholds(nm)

        y = real([-9.2d0, -4.4d0, 6.01d0, 7.2d0, 15.4d0, 0.2d0, 2.05d0, -1.992d0, 1.98d0, -5.976d0], c_double)
        x_i = noisemapper_decide_symbol(nm, y)
        call check(error, all(x_i .eq. [0, 1, 7, 7, 7, 4, 5, 3, 4, 1]))
    end subroutine test_decide_symbol


    ! +----------------------------------------+
    ! | HARD REVERSE reconciliation procedures |
    ! +----------------------------------------+

    ! +----------------------------------------+
    ! | SOFT REVERSE reconciliation procedures |
    ! +----------------------------------------+

    subroutine test_set_y_thresholds_uniform(error)
        type(error_type), allocatable, intent(out) :: error

        type(noisemapper_type) :: nm
        ! integer :: i

        nm = noisemapper_create(3)
        call noisemapper_update_N0_from_snrdb(nm, 0d0)
        call noisemapper_set_Fy_grids(nm)
        call noisemapper_set_y_thresholds_uniform(nm)

        call check(error, all(abs(nm%delta_Fy - 0.125d0) .lt. 1d-12))
        if (allocated(error)) then
            print *, nm%Fy_grid(1::1000)
            print *, "Thresholds: ", nm%y_thresholds
            print *, "CDF: ", nm%Fy_thresholds
            print *, "Deltas:", nm%delta_Fy
            return
        end if
    end subroutine test_set_y_thresholds_uniform
end module re2often_noisemapper_suite
