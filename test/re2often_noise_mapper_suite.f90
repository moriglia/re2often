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
    use stdlib_stats_distribution_normal, only: cdf_normal
    use re2often_noise_mapper
    use testdrive, only : new_unittest, unittest_type, error_type, check
    implicit none

    private

    public :: collect_suite

contains

    subroutine collect_suite(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [&
            new_unittest("Constructor", test_constructor),&
            new_unittest("Draw random symbols", test_random_symbol), &
            new_unittest("Convert symbol index to constellation point", test_symbol_index_to_value), &
            new_unittest("Convert symbol sequence to word", test_symbol_to_word), &
            new_unittest("Update N0", test_update_N0), &
            new_unittest("Test LAPPR construction for direct channel", test_y_to_lappr),&
            new_unittest("Test threshold setup and update", test_set_y_thresholds), &
            new_unittest("Test symbol decision", test_decide_symbol), &
            new_unittest("Test CDF of output channel", test_cdf_y), &
            new_unittest("Generation of soft metric", test_generate_soft_metric), &
            new_unittest("Reconstruction of channel output sample", test_generate_tentative_channel_sample), &
            new_unittest("LAPPR generation", test_generate_lappr)]
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


    subroutine test_random_symbol(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm
        integer :: x(100)

        nm = TNoiseMapper(2, 0.22d0, [0d0, 0d0, 1d0, 0d0])
        x = nm%random_symbols()
        call check(error, all(x==2))
        if (allocated(error)) return

        nm = TNoiseMapper(2, 0.22d0, [0d0, 0d0, 0d0, 1d0])
        x = nm%random_symbols()
        call check(error, all(x==3))
        if (allocated(error)) return

        nm = TNoiseMapper(2, 0.22d0, [1d0, 0d0, 0d0, 0d0])
        x = nm%random_symbols()
        call check(error, all(x==0))
        if (allocated(error)) return

        nm = TNoiseMapper(2, 0.22d0, [5d-1, 0d0, 2d-1, 3d-1])
        x = nm%random_symbols()
        call check(error, all(x/=1))
    end subroutine test_random_symbol


    subroutine test_symbol_index_to_value(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm
        integer :: x_i(100)
        double precision :: x(100), x_e(100)

        integer :: i

        nm = TNoiseMapper(3, 0.5d0)
        x_i = [(mod(i, 8), i = 1, 100)]
        x_e = -7d0 + real(x_i, dp)*2d0

        x = nm%symbol_index_to_value(x_i)
        call check(error, all(abs(x_e - x) < 1e-12))
    end subroutine test_symbol_index_to_value


    subroutine test_symbol_to_word(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm
        integer :: x(5)

        x = [0,3,3,1,2]

        nm = TNoiseMapper(bps=2, N0=1d0)
        call check(error, all(nm%symbol_to_word(x) .eqv. [&
            .false., .false., &
            .false., .true. , &
            .false., .true. , &
            .true. , .false., &
            .true. , .true.   ]))
    end subroutine test_symbol_to_word


    subroutine test_update_N0(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm

        nm = TNoiseMapper(2, 0.5d0)

        call check(error, nm%N0, 0.5d0, thr=1d-12)
        if (allocated(error)) return

        call nm%update_N0(0.75d0)
        call check(error, nm%N0, 0.75d0, thr=1d-12)
        if (allocated(error)) return

        call nm%update_N0(-0.15d0)
        call check(error, nm%N0, 0d0, thr=1d-12)
    end subroutine test_update_N0


    subroutine test_y_to_lappr(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm

        double precision :: y(10), lappr(20)

        nm = TNoiseMapper(1, 0.5d0)

        call random_number(y)
        ! even though they are not gaussian distributed,
        ! we are just testing the deterministic relation
        ! between channel output and lappr
        y = 4d0*y - 1

        lappr(:10) = nm%y_to_lappr(y)
        call check(error, all(abs(lappr(:10) + 4d0*y/0.5d0) < 1e-12))
        if (allocated(error)) return

        nm = TNoiseMapper(2, 0.5d0)
        lappr = nm%y_to_lappr(y)
        ! print *, lappr
        call check(error, all(abs(lappr(1::2) &
            + log(exp(-(y-1)**2/0.5d0) + exp(-(y+1)**2/0.5d0)) &
            - log(exp(-(y-3)**2/0.5d0) + exp(-(y+3)**2/0.5d0))) &
            .lt. 1d-12))
        if (allocated(error)) return

        call check(error, all(abs(lappr(2::2) &
            - log(exp(-(y+3)**2/0.5d0) + exp(-(y+1)**2/0.5d0)) &
            + log(exp(-(y-3)**2/0.5d0) + exp(-(y-1)**2/0.5d0))) &
            .lt. 1d-12))
    end subroutine test_y_to_lappr


    subroutine test_set_y_thresholds(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm

        nm = TNoiseMapper(bps=2, N0=1d0)

        call check(error, all(abs(nm%y_thresholds - real([-2, 0, 2], dp)) .lt. 1d-12))
        if (allocated(error)) return

        call nm%set_y_thresholds(real([-2.3, 0.15, 1.98], dp))

        call check(error, all(abs(nm%y_thresholds - real([-2.3, 0.15, 1.98], dp)) < 1d-12))
        if (allocated(error)) return
        call check(error, nm%Fy_thresholds(0) .lt. 1d-12)
        if (allocated(error)) return
        call check(error, nm%Fy_thresholds(nm%M) .gt. (1d0-1d-12))
        if (allocated(error)) return


        nm = TNoiseMapper(bps=1, N0=1d0)
        call check(error, nm%Fy_thresholds(1), 0.5d0, thr=1d-12)
        if (allocated(error)) return

        nm = TNoiseMapper(bps=1, N0=1d0, probabilities=real([0.1, 0.9], dp))
        call check(error, all(abs(nm%Fy_thresholds &
            - [0d0, &
            0.9d0 - 0.8d0*cdf_normal(x=1d0, loc=0d0, scale=sqrt(nm%N0/2d0)), &
            1d0]) .lt. 1d-9))
    end subroutine test_set_y_thresholds


    subroutine test_decide_symbol(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm

        nm = TNoiseMapper(bps=3, N0=1d0)
        call check(error, nm%decide_symbol(-9.3d0), 0)
        if (allocated(error)) return
        call check(error, nm%decide_symbol(-7.3d0), 0)
        if (allocated(error)) return
        call check(error, nm%decide_symbol(-6d0), 1)
        if (allocated(error)) return
        call check(error, nm%decide_symbol(-5.9d0), 1)
        if (allocated(error)) return
        call check(error, nm%decide_symbol(0d0), 4)
        if (allocated(error)) return
        call check(error, nm%decide_symbol(0.2d0), 4)
        if (allocated(error)) return
        call check(error, nm%decide_symbol(5.999d0), 6)
        if (allocated(error)) return
        call check(error, nm%decide_symbol(6.3d0), 7)
        if (allocated(error)) return

        call check(error, &
            all(nm%decide_symbol(real([-6.1, -2.3, 6.1], dp)) == [0, 2, 7]))
    end subroutine test_decide_symbol


    subroutine test_cdf_y(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm

        nm = TNoiseMapper(bps=1, N0=1d0, probabilities=real([0.1, 0.9], dp))
        call check(error, &
            nm%cdf_y(0.734d0), &
            0.9d0   * cdf_normal(x=.734d0, loc=1d0, scale=sqrt(nm%N0/2d0))   &
            + 0.1d0 * cdf_normal(x=.734d0, loc=-1d0, scale=sqrt(nm%N0/2d0)), &
            thr=1d-6)
        if (allocated(error)) return
        call check(error, &
            all(abs(nm%cdf_y([0.734d0, 0d0, 1d0]) - &
            [ 0.9d0   * cdf_normal(x=.734d0, loc=1d0, scale=sqrt(nm%N0/2d0))   &
            + 0.1d0 * cdf_normal(x=.734d0, loc=-1d0, scale=sqrt(nm%N0/2d0)),   &
            0.9d0 - 0.8d0*cdf_normal(x=1d0, loc=0d0, scale=sqrt(nm%N0/2d0)),   &
            0.9d0 * 0.5d0 + 0.1d0 * cdf_normal(x=1d0, loc=-1d0, scale=sqrt(nm%N0/2d0))])  &
            .lt. 1d-6))
    end subroutine test_cdf_y


    subroutine test_generate_soft_metric(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm

        double precision   :: y(8), nhat(8)
        integer            :: xhat_i(8), i

        nm = TNoiseMapper(bps=3, N0=1d0)

        y = real([-9.2, -5.2, -2.68, -1.12, 0.2, 2.14, 5.2, 6.47], dp)
        ! Each element is in  a different decision region


        ! Test elemental operation
        call nm%generate_soft_metric(y, nhat, xhat_i)

        ! Verify decisions
        do i = 1, 8
            call check(error, xhat_i(i), i-1)
            if (allocated(error)) return
        end do

        ! Verify nhat is in [0, 1]
        call check(error, all(nhat .le. 1))
        if (allocated(error)) return
        call check(error, all(nhat .ge. 0))
    end subroutine test_generate_soft_metric


    subroutine test_generate_tentative_channel_sample(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm

        double precision   :: y(8), nhat(8), yhat(8)
        integer            :: xhat_i(8), i

        nm = TNoiseMapper(bps=3, N0=1d0)

        y = real([-9.2, -5.2, -2.68, -1.12, 0.2, 2.14, 5.2, 6.47], dp)
        ! Each element is in  a different decision region

        ! Test elemental operation
        call nm%generate_soft_metric(y, nhat, xhat_i)

        ! Verify perfect reconstruction for correct symbols
        do i = 1, 8
            yhat(i) = nm%generate_tentative_channel_sample(nhat(i), xhat_i(i))
        end do

        call check(error, all(abs(yhat-y) .le. 1d-6))
        if (allocated(error)) return

        ! Add reconstruction tests for boundaries of decision regions
        call check(error, nm%generate_tentative_channel_sample(1d0, 0), nm%y_thresholds(1), thr=1d-6)
        if (allocated(error)) return
        do i = 1, 6
            ! In the following, the `iend()` addenda flip the monotonicity when needed
            ! The monotonicity is the default
            call check(error, nm%generate_tentative_channel_sample(0d0, i), &
                nm%y_thresholds(i + iand(i, 1)), thr=1d-6)
            if (allocated(error)) return
            call check(error, nm%generate_tentative_channel_sample(1d0, i), &
                nm%y_thresholds(i + 1 - iand(i, 1)), thr=1d-6)
            if (allocated(error)) return
        end do
        call check(error, nm%generate_tentative_channel_sample(1d0, 7), &
            nm%y_thresholds(7), thr=1d-6)
    end subroutine test_generate_tentative_channel_sample


    subroutine test_generate_lappr(error)
        type(error_type), allocatable, intent(out) :: error

        type(TNoiseMapper) :: nm

        double precision :: nhat
        double precision :: lappr_0(2)
        double precision :: lappr_1(2)
        double precision :: lappr_2(2)
        double precision :: lappr_3(2)

        nm = TNoiseMapper(bps=1, N0=0.5d0)
        nhat = 0.2 ! either y very low or very high
        call nm%generate_lappr(0, nhat, lappr_0)
        call nm%generate_lappr(1, nhat, lappr_1)
        call check(error, abs(lappr_0(1) + lappr_1(1)) .lt. 1d-9)
        if (allocated(error)) return
        call check(error, lappr_0(1) .gt. 0)
        if (allocated(error)) return
        call check(error, lappr_1(1) .lt. 0)
        if (allocated(error)) return

        nhat = 0.7
        call nm%generate_lappr(0, nhat, lappr_2)
        call nm%generate_lappr(1, nhat, lappr_3)
        call check(error, abs(lappr_2(1) + lappr_3(1)) .lt. 1d-9)
        if (allocated(error)) return
        call check(error, lappr_2(1) .gt. 0)
        if (allocated(error)) return
        call check(error, lappr_3(1) .lt. 0)
        if (allocated(error)) return

        ! Verify that the lappr values with nhat=0.7 have a lower
        ! absolute value then with nhat=0.2
        call check(error, lappr_0(1) .gt. lappr_2(1))
        if (allocated(error)) return
        call check(error, lappr_1(1) .lt. lappr_3(1))
        if (allocated(error)) return

        nm = TNoiseMapper(bps=2, N0=1d0)
        nhat=0.5
        call nm%generate_lappr(0, nhat, lappr_0)
        call nm%generate_lappr(1, nhat, lappr_1)
        call nm%generate_lappr(2, nhat, lappr_2)

        call check(error, all(lappr_0 .gt. 0))
        if (allocated(error)) return
        call check(error, (lappr_1(1) .lt. 0) .and. (lappr_1(2) .gt. 0))
        if (allocated(error)) return
        call check(error, all(lappr_2 .lt. 0))
        if (allocated(error)) return

        nhat=0.66
        call nm%generate_lappr(0, nhat, lappr_0)
        call nm%generate_lappr(1, nhat, lappr_1)
        call nm%generate_lappr(2, nhat, lappr_2)
        call nm%generate_lappr(3, nhat, lappr_3)

        ! Checks on one bit
        call check(error, abs(lappr_0(2) + lappr_3(2)) .lt. 1d-9)
        if (allocated(error)) return
        call check(error, abs(lappr_1(2) + lappr_2(2)) .lt. 1d-9)
        if (allocated(error)) return

        call check(error, abs(lappr_0(1) - lappr_3(1)) .lt. 1d-9)
        if (allocated(error)) return
        call check(error, abs(lappr_1(1) - lappr_2(1)) .lt. 1d-9)
    end subroutine test_generate_lappr
end module re2often_noise_mapper_suite
