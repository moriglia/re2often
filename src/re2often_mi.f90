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
module re2often_mi
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Mutual information functions
    use, intrinsic :: iso_c_binding
    use re2often_noisemapper
    use external_hermite
    use quadpack, only: dqags, dqagi
    implicit none

    private
    public :: I_soft_reverse_equidistant_th, I_soft_reverse_uniform_output_th
    public :: I_hard_reverse_equidistant_th, I_hard_reverse_uniform_output_th
    public :: I_direct

    type(noisemapper_type) :: nm

    public :: nm

    real(c_double), parameter :: sq2 = sqrt(2d0)

contains

    real(c_double) elemental function log2(arg) result(l)
        real(c_double), intent(in) :: arg

        l = log(arg)/log(2d0)
    end function log2

    ! +-----------------------------+
    ! | Soft Reverse Reconciliation |
    ! +-----------------------------+

    real(c_double) impure elemental function f_n_xhat_cond_x(n, xhat, x) result(pdf)
        !! PDF of \(N, \hat{X}|X\)
        real(c_double), intent(in) :: n
        !! Soft metric
        integer(c_int), intent(in) :: xhat
        !! Bob's decided symbol (index within constellation)
        integer(c_int), intent(in) :: x
        !! Alice's transmitted symbol (index within constellation)

        integer :: k
        real(c_double) :: a_j, two_y_i

        pdf = 0

        a_j = nm%constellation(x)
        two_y_i = 2*noisemapper_invert_soft_metric(nm, n, xhat)
        do k = 0, nm%M-1
            pdf = pdf + nm%probabilities(k) * &
                exp((nm%constellation(k) - a_j)*(two_y_i - nm%constellation(k) - a_j)/nm%N0)
        end do
        pdf = nm%delta_Fy(xhat) / pdf
    end function f_n_xhat_cond_x


    real(c_double) function f_soft_reverse(n) result(f)
        !! Argument of the finite integral
        real(c_double), intent(in) :: n
        !! Soft metric

        integer :: j, i
        real(c_double), allocatable :: f_n_xhat_cond_x_array(:,:)
        real(c_double) :: log2_f_n_cond_x

        allocate(f_n_xhat_cond_x_array(0:nm%M-1, 0:nm%M-1))
        ! 1st index : X
        ! 2nd index : Xhat

        do j = 0, nm%M-1
            f_n_xhat_cond_x_array(j, :) = f_n_xhat_cond_x(n, [(i, i=0, nm%M-1)], j)
        end do

        f = 0

        do j = 0, nm%M-1 ! X index
            log2_f_n_cond_x = log2(sum(f_n_xhat_cond_x_array(j, :)))
            do i = 0, nm%M-1 ! Xhat index
                f = f + nm%probabilities(j) * f_n_xhat_cond_x_array(j, i) * &
                    (log2(f_n_xhat_cond_x_array(j, i)) - log2_f_n_cond_x)
            end do
        end do
    end function f_soft_reverse

    real(c_double) function H_Xhat(nm) result(H)
        !! Entropy of the output symbols
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper

        H = sum(nm%delta_Fy * log2(nm%delta_Fy))
    end function H_Xhat


    real(c_double) function I_soft_reverse_equidistant_th(snrdb) result(I)
        !! Mutual information of the soft reverse reconciliation scheme
        real(c_double), intent(in) :: snrdb
        !! SNR [dB] at which to calculate the mutual information

        real(c_double) :: Abserr
        integer :: Neval, Ier, Limit, Lenw, Last

        integer :: Iwork(100)
        real(c_double) :: Work(400)
        Limit = 100
        Lenw = 400

        call noisemapper_update_N0_from_snrdb(nm, snrdb)
        call noisemapper_set_y_thresholds(nm)
        call noisemapper_set_Fy_grids(nm)

        call dqags(f_soft_reverse, 0d0, 1d0, 1d-12, 1d-6, &
            I, Abserr, Neval, Ier, &
            Limit, Lenw, Last, Iwork, Work)

        if (Ier /= 1) then
            print '("Error at ", f10.3, " [dB]: error ", i1)', snrdb, Ier
        end if

        I = I + H_Xhat(nm)
    end function I_soft_reverse_equidistant_th


    real(c_double) function I_soft_reverse_uniform_output_th(snrdb) result(I)
        !! Mutual information of the soft reverse reconciliation scheme
        real(c_double), intent(in) :: snrdb
        !! SNR [dB] at which to calculate the mutual information

        real(c_double) :: Abserr
        integer :: Neval, Ier, Limit, Lenw, Last

        integer :: Iwork(100)
        real(c_double) :: Work(400)
        Limit = 100
        Lenw = 400

        call noisemapper_update_N0_from_snrdb(nm, snrdb)
        call noisemapper_set_y_thresholds_uniform(nm)
        call noisemapper_set_Fy_grids(nm)

        call dqags(f_soft_reverse, 0d0, 1d0, 1d-12, 1d-6, &
            I, Abserr, Neval, Ier, &
            Limit, Lenw, Last, Iwork, Work)

        if (Ier /= 1) then
            print '("Error at ", f10.3, " [dB]: error ", i1)', snrdb, Ier
        end if

        I = I + H_Xhat(nm)
    end function I_soft_reverse_uniform_output_th

    ! +-----------------------------+
    ! | Hard reverse reconciliation |
    ! +-----------------------------+

    real(c_double) function H_Xhat_cond_X(nm) result(H)
        !! conditional entropy \(H(\hat{X}|X)\)
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper

        integer :: i, j
        H = 0

        do i = 0, nm%M-1
            do j = 0, nm%M-1
                H = nm%probabilities(j) * nm%fwd_probabilities(j, i) * &
                    log2(nm%fwd_probabilities(j, i))
            end do
        end do
        H = - H
    end function H_Xhat_cond_X


    real(c_double) function I_hard_reverse_equidistant_th(snrdb) result (I)
        !! Mutual information of the discrete Input and Output channel
        real(c_double), intent(in) :: snrdb
        !! SNR [dB] at which to evaluate the Mutual information

        call noisemapper_update_N0_from_snrdb(nm, snrdb)
        call noisemapper_set_y_thresholds(nm)
        call noisemapper_update_hard_reverse_tables(nm)

        I = H_Xhat(nm) - H_Xhat_cond_X(nm)
    end function I_hard_reverse_equidistant_th


    real(c_double) function I_hard_reverse_uniform_output_th(snrdb) result (I)
        !! Mutual information of the discrete Input and Output channel
        real(c_double), intent(in) :: snrdb
        !! SNR [dB] at which to evaluate the Mutual information

        call noisemapper_update_N0_from_snrdb(nm, snrdb)
        call noisemapper_set_y_thresholds_uniform(nm)
        call noisemapper_update_hard_reverse_tables(nm)

        I = H_Xhat(nm) - H_Xhat_cond_X(nm)
    end function I_hard_reverse_uniform_output_th


    ! +----------------------------+
    ! ! Soft direct reconciliation |
    ! +----------------------------+

    real(c_double) function f_integrand_GH(x) result(f)
        !! Function to be integrated with the Gauss-Hermite quadrature rule
        real(c_double), intent(in) :: x
        !! Scaled Y

        integer :: j, k
        real(c_double) :: log_arg, a_j, a_k
        f = 0

        do j = 0, nm%M-1
            log_arg = 0
            a_j = nm%constellation(j)
            do k = 0, nm%M-1
                a_k = nm%constellation(k)
                log_arg = log_arg + nm%probabilities(k) * &
                    exp((a_k-a_j)*(a_j-a_k + 2*x*sq2*nm%sigma)/nm%N0)
            end do
            f = f + nm%probabilities(j) * log2(log_arg)
        end do
        f = f / sqrt(acos(-1d0))
    end function f_integrand_GH


    real(c_double) function f_integrand_soft_direct(x) result(f)
        !! Integrand function for soft direct reconciliation
        real(c_double), intent(in) :: x
        !!

        f = - exp(-x**2) * f_integrand_GH(x)
    end function f_integrand_soft_direct


    real(c_double) function I_direct(snrdb) result(I)
        !! Mutual information of the direct reconciliation scheme
        real(c_double), intent(in) :: snrdb
        !! SNR [dB] at which to calculate the mutual information

        real(c_double) :: Abserr
        integer :: Ier
        ! integer :: Neval, Ier, Limit, Lenw, Last

        ! integer :: Iwork(100)
        ! real(c_double) :: Work(400)
        ! Limit = 100
        ! Lenw = 400

        call noisemapper_update_N0_from_snrdb(nm, snrdb)

        ! call dqagi(f_integrand_soft_direct, 0d0, 2, 1d-12, 1d-6, &
        !     I, Abserr, Neval, Ier, &
        !     Limit, Lenw, Last, Iwork, Work)

        I = - hermite(20, f_integrand_GH, Ier)

        if (Ier /= 0) then
            print '("Error at ", f10.3, " [dB]: error ", i1)', snrdb, Ier
        end if
    end function I_direct
end module re2often_mi
