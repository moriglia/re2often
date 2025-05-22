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
submodule (re2often_mi) re2often_mi_gmi
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Generalized mutual information

    integer(c_int) :: gmi_xtilde
    real(c_double) :: gmi_s
    logical :: gmi_soft_reverse_useDenominator
    procedure(q_soft_direct), pointer :: qfun_sd
    procedure(q_hard_hard_soft), pointer :: qfun_map_reverse_soft
contains
    ! function expected_metric_hard(q, p, s) result (E_of_q)
    !     !! Compute the expectation of the values given in `q`,
    !     !! when every value `q(i)` has probability `p(i)`.
    !     real(c_double), intent(in) :: q(:)
    !     !! Array of which to compute the expectation
    !     real(c_double), intent(in) :: p(size(q))
    !     !! Probability weights
    !     real(c_double), intent(in), optional :: s
    !     !! Exponent, if not present, 1 is assumed
    !     real(c_double) :: E_of_q
    !     !! Expectation of \(q^s\)

    !     if (present(s)) then
    !         E_of_q = sum(p*(q**s))
    !     else
    !         E_of_q = sum(p*q)
    !     end if
    ! end function expected_metric_hard


    ! I don't need this function right now. Better not to compile it
    ! function I_s_map_hard_direct(q, s) result(I_s)
    !     !! GMI for MAP criterion with hard information only, direct direction
    !     procedure(q_hard) :: q
    !     !! Function that computes the metric
    !     real(c_double), intent(in), optional :: s
    !     !! Parameter s of the GMI
    !     real(c_double) :: I_s
    !     !! GMI(s)

    !     integer(c_int) :: x, xhat
    !     real(c_double) :: q_vals(0:nm%M-1)
    !     real(c_double) :: q_expected

    !     I_s = 0

    !     do x = 0, nm%M-1
    !         do xhat = 0, nm%M-1
    !             q_vals(xhat) = q(x, xhat)
    !         end do

    !         if (present(s)) then
    !             q_expected = log0(sum( (q_vals**s) * nm%delta_Fy))
    !             q_vals = s * log0(q_vals)
    !         else
    !             q_expected = log0(sum(q_vals * nm%delta_Fy))
    !             q_vals = log0(q_vals)
    !         end if
    !         q_vals = q_vals - q_expected

    !         I_s = I_s + nm%probabilities(x) * sum(q_vals * nm%fwd_probabilities(x, :))
    !     end do

    !     I_s = I_s / log0(2d0)
    ! end function I_s_map_hard_direct


    ! +------------------+
    ! | MAP hard reverse |
    ! +------------------+

    module function q_map_hard_product(xhat, x) result (q)
        integer(c_int), intent(in) :: xhat
        integer(c_int), intent(in) :: x
        real(c_double) :: q

        real(c_double) :: Pb_given_x
        integer(c_int) :: k, i

        q = 1d0

        do k = 0, nm%bps-1
            Pb_given_x = 0
            do i = 0, nm%M-1
                if (nm%s_to_b(xhat, k) .eqv. nm%s_to_b(i, k)) then
                    Pb_given_x = Pb_given_x + nm%fwd_probabilities(x, i)
                end if
            end do
            q = q * Pb_given_x
        end do
    end function q_map_hard_product


    module function q_map_hard_opt(xhat, x) result(q)
        integer(c_int), intent(in) :: xhat
        integer(c_int), intent(in) :: x
        real(c_double) :: q

        q = nm%fwd_probabilities(x, xhat)
    end function q_map_hard_opt


    module function I_s_map_hard_reverse(q, s) result(I_s)
        !! GMI for MAP criterion with hard information only, reverse direction
        procedure(q_hard) :: q
        !! Function that computes the metric
        real(c_double), intent(in), optional :: s
        !! Positive parameter s of the GMI
        real(c_double) :: I_s
        !! GMI(s)
        integer(c_int) :: x, xhat
        real(c_double) :: q_vals(0:nm%M-1)
        real(c_double) :: q_expected

        I_s = 0

        do xhat = 0, nm%M-1
            do x = 0, nm%M-1
                q_vals(x) = q(xhat, x)
            end do

            if (present(s)) then
                q_expected = log0( sum( (q_vals**s) * nm%probabilities) )
                q_vals = s * log0( q_vals )
            else
                q_expected = log0( sum(q_vals * nm%probabilities) )
                q_vals = log0( q_vals )
            end if
            q_vals = q_vals - q_expected
            I_s = I_s + sum( nm%probabilities * nm%fwd_probabilities(:, xhat) * q_vals )
        end do

        I_s = I_s / log0(2d0)
    end function I_s_map_hard_reverse

    ! +-----------------+
    ! | MAP soft direct |
    ! +-----------------+

    function f_GH_map_expectation_qs(u) result(f)
        real(c_double), intent(in) :: u
        real(c_double) :: f

        integer(c_int) :: x

        f = 0

        do x = 0, nm%M-1
            f = f + nm%probabilities(x) &
                * qfun_sd(gmi_xtilde, nm%constellation(x) + sqrtN0*u)**gmi_s
        end do
    end function f_GH_map_expectation_qs


    function f_GH_map_expectation_ln_qs(u) result(f)
        real(c_double), intent(in) :: u
        real(c_double) :: f

        integer(c_int) :: x

        f = 0

        do x = 0, nm%M-1
            f = f + nm%probabilities(x) * &
                log0( qfun_sd(x, nm%constellation(x) + sqrtN0 * u) )
        end do
    end function f_GH_map_expectation_ln_qs


    module function I_s_map_soft_direct(q, s) result(I_s)
        !! GMI for MAP criterion with hard information only, reverse direction
        procedure(q_soft_direct) :: q
        !! Function that computes the metric
        real(c_double), intent(in), optional :: s
        !! Positive parameter s of the GMI
        real(c_double) :: I_s
        !! GMI(s)
        integer(c_int) :: x
        integer :: Ier

        qfun_sd => q

        I_s = hermite(20, f_GH_map_expectation_ln_qs, Ier)

        if (present(s)) then
            I_s = I_s * s
            gmi_s = s
        else
            gmi_s = 1
        end if

        I_s = I_s / sqrtPi

        do x = 0, nm%M-1
            gmi_xtilde = x
            I_s = I_s - nm%probabilities(x) * log0(hermite(20, f_GH_map_expectation_qs, Ier))
        end do
        I_s = I_s + log(sqrtPi)
        I_s = I_s / log(2d0)
    end function I_s_map_soft_direct


    module function q_map_soft_direct_prod(x, y) result (q)
        integer(c_int), intent(in) :: x
        real(c_double), intent(in) :: y
        real(c_double) :: q

        real(c_double) :: f_Y_joint_X(0:nm%M-1)
        real(c_double) :: P_Bl_joint_Y
        integer :: l, k

        f_Y_joint_X = nm%probabilities * pdf_normal(y, loc=nm%constellation, scale=nm%sigma)

        q = 1d0

        do l = 0, nm%bps-1
            P_Bl_joint_Y = 0
            do k = 0, nm%M-1
                if (nm%s_to_b(k,l) .eqv. nm%s_to_b(x, l)) then
                    P_Bl_joint_Y = P_Bl_joint_Y + f_Y_joint_X(k)
                end if
            end do
            q = q * P_Bl_joint_Y
        end do
        q = q / (sum(f_Y_joint_X)**nm%bps)
    end function q_map_soft_direct_prod


    ! +----------------+
    ! | ML hard direct |
    ! +----------------+
    module function I_s_ml_hard_direct(q, s) result (I_s)
        procedure(q_hard) :: q
        real(c_double), intent(in), optional :: s
        real(c_double) :: I_s

        real(c_double) :: q_buff(0:nm%M-1)
        real(c_double) :: tmp
        integer :: xhat, x

        real(c_double) :: s_local

        if (present(s)) then
            s_local = s
        else
            s_local = 1d0
        end if

        I_s = 0
        do xhat = 0, nm%M-1
            do x = 0, nm%M-1
                q_buff(x) = q(x, xhat)
            end do

            tmp = 0
            do x = 0, nm%M-1
                tmp = tmp + nm%fwd_probabilities(x, xhat)*nm%probabilities(x) &
                    * log0(q_buff(x))
            end do
            tmp = s_local * tmp &
                - nm%delta_Fy(xhat) * log0(sum(nm%probabilities * q_buff**s_local))
            I_s = I_s + tmp
        end do

        I_s = I_s / log(2d0)
    end function I_s_ml_hard_direct


    module function q_ml_hard_direct_prod(x, xhat) result(q)
        integer(c_int), intent(in) :: x
        integer(c_int), intent(in) :: xhat
        real(c_double) :: q

        integer :: x_local, l

        real(c_double) :: frac_n, frac_d

        q = 1d0

        do l = 0, nm%bps - 1
            frac_n = 0
            frac_d = 0
            do x_local = 0, nm%M-1
                if (nm%s_to_b(x_local, l) .eqv. nm%s_to_b(x, l)) then
                    frac_n = frac_n &
                        + nm%probabilities(x_local) &
                        * nm%fwd_probabilities(x_local, xhat)
                    frac_d = frac_d + nm%probabilities(x_local)
                end if
            end do
            q = q * frac_n / frac_d
        end do
    end function q_ml_hard_direct_prod

    ! +------------------+
    ! | MAP soft reverse |
    ! +------------------+
    function f_N_GH_map_expectation_qs(n) result (f)
        real(c_double), intent(in) :: n
        real(c_double)             :: f

        integer :: x, xhat
        real(c_double) :: tmp

        f = 0
        do x = 0, nm%M-1
            tmp = 0
            do xhat = 0, nm%M-1
                tmp = tmp + f_n_xhat_cond_x(n, xhat, x)
            end do
            f = f + tmp * nm%probabilities(x) * qfun_map_reverse_soft(x, gmi_xtilde, n)**gmi_s
        end do
    end function f_N_GH_map_expectation_qs


    function f_N_GH_map_expectation_log_qs(n) result(f)
        real(c_double), intent(in) :: n
        real(c_double)             :: f

        integer        :: x, xhat
        real(c_double) :: tmp

        f = 0
        do x = 0, nm%M-1
            tmp = 0
            do xhat = 0, nm%M-1
                tmp = tmp + f_n_xhat_cond_x(n, xhat, x) &
                    * log0(qfun_map_reverse_soft(x, xhat, n)**gmi_s)
            end do
            f = f + nm%probabilities(x) * tmp
        end do
    end function f_N_GH_map_expectation_log_qs


    module function I_s_map_soft_reverse(q, s, useDenominator) result(I_s)
        !! GMI-MAP
        procedure(q_hard_hard_soft)          :: q
        real(c_double), intent(in), optional :: s
        logical       , intent(in), optional :: useDenominator
        real(c_double)                       :: I_s

        ! Data for DQAGS
        real(c_double) :: Abserr
        integer :: Neval, Ier, Limit, Lenw, Last
        integer :: Iwork(100)
        real(c_double) :: Work(400)

        ! Further auxiliaries
        real(c_double) :: I_aux
        integer :: xhat

        Limit = 100
        Lenw = 400

        if (present(s)) then
            gmi_s = s
        else
            gmi_s = 1d0
        end if

        if (present(useDenominator)) then
            gmi_soft_reverse_useDenominator = useDenominator
        else
            gmi_soft_reverse_useDenominator = .true.
        end if

        qfun_map_reverse_soft => q

        call dqags(f_N_GH_map_expectation_log_qs, 0d0, 1d0, 1d-12, 1d-6, &
            I_s, Abserr, Neval, Ier, &
            Limit, Lenw, Last, Iwork, Work)

        if (Ier /= 0) then
            print '("DQAGS error in f_N_GH_map_expectation_log_qs ", i1)', Ier
        end if

        I_s = I_s * gmi_s

        do xhat = 0, nm%M-1
            gmi_xtilde = xhat
            call dqags(f_N_GH_map_expectation_qs, 0d0, 1d0, 1d-12, 1d-6, &
                I_aux, Abserr, Neval, Ier, &
                Limit, Lenw, Last, Iwork, Work)
            if (Ier /= 0) then
                print '("DQAGS error in f_N_GH_map_expectation_log_qs ", i1)', Ier
            end if
            I_s  = I_s - nm%probabilities(xhat) * log0(I_aux)
        end do
        I_s = I_s / log(2d0)
    end function I_s_map_soft_reverse


    module function q_map_soft_reverse_prod(x, xhat, n) result(q)
        integer(c_int), intent(in) :: x
        integer(c_int), intent(in) :: xhat
        real(c_double), intent(in) :: n
        real(c_double)             :: q

        real(c_double) :: f_N_given_X, tmp
        integer :: l, xhat_var

        q = 1d0

        f_N_given_X = 0

        if (gmi_soft_reverse_useDenominator) then
            do xhat_var = 0, nm%M-1
                f_N_given_X = f_N_given_X + f_n_xhat_cond_x(n, xhat_var, x)
            end do
        end if

        do l = 0, nm%bps-1
            tmp = 0
            do xhat_var = 0, nm%M-1
                if (nm%s_to_b(xhat, l) .eqv. nm%s_to_b(xhat_var, l)) then
                    tmp = tmp + f_n_xhat_cond_x(n, xhat_var, x)
                end if
            end do
            q = q * tmp
        end do
        if (gmi_soft_reverse_useDenominator) then
            q = q / (f_N_given_X**nm%bps)
        end if
    end function q_map_soft_reverse_prod
end submodule re2often_mi_gmi
