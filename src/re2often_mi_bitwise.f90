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
submodule (re2often_mi) re2often_mi_bitwise
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    !!
    !! Bitwise mutual information

    real(c_double), parameter :: pi = acos(-1d0)
    real(c_double), parameter :: sqrtPi = sqrt(pi)
    real(c_double), parameter :: twoSqrt2 = 2*sqrt(2d0)
contains
    ! +-----------------------------+
    ! | Soft reverse reconciliation |
    ! +-----------------------------+
    module subroutine compute_pdf_N_B1_cond_X(n, a_j, S, tau)
        !! Computes \(f_{N,B_l|X}(n, 1| a_j)\)
        !! for all \(l\in 1, \dots,\log_2M\
        !! and \(f_{N|X}(n| a_j)\)
        real(c_double), intent(in) :: n
        !! soft metric of RRS
        integer(c_int), intent(in) :: a_j
        !! transmitted symbol
        real(c_double), intent(out) :: S
        !! \(f_{N|X}(n| a_j)\)
        real(c_double), intent(out) :: tau(0:nm%bps-1)
        !! \(f_{N,B_l|X}(n, 1| a_j)\) where l is the index of tau

        integer(c_int) :: a_k, l
        real(c_double) :: pdf_tmp

        S = 0
        tau(:) = 0

        do a_k = 0, nm%M-1
            pdf_tmp = f_n_xhat_cond_x(n, a_k, a_j)
            S = S + pdf_tmp
            do l = 0, nm%bps
                if (nm%s_to_b(a_k, l)) then
                    tau(l) = tau(l) + pdf_tmp
                end if
            end do
        end do
    end subroutine compute_pdf_N_B1_cond_X


    module function f_bitwise(n) result(f)
        !! Argument for the finite integral for bit-wise MI.
        !! NOTE that it is not scaled to log2
        real(c_double), intent(in) :: n
        !! soft metric
        real(c_double) :: f

        integer(c_int) :: a_j
        real(c_double) :: tau(0:nm%bps-1)
        real(c_double) :: S

        f = 0

        do a_j = 0, nm%M-1
            call compute_pdf_N_B1_cond_X(n, a_j, S, tau)
            f = f + nm%probabilities(a_j) * &
                (sum(tau * log0(tau) + (S-tau) * log0(S-tau)) - nm%M*S*log0(S))
        end do
    end function f_bitwise


    module function H_Bl_reverse(nm) result(H)
        !! Entropy of the received bits
        !! NOTE that it is not scaled to log2
        type(noisemapper_type), intent(in) :: nm
        !! Noise mapper
        real(c_double) :: H(0:nm%bps-1)

        integer :: a_i, l

        H(:) = 0

        do a_i = 0, nm%M-1
            do l = 0, nm%bps-1
                if (nm%s_to_b(a_i, l)) then
                    H(l) = H(l) + nm%delta_Fy(a_i)
                end if
            end do
        end do

        ! Now H contains the probability that bit l is 1

        H = - H * log0(H) - (1-H) * log0(1-H)
    end function H_Bl_Reverse

    module function I_soft_reverse_bitwise(snrdb, thresholds, uf) result (I)
        real(c_double), intent(in) :: snrdb
        !! SNR [dB] at which to calculate the mutual information
        real(c_double), intent(in), optional :: thresholds(1:nm%M-1)
        !! Decision Thresholds, overrides `uf`
        logical, intent(in), optional :: uf
        !! Flag for uniform output thresholds
        real(c_double) :: I

        real(c_double) :: Abserr
        integer :: Neval, Ier, Limit, Lenw, Last

        integer :: Iwork(100)
        real(c_double) :: Work(400)

        call noisemapper_update_N0_from_snrdb(nm, snrdb)
        if (present(thresholds)) then
            call noisemapper_set_y_thresholds(nm, thresholds)
        elseif (present(uf)) then
            if (uf) then
                call noisemapper_set_y_thresholds_uniform(nm)
            else
                goto 100
            end if
        else
            goto 100
        end if
        goto 101

100     call noisemapper_set_y_thresholds(nm)
101     call noisemapper_set_Fy_grids(nm)

        Limit = 100
        Lenw = 400

        call dqags(f_bitwise, 0d0, 1d0, 1d-12, 1d-6, &
            I, Abserr, Neval, Ier, &
            Limit, Lenw, Last, Iwork, Work)

        if (Ier /= 0) then
            print '("Error at ", f10.3, " [dB]: error ", i1)', snrdb, Ier
        end if

        I = I + sum(H_Bl_reverse(nm))
        I = I/log(2d0)
    end function I_soft_reverse_bitwise


    ! +-----------------------+
    ! | Direct reconciliation |
    ! +-----------------------+
    module subroutine compute_pdf_Y_cond_B1(u, phi, tau)
        !! Computes \(\tilde{\tau}_l^{(1)}(u) :=
        !! \frac{1}{\sqrt{\pi}}\sum_{a_k \in \mathcal{A}_l^{(1)}}
        !! P_X(a_k) e^{\frac{2a_ku\sqrt{N_0} - a_k^2}{N_0}}\)
        real(c_double), intent(in) :: u
        !! Integration variable for the Gauss Hermite method
        real(c_double), intent(out) :: phi
        !! \(\tilde{\tau}_l^{(1)}(u) + \tilde{\tau}_l^{(1)}(u) \)
        !! (it does not depend on `l`)
        real(c_double), intent(out) :: tau(0:nm%bps-1)
        !! \(\tilde{\tau}_l^{(1)}(u)\) for each `l`

        integer :: k, l
        real(c_double) :: tmp, fourSigma

        fourSigma= twoSqrt2*nm%sigma

        phi = 0
        tau(:) = 0

        do k = 0, nm%M-1
            tmp = nm%probabilities(k) * &
                exp(nm%constellation(k)*(fourSigma*u -nm%constellation(k))/nm%N0)
            phi = phi + tmp
            do l = 0, nm%bps-1
                if (nm%s_to_b(k, l)) then
                    tau(l) = tau(l) + tmp
                end if
            end do
        end do
        phi = phi/sqrtPi
        tau = tau/sqrtPi
    end subroutine compute_pdf_Y_cond_B1


    module function f_integrand_GH_bitwise(u) result(f)
        real(c_double), intent(in) :: u
        real(c_double) :: f

        real(c_double) :: phi, tau(0:nm%bps-1)

        call compute_pdf_Y_cond_B1(u, phi, tau)

        f = sum(tau * log0(tau))
        tau = phi - tau
        f = f + sum(tau * log0(tau))
        f = f - nm%bps * phi * log0(phi)
    end function f_integrand_GH_bitwise


    module function H_Bl_direct() result(H)
        real(c_double) :: H(0:nm%bps-1)

        integer :: l, k

        H(:) = 0
        do k = 0, nm%M-1
            do l = 0, nm%bps-1
                if (nm%s_to_b(k, l)) then
                    H(l) = H(l) + nm%probabilities(k)
                end if
            end do
        end do
        H = - H * log0(H) - (1-H) * log0(1-H)
    end function H_Bl_direct


    module function I_direct_bitwise(snrdb) result(I)
        real(c_double), intent(in) :: snrdb
        !! SNR [dB] at which to calculate the mutual information
        real(c_double) :: I

        integer :: Ier

        call noisemapper_update_N0_from_snrdb(nm, snrdb)

        I = hermite(20, f_integrand_GH_bitwise, Ier)

        if (Ier /= 0) then
            print '("Error at ", f10.3, " [dB]: error ", i1)', snrdb, Ier
        end if

        I = (I + sum(H_Bl_direct())) / log(2d0)
    end function I_direct_bitwise
end submodule re2often_mi_bitwise
