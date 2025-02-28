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
program mi_biawgn
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    use iso_fortran_env, only: dp => real64
    use io_fortran_lib, only: to_file
    use external_hermite, only: hermite

    implicit none

    real(dp), parameter :: ln_of_2 = log(2d0)
    real(dp), parameter :: twoSqrtPi = 2*sqrt(acos(-1d0))

    real(dp) :: N_0
    real(dp), allocatable, target :: outdata(:,:)
    real(dp), pointer :: snrdb_array(:)
    real(dp), pointer :: mi_array(:)
    integer  :: i, Npts
    real(dp) :: start, stp

    start = -50
    stp   = 20
    Npts = int(ceiling((stp - start)/.1d0)) + 1

    allocate(outdata(Npts, 2))
    snrdb_array => outdata(:, 1)
    mi_array    => outdata(:, 2)

    snrdb_array = [(start + i*0.1d0, i=0, Npts-1)]
    mi_array(:) = 2 ! impossible value, so to detect errors...

    do i = 1, Npts
        mi_array(i) = mutual_information(snrdb_array(i))
    end do

    call to_file(outdata, file="mi_biawgn.csv", fmt="e")

contains
    real(dp) elemental function log_base_2(arg) result (l)
        real(dp), intent(in) :: arg

        l = log(arg)/ln_of_2
    end function log_base_2


    real(dp) function f_gh(x) result(f)
        real(dp), intent(in) :: x

        f = log_base_2(1 + 2*exp(-4/N_0)*cosh(4*x/sqrt(N_0)) + exp(-8/N_0))
    end function f_gh


    real(dp) function mutual_information(snrdb) result(mi)
        real(dp), intent(in) :: snrdb

        integer :: ier

        N_0 = 10_dp ** (-snrdb/10_dp)

        mi = 1 - hermite(20, f_gh, ier)/twoSqrtPi
    end function mutual_information

end program mi_biawgn
