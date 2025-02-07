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
module re2often_utils
    !! author: Marco Origlia
    !! license: GPL-3.0-or-later
    use io_fortran_lib, only: to_file
    implicit none

    public :: binsearch, save_data, make_directory_and_file_name, logical2integer

contains
    pure function binsearch(vector, y) result (index)
        !! Retrieve index of datum within a vector
        double precision, intent(in) :: vector(:)
        !! 1-based array \(\mathbf{v}\) to search data from
        double precision, intent(in) :: y
        !! Datum \(y\) to search
        integer :: index
        !! `index` \(i\) such that \(\mathbf{v}_i \leq y \lt \mathbf{v}_{i+1}\).
        !!  `index` is set to 0 if \(y \lt \mathbf{v}_{1}\).
        !!  `index` is set to `size(vector)` if \(y \geq \mathbf{v}_{\texttt{size(vector)}}\)

        integer :: index_u
        integer :: index_l

        if (y .lt. vector(1)) then
            index = 0
            return
        end if
        if (y .ge. vector(size(vector))) then
            index = size(vector)
            return
        end if

        index_u   = size(vector)
        index_l   = 1
        index = ishft(index_u + index_l, -1)
        do while ( (index_u - index_l) .gt. 1)
            if (y .ge. vector(index)) then
                index_l   = index
            else
                index_u   = index
            end if
            index = ishft(index_u + index_l, -1)
        end do
    end function binsearch


    subroutine save_data(data, root_dir, bps, isReverse, isHard, snr, nsnr, min_sim, max_sim, max_iter, min_ferr, alpha)
        !! Save data in the appropriate directory with a name consistend with simulation parameters
        double precision, intent(in) :: data(:,:)
        !! data to be saved in 3 columns: SNR [dB], BER, FER
        character(*), intent(in) :: root_dir
        !! Root directory of results without trailing "/"
        integer, intent(in) :: bps
        !! Bit per symbol
        logical, intent(in) :: isReverse
        !! Indicates whether the results refer to a reverse reconciliation problem
        logical, intent(in) :: isHard
        !! Indicates whether the results refer to a hard reverse reconciliation
        !! (ignored if `isReverse` is `.false.`)
        double precision, intent(in) :: snr(2)
        !! Start and stop SNR values
        integer, intent(in) :: nsnr
        !! Number of snr points
        integer, intent(in) :: min_sim
        !! minimum number of channel realizations
        integer, intent(in) :: max_sim
        !! maximum number of channel relizations
        integer, intent(in) :: max_iter
        !! maximum number of LDPC iterations
        integer, intent(in) :: min_ferr
        !! Number of frame errors to stop the simulation earlier
        double precision, optional, intent(in) :: alpha
        !! Scaling factor for the LAPPR

        character(len=250) :: dir
        character(len=250) :: fn
        character(len=500) :: file_name

        if (present(alpha)) then
            call make_directory_and_file_name(root_dir, bps, isReverse, isHard, &
                snr, nsnr, min_sim, max_sim, max_iter, min_ferr,                &
                dir, fn, alpha)
        else
            call make_directory_and_file_name(root_dir, bps, isReverse, isHard, &
                snr, nsnr, min_sim, max_sim, max_iter, min_ferr,                &
                dir, fn)
        end if

        call execute_command_line("mkdir -p " // trim(dir))

        write(file_name, '(A, "/", A, ".csv")') trim(dir), trim(fn)
        call to_file(x=data(:,:3), file=trim(file_name), &
            header=["SNR", "BER", "FER"], fmt="f")
    end subroutine save_data


    subroutine make_directory_and_file_name(root_dir, bps, isReverse, isHard, &
        snr, nsnr, min_sim, max_sim, max_iter, min_ferr, &
        output_dir, output_file, alpha)
        character(*), intent(in) :: root_dir
        !! Root directory of results without trailing "/"
        integer, intent(in) :: bps
        !! Bit per symbol
        logical, intent(in) :: isReverse
        !! Indicates whether the results refer to a reverse reconciliation problem
        logical, intent(in) :: isHard
        !! Indicates whether the results refer to a hard reverse reconciliation
        !! (ignored if `isReverse` is `.false.`)
        double precision, intent(in) :: snr(2)
        !! Start and stop SNR values
        integer, intent(in) :: nsnr
        !! Number of snr points
        integer, intent(in) :: min_sim
        !! minimum number of channel realizations
        integer, intent(in) :: max_sim
        !! maximum number of channel relizations
        integer, intent(in) :: max_iter
        !! maximum number of LDPC iterations
        integer, intent(in) :: min_ferr
        !! Number of frame errors to stop the simulation earlier
        character(len=250), intent(out) :: output_dir
        !! Output dir
        character(len=250), intent(out) :: output_file
        !! Output file
        double precision, optional, intent(in) :: alpha
        !! Scaling coefficient for the LAPPR

        output_dir = trim(root_dir) // "/"

        if (isReverse) then
            if (isHard) then
                output_dir = trim(output_dir) // "rev-h" ! for reverse hard
            else
                output_dir = trim(output_dir) // "rev-s" ! for reverse soft
            end if
        else
            output_dir = trim(output_dir) // "dir" ! for direct
        end if
        write(output_dir, '(A, "/bps", I0)') trim(output_dir), bps


        write(output_file, '("snr_", 2(SP, F0.3, "_"))') snr(1), snr(2)
        write(output_file, '(A, I0)') trim(output_file), nsnr
        write(output_file, '(A, "_it", I0)') trim(output_file), max_iter
        write(output_file, '(A, "_sim", I4.4, "_", I6.6, "_ferr", I4.4)') trim(output_file), &
            min_sim, max_sim, min_ferr

        if (present(alpha)) then
            output_dir =  trim(output_dir) // "/alpha"
            write(output_file, '(A, f8.6)') trim(output_file) // "_a", alpha
        end if
    end subroutine make_directory_and_file_name

    pure function logical2integer(arr) result(res)
        logical, intent(in) :: arr(0:)
        integer :: res

        integer :: i

        res = 0
        do i = 0, size(arr) - 1
            if (arr(i)) then
                res = res + ishft(1, i)
            end if
        end do
    end function logical2integer
end module re2often_utils
