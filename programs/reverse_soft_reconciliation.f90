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
program reverse_soft_reconciliation
  use sim_rrs, only: simulate_reverse_pam
  use io_fortran_lib, only: from_file, to_file
  use iso_c_binding, only: wp => c_double, ki => c_int, kl => c_bool
  implicit none

  integer :: argc
  character(len=500), allocatable :: argv(:)

  character(500) :: tanner_file
  character(500) :: output_file
  real(wp)    :: snr(2)         ! Signal to Noise Ratio in dB
  integer(ki) :: nsnr           ! Number of SNR points
  integer(ki) :: bps            ! Bits per symbol
  integer(ki) :: min_ferr       ! Minimum frame errors
  integer(ki) :: max_sim        ! Maximum simulation loops
  integer(ki) :: min_sim        ! Minimum simulation loops
  integer(ki) :: max_iter       ! Maximum number of LDPC iterations

  integer(ki), allocatable :: edge_definition(:,:)
  real(wp), allocatable :: outdata(:,:)
  
  integer :: i

  argc = command_argument_count()
  allocate(argv(argc))

  do i = 1, argc
     call get_command_argument(i, argv(i))
  end do
  
  i = 1
  do while(i <= argc)
     if (argv(i) == "--nsnr") then
        read(argv(i+1),*) nsnr
        i = i + 2
        ! print *, "nsnr", nsnr
     elseif (argv(i) == "--snr") then
        read(argv(i+1), *) snr(1)
        read(argv(i+2), *) snr(2)
        i = i+3
        ! print *, "snr", snr
     elseif (argv(i) == "--bps") then
        read(argv(i+1),*) bps
        i = i + 2
        ! print *, "bps", bps
     elseif (argv(i) == "--minframeerr") then
        read(argv(i+1), *) min_ferr
        i = i + 2
        ! print *, "ferr", min_ferr
     elseif (argv(i) == "--minsimloops") then
        read(argv(i+1), *) min_sim
        i = i + 2
        ! print *, "min_sim", min_sim
     elseif (argv(i) == "--maxsimloops") then
        read(argv(i+1), *) max_sim
        i = i + 2
        ! print *, "max_sim", max_sim
     elseif (argv(i) == "--maxldpciter") then
        read(argv(i+1), *) max_iter
        i = i + 2
        ! print *, "ldpc", max_iter
     elseif (argv(i) == "--tannerfile") then
        call get_command_argument(i+1, tanner_file)
        i = i + 2
        ! print*, tanner_file
     elseif (argv(i) == "--csvout") then
        call get_command_argument(i+1, output_file)
        i = i + 2
        ! print *, output_file
     else
        print *, "Unrecognized argument: ", argv(i)
        stop
     end if
  end do

  call from_file(file=tanner_file, into=edge_definition, header=.true.)

  allocate(outdata(nsnr,3))
  outdata(:, 1) = [(snr(1) + real(i, wp)*(snr(2)-snr(1))/real(nsnr - 1, wp), i = 0, nsnr-1)]
  outdata(:,2:) = 0.0_wp
  
  call simulate_reverse_pam(outdata(:,1), nsnr, bps, &
       min_ferr, max_sim, min_sim, &
       edge_definition(2:,2), & ! e_to_v is vid in the tanner file
       edge_definition(2:,3), & ! e_to_c is cid in the tanner file
       edge_definition(1,1), &  ! Number of edges is the first entry of EID in the tanner file
       max_iter, &
       outdata(:,2), outdata(:,3))

  if (this_image() == 1) then
     call to_file(x=outdata, file=output_file, header=["SNR", "BER", "FER"], fmt="f")
  end if

end program reverse_soft_reconciliation
