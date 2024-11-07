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
program direct_reconciliation
  use sim_direct_channel, only: simulate_pam
  use flap, only: command_line_interface
  use io_fortran_lib, only: from_file, to_file
  use iso_c_binding, only: wp => c_double, ki => c_int, kl => c_bool
  implicit none

  type(command_line_interface) :: cli

  character(1000) :: tanner_file
  character(1000) :: output_file
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


  call cli%init(progname="direct_reconciliation", &
       description="Evaluate bit error rate for the direct reconciliation with a PAM modulation", &
       license="GPL-3.0-or-later")

  call cli%add(switch="--snr", nargs='2', help="Signal to Noise Ratio range [dB]", required=.true.)
  call cli%add(switch="--nsnr", help="Number of SNR points", required=.true.)
  call cli%add(switch="--bps", help="Bits per Symbol", required=.true.)
  call cli%add(switch="--minframeerr", help="Minimum frame errors", required=.true.)
  call cli%add(switch="--maxsimloops", help="Maximum simulation loops", required=.true.)
  call cli%add(switch="--minsimloops", help="Minimum simulation loops", required=.true.)
  call cli%add(switch="--tannerfile", help="CSV file with Tanner graph definition", required=.true.)
  call cli%add(switch="--maxldpciter", help="Maximum number of LDPC iterations", required=.true.)
  call cli%add(switch="--csvout", help="Output file", required=.true.)

  call cli%get(val=snr, switch="--snr")
  call cli%get(val=nsnr, switch="--nsnr")
  call cli%get(val=bps, switch="--bps")
  call cli%get(val=min_ferr, switch="--minframeerr")
  call cli%get(val=max_sim, switch="--maxsimloops")
  call cli%get(val=min_sim, switch="--minsimloops")
  call cli%get(val=max_iter, switch="--maxldpciter")
  call cli%get(val=tanner_file, switch="--tannerfile")
  call cli%get(val=output_file, switch="--csvout")
  

  
  call from_file(file=tanner_file, into=edge_definition, header=.true.)

  allocate(outdata(nsnr,3))
  outdata(:, 1) = [(snr(1) + real(i, wp)*(snr(2)-snr(1))/real(nsnr - 1, wp), i = 0, nsnr-1)]
  outdata(:,2:) = 0.0_wp

  


  call simulate_pam(outdata(:,1), nsnr, bps, &
       min_ferr, max_sim, min_sim, &
       edge_definition(2:,2), & ! e_to_v is vid in the tanner file
       edge_definition(2:,3), & ! e_to_c is cid in the tanner file
       edge_definition(1,1), &  ! Number of edges is the first entry of EID in the tanner file
       max_iter, &
       outdata(:,2), outdata(:,3))

  if (this_image() == 1) then
     call to_file(x=outdata, file=output_file, header=["SNR", "BER", "FER"], fmt="f")
  end if

end program direct_reconciliation
