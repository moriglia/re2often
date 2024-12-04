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
  use io_fortran_lib, only: from_file, to_file
  use iso_c_binding, only: wp => c_double, ki => c_int, kl => c_bool
  use ldpc_decoder, only: TDecoder
  use noise_mapper, only: TNoiseMapper
  use forbear, only: bar_object
  use stdlib_random, only: stdlib_random_seed => random_seed
  implicit none

  integer :: argc
  character(len=500), allocatable :: argv(:)

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
  real(wp), pointer :: outdata(:,:)

  real(wp), pointer :: snrdb_array(:), ber(:), fer(:)
  
  integer :: i
  
  type(TDecoder) :: decoder
  type(TNoiseMapper) :: noisemapper
  integer :: time_seed
  integer:: seed(8)

  integer :: i_snr, i_frame
  integer :: N_symb, K

  integer, allocatable  :: x_i(:)
  real(wp), allocatable :: x(:)
  real(wp), allocatable :: y(:)
  logical, allocatable :: word(:)
  logical, allocatable :: synd(:)
  real(wp), allocatable :: lappr(:)
  real(wp), allocatable :: lappr_updated(:)

  integer :: new_errors, n_ldpc_it

  integer, allocatable :: b_error_count(:)[:]
  integer, allocatable :: f_error_count(:)[:]
  integer, allocatable :: f_count(:)[:]

  type(bar_object) :: progress_bar
  

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
     elseif (argv(i) == "--snr") then
        read(argv(i+1), *) snr(1)
        read(argv(i+2), *) snr(2)
        i = i+3
     elseif (argv(i) == "--bps") then
        read(argv(i+1),*) bps
        i = i + 2
     elseif (argv(i) == "--minframeerr") then
        read(argv(i+1), *) min_ferr
        i = i + 2
     elseif (argv(i) == "--minsimloops") then
        read(argv(i+1), *) min_sim
        i = i + 2
     elseif (argv(i) == "--maxsimloops") then
        read(argv(i+1), *) max_sim
        i = i + 2
     elseif (argv(i) == "--maxldpciter") then
        read(argv(i+1), *) max_iter
        i = i + 2
     elseif (argv(i) == "--tannerfile") then
        call get_command_argument(i+1, tanner_file)
        i = i + 2
     elseif (argv(i) == "--csvout") then
        call get_command_argument(i+1, output_file)
        i = i + 2
     else
        print *, "Unrecognized argument: ", argv(i)
        stop
     end if
  end do  

  
  call from_file(file=tanner_file, into=edge_definition, header=.true.)

  allocate(outdata(nsnr,3))

  snrdb_array => outdata(:,1)
  ber         => outdata(:,2)
  fer         => outdata(:,3)
  
  snrdb_array = [(snr(1) + real(i, wp)*(snr(2)-snr(1))/real(nsnr - 1, wp), i = 0, nsnr-1)]
  outdata(:,2:) = 0.0_wp
  
  decoder = TDecoder(&
       edge_definition(1 , 1), &
       edge_definition(2:, 2), &
       edge_definition(2:, 3))
  noisemapper = TNoiseMapper(B=bps, sigma=1.0_wp)
  
  K = decoder%vnum - decoder%cnum
  N_symb = decoder%vnum / bps
  allocate(x(N_symb))
  allocate(x_i(N_symb))
  allocate(y(N_symb))
  allocate(word(decoder%vnum))
  allocate(synd(decoder%cnum))
  allocate(lappr(decoder%vnum))
  allocate(lappr_updated(decoder%vnum))

  call random_seed(get=seed)
  call sleep(this_image())
  call system_clock(time_seed)
  seed(1) = mod(seed(1)*sum(seed(this_image():)), abs(time_seed) + 2) - 12399027
  seed(2) = 467738 * this_image() * (seed(3) - time_seed) + iand(sum(seed(:this_image())), time_seed)
  seed(3) = (time_seed + this_image()) * (2 + mod(abs(883 + seed(mod(time_seed, 8) + 1)), this_image())) + &
       mod(sum(seed), max(maxval(seed(:)), 37))
  seed(4) = mod(987654321*time_seed, abs(time_seed - this_image()*this_image()*this_image()) + 3)
  seed(5) = seed(7)/7 - time_seed*this_image() + this_image() ** mod(abs(seed(5)), 29)
  seed(6) = sum(seed(::2) ** abs(time_seed)) + 929812093 / sum(seed(1::2))
  seed(7) = this_image() * seed(6) * 661242 - this_image() / (seed(6) + time_seed**abs(sum(seed)))
  seed(8) = time_seed - seed(4)**this_image() + sum(seed(2::2)) ** abs(time_seed) + 329999999/time_seed
  call stdlib_random_seed(sum(seed), time_seed) ! Necessary for the random functions in the stdlib
  
  allocate(b_error_count(nsnr)[*])
  allocate(f_error_count(nsnr)[*])
  allocate(f_count(nsnr)[*])
  
  b_error_count(:) = 0
  f_error_count(:) = 0
  f_count(:) = 0
  
  sync all
  
  if (this_image() == 1) then
     call progress_bar%initialize(&
          filled_char_string='+', prefix_string='SNR points progress |',&
          suffix_string='| ', add_progress_percent=.true.)
     call progress_bar%start
  end if
  
  snr_loop : do i_snr = 1 , nsnr
     call noisemapper%update_noise_sigma_based_on_snr(snrdb_array(i_snr))
     
     frame_loop : do i_frame = 1, max_sim
        call noisemapper%random_symbol(x_i)
        word = noisemapper%symbol_to_grey_array(x_i)
        synd = decoder%word_to_synd(word)
        
        x = noisemapper%symbol_index_to_real(x_i)
        call noisemapper%add_noise(x, y)
        lappr = noisemapper%y_to_lappr_norm(y)

        n_ldpc_it = max_iter
        call decoder%decode(lappr, lappr_updated, synd, n_ldpc_it)
        new_errors = word_lappr_errors(word(:K), lappr_updated(:K))
        
        critical
          if (new_errors > 0) then
             b_error_count(i_snr)[1] = b_error_count(i_snr)[1] + new_errors
             f_error_count(i_snr)[1] = f_error_count(i_snr)[1] + 1
          end if
          f_count(i_snr)[1] = f_count(i_snr)[1] + 1
        end critical

        if ((f_error_count(i_snr)[1] > min_ferr) .and. &
             (f_count(i_snr)[1] > min_sim)) then
           exit frame_loop
        end if
     end do frame_loop

     if (this_image()==1) then
        call progress_bar%update(current=real(i_snr, 8)/real(Nsnr, 8))
     end if
  end do snr_loop

  sync all

  if (this_image() == 1) then
     ber = real(b_error_count, wp)/real(f_count*K, wp)
     fer = real(f_error_count, wp)/real(f_count, wp)

     call to_file(x=outdata, file=output_file, header=["SNR", "BER", "FER"], fmt="f")
  end if
  
  deallocate(outdata)

contains

  function word_lappr_errors(word, lappr) result (errors)
    logical,  intent(in) :: word(:)
    real(wp), intent(in) :: lappr(size(word))

    integer :: errors

    integer :: i
    errors = 0
    do i = 1, size(word)
       if (word(i) .neqv. (lappr(i) < 0)) then
          errors = errors + 1
       end if
    end do
  end function word_lappr_errors

end program direct_reconciliation
