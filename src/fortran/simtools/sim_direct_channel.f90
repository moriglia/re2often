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
module sim_direct_channel
  use iso_c_binding, only: wp => c_double, ki => c_int, kl => c_bool
  use ldpc_decoder, only: TDecoder
  use alpha_pam, only: TAlphaPAM
  use stdlib_stats_distribution_normal, only: rvs_normal
  implicit none

contains
  
  subroutine simulate_pam (snrdb_array, bps, &
       min_ferr, max_loops, min_loops, &
       e_to_v, e_to_c, ldpc_iterations, &
       ber, fer) bind (C)
    real(wp), intent(in) :: snrdb_array(:)
    integer(ki), intent(in)  :: bps, min_ferr, max_loops, min_loops
    integer(ki), intent(in)  :: e_to_v(:)
    integer(ki), intent(in)  :: e_to_c(size(e_to_v))
    integer(ki), intent(in)  :: ldpc_iterations
    real(wp), intent(out)    :: ber(size(snrdb_array))
    real(wp), intent(out)    :: fer(size(snrdb_array))

    type(TDecoder) :: decoder
    type(TAlphaPAM) :: alphabet
    integer :: time_seed, seed(8)

    integer :: i_snr, i_frame
    integer :: N_symb, K

    integer, allocatable :: x(:)
    real(wp), allocatable :: y(:)
    logical, allocatable :: word(:)
    logical, allocatable :: synd(:)
    real(wp), allocatable :: llr(:)
    real(wp), allocatable :: llr_updated(:)

    real(wp) :: N0, N0_half

    integer :: new_errors, n_ldpc_it

    integer, allocatable :: b_error_count(:)[:]
    integer, allocatable :: f_error_count(:)[:]
    integer, allocatable :: f_count(:)[:]
    
    decoder = TDecoder(size(e_to_v), e_to_v, e_to_c)
    alphabet = TAlphaPAM(bps)

    K = decoder%vnum - decoder%cnum
    N_symb = decoder%vnum / bps
    allocate(x(N_symb))
    allocate(y(N_symb))
    allocate(word(decoder%vnum))
    allocate(synd(decoder%cnum))
    allocate(llr(decoder%vnum))
    allocate(llr_updated(decoder%vnum))


    call random_init(.false., .true.)
    call random_seed(get=seed)
    print *, "Original seed: ", seed
    call system_clock(time_seed)
    print *, "Time seed: ", time_seed, "their XOR: ", ieor(time_seed, seed)
    call random_seed(put=[ieor(time_seed, seed), this_image()])

    allocate(b_error_count(size(snrdb_array))[*])
    allocate(f_error_count(size(snrdb_array))[*])
    allocate(f_count(size(snrdb_array))[*])
    
    b_error_count(:) = 0
    f_error_count(:) = 0
    f_count(:) = 0
    sync all

    snr_loop : do i_snr = 1 , size(snrdb_array)
       N0 = 10.0_wp**(-snrdb_array(i_snr)/10.0_wp) * alphabet%variance
       N0_half = 0.5_wp * N0
       
       frame_loop : do i_frame = 1, max_loops
          call alphabet%random_symbol(x)
          word = alphabet%symbol_to_grey_array(x)
          synd = decoder%word_to_synd(word)
          
          y = rvs_normal(loc=alphabet%symbol_index_to_real(x), scale=N0_half)
          llr = y_to_llr_grey_array(N0, alphabet, y)

          n_ldpc_it = ldpc_iterations
          call decoder%decode(llr, llr_updated, synd, n_ldpc_it)

          new_errors = word_llr_errors(word(:K), llr_updated(:K))

          critical
            if (new_errors > 0) then
               b_error_count(i_snr)[1] = b_error_count(i_snr)[1] + new_errors
               f_error_count(i_snr)[1] = f_error_count(i_snr)[1] + 1
            end if
            f_count(i_snr)[1] = f_count(i_snr)[1] + 1
          end critical
          
          if ((f_error_count(i_snr)[1] > min_ferr) .and. &
               (f_count(i_snr) > min_loops)) then
             exit frame_loop
          end if
       end do frame_loop
    end do snr_loop


    sync all

    if (this_image() == 1) then
       ber = real(b_error_count, wp)/real(f_count*K, wp)
       fer = real(f_error_count, wp)/real(f_count, wp)
    end if
    
  end subroutine simulate_pam


  function y_to_llr_grey(N0, pa, y) result (llr)
    real(wp), intent(in) :: N0 ! Note that this is 2\sigma^2
    type(TAlphaPAM), intent(in)  :: pa
    real(wp), intent(in) :: y
    real(wp) :: llr(pa%B)
    real(wp) :: D(pa%B)

    integer :: M, i, b, x
    real(wp) :: addendum
    M = ishft(1, pa%B)

    llr(:) = 0.0_wp
    D(:)   = 0.0_wp
    
    do i = 0, M-1
       x = i
       addendum = exp(-((y-pa%constellation(i))**2.0_wp)/N0)
       do b = 1, pa%B
          if (iand(x*(x+1), b'11') == 0) then
             llr(b) = llr(b) + addendum
          else
             D(b) = D(b) + addendum
          end if
          x = ishft(x, -1)
       end do
    end do

    llr = log(llr) - log(D)
  end function y_to_llr_grey


  function y_to_llr_grey_array(N0, pa, y) result(llr)
    real(wp), intent(in) :: N0 ! Note that this is 2\sigma^2
    type(TAlphaPAM), intent(in)  :: pa
    real(wp), intent(in) :: y(0:)
    real(wp) :: llr(0:size(y)*pa%B-1)

    integer :: i

    do i = 0, size(y)-1
       llr(i*pa%B : (i+1)*pa%B-1) = y_to_llr_grey(N0, pa, y(i))
    end do
  end function y_to_llr_grey_array


  function word_llr_errors(word, llr) result (errors)
    logical, intent(in) :: word(:)
    real(wp), intent(in) :: llr(size(word))
    integer :: errors

    integer :: i
    errors = 0
    do i = 1, size(word)
       if (xor(word(i), llr(i) < 0.0_wp)) then
          errors = errors + 1
       end if
    end do
  end function word_llr_errors


end module sim_direct_channel