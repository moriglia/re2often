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
  use forbear, only: bar_object
  implicit none

  public :: simulate_pam, y_to_llr_grey, y_to_llr_grey_array, word_llr_errors

contains
  
  subroutine simulate_pam (snrdb_array, Nsnr, bps, &
       min_ferr, max_loops, min_loops, &
       e_to_v, e_to_c, Ne, ldpc_iterations, &
       ber, fer) bind (C)
    integer(ki), intent(in)  :: Nsnr
    real(wp), intent(in)     :: snrdb_array(Nsnr)
    integer(ki), intent(in)  :: bps, min_ferr, max_loops, min_loops
    integer(ki), intent(in)  :: Ne
    integer(ki), intent(in)  :: e_to_v(Ne)
    integer(ki), intent(in)  :: e_to_c(Ne)
    integer(ki), intent(in)  :: ldpc_iterations
    real(wp), intent(out)    :: ber(Nsnr)
    real(wp), intent(out)    :: fer(Nsnr)

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

    real(wp) :: N0, N0_half, sigma

    integer :: new_errors, n_ldpc_it

    integer, allocatable :: b_error_count(:)[:]
    integer, allocatable :: f_error_count(:)[:]
    integer, allocatable :: f_count(:)[:]

    type(bar_object) :: progress_bar
    
    
    decoder = TDecoder(Ne, e_to_v, e_to_c)
    alphabet = TAlphaPAM(bps)

    K = decoder%vnum - decoder%cnum
    N_symb = decoder%vnum / bps
    allocate(x(N_symb))
    allocate(y(N_symb))
    allocate(word(decoder%vnum))
    allocate(synd(decoder%cnum))
    allocate(llr(decoder%vnum))
    allocate(llr_updated(decoder%vnum))


    ! call random_init(.false., .true.)
    call random_seed(get=seed)
    call sleep(this_image())
    call system_clock(time_seed)
    call random_seed(put=[ieor(time_seed, seed), this_image()])

    allocate(b_error_count(Nsnr)[*])
    allocate(f_error_count(Nsnr)[*])
    allocate(f_count(Nsnr)[*])

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
    
    snr_loop : do i_snr = 1 , Nsnr
       N0 = 10.0_wp**(-snrdb_array(i_snr)/10.0_wp) * alphabet%variance
       N0_half = 0.5_wp * N0
       sigma = sqrt(N0_half)
       
       frame_loop : do i_frame = 1, max_loops
          call alphabet%random_symbol(x)
          word = alphabet%symbol_to_grey_array(x)
          synd = decoder%word_to_synd(word)
          
          y = rvs_normal(loc=alphabet%symbol_index_to_real(x), scale=sigma)
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
       if (this_image()==1) then
          call progress_bar%update(current=real(i_snr, 8)/real(Nsnr, 8))
       end if
    end do snr_loop


    sync all

    if (this_image() == 1) then
       ber = real(b_error_count, wp)/real(f_count*K, wp)
       fer = real(f_error_count, wp)/real(f_count, wp)
    end if
    
  end subroutine simulate_pam


  pure function y_to_llr_grey(N0, pa, y) result (llr)
    real(wp), intent(in) :: N0 ! Note that this is 2\sigma^2
    type(TAlphaPAM), intent(in)  :: pa
    real(wp), intent(in) :: y
    real(wp) :: llr(0:pa%B-1)
    real(wp) :: D(0:pa%B-1)

    integer :: M, i, b
    real(wp) :: addendum
    M = ishft(1, pa%B)

    llr(:) = 0.0_wp
    D(:)   = 0.0_wp
    
    do i = 0, M-1
       addendum = pa%probabilities(i) * exp(-((y-pa%constellation(i))**2.0_wp)/N0)
       do b = 0, pa%B - 1
          if (pa%symbol_to_bit_map(i,b)) then
             D(b) = D(b) + addendum
          else
             llr(b) = llr(b) + addendum
          end if
       end do
    end do

    llr = log(llr) - log(D)
  end function y_to_llr_grey


  pure function y_to_llr_grey_array(N0, pa, y) result(llr)
    real(wp), intent(in) :: N0 ! Note that this is 2\sigma^2
    type(TAlphaPAM), intent(in)  :: pa
    real(wp), intent(in) :: y(0:)
    real(wp) :: llr(0:size(y)*pa%B-1)

    integer :: i

    do i = 0, size(y)-1
       llr(i*pa%B : (i+1)*pa%B-1) = y_to_llr_grey(N0, pa, y(i))
    end do
  end function y_to_llr_grey_array


  pure function word_llr_errors(word, llr) result (errors)
    logical, intent(in) :: word(:)
    real(wp), intent(in) :: llr(size(word))
    integer :: errors

    integer :: i
    errors = 0
    do i = 1, size(word)
       if (word(i) .neqv. (llr(i) < 0.0_wp)) then
          errors = errors + 1
       end if
    end do
  end function word_llr_errors


end module sim_direct_channel
