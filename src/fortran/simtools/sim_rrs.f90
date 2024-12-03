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
module sim_rrs
  use iso_c_binding, only: wp => c_double, ki => c_int, kl => c_bool
  use ldpc_decoder, only: TDecoder
  use noise_mapper, only: TNoiseMapper
  use stdlib_stats_distribution_normal, only: rvs_normal
  use forbear, only: bar_object
  implicit none

  public :: simulate_reverse_pam

contains
  
  subroutine simulate_reverse_pam (snrdb_array, Nsnr, bps, &
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
    type(TNoiseMapper) :: noisemapper
    integer :: time_seed, seed(8)

    integer :: i_snr, i_frame
    integer :: N_symb, K

    integer, allocatable :: x(:)
    real(wp), allocatable :: y(:)
    real(wp), allocatable :: nhat(:)
    integer, allocatable :: xhat(:)
    logical, allocatable :: word(:)
    logical, allocatable :: synd(:)
    real(wp), allocatable :: lappr(:)
    real(wp), allocatable :: lappr_updated(:)

    real(wp) :: N0, N0_half, sigma

    integer :: new_errors, n_ldpc_it

    integer, allocatable :: b_error_count(:)[:]
    integer, allocatable :: f_error_count(:)[:]
    integer, allocatable :: f_count(:)[:]

    type(bar_object) :: progress_bar
    
    
    decoder = TDecoder(Ne, e_to_v, e_to_c)
    noisemapper = TNoiseMapper(B=bps, sigma=1.0_wp)

    K = decoder%vnum - decoder%cnum
    N_symb = decoder%vnum / bps
    allocate(x(N_symb))
    allocate(y(N_symb))
    allocate(xhat(N_symb))
    allocate(nhat(N_symb))
    allocate(word(decoder%vnum))
    allocate(synd(decoder%cnum))
    allocate(lappr(decoder%vnum))
    allocate(lappr_updated(decoder%vnum))


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
       N0 = 10.0_wp**(-snrdb_array(i_snr)/10.0_wp) * noisemapper%variance
       N0_half = 0.5_wp * N0
       sigma = sqrt(N0_half)

       call noisemapper%update_noise_sigma(sigma)

       frame_loop : do i_frame = 1, max_loops
          call noisemapper%random_symbol(x)
          word = noisemapper%symbol_to_grey_array(x)
          synd = decoder%word_to_synd(word)
          
          y = rvs_normal(loc=noisemapper%symbol_index_to_real(x), scale=sigma)

          xhat = noisemapper%hard_decide_index(y)
          nhat = noisemapper%generate_soft_metric(y, xhat)
          lappr = noisemapper%demap_metric_to_lappr(nhat, x)
          
          n_ldpc_it = ldpc_iterations
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
               (f_count(i_snr)[1] > min_loops)) then
             exit frame_loop
          end if
       end do frame_loop
       if (this_image() == 1) then
          call progress_bar%update(current=real(i_snr, 8)/real(Nsnr, 8))
       end if
    end do snr_loop


    sync all

    if (this_image() == 1) then
       ber = real(b_error_count, wp)/real(f_count*K, wp)
       fer = real(f_error_count, wp)/real(f_count, wp)
    end if
    
  end subroutine simulate_reverse_pam


  pure function word_lappr_errors(word, lappr) result (errors)
    logical, intent(in) :: word(:)
    real(wp), intent(in) :: lappr(size(word))
    integer :: errors

    integer :: i
    errors = 0
    do i = 1, size(word)
       if (word(i) .neqv. (lappr(i) < 0.0_wp)) then
          errors = errors + 1
       end if
    end do
  end function word_lappr_errors


end module sim_rrs
