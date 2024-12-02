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
module noise_mapper
  use iso_c_binding, only: wp => c_double, kb => c_bool
  use alpha_pam, only: TAlphaPAM
  use stdlib_stats_distribution_normal, only: cdf_normal
  use stdlib_kinds, only: dp
  implicit none

  public :: TNoiseMapper, binsearch, interpolate

  type, public, extends(TAlphaPAM) :: TNoiseMapper
     ! noise parameters
     real(wp) :: sigma
     real(wp) :: sigmaSquare

     ! Threshold parameters
     real(wp), allocatable :: y_thresholds(:)
     real(wp), allocatable :: F_Y_thresholds(:)
     real(wp), allocatable :: delta_F_Y(:)

     ! Pre-calculated parameters
     real(wp), allocatable :: y_range(:)
     real(wp), allocatable :: F_Y_range(:)
     integer :: intervals_per_step

     ! Monotonicity Configuration
     logical(kb), allocatable :: monotonicity_config(:)

   contains
     procedure, pass(this), public :: set_y_thresholds

     procedure, pass, private :: y_to_cdf
     procedure, pass, private:: y_to_cdf_array
     generic, public :: CDF_Y => y_to_cdf, y_to_cdf_array

     procedure, pass, public :: hard_decide_index

     procedure, pass, private :: generate_soft_metric_single
     procedure, pass, private :: generate_soft_metric_array
     generic, public :: generate_soft_metric => &
          generate_soft_metric_single, generate_soft_metric_array
     procedure, pass, public :: reconstruct_sample_from_metric

     procedure, pass, private :: demap_metric_to_lappr_single_transmission
     procedure, pass, private :: demap_metric_to_lappr_array
     generic, public :: demap_metric_to_lappr => &
          demap_metric_to_lappr_single_transmission, demap_metric_to_lappr_array

     procedure, pass, public  :: update_noise_sigma
     procedure, pass, public :: free
  end type TNoiseMapper

  interface TNoiseMapper
     module procedure TNoiseMapperConstructor
  end interface TNoiseMapper
contains

  subroutine free(this)
    class(TNoiseMapper), intent(inout) :: this

    if (allocated(this%y_thresholds)) deallocate(this%y_thresholds)
    if (allocated(this%F_Y_thresholds)) deallocate(this%F_Y_thresholds)
    if (allocated(this%delta_F_Y)) deallocate(this%delta_F_Y)
    if (allocated(this%monotonicity_config)) deallocate(this%monotonicity_config)
    if (allocated(this%y_range))   deallocate(this%y_range) 
    if (allocated(this%F_Y_range)) deallocate(this%F_Y_range)
  end subroutine free

  
  function TNoiseMapperConstructor(B, probabilities, step, & ! Alphabet parameters
       sigma, monotonicity_config, intervals_per_step) &     ! Noise Mapper parameters
       result (this)
    integer, intent(in) :: B
    real(wp), intent(in), optional :: probabilities(0:ishft(1,B)-1)
    real(wp), intent(in), optional :: step
    
    real(wp), intent(in) :: sigma
    logical(kb), intent(in), optional :: monotonicity_config(0:ishft(1, B)-1)
    integer, intent(in), optional :: intervals_per_step

    type(TNoiseMapper) :: this

    integer :: i, n_range
    real(wp) :: y_low, y_high

    call this%TAlphaPAM%free
    call this%free
    

    if (present(probabilities)) then
       if (present(step)) then
          this%TAlphaPAM = TAlphaPAM(B, probabilities, step)
       else
          this%TAlphaPAM = TAlphaPAM(B, probabilities=probabilities)
       end if
    else
       if (present(step)) then
          this%TAlphaPAM = TAlphaPAM(B, step=step)
       else
          this%TAlphaPAM = TAlphaPAM(B)
       end if
    end if

    this%sigma = sigma
    this%sigmaSquare = sigma**2

    allocate(this%y_thresholds(1:this%M-1))
    allocate(this%F_Y_thresholds(0:this%M))
    allocate(this%delta_F_Y(0:this%M-1))

    call this%set_y_thresholds([(real(i+1-ishft(this%M, -1), wp)*this%step, i = 0, this%M-2)])

    allocate(this%monotonicity_config(0:this%M-1))
    if (present(monotonicity_config)) then
       this%monotonicity_config = monotonicity_config
    else
       this%monotonicity_config(0::2) = .false.
       this%monotonicity_config(1::2) = .true.
    end if


    y_high = this%constellation(this%M-1) + sigma*sqrt(-2.0_wp*log(0.01_wp))
    y_low  = - y_high
    if (present(intervals_per_step)) then
       this%intervals_per_step = intervals_per_step
       n_range = ceiling((y_high-y_low)*intervals_per_step/this%step)
    else
       this%intervals_per_step = 1000
       n_range = ceiling((y_high-y_low)*1000/this%step)
    end if
    allocate(this%y_range(0:n_range))
    allocate(this%F_Y_range(0:n_range))
    this%y_range = [(y_low + real(i,wp)*(y_high-y_low)/n_range, i=0, n_range)]
    this%F_Y_range = this%CDF_Y(this%y_range)
  end function TNoiseMapperConstructor


  subroutine set_y_thresholds(this, y_thresholds)
    class(TNoiseMapper) :: this
    real(wp), intent(in) :: y_thresholds(1:this%M-1)

    integer :: i
    
    this%y_thresholds = y_thresholds
    this%F_Y_thresholds(1:this%M-1) = this%CDF_Y(y_thresholds)

    this%F_Y_thresholds(0) = 0.0_wp
    this%F_Y_thresholds(this%M) = 1.0_wp
    
    this%delta_F_Y(1:this%M-1) = this%F_Y_thresholds(2:this%M) - this%F_Y_thresholds(1:this%M-1)
    this%delta_F_Y(0) = this%F_Y_thresholds(1)
  end subroutine set_y_thresholds


  subroutine update_noise_sigma(this, sigma)
    class(TNoiseMapper), intent(inout) :: this

    real(wp), intent(in) :: sigma

    real(wp) :: y_high, y_low
    integer  :: n_range, i

    y_high = this%constellation(this%M-1) + sigma*sqrt(-2.0_wp*log(0.01_wp))
    y_low  = - y_high
    n_range = ceiling((y_high-y_low)*this%intervals_per_step/this%step)

    if (allocated(this%y_range)) deallocate(this%y_range)
    if (allocated(this%F_Y_range)) deallocate(this%F_Y_range)

    allocate(this%y_range(0:n_range))
    allocate(this%F_Y_range(0:n_range))
    
    this%y_range = [(y_low + real(i,wp)*(y_high-y_low)/n_range, i=0, n_range)]
    this%F_Y_range = this%CDF_Y(this%y_range)

    call this%set_y_thresholds(this%y_thresholds)
  end subroutine update_noise_sigma

  function y_to_cdf(this, y) result (F)
    class(TNoiseMapper) :: this
    real(wp), intent(in) :: y

    real(wp) :: F

    F = sum(this%probabilities * cdf_normal(y, loc=this%constellation, scale=this%sigma))
  end function y_to_cdf


  function y_to_cdf_array(this, y) result (F)
    class(TNoiseMapper):: this
    real(wp), intent(in) :: y(:)

    real(wp) :: F(size(y))

    integer :: i

    F(:) = 0.0_wp

    do i = 1, size(y)
       F(i) = this%CDF_Y(y(i))
    end do
  end function y_to_cdf_array


  pure function binsearch(domain, val) result(ind)
    real(wp), intent(in) :: domain(0:)
    real(wp), intent(in) :: val
    integer :: ind

    integer :: upper, lower
    upper = size(domain) ! exclusive upper bound
    lower = 0 ! inclusive lower bound

    if (domain(upper-1) <= val) then
       ind = upper - 1
       return
    end if

    if (domain(0) > val) then
       ind = 0
       return
    end if
    
    ind = ishft(upper, -1)

    do while (upper > lower+1)
       if (val < domain(ind)) then
          upper = ind
       else
          lower = ind
       end if
       ind = (upper + lower)/2
    end do
  end function binsearch


  pure function interpolate(domain, codomain, d_val) result(c_val)
    real(wp), intent(in) :: domain(0:)
    real(wp), intent(in) :: codomain(0:)
    real(wp), intent(in) :: d_val

    real(wp) :: c_val

    c_val = codomain(binsearch(domain, d_val))
    ! Possible improvement: linear interpolation codomain(i) and (i+1)
  end function interpolate


  elemental function hard_decide_index(this, y) result (x_ind)
    class(TNoiseMapper), intent(in) :: this
    real(wp), intent(in) :: y
    integer :: x_ind

    if (y < this%y_thresholds(1)) then
       x_ind = 0
       return
    end if
    
    x_ind = binsearch(this%y_thresholds(:), y) + 1
  end function hard_decide_index


  function generate_soft_metric_single(this, y, i) result (nhat)
    class(TNoiseMapper), intent(in) :: this
    real(wp), intent(in) :: y
    integer, intent(in)  :: i

    real(wp) :: nhat

    if (this%monotonicity_config(i)) then
       nhat = (this%F_Y_thresholds(i+1) -  this%CDF_Y(y))/this%delta_F_Y(i)
    else
       nhat = (this%CDF_Y(y) - this%F_Y_thresholds(i))/this%delta_F_Y(i)
    end if
  end function generate_soft_metric_single


  function generate_soft_metric_array(this, y, i) result (nhat)
    class(TNoiseMapper), intent(in) :: this
    real(wp), intent(in) :: y(:)
    integer, intent(in)  :: i(size(y))

    real(wp) :: nhat(size(y))

    integer :: k

    do k = 1, size(y)
       nhat(k) = this%generate_soft_metric_single(y(k), i(k))
    end do
  end function generate_soft_metric_array

  
  function reconstruct_sample_from_metric(this, nhat, i) result(y)
    class(TNoiseMapper), intent(in) :: this
    real(wp), intent(in) :: nhat
    integer, intent(in)  :: i

    real(wp) :: y

    if (this%monotonicity_config(i)) then
       y = interpolate(this%F_Y_range, this%y_range, &
            this%F_Y_thresholds(i+1) -  nhat*this%delta_F_Y(i))
    else
       y = interpolate(this%F_Y_range, this%y_range, &
            this%F_Y_thresholds(i) +  nhat*this%delta_F_Y(i))
    end if
  end function reconstruct_sample_from_metric


  function demap_metric_to_lappr_single_transmission(this, nhat, j) result (lappr)
    class(TNoiseMapper) :: this
    integer, intent(in) :: j      ! transmitted symbol
    real(wp), intent(in):: nhat   ! softening metric

    real(wp) :: lappr(0:this%B-1)

    real(wp) :: N(0:this%B-1)
    real(wp) :: D(0:this%B-1)

    integer :: l ! bit index within symbol
    integer :: i ! possible received symbol
    integer :: k ! possible transmitted symbols

    real(wp) :: a_j, y_i, addendum, twoSigmaSquare

    N(:) = 0.0_wp
    D(:) = 0.0_wp

    a_j = this%constellation(j)
    twoSigmaSquare = this%sigmaSquare * 2.0_wp
    do i = 0, this%M
       y_i = this%reconstruct_sample_from_metric(nhat, i)
       addendum = 0.0_wp
       do k = 0, this%M-1
          addendum = addendum + this%probabilities(k) * &
               exp( (this%constellation(k) - a_j) * &
               (2*y_i - this%constellation(k) - a_j) / twoSigmaSquare)
       end do

       do l = 0, this%B-1
          if (this%symbol_to_bit_map(i, l)) then
             D(l) = D(l) + addendum
          else
             N(l) = N(l) + addendum
          end if
       end do
    end do

    lappr = log(N) - log(D)
  end function demap_metric_to_lappr_single_transmission


  function demap_metric_to_lappr_array(this, nhat, j) result (lappr)
    class(TNoiseMapper) :: this
    integer, intent(in) :: j(0:)     ! transmitted symbol
    real(wp), intent(in):: nhat(0:size(j)-1)   ! softening metric

    real(wp) :: lappr(0:size(j)*this%B-1)

    integer :: k

    do k = 0, size(j)-1
       lappr( k*this%B : (k+1)*this%B - 1 ) = this%demap_metric_to_lappr_single_transmission(nhat(k), j(k))
    end do
  end function demap_metric_to_lappr_array
end module noise_mapper
