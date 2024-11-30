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

  public :: TNoiseMapper

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

     ! Monotonicity Configuration
     logical(kb), allocatable :: monotonicity_config(:)

   contains
     procedure, pass(this), public :: set_y_thresholds

     procedure, pass, public :: y_to_cdf
     procedure, pass, public :: y_to_cdf_array
     generic, public :: CDF_Y => y_to_cdf, y_to_cdf_array
  end type TNoiseMapper

  interface TNoiseMapper
     module procedure TNoiseMapperConstructor
  end interface TNoiseMapper

  ! interface F_Y
  !    module procedure y_to_cdf
  !    module procedure y_to_cdf_array
  ! end interface F_Y

contains
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
    allocate(this%F_Y_thresholds(1:this%M-1))
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
       n_range = ceiling((y_high-y_low)*intervals_per_step/this%step)
    else
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
    this%F_Y_thresholds = this%CDF_Y(y_thresholds)

    this%delta_F_Y(1:this%M-2) = this%F_Y_thresholds(2:this%M-1) - this%F_Y_thresholds(1:this%M-2)
    this%delta_F_Y(0) = this%F_Y_thresholds(1)
    this%delta_F_Y(this%M-1) = 1.0_wp - this%F_Y_thresholds(this%M-1)
  end subroutine set_y_thresholds


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
  
end module noise_mapper
