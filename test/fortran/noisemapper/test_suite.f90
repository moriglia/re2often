module test_suite_noise_mapper
  use noise_mapper, only: TNoiseMapper, binsearch, interpolate
  use iso_c_binding, only: wp=>c_double
  use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
  implicit none

  private

  public :: collect_suite


contains

  subroutine collect_suite(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
         new_unittest("construction without optionals", test_construction_without_optionals),&
         new_unittest("binsearch", test_binsearch),&
         new_unittest("hard decision", test_hard_decide_index),&
         new_unittest("interpolation", test_interpolate),&
         new_unittest("soft metric generation and usage", test_soft_metric)&
         ]
  end subroutine collect_suite
  
  
  subroutine test_construction_without_optionals(error)
    type(error_type), allocatable, intent(out) :: error

    integer :: i
    type(TNoiseMapper) :: nm
    nm = TNoiseMapper(B=4, sigma=0.5_wp)

    ! PAM Alphabet inherited parameters
    call check(error, nm%M, 16)
    if (allocated(error)) return

    call check(error, nm%B, 4)
    if (allocated(error)) return
    call check(error, nm%step, 2.0_wp, thr=0.001_wp)
    if (allocated(error)) return
    do i = 0, 15
       call check(error, nm%constellation(i), &
            real(real(-15, wp) + real(i, wp)*2.0_wp, wp), thr=0.001_wp)
       if (allocated(error)) return
       call check(error, nm%probabilities(i), 0.0625_wp, thr=0.000001_wp)
       if (allocated(error)) return
    end do

    call check(error, nm%variance, 85.0_wp, thr=0.01_wp)
    if (allocated(error)) return

    call check(error, nm%sigma, 0.5_wp, thr=0.001_wp)
    if (allocated(error)) return

    call check(error, nm%sigmaSquare, 0.25_wp, thr=0.001_wp)
    if (allocated(error)) return
    
  end subroutine test_construction_without_optionals


  subroutine test_binsearch(error)
    type(error_type), allocatable, intent(out) :: error

    call check(error, binsearch(real([3,4,5,6], wp), 0._wp) .eq. 0)
    if (allocated(error)) return

    call check(error, binsearch(real([3,4,5,6], wp), 6.4_wp) .eq. 3)
    if (allocated(error)) return

    call check(error, binsearch(real([3,4,5,6], wp), 4.2_wp) .eq. 1)
  end subroutine test_binsearch


  subroutine test_hard_decide_index(error)
    type(error_type), allocatable, intent(out) :: error

    integer :: i
    type(TNoiseMapper) :: nm
    nm = TNoiseMapper(B=3, sigma=0.5_wp, step=2.0_wp)

    call check(error, nm%hard_decide_index(-6.4_wp) .eq. 0)
    if (allocated(error)) return

    call check(error, all(nm%hard_decide_index([-6.4_wp, 5.9_wp, 7.2_wp, 0.1_wp, -5.9_wp]) .eq. &
         [0, 6, 7, 4, 1]))
    if (allocated(error)) return
  end subroutine test_hard_decide_index


  subroutine test_interpolate(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: i

    call check(error, &
         interpolate(&
         real([0,1,2,3,4], wp),&
         real([.2, .5, .6, .67, .9], wp), -1.4_wp), &
         0.2_wp, &
         thr=0.001_wp)
    if (allocated(error)) return

    call check(error, &
         interpolate(&
         real([0,1,2,3,4], wp),&
         real([.2, .5, .6, .67, .9], wp), 0.9_wp), &
         0.2_wp, &
         thr=0.001_wp)
    if (allocated(error)) return

    call check(error, &
         interpolate(&
         real([0,1,2,3,4], wp),&
         real([.2, .5, .6, .67, .9], wp), 5.2_wp), &
         0.9_wp, &
         thr=0.001_wp)
    if (allocated(error)) return

    call check(error, &
         interpolate(&
         real([0,1,2,3,4], wp),&
         real([.2, .5, .6, .67, .9], wp), 2.4_wp), &
         0.6_wp, &
         thr=0.001_wp)
  end subroutine test_interpolate


  subroutine test_soft_metric(error)
    type(error_type), allocatable, intent(out) :: error

    integer :: i
    real(wp) :: y, nhat
    type(TNoiseMapper) :: nm
    nm = TNoiseMapper(B=2, sigma=0.5_wp, step=2.0_wp)

    y = -3.2_wp
    i = nm%hard_decide_index(y)
    nhat = nm%generate_soft_metric(y, i)
    call check(error, nm%reconstruct_sample_from_metric(nhat, i), y, thr=0.002_wp)
    if (allocated(error)) return
    
    y = 3.2_wp
    i = nm%hard_decide_index(y)
    nhat = nm%generate_soft_metric(y, i)
    call check(error, nm%reconstruct_sample_from_metric(nhat, i), y, thr=0.002_wp)
    if (allocated(error)) return
    
    y = 1.1_wp
    i = nm%hard_decide_index(y)
    nhat = nm%generate_soft_metric(y, i)
    call check(error, nm%reconstruct_sample_from_metric(nhat, i), y, thr=0.002_wp)
  end subroutine test_soft_metric
end module test_suite_noise_mapper
