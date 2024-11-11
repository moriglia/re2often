module test_direct_channel_suite
  use sim_direct_channel, only: y_to_llr_grey, y_to_llr_grey_array, word_llr_errors
  use alpha_pam, only: TAlphaPAM
  use iso_c_binding, only: wp=>c_double
  use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
  implicit none

  private

  public :: collect_suite


contains

  subroutine collect_suite(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
         new_unittest("y_to_llr_grey", test_y_to_llr_grey),&
         new_unittest("y_to_llr_grey_array", test_y_to_llr_grey_array),&
         new_unittest("word_llr_errors", test_word_llr_errors)&
         ]
  end subroutine collect_suite


  subroutine test_y_to_llr_grey(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: N0
    type(TAlphaPAM) :: alphabet
    real(wp) :: y
    real(wp) :: llr(3)
    integer :: bps

    bps = 3
    alphabet = TAlphaPAM(bps)
    N0 = 1.0_wp
    y = -10.0_wp
    
    llr = y_to_llr_grey(N0, alphabet, y)
    

    call check(error, llr(1), 16.0_wp, thr=0.1_wp)
    if (allocated(error)) return

    call check(error, llr(2), 40.0_wp, thr=0.1_wp)
    if (allocated(error)) return

    call check(error, llr(3), 112.0_wp, thr=0.1_wp)
    if (allocated(error)) return
    
  end subroutine test_y_to_llr_grey


  subroutine test_y_to_llr_grey_array(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: N0
    type(TAlphaPAM) :: alphabet
    real(wp) :: y(3)
    real(wp) :: llr(9)
    integer :: bps

    bps = 3
    alphabet = TAlphaPAM(bps)
    N0 = 1.0_wp
    y = [-10.0_wp, 0.0_wp, 3.1_wp]
    
    llr = y_to_llr_grey_array(N0, alphabet, y)
    

    call check(error, llr(1), 16.0_wp, thr=0.01_wp)
    if (allocated(error)) return

    call check(error, llr(2), 40.0_wp, thr=0.01_wp)
    if (allocated(error)) return

    call check(error, llr(3), 112.0_wp, thr=0.01_wp)
    if (allocated(error)) return

    call check(error, llr(4), 8.0_wp, thr=0.01_wp)
    if (allocated(error)) return

    call check(error, llr(5), -24.0_wp, thr=0.01_wp)
    if (allocated(error)) return

    call check(error, llr(6), 0.0_wp, thr=0.000001_wp)
    if (allocated(error)) return

    call check(error, llr(7), -4.4_wp, thr=0.05_wp)
    if (allocated(error)) return

    call check(error, llr(8), -3.6_wp, thr=0.05_wp)
    if (allocated(error)) return

    call check(error, llr(9), -16.8_wp, thr=0.05_wp)
    
  end subroutine test_y_to_llr_grey_array


  subroutine test_word_llr_errors(error)
    type(error_type), allocatable, intent(out) :: error

    real(wp) :: llr(10)
    logical  :: word(10)


    llr = [-0.2_wp, 0.3_wp, 9.2_wp, 0.44_wp, -3.2_wp, 0.22_wp, -9.23_wp, 32.0_wp, 3.1_wp, -4.4_wp]
    word = [.true., .false., .false., .true., .false., .true., .false., .true., .true., .false.]

    call check(error, word_llr_errors(word, llr), 7)
  end subroutine test_word_llr_errors

end module test_direct_channel_suite
