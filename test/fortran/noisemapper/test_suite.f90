module test_suite_noise_mapper
  use noise_mapper, only: TNoiseMapper
  use iso_c_binding, only: wp=>c_double
  use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
  implicit none

  private

  public :: collect_suite


contains

  subroutine collect_suite(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
         ! new_unittest("only bps", test_bps_only),&
         ! new_unittest("bps + probabilities", test_bps_probabilities),&
         ! new_unittest("bps + step", test_bps_step),&
         ! new_unittest("bps + step + probabilities", test_bps_step_probabilities),&
         ! new_unittest("symbol_to_bit_map", test_symbol_to_bit_map),&
         ! new_unittest("symbol_to_grey", test_symbol_to_grey),&
         ! new_unittest("symbol_to_grey_array", test_symbol_to_grey_array),&
         new_unittest("construction without optionals", test_construction_without_optionals)&
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
       call check(error, nm%constellation(i), real(real(-15, wp) + real(i, wp)*2.0_wp, wp), thr=0.001_wp)
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



  ! subroutine test_bps_probabilities(error)
  !   type(error_type), allocatable, intent(out) :: error

  !   integer :: i
  !   type(TAlphaPAM) :: pa
  !   pa = TAlphaPAM(2, [0.1_wp, 0.4_wp, 0.4_wp, 0.1_wp])

  !   call check(error, pa%M, 4)
  !   if (allocated(error)) return

  !   call check(error, pa%B, 2)
  !   if (allocated(error)) return
  !   call check(error, pa%step, 2.0_wp, thr=0.001_wp)
  !   if (allocated(error)) return
  !   do i = 0, 3
  !      call check(error, pa%constellation(i), real(real(-3, wp) + real(i, wp)*2.0_wp, wp), thr=0.001_wp)
  !      if (allocated(error)) return
  !   end do

  !   call check(error, pa%probabilities(0), 0.1_wp, thr=0.000001_wp)
  !   if (allocated(error)) return
  !   call check(error, pa%probabilities(1), 0.4_wp, thr=0.000001_wp)
  !   if (allocated(error)) return
  !   call check(error, pa%probabilities(2), 0.4_wp, thr=0.000001_wp)
  !   if (allocated(error)) return
  !   call check(error, pa%probabilities(3), 0.1_wp, thr=0.000001_wp)
  !   if (allocated(error)) return
    
  !   call check(error, pa%variance, 2.6_wp, thr=0.001_wp)
  !   if (allocated(error)) return
  ! end subroutine test_bps_probabilities


  ! subroutine test_bps_step(error)
  !   type(error_type), allocatable, intent(out) :: error

  !   integer :: i
  !   type(TAlphaPAM) :: pa
  !   pa = TAlphaPAM(2, step=1.0_wp)

  !   call check(error, pa%M, 4)
  !   if (allocated(error)) return

  !   call check(error, pa%B, 2)
  !   if (allocated(error)) return
  !   call check(error, pa%step, 1.0_wp, thr=0.001_wp)
  !   if (allocated(error)) return
  !   do i = 0, 3
  !      call check(error, pa%constellation(i), real(real(-3, wp)/2.0_wp + real(i, wp), wp), thr=0.001_wp)
  !      if (allocated(error)) return

  !      call check(error, pa%probabilities(i), 0.25_wp, thr=0.000001_wp)
  !      if (allocated(error)) return
  !   end do

    
  !   call check(error, pa%variance, 1.25_wp, thr=0.001_wp)
  !   if (allocated(error)) return
  ! end subroutine

  ! subroutine test_bps_step_probabilities(error)
  !   type(error_type), allocatable, intent(out) :: error

  !   integer :: i
  !   type(TAlphaPAM) :: pa
  !   pa = TAlphaPAM(2, [0.1_wp, 0.4_wp, 0.4_wp, 0.1_wp], 1.0_wp)

  !   call check(error, pa%M, 4)
  !   if (allocated(error)) return

  !   call check(error, pa%B, 2)
  !   if (allocated(error)) return
  !   call check(error, pa%step, 1.0_wp, thr=0.001_wp)
  !   if (allocated(error)) return
  !   do i = 0, 3
  !      call check(error, pa%constellation(i), real(-3, wp)/2.0_wp + real(i, wp), thr=0.001_wp)
  !      if (allocated(error)) return
  !   end do

  !   call check(error, pa%probabilities(0), 0.1_wp, thr=0.000001_wp)
  !   if (allocated(error)) return
  !   call check(error, pa%probabilities(1), 0.4_wp, thr=0.000001_wp)
  !   if (allocated(error)) return
  !   call check(error, pa%probabilities(2), 0.4_wp, thr=0.000001_wp)
  !   if (allocated(error)) return
  !   call check(error, pa%probabilities(3), 0.1_wp, thr=0.000001_wp)
  !   if (allocated(error)) return
    
  !   call check(error, pa%variance, 0.65_wp, thr=0.001_wp)
  !   if (allocated(error)) return
  ! end subroutine test_bps_step_probabilities


  ! subroutine test_symbol_to_bit_map(error)
  !   type(error_type), allocatable, intent(out) :: error
  !   logical :: table(0:7, 0:2)
  !   integer :: i, j
  !   type(TAlphaPAM) :: pa
  !   pa = TAlphaPAM(3)

    

  !   table(:3,2) = .false.
  !   table(4:,2) = .true.

  !   table(:, 1) = [.false., .false., .true., .true.,.true., .true., .false., .false.]
  !   table(:, 0) = [.false., .true., .true., .false.,.false., .true., .true., .false.]


  !   do i=0, pa%M-1
  !      do j = 0, pa%B-1
  !         call check(error, pa%symbol_to_bit_map(i, j) .eqv. table(i,j))
  !         if (allocated(error)) return
  !      end do
  !   end do
  ! end subroutine test_symbol_to_bit_map


  ! subroutine test_symbol_to_grey(error)
  !   type(error_type), allocatable, intent(out) :: error
  !   logical :: table(0:7, 0:2)
  !   integer :: i
  !   type(TAlphaPAM) :: pa
  !   pa = TAlphaPAM(3)

  !   table(:3,2) = .false.
  !   table(4:,2) = .true.

  !   table(:, 1) = [.false., .false., .true., .true.,.true., .true., .false., .false.]
  !   table(:, 0) = [.false., .true., .true., .false.,.false., .true., .true., .false.]


  !   do i=0, pa%M-1
  !      call check(error, all(pa%symbol_to_grey(i) .eqv. table(i,:)))
  !      if (allocated(error)) return
  !   end do
  ! end subroutine test_symbol_to_grey


  ! subroutine test_symbol_to_grey_array(error)
  !   type(error_type), allocatable, intent(out) :: error
  !   integer :: i
  !   integer :: symbols(3)
  !   type(TAlphaPAM) :: pa
  !   pa = TAlphaPAM(3)

  !   symbols = [0,1,7]

  !   do i=0, pa%M-1
  !      call check(error, all(&
  !           pa%symbol_to_grey_array(symbols) .eqv. &
  !           [.false., .false., .false., .true., .false., .false., .false., .false., .true.]))
  !      if (allocated(error)) then
  !         print*, pa%symbol_to_grey_array(symbols)
  !         return
  !      end if
       
  !   end do
  ! end subroutine test_symbol_to_grey_array


  ! subroutine test_symbol_index_to_real(error)
  !   type(error_type), allocatable, intent(out) :: error
  !   integer :: i
  !   integer :: symbols(8)
  !   type(TAlphaPAM) :: pa
  !   pa = TAlphaPAM(3)

  !   symbols = [0,1,3,4,2,6,7,5]

  !   call check(error, all(abs(&
  !        pa%symbol_index_to_real(symbols) - &
  !        [-7.0_wp, -5.0_wp, -1.0_wp, 1.0_wp, -3.0_wp, 5.0_wp, 7.0_wp, 3.0_wp]) < 0.001_wp))

  ! end subroutine test_symbol_index_to_real
  
end module test_suite_noise_mapper
