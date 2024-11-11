module test_ldpc_edge_list
  use ldpc_edge_list, only: TEdgeList
  use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
  implicit none

  private

  public :: collect_suite


contains

  subroutine collect_suite(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [&
         new_unittest("TEdgelistConstructor", test_ldpc_edge_list_construction)&
         ]
  end subroutine collect_suite


  subroutine test_ldpc_edge_list_construction(error)
    type(error_type), allocatable, intent(out) :: error

    type(TEdgeList) :: el

    el = TEdgeList(5)

    call check(error, el%N, 5)
    if (allocated(error)) return

    call check(error, allocated(el%data))
    if (allocated(error)) return

    el%data(5) = 3 ! Should not segfault


    el = TEdgeList(9)
    call check(error, el%N, 9)
    if (allocated(error)) return

    call check(error, allocated(el%data))
    if (allocated(error)) return

    el%data(9) = 0 ! Should not segfault
    
    
  end subroutine test_ldpc_edge_list_construction

end module test_ldpc_edge_list
