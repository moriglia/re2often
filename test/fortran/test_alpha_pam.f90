program test_alpha_pam
  use alpha_pam, only: TAlphaPAM
  use iso_c_binding, only: wp=>c_double
  implicit none
  
  type(TAlphaPAM) :: pa
  integer :: i
  integer :: x(10)


  pa = TAlphaPAM(2, [.2_wp, .3_wp, .3_wp, .2_wp], 2.0_wp)

  print *, pa%B, pa%M, pa%step, pa%variance, pa%probabilities

  print *, "symbol_to_grey output"
  do i = 0, 3
     print *, i, "=>", pa%symbol_to_grey(i)
  end do

  print *, "Map is:"
  do i = 0, 3
     print *, i, "=>", pa%symbol_to_bit_map(i,:)
  end do


  print *, ""
  print *, "Testing random symbols: "
  call pa%random_symbol(i)
  print *, i
  call pa%random_symbol(x)
  print *, x
  print *, pa%symbol_index_to_real(x)

  print *, "Bit sequence from random symbols: "
  print *, pa%symbol_to_grey_array(x)

  print *, "Default uniform probability"
  pa = TAlphaPAM(3, step=1.0_wp)
  print *, pa%M, pa%step, pa%variance, pa%probabilities

  print *, "All default but bps"
  pa = TAlphaPAM(4)
  print *, pa%M, pa%step, pa%variance, pa%probabilities
end program test_alpha_pam
