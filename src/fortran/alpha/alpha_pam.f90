module alpha_pam
  use iso_c_binding, only: wp => c_double
  
  type, public :: TAlphaPAM
     integer :: M ! modulation order
     integer :: B ! modulation bit / symbol
     real(wp) :: step
     real(wp), allocatable :: probabilities(:) ! 0 : M-1
     real(wp), allocatable :: constellation(:)
     real(wp) :: variance
     logical, allocatable :: symbol_to_bit_map(:,:) ! (0:M-1 , 0:B-1) 
   contains
     procedure, pass(this), public :: symbol_to_grey
     procedure, pass(this), public :: symbol_to_grey_array
     procedure, pass(this), public :: random_symbol
     procedure, pass(this), public :: symbol_index_to_real
  end type TAlphaPAM


  interface TAlphaPAM
     module procedure TAlphaPAMConstructor
  end interface TAlphaPAM


contains

  function TAlphaPAMConstructor(B, probabilities, step) result (this)
    integer, intent(in) :: B
    real(wp), intent(in) :: probabilities(0:ishft(1,B)-1)
    real(wp), intent(in), optional :: step
    type(TAlphaPAM) :: this

    integer :: i, M
    real(wp) :: base

    this%B = B
    M = ishft(1,B)
    this%M = M
    
    if (.not. present(step)) then
       this%step = 2.0_wp
    else
       this%step = step
    end if
    
    allocate(this%constellation(0:M-1))
    ! step*(real(((M_half-1)*M_half)/2 + ((M_half-1)*M_half*(M-1))/6, wp) + 0.25_wp)
    base = real(1-M, wp)*step*0.5_wp
    this%variance = 0.0_wp
    do i=0, M - 1
       this%constellation(i) = base + real(i, wp)*this%step
       this%variance = this%variance + probabilities(i) * this%constellation(i)**2
    end do

    allocate(this%probabilities(0:M-1))
    this%probabilities = probabilities

    allocate(this%symbol_to_bit_map(0:M-1,0:B-1))
    do i = 1, M
       this%symbol_to_bit_map(i, :) = this%symbol_to_grey(i-1)
    end do
  end function TAlphaPAMConstructor


  pure function symbol_to_grey(this, symbol) result (bits)
    ! Note that symbols go from 0 to M-1, rather than 1 to M
    class(TAlphaPAM), intent(in) :: this
    integer, intent(in) :: symbol
    logical :: bits(this%B)

    integer :: b, symb_tmp

    symb_tmp = symbol
    do b = 1, this%B
       bits(b) = iand(symb_tmp*(symb_tmp+1), b'11') /= 0
       symb_tmp = ishft(symb_tmp, -1)
    end do
  end function symbol_to_grey


  pure function symbol_to_grey_array(this, symbol_array) result (bits)
    class(TAlphaPAM), intent(in) :: this
    integer, intent(in) :: symbol_array(0:)
    logical :: bits(0:size(symbol_array)*this%B-1)

    integer :: i

    do i = 0,  size(symbol_array) - 1
       bits(i*this%B : (i+1)*this%B - 1) = this%symbol_to_bit_map(symbol_array(i),:)
    end do
  end function symbol_to_grey_array
    


  impure elemental subroutine random_symbol(this, x)
    class(TAlphaPAM), intent(in) :: this
    integer, intent(out) :: x

    real(wp) :: rnd
    integer :: i
    
    call random_number(rnd)
    
    x = this%M-1
    
    symbol_loop : do i = 0, this%M-2
       if (rnd < this%probabilities(i)) then
          x = i
          exit symbol_loop
       end if
       rnd = rnd - this%probabilities(i)
    end do symbol_loop
  end subroutine random_symbol

  
  elemental function symbol_index_to_real(this, idx) result (x)
    class(TAlphaPAM), intent(in) :: this
    integer, intent(in) :: idx

    real(wp) :: x

    x = this%constellation(idx)
  end function symbol_index_to_real
  
end module alpha_pam
