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
module ldpc_decoder
  use ldpc_edge_list, only: TEdgeList
  use iso_c_binding, only: wp => c_double
  implicit none

  public :: TDecoder

  type, public :: TDecoder
     integer :: cnum
     integer :: vnum
     integer :: Ne
     type(TEdgeList), allocatable :: c_to_e(:)
     type(TEdgeList), allocatable :: v_to_e(:)
     type(TEdgeList), allocatable :: c_to_v(:)
     type(TEdgeList), allocatable :: v_to_c(:)
     real(wp), allocatable :: B_buffer(:)
     real(wp), allocatable :: F_buffer(:)
   contains
     final :: destructor
     procedure, pass(this), public :: print => print_decoder
     procedure, pass(this), public :: check_llr
     procedure, pass(this), public :: word_to_synd
     procedure, pass(this), public :: decode
  end type TDecoder

  interface TDecoder
     module procedure TDecoderConstructor
  end interface TDecoder

contains
  ! ------------------------------
  ! --- Construction and setup --- 
  ! ------------------------------
  function TDecoderConstructor(N, e_to_v, e_to_c) result(decoder)
    type(TDecoder) :: decoder
    integer, intent(in) :: N
    integer, intent(in) :: e_to_c(N)
    integer, intent(in) :: e_to_v(N)
    integer :: j, i, max_buff

    decoder%cnum = maxval(e_to_c) + 1
    decoder%vnum = maxval(e_to_v) + 1
    decoder%Ne   = N

    allocate(decoder%c_to_e(decoder%cnum))
    allocate(decoder%v_to_e(decoder%vnum))
    call invert_table(N, e_to_c, decoder%cnum, decoder%c_to_e)
    call invert_table(N, e_to_v, decoder%vnum, decoder%v_to_e)

    allocate(decoder%c_to_v(decoder%cnum))
    allocate(decoder%v_to_c(decoder%vnum))
    call edge_to_node_convert_table(decoder%cnum, decoder%c_to_e, N, e_to_v, decoder%c_to_v)
    call edge_to_node_convert_table(decoder%vnum, decoder%v_to_e, N, e_to_c, decoder%v_to_c)

    max_buff = 0
    do i = 1, decoder%cnum
       max_buff = max(max_buff, decoder%c_to_e(i)%N)
    end do
    do i = 1, decoder%vnum
       max_buff = max(max_buff, decoder%v_to_e(i)%N)
    end do
    allocate(decoder%B_buffer(2:max_buff))
    allocate(decoder%F_buffer(max_buff-1))
  end function TDecoderConstructor

  subroutine invert_table(N, e_to_x, xnum, x_to_e)
    integer, intent(in) :: N
    integer, intent(in) :: e_to_x(N)
    integer, intent(in) :: xnum
    type(TEdgeList), intent(inout) :: x_to_e(xnum)

    integer :: xcount(xnum)
    integer :: i, j

    xcount(:) = 0
    do i = 1, N
       xcount(e_to_x(i)+1) = xcount(e_to_x(i)+1) + 1
    end do

    do i = 1, xnum
       x_to_e(i) = TEdgeList(xcount(i))
    end do

    do i = N, 1, -1
       j = e_to_x(i) + 1
       x_to_e(j)%data(xcount(j)) = i
       xcount(j) = xcount(j) - 1
    end do
    
  end subroutine invert_table

  
  subroutine edge_to_node_convert_table(N, x_to_e, Ne, e_to_y, x_to_y)
    integer, intent(in) :: N
    type(TEdgeList), intent(in) :: x_to_e(N)
    integer, intent(in) :: Ne
    integer, intent(in) :: e_to_y(Ne)
    type(TEdgeList), intent(inout) :: x_to_y(N)

    integer :: i, j

    do i = 1, N
       x_to_y(i) = TEdgeList(x_to_e(i)%N)
       do j = 1, x_to_y(i)%N
          x_to_y(i)%data(j) = 1 + e_to_y(x_to_e(i)%data(j))
       end do
    end do
  end subroutine edge_to_node_convert_table

  
  subroutine destructor(decoder)
    type(TDecoder) :: decoder
    
    if (allocated(decoder%c_to_e)  ) deallocate(decoder%c_to_e)
    if (allocated(decoder%v_to_e)  ) deallocate(decoder%v_to_e)
    if (allocated(decoder%v_to_c)  ) deallocate(decoder%v_to_c)
    if (allocated(decoder%c_to_v)  ) deallocate(decoder%c_to_v)
    if (allocated(decoder%F_buffer)) deallocate(decoder%F_buffer)
    if (allocated(decoder%B_buffer)) deallocate(decoder%B_buffer)
  end subroutine destructor

  subroutine print_decoder(this)
    class(TDecoder) :: this
    integer :: i
    print *, "TDecoder"
    print *, "    CNum =", this%cnum
    print *, "    VNum =", this%vnum
    print *, "    ENum =", this%Ne
  end subroutine print_decoder


  ! --------------------------------
  ! --- Message passing routines ---
  ! --------------------------------
  
  subroutine process_xnode(buffer_x_to_y, Ne, m_in, m_out, f_plus_kind, total, B_buffer, F_buffer)
    ! Generic processing function, it works both for variable- and for check- nodes
    ! - for variable nodes, f_plus_kind is a simple addition routine
    ! - for check nodes, f_plus_kind can be:
    !     + the "box plus" for exact message passing
    !     + the min abs with sign product for the min-sum algorithm
    !     + whatever you want for your approximation
    type(TEdgeList), intent(in) :: buffer_x_to_y
    integer, intent(in)       :: Ne
    real(wp), intent(in)      :: m_in(Ne)
    real(wp), intent(out)     :: m_out(Ne)
    interface
       pure function f_real_real(x, y) result (z)
         import wp
         real(wp), intent(in) :: x
         real(wp), intent(in) :: y
         real(wp)             :: z
       end function f_real_real
    end interface
    procedure(f_real_real) :: f_plus_kind
    real(wp), intent(out), optional  :: total
    real(wp), intent(inout), target, optional  :: B_buffer(2:buffer_x_to_y%N)
    real(wp), intent(inout), target, optional  :: F_buffer(buffer_x_to_y%N-1)


    real(wp), pointer :: buffer_fwd(:)
    real(wp), pointer :: buffer_bwd(:)

    integer :: i, j, N

    N = buffer_x_to_y%N
    
    if (present(B_buffer)) then
       buffer_bwd => B_buffer
    else
       allocate(buffer_bwd(2:N))
    end if
    
    if (present(F_buffer)) then
       buffer_fwd => F_buffer
    else
       allocate(buffer_fwd(N-1))
    end if
    
    buffer_fwd(1) = m_in(buffer_x_to_y%data(1))
    j = 1
    do i=2, N - 1
       buffer_fwd(i) = f_plus_kind(buffer_fwd(j), m_in(buffer_x_to_y%data(i)))
       j = i
    end do

    buffer_bwd(N) = m_in(buffer_x_to_y%data(N))
    j = N
    do i = N - 1, 2, -1
       buffer_bwd(i) = f_plus_kind(buffer_bwd(j), m_in(buffer_x_to_y%data(i)))
       j = i
    end do

    m_out(buffer_x_to_y%data(1)) = buffer_bwd(2)
    do i = 2, buffer_x_to_y%N-1
       m_out(buffer_x_to_y%data(i)) = f_plus_kind(buffer_fwd(i-1), buffer_bwd(i+1))
    end do
    m_out(buffer_x_to_y%data(N)) = buffer_fwd(N-1)


    if (present(total)) then
       total = f_plus_kind(buffer_fwd(N-1), m_in(buffer_x_to_y%data(N)))
    end if
    if (.not. present(B_buffer)) then
       deallocate(buffer_bwd)
    end if
    if (.not. present(F_buffer)) then
       deallocate(buffer_fwd)
    end if
  end subroutine process_xnode

  pure function f_plus_add(x, y) result (z)
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y
    real(wp)             :: z

    z = x + y
  end function f_plus_add


  pure function f_plus_box(x, y) result (z)
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y
    real(wp)             :: z

    z = sign(min(abs(x), abs(y)), x*y) + &
         log(1 + exp(-abs(x+y))) - &
         log(1 + exp(-abs(x-y)))
  end function f_plus_box

  subroutine process_vnode(buffer_v_to_e, llr_channel, llr_updated, Ne, m_c_to_v, m_v_to_c, B_buffer, F_buffer)
    type(TEdgeList), intent(in) :: buffer_v_to_e
    real(wp), intent(in)  :: llr_channel
    real(wp), intent(out) :: llr_updated
    integer, intent(in)   :: Ne
    real(wp), intent(in)  :: m_c_to_v(Ne)
    real(wp), intent(out) :: m_v_to_c(Ne)
    real(wp), intent(inout), target, optional :: B_buffer(2:buffer_v_to_e%N)
    real(wp), intent(inout), target, optional :: F_buffer(buffer_v_to_e%N-1)

    real(wp) :: total
    integer :: i

    if (buffer_v_to_e%N > 1) then
       if (present(B_buffer) .and. present(F_buffer)) then
          call process_xnode(buffer_v_to_e, Ne, m_c_to_v, m_v_to_c, f_plus_add, total, B_buffer, F_buffer)
       else
          call process_xnode(buffer_v_to_e, Ne, m_c_to_v, m_v_to_c, f_plus_add, total)
       end if
    else
       total = 0.0d0
    end if
    
    do i = 1, buffer_v_to_e%N
       m_v_to_c(buffer_v_to_e%data(i)) = m_v_to_c(buffer_v_to_e%data(i)) + llr_channel
    end do

    llr_updated = total + llr_channel
  end subroutine process_vnode

  subroutine process_cnode(buffer_c_to_e, s, Ne, m_v_to_c, m_c_to_v, B_buffer, F_buffer)
    type(TEdgeList), intent(in) :: buffer_c_to_e
    logical, intent(in)   :: s
    integer, intent(in)   :: Ne
    real(wp), intent(in)  :: m_v_to_c(Ne)
    real(wp), intent(out) :: m_c_to_v(Ne)
    real(wp), intent(inout), target, optional :: B_buffer(2:buffer_c_to_e%N)
    real(wp), intent(inout), target, optional :: F_buffer(buffer_c_to_e%N-1)


    integer :: i
    real(wp) :: sign_factor

    if (s) then
       sign_factor = -1.0_wp
    else
       sign_factor = 1.0_wp
    end if
    
    if (buffer_c_to_e%N < 2) then
       ! Very unlikely code design: it means a specific bit is known through the syndrome
       m_c_to_v(buffer_c_to_e%data(1)) = sign_factor * 1e300_wp
       return
    end if

    if (present(B_buffer) .and. present(F_buffer)) then
       call process_xnode(buffer_c_to_e, Ne, m_v_to_c, m_c_to_v, f_plus_box, B_buffer=B_buffer, F_buffer=F_buffer)
    else
       call process_xnode(buffer_c_to_e, Ne, m_v_to_c, m_c_to_v, f_plus_box)
    end if
    
    do i = 1, buffer_c_to_e%N
       m_c_to_v(buffer_c_to_e%data(i)) = m_c_to_v(buffer_c_to_e%data(i)) * sign_factor
    end do
  end subroutine process_cnode


  logical function check_llr(this, llr, synd)
    class(TDecoder) :: this
    real(wp), intent(in) :: llr(this%vnum)
    logical, intent(in) :: synd(this%cnum)
    logical p

    integer :: i, j

    check_llr = .true.
    
    checknode_loop : do i = 1, this%cnum
       p = .false.
       do j = 1, this%c_to_v(i)%N
          p = p .xor. (llr(this%c_to_v(i)%data(j)) < 0.0d0)
       end do
       
       if (xor(p, synd(i))) then
          check_llr = .false.
          exit checknode_loop
       end if
    end do checknode_loop
  end function check_llr

  ! --------------------------------------------------
  ! --- Helpers for orther programs and simulators ---
  ! --------------------------------------------------
  
  function word_to_synd(this, word) result(synd)
    class(TDecoder) :: this
    logical, intent(in) :: word(this%vnum)
    logical :: synd(this%cnum)

    integer :: i, j

    do i = 1, this%cnum
       synd(i) = .false.
       do j = 1, this%c_to_v(i)%N
          synd(i) = xor(synd(i), word(this%c_to_v(i)%data(j)))
       end do
    end do
  end function word_to_synd


  function llr_to_word(N, llr) result(word)
    integer, intent(in) :: N
    real(wp), intent(in):: llr(N)
    logical :: word(N)

    word = (llr < 0.0d0)
  end function llr_to_word


  ! ------------------------
  ! --- Decoding routine ---
  ! ------------------------
  subroutine decode(this, llr_channel, llr_updated, synd, N_iterations)
    class(TDecoder)        :: this
    real(wp), intent(in)   :: llr_channel(this%vnum)
    real(wp), intent(out)  :: llr_updated(this%vnum)
    logical, intent(in)    :: synd(this%cnum)
    integer, intent(inout) :: N_iterations
    ! N_iterations is used to return the actual number of iterations

    real(wp) :: m_c_to_v(this%Ne), m_v_to_c(this%Ne)
    integer :: it, i
    
    if (this%check_llr(llr_channel, synd)) then
       llr_updated = llr_channel
       N_iterations = 0
       return
    end if

    m_c_to_v(:) = 0.0_wp
    do i = 1, this%vnum
       ! First round to propagate the channel LLRs to the checknodes inputs
       call process_vnode(this%v_to_e(i), llr_channel(i), llr_updated(i), this%Ne, m_c_to_v, m_v_to_c)
    end do

    decoding_loop: do it = 1, N_iterations
       do i= 1, this%cnum
          call process_cnode(this%c_to_e(i), synd(i), this%Ne, m_v_to_c, m_c_to_v, &
               this%B_buffer, this%F_buffer)
       end do
       
       do i = 1, this%vnum
          call process_vnode(this%v_to_e(i), llr_channel(i), llr_updated(i), this%Ne, m_c_to_v, m_v_to_c, &
               this%B_buffer, this%F_buffer)
       end do

       if (this%check_llr(llr_updated, synd)) then
          ! Early stop if the current llr already satisfies the syndrome
          N_iterations = it
          exit decoding_loop
       end if
    end do decoding_loop
  end subroutine decode
end module ldpc_decoder
