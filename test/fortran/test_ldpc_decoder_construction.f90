program test_ldpc_decoder_construction
  use iso_c_binding, only: wp => c_double
  use ldpc_edge_list, only: TEdgeList
  use ldpc_decoder, only: TDecoder, process_xnode, process_vnode, process_cnode, llr_to_word
  implicit none
  integer :: i
  integer :: N
  integer, allocatable :: e_to_c(:), e_to_v(:)
  real(wp), allocatable :: m_c_to_v(:), m_v_to_c(:) 
  type(TEdgeList) :: buff
  type(TDecoder) :: decoder
  real(wp) :: llr_channel, llr_updated
  logical :: s

  real(wp), allocatable :: llr(:)
  logical, allocatable  :: synd(:), word(:)

  buff = TEdgeList(5)
  print *, buff%N
  buff%data(:) = [(i, i=1,5)]
  print *, buff%data(:)

  N = 5
  allocate(e_to_c(N))
  allocate(e_to_v(N))
  e_to_c(:) = [1, 0, 2, 1, 0]
  e_to_v(:) = [0, 1, 2, 3, 4]
  decoder = TDecoder(N, e_to_v, e_to_c)
  print *, "Decoder initialized"
  deallocate(e_to_c)
  deallocate(e_to_v)
  call decoder%print()

  N= 8
  allocate(e_to_c(N))
  allocate(e_to_v(N))
  e_to_c(:) = [0,3,4,2,2,1,3,0]
  e_to_v(:) = [0,1,2,3,4,5,6,5]
  print *, "New Decoder"
  decoder = TDecoder(8, e_to_v, e_to_c)
  print *, "New decoder done!"
  call decoder%print()
  deallocate(e_to_c)
  deallocate(e_to_v)


  N = 7
  allocate(e_to_c(N))
  allocate(e_to_v(N))
  e_to_v(:) = [0, 1, 2, 2, 3, 1, 4]
  e_to_c(:) = [0, 0, 0, 1, 1, 2, 2]
  decoder = TDecoder(N, e_to_v, e_to_c)
  call decoder%print()
  allocate(m_c_to_v(N))
  allocate(m_v_to_c(N))

  m_v_to_c(:) = [(0.1d0*i, i=1,7)]
  m_c_to_v(:) = 0
  call process_xnode(decoder%c_to_v(1), N, m_v_to_c, m_c_to_v, f_plus)

  print *, "Message list in : ", m_v_to_c(:)
  print *, "Message list out: ", m_c_to_v(:)

  deallocate(m_c_to_v)
  deallocate(m_v_to_c)
  deallocate(e_to_c)
  deallocate(e_to_v)

  
  
  N = 8
  allocate(e_to_c(N))
  allocate(e_to_v(N))
  e_to_v(:) = [0, 1, 2, 2, 3, 1, 4, 2]
  e_to_c(:) = [0, 0, 0, 1, 1, 2, 2, 2]
  decoder = TDecoder(N, e_to_v, e_to_c)
  call decoder%print()
  allocate(m_c_to_v(N))
  allocate(m_v_to_c(N))

  m_c_to_v    = [(0.1d0*i, i=1,N)]
  m_v_to_c(:) = 0
  llr_channel = 1.0d0
  call process_vnode(decoder%v_to_e(3), llr_channel, llr_updated, N, m_c_to_v, m_v_to_c)

  print *, "Message list C-V: ", m_c_to_v(:)
  print *, "Message list V-C: ", m_v_to_c(:)
  print *, llr_updated

  s = .true.
  ! just prepare some data for checknode
  call process_vnode(decoder%v_to_e(1), llr_channel, llr_updated, N, m_c_to_v, m_v_to_c)
  call process_vnode(decoder%v_to_e(2), llr_channel, llr_updated, N, m_c_to_v, m_v_to_c)
  ! test checknode
  call process_cnode(decoder%c_to_e(1), s, N, m_v_to_c, m_c_to_v)
  print *, "Message list V-C: ", m_v_to_c(:)
  print *, "Message list C-V: ", m_c_to_v(:)


  ! test check_llr
  allocate(llr(decoder%vnum))
  allocate(synd(decoder%cnum))

  llr = [-3.0, -1.0, 2.3, 0.8, 0.6]
  synd= [.false., .false., .true.]
  print*, "LLR checks to syndrome? ", decoder%check_llr(llr, synd)

  allocate(word(decoder%vnum))
  word = [.true., .false., .true., .false., .false.]
  synd = decoder%word_to_synd(word)

  print *, "word:",  word
  print *, "synd:",  synd
  print *, "llr=", llr, "      word=", llr_to_word(decoder%vnum, llr)

contains
  pure function f_plus(x,y) result(z)
    real(wp), intent(in) :: x
    real(wp), intent(in) :: y
    real(wp)             :: z

    z = x+y
  end function f_plus
end program test_ldpc_decoder_construction
