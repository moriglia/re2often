program test_ldpc_edge_list
  use ldpc_edge_list, only: TEdgeList

  type(TEdgeList) :: edge_list


  edge_list = TEdgeList(5)

  edge_list%data(:) = 6

  print *, edge_list%N
  print *, edge_list%data

  edge_list = TEdgeList(9)
  edge_list%data(:) = 3

  print *, edge_list%N
  print *, edge_list%data
  print *, "Test TEdgeList completed"

end program test_ldpc_edge_list
