  program test
    integer, dimension(10) :: test_arr
    test_arr=(/1.,2.,3.,4.,5.,6.,7.,8.,9.,10./)
    write(6,*) 'test_arr: ',test_arr
    write(6,*) 'shift test_arr 3: ',test_arr(3:2:1)
    write(6,*) 'shift test_arr 6: ',test_arr(6:5:1)
  end program test
