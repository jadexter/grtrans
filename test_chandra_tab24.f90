program test_chandra_tab24

use chandra_tab24

implicit none

real, dimension(5) :: testmu,testI,testdel
testmu=(/0.,1.,0.07,0.33,0.68/)

call load_chandra_tab24()
write(6,*) 'chandra vals: ',minval(ch_delta),maxval(ch_delta)
write(6,*) 'chandra vals: ',ch_I(1),ch_I(10)
write(6,*) 'chandra vals: ',ch_mu(20)

call interp_chandra_tab24(testmu,testI,testdel)
write(6,*) 'testmu: ',testI,testdel

call del_chandra_tab24

end program
