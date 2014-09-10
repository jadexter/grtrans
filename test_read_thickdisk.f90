program test_read_thickdisk

character(len=500) :: line
character(len=10) :: substr
integer :: val=11
open(unit=12,file='fieldline5206.bin',form='formatted')
read(12,'(A)') line
close(unit=12)
write(6,*) 'header: ',line
write(6,*) 'header: ',line(1:390)
write(6,*) 'header: ',line(1:392)
write(6,*) 'header: ',line(1:396)
write(6,*) 'header: ',line(1:398)
write(substr,fmt='(I4)') val
write(6,*) 'substr: ',substr,val

write(6,*) 'index: ',index(line,trim(adjustl(substr)),back=.true.)

end program test_read_thickdisk
