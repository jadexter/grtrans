

  program testgeo

  use geodesics
  use grtrans_inputs
  use class_four_vector

  implicit none

  type (geo) :: g,g2,g3
  type (geokerr_args) :: gargs
  integer :: i,gunit,status,j,l
  character(len=40) :: ifile='inputs.in'
  type (four_vector), dimension(:), allocatable :: kdifft
  real(kind=8) :: dummy, geodiff1, geodiff2, maxgeodiff, geodiff3

  gunit=12
  call read_inputs(ifile)

  i=1; j=1

  write(6,*) 'init args'
  call initialize_geokerr_args(gargs,nro*nphi)
  write(6,*) 'pixels'
!  write(6,*) 'args: ',use_geokerr,standard,mu0,j,mu0(j),spin,uout,uin, &
!       rcut,nrotype,a1,a2,b1,b2,nro,nphi,nup
  call initialize_pixels(gargs,use_geokerr,standard,mu0(j),spin,uout,uin, &
                  rcut,nrotype,a1,a2,b1,b2,nro,nphi,nup)
  write(6,*) 'geo',nup,nro,nphi
  call initialize_geodesic(g,gargs,i,status)
  write(6,*) 'status', status
!  allocate(kdifft(g%npts-1))

!  write(6,*) 'lambda', g%npts, g%gk%alpha(1), g%gk%beta(1), g%gk%a, g%gk%tpr, &
!       g%gk%mufill
!  do l=1,g%npts-1
!     write(6,*) cos(g%x(l)%data(3)), g%lambda(l+1)-g%lambda(l), g%tpmarr(l), g%tprarr(l)
!  enddo

! Read geokerr output:
  open(unit=gunit,file='geokerr_unit_test1.out')
  read(gunit,*) dummy
  write(6,*) 'dummy: ',dummy
!  read(gunit,*) dummy
  i=1
  call initialize_geodesic(g2,gargs,i,status,gunit)
  write(6,*) 'status: ',status
  close(unit=gunit)

  open(unit=gunit,file='geokerr_unit_test2.out')
  read(gunit,*) dummy
  write(6,*) 'dummy: ',dummy
!  read(gunit,*) dummy
  i=1
  call initialize_geodesic(g3,gargs,i,status,gunit)
  close(unit=gunit)

  geodiff1=maxval(abs(g%x%data(3)-g2%x%data(3)))
  geodiff2=maxval(abs(g%x%data(2)-g3%x%data(2)))
  
! Test difference between k^\mu, dx^\mu / d\lambda
 !     kdifft=(g%x(2:g%npts)-g%x(1:g%npts-1))
 !     kdifft%data(1)=kdifft%data(1)/(g%lambda(2:g%npts)-g%lambda(1:g%npts-1))
 !     kdifft%data(2)=kdifft%data(1)/(g%lambda(2:g%npts)-g%lambda(1:g%npts-1))
 !     kdifft%data(3)=kdifft%data(1)/(g%lambda(2:g%npts)-g%lambda(1:g%npts-1))
 !     kdifft%data(4)=kdifft%data(1)/(g%lambda(2:g%npts)-g%lambda(1:g%npts-1))
! Need to put this test in a different place -- can't do with npts = 1 which is the rest of the test.
 !     maxgeodiff=maxval((/(maxval(abs(g%k(2:g%npts)%data(1)-kdifft%data(1)))), &
 !          (maxval(abs(g%k(2:g%npts)%data(2)-kdifft%data(2)))),&
 !          (maxval(abs(g%k(2:g%npts)%data(3)-kdifft%data(3)))),&
 !          (maxval(abs(g%k(2:g%npts)%data(4)-kdifft%data(4))))/))

  open(unit=12,file='unit_test_geo.out')
  write(12,*) geodiff1, geodiff2, geodiff3, maxgeodiff
  close(unit=12)

  call del_geodesic(g)
  call del_geokerr_args(gargs)

!  deallocate(kdifft)
  end  program
