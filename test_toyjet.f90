

  program testffjet

  use class_four_vector
  use fluid_model_ffjet, only: initialize_ffjet_model, ffjet_vals, del_ffjet_data
  use kerr, only: calc_rms

  type (four_vector), dimension(:), allocatable :: x0,b,u
  real, dimension(:,:), allocatable :: r2d,z2d
  real, dimension(:), allocatable :: rho,p,r,z,rtemp,ztemp,bmag
  real :: a,rmin,rmax,zmin,zmax
  integer :: nr, i, k
  integer, dimension(:), allocatable :: indx

  nr=128; nz=128
  a=0.998; rmin=calc_rms(a); rmax=25.
  zmin=rmin; zmax=100.

  allocate(rtemp(nr)); allocate(ztemp(nz))
  allocate(z(nr*nz)); allocate(r(nr*nz)); allocate(rho(nr*nz)); allocate(p(nr*nz))
  allocate(r2d(nr,nz)); allocate(z2d(nr,nz)); allocate(bmag(nr*nz))

  allocate(x0(nr*nz)); allocate(b(nr*nz)); allocate(u(nr*nz))

  ztemp=(/(zmin*exp(real(k-1)/(nz-1)*log(zmax/zmin)),k=1,nz)/)
  rtemp=(/(rmin*exp(real(k-1)/(nr-1)*log(rmax/rmin)),k=1,nr)/)

  do i=1,nz 
     r2d(:,i)=rtemp
  enddo
  do i=1,nr
     z2d(i,:)=ztemp
  enddo

  r=reshape(r2d,(/nr*nz/))
  z=reshape(z2d,(/nr*nz/))

!  write(6,*) 'ztemp: ',zmin,zmax,ztemp
!  write(6,*) 'rtemp: ',rmin,rmax,rtemp
!  write(6,*) 'after: ',size(r),size(z),size(r2d),size(z2d)

  call initialize_ffjet_model(dble(a))
  
  write(6,*) 'after init'

! Set up x four-vector:
  x0%data(1)=0.; x0%data(2)=sqrt(r*r+z*z)
  x0%data(3)=acos(z/sqrt(r*r+z*z)); x0%data(4)=0.

  call ffjet_vals(x0,a,rho,p,b,u,bmag)
! Make proper arrays of rho, v, b to check with BL09

!  write(6,*) 'after vals rho', rho
!  write(6,*) 'after vals b', b*b

! Output so I can plot in IDL

  open(unit=12,file='test_ffjet.out')
  write(12,*) nr*nz, a
  write(12,*) r, z
  write(12,*) rho, p
  write(12,*) b*b
  close(unit=12)

  write(6,*) 'unit test ffjet: ',maxval(abs(u*u+1.)),maxval(u*b),minval(rho)
  open(unit=12,file='unit_test_ffjet.out')
  write(12,*) minval(rho), maxval(abs(u*u+1.)), maxval(abs(u*b))
  close(unit=12)

  call del_ffjet_data()

  deallocate(rtemp); deallocate(ztemp); deallocate(z); deallocate(r); deallocate(rho)
  deallocate(p); deallocate(b); deallocate(x0); deallocate(u); deallocate(bmag)

  end program
