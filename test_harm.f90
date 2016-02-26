
  program testharm

  use class_four_vector
  use fluid_model_harm, only: initialize_harm_model, harm_vals, del_harm_data
  use kerr, only: calc_rms

  type (four_vector), dimension(:), allocatable :: x0,b,u
  real, dimension(:,:), allocatable :: r2d,z2d
  real, dimension(:), allocatable :: rho,p,r,z,rtemp,ztemp,bmag,ub,norm,zero
  real :: a,rmin,rmax,zmin,zmax
  integer :: nr, i, k
  integer, dimension(:), allocatable :: indx

  nr=128; nz=128
  a=0.9375; rmin=calc_rms(a); rmax=25.
  zmin=-10.; zmax=10.

  allocate(rtemp(nr)); allocate(ztemp(nz))
  allocate(z(nr*nz)); allocate(r(nr*nz)); allocate(rho(nr*nz)); allocate(p(nr*nz))
  allocate(r2d(nr,nz)); allocate(z2d(nr,nz)); allocate(bmag(nr*nz))
  allocate(norm(nr*nz)); allocate(ub(nr*nz)); allocate(zero(nr*nz))

  zero=0.

  allocate(x0(nr*nz)); allocate(b(nr*nz)); allocate(u(nr*nz))

!  ztemp=(/(zmin*exp(real(k-1)/(nz-1)*log(zmax/zmin)),k=1,nz)/)
!  rtemp=(/(rmin*exp(real(k-1)/(nr-1)*log(rmax/rmin)),k=1,nr)/)
  ztemp=(/(zmin+(zmax-zmin)*real(k-1)/(nz-1),k=1,nz)/)
  rtemp=(/(rmin+(rmax-rmin)*real(k-1)/(nr-1),k=1,nr)/)
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

  call initialize_harm_model(dble(a))
  
  write(6,*) 'after init'

! Set up x four-vector:
  x0%data(1)=-100.; x0%data(2)=sqrt(r*r+z*z)
  x0%data(3)=acos(z/sqrt(r*r+z*z)); x0%data(4)=0.

  call harm_vals(x0,a,rho,p,b,u,bmag)
! Make proper arrays of rho, v, b to check with BL09

!  write(6,*) 'after vals rho', rho
!  write(6,*) 'after vals b', b*b

! Output so I can plot in IDL

  open(unit=12,file='test_harm.out')
  write(12,*) nr*nz, a
  write(12,*) r, z
  write(12,*) rho, p
  write(12,*) b*b
  close(unit=12)
  norm=merge(real(abs(u*u+1.)),zero,r.gt.(1.+sqrt(1.-a*a)))
  ub=merge(real(abs(u*b)),zero,r.gt.(1.+sqrt(1.-a*a)))
  write(6,*) 'unit test harm: ',maxval(norm),maxval(abs(ub)),minval(rho)
  open(unit=12,file='unit_test_harm.out')
  write(12,*) minval(rho), maxval(abs(norm)), maxval(abs(ub))
  close(unit=12)

  call del_harm_data()

  deallocate(rtemp); deallocate(ztemp); deallocate(z); deallocate(r); deallocate(rho)
  deallocate(p); deallocate(b); deallocate(x0); deallocate(u)
  deallocate(ub); deallocate(norm); deallocate(zero)

  end program
