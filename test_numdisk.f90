

  program testnumdisk

  use fluid_model_numdisk, only: numdisk_vals, initialize_numdisk_model, del_numdisk_data
!  use fluid_model_thindisk, only: thindisk_vals, init_thindisk
 ! use emissivity, only: bbemis

  real, dimension(:,:), allocatable :: fnu, fnuthin
  real, dimension(:,:), allocatable :: phi2, r2d
  real, dimension(:), allocatable :: omega, T, omegathin, rthin, r, phi
  integer :: nr, nphi, i, k
  real :: a, rmax, rmin, phmin, phmax, rtemp, phtemp

  nr=100; a=0.; rmin=3.; rmax=1000.; phmin=0.; phmax=2.*acos(-1.)
  nphi=100

  write(6,*) 'allocating test num disk'

  allocate(r2d(nr,nphi)); allocate(phi2(nr,nphi))
  allocate(omega(nr*nphi))
  allocate(T(nr*nphi))
  allocate(r(nr*nphi))
  allocate(phi(nphi*nr))

  write(6,*) '2d phi, r arrays'

  do i=1,nphi
     r2d(:,i)=rmin*exp((/(k-1,k=1,nr)/)*log(rmax/rmin)/(nr-1))
  enddo
  do i=1,nr
     phi2(i,:)=phmin+(phmax-phmin)*((/(k-1,k=1,nphi)/)+0.5)/nphi
  enddo
  phi=reshape(phi2,(/nr*nphi/))
  r=reshape(r2d,(/nr*nphi/))

  write(6,*) 'init numdisk model'

  rtemp=71.88*2.
  phtemp=.06

  call initialize_numdisk_model(dble(a))
  
!  write(6,*) 'numdisk vals', r2d(:,1)
!  write(6,*) 'numdisk vals', phi2(1,:)

  call numdisk_vals(r,phi,a,T,omega)

!  write(6,*) 'omega: ',omega
!  write(6,*) 'T: ',T

  open(unit=12,file='test_numdisk.out')
  write(12,*) nr, nphi
  write(12,*) r, phi, T
  close(unit=12)

  deallocate(omega)
  deallocate(T)
  deallocate(r)

  call del_numdisk_data()

  end program
