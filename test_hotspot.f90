  
  program testhotspot

  use class_four_vector
  use fluid_model_hotspot, only: hotspot_vals, init_hotspot, advance_hotspot_timestep
  use emissivity, only: polsynchpl

  type (four_vector), dimension(:), allocatable :: u,b,x,xx
  real, dimension(:), allocatable :: ph, r, rtemp, phtemp, n, minn, maxnorm, maxub, norm
  real, dimension(:,:), allocatable :: ph2d, r2d
  integer :: nr, i, nph, nt
  real :: a, rmax, rmin, phmax, phmin, dt

  nr=10; nph=50; a=0.998; rmax=10.; phmin=0.; phmax=2.*acos(-1.)
  rmin=2; rmax=12.
  allocate(x(nr*nph)); allocate(b(nr*nph)); allocate(u(nr*nph))
  allocate(r(nr*nph)); allocate(n(nr*nph)); allocate(xx(nr*nph))
  allocate(r2d(nr,nph)); allocate(ph2d(nr,nph))
  allocate(rtemp(nr)); allocate(phtemp(nph))
  allocate(ph(nr*nph)); allocate(norm(nr*nph))
  phtemp=(/(phmin+real(k-1)/(nph-1)*(phmax-phmin),k=1,nph)/)
  rtemp=(/(rmin+real(k-1)/(nr-1.)*(rmax-rmin),k=1,nr)/)


  dt=5.; nt=20
  allocate(minn(nt)); allocate(maxnorm(nt)); allocate(maxub(nt))
  write(6,*) 'test hotspot alloc'

  do i=1,nph
     r2d(:,i)=rtemp
  enddo
  do i=1,nr
     ph2d(i,:)=phtemp
  enddo
  r = reshape(r2d,(/nr*nph/))
  ph = reshape(ph2d,(/nr*nph/))

  open(unit=12,file='test_hotspot.out')
  write(12,*) nr, nph, nt
  write(12,*) r
  write(12,*) ph

!  write(6,*) 'r: ',r
!  write(6,*) 'phi: ',ph
  x%data(1)=0d0
  x%data(2)=r
  x%data(4)=ph
  x%data(3)=acos(0.)

  write(6,*) 'test hotspot init'

  call init_hotspot()
  write(6,*) 'test hotspot vals', size(x)
  do i=1,nt
     call hotspot_vals(x,a,n,b,u,xx)
     call advance_hotspot_timestep(dt)
     norm=u*u
     minn(i)=minval(n); maxnorm(i)=maxval(abs(norm+1.))
     maxub(i)=maxval(abs(u*b))
     write(12,*) n
     write(12,*) xx%data(2)
     write(12,*) xx%data(4)
  enddo
  close(unit=12)
  open(unit=12,file='unit_test_hotspot.out')
     write(12,*) minval(minn), maxval(maxnorm), maxval(maxub)
  close(unit=12)
  deallocate(r); deallocate(ph); deallocate(x); deallocate(b)
  deallocate(u); deallocate(n); deallocate(xx)
  deallocate(r2d); deallocate(ph2d)
  deallocate(phtemp); deallocate(rtemp)
  deallocate(norm); deallocate(minn); deallocate(maxnorm)
  deallocate(maxub)

  end program
