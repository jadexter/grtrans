

  program testpol

  use kerr, only: calc_polar_psi
  use geodesics, only: initialize_geodesic, del_geodesic, geokerr_args, geo, &
       initialize_geokerr_args, initialize_pixels, del_geokerr_args

  implicit none

  type (geo) :: g
  type (geokerr_args) :: gargs
  integer :: gunit,i,status
  double precision, dimension(1) :: rshift=1., polarpsi,cosne
  double precision :: a1,a2,b1,b2,rcut,a,uout,uin,mu0
  logical :: use_geokerr
  integer :: nro,nphi,nup,nrotype,standard

  i=1
  standard=2
  a=0.99d0

  write(6,*) 'init args'
  call initialize_geokerr_args(gargs,1)
  write(6,*) 'pixels'
  call initialize_pixels(gargs,.true.,2,.25d0,a,1d-4,1d0, &
       1d0,2,-10d0,-10d0,5d0,5d0,1,1,1)
  write(6,*) 'geo'
  call initialize_geodesic(g,gargs,gunit,i,status)
!  call comoving_ortho(g%x%data(2),g%x%data(3),g%gk%a, &
!           g%gk%alpha(1),g%gk%beta(1),g%gk%mu0,f%u,f%b,g%k,s2xi,c2xi, &
!           ang,rshift)
  write(6,*) 'k: ',g%k(1)%data
  write(6,*) 'x: ',g%x(1)%data
  write(6,*) 'tpm: ',g%tpmarr,g%gk%sm(1)

  call calc_polar_psi(g%x%data(2),cos(g%x%data(3)),g%gk%q2,g%gk%a, &
       g%gk%alpha,g%gk%beta,rshift,g%gk%mu0,g%k, &
       polarpsi,cosne)

!  write(6,*) 'psi: ',polarpsi
  call del_geodesic(g)
  call del_geokerr_args(gargs)

  write(6,*) 'polarpsi: ',polarpsi
  write(6,*) 'cosne: ',cosne

  end  program
