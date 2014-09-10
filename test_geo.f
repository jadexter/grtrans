
      include 'geodesics.f'
!      include 'fits.f'
      include 'camera.f'
      include 'read_inputs.f'
      program test_geo
      
      use ray_trace
      use geodesics
      use fits
      use grtrans_inputs
      implicit none
      type (geo) :: g
      type (ray_set) :: c
      type (geokerr_args) :: gargs
      integer :: unit, status, gunit, i
      integer, dimension(1) :: naxes
      integer :: tpm,tpr
!      double precision :: a1,a2,b1,b2,rcut,a
      double precision, dimension(:), allocatable :: alarr,
     &  bearr,q2arr,larr,suarr,smarr,ufarr,mufarr
      integer, dimension(:), allocatable :: tpmarr, tprarr
      double precision :: u0,uf,offset,uout
      double precision, dimension(:), allocatable :: 
     &                                tpri,tpmi,lambdai,
     &                                 dti,dphi,ufi,mufi
      integer :: usegeor,mufill,phit
      double precision, dimension(:), allocatable :: intensity
      character(len=20) :: outfile='test_geo.out'
 
      call read_inputs()
!      nvals=1
      write(6,*) 'nup: ',nup
      call initialize_raytrace_camera(c,nro,nphi,nvals)
      call initialize_pixels(gargs,.true.,standard,mu0,a,rcut,
     &  nrotype,a1,a2,b1,b2,nro,nphi,nup)
      allocate(intensity(nro*nphi))
      do i=1,nro*nphi
        call raytrace_camera_pixel_getnext(c)
        gargs%geonum=i
        call initialize_geodesic(g,gargs,gunit)
!        write(6,*) 'vals: ',real(g%lambda(g%npts))
        call save_raytrace_camera_pixel(c,(/real(gargs%alpha(i)),
     &              real(gargs%beta(i))/),(/real(g%lambda(g%npts))/))
!        intensity(i)=g%lambda(g%npts)
        call del_geodesic(g)
      end do
!      open(unit=12,file='test_geo.out')
!      write(12,*) gargs%alpha
!      write(12,*) gargs%beta
!      write(12,*) intensity
!      close(unit=12)
! Write FITS instead!
!      call create_fits(unit,outfile,status)
!      naxes(1)=nro*nphi
!      call write_fits_header(unit,16,1,naxes,status)
!      call write_fits_array(unit,1,1,nro*nphi,alarr,status)
!      call close_fits(unit,status)
!      call delete_fits(outfile,status)
 !     deallocate(intensity)
!      write(6,*) 'pixloc: ',c%pixloc(:,7),c%pixvals(:,7)
      write(6,*) 'c: ',c%nx, c%ny, c%nvals
      call write_raytrace_camera(c,12,outfile,0)
      call del_geokerr_args(gargs)
      call del_raytrace_camera(c)
      end program
