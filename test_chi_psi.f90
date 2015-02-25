
      program test_kerr

      use geodesics
      use class_four_Vector
      use kerr
      use grtrans_inputs
      use phys_constants, GC=>G
      use fluid_model

      implicit none

      real(kind=8), dimension(:), allocatable :: r, th, s2psi, c2psi, chi, cosne, pssave, &
           pcsave, cssave, ccsave, angsave, gsave, psi
      real(kind=8) :: q2,l,alpha,beta,su,sm
      real(kind=8), dimension(:), allocatable :: s2xi,c2xi,ang,rshift,vrl,vtl,vpl
      real(kind=8), dimension(:,:), allocatable :: aahat,bbhat, grshift, kbang
      type (four_Vector), dimension(:), allocatable :: uhat,bhat,khat,aa
      integer :: gunit, status, i
      type (geo) :: g
      type (geokerr_args) :: gargs
      type (fluid) :: f
      character(len=20) :: file='test_chi_psi.in'

      call read_inputs(file)

      write(6,*) 'inputs read'
      call load_fluid_model(fname,spin)
!      call initialize_pixels(gargs,.true.,standard,mu0,spin,uout,rcut,
!     & nrotype,a1,a2,b1,b2,nro,nphi,nup)
      write(6,*) 'fluid loaded', standard, mu0(1), spin, uout, uin, rcut, nrotype
      write(6,*) 'fluid: ',a1,a2,b1,b2,nro,nphi,nup
      call initialize_geokerr_args(gargs,nro*nphi)
!       write(6,*) 'init pixels: ',allocated(gargs%muf)
      call initialize_pixels(gargs,.true.,standard,mu0(1),spin,uout,uin, &
        rcut,nrotype,a1,a2,b1,b2,nro,nphi,nup)
      write(6,*) 'pixels'
      allocate(pssave(nro*nphi)); allocate(pcsave(nro*nphi))
      allocate(cssave(nro*nphi)); allocate(ccsave(nro*nphi))
      allocate(grshift(nup,nro*nphi)); allocate(kbang(nup,nro*nphi))
      allocate(gsave(nro*nphi)); allocate(angsave(nro*nphi))
      do i=1,nro*nphi
         call initialize_geodesic(g,gargs,gunit,i,status)
         !write(6,*) 'geo'
         call initialize_fluid_model(f,fname,spin,g%npts)
         !write(6,*) 'init fluid', size(g%x), g%gk%a
         call get_fluid_vars(g%x,g%k,g%gk%a,f)
         allocate(s2xi(g%npts)); allocate(c2xi(g%npts)); allocate(ang(g%npts))
         allocate(rshift(g%npts)); allocate(aa(g%npts)); allocate(aahat(g%npts,3))
         allocate(khat(g%npts)); allocate(bhat(g%npts))
         allocate(chi(g%npts)); allocate(psi(g%npts)); allocate(cosne(g%npts))
         allocate(vrl(g%npts)); allocate(s2psi(g%npts)); allocate(c2psi(g%npts))
         allocate(vtl(g%npts)); allocate(vpl(g%npts))
         ! assign f%b = fpar from Agol 1997 Ch 5
         f%b = calc_polvec(g%x%data(2),cos(g%x%data(3)),g%k,g%gk%a,asin(1d0))
!         kbang(:,i) = calc_kb_ang(g%k,f%b,f%u,g%x%data(2),g%x%data(3),g%gk%a)
         call comoving_ortho(g%x%data(2),g%x%data(3),g%gk%a,g%gk%alpha(1), &
              g%gk%beta(1),g%gk%mu0,f%u,f%b,g%k,s2xi,c2xi,ang,rshift,cosne)
         call calc_polar_psi(g%x%data(2),cos(g%x%data(3)), &
              g%gk%q2, g%gk%a, g%gk%alpha, g%gk%beta, rshift, g%gk%mu0, &
              g%k, c2psi, s2psi, cosne)
         ! other method for g
         call lnrf_frame(f%u%data(2)/f%u%data(1),f%u%data(3)/f%u%data(1),f%u%data(4)/f%u%data(1),g%x%data(2),g%gk%a, &
              g%x%data(3),vrl,vtl,vpl)
         grshift(:,i) = calcg(1d0/g%x%data(2),cos(g%x%data(3)),g%gk%q2,g%gk%l,g%gk%a,g%tpmarr,g%tprarr, &
              g%gk%su,g%gk%sm,vrl,vtl,vpl)
         ! other method for ang
         chi=atan(1.-c2xi,1.+c2xi)
         write(6,*) 'psi chi: ', c2psi, c2xi, s2psi, s2xi
!         write(6,*) 'ang: ',ang,kbang(:,i)
         pssave(i)=s2psi(1)
         pcsave(i)=c2psi(1)
         cssave(i)=s2xi(1)
         ccsave(i)=c2xi(1)
         angsave(i)=ang(1)
         gsave(i)=rshift(1)
         call del_geodesic(g); call del_fluid_model(f)
         deallocate(s2xi); deallocate(c2xi); deallocate(ang); deallocate(rshift)
         deallocate(aahat); deallocate(aa); deallocate(bhat); deallocate(khat)
         deallocate(chi); deallocate(psi); deallocate(cosne); deallocate(c2psi)
         deallocate(vrl); deallocate(vtl); deallocate(vpl); deallocate(s2psi)
      enddo
      open(unit=12,file='test_chi_psi.out')
      write(12,*) gargs%alpha
      write(12,*) gargs%beta
      write(12,*) pssave
      write(12,*) pcsave
      write(12,*) cssave
      write(12,*) ccsave
      write(12,*) angsave
!      write(12,*) kbang
      write(12,*) gsave
      write(12,*) grshift
      close(unit=12)
      deallocate(kbang); deallocate(grshift)
      deallocate(angsave); deallocate(gsave)
      deallocate(pssave); deallocate(pcsave); deallocate(cssave); deallocate(ccsave)
      end program test_kerr
