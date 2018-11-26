  program test_kerr

      use geodesics
      use class_four_Vector
      use kerr, only: calc_polar_psi, comoving_ortho_debug, comoving_ortho, calc_kappapw
      use grtrans_inputs
      use phys_constants, GC=>G
      use fluid_model

      implicit none

      real(kind=8), dimension(:), allocatable :: r, th
      real(kind=8) :: q2,l,alpha,beta,su,sm
      real(kind=8), dimension(:), allocatable :: s2xi,c2xi,c2psi,s2psi,ang,rshift,cosne, &
           aat,aar,aath,aaph,kht,khr,khth,khph,bhatt,bhatr,bhatth,bhatph, cosne2
      real(kind=8), dimension(:,:), allocatable :: aahat,bbhat
      type (four_Vector), dimension(:), allocatable :: uhat,bhat,khat,aa,kdifft
      real(kind=8) :: maxang, maxangdiff, maxgeodiff,maxkbdiff
      integer :: gunit, status
      type (geo) :: g
      type (geokerr_args) :: gargs
      type (fluid_args) :: fargs
      type (fluid) :: f
      character(len=100) :: file='inputs.in'

      call read_inputs(file)
      call assign_fluid_args(fargs,fdfile,fhfile,fgfile,fsim,fnt,findf,fnfiles,fjonfix, &
            fnw,fnfreq_tab,fnr,foffset,fdindf,fmagcrit,frspot,fr0spot,fn0spot,ftscl,frscl, &
            fwmin,fwmax,ffmin,ffmax,frmax,fsigt,ffcol,fmdot,mbh,fnscl,fnnthscl,fnnthp,fbeta,&
            fbl06,fnp,ftp,frin,frout,fthin,fthout,fphiin,fphiout,fscalefac)

      write(6,*) 'inputs read'
      call load_fluid_model(fname,spin,fargs)
!      call initialize_pixels(gargs,.true.,standard,mu0,spin,uout,rcut,
!     & nrotype,a1,a2,b1,b2,nro,nphi,nup)
      write(6,*) 'fluid loaded', standard, mu0(1), spin, uout, uin, rcut, nrotype
      write(6,*) 'fluid: ',a1,a2,b1,b2,nro,nphi,nup
      call initialize_geokerr_args(gargs,nro*nphi)
!       write(6,*) 'init pixels: ',allocated(gargs%muf)
      call initialize_pixels(gargs,.true.,standard,mu0(1),phi0,spin,uout,uin, &
        rcut,nrotype,a1,a2,b1,b2,nro,nphi,nup)
      write(6,*) 'pixels'
      call initialize_geodesic(g,gargs,6258,status)
      write(6,*) 'geo'
      call initialize_fluid_model(f,fname,spin,g%npts)
      write(6,*) 'init fluid', size(g%x), g%gk%a
      call get_fluid_vars(g%x,g%k,g%gk%a,f)
      write(6,*) 'kdotk: ',maxval(g%k*g%k)
      write(6,*) g%gk%alpha(1), g%gk%beta(1)
!      r=g%x%data(2); th=g%x%data(3)
      allocate(s2xi(g%npts)); allocate(c2xi(g%npts)); allocate(ang(g%npts))
      allocate(rshift(g%npts)); allocate(aahat(g%npts,3))
      allocate(kdifft(g%npts-1)); allocate(aa(g%npts))
      allocate(s2psi(g%npts)); allocate(c2psi(g%npts)); allocate(cosne(g%npts))
      allocate(aat(g%npts)); allocate(aar(g%npts)); allocate(cosne2(g%npts))
      allocate(aath(g%npts)); allocate(aaph(g%npts)); allocate(bhatt(g%npts))
      allocate(bhatr(g%npts)); allocate(bhatth(g%npts)); allocate(bhatph(g%npts))
      allocate(kht(g%npts)); allocate(khr(g%npts)); allocate(khth(g%npts))
      allocate(khph(g%npts)); allocate(bhat(g%npts)); allocate(khat(g%npts))
      call comoving_ortho_debug(g%x%data(2),g%x%data(3),g%gk%a,g%gk%alpha(1), &
       g%gk%beta(1),g%gk%mu0,f%u,f%b,g%k,s2xi,c2xi,ang,rshift,cosne2,aahat,aat,aar, &
       aath,aaph,kht,khr,khth,khph,bhatt,bhatr,bhatth,bhatph)
      khat%data(1)=kht; khat%data(2)=khr; khat%data(3)=khth; khat%data(4)=khph
      bhat%data(1)=bhatt; bhat%data(2)=bhatr; bhat%data(3)=bhatth; bhat%data(4)=bhatph
      aa%data(1)=aat; aa%data(2)=aar; aa%data(3)=aath; aa%data(4)=aaph
      call calc_polar_psi(g%x%data(2),cos(g%x%data(3)),g%gk%q2,g%gk%a,g%gk%alpha, &
           g%gk%beta,rshift,g%gk%mu0,g%k,c2psi,s2psi,cosne)
! Test difference between psi, chi:
      maxang=maxval(abs(c2psi-c2xi))
      write(6,*) 'test kerr xi: ',c2xi
      write(6,*) 'test kerr psi: ',c2psi
      maxkbdiff=maxval(abs(cosne2-cosne))
! Test chi:
      maxangdiff=maxval(abs(c2xi*c2xi+s2xi*s2xi-1d0))
      write(6,*) 'diffs: ',maxang,maxangdiff,maxval(abs(g%k*g%k))
!      write(6,*) 'chi: ', c2xi
!      write(6,*) 'psi: ',c2psi
!      write(6,*) 'ang: ',cos(ang)
!      write(6,*) 'cosne: ',cosne
!      write(6,*) 'cosne2: ',g%k*f%b / (g%k*f%u)
!      write(6,*) 'cosne2: ',g%k(1)%data
!      write(6,*) 'cosne2: ',f%b(1)%data
!      write(6,*) 'cosne2: ',f%u(1)%data
!      write(6,*) 'cosne3: ',sqrt(g%gk%q2)*rshift/g%x%data(2), g%x%data(2), rshift, g%gk%q2
 ! write to file
! the definition in calc_kappapw is K1 - i K2 while in my transport_perpk it is K1 - i K2, so compare -K2 to imaginary part
      write(6,*) 'kappa: ',calc_kappapw(g%gk%a,g%x%data(2),cos(g%x%data(3)),g%k,aa), &
           g%gk%alpha+g%gk%a*sqrt(1.-g%gk%mu0**2.), -(-g%gk%beta)
      open(unit=12,file='unit_test_kerr.out')
      write(12,*) maxang, maxangdiff, maxval(abs(g%k*g%k)), &
           maxval(abs(khat%data(2)*aahat(:,1)+khat%data(3)*aahat(:,2) & 
           +khat%data(4)*aahat(:,3))), &
           maxval(abs(g%gk%beta-aimag(calc_kappapw(g%gk%a, &
           g%x%data(2),cos(g%x%data(3)),g%k,aa))))
      close(unit=12)
      open(unit=12,file='test_kerr.out')
      write(12,*) g%gk%alpha(1), g%gk%beta(1), g%gk%a, g%gk%mu0, g%gk%sm, g%gk%su, g%npts
      write(12,*) g%x%data(1)
      write(12,*) g%x%data(2)
      write(12,*) g%x%data(3)
      write(12,*) g%x%data(4)
      write(12,*) g%tprarr
      write(12,*) g%tpmarr
      write(12,*) g%k%data(1)
      write(12,*) g%k%data(2)
      write(12,*) g%k%data(3)
      write(12,*) g%k%data(4)
      write(12,*) f%u%data(1)
      write(12,*) f%u%data(2)
      write(12,*) f%u%data(3)
      write(12,*) f%u%data(4)
      write(12,*) f%b%data(1)
      write(12,*) f%b%data(2)
      write(12,*) f%b%data(3)
      write(12,*) f%b%data(4)
      write(12,*) ang
      write(12,*) rshift
      write(12,*) s2xi
      write(12,*) c2xi
      write(12,*) aa%data(2)
      write(12,*) aa%data(3)
      write(12,*) aa%data(4)
      write(12,*) aahat(:,1)
      write(12,*) aahat(:,2)
      write(12,*) aahat(:,3)
      write(12,*) bhat%data(1)
      write(12,*) bhat%data(2)
      write(12,*) bhat%data(3)
      write(12,*) bhat%data(4)
      write(12,*) khat%data(1)
      write(12,*) khat%data(2)
      write(12,*) khat%data(3)
      write(12,*) khat%data(4)
      close(unit=12)
      call del_geodesic(g); call del_fluid_model(f)
      deallocate(s2xi); deallocate(c2xi); deallocate(ang); deallocate(rshift)
      deallocate(aahat); deallocate(aa); deallocate(bhat); deallocate(khat)
      deallocate(c2psi); deallocate(s2psi); deallocate(kdifft); deallocate(cosne)
      deallocate(aat); deallocate(aar); deallocate(aath); deallocate(aaph)
      deallocate(bhatt); deallocate(bhatr); deallocate(bhatth); deallocate(bhatph)
      deallocate(kht); deallocate(khr); deallocate(khth); deallocate(khph)
      deallocate(cosne2)

      end program test_kerr
