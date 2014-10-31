
       module grtrans

       use fluid_model, only: initialize_fluid_model, del_fluid_model, &
        get_fluid_vars, fluid, source_params, convert_fluid_vars, &
        initialize_source_params,del_source_params
       use class_rad_trans, only: initialize_rad_trans, del_rad_trans, &
        rad_trans 
       use geodesics, only: initialize_geodesic, del_geodesic, &
         geokerr_args, geo
       use emissivity
       use ray_trace, only: save_raytrace_camera_pixel, ray_set
!       use interpolate, only: get_weight, locate
       use phys_constants, only: GC => G, EC => E, CC => C, C2, MSUN
       use kerr, only: comoving_ortho, calc_polar_psi
       use chandra_tab24, only: interp_chandra_tab24
       use radtrans_integrate, only: integrate, del_radtrans_integrate_data, &
            init_radtrans_integrate_data, intensity

       implicit none

! grtrans object is made up of fluid, geodesic, emissivity, and rad trans objects

!        integer, parameter :: IS_LINEAR_STOKES=1
        integer, parameter :: WRITE_GEO=1
        integer :: EXTRA_QUANTS=0
!        integer, dimension(3) :: stats
!        double precision, parameter :: MAX_TAU = 10d0
        double precision :: LBH, MBH
!        double precision :: ortol=1.d-6, oatol=1.d-8, hmax
!        double precision, dimension(:), allocatable :: ss
         type (fluid) :: f
         type (rad_trans) :: r
         type (geo) :: g
         type (emis) :: e
!         type (emis_params) :: eparams
!         type (source_params) :: sparams
!$omp threadprivate(f,r,g,e)

!       interface grtrans_driver_jac_form
!         module procedure grtrans_driver_jac_form
!       end interface

!       interface grtrans_driver_rhs_form
!         module procedure grtrans_driver_rhs_form
!       end interface

       interface grtrans_driver
         module procedure grtrans_driver
       end interface

       contains

         subroutine grtrans_driver(gargs,gunit,c,i,l, &
          iname,ename,fname,sparams,eparams,nfreq,nparams,freqs,nup,extra)
! Driver routine for GRTrans. JAD 2/27/2011
         implicit none
         type (geokerr_args), intent(in) :: gargs
         type (ray_set), intent(inout), dimension(:) :: c
         type (emis_params), intent(in) :: eparams
         type (source_params), dimension(:), intent(in) :: sparams
         type (source_params) :: sp
         type (emis_params) :: ep
!         double precision, intent(in) :: mbh1,mdot
         double precision, intent(in), dimension(:) :: freqs
         integer, intent(in) :: gunit,i,l,nup, nfreq, nparams, extra
         character(len=20), intent(in) :: iname,fname,ename
         double precision :: fac
         double precision, dimension(:), allocatable :: s2xi,c2xi, &
          rshift,ang,nu,cosne,tau,tau_temp,intvals,dummy
         double precision, dimension(:,:), allocatable :: tau_arr
         integer :: k,npow,status,j,nf,nitems,m,kk,taudex,ii
         integer, dimension(:), allocatable :: inds
         npow=3; j=0; status=0; k=0
         ep=eparams
         if(extra==1) EXTRA_QUANTS=1
         ! these should all bes electron model inputs
!         eparams%gmin=100.; eparams%gmax=1.e5; eparams%p1=3.5
!         eparams%p2=3.5; eparams%mu=1./4.
!         write(6,*) 'mdot: ',mdot,i
! Need to make this a decision of load geodesics from file or generate new ones, and pass 
!gunit as last argument if loading.
!         write(6,*) 'grtrans_driver: ',i
         call initialize_geodesic(g,gargs,i,status)
!         write(6,*) 'geo: ',g%gk%alpha,g%gk%beta
!         write(6,*) 'geo r: ',g%x%data(2), g%x%data(3), g%x%data(4)
!         write(6,*) 'status: ',status
         if (WRITE_GEO==1) then
            nitems=18
! Dump geodesic information for debugging
            open(unit=8,file='geodebug.out',form='formatted')
              ! header with number of points, number of items, mdot, mbh, frequency
            write(8,*) g%npts, nitems, g%gk%alpha(1), g%gk%beta(1),g%gk%su,g%gk%sm, &
                          sparams%mdot,sparams%mbh, fac
! geodesic properties
            write(8,*) g%lambda
            write(8,*) g%x%data(1)
            write(8,*) g%x%data(2)
            write(8,*) g%x%data(3)
            write(8,*) g%x%data(4)
            write(8,*) g%k%data(1)
            write(8,*) g%k%data(2)
            write(8,*) g%k%data(3)
            write(8,*) g%k%data(4)
            write(8,*) g%tprarr
            write(8,*) g%tpmarr
         endif
         if(status==1.and.g%npts.gt.0) then
!            hmax=(g%lambda(1)-g%lambda(2))/3.
!            if(hmax.gt.10.) hmax=10.
!            hmax=0.1
!            write(6,*) 'hmax: ',hmax
!         write(6,*) 'after init', g%npts,c(1)%nvals
            allocate(nu(g%npts))
            allocate(s2xi(g%npts)); allocate(c2xi(g%npts))
            allocate(rshift(g%npts)); allocate(ang(g%npts))
            allocate(cosne(g%npts)); allocate(tau(g%npts))
            call initialize_fluid_model(f,fname,gargs%a,g%npts)
            call select_emissivity_values(e,ename)
!         write(6,*) 'evals: ',e%nk,e%neq
            ! Get fluid variables:
            call get_fluid_vars(g%x,g%k,g%gk%a,f)
!            write(6,*) 'fluid model: ',f%rho
            call comoving_ortho(g%x%data(2),g%x%data(3),g%gk%a, &
                 g%gk%alpha(1),g%gk%beta(1),g%gk%mu0,f%u,f%b,g%k,s2xi,c2xi, &
                 ang,rshift,cosne)
!            call init_radtrans_integrate_data(r%iflag,r%neq,g%npts,r%npts)
            do m=1,nparams
               sp=sparams(m)
               call initialize_source_params(sp,g%npts)
               call initialize_emis_params(ep,g%npts)
               call initialize_emissivity(e,g%npts,f%nfreq,rshift,ang,cosne)
!               write(6,*) 'm: ',m,nparams,sparams(m)%mdot,size(sparams)
!               allocate(sparams(m)%gmin(npts))
               MBH=sp%mbh
               LBH=GC*MBH*MSUN/C2
               call convert_model(sp)
!               write(6,*) 'gmin mu: ',sp%gmin
!               write(6,*) 'gmin mu: ',sp%mu
               ep%gmin=sp%gmin
               ep%mu=sp%mu
!              write(6,*) 'after convert'
!           write(6,*) 'emis model',FREQS
!               write(6,*) 'emis model'
               call emis_model(e,ep)
! Compute emissivity:
!               write(6,*) 'emis',NFREQ
               do k=1,NFREQ
                  nu=FREQS(k)/rshift
!                  write(6,*) 'init rad trans'
                  call initialize_rad_trans(r,iname,g%npts,c((l-1)*nfreq*nparams+(m-1)*nfreq+k)%nvals,extra)
!                  write(6,*) 'calc_emissivity',m,k,nfreq,nparams
!                  write(6,*) 'i: ',i
                  call calc_emissivity(nu,e)
!                  write(6,*) 'ej1: ',e%j(:,1)
!                  write(6,*) 'ej2: ',e%j(:,2)
!                  write(6,*) 'ej4: ',e%j(:,4)
                  if(any(isnan(e%j))) then
                     write(6,*) 'NaN in emissivity at i: ',i
                     write(6,*) 'fluid properties: ',e%ncgs,e%tcgs,e%bcgs
                  endif
!               write(6,*) 'emis j: ',e%j(:,1)
!               write(6,*) 'emis n: ',e%ncgsnth
!               write(6,*) 'emis b: ',e%bcgs
!               write(6,*) 'emis ang: ',e%incang
! don't integrate if there's no emission:
!                  write(6,*) 'ej: ',sum(e%j(:,1))
!               write(6,*) 'ortho: ',s2xi
                  if (sum(e%j(:,1)).ne.0d0) then
                     if(e%neq.eq.4) call rotate_emis(e,s2xi,c2xi)
! Call integration routine:
!                     write(6,*) 'integrate'
                     if(any(isnan(e%j))) then
                        write(6,*) 'NaN in emissivity after rotate at i: ',i
                        write(6,*) 'fluid properties: ',e%ncgs,e%tcgs,e%bcgs
                     endif
                     if(g%npts.ne.1) then
                        call invariant_emis(e,rshift)
!                     write(6,*) 'integrate: ',e%j(:,1)
!                     write(6,*) 'integrate: ',e%j(:,2)
!                     write(6,*) 'integrate: ',e%j(:,4)
!                     write(6,*) 'integrate: ',e%K(:,1)
!                     write(6,*) 'integrate: ',e%K(:,2)
!                     write(6,*) 'integrate: ',e%j(:,1)
                        fac=sum(e%j(:,1))
! scale emission close to 1 so that LSODA integrates accurately and put anu in cgs units
                        e%j=e%j/fac; e%K=e%K*LBH
!                        write(6,*) 'i: ',i
!                        call calc_opt_depth(g%lambda(1)-g%lambda(1:g%npts),e,tau,1)
!                        call calc_opt_depth(g%lambda(g%npts:1:-1),e,tau,1)
                        call calc_opt_depth(g%lambda(1:g%npts),e,tau,1)
                        tau = -tau
                        call init_radtrans_integrate_data(r%iflag,r%neq,g%npts,r%npts)
                        call integrate(g%lambda(1:g%npts),e%j,e%K,tau,r%npts)
! may be better to do init, del data steps separately and then can move to outside of some loops
                        r%I(1:r%neq,1:r%npts) = intensity(1:r%neq,1:r%npts);
                        call del_radtrans_integrate_data()
!                        write(6,*) 'integrate: ',g%npts,r%neq,r%npts,r%I(1,r%npts)
                        r%I=r%I*fac
                        !write(6,*) 'integrate'
                        if (EXTRA_QUANTS==1) then
                           allocate(tau_arr(g%npts,6))
                           allocate(tau_temp(g%npts)); allocate(intvals(g%npts))
                           allocate(dummy(g%npts)); allocate(inds(7))
                           inds=(/1,2,3,4,5,7/)
                           !            write(6,*) 'tau_arr'
                           do ii=1,6
                              call calc_opt_depth(g%lambda(1)-g%lambda(1:g%npts),e,tau_temp,inds(i))
                              tau_arr(:,i)=tau_temp
                           enddo
                           ! now calculate photosphere location, values there:
                           taudex=minloc(abs(tau-1.),1)
                           !            write(6,*) 'taudex: ',taudex,tau(g%npts)
                           r%tau(1:6)=tau_arr(taudex,:)
                           !            write(6,*) 'rtau',size(r%tau(1:6)),size(tau_arr(taudex,:))
                           ! now do emissivity-weighted ray averages of other quantities:
                           dummy=e%j(:,1)*exp(-tau(1:g%npts))
                           !            write(6,*) 'rtau', size(r%tau),size(r%tau(7:7))
                           !            write(6,*) 'tsum: ',tsum(g%lambda(1:taudex),dummy(1:taudex))
                           !            write(6,*) 'tsum: ',tsum(g%lambda(1:taudex),dummy(1:taudex)*g%x(1:taudex)%data(2))
                           intvals=1./tsum(g%lambda,dummy)*tsum(g%lambda,dummy*g%x%data(2))
                           r%tau(7)=intvals(taudex)
                           !            write(6,*) '7'
                           intvals=1./tsum(g%lambda,dummy)*tsum(g%lambda,dummy*g%x%data(3))
                           r%tau(8)=intvals(taudex)
                           intvals=1./tsum(g%lambda,dummy)*tsum(g%lambda,dummy*g%x%data(4))
                           r%tau(9)=intvals(taudex)
                           intvals=1./tsum(g%lambda,dummy)*tsum(g%lambda,dummy*e%ncgs)
                           r%tau(10)=intvals(taudex)
                           intvals=1./tsum(g%lambda,dummy)*tsum(g%lambda,dummy*e%tcgs)
                           r%tau(11)=intvals(taudex)
                           intvals=1./tsum(g%lambda,dummy)*tsum(g%lambda,dummy*e%bcgs)
                           r%tau(12)=intvals(taudex)
! this one is \beta...
                           intvals=1./tsum(g%lambda,dummy)*tsum(g%lambda,dummy*f%p*2./f%bmag**2.)
                           r%tau(13)=intvals(taudex)
                           deallocate(tau_arr); deallocate(tau_temp); deallocate(intvals); deallocate(dummy)
                           deallocate(inds)
                        endif
                     else
!             write(6,*) 'rshift :', size(rshift)
                        call invariant_emis(e,rshift,npow)
!             write(6,*) 'transpol'
!             if (e%neq==4) then
!                call transpol(rshift)
!                write(6,*) 'after emis', e%j(:,1),e%j(:,2),e%j(:,3)
!                call subtest()
!             endif
!             write(6,*) 'compute intensity'
                        call grtrans_compute_intensity()
                     endif
                  else
                     r%I=0d0; r%npts=1
!                     write(6,*) 'zero: ',status, size(r%I,1), size(r%I,2), r%I(:,r%npts)
                  endif
!           write(6,*) 'after intensity', i, r%npts, r%I(:,r%npts), sum(e%j(:,1))!, e%j(:,1), e%K(:,1)
                  if(mod(i,500).eq.0.and.mod(k,nfreq).eq.0.and.mod(m,nparams).eq.0) write(6,*) 'save ', i, r%I(:,r%npts),l,k,m
                  if(EXTRA_QUANTS==0) then
! this is normal situation
                     call save_raytrace_camera_pixel(c((l-1)*nfreq*nparams+(m-1)*nfreq+k), &
                       (/real(gargs%alpha(i)),real(gargs%beta(i))/),real(r%I(:,r%npts)),i)
                  else
! this is where we want to save additional info
!                     write(6,*) 'save: ',size(r%I(:,r%npts)),size(r%tau),&
!                          size((/real(r%I(:,r%npts)),real(r%tau)/))
!                     write(6,*) 'cam: ',size(c((l-1)*nfreq*nparams+(m-1)*nfreq+k)%pixvals)
                     call save_raytrace_camera_pixel(c((l-1)*nfreq*nparams+(m-1)*nfreq+k), &
                       (/real(gargs%alpha(i)),real(gargs%beta(i))/),(/real(r%I(:,r%npts)),real(r%tau)/),i)
!                     write(6,*) 'extra'
                  endif
!                  write(6,*) 'after save'
                  if (WRITE_GEO==1) then
! fluid properties
                     write(8,*) f%rho
                     write(8,*) f%p
                     write(8,*) f%bmag
                     write(8,*) f%u%data(1)
                     write(8,*) f%u%data(2)
                     write(8,*) f%u%data(3)
                     write(8,*) f%u%data(4)
                     write(8,*) f%b%data(1)
                     write(8,*) f%b%data(2)
                     write(8,*) f%b%data(3)
                     write(8,*) f%b%data(4)
! emission properties
                     write(8,*) e%ncgs
                     write(8,*) e%tcgs
                     write(8,*) e%bcgs
                     write(8,*) e%rshift
                     write(8,*) e%incang
                     write(8,*) e%j(:,1)
                     write(8,*) e%K(:,1)
! Faraday rotation, conversion coefs
                     write(8,*) e%j(:,2)
                     write(8,*) e%j(:,4)
                     write(8,*) e%K(:,2)
                     write(8,*) e%K(:,4)
                     write(8,*) e%K(:,5)
                     write(8,*) e%K(:,7)
                     write(8,*) r%I(1,:)/fac
                     write(8,*) r%I(2,:)/fac
                     write(8,*) r%I(3,:)/fac
                     write(8,*) r%I(4,:)/fac
                     write(8,*) tau
                     write(8,*) e%j(:,3)
                     write(8,*) e%K(:,3)
                     write(8,*) e%K(:,6)
                     close(unit=8)
                  endif
                  call del_rad_trans(r)
           !           write(6,*) 'after rad'
               enddo
!         write(6,*) 'before emis'
               call del_source_params(sp)
               call del_emis_params(ep)
               call del_emissivity(e)
        !         write(6,*) 'after emis'
!               deallocate(sparams(m)%gmin)
            enddo
!            call del_radtrans_integrate_data()
            call del_fluid_model(f)
            
     !         write(6,*) 'after fluid'
            deallocate(nu); deallocate(s2xi); deallocate(c2xi)
            deallocate(rshift); deallocate(ang); deallocate(cosne); deallocate(tau)
!         write(6,*) 'after xi'
!         if(g%npts==1.and.e%neq==4) then
!            deallocate(jout)
!            deallocate(jtemp)
!         endif
!         write(6,*) 'after polarpsi'
         else
! Geodesic failure. Save intensity as zero
!         write(6,*) 'status m1'
            do m=1,NPARAMS
               do k=1,NFREQ
                  call save_raytrace_camera_pixel(c((l-1)*nfreq*nparams+(m-1)*nfreq+k),(/real(gargs%alpha(i)), &
                       real(gargs%beta(i))/),(/(0.,j=1,c((l-1)*nfreq*nparams+(m-1)*nfreq+k)%nvals)/),i)
               enddo
            enddo
         endif
!         write(6,*) 'del geo: ',allocated(g%gk%muf)
         call del_geodesic(g)
  !        write(6,*) 'del geo: ',allocated(g%gk%muf)
         end subroutine grtrans_driver

         subroutine convert_model(sp)
!         type (fluid), intent(in) :: f
!         type (emis), intent(inout) :: e
         type (source_params), intent(inout) :: sp
         double precision, dimension(e%npts) :: ncgs,bcgs,tcgs,ncgsnth
         double precision, dimension(:), allocatable :: freqarr
         double precision, dimension(:,:), allocatable :: fnu
!         write(6,*) 'convert', size(f%fnu)
         call convert_fluid_vars(f,ncgs,ncgsnth,bcgs,tcgs,fnu,freqarr,sp)
!         write(6,*) 'sz: ',size(ncgs), size(bcgs), size(tcgs)
!         write(6,*) 'alloc: ',allocated(e%ncgs),allocated(e%bcgs),allocated(e%tcgs)

         call assign_emis_params(e,ncgs,ncgsnth,bcgs,tcgs,fnu,freqarr,f%nfreq)
         end subroutine convert_model

         subroutine subtest()
! Routine to test memory allocation
           double precision, dimension(:), allocatable ::  subtest1,subtest2, &
                subtest3, subtest4, subtest5, subtest6
           allocate(subtest1(1)); allocate(subtest2(2)); allocate(subtest3(4))
           allocate(subtest4(8)); allocate(subtest5(16)); allocate(subtest6(32))
           deallocate(subtest1); deallocate(subtest2); deallocate(subtest3)
           deallocate(subtest4); deallocate(subtest5); deallocate(subtest6)
         end subroutine subtest

         subroutine transpol(rshift)
! Subroutine to do CPS80 polarization for single points
! For now just uses Chandrasekhar table for degree, but could change that
! JAD 4/20/2012
         double precision, dimension(:), intent(in) :: rshift
         real, dimension(g%npts) :: interpI,interpdel
         double precision, dimension(g%npts) :: polarpsi,s2psi,c2psi,cosne
!         write(6,*) 'rmu: ',g%x%data(2), g%x%data(3)
           call calc_polar_psi(g%x%data(2),cos(g%x%data(3)), &
                  g%gk%q2,g%gk%a, &
                  g%gk%alpha, &
                  g%gk%beta,rshift,g%gk%mu0,g%k, &
                  c2psi,s2psi,cosne)
!         write(6,*) 'polar psi: ',cosne,polarpsi
           call interp_chandra_tab24(real(cosne),interpI,interpdel)
           e%j(:,1)=e%j(:,1)*dble(interpI)
!           e%j(:,2)=e%j(:,1)*cos(polarpsi)*dble(interpdel)
!           e%j(:,3)=e%j(:,1)*sin(polarpsi)*dble(interpdel)
           e%j(:,2)=e%j(:,1)*c2psi*dble(interpdel)
           e%j(:,3)=e%j(:,1)*s2psi*dble(interpdel)
           e%j(:,4)=0.
         end subroutine transpol

         subroutine grtrans_compute_intensity()
! Subroutine to compute intensity at single point rather than integrate
! JAD 10/13/2011
!         write(6,*) 'grtrans intensity: ',e%neq,e%npts,r%neq,g%npts
!         write(6,*) 'grtrans intensity: ',size(r%I), size(e%j)
         r%I(1:e%neq,1:e%npts)=transpose(e%j(1:e%npts,1:e%neq))
         end subroutine grtrans_compute_intensity

       end module grtrans
