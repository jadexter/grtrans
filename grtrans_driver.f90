
       module grtrans

       use omp_lib

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
       use kerr, only: comoving_ortho, calc_polar_psi, calc_kb_ang, calcg, lnrf_frame, &
            comoving_ortho_debug
       use chandra_tab24, only: interp_chandra_tab24
       use radtrans_integrate, only: integrate, del_radtrans_integrate_data, &
            init_radtrans_integrate_data, intensity
       use math, only: tsum

       implicit none

! grtrans object is made up of fluid, geodesic, emissivity, and rad trans objects

!        integer, parameter :: IS_LINEAR_STOKES=1
       integer :: EXTRA_QUANTS=0, WRITE_GEO=0
!        integer, dimension(3) :: stats
!        real(kind=8), parameter :: MAX_TAU = 10d0
       real(kind=8) :: LBH, MBH
!        real(kind=8) :: ortol=1.d-6, oatol=1.d-8, hmax
!        real(kind=8), dimension(:), allocatable :: ss
       type (fluid) :: f
       type (rad_trans) :: r
       type (geo) :: g
       type (emis) :: e
!         type (emis_params) :: eparams
!         type (source_params) :: sparams
!$omp threadprivate(f,r,g,e,EXTRA_QUANTS,WRITE_GEO,LBH,MBH)

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
          iname,ename,fname,sparams,eparams,nfreq,nparams,freqs,nup,extra,debug)
! Driver routine for GRTrans. JAD 2/27/2011
         implicit none
         type (geokerr_args), intent(in) :: gargs
         type (ray_set), intent(inout), dimension(:) :: c
         type (emis_params), intent(in) :: eparams
         type (source_params), dimension(:), intent(in) :: sparams
         type (source_params) :: sp
         type (emis_params) :: ep
!         real(kind=8), intent(in) :: mbh1,mdot
         real(kind=8), intent(in), dimension(:) :: freqs
         integer, intent(in) :: gunit,i,l,nup, nfreq, nparams, extra, debug
         character(len=100), intent(in) :: iname,fname,ename
         real(kind=8) :: fac
         real(kind=8), dimension(:), allocatable :: s2xi,c2xi,s2psi,c2psi, &
          rshift,ang,nu,cosne,tau,tau_temp,intvals,dummy,vrl,vtl,vpl,cosne2, &
         aat,aar,aath,aaph,kht,khr,khth,khph,bht,bhr,bhth,bhph,dlp
         real(kind=8), dimension(:,:), allocatable :: aahat
         real(kind=8), dimension(:,:), allocatable :: tau_arr
         integer :: k,npow,status,j,nf,nitems,m,kk,taudex,ii,nvals
         integer, dimension(:), allocatable :: inds
         npow=3; j=0; status=0; k=0
         ep=eparams
         if(extra==1) EXTRA_QUANTS=1
! Need to make this a decision of load geodesics from file or generate new ones, and pass 
!gunit as last argument if loading.
!         write(6,*) 'grtrans_driver: ',i,debug
         call initialize_geodesic(g,gargs,i,status)
         if(debug==1) then
            WRITE_GEO=1
         else
            WRITE_GEO=0
         endif
         if (WRITE_GEO==1) then
            nitems=18
! Dump geodesic information for debugging
            open(unit=9,file='geodebug.out',form='formatted')
              ! header with number of points, number of items, mdot, mbh, frequency
            write(9,*) g%npts, nitems, g%gk%alpha(1), g%gk%beta(1),g%gk%su,g%gk%sm, &
                          sparams%mdot,sparams%mbh, fac, g%gk%a, g%gk%mu0
! geodesic properties
            write(9,*) g%lambda
            write(9,*) g%x%data(1)
            write(9,*) g%x%data(2)
            write(9,*) g%x%data(3)
            write(9,*) g%x%data(4)
            write(9,*) g%k%data(1)
            write(9,*) g%k%data(2)
            write(9,*) g%k%data(3)
            write(9,*) g%k%data(4)
            write(9,*) g%tprarr
            write(9,*) g%tpmarr
         endif
         if(status==1.and.g%npts.gt.0) then
!            hmax=(g%lambda(1)-g%lambda(2))/3.
!            if(hmax.gt.10.) hmax=10.
!            hmax=0.1
!            write(6,*) 'hmax: ',hmax
!         write(6,*) 'after init', g%npts,c(1)%nvals
            allocate(nu(g%npts))
            allocate(s2xi(g%npts)); allocate(c2xi(g%npts))
            s2xi=0d0; c2xi=0d0
            allocate(rshift(g%npts)); allocate(ang(g%npts))
            allocate(cosne(g%npts)); allocate(tau(g%npts))
            call initialize_fluid_model(f,fname,gargs%a,g%npts)
            call select_emissivity_values(e,ename)
!         write(6,*) 'evals: ',e%nk,e%neq
            ! Get fluid variables:
            call get_fluid_vars(g%x,g%k,g%gk%a,f)
            !write(6,*) 'fluid model: '!,f%rho
            call comoving_ortho(g%x%data(2),g%x%data(3),g%gk%a, &
                 g%gk%alpha(1),g%gk%beta(1),g%gk%mu0,f%u,f%b,g%k,s2xi,c2xi, &
                 ang,rshift,cosne)
            !write(6,*) 'comoving ortho grtrans driver: '!,s2xi,c2xi,'\n'
            call initialize_rad_trans(r,iname,g%npts,c(1)%nvals,extra)
            call init_radtrans_integrate_data(r%iflag,r%neq,g%npts,r%npts)
            do m=1,nparams
               sp=sparams(m)
               !write(6,*) 'grtrans_driver sp: ',sp%nfac,sp%gminval
               call initialize_source_params(sp,g%npts)
               call initialize_emis_params(ep,g%npts)
               call initialize_emissivity(e,g%npts,f%nfreq,rshift,ang,cosne,f%nrelbin,f%bingammamin,f%bingammamax)!,emisargs)
               !write(6,*) 'f%nrelbin',f%nrelbin
               !               write(6,*) 'm: ',m,nparams,sparams(m)%mdot,size(sparams)
!               allocate(sparams(m)%gmin(npts))
               MBH=sp%mbh
               LBH=GC*MBH*MSUN/C2
!               write(6,*) 'grtrans_driver convert: ',sp%nfac,sp%gminval
               call convert_model(sp)
!               write(6,*) 'gmin mu: ',sp%gmin
!               write(6,*) 'gmin mu: ',sp%mu
               ep%gmin=sp%gmin
               ep%mu=sp%mu
              !write(6,*) 'after convert'
!              write(6,*) 'emis model',FREQS
!              write(6,*) 'emis model'
               call emis_model(e,ep)
! Compute emissivity:
!               write(6,*) 'emis',NFREQ
               do k=1,NFREQ
                  nu=FREQS(k)/rshift
                  call calc_emissivity(nu,e,ep)
                  if(any(isnan(e%j))) then
                     write(6,*) 'NaN in emissivity at i: ',i
                  endif

!                  if(any(isnan(e%K))) then
                     !write(6,*) 'NaN in opacity at i: ',i
!                     write(6,*) 'NaN in emissivity where: ',maxloc(isnan(e%j))
!                     write(6,*) 'NaN in emissivity where: ',maxloc(isnan(e%K))
                     !write(6,*) 'NaN in emissivity fac LBH: ',fac,LBH
                     !write(6,*) 'fluid properties: ',e%ncgs,e%tcgs,e%bcgs
!                  endif
   
!               write(6,*) 'emis i: ',i
!               write(6,*) 'emis j: ',minval(e%j(:,1))
!!               write(6,*) 'emis n: ',minval(e%ncgsnth)
!               write(6,*) 'emis b: ',minval(e%bcgs)
!               write(6,*) 'emis ang: ',e%incang
! don't integrate if there's no emission:
!                  write(6,*) 'ej: ',sum(e%j(:,1))
!               write(6,*) 'ortho: ',s2xi
!                  write(6,*) 'tiny: ',epsilon(sum(e%j(:,1))),g%npts*epsilon(sum(e%j(:,1)))
!                  if (sum(e%j(:,1)).gt.(g%npts*tiny(sum(e%j(:,1))))) then
!                  fac=sum(e%j(:,1))/g%npts
                  if (g%npts.gt.1) then
                     fac=sum(e%j(:,1)*(g%lambda(1)-g%lambda(1:g%npts)))
                  else
                     fac=sum(e%j(:,1)/g%npts)
                  endif
                  !write(6,*) 'i fac: ',i,fac,1d0/(1d0/fac)
!                  if(g%npts.gt.1) then
!                     fac=tsum(g%lambda(1)-g%lambda(1:g%npts),e%j(:,1))
!                  else
!                     fac=1d0
!                  endif
!                  if (1d0/(1d0/(fac/1d10)).gt.0d0) then
                  if (fac/1d16.gt.0d0) then
!                     write(6,*) 'emis i: ',i
!                     write(6,*) 'emis j: ',maxval(e%j(:,1))
!                     write(6,*) 'emis n: ',minval(e%ncgs)
!                     write(6,*) 'emis b: ',minval(e%bcgs)
!                     write(6,*) 'emis te: ',i,e%tcgs(maxloc(e%ncgs))
                     if(e%neq.eq.4) call rotate_emis(e,s2xi,c2xi)
! Call integration routine:
                     !write(6,*) 'integrate'
!                     if(any(isnan(e%j))) then
                        !write(6,*) 'NaN in emissivity after rotate at i: ',i
                        !write(6,*) 'fluid properties: ',e%ncgs,e%tcgs,e%bcgs
!                     endif
!                     if(any(isnan(e%K))) then
                        !write(6,*) 'NaN in opacity after rotate at i: ',i
                        !write(6,*) 'fluid properties: ',e%ncgs,e%tcgs,e%bcgs
!                     endif

                     if(g%npts.ne.1) then
                        call invariant_emis(e,rshift)
!                        fac=sum(e%j(:,1))/g%npts
! scale emission close to 1 so that LSODA integrates accurately and put anu in cgs units
                        e%j=e%j/fac; e%K=e%K*LBH
!                        write(6,*) 'i: ',i
                        call calc_opt_depth(g%lambda(1:g%npts),e,tau,1)
                        tau = -tau
!                        call init_radtrans_integrate_data(r%iflag,r%neq,g%npts,r%npts)
                        !write(6,*) 'maxtau: ',maxval(e%K(:,1)),maxval(abs(sqrt(e%K(:,5)**2.+e%K(:,6)**2.))),maxval(abs(tau))
                        call integrate(g%lambda(1:g%npts),e%j,e%K,tau,r%npts)
                        r%I(1:r%neq,1:r%npts) = intensity(1:r%neq,1:r%npts)
!                        call del_radtrans_integrate_data()
!                        write(6,*) 'integrate: ',g%npts,r%neq,r%npts,r%I(1,r%npts)
! added factor of LBH for cgs intensity units ergs / cm^2 / s / Hz / ster
                        r%I=r%I*fac*LBH
!                        write(6,*) 'grtrans driver extra quants: ',EXTRA_QUANTS,r%npts
                        if (EXTRA_QUANTS==1) then
!                           if(r%npts.gt.1) then
                           allocate(tau_arr(g%npts,6))
                           allocate(tau_temp(g%npts)); allocate(intvals(g%npts))
                           allocate(dummy(g%npts)); allocate(inds(6)); allocate(dlp(r%npts-1))
                           inds=(/1,2,3,4,5,7/)
!                           write(6,*) 'tau_arr', size(g%lambda), size(tau_temp)
                           do ii=1,6
                              call calc_opt_depth(g%lambda(1)-g%lambda(1:g%npts),e,tau_temp,inds(ii))
!                              write(6,*) 'extra after opt_depth ',ii, &
!                                   size(tau_arr,1), size(tau_arr,2), size(tau_temp)
                              tau_arr(:,ii)=tau_temp
                           enddo
                           ! now calculate photosphere location, values there:
                           taudex=minloc(abs(tau-1d0),1)
                           if(minval(tau).lt.1d0) taudex=g%npts
                           if(taudex.eq.1) taudex=2
!                           write(6,*) 'taudex: ',taudex,tau(1),tau(g%npts),tau(taudex)
!                           write(6,*) 'ej: ',e%j(:,1), fac, r%I, tau
                           r%tau(1:6)=tau_arr(taudex,:)
!                           write(6,*) 'rtau',size(r%tau(1:6)),size(tau_arr(taudex,:))
                           ! now do emissivity-weighted ray averages of other quantities
                           dummy=e%j(:,1)*exp(-tau(1:g%npts))
                           if(r%neq==4) then
                              dlp=sqrt(r%I(1,r%npts:2:-1)**2.+r%I(2,r%npts:2:-1)**2.)-&
                                sqrt(r%I(1,r%npts-1:1:-1)**2.+r%I(2,r%npts-1:1:-1)**2.)
                           else
                              dlp(:)=0.
                           endif
!                           write(6,*) 'dlp: ',minval(dlp),maxval(dlp)
!                           write(6,*) 'rtau', size(r%tau),size(r%tau(7:7))
!                           write(6,*) 'tsum: ',tsum(g%lambda(1:taudex),dummy(1:taudex))
!                           write(6,*) 'tsum: ',tsum(g%lambda(1:taudex),dummy(1:taudex)*g%x(1:taudex)%data(2))
                           intvals=1d0/tsum(g%lambda,dummy)*tsum(g%lambda,dummy*g%x%data(2))
!                           write(6,*) 'intvals 1: ',intvals,tsum(g%lambda,dummy),tsum(g%lambda,dummy*g%x%data(2))
!                           write(6,*) 'g lambda: ',g%lambda(1:g%npts),g%npts
                           r%tau(7)=intvals(taudex)
!                           write(6,*) '7'
                           intvals=1d0/tsum(g%lambda,dummy)*tsum(g%lambda,dummy*g%x%data(3))
                           r%tau(8)=intvals(taudex)
                           intvals=1d0/tsum(g%lambda,dummy)*tsum(g%lambda,dummy*g%x%data(4))
                           r%tau(9)=intvals(taudex)
                           intvals=1d0/tsum(g%lambda,dummy)*tsum(g%lambda,dummy*e%ncgs)
                           r%tau(10)=intvals(taudex)
                           intvals=1d0/tsum(g%lambda,dummy)*tsum(g%lambda,dummy*e%tcgs)
                           r%tau(11)=intvals(taudex)
                           intvals=1d0/tsum(g%lambda,dummy)*tsum(g%lambda,dummy*e%bcgs)
                           r%tau(12)=intvals(taudex)
! this one is \beta...
                           intvals=1d0/tsum(g%lambda,dummy)*tsum(g%lambda,dummy*f%p*2./f%bmag**2.)
                           r%tau(13)=intvals(taudex)
! new ones added for fraction of emission produced above/below midplane (e.g. counter-jet or forward jet images), and weighted r,th,\tau_FR,\tau_FC for LP
                           intvals=1d0/tsum(g%lambda,dummy)*tsum(g%lambda,dummy*sign(1d0,cos(g%x%data(3))))
                           r%tau(14)=intvals(taudex)
                           r%tau(15)=sum(dlp*g%x(1:r%npts-1)%data(2))/sum(dlp)
                           r%tau(16)=sum(dlp*g%x(1:r%npts-1)%data(3))/sum(dlp)
                           r%tau(17)=sum(dlp*tau_arr(1:r%npts-1,5))/sum(dlp)
                           r%tau(18)=sum(dlp*tau_arr(1:r%npts-1,6))/sum(dlp)
                           r%tau(19)=sum(dlp*sign(1d0,cos(g%x(1:r%npts-1)%data(3))))/sum(dlp)
                           deallocate(tau_arr); deallocate(tau_temp); deallocate(intvals); deallocate(dummy)
                           deallocate(inds); deallocate(dlp)
!                           else
                        endif
                     else
!             write(6,*) 'rshift :', size(rshift)
                        call invariant_emis(e,rshift,npow)
!             write(6,*) 'compute intensity'
                        call grtrans_compute_intensity()
                        if(EXTRA_QUANTS==1) then
!                           write(6,*) 'grtrans driver extra quants npts 1'
                           allocate(s2psi(g%npts)); s2psi=0d0
                           allocate(c2psi(g%npts)); c2psi=0d0
                           allocate(cosne2(g%npts)); cosne2=0d0
                           call transpol(rshift,s2psi,c2psi,cosne2)
                           r%tau(1)=s2xi(1)
                           r%tau(2)=s2psi(1)
                           r%tau(3)=cosne2(1)
                           r%tau(4)=ang(1)
                           r%tau(5)=cosne(1)
                           r%tau(6)=c2xi(1)
                           r%tau(7)=c2psi(1)
                           deallocate(s2psi); deallocate(cosne2); deallocate(c2psi)
                        endif
                     endif
                  else
                     r%I(:,:)=0d0; r%npts=1
                     if(EXTRA_QUANTS==1) r%tau(:)=0d0
!                     write(6,*) 'zero: ',status, size(r%I,1), size(r%I,2), r%I(:,r%npts)
                  endif
!           write(6,*) 'after intensity', i, r%npts, r%I(:,r%npts), sum(e%j(:,1))!, e%j(:,1), e%K(:,1)
!                  if(g%npts.gt.1.and.mod(i,500).eq.0.and.mod(k,nfreq).eq.0.and.mod(m,nparams).eq.0) &
!                       write(6,*) 'save ', i, r%I(:,r%npts),l,k,m
                  if(EXTRA_QUANTS==0) then
! this is normal situation
                     call save_raytrace_camera_pixel(c((l-1)*nfreq*nparams+(m-1)*nfreq+k), &
                       (/real(gargs%alpha(i)),real(gargs%beta(i))/),real(r%I(:,r%npts)),i)
                  else
! this is where we want to save additional info
!                     write(6,*) 'save: ',size(r%I(:,r%npts)),size(r%tau),&
!                          size((/real(r%I(:,r%npts)),real(r%tau)/))
!                     write(6,*) 'cam: ',size(c((l-1)*nfreq*nparams+(m-1)*nfreq+k)%pixvals)
!                     write(6,*) 'save cam extra: ',r%tau,s2xi(1),c2xi(1),status,g%npts,r%npts,EXTRA_QUANTS
!                     write(6,*) 'save cam extra: ',r%tau,real(r%tau)
!                     write(6,*) 'save cam extra: ',(/real(r%I(:,r%npts)),real(r%tau)/)
!                     write(6,*) 'save cam extra: ',sum(e%j(:,1))
                     call save_raytrace_camera_pixel(c((l-1)*nfreq*nparams+(m-1)*nfreq+k), &
                       (/real(gargs%alpha(i)),real(gargs%beta(i))/),(/real(r%I(:,r%npts)),&
                       real(r%tau)/),i)
!                     write(6,*) 'extra'
                  endif
                  !write(6,*) 'after save'
                  if (WRITE_GEO==1) then
! fluid properties
                     write(9,*) f%rho
                     write(9,*) f%p
                     write(9,*) f%bmag
                     write(9,*) f%u%data(1)
                     write(9,*) f%u%data(2)
                     write(9,*) f%u%data(3)
                     write(9,*) f%u%data(4)
                     write(9,*) f%b%data(1)
                     write(9,*) f%b%data(2)
                     write(9,*) f%b%data(3)
                     write(9,*) f%b%data(4)
! emission properties
                     write(9,*) e%ncgs
                     write(9,*) e%tcgs
                     write(9,*) e%bcgs
                     write(9,*) e%rshift
                     write(9,*) e%incang
                     write(9,*) e%j(:,1)
                     write(9,*) e%K(:,1)
! Faraday rotation, conversion coefs
                     write(9,*) e%j(:,2)
                     write(9,*) e%j(:,4)
                     write(9,*) e%K(:,2)
                     write(9,*) e%K(:,4)
                     write(9,*) e%K(:,5)
                     write(9,*) e%K(:,7)
! JAD 10/28/2019 removed /LBH/fac part here not sure ??
                     write(9,*) r%I(1,:)/fac/LBH
                     write(9,*) r%I(2,:)/fac/LBH
                     write(9,*) r%I(3,:)/fac/LBH
                     write(9,*) r%I(4,:)/fac/LBH
                     write(9,*) tau
                     write(9,*) e%j(:,3)
                     write(9,*) e%K(:,3)
                     write(9,*) e%K(:,6)
                     write(9,*) s2xi
                     write(9,*) c2xi
! calculate quantities with other methods:
                     write(9,*) calc_kb_ang(g%k,f%b,f%u,g%x%data(2),g%x%data(3),g%gk%a)
!                     call calc_polar_psi(g%x%data(2),cos(g%x%data(3)),g%gk%q2,g%gk%a, &
!                          g%gk%alpha,g%gk%beta,rshift,g%gk%mu0,g%k, &
!                          c2xi,s2xi,cosne)
!                     write(9,*) s2xi
!                     write(9,*) c2xi
                     allocate(vrl(g%npts)); allocate(vtl(g%npts)); allocate(vpl(g%npts))
                     call lnrf_frame(f%u%data(2)/f%u%data(1),f%u%data(3)/f%u%data(1), &
                          f%u%data(4)/f%u%data(1),g%x%data(2),g%gk%a,g%x%data(3), &
                          vrl,vtl,vpl)
                     write(9,*) calcg(1d0/g%x%data(2),cos(g%x%data(3)),g%gk%q2, &
                          g%gk%l,g%gk%a, &
                          g%tpmarr,g%tprarr,g%gk%su,g%gk%sm,vrl,vtl,vpl)
                     deallocate(vrl); deallocate(vtl); deallocate(vpl)
                     ! WRITE OTHER DEBUG QUANTITIES FROM COMOVING ORTHO
                     allocate(aahat(g%npts,3)); allocate(aat(g%npts))
                     allocate(bht(g%npts)); allocate(kht(g%npts))
                     allocate(bhr(g%npts)); allocate(khr(g%npts)); allocate(aar(g%npts))
                     allocate(bhth(g%npts)); allocate(khth(g%npts)); allocate(aath(g%npts))
                     allocate(bhph(g%npts)); allocate(khph(g%npts)); allocate(aaph(g%npts))
                     call comoving_ortho_debug(g%x%data(2),g%x%data(3),g%gk%a, &
                          g%gk%alpha(1),g%gk%beta(1),g%gk%mu0,f%u,f%b,g%k,s2xi,c2xi, &
                          ang,rshift,cosne,aahat,aat,aar,aath,aaph,kht,khr,khth,khph, &
                          bht,bhr,bhth,bhph)
                     write(9,*) s2xi
                     write(9,*) c2xi
                     write(9,*) aahat
                     write(9,*) aat
                     write(9,*) aar
                     write(9,*) aath
                     write(9,*) aaph
                     write(9,*) kht
                     write(9,*) khr
                     write(9,*) khth
                     write(9,*) khph
                     write(9,*) bht
                     write(9,*) bhr
                     write(9,*) bhth
                     write(9,*) bhph
                     write(9,*) fac
                     close(unit=9)
                     deallocate(aahat); deallocate(aat)
                     deallocate(bht); deallocate(kht)
                     deallocate(bhr); deallocate(khr); deallocate(aar)
                     deallocate(bhth); deallocate(khth); deallocate(aath)
                     deallocate(bhph); deallocate(khph); deallocate(aaph)
                  endif
!                  call del_rad_trans(r)
           !           write(6,*) 'after rad'
               enddo
!         write(6,*) 'before emis'
               call del_source_params(sp)
               call del_emis_params(ep)
               call del_emissivity(e)
        !         write(6,*) 'after emis'
!               deallocate(sparams(m)%gmin)
            enddo

            call del_rad_trans(r)
            call del_radtrans_integrate_data()
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
                       real(gargs%beta(i))/),(/(0.,j=1,c((l-1)*nfreq*nparams+(m-1)*nfreq+k)%nvals+c((l-1)* &
                       nfreq*nparams+(m-1)*nfreq+k)%nextra)/),i)
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
         real(kind=8), dimension(e%npts) :: ncgs,bcgs,tcgs,ncgsnth
         real(kind=8), dimension(e%npts,f%nrelbin) :: nnthcgs !AC f%nrelbin or e%nrelbin?
         real(kind=8), dimension(:), allocatable :: freqarr
         real(kind=8), dimension(:,:), allocatable :: fnu
         call convert_fluid_vars(f,ncgs,ncgsnth,nnthcgs,bcgs,tcgs,fnu,freqarr,sp)
!         write(6,*) 'CONVERTED'
!         write(6,*) 'sz: ',size(ncgs), size(bcgs), size(tcgs)
!         write(6,*) 'alloc: ',allocated(e%ncgs),allocated(e%bcgs),allocated(e%tcgs)
         call assign_emis_params(e,ncgs,ncgsnth,nnthcgs,bcgs,tcgs,fnu,freqarr,f%nfreq)
         end subroutine convert_model

         subroutine transpol(rshift,s2psi,c2psi,cosne)
! Subroutine to do CPS80 polarization for single points
! For now just uses Chandrasekhar table for degree, but could change that
! JAD 4/20/2012
         real(kind=8), dimension(:), intent(in) :: rshift
         real, dimension(g%npts) :: interpI,interpdel
         real(kind=8), dimension(g%npts) :: polarpsi
         real(kind=8), dimension(g%npts), intent(out) :: s2psi,cosne,c2psi
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
