
       module pgrtrans

       implicit none

       real(8), allocatable, dimension(:,:,:) :: ivals!,ab
       real(8), allocatable, dimension(:,:) :: ab
       real(8), allocatable, dimension(:) :: freqs!,mdots,mu0

       contains

          subroutine grtrans_main(standard,mumin,mumax,nmu,phi0,spin,&
            uout,uin, rcut, nrotype, gridvals, nn,i1,i2,fname, dt, nt, nload, &
            nmdot, mdotmin, mdotmax,ename, mbh, nfreq, fmin, fmax, muval,&
            gmin, gmax,p1, p2, jetalpha, stype,use_geokerr, nvals, iname,&
            cflag, extra, debug,outfile,fdfile,fhfile,fgfile,fsim,fnt,findf,fnfiles,fjonfix, &
            fnw,fnfreq_tab,fnr,foffset,fdindf,fmagcrit,frspot,fr0spot,fn0spot,ftscl,frscl, &
            fwmin,fwmax,ffmin,ffmax,frmax,fsigt,ffcol,fmdot,fnscl,fnnthscl,fnnthp,fbeta, &
            fbl06,fnp,ftp,frin,frout,fthin,fthout,fphiin,fphiout,epcoefindx,epotherargs,nepotherargs)
           
             use omp_lib
             use grtrans, only: grtrans_driver
             use ray_trace, only: ray_set, initialize_raytrace_camera, &
                  kwrite_raytrace_camera, del_raytrace_camera
             use fluid_model, only: load_fluid_model, unload_fluid_model, &
                  advance_fluid_timestep, source_params, assign_source_params_type, &
                  fluid_args, assign_fluid_args
             use emissivity, only: emis_params
             use geodesics, only: initialize_pixels, geokerr_args, &
                  del_geokerr_args,initialize_geokerr_args, initialize_geo_tabs
             use chandra_tab24, only: load_chandra_tab24, del_chandra_tab24

!INPUTS=====================
            integer, intent(in) :: standard,nrotype,nvals,nfreq,nmu,cflag,nt, &
                 nmdot,nload,extra,i1,i2,nepotherargs,debug
            integer :: nro,nphi,nup,tempi
            logical, intent(in) :: use_geokerr
            real(kind=8), intent(in) :: mumax,mumin,rcut,mbh,uout,uin, & 
                 fmin,fmax,dt,mdotmin,mdotmax,phi0,muval,gmin,gmax,p1,p2,jetalpha
            real(kind=8) :: a1,a2,b1,b2
            character(len=100), intent(in) :: ename,fname,iname,stype,outfile
            real(kind=8), dimension(:), allocatable :: mdots,mu0
!            real(kind=8), dimension(nfreq), intent(out) :: freqs
            real(kind=8), dimension(4),intent(in) :: gridvals
            real(kind=8) :: spin
            integer, dimension(3), intent(in) :: nn
! FLUID ARGUMENTS
            character(len=300), intent(in) :: fdfile,fhfile,fgfile,fsim
            integer, intent(in) :: fnt,findf,fnfiles,fjonfix,fnw,fnfreq_tab, &
                 fnr,foffset,fdindf,fmagcrit,fbl06
            real(8), intent(in) :: frspot,fr0spot,fn0spot,ftscl,frscl,fwmin,fwmax,ffmin, &
                 ffmax,frmax,fsigt,ffcol,fmdot,fnnthp,fnnthscl,fnscl,fbeta,ftp,fnp, &
                 frin,frout,fthin,fthout,fphiin,fphiout
            real(8), dimension(nepotherargs), intent(in) :: epotherargs
!            integer, intent(in) :: nepcoefindx
            integer, dimension(7), intent(in) :: epcoefindx
            !INPUTS====================
            integer :: nextra=0, inum, gunit, i, ncams, j, m, l, nparams, iii
            real(kind=8) :: wtime
            type (ray_set), dimension(:), allocatable :: c
            type (geokerr_args) :: gargs
            type (fluid_args) :: fargs
            type (source_params), dimension(:), allocatable :: sparams
            type (emis_params) :: eparams
            character(len=20), dimension(3) :: knames,kdescs
            !FITS output relevant information
            knames(1)='nx'; knames(2)='ny'; kdescs(1)='# x pixels'
            kdescs(2)='# y pixels'
            knames(3)='nu'; kdescs(3)='Frequency (Hz)'
            gunit=12
            !COMPUTE INPUTS==================================
            a1=dble(gridvals(1))
            a2=dble(gridvals(2))
            b1=dble(gridvals(3))
            b2=dble(gridvals(4))
            nro=nn(1)
            nphi=nn(2)
            nup=nn(3)
       !          write(6,*) 'read inputs nup: ',spin
            ! Compute frequencies w/ log spacing between and including fmin, fmax:
            allocate(freqs(nfreq)); allocate(mu0(nmu)); allocate(mdots(nmdot))
            if(nfreq==1) then
               freqs=(/fmin/)
            else
               freqs=fmin*exp((/(tempi-1,tempi=1,nfreq)/)*log(fmax/fmin)/(nfreq-1))
            endif
            if(nmdot==1) then
               mdots=(/mdotmin/)
            else
               mdots=mdotmin*exp((/(tempi-1,tempi=1,nmdot)/)*log(mdotmax/mdotmin)/(nmdot-1))
            endif
            if(nmu==1) then
               mu0=(/mumin/)
            else
               mu0=mumin+(mumax-mumin)/(nmu-1)*(/(tempi-1,tempi=1,nmu)/)
            endif
            !END COMPUTE INPUT AS IN READ_INPUTS================================
            
            !replace read_inputs with actual inputs
            if(extra==1) then
               if (nup.gt.1) then
                  nextra=19
               else
                  nextra=7
               endif
            endif
            ! these can later be added to a loop over emis parameter structures
            !       eparams%gmin=gmin; 
            eparams%gmax=gmax; eparams%p1=p1
            eparams%p2=p2;
            allocate(eparams%otherargs(nepotherargs))
            eparams%otherargs = epotherargs
            eparams%coefindx = epcoefindx
            !alwinnote 2015/04/05

            nparams=size(mdots)
            allocate(sparams(nparams))
            do iii=1,nparams
               sparams(iii)%mdot=mdots(iii); sparams(iii)%mbh=mbh
               sparams(iii)%nfac=ftscl; sparams(iii)%bfac=frscl
               sparams(iii)%jetalphaval=jetalpha
               sparams(iii)%gminval=gmin; sparams(iii)%muval=muval
               sparams(iii)%gmax=gmax
               sparams(iii)%p1=p1
               sparams(iii)%p2=p2
               call assign_source_params_type(sparams(iii),stype)
            enddo
            NCAMS=nparams*nfreq*nmu*nt
            allocate(c(NCAMS))
            if(outfile.ne."") write(6,*) 'outfile grtrans: ',outfile, NCAMS
            do m=1,NCAMS
               call initialize_raytrace_camera(c(m),nro,nphi,nvals,nextra)
            enddo
            
            call assign_fluid_args(fargs,fdfile,fhfile,fgfile,fsim,fnt,findf,fnfiles,fjonfix, &
            fnw,fnfreq_tab,fnr,foffset,fdindf,fmagcrit,frspot,fr0spot,fn0spot,ftscl,frscl, &
            fwmin,fwmax,ffmin,ffmax,frmax,fsigt,ffcol,fmdot,mbh,fnscl,fnnthscl,fnnthp,fbeta, &
            fbl06,fnp,ftp,frin,frout,fthin,fthout,fphiin,fphiout)
            call load_fluid_model(fname,spin,fargs)
            if(nup.eq.1.and.nvals.eq.4) call load_chandra_tab24()
            do j=1,nmu
               call initialize_geokerr_args(gargs,nro*nphi)
               call initialize_pixels(gargs,use_geokerr,standard,mu0(j),phi0,spin,uout,uin, &
                    rcut,nrotype,a1,a2,b1,b2,nro,nphi,nup)
               if(nload.gt.1) then
                  ! loop over geodesics calculating first point to get t0 for everything
                  gargs%nup=1
                  write(6,*) 'grtrans nload gt 1 loop',size(gargs%t0),allocated(gargs%t0),gargs%nup
                  !$omp parallel do private(i) shared(gargs)
                  do i=i1,i2
                     call initialize_geo_tabs(gargs,i)
                  enddo
                  !$omp end parallel do
                  write(6,*) 'grtrans after init geo tabs',minval(gargs%t0),maxval(gargs%t0)
                  where(gargs%t0.gt.0.)
                     gargs%t0=gargs%t0-minval(gargs%t0)
                  elsewhere
                     gargs%t0=0.
                  endwhere
                  gargs%nup=nup
                  write(6,*) 'grtrans after init geo tabs',minval(gargs%t0),maxval(gargs%t0)
               else
                  gargs%t0=0.
               endif
               wtime = omp_get_wtime()
               do l=1,nt
!$omp parallel do schedule(static,1) private(i) shared(gargs,gunit,c,j,nt,l,spin, &
!$omp& iname,ename,fname,sparams,eparams,nfreq,nparams,freqs,nup,i1,i2,extra,debug)
                  do i=i1,i2
!                  do i=15501,15999
                     call grtrans_driver(gargs,gunit,c,i,(j-1)*nt+l,iname,ename,fname, &
                     sparams,eparams,nfreq,nparams,freqs,nup,extra,debug)
                  enddo
!$omp end parallel do
                  if(l.lt.nt) call advance_fluid_timestep(fname,dt)
               enddo
               write(6,*) 'grtrans wall time elapsed: ', omp_get_wtime() - wtime
               spin=gargs%a
               call del_geokerr_args(gargs)
            enddo
            if(nup.eq.1.and.nvals.eq.4) call del_chandra_tab24()
            write(6,*) 'Write camera', nvals, nextra
            allocate(ivals(nvals+nextra,nro*nphi,NCAMS))
            ivals(:,:,:)=0d0
            allocate(ab(2,nro*nphi))
            ab=0d0
            ab(:,:) = dble(c(1)%pixloc)
            do m=1,ncams
               if(outfile.ne."") call kwrite_raytrace_camera(c(m),&                           
                    12,outfile,cflag,m,ncams,size(knames), &                                      
                    knames,kdescs,(/real(c(m)%nx),real(c(m)%ny),&                                 
                    real(FREQS(1+mod(m-1,nfreq)))/), &                                            
                    standard,mumin,mumax,nmu,phi0,spin,&                                          
                    uout,uin, rcut, nrotype, gridvals, nn, &                                      
                    fname, dt, nt, nload, nmdot, mdotmin, mdotmax, &                              
                    ename, mbh, nfreq, fmin, fmax, muval, gmin, gmax,&                            
                    p1, p2, jetalpha, stype, &                                                    
                    use_geokerr, nvals, iname, extra)
               ivals(:,:,m)=dble(c(m)%pixvals)
               call del_raytrace_camera(c(m))
            enddo
            call unload_fluid_model(fname)
            deallocate(c)
            deallocate(sparams); deallocate(eparams%otherargs)
            deallocate(mu0); deallocate(mdots)
            return
          end subroutine grtrans_main

          subroutine del_pgrtrans_data()
            deallocate(ivals); deallocate(ab); deallocate(freqs)
!            deallocate(mdots); deallocate(mu0)
          end subroutine del_pgrtrans_data
          
        end module pgrtrans
