      module fluid_model

      use class_four_vector
      use phys_constants, GC=>G
      use interpolate, only: interp, get_weight
      use kerr, only: kerr_metric, lnrf_frame, calc_rms, krolikc, calc_polvec, calc_u0, rms_vel
      use fluid_model_sphacc, only: sphacc_vals, init_sphacc, del_sphacc
      use fluid_model_sariaf, only: sariaf_vals, init_sariaf, del_sariaf
      use fluid_model_powerlaw, only: powerlaw_vals, init_powerlaw, del_powerlaw
      use fluid_model_ffjet, only: initialize_ffjet_model, del_ffjet_data, &
                                    ffjet_vals
      use fluid_model_phatdisk, only: phatdisk_vals, init_phatdisk, del_phatdisk, freq_tab
      use fluid_model_thindisk, only: thindisk_vals, init_thindisk
      use fluid_model_numdisk, only: initialize_numdisk_model, del_numdisk_data, &
           numdisk_vals
      use fluid_model_hotspot, only: hotspot_vals, init_hotspot, advance_hotspot_timestep
      use fluid_model_hotspot_schnittman, only: init_schnittman_hotspot, & 
           advance_schnittman_hotspot_timestep, schnittman_hotspot_vals
      use fluid_model_harm, only: initialize_harm_model, del_harm_data, harm_vals, &
           advance_harm_timestep
      use fluid_model_harm3d, only: initialize_harm3d_model, del_harm3d_data, harm3d_vals, &
          advance_harm3d_timestep
      use fluid_model_thickdisk, only: initialize_thickdisk_model, del_thickdisk_data, thickdisk_vals, &
           advance_thickdisk_timestep
      use fluid_model_mb09, only: initialize_mb09_model, del_mb09_data, mb09_vals, &
           advance_mb09_timestep
      use calc_gmin, only: calc_gmin_subroutine
      implicit none

      integer, parameter :: CONST=0,TAIL=1
      integer, parameter :: DUMMY=0,SPHACC=1,THINDISK=2,RIAF=3,HOTSPOT=4,PHATDISK=5,SCHNITTMAN=6
      integer, parameter :: COSMOS=10,MB=11,HARM=12,FFJET=13,NUMDISK=14,THICKDISK=15,MB09=16
      integer, parameter :: SARIAF=17,POWERLAW=18,HARM3D=19

      type fluid
        integer :: model, nfreq
        real :: rin
        real, dimension(:), allocatable :: rho,p,bmag,rho2
        real, dimension(:,:), allocatable :: fnu
        type (four_vector), dimension(:), allocatable :: u,b
      end type

      type fluid_args
         character(len=300) :: dfile,hfile,gfile,sim
         integer :: nt,indf,nfiles,jonfix,nw,nfreq_tab,nr,offset, &
              dindf,magcrit,bl06
         real(8) :: rspot,r0spot,n0spot,tscl,rscl,wmin,wmax,fmin, &
              fmax,rmax,sigt,fcol,mdot,mbh,nscl,nnthscl,nnthp,beta, &
              np,tp,rin,rout,thin,thout,phiin,phiout
      end type

      type source_params
!        real(kind=8), dimension(:), allocatable :: mdot,lleddeta,mu
        real(kind=8) :: nfac,bfac,mbh,mdot,p1,p2,gmax,gminval,jetalphaval,muval
        real(kind=8), dimension(:), allocatable :: gmin,jetalpha,mu
        integer :: type
      end type

!ultimately all source params stuff should probably go in own file / module
      interface initialize_source_params
         module procedure initialize_source_params
      end interface

      interface del_source_params
         module procedure del_source_params
      end interface
      
      interface assign_source_params
         module procedure assign_source_params
      end interface

      interface get_fluid_vars
        module procedure get_fluid_vars_single
        module procedure get_fluid_vars_arr
      end interface

      interface convert_fluid_vars
!       module procedure convert_fluid_vars_single
        module procedure convert_fluid_vars_arr
      end interface convert_fluid_vars

      interface initialize_fluid_model
        module procedure initialize_fluid_model
      end interface

      interface assign_fluid_args
         module procedure assign_fluid_args
      end interface
 
!      interface del_fluid_model
!        module procedure del_fluid_model
!      end interface

!      interface initialize_num_fluid_model
!        module procedure initialize_num_fluid_model
!      end interface

!      interface del_num_fluid_model
!        module procedure del_num_fluid_model
!      end interface

      interface advance_fluid_timestep
         module procedure advance_fluid_timestep
      end interface

      interface load_fluid_model
        module procedure load_fluid_model
      end interface

      interface unload_fluid_model
        module procedure unload_fluid_model
      end interface unload_fluid_model

      contains

        subroutine assign_fluid_args(fargs,dfile,hfile,gfile,sim,nt,indf,nfiles,jonfix, &
             nw,nfreq_tab,nr,offset,dindf,magcrit,rspot,r0spot,n0spot,tscl,rscl, &
             wmin,wmax,fmin,fmax,rmax,sigt,fcol,mdot,mbh,nscl,nnthscl,nnthp,beta,bl06,np,tp, &
             rin,rout,thin,thout,phiin,phiout)
          type (fluid_args), intent(inout) :: fargs
          character(len=100), intent(in) :: dfile,hfile,gfile,sim
          integer, intent(in) :: nt,indf,nfiles,jonfix,nw,nfreq_tab,nr,offset,dindf, &
               magcrit,bl06
          real(8), intent(in) :: rspot,r0spot,n0spot,tscl,rscl,wmin,wmax,fmin, &
               fmax,rmax,sigt,fcol,mdot,mbh,nscl,nnthscl,nnthp,beta,np,tp, &
               rin,rout,thin,thout,phiin,phiout
          fargs%dfile = dfile; fargs%hfile = hfile; fargs%gfile=gfile
          write(6,*) 'assign fluid args: ',fargs%dfile
          fargs%sim = sim; fargs%nt = nt; fargs%indf = indf; fargs%nfiles = nfiles
          fargs%jonfix = jonfix; fargs%nw = nw; fargs%nfreq_tab = nfreq_tab
          fargs%nr = nr; fargs%offset = offset; fargs%dindf = dindf
          fargs%magcrit = magcrit; fargs%rspot = rspot; fargs%r0spot = r0spot
          fargs%n0spot = n0spot; fargs%tscl = tscl; fargs%rscl = rscl
          fargs%wmin = wmin; fargs%wmax = wmax; fargs%fmin = fmin
          fargs%fmax = fmax; fargs%rmax = rmax; fargs%sigt = sigt
          fargs%mbh = mbh; fargs%fcol = fcol; fargs%mdot = mdot
          fargs%nscl = nscl; fargs%nnthscl = nnthscl; fargs%nnthp = nnthp
          fargs%beta = beta; fargs%bl06 = bl06; fargs%np = np; fargs%tp=tp
          fargs%rin = rin; fargs%rout = rout; fargs%thin = thin
          fargs%thout = thout; fargs%phiin = phiin; fargs%phiout = phiout
!          write(6,*) 'assign fluid args: ',jonfix,offset
        end subroutine assign_fluid_args

        subroutine load_fluid_model(fname,a,fargs)
        real(kind=8), intent(in) :: a
        character(len=20), intent(in) :: fname
        character(len=20) :: ifile
        type (fluid_args) :: fargs
        if(fname=='COSMOS') then
!          call initialize_cosmos_model(a,fargs)
        elseif(fname=='SARIAF') then
           call init_sariaf(real(fargs%nscl),real(fargs%tscl),real(fargs%nnthscl), &
                real(fargs%nnthp),real(fargs%beta),fargs%bl06) !alwinremark
!           call init_sariaf()
        elseif(fname=='POWERLAW') then
           call init_powerlaw(real(fargs%nscl),real(fargs%tscl),real(fargs%nnthscl), &
                real(fargs%nnthp),real(fargs%beta),real(fargs%np),real(fargs%tp), &
                real(fargs%rin),real(fargs%rout),real(fargs%thin),real(fargs%thout), &
                real(fargs%phiin),real(fargs%phiout)) !alwinremark
!           call init_powerlaw()
        elseif(fname=='MB') then
!          call intiialize_mb_model(a)
        elseif(fname=='THICKDISK') then
           call initialize_thickdisk_model(a,1,ifile,fargs%gfile,fargs%dfile,fargs%nt, &
                fargs%nfiles,fargs%indf,fargs%jonfix,fargs%offset,fargs%sim, &
                fargs%dindf,fargs%magcrit)
        elseif(fname=='MB09') then
           call initialize_mb09_model(a,1,ifile,fargs%gfile,fargs%dfile,fargs%nt, &
                fargs%nfiles,fargs%indf,fargs%jonfix,fargs%sim)
        elseif(fname=='HARM') then
          call initialize_harm_model(a,ifile,fargs%dfile,fargs%hfile,fargs%nt,fargs%indf)
        elseif(fname=='HARM3D') then
            call initialize_harm3d_model(a,ifile,fargs%dfile,fargs%hfile,fargs%gfile,fargs%nt,fargs%indf)
        elseif(fname=='SPHACC') then
           call init_sphacc()
        elseif(fname=='FFJET') then
    !      write(6,*) 'load'
          call initialize_ffjet_model(a,ifile,fargs%dfile)
        elseif(fname=='THINDISK') then
          call init_thindisk(real(a),ifile,real(fargs%mdot),real(fargs%mbh),real(fargs%rin),real(fargs%rout))
        elseif(fname=='PHATDISK') then
          call init_phatdisk(real(a),ifile,fargs%nw,real(fargs%wmin), &
               real(fargs%wmax),fargs%nfreq_tab,real(fargs%fmin), &
               real(fargs%fmax),fargs%nr, real(fargs%sigt), &
               real(fargs%fcol))
        elseif(fname=='NUMDISK') then
          call initialize_numdisk_model(ifile,fargs%dfile,real(fargs%tscl),&
               real(fargs%rscl))
        elseif(fname=='HOTSPOT') then
           call init_hotspot(ifile,real(fargs%rspot),real(fargs%r0spot), &
                real(fargs%n0spot))
        elseif(fname=='SCHNITTMAN') then
           call init_schnittman_hotspot(ifile,real(fargs%rspot),real(fargs%r0spot) &
                ,real(fargs%n0spot))
        endif
        end subroutine load_fluid_model

        subroutine advance_fluid_timestep(fname,dt)
        real(kind=8), intent(in) :: dt
        character(len=20), intent(in) :: fname
        if(fname=='COSMOS') then
!         call advance_cosmos_timestep(dt)
        elseif(fname=='MB') then
!           call advance_mb_timestep(dt)
        elseif(fname=='HARM') then
           call advance_harm_timestep(dt)
        elseif(fname=='HARM3D') then
           call advance_harm3d_timestep(dt)
        elseif(fname=='THICKDISK') then 
           call advance_thickdisk_timestep(dt)
        elseif(fname=='MB09') then 
           call advance_mb09_timestep(dt)
        elseif(fname=='HOTSPOT') then
           call advance_hotspot_timestep(real(dt))
        elseif(fname=='SCHNITTMAN') then
           call advance_schnittman_hotspot_timestep(real(dt))
        endif
        end subroutine advance_fluid_timestep

        subroutine initialize_fluid_model(f,fname,a,nup)
        character(len=20), intent(in) :: fname
        type (fluid), intent(out) :: f
        integer, intent(in) :: nup
        real(kind=8), intent(in) :: a
!        real(kind=8), intent(inout) :: uin
!        write(6,*) 'init: ',nup,a,fname
        if(fname=='PHATDISK') then
           f%model=PHATDISK; f%nfreq=size(freq_tab)
           allocate(f%u(nup)); allocate(f%fnu(nup,f%nfreq))
           allocate(f%b(nup))
        else
           allocate(f%u(nup))
! rho is used to store T for thindisk...
           allocate(f%rho(nup))
           allocate(f%b(nup))
           if(fname=='THINDISK') then
              f%model=THINDISK
           elseif(fname=='NUMDISK') then
              f%model=NUMDISK
           else
              allocate(f%bmag(nup))
              allocate(f%p(nup))
 ! !      write(6,*) 'sz: ',size(f%bmag)
              if(fname=='COSMOS') then
                 f%model=COSMOS
              elseif(fname=='MB') then
                 f%model=MB
              elseif(fname=='HARM') then
                 f%model=HARM
              elseif(fname=='HARM3D') then
                 f%model=HARM3D
              elseif(fname=='THICKDISK') then
                 f%model=THICKDISK
              elseif(fname=='MB09') then
                 f%model=MB09
              elseif(fname=='FFJET') then
                 f%model=FFJET
              elseif(fname=='SPHACC') then
                 f%model=SPHACC
              elseif(fname=='RIAF') then
                 f%model=RIAF
              elseif(fname=='HOTSPOT') then
                 f%model=HOTSPOT
              elseif(fname=='SCHNITTMAN') then
                 f%model=SCHNITTMAN
              elseif(fname=='SARIAF') then
                 allocate(f%rho2(nup))
                 f%model=SARIAF !alwinremark
              elseif(fname=='POWERLAW') then
                 allocate(f%rho2(nup))
                 f%model=POWERLAW
              else
                 write(6,*) 'WARNING: Unsupported fluid model -- using DUMMY'
                 f%model=DUMMY
              endif
           endif
        endif
!        f%model=model
!  !      write(6,*) 'fm: ',f%model, NUM_MODEL
!        if(f%model.gt.NUM_MODEL) call initialize_num_fluid_model(f,a)
        end subroutine initialize_fluid_model
 
        subroutine unload_fluid_model(fname)
        character(len=20), intent(in) :: fname
        if(fname=='COSMOS') then
!          call initialize_cosmos_model(a)
        elseif(fname=='MB') then
!          call intiialize_mb_model(a)
        elseif(fname=='HARM') then
           call del_harm_data()
        elseif(fname=='HARM3D') then
           call del_harm3d_data()
        elseif(fname=='THICKDISK') then
           call del_thickdisk_data()
        elseif(fname=='MB09') then
           call del_mb09_data()
        elseif(fname=='FFJET') then
    !      write(6,*) 'load'
           call del_ffjet_data()
        elseif(fname=='PHATDISK') then
           call del_phatdisk()
        elseif(fname=='NUMDISK') then
           call del_numdisk_data()
        elseif(fname=='SPHACC') then
           call del_sphacc()
        elseif(fname=='SARIAF') then
           call del_sariaf() !alwinremark does nothing so far
        elseif(fname=='POWERLAW') then
           call del_powerlaw()
        endif
        end subroutine unload_fluid_model

        subroutine del_fluid_model(f)
        type (fluid), intent(inout) :: f
        if(f%model==PHATDISK) then
           deallocate(f%u); deallocate(f%fnu)
           deallocate(f%b)
        else
           deallocate(f%u); deallocate(f%rho); deallocate(f%b)
!           write(6,*) 'del_fluid: ',f%model,THINDISK
           if(f%model.ne.THINDISK.and.f%model.ne.NUMDISK) then
              deallocate(f%p)
              deallocate(f%bmag)
           endif
           if(f%model==SARIAF.or.f%model==POWERLAW) deallocate(f%rho2)
        endif
        f%model=-1
        end subroutine del_fluid_model
       
        subroutine get_fluid_vars_single(x0,k0,a,f)
        type (four_vector), intent(in) :: x0, k0
        type (fluid), intent(inout) :: f
        real(kind=8), intent(in) :: a
 ! !      write(6,*) 'fluid: ',size(x0),f%model
        SELECT CASE(f%model)
          CASE (SPHACC)
!            call get_sphacc_fluidvars(x0,f)
          CASE (FFJET)
!      !      write(6,*) 'made it'
!            call get_ffjet_fluidvars(x0,real(a),f)
          CASE (THINDISK)
!            call initialize_thindisk_model(f)
          CASE (HOTSPOT) 
!            call initialize_hotspot_model(f)
          CASE (RIAF)
!            call initialize_ffjet_model(f)
          CASE (DUMMY)
        END SELECT
        end subroutine get_fluid_vars_single

        subroutine get_fluid_vars_arr(x0,k0,a,f)
        type (four_vector), intent(in), dimension(:) :: x0, k0
        type (fluid), intent(inout) :: f
        real(kind=8), intent(in) :: a
 ! !      write(6,*) 'fluid: ',size(x0),f%model
        SELECT CASE(f%model)
          CASE (SPHACC)
            call get_sphacc_fluidvars(x0,f)
          CASE (FFJET)
      !      write(6,*) 'made it'
            call get_ffjet_fluidvars(x0,real(a),f)
          CASE (THINDISK)
            call get_thindisk_fluidvars(x0,k0,real(a),f)
          CASE (PHATDISK)
            call get_phat_fluidvars(x0,k0,real(a),f)
          CASE (NUMDISK)
            call get_numdisk_fluidvars(x0,k0,real(a),f)
          CASE (HOTSPOT)
            call get_hotspot_fluidvars(x0,real(a),f)
          CASE (SCHNITTMAN)
            call get_schnittman_hotspot_fluidvars(x0,real(a),f)
          CASE (HARM)
             call get_harm_fluidvars(x0,real(a),f)
          CASE (HARM3D)
             call get_harm3d_fluidvars(x0,real(a),f)
          CASE (THICKDISK)
             call get_thickdisk_fluidvars(x0,real(a),f)
          CASE (MB09)
             call get_mb09_fluidvars(x0,real(a),f)
          CASE (SARIAF)
             call get_sariaf_fluidvars(x0,real(a),f) !alwinremark this exists below
          CASE (POWERLAW)
             call get_powerlaw_fluidvars(x0,real(a),f)
          CASE (DUMMY)
        END SELECT
        end subroutine get_fluid_vars_arr

        subroutine convert_fluid_vars_arr(f,ncgs,ncgsnth,bcgs,tcgs,fnuvals,freqvals,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(inout) :: sp
        real(kind=8), intent(out), dimension(size(f%rho)) :: ncgs,ncgsnth,bcgs,tcgs
        real(kind=8), intent(out), dimension(:), allocatable :: freqvals
        real(kind=8), intent(out), dimension(:,:), allocatable :: fnuvals
!        write(6,*) 'fluid convert: ',f%model,size(bcgs)
        SELECT CASE(f%model)
          CASE (SPHACC)
            call convert_fluidvars_sphacc(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (FFJET)
            call convert_fluidvars_ffjet(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (THINDISK)
            call convert_fluidvars_thindisk(f,tcgs,ncgs)
          CASE(PHATDISK)
             call convert_fluidvars_phatdisk(f,fnuvals,freqvals)
          CASE(NUMDISK)
             call convert_fluidvars_numdisk(f,tcgs,ncgs)
          CASE (HOTSPOT)
            call convert_fluidvars_hotspot(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (SCHNITTMAN)
            call convert_fluidvars_schnittman_hotspot(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (HARM)
             call convert_fluidvars_harm(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (HARM3D)
             call convert_fluidvars_harm3d(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (THICKDISK)
             call convert_fluidvars_thickdisk(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (MB09)
             call convert_fluidvars_mb09(f,ncgs,ncgsnth,bcgs,tcgs,sp)
          CASE (SARIAF)
             call convert_fluidvars_sariaf(f,ncgs,ncgsnth,bcgs,tcgs,sp) !alwinremark exists below
          CASE (POWERLAW)
             call convert_fluidvars_powerlaw(f,ncgs,ncgsnth,bcgs,tcgs,sp)
!          CASE (DUMMY)
        END SELECT
        call assign_source_params(sp,ncgs,tcgs,ncgsnth)

!        write(6,*) 'after convert'
        end subroutine convert_fluid_vars_arr

        subroutine get_thindisk_fluidvars(x0,k0,a,f)
        type (four_vector), intent(in), dimension(:) :: x0, k0
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
!        real :: rms,T0
        real, dimension(size(x0)) :: T,omega
        real, dimension(size(x0),10) :: metric
        call thindisk_vals(real(x0%data(2)),real(x0%data(3)),a,T,omega)
        f%rho=T
!        write(6,*) 'get fluidvars: ',size(f%rho), size(T)
        f%u%data(2)=0.; f%u%data(3)=0.
        f%b%data(1)=0.; f%b%data(2)=0.; f%b%data(3)=0.; f%b%data(4)=0.
!        write(6,*) 't: ',T,f%u%data(4),om,1./(r**3./2.+a),a*(r*r+a*a-d)/ &
        !  ((r*r+a*a)**2.-d*a*a*sin(x0%data(3))**2.)
        metric=kerr_metric(real(x0%data(2)),real(x0%data(3)),a)
!        f%b%data(2)=cos(x0%data(4))/metric(:,5)
!        f%b%data(4)=-sin(x0%data(4))/metric(:,10)
!        call assign_metric(f%b,transpose(metric))
        call assign_metric(f%u,transpose(metric))
        call assign_metric(f%b,transpose(metric))
        f%u%data(1)=sqrt(-1./(metric(:,1)+2.*metric(:,4)*omega+metric(:,10)* &
         omega*omega))
        f%u%data(4)=omega*f%u%data(1)
! Assign normal vector as magnetic field for comoving_ortho:
         f%b = calc_polvec(x0%data(2),cos(x0%data(3)),k0,dble(a),asin(1d0))
!         f%b%data(1)=-f%b%data(1)
!         f%b%data(2)=-f%b%data(2)
!         f%b%data(3)=-f%b%data(3)
!         f%b%data(4)=-f%b%data(4)
!        write(6,*) 'udotu: ',f%u*f%u,omega,f%u%data(1),r,rms
 !       write(6,*) 'thindisk fluidvars: ',allocated(f%u), &
!             allocated(f%b), allocated(f%rho)
        end subroutine get_thindisk_fluidvars

        subroutine get_numdisk_fluidvars(x0,k0,a,f)
        type (four_vector), intent(in), dimension(:) :: x0, k0
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
!        real :: rms,T0
        real, dimension(size(x0)) :: T,omega,phi
        real, dimension(size(x0),10) :: metric
        phi=x0%data(4); phi=phi+12.*acos(-1.)
        phi=mod(phi,(2.*acos(-1.)))
!        write(6,*) 'numdisk phi: ',minval(phi),maxval(phi)
        call numdisk_vals(real(x0%data(2)),phi,a,T,omega)
        f%rho=T
!        write(6,*) 'get fluidvars: ',size(f%rho), size(T)
        f%u%data(2)=0.; f%u%data(3)=0.
        f%b%data(1)=0.; f%b%data(2)=0.; f%b%data(3)=0.; f%b%data(4)=0.
!        write(6,*) 't: ',T,f%u%data(4),om,1./(r**3./2.+a),a*(r*r+a*a-d)/ &
        !  ((r*r+a*a)**2.-d*a*a*sin(x0%data(3))**2.)
        metric=kerr_metric(real(x0%data(2)),real(x0%data(3)),a)
        f%b%data(2)=cos(x0%data(4))/metric(:,5)
        f%b%data(4)=-sin(x0%data(4))/metric(:,10)
!        call assign_metric(f%b,transpose(metric))
        call assign_metric(f%u,transpose(metric))
        call assign_metric(f%b,transpose(metric))
        f%u%data(1)=sqrt(-1./(metric(:,1)+2.*metric(:,4)*omega+metric(:,10)* &
         omega*omega))
        f%u%data(4)=omega*f%u%data(1)
        f%b = calc_polvec(x0%data(2),cos(x0%data(3)),k0,dble(a),0.d0)
!        write(6,*) 'udotu: ',f%u*f%u,omega,f%u%data(1),r,rms
 !       write(6,*) 'numdisk fluidvars: ',allocated(f%u), &
!             allocated(f%b), allocated(f%rho)
        end subroutine get_numdisk_fluidvars

        subroutine get_phat_fluidvars(x0,k0,a,f)
        type (four_vector), intent(in), dimension(:) :: x0, k0
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
!        real :: rms,T0
        real, dimension(size(x0)) :: omega
        real, dimension(size(x0),10) :: metric
        real, dimension(size(x0),size(freq_tab)) :: fnu
!        write(6,*) 'phatdisk vals', size(x0%data(2))
        call phatdisk_vals(real(x0%data(2)),a,fnu,omega)
        f%fnu=fnu
        f%u%data(2)=0.; f%u%data(3)=0.
        f%b%data(1)=0.; f%b%data(2)=0.; f%b%data(3)=0.; f%b%data(4)=0.
!        write(6,*) 't: ',T,f%u%data(4),om,1./(r**3./2.+a),a*(r*r+a*a-d)/ &
        !  ((r*r+a*a)**2.-d*a*a*sin(x0%data(3))**2.)
        metric=kerr_metric(real(x0%data(2)),real(x0%data(3)),a)
!        call assign_metric(f%b,transpose(metric))
        call assign_metric(f%u,transpose(metric))
        f%u%data(1)=sqrt(-1./(metric(:,1)+2.*metric(:,4)*omega+metric(:,10)* &
         omega*omega))
        f%u%data(4)=omega*f%u%data(1)
        f%b = calc_polvec(x0%data(2),cos(x0%data(3)),k0,dble(a),0.d0)
!        write(6,*) 'udotu: ',f%u*f%u,omega,f%u%data(1),r,rms
        end subroutine get_phat_fluidvars

        subroutine get_ffjet_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f
        ! Computes properties of jet solution from Broderick & Loeb (2009)
        ! JAD 4/23/2010, fortran 3/30/2011
        call ffjet_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
!        write(6,*) 'ffjet u: ',f%u*f%u, f%b*f%b
        end subroutine get_ffjet_fluidvars

        subroutine get_harm_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f
        ! Computes properties of jet solution from Broderick & Loeb (2009)
        ! JAD 4/23/2010, fortran 3/30/2011
        call harm_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
!        write(6,*) 'harm u: ',f%u*f%u, f%b*f%b
        end subroutine get_harm_fluidvars

        subroutine get_harm3d_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f
        ! Computes properties of jet solution from Broderick & Loeb (2009)
        ! JAD 4/23/2010, fortran 3/30/2011
        call harm3d_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
!        write(6,*) 'harm u: ',f%u*f%u, f%b*f%b
        end subroutine get_harm3d_fluidvars

        subroutine get_thickdisk_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f
        ! Computes properties of jet solution from Broderick & Loeb (2009)
        ! JAD 4/23/2010, fortran 3/30/2011
        call thickdisk_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
!        write(6,*) 'thickdisk u: ',f%u*f%u, f%b*f%b
        end subroutine get_thickdisk_fluidvars

        subroutine get_mb09_fluidvars(x0,a,f)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        type (fluid), intent(inout) :: f
        ! Computes properties of jet solution from Broderick & Loeb (2009)
        ! JAD 4/23/2010, fortran 3/30/2011
        call mb09_vals(x0,a,f%rho,f%p,f%b,f%u,f%bmag)
!        write(6,*) 'mb09 u: ',f%u*f%u, f%b*f%b
        end subroutine get_mb09_fluidvars

! scale code to cgs units given mbh and simulation and cgs mdot values
        subroutine scale_sim_units(mbh,mdotcgs,mdot,rho,p,bmag, &
             ncgs,bcgs,tempcgs)
          real(kind=8), intent(in) :: mbh,mdot,mdotcgs
          real(kind=8) :: lcgs,tcgs
          real, dimension(:) :: rho,p,bmag
          real(kind=8), dimension(size(rho)) :: rhocgs,pcgs
          real(kind=8), dimension(size(rho)), intent(inout) :: ncgs,bcgs,tempcgs
        ! JAD 11/26/2012 adapted from IDL code
        ! Black hole mass sets time and length scales:
!        write(6,*) 'convert mbh: ',sp%mbh
        lcgs=GC*mbh*msun/c**2; tcgs=lcgs/c
        rhocgs=mdotcgs/mdot/lcgs**3*tcgs*rho; ncgs=rhocgs/mp
        ! Use this to convert pressure:
        pcgs=p*rhocgs/rho*c**2.
        ! Ideal gas temperature for single fluid (i.e., no separate e-/p):
        tempcgs=pcgs/ncgs/k
        ! And finally, bfield conversion is just square root of this:
        bcgs=bmag*sqrt(rhocgs/rho)*c
        ! Convert HL units to cgs:
        bcgs=bcgs*sqrt(4.*pi)
        end subroutine scale_sim_units

        subroutine monika_e(rho,p,b,beta_trans,rlow,rhigh,trat)
          real(kind=4), intent(in), dimension(:) :: rho,p,b
          real(kind=8), intent(in) :: rlow,rhigh,beta_trans
          real(kind=8), intent(inout), dimension(size(rho)) :: trat
          real(kind=8), dimension(size(rho)) :: beta,b2
          beta=p/(b*b)/0.5
!          beta_trans=1.d0
          b2=beta*beta
! defaults to 100, for now gmin even though should just be its own source param
!          rhigh=sp%gminval
!          rlow=1d0
!        Nhigh=1d0
!        Nlow=1d0
          where(b.gt.0d0)
             trat=rhigh*b2/(1d0+b2)+rlow/(1d0+b2)
!           nrat=Nhigh*b2/(1d0+b2)+Nlow/(1d0+b2)
          elsewhere
             trat=rhigh
!           nrat=Nhigh
          endwhere
! changed to add 1+trat to have maximum T_e = T_tot / 2
!          tempcgs=(p/rho)*mp*c*c/k/(1d0+trat)
      end subroutine monika_e

! non-thermal e- where jet energy density is high (e.g. Broderick & McKinney 2010, Dexter+2012)
        subroutine nonthermale_b2(alpha,gmin,p1,p2,bmagrho,bcgs,ncgsnth)
          real(kind=8), intent(in) :: alpha,gmin,p1,p2
          real(kind=8), intent(in), dimension(:) :: bmagrho,bcgs
          real(kind=8), intent(inout), dimension(:) :: ncgsnth
        where(bmagrho.gt.1.)
           ncgsnth=alpha*bcgs**2./8./pi/gmin*(p1-2.)/(p1-1.)/8.2e-7
        elsewhere
           ncgsnth=0.
        endwhere
        end subroutine nonthermale_b2

        subroutine convert_fluidvars_thickdisk(f,ncgs,ncgsnth,bcgs,tempcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tempcgs
        real(kind=8), dimension(size(f%rho)) :: trat
        real(kind=8) :: mdot,beta_trans
        mdot=0.0013; beta_trans=1d0
        call scale_sim_units(sp%mbh,sp%mdot,mdot,f%rho,f%p,f%bmag,ncgs, &
             bcgs,tempcgs)
        call monika_e(f%rho,f%p,f%bmag,beta_trans,1d0/sp%muval-1d0, &
             sp%gminval*(1d0/sp%muval-1d0),trat)
        tempcgs = tempcgs/(1d0+trat)
        call nonthermale_b2(sp%jetalphaval,sp%gminval,sp%p1,sp%p2, &
             f%bmag**2d0/f%rho,bcgs,ncgsnth)
        end subroutine convert_fluidvars_thickdisk

        subroutine convert_fluidvars_mb09(f,ncgs,ncgsnth,bcgs,tempcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tempcgs
        real(kind=8), dimension(size(f%rho)) :: trat
        real(kind=8) :: mdot,beta_trans
        mdot=.0013; beta_trans=1d0
        call scale_sim_units(sp%mbh,sp%mdot,mdot,f%rho,f%p,f%bmag,ncgs, &
             bcgs,tempcgs)
        call monika_e(f%rho,f%p,f%bmag,beta_trans,1d0/sp%muval-1d0, &
             sp%gminval*(1d0/sp%muval-1d0),trat)
        tempcgs = tempcgs/(1d0+trat)
! non-thermal e- by hand
        ncgsnth=ncgs
        end subroutine convert_fluidvars_mb09

        subroutine convert_fluidvars_harm(f,ncgs,ncgsnth,bcgs,tempcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tempcgs
        real(kind=8), dimension(size(f%rho)) :: trat
        real(kind=8) :: mdot,beta_trans
        mdot=.003; beta_trans=1d0
        call scale_sim_units(sp%mbh,sp%mdot,mdot,f%rho,f%p,f%bmag,ncgs, &
             bcgs,tempcgs)
! Moscibrodzka+2016 e- model with rlow = T_p / T_e from muval, rhigh = gmin*rlow
! reduces to T_p / T_e = const when gmin = 1 (should change input used for this)
! CHANGED TO USE MU INSTEAD OF TRAT to allow mu > 1/2 models to make sense 3/21/2017
        call monika_e(f%rho,f%p,f%bmag,beta_trans,sp%muval, &
             sp%muval/sp%gminval,trat)
        tempcgs = tempcgs*trat
        ncgsnth=ncgs
        end subroutine convert_fluidvars_harm

        subroutine convert_fluidvars_harm3d(f,ncgs,ncgsnth,bcgs,tempcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)), intent(out) :: ncgs,ncgsnth,bcgs,tempcgs
        real(kind=8), dimension(size(f%rho)) :: trat
        real(kind=8) :: mdot,beta_trans
        mdot=GC*sp%mbh*msun/c**3; beta_trans=1d0
        call scale_sim_units(sp%mbh,sp%mdot,mdot,f%rho,f%p,f%bmag,ncgs, &
             bcgs,tempcgs)
        call monika_e(f%rho,f%p,f%bmag,beta_trans,1d0/sp%muval-1d0, &
             sp%gminval*(1d0/sp%muval-1d0),trat)
        tempcgs = tempcgs/(1d0+trat)
        call nonthermale_b2(sp%jetalphaval,sp%gminval,sp%p1,sp%p2, &
             f%bmag**2d0/f%rho,bcgs,ncgsnth)
        end subroutine convert_fluidvars_harm3d

        subroutine convert_fluidvars_ffjet(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), &
          intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        type (source_params), intent(in) :: sp
        !real :: bfac=70., nfac=2.
!        write(6,*) 'sp: ',sp%nfac, sp%bfac
        ncgsnth=f%rho*sp%nfac; bcgs=f%bmag*sp%bfac; ncgs=0.; tcgs=0.
        end subroutine convert_fluidvars_ffjet

        subroutine convert_fluidvars_hotspot(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), &
          intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        type (source_params), intent(in) :: sp
        ncgs=f%rho; bcgs=f%bmag
        end subroutine convert_fluidvars_hotspot

        subroutine convert_fluidvars_schnittman_hotspot(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), &
          intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        type (source_params), intent(in) :: sp
        ncgs=f%rho; bcgs=1d0
        end subroutine convert_fluidvars_schnittman_hotspot

        subroutine convert_fluidvars_thindisk(f,tcgs,ncgs)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), intent(out) :: tcgs,ncgs
!        write(6,*) 'convert thindisk: ',size(f%rho),size(ncgs),size(tcgs)
        tcgs=f%rho
        ncgs=1.
        end subroutine convert_fluidvars_thindisk

        subroutine convert_fluidvars_numdisk(f,tcgs,ncgs)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(size(f%rho)), intent(out) :: tcgs,ncgs
!        write(6,*) 'convert numdisk: ',size(f%rho),size(ncgs),size(tcgs)
        tcgs=f%rho
        ncgs=1.
        end subroutine convert_fluidvars_numdisk

        subroutine convert_fluidvars_phatdisk(f,fnu,nu)
        type (fluid), intent(in) :: f
        real(kind=8), dimension(:,:), allocatable, intent(out) :: fnu
        real(kind=8), dimension(:), allocatable, intent(out) :: nu
        allocate(fnu(size(f%fnu,1),size(f%fnu,2)))
        allocate(nu(size(freq_tab)))
        fnu=f%fnu; nu=freq_tab
        end subroutine convert_fluidvars_phatdisk

        subroutine get_sphacc_fluidvars(x0,f)
        type (four_vector), intent(in), dimension(:) :: x0
        type (fluid), intent(inout) :: f
        real, dimension(size(x0)) :: u,grr,g00,n,B,T,ur
        u=1./x0%data(2)
        call sphacc_vals(u,n,B,T,ur)
       ! Equipartition B field
!        write(6,*) 'sphacc: ',B,size(x0)
        g00=-(1.-2.*u)
        grr=-1d0/g00
        f%u%data(2)=-ur
        f%u%data(3)=0d0
        f%u%data(4)=0d0
!        write(6,*) 'sphacc u: '!,-(grr*f%u%data(2)*f%u%data(2)+1)/g00
        f%u%data(1)=sqrt((-grr*f%u%data(2)*f%u%data(2)-1)/g00)
        f%b%data(3)=0d0
        f%b%data(4)=0d0
! use u dot b = 0, b dot b = B^2 to get components:
        f%b%data(1)=sqrt(f%u%data(2)**2*grr*B**2/ &
        (f%u%data(2)**2*g00*grr+f%u%data(1)**2*g00*g00))
        f%b%data(2)=-sqrt(B**2/grr-f%b%data(1)**2*g00/grr)
        f%rho=n
        f%p=T
        f%bmag=B
!        write(6,*) 'after sphacc u'
!        write(6,'((E9.4,2X))') u
!       write(6,*) 'T'
!       write(6,'((E9.4,2X))') T
!       write(6,*) 'n'
!       write(6,'((E9.4,2X))') n
!       write(6,*) 'B'
!       write(6,'((E9.4,2X))') B
        end subroutine get_sphacc_fluidvars

        subroutine convert_fluidvars_sphacc(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        type (source_params), intent(in) :: sp
        real(kind=8), dimension(size(f%rho)),  &
          intent(inout) :: ncgs,ncgsnth,bcgs,tcgs
!        write(6,*) 'convert', size(f%rho),size(ncgs)
!        write(6,*) 'convert',size(bcgs),size(tcgs)
!        write(6,*) 'convert',size(f%rho),size(f%bmag),size(f%p)
        ncgs=f%rho; bcgs=f%bmag; tcgs=f%p; ncgsnth=0.
!        write(6,*) 'after convert'
        end subroutine convert_fluidvars_sphacc

        subroutine get_hotspot_fluidvars(x0,a,f)
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
        type (four_vector), intent(in), dimension(:) :: x0
        type (four_Vector), dimension(size(x0)) :: x, xout
        real, dimension(size(x0)) :: n
        x=x0
! transform geodesic phi, t
        x%data(4) = -1*acos(0.) - x%data(4)
        x%data(1) = -1.*x%data(1)
!        write(6,*) 't hotspot: ',x%data(1)
        call hotspot_vals(x,a,n,f%b,f%u,xout)
        f%rho = dble(n)
        f%bmag = sqrt(f%b*f%b)
        end subroutine get_hotspot_fluidvars

        subroutine get_schnittman_hotspot_fluidvars(x0,a,f)
        type (fluid), intent(inout) :: f
        real, intent(in) :: a
        type (four_vector), intent(in), dimension(:) :: x0
        type (four_Vector), dimension(size(x0)) :: x, xout
        real, dimension(size(x0),10) :: metric
        real, dimension(size(x0)) :: n,omega,gfac,omt,ut,lc,hc,om,ar,d,safe
        real :: rms
        call schnittman_hotspot_vals(x0,a,n)
        !write(6,*) 'schnittman hotspot vars: ',n
        omega=1./(x0%data(2)**(3./2.)+a)
        f%rho = dble(n)
        f%u%data(2)=0.; f%u%data(3)=0.
        f%b%data(1)=0.; f%b%data(2)=0.; f%b%data(3)=0.; f%b%data(4)=0.
        metric=kerr_metric(real(x0%data(2)),real(x0%data(3)),a)
        f%u%data(1)=sqrt(-1./(metric(:,1)+2.*metric(:,4)*omega+metric(:,10)* &
         omega*omega))
        f%u%data(4)=omega*f%u%data(1)
        rms=calc_rms(a)
        d=x0%data(2)*x0%data(2)-2.*x0%data(2)+a*a
        lc=(rms*rms-2.*a*sqrt(rms)+a*a)/(rms**1.5-2.*sqrt(rms)+a)
        hc=(2.*x0%data(2)-a*lc)/d
        ar=(x0%data(2)*x0%data(2)+a*a)**2.-a*a*d*sin(x0%data(3))**2.
        om=2.*a*x0%data(2)/ar
        where(x0%data(2).gt.rms)
           omt=max(1./(x0%data(2)**(3./2.)+a),om)
        elsewhere
           omt=max((lc+a*hc)/(x0%data(2)*x0%data(2)+2.*x0%data(2)*(1.+hc)),om)
        endwhere
!        write(6,*) 'fluid hotspot assign xspot', omega, tspot, r0spot
        ut=sqrt(-1./(metric(:,1)+2.*metric(:,4)*omt+metric(:,10)* &
         omt*omt))
        safe=metric(:,1)+2.*metric(:,4)*omega+metric(:,10)* &
             omega*omega
        f%u%data(1)=merge(f%u%data(1),dble(ut),safe.lt.0d0)
        f%u%data(4)=merge(f%u%data(4),omt*f%u%data(1),safe.lt.0d0)
!        write(6,*) 'schnittman u safe: ',safe
!        write(6,*) 'schnittman u omt: ',omt
!        write(6,*) 'schnittman ut: ',f%u%data(1)
        gfac=1d0/sqrt((metric(:,10)*metric(:,1)-metric(:,4)*metric(:,4))* & 
           (metric(:,10)*f%u%data(4)*f%u%data(4)+f%u%data(1)* & 
           (2d0*metric(:,4)*f%u%data(4)+metric(:,1)*f%u%data(1))))
        f%b%data(1)=gfac*abs(metric(:,10)*f%u%data(4)+metric(:,4)*f%u%data(1))
        f%b%data(4)=-sign(1d0,metric(:,10)*f%u%data(4)+metric(:,4)*f%u%data(1)) &
           *(f%u%data(1)*metric(:,1)+metric(:,4)*f%u%data(4))*gfac
        call assign_metric(f%b,transpose(metric))
        call assign_metric(f%u,transpose(metric))
! Toroidal magnetic field:
        f%bmag = sqrt(f%b*f%b)
!        write(6,*) 'schnittman bmag u*u: ',f%u*f%u,f%bmag
        end subroutine get_schnittman_hotspot_fluidvars

        subroutine get_sariaf_fluidvars(x0,a,f)
!Semi-analytic RIAF model currently in progress. Inputs
        type (four_vector), intent(in), dimension(:) :: x0
        type (fluid), intent(inout) :: f
!inputs into sariaf_vals
        real, intent(in) :: a
        real, dimension(size(x0)) :: u,ctheta
!        real, dimension(size(x0)) :: riaf_u,riaf_neth,riaf_te,riaf_B
!outputs from sariaf_vals
        real, dimension(size(x0)) :: riaf_vr,riaf_vth,riaf_omega, &
             bmag,n,t,nnth
!        real, dimension(size(x0)) :: g00,grr,ur
!interim variables
        real, dimension(size(x0)) :: ub,aleph,bb
        real, dimension(size(x0)) :: rr, rho2,psi4,stheta
!rms related variables
        real :: rms
!r < rms intermediate variables
        real :: lambdae,game 
        real, dimension(size(x0)) :: delta,hhh 
        real(8), dimension(size(x0),10) :: metric
        real(8), dimension(size(x0)) :: hc,lc,d,ar,om,omtest,zero,vth,vphi,vr,u0!,one,two
        integer :: bl06
! debugging stuff
        real, dimension(size(x0)) :: rrcompare, ferret, ferretbb, ferretub
        integer :: i,alwingood,alwinbad,idlgood,idlbad
        real :: checkacc
        type (four_vector), dimension(size(x0)) :: fu

        rr = x0%data(2)
        rms = calc_rms(a)
        zero = 0d0

        lambdae = (rms**2. - 2.*a*sqrt(rms) + a**2.)/(rms**(3./2.)-2.*sqrt(rms)+a)
        game = sqrt(1.-2./3./rms)
        delta = rr*rr - 2.*rr + a*a
        hhh = (2.*rr-a*lambdae)/delta

        u=1./rr
        ctheta =cos(x0%data(3))
        stheta = sin(x0%data(3))
        rho2 = rr**2. + a**2. * (ctheta)**2.
        psi4 = 2.*rr/rho2
!kerr metric values copied from kerr.f90 kmetric_cov
!        gtt = -1.*(1.-psi4) !metric 1
!        gphi = (stheta*stheta) * (rho2 + a*a*(1.+2.*rr/rho2)*stheta*stheta) !metric 10
!        gtphi = -1.*psi4*a*stheta**2. !metric 4
        call sariaf_vals(a,ctheta,u,n,t,bmag,riaf_vr,riaf_vth,riaf_omega,nnth,bl06)
!       ! b = riaf_B
!        n = riaf_neth
!        t = riaf_te
!        u = riaf_u
!       ! Equipartition B field
!        write(6,*) 'sphacc: ',B,size(x0)
        
!        gtt= -(1.-2.*u)!alwinremark
!        gphi = !alwinremark
!        ub = gtt + riaf_omega*riaf_omega*gphi + 2.*riaf_omega*gtphi &
!             + grr*riaf_vr*riaf_vr + gtheta*riaf_vth*riaf_vth
!ub * u0**2. = -1
!        where(rr.lt.rms)
!           f%u%data(1) = game*(1.+2.*(1.+hhh)/rr)
!           f%u%data(2) = -1.*sqrt(2./3./rms)*(rms/rr-1.)**(3./2.)
!           f%u%data(3) = 0d0
!           f%u%data(4) = game*(lambdae + a*hhh)/rr/rr
!        elsewhere
!           f%u%data(1) =  sqrt(-1./ub) !what sign?
!           f%u%data(2) = 0d0 ! vr*f%u%data(1)
!           f%u%data(3) = 0d0 ! vth*f%u%data(1)
!           f%u%data(4) = riaf_omega * f%u%data(1)
!        endwhere

! construct four-velocity from three-velocity
        metric=kerr_metric(real(rr),real(x0%data(3)),a)
        lc = (rms**2d0-2d0*a*sqrt(rms)+a**2d0)/(rms**1.5d0-2d0*sqrt(rms)+a)
        d = rr**2d0-2d0*rr+a**2d0
        ar = (rr**2d0+a**2d0)**2d0-a**2d0*d*sin(x0%data(3))
        om = 2d0*a*rr/ar
        hc = (2d0*rr-a*lc)/d
!        if(bl06.eq.1) then
! stationary or free-fall inside ISCO, from bhimage.f
!           where(rr.lt.rms)
!              vr = zero
!              vr = sqrt(2d0*rr*(a**2d0+rr**2d0))*d/ar
!              vth = zero
!              vphi = om
!           elsewhere
!              vr = riaf_vr
!              vth = riaf_vth
!              vphi = riaf_omega
!           endwhere
!        else
!           where(rr.lt.rms)
! conserved quantities
!              vr = zero
!              vth = zero
!              omtest=(lc+a*hc)/(rr**2d0+2d0*rr*(1d0+hc))
!              vphi=merge(omtest,om,omtest.ge.om)
!           elsewhere
!              vr = riaf_vr
!              vth = riaf_vth
!              vphi = riaf_omega
!           endwhere
!        endif

        u0 = calc_u0(metric,dble(riaf_vr),dble(riaf_vth),dble(riaf_omega))
        fu = rms_vel(dble(a),x0%data(3),x0%data(2))
        where(rr.lt.rms)
           f%u%data(1) = fu%data(1)
           f%u%data(2) = fu%data(2)
           f%u%data(3) = fu%data(3)
           f%u%data(4) = fu%data(4)
        elsewhere
           f%u%data(1) = u0
           f%u%data(2) = riaf_vr*f%u%data(1)
           f%u%data(3) = riaf_vth*f%u%data(1)
           f%u%data(4) = riaf_omega*f%u%data(1)
        endwhere

!        f%u%data(1) = calc_u0(metric,vr,vth,vphi)
!        f%u%data(2) = vr*f%u%data(1)
!        f%u%data(3) = vth*f%u%data(1)
!        f%u%data(4) = vphi*f%u%data(1)

        aleph = -1.*(metric(:,4)*f%u%data(1)+metric(:,10)*f%u%data(4)) &
             /(metric(:,1)*f%u%data(1)+metric(:,4) * f%u%data(4))
!b0 = aleph*b_phi
        bb = metric(:,1)*aleph*aleph + metric(:,10) + 2.*metric(:,4)*aleph !I hope this is never negative
!bb*b_phi**2. = Bmag**2.
        f%b%data(4) = bmag/sqrt(bb) !what sign?
        f%b%data(3) = 0d0
        f%b%data(2) = 0d0
        f%b%data(1) = aleph * f%b%data(4)
        f%rho = n
        f%p = t
        f%bmag = bmag
        f%rho2 = nnth
!BEGIN TESTING 4 VECTOR ROUTINES
!testing dot product
!        where(rr.lt.rms)
!           rrcompare = 0.0
!        elsewhere
!           rrcompare = 1.0
!        endwhere
        checkacc = 1e-5
!        if(maxval(abs(gtt - metric(:,1))).gt.1e-5) then
!           write(6,*) 'ERROR: gtt Metric'
!        elseif(maxval(abs(gphi - metric(:,10))).gt.2e-4) then
!           write(6,*) 'ERROR: gphi Metric: ',maxval(abs(gphi-metric(:,10)))
!        elseif(maxval(abs(metric(:,4) - metric(:,4))).gt.1e-4) then
!           write(6,*) 'ERROR: metric(:,4) Metric'
!        else
!           write(6,*) 'METRIC GOOD'
!        endif
        call assign_metric(f%u,transpose(metric))
        call assign_metric(f%b,transpose(metric))
!        ferret = abs(f%u * f%u + 1.0)
!        ferretub = abs(f%u * f%b)
!        ferretbb = abs(f%b * f%b - bmag**2.)
!        alwinbad = 0
!        alwingood = 0
!        idlbad = 0
!        idlgood = 0
!        do i=1,size(x0)
!           if(rrcompare(i).gt.(0.5)) then
              !alwin stuff
!              if(ferret(i).gt.checkacc) then
!                 alwinbad = alwinbad + 1
!              else
!                 alwingood = alwingood + 1                
!                 if(ferretub(i).gt.checkacc) then
!                    write(6,*) 'WARNING: u dot b is wrong somewhere ',ferretub(i)
!                 elseif(ferretbb(i).gt.checkacc) then
!                       write(6,*) 'WARNING: b dot b is wrong somewhere ',ferretbb(i)
!                 endif
!              endif
!           else
!              if(ferret(i).gt.checkacc) then
!                 idlbad = idlbad + 1
!                 write(6,*) ferret(i)
!              else
!                 idlgood = idlgood + 1
!                 if(ferretub(i).gt.checkacc) then
!                    write(6,*) 'WARNING: u dot b is wrong somewhere ',ferretub(i)
!                 elseif(ferretbb(i).gt.checkacc) then
!                       write(6,*) 'WARNING: b dot b is wrong somewhere ',ferretbb(i)
!                 endif
!              endif
!           endif
!        enddo
!        if((idlgood+idlbad).gt.0) then
!        if((1.0*alwingood/(1.0*alwingood+1.0*alwinbad)).lt.(1.0*idlgood/(1.0*idlgood+1.0*idlbad))) then
!           write(6,*) 'ALWIN',alwingood,alwinbad
!        else
!           write(6,*) 'IDL',idlgood,idlbad
!        endif
!        elseif(alwinbad.gt.0) then
!           write(6,*) 'When 0 points inside RMS, ERROR: ', alwingood,alwinbad
!        endif
!END 4 VECTOR CHECKING ROUTINES
!        rrcompare = 1d0
!        write(6,*) 'u dot b: ',f%u*f%b
!        write(6,*) 'u dot u: ',f%u*f%u
!        if(maxval(abs(f%u * f%u + 1.0)).gt.checkacc) then
!           write(6,*) 'ALWIN WARNING: u dot u is wrong somewhere '
!           write(6,*) 'Error size: ',maxval(rrcompare*abs(f%u * f%u + 1.0))
!           write(6,*) 'u0: ',f%u%data(1)
!           write(6,*) 'vr: ',vr
!           write(6,*) 'vphi: ',vphi
!        endif
!        if(maxval((1.0-rrcompare)*abs(f%u * f%u + 1.0)).gt.checkacc) then
!           write(6,*) 'IDL WARNING: u dot u is wrong somewhere '
!           write(6,*) 'Error size: ',maxval((1.0-rrcompare)*abs(f%u * f%u + 1.0))
!        endif

!        if(maxval(abs(f%u * f%b)).gt.checkacc) then
!           write(6,*) 'WARNING: u dot b is wrong somewhere '
!           write(6,*) 'Error size: ',maxval(abs(f%u * f%b))
!           write(6,*) 'b0: ',f%b%data(1)
!           write(6,*) 'u0: ',f%u%data(1)
!        endif
!        if(maxval(abs(f%b * f%b - bmag**2.)).gt.checkacc) then
!           write(6,*) 'WARNING: b dot b is wrong somewhere ' 
!           write(6,*) 'Error size: ',maxval(abs(f%b * f%b - bmag**2.))
        end subroutine get_sariaf_fluidvars


        subroutine get_powerlaw_fluidvars(x0,a,f)
!powerlaw fluid model currently in progress (can also be one zone). Inputs
        type (four_vector), intent(in), dimension(:) :: x0
        type (fluid), intent(inout) :: f
!inputs into sariaf_vals
        real, intent(in) :: a
        real, dimension(size(x0)) :: u,ctheta
!        real, dimension(size(x0)) :: u,neth,te,B
!outputs from sariaf_vals
        real, dimension(size(x0)) :: vr,vth,omega, &
             bmag,n,t,nnth
!        real, dimension(size(x0)) :: g00,grr,ur
!interim variables
        real, dimension(size(x0)) :: ub,aleph,bb
        real, dimension(size(x0)) :: rr, rho2,psi4,stheta
!rms related variables
        real :: rms
!r < rms intermediate variables
        real :: lambdae,game 
        real, dimension(size(x0)) :: delta,hhh 
        real(8), dimension(size(x0),10) :: metric
        real(8), dimension(size(x0)) :: hc,lc,d,ar,om,omtest,zero,u0!,one,two
        type (four_vector), dimension(size(x0)) :: fu
        integer :: testindx

        rr = x0%data(2)
        rms = calc_rms(a)
        zero = 0d0

        lambdae = (rms**2. - 2.*a*sqrt(rms) + a**2.)/(rms**(3./2.)-2.*sqrt(rms)+a)
        game = sqrt(1.-2./3./rms)
        delta = rr*rr - 2.*rr + a*a
        hhh = (2.*rr-a*lambdae)/delta

        u=1./rr
        ctheta =cos(x0%data(3))
        stheta = sin(x0%data(3))
        rho2 = rr**2. + a**2. * (ctheta)**2.
        psi4 = 2.*rr/rho2

        call powerlaw_vals(a,ctheta,u,n,t,bmag,vr,vth,omega,nnth)

!        testindx = maxloc(bmag,1)
!        write(6,*) 'powerlaw fluidvars x0: ',x0(testindx)%data(2), &
!             x0(testindx)%data(3),x0(testindx)%data(4)

! construct four-velocity from three-velocity

        metric = kerr_metric(x0%data(2),x0%data(3),dble(a))
        u0 = calc_u0(metric,dble(vr),dble(vth),dble(omega))
   !     fu = rms_vel(dble(a),x0%data(3),x0%data(2))
!        where(rr.lt.rms)
!           f%u%data(1) = fu%data(1)
!           f%u%data(2) = fu%data(2)
!           f%u%data(3) = fu%data(3)
!           f%u%data(4) = fu%data(4)
!        elsewhere
!           f%u%data(1) = u0
!           f%u%data(2) = vr*f%u%data(1)
!           f%u%data(3) = vth*f%u%data(1)
!           f%u%data(4) = omega*f%u%data(1)
!        endwhere

        f%u%data(1) = u0
        f%u%data(2) = vr*u0
        f%u%data(3) = vth*u0
        f%u%data(4) = omega*u0

        aleph = -1.*(metric(:,4)*f%u%data(1)+metric(:,10)*f%u%data(4)) &
             /(metric(:,1)*f%u%data(1)+metric(:,4) * f%u%data(4))

        bb = metric(:,1)*aleph*aleph + metric(:,10) + 2.*metric(:,4)*aleph !I hope this is never negative

        f%b%data(4) = bmag/sqrt(bb) !what sign?
        f%b%data(3) = 0d0
        f%b%data(2) = 0d0
        f%b%data(1) = aleph * f%b%data(4)
        f%rho = n
        f%p = t
        f%bmag = bmag
        f%rho2 = nnth

        call assign_metric(f%u,transpose(metric))
        call assign_metric(f%b,transpose(metric))

        end subroutine get_powerlaw_fluidvars


        subroutine convert_fluidvars_sariaf(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8) :: riaf_n0, riaf_beta
!        real(kind=8) :: riaf_t0 !no longer used
        real(kind=8), dimension(size(f%rho)), &
             intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        type (source_params), intent(in) :: sp
!        riaf_n0 = sp%mdot !4.e7 
!        riaf_t0 = 1.6d11
! either non-unity here or in fluid_model_sariaf.f90
!riaf_t0 = 1.0 because for sariaf muval = t 
!and it gets changed in emis!sp%muval so not here
!to avoid double counting 
!Old value was 1.6e11 !these will change to be based on sp, I think.
!        riaf_beta = 1d1
!        ncgs= riaf_n0 * f%rho
        ncgs = f%rho
!        bcgs= sqrt(riaf_n0 / riaf_beta) * f%bmag
        bcgs = f%bmag
! new var rho2 is used to hold a second density, in this non-thermal e-
        ncgsnth= f%rho2
        tcgs= f%p !* riaf_t0
!        f%b%data(1) = sqrt(riaf_n0 / riaf_beta) * f%b%data(1)
!        f%b%data(4) = sqrt(riaf_n0 / riaf_beta) * f%b%data(4)
!        not necessary to scale f%b because f%b is only used for an angle
        end subroutine convert_fluidvars_sariaf

        subroutine convert_fluidvars_powerlaw(f,ncgs,ncgsnth,bcgs,tcgs,sp)
        type (fluid), intent(in) :: f
        real(kind=8) :: n0, beta, beta_trans
        real(kind=8), dimension(size(f%rho)), &
             intent(out) :: ncgs,ncgsnth,bcgs,tcgs
        real(kind=8), dimension(size(f%rho)) :: trat
        type (source_params), intent(in) :: sp
        beta_trans = 1d0
        ncgs = f%rho
        bcgs = f%bmag
        ncgsnth= f%rho2
        call monika_e(f%rho,f%rho,f%bmag,beta_trans,1d0/sp%muval-1d0, &
             sp%gminval*(1d0/sp%muval-1d0),trat)
        tcgs= f%p/(1d0+trat)
        end subroutine convert_fluidvars_powerlaw

! source param routines
        subroutine assign_source_params_type(sp,type)
          character(len=20), intent(in) :: type
          type (source_params), intent(inout) :: sp
          if(type=='const') then
             sp%type=CONST
          elseif(type=='tail') then
             sp%type=TAIL
          else
             write(6,*) 'ERROR in assign_source_params_type: type not recognized', type
          endif
        end subroutine assign_source_params_type

        subroutine initialize_source_params(sp,nup)
          type (source_params), intent(inout) :: sp
          integer, intent(in) :: nup
          allocate(sp%gmin(nup))
          allocate(sp%jetalpha(nup))
          allocate(sp%mu(nup))
        end subroutine initialize_source_params

        subroutine del_source_params(sp)
          type (source_params), intent(inout) :: sp
          deallocate(sp%gmin)
          deallocate(sp%jetalpha)
          deallocate(sp%mu)
        end subroutine del_source_params
        
        subroutine assign_source_params(sp,ncgs,tcgs,ncgsnth)
          type (source_params), intent(inout) :: sp
          real(kind=8), dimension(:), intent(in) :: ncgs,tcgs
          real(kind=8), dimension(:), intent(inout) :: ncgsnth
          real(kind=8), dimension(size(ncgs)) :: x,one,gmin,gmax,zero,factor
          zero=0d0
          one=1d0
          gmax=sp%gmax
!          write(6,*) 'fluid assign source params: ',sp%type,CONST,TAIL
          select case(sp%type)
             case (CONST)
                sp%gmin=sp%gminval
                sp%jetalpha=sp%jetalphaval
                sp%mu=sp%muval
!                write(6,*) 'fluid assign source params: ',sp%gminval,sp%gmin
             case (TAIL)
                sp%jetalpha=sp%jetalphaval
                sp%mu=sp%muval
!Fluid -> calc_gmin -> mu*tcgs -> emis means that calc_gmin is given a tcgs pre-mu correction and needs to be corrected here.                 
                call calc_gmin_subroutine(sp%p2,sp%mu*k*tcgs/m/c/c,sp%jetalpha,gmin,x)
!                sp%gmin=merge(merge(gmin,one,gmin.ge.1d0),gmax/2d0,gmin.le.gmax)
!                trust calc_gmin to merge gmin < 1 now
                sp%gmin=merge(gmin,gmax/2d0,gmin.le.gmax)
                factor=merge(one,(gmax/2d0/gmin)**(sp%p2 - 2.),gmin.le.gmax)
!when gmin is corrected for being too large, multiply ncgsnth by a corrective factor. The correction (1-p) is already applied, so the correction p-2 is needed.
                ncgsnth=factor * merge(x*ncgs*sp%gmin**(1.-sp%p2),zero,x.gt.0d0)
!ncgsnth proportional to gamma**-1
!                sp%gmin=gmin
!                where(ncgsnth.ne.ncgsnth)
!                   ncgsnth=-1e8
!                endwhere
!                write(6,*) 'gmin ncgsnth tail: ',minval(sp%gmin),maxval(sp%gmin),minval(ncgsnth),maxval(ncgsnth)
!                write(6,*) 'gmin ncgsnth tail: ',gmin
!                write(6,*) 'gmin ncgsnth tail x: ',x
!                write(6,*) 'gmin ncgsnth tail: ',ncgs
!                write(6,*) 'gmin ncgsnth tail: ',ncgsnth
          end select
        end subroutine assign_source_params

      end module fluid_model
