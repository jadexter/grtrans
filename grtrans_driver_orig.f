       include 'geodesics.f'
       include 'fluid.f'
       include 'rad_trans.f'
       include 'read_inputs.f'
       include 'camera.f'
       include 'emis.f'
       include 'odepack.f'

       module grtrans

       use fluid_model
       use class_rad_trans
       use geodesics
       use emissivity
       use odepack

       implicit none

! For now make it global till I figure out how to pass through ODEPACK...
! grtrans object is made up of a fluid_model, a geodesic and a rad_trans object:
!       type grtrans
        integer :: NK=7
        integer :: IS_LINEAR_STOKES=1 
        double precision :: ortol=1.d-8, oatol=1.d-10
         type (fluid) :: f
         type (rad_trans) :: r
         type (geo) :: g
         type (geokerr_args) :: gargs
         type (emis) :: e
!       end type grtrans

       interface grtrans_driver_jac_form
         module procedure grtrans_driver_jac_form
       end interface

       interface grtrans_driver_rhs_form
         module procedure grtrans_driver_rhs_form
       end interface

       contains

         subroutine grtrans_driver(gargs,gunit)
! Driver routine for GRTrans. JAD 2/27/2011
         implicit none
!         double precision, dimension(2) :: pixel
         type (geokerr_args), intent(in) :: gargs
         integer, intent(in) :: gunit
!         type (motion_constants) :: pixel
!         write(6,*) 'geodesic'
         call initialize_geodesic(g,gargs,gunit)
         call initialize_emissivity(e,r%neq,NK)
! Call integration routine:
!         write(6,*) 'integrate'
         if(g%npts.ne.1) then
           call grtrans_driver_integrate()
         else
           call grtrans_compute_intensity()
         endif
!         write(6,*) 'delete'
         call del_geodesic(g)
         call del_emissivity(e)
         return
         end subroutine grtrans_driver

         subroutine grtrans_compute_intensity()

         end subroutine grtrans_compute_intensity
 
         subroutine grtrans_driver_integrate()
! Subroutine to call integration routine
! JAD 2/28/2011
         if (r%iflag==0) call grtrans_driver_integrate_lsoda()
         return
         end subroutine grtrans_driver_integrate

         subroutine grtrans_driver_integrate_lsoda()
!         external grtrans_driver_calc_rhs
!         external grtrans_driver_calc_jac
!         call lsoda_wrapper(grtrans_driver_rhs_form,
!     &   r%neq,grtrans_driver_jac_form)
         double precision, dimension(4) :: I0
!         write(6,*) 'lsoda basic', size(r%I), size(g%lambda)
!         write(6,*) 'lambda: ',g%lambda, I0
         I0=0d0
!         write(6,*) 'I0: ',I0
         call lsoda_basic(grtrans_driver_lsoda_calc_rhs,I0,g%lambda,oatol,
     &    ortol,grtrans_driver_lsoda_calc_jac,r%I)
!         write(6,*) 'integrate: ',r%I(1,:)-g%lambda
         return
         end subroutine grtrans_driver_integrate_lsoda

         subroutine grtrans_driver_lsoda_calc_rhs(neq,lam,I,dIdlam)
! Compute RHS dIdlam for LSODA
         integer, intent(in) :: neq
         double precision, intent(in) :: lam
         double precision, intent(in), dimension(neq) :: I
         double precision, intent(out), dimension(neq) :: dIdlam
!         double precision, intent(out), dimension(neq,neq) :: jac
         double precision, dimension(neq) :: j
         double precision, dimension(NK) :: K
         call grtrans_driver_aux(neq,lam,j,K)
         call grtrans_driver_rhs_form(neq,j,K,dIdlam,I)
!         write(6,*) 'dIdlam: ',dIdlam
         return 
         end subroutine grtrans_driver_lsoda_calc_rhs

         subroutine grtrans_driver_rhs_form(neq,j,K,dIdlam,I)
         integer, intent(in) :: neq
         double precision, intent(in), dimension(neq) :: j
         double precision, intent(in), dimension(NK) :: K
         double precision, intent(out), dimension(neq) :: dIdlam
         double precision, intent(in), dimension(neq) :: I
         if (IS_LINEAR_STOKES==1) then
           dIdlam(1)=j(1)-(K(1)*I(1)+K(2)*I(2)+K(3)*I(3)+K(4)*I(4))
           dIdlam(2)=j(2)-(K(2)*I(1)+K(1)*I(2)+K(7)*I(3)-K(6)*I(4))
           dIdlam(3)=j(3)-(K(3)*I(1)-K(7)*I(2)+K(1)*I(3)+K(5)*I(4))
           dIdlam(4)=j(4)-(K(4)*I(1)+K(6)*I(2)-K(5)*I(3)+K(1)*I(4))
!           write(6,*) 'made it'
         endif
!         write(6,*) 'dIdlam ', j, K, I, dIdlam
         return
         end subroutine grtrans_driver_rhs_form

         subroutine grtrans_driver_jac_form(neq,j,K,nrowpd,pd)
         integer, intent(in) :: neq, nrowpd
         double precision, intent(in), dimension(neq) :: j
         double precision, intent(in), dimension(NK) :: K
         double precision, intent(out), dimension(nrowpd,neq) :: pd
         if (IS_LINEAR_STOKES==1) then
           pd(1,1)=K(1)
           pd(1,2)=K(2)
           pd(1,3)=K(3)
           pd(1,4)=K(4)
           pd(2,1)=K(2)
           pd(2,2)=K(1)
           pd(2,3)=K(7)
           pd(2,4)=-K(6)
           pd(3,1)=K(3)
           pd(3,2)=-K(7)
           pd(3,3)=K(1)
           pd(3,4)=K(5)
           pd(4,1)=K(4)
           pd(4,2)=K(6)
           pd(4,3)=-K(5)
           pd(4,4)=K(1)
           pd=-1d0*pd
         endif
!         write(6,*) 'pd: ',pd
         return
         end subroutine grtrans_driver_jac_form
       
         subroutine grtrans_driver_lsoda_calc_jac(neq,lam,ml
     &     ,mu,pd,nrowpd)
! Compute Jacobian for LSODA
         integer, intent(in) :: neq, nrowpd
         double precision, intent(in) :: lam
         double precision, intent(in) :: ml
         double precision, intent(in) :: mu
         double precision, intent(out), dimension(nrowpd,neq) :: pd
         double precision, dimension(neq) :: j
         double precision, dimension(NK) :: K
         call grtrans_driver_aux(neq,lam,j,K)
         call grtrans_driver_jac_form(neq,j,K,nrowpd,pd)
!         write(6,*) 'pd: ', pd
         return
         end subroutine grtrans_driver_lsoda_calc_jac

         subroutine grtrans_driver_aux(neq,lam,j,K)
         integer, intent(in) :: neq
         double precision, intent(in) :: lam
         double precision, intent(out), dimension(neq) :: j
         double precision, intent(out), dimension(NK) :: K
!         call interp_geo_point(lam,g)
!         call get_fluid_vars(g%x0,f)
!         call calc_emissivity(r%nu,f,e)
         call calc_emissivity(e)
         j=e%j
         K=e%K
!         call calc_rad_trans_coefs(e%j,e%K)
!         j=1d0
!         K=0d0
         return
         end subroutine grtrans_driver_aux

       end module grtrans

       subroutine grtrans_main()
       use grtrans_inputs
       use grtrans, only: grtrans_driver, r, f
       use ray_trace
       use fluid_model!, only: initialize_fluid_model,del_fluid_model
       use class_rad_trans!, only: initialize_rad_trans, del_rad_trans
       use geodesics!, only: initialize_pixels
       implicit none
       integer :: inum, gunit, i
       type (ray_set) :: c
!       type (rad_trans) :: r
!       type (fluid) :: f
       type (geokerr_args) :: gargs
       character(len=20) :: outfile='grtrans.out', ifile='inlist'
       call read_inputs()
!       write(6,*) 'after inputs'
       call initialize_raytrace_camera(c,nro,nphi,nvals)
!       write(6,*) 'camera'
       call initialize_pixels(gargs,use_geokerr,standard,mu0,a,rcut,
     & nrotype,a1,a2,b1,b2,nro,nphi,nup)
!       write(6,*) 'pixels'
       call initialize_fluid_model(f,fname)
!       write(6,*) 'fluid'
!       call initialize_rad_trans(r,iname)
!       write(6,*) 'rad_trans'
       do i=1,c%nx*c%ny
!         write(6,*) 'loop'
         call initialize_rad_trans(r,iname)
         call raytrace_camera_pixel_getnext(c)
         gargs%geonum=i
!         write(6,*) 'driver'
         call grtrans_driver(gargs,gunit)
!         write(6,*) 'save'
         call save_raytrace_camera_pixel(c,(/real(gargs%alpha(i)),
     &   real(gargs%beta(i))/),real(r%I(:,r%npts)))
         call del_rad_trans(r)
       enddo
!       write(6,*) 'after loop', c%pixvals
       call write_raytrace_camera(c,12,outfile,0)
!       call del_rad_trans(r)
       call del_fluid_model(f)
       call del_geokerr_args(gargs)
       call del_raytrace_camera(c)
       return
       end subroutine grtrans_main

       program call_grtrans
       implicit none
       call grtrans_main()
       end program
