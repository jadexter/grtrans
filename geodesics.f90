!        include 'geokerr_wrapper.f'

       module geodesics
       
       use class_four_vector
       use interpolate, only: get_weight
       use kerr, only: calc_nullp, kerr_metric, bl2ks
!       use class_geokerr, only: geokerr

       implicit none

       type geokerr_args
          real(kind=8) :: u0,uout,uin,offset,mu0,phi0
          integer :: usegeor,phit,mufill,ncase,kext,next, &
               nup
          real(kind=8) :: a
          real(kind=8), dimension(:), allocatable :: uf,muf, &
                                                     su,sm,alpha,beta, &
                                                     l,q2,t0
          integer, dimension(:), allocatable :: tpm,tpr
       end type geokerr_args

       type geo
          type (four_Vector) :: x0,k0
          type (geokerr_args) :: gk
          real(kind=8), dimension(:), allocatable :: lambda
          real(kind=8), dimension(:), allocatable :: tpmarr, tprarr
          integer :: lindx=25,npts
          type (four_Vector),  dimension(:), allocatable :: x,k
       end type geo

       interface construct_geodesic
         module procedure load_geodesic
         module procedure construct_geodesic
       end interface

       interface initialize_pixels
         module procedure initialize_pixels
       end interface
       
       interface initialize_geo_tabs
          module procedure initialize_geo_tabs
       end interface

       interface  get_pixel_locations
          subroutine initialize_camera_geokerr(standard,a1,a2, &
               b1,b2,rcut, &
               nrotype,nro,nphi,nup,u0,uout,mu0,a,offset, &
               ufarr,mufarr,alarr,bearr,q2arr,larr,smarr,suarr, &
               tpmarr,tprarr)
            integer, intent(in) :: nro,nphi,nup,nrotype,standard
            real(kind=8), intent(in) :: a1,a2,b1,b2,rcut,a,uout
            real(kind=8), intent(out), dimension(nro*nphi) :: alarr, &
                 bearr,ufarr,mufarr, &
                 q2arr,larr,smarr,suarr
            integer, intent(out), dimension(nro*nphi) :: tpmarr,tprarr
            real(kind=8), intent(inout) :: u0,mu0
            real(kind=8), intent(out) :: offset
          end subroutine initialize_camera_geokerr
       end interface get_pixel_locations

       contains

         subroutine initialize_pixels(geargs,use_geokerr,standard,mu0,phi0,a, &
                            uout,uin,rcut,nrotype,a1,a2,b1,b2,nro,nphi,nup)
         integer, intent(in) :: nro,nphi,nup,nrotype,standard
         type (geokerr_args), intent(inout) :: geargs
         real(kind=8), intent(in) :: a1,a2,b1,b2,rcut,a,uout,uin,mu0,phi0
         logical :: use_geokerr
!         write(6,*) 'a: ',a
         geargs%a=a
         geargs%mu0=mu0
         geargs%phi0=phi0
         geargs%nup=nup
         geargs%uout=uout
         if(USE_GEOKERR) then
            if(standard==1) then
              geargs%usegeor=0
              geargs%phit=1
              if(nup.gt.1) then
                 geargs%mufill=1
!                 geargs%kext=3
!                 geargs%next=60
! provisional idea to scale these with nup to improve accuracy for large nup
                 geargs%kext=ceiling(3d0*geargs%nup/400d0)
                 geargs%next=ceiling(60d0*geargs%nup/400d0)
                 write(6,*) 'kext next: ',geargs%kext,geargs%next
              else
                 geargs%mufill=0
                 geargs%kext=0
                 geargs%next=0
              endif
            else
              geargs%mufill=0
              geargs%phit=1
              geargs%usegeor=1
              geargs%kext=0
              geargs%next=0
            endif
!            write(6,*) 'get_pixel_loc', geargs%kext, geargs%next
!            write(6,*) 'geo pixel loc', a,standard,a1,a2,b1,b2
            call get_pixel_locations(standard,a1,a2,b1,b2,rcut, &
            nrotype,nro,nphi,nup,geargs%u0,geargs%uout,geargs%mu0,a &
            ,geargs%offset, &
            geargs%uf,geargs%muf &
            ,geargs%alpha,geargs%beta,geargs%q2,geargs%l, &
            geargs%sm,geargs%su,geargs%tpm,geargs%tpr)
!            write(6,*) 'pixel loc: ',a,geargs%mu0
!            write(6,*) 'after',geargs%alpha,geargs%beta,geargs%q2,geargs%l
         endif
         end subroutine initialize_pixels

         subroutine initialize_geo_tabs(gargs,i)
         type (geokerr_args), intent(inout) :: gargs
         integer, intent(in) :: i
         integer :: status
         type (geo) :: g
         real(8), dimension(1) :: t0
!         write(6,*) 'init geo tabs'
         gargs%kext=0; gargs%next=0; gargs%mufill=0
         call initialize_geodesic(g,gargs,i,status)
!         write(6,*) 'call get_geo_tabs'
         call get_geo_tabs(g,t0)
         gargs%t0(i)=t0(1)
!         write(6,*) 'gargs assign: ',i,gargs%t0(i)
         call del_geodesic(g)
!         write(6,*) 't0: ',gargs%t0
         end subroutine initialize_geo_tabs

         subroutine get_geo_tabs(g,t0,ks)
         ! Get the initial time    
         type (geo), intent(in) :: g
         integer, intent(in), optional :: ks
         real(8), dimension(1), intent(out) :: t0
           if(g%npts.eq.0) then
              t0=0.
           else
              if(present(ks)) then
                 t0=bl2ks(g%x%data(2),g%x%data(1),g%gk%a,1)
              else
                 t0=g%x%data(1)
              endif
           endif
         end subroutine get_geo_tabs

         subroutine construct_geodesic(g)
         type (geo), intent(inout) :: g
!         type (motion_constants), intent(inout) :: pixel
!         write(6,*) 'nup: ',g%gk%q2,g%gk%l
         call geokerr_wrapper(g)
   !      write(6,*) 'after wrapper'
         return
         end subroutine construct_geodesic 

         subroutine load_geodesic(g,gunit)
         type (geo), intent(inout) :: g
         integer, intent(in) :: gunit
         real(kind=8) :: ufi,mufi,dti,dphii,lambdai,tpmi,tpri
         real(kind=8) :: alpha,beta,nup
         integer :: k
         write(6,*) 'load geodesic: ',gunit
         read(gunit,*) alpha,beta,nup
         write(6,*) 'alphabeta: ',alpha,beta
         g%gk%alpha=alpha
         g%gk%beta=beta
         g%npts=nup
         write(6,*) 'load geo: ',g%npts,nup
         allocate(g%x(g%npts))
         allocate(g%lambda(g%npts))
         allocate(g%tpmarr(g%npts))
         allocate(g%tprarr(g%npts))
         allocate(g%k(g%npts))
         do k=1,g%npts
           read(gunit,*) ufi,mufi,dti,dphii,lambdai,tpmi,tpri
           write(6,*) 'after read: ',allocated(g%x)
           g%x(k)%data(2)=1./ufi
           g%x(k)%data(3)=acos(mufi)
           g%x(k)%data(1)=dti
           g%x(k)%data(4)=dphii
           g%lambda(k)=lambdai
           g%tpmarr(k)=tpmi
           g%tprarr(k)=tpri
         end do
!         g%gk%t0=g%x(1)%data(1)
         g%x(1)%data(1)=0d0
         return
         end subroutine load_geodesic

         subroutine initialize_geokerr_args(gk,n)
         integer, intent(in) :: n
         type (geokerr_args), intent(inout) :: gk
!         write(6,*) 'init geokerr args: ',allocated(gk%muf),allocated(gk%uf)
         allocate(gk%muf(n))
         allocate(gk%uf(n))
         allocate(gk%tpm(n))
         allocate(gk%tpr(n))
         allocate(gk%su(n))
         allocate(gk%sm(n))
         allocate(gk%alpha(n))
         allocate(gk%beta(n))
         allocate(gk%q2(n))
         allocate(gk%l(n))
         allocate(gk%t0(n))
         end subroutine initialize_geokerr_args

         subroutine initialize_geodesic(g,gargs,i,status,gunit)
!        type (motion_constants), intent(in) :: pixel
         type (geo), intent(inout) :: g
         type (geokerr_args), intent(in) :: gargs
         integer, intent(in) :: i
         integer, intent(out) :: status
         integer, intent(in), optional :: gunit
         status=1
!         g%m=pixel
!         write(6,*) 'initialize_geodesic from init geo', i,allocated(g%gk%muf)
         call initialize_geokerr_args(g%gk,1)
         g%gk%mu0=gargs%mu0
         g%gk%phi0=acos(-1.)*gargs%phi0
         g%gk%a=gargs%a
         g%gk%t0=gargs%t0(i)
!         write(6,*) 'a: ',gargs%a,g%gk%a
!         write(6,*) 'tpr: ', gargs%geonum
         g%gk%u0=gargs%u0
         g%gk%uf=gargs%uf(i)
         g%gk%uout=gargs%uout
         g%gk%uin=gargs%uin
         g%gk%alpha=gargs%alpha(i)
         g%gk%beta=gargs%beta(i)
         g%gk%su=gargs%su(i)
         g%gk%sm=gargs%sm(i)
         g%gk%l=gargs%l(i)
         g%gk%q2=gargs%q2(i)
         g%gk%tpm=gargs%tpm(i)
         g%gk%tpr=gargs%tpr(i)
         g%gk%usegeor=gargs%usegeor
         g%gk%mufill=gargs%mufill
 !        write(6,*) 'gargs: ',g%gk%mufill,gargs%nup,gargs%next,gargs%kext
         g%npts=gargs%nup+gargs%next-2*gargs%kext
!         write(6,*) 'init: ',gargs%nup,gargs%next,gargs%kext,g%npts
!         write(6,*) 'initialize geo: ',g%npts,gargs%nup,gargs%next,gargs%kext
         g%gk%nup=gargs%nup
         g%gk%offset=gargs%offset
         g%gk%phit=gargs%phit
         g%gk%muf=gargs%muf(i)
         g%gk%kext=gargs%kext
         g%gk%next=gargs%next
!         write(6,*) 'present: ',present(gunit)
         if(present(gunit)) then
            write(6,*) 'load geodesic: ',gunit
            call construct_geodesic(g,gunit)
         else
            call construct_geodesic(g)
         endif
!         call construct_geodesic(g)
         if(sum(g%x%data(2)).lt.0.) status=-1
         return
         end subroutine initialize_geodesic

         subroutine interp_geo_point(lam,g)
         type (geo), intent(inout) :: g
         real(kind=8), intent(in) :: lam
         real(kind=8) :: weight
         call get_weight(g%lambda,lam,g%lindx,weight)
         g%x0=(1d0-weight)*g%x(g%lindx)+weight*g%x(g%lindx+1)
!         write(6,*) 'after'
         g%k0=(1d0-weight)*g%k(g%lindx)+weight*g%k(g%lindx+1)
!         write(6,*) 'k'
         return
         end subroutine interp_geo_point

         subroutine geokerr_wrapper(g)
         type (geo), intent(inout) :: g
!         type (motion_constants), intent(in) :: pixel
         real(kind=8), dimension(:), allocatable :: &
                              ufi,mufi, &
                              dti,dphi,lambdai
         integer, dimension(:), allocatable :: tpmi,tpri
         integer :: nupdum,ii,i1,i2
         nupdum=g%gk%nup
!         write(6,*) 'geokerr nup: ',nupdum,g%gk%nup
         allocate(ufi(g%npts)); allocate(mufi(g%npts))
         allocate(dti(g%npts)); allocate(dphi(g%npts))
         allocate(lambdai(g%npts)); allocate(tpmi(g%npts))
         allocate(tpri(g%npts))
!         write(6,*) 'geokerr_wrapper npts: ',g%gk%nup,g%npts
         call geokerr(g%gk%u0,g%gk%uf,g%gk%uout,g%gk%mu0, &
          g%gk%muf,g%gk%a,g%gk%l,g%gk%q2, &
          g%gk%alpha,g%gk%beta,g%gk%TPM,g%gk%TPR, &
             g%gk%SU,g%gk%SM,g%gk%nup,g%gk%OFFSET, &
         g%gk%phit,g%gk%usegeor,g%gk%mufill,g%gk%ncase, &
         g%gk%kext,g%gk%next, &
         UFI,MUFI,DTI,DPHI,TPMI,TPRI,LAMBDAI)
!         write(6,*) 'after geokerr: ',TPMI
          g%npts=g%gk%nup
!         write(6,*) 'gnpts: ',g%npts
! remove locations with invalid solutions:
         i1=1; i2=g%npts
         if (g%gk%usegeor.eq.1.and.g%npts.gt.1) then
         if (ufi(1).lt.0d0) then
! Find first positive element:
            do ii=2,g%npts
               if (ufi(ii).ge.0d0) then
                  i1=ii
                  exit
               endif
            enddo
         else
            i1=1
         endif
!         i1=1
         if (ufi(g%npts).lt.0d0) then
! Find last non-zero element:
            do ii=1,g%npts-1
               if (ufi(g%npts-ii).ge.0d0) then
                  i2=g%npts-ii
                  exit
               endif
            enddo
         else
            i2=g%npts
         endif
         g%npts=i2-i1+1
         endif
!         write(6,*) 'i1: ',i1,i2,g%npts
!         write(6,*) 'uf: ',ufi
!         write(6,*) 'geokerr_wrapper: ',allocated(g%x),allocated(g%lambda),allocated(g%tpmarr)
!         write(6,*) 'geokerr_wrapper 2: ',allocated(g%tprarr),allocated(g%k)
         allocate(g%x(g%npts))
         allocate(g%lambda(g%npts))
         allocate(g%tpmarr(g%npts))
         allocate(g%tprarr(g%npts))
         allocate(g%k(g%npts))
         g%x%data(2)=1d0/UFI(i1:i2)
!         write(6,*) 'geokerr wrapper: ',g%npts,mufi,tpmi
         g%x%data(3)=ACOS(MUFI(i1:i2))
         g%x%data(4)=g%gk%phi0-DPHI(i1:i2)
! fixes for pole-on viewing and strange turning point behaviors from IDL code
! if this is a pole on case then need to correct phi:
         if (abs(g%gk%mu0).eq.1d0) then
            g%x%data(4)=g%x%data(4)+sign(1d0,g%gk%mu0)*atan2(g%gk%beta(1),g%gk%alpha(1))
         endif
! do differences forwards in time rather than backwards JAD 1/14/2013
         if (g%npts.ne.1) then 
           g%lambda=LAMBDAI(g%npts)-LAMBDAI(i1:i2)
           g%x%data(1)=DTI(1)-DTI(i1:i2)-g%gk%t0(1)
           g%tpmarr=TPMI(i1:i2)
           g%tprarr=TPRI(i1:i2)
         else
           g%lambda=LAMBDAI(1)
           g%x%data(1)=DTI(1)
           g%tpmarr=TPMI(1)
           g%tprarr=TPRI(1)
         endif
!         write(6,*) 'TPMI: ',MUFI,TPMI
!         g%tpmarr=TPMI(i1:i2)
!         g%tprarr=TPRI(i1:i2)
!         write(6,*) 'tprarr: ',g%tprarr
!         g%gk%t0=DTI(1)
!         write(6,*) 't0: ',g%gk%t0, 1./g%gk%u0
!         g%x(1)%data(1)=0d0
! COMPUTATION OF K FROM X
         g%k=calc_nullp(g%gk%q2(1),g%gk%l(1),g%gk%a, &
         g%x%data(2),cos(g%x%data(3)),g%gk%SU(1) &
        *((-1)**g%tprarr),(g%gk%SM(1)*(-1)**g%tpmarr),1,1)
!         write(6,*) 'after geokerr wrap'
         call assign_metric(g%x,transpose(kerr_metric( &
         g%x%data(2),g%x%data(3),g%gk%a)))
         call assign_metric(g%k,transpose(kerr_metric( &
         g%x%data(2),g%x%data(3),g%gk%a)))
 !        write(6,*) 'g%k: ',g%k(20)%data,g%k(20)%metric, &
 !         g%x(20)%data(2:3)
!         write(6,*) 'geokerr_wrapper: ',g%gk%alpha,g%gk%beta,1./UFI,MUFI
         deallocate(ufi); deallocate(mufi)
         deallocate(dti); deallocate(dphi)
         deallocate(lambdai); deallocate(tpmi)
         deallocate(tpri)
         end subroutine geokerr_wrapper
 
         subroutine del_geodesic(g)
         type (geo), intent(inout) :: g
         deallocate(g%x)
         deallocate(g%k)
         deallocate(g%lambda)
         deallocate(g%tpmarr)
         deallocate(g%tprarr)
!         write(6,*) 'del geokerr args',allocated(g%gk%muf)
         call del_geokerr_args(g%gk)
!         write(6,*) 'after del geokerr args'
         end subroutine del_geodesic

         subroutine del_geokerr_args(gk)
         type (geokerr_args), intent(inout) :: gk
!         write(6,*) 'del ge
         deallocate(gk%uf)
         deallocate(gk%muf)
         deallocate(gk%sm)
         deallocate(gk%su)
         deallocate(gk%tpr)
         deallocate(gk%tpm) 
         deallocate(gk%alpha)
         deallocate(gk%beta)
         deallocate(gk%q2)
         deallocate(gk%l)
         deallocate(gk%t0)
         end subroutine del_geokerr_args

       end module geodesics
