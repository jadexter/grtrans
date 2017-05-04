! wrappers for geokerr routines to allow for f2py use of code

module class_geokerr
  implicit none
   
    integer :: ncase
    real(kind=8) :: offset,u0
    real(kind=8), dimension(:), allocatable :: u,mu,dt,dphi,lam
    integer, dimension(:), allocatable :: tpmi,tpri
    integer, dimension(:), allocatable :: tpmarr,tprarr
    real(kind=8), dimension(:), allocatable :: aarr,barr, &
         ufarr,mufarr,q2arr,larr,smarr,suarr

!$omp threadprivate(u,mu,dt,dphi,lam,tpmi,tpri)

    interface geokerr
       subroutine geokerr(u0,uf,uout,mu0,muf,a,l,q2,alpha, &
            beta,tpm,tpr,su,sm,nup,offset,phit,usegeor,mufill, &
            ncase,kext,next,ufi,mufi,dti,dphi,tpmi,tpri,lambdai)
         real(kind=8), intent(in) :: u0,uf,uout,mu0,muf,a,l,q2,alpha, &
              beta,su,sm,offset
         integer, intent(in) :: kext,next,nup,tpm,tpr,phit,usegeor,mufill
         real(kind=8), intent(out), dimension(nup+next-2*kext) :: ufi,mufi, &
              dti,dphi,lambdai
         integer, intent(out), dimension(nup+next-2*kext) :: tpmi,tpri
         integer, intent(out) :: ncase
       end subroutine geokerr
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

    interface call_geokerr
       module procedure call_geokerr
    end interface

    interface init_pix_data
       module procedure init_pix_data
    end interface

    interface del_pix_data
       module procedure del_pix_data
    end interface

    interface init_geokerr_data
       module procedure init_geokerr_data
    end interface

    interface del_geokerr_data
       module procedure del_geokerr_data
    end interface del_geokerr_data

    contains

      subroutine call_geokerr(u0,uf,uout,mu0,muf,a,l,q2,alpha, &
           beta,tpm,tpr,su,sm,nup,offset,phit,usegeor,mufill, &
           kext,next)
        real(kind=8), intent(in) :: u0,uf,uout,mu0,muf,a,l,q2,alpha, &
              beta,su,sm,offset
        integer, intent(in) :: kext,next,nup,tpm,tpr,phit,usegeor,mufill
        real(kind=8), dimension(nup+next-2*kext) :: ufi,mufi, &
              dti,dphii,lambdai
        real(kind=8) :: tstart,tend
        integer, dimension(nup+next-2*kext) :: tpmii,tprii
!        write(6,*) 'calling geokerr inputs: ',nup+next-2*kext
!        call cpu_time(tstart)
        call geokerr(u0,uf,uout,mu0,muf,a,l,q2,alpha,beta, &
             tpm,tpr,su,sm,nup,offset,phit,usegeor,mufill,ncase, &
             kext,next,u,mu,dt,dphi,tpmi,tpri,lam)
!        call cpu_time(tend)
!        write(6,*) 'geokerr time elapsed: ',tend-tstart
!        write(6,*) 'after geokerr', ncase
!        write(6,*) 'sizes: ',size(ufi), size(u), size(mufi), size(mu)
!        write(6,*) 'sizes: ',size(dti), size(dt), size(lambdai), size(lam)
!        write(6,*) 'sizes: ',size(dphi),size(dphii),size(tpri),size(tprii)
!        write(6,*) 'sizes: ',size(tpmi), size(tpmii)
!        u=ufi; mu=mufi; dphi=dphii
!        write(6,*) 'past u mu phi'
!        dt=dti; lam=lambdai
!        write(6,*) 'past dt lam'
!        tpmi=tpmii; tpri=tprii
!        write(6,*) 'past tpmi tpmri'
      end subroutine call_geokerr

      subroutine  init_pix_data(npix)
        integer, intent(in) :: npix
        allocate(aarr(npix)); allocate(barr(npix)); allocate(q2arr(npix))
        allocate(larr(npix)); allocate(ufarr(npix)); allocate(mufarr(npix))
        allocate(tprarr(npix)); allocate(tpmarr(npix)); allocate(smarr(npix))
        allocate(suarr(npix))
      end subroutine init_pix_data

      subroutine del_pix_data()
        deallocate(aarr); deallocate(barr); deallocate(q2arr)
        deallocate(larr); deallocate(ufarr); deallocate(mufarr)
        deallocate(tprarr); deallocate(tpmarr); deallocate(smarr)
        deallocate(suarr)
      end subroutine del_pix_data


      subroutine init_geokerr_data(npts)
        integer, intent(in) :: npts
        allocate(u(npts)); allocate(mu(npts)); allocate(dt(npts))
        allocate(dphi(npts)); allocate(lam(npts)); allocate(tpmi(npts))
        allocate(tpri(npts))
      end subroutine init_geokerr_data

      subroutine del_geokerr_data()
        deallocate(u); deallocate(mu); deallocate(dt)
        deallocate(dphi); deallocate(lam); deallocate(tpmi)
        deallocate(tpri)
      end subroutine del_geokerr_data

      subroutine pixel(standard,a1,a2,b1,b2,rcut,nrotype,nro,nphi,nup, &
           uout,mu0,a)
         integer, intent(in) :: nro,nphi,nup,nrotype,standard
         real(kind=8), intent(in) :: a1,a2,b1,b2,rcut,a,uout
         real(kind=8), intent(inout) :: mu0
         real(kind=8), dimension(nro*nphi) :: alarr, &
              bearr,uf,muf, &
              q2,l,sm,su
         integer, dimension(nro*nphi) :: tpm,tpr

         call get_pixel_locations(standard,a1,a2,b1,b2,rcut, &
            nrotype,nro,nphi,nup,u0,uout,mu0,a &
            ,offset,ufarr,mufarr,aarr,barr,q2arr,larr, &
            smarr,suarr,tpmarr,tprarr)
!         aarr=alarr; barr=bearr; q2arr=q2; larr=l
!         ufarr=uf; mufarr=muf; tpmarr=tpm; tprarr=tpr
       end subroutine pixel

    end module class_geokerr
