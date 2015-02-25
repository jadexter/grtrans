   module fluid_model_numdisk

!      use class_four_vector
      use interpolate, only: interp
      use fluid_model_thindisk, only: thindisk_vals
!      use kerr, only: kerr_metric, lnrf_frame
      implicit none

      namelist /numdisk/  dfile,tscl,rscl

      character(len=40) :: dfile
      integer :: nr, nphi
      real :: tscl,rscl
      real, dimension(:), allocatable :: rc_arr, phc_arr
      real, dimension(:), allocatable :: T_arr
!      real, dimension(:), allocatable :: p_arr
!      real, dimension(:), allocatable :: u0_arr
!      real, dimension(:), allocatable :: vr_arr
!      real, dimension(:), allocatable :: vth_arr
!      real, dimension(:), allocatable :: vph_arr
!      real, dimension(:), allocatable :: b0_arr
!      real, dimension(:), allocatable :: br_arr
!      real, dimension(:), allocatable :: bth_arr
!      real, dimension(:), allocatable :: bph_arr
!      real, dimension(:), allocatable :: aux1_arr, aux2_arr

      interface init_numdisk_data
        module procedure init_numdisk_data
      end interface

      interface del_numdisk_data
        module procedure del_numdisk_data
      end interface

      interface initialize_numdisk_model
        module procedure initialize_numdisk_model
      end interface
 
      interface numdisk_vals
        module procedure numdisk_vals
      end interface

      contains

        subroutine numdisk_vals(r,phi,a,T,omega)
!        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real, intent(in), dimension(:) :: r,phi
!        real(kind=8), dimension(size(r)) :: done
        real, dimension(size(r)) :: x2,x1,fone, &
         rd,pd,zr,dzero,dummy
        real, dimension(int(sqrt(real(size(rc_arr))))) :: uniqx1 &
          ,uniqx2,uniqr,uniqp,omega
        real, dimension(size(r),4) :: Ti
        integer, dimension(size(r)*4) :: indx
        integer, dimension(size(r)) :: lx1,lx2,ux1,ux2,x1l,x1u,x2l,x2u, &
         umax,one
        integer :: npts,nx1,nx2
        real :: dph
        real, dimension(size(r)), intent(out) :: T
!        type (four_Vector), intent(out), dimension(size(r)) :: u
        ! Load numerical 2D disk model (r,phi)
        ! JAD 4/24/2012
!        write(6,*) 'numdisk: ',size(r)
        x1=log(r); x2=phi
!        write(6,*) 'numdisk vals coords: ',x1,x2
        nx1=nr; nx2=nphi; umax=1; one=1
        dzero=0d0; fone=1.
!        write(6,*) 'sz: ',size(rc_arr(1:nx1)),nx1,size(thc_arr),nx2,size(uniqr)
!        write(6,*) 'sz2: ',rc_arr(1)
        uniqr=rc_arr(1:nx1)
!        write(6,*) 'between'
        uniqx2=phc_arr(1:(nx1-1)*nx2+1:nx1)
!        write(6,*) 'uniq: ',uniqr, uniqx2
!        write(6,*) 'between', size(uniqp), size(uniqx2), &
!             size(uniqr), size(uniqx1)
        uniqp=uniqx2 ; uniqx1=log(uniqr)
!        write(6,*) 'uniq: ',minval(rc_arr),maxval(rc_arr), &
!        minval(uniqx1),maxval(uniqx1)
!        write(6,*) 'uniq: ',minval(uniqx2),maxval(uniqx2)
!        x2=x0%data(4); phi=x2
!        zm=cos(phi)
!        x2=acos(abs(zm)); phi=x2
!        zr=x0%data(2)
!        write(6,*) 'zr: ',zr, zm
        lx1=int((x1-uniqx1(1))/(uniqx1(nx1)-uniqx1(1))*(nx1-1))+1
        ux1=lx1+1
        lx2=int((x2-uniqx2(1))/(uniqx2(nx2)-uniqx2(1))*(nx2-1))+1
!        write(6,*) 'lx1: ',x1,minval(lx1),maxval(lx1),uniqx1(1)
!        write(6,*) 'lx2: ',x2,zm,minval(lx2),maxval(lx2),size(uniqp)
        lx2=merge(merge(lx2,one,lx2.ge.1),umax,lx2.le.(nx2-1))
        ux2=lx2+1
        dph=uniqp(2)-uniqp(1)
        pd=(phi-uniqp(lx2))/dph
        rd=(zr-uniqr(lx1))/(uniqr(ux1)-uniqr(lx1))
!  !      write(6,*) 'weights: ',minval(lx1),minval(lx2),maxval(lx1), &
!         maxval(lx2),minval(pd),maxval(pd),minval(rd),maxval(rd)
!  !      write(6,*) 'th(lx2): ',zr(37),uniqr(lx1(37)),uniqr(ux1(37))
        !  r is fastest changing index
        x1l=lx1-1 ; x1u=ux1-1
        x2l=(lx2-1)*nx1 ; x2u=(ux2-1)*nx1
!        write(6,*) 'indx: ',lx1, lx2
        indx=(/(/x1l+x2l+1/),(/x1l+x2u+1/),(/x1u+x2l+1/),(/x1u+x2u+1/)/)
!  !      write(6,*) 'indx: ',minval(indx),maxval(indx),size(indx), &
!         size(rho_arr)
!  !      write(6,*) 'indx: ',rc_arr(x1l(36)+x2l(36)+1),zr(36), &
!         rc_arr(x1u(36)+x2l(36)+1),rd(36)
!  !      write(6,*) 'indx: ',thc_arr(x2l(17)+x1l(17)+1),phi(17), &
!         thc_arr(x2u(17)+x1l(17)+1),pd(17)
!  !      write(6,*) 'indx: ',indx(51:100), npts
!  !      write(6,*) 'rhoindx: ',rho_arr(indx)
!        write(6,*) 'reshape'
        npts=size(r)
        ti=reshape(t_arr(indx),(/npts,4/))
!        vrli=reshape(vr_arr(indx),(/npts,4/))
!        vtli=reshape(vth_arr(indx),(/npts,4/))
!        vpli=reshape(vph_arr(indx),(/npts,4/))
!        ppi=reshape(p_arr(indx),(/npts,4/))
!        b0i=reshape(b0_arr(indx),(/npts,4/))
!        bri=reshape(br_arr(indx),(/npts,4/))
!        bthi=reshape(bth_arr(indx),(/npts,4/))
!        bphi=reshape(bph_arr(indx),(/npts,4/))
!        u0i=reshape(u0_arr(indx),(/npts,4/))
!        write(6,*) 'after reshape'
!        rho=merge(interp(rhoi,rd,pd),dzero,x1.gt.uniqx1(1))
        T=merge(ti(:,1),fone,x1.gt.uniqx1(1))
!        write(6,*) 'rho: ', rho, p
!        vrl0=merge(interp(vrli,rd,pd),dzero,x1.gt.uniqx1(1))
!        vtl0=merge(interp(vtli,rd,pd),fone,x1.gt.uniqx1(1))
!        write(6,*) 'v'
!        b%data(1)=merge(interp(b0i,rd,pd),fone,x1.gt.uniqx1(1))
!        b%data(2)=merge(interp(bri,rd,pd),fone,x1.gt.uniqx1(1))
!        write(6,*) 'b'
!        b%data(3)=merge(interp(bthi,rd,pd),fone,x1.gt.uniqx1(1))
!        b%data(4)=merge(interp(bphi,rd,pd),fone,x1.gt.uniqx1(1))
!        write(6,*) 'u'
!        u%data(1)=merge(dble(interp(u0i,rd,pd)),done,x1.gt.uniqx1(1))
!        vpl0=merge(interp(vpli,rd,pd),dzero,x1.gt.uniqx1(1))
  !      write(6,*) 'after assign'
        ! Protect azimuthal velocities at poles
        ! Compute magnitude of interpolated b-field and force b^2 > 0 (need to look into numerical issues here):
!        call assign_metric(b,transpose(kerr_metric(zr, &
!         real(x0%data(3)),a)))
!        bmag=b*b; bmag=merge(bmag,dzero,bmag.ge.0d0)
!        bmag=sqrt(bmag)
!        write(6,*) 'maxr: ',maxval(uniqr),maxval(zr)
!        write(6,*) 'bmag: ',maxval(bmag),minval(bmag)
!        write(6,*) 'bmag: ',bmag
!        bmag0=sqrt(max(((normu(zfb0,zfb0,((zr),(acos(zm))),a,/bl)),(fltarr(size(zm)))),dim=2))
        ! Correct \phi, \phi components for reflecting sol'n:
!  !      write(6,*) 'before sign: ',vtl0
!        vtl0=sign(vtl0,zm) ; vpl0=sign(vpl0,zm)
!        b%data(3)=sign(b%data(3),dble(zm))
!        b%data(4)=sign(b%data(4),dble(zm))
        ! Get four-velocity:
!        call lnrf_frame(vrl0,vtl0,vpl0,zr,a,real(x0%data(3)) &
!         ,vr0,vth0,vph0,1)
!         write(6,*) 'lnrf', size(u%data(1)), size(vr0)
!        u%data(2)=u%data(1)*dble(vr0)
!        u%data(3)=u%data(1)*dble(vth0)
!        u%data(4)=u%data(1)*dble(vph0)
!        write(6,*) 'u'
        call thindisk_vals(r,dzero+acos(-fone)/2.,a,dummy,omega)
!        call assign_metric(u,transpose(kerr_metric(zr,real(x0%data(3)) &
!        ,a)))
!        write(6,*) 'leaving numdisk vals',bmag
        end subroutine numdisk_vals

        subroutine read_numdisk_inputs(ifile)
        character(len=20), intent(in) :: ifile
        open(unit=8,file=ifile,form='formatted',status='old')
        read(8,nml=numdisk)
!        write(6,*) 'read: ',dfile
        close(unit=8)
        end subroutine read_numdisk_inputs

        subroutine initialize_numdisk_model(ifile,df,ts,rs)
!        real(kind=8), intent(in) :: a
        real :: aa
        real, dimension(:), allocatable :: b
        character(len=20), intent(in), optional :: ifile
        character(len=20) :: default_ifile='numdisk.in', status
        real, intent(in), optional :: ts,rs
        character(len=40), intent(in), optional :: df
        if (present(df)) then
           dfile = df
           tscl = ts
           rscl = rs
        else
           if (present(ifile)) then
              call read_numdisk_inputs(ifile)
           else
              call read_numdisk_inputs(default_ifile)
           endif
        endif
!        write(6,*) 'init numdisk ', dfile
        open(unit=8,file=dfile,form='unformatted',status='old')
        write(6,*) 'file open'
        read(8) nr
        read(8) nphi
        write(6,*) 'nx: ',nr,nphi,nr*nphi,kind(nr)
        call init_numdisk_data()
        read(8) rc_arr, phc_arr, T_arr
        write(6,*) 'rc: ',minval(rc_arr),maxval(rc_arr)
        write(6,*) 'phc: ',minval(phc_arr),maxval(phc_arr)
        write(6,*) 't: ',minval(t_arr),maxval(t_arr)
        close(8)
! Use scalings for T0, r0 from numdisk.in:
        T_arr=T_arr*tscl
        rc_arr=rc_arr*rscl
 !       write(6,*) 'aa: ',aa,nx,size(rho_arr)
        end subroutine initialize_numdisk_model
 
        subroutine init_numdisk_data()
!        integer, intent(in), optional :: aux,pd,auxn
!        write(6,*) 'init numdisk data: '
        allocate(T_arr(nr*nphi))
        allocate(rc_arr(nr*nphi)); allocate(phc_arr(nr*nphi))
!        write(6,*) 'init numdisk data: ',nx,n,aux,auxn
!        if(present(pd)) then
!          allocate(phc_arr(nx))
!        endif
!        if(present(aux)) then
!          allocate(aux1_arr(n)); allocate(aux2_arr(n))
!        endif
        end subroutine init_numdisk_data

        subroutine del_numdisk_data()
        deallocate(T_arr)
        deallocate(rc_arr); deallocate(phc_arr)
 !       if(allocated(phc_arr)) deallocate(phc_arr)
 !       if(allocated(aux1_arr)) then
 !         deallocate(aux1_arr); deallocate(aux2_arr)
 !       endif
        end subroutine del_numdisk_data

      end module fluid_model_numdisk
