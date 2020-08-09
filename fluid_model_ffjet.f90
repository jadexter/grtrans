   module fluid_model_ffjet

      use class_four_vector
      use interpolate, only: interp
      use kerr, only: kerr_metric, lnrf_frame
      implicit none

      namelist /ffjet/  dfile

      character(len=40) :: dfile
      real, dimension(:), allocatable :: rc_arr, thc_arr, phc_arr
      real, dimension(:), allocatable :: rho_arr
      real, dimension(:), allocatable :: p_arr
      real, dimension(:), allocatable :: u0_arr
      real, dimension(:), allocatable :: vr_arr
      real, dimension(:), allocatable :: vth_arr
      real, dimension(:), allocatable :: vph_arr
      real, dimension(:), allocatable :: b0_arr
      real, dimension(:), allocatable :: br_arr
      real, dimension(:), allocatable :: bth_arr
      real, dimension(:), allocatable :: bph_arr
      real, dimension(:), allocatable :: aux1_arr, aux2_arr

      interface init_ffjet_data
        module procedure init_ffjet_data
      end interface

      interface del_ffjet_data
        module procedure del_ffjet_data
      end interface

      interface initialize_ffjet_model
        module procedure initialize_ffjet_model
      end interface
 
      interface ffjet_vals
        module procedure ffjet_vals
      end interface

      contains

        subroutine ffjet_vals(x0,a,rho,p,b,u,bmag)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real(kind=8), dimension(size(x0)) :: done
        real, dimension(size(x0)) :: x2,x1,zm,theta,fone, &
         vpl0,vrl0,vtl0,p0,rho0,bph0,bth0,rd,td,zr,dzero,vr0,vth0,vph0,bth,dummy
        real, dimension(int(sqrt(real(size(rc_arr))))) :: uniqx1 &
          ,uniqx2,uniqr,uniqth
        real, dimension(size(x0),4) :: ppi,rhoi,vrli,vtli, &
         vpli,bri,bthi,bphi,b0i,u0i
        integer, dimension(size(x0)*4) :: indx
        integer, dimension(size(x0)) :: lx1,lx2,ux1,ux2,x1l,x1u,x2l,x2u, &
         umax,one
        integer :: npts,nx1,nx2
        real, dimension(size(x0)), intent(out) :: rho,p,bmag
        type (four_Vector), intent(out), dimension(size(x0)) :: u,b
        ! Computes properties of jet solution from Broderick & Loeb (2009)
        ! JAD 4/23/2010, fortran 3/30/2011
!        write(6,*) 'ffjet: ',size(x0)
        nx1=int(sqrt(real(size(rc_arr)))); nx2=nx1; umax=nx2-1; one=1
        dzero=0d0; done=1d0; fone=1.
        ! Jet solution is symmetric about xy plane
!        write(6,*) 'sz: ',size(rc_arr(1:nx1)),nx1,size(thc_arr),nx2,size(uniqr)
!        write(6,*) 'sz2: ',rc_arr(1)
        uniqr=rc_arr(1:nx1)
!        write(6,*) 'between'
        uniqx2=thc_arr(1:(nx1-1)*nx2+1:nx1)
!        write(6,*) 'between', size(uniqth), size(uniqx2), &
!             size(uniqr), size(uniqx1)
        uniqth=uniqx2 ; uniqx1=log(uniqr)
!        write(6,*) 'uniq: ',minval(rc_arr),maxval(rc_arr), &
!        minval(uniqx1),maxval(uniqx1)
!        write(6,*) 'uniq: ',minval(uniqx2),maxval(uniqx2)
        npts=size(x0)
        theta=x0%data(3)
        zm=cos(theta)
        x2=acos(abs(zm)); theta=x2
        zr=x0%data(2)
!        write(6,*) 'zr: ',zr, zm
        x1=log(zr)
        lx1=int((x1-uniqx1(1))/(uniqx1(nx1)-uniqx1(1))*(nx1-1))+1
        ux1=lx1+1
        lx2=int((x2-uniqx2(1))/(uniqx2(nx2)-uniqx2(1))*(nx2-1))+1
        !write(6,*) 'lx1: ',minval(lx1),maxval(lx1),uniqx1(1)
        !write(6,*) 'lx2: ',minval(lx2),maxval(lx2),size(uniqth)
        lx2=merge(merge(lx2,one,lx2.ge.1),umax,lx2.le.(nx2-1))
        ux2=lx2+1
        td=(theta-uniqth(lx2))/(uniqth(ux2)-uniqth(lx2))
        rd=(zr-uniqr(lx1))/(uniqr(ux1)-uniqr(lx1))
!  !      write(6,*) 'weights: ',minval(lx1),minval(lx2),maxval(lx1), &
!         maxval(lx2),minval(td),maxval(td),minval(rd),maxval(rd)
!  !      write(6,*) 'th(lx2): ',zr(37),uniqr(lx1(37)),uniqr(ux1(37))
        !  r is fastest changing index
!        write(6,*) 'rdtd: ',rd(295),td(295),zr(295),zm(295),lx1(295),lx2(295)
        x1l=lx1-1 ; x1u=ux1-1
        x2l=(lx2-1)*nx1 ; x2u=(ux2-1)*nx1
        indx=(/(/x1l+x2l+1/),(/x1l+x2u+1/),(/x1u+x2l+1/),(/x1u+x2u+1/)/)
        !write(6,*) 'indx: ',minval(indx),maxval(indx),size(indx), &
!         size(rho_arr)
!  !      write(6,*) 'indx: ',rc_arr(x1l(36)+x2l(36)+1),zr(36), &
!         rc_arr(x1u(36)+x2l(36)+1),rd(36)
!  !      write(6,*) 'indx: ',thc_arr(x2l(17)+x1l(17)+1),theta(17), &
!         thc_arr(x2u(17)+x1l(17)+1),td(17)
!  !      write(6,*) 'indx: ',indx(51:100), npts
!  !      write(6,*) 'rhoindx: ',rho_arr(indx)
!        write(6,*) 'reshape'
        rhoi=reshape(rho_arr(indx),(/npts,4/))
        vrli=reshape(vr_arr(indx),(/npts,4/))
        vtli=reshape(vth_arr(indx),(/npts,4/))
        vpli=reshape(vph_arr(indx),(/npts,4/))
        ppi=reshape(p_arr(indx),(/npts,4/))
        b0i=reshape(b0_arr(indx),(/npts,4/))
        bri=reshape(br_arr(indx),(/npts,4/))
        bthi=reshape(bth_arr(indx),(/npts,4/))
        bphi=reshape(bph_arr(indx),(/npts,4/))
        u0i=reshape(u0_arr(indx),(/npts,4/))
!        write(6,*) 'after reshape', vpli(294,:)
!        write(6,*) 'after reshape', vpli(295,:)
!        write(6,*) 'after reshape', vpli(296,:)
        rho=merge(interp(rhoi,rd,td),dzero,x1.gt.uniqx1(1))
        p=merge(interp(ppi,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'rho: ', rho, p
        vrl0=merge(interp(vrli,rd,td),dzero,x1.gt.uniqx1(1))
        vtl0=merge(interp(vtli,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'v'
        b%data(1)=merge(interp(b0i,rd,td),fone,x1.gt.uniqx1(1))
        b%data(2)=merge(interp(bri,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'b'
        b%data(3)=merge(interp(bthi,rd,td),fone,x1.gt.uniqx1(1))
        b%data(4)=merge(interp(bphi,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'u'
        u%data(1)=merge(dble(interp(u0i,rd,td)),done,x1.gt.uniqx1(1))
        vpl0=merge(interp(vpli,rd,td),dzero,x1.gt.uniqx1(1))
  !      write(6,*) 'after assign'
        ! Protect azimuthal velocities at poles
        ! Compute magnitude of interpolated b-field and force b^2 > 0 (need to look into numerical issues here):
        call assign_metric(b,transpose(kerr_metric(zr, &
         real(x0%data(3)),a)))
        bmag=b*b; bmag=merge(bmag,dzero,bmag.ge.0d0)
        bmag=sqrt(bmag)
!        write(6,*) 'maxr: ',maxval(uniqr),maxval(zr)
!        write(6,*) 'bmag: ',maxval(bmag),minval(bmag)
!        write(6,*) 'bmag: ',bmag
!        bmag0=sqrt(max(((normu(zfb0,zfb0,((zr),(acos(zm))),a,/bl)),(fltarr(size(zm)))),dim=2))
        ! Correct \theta, \phi components for reflecting sol'n
        ! assume reflection is that \hat{z} stays same for field (flux),
        ! then assume that \hat{\phi} stays unchanged (L_z)
        ! flips for velocity. so need to flip v^th and b^r. 
!        vtl0=sign(fone,zm)*vtl0
!        b%data(3)=sign(fone,zm)*b%data(3)
!        b%data(2)=sign(fone,zm)*b%data(2)
!         vpl0=sign(fone,zm)*vpl0
!        b%data(4)=sign(fone,zm)*b%data(4)
        ! Get four-velocity:
        call lnrf_frame(vrl0,vtl0,vpl0,zr,a,real(x0%data(3)) &
         ,vr0,vth0,vph0,1)
!        call lnrf_frame(vr_arr,vth_arr,vph_arr,rc_arr,a,thc_arr, &
             
       !  write(6,*) 'lnrf', size(u%data(1)), size(vr0)
        u%data(2)=u%data(1)*dble(vr0)
        u%data(3)=u%data(1)*dble(vth0)
        u%data(4)=u%data(1)*dble(vph0)
!        write(6,*) 'u'
        call assign_metric(u,transpose(kerr_metric(zr,real(x0%data(3)) &
        ,a)))
!        write(6,*) 'leaving ffjet vals',bmag
        !write(6,*) 'udotu: ',maxval(abs(u*u+1d0))
 ! Cut off emission from below eq. plane:
!        rho=merge(rho,rho*0.,zm.ge.0)
        end subroutine ffjet_vals

        subroutine read_ffjet_inputs(ifile)
        character(len=20), intent(in) :: ifile
        open(unit=8,file=ifile,form='formatted',status='old')
        read(8,nml=ffjet)
!        write(6,*) 'read: ',dfile
        close(unit=8)
        end subroutine read_ffjet_inputs

        subroutine initialize_ffjet_model(a,ifile,df)
        real(kind=8), intent(in) :: a
        real :: aa
        integer :: nx, status
        real, dimension(:), allocatable :: b
        character(len=20), intent(in), optional :: ifile
        character(len=40), intent(in), optional :: df
        character(len=20) :: default_ifile='ffjet.in'
        if (present(df)) then
           dfile = df
        else
           if (present(ifile)) then
              call read_ffjet_inputs(ifile)
           else
              call read_ffjet_inputs(default_ifile)
           endif
        endif
!        dfile='m87bl09rfp10xi5a998fluidvars.bin'
       write(6,*) 'init ffjet', dfile
        open(unit=8,file=dfile,form='unformatted',status='old')
        read(8) aa,nx
        call init_ffjet_data(nx,nx,1,nx); allocate(b(nx))
        read(8) rc_arr, thc_arr, rho_arr
        read(8) b, b0_arr, br_arr, bth_arr, bph_arr
        read(8) u0_arr, vr_arr, vth_arr, vph_arr
        read(8) aux1_arr, aux2_arr
        close(8)
 !       write(6,*) 'aa: ',aa,nx,size(rho_arr)
        deallocate(b)
        end subroutine initialize_ffjet_model
 
        subroutine init_ffjet_data(nx,n,aux,auxn,td)
        integer, intent(in) :: n,nx
        integer, intent(in), optional :: aux,td,auxn
        allocate(rho_arr(n)); allocate(p_arr(n)); allocate(u0_arr(n))
        allocate(vr_arr(n)); allocate(vth_arr(n)); allocate(vph_arr(n))
        allocate(b0_arr(n)); allocate(br_arr(n)); allocate(bth_arr(n))
        allocate(bph_arr(n))
        allocate(rc_arr(nx)); allocate(thc_arr(nx))
        write(6,*) 'init ffjet data: ',nx,n,aux,auxn
        if(present(td)) then
          allocate(phc_arr(nx))
        endif
        if(present(aux)) then
          allocate(aux1_arr(n)); allocate(aux2_arr(n))
        endif
        end subroutine init_ffjet_data

        subroutine del_ffjet_data()
        deallocate(rho_arr); deallocate(p_arr); deallocate(u0_arr)
        deallocate(vr_arr); deallocate(vth_arr); deallocate(vph_arr)
        deallocate(b0_arr); deallocate(br_arr); deallocate(bth_arr)
        deallocate(bph_arr); deallocate(rc_arr); deallocate(thc_arr)
        if(allocated(phc_arr)) deallocate(phc_arr)
        if(allocated(aux1_arr)) then
          deallocate(aux1_arr); deallocate(aux2_arr)
        endif
        end subroutine del_ffjet_data

      end module fluid_model_ffjet
