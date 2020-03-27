   module fluid_model_iharm

      use class_four_vector
      use interpolate, only: interp
      use kerr, only: kerr_metric, lnrf_frame, uks2ubl, bl2ks
      use phys_constants, only: pi
      use math, only: zbrent
      implicit none

      namelist /harm/  dfile, hfile, nt, indf, gfile

      character(len=100) :: dfile, hfile, gfile
      integer :: n, ndumps, nt, indf, dlen, nhead
      integer :: ndim=3, IS_MMKS
      integer, dimension(:), allocatable :: dumps
      real(kind=8) :: tstep=20d0, &
           toffset=0d0,tcur,dx1,dx2, &
           dx3,a,gam,Rin,Rout,hslope,R0,ti,tj,tk,startx1, startx2,rdump_cnt,rdump01_cnt, &
           startx3,game,lim,game4, game5, tf, &
           tavgstart,SMALL=1e-20, poly_xt, poly_alpha, mks_smooth, Lunit, Munit
      integer :: SDUMP,eHEAT,eCOND,nx1,nx2,nx3,metric,n_prim,eRAD,BL
      real :: asim
      real, dimension(:), allocatable :: t
      real, dimension(:), allocatable :: x1_arr, x2_arr, x3_arr, r_arr, th_arr, ph_arr
      real, dimension(:), allocatable :: rho_arr, p_arr, u0_arr, vrl_arr, &
        vpl_arr, vtl_arr, temp_arr, kela_arr
! global variables for BL3 but should have these read in from header files rather than hard-coded
      real, dimension(:), allocatable :: b0_arr, br_arr, bth_arr, bph_arr, gdet
      real, dimension(:,:), allocatable :: drdx,metric_cov,metric_con

      interface init_iharm_data
        module procedure init_iharm_data
      end interface

      interface del_iharm_data
        module procedure del_iharm_data
      end interface

      interface initialize_iharm_model
        module procedure initialize_iharm_model
      end interface

      interface transformbl2mmks
         module procedure transformbl2mmks
      end interface

      interface ummks2uks
         module procedure ummks2uks
      end interface

      interface iharm_vals
        module procedure iharm_vals
      end interface

!      interface calcrmks
!         module procedure calcrmks
!         module procedure calcrmks_single
!      end interface

!      interface findx1mks
!         module procedure findx1mks
!         module procedure findx1mks_array
!      end interface

      interface findx2iharm
         module procedure findx2iharm
         module procedure findx2iharm_array
      end interface

      contains

        function findx2iharm(x2,args) result(diff)
        real(kind=8), intent(in) :: x2
        real(kind=8), intent(in), dimension(:) :: args
        real(kind=8) :: h, theta, diff
!        real(kind=8), dimension(size(x2)) :: thetacur
!        h=args(2); theta=args(1)
!        thetacur=calcthmksh(x2)
!        diff=theta-thetacur
        end function findx2iharm

        function calcthmksh(x2) result(th)
          real(kind=8), intent(in), dimension(:) :: x2
          real(kind=8), dimension(size(x2)) :: th
!          th=pi/2.*(1.+x2)+((1.-hslope)/2.)*sin(pi*(1.+x2))
          th=pi*x2+(1.-hslope)/2.*sin(2.*pi*x2)
        end function calcthmksh

        function findx2iharm_array(x2,args) result(diff)
        real(kind=8), intent(in), dimension(:) :: x2
        real(kind=8), intent(in), dimension(:,:) :: args
        real(kind=8) :: h
        real(kind=8), dimension(size(x2)) :: theta, diff
        h=args(1,2); theta=args(:,1)
        diff=theta-(pi/2d0*(1.+x2)+((1d0-h)/2d0)*sin(pi*(1d0+x2)))
        end function findx2iharm_array

        function  findx2mksh(x2,args) result(diff)
        real(kind=8), intent(in), dimension(:) :: x2
        real(kind=8), intent(in), dimension(:,:) :: args
        real(kind=8) :: h
        real(kind=8), dimension(size(x2)) :: theta, diff, thetacur
        theta=args(:,1)
        thetacur=calcthmksh(x2)
        diff=theta-thetacur
        end function findx2mksh

        function calcthmmks(x2,x1) result(theta)
          real, intent(in), dimension(:) :: x2, x1
          real, dimension(size(x2)) :: theta,thetag,thetaj,D
          real :: A,B,C
          A=mks_smooth
          B=poly_xt
          C=poly_alpha
          D=pi/(2.+2./(B**C*(1+C)))
          thetag=pi*x2+(1.-hslope)/2.*sin(2.*pi*x2)
          thetaj=D*(2.*x2-1.)*(1.+((2.*x2-1.)/B)**C/(1.+C))+pi/2.
          theta=thetag+exp(-A*(x1-startx1))*(thetaj-thetag)
        end function calcthmmks

        function findx2mmks(x2,args) result(diff)
          ! Calculates \theta for iharm
          ! JAD 5/24/2019
          real(kind=8), intent(in), dimension(:) :: x2
          real(kind=8), intent(in), dimension(:,:) :: args
          real(kind=8), dimension(size(x2)) :: diff,th,theta,x1
          th=args(:,2) ; x1=args(:,1)
          theta=calcthmmks(real(x2),real(x1))
          diff=th-theta
        end function findx2mmks

        subroutine transformbl2mmks(r,th,ph,x1,x2,x3)
        ! transform Boyer-Lindquist coordinates to modified Kerr-Schild coordinates used by thickdisk
        real, intent(in), dimension(:) :: r,th,ph
        real, intent(out), dimension(size(r)) :: x1,x2,x3
        real(kind=8), dimension(5) :: args
        real(kind=8), dimension(size(r),2) :: thargs
        real(kind=8), dimension(size(r)) :: mone,one
        integer :: i
        mone(:)=-1d0; one(:)=1d0
        x1=log(r)
        thargs(:,1)=dble(x1)
        thargs(:,2)=dble(th)
        x2=zbrent(findx2mmks,mone,one,thargs,1d-6)
        x3=ph
        end subroutine transformbl2mmks

        subroutine transformbl2mksh(r,th,ph,x1,x2,x3)
        ! transform Boyer-Lindquist coordinates to modified Kerr-Schild coordinates used by thickdisk
        real, intent(in), dimension(:) :: r,th,ph
        real, intent(out), dimension(size(r)) :: x1,x2,x3
        real(kind=8), dimension(5) :: args
        real(kind=8), dimension(size(r),1) :: thargs
        real(kind=8), dimension(size(r)) :: mone,one,zero
        integer :: i
        mone(:)=-1d0; one(:)=1d0; zero=0d0
        x1=log(r)
!        thargs(:,1)=dble(x1)
        thargs(:,1)=dble(th)
! X2 IS ALWAYS [0,1] here
        x2=zbrent(findx2mksh,zero,one,thargs,1d-6)
        x3=ph
        end subroutine transformbl2mksh

        function ummks2uks(fum,x1,x2) result(uks)
          real(kind=8), intent(in), dimension(:) :: x1,x2
          type (four_vector), dimension(:), intent(in) :: fum
!          real(kind=8), intent(in) :: xbr
          type (four_vector), dimension(size(fum)) :: uks
          real(kind=8), dimension(size(x1)) :: r,dx2,dx1,drdx1,dthdx1,dthdx2
          ! Convert four-vectors from MMKS to KS numerically using central differences.
          r=exp(x1)
          dx2=1d-6*x2; dx1=1d-4*x1
          drdx1=r
          dthdx1=(calcthmmks(real(x2),real(x1+.5*dx1))- &
               calcthmmks(real(x2),real(x1-.5*dx1)))/dx1
          dthdx2=(calcthmmks(real(x2+.5*dx2),real(x1))-calcthmmks(real(x2-.5*dx2),real(x1)))/dx2
          uks%data(2)=drdx1*fum%data(2)
          uks%data(3)=dthdx1*fum%data(2)+dthdx2*fum%data(3)
          ! phi doesn't change
          uks%data(4)=fum%data(4)
          ! time component doesn't change
          uks%data(1)=fum%data(1)
        end function ummks2uks

        function umksh2uks(fum,r,x2) result(uks)
          ! Converts a Kerr-Schild spherical 4-velocity to a Boyer-Lindquist one.
          ! JAD 11/10/2008
          ! Do simple transformation first:
          ! to 3D on 19/11/15 JD
          real, dimension(:), intent(in) :: r, x2
          type (four_vector), dimension(:), intent(in) :: fum
          type (four_vector), dimension(size(r)) :: uks
!          real, intent(in), optional :: hin
!          real :: hval
          real(kind=8), dimension(size(r)) :: ur, uth, uph
          real, dimension(size(r)) :: dthdx2
          write(6,*) 'read harm umksh2uks present h'
!          if (present(hin)) then
!             hval=hin
!          else
!             hval=0.3
!          endif
          write(6,*) 'hslope: ',hslope
          ur=dble(r)*fum%data(2)
          ! Compute partial derivatives:
!          dthdx2=pi
!          dthdx2=pi*(1.+(1.-hval)*cos(2.*pi*x2))
!          dthdx2=pi/2.*(1.+(1.-hslope)*cos(pi*(1.+x2)))
          dthdx2=pi*(1.+(1.-hslope)*cos(2*pi*x2))
! new version dth/dx2 = pi * (1 + (1-h) * cos(2 pi x2))
! matches above i'd already commented out...
          uth=fum%data(3)*dble(dthdx2)
          uph=fum%data(4)
          uks%data(1)=fum%data(1)
          uks%data(2)=ur
          uks%data(3)=uth
          uks%data(4)=uph
          write(6,*) 'done with umksh2uks'
        end function umksh2uks

        subroutine iharm_vals(x0,a,rho,p,b,u,bmag,kela,gamma)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real(kind=8), dimension(size(x0)) :: done,pfac,nfac
        real, dimension(size(x0)) :: x3,x2,x1,zm,theta,fone, &
         vpl0,vrl0,vtl0,rd,td, pd,rttd,zr,dzero, &
         vr0,vth0,vph0,bth,dummy,tt,ttd,zt,zphi,zpp
        real, dimension(nx1) :: uniqx1,uniqr
        real, dimension(nx2) :: uniqx2,uniqth
        real, dimension(nx3) :: uniqx3,uniqph
        real, dimension(size(x0),2**(ndim+1)) :: ppi,rhoi,vrli,vtli, &
         vpli,bri,bthi,bphi,b0i,u0i,ri,thi,phii,kelai
!        integer, dimension(size(x0)*(2**ndim)) :: sindx
        real, dimension(size(x0)) :: dth, minph
        integer, dimension(size(x0)*2*(2**ndim)) :: indx
        integer, dimension(size(x0)) :: lx1,lx2, lx3, ux1,ux2,ux3,x1l,x1u,x2l,x2u,x3l,x3u, &
         umax,one,tindx
        integer :: npts,i,maxtindx
        real, dimension(size(x0)), intent(out) :: rho,p,bmag,kela
        real, intent(out) :: gamma
        type (four_Vector), intent(out), dimension(size(x0)) :: u,b
        ! Interpolates HARM data to input coordinates
        ! JAD 3/20/2009, fortran 11/12/2012
!        write(6,*) 'harm: ',size(x0)
        !nx1=floor(sqrt(real(size(x1_arr)))); nx2=nx1
!        write(6,*) 'iharm vals: ',nx1,nx2,nx3,size(x1_arr),size(x2_arr),size(x3_arr)
!        write(6,*) 'iharm vals'
        umax=nx2-1; one=1
        dzero=0d0; done=1d0; fone=1.
        uniqx3=x3_arr(1:nx3)
        uniqx2=x2_arr(1:nx3*(nx2-1)+1:nx3)
        uniqx1=x1_arr(1:nx3*(nx2-1)*nx1+1:nx3*nx2)
        uniqr=exp(uniqx1)
        if(BL.eq.1) uniqth=pi/2.*(1.+uniqx2)+((1.-hslope)/2.)*sin(pi*(1.+uniqx2))
        uniqph=uniqx3
        npts=size(x0)
        theta=x0%data(3)
        zr=x0%data(2)
        zt=x0%data(1)
        tt=bl2ks(dble(zr),dble(zt),dble(a),1)-bl2ks(dble(zr(1)),0d0,dble(a),1)
        tt=-tt
        zpp=x0%data(4)
        zphi=bl2ks(dble(zr),dble(zpp),dble(a))
        zphi=mod(zphi,(2.*pi))
        where(zphi.lt.0.)
            zphi=zphi+2.*pi
        endwhere
!        write(6,*) 'iharm vals transform'
!        if(BL.eq.3) then
        if(IS_MMKS.eq.1) then
           call transformbl2mmks(zr,theta,zphi,x1,x2,x3)
        else
           call transformbl2mksh(zr,theta,zphi,x1,x2,x3)
        endif
!        write(6,*) 'transform r: ',maxval(zr), minval(zr), maxval(x1), minval(x1)
!        write(6,*) 'transform th: ',maxval(theta), minval(theta), minval(x2), maxval(x2)
!        write(6,*) 'transform phi: ',maxval(zphi), minval(zphi), minval(x3), maxval(x3)
!        write(6,*) 'dif x3: ',(uniqx3(nx3)-uniqx3(1)),x3,uniqx3
        lx1=floor((x1-uniqx1(1))/(uniqx1(nx1)-uniqx1(1))*(nx1-1))+1
        ux1=lx1+1
        lx2=floor((x2-uniqx2(1))/(uniqx2(nx2)-uniqx2(1))*(nx2-1))+1
        ux2=lx2+1
        lx3=floor((x3-uniqx3(1))/(uniqx3(nx3)-uniqx3(1))*(nx3-1))+1
        ux3=lx3+1
!        lx2=merge(merge(lx2,one,lx2.ge.1),umax,lx2.le.(nx2-1))
        lx2=merge(lx2,one,lx2.ge.1)
        lx2=merge(lx2,nx2,lx2.le.nx2)
        ux2=merge(ux2,one,ux2.ge.1)
        ux2=merge(ux2,nx2,ux2.le.nx2)
        where(ux3.gt.nx3)
            ux3=1
        endwhere
        where(lx3.lt.1)
            lx3=nx3
            minph=uniqph(lx3)-2.*pi
        endwhere
! use nearest neighbor rather than interpolate
        where(uniqr(lx1).le.(1.+sqrt(1.-a**2.)))
           rd=1.
           pfac=1e-3
           nfac=1e-3
        elsewhere
           pfac=1.
           nfac=1.
        endwhere
!  nearest neighbor
        rd(:)=1.; td(:)=1.; pd(:)=1.
!        write(6,*) 'coords: ',
!        write(6,*) 'rd td pd: ',minval(rd),maxval(rd),minval(td),maxval(td),minval(pd),maxval(pd)
!        write(6,*) 'iharm ux lx: ',minval(lx1),maxval(ux1),minval(lx2),maxval(ux2),minval(lx3),maxval(ux3)
        ! th is fastest changing index
        x3l=lx3-1; x3u=ux3-1
        x2l=(lx2-1)*nx3 ; x2u=(ux2-1)*nx3
        x1l=(lx1-1)*nx2*nx3; x1u=(ux1-1)*nx2*nx3
! Indices for nearest neighbors in 2D
!        write(6,*) 'sindx', size(sindx)
!        sindx=(/(/x1l+x2l+1/),(/x1l+x2u+1/),(/x1u+x2l+1/),(/x1u+x2u+1/)/)
!        write(6,*) 'sindx', maxval(sindx), sindx(npts+1), x1l(1)+x2u(1)+1, size(sindx)
! Now find timestep indices:
        maxtindx=nt-2
        where(floor(tt/tstep).le.maxtindx)
           tindx(:)=floor(tt/tstep)
           ttd=(tt-tindx*tstep)/tstep
        elsewhere
           tindx=maxtindx+1
           ttd=0.
        endwhere
!        write(6,*) 'tindx: ',tindx
!        write(6,*) 'tindx: ',maxtindx,floor(tt/tstep)
! Create index array across times and 2D nearest neighbors (first half at first time, second half at second time):
! Need to make sure there aren't integer problems when these indices get very large
!        write(6,*) 'indx', size(indx), 2**ndim, maxval(tindx), minval(tindx), tstep
        indx=(/(/x1l+x2l+x3l+1+n*tindx/),(/x1l+x2u+x3l+1+n*tindx/),(/x1u+x2l+x3l+1+n*tindx/), &
            (/x1u+x2u+x3l+1+n*tindx/),(/x1l+x2l+x3u+1+n*tindx/), &
            (/x1l+x2u+x3u+1+n*tindx/),(/x1u+x2l+x3u+1+n*tindx/),(/x1u+x2u+x3u+1+n*tindx/), &
            (/x1l+x2l+x3l+1+n*(tindx+1)/),(/x1l+x2u+x3l+1+n*(tindx+1)/),(/x1u+x2l+x3l+1+n*(tindx+1)/), &
            (/x1u+x2u+x3l+1+n*(tindx+1)/),(/x1l+x2l+x3u+1+n*(tindx+1)/),(/x1l+x2u+x3u+1+n*(tindx+1)/), &
            (/x1u+x2l+x3u+1+n*(tindx+1)/),(/x1u+x2u+x3u+1+n*(tindx+1)/)/)
! keep indx from going out of bounds when loading one time slice
        where(indx.gt.size(rho_arr))
           indx=indx-n
        endwhere
!        do i=1,2**ndim
!           indx((i-1)*npts+1:i*npts)=sindx((i-1)*npts+1)+npts*tindx
!           indx((2**ndim+i-1)*npts+1:(i+2**ndim)*npts)=sindx((i-1)*npts+1)+npts*(tindx+1)
!        end do
!        write(6,*) 'after indx', maxval(indx), size(indx),npts,ndim!, maxval(sindx), indx(npts+1), n*maxval(tindx+1),size(rho_arr),maxloc(indx)!size(rho_arr(indx))
!        write(6,*) 'after indx: ',indx,size(x2_arr)
        rhoi=reshape(rho_arr(indx),(/npts,2**(ndim+1)/))
        vrli=reshape(vrl_arr(indx),(/npts,2**(ndim+1)/))
        vtli=reshape(vtl_arr(indx),(/npts,2**(ndim+1)/))
        vpli=reshape(vpl_arr(indx),(/npts,2**(ndim+1)/))
        ppi=reshape(p_arr(indx),(/npts,2**(ndim+1)/))
!        write(6,*) 'after ppi: ',npts,ndim
        b0i=reshape(b0_arr(indx),(/npts,2**(ndim+1)/))
        bri=reshape(br_arr(indx),(/npts,2**(ndim+1)/))
        bthi=reshape(bth_arr(indx),(/npts,2**(ndim+1)/))
        bphi=reshape(bph_arr(indx),(/npts,2**(ndim+1)/))
        u0i=reshape(u0_arr(indx),(/npts,2**(ndim+1)/))
        if(eHEAT.eq.1) then
           kelai=reshape(kela_arr(indx),(/npts,2**(ndim+1)/))
        endif
! coordinates for debugging
        ri=reshape(r_arr(indx),(/npts,2**(ndim+1)/))
        thi=reshape(th_arr(indx),(/npts,2**(ndim+1)/))
        phii=reshape(ph_arr(indx),(/npts,2**(ndim+1)/))
        rttd=0.
        rho=merge(interp(rhoi,rttd,pd,rd,td),dzero,x1.gt.uniqx1(1))*nfac
        p=merge(interp(ppi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))*pfac
        if(eHEAT.eq.1) then
           kela=merge(interp(kelai,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))*pfac
        endif
!        write(6,*) 'rho: ', rho, p
        vrl0=merge(interp(vrli,rttd,pd,rd,td),dzero,x1.gt.uniqx1(1))
        vtl0=merge(interp(vtli,rttd,pd,rd,td),dzero,x1.gt.uniqx1(1))
!        write(6,*) 'v'
        b%data(1)=merge(interp(b0i,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
        b%data(2)=merge(interp(bri,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'b'
        b%data(3)=merge(interp(bthi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
        b%data(4)=merge(interp(bphi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
 !       write(6,*) 'u'
        u%data(1)=merge(dble(interp(u0i,rttd,pd,rd,td)),done,x1.gt.uniqx1(1))
        vpl0=merge(interp(vpli,rttd,pd,rd,td),dzero,x1.gt.uniqx1(1))
!        write(6,*) 'after assign'
        ! Protect azimuthal velocities at poles
        ! Compute magnitude of interpolated b-field and force b^2 > 0 (need to look into numerical issues here):
        call assign_metric(b,transpose(kerr_metric(zr, &
         real(x0%data(3)),a)))
        bmag=b*b; bmag=merge(bmag,dzero,bmag.ge.0d0)
        bmag=sqrt(bmag)
        ! Get four-velocity:
        call lnrf_frame(vrl0,vtl0,vpl0,zr,a,real(x0%data(3)) &
         ,vr0,vth0,vph0,1)
        u%data(2)=u%data(1)*dble(vr0)
        u%data(3)=u%data(1)*dble(vth0)
        u%data(4)=u%data(1)*dble(vph0)
        call assign_metric(u,transpose(kerr_metric(zr,real(x0%data(3)) &
        ,a)))
! ASSIGNING EOS GAMMA as output to have enough information for both p, u
        gamma=gam
!        write(6,*) 'min vals u: ', minval(u%data(1)), maxval(abs(u*u+1))
!        write(6,*) 'min vals p: ',minval(rhoi),minval(ppi)
!        write(6,*) 'harm vals r: ',zr
!        write(6,*) 'harm vals rho: ',rho
!        write(6,*) 'harm vals rd: ',rd
!        write(6,*) 'harm vals ux: ',nx1,ux1
!        write(6,*) 'harm vals udotu: ',(abs(u*u+1))(minloc(zr))
!        write(6,*) 'harm vals minloc: ',minloc(zr),rd(minloc(zr)), &
!             td(minloc(zr)),x0(minloc(zr))%data(3),uniqth(lx2(minloc(zr))), &
!             lx2(minloc(zr)),pd(minloc(zr)),lx1(minloc(zr)),ux1(minloc(zr)), &
!             uniqx1(1),x1(minloc(zr))
!        write(6,*) 'harm vals minloc b0i: ',b0i(minloc(zr),:)
!        write(6,*) 'harm vals minloc bri: ',bri(minloc(zr),:)
!        write(6,*) 'harm vals minloc bthi: ',bthi(minloc(zr),:)
!        write(6,*) 'harm vals minloc bphi: ',bphi(minloc(zr),:)
!        write(6,*) 'harm vals interp bi: ',b(minloc(zr))%data(1),b(minloc(zr))%data(2), &
!             b(minloc(zr))%data(3),b(minloc(zr))%data(4)
!        write(6,*) 'harm vals minloc x0: ', ri(minloc(zr),:),thi(minloc(zr),:),&
!             phii(minloc(zr),:)
!        write(6,*) 'harm coords: ',zr(minloc(p)),theta(minloc(p)),zphi(minloc(p))
!        write(6,*) 'harm vals rhoi: ',rhoi(minloc(p),:)
        end subroutine iharm_vals

! this is unchanged for plain text / binary version as long as you provide the right header file name
        subroutine read_iharm_data_header()
!          integer, intent(in) :: nhead
          integer :: n,nheadmax=200,status
          real(kind=8), dimension(:), allocatable :: header
!          character(len=100) :: metric
          allocate(header(nheadmax)); header(:)=-1d20
          write(6,*) 'read iharm header: ',hfile
          open(unit=8,file=hfile)
          read(8,*,iostat=status) header
          close(unit=8)
          write(6,*) 'iharm header: ',header(1),header(2)
          write(6,*) 'iharm header: ',header(3),header(4)
          write(6,*) 'iharm header: ',header(5),header(6)
! now read variables from header
          tcur=header(1)
          nx1=int(header(2))
          nx2=int(header(3))
          nx3=int(header(4))
          asim=header(5)
          hslope=header(6)
          gam=header(7)
          mks_smooth=header(8)
          poly_xt=header(9)
          poly_alpha=header(10)
          startx1=header(11)
          metric=header(12)
! flag for MMKS: metric=0 for MKS, metric=1 for MMKS (=FMKS)
          IS_MMKS=metric
          eHEAT=header(13)
          eRAD=header(14)
          n_prim=header(15)
          Rin=header(16)
          Rout=header(17)
! for radiation, get Lunit and Munit to compare with M_BH, Mdot
          if(eRAD.eq.1) then
             Lunit=header(18)
             Munit=header(19)
          endif
! TEMPORARY set so that grid information is read from same file
          SDUMP=0
          write(6,*) 'done reading header: ',gam,mks_smooth, &
               poly_xt,poly_alpha,startx1
        end subroutine read_iharm_data_header

        subroutine read_iharm_inputs(ifile)
        character(len=20), intent(in) :: ifile
        open(unit=8,file=ifile,form='formatted',status='old')
        read(8,nml=harm)
        write(6,*) 'read: ',nt
        close(unit=8)
        end subroutine read_iharm_inputs


! construct grid from iharm data startx1,startx2,startx3,dx1,dx2,dx3
        subroutine read_iharm_grid_file(grid_file)
          character(len=100), intent(in) :: grid_file
          character :: header_byte
          integer :: nheader_bytes,n,grid_dlen,i,j,k,indx
! should be x3 changing fastest, x1 slowest
          do i=1,nx1
             do j=1,nx2
                do k=1,nx3
                   indx=i*nx2*nx3+j*nx3+k
                   x1_arr(indx)=startx1+dx1*i
                   x2_arr(indx)=startx2+dx2*j
                   x3_arr(indx)=startx3+dx3*k
                enddo
             enddo
          enddo
          r_arr=exp(x1_arr)
! the 2*pi factor is from my conversion choice from HDF5
          ph_arr=x3_arr*2.*pi
          if(IS_MMKS.eq.1) then
             th_arr=calcthmmks(x2_arr,x1_arr)
          else
             th_arr=real(calcthmksh(dble(x2_arr)))
          endif
          write(6,*) 'iharm grid: ',maxval(x1_arr),maxval(x2_arr),minval(th_arr)
        end subroutine read_iharm_grid_file

        subroutine read_iharm_data_file(data_file,rho,p,u,b,kela,mdot)
          character(len=100), intent(in) :: data_file
!          character(len=100) :: data_file_app
          character :: header_byte
          real(8), dimension(:), allocatable :: header
          integer :: nhead,nheader_bytes
          integer :: rhopos,ppos,vpos,bpos,gdetpos,kelapos
!          logical, intent(in), optional :: gridonly
          integer, intent(in), optional :: mdot
!          logical :: gridread
!          real(8), intent(out) :: tcur
          real(8), dimension(:), allocatable, intent(out) :: p,rho,kela
          real(8), dimension(:), allocatable :: udotu,bdotu,bdotb,pmin
          real(8), dimension(:), allocatable :: alpha,gamma,zero
          type (four_vector), dimension(:), allocatable, intent(out) :: u,b
          type (four_vector), dimension(:), allocatable :: uks,bks,beta
          real(8), dimension(:,:), allocatable :: grid
          real(4), dimension(:,:), allocatable :: data
          integer :: i,j, nelem, nelem2
          rhopos=4; ppos=5; vpos=6; bpos=10
          kelapos=14;
          nhead=6
          allocate(header(nhead))
!          end if
! FOR IHARM FORMAT FROM CONVERT_HARM_EBHLIGHT_HDF5_TO_GRTRANS.ipynb
          dlen=13+eHEAT; SDUMP=0
          write(6,*) 'rhopos: ',rhopos,ppos,vpos,bpos
          write(6,*) 'data file: ',data_file
          write(6,*) 'iharm data size: ',dlen,nx1,nx2,nx3
!          write(6,*) 'iharm data file: ',data_file
! reversing this for my iharm format from np.savetxt
          allocate(data(dlen,nx1*nx2*nx3))
          open(unit=8,file=data_file,form='unformatted',access='stream',status='old',action='read')
!          open(unit=8,file=data_file,form='formatted',status='old',action='read')
!          read(8,*) header
!          read(8,*) data
          read(8) data
          close(8)
          write(6,*) 'iharm data: ',data(:,1)
          deallocate(header)
          if(SDUMP.eq.0) then
             allocate(grid(nx1*nx2*nx3,3))
             grid=transpose(data(1:3,:))
             x1_arr=grid(:,1); x2_arr=grid(:,2); x3_arr=grid(:,3)
!             r_arr=grid(:,4); th_arr=grid(:,5); ph_arr=grid(:,6)
             deallocate(grid)
             write(6,*) 'read iharm grid sizes', size(x1_arr), size(x2_arr), &
                  size(x3_arr)
! case where we are using a regular dump file but only to get the grid
! replacing with gdump for metric since sdump files use three-vector v^i, B^i
!             if(gridread) return
          endif
          r_arr=exp(x1_arr)
          ph_arr=x3_arr
          if(IS_MMKS.eq.1) then
             th_arr=calcthmmks(x2_arr,x1_arr)
          else
             th_arr=real(calcthmksh(dble(x2_arr)))
          endif
          allocate(p(nx1*nx2*nx3))
          allocate(rho(nx1*nx2*nx3)); allocate(pmin(nx1*nx2*nx3))
          allocate(u(nx1*nx2*nx3)); allocate(b(nx1*nx2*nx3))
!          allocate(gdet(nx1*nx2*nx3))
          allocate(uks(nx1*nx2*nx3)); allocate(bks(nx1*nx2*nx3))
          rho=data(rhopos,:)
          p=data(ppos,:)
          if(eHEAT.eq.1) then
             allocate(kela(nx1*nx2*nx3))
             kela=data(kelapos,:)
          endif
! for small dumps these are 3-vector velocities that need to be converted using metric from gdump
!          allocate(alpha(nx1*nx2*nx3))
!          allocate(beta(nx1*nx2*nx3))
!          allocate(gamma(nx1*nx2*nx3))
!          allocate(zero(nx1*nx2*nx3))
!          zero(:)=0d0
! UPDATED ASSUMING WE ARE SAVING MKS FOUR VECTORS
          b%data(1)=data(bpos,:)
          b%data(2)=data(bpos+1,:)
          b%data(3)=data(bpos+2,:)
          b%data(4)=data(bpos+3,:)
          u%data(1)=data(vpos,:)
          u%data(2)=data(vpos+1,:)
          u%data(3)=data(vpos+2,:)
          u%data(4)=data(vpos+3,:)
          deallocate(data)
!          write(6,*) 'iharm data done assigning rho,p,b,u'
! should check on this code for getting b^t, u^t
! minimally would need to calculate the metric
!          alpha = 1d0/sqrt(-metric_con(1,:))
!          beta%data(1)=1d0
!          beta%data(2)=metric_con(2,:)*alpha*alpha
!          beta%data(3)=metric_con(3,:)*alpha*alpha
!          beta%data(4)=metric_con(4,:)*alpha*alpha
!          call assign_metric(u,metric_cov)
!          call assign_metric(b,metric_cov)
!          gamma=u*u
!          gamma=sqrt(1d0+merge(gamma,zero,gamma.gt.0d0))
!          u%data(1) = gamma/alpha
!          u%data(2) = u%data(2)-gamma*beta%data(2)/alpha
!          u%data(3) = u%data(3)-gamma*beta%data(3)/alpha
!          u%data(4) = u%data(4)-gamma*beta%data(4)/alpha
!          deallocate(alpha); deallocate(beta); deallocate(gamma)
!          deallocate(zero)
!          b%data(1) = u*b
!          b%data(2) = (b%data(2)+b%data(1)*u%data(2))/u%data(1)
!          b%data(3) = (b%data(3)+b%data(1)*u%data(3))/u%data(1)
!          b%data(4) = (b%data(4)+b%data(1)*u%data(4))/u%data(1)
! test code for these transformations
!          allocate(bdotu(n)); allocate(udotu(n)); allocate(bdotb(n))
!          udotu=abs(u*u+1.)
!          write(6,*) 'udotu: ',maxval(udotu),udotu(30*nx1*nx2+50)
!          bdotb = b*b; bdotu = b*u
!          write(6,*) 'bdotu: ',maxval(abs(bdotu)),bdotu(30*nx1*nx2+50), &
!               bdotb(30*nx1*nx2+50)
!          deallocate(bdotu); deallocate(udotu); deallocate(bdotb)
!          p=merge(p,pmin,p.gt.pmin)
!          write(6,*) "min vals rho", minval(rho)
!          write(6,*) "min vals p", minval(p)
!          write(6,*) "min vals u", minval(u%data(1)), minval(u%data(2)), minval(u%data(3)), minval(u%data(4))
!          write(6,*) "max vals u", maxval(u%data(1)), maxval(u%data(2)), maxval(u%data(3)), maxval(u%data(4))
!          write(6,*) 'read harm transform coords r ', minval(r_arr), maxval(r_arr), asim
          write(6,*) 'read harm transform coords th ',minval(x2_arr), maxval(x2_arr),minval(th_arr),maxval(th_arr)
          write(6,*) 'read harm transform coords phi',minval(x3_arr), maxval(x3_arr),minval(ph_arr),maxval(ph_arr)
          write(6,*) 'read harm transform coords u ',minval(u%data(1)),minval(u%data(2))
          write(6,*) 'read harm transform coords u size ',size(u),hslope
!          if(BL.eq.3) then
          if(IS_MMKS.eq.1) then
             uks = ummks2uks(u,dble(x1_arr),dble(x2_arr))
             bks = ummks2uks(b,dble(x1_arr),dble(x2_arr))
          else
             uks = umksh2uks(u,r_arr,x2_arr)
             bks = umksh2uks(b,r_arr,x2_arr)
          end if
          write(6,*) 'after uks ',minval(uks%data(1))
          u = uks2ubl(uks,dble(r_arr),dble(asim))
          write(6,*) 'read harm transform coords b', minval(u%data(3)),maxval(u%data(2)),&
               minval(b%data(1)),maxval(b%data(4))
          write(6,*) "min vals u", minval(u%data(1)), minval(u%data(2)), minval(u%data(3)), minval(u%data(4))
          write(6,*) "max vals u", maxval(u%data(1)), maxval(u%data(2)), maxval(u%data(3)), maxval(u%data(4))
          b   = uks2ubl(bks,dble(r_arr),dble(asim))
! test code...
          allocate(bdotu(n)); allocate(udotu(n)); allocate(bdotb(n))
          call assign_metric(u,transpose(kerr_metric(r_arr,th_arr,real(asim))))
          call assign_metric(b,transpose(kerr_metric(r_arr,th_arr,real(asim))))
          where(r_arr.gt.(1.+sqrt(1.-asim**2.)))
             udotu=abs(u*u+1.)
             bdotu=abs(u*b)
             bdotb=abs(b*b)
          elsewhere
             udotu=0.
             bdotu=0.
             bdotb=0.
          endwhere
          write(6,*) 'after transform', rho(1:4), p(1:4)
          write(6,*) 'after transform',p(5*nx1+23), rho(7*nx1+44), &
               maxval(abs(udotu+1.)), maxval(abs(bdotu)),minval(bdotb)
          write(6,*) 'udotu: ',udotu(30*nx1*nx2+50:30*nx1*nx2+60)
          write(6,*) 'bdotb: ',bdotb(30*nx1*nx2+50:30*nx1*nx2+60)
          write(6,*) 'bdotu: ',bdotu(30*nx1*nx2+50:30*nx1*nx2+60)
          write(6,*) 'x1: ',x1_arr(30*nx1*nx2+50),x2_arr(30*nx1*nx2+60)
!          write(6,*) 'bmag: ',bmag(123),bmag(244)
          deallocate(udotu); deallocate(bdotu); deallocate(bdotb)
          deallocate(uks); deallocate(bks); deallocate(pmin)
        end subroutine read_iharm_data_file

        subroutine initialize_iharm_model(a,ifile,df,hf,gf,ntt,indft)
        real(kind=8), intent(in) :: a
        integer :: nx, status, nhead
        character(len=20), intent(in), optional :: ifile
        character(len=20) :: default_ifile='harm.in'
        character(len=100), intent(in), optional :: df,hf,gf
        integer, intent(in), optional :: ntt,indft
!        if(SDUMP) then
!           nhead=46
!        else
!           nhead=70
!        endif
! if all inputs are given, assign variables rather than reading from file
        if (present(df)) then
           dfile = df
           hfile = hf
           gfile = gf
           nt = ntt
           indf = indft
        else
           if (present(ifile)) then
              call read_iharm_inputs(ifile)
           else
              call read_iharm_inputs(default_ifile)
           endif
        endif
        write(6,*) 'iharm fluid model: ',hfile,dfile,gfile,nt,indf
        call read_iharm_data_header()
        if (abs(asim-a).gt.1e-4) then 
           write(6,*) 'ERROR -- Different simulation and grtrans spin values!', asim, a
           return
        endif
        write(6,*) 'read header', SDUMP, BL, eCOND, asim
        n=nx1*nx2*nx3
        call init_iharm_data(n,n*nt)
        write(6,*) 'init data', nt, indf, hfile, dfile
        call load_iharm_data(nt)
        end subroutine initialize_iharm_model

        subroutine load_iharm_data(nt)
        real(8), dimension(:), allocatable :: rho,p,kela
        real, dimension(:), allocatable :: vrl, vtl, vpl
        type (four_vector), dimension(:), allocatable :: u, b
        integer, intent(in) :: nt
        character(len=20) :: append,ending
        character(len=100) :: data_file
        integer :: k
        real(8) :: tstep_test
        allocate(rho(n)); allocate(p(n))
        allocate(vrl(n)); allocate(vtl(n)); allocate(vpl(n))
        allocate(u(n)); allocate(b(n))
        if(eHEAT.eq.1) then
           allocate(kela(n))
        endif
! if we're using small dumps then get grid and metric from a regular dump file
        write(6,*) 'eHEAT: ',eHEAT
! my current scheme saves grid and data together. probably should 
! just switch to an sdump scheme right away though with library, ebhlight
! in mind with many files
        if(SDUMP.eq.1) then
           call init_iharm_grid_data(nx1*nx2*nx3)
           call read_iharm_grid_file(gfile)
        endif
        write(6,*) 'grid: ',SDUMP,nt
! now loop over and load data files
!        ending='.out'
        ending='.bin'
        do k=1,nt
           write(append, fmt='(I5.3)') indf-(k-1)
           data_file = trim(dfile) // trim(adjustl(append)) // trim(ending)
           write(6,*) 'data_file: ',indf-(k-1),append,data_file
           call read_iharm_data_file(data_file,rho,p,u,b,kela)
           t(k)=tcur
           write(6,*) 'after iharm data: ',tcur,IS_MMKS
           call lnrf_frame(real(u%data(2)/u%data(1)),real(u%data(3)/u%data(1)), & 
                real(u%data(4)/u%data(1)),r_arr,asim,th_arr,vrl,vtl,vpl)
!           write(6,*) 'lnrf transform', size(b0_arr), size(b%data(1)), size(b0_arr((k-1)*n+1:k*n))
!           write(6,*) 'lnrf transform', size(vrl), size(vtl), size(vpl), size(u0_arr)
           ! now assign to data arrays
           rho_arr((k-1)*n+1:k*n)=rho
! conversion between internal energy and pressure
! REMOVING THIS TO USE INTERNAL ENERGY IN CALCULATING T_E FOR IHARM!
! f%p will be UU!
!           p_arr((k-1)*n+1:k*n)=p*(real(gam)-1.)
           p_arr((k-1)*n+1:k*n)=p
           b0_arr((k-1)*n+1:k*n)=b%data(1)
           br_arr((k-1)*n+1:k*n)=b%data(2)
           bth_arr((k-1)*n+1:k*n)=b%data(3)
           bph_arr((k-1)*n+1:k*n)=b%data(4)
!           write(6,*) 'v assign'
           u0_arr((k-1)*n+1:k*n)=u%data(1)
!           write(6,*) 'vrl assign'
           vrl_arr((k-1)*n+1:k*n)=real(vrl)
!           write(6,*) 'after assign', size(vtl_arr((k-1)*n+1:k*n)), size(vtl)
           vpl_arr((k-1)*n+1:k*n)=real(vpl)
!           write(6,*) 'after vpl', vtl(1), vtl(n), vtl
           vtl_arr((k-1)*n+1:k*n)=real(vtl)
           if(eHEAT.eq.1) then
              kela_arr((k-1)*n+1:k*n)=kela
           endif
        end do
!        tstep=10.
! Need to be able to do light curves even when nload (nt) = 1... Maybe read headers for first two files.
        if (nt.gt.1) then
! use default sim time step unless told otherwise
           tstep=t(1)-t(2)
           tstep_test=t(nt-1)-t(nt)
           if (abs((tstep-tstep_test)/tstep).gt.0.1) then 
              write(6,*) 'WARNING -- Time step changes over dumps!'
           endif
        endif
!        write(6,*) 'after loop', maxval(abs(u*u+1.))
        if(SDUMP.eq.1) call del_iharm_grid_data()
        deallocate(rho); deallocate(p)
        deallocate(vrl); deallocate(vtl); deallocate(vpl)
        deallocate(u); deallocate(b)
        if(eHEAT.eq.1) then
           deallocate(kela)
        endif
!        write(6,*) 'read data file', a
        write(6,*) maxval(vrl_arr**2.+vtl_arr**2.+vpl_arr**2.)
!        write(6,*) 'maxmin temp: ',minval(p_arr/rho_arr*1.67e-24*9e20/1.38e-16), &
!             maxval(p_arr/rho_arr*1.67e-24*9e20/1.38e-16)
        end subroutine  load_iharm_data

        subroutine advance_iharm_timestep(dtobs)
        real(8), intent(in) :: dtobs
        integer :: nupdate
        nupdate = floor(dtobs / tstep)
        toffset = toffset+dtobs-tstep*nupdate
! after a couple stpes might be accumulating offset
        nupdate=nupdate+floor(toffset/tstep)
        toffset=toffset-floor(toffset/tstep)*tstep
        indf=indf+nupdate
        write(6,*) 'advance harm timestep: ',dtobs,tstep,toffset,indf,nupdate
        call update_iharm_data(nupdate)
        end subroutine advance_iharm_timestep

        subroutine update_iharm_data(nupdate)
        integer, intent(in) :: nupdate
        integer :: nshift
        if(nupdate.gt.nt) then
! this is case where we just load all new data
           call load_iharm_data(nt)
        else
! this is case where we shift existing data back and then load more
           nshift=n*nupdate+1
!           rho_arr=rho_arr(nshift:nshift-1:1)
!           allocate(temp_arr(n))
           rho_arr(nshift:n*nt)=rho_arr(1:n*nt-nshift)
           p_arr(nshift:n*nt)=p_arr(1:n*nt-nshift)
           br_arr(nshift:n*nt)=br_arr(1:n*nt-nshift)
           bth_arr(nshift:n*nt)=bth_arr(1:n*nt-nshift)
           bph_arr(nshift:n*nt)=bph_arr(1:n*nt-nshift)
           b0_arr(nshift:n*nt)=b0_arr(1:n*nt-nshift)
           vrl_arr(nshift:n*nt)=vrl_arr(1:n*nt-nshift)
           vtl_arr(nshift:n*nt)=vtl_arr(1:n*nt-nshift)
           vpl_arr(nshift:n*nt)=vpl_arr(1:n*nt-nshift)
           if(eHEAT.eq.1) then
              kela_arr(nshift:n*nt)=kela_arr(1:n*nt-nshift)
           endif
           call load_iharm_data(nupdate)
        endif
        end subroutine update_iharm_data

        subroutine init_iharm_data(nx,n)
        integer, intent(in) :: n,nx
        allocate(rho_arr(n)); allocate(p_arr(n)); allocate(u0_arr(n))
        allocate(vrl_arr(n)); allocate(vtl_arr(n)); allocate(vpl_arr(n))
        allocate(b0_arr(n)); allocate(br_arr(n)); allocate(bth_arr(n))
        allocate(bph_arr(n))
        allocate(x1_arr(nx)); allocate(x2_arr(nx)); allocate(r_arr(nx))
        allocate(th_arr(nx)); allocate(t(nt)); allocate(x3_arr(nx))
        allocate(ph_arr(nx))
        if(eHEAT.eq.1) allocate(kela_arr(n))
        end subroutine init_iharm_data

        subroutine del_iharm_data()
        deallocate(rho_arr); deallocate(p_arr); deallocate(u0_arr)
        deallocate(vrl_arr); deallocate(vtl_arr); deallocate(vpl_arr)
        deallocate(b0_arr); deallocate(br_arr); deallocate(bth_arr)
        deallocate(bph_arr); deallocate(x1_arr); deallocate(x2_arr)
        deallocate(r_arr); deallocate(th_arr); deallocate(t)
        deallocate(x3_arr)
        deallocate(ph_arr)
        if(eHEAT.eq.1) deallocate(kela_arr)
        end subroutine del_iharm_data

        subroutine init_iharm_grid_data(nx)
          integer, intent(in) :: nx
          allocate(gdet(nx)); allocate(metric_cov(10,n))
          allocate(metric_con(10,n)); allocate(drdx(16,n))
        end subroutine init_iharm_grid_data
        
        subroutine del_iharm_grid_data()
          deallocate(metric_cov)
          deallocate(metric_con)
          deallocate(gdet)
          deallocate(drdx)
        end subroutine del_iharm_grid_data

      end module fluid_model_iharm
