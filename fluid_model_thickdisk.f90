   module fluid_model_thickdisk

      use class_four_vector
      use interpolate, only: interp
      use kerr, only: kerr_metric, lnrf_frame, uks2ubl, bl2ks, surf_integral
      use phys_constants, only: pi
      use math, only: zbrent
      implicit none

      namelist /thickdisk/  dfile, gfile, nt, nfiles, indf, jonfix, offset, sim, dindf, magcrit
      integer :: nx1, nx2, nx3, ndumps, nt, nfiles, indf, defcoord, &
           jonfix, dlen, offset, &
           header_length, field_header_length, dindf, magcrit
      integer, parameter  :: file_len = 300, sim_len=40, hstr_len=500
      character(len=file_len) :: dfile, gfile
      character(len=40) :: sim
      integer :: n, readthickgrid
      integer :: ndim=3, nhead=30, test=1, binary=1, DEBUG=0
      integer, dimension(:), allocatable :: dumps
      real :: tstep, h, dx1, dx2, dx3, gam, startx1, startx2, startx3, dt, &
           r0, rin, rout
      real(kind=8) :: xbr,toffset=0.,asim,tcur
      real, dimension(:), allocatable :: t
      real(kind=8), dimension(:), allocatable :: x1_arr, x2_arr, x3_arr, r_arr, th_arr, ph_arr
      real, dimension(:), allocatable :: rho_arr, p_arr, u0_arr, vrl_arr, &
        vpl_arr, vtl_arr
      real, dimension(:), allocatable :: b0_arr, br_arr, bth_arr, bph_arr

      interface init_thickdisk_data
        module procedure init_thickdisk_data
      end interface

      interface calcthmks
         module procedure calcthmks6
         module procedure calcthmks6_single
         module procedure calcthmks6_8
         module procedure calcthmks6_8single
      end interface
      
      interface del_thickdisk_data
        module procedure del_thickdisk_data
      end interface

      interface initialize_thickdisk_model
        module procedure initialize_thickdisk_model
      end interface

      interface transformbl2mks
         module procedure transformbl2mks
      end interface

      interface umks2uks
         module procedure umks2uks
      end interface

      interface thickdisk_vals
        module procedure thickdisk_vals
      end interface

      interface coord_def
         module procedure coord_def
      end interface

      interface calcrmks
         module procedure calcrmks
         module procedure calcrmks_single
      end interface

!      interface 
!         module procedure get_header_lengths
!      end interface

      contains

        subroutine coord_def
! routine to set up global coordinate system variables if desired
        end subroutine coord_def

        function calcrmks(x1,xbr) result(r)
          ! Compute r given x1 for Jon's simulations                                                                 
          ! JAD 7/24/2011
          real(kind=8), intent(in), dimension(:) :: x1
          real(kind=8), intent(in) :: xbr
          real(kind=8) :: npow2
          real(kind=8), dimension(size(x1)) :: r, xi
          npow2=10d0
          where(x1.gt.xbr)
             xi=x1+(x1-xbr)**npow2
          elsewhere
             xi=x1
          endwhere
          r=exp(xi)
        end function calcrmks

        function calcrmks_single(x1,xbr) result(r)
          ! Compute r given x1 for Jon's simulations                                                                 
          ! JAD 7/24/2011
          real(kind=8), intent(in) :: x1,xbr
          real(kind=8) :: npow2
          real(kind=8) :: r, xi
          npow2=10d0
          if(x1.gt.xbr) then
             xi=x1+(x1-xbr)**npow2
          else
             xi=x1
          endif
          r=exp(xi)
        end function calcrmks_single

        function umks2uks(fum,x1,x2,xbr) result(uks)
          real(kind=8), intent(in), dimension(:) :: x1,x2
          type (four_vector), dimension(:), intent(in) :: fum
          real(kind=8), intent(in) :: xbr
          type (four_vector), dimension(size(fum)) :: uks
          real(kind=8), dimension(size(x1)) :: r,dr,dx2,dx1,drdx1,dthdrnum,dthdx2num
          ! Convert four-vectors from MKS to KS numerically using central differences.
          ! JAD 5/24/2010
          ! These parameters are from comparing numerical and analytical values 
          ! using the jetcoords3 grid.
          r=calcrmks(x1,xbr)
          dr=1d-4*r; dx2=1d-6*x2; dx1=1d-4*x1
          ! r numerically:
          drdx1=(calcrmks(x1+.5d0*dx1,xbr)-calcrmks(x1-.5d0*dx1,xbr))/dx1
          uks%data(2)=drdx1*fum%data(2)
          ! Now do \theta numerically:
          dthdrnum=(calcthmks(x2,r+.5d0*dr)-calcthmks(x2,r-.5d0*dr))/dr
          dthdx2num=(calcthmks(x2+.5d0*dx2,r)-calcthmks(x2-.5d0*dx2,r))/dx2
          uks%data(3)=fum%data(3)*dthdx2num+uks%data(2)*dthdrnum
          ! simple phi
          uks%data(4)=fum%data(4)*2.*pi
          ! time component doesn't change
          uks%data(1)=fum%data(1)
!          write(6,*) 'umks2uks r: ',r(1),x1(1),x2(1),xbr
!          write(6,*) 'umks2uks dr: ',dr(1),dx1(1),dx2(1)
!          write(6,*) 'umks2uks drdx: ',calcrmks(x1(1:2)+0.5d0*dx1(1:2),xbr)
!          write(6,*) 'umks2uks umks: ',fum(1)%data(1:4)
!          write(6,*) 'umks2uks drdx1: ',drdx1(1)
 !         write(6,*) 'umks2uks dthdrnum: ',dthdrnum(1)
 !         write(6,*) 'umks2uks dthdx2num: ',dthdx2num(1)
 !         write(6,*) 'umks2uks uks: ',uks(1)%data
        end function umks2uks

        function calcthmks6(x2,r) result(theta)
          real, intent(in), dimension(:) :: x2, r
          real, dimension(size(x2)) :: theta,myhslope,th2, &
               th0,switch0,switch2,theta1,theta2,arctan2
          real :: r0r,r1jet,njet,r0jet,rsjet,qjet, &
               rs,r0,r0jet3,rsjet3,h0,ntheta,htheta,rsjet2,r0jet2
          r0r=0.
          r1jet=2.8
          njet=0.3
          r0jet=15.
          rsjet=40.
          qjet=1.3
          rs=40.
          r0=20.
          r0jet3=20.
          rsjet3=0.
          h0=0.3
          njet=1.0
          ntheta=5.
          htheta=0.15
          rsjet2=5.0
          r0jet2=2.0
          myhslope=h0+((r-rsjet3)/r0jet3)**njet
          th2=0.5*pi*(1.+atan(myhslope*(x2-0.5))/atan(myhslope*.5))
          myhslope=2.0-qjet*(r/r1jet)**(-njet*(0.5+1.0/pi*atan(r/r0jet-rsjet/r0jet)))
          th0=pi*x2+((1.-myhslope)*0.5)*sin(2.*pi*x2)
          switch0=0.5+1.0/pi*atan((r-rs)/r0)
          switch2=0.5-1.0/pi*atan((r-rs)/r0)
          theta1=th0*switch2+th2*switch0
          theta2=pi*.5*(htheta*(2.*x2-1.)+(1.-htheta)*(2.*x2-1.)**ntheta+1.)
          arctan2=.5+1./pi*(atan((r-rsjet2)/r0jet2))
          theta=theta2+arctan2*(theta1-theta2)
        end function calcthmks6

        function calcthmks6_single(x2,r) result(theta)
          ! Calculates \theta for Jon's defcoord = 1401
          ! JAD 5/14/2010, fortran 12/28/2012
          real, intent(in) :: x2,r
          real :: th,r0r,r1jet,njet,r0jet,rsjet,qjet, &
               rs,r0,r0jet3,rsjet3,h0,ntheta,htheta,rsjet2,r0jet2,myhslope,th2, &
               th0,switch0,switch2,theta1,theta2,arctan2,theta
          r0r=0.
          r1jet=2.8
          njet=0.3
          r0jet=15.
          rsjet=40.
          qjet=1.3
          rs=40.
          r0=20.
          r0jet3=20.
          rsjet3=0.
          h0=0.3
          njet=1.0
          ntheta=5.
          htheta=0.15
          rsjet2=5.0
          r0jet2=2.0
          myhslope=h0+((r-rsjet3)/r0jet3)**njet
          th2=0.5*pi*(1.+atan(myhslope*(x2-0.5))/atan(myhslope*.5))
          myhslope=2.0-qjet*(r/r1jet)**(-njet*(0.5+1.0/pi*atan(r/r0jet-rsjet/r0jet)))
          th0=pi*x2+((1.-myhslope)*0.5)*sin(2.*pi*x2)
          switch0=0.5+1.0/pi*atan((r-rs)/r0)
          switch2=0.5-1.0/pi*atan((r-rs)/r0)
          theta1=th0*switch2+th2*switch0
          theta2=pi*.5*(htheta*(2.*x2-1.)+(1.-htheta)*(2.*x2-1.)**ntheta+1.)
          arctan2=.5+1./pi*(atan((r-rsjet2)/r0jet2))
          theta=theta2+arctan2*(theta1-theta2)
        end function calcthmks6_single

        function calcthmks6_8(x2,r) result(theta)
          real(kind=8), intent(in), dimension(:) :: x2, r
          real(kind=8), dimension(size(x2)) :: theta,myhslope,th2, &
               th0,switch0,switch2,theta1,theta2,arctan2
          real(kind=8) :: r0r,r1jet,njet,r0jet,rsjet,qjet, &
               rs,r0,r0jet3,rsjet3,h0,ntheta,htheta,rsjet2,r0jet2
          r0r=0.
          r1jet=2.8
          njet=0.3
          r0jet=15.
          rsjet=40.
          qjet=1.3
          rs=40.
          r0=20.
          r0jet3=20.
          rsjet3=0.
          h0=0.3
          njet=1.0
          ntheta=5.
          htheta=0.15
          rsjet2=5.0
          r0jet2=2.0
          myhslope=h0+((r-rsjet3)/r0jet3)**njet
          th2=0.5*pi*(1.+atan(myhslope*(x2-0.5))/atan(myhslope*.5))
          myhslope=2.0-qjet*(r/r1jet)**(-njet*(0.5+1.0/pi*atan(r/r0jet-rsjet/r0jet)))
          th0=pi*x2+((1.-myhslope)*0.5)*sin(2.*pi*x2)
          switch0=0.5+1.0/pi*atan((r-rs)/r0)
          switch2=0.5-1.0/pi*atan((r-rs)/r0)
          theta1=th0*switch2+th2*switch0
          theta2=pi*.5*(htheta*(2.*x2-1.)+(1.-htheta)*(2.*x2-1.)**ntheta+1.)
          arctan2=.5+1./pi*(atan((r-rsjet2)/r0jet2))
          theta=theta2+arctan2*(theta1-theta2)
        end function calcthmks6_8


        function calcthmks6_8single(x2,r) result(theta)
          ! Calculates \theta for Jon's defcoord = 1401
          ! JAD 5/14/2010, fortran 12/28/2012    
          real(kind=8), intent(in) :: x2,r
          real(kind=8) :: th,r0r,r1jet,njet,r0jet,rsjet,qjet, &
               rs,r0,r0jet3,rsjet3,h0,ntheta,htheta,rsjet2,r0jet2,myhslope,th2, &
               th0,switch0,switch2,theta1,theta2,arctan2,theta
          r0r=0.
          r1jet=2.8
          njet=0.3
          r0jet=15.
          rsjet=40.
          qjet=1.3
          rs=40.
          r0=20.
          r0jet3=20.
          rsjet3=0.
          h0=0.3
          njet=1.0
          ntheta=5.
          htheta=0.15
          rsjet2=5.0
          r0jet2=2.0
          myhslope=h0+((r-rsjet3)/r0jet3)**njet
          th2=0.5*pi*(1.+atan(myhslope*(x2-0.5))/atan(myhslope*.5))
          myhslope=2.0-qjet*(r/r1jet)**(-njet*(0.5+1.0/pi*atan(r/r0jet-rsjet/r0jet)))
          th0=pi*x2+((1.-myhslope)*0.5)*sin(2.*pi*x2)
          switch0=0.5+1.0/pi*atan((r-rs)/r0)
          switch2=0.5-1.0/pi*atan((r-rs)/r0)
          theta1=th0*switch2+th2*switch0
          theta2=pi*.5*(htheta*(2.*x2-1.)+(1.-htheta)*(2.*x2-1.)**ntheta+1.)
          arctan2=.5+1./pi*(atan((r-rsjet2)/r0jet2))
          theta=theta2+arctan2*(theta1-theta2)
        end function calcthmks6_8single


        function findx2mks6(x2,args) result(diff)
          ! Calculates \theta for Jon's defcoord = 1401
          ! JAD 5/14/2010, fortran 12/28/2012    
          real(kind=8), intent(in) :: x2
          real(kind=8), intent(in), dimension(:) :: args
          real(kind=8) :: diff,th,r,r0r,r1jet,njet,r0jet,rsjet,qjet, &
               rs,r0,r0jet3,rsjet3,h0,ntheta,htheta,rsjet2,r0jet2,myhslope,th2, &
               th0,switch0,switch2,theta1,theta2,arctan2,theta
          th=args(2) ; r=args(1)
          theta=calcthmks(real(x2),real(r))
          diff=th-theta
        end function findx2mks6

        function findx2mks(x2,args) result(diff)
          ! Calculate theta from mks coordinate x2 (McKinney & Gammie (2006a)) and bl/ks coordinate r.             
          ! JAD 4/9/2009                                              
          ! These parameter definitions are based on Jon's defcoord=9:
          real(kind=8), intent(in), dimension(:) :: args
          real(kind=8), intent(in) :: x2
          real(kind=8) :: diff,th,r,rj,r0j,nj,rsj,q,g,h
          th=args(2); r=args(1)
          rj=2.8
          nj=0.3
          r0j=20.
          rsj=80.
          q=1.3
          g=-nj*(.5+1./pi*atan((r-rsj)/r0j))
          h=2.-q*(r/rj)**g
          if(x2.lt.0.5) then
             diff=th-((pi*x2+(1.-h)/2.*sin(2.*pi*x2)))
          else
             diff=th-(pi*x2-(1.-h)/2.*sin(2.*pi*(1.-x2)))
          endif
        end function findx2mks
        
        function findx1mks(x1,args) result(diff)
          real(kind=8), intent(in), dimension(:) :: args
          real(kind=8), intent(in) :: x1
          real(kind=8) :: xbr,npow2,r,xi,diff
        ! Compute x1 given r for Jon's simulations 
          xbr=args(2)
          npow2=10d0
          r=args(1)
          diff=r-calcrmks(x1,xbr)
        end function findx1mks

        subroutine transformbl2mks(r,th,ph,x1,x2,x3)
        ! transform Boyer-Lindquist coordinates to modified Kerr-Schild coordinates used by thickdisk
        real, intent(in), dimension(:) :: r,th,ph
        real, intent(out), dimension(size(r)) :: x1,x2,x3
        real, dimension(2) :: args
        integer :: i
        do i=1,size(r)
           args(1)=r(i); args(2)=xbr
           x1(i)=zbrent(findx1mks,0d0,10d0,dble(args),1d-6)
           args(2)=th(i)
           x2(i)=zbrent(findx2mks6,0d0,1d0,dble(args),1d-6)
        enddo
        x3=ph/2./pi
        end subroutine transformbl2mks

        subroutine thickdisk_vals(x0,a,rho,p,b,u,bmag)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real(kind=8), dimension(size(x0)) :: done
        real, dimension(size(x0)) :: x2,x1,x3,zm,theta,fone, &
         vpl0,vrl0,vtl0,rd,td,pd,rttd,zr,dzero, &
         vr0,vth0,vph0,bth,dummy,tt,ttd,zt,zphi,zpp
        real(kind=8), dimension(nx1) :: uniqx1,uniqr
        real(kind=8), dimension(nx2) :: uniqx2
        real(kind=8), dimension(nx3) :: uniqx3,uniqp
        real, dimension(size(x0),2**(ndim+1)) :: ppi,rhoi,vrli,vtli, &
         vpli,bri,bthi,bphi,b0i,u0i
        real, dimension(size(x0)) :: dth,minp
        real :: dph
        integer, dimension(size(x0)*2*(2**ndim)) :: indx
        integer, dimension(size(x0)) :: x1l,x1u,x2l,x2u,x3l,x3u,tindx
        integer, dimension(size(x0)) :: lx1,lx2,ux1,ux2, &
         umax,one,lx3,ux3
        integer :: npts,i,maxtindx
        real, dimension(size(x0)), intent(out) :: rho,p,bmag
        type (four_Vector), intent(out), dimension(size(x0)) :: u,b
        ! Interpolates THICKDISK data to input coordinates
        ! JAD 3/20/2009, fortran 12/28/2012
!        write(6,*) 'thickdisk_vals'
        one=1; dzero=0d0; done=1d0; fone=1.
!        write(6,*) 'before uniq: ',nx1,nx2,nx3
        uniqx1=x1_arr(1:nx1)
!        write(6,*) 'between', minval(x0%data(4)), maxval(x0%data(4))
        uniqx2=x2_arr(1:nx1*(nx2-1)+1:nx1)
!        write(6,*) 'uniqx3: ',size(x3_arr(1:nx1*nx2*(nx3-1)+1:nx1*nx2)),nx3
!        write(6,*) 'uniqx3: ',size(x3_arr),nx1*nx2*nx3,nx1*nx2*(nx3-1)+1
        uniqx3=x3_arr(1:nx1*nx2*(nx3-1)+1:nx1*nx2)
! put case statement for coordinate systems here? Then do all transformations?
!        write(6,*) 'after uniqx3'
        uniqr=calcrmks(uniqx1,xbr)
        uniqp=2.*pi*uniqx3
!        write(6,*) 'uniqx1: ',maxval(uniqx1)
!        write(6,*) 'uniqx2: ',uniqx2
!        write(6,*) 'uniqx3: ',uniqp
        npts=size(x0)
        theta=x0%data(3)
        zr=x0%data(2)
        zt=x0%data(1)
!        write(6,*) 'zt: ',zt(1),bl2ks(dble(zr(1)),0d0,dble(a),1)
!        write(6,*) 'zt: ',zt
        tt=bl2ks(dble(zr),dble(zt),dble(a),1)-bl2ks(dble(zr(1)),0d0,dble(a),1)
!        write(6,*) 'tt: ',tt(1:5)
! after converting to KS switch back to Kerr-Schild
        tt=-tt
!        write(6,*) 'tt: ',minval(x0%data(1)),maxval(x0%data(1)),minval(tt),maxval(tt)
        zpp=x0%data(4)
        zphi=bl2ks(dble(zr),dble(zpp),dble(a))
! put thickdisk on [0,2pi]
        zphi=mod(zphi,(2.*pi))
        where(zphi.lt.0.)
           zphi=zphi+2.*pi
        endwhere
        call transformbl2mks(zr,theta,zphi,x1,x2,x3)
        lx1=floor((x1-uniqx1(1))/(uniqx1(nx1)-uniqx1(1))*(nx1-1))+1
        ux1=lx1+1
        lx2=floor((x2-uniqx2(1))/(uniqx2(nx2)-uniqx2(1))*(nx2-1))+1
        ux2=lx2+1
        lx3=floor((x3-uniqx3(1))/(uniqx3(nx3)-uniqx3(1))*(nx3-1))+1
        ux3=lx3+1
        lx2=merge(lx2,one,lx2.ge.1)
        ux2=merge(ux2,nx2,ux2.le.nx2)
        lx3=merge(lx3,one,lx3.ge.1)
        ux3=merge(ux3,nx3,ux3.le.nx3)
! Deal with poles
!        write(6,*) 'poles',size(zr),size(uniqx2(lx2)),size(uniqx2(ux2))
        where(ux2.ne.lx2)
           dth=calcthmks(uniqx2(ux2),dble(zr))-calcthmks(uniqx2(lx2),dble(zr))
        elsewhere
           dth=calcthmks(uniqx2(1)+0d0*dble(zr),dble(zr))
        endwhere
! periodic in phi
        minp=uniqp(lx3)
        where(ux3.gt.nx3)
           ux3=1
        endwhere
        where(lx3.lt.1)
           lx3=nx3
           minp=uniqp(lx3)-2.*pi
        endwhere
! uniform in phi
        pd=(zphi-minp)/(uniqp(2)-uniqp(1))
!        write(6,*) 'pd: ',x3(maxloc(abs(pd))), zphi(maxloc(abs(pd))), uniqp(lx3(maxloc(abs(pd)))), &
!             lx3(maxloc(abs(pd))), pd(maxloc(abs(pd)))
        td=abs(theta-calcthmks(uniqx2(lx2),dble(zr)))/dth
!        write(6,*) 'lx2: ',minval(lx2),maxval(lx2),nx2
        rd=(zr-uniqr(lx1))/(uniqr(ux1)-uniqr(lx1))
! When the geodesic is inside the innermost zone center located outside the horizon, 
! use nearest neighbor rather than interpolate
!        write(6,*) 'fluid model thickdisk a: ',a,(1.+sqrt(1.-a**2.))
!        write(6,*) 'fluid model thickdisk uniqr: ',uniqr(lx1)
        where(uniqr(lx1).le.(1.+sqrt(1.-a**2.)))
           rd=1.
        endwhere
!        write(6,*) 'rd td: ',minval(rd),maxval(rd),minval(td),maxval(td),minval(pd),maxval(pd)
!        write(6,*) 'ux lx: ',minval(lx1),maxval(ux1),minval(lx2),maxval(ux2)
        ! th is fastest changing index
        x1l=lx1-1; x1u=ux1-1
        x2l=(lx2-1)*nx1 ; x2u=(ux2-1)*nx1
        x3l=(lx3-1)*nx1*nx2; x3u=(ux3-1)*nx1*nx2
! Indices for nearest neighbors in 3D
!        write(6,*) 'sindx', maxval(sindx), sindx(npts+1), x1l(1)+x2u(1)+1, size(sindx)
! Now find timestep indices:
        maxtindx=nt-2
        where(floor(tt/tstep).le.maxtindx)
           tindx=floor(tt/tstep)
           ttd=(tt-tindx*tstep)/tstep
        elsewhere
           tindx=maxtindx+1
           ttd=0.
        endwhere
!        write(6,*) 'tindx: ',tindx,n*tindx
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
!        write(6,*) 'after indx: ', minloc(indx), x1(minloc(indx)), x2(minloc(indx)), x3(minloc(indx))
!        write(6,*) 'after indx: ',lx1(minloc(indx)),lx2(minloc(indx)),lx3(minloc(indx))
!        write(6,*) 'after indx: ',x1u(minloc(indx)),x2u(minloc(indx)),x3u(minloc(indx))
!        write(6,*) 'after indx: ',x1l(minloc(indx)),x2l(minloc(indx)),x3l(minloc(indx))
!        write(6,*) 'after indx', maxval(indx),minval(indx),size(rho_arr)!,size(rho_arr(indx))
!        write(6,*) 'after indx',npts,2**(ndim+1),size(rhoi)
        rhoi=reshape(rho_arr(indx),(/npts,2**(ndim+1)/))
!        write(6,*) 'rhoi',size(rho_arr(indx)),npts,2**(ndim+1),size(rhoi)
        vrli=reshape(vrl_arr(indx),(/npts,2**(ndim+1)/))
        vtli=reshape(vtl_arr(indx),(/npts,2**(ndim+1)/))
        vpli=reshape(vpl_arr(indx),(/npts,2**(ndim+1)/))
        ppi=reshape(p_arr(indx),(/npts,2**(ndim+1)/))
        b0i=reshape(b0_arr(indx),(/npts,2**(ndim+1)/))
        bri=reshape(br_arr(indx),(/npts,2**(ndim+1)/))
        bthi=reshape(bth_arr(indx),(/npts,2**(ndim+1)/))
        bphi=reshape(bph_arr(indx),(/npts,2**(ndim+1)/))
        u0i=reshape(u0_arr(indx),(/npts,2**(ndim+1)/))
!        rttd=ttd
        rttd=0.
        rho=merge(interp(rhoi,rttd,pd,rd,td),dzero,x1.gt.uniqx1(1))
!        write(6,*) 'interp test u0i: ',u0i(7041,:)
!        write(6,*) 'interp test u0i: ',vrli(7041,:)
!        write(6,*) 'interp test u0i: ',vtli(7041,:)
!        write(6,*) 'interp test u0i: ',vpli(7041,:)
!        write(6,*) 'interp test pd: ',pd(7041),rd(7041),td(7041)
!        write(6,*) 'interp test rho: ',u0(7041)
        p=merge(interp(ppi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'rho: ', rho, p
        vrl0=merge(interp(vrli,rttd,pd,rd,td),dzero,x1.gt.uniqx1(1))
        vtl0=merge(interp(vtli,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'v'
        b%data(1)=merge(interp(b0i,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
        b%data(2)=merge(interp(bri,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'b'
        b%data(3)=merge(interp(bthi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
        b%data(4)=merge(interp(bphi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'u'
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
!        write(6,*) 'interp test u: ',vrl0(7041),vtl0(7041),vpl0(7041)
!        write(6,*) 'interp test u: ',u(7041)%data
!        write(6,*) 'interp test u: ',vrl0(7041)**2+vtl0(7041)**2+vpl0(7041)**2
!        write(6,*) 'interp test u: ',x1_arr(indx(7041)),x2_arr(indx(7041)), &
!             x3_arr(indx(7041))
!        write(6,*) 'interp test u: ',x1_arr(indx(7041+size(x0)*4)),x2_arr(indx(7041+size(x0)*4)), &
!             x3_arr(indx(7041+size(x0)*4))
!        write(6,*) 'interp test u: ',x1(7041),x2(7041),x3(7041)
        call assign_metric(u,transpose(kerr_metric(zr,real(x0%data(3)) &
        ,a)))
!        write(6,*) 'fluid model thickdisk udotu: ',u*u
!        write(6,*) '7041: ',u(7041)*u(7041)
        if(DEBUG.eq.1) then
        ! output quantities for debugging
           write(6,*) 'debug!'
           open(unit=8,file='fluid_model_thickdisk_debug.out')
           write(8,*) nx1,nx2,nx3,size(x0%data(1)),xbr
           write(8,*) rho_arr
           write(8,*) x0%data(1), x0%data(2), x0%data(3), x0%data(4)
           write(8,*) tt,zr,theta,zphi
           write(8,*) x1,x2,x3
           write(8,*) uniqx1,uniqx2,uniqx3
           write(8,*) x1_arr, x2_arr, x3_arr
           write(8,*) lx1,lx2,lx3,ux1,ux2,ux3
           write(8,*) rd,td,pd
           write(8,*) x1l,x1u,x2l,x2u,x3l,x3u
           write(8,*) indx
!           write(6,*) 'test u0i vrli: ',maxval(abs(u0i-vrli))
           write(8,*) u0i
           write(8,*) vrli
           write(8,*) vtli
           write(8,*) vpli
           write(8,*) u%data(1), u%data(2), u%data(3), u%data(4)
           close(unit=8)
        endif
!        write(6,*) 'thickdisk vals rho: ',rho
!        write(6,*) 'thickdisk vals minloc: ',minloc(rho),rd(minloc(rho)), &
!             td(minloc(rho)),x0(minloc(rho))%data(3),uniqth(lx2(minloc(rho))), &
!             lx2(minloc(rho))
!        write(6,*) 'thickdisk vals rhoi: ',rhoi(minloc(rho),:)
        end subroutine thickdisk_vals

        subroutine read_thickdisk_inputs(ifile)
        character(len=20), intent(in) :: ifile
        open(unit=8,file=ifile,form='formatted',status='old')
        read(8,nml=thickdisk)
        write(6,*) 'read: ',nt,nfiles,dindf
        close(unit=8)
        end subroutine read_thickdisk_inputs

        subroutine read_thickdisk_data_header(hfile,header_len)
        character(len=file_len), intent(in) :: hfile
        character(len=10) :: dstr
        character(len=hstr_len) :: header_str
        integer, intent(out) :: header_len
!        integer, intent(in) :: nhead
        real(kind=8), dimension(nhead) :: header
   real(kind=8)      :: tcur,mbh,qbh
        ! Read thickdisk header file
        ! JAD 11/24/2012 based on previous IDL codes
        ! header format from thickdisk dump.c
        open(unit=8,file=hfile,form='formatted',status='old')
        write(6,*) 'nhead: ',nhead,size(header)
        read(8,*) header
        nx1=header(2)
        nx2=header(3)
        nx3=header(4)
        startx1=header(5)
        startx2=header(6)
        startx3=header(7)
        dx1=header(8)
        dx2=header(9)
        dx3=header(10)
        asim=header(13)
        gam=header(12)
        r0=header(14)
        rin=header(15)
        rout=header(16)
        h=header(17)
        dt=header(18)
        tcur=header(1)
        mbh=header(20)
        qbh=header(21)
! is this true?
        dlen=header(nhead)
        if ((qbh.ne.0.).or.(mbh.ne.1.)) write(6,*) 'ERROR! Unexpected values in thickdisk header' 
        defcoord=header(19)
        write(6,*) 'header: ',defcoord,asim,tcur,qbh,gam
        close(unit=8)
! now re-open to find header length
        open(unit=8,file=hfile,form='formatted',status='old')
        read(8,'(A)') header_str
        close(unit=8)
        write(dstr,fmt='(I4)') dlen
        write(6,*) 'header_str: ',header_str
        write(6,*) 'header_str: ',trim(adjustl(dstr)),dlen
        header_len=index(header_str,trim(adjustl(dstr)),back=.true.)+2
        write(6,*) 'header len: ',header_len
        end subroutine read_thickdisk_data_header

        subroutine read_thickdisk_grid_file(grid_file,tcur,binary,mdot)
        character(len=file_len), intent(in) :: grid_file
        integer, intent(in), optional :: binary
!        integer :: dlen
        integer :: gdetpos,rhopos,ppos,vpos,bpos,glen
        integer, intent(in), optional :: mdot
        real(kind=8), dimension(:), allocatable :: gdet, header, udotu, bdotu, rho, p
        real(kind=8), dimension(:,:), allocatable :: grid, data, tmetric
        real(kind=8), intent(out) :: tcur
!        real(kind=8), dimension(:), allocatable :: tempval
        real(kind=8) :: tsteppartf,realstartx1, &
             realstartx2,realstartx3,dx1,dx2,dx3,gamval,aval,r0val,rinval,routval, &
             hslope,localdt,mbhval,qbhval
        integer(4) :: is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumns,defcoordval, &
             realtotalsize1,realtotalsize2,realtotalsize3,localrealnstep
        character(len=300) :: empty,nstr
        character(len=2) :: endl
        character(len=header_length) :: dummy
        type (four_vector), dimension(:), allocatable :: u,b
        integer :: i, nelem
          ! Read thickdisk data from a file of given line length, 
          ! positions of desired variables
          ! JAD 11/24/2012 based on previous IDL codes
          ! Set grid parameters so we can generalize later:
          !nhead=21 ; dlen=72
          !dlen=34; nhead=26
!          rhopos=4; ppos=5; vpos=13; bpos=21; gdetpos=33
        rhopos=10; ppos=18; vpos=30; bpos=38; gdetpos=52
        write(6,*) 'nhead: ',nhead
          allocate(header(nhead))
          write(6,*) 'grid file: ',grid_file
!          write(6,*) 'data file: ',grid_file
          allocate(grid(nx1*nx2*nx3,6));
          allocate(p(nx1*nx2*nx3)); allocate(rho(nx1*nx2*nx3))
          allocate(u(nx1*nx2*nx3)); allocate(b(nx1*nx2*nx3)); allocate(gdet(nx1*nx2*nx3))
!          allocate(uks(nx1*nx2*nx3)); allocate(bks(nx1*nx2*nx3))
!          write(6,*) 'read thickdisk sizes', dlen, nx2, size(data)
          open(unit=8,file=grid_file,status='old',action='read',form='formatted')
          read(8,*) header
          write(6,*) 'header: ',header
          glen=header(nhead)
          tcur=header(1)
          close(unit=8)
          if(present(binary)) then
             open(unit=8,file=grid_file,status='old',action='read',form='unformatted',access='stream')
!             read(8) tsteppartf,realtotalsize1,realtotalsize2,realtotalsize3,realstartx1,realstartx2, &
!                  realstartx3, dx1, dx2, dx3, localrealnstep, gam, aval, r0val,rinval,routval,hslope,localdt, &
!                  defcoordval,mbhval,qbhval,is,ie,js,je,ks,ke,whichdump,whichdumpversion,numcolumns,endl
             read(8) dummy
!             write(6,*) 'header vals: ',tsteppartf,realtotalsize1,realtotalsize2,realtotalsize3
!             write(6,*) 'header vals: ',numcolumns
             write(6,*) 'header vals: ',dummy
!             allocate(tempval(73))
!             write(nstr,fmt='(I8)') glen*nx1*nx2*nx3
!             write(6,*) 'nstr: ',nstr, glen*nx1*nx2*nx3, nx1, nx2, nx3
!             do i=1,nx1*nx2*nx3
!                read(8,'(73A8)') tempval
!                data(:,i)=tempval
!             enddo
!             read(8,'('//trim(nstr)//'A)') data
!             read(8,'(A)') tempval
!             read(8,'(A)') tempval
!             read(8,'(A8)') data
             allocate(data(glen,nx1*nx2))
             write(6,*) 'data size: ',glen,nx1*nx2*nx3,size(data)
             nelem=nx1*nx2
             do i=0,nx3-1
!                write(6,*) 'i: ',i,nx1,nx2,nx3
                read(8) data
                p(1+i*nelem:(i+1)*nelem)=data(ppos+1,:)
                rho(1+i*nelem:(i+1)*nelem)=data(rhopos+1,:)
                b(1+i*nelem:(i+1)*nelem)%data(1)=(data(bpos+1,:))
                b(1+i*nelem:(i+1)*nelem)%data(2)=(data(bpos+2,:))
                b(1+i*nelem:(i+1)*nelem)%data(3)=(data(bpos+3,:))
                b(1+i*nelem:(i+1)*nelem)%data(4)=(data(bpos+4,:))
                !             write(6,*) 'sizes: ',size(b(1+i*nelem:(i+1)*nelem)%data(4)),size(data(bpos+4,:))
                u(1+i*nelem:(i+1)*nelem)%data(1)=(data(vpos+1,:))
                u(1+i*nelem:(i+1)*nelem)%data(2)=(data(vpos+2,:))
                u(1+i*nelem:(i+1)*nelem)%data(3)=(data(vpos+3,:))
                u(1+i*nelem:(i+1)*nelem)%data(4)=(data(vpos+4,:))
!             write(6,*) 'sizes: ',size(u(1+i*nelem:(i+1)*nelem)%data(4)),size(data(vpos+4,:))
!             write(6,*) 'sizes: ',size(grid(1+i*nelem:(i+1)*nelem,:)),size(data(1:4,:))
                grid(1+i*nelem:(i+1)*nelem,:)=transpose(data(4:9,:))
!             write(6,*) 'sizes: ',size(gdet(1+i*nelem:(i+1)*nelem)),size(data(gdetpos+1,:))
                gdet(1+i*nelem:(i+1)*nelem)=(data(gdetpos+1,:))
             enddo
!             read(8) data
!             rho=data(rhopos+1,:)
!             p=data(ppos+1,:)
!             u%data(1)=data(vpos+1,:)
!             u%data(2)=data(vpos+2,:)
!             u%data(3)=data(vpos+3,:)
!             u%data(4)=data(vpos+4,:)
!             b%data(1)=data(bpos+1,:)
!             b%data(2)=data(bpos+2,:)
!             b%data(3)=data(bpos+3,:)
!             b%data(4)=data(bpos+4,:)
!             grid=transpose(data(4:9,:))
!             gdet=data(gdetpos+1,:)
             close(unit=8)
           else
             open(unit=8,file=grid_file,status='old',action='read',form='formatted')
             read(8,*) header
!             write(6,*) 'header: ',header
             glen=header(nhead)
             tcur=header(1)
             allocate(data(dlen,nx2)); nelem=nx1
             do i=0,(nx1*nx2*nx3)/nelem-1
                read(8,*) data
!             if (i.eq.0) write(6,*) 'thickdisk data: ', data(:,1)
!             write(6,*) 'data loop: ',i,size(rho),size(p),size(b)
!             write(6,*) 'data vals: ',data(:,1),data(:,gdetpos+1)
                rho(1+i*nelem:(i+1)*nelem)=data(rhopos+1,:)
!             write(6,*) 'sizes: ',size(rho(1+i*nelem:(i+1)*nelem)),size(data(rhopos+1,:))
                p(1+i*nelem:(i+1)*nelem)=data(ppos+1,:)
                b(1+i*nelem:(i+1)*nelem)%data(1)=(data(bpos+1,:))
                b(1+i*nelem:(i+1)*nelem)%data(2)=(data(bpos+2,:))
                b(1+i*nelem:(i+1)*nelem)%data(3)=(data(bpos+3,:))
                b(1+i*nelem:(i+1)*nelem)%data(4)=(data(bpos+4,:))
                !             write(6,*) 'sizes: ',size(b(1+i*nelem:(i+1)*nelem)%data(4)),size(data(bpos+4,:))
                u(1+i*nelem:(i+1)*nelem)%data(1)=(data(vpos+1,:))
                u(1+i*nelem:(i+1)*nelem)%data(2)=(data(vpos+2,:))
                u(1+i*nelem:(i+1)*nelem)%data(3)=(data(vpos+3,:))
                u(1+i*nelem:(i+1)*nelem)%data(4)=(data(vpos+4,:))
!             write(6,*) 'sizes: ',size(u(1+i*nelem:(i+1)*nelem)%data(4)),size(data(vpos+4,:))
!             write(6,*) 'sizes: ',size(grid(1+i*nelem:(i+1)*nelem,:)),size(data(1:4,:))
                grid(1+i*nelem:(i+1)*nelem,:)=transpose(data(4:9,:))
!             write(6,*) 'sizes: ',size(gdet(1+i*nelem:(i+1)*nelem)),size(data(gdetpos+1,:))
                gdet(1+i*nelem:(i+1)*nelem)=(data(gdetpos+1,:))
                !             write(6,*) 'data loop'
             end do
          endif
!          write(6,*) 'after read loop'
          close(unit=8)
          deallocate(data); deallocate(header)
!          write(6,*) 'grid sizes', size(x1_arr), size(x2_arr), size(r_arr), size(th_arr)
          x1_arr=grid(:,1); x2_arr=grid(:,2); x3_arr=grid(:,3)
          r_arr=grid(:,4); th_arr=grid(:,5); ph_arr=grid(:,6)
          write(6,*) 'assign grid', r_arr(1), th_arr(1), ph_arr(1)
          write(6,*) 'assign grid', x1_arr(1), x2_arr(1), x3_arr(1)
!          deallocate(grid)
!          if (present(mdot)) then
             ! Calculate accretion rate in code units:
             !  nx1=n_elements(uniqx1) ; nx2=n_elements(uniqx2) ; nz=n_elements(uniqx3) 
!             dx2=uniqx2(2)-uniqx2(1) ; dx3=uniqx3(2)-uniqx3(1)
!             mdotarr=-1.*sum(sum(reform(gdet*rho*v(:,1),nx1,nx2,nz),3),2)*dx2*dx3
!          mdot=surf_integral(-u%data(2)*rho,r_arr,th_arr,ph_arr,a,gdet)
!          endif
          ! Transform velocities, magnetic fields from MKS to KS and then BL:
 !         write(6,*) 'transform coords'
!          uks = umks2uks(u,r_arr,x2_arr,xbr)
!          u = uks2ubl(uks,dble(r_arr),dble(asim))
!          bks = umks2uks(b,r_arr,x2_arr,xbr)
!          b   = uks2ubl(bks,dble(r_arr),dble(asim))
! test code...
!          allocate(bdotu(n)); allocate(udotu(n))
!          call assign_metric(u,transpose(kerr_metric(r_arr,th_arr,real(asim))))
!          call assign_metric(b,transpose(kerr_metric(r_arr,th_arr,real(asim))))
!          where(r_arr.gt.(1.+sqrt(1.-asim**2.)))
!             udotu=abs(u*u+1.)
!             bdotu=abs(u*b)
!          elsewhere
!             udotu=0.
!             bdotu=0.
!          endwhere
!          write(6,*) 'after transform', rho(1:4), p(1:4)
!          write(6,*) 'after transform',p(5*nx1+23), rho(7*nx1+131), &
!               maxval(udotu), maxval(bdotu)
!          deallocate(udotu); deallocate(bdotu)
        deallocate(grid)
        deallocate(p); deallocate(rho); deallocate(gdet); deallocate(u); deallocate(b)
        end subroutine read_thickdisk_grid_file

        subroutine read_mb09_grid_file(nt)
        integer, intent(in) :: nt
        integer :: n
          open(unit=8,file=gfile,form='unformatted',status='old')
          read(8) nx1,nx2,nx3
          n=nx1*nx2*nx3
          write(6,*) 'read grid: ',size(x1_arr),size(x2_arr),size(x3_arr)
          read(8) x1_arr
          read(8) x2_arr
          read(8) x3_arr
          close(unit=8)
        end subroutine read_mb09_grid_file

        subroutine read_thickdisk_fieldline_file(data_file,tcur,rho,p,u,b,jonfix,binary,test)
        character(len=file_len), intent(in) :: data_file
        type (four_vector), intent(inout), dimension(:) :: u,b
        real(kind=8), dimension(:), allocatable, intent(inout) :: rho,p
        integer, intent(in), optional :: jonfix,test,binary
        real(kind=8), intent(out) :: tcur
        real(kind=8), dimension(:), allocatable :: header,b0,bsqorho,condmaxbsqorhorhs, &
             one,zero,rinterp
        real(4), dimension(:,:), allocatable :: data
        type (four_vector), dimension(:), allocatable :: uks,bks
        real(kind=8) :: maxbsqorhofar,maxbsqorhohigh,maxbsqorhonear,maxuu0high
        integer :: n,rhopos,ppos,vpos,bpos,i,nelem
        character(len=field_header_length) :: dummy
!           dlen=11
           rhopos=1 ; ppos=2 ; vpos=5 ; bpos=9
!           allocate(header(nhead))
!           t=header(1)
!           a=header(12)
!           nx=header(1) ; ny=header(2) ; nz=header(3)
!           deallocate(header)
           n=nx1*nx2*nx3; nelem=nx1
           ! Create desired data arrays
!           allocate(b(n)); allocate(u(n)); allocate(rho(n))
!           allocate(p(n)); allocate(bks(n)); allocate(uks(n)); allocate(b0(n))
           allocate(uks(n)); allocate(b0(n)); allocate(bks(n))
!           defcoord=header(18); rout=header(15)
           if(rout.gt.1.e3) then
              xbr=log(500.)
           else 
              xbr=log(1.e5)
           endif
!           select case (defcoord)
!              case (9)
!                 thfunc='calcthmks'
!              case(1401)
!                 thfunc='calcthmks6'
!              case default
!                 write(6,*) 'WARNING -- Coordinate system not recognized!'
!            end select
!           if n_elements(cfile) ne 0 then
!           ! Read in coordinate parameter file:
!           openr, lunc, cfile, /get_lun
!           coords=''
!           readf,lunc,coords
!           coords=strsplit(coords,' ',/extract) ; ncoords=n_elements(coords)
!           free_lun, lunc
!        endif
           th_arr=calcthmks(x2_arr,calcrmks(x1_arr,xbr))
           b%data(1)=0.
           if(present(binary)) then
              open(unit=8,file=data_file,form='unformatted',access='stream')
              read(8) dummy
              write(6,*) 'dummy: ',dummy
              write(6,*) 'dlen: ',dlen,n
              allocate(data(dlen,n))
              read(8) data
              write(6,*) 'field data: ',data(:,1)
              rho=data(rhopos,:)
              p=data(ppos,:)
              b%data(2)=data(bpos,:)
              b%data(3)=data(bpos+1,:)
              b%data(4)=data(bpos+2,:)
              u%data(1)=data(vpos,:)
              u%data(2)=data(vpos+1,:)
              u%data(3)=data(vpos+2,:)
              u%data(4)=data(vpos+3,:)
           else
              allocate(header(nhead))
              open(unit=8,file=data_file,form='formatted')
              read(8) header
              deallocate(header)
              allocate(data(dlen,nelem))
              do i=0,n/nelem-1
                 read(8) data
                 rho(i*nelem+1:(i+1)*nelem)=data(rhopos,:)
                 p(i*nelem+1:(i+1)*nelem)=data(ppos,:)
                 b(i*nelem+1:(i+1)*nelem)%data(2)=data(bpos,:)
                 b(i*nelem+1:(i+1)*nelem)%data(3)=data(bpos+1,:)
                 b(i*nelem+1:(i+1)*nelem)%data(4)=data(bpos+2,:)
!                 b(i*nelem+1:i*nelem+nelem,:)=transpose(data(bpos:bpos+2,:))
                 u(i*nelem+1:(i+1)*nelem)%data(1)=data(vpos,:)
                 u(i*nelem+1:(i+1)*nelem)%data(2)=data(vpos+1,:)
                 u(i*nelem+1:(i+1)*nelem)%data(3)=data(vpos+2,:)
                 u(i*nelem+1:(i+1)*nelem)%data(4)=data(vpos+3,:)
!                 v(i*nelem+1:i*nelem+nelem,:)=transpose(data(vpos:vpos+3,:))
              enddo
            endif
            close(unit=8)
            write(6,*) 'th: ',th_arr(1),calcthmks(x2_arr(1:2),calcrmks(x1_arr(1:2),xbr))
            write(6,*) 'data: ',data(vpos,2223), u(2223)%data(1)
            deallocate(data)
            ! Convert energy density to pressure:
            !gamma=header(11)
            p=(gam-1.)*p
            ! Get four velocity in MKS coords from transport velocity:
            !v(*,1)=v(*,1)*v(*,0) ; v(*,2)=v(*,2)*v(*,0) ; v(*,3)=v(*,3)*v(*,0)
            u%data(2)=u%data(2)*u%data(1)
            u%data(3)=u%data(3)*u%data(1)
            u%data(4)=u%data(4)*u%data(1)
            write(6,*) 'umks: ',u(1)%data
            ! Transform from MKS to BL via KS:
            uks=umks2uks(u,x1_arr,x2_arr,xbr)
            u=uks2ubl(uks,dble(calcrmks(x1_arr,xbr)),dble(asim))
            ! First convert what we've got to KS:
            write(6,*) 'bmks: ',b(1)%data, maxval(abs(r_arr-calcrmks(x1_arr,xbr))/r_arr)
            bks=umks2uks(b,x1_arr,x2_arr,xbr)
!            write(6,*) 'bks: ',bks(1)%data,uks(1)%data,r_arr(1),x2_arr(1)
            ! Recover time-component of four-vector magnetic field:
            call assign_metric(uks,transpose(kerr_metric(r_arr,th_arr,ph_arr,asim)))
            write(6,*) 'uksdotuks: ',maxval(abs(uks*uks+1.))
            write(6,*) 'uks: ',uks(1)%data
            write(6,*) 'uks: ',uks(1)%metric
            write(6,*) 'ks metric: ',r_arr(1:2),th_arr(1:2),ph_arr(1:2) &
                 ,kerr_metric(r_arr(1:2),th_arr(1:2),ph_arr(1:2),asim)
            call assign_metric(bks,transpose(kerr_metric(r_arr,th_arr,ph_arr,asim)))
            b0=bks*uks
            write(6,*) 'bks orig: ',bks(449827)%data
            bks%data(1)=b0; bks%data(2)=(bks%data(2)+b0*uks%data(2))/uks%data(1)
            bks%data(3)=(bks%data(3)+b0*uks%data(3))/uks%data(1)
            bks%data(4)=(bks%data(4)+b0*uks%data(4))/uks%data(1)
            write(6,*) 'b0: ',b0(449827)
            write(6,*) 'b0 uks: ',uks(449827)%data
            write(6,*) 'bks: ',bks(449827)%data
            deallocate(uks)
            b=uks2ubl(bks,dble(calcrmks(x1_arr,xbr)),dble(asim)); deallocate(b0); deallocate(bks)
            call assign_metric(b,transpose(kerr_metric(r_arr,th_arr,asim)))
            call assign_metric(u,transpose(kerr_metric(r_arr,th_arr,asim)))
            ! Apply Jon's approximate fix for the density/energy floors if desired:
            write(6,*) 'fluid model thickdisk jonfix: ',jonfix
            if(jonfix.eq.1) then
               write(6,*) 'fluid model thickdisk performing jonfix'
               maxbsqorhohigh=45.
               maxbsqorhonear=30.
               maxbsqorhofar=10.
               maxuu0high=50.
               allocate(rinterp(n)); allocate(one(n)); allocate(zero(n))
               one=1d0; zero=0d0
               allocate(condmaxbsqorhorhs(n)); allocate(bsqorho(n))
               rinterp=(calcrmks(x1_arr,xbr)-9d0)*(1d0-0d0)/(0d0-9d0) ! gives 0 for use near 9   gives 1 for use near 0
!               write(6,*) 'fluid model thickdisk rmks: ', calcrmks(x1_arr,xbr), xbr
!               write(6,*) 'fluid model thickdisk rinterp: ',rinterp
               rinterp = merge(merge(rinterp,one,rinterp.le.one),zero,rinterp.ge.zero)
               write(6,*) 'rinterp: ',minval(rinterp),maxval(rinterp)
               !    rinterp(rinterp>1.0)=1.0
               !    rinterp(rinterp<0.0)=0.0
               condmaxbsqorhorhs=rinterp*maxbsqorhonear + (1d0-rinterp)*maxbsqorhofar
               bsqorho=(b*b)/rho
               write(6,*) 'bsqrho: ',maxval(bsqorho),minval(bsqorho)
               where((bsqorho.gt.maxbsqorhonear).or.(bsqorho.ge.condmaxbsqorhorhs))
                  rho=1e-18
                  p=1e-18
               endwhere
               deallocate(rinterp); deallocate(one); deallocate(zero); deallocate(condmaxbsqorhorhs)
               deallocate(bsqorho)
            endif
            if(present(test)) then
               write(6,*) 'vals: ',xbr,rout,b(1)%data,u(1)%data,x1_arr(1),x2_arr(1),x3_arr(1)
               write(6,*) 'bmin: ',minval(b*b),maxval(b*b), b(1)*b(1)
               write(6,*) 'fluid model thickdisk udotu: ',maxval(abs(u*u+1.)), u(1)*u(1)+1.
            endif
        end subroutine read_thickdisk_fieldline_file

        subroutine calc_thickdisk_mdot(mdot)
        real(kind=8), dimension(nx1), intent(out) :: mdot
        real(kind=8) :: dph
        real(kind=8), dimension(nx1,nx2) :: dth
        real(kind=8), dimension(nx1,nx2,nx3) :: gdet
        mdot=surf_integral(reshape(dble(-rho_arr*vrl_arr),(/nx1,nx2,nx3/)),&
             reshape(r_arr,(/nx1,nx2,nx3/)),&
             reshape(th_arr,(/nx1,nx2,nx3/)),reshape(ph_arr,(/nx1,nx2,nx3/)),asim)
        end subroutine

        subroutine calc_thickdisk_jet_power(jet_power,wind_power,therm,kinetic,jet_magenergy, &
             jet_vel,jetmom,jet_totpower,gaminf_az,wind_therm,wind_kinetic,jetpower2d,therm2d, &
             thermout2d, jetvel2d, jetmom2d, totpower2d,jetpower2dth,therm2dth,thermout2dth, &
             totpower2dth)
        ! calculate jet and wind power using definitions from mckinney et al. 2012
        ! this assumes data loaded is four-velocity rather than lnrf vels
        real(kind=8), dimension(nx1), intent(out) :: jet_power, wind_power, therm, kinetic, &
             jet_magenergy,jet_vel,jetmom,jet_totpower,wind_therm,wind_kinetic
        real(kind=8), dimension(nx1,nx2), intent(out) :: gaminf_az,jetpower2d,therm2d,thermout2d, &
             jetvel2d, jetmom2d, totpower2d,therm2dth,thermout2dth,totpower2dth,jetpower2dth
        integer, dimension(nx1,nx2,nx3) :: wjet, wwind
        type (four_vector), dimension(nx1*nx2*nx3) :: ucov,bcov,u,b
        real(kind=8), dimension(nx1,nx2,nx3) :: trtem,trtkin,trttherm,trttot,trtpa, &
             r3,t3,p3,massflux,specenth,trtemth,trtkinth,trtthermth,trttotth, &
             trtpath
        ! the above four_vector arrays use a huge amount of memory with high res simulations because of 4x metric
        wjet=0; wwind=0
        u%data(1)=u0_arr; u%data(2)=vrl_arr; u%data(3)=vtl_arr; u%data(4)=vpl_arr
        b%data(1)=b0_arr; b%data(2)=br_arr; b%data(3)=bth_arr; b%data(4)=bph_arr
        call assign_metric(b,transpose(kerr_metric(r_arr,th_arr,asim)))
        call assign_metric(u,transpose(kerr_metric(r_arr,th_arr,asim)))
  
        r3=reshape(r_arr,(/nx1,nx2,nx3/))
        t3=reshape(th_arr,(/nx1,nx2,nx3/))
        p3=reshape(ph_arr,(/nx1,nx2,nx3/))

        ucov=lower(u); bcov=lower(b)

        if(magcrit.eq.1) then
           where(reshape(b*b/rho_arr, &
                (/nx1,nx2,nx3/)).ge.1)
              wjet=1
           endwhere
           where(reshape(b*b/rho_arr, &
                (/nx1,nx2,nx3/)).lt.1. & 
                .and.reshape(p_arr(1:nx1*nx2*nx3)/(b(1:nx1*nx2*nx3)*b(1:nx1*nx2*nx3)/2.), &
                (/nx1,nx2,nx3/)).lt.2.)
              wwind=1
           endwhere
        else
!           where(abs(cos(t3)).gt.0.98.and.reshape(ucov%data(1)*u%data(2),(/nx1,nx2,nx3/)).lt.0.)
!              wjet=1
!           endwhere
!           where(abs(cos(t3)).gt.0.866.and. &
!                abs(cos(t3)).lt.0.98.and.reshape(ucov%data(1)*u%data(2),(/nx1,nx2,nx3/)).lt.0.)
!              wwind=1
!           endwhere
           specenth=reshape(ucov%data(1)*(4.*p_arr+rho_arr)/rho_arr, &
                (/nx1,nx2,nx3/))
           where(specenth.lt.(-1.3).and.reshape(ucov%data(1)*u%data(2), &
                (/nx1,nx2,nx3/)).lt.0.)
              wjet=1
           endwhere
           where(specenth.lt.(-1.0).and.specenth.ge.(-1.3).and.reshape(ucov%data(1)*u%data(2), &
                (/nx1,nx2,nx3/)).lt.0.)
              wwind=1
           endwhere
        endif
           write(6,*) 'wjet, wwind: ',sum(wjet),sum(wwind)
        write(6,*) 'utest: ',maxval(abs(u*u+1.)), &
             maxval(abs(ucov%data(1)*u%data(1)+ucov%data(2)*u%data(2) &
             +u%data(3)*ucov%data(3)+u%data(4)*ucov%data(4)+1.))
        trtem=reshape(b*b*ucov%data(1)*u%data(2)-bcov%data(1)*b%data(2),(/nx1,nx2,nx3/))
        trtkin=reshape((ucov%data(1)+1.)*u%data(2)*rho_arr,(/nx1,nx2,nx3/))
        trtpa=reshape(rho_arr*ucov%data(1)*u%data(2),(/nx1,nx2,nx3/))
        trttherm=reshape(4.*p_arr*ucov%data(1)*u%data(2),(/nx1,nx2,nx3/))
        trttot=trtpa+trttherm+trtem
        massflux=reshape(u%data(2)*rho_arr,(/nx1,nx2,nx3/))
        wind_power=surf_integral(trtem*wwind,r3,t3,p3,asim)
        jet_power=surf_integral(trtem*wjet,r3,t3,p3,asim)
! for now kinetic, thermal terms include both wind & jet
        kinetic=surf_integral(trtkin*(wjet+wwind),r3,t3,p3,asim)
        therm=surf_integral(trttherm*(wjet+wwind),r3,t3,p3,asim)
        wind_kinetic=surf_integral(trtkin*wwind,r3,t3,p3,asim)
        wind_therm=surf_integral(trttherm*wwind,r3,t3,p3,asim)
        jet_magenergy=surf_integral(reshape((b*b)/2.,(/nx1,nx2,nx3/)) &
             *(wwind+wjet),r3,t3,p3,asim)
        jet_vel=surf_integral(reshape(u%data(2)/u%data(1),(/nx1,nx2,nx3/)) &
             *(wwind+wjet),r3,t3,p3,asim)/ &
             surf_integral(1d0*(wwind+wjet),r3,t3,p3,asim)
        jetmom=surf_integral(massflux*(wwind+wjet),r3,t3,p3,asim)
        jet_totpower=surf_integral(trttot*(wwind+wjet),r3,t3,p3,asim)
        write(6,*) 'gaminf: ',minval(trttot),maxval(trttot),minval(massflux),maxval(massflux), &
             minval(trttot/massflux),maxval(trttot/massflux)
        gaminf_az=surf_integral((-trttot/massflux)*(wwind+wjet),r3,t3,p3,asim,1,1)
        jetvel2d=surf_integral(reshape(u%data(2)/u%data(1),(/nx1,nx2,nx3/)),r3,t3,p3,asim,1,1)
        jetmom2d=surf_integral(massflux,r3,t3,p3,asim,1,1)
        totpower2d=surf_integral(trttot,r3,t3,p3,asim,1,1)
        therm2d=surf_integral(trttherm,r3,t3,p3,asim,1,1)
        jetpower2d=surf_integral(trtem,r3,t3,p3,asim,1,1)
        thermout2d=surf_integral((wjet+wwind)*trttherm,r3,t3,p3,asim,1,1)
! extra 2d \theta quantities:
        trtemth=reshape(b*b*ucov%data(1)*u%data(3)-bcov%data(1)*b%data(3),(/nx1,nx2,nx3/))
        trtkinth=reshape((ucov%data(1)+1.)*u%data(3)*rho_arr,(/nx1,nx2,nx3/))
        trtpath=reshape(rho_arr*ucov%data(1)*u%data(3),(/nx1,nx2,nx3/))
        trtthermth=reshape(4.*p_arr*ucov%data(1)*u%data(3),(/nx1,nx2,nx3/))
        trttotth=trtpath+trtthermth+trtemth
        totpower2dth=surf_integral(trttotth,r3,t3,p3,asim,1,1)
        therm2dth=surf_integral(trtthermth,r3,t3,p3,asim,1,1)
        jetpower2dth=surf_integral(trtemth,r3,t3,p3,asim,1,1)
        thermout2dth=surf_integral(wjet*trtthermth,r3,t3,p3,asim,1,1)
        write(6,*) 'trtem: ',sum(trtem)
        write(6,*) 'trtkin: ',sum(trtkin)
        write(6,*) 'trttherm: ',sum(trttherm)
        write(6,*) 'jet_vel: ',sum(jet_vel)
        write(6,*) 'jet_magenergy: ',sum(jet_magenergy)
        write(6,*) 'ucov: ',minval(ucov%data(1)), maxval(ucov%data(1))
        end subroutine calc_thickdisk_jet_power

!        subroutine calc_thickdisk_inner_edge(rin)
        ! calculate disk inner edge using definitions from hawley & krolik 2001
!        rin=
!        end subroutine calc_thickdisk_inner_edge

        subroutine calc_thickdisk_shellavg_properties(vr,ur,smri,pmag,pmagnowt,pgas,pgasnowt,beta,thetad,bz, &
             rhoint,vr2,alphamag,omega,bzcov,bphi,bphicov,pgas2d,pmag2d)
        real(kind=8), dimension(nx1), intent(out) :: vr, ur, pmag, &
             pgas, beta, thetad,  smri, bz, rhoint, vr2,alphamag,omega, &
             pmagnowt, pgasnowt, bzcov, bphi, bphicov
        real(kind=8), dimension(nx1) :: vtha, lammri
        real(kind=8), dimension(nx1,nx3) :: theta0,rhoth
        real(kind=8), dimension(nx1,nx2) :: pgas2d, pmag2d
        real(kind=8), dimension(nx1*nx2*nx3) :: theta03
        real(kind=8), dimension(nx1,nx2,nx3) :: theta3,r3,t3,p3,thetadd,lamdd
        integer :: i,j
        type (four_vector), dimension(nx1*nx2*nx3) :: ucov,bcov,u,b
        write(6,*) 'start of thickdisk shellavg'
        u%data(1)=u0_arr; u%data(2)=vrl_arr; u%data(3)=vtl_arr; u%data(4)=vpl_arr
        b%data(1)=b0_arr; b%data(2)=br_arr; b%data(3)=bth_arr; b%data(4)=bph_arr
        call assign_metric(b,transpose(kerr_metric(r_arr,th_arr,asim)))
        call assign_metric(u,transpose(kerr_metric(r_arr,th_arr,asim)))
        ucov=lower(u); bcov=lower(b)
! calculate some density weighted shell averages
        r3=reshape(r_arr,(/nx1,nx2,nx3/))
        t3=reshape(th_arr,(/nx1,nx2,nx3/))
        p3=reshape(ph_arr,(/nx1,nx2,nx3/))
        rhoint=surf_integral(reshape(1d0*rho_arr,(/nx1,nx2,nx3/)),r3,t3,p3,asim)
        write(6,*) 'shell rhoint', minval(rhoint), maxval(rhoint), rhoint(1:10), &
             rho_arr(1:10),1d0*rho_arr(1:10)
        rhoth=surf_integral(reshape(1d0*rho_arr,(/nx1,nx2,nx3/)),r3,t3,p3,asim,1)
        write(6,*) 'shell rhoth'
        ur=surf_integral(reshape(rho_arr*u%data(2),(/nx1,nx2,nx3/)),r3,t3,p3,asim)/rhoint
        vr=surf_integral(reshape(rho_arr*u%data(2)/u%data(1),(/nx1,nx2,nx3/)),r3,t3,p3,asim)/rhoint
        vr2=surf_integral(reshape(rho_arr*(u%data(2)/u%data(1))**2., &
             (/nx1,nx2,nx3/)),r3,t3,p3,asim)/rhoint
        write(6,*) 'shell ur', ur(1:10)
        omega=surf_integral(reshape(1d0*rho_arr*(u%data(4)/u%data(1)),(/nx1,nx2,nx3/)), &
             r3,t3,p3,asim)/rhoint
        write(6,*) 'shell omega'
        beta=surf_integral(reshape(1d0*rho_arr*p_arr/(b*b/2.),(/nx1,nx2,nx3/)),r3,t3, &
            p3,asim)/rhoint
        write(6,*) 'shell beta'
        pmag=surf_integral(reshape(rho_arr*(b*b)/2.,(/nx1,nx2,nx3/)),r3,t3, &
            p3,asim)/rhoint
        pmagnowt=surf_integral(reshape((b*b)/2.,(/nx1,nx2,nx3/)),r3,t3,p3,asim)
        pgasnowt=surf_integral(reshape(1d0*p_arr,(/nx1,nx2,nx3/)),r3,t3,p3,asim)
        pmag2d=surf_integral(reshape((b*b)/2.,(/nx1,nx2,nx3/)),r3,t3,p3,asim,1,1)
        pgas2d=surf_integral(reshape(1d0*p_arr,(/nx1,nx2,nx3/)),r3,t3,p3,asim,1,1)
        alphamag=surf_integral(reshape(rho_arr*1./2./p_arr*(u%data(2)*u%data(4)*(b*b) &
             -b%data(2)*b%data(4)),(/nx1,nx2,nx3/)),r3,t3, &
            p3,asim)/rhoint
        write(6,*) 'shell pmag'
        pgas=surf_integral(reshape(1d0*rho_arr*p_arr,(/nx1,nx2,nx3/)),r3,t3, &
             p3,asim)/rhoint
        write(6,*) 'shell pgas'
        bz=surf_integral(reshape(rho_arr*b%data(3),(/nx1,nx2,nx3/)),r3, &
             t3,p3,asim)/rhoint
        bzcov=surf_integral(reshape(rho_arr*bcov%data(3),(/nx1,nx2,nx3/)),r3, &
             t3,p3,asim)/rhoint
        bphicov=surf_integral(reshape(rho_arr*bcov%data(4),(/nx1,nx2,nx3/)),r3, &
             t3,p3,asim)/rhoint
        bphi=surf_integral(reshape(rho_arr*b%data(4),(/nx1,nx2,nx3/)),r3, &
             t3,p3,asim)/rhoint
        write(6,*) 'shell bz'
        theta0=pi/2.+surf_integral(reshape(rho_arr*(th_arr-pi/2.),(/nx1,nx2,nx3/)),r3, &
             t3,p3,asim,1)/rhoth
        write(6,*) 'shell theta0', size(theta0,1), size(theta0,2), size(theta3,1), size(theta3,3)
        do i=1,nx2
           theta3(:,i,:)=theta0
        enddo
        write(6,*) 'shell theta3'
        theta03=reshape(theta3,(/nx1*nx2*nx3/))
        write(6,*) 'shell theta03'!,size(sum(theta0,2))!, &
!             size(sum(theta0,2),2)!, shape(thetad),shape(sum(theta0,2))
        thetad=sum(sqrt(surf_integral(reshape(rho_arr*(th_arr-theta03)**2.,(/nx1,nx2,nx3/)), &
             r3,t3,p3,asim,1)/rhoth),2)/nx3
        write(6,*) 'shell thetad'
        vtha=surf_integral(reshape(rho_arr*sqrt(b%data(3)*bcov%data(3)/(b*b+rho_arr+4.*p_arr)),(/nx1,nx2,nx3/)), &
             r3,t3,p3,asim)/rhoint
        write(6,*) 'shell vtha'
        lammri=2.*pi*abs(vtha/omega)
        write(6,*) 'shell lammri'
        do i=1,nx2
           do j=1,nx3
              thetadd(:,i,j)=thetad
              lamdd(:,i,j)=lammri
           enddo
        enddo
        smri=surf_integral(reshape(rho_arr*2.*r_arr,(/nx1,nx2,nx3/))/lamdd*thetadd, &
             r3,t3,p3,asim)/rhoint
        write(6,*) 'shell smri'
        end subroutine calc_thickdisk_shellavg_properties

        subroutine initialize_thickdisk_model(a,transform,ifile,gf,df,ntt,nft,indft,jf, &
             off,simt,dindft,mct)
        real(kind=8), intent(in) :: a
        real(kind=8) :: tcur, tstep_test
        integer :: nx, status, writeall
        integer, intent(in) :: transform
        character(len=20), intent(in), optional :: ifile
        character(len=100), intent(in), optional :: simt,gf,df
        integer, intent(in), optional :: nft,indft,jf,off,ntt,dindft,mct
        character(len=20) :: default_ifile='thickdisk.in', append
        character(len=file_len) :: header_file
        writeall=0
        if (present(gf)) then
           gfile = gf
           dfile = df
           nt = ntt
           nfiles = nft
           indf = indft
           jonfix = jf
        else
           if (present(ifile)) then
              call read_thickdisk_inputs(ifile)
           else
              call read_thickdisk_inputs(default_ifile)
           endif
        endif
        write(6,*) 'inputs read: ',gfile,dfile,nt,nfiles,indf,jonfix
!        dfile='m87bl09rfp10xi5a998fluidvars.bin'
!        write(6,*) 'init thickdisk'
! length of header at start of thickdisk files (these should be determined by Python)
!        header_length=390+offset
!        field_header_length=header_length+6
! should change fmt to be more general -- figure out max # files, assign fmt accordingly
        write(append, fmt='(I4.4)') indf
        header_file=trim(dfile) // trim(append) // trim('.bin')
        write(6,*) 'nhead: ',nhead
        call read_thickdisk_data_header(header_file,header_length)
!        call read_thickdisk_data_header(header_file,field_header_length)
        write(6,*) 'header lengths: ',header_length,field_header_length
        write(6,*) 'rin: ',rin,rout,calcrmks(10d0,xbr)
        write(6,*) 'nhead: ',nhead,allocated(ph_arr)
        if(abs(asim-a).gt.1e-4) then 
           write(6,*) 'ERROR -- Different simulation and grtrans spin values!', asim, a
           return
        endif
        n=nx1*nx2*nx3
        call init_thickdisk_data(n,n*nt)
! add option to mb09 style or not depending on sim name
        if(trim(gfile)=='thickdisk7grid.bin') then
! need to have xbr, rout set here
           call read_mb09_grid_file(nt)
           if(rout.gt.1.e3) then
              xbr=log(500.)
           else 
              xbr=log(1.e5)
           endif
           r_arr = calcrmks(x1_arr,xbr)
           th_arr = calcthmks(x2_arr,calcrmks(x1_arr,xbr))
           ph_arr = x3_arr*2.*pi
        else
           call read_thickdisk_data_header(gfile,header_length)
           write(6,*) 'updated header_length for gfile: ',header_length
           call read_thickdisk_grid_file(gfile,tcur,binary)
        endif
        write(6,*) 'nhead: ',nhead,allocated(ph_arr)
        write(6,*) 'read header', nx1, nx2, allocated(ph_arr), allocated(r_arr)
        write(6,*) 'init data', nt, indf, dfile
        call load_thickdisk_data(nt,transform)
        if(writeall==1) then
! write all BL coordinate data to a file
           open(unit=10,file='thickdisk3ddata.out',form='unformatted')
           write(10) r_arr
           write(10) th_arr
           write(10) ph_arr
           write(10) rho_arr
           write(10) p_arr
           write(10) b0_arr
           write(10) br_arr
           write(10) bth_arr
           write(10) bph_arr
           write(10) u0_arr
           write(10) vrl_arr
           write(10) vtl_arr
           write(10) vpl_arr
           close(10)
        endif
        end subroutine initialize_thickdisk_model

        subroutine load_thickdisk_data(nt,transform)
        real(kind=8), dimension(:), allocatable :: rho,p
        real, dimension(:), allocatable :: vrl, vtl, vpl
        type (four_vector), dimension(:), allocatable :: u, b
        integer, intent(in) :: nt, transform
        character(len=file_len) :: data_file
        character(len=20) :: append
        integer :: k
        real(kind=8) :: tcur, tstep_test
        allocate(rho(n)); allocate(p(n))
        allocate(vrl(n)); allocate(vtl(n)); allocate(vpl(n))
        allocate(u(n)); allocate(b(n))
        do k=1,nt
           write(append, fmt='(I4.4)') indf-(k-1)
           data_file = trim(dfile) // trim(append) // '.bin'
!           write(6,*) 'data_file: ',indf-(k-1),append,data_file
           call read_thickdisk_data_header(data_file,field_header_length)
           call read_thickdisk_fieldline_file(data_file,tcur,rho,p,u,b,jonfix,binary,test)
           if(k.eq.1) t(k)=tcur
           ! transform velocities to LNRF frame
!           write(6,*) 'lnrf sizes: ',size(u%data(2)), size(r_arr), size(th_arr), size(vrl), size(vtl)
!           write(6,*) 'lnrf transform', size(b0_arr), size(b%data(1)), size(b0_arr((k-1)*n+1:k*n))
!           write(6,*) 'lnrf transform', size(vrl), size(vtl), size(vpl), size(u0_arr)
           ! now assign to data arrays
           rho_arr((k-1)*n+1:k*n)=rho
! conversion between internal energy and pressure is done in fieldline read
           p_arr((k-1)*n+1:k*n)=p
           b0_arr((k-1)*n+1:k*n)=b%data(1)
           br_arr((k-1)*n+1:k*n)=b%data(2)
           bth_arr((k-1)*n+1:k*n)=b%data(3)
           bph_arr((k-1)*n+1:k*n)=b%data(4)
           write(6,*) 'fluid model thickdisk v assign: ',k,n,(k-1)*n+1,k*n,size(u%data(1))
           u0_arr((k-1)*n+1:k*n)=u%data(1)
           write(6,*) 'fluid model thickdisk u0_arr min max udotu: ', minval(u0_arr),maxval(u0_arr),maxval(abs(u*u+1d0))
!           write(6,*) 'vrl assign'
! for some applications want lnrf vels rather than four-velocity
           if(transform.eq.1) then
              call lnrf_frame(real(u%data(2)/u%data(1)),real(u%data(3)/u%data(1)), & 
                real(u%data(4)/u%data(1)),real(r_arr),real(asim),real(th_arr),vrl,vtl,vpl)
              vrl_arr((k-1)*n+1:k*n)=real(vrl)
              !           write(6,*) 'after assign', size(vtl_arr((k-1)*n+1:k*n)), size(vtl)
              vpl_arr((k-1)*n+1:k*n)=real(vpl)
              !           write(6,*) 'after vpl', vtl(1), vtl(n), vtl
              vtl_arr((k-1)*n+1:k*n)=real(vtl)
              !           write(6,*) 'after vtl'
              !           write(6,*) 'assign'
           else
              vrl_arr((k-1)*n+1:k*n)=real(u%data(2))
              vtl_arr((k-1)*n+1:k*n)=real(u%data(3))
              vpl_arr((k-1)*n+1:k*n)=real(u%data(4))
           endif
        end do
        tstep=2.
        if (nt.gt.1) then
           tstep=t(1)-t(2)
           tstep_test=t(nt-1)-t(nt)
           if (abs((tstep-tstep_test)/tstep).gt.0.1) then 
              write(6,*) 'WARNING -- Time step changes over dumps!'
           endif
        endif
!        write(6,*) 'after loop', maxval(abs(u*u+1.))
        deallocate(rho); deallocate(p); deallocate(u); deallocate(b)
        deallocate(vrl); deallocate(vtl); deallocate(vpl)
!        write(6,*) 'read data file', a
        write(6,*) 'fluid model thickdisk vmax: ',maxval(vrl_arr**2.+vtl_arr**2.+vpl_arr**2.)
!        write(6,*) 'maxmin temp: ',minval(p_arr/rho_arr*1.67e-24*9e20/1.38e-16), &
!             maxval(p_arr/rho_arr*1.67e-24*9e20/1.38e-16)
        end subroutine load_thickdisk_data

        subroutine advance_thickdisk_timestep(dtobs)
        real(kind=8), intent(in) :: dtobs
        integer :: nupdate
        nupdate = floor(dtobs / tstep)
        toffset = toffset+dtobs-tstep*nupdate
! after a couple stpes might be accumulating offset
        nupdate=nupdate+floor(toffset/tstep)
        toffset=toffset-floor(toffset/tstep)*tstep
        indf=indf+nupdate
        write(6,*) 'advance thickdisk timestep: ',dtobs,tstep,toffset,indf,nupdate
        call update_thickdisk_data(nupdate)
        end subroutine advance_thickdisk_timestep

        subroutine update_thickdisk_data(nupdate)
        integer, intent(in) :: nupdate
        integer :: nshift
        if(nupdate.gt.nt) then
! this is case where we just load all new data
           call load_thickdisk_data(nt,1)
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
           call load_thickdisk_data(nupdate,1)
        end if
        end subroutine update_thickdisk_data
 
        subroutine init_thickdisk_data(nx,n)
        integer, intent(in) :: n,nx
        allocate(rho_arr(n)); allocate(p_arr(n)); allocate(u0_arr(n))
        allocate(vrl_arr(n)); allocate(vtl_arr(n)); allocate(vpl_arr(n))
        allocate(b0_arr(n)); allocate(br_arr(n)); allocate(bth_arr(n))
        allocate(bph_arr(n)); allocate(ph_arr(nx))
        allocate(x1_arr(nx)); allocate(x2_arr(nx)); allocate(r_arr(nx))
        allocate(th_arr(nx)); allocate(t(nt)); allocate(x3_arr(nx))
        end subroutine init_thickdisk_data

        subroutine del_thickdisk_data()
        deallocate(rho_arr); deallocate(p_arr); deallocate(u0_arr)
        deallocate(vrl_arr); deallocate(vtl_arr); deallocate(vpl_arr)
        deallocate(b0_arr); deallocate(br_arr); deallocate(bth_arr)
        deallocate(bph_arr); deallocate(x1_arr); deallocate(x2_arr)
        deallocate(r_arr); deallocate(th_arr); deallocate(ph_arr)
        deallocate(t); deallocate(x3_arr)
        end subroutine del_thickdisk_data

      end module fluid_model_thickdisk
