   module fluid_model_mb09

      use class_four_vector
      use interpolate, only: interp
      use kerr, only: kerr_metric, lnrf_frame, uks2ubl, bl2ks, surf_integral
      use phys_constants, only: pi
      use math, only: zbrent
      implicit none

      namelist /mb09/  dfile, gfile, nt, nfiles, indf, jonfix, offset, sim
      integer :: nx1, nx2, nx3, ndumps, nt, nfiles, indf, defcoord, &
           jonfix, dlen, offset, header_length, field_header_length
      character(len=100) :: dfile, gfile
      character(len=100) :: sim
      integer :: n
      integer :: ndim=3, nhead=30, test=1, binary=1, DEBUG=0
      integer, dimension(:), allocatable :: dumps
      real :: tstep, h, dx1, dx2, dx3, gam, startx1, startx2, startx3, dt, &
           r0, rin, rout
      real(8) :: xbr,toffset=0.,asim
      real, dimension(:), allocatable :: t
      real(8), dimension(:), allocatable :: x1_arr, x2_arr, x3_arr, r_arr, th_arr, ph_arr
      real, dimension(:), allocatable :: rho_arr, p_arr, u0_arr, vrl_arr, &
        vpl_arr, vtl_arr
      real, dimension(:), allocatable :: b0_arr, br_arr, bth_arr, bph_arr

      interface init_mb09_data
        module procedure init_mb09_data
      end interface

      interface calcthmks
         module procedure calcthmks
         module procedure calcthmks_single
         module procedure calcthmks8
         module procedure calcthmks8_single
      end interface
      
      interface del_mb09_data
        module procedure del_mb09_data
      end interface

      interface initialize_mb09_model
        module procedure initialize_mb09_model
      end interface

      interface transformbl2mks
         module procedure transformbl2mks
      end interface

      interface umks2uks
         module procedure umks2uks
      end interface

      interface mb09_vals
        module procedure mb09_vals
      end interface

!      interface coord_def
!         module procedure coord_def
!      end interface

      interface calcrmks
         module procedure calcrmks
         module procedure calcrmks_single
      end interface

      contains

        function calcrmks(x1,xbr) result(r)
          ! Compute r given x1 for Jon's simulations
          ! JAD 7/24/2011
          real(8), intent(in), dimension(:) :: x1
          real(8), intent(in) :: xbr
          real(8) :: npow2
          real(8), dimension(size(x1)) :: r, xi
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
          real(8), intent(in) :: x1,xbr
          real(8) :: npow2
          real(8) :: r, xi
          npow2=10d0
          if(x1.gt.xbr) then
             xi=x1+(x1-xbr)**npow2
          else
             xi=x1
          endif
          r=exp(xi)
        end function calcrmks_single

        function umks2uks(fum,x1,x2,xbr) result(uks)
          real(8), intent(in), dimension(:) :: x1,x2
          type (four_vector), dimension(:), intent(in) :: fum
          real(8), intent(in) :: xbr
          type (four_vector), dimension(size(fum)) :: uks
          real(8), dimension(size(x1)) :: r,dr,dx2,dx1,drdx1,dthdrnum,dthdx2num
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
          write(6,*) 'umks2uks r: ',r(1),x1(1),x2(1),xbr
          write(6,*) 'umks2uks dr: ',dr(1),dx1(1),dx2(1)
          write(6,*) 'umks2uks drdx: ',calcrmks(x1(1:2)+0.5d0*dx1(1:2),xbr)
          write(6,*) 'umks2uks umks: ',fum(1)%data(1:4)
          write(6,*) 'umks2uks drdx1: ',drdx1(1)
          write(6,*) 'umks2uks dthdrnum: ',dthdrnum(1)
          write(6,*) 'umks2uks dthdx2num: ',dthdx2num(1)
          write(6,*) 'umks2uks uks: ',uks(1)%data
        end function umks2uks

        function calcthmks(x2,r) result(th)
!          Calculate theta from mks coordinate x2 (McKinney & Gammie (2006a)) and bl/ks coordinate r.
!          JAD 4/9/2009
!          These parameter definitions are based on Jon's defcoord=9:
          real, dimension(:), intent(in) :: x2,r
          real, dimension(size(x2)) :: th
          real, dimension(size(r)) :: g,h
          real :: rj,nj,r0j,rsj,q
          rj=2.8
          nj=0.3
          r0j=20.
          rsj=80.
          q=1.3
          g=-nj*(.5+1./pi*atan((r-rsj)/r0j))
          h=2.-q*(r/rj)**g
          where(x2.lt.0.5)
             th=(pi*x2+(1.-h)/2.*sin(2.*pi*x2))
          elsewhere
             th=(pi*x2-(1.-h)/2.*sin(2.*pi*(1.-x2)))
          endwhere
        end function calcthmks

        function calcthmks_single(x2,r) result(th)
!          Calculate theta from mks coordinate x2 (McKinney & Gammie (2006a)) and bl/ks coordinate r.
!          JAD 4/9/2009
!          These parameter definitions are based on Jon's defcoord=9:
          real, intent(in) :: x2,r
          real :: th
          real :: rj,nj,r0j,rsj,q,g,h
          rj=2.8
          nj=0.3
          r0j=20.
          rsj=80.
          q=1.3
          g=-nj*(.5+1./pi*atan((r-rsj)/r0j))
          h=2.-q*(r/rj)**g
          if(x2.lt.0.5) then
             th=(pi*x2+(1.-h)/2.*sin(2.*pi*x2))
          else
             th=(pi*x2-(1.-h)/2.*sin(2.*pi*(1.-x2)))
          endif
        end function calcthmks_single

        function calcthmks8(x2,r) result(th)
!          Calculate theta from mks coordinate x2 (McKinney & Gammie (2006a)) and bl/ks coordinate r.
!          JAD 4/9/2009
!          These parameter definitions are based on Jon's defcoord=9:
          real(8), dimension(:), intent(in) :: x2,r
          real(8), dimension(size(x2)) :: th
          real(8), dimension(size(r)) :: g,h
          real(8) :: rj,nj,r0j,rsj,q
          rj=2.8
          nj=0.3
          r0j=20.
          rsj=80.
          q=1.3
          g=-nj*(.5+1./pi*atan((r-rsj)/r0j))
          h=2.-q*(r/rj)**g
          where(x2.lt.0.5)
             th=(pi*x2+(1.-h)/2.*sin(2.*pi*x2))
          elsewhere
             th=(pi*x2-(1.-h)/2.*sin(2.*pi*(1.-x2)))
          endwhere
        end function calcthmks8

        function calcthmks8_single(x2,r) result(th)
!          Calculate theta from mks coordinate x2 (McKinney & Gammie (2006a)) and bl/ks coordinate r.
!          JAD 4/9/2009
!          These parameter definitions are based on Jon's defcoord=9:
          real(8), intent(in) :: x2,r
          real(8) :: th
          real(8) :: rj,nj,r0j,rsj,q,g,h
          rj=2.8
          nj=0.3
          r0j=20.
          rsj=80.
          q=1.3
          g=-nj*(.5+1./pi*atan((r-rsj)/r0j))
          h=2.-q*(r/rj)**g
          if(x2.lt.0.5) then
             th=(pi*x2+(1.-h)/2.*sin(2.*pi*x2))
          else
             th=(pi*x2-(1.-h)/2.*sin(2.*pi*(1.-x2)))
          endif
        end function calcthmks8_single


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
          real(8), intent(in), dimension(:) :: x2, r
          real(8), dimension(size(x2)) :: theta,myhslope,th2, &
               th0,switch0,switch2,theta1,theta2,arctan2
          real(8) :: r0r,r1jet,njet,r0jet,rsjet,qjet, &
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
          real(8), intent(in) :: x2,r
          real(8) :: th,r0r,r1jet,njet,r0jet,rsjet,qjet, &
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
          real(8), intent(in), dimension(:) :: args
          real(8), intent(in) :: x1
          real(8) :: xbr,npow2,r,xi,diff
        ! Compute x1 given r for Jon's simulations 
          xbr=args(2)
          npow2=10d0
          r=args(1)
          diff=r-calcrmks(x1,xbr)
        end function findx1mks

        subroutine transformbl2mks(r,th,ph,x1,x2,x3)
        ! transform Boyer-Lindquist coordinates to modified Kerr-Schild coordinates used by mb09
        real, intent(in), dimension(:) :: r,th,ph
        real, intent(out), dimension(size(r)) :: x1,x2,x3
        real, dimension(2) :: args
        integer :: i
        do i=1,size(r)
           args(1)=r(i); args(2)=xbr
           x1(i)=zbrent(findx1mks,0d0,10d0,dble(args),1d-6)
           args(2)=th(i)
           x2(i)=zbrent(findx2mks,0d0,1d0,dble(args),1d-6)
        enddo
        x3=ph/2./pi
        end subroutine transformbl2mks

        subroutine mb09_vals(x0,a,rho,p,b,u,bmag)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real(8), dimension(size(x0)) :: done
        real, dimension(size(x0)) :: x2,x1,x3,zm,theta,fone, &
         vpl0,vrl0,vtl0,rd,td,pd,rttd,zr,dzero, &
         vr0,vth0,vph0,bth,dummy,tt,ttd,zt,zphi,zpp
        real(8), dimension(nx1) :: uniqx1,uniqr
        real(8), dimension(nx2) :: uniqx2
        real(8), dimension(nx3) :: uniqx3,uniqp
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
        ! Interpolates MB09 data to input coordinates
        ! JAD 3/20/2009, fortran 12/28/2012
!        write(6,*) 'mb09_vals'
        one=1; dzero=0d0; done=1d0; fone=1.
!        write(6,*) 'before uniq: ',nx1,nx2,nx3
        uniqx1=x1_arr(1:nx1)
!        write(6,*) 'between', minval(x0%data(4)), maxval(x0%data(4))
        uniqx2=x2_arr(1:nx1*(nx2-1)+1:nx1)
!        write(6,*) 'uniqx3: ',size(x3_arr(1:nx1*nx2*(nx3-1)+1:nx1*nx2)),nx3
!        write(6,*) 'uniqx3: ',size(x3_arr),nx1*nx2*nx3,nx1*nx2*(nx3-1)+1
        uniqx3=x3_arr(1:nx1*nx2*(nx3-1)+1:nx1*nx2)
! put case statement for coordinate systems here? Then do all transformations?
 !       write(6,*) 'after uniqx3',minval(uniqx1),maxval(uniqx1),xbr
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
! put mb09 on [0,2pi]
        zphi=mod(zphi,(2.*pi))
        where(zphi.lt.0.)
           zphi=zphi+2.*pi
        endwhere
!        where(zphi.gt.(2.*pi))
!           zphi=zphi-2.*pi
!        endwhere
!        write(6,*) 'zphi: ',minval(zphi),maxval(zphi)
!        write(6,*) 'transform: ',xbr,minval(zr),maxval(zr)
        call transformbl2mks(zr,theta,zphi,x1,x2,x3)
!        write(6,*) 'transform: ',maxval(zr), minval(zr), maxval(x1), minval(x1)
!        write(6,*) 'transform: ',maxval(theta), minval(theta), minval(x2), maxval(x2)
        lx1=floor((x1-uniqx1(1))/(uniqx1(nx1)-uniqx1(1))*(nx1-1))+1
        ux1=lx1+1
        lx2=floor((x2-uniqx2(1))/(uniqx2(nx2)-uniqx2(1))*(nx2-1))+1
        ux2=lx2+1
        lx3=floor((x3-uniqx3(1))/(uniqx3(nx3)-uniqx3(1))*(nx3-1))+1
        ux3=lx3+1
        lx2=merge(lx2,one,lx2.ge.1)
        ux2=merge(ux2,nx2,ux2.le.nx2)
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
!        write(6,*) '7041: ',u(7041)*u(7041)
        if(DEBUG.eq.1) then
        ! output quantities for debugging
           write(6,*) 'debug!'
           open(unit=8,file='fluid_model_mb09_debug.out')
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
           write(8,*) u0i, vrli, vtli, vpli
           write(8,*) u%data(1), u%data(2), u%data(3), u%data(4)
        endif
!        write(6,*) 'mb09 vals rho: ',rho
!        write(6,*) 'mb09 vals minloc: ',minloc(rho),rd(minloc(rho)), &
!             td(minloc(rho)),x0(minloc(rho))%data(3),uniqth(lx2(minloc(rho))), &
!             lx2(minloc(rho))
!        write(6,*) 'mb09 vals rhoi: ',rhoi(minloc(rho),:)
        end subroutine mb09_vals

        subroutine read_mb09_inputs(ifile)
        character(len=20), intent(in) :: ifile
        open(unit=8,file=ifile,form='formatted',status='old')
        read(8,nml=mb09)
        write(6,*) 'read: ',nt
        close(unit=8)
        xbr=25.
        end subroutine read_mb09_inputs

        subroutine calc_mb09_jet_power(jet_power,wind_power)
        ! calculate jet and wind power using definitions from mckinney et al. 2012
        ! this assumes data loaded is four-velocity rather than lnrf vels
        real(8), dimension(nx1), intent(out) :: jet_power, wind_power
        integer, dimension(nx1,nx2,nx3) :: wjet, wwind
        type (four_vector), dimension(nx1*nx2*nx3) :: ucov,bcov,u,b
        real(8), dimension(nx1,nx2,nx3) :: trtem
        wjet=0; wwind=0
        u%data(1)=u0_arr; u%data(2)=vrl_arr; u%data(3)=vtl_arr; u%data(4)=vpl_arr
        b%data(1)=b0_arr; b%data(2)=br_arr; b%data(3)=bth_arr; b%data(4)=bph_arr
        call assign_metric(b,transpose(kerr_metric(r_arr,th_arr,asim)))
        call assign_metric(u,transpose(kerr_metric(r_arr,th_arr,asim)))
        where(reshape(b(1:nx1*nx2*nx3)*b(1:nx1*nx2*nx3)/rho_arr(1:nx1*nx2*nx3), &
             (/nx1,nx2,nx3/)).gt.1)
           wjet=1
        endwhere
        where(reshape(b(1:nx1*nx2*nx3)*b(1:nx1*nx2*nx3)/rho_arr(1:nx1*nx2*nx3), &
             (/nx1,nx2,nx3/)).lt.1. & 
             .and.reshape(p_arr(1:nx1*nx2*nx3)/(b(1:nx1*nx2*nx3)*b(1:nx1*nx2*nx3)/8./pi), &
             (/nx1,nx2,nx3/)).lt.2.)
           wwind=1
        endwhere
        write(6,*) 'wjet, wwind: ',sum(wjet),sum(wwind)
        ucov=lower(u); bcov=lower(b)
        trtem=reshape(b*b*ucov%data(1)*u%data(2)-bcov%data(1)*b%data(2),(/nx1,nx2,nx3/))
        wind_power=surf_integral(trtem*wwind,reshape(r_arr,(/nx1,nx2,nx3/)),reshape(th_arr,(/nx1,nx2,nx3/)),&
             reshape(ph_arr,(/nx1,nx2,nx3/)),asim)
        jet_power=surf_integral(trtem*wjet,reshape(r_arr,(/nx1,nx2,nx3/)),reshape(th_arr,(/nx1,nx2,nx3/)),&
             reshape(ph_arr,(/nx1,nx2,nx3/)),asim)
        end subroutine calc_mb09_jet_power

!        subroutine calc_mb09_inner_edge(rin)
        ! calculate disk inner edge using definitions from hawley & krolik 2001
!        rin=
!        end subroutine calc_mb09_inner_edge

        subroutine initialize_mb09_model(a,transform,ifile,gf,df,ntt,nft, &
             indft,jf,simt)
        real(8), intent(in) :: a
        real(8) :: tcur, tstep_test
        integer :: nx, status
        integer, intent(in) :: transform
        character(len=20), intent(in), optional :: ifile,simt
        character(len=20) :: default_ifile='mb09.in', append
        character(len=100), intent(in), optional :: gf,df
        integer, intent(in), optional :: ntt,nft,indft,jf
        character(len=100) :: header_file
        if (present(gf)) then
           dfile = df
           gfile = gf
           nt = ntt
           nfiles = nft
           indf = indft
           jonfix = jf
           sim = simt
        else
           if (present(ifile)) then
              call read_mb09_inputs(ifile)
           else
              call read_mb09_inputs(default_ifile)
           endif
        endif
! hard coded to be large and not be used.
        xbr = 25.
        write(6,*) 'inputs read: ',allocated(ph_arr),a
!        write(6,*) 'init mb09'
! length of header at start of mb09 files (these should be determined by Python)
!        header_length=390+offset
!        field_header_length=header_length+6
! should change fmt to be more general -- figure out max # files, assign fmt accordingly
        write(append, fmt='(I4.4)') indf
!        header_file=trim(dfile) // trim(append) // trim('.bin')
        write(6,*) 'nhead: ',nhead
!        call read_mb09_data_header(gfile,header_length)
!        call read_mb09_data_header(header_file,field_header_length)
!        write(6,*) 'header lengths: ',header_length,field_header_length
!        write(6,*) 'rin: ',rin,rout,calcrmks(10d0,xbr)
        write(6,*) 'nhead: ',nhead,allocated(ph_arr)
!        if(abs(asim-a).gt.1e-4) then 
!           write(6,*) 'ERROR -- Different simulation and grtrans spin values!', asim, a
!           return
!        endif
        asim=a
        n=nx1*nx2*nx3
!        call init_mb09_data(n,n*nt)
        write(6,*) 'init: ',allocated(rho_arr),size(rho_arr),n,nt
        call read_mb09_grid_file(nt)
        write(6,*) 'read grid: ',minval(x1_arr),maxval(x1_arr)
! transform to get others
        ph_arr=2.*pi*x3_arr
        r_arr=calcrmks(x1_arr,xbr)
        th_arr=calcthmks(x2_arr,r_arr)
        write(6,*) 'r th ph: ',minval(r_arr),maxval(th_arr)
        write(6,*) 'nhead: ',nhead,allocated(ph_arr)
        write(6,*) 'read header', nx1, nx2, allocated(ph_arr), allocated(r_arr)
        write(6,*) 'init data', nt, indf, dfile
        call load_mb09_data(nt,1)
        end subroutine initialize_mb09_model

        subroutine read_mb09_grid_file(nt)
        integer, intent(in) :: nt
        integer :: n
          open(unit=8,file=gfile,form='unformatted',status='old')
          read(8) nx1,nx2,nx3
          n=nx1*nx2*nx3
          ! allocate data
          call init_mb09_data(n,n*nt)
          write(6,*) 'read grid: ',size(x1_arr),size(x2_arr),size(x3_arr)
          read(8) x1_arr
          read(8) x2_arr
          read(8) x3_arr
          close(unit=8)
        end subroutine read_mb09_grid_file

        subroutine read_mb09_data(dfile,rho,p,u0,vr,vth,vph,b0,br,bth,bph)
        integer :: nx, status
        real, dimension(:), allocatable :: data
        real, intent(inout), dimension(:) :: rho,p,u0,vr,vth,vph,b0,br,bth,bph
        real, dimension(size(rho),10) :: metric
        real, dimension(size(rho)) :: ui2
        character(len=100), intent(in) :: dfile
        open(unit=8,file=dfile,form='unformatted',status='old')
        read(8) nx
        if(nx.ne.9*nx1*nx2*nx3) then
           write(6,*) 'ERROR: INCORRECT DATA SIZE IN READ_MB09_DATA: ',nx,nx1,nx2,nx3
           return
        endif
        allocate(data(nx))
        read(8) data
        close(8)
        rho=data(1:nx1*nx2*nx3)
        p=data(nx1*nx2*nx3+1:2*nx1*nx2*nx3)
        vr=data(2*nx1*nx2*nx3+1:3*nx1*nx2*nx3)
        vth=data(3*nx1*nx2*nx3+1:4*nx1*nx2*nx3)
        vph=data(4*nx1*nx2*nx3+1:5*nx1*nx2*nx3)
        b0=data(5*nx1*nx2*nx3+1:6*nx1*nx2*nx3)
        br=data(6*nx1*nx2*nx3+1:7*nx1*nx2*nx3)
        bth=data(7*nx1*nx2*nx3+1:8*nx1*nx2*nx3)
        bph=data(8*nx1*nx2*nx3+1:9*nx1*nx2*nx3)
        deallocate(data)
! calculate u0 from others (from IDL grtrans):
        metric=kerr_metric(r_arr,th_arr,asim)
        ui2=metric(:,1)+2.*metric(:,4)*vph+vr**2*metric(:,5)+metric(:,8)*vth**2+metric(:,10)*vph**2
        u0=1./sqrt(-ui2)
        end subroutine read_mb09_data

! routine to load a single simh52binary.pro binary file
        subroutine load_mb09_data(nt,transform)
        integer :: nx, k
        integer, intent(in) :: nt,transform
        real, dimension(:), allocatable :: rho,p,u0,vr,vth,vph,b0,br,bth,bph,vrl,vtl,vpl
        real :: tcur
!        character(len=100), intent(in) :: dfile
        character(len=100) :: data_file
        character(len=20) :: append
        tstep=2.
        write(6,*) 'load mb09', dfile
        n=nx1*nx2*nx3
        allocate(rho(n)); allocate(p(n)); allocate(u0(n))
        allocate(vph(n)); allocate(vth(n)); allocate(b0(n))
        allocate(br(n)); allocate(bth(n)); allocate(bph(n))
        allocate(vr(n))
        if(transform.eq.1) then
           allocate(vrl(n)); allocate(vtl(n)); allocate(vpl(n))
        endif
        do k=1,nt
           write(append, fmt='(I4.4)') indf-(k-1)
           data_file = trim(dfile) // trim(append) // '.bin'
           write(6,*) 'data_file: ',indf-(k-1),append,data_file
           call read_mb09_data(data_file,rho,p,u0,vr,vth,vph,b0,br,bth,bph)
           write(6,*) 'read mb09 data: ',minval(u0),maxval(u0),minval(vph),maxval(vph)
           write(6,*) 'read mb09 data: ',minval(vr),maxval(vr),minval(vth),maxval(vth)
           write(6,*) 'read mb09 data: ',asim,minval(r_arr),maxval(r_arr),minval(th_arr),maxval(th_arr)
           write(6,*) 'read mb09 data: ',minval(b0),maxval(b0),nx1,nx2,nx3
           write(6,*) 'after read mb09 data', n, size(rho), size(rho_arr)
           write(6,*) 'after read mb09 data', minval(vph), maxval(vph), minval(vr), maxval(vr)
           if(k.eq.1) t(k)=tcur
           ! assign to data arrays
           rho_arr((k-1)*n+1:k*n)=rho
           write(6,*) 'after assign rho'
           p_arr((k-1)*n+1:k*n)=p
           b0_arr((k-1)*n+1:k*n)=b0
           br_arr((k-1)*n+1:k*n)=br
           bth_arr((k-1)*n+1:k*n)=bth
           bph_arr((k-1)*n+1:k*n)=bph
!           write(6,*) 'v assign'
           u0_arr((k-1)*n+1:k*n)=u0
           write(6,*) 'after all assign',minval(p_arr),maxval(rho_arr)
! for some applications want lnrf vels rather than four-velocity
           if(transform.eq.1) then
              call lnrf_frame(vr,vth,vph,real(r_arr),real(asim),real(th_arr),vrl,vtl,vpl)
              vrl_arr((k-1)*n+1:k*n)=real(vrl)
              !           write(6,*) 'after assign', size(vtl_arr((k-1)*n+1:k*n)), size(vtl)
              vpl_arr((k-1)*n+1:k*n)=real(vpl)
              !           write(6,*) 'after vpl', vtl(1), vtl(n), vtl
              vtl_arr((k-1)*n+1:k*n)=real(vtl)
                         write(6,*) 'after vtl',minval(vrl_arr),maxval(vrl_arr),minval(vpl_arr)
              !           write(6,*) 'assign'
           else
              vrl_arr((k-1)*n+1:k*n)=vr
              vpl_arr((k-1)*n+1:k*n)=vph
              vtl_arr((k-1)*n+1:k*n)=vth
           endif
        end do
        deallocate(rho); deallocate(p); deallocate(u0)
        deallocate(vph); deallocate(vth); deallocate(b0)
        deallocate(br); deallocate(bth); deallocate(bph)
        deallocate(vr)
        if(transform.eq.1) then
           deallocate(vrl); deallocate(vtl); deallocate(vpl)
        endif
           
        end subroutine load_mb09_data

        subroutine advance_mb09_timestep(dtobs)
        real(8), intent(in) :: dtobs
        integer :: nupdate
        nupdate = floor(dtobs / tstep)
        toffset = toffset+dtobs-tstep*nupdate
! after a couple stpes might be accumulating offset
        nupdate=nupdate+floor(toffset/tstep)
        toffset=toffset-floor(toffset/tstep)*tstep
        indf=indf+nupdate
        write(6,*) 'advance mb09 timestep: ',dtobs,tstep,toffset,indf,nupdate
        call update_mb09_data(nupdate)
        end subroutine advance_mb09_timestep

        subroutine update_mb09_data(nupdate)
        integer, intent(in) :: nupdate
        integer :: nshift
        if(nupdate.gt.nt) then
! this is case where we just load all new data
           call load_mb09_data(nt,1)
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
           call load_mb09_data(nupdate,1)
        end if
        end subroutine update_mb09_data
 
        subroutine init_mb09_data(nx,n)
        integer, intent(in) :: n,nx
        allocate(rho_arr(n)); allocate(p_arr(n)); allocate(u0_arr(n))
        allocate(vrl_arr(n)); allocate(vtl_arr(n)); allocate(vpl_arr(n))
        allocate(b0_arr(n)); allocate(br_arr(n)); allocate(bth_arr(n))
        allocate(bph_arr(n)); allocate(ph_arr(nx))
        allocate(x1_arr(nx)); allocate(x2_arr(nx)); allocate(r_arr(nx))
        allocate(th_arr(nx)); allocate(t(nt)); allocate(x3_arr(nx))
        end subroutine init_mb09_data

        subroutine del_mb09_data()
        deallocate(rho_arr); deallocate(p_arr); deallocate(u0_arr)
        deallocate(vrl_arr); deallocate(vtl_arr); deallocate(vpl_arr)
        deallocate(b0_arr); deallocate(br_arr); deallocate(bth_arr)
        deallocate(bph_arr); deallocate(x1_arr); deallocate(x2_arr)
        deallocate(r_arr); deallocate(th_arr); deallocate(ph_arr)
        deallocate(t); deallocate(x3_arr)
        end subroutine del_mb09_data

      end module fluid_model_mb09
