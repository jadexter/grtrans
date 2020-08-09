   module fluid_model_harmpi

      use class_four_vector
      use interpolate, only: interp
      use kerr, only: kerr_metric, lnrf_frame, uks2ubl, bl2ks
      use phys_constants, only: pi
      use math, only: zbrent
      implicit none

      namelist /harm/  dfile, hfile, nt, indf, gfile

      character(len=100) :: dfile, hfile, gfile
      integer :: n, ndumps, nt, indf, dlen, nhead
      integer :: ndim=3
      integer, dimension(:), allocatable :: dumps
      real(kind=8) :: tstep=1d0, & !h, dx1, dx2, dx3, gam, startx1, startx2, startx3, &
           toffset=0d0,tcur,dx1,dx2, &
           dx3,a,gam,Rin,Rout,hslope,R0,ti,tj,tk,startx1, startx2,rdump_cnt,rdump01_cnt, &
           startx3, x10, x20, game, cour,dt,DTd,DTi,DTl,DTr,DTr01,dump_cnt,failed,xbr,&
           image_cnt,lim,game4, game5, tf, &
           fractheta,fracphi,rbr,npow2,cpow2,global_x10,global_x20,global_fracdisk, &
           global_fracjet,global_r0disk,global_rdiskend,global_r0jet,global_rjetend, &
           global_jetnu1, global_jetnu2, global_disknu1, global_disknu2,global_rsjet, &
           global_r0grid,BL,DOQDOTAVG, tavgstart,SMALL=1e-20
      integer :: SDUMP,eHEAT,eCOND,DONUCLEAR,DOFLR,DOCYLINDRIFYCOORDS,EVOLVEVPOT,N1,N2,N3,N1G, &
           N2G,N3G,nx1,nx2,nx3,nstep,starti,startj,startk,DOKTOT,NPR,DOPARTICLES,myNp,NPTOT
      real :: asim
      real, dimension(:), allocatable :: t
      real, dimension(:), allocatable :: x1_arr, x2_arr, x3_arr, r_arr, th_arr, ph_arr
      real, dimension(:), allocatable :: rho_arr, p_arr, u0_arr, vrl_arr, &
        vpl_arr, vtl_arr, temp_arr, kela_arr, kelb_arr, kelc_arr, keld_arr
! global variables for BL3 but should have these read in from header files rather than hard-coded
      real, dimension(:), allocatable :: b0_arr, br_arr, bth_arr, bph_arr, gdet
      real, dimension(:,:), allocatable :: drdx,metric_cov,metric_con
      
      interface read_bl3_vars
         module procedure read_bl3_vars
      end interface

      interface init_harmpi_data
        module procedure init_harmpi_data
      end interface

      interface del_harmpi_data
        module procedure del_harmpi_data
      end interface

      interface initialize_harmpi_model
        module procedure initialize_harmpi_model
      end interface

      interface transformbl2mksh
         module procedure transformbl2mksh
      end interface

      interface umksh2uks
         module procedure umksh2uks
      end interface

      interface harmpi_vals
        module procedure harmpi_vals
      end interface

      interface calcrmks
         module procedure calcrmks
         module procedure calcrmks_single
      end interface

      interface findx1mks
         module procedure findx1mks
         module procedure findx1mks_array
      end interface

      interface findx2harmpi
         module procedure findx2harmpi
         module procedure findx2harmpi_array
      end interface

      interface mins
         module procedure mins
         module procedure mins_single
      end interface

      interface limlin
         module procedure limlin
         module procedure limlin_array
      end interface

      contains

        function findx2harmpi(x2,args) result(diff)
        real(kind=8), intent(in) :: x2
        real(kind=8), intent(in), dimension(:) :: args
        real(kind=8) :: h, theta, diff
        h=args(2); theta=args(1)
        diff=theta-(pi/2.*(1.+x2)+((1.-h)/2.)*sin(pi*(1.+x2)))
        end function findx2harmpi

        function findx2harmpi_array(x2,args) result(diff)
        real(kind=8), intent(in), dimension(:) :: x2
        real(kind=8), intent(in), dimension(:,:) :: args
        real(kind=8) :: h
        real(kind=8), dimension(size(x2)) :: theta, diff
        h=args(1,2); theta=args(:,1)
        diff=theta-(pi/2d0*(1.+x2)+((1d0-h)/2d0)*sin(pi*(1d0+x2)))
        end function findx2harmpi_array

        subroutine transformbl2mksh(r,th,phi,x1,x2,x3)
        ! to 3D on 19/11/15 JD
        ! transform Boyer-Lindquist coordinates to modified Kerr-Schild coordinates used by HARM
        real, intent(in), dimension(:) :: r,th,phi
        real, intent(out), dimension(size(r)) :: x1,x2,x3
!        real(kind=8), intent(in) :: hslope
        real(kind=8) :: hval
        real(kind=8), dimension(2) :: args
        integer :: i
! temporary cheat like i used to do assuming all r < rbr
        x1=log(r)
        x3=phi
!        x2=th/pi
! CHECK THAT THIS WORKS WITH HSLOPE
!        hval=1d0
        hval=hslope
        args=(/th(1),real(hval)/)
        do i=1,size(r)
! put r part here           
           args(1)=th(i)
           x2(i)=zbrent(findx2harmpi,-1d0,1d0,dble(args),1d-6)
        enddo
        end subroutine transformbl2mksh

        function Ftr(x) result(res)
          real(kind=8), intent(in), dimension(:) :: x
          real(kind=8), dimension(size(x)) :: res
          real(kind=8), dimension(size(x)) :: zero,one
          zero=0.; one=1.
          res=merge(zero,(64+cos(5*pi*x)+70*sin((pi*(-1+2*x))/2.)+5*sin((3*pi*(-1+2*x))/2.))/128.,x.le.0.)
          res=merge(res,one,x.lt.1.)
        end function Ftr

!        function Ftrgenlin( x, xa, xb, ya, yb ) result(res)
!          real(kind=8), intent(in), dimension(:) :: x,
!          res = (x*ya)/xa + (-((x*ya)/xa) + ((x - xb)*(1 - yb))/(1 - xb) + yb)*Ftr((x - xa)/(-xa + xb));

        function Ftrgen( x, xa, xb, ya, yb ) result(res)
          real(kind=8), intent(in), dimension(:) :: x
          real(kind=8), intent(in) :: xa,xb,ya,yb
          real(kind=8), dimension(size(x)) :: res
          res = ya + (yb-ya)*Ftr( (x-xa)/(xb-xa) )
        end function Ftrgen

        function Fangle( x ) result(res)
          real(kind=8), intent(in), dimension(:) :: x
          real(kind=8), dimension(size(x)) :: res
          real(kind=8), dimension(size(x)) :: zero,val
          zero=0d0
          res=zero
          where(x.ge.-1.and.x.lt.1.)
             res=(1 + x + (-140*sin((pi*(1 + x))/2d0) + &
                  (10*sin((3*pi*(1 + x))/2d0))/3d0 + &
                  (2*sin((5*pi*(1 + x))/2d0))/5d0)/ &
                  (64d0*pi))/2d0
          endwhere
          where(x.gt.1d0)
             res=x
          endwhere
        end function Fangle

        function limlin( x, x0, dx, y0 ) result(res)
          real(kind=8), intent(in), dimension(:) :: x
          real(kind=8), intent(in) :: x0,dx,y0
          real(kind=8), dimension(size(x)) :: res
          res = y0 - dx * Fangle(-(x-x0)/dx)
        end function limlin

        function limlin_array( x, x0, dx, y0 ) result(res)
          real(kind=8), intent(in), dimension(:) :: x,x0,dx,y0
!          real(kind=8), intent(in) :: x0,dx,y0
          real(kind=8), dimension(size(x)) :: res
          res = y0 - dx * Fangle(-(x-x0)/dx)
        end function limlin_array
        
!        function limlin_single( x, x0, dx, y0 ) result(res)
!          real(kind=8), intent(in), dimension(:) :: x
!          real(kind=8), intent(in) :: x0,dx,y0
!          real(kind=8), dimension(size(x)) :: res
!          res = y0 - dx * Fangle(-(x-x0)/dx)
!        end function limlin_single

        function mins( f1, f2, df ) result(res)
          real(kind=8), intent(in), dimension(:) :: f1,f2,df
          real(kind=8), dimension(size(f1)) :: res
          res = limlin(f1, f2, df, f2)
        end function mins

        function mins_single( f1, f2, df ) result(res)
          real(kind=8), intent(in), dimension(:) :: f1
          real(kind=8), intent(in) :: f2,df
          real(kind=8), dimension(size(f1)) :: res
          res = limlin(f1, f2, df, f2)
        end function mins_single

        function maxs( f1, f2, df ) result(res)
          real(kind=8), intent(in), dimension(:) :: f1,f2,df
          real(kind=8), dimension(size(f1)) :: res
          res = -mins(-f1, -f2, df)
        end function maxs

        function minmaxs( f1, f2, df, dir ) result(res)
          real(kind=8), intent(in), dimension(:) :: f1,f2,df,dir
          real(kind=8), dimension(size(f1)) :: res,zero
          zero=0.
          res = zero
          where(dir.gt.0.)
             res = maxs(f1, f2, df)
          elsewhere
             res = mins(f1, f2, df)
          endwhere
        end function minmaxs

        subroutine to1stquadrant( x2in, x2out, ismirrored )
          real(kind=8), intent(in), dimension(:) :: x2in
          real(kind=8), intent(out), dimension(size(x2in)) :: x2out
          integer, intent(out), dimension(size(x2in)) :: ismirrored
          real(kind=8), dimension(size(x2in)) :: ntimes
      
          x2out = x2in
          ismirrored=0
          ntimes = floor( (x2in+2.0)/4.0 )
          x2out = x2out - 4.*ntimes

          where(x2out.gt.0.)
             x2out = -x2out
             ismirrored = 1-ismirrored
          endwhere
          where(x2out.lt.-1.)
             x2out = -2.-x2out
             ismirrored = 1-ismirrored
          endwhere
        end subroutine to1stquadrant

        function func1( r,x2 ) result(res)
          real(kind=8), intent(in), dimension(:) :: r,x2
          real(kind=8), dimension(size(x2)) :: res
!          external :: calcth
          res = sin(calcthmksbl3( x2,r ))
        end function func1

        function func2( r0, r, x20, x2) result(res)
          real(kind=8), intent(in), dimension(:) :: r,x2
          real(kind=8), intent(in) :: r0,x20
          real(kind=8), dimension(size(r)) :: res
          real(kind=8), dimension(size(r)) :: sth1in,sth2in,sth1inaxis,sth2inaxis,mone
!          external :: calcth
          real(kind=8) :: SMALL
          SMALL = 1e-20
          mone=-1.
          
!          th = calcth( r,x2 )

!          write(6,*) 'in func2: ',r0,x20
!          write(6,*) 'in func2: ',r,x2
          sth1in = sinth1in( r0, r, x20, x2)
          sth2in = sin( th2in(r0,r,x20,x2) )

          sth1inaxis = sinth1in( r0,r,x20,mone)
          sth2inaxis = sin( th2in(r0,r,x20,mone) )

          res = minmaxs( sth1in, sth2in, abs(sth2inaxis-sth1inaxis)+SMALL, r-r0 )
!          write(6,*) 'func2: ',res,sth1in,sth2in
!          write(6,*) 'func2: ',sth2inaxis,sth1inaxis
        end function func2

        function sinth0( x20, r0, r ) result(res)
          real(kind=8), intent(in), dimension(:) :: x20,r0,r
          real(kind=8),  dimension(size(r)) :: res
          real(kind=8), dimension(size(x20)) :: th0
!          external :: calcth
          th0 = calcthmksbl3( x20, r0 )
          res = r0(1) * sin(th0(1)) / r
        end function sinth0

        function th2in( r0,r,x20,x2 ) result(res)
          real(kind=8), intent(in), dimension(:) :: x2,r!,x20,r0
          real(kind=8), intent(in) :: x20,r0
!          external :: calcth
          real(kind=8), dimension(size(x2)) :: res
          real(kind=8), dimension(size(r)) :: zero,th0,thetamid, &
               theta,thetac
!          real(kind=8), dimension(1) :: theta0,thetac

          zero(:) = 0d0
          thetac = calcthmksbl3(x20+zero,r)
          thetamid = calcthmksbl3(zero,r)
!          theta0 = calcthmksbl3((/r0/), (/x20/))
          theta = calcthmksbl3(x2,r)

          th0 = asin( sinth0(x20+zero, r0+zero, r) )
          res = (theta - thetac)/(thetamid - thetac) * (thetamid-th0) + th0
!          write(6,*) 'in thetac: ',x20,r,r0,x2
!          write(6,*) 'in thetac: ',thetac,thetamid,theta
!          write(6,*) 'in thetac: ',th0,res
        end function th2in

        function sinth1in( r0, r, x20, x2) result(res)
          real(kind=8), intent(in), dimension(:) :: r,x2
          real(kind=8), intent(in) :: r0,x20
          real(kind=8), dimension(size(x2)) :: res,r0t
          real(kind=8), dimension(size(x2)) :: thc
!          external :: calcth
! do we also need x1 variables here? don't think so
          r0t(:)=r0
          thc = calcthmksbl3( x2,r0t )
          res =  r0 * sin(thc) / r
!          write(6,*) 'sinth1in: ',r0,thc,r
        end function sinth1in

        function thetaofx2(x2, ror0nu) result(res)
          real(kind=8), intent(in), dimension(:) ::  x2,ror0nu
!          real(kind=8), dimension(size(x2)) :: res
          real(kind=8), dimension(size(x2)) :: theta,theta2,theta3,res
          theta = atan( tan((x2+1)*pi/2.)/ror0nu )
          theta2 = (pi    + atan( tan((x2-1)*pi/2.)/ror0nu ))
          theta3 = (pi/2. + atan( tan(x2*pi/2.)*ror0nu ))
          res = merge(merge(theta,theta2,x2.lt.-0.5),theta3,x2.gt.0.5)
!          res = theta+theta2+theta3
        end function thetaofx2

        function calcth_cylindrified(x2in,rin) result(theta)
          real(kind=8), intent(in), dimension(:) :: rin,x2in
!          real(kind=8), dimension(size(rin)) :: theta
          real(kind=8), dimension(size(rin)) :: theta,x2mirror,thmirror, &
               thorig,f1,f2,sinth,th,dftr,rtr,x1tr,thtr
          integer, dimension(size(rin)) :: ismirrored
          real(kind=8) :: r0,x20
!          external :: calcth,calcr
!          SMALL=1e-20
          thorig = calcthmksbl3(x2in,rin)
          call to1stquadrant( x2in,x2mirror,ismirrored )
          thmirror = calcthmksbl3(x2mirror,rin)
          
          r0 = calcrmks(global_x10)
          x20 = global_x20

          x1tr = log( 0.5*( exp(global_x10)+exp(startx1) ) )   !always bound to be between startx[1] and r0
          rtr = calcrmks(x1tr)
          thtr = calcthmksbl3( x2mirror,rtr )
!          write(6,*) 'x1tr: ',x2mirror,x1tr,rtr,thtr
          
          f1 = func1( rin,x2mirror )
          f2 = func2( r0, rin, x20, x2mirror )
!          write(6,*) 'after f1f2: ',f1,f2
          dftr = func2( r0, rtr, x20, x2mirror ) - func1( rtr,x2mirror )
          sinth = maxs( rin*f1, rin*f2, rtr*abs(dftr)+SMALL ) / rin

!          write(6,*) 'after sinth: ',r0,rtr,x1tr
!          write(6,*) 'after sinth: ',sinth,dftr
          th = asin( sinth )
          theta = thorig
! presumably this can be replaced with a single line and sign() call...
          where(ismirrored.eq.0)
             theta = thorig + (th - thmirror)
          endwhere
          where(ismirrored.ne.0.)
             theta  = thorig  - (th  - thmirror )
          endwhere
        end function calcth_cylindrified

        subroutine read_bl3_vars()
          SMALL= 1e-20
          rbr = 400.
          R0 = 0.0
          Rin = 0.87*(1.+sqrt(1.-0.5**2.))
          xbr = log(rbr-R0)
          npow2 = 4.0
          cpow2 = 1.0
          hslope = 0.3
          startx1 = log(Rin-R0)
          
          fractheta = 1.
          fracphi = 1.
          
          global_fracdisk = 0.25 !fraction of resolution going into the disk
          global_fracjet = 0.40 !fraction of resolution going into the jets
          global_disknu1 = -2
          global_disknu2 = 0.75
          global_jetnu1 = -2  !the nu-parameter that determines jet shape at small radii (r < min(r0jet,r0disk))
          global_jetnu2 = 0.75  !the nu-parameter that determines jet shape at large radii (r > r0jet)
          global_rsjet = 0.0
          global_r0grid = Rin
          global_r0jet = 2*Rin
          global_rjetend = 1e3
          global_r0disk = 2*Rin
          global_rdiskend = 5*Rin
          global_x10 = 5.;  !radial distance in MCOORD until which the innermost angular cell is cylinrdical
          global_x20 = -1. + 1./256.
        end subroutine read_bl3_vars

        function calcthmksbl3(x2,r) result(theta)
          real(kind=8), intent(in), dimension(:) :: x2, r
          real(kind=8), dimension(size(x2)) :: fac,r1disk,r2disk,dr,r1jet,r2jet, &
               ror0nudisk,ror0nujet,thetadisk,thetajet
          real(kind=8), dimension(size(x2)) :: theta

          fac = Ftrgen( abs(x2), global_fracdisk, 1-global_fracjet, 0d0, 1d0 )
          r1disk = mins( r/global_r0disk, 1d0 , 0.5d0 ) * (global_r0disk/global_r0grid)
          r2disk = r/(r1disk*global_r0grid)
          dr = global_rdiskend/global_r0disk
          r2disk = mins( r2disk, dr, 0.5d0*dr )

          r1jet = mins( r/global_r0jet, 1d0 , 0.5d0 ) * (global_r0jet/global_r0grid)
          r2jet = r/(r1jet*global_r0grid)
          dr = global_rjetend/global_r0jet
          r2jet = mins( r2jet, dr, 0.5d0*dr )

          ror0nudisk = r1disk**(0.5d0*global_disknu1) * r2disk**(0.5d0*global_disknu2)
          ror0nujet = r1jet**(0.5d0*global_jetnu1) * r2jet**(0.5d0*global_jetnu2)

          thetadisk = thetaofx2( x2, ror0nudisk )
          thetajet = thetaofx2( x2, ror0nujet )
          theta = fac*thetajet + (1 - fac)*thetadisk

        end function calcthmksbl3

        function findx2mksbl3(x2,args) result(diff)
          ! Calculates \theta for harmpi BL=3 coordinates
          ! JAD 1/5/2018
          real(kind=8), intent(in), dimension(:) :: x2
          real(kind=8), intent(in), dimension(:,:) :: args
          real(kind=8), dimension(size(x2)) :: r,th,diff,theta
!          real(kind=8) :: r0r,r1jet,njet,r0jet,rsjet,qjet, &
!               rs,r0,r0jet3,rsjet3,h0,ntheta,htheta,rsjet2,r0jet2,myhslope,th2, &
!               th0,switch0,switch2,theta1,theta2,arctan2
          th=args(:,2); r=args(:,1)
!          theta=calcth_cylindrified(x2,r)
          theta=calcthmksbl3(x2,r)
          diff=th-theta
!          write(6,*) 'findx2mksbl3: ',th(1),theta(1),x2(1),r(1)
!          write(6,*) 'findx2 diff: ',diff
        end function findx2mksbl3

        function calcrmks(x1) result(r)
          ! Compute r given x1 for Jon's simulations
          ! JAD 7/24/2011
          real(kind=8), intent(in), dimension(:) :: x1
!          real(kind=8), intent(in) :: xbr,npow2,cpow2,R0
!          real(kind=8) :: npow2,cpow2,R0
          real(kind=8), dimension(size(x1)) :: r, xi
!          npow2=10d0
          where(x1.gt.xbr)
             xi=x1+cpow2*(x1-xbr)**npow2
          elsewhere
             xi=x1
          endwhere
!          write(6,*) 'calcrmks: ',R0,cpow2,npow2,xbr
          r=R0+exp(xi)
        end function calcrmks

        function calcrmks_single(x1) result(r)
          ! Compute r given x1 for Jon's simulations
          ! JAD 7/24/2011
          real(kind=8), intent(in) :: x1
!          real(kind=8), intent(in) :: xbr,npow2,cpow2,R0
!          real(kind=8) :: npow2,cpow2,R0
          real(kind=8)  :: r, xi
!          npow2=10d0
          if(x1.gt.xbr) then
             xi=x1+cpow2*(x1-xbr)**npow2
          else
             xi=x1
          endif
          r=R0+exp(xi)
        end function calcrmks_single

        function findx1mks(x1,args) result(diff)
          real(kind=8), intent(in), dimension(:) :: args
          real(kind=8), intent(in) :: x1
          real(kind=8) :: diff,r
!          real(kind=8) :: xbr,npow2,r,xi,diff,R0,cpow2
        ! Compute x1 given r for Jon's simulations 
          R0=args(2)
          xbr=args(3)
          npow2=args(4)
          cpow2=args(5)
          r=args(1)
          diff=r-calcrmks(x1)
        end function findx1mks

        function findx1mks_array(x1,args) result(diff)
          real(kind=8), intent(in), dimension(:,:) :: args
          real(kind=8), intent(in), dimension(:) :: x1
          real(kind=8), dimension(size(x1)) :: diff,r
!          real(kind=8) :: xbr,npow2,r,xi,diff,R0,cpow2
        ! Compute x1 given r for Jon's simulations 
!          R0=args(2)
!          xbr=args(3)
!          npow2=args(4)
!          cpow2=args(5)
          r=args(:,1)
          diff=r-calcrmks(x1)
        end function findx1mks_array

        subroutine transformbl2mksbl3(r,th,ph,x1,x2,x3)
        ! transform Boyer-Lindquist coordinates to modified Kerr-Schild coordinates used by thickdisk
        real, intent(in), dimension(:) :: r,th,ph
        real, intent(out), dimension(size(r)) :: x1,x2,x3
        real(kind=8), dimension(5) :: args
        real(kind=8), dimension(size(r),2) :: thargs
        real(kind=8), dimension(size(r)) :: mone,one
        integer :: i
        mone(:)=-1d0; one(:)=1d0
!        real, intent(in) :: xbr,cpow2,npow2,R0
! this should be init somewhere once and then left alone. minimal cost but no need to re-assign global vars over and over. should really just be part of reading the header so move it there.
!        call read_bl3_vars()
        args(2)=R0; args(3)=xbr; args(4)=npow2; args(5)=cpow2
        do i=1,size(r)
           args(1)=r(i)
           x1(i)=zbrent(findx1mks,0d0,10d0,args,1d-6)
!           args(2)=th(i)
!           x2(i)=zbrent(findx2mksbl3,0d0,1d0,dble(args),1d-6)
        enddo
! test call to zbrent_array so that i don't have to write *_single versions of every sasha BL3 routine
!        if(BL==3) then
        thargs(:,1)=dble(r)
        thargs(:,2)=dble(th)
        x2=zbrent(findx2mksbl3,mone,one,thargs,1d-6)
!        else
!           thargs(:,1)=dble(hval)
!           thargs(:,2)=dble(th)
!           x2=zbrent(findx2harmpi,mone,one,thargs,1d-6)
!        end if
        x3=ph
        end subroutine transformbl2mksbl3

        function umks2uksbl3(fum,x1,x2) result(uks)
          real(kind=8), intent(in), dimension(:) :: x1,x2
          type (four_vector), dimension(:), intent(in) :: fum
!          real(kind=8), intent(in) :: xbr
          type (four_vector), dimension(size(fum)) :: uks
          real(kind=8), dimension(size(x1)) :: r,dr,dx2,dx1,drdx1,dthdx1,dthdx2
          ! Convert four-vectors from MKS to KS numerically using central differences.
          ! JAD 5/24/2010 for harmpi Sasha BL=3 grid 12/30/2017
          ! alternatively just use gdump file
!  CHANGED to use dth/dx1 instead of dth/dr since that's what i have from gdump file
! i think the old way with dthdr = dthdx1 / drdx1 was fine b/c r(x1) is 
!  independent of other coordinates but not sure
          if(allocated(drdx)) then
             drdx1=drdx(6,:)
             dthdx1=drdx(10,:)
             dthdx2=drdx(11,:)
          else
          ! parameters are from comparing numerical and analytical values 
          ! using the jetcoords3 grid of Jon's but now using for harmpi.
             dr=1d-4*r; dx2=1d-6*x2; dx1=1d-4*x1
             drdx1=(calcrmks(x1+.5*dx1)-calcrmks(x1-.5d0*dx1))/dx1
!             dthdr=(calcthmksbl3(x2,r+.5d0*dr)-calcthmksbl3(x2,r-.5d0*dr))/dr
             dthdx1=(calcth_cylindrified(x2,calcrmks(x1+.5*dx1))- &
                  calcth_cylindrified(x2,calcrmks(x1-.5*dx1)))/dx1
             dthdx2=(calcth_cylindrified(x2+.5d0*dx2,r)-calcth_cylindrified(x2-.5d0*dx2,r))/dx2
!             dthdr=(calcthmksbl3(x2,r+.5d0*dr)-calcthmksbl3(x2,r-.5d0*dr))/dr
!             dthdx1=(calcthmksbl3(x2,calcrmks(x1+.5*dx1))- &
!                  calcthmksbl3(x2,calcrmks(x1-.5*dx1)))/dx1
!             dthdx2=(calcthmksbl3(x2+.5d0*dx2,r)-calcthmksbl3(x2-.5d0*dx2,r))/dx2
          end if
          r=calcrmks(x1)!,xbr)
          uks%data(2)=drdx1*fum%data(2)
          uks%data(3)=dthdx1*fum%data(2)+dthdx2*fum%data(3)
!          uks%data(3)=fum%data(3)*dthdx2+uks%data(2)*dthdr
          ! phi doesn't change
          uks%data(4)=fum%data(4)
          ! time component doesn't change
          uks%data(1)=fum%data(1)
        end function umks2uksbl3

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
          dthdx2=pi/2.*(1.+(1.-hslope)*cos(pi*(1.+x2)))
          uth=fum%data(3)*dble(dthdx2)
          uph=fum%data(4)
          uks=fum
          uks%data(2)=ur
          uks%data(3)=uth
          uks%data(4)=uph
        end function umksh2uks

        subroutine harmpi_vals(x0,a,rho,p,b,u,bmag,kela,kelb, &
             kelc,keld)
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
         vpli,bri,bthi,bphi,b0i,u0i,ri,thi,phii,kelai,kelbi,kelci,keldi
!        integer, dimension(size(x0)*(2**ndim)) :: sindx
        real, dimension(size(x0)) :: dth, minph
        integer, dimension(size(x0)*2*(2**ndim)) :: indx
        integer, dimension(size(x0)) :: lx1,lx2, lx3, ux1,ux2,ux3,x1l,x1u,x2l,x2u,x3l,x3u, &
         umax,one,tindx
        integer :: npts,i,maxtindx
        real, dimension(size(x0)), intent(out) :: rho,p,bmag,kela,kelb, &
             kelc,keld
        type (four_Vector), intent(out), dimension(size(x0)) :: u,b
        ! Interpolates HARM data to input coordinates
        ! JAD 3/20/2009, fortran 11/12/2012
!        write(6,*) 'harm: ',size(x0)
        !nx1=floor(sqrt(real(size(x1_arr)))); nx2=nx1
!        write(6,*) 'harmpi vals: ',nx1,nx2,nx3,size(x1_arr),size(x2_arr),size(x3_arr)
!        write(6,*) 'harmpi vals'
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
! JAD 10/28/2019 ADDING ANOTHER ROUND FOR < -2\pi!
        where(zphi.lt.0.)
           zphi=zphi+2.*pi
        endwhere
!        write(6,*) 'harmpi vals transform'
        if(BL.eq.3) then
           call transformbl2mksbl3(zr,theta,zphi,x1,x2,x3)
        else
           call transformbl2mksh(zr,theta,zphi,x1,x2,x3)
        end if

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
!        write(6,*) 'nx1 nx2 nx3: ',nx1,nx2,nx3,minval(zphi),maxval(zphi),minval(zr),maxval(zr),minval(theta), &
!             maxval(theta)
!       write(6,*) 'ux lx: ',minval(lx1),maxval(ux1),minval(lx2),maxval(ux2), minval(lx3), maxval(ux3)

! commenting out below while interpolation is off
! Deal with poles
!        write(6,*) 'poles'
!        if(BL.eq.1) then
!           where(ux2.ne.lx2)
!              dth=uniqth(ux2)-uniqth(lx2)
!           elsewhere
!              dth=uniqth(1)
!           endwhere
!        else
!!           where(ux2.ne.lx2)
!!              dth=calcth_cylindrified(dble(uniqx2(ux2)),dble(zr))-calcth_cylindrified(dble(uniqx2(lx2)),dble(zr))
!!           elsewhere
!!              dth=calcth_cylindrified(dble(uniqx2(1))+0d0*dble(zr),dble(zr))
!!           endwhere
!
!           where(ux2.ne.lx2)
!              dth=calcthmksbl3(dble(uniqx2(ux2)),dble(zr))-calcthmksbl3(dble(uniqx2(lx2)),dble(zr))
!           elsewhere
!              dth=calcthmksbl3(dble(uniqx2(1))+0d0*dble(zr),dble(zr))
!           endwhere
!        end if
! periodic in phi
!        minph=uniqph(lx3)
        where(ux3.gt.nx3)
            ux3=1
        endwhere
        where(lx3.lt.1)
            lx3=nx3
            minph=uniqph(lx3)-2.*pi
        endwhere
!        write(6,*) 'uniform phi'
! uniform in phi
!        pd=(zphi-minph)/(uniqph(2)-uniqph(1))
!        if(BL.eq.3) then
!!           td=abs(theta-calcth_cylindrified(dble(uniqx2(lx2)),dble(zr)))/dth
!           td=abs(theta-calcthmksbl3(dble(uniqx2(lx2)),dble(zr)))/dth
!        else
!           td=abs(theta-uniqth(lx2))/dth
 !       end if
 !       rd=(zr-uniqr(lx1))/(uniqr(ux1)-uniqr(lx1))
! When the geodesic is inside the innermost zone center located outside the horizon, 
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
!        write(6,*) 'ux lx: ',minval(lx1),maxval(ux1),minval(lx2),maxval(ux2),minval(lx3),maxval(ux3)
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
!        write(6,*) 'before keli', eCOND,eHEAT,allocated(kel_arr)
        if(eCOND.eq.1.or.eHEAT.eq.1) then
           kelai=reshape(kela_arr(indx),(/npts,2**(ndim+1)/))
           kelbi=reshape(kelb_arr(indx),(/npts,2**(ndim+1)/))
           kelci=reshape(kelc_arr(indx),(/npts,2**(ndim+1)/))
           keldi=reshape(keld_arr(indx),(/npts,2**(ndim+1)/))
        endif
! coordinates for debugging
        ri=reshape(r_arr(indx),(/npts,2**(ndim+1)/))
        thi=reshape(th_arr(indx),(/npts,2**(ndim+1)/))
        phii=reshape(ph_arr(indx),(/npts,2**(ndim+1)/))
        rttd=0.
        rho=merge(interp(rhoi,rttd,pd,rd,td),dzero,x1.gt.uniqx1(1))*nfac
        p=merge(interp(ppi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))*pfac
        if(eCOND.eq.1.or.eHEAT.eq.1) then
           kela=merge(interp(kelai,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))*pfac
           kelb=merge(interp(kelbi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))*pfac
           kelc=merge(interp(kelci,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))*pfac
           keld=merge(interp(keldi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))*pfac
        else
           kela=0d0; kelb=0d0; kelc=0d0; keld=0d0
        end if
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
!        write(6,*) 'min vals u', minval(u%data(1)), maxval(abs(u*u+1))
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
        end subroutine harmpi_vals

        subroutine read_harmpi_inputs(ifile)
        character(len=20), intent(in) :: ifile
        open(unit=8,file=ifile,form='formatted',status='old')
        read(8,nml=harm)
        write(6,*) 'read: ',nt
        close(unit=8)
        end subroutine read_harmpi_inputs

!        subroutine read_harmpi_data_header(nhead)
!        integer, intent(in) :: nhead
!        real(8), dimension(nhead) :: header
!        real(8) :: tcur
        ! Read HARM header file
        ! JAD 11/24/2012 based on previous IDL codes
        ! header format from HARM dump.c
!        write(6,*) 'read harm header: ',nhead
!        open(unit=8,file=hfile,form='formatted',status='old')
!        read(8,*) header
!        nx1=header(2)
!        nx2=header(3)
!        nx3=header(4)
!        startx1=header(11)
!        startx2=header(12)
!        startx3=header(13)
!        dx1=header(14)
!        dx2=header(15)
!        dx3=header(16)
!        asim=header(19)
!        gam=header(20)
!        h=header(36)
! CHANGE READ MORE BASED ON VALUE OF BL
!        tcur=header(1)
!        write(6,*) 'harmpi header: ',tcur,nx1,nx2,nx3,startx1,startx2
!        write(6,*) 'harmpi header: ',dx3,asim,gam,h,header(34:38)
! is this true for HARM?
!        defcoord=header(19)
!        close(unit=8)
!        end subroutine read_harmpi_data_header

        subroutine read_harmpi_data_header()
!          integer, intent(in) :: nhead
          integer :: n,nheadmax=200,status
          real(kind=8), dimension(:), allocatable :: header
          nhead=70
          allocate(header(nheadmax)); header(:)=-1d20
          open(unit=8,file=hfile)
          read(8,*,iostat=status) header
          close(unit=8)
!          write(6,*) 'header: ',header
          n = 1
!          nhead=70
!          write(6,*) 'minloc: ',minloc(header),size(header),size(minloc(header))
          nhead=minloc(header,1)-1
!          nhead=nhead-1
          write(6,*) 'header: ',nhead,header(nhead),header(nhead+1)
          tcur = dble(header(n)); n=n+1
!per tile resolution
          N1 = int(header(n)); n=n+1
          N2 = int(header(n)); n=n+1
          N3 = int(header(n)); n=n+1
        !total resolution
          nx1 = int(header(n)); n=n+1
          nx2 = int(header(n)); n=n+1
          nx3 = int(header(n)); n=n+1
        !numbers of ghost cells
          N1G = int(header(n)); n=n+1
          N2G = int(header(n)); n=n+1
          N3G = int(header(n)); n=n+1
          startx1 = dble(header(n)); n=n+1
          startx2 = dble(header(n)); n=n+1
          startx3 = dble(header(n)); n=n+1
          dx1=dble(header(n)); n=n+1
          dx2=dble(header(n)); n=n+1
          dx3=dble(header(n)); n=n+1
          tf=dble(header(n)); n=n+1
          nstep=dble(header(n)); n=n+1
          asim=dble(header(n)); n=n+1
          gam=dble(header(n)); n=n+1
          cour=dble(header(n)); n=n+1
          DTd=dble(header(n)); n=n+1
          DTl=dble(header(n)); n=n+1
          DTi=dble(header(n)); n=n+1
          DTr=dble(header(n)); n=n+1
          DTr01=dble(header(n)); n=n+1
          dump_cnt=dble(header(n)); n=n+1
          image_cnt=dble(header(n)); n=n+1
          rdump_cnt=dble(header(n)); n=n+1
          rdump01_cnt=dble(header(n)); n=n+1
          dt=dble(header(n)); n=n+1
          lim=dble(header(n)); n=n+1
          failed=dble(header(n)); n=n+1
          Rin=dble(header(n)); n=n+1
          Rout=dble(header(n)); n=n+1
          hslope=dble(header(n)); n=n+1
          R0=dble(header(n)); n=n+1
! up until here we are compatible with both public / private harmpi but then we split
! now we separate for public version (nhead=46)
          if (nhead==46) then
             NPR=int(header(n)); n=n+1
             DOKTOT=int(header(n)); n=n+1
             fractheta = dble(header(n)); n=n+1
             fracphi   = dble(header(n)); n=n+1
             rbr       = dble(header(n)); n=n+1
             xbr = log(rbr-R0)
             npow2     = dble(header(n)); n=n+1
             cpow2     = dble(header(n)); n=n+1
             BL = dble(header(n)); n=n+1
             eHEAT=-1
             eCOND=-1
             DONUCLEAR=0
             DOFLR = 0
             DOCYLINDRIFYCOORDS = 0
          else
             if (n<nhead) then
                NPR=int(header(n)); n=n+1
                DOKTOT=int(header(n)); n=n+1
                eHEAT=int(header(n)); n=n+1
                eCOND=int(header(n)); n=n+1
                DONUCLEAR=int(header(n)); n=n+1
!             nheader = 42
             else
                NPR=-1
                DOKTOT=-1
                eHEAT=-1
                eCOND=-1
                DONUCLEAR=0
                DOFLR = 0
             end if
             if (n<nhead) then
                DOFLR=int(header(n)); n=n+1
                !             nheader = 43
             else
                DOFLR = 0
             endif
             if (n<nhead) then
                !             nheader = 60
                DOCYLINDRIFYCOORDS = dble(header(n)); n=n+1
                fractheta = dble(header(n)); n=n+1
                fracphi   = dble(header(n)); n=n+1
                rbr       = dble(header(n)); n=n+1
                xbr = log(rbr-R0)
                npow2     = dble(header(n)); n=n+1
                cpow2     = dble(header(n)); n=n+1
                global_x10 = dble(header(n)); n=n+1
                global_x20 = dble(header(n)); n=n+1
                global_fracdisk   = dble(header(n)); n=n+1
                global_fracjet    = dble(header(n)); n=n+1
                global_r0disk     = dble(header(n)); n=n+1
                global_rdiskend   = dble(header(n)); n=n+1
                global_r0jet      = dble(header(n)); n=n+1
                global_rjetend    = dble(header(n)); n=n+1
                global_jetnu2      = dble(header(n)); n=n+1
                global_rsjet      = dble(header(n)); n=n+1
                global_r0grid     = dble(header(n)); n=n+1
             endif
             if (n<nhead) then
!             nheader = 61
                BL = dble(header(n)); n=n+1
             endif
             if (n<nhead) then
!             nheader = 62
                EVOLVEVPOT = int(header(n)); n=n+1
             else
                EVOLVEVPOT = 0
             endif
             if (n<nhead) then
                !             nheader = 63
                global_jetnu1 = dble(header(n)); n=n+1
             else
                global_jetnu1 = 0
             end if
             if (n<nhead) then
                !             nheader = 64
                global_disknu1 = dble(header(n)); n=n+1
             else
                global_disknu1 = 0
             end if
             if (n<nhead) then
                !             nheader = 65
                global_disknu2 = dble(header(n)); n=n+1
             else
                global_disknu2 = 0
             end if
             if (n<nhead) then
                myNp = int(header(n)); n=n+1
                NPTOT = int(header(n)); n=n+1
                DOPARTICLES = 1
                !            nheader = 67
             end if
             if (n<nhead) then
                SDUMP = int(header(n)); n=n+1
                !            nheader = 68
             else
                SDUMP = 0
             end if
             if (n<nhead) then
                DOQDOTAVG = int(header(n)); n=n+1
                tavgstart = dble(header(n)); n=n+1
                !            nheader = 70
             else
                DOQDOTAVG = 0
                tavgstart = -1d0
             end if
          endif
          deallocate(header)
          if(SDUMP.eq.1) then
             if(eCOND.eq.1.or.eHEAT.eq.1) then
                dlen=NPR+3
             else
                dlen=NPR
             endif
          else
! these are untested so be careful
             if(eCOND.eq.1.or.eHEAT.eq.1) then
                dlen=58-19+NPR
             else
                dlen=42
             end if
          endif
          write(6,*) 'header: ',eHEAT,SDUMP,dlen
        end subroutine read_harmpi_data_header

! read gdump file after Sasha's harm_script.py
        subroutine read_harmpi_grid_file(grid_file)
          character(len=100), intent(in) :: grid_file
          character :: header_byte
          integer :: nheader_bytes,n,grid_dlen
          real(8), dimension(:,:), allocatable :: grid
          real(4), dimension(:,:), allocatable :: data
          grid_dlen=58
          allocate(data(grid_dlen,nx1*nx2*nx3))
          open(unit=8,file=grid_file,form='unformatted',access='stream',status='old',action='read')
          nheader_bytes=1
          read(8) header_byte
          do while(header_byte.ne.char(10)) 
             read(8) header_byte
             nheader_bytes=nheader_bytes+1
          end do
          read(8) data
          close(unit=8)
! first read grid
          allocate(grid(nx1*nx2*nx3,6))
          grid=transpose(data(4:9,:))
          x1_arr=grid(:,1); x2_arr=grid(:,2); x3_arr=grid(:,3)
          r_arr=grid(:,4); th_arr=grid(:,5); ph_arr=grid(:,6)
          deallocate(grid)
          write(6,*) 'harmpi grid: ',maxval(x1_arr),maxval(x2_arr),minval(th_arr)
! now metric
! want only unique components to use for your four vector routines which are
! 1 2 3 4 (5) 6 7 8 (9 10) 11 12 (13 14 15) 16
! 10:25, 26:41
          metric_cov=data((/10,11,12,13,15,16,17,20,21,25/),:)
          metric_con=data((/26,27,28,29,31,32,33,36,37,41/),:)
          gdet=data(42,:)
          drdx=data(43:58,:)
          write(6,*) 'harmpi grid gdet: ',minval(gdet),maxval(gdet)
          deallocate(data)
        end subroutine read_harmpi_grid_file

        subroutine read_harmpi_data_file(data_file,rho,p,u,b,kela, &
             kelb,kelc,keld,mdot)
          character(len=100), intent(in) :: data_file
!          character(len=100) :: data_file_app
          character :: header_byte
          integer :: nhead,nheader_bytes
          integer :: rhopos,ppos,vpos,bpos,gdetpos,kelapos, &
               kelbpos,kelcpos,keldpos
!          logical, intent(in), optional :: gridonly
          integer, intent(in), optional :: mdot
!          logical :: gridread
!          real(8), intent(out) :: tcur
          real(8), dimension(:), allocatable, intent(out) :: p,rho,kela, &
               kelb,kelc,keld
          real(8), dimension(:), allocatable :: udotu,bdotu,bdotb,pmin
          real(8), dimension(:), allocatable :: alpha,gamma,zero
          type (four_vector), dimension(:), allocatable, intent(out) :: u,b
          type (four_vector), dimension(:), allocatable :: uks,bks,beta
          real(8), dimension(:,:), allocatable :: grid
          real(4), dimension(:,:), allocatable :: data
          integer :: i,j, nelem, nelem2
          ! Read HARM data from a file of given line length, 
          ! positions of desired variables
          ! FIX these data lengths must depend on eHEAT and other arguments
! SET DLEN IN READ HEADER AS A GLOBAL VAR FOR DUMP, SDUMP, GDUMP
!          allocate(header(nhead))
          write(6,*) 'in read harmpi data file: ',SDUMP,data_file,dlen, &
               nx1*nx2*nx3
          if(SDUMP.eq.0) then
             rhopos=9; ppos=10; kelapos=rhopos+8; kelbpos=kelapos+1
             kelcpos=kelapos+2; keldpos=kelapos+3
             vpos=18+DOKTOT; bpos=vpos+8
          else
             rhopos=0; ppos=1; vpos=2; bpos=5; kelapos=rhopos+10
             kelbpos=kelapos+1; kelcpos=kelapos+2; keldpos=kelapos+3
! this should be kel4a but can check in your new ones with no conduction
          end if
          write(6,*) 'rhopos: ',rhopos,ppos,vpos,bpos,kelapos
          write(6,*) 'data file: ',data_file
          write(6,*) 'harmpi data size: ',dlen,nx1,nx2,nx3
!          write(6,*) 'harmpi data file: ',data_file
          allocate(data(dlen,nx1*nx2*nx3))
          open(unit=8,file=data_file,form='unformatted',access='stream',status='old',action='read')
          nheader_bytes=1
          read(8) header_byte
          do while(header_byte.ne.char(10)) 
             read(8) header_byte
             nheader_bytes=nheader_bytes+1
          end do
          write(6,*) 'read header: ',nheader_bytes,dlen,SDUMP
          read(8) data
          close(8)
          write(6,*) 'harmpi read done with data'
          if(SDUMP.eq.0) then
             allocate(grid(nx1*nx2*nx3,6))
             grid=transpose(data(4:9,:))
             x1_arr=grid(:,1); x2_arr=grid(:,2); x3_arr=grid(:,3)
             r_arr=grid(:,4); th_arr=grid(:,5); ph_arr=grid(:,6)
             deallocate(grid)
             write(6,*) 'read harm grid sizes', size(x1_arr), size(x2_arr), &
                  size(x3_arr), size(r_arr), size(th_arr), size(ph_arr)
! case where we are using a regular dump file but only to get the grid
! replacing with gdump for metric since sdump files use three-vector v^i, B^i
!             if(gridread) return
          endif
          allocate(p(nx1*nx2*nx3))
          allocate(rho(nx1*nx2*nx3)); allocate(pmin(nx1*nx2*nx3))
          allocate(u(nx1*nx2*nx3)); allocate(b(nx1*nx2*nx3))
!          allocate(gdet(nx1*nx2*nx3))
          allocate(uks(nx1*nx2*nx3)); allocate(bks(nx1*nx2*nx3))
          rho=data(rhopos+1,:)
          p=data(ppos+1,:)
          if(eCOND.eq.1.or.eHEAT.eq.1) then
             allocate(kela(nx1*nx2*nx3))
             allocate(kelb(nx1*nx2*nx3))
             allocate(kelc(nx1*nx2*nx3))
             allocate(keld(nx1*nx2*nx3))
             kela=data(kelapos+1,:)
             kelb=data(kelbpos+1,:)
             kelc=data(kelcpos+1,:)
             keld=data(keldpos+1,:)
          end if
          if(SDUMP.eq.0) then
             b%data(1)=data(bpos+1,:)
             b%data(2)=data(bpos+2,:)
             b%data(3)=data(bpos+3,:)
             b%data(4)=data(bpos+4,:)
!                write(6,*) 'sizes b: ',size(b%data(4)),size(data(bpos+4,:))
             u%data(1)=data(vpos+1,:)
             u%data(2)=data(vpos+2,:)
             u%data(3)=data(vpos+3,:)
             u%data(4)=data(vpos+4,:)
          else
! for small dumps these are 3-vector velocities that need to be converted using metric from gdump
             allocate(alpha(nx1*nx2*nx3))
             allocate(beta(nx1*nx2*nx3))
             allocate(gamma(nx1*nx2*nx3))
             allocate(zero(nx1*nx2*nx3))
             zero(:)=0d0
             b%data(1)=0d0
             b%data(2)=data(bpos+1,:)
             b%data(3)=data(bpos+2,:)
             b%data(4)=data(bpos+3,:)
             u%data(1)=0d0
             u%data(2)=data(vpos+1,:)
             u%data(3)=data(vpos+2,:)
             u%data(4)=data(vpos+3,:)
             alpha = 1d0/sqrt(-metric_con(1,:))
             beta%data(1)=1d0
             beta%data(2)=metric_con(2,:)*alpha*alpha
             beta%data(3)=metric_con(3,:)*alpha*alpha
             beta%data(4)=metric_con(4,:)*alpha*alpha
             call assign_metric(u,metric_cov)
             call assign_metric(b,metric_cov)
             gamma=u*u
             gamma=sqrt(1d0+merge(gamma,zero,gamma.gt.0d0))
!             gamma = (metric_cov(6,:)*u%data(2)**2.+metric_cov(11,:)*u%data(3)**2.+&
!                  metric_cov(16,:)*u%data(4)**2.+2.*(metric_cov(7,:)*u%data(2)*u%data(3)+&
!                  metric_cov(8,:)*u%data(2)*u%data(4)+metric_cov(9,:)*u%data(3)*u%data(4)))
             u%data(1) = gamma/alpha
             u%data(2) = u%data(2)-gamma*beta%data(2)/alpha
             u%data(3) = u%data(3)-gamma*beta%data(3)/alpha
             u%data(4) = u%data(4)-gamma*beta%data(4)/alpha
             deallocate(alpha); deallocate(beta); deallocate(gamma)
             deallocate(zero)
             b%data(1) = u*b
             b%data(2) = (b%data(2)+b%data(1)*u%data(2))/u%data(1)
             b%data(3) = (b%data(3)+b%data(1)*u%data(3))/u%data(1)
             b%data(4) = (b%data(4)+b%data(1)*u%data(4))/u%data(1)
          endif
          deallocate(data)
! test code for these transformations
!          allocate(bdotu(n)); allocate(udotu(n)); allocate(bdotb(n))
!          udotu=abs(u*u+1.)
!          write(6,*) 'udotu: ',maxval(udotu),udotu(30*nx1*nx2+50)
!          bdotb = b*b; bdotu = b*u
!          write(6,*) 'bdotu: ',maxval(abs(bdotu)),bdotu(30*nx1*nx2+50), &
!               bdotb(30*nx1*nx2+50)
!          deallocate(bdotu); deallocate(udotu); deallocate(bdotb)
          p=merge(p,pmin,p.gt.pmin)
          write(6,*) "min vals rho", minval(rho)
          write(6,*) "min vals p", minval(p)
          if(allocated(kela)) write(6,*) 'min vals kel: ',minval(kela)
          write(6,*) "min vals u", minval(u%data(1)), minval(u%data(2)), minval(u%data(3)), minval(u%data(4))
          write(6,*) "max vals u", maxval(u%data(1)), maxval(u%data(2)), maxval(u%data(3)), maxval(u%data(4))
!          if (present(mdot)) then
             ! Calculate accretion rate in code units:
             !  nx1=n_elements(uniqx1) ; nx2=n_elements(uniqx2) ; nx3=n_elements(uniqx3) 
!             dx2=uniqx2(2)-uniqx2(1) ; dx3=uniqx3(2)-uniqx3(1)
!             mdotarr=-1.*sum(sum(reform(gdet*rho*v(:,1),nx1,nx2,nx3),3),2)*dx2*dx3
!          endif
          ! Transform velocities, magnetic fields from MKS to KS and then BL:
          write(6,*) 'read harm transform coords r ', minval(r_arr), maxval(r_arr), asim
          write(6,*) 'read harm transform coords th ',minval(x2_arr), maxval(x2_arr)
          write(6,*) 'read harm transform coords phi',minval(x3_arr), maxval(x3_arr)
          write(6,*) 'read harm transform coords u ',minval(u%data(1)),minval(u%data(2))
          write(6,*) 'read harm transform coords u size ',size(u),hslope
          if(BL.eq.3) then
             uks = umks2uksbl3(u,dble(x1_arr),dble(x2_arr))
             bks = umks2uksbl3(b,dble(x1_arr),dble(x2_arr))
          else
             uks = umksh2uks(u,x1_arr,x2_arr)
             bks = umksh2uks(b,x1_arr,x2_arr)
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
          write(6,*) 'after transform',p(5*nx1+23), rho(7*nx1+131), &
               maxval(abs(udotu)), maxval(abs(bdotu)),minval(bdotb)
          write(6,*) 'udotu: ',udotu(30*nx1*nx2+50:30*nx1*nx2+60)
          write(6,*) 'bdotb: ',bdotb(30*nx1*nx2+50:30*nx1*nx2+60)
          write(6,*) 'bdotu: ',bdotu(30*nx1*nx2+50:30*nx1*nx2+60)
          write(6,*) 'x1: ',x1_arr(30*nx1*nx2+50),x2_arr(30*nx1*nx2+60)
          deallocate(udotu); deallocate(bdotu); deallocate(bdotb)
          deallocate(uks); deallocate(bks); deallocate(pmin)
        end subroutine read_harmpi_data_file

        subroutine initialize_harmpi_model(a,ifile,df,hf,gf,ntt,indft)
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
              call read_harmpi_inputs(ifile)
           else
              call read_harmpi_inputs(default_ifile)
           endif
        endif
        call read_harmpi_data_header()
        if (abs(asim-a).gt.1e-4) then 
           write(6,*) 'ERROR -- Different simulation and grtrans spin values!', asim, a
           return
        endif
        write(6,*) 'read header', SDUMP, BL, eCOND, asim
        n=nx1*nx2*nx3
        call init_harmpi_data(n,n*nt)
        write(6,*) 'init data', nt, indf, hfile, dfile
        call load_harmpi_data(nt)
        end subroutine initialize_harmpi_model

        subroutine load_harmpi_data(nt)
        real(8), dimension(:), allocatable :: rho,p,kela, &
             kelb,kelc,keld
        real, dimension(:), allocatable :: vrl, vtl, vpl
        type (four_vector), dimension(:), allocatable :: u, b
        integer, intent(in) :: nt
        character(len=20) :: append
        character(len=100) :: data_file
        integer :: k
        real(8) :: tcur, tstep_test
        allocate(rho(n)); allocate(p(n))
        allocate(vrl(n)); allocate(vtl(n)); allocate(vpl(n))
        allocate(u(n)); allocate(b(n))
        if(eCOND.eq.1.or.eHEAT.eq.1) then
           allocate(kela(n)); allocate(kelb(n))
           allocate(kelc(n)); allocate(keld(n))
        endif
! if we're using small dumps then get grid and metric from a regular dump file
        write(6,*) 'eHEAT: ',eHEAT,allocated(kela)
        if(SDUMP.eq.1) then
           call init_harmpi_grid_data(nx1*nx2*nx3)
           call read_harmpi_grid_file(gfile)
        endif
        write(6,*) 'grid: ',SDUMP,nt
! now loop over and load data files
        do k=1,nt
           write(append, fmt='(I5.3)') indf-(k-1)
           data_file = trim(dfile) // trim(adjustl(append))
           write(6,*) 'data_file: ',indf-(k-1),append,data_file
           call read_harmpi_data_file(data_file,rho,p,u,b,kela, &
                kelb,kelc,keld)
           t(k)=tcur
           write(6,*) 'after harmpi data: ',tcur,allocated(kela)
           call lnrf_frame(real(u%data(2)/u%data(1)),real(u%data(3)/u%data(1)), & 
                real(u%data(4)/u%data(1)),r_arr,asim,th_arr,vrl,vtl,vpl)
!           write(6,*) 'lnrf transform', size(b0_arr), size(b%data(1)), size(b0_arr((k-1)*n+1:k*n))
!           write(6,*) 'lnrf transform', size(vrl), size(vtl), size(vpl), size(u0_arr)
           ! now assign to data arrays
           rho_arr((k-1)*n+1:k*n)=rho
! conversion between internal energy and pressure
           p_arr((k-1)*n+1:k*n)=p*(gam-1.)
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
           if(eHEAT.eq.1.or.eCOND.eq.1) then
              kela_arr((k-1)*n+1:k*n)=kela
              kelb_arr((k-1)*n+1:k*n)=kelb
              kelc_arr((k-1)*n+1:k*n)=kelc
              keld_arr((k-1)*n+1:k*n)=keld
           endif
!           write(6,*) 'after vtl'
!           write(6,*) 'assign'
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
! HARD CODING FOR NOW SMALL CHANGES BETWEEN DUMPS CAN CHANGE THINGS           
        tstep=1d0
        write(6,*) 'tstep slow slight: ',tstep,tstep_test
!        write(6,*) 'after loop', maxval(abs(u*u+1.))
        if(SDUMP.eq.1) call del_harmpi_grid_data()
        deallocate(rho); deallocate(p)
        deallocate(vrl); deallocate(vtl); deallocate(vpl)
        deallocate(u); deallocate(b)
        if(eCOND.eq.1.or.eHEAT.eq.1) then
           deallocate(kela); deallocate(kelb)
           deallocate(kelc); deallocate(keld)
        endif
!        write(6,*) 'read data file', a
        write(6,*) maxval(vrl_arr**2.+vtl_arr**2.+vpl_arr**2.)
!        write(6,*) 'maxmin temp: ',minval(p_arr/rho_arr*1.67e-24*9e20/1.38e-16), &
!             maxval(p_arr/rho_arr*1.67e-24*9e20/1.38e-16)
        end subroutine  load_harmpi_data

        subroutine advance_harmpi_timestep(dtobs)
        real(8), intent(in) :: dtobs
        integer :: nupdate
        nupdate = floor(dtobs / tstep)
        toffset = toffset+dtobs-tstep*nupdate
! after a couple stpes might be accumulating offset
        nupdate=nupdate+floor(toffset/tstep)
        toffset=toffset-floor(toffset/tstep)*tstep
        indf=indf+nupdate
        write(6,*) 'advance harm timestep: ',dtobs,tstep,toffset,indf,nupdate
        call update_harmpi_data(nupdate)
        end subroutine advance_harmpi_timestep

        subroutine update_harmpi_data(nupdate)
        integer, intent(in) :: nupdate
        integer :: nshift
        if(nupdate.gt.nt) then
! this is case where we just load all new data
           call load_harmpi_data(nt)
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
           if(eCOND.eq.1.or.eHEAT.eq.1) then
              kela_arr(nshift:n*nt)=kela_arr(1:n*nt-nshift)
              kelb_arr(nshift:n*nt)=kelb_arr(1:n*nt-nshift)
              kelc_arr(nshift:n*nt)=kelc_arr(1:n*nt-nshift)
              keld_arr(nshift:n*nt)=keld_arr(1:n*nt-nshift)
           endif
           call load_harmpi_data(nupdate)
        endif
        end subroutine update_harmpi_data

        subroutine init_harmpi_data(nx,n)
        integer, intent(in) :: n,nx
        allocate(rho_arr(n)); allocate(p_arr(n)); allocate(u0_arr(n))
        allocate(vrl_arr(n)); allocate(vtl_arr(n)); allocate(vpl_arr(n))
        allocate(b0_arr(n)); allocate(br_arr(n)); allocate(bth_arr(n))
        allocate(bph_arr(n))
        allocate(x1_arr(nx)); allocate(x2_arr(nx)); allocate(r_arr(nx))
        allocate(th_arr(nx)); allocate(t(nt)); allocate(x3_arr(nx))
        allocate(ph_arr(nx))
        if(eCOND.eq.1.or.eHEAT.eq.1) then
           allocate(kela_arr(n))
           allocate(kelb_arr(n))
           allocate(kelc_arr(n))
           allocate(keld_arr(n))
        endif
        end subroutine init_harmpi_data

        subroutine del_harmpi_data()
        deallocate(rho_arr); deallocate(p_arr); deallocate(u0_arr)
        deallocate(vrl_arr); deallocate(vtl_arr); deallocate(vpl_arr)
        deallocate(b0_arr); deallocate(br_arr); deallocate(bth_arr)
        deallocate(bph_arr); deallocate(x1_arr); deallocate(x2_arr)
        deallocate(r_arr); deallocate(th_arr); deallocate(t)
        deallocate(x3_arr)
        deallocate(ph_arr)
        if(eCOND.eq.1.or.eHEAT.eq.1) then
           deallocate(kela_arr)
           deallocate(kelb_arr)
           deallocate(kelc_arr)
           deallocate(keld_arr)
        end if
        end subroutine del_harmpi_data

        subroutine init_harmpi_grid_data(nx)
          integer, intent(in) :: nx
          allocate(gdet(nx)); allocate(metric_cov(10,n))
          allocate(metric_con(10,n)); allocate(drdx(16,n))
        end subroutine init_harmpi_grid_data
        
        subroutine del_harmpi_grid_data()
          deallocate(metric_cov)
          deallocate(metric_con)
          deallocate(gdet)
          deallocate(drdx)
        end subroutine del_harmpi_grid_data

      end module fluid_model_harmpi
