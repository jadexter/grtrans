module sashacoords

  use phys_constants, only: pi

  integer :: N1,N2,N3
  real :: a
  real :: deltaphi,R0,Rout,fracdisk,fracjet,fraccorona,nu,rsjet,rhor,Rin,r0grid, &
       r0jet,rjetend,r0disk,hoverr,rmaxdisk,rdiskend,logfac,rinsidediskmax,rbr, &
       npow2,cpow2,x10,x20,x1min,x2min

  contains

    subroutine assign_inputs(N1in,N2in,N3in,ain)
      integer, intent(in) :: N1in,N2in,N3in
      real, intent(in) :: ain
      a = ain
      N1 = N1in
      N2 = N2in
      N3 = N3in
      deltaphi = 2.*pi
      R0 = 0.3
      Rout = 1e5
      
      fracdisk = 0.36
      fracjet = 0.3
      fraccorona = 1. - fracdisk - fracjet
      nu = 0.75
      rsjet = 0.
      a = 0.99
      rhor = 1.+sqrt(1.-a*a)
      Rin = 0.83*rhor
      r0grid = 5.*Rin
      r0jet = Rin
      rjetend = 5.
      r0disk = Rin
      hoverr = 0.1
      rmaxdisk=34.
      rdiskend = 300.
      logfac = 0.
      rinsidediskmax = (1.+logfac*log10(rdiskend/rmaxdisk))**(2./nu)
      rbr = 1000.
      npow2 = 4.
      cpow2 = 1.
      x10 = 2.4
      x20 = -1.+1./N2
      x1min = log(Rin-R0)
      x1br = log(rbr-R0)
    end subroutine assign_inputs

    function step(x) result(s)
      real, dimension(:), intent(in) :: x
      real, dimension(size(x)) :: s
      s = 0.5*(sign(x)+1.)
    end function step

    function fcore(x) result(fc)
      real, dimension(:), intent(in) :: x
      real, dimension(size(x)) :: fc
      fc = 0.5*((x+1.)+2./128./pi*(-140.*sin(pi*(x+1.)/2.)+10./3.*sin(3.*pi*(x+1.)/2.)+2./5.*sin(5.*pi*(x+1.)/2.)))
    end function fcore

    function f(x) result(r)
      real, dimension(:), intent(in) :: x
      real, dimension(size(x)) :: r
      !    x = array([xx])
      !    zero = zeros(size(x))
      r = merge(merge(fcore(x),x,x.lt.1.),0.,x.gt.-1.)
    end function f

    function limlin(x,x0,dx) result(ll)
      real, dimension(:), intent(in) :: x,x0,dx
      real, dimension(size(x)) :: ll
      ll = x0-dx*f(-(x-x0)/dx)
    end function limlin

    function minlin(x,x0,dx,y0) result(ml)
      real, dimension(:), intent(in) :: x,x0,dx,y0
      real, dimension(size(x)) :: ml
      ml = y0+dx*f((x-x0)/dx)
    end function minlin
  
    function mins(f1,f2,df) result(mms)
      real, dimension(:), intent(in) :: f1,f2,df
      real, dimension(size(f1)) :: mms
      mms = limlin(f1,f2,df)
    end function mins

    function maxs(f1,f2,df) result(ms)
      real, dimension(:), intent(in) :: f1,f2,df
      real, dimension(size(f1)) :: ms
      ms = -limlin(-f1,-f2,df)
    end function maxs

    function minmaxs(f1,f2,df,dir) result(mms)
      real, dimension(:), intent(in) :: f1,f2,df,dir
      real, dimension(size(f1)) :: mms
      s = step(dir)
      mms = s*maxs(f1,f2,df)+(1.-s)*mins(f1,f2,df)
    end function minmaxs

    function calc_r(x1) result(r)
      complex, dimension(:), intent(in) :: x1
      complex, dimension(size(x1)) :: r
      r = R0+exp(x1+cpow2*(step(x1-x1br)*(x1-x1br))**npow2)
    end function calc_r
    
    function rinsidedisk(x1) result(rid)
      complex, dimension(:), intent(in) :: r
      complex, dimension(size(x1)) :: rid
      complex, dimension(size(x1)) :: half,zero,one
      half = 0.5; zero = 0.; one = 1.
      r = calc_r(x1)
      rid = (1.+logfac*log10(mins(maxs(one,r/rmaxdisk,half), &
           rdiskend/rmaxdisk+zero,half*rdiskend/rmaxdisk)))**(2./nu)
    end function rinsidedisk

    function rbeforedisk(x1) result(rbd)
      complex, dimension(:), intent(in) :: x1
      complex, dimension(size(x1)) :: rbd
      complex, dimension(size(x1)) :: one
      one = 1.
      rbd = mins(calc_r(x1),r0disk*one,0.5*r0disk*one)
    end function rbeforedisk
    
    function rafterdisk(x1) result(rad)
      complex, dimension(:), intent(in) :: x1
      complex, dimension(size(x1)) :: rad
      complex, dimension(size(x1)) :: half
      half = 0.5
      rad = maxs(one,1.+(calc_r(x1)-rdiskend)*r0jet/(rjetend*r0disk*rinsidediskmax), &
           half*rdiskend*r0jet/(rjetend*r0disk*rinsidediskmax))
    end function rafterdisk

    function fakerdisk(x1,x2) result(frd)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: frd
      frd = rbeforedisk(x1)*rinsidedisk(x1)*rafterdisk(x1)
    end function fakerdisk

    function fakerjet(x1,x2) result(frj)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: frj
      complex, dimension(size(x1)) :: one
      one = 1.
      r = calc_r(x1)
      frj = mins(r,r0jet*one,one*0.5*r0jet)*maxs(one,r/rjetend,one*0.5)
    end function fakerjet

    function Ftr1(x) result(f)
      complex, dimension(:), intent(in) :: x
      complex, dimension(size(x)) :: f
      f = 1./128.*(64.+cos(5./2.*pi*(1.+x))+70.*sin(pi*x/2.)+5.*sin(3.*pi*x/2.))
    end function Ftr1
    
    function Ftr(x) result(f)
      complex, dimension(:), intent(in) :: x
      complex, dimension(size(x)) :: f
      !    x = array([xx])
      !    zero = zeros(size(x))
      !    one = zeros(size(x))
      f = merge(0.,merge(1.,Ftr1(2.*x-1.),x.gt.1.),x.lt.0.)
    end function Ftr

    function f0(x,xa,xb,r0a,r0b) result(f00)
      complex, dimension(:), intent(in) :: x,xa,xb,r0a,r0b
      complex, dimension(size(x)) :: f00
      f00 = r0a+(r0b-r0a)*Ftr((x-xa)/(xb-xa))
    end function f0
    
    function fac(x2) result(fc)
      complex, dimension(:), intent(in) :: x2
      complex, dimension(size(x2)) :: fc
      complex, dimension(size(x2)) :: one,zero
      one=1.; zero=0.
      fc = f0(abs(x2),fracdisk+zero,one-fracjet,zero,one)
    end function fac

    function faker(x1,x2) result(fkr)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: fkr
      fkr = fakerdisk(x1,x2)*(1.-fac(x2))+fakerjet(x1,x2)*fac(x2)-rsjet*Rin
    end function faker

    function prefact(x1,x2) result(pref)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: pref
      pref = (faker(x1,x2)/r0grid)**(nu/2.)
    end function prefact
    !function thori(x1,x2):
!    th = pi/2.+arctan(tan(x2*pi/2.)*prefact(x1,x2))
    !    return th

    function prefactdisk(x1,x2) result(pred)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: pred
      pred = (fakerdisk(x1,x2)/r0grid)**(nu/2.)
    end function prefactdisk

    function prefactjet(x1,x2) result(prej)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: prej
      prej = (fakerjet(x1,x2)/r0grid)**(nu/2.)
    end function prefactjet

    function thdisk(x1,x2) result(thd)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: thd
      thd = pi/2.+arctan(tan(x2*pi/2.)*prefactdisk(x1,x2))
    end function thdisk

    function thjet(x1,x2) result(thj)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: thj
      thj = pi/2.+arctan(tan(x2*pi/2.)*prefactjet(x1,x2))
    end function thjet

    function thori(x1,x2) result(thor)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: thor
      thor = thdisk(x1,x2)*(1.-fac(x2))+thjet(x1,x2)*fac(x2)-rsjet*Rin
    end function thori

    function th0(x1,x2) result(tho)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: tho
      tho = arcsin(calc_r(x10)*sin(thori(x10,x20))/calc_r(x1))
    end function th0

    function th1in(x1,x2) result(th1)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: th1
      th1 = arcsin(calc_r(x10)*sin(thori(x10,x2))/calc_r(x1))
    end function th1in

    function th2in(x1,x2) result(th2)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: th2
      th012 = th0(x1,x2)
      th2 = ((thori(x1,x2)-thori(x1,x20))/(thori(x1,0.)-thori(x1,x20)))*(thori(x1,0.)-th012)+th012
    end function th2in

    function func1(x1,x2) result(f1)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: f1
      f1=sin(thori(x1,x2))
    end function func1

    function func2(x1,x2) result(f2)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: f2
      complex, dimension(size(x1)) :: m1
      m1 = -1.
      arg1 = real(sin(th1in(x1,x2)))
      arg2 = real(sin(th2in(x1,x2)))
      arg3 = abs(real(sin(th2in(x1,m1)))-real(sin(th1in(x1,m1))))
      arg4 = x1 - x10
      f2 = minmaxs(arg1,arg2,arg3,arg4)
    end function func2

    function dfunc(x1,x2) result(df)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: df
      df = func2(x1,x2)-func1(x1,x2)
    end function dfunc

    function calc_th(x1,x2) result(th)
      complex, dimension(:), intent(in) :: x1,x2
      complex, dimension(size(x1)) :: th
      r = calc_r(x1)
      argmax = calc_r(log(0.5*exp(x10)+0.5*exp(x1min)))*(abs(dfunc(log(0.5*exp(x10)+0.5*exp(x1min)),x2)))
      th = arcsin(maxs(r*func1(x1,x2),r*func2(x1,x2),argmax)/r)
    end function calc_th
  
    function x1func(i,N,x1min,x1max) result(x1)
      integer, dimension(:), intent(in) :: i
      integer, intent(in) :: N
      real, intent(in) :: x1min,x1max
      x1 = x1min+(x1max-x1min)*(i+0.5)/N
    end function x1func

    subroutine construct_grid(x1,x2,x3,r,th,phi)
      ! r grid
      complex, dimension(:,:,:), allocatable, intent(out) :: x1,x2,x3
      complex, dimension(:), allocatable, intent(out) :: r,th,phi
      integer :: k,j,l
      complex, dimension(:), allocatable :: x1d,x2d,x3d
      integer, dimension(:), allocatable :: i
      ! grid and BH spin parameters
      call assign_inputs(288,128,64,0.99)
      allocate(x1d(N1)); allocate(x2d(N2)); allocate(x3d(N3))
      allocate(x1(N1,N2,N3)); allocate(x2(N1,N2,N3))
      allocate(x3(N1,N2,N3)); allocate(r(N1,N2,N3))
      allocate(th(N1*N2*N3)); allocate(phi(N1*N2*N3))
      allocate(i(N1))
      i = (/(k,k=1,N1)/)
      x1max = 8.25133
      !    x1max = findx1(Rout)
      x1d = x1func(i,N1,x1min,x1max)
      x2d = (arange(N2)+0.5)/N2
      x3d = (arange(N3)+0.5)/N3
      do l=1,N1
         do j=1,N2
            x3(l,j,:) = x3d
         enddo
      enddo
      do l=1,N1
         do j=1,N3
            x2(l,:,j) = x2d
         enddo
      enddo
      do l=1,N2
         do j=1,N3
            x1(:,l,j)=x1d
         enddo
      enddo
      
      r = calc_r(reshape(x1,(/N1*N2*N3/)))
      th = calc_th(reshape(x1,(/N1*N2*N3/)),reshape(x2,(/N1*N2*N3/)))
      phi = x3*2.*pi
      deallocate(x1d); deallocate(x2d); deallocate(x3d)
      deallocate(i)
    end subroutine construct_grid
  
  end module sashacoords
