  module fluid_model_hotspot
! Hotspot model from Broderick & Loeb (2005, 2006)

  use class_four_vector
  use kerr, only: calc_rms, ledd, kerr_metric
  use phys_constants, only: G, c2, pi, sigb, msun

  implicit none

  namelist /hotspot/ rspot, r0spot, n0spot
  
  real :: rspot,r0spot,n0spot,tspot

  interface hotspot_vals
    module procedure hotspot_vals
  end interface

  interface init_hotspot
     module procedure init_hotspot
  end interface

  interface advance_hotspot_timestep
     module procedure advance_hotspot_timestep
  end interface

  contains

    subroutine read_hotspot_inputs(ifile)
    character(len=20), intent(in) :: ifile
    open(unit=8,file=ifile,form='formatted',status='old')
    read(8,nml=hotspot)
    close(unit=8)
!    write(6,*) 'hotspot: ',r0spot,rspot,n0spot
    end subroutine read_hotspot_inputs

    subroutine init_hotspot(ifile,rs,r0,n0)
    character(len=20), intent(in), optional :: ifile
    character(len=20) :: default_ifile='hotspot.in'
    real, intent(in), optional :: rs,r0,n0
    if (present(rs)) then
       rspot = rs
       r0spot = r0
       n0spot = n0
    else
       if (present(ifile)) then
          call read_hotspot_inputs(ifile)
       else
          call read_hotspot_inputs(default_ifile)
       endif
    endif
    tspot=0d0
    end subroutine init_hotspot

    subroutine advance_hotspot_timestep(dt)
    real, intent(in) :: dt
    tspot=tspot+dt
    end subroutine advance_hotspot_timestep

    subroutine hotspot_vals(x0,a,n,b,u,x)
      type (four_vector), intent(in), dimension(:) :: x0
      real, intent(in) :: a
      real, intent(out), dimension(size(x0)) :: n
      type (four_vector), intent(out), dimension(size(x0)) :: b,u
      type (four_vector), dimension(size(x0)), intent(out) :: x
      type (four_vector), dimension(size(x0)) :: dx
      real :: rms,dcrit,bmin
      real, dimension(size(x0)) :: kc, lc, hc, ar, d, om, safe
      real, dimension(1) :: omega, dterm1
      type (four_vector), dimension(1) :: xspot, uspot
      real, dimension(size(x0)) :: bmag,r,dnorm, &
           gfac,dterm2,dterm3,dterm4,zero,one,omt,ut
      real, dimension(size(x0),10) :: metric
      real, dimension(1,10) :: metrics
      real, dimension(10,1) :: tmetrics
      real, dimension(10,size(x0)) :: tmetric
!      write(6,*) 'hotspot sizes: ',size(x0), size(u)
      zero=0d0; one=1d0; bmin=1e-6
      dcrit=8d0
      r=x0%data(2)
      x=x0
      rms=calc_rms(a)
      d=x%data(2)*x%data(2)-2.*x%data(2)+a*a
      lc=(rms*rms-2.*a*sqrt(rms)+a*a)/(rms**1.5-2.*sqrt(rms)+a)
      hc=(2.*x%data(2)-a*lc)/d
      ar=(x%data(2)*x%data(2)+a*a)**2.-a*a*d*sin(x%data(3))**2.
      om=2.*a*x%data(2)/ar
! calculate hotspot (constant) omega:
!      where(r0spot.gt.rms)
        omega=1./(r0spot**(3./2.)+a)
!      elsewhere
!        omega=max((lc+a*hc)/(r0spot*r0spot+2.*r0spot*(1.+hc)),om)
!      endwhere
!        write(6,*) 'fluid hotspot assign xspot', omega, tspot, r0spot
      xspot%data(1)=0d0!tspot
      xspot%data(2)=r0spot
      xspot%data(3)=acos(0.)
! phi of spot is at 0 so need to rotate to those coords
      x%data(4) = x%data(4) - (tspot+x%data(1))*omega(1)
!      write(6,*) 'x1: ',x(1)%data(1),x(100)%data(1)
! make sure spot phi is between -pi, pi
      x%data(4)=atan2(sin(x%data(4)),cos(x%data(4)))
!      write(6,*) 'ps: ',x(10)%data(2),x(10)%data(4)
! put spot at phi=0
!      xspot%data(4)=atan(sin(tspot*omega),cos(tspot*omega))
      xspot%data(4)=0
!      write(6,*) 'fluid hotspot metrics'
      metric=real(kerr_metric(dble(r),x%data(3),dble(a))); tmetric=transpose(metric)
      metrics=real(kerr_metric(xspot%data(2),xspot%data(3),dble(a))); tmetrics=transpose(metrics)
!      write(6,*) 'fluid hotspot size x', size(x), size(xspot), size(u)
      call assign_metric(x,dble(tmetric))
!      write(6,*) 'fluid hotspot size x',size(x), size(xspot), size(u)
      call assign_metric(xspot,dble(tmetrics))
!      write(6,*) 'fluid hotspot u4'
!      write(6,*) 'fluid hotspot dnorm'
      dterm1 = xspot(1)*xspot(1)
      dterm2 = x*x
      dterm3 = x*xspot(1)
      dterm4 = xspot(1)*x
      dx=xspot(1)-x
      dx%data(1)=0d0
!      write(6,*) 'dx2: ',dx%data(2)
!      write(6,*) 'dx3: ',dx%data(3)
!      write(6,*) 'dx4: ',dx%data(4)
      call assign_metric(dx,dble(tmetrics(:,1)))
!      write(6,*) 'dx: ',dx(1)%metric
!      write(6,*) 'dx: ',dx(2)%metric
!      write(6,*) 'dx: ',dx(3)%metric
      !dnorm=dterm1+dterm2-dterm3-dterm4!+(x*u-xspot(1)*u)**2.
      uspot%data(1)=sqrt(-1./(metrics(1,1)+2.*metrics(1,4)*omega(1)+metrics(1,10)* &
         omega(1)*omega(1)))
      uspot%data(4)=omega(1)*uspot%data(1)
      uspot%data(2)=0d0; uspot%data(3)=0d0
!      write(6,*) 'assign metric', size(u), size(tmetric,2)
      call assign_metric(uspot,dble(tmetrics))
      dnorm=dx*dx+(dx*uspot(1))**2.
!      write(6,*) 'dnorm: ',dnorm(1:3)
!      write(6,*) 'fluid hotspot coords'
      u%data(1)=sqrt(-1./(metric(:,1)+2.*metric(:,4)*omega(1)+metric(:,10)* &
         omega(1)*omega(1)))
      u%data(4)=omega(1)*u%data(1)
      u%data(2)=0d0; u%data(3)=0d0
!      write(6,*) 'assign metric', size(u), size(tmetric,2)
      call assign_metric(uspot,dble(tmetrics))
      n=n0spot*exp(-dnorm/2d0/rspot**2.)
      bmag=sqrt(0.1*8.*acos(-1.)*n*1.67e-24/2.*9e20/r)
!      write(6,*) 'fluid hotspot b'x
      ! toroidal field
      !metric=kerr_metric(x%data(2),x%data(3),a); tmetric=transpose(metric)
      where(x%data(2).gt.rms)
        omt=max(1./(x%data(2)**(3./2.)+a),om)
      elsewhere
        omt=max((lc+a*hc)/(x%data(2)*x%data(2)+2.*x%data(2)*(1.+hc)),om)
      endwhere
!        write(6,*) 'fluid hotspot assign xspot', omega, tspot, r0spot
      ut=sqrt(-1./(metric(:,1)+2.*metric(:,4)*omt+metric(:,10)* &
         omt*omt))
!      u%data(1)=ut; u%data(4)=omt*u%data(1)
!      write(6,*) 'omega: ',om,omt,omega
!      write(6,*) 'u merge: ',dnorm/2d0/rspot**2., ut, u%data(1), metric, x%data(2), x%data(3)
      safe=metric(:,1)+2.*metric(:,4)*omega(1)+metric(:,10)* &
         omega(1)*omega(1)
!      u%data(1)=merge(u%data(1),dble(ut),dnorm/2d0/rspot**2..lt.dcrit)
!      u%data(4)=merge(u%data(4),omt*u%data(1),dnorm/2d0/rspot**2..lt.dcrit)
      u%data(1)=merge(u%data(1),dble(ut),safe.lt.0d0)
      u%data(4)=merge(u%data(4),omt*u%data(1),safe.lt.0d0)
      gfac=1d0/sqrt((metric(:,10)*metric(:,1)-metric(:,4)*metric(:,4))* & 
           (metric(:,10)*u%data(4)*u%data(4)+u%data(1)* & 
           (2d0*metric(:,4)*u%data(4)+metric(:,1)*u%data(1))))
      bmag=merge(bmag,one,dnorm/2d0/rspot**2..lt.dcrit)
!      bmag=bmag+bmin
      b%data(1)=bmag*gfac*abs(metric(:,10)*u%data(4)+metric(:,4)*u%data(1))
      b%data(4)=-bmag*sign(1d0,metric(:,10)*u%data(4)+metric(:,4)*u%data(1)) &
           *(u%data(1)*metric(:,1)+metric(:,4)*u%data(4))*gfac
      b%data(2)=0d0
      b%data(3)=0d0
      call assign_metric(b,tmetric); call assign_metric(u,tmetric)
      n=merge(n,zero,dnorm/2d0/rspot**2..lt.dcrit)
!      write(6,*) 'fluid hotspot tests: ',maxval(abs(u*u+1d0)), maxval(abs(b*u)),&
!           maxval(abs(b*b-bmag*bmag)/bmag/bmag)
!      write(6,*) 'r: ',r
!      write(6,*) 'dnorm: ',dnorm
!      write(6,*) 'n: ',n
!      write(6,*) 'bmag: ',bmag
!      write(6,*) 'b1: ',b%data(1)
!      write(6,*) 'b4: ',b%data(4)
!      write(6,*) 'u1: ',u%data(1)
!!      write(6,*) 'u4: ',u%data(4)
!      write(6,*) 'ut: ',ut
!      write(6,*) 'omega: ',omega
!      write(6,*) 'omt: ',omt
!      write(6,*) 'om: ',om
!      write(6,*) 'dx: ',dx(1)%data, uspot(1)%data
!      write(6,*) 'dx: ',dx(1)*dx(1), dx(1)*uspot(1)
     end subroutine hotspot_vals

  end module fluid_model_hotspot
