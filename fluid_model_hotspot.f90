  module fluid_model_hotspot
! Hotspot model from Broderick & Loeb (2005, 2006)

  use class_four_vector
  use kerr, only: calc_rms, ledd, kerr_metric, calc_polvec
  use phys_constants, only: G, c2, pi, sigb, msun

  implicit none

  namelist /hotspot/ rspot, r0spot, n0spot, bl06
  
  real :: rspot,r0spot,n0spot,tspot
  integer :: bl06

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

    subroutine init_hotspot(ifile,rs,r0,n0,bl)
    character(len=20), intent(in), optional :: ifile
    character(len=20) :: default_ifile='hotspot.in'
    real, intent(in), optional :: rs,r0,n0
    integer, intent(in), optional :: bl
    if (present(rs)) then
       rspot = rs
       r0spot = r0
       n0spot = n0
       bl06 = bl
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
    tspot=tspot-dt
    end subroutine advance_hotspot_timestep

    subroutine hotspot_vals(x0,k0,a,n,b,u,x)
      type (four_vector), intent(in), dimension(:) :: x0,k0
      real(kind=8), intent(in) :: a
      real(kind=8), intent(out), dimension(size(x0)) :: n
      type (four_vector), intent(out), dimension(size(x0)) :: b,u
      type (four_vector), dimension(size(x0)), intent(out) :: x
      type (four_vector), dimension(size(x0)) :: dx
      real(kind=8) :: rms,dcrit,bmin
      real(kind=8), dimension(size(x0)) :: kc, lc, hc, ar, d, om, safe
      real(kind=8), dimension(1) :: omega
      type (four_vector), dimension(1) :: xspot, uspot
      real(kind=8), dimension(size(x0)) :: bmag,r,dnorm, &
           gfac,zero,one,omt,ut
      real(kind=8), dimension(size(x0),10) :: metric
      real(kind=8), dimension(1,10) :: metrics
      real(kind=8), dimension(10,1) :: tmetrics
      real(kind=8), dimension(10,size(x0)) :: tmetric
!      write(6,*) 'hotspot sizes: ',size(x0), size(u)
      zero=0d0; one=1d0; bmin=1e-6
      dcrit=8d0
      r=x0%data(2)
      x=x0
      rms=calc_rms(real(a))
      d=x%data(2)*x%data(2)-2.*x%data(2)+a*a
      lc=(rms*rms-2d0*a*sqrt(rms)+a*a)/(rms**1.5d0-2d0*sqrt(rms)+a)
      hc=(2.*x%data(2)-a*lc)/d
      ar=(x%data(2)*x%data(2)+a*a)**2.-a*a*d*sin(x%data(3))**2.
      om=2.*a*x%data(2)/ar
! calculate hotspot (constant) omega:
      omega=1d0/(r0spot**(3d0/2d0)+a)
      xspot%data(1)=0d0!tspot
      xspot%data(2)=r0spot
      xspot%data(3)=acos(0d0)
! put spot at phi=0, geodesic -pi to pi to avoid phi wrapping problems
      xspot%data(4)=0d0
! phi of spot is at 0 so need to rotate to those coords
      x%data(4) = x%data(4) - (tspot+x%data(1))*omega(1)
      x%data(4) = atan2(sin(x%data(4)),cos(x%data(4)))
      metric=(kerr_metric(dble(r),x%data(3),dble(a))); tmetric=transpose(metric)
      metrics=(kerr_metric(xspot%data(2),xspot%data(3),dble(a))); tmetrics=transpose(metrics)
      call assign_metric(x,tmetric)
      call assign_metric(xspot,tmetrics)
      dx=xspot(1)-x
      dx%data(1)=0d0
      call assign_metric(dx,dble(tmetrics(:,1)))
      uspot%data(1)=sqrt(-1d0/(metrics(1,1)+2d0*metrics(1,4)*omega(1)+metrics(1,10)* &
         omega(1)*omega(1)))
      uspot%data(4)=omega(1)*uspot%data(1)
      uspot%data(2)=0d0; uspot%data(3)=0d0
      call assign_metric(uspot,dble(tmetrics))
      dnorm=dx*dx+(dx*uspot(1))**2d0
! TEST CASE no velocity effects when bl06 < 0
      if (bl06.lt.0) omega(1)=0d0
      u%data(1)=sqrt(-1d0/(metric(:,1)+2d0*metric(:,4)*omega(1)+metric(:,10)* &
         omega(1)*omega(1)))
      u%data(4)=omega(1)*u%data(1)
      u%data(2)=0d0; u%data(3)=0d0
      n=n0spot*exp(-dnorm/2d0/rspot**2.)
! this bmag is equipartition * some factor nspot / ntot = 100 for now
      bmag=sqrt(0.1d0*8d0*acos(-1d0)*n*100d0*1.67d-24/2d0*9d20/r)
      where(x%data(2).gt.rms)
        omt=max(1d0/(x%data(2)**(3d0/2d0)+a),om)
      elsewhere
        omt=max((lc+a*hc)/(x%data(2)*x%data(2)+2d0*x%data(2)*(1d0+hc)),om)
      endwhere
      ut=sqrt(-1d0/(metric(:,1)+2d0*metric(:,4)*omt+metric(:,10)* &
         omt*omt))
      safe=metric(:,1)+2d0*metric(:,4)*omega(1)+metric(:,10)* &
         omega(1)*omega(1)
      u%data(1)=merge(u%data(1),dble(ut),safe.lt.0d0)
      u%data(4)=merge(u%data(4),omt*u%data(1),safe.lt.0d0)
! floor on bmag and zero n when far away from spot
      bmag=merge(bmag,one,dnorm/2d0/rspot**2d0.lt.dcrit)
      n=merge(n,zero,dnorm/2d0/rspot**2d0.lt.dcrit)
! change so that Omega = 0 is bl06 < 0 others are all > 0 with different values
! toroidal field model bl06 = 1
! poloidal field model bl06 = 0
! vertical field model bl06 = 2
      if(abs(bl06).eq.1) then
         gfac=1d0/sqrt((metric(:,10)*metric(:,1)-metric(:,4)*metric(:,4))* & 
           (metric(:,10)*u%data(4)*u%data(4)+u%data(1)* & 
           (2d0*metric(:,4)*u%data(4)+metric(:,1)*u%data(1))))
         b%data(1)=bmag*gfac*abs(metric(:,10)*u%data(4)+metric(:,4)*u%data(1))
         b%data(4)=-bmag*sign(1d0,metric(:,10)*u%data(4)+metric(:,4)*u%data(1)) &
           *(u%data(1)*metric(:,1)+metric(:,4)*u%data(4))*gfac
         b%data(2)=0d0
         b%data(3)=0d0
      else if(bl06.eq.0) then
! poloidal field model
         b%data(1)=0d0
         b%data(2)=0d0
! bth = bmag/sqrt(gthth)
         b%data(3)=bmag/sqrt(metric(:,8))
         b%data(4)=0d0
      else if(abs(bl06).eq.2) then
! vertical field model
         b%data(1)=0d0
         b%data(4)=0d0
! using b^2 = bmag^2, bt = 0, sqrt(grr)br/sqrt(gthth)bth = -1/tan(theta)
         b%data(3)=bmag/sqrt(metric(:,8))*sin(x0%data(3))
         b%data(2)=-bmag/sqrt(metric(:,5))*cos(x0%data(3))
! calc_polvec case of vertical field in the comoving frame
      else if(abs(bl06).eq.3) then
         b=calc_polvec(x0%data(2),cos(x0%data(3)),k0,a,asin(1d0))
      else
         write(6,*) 'ERROR: unrecognized hotspot magnetic field model ',bl06
      endif
      call assign_metric(b,tmetric); call assign_metric(u,tmetric)
     end subroutine hotspot_vals

  end module fluid_model_hotspot
