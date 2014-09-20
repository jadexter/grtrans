  module fluid_model_hotspot_schnittman
! Hotspot model from Schnittman & Bertschinger (2004)

  use class_four_vector
  use kerr, only: calc_rms, ledd, kerr_metric
  use phys_constants, only: pi

  implicit none

  namelist /hotspot/ rspot, r0spot, n0spot
  
  real :: rspot,r0spot,n0spot,tspot

  interface schnittman_hotspot_vals
    module procedure schnittman_hotspot_vals
  end interface

  interface init_schnittman_hotspot
     module procedure init_schnittman_hotspot
  end interface

  interface advance_schnittman_hotspot_timestep
     module procedure advance_schnittman_hotspot_timestep
  end interface

  contains

    subroutine read_schnittman_hotspot_inputs(ifile)
    character(len=20), intent(in) :: ifile
    open(unit=8,file=ifile,form='formatted',status='old')
    read(8,nml=hotspot)
    close(unit=8)
!    write(6,*) 'hotspot: ',r0spot,rspot,n0spot
    end subroutine read_schnittman_hotspot_inputs

    subroutine init_schnittman_hotspot(ifile,rs,r0,n0)
    character(len=20), intent(in), optional :: ifile
    character(len=20) :: default_ifile='hotspot.in'
    real, intent(in), optional :: rs,r0,n0
    if (present(rs)) then
       rspot = rs
       r0spot = r0
       n0spot = n0
    else
       if (present(ifile)) then
          call read_schnittman_hotspot_inputs(ifile)
       else
          call read_schnittman_hotspot_inputs(default_ifile)
       endif
    endif
    tspot=0d0
    end subroutine init_schnittman_hotspot

    subroutine advance_schnittman_hotspot_timestep(dt)
    real, intent(in) :: dt
    tspot=tspot+dt
    end subroutine advance_schnittman_hotspot_timestep

    subroutine schnittman_hotspot_vals(x0,a,n)
      type (four_vector), intent(in), dimension(:) :: x0
      real, intent(in) :: a
      real, intent(out), dimension(size(x0)) :: n
      real, dimension(size(x0)) :: x,y,z,d2,t,phispot,xspot, & 
           yspot,r,phi,th,zero
      real :: dcrit,omega
!      write(6,*) 'hotspot sizes: ',size(x0), size(u)
      zero=0d0
      r = x0%data(2)
      phi = x0%data(4)
      t = x0%data(1)
      th = x0%data(3)
      omega = 1d0/(r0spot**(3d0/2d0)+a)
! shift phi & t: (THIS SHOULD NOW BE DONE ELSEWHERE)
!      phi = -phi-pi/2d0
!      t = -t
! "cartesian" coords:
      x = r*sin(th)*cos(phi)
      y = r*sin(th)*sin(phi)
      z = r*cos(th)
      phispot = omega*(t+tspot)
      xspot = r0spot*cos(phispot)
      yspot = r0spot*sin(phispot)
! "distance" between geodesics & spot:
      d2 = (x-xspot)**2.+(y-yspot)**2.+z**2.
      dcrit = 16d0*rspot**2.
      n = merge(exp(-d2/2./rspot/rspot),zero,d2.lt.dcrit)
!      write(6,*) 'schnittman hotspot r: ',r
!      write(6,*) 'schnittman hotspot t: ',t
!      write(6,*) 'schnittman hotspot th: ',th
!      write(6,*) 'schnittman hotspot phi: ',phi
!      write(6,*) 'schnittman hotspot d2: ',d2
!      write(6,*) 'schnittman hotspot n: ',n
     end subroutine schnittman_hotspot_vals

  end module fluid_model_hotspot_schnittman
