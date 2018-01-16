  module fluid_model_thindisk
! Standard NT73 thin disk model

  use class_four_vector
  use kerr, only: calc_rms, krolikc, ledd
  use phys_constants, only: G, c2, pi, sigb, msun

  implicit none

  namelist /thindisk/ mbh, mdot
  
  real :: mbh, mdot, rin, rout

  interface thindisk_vals
    module procedure thindisk_vals
  end interface

  interface init_thindisk
     module procedure init_thindisk
  end interface

  contains

    subroutine read_thindisk_inputs(ifile)
    character(len=20), intent(in) :: ifile
    open(unit=8,file=ifile,form='formatted',status='old')
    read(8,nml=thindisk)
    close(unit=8)
    end subroutine read_thindisk_inputs

    subroutine init_thindisk(a,ifile,mdott,mbht,rint,routt)
    real, intent(in) :: a
    character(len=20), intent(in), optional :: ifile
    character(len=20) :: default_ifile='thindisk.in'
    real, intent(in), optional :: mdott,mbht,rint,routt
!    write(6,*) 'init thindisk spin: ',a
    if (present(mdott)) then
       mdot = mdott
       mbh = mbht
       rin = rint
       rout = routt
    else
       if (present(ifile)) then
          call read_thindisk_inputs(ifile)
       else
          call read_thindisk_inputs(default_ifile)
       endif
    endif
    end subroutine init_thindisk

    subroutine thindisk_vals(r,th,a,T,omega)
      real, intent(in), dimension(:) :: r,th
      real, intent(in) :: a
      real, intent(out), dimension(size(r)) :: T, omega
      real :: rms,T0,lbh, Mdotedd
      real, dimension(size(r)) :: b,kc,d,ar,lc,hc,om
      real, dimension(size(r),10) :: metric
      lbh=MBH*msun*G/c2
      b=1.-3./r+2.*a/r**(3./2.)
      rms=calc_rms(a)
      rin=max(rms,rin)
      kc=krolikc(r,a)
      d=r*r-2.*r+a*a
      lc=(rms*rms-2.*a*sqrt(rms)+a*a)/(rms**1.5-2.*sqrt(rms)+a)
      hc=(2.*r-a*lc)/d
      ar=(r*r+a*a)**2.-a*a*d*sin(th)**2.
      om=2.*a*r/ar
      Mdotedd=ledd(MBH)/c2
!      write(6,*) 'T0: ',pi,G,MBH,msun,Mdot,Mdotedd,lbh,sigb,real(MBH*msun)
!      write(6,*) 'mbh: ',mbh,msun,mdot,mdotedd,lbh
      T0=(3.0/8.0/pi*G*MBH*msun*Mdot*Mdotedd/lbh/lbh/lbh/sigb)**(.25)
!      write(6,*) 'T0val: ',T0
      where(r.gt.rms)
        omega=max(1./(r**(3./2.)+a),om)
      elsewhere
        omega=max((lc+a*hc)/(r*r+2.*r*(1.+hc)),om)
      endwhere
      where(r.gt.rin.and.r.lt.rout)
         T=T0*(kc/b/r**3.)**(1./4.)
      elsewhere
         T=T0/1d5
      end where
!        write(6,*) 'T: ',r,rms,T,T0,MBH,Mdot
     end subroutine thindisk_vals


  end module fluid_model_thindisk
