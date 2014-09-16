  module fluid_model_sariaf
! Semi-analytic RIAF model from Yuan, Quataert, Narayan (2003)
! And Broderick+ (2009)

!  use class_four_vector
!  use kerr, only: calc_rms, ledd, kerr_metric
  use phys_constants, only: c2, pi, mp

  implicit none

  namelist /sariaf/ rspot, r0spot, n0spot !confused
  
  real :: rspot,r0spot,n0spot,tspot
  interface read_sariaf_inputs
     module procedure read_sariaf_inputs
  end interface

  interface sariaf_vals
    module procedure sariaf_vals
  end interface

  interface init_sariaf
     module procedure init_sariaf
  end interface

  interface del_sariaf
     module procedure del_sariaf
  end interface

! FROM HOTSPOT CODE
!  interface advance_hotspot_timestep
!     module procedure advance_hotspot_timestep
!  end interface
!
  contains

    subroutine read_sariaf_inputs(ifile)
    character(len=20), intent(in) :: ifile
    open(unit=8,file=ifile,form='formatted',status='old')
    read(8,nml=sariaf) !change from nml = hotspot
    close(unit=8)
!    write(6,*) 'hotspot: ',r0spot,rspot,n0spot !change
    end subroutine read_sariaf_inputs

    subroutine init_sariaf(ifile)
!do nothing
    character(len=20), intent(in), optional :: ifile
!    character(len=20) :: default_ifile='sariaf.in'
!    if (present(ifile)) then
!       call read_sariaf_inputs(ifile)
!    else
!       call read_sariaf_inputs(default_ifile)
!    endif
!    tspot=0d0 !change
    end subroutine init_sariaf
!    FROM HOTSPOT  
!    subroutine advance_hotspot_timestep(dt)
!    real, intent(in) :: dt
!    tspot=tspot+dt
!    end subroutine advance_hotspot_timestep
!
    subroutine sariaf_vals(a,mu,riaf_u,riaf_neth,riaf_te,riaf_B,riaf_vr,riaf_vth,riaf_omega)
      real :: riaf_n0, riaf_t0,riaf_beta
      real, intent(in) :: a
      real, dimension(:), intent(in) :: riaf_u,mu
      real, dimension(size(riaf_u)) :: riaf_r, riaf_rs,riaf_z,riaf_a2
!      real, dimension(size(riaf_u)) :: riaf_rms, riaf_lambdae, riaf_delta, riaf_game, riaf_H, riaf_compare
!      real, dimension(size(riaf_u)) :: riaf_urimsarg,riaf_urims,riaf_utims,riaf_uphims
      real :: riaf_a
      real, dimension(size(mu)) :: riaf_mu
      real, dimension(size(riaf_u)), intent(out) :: riaf_neth, riaf_te, riaf_B, riaf_vr, riaf_omega, riaf_vth

!      write(6,*) 'hotspot sizes: ',size(x0), size(u)
!r,mu,a,z,a2,rms,lambdae,delta,game,H,urimsarg,urims,utims,uphims,n0,t0,beta,rs
!compare,neth,te,B,vr
!c2 and mp from phys_constants      
!parameters
      riaf_n0 = 1.0 !4.e7
      riaf_t0 = 1.0 !1.6e11
      riaf_beta = 1.0 !10.

!seems like x0 is only useful for the r.
      riaf_r = 1./riaf_u!what is geoData.r?
      riaf_rs = riaf_r/2. !scaled r/r_s
      riaf_mu = mu !what is geoData.mu?
      riaf_a = a !what is geoData.a?
      riaf_z=riaf_r*riaf_mu
      riaf_a2 = sqrt(riaf_r**2.-riaf_z**2.)
!kerr metric
!
!      riaf_Z1 = 1.+(1.-riaf_a**2.)**(1./3.)*((1.+riaf_a)**(1./3.)+(1.-a)**(1./3.))
!      riaf_Z2 = sqrt(3.*(riaf_a**2.)+riaf_Z1**2.)


!      riaf_rms = calc_rms(riaf_a)
!      riaf_lambdae = (riaf_rms**2.-2.*riaf_a*sqrt(riaf_rms)+riaf_a**2.)/(riaf_rms**(3./2.)-2.*sqrt(riaf_rms)+riaf_a)
!      riaf_delta = riaf_r**2. -2.*riaf_r + riaf_a**2.
!      riaf_game = sqrt(1.-2./3./riaf_rms)
!      riaf_H = (2.*riaf_r-riaf_a*riaf_lambdae)/riaf_delta
!      riaf_compare = 0*riaf_rms
!      where(riaf_r.lt.riaf_rms)
!         riaf_compare = 1.
!      endwhere
!      riaf_urimsarg = (riaf_rms/riaf_r - 1.)*riaf_compare 
!      riaf_urims = -1.*sqrt(2./3./riaf_rms)*(riaf_urimsarg**(3./2.))
!      riaf_utims = riaf_game*(1.+(2./riaf_r)*(1.+riaf_H))
!      riaf_uphims = (riaf_game/riaf_r**2.)*(riaf_lambdae+riaf_a*riaf_H)
!From Broderick et al. (2009)
      riaf_neth = riaf_n0 * ((riaf_rs)**(-1.1))*exp(-.5*(riaf_z/riaf_a2)**2.) !will return rho0 eq
      riaf_te = riaf_t0 * (riaf_rs)**(-0.84) !will return p0 eq
      riaf_B = sqrt(8.*pi* riaf_neth * mp * c2 / riaf_rs / 12. / riaf_beta)
      riaf_vr = 0.*riaf_r!riaf_compare*riaf_urims/riaf_utims  !what is vr?
      riaf_vth = 0.*riaf_r
      riaf_omega = 1./(riaf_r**(3./2.)+riaf_a) !*(1.-riaf_compare) + riaf_compare*riaf_uphims/riaf_utims
!B and v are four vectors now?
     end subroutine sariaf_vals

     subroutine del_sariaf
       !nothing so far
     end subroutine del_sariaf

  end module fluid_model_sariaf
