  module fluid_model_powerlaw
! power law model with inputs in B, T, n, assumes toroidal magnetic field 
! and velocity vectors 

  use phys_constants, only: c2, pi, mp

  implicit none

  real :: pl_pnth, pl_n0, pl_t0, pl_nnth0, pl_beta, pl_pn, pl_pt, &
       pl_rin, pl_rout, pl_thin, pl_thout, pl_phiin, pl_phiout
  
  interface powerlaw_vals
    module procedure powerlaw_vals
  end interface

  interface init_powerlaw
     module procedure init_powerlaw
  end interface

  interface del_powerlaw
     module procedure del_powerlaw
  end interface

  contains

    subroutine init_powerlaw(n0,t0,nnth0,pnth,beta,pn,pt,rin,rout,thin,thout, &
         phiin,phiout)
! assign inputs, defaults for one zone model w/ Gaussian height
      real, intent(in), optional :: n0,t0,nnth0,pnth,beta,pn,pt, &
           rin,rout,thin,thout,phiin,phiout
      write(6,*) 'init powerlaw: ',present(n0),present(t0),present(nnth0),&
           present(pnth),present(beta)
      if(present(pnth)) then
         pl_pnth = pnth
      else
         pl_pnth = 0.
      endif
      if(present(pt)) then
         pl_pt = pt
      else
         pl_pt = 0.
      endif
      if(present(pn)) then
         pl_pn = pn
      else
         pl_pn = 0.
      endif
      if(present(n0)) then
         pl_n0 = n0
      else
         pl_n0 = 3e7
      endif
      if(present(t0)) then
         pl_t0 = t0
      else
         pl_t0 = 6e10
      endif
      if(present(nnth0)) then
         pl_nnth0 = nnth0
      else
         pl_nnth0 = 8e4
      endif
      if(present(beta)) then
         pl_beta = beta
      else
         pl_beta = 10.
      endif
      if(present(rin)) then
         pl_rin = rin
      else
         pl_rin = 0.
      endif
      if(present(rout)) then
         pl_rout = rout
      else
         pl_rout = 1e8
      endif
      if(present(thin)) then
         pl_thin = thin
      else
         pl_thin = -10.
      endif
      if(present(thout)) then
         pl_thout = thout
      else
         pl_thout = 10.
      endif
      if(present(phiin)) then
         pl_phiin = phiin
      else
         pl_phiin = 0.
      endif
      if(present(phiout)) then
         pl_phiout = phiout
      else
         pl_phiout = 1e4
      endif
      write(6,*) 'fluid model powerlaw inputs: ',pl_pnth,pl_n0,pl_t0,pl_nnth0,pl_beta,pl_pn,pl_pt
      write(6,*) 'fluid model powerlaw inputs cuts: ',pl_rin,pl_rout,pl_thin,pl_thout
    end subroutine init_powerlaw

    subroutine powerlaw_vals(a,mu,u,neth,te,B,vr,vth,omega, &
         nenth)
      real, intent(in) :: a
      real, dimension(:), intent(in) :: u,mu
      real, dimension(size(u)) :: r, rs,z,a2
      real, dimension(size(u)), intent(out) :: neth, te, B, vr, omega, &
           vth, nenth
      r = 1./u
      rs = r/2.
      z=r*mu
      a2 = sqrt(r**2.-z**2.)
      neth = pl_n0 * ((rs)**(-pl_pn))!*exp(-.5*(z/a2)**2.)
      nenth = pl_nnth0 * ((rs)**(-pl_pnth))!*exp(-.5*(z/a2)**2.)
      te = pl_t0 * (rs)**(-pl_pt)
      ! roughly equipartition from B+09
!      B = sqrt(8.*pi* neth * mp * c2 / rs / 12. / pl_beta)
      vr = 0.*r
      vth = 0.*r
! velocity scaling in terms of c for equatorial plane, default is stationary to mimic non-rel one zone model
      omega = pl_phiin/r
      where(r.gt.pl_rout.or.r.lt.pl_rin.or.mu.lt.pl_thin.or. &
           mu.gt.pl_thout)
         neth = 0.
         nenth = 0.
      endwhere
! changed to use constant rs for constant value of B corresponding to R = 20 M
      B = sqrt(8.*pi* neth * mp * c2 / 10. / 12. / pl_beta)
!      write(6,*) 'powerlaw_vals vals: ',r,z,neth,nenth,te
     end subroutine powerlaw_vals

     subroutine del_powerlaw
       !nothing so far
     end subroutine del_powerlaw

  end module fluid_model_powerlaw
