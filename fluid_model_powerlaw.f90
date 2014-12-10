  module fluid_model_powerlaw
! Semi-analytic RIAF model from Yuan, Quataert, Narayan (2003)
! And Broderick+ (2009)

  use phys_constants, only: c2, pi, mp

  implicit none

  real :: pl_pnth, pl_n0, pl_t0, pl_nnth0, pl_beta, pl_pn, pl_pt
  
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

    subroutine init_powerlaw(n0,t0,nnth0,pnth,beta,pn,pt)
! assign inputs, defaults for one zone model w/ Gaussian height
      real, intent(in), optional :: n0,t0,nnth0,pnth,beta,pn,pt
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
      write(6,*) 'fluid model powerlaw inputs: ',pl_pnth,pl_n0,pl_t0,pl_nnth0,pl_beta,pl_pn,pl_pt
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
      neth = pl_n0 * ((rs)**(-pl_pn))*exp(-.5*(z/a2)**2.)
      nenth = pl_nnth0 * ((rs)**(-pl_pnth))*exp(-.5*(z/a2)**2.)
      te = pl_t0 * (rs)**(-pl_pt)
      ! roughly equipartition from B+09
      B = sqrt(8.*pi* neth * mp * c2 / rs / 12. / pl_beta)
      vr = 0.*r
      vth = 0.*r
      ! stationary vel to mimic non-rel
      omega = 0.*r
!      write(6,*) 'powerlaw_vals vals: ',r,z,neth,nenth,te
     end subroutine powerlaw_vals

     subroutine del_powerlaw
       !nothing so far
     end subroutine del_powerlaw

  end module fluid_model_powerlaw
