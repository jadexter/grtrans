  module fluid_model_toy
! toy fluid model used for Roman's code comparison tests

  use phys_constants, only: c2, pi, mp

  implicit none

  real(kind=8) :: toy_l0, toy_n0, toy_h
  integer :: toy_which

  interface toy_vals
    module procedure toy_vals
  end interface

  interface init_toy
     module procedure init_toy
  end interface

  interface del_toy
     module procedure del_toy
  end interface

  contains

    subroutine init_toy(n0,h,l0,which)
! assign inputs, which is which problem we are doing
      real(kind=8), intent(in) :: n0,h,l0
      integer, intent(in) :: which
      toy_n0 = n0
      toy_h = h
      toy_l0 = l0
      toy_which = which
      write(6,*) 'init_toy: ',toy_n0,toy_h,toy_l0,toy_which
    end subroutine init_toy

! code comparison paper equations 1 and 2
    subroutine toy_vals(a,mu,r,n,l)
      real(kind=8), intent(in) :: a
      real(kind=8), dimension(:), intent(in) :: r,mu
      real(kind=8), dimension(size(r)), intent(out) :: n,l
      real(kind=8), dimension(size(r)) :: rcyl,z,dist
      real(kind=8) :: q
      q=0.5
      rcyl = r*sqrt(1.-mu**2.)
      l=toy_l0/(1.+rcyl)*rcyl**(1.+q)
      z=toy_h*mu
!      l=0.; z=0.
      dist=((r/10.)**2.+z**2.)
      where(dist.lt.20.)
         n=toy_n0*exp(-dist/2.)
      elsewhere
         n=0.
      endwhere
!      toy_omega = 1./(toy_r**(3./2.)+toy_a)
     end subroutine toy_vals

     subroutine del_toy
       !nothing so far
     end subroutine del_toy

  end module fluid_model_toy
