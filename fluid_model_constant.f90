      module fluid_model_constant

      implicit none

      contains

        subroutine init_constant()
! dummy routine for grtrans
          return
        end subroutine init_constant

        subroutine del_constant()
! dummy routine for grtrans
          return
        end subroutine del_constant

        subroutine constant_vals(u,n,B,T)
        real, dimension(:), intent(in) :: u
        real, dimension(size(u)), intent(out) :: n,B,T
        ! should update to read a file for input values so that i can easily test different regimes
        T   = 1.0d12  ! In Kelvin
        n = 1e5
        b = 10.
        end subroutine constant_vals

      end module fluid_model_constant
