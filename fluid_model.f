! Module for fluid_model objects that load and interpolate fluid variable data to
! four-vector positions.
       module fluid_model
       real, allocatable, dimension(:), rho_arr, p_arr, b0_arr, br_arr, bth_arr
       real, allocatable, dimension(:), bph_arr, u0_arr, vr_arr, vth_arr, vph_arr
       type fluid_model
         double precision :: rho,p,b, 
       end type fluid_model

       interface get_fluid_vars
         module procedure get_fluid_vars
       end interface

       interface initialize_fluid_model
         module procedure intialize_num_fluid_model
         module procedure initialize_anal_fluid_model
       end interface

       contains
         
         subroutine initialize_num_fluid_model(nflag,f,dheader)
         type (fluid_model), intent(inout) :: f
         character(len=20), intent(in) :: dheader
         integer, intent(in) :: nflag
!        read_fluid_header(dheader)
!        allocate(data_arrays(nel))
         return
         end subroutine initialize_num_fluid_model

         subroutine initialize_anal_fluid_model(aflag,f)
         type (fluid_model), intent(inout) :: f
         integer, intent(in) :: aflag
         return
         end subroutine initialize_anal_fluid_model

         subroutine get_fluid_vars(x,f)
         type (four_vector) :: x
         type (fluid_model) :: f
! Populate fluid variables at position x^\mu for fluid model object f
                  
         return
         end subroutine get_fluid_vars

       end module fluid_model
