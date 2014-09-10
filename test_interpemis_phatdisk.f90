      include 'emis.f90'
      include 'fluid.f90'
      include 'fluid_model_thindisk.f90'
      include 'math.f90'
      include 'interpolate.f90'
      include 'kerr.f90'
      include 'class_four_vector.f90'

      program testemis

      use emissivity
      use fluid_model
      use class_four_Vector
      implicit none

      type (emis) :: ee
      type (fluid) :: f
      type (four_Vector), dimension(1) :: x
      character(len=20) :: ff='SPHACC'
!      real, dimension(100,20) :: 

!      call initialize_fluid_model(f,ff,0d0,1)

      call initialize_emissivity(ee,'',4,7,1)
      write(6,*) 'init: ',ee%npts,ee%neq
      call calc_emissivity(2.3d11,f,ee)

      write(6,*) ee%j, ee%K
      

      call get_fluid_vars(x,f)
      write(6,*) f%ncgs, f%bcgs, f%tcgs

      call del_fluid_model(f)
      call del_emissivity(ee)
      end program testemis
