      include 'emis.f'
      include 'fluid.f'
      include 'rad_trans.f'
      include 'math.f'
      include 'interpolate.f'
      program testemis

      use emissivity
      use fluid_model
      use class_rad_trans
      use class_four_Vector
      implicit none

      type (emis) :: ee
      type (fluid) :: f
      type (four_Vector), dimension(1) :: x
      character(len=20) :: ff='SPHACC'
!      real, dimension(100,20) :: 

      x(1)%data=(/1d0,5d0,0d0,1d0/)
      call initialize_fluid_model(f,ff,0d0,1)

      f%ncgs=1d10; f%tcgs=1d11; f%bcgs=100d0; f%incang=1d0
         

      call initialize_emissivity(ee,'',4,7,1)
      write(6,*) 'init: ',ee%npts,ee%neq
      call calc_emissivity(2.3d11,f,ee)

      write(6,*) ee%j, ee%K
      

      call get_fluid_vars(x,f)
      write(6,*) f%ncgs, f%bcgs, f%tcgs

      call del_fluid_model(f)
      call del_emissivity(ee)
      end program testemis
