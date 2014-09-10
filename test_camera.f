      program testcamera
      
      use ray_trace
  
      integer :: nvals=3
      integer :: nx=50,ny=50
      integer :: test1,test2,test3
      double precision, dimension(50,50) :: pixtest
      double precision, dimension(3,50,50) :: pixvaltest
      character(len=20) :: outfile
      double precision, dimension(3) :: vals
      type (ray_set) :: c
      outfile='raytrace.out'
      call initialize_raytrace_camera(c,nx,ny,nvals)
      write (6,*) size(c%pixloc), size(c%pixvals)
      write(6,*) 'initialized'
      c%pixloc=0d0
      write(6,*) 'pix: ',c%pixloc(33,27:29)
      c%pixloc(33,27)=4.5
      c%pixloc(33,29)=-2.
      c%pixloc(1,1)=-8.
      write(6,*) 'pixloc: ',c%pixloc(33,27:29)
      call raytrace_camera_pixel_getnext(c)
      c%next_pixel=(/50,1/)
      call raytrace_camera_pixel_getnext(c)
      write(6,*) 'pixels: ',c%current_pixel,c%next_pixel
      vals=(/-1.,1.,7./)
      call save_raytrace_camera_pixel(c,vals)
      write(6,*) 'save: ',c%current_pixel,
     &     c%pixvals(:,c%current_pixel(1),
     &                            c%current_pixel(2))
      call write_raytrace_camera(c,12,outfile)
      write(6,*) 'write'
      call del_raytrace_camera(c)
      open(unit=12,file=outfile)
      read(12,*) test1,test2,test3
      read(12,*) pixtest
      read(12,*) pixvaltest
      write(6,*) 'test: ',test1,test2,test3
      write(6,*) pixtest(33,27:29)
      write(6,*) pixvaltest(:,50,1)
      end program
