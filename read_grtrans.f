      include 'camera.f'

      program read_grtrans
      use ray_trace
      implicit none
      type (ray_set) :: c
      integer :: unit,maxdim,bitpix,naxis
      integer, dimension(:) :: naxes
      integer :: pcount, gcount
      logical :: simple,extend
      integer :: status=0
      call read_fits_header(unit,maxdim,simple,bitpix,naxis,naxes,
     & pcount,gcount,extend,status)
      call read_fits_array()
      call read_fits_array()
      call read_fits_array()
      
      end program
      
