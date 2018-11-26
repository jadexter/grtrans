       module fits

       implicit none

       interface create_fits
         module procedure create_fits
       end interface
 
       interface delete_fits
         module procedure delete_fits
       end interface

       interface open_fits
         module procedure open_fits
       end interface
  
       interface close_fits
         module procedure close_fits
       end interface

       interface write_fits_array
         module procedure write_fits_array_int
         module procedure write_fits_image_int
         module procedure write_fits_array_dp
         module procedure write_fits_image_dp
         module procedure write_fits_array_real
         module procedure write_fits_image_real
       end interface

       interface read_fits_array
         module procedure read_fits_array_int
         module procedure read_fits_image_int
         module procedure read_fits_array_dp
         module procedure read_fits_image_dp
         module procedure read_fits_array_real
         module procedure read_fits_image_real
       end interface

       interface new_fits_extension
         subroutine ftcrhd(unit,status)
         integer, intent(in) :: unit
         integer, intent(inout) :: status
         end subroutine ftcrhd
       end interface

       interface next_fits_extension
         
       end interface

       interface write_fits_extra_keyword
         module procedure write_fits_extra_keyword_int
         subroutine ftpkye(unit,kname,kval,dec,kdesc,status)
         integer, intent(in) :: unit, dec
         integer, intent(inout) :: status
         real, intent(in) :: kval
         character(len=20), intent(in) :: kname,kdesc
         end subroutine ftpkye
         subroutine ftpkyd(unit,kname,kval,dec,kdesc,status)
         integer, intent(in) :: unit, dec
         real(kind=8), intent(in) :: kval
         integer, intent(inout) :: status
         character(len=20), intent(in) :: kname,kdesc
         end subroutine ftpkyd
       end interface

       interface write_fits_header
         subroutine FTPHPS(unit,bitpix,naxis,naxes,status)
         integer, intent(in) :: unit,bitpix,naxis
         integer, intent(in), dimension(naxis) :: naxes
         integer, intent(inout) :: status
         end subroutine ftphps
         subroutine FTPHPR(unit,simple,bitpix,naxis,naxes,pcount, &
         gcount,extend,status)
         integer, intent(in) :: unit,bitpix,naxis, &
                                pcount,gcount
         integer, intent(in), dimension(naxis) :: naxes
         integer, intent(inout) :: status
         logical, intent(in) :: simple,extend
         end subroutine ftphpr
       end interface

       interface read_fits_header
         subroutine FTGHPR(unit,maxdim,simple,bitpix,naxis,naxes, &
         pcount,gcount,extend,status)
         integer, intent(in) :: unit,maxdim
         integer, intent(out) :: bitpix,naxis
         integer, intent(out), dimension(maxdim) :: naxes
         integer, intent(out) :: pcount,gcount
         logical, intent(out) :: simple,extend
         integer, intent(inout) :: status
         end subroutine ftghpr
       end interface

       contains 
       
         subroutine create_fits(unit,filename,status)
         integer, intent(in) :: unit
         character(len=100), intent(in) :: filename
         integer, intent(out) :: status
         integer :: blocksize=1
         status=0
         call ftgiou(unit,status)
         call ftinit(unit,filename,blocksize,status)
         end subroutine create_fits

         subroutine delete_fits(filename,status)
         character(len=100), intent(in) :: filename
         integer, intent(inout) :: status
         integer :: unit
         call open_fits(unit,filename,1,status)
         call ftdelt(unit,status)
         call close_fits(unit,status)
         end subroutine delete_fits

         subroutine open_fits(unit,filename,readwrite,status)
         integer, intent(in) :: unit, readwrite
         character(len=100), intent(in) :: filename
         integer, intent(out) :: status
         integer :: blocksize=1
         status=0
         call ftgiou(unit,status)
         call ftopen(unit,filename,readwrite,blocksize,status)
         end subroutine open_fits

         subroutine close_fits(unit,status)
         integer, intent(in) :: unit
         integer, intent(inout) :: status         
         call ftclos(unit,status)
         call ftfiou(unit,status)
         end subroutine close_fits

         subroutine write_fits_array_int(unit,group,fpixel,nelements, &
         values,status)
         integer, intent(in) :: unit,group,fpixel,nelements
         integer, intent(inout) :: status
         integer, intent(in), dimension(:) :: values
         call ftpprj(unit,group,fpixel,nelements,values,status)
         end subroutine write_fits_array_int

         subroutine write_fits_image_int(unit,group,fpixel,nelements, &
         values,status)
         integer, intent(in) :: unit,group,fpixel,nelements
         integer, intent(inout) :: status
         integer, intent(in), dimension(:,:) :: values
         call ftpprj(unit,group,fpixel,nelements,values,status)
         end subroutine write_fits_image_int

         subroutine read_fits_array_real(unit,group,fpixel,nelements, &
         nullval, &
         values,anyf,status)
         integer, intent(in) :: unit,group,fpixel,nelements,nullval
         real, intent(out), dimension(:) :: values
         integer, intent(inout) :: status
         integer, intent(out) :: anyf
         call ftgpve(unit,group,fpixel,nelements,nullval,values,anyf, &
          status)
         end subroutine read_fits_array_real

         subroutine read_fits_image_real(unit,group,fpixel,nelements &
         ,nullval, &
         values,anyf,status)
         integer, intent(in) :: unit,group,fpixel,nelements,nullval
         real, intent(out), dimension(:,:) :: values
         integer, intent(inout) :: status
         integer, intent(out) :: anyf
         call ftgpve(unit,group,fpixel,nelements,nullval,values,anyf, &
         status)
         end subroutine read_fits_image_real

         subroutine write_fits_array_real(unit,group,fpixel,nelements, &
         values,status)
         integer, intent(in) :: unit,group,fpixel,nelements
         integer, intent(inout) :: status
         real, intent(in), dimension(:) :: values
         call ftppre(unit,group,fpixel,nelements,values,status)
         end subroutine write_fits_array_real

         subroutine write_fits_image_real(unit,group,fpixel,nelements, &
         values,status)
         integer, intent(in) :: unit,group,fpixel,nelements
         integer, intent(inout) :: status
         real, intent(in), dimension(:,:) :: values
         call ftppre(unit,group,fpixel,nelements,values,status)
         end subroutine write_fits_image_real

         subroutine read_fits_array_int(unit,group,fpixel, &
         nelements,nullval, &
         values,anyf,status)
         integer, intent(in) :: unit,group,fpixel,nelements,nullval
         integer, intent(out), dimension(:) :: values
         integer, intent(inout) :: status
         integer, intent(out) :: anyf
         call ftgpvj(unit,group,fpixel,nelements,nullval,values,anyf, &
          status)
         end subroutine read_fits_array_int

         subroutine read_fits_image_int(unit,group,fpixel, &
         nelements,nullval, &
         values,anyf,status)
         integer, intent(in) :: unit,group,fpixel,nelements,nullval
         integer, intent(out), dimension(:,:) :: values
         integer, intent(inout) :: status
         integer, intent(out) :: anyf
         call ftgpvj(unit,group,fpixel,nelements,nullval,values,anyf, &
          status)
         end subroutine read_fits_image_int

         subroutine read_fits_array_dp(unit,group,fpixel,nelements, &
          nullval, &
         values,anyf,status)
         integer, intent(in) :: unit,group,fpixel,nelements,nullval
         real(kind=8), intent(out), dimension(:) :: values
         integer, intent(inout) :: status
         integer, intent(out) :: anyf
         call ftgpvd(unit,group,fpixel,nelements,nullval,values,anyf, &
          status)
         end subroutine read_fits_array_dp

         subroutine read_fits_image_dp(unit,group,fpixel,nelements &
         ,nullval, &
         values,anyf,status)
         integer, intent(in) :: unit,group,fpixel,nelements,nullval
         real(kind=8), intent(out), dimension(:,:) :: values
         integer, intent(inout) :: status
         integer, intent(out) :: anyf
         call ftgpvd(unit,group,fpixel,nelements,nullval,values,anyf, &
          status)
         end subroutine read_fits_image_dp

         subroutine write_fits_array_dp(unit,group,fpixel,nelements, &
         values,status)
         integer, intent(in) :: unit,group,fpixel,nelements
         integer, intent(inout) :: status
         real(kind=8), intent(in), dimension(:) :: values
         call ftpprd(unit,group,fpixel,nelements,values,status)
         end subroutine write_fits_array_dp

         subroutine write_fits_image_dp(unit,group,fpixel,nelements, &
         values,status)
         integer, intent(in) :: unit,group,fpixel,nelements
         integer, intent(inout) :: status
         real(kind=8), intent(in), dimension(:,:) :: values
         call ftpprd(unit,group,fpixel,nelements,values,status)
         end subroutine write_fits_image_dp

         subroutine write_fits_extra_keyword_int(unit,kname,kval, &
         decimal,kdesc,status)
         integer, intent(in) :: unit,kval,decimal
         integer, intent(inout) :: status
         character(len=20), intent(in) :: kname,kdesc
         call ftpkyj(unit,kname,kval,kdesc,status)
         end subroutine write_fits_extra_keyword_int

         subroutine aftpkyj(unit,kname,kval, &
         kdesc,status)
         integer, intent(in) :: unit,kval
         integer, intent(inout) :: status
         character(len=20), intent(in) :: kname,kdesc
         call ftpkyj(unit,kname,kval,kdesc,status)
         end subroutine aftpkyj

         subroutine aftpkyd(unit,kname,kval, &
         decimal,kdesc,status)
         integer, intent(in) :: unit,decimal
         real(kind=8), intent(in) :: kval
         integer, intent(inout) :: status
         character(len=20), intent(in) :: kname,kdesc
         call ftpkyd(unit,kname,kval,decimal,kdesc,status)
         end subroutine aftpkyd

         subroutine aftpkys(unit,kname,kval, &
         kdesc,status)
         integer, intent(in) :: unit
         character(len=20), intent(in) :: kval
         integer, intent(inout) :: status
         character(len=20), intent(in) :: kname,kdesc
         call ftpkys(unit,kname,kval,kdesc,status)
         end subroutine aftpkys

       end module
