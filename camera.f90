
      module ray_trace

      use fits

      implicit none

      type ray_set
      integer :: nx,ny,nvals,nextra
      real, dimension(:,:), allocatable :: pixloc
      real, dimension(:,:), allocatable :: pixvals
!      integer :: current_pixel
!      integer :: next_pixel
      end type

      interface write_raytrace_camera
        module procedure write_raytrace_camera
!        module procedure write_raytrace_camera_fits
      end interface

      interface kwrite_raytrace_camera
        module procedure kwrite_raytrace_camera
!        module procedure write_raytrace_camera_fits
      end interface
 
      contains

        subroutine save_raytrace_camera_pixel(c,pixel,vals,pnum)
        type (ray_set), intent(inout) :: c
        real, intent(in), dimension(2) :: pixel
        real, intent(in), dimension(c%nvals+c%nextra) :: vals
        integer, intent(in) :: pnum
!        write(6,*) 'in save camera: ', pixel
!        write(6,*) 'in save camera: ',vals
        c%pixvals(:,pnum)=vals
        c%pixloc(:,pnum)=pixel
!        write(6,*) 'save camera end of save'
        return
        end subroutine save_raytrace_camera_pixel

!        subroutine raytrace_camera_pixel_getnext(c)
!        type (ray_set), intent(inout) :: c
!        integer :: pixnum
!        c%current_pixel=c%next_pixel
!        write(6,*) 'npix: ',c%current_pixel(2)+1,c%ny+1
!        pixnum=c%nx*c%current_pixel(2)+c%current_pixel(1)
!        c%next_pixel=(/mod(pixnum+1,c%nx),(pixnum+1)/c%ny/)
!        c%next_pixel=c%current_pixel+1
!        return
!        end subroutine raytrace_camera_pixel_getnext

        subroutine initialize_raytrace_camera(c,nx,ny,nvals,nextra)
        type (ray_set), intent(out) :: c
        integer, intent(in) :: nx,ny,nvals,nextra
        c%nx=nx
        c%ny=ny
        c%nvals=nvals
        c%nextra=nextra
!        write(6,*) 'init cam: ',nvals,nextra
        allocate(c%pixloc(2,c%nx*c%ny))
        allocate(c%pixvals(c%nvals+c%nextra,c%nx*c%ny))
        c%pixvals(:,:)=0.
        return
        end subroutine initialize_raytrace_camera

        subroutine del_raytrace_camera(c)
        type (ray_set), intent(inout) :: c
        deallocate(c%pixloc)
        deallocate(c%pixvals)
        return
        end subroutine del_raytrace_camera

        subroutine write_raytrace_camera(c,outunit,outfile,cflag, &
         cnum,ncams,nkey,knames,kdescs,kvals)
        type (ray_set), intent(in) :: c
        integer, intent(in) :: outunit, cflag, cnum, ncams
        integer, intent(in), optional :: nkey
        character(len=20), dimension(:), intent(in), optional :: &
         kdescs,knames
        real, dimension(:), intent(in), optional:: kvals
        character(len=100), intent(in) :: outfile
        integer, dimension(2) :: naxes
        integer :: unit, status, k, nkeyz
        if(cflag==1) then
! use FITS
          write(6,*) cnum, ncams, outunit, size(c%pixloc), &
       size(c%pixvals)
          write(6,*) 'FITS ', outfile
          if(cnum==1) then
            call create_fits(unit,outfile,status)
            if(status.ne.0) write(6,*) 'Create error', status
            naxes(1)=c%nx*c%ny*2
            write(6,*) 'nvals: ',c%nvals
            call write_fits_header(unit,.true.,-32,1,(/c%nx*c%ny*2/), &
            0,1,.true.,status)
            if(status.ne.0) write(6,*) 'Header error', status
            call write_fits_array(unit,1,1,c%nx*c%ny*2,c%pixloc,status)
            if(status.ne.0) write(6,*) 'Array error', status
            call new_fits_extension(unit,status)
            if(status.ne.0) write(6,*) 'Ext error', status
          else
            call open_fits(unit,outfile,1,status)
            call new_fits_extension(unit,status)
          endif
          call write_fits_header(unit,-32,1,(/c%nx*c%ny*(c%nvals+c%nextra)/), &
          status)
          if(status.ne.0) write(6,*) 'Header error 2', status, cnum, &
               c%nx*c%ny*(c%nvals*c%nextra)
          if(present(nkey)) then
            do k=1,nkey
              write(6,*) 'k: ',k,knames(k),kvals(k),kdescs(k)
              call write_fits_extra_keyword(unit,knames(k),kvals(k), &
              6,kdescs(k),status)
              if(status.ne.0) write(6,*) 'Keyword error', k, status
            enddo
          endif
!          write(6,*) 'made it'
          write(6,*) 'write_raytrace_camera: ',size(c%pixvals)
          call write_fits_array(unit,1,1,c%nx*c%ny*(c%nvals+c%nextra),c%pixvals &
                                                        ,status)
          if(status.ne.0) write(6,*) 'Array error 2', status
          call close_fits(unit,status) 
          if(status.ne.0) write(6,*) 'Close error', status
!          call delete_fits(outfile,status)
        else
          if(cnum==1) then
            open(unit=outunit,file=outfile,FORM="unformatted")
          else
            open(unit=outunit,file=outfile,form="unformatted",access='APPEND')
          endif
!          write(6) 'write: ',c%pixloc
!          write(6,*) 'vals: ',c%pixvals
          write(outunit) c%nx,c%ny,c%nvals+c%nextra
          if(present(nkey)) then
             write(outunit) nkey
             do k=1,nkey 
                write(outunit) kvals(k)
             enddo
          else
             nkeyz=0
             write(outunit) nkeyz
          endif
          write(outunit) c%pixloc
          write(outunit) c%pixvals
          close(outunit)
        endif
        return
        end subroutine write_raytrace_camera







        subroutine kwrite_raytrace_camera(c,outunit,outfile,cflag, &
         cnum,ncams,nkey,knames,kdescs,kvals, &
         standard,mumin,mumax,nmu,phi0,spin,&
         uout,uin, rcut, nrotype, gridvals, nn, &
         fname, dt, nt, nload, nmdot, mdotmin, mdotmax, &
         ename, mbh, nfreq, fmin, fmax, muval, gmin, gmax,&
         p1, p2, jetalpha, stype, &
         use_geokerr, nvals, iname, extra)


        type (ray_set), intent(in) :: c
!MORE INPUTS
        integer, intent(in) :: standard,nrotype,nvals,nfreq,nmu,cflag, nt,nmdot,nload,extra
        logical, intent(in) :: use_geokerr
        real(kind=8), intent(in) :: mumax,mumin,spin,rcut,mbh,uout,uin, & 
           fmin,fmax,dt,mdotmin,mdotmax,phi0,muval,gmin,gmax,p1,p2,jetalpha
        character(len=100), intent(in) :: ename,fname,iname,stype
        real(kind=8), dimension(4),intent(in) :: gridvals
        integer, dimension(3), intent(in) :: nn
        character(len=20) :: tch,tcha
!End INPUTS
        integer, intent(in) :: outunit, cnum, ncams
        integer, intent(in), optional :: nkey
        character(len=20), dimension(:), intent(in), optional :: &
         kdescs,knames
        real, dimension(:), intent(in), optional:: kvals
        character(len=300), intent(in) :: outfile
        integer, dimension(2) :: naxes
        integer :: unit, status, k, nkeyz
        if(cflag==1) then
! use FITS
          write(6,*) cnum, ncams, outunit, size(c%pixloc), &
       size(c%pixvals)
          write(6,*) 'FITS ', outfile
          if(cnum==1) then
            call create_fits(unit,outfile,status)
            if(status.ne.0) write(6,*) 'Create error', status
            naxes(1)=c%nx*c%ny*2
            write(6,*) 'nvals: ',c%nvals
            call write_fits_header(unit,.true.,-32,1,(/c%nx*c%ny*2/), &
            0,1,.true.,status)
            if(status.ne.0) write(6,*) 'Header error', status
            call write_fits_array(unit,1,1,c%nx*c%ny*2,c%pixloc,status)
            if(status.ne.0) write(6,*) 'Array error', status
            call new_fits_extension(unit,status)
            if(status.ne.0) write(6,*) 'Ext error', status
          else
            call open_fits(unit,outfile,1,status)
            call new_fits_extension(unit,status)
          endif
          call write_fits_header(unit,-32,1,(/c%nx*c%ny*(c%nvals+c%nextra)/), &
          status)
          if(status.ne.0) write(6,*) 'Header error 2', status, cnum, &
               c%nx*c%ny*(c%nvals*c%nextra)
          !call ftpdat(unit)
          if(present(nkey)) then
            do k=1,nkey
              write(6,*) 'k: ',k,knames(k),kvals(k),kdescs(k)
              call write_fits_extra_keyword(unit,knames(k),kvals(k), &
              6,kdescs(k),status)
              if(status.ne.0) write(6,*) 'Keyword error', k, status
            enddo
          endif
          tch = 'standard'
          call aftpkyj(unit,tch,standard,tch,status)
          tch = 'mumin'
          call aftpkyd(unit,tch,mumin,6,tch,status)
          tch = 'mumax'
          call aftpkyd(unit,tch,mumax,6,tch,status)
          tch = 'nmu'
          call aftpkyj(unit,tch,nmu,tch,status)
          tch = 'phi0'
          call aftpkyd(unit,tch,phi0,6,tch,status)
          tch = 'spin'
          call aftpkyd(unit,tch,spin,6,tch,status)
          tch = 'uout'
          call aftpkyd(unit,tch,uout,6,tch,status)
          tch = 'uin'
          call aftpkyd(unit,tch,uin,6,tch,status)
          tch = 'rcut'
          call aftpkyd(unit,tch,rcut,6,tch,status)
          tch = 'nrotype'
          call aftpkyj(unit,tch,nrotype,tch,status)
          tch = 'gridvals(1)'
          call aftpkyd(unit,tch,gridvals(1),6,tch,status)
          tch = 'gridvals(2)'
          call aftpkyd(unit,tch,gridvals(2),6,tch,status)
          tch = 'gridvals(3)'
          call aftpkyd(unit,tch,gridvals(3),6,tch,status)
          tch = 'gridvals(4)'
          call aftpkyd(unit,tch,gridvals(4),6,tch,status)
          tch = 'nn(1)'
          call aftpkyj(unit,tch,nn(1),tch,status)
          tch = 'nn(2)'
          call aftpkyj(unit,tch,nn(2),tch,status)
          tch = 'nn(3)'
          call aftpkyj(unit,tch,nn(3),tch,status)
          tch = 'fname'
          call aftpkys(unit,tch,fname,tch,status)
          tch = 'dt'
          call aftpkyd(unit,tch,dt,6,tch,status)
          tch = 'nt'
          call aftpkyj(unit,tch,nt,tch,status)
          tch = 'nload'
          call aftpkyj(unit,tch,nload,tch,status)
          tch = 'nmdot'
          call aftpkyj(unit,tch,nmdot,tch,status)
          tch = 'mdotmin'
          call aftpkyd(unit,tch,mdotmin,6,tch,status)
          tch = 'mdotmax'
          call aftpkyd(unit,tch,mdotmax,6,tch,status)
          tch = 'ename'
          call aftpkys(unit,tch,ename,tch,status)
          tch = ' mbh'
          call aftpkyd(unit,tch, mbh,6,tch,status)
          tch = 'nfreq'
          call aftpkyj(unit,tch,nfreq,tch,status)
          tch = 'fmin'
          call aftpkyd(unit,tch,fmin,6,tch,status)
          tch = 'fmax'
          call aftpkyd(unit,tch,fmax,6,tch,status)
          tch = 'muval'
          call aftpkyd(unit,tch,muval,6,tch,status)
          tch = 'gmin'
          call aftpkyd(unit,tch,gmin,6,tch,status)
          tch = 'gmax'
          call aftpkyd(unit,tch,gmax,6,tch,status)
          tch = 'p1'
          call aftpkyd(unit,tch,p1,6,tch,status)
          tch = 'p2'
          call aftpkyd(unit,tch,p2,6,tch,status)
          tch = 'jetalpha'
          call aftpkyd(unit,tch,jetalpha,6,tch,status)
          tch = 'stype'
          call aftpkys(unit,tch,stype,tch,status)
          tch = 'use_geokerr'
          if(USE_GEOKERR) then
             tcha = 'T'
          else
             tcha = 'F'
          endif
          call aftpkys(unit,tch,tcha,tch,status)
          tch = 'nvals'
          call aftpkyj(unit,tch,nvals,tch,status)
          tch = 'iname'
          call aftpkys(unit,tch,iname,tch,status)
          tch = 'cflag'
          call aftpkyj(unit,tch,cflag,tch,status)
          tch = 'extra'
          call aftpkyj(unit,tch,extra,tch,status)
!GEOKERR IS LOGICAL






!          write(6,*) 'made it'
          write(6,*) 'write_raytrace_camera: ',size(c%pixvals)
          call write_fits_array(unit,1,1,c%nx*c%ny*(c%nvals+c%nextra),c%pixvals &
                                                        ,status)
          if(status.ne.0) write(6,*) 'Array error 2', status
          call close_fits(unit,status) 
          if(status.ne.0) write(6,*) 'Close error', status
!          call delete_fits(outfile,status)
        else
          if(cnum==1) then
            open(unit=outunit,file=outfile,FORM="unformatted")
          else
            open(unit=outunit,file=outfile,form="unformatted",access='APPEND')
          endif
!          write(6) 'write: ',c%pixloc
!          write(6,*) 'vals: ',c%pixvals
          write(outunit) c%nx,c%ny,c%nvals+c%nextra
          if(present(nkey)) then
             write(outunit) nkey
             do k=1,nkey 
                write(outunit) kvals(k)
             enddo
          else
             nkeyz=0
             write(outunit) nkeyz
          endif
          write(outunit) c%pixloc
          write(outunit) c%pixvals
          close(outunit)
        endif
        return
        end subroutine kwrite_raytrace_camera

      end module ray_trace
