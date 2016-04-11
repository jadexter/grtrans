       ! Module for radiative transfer objects including intensities
       ! and methods for computing emissivities and for 
       ! integrating the radiative transfer equation
       module class_rad_trans
       use grtrans_inputs
       implicit none
       type rad_trans
         real(kind=8), dimension(:,:), allocatable :: I
         real(kind=8), dimension(:), allocatable :: tau
         real(kind=8) :: nu
         integer :: neq,iflag,onpts,npts,ename,idex
       end type rad_trans
        
       interface initialize_rad_trans
         module procedure initialize_rad_trans
       end interface

       interface del_rad_trans
         module procedure del_rad_trans
       end interface

       contains
         
         subroutine initialize_rad_trans(r,iname,npts,neq,extra)
         type (rad_trans), intent(out) :: r
         integer, intent(in) :: npts, neq,extra
         character(len=20), intent(in) :: iname
         integer :: nextra
         if(iname=='lsoda') then
            r%iflag=0
         else if (iname=='delo') then
            r%iflag=1
         else if (iname=='formal') then
            r%iflag=2
         else if (iname=='lsodasph') then
            r%iflag=3
         endif
!         write(6,*) 'rneq: ',neq
         r%neq=neq
! Determine which integration routine we're using:
         r%onpts=npts
! Allocate intensity array including intermediate points:
         allocate(r%I(r%neq,r%onpts))
         r%I(:,:)=0d0; r%npts=r%onpts
         if(extra.eq.1) then
            if(npts.ne.1) then
               nextra=13
            else
               nextra=7
            endif
            allocate(r%tau(nextra)); r%tau(:)=0d0
         endif
!         write(6,*) 'r%tau, r%I init: ',minval(r%tau),maxval(r%tau),minval(r%I),maxval(r%I)
         r%ename=0
         return
         end subroutine initialize_rad_trans

         subroutine del_rad_trans(r)
         type (rad_trans), intent(inout) :: r
         deallocate(r%I)
         if(allocated(r%tau)) deallocate(r%tau)
         r%ename=0
         r%iflag=0
         end subroutine del_rad_trans

       end module class_rad_trans
