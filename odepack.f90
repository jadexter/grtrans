!      include 'opkda2.f'
!      include 'opkda1.f'
!      include 'opkdmain.f'

      module odepack
      implicit none
      
      real(kind=8), dimension(:), allocatable :: rwork
      integer, dimension(:), allocatable :: iwork
      integer :: jt, itol, iopt, itask, lrw, liw, istate, neq
      real(kind=8) :: rtol,atol
      real(kind=8) :: tsave
!$omp threadprivate(rwork,iwork,jt,itol,iopt,itask,lrw,liw,istate,neq,rtol,atol,tsave)

      interface lsoda_init
        module procedure lsoda_init
      end interface

      interface lsoda_wrapper
        module procedure lsoda_wrapper
      end interface

      interface lsoda
        subroutine dlsoda(f,neq,y,t,tout,itol,rtol,atol,itask, &
        istate,iopt,rwork,lrw,iwork,liw,jac,jt)
        external f, jac
        integer, intent(in) :: neq, itol, itask, liw, lrw,jt,iopt
        integer, intent(inout) :: istate
        real(kind=8), intent(in) :: tout
        real(kind=8), intent(in) :: atol,rtol
        real(kind=8), intent(inout) :: t
        real(kind=8), intent(inout), dimension(neq) :: y
        integer, intent(inout), dimension(liw) :: iwork
        real(kind=8), intent(inout), dimension(lrw) :: rwork
        end subroutine dlsoda
!        module procedure lsoda_basic
      end interface

      contains

        subroutine lsoda_basic(f,y,t,oatol,ortol,jac,yout,mt,mxs,si, &
        hmin,hmax)
        external f, jac
        real(kind=8), intent(inout), dimension(:) :: y
        real(kind=8), intent(in), dimension(:) :: t
!        real(kind=8), intent(out), dimension(size(t)) :: tarr
        real(kind=8), intent(out), dimension(size(y),size(t)) :: &
                                                             yout
        real(kind=8), intent(in) :: oatol, ortol
        integer, intent(out), optional, dimension(4) :: si
        integer, intent(in), optional :: mt,mxs
        integer :: oiopt,oitask,ojt,oneq,npts,i
        real(kind=8) :: tout, t0
        real(kind=8), intent(in), optional :: hmin, hmax
        character(len=50) :: errmsg
        oneq=size(y)
        npts=size(t)
        ojt=1; oitask=1; oiopt=0
        t0=t(1)
!        write(6,*) 'init: ',y,npts,t
!        write(6,*) 't0: ',t0,oatol,ortol,oneq
        if(present(mt)) oitask=4
        if(present(mxs)) oiopt=1
!        write(6,*) 'yout: ',size(yout),y,t,oatol,ortol
        call lsoda_init(y,t0,oatol,ortol, &
        oiopt,oitask,ojt,oneq)
!        iwork=0
!        rwork=0d0
        if(present(mt)) rwork(1)=maxval(t)
        if(present(mxs)) then
          iwork(6)=mxs
          rwork(5:7)=0d0
          iwork(5)=0
          iwork(7:9)=0
        endif
        if(present(hmin)) rwork(7)=hmin
        if(present(hmax)) rwork(6)=hmax
        yout(:,1)=y
        do i=2,npts
!          write(6,*) 'sizet: ',size(t)
          tout=t(i)
!          write(6,*) 'tout: ',tout, size(y), y, jt, oneq
          call lsoda_wrapper(f,y,tout,jac,errmsg)
!          write(6,*) 'tout: ',y
!         tarr(i)=tsave
!          write(6,*) 'yout: ',i,size(yout(:,i)),size(y)
          yout(:,i)=y
!          write(6,*) 'tout: ',tout,size(y),y,(yout(1,i)-yout(1,i-1))
!     &       /(t(i)-t(i-1)),jt,oneq
!          write(6,*) 'made it'
          if(istate.lt.0) then
! error: set this pixel to zero and break
!             write(6,*) 'istate lt 0 in lsoda: ',i,istate,errmsg
             yout(:,:)=0d0
             exit
          endif
!             istate=2   ! Reset istate in case of error
        enddo
! Write some statistics:
!        write(6,*) 'Last method switch: ',rwork(15)
!        write(6,*) 'N_f: ',iwork(12)
!        write(6,*) 'N_jac: ',iwork(13)
!        write(6,*) 'Method: ',iwork(19)
        if(present(si)) then
          si(1)=iwork(12); si(2)=iwork(13); si(3)=iwork(19); si(4)=iwork(11)
        endif
        call lsoda_destroy()
        end subroutine lsoda_basic

        subroutine lsoda_wrapper(f,y,t,jac,errmsg)
        external f, jac
        real(kind=8), intent(inout), dimension(:) :: y
        real(kind=8), intent(inout) :: t
        character(len=50), optional, intent(out) :: errmsg
!        write(6,*) 'lsoda: ',neq,y(1:2),tsave,t
!        write(6,*) 'args: ',itol,rtol,atol,itask
!        write(6,*) 'more: ',istate,iopt,lrw,liw,jt
!        write(6,*) 'tsave: ',tsave, t, y
        call lsoda(f,neq,y,tsave,t,itol,rtol,atol,itask, &
        istate,iopt,rwork,lrw,iwork,liw,jac,jt)
!        write(6,*) 'tsave: ',tsave-t
        call lsoda_err(errmsg)
!        write(6,*) errmsg
!        istate=2   ! Reset istate to continue run
        tsave=t
        end subroutine lsoda_wrapper

        subroutine lsoda_init(y,t,oatol,ortol, &
        oiopt,oitask,ojt,oneq)
!        external f, jac
        real(kind=8), intent(inout), dimension(:) :: y
        real(kind=8), intent(inout) :: t
        integer, optional :: oiopt, oitask, ojt, oneq
!        integer :: oiopt, oitask, ojt
!        real(kind=8), intent(in), dimension(:) :: oatol, ortol
        real(kind=8), intent(in) :: oatol, ortol
!        character(len=50), optional, intent(out) :: errmsg
        if (present(oneq)) then
          neq=oneq
        else
          neq=size(y)
        endif
        if (present(oiopt)) then
          iopt=oiopt
        else
          iopt=0
        endif
        if (present(oitask)) then
          itask=oitask
        else
          itask=0
        endif
        if (present(ojt)) then
          jt=ojt
        else
          jt=1
        endif
!        jt=ojt; itask=oitask; iopt=oiopt
        istate=1
        itol=1
        liw=20+neq; lrw=max(20+16*neq,22+9*neq+neq*neq);
        allocate(iwork(liw)); allocate(rwork(lrw))
        rtol=ortol; atol=oatol
!        write(6,*) 'init: ',neq,y(1:2),t
!        write(6,*) 'args: ',itol,rtol,atol,itask
!        write(6,*) 'more: ',istate,iopt,lrw,liw,jt
!        call lsoda(f,neq,y,t,tout,itol,rtol,atol,itask, &
!        istate,iopt,rwork,lrw,iwork,liw,jac,jt)
!        call lsoda_err(errmsg); write(6,*) istate, errmsg
        tsave=t
        end subroutine lsoda_init
 
        subroutine lsoda_err(errmsg)
        character(len=50), optional, intent(out) :: errmsg
! Check for ISTATE errors:
        if (present (errmsg)) then
        SELECT CASE (istate)
          CASE (1)
            errmsg='LSODA not initialized before call'
          CASE (2)
            errmsg='Integration successful'
          CASE (-1)
            errmsg='Excessive work done this step'
          CASE (-2)
            errmsg='Too much accuracy requested'
          CASE (-3)
            errmsg='Illegal input detected'
          CASE (-4)
            errmsg='Inappropriate input or singularity'
          CASE (-5) 
            errmsg='Repeated convergence test failures'
          CASE (-6)
            errmsg='Error weight became zero'
          CASE (-7)
            errmsg='Work arrays too small for new method'
        END SELECT
        else
        select case (istate)
          CASE (1)
            write(6,*) 'LSODA not initialized before call'
          CASE (2)
            write(6,*) 'Integration successful'
          CASE (-1)
            write(6,*) 'Excessive work done this step'
          CASE (-2)
            write(6,*) 'Too much accuracy requested'
          CASE (-3)
            write(6,*) 'Illegal input detected'
          CASE (-4)
            write(6,*) 'Inappropriate input or singularity'
          CASE (-5) 
            write(6,*) 'Repeated convergence test failures'
          CASE (-6)
            write(6,*) 'Error weight became zero'
          CASE (-7)
            write(6,*) 'Work arrays too small for new method'
        end select
        endif
        end subroutine lsoda_err

        subroutine lsoda_destroy()
        deallocate(iwork)
        deallocate(rwork)
!        deallocate(rtol)
!        deallocate(atol)
        end subroutine lsoda_destroy

      end module
