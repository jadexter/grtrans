  module fluid_model_phatdisk
! Standard NT73 thin disk model

  use class_four_vector
  use kerr, only: calc_rms, krolikc, ledd
  use phys_constants, only: G, c2, pi, sigb, msun, h, kb => k
  use fluid_model_thindisk, only: init_thindisk, thindisk_vals, MBH, Mdot
  use math, only: tsum
  use interpolate, only: get_weight

  implicit none

  namelist /phatdisk/ nw, wmin, wmax, nfreq_tab, fmin, fmax, rmax, nr, sigt, fcol
  
  real :: fmin,fmax, wmin, wmax, rmax, sigt, fcol
  integer :: nw, nfreq_tab, nr
  real, dimension(:), allocatable :: freq_tab, w, r_tab, om_tab
  real(kind=8), dimension(:,:), allocatable :: fnu_tab

  interface phatdisk_vals
    module procedure phatdisk_vals
  end interface

  interface init_phatdisk
     module procedure init_phatdisk
  end interface

  contains

    subroutine read_phatdisk_inputs(ifile)
    character(len=20), intent(in) :: ifile
    integer :: i
    open(unit=8,file=ifile,form='formatted',status='old')
    read(8,nml=phatdisk)
    close(unit=8)
    write(6,*) 'read phat: ',nfreq_tab
    allocate(freq_tab(nfreq_tab)); allocate(w(nw)); allocate(r_tab(nr))
    allocate(om_tab(nr)); allocate(fnu_tab(nr,nfreq_tab))
    if(nfreq_tab==1) then
       freq_tab=(/fmin/)
    else
       freq_tab=fmin*exp((/(i-1,i=1,nfreq_tab)/)*log(fmax/fmin)/(nfreq_tab-1))
    endif
    if(nw==1) then
       w=(/wmin/)
    else
       w=wmin*exp((/(i-1,i=1,nw)/)*log(wmax/wmin)/(nw-1))
    endif
    end subroutine read_phatdisk_inputs

    subroutine init_phatdisk(a,ifile,nwt,wmint,wmaxt,nft,fmint,fmaxt, &
         nrt,sigtt,fcolt)
    real, intent(in) :: a
    character(len=20), intent(in), optional :: ifile
    character(len=20) :: default_ifile='phatdisk.in'
    integer, intent(in), optional :: nwt,nft,nrt
    real, intent(in), optional :: wmint,wmaxt,fmint,fmaxt,sigtt,fcolt
    real :: rh,l10sigt
    real, dimension(:), allocatable :: x, fw, zero, denom
    real, dimension(:), allocatable :: z, T, tsum_arr,igrand
    integer :: i,k
!    write(6,*) 'init phat'
!    write(6,*) 'init phat: ',a
    call init_thindisk(a,ifile)
    if (present(nwt)) then
       nw = nwt
       wmin = wmint
       wmax = wmaxt
       nfreq_tab = nft
       fmin = fmint
       fmax = fmaxt
       nr = nrt
       sigt = sigtt
       fcol = fcolt
    else
       if (present(ifile)) then
          call read_phatdisk_inputs(ifile)
       else
          call read_phatdisk_inputs(default_ifile)
       end if
    endif
    allocate(z(nr)); allocate(T(nr)); allocate(x(nw)); allocate(fw(nw))
    allocate(tsum_arr(nw)); allocate(zero(nw)); allocate(igrand(nw))
    allocate(denom(nw))
    zero=0.
!    write(6,*) 'after phat inputs: ',sigt,nw,a,nr,rmax
    l10sigt=log(10.)*sigt
    x=log(w)
! Set up tables of values
!    rms=calc_rms(a)
    rh=1.+sqrt(1.-a*a)
!    write(6,*) 'r_tab: ',rh,nr,maxval((/(i,i=1,nr)/)),rmax
!    write(6,*) 'array: ',(/(i,i=1,nr)/)
    r_tab=rh*exp((/(i,i=1,nr)/)/real(nr-1)*alog(rmax/rh))
!    write(6,*) 'r_tab: ',maxval(r_tab)
! Thin disk solution:
!    write(6,*) 'before thindisk_vals', size(r_tab), size(T)
    call thindisk_vals(r_tab,r_tab*0.,a,T,om_tab)
!    write(6,*) 'after thindisk: ',r_tab
!    write(6,*) 'after thindisk: ',om_tab
!    write(6,*) 'after thindisk: ',T
! Now calculate log-normal flux:
!    fw=exp(-(log(w)+2.*sigf*sigf)**2/(2.*sigf*sigf))/sigf/sqrt(2.*pi)/w
!    write(6,*) 'l10sigt: ',l10sigt**2., x
    fw=exp(-(x+l10sigt**2.)**2./l10sigt**2.)/l10sigt/sqrt(pi)
!    write(6,*) 'fw: ',fw
!    fw=exp(-(log(w)+(l10sigt)**2.)**2./l10sigt**2.)/w/l10sigt/sqrt(pi)
    write(6,*) 'freq_tab: ',fcol, x(2)-x(1), x(nw)-x(nw-1)
    do k=1,nfreq_tab
       z=h*freq_tab(k)/kb/T/fcol
!       write(6,*) 'z: ',minval(z), maxval(z)
       do i=1,nr
! Set integrand to zero where fw=0
          denom=merge(exp(z(i)/exp(x))-1.,z(i)/exp(x),z(i)/exp(x) > 1.e-4)
          igrand=merge(fw/denom,zero,fw > 0.)
          tsum_arr=tsum(dble(x),dble(igrand))
!          tsum_arr=sum((x(2)-x(1))*igrand)
          fnu_tab(i,k)=fcol**(-4.)*2.*pi*z(i)**3*(kb*fcol*T(i))**3/h/h/c2*tsum_arr(nw)
       enddo
    enddo
!    write(6,*) 'fnu_tab: ',minval(fnu_tab),maxval(fnu_tab),fnu_tab(:,nfreq_tab)
    deallocate(z); deallocate(T); deallocate(x); deallocate(fw)
    deallocate(tsum_arr); deallocate(zero)
    end subroutine init_phatdisk

    subroutine phatdisk_vals(r,a,fnu,omega)
    real, intent(in), dimension(:) :: r
    real, intent(in) :: a
    real, intent(out), dimension(size(r)) :: omega
    real, intent(out), dimension(size(r),nfreq_tab) :: fnu
    integer :: i
    integer, dimension(size(r)) :: indx
    real, dimension(size(r)) :: weight
    call get_weight(r_tab,r,indx,weight)
!    write(6,*) 'weight: ',size(r_tab),size(r)
    do i=1,nfreq_tab
       fnu(:,i)=(1.-weight)*fnu_tab(indx,i)+weight*fnu_tab(indx+1,i)
    enddo
    omega=(1.-weight)*om_tab(indx)+weight*om_tab(indx+1)
    end subroutine phatdisk_vals

    subroutine del_phatdisk()
    deallocate(w); deallocate(freq_tab)
    deallocate(fnu_tab); deallocate(r_tab); deallocate(om_tab)
    end subroutine del_phatdisk

  end module fluid_model_phatdisk
