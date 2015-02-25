

  program testphatdisk

  use fluid_model_phatdisk, only: phatdisk_vals, init_phatdisk, del_phatdisk, nfreq_tab, freq_tab
  use fluid_model_thindisk, only: thindisk_vals, init_thindisk
  use emissivity, only: bbemis

  real, dimension(:,:), allocatable :: fnu, fnuthin
  real(kind=8), dimension(:,:), allocatable :: K
  real, dimension(:), allocatable :: omega, T, omegathin, rthin, r, ftemp
  integer :: nr, i
  real :: a, rmax

  nr=1000; a=0.; athin=0.9; rmax=1000.
  allocate(omega(nr)); allocate(omegathin(nr))
  allocate(T(nr)); allocate(K(nr,11)); allocate(rthin(nr))
  allocate(r(nr)); allocate(ftemp(nr))

  rthin=1.2*exp((/(i,i=1,nr)/)*8./(nr-1)); r=3.*exp((/(i,i=1,nr)/)*7./(nr-1))

!  write(6,*) 'testphatdisk', r

  call init_phatdisk(a)
  
!  write(6,*) 'after init', nfreq_tab
  
  allocate(fnuthin(nr,nfreq_tab)); allocate(fnu(nr,nfreq_tab))

  call phatdisk_vals(r,a,fnu,omega)

!  write(6,*) 'omega: ',omega
!  write(6,*) 'fnu: ',fnu

!  call del_phatdisk()

!  write(6,*) 'after delete'

  call init_thindisk(athin)
  call thindisk_vals(rthin,rthin*0.,athin,T,omegathin)
  do i=1,nfreq_tab
     ftemp=freq_tab(i)
     call bbemis(dble(ftemp),dble(T),K)
     fnuthin(:,i)=K(:,1)
  enddo

  open(unit=12,file='test_phatdisk.out')
  write(12,*) nr, nfreq_tab
  write(12,*) freq_tab
  write(12,*) r
  write(12,*) fnu
  write(12,*) rthin
  write(12,*) fnuthin
  close(unit=12)


  deallocate(omega); deallocate(fnu); deallocate(omegathin)
  deallocate(T); deallocate(rthin); deallocate(fnuthin)
  deallocate(r)

  call del_phatdisk()

  end program
