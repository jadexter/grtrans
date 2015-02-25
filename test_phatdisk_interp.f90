

  program testphatdisk

!  use fluid_model_phatdisk, only: phatdisk_vals, init_phatdisk, del_phatdisk, nfreq_tab, freq_tab
  use fluid_model, only: initialize_fluid_model, get_fluid_vars, del_fluid_model, convert_fluid_vars, &
       fluid, source_params, load_fluid_model
  use fluid_model_phatdisk, only: freq_tab, nfreq_tab, fnu_tab, r_tab
  use emissivity, only: interpemis, assign_emis_vars, emis, emis_params, calc_emissivity, &
       select_emissivity_values, initialize_emissivity, emis_model, del_emissivity
  use class_four_vector

  implicit none

  real, dimension(:,:), allocatable :: fnu,anu
  real(kind=8), dimension(:,:), allocatable :: K, fnuvals
  real, dimension(:), allocatable :: omega, r
  real(kind=8), dimension(:), allocatable :: ncgs,bcgs,tcgs,freqvals, nu
  type (fluid) :: f
  type (emis) :: e
  type (four_vector), dimension(:), allocatable :: x0
  type (source_params) :: sp
  type (emis_params) :: eparams
  integer :: nr, i
  real :: a, rmax
  character(len=20) :: fname="PHATDISK", ename="INTERP"

  write(6,*) 'start'
  nr=10; a=0.; rmax=1000.
  allocate(omega(nr)); allocate(nu(nr))
  allocate(K(nr,11)); allocate(x0(nr))
  allocate(r(nr)); allocate(bcgs(nr)); allocate(tcgs(nr)); allocate(ncgs(nr))

  write(6,*) 'r'
  r=3.*exp((/(i,i=1,nr)/)*7./(nr-1))

  nu=2.41e17
  write(6,*) 'x0'
  x0%data(1)=0.
  x0%data(2)=r
  x0%data(3)=0.
  x0%data(4)=0.

!  write(6,*) 'testphatdisk', r

!  call init_phatdisk(a)
  call load_fluid_model(fname,dble(a))
  write(6,*) 'init fluid'
  call initialize_fluid_model(f,fname,dble(a),nr)
  write(6,*) 'select emis'
  call select_emissivity_values(e,ename)
  
  write(6,*) 'after init', nfreq_tab

  write(6,*) 'get fluid vars'
  call get_fluid_vars(x0,dble(a),f)

  write(6,*) 'init emis',size(f%fnu,1),size(f%fnu,2)
  call initialize_emissivity(e,nr,f%nfreq)

  write(6,*) 'convert fluid vars'
  call convert_fluid_vars(f,ncgs,bcgs,tcgs,fnuvals,freqvals,sp)
  write(6,*) 'assign emis vars',fnuvals
  call assign_emis_vars(e,ncgs,bcgs,tcgs,fnuvals,freqvals,nfreq_tab)
  write(6,*) 'assign emis vars',e%fnu
  write(6,*) 'emis model'
  call emis_model(e,eparams)

  write(6,*) 'calc emis',freq_tab,nfreq_tab
  call calc_emissivity(nu,e)

  write(6,*) 'fnu: ',size(e%fnu)
  write(6,*) 'interp: ',e%j(:,1)
  write(6,*) 'check: ',fnu_tab(:,1)

! Write to file to check

  open(unit=12,file='phatdiskinterp.out')
  write(12,*) nr, size(r_tab)
  write(12,*) r, e%j(:,1)
  write(12,*) r_tab, fnu_tab(:,1)
  close(unit=12)

!  write(6,*) 'after delete'
  deallocate(omega)
  deallocate(r)
  deallocate(ncgs); deallocate(bcgs); deallocate(tcgs)

  call del_emissivity(e)
  call del_fluid_model(f)

  end program
