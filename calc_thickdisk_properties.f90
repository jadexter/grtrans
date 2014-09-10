
  program calc_thickdisk_properties

  use grtrans_inputs, only: read_inputs, spin
  use class_four_vector
  use fluid_model_thickdisk, only: initialize_thickdisk_model, thickdisk_vals, del_thickdisk_data, &
       calc_thickdisk_mdot,calc_thickdisk_jet_power,nx1,nx2, nfiles, magcrit, nt,indf,dindf,load_thickdisk_data,r_arr, &
       advance_thickdisk_timestep,tcur, calc_thickdisk_shellavg_properties, sim, th_arr!, calc_thickdisk_inner_edge
  use kerr, only: calc_rms!, acc_efficiency

  implicit none

  real(8), dimension(:), allocatable :: mdot_arr, jet_power_arr, disk_power_arr, mdot, &
       jet_power, disk_power, wind_power, kinetic, therm, therm_arr,kinetic_arr, wind_power_arr, &
       vr, ur, smri, pmag, pmagnowt, pgasnowt, pgas, beta, thetad, bz, rho, alphamag, omega, vr2,jet_magenergy_arr, &
       jet_vel_arr,jet_magenergy,jet_vel,jetmom,jetmom_arr,jet_totpower, &
       jet_totpower_arr, gaminf, wind_kinetic, wind_therm, wind_kinetic_arr, wind_therm_arr, bzcov, bphi, bphicov
  real(8), dimension(:,:), allocatable :: gaminf_az, pgas2d, pmag2d, jetpower2d, jetvel2d, jetmom2d, therm2d, &
       thermout2d, totpower2d,jetpower2dth,therm2dth,thermout2dth,totpower2dth
  real(8) :: rin
  real :: a
  real(8) :: tstep=2.
  character(len=40) :: fstring, ifile
  character(len=10) :: append
  integer :: rdex,hdex,i,mdex

  ifile='inputs.in'
  call read_inputs(ifile)
  a=spin
  call initialize_thickdisk_model(dble(a),0)
  rdex=minloc(abs(r_arr(1:nx1)-50.),1)
  hdex=minloc(abs(r_arr(1:nx1)-(1.+sqrt(1.-a**2.))),1)+1
  mdex=minloc(abs(r_arr(1:nx1)-3.),1)
  write(6,*) 'rdex: ',rdex,hdex

  allocate(pmag2d(nx1,nx2)); allocate(pgas2d(nx1,nx2))
  allocate(totpower2d(nx1,nx2)); allocate(jetpower2d(nx1,nx2))
  allocate(therm2d(nx1,nx2)); allocate(thermout2d(nx1,nx2))
  allocate(totpower2dth(nx1,nx2)); allocate(jetpower2dth(nx1,nx2))
  allocate(therm2dth(nx1,nx2)); allocate(thermout2dth(nx1,nx2))
  allocate(jetmom2d(nx1,nx2)); allocate(jetvel2d(nx1,nx2))

  allocate(gaminf_az(nx1,nx2)); allocate(gaminf(nx1))
  allocate(mdot(nx1)); allocate(mdot_arr(nfiles))
  allocate(jet_power(nx1)); allocate(jet_power_arr(nfiles))
  allocate(jet_magenergy(nx1)); allocate(jet_magenergy_arr(nfiles))
  allocate(jet_vel(nx1)); allocate(jet_vel_arr(nfiles))
  allocate(jetmom(nx1)); allocate(jetmom_arr(nfiles))
  allocate(wind_power(nx1)); allocate(disk_power_arr(nfiles))
  allocate(disk_power(1)); allocate(therm(nx1)); allocate(kinetic(nx1))
  allocate(wind_power_arr(nfiles)); allocate(therm_arr(nfiles))
  allocate(kinetic_arr(nfiles))
  allocate(jet_totpower(nx1)); allocate(jet_totpower_arr(nfiles))
  allocate(wind_kinetic(nx1)); allocate(wind_kinetic_arr(nfiles))
  allocate(wind_therm(nx1)); allocate(wind_therm_arr(nfiles))

  allocate(vr(nx1)); allocate(ur(nx1)); allocate(rho(nx1))
  allocate(smri(nx1)); allocate(pmag(nx1)); allocate(pgas(nx1))
  allocate(beta(nx1)); allocate(thetad(nx1)); allocate(bz(nx1))
  allocate(bzcov(nx1)); allocate(bphi(nx1)); allocate(bphicov(nx1))
  allocate(omega(nx1)); allocate(vr2(nx1)); allocate(alphamag(nx1))
  allocate(pmagnowt(nx1)); allocate(pgasnowt(nx1))

  do i=1,nfiles
! now calculate mdot, jet/disk power
     write(6,*) 'pre mdot'
     call calc_thickdisk_mdot(mdot)
     write(6,*) 'after mdot pre jet'
     call calc_thickdisk_jet_power(jet_power,wind_power,therm,kinetic,jet_magenergy, &
          jet_vel,jetmom,jet_totpower,gaminf_az,wind_kinetic,wind_therm,jetpower2d,therm2d, &
          thermout2d, jetvel2d, jetmom2d, totpower2d,jetpower2dth,therm2dth,thermout2dth, &
          totpower2dth)
     write(6,*) 'pre shell',minval(gaminf_az),maxval(gaminf_az),minval(jetmom),maxval(jetmom)
     call calc_thickdisk_shellavg_properties(vr,ur,smri,pmag,pmagnowt,pgas,pgasnowt,beta &
          ,thetad,bz,rho,vr2,alphamag,omega,bzcov,bphi,bphicov,pgas2d,pmag2d)
     write(6,*) 'after jet'
!     call calc_thickdisk_inner_edge(rin)
     rin=calc_rms(a)
! assumed thindisk efficiency from rin to \infty for spin a
!     disk_power = acc_efficiency(rin,a)
     disk_power=1./2./calc_rms(a)
     !write(6,*) 'disk power: ',disk_power
     disk_power_arr(i)=disk_power(1)
     gaminf=gaminf_az(rdex,:)
     jet_power_arr(i)=jet_power(rdex)
     jet_totpower_arr(i)=jet_totpower(rdex)
     jet_magenergy_arr(i)=jet_magenergy(rdex)
     jet_vel_arr(i)=jet_vel(rdex)
     jetmom_arr(i)=jetmom(rdex)
     wind_power_arr(i)=wind_power(rdex)
     therm_arr(i)=therm(rdex)
     kinetic_arr(i)=kinetic(rdex)
     wind_kinetic_arr(i)=wind_kinetic(rdex)
     wind_therm_arr(i)=wind_therm(rdex)
!     mdot_arr(i)=sum(mdot(hdex:hdex+9))/10.
     mdot_arr(i)=mdot(mdex)
     write(6,*) 'after assigns', mdot_arr(i)
     write(6,*) mdot(1:40)
     write(6,*) disk_power_arr(i), jet_power_arr(i), &
          wind_power_arr(i), therm_arr(i), kinetic_arr(i)

     write(append,fmt='(I04)') indf
     fstring='thickdisk_shellavg_' // trim(append) // '.out'
     open(unit=12,file=fstring)
     write(12,*) nx1, nx2
     write(12,*) vr, ur, smri, pgas, pmag, pmagnowt, beta, thetad, bz, mdot, rho, r_arr(1:nx1), &
          vr2, alphamag, omega, jet_power, therm, kinetic, jet_magenergy, jet_vel,jetmom, &
          jet_totpower,gaminf,wind_kinetic,wind_therm,wind_power,pgasnowt,bzcov,bphi,bphicov
     close(unit=12)

     fstring='thickdisk_azavg_' // trim(append) // '.out'
     open(unit=12,file=fstring)
     write(12,*) nx1, nx2
     write(12,*) r_arr(1:nx1*nx2), th_arr(1:nx1*nx2), pmag2d, pgas2d,&
          jetpower2d, therm2d, thermout2d, jetvel2d,jetmom2d,totpower2d, &
          jetpower2dth,therm2dth,thermout2dth,totpower2dth
     close(unit=12)


! advance time step
     if(i.lt.nfiles) then
        indf=indf+dindf
        write(6,*) 'indf: ',indf
        call load_thickdisk_data(nt,0)
     endif
  enddo

  call del_thickdisk_data()

! write results to file
  open(unit=12,file='thickdisk_properties.out')
  write(12,*) indf, tcur, mdot_arr, jet_power_arr, wind_power_arr, &
       disk_power_arr, therm_arr, kinetic_arr, jet_magenergy_arr, jet_vel_arr, &
       jetmom_arr,jet_totpower_arr,wind_kinetic_arr,wind_therm_arr
  close(unit=12)

  deallocate(jet_totpower); deallocate(jet_totpower_arr)
  deallocate(wind_kinetic); deallocate(wind_kinetic_arr)
  deallocate(wind_therm); deallocate(wind_therm_arr)
  deallocate(gaminf_az); deallocate(gaminf)
  deallocate(jetmom); deallocate(jetmom_arr)
  deallocate(mdot); deallocate(mdot_arr); deallocate(jet_power_arr)
  deallocate(disk_power_arr); deallocate(jet_power); deallocate(wind_power)
  deallocate(disk_power); deallocate(therm); deallocate(therm_arr)
  deallocate(kinetic); deallocate(kinetic_arr); deallocate(wind_power_arr)
  deallocate(vr); deallocate(ur); deallocate(smri)
  deallocate(pmag); deallocate(pgas); deallocate(rho)
  deallocate(beta); deallocate(thetad); deallocate(bz)
  deallocate(bzcov); deallocate(bphicov); deallocate(bphi)
  deallocate(vr2); deallocate(alphamag); deallocate(omega)
  deallocate(jet_magenergy); deallocate(jet_magenergy_arr)
  deallocate(jet_vel); deallocate(jet_vel_arr)
  deallocate(pmagnowt); deallocate(pgasnowt)
  deallocate(jetpower2d); deallocate(jetmom2d)
  deallocate(jetvel2d); deallocate(therm2d)
  deallocate(totpower2d); deallocate(totpower2dth)
  deallocate(thermout2d); deallocate(thermout2dth)
  deallocate(jetpower2dth); deallocate(therm2dth)

  end program
