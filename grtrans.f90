
       subroutine grtrans_run(ifile,outfile)
       use omp_lib
       use grtrans_inputs
       use pgrtrans
       use grtrans, only: grtrans_driver
       use ray_trace, only: ray_set,initialize_raytrace_camera, &
            kwrite_raytrace_camera, del_raytrace_camera
       use fluid_model, only: load_fluid_model, unload_fluid_model, &
            advance_fluid_timestep, source_params, assign_source_params_type, &
            fluid_args, assign_fluid_args
       use emissivity, only: emis_params
       use geodesics, only: initialize_pixels, geokerr_args, &
         del_geokerr_args,initialize_geokerr_args, initialize_geo_tabs
       use chandra_tab24, only: load_chandra_tab24, del_chandra_tab24

       implicit none
       character(len=100), intent(in) :: outfile,ifile
!       character(len=40) :: outfile,ifile
       integer :: nextra=0, inum, gunit, i, indx, ncams, j, m, l, nparams, &
            nthreads, threadnum, iii
       real (kind=8) :: wtime
!       integer, dimension(:), allocatable :: indx
       type (ray_set), dimension(:), allocatable :: c
       type (geokerr_args) :: gargs
       type (fluid_args) :: fargs
       type (source_params), dimension(:), allocatable :: sparams
       type (emis_params) :: eparams
       character(len=20), dimension(3) :: knames,kdescs
       knames(1)='nx'; knames(2)='ny'; kdescs(1)='# x pixels'
       kdescs(2)='# y pixels'
       knames(3)='nu'; kdescs(3)='Frequency (Hz)'
       gunit=12
       call read_inputs(ifile)
       write(6,*) 'dfile: ',fdfile
! overlapping code with pgrtrans
! call grtrans_main here with correct argument list
       call grtrans_main(standard,mumin,mumax,nmu,phi0,spin,&
            uout,uin, rcut, nrotype, gridvals, nn,i1,i2,fname, dt, nt, nload, &
            nmdot, mdotmin, mdotmax,ename, mbh, nfreq, fmin, fmax, muval,&
            gmin, gmax,p1, p2, jetalpha, stype,use_geokerr, nvals, iname,&
            cflag, extra, debug,outfile,fdfile,fhfile,fgfile,fsim,fnt,findf,fnfiles,fjonfix, &
            fnw,fnfreq_tab,fnr,foffset,fdindf,fmagcrit,frspot,fr0spot,fn0spot,ftscl,frscl, &
            fwmin,fwmax,ffmin,ffmax,frmax,fsigt,ffcol,fmdot,fnscl,fnnthscl,fnnthp,fbeta, &
            fbl06,fnp,ftp,frin,frout,fthin,fthout,fphiin,fphiout,coefindx, &
            epotherargs,nepotherargs)
! if you want to use pgrtrans ivals, ab, freqs then do so before here
       call del_pgrtrans_data()
! this is now done directly in pgrtrans
!       call delete_inputs()
       return
       end subroutine grtrans_run
