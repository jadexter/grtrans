      module grtrans_inputs

! for future with multiple emis names
!      use strings, only: countsubstring

      implicit none

      namelist /geodata/   standard,mumin,mumax,nmu,phi0,spin, uout,uin, rcut, &
       nrotype, gridvals, nn
      namelist /fluiddata/ fname, dt, nt, nload, nmdot, mdotmin, mdotmax
      namelist /emisdata/  ename, mbh, nfreq, fmin, fmax, muval, gmin, gmax, p1, p2, jetalpha, stype
      namelist /general/   use_geokerr, nvals, iname, cflag, extra
! namelists for fluid models
      namelist /harm/  fdfile, fgfile, fhfile, fnt, fnfiles, findf, fjonfix, &
           foffset, fsim, fdindf, fmagcrit
      namelist /analytic/ fnw, fwmin, fwmax, fnfreq_tab, ffmin, ffmax, frmax, fnr, fsigt, ffcol, &
           frspot,fr0spot,fn0spot,ftscl,frscl,fmdot
      
      integer :: standard,nrotype,nro,nphi,nup,nvals,nfreq,nmu,cflag, nt,nmdot,nload,extra
      logical :: use_geokerr
      double precision :: mumax,mumin,spin,rcut,a1,a2,b1,b2,mbh,uout,uin, & 
           fmin,fmax,dt,mdotmin,mdotmax,phi0,muval,gmin,gmax,p1,p2,jetalpha
      character(len=100) :: ename,fname,iname,stype
      double precision, dimension(:), allocatable :: freqs,mdots,mu0
      double precision, dimension(4) :: gridvals
      integer, dimension(3) :: nn
      ! fluid arguments
      character(len=40) :: fdfile,fhfile,fgfile,fsim
      integer :: fnt,findf,fnfiles,fjonfix,fnw,fnfreq_tab, &
           fnr,foffset,fdindf,fmagcrit
      real(8) :: frspot,fr0spot,fn0spot,ftscl,frscl,fwmin,fwmax,ffmin, &
           ffmax,frmax,fsigt,ffcol,fmdot

      interface read_inputs
        module procedure read_inputs
      end interface
 
      contains

        subroutine read_inputs(file)
        character(len=100), intent(in) :: file
        integer :: i
          open(unit=8,file=file)
          read(8,nml=geodata)
!          write(6,*) 'n: ',n
!          read(8,nml=emisdata)
          read(8,nml=fluiddata)
!          write(6,*) 'model: ',fname
          read(8,nml=emisdata)
          read(8,nml=general)
          read(8,nml=harm)
          read(8,nml=analytic)
!          write(6,*) 'general: ',iname
          close(unit=8)
          a1=dble(gridvals(1))
          a2=dble(gridvals(2))
          b1=dble(gridvals(3))
          b2=dble(gridvals(4))
          nro=nn(1)
          nphi=nn(2)
          nup=nn(3)
!          write(6,*) 'read inputs nup: ',spin
! Compute frequencies w/ log spacing between and including fmin, fmax:
          allocate(freqs(nfreq)); allocate(mu0(nmu)); allocate(mdots(nmdot))
          if(nfreq==1) then
            freqs=(/fmin/)
          else
            freqs=fmin*exp((/(i-1,i=1,nfreq)/)*log(fmax/fmin)/(nfreq-1))
          endif
          if(nmdot==1) then
            mdots=(/mdotmin/)
          else
            mdots=mdotmin*exp((/(i-1,i=1,nmdot)/)*log(mdotmax/mdotmin)/(nmdot-1))
          endif
          if(nmu==1) then
            mu0=(/mumin/)
          else
            mu0=mumin+(mumax-mumin)/(nmu-1)*(/(i-1,i=1,nmu)/)
          endif
!          write(6,*) 'nup: ',nup,freqs,fmax,fmin,fmin*exp(log(fmax/fmin))
        end subroutine read_inputs

        subroutine delete_inputs()
! deallocate input variables
        deallocate(mu0); deallocate(freqs); deallocate(mdots)
        end subroutine delete_inputs

      end module grtrans_inputs
