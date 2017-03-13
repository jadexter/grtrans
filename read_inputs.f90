      module grtrans_inputs

! for future with multiple emis names
!      use strings, only: countsubstring

      implicit none

      namelist /geodata/   standard,mumin,mumax,nmu,phi0,spin, uout,uin, rcut, &
       nrotype, gridvals, nn, i1, i2, extra, debug
      namelist /fluiddata/ fname, dt, nt, nload, nmdot, mdotmin, mdotmax
      namelist /emisdata/  ename, mbh, nfreq, fmin, fmax, muval, gmin, gmax, p1, p2, jetalpha, &
           stype, delta, nweights, coefindx
      namelist /general/   use_geokerr, nvals, iname, cflag
! namelists for fluid models
      namelist /harm/  fdfile, fgfile, fhfile, fnt, fnfiles, findf, fjonfix, &
           foffset, fsim, fdindf, fmagcrit
      namelist /analytic/ fnw, fwmin, fwmax, fnfreq_tab, ffmin, ffmax, frmax, fnr, fsigt, ffcol, &
           frspot,fr0spot,fn0spot,ftscl,frscl,fmdot,fnscl,fnnthscl,fnnthp,fbeta,fbl06,fnp,ftp, &
           frin,frout,fthin,fthout,fphiin,fphiout
      
      integer :: standard,nrotype,nro,nphi,nup,nvals,nfreq,nmu,cflag, nt,nmdot, &
           nload,extra,i1,i2,debug
      logical :: use_geokerr
      real(kind=8) :: mumax,mumin,spin,rcut,a1,a2,b1,b2,mbh,uout,uin, & 
           fmin,fmax,dt,mdotmin,mdotmax,phi0,muval,gmin,gmax,p1,p2, &
           jetalpha, delta
      character(len=100) :: ename,fname,iname,stype
      real(kind=8), dimension(:), allocatable :: freqs,mdots,mu0
      real(kind=8), dimension(4) :: gridvals
      integer, dimension(3) :: nn
      ! fluid arguments
      character(len=300) :: fdfile,fhfile,fgfile,fsim
      integer :: fnt,findf,fnfiles,fjonfix,fnw,fnfreq_tab, &
           fnr,foffset,fdindf,fmagcrit,fbl06,nweights,nepotherargs
      integer, dimension(7) :: coefindx
      real(8) :: frspot,fr0spot,fn0spot,ftscl,frscl,fwmin,fwmax,ffmin, &
           ffmax,frmax,fsigt,ffcol,fmdot,fnscl,fnnthscl,fnnthp,fbeta,fnp,ftp, &
           frin,frout,fthin,fthout,fphiin,fphiout
      real(kind=8), dimension(:), allocatable :: epotherargs
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
! new stuff for epotherargs:
          if(nweights==1) then
             allocate(epotherargs(2))
             epotherargs(1)=delta
             epotherargs(2)=1.
          else
             allocate(epotherargs(7))
             epotherargs(1)=delta
             epotherargs(2:7)=(/0.10714295,0.4459917,1.34219997, &
                  2.17450422,1.36515333,0.52892829/)
          endif
          nepotherargs=nweights+1
!          write(6,*) 'nup: ',nup,freqs,fmax,fmin,fmin*exp(log(fmax/fmin))
        end subroutine read_inputs

!        subroutine delete_inputs()
! deallocate input variables
!        deallocate(mu0); deallocate(freqs); deallocate(mdots); deallocate(epotherargs)
!        end subroutine delete_inputs

      end module grtrans_inputs
