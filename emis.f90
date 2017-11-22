
       module emissivity
       use math
       use polsynchemis, only: initialize_polsynchpl,polsynchpl, &
        polsynchth,del_polsynchpl,synchpl,bnu,synchemis,synchemisnoabs
       use chandra_tab24, only: interp_chandra_tab24
       use calc_maxjutt, only: calc_maxjutt_subroutine
       use calc_maxcomp, only: calc_maxcomp_subroutine
       implicit none

! Global constants for selecting emissivity
       integer,parameter :: ELAMBDA=0,ESYNCHTH=1,EPOLSYNCHTH=2, &
                  ESYNCHPL=3,EPOLSYNCHPL=4,EBREMS=5,ELINE=6, &
                  EIRON=7,EBB=8,EBBPOL=9,EINTERP=10,EFBB=11,ERHO=12, &
                  ESYNCHTHAV=13, ESYNCHTHAVNOABS=14, EHYBRIDTHPL = 20, &
                  EHYBRIDTH=21, EHYBRIDPL=22, EMAXJUTT=23, EMAXCOMP=24

       type emis
         real(kind=8), dimension(:,:), allocatable :: j,K
         real(kind=8), dimension(:), allocatable :: ncgs,bcgs, &
          tcgs,incang,p,ncgsnth,rshift,freqarr,fnu,deltapol,psipol
         real(kind=8) :: gmax,fcol
         real(kind=8), dimension(:), allocatable :: cosne,gmin
!         real(kind=8), dimension(:), allocatable :: args
         integer :: neq,nk,npts,type,nfreq
       end type emis

       type emis_params
!         real :: gmin,mu
         real :: gmax,p1,p2
         integer :: nfreq_tab
         real(8), dimension(:), allocatable :: otherargs
         real, dimension(:), allocatable :: freq_tab,gmin,mu
         integer, dimension(:), allocatable :: usecoefs
         integer, dimension(7) :: coefindx
       end type

       interface initialize_emis_params
          module procedure initialize_emis_params
       end interface

       interface del_emis_params
          module procedure del_emis_params
       end interface

       interface assign_emis_params
          module procedure assign_emis_params
       end interface

       interface emis_model 
          module procedure emis_model
!         module procedure emis_model_lambda
!         module procedure emis_model_synchpl
!         module procedure emis_model_synchth
       end interface

       interface calc_emissivity
         module procedure calc_emissivity
       end interface

       interface rotate_emis
          module procedure rotate_emis
       end interface

       interface invariant_emis
          module procedure invariant_emis
          module procedure invariant_intensity
       end interface

       interface calc_opt_depth
         module procedure calc_opt_depth
       end interface

       contains

         subroutine interpemis(nu,freqarr,jnu,anu,K)
           ! Interpolate jnu & anu tabulated over freqarr to locations nu
           ! Adapted from IDL version of interpemisnew
           ! JAD 11/29/2011
           real(kind=8), dimension(:), intent(in) :: nu,freqarr
           real(kind=8), dimension(:), intent(in) :: jnu, anu
           real(kind=8), dimension(size(jnu)) :: logj, loga,minarr
           real(kind=8), dimension(:,:), intent(inout) :: K
           real :: minval
           integer :: nf,i
           real, dimension(size(nu)) :: yy,aa,xix,xix1,xx,aslope,jslope,yix,yix1
           integer, dimension(size(nu)) :: iff,i1,zero,nf1
           nf=size(freqarr) ; minval=-500.; nf1=nf-1; zero=0.
           minarr=minval
           logj=merge(log(jnu),minarr,jnu.gt.0)
           loga=merge(log(anu),minarr,anu.gt.0)
           iff=floor((log(nu)-log(freqarr(1)))/(log(freqarr(nf))- &
                log(freqarr(1)))*(nf-1))+1
           write(6,*) 'iff: ',log(nu),log(freqarr),nf,iff
           iff = merge(merge(iff,zero,iff.gt.0),nf-1,iff.lt.(nf-1))
           i1=(/(i,i=1,size(nu))/)
           xx=real(log(nu))
           yix=logj((i1-1)*size(freqarr)+iff) ; yix1=logj((i1-1)*size(freqarr)+iff+1)
!           write(6,*) 'yix: ',(i1-1)*size(freqarr)+iff,yix,yix1
           xix1=log(freqarr(iff+1)) ; xix=log(freqarr(iff))
           jslope=(yix1-yix)/(xix1-xix)
           yy=yix+jslope*(xx-xix)
           yix=loga((i1-1)*size(freqarr)+iff) ; yix1=loga((i1-1)*size(freqarr)+iff+1)
           xix1=log(freqarr(iff+1)) ; xix=log(freqarr(iff))
           aslope=(yix1-yix)/(xix1-xix)
           aa=yix+aslope*(xx-xix)
           K(:,1)=exp(yy); K(:,5)=exp(aa)
         end subroutine interpemis

         subroutine interpemis_noabs(nu,freqarr,jnu,K)
           ! Interpolate jnu tabulated over freqarr to locations nu
           ! Adapted from IDL version of interpemisnew
           ! JAD 11/29/2011
           real(kind=8), dimension(:), intent(in) :: nu,freqarr
           real(kind=8), dimension(:), intent(in) :: jnu
           real(kind=8), dimension(:,:), intent(inout) :: K
           real(kind=8), dimension(size(jnu)) :: logj,minarr
           real(kind=8) :: minval
           integer :: nf,i
           real, dimension(size(nu)) :: yy,xix,xix1,xx,jslope,yix,yix1
!           real, dimension(size(jnu,1),size(jnu,2)) :: yix1, yix
           integer, dimension(size(nu)) :: iff,i1,zero,nf1
           nf=size(freqarr) ; minval=-500.d0; nf1=nf-1; zero=0.
           minarr=minval
           logj=merge(log(jnu),minarr,jnu.gt.0)
           iff=floor((log(nu)-log(freqarr(1)))/(log(freqarr(nf))- &
                log(freqarr(1)))*(nf-1))+1
           iff = merge(merge(iff,zero,iff.gt.0),nf-1,iff.lt.(nf-1))
           i1=(/(i,i=1,size(nu))/)
           xx=real(log(nu))
           yix=logj((i1-1)*size(freqarr)+iff) ; yix1=logj((i1-1)*size(freqarr)+iff+1)
!           write(6,*) 'yix: ',(i1-1)*size(freqarr)+iff,yix
!           write(6,*) 'yix: ',i1,iff
           xix1=log(freqarr(iff+1)) ; xix=log(freqarr(iff))
           jslope=(yix1-yix)/(xix1-xix)
           yy=yix+jslope*(xx-xix)
           K(:,1)=exp(yy)
         end subroutine interpemis_noabs

         subroutine rhoemis(rho,rshift,K)
         real(kind=8), dimension(:), intent(in) :: rho,rshift
         real(kind=8), dimension(:,:), intent(inout) :: K
         K=0d0
! extra factor of rshift to make up for frequency integration
         K(:,1)=rho*rshift
!         write(6,*) 'rhoemis: ',rho
         end subroutine rhoemis

         subroutine fbbemis(nu,T,f,K)
         real(kind=8), dimension(:), intent(in) :: nu,T
         real(kind=8), dimension(:,:), intent(inout) :: K
         real(kind=8) :: f
!         write(6,*) 'bbemis', size(K,1), size(T), size(nu)
         K=0.d0
         K(:,1)=f**(-4d0)*bnu(T*f,nu)
         end subroutine fbbemis

         subroutine bbemis(nu,T,K)
         real(kind=8), dimension(:), intent(in) :: nu,T
         real(kind=8), dimension(:,:), intent(inout) :: K
!         write(6,*) 'bbemis', size(K,1), size(T), size(nu)
         K=0.d0
         K(:,1)=bnu(T,nu)
         end subroutine bbemis

         subroutine fbbpolemis(nu,T,f,cosne,K)
         real(kind=8), dimension(:), intent(in) :: nu,T,cosne
         real(kind=8), dimension(:,:), intent(inout) :: K
         real, dimension(size(cosne)) :: interpI, interpdel
         real(kind=8) :: f
         K=0.d0
!         write(6,*) 'K: ',size(K,1),size(K,2),size(T)
         f=1.8d0
!         write(6,*) 'T: ',T,nu
         K(:,1)=f**(-4d0)*bnu(T*f,nu)
! assumes Chandrasekhar electron scattering from semi-infinite atmosphere
         call interp_chandra_tab24(real(cosne),interpI,interpdel)
         K(:,1)=K(:,1)*dble(interpI)
         K(:,2)=K(:,1)*dble(interpdel)
         end subroutine fbbpolemis

         subroutine select_emissivity_values(e,ename)
         type (emis), intent(out) :: e
         character(len=20), intent(in) :: ename
         if(ename=='POLSYNCHTH') then
            e%type=EPOLSYNCHTH
            e%neq=4
         elseif(ename=='MAXJUTT') then
            e%type=EMAXJUTT
            e%neq=4
         elseif(ename=='MAXCOMP') then
            e%type=EMAXCOMP
            e%neq=4
         elseif(ename.eq.'HYBRIDTH') then
            e%type=EHYBRIDTH 
            e%neq=4
         elseif(ename.eq.'HYBRIDPL') then
            e%type=EHYBRIDPL 
            e%neq=4
         elseif(ename=='HYBRIDTHPL') then
            e%type=EHYBRIDTHPL 
            e%neq=4
         elseif(ename=='SYNCHTHAV') then
            e%type=ESYNCHTHAV
            e%neq=1
         elseif(ename=='SYNCHTHAVNOABS') then
            e%type=ESYNCHTHAVNOABS
            e%neq=1
         elseif(ename=='POLSYNCHPL') then
           e%neq=4
           e%type=EPOLSYNCHPL
         elseif(ename=='SYNCHPL') then
           e%neq=1
           e%type=ESYNCHPL
         elseif(ename=='lambda') then
           e%neq=1
           e%type=ELAMBDA
         elseif(ename=='BB') then
           e%neq=1
           e%type=EBB
         elseif(ename=='FBB') then
           e%neq=1
           e%type=EFBB
         elseif(ename=='RHO') then
            e%neq=1
            e%type=ERHO
         elseif(ename=='BBPOL') then
           e%neq=4
           e%type=EBBPOL
         elseif(ename=='INTERP') then
           e%neq=1
           e%type=EINTERP
         elseif(ename=='INTERPPOL') then
           e%neq=4
           e%type=EINTERP
         else
           write(6,*) 'WARNING -- Emissivity not recognized!'
         endif
         e%nk=1+e%neq*(e%neq-1)/2
         end subroutine select_emissivity_values

         subroutine initialize_emissivity(e,npts,nfreq,rshift,ang,cosne)
         type (emis), intent(inout) :: e
         integer, intent(in) :: npts, nfreq
         real(kind=8), dimension(:), intent(in) :: rshift,ang,cosne
 !        write(6,*) 'init emis: ',npts,e%nk,e%neq
         allocate(e%j(npts,e%neq)); allocate(e%K(npts,e%nk))
         e%j=0d0; e%K=0d0
         allocate(e%rshift(npts)); allocate(e%gmin(npts))
         allocate(e%incang(npts)); allocate(e%cosne(npts))
         e%rshift=rshift; e%incang=ang; e%cosne=cosne
         SELECT CASE (e%type)
            
           CASE (EHYBRIDTHPL)
             allocate(e%ncgs(npts)); allocate(e%tcgs(npts))
             allocate(e%bcgs(npts))
             allocate(e%ncgsnth(npts)); allocate(e%p(npts))
             call initialize_polsynchpl(e%neq) !I guess this is right?
           CASE (EHYBRIDTH) !Thermal decomposition of the Hybrid Image
             allocate(e%ncgs(npts)); allocate(e%tcgs(npts))
             allocate(e%bcgs(npts))
             allocate(e%ncgsnth(npts)); allocate(e%p(npts))
             call initialize_polsynchpl(e%neq) !I guess this is right?
           CASE (EHYBRIDPL) !Power Law decomposition of the Hybrid Image
             allocate(e%ncgs(npts)); allocate(e%tcgs(npts))
             allocate(e%bcgs(npts))
             allocate(e%ncgsnth(npts)); allocate(e%p(npts))
             call initialize_polsynchpl(e%neq) !I guess this is right?
           CASE (EMAXJUTT) !alwinnote 2015/03/05
             allocate(e%tcgs(npts)); allocate(e%ncgs(npts))
             allocate(e%bcgs(npts))!; allocate(e%args(npts))
           CASE (EMAXCOMP) !alwinnote 2015/03/05
             allocate(e%tcgs(npts)); allocate(e%ncgs(npts))
             allocate(e%bcgs(npts))!; allocate(e%args(npts))
!             e%args=emisargs
           CASE (EPOLSYNCHTH)
             allocate(e%tcgs(npts)); allocate(e%ncgs(npts))
             allocate(e%bcgs(npts))
           CASE (EPOLSYNCHPL)
             allocate(e%ncgs(npts)); allocate(e%tcgs(npts))
             allocate(e%bcgs(npts))
             allocate(e%ncgsnth(npts)); allocate(e%p(npts))
             call initialize_polsynchpl(e%neq)
           CASE (ESYNCHPL)
!             write(6,*) 'init epl: ',allocated(e%ncgs),allocated(e%tcgs)
             allocate(e%ncgs(npts)); allocate(e%tcgs(npts))
             allocate(e%bcgs(npts))
             allocate(e%ncgsnth(npts)); allocate(e%p(npts))
             call initialize_polsynchpl(e%neq)
           CASE (ESYNCHTHAV)
             allocate(e%tcgs(npts)); allocate(e%ncgs(npts))
             allocate(e%bcgs(npts))
           CASE (ESYNCHTHAVNOABS)
             allocate(e%tcgs(npts)); allocate(e%ncgs(npts))
             allocate(e%bcgs(npts))
           CASE (EBB)
             allocate(e%tcgs(npts))
           CASE (EFBB)
              allocate(e%tcgs(npts))
           CASE (ERHO)
              allocate(e%ncgs(npts))
           CASE (EBBPOL)
             allocate(e%tcgs(npts))
             e%incang=0.
           CASE (EINTERP)
             allocate(e%fnu(npts*nfreq)); allocate(e%freqarr(nfreq))
         END SELECT
         e%npts=npts
         e%nfreq=nfreq
         end subroutine initialize_emissivity

         subroutine calc_emissivity(nu,e,ep)
         type (emis), intent(inout) :: e
         type (emis_params), intent(in) :: ep
         real(kind=8), intent(in), dimension(:) :: nu
         real(kind=8), dimension(e%npts,11) :: K, Kth, Kpl

         integer :: i
         select  case(e%type)
           case(emaxjutt)
              call calc_maxjutt_subroutine(nu,e%ncgs,e%bcgs,e%tcgs,e%incang,ep%otherargs,K)
           case(emaxcomp)
              call calc_maxcomp_subroutine(nu,e%ncgs,e%bcgs,e%tcgs,e%incang,ep%otherargs,K)
           case(ehybridthpl)
              call polsynchth(nu,e%ncgs,e%bcgs,e%tcgs,e%incang,Kth)
              call polsynchpl(nu,e%ncgsnth,e%bcgs,e%incang,e%p,e%gmin, &
               e%gmax,Kpl)
              K = Kth + Kpl
           case(ehybridth)
              call polsynchth(nu,e%ncgs,e%bcgs,e%tcgs,e%incang,Kth)
              call polsynchpl(nu,e%ncgsnth,e%bcgs,e%incang,e%p,e%gmin, &
               e%gmax,Kpl)
              K = Kth !get thermal emission
              K(:,5:11) = Kth(:,5:11) + Kpl(:,5:11) !and total absorption
           case(ehybridpl)
              call polsynchth(nu,e%ncgs,e%bcgs,e%tcgs,e%incang,Kth)
              call polsynchpl(nu,e%ncgsnth,e%bcgs,e%incang,e%p,e%gmin, &
               e%gmax,Kpl)
              K = Kpl !get power law emission
              K(:,5:11) = Kth(:,5:11) + Kpl(:,5:11) !and total absorption
           case(epolsynchth)
!              write(6,*) 'polsynchth ncgs: ',allocated(e%ncgs),size(e%ncgs)
!              write(6,*) 'polsynchth bcgs: ',e%bcgs
!              write(6,*) 'polsynchth tcgs: ',e%tcgs
             call polsynchth(nu,e%ncgs,e%bcgs,e%tcgs,e%incang,K)
           case(epolsynchpl)
!              write(6,*) 'ppl nu', nu
!             write(6,*) 'ppl n', e%ncgsnth
!             write(6,*) 'ppl b',e%bcgs
!             write(6,*) 'ppl i',e%incang
!             write(6,*) 'ppl p',e%p
             call polsynchpl(nu,e%ncgsnth,e%bcgs,e%incang,e%p,e%gmin, &
              e%gmax,K)
             if(any(isnan(K))) then
                write(6,*) 'NaN in polsynch sizes: ', &
                     e%npts,size(nu), &
                     size(e%ncgsnth),size(e%bcgs),size(e%incang),size(e%p), &
                     size(e%gmin),size(K,1),size(K,2)
                write(6,*) 'NaN in polsynchpl K8: ',K(:,8)
                write(6,*) 'NaN in polsynchpl K9: ',K(:,9)
                write(6,*) 'NaN in polsynchpl nu: ',nu
                write(6,*) 'NaN in polsynchpl ang: ',e%incang
                write(6,*) 'NaN in polsynchpl n: ',e%ncgsnth
                write(6,*) 'NaN in polsynchpl B: ',e%bcgs
                write(6,*) 'NaN in polsynchpl p: ',e%p
                write(6,*) 'NaN in polsynchpl gmin: ',e%gmin
             endif
           case(esynchpl)
             call synchpl(nu,e%ncgsnth,e%bcgs,e%incang,e%p,e%gmin, &
              e%gmax,K)
!             write(6,*) 'synchpl: ',e%ncgsnth,e%bcgs,K(:,1)
!             write(6,*) 'synchpl sum: ',sum(K(:,1))
           case(esynchthav)
              call synchemis(nu,e%ncgs,e%bcgs,e%tcgs,K)
           case(esynchthavnoabs)
              call synchemisnoabs(nu,e%ncgs,e%bcgs,e%tcgs,K)
           case(ebb)
             call bbemis(nu,e%tcgs,K)
           case(ebbpol)
!              write(6,*) 'call bbpolemis: ',size(e%tcgs),size(K,1),size(K,2), size(nu)
             call fbbpolemis(nu,e%tcgs,e%fcol,e%cosne,K)
           case(erho)
              call rhoemis(e%ncgs,e%rshift,K)
           case(einterp)
!              write(6,*) 'interpemis no abs: ',size(nu), size(e%freqarr), size(e%fnu)
!              write(6,*) 'interpemis K: ',size(K,1), size(K,2)
             call interpemis_noabs(nu,e%freqarr,e%fnu,K)
         end select
!         write(6,*) 'afterpolsynch', e%bcgs, e%ncgs, e%tcgs, K(1:e%npts,1)
         if (e%neq==4) then
           e%j(1:e%npts,1:e%neq)=K(1:e%npts,1:e%neq)
           if (allocated(ep%usecoefs)) &
                e%K(1:e%npts,ep%usecoefs)=K(1:e%npts,ep%usecoefs+e%neq)
!           e%K(1:e%npts,1:e%nk)=K(1:e%npts,e%neq+1:e%nk+e%neq)
!           if(allocated(ep%usecoefs)) then
!              e%K(1:e%npts,ep%usecoefs)=K(1:e%npts,ep%usecoefs+e%neq)
!           else
!              e%K(1:e%npts,1:e%nk)=K(1:e%npts,e%neq+1:e%nk+e%neq)
!           endif
!           write(6,*) 'pol assign e ',size(e%j),size(e%K),size(K(:,1:4)),size(K(:,5:11))
         else
           e%j(1:e%npts,1)=K(1:e%npts,1)
           e%K(1:e%npts,1)=K(1:e%npts,5)
!           write(6,*) 'assign e', e%j,e%npts
         endif
         end subroutine calc_emissivity

         subroutine assign_emis_params(e,ncgs,ncgsnth,bcgs,tcgs,fnu,freqarr,nf)
         type (emis), intent(inout) :: e
         real(kind=8), dimension(:), intent(in) :: bcgs,ncgs,tcgs,ncgsnth
         real(kind=8), dimension(:), intent(in) :: freqarr
         real(kind=8), dimension(:,:), intent(in) :: fnu
         integer, intent(in) :: nf
         SELECT CASE (e%type)
           CASE (EHYBRIDTHPL)
             e%ncgs=ncgs; e%bcgs=bcgs; e%tcgs=tcgs; e%ncgsnth=ncgsnth 
           CASE (EHYBRIDTH)
             e%ncgs=ncgs; e%bcgs=bcgs; e%tcgs=tcgs; e%ncgsnth=ncgsnth 
           CASE (EHYBRIDPL)
             e%ncgs=ncgs; e%bcgs=bcgs; e%tcgs=tcgs; e%ncgsnth=ncgsnth 
           CASE (EMAXJUTT) !alwinnote 2015/03/05
             e%ncgs=ncgs; e%bcgs=bcgs; e%tcgs=tcgs
           CASE (EMAXCOMP) !alwinnote 2015/03/05
             e%ncgs=ncgs; e%bcgs=bcgs; e%tcgs=tcgs
           CASE (EPOLSYNCHTH)
             e%ncgs=ncgs; e%bcgs=bcgs; e%tcgs=tcgs
           CASE (EPOLSYNCHPL)
             e%ncgs=ncgs; e%bcgs=bcgs; e%tcgs=tcgs; e%ncgsnth=ncgsnth
           CASE (ESYNCHPL)
             e%ncgs=ncgs; e%bcgs=bcgs; e%tcgs=tcgs; e%ncgsnth=ncgsnth
           CASE (ESYNCHTHAV)
             e%ncgs=ncgs; e%bcgs=bcgs; e%tcgs=tcgs
           CASE (ESYNCHTHAVNOABS)
             e%ncgs=ncgs; e%bcgs=bcgs; e%tcgs=tcgs
           CASE (EBB)
             e%tcgs=tcgs
           CASE (EFBB)
             e%tcgs=tcgs
           CASE (ERHO)
             e%ncgs=ncgs
           CASE (EBBPOL)
             e%tcgs=tcgs
           CASE (EINTERP)
!             write(6,*) 'assign emis: ',size(fnu),size(e%fnu)
!             write(6,*) 'assign emis: ',size(e%freqarr),size(freqarr)
             e%fnu=reshape(transpose(fnu),(/size(fnu)/)); e%freqarr=freqarr
         END SELECT
         end subroutine assign_emis_params

         subroutine initialize_emis_params(ep,n)
           integer, intent(in) :: n
           type (emis_params), intent(inout) :: ep
           integer :: nvals,k,i
           allocate(ep%gmin(n))
           allocate(ep%mu(n))
           nvals=sum(ep%coefindx)
           if(nvals.gt.0) then
              allocate(ep%usecoefs(nvals))
              k=1
              do i=1,size(ep%coefindx)
                 if(ep%coefindx(i).eq.1) then
                    ep%usecoefs(k)=i
                    k=k+1
                 endif
              enddo
           endif
!           write(6,*) 'emis coefindx: ',ep%coefindx,ep%usecoefs
!           allocate(ep%args(n))
         end subroutine initialize_emis_params

         subroutine del_emis_params(ep)
           type (emis_params), intent(inout) :: ep
           deallocate(ep%gmin)
           deallocate(ep%mu)
           if(allocated(ep%usecoefs)) deallocate(ep%usecoefs)
!           deallocate(ep%args)
         end subroutine del_emis_params

         subroutine polsynchemis_wrapper(nu,e)
         type (emis), intent(inout) :: e
         real(kind=8), intent(in), dimension(:) :: nu
         real(kind=8), dimension(e%npts,11) :: K, Kth, Kpl
!         write(6,*) 'polsynch', size(K), size(e%bcgs), size(nu)
   !      write(6,*) size(e%ncgs), size(e%tcgs), size(e%incang),e%npts
         select case(e%type)
           case(ehybridthpl)
             call polsynchth(nu,e%ncgs,e%bcgs,e%tcgs,e%incang,Kth)
             call polsynchpl(nu,e%ncgsnth,e%bcgs,e%incang,e%p,e%gmin, &
              e%gmax,Kpl)
             K = Kth + Kpl
           case(ehybridth)
             call polsynchth(nu,e%ncgs,e%bcgs,e%tcgs,e%incang,Kth)
             call polsynchpl(nu,e%ncgsnth,e%bcgs,e%incang,e%p,e%gmin, &
              e%gmax,Kpl)
             K = Kth
             K(:,5:11) = Kth(:,5:11) + Kpl(:,5:11)

           case(ehybridpl)
             call polsynchth(nu,e%ncgs,e%bcgs,e%tcgs,e%incang,Kth)
             call polsynchpl(nu,e%ncgsnth,e%bcgs,e%incang,e%p,e%gmin, &
              e%gmax,Kpl)
             K = Kpl
             K(:,5:11) = Kth(:,5:11) + Kpl(:,5:11)

           case(epolsynchth)
              write(6,*) 'polsynchth nu: ',nu
              write(6,*) 'polsynchth ncgs: ',e%ncgs
              write(6,*) 'polsynchth bcgs: ',e%bcgs
              write(6,*) 'polsynchth tcgs: ',e%tcgs
             call polsynchth(nu,e%ncgs,e%bcgs,e%tcgs,e%incang,K)
           case(epolsynchpl)
       !      write(6,*) 'ppl', e%ncgsnth, e%bcgs, e%incang, e%p
             call polsynchpl(nu,e%ncgsnth,e%bcgs,e%incang,e%p,e%gmin, &
              e%gmax,K)
           case(esynchpl)
             call synchpl(nu,e%ncgsnth,e%bcgs,e%incang,e%p,e%gmin, &
              e%gmax,K)
         end select
         write(6,*) 'afterpolsynch'
         if (e%neq==4) then
           e%j=K(:,1:4)
           e%K=K(:,5:11)
!           write(6,*) 'pol assign e ',e%K
         else
           e%j(:,1)=K(:,1)
           e%K(:,1)=K(:,5)
!           write(6,*) 'assign e', e%K
         endif
         end subroutine polsynchemis_wrapper

         subroutine lambda(e)
         type (emis), intent(inout) :: e
         e%j=1d0; e%K=0d0
         end subroutine lambda

         subroutine del_emissivity(e)
         type (emis), intent(inout) :: e
         deallocate(e%j); deallocate(e%K); deallocate(e%rshift)
         deallocate(e%cosne); deallocate(e%gmin)
         SELECT CASE (e%type)
           CASE (EHYBRIDTHPL)
             deallocate(e%ncgs); deallocate(e%tcgs)
             deallocate(e%bcgs); deallocate(e%incang)
             deallocate(e%ncgsnth); deallocate(e%p)
             call del_polsynchpl(e%neq)              
           CASE (EHYBRIDTH)
             deallocate(e%ncgs); deallocate(e%tcgs)
             deallocate(e%bcgs); deallocate(e%incang)
             deallocate(e%ncgsnth); deallocate(e%p)
             call del_polsynchpl(e%neq)              
           CASE (EHYBRIDPL)
             deallocate(e%ncgs); deallocate(e%tcgs)
             deallocate(e%bcgs); deallocate(e%incang)
             deallocate(e%ncgsnth); deallocate(e%p)
             call del_polsynchpl(e%neq)              
           CASE (EMAXJUTT) !alwinnote 2015/03/05
             deallocate(e%tcgs); deallocate(e%ncgs)
             deallocate(e%bcgs); deallocate(e%incang)
           CASE (EMAXCOMP) !alwinnote 2015/07/22
             deallocate(e%tcgs); deallocate(e%ncgs)
             deallocate(e%bcgs); deallocate(e%incang)
!             deallocate(e%args)
           CASE (EPOLSYNCHTH)
             deallocate(e%tcgs); deallocate(e%ncgs)
             deallocate(e%bcgs); deallocate(e%incang)
           CASE (EPOLSYNCHPL)
             deallocate(e%ncgs); deallocate(e%tcgs)
             deallocate(e%bcgs); deallocate(e%incang)
             deallocate(e%ncgsnth); deallocate(e%p)
             call del_polsynchpl(e%neq)
           CASE (ESYNCHPL)
             deallocate(e%ncgs); deallocate(e%tcgs)
             deallocate(e%bcgs); deallocate(e%incang)
             deallocate(e%ncgsnth); deallocate(e%p)
             call del_polsynchpl(e%neq)
           CASE (ESYNCHTHAV)
             deallocate(e%tcgs); deallocate(e%ncgs)
             deallocate(e%bcgs); deallocate(e%incang)
           CASE (ESYNCHTHAVNOABS)
             deallocate(e%tcgs); deallocate(e%ncgs)
             deallocate(e%bcgs); deallocate(e%incang)
           CASE (EBB)
             deallocate(e%tcgs); deallocate(e%incang)
           CASE (EFBB)
             deallocate(e%tcgs); deallocate(e%incang)
           CASE (ERHO)
             deallocate(e%ncgs); deallocate(e%incang)
           CASE (EBBPOL)
             deallocate(e%tcgs); deallocate(e%incang)
           CASE (EINTERP)
              deallocate(e%fnu)
         END SELECT
         e%npts=0
         end subroutine del_emissivity

         subroutine calc_opt_depth(s,e,tau,indx)
         real(kind=8), intent(in), dimension(:) :: s
         type (emis), intent(in) :: e
         real(kind=8), intent(out), dimension(size(s)) :: tau
         integer, intent(in) :: indx
!         write(6,*) 'opt depth: ',size(s),size(e%K(:,1))
!         write(6,*) 's: ', s
!         write(6,*) 'K: ',e%K(:,1)
         tau=tsum(s,abs(e%K(:,indx)))
         end subroutine calc_opt_depth

         subroutine rotate_emis(e,s2xi,c2xi)
         type (emis), intent(inout) :: e
         real(kind=8), dimension(:), intent(in) :: s2xi,c2xi
         real(kind=8), dimension(size(s2xi)) :: ju,jq,rhoq,rhou,aq,au
! Based on Shcherbakov & Huang (2011) and Mathematica rotate_emis.nb
         jq=e%j(:,2); ju=e%j(:,3); aq=e%K(:,2); au=e%K(:,3)
         rhoq=e%K(:,5); rhou=e%K(:,6)
!         write(6,*) 'c2xi: ',c2xi
!         write(6,*) 's2xi: ',s2xi
! if bcgs = 0, get NaN s2xi,c2xi. set to arbitrary angle.
!         where(e%bcgs.eq.0d0)
!            c2xi=1d0
!            s2xi=0d0
!         endwhere
         e%j(:,2)=c2xi*jq-s2xi*ju
         e%j(:,3)=s2xi*jq+c2xi*ju
         e%K(:,2)=c2xi*aq-s2xi*au
         e%K(:,3)=s2xi*aq+c2xi*au
         e%K(:,5)=c2xi*rhoq-s2xi*rhou
         e%K(:,6)=s2xi*rhoq+c2xi*rhou
!         e%j(:,3)=s2xi*e%j(:,2); e%j(:,2)=c2xi*e%j(:,2)
!         e%K(:,3)=s2xi*e%K(:,2)+c2xi*e%K(:,3); e%K(:,2)=c2xi*e%K(:,2)-s2xi*e%K(:,3)
!         e%K(:,6)=s2xi*e%K(:,5); e%K(:,5)=c2xi*e%K(:,5)
         if(any(isnan(e%j).or.any(isnan(e%K)))) then
            write(6,*) 'NaN in rotate_emis!'
!            write(6,*) 's2xi: ',s2xi
!            write(6,*) 'c2xi: ',c2xi
!            write(6,*) 'jq: ',jq
!            write(6,*) 'ej2: ',e%j(:,2)
!            write(6,*) 'ang: ',e%incang
!            write(6,*) 'ek5: ',e%K(:,5)
         endif
         end subroutine rotate_emis

         subroutine invariant_emis(e,g)
         type (emis), intent(inout) :: e
         integer :: i
         real(kind=8), dimension(:), intent(in) :: g
!         e%rshift=g
         do i=1,e%neq; e%j(:,i)=e%j(:,i)*g*g; enddo
         do i=1,e%nk; e%K(:,i)=e%K(:,i)/g; enddo
         end subroutine invariant_emis

         subroutine invariant_intensity(e,g,npow)
         type (emis), intent(inout) :: e
         real(kind=8), dimension(:), intent(in) :: g
         integer :: i
         integer, intent(in) :: npow
!         e%rshift=g
         do i=1,e%neq; e%j(:,i)=e%j(:,i)*g**npow; enddo
         end subroutine invariant_intensity

         subroutine emis_model(e,ep)
         type (emis), intent(inout) :: e
         type (emis_params), intent(in) :: ep
         SELECT CASE (e%type)
           CASE (EHYBRIDTHPL)
             call emis_model_synchth(e,ep%mu)
             call emis_model_synchpl(e,ep%gmin,ep%gmax,ep%p1)
           CASE (EHYBRIDTH) !alwinnote
             call emis_model_synchth(e,ep%mu)
             call emis_model_synchpl(e,ep%gmin,ep%gmax,ep%p1)
           CASE (EHYBRIDPL) !alwinnote
             call emis_model_synchth(e,ep%mu)
             call emis_model_synchpl(e,ep%gmin,ep%gmax,ep%p1)
           CASE (EMAXJUTT) !alwinnote 2015/03/05
             call emis_model_synchth(e,ep%mu)
           CASE (EMAXCOMP) !alwinnote 2015/07/22
             call emis_model_synchth(e,ep%mu)
           CASE (EPOLSYNCHTH)
             call emis_model_synchth(e,ep%mu)
           CASE (EPOLSYNCHPL)
             call emis_model_synchpl(e,ep%gmin,ep%gmax,ep%p1)
           CASE (ESYNCHPL)
             call emis_model_synchpl(e,ep%gmin,ep%gmax,ep%p1)
           CASE (EINTERP)
             call emis_model_interp(e,ep%nfreq_tab,ep%freq_tab)
         END SELECT
         end subroutine emis_model

         subroutine emis_model_synchpl(e,gmin,gmax,p)
         type (emis), intent(inout) :: e
         real, intent(in) :: gmax,p
         real, intent(in), dimension(:) :: gmin
 !        write(6,*) 'em: ',gmin,gmax,p
         e%p=p; e%gmin=gmin; e%gmax=gmax
!         e%ncgsnth=e%ncgs
         end subroutine emis_model_synchpl

         subroutine emis_model_interp(e,nfreq,freq)
         type (emis), intent(inout) :: e
         integer, intent(in) :: nfreq
         real, dimension(:), intent(in) :: freq
         end subroutine emis_model_interp

         subroutine emis_model_synchth(e,mu)
         type (emis), intent(inout) :: e
         real, dimension(:), intent(in) :: mu
! This is constant Ti/Te electron model
! UPDATE 3/13/2017 don't want this here since e- temp can depend on other parameters
! for now only uses things in e, but could e.g. use polar angle or who knows what else.
!         e%tcgs=e%tcgs*mu
         end subroutine emis_model_synchth

         subroutine emis_model_lambda(e)
         type (emis), intent(inout) :: e
         end subroutine emis_model_lambda 

       end module emissivity
