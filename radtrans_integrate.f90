  module radtrans_integrate

    use odepack, only: lsoda_basic
    use interpolate, only: get_weight, locate
    use math, only: tsum

    implicit none

    ! global data to avoid needing to pass arguments through external integration routines (e.g., lsoda)
    integer :: lindx=25,nequations,nptstot,npts,nptsout,iflag
    real(kind=8), dimension(:,:), allocatable :: KK,jj,intensity,PP
    real(kind=8), dimension(:), allocatable :: s,ss,s0,tau,stokesq,stokesu,stokesv
    real(kind=8), dimension(:,:,:), allocatable :: QQ,imm,OO

    integer :: IS_LINEAR_STOKES = 1
    integer(kind=4) :: maxsteps = 100000
    integer, dimension(4) :: stats
    real(kind=8) :: ortol = 1d-6, oatol = 1d-8, hmax = 10, MAX_TAU, MAX_TAU_DEFAULT = 10d0, thin = 1d-2, sphtolfac=1d0

!$omp threadprivate(ss,tau,s0,jj,KK,intensity,lindx,nptstot,npts,nptsout,s,stats,stokesq,stokesu,stokesv)
!$omp threadprivate(ortol,oatol,hmax,thin,MAX_TAU,QQ,PP,imm,OO,IS_LINEAR_STOKES,sphtolfac)

    interface integrate
       module procedure integrate
    end interface

    interface calc_delo_Q
       module procedure calc_delo_Q_single
    end interface
    
    interface calc_delo_P
       module procedure calc_delo_P_single
    end interface

    interface imatrix_4
       module procedure imatrix_4_single
    end interface

    interface invert_delo_matrix
       module procedure invert_delo_matrix_single
    end interface

    interface calc_delo_Q_thin
       module procedure calc_delo_Q_thin_single
    end interface
    
    interface calc_delo_P_thin
       module procedure calc_delo_P_thin_single
    end interface

    interface invert_delo_matrix_thin
       module procedure invert_delo_matrix_thin_single
    end interface

    interface calc_O
       module procedure calc_O
    end interface

    interface radtrans_integrate_formal
       module procedure radtrans_integrate_formal
    end interface

    contains

      subroutine init_radtrans_integrate_data(riflag,rneq,gnpts,rnpts,maxtau,hm,oa,ort,th,maxst)
        integer, intent(in) :: riflag,rneq,gnpts,rnpts
!        real(kind=8), intent(in), optional :: maxtau,hm,oa,ort
        real(kind=8), intent(in), optional :: maxtau,hm,oa,ort,th
        integer(kind=4), intent(in), optional :: maxst
        nequations = rneq; npts = gnpts; iflag = riflag; nptsout = rnpts
        if(present(maxtau)) then
!        MAX_TAU = maxtau; hmax = hm; oatol = oa; ortol = ort
           MAX_TAU = maxtau
        else
           MAX_TAU=MAX_TAU_DEFAULT
        endif
        if(present(hm)) then
           hmax = hm
        else
           hmax = 0.1d0
        endif
        if(present(ort)) then 
           ortol = ort
        else
           ortol = 1d-6
        endif
        if(present(oa)) then
           oatol = oa
        else
           oatol = 1d-8
        endif
        if(present(th)) then
           thin = th
        else
           thin = 1d-2
        endif
        if(present(maxst)) then
           maxsteps = maxst
        else
           maxsteps = 100000
        endif
        allocate(tau(npts))
        allocate(s0(npts))
        allocate(jj(npts,nequations))
        allocate(KK(npts,1+(nequations)*(nequations-1)/2))
        allocate(intensity(nequations,npts))
        allocate(s(npts)); allocate(ss(npts))
        s(:) = 0d0
        ss(:) = 0d0
        intensity(:,:) = 0d0
        jj(:,:) = 0d0
        tau(:) = 0d0
        KK(:,:) = 0d0
        s0(:) = 0d0
        if(iflag==1) then
           allocate(QQ(npts,nequations,nequations)); allocate(PP(npts,nequations))
           allocate(imm(npts,nequations,nequations))
           QQ = 0d0; PP = 0d0; imm = 0d0
        elseif(iflag==2) then
           allocate(OO(nequations,nequations,npts))
           OO = 0d0
        elseif(iflag==3) then
           allocate(stokesq(npts))
           allocate(stokesu(npts))
           allocate(stokesv(npts))
        endif
      end subroutine init_radtrans_integrate_data

      subroutine integrate(sinput,ej,eK,tauin,rnpts)
! these array sizes are really all known in terms of npts, nequations...
        integer, intent(inout) :: rnpts
!        real(kind=8), dimension(:,:), intent(inout) :: rI
        real(kind=8), dimension(:,:), intent(in) :: ej,eK
        real(kind=8), dimension(:), intent(in) :: sinput
        real(kind=8), dimension(:), intent(in) :: tauin
!        real(kind=8), dimension(:), intent(inout) :: I0
!        real(kind=8), intent(in), optional :: maxt
    ! Subroutine to call integration routine for RT equation
    ! JAD 2/28/2011 moved from grtrans_driver 8/7/2014
        s0=sinput; s=sinput; ss=sinput; intensity=0d0
        tau = tauin
        jj=ej; KK=eK
        if(any(isnan(jj))) then
           write(6,*) 'nan in radtrans_integrate j'
!           write(6,*) 'radtrans_integrate j1: ',jj(:,1)
!           write(6,*) 'radtrans_integrate j2: ',jj(:,2)
        endif
        if(any(isnan(KK))) then
           write(6,*) 'nan in radtrans_integrate K'
!           write(6,*) 'radtrans_integrate K1: ',KK(:,1)
!           write(6,*) 'radtrans_integrate K2: ',KK(:,2)
!           write(6,*) 'radtrans_integrate K4: ',KK(:,4)
!           write(6,*) 'radtrans_integrate K5: ',KK(:,5)
!           write(6,*) 'radtrans_integrate K7: ',KK(:,7)
        endif
! spherical stokes case
        if (iflag==0) then
           IS_LINEAR_STOKES=1
           
           call radtrans_integrate_lsoda()
        elseif (iflag==3) then
           IS_LINEAR_STOKES=0
           call radtrans_integrate_lsoda()
        elseif (iflag==1) then
           if(nequations==4) then
              call radtrans_integrate_delo(s,tau,jj,KK(:,1:4),KK(:,5:7),PP,QQ,imm)
           else
              call radtrans_integrate_quadrature(s,jj(:,1),KK(:,1))
           endif
        elseif (iflag==2) then
           if(nequations==4) then
              call radtrans_integrate_formal(s,jj,KK(:,1:4),KK(:,5:7),OO)
           else
              call radtrans_integrate_quadrature(s,jj(:,1),KK(:,1))
           endif
        endif
!        write(6,*) 'assign intensity', size(rI,1), size(rI,2), size(intensity,1), size(intensity,2), nptsout
!        rI(1:nequations,1:nptsout) = intensity(1:nequations,1:nptsout);
        rnpts = nptsout
!        write(6,*) 'intensity: ', rnpts, rI(rnpts,1)
        if(isnan(intensity(1,rnpts))) then
           write(6,*) 'NaN in integrate ej: '!,ej
!           write(6,*) 'NaN in integrate jj 2: ',jj
!           write(6,*) 'NaN in integrate eK: ',eK
        endif
        return
      end subroutine integrate
  
      subroutine radtrans_integrate_lsoda()
        real(kind=8), dimension(4) :: I0
        real(kind=8), dimension(npts) :: dummy
        real(kind=8), dimension(:,:), allocatable :: tau_arr
        real(kind=8), dimension(:), allocatable :: tau_temp,intvals,q,u,v
        integer, dimension(:), allocatable :: inds
        real(kind=8) :: weight
        integer :: lamdex,i,ii,i1,i2,taudex
        i1=1; i2=nptstot
        I0 = 0d0
        if(maxval(tau).le.MAX_TAU) then
           lamdex=npts
        else
           call locate(tau,MAX_TAU,lamdex)
           lamdex=lamdex+1
         !           write(6,*) 'locate', lamdex
!           if(lamdex==2) then
!              call get_weight(tau,MAX_TAU,lamdex,weight)
! need to re-work this. just trying to interpolate between s0(1) and s0(2).
!             write(6,*) 'weight: ',lamdex,weight,s(lamdex),s(lamdex+1)
!              s=(/(s(1)-s(2))*(1d0-weight),0d0/)
!              lamdex=lamdex+1
!             write(6,*) 'radtrans integrate s: ',s(1:lamdex),lamdex
!           endif
        endif
! Only use parts of ray where emissivity is non-zero:
        if (jj(1,1).eq.0d0) then
           ! Find first non-zero element:
           do ii=2,npts
              if (jj(ii,1).ne.0d0) then
                 i1=ii
                 exit
              endif
           enddo
        else
           i1=1
        endif
      !         i1=1
        if (jj(lamdex,1).eq.0d0) then
! Find last non-zero element:
           do ii=1,lamdex-1
              if (jj(lamdex-ii,1).ne.0d0) then
                 i2=lamdex-ii
                 exit
              endif
           enddo
        else
           i2=lamdex
        endif
!         write(6,*) 'i1i2: ',i1,i2,npts, jj(i2,1), jj(i2-1,1)
        nptsout=i2
!         write(6,*) 's: ',s(lamdex),neq,npts,tau(lamdex),lamdex
   !      write(6,*) 'I0: ',s(npts), tau(npts)
 !        write(6,*) 'lam: ',s0
! try to figure out hmax:
        ss(1:npts)=s0(npts:1:-1)
        s(i1:i2)=s(i2:i1:-1)
!        write(6,*) 'integrate s: ',minval(s), maxval(s), i1, i2, lamdex
! should make this the start of grtrans_integrate_lsoda and put everything else in general integrate subroutine since it's generic to all methods
        if(nequations==4) then
           if(IS_LINEAR_STOKES==1) then
              call lsoda_basic(radtrans_lsoda_calc_rhs,I0(1:nequations), &
                   s(i1:i2),oatol, &
                   ortol,radtrans_lsoda_calc_jac,intensity(:,i1:i2), &
                   1,maxsteps,stats,hmax=hmax)
           else
! spherical stokes under development
!              I0=I0+I0sphmin
!              I0(1)=1d-64; I0(2)=1d-64
              I0(1)=1d-8; I0(2)=1d-8
              I0(3)=0.1; I0(4)=-0.1
!              I0=0d0
! testing higher error tolerance for sph stokes to converge on difficult problems
              call lsoda_basic(radtrans_lsoda_calc_rhs_sph,I0(1:nequations), &
                   s(i1:i2),oatol*sphtolfac, &
                   ortol*sphtolfac,radtrans_lsoda_calc_jac_sph,intensity(:,i1:i2), &
                   1,maxsteps,stats,hmax=hmax)
! convert to linear stokes parameters
!              allocate(q(npts));allocate(u(npts))
!              allocate(v(npts))
!              write(6,*) 'sph stokes intensity: ',maxval(intensity(1,i1:i2)), &
!                   maxval(intensity(2,i1:i2))
              stokesq(i1:i2)=intensity(2,i1:i2)*sin(intensity(4,i1:i2)) &
                   *cos(intensity(3,i1:i2))
              stokesu(i1:i2)=intensity(2,i1:i2)*sin(intensity(4,i1:i2)) &
                   *sin(intensity(3,i1:i2))
              stokesv(i1:i2)=intensity(2,i1:i2)*cos(intensity(4,i1:i2))
              intensity(2,i1:i2)=stokesq(i1:i2); intensity(3,i1:i2)=stokesu(i1:i2)
              intensity(4,i1:i2)=stokesv(i1:i2)
!              write(6,*) 'sph stokes intensity: ',maxval(q),maxval(v)
!              write(6,*) 'sph stokes intensity: ',maxval(intensity(2,i1:i2))
!              write(6,*) 'sph stokes intensity: ',maxval(intensity(1,i1:i2))
!              deallocate(q); deallocate(u); deallocate(v)
           endif
        else
           call lsoda_basic(radtrans_lsoda_calc_rhs_npol, &
              I0(1:nequations),s(i1:i2),oatol, &
              ortol,radtrans_lsoda_calc_jac_npol,intensity(:,i1:i2), &
              1,maxsteps,stats,hmax=hmax)
        endif
        if(isnan(intensity(1,i2))) then
!           write(6,*) 'NaN in integrate: ',i1,i2,s(i1:i2)
!           write(6,*) 'NaN in integrate intensity: ',intensity(1,i1:i2)
!           write(6,*) 'NaN in integrate j: ',jj(i1:i2,:)
!           write(6,*) 'NaN in integrate K: ',KK(i1:i2,:)
        endif
        return
      end subroutine radtrans_integrate_lsoda

!      subroutine convert_sph_lin_stokes()
! convert spherical stokes output to normal (I,Q,U,V)
!        intensity(2,i1:i2

      subroutine radtrans_lsoda_calc_rhs(neq,lam,I,dIdlam)
! Compute RHS dIdlam for LSODA
        integer, intent(in) :: neq
        real(kind=8), intent(in) :: lam
        real(kind=8), intent(in), dimension(neq) :: I
        real(kind=8), intent(out), dimension(neq) :: dIdlam
!         real(kind=8), intent(out), dimension(neq,neq) :: jac
        real(kind=8), dimension(neq) :: j
        real(kind=8), dimension(1+neq*(neq-1)/2) :: K
        call radtrans_aux(neq,lam,j,K)
        call radtrans_rhs_form(neq,j,K,dIdlam,I)
!         write(6,*) 'dIdlam: ',lam,dIdlam
!         write(6,*) 'jk: ',jj(1),jj(2),j(4)
!         write(6,*) 'K: ',K(1),K(4),
!     &    K(5),K(7)
        return 
      end subroutine radtrans_lsoda_calc_rhs

      subroutine radtrans_lsoda_calc_rhs_npol(neq,lam,I,dIdlam)
! Compute RHS dIdlam for LSODA
        integer, intent(in) :: neq
        real(kind=8), intent(in) :: lam
        real(kind=8), intent(in), dimension(neq) :: I
        real(kind=8), intent(out), dimension(neq) :: dIdlam
!         real(kind=8), intent(out), dimension(neq,neq) :: jac
        real(kind=8), dimension(neq) :: j
        real(kind=8), dimension(1+neq*(neq-1)/2) :: K
        call radtrans_aux(neq,lam,j,K)
        call radtrans_rhs_form_npol(neq,j,K,dIdlam,I)
!         didlam=didlam*1e25
!         write(6,*) 'dIdlam: ',lam,dIdlam
!         write(6,*) 'jk: ',j(1),j(2),j(4)
!         write(6,*) 'K: ',K(1),K(4),
!     &    K(5),K(7)
      end subroutine radtrans_lsoda_calc_rhs_npol
         
      subroutine radtrans_rhs_form(neq,j,K,dIdlam,I)
        integer, intent(in) :: neq
        real(kind=8), intent(in), dimension(neq) :: j
        real(kind=8), intent(in), dimension(1+neq*(neq-1)/2) :: K
        real(kind=8), intent(out), dimension(neq) :: dIdlam
        real(kind=8), intent(in), dimension(neq) :: I
      !         write(6,*) 'rhs: ',IS_LINEAR_STOKES,size(I),size(K),size(J)
!        if (IS_LINEAR_STOKES==1) then
           dIdlam(1)=j(1)-(K(1)*I(1)+K(2)*I(2)+K(3)*I(3)+K(4)*I(4))
           dIdlam(2)=j(2)-(K(2)*I(1)+K(1)*I(2)+K(7)*I(3)-K(6)*I(4))
           dIdlam(3)=j(3)-(K(3)*I(1)-K(7)*I(2)+K(1)*I(3)+K(5)*I(4))
           dIdlam(4)=j(4)-(K(4)*I(1)+K(6)*I(2)-K(5)*I(3)+K(1)*I(4))
!        endif
      end subroutine radtrans_rhs_form

      subroutine radtrans_rhs_form_npol(neq,j,K,dIdlam,I)
        integer, intent(in) :: neq
        real(kind=8), intent(in), dimension(neq) :: j
        real(kind=8), intent(in), dimension(1+neq*(neq-1)/2) :: K
        real(kind=8), intent(out), dimension(neq) :: dIdlam
        real(kind=8), intent(in), dimension(neq) :: I
!         write(6,*) 'rhs npol: ',size(I),size(K),size(J)
!         dIdlam(1)=maxval((/j(1)-K(1)*I(1),0d0/))
        dIdlam(1)=j(1)-K(1)*I(1)
!         write(6,*) 'rhs: ',dIdlam(1),j(1),K(1),I(1)
        return
      end subroutine radtrans_rhs_form_npol

      subroutine radtrans_jac_form_npol(neq,j,K,nrowpd,pd)
        integer, intent(in) :: neq, nrowpd
        real(kind=8), intent(in), dimension(neq) :: j
        real(kind=8), intent(in), dimension(1+neq*(neq-1)/2) :: K
        real(kind=8), intent(out), dimension(nrowpd,neq) :: pd
   !      write(6,*) 'jac: ',nrowpd,neq,size(K),K(1)
        pd(1,1)=-K(1)
!         write(6,*) 'pd: ',pd
        return
      end subroutine radtrans_jac_form_npol

      subroutine radtrans_jac_form(neq,j,K,nrowpd,pd)
        integer, intent(in) :: neq, nrowpd
        real(kind=8), intent(in), dimension(neq) :: j
        real(kind=8), intent(in), dimension(1+neq*(neq-1)/2) :: K
        real(kind=8), intent(out), dimension(nrowpd,neq) :: pd
        !      write(6,*) 'jac: ',nrowpd,neq,size(K)
!        if (IS_LINEAR_STOKES==1) then
           pd(1,1)=K(1)
           pd(1,2)=K(2)
           pd(1,3)=K(3)
           pd(1,4)=K(4)
           pd(2,1)=K(2)
           pd(2,2)=K(1)
           pd(2,3)=K(7)
           pd(2,4)=-K(6)
           pd(3,1)=K(3)
           pd(3,2)=-K(7)
           pd(3,3)=K(1)
           pd(3,4)=K(5)
           pd(4,1)=K(4)
           pd(4,2)=K(6)
           pd(4,3)=-K(5)
           pd(4,4)=K(1)
           pd=-1d0*pd
!        endif
!         write(6,*) 'pd: ',pd
        return
      end subroutine radtrans_jac_form
       
      subroutine radtrans_lsoda_calc_jac(neq,lam,I,ml &
           ,mu,pd,nrowpd)
! Compute Jacobian for LSODA
        integer, intent(in) :: neq, nrowpd
        real(kind=8), intent(in) :: lam
        real(kind=8), intent(in), dimension(neq) :: I
        real(kind=8), intent(in) :: ml
        real(kind=8), intent(in) :: mu
        real(kind=8), intent(out), dimension(nrowpd,neq) :: pd
        real(kind=8), dimension(neq) :: j
        real(kind=8), dimension(1+neq*(neq-1)/2) :: K
!         write(6,*) 'jac: ',nrowpd
        call radtrans_aux(neq,lam,j,K)
        call radtrans_jac_form(neq,j,K,nrowpd,pd)
        !         write(6,*) 'pd: ', pd
        return
      end subroutine radtrans_lsoda_calc_jac

      subroutine radtrans_lsoda_calc_jac_npol(neq,lam,I,ml &
           ,mu,pd,nrowpd)
! Compute Jacobian for LSODA
        integer, intent(in) :: neq, nrowpd
        real(kind=8), intent(in) :: lam
        real(kind=8), intent(in), dimension(neq) :: I
        real(kind=8), intent(in) :: ml
        real(kind=8), intent(in) :: mu
        real(kind=8), intent(out), dimension(nrowpd,neq) :: pd
        real(kind=8), dimension(neq) :: j
        real(kind=8), dimension(1+neq*(neq-1)/2) :: K
!         write(6,*) 'jac: ',nrowpd
        call radtrans_aux(neq,lam,j,K)
        call radtrans_jac_form_npol(neq,j,K,nrowpd,pd)
!         write(6,*) 'pd: ', pd
        return
      end subroutine radtrans_lsoda_calc_jac_npol

      subroutine radtrans_aux(neq,lam,j,K)
        integer, intent(in) :: neq
        real(kind=8), intent(in) :: lam
        real(kind=8), intent(out), dimension(neq) :: j
        real(kind=8), intent(out), dimension(1+neq*(neq-1)/2) :: K
        real(kind=8) :: weight
        integer :: indx,uindx
! here is where we would need to change things: have a new routine that calculates new j, K with new rel. factors at point lam
        call get_weight(ss,lam,lindx,weight)
        indx=npts-lindx+1; uindx=minval((/indx+1,npts/))
        j=(1d0-weight)*jj(indx,:)+weight*jj(uindx,:)
        K=(1d0-weight)*KK(indx,:)+weight*KK(uindx,:)
      end subroutine radtrans_aux

! spherical stokes stuff in development
      subroutine radtrans_lsoda_calc_rhs_sph(neq,lam,I,dIdlam)
! Compute RHS dIdlam for LSODA
        integer, intent(in) :: neq
        real(kind=8), intent(in) :: lam
        real(kind=8), intent(in), dimension(neq) :: I
        real(kind=8), intent(out), dimension(neq) :: dIdlam
!         real(kind=8), intent(out), dimension(neq,neq) :: jac
        real(kind=8), dimension(neq) :: j
        real(kind=8), dimension(1+neq*(neq-1)/2) :: K
        call radtrans_aux(neq,lam,j,K)
        call radtrans_rhs_form_sph(neq,j,K,dIdlam,I)
!         write(6,*) 'dIdlam: ',lam,dIdlam
!         write(6,*) 'jk: ',jj(1),jj(2),j(4)
!         write(6,*) 'K: ',K(1),K(4),
!     &    K(5),K(7)
        return 
      end subroutine radtrans_lsoda_calc_rhs_sph

      subroutine radtrans_rhs_form_sph(neq,j,K,dIdlam,I)
        integer, intent(in) :: neq
        real(kind=8), intent(in), dimension(neq) :: j
        real(kind=8), intent(in), dimension(1+neq*(neq-1)/2) :: K
        real(kind=8), intent(out), dimension(neq) :: dIdlam
        real(kind=8), intent(in), dimension(neq) :: I
        real(kind=8) :: sphi,spsi,cphi,cpsi
      !         write(6,*) 'rhs: ',IS_LINEAR_STOKES,size(I),size(K),size(J)
!        if (IS_LINEAR_STOKES==1) then
!           dIdlam(1)=j(1)-(K(1)*I(1)+K(2)*I(2)+K(3)*I(3)+K(4)*I(4))
!           dIdlam(2)=j(2)-(K(2)*I(1)+K(1)*I(2)+K(7)*I(3)-K(6)*I(4))
!           dIdlam(3)=j(3)-(K(3)*I(1)-K(7)*I(2)+K(1)*I(3)+K(5)*I(4))
!           dIdlam(4)=j(4)-(K(4)*I(1)+K(6)*I(2)-K(5)*I(3)+K(1)*I(4))
        sphi=sin(I(3)); cphi=cos(I(3)); spsi=sin(I(4)); cpsi=cos(I(4))
        dIdlam(1)=j(1)-K(1)*I(1)-(cphi*spsi*K(2)+sphi*spsi*K(3)+cpsi*K(4))*I(2)
        dIdlam(2)=-K(1)*I(2)-cphi*spsi*K(2)*I(1)-sphi*spsi*K(3)*I(1)+ &
             spsi*(cphi*j(2)+sphi*j(3))+cpsi*(-I(1)*K(4)+j(4))
!        dIdlam(3)=1d0/I(2)*(-cphi/spsi*(I(1)*K(3)-j(3)+I(2)*cpsi*K(5))- &
!             sphi/spsi*(-I(1)*K(2)+j(2)+I(2)*cpsi*K(6))+I(2)*K(7))
        dIdlam(3)=1d0/I(2)/spsi*(cphi*(j(3)-I(1)*K(3))+sphi*(-j(2)+I(1)*K(2))) &
             -cpsi/spsi*(cphi*K(5)+sphi*K(6))+K(7)
!        dIdlam(4)=1d0/I(2)*(-I(1)*cphi*cpsi*K(2)+spsi*(-j(4)+I(1)*K(4))- &
!             sphi*(I(1)*cpsi*K(3)-cpsi*j(3)+I(2)*K(5))+cphi*cpsi*j(2))+cphi*K(6)
        dIdlam(4)=1d0/I(2)*(spsi*(-j(4)+I(1)*K(4))+cpsi*(cphi*(j(2)-K(2)*I(1)) &
             +sphi*(-I(1)*K(3)+j(3))))+cphi*K(6)-sphi*K(5)
!        endif
      end subroutine radtrans_rhs_form_sph

      subroutine radtrans_jac_form_sph(neq,j,K,I,nrowpd,pd)
        integer, intent(in) :: neq, nrowpd
        real(kind=8), intent(in), dimension(neq) :: j,I
        real(kind=8), intent(in), dimension(1+neq*(neq-1)/2) :: K
        real(kind=8), intent(out), dimension(nrowpd,neq) :: pd
        real(kind=8) :: cphi,sphi,cpsi,spsi
        !      write(6,*) 'jac: ',nrowpd,neq,size(K)
!        if (IS_LINEAR_STOKES==1) then
        cphi=cos(I(3)); sphi=sin(I(3)); cpsi=cos(I(4)); spsi=sin(I(4))
        pd(1,1)=-K(1)
        pd(1,2)=-cphi*spsi*K(2)-sphi*spsi*K(3)-cpsi*K(4)
        pd(1,3)=-I(2)*(-sphi*spsi*K(2)+cphi*spsi*K(3))
        pd(1,4)=-I(2)*(cphi*cpsi*K(2)+cpsi*sphi*K(3)-spsi*K(4))
        pd(2,1)=pd(1,2)
        pd(2,2)=-K(1)
        pd(2,3)=I(1)*spsi*(sphi*K(2)-cphi*K(3))+spsi*(cphi*j(3)-sphi*j(2))
        pd(2,4)=-I(1)*cpsi*(cphi*K(2)+sphi*K(3))+cpsi*(cphi*j(2)+sphi*j(3))-spsi*(j(4)-I(1)*K(4))
!        pd(3,1)=1d0/I(2)*(sphi/spsi*K(2)-cphi/spsi*K(3))
        pd(3,1)=1d0/spsi/I(2)*(sphi*K(2)-cphi*K(3))
!        pd(3,2)=(-cphi*cpsi/spsi*K(5)-cpsi/spsi*sphi*K(6)+K(7))/I(2)-(-cphi/spsi*(I(1)*K(3)-& 
!             j(3)+I(2)*cpsi*K(5))-1d0/spsi*sphi*(-I(1)*K(2)+j(2)+I(2)*cpsi*K(6))+I(2)*K(7))/I(2)**2.
        pd(3,2)=1d0/spsi/I(2)**2d0*(cphi*(I(1)*K(3)-j(3))+sphi*(j(2)-I(1)*K(2)))
!        pd(3,3)=(1d0/spsi*sphi*(I(1)*K(3)-j(3)+I(2)*cpsi*K(5))-cphi/spsi*(-I(1)*K(2)+j(2)+I(2)*cpsi*K(6)))/I(2)
        pd(3,3)=1d0/I(2)/spsi*(sphi*(I(1)*K(3)-j(3))+cphi*(I(1)*K(2)-j(2)))+1d0/spsi*(sphi*K(5)-cphi*K(6))
!        pd(3,4)=(I(2)*cphi*K(5)+cphi*cpsi/spsi/spsi*(I(1)*K(3)-j(3)+I(2)*cpsi*K(5))+I(2)*sphi*K(6)+ &
!             cpsi/spsi/spsi*sphi*(-I(1)*K(2)+j(2)+I(2)*cpsi*K(6)))/I(2)
        pd(3,4)=cpsi/spsi/spsi/I(2)*(cphi*(I(1)*K(3)-j(3))-sphi*(I(1)*K(2)-j(2)))+1d0/spsi/spsi*(cphi*K(5)+sphi*K(6))
        pd(4,1)=(-cphi*cpsi*K(2)-cpsi*sphi*K(3)+spsi*K(4))/I(2)        
!        pd(4,2)=(-sphi*K(5)+cphi*K(6))/I(2)-(-I(1)*cphi*cpsi*K(2)+spsi*(I(1)*K(4)-j(4))-sphi*(I(1)*cpsi*K(3) &
!             -cpsi*j(3)+I(2)*K(5))+cphi*(cpsi*j(2)+I(2)*K(6)))/I(2)**2d0
        pd(4,2)=1d0/I(2)**2d0*(cphi*cpsi*(I(1)*K(2)-j(2))-spsi*(I(1)*K(4)-j(4))+sphi*cpsi*(I(1)*K(3)-j(3)))
!        pd(4,3)=(I(1)*cpsi*sphi*K(2)-cphi*(I(1)*cpsi*K(3)-cpsi*j(3)+I(2)*K(5))-sphi*(cpsi*j(2)+I(2)*K(6)))/I(2)
        pd(4,3)=1d0/I(2)*(cpsi*sphi*(I(1)*K(2)-j(2))-cpsi*cphi*(I(1)*K(3)-j(3)))-cphi*K(5)-sphi*K(6)
!        pd(4,4)=(I(1)*cphi*spsi*K(2)-cphi*spsi*j(2)-sphi*(-I(1)*spsi*K(3)+spsi*j(3))+cpsi*(I(1)*K(4)-j(4)))/I(2)
        pd(4,4)=1d0/I(2)*(cphi*spsi*(I(1)*K(2)-j(2))+sphi*spsi*(I(1)*K(3)-j(3))+cpsi*(I(1)*K(4)-j(4)))
!           pd=-1d0*pd
!        endif
!         write(6,*) 'pd: ',pd
        return
      end subroutine radtrans_jac_form_sph

      subroutine radtrans_jac_form_sph_old(neq,j,K,I,nrowpd,pd)
        integer, intent(in) :: neq, nrowpd
        real(kind=8), intent(in), dimension(neq) :: j,I
        real(kind=8), intent(in), dimension(1+neq*(neq-1)/2) :: K
        real(kind=8), intent(out), dimension(nrowpd,neq) :: pd
        real(kind=8) :: cphi,sphi,cpsi,spsi
        !      write(6,*) 'jac: ',nrowpd,neq,size(K)
!        if (IS_LINEAR_STOKES==1) then
        cphi=cos(I(3)); sphi=sin(I(3)); cpsi=cos(I(4)); spsi=sin(I(4))
        pd(1,1)=-K(1)
        pd(1,2)=-cphi*spsi*K(2)-sphi*spsi*K(3)-cpsi*K(4)
        pd(1,3)=-I(2)*(-sphi*spsi*K(2)+cphi*spsi*K(3))
        pd(1,4)=-I(2)*(cphi*cpsi*K(2)+cpsi*sphi*K(3)-spsi*K(4))
        pd(2,1)=pd(1,2)
        pd(2,2)=-K(1)
        pd(2,3)=I(1)*spsi*(sphi*K(2)-cphi*K(3))+spsi*(cphi*j(3)-sphi*j(2))
        pd(2,4)=-I(1)*cpsi*(cphi*K(2)+sphi*K(3))+cpsi*(cphi*j(2)+sphi*j(3))-spsi*(j(4)-I(1)*K(4))
        pd(3,1)=1d0/I(2)*(sphi/spsi*K(2)-cphi/spsi*K(3))
!        pd(3,1)=1d0/spsi/I(2)*(sphi*K(2)-cphi*K(3))
        pd(3,2)=(-cphi*cpsi/spsi*K(5)-cpsi/spsi*sphi*K(6)+K(7))/I(2)-(-cphi/spsi*(I(1)*K(3)-& 
             j(3)+I(2)*cpsi*K(5))-1d0/spsi*sphi*(-I(1)*K(2)+j(2)+I(2)*cpsi*K(6))+I(2)*K(7))/I(2)**2.
!        pd(3,2)=1d0/spsi/I(2)**2d0*(cphi*(I(1)*K(3)-j(3))+sphi*(j(2)-I(1)*K(2)))
        pd(3,3)=(1d0/spsi*sphi*(I(1)*K(3)-j(3)+I(2)*cpsi*K(5))-cphi/spsi*(-I(1)*K(2)+j(2)+I(2)*cpsi*K(6)))/I(2)
!        pd(3,3)=1d0/I(2)/spsi*(sphi*(I(1)*K(3)-j(3))+cphi*(I(1)*K(2)-j(2)))+1d0/spsi*(sphi*K(5)-cphi*K(6))
        pd(3,4)=(I(2)*cphi*K(5)+cphi*cpsi/spsi/spsi*(I(1)*K(3)-j(3)+I(2)*cpsi*K(5))+I(2)*sphi*K(6)+ &
             cpsi/spsi/spsi*sphi*(-I(1)*K(2)+j(2)+I(2)*cpsi*K(6)))/I(2)
!        pd(3,4)=cpsi/spsi/spsi/I(2)*(cphi*(I(1)*K(3)-j(3))-sphi*(I(1)*K(2)-j(2)))+1d0/spsi/spsi*(cphi*K(5)+sphi*K(6))
        pd(4,1)=(-cphi*cpsi*K(2)-cpsi*sphi*K(3)+spsi*K(4))/I(2)
        pd(4,2)=(-sphi*K(5)+cphi*K(6))/I(2)-(-I(1)*cphi*cpsi*K(2)+spsi*(I(1)*K(4)-j(4))-sphi*(I(1)*cpsi*K(3) &
             -cpsi*j(3)+I(2)*K(5))+cphi*(cpsi*j(2)+I(2)*K(6)))/I(2)**2d0
!        pd(4,2)=1d0/I(2)**2d0*(cphi*cpsi*(I(1)*K(2)-j(2))-spsi*(I(1)*K(4)-j(4))+sphi*cpsi*(I(1)*K(3)-j(3)))
        pd(4,3)=(I(1)*cpsi*sphi*K(2)-cphi*(I(1)*cpsi*K(3)-cpsi*j(3)+I(2)*K(5))-sphi*(cpsi*j(2)+I(2)*K(6)))/I(2)
!        pd(4,3)=1d0/I(2)*(cpsi*sphi*(I(1)*K(2)-j(2))-cpsi*cphi*(I(1)*K(3)-j(3)))-cphi*K(5)-sphi*K(6)
        pd(4,4)=(I(1)*cphi*spsi*K(2)-cphi*spsi*j(2)-sphi*(-I(1)*spsi*K(3)+spsi*j(3))+cpsi*(I(1)*K(4)-j(4)))/I(2)
!        pd(4,4)=1d0/I(2)*(cphi*spsi*(I(1)*K(2)-j(2))+sphi*spsi*(I(1)*K(3)-j(3))+cpsi*(I(1)*K(4)-j(4)))
!           pd=-1d0*pd
!        endif
!         write(6,*) 'pd: ',pd
        return
      end subroutine radtrans_jac_form_sph_old

      subroutine radtrans_lsoda_calc_jac_sph(neq,lam,I,ml &
           ,mu,pd,nrowpd)
! Compute Jacobian for LSODA
        integer, intent(in) :: neq, nrowpd
        real(kind=8), intent(in) :: lam
        real(kind=8), intent(in), dimension(neq) :: I
        real(kind=8), intent(in) :: ml
        real(kind=8), intent(in) :: mu
        real(kind=8), intent(out), dimension(nrowpd,neq) :: pd
        real(kind=8), dimension(neq) :: j
        real(kind=8), dimension(1+neq*(neq-1)/2) :: K
!         write(6,*) 'jac: ',nrowpd
        call radtrans_aux(neq,lam,j,K)
! changed from above to add I as argument since not linear eqns any more
        call radtrans_jac_form_sph(neq,j,K,I,nrowpd,pd)
        !         write(6,*) 'pd: ', pd
        return
      end subroutine radtrans_lsoda_calc_jac_sph

      subroutine calc_O(a,rho,dx,identity,O,M1,M2,M3,M4)
        real(kind=8), dimension(4), intent(in) :: a
        real(kind=8), dimension(3), intent(in) :: rho
        real(kind=8), dimension(4,4) :: identity,onopol!,M1,M2,M3,M4
        real(kind=8), dimension(4,4), intent(out) :: M1,M2,M3,M4,O
        real(kind=8) :: lam1,lam2,ap,theta,sig,a2,p2
        real(kind=8) :: aq,au,av,rhoq,rhou,rhov
        real(kind=8) :: dx,coshlam1,coslam2,sinhlam1,sinlam2
        onopol = exp(-a(1)*dx)
        aq = a(2); au = a(3); av = a(4)
        rhoq = rho(1); rhou = rho(2); rhov = rho(3)
        a2 = aq**2d0+au**2d0+av**2d0
        p2 = rhoq**2d0+rhou**2d0+rhov**2d0
        if(a2.eq.0d0.and.p2.eq.0d0) then
           O = identity*onopol
        else
           ap = aq*rhoq+au*rhou+av*rhov
           lam1 = sqrt(sqrt((a2-p2)**2d0/4d0+ap**2d0)+(a2-p2)/2d0)
           lam2 = sqrt(sqrt((a2-p2)**2d0/4d0+ap**2d0)-(a2-p2)/2d0)
           theta = lam1**2d0+lam2**2d0
           sig = sign(1d0,ap)
           M1 = identity
           M2(:,1) = (/0d0,lam2*aq-sig*lam1*rhoq,lam2*au-sig*lam1*rhou,lam2*av-sig*lam1*rhov/)
           M2(:,2) = (/lam2*aq-sig*lam1*rhoq,0d0,-sig*lam1*av-lam2*rhov,sig*lam1*au+lam2*rhou/)
           M2(:,3) = (/lam2*au-sig*lam1*rhou,sig*lam1*av+lam2*rhov,0d0,-sig*lam1*aq-lam2*rhoq/)
           M2(:,4) = (/lam2*av-sig*lam1*rhov,-sig*lam1*au-lam2*rhou,sig*lam1*aq+lam2*rhoq,0d0/)
           M2 = 1d0/theta*M2
           M3(:,1) = (/0d0,lam1*aq+sig*lam2*rhoq,lam1*au+sig*lam2*rhou,lam1*av+sig*lam2*rhov/)
           M3(:,2) = (/lam1*aq+sig*lam2*rhoq,0d0,sig*lam2*av-lam1*rhov,-sig*lam2*au+lam1*rhou/)
           M3(:,3) = (/lam1*au+sig*lam2*rhou,-sig*lam2*av+lam1*rhov,0d0,sig*lam2*aq-lam1*rhoq/)
           M3(:,4) = (/lam1*av+sig*lam2*rhov,sig*lam2*au-lam1*rhou,-sig*lam2*aq+lam1*rhoq,0d0/)
           M3=1d0/theta*M3
           M4(:,1) = (/(a2+p2)/2d0,au*rhov-av*rhou,av*rhoq-aq*rhov,aq*rhou-au*rhoq/)
           M4(:,2) = (/av*rhou-au*rhov,aq*aq+rhoq*rhoq-(a2+p2)/2d0,aq*au+rhoq*rhou,av*aq+rhov*rhoq/)
           M4(:,3) = (/aq*rhov-av*rhoq,aq*au+rhoq*rhou,au*au+rhou*rhou-(a2+p2)/2d0,au*av+rhou*rhov/)
           M4(:,4) = (/au*rhoq-aq*rhou,av*aq+rhov*rhoq,au*av+rhou*rhov,av*av+rhov*rhov-(a2+p2)/2d0/)
           M4=2d0/theta*M4
           coshlam1=cosh(lam1*dx)
           coslam2=cos(lam2*dx)
           sinhlam1=sinh(lam1*dx)
           sinlam2=sin(lam2*dx)
           O = onopol*(1d0/2d0*(coshlam1+coslam2)*M1 - sinlam2*M2-sinhlam1*M3+1d0/2d0 &
                *(coshlam1-coslam2)*M4)
        endif
      end subroutine calc_O

      subroutine imatrix_4_single(a,b)
        real(kind=8), dimension(4,4), intent(in) :: a
        real(kind=8), dimension(4,4), intent(out) :: b
        real(kind=8) :: detA,a11,a21,a31,a41,a12,a22, &
             a32,a42,a13,a23,a33,a43,a14,a24,a34,a44

        b = 0d0

        a11 = a(1,1); a12 = a(1,2); a13 = a(1,3); a14 = a(1,4)
        a21 = a(2,1); a22 = a(2,2); a23 = a(2,3); a24 = a(2,4)
        a31 = a(3,1); a32 = a(3,2); a33 = a(3,3); a34 = a(3,4)
        a41 = a(4,1); a42 = a(4,2); a43 = a(4,3); a44 = a(4,4)
        
        b(1,1) = a22*a33*a44 + a23*a34*a42 + a24*a32*a43  &
             - a22*a34*a43 - a23*a32*a44 - a24*a33*a42
        b(1,2) = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 &
             - a12*a33*a44 - a13*a34*a42 - a14*a32*a43
        b(1,3) = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 &
             - a12*a24*a43 - a13*a22*a44 - a14*a23*a42
        b(1,4) = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 &
             - a12*a23*a34 - a13*a24*a32 - a14*a22*a33
        b(2,1) = a21*a34*a43 + a23*a31*a44 + a24*a33*a41 &
             - a21*a33*a44 - a23*a34*a41 - a24*a31*a43
        b(2,2) = a11*a33*a44 + a13*a34*a41 + a14*a31*a43 &
             - a11*a34*a43 - a13*a31*a44 - a14*a33*a41
        b(2,3) = a11*a24*a43 + a13*a21*a44 + a14*a23*a41 &
             - a11*a23*a44 - a13*a24*a41 - a14*a21*a43
        b(2,4) = a11*a23*a34 + a13*a24*a31 + a14*a21*a33 &
             - a11*a24*a33 - a13*a21*a34 - a14*a23*a31
        b(3,1) = a21*a32*a44 + a22*a34*a41 + a24*a31*a42 &
             - a21*a34*a42 - a22*a31*a44 - a24*a32*a41
        b(3,2) = a11*a34*a42 + a12*a31*a44 + a14*a32*a41 &
             - a11*a32*a44 - a12*a34*a41 - a14*a31*a42
        b(3,3) = a11*a22*a44 + a12*a24*a41 + a14*a21*a42 &
             - a11*a24*a42 - a12*a21*a44 - a14*a22*a41
        b(3,4) = a11*a24*a32 + a12*a21*a34 + a14*a22*a31 &
             - a11*a22*a34 - a12*a24*a31 - a14*a21*a32
        b(4,1) = a21*a33*a42 + a22*a31*a43 + a23*a32*a41 &
             - a21*a32*a43 - a22*a33*a41 - a23*a31*a42
        b(4,2) = a11*a32*a43 + a12*a33*a41 + a13*a31*a42 &
             - a11*a33*a42 - a12*a31*a43 - a13*a32*a41
        b(4,3) = a11*a23*a42 + a12*a21*a43 + a13*a22*a41 &
             - a11*a22*a43 - a12*a23*a41 - a13*a21*a42
        b(4,4) = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 &
             - a11*a23*a32 - a12*a21*a33 - a13*a22*a31
        detA = a(1,1)*b(1,1) + a(2,1)*b(1,2) + &
             a(3,1)*b(1,3) + a(4,1)*b(1,4)
        b = b/detA
      end subroutine imatrix_4_single
      
      subroutine opacity_matrix(a,p,Karr)
      real(kind=8), dimension(:,:), intent(in) :: a
      real(kind=8), dimension(:,:), intent(in) :: p
      real(kind=8), dimension(size(a,1),4,4), intent(out) :: Karr
      Karr(:,1,1) = a(:,1); Karr(:,2,2) = a(:,1); Karr(:,3,3) = a(:,1); Karr(:,4,4) = a(:,1)
      Karr(:,2,1) = a(:,2); Karr(:,3,1) = a(:,3); Karr(:,4,1) = a(:,4)
      Karr(:,1,2) = a(:,2); Karr(:,1,3) = a(:,3); Karr(:,1,4) = a(:,4)
      Karr(:,3,2) = -p(:,3); Karr(:,4,2) = p(:,2); Karr(:,2,3) = p(:,3); Karr(:,2,4) = -p(:,2)
      Karr(:,4,3) = -p(:,1); Karr(:,3,4) = p(:,1)
    end subroutine opacity_matrix

      subroutine invert_delo_matrix_thin_single(dx,K0,ki,delta,identity,matrix,imatrix)
        real(kind=8), intent(in) :: dx,ki,delta
        real(kind=8), intent(in), dimension(4,4) :: identity,K0
        real(kind=8), intent(out), dimension(4,4) :: matrix,imatrix
        matrix = (1d0-delta/2d0+delta**2d0/6d0)*identity & 
             +(0.5d0*dx-1d0/6d0*dx**2d0*ki)*K0
        call imatrix_4(matrix,imatrix)
      end subroutine invert_delo_matrix_thin_single

      subroutine calc_delo_P_thin_single(imatrix,dx,j,j1,ki,ki1,P)
        real(kind=8), intent(in), dimension(4,4) :: imatrix
        real(kind=8), intent(in) :: dx,ki,ki1
        real(kind=8), intent(in), dimension(4) :: j,j1
        real(kind=8), intent(out), dimension(4) :: P
        P = matmul(imatrix,(0.5d0*dx*j-1d0/6d0*dx**2d0*ki*j)+ &
             (0.5d0*j1*dx-1d0/3d0*dx**2d0*ki*j1))
      end subroutine calc_delo_P_thin_single

      subroutine calc_delo_Q_thin_single(imatrix,dx,ki,ki1,K1,identity,Q)
        real(kind=8), intent(in), dimension(4,4) :: imatrix,K1,identity
        real(kind=8), intent(in) :: dx,ki,ki1
        real(kind=8), intent(out), dimension(4,4) :: Q
        Q = matmul(imatrix,identity*(1d0-0.5d0*dx*ki+1d0/6d0*dx**2d0*ki**2d0) &
             -(0.5d0*dx-1d0/3d0*dx**2d0)*K1)
      end subroutine calc_delo_Q_thin_single

      subroutine invert_delo_matrix_single(F,G,Kp,identity,matrix,imatrix)
        real(kind=8), intent(in) :: F,G
        real(kind=8), intent(in), dimension(4,4) :: Kp,identity
        real(kind=8), intent(out), dimension(4,4) :: matrix,imatrix
        matrix = identity+(F-G)*Kp
        call imatrix_4(matrix,imatrix)
      end subroutine invert_delo_matrix_single

      subroutine calc_delo_P_single(imatrix,F,G,Sp,Sp1,P)
        real(kind=8), dimension(4,4), intent(in) :: imatrix
        real(kind=8), intent(in) :: F,G
        real(kind=8), intent(in), dimension(4) :: Sp,Sp1
        real(kind=8), intent(out), dimension(4) :: P
        P = matmul(imatrix,(F-G)*Sp+G*Sp1)
      end subroutine calc_delo_P_single

      subroutine calc_delo_Q_single(imatrix,E,F,G,Kp1,identity,Q)
        real(kind=8), dimension(4,4), intent(in) :: imatrix,Kp1,identity
        real(kind=8), intent(in) :: E,F,G
        real(kind=8), dimension(4,4), intent(out) :: Q
        Q = matmul(imatrix,identity*E-G*Kp1)
      end subroutine calc_delo_Q_single

      subroutine radtrans_integrate_delo(x,tau,j,a,rho,P,Q,im)
        real(kind=8), dimension(:), intent(in) :: x,tau
        real(kind=8), dimension(:,:), intent(in) :: j,a,rho
        real(kind=8), dimension(size(x),4), intent(inout) :: P
        real(kind=8), dimension(size(x),4,4), intent(inout) :: Q,im
        real(kind=8), dimension(size(x)-1) :: delta,dx
        real(kind=8) :: E,F,G
        real(kind=8), dimension(4) :: iprev,Sp,Sp1,pt,I0
        real(kind=8), dimension(4,4) :: identity,Kp,Kp1,K0,K1,qt,imatrix,matrix
        real(kind=8), dimension(size(x),4,4) :: Karr
        integer :: k
        I0=0d0
        identity = reshape((/1d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,0d0, &
             1d0,0d0,0d0,0d0,0d0,1d0/),(/4,4/))
        delta = tau(2:npts) - tau(1:npts-1)
        call opacity_matrix(a,rho,Karr)
        dx = x(1:npts-1) - x(2:npts)
        !        write(6,*) 'delo integrate: ',thin,minval(dx),maxval(tau),minval(delta),maxval(delta)
! integration is from deepest point out for starting intensity I0
        intensity(:,1) = I0; iprev = I0
! alternative way to write this: where(delta.gt.thin) and calculate pt, qt, elsewhere() and calculate ptt, qtt and then do loop free of if statements
        do k=npts-1,1,-1
           K0 = Karr(k,:,:)
           K1 = Karr(k+1,:,:)
           if (delta(k).gt.thin) then
              E = exp(-delta(k))
              F = 1d0-E
              G = (1d0-(1d0+delta(k))*E)/delta(k)
              Sp = j(k,:)/a(k,1)
              Sp1 = j(k+1,:)/a(k+1,1)
              Kp = K0/a(k,1)-identity; Kp1 = K1/a(k+1,1)-identity
              call invert_delo_matrix(F,G,Kp,identity,matrix,imatrix)
              call calc_delo_P(imatrix,F,G,Sp,Sp1,pt)
              call calc_delo_Q(imatrix,E,F,G,Kp1,identity,qt)
           else
              call invert_delo_matrix_thin(dx(k),K0,a(k,1),delta(k),identity,matrix,imatrix)
              call calc_delo_P_thin(imatrix,dx(k),j(k,:),j(k+1,:),a(k,1),a(k+1,1),pt)
              call calc_delo_Q_thin(imatrix,dx(k),a(k,1),a(k+1,1),K1,identity,qt)
           endif
           intensity(:,npts-k+1) = pt + matmul(qt,iprev)
           P(k,:) = pt
           Q(k,:,:) = qt
           im(k,:,:) = imatrix
           iprev = intensity(:,npts-k+1)
        end do
!        write(6,*) 'delo intensity: ',intensity(1,:)
!        write(6,*) 'delo intensity P: ',P(:,1)
      end subroutine radtrans_integrate_delo

      subroutine radtrans_integrate_formal(x,j,a,rho,O)
        real(kind=8), dimension(:), intent(in) :: x
        real(kind=8), dimension(size(x),4), intent(in) :: j,a
        real(kind=8), dimension(size(x),3), intent(in) :: rho
        real(kind=8), dimension(4,4,size(x)), intent(inout) :: O
        real(kind=8), dimension(size(x)-1) :: dx
        real(kind=8), dimension(4) :: I0,iprev
        real(kind=8), dimension(4,4) :: identity,M1,M2,M3,M4,Ot
        integer :: k
        I0=0d0
        identity = reshape((/1d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,0d0, &
             1d0,0d0,0d0,0d0,0d0,1d0/),(/4,4/))
        dx = x(1:npts-1) - x(2:npts)
        intensity(:,1) = I0; iprev = I0
        do k=npts-1,1,-1
!           call calc_O(a(k,:),rho(k,:),dx(k),identity,Ot,M1,M2,M3,M4)
!           intensity(:,npts-k+1) = matmul(Ot,j(k,:))*dx(k)+matmul(Ot,iprev)
           call calc_O(a(k,:),rho(k,:),dx(k),identity,Ot,M1,M2,M3,M4)
           intensity(:,npts-k+1) = matmul(Ot,j(k,:))*dx(k)+matmul(Ot,iprev)
           iprev = intensity(:,npts-k+1)
           O(:,:,k)=Ot
        end do
        return
      end subroutine radtrans_integrate_formal

      subroutine radtrans_integrate_quadrature(s,j,K)
        real(kind=8), dimension(:), intent(in) :: s,j,K
! quadrature solution assuming I0=0, - sign from s,j,tau being tabulated backwards:
        intensity(1,:)=-tsum(s,j*exp(-tau))
      end subroutine radtrans_integrate_quadrature

      subroutine del_radtrans_integrate_data()
        deallocate(jj); deallocate(KK)
        deallocate(s); deallocate(ss); deallocate(s0); deallocate(tau)
        deallocate(intensity)
        if(iflag==1) then
           deallocate(PP); deallocate(QQ); deallocate(imm)
        elseif(iflag==2) then
           deallocate(OO)
        elseif(iflag==3) then
           deallocate(stokesq);deallocate(stokesu);deallocate(stokesv)
        endif
      end subroutine del_radtrans_integrate_data
   
    end module radtrans_integrate
