  module radtrans_integrate

    use odepack, only: lsoda_basic
    use interpolate, only: get_weight, locate

    implicit none

    ! global data to avoid needing to pass arguments through external integration routines (e.g., lsoda)
    integer :: lindx=25,nequations,nptstot,npts,nptsout,iflag
    double precision, dimension(:,:), allocatable :: KK,jj,intensity,PP
    double precision, dimension(:), allocatable :: s,ss,s0,tau
    double precision, dimension(:,:,:), allocatable :: QQ,imm

    integer :: IS_LINEAR_STOKES = 1
    integer, dimension(3) :: stats
    double precision :: ortol = 1d-6, oatol = 1d-8, hmax = 10, MAX_TAU, MAX_TAU_DEFAULT = 10d0, thin = 1d-2

!$omp threadprivate(ss,tau,s0,jj,KK,intensity,lindx,nptstot,npts,nptsout,s,stats)
!$omp threadprivate(ortol,oatol,hmax,thin,MAX_TAU,QQ,PP,imm)

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

    contains

      subroutine init_radtrans_integrate_data(riflag,rneq,gnpts,rnpts,maxtau,hm,oa,ort,th)
        integer, intent(in) :: riflag,rneq,gnpts,rnpts
!        double precision, intent(in), optional :: maxtau,hm,oa,ort
        double precision, intent(in), optional :: maxtau,hm,oa,ort,th
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
        allocate(tau(npts))
        allocate(s0(npts))
        allocate(jj(npts,nequations))
        allocate(KK(npts,1+(nequations)*(nequations-1)/2))
        allocate(intensity(nequations,npts))
        allocate(s(npts)); allocate(ss(npts))
        s = 0d0
        ss = 0d0
        intensity = 0d0
        jj = 0d0
        tau = 0d0
        KK = 0d0
        s0 = 0d0
        if(iflag==1) then
           allocate(QQ(npts,nequations,nequations)); allocate(PP(npts,nequations))
           allocate(imm(npts,nequations,nequations))
           QQ = 0d0; PP = 0d0; imm = 0d0
        endif
      end subroutine init_radtrans_integrate_data

      subroutine integrate(sinput,ej,eK,tauin,rnpts)
! these array sizes are really all known in terms of npts, nequations...
        integer, intent(inout) :: rnpts
!        double precision, dimension(:,:), intent(inout) :: rI
        double precision, dimension(:,:), intent(in) :: ej,eK
        double precision, dimension(:), intent(in) :: sinput
        double precision, dimension(:), intent(in) :: tauin
!        double precision, dimension(:), intent(inout) :: I0
!        double precision, intent(in), optional :: maxt
    ! Subroutine to call integration routine for RT equation
    ! JAD 2/28/2011 moved from grtrans_driver 8/7/2014
     ! assign global data
!        call init_radtrans_integrate_data()
!        allocate(tau(gnpts))
!        allocate(s0(size(sinput)))
 !       allocate(jj(size(ej,1),size(ej,2)))
!        allocate(KK(size(eK,1),size(eK,2)))
!        allocate(intensity(rneq,rnpts))
!        allocate(s(size(sinput))); allocate(ss(size(sinput)))
!        write(6,*) 'radtrans_integrate: ',sinput
        s0=sinput; s=sinput; ss=sinput; intensity=0d0
!        write(6,*) 'radtrans_integrate: ',s
        tau = tauin
!        write(6,*) 'made it to integrate'
        jj=ej; KK=eK
        if (iflag==0) then
           call radtrans_integrate_lsoda()
        elseif (iflag==1) then
           call radtrans_integrate_delo(s,tau,jj,KK(:,1:4),KK(:,5:7),PP,QQ,imm)
        endif
!        write(6,*) 'assign intensity', size(rI,1), size(rI,2), size(intensity,1), size(intensity,2), nptsout
!        rI(1:nequations,1:nptsout) = intensity(1:nequations,1:nptsout);
        rnpts = nptsout
!        write(6,*) 'intensity: ', rnpts, rI(rnpts,1)
        if(isnan(intensity(1,rnpts))) then
           write(6,*) 'NaN in integrate ej: ',ej
           write(6,*) 'NaN in integrate jj 2: ',jj
           write(6,*) 'NaN in integrate eK: ',eK
           write(6,*) 'NaN in integrate KK: ',KK
        endif
!        deallocate(jj); deallocate(KK)
!        deallocate(s); deallocate(ss); deallocate(s0); deallocate(tau)
!        deallocate(intensity)
        return
      end subroutine integrate
  
      subroutine radtrans_integrate_lsoda()
!        double precision, dimension(:), intent(inout) :: I0
        double precision, dimension(4) :: I0
        double precision, dimension(npts) :: dummy
        double precision, dimension(:,:), allocatable :: tau_arr
        double precision, dimension(:), allocatable :: tau_temp,intvals
        integer, dimension(:), allocatable :: inds
        double precision :: weight
        integer :: lamdex,i,ii,i1,i2,taudex
      !         
!         write(6,*) 'Kmax: ',size(KK), size(jj)
!         write(6,*) 'tau: ',tau(2),maxval(tau),MAX_TAU
!         fac=max(1d0,tau(npts)); jj=jj*fac
        i1=1; i2=nptstot
        I0 = 0d0
        if(maxval(tau).le.MAX_TAU) then
           lamdex=npts
        else
!           write(6,*) 'locate', size(tau), maxval(tau), tau(2)
           call locate(tau,MAX_TAU,lamdex)
           lamdex=lamdex+1
         !           write(6,*) 'locate', lamdex
           if(lamdex==2) then
              call get_weight(tau,MAX_TAU,lamdex,weight)
! need to re-work this. just trying to interpolate between s0(1) and s0(2).
!             write(6,*) 'weight: ',lamdex,weight,s(lamdex),s(lamdex+1)
              s=(/(s(1)-s(2))*(1d0-weight),0d0/)
              lamdex=lamdex+1
!             write(6,*) 's: ',s(1:lamdex),lamdex
           endif
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
         !           s(i1:i2)=s(i2:i1:-1)
!            write(6,*) 'lsoda basic'
           ! test I0:
!           I0(1) = jj(npts-i1,1)/KK(npts-i1,1)
           call lsoda_basic(radtrans_lsoda_calc_rhs,I0(1:nequations), &
                s(i1:i2),oatol, &
                ortol,radtrans_lsoda_calc_jac,intensity(:,i1:i2), &
                1,100000,stats,hmax=hmax)
!           write(6,*) 'after lsoda'
        else
!           write(6,*) 'rneq: ',neq,size(I)
!           write(6,*) 'si: ',s(lamdex),jj(lamdex,1),KK(lamdex,1)
!           write(6,*) 'si: ',bcgs(lamdex),ncgs(lamdex)
!           s(i1:i2)=s(npts)-s(i2:i1:-1)
!            s(i1:i2)=s(i2:i1:-1)
!           write(6,*) 's: ',s,i1,i2
           call lsoda_basic(radtrans_lsoda_calc_rhs_npol, &
              I0(1:nequations),s(i1:i2),oatol, &
              ortol,radtrans_lsoda_calc_jac_npol,intensity(:,i1:i2), &
              1,100000,stats,hmax=hmax)
!           call lsoda_basic(radtrans_lsoda_calc_rhs_npol,
!     &      I0(1:neq),s(npts)-s(lamdex:1:-1),oatol,
!     &      ortol,radtrans_lsoda_calc_jac_npol,intensity(:,lamdex:1:-1),
!     &      1,100000,stats,hmax=10d0)
!          write(6,*) 'after lsoda',intensity(:,lamdex)
!          npts=lamdex
! quadrature solution:
           !dummy=tsum(s0(1:npts),jj(:,1)*exp(-tau(1:npts)))
!           do i=1,lamdex; write(6,*) 'tsum: ',s0(i),tau(i),dummy(i); enddo 
!           write(6,*) 'dummy: ',dummy(npts),tau(npts)
!            write(6,*) 'rshift: ',rshift
!           write(6,*) 'occu: ',intensity(1,i1:i2)/rshift(i1:i2)**3
            
!           write(6,*) 'tsum: ',tsum(s0(1:npts),jj(:,1)*exp(-tau(1:npts)))
!           write(6,*) 'tsum tau: ',tsum(tau(i1:i2),jj(:,1)/KK(:,1)*exp(-tau(i1:i2)))
!           write(6,*) 'tsum diff: ',(intensity(1,i1:i2)-tsum(s0(1:npts),
!     &     jj(:,1)*exp(-tau(1:npts))))/tsum(s0(1:npts),jj(:,1)*
!     &      exp(-tau(1:npts)))
!           write(6,*) 'didlam: ',jj(i1:i2,1)-KK(i1:i2,1)*
!     &        tsum(s0(i1:i2),
!     &       jj(i1:i2,1)*exp(-tau(i1:i2)))
        endif
!         do while((intensity(1,npts).lt.0d0).or.isnan(intensity(1,npts)))
!           npts=max(npts-1,1)
!           if(npts==1) then
!             npts=lamdex; I0=0d0; s=s/2d0
!              write(6,*) 'repeat: ',s
!             if(neq==4) then
!               call lsoda_basic(radtrans_lsoda_calc_rhs,
!     &          I0(1:neq),s(i1:i2),oatol,
!     &          ortol,radtrans_lsoda_calc_jac,intensity(:,i1:i2),
!     &          1,100000,stats,hmax=10d0)
!             else
!               write(6,*) 'lsoda basic: ',neq
!               call lsoda_basic(radtrans_lsoda_calc_rhs_npol,
!     &          I0(1:neq),s(i1:i2),oatol,
!     &          ortol,radtrans_lsoda_calc_jac_npol,
!     &          intensity(:,i1:i2),
!     &          1,100000,stats,hmax=10d0)
!             endif                   
!           endif
!         enddo
! Negative total intensity
!         write(6,*) 'made it'
!         write(6,*) 'r: ',size(I), maxval(I)
!         write(6,*) 'integrate: ',intensity(1,:)-s0
! Now make intensity from occupation number:
 !        do i=1,neq
 !          intensity(i,i1:i2)=intensity(i,i1:i2)/
 !    &           rshift(i1:i2)**3
         !         enddo
!        deallocate(ss)
        if(isnan(intensity(1,i2))) then
           write(6,*) 'NaN in integrate: ',i1,i2,s(i1:i2)
           write(6,*) 'NaN in integrate intensity: ',intensity(1,i1:i2)
           write(6,*) 'NaN in integrate j: ',jj(i1:i2,:)
           write(6,*) 'NaN in integrate K: ',KK(i1:i2,:)
        endif
        return
      end subroutine radtrans_integrate_lsoda

      subroutine radtrans_lsoda_calc_rhs(neq,lam,I,dIdlam)
! Compute RHS dIdlam for LSODA
        integer, intent(in) :: neq
        double precision, intent(in) :: lam
        double precision, intent(in), dimension(neq) :: I
        double precision, intent(out), dimension(neq) :: dIdlam
!         double precision, intent(out), dimension(neq,neq) :: jac
        double precision, dimension(neq) :: j
        double precision, dimension(1+neq*(neq-1)/2) :: K
        double precision :: rshift
        call radtrans_aux(neq,lam,j,K,rshift)
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
        double precision, intent(in) :: lam
        double precision, intent(in), dimension(neq) :: I
        double precision, intent(out), dimension(neq) :: dIdlam
!         double precision, intent(out), dimension(neq,neq) :: jac
        double precision, dimension(neq) :: j
        double precision, dimension(1+neq*(neq-1)/2) :: K
        double precision :: rshift
        call radtrans_aux(neq,lam,j,K,rshift)
        call radtrans_rhs_form_npol(neq,j,K,rshift,dIdlam,I)
!         didlam=didlam*1e25
!         write(6,*) 'dIdlam: ',lam,dIdlam
!         write(6,*) 'jk: ',j(1),j(2),j(4)
!         write(6,*) 'K: ',K(1),K(4),
!     &    K(5),K(7)
      end subroutine radtrans_lsoda_calc_rhs_npol
         
      subroutine radtrans_rhs_form(neq,j,K,dIdlam,I)
        integer, intent(in) :: neq
        double precision, intent(in), dimension(neq) :: j
        double precision, intent(in), dimension(1+neq*(neq-1)/2) :: K
        double precision, intent(out), dimension(neq) :: dIdlam
        double precision, intent(in), dimension(neq) :: I
      !         write(6,*) 'rhs: ',IS_LINEAR_STOKES,size(I),size(K),size(J)
        if (IS_LINEAR_STOKES==1) then
           dIdlam(1)=j(1)-(K(1)*I(1)+K(2)*I(2)+K(3)*I(3)+K(4)*I(4))
           dIdlam(2)=j(2)-(K(2)*I(1)+K(1)*I(2)+K(7)*I(3)-K(6)*I(4))
           dIdlam(3)=j(3)-(K(3)*I(1)-K(7)*I(2)+K(1)*I(3)+K(5)*I(4))
           dIdlam(4)=j(4)-(K(4)*I(1)+K(6)*I(2)-K(5)*I(3)+K(1)*I(4))
        endif
      end subroutine radtrans_rhs_form

      subroutine radtrans_rhs_form_npol(neq,j,K,rshift,dIdlam,I)
        integer, intent(in) :: neq
        double precision, intent(in), dimension(neq) :: j
        double precision, intent(in), dimension(1+neq*(neq-1)/2) :: K
        double precision, intent(out), dimension(neq) :: dIdlam
        double precision, intent(in), dimension(neq) :: I
        double precision, intent(in) :: rshift
!         write(6,*) 'rhs npol: ',size(I),size(K),size(J)
!         dIdlam(1)=maxval((/j(1)-K(1)*I(1),0d0/))
        dIdlam(1)=j(1)-K(1)*I(1)
!         write(6,*) 'rhs: ',dIdlam(1),j(1),K(1),I(1)
        return
      end subroutine radtrans_rhs_form_npol

      subroutine radtrans_jac_form_npol(neq,j,K,nrowpd,pd)
        integer, intent(in) :: neq, nrowpd
        double precision, intent(in), dimension(neq) :: j
        double precision, intent(in), dimension(1+neq*(neq-1)/2) :: K
        double precision, intent(out), dimension(nrowpd,neq) :: pd
   !      write(6,*) 'jac: ',nrowpd,neq,size(K),K(1)
        pd(1,1)=-K(1)
!         write(6,*) 'pd: ',pd
        return
      end subroutine radtrans_jac_form_npol

      subroutine radtrans_jac_form(neq,j,K,nrowpd,pd)
        integer, intent(in) :: neq, nrowpd
        double precision, intent(in), dimension(neq) :: j
        double precision, intent(in), dimension(1+neq*(neq-1)/2) :: K
        double precision, intent(out), dimension(nrowpd,neq) :: pd
        !      write(6,*) 'jac: ',nrowpd,neq,size(K)
        if (IS_LINEAR_STOKES==1) then
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
        endif
!         write(6,*) 'pd: ',pd
        return
      end subroutine radtrans_jac_form
       
      subroutine radtrans_lsoda_calc_jac(neq,lam,I,ml &
           ,mu,pd,nrowpd)
! Compute Jacobian for LSODA
        integer, intent(in) :: neq, nrowpd
        double precision, intent(in) :: lam
        double precision, intent(in), dimension(neq) :: I
        double precision, intent(in) :: ml
        double precision, intent(in) :: mu
        double precision, intent(out), dimension(nrowpd,neq) :: pd
        double precision, dimension(neq) :: j
        double precision, dimension(1+neq*(neq-1)/2) :: K
        double precision :: rshift
!         write(6,*) 'jac: ',nrowpd
        call radtrans_aux(neq,lam,j,K,rshift)
        call radtrans_jac_form(neq,j,K,nrowpd,pd)
        !         write(6,*) 'pd: ', pd
        return
      end subroutine radtrans_lsoda_calc_jac

      subroutine radtrans_lsoda_calc_jac_npol(neq,lam,I,ml &
           ,mu,pd,nrowpd)
! Compute Jacobian for LSODA
        integer, intent(in) :: neq, nrowpd
        double precision, intent(in) :: lam
        double precision, intent(in), dimension(neq) :: I
        double precision, intent(in) :: ml
        double precision, intent(in) :: mu
        double precision, intent(out), dimension(nrowpd,neq) :: pd
        double precision, dimension(neq) :: j
        double precision, dimension(1+neq*(neq-1)/2) :: K
        double precision :: rshift
!         write(6,*) 'jac: ',nrowpd
        call radtrans_aux(neq,lam,j,K,rshift)
        call radtrans_jac_form_npol(neq,j,K,nrowpd,pd)
!         write(6,*) 'pd: ', pd
        return
      end subroutine radtrans_lsoda_calc_jac_npol

      subroutine radtrans_aux(neq,lam,j,K,rshift)
        integer, intent(in) :: neq
        double precision, intent(in) :: lam
        double precision, intent(out) :: rshift
        double precision, intent(out), dimension(neq) :: j
        double precision, intent(out), dimension(1+neq*(neq-1)/2) :: K
        double precision :: weight
      !         double precision, dimension(size(s0)) :: dummy
        integer :: indx,uindx
!         call interp_geo_point(lam,g)
!         call get_fluid_vars(x0,f)
!         call calc_emissivity(nu,f,e)
!         call calc_emissivity(e)
 !        call locate(s0,lam,lindx)
!         write(6,*) 'aux: ',lam,size(s0)
!         write(6,*) 'glam: ',s0
 ! NEED TO ADD ASSIGNMENT OF SS
        call get_weight(ss,lam,lindx,weight)
!         call get_weight(s0(1:npts),lam,lindx,weight)
 !        write(6,*) 'sz: ',size(j),size(K),size(j,2),size(K,2)
!          j=exp((1d0-weight)*log(j(lindx,:))+weight*log(j(lindx+1,:)))
!          K=exp((1d0-weight)*log(K(lindx,:))+weight*log(K(lindx+1,:)))
!          rshift=exp((1d0-weight)*log(rshift(lindx))+ &
!               weight*log(rshift(lindx+1)))
        indx=npts-lindx+1; uindx=minval((/indx+1,npts/))
        !          write(6,*) 'j: ',j
        j=(1d0-weight)*jj(indx,:)+weight*jj(uindx,:)
        K=(1d0-weight)*KK(indx,:)+weight*KK(uindx,:)
!          write(6,*) 'j: ',j,weight,indx,uindx,jj(indx,:),jj(uindx,:)
!          write(6,*) 'j: ',lam,lindx,jj
!          j=j(lindx,:)
!          K=K(lindx,:)
!         write(6,*) 'lindx: ',lam,indx,lindx,weight,size(K),size(j)
 !        write(6,*) 'lindx2: ',jj(lindx,1),jj(lindx+1,1),K(lindx,1),K(lindx+1,1),j,K
     !    write(6,*) 'ek: ',K(:,1)
!         evals=interp_emis(lam)
!         j=j
!         K=K
!         call calc_rad_trans_coefs(j,K)
!         j=1d0
!         K=0d0
      end subroutine radtrans_aux

!      subroutine imatrix_4(a,b)
!        double precision, dimension(:,:,:), intent(in) :: a
!        double precision, dimension(size(a,1),4,4), intent(out) :: b
!        double precision, dimension(size(a,1)) :: detA
!        integer :: i
!        double precision, dimension(size(a,1)) :: a11,a21,a31,a41,a12,&
!             a22,a32,a42,a13,a23,a33,a43,a14,a24,a34,a44
!        a11 = a(:,1,1); a12 = a(:,2,1); a13 = a(:,3,1); a14 = a(:,4,1)
!        a21 = a(:,1,2); a22 = a(:,2,2); a23 = a(:,3,2); a24 = a(:,4,2)
!        a31 = a(:,1,3); a32 = a(:,2,3); a33 = a(:,3,3); a34 = a(:,4,3)
!        a41 = a(:,1,4); a42 = a(:,2,4); a43 = a(:,3,4); a44 = a(:,4,4)
        
!        b(:,1,1) = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42
!        b(:,2,1) = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43
!        b(:,3,1) = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42
!        b(:,4,1) = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33
!        b(:,1,2) = a21*a34*a43 + a23*a32*a44 + a24*a33*a41 - a21*a33*a44 - a23*a34*a41 - a24*a31*a43
!        b(:,2,2) = a11*a33*a44 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 - a14*a33*a41
!        b(:,3,2) = a11*a24*a43 + a13*a21*a44 + a14*a23*a41 - a11*a23*a44 - a13*a24*a41 - a14*a21*a43
!        b(:,4,2) = a11*a23*a34 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 - a14*a23*a31
!        b(:,1,3) = a21*a32*a44 + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 - a24*a32*a41
!        b(:,2,3) = a11*a34*a42 + a12*a31*a44 + a14*a32*a41 - a11*a32*a44 - a12*a34*a41 - a14*a31*a42
!        b(:,3,3) = a11*a22*a44 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 - a14*a22*a41
!        b(:,4,3) = a11*a24*a32 + a12*a21*a34 + a14*a22*a31 - a11*a22*a34 - a12*a24*a31 - a14*a21*a32
!        b(:,1,4) = a21*a33*a42 + a22*a31*a43 + a23*a32*a41 - a21*a32*a43 - a22*a33*a41 - a23*a31*a42
!        b(:,2,4) = a11*a32*a43 + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 - a13*a32*a41
!        b(:,3,4) = a11*a23*a42 + a12*a21*a43 + a13*a22*a41 - a11*a22*a43 - a12*a23*a41 - a13*a21*a42
!        b(:,4,4) = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31

!        detA = a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 + a12*a21*a34*a43 + a12*a23*a31*a44 &
!             + a12*a24*a33*a41 + a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 + a14*a21*a33*a42 &
!             + a14*a22*a31*a43 + a14*a23*a32*a41 - a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42 &
!             - a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 - a13*a21*a34*a42 - a13*a22*a31*a44 &
!             - a13*a24*a32*a41 - a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42

!        detA = a(:,1,1)*b(:,1,1) + a(:,2,1)*b(:,1,2) + a(:,3,1)*b(:,1,3) + a(:,4,1)*b(:,1,4)
!        do i=1,size(a,1)
!           b(i,:,:) = 1./detA(i)*b(i,:,:)
!        enddo
!      end subroutine imatrix_4

      subroutine imatrix_4_single(a,b)
        double precision, dimension(4,4), intent(in) :: a
        double precision, dimension(4,4), intent(out) :: b
        double precision :: detA,a11,a21,a31,a41,a12,a22, &
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

!        detA = a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 + a12*a21*a34*a43 + a12*a23*a31*a44 &
!             + a12*a24*a33*a41 + a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42 + a14*a21*a33*a42 &
!             + a14*a22*a31*a43 + a14*a23*a32*a41 - a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42 &
!             - a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43 - a13*a21*a34*a42 - a13*a22*a31*a44 &
!             - a13*a24*a32*a41 - a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42

        detA = a(1,1)*b(1,1) + a(2,1)*b(1,2) + &
             a(3,1)*b(1,3) + a(4,1)*b(1,4)
        b = b/detA
      end subroutine imatrix_4_single
      
      subroutine opacity_matrix(a,p,Karr)
      double precision, dimension(:,:), intent(in) :: a
      double precision, dimension(:,:), intent(in) :: p
      double precision, dimension(size(a,1),4,4), intent(out) :: Karr
!      K = np.array(((a(0),a(1),a(2),a(3)),(a(1),a(0),p(2),-p(1)),(a(2),-p(2),a(0),p(0)),(a(3),p(1),-p(0),a(0))))
      Karr(:,1,1) = a(:,1); Karr(:,2,2) = a(:,1); Karr(:,3,3) = a(:,1); Karr(:,4,4) = a(:,1)
      Karr(:,2,1) = a(:,2); Karr(:,3,1) = a(:,3); Karr(:,4,1) = a(:,4)
      Karr(:,1,2) = a(:,2); Karr(:,1,3) = a(:,3); Karr(:,1,4) = a(:,4)
!      Karr(:,3,2) = p(:,3); Karr(:,4,2) = -p(:,2); Karr(:,2,3) = -p(:,3); Karr(:,2,4) = p(:,2)
!      Karr(:,4,3) = p(:,1); Karr(:,3,4) = -p(:,1)
      Karr(:,3,2) = -p(:,3); Karr(:,4,2) = p(:,2); Karr(:,2,3) = p(:,3); Karr(:,2,4) = -p(:,2)
      Karr(:,4,3) = -p(:,1); Karr(:,3,4) = p(:,1)
    end subroutine opacity_matrix

      subroutine invert_delo_matrix_thin_single(dx,K0,ki,delta,identity,matrix,imatrix)
        double precision, intent(in) :: dx,ki,delta
        double precision, intent(in), dimension(4,4) :: identity,K0
        double precision, intent(out), dimension(4,4) :: matrix,imatrix
        matrix = (1d0-delta/2d0+delta**2d0/6d0)*identity & 
             +(0.5d0*dx-1d0/6d0*dx**2d0*ki)*K0
        call imatrix_4(matrix,imatrix)
      end subroutine invert_delo_matrix_thin_single

      subroutine calc_delo_P_thin_single(imatrix,dx,j,j1,ki,ki1,P)
        double precision, intent(in), dimension(4,4) :: imatrix
        double precision, intent(in) :: dx,ki,ki1
        double precision, intent(in), dimension(4) :: j,j1
        double precision, intent(out), dimension(4) :: P
        P = matmul(imatrix,(0.5d0*dx*j-1d0/6d0*dx**2d0*ki*j)+ &
             (0.5d0*j1*dx-1d0/3d0*dx**2d0*ki*j1))
      end subroutine calc_delo_P_thin_single

      subroutine calc_delo_Q_thin_single(imatrix,dx,ki,ki1,K1,identity,Q)
        double precision, intent(in), dimension(4,4) :: imatrix,K1,identity
        double precision, intent(in) :: dx,ki,ki1
        double precision, intent(out), dimension(4,4) :: Q
        Q = matmul(imatrix,identity*(1d0-0.5d0*dx*ki+1d0/6d0*dx**2d0*ki**2d0) &
             -(0.5d0*dx-1d0/3d0*dx**2d0)*K1)
      end subroutine calc_delo_Q_thin_single

      subroutine invert_delo_matrix_single(F,G,Kp,identity,matrix,imatrix)
        double precision, intent(in) :: F,G
        double precision, intent(in), dimension(4,4) :: Kp,identity
        double precision, intent(out), dimension(4,4) :: matrix,imatrix
        matrix = identity+(F-G)*Kp
        call imatrix_4(matrix,imatrix)
      end subroutine invert_delo_matrix_single

      subroutine calc_delo_P_single(imatrix,F,G,Sp,Sp1,P)
        double precision, dimension(4,4), intent(in) :: imatrix
        double precision, intent(in) :: F,G
        double precision, intent(in), dimension(4) :: Sp,Sp1
        double precision, intent(out), dimension(4) :: P
        P = matmul(imatrix,(F-G)*Sp+G*Sp1)
      end subroutine calc_delo_P_single

      subroutine calc_delo_Q_single(imatrix,E,F,G,Kp1,identity,Q)
        double precision, dimension(4,4), intent(in) :: imatrix,Kp1,identity
        double precision, intent(in) :: E,F,G
        double precision, dimension(4,4), intent(out) :: Q
        Q = matmul(imatrix,identity*E-G*Kp1)
      end subroutine calc_delo_Q_single

      subroutine radtrans_integrate_delo(x,tau,j,a,rho,P,Q,im)
        double precision, dimension(:), intent(in) :: x,tau
        double precision, dimension(:,:), intent(in) :: j,a,rho
!        double precision, intent(out), dimension(size(j,1),size(j,2)) :: P
!        double precision, intent(out), dimension(size(x),4,4) :: Q,im
!        double precision, intent(out), dimension(size(x)) :: delta
        double precision, dimension(size(x),4), intent(inout) :: P
        double precision, dimension(size(x),4,4), intent(inout) :: Q,im
        double precision, dimension(size(x)-1) :: delta,dx
        double precision :: E,F,G
        double precision, dimension(4) :: iprev,Sp,Sp1,pt,I0
        double precision, dimension(4,4) :: identity,Kp,Kp1,K0,K1,qt,imatrix,matrix
        double precision, dimension(size(x),4,4) :: Karr
        integer :: k
        I0=0d0
        identity = reshape((/1d0,0d0,0d0,0d0,0d0,1d0,0d0,0d0,0d0,0d0, &
             1d0,0d0,0d0,0d0,0d0,1d0/),(/4,4/))
        delta = tau(2:npts) - tau(1:npts-1)
!        E = np.exp(-delta)
!        F = 1.-E
!        G = (1.-(1.+delta)*E)/delta
        call opacity_matrix(a,rho,Karr)
        dx = x(1:npts-1) - x(2:npts)
!        dx = x(2:npts) - x(1:npts-1)
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
!              write(6,*) 'invert delo matrix call: ',k,dx(k),K0,a(k,1),delta(k)
!              write(6,*) 'invert delo matrix call identity: ',identity
              call invert_delo_matrix_thin(dx(k),K0,a(k,1),delta(k),identity,matrix,imatrix)
!              write(6,*) 'invert delo matrix after call matrix: ',matrix
!              write(6,*) 'invert delo matrix after call imatrix: ',imatrix
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

      subroutine del_radtrans_integrate_data()
        deallocate(jj); deallocate(KK)
        deallocate(s); deallocate(ss); deallocate(s0); deallocate(tau)
        deallocate(intensity)
        if(iflag.eq.1)then
           deallocate(PP); deallocate(QQ); deallocate(imm)
        endif
      end subroutine del_radtrans_integrate_data
   
    end module radtrans_integrate
