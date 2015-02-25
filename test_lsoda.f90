
      module common_data
      real(kind=8) :: nu=2.3d11
      real(kind=8), dimension(454) :: s,ncgs,tcgs,
     &  bcgs,ang
      real(kind=8) :: fac=1d15
      end module

      program lsoda_radtrans

      use odepack

      implicit none
      external dIdspol, jacpol, dIds, jac
      integer :: n=454, k
      real(kind=8), dimension(1,454) :: Iout, Iexact, tau
      real(kind=8), dimension(1) :: I
      real(kind=8), dimension(4,454) :: Ioutpol
      real(kind=8), dimension(4) :: Ipol
      real(kind=8), dimension(454) :: sout
      real(kind=8) :: oatol,ortol,temp1,temp2,temp3,temp4,temp5
      integer, dimension(3) :: stats
      type (fluid) :: vars
      type (emis) :: e
      character(len=40) :: FMT1="(F15.12,E15.8,E15.8,F15.12,F15.12)"
      ortol=1d-7; oatol=1d-7
      Ipol=0d0; I=0d0
!      allocate(s(n)); allocate(ncgs(n)); allocate(tcgs(n))
!      allocate(bcgs(n)); allocate(ang(n))
      open(unit=8,file='mb09m87_lsoda.inp',form='unformatted')
!      do k=1,n
        read(8) s
!        write(6,*) 's: ',s
        read(8) ncgs
!        write(6,*) 'ncgs: ',ncgs
        read(8) tcgs
!        write(6,*) 'tcgs: ',tcgs
        read(8) bcgs
!        write(6,*) 'bcgs: ', bcgs
        read(8) ang
  !    enddo
!      vars%ncgs=ncgs; vars%bcgs=bcgs; vars%tcgs=tcgs; vars%incang=ang
      ncgs=ncgs/83893.338d0*10d0
      write(6,*) maxval(ncgs), minval(tcgs), maxval(bcgs), minval(ang)
      write(6,*) maxval(s), size(s)
! First do non-pol:
      call lsoda_basic(dIds,I,(/s(1),s(454)/),
     & oatol,ortol,jac,Iout,1,100000,stats)
      Iexact(1,:)=(1d0-exp(-s*1d0))/1d0
!      write(6,*) 'non-pol: ',(Iout-Iexact)/Iexact
!      write(6,*) 'exact: ',Iexact
!      write(6,*) 't: ',tsave-s(25)
      write(6,*) 'unpol stats'
      write(6,*) 'N_f: ',stats(1)
      write(6,*) 'N_jac: ',stats(2)
      write(6,*) 'Method: ',stats(3)
! Now do pol:
!      do k=1,22500
      call lsoda_basic(dIdspol,Ipol,(/s(1),s(454)/)
     & ,oatol,
     & ortol,jacpol,Ioutpol,1,100000,stats)
 !     enddo
      write(6,*) 'pol stats'
      write(6,*) 'N_f: ',stats(1)
      write(6,*) 'N_jac: ',stats(2)
      write(6,*) 'Method: ',stats(3)
      write(6,*) 'pol: ',Ioutpol(:,2)
      write(6,*) 'nonpol: ',Iout(:,2)
!      deallocate(ang); deallocate(bcgs); deallocate(ncgs)
!      deallocate(s); deallocate(tcgs)
      close(unit=8)
      end
!      contains

      subroutine calc_coefs(s0,j,K)
      use fluid_model
      use emissivity
      use common_data
      real(kind=8), intent(in) :: s0
      real(kind=8), intent(out), dimension(4) :: j
      real(kind=8), intent(out), dimension(7) :: K
      integer, save :: lindx=1
      type (fluid) :: f
      type (emis) :: e
!      write(6,*) 'calc_coefs'
!      j=1d0
!      K(1)=1000.d0; K(2)=400.0d0; K(3)=30.3d0
!      K(4)=0.5d0; K(5)=0.5d0; K(6)=0.35d0; K(7)=-0.2d0
!      K(5:6)=0.01d0; K(7)=0.05d0
!      K=0d0
      call get_f(s0,f,lindx)
      call calc_emissivity(nu,f,e)
      j=e%j
      K=e%K
      end subroutine calc_coefs
 
      subroutine get_f(s0,f,lindx)
      use fluid_model
      use common_data
      use interpolate
      real(kind=8), intent(in) :: s0
      real(kind=8) :: weight
      type (fluid), intent(out) :: f
      integer, intent(inout) :: lindx
!      write(6,*) 'get_f: ',s0,size(s)
      call locate(s,s0,lindx)
!      write(6,*) 'hunt: ',lindx
      weight=(s0-s(lindx))/(s(lindx+1)-s(lindx))
      f%ncgs=(1d0-weight)*ncgs(lindx+1)+weight*ncgs(lindx)
      f%bcgs=(1d0-weight)*bcgs(lindx+1)+weight*bcgs(lindx)
      f%tcgs=(1d0-weight)*tcgs(lindx+1)+weight*tcgs(lindx)
      f%incang=(1d0-weight)*ang(lindx+1)+weight*ang(lindx)
!      write(6,*) 'f: ',f%ncgs,f%incang
      end subroutine get_f

!      end

      subroutine dIds(neq,s,I,Idot)
      integer, intent(in) :: neq
      real(kind=8), intent(in) :: s, I
      real(kind=8), intent(out) :: Idot
      real(kind=8), dimension(4) :: j
      real(kind=8), dimension(7) :: K 
      call calc_coefs(s,j,K)
!      K(5:7)=0d0
      Idot=j(1)*1d15-K(1)*I*9.41d14
      end subroutine dIds
 
      subroutine jac(neq,s,I,ml,mu,pd,nrowpd)
      integer, intent(in) :: neq, ml, mu, nrowpd
      real(kind=8), intent(in) :: I, s
      real(kind=8), intent(out), dimension(nrowpd,neq) :: pd
      real(kind=8), dimension(4) :: j
      real(kind=8), dimension(7) :: K
 !     write(6,*) 'jac'
      call calc_coefs(s,j,K)
!      K(5:7)=0d0
      pd(1,1)=-K(1)
      end subroutine jac
 
      subroutine dIdspol(neq,s,I,Idot)
      integer, intent(in) :: neq
      real(kind=8), intent(in) :: s
      real(kind=8), dimension(neq), intent(in) :: I
      real(kind=8), dimension(neq), intent(out) :: Idot
      real(kind=8), dimension(4) :: j
      real(kind=8), dimension(7) :: K
      call calc_coefs(s,j,K)
!      K(5:7)=0d0
!      write(6,*) 'dids'
      Idot(1)=j(1)*1d15-(K(1)*I(1)+K(2)*I(2)+K(3)*I(3)+K(4)*I(4))*
     &9.41d14
      Idot(2)=j(2)*1d15-(K(2)*I(1)+K(1)*I(2)+K(7)*I(3)-K(6)*I(4))*
     &9.41d14
      Idot(3)=j(3)*1d15-(K(3)*I(1)-K(7)*I(2)+K(1)*I(3)+K(5)*I(4))*
     &9.41d14
      Idot(4)=j(4)*1d15-(K(4)*I(1)+K(6)*I(2)-K(5)*I(3)+K(1)*I(4))*
     &9.41d14
      end subroutine dIdspol

      subroutine jacpol(neq,s,I,ml,mu,pd,nrowpd)
      integer, intent(in) :: neq, ml, mu, nrowpd
      real(kind=8), intent(in) :: s
      real(kind=8), intent(in), dimension(neq) :: I
      real(kind=8), intent(out), dimension(nrowpd,neq) :: pd
      real(kind=8), dimension(4) :: j
      real(kind=8), dimension(7) :: K
      call calc_coefs(s,j,K)   
 !     K(5:7)=0d0  
 !     write(6,*) 'jac' 
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
!      write(6,*) 'pd: ',nrowpd,neq,pd
      end subroutine jacpol
