      module kerr
      
      use class_four_Vector
      use math, only: dot_product
      use phys_constants, only: msun, G, c, sigt, pi, mp

      implicit none

      interface surf_integral
         module procedure surf_integral
         module procedure th_integral
         module procedure ph_integral
      end interface

      interface kerr_metric
        module procedure blmetric_cov
        module procedure blmetric_con
        module procedure blmetric_cov_real
        module procedure blmetric_con_real
        module procedure ksmetric_cov
!        module procedure ksmetric_con
!        module procedure kstmetric_sph
!        module procedure kstmetric_car
      end interface

      interface lnrf_frame
        module procedure lnrf_frame
        module procedure lnrf_frame_inv
        module procedure lnrf_frame_real
        module procedure lnrf_frame_inv_real
      end interface

      interface bl2ks
         module procedure bl2ks_phi
         module procedure bl2ks_time
         module procedure bl2ks_phi_single
         module procedure bl2ks_time_single
      end interface

      interface uks2ubl
         module procedure uks2ubl
      end interface

!      interface calc_kb_ang
!         module procedure calc_kb_ang
!      end interface
      
      interface calcg
         module procedure calcg
      end interface

      interface comoving_ortho
        module procedure comoving_ortho
      end interface

      interface calc_rms
        module procedure calc_rms
      end interface

      interface krolikc
        module procedure krolikc
      end interface

      interface ledd
         module procedure ledd
      end interface

      interface calc_polar_psi
         module procedure calc_polar_psi
      end interface

      interface calc_polvec
         module procedure calc_polvec
      end interface
         
      interface calc_kappapw
         module procedure calc_kappapw
      end interface

      interface calc_u0
         module procedure calc_u0
      end interface

      interface rms_vel
         module procedure rms_vel
      end interface

      interface calc_plunging_vel
         module procedure calc_plunging_vel
      end interface calc_plunging_vel

      contains

        function ledd(m) result(leddval)
        real(kind=8) :: leddval
        real, intent(in) :: m
! Eddington luminosity
        leddval=4.*pi*G*m*msun*mp*c/sigt
        end function ledd

        function calc_rms(a) result(rms)
        real, intent(in) :: a
        real :: rms,z1,z2
        z1=1.+(1.-a*a)**(1./3.)*((1.+a)**(1./3.)+(1.-a)**(1./3.))
        z2=sqrt(3.*a*a+z1*z1)
        rms=3.+z2-sign(sqrt((3.-z1)*(3.+z1+2.*z2)),a)
        end function calc_rms
 
        function krolikc(r,a) result(kc)
        real, intent(in) :: a
        real :: rms, pi
        real, dimension(:), intent(in) :: r
        real, dimension(size(r)) :: kc,y,yms,y1,y2,y3,arg1,arg2,arg3,arg4
        pi=acos(-1.)
        rms=calc_rms(a)
        y=sqrt(r)
        yms=sqrt(rms)
        y1=2d0*cos(1./3.*(acos(a)-pi))
        y2=2d0*cos(1./3.*(acos(a)+pi))
        y3=-2d0*cos(1./3.*acos(a))
        arg1=3.*a/(2d0*y)
        arg2=3.*(y1-a)**2d0/(y*y1*(y1-y2)*(y1-y3))
        arg3=3.*(y2-a)**2d0/(y*y2*(y2-y1)*(y2-y3))
        arg4=3.*(y3-a)**2d0/(y*y3*(y3-y1)*(y3-y2))
        kc=1.-yms/y-arg1*log(y/yms)-arg2*log((y-y1)/(yms-y1))- &
          arg3*log((y-y2)/(yms-y2))-arg4*log((y-y3)/(yms-y3))
!        write(6,*) 'krolikc: ',kc,arg1,arg2,arg3,arg4,y1,y2,y3
!        write(6,*) 'krolikc: ',a,yms,(y-y2)/(yms-y2), (y-y1)/(yms-y1),(y-y3)/(yms-y3)
        end function krolikc

        function uks2ubl(fut,r,a) result(fu)
        ! Converts a Kerr-Schild spherical 4-velocity to a Boyer-Lindquist one.
        ! JAD 11/10/2008
        type (four_vector), dimension(:), intent(in) :: fut
        real(kind=8), dimension(:), intent(in) :: r
        real(kind=8), intent(in) :: a
        real(kind=8), dimension(size(r)) :: delta
        type (four_vector), dimension(size(fut)) :: fu
        delta=r**2d0-2d0*r+a**2
        fu=fut
! from Font+1999
        fu%data(1)=fut%data(1)-2d0*r/delta*fut%data(2)
        fu%data(4)=fut%data(4)-a/delta*fut%data(2)
        end function uks2ubl


        function bl2ks_time(r,x,a,time) result(xtilde)
        real(kind=8), dimension(:), intent(in) :: r,x
        real(kind=8), intent(in) :: a
        real(kind=8), dimension(size(x)) :: xtilde
        integer, intent(in) :: time
        ! Function to convert t/phi from Boyer-Lindquist coordinates to Kerr-Schild coordinates as shown in PF98 and FIP99
        xtilde=x+log(r**2-2*r+a**2)+1./(2d0*sqrt(1.-a**2))*log((r-1.-sqrt(1.-a**2))/(r-1+sqrt(1.-a**2))) 
        end function bl2ks_time

        function bl2ks_phi(r,x,a) result(xtilde)
        real(kind=8), dimension(:), intent(in) :: r,x
        real(kind=8), intent(in) :: a
        real(kind=8), dimension(size(x)) :: xtilde
        ! Function to convert t/phi from Boyer-Lindquist coordinates to Kerr-Schild coordinates as shown in PF98 and FIP99
        xtilde=x+a/(2d0*sqrt(1.-a**2))*log((r-1.-sqrt(1.-a**2))/(r-1.+sqrt(1.-a**2)))
        end function bl2ks_phi

        function bl2ks_time_single(r,x,a,time) result(xtilde)
        real(kind=8), intent(in) :: r,x
        real(kind=8), intent(in) :: a
        real(kind=8) :: xtilde
        integer, intent(in) :: time
        ! Function to convert t/phi from Boyer-Lindquist coordinates to Kerr-Schild coordinates as shown in PF98 and FIP99
        xtilde=x+log(r**2-2*r+a**2)+1./(2d0*sqrt(1.-a**2))*log((r-1.-sqrt(1.-a**2))/(r-1+sqrt(1.-a**2))) 
        end function bl2ks_time_single

        function bl2ks_phi_single(r,x,a) result(xtilde)
        real(kind=8), intent(in) :: r,x
        real(kind=8), intent(in) :: a
        real(kind=8) :: xtilde
        ! Function to convert t/phi from Boyer-Lindquist coordinates to Kerr-Schild coordinates as shown in PF98 and FIP99
        xtilde=x+a/(2d0*sqrt(1.-a**2))*log((r-1.-sqrt(1.-a**2))/(r-1.+sqrt(1.-a**2)))
        end function bl2ks_phi_single

        FUNCTION CALCG(U,MU,Q2,L,A,TPM,TPR,SU,SM,VRL,VTL,VPL) result(g)
          real(kind=8), dimension(:), intent(in) :: su,sm,u,mu,tpm,tpr,vrl,vtl,vpl,q2,l
          real(kind=8), intent(in) :: a
          real(kind=8), dimension(size(u)) :: zero,r,z1,z2,rms,d,ar,rho,enu, &
               emu1,emu2,epsi,om,sr,st,omega,gam,rr,tt,g
          real(kind=8) :: one,two,tres,third
          ! see 2/25/99 notes, p.17
          zero=0.d0
          ONE=1.d0
          TWO=2D0
          TRES=3.d0
          THIRD=ONE/TRES
          R=ONE/U
          Z1=ONE+(ONE-A*A)**THIRD*((ONE+A)**THIRD+(ONE-A)**THIRD)
          Z2=SQRT(TRES*A*A+Z1*Z1)
          RMS=TRES+Z2-SIGN(1.d0,A)*SQRT((TRES-Z1)*(TRES+Z1+TWO*Z2))
          D=R*R-TWO*R+A*A
          AR=(R*R+A*A)**2-A*A*D*(ONE-MU*MU)
          RHO=R*R+A*A*MU*MU
          ENU=SQRT(D*RHO/AR)
          EMU1=SQRT(RHO/D)
          EMU2=SQRT(RHO)
          EPSI=SQRT(ONE-MU*MU)*SQRT(AR/RHO)
          OM=TWO*A*R/AR
          SR=(-ONE)**TPR*SU
          ST=-(-ONE)**TPM*SM
          ! Now, compute velocities in LNRF
          OMEGA=ENU/EPSI*VPL+OM
          omega=merge(omega,zero,epsi.ne.0d0)
          GAM=ONE/SQRT(ONE-(VRL*VRL+VTL*VTL+VPL*VPL))
! Now, compute R(r) and Theta(theta):
          RR=-A*A*Q2*U**4+TWO*U**3*(Q2+(A-L)**2)+U*U*(A*A-Q2-L*L)+ONE
          TT=(Q2+MU*MU*(A*A-L*L-Q2)-A*A*MU**4)/(ONE-MU*MU)
          tt=merge(sqrt(tt),zero,tt.ge.0d0)
          rr=merge(sqrt(rr)*r*r,zero,rr.ge.0d0)
! Now, we can compute the redshift:
          G=ENU/GAM/(ONE-L*OMEGA-EMU1*ENU*VRL/RHO*SR*RR-EMU2*ENU*VTL/RHO*ST*TT)
        END FUNCTION CALCG

        function calc_kb_ang(k,b,u,r,th,a) result(ang)
          ! Calculate angle between wave-vector and magnetic field in fluid                                
          ! frame based on Broderick (2004) covariant method                                      
          ! JAD 2/8/2011
          type (four_vector), dimension(:), intent(in) :: k,b,u
          real(kind=8), dimension(:), intent(in) :: r,th
          real(kind=8), intent(in) :: a
          real(kind=8), dimension(size(r)) :: bdotk,bdotb,kdotk,om,angnorm,ang,zero,one, &
               angmin,angmax, cdot2
          real(kind=8), dimension(size(r),10) :: metric
          real(kind=8), dimension(10,size(r)) :: tmetric
          zero=0.d0; one=1.d0; angmin=-.9999d0; angmax=0.9999d0
          metric=kerr_metric(r,th,a); tmetric=transpose(metric)
          call assign_metric(b,tmetric)
          call assign_metric(u,tmetric)
          call assign_metric(k,tmetric)
!          write(6,*) 'metric: ',b(1)%metric, u(1)%metric,k(1)%metric
!          write(6,*) 'vectors: ',b(1)%data,u(1)%data,k(1)%data
          bdotk=b*k
          bdotb=b*b
          kdotk=k*k
          om=-k*u
!          angnorm=bdotk/sqrt(om*om)/sqrt(bdotb)
!          ang=acos(merge(merge(angnorm,angmax,angnorm.le.0.9999d0),angmin, &
!               angnorm.ge.-.9999d0))
        
          cdot2=bdotk*bdotk/bdotb/(kdotk+om*om)
!          write(6,*) 'kb: ',r,th,a
!          write(6,*) 'kb: ',k(1)%data, b(1)%data, u(1)%data
!          write(6,*) 'cdot2: ',cdot2,bdotk,bdotb,kdotk,om*om
          cdot2=merge(merge(cdot2,zero,cdot2.ge.0.d0),one,cdot2.le.1.d0)
          ang=acos(sign(1d0,bdotk)*sqrt(cdot2))
        end function calc_kb_ang


        function calc_nullp(q2,l,a,r,mu,su,smu,rcomp,thcomp) result(k)
        real(kind=8), intent(in) :: q2,l,a
        real(kind=8), dimension(:), intent(in) :: r,mu,su,smu
        integer, intent(in), optional :: rcomp,thcomp
        real(kind=8), dimension(size(r)) :: zero,Mfunc,Ufunc, &
          rho2,delta,u,pu,pt,pmu,pphi,Rfunc
        type (four_Vector), dimension(size(r)) :: k
        zero=0d0
        u=1d0/r
        rho2=r*r+a*a*mu*mu
        delta=r*r-2.d0*r+a*a
        Mfunc=q2+(a*a-q2-l*l)*mu*mu-a*a*mu*mu*mu*mu
        Mfunc=merge(Mfunc,zero,Mfunc.gt.0d0)
        if(present(thcomp)) then
! This is negative of RB94 to account for dlambda -> -dlambda for forward in time.
! see 12/11/12 notes
          pmu=1d0/rho2*smu*sqrt(Mfunc/(1.d0-mu*mu))
!          write(6,*) 'thcomp'
        else
          pmu=-smu*sqrt(Mfunc)/rho2
        endif
        Ufunc=1.d0+(a*a-q2-l*l)*u*u+2.d0*((a-l)**2+q2)*u*u*u- &
              a*a*q2*u*u*u*u
        Ufunc=merge(Ufunc,zero,Ufunc.gt.0d0)
        if(present(rcomp)) then
          Rfunc=R*R*sqrt(Ufunc)
          pu=su*Rfunc/rho2
!          write(6,*) 'rcomp'
        else
          pu=-su*sqrt(Ufunc)/rho2
        endif
        pt=(-a*(a*(1.d0-mu*mu)-l)+(r*r+a*a)/delta*(r*r+a*a-a*l))/rho2
        pphi=(-a+l/(1.d0-mu*mu)+a/delta*(r*r+a*a-a*l))/rho2
        k%data(1)=pt; k%data(2)=pu; k%data(3)=pmu; k%data(4)=pphi
!        write(6,*) 'nullp: ',pt,pu,pmu,pphi
        end function calc_nullp

        function blmetric_con_real(r,th,a,con) result(blmetric)
        ! Returns boyer-lindquist covariant metric components
        ! JAD 10/1/2008, fortran 3/28/2011
        real, dimension(:), intent(in) :: r,th
        real, intent(in) :: a
        integer, intent(in) :: con
        real, dimension(size(r),10) :: blmetric
!        real, dimension(10,size(r)) :: gmunu
        real, dimension(size(r)) :: cth,sth,delta, &
        rho2,sigma
        cth=cos(th); sth=sin(th)
        delta=r*r-2.*r+a*a
        rho2=r*r+a*a*cth*cth
        sigma=(r*r+a*a)**2-a*a*delta*sth*sth
        blmetric=0d0
        blmetric(:,1)=-((r*r+a*a)**2-a*a*delta*sth*sth)/rho2/delta
        blmetric(:,4)=-2*a*r/rho2/delta
        blmetric(:,5)=delta/rho2
        blmetric(:,8)=1./rho2
        blmetric(:,10)=(delta-a*a*sth*sth)/delta/rho2/sth/sth
!        gmunu=transpose(blmetric)
        end function blmetric_con_real

        function ksmetric_cov(r,th,ph,a) result(ksmetric)
        ! Returns boyer-lindquist covariant metric components
        ! JAD 10/1/2008, fortran 3/28/2011
        ! Doesn't actually depend on ph unless you're doing tilted metric!
        real(8), dimension(:), intent(in) :: r,th,ph
        real(8), intent(in) :: a
        real(8), dimension(size(r),10) :: ksmetric
        real(8), dimension(size(r)) :: ctheta,stheta,rho2,psi4
        ctheta=cos(th); stheta=sin(th)
        rho2=r**2+a**2*ctheta**2
        psi4=2d0*r/rho2
        ksmetric=0.
        ksmetric(:,1)=-(1.-psi4)
        ksmetric(:,2)=psi4
        ksmetric(:,4)=-a*stheta**2*psi4
        ksmetric(:,5)=1.+psi4
        ksmetric(:,7)=-a*stheta**2*(1.+psi4)
        ksmetric(:,8)=rho2
        ksmetric(:,10)=stheta**2*(rho2+a**2*(1.+2d0*r/rho2)*stheta**2)
!        write(6,*) 'ksmetric kerr: ',ksmetric(1,:)
      end function ksmetric_cov

        function blmetric_con(r,th,a,con) result(blmetric)
        ! Returns boyer-lindquist covariant metric components
        ! JAD 10/1/2008, fortran 3/28/2011
        real(kind=8), dimension(:), intent(in) :: r,th
        real(kind=8), intent(in) :: a
        integer, intent(in) :: con
        real(kind=8), dimension(size(r),10) :: blmetric
!        real(kind=8), dimension(10,size(r)) :: gmunu
        real(kind=8), dimension(size(r)) :: cth,sth,delta, &
        rho2,sigma
        cth=cos(th); sth=sin(th)
        delta=r*r-2d0*r+a*a
        rho2=r*r+a*a*cth*cth
        sigma=(r*r+a*a)**2-a*a*delta*sth*sth
        blmetric=0d0
        blmetric(:,1)=-((r*r+a*a)**2-a*a*delta*sth*sth)/rho2/delta
        blmetric(:,4)=-2d0*a*r/rho2/delta
        blmetric(:,5)=delta/rho2
        blmetric(:,8)=1d0/rho2
        blmetric(:,10)=(delta-a*a*sth*sth)/delta/rho2/sth/sth
!        gmunu=transpose(blmetric)
        end function blmetric_con

        function blmetric_cov_real(r,th,a) result(blmetric)
        ! Returns boyer-lindquist covariant metric components
        ! JAD 10/1/2008, fortran 3/28/2011
        real, dimension(:), intent(in) :: r,th
        real, intent(in) :: a
        real, dimension(size(r),10) :: blmetric
!        real, dimension(10,size(r)) :: gmunu
        real, dimension(size(r)) :: cth,sth,delta, &
        rho2,sigma
        cth=cos(th); sth=sin(th)
        delta=r*r-2.*r+a*a
        rho2=r*r+a*a*cth*cth
        sigma=(r*r+a*a)**2-a*a*delta*sth*sth
        blmetric=0d0
        blmetric(:,1)=-(delta-a*a*sth*sth)/rho2
        blmetric(:,4)=-2.*a*r*sth*sth/rho2
        blmetric(:,5)=rho2/delta
        blmetric(:,8)=rho2
        blmetric(:,10)=sigma/rho2*sth*sth
        end function blmetric_cov_real

        function blmetric_cov(r,th,a) result(blmetric)
        ! Returns boyer-lindquist covariant metric components
        ! JAD 10/1/2008, fortran 3/28/2011
        real(kind=8), dimension(:), intent(in) :: r,th
        real(kind=8), intent(in) :: a
        real(kind=8), dimension(size(r),10) :: blmetric
!        real(kind=8), dimension(10,size(r)) :: gmunu
        real(kind=8), dimension(size(r)) :: cth,sth,delta, &
        rho2,sigma
        cth=cos(th); sth=sin(th)
        delta=r*r-2d0*r+a*a
        rho2=r*r+a*a*cth*cth
        sigma=(r*r+a*a)**2-a*a*delta*sth*sth
        blmetric=0d0
        blmetric(:,1)=-(delta-a*a*sth*sth)/rho2
        blmetric(:,4)=-2d0*a*r*sth*sth/rho2
        blmetric(:,5)=rho2/delta
        blmetric(:,8)=rho2
        blmetric(:,10)=sigma/rho2*sth*sth
        end function blmetric_cov

        subroutine lnrf_frame(vr,vt,omega,r,a,th,vrl,vtl,vpl)
        real(kind=8), dimension(:), intent(in) :: vr,vt,omega,r,th
        real(kind=8), intent(in) :: a
        real(kind=8), dimension(size(r)), intent(out) :: vrl,vtl,vpl
        real(kind=8), dimension(size(r)) :: zero,d,ar,rho,enu,emu1, &
         emu2,epsi,om,mu
        zero=0d0
        mu=cos(th)
        D=R*R-2d0*R+A*A
        AR=(R*R+A*A)**2-A*A*D*(1.-MU*MU)
        RHO=R*R+A*A*MU*MU
        ENU=SQRT(D*RHO/AR)
        EMU1=SQRT(RHO/D)
        EMU2=SQRT(RHO)
        EPSI=SQRT(1.-MU*MU)*SQRT(AR/RHO)
        OM=2d0*A*R/AR
        VRL=EMU1/ENU*VR
        VTL=EMU2/ENU*VT
        VPL=EPSI/ENU*(OMEGA-OM)
        VRL=MERGE(VRL,ZERO,D.GT.0d0)
        VTL=MERGE(VTL,ZERO,D.GT.0d0)
!       VPL=MERGE(VPL,ZERO,(D.GT.0d0.AND.ABS(MU).LT.1d0))
        VPL=MERGE(VPL,ZERO,D.GT.0d0)
        end subroutine lnrf_frame

        subroutine lnrf_frame_real(vr,vt,omega,r,a,th,vrl,vtl,vpl)
        real, dimension(:), intent(in) :: vr,vt,omega,r,th
        real, intent(in) :: a
        real, dimension(size(r)), intent(out) :: vrl,vtl,vpl
        real, dimension(size(r)) :: zero,d,ar,rho,enu,emu1, &
         emu2,epsi,om,mu
        zero=0d0
        mu=cos(th)
        D=R*R-2.*R+A*A
        AR=(R*R+A*A)**2-A*A*D*(1.-MU*MU)
        RHO=R*R+A*A*MU*MU
        ENU=SQRT(D*RHO/AR)
        EMU1=SQRT(RHO/D)
        EMU2=SQRT(RHO)
        EPSI=SQRT(1.-MU*MU)*SQRT(AR/RHO)
        OM=2.*A*R/AR
        VRL=EMU1/ENU*VR
        VTL=EMU2/ENU*VT
        VPL=EPSI/ENU*(OMEGA-OM)
        VRL=MERGE(VRL,ZERO,D.GT.0d0)
        VTL=MERGE(VTL,ZERO,D.GT.0d0)
        VPL=MERGE(VPL,ZERO,D.GT.0d0)
        end subroutine lnrf_frame_real

        subroutine lnrf_frame_inv(vr,vt,omega,r,a,th,vrl,vtl,vpl,inv)
        real(kind=8), dimension(:), intent(in) :: vr,vt,omega,r,th
        real(kind=8), intent(in) :: a
        integer, intent(in) :: inv
        real(kind=8), dimension(size(r)), intent(out) :: vrl,vtl,vpl
        real(kind=8), dimension(size(r)) :: zero,d,ar,rho,enu,emu1, &
        emu2,epsi,om,mu
        zero=0d0
        mu=cos(th)
        D=R*R-2.*R+A*A
        AR=(R*R+A*A)**2-A*A*D*(1.-MU*MU)
        RHO=R*R+A*A*MU*MU
        ENU=SQRT(D*RHO/AR)
        EMU1=SQRT(RHO/D)
        EMU2=SQRT(RHO)
        EPSI=SQRT(1.-MU*MU)*SQRT(AR/RHO)
        OM=2.*A*R/AR
        VRL=ENU/EMU1*VR
        VTL=ENU/EMU2*VT
        VPL=ENU/EPSI*OMEGA+OM
        VRL=MERGE(VRL,ZERO,D.GT.0d0)
        VTL=MERGE(VTL,ZERO,D.GT.0d0)
        VPL=MERGE(VPL,ZERO,D.GT.0d0)
        end subroutine lnrf_frame_inv

        subroutine lnrf_frame_inv_real(vr,vt,omega,r,a,th,vrl,vtl, &
         vpl,inv)
        real, dimension(:), intent(in) :: vr,vt,omega,r,th
        real, intent(in) :: a
        integer, intent(in) :: inv
        real, dimension(size(r)), intent(out) :: vrl,vtl,vpl
        real, dimension(size(r)) :: zero,d,ar,rho,enu,emu1, &
        emu2,epsi,om,mu
        zero=0d0
        mu=cos(th)
        D=R*R-2.*R+A*A
        AR=(R*R+A*A)**2-A*A*D*(1.-MU*MU)
        RHO=R*R+A*A*MU*MU
        ENU=SQRT(D*RHO/AR)
        EMU1=SQRT(RHO/D)
        EMU2=SQRT(RHO)
        EPSI=SQRT(1.-MU*MU)*SQRT(AR/RHO)
        OM=2.*A*R/AR
        VRL=ENU/EMU1*VR
        VTL=ENU/EMU2*VT
        VPL=ENU/EPSI*OMEGA+OM
        VRL=MERGE(VRL,ZERO,D.GT.0d0)
        VTL=MERGE(VTL,ZERO,D.GT.0d0)
        VPL=MERGE(VPL,ZERO,D.GT.0d0)
        end subroutine lnrf_frame_inv_real

        subroutine transport_perpk(kvec,r,th,a,metric,Kap1,Kap2,f1, &
         f2,f3)
        ! Calculate components of parallel transported vector with
        ! Walker-Penrose constants Kap1, Kap2 at locations 
        !r,th in Boyer-Lindquist coordinates assuming f0=0.
        ! JAD 2/22/2011 from notes
        real(kind=8), intent(in), dimension(:) :: r,th
        type (four_Vector), dimension(:), intent(in) :: kvec
        real(kind=8), intent(in) :: a,Kap1,Kap2
        real(kind=8), intent(in), dimension(:,:) :: metric
        real(kind=8), intent(out), dimension(size(r)) :: f1,f2,f3
        real(kind=8), dimension(size(r)) :: g03,g11,g22,g33, &
         cth,sth,k0,k1,k2,k3,gam1,gam2,gam3,del1,del2,del3,denom
!  !      write(6,*) 'perpk', size(metric,1), size(metric,2)
        g03=metric(:,4); g11=metric(:,5)
        g22=metric(:,8); g33=metric(:,10)
!  !      write(6,*) 'perpk metric'
        cth=cos(th); sth=sin(th)
 ! !      write(6,*) 'perpk kk: ',kvec%data(1),
 !    &   kvec%data(2),kvec%data(3),
 !    &   kvec%data(4),kvec*kvec
        k0=kvec%data(1); k1=kvec%data(2)
        k2=kvec%data(3); k3=kvec%data(4)
!  !      write(6,*) 'perpk k'
        gam1=a*cth*k0-a*a*cth*sth*sth*k3
        gam2=r*(r*r+a*a)*sth*k3-a*r*sth*k0
        gam3=a*a*cth*sth*sth*k1-r*(r*r+a*a)*sth*k2
        del1=r*k0-r*a*sth*sth*k3
!        del2=a*cth*sth*(r*r+a*a)*k3-a*a*sth*cth*k0
! sign changed 8/24/2015 see notes! leads to small changes ~% usually
        del2=-a*cth*sth*(r*r+a*a)*k3+a*a*sth*cth*k0
        del3=r*a*sth*sth*k1+a*cth*sth*(r*r+a*a)*k2
! new version based on mathematica notebook 2/23/2015 trying to avoid denominator problems
! changed back to version in mathematica notebook now that i am using correct initial Kap1, Kap2 and not swapped
        denom=(gam2*del1-gam1*del2)*(g33*k3+g03*k0)+(gam3*del2-gam2*del3)*g11*k1-(gam3*del1-gam1*del3)*g22*k2
        f1=(gam2*Kap1-del2*Kap2)*(g33*k3+g03*k0)-g22*k2*(gam3*Kap1-del3*Kap2)
        f2=(del1*Kap2-gam1*Kap1)*(g33*k3+g03*k0)+g11*k1*(gam3*Kap1-del3*Kap2)
        f3=g22*k2*(gam1*Kap1-del1*Kap2)-g11*k1*(gam2*Kap1-del2*Kap2)
!        f1=(gam2*Kap2-del2*Kap1)*(g33*k3+g03*k0)-g22*k2*(gam3*Kap2-del3*Kap1)
!        f2=(del1*Kap1-gam1*Kap2)*(g33*k3+g03*k0)+g11*k1*(gam3*Kap2-del3*Kap1)
!        f3=g22*k2*(gam1*Kap2-del1*Kap1)-g11*k1*(gam2*Kap2-del2*Kap1)
        where(abs(denom).gt.0d0)
           f1=f1/denom
           f2=f2/denom
           f3=f3/denom
        endwhere
        end subroutine transport_perpk

        subroutine comoving_ortho(r,th,a,alpha,beta,mus,u,b,k, &
             s2xi,c2xi,ang,g,cosne)
        type (four_Vector), dimension(:), intent(in) :: u,b,k
        real(kind=8), dimension(:), intent(in) :: r,th
        real(kind=8), intent(in) :: a,alpha,beta,mus
        real(kind=8), intent(out), dimension(size(r)) :: s2xi,c2xi, &
             ang,g,cosne
        type (four_vector), dimension(size(r)) :: bhat,khat,aa
        real(kind=8), dimension(size(r),3) :: aahat
        call comoving_ortho_core(r,th,a,alpha,beta,mus,u,b,k,s2xi, &
             c2xi,ang,g,cosne,aa,aahat,bhat,khat)
        return
        end subroutine comoving_ortho

        subroutine comoving_ortho_debug(r,th,a,alpha,beta,mus,u,b,k, &
             s2xi,c2xi,ang,g,cosne,aahat,aat,aar,aath,aaph,kht,khr,khth,khph, &
             bht,bhr,bhth,bhph)
        type (four_Vector), dimension(:), intent(in) :: u,b,k
        real(kind=8), dimension(:), intent(in) :: r,th
        real(kind=8), intent(in) :: a,alpha,beta,mus
        real(kind=8), intent(out), dimension(size(r)) :: s2xi,c2xi, &
             ang,g,cosne
        type (four_vector), dimension(size(r)) :: bhat,khat,aa
        real(kind=8), dimension(size(r),3), intent(out) :: aahat
        real(kind=8), dimension(size(r)), intent(out) :: bht,bhr,bhth,bhph, &
             aat,aar,aath,aaph,kht,khr,khth,khph
        call comoving_ortho_core(r,th,a,alpha,beta,mus,u,b,k,s2xi, &
             c2xi,ang,g,cosne,aa,aahat,bhat,khat)
        bht=bhat%data(1); bhr=bhat%data(2); bhth=bhat%data(3)
        bhph=bhat%data(4); kht=khat%data(1); khr=khat%data(2)
        khth=khat%data(3); khph=khat%data(4)
        aat=aa%data(1); aar=aa%data(2); aath=aa%data(3)
        aaph=aa%data(4)
        return
      end subroutine comoving_ortho_debug

        subroutine comoving_ortho_core(r,th,a,alpha,beta,mus,u,b,k, &
         s2xi,c2xi,ang,g,cosne,aa,aahat,bhat,khat)
! Transform four-velocity u, magnetic field b and wave-vector k to the
! co-moving orthonormal frame following Beckwith et al (2008) and
! Shcherbakov & Huang (2011)
! JAD 2/14/2011, fortran 3/27/2011
!        implicit none
        type (four_Vector), dimension(:), intent(in) :: u,b,k
        real(kind=8), dimension(:), intent(in) :: r,th
        real(kind=8), intent(in) :: a,alpha,beta,mus
        real(kind=8), intent(out), dimension(size(r)) :: s2xi,c2xi, &
             ang,g,cosne
        type (four_Vector), dimension(size(r)) :: uhat!, &
!         bhat,khat
        type (four_vector), dimension(size(r)), intent(out) :: bhat,khat,aa
        real(kind=8), dimension(size(r),10) :: metric
        real(kind=8), dimension(10,size(r)) :: tmetric
        real(kind=8), dimension(size(r)) :: gtt,gtp,grr,gmm,gpp,ut, &
             ur,um,up,urc,umc,utc,upc,np2,nm2,nr2,uhatt,uhatr,uhatm, &
             bhatt,bhatr,bhatm,bhatp,khatt,khatr,khatm,khatp,aahatt, &
             aahatm,aahatr,aahatp,delta,al1,al2,al3,uhatp,ahatt,ahatm, &
             ahatr,ahatp,bdotb,bdotk,kdotk,om,om2,knorm,anorm,bpdotbp, &
             bpdotbb,aadotbp,sxi,cxi,eps,one,mone,angnorm,angmin,angmax, &
             cxitest,sxitest,xi,be1,be2
        real(kind=8), dimension(size(r),3) :: bbhat,bphat
!        real(kind=8), dimension(size(r),3) :: aahat &
!           ,bbhat
        real(kind=8), dimension(size(r),3), intent(out) :: aahat
        type (four_Vector), dimension(size(r)) :: &
             ekt,ekr,ekm,ekp,stb,srb,spb,smb,ahat
        integer :: i, nr, testindx
        real(kind=8) :: kap1,kap2
        one=1d0; mone=-one; angmin=-0.99d0; angmax=0.99d0
        eps=epsilon(eps)
        metric=kerr_metric(r,th,a); tmetric=transpose(metric)
! Co-variant metric components
        gtt=metric(:,1) ; gtp=metric(:,4) ; grr=metric(:,5)
        gmm=metric(:,8) ; gpp=metric(:,10)
! Gather all possible four-velocity components to simplify things:
        ut=u%data(1) ; ur=u%data(2) ; um=u%data(3) ; up=u%data(4)
        utc=gtt*ut+gtp*up; upc=gpp*up+gtp*ut ; urc=grr*ur ; umc=um*gmm
        nr=size(r)
! Find parallel transported basis e^\th, e^ph everywhere on
! ray from complex Walker-Penrose constant
! (see 2/22/2011 notes, Chandrasekhar 1983 and CS77, CPS80)
!        Kap1=beta
!        Kap2=-alpha-a*sqrt(1.-mus*mus)
! JAD 10/5/2015
! Kap1, Kap2 changed to be consistent with my math. accompanying change in transport_perpk leaves everything unchanged and makes it more clear.
        Kap1=alpha+a*sqrt(1.-mus*mus)
        Kap2=-beta
        if((Kap1.ne.0d0).or.(Kap2.ne.0d0)) then
          call transport_perpk(k,r,th,a,metric,Kap1,Kap2,al1,al2,al3)
        else
          al1=0d0; al2=0d0; al3=1d0/sqrt(metric(:,10))
        endif
        aa%data(1)=0d0;aa%data(2)=al1;aa%data(3)=al2;aa%data(4)=al3
        call assign_metric(aa,tmetric)
! Kulkarni+2011 normalizations:
        delta=r*r+a*a-2.*r
        nr2=-grr*(utc*ut+upc*up)*(1.+umc*um)
        nm2=gmm*(1.+umc*um)
        np2=-(utc*ut+upc*up)*delta*sin(th)*sin(th)
        ekt=(-1d0)*u
        ekr%data(1)=(urc*ut)/sqrt(nr2); ekr%data(2)=-(utc*ut+upc*up) &
            /sqrt(nr2)
        ekr%data(3)=0d0; ekr%data(4)=(urc*up)/sqrt(nr2)
        ekm%data(1)=(umc*ut)/sqrt(nm2); ekm%data(2)=(umc*ur)/sqrt(nm2)
        ekm%data(3)=(1+umc*um)/sqrt(nm2); ekm%data(4)=(umc*up)/sqrt(nm2)
        ekp%data(1)=upc/sqrt(np2); ekp%data(2)=0d0; ekp%data(3)=0d0
        ekp%data(4)=-utc/sqrt(np2)
        call assign_metric(ekt,tmetric); call assign_metric(ekm,tmetric)
        call assign_metric(ekp,tmetric); call assign_metric(ekr,tmetric)
        uhatt=ekt*u; uhatr=ekr*u; uhatm=ekm*u; uhatp=ekp*u
        bhatt=ekt*b; bhatr=ekr*b; bhatm=ekm*b; bhatp=ekp*b
        khatt=ekt*k; khatr=ekr*k; khatm=ekm*k; khatp=ekp*k
        ahatt=ekt*aa; ahatr=ekr*aa; ahatm=ekm*aa; ahatp=ekp*aa
        uhat%data(1)=uhatt; uhat%data(2)=uhatr; uhat%data(3)=uhatm; uhat%data(4)=uhatp
! now transform b, k:
        bhat%data(1)=bhatt; bhat%data(2)=bhatr; bhat%data(3)=bhatm; bhat%data(4)=bhatp
        khat%data(1)=khatt; khat%data(2)=khatr; khat%data(3)=khatm; khat%data(4)=khatp
        ahat%data(1)=ahatt; ahat%data(2)=ahatr; ahat%data(3)=ahatm; ahat%data(4)=ahatp
        call assign_metric(b,tmetric); call assign_metric(u,tmetric)
        kdotk=khatr*khatr+khatm*khatm+khatp*khatp
!        write(6,*) 'metric ortho: ',b(1)%metric,u(1)%metric,k(1)%metric
!        write(6,*) 'vector ortho: ',b(1)%data,u(1)%data,k(1)%data
        bdotk=b*k; bdotb=b*b
        om=-k*u; om2=om*om
        aahat(:,1)=ahatr-khat%data(2)*ahatt/khat%data(1)
        aahat(:,2)=ahatm-khat%data(3)*ahatt/khat%data(1)
        aahat(:,3)=ahatp-khat%data(4)*ahatt/khat%data(1)
!        write(6,*) 'comoving ortho aa 2: ',aa%data(2)
!        write(6,*) 'comoving ortho aa 3: ',aa%data(3)
!        write(6,*) 'comoving ortho aa 4: ',aa%data(4)
!        write(6,*) 'comoving ortho ahatt: ',ahatt
!        write(6,*) 'comoving ortho khatt: ',khat%data(1)
!        write(6,*) 'comoving ortho aahat 1: ',aahat(:,1)
!        write(6,*) 'comoving ortho aahat 2: ',aahat(:,2)
!        write(6,*) 'comoving ortho aahat 3: ',aahat(:,3)
!        anorm=dot_product(aahat,aahat)
!        write(6,*) 'anorm: ',anorm
!        aahat(:,1)=aahat(:,1)/sqrt(anorm)
!        aahat(:,2)=aahat(:,2)/sqrt(anorm)
!        aahat(:,3)=aahat(:,3)/sqrt(anorm)
        knorm=khatr*khatr+khatm*khatm+khatp*khatp
        bbhat(:,1)=-(aahat(:,2)*khatp-aahat(:,3)*khatm)/sqrt(knorm)
        bbhat(:,2)=-(aahat(:,3)*khatr-aahat(:,1)*khatp)/sqrt(knorm)
        bbhat(:,3)=-(aahat(:,1)*khatm-aahat(:,2)*khatr)/sqrt(knorm)
        bdotk=bhat%data(2)*khatr+bhat%data(3)*khatm+bhat%data(4)*khatp
!        bphat(:,1)=bhat%data(2)-khatr*bdotk/knorm
!        bphat(:,2)=bhat%data(3)-khatm*bdotk/knorm
!        bphat(:,3)=bhat%data(4)-khatp*bdotk/knorm
!        aadotbp=dot_product(aahat,bphat)
!        bpdotbp=dot_product(bphat,bphat)
!        bpdotbb=dot_product(bphat,bbhat)
!        write(6,*) 'test aadotbp: ',aadotbp/sqrt(bpdotbp)
!        write(6,*) 'test aadotbp 2: ',(aahat(:,1)*bhat%data(2)+aahat(:,2)*bhat%data(3)+bhat%data(4)*aahat(:,3))/ &
!             sqrt((aahat(:,1)*bhat%data(2)+aahat(:,2)*bhat%data(3)+aahat(:,3)*bhat%data(4))**2.+(bbhat(:,1) &
!             *bhat%data(2)+bbhat(:,2)*bhat%data(3)+bbhat(:,3)*bhat%data(4))**2.)
!        write(6,*) 'test aadotbp max: ',maxval(abs(aadotbp-aahat(:,1)*bhat%data(2)+aahat(:,2)*bhat%data(3)+bhat%data(4)*aahat(:,3)))
! Now compute angle between polarization basis and magnetic field:
! if field is zero, emissivity = 0 so arbitrary angles to prevent NaNs
        where(bdotb.gt.0d0)
           aadotbp=bhat%data(2)*aahat(:,1)+bhat%data(3)*aahat(:,2)+bhat%data(4)*aahat(:,3)
           bpdotbb=bhat%data(2)*bbhat(:,1)+bhat%data(3)*bbhat(:,2)+bhat%data(4)*bbhat(:,3)
           s2xi=-2d0*aadotbp*bpdotbb/(aadotbp**2d0+bpdotbb**2d0)
           c2xi=(bpdotbb*bpdotbb-aadotbp*aadotbp)/(aadotbp**2d0+bpdotbb**2d0)
           angnorm=bdotk/sqrt(knorm)/sqrt(bdotb)
        elsewhere
           s2xi=0d0
           c2xi=one
           angnorm=0.5d0
        endwhere
!        testindx=maxloc(bdotb,1)
!        write(6,*) 'bperp bp: ',bpdotbp(testindx)
!        write(6,*) 'bperp bdotb: ',bdotb(testindx)
!        write(6,*) 'bperp aadotbp: ',aadotbp(testindx)
!        write(6,*) 'bperp s2xi: ',s2xi(testindx)
!        write(6,*) 'bperp c2xi: ',c2xi(testindx)
!        write(6,*) 'unit vectors: ',dot_product(bbhat,bbhat),dot_product(aahat,aahat),sqrt(knorm)
!        write(6,*) 'orthogonal: ',dot_product(aahat,bbhat),(aahat(:,1)*khatr+aahat(:,2)*khatm+aahat(:,3)*khatp)/sqrt(knorm), &
!             (bbhat(:,1)*khatr+bbhat(:,2)*khatm+bbhat(:,3)*khatp)/sqrt(knorm)
!        cxitest=-(bbhat(:,1)*bhat%data(2)+bbhat(:,2)*bhat%data(3)+bbhat(:,3)*bhat%data(4))/sqrt(bdotb) + &
!             (bbhat(:,1)*khatr+bbhat(:,2)*khatm+bbhat(:,3)*khatp)*(khat*bhat)/knorm/sqrt(bdotb)
!        sxitest=(aahat(:,1)*bhat%data(2)+aahat(:,2)*bhat%data(3)+aahat(:,3)*bhat%data(4))/sqrt(bdotb) - &
!             (aahat(:,1)*khatr+aahat(:,2)*khatm+aahat(:,3)*khatp)*(khat*bhat)/knorm/sqrt(bdotb)
        !write(6,*) 'cxitest: ',cxi,cxitest
        !write(6,*) 'sxitest: ',sxi,sxitest
!        s2xi=2d0*sxi*cxi; c2xi=cxi*cxi-sxi*sxi
! roman version manifestly has s2xi^2+c2xi^2 = 1, doesn't need bperp.
!        write(6,*) 's2xi: ',s2xi
!        Be1=bhat%data(2)*ahat%data(2)+bhat%data(3)*ahat%data(3)+bhat%data(4)*ahat%data(4)
!        Be2=bhat%data(2)*bbhat(:,1)+bhat%data(3)*bbhat(:,2)+bhat%data(4)*bbhat(:,3)
!        write(6,*) 's2xi roman: ',-2*Be1*Be2/(Be1*Be1+Be2*Be2)
!        write(6,*) 'c2xi: ',c2xi
!        write(6,*) 'c2xi roman: ',(-Be1*Be1+Be2*Be2)/(Be1*Be1+Be2*Be2)
!        write(6,*) 'c2xi bp: ',(-aadotbp*aadotbp+bpdotbb*bpdotbb)/(aadotbp*aadotbp+bpdotbb*bpdotbb)
!        s2xi=sin(2d0*xi); c2xi=cos(2d0*xi)
        !angnorm=bdotk/sqrt(knorm)/sqrt(bdotb)
!        write(6,*) 'cdot: ',b(1)%data,u(1)%data,k(1)%data
!        write(6,*) 'cdot: ',bdotk,bdotb,om*om
!        write(6,*) 'ang: ',angnorm, bdotk/sqrt(om2)/sqrt(bdotb)
        ang=acos(merge(merge(angnorm,angmax,angnorm.le.angmax),angmin, &
         angnorm.ge.angmin))
        ! CHANGED SIGN following \nu = -u^\mu k_\mu = -u^t k_t = k^t in comoving ortho so g = 1 / k^t
        g=1d0/khat%data(1)
!        write(6,*) 'comoving ortho ang g: ',ang(testindx),g(testindx)
!        write(6,*) 'comoving ortho r th: ',r(testindx),th(testindx)
!        write(6,*) 'comoving ortho vectors a: ',aahat(testindx,1), &
!             aahat(testindx,2),aahat(testindx,3)
!        write(6,*) 'comoving ortho vectors k: ',khat(testindx)%data(2),&
!             khat(testindx)%data(3),khat(testindx)%data(4)
!        write(6,*) 'comoving ortho vectors b: ',bbhat(testindx,1),&
!             bbhat(testindx,2),bbhat(testindx,3)
!        write(6,*) 'comoving ortho vectors B: ',bhat(testindx)%data(2),&
!             bhat(testindx)%data(3),bhat(testindx)%data(4)
!        write(6,*) 'comoving ortho: ',u*u
!        write(6,*) 'comoving ortho: ',k*k
!        write(6,*) 'comoving ortho khat: ',khatr
!        write(6,*) 'comoving ortho khat: ',khat%data(1)
        cosne=g*sqrt(beta*beta+mus*mus*(alpha*alpha-a*a))/r
        if(any(isnan(s2xi)).or.any(isnan(c2xi))) then
!           write(6,*) 'NaN in comoving ortho aadotbp: ',aadotbp
!           write(6,*) 'NaN in comoving ortho bpdotbb: ',bpdotbb
!           write(6,*) 'NaN in comoving ortho ahat: ',aa%data(2),aa%data(3),aa%data(4)
!           write(6,*) 'NaN in comoving ortho Kpw: ',Kap1,Kap2
        endif
      end subroutine comoving_ortho_core

      subroutine comoving_ortho_debug_old(r,th,a,alpha,beta,mus,u,b,k, &
         s2xi,c2xi,ang,g,aa,aahat,bhat,khat)
! Transform four-velocity u, magnetic field b and wave-vector k to the
! co-moving orthonormal frame following Beckwith et al (2008) and
! Shcherbakov & Huang (2011)
! JAD 2/14/2011, fortran 3/27/2011
!        implicit none
        type (four_Vector), dimension(:), intent(in) :: u,b,k
        real(kind=8), dimension(:), intent(in) :: r,th
        real(kind=8), intent(in) :: a,alpha,beta,mus
        real(kind=8), intent(out), dimension(size(r)) :: s2xi,c2xi, &
             ang,g
        type (four_Vector), dimension(size(r)) :: uhat, &
         bhat,khat
        real(kind=8), dimension(size(r),10) :: metric
        real(kind=8), dimension(10,size(r)) :: tmetric
        real(kind=8), dimension(size(r)) :: gtt,gtp,grr,gmm,gpp,ut, &
             ur,um,up,urc,umc,utc,upc,np2,nm2,nr2,uhatt,uhatr,uhatm, &
             bhatt,bhatr,bhatm,bhatp,khatt,khatr,khatm,khatp,aahatt, &
             aahatm,aahatr,aahatp,delta,al1,al2,al3,uhatp,ahatt,ahatm, &
             ahatr,ahatp,bdotb,bdotk,kdotk,om,om2,knorm,bpdotbp, &
             bpdotbb,aadotbp,sxi,cxi,eps,one,mone,angnorm,angmin,angmax 
        real(kind=8), dimension(size(r),3) :: aahat &
           ,bbhat
        real(kind=8), dimension(size(r),3) :: bphat
        type (four_Vector), dimension(size(r)) :: &
             ekt,ekr,ekm,ekp,stb,srb,spb,smb,aa,ahat
        integer :: i, nr
        real(kind=8) :: kap1,kap2
        eps=epsilon(eps)
        metric=kerr_metric(r,th,a); tmetric=transpose(metric)
! Co-variant metric components
!        call assign_metric(u,tmetric); call assign_metric(b,tmetric)
!        call assign_metric(k,tmetric)
        gtt=metric(:,1) ; gtp=metric(:,4) ; grr=metric(:,5)
        gmm=metric(:,8) ; gpp=metric(:,10)
  !      write(6,*) 'ut ',size(r),size(th),a,alpha,beta,mus,size(u), &
!         size(b), size(k)
! Gather all possible four-velocity components to simplify things:
        ut=u%data(1) ; ur=u%data(2) ; um=u%data(3) ; up=u%data(4)
        utc=gtt*ut+gtp*up; upc=gpp*up+gtp*ut ; urc=grr*ur ; umc=um*gmm
        nr=size(r)
!        write(6,*) 'a: ',a
!  !      write(6,*) 'kap', nr
! Find parallel transported basis e^\th, e^ph everywhere on
! ray from complex Walker-Penrose constant
! (see 2/22 notes, Chandrasekhar 1983 and CS77, CPS80)
        Kap1=beta
        Kap2=-alpha-a*sqrt(1.-mus*mus)
  !      write(6,*) 'perpk call',Kap1,Kap2,alpha,beta
        if((Kap1.ne.0d0).or.(Kap2.ne.0d0)) then
          call transport_perpk(k,r,th,a,metric,Kap1,Kap2,al1,al2,al3)
          write(6,*) 'out', metric(1,:)
          write(6,*) 'out', al1(1)*sqrt(metric(1,5)), sqrt(metric(1,8))*al2(1), sqrt(metric(1,10))*al3(1)
        else
          al1=0d0; al2=0d0; al3=1d0/sqrt(metric(:,10))
        endif
        aa%data(1)=0d0;aa%data(2)=al1;aa%data(3)=al2;aa%data(4)=al3
        call assign_metric(aa,tmetric)
!  !      write(6,*) 'kdota: ',maxval(abs(k*aa)),maxval(abs(aa*aa-1))
!  !      write(6,*) 'perpk'
! Kulkarni+2011 normalizations:
        delta=r*r+a*a-2d0*r
!        write(6,*) 'delta: ',r,a,delta,grr
        nr2=-grr*(utc*ut+upc*up)*(1.+umc*um)
!  !      write(6,*) 'nr2: ',utc,ut,upc,up,umc,um,nr2
        nm2=gmm*(1.+umc*um)
        np2=-(utc*ut+upc*up)*delta*sin(th)*sin(th)
        ekt=u
        ekr%data(1)=(urc*ut)/sqrt(nr2); ekr%data(2)=-(utc*ut+upc*up) &
            /sqrt(nr2)
        ekr%data(3)=0d0; ekr%data(4)=(urc*up)/sqrt(nr2)
        ekm%data(1)=(umc*ut)/sqrt(nm2); ekm%data(2)=(umc*ur)/sqrt(nm2)
        ekm%data(3)=(1+umc*um)/sqrt(nm2); ekm%data(4)=(umc*up)/sqrt(nm2)
        ekp%data(1)=upc/sqrt(np2); ekp%data(2)=0d0; ekp%data(3)=0d0
        ekp%data(4)=-utc/sqrt(np2)
  !      write(6,*) 'ek'!, urc, up, umc
        call assign_metric(ekt,tmetric); call assign_metric(ekm,tmetric)
        call assign_metric(ekp,tmetric); call assign_metric(ekr,tmetric)
!        write(6,*) 'ek', maxval(abs(ekt*ekt+1)), maxval(abs(ekp*ekp-1))
!        write(6,*) 'ek', maxval(abs(ekr*ekr-1)), maxval(abs(ekm*ekm-1))
!        write(6,*) 'ek', maxval(abs(ekt*ekp)), maxval(abs(ekt*ekr)),  &
!        maxval(abs(ekt*ekm)), maxval(abs(ekr*ekm)),maxval(abs(ekr*ekp))
        uhatt=ekt*u
        uhatr=ekr*u
        uhatm=ekm*u
        uhatp=ekp*u
        bhatt=ekt*b
        bhatr=ekr*b
        bhatm=ekm*b
        bhatp=ekp*b
!        write(6,*) 'bt: ',maxval(abs(bhatt)),maxval(abs(uhatt+1))
!        write(6,*) 'bhat: ',maxval(abs(bhatr)),maxval(abs(bhatm)), &
!         maxval(abs(bhatp))
!        write(6,*) 'uhatr: ',uhatr,uhatp
        khatt=ekt*k
        khatr=ekr*k
        khatm=ekm*k
        khatp=ekp*k
        ahatt=ekt*aa
        ahatr=ekr*aa
        ahatm=ekm*aa
        ahatp=ekp*aa
!        write(6,*) 'khatt: ', &
!        maxval(abs(khatt*khatt-khatr*khatr-khatm*khatm-khatp*khatp))
 ! !      write(6,*) 'ahatt: ',maxval(abs(ahatt)),maxval(abs(ahatr*ahatr+ &
 !       ahatm*ahatm+ahatp*ahatp-1))
        uhat%data(1)=uhatt
        uhat%data(2)=uhatr
        uhat%data(3)=uhatm
        uhat%data(4)=uhatp
  !      write(6,*) 'uhat: ',maxval(abs(uhatt+1))
! now transform b, k:
        bhat%data(1)=bhatt
        bhat%data(2)=bhatr
        bhat%data(3)=bhatm
        bhat%data(4)=bhatp
        khat%data(1)=khatt
        khat%data(2)=khatr
        khat%data(3)=khatm
        khat%data(4)=khatp
        ahat%data(1)=ahatt
        ahat%data(2)=ahatr
        ahat%data(3)=ahatm
        ahat%data(4)=ahatp
  !      write(6,*) 'ahatt: ',maxval(abs(ahatt))
        call assign_metric(b,tmetric)
        call assign_metric(u,tmetric)
!        write(6,*) 'udotu: ',u*u
!        write(6,*) 'fb: ',maxval(abs(b*u)),maxval(sqrt(b*b)), &
!         maxval(abs(1+u*u)),maxval(abs(k*k))
!        write(6,*) 'bhat: ',b*b-(bhatt**2+bhatr**2+bhatm**2+bhatp**2)
!        call assign_metric(khat,metric)
!        write(6,*) 'uhatr: ',uhatr
!        write(6,*) 'uhatm: ',uhatm
!        write(6,*) 'uhatp: ',uhatp
        kdotk=khatr*khatr+khatm*khatm+khatp*khatp
        bdotk=b*k
        bdotb=b*b
!!        write(6,*) 'kk: ',k*k
  !      write(6,*) 'bdotb: ',bdotb!,bdotk,sqrt(bdotb),maxval(bdotk/sqrt(bdotb))
        om=-k*u; om2=om*om
!        write(6,*) 'khat: ',-khatt*uhatt,om
!        write(6,*) 'uhat: ',uhatt
!  !      write(6,*) 'om: ',om2, om, om2-kdotk
        aahat(:,1)=ahatr-khat%data(2)*ahatt/khat%data(1)
        aahat(:,2)=ahatm-khat%data(3)*ahatt/khat%data(1)
        aahat(:,3)=ahatp-khat%data(4)*ahatt/khat%data(1)
!  !      write(6,*) 'sane: ',maxval(abs(ahatt-khat%data(1)*ahatt/khat &
!        %data(1)))
!  !      write(6,*) 'adota: ',maxval(abs(aahat(:,1)**2+aahat(:,2)**2 &
!         +aahat(:,3)**2-1))
        knorm=khatr*khatr+khatm*khatm+khatp*khatp
!  !      write(6,*) 'knorm: ',knorm
        bbhat(:,1)=(aahat(:,2)*khatp-aahat(:,3)*khatm)/sqrt(knorm)
        bbhat(:,2)=(aahat(:,3)*khatr-aahat(:,1)*khatp)/sqrt(knorm)
        bbhat(:,3)=(aahat(:,1)*khatm-aahat(:,2)*khatr)/sqrt(knorm)
 ! !      write(6,*) 'aabb: ',maxval(abs(aahat(:,1))),
 !    &  maxval(abs(aahat(:,2))),maxval(abs(aahat(:,3))),
 !    &  maxval(abs(bbhat(:,1))),maxval(abs(bbhat(:,2))),
 !    &  maxval(abs(bbhat(:,3)))
!  !      write(6,*) 'aabb2: ',aahat
!  !      write(6,*) 'khat: ',khatr,khatm,khatp
!  !      write(6,*) 'aadotbb: ',bbhat(:,1)*aahat(:,1)+aahat(:,2)*bbhat(
!         :,2)+aahat(:,3)*bbhat(:,3)
!        write(6,*) 'bbdotbb: ',maxval(abs(bbhat(:,1)**2+bbhat(:,2)**2
!          +bbhat(:,3)**2-1))
! Find perpendicular magnetic field:
!  !      write(6,*) 'bdotk: ',maxval(abs(khatr*bhatr+khatm*bhatm+
!         khatp*khatp-bdotk))
!  !      write(6,*) 'abperp: ',maxval(abs(aahat(:,1)*bbhat(:,1)+
!         aahat(:,2)*bbhat(:,2)+aahat(:,3)*bbhat(:,3))),
!        maxval(abs(aahat(:,1)*khatr+aahat(:,2)*khatm+aahat(:,3)*khatp)),
!        maxval(abs(bbhat(:,1)*khatr+bbhat(:,2)*khatm+bbhat(:,3)*khatp)),
!        maxval(abs(dot_product(aahat,bbhat)))
        bdotk=bhat%data(2)*khatr+bhat%data(3)*khatm+bhat%data(4)*khatp
  !      write(6,*) 'bphat: ',maxval(abs(bhat%data(1)))
  !      write(6,*) 'bphat: ',maxval(abs(bhat%data(2)))
        bphat(:,1)=bhat%data(2)-khatr*bdotk/knorm
        bphat(:,2)=bhat%data(3)-khatm*bdotk/knorm
        bphat(:,3)=bhat%data(4)-khatp*bdotk/knorm
  !      write(6,*) 'bphat: ',maxval(abs(bphat(:,1)*khatr+bphat(:,2)*
  !   &   khatm+bphat(:,3)*khatp))
        aadotbp=dot_product(aahat,bphat)
        bpdotbp=dot_product(bphat,bphat)
        bpdotbb=dot_product(bphat,bbhat)
        write(6,*) 'aabb: ',maxval(abs(dot_product(aahat,aahat)-1)), &
        maxval(abs(dot_product(bbhat,bbhat)-1))
        write(6,*) 'dot: ',maxval(aadotbp),bpdotbp,bpdotbb
! Now compute angle between polarization basis and magnetic field:
        one=1d0; mone=-one; angmin=-.9999d0; angmax=-1d0*angmin
        where(bpdotbp.gt.0d0)
           aadotbp=bhat%data(2)*aahat(:,1)+bhat%data(3)*aahat(:,2)+bhat%data(4)*aahat(:,3)
           bpdotbb=bhat%data(2)*bbhat(:,1)+bhat%data(3)*bbhat(:,2)+bhat%data(4)*bbhat(:,3)
           s2xi=-2d0*aadotbp*bpdotbb/(aadotbp**2d0+bpdotbb**2d0)
           c2xi=(bpdotbb*bpdotbb-aadotbp*aadotbp)/(aadotbp**2d0+bpdotbb**2d0)
           sxi=aadotbp/sqrt(bpdotbp)
           cxi=-bpdotbb/sqrt(bpdotbp)
           angnorm=bdotk/sqrt(knorm)/sqrt(bdotb)
!        where(
        elsewhere
           sxi=0d0
           cxi=one
           angnorm=0.5d0
        endwhere
!        write(6,*) 'cxi: ',cxi
!        s2xi=2.*sxi*cxi; c2xi=cxi*cxi-sxi*sxi
!        write(6,*) 'angnorm', angnorm
!        write(6,*) 'bdotk: ',bdotk/sqrt(knorm)/sqrt(bdotb)
!        write(6,*) 'knorm: ',knorm
!        write(6,*) 'bdotb: ',bdotb
!        write(6,*) 'bdotk: ',(b*k)**2/om2/bdotb
  !      write(6,*) 'mone',merge(merge(angnorm,one,angnorm.le.1d0),mone,
 !    &   angnorm.ge.-1d0)
!  !      write(6,*) 'mone',merge(angnorm,mone,angnorm.ge.-1d0)
        ang=acos(merge(merge(angnorm,angmax,angnorm.le..9999d0),angmin, &
         angnorm.ge.-.9999d0))
!        write(6,*) 'sc: ',cxi*cxi+sxi*sxi
!        write(6,*) 'ang: ',ang
        g=-1./khat%data(1)
      end subroutine comoving_ortho_debug_old

      subroutine calc_polar_psi(r,muf,q2,a,alpha,beta,rshift,mus,p, &
           c2psi,s2psi,cosne)
      real(kind=8), dimension(:), intent(in) :: r,muf,q2, &
           alpha,beta,rshift
      type (four_vector), dimension(:), intent(in) :: p
      real(kind=8), intent(in) :: a, mus
      type (four_vector), dimension(size(r)) :: f
      real(kind=8), dimension(size(r)), intent(out) :: cosne,s2psi,c2psi
      real(kind=8), dimension(size(r)) :: kappa1, kappa2, gammac, num, denom, polarpsi 
      complex, dimension(size(r)) :: kappa_pw
!      write(6,*) 'polar psi: ',size(r)
! Components of vector in plane of disk
      f=calc_polvec(r,muf,p,a,0.d0)
! penrose-walker constant for that vector
      kappa_pw=calc_kappapw(a,r,muf,p,f)
      kappa2=real(kappa_pw)
      kappa1=-aimag(kappa_pw)
      gammac=-alpha-a*(1.d0-mus**2)
! CPS80 but I think these are relative to e^\theta, e.g., \beta axis
!      c2psi=(gammac*kappa1-beta*kappa2)/(gammac**2d0+beta**2d0)
!      s2psi=(-gammac*kappa2-beta*kappa1)/(gammac**2d0+beta**2d0)
! Eric's (relative to \alpha axis)
      denom=(beta*kappa2-gammac*kappa1)!/(gammac**2d0+beta**2d0)
      num=(-beta*kappa1-gammac*kappa2)!/(gammac**2d0+beta**2d0)
! Eric's equation for polarpsi (seems to disagree with CPS80...)
      polarpsi=atan2(denom,num)
      s2psi=sin(2d0*polarpsi); c2psi=cos(2d0*polarpsi)
! polarization angle for that vector transported to infinity
!      num=-gammac*kappa2-beta*kappa1
!      denom=gammac*kappa1-beta*kappa2
!      polarpsi=0.5d0*atan(num/denom)
!      s2psi=2d0*num*denom/(gammac**2d0+beta**2d0)**2d0
!      s2psi=2d0*kappa1*kappa2*(beta**2d0-gammac**2d0)/(beta**2d0+gammac**2d0)**2d0
!      c2psi=(denom**2d0-num**2d0)/(gammac**2d0+beta**2d0)**2d0
!      c2psi=((gammac**2d0-beta**2d0)*kappa1**2d0-4d0*beta*gammac*kappa1*kappa2+ &
!           (beta**2d0-gammac**2d0)*kappa2**2d0)/(gammac**2d0+beta**2d0)**2d0
!      c2psi=sqrt(1d0-s2psi**2d0)
!      c2psi=2d0*denom**2d0/(gammac**2d0+beta**2d0)**2d0-1d0
!      c2psi=cos(2d0*polarpsi); s2psi=sin(2d0*polarpsi)
!      write(6,*) 's2psi: ',2d0*num*denom/(gammac**2d0+beta**2d0)**2d0, sin(2d0*polarpsi)
!      write(6,*) 'c2psi: ',(denom**2d0-num**2d0)/(gammac**2d0+beta**2d0)**2d0, cos(2d0*polarpsi)
      cosne=rshift*sqrt(q2)/r
      end subroutine

      function calc_polvec(r,mu,p,a,psi) result(fourf)
      real(kind=8), dimension(:), intent(in) :: r, mu
      real(kind=8), intent(in) :: a, psi
      type (four_vector), dimension(:), intent(in) :: p
      real(kind=8), dimension(size(r)) :: ptt, prt, ptht, ppht, &
           vel
      real(kind=8), dimension(10,size(r)) :: tmetric
      real(kind=8), dimension(size(r),10) :: metric
      type (four_vector), dimension(size(r)) :: fourf
      real(kind=8), dimension(size(r)) :: delta,ar,om,rho,frl,fthl, &
           fphl,frp,fthp,fphp,fr,fth,fph,normf,epsi,enu
! Calculate polarization vector assuming convention f^0=0, using Agol (1997)
! LNRF momentum:
      delta=r**2-2.d0*r+a**2
      AR=(R*R+A*A)**2-A*A*delta*(1.d0-mu*mu)
      OM=2.d0*A*R/AR; rho=r**2+a**2*mu**2
      ptt=r*sqrt(delta/AR)*p%data(1)
      prt=r/sqrt(delta)*p%data(2)
      ptht=r*p%data(3)
      ppht=sqrt(AR)/r*(p%data(4)-om*p%data(1))
! LNRF velocity for a Keplerian disk:
      vel=1.d0/(r**(3.d0/2.d0)+a) 
      EPSI=SQRT(1.d0-mu*mu)*SQRT(AR/RHO) 
      ENU=SQRT(delta*RHO/AR)
      vel=EPSI/ENU*(vel-om)
      frl=sqrt(delta)/r*(vel*(ptt-prt**2/ptt)-ppht)
      fthl=-vel*prt*ptht/ptt/r
      fphl=r*prt/sqrt(AR)*(1.d0-vel*ppht/ptt)
      frp=sqrt(delta)*ptht*prt/r*(-1.d0+vel*ppht/ptt)
      fthp=1.d0/r*(prt**2+(1.d0+vel**2)*ppht**2-2.d0*vel*ppht*ptt+vel*ptht**2*ppht/ptt)
      fphp=r*ptht/sqrt(AR)*(-(1.d0+vel**2)*ppht+vel*ptt+vel*ppht**2/ptt)
      fr=cos(psi)*frl+sin(psi)*frp
      fth=cos(psi)*fthl+sin(psi)*fthp
      fph=cos(psi)*fphl+sin(psi)*fphp
      fourf%data(1)=0.*fr
      fourf%data(2)=fr
      fourf%data(3)=fth
      fourf%data(4)=fph
      ! Now normalize:
      metric=kerr_metric(r,acos(mu),a); tmetric=transpose(metric)
      call assign_metric(fourf,tmetric)
      normf=fourf*fourf
      fourf%data(1)=fourf%data(1)/sqrt(normf)
      fourf%data(2)=fourf%data(2)/sqrt(normf)
      fourf%data(3)=fourf%data(3)/sqrt(normf)
      fourf%data(4)=fourf%data(4)/sqrt(normf)
      return
      end function calc_polvec

      function calc_kappapw(a,r,mu,p,f) result(kappapw)
      real(kind=8), intent(in), dimension(:) :: r,mu
      type (four_vector), intent(in), dimension(:) :: p
      type (four_vector), intent(in), dimension(:) :: f
      complex, dimension(size(r)) :: kappapw
      real(kind=8), intent(in) :: a
      real(kind=8), dimension(size(r)) :: alpha,beta
      ! Calculates the complex Penrose-Walker constant for some four-vector
      ! f perpendicular to p
      ! JAD 10/23/2008
      ! Calculate complex Penrose-Walker constant
      alpha=(p%data(1)*f%data(2)-p%data(2)*f%data(1))+a*(1.d0-mu**2)* &
           (p%data(2)*f%data(4)-p%data(4)*f%data(2))
      beta=(r**2+a**2)*sqrt(1.d0-mu**2)*(p%data(4)*f%data(3)-p%data(3)*f%data(4))& 
           -a*sqrt(1.d0-mu**2)*(p%data(1)*f%data(3)-p%data(3)*f%data(1))
      kappapw=cmplx(alpha,-beta)*cmplx(r,-a*mu)
      return
      end function calc_kappapw

      function surf_integral(x,r,th,ph,a) result (s)
      ! calculate radial profile of angle-integrated quantity x assuming uniform ph
      ! JAD 2/23/2013
      real(8), intent(in), dimension(:,:,:) :: x,r,th,ph
      real(8), intent(in) :: a
      real(8), dimension(size(r,1),size(r,2),size(r,3)) :: dth, gdet
      integer :: nx1,nx2,nx3
      real(8) :: dph
      real(8), dimension(size(r,1)) :: s
      dph=ph(1,1,2)-ph(1,1,1)
      gdet=(r**2d0+a**2d0*cos(th)**2d0)*sin(th)
      nx1=size(r,1); nx2=size(r,2); nx3=size(r,3)
      dth(:,2:nx2,:)=th(:,2:nx2,:)-th(:,1:nx2-1,:)
      dth(:,1,:)=dth(:,nx2,:)
      s=sum(sum(x*gdet*dth*dph,3),2)
      end function surf_integral

      function th_integral(x,r,th,ph,a,onlyth) result (s)
      ! calculate radial profile of angle-integrated quantity x assuming uniform ph
      ! JAD 2/23/2013
      real(8), intent(in), dimension(:,:,:) :: x,r,th,ph
      real(8), intent(in) :: a
      real(8), dimension(size(r,1),size(r,2),size(r,3)) :: dth, gdet
      integer :: nx1,nx2,nx3
      real(8) :: dph
      real(8), dimension(size(r,1),size(r,3)) :: s
      integer, intent(in) :: onlyth
      dph=ph(1,1,2)-ph(1,1,1)
      gdet=(r**2d0+a**2d0*cos(th)**2d0)*sin(th)
      nx1=size(r,1); nx2=size(r,2); nx3=size(r,3)
      dth(:,2:nx2,:)=th(:,2:nx2,:)-th(:,1:nx2-1,:)
      dth(:,1,:)=dth(:,nx2,:)
      s=sum(x*gdet*dth*dph,2)*nx3
      end function th_integral

      function ph_integral(x,r,th,ph,a,onlyth,onlyph) result (s)
      ! calculate radial profile of angle-integrated quantity x assuming uniform ph
      ! JAD 2/23/2013
      real(8), intent(in), dimension(:,:,:) :: x,r,th,ph
      real(8), intent(in) :: a
      real(8), dimension(size(r,1),size(r,2),size(r,3)) :: dth, gdet
      integer :: nx1,nx2,nx3
      real(8) :: dph
      real(8), dimension(size(r,1),size(r,2)) :: s
      integer, intent(in) :: onlyth, onlyph
      dph=ph(1,1,2)-ph(1,1,1)
      gdet=(r**2d0+a**2d0*cos(th)**2d0)*sin(th)
      nx1=size(r,1); nx2=size(r,2); nx3=size(r,3)
      dth(:,2:nx2,:)=th(:,2:nx2,:)-th(:,1:nx2-1,:)
      dth(:,1,:)=dth(:,nx2,:)
      s=sum(x*gdet*dth*dph,3)*nx2
!      write(6,*) 'ph integral: ',minval(s),maxval(s),minval(x),maxval(x)
      end function ph_integral

      function calc_u0(metric,vr,vth,vph) result(u0)
! calculate u0 from three-velocity in BL coordinates
        real(8), intent(in), dimension(:) :: vr,vth,vph
        real(8), intent(in), dimension(:,:) :: metric
        real(8), dimension(size(vr)) :: u0
        u0 = sqrt(-1d0/(metric(:,1) + metric(:,5)*vr**2d0 + metric(:,8)*vth**2d0 &
             + metric(:,10)*vph**2d0 + 2d0*metric(:,4)*vph))
      end function calc_u0

      subroutine calc_rms_constants(a,ems,lms,rms)
        real(8), intent(in) :: a
        real(8), intent(out) :: ems,lms,rms
        real(8) :: v
        rms = dble(calc_rms(real(a)))
        v = 1d0/sqrt(rms)
        ems = (1d0-2d0*v*v+a*v*v*v)/sqrt(1d0-3d0*v*v+2d0*a*v*v*v)
        lms = rms*v*(1d0-2d0*a*v*v*v+a*a*v*v*v*v)/sqrt(1d0-3d0*v*v+2d0*a*v*v*v)
        return
      end subroutine calc_rms_constants

      function calc_plunging_vel(a,r) result(fu)
        real(8), intent(in) :: a
        real(8), dimension(:), intent(in) :: r
        type (four_vector), dimension(size(r)) :: fu
        real(8) :: ems,lms,rms,rh
        real(8), dimension(size(r)) :: pt,pr,pphi,th,denom
        real(8), dimension(size(r),10) :: metriccon
! this is way based on Hughes (2000, 2001), Schnittman PhD. equivalently there is analytic way from old code that seems to give ~same result.
        call calc_rms_constants(a,ems,lms,rms)
!        write(6,*) 'kerr plunging vel rms constants: ',ems,lms,rms,minval(r),maxval(r)
        th = pi/2d0
        metriccon = kerr_metric(r,th,a,1)
        pt = -metriccon(:,1)*ems + metriccon(:,4)*lms
        denom = -metriccon(:,5)*(1d0+metriccon(:,1)*ems*ems-2d0*metriccon(:,4)*ems*lms+metriccon(:,10)*lms*lms)
        where(denom.gt.0.)
           pr = -sqrt(denom)
        elsewhere
           pr = 0d0
        endwhere
!        pr = -sqrt(-metriccon(:,5)*(1d0+metriccon(:,1)*ems*ems-2d0*metriccon(:,4)*ems*lms+metriccon(:,10)*lms*lms))
        pphi = -metriccon(:,4)*ems+metriccon(:,10)*lms
        fu%data(1)=pt
        fu%data(2)=pr
        fu%data(3)=0d0
        fu%data(4)=pphi
!        write(6,*) 'kerr plunging vel fu: ',minval(pt),maxval(pt),minval(pr),maxval(pr),minval(pphi),maxval(pphi)
      end function calc_plunging_vel

      function rms_vel(a,th,r) result(fu)
        real(8), intent(in) :: a
        real(8), intent(in), dimension(:) :: th,r
        type (four_vector), dimension(size(r)) :: fu,fueq
        real(8), dimension(size(r)) :: vrl,vtl,vpl,vrr,vtt,vpp,u0,theq
        real(8), dimension(size(r),10) :: metric
! plunging four-velocity in equatorial plane
        fueq = calc_plunging_vel(a,r)
        theq = pi/2d0
! transform to lnrf
        call lnrf_frame(fueq%data(2)/fueq%data(1),fueq%data(3)/fueq%data(1), &
             fueq%data(4)/fueq%data(1),r,a,theq,vrl,vtl,vpl)
!        write(6,*) 'kerr rms_vel lnrf 1: ',maxval(vrl),minval(vrl),maxval(vpl),minval(vpl)
! invert but at real \theta outside equatorial plane
        call lnrf_frame(vrl,vtl,vpl,r,a,th,vrr,vtt,vpp,1)
        metric=kerr_metric(r,th,a)
        u0=calc_u0(metric,vrr,vtt,vpp)
!        write(6,*) 'kerr rms_vel u0: ',minval(u0),maxval(u0)
        fu%data(1)=u0; fu%data(2)=u0*vrr; fu%data(3)=u0*vtt; fu%data(4)=u0*vpp
! now test
        call assign_metric(fu,transpose(metric))
!        write(6,*) 'kerr rms_vel udotu: ',maxval(abs(fu*fu+1d0))
      end function rms_vel

    end module kerr
