   module fluid_model_harm3d

      use class_four_vector
      use interpolate, only: interp
      use kerr, only: kerr_metric, lnrf_frame, uks2ubl, bl2ks
      use phys_constants, only: pi
      use math, only: zbrent
      implicit none

      namelist /harm/  dfile, hfile, nt, indf, gfile

      character(len=100) :: dfile, hfile, gfile
      integer :: nx1, nx2, nx3, n, ndumps, nt, indf, dlen, nhead
      integer :: ndim=3
      integer, dimension(:), allocatable :: dumps
      real :: tstep=20., h, asim, dx1, dx2, dx3, gam, startx1, startx2, startx3, toffset=0.
      real, dimension(:), allocatable :: t
      real, dimension(:), allocatable :: x1_arr, x2_arr, x3_arr, r_arr, th_arr, ph_arr
      real, dimension(:), allocatable :: rho_arr, p_arr, u0_arr, vrl_arr, &
        vpl_arr, vtl_arr, temp_arr
      real, dimension(:), allocatable :: b0_arr, br_arr, bth_arr, bph_arr

      interface init_harm3d_data
        module procedure init_harm3d_data
      end interface

      interface del_harm3d_data
        module procedure del_harm3d_data
      end interface

      interface initialize_harm3d_model
        module procedure initialize_harm3d_model
      end interface

      interface transformbl2mksh
         module procedure transformbl2mksh
      end interface

      interface umksh2uks
         module procedure umksh2uks
      end interface

      interface harm3d_vals
        module procedure harm3d_vals
      end interface

      contains

        function findx2harm(x2,args) result(diff)
        real(kind=8), intent(in) :: x2
        real(kind=8), intent(in), dimension(:) :: args
        real(kind=8) :: h, theta, diff
        h=args(2); theta=args(1)
        diff=theta-pi*x2
        end function findx2harm

        subroutine transformbl2mksh(r,th,phi,x1,x2,x3,h)
        ! to 3D on 19/11/15 JD
        ! transform Boyer-Lindquist coordinates to modified Kerr-Schild coordinates used by HARM
        real, intent(in), dimension(:) :: r,th,phi
        real, intent(out), dimension(size(r)) :: x1,x2,x3
        real, intent(in) :: h
        real :: hval
        real, dimension(2) :: args
        integer :: i
        x1=log(r)
        x3=phi
        x2=th/pi
!        hval=1d0
!        args=(/th(1),hval/)
!        do i=1,size(r)
!           args(1)=th(i)
!           x2(i)=zbrent(findx2harm,0d0,1d0,dble(args),1d-6)
!        enddo
        end subroutine transformbl2mksh

        function umksh2uks(fum,r,x2,hin) result(uks)
          ! Converts a Kerr-Schild spherical 4-velocity to a Boyer-Lindquist one.
          ! JAD 11/10/2008
          ! Do simple transformation first:
          ! to 3D on 19/11/15 JD
          real, dimension(:), intent(in) :: r, x2
          type (four_vector), dimension(:), intent(in) :: fum
          type (four_vector), dimension(size(r)) :: uks
          real, intent(in), optional :: hin
          real :: hval
          real(kind=8), dimension(size(r)) :: ur, uth, uph
          real, dimension(size(r)) :: dthdx2
          write(6,*) 'read harm umksh2uks present h'
          if (present(hin)) then
             hval=hin
          else
             hval=0.3
          endif
          write(6,*) 'hval: ',present(hin), hval
          ur=dble(r)*fum%data(2)
          ! Compute partial derivatives:
          dthdx2=pi
          uth=fum%data(3)*dble(dthdx2)
          uph=fum%data(4)
          uks=fum
          uks%data(2)=ur
          uks%data(3)=uth
          uks%data(4)=uph
        end function umksh2uks

        subroutine harm3d_vals(x0,a,rho,p,b,u,bmag)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real(kind=8), dimension(size(x0)) :: done,pfac,nfac
        real, dimension(size(x0)) :: x3,x2,x1,zm,theta,fone, &
         vpl0,vrl0,vtl0,rd,td, pd,rttd,zr,dzero, &
         vr0,vth0,vph0,bth,dummy,tt,ttd,zt,zphi,zpp
        real, dimension(nx1) :: uniqx1,uniqr
        real, dimension(nx2) :: uniqx2,uniqth
        real, dimension(nx3) :: uniqx3,uniqph
        real, dimension(size(x0),2**(ndim+1)) :: ppi,rhoi,vrli,vtli, &
         vpli,bri,bthi,bphi,b0i,u0i,ri,thi,phii
!        integer, dimension(size(x0)*(2**ndim)) :: sindx
        real, dimension(size(x0)) :: dth, minph
        integer, dimension(size(x0)*2*(2**ndim)) :: indx
        integer, dimension(size(x0)) :: lx1,lx2, lx3, ux1,ux2,ux3,x1l,x1u,x2l,x2u,x3l,x3u, &
         umax,one,tindx
        integer :: npts,i,maxtindx
        real, dimension(size(x0)), intent(out) :: rho,p,bmag
        type (four_Vector), intent(out), dimension(size(x0)) :: u,b
        ! Interpolates HARM data to input coordinates
        ! JAD 3/20/2009, fortran 11/12/2012
!        write(6,*) 'harm: ',size(x0)
        !nx1=floor(sqrt(real(size(x1_arr)))); nx2=nx1
!        write(6,*) 'harm3d vals: ',nx1,nx2,nx3,size(x1_arr),size(x2_arr),size(x3_arr)
!        write(6,*) 'harm3d vals'
        umax=nx2-1; one=1
        dzero=0d0; done=1d0; fone=1.
        uniqx3=x3_arr(1:nx3)
        uniqx2=x2_arr(1:nx3*(nx2-1)+1:nx3)
!       write(6,*) 'between'
        uniqx1=x1_arr(1:nx3*(nx2-1)*nx1+1:nx3*nx2)
        uniqr=exp(uniqx1)
        uniqth=pi*uniqx2
        uniqph=uniqx3
!        write(6,*) 'uniqr: ',minval(uniqr), maxval(uniqr)
!        write(6,*) 'uniqth: ',minval(uniqth), maxval(uniqth)
!        write(6,*) 'uniqph: ',minval(uniqph), maxval(uniqph)
        npts=size(x0)
        theta=x0%data(3)
        zr=x0%data(2)
        zt=x0%data(1)
        tt=bl2ks(dble(zr),dble(zt),dble(a),1)-bl2ks(dble(zr(1)),0d0,dble(a),1)
        tt=-tt
!        write(6,*) 'tt: ',minval(x0%data(1)),maxval(x0%data(1)),minval(tt),maxval(tt)
!        ztt=-ztt
!        tt=ztt-min(ztt)+tabs0
        zpp=x0%data(4)
        zphi=bl2ks(dble(zr),dble(zpp),dble(a))
        zphi=mod(zphi,(2.*pi))
        where(zphi.lt.0.)
            zphi=zphi+2.*pi
        endwhere
!        write(6,*) 'harm3d vals transform'
        call transformbl2mksh(zr,theta,zphi,x1,x2,x3,h)

!        write(6,*) 'transform r: ',maxval(zr), minval(zr), maxval(x1), minval(x1)
!        write(6,*) 'transform th: ',maxval(theta), minval(theta), minval(x2), maxval(x2)
!        write(6,*) 'transform phi: ',maxval(zphi), minval(zphi), minval(x3), maxval(x3)
!        write(6,*) 'dif x3: ',(uniqx3(nx3)-uniqx3(1)),x3,uniqx3
        lx1=floor((x1-uniqx1(1))/(uniqx1(nx1)-uniqx1(1))*(nx1-1))+1
        ux1=lx1+1
        lx2=floor((x2-uniqx2(1))/(uniqx2(nx2)-uniqx2(1))*(nx2-1))+1
        ux2=lx2+1
        lx3=floor((x3-uniqx3(1))/(uniqx3(nx3)-uniqx3(1))*(nx3-1))+1
        ux3=lx3+1
!        lx2=merge(merge(lx2,one,lx2.ge.1),umax,lx2.le.(nx2-1))
        lx2=merge(lx2,one,lx2.ge.1)
        lx2=merge(lx2,nx2,lx2.le.nx2)
        ux2=merge(ux2,one,ux2.ge.1)
        ux2=merge(ux2,nx2,ux2.le.nx2)
!        write(6,*) 'nx1 nx2 nx3: ',nx1,nx2,nx3,minval(zphi),maxval(zphi),minval(zr),maxval(zr),minval(theta), &
!             maxval(theta)
!       write(6,*) 'ux lx: ',minval(lx1),maxval(ux1),minval(lx2),maxval(ux2), minval(lx3), maxval(ux3)
! Deal with poles
!        write(6,*) 'poles'
        where(ux2.ne.lx2)
           dth=uniqth(ux2)-uniqth(lx2)
        elsewhere
            dth=uniqth(1)*2
        endwhere
! periodic in phi
        minph=uniqph(lx3)
        where(ux3.gt.nx3)
            ux3=1
        endwhere
        where(lx3.lt.1)
            lx3=nx3
            minph=uniqph(lx3)-2.*pi
        endwhere
!        write(6,*) 'uniform phi'
! uniform in phi
        pd=(zphi-minph)/(uniqph(2)-uniqph(1))
        td=abs(theta-uniqth(lx2))/dth
        rd=(zr-uniqr(lx1))/(uniqr(ux1)-uniqr(lx1))
! When the geodesic is inside the innermost zone center located outside the horizon, 
! use nearest neighbor rather than interpolate
        where(uniqr(lx1).le.(1.+sqrt(1.-a**2.)))
           rd=1.
           pfac=1e-3
           nfac=1e-3
        elsewhere
           pfac=1.
           nfac=1.
        endwhere
!        write(6,*) 'coords: ',
!        write(6,*) 'rd td pd: ',minval(rd),maxval(rd),minval(td),maxval(td),minval(pd),maxval(pd)
!        write(6,*) 'ux lx: ',minval(lx1),maxval(ux1),minval(lx2),maxval(ux2),minval(lx3),maxval(ux3)
        ! th is fastest changing index
        x3l=lx3-1; x3u=ux3-1
        x2l=(lx2-1)*nx3 ; x2u=(ux2-1)*nx3
        x1l=(lx1-1)*nx2*nx3; x1u=(ux1-1)*nx2*nx3
! Indices for nearest neighbors in 2D
!        write(6,*) 'sindx', size(sindx)
!        sindx=(/(/x1l+x2l+1/),(/x1l+x2u+1/),(/x1u+x2l+1/),(/x1u+x2u+1/)/)
!        write(6,*) 'sindx', maxval(sindx), sindx(npts+1), x1l(1)+x2u(1)+1, size(sindx)
! Now find timestep indices:
        maxtindx=nt-2
        where(floor(tt/tstep).le.maxtindx)
           tindx(:)=floor(tt/tstep)
           ttd=(tt-tindx*tstep)/tstep
        elsewhere
           tindx=maxtindx+1
           ttd=0.
        endwhere
!        write(6,*) 'tindx: ',minval(tindx),maxval(tindx)
! Create index array across times and 2D nearest neighbors (first half at first time, second half at second time):
! Need to make sure there aren't integer problems when these indices get very large
!        write(6,*) 'indx', size(indx), 2**ndim, maxval(tindx), minval(tindx), tstep
        indx=(/(/x1l+x2l+x3l+1+n*tindx/),(/x1l+x2u+x3l+1+n*tindx/),(/x1u+x2l+x3l+1+n*tindx/), &
            (/x1u+x2u+x3l+1+n*tindx/),(/x1l+x2l+x3u+1+n*tindx/), &
            (/x1l+x2u+x3u+1+n*tindx/),(/x1u+x2l+x3u+1+n*tindx/),(/x1u+x2u+x3u+1+n*tindx/), &
            (/x1l+x2l+x3l+1+n*(tindx+1)/),(/x1l+x2u+x3l+1+n*(tindx+1)/),(/x1u+x2l+x3l+1+n*(tindx+1)/), &
            (/x1u+x2u+x3l+1+n*(tindx+1)/),(/x1l+x2l+x3u+1+n*(tindx+1)/),(/x1l+x2u+x3u+1+n*(tindx+1)/), &
            (/x1u+x2l+x3u+1+n*(tindx+1)/),(/x1u+x2u+x3u+1+n*(tindx+1)/)/)
! keep indx from going out of bounds when loading one time slice
        where(indx.gt.size(rho_arr))
           indx=indx-n
        endwhere
!        do i=1,2**ndim
!           indx((i-1)*npts+1:i*npts)=sindx((i-1)*npts+1)+npts*tindx
!           indx((2**ndim+i-1)*npts+1:(i+2**ndim)*npts)=sindx((i-1)*npts+1)+npts*(tindx+1)
!        end do
!        write(6,*) 'after indx', maxval(indx), maxval(sindx), indx(npts+1), n*maxval(tindx+1),size(rho_arr),maxloc(indx)!size(rho_arr(indx))
        rhoi=reshape(rho_arr(indx),(/npts,2**(ndim+1)/))
        vrli=reshape(vrl_arr(indx),(/npts,2**(ndim+1)/))
        vtli=reshape(vtl_arr(indx),(/npts,2**(ndim+1)/))
        vpli=reshape(vpl_arr(indx),(/npts,2**(ndim+1)/))
        ppi=reshape(p_arr(indx),(/npts,2**(ndim+1)/))
        b0i=reshape(b0_arr(indx),(/npts,2**(ndim+1)/))
        bri=reshape(br_arr(indx),(/npts,2**(ndim+1)/))
        bthi=reshape(bth_arr(indx),(/npts,2**(ndim+1)/))
        bphi=reshape(bph_arr(indx),(/npts,2**(ndim+1)/))
        u0i=reshape(u0_arr(indx),(/npts,2**(ndim+1)/))
! coordinates for debugging
        ri=reshape(r_arr(indx),(/npts,2**(ndim+1)/))
        thi=reshape(th_arr(indx),(/npts,2**(ndim+1)/))
        phii=reshape(ph_arr(indx),(/npts,2**(ndim+1)/))
        rttd=0.
        rho=merge(interp(rhoi,rttd,pd,rd,td),dzero,x1.gt.uniqx1(1))*nfac
        p=merge(interp(ppi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))*pfac
!        write(6,*) 'rho: ', rho, p
        vrl0=merge(interp(vrli,rttd,pd,rd,td),dzero,x1.gt.uniqx1(1))
        vtl0=merge(interp(vtli,rttd,pd,rd,td),dzero,x1.gt.uniqx1(1))
!        write(6,*) 'v'
        b%data(1)=merge(interp(b0i,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
        b%data(2)=merge(interp(bri,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'b'
        b%data(3)=merge(interp(bthi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
        b%data(4)=merge(interp(bphi,rttd,pd,rd,td),fone,x1.gt.uniqx1(1))
 !       write(6,*) 'u'
        u%data(1)=merge(dble(interp(u0i,rttd,pd,rd,td)),done,x1.gt.uniqx1(1))
        vpl0=merge(interp(vpli,rttd,pd,rd,td),dzero,x1.gt.uniqx1(1))
!        write(6,*) 'after assign'
        ! Protect azimuthal velocities at poles
        ! Compute magnitude of interpolated b-field and force b^2 > 0 (need to look into numerical issues here):
        call assign_metric(b,transpose(kerr_metric(zr, &
         real(x0%data(3)),a)))
        bmag=b*b; bmag=merge(bmag,dzero,bmag.ge.0d0)
        bmag=sqrt(bmag)
        ! Get four-velocity:
        call lnrf_frame(vrl0,vtl0,vpl0,zr,a,real(x0%data(3)) &
         ,vr0,vth0,vph0,1)
        u%data(2)=u%data(1)*dble(vr0)
        u%data(3)=u%data(1)*dble(vth0)
        u%data(4)=u%data(1)*dble(vph0)
        call assign_metric(u,transpose(kerr_metric(zr,real(x0%data(3)) &
        ,a)))
!        write(6,*) 'min vals u', minval(u%data(1)), maxval(abs(u*u+1))
!        write(6,*) 'harm vals r: ',zr
!        write(6,*) 'harm vals rho: ',rho
!        write(6,*) 'harm vals rd: ',rd
!        write(6,*) 'harm vals ux: ',nx1,ux1
!        write(6,*) 'harm vals udotu: ',(abs(u*u+1))(minloc(zr))
!        write(6,*) 'harm vals minloc: ',minloc(zr),rd(minloc(zr)), &
!             td(minloc(zr)),x0(minloc(zr))%data(3),uniqth(lx2(minloc(zr))), &
!             lx2(minloc(zr)),pd(minloc(zr)),lx1(minloc(zr)),ux1(minloc(zr)), &
!             uniqx1(1),x1(minloc(zr))
!        write(6,*) 'harm vals minloc b0i: ',b0i(minloc(zr),:)
!        write(6,*) 'harm vals minloc bri: ',bri(minloc(zr),:)
!        write(6,*) 'harm vals minloc bthi: ',bthi(minloc(zr),:)
!        write(6,*) 'harm vals minloc bphi: ',bphi(minloc(zr),:)
!        write(6,*) 'harm vals interp bi: ',b(minloc(zr))%data(1),b(minloc(zr))%data(2), &
!             b(minloc(zr))%data(3),b(minloc(zr))%data(4)
!        write(6,*) 'harm vals minloc x0: ', ri(minloc(zr),:),thi(minloc(zr),:),&
!             phii(minloc(zr),:)
!        write(6,*) 'harm coords: ',zr(minloc(p)),theta(minloc(p)),zphi(minloc(p))
!        write(6,*) 'harm vals rhoi: ',rhoi(minloc(p),:)
        end subroutine harm3d_vals

        subroutine read_harm3d_inputs(ifile)
        character(len=20), intent(in) :: ifile
        open(unit=8,file=ifile,form='formatted',status='old')
        read(8,nml=harm)
        write(6,*) 'read: ',nt
        close(unit=8)
        end subroutine read_harm3d_inputs

        subroutine read_harm3d_data_header(nhead)
        integer, intent(in) :: nhead
        real(8), dimension(nhead) :: header
        real(8) :: tcur
        ! Read HARM header file
        ! JAD 11/24/2012 based on previous IDL codes
        ! header format from HARM dump.c
        write(6,*) 'read harm header: ',nhead
        open(unit=8,file=hfile,form='formatted',status='old')
        read(8,*) header
        nx1=header(2)
        nx2=header(3)
        nx3=header(4)
        startx1=header(5)
        startx2=header(6)
        startx3=header(7)
        dx1=header(8)
        dx2=header(9)
        dx3=header(10)
        asim=header(13)
        gam=header(14)
        h=header(nhead-1)
        tcur=header(1)
! is this true for HARM?
!        defcoord=header(19)
        close(unit=8)
        end subroutine read_harm3d_data_header

        subroutine read_harm3d_data_file(data_file,tcur,rho,p,u,b,mdot)
        character(len=100), intent(in) :: data_file
        character(len=100) :: data_file_app
        integer :: nhead,dlen
        integer :: rhopos,ppos,vpos,bpos,gdetpos
        integer, intent(in), optional :: mdot
        real(8), intent(out) :: tcur
        real(8), dimension(:), allocatable, intent(out) :: p,rho
        real(8), dimension(:), allocatable :: gdet, header, udotu, bdotu,bdotb, pmin
        type (four_vector), dimension(:), allocatable, intent(out) :: u,b
        type (four_vector), dimension(:), allocatable :: uks,bks
        real(8), dimension(:,:), allocatable :: grid, data, tmetric
        integer :: i,j, nelem, nelem2
        integer :: BINARY=1
! option to write an ASCII file read here back out as a binary
          ! Read HARM data from a file of given line length, 
          ! positions of desired variables
          ! JAD 11/24/2012 based on previous IDL codes
          ! Set grid parameters so we can generalize later:
!          nhead=21 ; dlen=72
          dlen=31; nhead=29
          allocate(header(nhead))
          rhopos=6; ppos=7; vpos=14; bpos=22; gdetpos=30
          allocate(grid(nx1*nx2*nx3,6)); allocate(p(nx1*nx2*nx3)); allocate(rho(nx1*nx2*nx3)); allocate(pmin(nx1*nx2*nx3))
          allocate(u(nx1*nx2*nx3)); allocate(b(nx1*nx2*nx3)); allocate(gdet(nx1*nx2*nx3))
          allocate(uks(nx1*nx2*nx3)); allocate(bks(nx1*nx2*nx3))
          pmin=1e-18
          write(6,*) 'data file: ',data_file
          if(BINARY.ne.1) then
             open(unit=8,file=data_file,form='formatted',status='old',action='read')
             allocate(data(dlen,nx3)); nelem=nx3
             write(6,*) 'read harm sizes', dlen, nx2, size(data)
             read(8,*) header
             tcur=header(1)
             deallocate(header)
             write(6,*) 'read harm past header: ',nx1,nx2,nx3,nelem,dlen
             do i=0,(nx1*nx2*nx3)/(nelem)-1
!                write(6,*) 'i: ',i
                read(8,*) data
!                write(6,*) 'harm data: ', data(:,1)
!             write(6,*) 'data loop: ',i,size(rho),size(p),size(b)
!             write(6,*) 'data vals: ',data(:,1),data(:,gdetpos+1)
                rho(1+i*nelem:(1+i)*nelem)=data(rhopos+1,:)
!                write(6,*) 'sizes rho: ',size(rho),size(data(rhopos+1,:))
                p(1+i*nelem:(1+i)*nelem)=data(ppos+1,:)
                b(1+i*nelem:(1+i)*nelem)%data(1)=(data(bpos+1,:))
                b(1+i*nelem:(1+i)*nelem)%data(2)=(data(bpos+2,:))
                b(1+i*nelem:(1+i)*nelem)%data(3)=(data(bpos+3,:))
                b(1+i*nelem:(1+i)*nelem)%data(4)=(data(bpos+4,:))
!                write(6,*) 'sizes b: ',size(b%data(4)),size(data(bpos+4,:))
                u(1+i*nelem:(1+i)*nelem)%data(1)=(data(vpos+1,:))
                u(1+i*nelem:(1+i)*nelem)%data(2)=(data(vpos+2,:))
                u(1+i*nelem:(1+i)*nelem)%data(3)=(data(vpos+3,:))
                u(1+i*nelem:(1+i)*nelem)%data(4)=(data(vpos+4,:))
!                write(6,*) 'sizes u: ',size(u%data(4)),size(data(vpos+4,:))
!               write(6,*) 'sizes: ',size(grid(1+i*nelem:(i+1)*nelem,:)),size(data(1:4,:))
 !               write(6,*) 'grid '
                grid(1+i*nelem:(1+i)*nelem,:)=transpose(data(1:6,:))
!                write(6,*) 'sizes grid: ',size(grid(1+i*nelem+j*nelem2:(1+i)*nelem+j*nelem2:(1+i)*nelem+(1+j)*nelem2,:)),size(data(1:6,:))
!            write(6,*) 'sizes: ',size(gdet),size(data(gdetpos+1,:))
                gdet(1+i*nelem:(1+i)*nelem)=(data(gdetpos+1,:))
                !             write(6,*) 'data loop'
             end do
!          write(6,*) 'after read loop'
             close(unit=8)
             deallocate(data)
             x1_arr=grid(:,1); x2_arr=grid(:,2); x3_arr=grid(:,3); r_arr=grid(:,4); th_arr=grid(:,5); ph_arr=grid(:,6)
          else
             data_file_app = trim(data_file) // '.bin'
             call read_harm3d_data(data_file_app,rho,p,u%data(1),u%data(2),u%data(3),u%data(4),b%data(1), &
                  b%data(2),b%data(3),b%data(4))
             call read_harm3d_grid_file()
             r_arr=exp(x1_arr); th_arr=x2_arr*pi; ph_arr=x3_arr
          endif
          p=merge(p,pmin,p.gt.pmin)
          write(6,*) "min vals p", minval(p)
          write(6,*) "min vals u", minval(u%data(1)), minval(u%data(2)), minval(u%data(3)), minval(u%data(4))
          write(6,*) "max vals u", maxval(u%data(1)), maxval(u%data(2)), maxval(u%data(3)), maxval(u%data(4))
          write(6,*) 'read harm grid sizes', size(x1_arr), size(x2_arr), size(x3_arr), size(r_arr), size(th_arr), size(ph_arr)
          write(6,*) 'read harm assign grid'
          deallocate(grid)
!          if (present(mdot)) then
             ! Calculate accretion rate in code units:
             !  nx1=n_elements(uniqx1) ; nx2=n_elements(uniqx2) ; nz=n_elements(uniqx3) 
!             dx2=uniqx2(2)-uniqx2(1) ; dx3=uniqx3(2)-uniqx3(1)
!             mdotarr=-1.*sum(sum(reform(gdet*rho*v(:,1),nx1,nx2,nz),3),2)*dx2*dx3
!          endif
          ! Transform velocities, magnetic fields from MKS to KS and then BL:
          write(6,*) 'read harm transform coords r ', minval(r_arr), maxval(r_arr), asim
          write(6,*) 'read harm transform coords th ',minval(x2_arr), maxval(x2_arr), h
          write(6,*) 'read harm transform coords phi',minval(x3_arr), maxval(x3_arr)
          write(6,*) 'read harm transform coords u ',minval(u%data(1)),minval(u%data(2))
!          write(6,*) 'read harm transform coords u ',u%data(1)
!          write(6,*) 'read harm transform coords u ',u%data(2)
!          write(6,*) 'read harm transform coords u ',u%data(3)
!          write(6,*) 'read harm transform coords u ',u%data(4)
          write(6,*) 'read harm transform coords u size ',size(u),h
          uks = umksh2uks(u,r_arr,x2_arr,h)
          write(6,*) 'after uks ',minval(uks%data(1))
          u = uks2ubl(uks,dble(r_arr),dble(asim))
          write(6,*) 'read harm transform coords b', minval(u%data(3)),maxval(u%data(2)),&
               minval(b%data(1)),maxval(b%data(4))
          write(6,*) "min vals u", minval(u%data(1)), minval(u%data(2)), minval(u%data(3)), minval(u%data(4))
          write(6,*) "max vals u", maxval(u%data(1)), maxval(u%data(2)), maxval(u%data(3)), maxval(u%data(4))
          bks = umksh2uks(b,r_arr,x2_arr,h)
          b   = uks2ubl(bks,dble(r_arr),dble(asim))
! test code...
          allocate(bdotu(n)); allocate(udotu(n)); allocate(bdotb(n))
          call assign_metric(u,transpose(kerr_metric(r_arr,th_arr,real(asim))))
          call assign_metric(b,transpose(kerr_metric(r_arr,th_arr,real(asim))))
          where(r_arr.gt.(1.+sqrt(1.-asim**2.)))
             udotu=abs(u*u+1.)
             bdotu=abs(u*b)
             bdotb=abs(b*b)
          elsewhere
             udotu=0.
             bdotu=0.
             bdotb=0.
          endwhere
          write(6,*) 'after transform', rho(1:4), p(1:4)
          write(6,*) 'after transform',p(5*nx1+23), rho(7*nx1+131), &
               maxval(abs(udotu)), maxval(abs(bdotu)),minval(bdotb)
          deallocate(udotu); deallocate(bdotu); deallocate(bdotb)
        end subroutine read_harm3d_data_file

        subroutine read_harm3d_grid_file()
!        integer, intent(in) :: nt
!        integer :: n
        open(unit=8,file=gfile,form='unformatted',status='old')
!        read(8) nx1,nx2,nx3
!        n=nx1*nx2*nx3
        ! allocate data
        !          call init_harm3d_data(n,n*nt)
        write(6,*) 'read grid: ',size(x1_arr),size(x2_arr),size(x3_arr)
        read(8) x1_arr
        read(8) x2_arr
        read(8) x3_arr
        close(unit=8)
        end subroutine read_harm3d_grid_file

        subroutine read_harm3d_data(dfile,rho,p,u0,vr,vth,vph,b0,br,bth,bph)
        integer :: nx, status
        real, dimension(:), allocatable :: data
        real(8), intent(inout), dimension(:) :: rho,p,u0,vr,vth,vph,b0,br,bth,bph
        real, dimension(size(rho),10) :: metric
        real, dimension(size(rho)) :: ui2
        character(len=100), intent(in) :: dfile
        open(unit=8,file=dfile,form='unformatted',status='old')
        read(8) nx
        if(nx.ne.10*nx1*nx2*nx3) then
           write(6,*) 'ERROR: INCORRECT DATA SIZE IN READ_HARM3D_DATA: ',nx,nx1,nx2,nx3
           return
        endif
        allocate(data(nx))
        read(8) data
        close(8)
        rho=data(1:nx1*nx2*nx3)
        p=data(nx1*nx2*nx3+1:2*nx1*nx2*nx3)
        u0=data(2*nx1*nx2*nx3+1:3*nx1*nx2*nx3)
        vr=data(3*nx1*nx2*nx3+1:4*nx1*nx2*nx3)
        vth=data(4*nx1*nx2*nx3+1:5*nx1*nx2*nx3)
        vph=data(5*nx1*nx2*nx3+1:6*nx1*nx2*nx3)
        b0=data(6*nx1*nx2*nx3+1:7*nx1*nx2*nx3)
        br=data(7*nx1*nx2*nx3+1:8*nx1*nx2*nx3)
        bth=data(8*nx1*nx2*nx3+1:9*nx1*nx2*nx3)
        bph=data(9*nx1*nx2*nx3+1:10*nx1*nx2*nx3)
        deallocate(data)
! calculate u0 from others (from IDL grtrans):
!        metric=kerr_metric(r_arr,th_arr,asim)
!        ui2=metric(:,1)+2.*metric(:,4)*vph+vr**2*metric(:,5)+metric(:,8)*vth**2+metric(:,10)*vph**2
!        u0=1./sqrt(-ui2)
        end subroutine read_harm3d_data

        subroutine initialize_harm3d_model(a,ifile,df,hf,gf,ntt,indft)
        real(kind=8), intent(in) :: a
        integer :: nx, status, nhead
        character(len=20), intent(in), optional :: ifile
        character(len=20) :: default_ifile='harm.in'
        character(len=100), intent(in), optional :: df,hf,gf
        integer, intent(in), optional :: ntt,indft
        nhead=29
! if all inputs are given, assign variables rather than reading from file
        if (present(df)) then
           dfile = df
           hfile = hf
           gfile = gf
           nt = ntt
           indf = indft
        else
           if (present(ifile)) then
              call read_harm3d_inputs(ifile)
           else
              call read_harm3d_inputs(default_ifile)
           endif
        endif
!        dfile='m87bl09rfp10xi5a998fluidvars.bin'
!        write(6,*) 'init harm'
        call read_harm3d_data_header(nhead)
        if (abs(asim-a).gt.1e-4) then 
           write(6,*) 'ERROR -- Different simulation and grtrans spin values!'
           return
        endif
!        write(6,*) 'read header', nx1, nx2
        n=nx1*nx2*nx3
        call init_harm3d_data(n,n*nt)
        write(6,*) 'init data', nt, indf, hfile, dfile
        call load_harm3d_data(nt)
        end subroutine initialize_harm3d_model

        subroutine load_harm3d_data(nt)
        real(8), dimension(:), allocatable :: rho,p
        real, dimension(:), allocatable :: vrl, vtl, vpl
        type (four_vector), dimension(:), allocatable :: u, b
        integer, intent(in) :: nt
        character(len=20) :: append
        character(len=100) :: data_file
        integer :: k
        real(8) :: tcur, tstep_test
        allocate(rho(n)); allocate(p(n))
        allocate(vrl(n)); allocate(vtl(n)); allocate(vpl(n))
        allocate(u(n)); allocate(b(n))
        do k=1,nt
           write(append, fmt='(I3.3)') indf-(k-1)
           data_file = trim(dfile) // append
!           write(6,*) 'data_file: ',indf-(k-1),append,data_file
           call read_harm3d_data_file(data_file,tcur,rho,p,u,b)
           t(k)=tcur
           ! transform velocities to LNRF frame
!           write(6,*) 'lnrf sizes: ',size(u%data(2)), size(r_arr), size(th_arr), size(vrl), size(vtl)
           call lnrf_frame(real(u%data(2)/u%data(1)),real(u%data(3)/u%data(1)), & 
                real(u%data(4)/u%data(1)),r_arr,asim,th_arr,vrl,vtl,vpl)
!           write(6,*) 'lnrf transform', size(b0_arr), size(b%data(1)), size(b0_arr((k-1)*n+1:k*n))
!           write(6,*) 'lnrf transform', size(vrl), size(vtl), size(vpl), size(u0_arr)
           ! now assign to data arrays
           rho_arr((k-1)*n+1:k*n)=rho
! conversion between internal energy and pressure
           p_arr((k-1)*n+1:k*n)=p*(gam-1.)
           b0_arr((k-1)*n+1:k*n)=b%data(1)
           br_arr((k-1)*n+1:k*n)=b%data(2)
           bth_arr((k-1)*n+1:k*n)=b%data(3)
           bph_arr((k-1)*n+1:k*n)=b%data(4)
!           write(6,*) 'v assign'
           u0_arr((k-1)*n+1:k*n)=u%data(1)
!           write(6,*) 'vrl assign'
           vrl_arr((k-1)*n+1:k*n)=real(vrl)
!           write(6,*) 'after assign', size(vtl_arr((k-1)*n+1:k*n)), size(vtl)
           vpl_arr((k-1)*n+1:k*n)=real(vpl)
!           write(6,*) 'after vpl', vtl(1), vtl(n), vtl
           vtl_arr((k-1)*n+1:k*n)=real(vtl)
!           write(6,*) 'after vtl'
!           write(6,*) 'assign'
        end do
!        tstep=10.
! Need to be able to do light curves even when nload (nt) = 1... Maybe read headers for first two files.
        if (nt.gt.1) then
! use default sim time step unless told otherwise
           tstep=t(1)-t(2)
           tstep_test=t(nt-1)-t(nt)
           if (abs((tstep-tstep_test)/tstep).gt.0.1) then 
              write(6,*) 'WARNING -- Time step changes over dumps!'
           endif
        endif
!        write(6,*) 'after loop', maxval(abs(u*u+1.))
        deallocate(rho); deallocate(p)
        deallocate(vrl); deallocate(vtl); deallocate(vpl)
        deallocate(u); deallocate(b)
!        write(6,*) 'read data file', a
        write(6,*) maxval(vrl_arr**2.+vtl_arr**2.+vpl_arr**2.)
!        write(6,*) 'maxmin temp: ',minval(p_arr/rho_arr*1.67e-24*9e20/1.38e-16), &
!             maxval(p_arr/rho_arr*1.67e-24*9e20/1.38e-16)
        end subroutine  load_harm3d_data

        subroutine advance_harm3d_timestep(dtobs)
        real(8), intent(in) :: dtobs
        integer :: nupdate
        nupdate = floor(dtobs / tstep)
        toffset = toffset+dtobs-tstep*nupdate
! after a couple stpes might be accumulating offset
        nupdate=nupdate+floor(toffset/tstep)
        toffset=toffset-floor(toffset/tstep)*tstep
        indf=indf+nupdate
        write(6,*) 'advance harm timestep: ',dtobs,tstep,toffset,indf,nupdate
        call update_harm3d_data(nupdate)
        end subroutine advance_harm3d_timestep

        subroutine update_harm3d_data(nupdate)
        integer, intent(in) :: nupdate
        integer :: nshift
        if(nupdate.gt.nt) then
! this is case where we just load all new data
           call load_harm3d_data(nt)
        else
! this is case where we shift existing data back and then load more
           nshift=n*nupdate+1
!           rho_arr=rho_arr(nshift:nshift-1:1)
!           allocate(temp_arr(n))
           rho_arr(nshift:n*nt)=rho_arr(1:n*nt-nshift)
           p_arr(nshift:n*nt)=p_arr(1:n*nt-nshift)
           br_arr(nshift:n*nt)=br_arr(1:n*nt-nshift)
           bth_arr(nshift:n*nt)=bth_arr(1:n*nt-nshift)
           bph_arr(nshift:n*nt)=bph_arr(1:n*nt-nshift)
           b0_arr(nshift:n*nt)=b0_arr(1:n*nt-nshift)
           vrl_arr(nshift:n*nt)=vrl_arr(1:n*nt-nshift)
           vtl_arr(nshift:n*nt)=vtl_arr(1:n*nt-nshift)
           vpl_arr(nshift:n*nt)=vpl_arr(1:n*nt-nshift)
           call load_harm3d_data(nupdate)
        end if
        end subroutine update_harm3d_data

        subroutine init_harm3d_data(nx,n)
        integer, intent(in) :: n,nx
!        integer, intent(in), optional :: aux,td,auxn
        allocate(rho_arr(n)); allocate(p_arr(n)); allocate(u0_arr(n))
        allocate(vrl_arr(n)); allocate(vtl_arr(n)); allocate(vpl_arr(n))
        allocate(b0_arr(n)); allocate(br_arr(n)); allocate(bth_arr(n))
        allocate(bph_arr(n))
        allocate(x1_arr(nx)); allocate(x2_arr(nx)); allocate(r_arr(nx))
        allocate(th_arr(nx)); allocate(t(nt)); allocate(x3_arr(nx))
        allocate(ph_arr(nx))
        end subroutine init_harm3d_data

        subroutine del_harm3d_data()
        deallocate(rho_arr); deallocate(p_arr); deallocate(u0_arr)
        deallocate(vrl_arr); deallocate(vtl_arr); deallocate(vpl_arr)
        deallocate(b0_arr); deallocate(br_arr); deallocate(bth_arr)
        deallocate(bph_arr); deallocate(x1_arr); deallocate(x2_arr)
        deallocate(r_arr); deallocate(th_arr); deallocate(t)
        deallocate(x3_arr)
        deallocate(ph_arr)
        end subroutine del_harm3d_data

      end module fluid_model_harm3d
