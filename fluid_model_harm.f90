   module fluid_model_harm

      use class_four_vector
      use interpolate, only: interp
      use kerr, only: kerr_metric, lnrf_frame, uks2ubl, bl2ks
      use phys_constants, only: pi
      use math, only: zbrent
      implicit none

      namelist /harm/  dfile, hfile, nt, indf

      character(len=40) :: dfile, hfile
      integer :: nx1, nx2, n, ndumps, nt, indf, dlen, nhead
      integer :: ndim=2
      integer, dimension(:), allocatable :: dumps
      real :: tstep=20., h, asim, dx1, dx2, gam, startx1, startx2, toffset=0.
      real, dimension(:), allocatable :: t
      real, dimension(:), allocatable :: x1_arr, x2_arr, r_arr, th_arr
      real, dimension(:), allocatable :: rho_arr, p_arr, u0_arr, vrl_arr, &
        vpl_arr, vtl_arr, temp_arr
      real, dimension(:), allocatable :: b0_arr, br_arr, bth_arr, bph_arr

      interface init_harm_data
        module procedure init_harm_data
      end interface

      interface del_harm_data
        module procedure del_harm_data
      end interface

      interface initialize_harm_model
        module procedure initialize_harm_model
      end interface

      interface transformbl2mksh
         module procedure transformbl2mksh
      end interface

      interface umksh2uks
         module procedure umksh2uks
      end interface

      interface harm_vals
        module procedure harm_vals
      end interface

      contains

        function findx2harm(x2,args) result(diff)
        real(kind=8), intent(in) :: x2
        real(kind=8), intent(in), dimension(:) :: args
        real(kind=8) :: h, theta, diff
        h=args(2); theta=args(1)
        diff=theta-pi*x2-.5*(1-h)*sin(2.*pi*x2)
        end function findx2harm


        subroutine transformbl2mksh(r,th,x1,x2,h)
        ! transform Boyer-Lindquist coordinates to modified Kerr-Schild coordinates used by HARM
        real, intent(in), dimension(:) :: r,th
        real, intent(out), dimension(size(r)) :: x1,x2
        real, intent(in) :: h
        real, dimension(2) :: args
        integer :: i
        x1=log(r)
        args=(/th(1),h/)
        do i=1,size(r)
           args(1)=th(i)
           x2(i)=zbrent(findx2harm,0d0,1d0,dble(args),1d-6)
        enddo
        end subroutine transformbl2mksh

        function umksh2uks(fum,r,x2,hin) result(uks)
          ! Converts a Kerr-Schild spherical 4-velocity to a Boyer-Lindquist one.
          ! JAD 11/10/2008
          ! Do simple transformation first:
          real, dimension(:), intent(in) :: r, x2
          type (four_vector), dimension(:), intent(in) :: fum
          type (four_vector), dimension(size(r)) :: uks
          real, intent(in), optional :: hin
          real :: hval
          real(kind=8), dimension(size(r)) :: ur, uth
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
          dthdx2=pi*(1.+(1.-hval)*cos(2.*pi*x2))
          uth=fum%data(3)*dble(dthdx2)
          uks=fum
          uks%data(2)=ur
          uks%data(3)=uth
        end function umksh2uks

        subroutine harm_vals(x0,a,rho,p,b,u,bmag)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real(kind=8), dimension(size(x0)) :: done
        real, dimension(size(x0)) :: x2,x1,zm,theta,fone, &
         vpl0,vrl0,vtl0,rd,td,rttd,zr,dzero, &
         vr0,vth0,vph0,bth,dummy,tt,ttd,zt,zphi,zpp
        real, dimension(nx1) :: uniqx1,uniqr
        real, dimension(nx2) :: uniqx2,uniqth
        real, dimension(size(x0),2**(ndim+1)) :: ppi,rhoi,vrli,vtli, &
         vpli,bri,bthi,bphi,b0i,u0i
!        integer, dimension(size(x0)*(2**ndim)) :: sindx
        real, dimension(size(x0)) :: dth
        integer, dimension(size(x0)*2*(2**ndim)) :: indx
        integer, dimension(size(x0)) :: lx1,lx2,ux1,ux2,x1l,x1u,x2l,x2u, &
         umax,one,tindx
        integer :: npts,i,maxtindx
        real, dimension(size(x0)), intent(out) :: rho,p,bmag
        type (four_Vector), intent(out), dimension(size(x0)) :: u,b
        ! Interpolates HARM data to input coordinates
        ! JAD 3/20/2009, fortran 11/12/2012
!        write(6,*) 'harm: ',size(x0)
        !nx1=floor(sqrt(real(size(x1_arr)))); nx2=nx1
        umax=nx2-1; one=1
        dzero=0d0; done=1d0; fone=1.
        uniqx2=x2_arr(1:nx2)
!        write(6,*) 'between'
        uniqx1=x1_arr(1:(nx2-1)*nx1+1:nx2)
        uniqr=exp(uniqx1)
        uniqth=pi*uniqx2+(1.-h)/2.*sin(2.*pi*uniqx2)
!        write(6,*) 'uniqx1: ',uniqx1
!        write(6,*) 'uniqx2: ',uniqx2
        npts=size(x0)
        theta=x0%data(3)
        zr=x0%data(2)
        zt=x0%data(1)
        tt=bl2ks(dble(zr),dble(zt),dble(a),1)-bl2ks(dble(zr(1)),0d0,dble(a),1)
        tt=-tt+1d-8
!        write(6,*) 'tt: ',minval(x0%data(1)),maxval(x0%data(1)),minval(tt),maxval(tt)
!        ztt=-ztt
!        tt=ztt-min(ztt)+tabs0
        zpp=x0%data(4)
        zphi=bl2ks(dble(zr),dble(zpp),dble(a))
        call transformbl2mksh(zr,theta,x1,x2,h)
!        write(6,*) 'transform: ',maxval(zr), minval(zr), maxval(x1), minval(x1)
!        write(6,*) 'transform: ',maxval(theta), minval(theta), minval(x2), maxval(x2)
        lx1=floor((x1-uniqx1(1))/(uniqx1(nx1)-uniqx1(1))*(nx1-1))+1
        ux1=lx1+1
        lx2=floor((x2-uniqx2(1))/(uniqx2(nx2)-uniqx2(1))*(nx2-1))+1
        ux2=lx2+1
!        lx2=merge(merge(lx2,one,lx2.ge.1),umax,lx2.le.(nx2-1))
        lx2=merge(lx2,one,lx2.ge.1)
        ux2=merge(ux2,nx2,ux2.le.nx2)
! Deal with poles
        where(ux2.ne.lx2)
           dth=uniqth(ux2)-uniqth(lx2)
        elsewhere
           dth=uniqth(1)
        endwhere
        td=abs(theta-uniqth(lx2))/dth
        rd=(zr-uniqr(lx1))/(uniqr(ux1)-uniqr(lx1))
! When the geodesic is inside the innermost zone center located outside the horizon, 
! use nearest neighbor rather than interpolate
        where(uniqr(lx1).le.(1.+sqrt(1.-a**2.)))
           rd=1.
        endwhere
!        write(6,*) 'rd td: ',minval(rd),maxval(rd),minval(td),maxval(td)
!        write(6,*) 'ux lx: ',minval(lx1),maxval(ux1),minval(lx2),maxval(ux2)
        ! th is fastest changing index
        x2l=lx2-1 ; x2u=ux2-1
        x1l=(lx1-1)*nx2 ; x1u=(ux1-1)*nx2
! Indices for nearest neighbors in 2D
!        write(6,*) 'sindx', size(sindx)
!        sindx=(/(/x1l+x2l+1/),(/x1l+x2u+1/),(/x1u+x2l+1/),(/x1u+x2u+1/)/)
!        write(6,*) 'sindx', maxval(sindx), sindx(npts+1), x1l(1)+x2u(1)+1, size(sindx)
! Now find timestep indices:
        maxtindx=nt-2
        where(floor(tt/tstep).le.maxtindx)
           tindx=floor(tt/tstep)
           ttd=(tt-tindx*tstep)/tstep
        elsewhere
           tindx=maxtindx+1
           ttd=0.
        endwhere
!        write(6,*) 'tindx: ',minval(tindx),maxval(tindx)
! Create index array across times and 2D nearest neighbors (first half at first time, second half at second time):
! Need to make sure there aren't integer problems when these indices get very large
!        write(6,*) 'indx', size(indx), 2**ndim, maxval(tindx), minval(tindx), tstep
        indx=(/(/x1l+x2l+1+n*tindx/),(/x1l+x2u+1+n*tindx/),(/x1u+x2l+1+n*tindx/),(/x1u+x2u+1+n*tindx/), &
             (/x1l+x2l+1+n*(tindx+1)/),(/x1l+x2u+1+n*(tindx+1)/),(/x1u+x2l+1+n*(tindx+1)/),(/x1u+x2u+1+n*(tindx+1)/)/)
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
!        write(6,*) 'after reshape', minval(ttd), maxval(ttd)
!        rttd=ttd
        rttd=0.
        rho=merge(interp(rhoi,rttd,rd,td),dzero,x1.gt.uniqx1(1))
        p=merge(interp(ppi,rttd,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'rho: ', rho(1),rhoi(1,:),rd(1),td(1)
!        write(6,*) 'rho: ',zr(1),x0(1)%data,theta(1),x1(1),x2(1)
!        write(6,*) 'rho: ',lx1(1),lx2(1),ux1(1),ux2(1)
!        write(6,*) 'rho: ',size(indx),size(x1l)*8
!        write(6,*) 'rho: ',size(((/x1l+x2l+1+n*tindx/),(/x1l+x2u+1+n*tindx/),(/!x1u+x\
!2l+1+n*tindx/),(/x1u+x2u+1+n*tindx/), &                                                
!             (/x1l+x2l+1+n*(tindx+1)/),(/x1l+x2u+1+n*(tindx+1)/),(/x1u+x2l+1+n*!(tind\
!x+1)/), &                                                                              
!             (/x1u+x2u+1+n*(tindx+1)/)/))                                             
!        write(6,*) 'rho: ',rttd(1)
!        write(6,*) 'rho: ',rttd(1),maxtindx,tindx(1),tindx(2),tt(1),tstep, &
!             floor(tt(1)/tstep), floor(tt(2)/tstep),tt(2)
!        write(6,*) 'rho: ',indx((/1,size(x1u)+1,2*size(x1u)+1,3*size(x1u)+1, &
!             4*size(x1u)+1,5*size(x1u)+1,6*size(x1u)+1,7*size(x1u)+1/))
!        write(6,*) 'rho: ',rho_arr(indx((/1,size(x1u)+1,2*size(x1u)+1,3*size(x1u)+1, &
!             4*size(x1u)+1,5*size(x1u)+1,6*size(x1u)+1,7*size(x1u)+1/)))
!        write(6,*) 'rho: ', rho, p
        vrl0=merge(interp(vrli,rttd,rd,td),dzero,x1.gt.uniqx1(1))
        vtl0=merge(interp(vtli,rttd,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'v'
        b%data(1)=merge(interp(b0i,rttd,rd,td),fone,x1.gt.uniqx1(1))
        b%data(2)=merge(interp(bri,rttd,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'b'
        b%data(3)=merge(interp(bthi,rttd,rd,td),fone,x1.gt.uniqx1(1))
        b%data(4)=merge(interp(bphi,rttd,rd,td),fone,x1.gt.uniqx1(1))
!        write(6,*) 'u'
        u%data(1)=merge(dble(interp(u0i,rttd,rd,td)),done,x1.gt.uniqx1(1))
        vpl0=merge(interp(vpli,rttd,rd,td),dzero,x1.gt.uniqx1(1))
  !      write(6,*) 'after assign'
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
!        write(6,*) 'harm vals r: ',zr
!        write(6,*) 'harm vals rho: ',rho
!        write(6,*) 'harm vals rd: ',rd
!        write(6,*) 'harm vals ux: ',nx1,ux1
!        write(6,*) 'harm vals minloc: ',minloc(rho),rd(minloc(rho)), &
!             td(minloc(rho)),x0(minloc(rho))%data(3),uniqth(lx2(minloc(rho))), &
!             lx2(minloc(rho))
!        write(6,*) 'harm vals rhoi: ',rhoi(minloc(rho),:)
        end subroutine harm_vals

        subroutine read_harm_inputs(ifile)
        character(len=20), intent(in) :: ifile
        open(unit=8,file=ifile,form='formatted',status='old')
        read(8,nml=harm)
        write(6,*) 'read: ',nt
        close(unit=8)
        end subroutine read_harm_inputs

        subroutine read_harm_data_header(nhead)
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
        startx1=header(4)
        startx2=header(5)
        dx1=header(6)
        dx2=header(7)
        asim=header(10)
        gam=header(11)
        h=header(nhead-1)
        tcur=header(1)
! is this true for HARM?
!        defcoord=header(19)
        close(unit=8)
        end subroutine read_harm_data_header

        subroutine read_harm_data_file(data_file,tcur,rho,p,u,b,mdot)
        character(len=40), intent(in) :: data_file
        integer :: nhead,dlen
        integer :: rhopos,ppos,vpos,bpos,gdetpos
        integer, intent(in), optional :: mdot
        real(8), intent(out) :: tcur
        real(8), dimension(:), allocatable, intent(out) :: p,rho
        real(8), dimension(:), allocatable :: gdet, header, udotu, bdotu
        type (four_vector), dimension(:), allocatable, intent(out) :: u,b
        type (four_vector), dimension(:), allocatable :: uks,bks
        real(8), dimension(:,:), allocatable :: grid, data, tmetric
        integer :: i, nelem
          ! Read HARM data from a file of given line length, 
          ! positions of desired variables
          ! JAD 11/24/2012 based on previous IDL codes
          ! Set grid parameters so we can generalize later:
!          nhead=21 ; dlen=72
          dlen=34; nhead=26
          allocate(header(nhead))
          rhopos=4; ppos=5; vpos=13; bpos=21; gdetpos=33
          write(6,*) 'data file: ',data_file
          open(unit=8,file=data_file,form='formatted',status='old',action='read')
          allocate(data(dlen,nx2)); nelem=nx2
          allocate(grid(nx1*nx2,4)); allocate(p(nx1*nx2)); allocate(rho(nx1*nx2))
          allocate(u(nx1*nx2)); allocate(b(nx1*nx2)); allocate(gdet(nx1*nx2))
          allocate(uks(nx1*nx2)); allocate(bks(nx1*nx2))
          write(6,*) 'read harm sizes', dlen, nx2, size(data)
          read(8,*) header
          tcur=header(1)
          deallocate(header)
          write(6,*) 'read harm past header: ',nx1,nx2,nelem,dlen
          do i=0,(nx1*nx2)/nelem-1
             read(8,*) data
!             if (i.eq.0) write(6,*) 'harm data: ', data(:,1)
!             write(6,*) 'data loop: ',i,size(rho),size(p),size(b)
!             write(6,*) 'data vals: ',data(:,1),data(:,gdetpos+1)
             rho(1+i*nelem:(i+1)*nelem)=data(rhopos+1,:)
!             write(6,*) 'sizes: ',size(rho(1+i*nelem:(i+1)*nelem)),size(data(rhopos+1,:))
             p(1+i*nelem:(i+1)*nelem)=data(ppos+1,:)
             b(1+i*nelem:(i+1)*nelem)%data(1)=(data(bpos+1,:))
             b(1+i*nelem:(i+1)*nelem)%data(2)=(data(bpos+2,:))
             b(1+i*nelem:(i+1)*nelem)%data(3)=(data(bpos+3,:))
             b(1+i*nelem:(i+1)*nelem)%data(4)=(data(bpos+4,:))
!             write(6,*) 'sizes: ',size(b(1+i*nelem:(i+1)*nelem)%data(4)),size(data(bpos+4,:))
             u(1+i*nelem:(i+1)*nelem)%data(1)=(data(vpos+1,:))
             u(1+i*nelem:(i+1)*nelem)%data(2)=(data(vpos+2,:))
             u(1+i*nelem:(i+1)*nelem)%data(3)=(data(vpos+3,:))
             u(1+i*nelem:(i+1)*nelem)%data(4)=(data(vpos+4,:))
!             write(6,*) 'sizes: ',size(u(1+i*nelem:(i+1)*nelem)%data(4)),size(data(vpos+4,:))
!             write(6,*) 'sizes: ',size(grid(1+i*nelem:(i+1)*nelem,:)),size(data(1:4,:))
             grid(1+i*nelem:(i+1)*nelem,:)=transpose(data(1:4,:))
!             write(6,*) 'sizes: ',size(gdet(1+i*nelem:(i+1)*nelem)),size(data(gdetpos+1,:))
             gdet(1+i*nelem:(i+1)*nelem)=(data(gdetpos+1,:))
!             write(6,*) 'data loop'
          end do
!          write(6,*) 'after read loop'
          close(unit=8)
          deallocate(data)
          write(6,*) 'read harm grid sizes', size(x1_arr), size(x2_arr), size(r_arr), size(th_arr)
          x1_arr=grid(:,1); x2_arr=grid(:,2); r_arr=grid(:,3); th_arr=grid(:,4)
          write(6,*) 'read harm assign grid'
          deallocate(grid)
          if (present(mdot)) then
             ! Calculate accretion rate in code units:
             !  nx1=n_elements(uniqx1) ; nx2=n_elements(uniqx2) ; nz=n_elements(uniqx3) 
!             dx2=uniqx2(2)-uniqx2(1) ; dx3=uniqx3(2)-uniqx3(1)
!             mdotarr=-1.*sum(sum(reform(gdet*rho*v(:,1),nx1,nx2,nz),3),2)*dx2*dx3
          endif
          ! Transform velocities, magnetic fields from MKS to KS and then BL:
          write(6,*) 'read harm transform coords u ', minval(r_arr), maxval(r_arr), asim
          write(6,*) 'read harm transform coords u ', minval(th_arr), maxval(th_arr), asim
!          write(6,*) 'read harm transform coords u ',minval(x2_arr), maxval(x2_arr), h
!          write(6,*) 'read harm transform coords u ',minval(u%data(1)),maxval(u%data(1))
!          write(6,*) 'read harm transform coords b ',minval(b%data(1)),maxval(b%data(1))

!          write(6,*) 'read harm transform coords u ',u%data(1)
!          write(6,*) 'read harm transform coords u ',u%data(2)
!          write(6,*) 'read harm transform coords u ',u%data(3)
!          write(6,*) 'read harm transform coords u ',u%data(4)
!          write(6,*) 'read harm transform coords u size ',size(u),h
          uks = umksh2uks(u,r_arr,x2_arr,h)
          write(6,*) 'after uks ',minval(uks%data(1))
          u = uks2ubl(uks,dble(r_arr),dble(asim))
          write(6,*) 'read harm transform cords b', minval(u%data(3)),maxval(u%data(2)),&
               minval(b%data(1)),maxval(b%data(4))
          bks = umksh2uks(b,r_arr,x2_arr,h)
          b   = uks2ubl(bks,dble(r_arr),dble(asim))

          write(6,*) 'read harm transform coords u ',maxval(u%data(1))
          write(6,*) 'read harm transform coords u ',maxval(u%data(2))
          write(6,*) 'read harm transform coords u ',maxval(u%data(3))
          write(6,*) 'read harm transform coords u ',maxval(u%data(4))

! test code...
          allocate(bdotu(n)); allocate(udotu(n))
          call assign_metric(u,transpose(kerr_metric(r_arr,th_arr,real(asim))))
          call assign_metric(b,transpose(kerr_metric(r_arr,th_arr,real(asim))))
          where(r_arr.gt.(1.+sqrt(1.-asim**2.)))
             udotu=abs(u*u+1.)
             bdotu=abs(u*b)
          elsewhere
             udotu=0.
             bdotu=0.
          endwhere
          write(6,*) 'after transform', rho(1:4), p(1:4)
          write(6,*) 'after transform',p(5*nx1+23), rho(7*nx1+131), &
               maxval(udotu), maxval(bdotu)
          deallocate(udotu); deallocate(bdotu)
        end subroutine read_harm_data_file


        subroutine initialize_harm_model(a,ifile,df,hf,ntt,indft)
        real(kind=8), intent(in) :: a
        integer :: nx, status, nhead
        character(len=20), intent(in), optional :: ifile
        character(len=20) :: default_ifile='harm.in'
        character(len=40), intent(in), optional :: df,hf
        integer, intent(in), optional :: ntt,indft
        nhead=26
! if all inputs are given, assign variables rather than reading from file
        if (present(df)) then
           dfile = df
           hfile = hf
           nt = ntt
           indf = indft
        else
           if (present(ifile)) then
              call read_harm_inputs(ifile)
           else
              call read_harm_inputs(default_ifile)
           endif
        endif
!        dfile='m87bl09rfp10xi5a998fluidvars.bin'
!        write(6,*) 'init harm'
        call read_harm_data_header(nhead)
        if (abs(asim-a).gt.1e-4) then 
           write(6,*) 'ERROR -- Different simulation and grtrans spin values!'
           return
        endif
!        write(6,*) 'read header', nx1, nx2
        n=nx1*nx2
        call init_harm_data(n,n*nt)
!        write(6,*) 'init data', nt, indf, hfile, dfile
        call load_harm_data(nt)
        end subroutine initialize_harm_model

        subroutine load_harm_data(nt)
        real(8), dimension(:), allocatable :: rho,p
        real, dimension(:), allocatable :: vrl, vtl, vpl
        type (four_vector), dimension(:), allocatable :: u, b
        integer, intent(in) :: nt
        character(len=20) :: append
        character(len=40) :: data_file
        integer :: k
        real(8) :: tcur, tstep_test
        allocate(rho(n)); allocate(p(n))
        allocate(vrl(n)); allocate(vtl(n)); allocate(vpl(n))
        allocate(u(n)); allocate(b(n))

        !AC loop is backward?
        do k=1,nt
           write(append, fmt='(I3.3)') indf-(k-1)
           data_file = trim(dfile) // append
!           write(6,*) 'data_file: ',indf-(k-1),append,data_file
           call read_harm_data_file(data_file,tcur,rho,p,u,b)
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
        end subroutine  load_harm_data

        subroutine advance_harm_timestep(dtobs)
        real(8), intent(in) :: dtobs
        integer :: nupdate
        nupdate = floor(dtobs / tstep)
        toffset = toffset+dtobs-tstep*nupdate
! after a couple stpes might be accumulating offset
        nupdate=nupdate+floor(toffset/tstep)
        toffset=toffset-floor(toffset/tstep)*tstep
        indf=indf+nupdate
        write(6,*) 'advance harm timestep: ',dtobs,tstep,toffset,indf,nupdate
        call update_harm_data(nupdate)
        end subroutine advance_harm_timestep

        subroutine update_harm_data(nupdate)
        integer, intent(in) :: nupdate
        integer :: nshift
        if(nupdate.gt.nt) then
! this is case where we just load all new data
           call load_harm_data(nt)
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
           call load_harm_data(nupdate)
        end if
        end subroutine update_harm_data

        subroutine init_harm_data(nx,n)
        integer, intent(in) :: n,nx
!        integer, intent(in), optional :: aux,td,auxn
        allocate(rho_arr(n)); allocate(p_arr(n)); allocate(u0_arr(n))
        allocate(vrl_arr(n)); allocate(vtl_arr(n)); allocate(vpl_arr(n))
        allocate(b0_arr(n)); allocate(br_arr(n)); allocate(bth_arr(n))
        allocate(bph_arr(n))
        allocate(x1_arr(nx)); allocate(x2_arr(nx)); allocate(r_arr(nx))
        allocate(th_arr(nx)); allocate(t(nt))
        end subroutine init_harm_data

        subroutine del_harm_data()
        deallocate(rho_arr); deallocate(p_arr); deallocate(u0_arr)
        deallocate(vrl_arr); deallocate(vtl_arr); deallocate(vpl_arr)
        deallocate(b0_arr); deallocate(br_arr); deallocate(bth_arr)
        deallocate(bph_arr); deallocate(x1_arr); deallocate(x2_arr)
        deallocate(r_arr); deallocate(th_arr); deallocate(t)
        end subroutine del_harm_data

      end module fluid_model_harm
