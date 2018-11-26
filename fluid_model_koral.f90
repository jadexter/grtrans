   !AC -- read KORAL files 2/2/17
   !   -- modified from fluid_model_harm.f90
   module fluid_model_koral

      use class_four_vector
      use interpolate, only: interp
      use kerr, only: kerr_metric, lnrf_frame, uks2ubl, bl2ks
      use phys_constants, only: pi
      use math, only: zbrent
      implicit none

      !AC koral will read modified harm .in files for now
      !AC should change to koral namelist
      namelist /harm/  dfile, hfile, nt, indf

      logical :: doingkoralnth
      character(len=100) :: dfile, hfile
      integer :: nx1, nx2, n, ndumps, nt, indf, dlen, nhead
      integer :: nrelbin=0
      integer :: ndim=2
      integer :: minpolecell=4
      integer, dimension(:), allocatable :: dumps
      real :: tstep=20., h, r0, aa, bb, pp, asim, dx1, dx2, gam, startx1, startx2, toffset=0.
      real :: relgammamin=1., relgammamax=1.
      real(8) :: scalefac
      real, dimension(:), allocatable :: t
      real, dimension(:), allocatable :: x1_arr, x2_arr, r_arr, th_arr
      real, dimension(:), allocatable :: rho_arr, p_arr, u0_arr, vrl_arr, &
           vpl_arr, vtl_arr, temp_arr
      real, dimension(:,:), allocatable :: nth_arr
      real, dimension(:), allocatable :: b0_arr, br_arr, bth_arr, bph_arr


      interface init_koral_data
        module procedure init_koral_data
      end interface

      interface del_koral_data
        module procedure del_koral_data
      end interface

      interface initialize_koral_model
        module procedure initialize_koral_model
      end interface

      interface transformbl2mksh3
         module procedure transformbl2mksh3
      end interface

      interface transformmksh32bl
         module procedure transformmksh32bl
      end interface

      !interface umksh2uks
      !   module procedure umksh2uks
      !end interface

      interface koral_vals
        module procedure koral_vals
      end interface

      contains

        subroutine transformbl2mksh3(th,r,x2,h,a,b,p)
        ! transform Boyer-Lindquist theta coordinates to modified Kerr-Schild x2 at a corresponding array of r 
        real, intent(in), dimension(:) :: th,r
        real, intent(out), dimension(size(th)) :: x2
        real, intent(in) :: h,a,b,p

        x2=0.5*(1 + ((r**p)/(h*pi)) * (atan(tan(0.5*h*pi)*(1-2*th/pi))/((b-a)*(2**p)+(a-0.5)*(r**p))))
        end subroutine transformbl2mksh3


        subroutine transformmksh32bl(x2,r,th,h,a,b,p)
        ! transform modified Kerr-Schild x2 to Boyer-Lindquist theta at a corresponding array of r=r0+exp(x1)
        real, intent(in), dimension(:) :: x2,r
        real, intent(out), dimension(size(x2)) :: th
        real, intent(in) :: h,a,b,p
        
        th = 0.5*pi*(1 + tan(h*pi*(-0.5+x2+(1-2*x2)*(a+(2**p)*(b-a)/(r**p))))/tan(0.5*h*pi))
        !write(6,*) 'x2len: ', size(x2)
        end subroutine transformmksh32bl

        subroutine koral_vals(x0,a,rho,p,b,u,bmag,nth)
        type (four_Vector), intent(in), dimension(:) :: x0
        real, intent(in) :: a
        real(kind=8), dimension(size(x0)) :: done
        real, dimension(size(x0)) :: x2,x1,zm,theta,fone, &
         vpl0,vrl0,vtl0,rd,td,rttd,zr,dzero, &
         vr0,vth0,vph0,bth,dummy,tt,ttd,zt,zphi,zpp
        real, dimension(nx1) :: uniqx1,uniqr
        real, dimension(nx2) :: uniqx2
        real, dimension(size(x0),2**(ndim+1)) :: ppi,rhoi,vrli,vtli, &
             vpli,bri,bthi,bphi,b0i,u0i, nth_bini

        real, dimension(size(x0)) :: dth,thu,thl,thone
        integer, dimension(size(x0)*2*(2**ndim)) :: indx
        integer, dimension(size(x0)) :: lx1,lx2,ux1,ux2,x1l,x1u,x2l,x2u, &
         umax,one,tindx,lx2l,lx2u,ux2l,ux2u,low,high
        real :: x2mintrust, x2maxtrust
        integer :: npts,i,j,maxtindx,ie
        real, dimension(size(x0)) :: b1tmp,b2tmp,b3tmp,b4tmp,u1tmp,ntmp
        real, dimension(size(x0)), intent(out) :: rho,p,bmag
        real, dimension(size(x0),nrelbin), intent(out) :: nth
        type (four_Vector), intent(out), dimension(size(x0)) :: u,b
        
        ! Interpolates HARM data to input coordinates
        ! JAD 3/20/2009, fortran 11/12/2012
       
        umax=nx2-1; one=1
        dzero=0d0; done=1d0; fone=1
         
        uniqx2=x2_arr(1:nx2)
        uniqx1=x1_arr(1:(nx2-1)*nx1+1:nx2)
        


        npts=size(x0)
        theta=x0%data(3)
        zr=x0%data(2)
        zt=x0%data(1)
        tt=bl2ks(dble(zr),dble(zt),dble(a),1) - bl2ks(dble(zr(1)),0d0,dble(a),1)
        tt=-tt+1d-8

        x1=log(dble(zr)-r0)
        call transformbl2mksh3(theta,zr,x2,h,aa,bb,pp)

        lx1=floor((x1-uniqx1(1))/(uniqx1(nx1)-uniqx1(1))*(nx1-1))+1
        ux1=lx1+1
        lx1=merge(lx1,one,lx1.ge.1)
        ux1=merge(ux1,nx1,ux1.le.nx1)

        lx2=floor((x2-uniqx2(1))/(uniqx2(nx2)-uniqx2(1))*(nx2-1))+1
        ux2=lx2+1
        lx2=merge(lx2,one,lx2.ge.1)
        ux2=merge(ux2,nx2,ux2.le.nx2)

!AC do x1,x2 distances instead of r, th distances -- ok??
! r-distance
        uniqr=r0+exp(uniqx1)
        rd=(zr-uniqr(lx1))/(uniqr(ux1)-uniqr(lx1))

! When the geodesic is inside the innermost zone center located outside the horizon, 
! use nearest neighbor rather than interpolate
        where((uniqr(lx1).le.(1.+sqrt(1.-a**2.))).or.(lx1.eq.1))
           rd=1.
        endwhere

! theta-distance
! Deal with poles 

        !AC for interpolation, choose rl for theta distance? 
        call transformmksh32bl(uniqx2(lx2),uniqr(lx1),thl,h,aa,bb,pp)
        call transformmksh32bl(uniqx2(ux2),uniqr(lx1),thu,h,aa,bb,pp)
        call transformmksh32bl(uniqx2(one),uniqr(lx1),thone,h,aa,bb,pp)
        
        !write(6,*) minval(thone),maxval(thone)
        !write(6,*) minval(thu-thl),maxval(thu-thl)
        !write(6,*) minval(thl),maxval(thl)
        !write(6,*) minval(thu),maxval(thu)

        where(ux2.ne.lx2)
           td=abs((theta-thl)/(thu-thl))
        elsewhere
           td=abs((theta-thl)/thone)
        endwhere
        
        !where(ux2.ne.lx2)
        !   td=abs((x2-uniqx2(lx2))/(uniqx2(ux2)-uniqx2(lx2)))
        !elsewhere
        !   td=(x2-uniqx2(lx2))/uniqx2(1)
        !endwhere 
        !write(6,*) size(td), minval(td), maxval(td)

! th is fastest changing index
        x2l=lx2-1 
        x2u=ux2-1
        x1l=(lx1-1)*nx2 
        x1u=(ux1-1)*nx2

! Indices for nearest neighbors in 2D
! Now find timestep indices:
        maxtindx=nt-2
        where(floor(tt/tstep).le.maxtindx)
           tindx=floor(tt/tstep)
           ttd=(tt-tindx*tstep)/tstep
        elsewhere
           tindx=maxtindx+1
           ttd=0.
        endwhere

! Create index array across times and 2D nearest neighbors (first half at first time, second half at second time):
! Need to make sure there aren't integer problems when these indices get very large
        indx=(/(/x1l+x2l+1+n*tindx/),(/x1l+x2u+1+n*tindx/),(/x1u+x2l+1+n*tindx/),(/x1u+x2u+1+n*tindx/), &
             (/x1l+x2l+1+n*(tindx+1)/),(/x1l+x2u+1+n*(tindx+1)/),(/x1u+x2l+1+n*(tindx+1)/),(/x1u+x2u+1+n*(tindx+1)/)/)
! keep indx from going out of bounds when loading one time slice
        where(indx.gt.size(rho_arr))
           indx=indx-n
        endwhere
     
        !write (6,*) 'index :', shape(indx) ,2**(ndim+1), shape(vpl_arr)
        !write (6,*) 'lx1 min/max:', minval(lx1) , maxval(lx1)
        !write (6,*) 'ux1 min/max:', minval(ux1) , maxval(ux1)
        !write (6,*) 'lx2 min/max:', minval(lx2) , maxval(lx2)
        !write (6,*) 'ux2 min/max:', minval(ux2) , maxval(ux2)
        !write (6,*) 'index min/max:', minval(indx) , maxval(indx)
        !write (6,*) 'vtl', shape(vtl_arr(indx))
        !write (6,*) 'vpl', shape(vpl_arr(indx))
      
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

        !interpolate data
        rttd=0.
        rho=interp(rhoi,rttd,rd,td)
        p=interp(ppi,rttd,rd,td)
        vrl0=interp(vrli,rttd,rd,td)
        vtl0=interp(vtli,rttd,rd,td)
        b1tmp=interp(b0i,rttd,rd,td)
        b2tmp=interp(bri,rttd,rd,td)
        b3tmp=interp(bthi,rttd,rd,td)
        b4tmp=interp(bphi,rttd,rd,td)
        u1tmp=interp(u0i,rttd,rd,td)
        vpl0=interp(vpli,rttd,rd,td)

        !mask data outside trusted region
        rho=merge(rho,dzero,x1.gt.uniqx1(1))
        p=merge(p,fone,x1.gt.uniqx1(1)) 
        vrl0=merge(vrl0,dzero,x1.gt.uniqx1(1))
        vtl0=merge(vtl0,fone,x1.gt.uniqx1(1))
        b1tmp=merge(b1tmp,fone,x1.gt.uniqx1(1))
        b2tmp=merge(b2tmp,fone,x1.gt.uniqx1(1))
        b3tmp=merge(b3tmp,fone,x1.gt.uniqx1(1))
        b4tmp=merge(b4tmp,fone,x1.gt.uniqx1(1))
        u1tmp=merge(dble(u1tmp),done,x1.gt.uniqx1(1))
        vpl0=merge(vpl0,dzero,x1.gt.uniqx1(1))

        !AC minimum polar cell - change!
        !AC this is the minimum polar coordinate x2 we trust. 
        x2mintrust = uniqx2(minpolecell)
        x2maxtrust = uniqx2(size(uniqx2)-(minpolecell-1))

        rho=merge(rho,dzero,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        p=merge(p,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust)) 
        vrl0=merge(vrl0,dzero,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        vtl0=merge(vtl0,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        b%data(1)=merge(b1tmp,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        b%data(2)=merge(b2tmp,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        b%data(3)=merge(b3tmp,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        b%data(4)=merge(b4tmp,fone,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        u%data(1)=merge(dble(u1tmp),done,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))
        vpl0=merge(vpl0,dzero,(x2.gt.x2mintrust).and.(x2.lt.x2maxtrust))

        !AC do we have to do this in a loop? 
        if(doingkoralnth) then
           do ie=1, nrelbin
              nth_bini=reshape(nth_arr(indx,ie),(/npts,2**(ndim+1)/))
              ntmp=interp(nth_bini,rttd,rd,td)
              ntmp=merge(ntmp,dzero,x1.gt.uniqx1(1))
              nth(:,ie)=merge(ntmp,dzero,x2.gt.x2mintrust)
           end do
        end if
        
        ! Protect azimuthal velocities at poles
        ! Compute magnitude of interpolated b-field and force b^2 > 0 (need to look into numerical issues here):
        call assign_metric(b,transpose(kerr_metric(zr, real(x0%data(3)),a)))
        bmag=b*b; bmag=merge(bmag,dzero,bmag.ge.0d0)
        bmag=sqrt(bmag)

        call lnrf_frame(vrl0,vtl0,vpl0,zr,a,real(x0%data(3)),vr0,vth0,vph0,1)

        u%data(2)=u%data(1)*dble(vr0)
        u%data(3)=u%data(1)*dble(vth0)
        u%data(4)=u%data(1)*dble(vph0)

        call assign_metric(u,transpose(kerr_metric(zr,real(x0%data(3)) ,a)))

        end subroutine koral_vals

        subroutine read_koral_inputs(ifile)
        character(len=20), intent(in) :: ifile
        open(unit=8,file=ifile,form='formatted',status='old')
        read(8,nml=harm) !AC koral will share the harm input file structure for now 
        write(6,*) 'read: ',nt
        close(unit=8)
        end subroutine read_koral_inputs

        !AC MODIFY
        subroutine read_koral_data_header(nhead)
        integer, intent(in) :: nhead
        real(8), dimension(nhead) :: header
        real(8), dimension(3) :: header2
        real(8) :: tcur
     
        write(6,*) 'read header: ', nhead
        open(unit=8,file=hfile,form='formatted',status='old')
        read(8,*) header
        
        tcur=header(1)
        nx1=header(2)
        nx2=header(3)
        asim=header(4)
        !mbh=header(5)
        r0=header(6)
        h=header(7)
        aa=header(8)
        bb=header(9)
        pp=header(10)

        !AC nonthermal
        if(doingkoralnth) then
           read(8,*) header2
           nrelbin=header2(1)
           relgammamin=header2(2)
           relgammamax=header2(3)
           write(6,*) 'read nth header vals ', nrelbin, relgammamin, relgammamax

        end if
        
        write(6,*) 'read header vals: ', asim, r0, h, aa, bb, pp
        close(unit=8)
        end subroutine read_koral_data_header

        subroutine read_koral_data_file(data_file,tcur,rho,p,u,b,nnth)
        character(len=100), intent(in) :: data_file
        integer :: nhead, nhead2, dlen
        integer :: rhopos,ppos,vpos,bpos,gdetpos,gridpos,nthpos
        real(8), intent(out) :: tcur
        real(8), dimension(:), allocatable, intent(out) :: p,rho
        real(8), dimension(:,:), allocatable, intent(out) :: nnth
        real(8), dimension(:), allocatable :: gdet, header, header2
        type (four_vector), dimension(:), allocatable, intent(out) :: u,b
        real(8), dimension(:,:), allocatable :: grid, data, tmetric
        integer :: i, nelem

          !AC change this for new fields!!
          if(nrelbin>0) then
             dlen=41+nrelbin;
          else
             dlen = 38
          endif
           
          nhead=10; nhead2=3
          allocate(header(nhead))
          allocate(header2(nhead2))
          
          rhopos=10; ppos=33; vpos=12; bpos=25; gdetpos=16; gridpos=4; nthpos=43
          !rhopos=7; ppos=30; vpos=9; bpos=22; gdetpos=13; gridpos=3; nthpos=40

          write(6,*) 'data file: ',data_file
          open(unit=8,file=data_file,form='formatted',status='old',action='read')
          
          allocate(data(dlen,nx2)); nelem=nx2
          allocate(grid(nx1*nx2,4)); !AC 4 (x1,x2,r,theta)
          allocate(p(nx1*nx2)); allocate(rho(nx1*nx2))
          allocate(u(nx1*nx2)); allocate(b(nx1*nx2)); allocate(gdet(nx1*nx2))
          allocate(nnth(nx1*nx2, nrelbin));

          write(6,*) 'read koral sizes', dlen, nx1, nx2, size(data)
          read(8,*) header
          tcur=header(1)
          
          if(doingkoralnth) then
             read(8,*) header2
          end if 

          deallocate(header)
          deallocate(header2)
          
          write(6,*) 'read koral past header: ', nx1, nx2, nelem, dlen
       
          !AC read the data
          do i=0, (nx1*nx2)/nelem-1
             read(8,*) data
             rho(1+i*nelem:(i+1)*nelem)=data(rhopos,:)
             p(1+i*nelem:(i+1)*nelem)=data(ppos,:) !AC p is actually Te
             
             b(1+i*nelem:(i+1)*nelem)%data(1)=(data(bpos,:))
             b(1+i*nelem:(i+1)*nelem)%data(2)=(data(bpos+1,:))
             b(1+i*nelem:(i+1)*nelem)%data(3)=(data(bpos+2,:))
             b(1+i*nelem:(i+1)*nelem)%data(4)=(data(bpos+3,:))
             u(1+i*nelem:(i+1)*nelem)%data(1)=(data(vpos,:))
             u(1+i*nelem:(i+1)*nelem)%data(2)=(data(vpos+1,:))
             u(1+i*nelem:(i+1)*nelem)%data(3)=(data(vpos+2,:))
             u(1+i*nelem:(i+1)*nelem)%data(4)=(data(vpos+3,:))

             !AC grid is x1,x2,r,theta
             grid(1+i*nelem:(i+1)*nelem,:)=transpose(data(gridpos:gridpos+3,:)) 
             gdet(1+i*nelem:(i+1)*nelem)=(data(gdetpos,:))

             if(doingkoralnth) then
                nnth(1+i*nelem:(i+1)*nelem, : ) = transpose(data(nthpos:nthpos+nrelbin-1, : ))
             endif

          end do
          write(6,*) 'after read loop'
          
          close(unit=8)
          deallocate(data)
          write(6,*) 'read koral grid sizes', size(x1_arr), size(x2_arr), size(r_arr), size(th_arr)
          x1_arr=grid(:,1); x2_arr=grid(:,2); r_arr=grid(:,3); th_arr=grid(:,4)
          write(6,*) 'read koral assign grid'
          deallocate(grid)

          ! Transform velocities, magnetic fields from MKS to KS and then BL:
          if(doingkoralnth) then
             write(6,*) 'read koral nnth 1 ', minval(nnth(:,1)), maxval(nnth(:,1))
          endif
          write(6,*) 'read koral Te ', minval(p), maxval(p)
          write(6,*) 'read koral gdet ', minval(gdet), maxval(gdet)
          write(6,*) 'read koral rho ', minval(rho), maxval(rho)
          write(6,*) 'read koral b0', minval(b%data(1)), maxval(b%data(1))
          write(6,*) 'read koral u0 ', minval(u%data(1)), maxval(u%data(1))
!          write(6,*) 'read koral u4 ', minval(u%data(4)), maxval(u%data(4))

          write(6,*) 'read koral transform coords r ', minval(r_arr), maxval(r_arr), asim
          write(6,*) 'read koral transform coords th ', minval(th_arr), maxval(th_arr), asim
          write(6,*) 'WARNING: ignoring minpolecell=',minpolecell, ' cells on polar axis!!'
!          write(6,*) 'read koral transform coords  ', minval(x2_arr), maxval(x2_arr), h
!          write(6,*) 'read koral transform coords  ', minval(x1_arr), maxval(x1_arr), h

          
        end subroutine read_koral_data_file

        subroutine initialize_koral_model(a,ifile,doingnth,df,hf,ntt,indft,sfac,nrelbin_bk,relgammamin_bk,relgammamax_bk)
        real(kind=8), intent(in) :: a
        integer,intent(out),optional :: nrelbin_bk
        real,intent(out),optional :: relgammamin_bk, relgammamax_bk
        integer :: nx, status, nhead
        logical, intent(in) :: doingnth
        real(8),intent(in),optional :: sfac
        character(len=20), intent(in), optional :: ifile
        character(len=20) :: default_ifile='koral.in' !AC
        character(len=100), intent(in), optional :: df,hf
        integer, intent(in), optional :: ntt,indft
        nhead=10

        doingkoralnth = doingnth !ANDREW AA AC will there be issues here?
        
! if all inputs are given, assign variables rather than reading from file
        if (present(df)) then
           dfile = df
           hfile = hf
           nt = ntt
           indf = indft
        else
           if (present(ifile)) then
              call read_koral_inputs(ifile)
           else
              call read_koral_inputs(default_ifile)
           endif
        endif

        if (present(sfac)) then
            scalefac=sfac
        else
            scalefac=1.
        endif

        call read_koral_data_header(nhead)

        nrelbin_bk = nrelbin
        relgammamin_bk = relgammamin
        relgammamax_bk = relgammamax

        if (abs(asim-a).gt.1e-4) then 
           write(6,*) 'ERROR -- Different simulation and grtrans spin values!'
           return
        endif
        n=nx1*nx2
        call init_koral_data(n,n*nt)
        call load_koral_data(nt)
        write(6,*) 'nrelbin',nrelbin
        end subroutine initialize_koral_model

        subroutine load_koral_data(nt)
        real(8), dimension(:), allocatable :: rho,p
        real(8), dimension(:,:), allocatable :: nnth  
        real, dimension(:), allocatable :: vrl, vtl, vpl
        type (four_vector), dimension(:), allocatable :: u, b
        integer, intent(in) :: nt
        character(len=250) :: append
        character(len=100) :: data_file
        integer :: k
        real(8) :: tcur, tstep_test

        allocate(nnth(n,nrelbin))
        allocate(rho(n)); allocate(p(n))
        allocate(vrl(n)); allocate(vtl(n)); allocate(vpl(n))
        allocate(u(n)); allocate(b(n))
        do k=1,nt
           write(append, fmt='(I3.3)') indf-(k-1)
           data_file = trim(dfile) // append
!           write(6,*) 'data_file: ',indf-(k-1),append,data_file
           call read_koral_data_file(data_file,tcur,rho,p,u,b,nnth)
           t(k)=tcur
           ! transform velocities to LNRF frame
!           write(6,*) 'lnrf sizes: ',size(u%data(2)), size(r_arr), size(th_arr), size(vrl), size(vtl)
           call lnrf_frame(real(u%data(2)/u%data(1)),real(u%data(3)/u%data(1)), & 
                           real(u%data(4)/u%data(1)), &
                           r_arr,asim,th_arr,vrl,vtl,vpl)
           write(6,*) 'lnrf transform', size(b0_arr), size(b%data(1)), size(b0_arr((k-1)*n+1:k*n))
           write(6,*) 'lnrf transform', size(vrl), size(vtl), size(vpl), size(u0_arr), size(nth_arr(:,1))

           ! now assign to data arrays
           write(6,*) scalefac,sqrt(scalefac)
           rho_arr((k-1)*n+1:k*n)=rho*scalefac

!AC p_arr is actually the ELECTRON TEMPERATURE straight from koral -- no scalefac!!
           p_arr((k-1)*n+1:k*n)=p!*scalefac
           b0_arr((k-1)*n+1:k*n)=b%data(1)*sqrt(scalefac)
           br_arr((k-1)*n+1:k*n)=b%data(2)*sqrt(scalefac)
           bth_arr((k-1)*n+1:k*n)=b%data(3)*sqrt(scalefac)
           bph_arr((k-1)*n+1:k*n)=b%data(4)*sqrt(scalefac)
           u0_arr((k-1)*n+1:k*n)=u%data(1)!*scalefac !??TODO why would u^0 be scaled?
           vrl_arr((k-1)*n+1:k*n)=real(vrl)
           vpl_arr((k-1)*n+1:k*n)=real(vpl)
           vtl_arr((k-1)*n+1:k*n)=real(vtl)

           if(doingkoralnth) then
              nth_arr((k-1)*n+1:k*n, :) = nnth*scalefac;
           endif
        end do
        !write(6,*) nt
! Need to be able to do light curves even when nload (nt) = 1... Maybe read headers for first two files.
        if (nt.gt.1) then
! use default sim time step unless told otherwise
           tstep=t(1)-t(2)
           tstep_test=t(nt-1)-t(nt)
           if (abs((tstep-tstep_test)/tstep).gt.0.1) then 
              write(6,*) 'WARNING -- Time step changes over dumps!'
           endif
        endif

        deallocate(rho)
        deallocate(p)
        deallocate(vrl); deallocate(vtl); deallocate(vpl)
        deallocate(u); deallocate(b)
        deallocate(nnth)
        end subroutine  load_koral_data

        subroutine advance_koral_timestep(dtobs)
        real(8), intent(in) :: dtobs
        integer :: nupdate
        nupdate = floor(dtobs / tstep)
        toffset = toffset+dtobs-tstep*nupdate
! after a couple stpes might be accumulating offset
        nupdate=nupdate+floor(toffset/tstep)
        toffset=toffset-floor(toffset/tstep)*tstep
        indf=indf+nupdate
        write(6,*) 'advance koral timestep: ',dtobs,tstep,toffset,indf,nupdate !AC
        call update_koral_data(nupdate)
        end subroutine advance_koral_timestep

        subroutine update_koral_data(nupdate)
        integer, intent(in) :: nupdate
        integer :: nshift
        if(nupdate.gt.nt) then
! this is case where we just load all new data
           call load_koral_data(nt)
        else
! this is case where we shift existing data back and then load more
           nshift=n*nupdate+1
           rho_arr(nshift:n*nt)=rho_arr(1:n*nt-nshift)
           p_arr(nshift:n*nt)=p_arr(1:n*nt-nshift)
           br_arr(nshift:n*nt)=br_arr(1:n*nt-nshift)
           bth_arr(nshift:n*nt)=bth_arr(1:n*nt-nshift)
           bph_arr(nshift:n*nt)=bph_arr(1:n*nt-nshift)
           b0_arr(nshift:n*nt)=b0_arr(1:n*nt-nshift)
           vrl_arr(nshift:n*nt)=vrl_arr(1:n*nt-nshift)
           vtl_arr(nshift:n*nt)=vtl_arr(1:n*nt-nshift)
           vpl_arr(nshift:n*nt)=vpl_arr(1:n*nt-nshift)
           if(doingkoralnth) then
              nth_arr(nshift:n*nt,:)=nth_arr(1:n*nt-nshift,:)
           endif
           call load_koral_data(nupdate)
        end if
        end subroutine update_koral_data

        subroutine init_koral_data(nx,n)
          integer, intent(in) :: n,nx
          integer :: i
          real :: logbinspace, binspace
!        integer, intent(in), optional :: aux,td,auxn
        allocate(rho_arr(n)); allocate(p_arr(n)); allocate(u0_arr(n))
        allocate(vrl_arr(n)); allocate(vtl_arr(n)); allocate(vpl_arr(n))
        allocate(b0_arr(n)); allocate(br_arr(n)); allocate(bth_arr(n))
        allocate(bph_arr(n))
        allocate(x1_arr(nx)); allocate(x2_arr(nx)); allocate(r_arr(nx))
        allocate(th_arr(nx)); allocate(t(nt))

        !AC nonthermal
        !AC probably don't need the relel_gammas here - move to emis
        if(doingkoralnth) then
           allocate(nth_arr(n,nrelbin))
        end if
        
        end subroutine init_koral_data

        subroutine del_koral_data()
        deallocate(rho_arr); deallocate(p_arr); deallocate(u0_arr)
        deallocate(vrl_arr); deallocate(vtl_arr); deallocate(vpl_arr)
        deallocate(b0_arr); deallocate(br_arr); deallocate(bth_arr)
        deallocate(bph_arr); deallocate(x1_arr); deallocate(x2_arr)
        deallocate(r_arr); deallocate(th_arr); deallocate(t)
        if(doingkoralnth) then
           deallocate(nth_arr)
        end if
        end subroutine del_koral_data

      end module fluid_model_koral
