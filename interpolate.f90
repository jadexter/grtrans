       module interpolate
         
       interface interp_1d
         subroutine polyvl(NDER,XX,YFIT,YP,N,X,C,WORK,IERR)
         integer, intent(in) :: NDER,N
         integer, intent(out) :: IERR
         real(kind=8), intent(in), dimension(N) :: X,C
         real(kind=8), intent(in) :: XX
         real(kind=8), intent(out) :: YFIT
         real(kind=8), intent(in), dimension(:) :: YP
         real(kind=8), intent(in), dimension(:) :: WORK
         end subroutine polyvl
       end interface

       interface interp
         module procedure bilininterp
         module procedure trilininterp
         module procedure quadlininterp
         module procedure dbilininterp
         module procedure dtrilininterp
         module procedure dquadlininterp
       end interface

       interface polint
         subroutine polint(N,X,Y,C)
         integer, intent(in) :: N
         real(kind=8), intent(in), dimension(N) :: X,Y
         real(kind=8), intent(out), dimension(N) :: C
         end subroutine polint
       end interface

       interface get_weight
         module procedure dget_weight
         module procedure get_weight
         module procedure get_weight_indx
         module procedure get_weight_arr
       end interface

       interface locate
         subroutine dlocate(xx,x,indx)
         real(kind=8), dimension(:), intent(in) :: xx
         real(kind=8), intent(in) :: x
         integer, intent(out) :: indx
         end subroutine dlocate
         subroutine locate(xx,x,indx)
         real, dimension(:), intent(in) :: xx
         real, intent(in) :: x
         integer, intent(out) :: indx
         end subroutine locate
       end interface
 
       interface hunt
         subroutine dhunt(xx,x,jlo)
         integer, intent(inout) :: jlo
         real(kind=8), intent(in) :: x
         real(kind=8), dimension(:), intent(in) :: xx
         end subroutine dhunt
         subroutine hunt(xx,x,jlo)
         integer, intent(inout) :: jlo
         real, intent(in) :: x
         real, dimension(:), intent(in) :: xx
         end subroutine hunt
       end interface

       contains

         subroutine dget_weight(xx,x,jlo,weight)
         real(kind=8), intent(in) :: x
         real(kind=8), dimension(:), intent(in) :: xx
         integer, intent(inout) :: jlo
         real(kind=8), intent(out) :: weight
         call hunt(xx,x,jlo)
         weight=(x-xx(jlo))/(xx(jlo+1)-xx(jlo))
         end subroutine dget_weight

         subroutine get_weight(xx,x,jlo,weight)
         real, intent(in) :: x
         real, dimension(:), intent(in) :: xx
         integer, intent(inout) :: jlo
         real, intent(out) :: weight
         call hunt(xx,x,jlo)
         weight=(x-xx(jlo))/(xx(jlo+1)-xx(jlo))
         end subroutine get_weight

         subroutine get_weight_indx(xx,x,weight)
         real, intent(in) :: x
         real, dimension(:), intent(in) :: xx
         integer :: jlo
         real, intent(out) :: weight
         call locate(xx,x,jlo)
         weight=(x-xx(jlo))/(xx(jlo+1)-xx(jlo))+jlo
         end subroutine get_weight_indx

         subroutine get_weight_arr(xx,x,indx,weight)
         real, dimension(:), intent(in) :: x,xx
         real, dimension(size(x)), intent(out) :: weight
         integer :: i,jlo
         integer, dimension(size(x)), intent(out) :: indx
         call locate(xx,x(1),jlo)
         indx(1)=jlo
         do i=2,size(x)
           call hunt(xx,x(i),jlo)
           indx(i)=jlo
         enddo
         weight=(x-xx(indx))/(xx(indx+1)-xx(indx))
         end subroutine get_weight_arr

        function dbilininterp(vv,xd,yd) result(y)
        ! Trilinearly interpolate given points (xd,yd,zd) assumed to be normalized so grid spacing is 1 in each direction.
        ! Based on formula from wikipedia.
        ! JAD 9/8/2008
        real(kind=8), intent(in), dimension(:) :: xd,yd
        real(kind=8), intent(in), dimension(size(xd),4) :: vv
        real(kind=8), dimension(size(xd)) :: y,w1,w2
        ! Interpolate along y
        w1=vv(:,1)*(1.-yd)+vv(:,2)*yd
        w2=vv(:,3)*(1.-yd)+vv(:,4)*yd
        ! Finally, interpolate y weights along x
        y=w1*(1.-xd)+w2*xd
        end function dbilininterp

        function dtrilininterp(vv,xd,yd,zd) result(y)
        ! Trilinearly interpolate given points (xd,yd,zd) assumed to be normalized so grid spacing is 1 in each direction.
        ! Based on formula from wikipedia.
        ! JAD 9/8/2008
        real(kind=8), intent(in), dimension(:) :: xd,yd,zd
        real(kind=8), intent(in), dimension(size(xd),8) :: vv
        real(kind=8), dimension(size(xd)) :: y,i1,i2,j1,j2,w1,w2
!        write(6,*) 'dtrilin'
        ! First interpolate along z:
        i1=vv(:,1)*(1.-zd)+vv(:,2)*zd
        i2=vv(:,3)*(1.-zd)+vv(:,4)*zd
        j1=vv(:,5)*(1.-zd)+vv(:,6)*zd
        j2=vv(:,7)*(1.-zd)+vv(:,8)*zd
        ! Then interpolate these z-weights along y
        w1=i1*(1.-yd)+i2*yd
        w2=j1*(1.-yd)+j2*yd
        ! Finally, interpolate y,z weights along x
        y=w1*(1.-xd)+w2*xd
        end function dtrilininterp

        function dquadlininterp(vv,td,xd,yd,zd) result(y)
        ! Quadrilinearly interpolate given points (td,xd,yd,zd) assumed to be normalized so grid spacing is 1 in each direction.
        ! Based on trilinear interpolation formula from wikipedia.
        ! JAD 9/8/2008
        real(kind=8), intent(in), dimension(:) :: td,xd,yd,zd
        real(kind=8), intent(in), dimension(size(td),16) :: vv
        real(kind=8), dimension(size(td)) :: y,i1,i2,j1,j2,k1, &
         k2,l1,l2,m1,m2,n1,n2,w1,w2
        ! First interpolate along z:
        i1=vv(:,1)*(1.-zd)+vv(:,2)*zd
        i2=vv(:,3)*(1.-zd)+vv(:,4)*zd
        j1=vv(:,5)*(1.-zd)+vv(:,6)*zd
        j2=vv(:,7)*(1.-zd)+vv(:,8)*zd
        k1=vv(:,9)*(1.-zd)+vv(:,10)*zd
        k2=vv(:,11)*(1.-zd)+vv(:,12)*zd
        l1=vv(:,13)*(1.-zd)+vv(:,14)*zd
        l2=vv(:,15)*(1.-zd)+vv(:,16)*zd
        ! Then interpolate these z-weights along y
        m1=i1*(1.-yd)+i2*yd
        m2=j1*(1.-yd)+j2*yd
        n1=k1*(1.-yd)+k2*yd
        n2=l1*(1.-yd)+l2*yd
        ! Then interpolate y,z weights along x:
        w1=m1*(1.-xd)+m2*xd
        w2=n1*(1.-xd)+n2*xd
        ! Finally, interpolate x,y,z weights along t:
        y=w1*(1.-td)+w2*td
        end function dquadlininterp

        function bilininterp(vv,xd,yd) result(y)
        ! Trilinearly interpolate given points (xd,yd,zd) assumed to be normalized so grid spacing is 1 in each direction.
        ! Based on formula from wikipedia.
        ! JAD 9/8/2008
        real, intent(in), dimension(:) :: xd,yd
        real, intent(in), dimension(:,:) :: vv
        real, dimension(size(xd)) :: y,w1,w2
        ! Interpolate along y
        w1=vv(:,1)*(1.-yd)+vv(:,2)*yd
        w2=vv(:,3)*(1.-yd)+vv(:,4)*yd
        ! Finally, interpolate y weights along x
        y=w1*(1.-xd)+w2*xd
        end function bilininterp

        function trilininterp(vv,xd,yd,zd) result(y)
        ! Trilinearly interpolate given points (xd,yd,zd) assumed to be normalized so grid spacing is 1 in each direction.
        ! Based on formula from wikipedia.
        ! JAD 9/8/2008
        real, intent(in), dimension(:) :: xd,yd,zd
        real, intent(in), dimension(:,:) :: vv
        real, dimension(size(xd)) :: y,i1,i2,j1,j2,w1,w2
!        write(6,*) 'trilin'
        ! First interpolate along z:
        i1=vv(:,1)*(1.-zd)+vv(:,2)*zd
        i2=vv(:,3)*(1.-zd)+vv(:,4)*zd
        j1=vv(:,5)*(1.-zd)+vv(:,6)*zd
        j2=vv(:,7)*(1.-zd)+vv(:,8)*zd
        ! Then interpolate these z-weights along y
        w1=i1*(1.-yd)+i2*yd
        w2=j1*(1.-yd)+j2*yd
        ! Finally, interpolate y,z weights along x
        y=w1*(1.-xd)+w2*xd
        end function trilininterp

        function quadlininterp(vv,td,xd,yd,zd) result(y)
        ! Quadrilinearly interpolate given points (td,xd,yd,zd) assumed to be normalized so grid spacing is 1 in each direction.
        ! Based on trilinear interpolation formula from wikipedia.
        ! JAD 9/8/2008
        real, intent(in), dimension(:) :: td,xd,yd,zd
        real, intent(in), dimension(:,:) :: vv
        real, dimension(size(td)) :: y,i1,i2,j1,j2,k1, &
         k2,l1,l2,m1,m2,n1,n2,w1,w2
        ! First interpolate along z:
        i1=vv(:,1)*(1.-zd)+vv(:,2)*zd
        i2=vv(:,3)*(1.-zd)+vv(:,4)*zd
        j1=vv(:,5)*(1.-zd)+vv(:,6)*zd
        j2=vv(:,7)*(1.-zd)+vv(:,8)*zd
        k1=vv(:,9)*(1.-zd)+vv(:,10)*zd
        k2=vv(:,11)*(1.-zd)+vv(:,12)*zd
        l1=vv(:,13)*(1.-zd)+vv(:,14)*zd
        l2=vv(:,15)*(1.-zd)+vv(:,16)*zd
        ! Then interpolate these z-weights along y
        m1=i1*(1.-yd)+i2*yd
        m2=j1*(1.-yd)+j2*yd
        n1=k1*(1.-yd)+k2*yd
        n2=l1*(1.-yd)+l2*yd
        ! Then interpolate y,z weights along x:
        w1=m1*(1.-xd)+m2*xd
        w2=n1*(1.-xd)+n2*xd
        ! Finally, interpolate x,y,z weights along t:
        y=w1*(1.-td)+w2*td
        end function quadlininterp

      end module interpolate
