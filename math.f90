       module math
 
       interface tsum
         module procedure tsum
       end interface

       interface sum
         module procedure cum_sum
       end interface

       interface dot_product
         module procedure dot_product_arr
       end interface

       interface zbrent
          module procedure zbrent
       end interface

       contains

       pure function dot_product_arr(a,b) result(dot)
       real(kind=8), dimension(:,:), intent(in) :: a,b
       real(kind=8), dimension(size(a,1)) :: dot
       integer :: i
       dot=0d0
       do i=1,size(a,2); dot=dot+a(:,i)*b(:,i); enddo
       end function dot_product_arr

       function tsum(x,y)
       real(kind=8), intent(in), dimension(:) :: x
       real(kind=8), intent(in), dimension(:) :: y
       integer :: n
       real(kind=8), dimension(size(x)) :: tsum,tsumint
       real(kind=8), dimension(size(x)-1) :: xdif, yavg
       n=size(x)
       xdif=x(2:n)-x(1:n-1)
       yavg=(y(1:n-1)+y(2:n))/2d0
       tsumint(1)=0d0
       tsumint(2:n)=xdif*yavg
!       write(6,*) 'y: ',y
!       write(6,*) 'xdif yavg: ',xdif(size(x)-1), yavg(size(x)-1)
       tsum=sum(tsumint,1)
       end function tsum

       function cum_sum(x,cumu)
       real(kind=8), intent(in), dimension(:) :: x
       real(kind=8), dimension(size(x)) :: cum_sum
       integer :: i
       integer, intent(in) :: cumu
       do i=1,size(x)
         cum_sum(i)=sum(x(1:i))
       enddo
       end function cum_sum

       function zbrent(func,x1,x2,args,tol)
         implicit none
         interface
            function func(x,args)
              implicit none
              real(8), intent(in) :: x
              real(8), intent(in), dimension(:) :: args
              real(8)             :: func
            end function func
         end interface
         real(8), intent(in), dimension(:) :: args
         real(8), intent(in) :: x1,x2,tol
         real(8)             :: zbrent
         integer, parameter  :: itmax=100
         real(8), parameter  :: eps=epsilon(x1)
         integer             :: iter
         real(8)             :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
         a=x1       ; b=x2
         fa=func(a,args) ; fb=func(b,args)
         if ((fa > 0.0 .AND. fb > 0.0) .OR. (fa < 0.0 .AND. fb < 0.0))then
            write(*,"(A)")'ROOT MUST BE BRACKETED FOR ZBRENT'
            write(*,*) x1,x2,fa,fb
            stop
         endif

         c=b ; fc=fb
         do iter=1,itmax
            if ((fb > 0.0 .AND. fc > 0.0) .OR. (fb < 0.0 .AND. fc < 0.0)) then
               c=a
               fc=fa
               d=b-a
               e=d
            end if
            if (abs(fc) < abs(fb)) then
               a=b
               b=c
               c=a
               fa=fb
               fb=fc
               fc=fa
            end if
            tol1=2.0d0*eps*abs(b)+0.5d0*tol
            xm=0.5d0*(c-b)
            if (abs(xm) <= tol1 .or. fb == 0.0) then
               zbrent=b
               return
            end if
            !
            if (abs(e) >= tol1 .AND. abs(fa) > abs(fb)) then
               s=fb/fa
               if (a == c) then
                  p=2.0d0*xm*s
                  q=1.0d0-s
               else
                  q=fa/fc
                  r=fb/fc
                  p=s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
                  q=(q-1.0d0)*(r-1.0d0)*(s-1.0d0)
               end if
               if (p > 0.0) q=-q
               p=abs(p)
               if (2.0d0*p  <  min(3.0d0*xm*q-abs(tol1*q),abs(e*q))) then
                  e=d
                  d=p/q
               else
                  d=xm
                  e=d
               end if
            else
               d=xm
               e=d
            end if
            a=b
            fa=fb
            b=b+merge(d,sign(tol1,xm), abs(d) > tol1 )
            fb=func(b,args)
         end do
         write(*,"(A)")'zbrent: exceeded maximum iterations'
         zbrent=b
       end function zbrent

       end module math
