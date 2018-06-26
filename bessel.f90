       module bessel
         implicit none
         
       interface besselk0
          module procedure besselk0
       end interface

       interface besselk1
          module procedure besselk1
       end interface

       interface besselk
          module procedure besselk
       end interface

       contains

         function besselk(n,x)
           integer, intent(in) :: n
           real(8), dimension(:), intent(in) :: x
           real(8), dimension(size(x)) :: besselk
           integer :: j
           real(8), dimension(size(x)) :: bk,bkm,bkp,tox
           tox=2d0/x
           bkm=besselk0(x)
           bk=besselk1(x)
!           if(n.eq.0) then
!           besselk=bkm
!           else if(n.eq.1) then
!              besselk=bk
!           else
           do j=1,n-1
              bkp=bkm+j*tox*bk
              bkm=bk
              bk=bkp
           end do
           besselk=bk
!           end if
         end function besselk

         function besselk0(x)
           real(8), dimension(:), intent(in) :: x
           real(8), dimension(size(x)) :: besselk0
           real(8), dimension(size(x)) :: y
           logical, dimension(size(x)) :: mask
           real(8), dimension(7) :: p = (/-0.57721566d0,0.42278420d0,&
                0.23069756d0,0.3488590d-1,0.262698d-2,0.10750d-3,&
                0.74d-5/)
           real(8), dimension(7) :: q = (/1.25331414d0,-0.7832358d-1,&
                0.2189568d-1,-0.1062446d-1,0.587872d-2,&
                -0.251540d-2,0.53208d-3/)
!           write(6,*) 'bessel k0 before mask'
           mask = (x <= 2.0d0)
!           write(6,*) 'before where mask'
           where (mask)
              y=x*x/4.0d0
              besselk0=(-log(x/2.0d0)*besseli0(x))+poly_msk(y,p,mask)
           elsewhere
              y=(2.0d0/x)
              besselk0=(exp(-x)/sqrt(x))*poly_msk(y,q,.not.mask)
           end where
         end function besselk0
         
         function besselk1(x)
           real(8), dimension(:), intent(in) :: x
           real(8), dimension(size(x)) :: besselk1
           real(8), dimension(size(x)) :: y
           logical, dimension(size(x)) :: mask
           real(8), dimension(7) :: p = (/1.0d0,0.15443144d0,&
                -0.67278579d0,-0.18156897d0,-0.1919402d-1,&
                -0.110404d-2,-0.4686d-4/)
           real(8), dimension(7) :: q = (/1.25331414d0,0.23498619d0,&
                -0.3655620d-1,0.1504268d-1,-0.780353d-2,&
                0.325614d-2,-0.68245d-3/)
!           write(6,*) 'bessel k1 before mask'
           mask = (x <= 2.0d0)
!           write(6,*) 'bessel k1 after mask'
           where (mask)
              y=x*x/4.0d0
              besselk1=(log(x/2.0d0)*besseli1(x))+(1.0d0/x)*poly_msk(y,p,mask)
           elsewhere
              y=2.0d0/x
              besselk1=(exp(-x)/sqrt(x))*poly_msk(y,q,.not.mask)
           end where
         end function besselk1
         
         function besseli0(x)
           real(8), dimension(:), intent(in) :: x
           real(8), dimension(size(x)) :: besseli0
           real(8), dimension(size(x)) :: ax
           real(8), dimension(size(x)) :: y
           logical, dimension(size(x)) :: mask
           real(8), dimension(7) :: p = (/1.00d0,3.51562290d0,&
                3.08994240d0,1.20674920d0,0.26597320d0,0.360768d-10,&
                0.45813d-20/)
           real(8), dimension(9) :: q = (/0.398942280d0,0.1328592d-10,&
                0.225319d-20,-0.157565d-20,0.916281d-20,&
                -0.2057706d-10,0.2635537d-10,-0.1647633d-10,&
                0.392377d-20/)
           ax=abs(x)
           mask = (ax < 3.75d0)
           where (mask)
              besseli0=poly_msk((x/3.750d0)**2,p,mask)
           elsewhere
            y=3.750d0/ax
            besseli0=(exp(ax)/sqrt(ax))*poly_msk(y,q,.not.mask)
         end where
       end function besseli0
       
       function besseli1(x)
         real(8), dimension(:), intent(in) :: x
         real(8), dimension(size(x)) :: besseli1
         real(8), dimension(size(x)) :: ax
         real(8), dimension(size(x)) :: y
         logical, dimension(size(x)) :: mask
         real(8), dimension(7) :: p = (/0.50d0,0.878905940d0,&
              0.514988690d0,0.150849340d0,0.2658733d-10,&
              0.301532d-20,0.32411d-30/)
          real(8), dimension(9) :: q = (/0.398942280d0,-0.3988024d-10,&
               -0.362018d-20,0.163801d-20,-0.1031555d-10,&
               0.2282967d-10,-0.2895312d-10,0.1787654d-10,&
               -0.420059d-20/)
          ax=abs(x)
          mask = (ax < 3.75d0)
          where (mask)
             besseli1=ax*poly_msk((x/3.750d0)**2,p,mask)
          elsewhere
             y=3.750d0/ax
             besseli1=(exp(ax)/sqrt(ax))*poly_msk(y,q,.not.mask)
          end where
          where (x < 0.0d0) besseli1=-besseli1
        end function besseli1
        
        function poly_msk(x,coeffs,mask)
          real(8), dimension(:), intent(in) :: coeffs,x
          logical, dimension(:), intent(in) :: mask
          real(8), dimension(size(x)) :: poly_msk
          poly_msk=unpack(poly(pack(x,mask),coeffs),mask,0.0d0)
        end function poly_msk
        
        function poly(x,coeffs)
          real(8), dimension(:), intent(in) :: coeffs,x
          real(8), dimension(size(x)) :: poly
          integer :: i,m
          m=size(coeffs)
!          n=size(x)
          if (m <= 0) then
             poly=0.0d0
          else !if (m < n .or. m < NPAR_POLY) then
             poly=coeffs(m)
             do i=m-1,1,-1
                poly=x*poly+coeffs(i)
             end do
          endif!else
!             do i=1,n
!                poly(i)=poly_rr(x(i),coeffs)
!             end do
!          end if
        end function poly

        end module bessel
