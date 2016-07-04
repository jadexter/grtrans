      module phys_constants
      
      real(kind=8), parameter :: h=6.626d-27,k=1.38d-16, &
       c=2.99792458d10,e=4.8032d-10,G=6.67d-8, &
       m=9.10938188d-28,mp=1.67262158d-24,pi=acos(-1d0),c2=c*c,sigb=5.6704d-5, &
       msun=1.998d33, sigt=6.6523d-25

      contains
        
        subroutine phys_constants_dummy()
        end subroutine phys_constants_dummy

      end module
