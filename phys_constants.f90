      module phys_constants
      
      double precision, parameter :: h=6.626d-27,k=1.38d-16, &
       c=3d10,e=4.8d-10,G=6.67d-8, &
       m=9.11d-28,mp=1.67d-24,pi=acos(-1d0),c2=c*c,sigb=5.67d-5, &
       msun=1.998d33, sigt=6.65d-25

      contains
        
        subroutine phys_constants_dummy()
        end subroutine phys_constants_dummy

      end module
