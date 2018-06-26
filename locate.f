      SUBROUTINE dlocate(xx,x,indx)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: xx
      DOUBLE PRECISION, INTENT(IN) :: x
      INTEGER, INTENT(OUT) :: indx
      INTEGER :: n,jl,jm,ju
      LOGICAL :: ascnd
!      write(6,*) 'xx: ',xx,x,size(xx)
      n=size(xx)
      ascnd = (xx(n) >= xx(1))
      jl=0
      ju=n+1
      do
        if (ju-jl <= 1) exit
        jm=(ju+jl)/2
        if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
        else
          ju=jm
        end if
      end do
      if (x == xx(1)) then
        indx=1
      else if (x == xx(n)) then
        indx=n-1
      else
        indx=jl
      end if
      END SUBROUTINE dlocate

      SUBROUTINE locate(xx,x,indx)
      IMPLICIT NONE
      REAL, DIMENSION(:), INTENT(IN) :: xx
      REAL, INTENT(IN) :: x
      INTEGER, INTENT(OUT) :: indx
      INTEGER :: n,jl,jm,ju
      LOGICAL :: ascnd
!      write(6,*) 'xx: ',xx,x,size(xx)
      n=size(xx)
      ascnd = (xx(n) >= xx(1))
      jl=0
      ju=n+1
      do
        if (ju-jl <= 1) exit
        jm=(ju+jl)/2
        if (ascnd .eqv. (x >= xx(jm))) then
          jl=jm
        else
          ju=jm
        end if
      end do
      if (x == xx(1)) then
        indx=1
      else if (x == xx(n)) then
        indx=n-1
      else
        indx=jl
      end if
      END SUBROUTINE locate

