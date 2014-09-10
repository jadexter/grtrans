!     Testing four_Vector Class Constructors & Operators 
      include 'class_vector.f' ! see previous figure 
      program check_vector_class 
      use class_four_Vector 
      type (four_Vector) :: x, y, z
      type (four_Vector), dimension(50) :: u,v,w
      double precision, dimension(4,50) :: t
      double precision, dimension(50) :: s
      real, dimension(30) :: a
      integer :: j
      a=4.; a=a+1.
      write(6,*) 'a: ',a(20)
      do j=1,50
        t(:,j)=(/3.,1.*j,1.*j-1.,1.*j+1./)
      enddo
      u=t
      write(6,*) 'u=t: ',u(25)%data-t(:,25),u(25)%data
      u=u*3d0
      u=u*.3333d0
      v=u
      write(6,*) v%data(1)-u%data(1)
      call assign_metric(u,(/-1.,0.,0.,0.,1.,0.,0.,1.,0.,1./))
      call assign_metric(v,(/-1.,0.,0.,0.,1.,0.,0.,1.,0.,1./))
      s=u*v
      write(6,*) 'u*v: ',s(17),u(17)*v(17)
      write (*,'("array of vectors: ")', advance='no')
 !     call list(t(4))
 !     write(6,*) 'add: '; call list(t(7)-3.)
!      write(6,*) dot_product(t(7)%data,t(4)%data)
 !     t(17)=(t(4)*t(7))
 !     call list(t(17))
!     test optional constructors: assign, and copy 
!      x = make_four_Vector ()        ! single scalar zero 
 !     write (*,'("made scalar x = ")', advance='no'); call list (x) 
!      call delete_four_Vector (x) ; y = make_four_Vector (4) ! 4 zero values 
 !     write (*,'("made null y = ")',   advance='no'); call list (y) 
!      z = make_four_Vector (4, (/11., 12., 13., 14./) ) ! 4 non-zero values 
 !     write (*,'("made full z = ")',   advance='no'); call list (z) 
 !     write (*,'("assign [ 31., 32., 33., 34. ] to x")') 
!      x = assign( (/31., 32., 33., 34./) ) ! (4) non-zeros 
 !     write (*,'("assigned  x = ")',   advance='no'); call list (x) 
!      x = four_Vector_(4, (/31., 32., 33., 34./) ) ! 4 non-zero values 
 !     write (*,'("public    x = ")',   advance='no'); call list (x) 
 !     write (*,'("copy x to y =")', advance='no')
!      y = copy_four_Vector (x) ; call list (y) ! copy 
!     test overloaded operators 
!      write (*,'("z * x gives ")', advance='no'); print *, z*x ! dot 
!      write (*,'("z + x gives ")', advance='no'); call list (z+x) ! add 
!      y = 25.6                  ! real to vector 
      write (*,'("y = 25.6 gives ")',   advance='no'); call list (y) 
!      y = z                     ! equality 
!      write (*,'("y = z gives y as ")',   advance='no'); call list (y) 
!      write (*,'("logic y == x gives ")', advance='no'); print *, y==x 
!      write (*,'("logic y == z gives ")', advance='no'); print *, y==z 
!     test destructor, accessors 
!      call delete_four_Vector (y)    ! destructor 
 !     write (*,'("deleting y gives y = ")', advance='no'); call list (y) 
 !     print *, "size of x is ", length (x) ! accessor 
 !     print *, "data in x are [", values (x), "]" ! accessor 
 !     write (*,'("2. times x is ")',  advance='no'); call list (2.0*x) 
 !     write (*,'("x times 2. is ")',  advance='no'); call list (x*2.0) 
 !     call delete_four_Vector (x); call delete_four_Vector (z) ! clean up 
      end program check_vector_class 
!     Running gives the output:       ! made scalar x = [0.] 
!     made null y = [0., 0., 0., 0.]   ! made full z = [11., 12., 13., 14.] 
!     assign [31., 32., 33., 34.] to x ! assigned x  = [31., 32., 33., 34.] 
!     public  x = [31., 32., 33., 34.] ! copy x to y = [31., 32., 33., 34.] 
!     z * x gives 1630.                ! z + x gives   [42., 44., 46., 48.] 
!     y = 25.6 gives [25.6000004]      ! y = z,   y =  [11., 12., 13., 14.] 
!     logic y == x gives F             ! logic y == z gives T 
!     deleting y gives y = []          ! size of x is 4 
!     data in x : [31., 32., 33., 34.] ! 2. times x is [62., 64., 66., 68.] 
!     x times 2. is [62., 64., 66., 68.] 
