      module class_Vector

!     filename: class_Vector.inc 
!     public, everything by default, but can specify any 
      type Vector 
        integer                     :: size ! vector length 
        real, pointer, dimension(:) :: data ! component values 
      end type Vector
 
!     Overload common operators 
      interface operator (+)    ! add others later 
      module procedure add_Vector, add_Real_to_Vector ; end interface 

      interface operator (-)    ! add unary versions later 
      module procedure subtract_Vector, subtract_Real ; end interface 

      interface operator (*)    ! overload * 
      module procedure dot_Vector, real_mult_Vector, Vector_mult_real 
      end interface 

      interface assignment (=)  ! overload = 
      module procedure equal_Real ; end interface 

      interface operator (==)   ! overload == 
      module procedure is_equal_to ; end interface 

      contains                                  ! functions & operators 

        function add_Real_to_Vector (v, r) result (new) ! overload + 
        type (Vector), intent(in) :: v 
        real,          intent(in) :: r 
        type (Vector)             :: new ! new = v + r 
        if ( v%size < 1 ) stop "No sizes in add_Real_to_Vector" 
        allocate ( new%data(v%size) ) ; new%size = v%size 
!  new%data = v%data + r ! as array operation, or use implied loop 
        new%data(1:v%size) = v%data(1:v%size) + r ; end function 
    
        function add_Vector (a, b) result (new) ! vector + vector 
        type (Vector), intent(in) :: a, b 
        type (Vector)             :: new ! new = a + b 
        if ( a%size /= b%size ) stop "Sizes differ in add_Vector" 
        allocate ( new%data(a%size) ) ; new%size = a%size 
        new%data = a%data + b%data    ; end function add_Vector 
     
        function assign (values) result (name) ! array to vector constructor
        real, intent(in) :: values(:) ! given rank 1 array 
        integer          :: length ! array size 
        type (Vector)    :: name  ! Vector to create 
        length = size(values); allocate ( name%data(length) ) 
        name % size = length ; name % data = values; end function assign 
        function copy_Vector (name) result (new) 
        type (Vector), intent(in) :: name 
        type (Vector)             :: new 
        allocate ( new%data(name%size) ) ; new%size = name%size 
        new%data = name%data             ; end function copy_Vector 
    
        subroutine  delete_Vector (name) ! deallocate allocated items 
        type (Vector), intent(inout) :: name 
        integer                      :: ok ! check deallocate status 
        deallocate (name%data, stat = ok ) 
        if ( ok /= 0 ) stop "Vector not allocated in delete_Vector" 
        name%size = 0 ; end subroutine delete_Vector 

        function dot_Vector (a, b) result (c) ! overload * 
        type (Vector), intent(in) :: a, b 
        real                      :: c 
        if ( a%size /= b%size ) stop "Sizes differ in dot_Vector" 
        c = dot_product (a%data, b%data) ; end function dot_Vector 

        subroutine equal_Real (new, R) ! overload =, real to vector 
        type (Vector), intent(inout) :: new 
        real,          intent(in)    :: R 
        if ( associated (new%data) ) deallocate (new%data) 
        allocate ( new%data(1) ); new%size = 1 
        new%data = R            ; end subroutine equal_Real 

        logical function is_equal_to (a, b) result (t_f) ! overload == 
        type (Vector), intent(in) :: a, b ! left & right of == 
        t_f = .false.             ! initialize 
        if ( a%size /= b%size ) return ! same size ? 
        t_f = all ( a%data == b%data ) ! and all values match 
        end function is_equal_to 

        function length (name) result (n) ! accessor member 
        type (Vector), intent(in) :: name 
        integer                  :: n 
        n = name % size ; end function length 

        subroutine list (name)    ! accessor member 
        type (Vector), intent(in) :: name 
        print *,"[", name % data(1:name%size), "]"; end subroutine list 
     
        function make_Vector (len, values) result(v) ! Optional Constructor 
        integer, optional, intent(in) :: len ! number of values 
        real,    optional, intent(in) :: values(:) ! given values 
        type (Vector)                 :: v 
        if ( present (len) ) then ! create vector data 
          v%size = len ; allocate ( v%data(len) ) 
          if ( present (values)) then ; v%data = values ! vector 
          else                      ; v%data = 0.d0 ! null vector 
          end if                 ! values present
        else                      ! scalar constant 
          v%size = 1                  ; allocate ( v%data(1) ) ! default 
          if ( present (values)) then ; v%data(1) = values(1) ! scalar 
          else                      ; v%data(1) = 0.d0 ! null 
          end if                 ! value present 
        end if                    ! len present 
        end function make_Vector 

        function normalize_Vector (name)  result (new) 
        type (Vector), intent(in) :: name 
        type (Vector)             :: new 
        real                      :: total, nil = epsilon(nil) ! tolerance 
        allocate ( new%data(name%size) ) ; new%size = name%size 
        total = sqrt ( sum ( name%data**2 ) ) ! intrinsic functions 
        if ( total < nil ) then ; new%data = 0.d0 ! avoid division by 0 
        else                  ; new%data = name%data/total 
        end if                  ; end function normalize_Vector 

        subroutine read_Vector (name) ! read array, assign 
        type (Vector), intent(inout) :: name 
        integer, parameter           :: max = 999 
        integer                      :: length 
        read (*,'(i1)', advance = 'no') length 
        if ( length <= 0 )   stop "Invalid length in read_Vector" 
        if ( length >= max ) stop "Maximum length in read_Vector" 
        allocate ( name % data(length) ) ; name % size = length 
        read *, name % data(1:length)    ; end subroutine read_Vector 
        
        function real_mult_Vector (r, v) result (new) ! overload * 
        real,          intent(in) :: r 
        type (Vector), intent(in) :: v 
        type (Vector)             :: new ! new = r * v 
        if ( v%size < 1 ) stop "Zero size in real_mult_Vector" 
        allocate ( new%data(v%size) ) ; new%size = v%size 
        new%data = r * v%data         ; end function real_mult_Vector 

        function size_Vector (name) result (n) ! accessor member 
        type (Vector), intent(in) :: name 
        integer                   :: n 
        n = name % size ; end function size_Vector 
      
        function subtract_Real (v, r) result (new) ! vector-real, overload - 
        type (Vector), intent(in) :: v 
        real,          intent(in) :: r 
        type (Vector)             :: new ! new = v + r 
        if ( v%size < 1 ) stop "Zero length in subtract_Real" 
        allocate ( new%data(v%size) ) ; new%size = v%size 
        new%data = v%data - r         ; end function subtract_Real 

        function subtract_Vector (a, b) result (new) ! overload - 
        type (Vector), intent(in) :: a, b 
        type (Vector)             :: new 
        if ( a%size /= b%size ) stop "Sizes differ in subtract_Vector" 
        allocate ( new%data(a%size) ) ; new%size = a%size 
        new%data = a%data - b%data    ; end function subtract_Vector 

        function values (name) result (array) ! accessor member

        type (Vector), intent(in) :: name 
        real                      :: array(name%size) 
        array = name % data ; end function values 

        function Vector_ (length, values) result(name) ! Public constructor 
        integer,      intent(in) :: length ! array size 
        real, target, intent(in) :: values(length) ! given array 
        real, pointer            :: pt_to_val(:) ! pointer to array 
        type (Vector)            :: name ! Vector to create 
        integer                  :: get_m ! allocate flag 
        allocate ( pt_to_val (length), stat = get_m ) ! allocate 
        if ( get_m /= 0 ) stop 'allocate error' ! check 
        pt_to_val = values        ! dereference values 
        name      = Vector(length, pt_to_val) ! intrinsic constructor 
        end function Vector_ 

        function Vector_max_value (a) result (v) ! accessor member 
        type (Vector), intent(in) :: a 
        real                      :: v 
        v = maxval ( a%data(1:a%size) ); end function Vector_max_value 

        function Vector_min_value (a) result (v) ! accessor member 
        type (Vector), intent(in) :: a 
        real                      :: v 
        v = minval ( a%data(1:a%size) ) ; end function Vector_min_value 

        function Vector_mult_real (v, r) result (new) ! vector*real, overload * 
        type (Vector), intent(in) :: v 
        real,          intent(in) :: r 
        type (Vector)             :: new            ! new = v * r 
        if ( v%size < 1 ) stop "Zero size in Vector_mult_real" 
        new = Real_mult_Vector (r, v) ; end function Vector_mult_real 

      end module class_Vector 
