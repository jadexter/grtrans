      module class_four_Vector

!     filename: class_Vector.inc 
!     public, everything by default, but can specify any 
      type four_Vector
!        private
        real(kind=8), dimension(4) :: data ! component values
        real(kind=8), dimension(10) :: metric ! metric
      end type four_Vector
 
      interface size
        module procedure four_Vector_size
      end interface

      interface lower
         module procedure lower
         module procedure lower_arr
      end interface

      interface assign_metric
        module procedure assign_metric, assign_metric_arr
        module procedure assign_metric_real, assign_metric_arr_real
        module procedure assign_metric_arr_real_same
        module procedure assign_metric_arr_same
      end interface

!     Overload common operators 
      interface operator (+)    ! add others later 
      module procedure add_four_Vector, add_Double_to_four_Vector
      module procedure add_four_Vector_arr
      module procedure add_Double_to_four_Vector_arr
      module procedure add_real_to_four_Vector_arr
      module procedure add_real_to_four_Vector
      end interface 

      interface operator (-)    ! add unary versions later 
      module procedure subtract_four_Vector, subtract_Double
      module procedure subtract_four_Vector_arr, subtract_Double_arr
      module procedure subtract_Real,subtract_Real_arr
      module procedure subtract_four_Vector_arr_single
      module procedure subtract_four_Vector_single_arr
      end interface 

      interface operator (*)    ! overload * 
      module procedure dot_four_Vector, double_mult_four_Vector
      module procedure four_Vector_mult_double 
      module procedure dot_four_Vector_arr, double_mult_four_Vector_arr
      module procedure dot_four_vector_arr_single
      module procedure dot_four_Vector_single_arr
      module procedure four_Vector_mult_double_arr
      module procedure four_Vector_mult_real_arr
      module procedure real_mult_four_Vector_arr
      module procedure real_mult_four_Vector,four_Vector_mult_real
      end interface 

      interface assignment (=)  ! overload = 
      module procedure equal_Double, equal_Double_arr
      module procedure equal_Real, equal_Real_arr
      module procedure assign, assign_arr
      end interface 

      interface operator (==)   ! overload == 
      module procedure is_equal_to, is_equal_to_arr
      end interface 

      contains                     ! functions & operators 

        pure function four_Vector_size(name) result (n) ! accessor member
        type (four_Vector), intent(in), dimension(:) :: name
        integer :: n
        n = size(name%data(1)) ; end function four_Vector_size 

        function add_Double_to_four_Vector (v, r) result (new) ! overload + 
        type (four_Vector), intent(in) :: v 
        real(kind=8),          intent(in) :: r 
        type (four_Vector)             :: new ! new = v + r 
        new%data = v%data + r ; end function 
    
        function add_four_Vector (a, b) result (new) ! vector + vector 
        type (four_Vector), intent(in) :: a, b 
        type (four_Vector)             :: new ! new = a + b 
!        if ( a%size /= b%size ) stop "Sizes differ in add_four_Vector" 
!        allocate ( new%data(a%size) ) ; new%size = a%size 
        new%data = a%data + b%data    ; end function add_four_Vector 

        subroutine assign (name, values)
        real(kind=8), intent(in) :: values(4)
        type (four_Vector), intent(inout) :: name
        name%data=values
        end subroutine assign
     
!        function assign (values) result (name) ! array to vector constructor
!        real(kind=8), intent(in) :: values(4) ! given rank 1 array 
!        integer          :: length ! array size 
!        type (four_Vector)    :: name  ! four_Vector to create 
!        length = size(values); allocate ( name%data(length) ) 
!        name % data = values; end function assign 

        function copy_four_Vector (name) result (new) 
        type (four_Vector), intent(in) :: name 
        type (four_Vector)             :: new 
        new%data = name%data; new%metric = name%metric
        end function copy_four_Vector 
    
        subroutine  delete_four_Vector (name) ! deallocate allocated items 
        type (four_Vector), intent(inout) :: name 
!        integer                      :: ok ! check deallocate status 
!        deallocate (name%data, stat = ok ) 
!        if ( ok /= 0 ) stop "four_Vector &
!          not allocated in delete_four_Vector" 
!        name%size = 0 ;
        name % data = 0; name % metric = 0
        end subroutine delete_four_Vector 
 
        function dot_four_Vector (a, b) result (c) ! overload * 
        type (four_Vector), intent(in) :: a, b
        type (four_Vector) :: d
        real(kind=8)                      :: c 
!        if ( a%size /= b%size ) stop "Sizes differ in dot_four_Vector"
        d=lower(a) ! change a to co-variant
        c = dot_product (d%data, b%data) ; end function dot_four_Vector 

        function lower(a) result(b)
        type (four_Vector), intent(in) :: a
        type (four_Vector) :: b
        b%data=(/a%metric(1)*a%data(1)+a%metric(2)*a%data(2)+a%metric(3) &
        *a%data(3)+ &
        a%metric(4)*a%data(4),a%metric(2)*a%data(1)+a%metric(5)* &
        a%data(2)+a%metric(6)* &
        a%data(3)+a%metric(7)*a%data(4),a%metric(3)*a%data(1)+ &
        a%metric(6)*a%data(2)+ &
        a%metric(8)*a%data(3)+a%metric(9)*a%data(4),a%metric(4)* &
        a%data(1)+a%metric(7)* &
        a%data(2)+a%metric(9)*a%data(3)+a%metric(10)*a%data(4)/) 
        end function lower

        function lower_arr(a) result(b)
        type (four_Vector), dimension(:), intent(in) :: a
        type (four_Vector), dimension(size(a)) :: b
        integer :: i
        do i=1,size(a)
           b(i)=lower(a(i))
        end do
        end function lower_arr

        subroutine equal_Double (new, R) ! overload =, real(kind=8) to vector 
        type (four_Vector), intent(inout) :: new 
        real(kind=8),          intent(in)    :: R 
!        if ( associated (new%data) ) deallocate (new%data) 
!        allocate ( new%data(1) ); new%size = 1 
        new%data = R            ; end subroutine equal_Double 

        logical function is_equal_to (a, b) result (t_f) ! overload == 
        type (four_Vector), intent(in) :: a, b ! left & right of == 
        t_f = .false.             ! initialize 
!        if ( a%size /= b%size ) return ! same size ? 
        t_f = all ( a%data == b%data ) ! and all values match 
        end function is_equal_to 

!        function length (name) result (n) ! accessor member 
!        type (four_Vector), intent(in) :: name 
!        integer                  :: n 
!        n = name % size ; end function length 

        subroutine list (name)    ! accessor member 
        type (four_Vector), intent(in) :: name 
        print *,"[", name % data, "]"; end subroutine list 
     
        function make_four_Vector (values, m_values) result(v) ! Optional Constructor 
!        integer, optional, intent(in) :: len ! number of values 
        real(kind=8), optional, intent(in) :: values(4)    ! given values
        real(kind=8), optional, intent(in) :: m_values(10) ! metric values 
        type (four_Vector)                 :: v 
!        if ( present (len) ) then ! create vector data 
!          v%size = len ; allocate ( v%data(len) ) 
          if ( present (values)) then ; v%data = values ! vector 
          else                      ; v%data = 0.d0 ! null vector 
          end if                 ! values present
          if (present(m_values)) then ; v%metric = m_values
          else ; v% metric=(/-1.,0.,0.,0.,1.,0.,0.,1.,0.,1./)
          end if
        end function make_four_Vector 

        function normalize_four_Vector (name)  result (new) 
        type (four_Vector), intent(in) :: name 
        type (four_Vector)             :: new 
        real(kind=8)               :: total, nil = epsilon(nil) ! tolerance
        total = sqrt ( name*name ) ! use raise/lower op
        if ( total < nil ) then ; new%data = 0.d0 ! avoid division by 0 
        else                  ; new%data = name%data/total 
        end if                  ; end function normalize_four_Vector 
        
        function double_mult_four_Vector (r, v) result (new) ! overload * 
        real(kind=8),          intent(in) :: r 
        type (four_Vector), intent(in) :: v 
        type (four_Vector)             :: new ! new = r * v 
        new%data = r * v%data         
        end function double_mult_four_Vector 
      
        function subtract_Double (v, r) result (new) ! vector-real(kind=8), overload - 
        type (four_Vector), intent(in) :: v 
        real(kind=8),          intent(in) :: r 
        type (four_Vector)             :: new ! new = v + r 
        new%data = v%data - r         ; end function subtract_double 

        function subtract_four_Vector (a, b) result (new) ! overload - 
        type (four_Vector), intent(in) :: a, b 
        type (four_Vector)             :: new 
        new%data = a%data - b%data ; end function subtract_four_Vector 

        function values (name) result (array) ! accessor member

        type (four_Vector), intent(in) :: name 
        real(kind=8)                      :: array(4) 
        array = name % data ; end function values 

        function four_Vector_ (values, metric) result(name) ! Public constructor 
!        integer,      intent(in) :: length ! array size 
        real(kind=8), intent(in) :: values(4) ! given
        real(kind=8), optional, intent(in) :: metric(10)
!        real(kind=8), pointer            :: pt_to_val(:) ! pointer to array 
        type (four_Vector)            :: name ! four_Vector to create 
        if(.not.(present(metric))) then 
          name=four_Vector(values,(/-1.,0.,0.,0.,1.,0.,0.,1.,0.,1./))
        else
        name      = four_Vector(values,metric) ! intrinsic constructor 
        end if
        end function four_Vector_

        function four_Vector_max_value (a) result (v) ! accessor member 
        type (four_Vector), intent(in) :: a 
        real(kind=8)                      :: v 
        v = maxval ( a%data ); end function four_Vector_max_value 

        function four_Vector_min_value (a) result (v) ! accessor member 
        type (four_Vector), intent(in) :: a 
        real(kind=8)                      :: v 
        v = minval ( a%data ) ; end function four_Vector_min_value 

        function four_Vector_mult_double (v, r) result (new) ! vector*real(kind=8), overload * 
        type (four_Vector), intent(in) :: v 
        real(kind=8),          intent(in) :: r 
        type (four_Vector)             :: new            ! new = v * r 
!        if ( v%size < 1 ) stop "Zero size in four_Vector_mult_real(kind=8)" 
        new = Double_mult_four_Vector (r, v)
        end function four_Vector_mult_double 

        function add_Double_to_four_Vector_arr (v, r) result (new) ! overload + 
        type (four_Vector), intent(in), dimension(:) :: v 
        real(kind=8),          intent(in)  :: r 
        type (four_Vector), dimension(size(v)) :: new ! new = v + r 
        new%data(1) = v%data(1) + r
        new%data(2) = v%data(2) + r
        new%data(3) = v%data(3) + r
        new%data(4) = v%data(4) + r
        end function 
    
        function add_four_Vector_arr (a, b) result (new) ! vector + vector 
        type (four_Vector), intent(in), dimension(:) :: a, b 
        type (four_Vector), dimension(size(a)) :: new ! new = a + b 
        integer :: i
        do i=1,size(a); new(i)%data=a(i)%data+b(i)%data; enddo
        end function add_four_Vector_arr 
     
        subroutine assign_arr (name,values)
! array to vector constructor
        real(kind=8), intent(in) :: values(:,:) ! given rank 1 array 
        integer          :: i ! array size 
        type (four_Vector), dimension(size(values,2)), intent(inout) &
           :: name  
! four_Vector to create 
        do i=1,size(name); name(i)%data=values(:,i); enddo
        end subroutine assign_arr

        subroutine assign_metric(name,values)
        real(kind=8), intent(in), dimension(10) :: values
        type (four_Vector) :: name
        name%metric=values
        end subroutine assign_metric

        subroutine assign_metric_real(name,values)
        real, intent(in), dimension(10) :: values
        type (four_Vector) :: name
        name%metric=values
        end subroutine assign_metric_real

        subroutine assign_metric_arr_real(name,values)
        real, intent(in), dimension(:,:) :: values
        type (four_Vector), dimension(:) :: name
        integer :: i
        do i=1,size(values,2); name(i)%metric=values(:,i); enddo
        end subroutine assign_metric_arr_real

        subroutine assign_metric_arr_real_same(name,values)
        real, intent(in), dimension(10) :: values
        type (four_Vector), dimension(:) :: name
        integer :: i
        do i=1,size(name); name(i)%metric=values; enddo
        end subroutine assign_metric_arr_real_same

        subroutine assign_metric_arr_same(name,values)
        real(kind=8), intent(in), dimension(10) :: values
        type (four_Vector), dimension(:) :: name
        integer :: i
        do i=1,size(name); name(i)%metric=values; enddo
        end subroutine assign_metric_arr_same

        subroutine assign_metric_arr(name,values)
        real(kind=8), intent(in), dimension(:,:) :: values
        type (four_Vector), dimension(:) :: name
        integer :: i
        do i=1,size(values,2); name(i)%metric=values(:,i); enddo
        end subroutine assign_metric_arr

        function copy_four_Vector_arr (name) result (new) 
        type (four_Vector), intent(in), dimension(:) :: name 
        type (four_Vector), dimension(size(name)) :: new 
        integer :: i
        do i=1,size(name)
          new(i)%data = name(i)%data; new(i)%metric=name(i)%metric
        enddo
        end function copy_four_Vector_arr 
    
        subroutine  delete_four_Vector_arr (name) ! deallocate allocated items
        type (four_Vector), intent(inout), dimension(:) :: name 
        integer                      :: i ! check deallocate status 
        do i=1,size(name); name(i)%data=0;name(i)%metric=0; enddo
        end subroutine delete_four_Vector_arr 
 
        function dot_four_Vector_arr (a, b) result (c) ! overload * 
        type (four_Vector), intent(in), dimension(:) :: a, b
        type (four_Vector) :: d
        real(kind=8),dimension(size(a)) :: c
        integer :: i
        do i=1,size(a)
          d=lower(a(i)) ! change a to co-variant
          c(i)=dot_product(d%data,b(i)%data)
        enddo
        end function dot_four_Vector_arr

        function dot_four_Vector_arr_single (a, b) result (c) ! overload * 
        type (four_Vector), intent(in), dimension(:) :: a
        type (four_Vector), intent(in) :: b
        type (four_Vector) :: d
        real(kind=8),dimension(size(a)) :: c
        integer :: i
        do i=1,size(a)
          d=lower(a(i)) ! change a to co-variant
          c(i)=dot_product(d%data,b%data)
        enddo
        end function dot_four_Vector_arr_single

        function dot_four_Vector_single_arr (a, b) result (c) ! overload * 
        type (four_Vector), intent(in), dimension(:) :: b
        type (four_vector), intent(in) :: a
        type (four_Vector) :: d
        real(kind=8),dimension(size(b)) :: c
        integer :: i
        d=lower(a)
        do i=1,size(b)
          !d=lower(a(i)) ! change a to co-variant
          c(i)=dot_product(d%data,b(i)%data)
        enddo
        end function dot_four_Vector_single_arr

!        function lower_arr(a) result(b)
!        type (four_Vector), intent(in) :: a
!        type (four_Vector) :: b
!!        b%data=(/a%metric(1)*a%data(1)+a%metric(2)*a%data(2)+a%metric(3) & &
!        *a%data(3)+ &
!        a%metric(4)*a%data(4),a%metric(2)*a%data(1)+a%metric(5)* &
!        a%data(2)+a%metric(6)* &
!        a%data(3)+a%metric(7)*a%data(4),a%metric(2)*a%data(1)+ &
!!        a%metric(6)*a%data(2)+ &
!        a%metric(8)*a%data(3)+a%metric(9)*a%data(4),a%metric(4)* &
!        a%data(1)+a%metric(7)* &
!        a%data(2)+a%metric(9)*a%data(3)+a%metric(10)*a%data(4)/)
!        end function lower_arr

        subroutine equal_Double_arr (new, R) ! overload =, real(kind=8) to vector 
        type (four_Vector), intent(inout), dimension(:) :: new 
        real(kind=8), intent(in) :: R 
        do i=1,size(new);new(i)%data=R;enddo 
        end subroutine equal_Double_arr 

        logical function is_equal_to_arr (a, b) result (t_f) ! overload == 
        type (four_Vector), intent(in), dimension(:) :: a, b ! left & right of == 
        integer :: i
        t_f = .false. ! initialize
        do i=1,size(a)
          t_f = all( a(i)%data == b(i)%data ) ! and all values match
          if(.not.t_f) exit
        enddo
        end function is_equal_to_arr

        subroutine list_arr (name)    ! accessor member 
        type (four_Vector), intent(in) :: name 
        print *,"[", name % data, "]"; end subroutine list_arr 
     
        function make_four_Vector_arr (values, m_values) result(v) ! Optional Constructor 
!        integer, optional, intent(in) :: len ! number of values 
        real(kind=8), intent(in) :: values(:,:)    
! given values
        real(kind=8), intent(in) :: m_values(:,:) 
! metric values 
        type (four_Vector), dimension(size(values,2)) :: v 
        integer :: i
!        if ( present (len) ) then ! create vector data 
!          v%size = len ; allocate ( v%data(len) )
        do i=1,size(values,2)
          v(i)%data = values(:,i) ! vector 
          v(i)%metric = m_values(:,i)
        enddo
        end function make_four_Vector_arr 

        function normalize_four_Vector_arr (name)  result (new) 
        type (four_Vector), intent(in), dimension(:) :: name 
        type (four_Vector), dimension(size(name))    :: new 
        real(kind=8)          :: total, nil=epsilon(nil) 
        integer :: i
! tolerance
        do i=1,size(name)
        total = sqrt ( name(i)*name(i) ) ! use raise/lower op
        if ( total < nil ) then ; new(i)%data = 0.d0 ! avoid division by 0 
        else ; new(i)%data = name(i)%data/total
        end if; enddo; end function normalize_four_Vector_arr  
        
        function double_mult_four_Vector_arr (r, v) result (new) ! overload * 
        real(kind=8), intent(in) :: r 
        type (four_Vector), intent(in), dimension(:) :: v 
        type (four_Vector), dimension(size(v)) :: new ! new = r * v 
        integer :: i
        do i=1,size(v); new(i)%data = r*v(i)%data; enddo         
        end function double_mult_four_Vector_arr 

        function size_four_Vector (name) result (n) ! accessor member 
        type (four_Vector), intent(in), dimension(:) :: name 
        integer                   :: n 
        n = size(name) ; end function size_four_Vector 
      
        function subtract_Double_arr (v, r) result (new) ! vector-real(kind=8), overload - 
        type (four_Vector), intent(in), dimension(:) :: v 
        real(kind=8), intent(in) :: r 
        type (four_Vector), dimension(size(v)) :: new ! new = v + r 
        integer :: i
        do i=1,size(v); new(i)%data = v(i)%data - r; enddo
        end function subtract_double_arr 

        function subtract_four_Vector_arr (a, b) result (new) ! overload -
        type (four_Vector), intent(in), dimension(:) :: a, b 
        type (four_Vector), dimension(size(a)) :: new 
        integer :: i
        do i=1,size(a); new(i)%data = a(i)%data - b(i)%data; enddo
        end function subtract_four_Vector_arr 
        
        function subtract_four_Vector_single_arr (a, b) result (new) ! overload -
        type (four_Vector), intent(in), dimension(:) :: b
        type (four_Vector), intent(in) :: a
        type (four_Vector), dimension(size(b)) :: new 
        integer :: i
        do i=1,size(b); new(i)%data = a%data - b(i)%data; enddo
        end function subtract_four_Vector_single_arr

        function subtract_four_Vector_arr_single (a, b) result (new) ! overload -
        type (four_Vector), intent(in), dimension(:) :: a
        type (four_Vector), intent(in) :: b
        type (four_Vector), dimension(size(a)) :: new
        integer :: i
        do i=1,size(a); new(i)%data = a(i)%data - b%data; enddo
        end function subtract_four_Vector_arr_single

        function values_arr (name) result (array) ! accessor member
        type (four_Vector), intent(in), dimension(:) :: name 
        real(kind=8) :: array(4,size(name)) 
        integer :: i
        do i=1,size(name); array(:,i)=name(i)%data ; enddo
        end function values_arr 

        function four_Vector_mult_double_arr (v, r) result (new) ! vector*real(kind=8), overload * 
        type (four_Vector), intent(in), dimension(:) :: v 
        real(kind=8), intent(in) :: r 
        type (four_Vector), dimension(size(v)) :: new  ! new = v * r 
        new = Double_mult_four_Vector_arr (r, v)
        end function four_Vector_mult_double_arr


        function add_real_to_four_Vector (v, r) result (new) ! overload + 
        type (four_Vector), intent(in) :: v 
        real,          intent(in) :: r 
        type (four_Vector)             :: new ! new = v + r 
!        if ( v%size < 1 ) stop "No sizes in add_real Precision_to_four_Vector" 
!        allocate ( new%data(4) )
!  new%data = v%data + r ! as array operation, or use implied loop 
        new%data = v%data + r ; end function 
     
        subroutine assign_real (name,values) ! array to vector constructor
        real, intent(in) :: values(4) ! given rank 1 array 
        type (four_Vector)    :: name  ! four_Vector to create 
        name % data = values; end subroutine assign_real
 
        function dot_four_Vector_real (a, b) result (c) ! overload * 
        type (four_Vector), intent(in) :: a, b
        type (four_Vector) :: d
        real :: c
!        if ( a%size /= b%size ) stop "Sizes differ in dot_four_Vector"
        d=lower(a) ! change a to co-variant
        c = dot_product (d%data, b%data)
        end function dot_four_Vector_real 

        subroutine equal_real (new, R) ! overload =, real to vector 
        type (four_Vector), intent(inout) :: new 
        real,          intent(in)    :: R 
!        if ( associated (new%data) ) deallocate (new%data) 
!        allocate ( new%data(1) ); new%size = 1 
        new%data = R            ; end subroutine equal_real 
     
        function make_four_Vector_real (values, m_values) result(v) ! Optional Constructor 
!        integer, optional, intent(in) :: len ! number of values 
        real, optional, intent(in) :: values(4)    ! given values
        real, optional, intent(in) :: m_values(10) ! metric values 
        type (four_Vector)                 :: v 
!        if ( present (len) ) then ! create vector data 
!          v%size = len ; allocate ( v%data(len) ) 
          if ( present (values)) then ; v%data = values ! vector 
          else                      ; v%data = 0.d0 ! null vector 
          end if                 ! values present
          if (present(m_values)) then ; v%metric = m_values
          else ; v%metric=(/-1.,0.,0.,0.,1.,0.,0.,1.,0.,1./)
          end if
        end function make_four_Vector_real

        function normalize_four_Vector_real (name)  result (new) 
        type (four_Vector), intent(in) :: name 
        type (four_Vector)             :: new 
        real               :: total, nil = epsilon(nil) ! tolerance
        total = sqrt ( name*name ) ! use raise/lower op
        if ( total < nil ) then ; new%data = 0.d0 ! avoid division by 0 
        else                  ; new%data = name%data/total 
        end if     ; end function normalize_four_Vector_real 
        
        function real_mult_four_Vector (r, v) result (new) ! overload * 
        real,          intent(in) :: r 
        type (four_Vector), intent(in) :: v 
        type (four_Vector)             :: new ! new = r * v 
        new%data = r * v%data         
        end function real_mult_four_Vector 
      
        function subtract_real (v, r) result (new) ! vector-real, overload
        type (four_Vector), intent(in) :: v 
        real,          intent(in) :: r 
        type (four_Vector)             :: new ! new = v + r 
        new%data = v%data - r         ; end function subtract_real 

        function values_real (name) result (array) ! accessor member
        type (four_Vector), intent(in) :: name 
        real                      :: array(4) 
        array = name % data ; end function values_real 

        function four_Vector_real_ (values, metric) result(name) ! Public constructor 
!        integer,      intent(in) :: length ! array size 
        real, intent(in) :: values(4) ! given
        real, optional, intent(in) :: metric(10)
!        real, pointer            :: pt_to_val(:) ! pointer to array 
        type (four_Vector)            :: name ! four_Vector to create 
        if(.not.(present(metric))) then 
          name=four_Vector(values,(/-1.,0.,0.,0.,1.,0.,0.,1.,0.,1./))
        else
        name      = four_Vector(values,metric) ! intrinsic constructor 
        end if
        end function four_Vector_real_

        function four_Vector_max_value_real (a) result (v) ! accessor member 
        type (four_Vector), intent(in) :: a 
        real                      :: v 
        v = maxval ( a%data ); end function four_Vector_max_value_real 

        function four_Vector_min_value_real (a) result (v) ! accessor member 
        type (four_Vector), intent(in) :: a 
        real                      :: v 
        v = minval ( a%data ) ; end function four_Vector_min_value_real 

        function four_Vector_mult_real (v, r) result (new) ! vector*real, overload * 
        type (four_Vector), intent(in) :: v 
        real,          intent(in) :: r 
        type (four_Vector)             :: new            ! new = v * r 
!        if ( v%size < 1 ) stop "Zero size in four_Vector_mult_real" 
        new = real_mult_four_Vector (r, v)
        end function four_Vector_mult_real 

        function add_real_to_four_Vector_arr (v, r) result (new) ! overload + 
        type (four_Vector), intent(in), dimension(:) :: v 
        real,          intent(in)  :: r 
        type (four_Vector), dimension(size(v)) :: new ! new = v + r 
        integer :: i
        do i=1,size(v); new(i)%data=v(i)%data+r; enddo
        end function 
     
        subroutine assign_arr_real (name,values) 
! array to vector constructor
        real, intent(in) :: values(:,:) ! given rank 2 array 
        integer          :: i ! array size 
        type (four_Vector), dimension(size(values,1)) :: name 
! four_Vector to create 
        do i=1,size(name); name(i)%data=values(:,i); enddo
        end subroutine assign_arr_real
 
        function dot_four_Vector_arr_real (a, b) result (c) ! overload * 
        type (four_Vector), intent(in), dimension(:) :: a, b
        type (four_Vector) :: d
        real,dimension(size(a)) :: c
        integer :: i
!        if ( a%size /= b%size ) stop "Sizes differ in dot_four_Vector"
        do i=1,size(a)
          d=lower(a(i)) ! change a to co-variant
          c(i)=dot_product(d%data,b(i)%data)
        enddo
        end function dot_four_Vector_arr_real

        subroutine equal_real_arr (new, R) ! overload =, real to vector 
        type (four_Vector), intent(inout), dimension(:) :: new 
        real, intent(in) :: R 
        do i=1,size(new);new(i)%data=R;enddo 
        end subroutine equal_real_arr
     
        function make_four_Vector_arr_real (values, m_values) result(v) ! Optional Constructor 
!        integer, optional, intent(in) :: len ! number of values 
        real, intent(in) :: values(:,:)    
! given values
        real, intent(in) :: m_values(:,:) 
! metric values 
        type (four_Vector), dimension(size(values,2)) :: v 
        integer :: i
!        if ( present (len) ) then ! create vector data 
!          v%size = len ; allocate ( v%data(len) )
        do i=1,size(values,2)
          v(i)%data = values(:,i) ! vector 
          v(i)%metric = m_values(:,i)
        enddo
        end function make_four_Vector_arr_real

        function normalize_four_Vector_arr_real (name)  result (new) 
        type (four_Vector), intent(in), dimension(:) :: name 
        type (four_Vector), dimension(size(name))    :: new 
        real          :: total, nil=epsilon(nil) 
        integer :: i
! tolerance
        do i=1,size(name)
        total = sqrt ( name(i)*name(i) ) ! use raise/lower op
        if ( total < nil ) then ; new(i)%data = 0.d0 ! avoid division by 0 
        else ; new(i)%data = name(i)%data/total
        end if; enddo; end function normalize_four_Vector_arr_real  
        
        function real_mult_four_Vector_arr (r, v) result (new) ! overload * 
        real, intent(in) :: r 
        type (four_Vector), intent(in), dimension(:) :: v 
        type (four_Vector), dimension(size(v)) :: new ! new = r * v 
        integer :: i
        do i=1,size(v); new(i)%data = r*v(i)%data; enddo         
        end function real_mult_four_Vector_arr 
      
        function subtract_real_arr (v, r) result (new) ! vector-real, overload - 
        type (four_Vector), intent(in), dimension(:) :: v 
        real, intent(in) :: r
        type (four_Vector), dimension(size(v)) :: new ! new = v + r 
        integer :: i
        do i=1,size(v); new(i)%data = v(i)%data - r; enddo
        end function subtract_real_arr 

        function values_arr_real (name) result (array) ! accessor member
        type (four_Vector), intent(in), dimension(:) :: name 
        real :: array(4,size(name)) 
        integer :: i
        do i=1,size(name); array(:,i)=name(i)%data ; enddo
        end function values_arr_real

        function four_Vector_mult_real_arr (v, r) result (new) ! vector*real, overload * 
        type (four_Vector), intent(in), dimension(:) :: v 
        real, intent(in) :: r 
        type (four_Vector), dimension(size(v)) :: new  ! new = v * r 
        new = real_mult_four_Vector_arr (r, v)
        end function four_Vector_mult_real_arr

      end module class_four_Vector
