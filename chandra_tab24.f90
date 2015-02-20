module chandra_tab24

use interpolate, only: get_weight

implicit none

real, dimension(:), allocatable :: ch_mu,ch_I,ch_delta

contains
  
  subroutine load_chandra_tab24()
  character(len=20) :: file='ch24_vals.txt'
  integer :: npts
  open(unit=8,file=file)
  read(8,*) npts
  allocate(ch_mu(npts)); allocate(ch_I(npts)); allocate(ch_delta(npts))
  read(8,*) ch_mu
  read(8,*) ch_I
  read(8,*) ch_delta
  close(unit=8)
!  write(6,*) 'load chandra'
  end subroutine load_chandra_tab24

  subroutine del_chandra_tab24()
!    write(6,*) 'del chandra'
  deallocate(ch_mu); deallocate(ch_I); deallocate(ch_delta)
  end subroutine del_chandra_tab24

  subroutine interp_chandra_tab24(mu,I,del)
  real, dimension(:), intent(in) :: mu
  real, dimension(size(mu)), intent(out) :: I, del
  integer, dimension(size(mu)) :: indx
  real, dimension(size(mu)) :: weight
  call get_weight(ch_mu,mu,indx,weight)
  I=(1.-weight)*ch_I(indx)+weight*ch_I(indx+1)
  del=(1.-weight)*ch_delta(indx)+weight*ch_delta(indx+1)
!  write(6,*) 'interp: ',size(mu),mu,I,del
  end subroutine interp_chandra_tab24

end module
