       module calc_maxcomp

        use polsynchemis, only: polsynchth

       contains

!         subroutine calc_maxjutt_subroutine(nu,ncgs,bcgs,tcgs,usim,incang,maxjutt_args,ktotal)
         subroutine calc_maxcomp_subroutine(nu,ncgs,bcgs,tcgs,incang,maxjutt_args,ktotal)
           double precision, dimension(:), intent(in) :: nu
           double precision, dimension(:), intent(in) :: ncgs,bcgs,tcgs,incang
           double precision, dimension(:), intent(in) :: maxjutt_args
!6 weights, 1 delta, 1 selection option
!           double precision, dimension(:), intent(in) :: usim
!2,3,5,7,11,13,17
           double precision, dimension(size(maxjutt_args) - 2) :: weights_arr, i_arr, delta_arr
           double precision, dimension(size(ncgs)) :: tcgs_min
!           double precision, dimension(size(ncgs)) :: theta_min
           double precision :: delta, a, dwsum

           double precision, dimension(size(ncgs),11) :: ktemp
           double precision, dimension(size(ncgs),11),intent(out) :: ktotal

           integer :: isel

!           ones = 1d0
!           write(6,*) 'maxjutt: ',size(maxjutt_args),size(weights_arr)
           a = 3d0
!Note derivation assumes relativistic a = 3.0
           delta = maxjutt_args(1) !first argument is delta
           selection = maxjutt_args(2) !second argument selects the decomposition
           weights_arr = maxjutt_args(3:)/SUM(maxjutt_args(3:)) !normalize
!           write(6,*) 'maxjutt weights: ',weights_arr,delta
           i_arr = (/ (I, I=0,size(weights_arr)-1) /) !indices
           delta_arr = delta**i_arr !all deltas to be multiplied by thetamin
           dwsum = SUM(weights_arr * delta_arr)
           tcgs_min = tcgs/dwsum

!pseudocode
!ktotal = whatever
!convert units to cgs?
!for i in (1,size(weights_arr)):
!    ktemp  = polsynchth(nu,weights_arr(i)*ncgs,bcgs,thetamin*delta_arr(i),e%incang,Kth)
!    ktotal = ktotal + ktemp
!total absorption
!partial absorption from selection
!selection = 6 returns the top maxwellian

!return K to emis.f90
           ktotal = 0d0
           do i=1,size(weights_arr)
              ktemp = 0d0
!              write(6,*) 'calc_maxjutt call polsynchth: ',weights_arr(i),ncgs(1),bcgs(1),tcgs_min(1)*delta_arr(1),incang(1)
              call polsynchth(nu,weights_arr(i)*ncgs,bcgs,tcgs_min*delta_arr(i),incang,ktemp)
              ktotal = ktotal + ktemp
           enddo


           isel = INT(selection)
           if(isel.gt.0) then
              if(isel.le.(size(weights_arr))) then
                 ktemp = 0d0
                 call polsynchth(nu,weights_arr(isel)*ncgs,bcgs,tcgs_min*delta_arr(isel),incang,ktemp)
                 ktotal(:,1:4) = ktemp(:,1:4)
              endif
           endif

!           write(6,*) 'maxjutt loop: ',tcgs_min,delta_arr,i,size(weights_arr)
         end subroutine calc_maxcomp_subroutine

         end module calc_maxcomp
