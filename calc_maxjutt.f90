       module calc_maxjutt
!       use math
!        use polsynchemis, only: initialize_polsynchpl,polsynchpl, &
!        polsynchth,del_polsynchpl,synchpl,bnu,synchemis,synchemisnoabs
        use polsynchemis, only: polsynchth

       contains
!         subroutine calc_maxjutt_subroutine(nsim,usim,maxjutt_args,maxjutt_thetas,mj_ncgs,mj_tcgs)
         subroutine calc_maxjutt_subroutine(nsim,usim,nu,bcgs,incang,maxjutt_args,ktotal)
           double precision, dimension(:), intent(in) :: nsim,usim
           double precision, dimension(:), intent(in) :: nu,bcgs,incang
           double precision, dimension(:), intent(in) :: maxjutt_args
!6 weights, 1 delta

!           double precision, intent(in) :: nu,bcgs,incang
!           double precision, dimension(size(nsim)) :: ones
           double precision, dimension(size(maxjutt_args) - 1) :: weights_arr, i_arr, delta_arr

           double precision, dimension(size(nsim)) :: theta_min
           double precision :: delta, a, dwsum

           double precision, dimension(size(nsim),11) :: ktemp
           double precision, dimension(size(nsim),11),intent(out) :: ktotal
!           ones = 1d0
           a = 3d0
           delta = maxjutt_args(1) !first argument is delta
           weights_arr = maxjutt_args(2:)/SUM(maxjutt_args(2:)) !normalize
           i_arr = (/ (I, I=0,size(weights_arr)-1) /) !indices
           delta_arr = delta**i_arr !all deltas to be multiplied by thetamin
           dwsum = SUM(weights_arr * delta_arr)
           theta_min = (usim/nsim)/a/dwsum 
!           write(6,*) weights_arr
!           write(6,*) delta_arr
!N array of theta_mins calculated to satisfy sum energy = usim
!3/5/2015 Checking thetamin calculation. Weight calculation is straightforward. Comparing to edf.py
!Comparison revealed a = 3 discrepancy, now in agreement
!weights_arr and delta_arr working as intended

!pseudocode
!ktotal = whatever
!convert units to cgs?
!for i in (1,size(weights_arr)):
!    ktemp  = polsynchth(nu,weights_arr(i)*ncgs,bcgs,thetamin*delta_arr(i),e%incang,Kth)           
!    ktotal = ktotal + ktemp
!return K to emis.f90

           do i=1,size(weights_arr)
              ktemp = 0d0
              call polsynchth(nu,weights_arr(i)*nsim,bcgs,theta_min*delta_arr(i),incang,ktemp)
              ktotal = ktotal + ktemp
!ktotal = ktotal + ktemp
           enddo
         end subroutine calc_maxjutt_subroutine
         end module calc_maxjutt
