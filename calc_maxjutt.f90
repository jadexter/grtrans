       module calc_maxjutt

        use polsynchemis, only: polsynchth

       contains

!         subroutine calc_maxjutt_subroutine(nu,ncgs,bcgs,tcgs,usim,incang,maxjutt_args,ktotal)
         subroutine calc_maxjutt_subroutine(nu,ncgs,bcgs,tcgs,incang,maxjutt_args,ktotal)
           double precision, dimension(:), intent(in) :: nu
           double precision, dimension(:), intent(in) :: ncgs,bcgs,tcgs,incang
           double precision, dimension(:), intent(in) :: maxjutt_args
!6 weights, 1 delta
!           double precision, dimension(:), intent(in) :: usim

           double precision, dimension(size(maxjutt_args) - 1) :: weights_arr, i_arr, delta_arr
           double precision, dimension(size(ncgs)) :: tcgs_min
!           double precision, dimension(size(ncgs)) :: theta_min
           double precision :: delta, a, dwsum

           double precision, dimension(size(ncgs),11) :: ktemp
           double precision, dimension(size(ncgs),11),intent(out) :: ktotal
!           ones = 1d0
!           write(6,*) 'maxjutt: ',size(maxjutt_args),size(weights_arr)
           a = 3d0
!Note derivation assumes relativistic a = 3.0
           delta = maxjutt_args(1) !first argument is delta
           weights_arr = maxjutt_args(2:)/SUM(maxjutt_args(2:)) !normalize
!           write(6,*) 'maxjutt weights: ',weights_arr,delta
           i_arr = (/ (I, I=0,size(weights_arr)-1) /) !indices
           delta_arr = delta**i_arr !all deltas to be multiplied by thetamin
           dwsum = SUM(weights_arr * delta_arr)

!           theta_min = (usim/ncgs)/a/dwsum !USIM

!            u_e = 3 ncgs k tcgs
!            theta_min = 3 k tcgs / a / dwsum / m / c / c
!            k tcgs_min = theta_min * m * c * c = 
!            = 3 k tcgs / a / dwsum
!            = k tcgs / dwsum
           tcgs_min = tcgs/dwsum
!check if that second a is necessary

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
           ktotal = 0d0
           do i=1,size(weights_arr)
              ktemp = 0d0
!              write(6,*) 'calc_maxjutt call polsynchth: ',weights_arr(i),ncgs(1),bcgs(1),tcgs_min(1)*delta_arr(1),incang(1)
              call polsynchth(nu,weights_arr(i)*ncgs,bcgs,tcgs_min*delta_arr(i),incang,ktemp)
              ktotal = ktotal + ktemp
           enddo
!           write(6,*) 'maxjutt loop: ',tcgs_min,delta_arr,i,size(weights_arr)
         end subroutine calc_maxjutt_subroutine

         end module calc_maxjutt
