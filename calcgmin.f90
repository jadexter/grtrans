
       module calc_gmin
!       use math

       contains

         subroutine calc_gmin_subroutine(gmin_p,gmin_theta,gmin_etas,gmin_gmins,gmin_nfracs)
           real(kind=8), intent(in) :: gmin_p
           real(kind=8) :: gmin_lin_cons, gmin_lin_coeff, gmin_lin_power
           real(kind=8) :: gmin_inv_cons, gmin_inv_coeff, gmin_inv_power
           real(kind=8) :: gmin_inv_sin_coeff, gmin_inv_sin_freq, gmin_inv_sin_delay
           real(kind=8) :: gmin_lin_sin_coeff, gmin_lin_sin_freq, gmin_lin_sin_delay
           real(kind=8), dimension(:), intent(in) :: gmin_theta, gmin_etas
           real(kind=8), dimension(size(gmin_theta)) :: one
           real(kind=8), dimension(size(gmin_theta)) :: gmin_linear_constant, gmin_inverse_constant, gmin_atheta
           real(kind=8), dimension(size(gmin_theta)) :: gmin_raw_gmins
           real(kind=8) :: acenter, anormal, azero, astwo
!          intermediate variables
           real(kind=8), dimension(size(gmin_theta)), intent(out) :: gmin_gmins, gmin_nfracs
           one = 1d0
           acenter = 0.5668090982352612
           anormal = 0.52624783
           azero = 3d0 / sqrt(2d0)
           astwo = LOG(SQRT(2d0))


           if(gmin_p==3.5) then
              gmin_lin_power =0.276589155355
              gmin_lin_coeff =-13.5593749125
              gmin_lin_cons =16.0797900684
              gmin_inv_power =6.53997654139
              gmin_inv_coeff =151.597731214
              gmin_inv_cons =0.722506578136
              gmin_lin_sin_coeff =-0.143307918052
              gmin_lin_sin_freq =-0.143307918052
              gmin_lin_sin_delay =-0.143307918052
              gmin_lin_sin_coeff =0.121815691108
              gmin_lin_sin_freq =0.121815691108
              gmin_lin_sin_delay =0.121815691108

           else
              gmin_lin_cons = 21.38307186
              gmin_lin_coeff = -16.7811712
              gmin_lin_power = 0.15128533
!           gmin_inv_cons = 0.53189121
!           gmin_inv_coeff = 0.91485025
!           gmin_inv_power = 1.14368368
!           gmin_inv_sin_coeff = .0111806594
!           gmin_inv_sin_freq = -17.0100822
!           gmin_inv_sin_delay = 4.03979965
              gmin_inv_cons = 0.74798712
              gmin_inv_coeff = 0.62609462
              gmin_inv_power = 0.81567379
              gmin_inv_sin_coeff = .00638946501
              gmin_inv_sin_freq = -16.8034428
              gmin_inv_sin_delay = 3.72208398
              gmin_lin_sin_coeff = 0.
              gmin_lin_sin_freq = 0.
              gmin_lin_sin_delay = 0.
           endif

           gmin_linear_constant = gmin_lin_cons + gmin_lin_coeff*(gmin_etas**gmin_lin_power)&
                + gmin_lin_sin_coeff*SIN(gmin_etas*gmin_lin_sin_freq + gmin_lin_sin_delay)
           gmin_inverse_constant = gmin_inv_cons + gmin_inv_coeff*(gmin_etas**gmin_inv_power)&
                + gmin_inv_sin_coeff*SIN(gmin_etas*gmin_inv_sin_freq + gmin_inv_sin_delay)


           gmin_raw_gmins = gmin_theta* gmin_linear_constant + gmin_inverse_constant  
           gmin_gmins = merge(gmin_raw_gmins,one,gmin_raw_gmins.ge.1d0)
!I decided to combine a(theta) *theta into a single function atheta
           gmin_atheta = gmin_theta * azero * EXP(astwo * TANH(anormal*LOG(gmin_theta / acenter)))
!azero is moving tanh down to the midpoint of 1.5 and 3.0
!astwo is scaling the amplitude of tanh in loglog space to sqrt(2)
!           gmin_atheta = -0.918344120 + 2.99892343*gmin_theta

           gmin_nfracs = gmin_etas*gmin_atheta*((gmin_p - 2.)/(gmin_p - 1.)) * (gmin_gmins**(gmin_p - 2.))

         end subroutine calc_gmin_subroutine
         end module calc_gmin
