SUBROUTINE Apply_alpha_model_vars(state1, state2)

! Use the same alpha fields for all variables, when inter-variable localisation is disabled

USE DefConsTypes, ONLY :     &
    ABC_type,                &
    nlongs, nlevs

IMPLICIT NONE

TYPE(ABC_type), INTENT(INOUT)   :: state1    ! Running total
TYPE(ABC_type), INTENT(IN)      :: state2


  state1 % u(0:nlongs+1,0:nlevs+1)           = state2 % u(0:nlongs+1,0:nlevs+1) + state2 % v(0:nlongs+1,0:nlevs+1) + &
                                               state2 % w(0:nlongs+1,0:nlevs+1) + state2 % r(0:nlongs+1,0:nlevs+1) + &
                                               state2 % b(0:nlongs+1,0:nlevs+1)
  state1 % v(0:nlongs+1,0:nlevs+1)           = state2 % u(0:nlongs+1,0:nlevs+1) + state2 % v(0:nlongs+1,0:nlevs+1) + &
                                               state2 % w(0:nlongs+1,0:nlevs+1) + state2 % r(0:nlongs+1,0:nlevs+1) + &
                                               state2 % b(0:nlongs+1,0:nlevs+1)
  state1 % w(0:nlongs+1,0:nlevs+1)           = state2 % u(0:nlongs+1,0:nlevs+1) + state2 % v(0:nlongs+1,0:nlevs+1) + &
                                               state2 % w(0:nlongs+1,0:nlevs+1) + state2 % r(0:nlongs+1,0:nlevs+1) + &
                                               state2 % b(0:nlongs+1,0:nlevs+1)
  state1 % r(0:nlongs+1,0:nlevs+1)           = state2 % u(0:nlongs+1,0:nlevs+1) + state2 % v(0:nlongs+1,0:nlevs+1) + &
                                               state2 % w(0:nlongs+1,0:nlevs+1) + state2 % r(0:nlongs+1,0:nlevs+1) + &
                                               state2 % b(0:nlongs+1,0:nlevs+1)
  state1 % b(0:nlongs+1,0:nlevs+1)           = state2 % u(0:nlongs+1,0:nlevs+1) + state2 % v(0:nlongs+1,0:nlevs+1) + &
                                               state2 % w(0:nlongs+1,0:nlevs+1) + state2 % r(0:nlongs+1,0:nlevs+1) + &
                                               state2 % b(0:nlongs+1,0:nlevs+1)


END SUBROUTINE Apply_alpha_model_vars
