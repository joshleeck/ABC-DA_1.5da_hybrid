SUBROUTINE Add_pert_model_vars(state_out, state1_in, state2_in, weight)

! Computes the following
!   state_out = state1_in + weight * state2_in

USE DefConsTypes, ONLY :     &
    ZREAL8,                  &
    ABC_type,                 &
    nlongs, nlevs

IMPLICIT NONE

TYPE(ABC_type),  INTENT(INOUT)   :: state_out
TYPE(ABC_type),  INTENT(IN)      :: state1_in
TYPE(ABC_type),  INTENT(IN)      :: state2_in
REAL(ZREAL8),   INTENT(IN)      :: weight

state_out % u(0:nlongs+1,0:nlevs+1) = state1_in % u(0:nlongs+1,0:nlevs+1) + &
                                       weight * state2_in % u(0:nlongs+1,0:nlevs+1)

state_out % v(0:nlongs+1,0:nlevs+1) = state1_in % v(0:nlongs+1,0:nlevs+1) + &
                                       weight * state2_in % v(0:nlongs+1,0:nlevs+1)

state_out % w(0:nlongs+1,0:nlevs+1) = state1_in % w(0:nlongs+1,0:nlevs+1) + &
                                       weight * state2_in % w(0:nlongs+1,0:nlevs+1)

state_out % r(0:nlongs+1,0:nlevs+1) = state1_in % r(0:nlongs+1,0:nlevs+1) + &
                                       weight * state2_in % r(0:nlongs+1,0:nlevs+1)

state_out % b(0:nlongs+1,0:nlevs+1) = state1_in % b(0:nlongs+1,0:nlevs+1) + &
                                       weight * state2_in % b(0:nlongs+1,0:nlevs+1)


END SUBROUTINE Add_pert_model_vars
