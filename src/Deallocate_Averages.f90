SUBROUTINE Deallocate_Averages(state)

USE DefConsTypes, ONLY :     &
    Averages_type

IMPLICIT NONE

TYPE(Averages_type), INTENT(INOUT)   :: state

DEALLOCATE(state % u_1, state % u_2, state % u_m, state % w_1, state % w_2, state % w_m)

END SUBROUTINE Deallocate_Averages
