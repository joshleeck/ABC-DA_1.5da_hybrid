SUBROUTINE Deallocate_um_data (state)

USE DefConsTypes, ONLY :   &
    UM_type

IMPLICIT NONE

TYPE(UM_type), INTENT(INOUT)   :: state

DEALLOCATE(state % longs_u, state % longs_v, state % half_levs, state % full_levs)
DEALLOCATE(state % u, state % v, state % w)
DEALLOCATE(state % density, state % theta, state % exner_pressure, state % orog_height)

END SUBROUTINE Deallocate_um_data
