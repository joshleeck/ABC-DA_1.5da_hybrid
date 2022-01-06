SUBROUTINE Deallocate_model_vars(state)

USE DefConsTypes, ONLY :     &
    ABC_type

IMPLICIT NONE

TYPE(ABC_type), INTENT(INOUT)   :: state

DEALLOCATE(state % u, state % v, state % w, state % r, state % b, state % tracer)
DEALLOCATE(state % rho, state % b_ef, state % hydro_imbal, state % geost_imbal)
DEALLOCATE(state % vert_mom_source, state % horiz_div, state % horiz_vort)

END SUBROUTINE Deallocate_model_vars
