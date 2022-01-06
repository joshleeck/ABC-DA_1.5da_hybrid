SUBROUTINE Deallocate_CVs (state)

USE DefConsTypes, ONLY :     &
    CV_type

IMPLICIT NONE

TYPE(CV_type), INTENT(INOUT)   :: state

DEALLOCATE(state % v1, state % v2, state % v3, state % v4, state % v5, state % v6)

END SUBROUTINE Deallocate_CVs
