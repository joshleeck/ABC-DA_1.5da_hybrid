SUBROUTINE Mul_CV_cons(state, mul_cons)

! Multiply a control variable by a constant

USE DefConsTypes, ONLY :     &
    ZREAL8,                  &
    CV_type,                 &
    nlongs, nlevs

IMPLICIT NONE

TYPE(CV_type),  INTENT(INOUT)   :: state
REAL(ZREAL8),   INTENT(IN)      :: mul_cons

state % v1(0:nlongs+1,0:nlevs+1) = state % v1(0:nlongs+1,0:nlevs+1) * mul_cons
state % v2(0:nlongs+1,0:nlevs+1) = state % v2(0:nlongs+1,0:nlevs+1) * mul_cons
state % v3(0:nlongs+1,0:nlevs+1) = state % v3(0:nlongs+1,0:nlevs+1) * mul_cons
state % v4(0:nlongs+1,0:nlevs+1) = state % v4(0:nlongs+1,0:nlevs+1) * mul_cons
state % v5(0:nlongs+1,0:nlevs+1) = state % v5(0:nlongs+1,0:nlevs+1) * mul_cons
state % v6(0:nlongs+1,0:nlevs+1) = state % v6(0:nlongs+1,0:nlevs+1) * mul_cons

END SUBROUTINE Mul_CV_cons
