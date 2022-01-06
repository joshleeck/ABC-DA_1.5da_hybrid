SUBROUTINE Deallocate_dims (state)

USE DefConsTypes, ONLY :   &
    dims_type

IMPLICIT NONE

TYPE(dims_type), INTENT(INOUT)   :: state

DEALLOCATE (state % longs_u, state % longs_v)
DEALLOCATE (state % half_levs, state % full_levs)
DEALLOCATE (state % a1, state % b1, state % a2, state % b2)
DEALLOCATE (state % recip_half_kp1_k, state % recip_half_k_km1)
DEALLOCATE (state % recip_full_kp1_k, state % recip_full_k_km1)

END SUBROUTINE Deallocate_dims
