SUBROUTINE Initialise_dims (state)

! Initialise dimension variables

USE DefConsTypes, ONLY :   &
    dims_type,             &
    nlongs, nlevs

IMPLICIT NONE

TYPE(dims_type), INTENT(INOUT)   :: state

IF (.NOT.ALLOCATED(state % longs_u))          ALLOCATE (state % longs_u(0:nlongs+1))
IF (.NOT.ALLOCATED(state % longs_v))          ALLOCATE (state % longs_v(0:nlongs+1))
IF (.NOT.ALLOCATED(state % half_levs))        ALLOCATE (state % half_levs(0:nlevs+1))
IF (.NOT.ALLOCATED(state % full_levs))        ALLOCATE (state % full_levs(0:nlevs+1))
IF (.NOT.ALLOCATED(state % a1))               ALLOCATE (state % a1(1:nlevs+1)) 
IF (.NOT.ALLOCATED(state % b1))               ALLOCATE (state % b1(1:nlevs+1))
IF (.NOT.ALLOCATED(state % a2))               ALLOCATE (state % a2(0:nlevs))
IF (.NOT.ALLOCATED(state % b2))               ALLOCATE (state % b2(0:nlevs))
IF (.NOT.ALLOCATED(state % recip_half_kp1_k)) ALLOCATE (state % recip_half_kp1_k(0:nlevs))
IF (.NOT.ALLOCATED(state % recip_half_k_km1)) ALLOCATE (state % recip_half_k_km1(1:nlevs+1))
IF (.NOT.ALLOCATED(state % recip_full_kp1_k)) ALLOCATE (state % recip_full_kp1_k(0:nlevs))
IF (.NOT.ALLOCATED(state % recip_full_k_km1)) ALLOCATE (state % recip_full_k_km1(1:nlevs+1))

state % longs_u(0:nlongs+1)         = 0.0
state % longs_v(0:nlongs+1)         = 0.0
state % half_levs(0:nlevs+1)        = 0.0
state % full_levs(0:nlevs+1)        = 0.0
state % a1(1:nlevs+1)               = 0.0
state % b1(1:nlevs+1)               = 0.0
state % a2(0:nlevs)                 = 0.0
state % b2(0:nlevs)                 = 0.0
state % recip_half_kp1_k(0:nlevs)   = 0.0
state % recip_half_k_km1(1:nlevs+1) = 0.0
state % recip_full_kp1_k(0:nlevs)   = 0.0
state % recip_full_k_km1(1:nlevs+1) = 0.0

END SUBROUTINE Initialise_dims
