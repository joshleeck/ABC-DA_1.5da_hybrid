SUBROUTINE Initialise_um_data (state)

! Initialise um-data variables

USE DefConsTypes, ONLY :   &
    UM_type,               &
    nlongs, nlevs

IMPLICIT NONE

TYPE(UM_type), INTENT(INOUT)   :: state

IF (.NOT.ALLOCATED(state % longs_u))        ALLOCATE(state % longs_u(1:nlongs))
IF (.NOT.ALLOCATED(state % longs_v))        ALLOCATE(state % longs_v(1:nlongs))
IF (.NOT.ALLOCATED(state % half_levs))      ALLOCATE(state % half_levs(1:nlevs+1))
IF (.NOT.ALLOCATED(state % full_levs))      ALLOCATE(state % full_levs(0:nlevs+1))
IF (.NOT.ALLOCATED(state % u))              ALLOCATE(state % u(1:nlongs,1:nlevs))
IF (.NOT.ALLOCATED(state % v))              ALLOCATE(state % v(1:nlongs,1:nlevs))
IF (.NOT.ALLOCATED(state % w))              ALLOCATE(state % w(1:nlongs,0:nlevs))
IF (.NOT.ALLOCATED(state % density))        ALLOCATE(state % density(1:nlongs,1:nlevs))
IF (.NOT.ALLOCATED(state % theta))          ALLOCATE(state % theta(1:nlongs,1:nlevs))
IF (.NOT.ALLOCATED(state % exner_pressure)) ALLOCATE(state % exner_pressure(1:nlongs,1:nlevs+1))
IF (.NOT.ALLOCATED(state % orog_height))    ALLOCATE(state % orog_height(1:nlongs))

state % longs_u(1:nlongs)                  = 0.0
state % longs_v(1:nlongs)                  = 0.0
state % half_levs(1:nlevs+1)               = 0.0
state % full_levs(0:nlevs+1)               = 0.0
state % u(1:nlongs,1:nlevs)                = 0.0
state % v(1:nlongs,1:nlevs)                = 0.0
state % w(1:nlongs,0:nlevs)                = 0.0
state % density(1:nlongs,1:nlevs)          = 0.0
state % theta(1:nlongs,1:nlevs)            = 0.0
state % exner_pressure(1:nlongs,1:nlevs+1) = 0.0
state % orog_height(1:nlongs)              = 0.0

END SUBROUTINE Initialise_um_data
