SUBROUTINE Initialise_model_vars(state, random_init)

! Initialise model variables

USE DefConsTypes, ONLY :     &
    ABC_type,                &
    nlongs, nlevs

IMPLICIT NONE

INCLUDE "Boundaries.interface"

TYPE(ABC_type), INTENT(INOUT)   :: state
LOGICAL,        INTENT(IN)      :: random_init



IF (.NOT.ALLOCATED(state % u))               ALLOCATE(state % u(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % v))               ALLOCATE(state % v(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % w))               ALLOCATE(state % w(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % r))               ALLOCATE(state % r(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % b))               ALLOCATE(state % b(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % tracer))          ALLOCATE(state % tracer(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % rho))             ALLOCATE(state % rho(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % b_ef))            ALLOCATE(state % b_ef(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % hydro_imbal))     ALLOCATE(state % hydro_imbal(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % geost_imbal))     ALLOCATE(state % geost_imbal(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % vert_mom_source)) ALLOCATE(state % vert_mom_source(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % horiz_div))       ALLOCATE(state % horiz_div(0:nlongs+1,0:nlevs+1))
IF (.NOT.ALLOCATED(state % horiz_vort))      ALLOCATE(state % horiz_vort(0:nlongs+1,0:nlevs+1))

IF (random_init) THEN

  CALL RANDOM_NUMBER (state % u(1:nlongs,1:nlevs))
  CALL RANDOM_NUMBER (state % v(1:nlongs,1:nlevs))
  CALL RANDOM_NUMBER (state % w(1:nlongs,1:nlevs))
  CALL RANDOM_NUMBER (state % r(1:nlongs,1:nlevs))
  CALL RANDOM_NUMBER (state % b(1:nlongs,1:nlevs))
  CALL RANDOM_NUMBER (state % tracer(1:nlongs,1:nlevs))
  state % rho(1:nlongs,1:nlevs) = 1.0 + state % r(1:nlongs,1:nlevs)

  CALL Boundaries (state)

ELSE

  state % u(0:nlongs+1,0:nlevs+1)           = 0.0
  state % v(0:nlongs+1,0:nlevs+1)           = 0.0
  state % w(0:nlongs+1,0:nlevs+1)           = 0.0
  state % r(0:nlongs+1,0:nlevs+1)           = 0.0     ! density perturbation
  state % b(0:nlongs+1,0:nlevs+1)           = 0.0     ! buoyancy perturbation
  state % rho(0:nlongs+1,0:nlevs+1)         = 0.0     ! density full field
  state % b_ef(0:nlongs+1,0:nlevs+1)        = 0.0     ! Effective buoyancy
  state % tracer(0:nlongs+1,0:nlevs+1)      = 0.0
  state % hydro_imbal(0:nlongs+1,0:nlevs+1) = 0.0
  state % geost_imbal(0:nlongs+1,0:nlevs+1) = 0.0
  state % Kinetic_Energy                    = 0.0
  state % Buoyant_Energy                    = 0.0
  state % Elastic_Energy                    = 0.0
  state % Total_Energy                      = 0.0

END IF

END SUBROUTINE Initialise_model_vars
