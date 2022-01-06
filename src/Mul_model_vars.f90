SUBROUTINE Mul_model_vars(state, mul_state, core_only)

! Multiply states element-wise, essentially a Schur product

USE DefConsTypes, ONLY :     &
    ABC_type,                &
    nlongs, nlevs

IMPLICIT NONE

TYPE(ABC_type), INTENT(INOUT)   :: state    ! Running total
TYPE(ABC_type), INTENT(IN)      :: mul_state
LOGICAL,        INTENT(IN)      :: core_only


  state % u(0:nlongs+1,0:nlevs+1)           = state % u(0:nlongs+1,0:nlevs+1) * &
                                              mul_state % u(0:nlongs+1,0:nlevs+1)
  state % v(0:nlongs+1,0:nlevs+1)           = state % v(0:nlongs+1,0:nlevs+1) * &
                                              mul_state % v(0:nlongs+1,0:nlevs+1)
  state % w(0:nlongs+1,0:nlevs+1)           = state % w(0:nlongs+1,0:nlevs+1) * &
                                              mul_state % w(0:nlongs+1,0:nlevs+1)
  state % r(0:nlongs+1,0:nlevs+1)           = state % r(0:nlongs+1,0:nlevs+1) * &
                                              mul_state % r(0:nlongs+1,0:nlevs+1)
  state % b(0:nlongs+1,0:nlevs+1)           = state % b(0:nlongs+1,0:nlevs+1) * &
                                              mul_state % b(0:nlongs+1,0:nlevs+1)

  IF (.NOT.core_only) THEN
    state % rho(0:nlongs+1,0:nlevs+1)         = state % rho(0:nlongs+1,0:nlevs+1) * &
                                                mul_state % rho(0:nlongs+1,0:nlevs+1)
    state % b_ef(0:nlongs+1,0:nlevs+1)        = state % b_ef(0:nlongs+1,0:nlevs+1) * &
                                                mul_state % b_ef(0:nlongs+1,0:nlevs+1)
    state % tracer(0:nlongs+1,0:nlevs+1)      = state % tracer(0:nlongs+1,0:nlevs+1) * &
                                                mul_state % tracer(0:nlongs+1,0:nlevs+1)
    state % hydro_imbal(0:nlongs+1,0:nlevs+1) = state % hydro_imbal(0:nlongs+1,0:nlevs+1) * &
                                                mul_state % hydro_imbal(0:nlongs+1,0:nlevs+1)
    state % geost_imbal(0:nlongs+1,0:nlevs+1) = state % geost_imbal(0:nlongs+1,0:nlevs+1) * &
                                                mul_state % geost_imbal(0:nlongs+1,0:nlevs+1)
    state % Kinetic_Energy                    = state % Kinetic_Energy * &
                                                mul_state % Kinetic_Energy
    state % Buoyant_Energy                    = state % Buoyant_Energy * &
                                                mul_state % Buoyant_Energy
    state % Elastic_Energy                    = state % Elastic_Energy * &
                                                mul_state % Elastic_Energy
    state % Total_Energy                      = state % Total_Energy * &
                                                mul_state % Total_Energy
  END IF

END SUBROUTINE Mul_model_vars
