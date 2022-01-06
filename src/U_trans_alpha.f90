SUBROUTINE U_trans_alpha (LS, ModelIn, ModelOut, CVTdata, dims, L_alpha)

! Code to perform the localisation transforms: ModelVars = U ControlVars

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  CV_type,                &
  CVT_type,               &
  aCVT_type,              &
  dims_type,              &
  nlongs,                 &
  nlevs


IMPLICIT NONE



TYPE(ABC_type),  INTENT(IN)    :: LS
TYPE(ABC_type),  INTENT(IN)    :: ModelIn
TYPE(ABC_type),  INTENT(INOUT) :: ModelOut
TYPE(CVT_type),  INTENT(IN)    :: CVTdata
TYPE(aCVT_type), INTENT(IN)    :: L_alpha
TYPE(dims_type), INTENT(IN)    :: dims

TYPE(ABC_type)                 :: Interim1


CALL Initialise_model_vars (Interim1, .FALSE.)

CALL U_h_alpha (ModelIn, Interim1, L_alpha)

CALL U_v_alpha (Interim1, ModelOut, L_alpha) 

CALL Deallocate_model_vars (Interim1)

END SUBROUTINE U_trans_alpha
