SUBROUTINE U_trans_alpha_adj (LS, ModelOut, ModelIn, CVTdata, dims, L_alpha)

! Code to perform the adjoint of the localisation transform: ModelOut = U* ModelIn

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


TYPE(ABC_type),  INTENT(IN)     :: LS
TYPE(ABC_type),  INTENT(INOUT)  :: ModelOut
TYPE(ABC_type),  INTENT(IN)     :: ModelIn
TYPE(CVT_type),  INTENT(IN)     :: CVTdata
TYPE(aCVT_type), INTENT(IN)     :: L_alpha
TYPE(dims_type), INTENT(IN)     :: dims

TYPE(ABC_type)                  :: Interim1


CALL Initialise_model_vars (Interim1, .FALSE.)

CALL U_v_alpha_adj (Interim1, ModelIn, L_alpha)

CALL U_h_alpha_adj (ModelOut, Interim1, L_alpha)

CALL Deallocate_model_vars (Interim1)

END SUBROUTINE U_trans_alpha_adj
