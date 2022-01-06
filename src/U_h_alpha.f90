SUBROUTINE U_h_alpha (ControlVar1, ControlVar2, L_alpha)

! Code to perform the horizontal localisation: ControlVar2 = U_h_alpha* ControlVar1

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  aCVT_type,              &
  nlongs,                 &
  nlevs 


IMPLICIT NONE

TYPE(ABC_type),  INTENT(IN)    :: ControlVar1
TYPE(ABC_type),  INTENT(INOUT) :: ControlVar2
TYPE(aCVT_type), INTENT(IN)    :: L_alpha

INTEGER                        :: z
TYPE(ABC_type)                 :: Interim

CALL Initialise_model_vars (Interim, .FALSE.)

Interim % u(1:nlongs,1:nlevs) = ControlVar1 % u(1:nlongs,1:nlevs)
Interim % v(1:nlongs,1:nlevs) = ControlVar1 % v(1:nlongs,1:nlevs)
Interim % w(1:nlongs,1:nlevs) = ControlVar1 % w(1:nlongs,1:nlevs)
Interim % r(1:nlongs,1:nlevs) = ControlVar1 % r(1:nlongs,1:nlevs)
Interim % b(1:nlongs,1:nlevs) = ControlVar1 % b(1:nlongs,1:nlevs)


! Multiply by the evs (these are actually the square-roots)
! These are factors in horiz mode space
DO z = 1, nlevs
  Interim % u(1:nlongs,z) = L_alpha % HorizEV1(1:nlongs,1) * Interim % u(1:nlongs,z)
  Interim % v(1:nlongs,z) = L_alpha % HorizEV2(1:nlongs,1) * Interim % v(1:nlongs,z)
  Interim % w(1:nlongs,z) = L_alpha % HorizEV3(1:nlongs,1) * Interim % w(1:nlongs,z)
  Interim % r(1:nlongs,z) = L_alpha % HorizEV4(1:nlongs,1) * Interim % r(1:nlongs,z)
  Interim % b(1:nlongs,z) = L_alpha % HorizEV5(1:nlongs,1) * Interim % b(1:nlongs,z)
END DO

! Do transform from horiz mode space to gridpoint space
DO z = 1, nlevs
  ControlVar2 % u(1:nlongs,z) = MATMUL(L_alpha % HorizMode1(1:nlongs,1:nlongs,1), Interim % u(1:nlongs,z))
  ControlVar2 % v(1:nlongs,z) = MATMUL(L_alpha % HorizMode2(1:nlongs,1:nlongs,1), Interim % v(1:nlongs,z))
  ControlVar2 % w(1:nlongs,z) = MATMUL(L_alpha % HorizMode3(1:nlongs,1:nlongs,1), Interim % w(1:nlongs,z))
  ControlVar2 % r(1:nlongs,z) = MATMUL(L_alpha % HorizMode4(1:nlongs,1:nlongs,1), Interim % r(1:nlongs,z))
  ControlVar2 % b(1:nlongs,z) = MATMUL(L_alpha % HorizMode5(1:nlongs,1:nlongs,1), Interim % b(1:nlongs,z))
END DO

CALL Deallocate_model_vars (Interim) 

END SUBROUTINE U_h_alpha
