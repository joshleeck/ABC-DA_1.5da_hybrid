SUBROUTINE U_v_alpha_adj (ControlVar1, ControlVar2, L_alpha)

! Code to perform the adjoint vertical localisation: ControlVar1 = U_v_alpha* ControlVar2

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,                &
  aCVT_type,               &
  nlongs,                 &
  nlevs


IMPLICIT NONE

TYPE(ABC_type),  INTENT(INOUT) :: ControlVar1
TYPE(ABC_type),  INTENT(IN)    :: ControlVar2
TYPE(aCVT_type), INTENT(IN)    :: L_alpha
INTEGER                        :: x
TYPE(ABC_type)                 :: Interim

CALL Initialise_model_vars (Interim, .FALSE.)

! Do transform from levels to vertical mode space
DO x = 1, nlongs
  Interim % u(x,1:nlevs) = MATMUL(TRANSPOSE(L_alpha % VertMode1(1:nlevs,1:nlevs,1)), ControlVar2 % u(x,1:nlevs))
  Interim % v(x,1:nlevs) = MATMUL(TRANSPOSE(L_alpha % VertMode2(1:nlevs,1:nlevs,1)), ControlVar2 % v(x,1:nlevs))
  Interim % w(x,1:nlevs) = MATMUL(TRANSPOSE(L_alpha % VertMode3(1:nlevs,1:nlevs,1)), ControlVar2 % w(x,1:nlevs))
  Interim % r(x,1:nlevs) = MATMUL(TRANSPOSE(L_alpha % VertMode4(1:nlevs,1:nlevs,1)), ControlVar2 % r(x,1:nlevs))
  Interim % b(x,1:nlevs) = MATMUL(TRANSPOSE(L_alpha % VertMode5(1:nlevs,1:nlevs,1)), ControlVar2 % b(x,1:nlevs))
END DO

! Multiply by the evs (these are actually the square-roots)
! These are factors in vert mode space
DO x = 1, nlongs
  Interim % u(x,1:nlevs) = Interim % u(x,1:nlevs) * L_alpha % VertEV1(1:nlevs,1)
  Interim % v(x,1:nlevs) = Interim % v(x,1:nlevs) * L_alpha % VertEV2(1:nlevs,1)
  Interim % w(x,1:nlevs) = Interim % w(x,1:nlevs) * L_alpha % VertEV3(1:nlevs,1)
  Interim % r(x,1:nlevs) = Interim % r(x,1:nlevs) * L_alpha % VertEV4(1:nlevs,1)
  Interim % b(x,1:nlevs) = Interim % b(x,1:nlevs) * L_alpha % VertEV5(1:nlevs,1)
END DO

! Do non-symmetric transform - keep in vert mode space
ControlVar1 % u(1:nlongs,1:nlevs) = Interim % u(1:nlongs,1:nlevs)
ControlVar1 % v(1:nlongs,1:nlevs) = Interim % v(1:nlongs,1:nlevs)
ControlVar1 % w(1:nlongs,1:nlevs) = Interim % w(1:nlongs,1:nlevs)
ControlVar1 % r(1:nlongs,1:nlevs) = Interim % r(1:nlongs,1:nlevs)
ControlVar1 % b(1:nlongs,1:nlevs) = Interim % b(1:nlongs,1:nlevs)

CALL Deallocate_model_vars (Interim)

END SUBROUTINE U_v_alpha_adj
