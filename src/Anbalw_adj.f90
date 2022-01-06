SUBROUTINE Anbalw_adj (LS, du, drhop, w_b, dims)

! Code to compute the anelastically balanced component of w


USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  dims_type,              &
  ABC_type,               &
  nlongs,                 &
  nlevs,                  &
  recip2dx
  

IMPLICIT NONE

! Subroutine parameters
TYPE(ABC_type),  INTENT(IN)    :: LS
REAL(ZREAL8),    INTENT(INOUT) :: du(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),    INTENT(INOUT) :: drhop(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),    INTENT(IN)    :: w_b(1:nlongs,1:nlevs)
TYPE(dims_type), INTENT(IN)    :: dims

! Local variables
INTEGER                        :: x, z
REAL(ZREAL8), ALLOCATABLE      :: term1(:), term2(:), term3(:), diff(:), rhodw(:)
REAL(ZREAL8)                   :: recipdz

! Function
REAL(ZREAL8)                   :: INT_HF



ALLOCATE (term1(1:nlevs))
ALLOCATE (term2(1:nlevs))
ALLOCATE (term3(1:nlevs))
ALLOCATE (diff(1:nlevs+1))
ALLOCATE (rhodw(0:nlevs))


DO x = 1, nlongs

  ! Calculate terms on half levels

  ! Compute balanced w by dividing by density
  ! -----------------------------------------
  DO z = 1, nlevs
    rhodw(z) = w_b(x,z) / (INT_HF(LS % rho(x,z), LS % rho(x,z+1), z, dims))
  END DO

  ! Integrate from the bottom of the domain
  ! ---------------------------------------
  ! Adjoint integrates from top to bottom
  diff(nlevs+1) = 0.0
  DO z = nlevs, 1, -1
    diff(z) = diff(z+1) + rhodw(z) * (dims % full_levs(z) - dims % full_levs(z-1))
  END DO

  ! Add these (and negate)
  ! ----------------------
  term1(1:nlevs) = -1.0 * diff(1:nlevs)
  term2(1:nlevs) = -1.0 * diff(1:nlevs)
  term3(1:nlevs) = -1.0 * diff(1:nlevs)

  ! Term 3: d/dx (rho du)
  ! ---------------------
  DO z = 1, nlevs
    du(x,z)   = du(x,z) + (LS % rho(x,z) + LS % rho(x+1,z)) * recip2dx * term3(z)
    du(x-1,z) = du(x-1,z) - (LS % rho(x-1,z) + LS % rho(x,z)) * recip2dx * term3(z)
  END DO

  ! Term 2: d/dz (deltarho' w)
  ! --------------------------
  DO z = 1, nlevs
    ! This is the forward part of the code, evaluating the function INT_HF
    ! This will help construction of the adjoint
    !term2(z) = ((dims % a2(z) * drhop(x,z+1) + dims % b2(z) * drhop(x,z)) * LS % w(x,z) -      &
    !            (dims % a2(z-1) * drhop(x,z) + dims % b2(z-1) * drhop(x,z-1)) * LS % w(x,z-1)) / &
    !           (dims % full_levs(z) - dims % full_levs(z-1))

    recipdz = 1.0 / (dims % full_levs(z) - dims % full_levs(z-1))
    drhop(x,z+1) = drhop(x,z+1) + dims % a2(z) * LS % w(x,z)     * recipdz * term2(z)
    drhop(x,z)   = drhop(x,z)   + dims % b2(z) * LS % w(x,z)     * recipdz * term2(z)
    drhop(x,z)   = drhop(x,z)   - dims % a2(z-1) * LS % w(x,z-1) * recipdz * term2(z)
    drhop(x,z-1) = drhop(x,z-1) - dims % b2(z-1) * LS % w(x,z-1) * recipdz * term2(z)

  END DO

  ! Term 1: d/dx (deltarho' u)
  ! --------------------------
  DO z = 1, nlevs
    drhop(x,z)   = drhop(x,z)   + LS % u(x,z)   * recip2dx * term1(z)
    drhop(x+1,z) = drhop(x+1,z) + LS % u(x,z)   * recip2dx * term1(z)
    drhop(x-1,z) = drhop(x-1,z) - LS % u(x-1,z) * recip2dx * term1(z)
    drhop(x,z)   = drhop(x,z)   - LS % u(x-1,z) * recip2dx * term1(z)
  END DO

END DO


DEALLOCATE (term1, term2, term3, diff, rhodw)

END SUBROUTINE Anbalw_adj
