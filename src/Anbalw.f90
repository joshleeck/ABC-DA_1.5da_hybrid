SUBROUTINE Anbalw (LS, du, drhop, w_b, dims)

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
REAL(ZREAL8),    INTENT(IN)    :: du(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),    INTENT(IN)    :: drhop(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),    INTENT(INOUT) :: w_b(1:nlongs,1:nlevs)
TYPE(dims_type), INTENT(IN)    :: dims

! Local variables
INTEGER                        :: x, z
REAL(ZREAL8), ALLOCATABLE      :: term1(:), term2(:), term3(:), diff(:), rhodw(:)

! Function
REAL(ZREAL8)                   :: INT_HF



ALLOCATE (term1(1:nlevs))
ALLOCATE (term2(1:nlevs))
ALLOCATE (term3(1:nlevs))
ALLOCATE (diff(1:nlevs))
ALLOCATE (rhodw(0:nlevs))


DO x = 1, nlongs

  ! Calculate terms on half levels

  ! Term 1: d/dx (deltarho' u)
  ! --------------------------
  DO z = 1, nlevs
    term1(z) = ((drhop(x,z)   + drhop(x+1,z)) * LS % u(x,z) -      &
                (drhop(x-1,z) + drhop(x,z))   * LS % u(x-1,z)) *   &
               recip2dx
  END DO

  ! Term 2: d/dz (deltarho' w)
  ! --------------------------
  DO z = 1, nlevs
    term2(z) = (INT_HF(drhop(x,z),   drhop(x,z+1), z,   dims) * LS % w(x,z) -      &
                INT_HF(drhop(x,z-1), drhop(x,z),   z-1, dims) * LS % w(x,z-1)) / &
               (dims % full_levs(z) - dims % full_levs(z-1))
  END DO

  ! Term 3: d/dx (rho du)
  ! ---------------------
  DO z = 1, nlevs
    term3(z) = ((LS % rho(x,z) + LS % rho(x+1,z)) * du(x,z) -      &
                (LS % rho(x-1,z) + LS % rho(x,z)) * du(x-1,z)) *   &
               recip2dx
  END DO

  ! Add these (and negate)
  ! ----------------------
  diff(1:nlevs) = -1.0 * term1(1:nlevs) - term2(1:nlevs) - term3(1:nlevs)
  ! This is held on half levels


  ! Integrate from the bottom of the domain
  ! ---------------------------------------
  rhodw(0) = 0.0
  DO z = 1, nlevs
    rhodw(z) = rhodw(z-1) + diff(z) * (dims % full_levs(z) - dims % full_levs(z-1))
  END DO

  ! Compute balanced w by dividing by density
  ! -----------------------------------------
  DO z = 1, nlevs
    w_b(x,z) = rhodw(z) / (INT_HF(LS % rho(x,z), LS % rho(x,z+1), z, dims))
  END DO

END DO


DEALLOCATE (term1, term2, term3, diff, rhodw)

END SUBROUTINE Anbalw
