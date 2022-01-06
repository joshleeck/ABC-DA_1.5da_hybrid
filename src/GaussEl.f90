! -------------------------------------------------------------------------------
SUBROUTINE GaussEl ( N,      & ! Size of system
                     A,      & ! Matrix (input)
                     y,      & ! RHS (input)
                     x )       ! Solution (output)
! Gaussian elimination for a full square matrix
! Declare local variables

IMPLICIT NONE

INTEGER,      PARAMETER   :: ZREAL8 = SELECTED_REAL_KIND(15,307)

! Subroutine parameters
INTEGER,      INTENT(IN)  :: N
REAL(ZREAL8), INTENT(IN)  :: A(1:N, 1:N)
REAL(ZREAL8), INTENT(IN)  :: y(1:N)
REAL(ZREAL8), INTENT(OUT) :: x(1:N)

! Local variables
INTEGER                   :: i, j, jp
REAL(ZREAL8)              :: recip_Aii, Aii, Aji, partialsum
REAL(ZREAL8), ALLOCATABLE :: Acopy(:,:)
REAL(ZREAL8), ALLOCATABLE :: ycopy(:)


ALLOCATE (Acopy(1:N, 1:N))
ALLOCATE (ycopy(1:N))

Acopy(1:N, 1:N) = A(1:N, 1:N)
ycopy(1:N)      = y(1:N)

! First go through the matrix
DO i = 1, N-1
  recip_Aii = 1.0 / Acopy(i,i)
  ! Modify later rows accordingly
  DO j = i+1, N
    Aji = Acopy(j,i)
    ! Modify RHS
    ycopy(j) = ycopy(j) - Aji * ycopy(i) * recip_Aii;
    ! Modify matrix
    DO jp = i+1, N
      Acopy(j,jp) = Acopy(j,jp) - Aji * Acopy(i,jp) * recip_Aii
    END DO
  END DO
END DO

! Now solve the matrix backwards
DO i = N, 1, -1
  IF (i < N) THEN
    ! Compute partial sum
    partialsum = SUM(Acopy(i,i+1:N) * x(i+1:N))
  ELSE
    partialsum = 0.0;
  END IF
  Aii       = Acopy(i,i)
  recip_Aii = 1.0 / Aii
  x(i)      = (ycopy(i) - partialsum) * recip_Aii;
END DO


DEALLOCATE (Acopy, ycopy)

END SUBROUTINE GaussEl
