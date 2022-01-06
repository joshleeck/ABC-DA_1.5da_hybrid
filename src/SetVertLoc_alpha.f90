SUBROUTINE SetVertLoc_alpha (vScale_alpha, dims, L_alpha)

! Code to compute the eigenvectors and eigenvalues for the vertical localisation matrix, write to L_alpha to store

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  aCVT_type,               &
  dims_type,              &
  nlongs,                 &
  nlevs,                  &
  small, ConditionFudge


IMPLICIT NONE

REAL(ZREAL8),    INTENT(IN)     :: vScale_alpha
TYPE(dims_type), INTENT(IN)     :: dims
TYPE(aCVT_type),  INTENT(INOUT)  :: L_alpha

INTEGER                      :: z1, z2, lim, ierr, k
REAL(ZREAL8)                 :: increment, zdiff
REAL(ZREAL8), ALLOCATABLE    :: Work(:), Lv(:,:)
REAL(ZREAL8), EXTERNAL       :: Fn_GC_corr

ALLOCATE (Lv(1:nlevs, 1:nlevs))
ALLOCATE (Work(1:3*nlevs))

! Compute Gaspari Cohn localisation function
DO z1 = 1, nlevs 
  DO z2 = 1, nlevs
    zdiff = ABS(dims % full_levs(z2) &
               - dims % full_levs(z1))

    Lv(z1,z2) = Fn_GC_corr(zdiff, vScale_alpha)
  END DO
END DO

! Note here we assume uniform inner product, i.e. just decompose L_alpha
! Also same localisation matrix for each control variable
! Horizontal localisation before vertical, so only 1 covariance matrix
lim = 1
DO k = 1, lim

  L_alpha % VertMode1(1:nlevs,1:nlevs,k) = Lv(1:nlevs,1:nlevs)
  L_alpha % VertMode2(1:nlevs,1:nlevs,k) = Lv(1:nlevs,1:nlevs)
  L_alpha % VertMode3(1:nlevs,1:nlevs,k) = Lv(1:nlevs,1:nlevs)
  L_alpha % VertMode4(1:nlevs,1:nlevs,k) = Lv(1:nlevs,1:nlevs)
  L_alpha % VertMode5(1:nlevs,1:nlevs,k) = Lv(1:nlevs,1:nlevs)


  CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
             'U',                                  & ! Upper triangular of matrix specified
              nlevs,                                & ! Order of the matrix to be diagonalised
              L_alpha % VertMode1(1:nlevs, 1:nlevs, k), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
              nlevs,                                & ! Leading dimension of the matrix
              L_alpha % VertEV1(1:nlevs,k),             & ! Eigenvalue array
              Work(1:3*nlevs),                      & ! Work array
              3*nlevs,                              & ! Size of work array
              ierr)

  DO z1 = 1, nlevs
    IF (L_alpha % VertEV1(z1,k) < 0) THEN
      L_alpha % VertEV1(z1,k) = 0
    END IF
  END DO
  ! Condition for zero values
  !increment = SUM(L_alpha % VertEV1(1:nlevs,k)) * ConditionFudge
  !IF (increment < small) increment = small
  !  L_alpha % VertEV1(1:nlevs,k) = L_alpha % VertEV1(1:nlevs,k) + increment
  L_alpha % VertEV1(1:nlevs,k) = SQRT(L_alpha % VertEV1(1:nlevs,k))

  CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
             'U',                                  & ! Upper triangular of matrix specified
              nlevs,                                & ! Order of the matrix to be diagonalised
              L_alpha % VertMode2(1:nlevs, 1:nlevs, k), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
              nlevs,                                & ! Leading dimension of the matrix
              L_alpha % VertEV2(1:nlevs,k),             & ! Eigenvalue array
              Work(1:3*nlevs),                      & ! Work array
              3*nlevs,                              & ! Size of work array
              ierr)

  DO z1 = 1, nlevs
    IF (L_alpha % VertEV2(z1,k) < 0) THEN
      L_alpha % VertEV2(z1,k) = 0
    END IF
  END DO
  ! Condition for zero values
  !increment = SUM(L_alpha % VertEV2(1:nlevs,k)) * ConditionFudge
  !IF (increment < small) increment = small
  !  L_alpha % VertEV2(1:nlevs,k) = L_alpha % VertEV2(1:nlevs,k) + increment
  L_alpha % VertEV2(1:nlevs,k) = SQRT(L_alpha % VertEV2(1:nlevs,k))

  CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
             'U',                                  & ! Upper triangular of matrix specified
              nlevs,                                & ! Order of the matrix to be diagonalised
              L_alpha % VertMode3(1:nlevs, 1:nlevs, k), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
              nlevs,                                & ! Leading dimension of the matrix
              L_alpha % VertEV3(1:nlevs,k),             & ! Eigenvalue array
              Work(1:3*nlevs),                      & ! Work array
              3*nlevs,                              & ! Size of work array
              ierr)

  DO z1 = 1, nlevs
    IF (L_alpha % VertEV3(z1,k) < 0) THEN
      L_alpha % VertEV3(z1,k) = 0
    END IF
  END DO
  ! Condition for zero values
  !increment = SUM(L_alpha % VertEV3(1:nlevs,k)) * ConditionFudge
  !IF (increment < small) increment = small
  !  L_alpha % VertEV3(1:nlevs,k) = L_alpha % VertEV3(1:nlevs,k) + increment
  L_alpha % VertEV3(1:nlevs,k) = SQRT(L_alpha % VertEV3(1:nlevs,k))

  CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
             'U',                                  & ! Upper triangular of matrix specified
              nlevs,                                & ! Order of the matrix to be diagonalised
              L_alpha % VertMode4(1:nlevs, 1:nlevs, k), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
              nlevs,                                & ! Leading dimension of the matrix
              L_alpha % VertEV4(1:nlevs,k),             & ! Eigenvalue array
              Work(1:3*nlevs),                      & ! Work array
              3*nlevs,                              & ! Size of work array
              ierr)

  DO z1 = 1, nlevs
    IF (L_alpha % VertEV4(z1,k) < 0) THEN
      L_alpha % VertEV4(z1,k) = 0
    END IF
  END DO
  ! Condition for zero values
  !increment = SUM(L_alpha % VertEV4(1:nlevs,k)) * ConditionFudge
  !IF (increment < small) increment = small
  !  L_alpha % VertEV4(1:nlevs,k) = L_alpha % VertEV4(1:nlevs,k) + increment
  L_alpha % VertEV4(1:nlevs,k) = SQRT(L_alpha % VertEV4(1:nlevs,k))

  CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
             'U',                                  & ! Upper triangular of matrix specified
              nlevs,                                & ! Order of the matrix to be diagonalised
              L_alpha % VertMode5(1:nlevs, 1:nlevs, k), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
              nlevs,                                & ! Leading dimension of the matrix
              L_alpha % VertEV5(1:nlevs,k),             & ! Eigenvalue array
              Work(1:3*nlevs),                      & ! Work array
              3*nlevs,                              & ! Size of work array
              ierr)

  DO z1 = 1, nlevs
    IF (L_alpha % VertEV5(z1,k) < 0) THEN
      L_alpha % VertEV5(z1,k) = 0
    END IF
  END DO
  ! Condition for zero values
  !increment = SUM(L_alpha % VertEV5(1:nlevs,k)) * ConditionFudge
  !IF (increment < small) increment = small
  !  L_alpha % VertEV5(1:nlevs,k) = L_alpha % VertEV5(1:nlevs,k) + increment
  L_alpha % VertEV5(1:nlevs,k) = SQRT(L_alpha % VertEV5(1:nlevs,k))

END DO

DEALLOCATE (Lv)
DEALLOCATE (Work)

END SUBROUTINE SetVertLoc_alpha
