SUBROUTINE SetHorizLoc_alpha (hScale_alpha, dims, L_alpha)

! Code to compute the eigenvectors and eigenvalues for the horizontal localisation matrix, write to L_alpha to store

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  aCVT_type,               &
  dims_type,              &
  nlongs,                 &
  nlevs, dx,              &
  ConditionFudge,         &
  small


IMPLICIT NONE

REAL(ZREAL8),    INTENT(IN)     :: hScale_alpha
TYPE(dims_type), INTENT(IN)     :: dims
TYPE(aCVT_type),  INTENT(INOUT)  :: L_alpha

INTEGER                      :: k, x1, x2, lev, lim, ierr
REAL(ZREAL8)                 :: alpha, increment, xdiff1, xdiff2, xdiff
REAL(ZREAL8), ALLOCATABLE    :: Work(:), Func(:,:), Interim(:,:), Lh(:,:)
REAL(ZREAL8), EXTERNAL       :: Fn_GC_corr

ALLOCATE (Func(1:nlongs, 1:nlevs))
ALLOCATE (Interim(1:nlongs, 1:nlevs))
ALLOCATE (Lh(1:nlongs, 1:nlongs))
ALLOCATE (Work(1:3*nlongs))

! Compute Gaspari Cohn localisation function
DO x1 = 1, nlongs
  DO x2 = 1, nlongs
    ! Need to account for 'two distances' between any two points because of periodic BC and construct circulant Lh 
    xdiff1 = ABS(x2 - x1)*dx
    IF (xdiff1 <= 0.5*dims % longs_u(nlongs)) THEN
      Lh(x1,x2) = Fn_GC_corr(xdiff1, hScale_alpha)
      !Lh(x1,x2) = EXP(-0.5*(xdiff1/hScale_alpha)**2)
    ELSE
      IF (x2 > x1) THEN
        xdiff2 = (x1 + (nlongs-x2))*dx
      ELSE
        xdiff2 = (x2 + (nlongs-x1))*dx
      END IF
    Lh(x1,x2) = Fn_GC_corr(xdiff2, hScale_alpha)
    !Lh(x1,x2) = EXP(-0.5*(xdiff2/hScale_alpha)**2)
    END IF
  END DO
END DO


! Check for symmetry
!DO x1 = 1, nlongs
!  DO x2 = 1, nlongs
!    IF (Lh(x1,x2) .NE. Lh(x2,x1)) THEN
!      PRINT *, Lh(x1,x2)
!    END IF
!  END DO
!END DO

! Note here we assume uniform inner product, i.e. just decompose L_alpha
! Also same localisation matrix for each control variable
! Horizontal localisation before vertical, so only 1 covariance matrix
lim = 1
DO lev = 1, lim

  L_alpha % HorizMode1(1:nlongs,1:nlongs,lev) = Lh(1:nlongs,1:nlongs)
  L_alpha % HorizMode2(1:nlongs,1:nlongs,lev) = Lh(1:nlongs,1:nlongs)
  L_alpha % HorizMode3(1:nlongs,1:nlongs,lev) = Lh(1:nlongs,1:nlongs)
  L_alpha % HorizMode4(1:nlongs,1:nlongs,lev) = Lh(1:nlongs,1:nlongs)
  L_alpha % HorizMode5(1:nlongs,1:nlongs,lev) = Lh(1:nlongs,1:nlongs)


  CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
             'U',                                  & ! Upper triangular of matrix specified
              nlongs,                                & ! Order of the matrix to be diagonalised
              L_alpha % HorizMode1(1:nlongs, 1:nlongs, lev), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
              nlongs,                                & ! Leading dimension of the matrix
              L_alpha % HorizEV1(1:nlongs,lev),             & ! Eigenvalue array
              Work(1:3*nlongs),                      & ! Work array
              3*nlongs,                              & ! Size of work array
              ierr)

  !PRINT *, L_alpha % HorizEV1(:, lev)
  ! As hScale_alpha tends to infinity, circulant Lh is no longer positive semi-definite, need to truncate negative eigenvalues
  DO x1 = 1, nlongs
    IF (L_alpha % HorizEV1(x1,lev) < 0) THEN
      L_alpha % HorizEV1(x1,lev) = 0
    END IF
  END DO
  ! Condition for zero values
  !increment = SUM(L_alpha % HorizEV1(1:nlongs,lev)) * ConditionFudge
  !IF (increment < small) increment = small
  !  L_alpha % HorizEV1(1:nlongs,lev) = L_alpha % HorizEV1(1:nlongs,lev) + increment
  L_alpha % HorizEV1(1:nlongs,lev) = SQRT(L_alpha % HorizEV1(1:nlongs,lev))

  CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
             'U',                                  & ! Upper triangular of matrix specified
              nlongs,                                & ! Order of the matrix to be diagonalised
              L_alpha % HorizMode2(1:nlongs, 1:nlongs, lev), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
              nlongs,                                & ! Leading dimension of the matrix
              L_alpha % HorizEV2(1:nlongs,lev),             & ! Eigenvalue array
              Work(1:3*nlongs),                      & ! Work array
              3*nlongs,                              & ! Size of work array
              ierr)

  ! As hScale_alpha tends to infinity, Lh is no longer positive semi-definite, need to truncate negative eigenvalues
  DO x1 = 1, nlongs
    IF (L_alpha % HorizEV2(x1,lev) < 0) THEN
      L_alpha % HorizEV2(x1,lev) = 0
    END IF
  END DO
  ! Condition for zero values
  !increment = SUM(L_alpha % HorizEV2(1:nlongs,lev)) * ConditionFudge
  !IF (increment < small) increment = small
  !  L_alpha % HorizEV2(1:nlongs,lev) = L_alpha % HorizEV2(1:nlongs,lev) + increment
  L_alpha % HorizEV2(1:nlongs,lev) = SQRT(L_alpha % HorizEV2(1:nlongs,lev))

  CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
             'U',                                  & ! Upper triangular of matrix specified
              nlongs,                                & ! Order of the matrix to be diagonalised
              L_alpha % HorizMode3(1:nlongs, 1:nlongs, lev), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
              nlongs,                                & ! Leading dimension of the matrix
              L_alpha % HorizEV3(1:nlongs,lev),             & ! Eigenvalue array
              Work(1:3*nlongs),                      & ! Work array
              3*nlongs,                              & ! Size of work array
              ierr)

  ! As hScale_alpha tends to infinity, Lh is no longer positive semi-definite, need to truncate negative eigenvalues
  DO x1 = 1, nlongs
    IF (L_alpha % HorizEV3(x1,lev) < 0) THEN
      L_alpha % HorizEV3(x1,lev) = 0
    END IF
  END DO
  ! Condition for zero values
  !increment = SUM(L_alpha % HorizEV3(1:nlongs,lev)) * ConditionFudge
  !IF (increment < small) increment = small
  !  L_alpha % HorizEV3(1:nlongs,lev) = L_alpha % HorizEV3(1:nlongs,lev) + increment
  L_alpha % HorizEV3(1:nlongs,lev) = SQRT(L_alpha % HorizEV3(1:nlongs,lev))

  CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
             'U',                                  & ! Upper triangular of matrix specified
              nlongs,                                & ! Order of the matrix to be diagonalised
              L_alpha % HorizMode4(1:nlongs, 1:nlongs, lev), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
              nlongs,                                & ! Leading dimension of the matrix
              L_alpha % HorizEV4(1:nlongs,lev),             & ! Eigenvalue array
              Work(1:3*nlongs),                      & ! Work array
              3*nlongs,                              & ! Size of work array
              ierr)

  ! As hScale_alpha tends to infinity, Lh is no longer positive semi-definite, need to truncate negative eigenvalues
  DO x1 = 1, nlongs
    IF (L_alpha % HorizEV4(x1,lev) < 0) THEN
      L_alpha % HorizEV4(x1,lev) = 0
    END IF
  END DO
  ! Condition for zero values
  !increment = SUM(L_alpha % HorizEV4(1:nlongs,lev)) * ConditionFudge
  !IF (increment < small) increment = small
  !  L_alpha % HorizEV4(1:nlongs,lev) = L_alpha % HorizEV4(1:nlongs,lev) + increment
  L_alpha % HorizEV4(1:nlongs,lev) = SQRT(L_alpha % HorizEV4(1:nlongs,lev))

  CALL DSYEV('V',                                  & ! Eigenvalues and vectors to be computed
             'U',                                  & ! Upper triangular of matrix specified
              nlongs,                                & ! Order of the matrix to be diagonalised
              L_alpha % HorizMode5(1:nlongs, 1:nlongs, lev), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
              nlongs,                                & ! Leading dimension of the matrix
              L_alpha % HorizEV5(1:nlongs,lev),             & ! Eigenvalue array
              Work(1:3*nlongs),                      & ! Work array
              3*nlongs,                              & ! Size of work array
              ierr)

  ! As hScale_alpha tends to infinity, Lh is no longer positive semi-definite, need to truncate negative eigenvalues
  DO x1 = 1, nlongs
    IF (L_alpha % HorizEV5(x1,lev) < 0) THEN
      L_alpha % HorizEV5(x1,lev) = 0
    END IF
  END DO
  ! Condition for zero values
  !increment = SUM(L_alpha % HorizEV5(1:nlongs,lev)) * ConditionFudge
  !IF (increment < small) increment = small
  !  L_alpha % HorizEV5(1:nlongs,lev) = L_alpha % HorizEV5(1:nlongs,lev) + increment
  L_alpha % HorizEV5(1:nlongs,lev) = SQRT(L_alpha % HorizEV5(1:nlongs,lev))

END DO

DEALLOCATE (Func)
DEALLOCATE (Interim)
DEALLOCATE (Lh)
DEALLOCATE (Work)

END SUBROUTINE SetHorizLoc_alpha
