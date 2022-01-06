SUBROUTINE ConvScaleDens (u, v, w, rpls, rpcs, rhols, bpcs,       &   ! Reference state
                          du, dv, dw, drpls, dbpcs,               &   ! Input increments
                          drpcs,                                  &   ! Output increment
                          dims )                                      ! Meta data

! Code to find the convective-scale density perturbation.
! Involves solving an elliptic equation.
! Ross Bannister, National Centre for Earth Observation, May 2020.

! This code should be viewed in combination with the documentation
! (AnelasticBal.pdf)

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  nlongs,                 &
  nlevs,                  &
  dims_type,              &
  f, B, C, half, recipdx, &
  recip2dx, pi, dx

IMPLICIT NONE

! Variable naming convection:
! ---------------------------
! No 'd' prefix means reference state variable
! 'p' is for primed variable
! 'ls' suffix is for large scale
! 'cs' suffix is for convective scale
! 'zm' suffix is for zonal mean

! Delcare input reference state variables
REAL(ZREAL8),   INTENT(IN)    :: u(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),   INTENT(IN)    :: v(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),   INTENT(IN)    :: w(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),   INTENT(IN)    :: rpls(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),   INTENT(IN)    :: rpcs(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),   INTENT(IN)    :: rhols(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),   INTENT(IN)    :: bpcs(0:nlongs+1, 0:nlevs+1)

! Declare input pert variables
REAL(ZREAL8),   INTENT(IN)    :: du(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),   INTENT(IN)    :: dv(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),   INTENT(IN)    :: dw(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),   INTENT(IN)    :: drpls(0:nlongs+1, 0:nlevs+1)
REAL(ZREAL8),   INTENT(IN)    :: dbpcs(0:nlongs+1, 0:nlevs+1)

! Declare the output pert variable
REAL(ZREAL8),   INTENT(OUT)   :: drpcs(0:nlongs+1, 0:nlevs+1)

! Declare other parameters
TYPE(dims_type),INTENT(IN)    :: dims

! Declare local variables
INTEGER                       :: x, z, zp, k
REAL(ZREAL8), ALLOCATABLE     :: level(:), column(:)
REAL(ZREAL8), ALLOCATABLE     :: field(:,:)
REAL(ZREAL8), ALLOCATABLE     :: drpls_fl(:,:)
REAL(ZREAL8), ALLOCATABLE     :: rpls_fl(:, :)
REAL(ZREAL8), ALLOCATABLE     :: rhols_fl(:,:)
REAL(ZREAL8), ALLOCATABLE     :: rhols_fl_zm(:)
REAL(ZREAL8), ALLOCATABLE     :: rhols_zm(:)
REAL(ZREAL8)                  :: BdivC, fdivC, specfac
REAL(ZREAL8), ALLOCATABLE     :: Term1a(:,:), Term1b(:,:)
REAL(ZREAL8), ALLOCATABLE     :: Term2a(:,:), Term2b(:,:)
REAL(ZREAL8), ALLOCATABLE     :: Term3a(:,:), Term3b(:,:)
REAL(ZREAL8), ALLOCATABLE     :: Term4a(:,:), Term4b(:,:)
REAL(ZREAL8), ALLOCATABLE     :: Term5a(:,:), Term5b(:,:)
REAL(ZREAL8), ALLOCATABLE     :: Term6a(:,:), Term6b(:,:)
REAL(ZREAL8), ALLOCATABLE     :: Term7a(:,:), Term7b(:,:)
REAL(ZREAL8), ALLOCATABLE     :: Term8a(:,:), Term8b(:,:)
REAL(ZREAL8), ALLOCATABLE     :: Term9a(:,:), Term9b(:,:)
REAL(ZREAL8), ALLOCATABLE     :: Term10a(:,:), Term10b(:,:)
REAL(ZREAL8), ALLOCATABLE     :: RHSr(:,:), RHSs(:,:)
REAL(ZREAL8), ALLOCATABLE     :: P(:,:), Pinv(:,:,:)
!REAL(ZREAL8), ALLOCATABLE     :: Check(:)
REAL(ZREAL8), ALLOCATABLE     :: y(:)
REAL(ZREAL8), ALLOCATABLE     :: Solutions(:,:)
INTEGER                       :: real_index, imag_index

! Declare functions
REAL(ZREAL8)                  :: INT_HF, INT_FH




BdivC   = B / C
fdivC   = f / C
specfac = -4.0 * pi*pi / (REAL(nlongs)*REAL(nlongs)*dx*dx)


PRINT *, 'Computing a selection of data on full levels and zonal means'


! Calculate drpls on full levels
! ------------------------------
ALLOCATE (drpls_fl(0:nlongs+1, 0:nlevs+1))
DO x = 1, nlongs
  DO z = 0, nlevs
    drpls_fl(x,z) = INT_HF(drpls(x,z), drpls(x,z+1), z, dims)
  END DO
  drpls_fl(x,nlevs+1) = drpls_fl(x,nlevs)
END DO

! Calculate rpls on full levels
! -----------------------------
ALLOCATE (rpls_fl(0:nlongs+1, 0:nlevs+1))
DO x = 1, nlongs
  DO z = 0, nlevs
    rpls_fl(x,z) = INT_HF(rpls(x,z), rpls(x,z+1), z, dims)
  END DO
  rpls_fl(x,nlevs+1) = rpls_fl(x,nlevs)
END DO

! Calculate rhols on full levels
! ------------------------------
ALLOCATE (rhols_fl(0:nlongs+1, 0:nlevs+1))
DO x = 1, nlongs
  DO z = 0, nlevs
    rhols_fl(x,z) = INT_HF(rhols(x,z), rhols(x,z+1), z, dims)
  END DO
  rhols_fl(x,nlevs+1) = rhols_fl(x,nlevs)
END DO

! Calculate the zonal mean of rhols on full levels
! ------------------------------------------------
ALLOCATE (rhols_fl_zm(0:nlevs+1))
DO z = 0, nlevs+1
  rhols_fl_zm(z) = SUM(rhols_fl(1:nlongs,z)) / REAL(nlongs)
END DO

! Calculate the zonal mean of rhols
! ---------------------------------
ALLOCATE (rhols_zm(0:nlevs+1))
DO z = 0, nlevs+1
  rhols_zm(z) = SUM(rhols(1:nlongs,z)) / REAL(nlongs)
END DO



ALLOCATE (level(0:nlongs+1))
ALLOCATE (column(0:nlevs+1))
ALLOCATE (field(0:nlongs+1, 0:nlevs+1))

! Term 1a
! -------
ALLOCATE (Term1a(1:nlongs, 1:nlevs))
PRINT *, 'Term 1a'
DO z = 1, nlevs
  DO x = 1, nlongs
    level(x) = (rpcs(x+1,z) - rpcs(x,z)) * recipdx * (drpls(x,z) + drpls(x+1,z)) * half
  END DO
  level(0)        = level(nlongs)
  level(nlongs+1) = level(1)
  DO x = 1, nlongs
    Term1a(x,z) = -1.0 * (level(x) - level(x-1)) * recipdx
  END DO
END DO

! Term 1b
! -------
ALLOCATE (Term1b(1:nlongs, 1:nlevs))
PRINT *, 'Term 1b'
DO x = 1, nlongs
  DO z = 1, nlevs
    column(z) = (rpcs(x,z+1) - rpcs(x,z)) / (dims % half_levs(z+1) - dims % half_levs(z)) * &
                drpls_fl(x,z)
  END DO
  column(0)       = column(1)
  column(nlevs+1) = column(nlevs)
  DO z = 1, nlevs
    Term1b(x,z)   = -1.0 * (column(z) - column(z-1)) / (dims % full_levs(z) - dims % full_levs(z-1))
  END DO
END DO

! Term 2a
! -------
ALLOCATE (Term2a(1:nlongs, 1:nlevs))
PRINT *, 'Term 2a'
DO z = 1, nlevs
  DO x = 1, nlongs
    level(x) = (u(x,z) - u(x-1,z)) * recipdx * (u(x,z) + u(x-1,z)) * half * drpls(x,z)
  END DO
  level(0)        = level(nlongs)
  level(nlongs+1) = level(1)
  DO x = 1, nlongs
    Term2a(x,z)     = -1.0 * BdivC * (level(x+1) - level(x-1)) * recip2dx
  END DO
END DO

! Term 2b
! -------
ALLOCATE (Term2b(1:nlongs, 1:nlevs))
PRINT *, 'Term 2b'
DO z = 0, nlevs
  DO x = 1, nlongs
    level(x) = (u(x,z+1) - u(x,z)) / (dims % half_levs(z+1) - dims % half_levs(z)) * &
               (w(x,z) + w(x+1,z)) * half * (drpls_fl(x,z) + drpls_fl(x+1,z)) * half
  END DO
  level(0)        = level(nlongs)
  level(nlongs+1) = level(1)
  DO x = 1, nlongs
    field(x,z) = (level(x) - level(x-1)) * recipdx
  END DO
END DO
DO x = 1, nlongs
  field(x,nlevs+1) = field(x,nlevs)
END DO

DO z = 1, nlevs
  DO x = 1, nlongs
    Term2b(x,z) = -1.0 * BdivC * INT_FH(field(x,z-1), field(x,z), z, dims)
  END DO
END DO

! Term 3a (this is just an alteration of 2a)
! ------------------------------------------
ALLOCATE (Term3a(1:nlongs, 1:nlevs))
PRINT *, 'Term 3a'
DO z = 1, nlevs
  DO x = 1, nlongs
    level(x) = (u(x,z) - u(x-1,z)) * recipdx * (du(x,z) + du(x-1,z)) * half * rpls(x,z)
  END DO
  level(0)        = level(nlongs)
  level(nlongs+1) = level(1)
  DO x = 1, nlongs
    Term3a(x,z)     = -1.0 * BdivC * (level(x+1) - level(x-1)) * recip2dx
  END DO
END DO

! Term 3b (this is just an alteration of 2b)
! ------------------------------------------
ALLOCATE (Term3b(1:nlongs, 1:nlevs))
PRINT *, 'Term 3b'
DO z = 0, nlevs
  DO x = 1, nlongs
    level(x) = (u(x,z+1) - u(x,z)) / (dims % half_levs(z+1) - dims % half_levs(z)) * &
               (dw(x,z) + dw(x+1,z)) * half * (rpls_fl(x,z) + rpls_fl(x+1,z)) * half
  END DO
  level(0)        = level(nlongs)
  level(nlongs+1) = level(1)
  DO x = 1, nlongs
    field(x,z) = (level(x) - level(x-1)) * recipdx
  END DO
END DO
DO x = 1, nlongs
  field(x,nlevs+1) = field(x,nlevs)
END DO

DO z = 1, nlevs
  DO x = 1, nlongs
    Term3b(x,z) = -1.0 * BdivC * INT_FH (field(x,z-1), field(x,z), z, dims)
  END DO
END DO

! Term 4a (this is just an alteration of 2a)
! ------------------------------------------
ALLOCATE (Term4a(1:nlongs, 1:nlevs))
PRINT *, 'Term 4a'
DO z = 1, nlevs
  DO x = 1, nlongs
    level(x) = (du(x,z) - du(x-1,z)) * recipdx * (u(x,z) + u(x-1,z)) * half * rpls(x,z)
  END DO
  level(0)        = level(nlongs)
  level(nlongs+1) = level(1)
  DO x = 1, nlongs
    Term4a(x,z)     = -1.0 * BdivC * (level(x+1) - level(x-1)) * recip2dx
  END DO
END DO

! Term 4b (this is just an alteration of 2b)
! ------------------------------------------
ALLOCATE (Term4b(1:nlongs, 1:nlevs))
PRINT *, 'Term 4b'
DO z = 0, nlevs
  DO x = 1, nlongs
    level(x) = (du(x,z+1) - du(x,z)) / (dims % half_levs(z+1) - dims % half_levs(z)) * &
               (w(x,z) + w(x+1,z)) * half * (rpls_fl(x,z) + rpls_fl(x+1,z)) * half
  END DO
  level(0)        = level(nlongs)
  level(nlongs+1) = level(1)
  DO x = 1, nlongs
    field(x,z) = (level(x) - level(x-1)) * recipdx
  END DO
END DO
DO x = 1, nlongs
  field(x,nlevs+1) = field(x,nlevs)
END DO

DO z = 1, nlevs
  DO x = 1, nlongs
    Term4b(x,z) = -1.0 * BdivC * INT_FH(field(x,z-1), field(x,z), z, dims)
  END DO
END DO

! Term 5a
! -------
ALLOCATE (Term5a(1:nlongs, 1:nlevs))
PRINT *, 'Term 5a'
DO x = 1, nlongs
  DO z = 0, nlevs
    column(z) = (w(x+1,z) - w(x,z)) * recipdx * INT_HF (u(x,z), u(x,z+1), z, dims) *    &
                (drpls_fl(x,z) + drpls_fl(x+1,z)) * half
  END DO
  column(nlevs+1) = column(z)

  DO z = 1, nlevs+1
    field(x,z) = (column(z) - column(z-1)) / (dims % full_levs(z) - dims % full_levs(z-1))
  END DO
END DO
DO z = 1, nlevs
  field(0,z)       = field(nlevs,z)
  field(nlevs+1,z) = field(1,z)
END DO
DO x = 1, nlongs
  DO z = 1, nlevs
    Term5a(x,z) = -1.0 * BdivC * (field(x-1,z) + field(x,z)) * half
  END DO
END DO

! Term 5b
! -------
ALLOCATE (Term5b(1:nlongs, 1:nlevs))
PRINT *, 'Term 5b'
DO x = 1, nlongs
  DO z = 1, nlevs+1
    column(z) = (w(x,z) - w(x,z-1)) / (dims % full_levs(z) - dims % full_levs(z-1)) *    &
                INT_FH (w(x,z-1), w(x,z), z, dims) * drpls(x,z)
  END DO
  column(0) = column(1)
  DO z = 1, nlevs
    Term5b(x,z) = -1.0 * BdivC * (column(z+1) - column(z-1)) / (dims % half_levs(z+1) - dims % half_levs(z-1))
  END DO
END DO

! Term 6a (this is just an alteration of 5a)
! ------------------------------------------
ALLOCATE (Term6a(1:nlongs, 1:nlevs))
PRINT *, 'Term 6a'
DO x = 1, nlongs
  DO z = 0, nlevs
    column(z) = (w(x+1,z) - w(x,z)) * recipdx * INT_HF (du(x,z), du(x,z+1), z, dims) *    &
                (rpls_fl(x,z) + rpls_fl(x+1,z)) * half
  END DO
  column(nlevs+1) = column(z)

  DO z = 1, nlevs+1
    field(x,z) = (column(z) - column(z-1)) / (dims % full_levs(z) - dims % full_levs(z-1))
  END DO
END DO
DO z = 1, nlevs
  field(0,z)       = field(nlevs,z)
  field(nlevs+1,z) = field(1,z)
END DO
DO x = 1, nlongs
  DO z = 1, nlevs
    Term6a(x,z) = -1.0 * BdivC * (field(x-1,z) + field(x,z)) * half
  END DO
END DO

! Term 6b (this is just an alteration of 5b)
! ------------------------------------------
ALLOCATE (Term6b(1:nlongs, 1:nlevs))
PRINT *, 'Term 6b'
DO x = 1, nlongs
  DO z = 1, nlevs+1
    column(z) = (w(x,z) - w(x,z-1)) / (dims % full_levs(z) - dims % full_levs(z-1)) *    &
                INT_FH (dw(x,z-1), dw(x,z), z, dims) * rpls(x,z)
  END DO
  column(0) = column(1)
  DO z = 1, nlevs
    Term6b(x,z) = -1.0 * BdivC * (column(z+1) - column(z-1)) / (dims % half_levs(z+1) - dims % half_levs(z-1))
  END DO
END DO

! Term 7a (this is just an alteration of 5a)
! ------------------------------------------
ALLOCATE (Term7a(1:nlongs, 1:nlevs))
PRINT *, 'Term 7a'
DO x = 1, nlongs
  DO z = 0, nlevs
    column(z) = (dw(x+1,z) - dw(x,z)) * recipdx * INT_HF (u(x,z), u(x,z+1), z, dims) *    &
                (rpls_fl(x,z) + rpls_fl(x+1,z)) * half
  END DO
  column(nlevs+1) = column(z)

  DO z = 1, nlevs+1
    field(x,z) = (column(z) - column(z-1)) / (dims % full_levs(z) - dims % full_levs(z-1))
  END DO
END DO
DO z = 1, nlevs
  field(0,z)       = field(nlevs,z)
  field(nlevs+1,z) = field(1,z)
END DO
DO x = 1, nlongs
  DO z = 1, nlevs
    Term7a(x,z) = -1.0 * BdivC * (field(x-1,z) + field(x,z)) * half
  END DO
END DO

! Term 7b (this is just an alteration of 5b)
! ------------------------------------------
ALLOCATE (Term7b(1:nlongs, 1:nlevs))
PRINT *, 'Term 7b'
DO x = 1, nlongs
  DO z = 1, nlevs+1
    column(z) = (dw(x,z) - dw(x,z-1)) / (dims % full_levs(z) - dims % full_levs(z-1)) *    &
                INT_FH (w(x,z-1), w(x,z), z, dims) * rpls(x,z)
  END DO
  column(0) = column(1)
  DO z = 1, nlevs
    Term7b(x,z) = -1.0 * BdivC * (column(z+1) - column(z-1)) / (dims % half_levs(z+1) - dims % half_levs(z-1))
  END DO
END DO

! Term 8a
! -------
ALLOCATE (Term8a(1:nlongs, 1:nlevs))
PRINT *, 'Term 8a'
DO z = 1, nlevs
  DO x = 1, nlongs
    level(x) = (rpls(x+1,z) - rpls(x,z)) * recipdx * (drpls(x,z) + drpls(x+1,z)) * half
  END DO
  level(0)        = level(nlongs)
  level(nlongs+1) = level(1)
  DO x = 1, nlongs
    Term8a(x,z)     = -1.0 * (level(x) - level(x-1)) * recipdx
  END DO
END DO

! Term 8b (this is just an alteration of 8a)
! ------------------------------------------
ALLOCATE (Term8b(1:nlongs, 1:nlevs))
PRINT *, 'Term 8b'
DO z = 1, nlevs
  DO x = 1, nlongs
    level(x) = (drpls(x+1,z) - drpls(x,z)) * recipdx * (rpls(x,z) + rpls(x+1,z)) * half
  END DO
  level(0)        = level(nlongs)
  level(nlongs+1) = level(1)
  DO x = 1, nlongs
    Term8b(x,z)     = -1.0 * (level(x) - level(x-1)) * recipdx
  END DO
END DO

! Term 9a
! -------
ALLOCATE (Term9a(1:nlongs, 1:nlevs))
PRINT *, 'Term 9a'
DO z = 1, nlevs
  DO x = 1, nlongs
    Term9a(x,z) = fdivC * (v(x+1,z) * drpls(x+1,z) - v(x-1,z) * drpls(x-1,z)) * recip2dx
  END DO
END DO

! Term 9b (this is just an alteration of 9a)
! ------------------------------------------
ALLOCATE (Term9b(1:nlongs, 1:nlevs))
PRINT *, 'Term 9b'
DO z = 1, nlevs
  DO x = 1, nlongs
    Term9b(x,z) = fdivC * (dv(x+1,z) * rpls(x+1,z) - dv(x-1,z) * rpls(x-1,z)) * recip2dx
  END DO
END DO

! Term 10a
! --------
ALLOCATE (Term10a(1:nlongs, 1:nlevs))
PRINT *, 'Term 10a'
DO z = 1, nlevs
  DO x = 1, nlongs
    Term10a(x,z) = (drpls_fl(x,z) * bpcs(x,z) - drpls_fl(x,z-1) * bpcs(x,z-1)) /&
                   (dims % full_levs(z) - dims % full_levs(z-1)) / C
  END DO
END DO

! Term 10b (this is just an alteration of 10a)
! --------------------------------------------
ALLOCATE (Term10b(1:nlongs, 1:nlevs))
PRINT *, 'Term 10b'
DO z = 1, nlevs
  DO x = 1, nlongs
    Term10b(x,z) = (rpls_fl(x,z) * dbpcs(x,z) - rpls_fl(x,z-1) * dbpcs(x,z-1)) /&
                   (dims % full_levs(z) - dims % full_levs(z-1)) / C
  END DO
END DO

DEALLOCATE (level, column, field)



! Find the total RHS field
! ------------------------
ALLOCATE (RHSr(1:nlongs, 1:nlevs))
PRINT *, 'Calculating total RHS'
RHSr(1:nlongs, 1:nlevs) = Term1a (1:nlongs, 1:nlevs) + Term1b (1:nlongs, 1:nlevs) + &
                          Term2a (1:nlongs, 1:nlevs) + Term2b (1:nlongs, 1:nlevs) + &
                          Term3a (1:nlongs, 1:nlevs) + Term3b (1:nlongs, 1:nlevs) + &
                          Term4a (1:nlongs, 1:nlevs) + Term4b (1:nlongs, 1:nlevs) + &
                          Term5a (1:nlongs, 1:nlevs) + Term5b (1:nlongs, 1:nlevs) + &
                          Term6a (1:nlongs, 1:nlevs) + Term6b (1:nlongs, 1:nlevs) + &
                          Term7a (1:nlongs, 1:nlevs) + Term7b (1:nlongs, 1:nlevs) + &
                          Term8a (1:nlongs, 1:nlevs) + Term8b (1:nlongs, 1:nlevs) + &
                          Term9a (1:nlongs, 1:nlevs) + Term9b (1:nlongs, 1:nlevs) + &
                          Term10a(1:nlongs, 1:nlevs) + Term10b(1:nlongs, 1:nlevs)

CALL Write_one_field ('Term1a.nc', nlongs, nlevs, Term1a(1:nlongs,1:nlevs), 'Term1a')
CALL Write_one_field ('Term1b.nc', nlongs, nlevs, Term1b(1:nlongs,1:nlevs), 'Term1b')
CALL Write_one_field ('Term2a.nc', nlongs, nlevs, Term2a(1:nlongs,1:nlevs), 'Term2a')
CALL Write_one_field ('Term2b.nc', nlongs, nlevs, Term2b(1:nlongs,1:nlevs), 'Term2b')
CALL Write_one_field ('Term3a.nc', nlongs, nlevs, Term3a(1:nlongs,1:nlevs), 'Term3a')
CALL Write_one_field ('Term3b.nc', nlongs, nlevs, Term3b(1:nlongs,1:nlevs), 'Term3b')
CALL Write_one_field ('Term4a.nc', nlongs, nlevs, Term4a(1:nlongs,1:nlevs), 'Term4a')
CALL Write_one_field ('Term4b.nc', nlongs, nlevs, Term4b(1:nlongs,1:nlevs), 'Term4b')
CALL Write_one_field ('Term5a.nc', nlongs, nlevs, Term5a(1:nlongs,1:nlevs), 'Term5a')
CALL Write_one_field ('Term5b.nc', nlongs, nlevs, Term5b(1:nlongs,1:nlevs), 'Term5b')
CALL Write_one_field ('Term6a.nc', nlongs, nlevs, Term6a(1:nlongs,1:nlevs), 'Term6a')
CALL Write_one_field ('Term6b.nc', nlongs, nlevs, Term6b(1:nlongs,1:nlevs), 'Term6b')
CALL Write_one_field ('Term7a.nc', nlongs, nlevs, Term7a(1:nlongs,1:nlevs), 'Term7a')
CALL Write_one_field ('Term7b.nc', nlongs, nlevs, Term7b(1:nlongs,1:nlevs), 'Term7b')
CALL Write_one_field ('Term8a.nc', nlongs, nlevs, Term8a(1:nlongs,1:nlevs), 'Term8a')
CALL Write_one_field ('Term8b.nc', nlongs, nlevs, Term8b(1:nlongs,1:nlevs), 'Term8b')
CALL Write_one_field ('Term9a.nc', nlongs, nlevs, Term9a(1:nlongs,1:nlevs), 'Term9a')
CALL Write_one_field ('Term9b.nc', nlongs, nlevs, Term9b(1:nlongs,1:nlevs), 'Term9b')
CALL Write_one_field ('Term10a.nc', nlongs, nlevs, Term10a(1:nlongs,1:nlevs), 'Term10a')
CALL Write_one_field ('Term10b.nc', nlongs, nlevs, Term10b(1:nlongs,1:nlevs), 'Term10b')
CALL Write_one_field ('RHS.nc', nlongs, nlevs, RHSr(1:nlongs,1:nlevs), 'RHS')

DEALLOCATE (Term1a, Term1b, Term2a, Term2b, Term3a, Term3b, Term4a, Term4b, Term5a, Term5b)
DEALLOCATE (Term6a, Term6b, Term7a, Term7b, Term8a, Term8b, Term9a, Term9b, Term10a, Term10b)
DEALLOCATE (drpls_fl, rpls_fl, rhols_fl)


! Fourier transform the RHS field
! -------------------------------
ALLOCATE (RHSs(1:nlongs, 1:nlevs))
PRINT *, 'Converting the total RHS to spectral space'
CALL fft_real2spec ( RHSr(1:nlongs,1:nlevs), RHSs(1:nlongs,1:nlevs) )
DEALLOCATE (RHSr)



! Construct the vertical operators and their inverses
! ---------------------------------------------------
ALLOCATE (P(1:nlevs, 1:nlevs))
ALLOCATE (Pinv(1:nlevs, 1:nlevs, 0:nlongs/2))
ALLOCATE (y(1:nlevs))
PRINT *, 'Constructing vertical operators and inverses'
! Loop over wavenumbers
DO k = 0, nlongs/2
  PRINT *, 'Wavenumber ', k
  P(1:nlevs, 1:nlevs) = 0.0
  ! Loop over each row of P
  DO z = 1, nlevs
    IF (z > 1) &
      P(z,z-1) = rhols_fl_zm(z-1) * dims % recip_full_k_km1(z) * dims % recip_half_k_km1(z)
    P(z,z)     = specfac * REAL(k*k) * rhols_zm(z) -                                         &
                 rhols_fl_zm(z) * dims % recip_full_k_km1(z) * dims % recip_half_kp1_k(z) -  &
                 rhols_fl_zm(z-1) * dims % recip_full_k_km1(z) * dims % recip_half_k_km1(z)
    IF (z < nlevs) &
      P(z,z+1) = rhols_fl_zm(z) * dims % recip_full_k_km1(z) * dims % recip_half_kp1_k(z)
  END DO

  ! Construct the inverse of each vertical operator
  ! Loop over each column of the inverse Pinv
  PRINT *, '  inverse'
  y(1:nlevs) = 0.0
  DO z = 1, nlevs
    y(z) = 1.0
    CALL GaussEl ( nlevs,              &
                   P(1:nlevs,1:nlevs), &
                   y(1:nlevs),         &
                   Pinv(1:nlevs,z,k) )
    y(z) = 0.0
  END DO


!  ! Output a selection of operators as a check
!  IF ((k == 0) .OR. (k == nlongs/2)) THEN
!    ALLOCATE (Check(1:nlevs))
!    PRINT *, 'Here is the operator P for wavenumber k = ', k
!    DO z = 1, nlevs
!      PRINT *, z, P(z,1:nlevs)
!    END DO
!    PRINT *, 'Here is the operator Pinv for wavenumber k = ', k
!    DO z = 1, nlevs
!      PRINT *, Pinv(z,1:nlevs,k)
!    END DO
!    PRINT *, 'Here is the result of P Pinv for wavenumber k = ', k
!    DO z = 1, nlevs
!      DO zp = 1, nlevs
!        Check(zp) = SUM(P(z,1:nlevs) * Pinv(1:nlevs,zp,k))
!      END DO
!      PRINT *, z, Check(1:nlevs)
!    END DO
!    PRINT *, 'Here is the result of Pinv P for wavenumber k = ', k
!    DO z = 1, nlevs
!      DO zp = 1, nlevs
!        Check(zp) = SUM(Pinv(z,1:nlevs,k) * P(1:nlevs,zp))
!      END DO
!      PRINT *, z, Check(1:nlevs)
!    END DO
!    DEALLOCATE (Check)
!  END IF

END DO
DEALLOCATE (P, y)
DEALLOCATE (rhols_fl_zm, rhols_zm)


!!!!!! Note: Pinv will eventually need to be stored in the CVT structure





! Operate with the inverse of the vertical operators to construct the solution in spectral space
! ----------------------------------------------------------------------------------------------
ALLOCATE (Solutions(1:nlongs, 1:nlevs))

PRINT *, 'Computing the solution in spectral space'
! The largest scale (spectral index 1, corresponding to k=0)
real_index = 1
DO z = 1, nlevs
  Solutions(real_index,z) = SUM(Pinv(z,1:nlevs,0) * RHSs(real_index,1:nlevs))
END DO

! The intermediate scales
DO k = 1, nlongs/2 - 1
  real_index = 2*k
  imag_index = 2*k + 1
  DO z = 1, nlevs
    Solutions(real_index,z) = SUM(Pinv(z,1:nlevs,k) * RHSs(real_index,1:nlevs))
    Solutions(imag_index,z) = SUM(Pinv(z,1:nlevs,k) * RHSs(imag_index,1:nlevs))
  END DO
END DO

! The smallest scale (spectral index nlongs, corresponding to k=nlongs/2)
real_index = nlongs
DO z = 1, nlevs
  Solutions(real_index,z) = SUM(Pinv(z,1:nlevs,nlongs/2) * RHSs(real_index,1:nlevs))
END DO
DEALLOCATE (RHSs, Pinv)

! Find the solution in real space
! -------------------------------

PRINT *, 'Converting the solution to real space'
CALL fft_spec2real(Solutions(1:nlongs, 1:nlevs), drpcs(1:nlongs, 1:nlevs))
DEALLOCATE (Solutions)


! Horizontal boundaries
PRINT *, 'Dealing with horizontal boundaries for the solution'
drpcs(0, 1:nlevs)          = drpcs(nlongs, 1:nlevs)
drpcs(nlongs+1, 1:nlevs)   = drpcs(1, 1:nlevs)
! Vertical boundaries
PRINT *, 'Dealing with vertical boundaries for the solution'
drpcs(0:nlongs+1,0)        = drpcs(0:nlongs+1,1)
drpcs(0:nlongs+1,nlevs+1)  = drpcs(0:nlongs+1,nlevs)
PRINT *, 'Finished dealing with boundaries'
PRINT *, 'End of subroutine'

END SUBROUTINE ConvScaleDens
