PROGRAM Util_RawCor

!*****************************************************
!*   Code to compute a selection of raw              *
!*   background error correlations.                  *
!*                                                   *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    nlongs, nlevs,               &
    dims_type,                   &
    ABC_type,                    &
    datadirRawCov,               &
    datadirABCperts,             &
    ImplCov_npoints,             &
    longindex,                   &
    levindex,                    &
    Npointsmax,                  &
    NEnsmax, NEns,               &
    NEnsMems,                    &
    NNMCmax, NNMC,               &
    Nlats


IMPLICIT NONE

! Declare variables
!==========================
TYPE(dims_type)          :: dims
TYPE(ABC_type)           :: ABC_member, Cov, StdDev
REAL(ZREAL8)             :: Neffmems_inv, source, sourcestddev
INTEGER                  :: x_source, z_source
INTEGER                  :: point, var, ens, mem, lat, item
CHARACTER(LEN=6)         :: varname(1:6)
CHARACTER(LEN=320)       :: ABCfile, RawCov_filename


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_RawCov'
PRINT*, '*************************************************************************'

! Read namelist
CALL SetOptions

CALL Initialise_dims(dims)
CALL Initialise_model_vars(ABC_member, .FALSE.)
CALL Initialise_model_vars(Cov, .FALSE.)
CALL Initialise_model_vars(StdDev, .FALSE.)

IF ((ImplCov_npoints < 1) .OR. (ImplCov_npoints > Npointsmax)) THEN
  PRINT *, 'ImplCov_npoints must be between 1 and ', Npointsmax
  PRINT *, 'Either run with different ImplCov_npoints, or recompile changing Npointsmax'
  STOP
END IF

varname(1) = 'u'
varname(2) = 'v'
varname(3) = 'w'
varname(4) = 'r'
varname(5) = 'b'
varname(6) = 'tracer'


Neffmems_inv = 1.0 / REAL(NEns * NEnsMems * Nlats)

IF (ImplCov_npoints > 0) THEN
  PRINT *, '===== Running some point-by-point covariances ====='


  ! Compute standard deviations
  CALL Initialise_model_vars (StdDev, .FALSE.)
  IF (NEns > 0) THEN
    DO ens = 1, NEns
      item = 0
      DO lat = 1, Nlats
        DO mem = 1, NEnsMems
          item = item + 1
          ! Read-in this ensemble member
          WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
          PRINT *, 'Reading file ', TRIM(ABCfile)
          CALL Read_state_2d (ABCfile, ABC_member, dims, 1, .TRUE.)

          IF ((point ==1).AND.(var == 1).AND.(ens == 1).AND.(mem == 1)) THEN
            ! Set-up model grid with info in file read-in
            CALL Set_ht_dep_cons (dims)
          END IF

          StdDev % u(1:nlongs,1:nlevs) = StdDev % u(1:nlongs,1:nlevs) + &
                                         ABC_member % u(1:nlongs,1:nlevs) * ABC_member % u(1:nlongs,1:nlevs)
          StdDev % v(1:nlongs,1:nlevs) = StdDev % v(1:nlongs,1:nlevs) + &
                                         ABC_member % v(1:nlongs,1:nlevs) * ABC_member % v(1:nlongs,1:nlevs)
          StdDev % w(1:nlongs,1:nlevs) = StdDev % w(1:nlongs,1:nlevs) + &
                                         ABC_member % w(1:nlongs,1:nlevs) * ABC_member % w(1:nlongs,1:nlevs)
          StdDev % r(1:nlongs,1:nlevs) = StdDev % r(1:nlongs,1:nlevs) + &
                                         ABC_member % r(1:nlongs,1:nlevs) * ABC_member % r(1:nlongs,1:nlevs)
          StdDev % b(1:nlongs,1:nlevs) = StdDev % b(1:nlongs,1:nlevs) + &
                                         ABC_member % b(1:nlongs,1:nlevs) * ABC_member % b(1:nlongs,1:nlevs)
        END DO
      END DO
    END DO

    ! Normalize
    StdDev % u(1:nlongs,1:nlevs) = SQRT(StdDev % u(1:nlongs,1:nlevs) * Neffmems_inv)
    StdDev % v(1:nlongs,1:nlevs) = SQRT(StdDev % v(1:nlongs,1:nlevs) * Neffmems_inv)
    StdDev % w(1:nlongs,1:nlevs) = SQRT(StdDev % w(1:nlongs,1:nlevs) * Neffmems_inv)
    StdDev % r(1:nlongs,1:nlevs) = SQRT(StdDev % r(1:nlongs,1:nlevs) * Neffmems_inv)
    StdDev % b(1:nlongs,1:nlevs) = SQRT(StdDev % b(1:nlongs,1:nlevs) * Neffmems_inv)
  END IF


  ! Loop over source points
  DO point = 1, ImplCov_npoints

    x_source = longindex(point)
    z_source = levindex(point)

    ! Repeat for each source field
    DO var = 1, 6
      CALL Initialise_model_vars (Cov, .FALSE.)
      IF (NEns > 0) THEN
        DO ens = 1, NEns
          item = 0
          DO lat = 1, Nlats
            DO mem = 1, NEnsMems
              item = item + 1
              ! Read-in this ensemble member
              WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
              PRINT *, 'Reading file ', TRIM(ABCfile)
              CALL Read_state_2d (ABCfile, ABC_member, dims, 1, .TRUE.)

              IF ((point ==1).AND.(var == 1).AND.(ens == 1).AND.(mem == 1)) THEN
                ! Set-up model grid with info in file read-in
                CALL Set_ht_dep_cons (dims)
              END IF

              ! Source point
              SELECT CASE (var)
              CASE (1)
                source = ABC_member % u(x_source,z_source)
              CASE (2)
                source = ABC_member % v(x_source,z_source)
              CASE (3)
                source = ABC_member % w(x_source,z_source)
              CASE (4)
                source = ABC_member % r(x_source,z_source)
              CASE (5)
                source = ABC_member % b(x_source,z_source)
              CASE (6)
                source = ABC_member % tracer(x_source,z_source)
              END SELECT
              Cov % u(1:nlongs,1:nlevs)      = Cov % u(1:nlongs,1:nlevs)      + source * ABC_member % u(1:nlongs,1:nlevs)
              Cov % v(1:nlongs,1:nlevs)      = Cov % v(1:nlongs,1:nlevs)      + source * ABC_member % v(1:nlongs,1:nlevs)
              Cov % w(1:nlongs,1:nlevs)      = Cov % w(1:nlongs,1:nlevs)      + source * ABC_member % w(1:nlongs,1:nlevs)
              Cov % r(1:nlongs,1:nlevs)      = Cov % r(1:nlongs,1:nlevs)      + source * ABC_member % r(1:nlongs,1:nlevs)
              Cov % b(1:nlongs,1:nlevs)      = Cov % b(1:nlongs,1:nlevs)      + source * ABC_member % b(1:nlongs,1:nlevs)
              Cov % tracer(1:nlongs,1:nlevs) = Cov % tracer(1:nlongs,1:nlevs) + source * ABC_member % tracer(1:nlongs,1:nlevs)
            END DO
          END DO
        END DO

        ! Normalize
        Cov % u(1:nlongs,1:nlevs)      = Cov % u(1:nlongs,1:nlevs)      * Neffmems_inv
        Cov % v(1:nlongs,1:nlevs)      = Cov % v(1:nlongs,1:nlevs)      * Neffmems_inv
        Cov % w(1:nlongs,1:nlevs)      = Cov % w(1:nlongs,1:nlevs)      * Neffmems_inv
        Cov % r(1:nlongs,1:nlevs)      = Cov % r(1:nlongs,1:nlevs)      * Neffmems_inv
        Cov % b(1:nlongs,1:nlevs)      = Cov % b(1:nlongs,1:nlevs)      * Neffmems_inv
        Cov % tracer(1:nlongs,1:nlevs) = Cov % tracer(1:nlongs,1:nlevs) * Neffmems_inv

        ! Convert covariances into correlations
        ! Source point standard deviation
        SELECT CASE (var)
        CASE (1)
          sourcestddev = StdDev % u(x_source,z_source)
        CASE (2)
          sourcestddev = StdDev % v(x_source,z_source)
        CASE (3)
          sourcestddev = StdDev % w(x_source,z_source)
        CASE (4)
          sourcestddev = StdDev % r(x_source,z_source)
        CASE (5)
          sourcestddev = StdDev % b(x_source,z_source)
        CASE (6)
          sourcestddev = StdDev % tracer(x_source,z_source)
        END SELECT

        Cov % u(1:nlongs,1:nlevs) = Cov % u(1:nlongs,1:nlevs) / (StdDev % u(1:nlongs,1:nlevs) * sourcestddev)
        Cov % v(1:nlongs,1:nlevs) = Cov % v(1:nlongs,1:nlevs) / (StdDev % v(1:nlongs,1:nlevs) * sourcestddev)
        Cov % w(1:nlongs,1:nlevs) = Cov % w(1:nlongs,1:nlevs) / (StdDev % w(1:nlongs,1:nlevs) * sourcestddev)
        Cov % r(1:nlongs,1:nlevs) = Cov % r(1:nlongs,1:nlevs) / (StdDev % r(1:nlongs,1:nlevs) * sourcestddev)
        Cov % b(1:nlongs,1:nlevs) = Cov % b(1:nlongs,1:nlevs) / (StdDev % b(1:nlongs,1:nlevs) * sourcestddev)

        ! Output these covariance fields
        WRITE (RawCov_filename, '(A,A,I0.3,A,A,A)') TRIM(datadirRawCov), '/Point_', point, &
               '_delta', TRIM(varname(var)), '.nc'
        CALL Write_state_2d (RawCov_filename, Cov, dims, 1, 0, 0, .TRUE.)


      END IF
    END DO

  END DO
END IF

! Tidy up
CALL Deallocate_dims(dims)
CALL Deallocate_model_vars(ABC_member)
CALL Deallocate_model_vars(Cov)
CALL Deallocate_model_vars(StdDev)

END PROGRAM Util_RawCor
