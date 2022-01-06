PROGRAM Util_ImpliedCor

!*****************************************************
!*   Code to compute a selection of implied          *
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
    CV_type,                     &
    CVT_type,                    &
    datadirImpliedCov,           &
    datadirABCfcs,               &
    LS_file,                     &
    datadirCVT,                  &
    CVT_file,                    &
    ImplCov_npoints,             &
    longindex,                   &
    levindex,                    &
    Npointsmax,                  &
    fft_wsave_x,                 &
    fft_work_x


IMPLICIT NONE

! Declare variables
!==========================
TYPE(dims_type)           :: dims
TYPE(ABC_type)            :: Model_data1, Model_data2, LS
TYPE(CV_type)             :: CV_data
TYPE(CVT_type)            :: CVT
INTEGER                   :: point, x, z
REAL(ZREAL8), ALLOCATABLE :: StdDevu(:)
REAL(ZREAL8), ALLOCATABLE :: StdDevv(:)
REAL(ZREAL8), ALLOCATABLE :: StdDevw(:)
REAL(ZREAL8), ALLOCATABLE :: StdDevr(:)
REAL(ZREAL8), ALLOCATABLE :: StdDevb(:)

CHARACTER(LEN=320)       :: LS_filename, CVT_filename, ImplCov_filename


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_ImpliedCov'
PRINT*, '*************************************************************************'

! Read namelist
CALL SetOptions

CALL Initialise_dims(dims)
CALL Initialise_model_vars(Model_data1, .FALSE.)
CALL Initialise_model_vars(Model_data2, .FALSE.)
CALL Initialise_model_vars(LS, .FALSE.)
CALL Initialise_CVs(CV_data, .FALSE.)
CALL Initialise_CVT(CVT)
ALLOCATE (StdDevu(1:nlevs))
ALLOCATE (StdDevv(1:nlevs))
ALLOCATE (StdDevw(1:nlevs))
ALLOCATE (StdDevr(1:nlevs))
ALLOCATE (StdDevb(1:nlevs))

IF ((ImplCov_npoints < 1) .OR. (ImplCov_npoints > Npointsmax)) THEN
  PRINT *, 'ImplCov_npoints must be between 1 and ', Npointsmax
  PRINT *, 'Either run with different ImplCov_npoints, or recompile changing Npointsmax'
  STOP
END IF


LS_filename          = TRIM(datadirABCfcs)   // '/' // TRIM(LS_file)
CVT_filename         = TRIM(datadirCVT)      // '/' // TRIM(CVT_file)

! Read in LS
PRINT*, 'Reading in linearisation state ...'
CALL Read_state_2d (LS_filename, LS, dims, 1, .TRUE.)
PRINT*, '-- done'

! Set some commonly-used constants
CALL Set_ht_dep_cons (dims)

! Read in CVT data prepared beforehand
PRINT*, 'Reading in CVT data ...'
CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
CALL Read_Covs (CVT_filename, CVT,              &
                .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.)
PRINT *, '-- done'



IF (ImplCov_npoints > 0) THEN


  ! Compute standard devs as a function of height (in order to allow computation of correlations)
  CALL Initialise_model_vars (Model_data1, .FALSE.)

  ! u field
  DO z = 1, nlevs
    PRINT *, 'Computing implied std for u level', z
    ! Put a point at an arbitrary horiz point in the u-field at level z
    Model_data1 % u(longindex(1), z) = 1.0
    ! Operate with U transpose
    CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    ! Store the StdDev
    StdDevu(z) = SQRT(Model_data2 % u(longindex(1), z))
    ! Replace with zero
    Model_data1 % u(longindex(1), z) = 0.0
  END DO

  ! v field
  DO z = 1, nlevs
    PRINT *, 'Computing implied std for v level', z
    ! Put a point at an arbitrary horiz point in the v-field at level z
    Model_data1 % v(longindex(1), z) = 1.0
    ! Operate with U transpose
    CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    ! Store the StdDev
    StdDevv(z) = SQRT(Model_data2 % v(longindex(1), z))
    ! Replace with zero
    Model_data1 % v(longindex(1), z) = 0.0
  END DO

  ! w field
  DO z = 1, nlevs
    PRINT *, 'Computing implied std for w level', z
    ! Put a point at an arbitrary horiz point in the w-field at level z
    Model_data1 % w(longindex(1), z) = 1.0
    ! Operate with U transpose
    CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    ! Store the StdDev
    StdDevw(z) = SQRT(Model_data2 % w(longindex(1), z))
    ! Replace with zero
    Model_data1 % w(longindex(1), z) = 0.0
  END DO

  ! r field
  DO z = 1, nlevs
    PRINT *, 'Computing implied std for r level', z
    ! Put a point at an arbitrary horiz point in the r-field at level z
    Model_data1 % r(longindex(1), z) = 1.0
    ! Operate with U transpose
    CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    ! Store the StdDev
    StdDevr(z) = SQRT(Model_data2 % r(longindex(1), z))
    ! Replace with zero
    Model_data1 % r(longindex(1), z) = 0.0
  END DO

  ! b field
  DO z = 1, nlevs
    PRINT *, 'Computing implied std for b level', z
    ! Put a point at an arbitrary horiz point in the b-field at level z
    Model_data1 % b(longindex(1), z) = 1.0
    ! Operate with U transpose
    CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    ! Store the StdDev
    StdDevb(z) = SQRT(Model_data2 % b(longindex(1), z))
    ! Replace with zero
    Model_data1 % b(longindex(1), z) = 0.0
  END DO




  PRINT *, '===== Running some point-by-point covariances ====='

  DO point = 1, ImplCov_npoints

    CALL Initialise_model_vars (Model_data1, .FALSE.)

    ! Put delta function in the u field
    ! ---------------------------------
    PRINT *, '  Point', point, ' delta u'
    Model_data1 % u(longindex(point), levindex(point)) = 1.0
    ! Operate with U transpose
    CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    ! Convert covariances into correlations
    DO x = 1, nlongs
      Model_data2 % u(x, 1:nlevs) = Model_data2 % u(x, 1:nlevs) / &
                                    (StdDevu(1:nlevs) * StdDevu(levindex(point)))
      Model_data2 % v(x, 1:nlevs) = Model_data2 % v(x, 1:nlevs) / &
                                    (StdDevv(1:nlevs) * StdDevu(levindex(point)))
      Model_data2 % w(x, 1:nlevs) = Model_data2 % w(x, 1:nlevs) / &
                                    (StdDevw(1:nlevs) * StdDevu(levindex(point)))
      Model_data2 % r(x, 1:nlevs) = Model_data2 % r(x, 1:nlevs) / &
                                    (StdDevr(1:nlevs) * StdDevu(levindex(point)))
      Model_data2 % b(x, 1:nlevs) = Model_data2 % b(x, 1:nlevs) / &
                                    (StdDevb(1:nlevs) * StdDevu(levindex(point)))
    END DO
    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltau.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)


    ! Put delta function in the v field
    ! ---------------------------------
    PRINT *, '  Point', point, ' delta v'
    Model_data1 % u(longindex(point), levindex(point)) = 0.0
    Model_data1 % v(longindex(point), levindex(point)) = 1.0
    ! Operate with U transpose
    CALL Initialise_CVs (CV_data, .FALSE.)
    CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    ! Operate with U
    CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    ! Convert covariances into correlations
    DO x = 1, nlongs
      Model_data2 % u(x, 1:nlevs) = Model_data2 % u(x, 1:nlevs) / &
                                    (StdDevu(1:nlevs) * StdDevv(levindex(point)))
      Model_data2 % v(x, 1:nlevs) = Model_data2 % v(x, 1:nlevs) / &
                                    (StdDevv(1:nlevs) * StdDevv(levindex(point)))
      Model_data2 % w(x, 1:nlevs) = Model_data2 % w(x, 1:nlevs) / &
                                    (StdDevw(1:nlevs) * StdDevv(levindex(point)))
      Model_data2 % r(x, 1:nlevs) = Model_data2 % r(x, 1:nlevs) / &
                                    (StdDevr(1:nlevs) * StdDevv(levindex(point)))
      Model_data2 % b(x, 1:nlevs) = Model_data2 % b(x, 1:nlevs) / &
                                    (StdDevb(1:nlevs) * StdDevv(levindex(point)))
    END DO
    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltav.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)


    ! Put delta function in the w field
    ! ---------------------------------
    PRINT *, '  Point', point, ' delta w'
    Model_data1 % v(longindex(point), levindex(point)) = 0.0
    Model_data1 % w(longindex(point), levindex(point)) = 1.0
    ! Operate with U transpose
    CALL Initialise_CVs (CV_data, .FALSE.)
    CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    ! Operate with U
    CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    ! Convert covariances into correlations
    DO x = 1, nlongs
      Model_data2 % u(x, 1:nlevs) = Model_data2 % u(x, 1:nlevs) / &
                                    (StdDevu(1:nlevs) * StdDevw(levindex(point)))
      Model_data2 % v(x, 1:nlevs) = Model_data2 % v(x, 1:nlevs) / &
                                    (StdDevv(1:nlevs) * StdDevw(levindex(point)))
      Model_data2 % w(x, 1:nlevs) = Model_data2 % w(x, 1:nlevs) / &
                                    (StdDevw(1:nlevs) * StdDevw(levindex(point)))
      Model_data2 % r(x, 1:nlevs) = Model_data2 % r(x, 1:nlevs) / &
                                    (StdDevr(1:nlevs) * StdDevw(levindex(point)))
      Model_data2 % b(x, 1:nlevs) = Model_data2 % b(x, 1:nlevs) / &
                                    (StdDevb(1:nlevs) * StdDevw(levindex(point)))
    END DO
    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltaw.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)


    ! Put delta function in the r field
    ! ---------------------------------
    PRINT *, '  Point', point, ' delta r'
    Model_data1 % w(longindex(point), levindex(point)) = 0.0
    Model_data1 % r(longindex(point), levindex(point)) = 1.0
    ! Operate with U transpose
    CALL Initialise_CVs (CV_data, .FALSE.)
    CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    ! Operate with U
    CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    ! Convert covariances into correlations
    DO x = 1, nlongs
      Model_data2 % u(x, 1:nlevs) = Model_data2 % u(x, 1:nlevs) / &
                                    (StdDevu(1:nlevs) * StdDevr(levindex(point)))
      Model_data2 % v(x, 1:nlevs) = Model_data2 % v(x, 1:nlevs) / &
                                    (StdDevv(1:nlevs) * StdDevr(levindex(point)))
      Model_data2 % w(x, 1:nlevs) = Model_data2 % w(x, 1:nlevs) / &
                                    (StdDevw(1:nlevs) * StdDevr(levindex(point)))
      Model_data2 % r(x, 1:nlevs) = Model_data2 % r(x, 1:nlevs) / &
                                    (StdDevr(1:nlevs) * StdDevr(levindex(point)))
      Model_data2 % b(x, 1:nlevs) = Model_data2 % b(x, 1:nlevs) / &
                                    (StdDevb(1:nlevs) * StdDevr(levindex(point)))
    END DO
    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltar.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)


    ! Put delta function in the b field
    ! ---------------------------------
    PRINT *, '  Point', point, ' delta b'
    Model_data1 % r(longindex(point), levindex(point)) = 0.0
    Model_data1 % b(longindex(point), levindex(point)) = 1.0
    ! Operate with U transpose
    CALL Initialise_CVs (CV_data, .FALSE.)
    CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    ! Operate with U
    CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    ! Convert covariances into correlations
    DO x = 1, nlongs
      Model_data2 % u(x, 1:nlevs) = Model_data2 % u(x, 1:nlevs) / &
                                    (StdDevu(1:nlevs) * StdDevb(levindex(point)))
      Model_data2 % v(x, 1:nlevs) = Model_data2 % v(x, 1:nlevs) / &
                                    (StdDevv(1:nlevs) * StdDevb(levindex(point)))
      Model_data2 % w(x, 1:nlevs) = Model_data2 % w(x, 1:nlevs) / &
                                    (StdDevw(1:nlevs) * StdDevb(levindex(point)))
      Model_data2 % r(x, 1:nlevs) = Model_data2 % r(x, 1:nlevs) / &
                                    (StdDevr(1:nlevs) * StdDevb(levindex(point)))
      Model_data2 % b(x, 1:nlevs) = Model_data2 % b(x, 1:nlevs) / &
                                    (StdDevb(1:nlevs) * StdDevb(levindex(point)))
    END DO
    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltab.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)


    ! Put delta function in the tracer field
    ! ---------------------------------
    PRINT *, '  Point', point, ' delta b'
    Model_data1 % b(longindex(point), levindex(point)) = 0.0
    Model_data1 % tracer(longindex(point), levindex(point)) = 1.0
    ! Operate with U transpose
    CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltatracer.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)

  END DO
END IF

CALL Deallocate_dims(dims)
CALL Deallocate_model_vars(Model_data1)
CALL Deallocate_model_vars(Model_data2)
CALL Deallocate_model_vars(LS)
CALL Deallocate_CVs(CV_data)
CALL Deallocate_CVT(CVT)
DEALLOCATE (StdDevu, StdDevv, StdDevw, StdDevr, StdDevb)
IF (ALLOCATED(fft_wsave_x)) DEALLOCATE(fft_wsave_x)
IF (ALLOCATED(fft_work_x)) DEALLOCATE(fft_work_x)

END PROGRAM Util_ImpliedCor
