PROGRAM Master_ImpliedCov

!*****************************************************
!*   Code to compute a selection of implied          *
!*   background error covariances.                   *
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
    aCVT_type,                   &
    datadirImpliedCov,           &
    datadirABCfcs,               &
    Hybrid_opt,                  &
    NEnsMems,                    &
    Nlats,                       &
    Cov_WeightE,                 &
    Cov_WeightC,                 &
    VarLoc,                      &
    Var_dep_loc,                 &
    hScale_alpha,                &
    vScale_alpha,                &
    datadirEM,                   &
    LS_file,                     &
    datadirCVT,                  &
    CVT_file,                    &
    ImplCov_npoints,             &
    longindex,                   &
    levindex,                    &
    Npointsmax,                  &
    fft_wsave_x,                 &
    fft_work_x, Nvariable


IMPLICIT NONE

! Declare variables
!==========================
TYPE(dims_type)          :: dims
TYPE(ABC_type)           :: Model_data1, Model_data2, LS, dxb, dxe
TYPE(CV_type)            :: CV_data
TYPE(CVT_type)           :: CVT
TYPE(aCVT_type)          :: L_alpha
INTEGER                  :: point, Nacv, n, Neffmems, item
LOGICAL                  :: Use_EOTD, exist

CHARACTER(LEN=320)          :: LS_filename, CVT_filename, ImplCov_filename, ABCfilename
TYPE(ABC_type), ALLOCATABLE :: EM_x(:), EM_x_copy(:), tmp_Uv(:)
TYPE(ABC_type), ALLOCATABLE :: CV_alpha_data(:)
REAL(ZREAL8)                :: denom, beta_e, beta_c
REAl(ZREAL8), ALLOCATABLE   :: Work(:,:,:)


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_ImpliedCov'
PRINT*, '*************************************************************************'

! Read namelist
CALL SetOptions

CALL Initialise_dims (dims)
CALL Initialise_model_vars (Model_data1, .FALSE.)
CALL Initialise_model_vars (Model_data2, .FALSE.)
CALL Initialise_model_vars (LS, .FALSE.)
CALL Initialise_CVs (CV_data, .FALSE.)
CALL Initialise_CVT (CVT)

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

! Check if hybrid is used and allocate appropriate arrays
IF (Hybrid_opt == 2 .OR. Hybrid_opt == 3) THEN
  Use_EOTD = .TRUE.
  Neffmems = NEnsMems
  Nacv = Neffmems

  ! Ensure number of effective ensemble members are more than 1
  IF (Neffmems > 1) THEN
    PRINT*, 'Using', Neffmems, 'ensemble error modes'
  ELSE
    PRINT*, 'EnVar enabled but invalid number of ensemble error modes'
    STOP
  END IF

  ! Set weighting factors
  beta_c = SQRT(Cov_WeightC/100)
  beta_e = SQRT(Cov_WeightE/100)

  IF (Hybrid_opt == 2) THEN
    ! Pure EnVar so forcing static weighting to 0
    ! Set weighting factors
    beta_c = 0
  END IF

  PRINT *, 'Beta (climatological) is ', beta_c
  PRINT *, 'Beta (ensemble) is ', beta_e
ELSE
  Use_EOTD = .FALSE.
  ! Set Nacv with a dummy value of 1 when hybrid is not enabled
  Nacv = 1
END IF

! Allocate arrays
ALLOCATE (EM_x(1:Nacv))
ALLOCATE (EM_x_copy(1:Nacv))
ALLOCATE (CV_alpha_data(1:Nacv))
ALLOCATE (tmp_Uv(1:Nacv))
DO n = 1, Nacv
  CALL Initialise_model_vars (CV_alpha_data(n), .FALSE.)
  CALL Initialise_model_vars (tmp_Uv(n), .FALSE.)
END DO

! Set up localisation matrix
IF (Use_EOTD) THEN
  CALL Initialise_aCVT(L_alpha)
  IF (Var_dep_loc) THEN
    CALL SetHorizLoc_alpha_vdl (hScale_alpha, dims, L_alpha)
    CALL SetVertLoc_alpha_vdl (vScale_alpha, dims, L_alpha)
  ELSE
    CALL SetHorizLoc_alpha (hScale_alpha, dims, L_alpha)
    CALL SetVertLoc_alpha (vScale_alpha, dims, L_alpha)
  END IF
END IF

! Count number of error mode files and ensure it matches number of effective ensemble members
IF (Use_EOTD) THEN
  item = 0
  DO n = 1, Nacv
    ! Check if file exists in directory
    WRITE (ABCfilename, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirEM), '/PertABC_Item', n, '.nc'
    INQUIRE(file=ABCfilename, exist=exist)
    IF (exist) THEN
      ! Read in error mode states to a model_state_array if file exists
      CALL Initialise_model_vars (EM_x(n), .FALSE.)
      CALL Read_state_2d (ABCfilename, EM_x(n), dims, -1, .TRUE.)
      ! Divide by number of Nacv-1
      denom = SQRT(FLOAT(Nacv-1))
      CALL Div_model_cons (EM_x(n), denom, .TRUE.)
    ELSE
      EXIT
    END IF
    item = item + 1
  END DO 
END IF

IF (Use_EOTD) THEN
  IF (Nacv .NE. item) THEN
        ! Usually due to missing indexes
    PRINT*, 'Number of error modes files available does not match number of effective ensemble members'
    STOP
  END IF
END IF


IF (ImplCov_npoints > 0) THEN
  PRINT *, '===== Running some point-by-point covariances ====='

  DO point = 1, ImplCov_npoints

    CALL Initialise_model_vars (Model_data1, .FALSE.)

    ! Put delta function in the u field
    ! ---------------------------------
    PRINT *, '  Point', point, ' delta u'
    Model_data1 % u(longindex(point), levindex(point)) = 1.0

    IF (Var_dep_loc) THEN
      ALLOCATE (Work(1:Nvariable*nlongs,1:Nvariable*nlevs,1:Nacv))
      Work(1:Nvariable*nlongs,1:Nvariable*nlevs,1:Nacv) = 0
    END IF

    ! Operate with U transpose
    ! Hybrid implied covariances
    IF (Use_EOTD) THEN
      DO n = 1, Nacv
        ! Copy error modes data so original does not get modified
        CALL Initialise_model_vars (EM_x_copy(n), .FALSE.)
        CALL Add_model_vars (EM_x_copy(n), EM_x(n), .FALSE.)

        CALL Mul_model_vars (EM_x_copy(n), Model_data1, .TRUE.)

        IF (Var_dep_loc) THEN ! Var_dep_loc deals with full control vector explicitly
          CALL U_trans_alpha_vdl_adj (LS, Work(:,:,n), EM_x_copy(n), CVT, dims, L_alpha)

        ELSE
          CALL U_trans_alpha_adj (LS, CV_alpha_data(n), EM_x_copy(n), CVT, dims, L_alpha)
          IF (VarLoc == 2) THEN
            ! Same alpha fields applied on each variable
            ! Extracting the column of the implied covariance associated with u at a point
            CV_alpha_data(n) % v (:,:) = CV_alpha_data(n) % u (:,:)
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % u (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % u (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % u (:,:)
          ELSE IF (VarLoc == 3) THEN          
            ! Same alpha fields applied on u and v
            ! Extracting the column of the implied covariance associated with u at a point
            CV_alpha_data(n) % v (:,:) = CV_alpha_data(n) % u (:,:)
          ELSE IF (VarLoc == 4) THEN
            ! Same alpha fields applied on each variable (only v excluded)
            ! Extracting the column of the implied covariance associated with u at a point
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % u (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % u (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % u (:,:)
          ELSE IF (VarLoc == 5) THEN
            ! Same alpha fields applied on each variable (only r excluded)
            ! Extracting the column of the implied covariance associated with u at a point
            CV_alpha_data(n) % v (:,:) = CV_alpha_data(n) % u (:,:)
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % u (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % u (:,:)
          ELSE IF (VarLoc == 6) THEN
            ! Same alpha fields applied on each variable (only v and w excluded)
            ! Extracting the column of the implied covariance associated with u at a point
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % u (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % u (:,:)
          !ELSE IF (VarLoc == 7) THEN
            ! Same alpha fields applied on each variable (only u excluded)
            ! Extracting the column of the implied covariance associated with u at a point
          END IF
        END IF

        ! Apply weighting
        IF (Var_dep_loc) THEN
          Work(:,:,n) = Work(:,:,n)*beta_e
        ELSE 
          CALL Mul_model_cons (CV_alpha_data(n), beta_e, .TRUE.)
        END IF
      END DO
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      
      CALL Mul_CV_cons (CV_data, beta_c)
    ! Normal implied covariances
    ELSE
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    END IF

    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    ! Hybrid implied covariances
    IF (Use_EOTD) THEN
      CALL Initialise_model_vars (dxb, .FALSE.)
      CALL Initialise_model_vars (dxe, .FALSE.)
      
      DO n = 1, Nacv
        IF (Var_dep_loc) THEN
          CALL U_trans_alpha_vdl (LS, Work(:,:,n), tmp_Uv(n), CVT, dims, L_alpha)
        ELSE
          CALL U_trans_alpha (LS, CV_alpha_data(n), tmp_Uv(n), CVT, dims, L_alpha)
        END IF
        ! Schur product multiplication between error modes and alpha control variables in model space

        CALL Mul_model_vars(tmp_Uv(n), EM_x(n), .TRUE.)
        CALL Add_model_vars(dxe, tmp_Uv(n), .TRUE.)          
      END DO
      CALL Mul_model_cons(dxe, beta_e, .TRUE.)

      CALL U_trans (LS, CV_data, dxb, CVT, dims)
      CALL Mul_model_cons (dxb, beta_c, .TRUE.)

      CALL Add_model_vars (Model_data2, dxb, .TRUE.)
      CALL Add_model_vars (Model_data2, dxe, .TRUE.)
    ! Normal implied covariances
    ELSE
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    END IF

    PRINT *, Model_data2 % u(longindex(point), levindex(point))

    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltau.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)

    ! Put delta function in the v field
    ! ---------------------------------
    PRINT *, '  Point', point, ' delta v'
    Model_data1 % u(longindex(point), levindex(point)) = 0.0
    Model_data1 % v(longindex(point), levindex(point)) = 1.0

    IF (Var_dep_loc) THEN
      Work(1:Nvariable*nlongs,1:Nvariable*nlevs,1:Nacv) = 0
    END IF

    ! Operate with U transpose  
    ! Hybrid implied covariances
    IF (Use_EOTD) THEN
      DO n = 1, Nacv
        ! Copy error modes data so original does not get modified
        CALL Initialise_model_vars (EM_x_copy(n), .FALSE.)
        CALL Add_model_vars (EM_x_copy(n), EM_x(n), .FALSE.)

        CALL Mul_model_vars (EM_x_copy(n), Model_data1, .TRUE.)

        IF (Var_dep_loc) THEN ! Var_dep_loc deals with full control vector explicitly
          CALL U_trans_alpha_vdl_adj (LS, Work(:,:,n), EM_x_copy(n), CVT, dims, L_alpha)

        ELSE
          CALL U_trans_alpha_adj (LS, CV_alpha_data(n), EM_x_copy(n), CVT, dims, L_alpha)
          IF (VarLoc == 2) THEN
            ! Same alpha fields applied on each variable
            ! Extracting the column of the implied covariance associated with v at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % v (:,:)
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % v (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % v (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % v (:,:)
          ELSE IF (VarLoc == 3) THEN
            ! Same alpha fields applied on u and v
              ! Extracting the column of the implied covariance associated with u at
            ! a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % v (:,:)
          !ELSE IF (VarLoc == 4) THEN
            ! Same alpha fields applied on each variable (only v excluded)
            ! Extracting the column of the implied covariance associated with v at a point
          ELSE IF (VarLoc == 5) THEN
            ! Same alpha fields applied on each variable (only r excluded)
            ! Extracting the column of the implied covariance associated with v at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % v (:,:)
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % v (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % v (:,:)
          !ELSE IF (VarLoc == 6) THEN
            ! Same alpha fields applied on each variable (only v and w excluded)
            ! Extracting the column of the implied covariance associated with v at a point
          ELSE IF (VarLoc == 7) THEN
            ! Same alpha fields applied on each variable (only u excluded)
            ! Extracting the column of the implied covariance associated with v at a point
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % v (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % v (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % v (:,:)
          END IF
        END IF
   
        ! Apply weighting
        IF (Var_dep_loc) THEN
          Work(:,:,n) = Work(:,:,n)*beta_e
        ELSE
          CALL Mul_model_cons (CV_alpha_data(n), beta_e, .TRUE.)
        END IF
      END DO
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      CALL Mul_CV_cons (CV_data, beta_c)
    ! Normal implied covariances
    ELSE
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    END IF

    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    ! Hybrid implied covariances
    IF (Use_EOTD) THEN
      CALL Initialise_model_vars (dxb, .FALSE.)
      CALL Initialise_model_vars (dxe, .FALSE.)
      DO n = 1, Nacv
        IF (Var_dep_loc) THEN
          CALL U_trans_alpha_vdl (LS, Work(:,:,n), tmp_Uv(n), CVT, dims, L_alpha)
        ELSE
          CALL U_trans_alpha (LS, CV_alpha_data(n), tmp_Uv(n), CVT, dims, L_alpha)
        END IF
        ! Schur product multiplication between error modes and alpha control variables in model space

        CALL Mul_model_vars(tmp_Uv(n), EM_x(n), .TRUE.)
        CALL Add_model_vars(dxe, tmp_Uv(n), .TRUE.)          
      END DO
      CALL Mul_model_cons(dxe, beta_e, .TRUE.)

      CALL U_trans (LS, CV_data, dxb, CVT, dims)
      CALL Mul_model_cons (dxb, beta_c, .TRUE.)

      CALL Add_model_vars (Model_data2, dxb, .TRUE.)
      CALL Add_model_vars (Model_data2, dxe, .TRUE.)
    ! Normal implied covariances
    ELSE
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    END IF

    PRINT *, Model_data2 % v(longindex(point), levindex(point))

    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltav.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)


    ! Put delta function in the w field
    ! ---------------------------------
    PRINT *, '  Point', point, ' delta w'
    Model_data1 % v(longindex(point), levindex(point)) = 0.0
    Model_data1 % w(longindex(point), levindex(point)) = 1.0

    ! Operate with U transpose
    ! Hybrid implied covariances
    IF (Use_EOTD) THEN
      DO n = 1, Nacv
        ! Copy error modes data so original does not get modified
        CALL Initialise_model_vars (EM_x_copy(n), .FALSE.)
        CALL Add_model_vars (EM_x_copy(n), EM_x(n), .FALSE.)

        CALL Mul_model_vars (EM_x_copy(n), Model_data1, .TRUE.)

        IF (Var_dep_loc) THEN ! Var_dep_loc deals with full control vector explicitly
          CALL U_trans_alpha_vdl_adj (LS, Work(:,:,n), EM_x_copy(n), CVT, dims, L_alpha)

        ELSE
          CALL U_trans_alpha_adj (LS, CV_alpha_data(n), EM_x_copy(n), CVT, dims, L_alpha)
          IF (VarLoc == 2) THEN
            ! Same alpha fields applied on each variable
            ! Extracting the column of the implied covariance associated with w at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % w (:,:)
            CV_alpha_data(n) % v (:,:) = CV_alpha_data(n) % w (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % w (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % w (:,:)
          ELSE IF (VarLoc == 3) THEN
            ! Same alpha fields applied on w,r,b
            ! Extracting the column of the implied covariance associated with w at
            ! a point
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % w (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % w (:,:)
          ELSE IF (VarLoc == 4) THEN
            ! Same alpha fields applied on each variable (only v excluded)
            ! Extracting the column of the implied covariance associated with w at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % w (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % w (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % w (:,:)
          ELSE IF (VarLoc == 5) THEN
            ! Same alpha fields applied on each variable (only r excluded)
            ! Extracting the column of the implied covariance associated with w at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % w (:,:)
            CV_alpha_data(n) % v (:,:) = CV_alpha_data(n) % w (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % w (:,:)
          !ELSE IF (VarLoc == 6) THEN
            ! Same alpha fields applied on each variable (only v and w excluded)
            ! Extracting the column of the implied covariance associated with w at a point
          ELSE IF (VarLoc == 7) THEN
            ! Same alpha fields applied on each variable (only u excluded)
            ! Extracting the column of the implied covariance associated with w at a point
            CV_alpha_data(n) % v (:,:) = CV_alpha_data(n) % w (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % w (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % w (:,:)
          END IF
        END IF

        ! Apply weighting
        IF (Var_dep_loc) THEN
          Work(:,:,n) = Work(:,:,n)*beta_e
        ELSE
          CALL Mul_model_cons (CV_alpha_data(n), beta_e, .TRUE.)
        END IF
      END DO
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      CALL Mul_CV_cons (CV_data, beta_c)
    ! Normal implied covariances
    ELSE
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    END IF

    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    ! Hybrid implied covariances
    IF (Use_EOTD) THEN
      CALL Initialise_model_vars (dxb, .FALSE.)
      CALL Initialise_model_vars (dxe, .FALSE.)
      DO n = 1, Nacv
        IF (Var_dep_loc) THEN
          CALL U_trans_alpha_vdl (LS, Work(:,:,n), tmp_Uv(n), CVT, dims, L_alpha)
        ELSE
          CALL U_trans_alpha (LS, CV_alpha_data(n), tmp_Uv(n), CVT, dims, L_alpha)
        END IF

        ! Schur product multiplication between error modes and alpha control variables in model space
        CALL Mul_model_vars(tmp_Uv(n), EM_x(n), .TRUE.)
        CALL Add_model_vars(dxe, tmp_Uv(n), .TRUE.)          
      END DO
      CALL Mul_model_cons(dxe, beta_e, .TRUE.)

      CALL U_trans (LS, CV_data, dxb, CVT, dims)
      CALL Mul_model_cons (dxb, beta_c, .TRUE.)

      CALL Add_model_vars (Model_data2, dxb, .TRUE.)
      CALL Add_model_vars (Model_data2, dxe, .TRUE.)
    ! Normal implied covariances
    ELSE
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    END IF

    PRINT *, Model_data2 % w(longindex(point), levindex(point))

    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltaw.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)


    ! Put delta function in the r field
    ! ---------------------------------
    PRINT *, '  Point', point, ' delta r'
    Model_data1 % w(longindex(point), levindex(point)) = 0.0
    Model_data1 % r(longindex(point), levindex(point)) = 1.0

    ! Operate with U transpose
    ! Hybrid implied covariances
    IF (Use_EOTD) THEN
      DO n = 1, Nacv
        ! Copy error modes data so original does not get modified
        CALL Initialise_model_vars (EM_x_copy(n), .FALSE.)
        CALL Add_model_vars (EM_x_copy(n), EM_x(n), .FALSE.)

        CALL Mul_model_vars (EM_x_copy(n), Model_data1, .TRUE.)

        IF (Var_dep_loc) THEN ! Var_dep_loc deals with full control vector explicitly
          CALL U_trans_alpha_vdl_adj (LS, Work(:,:,n), EM_x_copy(n), CVT, dims, L_alpha)

        ELSE
          CALL U_trans_alpha_adj (LS, CV_alpha_data(n), EM_x_copy(n), CVT, dims, L_alpha)
          IF (VarLoc == 2) THEN
            ! Same alpha fields applied on each variable
            ! Extracting the column of the implied covariance associated with r at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % r (:,:)
            CV_alpha_data(n) % v (:,:) = CV_alpha_data(n) % r (:,:)
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % r (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % r (:,:)
          ELSE IF (VarLoc == 3) THEN
            ! Same alpha fields applied on w,r,b
            ! Extracting the column of the implied covariance associated with r at a point
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % r (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % r (:,:)
          ELSE IF (VarLoc == 4) THEN
            ! Same alpha fields applied on each variable (only v excluded)
            ! Extracting the column of the implied covariance associated with r at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % r (:,:)
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % r (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % r (:,:)
            !ELSE IF (VarLoc == 5) THEN
            ! Same alpha fields applied on each variable (only r excluded)
            ! Extracting the column of the implied covariance associated with r at a point
          ELSE IF (VarLoc == 6) THEN
            ! Same alpha fields applied on each variable (only v and w excluded)
            ! Extracting the column of the implied covariance associated with r at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % r (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % r (:,:)
          ELSE IF (VarLoc == 7) THEN
            ! Same alpha fields applied on each variable (only u excluded)
            ! Extracting the column of the implied covariance associated with r at a point
            CV_alpha_data(n) % v (:,:) = CV_alpha_data(n) % r (:,:)
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % r (:,:)
            CV_alpha_data(n) % b (:,:) = CV_alpha_data(n) % r (:,:)
          END IF
        END IF

        ! Apply weighting
        IF (Var_dep_loc) THEN
          Work(:,:,n) = Work(:,:,n)*beta_e
        ELSE
          CALL Mul_model_cons (CV_alpha_data(n), beta_e, .TRUE.)
        END IF
      END DO
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      CALL Mul_CV_cons (CV_data, beta_c)
    ! Normal implied covariances
    ELSE
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    END IF

    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    ! Hybrid implied covariances
    IF (Use_EOTD) THEN
      CALL Initialise_model_vars (dxb, .FALSE.)
      CALL Initialise_model_vars (dxe, .FALSE.)
      DO n = 1, Nacv
        CALL U_trans_alpha (LS, CV_alpha_data(n), tmp_Uv(n), CVT, dims, L_alpha)
       
        IF (Var_dep_loc) THEN
          CALL U_trans_alpha_vdl (LS, Work(:,:,n), tmp_Uv(n), CVT, dims, L_alpha)
        ELSE
          CALL U_trans_alpha (LS, CV_alpha_data(n), tmp_Uv(n), CVT, dims, L_alpha)
        END IF

        ! Schur product multiplication between error modes and alpha control variables in model space 
        CALL Mul_model_vars(tmp_Uv(n), EM_x(n), .TRUE.)
        CALL Add_model_vars(dxe, tmp_Uv(n), .TRUE.)          
      END DO
      CALL Mul_model_cons(dxe, beta_e, .TRUE.)

      CALL U_trans (LS, CV_data, dxb, CVT, dims)
      CALL Mul_model_cons (dxb, beta_c, .TRUE.)

      CALL Add_model_vars (Model_data2, dxb, .TRUE.)
      CALL Add_model_vars (Model_data2, dxe, .TRUE.)
    ! Normal implied covariances
    ELSE
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    END IF

    PRINT *, Model_data2 % r(longindex(point), levindex(point))

    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltar.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)


    ! Put delta function in the b field
    ! ---------------------------------
    PRINT *, '  Point', point, ' delta b'
    Model_data1 % r(longindex(point), levindex(point)) = 0.0
    Model_data1 % b(longindex(point), levindex(point)) = 1.0

    ! Operate with U transpose
    ! Hybrid implied covariances
    IF (Use_EOTD) THEN
      DO n = 1, Nacv
        ! Copy error modes data so original does not get modified
        CALL Initialise_model_vars (EM_x_copy(n), .FALSE.)
        CALL Add_model_vars (EM_x_copy(n), EM_x(n), .FALSE.)

        CALL Mul_model_vars (EM_x_copy(n), Model_data1, .TRUE.)
        
        IF (Var_dep_loc) THEN ! Var_dep_loc deals with full control vector explicitly
          CALL U_trans_alpha_vdl_adj (LS, Work(:,:,n), EM_x_copy(n), CVT, dims, L_alpha)

        ELSE
          CALL U_trans_alpha_adj (LS, CV_alpha_data(n), EM_x_copy(n), CVT, dims, L_alpha)
          IF (VarLoc == 2) THEN
            ! Same alpha fields applied on each variable
            ! Extracting the column of the implied covariance associated with b at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % b (:,:)
            CV_alpha_data(n) % v (:,:) = CV_alpha_data(n) % b (:,:)
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % b (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % b (:,:)
          ELSE IF (VarLoc == 3) THEN
            ! Same alpha fields applied on w,r,b
            ! Extracting the column of the implied covariance associated with b at a point
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % b (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % b (:,:)
          ELSE IF (VarLoc == 4) THEN
            ! Same alpha fields applied on each variable (only v excluded)
            ! Extracting the column of the implied covariance associated with b at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % b (:,:)
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % b (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % b (:,:)
          ELSE IF (VarLoc == 5) THEN
            ! Same alpha fields applied on each variable (only r excluded)
            ! Extracting the column of the implied covariance associated with b at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % b (:,:)
            CV_alpha_data(n) % v (:,:) = CV_alpha_data(n) % b (:,:)
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % b (:,:)
          ELSE IF (VarLoc == 6) THEN
            ! Same alpha fields applied on each variable (only v and w excluded)
            ! Extracting the column of the implied covariance associated with b at a point
            CV_alpha_data(n) % u (:,:) = CV_alpha_data(n) % b (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % b (:,:) 
          ELSE IF (VarLoc == 7) THEN
            ! Same alpha fields applied on each variable (only u excluded)
            ! Extracting the column of the implied covariance associated with b at a point
            CV_alpha_data(n) % v (:,:) = CV_alpha_data(n) % b (:,:)
            CV_alpha_data(n) % w (:,:) = CV_alpha_data(n) % b (:,:)
            CV_alpha_data(n) % r (:,:) = CV_alpha_data(n) % b (:,:)
          END IF
        END IF

        ! Apply weighting
        IF (Var_dep_loc) THEN
          Work(:,:,n) = Work(:,:,n)*beta_e
        ELSE
          CALL Mul_model_cons (CV_alpha_data(n), beta_e, .TRUE.)
        END IF
      END DO
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      CALL Mul_CV_cons (CV_data, beta_c)
    ! Normal implied covariances
    ELSE
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    END IF

    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
    ! Hybrid implied covariances
    IF (Use_EOTD) THEN
      CALL Initialise_model_vars (dxb, .FALSE.)
      CALL Initialise_model_vars (dxe, .FALSE.)
      DO n = 1, Nacv
        
        IF (Var_dep_loc) THEN
          CALL U_trans_alpha_vdl (LS, Work(:,:,n), tmp_Uv(n), CVT, dims, L_alpha)
        ELSE
          CALL U_trans_alpha (LS, CV_alpha_data(n), tmp_Uv(n), CVT, dims, L_alpha)
        END IF

        ! Schur product multiplication between error modes and alpha control variables in model space
        CALL Mul_model_vars(tmp_Uv(n), EM_x(n), .TRUE.)
        CALL Add_model_vars(dxe, tmp_Uv(n), .TRUE.)          
      END DO
      CALL Mul_model_cons(dxe, beta_e, .TRUE.)

      CALL U_trans (LS, CV_data, dxb, CVT, dims)
      CALL Mul_model_cons (dxb, beta_c, .TRUE.)

      CALL Add_model_vars (Model_data2, dxb, .TRUE.)
      CALL Add_model_vars (Model_data2, dxe, .TRUE.)
    ! Normal implied covariances
    ELSE
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    END IF

    PRINT *, Model_data2 % b(longindex(point), levindex(point))

    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltab.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)


    ! Put delta function in the tracer field
    ! ---------------------------------
    ! Many features not implemented in tracer field
    PRINT *, '  Point', point, ' delta tracer'
    Model_data1 % b(longindex(point), levindex(point)) = 0.0
    Model_data1 % tracer(longindex(point), levindex(point)) = 1.0

    ! Operate with U transpose
	! Hybrid implied covariances
    IF (Use_EOTD) THEN
      DO n = 1, Nacv
        ! Copy error modes data so original does not get modified
        CALL Initialise_model_vars (EM_x_copy(n), .FALSE.)
        CALL Add_model_vars (EM_x_copy(n), EM_x(n), .FALSE.)

        CALL Mul_model_vars (EM_x_copy(n), Model_data1, .TRUE.)
        CALL U_trans_alpha_adj (LS, CV_alpha_data(n), EM_x_copy(n), CVT, dims, L_alpha)
        ! Apply weighting
        CALL Mul_model_cons (CV_alpha_data(n), beta_e, .TRUE.)
      END DO
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      CALL Mul_CV_cons (CV_data, beta_c)
    ! Normal implied covariances
    ELSE
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
    END IF

    ! Operate with U
    CALL Initialise_model_vars (Model_data2, .FALSE.)
	! Hybrid implied covariances
    IF (Use_EOTD) THEN
      CALL Initialise_model_vars (dxb, .FALSE.)
      CALL Initialise_model_vars (dxe, .FALSE.)
      DO n = 1, Nacv
        CALL U_trans_alpha (LS, CV_alpha_data(n), tmp_Uv(n), CVT, dims, L_alpha)
        ! Schur product multiplication between error modes and alpha control variables in model space
        CALL Mul_model_vars(tmp_Uv(n), EM_x(n), .TRUE.)
        CALL Add_model_vars(dxe, tmp_Uv(n), .TRUE.)          
      END DO
      CALL Mul_model_cons(dxe, beta_e, .TRUE.)

      CALL U_trans (LS, CV_data, dxb, CVT, dims)
      CALL Mul_model_cons (dxb, beta_c, .TRUE.)

      CALL Add_model_vars (Model_data2, dxb, .TRUE.)
      CALL Add_model_vars (Model_data2, dxe, .TRUE.)
    ! Normal implied covariances
    ELSE
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
    END IF

    ! Output the result
    WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltatracer.nc'
    CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0, .TRUE.)

  END DO
END IF

! Deallocate
PRINT*, 'Deallocating ...'
DEALLOCATE (EM_x)
DEALLOCATE (EM_x_copy)
DEALLOCATE (CV_alpha_data)
DEALLOCATE (tmp_Uv)
CALL Deallocate_dims (dims)
CALL Deallocate_model_vars (Model_data1)
CALL Deallocate_model_vars (Model_data2)
CALL Deallocate_model_vars (LS)
CALL Deallocate_CVs (CV_data)
CALL Deallocate_CVT (CVT)
IF (ALLOCATED(fft_wsave_x)) DEALLOCATE(fft_wsave_x)
IF (ALLOCATED(fft_work_x)) DEALLOCATE(fft_work_x)
IF (Var_dep_loc) DEALLOCATE (Work)

END PROGRAM Master_ImpliedCov
