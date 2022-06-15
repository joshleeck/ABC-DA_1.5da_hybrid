PROGRAM Master_Assimilate

!*****************************************************
!*   Master code to do data assimilation with ABC    *
!*   model                                           *
!*                                                   *
!*   Ross Bannister, r.n.bannister@reading.ac.uk     *
!*   25/04/18                                        *
!*   Joshua Lee,     joshua_lee@nea.gov.sg           *
!*   06/06/21                                        *
!*       - added hybrid algorithm                    *
!*                                                   *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    nlongs, nlevs,               &
    Hybrid_opt,                  &
    NEnsMems,                    &
    Cov_WeightE,                 &
    Cov_WeightC,                 &
    VarLoc,                      &
    hScale_alpha,                &
    vScale_alpha,                &
    datadirEM,                   &
    Vartype,                     &
!    LS_file,                     &
    datadirCVT,                  &
    CVT_file,                    &
    datadirAnal,                 &
    anal_file,                   &
    analinc_file,                &
    datadir_Bg,                  &
    Bg_file,                     &
    datadir_Obs,                 &
    Obs_file,                    &
    N_outerloops,                &
    N_innerloops_max,            &
    dt, dt_da, t0,               &
    mu, minus_mu,                &
    crit_inner,                  &
    dims_type,                   &
    ABC_type,                    &
    CV_type,                     &
    CVT_type,                    &
    aCVT_type,                   &
    diagnostics_file,            &
    Obs_type,                    &
    fft_wsave_x,                 &
    fft_work_x




IMPLICIT NONE

INCLUDE "InnerProdControlSpace.interface"
INCLUDE "InnerProdModelSpace.interface"
INCLUDE "Read_Obs.interface"
INCLUDE "DeAllocate_Obs.interface"
INCLUDE "PenAndGrad.interface"


! Declare variables
!==========================
CHARACTER(LEN=320)          :: Bg_filename, CVT_filename, Obs_filename
CHARACTER(LEN=320)          :: Scalar_diags_filename, Anal_filename, Analinc_filename, ABCfilename
TYPE(ABC_type)              :: Bg0, dBg0, dx0, anal_inc
TYPE(ABC_type), ALLOCATABLE :: LSfc(:)
TYPE(dims_type)             :: dims
TYPE(CV_type)               :: dchiB0, dchi0_i, dchi0_ip1
TYPE(CV_type)               :: r_i, r_ip1, p_i, p_ip1, pert
TYPE(Obs_type), POINTER     :: Observations
INTEGER                     :: timesteps, DAtimesteps, maxtime, Neffmems, Nacv
TYPE(CVT_type)              :: CVT
TYPE(aCVT_type)             :: L_alpha
INTEGER                     :: k, i, item, n, t
LOGICAL                     :: dchiB0zero, dchi0zero, Converged_inner, ForceNLmodel, Use_EOTD, exist
REAL(ZREAL8)                :: Jb, Jo, Je, J, Jb_plus, Jo_plus, J_plus, Jb_minus, Jo_minus, J_minus, Je_minus, Je_plus
REAL(ZREAL8)                :: magnitude2_r_i, magnitude2_r_0, magnitude_adjust_inner, magnitude2_r_ip1
REAL(ZREAL8)                :: magnitude2_r_i_alpha, magnitude2_r_ip1_alpha
REAL(ZREAL8)                :: beta_i, alpha_i_i, beta_c, beta_e, denom
REAL(ZREAL8)                :: m1 = -1.0, N_p = 5.0

TYPE(ABC_type)              :: dxb0, dxe0
TYPE(ABC_type), ALLOCATABLE :: chi0_i_alpha(:), r_i_alpha(:), p_i_alpha(:), pert_alpha(:), r_i_alpha_copy(:)
TYPE(ABC_type), ALLOCATABLE :: chi0_ip1_alpha(:), r_ip1_alpha(:), p_ip1_alpha(:), r_ip1_alpha_copy(:)
TYPE(ABC_type), ALLOCATABLE :: EM_x(:), tmp_Uv(:)



PRINT*, '*************************************************************************'
PRINT*, 'Running Master_Assimilate'
PRINT*, '*************************************************************************'

! Read namelist
CALL SetOptions

CALL Initialise_dims (dims)
CALL Initialise_model_vars (Bg0, .FALSE.)
CALL Initialise_model_vars (dBg0, .FALSE.)
CALL Initialise_model_vars (dx0, .FALSE.)
CALL Initialise_model_vars (anal_inc, .FALSE.)
CALL Initialise_CVs (dchi0_ip1, .FALSE.)
CALL Initialise_CVs (r_i, .FALSE.)
CALL Initialise_CVs (r_ip1, .FALSE.)
CALL Initialise_CVs (p_i, .FALSE.)
CALL Initialise_CVs (p_ip1, .FALSE.)
CALL Initialise_CVs (pert, .FALSE.)
CALL Initialise_CVT (CVT)


! Main filenames (other files diagnostic files are used within the DA loops)
! Background state (input)
Bg_filename             = TRIM(datadir_Bg)  // '/' // TRIM(Bg_file)
! CVT data (intput)
CVT_filename            = TRIM(datadirCVT)  // '/' // TRIM(CVT_file)
! Obs data (input)
Obs_filename            = TRIM(datadir_Obs) // '/' // TRIM(Obs_file)
! Diagnostics (output)
Scalar_diags_filename   = TRIM(datadirAnal) // '/' // TRIM(diagnostics_file)
! Analysis (output)
Anal_filename           = TRIM(datadirAnal) // '/' // TRIM(anal_file)
! Analysis increment (output)
Analinc_filename        = TRIM(datadirAnal) // '/' // TRIM(analinc_file)

OPEN (12, file=Scalar_diags_filename)
WRITE (12,'(2A7, 9A15)') '#outer', 'inner', 'Jb', 'Jo', 'Je', 'J', '|grad|', 'KE', 'BE', 'EE', 'TE'
WRITE (12,'(2A7, 9A15)') '#-----', '-----', '--', '--', '--', '-', '------', '--', '--', '--', '--'

PRINT*, 'Reading in CVT data ...'
CALL Read_Covs (CVT_filename, CVT,              &
                .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.)
PRINT *, '-- done'

! Read in the background state (-1 means the last item in the input file)
PRINT*, 'Reading in background state ...'
CALL Read_state_2d (Bg_filename, Bg0, dims, -1, .TRUE.)
PRINT*, '-- done'

! Set some commonly-used constants
CALL Set_ht_dep_cons (dims)

! Compute balance and energy diagnostics of the background
CALL Calc_geost (Bg0)
CALL Calc_hydro (Bg0, dims)
CALL Energy (Bg0, dims)

! Initialise the observation structure (a linked list)
NULLIFY(Observations)

! Read-in the observations
PRINT*, 'Reading in observations ...'
CALL Read_Obs (Observations, Obs_filename, dt, timesteps, dt_da, DAtimesteps, maxtime)
PRINT *, 'Model timestep ', dt
PRINT *, 'DA timestep    ', dt_da
PRINT *, '-- done'

!Allocate the reference state array over the DA time steps
ALLOCATE (LSfc(0:DAtimesteps))
DO t = 0, DAtimesteps
  CALL Initialise_model_vars (LSfc(t), .FALSE.)
END DO

! Error checking for hybrid options (pure EnVar and hybrid)
IF (Hybrid_opt == 2 .OR. Hybrid_opt == 3) THEN
  Neffmems   = NEnsMems

  ! Ensure number of effective ensemble members are more than 1
  IF (Neffmems > 1) THEN
    PRINT*, 'Using', Neffmems, 'ensemble error modes'
  ELSE
    PRINT*, 'EnVar enabled but invalid number of ensemble error modes'
    STOP
  END IF
END IF

! Count number of alpha control variable fields required; based on number of ensemble members
! Set conditional to call hybrid subroutines instead
Use_EOTD = .FALSE.
IF (Hybrid_opt == 3) THEN
  Nacv    = Neffmems
  Use_EOTD = .TRUE.

  ! Set weighting factors
  beta_c = SQRT(Cov_WeightC/100)
  beta_e = SQRT(Cov_WeightE/100)

! Temporary code which uses hybrid algorithm but for a special case - full weight to static B
ELSE IF (Hybrid_opt == 2) THEN
  ! Pure EnVar so forcing static weighting to 0
  ! Set weighting factors
  beta_c = 0
  beta_e = SQRT(Cov_WeightE/100)
  Use_EOTD = .TRUE.

! Set Nacv with a dummy value of 1 when EnVar is not enabled
ELSE
  Nacv = 1
END IF

IF (Use_EOTD) THEN
  PRINT *, 'Beta (climatological) is ', beta_c
  PRINT *, 'Beta (ensemble) is ', beta_e
  ! Set up localisation matrix
  CALL Initialise_aCVT(L_alpha)
  CALL SetHorizLoc_alpha (hScale_alpha, dims, L_alpha)
  CALL SetVertLoc_alpha (vScale_alpha, dims, L_alpha)
ELSE
  CALL Initialise_aCVT(L_alpha)
END IF

! Allocate the number of alpha control variables based on above conditions, dummy value of 1
ALLOCATE (chi0_i_alpha(1:Nacv))
ALLOCATE (r_i_alpha(1:Nacv))
ALLOCATE (r_ip1_alpha(1:Nacv))
ALLOCATE (EM_x(1:Nacv))
ALLOCATE (tmp_Uv(1:Nacv))
ALLOCATE (pert_alpha(1:Nacv))
ALLOCATE (chi0_ip1_alpha(1:Nacv))
ALLOCATE (p_i_alpha(1:Nacv))
ALLOCATE (p_ip1_alpha(1:Nacv))
IF (VarLoc .NE. 1 .AND. VarLoc .NE. 2) THEN
  ALLOCATE (r_i_alpha_copy(1:Nacv))
  ALLOCATE (r_ip1_alpha_copy(1:Nacv))
END IF
! Initialise fields, but not used if EnVar is not enabled
DO n = 1, Nacv
  CALL Initialise_model_vars (r_i_alpha(n), .FALSE.)
  CALL Initialise_model_vars (r_ip1_alpha(n), .FALSE.)
  CALL Initialise_model_vars (pert_alpha(n), .FALSE.)
  CALL Initialise_model_vars (chi0_ip1_alpha(n), .FALSE.)
  CALL Initialise_model_vars (p_i_alpha(n), .FALSE.)
  CALL Initialise_model_vars (p_ip1_alpha(n), .FALSE.)
END DO

! Count number of error mode files and ensure it matches number of effective ensemble members
IF (Use_EOTD) THEN
  item = 0
  DO n = 1, Nacv
    ! Check if file exists in directory, files should follow standard template
    WRITE (ABCfilename, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirEM), '/PertABC_Item', n, '.nc'
    INQUIRE(file=ABCfilename, exist=exist)
    IF (exist) THEN
      ! Read in error mode states to a appropriate model state array if file exists
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

! ===== Start the outer loop ================================================
outerloop: DO k = 1, N_outerloops + 1

  PRINT *, '===== Starting outer loop', k


  ! Sort out the reference state and background perts
  PRINT *, 'Sorting out the initializations ...'
  IF (k == 1) THEN
    ! Set the reference state to the background
    LSfc(0) = Bg0
    ! This means that dBg0 and dchiB0 are both zero
    CALL Initialise_model_vars (dBg0, .FALSE.)
    CALL Initialise_CVs (dchiB0, .FALSE.)
    dchiB0zero = .TRUE.
  ELSE
    ! Compute the background as a perturbation
    dBg0 = Bg0
    CALL Subtract_model_vars (dBg0, LSfc(0), .TRUE.)
    ! Transform to control space
    CALL U_trans_inv (LSfc(0), dchiB0, dBg0, CVT, dims)

    dchiB0zero = .FALSE.
  END IF
  PRINT *, '-- done'

  ! In the inner loop, the first guess of the state is the reference
  ! This means that dchi0 is zero
  CALL Initialise_CVs (dchi0_i, .FALSE.)
  dchi0zero = .TRUE.

  ! Need to initialise alpha control variables in extended control vector if EnVar enabled
  ! chi0_i_alpha is an array of length Nacv; each index contains a model_vars structure (used in Schur product)
  IF (Use_EOTD) THEN
    DO n = 1, Nacv
      CALL Initialise_model_vars (chi0_i_alpha(n), .FALSE.)
    END DO
  ! Dummy fields, not eventually used in PenAndGrad
  ELSE
    DO n = 1, Nacv
      CALL Initialise_model_vars (chi0_i_alpha(n), .FALSE.)
      CALL Initialise_model_vars (EM_x(n), .FALSE.)
    END DO
  END IF

  ! Set to run the NL model in the PenAndGrad routine, despite assimilation mode
  ! This is used to run the NL model from the analysis (to find the next background) at the end of the DA
  ForceNLmodel = k > N_outerloops

  ! Compute the initial gradient of the cost function
  PRINT *, '===== Starting initial PenAndGrad routine'
  CALL PenAndGrad (VarType, DAtimesteps, timesteps, dims, t0, Nacv,                   &
                   Observations, dchi0_i, dchi0zero, dchiB0, dchiB0zero,              &
                   CVT, LSfc, Jb, Jo, Je, J, .TRUE., r_i, r_i_alpha,                  &
                   .TRUE., ForceNLmodel, datadirAnal, k, 0, chi0_i_alpha, Use_EOTD, EM_x, L_alpha)
  PRINT *, '-- done'

  IF (k > N_outerloops) THEN
    WRITE (12,'(2I7, 9F15.5)') k, 0, Jb, Jo, Je, J, 99999.0,                           &
                               LSfc(0) % Kinetic_Energy, LSfc(0) % Buoyant_Energy, &
                               LSfc(0) % Elastic_Energy, LSfc(0) % Total_Energy
  ELSE


    ! CALL Write_CV ('r_0_cv.nc', r_i, 14, CVT, dims % longs_v(1:nlongs), dims % full_levs(1:nlevs))


    ! The residual is minus this
    CALL Minus_CVs (r_i)
    magnitude2_r_i = REAL(InnerProdControlSpace(r_i, r_i, ComplexSpace=.TRUE.))

    ! Add magnitude contributions from residual of ensemble portions
    IF (Use_EOTD) THEN
      DO n = 1, Nacv
        CALL Mul_model_cons (r_i_alpha(n), m1, .TRUE.)

        ! Without inter-variable localisation, we treat use the same alpha field for each model variable
        ! or certain variables, so need to account for any repeated inner product computations
        IF (VarLoc == 1) THEN ! Full inter-variable localisation
          magnitude2_r_i_alpha = InnerProdModelSpace(r_i_alpha(n), r_i_alpha(n))
          magnitude2_r_i = magnitude2_r_i + magnitude2_r_i_alpha
        
        ELSE IF (VarLoc == 2) THEN ! Retain all inter-variable covariances
          magnitude2_r_i_alpha = InnerProdModelSpace(r_i_alpha(n), r_i_alpha(n))
          magnitude2_r_i = magnitude2_r_i + magnitude2_r_i_alpha/N_p
        
        ELSE IF (VarLoc == 4) THEN ! Retain all inter-variable covariances except with v
          ! Add components associated with all other variables
          CALL Initialise_model_vars (r_i_alpha_copy(n), .FALSE.)
          CALL Add_model_vars (r_i_alpha_copy(n), r_i_alpha(n), .TRUE.)
          r_i_alpha_copy(n) % v(0:nlongs+1,0:nlevs+1) = 0 
          magnitude2_r_i_alpha = InnerProdModelSpace(r_i_alpha_copy(n), r_i_alpha_copy(n))
          magnitude2_r_i = magnitude2_r_i + magnitude2_r_i_alpha/4.0

          ! Add components associated with v
          CALL Initialise_model_vars (r_i_alpha_copy(n), .FALSE.)
          CALL Add_model_vars (r_i_alpha_copy(n), r_i_alpha(n), .TRUE.)
          r_i_alpha_copy(n) % u(0:nlongs+1,0:nlevs+1) = 0
          r_i_alpha_copy(n) % w(0:nlongs+1,0:nlevs+1) = 0
          r_i_alpha_copy(n) % r(0:nlongs+1,0:nlevs+1) = 0
          r_i_alpha_copy(n) % b(0:nlongs+1,0:nlevs+1) = 0
          magnitude2_r_i_alpha = InnerProdModelSpace(r_i_alpha_copy(n), r_i_alpha_copy(n))
          magnitude2_r_i = magnitude2_r_i + magnitude2_r_i_alpha
 
        END IF
      END DO
    END IF

    ! Store this value as the initial magnitude squared of the residual
    magnitude2_r_0 = magnitude2_r_i

    ! Output initial diagnostics
    WRITE (12,'(2I7, 9F15.5)') k, 0, Jb, Jo, Je, J, SQRT(magnitude2_r_i),              &
                               LSfc(0) % Kinetic_Energy, LSfc(0) % Buoyant_Energy, &
                               LSfc(0) % Elastic_Energy, LSfc(0) % Total_Energy

    ! The initial search direction is the initial residual
    p_i = r_i

    ! Repeat and get initial search direction for each alpha control variable
    IF (Use_EOTD) THEN
      DO n = 1, Nacv
        p_i_alpha(n) = r_i_alpha(n)
      END DO
    END IF

    ! ===== Start the inner loop ================================================
    i = 1
    innerloop: DO

      PRINT *, 'Starting inner loop:', i

      ! Do a line minimization in the search direction
      ! Compute J at a positively perturbed state
      CALL Add_pert_CVs (pert, dchi0_i, p_i, mu)
      ! Repeat for each alpha control variable
      IF (Use_EOTD) THEN
        DO n = 1, Nacv
          CALL Add_pert_model_vars (pert_alpha(n), chi0_i_alpha(n), p_i_alpha(n), mu)
        END DO
      END IF   
      PRINT *, 'Starting PenAndGrad routine for J+'
      CALL PenAndGrad ( VarType, DAtimesteps, timesteps, dims, t0, Nacv,                        &
                        Observations, pert, .FALSE., dchiB0, dchiB0zero,                        &
                        CVT, LSfc, Jb_plus, Jo_plus, Je_plus, J_plus, .FALSE., r_i, r_i_alpha,  &
                        .FALSE., .FALSE., datadirAnal, k, i, pert_alpha, Use_EOTD, EM_x, L_alpha )
      PRINT *, 'Positive pert (Jb+, Jo+, Je+, J+):', Jb_plus, Jo_plus, Je_plus, J_plus


      ! Compute J at a negatively perturbed state
      CALL Add_pert_CVs (pert, dchi0_i, p_i, minus_mu)
      ! Repeat for each alpha control variable
      IF (Use_EOTD) THEN
        DO n = 1, Nacv
          CALL Add_pert_model_vars (pert_alpha(n), chi0_i_alpha(n), p_i_alpha(n), minus_mu)
        END DO
      END IF
      PRINT *, 'Starting PenAndGrad routine for J-'
      CALL PenAndGrad ( VarType, DAtimesteps, timesteps, dims, t0, Nacv,                        &
                        Observations, pert, .FALSE., dchiB0, dchiB0zero,                        &
                        CVT, LSfc, Jb_minus, Jo_minus, Je_minus, J_minus, .FALSE., r_i, r_i_alpha,  &
                        .FALSE., .FALSE., datadirAnal, k, i, pert_alpha, Use_EOTD, EM_x, L_alpha )
      PRINT *, 'Negative pert (Jb-, Jo-, Je-, J-):', Jb_minus, Jo_minus, Je_minus, J_minus

      beta_i = mu * (J_minus - J_plus) / (2.0 * (J_minus + J_plus - 2.0 * J))

      ! The adjustment to the state is beta_i * p_i
      ! What is the size of this adjustment?
      magnitude_adjust_inner = beta_i * SQRT(REAL(InnerProdControlSpace(p_i, p_i, ComplexSpace=.TRUE.)))
      !PRINT *, 'Change in dchi:', k, i, magnitude_adjust_inner

      ! Can add a similar diagnostic for alpha control variables here

      ! Update the state
      CALL Add_pert_CVs (dchi0_ip1, dchi0_i, p_i, beta_i)

      ! Repeat for each alpha control variable
      IF (Use_EOTD) THEN
        DO n = 1, Nacv
          CALL Add_pert_model_vars (chi0_ip1_alpha(n), chi0_i_alpha(n), p_i_alpha(n), beta_i)
        END DO
      END IF
     

      dchi0zero = .FALSE.

      ! Compute a new gradient
      PRINT *, 'Starting main PenAndGrad routine for inner loop', i
      CALL PenAndGrad ( VarType, DAtimesteps, timesteps, dims, t0, Nacv,                      &
                        Observations, dchi0_ip1, dchi0zero, dchiB0, dchiB0zero,               &
                        CVT, LSfc, Jb, Jo, Je, J, .TRUE., r_ip1, r_ip1_alpha,                     &
                        .TRUE., .FALSE., datadirAnal, k, i, chi0_ip1_alpha, Use_EOTD, EM_x, L_alpha ) 

      ! The residual is minus this
      CALL Minus_CVs (r_ip1)
      magnitude2_r_ip1 = REAL(InnerProdControlSpace(r_ip1, r_ip1, ComplexSpace=.TRUE.))

      ! Add magnitude contributions from residual of ensemble portions
      IF (Use_EOTD) THEN
        DO n = 1, Nacv
          CALL Mul_model_cons (r_ip1_alpha(n), m1, .TRUE.)

          ! Without inter-variable localisation, we treat use the same alpha field for each model variable
          ! or certain variables, so need to account for any repeated inner product computations
          IF (VarLoc == 1) THEN ! Full inter-variable localisation
            magnitude2_r_ip1_alpha = InnerProdModelSpace(r_ip1_alpha(n), r_ip1_alpha(n))
            magnitude2_r_ip1 = magnitude2_r_ip1 + magnitude2_r_ip1_alpha
          
          ELSE IF (VarLoc == 2) THEN ! Retain all inter-variable covariances
            magnitude2_r_ip1_alpha = InnerProdModelSpace(r_ip1_alpha(n), r_ip1_alpha(n))
            magnitude2_r_ip1 = magnitude2_r_ip1 + magnitude2_r_ip1_alpha/N_p
         
          ELSE IF (VarLoc == 4) THEN ! Retain all inter-variable covariances except with v
            ! Add components associated with all other variables
            CALL Initialise_model_vars (r_ip1_alpha_copy(n), .FALSE.)
            CALL Add_model_vars (r_ip1_alpha_copy(n), r_ip1_alpha(n), .TRUE.)
            r_ip1_alpha_copy(n) % v(0:nlongs+1,0:nlevs+1) = 0
            magnitude2_r_ip1_alpha = InnerProdModelSpace(r_ip1_alpha_copy(n), r_ip1_alpha_copy(n))
            magnitude2_r_ip1 = magnitude2_r_ip1 + magnitude2_r_ip1_alpha/4.0

            ! Add components associated with v
            CALL Initialise_model_vars (r_ip1_alpha_copy(n), .FALSE.)
            CALL Add_model_vars (r_ip1_alpha_copy(n), r_ip1_alpha(n), .TRUE.)
            r_ip1_alpha_copy(n) % u(0:nlongs+1,0:nlevs+1) = 0
            r_ip1_alpha_copy(n) % w(0:nlongs+1,0:nlevs+1) = 0
            r_ip1_alpha_copy(n) % r(0:nlongs+1,0:nlevs+1) = 0
            r_ip1_alpha_copy(n) % b(0:nlongs+1,0:nlevs+1) = 0
            magnitude2_r_ip1_alpha = InnerProdModelSpace(r_ip1_alpha_copy(n), r_ip1_alpha_copy(n))
            magnitude2_r_ip1 = magnitude2_r_ip1 + magnitude2_r_ip1_alpha

          END IF
        END DO
      END IF

      ! Output diagnostics
      WRITE (12,'(2I7, 9F15.5)') k, i, Jb, Jo, Je, J, SQRT(magnitude2_r_i),         &
                                 0.0, 0.0, 0.0, 0.0

      !PRINT *, 'Next residual:', k, i, SQRT(magnitude2_r_ip1)

      ! Have we converged?
      Converged_inner = ( SQRT(magnitude2_r_ip1 / magnitude2_r_0) < crit_inner )

      IF (.NOT.Converged_inner) THEN
        alpha_i_i = magnitude2_r_ip1 / magnitude2_r_i

        ! Compute the new search direction
        CALL Add_pert_CVs (p_ip1, r_ip1, p_i, alpha_i_i)

        ! Repeat for each alpha control variable
        IF (Use_EOTD) THEN
          DO n = 1, Nacv
            CALL Add_pert_model_vars (p_ip1_alpha(n), r_ip1_alpha(n), p_i_alpha(n), alpha_i_i)
          END DO
        END IF

        ! Increment inner loop counter, and shift variables
        i              = i + 1
        r_i            = r_ip1
        magnitude2_r_i = magnitude2_r_ip1
        p_i            = p_ip1
        dchi0_i        = dchi0_ip1

        ! Repeat for ensemble related variables, these are empty and unchanged if Use_EOTD is .FALSE.
        r_i_alpha      = r_ip1_alpha
        p_i_alpha      = p_ip1_alpha
        chi0_i_alpha   = chi0_ip1_alpha

        IF (i > N_innerloops_max) THEN
          EXIT innerloop
        END IF
      ELSE
        EXIT innerloop
      END IF
    END DO innerloop

    ! Compute the increment in model space
    CALL U_trans (LSfc(0), dchi0_ip1, dx0, CVT, dims)

    ! Add ensemble contributions, localise in model space and perform weighting
    IF (Use_EOTD) THEN       
      ! Weight background contribution
      CALL Initialise_model_vars (dxb0, .FALSE.)
      CALL Add_model_vars(dxb0, dx0, .TRUE.)
      CALL Mul_model_cons(dxb0, beta_c, .TRUE.)

      ! Weight ensemble contribution
      CALL Initialise_model_vars (dxe0, .FALSE.)
      DO n = 1, Nacv
        CALL Initialise_model_vars (tmp_Uv(n), .FALSE.)
        ! Localisation of alpha control variables
        ! Note that division factor of N_p (or otherwise) is not required because 
        ! we do not repeat this computation N_p times (and sum)
        ! unlike in true matrix form of U_alpha which has N_g x N_p cols
        CALL U_trans_alpha (LSfc(0), chi0_ip1_alpha(n), tmp_Uv(n), CVT, dims, L_alpha)
        
        ! Schur product multiplication between error modes and alpha control variables in model space
        CALL Mul_model_vars(tmp_Uv(n), EM_x(n), .TRUE.)
        CALL Add_model_vars(dxe0, tmp_Uv(n), .TRUE.)
      END DO

      CALL Mul_model_cons(dxe0, beta_e, .TRUE.)
    
      ! Increment is sum of weighted contributions
      CALL Initialise_model_vars (dx0, .FALSE.)
      CALL Add_model_vars(dx0, dxb0, .TRUE.)
      CALL Add_model_vars(dx0, dxe0, .TRUE.)
    END IF  

    ! Add to reference state
    CALL Add_model_vars (LSfc(0), dx0, .TRUE.)

  END IF

END DO outerloop


! Compute the overall analysis increment
anal_inc = LSfc(0)
CALL Subtract_model_vars (anal_inc, Bg0, .FALSE.)

! Output the analysis and analysis increment
PRINT*, 'Writing the analysis ...'
CALL Write_state_2d (Anal_filename, LSfc(0), dims, 1, 0, 0, .TRUE.)
PRINT*, '  -- done'
PRINT*, 'Writing the analysis increment ...'
CALL Write_state_2d (Analinc_filename, anal_inc, dims, 1, 0, 0, .TRUE.)
PRINT*, '  -- done'

CLOSE (12)

! Deallocate
PRINT*, 'Deallocating ...'
CALL DeAllocate_Obs (Observations)
DO t = 0, DAtimesteps
  CALL Deallocate_model_vars (LSfc(t))
END DO
DEALLOCATE (LSfc)
DEALLOCATE (chi0_i_alpha)
DEALLOCATE (r_i_alpha)
DEALLOCATE (r_ip1_alpha)
DEALLOCATE (EM_x)
DEALLOCATE (tmp_Uv)
DEALLOCATE (pert_alpha)
DEALLOCATE (chi0_ip1_alpha)
DEALLOCATE (p_i_alpha)
DEALLOCATE (p_ip1_alpha)
IF (VarLoc .NE. 1 .AND. VarLoc .NE. 2) THEN
  DEALLOCATE (r_i_alpha_copy)
  DEALLOCATE (r_ip1_alpha_copy)
END IF
CALL Deallocate_dims (dims)
CALL Deallocate_model_vars (Bg0)
CALL Deallocate_model_vars (dBg0)
CALL Deallocate_model_vars (dx0)
CALL Deallocate_model_vars (anal_inc)
CALL Deallocate_CVs (dchiB0)
CALL Deallocate_CVs (dchi0_i)
CALL Deallocate_CVs (dchi0_ip1)
CALL Deallocate_CVs (r_i)
CALL Deallocate_CVs (r_ip1)
CALL Deallocate_CVs (p_i)
CALL Deallocate_CVs (p_ip1)
CALL Deallocate_CVs (pert)
CALL Deallocate_CVT (CVT)
CALL Deallocate_aCVT (L_alpha)
IF (ALLOCATED(fft_wsave_x)) DEALLOCATE(fft_wsave_x)
IF (ALLOCATED(fft_work_x)) DEALLOCATE(fft_work_x)

PRINT*, '  -- done'

END PROGRAM Master_Assimilate
