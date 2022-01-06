SUBROUTINE PenAndGrad ( VarType, dasteps, steps, dims, t0, Nacv,                        &
                        Obs, dchi0, dchi0zero, dchib0, dchib0zero,                      &
                        CVTdata, LSfc, Jb, Jo, Je, J, compute_grad, grad0, grad0_alpha, &
                        outputdiags_flag, forceNLmodel, outputdir, outerloop,           &
                        innerloop, chi0_alpha, Use_EOTD, EM_x, L_alpha )

!*********************************************************************************
!*                                                                               *
!*  Compute penalty and gradient of the cost function                            *
!*  The gradient is in control space                                             *
!*                                                                               *
!*  INPUTS                                                                       *
!*  VarType                - 3=3DVar, 35=3D-FGAT, 4=4DVar                        *
!*  dasteps                - number of time states for the da (excluding 0)      *
!*  steps                  - number of model steps                               *
!*  dims                   - dimension data                                      *
!*  t0                     - time corresponding to time step 0 (seconds)         *
!*  Obs                    - pointer to the start of the observation linked list *
!*  dchi0                  - control variable at start of window                 *
!*  dchi0zero              - set if dchi0=0                                      *
!*  dchib0                 - background increment in control space               *
!*  dchib0zero             - set if dchib0=0                                     *
!*  CVTdata                - data that defines the CVT                           *
!*  INPUT/OUTPUT                                                                 *
!*  LSfc                   - Linearization state trajectory                      *
!*  OUTPUTS                                                                      *
!*  Jb                     - Background part of cost function                    *
!*  Jo                     - Observation part of cost function                   *
!*  Je                     - Ensemble part of cost function (0 if not Use_EOTD)  *
!*  J                      - Total cost function                                 *
!*  compute_grad           - Set to compute the gradient (otherwise only cost    *
!*                           function is found, leaving grad0 unchanged)         *
!*  grad0                  - Gradient of cost function in control space          *
!*  grad0_alpha            - Gradient of cost function in model space for        *
!*                           the alpha control variables                         *
!*  INPUTS                                                                       *
!*  outputdiags_flag       - Set to output various optional diagnostics          *
!*  forceNLmodel           - Set to force run of NL model in ref state           *
!*  outputdir              - Output directory                                    *
!*  outerloop              - outer loop index (for diagnostic file naming)       *
!*  innerloop              - inner loop index (for diagnostic file naming)       *
!*  Use_EOTD               - Logical which controls EnVar subroutines            *
!*  INPUTS                                                                       *
!*  Nacv                   - Number of alpha control variables (set to dummy     *
!*                           value of 1 if Hybrid_opt is not 2 or 3              *
!*  chi0_alpha             - alpha control variable at start of window           *
!*  EM_x                   - array containing error modes in model space         *
!*  L_alpha                - data that defines the localisation matrix           *
!*                                                                               *
!*   R. Bannister, 1.4da        28-02-2018                                       *
!*   J. Lee      , 1.5da_hybrid 06-05-2021 - added hybrid DA components          *
!*                                                                               *
!*********************************************************************************


USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  dims_type,              &
  Hybrid_opt,             &
  nlongs,                 &
  nlevs,                  &
  Obs_type,               &
  CV_type,                &
  CVT_type,               &
  aCVT_type,              &
  dt, dt_da,              &
  Cov_WeightE,            &
  Cov_WeightC,            &
  InterVarLoc             

IMPLICIT NONE

INCLUDE "InnerProduct_obsself.interface"
INCLUDE "InnerProdControlSpace.interface"
INCLUDE "InnerProdModelSpace.interface"
INCLUDE "ModelObservations.interface"
INCLUDE "ModelObservations_ZeroPert.interface"
INCLUDE "ModelObservations_linear.interface"
INCLUDE "ModelObservations_adj.interface"
INCLUDE "Write_Obs.interface"

! Parameters
!-----------
INTEGER,                 INTENT(IN)    :: Vartype
INTEGER,                 INTENT(IN)    :: dasteps
INTEGER,                 INTENT(IN)    :: steps
TYPE(dims_type),         INTENT(IN)    :: dims
INTEGER,                 INTENT(IN)    :: t0
TYPE(Obs_type), POINTER, INTENT(IN)    :: Obs
TYPE(CV_type),           INTENT(IN)    :: dchi0
LOGICAL,                 INTENT(IN)    :: dchi0zero
TYPE(CV_type),           INTENT(IN)    :: dchib0
LOGICAL,                 INTENT(IN)    :: dchib0zero
TYPE(CVT_type),          INTENT(IN)    :: CVTdata
TYPE(aCVT_type),         INTENT(IN)    :: L_alpha
TYPE(ABC_type),          INTENT(INOUT) :: LSfc(0:steps)
REAL(ZREAL8),            INTENT(OUT)   :: Jb
REAL(ZREAL8),            INTENT(OUT)   :: Jo
REAL(ZREAL8),            INTENT(OUT)   :: Je
REAL(ZREAL8),            INTENT(OUT)   :: J
LOGICAL,                 INTENT(IN)    :: compute_grad
TYPE(CV_type),           INTENT(INOUT) :: grad0
LOGICAL,                 INTENT(IN)    :: outputdiags_flag
LOGICAL,                 INTENT(IN)    :: forceNLmodel
CHARACTER(LEN=100),      INTENT(IN)    :: outputdir
INTEGER,                 INTENT(IN)    :: innerloop
INTEGER,                 INTENT(IN)    :: outerloop

LOGICAL,                 INTENT(IN)    :: Use_EOTD
INTEGER,                 INTENT(IN)    :: Nacv
TYPE(ABC_type),          INTENT(IN)    :: chi0_alpha(1:Nacv)
TYPE(ABC_type),          INTENT(INOUT) :: grad0_alpha(1:Nacv)
TYPE(ABC_type),          INTENT(IN)    :: EM_x(1:Nacv)

! Local variables
! ---------------
TYPE(ABC_type), ALLOCATABLE            :: deltax(:)
TYPE(ABC_type), ALLOCATABLE            :: deltax_hat(:)
TYPE(ABC_type), ALLOCATABLE            :: deltax_b(:)
TYPE(ABC_type), ALLOCATABLE            :: deltax_e(:)
INTEGER                                :: t, ModelStepsPerDAStep, n
COMPLEX(ZREAL8)                        :: Jb_complex
TYPE(CV_type)                          :: diffcv
CHARACTER(LEN=320)                     :: Diag_filename
TYPE(ABC_type), ALLOCATABLE            :: diffcv_alpha(:), tmp_Uv(:), EM_x_copy(:), tmp_alpha(:)
REAL(ZREAL8)                           :: beta_c, beta_e


! Prepare memory
ALLOCATE (deltax(0:dasteps))
ALLOCATE (deltax_b(0:dasteps))
ALLOCATE (deltax_e(0:dasteps))
DO t = 0, dasteps
  CALL Initialise_model_vars (deltax(t), .FALSE.)
END DO
CALL Initialise_CVs (diffcv, .FALSE.)
ALLOCATE (diffcv_alpha(1:Nacv))
ALLOCATE (tmp_Uv(1:Nacv))
ALLOCATE (EM_x_copy(1:Nacv))
ALLOCATE (tmp_alpha(1:Nacv))
DO n = 1, Nacv
  CALL Initialise_model_vars (diffcv_alpha(n), .FALSE.)
  CALL Initialise_model_vars (tmp_Uv(n), .FALSE.)
  CALL Initialise_model_vars (EM_x_copy(n), .FALSE.)
  CALL Initialise_model_vars (tmp_alpha(n), .FALSE.)
END DO


ModelStepsPerDAStep = INT(dt_da/dt)

IF (Use_EOTD) THEN
  beta_c = SQRT(Cov_WeightC/100)
  beta_e = SQRT(Cov_WeightE/100)
  IF (Hybrid_opt == 2) THEN
    beta_c = 0
  END IF
END IF

! ===================================
! --- FORWARD PART OF CALCULATION ---
! ===================================
IF (dchi0zero) THEN
  ! This means that we are at the start of a new outer loop
  PRINT *, 'dchi is zero'
  IF (outputdiags_flag) THEN
    PRINT *, 'Computing LS balance and energy diagnostics for t=0'
    CALL Calc_geost (LSfc(0))
    CALL Calc_hydro (LSfc(0), dims)
    CALL Energy (LSfc(0), dims)
  END IF
  IF ((Vartype == 35) .OR. (Vartype == 4) .OR. forceNLmodel) THEN
    ! 3DFGAT or 4DVar, so compute the non-linear model trajectory
    PRINT *, '3DFGAT, 4DVar, or just propagating analysis throughout time window.'
    PRINT *, 'About to run nonlinear model for ', steps, ' time steps'
    CALL ABC_NL_ModelDriver_DA (LSfc(0:dasteps), dims, steps, dasteps)
    PRINT *, ' -- done'
    IF (outputdiags_flag) THEN
      DO t = 1, dasteps
        PRINT *, 'Computing LS balance diagnostics for t=', t
        CALL Calc_geost (LSfc(t))
        CALL Calc_hydro (LSfc(t), dims)
        CALL Energy (LSfc(t), dims)
      END DO
    END IF
  ELSE

    ! 3DVar, so just copy t=0 state to all times
    PRINT *, '3DVar.  About to copy states'
    DO t = 1, dasteps
      LSfc(t) = LSfc(0)
    END DO
    PRINT *, ' -- done'
  END IF

  IF (outputdiags_flag) THEN
    ! Output the LS state throughout the window
    DO t = 0, dasteps
      WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/LS_Oloop', outerloop, '_Iloop', innerloop, '.nc'
      ! PRINT *, 'Time', t, ':  Outputting LS data to ', TRIM(Diag_filename)
      CALL Write_state_2d (TRIM(Diag_filename), LSfc(t), dims, dasteps+1, t, ModelStepsPerDAStep, .TRUE.)
    END DO
  END IF


  ! Compute the reference model observations
  PRINT *, 'Computing model observations ...'
  CALL ModelObservations (dasteps, LSfc(0:dasteps), dims, dt_da, t0, &
                          Obs, .FALSE.)
  PRINT *, ' -- done'

  ! Set the remaining observation structure elements for this ref state and zero pert
  PRINT *, 'Setting model perts to zero ...'
  CALL ModelObservations_ZeroPert (Obs)
  PRINT *, ' -- done'


  IF (outputdiags_flag) THEN
    ! Output the observation structure
    WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/Obs_', outerloop, '_Iloop', innerloop, '.dat'
    PRINT *, 'Outputting obs structure ', TRIM(Diag_filename)
    CALL Write_Obs ( Diag_filename, 0, dt, 0, dt_da, 0, Obs )
  END IF

ELSE

  ! The reference values are already known
  PRINT *, 'dchi is non-zero'

  ! Convert the control variable pert to model space
  CALL U_trans (LSfc(0), dchi0, deltax(0), CVTdata, dims)

  ! Add EOTD contributions, localise in model space and perform weighting
  IF (Use_EOTD) THEN
    ! Weight contribution from background
    ! PRINT *, 'Weighting background'
    DO t = 0, dasteps
      CALL Initialise_model_vars (deltax_b(t), .FALSE.)
    END DO
    CALL Add_model_vars(deltax_b(0), deltax(0), .TRUE.)
    CALL Mul_model_cons(deltax_b(0), beta_c, .TRUE.)

    ! Weight contribution from ensemble
    ! PRINT *, 'Weighting ensemble'
    DO t = 0, dasteps
      CALL Initialise_model_vars (deltax_e(t), .FALSE.)
    END DO
    ! Add linear combination of error modes weighted by alpha control variable

    DO n = 1, Nacv
      ! Localisation of alpha control variables in model space
      CALL U_trans_alpha (LSfc(0), chi0_alpha(n), tmp_Uv(n), CVTdata, dims, L_alpha)

      ! Schur product multiplication between error modes and alpha control variables in model space
      CALL Mul_model_vars(tmp_Uv(n), EM_x(n), .TRUE.)
      CALL Add_model_vars(deltax_e(0), tmp_Uv(n), .TRUE.)
    END DO

    CALL Mul_model_cons(deltax_e(0), beta_e, .TRUE.)
    
    ! Increment is sum of weighted contributions
    CALL Initialise_model_vars (deltax(0), .FALSE.)

    CALL Add_model_vars(deltax(0), deltax_b(0), .TRUE.)
    CALL Add_model_vars(deltax(0), deltax_e(0), .TRUE.)

  END IF  

  IF (outputdiags_flag) THEN
    PRINT *, 'Computing pert balance diagnostics for t=0'
    CALL Calc_geost (deltax(0))
    CALL Calc_hydro (deltax(0), dims)
  END IF

  ! Propagate the initial perturbation throughout the time window
  IF ((VarType == 3) .OR. (VarType == 35)) THEN
    ! 3DVar or 3DFGAT, so initial perturbation does not change
    DO t = 1, dasteps
      deltax(t) = deltax(0)
    END DO
  ELSE
    ! 4DVar, so propagate perturbation
    ! Need to write tangent linear model
    ! Need to call balance diagnostics as above
    PRINT*, 'Error - 4DVar not yet implemented.'
    STOP
  END IF

!  IF (outputdiags_flag) THEN
!    ! Output the perturbation throughout the window
!    DO t = 0, dasteps
!      WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/dx_Oloop', outerloop, &
!                            '_Iloop', innerloop, '.nc'
!      PRINT *, 'Time', t, ':  Outputting pert data to ', TRIM(Diag_filename)
!      CALL Write_state_2d (TRIM(Diag_filename), deltax(t), dims, dasteps+1, t, ModelStepsPerDAStep, .TRUE.)
!    END DO
!  END IF

  ! Act with the linear observation operator
  CALL ModelObservations_linear (dasteps, LSfc(0:dasteps), deltax(0:dasteps), dims, &
                                 dt_da, t0, Obs )
    

END IF




! To do with the background term
! Compute the difference dchi0 - dchib0
IF (.NOT.dchi0zero) THEN
  CALL Add_CVs(diffcv, dchi0)
END IF
IF (.NOT.dchib0zero) THEN
  CALL Subtract_CVs(diffcv, dchib0)
END IF

! To do with ensemble terms
! Set diffcv_alpha(:) = chi0_alpha(:) for easier referencing
IF (.NOT.dchi0zero) THEN
  IF (Use_EOTD) THEN
    DO n = 1, Nacv
      CALL Add_model_vars(diffcv_alpha(n), chi0_alpha(n), .TRUE.)
    END DO
  END IF
END IF

!IF (outputdiags_flag) THEN
!  ! Output the contribution to the gradient from the background term
!  WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/GradBg_Oloop', outerloop, &
!                        '_Iloop', innerloop, '.nc'
!  PRINT *, 'Outputting bg grad data to ', TRIM(Diag_filename)
!  CALL Write_CV (TRIM(Diag_filename), diffcv, 1, CVTdata, dims % longs_v(1:nlongs), dims % full_levs(1:nlevs))
!END IF


IF (compute_grad) THEN

  ! ===================================
  ! --- ADJOINT PART OF CALCULATION ---
  ! ===================================

  ! Set the adjoint states to zero
  ALLOCATE (deltax_hat(0:dasteps))
  DO t = 0, dasteps
    CALL Initialise_model_vars (deltax_hat(t), .FALSE.)
  END DO

  ! Perform the adjoint of the model observations
  CALL ModelObservations_adj (dasteps, LSfc(0:dasteps), deltax_hat(0:dasteps), dims,  &
                              dt_da, t0, Obs, 1)
  ! At this stage the deltax_hat states are the Delta vectors in the documentation

  IF (outputdiags_flag .AND. (outerloop == 1) .AND. (innerloop == 1)) THEN
    ! Output the Delta throughout the window
    DO t = 0, dasteps
      WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/Delta_Oloop', outerloop, &
                            '_Iloop', innerloop, '.nc'
      ! PRINT *, 'Time', t, ':  Outputting Delta data to ', TRIM(Diag_filename)
      CALL Write_state_2d (TRIM(Diag_filename), deltax_hat(t), dims, dasteps+1, t, ModelStepsPerDAStep, .TRUE.)
    END DO
  END IF

  ! Integrate the adjoint states backwards in time
  IF ((VarType == 3) .OR. (VarType == 35)) THEN
    ! 3DVar or 3DFGAT, so no time evolution of the adjoints
    DO t = dasteps-1, 0, -1
      CALL Add_model_vars(deltax_hat(t), deltax_hat(t+1), .TRUE.)
    END DO
  ELSE
    ! 4DVar, so propagate perturbation
    PRINT*, 'Error - 4DVar not yet implemented.'
    STOP
  END IF

  IF (outputdiags_flag .AND. (outerloop == 1) .AND. (innerloop == 1)) THEN
    ! Output the gradient of Jo throughout the window
    DO t = 0, dasteps
      WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/GradJo_Oloop', outerloop, &
                            '_Iloop', innerloop, '.nc'
      ! PRINT *, 'Time', t, ':  Outputting Delta data to ', TRIM(Diag_filename)
      CALL Write_state_2d (TRIM(Diag_filename), deltax_hat(t), dims, dasteps+1, t, ModelStepsPerDAStep, .TRUE.)
    END DO
  END IF
  
  ! Operate with the adjoint of the CVT to obtain the gradient (observation contribution) in control space
  CALL U_trans_adj (LSfc(0), grad0, deltax_hat(0), CVTdata, dims)

  ! Apply weighting
  IF (Use_EOTD) THEN
    CALL Mul_CV_cons (grad0, beta_c)
  END IF

  ! Add the background contribution to the gradient
  CALL Add_CVs(grad0, diffcv)

  ! Compute the ensemble part of the gradient
  IF (Use_EOTD) THEN

    ! Initialise empty gradient vector
    DO n = 1, Nacv
      CALL Initialise_model_vars (grad0_alpha(n), .FALSE.)
      CALL Initialise_model_vars (tmp_alpha(n), .FALSE.)

      ! Copy error modes data so original does not get modified
      CALL Add_model_vars (EM_x_copy(n), EM_x(n), .TRUE.)

      ! Schur product with same gradient contribution from observations for each n
      CALL Mul_model_vars( EM_x_copy(n), deltax_hat(0), .TRUE. )
   
      ! Weight ensemble contribution to gradient
      ! Apply adjoint of localisation on error modes
      CALL U_trans_alpha_adj( LSfc(0), tmp_alpha(n), EM_x_copy(n), CVTdata, dims, L_alpha )

      IF (.NOT. InterVarLoc) THEN
        ! Same alpha fields applied on each variable from each variable
        CALL Apply_alpha_model_vars (grad0_alpha(n), tmp_alpha(n))
      ELSE
        grad0_alpha(n) = tmp_alpha(n)
      END IF

      ! Apply weighting
      CALL Mul_model_cons (grad0_alpha(n), beta_e, .TRUE.)

      ! Add gradient contributions from both ensemble parts
      CALL Add_model_vars(grad0_alpha(n), diffcv_alpha(n), .TRUE.)
    END DO

  ELSE
    DO n = 1, Nacv
      ! Initialise empty gradient vector which is not changed
      CALL Initialise_model_vars (grad0_alpha(n), .FALSE.)
    END DO
  END IF

  ! Tidy up

  DO t = 0, dasteps
    CALL Deallocate_model_vars (deltax_hat(t))
  END DO
  DEALLOCATE (deltax_hat)
PRINT *, 'DONE'

END IF



! ===================================
! --- VALUES OF COST FUNCTION ---
! ===================================

! The background part of the cost function
Jb_complex = InnerProdControlSpace (diffcv, diffcv, ComplexSpace=.TRUE., ignore_halos=.TRUE.) / 2.0
PRINT *, 'Jb = ', Jb_complex
Jb         = REAL(Jb_complex)
PRINT *, 'Jb = ', Jb

! The observation part of the cost function
Jo         = InnerProduct_obsself (Obs, .TRUE., 2) / 2.0
PRINT *, 'Jo = ', Jo

! The ensemble part of the cost function
Je = 0
IF (Use_EOTD) THEN
  DO n = 1, Nacv
    IF (.NOT. InterVarLoc) THEN
      Je = Je + InnerProdModelSpace (diffcv_alpha(n), diffcv_alpha(n), ignore_halos=.TRUE.) / (2.0*5.0)
    ELSE
      Je = Je + InnerProdModelSpace (diffcv_alpha(n), diffcv_alpha(n), ignore_halos=.TRUE.) / 2.0
    END IF
  END DO
END IF
PRINT *, 'Je = ', Je

! The total cost function
J          = Jb + Jo + Je
PRINT *, 'J = ', J



! Tidy up
PRINT *, 'Tidying up for deltax'
DO t = 0, dasteps
PRINT *, t
  CALL Deallocate_model_vars (deltax(t))
END DO
PRINT *, 'Final tidy up for deltax'
DEALLOCATE (deltax)
PRINT *, 'Final tidy up for diffcv'
CALL Deallocate_CVs (diffcv)
PRINT *, 'Final tidy up for ensemble-related arrays'
DO n = 1, Nacv
  CALL Deallocate_model_vars (diffcv_alpha(n))
  CALL Deallocate_model_vars (tmp_Uv(n))
  CALL Deallocate_model_vars (EM_x_copy(n))
  CALL Deallocate_model_vars (tmp_alpha(n))
END DO
PRINT *, 'DONE'



END SUBROUTINE PenAndGrad
