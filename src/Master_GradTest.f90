PROGRAM Master_GradTest

!*****************************************************
!*   Master code run gradient test                   *
!*                                                   *
!*   Ross Bannister, r.n.bannister@reading.ac.uk     *
!*   25/04/18                                        *
!*   Joshua Lee,     joshua_lee@nea.gov.sg           *
!*   06/06/21                                        *
!*       - added portion for when hybrid is enabled  *
!*   Based on Master_Assimilate.f90                  *
!*                                                   *
!*   Phi = J(chi + alpha h) - J(chi)                 *
!*         -------------------------                 *
!*             alpha h^T grad                        *
!*                                                   *
!*   We choose h = grad / |grad|                     *
!*                                                   *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    nlongs, nlevs,               &
    Vartype,                     &
    datadirCVT,                  &
    Hybrid_opt,                  &
    NEnsMems,                    &
    Cov_WeightE,                 &
    Cov_WeightC,                 &
    hScale_alpha,                &
    vScale_alpha,                &
    datadirEM,                   &
    CVT_file,                    &
    datadirTestDA,               &
    datadir_Bg,                  &
    Bg_file,                     &
    datadir_Obs,                 &
    Obs_file,                    &
    dt, dt_da, t0,               &
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
INCLUDE "Read_Obs.interface"
INCLUDE "DeAllocate_Obs.interface"
INCLUDE "PenAndGrad.interface"


! Declare variables
!==========================
CHARACTER(LEN=320)          :: Bg_filename, CVT_filename, Obs_filename
CHARACTER(LEN=320)          :: Scalar_diags_filename, ABCfilename
TYPE(ABC_type)              :: Bg0
TYPE(ABC_type), ALLOCATABLE :: LSfc(:)
TYPE(dims_type)             :: dims
TYPE(CV_type)               :: dchiB0, chi0, chi0_pert, grad_chi, h
TYPE(Obs_type), POINTER     :: Observations
INTEGER                     :: timesteps, DAtimesteps, maxtime, power, mul, t
INTEGER                     :: Neffmems, Nacv, n, item
LOGICAL                     :: Use_EOTD, exist
TYPE(CVT_type)              :: CVT
TYPE(aCVT_type)             :: L_alpha
REAL(ZREAL8)                :: Jb, Jo, Je, J, J0, modgrad, Phi, alpha
REAL(ZREAL8)                :: multiply(1:4), beta_c, beta_e, denom
TYPE(ABC_type), ALLOCATABLE :: chi0_alpha(:), chi0_pert_alpha(:), grad_chi_alpha(:)
TYPE(ABC_type), ALLOCATABLE :: EM_x(:), h_alpha(:)



PRINT*, '*************************************************************************'
PRINT*, 'Running Master_GradTest'
PRINT*, '*************************************************************************'


! Read namelist
CALL SetOptions

CALL Initialise_dims (dims)
CALL Initialise_model_vars (Bg0, .FALSE.)
CALL Initialise_CVs (dchiB0, .FALSE.)
CALL Initialise_CVs (chi0, .FALSE.)
CALL Initialise_CVs (chi0_pert, .FALSE.)
CALL Initialise_CVs (grad_chi, .FALSE.)
CALL Initialise_CVs (h, .FALSE.)
CALL Initialise_CVT (CVT)

! Main filenames (other files diagnostic files are used within the DA loops)
! Background state (input)
Bg_filename             = TRIM(datadir_Bg)  // '/' // TRIM(Bg_file)
! CVT data (intput)
CVT_filename            = TRIM(datadirCVT)  // '/' // TRIM(CVT_file)
! Obs data (input)
Obs_filename            = TRIM(datadir_Obs) // '/' // TRIM(Obs_file)
! Diagnostics (output)
Scalar_diags_filename   = TRIM(datadirTestDA) // '/' // TRIM(diagnostics_file)

OPEN (12, file=Scalar_diags_filename)
WRITE (12,'(3A20)') '#alpha', 'Phi', '1-Phi'
WRITE (12,'(3A20)') '#-----', '-----', '-----'

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
  ! These are computed directly in PenAndGrad, but included here for diagnostic purposes
  beta_c = SQRT(Cov_WeightC/100)
  beta_e = SQRT(Cov_WeightE/100)

! Temporary code which uses hybrid algorithm but for a special case - full weight to static B
ELSE IF (Hybrid_opt == 2) THEN
  ! Pure EnVar so forcing static weighting to 0
  ! Set weighting factors
  beta_c = 0
  beta_e = SQRT(Cov_WeightE/100)

! Set Nacv with a dummy value of 1 when hybrid is not enabled
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
ALLOCATE (chi0_alpha(1:Nacv))
ALLOCATE (chi0_pert_alpha(1:Nacv))
ALLOCATE (grad_chi_alpha(1:Nacv))
ALLOCATE (EM_x(1:Nacv))
ALLOCATE (h_alpha(1:Nacv))

DO n = 1, Nacv
  CALL Initialise_model_vars (grad_chi_alpha(n), .FALSE.)
  CALL Initialise_model_vars (h_alpha(n), .FALSE.)
  CALL Initialise_model_vars (chi0_pert_alpha(n), .FALSE.)
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

! Set the reference state to the background
LSfc(0) = Bg0

! Choose a random state in control space (this is chi as in the header notes of this code)
CALL Initialise_CVs (chi0, .TRUE.)

! Also initialise random states of ensemble portions of the control vector
! This is not used if Use_EOTD is .FALSE.
DO n = 1, Nacv
  CALL Initialise_model_vars (chi0_alpha(n), .TRUE.)
END DO

! Choose to be at the background state (dchiB0=0)
CALL Initialise_CVs (dchiB0, .FALSE.)

! Compute the gradient, grad_chi (this is grad in the header notes of this code)
CALL PenAndGrad (VarType, DAtimesteps, timesteps, dims, t0, Nacv,                    &
                 Observations, chi0, .FALSE., dchiB0, .TRUE.,                        &
                 CVT, LSfc, Jb, Jo, Je, J, .TRUE., grad_chi, grad_chi_alpha,         &
                 .FALSE., .FALSE., datadirTestDA, 0, 0, chi0_alpha, Use_EOTD, EM_x, L_alpha)

! Store the value of J at this state J0
J0 = J

! Choose h (in the header notes of this code) to be grad / |grad|
! Compute the size of the gradient, |grad|
modgrad = SUM(grad_chi % v1(1:nlongs,1:nlevs) * grad_chi % v1(1:nlongs,1:nlevs)) +   &
          SUM(grad_chi % v2(1:nlongs,1:nlevs) * grad_chi % v2(1:nlongs,1:nlevs)) +   &
          SUM(grad_chi % v3(1:nlongs,1:nlevs) * grad_chi % v3(1:nlongs,1:nlevs)) +   &
          SUM(grad_chi % v4(1:nlongs,1:nlevs) * grad_chi % v4(1:nlongs,1:nlevs)) +   &
          SUM(grad_chi % v5(1:nlongs,1:nlevs) * grad_chi % v5(1:nlongs,1:nlevs))

! Add size of gradient from ensemble portions |grad_alpha|, for computing h_alpha
IF (Use_EOTD) THEN
  DO n = 1, Nacv
    modgrad = modgrad + &
          SUM(grad_chi_alpha(n) % u(1:nlongs,1:nlevs) * grad_chi_alpha(n) % u(1:nlongs,1:nlevs)) + &
          SUM(grad_chi_alpha(n) % v(1:nlongs,1:nlevs) * grad_chi_alpha(n) % v(1:nlongs,1:nlevs)) + &
          SUM(grad_chi_alpha(n) % w(1:nlongs,1:nlevs) * grad_chi_alpha(n) % w(1:nlongs,1:nlevs)) + &
          SUM(grad_chi_alpha(n) % r(1:nlongs,1:nlevs) * grad_chi_alpha(n) % r(1:nlongs,1:nlevs)) + &
          SUM(grad_chi_alpha(n) % b(1:nlongs,1:nlevs) * grad_chi_alpha(n) % b(1:nlongs,1:nlevs))
  END DO
END IF

modgrad = SQRT(modgrad)

! Divide grad_chi by modgrad to produce h
h % v1(1:nlongs,1:nlevs) = grad_chi % v1(1:nlongs,1:nlevs) / modgrad
h % v2(1:nlongs,1:nlevs) = grad_chi % v2(1:nlongs,1:nlevs) / modgrad
h % v3(1:nlongs,1:nlevs) = grad_chi % v3(1:nlongs,1:nlevs) / modgrad
h % v4(1:nlongs,1:nlevs) = grad_chi % v4(1:nlongs,1:nlevs) / modgrad
h % v5(1:nlongs,1:nlevs) = grad_chi % v5(1:nlongs,1:nlevs) / modgrad

! Repeat for ensemble portions (grad_chi_alpha)
IF (Use_EOTD) THEN
  DO n = 1, Nacv
    h_alpha(n) % u(1:nlongs,1:nlevs) = grad_chi_alpha(n) % u(1:nlongs,1:nlevs) / modgrad
    h_alpha(n) % v(1:nlongs,1:nlevs) = grad_chi_alpha(n) % v(1:nlongs,1:nlevs) / modgrad
    h_alpha(n) % w(1:nlongs,1:nlevs) = grad_chi_alpha(n) % w(1:nlongs,1:nlevs) / modgrad
    h_alpha(n) % r(1:nlongs,1:nlevs) = grad_chi_alpha(n) % r(1:nlongs,1:nlevs) / modgrad
    h_alpha(n) % b(1:nlongs,1:nlevs) = grad_chi_alpha(n) % b(1:nlongs,1:nlevs) / modgrad
  END DO
END IF

multiply(1) = 8.0
multiply(2) = 5.0
multiply(3) = 2.0
multiply(4) = 1.0

! ===== Loop over different alphas ================================================
DO power = 1, -16, -1
  DO mul = 1, 4
    alpha = multiply(mul) * 10.0 ** REAL(power)
    PRINT *, 'Running grad test with alpha = ', alpha

    ! Add alpha * r_i to dchi0_i
    CALL Add_pert_CVs(chi0_pert, chi0, h, alpha)

    ! Repeat for ensemble portions
    IF (Use_EOTD) THEN
      DO n = 1, Nacv
        CALL Add_pert_model_vars(chi0_pert_alpha(n), chi0_alpha(n), h_alpha(n), alpha)
      END DO     
    END IF

    ! Compute the perturbed value of J
    CALL PenAndGrad (VarType, DAtimesteps, timesteps, dims, t0, Nacv,                &
                     Observations, chi0_pert, .FALSE., dchiB0, .TRUE.,               &
                     CVT, LSfc, Jb, Jo, Je, J, .FALSE., grad_chi, grad_chi_alpha,    &
                     .FALSE., .FALSE., datadirTestDA, 0, 0, chi0_pert_alpha, Use_EOTD, EM_x, L_alpha)

    ! Compute and output the required diagnostic
    Phi = (J - J0) / (alpha * modgrad)
    WRITE (12, '(3E20.10)') alpha, Phi, 1.0-Phi

  END DO
END DO


CLOSE (12)

! Deallocate
PRINT*, 'Deallocating ...'
CALL DeAllocate_Obs (Observations)
DO t = 0, DAtimesteps
  CALL Deallocate_model_vars (LSfc(t), .FALSE.)
END DO
DEALLOCATE (LSfc)
DEALLOCATE (chi0_alpha)
DEALLOCATE (chi0_pert_alpha)
DEALLOCATE (grad_chi_alpha)
DEALLOCATE (EM_x)
DEALLOCATE (h_alpha)
CALL Deallocate_dims (dims)
CALL Deallocate_model_vars (Bg0)
CALL Deallocate_CVs (dchiB0)
CALL Deallocate_CVs (chi0)
CALL Deallocate_CVs (chi0_pert)
CALL Deallocate_CVs (grad_chi)
CALL Deallocate_CVs (h)
CALL Deallocate_CVT (CVT)
IF (ALLOCATED(fft_wsave_x)) DEALLOCATE(fft_wsave_x)
IF (ALLOCATED(fft_work_x)) DEALLOCATE(fft_work_x)

PRINT*, '  -- done'

END PROGRAM Master_GradTest
