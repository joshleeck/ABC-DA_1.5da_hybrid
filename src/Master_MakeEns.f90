PROGRAM Master_MakeEns

!*****************************************************
!*   Code to run to generate ensemble members        *
!*                                                   *
!*   J. Lee,       1.5_hybrid:  03-06-2021           *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    dims_type,                   &
    ZREAL8,                      &
    Obs_type,                    &
    ABC_type,                    &
    NEnsMems,                    &
    datadirABCEns_init,          &
    datadirABCEns_bg,            &
    datadirABCEns_anal,          &
    fullpath_Obs,                &
    ABC_init_ctrl_file,          &
    ABC_bg_ctrl_file,            &
    ABC_anal_ctrl_file,          &
    datadirABC_init,             &
    datadirABC_anal,             &
    datadirABC_bg,               &
    Ens_opt,                     &
    RF_tune,                     &
    EBV_tune,                    &
    EnSRF_tune,                  &
    RF_file,                     &
    uncorr_thresh,               &
    dt, dt_da, t0,               &
    nlongs,                      &
    nlevs,                       &
    epsilon_diagnostic,          &
    epsilon_file


IMPLICIT NONE

!NetCDF library (file format used to read/write data)
!----------------------------------------------------
INCLUDE '/scratch/singadm/pkg_hera/netcdf/4.3.3/gnu/4.9.2/include/netcdf.inc'

INCLUDE "InnerProdModelSpace.interface"
INCLUDE "ModelObservations.interface"
INCLUDE "ModelObservations_ZeroPert.interface"
INCLUDE "Write_Obs_hx.interface"
INCLUDE "Read_Obs.interface"
INCLUDE "Count_Obs.interface"
INCLUDE "DeAllocate_Obs.interface"

! Declare variables
!==========================
CHARACTER(LEN=256)          :: ABCfilename, Diag_filename, input_filename
INTEGER                     :: n, rint1, rint2, count, ncid, ierr, ntimes, dimid_time, DAtimesteps, timesteps, maxtime, t, obc, n2
REAL(ZREAL8)                :: m1 = -1.0, c, c_star, denom, fac
REAL(ZREAL8)                :: total, magnitude2_init, magnitude2_bg, epsilon_init, max_norm, r1, r2
LOGICAL                     :: exist
TYPE(dims_type)             :: dims
TYPE(Obs_type), POINTER     :: Obs
TYPE(ABC_type)              :: ABC_anal0_data, ABC_bg0_data, ABC_init0_data, RF_data1, RF_data2, tmp_state
TYPE(ABC_type), ALLOCATABLE :: diff_init_pert(:), diff_bg_pert(:)
TYPE(ABC_type), ALLOCATABLE :: ABC_init_data(:), ABC_bg_data(:), ABC_anal_data(:)
TYPE(ABC_type), ALLOCATABLE :: ABC_bg_data_t(:,:)
REAL(ZREAL8), ALLOCATABLE   :: hx(:,:), Yb_pert(:,:), R(:,:), YbYbT(:,:), Gamma(:,:), Yb_pert_T(:,:), YbTGamYb(:,:)
REAL(ZREAL8), ALLOCATABLE   :: Yb_mean(:), IPIV(:)

CHARACTER(LEN=320)          :: ABC_init0_file, ABC_anal0_file, ABC_bg0_file, Obs_filename


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_MakeEns'
PRINT*, '*************************************************************************'


! Read namelist
CALL SetOptions

! Prepare memory
ALLOCATE (diff_init_pert(1:NEnsMems))
ALLOCATE (diff_bg_pert(1:NEnsMems))
ALLOCATE (ABC_bg_data(1:NEnsMems))
ALLOCATE (ABC_init_data(1:NEnsMems))
ALLOCATE (ABC_anal_data(1:NEnsMems))
DO n = 1, NEnsMems
  CALL Initialise_model_vars (ABC_bg_data(n), .FALSE.)
  CALL Initialise_model_vars (ABC_init_data(n), .FALSE.)
  CALL Initialise_model_vars (ABC_anal_data(n), .FALSE.)
  CALL Initialise_model_vars (diff_bg_pert(n), .FALSE.)
  CALL Initialise_model_vars (diff_init_pert(n), .FALSE.)
END DO

CALL Initialise_dims (dims)
CALL Initialise_model_vars (ABC_anal0_data, .FALSE.)
CALL Initialise_model_vars (ABC_bg0_data, .FALSE.)
CALL Initialise_model_vars (ABC_init0_data, .FALSE.)
CALL Initialise_model_vars (RF_data1, .FALSE.)
CALL Initialise_model_vars (RF_data2, .FALSE.)


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_MakeEns'
PRINT*, '*************************************************************************'

ABC_init0_file   = TRIM(datadirABC_init) // '/' // TRIM(ABC_init_ctrl_file)
ABC_anal0_file   = TRIM(datadirABC_anal)  // '/' // TRIM(ABC_anal_ctrl_file)
ABC_bg0_file     = TRIM(datadirABC_bg)  // '/' // TRIM(ABC_bg_ctrl_file)


! Read in ABC control analysis data
PRINT*, 'Reading in control analysis ...'
CALL Read_state_2d (ABC_anal0_file, ABC_anal0_data, dims, -1, .TRUE.)
PRINT*, '-- done'

! Set some commonly-used constants
CALL Set_ht_dep_cons (dims)


IF (Ens_opt == 1 .OR. Ens_opt == 2) THEN
  ! Read in ABC control initial data (previous cycle analysis)
  PRINT*, 'Reading in control initial state ...'
  CALL Read_state_2d (ABC_init0_file, ABC_init0_data, dims, -1, .TRUE.)
  PRINT*, '-- done'

  ! Read in ABC control background data (usually the unperturbed forecast from previous analysis)
  PRINT*, 'Reading in control background ...'
  CALL Read_state_2d (ABC_bg0_file, ABC_bg0_data, dims, -1, .TRUE.)
  PRINT*, '-- done'

  ! Prepare subtraction of control background from each ensemble background
  CALL Mul_model_cons (ABC_bg0_data, m1, .TRUE.)
  ABC_bg0_data % rho(0:nlongs+1,0:nlevs+1) = ABC_bg0_data % rho(0:nlongs+1,0:nlevs+1)*m1

  ! Prepare subtraction of control initial from each ensemble initial state
  CALL Mul_model_cons (ABC_init0_data, m1, .TRUE.)
  ABC_init0_data % rho(0:nlongs+1,0:nlevs+1) = ABC_init0_data % rho(0:nlongs+1,0:nlevs+1)*m1

  ! Read in ABC ensemble initial data (previous cycle perturbed analyses)
  DO n = 1, NEnsMems
    WRITE (ABCfilename, '(A,A,I0.3,A)') TRIM(datadirABCEns_init), '/ABC_anal_Ens', n, '.nc'
    INQUIRE(file=ABCfilename, exist=exist)
    IF (exist) THEN
      ! Read in ensemble forecast states to a appropriate model state array if file exists
      CALL Read_state_2d (ABCfilename, ABC_init_data(n), dims, -1, .TRUE.)
    ELSE
      PRINT *, 'Missing ensemble analyses, or mismatch between NEnsMems and number of files'
      STOP
    END IF
  END DO

  ! Read in ABC ensemble background data
  DO n = 1, NEnsMems
    WRITE (ABCfilename, '(A,A,I0.3,A)') TRIM(datadirABCEns_bg), '/ABC_bg_Ens', n, '.nc'
    INQUIRE(file=ABCfilename, exist=exist)
    IF (exist) THEN
      ! Read in ensemble forecast states to a appropriate model state array if file exists
      CALL Read_state_2d (ABCfilename, ABC_bg_data(n), dims, -1, .TRUE.)
    ELSE
      PRINT *, 'Missing ensemble background, or mismatch between NEnsMems and number of files'
      STOP
    END IF
  END DO
END IF

IF (Ens_opt == 4) THEN
  ! Read observation file, Obs data (input)
  Obs_filename = TRIM(fullpath_Obs)

  ! Initialise the observation structure (a linked list)
  NULLIFY(Obs)

  ! Read-in the observations
  PRINT*, 'Reading in observations ...'
  CALL Read_Obs (Obs, Obs_filename, dt, timesteps, dt_da, DAtimesteps, maxtime)
  PRINT *, 'Model timestep ', dt
  PRINT *, 'DA timestep    ', dt_da
  PRINT *, '-- done'

  ! Count number of observations
  PRINT *, 'Counting observations to determine observation space...'
  CALL Count_Obs (Obs, count)
  PRINT *, '-- done'

  ! Now prepare memory and initialise arrays based on requirements
  ALLOCATE (ABC_bg_data_t(1:NEnsMems,0:DAtimesteps))
  ALLOCATE (hx(1:count,1:NEnsMems))
  ALLOCATE (Yb_pert(1:count,1:NEnsMems))
  ALLOCATE (Yb_mean(1:count))
  ALLOCATE (R(1:count,1:count))
  ALLOCATE (YbYbT(1:count,1:count))
  ALLOCATE (Gamma(1:count,1:count))
  ALLOCATE (Yb_pert_T(1:NEnsMems,1:count))
  ALLOCATE (YbTGamYb(1:NEnsMems,1:NEnsMems))
  
  DO n = 1, NEnsMems
    DO t = 0, DAtimesteps
      CALL Initialise_model_vars (ABC_bg_data_t(n,t), .FALSE.)
    END DO
  END DO


  ! Read in full ABC ensemble background data (i.e. starting from previous analysis time)
  DO n = 1, NEnsMems
    WRITE (ABCfilename, '(A,A,I0.3,A)') TRIM(datadirABCEns_bg), '/ABC_bg_Ens', n, '.nc'
    INQUIRE(file=ABCfilename, exist=exist)
    IF (exist) THEN
      DO t = 0, DAtimesteps
        ! Read in ensemble forecast states to a appropriate model state array if file exists
        ! Previous analysis time is the (t+1)th state
        CALL Read_state_2d (ABCfilename, ABC_bg_data_t(n,t), dims, t+1, .TRUE.)
      END DO
    ELSE
      PRINT *, 'Missing ensemble background, or mismatch between NEnsMems and number of files'
      STOP
    END IF
  END DO
  
  IF (NEnsMems == 1) THEN
    PRINT *, 'Number of ensemble members is 1, which does not make sense for ensemble square root filter'
    STOP
  END IF
END IF



! Code branches to perform various ensemble generation methods from here
! ===================================
! --- BRED VECTORS ---
! ===================================
IF (Ens_opt == 1) THEN
  ! Global scaling of all ensemble members
  ! Compute scaling factor using RMS
  magnitude2_bg   = 0.0
  magnitude2_init = 0.0
  DO n = 1, NEnsMems
    ! Perform dx_k = x_k - x_0 for ensemble background states
    CALL Add_model_vars (diff_bg_pert(n), ABC_bg_data(n), .TRUE.)
    CALL Add_model_vars (diff_bg_pert(n), ABC_bg0_data, .TRUE.)
    ! Compute RMS (does not include rho)
    ! magnitude2_bg = magnitude2_bg + InnerProdModelSpace(diff_bg_pert(n), diff_bg_pert(n))
    
    ! Compute energy norm
    diff_bg_pert(n) % rho(0:nlongs+1,0:nlevs+1) = ABC_bg_data(n) % rho(0:nlongs+1,0:nlevs+1)
    CALL Energy (diff_bg_pert(n), dims)
    magnitude2_bg = magnitude2_bg + diff_bg_pert(n) % Total_Energy

    ! Perform dx_k = x_k - x_0 for ensemble initial states
    CALL Add_model_vars (diff_init_pert(n), ABC_init_data(n), .TRUE.)
    CALL Add_model_vars (diff_init_pert(n), ABC_init0_data, .TRUE.)
    ! Compute RMS (does not include rho)
    ! magnitude2_init = magnitude2_init + InnerProdModelSpace(diff_init_pert(n), diff_init_pert(n))

    ! Compute energy norm
    diff_init_pert(n) % rho(0:nlongs+1,0:nlevs+1) = ABC_init_data(n) % rho(0:nlongs+1,0:nlevs+1)
    CALL Energy (diff_init_pert(n), dims)
    magnitude2_init = magnitude2_init + diff_init_pert(n) % Total_Energy

  END DO

  ! Magnitudes are divided by same total (NEnsMems*nlongs*nlevs*number of variables)
  c = SQRT(magnitude2_init)/SQRT(magnitude2_bg)


! ===================================
! --- ENSEMBLE BRED VECTORS ---
! ===================================
ELSE IF (Ens_opt == 2) THEN
  ! Ensemble bred vectors following Balci et al. (2011)
  max_norm = 0.0

  ! Initial epsilon in terms of average of inner product norm of ensemble members, this caused filter divergence
  ! Since average norm is always less than max norm, perturbations will be scaled down substantially each cycle
  magnitude2_init = 0.0

  DO n = 1, NEnsMems
    ! Perform dx_k = x_k - x_0 for ensemble initial states
    CALL Add_model_vars (diff_init_pert(n), ABC_init_data(n), .TRUE.)
    CALL Add_model_vars (diff_init_pert(n), ABC_init0_data, .TRUE.)
    ! Compute inner product norm (does not include rho)
    !magnitude2_init = magnitude2_init + InnerProdModelSpace(diff_init_pert(n), diff_init_pert(n))


    ! Compute energy norm
    diff_init_pert(n) % rho(0:nlongs+1,0:nlevs+1) = ABC_init_data(n) % rho(0:nlongs+1,0:nlevs+1)
    CALL Energy (diff_init_pert(n), dims)
    magnitude2_init = magnitude2_init + diff_init_pert(n) % Total_Energy 


    ! Perform dx_k = x_k - x_0 for ensemble background states
    CALL Add_model_vars (diff_bg_pert(n), ABC_bg_data(n), .TRUE.)
    CALL Add_model_vars (diff_bg_pert(n), ABC_bg0_data, .TRUE.)
    ! Find maximum size of perturbation (in terms of inner product norm)
    !magnitude2_bg = InnerProdModelSpace(diff_bg_pert(n), diff_bg_pert(n))
   
    ! Find maximum size of perturbation (in terms of energy norm)
    diff_bg_pert(n) % rho(0:nlongs+1,0:nlevs+1) = ABC_bg_data(n) % rho(0:nlongs+1,0:nlevs+1)
    CALL Energy (diff_bg_pert(n), dims)
    magnitude2_bg = diff_bg_pert(n) % Total_Energy

    ! Find largest norm from perturbations
    IF (SQRT(magnitude2_bg) .GT. max_norm) THEN
      max_norm = SQRT(magnitude2_bg)
    END IF
  END DO

  total = 1.0*NEnsMems
  epsilon_init = SQRT(magnitude2_init/total)

  ! Write mean norm as for diagnostic purposes
  Diag_filename   = TRIM(datadirABCEns_anal) // '/' // TRIM(epsilon_diagnostic)

  OPEN (10, file=Diag_filename)
  WRITE (10,*) epsilon_init
  CLOSE (10)

  ! Alternatively, can keep total energy norm constant based on input file to prevent filter divergence
  ! Overwrite epsilon_init with constant from input file
  input_filename = TRIM(epsilon_file)
  OPEN (9, file=input_filename)
  READ (9,*) epsilon_init
  CLOSE (9)

  c = epsilon_init/max_norm


  ! Scale perturbations by 1/SQRT(2) because estimated background error variances are larger than
  ! true background error by factor 2 for state minus state estimates of background error (Var[X-Y] = Var[X] + Var[Y])
  fac = 1/SQRT(2.0)
  c = c*EBV_tune*fac

! ===================================
! --- RANDOM FIELD PERTURBATIONS ---
! ===================================
ELSE IF (Ens_opt == 3) THEN
  ! Random field perturbation following Magnusson et al. (2009)
 
  ! Get 2 random state indexes from a very long simulation
  ! These mimic drawing 2 random analyses states from a time period
  ierr = NF_OPEN(RF_file, NF_NOWRITE, ncid)
  IF ( ierr .NE. 0 ) THEN
    PRINT*, ' *** Error opening file ***'
    PRINT*, ierr, NF_STRERROR(ierr)
    PRINT*, 'FILE :: ', RF_file
    STOP
  ENDIF

  ierr = NF_INQ_DIMID(ncid, 'time', dimid_time)
  IF ( ierr .NE. 0 ) THEN
    PRINT*,'ierr', ierr, NF_STRERROR(ierr) 
    STOP
  END IF

  ierr = NF_INQ_DIMLEN (ncid, dimid_time, ntimes)
  IF (ierr /= 0) THEN
    PRINT*, '***Error getting length of time dimension ***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  END IF

  PRINT *, 'Number of states in RF_file: ', ntimes

  ! Get ensemble pertubations based on difference between 2 random states
  magnitude2_bg = 0
  DO n = 1, NEnsMems
    rint1 = 0
    rint2 = 0
    count = 0
    ! Ensure that states are relatively far apart in time so they are assumed to be uncorrelated
    DO WHILE (ABS(rint1 - rint2) .LT. uncorr_thresh)
      CALL RANDOM_NUMBER(r1)
      rint1 = 1 + FLOOR((ntimes+1-n)*r1)
      CALL RANDOM_NUMBER(r2)
      rint2 = 1 + FLOOR((ntimes+1-n)*r2)
      count = count + 1
      IF (count .GT. 20) THEN
        PRINT *, 'Warning: More than 20 tries to sample uncorrelated random states based on uncorr_thresh'
        PRINT *, 'Consider lowering uncorr_thresh or increasing number of states in RF_file'
        STOP
      END IF
    END DO
    ! Read two random states
    CALL Read_state_2d (RF_file, RF_data1, dims, rint1, .TRUE.)
    CALL Read_state_2d (RF_file, RF_data2, dims, rint2, .TRUE.)

    ! Perform subtraction following equation (8) of Magnusson et al. (2009)
    CALL Add_model_vars (diff_bg_pert(n), RF_data1, .TRUE.)
    CALL Mul_model_cons (RF_data2, m1, .TRUE.)
    RF_data2 % rho(0:nlongs+1,0:nlevs+1) = RF_data2 % rho(0:nlongs+1,0:nlevs+1)*m1
    CALL Add_model_vars (diff_bg_pert(n), RF_data2, .TRUE.)
    PRINT *, 'Perturbation generated using state at time', rint1, 'and', rint2
   
    ! Calculate inner product norm (alternatively, can use energy norm)
    !magnitude2_bg = magnitude2_bg + InnerProdModelSpace(diff_bg_pert(n), diff_bg_pert(n))
    
    ! Calculate energy norm
    diff_bg_pert(n) % rho(0:nlongs+1,0:nlevs+1) = RF_data1 % rho(0:nlongs+1,0:nlevs+1)
    CALL Energy (diff_bg_pert(n), dims)
    magnitude2_bg = magnitude2_bg + diff_bg_pert(n) % Total_Energy

  END DO

  ! Find average of inner product/energy norm from perturbations
  total = 1.0*NEnsMems
  epsilon_init = magnitude2_bg/total

  ! Scale pertubations size to average of inner product norm (or energy norm)
  DO n = 1, NEnsMems
    !magnitude2_bg = InnerProdModelSpace(diff_bg_pert(n), diff_bg_pert(n))
    magnitude2_bg = diff_bg_pert(n) % Total_Energy

    c_star = epsilon_init/magnitude2_bg
    CALL Mul_model_cons (diff_bg_pert(n), c_star, .TRUE.)
  END DO

  ! Further scaling by c is based on tuning factor
  c = RF_tune  

  ! Scale perturbations by 1/SQRT(2) because estimated background error variances are larger than
  ! true background error by factor 2 for state minus state estimates of background error (Var[X-Y] = Var[X] + Var[Y])
  fac = 1/SQRT(2.0)
  c = c*fac

! ===================================
! --- ENSEMBLE SQUARE ROOT FILTER ---
! ===================================
ELSE IF (Ens_opt == 4) THEN
  ! Ensemble square root filter adapted from Sakov and Oke (2009)


  ! Compute Y^b's, ensemble background in observation space using ensemble background
  ! Also compute perturbations of Y^b
  Yb_mean(1:count) = 0.0
  R(1:count,1:count) = 0.0
  YbYbT(1:count,1:count) = 0.0
  hx(1:count,1:NEnsMems) = 0.0

  PRINT *, 'Computing ensemble background in observation space ...'
  DO n = 1, NEnsMems
    ! First pass through observation file - compute ensemble background in observation space
    CALL ModelObservations (DAtimesteps, ABC_bg_data_t(n,0:DAtimesteps), dims, dt_da, t0, &
                          Obs, .FALSE.)
    CALL ModelObservations_ZeroPert (Obs)
    PRINT *, ' -- done'

    ! Second pass through observation file - write out diagnostic file
    ! Also stores values in array hx and R
    WRITE (Diag_filename, '(A,A,I0.3,A)') TRIM(datadirABCEns_anal), '/Obs_Yb_', n, '.dat'
    PRINT *, 'Outputting ensemble background values in observation space ', TRIM(Diag_filename)
    CALL Write_Obs_hx ( TRIM(Diag_filename), Obs, n, NEnsMems, count, hx, R )

    Yb_mean = Yb_mean + hx(:,n)
  END DO
  ! Compute mean of Y^b's 
  denom = FLOAT(NEnsMems)
  Yb_mean = Yb_mean/denom

  ! Compute perturbations of Y^b
  DO n = 1, NEnsMems
    Yb_pert(:,n) = hx(:,n) - Yb_mean
  END DO

  ! Store Y^b^T for use later
  Yb_pert_T = TRANSPOSE(Yb_pert)

  ! Compute Y^b Y^b^T and add R
  ! Another way to compute Y^b Y^b^T
  !DO obc = 1, count
  !  DO obc2 = 1, count
  !    total = 0.0
  !    DO n = 1, NEnsMems
  !      total = total + Yb_pert(obc,n) * Yb_pert(obc2,n)
  !    END DO
  !
  !    YbYbT(obc,obc2) = total
  !  END DO
  !END DO

  YbYbT = MATMUL(Yb_pert,Yb_pert_T)
  denom = FLOAT(NEnsMems-1)
  Gamma = YbYbT/denom + R


  ! Compute ensemble mean and prepare for subtraction from each ensemble background
  ! This is performed for each discrete timestep, although only the final time index is eventually used later
  DO t = 0, DAtimesteps
    CALL Initialise_model_vars (tmp_state, .FALSE.)
    DO n = 1, NEnsMems
      CALL Add_model_vars (tmp_state, ABC_bg_data_t(n,t), .TRUE.)
    END DO
    denom = FLOAT(NEnsMems)
    CALL Div_model_cons (tmp_state, denom, .TRUE.)
    CALL Mul_model_cons (tmp_state, m1, .TRUE.)
    tmp_state % rho(0:nlongs+1,0:nlevs+1) = tmp_state % rho(0:nlongs+1,0:nlevs+1)*m1

    ! Compute perturbations after subtracting ensemble mean
    DO n = 1, NEnsMems
      CALL Add_model_vars (ABC_bg_data_t(n,t), tmp_state, .TRUE.)
    END DO
  END DO

  ! Solve linear system of equations for AX = B; get X=A^-1 B from the solve  

  ALLOCATE (IPIV(1:count)) 

  CALL DGESV(count,                                  & ! The order of A; number of linear equations (in)
             NEnsMems,                               & ! The number of columns of B; number of right-hand sides (in)
             Gamma, count,                           & ! The coefficient matrix A (inout LU from PLU) and leading dimension of A (out)
             IPIV,                                   & ! Pivot indices that define permutation matrix P (out)
             Yb_pert, count,                         & ! The matrix B (inout X) and the leading dimension of B (out)
             ierr)                                     ! got error? 0 => no error (out)


  
  ! Matrix multiplication to get Y^b^T Gamma^-1 Y^b
  YbTGamYb = MATMUL(Yb_pert_T,Yb_pert)

  !PRINT *, 'Yb_pert_T Gamma^-1 Yb_pert'
  !DO n = 1, NEnsMems
  !  PRINT *, YbTGamYb(n,:)
  !END DO

  ! This performs X^b(finaltime) YbTGamYb
  DO n = 1, NEnsMems 
    DO n2 = 1, NEnsMems
      CALL Initialise_model_vars (tmp_state, .FALSE.)
      CALL Add_model_vars (tmp_state, ABC_bg_data_t(n2,DAtimesteps), .TRUE.)

      ! Multiply by correct coefficient in YbTGamYb - each variable is multiplied by same coefficient
      c_star = YbTGamYb(n2, n)
      CALL Mul_model_cons (tmp_state, c_star, .TRUE.)
      CALL Add_model_vars (diff_bg_pert(n), tmp_state, .TRUE.)
    END DO
  END DO
  ! This is the scaling factor that is generalised by Le Duc
  c = EnSRF_tune

ELSE
  PRINT *, 'Ensemble option is invalid'
  STOP
END IF

PRINT *, 'Scaling factor is: ', c

! Scale perturbations before adding to analysis after DA
DO n = 1, NEnsMems
  CALL Mul_model_cons (diff_bg_pert(n), c, .TRUE.)

  ! For EnSRF, need analysis perturbations are computed after scaling
  IF (Ens_opt == 4) THEN
    ! Compute Xb(finaltime) - (c/(N_e-1))*X^b(finaltime) YbTGamYb (c was multiplied above)
    denom = FLOAT(NEnsMems-1)
    CALL Div_model_cons (diff_bg_pert(n), denom, .TRUE.)
    CALL Mul_model_cons (diff_bg_pert(n), m1, .TRUE.)
    CALL Add_model_vars (diff_bg_pert(n), ABC_bg_data_t(n, DAtimesteps), .TRUE.)
  END IF

  ! Set rho field to follow unperturbed analysis
  diff_bg_pert(n) % rho(0:nlongs+1,0:nlevs+1) = 0
  CALL Add_model_vars (ABC_anal_data(n), ABC_anal0_data, .TRUE.)
  CALL Add_model_vars (ABC_anal_data(n), diff_bg_pert(n), .TRUE.)

  ! Write out perturbed analyses to be used for next cycle
  WRITE (ABCfilename, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCEns_anal), '/ABC_anal_Ens', n, '.nc'
  CALL Write_state_2d (ABCfilename, ABC_anal_data(n), dims, 1, 0, 0, .TRUE.)
    
END DO


CALL Deallocate_dims (dims)
CALL Deallocate_model_vars (ABC_anal0_data)
CALL Deallocate_model_vars (ABC_bg0_data)
CALL Deallocate_model_vars (ABC_init0_data)
CALL Deallocate_model_vars (RF_data1)
CALL Deallocate_model_vars (RF_data2)
PRINT *, 'Final tidy up for ensemble-related arrays'
DO n = 1, NEnsMems
  CALL Deallocate_model_vars (ABC_init_data(n))
  CALL Deallocate_model_vars (ABC_bg_data(n))
  CALL Deallocate_model_vars (diff_init_pert(n))
  CALL Deallocate_model_vars (diff_bg_pert(n))
  CALL Deallocate_model_vars (ABC_anal_data(n))
END DO
IF (Ens_opt == 4) THEN
  DO n = 1, NEnsMems
    DO t = 0, DAtimesteps
      CALL Deallocate_model_vars (ABC_bg_data_t(n,t))
    END DO
  END DO
  DEALLOCATE (hx)
  DEALLOCATE (Yb_mean)
  DEALLOCATE (Yb_pert)
  DEALLOCATE (R)
  DEALLOCATE (YbYbT)
  DEALLOCATE (Gamma)
  DEALLOCATE (Yb_pert_T)
  DEALLOCATE (IPIV)
  DEALLOCATE (YbTGamYb)
  CALL Deallocate_model_vars (tmp_state)
  CALL DeAllocate_Obs (Obs)
END IF
PRINT *, 'DONE'

END PROGRAM Master_MakeEns
