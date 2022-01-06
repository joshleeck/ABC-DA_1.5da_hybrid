! Define namelist
! ---------------
&UserOptions

! Reading and processing UM data
! ------------------------------
  nlongs                   = 364,
  nlevs                    = 60,
  datadirImpliedCov        = '.',
  datadirCVT               = '../Master_Calibration/rf_noVR',
  Hybrid_opt               = 1,                  !Type of hybrid (or if pure Var)
  NEnsMems                 = 30,
  Cov_weightE              = 0.0,
  Cov_weightC              = 100.0,
  hScale_alpha             = 100000.0,
  vScale_alpha             = 5000.0,
  InterVarLoc              = .FALSE.
  datadirEM                = '../Master_Calibration/Master_Calibration_stage2/hybrid',
  CVT_file                 = 'CVT_rf_noVR.nc',
  datadirABCfcs            = '../Master_Calibration/Master_Calibration_stage1',
  LS_file                  = 'FC_Ens001_Item001.nc',
  ImplCov_npoints          = 1
  longindex(1)             = 182
  levindex(1)              = 30
/
