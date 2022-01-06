! Define namelist
! ---------------
&UserOptions

! Run the data assimilation
! -------------------------
  Vartype          = 3,                  !3=3DVar, 35=3D-FGAT, 4=4DVar
  Hybrid_opt       = 1,                  !Type of hybrid (or if pure Var)
  NEnsMems         = 60,
  nlongs           = 364,
  nlevs            = 60
  datadir_Bg       = '../Master_MakeBgObs/MakeBg',
  Bg_file          = 'Bg.nc',
  datadirCVT       = '../Master_Calibration',
  CVT_file         = 'CVT.nc',
  datadir_Obs      = '../Master_MakeBgObs/MakeObs',
  Obs_file         = 'Obs.dat',
  t0               = 0,                  !Time of start of this DA cycle
  N_outerloops     = 1,
  N_innerloops_max = 50,
  crit_inner       = 0.01,
  Cov_weightE      = 0.0,
  Cov_weightC      = 100.0,
  hScale_alpha     = 200000.0,
  vScale_alpha     = 8000.0,
  InterVarLoc      = .FALSE.
  datadirEM        = '../Master_Calibration//Master_Calibration_stage2/hybrid'
  datadirAnal      = '.',
  anal_file        = 'Anal.nc',
  analinc_file     = 'AnalInc.nc',
  diagnostics_file = 'diagnostics.dat'
/
