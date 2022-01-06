&UserOptions
! Generate observations
! ------------------------------
  Generate_mode   = 2
  datadir_ObsSpec = '../ObsSpecification'
  ObsSpec_file    = 'ObsSpec.dat'
  datadir_Obs     = '.'
  Obs_file        = 'Obs.dat'
  datadirABC_in   = '../../Master_PrepareABC_InitState'
  init_ABC_file   = 'ABC_InitialConds_SINGV.nc'
  output_ABC_file = 'Truth_run.nc'
  dt_da           = 600.0
  t0              = 0
  runlength       = 3600.0
  nlongs          = 364
  nlevs           = 60
  f               = 1.0E-5
  A               = 0.02
  B               = 0.01
  C               = 1.0E4
/
