! Define namelist
! ---------------
&UserOptions

! Running ABC model
! ------------------------------
  datadirABC_in            = '../Master_PrepareABC_InitState'
  init_ABC_file            = 'ABC_InitialConds_SINGV.nc'
  datadirABC_out           = '.'
  output_ABC_file          = 'ABC_ModelRun_SINGV_for_RF.nc'
  diagnostics_file         = 'ABC_Diagnostics.dat'
  runlength                = 4320000.0
  ndumps                   = 1200
  dt                       = 1.0
  nlongs                   = 364
  nlevs                    = 60
  Lengthscale_diagnostics  = .FALSE.
  f                        = 1.0E-5
  A                        = 0.02
  B                        = 0.01
  C                        = 1.0E4
  Adv_tracer               = .FALSE.
/
