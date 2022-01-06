&UserOptions
! Generate background
! ------------------------------
  Generate_mode = 3
  datadirABC_in = '../../Master_PrepareABC_InitState'
  init_ABC_file = 'ABC_InitialConds_SINGV.nc'
  datadirCVT    = '../../Master_Calibration'
  CVT_file      = 'CVT.nc'
  datadir_Bg    = '.'
  Pert_file     = 'Bgerr.nc'
  Bg_file       = 'Bg.nc'
  nlongs        = 364
  nlevs         = 60
/
