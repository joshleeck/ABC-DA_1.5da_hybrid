! Define namelist
! ---------------
&UserOptions

! Master_PrepareABC_InitState
! ------------------------------
  Init_ABC_opt             = 1
  datadirUM                = 'test/data'
  init_um_file             = 'SINGV_20190607T0100Z.nc'
  datadirABC_out           = '.'
  init_ABC_file            = 'ABC_InitialConds_SINGV.nc'
  latitude                 = 530
  nlongs                   = 364
  nlevs                    = 60
  Regular_vert_grid        = .TRUE.
  Adv_tracer               = .FALSE.
  gravity_wave_switch      = .FALSE.
  f                        = 1.0E-5
  A                        = 0.02
  B                        = 0.01
  C                        = 1.0E4
  BoundSpread              = 100.
/
