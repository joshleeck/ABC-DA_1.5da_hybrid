! Define namelist
! ---------------
&UserOptions

! Specification of observations
! ------------------------------
  Generate_mode            = 1
  ObsSpec%year0            = 2010
  ObsSpec%month0           = 1
  ObsSpec%day0             = 1
  ObsSpec%hour0            = 0
  ObsSpec%min0             = 0
  ObsSpec%sec0             = 0
  ObsSpec%NumBatches       = 5
  datadir_ObsSpec          = '.'
  ObsSpec_file             = 'ObsSpec.dat'
! ------------------------------
  ObsSpec%batch(1)         = 1
  ObsSpec%seconds(1)       = 0
  ObsSpec%ob_of_what(1)    = 1       ! u
  ObsSpec%NumObs_long(1)   = 20
  ObsSpec%NumObs_height(1) = 10
  ObsSpec%long_min(1)      = 50000.0
  ObsSpec%long_max(1)      = 500000.0
  ObsSpec%height_min(1)    = 4000.0
  ObsSpec%height_max(1)    = 14000.0
  ObsSpec%stddev(1)        = 0.5
! ------------------------------
  ObsSpec%batch(2)         = 1
  ObsSpec%seconds(2)      = 0
  ObsSpec%ob_of_what(2)    = 2       ! v
  ObsSpec%NumObs_long(2)   = 20
  ObsSpec%NumObs_height(2) = 10
  ObsSpec%long_min(2)      = 50000.0
  ObsSpec%long_max(2)      = 500000.0
  ObsSpec%height_min(2)    = 4000.0
  ObsSpec%height_max(2)    = 14000.0
  ObsSpec%stddev(2)        = 0.08
! ------------------------------
  ObsSpec%batch(3)         = 1
  ObsSpec%seconds(3)       = 0
  ObsSpec%ob_of_what(3)    = 3       ! w
  ObsSpec%NumObs_long(3)   = 20
  ObsSpec%NumObs_height(3) = 10
  ObsSpec%long_min(3)      = 50000.0
  ObsSpec%long_max(3)      = 500000.0
  ObsSpec%height_min(3)    = 4000.0
  ObsSpec%height_max(3)    = 14000.0
  ObsSpec%stddev(3)        = 0.00010000
! ------------------------------
  ObsSpec%batch(4)         = 1
  ObsSpec%seconds(4)       = 0
  ObsSpec%ob_of_what(4)    = 4       ! r
  ObsSpec%NumObs_long(4)   = 20
  ObsSpec%NumObs_height(4) = 10
  ObsSpec%long_min(4)      = 50000.0
  ObsSpec%long_max(4)      = 500000.0
  ObsSpec%height_min(4)    = 4000.0
  ObsSpec%height_max(4)    = 14000.0
  ObsSpec%stddev(4)        = 0.00100000
! ------------------------------
  ObsSpec%batch(5)         = 1
  ObsSpec%seconds(5)       = 0
  ObsSpec%ob_of_what(5)    = 5       ! b
  ObsSpec%NumObs_long(5)   = 20
  ObsSpec%NumObs_height(5) = 10
  ObsSpec%long_min(5)      = 50000.0
  ObsSpec%long_max(5)      = 500000.0
  ObsSpec%height_min(5)    = 4000.0
  ObsSpec%height_max(5)    = 14000.0
  ObsSpec%stddev(5)        = 0.00500000
/
