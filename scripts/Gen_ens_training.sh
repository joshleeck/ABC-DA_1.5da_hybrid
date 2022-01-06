#!/bin/bash

# ========================================================
# To generate the ensemble forecasts as training data for use
# in Master_Calibration_stage2
# Written by Joshua Lee, 08/06/2021
# ========================================================


# ===== USER-SPECIFIED VARIABLES =========================

# The size of the system
NLONGS=364
NLEVS=60
DX=1500.0

# The directory at the base of this cycling (output data - MUST BE FULL PATH)
BASE_DIR=/scratch/leeck/ABC-DA_1.5da_leeck_364/test

# Data assimilation window length (seconds)
DA_WINDOW=3600.0
DA_WINDOW_INT=3600

# The directory containing the main source files
CODE_DIR=/scratch/leeck/ABC-DA_1.5da_leeck_364/src

# Set the number of ensemble members excluding control/deterministic member
N_ENS=30

# Set to generate ensemble members and forecast for cold start (1 or 0) - set to 0 to skip
CALCS_ENS=1

# The file containing the initial analyses
INIT_ANA=ABC_InitialConds_SINGV.nc

# The directory containing initial analyses
INIT_ANA_DIR=/scratch/leeck/ABC-DA_1.5da_leeck_364/examples/Master_PrepareABC_InitState/

# The file containing long truth run for drawing random fields from
RF_FILE=/scratch/leeck/ABC-DA_1.5da_leeck_364/examples/Master_RunNLModel/ABC_ModelRun_SINGV_for_RF.nc

# Set the model parameters
f=1.0E-5
A=0.02
B=0.01
C=1.0E4

# ===== END OF USER-SPECIFIED VARIABLES ==================


echo "=============================================" > $BASE_DIR/SuiteOut
echo "Running ABC_DA suite" >> $BASE_DIR/SuiteOut
date >> $BASE_DIR/SuiteOut
echo "=============================================" >> $BASE_DIR/SuiteOut


# =====================================================================================
# GENERATE ENSEMBLE MEMBERS AND THEN PERFORM FORECASTS
# =====================================================================================

if [ $CALCS_ENS == 1 ]; then
  # Perform cold start to generate the initial ensemble members (cycle 0)
  # Control member analysis and forecast is assumed to be already generated
  # ------------------------------------
  # Ensemble forecasts also provide the error modes for hybrid EnVar for future cycles
  echo "===== COLD START GENERATION OF ENSEMBLE =====" >> $BASE_DIR/SuiteOut

  mkdir -p $BASE_DIR/InitEns
  # Rename file to convention
  cp $INIT_ANA_DIR/$INIT_ANA $BASE_DIR/InitEns/ABC_anal_control.nc

  echo "Generating initial conditions for $N_ENS ensemble members" >> $BASE_DIR/SuiteOut
  #Generate namelist (overwrites previous namelist)
  cat > $BASE_DIR/InitEns/UserOptions.nl << EOF
&UserOptions
! Make ensemble of perturbed analyses
! ------------------------------
  ABC_anal_ctrl_file       = 'ABC_anal_control.nc'
  datadirABC_anal          = '$BASE_DIR/InitEns'
  datadirABCEns_anal       = '$BASE_DIR/InitEns'
  NEnsMems                 = $N_ENS
  nlongs                   = $NLONGS
  nlevs                    = 60
  Breed_opt                = 3
  RF_tune                  = 1.0
  RF_file                  = '$RF_FILE'
  uncorr_thresh            = 100
/
EOF
  cd $BASE_DIR/InitEns
  $CODE_DIR/Master_MakeEns.out > stdout 2> stderr
  mv stdout MakeEns_stdout
  mv stderr MakeEns_stderr
  echo "  ... done" >> $BASE_DIR/SuiteOut

  for ENS in $(seq 1 $N_ENS)
  do  
    ENS_NUM=`printf "%03d\n" $ENS`
    # Perform forecasts from initial conditions
    echo "Performing forecast for ensemble $ENS" >> $BASE_DIR/SuiteOut
    # Generate namelist (overwrites previous namelist)
    cat > $BASE_DIR/InitEns/UserOptions.nl << EOF
&UserOptions
! Running ABC model
! ------------------------------
  datadirABC_in            = '.'
  init_ABC_file            = 'ABC_anal_Ens${ENS_NUM}.nc'
  datadirABC_out           = '.'
  output_ABC_file          = 'ABC_bg_Ens${ENS_NUM}.nc'
  diagnostics_file         = 'tmp.dat'
  runlength                = $DA_WINDOW
  ndumps                   = 6
  dt                       = 1.0
  nlongs                   = $NLONGS
  nlevs                    = $NLEVS
  Lengthscale_diagnostics  = .FALSE.
  f                        = $f
  A                        = $A
  B                        = $B
  C                        = $C
  Adv_tracer               = .FALSE.
/
EOF
    $CODE_DIR/Master_RunNLModel.out > stdout 2> stderr
    echo "  ... done" >> $BASE_DIR/SuiteOut
  done
 
  # Generate forecast for control member
  echo "Performing forecast for control member" >> $BASE_DIR/SuiteOut
  # Generate namelist (overwrites previous namelist)

    cat > $BASE_DIR/InitEns/UserOptions.nl << EOF
&UserOptions
! Running ABC model
! ------------------------------
  datadirABC_in            = '.'
  init_ABC_file            = 'ABC_anal_control.nc'
  datadirABC_out           = '.'
  output_ABC_file          = 'ABC_bg_control.nc'
  diagnostics_file         = 'tmp.dat'
  runlength                = $DA_WINDOW
  ndumps                   = 6
  dt                       = 1.0
  nlongs                   = $NLONGS
  nlevs                    = $NLEVS
  Lengthscale_diagnostics  = .FALSE.
  f                        = $f
  A                        = $A
  B                        = $B
  C                        = $C
  Adv_tracer               = .FALSE.
/
EOF
    $CODE_DIR/Master_RunNLModel.out > stdout 2> stderr
    echo "  ... done" >> $BASE_DIR/SuiteOut

fi

echo "Script finished"
