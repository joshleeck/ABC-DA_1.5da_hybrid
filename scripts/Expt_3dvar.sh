#!/bin/bash

# ========================================================
# To run an ABCvn1.5da cycling experiment
# Ross Bannister
# Adapted by Joshua Lee, 08/06/2021
# ========================================================
# This script assumes that the following already exist:
#   CVT file, truth file, observation network specification file
# ========================================================


# ===== USER-SPECIFIED VARIABLES =========================

# The size of the system
NLONGS=364
NLEVS=60
DX=1500.0

# The directory at the base of this cycling (output data - MUST BE FULL PATH)
BASE_DIR=/scratch/leeck/ABC-DA_1.5da_leeck_364/test

# The number of data assimilation cycles
N_CYCLES=50

# The number of outer loops
N_OUTER_LOOPS=1

# The maximum number of inner loops
N_INNER_LOOPS_MAX=75

# Data assimilation window length (seconds)
DA_WINDOW=3600.0
DA_WINDOW_INT=3600

# The directory containing the file that specifies which observations to make (input data only)
OBS_NETWORK_DIR=$BASE_DIR/ObsConfig

# The directory containing the true state at the initial time (input data only)
INITIAL_TRUTH_DIR=$BASE_DIR/RunNLModel

# The directory containing the CVT file (describing the background error covariances) (input data only)
CVT_DIR=$BASE_DIR/CVT

# The directory containing the main source files
CODE_DIR=/scratch/leeck/ABC-DA_1.5da_leeck_364/src

# The directory containing the plot source files
PLOT_CODE_DIR=/scratch/leeck/ABC-DA_1.5da_leeck_364/graphics

# Set if the background state at the start has already been computed from the truth (1 or 0)
BACKGROUND_ALREADY_COMPUTED=1

# Set the Var assim type 3=3DVar, 35=3D-FGAT, 4=4DVar
VAR_TYPE=35

# Set the hybrid option 1=standard B, 2=pure Envar, 3=hybrid Envar
HYBRID_OPT=1

# Set the number of ensemble members excluding control/deterministic member
N_ENS=30

# Set initial control analysis and its forecast for generation of ensemble (if required)
INIT_ANA=/scratch/leeck/ABC-DA_1.5da_leeck_364/examples/Master_PrepareABC_InitState/ABC_InitialConds_SINGV.nc
INIT_BG=/scratch/leeck/ABC-DA_1.5da_leeck_364/examples/Master_RunNLModel/ABC_ModelRun_SINGV.nc

# Set to generate ensemble members and forecast for cold start (1 or 0) - set to 0 to skip
CALCS_ENS=0

# Set to generate free background run
CALCS_FREE=0

# Set the model parameters
f=1.0E-5
A=0.02
B=0.01
C=1.0E4

# Set to perform plotting of results (1 or 0)
PLOT=1

# Set to perform plotting for ensemble (1 or 0)
PLOT_ENS=0

# Set to run the ABC programs (1 or 0) - set to 0 to do plotting only
CALCS=1

# ===== END OF USER-SPECIFIED VARIABLES ==================

# Set the name of the (N_OUTER_LOOPS+1)th LS state produced by the DA
# The last time state in this file is the background for the next cycle
# The first time state in this file is the analysis
let N_OUTER_LOOPSP1=$N_OUTER_LOOPS+1
OUTER_LOOP_FORM=`printf "%03u" $N_OUTER_LOOPSP1`
USUAL_BACKGROUND_FILE=LS_Oloop${OUTER_LOOP_FORM}_Iloop000.nc



echo "=============================================" > $BASE_DIR/SuiteOut
echo "Running ABC_DA suite" >> $BASE_DIR/SuiteOut
date >> $BASE_DIR/SuiteOut
echo "=============================================" >> $BASE_DIR/SuiteOut

echo "List of directories used in this cycling" > $BASE_DIR/ExpList.dat


# =====================================================================================
# START THE CYCLING
# =====================================================================================

if [ $CALCS_ENS == 1 ]; then
  # Perform cold start to generate the initial ensemble members (cycle 0)
  # Control member analysis and forecast is assumed to be already generated
  # ------------------------------------
  # Ensemble forecasts also provide the error modes for hybrid EnVar for future cycles
  echo "===== COLD START GENERATION OF ENSEMBLE =====" >> $BASE_DIR/SuiteOut

  mkdir -p $BASE_DIR/InitEns
  echo "Generating initial conditions for $N_ENS ensemble members" >> $BASE_DIR/SuiteOut
  
  cp $INIT_ANA $BASE_DIR/InitEns/ABC_anal_control.nc
  cp $INIT_BG $BASE_DIR/InitEns/ABC_bg_control.nc

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
  Ens_opt                  = 3
  RF_tune                  = 1.0
  RF_file                  = '/scratch/leeck/ABC-DA_1.5da_leeck_364/examples/Master_RunNLModel/ABC_ModelRun_SINGV_for_RF.nc'
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

  # Compute error modes by subtracting control member forecast fields
  # Scaling is done by Master_Assimilate
  cdo ensmean ABC_bg_Ens*.nc Ensmean.nc

  for ENS in $(seq 1 $N_ENS)
  do
    ENS_NUM=`printf "%03d\n" $ENS`
    #cdo sub ABC_bg_Ens${ENS_NUM}.nc ABC_bg_control.nc PertABC_Item${ENS_NUM}.nc
    cdo sub ABC_bg_Ens${ENS_NUM}.nc Ensmean.nc PertABC_Item${ENS_NUM}.nc
  done

  # Tidy up
  rm -f tmp.dat
  
fi



for CYCLE in $(seq 1 $N_CYCLES)
do
  CYCLE_FORM=`printf "%04u" $CYCLE`
  echo "===== DATA ASSIMILATION CYCLE NUMBER" $CYCLE_FORM "=====" >> $BASE_DIR/SuiteOut
  CYCLE_DIR=$BASE_DIR/da_cycle_$CYCLE_FORM

  # Make the directory for this da cycle
  mkdir -p $CYCLE_DIR

  echo $CYCLE_DIR >> $BASE_DIR/ExpList.dat


  if [ $CYCLE == 1 ]; then
    # Create a background state by perturbing the truth
    # -------------------------------------------------
    echo "  First cycle - need a background state" >> $BASE_DIR/SuiteOut
    if [ $BACKGROUND_ALREADY_COMPUTED == 0 ]; then
      mkdir -p $CYCLE_DIR/InitBg
      # Generate the namelist
      if [ $CALCS == 1 ]; then
        cat > $CYCLE_DIR/InitBg/UserOptions.nl << EOF
&UserOptions
! Generate background
! ------------------------------
  nlongs        = $NLONGS
  nlevs         = $NLEVS
  dx            = $DX
  Generate_mode = 3
  datadirABC_in = '$INITIAL_TRUTH_DIR'
  init_ABC_file = 'Truth_run.nc'
  datadirCVT    = '$CVT_DIR'
  CVT_file      = 'CVT.nc'
  Bg_inflation  = 1.0
  datadir_Bg    = '.'
  Pert_file     = 'Bgerr.nc'
  Bg_file       = 'Bg.nc'
/
EOF
      fi

      # Go into this directory and run the code to generate the background
      cd $CYCLE_DIR/InitBg
      if [ $CALCS == 1 ]; then
        echo "  Generating the first background state ..." >> $BASE_DIR/SuiteOut
        $CODE_DIR/Master_MakeBgObs.out > stdout 2> stderr
        echo "  ... done" >> $BASE_DIR/SuiteOut
      fi

    fi

    # Set where the truth is for generating the obs
    TRUTH_DIR=$INITIAL_TRUTH_DIR
    TRUTH_FILE=Truth_run.nc

    # Make sure that the first DA cycle can find the background file
    CYCLE_DIR_PREV=$CYCLE_DIR/InitBg
    BACK_FILE=Bg.nc
    EM_DIR=$BASE_DIR/InitEns

  else
    # We have moved beyond the first cycle
    echo "  Background state to come from previous cycle" >> $BASE_DIR/SuiteOut
    # Set where the truth is for generating the obs
    TRUTH_DIR=$CYCLE_DIR_PREV/Obs+Truth
    TRUTH_FILE=Truth_run.nc
    # Set the name of the background file
    BACK_FILE=$USUAL_BACKGROUND_FILE
    EM_DIR=$CYCLE_DIR_PREV/Ens
  fi

  echo "  This cycle's dir :" $CYCLE_DIR >> $BASE_DIR/SuiteOut
  echo "  Truth            :" $TRUTH_DIR/$TRUTH_FILE >> $BASE_DIR/SuiteOut
  echo "  Background       :" $CYCLE_DIR_PREV/$BACK_FILE >> $BASE_DIR/SuiteOut






  # Make the observations for this cycle
  # ------------------------------------
  # This also generates the truth state over the assimilation window
  mkdir -p $CYCLE_DIR/Obs+Truth
  if [ $CALCS == 1 ]; then
    cat > $CYCLE_DIR/Obs+Truth/UserOptions.nl << EOF
&UserOptions
! Generate observations
! ------------------------------
  nlongs          = $NLONGS
  nlevs           = $NLEVS
  dx              = $DX
  Generate_mode   = 2
  datadir_ObsSpec = '$OBS_NETWORK_DIR'
  ObsSpec_file    = 'ObsSpec.dat'
  datadir_Obs     = '.'
  Obs_file        = 'Obs.dat'
  datadirABC_in   = '$TRUTH_DIR'
  init_ABC_file   = '$TRUTH_FILE'
  output_ABC_file = 'Truth_run.nc'
  dt_da           = 600.0
  t0              = 0
  runlength       = $DA_WINDOW
  f               = $f
  A               = $A
  B               = $B
  C               = $C
/
EOF
  fi

  # Go into this directory and run the code to generate the observations
  cd $CYCLE_DIR/Obs+Truth
  if [ $CALCS == 1 ]; then
    echo "  Generating the synthetic obs ..." >> $BASE_DIR/SuiteOut
    $CODE_DIR/Master_MakeBgObs.out > stdout 2> stderr
    echo "  ... done" >> $BASE_DIR/SuiteOut
  fi




  # Run the data assimilation for this cycle
  # ----------------------------------------
  # At the end this also produces the background for the next cycle
  if [ $CALCS == 1 ]; then
    cat > $CYCLE_DIR/UserOptions.nl << EOF
&UserOptions
! Run the data assimilation
! -------------------------
  Vartype          = $VAR_TYPE                 !3=3DVar, 35=3D-FGAT, 4=4DVar
  Hybrid_opt       = $HYBRID_OPT        !Type of hybrid (or if pure Var)
  NEnsMems         = $N_ENS
  nlongs           = $NLONGS
  nlevs            = $NLEVS
  datadir_Bg       = '$CYCLE_DIR_PREV'
  Bg_file          = '$BACK_FILE'
  datadirCVT       = '$CVT_DIR'
  CVT_file         = 'CVT.nc'
  datadir_Obs      = 'Obs+Truth'
  Obs_file         = 'Obs.dat'
  t0               = 0                  !Time of start of this DA cycle
  N_outerloops     = $N_OUTER_LOOPS
  N_innerloops_max = $N_INNER_LOOPS_MAX
  crit_inner       = 0.001
  Cov_weightE      = 0.0
  Cov_weightC      = 100.0
  hScale_alpha     = 20000.0           !10000 20000 5000 40000
  vScale_alpha     = 2000000.0             !1000.0 2000.0 500 4000
  InterVarLoc      = .FALSE.
  datadirEM        = '$EM_DIR'
  datadirAnal      = '.'
  anal_file        = 'Anal.nc',
  analinc_file     = 'AnalInc.nc'
  diagnostics_file = 'diagnostics.dat'
/
EOF
  fi

  # Go into this directory and run the code to do the data assimilation
  cd $CYCLE_DIR
  if [ $CALCS == 1 ]; then
    echo "  Assimilation ..." >> $BASE_DIR/SuiteOut
    $CODE_DIR/Master_Assimilate.out > stdout 2> stderr
    echo "  ... done" >> $BASE_DIR/SuiteOut
  fi


  # Generate new ensemble forecast for next cycle using error breeding, centred on analysis
  mkdir -p $CYCLE_DIR/Ens
  cd $CYCLE_DIR/Ens
  # Use cold start directory for the first cycle
  if [ $CYCLE == 1 ]; then
    ENS_DIR=$BASE_DIR/InitEns
    OBS_DIR=$BASE_DIR/Obs+Truth
  # For other cycles, use previous cycle directory
  else
    ENS_DIR=$CYCLE_DIR_PREV/Ens
    OBS_DIR=$CYCLE_DIR_PREV/Obs+Truth
  fi
    
  # Generate namelist
  cat > $CYCLE_DIR/Ens/UserOptions.nl << EOF
&UserOptions
! Generate ensemble of perturbed analyses
! ------------------------------
  ABC_init_ctrl_file       = 'ABC_anal_control.nc'
  ABC_anal_ctrl_file       = 'Anal.nc'
  ABC_bg_ctrl_file         = 'ABC_bg_control.nc'
  datadirABC_init          = '$ENS_DIR'
  datadirABC_anal          = '$CYCLE_DIR'
  datadirABC_bg            = '$ENS_DIR'
  datadirABCEns_init       = '$ENS_DIR'
  datadirABCEns_bg         = '$ENS_DIR'
  datadirABCEns_anal       = '$CYCLE_DIR/Ens'
  NEnsMems                 = $N_ENS
  Ens_opt                  = 2
  epsilon_diagnostic       = 'epsilon_diag.dat'
  epsilon_file             = '$BASE_DIR/da_cycle_0001/Ens/epsilon_diag.dat'
  nlongs                   = $NLONGS
  nlevs                    = $NLEVS
/
EOF

  if [ $CALCS_ENS == 1 ]; then
    echo "  Producing ensemble members using EnSRF ..." >> $BASE_DIR/SuiteOut
    $CODE_DIR/Master_MakeEns.out > stdout 2> stderr
    mv stdout MakeEns_stdout
    mv stderr MakeEns_stderr
    echo "  ... done" >> $BASE_DIR/SuiteOut

    # Perform forecast and compute error modes
    for ENS in $(seq 1 $N_ENS)
    do
      ENS_NUM=`printf "%03d\n" $ENS`

      # Perform forecasts from initial conditions
      echo "  Performing forecast from perturbed analysis for ensemble $ENS ..." >> $BASE_DIR/SuiteOut
      # Generate namelist (overwrites previous namelist)
      cat > $CYCLE_DIR/Ens/UserOptions.nl << EOF
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

    # Compute error modes by subtracting control member forecast fields
    # Scaling is done by Master_Assimilate
    cdo ensmean ABC_bg_Ens*.nc Ensmean.nc
    
    for ENS in $(seq 1 $N_ENS)
    do
      ENS_NUM=`printf "%03d\n" $ENS`
      #cdo sub ABC_bg_Ens${ENS_NUM}.nc $CYCLE_DIR/$USUAL_BACKGROUND_FILE PertABC_Item${ENS_NUM}.nc
      cdo sub ABC_bg_Ens${ENS_NUM}.nc Ensmean.nc PertABC_Item${ENS_NUM}.nc
    done

    # Tidy up, copy relevant files to be used for next cycle
    cp $CYCLE_DIR/$USUAL_BACKGROUND_FILE ABC_bg_control.nc
    cp $CYCLE_DIR/Anal.nc ABC_anal_control.nc 
    rm -f tmp.dat
  fi

  # Run plotting code if required
  if [ $PLOT == 1 ]; then
    echo "  Plotting/processing ..." >> $BASE_DIR/SuiteOut
    python $PLOT_CODE_DIR/PlotAssimDiags.py $CYCLE_DIR
    echo "  ... done" >> $BASE_DIR/SuiteOut
  fi

  # Run plotting code for ensemble members if required
  if [ $PLOT_ENS == 1 ]; then
    echo " Plotting/processing ..." >> $BASE_DIR/SuiteOut
    python $PLOT_CODE_DIR/PlotAssimDiagsEns.py $CYCLE_DIR $N_ENS
    echo "  ... done"
  fi

  # Compute ensemble spread for each cycle
  if [ $PLOT_ENS == 1 ]; then

    arrFiles=()
    for ENS in $(seq 1 $N_ENS)
    do
      ENS_NUM=`printf "%03d\n" $ENS`
      arrFiles+=("ABC_bg_Ens${ENS_NUM}.nc")
    done
    #echo ${arrFiles[@]}
    cdo ensstd1 ${arrFiles[@]} ABC_bg_ensstd1.nc
  fi

  # Prepare for the next cycle
  CYCLE_DIR_PREV=$CYCLE_DIR

done

# Run a free forecast from the first background
# ---------------------------------------------
echo "Preparing to run a free forecast from the initial background" >> $BASE_DIR/SuiteOut
let FULL_RUN_LENGTH=$N_CYCLES*$DA_WINDOW_INT
mkdir -p $BASE_DIR/Master_RunNLModel_Fullbg
if [ $CALCS_FREE == 1 ]; then
  cat > $BASE_DIR/Master_RunNLModel_Fullbg/UserOptions.nl << EOF
! Define namelist
! ---------------
&UserOptions
! Running ABC model
! ------------------------------
  nlongs                   = $NLONGS
  nlevs                    = $NLEVS
  dx                       = $DX
  datadirABC_in            = '../da_cycle_0001/InitBg'
  init_ABC_file            = 'Bg.nc'
  datadirABC_out           = '.'
  output_ABC_file          = 'BgFc.nc'
  diagnostics_file         = 'ABC_Diagnostics.dat'
  runlength                = $FULL_RUN_LENGTH.0
  ndumps                   = $N_CYCLES
  dt                       = 1.0
  Lengthscale_diagnostics  = .FALSE.
  f                        = $f
  A                        = $A
  B                        = $B
  C                        = $C
  Adv_tracer               = .TRUE.
/
EOF
fi

# Go into this directory and run the forecast
cd $BASE_DIR/Master_RunNLModel_Fullbg
if [ $CALCS_FREE == 1 ]; then
  echo "Running a free forecast from initial background ..." >> $BASE_DIR/SuiteOut
  $CODE_DIR/Master_RunNLModel.out > stdout 2> stderr
  echo "... done" >> $BASE_DIR/SuiteOut
fi



# Run final plotting code to show aggregated diagnostics for the whole cycle
if [ $PLOT == 1 ]; then
  echo "Final plotting/processing ..." >> $BASE_DIR/SuiteOut
  python $PLOT_CODE_DIR/PlotMultiCycleErrors.py $BASE_DIR
  echo "... done" >> $BASE_DIR/SuiteOut
fi

# Run final plotting code to show aggregated diagnostics for whole cycle for each member
if [ $PLOT_ENS == 1 ]; then  
  echo "Final plotting/processing for ensemble ..." >> $BASE_DIR/SuiteOut
  python $PLOT_CODE_DIR/PlotMultiCycleErrorsEns.py $BASE_DIR $N_ENS
  echo "... done" >> $BASE_DIR/SuiteOut
fi

# Run final code to get aggregated diagnostics for the error and ensemble spread
if [ $PLOT_ENS == 1 ]; then
  echo "Final plotting/processing for ensemble spread ..." >> $BASE_DIR/SuiteOut
  python $PLOT_CODE_DIR/PlotMultiCycleSpreadError.py $BASE_DIR
  echo "... done" >> $BASE_DIR/SuiteOut
fi

echo "Script finished"
