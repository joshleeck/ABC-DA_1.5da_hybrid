#!/bin/bash

# ========================================================
# To run an ABCvn1.5da data assimilation cycle
# Author, Ross Bannister
# 01/08/2018
# ========================================================
# This script assumes that the following already exist:
#   CVT file, truth file, observation network specification file
# ========================================================


# ===== USER-SPECIFIED VARIABLES =========================

# The size of the system
NLONGS=3600
NLEVS=150
DX=150.0

# The directory at the base of this cycling (output data - MUST BE FULL PATH)
OBS_EXP=vrhob
BAL_EXP=Exp+GB+HB-AB
TRU_EXP=Truth_011
BASE_DIR=/scratch2/sws98rnb/HiRes/Assim/AtmosParams/With$TRU_EXP/$BAL_EXP/Assim_Obs_$OBS_EXP

# The number of data assimilation cycles
N_CYCLES=30

# The number of outer loops
N_OUTER_LOOPS=1

# The maximum number of inner loops
N_INNER_LOOPS_MAX=100

# Data assimilation window length (seconds)
DA_WINDOW=3600.0
DA_WINDOW_INT=3600

# The directory containing the file that specifies which observations to make (input data only)
OBS_NETWORK_DIR=/scratch2/sws98rnb/HiRes/Assim/ObsConfigs/$OBS_EXP

# The directory containing the true state at the initial time (input data only)
INITIAL_TRUTH_DIR=/scratch2/sws98rnb/HiRes/Truths/AtmosParams/$TRU_EXP/RunNLModel

# The directory containing the CVT file (describing the background error covariances) (input data only)
CVT_DIR=/scratch2/sws98rnb/HiRes/Ens+Calib/AtmosParams/$BAL_EXP

# The directory containing the main source files
CODE_DIR=/home/users/sws98rnb/DataAssim/RuthsModel/ABC_vn1.5da/src

# The directory containing the plot source files
PLOT_CODE_DIR=/home/users/sws98rnb/DataAssim/RuthsModel/ABC_vn1.5da/graphics

# Set if the background state at the start has already been computed from the truth (1 or 0)
BACKGROUND_ALREADY_COMPUTED=0

# Set the Var assim type 3=3DVar, 35=3D-FGAT, 4=4DVar
VAR_TYPE=3

# Set the model parameters
A=0.01
B=0.01
C=1.0E5

# Set to perform plotting of results (1 or 0)
PLOT=1

# Set to run the ABC programs (1 or 0) - set to 0 to do plotting only
CALCS=1

# ===== END OF USER-SPECIFIED VARIABLES ==================



# Set the name of the (N_OUTER_LOOPS+1)th LS state produced by the DA
# (the last time state in this file is the background for the next cycle)
let N_OUTER_LOOPSP1=$N_OUTER_LOOPS+1
OUTER_LOOP_FORM=`printf "%03u" $N_OUTER_LOOPSP1`
USUAL_BACKGROUND_FILE=LS_Oloop${OUTER_LOOP_FORM}_Iloop000.nc



echo "=============================================" > $BASE_DIR/SuiteOut
echo "Running ABC_DA suite" >> $BASE_DIR/SuiteOut
date >> $BASE_DIR/SuiteOut
echo "=============================================" >> $BASE_DIR/SuiteOut

echo "List of directories used in this cycling" > $BASE_DIR/ExpList.dat

# Set-up python
module load python2.7/canopy/2.1.8


# =====================================================================================
# START THE CYCLING
# =====================================================================================
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
      # Generate the namelist (see Sect. 4.8 of documentation, under Generate_mode=3)
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
  init_ABC_file = 'Truth.nc'
  datadirCVT    = '$CVT_DIR'
  CVT_file      = 'CVT.nc'
  Bg_inflation  = 5.0
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
    TRUTH_FILE=Truth.nc

    # Make sure that the first DA cycle can find the background file
    CYCLE_DIR_PREV=$CYCLE_DIR/InitBg
    BACK_FILE=Bg.nc

  else
    # We have moved beyond the first cycle
    echo "  Background state to come from previous cycle" >> $BASE_DIR/SuiteOut
    # Set where the truth is for generating the obs
    TRUTH_DIR=$CYCLE_DIR_PREV/Obs+Truth
    TRUTH_FILE=Truth.nc
    # Set the name of the background file
    BACK_FILE=$USUAL_BACKGROUND_FILE
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
  output_ABC_file = 'Truth.nc'
  dt_da           = 600.0
  t0              = 0
  runlength       = $DA_WINDOW
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
  nlongs           = $NLONGS
  nlevs            = $NLEVS
  dx               = $DX
  Vartype          = $VAR_TYPE,          !3=3DVar, 35=3D-FGAT, 4=4DVar
  Hybrid_opt       = 1,                  !Type of hybrid (or if pure Var)
  datadir_Bg       = '$CYCLE_DIR_PREV',
  Bg_file          = '$BACK_FILE',
  datadirCVT       = '$CVT_DIR',
  CVT_file         = 'CVT.nc',
  datadir_Obs      = 'Obs+Truth',
  Obs_file         = 'Obs.dat',
  t0               = 0,                  !Time of start of this DA cycle
  N_outerloops     = $N_OUTER_LOOPS,
  N_innerloops_max = $N_INNER_LOOPS_MAX,
  crit_inner       = 0.0001,
  ForceCor         = .FALSE.,
  datadirAnal      = '.',
  anal_file        = 'Anal.nc',
  analinc_file     = 'AnalInc.nc',
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




  # Run plotting code if required
  if [ $PLOT == 1 ]; then
    echo "  Plotting/processing ..." >> $BASE_DIR/SuiteOut
    python $PLOT_CODE_DIR/PlotAssimDiags.py $CYCLE_DIR
    echo "  ... done" >> $BASE_DIR/SuiteOut
  fi

  # Prepare for the next cycle
  CYCLE_DIR_PREV=$CYCLE_DIR

done


# Run a free forecast from the first background
# ---------------------------------------------
echo "Preparing to run a free forecast from the initial background" >> $BASE_DIR/SuiteOut
let FULL_RUN_LENGTH=$N_CYCLES*$DA_WINDOW_INT
mkdir -p $BASE_DIR/Master_RunNLModel_Fullbg
if [ $CALCS == 1 ]; then
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
  Lengthscale_diagnostics  = .FALSE.
  A                        = $A
  B                        = $B
  C                        = $C
  Adv_tracer               = .TRUE.
/
EOF
fi

# Go into this directory and run the forecast
cd $BASE_DIR/Master_RunNLModel_Fullbg
if [ $CALCS == 1 ]; then
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

echo "Script finished"
