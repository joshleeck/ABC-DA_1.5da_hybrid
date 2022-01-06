#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read and plot assimilation diagnostics for the ensemble
#
# Please edit the input details (e.g. location of data files)
# at the start of the main part of the code.
# This is located after the function definitions below.
#
# Joshua Lee, September 2021
# -------------------------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import colors, cm
import matplotlib
import os
import sys
from Routines4PlotAssimDiags import *
from ast import literal_eval

# Set variables
# The data directory will be specified by the command line argument if present
if len(sys.argv) > 1:
  data_dir  = sys.argv[1]
else:
  data_dir  = '/home/data'
print ('data_dir = ', data_dir)

try:
  nens        = int(sys.argv[2])
except:
  print "NEns needs to be an integer"
  sys.exit()

# Bound for taking spatial mean over
try:
  bound_tuple = sys.argv[3]
  bound_tuple = literal_eval(bound_tuple)
except:
  print "Need to input proper bounds for taking spatial mean"
  sys.exit()

#=================================================================================

plot_dir    = data_dir + '/Ens'
ScalarDiags = 'diagnostics.dat'
Truth_dir   = data_dir + '/Obs+Truth'
Truth_traj  = 'Truth_run.nc'

# Which routines do we run?
PlotAE_traj = True  # Analysis error trajectory

# Plot all times "AllTimes" or just first and last times "FirstLast" in above options?
# The option "None" is also used to turn-off plotting altogether (and just collect data).
TimeOutput  = "AllTimes"

# Set this to output only a compact (mean error and rms error) time sequence for this da cycle
# If set to true, many of above Plotxxx and TimeOutput settings will be over-ridden.
Output_only_compact_data = True

# Set the domain dimensions
nlongs      = 364
nlevs       = 60
C_param     = 10000.0
f           = 0.00001

# Output type - not yet implemented - will always output to web
# As there could be a large number of plots, there is an option to output results on a web page
# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'web'


# Over-ride settings if necessary
if Output_only_compact_data:
  PlotAE_traj = True
  TimeOutput  = "None"


#os.system('mkdir -p ' + plot_dir + '/Plots')
if (output_type == 'web'):
  filesuffix = '.png'
else:
  filesuffix = '.eps'


if (output_type == 'web'):
  # Set-up the html frames file
  html_file = open (plot_dir + '/Frames.html', 'w')
  html_file.write ('<html>\n')
  html_file.write ('<FRAMESET ROWS=50%,50%>\n')
  html_file.write ('<FRAME src=Plots.html></FRAME>\n')
  html_file.write ('<FRAME src=Plots.html></FRAME>\n')
  html_file.write ('</FRAMESET>\n')
  html_file.write ('</html>')
  html_file.close

  # Set-up the html file
  html_file = open (plot_dir + '/Plots.html', 'w')
  html_file.write ('<html>\n')
  html_file.write ('<h1>Assimilation output</h1>\n')
  html_file.write (data_dir + '/Ens')
else:
  html_file = 0

# Get RMS values, error values of control analysis over cropped domain (plot_diff_fields_custom)
OuterLoops = 1
Anal_traj   = 'LS_Oloop%03i_Iloop000.nc'% (OuterLoops+1)

if PlotAE_traj:
  if (output_type == 'web'):
    html_file.write ('<h2>Analysis error trajectory</h2>\n')
  times, QuantityNames, TRUTH_mean, ANAL_mean, ANAL_err, TRUTH_RMS, ANAL_RMS, ANAL_RMSE = \
        plot_diff_fields_custom (Truth_dir + '/' + Truth_traj,
                          data_dir + '/' + Anal_traj,
                          'Ens_ctrl_err', output_type, html_file, plot_dir,
                          C_param, f, TimeOutput, bound_tuple)

  # Dump the time sequence of error data to a file inside the data directory
  dump_scalar_time_seq (times, QuantityNames, 'analctrl', TRUTH_mean, ANAL_mean, ANAL_err, TRUTH_RMS, ANAL_RMS, ANAL_RMSE, data_dir + '/Ens')


# Get RMS values, error values for each ensemble member over cropped domain (plot_diff_fields_custom)
for ens in range(int(nens)):

  Anal_traj   = 'ABC_bg_Ens'+str(ens+1).zfill(3)+'.nc'

  # Plot the analysis error trajectory fields
  # --------------------------------------------
  # Also return the mean analysis error, and RMS analysis error for all quantities
  if PlotAE_traj:
    if (output_type == 'web'):
      html_file.write ('<h2>Analysis error trajectory</h2>\n')
    times, QuantityNames, TRUTH_mean, ANAL_mean, ANAL_err, TRUTH_RMS, ANAL_RMS, ANAL_RMSE = \
          plot_diff_fields_custom (Truth_dir + '/' + Truth_traj,
                            data_dir + '/Ens/' + Anal_traj,
                            'Ens_err', output_type, html_file, plot_dir,
                            C_param, f, TimeOutput, bound_tuple)

    # Dump the time sequence of error data to a file inside the data directory
    dump_scalar_time_seq (times, QuantityNames, 'anal'+str(ens+1).zfill(3), TRUTH_mean, ANAL_mean, ANAL_err, TRUTH_RMS, ANAL_RMS, ANAL_RMSE, data_dir + '/Ens')

    # Plot the scalar means
    if TimeOutput  != "None":
      plot_scalar_time_seq (times, [0, times[-1]], QuantityNames, 'anal'+str(ens+1).zfill(3), TRUTH_mean, ANAL_mean, ANAL_err, TRUTH_RMS, ANAL_RMS, ANAL_RMSE, output_type, html_file, plot_dir)

if (output_type == 'web'):
  html_file.write ('</html>')
  html_file.close
