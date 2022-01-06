#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read the summary files from four experiments and
# plot them together (scale-dependent errors)
#
# Ross Bannister, July 2020
# -------------------------------------------------------------------


def PlotErrs (scales, data, nExps, quantity_short, quantity_long, data_type, LineStyles, LineColours, Labels, PlotTitle, MainDir):
  import matplotlib.pyplot as plt
  from matplotlib import colors, cm
  import matplotlib

  # Plot the errors for all experiments as a function of scale
  # data_type is bg or anal
  print 'Plotting data for ', quantity_short
  fig, ax = plt.subplots()
  ax.set_xlabel('scale (m)', fontsize=28)
  ax.set_ylabel('RMSE', fontsize=28)
  plt.xticks(fontsize=24)
  plt.yticks(fontsize=24)
  fig.set_size_inches(14.0, 7.0)
  ax.set_xscale('log')
  ax.set_yscale('log')
  for experiment in range(nExps):
    ax.plot(scales[:], data[experiment][:],  linewidth=2, ls=LineStyles[experiment], color=LineColours[experiment], label=Labels[experiment])
  ax.legend(loc='best', fontsize=24) #(loc='upper right')
  plt.title(quantity_long + ' ' + quantity_short + ' (' + data_type + ') for ' + PlotTitle, fontsize=28)
  plt.savefig(MainDir + '/Summary_spectral_' + quantity_short + '_' + data_type + '.eps', bbox_inches='tight')
  plt.close('all')



# ==========================================================================
# ==========================================================================

import numpy as np
from netCDF4 import Dataset

PlotTitle   = 'hi-res (1.5km grid)'
#PlotTitle   = 'very hi-res (150m grid)'
MainDir     = '/media/ross/1297-5336/MeteororologyWork/DataAssim/RuthsModel/StdRes/Assim/AtmosParams/WithTruth_011'

# Plot for analysis (anal) or background (bg)
data_type = 'anal'

ExpNames    = ['alpha=0, beta=0', 'alpha=0, beta=1', 'alpha=1, beta=0', 'alpha=1, beta=1']
ExpFiles    = []
ExpFiles.append (MainDir + '/Exp-GB-HB-AB/Assim_Obs_vrhob/')
ExpFiles.append (MainDir + '/Exp-GB+HB-AB/Assim_Obs_vrhob/')
ExpFiles.append (MainDir + '/Exp+GB-HB-AB/Assim_Obs_vrhob/')
ExpFiles.append (MainDir + '/Exp+GB+HB-AB/Assim_Obs_vrhob/')
VarNames    = ['u', 'v', 'w', 'rp', 'bp']
nExps       = len(ExpFiles)

# Define the line characteristics
LineColours = ['black', 'red', 'blue', 'purple']
LineStyles  = ['solid', 'dashed', 'dotted', 'dashdot']



# Get some information from the first data file
input_file  = open (ExpFiles[0] + data_type + '_err_u_spec.dat', 'r')
line        = input_file.readline()
# Get the times
line        = input_file.readline()
times       = map(float, line.split())
ntimes      = len(times)
line        = input_file.readline()
line        = input_file.readline()
# Get the scales
line        = input_file.readline()
line        = input_file.readline()
scales      = map(float, line.split())
nscales     = len(scales)
input_file.close()

choosetime  = ntimes / 2


data_u  = np.zeros((nExps,nscales))
data_v  = np.zeros((nExps,nscales))
data_w  = np.zeros((nExps,nscales))
data_rp = np.zeros((nExps,nscales))
data_bp = np.zeros((nExps,nscales))


# Read the data
for experiment in range(nExps):

  # Read-in data for u for this experiment
  print 'Reading data for u'
  input_file = open (ExpFiles[experiment] + data_type + '_err_u_spec.dat', 'r')
  # Read the blank lines
  for l in range(7):
    line = input_file.readline()
  # Loop over each scale
  for scale in range(nscales):
    line       = input_file.readline()
    line       = input_file.readline()
    line_float = map(float, line.split())
    # Extract data for the particular time that is of interest
    data_u[experiment,scale] = line_float[choosetime-1]
  input_file.close()

  # Read-in data for v for this experiment
  print 'Reading data for v'
  input_file = open (ExpFiles[experiment] + data_type + '_err_v_spec.dat', 'r')
  # Read the blank lines
  for l in range(7):
    line = input_file.readline()
  # Loop over each scale
  for scale in range(nscales):
    line       = input_file.readline()
    line       = input_file.readline()
    line_float = map(float, line.split())
    # Extract data for the particular time that is of interest
    data_v[experiment,scale] = line_float[choosetime-1]
  input_file.close()

  # Read-in data for w for this experiment
  print 'Reading data for w'
  input_file = open (ExpFiles[experiment] + data_type + '_err_w_spec.dat', 'r')
  # Read the blank lines
  for l in range(7):
    line = input_file.readline()
  # Loop over each scale
  for scale in range(nscales):
    line       = input_file.readline()
    line       = input_file.readline()
    line_float = map(float, line.split())
    # Extract data for the particular time that is of interest
    data_w[experiment,scale] = line_float[choosetime-1]
  input_file.close()

  # Read-in data for rp for this experiment
  print 'Reading data for rp'
  input_file = open (ExpFiles[experiment] + data_type + '_err_rp_spec.dat', 'r')
  # Read the blank lines
  for l in range(7):
    line = input_file.readline()
  # Loop over each scale
  for scale in range(nscales):
    line       = input_file.readline()
    line       = input_file.readline()
    line_float = map(float, line.split())
    # Extract data for the particular time that is of interest
    data_rp[experiment,scale] = line_float[choosetime-1]
  input_file.close()
  
  # Read-in data for bp for this experiment
  print 'Reading data for bp'
  input_file = open (ExpFiles[experiment] + data_type + '_err_bp_spec.dat', 'r')
  # Read the blank lines
  for l in range(7):
    line = input_file.readline()
  # Loop over each scale
  for scale in range(nscales):
    line       = input_file.readline()
    line       = input_file.readline()
    line_float = map(float, line.split())
    # Extract data for the particular time that is of interest
    data_bp[experiment,scale] = line_float[choosetime-1]
  input_file.close()


# Plot the u errors
PlotErrs (scales, data_u, nExps, 'u', 'Zonal wind', data_type, LineStyles, LineColours, ExpNames, PlotTitle, MainDir)

# Plot the v errors
PlotErrs (scales, data_v, nExps, 'v', 'Meridional wind', data_type, LineStyles, LineColours, ExpNames, PlotTitle, MainDir)

# Plot the w errors
PlotErrs (scales, data_w, nExps, 'w', 'Vertical wind', data_type, LineStyles, LineColours, ExpNames, PlotTitle, MainDir)

# Plot the rp errors
PlotErrs (scales, data_rp, nExps, 'rp', 'Scaled density', data_type, LineStyles, LineColours, ExpNames, PlotTitle, MainDir)

# Plot the bp errors
PlotErrs (scales, data_bp, nExps, 'bp', 'Buoyancy', data_type, LineStyles, LineColours, ExpNames, PlotTitle, MainDir)
