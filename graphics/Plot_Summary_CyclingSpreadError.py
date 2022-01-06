#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read the summary file from an experiment and
# plot the ensemble spread, RMSE and B_c implied covariances
#
# Joshua Lee, October 2021
# -------------------------------------------------------------------


def PlotErrs (times_bound, times, data, nlines, quantity_short, quantity_long, data_type, units, LineStyles, LineColours, ylow, yhigh, PlotTitle, MainDir, ensspr_data, cov_c_data):
  import matplotlib.pyplot as plt
  from matplotlib import colors, cm
  import matplotlib

  # Plot the errors for all experiments as a time sequence
  # data_type is bg or anal
  fig, ax = plt.subplots()
  ax.set_xlabel('time (h)', fontsize=20)
  ax.set_ylabel('RMSE/Spread/S.D. (' + units + ')', fontsize=20)
  plt.xticks(fontsize=20)
  plt.yticks(fontsize=20)
  fig.set_size_inches(14.0, 7.0)
  for cycle in times_bound:
    plt.axvline(x=cycle, color='yellow')
  for experiment in range(nlines):
    ax.plot(times[:], data[experiment][:],  linewidth=1, ls='solid', color='black', label='RMSE')
    print quantity_short, 'RMSE', np.mean(data[experiment][:])

  try:
    ax.plot(times[:], ensspr_data,  linewidth=1, ls=LineStyles, color=LineColours, label='Spread')
    print quantity_short, 'EnsSpr', np.mean(ensspr_data)
  except:
    print 'No EnsSpr data provided'

  try:
    ax.plot(times[:], cov_c_data,  linewidth=1, ls='solid', color='blue', label='Bc S.D.')
    print quantity_short, 'B_c', np.mean(cov_c_data)
  except:
    print 'No B_c data provided'

  plt.ylim((ylow,yhigh))
  ax.legend(loc='upper right', fontsize=20) #(loc='upper right') #'best'
  plt.title(quantity_short + ' analysis', fontsize=20)
  plt.savefig(MainDir + 'Plots/Summary_' + quantity_short + '_' + data_type + '.png', bbox_inches='tight')
  plt.close('all')



# ==========================================================================
# ==========================================================================

import numpy as np
from netCDF4 import Dataset

# Plot for analysis (anal) or background (bg)
data_type = 'rmse_spr'
Res = 1500
PlotTitle = 'Normal'

MainDir = '/scratch/leeck/ABC-DA_1.5da_leeck_364/'

ensspr = True
cov_c = True

if data_type == 'rmse_spr':
  skiplines = 10
  if (Res == 1500):
    ylow_u  = 0;     yhigh_u  = 3.5
    ylow_v  = 0;     yhigh_v  = 0.6
    ylow_w  = 0.0;     yhigh_w  = 0.08
    ylow_rp = 0.0; yhigh_rp = 0.003
    ylow_bp = 0.0;    yhigh_bp = 0.03



ExpFiles    = []

ExpFiles.append (MainDir + '/EnSRFd/SummaryEnsSpread.dat')

impliedcovfile = '/scratch/leeck/ABC-DA_1.5da_leeck_364/examples/Master_ImpliedCov/ImpliedCov.dat'

nlines      = len(ExpFiles)

# Define the line characteristics
#LineColours = ['red', 'magenta', 'purple', 'blue']
#LineStyles  = ['solid', 'solid', 'solid', 'solid']
LineColours = 'red'
LineStyles = 'solid'


# Get some information from the first data file
input_file  = open (ExpFiles[0], 'r')
line        = input_file.readline()
line        = input_file.readline()
line        = input_file.readline()
line        = input_file.readline()
line        = input_file.readline()
# Get the cycle bounary times
line        = input_file.readline()
times_bound = map(float, line.split())
line        = input_file.readline()
# Get the data times
line        = input_file.readline()
times_data  = map(float, line.split())
input_file.close()


data_u  = []
data_v  = []
data_w  = []
data_rp = []
data_bp = []
# Read the data
for experiment in range(nlines):
  input_file = open (ExpFiles[experiment], 'r')

  # Read the blank lines
  for l in range(skiplines):
    line = input_file.readline()

  # Read-in data for u
  for l in range(4):
    line = input_file.readline()
  data_u.append(map(float, line[7:].split()))
  # Read-in data for v
  for l in range(4):
    line = input_file.readline()
  data_v.append(map(float, line[7:].split()))
  # Read-in data for w
  for l in range(4):
    line = input_file.readline()
  data_w.append(map(float, line[7:].split()))
  # Read-in data for rp
  for l in range(4):
    line = input_file.readline()
  data_rp.append(map(float, line[7:].split()))
  # Read-in data for bp
  for l in range(4):
    line = input_file.readline()
  data_bp.append(map(float, line[7:].split()))
  input_file.close()

if ensspr:
  input_file = open (ExpFiles[0], 'r')
  for l in range(36):
    line = input_file.readline()
  # Read-in data for u
  for l in range(3):
    line = input_file.readline()
  ensspr_u = map(float, line[7:].split())
  # Read-in data for v
  for l in range(3):
    line = input_file.readline()
  ensspr_v = map(float, line[7:].split())
  # Read-in data for w
  for l in range(3):
    line = input_file.readline()
  ensspr_w = map(float, line[7:].split())
  # Read-in data for rp
  for l in range(3):
    line = input_file.readline()
  ensspr_rp = map(float, line[7:].split())
  # Read-in data for bp
  for l in range(3):
    line = input_file.readline()
  ensspr_bp = map(float, line[7:].split())
  input_file.close()

if cov_c:
  input_file = open (impliedcovfile, 'r')
  for l in range(4):
    line = input_file.readline()
  cov_c_u = [float(line)]*len(times_data)
  for l in range(4):
    line = input_file.readline()
  cov_c_v = [float(line)]*len(times_data)  
  for l in range(4):
    line = input_file.readline()
  cov_c_w = [float(line)]*len(times_data)
  for l in range(4):
    line = input_file.readline()
  cov_c_rp = [float(line)]*len(times_data)
  for l in range(4):
    line = input_file.readline()
  cov_c_bp = [float(line)]*len(times_data)
  input_file.close()


# Plot the u errors
PlotErrs (times_bound, times_data, data_u, nlines, 'u', 'Zonal wind', data_type, 'm/s', LineStyles, LineColours, ylow_u, yhigh_u, PlotTitle, MainDir, ensspr_u, cov_c_u)

# Plot the v errors
PlotErrs (times_bound, times_data, data_v, nlines, 'v', 'Meridonal wind', data_type, 'm/s', LineStyles, LineColours, ylow_v, yhigh_v, PlotTitle, MainDir, ensspr_v, cov_c_v)

# Plot the w errors
PlotErrs (times_bound, times_data, data_w, nlines, 'w', 'Vertical wind', data_type, 'm/s', LineStyles, LineColours, ylow_w, yhigh_w, PlotTitle, MainDir, ensspr_w, cov_c_w)

# Plot the rp errors
PlotErrs (times_bound, times_data, data_rp, nlines, 'r_prime', 'Scaled density', data_type, 'dimensionless', LineStyles, LineColours, ylow_rp, yhigh_rp, PlotTitle, MainDir, ensspr_rp, cov_c_rp)

# Plot the bp errors
PlotErrs (times_bound, times_data, data_bp, nlines, 'b_prime', 'Buoyancy', data_type, 'm/s$^2$', LineStyles, LineColours, ylow_bp, yhigh_bp, PlotTitle, MainDir, ensspr_bp, cov_c_bp)
