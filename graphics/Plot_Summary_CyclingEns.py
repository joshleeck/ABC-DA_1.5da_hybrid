#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read the summary files from an experiment and
# plot the ensemble trajectories, mean and truth
#
# Joshua Lee, October 2021
# -------------------------------------------------------------------


def PlotErrs (times_bound, times, data, nlines, quantity_short, quantity_long, data_type, units, LineStyles, LineColours, ylow, yhigh, Labels, PlotTitle, MainDir, bg_data, truth_data, ensmean_data):
  import matplotlib.pyplot as plt
  from matplotlib import colors, cm
  import matplotlib

  # Plot the errors for all experiments as a time sequence
  # data_type is bg or anal
  fig, ax = plt.subplots()
  ax.set_xlabel('time (h)', fontsize=20)
  ax.set_ylabel('Cropped Domain Mean (' + units + ')', fontsize=20)
  plt.xticks(fontsize=20)
  plt.yticks(fontsize=20)
  fig.set_size_inches(14.0, 7.0)
  for cycle in times_bound:
    plt.axvline(x=cycle, color='yellow')
  for experiment in range(nlines):
    ax.plot(times[:], data[experiment][:],  linewidth=0.5, ls=LineStyles, color=LineColours, alpha=0.5)
    print quantity_short, Labels[experiment], np.mean(data[experiment][:])

  try:
    ax.plot(times[::7], bg_data, linewidth=1, ls='solid', color='blue', label='FreeBG')
    print quantity_short, 'FreeBG', np.mean(bg_data)
  except:
    print 'No FreeBG provided'

  try:
    ax.plot(times[:], truth_data,  linewidth=1, ls='solid', color='black', label='Truth')
    print quantity_short, 'Truth', np.mean(truth_data)
  except:
   print 'No Truth provided'

  try:
    ax.plot(times[:], ensmean_data,  linewidth=1, ls='solid', color='red', label='EnsMean')
    print quantity_short, 'EnsMean', np.mean(ensmean_data)
  except:
   print 'No EnsMean provided'

  plt.ylim((ylow,yhigh))
  ax.legend(loc='upper right', fontsize=20) #(loc='best')
  plt.title(quantity_short + ' analysis', fontsize=20)
  plt.savefig(MainDir + 'Plots/Summary_' + quantity_short + '_' + data_type + '.png', bbox_inches='tight')
  plt.close('all')



# ==========================================================================
# ==========================================================================

import numpy as np
from netCDF4 import Dataset

# Plot for analysis (anal) or background (bg)
data_type = 'anal'
Res = 1500
PlotTitle = 'Normal'

MainDir = '/scratch/leeck/ABC-DA_1.5da_leeck_364/'

freebg = True
truth = True
ensmean = True

if data_type == 'anal':
  skiplines = 10
  if (Res == 1500):
    ylow_u  = -6;     yhigh_u  = 10
    ylow_v  = 1;     yhigh_v  = 2.5
    ylow_w  = -0.2;     yhigh_w  = 0.2
    ylow_rp = -0.008; yhigh_rp = 0.008
    ylow_bp = -0.03;    yhigh_bp = 0.03


nens = 30
ExpNames   = ['Ens' + str(x+1).zfill(3) for x in range(nens)]
ExpFiles    = []


for ens in range(nens):
  ExpFiles.append (MainDir + '/EBVd/SummaryEns' + str(ens+1).zfill(3) + '.dat')

nlines      = len(ExpNames)

# Define the line characteristics
#LineColours = ['red', 'magenta', 'purple', 'blue']
#LineStyles  = ['solid', 'solid', 'solid', 'solid']
LineColours = 'grey'
LineStyles = 'solid'


# Get some information from the first data file
input_file  = open (ExpFiles[0], 'r')
line        = input_file.readline()
line        = input_file.readline()
line        = input_file.readline()
line        = input_file.readline()
line        = input_file.readline()
# Get the cycle boundary times
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

if freebg: 
  input_file = open (ExpFiles[0], 'r')
  for l in range(36):
    line = input_file.readline()
  # Read-in data for u
  for l in range(4):
    line = input_file.readline()
  bg_u = map(float, line[7:].split())
  # Read-in data for v
  for l in range(4):
    line = input_file.readline()
  bg_v = map(float, line[7:].split())
  # Read-in data for w
  for l in range(4):
    line = input_file.readline()
  bg_w = map(float, line[7:].split())
  # Read-in data for rp
  for l in range(4):
    line = input_file.readline()
  bg_rp = map(float, line[7:].split())
  # Read-in data for bp
  for l in range(4):
    line = input_file.readline()
  bg_bp = map(float, line[7:].split())
  input_file.close()

if truth:
  input_file = open (ExpFiles[0], 'r')
  for l in range(63):
    line = input_file.readline()
  # Read-in data for u
  for l in range(3):
    line = input_file.readline()
  truth_u = map(float, line[7:].split())
  # Read-in data for v
  for l in range(3):
    line = input_file.readline()
  truth_v = map(float, line[7:].split())
  # Read-in data for w
  for l in range(3):
    line = input_file.readline()
  truth_w = map(float, line[7:].split())
  # Read-in data for rp
  for l in range(3):
    line = input_file.readline()
  truth_rp = map(float, line[7:].split())
  # Read-in data for bp
  for l in range(3):
    line = input_file.readline()
  truth_bp = map(float, line[7:].split())
  input_file.close()

if ensmean:
  input_file = open (ExpFiles[0], 'r')
  for l in range(83):
    line = input_file.readline()
  # Read-in data for u
  for l in range(3):
    line = input_file.readline()
  ensmean_u = map(float, line[7:].split())
  # Read-in data for v
  for l in range(3):
    line = input_file.readline()
  ensmean_v = map(float, line[7:].split())
  # Read-in data for w
  for l in range(3):
    line = input_file.readline()
  ensmean_w = map(float, line[7:].split())
  # Read-in data for rp
  for l in range(3):
    line = input_file.readline()
  ensmean_rp = map(float, line[7:].split())
  # Read-in data for bp
  for l in range(3):
    line = input_file.readline()
  ensmean_bp = map(float, line[7:].split())
  input_file.close()

# Plot the u errors
PlotErrs (times_bound, times_data, data_u, nlines, 'u', 'Zonal wind', data_type, 'm/s', LineStyles, LineColours, ylow_u, yhigh_u, ExpNames, PlotTitle, MainDir, bg_u, truth_u, ensmean_u)

# Plot the v errors
PlotErrs (times_bound, times_data, data_v, nlines, 'v', 'Meridonal wind', data_type, 'm/s', LineStyles, LineColours, ylow_v, yhigh_v, ExpNames, PlotTitle, MainDir, bg_v, truth_v, ensmean_v)

# Plot the w errors
PlotErrs (times_bound, times_data, data_w, nlines, 'w', 'Vertical wind', data_type, 'm/s', LineStyles, LineColours, ylow_w, yhigh_w, ExpNames, PlotTitle, MainDir, bg_w, truth_w, ensmean_w)

# Plot the rp errors
PlotErrs (times_bound, times_data, data_rp, nlines, 'r_prime', 'Scaled density', data_type, 'dimensionless', LineStyles, LineColours, ylow_rp, yhigh_rp, ExpNames, PlotTitle, MainDir, bg_rp, truth_rp, ensmean_rp)

# Plot the bp errors
PlotErrs (times_bound, times_data, data_bp, nlines, 'b_prime', 'Buoyancy', data_type, 'm/s$^2$', LineStyles, LineColours, ylow_bp, yhigh_bp, ExpNames, PlotTitle, MainDir, bg_bp, truth_bp, ensmean_bp)
