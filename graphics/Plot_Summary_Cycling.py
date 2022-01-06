#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read the summary files from experiments and
# plot them together
#
# Ross Bannister, July 2020
# Adapted Joshua Lee, October 2021
# -------------------------------------------------------------------


def PlotErrs (times_bound, times, data, nlines, quantity_short, quantity_long, data_type, units, LineStyles, LineColours, ylow, yhigh, Labels, PlotTitle, MainDir, bg):
  import matplotlib.pyplot as plt
  from matplotlib import colors, cm
  import matplotlib

  # Plot the errors for all experiments as a time sequence
  # data_type is bg or anal
  fig, ax = plt.subplots()
  ax.set_xlabel('time (h)', fontsize=20)
  ax.set_ylabel('RMSE (' + units + ')', fontsize=20)
  plt.xticks(fontsize=20)
  plt.yticks(fontsize=20)
  fig.set_size_inches(14.0, 7.0)
  for cycle in times_bound:
    plt.axvline(x=cycle, color='yellow')
  for experiment in range(nlines):
    ax.plot(times[:], data[experiment][:],  linewidth=1, ls=LineStyles[experiment], color=LineColours[experiment], label=Labels[experiment])
    print quantity_short, Labels[experiment], np.mean(data[experiment][:])

  ax.plot(times[::7],bg, linewidth=1, ls='solid', color='black', label='FreeBG')
  print quantity_short, 'FreeBG', np.mean(bg)
  plt.ylim((ylow,yhigh))
  ax.legend(loc='upper right', fontsize=20) #(loc='best')
  plt.title(quantity_short + ' analysis errors', fontsize=20)
  plt.savefig(MainDir + 'Plots/Summary_' + quantity_short + '_' + data_type + '.png', bbox_inches='tight')
  plt.close('all')



# ==========================================================================
# ==========================================================================

import numpy as np
from netCDF4 import Dataset

# Set resolution of the model (150 or 1500) m grid
Res = 1500
# Plot for analysis (anal) or background (bg)
data_type = 'anal'

if (Res == 1500):
  # This is standard resolution (called high resolution in paper)
  PlotTitle = 'Normal'

MainDir = '/scratch/leeck/ABC-DA_1.5da_leeck_364/'

freebg = True

if data_type == 'anal':
  skiplines = 36
  if (Res == 1500):
    ylow_u  = 0.2;     yhigh_u  = 1.4
    ylow_v  = 0.2;     yhigh_v  = 1.3
    ylow_w  = 0.0;     yhigh_w  = 0.05
    ylow_rp = 0.0002; yhigh_rp = 0.001
    ylow_bp = 0.002;    yhigh_bp = 0.007
  else:
    ylow_u  = 1.5;     yhigh_u  = 3.5
    ylow_v  = 1.7;     yhigh_v  = 2.8
    ylow_w  = 0.0;     yhigh_w  = 0.05
    ylow_rp = 0.0005; yhigh_rp = 0.0009
    ylow_bp = 0.01;    yhigh_bp = 0.03
else:
  skiplines = 10
  if (Res == 1500):
    ylow_u  = 1.5;     yhigh_u  = 3.5
    ylow_v  = 1.8;     yhigh_v  = 4.0
    ylow_w  = 0.0;     yhigh_w  = 0.6
    ylow_rp = 0.00005; yhigh_rp = 0.0005
    ylow_bp = 0.01;    yhigh_bp = 0.05
  else:
    ylow_u  = 1.5;     yhigh_u  = 3.5
    ylow_v  = 1.7;     yhigh_v  = 6.0
    ylow_w  = 0.0;     yhigh_w  = 0.6
    ylow_rp = 0.00005; yhigh_rp = 0.0008
    ylow_bp = 0.01;    yhigh_bp = 0.06

#ExpNames    = ['100% Be', '80% Be, 20% Bc', '50% Be, 50% Be', '100% Bc']
#ExpNames     = ['100% Be-EBV', '100% Be-EnSRF', '100% Bc']
ExpNames     = ['100% Be-30mem', '100% Be-20mem', '100% Be-10mem', '100% Bc']
ExpFiles    = []
ExpFiles.append (MainDir + '/EBVd/Summary.dat')
#ExpFiles.append (MainDir + '/EnSRFd/Summary.dat')
#ExpFiles.append (MainDir + '/EnSRFc/Summary.dat')
#ExpFiles.append (MainDir + '/EnSRFb/Summary.dat')
ExpFiles.append (MainDir + '/EBVd-20mem/Summary.dat')
ExpFiles.append (MainDir + '/EBVd-10mem/Summary.dat')
ExpFiles.append (MainDir + '/3DVAR/Summary.dat')

nlines      = len(ExpNames)

# Define the line characteristics
#LineColours = ['red', 'magenta', 'purple', 'blue']
LineColours = ['red', 'cyan', 'brown', 'blue'] 
LineStyles  = ['solid', 'solid', 'solid', 'solid']
#LineColours = ['green','orange','blue']
#LineStyles = ['solid','solid','solid']


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
  data_u.append(map(float, line[12:].split()))
  # Read-in data for v
  for l in range(4):
    line = input_file.readline()
  data_v.append(map(float, line[12:].split()))
  # Read-in data for w
  for l in range(4):
    line = input_file.readline()
  data_w.append(map(float, line[12:].split()))
  # Read-in data for rp
  for l in range(4):
    line = input_file.readline()
  data_rp.append(map(float, line[12:].split()))
  # Read-in data for bp
  for l in range(4):
    line = input_file.readline()
  data_bp.append(map(float, line[12:].split()))
  input_file.close()

if freebg: 
  input_file = open (ExpFiles[0], 'r')
  for l in range(62):
    line = input_file.readline()
  # Read-in data for u
  for l in range(4):
    line = input_file.readline()
  bg_u = map(float, line[12:].split())
  # Read-in data for v
  for l in range(4):
    line = input_file.readline()
  bg_v = map(float, line[12:].split())
  # Read-in data for w
  for l in range(4):
    line = input_file.readline()
  bg_w = map(float, line[12:].split())
  # Read-in data for rp
  for l in range(4):
    line = input_file.readline()
  bg_rp = map(float, line[12:].split())
  # Read-in data for bp
  for l in range(4):
    line = input_file.readline()
  bg_bp = map(float, line[12:].split())
  input_file.close()


# Plot the u errors
PlotErrs (times_bound, times_data, data_u, nlines, 'u', 'Zonal wind', data_type, 'm/s', LineStyles, LineColours, ylow_u, yhigh_u, ExpNames, PlotTitle, MainDir, bg_u)

# Plot the v errors
PlotErrs (times_bound, times_data, data_v, nlines, 'v', 'Meridonal wind', data_type, 'm/s', LineStyles, LineColours, ylow_v, yhigh_v, ExpNames, PlotTitle, MainDir, bg_v)

# Plot the w errors
PlotErrs (times_bound, times_data, data_w, nlines, 'w', 'Vertical wind', data_type, 'm/s', LineStyles, LineColours, ylow_w, yhigh_w, ExpNames, PlotTitle, MainDir, bg_w)

# Plot the rp errors
PlotErrs (times_bound, times_data, data_rp, nlines, 'r_prime', 'Scaled density', data_type, 'dimensionless', LineStyles, LineColours, ylow_rp, yhigh_rp, ExpNames, PlotTitle, MainDir, bg_rp)

# Plot the bp errors
PlotErrs (times_bound, times_data, data_bp, nlines, 'b_prime', 'Buoyancy', data_type, 'm/s$^2$', LineStyles, LineColours, ylow_bp, yhigh_bp, ExpNames, PlotTitle, MainDir, bg_bp)
