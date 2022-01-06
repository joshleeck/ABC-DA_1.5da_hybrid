#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read spread and error diagnostics from each cycle and plot them
#
# Joshua Lee, October 2021
# -------------------------------------------------------------------


import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import matplotlib
from Routines4PlotAssimDiags import plot_scalar_time_seq, plot_scalar_time_seq_custom
import os
import sys
from ast import literal_eval

# The base directory of the cycling will be specified by the command line argument if present
if len(sys.argv) > 1:
  print ('Base directory specified on the command line')
  Base_dir  = sys.argv[1]
else:
  print ('Base directory specified in the python script')
  Base_dir  = '/home/data'
print ('Base_dir = ', Base_dir)

# Bound for taking spatial mean over
try:
  bound_tuple = sys.argv[2]
  z1, z2, x1, x2 = literal_eval(bound_tuple)
except:
  print "Need to input proper bounds for taking spatial mean"
  sys.exit()


# The file containing the list of sequential DA runs performed by the script
Exp_list    = 'ExpList.dat'

# Also plot data from the ensemble spread?
Plot_ensspr = True

# The location of the plots
plot_dir    = Base_dir

# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'web'


os.system('mkdir -p ' + plot_dir + '/Plots')
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
  html_file.write ('<h1>Multi-cycle assimilation error time sequences</h1>\n')
  html_file.write (Base_dir)


# Read-in the list of files
input_file = open (Base_dir + '/' + Exp_list, 'r')
cycle_list = []
line = input_file.readline()
line = input_file.readline()
while (line != ''):
  cycle_list.append(line[:-1])  # Remove newline character
  line = input_file.readline()
input_file.close()
print ('There are ', len(cycle_list), ' cycles in this DA run')


# Get the analysis error data from the DA cycles (to compare with spread)
data_type = 'analctrl'


print ('=========================')
print ('Dealing with ', data_type, 'errors')


if (output_type == 'web'):
  html_file.write ('<h2>Plots of ' + data_type + 'errors</h2>\n')

# Analysis data
print ('  Reading from files')

# The master time sequencies for all cycles concatenated together
data_mean       = []    # Many times per cycle
data_mean_err   = []    # Many times per cycle
data_rms        = []    # Many times per cycle
data_rms_err    = []    # Many times per cycle

ensspr_mean     = []    # Mant times per cycle

cycle_count     = -1

for cycle_dir in cycle_list:
    print ('  Dealing with directory : ', cycle_dir)
    cycle_count += 1
    input_file = open (cycle_dir + '/Ens/' + data_type + '.dat', 'r')
    line = input_file.readline()
    line = input_file.readline()
    # Get the times
    line = input_file.readline()
    line = input_file.readline()
    times_cycle = map(float, line.split())
    ntimes_eachcycle = len(times_cycle)
    #print ('    times_cycle:', times_cycle)

    # Get quantities that have mean and rms values and errors
    # Put into a generic structure first to discover what is available
    quantities          = []
    data_mean_cycle     = []
    data_mean_err_cycle = []
    data_rms_cycle      = []
    data_rms_err_cycle  = []

    for quantity in range(6):
      line = input_file.readline()
      quantities.append(line[:-1])  # -1 to remove newline
      line = input_file.readline()
      line = input_file.readline()
      data_mean_cycle.append(map(float, line.split()))
      line = input_file.readline()
      data_mean_err_cycle.append(map(float, line.split()))
      line = input_file.readline()
      line = input_file.readline()
      data_rms_cycle.append(map(float, line.split()))
      line = input_file.readline()
      data_rms_err_cycle.append(map(float, line.split())) 

    input_file.close()


    # Create the master structures if this is the first cycle through
    # This creates lists of empty time sequences
    if len(data_mean) == 0:
      for quantity in range(len(data_mean_cycle)):
        data_mean.append([])
        data_mean_err.append([])
    if len(data_rms) == 0:
      for quantity in range(len(data_rms_cycle)):
        data_rms.append([])
        data_rms_err.append([])
      if (Plot_ensspr):
        for quantity in range(len(data_rms_cycle)):
          ensspr_mean.append([])
     # Append these data to the master time sequencies
    for quantity in range(len(data_mean_cycle)):
      for time in data_mean_cycle[quantity]:
        data_mean[quantity].append(time)
      for time in data_mean_err_cycle[quantity]:
        data_mean_err[quantity].append(time)
    for quantity in range(len(data_rms_cycle)):
      for time in data_rms_cycle[quantity]:
        data_rms[quantity].append(time)
      for time in data_rms_err_cycle[quantity]:
        data_rms_err[quantity].append(time)

    if (Plot_ensspr):
      # Ensemble variance file
      ensvar_file = cycle_dir + '/Ens/ABC_bg_ensvar.nc'
      nc_file_ensvar = Dataset(ensvar_file)

      for quantity in range(6):
        n_times = nc_file_ensvar.variables[quantities[quantity]].shape[0]

        for time in range(n_times):
          ensvar = nc_file_ensvar.variables[quantities[quantity]][time,z1:z2,x1:x2]
          # Mean of ensemble spread at this time
          ensspr_mean[quantity].append(np.sqrt(np.mean(ensvar)))

      nc_file_ensvar.close


print ('Storing the analysis data')
anal_mean_err = np.asarray(data_mean_err)
anal_rms_err  = np.asarray(data_rms_err)

# Set-up time data - remember that the end time of one cycle is the start time of the next
times = []
dt    = (times_cycle[1]  - times_cycle[0]) / 3600.0
Dt    = dt * float(ntimes_eachcycle - 1)
print ('Data assimilation time step is ', dt)
print ('Cycle length is                ', Dt)
n_cycles = len(cycle_list)
print ('There are ', n_cycles, ' cycles')
print ('There are ', ntimes_eachcycle, ' data assimilation steps per cycle')

for cycle in range(n_cycles):
  for step in range(ntimes_eachcycle):
    times.append(float(cycle)*Dt + float(step)*dt)
#print ('The times are')
#print (times)

# Set-up the cycle boundaries
cycle_bound_times = []
for cycle in range(n_cycles+1):
  cycle_bound_times.append(float(cycle)*Dt)
    

if (output_type == 'web'):
  html_file.write ('</html>')
  html_file.close


# Finally, write-out all data for possible analysis in another program
# --------------------------------------------------------------------
summary_file = open (plot_dir + '/SummaryEnsSpread.dat', 'w')
summary_file.write ('=====================================================================\n')
summary_file.write ('This file contains analysis errors, and ensemble spread (domain mean) data\n')
summary_file.write (plot_dir + '\n')
summary_file.write ('=====================================================================\n')
summary_file.write ('These are the boundary times\n')
output_string = ''
for item in cycle_bound_times:
  output_string += str(item) + ' '
summary_file.write (output_string + '\n')
summary_file.write ('These are the times of the inter-cycle data\n')
output_string = ''
for item in times:
  output_string += str(item) + ' '
summary_file.write (output_string + '\n')

summary_file.write ('=====================================================================\n')
summary_file.write ('Data for the analysis errors\n')
for quantity in range(6):
  summary_file.write ('---------------------------------------------------------------------\n')
  summary_file.write ('Quantity   : ' + quantities[quantity] + '\n')
  output_string = ''
  for item in anal_mean_err[quantity]:
    output_string += str(item) + ' '
  summary_file.write ('mean error : ' + output_string + '\n')
  output_string = ''
  for item in data_rms_err[quantity]:
    output_string += str(item) + ' '
  summary_file.write ('rms  : ' + output_string + '\n')

summary_file.write ('=====================================================================\n')
summary_file.write ('Data for the ensemble spread (domain mean)\n')
for quantity in range(6):
  summary_file.write ('---------------------------------------------------------------------\n')
  summary_file.write ('Quantity   : ' + quantities[quantity] + '\n')
  output_string = ''
  for item in ensspr_mean[quantity]:
    output_string += str(item) + ' '
  summary_file.write ('mean : ' + output_string + '\n')

summary_file.close
