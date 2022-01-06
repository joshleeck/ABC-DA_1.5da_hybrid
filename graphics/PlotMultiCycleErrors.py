#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read error diagnostics from each cycle, splice them, and then plot
#
# Ross Bannister, August 2018
# -------------------------------------------------------------------


import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import matplotlib
from Routines4PlotAssimDiags import plot_scalar_time_seq, plot_scalar_time_seq_custom
import os
import sys

# The base directory of the cycling will be specified by the command line argument if present
if len(sys.argv) > 1:
  print ('Base directory specified on the command line')
  Base_dir  = sys.argv[1]
else:
  print ('Base directory specified in the python script')
  Base_dir  = '/home/data'
print ('Base_dir = ', Base_dir)

# The file containing the list of sequential DA runs performed by the script
Exp_list    = 'ExpList.dat'

# Also plot data from the free backgound run?
Plot_freebg = True
# The full path of the file containing the forecast from the first background (relevant only if Plot_freebg = True)
FullBg_file = Base_dir + '/Master_RunNLModel_Fullbg/BgFc.nc'

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




# If requested, open the free backgound run netcdf file
if Plot_freebg:
  nc_file_bg = Dataset(FullBg_file)


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


# Get the error data from the DA cycles
files = ['bg', 'bg_fbg_truth', 'anal', 'anal_fbg_truth', 'anal-bg']

for data_type in files:
  print ('=========================')
  print ('Dealing with ', data_type, ' data and errors')


  if (output_type == 'web'):
    html_file.write ('<h2>Plots of ' + data_type + ' data and errors</h2>\n')

  if (data_type == 'bg') or (data_type == 'anal'):
    # Either background or analysis data
    print ('  Reading from files')

    # The master time sequencies for all cycles concatenated together
    truth_mean      = []    # Many times per cycle
    data_mean       = []    # Many times per cycle
    data_mean_err   = []    # Many times per cycle

    truth_rms       = []    # Many times per cycle
    data_rms        = []    # Many times per cycle
    data_rms_err    = []    # Many times per cycle

    if data_type == 'bg':
      # Only need to initialise when dealing with background data
      freebg_mean     = []    # Once per cycle
      freebg_mean_err = []    # Once per cycle
      freebg_rms      = []    # Once per cycle
      freebg_rms_err  = []    # Once per cycle

    cycle_count     = -1

    for cycle_dir in cycle_list:
      print ('  Dealing with directory : ', cycle_dir)
      cycle_count += 1
      input_file = open (cycle_dir + '/' + data_type + '.dat', 'r')
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
      truth_mean_cycle    = []
      data_mean_cycle     = []
      data_mean_err_cycle = []
      truth_rms_cycle     = []
      data_rms_cycle      = []
      data_rms_err_cycle  = []
      for quantity in range(6):
        line = input_file.readline()
        quantities.append(line[:-1])  # -1 to remove newline
        line = input_file.readline()
        truth_mean_cycle.append(map(float, line.split()))
        line = input_file.readline()
        data_mean_cycle.append(map(float, line.split()))
        line = input_file.readline()
        data_mean_err_cycle.append(map(float, line.split()))
        line = input_file.readline()
        truth_rms_cycle.append(map(float, line.split()))
        line = input_file.readline()
        data_rms_cycle.append(map(float, line.split()))
        line = input_file.readline()
        data_rms_err_cycle.append(map(float, line.split()))

      '''
      # Get balance quantities that are to do with scale bands
      # These have only rms values (and we don't yet know how many bands there are)
      line = input_file.readline()
      while (line != ''):
        quantities.append(line[:-1])  # -1 to remove newline
        line = input_file.readline()
        truth_rms_cycle.append(map(float, line.split()))
        line = input_file.readline()
        data_rms_cycle.append(map(float, line.split()))
        line = input_file.readline()
        data_rms_err_cycle.append(map(float, line.split()))
        line = input_file.readline()
      Nscales = (len(data_rms_cycle) - len(data_mean_cycle)) / 2
      print ('  There are ', Nscales, ' scale bands in the file')
      '''

      input_file.close()


      # Get data for the free background run if requested
      if (Plot_freebg and (data_type == 'bg')):
        freebg_mean_cycle     = []
        freebg_mean_err_cycle = []
        freebg_rms_cycle      = []
        freebg_rms_err_cycle  = []

        # Truth file
        truth_file = cycle_dir + '/Obs+Truth/Truth_run.nc'
        nc_file_tr = Dataset(truth_file)

        # Read and compute for each quantity
        for quantity in range(6):
          freebg  = nc_file_bg.variables[quantities[quantity]][cycle_count,:,:]
          truth   = nc_file_tr.variables[quantities[quantity]][0,:,:]
          error   = freebg - truth
          # Mean of free background at this time
          freebg_mean_cycle.append(np.mean(freebg))
          # Mean error of free background at this time
          freebg_mean_err_cycle.append(np.mean(error))
          # RMS of free background at this time
          freebg_rms_cycle.append(np.sqrt(np.mean(freebg*freebg)))
          # RMS of free background error at this time
          freebg_rms_err_cycle.append(np.sqrt(np.mean(error*error)))

        nc_file_tr.close


      # Create the master structures if this is the first cycle through
      # This creates lists of empty time sequences
      if len(data_mean) == 0:
        for quantity in range(len(data_mean_cycle)):
          truth_mean.append([])
          data_mean.append([])
          data_mean_err.append([])
      if len(data_rms) == 0:
        for quantity in range(len(data_rms_cycle)):
          truth_rms.append([])
          data_rms.append([])
          data_rms_err.append([])
        if (Plot_freebg and (data_type == 'bg')):
          for quantity in range(6):
            freebg_mean.append([])
            freebg_mean_err.append([])
            freebg_rms.append([])
            freebg_rms_err.append([])

      # Append these data to the master time sequencies
      for quantity in range(len(data_mean_cycle)):
        for time in truth_mean_cycle[quantity]:
          truth_mean[quantity].append(time)
        for time in data_mean_cycle[quantity]:
          data_mean[quantity].append(time)
        for time in data_mean_err_cycle[quantity]:
          data_mean_err[quantity].append(time)
      for quantity in range(len(data_rms_cycle)):
        for time in truth_rms_cycle[quantity]:
          truth_rms[quantity].append(time)
        for time in data_rms_cycle[quantity]:
          data_rms[quantity].append(time)
        for time in data_rms_err_cycle[quantity]:
          data_rms_err[quantity].append(time)
      if (Plot_freebg and (data_type == 'bg')):
        for quantity in range(6):
          freebg_mean[quantity].append(freebg_mean_cycle[quantity])
          freebg_mean_err[quantity].append(freebg_mean_err_cycle[quantity])
          freebg_rms[quantity].append(freebg_rms_cycle[quantity])
          freebg_rms_err[quantity].append(freebg_rms_err_cycle[quantity])

    # Close the free background run netcdf file
    if (Plot_freebg and (data_type == 'bg')):
      nc_file_bg.close


    if data_type == 'bg':
      # Store the background data
      print ('Storing the background data')
      bg_mean_err = np.asarray(data_mean_err)
      bg_rms_err  = np.asarray(data_rms_err)
    else:
      print ('Storing the analysis data')
      anal_mean_err = np.asarray(data_mean_err)
      anal_rms_err  = np.asarray(data_rms_err)

    if data_type == 'bg':
      # For the first loop, set-up some stuff for plotting
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

    # Plot the data
    print ('Plotting for ', data_type)
    plot_scalar_time_seq (times, cycle_bound_times, quantities, data_type, truth_mean, data_mean, data_mean_err, truth_rms, data_rms, freebg_rms, data_rms_err, freebg_rms_err, output_type, html_file, plot_dir)

  elif (data_type == 'anal-bg'):
    # The third item in the loop will not read data, but find the difference between the analysis and background errors
    # Positive is then 'bad' and negative is 'good'
    print ('Taking the difference between analysis and background errors')
    data_mean_err_diff = np.absolute(anal_mean_err) - np.absolute(bg_mean_err)
    data_rms_err_diff  = np.absolute(anal_rms_err) - np.absolute(bg_rms_err)

    # Plot the difference data
    print ('Plotting for ', data_type)
    plot_scalar_time_seq (times, cycle_bound_times, quantities, data_type, truth_mean, data_mean, data_mean_err_diff, truth_rms, data_rms, freebg_rms, data_rms_err_diff, freebg_rms_err, output_type, html_file, plot_dir)
   
  elif (data_type == 'bg_fbg_truth'):
    print ('Plotting for ', data_type)
    plot_scalar_time_seq_custom (times, cycle_bound_times, quantities, data_type, truth_mean, data_mean, freebg_mean, truth_rms, data_rms, freebg_rms, data_rms_err, freebg_rms_err, output_type, html_file, plot_dir)
  elif (data_type == 'anal_fbg_truth'):
    print ('Plotting for ', data_type)
    plot_scalar_time_seq_custom (times, cycle_bound_times, quantities, data_type, truth_mean, data_mean, freebg_mean, truth_rms, data_rms, freebg_rms, data_rms_err, freebg_rms_err, output_type, html_file, plot_dir)
    



if (output_type == 'web'):
  html_file.write ('</html>')
  html_file.close


# Finally, write-out all data for possible analysis in another program
# --------------------------------------------------------------------
summary_file = open (plot_dir + '/Summary.dat', 'w')
summary_file.write ('=====================================================================\n')
summary_file.write ('This file contains analysis, background, and free bg forecast data\n')
summary_file.write (plot_dir + '\n')
summary_file.write ('=====================================================================\n')
summary_file.write ('These are the boundary times (and the times for the free bg forecast)\n')
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
summary_file.write ('Data for the background states\n')
for quantity in range(6):
  summary_file.write ('---------------------------------------------------------------------\n')
  summary_file.write ('Quantity   : ' + quantities[quantity] + '\n')
  output_string = ''
  for item in bg_mean_err[quantity]:
    output_string += str(item) + ' '
  summary_file.write ('mean error : ' + output_string + '\n')
  output_string = ''
  for item in bg_rms_err[quantity]:
    output_string += str(item) + ' '
  summary_file.write ('rms error  : ' + output_string + '\n')

summary_file.write ('=====================================================================\n')
summary_file.write ('Data for the analysis states\n')
for quantity in range(6):
  summary_file.write ('---------------------------------------------------------------------\n')
  summary_file.write ('Quantity   : ' + quantities[quantity] + '\n')
  output_string = ''
  for item in anal_mean_err[quantity]:
    output_string += str(item) + ' '
  summary_file.write ('mean error : ' + output_string + '\n')
  output_string = ''
  for item in anal_rms_err[quantity]:
    output_string += str(item) + ' '
  summary_file.write ('rms error  : ' + output_string + '\n')

summary_file.write ('=====================================================================\n')
summary_file.write ('Data for the free background\n')
for quantity in range(6):
  summary_file.write ('---------------------------------------------------------------------\n')
  summary_file.write ('Quantity   : ' + quantities[quantity] + '\n')
  output_string = ''
  for item in freebg_mean_err[quantity]:
    output_string += str(item) + ' '
  summary_file.write ('mean error : ' + output_string + '\n')
  output_string = ''
  for item in freebg_rms_err[quantity]:
    output_string += str(item) + ' '
  summary_file.write ('rms error  : ' + output_string + '\n')
summary_file.write ('=====================================================================\n')
summary_file.close
