#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read error diagnostics from each cycle, splice them, and then plot.
# This is done in spectral space in order to understand the characteristics of the
# assimilation as a function of scale.
#
# Ross Bannister, July 2020
# -------------------------------------------------------------------

# ===================================================================
def process_field (field_realspace):
  import numpy as np
  # Subroutine to take field, fourier transform it and return the level averaged spec dens

  nlevs            = field_realspace.shape[0]
  nlongs           = field_realspace.shape[1]
  nscales          = nlongs/2 + 1
  spec_dens        = np.zeros((nlongs))
  spec_dens_merged = np.zeros((nlongs/2+1))
  av_spec_dens     = np.zeros((nlongs/2+1))

  # Loop over levels to find average spectral error density
  for z in range(nlevs):
    # Do FFT of the err field
    fourier = np.fft.fft(field_realspace[z,:])
    # Find the density of the spectrum
    for element in range(nlongs):
      spec_dens[element] = abs(fourier[element])
    # Merge waves of equivalent scale
    for element in range(nlongs/2 + 1):
      if (element == 0):
        spec_dens_merged[element] = spec_dens[0]
      else:
        spec_dens_merged[element] = (spec_dens[element] + spec_dens[nlongs-element]) / 2.0

    # Add-on the contribution from this level
    av_spec_dens[:] += spec_dens_merged[:]
  # Normalise
  av_spec_dens[:] /= float(nlevs)

  return av_spec_dens

# ===================================================================
def plot_field (wavenumbers, scales, times, spectrum, quantity, data_type, plot_dir, filesuffix, html_file):

  import numpy as np
  import matplotlib.pyplot as plt
  from matplotlib import colors, cm
  import matplotlib

  # Plot the spectrum (wavenumber vs time) and also the spectrum at the mid-way timestep
  # data_type is either bg or anal

  ntimes  = len(times)
  nscales = len(scales)

  matplotlib.rc('xtick', labelsize=12)
  matplotlib.rc('ytick', labelsize=12)
  fig, ax = plt.subplots()
  cmap = cm.get_cmap('seismic', 11)
  cax  = ax.contourf(wavenumbers[0:nscales], times[0:ntimes], np.log(spectrum[0:ntimes,0:nscales]), cmap=cmap)
  cax  = ax.contour(wavenumbers[0:nscales], times[0:ntimes], np.log(spectrum[0:ntimes,0:nscales]), colors='k')
  cbar = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  # Labels etc
  ax.set_title('Error spectrum for ' + quantity, fontsize=14)
  ax.set_xlabel('wavenumber', fontsize=12)
  ax.set_ylabel('time (s)', fontsize=12)
  #plt.show()
  filename_core = data_type + '_err_' + quantity + '_spec'
  plt.savefig(plot_dir + '/Plots/' + filename_core + filesuffix, bbox_inches='tight')
  plt.close('all')


  # Plot the spectrum for one time
  matplotlib.rc('xtick', labelsize=12)
  matplotlib.rc('ytick', labelsize=12)
  fig, ax = plt.subplots()
  ax.set_xlabel('scale (m)', fontsize=12)
  ax.set_ylabel('spectral density', fontsize=12)
  ax.set_xscale('log')

  plt.title('Error spectrum for ' + quantity + ' at central time', fontsize=14)

  ax.plot(scales[1:nscales], spectrum[ntimes/2,1:nscales], linewidth=2, color='red')

  # Show the legend
  #ax.legend()

  # plt.show()             # Comment out to not display plot on the screen
  plt.savefig(plot_dir + '/Plots/' + filename_core + '1time' + filesuffix, bbox_inches='tight')  # Comment out to not save the file
  plt.close('all')


  if (filesuffix == '.png'):
    html_file.write ('<h2>Plots of ' + data_type + ' ' + quantity + ' spectral error</h2>\n')
    html_file.write ('<img src=Plots/' + filename_core  + '.png width=300></td>\n')
    html_file.write ('<img src=Plots/' + filename_core  + '1time' + '.png width=300></td>\n')



# ===================================================================
def output_field (wavenumbers, scales, times, spectrum, quantity, data_type, plot_dir):

  # Save the data for possible future use
  output_file = open (plot_dir + '/' + data_type + '_err_' + quantity + '_spec' + '.dat', 'w')
  output_file.write ('Times\n')
  for element in times:
    output_file.write (str(element) + '  ')
  output_file.write ('\nWavenumbers\n')
  for element in wavenumbers:
    output_file.write (str(element) + '  ')
  output_file.write ('\nScales\n')
  for element in scales:
    output_file.write (str(element) + '  ')
  output_file.write ('\nSpectral density data\n')
  for element in range(nscales):
    output_file.write ('Scale number ' + str(element) + '\n')
    for item in range(ntimes_total):
      output_file.write (str(spectrum[item,element]) + '  ')
    output_file.write ('\n')
  output_file.close()



# ==============================================================================
# ==============================================================================
# ==============================================================================
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import matplotlib
import os
import sys

# The base directory of the cycling will be specified by the command line argument if present
if len(sys.argv) > 1:
  print ('Base directory specified on the command line')
  Base_dir  = sys.argv[1]
else:
  print ('Base directory specified in the python script')
  Base_dir  = '/scratch2/sws98rnb/HiRes/Assim/AtmosParams/WithTruth_011/Exp-GB+HB-AB/Assim_Obs_vrhob'
print ('Base_dir = ', Base_dir)

# The number of data assimilation cycles
ncycles = 30
dx      = 1500.0


# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'web'

plot_dir = Base_dir
os.system('mkdir -p ' + plot_dir + '/Plots')

if (output_type == 'web'):
  # Set-up the html file
  html_file = open (plot_dir + '/PlotsScales.html', 'w')
  html_file.write ('<html>\n')
  html_file.write ('<h1>Multi-cycle assimilation error time sequences</h1>\n')
  html_file.write (Base_dir)
  filesuffix = '.png'
else:
  filesuffix = '.eps'


# The error data types considered over the DA cycles
data_types = ['bg', 'anal']


# Loop over background and analyis
for data_type in data_types:
  print ('=========================')
  print ('Dealing with ', data_type, ' data and errors')


  # Loop over DA cycles
  for cycle0 in range(ncycles):
    cycle = cycle0 + 1

    # Find the filename of the truth file for this cycle and open the netcdf file
    file_truth = Base_dir + '/da_cycle_%04i/Obs+Truth/Truth.nc'% (cycle)
    print 'Truth file : ', file_truth
    nc_file_true = Dataset(file_truth)

    # Find the filename of the background or analysis and open the netcdf file
    if data_type == 'bg':
      file_data = Base_dir + '/da_cycle_%04i/LS_Oloop001_Iloop000.nc'% (cycle)
    else:
      file_data = Base_dir + '/da_cycle_%04i/LS_Oloop002_Iloop000.nc'% (cycle)
    print 'Data file  : ', file_data
    nc_file_data = Dataset(file_data)

    if (data_type == 'bg') and (cycle == 1):
      # Determine the number of times
      subtimes    = nc_file_true.variables['time'][:]
      n_sub_times = len(subtimes)
      dt          = subtimes[1] - subtimes[0]
      Dt          = subtimes[n_sub_times-1] - subtimes[0]   # This is the time of the DA cycling


    # Loop over the subtimes
    for t in range(n_sub_times):

      # Get data for u (truth and background/analysis)
      true       = nc_file_true.variables['u'][t,:,:]
      data       = nc_file_data.variables['u'][t,:,:]
      err_raw_u  = data[:,:] - true[:,:]
      # Get data for v (truth and background/analysis)
      true       = nc_file_true.variables['v'][t,:,:]
      data       = nc_file_data.variables['v'][t,:,:]
      err_raw_v  = data[:,:] - true[:,:]
      # Get data for w (truth and background/analysis)
      true       = nc_file_true.variables['w'][t,:,:]
      data       = nc_file_data.variables['w'][t,:,:]
      err_raw_w  = data[:,:] - true[:,:]
      # Get data for rp (truth and background/analysis)
      true       = nc_file_true.variables['r_prime'][t,:,:]
      data       = nc_file_data.variables['r_prime'][t,:,:]
      err_raw_rp = data[:,:] - true[:,:]
      # Get data for bp (truth and background/analysis)
      true       = nc_file_true.variables['b_prime'][t,:,:]
      data       = nc_file_data.variables['b_prime'][t,:,:]
      err_raw_bp = data[:,:] - true[:,:]

      nc_file_true.close
      nc_file_data.close


      # At this point we know the size of the arrays that we are dealing with
      # So we can create the arrays needed
      if (cycle == 1) and (data_type == 'bg') and (t == 0):
        nlevs            = err_raw_u.shape[0]
        nlongs           = err_raw_u.shape[1]
        nscales          = nlongs/2 + 1
        ntimes_total     = ncycles * n_sub_times
        time             = np.zeros((ntimes_total))
        err_u            = np.zeros((ntimes_total, nscales))
        err_v            = np.zeros((ntimes_total, nscales))
        err_w            = np.zeros((ntimes_total, nscales))
        err_rp           = np.zeros((ntimes_total, nscales))
        err_bp           = np.zeros((ntimes_total, nscales))

        print 'The number of sub times ', n_sub_times
        print 'The total No times      ', ntimes_total
        print 'The number of scales    ', nscales
        print 'dt (small time step)  = ', dt
        print 'Dt (cycle time step)  = ', Dt
        print 'nlevs                 = ', nlevs
        print 'nlongs                = ', nlongs

        # Compute the wavenumbers and the equivalent scales
        wavenumber = range(nscales)
        scales     = np.zeros((nscales))
        for element in range(1,nscales):
          scales[element] = float(nlongs) * dx / float(element)
        scales[0] = scales[1]

      # Store the time
      time[cycle0*n_sub_times + t] = float(cycle0) * Dt + float(t) * dt
      if (t == n_sub_times-1):
        # At the end time of this cycle.  Take off a small amount from the time label (as it is repeated for the start of the next cycle)
        time[cycle0*n_sub_times + t] -= dt/2.0

      # Extract the av spectral density of the error fields
      err_u[cycle0*n_sub_times + t,:]  = process_field(err_raw_u)[:]
      err_v[cycle0*n_sub_times + t,:]  = process_field(err_raw_v)[:]
      err_w[cycle0*n_sub_times + t,:]  = process_field(err_raw_w)[:]
      err_rp[cycle0*n_sub_times + t,:] = process_field(err_raw_rp)[:]
      err_bp[cycle0*n_sub_times + t,:] = process_field(err_raw_bp)[:]




  # Plot the u data
  plot_field (wavenumber, scales, time, err_u, 'u', data_type, plot_dir, filesuffix, html_file)
  # Output the data to file
  output_field (wavenumber, scales, time, err_u, 'u', data_type, plot_dir)

  # Plot the v data
  plot_field (wavenumber, scales, time, err_v, 'v', data_type, plot_dir, filesuffix, html_file)
  # Output the data to file
  output_field (wavenumber, scales, time, err_v, 'v', data_type, plot_dir)

  # Plot the w data
  plot_field (wavenumber, scales, time, err_w, 'w', data_type, plot_dir, filesuffix, html_file)
  # Output the data to file
  output_field (wavenumber, scales, time, err_w, 'w', data_type, plot_dir)

  # Plot the rp data
  plot_field (wavenumber, scales, time, err_rp, 'rp', data_type, plot_dir, filesuffix, html_file)
  # Output the data to file
  output_field (wavenumber, scales, time, err_rp, 'rp', data_type, plot_dir)

  # Plot the bp data
  plot_field (wavenumber, scales, time, err_bp, 'bp', data_type, plot_dir, filesuffix, html_file)
  # Output the data to file
  output_field (wavenumber, scales, time, err_bp, 'bp', data_type, plot_dir)



if (output_type == 'web'):
  html_file.write ('</html>')
  html_file.close
