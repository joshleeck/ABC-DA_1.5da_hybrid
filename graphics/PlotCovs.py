#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to:
#   -read and plot the raw or implied forecast error covariance fields.
#
# Please edit the input details (e.g. location of data files)
# at the start of the main part of the code.
# This is located after any function definitions below.
#
# Ross Bannister, Dec 2016, Mar 2018
# -------------------------------------------------------------------

# ===================================================================
def max_field (field):
  # Subroutine to return maxmimum absolute value of field
  maxlev  = []
  for lev in field:
    maxlev.append(max(lev))
  globalmax = max(maxlev)
  return globalmax

def min_field (field):
  # Subroutine to return maxmimum absolute value of field
  minlev  = []
  for lev in field:
    minlev.append(min(lev))
  globalmin = min(minlev)
  return globalmin

def remove_extremes (field, low, high):
  # Set values less than low to low
  # Set values larger than high to high
  new_field = []
  for line in field:
    new_line = []
    for val in line:
      if val < low: 
        new_line.append(low)
      else:
        if val > high:
          new_line.append(high)
        else:
          new_line.append(val)
    new_field.append(new_line)
  return new_field
 
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import colors, cm
import matplotlib
import os
 
# ===================================================================

class MidpointNormalize(matplotlib.colors.Normalize):
  def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
    self.midpoint=midpoint
    matplotlib.colors.Normalize.__init__(self,vmin,vmax,clip)
  def __call__(self, value, clip=None):
    x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
    return np.ma.masked_array(np.interp(value, x, y))

# ===================================================================
  
# ===================================================================

# What to plot
plot_raw_covs     = False
plot_implied_covs = True

# Set directory of the input ensemble
data_dir_RawCovs     = 'test/examples/Master_RawCov/'
data_dir_ImpliedCovs = 'test/examples/Master_ImpliedCov/'

# Set the domain dimensions
nlongs   = 364
nlevs    = 60
tstep    = 0

# How many ensembles? (Used only if plot_ensemble_covs is true)
Nens       = 1
NEnsMems   = 30  # This is the number of ensembles * number of latitudes used in Fortran code

# What locations of source points? (Lists must match namelist for implied/raw covs run)
longindex = [182]
levindex  = [30]
points    = len(longindex)
if (len(levindex) != points):
  print ('Error - there must be the same number of longindex and levindex entries')

# Output type
# As there could be a large number of plots, there is an option to output results on a web page
# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'web'


if (output_type == 'web'):
  filesuffix = '.png'
else:
  filesuffix = '.eps'


varfilenames = ['u', 'v', 'w', 'r', 'b', 'tracer']
varnames     = ['u', 'v', 'w', 'r_prime', 'b_prime', 'tracer']





# ===================================================================
for type_of_cov in range(2):
  if (type_of_cov == 0):
    go       = plot_raw_covs
    data_dir = data_dir_RawCovs
    htmlfile = 'Plots_RawCovs.html'
    print ('Plotting raw covariances ...')
  else:
    go       = plot_implied_covs
    data_dir = data_dir_ImpliedCovs
    htmlfile = 'Plots_ImpliedCovs.html'
    print ('Plotting implied covariances ...')
  plotdir = data_dir


  if (go):

    os.system('mkdir -p ' + plotdir + '/Plots')

    if (output_type == 'web'):
      html_file = open (plotdir + '/' + htmlfile, 'w')
      html_file.write ('<html>\n')
      if (type_of_cov == 0):
        html_file.write ('<h1>Raw covariances for ABC model</h1>\n')
      else:
        html_file.write ('<h1>Implied covariances for ABC model</h1>\n')

    for point in range(points):
      print ('Plotting data for spatial point ', point+1)
      val_lim = [[9, 16, 0.6, 0.01, 0.12, 0.5],[0.20, 0.20, 0.0025, 0.00006, 0.0025, 0.5],[0.5, 0.8, 0.04, 0.0004, 0.006, 0.5],[0.005, 0.00025, 0.0003, 0.000005, 0.00004, 0.5],[0.05, 0.08, 0.004, 0.00004, 0.0006, 0.5],[1,1,1,1,1,1]]
      for var_source in range(6):
        print ('Dealing with source model variable ', var_source+1, ': ', varnames[var_source])
        if (output_type == 'web'):
          html_file.write ('<h2>Source point ' + str(point+1) + ', source variable ' + varnames[var_source] + '</h2>\n')
          html_file.write ('<table cols=6 border=0>\n')
          html_file.write ('<tr>\n')

        filename = 'Point_' + str(point+1).zfill(3) + '_delta' + varfilenames[var_source] + '.nc'
        # Open the file for this source point/variable
        nc_file = Dataset(data_dir + '/' + filename)

        # Get axis data
        if ((point == 0) and (var_source==0)):
          longs_u    = nc_file.variables['longs_u'][:] / 1000.0
          longs_v    = nc_file.variables['longs_v'][:] / 1000.0
          full_level = nc_file.variables['full_level'][:] / 1000.0
          half_level = nc_file.variables['half_level'][:] / 1000.0

        maxvalue_list = val_lim[var_source]
        # Loop round each field variable
        for var_field in range(6):
          print ('Dealing with field model variable ', var_field+1, ': ', varnames[var_field])
          field        = nc_file.variables[varnames[var_field]][tstep,:,:]
          minvalue     = np.abs(np.min(field))
          maxvalue     = np.abs(np.max(field))
          #maxvalue_use = max(minvalue,maxvalue)
          #minvalue_use = -maxvalue_use
          levels = np.linspace(-val_lim[var_source][var_field],val_lim[var_source][var_field],11)
          mid_val = levels[4]
          levels = np.delete(levels, 5) 
          levels = np.insert(levels, 5, -mid_val/2.)
          levels = np.insert(levels, 5, -mid_val/8.)
          levels = np.insert(levels, 5, -mid_val/20.)
          levels = np.insert(levels, 5, -mid_val/100.)
          levels = np.insert(levels, 5, mid_val/2.)
          levels = np.insert(levels, 6, mid_val/8.)
          levels = np.insert(levels, 7, mid_val/20.)
          levels = np.insert(levels, 8, mid_val/100.)
          #print levels
          mn           = np.mean(field)
          err          = field - mn
          rms          = np.sqrt(np.mean(err * err))
          if (rms == 0.0):
            minvalue_use = -1.0
            maxvalue_use = 1.0

          # Plot the implied covariance field
          matplotlib.rc('xtick', labelsize=16)
          matplotlib.rc('ytick', labelsize=16)
          fig, ax = plt.subplots()
          #cmap    = cm.get_cmap('Greys', 11)
          cmap    = cm.get_cmap('seismic')

          if (var_field+1 == 1):
            # Zonal wind
            cax = ax.contourf(longs_u, half_level, field, norm=MidpointNormalize(midpoint=0.), cmap=cmap, levels=levels)
            if (rms > 0.0):
              cax1 = ax.contour (longs_u, half_level, field, colors='k', levels=levels)
            plt.scatter(longs_u[longindex[point]], half_level[levindex[point]], marker='x', c='Y', s=400.0, linewidth=4)
          elif (var_field+1 == 2):
            # Meridional wind
            cax = ax.contourf(longs_v, half_level, field, norm=MidpointNormalize(midpoint=0.), cmap=cmap, levels=levels)
            if (rms > 0.0):
              cax1 = ax.contour (longs_v, half_level, field, colors='k', levels=levels)
            plt.scatter(longs_v[longindex[point]], half_level[levindex[point]], marker='x', c='Y', s=400.0, linewidth=4)
          elif (var_field+1 == 3):
            # Vertical wind
            cax = ax.contourf(longs_v, full_level, field, norm=MidpointNormalize(midpoint=0.), cmap=cmap, levels=levels)
            if (rms > 0.0):
              cax1 = ax.contour (longs_v, full_level, field, colors='k', levels=levels)
            plt.scatter(longs_v[longindex[point]], full_level[levindex[point]], marker='x', c='Y', s=400.0, linewidth=4)
          elif (var_field+1 == 4):
            # r_prime
            cax = ax.contourf(longs_v, half_level, field, norm=MidpointNormalize(midpoint=0.), cmap=cmap, levels=levels)
            if (rms > 0.0):
              cax1 = ax.contour (longs_v, half_level, field, colors='k', levels=levels)
            plt.scatter(longs_v[longindex[point]], half_level[levindex[point]], marker='x', c='Y', s=400.0, linewidth=4)
          elif (var_field+1 == 5):
            # b_prime
            cax = ax.contourf(longs_v, full_level, field, norm=MidpointNormalize(midpoint=0.), cmap=cmap, levels=levels)
            if (rms > 0.0):
              cax1 = ax.contour (longs_v, full_level, field, colors='k', levels=levels)
            plt.scatter(longs_v[longindex[point]], full_level[levindex[point]], marker='x', c='Y', s=400.0, linewidth=4)
          elif (var_field+1 == 6):
            # tracer
            cax = ax.contourf(longs_v, full_level, field, cmap=cmap, levels=levels)
            if (rms > 0.0):
              cax1 = ax.contour (longs_v, full_level, field, colors='k', levels=levels)
            plt.scatter(longs_v[longindex[point]], full_level[levindex[point]], marker='x', c='Y', s=400.0, linewidth=4)

          # Labels etc
          cbar = fig.colorbar(cax, orientation='vertical', cmap=cmap, ticks=levels)
          ax.set_title(varnames[var_source] + ' - ' + varnames[var_field], fontsize=16)
          ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
          ax.set_ylabel('Vertical distance (km)', fontsize=16)

          #plt.show()
          figfilename = 'Plots/Point' + str(point+1).zfill(3) + '__' + varnames[var_source] + '-' + varnames[var_field] + filesuffix
          plt.savefig(plotdir + '/' + figfilename, bbox_inches='tight')
          plt.close('all')

          if (output_type == 'web'):
            html_file.write ('<td><img src=' + figfilename + ' width=300></td>\n')

        if (output_type == 'web'):
          html_file.write ('</tr>\n')
          html_file.write ('</table>\n')


        nc_file.close


    if (output_type == 'web'):
      html_file.write ('</table>\n')
      html_file.write ('</html>')
      html_file.close()
      print ('An html file has been created to view the figures')
      print ('Please view the following html file with your browser')
      print (plotdir + '/Plots.html')
      print ()
