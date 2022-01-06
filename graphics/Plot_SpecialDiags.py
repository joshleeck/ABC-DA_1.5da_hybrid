#!/usr/bin/env python
# -------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import colors, cm
import matplotlib
import os

# Set the location of the data
#DataDir = '/media/ross/1297-5336/MeteororologyWork/DataAssim/RuthsModel/StdRes/SpecialDiags/AtmosParams/BalCors'
#DataDir = '/home/ross/tmp/StdRes'
#DataDir = '/scratch2/sws98rnb/StdRes/SpecialDiags/AtmosParams/Exp-GB-HB-AB'
DataDir = '/media/ross/1297-5336/MeteororologyWork/DataAssim/RuthsModel/StdRes/SpecialDiags/AtmosParams/Exp+GB+HB-AB'

# Std or Hi resolution?
res = 'Std'
#res = 'Hi'


# Set the location of the plots
PlotDir = DataDir

# Create the plotting directory
os.system('mkdir -p ' + PlotDir + '/Plots')


# As there could be a large number of plots, there is an option to output results on a web page
# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'eps'

if (output_type == 'web'):
  filesuffix = '.png'
else:
  filesuffix = '.eps'


# Which set of special diagnostics to plot?
DiagSwitch = []
# Real space correlations between control parameters
DiagSwitch.append(False)

# Spectral space correlations between control parameters
DiagSwitch.append(True)

# Actual control variables
DiagSwitch.append(False)

# Compute average KE spectrum
DiagSwitch.append(False)

# Strength of balance correlations with scale
DiagSwitch.append(False)

# Variances of scaled density (total, balanced, unbalanced)
DiagSwitch.append(False)








if (output_type == 'web'):
  # Set-up the html file
  html_file = open (PlotDir + '/Plots.html', 'w')
  html_file.write ('<html>\n')
  html_file.write ('<h1>Special Diagnostics</h1>\n')
  html_file.write (DataDir + '\n')
  html_file.write ('<br>')




if (DiagSwitch[0]):

  # --------------------------------------------------
  # Real space correlations between control parameters
  # --------------------------------------------------

  print ('Plotting real space correlations between control parameters')
  if (output_type == 'web'):
    html_file.write ('<h1>Real space correlations between control parameters</h1>\n')

  # Set the input filenames and the variable names
  InputFiles = []
  InputFiles.append('RealSpWith_psi.nc')
  InputFiles.append('RealSpWith_chi.nc')
  InputFiles.append('RealSpWith_ur.nc')
  InputFiles.append('RealSpWith_ub.nc')
  InputFiles.append('RealSpWith_uw.nc')
  FieldNames = []
  FieldNames.append('psi')
  FieldNames.append('chi')
  FieldNames.append('unbal_r')
  FieldNames.append('unbal_b')
  FieldNames.append('unbal_w')
  ncVarNames = []
  ncVarNames.append('v_1')
  ncVarNames.append('v_2')
  ncVarNames.append('v_3')
  ncVarNames.append('v_4')
  ncVarNames.append('v_5')


  if (output_type == 'web'):
    html_file.write ('<table>\n')


  # Loop around each variable
  for varA in range(5):

    if (output_type == 'web'):
      html_file.write ('<tr>\n')

    # Open the correlation file for this variable
    print ('Opening ', InputFiles[varA])
    nc_file = Dataset(DataDir + '/' + InputFiles[varA])
    longs   = nc_file.variables['longs'][:] / 1000.0
    level   = nc_file.variables['level'][:] / 1000.0

    nlongs  = len(longs)
    nlevs   = len(level)

    # Loop around each variable
    for varB in range(5):

      # Read-in the correlation field for this variable
      corfield = nc_file.variables[ncVarNames[varB]][:]

      # Compute its root mean square
      rms      = np.sqrt(np.mean(corfield[0:nlevs-1,0:nlongs] * corfield[0:nlevs-1,0:nlongs]))

      # Set a corner to be -1 and +1 (to set range)
      corfield[0][0] = -1.0
      corfield[1][0] = 1.0

      if (output_type == 'web'):
        html_file.write ('<td>' + FieldNames[varA] + ' with ' + FieldNames[varB] + '<br>\n')
        html_file.write ('RMS = ' + str(rms) + '<br>')


      matplotlib.rc('xtick', labelsize=16)
      matplotlib.rc('ytick', labelsize=16)
      fig, ax = plt.subplots()
      # cmap    = cm.get_cmap('Greys', 11)
      cmap    = cm.get_cmap('seismic', 11)
      cax     = ax.contourf(longs, level, corfield, cmap=cmap, vmin=-1.0, vmax=1.0)
      cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
      cax     = ax.contour (longs, level, corfield, colors='k', vmin=-1.0, vmax=1.0)
      # Labels etc
      ax.set_title('Corr ' + FieldNames[varA] + ' with ' + FieldNames[varB], fontsize=16)
      ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
      ax.set_ylabel('Vertical distance (km)', fontsize=16)
      #plt.show()
      graphics_file_name = 'Plots/' + 'RealSp_' + FieldNames[varA] + '_' + FieldNames[varB] + filesuffix
      plt.savefig(PlotDir + '/' + graphics_file_name, bbox_inches='tight')
      plt.close('all')

      if (output_type == 'web'):
        html_file.write ('<img src=' + graphics_file_name + ' width=300></td>\n')


    nc_file.close
    if (output_type == 'web'):
      html_file.write ('</tr>\n')

  if (output_type == 'web'):
    html_file.write ('</table>\n')





if (DiagSwitch[1]):

  # ------------------------------------------------------
  # Spectral space correlations between control parameters
  # ------------------------------------------------------

  if (res == 'Std'):
    # Average between 3 and 9 km (as in Bannister et al., 2018, Fig 9)
    levA = 12   # For standard resolution
    levB = 35   # For standard resolution
  else:
    levA = 30   # For high resolution
    levB = 88   # For high resolution


  print ('Plotting spectral space correlations between control parameters')
  if (output_type == 'web'):
    html_file.write ('<h1>Spectral space correlations between control parameters</h1>\n')

  # Set the input filenames and the variable names
  InputFiles = []
  InputFiles.append('psi')
  InputFiles.append('chi')
  InputFiles.append('ru')
  InputFiles.append('bu')
  InputFiles.append('w')
  FieldNames = []
  FieldNames.append('psi')
  FieldNames.append('chi')
  FieldNames.append('unbal_r')
  FieldNames.append('unbal_b')
  FieldNames.append('unbal_w')

  if (output_type == 'web'):
    html_file.write ('<table>\n')


  # Loop around each variable
  for varA in range(5):

    if (output_type == 'web'):
      html_file.write ('<tr>\n')

    # Loop around each variable
    for varB in range(5):

      # Open the correlation files for this variable
      variablename = 'SpecSpCor_' + InputFiles[varA] + '_' + InputFiles[varB]
      print ('Opening ', variablename, '.nc')
      nc_file = Dataset(DataDir + '/' + variablename + '.nc')
      if (varA == 0) and (varB == 0):
        wns     = nc_file.variables['x'][:] - 1.0
        level   = nc_file.variables['z'][:]
        nwns    = len(wns)
        nlevs   = len(level)
        scale   = np.zeros((nwns))
        for wn in range(1,nwns,1):
          scale[wn] = 540.0 / float(wn)
        scale[0] = scale[1]

      # Read-in the correlation field for this variable
      corfield = nc_file.variables[variablename][:]

      # Close the file
      nc_file.close

      # Compute its root mean square
      rms      = np.sqrt(np.mean(corfield[0:nlevs-1,0:nwns] * corfield[0:nlevs-1,0:nwns]))

      # Set a corner to be -1 and +1 (to set range)
      corfield[0][0] = -1.0
      corfield[1][0] = 1.0

      if (output_type == 'web'):
        html_file.write ('<td>' + FieldNames[varA] + ' with ' + FieldNames[varB] + '<br>\n')
        html_file.write ('RMS = ' + str(rms) + '<br>')


      #### Plot correlations wavenumber vs level
      matplotlib.rc('xtick', labelsize=16)
      matplotlib.rc('ytick', labelsize=16)
      fig, ax = plt.subplots()
      # cmap    = cm.get_cmap('Greys', 11)
      cmap    = cm.get_cmap('seismic', 11)
      cax     = ax.contourf(wns, level, corfield, cmap=cmap, vmin=-1.0, vmax=1.0)
      cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
      cax     = ax.contour (wns, level, corfield, colors='k', vmin=-1.0, vmax=1.0)
      # Labels etc
      ax.set_title('Corr ' + FieldNames[varA] + ' with ' + FieldNames[varB], fontsize=16)
      ax.set_xlabel('Wavenumber', fontsize=16)
      ax.set_ylabel('Vertical distance (km)', fontsize=16)
      #plt.show()
      graphics_file_name = 'Plots/' + 'SpecSp_' + FieldNames[varA] + '_' + FieldNames[varB] + filesuffix
      plt.savefig(PlotDir + '/' + graphics_file_name, bbox_inches='tight')
      plt.close('all')

      if (output_type == 'web'):
        html_file.write ('<img src=' + graphics_file_name + ' width=300>\n')


      #### Average the correlations vertically
      corfield_av = np.zeros((nwns))
      for wn in range(nwns):
        corfield_av[wn] = np.mean(corfield[levA:levB,wn])

      fig, ax = plt.subplots()
      ax.set_title('Corr ' + FieldNames[varA] + ' with ' + FieldNames[varB] + ' vert av', fontsize=16)
      ax.set_xlabel('scale (km)')
      ax.set_ylabel('correlation')
      ax.set_xscale('log')
      # Set the y plotting limits
      plt.ylim(-1.0, 1.0)
      ax.plot(scale, corfield_av, '.', linewidth=2, color='k')

      graphics_file_name = 'Plots/' + 'SpecSp_' + FieldNames[varA] + '_' + FieldNames[varB] + 'VertAv' + filesuffix
      plt.savefig(PlotDir + '/' + graphics_file_name, bbox_inches='tight')
      plt.close('all')

      if (output_type == 'web'):
        html_file.write ('<br><img src=' + graphics_file_name + ' width=300></td>\n')

    if (output_type == 'web'):
      html_file.write ('</tr>\n')

  if (output_type == 'web'):
    html_file.write ('</table>\n')







if (DiagSwitch[2]):

  # -----------------------------------------------------
  # Spectral space correlations between control variables
  # -----------------------------------------------------

  # Average over vertical modes
  modeA = 1
  modeB = 3


  print ('Plotting correlations between control variables')
  if (output_type == 'web'):
    html_file.write ('<h1>Correlations between control variables</h1>\n')

  # Set the input filenames and the variable names
  InputFiles = []
  InputFiles.append('psi')
  InputFiles.append('chi')
  InputFiles.append('ru')
  InputFiles.append('bu')
  InputFiles.append('w')
  FieldNames = []
  FieldNames.append('psi')
  FieldNames.append('chi')
  FieldNames.append('unbal_r')
  FieldNames.append('unbal_b')
  FieldNames.append('unbal_w')


  if (output_type == 'web'):
    html_file.write ('<table>\n')


  # Loop around each variable
  for varA in range(5):

    if (output_type == 'web'):
      html_file.write ('<tr>\n')

    # Loop around each variable
    for varB in range(5):

      # Open the correlation files for this variable
      variablename = 'ConVar_' + InputFiles[varA] + '_' + InputFiles[varB]
      print ('Opening ', variablename, '.nc')
      nc_file = Dataset(DataDir + '/' + variablename + '.nc')
      if (varA == 0) and (varB == 0):
        wns     = nc_file.variables['x'][:] - 1.0
        level   = nc_file.variables['z'][:]
        nwns    = len(wns)
        nlevs   = len(level)
        nmodes  = nlevs
        scale   = np.zeros((nwns))
        for wn in range(1,nwns,1):
          scale[wn] = 540.0 / float(wn)
        scale[0] = scale[1]
        modes = range(1,nmodes+1)

      # Read-in the correlation field for this variable
      corfield = nc_file.variables[variablename][:]
      
      # Close the file
      nc_file.close

      # Compute its root mean square
      rms = np.sqrt(np.mean(corfield[1:nmodes,0:nwns] * corfield[1:nmodes,0:nwns]))

      # Set a corner to be -1 and +1 (to set range)
      corfield[0][0] = -1.0
      corfield[1][0] = 1.0

      if (output_type == 'web'):
        html_file.write ('<td>' + FieldNames[varA] + ' with ' + FieldNames[varB] + '<br>\n')
        html_file.write ('RMS = ' + str(rms) + '<br>')


      matplotlib.rc('xtick', labelsize=16)
      matplotlib.rc('ytick', labelsize=16)
      fig, ax = plt.subplots()
      # cmap    = cm.get_cmap('Greys', 11)
      cmap    = cm.get_cmap('seismic', 11)
      cax     = ax.contourf(wns, modes, corfield, cmap=cmap, vmin=-1.0, vmax=1.0)
      cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
      cax     = ax.contour (wns, modes, corfield, colors='k', vmin=-1.0, vmax=1.0)
      # Labels etc
      ax.set_title('Corr ' + FieldNames[varA] + ' with ' + FieldNames[varB], fontsize=16)
      ax.set_xlabel('Wavenumber', fontsize=16)
      ax.set_ylabel('Vertical mode', fontsize=16)
      #plt.show()
      graphics_file_name = 'Plots/' + 'ConVar_' + FieldNames[varA] + '_' + FieldNames[varB] + filesuffix
      plt.savefig(PlotDir + '/' + graphics_file_name, bbox_inches='tight')
      plt.close('all')

      if (output_type == 'web'):
        html_file.write ('<img src=' + graphics_file_name + ' width=300>\n')


      #### Average the correlations over a number of vertical modes
      corfield_av = np.zeros((nwns))
      for wn in range(nwns):
        corfield_av[wn] = np.mean(corfield[modeA:modeB,wn])

      fig, ax = plt.subplots()
      ax.set_title('Corr ' + FieldNames[varA] + ' with ' + FieldNames[varB] + ' mode av', fontsize=16)
      ax.set_xlabel('scale (km)')
      ax.set_ylabel('correlation')
      ax.set_xscale('log')
      # Set the y plotting limits
      plt.ylim(-1.0, 1.0)
      ax.plot(scale, corfield_av, '.', linewidth=2, color='k')

      graphics_file_name = 'Plots/' + 'ConVar_' + FieldNames[varA] + '_' + FieldNames[varB] + 'ModeAv' + filesuffix
      plt.savefig(PlotDir + '/' + graphics_file_name, bbox_inches='tight')
      plt.close('all')

      if (output_type == 'web'):
        html_file.write ('<br><img src=' + graphics_file_name + ' width=300></td>\n')

    if (output_type == 'web'):
      html_file.write ('</tr>\n')

  if (output_type == 'web'):
    html_file.write ('</table>\n')







if (DiagSwitch[3]):

  # -----------
  # KE spectrum
  # -----------

  if (res == 'Std'):
    # Average between 3 and 9 km (as in Bannister et al., 2018, Fig 9)
    levA = 12   # For standard resolution
    levB = 35   # For standard resolution
  else:
    levA = 30   # For high resolution
    levB = 88   # For high resolution

  print ('Plotting KE spectrum')
  if (output_type == 'web'):
    html_file.write ('<h1>KE spectrum</h1>\n')

  print ('Reading from data file')
  input_file = open (DataDir + '/KEspec.dat')
  line  = input_file.readline()
  line  = input_file.readline()
  nwns  = int(line)
  line  = input_file.readline()
  line  = input_file.readline()
  nlevs = int(line)
  line  = input_file.readline()
  line  = input_file.readline()

  KEspecs = np.zeros((nlevs,nwns))

  for lev in range(nlevs):
    line  = input_file.readline()
    KEspecs[lev][0:nwns] = line.split()[0:nwns]
  input_file.close()

  # Remove 1 from the number of wavenumbers (to avoid bad data)
  nwns = nwns - 1

  # Average between vertical levels
  KEspec_vertav = np.zeros((nwns))
  for wn in range(nwns):
    KEspec_vertav[wn] = np.mean(KEspecs[levA:levB, wn])

  # Generate the wavenumbers
  wns = np.linspace(0.0, nwns, nwns)

  # Generate the lengthscales (km)
  lscales = np.zeros((nwns))
  for wn in range(1,nwns):
    lscales[wn] = 1.5 * 360.0 / float(wn)

  fig, ax = plt.subplots()
  ax.set_yscale('log')
  ax.set_xscale('log')
  ax.set_xlabel('Wavenumber')
  ax.set_ylabel('Spectral KE')
  ax.plot(wns[1:], KEspec_vertav[1:], '.', linewidth=2, color='k')
  x1,x2 = ax.get_xlim()

  C = 5.0
  ax.plot(wns[1:], C * wns[1:]**(-3.0), color='black')
  ax.plot(wns[1:], C * wns[1:]**(-5.0/3.0), color='black', linestyle='dashed')

  ax.annotate('k^{-3}', xy=(200.0,1e-6), xytext=(200.0,1e-6), size='large', color='black')
  ax.annotate('k^{-5/3}', xy=(200.0,1e-3), xytext=(200.0,1e-3), size='large', color='black')

  ax2=plt.twiny()
  ax2.set_xlim(lscales[1], lscales[-1])
  ax2.set_xscale('log') # or use ax2.set_yticks
  ax2.set_xlabel('Wavelength (km)')


  graphics_file_name = 'Plots/' + 'KEspec' + filesuffix
  plt.savefig(PlotDir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')

  if (output_type == 'web'):
    html_file.write ('<img src=' + graphics_file_name + ' width=300></td>\n')

  # Output the spectrum as a text file
  output_file = open (PlotDir + '/' + 'KEspecPythonAv.dat', 'w')
  for wn in range(nwns):
    output_file.write (str(wn) + '  ' + str(lscales[wn]) + '  ' + str(KEspec_vertav[wn]) + '\n')
  output_file.close()





if (DiagSwitch[4]):

  # -------------------------------------------
  # Strength of balance correlations with scale
  # -------------------------------------------

  if (res == 'Std'):
    # Average between 3 and 9 km (as in Bannister et al., 2018, Fig 9)
    levA  = 12   # For standard resolution
    levB  = 35   # For standard resolution
    nlevs = 60  # For standard resolution
  else:
    levA  = 30   # For high resolution
    levB  = 88   # For high resolution
    nlevs = 150 # For high resolution

  print ('Plotting strength of balance correlations with scale')
  if (output_type == 'web'):
    html_file.write ('<h1>Strength of balance correlations with scale</h1>\n')

  print ('Reading from data file (GeoBalCorr.dat)')
  input_file = open (DataDir + '/GeoBalCorr.dat')
  line  = input_file.readline()
  line  = input_file.readline()
  SmoothingGridBoxes = line.split()
  nscales = len(SmoothingGridBoxes)
  print ('There are ', nscales, 'different smoothing scales')
  GeoBalCorr = np.zeros((nscales, nlevs))
  line  = input_file.readline()
  for z in range(nlevs):
    line  = input_file.readline()
    GeoBalCorr[0:nscales,z] = line.split()[1:nscales+1]
  input_file.close()

  print ('Reading from data file (HydBalCorr.dat)')
  input_file = open (DataDir + '/HydBalCorr.dat')
  line  = input_file.readline()
  line  = input_file.readline()
  SmoothingGridBoxes = line.split()
  nscales = len(SmoothingGridBoxes)
  print ('There are ', nscales, 'different smoothing scales')
  for item in range(nscales):
    SmoothingGridBoxes[item] = float(SmoothingGridBoxes[item])

  HydBalCorr = np.zeros((nscales, nlevs))
  line  = input_file.readline()
  for z in range(nlevs):
    line  = input_file.readline()
    HydBalCorr[0:nscales,z] = line.split()[1:nscales+1]
  input_file.close()


  # Average over required levels
  GeoBal = np.zeros((nscales))
  for scale in range(nscales):
    GeoBal[scale] = np.abs(np.mean(GeoBalCorr[scale,levA:levB]))

  HydBal = np.zeros((nscales))
  for scale in range(nscales):
    HydBal[scale] = np.abs(np.mean(HydBalCorr[scale,levA:levB]))

  # Plot the results
  fig, ax1 = plt.subplots()
  ax2      = ax1.twinx()

  print ('Plotting geostrophic correlations')
  ax1.set_ylim([0.0,0.4])
  ax1.plot(SmoothingGridBoxes, GeoBal, color='black', linewidth='2')
  # ax1.legend(loc='upper left')
  # ax1.legend(loc=(0.3, legypos[dirno]))
  ax1.set_xlabel('Scale (grid boxes)', color='black', fontsize=16)
  ax1.set_ylabel('|geostrophic correlation| (black)', color='black', fontsize=16)
  #for tl in ax1.get_yticklabels():
  #  tl.set_color('blue')

  print ('Plotting hydrostatic correlations')
  ###ax2.set_ylim([0.85,1.0])     # For standard resolution
  ax2.set_ylim([0.997,1.0])    # For high resolution
  ax2.plot(SmoothingGridBoxes, HydBal, color='grey', linewidth='2')
  #ax2.legend(loc='lower right')
  #ax2.legend(loc=(0.6, legypos[dirno]))
  ax2.set_ylabel('|hydrostatic correlation| (grey)', color='black', fontsize=16)
  #for tl in ax2.get_yticklabels():
  #  tl.set_color('red')

  plt.title('Geostrophic and hydrostatic correlations', color='black', fontsize=16)

  graphics_file_name = 'Plots/' + 'BalCor' + filesuffix
  plt.savefig(PlotDir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')

  if (output_type == 'web'):
    html_file.write ('<img src=' + graphics_file_name + ' width=300></td>\n')








if (DiagSwitch[5]):

  # ------------------------------------------------------
  # Variances of scaled density (total, balanced, unbalanced)
  # in spectral space
  # ------------------------------------------------------

  if (res == 'Std'):
    # Average between 3 and 9 km (as in Bannister et al., 2018, Fig 9)
    #levA = 12   # For standard resolution
    #levB = 35   # For standard resolution  
    # Around 4 km
    levA = 16   # For standard resolution
    levB = 17   # For standard resolution
  else:
    # Average between 3 and 9 km (as in Bannister et al., 2018, Fig 9)
    #levA = 30   # For high resolution
    #levB = 88   # For high resolution

    # Around 4 km
    levA = 42   # For high resolution
    levB = 44   # For high resolution


  print ('Plotting spectral space variances of r, rb, ru')
  if (output_type == 'web'):
    html_file.write ('<h1>Spectral space variances of r, rb, ru</h1>\n')

  # Open the files and read-in the variances
  variablename = 'TotrVarSpec'
  filename     = DataDir + '/' + variablename + '.nc'
  print ('Opening ', filename)
  nc_file = Dataset(filename)
  # Read-in axis information
  wns     = nc_file.variables['x'][:] - 1.0
  level   = nc_file.variables['z'][:]
  nwns    = len(wns)
  nlevs   = len(level)
  scale   = np.zeros((nwns))
  for wn in range(1,nwns,1):
    scale[wn] = 540.0 / float(wn)
  scale[0] = scale[1] * 2.0
  # Read-in the variance field
  totrvar = nc_file.variables[variablename][:]
  nc_file.close

  variablename = 'BalrVarSpec'
  filename     = DataDir + '/' + variablename + '.nc'
  print ('Opening ', filename)
  nc_file = Dataset(filename)
  # Read-in the variance field
  balrvar = nc_file.variables[variablename][:]
  nc_file.close

  variablename = 'UnbrVarSpec'
  filename     = DataDir + '/' + variablename + '.nc'
  print ('Opening ', filename)
  nc_file = Dataset(filename)
  # Read-in the variance field
  unbrvar = nc_file.variables[variablename][:]
  nc_file.close


  #### Average the variances vertically
  totrvar_av = np.zeros((nwns))
  balrvar_av = np.zeros((nwns))
  unbrvar_av = np.zeros((nwns))
  for wn in range(nwns):
    totrvar_av[wn] = np.mean(totrvar[levA:levB,wn])
    balrvar_av[wn] = np.mean(balrvar[levA:levB,wn])
    unbrvar_av[wn] = np.mean(unbrvar[levA:levB,wn])

  fig, ax = plt.subplots()
  ax.set_title('Vertically averaged scaled den spec', fontsize=16)
  ax.set_xlabel('scale (km)', fontsize=16)
  ax.set_ylabel('correlation', fontsize=16)
  ax.set_xscale('log')
  ax.set_yscale('log')
  # Set the y plotting limits
  ax.plot(scale, totrvar_av, ls='solid', linewidth=2, color='k', label='Total')
  ax.plot(scale, balrvar_av, ls='dotted', linewidth=2, color='k', label='Bal')
  ax.plot(scale, totrvar_av, ls='dashdot', linewidth=2, color='k', label='Unbal')
  # Show the legend
  ax.legend(fontsize=16)

  graphics_file_name = 'Plots/ScaledDenVar_VertAv' + filesuffix
  plt.savefig(PlotDir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')

  if (output_type == 'web'):
    html_file.write ('<img src=' + graphics_file_name + ' width=300>\n')


if (output_type == 'web'):
  html_file.close
