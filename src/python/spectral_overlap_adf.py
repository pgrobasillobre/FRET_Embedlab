
from __future__ import division

import sys
import re
import numpy as np
import time
import matplotlib.pyplot as plt


# ==================================================
#                                                  #
#   EXECUTION details: python3 overlap_adf.py -h   #
#                                                  #
# ==================================================

# --> Initialize variables
min_energy    = 0.0   # Minimum Energy for the convoluted spectra in eV
max_energy    = 5.0   # Maximum Energy for the convoluted spectra in eV
grid_points   = 1000  # Number of energy points in the spectral range
fwhm          = 0.8   # Full Width at Half Maximum (FWHM) in eV 

# --> Conversion factors
ev_to_hartree = 1.0/27.211399


# ==============================
#           FUNCTIONS
# ==============================

def read_command_line(command_line):

   if len(sys.argv) < 2:
      print('')
      print('')
      print('   Please provide inputs files:')
      print('')
      print('      For more details --> python3 fret_coulomb.py -h')
      print('')
      print('')
      sys.exit()

   elif (len(sys.argv) > 6):
      print('')
      print('')
      print('   Too many arguments in input')
      print('')
      print('      For more details --> python3 fret_coulomb.py -h')
      print('')
      print('')
      sys.exit()
 
   elif sys.argv[1] == '-h' or sys.argv[1] == '-help':
      print('')
      print('')
      print(' Execution --> python3 overlap_adf.py adf1.log #state adf2.log #state plot_gaussians[optional]')
      print('')
      print('')
      sys.exit()

   else:
      plot_gaus = False
      if (len(sys.argv) == 6): plot_gaus = True
      return(sys.argv[1],int(sys.argv[2]),sys.argv[3],int(sys.argv[4]),plot_gaus)
 

# -----

def read_adf_tddft_log(infile,state):

   au_to_ev = 27.211324570273

   counter  = 0
   n_states = 0
   exc      = 0.0
   osc      = 0.0

   n_states_found    = False
   state_found       = False
   excitations_found = False

   with open(infile,'r') as f:

      lines    = f.readlines()

      for line in lines:

          if (re.search('lowest',line,re.IGNORECASE)): # Check the number of states used in the calculation
             n_states_found    = True
             n_states = int(line.split()[1])

          if ' no.  E/a.u.        E/eV      f           dE/a.u.' in line: excitations_found = True # Check that we have the excitations

          # Checkpoint
          if (excitations_found and not n_states_found):
             print(' ')
             print(' ')
             print('    Error: calculated TD-DFT states not found in file ' + str(infile))
             print(' ')
             print(' ')
             sys.exit()

          if (excitations_found): counter +=1 # To know when to read the calculated TD-DFT states 

          if (excitations_found and counter > 2 and counter-2==int(state) and counter-2 <= n_states): # Read the found TD-DFT state of interest

             state_found = True 

             if (int(line.split()[0]) == int(state)):
                exc = float(line.split()[1]) * au_to_ev 
                osc = float(line.split()[3]) * au_to_ev

             break # If we have found the state exit the for loop

          elif (counter-2 > n_states): # If the counter overpass the total number of states then there is a problem
                print(' ')
                print(' ')
                print('    Error: state ' + str(int(state)) + ' not found in file ' + str(infile))
                print(' ')
                print(' ')
                sys.exit()

   return (exc,osc)


# -----

def single_gaussian(grid_points, exc, osc, fwhm, min_energy, max_energy):

    # Shape of the gaussian: f(x) = (1 / (σ * √(2π))) * exp(-((x - μ)^2) / (2σ^2))
    # (weighted by oscillator strength)

    # Create an array of energy values in the specified range
    energies = np.linspace(min_energy, max_energy, grid_points)

    # Calculate the standard deviation (sigma) from the FWHM
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))

    # Calculate the Gaussian convolution
    gaussian = osc / (sigma * np.sqrt(2 * np.pi)) * np.exp(-((energies - exc) ** 2) / (2 * sigma**2))

    return energies, gaussian


          
# ==============================
#         MAIN PROGRAM 
# ==============================


# START TIMER
start = time.time()


# --> Read command line, excitation energies and oscillator strengths

adf_1, adf_1_state, adf_2, adf_2_state, plot_gaussians = read_command_line(sys.argv)

exc_1, osc_1 = read_adf_tddft_log(adf_1, adf_1_state)
exc_2, osc_2 = read_adf_tddft_log(adf_2, adf_2_state)


# --> Create Gaussian convolution

energies_1, gaussian_1 = single_gaussian(grid_points, exc_1, osc_1, fwhm, min_energy, max_energy)
energies_2, gaussian_2 = single_gaussian(grid_points, exc_2, osc_2, fwhm, min_energy, max_energy)


# --> Calculate the overlap using numerical integration (area under the product of the two Gaussians)

overlap = np.trapz(gaussian_1 * gaussian_2, energies_1) 


# END TIMER
end = time.time()


print('')
print('   Spectral Overlap (a.u.): ' + str(overlap * ev_to_hartree))
print('                            ')
print('   Omega_0 (a.u.): ' + str((abs(exc_1 - exc_2)/2.0) * ev_to_hartree))
print('')
print('   ----------------------------------')
print('   NORMAL TERMINATION')
print('')
print('   COMPUTATIONAL TIME: ', str(round(end-start,2)), ' s')
print('')
print('')


# --> Plot the convolution [OPTIONAL]

if (plot_gaussians): 

   # Figure size
   plt.figure(figsize=(9.0, 7.0))

   # Fonts and labels
   plt.rcParams['font.family'] = 'Times New Roman'
   plt.rcParams['mathtext.fontset'] = 'custom'
   plt.rcParams['mathtext.rm'] = 'Times New Roman'
   plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
   plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'
   plt.rcParams['mathtext.sf'] = 'Times New Roman'
   plt.rcParams['mathtext.default'] = 'regular'

   fontsize_axes   = 20
   fontsize_labels = 22
   fontsize_title   = 25

   # Ticks, labels and title
   plt.tick_params(axis='x', which='major', labelsize=fontsize_axes)
   plt.tick_params(axis='y', which='major', labelsize=fontsize_axes)

   plt.xlabel("Energy (eV)", fontsize=fontsize_labels, labelpad = 8.0)
   plt.ylabel("Osc. Strength (eV)", fontsize=fontsize_labels, labelpad=18.0)

   plt.title("Gaussian Convolution\nof TDDFT Excitations",fontsize=fontsize_title, pad=12.0)

   # Plot
   plt.plot(energies_1, gaussian_1, color='red',  label = 'State ' + str(adf_1_state))
   plt.plot(energies_2, gaussian_2, color='blue', label = 'State ' + str(adf_2_state))

   plt.legend(framealpha=1,fontsize=fontsize_labels)
   plt.grid(True)

   ### Calculate the maximum y-values for each Gaussian
   ##max_y1 = max(gaussian_1)
   ##max_y2 = max(gaussian_2)
   ##
   ### Plot vertical lines (pulses) from the bottom (0.0) to the maximum y-values
   ##plt.vlines(exc_1, 0.0, max_y1, colors='red',  linestyles='dashed', label='Pulse for Gaussian 1', alpha=0.2)
   ##plt.vlines(exc_2, 0.0, max_y2, colors='blue', linestyles='dashed', label='Pulse for Gaussian 2', alpha=0.2)

   plt.show()


# ==============================
#         MAIN PROGRAM 
# ==============================








