
from __future__ import division

import sys
import os
import re
import glob
import math
import numpy as np
import time
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D


# ========================================================
#                                                        #
#   EXECUTION details: python3 field_on_xyz_grid.py -h   #
#                                                        #
# ========================================================

# --> Initialize variables
incident_field_intenstiy = 0.1



# ==============================
#           FUNCTIONS
# ==============================

# ------------------------------- #
# ------ Read command line ------ #

def read_command_line(command_line):

   if sys.argv[1] == '-h' or sys.argv[1] == '-help':
      print('')
      print('')
      print(' Execution --> python3 field_on_xyz_grid.py tar_file grid(xyz format) dir{x/y/z}')
      print('')
      print('')
      sys.exit()

   elif len(sys.argv) < 4:
      print('')
      print('')
      print('   Please provide all inputs files:')
      print('')
      print('      For more details --> python3 field_on_xyz_grid.py -h')
      print('')
      print('')
      sys.exit()

   elif (len(sys.argv) > 4):
      print('')
      print('')
      print('   Too many arguments in input')
      print('')
      print('      For more details --> python3 field_on_xyz_grid.py -h')
      print('')
      print('')
      sys.exit()
 
   else:
      tar_file  = sys.argv[1]
      grid_file = sys.argv[2]
      direction = sys.argv[3]

      return(tar_file,grid_file,direction)


# ------------------------------- #
# --------- Check inputs -------- #

def check_inputs(tar_file,grid_file,direction):
#
   """
   Function to check tar and grid files

   :tar_file : tar.gz file from nanofq (bem_tommaso_ief_branch)
   :grid_file: grid file in xyz format
   :direction: field direction
   """
#
   if not os.path.exists(tar_file):
      print(' ')
      print('  STOP: tar.gz nanofq file "' + tar_file + '" is not in the current folder')
      print(' ')
      sys.exit()

   if not os.path.exists(grid_file):
      print(' ')
      print('  STOP: grid file "' + grid_file + '" is not in the current folder')
      print(' ')
      sys.exit()

   if tar_file[-6:] != 'tar.gz':
      print(' ')
      print('  STOP: tar.gz extension not found in "' + tar_file + '" file')
      print(' ')
      sys.exit()

   if grid_file[-4:] != '.xyz':
      print(' ')
      print('  STOP: .xyz extension not found in "' + grid_file + '" file')
      print(' ')
      sys.exit()

   if direction !='x' and direction !='y' and direction !='z':
      print(' ')
      print('  STOP: direction "' + direction + '" not recognised')
      print(' ')
      print('    Options:')
      print(' ')
      print('      x, y, z')
      print(' ')
      sys.exit()


# ------------------------------- #
# ------- Read geometries ------- #

def read_geom(input_file):
#
   """
   Function to read geometries from xyz file

   :input_file: file with geometry in xyz format
   """
#
   with open(input_file, 'r') as infile:
      n_atoms = int(infile.readline())
      infile.readline()
#   
      atoms = []
      x = [] 
      y = [] 
      z = []
#   
      for line in infile:
         atoms.append(line.split()[0])
         x.append(float(line.split()[1]))
         y.append(float(line.split()[2]))
         z.append(float(line.split()[3]))
#   
   return(n_atoms,atoms,x,y,z)


# ------------------------------- #
# ------ Untar nanofq file ------ #

def untar_nanofq_file(tar_file):
#
   """
   Function to untar nanofq file 

   :tar_file: tar nanofq file in bem_tommaso_ief format
   """
#
#  Untar the file
   if os.path.exists('tmp_analysis'): os.system('rm -rf tmp_analysis')
   os.system('mkdir tmp_analysis')
#   
   os.system('cp ' + tar_file + ' tmp_analysis/')
   os.system('tar -xf ./tmp_analysis/' + tar_file + ' -C tmp_analysis') 
#
#  Check number of frequency files
   freq_files = glob.glob('tmp_analysis/*freq')
   if (len(freq_files)) > 1:
      print(' ')
      print('  STOP: only 1 frequency file allowed.\n\n'  +
            '        There are ' + str(len(freq_files)) + ' frequency files in "'+ tar_file + '"')
      print(' ')
      sys.exit()


# ------------------------------- #
# ------ Read untared file ------ #

def read_untared_nanofq_file():
#
   """
   Function to untar nanofq file

   :tar_file: tar nanofq file in bem_tommaso_ief format
   """
#
#  Initialize variables
#
   read_x = False
   read_y = False
   read_z = False
#
   lcharges_and_dipoles = False
#
#  Read frequency file
#
   freq_file = glob.glob('tmp_analysis/*freq')[0]
#   
   with open(freq_file,'r') as f:
      # Take 1st line and check if we have charges or charges+dipoles
      header  = f.readline().strip()
#
      if (header.startswith('complex Charges/Dipoles:')):
         lcharges_and_dipoles = True
      else:
         print(' ')
         print('  STOP: only charges + dipoles (wFQFMu) is allowed.')
         print(' ')
         sys.exit()
#   
      # Rewind back to the start of the file
      f.seek(0)
#   
      # Check number of lines and guess number of atoms
      lines = f.readlines()
      n_lines = len(lines)
#
      if (lcharges_and_dipoles): 
         nAtoms = int((n_lines - 3) / (3*4))
      else:
         nAtoms = int((n_lines - 3) / (3))
#   
      # Read charges or charges + dipoles
      if (lcharges_and_dipoles):
#
         data_re_x = []
         data_im_x = []
#
         data_re_y = []
         data_im_y = []
#
         data_re_z = []
         data_im_z = []
#
         for line in lines:
            if (len(line.split())!=3 and (len(line.split())!=6)):
               print(' ')
               print('  STOP: freq file "' + freq_file + '" corrupt')
               print(' ')
               sys.exit()
#
            if (read_x):
               if len(line.split()) !=6:
                  data_re_x.append(float(line.split()[1]))
                  data_im_x.append(float(line.split()[2]))

            elif (read_y):
               if len(line.split()) !=6:
                  data_re_y.append(float(line.split()[1]))
                  data_im_y.append(float(line.split()[2]))

            elif (read_z):
               if len(line.split()) !=6:
                  data_re_z.append(float(line.split()[1]))
                  data_im_z.append(float(line.split()[2]))
#
            if line.startswith('complex Charges/Dipoles: x. freq'):
               read_x = True
               read_y = False
               read_z = False
            elif line.startswith('complex Charges/Dipoles: y. freq'):
               read_x = False
               read_y = True
               read_z = False
            elif line.startswith('complex Charges/Dipoles: z. freq'):
               read_x = False
               read_y = False
               read_z = True
#
      else: 
         print(' ')
         print('  STOP: only charges + dipoles (wFQFMu) is allowed.')
         print(' ')
         sys.exit()

#
   return(lcharges_and_dipoles,nAtoms,data_re_x,data_im_x,data_re_y,data_im_y,data_re_z,data_im_z)
   
   
# ------------------------------- #
# ------ Read np geometry ------- #

def read_np_geometry_and_params(nAtoms):
#
   """
   Function to read np geometry

   nAtoms: number of atoms
   """
#
#  Initialize variables
#
   x_np  = []
   y_np  = []
   z_np  = []
   Chi   = []
   Eta   = []
   Alpha = []
   R_q   = []
   R_mu  = []
#
   read_geometry     = False
#
#  Read info file
#
   info_file = glob.glob('tmp_analysis/*info')[0]
#
   with open(info_file,'r') as f:
      for line in f:
#
         if line.startswith('Parameters:'): read_geometry = False 
#
         if read_geometry:
            x_np.append(float(line.split()[2]))
            y_np.append(float(line.split()[3]))
            z_np.append(float(line.split()[4]))
#
         if line.startswith('Input Geometry (angstrom)'): read_geometry = True
#
      # Check
      if len(x_np) != nAtoms:
         print(' ')
         print('  STOP: FREQ filenumber of atoms ' + str(nAtoms) + 'different from INFO file number of atoms' + str(len(x_np)))
         print(' ')
         sys.exit()
#
      # Rewind back to the start of the file
      f.seek(0)
#
      found_param = False
      for line in f:
#
         if found_param:
            if line.startswith('Number of atoms per IMol'): 
               break
            elif len(line.split())>0:
               Chi.append(float(line.split()[1])) 
               Eta.append(float(line.split()[2]))
               Alpha.append(float(line.split()[3]))
               R_q.append(float(line.split()[4]))
               R_mu.append(float(line.split()[5]))
# 
         if line.startswith('Parameters:'):
            model = line.split()[1]
            found_param = True
            if model != 'fqfmu_pqeq':
               print(' ')
               print('  STOP: only "fqfmu_pqeq" supported. Found "' + model + '" in ' + info_file)
               print(' ')
               sys.exit()
#
      return(Chi,Eta,Alpha,R_q,R_mu,x_np,y_np,z_np)


# ------------------------------- #
# ------ Calculate e_field ------ #

def calculate_electric_field(incident_field_intenstiy,lcharges_and_dipoles,nAtoms,direction,R_q,R_mu,x_np,y_np,z_np,x_grid,y_grid,z_grid,omega_re,omega_im):
#
   """
   Function to calculate electric field

   :lcharges_and_dipoles: False = Charges; True = Charges + Dipoles
   :nAtoms              : Number of atoms
   :direction           : direction of incident electric field
   :omega_re/im_x/y/z   : charges or charges + dipoles (real/imaginary, with incident field along x/y/z
   """
#
#  Initialize variables
#
   one   = 1.0
   two   = 2.0
   three = 3.0
   four  = 4.0
   five  = 5.0
   to_bohr = 1.8897261254578281 # Angstroms to bohr
#
   e_field   = []
#
#  Read charges (dipoles)
#
   ChargeRe = [omega_re[i] for i in range(nAtoms)]
   ChargeIm = [omega_im[i] for i in range(nAtoms)]
   if(lcharges_and_dipoles):
      DipoleRe = [omega_re[i+nAtoms] for i in range(3*nAtoms)]
      DipoleIm = [omega_im[i+nAtoms] for i in range(3*nAtoms)]
#
#  Calculate electric field as function of the direction
#
   for x,y,z in zip(x_grid,y_grid,z_grid):
#
      EFieldX_Im = 0.0
      EFieldX_Re = 0.0
      EFieldY_Im = 0.0
      EFieldY_Re = 0.0
      EFieldZ_Im = 0.0
      EFieldZ_Re = 0.0
#
      for l in range(nAtoms):
#   
         x_IJ   = (x_np[l]-x)*to_bohr
         y_IJ   = (y_np[l]-y)*to_bohr
         z_IJ   = (z_np[l]-z)*to_bohr
#
         distIJ = math.sqrt(x_IJ**2 + y_IJ**2 + z_IJ**2)
         distIJ_2 = distIJ**two
         distIJ_3 = distIJ**three
         distIJ_5 = distIJ**five
#
         R_Q_I  = R_q[l]
         R_Q_I_2 = R_Q_I**two
#
         const = two/(math.sqrt(np.pi)*R_Q_I)
#
         factor = math.erf(distIJ/R_Q_I)-const*distIJ*math.exp(-distIJ_2/R_Q_I_2)
#
         if(distIJ < 0.1):
            print(' ')
            print('  STOP: a NP atom is too close to a point grid.')
            print(' ')
            sys.exit()
#
         if (direction=='x'):
            EFieldX_Im = EFieldX_Im - factor*x_IJ*ChargeIm[l]/distIJ_3
            EFieldX_Re = EFieldX_Re - factor*x_IJ*ChargeRe[l]/distIJ_3
         elif (direction=='y'):
            EFieldY_Im = EFieldY_Im - factor*y_IJ*ChargeIm[l]/distIJ_3
            EFieldY_Re = EFieldY_Re - factor*y_IJ*ChargeRe[l]/distIJ_3
         elif (direction=='z'):
            EFieldZ_Im = EFieldZ_Im - factor*z_IJ*ChargeIm[l]/distIJ_3
            EFieldZ_Re = EFieldZ_Re - factor*z_IJ*ChargeRe[l]/distIJ_3
#
         if(lcharges_and_dipoles):
#
            R_Mu_I   = R_mu[l]
            R_Mu_I_2 = R_Mu_I**two
            R_Mu_I_3 = R_Mu_I**three
#
            fD1 = math.erf(distIJ/R_Mu_I) - ((two*distIJ)/(math.sqrt(np.pi)*R_Mu_I))*math.exp(-distIJ_2/R_Mu_I_2)
            fD2 = (four/(math.sqrt(np.pi)*R_Mu_I_3)) * math.exp(-distIJ_2/R_Mu_I_2)
#
            index_1 = 3*l
#
            if direction=='x':
#
               EFIeldX_Im = EFieldX_Im + ( ((-three*x_IJ*x_IJ)/distIJ_5 + one/distIJ_3)*fD1*DipoleIm[index_1] +
                                              fD2*(x_IJ*x_IJ*DipoleIm[index_1])/distIJ_2                      +
                                           ((-three*x_IJ*y_IJ)/distIJ_5 )*fD1*DipoleIm[index_1+1]             +
                                              fD2*(x_IJ*y_IJ*DipoleIm[index_1+1])/distIJ_2                    +
                                           ((-three*x_IJ*z_IJ)/distIJ_5 )*fD1*DipoleIm[index_1+2]             +
                                              fD2*(x_IJ*z_IJ*DipoleIm[index_1+2])/distIJ_2 )
#
               EFIeldX_Re = EFieldX_Re + ( ((-three*x_IJ*x_IJ)/distIJ_5 + one/distIJ_3)*fD1*DipoleRe[index_1] +
                                              fD2*(x_IJ*x_IJ*DipoleRe[index_1])/distIJ_2                      +
                                           ((-three*x_IJ*y_IJ)/distIJ_5 )*fD1*DipoleRe[index_1+1]             +
                                              fD2*(x_IJ*y_IJ*DipoleRe[index_1+1])/distIJ_2                    +
                                           ((-three*x_IJ*z_IJ)/distIJ_5 )*fD1*DipoleRe[index_1+2]             +
                                              fD2*(x_IJ*z_IJ*DipoleRe[index_1+2])/distIJ_2 )
#
            elif direction=='y':
#
               EFIeldY_Im = EFieldY_Im + ( ((-three*y_IJ*y_IJ)/distIJ_5 + one/distIJ_3)*fD1*DipoleIm[index_1+1] +
                                              fD2*(y_IJ*y_IJ*DipoleIm[index_1+1])/distIJ_2                      +
                                           ((-three*y_IJ*x_IJ)/distIJ_5 )*fD1*DipoleIm[index_1]                 +
                                              fD2*(y_IJ*x_IJ*DipoleIm[index_1])/distIJ_2                        +
                                           ((-three*y_IJ*z_IJ)/distIJ_5 )*fD1*DipoleIm[index_1+2]               +
                                              fD2*(y_IJ*z_IJ*DipoleIm[index_1+2])/distIJ_2 )
#
               EFIeldY_Re = EFieldY_Re + ( ((-three*y_IJ*y_IJ)/distIJ_5 + one/distIJ_3)*fD1*DipoleRe[index_1+1] +
                                              fD2*(y_IJ*y_IJ*DipoleRe[index_1+1])/distIJ_2                      +
                                           ((-three*y_IJ*x_IJ)/distIJ_5 )*fD1*DipoleRe[index_1]                 +
                                              fD2*(y_IJ*x_IJ*DipoleRe[index_1])/distIJ_2                        +
                                           ((-three*y_IJ*z_IJ)/distIJ_5 )*fD1*DipoleRe[index_1+2]               +
                                              fD2*(y_IJ*z_IJ*DipoleRe[index_1+2])/distIJ_2 )
#
            elif direction=='z':
#
               EFIeldZ_Im = EFieldZ_Im + ( ((-three*z_IJ*z_IJ)/distIJ_5 + one/distIJ_3)*fD1*DipoleIm[index_1+2] +
                                              fD2*(z_IJ*z_IJ*DipoleIm[index_1+2])/distIJ_2                      +
                                           ((-three*z_IJ*x_IJ)/distIJ_5 )*fD1*DipoleIm[index_1]                 +
                                              fD2*(z_IJ*x_IJ*DipoleIm[index_1])/distIJ_2                        +
                                           ((-three*z_IJ*y_IJ)/distIJ_5 )*fD1*DipoleIm[index_1+1]               +
                                              fD2*(z_IJ*y_IJ*DipoleIm[index_1+1])/distIJ_2 )
#
               EFIeldZ_Re = EFieldZ_Re + ( ((-three*z_IJ*z_IJ)/distIJ_5 + one/distIJ_3)*fD1*DipoleRe[index_1+2] +
                                              fD2*(z_IJ*z_IJ*DipoleRe[index_1+2])/distIJ_2                      +
                                           ((-three*z_IJ*x_IJ)/distIJ_5 )*fD1*DipoleRe[index_1]                 +
                                              fD2*(z_IJ*x_IJ*DipoleRe[index_1])/distIJ_2                        +
                                           ((-three*z_IJ*y_IJ)/distIJ_5 )*fD1*DipoleRe[index_1+1]               +
                                              fD2*(z_IJ*y_IJ*DipoleRe[index_1+1])/distIJ_2 )
#
      e_field.append(math.sqrt(EFieldX_Im**2 + EFieldX_Re**2 +
                               EFieldY_Im**2 + EFieldY_Re**2 +
                               EFieldZ_Im**2 + EFieldZ_Re**2)/incident_field_intenstiy)
#
   return(e_field)


# ------------------------------- #
# -------- Plot e_field --------- #

def plot_electric_field(x,y,z,x_np,y_np,z_np,e_field):
#
   """
   Function to calculate electric field

   :x,y,z   : Coordinates at which the electric field was calculated
   :e_field : Electric field at each coordinate point
   """
#
   # Creating a 3D scatter plot
   fig = plt.figure()
   ax = fig.add_subplot(111, projection='3d')
   
   # Scatter plot using x, y, z for position and e_field for color
   scatter = ax.scatter(x_grid, y_grid, z_grid, c=e_field)

   ## Adding spheres
   #for x_np_, y_np_, z_np_ in zip(x_np, y_np, z_np):
   #    ax.scatter(x_np_, y_np_, z_np_, color='yellow') #, s=1)  # s is the size of the point

   
   # Adding a color bar
   plt.colorbar(scatter, ax=ax, label='Electric Field Enhancement')
   
   # Adding labels
   ax.set_xlabel('X axis')
   ax.set_ylabel('Y axis')
   ax.set_zlabel('Z axis')
   
   # Show the plot
   plt.show()


# -----

          
# ==============================
#         MAIN PROGRAM 
# ==============================


# START TIMER
start = time.time()


# --> Read command line and check inputs

tar_file, grid_file,direction = read_command_line(sys.argv)

check_inputs(tar_file,grid_file,direction)

# --> Untar tar file and read grid_file

untar_nanofq_file(tar_file)

n_points, grid_atoms, x_grid,y_grid,z_grid = read_geom(grid_file)

# --> Read charges/dipoles and nanoparticle geometry + model parameters

lcharges_and_dipoles,nAtoms,omega_re_x,omega_im_x,omega_re_y,omega_im_y,omega_re_z,omega_im_z = read_untared_nanofq_file()

Chi,Eta,Alpha,R_q,R_mu,x_np,y_np,z_np = read_np_geometry_and_params(nAtoms)

# --> Calculate Electric Field Enhancement

if direction=='x': e_field = calculate_electric_field(incident_field_intenstiy,lcharges_and_dipoles,
                                                      nAtoms,direction,R_q,R_mu,x_np,y_np,z_np,
                                                      x_grid,y_grid,z_grid,omega_re_x,omega_im_x)

if direction=='y': e_field = calculate_electric_field(incident_field_intenstiy,lcharges_and_dipoles,
                                                      nAtoms,direction,R_q,R_mu,x_np,y_np,z_np,
                                                      x_grid,y_grid,z_grid,omega_re_y,omega_im_y)

if direction=='z': e_field = calculate_electric_field(incident_field_intenstiy,lcharges_and_dipoles,
                                                      nAtoms,direction,R_q,R_mu,x_np,y_np,z_np,
                                                      x_grid,y_grid,z_grid,omega_re_z,omega_im_z)

# --> 3D plot of electric field

plot_electric_field(x_grid,y_grid,z_grid,x_np,y_np,z_np,e_field)





# END TIMER
end = time.time()


#print('')
#print('   ----------------------------------')
#print('   NORMAL TERMINATION')
#print('')
#print('   COMPUTATIONAL TIME: ', str(round(end-start,2)), ' s')
#print('')
#print('')



# ==============================
#         MAIN PROGRAM 
# ==============================








