#!/bin/bash

# ================ IMPORTANT NOTES ================
#
# We assume the reference initial position of the NP
# align with the Y axis.
#
# This means that initially the molecule is moved up
# or down respect that axis, and then a set of 
# translations and rotations are applied, according 
# to the options given (see below)
#
# EXECUTION: ./create_structures_full_grid.sh > log
#
# ERRORS: to check errors do (any match is an error):
#         $ grep -irs 'error' log
#         $ grep -irs 'stop' log
#         $ grep -irs 'traceback' log
#
# RESULTS: grid/up/grid-up.xyz 
#          grid/down/grid-down.xyz
#
# =================================================
# ==================== OPTIONS ====================
# =================================================

# WHERE IS THE GEOM.PY CODE? --> FULL PATH
geom=/Users/pablogrobasillobre/Dropbox/3rd-year/research/sveva-fret/programs/fret_embedlab_fortran/src/python/geom.py
#
np_name=cto_r30 # the np geometry cannot be stored in a folder
#
molecule_up=init-geom/up-y.xyz # do not change the name of the folder, you can only change the name of the xyz file
molecule_down=init-geom/down-y.xyz # do not change the name of the folder, you can only change the name of the xyz file
#
minimum_molecule_np_distance=5.0
#
angle_step=45 # Define the step size for angle

#
# =================================================
# =================== BASH CODE ===================
# =================================================

if [ -d "grid" ]; then
    rm -rf "grid"
fi
mkdir grid


# =================================================
# ================== MOLECULE UP ==================
# =================================================

echo " "
echo "   MOLECULE UP"
echo " "
#
mkdir grid/up

# Extract number of lines in the file specified by $molecule_up
nlines_molecule_up=$(wc -l < "$molecule_up")
natoms_molecule_up=$((nlines_molecule_up - 2))

# Define input for minimum molecule to nanoparticle distance
if [ -d "input_translate" ]; then
    rm -rf "input_translate"
fi
echo "$minimum_molecule_np_distance" > input_translate


# Generate angles for Y axis
if [ -d "input_angle_y" ]; then
    rm -rf "input_angle_y"
fi
if [ -d "input_angle_z" ]; then
    rm -rf "input_angle_z"
fi

angles_y=()
for ((i=0; i<360; i+=angle_step)); do
    angles_y+=($i)
done

# Generate angles for Z axis
angles_z=()
for ((i=0; i<180; i+=angle_step)); do
    angles_z+=($i)
done

# Write angles to files
for angle in "${angles_y[@]}"; do
    echo "$angle" >> input_angle_y
done

for angle in "${angles_z[@]}"; do
    echo "$angle" >> input_angle_z
done



# Rotate the np around the y axis
python3 $geom -r input_angle_z $np_name.xyz origin_CM_no +y > tmp
rm -rf tmp
#
echo " "
echo "   Y-NANOPARTICLE ROTATED"
echo " "
#
mv results_rotate rotated_np_y
#
counter=0
for rotated_np_y in rotated_np_y/*.xyz; do
#
   counter=$(echo "$counter + 1" | bc)
   mkdir "grid/up/"${counter}
#
   # Rotate the np around z axis
   python3 $geom -r input_angle_y $rotated_np_y origin_CM_no +z > tmp
   
   rm -rf tmp
   
   echo " "
   echo "   Z-NANOPARTICLE ROTATED"
   echo " "
   
   mv results_rotate rotated_np_y_z
   
#  ok 1
   
   # Traslate all the fret aceptor and donor files respect the rotated np
   # we expect to have ONE MOLECULE UP IN THE Y AXIS AND ONE MOLECULE
   # DOWN IN THE Y AXIS
   if [ -d "molecules" ]; then
       rm -rf "molecules"
   fi
   mkdir molecules
   
  
   # MANAGE FIRST MOLECULE --> UP IN THE Y AXIS
   
   for file in rotated_np_y_z/*.xyz; do

       # Translate
       python3 $geom -t input_translate $file origin_CM_1_no $molecule_up origin_CM_2_no +y verbose_yes

       translated_molecule_up=${file:14}

       mv results_translate/*xyz molecules/$translated_molecule_up

       rm -rf results_translate
   
   done
#  ok 2

   # Now rotate the translated molecules around Z
   if [ -d "final_molecules_z" ]; then
       rm -rf "final_molecules_z"
   fi
   mkdir final_molecules_z
   for file in molecules/*.xyz; do
   
       #Extract the part of the filename before .xyz
       filename_without_extension="${file%.xyz}"
   
       # Extract the number after the last underscore
       last_number="${filename_without_extension##*_}"
   
       # Print the extracted number
       echo "$last_number" > tmp_angle
   
       # Rotate the np
       python3 $geom -r tmp_angle $file origin_CM_no +z
   
       mv results_rotate/*xyz final_molecules_z
       rm -rf results_rotate
       rm -rf tmp_angle
   
   done
#  ok 3

   # Now rotate the translated molecules around Y
   if [ -d "final_molecules_z_y" ]; then
       rm -rf "final_molecules_z_y"
   fi
   mkdir final_molecules_z_y
   for file in final_molecules_z/*.xyz; do
   

       # Extract the part after "${np_name}_+y_degree_" and before the next underscore
       prefix="${np_name}_+y_degree_"
       temp=${file#*$prefix}
       y_degree=${temp%%_*}

       # Do whatever you need with y_degree
       echo "Extracted degree: $y_degree"

       # Print the extracted number
       echo "$y_degree" > tmp_angle

       # Rotate the np
       python3 $geom -r tmp_angle $file origin_CM_no +y

       mv results_rotate/*xyz final_molecules_z_y
       rm -rf results_rotate
       rm -rf tmp_angle
   
   done

   
   mv final_molecules_z_y/*xyz "grid/up/${counter}"

   rm -rf final_molecules_z rm -rf final_molecules_z_y molecules rotated_np_y_z

done # rotated np y

rm -rf rotated_np_y input_*

# Only consider the initial position and position at 180 degrees (Z axis), 
# otherwise we have overlapping issues within that specific point

cd grid/up
for i in *; do
   if [ "$i" -gt 1 ]; then  # Check if the folder name as a number is greater than 1
      find "$i" -type f -name '*degree*_0.0_*degree*_0.0_*' -exec rm {} +
      find "$i" -type f -name '*degree*_180.0_*degree*_180.0_*' -exec rm {} +
   fi
done

mkdir ../all
cp */* ../all 
mv ../all .

cd all

for i in *; do tail -$natoms_molecule_up  $i | cat >> tmp_xyz ; done 

wc -l tmp_xyz | awk '{print $1}' > tmp_n_points

echo "XYZ Grid Points" > tmp_title

cat tmp_n_points tmp_title tmp_xyz > grid-up.xyz

mv grid-up.xyz ../
rm -rf tmp*

cd ../../../


echo " "
echo " GRID CREATED --> check grid/up/grid-up.xyz"
echo " "



# =================================================
# ================= MOLECULE DOWN =================
# =================================================


echo " "
echo "   MOLECULE DOWN"
echo " "
#
mkdir grid/down

# Extract number of lines in the file specified by $molecule_down
nlines_molecule_down=$(wc -l < "$molecule_down")
natoms_molecule_down=$((nlines_molecule_down - 2))

# Define input for minimum molecule to nanoparticle distance
if [ -d "input_translate" ]; then
    rm -rf "input_translate"
fi
echo "$minimum_molecule_np_distance" > input_translate


# Generate angles for Y axis
if [ -d "input_angle_y" ]; then
    rm -rf "input_angle_y"
fi
if [ -d "input_angle_z" ]; then
    rm -rf "input_angle_z"
fi

angles_y=()
for ((i=0; i<360; i+=angle_step)); do
    angles_y+=($i)
done

# Generate angles for Z axis
angles_z=()
for ((i=0; i<180; i+=angle_step)); do
    angles_z+=($i)
done

# Write angles to files
for angle in "${angles_y[@]}"; do
    echo "$angle" >> input_angle_y
done

for angle in "${angles_z[@]}"; do
    echo "$angle" >> input_angle_z
done



# Rotate the np around the y axis
python3 $geom -r input_angle_z $np_name.xyz origin_CM_no +y > tmp
rm -rf tmp
#
echo " "
echo "   Y-NANOPARTICLE ROTATED"
echo " "
#
mv results_rotate rotated_np_y
#
counter=0
for rotated_np_y in rotated_np_y/*.xyz; do
#
   counter=$(echo "$counter + 1" | bc)
   mkdir "grid/down/"${counter}
#
   # Rotate the np around z axis
   python3 $geom -r input_angle_y $rotated_np_y origin_CM_no +z > tmp
   
   rm -rf tmp
   
   echo " "
   echo "   Z-NANOPARTICLE ROTATED"
   echo " "
   
   mv results_rotate rotated_np_y_z
   
#  ok 1
   
   # Traslate all the fret aceptor and donor files respect the rotated np
   # we expect to have ONE MOLECULE down IN THE Y AXIS AND ONE MOLECULE
   # DOWN IN THE Y AXIS
   if [ -d "molecules" ]; then
       rm -rf "molecules"
   fi
   mkdir molecules
   
  
   # MANAGE FIRST MOLECULE --> down IN THE Y AXIS
   
   for file in rotated_np_y_z/*.xyz; do

       # Translate
       python3 $geom -t input_translate $file origin_CM_1_no $molecule_down origin_CM_2_no -y verbose_yes

       translated_molecule_down=${file:14}

       mv results_translate/*xyz molecules/$translated_molecule_down

       rm -rf results_translate
   
   done
#  ok 2

   # Now rotate the translated molecules around Z
   if [ -d "final_molecules_z" ]; then
       rm -rf "final_molecules_z"
   fi
   mkdir final_molecules_z
   for file in molecules/*.xyz; do
   
       #Extract the part of the filename before .xyz
       filename_without_extension="${file%.xyz}"
   
       # Extract the number after the last underscore
       last_number="${filename_without_extension##*_}"
   
       # Print the extracted number
       echo "$last_number" > tmp_angle
   
       # Rotate the np
       python3 $geom -r tmp_angle $file origin_CM_no +z
   
       mv results_rotate/*xyz final_molecules_z
       rm -rf results_rotate
       rm -rf tmp_angle
   
   done
#  ok 3

   # Now rotate the translated molecules around Y
   if [ -d "final_molecules_z_y" ]; then
       rm -rf "final_molecules_z_y"
   fi
   mkdir final_molecules_z_y
   for file in final_molecules_z/*.xyz; do
   

       # Extract the part after "${np_name}_+y_degree_" and before the next underscore
       prefix="${np_name}_+y_degree_"
       temp=${file#*$prefix}
       y_degree=${temp%%_*}

       # Do whatever you need with y_degree
       echo "Extracted degree: $y_degree"

       # Print the extracted number
       echo "$y_degree" > tmp_angle

       # Rotate the np
       python3 $geom -r tmp_angle $file origin_CM_no +y

       mv results_rotate/*xyz final_molecules_z_y
       rm -rf results_rotate
       rm -rf tmp_angle
   
   done

   
   mv final_molecules_z_y/*xyz "grid/down/${counter}"

   rm -rf final_molecules_z rm -rf final_molecules_z_y molecules rotated_np_y_z

done # rotated np y

rm -rf rotated_np_y input_*

# Only consider the initial position and position at 180 degrees (Z axis), 
# otherwise we have overlapping issues within that specific point

cd grid/down
for i in *; do
   if [ "$i" -gt 1 ]; then  # Check if the folder name as a number is greater than 1
      find "$i" -type f -name '*degree*_0.0_*degree*_0.0_*' -exec rm {} +
      find "$i" -type f -name '*degree*_180.0_*degree*_180.0_*' -exec rm {} +
   fi
done

mkdir ../all
cp */* ../all 
mv ../all .

cd all

for i in *; do tail -$natoms_molecule_down  $i | cat >> tmp_xyz ; done 

wc -l tmp_xyz | awk '{print $1}' > tmp_n_points

echo "XYZ Grid Points" > tmp_title

cat tmp_n_points tmp_title tmp_xyz > grid-down.xyz

mv grid-down.xyz ../
rm -rf tmp*

cd ../../../


echo " "
echo " GRID CREATED --> check grid/down/grid-down.xyz"
echo " "




