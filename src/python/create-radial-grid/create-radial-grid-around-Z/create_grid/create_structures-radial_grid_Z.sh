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
geom=/home/pablo/Dropbox/3rd-year/research/sveva-fret/programs/fret_embedlab_fortran/src/python/geom.py
#
np_name=cto_r30 # the np geometry cannot be stored in a folder
#
molecule_up=init-geom/up-y.xyz # do not change the name of the folder, you can only change the name of the xyz file
molecule_down=init-geom/down-y.xyz # do not change the name of the folder, you can only change the name of the xyz file
#
minimum_molecule_np_distance=5.0
#
angle_step=60 # Define the step size for angle

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


# Generate angles for Z axis
if [ -d "input_angle_z" ]; then
    rm -rf "input_angle_z"
fi
angles_z=()
for ((i=0; i<360; i+=angle_step)); do
    angles_z+=($i)
done

# Write angles to files
for angle in "${angles_z[@]}"; do
    echo "$angle" >> input_angle_z
done


# Rotate the np around the z axis
python3 $geom -r input_angle_z $np_name.xyz origin_CM_no +z > tmp
rm -rf tmp
#
echo " "
echo "   Z-NANOPARTICLE ROTATED"
echo " "
#
mv results_rotate rotated_np_z
   
#  ok 1
   
# Traslate all the fret aceptor and donor files respect the rotated np
# we expect to have ONE MOLECULE UP IN THE Y AXIS AND ONE MOLECULE
# DOWN IN THE Y AXIS
if [ -d "molecules" ]; then
    rm -rf "molecules"
fi
mkdir molecules


# MANAGE FIRST MOLECULE --> UP IN THE Y AXIS

for file in rotated_np_z/*.xyz; do

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

mv final_molecules_z/*xyz "grid/up/${counter}"

rm -rf final_molecules_z

rm -rf rotated_np_* input_*

cd grid/up

mkdir all
cp *xyz all 
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


# Generate angles for Z axis
if [ -d "input_angle_z" ]; then
    rm -rf "input_angle_z"
fi
angles_z=()
for ((i=0; i<360; i+=angle_step)); do
    angles_z+=($i)
done
angles_z+=(360)  # Add the final angle 180

# Write angles to files
for angle in "${angles_z[@]}"; do
    echo "$angle" >> input_angle_z
done


# Rotate the np around the z axis
python3 $geom -r input_angle_z $np_name.xyz origin_CM_no +z > tmp
rm -rf tmp
#
echo " "
echo "   Z-NANOPARTICLE ROTATED"
echo " "
#
mv results_rotate rotated_np_z
   
#  ok 1
   
# Traslate all the fret aceptor and donor files respect the rotated np
# we expect to have ONE MOLECULE UP IN THE Y AXIS AND ONE MOLECULE
# DOWN IN THE Y AXIS
if [ -d "molecules" ]; then
    rm -rf "molecules"
fi
mkdir molecules


# MANAGE SECOND MOLECULE --> DOWN IN THE Y AXIS

for file in rotated_np_z/*.xyz; do

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

mv final_molecules_z/*xyz "grid/down/${counter}"

rm -rf final_molecules_z

rm -rf rotated_np_* input_*

cd grid/down

mkdir all
cp *xyz all 
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



