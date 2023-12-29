#!/bin/bash


# Rotate the np
python3 /Users/pablogrobasillobre/Dropbox/3rd-year/research/sveva-fret/programs/fret_embedlab_fortran/src/python/geom.py -r input_angle cto_r30.xyz origin_CM_no +y

mv results_rotate rotated_np


# Traslate all the fret aceptor (or donor) files respect the rotated np
if [ -d "molecules" ]; then
    rm -rf "molecules"
fi
mkdir molecules

counter=0
if [ -d "grid" ]; then
    rm -rf "grid"
fi
mkdir grid

for molecule in 1_create_vertical_grid/results_translate/*xyz; do

   counter=$(echo "$counter + 1" | bc)
   mkdir "grid/"${counter}

   for file in rotated_np/*.xyz; do

       # Translate
       python3 /Users/pablogrobasillobre/Dropbox/3rd-year/research/sveva-fret/programs/fret_embedlab_fortran/src/python/geom.py -t input_translate $file origin_CM_1_no $molecule origin_CM_2_no -x verbose_yes
   
       rotated_molecule=${file:11}
   
       mv results_translate/*xyz molecules/$rotated_molecule
   
       rm -rf results_translate
   
   done
   
   
   # Now rotate the translated molecules
   if [ -d "final_molecules" ]; then
       rm -rf "final_molecules"
   fi
   mkdir final_molecules
   for file in molecules/*.xyz; do
   
       #Extract the part of the filename before .xyz
       filename_without_extension="${file%.xyz}"
   
       # Extract the number after the last underscore
       last_number="${filename_without_extension##*_}"
   
       # Print the extracted number
       echo "$last_number" > tmp_angle
   
       # Rotate the np
       python3 /Users/pablogrobasillobre/Dropbox/3rd-year/research/sveva-fret/programs/fret_embedlab_fortran/src/python/geom.py -r tmp_angle $file origin_CM_no +y
   
       mv results_rotate/*xyz final_molecules
       rm -rf results_rotate
       rm -rf tmp_angle
   
   done

   mv final_molecules/*xyz "grid/${counter}"

done

rm -rf final_molecules molecules rotated_np


# END ORGANIZING THE RESULTS
cd grid

mkdir ../raw-data
mv * ../raw-data
mv ../raw-data .
cd raw-data

for i in *; do for j in $i/*; do tail -1  $j | cat >> tmp_xyz ; done ; done

wc -l tmp_xyz | awk '{print $1}' > tmp_n_points

echo "XYZ Grid Points" > tmp_title

cat tmp_n_points tmp_title tmp_xyz > grid.xyz

mv grid.xyz ../
rm -rf tmp*

cd ../../


echo " "
echo " GRID CREATED --> check grid/grid.xyz"
echo " "



