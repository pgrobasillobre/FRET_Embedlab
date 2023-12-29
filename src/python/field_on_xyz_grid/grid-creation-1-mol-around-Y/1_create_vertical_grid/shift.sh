#!/bin/bash
counter=0
for i in $(seq 14 -2.0 -16.0); do
   python3 /Users/pablogrobasillobre/Dropbox/3rd-year/research/sveva-fret/programs/fret_embedlab_fortran/src/python/geom.py -t1 $counter  h1-up.xyz origin_CM_no -y
   counter=$(echo "$counter + 2.0" | bc)
done
