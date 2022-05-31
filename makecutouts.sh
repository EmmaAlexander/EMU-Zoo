#!/bin/bash 
for i in 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150
do
   echo $i
   python emu-zoo_make_alpha_cutouts.py 0 $i
done