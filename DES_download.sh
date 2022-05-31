#!/bin/bash

# Input and output files
DES_FILE_INDEX_LOC='/Volumes/TARDIS/Work/EMUzoo/DES_tiles.txt'
DES_FILE_SUB_INDEX_LOC='/Volumes/TARDIS/Work/EMUzoo/DES_tiles_downloadtest.txt'
DES_DATA_LOC='/Volumes/TARDIS/Work/EMUzoo/data/DES/'

# Limits of sky area
# Take care when using negative declinations!! 
# Also remember to not use spaces
RA_MIN='2030' #HHMM
RA_MAX='2040' #HHMM
DEC_MIN='-5100' #DDMM 
DEC_MAX='-5030' #DDMM

# define DES bands to download 
#choose all or subset of ('Y' 'g' 'i' 'r' 'z')
declare -a bands=('g' 'i' 'r')

# make the file list
python DES_file_select.py $DES_FILE_INDEX_LOC $RA_MIN $RA_MAX $DEC_MIN $DEC_MAX $DES_FILE_SUB_INDEX_LOC
# /mnt/shared/home/ealexander/scripts/DES_file_select.py

# download the files
input=$DES_FILE_SUB_INDEX_LOC
#loop over DES tiles by reading file line by line
while IFS= read -r line
do
  echo "Downloading $line files"
  # loop over bands
  for b in "${bands[@]}"
  do
	  echo "$b"
	  curl -f "https://desdr-server.ncsa.illinois.edu/despublic/dr2_tiles/${line}/${line}_r4575p02_${b}.fits.fz" -L -o "${DES_DATA_LOC}${line}_r4575p02_${b}.fits.fz"
	  curl -f "https://desdr-server.ncsa.illinois.edu/despublic/dr2_tiles/${line}/${line}_r4575p01_${b}.fits.fz" -L -o "${DES_DATA_LOC}${line}_r4575p01_${b}.fits.fz"
  done
done < "$input"