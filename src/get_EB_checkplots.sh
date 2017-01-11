#!/bin/bash

DSP_lim=20

mapfile -t fieldIDs < INPUT_field_ID_numbers.txt

for id in "${fieldIDs[@]}"
do
  outfile="LOGS/G"$id"_DSP"$DSP_lim".out"
  echo "Running G"$id". Output: "$outfile

  printf "INPUT_field_ID_numbers.txt:\n" > $outfile
  cat INPUT_field_ID_numbers.txt >> $outfile

  printf "\nget_EB_checkplots.sh:\n" >> $outfile
  cat get_EB_checkplots.sh >> $outfile
  
  printf "\nExecute:" >> $outfile
  printf "\npython 02_get_EB_checkplots.py $id $DSP_lim >> $outfile &" >> $outfile
  printf "\n==============================\n" >> $outfile

  python 02_get_EB_checkplots.py $id $DSP_lim >> $outfile &

  wait $(jobs -p) # Wait for just-submitted job to complete before proceeding.
done
