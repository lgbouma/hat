#!/bin/bash

echo " "
echo "Finding periods:"
echo " "
                
for file in $( ls *dat ); do
  echo "item: $file"
  ./periodS2 $file periods.txt periods.err
  cut -d' ' -f1,2 periods.txt > list.txt
  ./debil list.txt list.dbl list.err
done

#for f in $files
#do
#  echo $file
#  ./periodS2 $file periods.txt periods.err
#
#  cut -d' ' -f1,2 periods.txt >! list.txt
#
#  echo " "
#  echo "Running DEBiL:"
#  echo " "
#                  
#  ./debil_OGLE list.txt list.dbl list.err
#
#  echo " "
#  echo "Results:"
#  echo " "
#  cat list.dbl
#done
