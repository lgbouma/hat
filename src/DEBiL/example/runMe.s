#!/bin/csh -f

echo " "
echo "Finding periods:"
echo " "
                
foreach file (*.dat)
echo $file
periodS2 $file periods.txt periods.err
end

cut -d' ' -f1,2 periods.txt >! list.txt

echo " "
echo "Running DEBiL:"
echo " "
                
debil_OGLE list.txt list.dbl list.err

echo " "
echo "Results:"
echo " "
cat list.dbl
