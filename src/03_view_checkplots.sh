#!/bin/bash

#Once LCs and periodograms have been created, find some detached eclipsing
#binaries through direct application of the human brain's pattern-recognition
#algorithm ("visual inspection").
#
#Usage:
#Manually enter field_id (e.g., as listed in 03_DATA_fields_todo.txt) in
#src/03_view_checkplots.sh.
#Then:
#   >>> ./03_view_checkplots.sh
#
#In other words, this script never get directly called.
#
#Once happy with visual inspection of a given field, manually write in the
#processed field ID to 04_DATA_processed_fields.txt


echo "Generating checkplotlist for server to read..."
python 03_look_at_LCs.py 116

wait $(jobs -p)

echo "Launching checkplotserver..."
#want to launch from proj/hat/src/

python checkplotserver.py
