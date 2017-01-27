#!/bin/bash

echo "Generating checkplotlist for server to read..."
python 03_look_at_LCs.py

wait $(jobs -p)

echo "Launching checkplotserver..."
#want to launch from proj/hat/src/

python checkplotserver.py
