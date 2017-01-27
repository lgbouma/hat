'''
Once LCs and periodograms have been created, find some detached eclipsing
binaries through direct application of the human brain's pattern-recognition
algorithm ("visual inspection").

Usage:
Manually enter field_id (if many, list them in 03_DATA_fields_todo.txt).
Then:
   >>> ./03_view_checkplots.sh

Then manually write in the processed field ID to 04_DATA_processed_fields.txt
'''

from astrobase import checkplotlist as cpl
import os
from numpy import all as npall, array as nparray

field_id = '093'

arg0 = 'checkplotlist.py'
arg1 = 'pkl'
cpdir = '../data/CPs_cut/G'+field_id+'_20'
arg2 = cpdir

cpnames = [f for f in os.listdir(cpdir) if 'checkplot' in f]
assert npall(nparray([cpn.endswith('.pkl') for cpn in cpnames])), \
        'extract checkplots to .pkl files!'

cpl.main([arg0,arg1,arg2])
