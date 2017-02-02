'''
Once LCs and periodograms have been created, find some detached eclipsing
binaries through direct application of the human brain's pattern-recognition
algorithm ("visual inspection").

Usage:
Manually enter field_id (e.g., as listed in 03_DATA_fields_todo.txt) in
src/03_view_checkplots.sh.
Then:
   >>> ./03_view_checkplots.sh

In other words, this script never get directly called.

Once happy with visual inspection of a given field, manually write in the
processed field ID to 04_DATA_processed_fields.txt
'''

from astrobase import checkplotlist as cpl
import os, subprocess, argparse
from numpy import all as npall, array as nparray

def run_script(script, stdin=None):
    '''
    Run a bash script specified by the `script` string.
    Returns (stdout, stderr), raises error on non-zero return code. Spaces,
    quotes, and newlines should not raise problems (subprocess.Popen parses
    everything).
    '''
    proc = subprocess.Popen(['bash', '-c', script],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        stdin=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    return stdout, stderr

parser = argparse.ArgumentParser()
parser.add_argument('args', nargs='+')
ao = parser.parse_args()
assert len(ao.args) == 1, 'one arg allowed: field_id'

field_id = ao.args[0]

arg0 = 'checkplotlist.py'
arg1 = 'pkl'
if field_id == 'weird':
    cpdir = '../data/weirdpkls'
elif field_id == 'deb':
    cpdir = '../data/debpkls'
else:
    cpdir = '../data/CPs_cut/G'+field_id+'_20'
arg2 = cpdir

cpnames = [f for f in os.listdir(cpdir) if 'checkplot' in f]

#Extract checkplots .pkl.gz to .pkl files:
if not npall(nparray([cpn.endswith('.pkl') for cpn in cpnames])):
    print('gunzip '+cpdir+'/*pkl.gz')
    stdout, stderr = run_script('gunzip '+cpdir+'/*pkl.gz')

cpl.main([arg0,arg1,arg2])
