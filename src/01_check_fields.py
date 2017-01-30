'''
We're going to have lots of HATN/S fields given to us from Waqas. We need to
make sure that everything we're given (listed in 01_DATA_waqas_fields.txt) wind
up making it to /data/sqlitecurve_headers.

Usage:
1) Enter field numbers given from Waqas Bhatti.
2) Run:
    >>> python 01_check_fields.py
'''
import numpy as np
import os

f = open('01_DATA_waqas_fields.txt')
lines = f.readlines()
f.close()
given_fields = np.sort(np.unique([l.strip('\n') for l in lines[1:]]))

sqlitecurve_fields = np.sort(np.array([f.strip('-sqlitecurves.csv')[1:] for f
    in os.listdir('../data/sqlitecurve_headers/')
    if f.endswith('-sqlitecurves.csv')]))

nonint = list(set(given_fields).symmetric_difference(sqlitecurve_fields))

if len(nonint)>=1:
    print('The following fields never made it to sqlitecurve_fields:')
    print(nonint)
else:
    np.testing.assert_array_equal(sqlitecurve_fields, given_fields)
    print('All fields entered in 01_DATA_waqas_fields.txt are in'+\
        'sqlitecurve_fields')

print('unique field are:')
print(sqlitecurve_fields)
print('with length:')
print(len(sqlitecurve_fields))
