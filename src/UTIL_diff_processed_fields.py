'''
Get lists of the fields that I've processed, fields that I've received from
Waqas, and the difference between them.
'''

import numpy as np, os, pandas as pd

with open('01_DATA_waqas_fields.txt') as f:
    receivedfields = np.unique(np.array([l.strip('\n') for l in f.readlines() \
            if not l.startswith('#')]))

with open('04_DATA_processed_fields.txt') as f:
    procfields = np.array([l.strip('\n') for l in f.readlines() \
            if not l.startswith('#')])

assert np.all(np.unique(procfields) == np.sort(procfields)), \
        '04_DATA_processed_fields.txt should not have duplicates.'

#Unit test on processed fields
mvdfields = np.array([d[1:4] for d in \
                os.listdir('/media/luke/LGB_tess_data/hat/CPs_cut/')])

difffields = np.setdiff1d(mvdfields, procfields)

if len(difffields)>=1:
    print(difffields)

    assert np.all(np.sort(procfields) == np.sort(mvdfields)), \
            'Fields listed in 04_DATA_processed_fields.txt should match those'+\
            'actually in the external harddrive.'

else:
    print(
    '''
    Fields on external harddrive are same as those in
    04_DATA_processed_fields.txt
    ''')


#What fields still need to be processed?
waitingfields = np.setdiff1d(receivedfields, procfields)

print('The following fields have yet to be processed:')
for ix, wf in enumerate(waitingfields):
    print('{:d}.\t'.format(ix+1),wf)

print('\n{:d} fields have been processed, {:d} have been received from WB.'.format(\
        len(procfields), len(receivedfields)))
