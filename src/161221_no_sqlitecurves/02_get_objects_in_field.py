'''
Having copied select data from our chosen HAT field (with BLS already run on 
them) to a local directory, we want to get the names of those objects, and 
then output a text file with the HAT IDs and corresponding 2MASS IDs.
'''

import numpy as np
import os
import json, requests

data_dir = '../data/HATpipe/'
field = 'G257' #FIXME generalize for any field

all_bls_files = os.listdir(data_dir+field+'/BLS')

### Check that tails are what we expect
tails = [f[-8:] for f in all_bls_files]
tails_u = np.unique(tails)
expected_tails = ['.blsanal', '.blsspec']
for el in tails_u:
    assert el in expected_tails, \
    '{:s} is unexpected filetype in {:s}/BLS'.format(el, field)

# Get list of HAT object IDs in field
hat_ids = np.unique([f[:-8] for f in all_bls_files])

# Find corresponding 2MASS IDs.
hat_url = 'https://hatsurveys.org/lightcurves/object/'
api_key = '?apikey=ZjAwYWY1NWU3NDlhYmU0ZWQyZTU3NGM0YmM4ZmExZDYyOTEzNmE5MTQ'+\
          'yZTdkN2JlMDM3YzkzMDUyNmVjYTU0Zg'
hat_url = 'https://hatsurveys.org/lightcurves/object/'

twomass_ids = []
for ix, hat_ID in enumerate(hat_ids):
    query_url = hat_url + hat_ID + api_key
    resp = requests.get(url=query_url)
    data = json.loads(resp.text)

    twomass_id = data['results']['twomassid']

    twomass_ids.append(twomass_id)
    if ix % 500 == 0:
        print(ix)

out = np.array([hat_ids, twomass_ids])
data_path = '../data/intermed/'
out_fname = '02_'+field+'_hatid_twomassid.out'
np.savetxt(data_path+out_fname, out.T, fmt='%s', delimiter=',')
