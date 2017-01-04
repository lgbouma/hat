'''
Given a sqlitecurve "header" (csv file with hatid, twomassid, lcfpath, and ndet
from the HATNet/HATS data), download the sqlitecurve to /data/LCs/ from HAT
servers.
'''

import pandas as pd
import requests, os, time

field_name = 'G199' #FIXME generalize for arbitary input args
head_base = '../data/sqlitecurves/headers/'
head_path = head_base+field_name+'-sqlitecurves.csv' #TODO get WB's full list

head_data = pd.read_csv(head_path)

# Make the directory to write sqlite.gz light curves
write_dir = '../data/LCs/'+field_name
if not os.path.exists(write_dir):
    os.makedirs(write_dir)

url_base = 'http://data.hatsurveys.org/archive/'
k = 0
for lcfpath in head_data['lcfpath']:
    # If you already have the LC, don't redownload it. 
    if os.path.exists(write_dir+'/'+lcfpath.split('/')[-1]):
        continue
    else:
        # Download and write the LC
        this_url = url_base+lcfpath
        r = requests.get(this_url)
        with open(write_dir+'/'+lcfpath.split('/')[-1], 'wb') as data:
                data.write(r.content)
        # Time delay; don't overload the HAT servers
        if k % 100 == 0:
            print(k)
            time.sleep(0.1)
        k += 1
