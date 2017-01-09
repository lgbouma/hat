'''
Given a sqlitecurve "header" (csv file with hatid, twomassid, lcfpath, and ndet
from the HATNet/HATS data), download the sqlitecurve to /data/LCs/ from HAT
servers.

Usage:
    python 02_retrive_LCs_given_sqlitecurve_head.py field_id
Args:
    field_id: integer specifying HAT field ID number (HAT-???-1234567).
'''

import pandas as pd
import requests, os, time
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('args', nargs='+')
    ao = parser.parse_args()
    assert len(ao.args) == 1, 'only one arg: HAT field_id'

    field_name = 'G' + ao.args[0]
    head_base = '../data/sqlitecurve_headers/'
    head_path = head_base+field_name+'-sqlitecurves.csv' #TODO get WB's full list

    head_data = pd.read_csv(head_path)

    # Make the directory to write sqlite.gz light curves
    write_dir = '../data/LCs/'+field_name
    if not os.path.exists(write_dir):
        os.makedirs(write_dir)

    # Let user know what they're doing
    print('Downloading sqlitecurves for {:s} ({:d} total). Update every 100'.\
            format(field_name, len(head_data)-1))

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

if __name__ == '__main__':
    main()
