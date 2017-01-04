def split_HAT_format_to_pandas_readable(lc_paths):
    '''
    We want to go from *rlc, *tfalc, *epdlc data (likely the last,
    since it has complete columns) to something python-readable that
    we can easily plot.
    This follows the format of data from https://hatsurveys.org/lightcurves,
    which I believe is valid after ~2013 (before that, columns change).
    '''
    for lc in lc_paths:
        print lc
        file_contents = []
        for line in open(lc, 'r'):
            file_contents.append(line)
        comments = [l for l in file_contents if ('#' in l) or ('\n' == l)]
        assert file_contents[len(comments)-1] == '\n', 'standard LC format has new line'+\
                ' before data dump'

        lc_data = file_contents[len(comments):]

        col_names = ['KEY', 'OFIELD', 'BJDc', 'MRAW0', 'MRAWERR0', 'PHOTFLAG0',
                     'MRAW1', 'MRAWERR1', 'PHOTFLAG1', 'MRAW2', 'MRAWERR2', 'PHOTFLAG2',
                     'MFIT0', 'MFIT1', 'MFIT2', 'MEPD0', 'MEPD1', 'MEPD2', 'MTFA0', 'MTFA1',
                     'MTFA2', 'X', 'Y', 'BG', 'BGERR', 'S', 'D', 'K', 'HA', 'Z', 'JDc']

        header_str = ''
        for ix, c in enumerate(col_names):
            if ix != len(col_names)-1:
                header_str += c+'\t'
            if ix == len(col_names)-1:
                header_str += c+'\n'

        data_fname = lc.split('/')[1]
        lc_data.insert(0, header_str)

        with open('pandas_readable/'+data_fname, 'w') as f:
            f.writelines(lc_data)
            f.close()

import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os

if __name__ == '__main__':
    '''
    Enter a directory name that contains your LC data, perhaps rsync'd from a
    phn??@astro.princeton.edu server.
    '''
    lc_dir = 'field_257/'
    lc_names = np.sort(os.listdir(lc_dir))
    # Assumes you want .tfalc
    lc_paths = [lc_dir+name for name in lc_names if '.tfalc' in name]

    split_HAT_format_to_pandas_readable(lc_paths)
