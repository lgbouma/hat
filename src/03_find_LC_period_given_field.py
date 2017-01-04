'''
Usage (from within hat/src/):
> python 03_find_LC_period_given_field.py

Non-understood bugs:
    1) 'out' astropy table has slightly-different versions of the same EBs at
    slightly-different periods.
'''

import subprocess, os, pipes
import pandas as pd, numpy as np
from astropy.io import ascii
from astropy.table import vstack
from astrobase import hatlc

def exists_remote(host, path, is_dir=False):
    '''
    Test if a file exists at path on a host accessible with SSH.
    Example:
    > exists_remote('lbouma@phn12.astro.princeton.edu',
    >               '/H/BIGPROJ/hatuser/2007_hatnet_phot/',
    >               is_dir=True)
    kwargs:
    is_dir: boolean for if path being a directory. By default, assumes a file.
    '''
    test_flag = '-d' if is_dir else '-f'
    status = subprocess.call(
        ['ssh', host, 'test '+test_flag+' {}'.format(pipes.quote(path))])
    if status == 0:
        return True
    if status == 1:
        return False
    raise Exception('SSH failed')

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

def get_blsanalsum_output():
    '''
    This script copies BLS results from a run of HATpipe (hosted on the HAT 
    servers) locally for processing.
    '''
    host_name = 'lbouma@phn12.astro.princeton.edu'
    base_dir='/H/BIGPROJ/hatuser/2007_hatnet_phot/'
    #FIXME: generalize identification of neighboring fields
    nbhr_fields = ['199','154','198','155','197','153','115','292','242']
    nbhr_field_subdirs = list(map(lambda x:'G'+x+'/', nbhr_fields))
    nbhr_field_paths = [base_dir+sd for sd in nbhr_field_subdirs]
    write_path = '../data/HATpipe/blsanalsums/'

    for nbhr_field_path in nbhr_field_paths:
        write_str = write_path+nbhr_field_path.split('/')[-2]
        # First check if data is already exists locally.
        if os.path.exists(write_str+'_Cand.txt') or \
            os.path.exists(write_str+'_Old.txt'):
            print('{} is already local. Continue.'.format(write_str))
            continue

        # If no /Cand subdir (old photometry) try upper-dir:
        elif not exists_remote(host_name, nbhr_field_path+'Cand', is_dir=True):
            print('{} not found. Trying {}...'.format(nbhr_field_path+'Cand',
                nbhr_field_path+'blsanalsum.txt'))
            if exists_remote(host_name, nbhr_field_path+'blsanalsum.txt'):
                write_str += '_Old.txt'
                pull_str = host_name+':'+nbhr_field_path+'blsanalsum.txt '
                stdout, stderr = run_script('scp '+pull_str+write_str)
            else:
                print('{} not found. Continue.'.format(\
                        nbhr_field_path+'blsanalsum.txt'))
                continue

        # Try /Cand/blsanalsum.txt. If it exists, post-pend "_Cand.txt"
        elif exists_remote(host_name, nbhr_field_path+'Cand', is_dir=True):
            if exists_remote(host_name, nbhr_field_path+'Cand/blsanalsum.txt'):
                write_str += '_Cand.txt'
                pull_str = host_name+':'+nbhr_field_path+'Cand/blsanalsum.txt '
                # scp it over (via subprocess)
                stdout, stderr = run_script('scp '+pull_str+write_str)
                print('scp '+pull_str+write_str)

def parse_blsanalsum_for_BLS_peaks(max_DSP=50):
    '''
    Take blsanalsum.txt output, parse it into text file of:
    HATID, ROWID, PERIOD, Q, EPOCH, SNR, DSP, NTR, OOTSIG.
    '''
    max_DSP = str(max_DSP)
    data_path = '../data/HATpipe/blsanalsums/'
    write_path = '../data/HATpipe/blsanalsums/cuts/'
    existing_blsanalsums = os.listdir(data_path)
    blsanalsum_paths = [data_path+eb for eb in existing_blsanalsums if '.txt' \
            in eb]

    for blsanalsum_path in blsanalsum_paths:
        this_file = blsanalsum_path.split('/')[-1]
        pre_str = blsanalsum_path.split('/')[-1].split('.')[0]
        out_path = '{:s}{:s}_DSP_cut_{:s}.txt'.format(write_path, \
            pre_str, max_DSP)
        if os.path.exists(out_path):
            print('{:s} exists; continue.'.format(out_path))
        else:
            gawk_call='gawk \'$16 > {:s} {{print}}\' {:s} '.format(\
                max_DSP, blsanalsum_path)+\
                '| gawk \'{printf ("%15s\\t%d\\t%.7f\\t%.5f\\t%.7'+\
                'f\\t%.5e\\t%.5e\\t%d\\t%.4e\\n", $1, $2, $4, $5, $6,'+\
                ' $15, $16, $25, $30)}\''+' > {:s}'.format(out_path)

            run_script(gawk_call)
            print('Parsed {:s}.'.format(out_path))

def sqlitecurve_matching_and_LC_processing(DSP_lim, field='G199'):
    '''
    Given a specified field (e.g., G199), previous steps have created
    data/HATpipe/blsanalsums/cuts for that field and neighbors. Now match 
    with sqlitecurve file-IDs. Then create a LC file with format: time, 
    mag, err (space delimited).
    Pipe these as txt files to be processed w/ DEBiL.

    Take ~1 sec per LC (slow in astrobase's hatlc.read*).

    Params:
    DSP_lim: int. Limit on DSP used to find LCs w/ strong periodicity. In
        practice, used to ID blsanalsums/cuts to crossmatch w/ LCs.

    Optional TODO:
    * make a few phase-folded LCs (a la astrobase plotting tools) (Not
        quite necessary yet).
    * add / use 2MASS IDs in later LC creation
    '''

    cut_dpath = '../data/HATpipe/blsanalsums/cuts/'
    cut_paths = [cut_dpath+p for p in os.listdir(cut_dpath) if
            'cut_'+str(DSP_lim)+'.txt' in p]
    col_names = ('HATID', 'ROWID', 'PERIOD', 'Q', 'EPOCH', 'SNR', 'DSP', 'NTR',
            'OOTSIG')

    # Get astropy table of blsanalsum for HAT objects with DSP>DSP_lim
    for ix, cut_path in enumerate(cut_paths):
        if ix == 0:
            tab = ascii.read(cut_path, delimiter='\t', names=col_names)
        else:
            end = ascii.read(cut_path, delimiter='\t', names=col_names)
            tab = vstack([tab, end])

    # `tab` df contains all fields. `out` df has those in `field` (e.g., G199).
    tab = tab.to_pandas()
    tab = tab.set_index('HATID')
    field_substr = '-'+field.split('G')[-1]+'-'
    out = tab[tab.index.str.contains(field_substr)]

    # Cross-match with downloaded sqlitecurves. First, get intersection of out
    # ['HATID'] and sqlc_names.
    intersxn = [f[:15] for f in os.listdir('../data/LCs/'+field) if f[:15] in
            out.index]
    out['has_sqlc'] = out.index.map(lambda x: True if x in intersxn else False)
    np.testing.assert_array_equal(out['has_sqlc'], np.ones_like(out.index))

    # File name format: HAT-199-0025234-V0-DR0-hatlc.sqlite.gz
    LC_path = '../data/LCs/'+field+'/'
    tail_str = '-V0-DR0-hatlc.sqlite.gz'

    for hatid in out.index:
        # Write in DEBiL-readable format.
        write_path = '../data/DEBiL_LCs/'+field
        f_path = write_path+'/'+hatid+'.txt'
        if not os.path.isdir(write_path):
            os.makedirs(write_path)

        if not os.path.exists(f_path):
            # Get sqlitecurve data.
            obj_path = LC_path+hatid+tail_str
            lcd, msg = hatlc.read_and_filter_sqlitecurve(obj_path)
            # make sure all observations are at the same zero-point. required
            # for very long time-base HATNet LCs since these include different
            # instrument combinations over the decade of survey operations.
            normlcd = hatlc.normalize_lcdict(lcd)

            # FIXME: by default, take smallest EPD aperture, & only with "G" flag.
            time = normlcd['rjd'][normlcd['aiq_000']=='G']
            mag = normlcd['aep_000'][normlcd['aiq_000']=='G']
            err = normlcd['aie_000'][normlcd['aiq_000']=='G']

            f_id = open(f_path, 'wb+') # "+" creates if non-existent
            data = np.array([time, mag, err])
            np.savetxt(f_id, data.T, fmt=['%.7f','%.6f','%.6f'])
            f_id.close()

        else:
            continue

    # Write DEBiL "input list" of HAT-IDs and periods.
    write_path = '../data/DEBiL_heads/'+field+'_DSP'+str(DSP_lim)+'.txt'
    if not os.path.exists(write_path):
        f_id = open(write_path, 'wb+')
        data = np.array([out.index, out['PERIOD']])
        np.savetxt(f_id, data.T, fmt=['%15s', '%.6f'])
        f_id.close()


if __name__ == '__main__':
    DSP_lim = 30
    get_blsanalsum_output()
    parse_blsanalsum_for_BLS_peaks(max_DSP=DSP_lim)
    sqlitecurve_matching_and_LC_processing(DSP_lim, field='G199')
