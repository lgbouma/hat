'''
Given field, DSP_lim, and LCs, find first-pass at interesting EB candidates by
periodicity analysis.

Usage (from within hat/src/):
> python 03_find_LC_period_given_field.py field_id DSP_lim
Args:
    field_id: integer specifying HAT field ID number (HAT-???-1234567).

Non-understood bugs:
    1) 'out' astropy table has slightly-different versions of the same EBs at
    slightly-different periods.
'''

import subprocess, os, pipes
import pandas as pd, numpy as np
from astropy.io import ascii
from astropy.table import vstack
from astrobase import hatlc, periodbase, plotbase
from shutil import copyfile
import argparse

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

def get_blsanalsum_output(field_id=None):
    '''
    This script copies BLS results (blsanalsum.txt files) from a run of
    HATpipe, hosted on the HAT servers, locally for processing.
    '''
    assert type(field_id) == str, 'call get_blsanalsum_output with field_id'

    # Given neighboring fields, get paths on HATNet servers to download BLS
    # pipeline output
    host_name = 'lbouma@phn12.astro.princeton.edu'
    base_dir='/H/BIGPROJ/hatuser/2007_hatnet_phot/'

    # 2 lines of `#`-prepended comments, then lines of neighboring fields.
    # Each field is 3 digits, so narrow indexing below works.
    fo = open('DATA_nbhr_fields.txt', 'r')
    nbhr_lines = fo.readlines()
    fo.close()
    nbhr_fields_line = [l for l in nbhr_lines[2:] \
            if int(l[:3]) == int(field_id)]
    assert len(nbhr_fields_line) == 1, 'nbhr_fields_line should be singleton'
    nbhr_fields = nbhr_fields_line[0][4:-1].split(',')
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

        # Default: try /Cand/blsanalsum.txt. If it exists, post-pend "_Cand.txt"
        elif exists_remote(host_name, nbhr_field_path+'Cand', is_dir=True):
            if exists_remote(host_name, nbhr_field_path+'Cand/blsanalsum.txt'):
                write_str += '_Cand.txt'
                pull_str = host_name+':'+nbhr_field_path+'Cand/blsanalsum.txt '
                # scp it over (via subprocess)
                stdout, stderr = run_script('scp '+pull_str+write_str)
                print('scp '+pull_str+write_str)

        # If no /Cand subdir (old photometry) try upper-dir. If forced here,
        # post-pend "_Old.txt"
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


def parse_blsanalsum_for_BLS_peaks(max_DSP=None):
    '''
    Take blsanalsum.txt output, parse it into text file of:
    HATID, ROWID, PERIOD, Q, EPOCH, SNR, DSP, NTR, OOTSIG. Save it in
    /data/HATpipe/blsanalsums/cuts/
    '''
    assert type(max_DSP) == int, 'call parse_blsanalsum_for_BLS_peaks with max_DSP'

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

def sqlitecurve_matching_and_LC_processing(
        DSP_lim=None,
        field_id=None,
        DEBiL_write=False):
    '''
    Given a specified field (e.g., G199), previous steps have created
    data/HATpipe/blsanalsums/cuts for that field and neighbors. Now match
    with sqlitecurve file-IDs. Then create a LC file with format: time,
    mag, err (space delimited).
    Pipe these as txt files to /data/LCs_cut.

    Takes ~1 sec per LC (slow in astrobase's hatlc.read*).

    Optional TODO:
    * make a few phase-folded LCs (a la astrobase plotting tools) (Not
        quite necessary yet).
    * add / use 2MASS IDs in later LC creation
    '''
    assert type(field_id) == str

    field_name = 'G'+field_id
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
    field_substr = '-'+field_id+'-'
    out = tab[tab.index.str.contains(field_substr)]

    # Cross-match with downloaded sqlitecurves. First, get intersection of out
    # ['HATID'] and sqlc_names.
    intersxn = [f[:15] for f in os.listdir('../data/LCs/'+field_name)
            if f[:15] in out.index]
    out['has_sqlc'] = out.index.map(lambda x: True if x in intersxn else False)

    # File name format: HAT-199-0025234-V0-DR0-hatlc.sqlite.gz
    LC_read_path = '../data/LCs/'+field_name+'/' # where sqlitecurves already exist
    tail_str = '-V0-DR0-hatlc.sqlite.gz'
    # paths for LCs and EB checkplots
    LC_write_path = '../data/LCs_cut/'+field_name+'_'+str(DSP_lim)
    CP_write_path = '../data/CPs_cut/'+field_name+'_'+str(DSP_lim)
    for outpath in [LC_write_path, LC_write_path+'/periodcut',
            CP_write_path, CP_write_path+'/periodcut']:
        if not os.path.isdir(outpath):
            os.makedirs(outpath)

    for hatid in out.index:
        if np.all(out.ix[hatid]['has_sqlc']):
            LC_cut_path = LC_write_path+'/'+hatid+tail_str
            LC_periodcut_path = LC_write_path+'/periodcut/'+hatid+tail_str
            CP_cut_path = CP_write_path+'/'+hatid+'.png'
            CP_periodcut_path = CP_write_path+'/periodcut/'+hatid+'.png'

            if (not os.path.exists(CP_cut_path)) and \
            (not os.path.exists(CP_periodcut_path)):
                # Get sqlitecurve data.
                obj_path = LC_read_path+hatid+tail_str
                lcd, msg = hatlc.read_and_filter_sqlitecurve(obj_path)
                # Make sure all observations are at the same zero-point.
                normlcd = hatlc.normalize_lcdict(lcd)
                # Select recommended EPD aperture with 'G' flag. (Alternate 
                # approach: take the smallest aperture to minimize crowding).
                ap = next(iter(lcd['lcbestaperture']['ap']))
                times = normlcd['rjd'][normlcd['aiq_'+ap]=='G']
                mags = normlcd['aep_'+ap][normlcd['aiq_'+ap]=='G']
                errs = normlcd['aie_'+ap][normlcd['aiq_'+ap]=='G']

                # Period analysis: Stellingwerf phase dispersion minimization
                # and rerun Box-Least-Squares. Range of interesting periods:
                # 0.5days-25days.
                smallest_p, biggest_p = 0.5, 100.
                spdmp = periodbase.stellingwerf_pdm(times,mags,errs,
                    autofreq=True,
                    startp=smallest_p,
                    endp=biggest_p,
                    normalize=False,
                    stepsize=1.0e-4,
                    phasebinsize=0.05,
                    mindetperbin=9,
                    nbestpeaks=5,
                    periodepsilon=0.1, # 0.1days
                    sigclip=None, # no sigma clipping
                    nworkers=None)

                blsp = periodbase.bls_parallel_pfind(times,mags,errs,
                    startp=smallest_p,
                    endp=biggest_p, # don't search full timebase
                    stepsize=1.0e-5,
                    mintransitduration=0.01, # minimum transit length in phase
                    maxtransitduration=0.8,  # maximum transit length in phase
                    nphasebins=200,
                    autofreq=False, # figure out f0, nf, and df automatically
                    nbestpeaks=5,
                    periodepsilon=0.1, # 0.1
                    nworkers=None,
                    sigclip=None)

                # Make and save checkplot to be looked at.
                cp = plotbase.make_eb_checkplot(spdmp,blsp,times,mags,errs,
                       objectinfo=normlcd['objectinfo'],
                       findercmap='gray_r',
                       normto='globalmedian',
                       normmingap=4.0,
                       outfile=CP_cut_path,
                       sigclip=None,
                       varepoch='min',
                       phasewrap=True,
                       phasesort=True,
                       phasebin=0.002,
                       plotxlim=[-0.6,0.6])

                # Copy LCs with DSP>DSP_lim to /data/LCs_cut
                if not os.path.exists(LC_cut_path):
                    copyfile(obj_path, LC_cut_path)

                # Cut: if _all_ nbestperiods > 25 days, likely Cepheid or other
                # pulsation periodicity that is not of interest
                bestperiods = spdmp['nbestperiods']+blsp['nbestperiods']
                print(bestperiods)
                period_cut = 25. # days
                if np.all(np.array(bestperiods) > period_cut):
                    # Move the eb_checkplot, and the LC to subdirs
                    os.rename(LC_cut_path, LC_periodcut_path)
                    os.rename(CP_cut_path, CP_periodcut_path)


        else:
            continue

    if DEBiL_write:
        # Write DEBiL "input list" of HAT-IDs and periods.
        write_path = '../data/DEBiL_heads/'+field+'_DSP'+str(DSP_lim)+'.txt'
        if not os.path.exists(write_path):
            f_id = open(write_path, 'wb+')
            data = np.array([out.index, out['PERIOD']])
            np.savetxt(f_id, data.T, fmt=['%15s', '%.6f'])
            f_id.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('args', nargs='+')
    ao = parser.parse_args()
    assert len(ao.args) == 2, 'two args allowed: HAT field_id and DSP lim'

    # Params needed to identify LCs of interest for periodicity analysis
    field_id = ao.args[0] # keep as string
    DSP_lim = int(ao.args[1])

    # Periodicity analysis
    get_blsanalsum_output(field_id=field_id)
    parse_blsanalsum_for_BLS_peaks(max_DSP=DSP_lim)
    sqlitecurve_matching_and_LC_processing(DSP_lim=DSP_lim, field_id=field_id)
