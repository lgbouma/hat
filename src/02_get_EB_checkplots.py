'''
Given field, DSP_lim, and LCs, find first-pass at interesting EB candidates by
periodicity analysis. Return EB checkplots.

Usage (from within hat/src/):
> python 02_get_EB_checkplots.py field_id DSP_lim > G???_DSP??.out

Args:
    field_id: integer specifying HAT field ID number (HAT-???-1234567).

Non-understood bugs:
    1) 'out' astropy table has slightly-different versions of the same EBs at
    slightly-different periods.
'''

import subprocess, os, pipes
import requests, os, time
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

def get_neighboring_fields(field_id):
    '''Given a field ID, get paths of neighboring fields, return a list of the
    fields, and a list of their paths on the HAT computers. '''

    # 2 lines of `#`-prepended comments, then lines of neighboring fields.
    # Each field is 3 digits, so narrow indexing below works.
    fo = open('DATA_nbhr_fields.txt', 'r')
    nbhr_lines = fo.readlines()
    fo.close()
    nbhr_fields_line = [l for l in nbhr_lines[2:]
            if int(l[:3]) == int(field_id)]
    assert len(nbhr_fields_line) == 1, 'nbhr_fields_line should be singleton'
    nbhr_fields = nbhr_fields_line[0][4:-1].split(',')
    nbhr_field_subdirs = list(map(lambda x:'G'+x+'/', nbhr_fields))

    base_dir='/H/BIGPROJ/hatuser/2007_hatnet_phot/'
    nbhr_field_paths = [base_dir+sd for sd in nbhr_field_subdirs]

    return nbhr_fields, nbhr_field_paths

def scp_blsanalsum_output(field_id=None):
    '''
    This script copies BLS results (blsanalsum.txt files) from a run of
    HATpipe, hosted on the HAT servers, locally for processing.
    Note: it copies the BLS results for your named field (`field_id`), as well
    as any neighboring fields.
    '''
    assert type(field_id) == str, 'call scp_blsanalsum_output with field_id'
    print('\nBeginning scp_blsanalsum_output...\n')

    # Given neighboring fields, get paths on HATNet servers to download BLS
    # pipeline output
    host_name = 'lbouma@phn12.astro.princeton.edu'

    nbhr_fields, nbhr_field_paths = get_neighboring_fields(field_id)
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
        else:
            print('{} not found. Trying {}...'.format(nbhr_field_path+'Cand',
                nbhr_field_path+'blsanalsum.txt'))
            if exists_remote(host_name, nbhr_field_path+'blsanalsum.txt'):
                write_str += '_Old.txt'
                pull_str = host_name+':'+nbhr_field_path+'blsanalsum.txt '
                stdout, stderr = run_script('scp '+pull_str+write_str)
            else:
                print('{} not found. Continuing.'.format(\
                        nbhr_field_path+'blsanalsum.txt'))
                continue


def parse_blsanalsum_for_BLS_peaks(max_DSP=None,Ntra_min=3):
    '''
    Take blsanalsum.txt output, parse it into text file of:
    HATID, ROWID, PERIOD, Q, EPOCH, SNR, DSP, NTR, OOTSIG. Save it in
    /data/HATpipe/blsanalsums/cuts/.
    Impose cuts on minimum number of transits (NTR) and their quality (NTV).

    Do this for every blsanalsum.txt results file found in
    "/data/HATpipe/blsanalsums/" (irrespective of field_id!).

    This simplifies subsequent analysis (and if you change Ntra_min, it will
    break the file naming convention).
    '''

    assert type(max_DSP) == int, 'call parse_blsanalsum_for_BLS_peaks with max_DSP'

    max_DSP = str(max_DSP)
    # Joel's "NTV" column example: 
    # "53200" means 2nd best transit has btwn 50-60% of covered pts as best 
    # transit. 3rd best has 30-40%. 4th 20-30%. Rest 0-10%.
    NTV_first_char_atleast = 6
    NTV_second_char_atleast = 4
    data_path = '../data/HATpipe/blsanalsums/'
    write_path = '../data/HATpipe/blsanalsums/cuts/'
    existing_blsanalsums = os.listdir(data_path)
    blsanalsum_paths = [data_path+eb for eb in existing_blsanalsums if '.txt' \
            in eb]

    print('\nBeginning parse_blsanalsum_for_BLS_peaks with:'+\
        '\nmax_DSP: {:s}'.format(max_DSP)+\
        '\nNtra_min: {:d}'.format(Ntra_min)+\
        '\nNTV_first_char_atleast: {:d}'.format(NTV_first_char_atleast)+\
        '\nNTV_second_char_atleast: {:d}\n'.format(NTV_second_char_atleast))

    for blsanalsum_path in blsanalsum_paths:
        this_file = blsanalsum_path.split('/')[-1]
        pre_str = blsanalsum_path.split('/')[-1].split('.')[0]
        out_path = '{:s}{:s}_DSP_cut_{:s}.txt'.format(write_path, \
            pre_str, max_DSP)
        if os.path.exists(out_path):
            print('{:s} exists; continue.'.format(out_path))
        else:
            # Apply gawk powers
            gawk_call='gawk \'($16 > {:s}) && ($25 >= {:d}) && '.format(\
                max_DSP, Ntra_min)+\
                '(substr($27,1,1) >= {:d}) && (substr($27,2,1) >= {:d})'.format(\
                NTV_first_char_atleast, NTV_second_char_atleast)+\
                ' {{print}}\' {:s} '.format(blsanalsum_path)+\
                '| gawk \'{printf ("%15s\\t%d\\t%.7f\\t%.5f\\t%.7'+\
                'f\\t%.5e\\t%.5e\\t%d\\t%.4e\\n", $1, $2, $4, $5, $6,'+\
                ' $15, $16, $25, $30)}\''+' > {:s}'.format(out_path)

            run_script(gawk_call)
            print('Parsed {:s}.'.format(out_path))


def download_parsed_LCs(DSP_lim=None,field_id=None):
    '''
    For LCs that met the cuts imposed in the blsanalsum parsing step (from
    Joel's big legacy BLS run), download their corresponding LCs from the HAT
    servers for subsequent periodicity analysis. The LC format is sqlitecurve.
    Note this is done for all stars in the named field (`field_id`) and its
    neighboring (on-sky-field), if the neighboring "Cand" dirs have stars that
    are in fact in the named field.

    Args:

    `field_id`: HAT-???-1234567

    Returns:
    `out`: astropy table with (ROWID, PERIOD, Q, EPOCH, SN, DSP, NTR, OOTSIG,
        and has_sqlc) columns.

    (Previously 02_retrive_LCs_given_sqlitecurve_head.py.)
    '''
    assert type(field_id)==str and type(DSP_lim)==int
    print('\nBeginning download_parsed_LCs...\n')

    field_name = 'G' + field_id # e.g., 'G081'
    head_base = '../data/sqlitecurve_headers/'
    head_path = head_base+field_name+'-sqlitecurves.csv'
    head_data = pd.read_csv(head_path)
    # Make the directory to write sqlite.gz light curves
    write_dir = '../data/LCs/'+field_name
    if not os.path.exists(write_dir):
        os.makedirs(write_dir)

    # Cross-reference lcfpaths with /data/HATpipe/blsanalsums/cuts HATIDs
    names = 'HATID,ROWID,PERIOD,Q,EPOCH,SNR,DSP,NTR,OOTSIG'
    col_names = tuple([n for n in names.split(',')])
    # Get neighboring fields, and paths to them on HAT computers.
    nbhr_fields, nbhr_field_paths = get_neighboring_fields(field_id)

    # For each blsanalsum in this field and neighboring fields, get astropy
    # table of blsanalsum data for HAT objects with DSP>DSP_lim, Ntra>Ntra_min,
    # and any other previous cuts.
    blsanalsum_path = '../data/HATpipe/blsanalsums/cuts/'
    blsdat_paths = []

    for field in nbhr_fields:
        thisbls_f = '{:s}_Cand_DSP_cut_{:d}.txt'.format(field_name, DSP_lim)
        blsdat_paths.append(blsanalsum_path+thisbls_f)

    for ix, blsdat_path in enumerate(blsdat_paths):
        if ix == 0:
            tab = ascii.read(blsdat_path, delimiter='\t', names=col_names)
        else:
            end = ascii.read(blsdat_path, delimiter='\t', names=col_names)
            tab = vstack([tab, end])

    # `tab` df contains all fields. `out` df has those in `field_id` (e.g., G199).
    # (The specific one being analyzed)
    tab = tab.to_pandas()
    tab = tab.set_index('HATID')
    field_substr = '-'+field_id+'-'
    out = tab[tab.index.str.contains(field_substr)]

    # Cross-match with downloaded sqlitecurves. First, get intersection of
    # `out` ['HATID'] and sqlc_names.
    intersxn = [hatid for hatid in out.index if hatid in
            np.array(head_data.hatid)]
    out['has_sqlc'] = out.index.map(lambda x: True if x in intersxn else False)

    # Download the sqlitecurves for our field that meet our cuts.
    url_base = 'http://data.hatsurveys.org/archive/'

    lcfpaths = ['hatnet/DR0/{:s}/'.format(field_id)+hatid+\
            '-V0-DR0-hatlc.sqlite.gz' for hatid in out.index \
            if np.all(out.ix[hatid].has_sqlc)]
    print('\nBeginning sqlitecurve download for {:s} ({:d} total).\n'.\
            format(field_name, len(np.unique(lcfpaths))))
    k = 0
    for ix, lcfpath in enumerate(np.unique(lcfpaths)):
        # If you already have the LC, don't redownload it. 
        if os.path.exists(write_dir+'/'+lcfpath.split('/')[-1]):
            print('{:d}: {:s} exists, continue.'.format(\
                    ix, write_dir+'/'+lcfpath.split('/')[-1]))
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

    print('\nDone downloading LCs that meet cuts.\n\n')

    return out




def periodicity_analysis(out,
        DSP_lim=None,
        field_id=None,
        DEBiL_write=False):
    '''
    Given a specified field (e.g., G199), previous steps have created
    data/HATpipe/blsanalsums/cuts for that field and neighbors (imposing cuts
    on DSP_lim, Ntra_min, and NTV). They've also downloaded the appropriate
    LCs.
    Now rerun the periodicity analysis for these LCs (Box-Least-Squares and
    Stellingwerf Phase Dispersion Minimization), and make eb_checkplots for
    subsequent looking-at ("visual inspection").

    Args:
       DEBiL_write (bool): whether to write a "name and best BLS period" file
           (in basically all use cases, not necessary).
    '''
    assert type(field_id) == str
    print('\nBeginning periodicity analysis...\n\n')

    # File name format: HAT-199-0025234-V0-DR0-hatlc.sqlite.gz
    field_name = 'G' + field_id # e.g., 'G081'
    LC_read_path = '../data/LCs/'+field_name+'/' # where sqlitecurves already exist
    tail_str = '-V0-DR0-hatlc.sqlite.gz'
    # paths for LCs and EB checkplots
    LC_write_path = '../data/LCs_cut/'+field_name+'_'+str(DSP_lim)
    CP_write_path = '../data/CPs_cut/'+field_name+'_'+str(DSP_lim)
    for outpath in [LC_write_path, LC_write_path+'/periodcut',
            CP_write_path, CP_write_path+'/periodcut']:
        if not os.path.isdir(outpath):
            os.makedirs(outpath)

    for ix, hatid in enumerate(out.index):
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

                # Cut: if _all_ nbestperiods > 30 days, likely Cepheid or other
                # pulsation periodicity that is not of interest.
                # While there may legitimately be detached binaries at > 30
                # days, we will struggle to find them, and their lower
                # harmonics should be in their periodograms.
                # Additional case to hide: if the only good peaks below 30 days
                # are multiples of 1 day (within 0.01days, absolute)
                bestperiods = spdmp['nbestperiods']+blsp['nbestperiods']
                best3periods = spdmp['nbestperiods'][:3]+\
                        blsp['nbestperiods'][:3]
                maxperiod = 30. # days
                minperiod = 0.502 # days; else this harmonic of 1d happens
                proxto1dmult = 0.01 # days
                bparr, b3parr = np.array(bestperiods), np.array(best3periods)
                if np.all(bparr > maxperiod) or \
                np.all(np.isclose(\
                (abs(b3parr[(b3parr<maxperiod)&(b3parr>minperiod)])-1.)%1.,\
                0., atol=1e-2)):
                    # Move the eb_checkplot, and the LC to subdirs
                    os.rename(LC_cut_path, LC_periodcut_path)
                    os.rename(CP_cut_path, CP_periodcut_path)


        else:
            print('{:d}: {:s} or LC counterpart exists; continue.'.\
                    format(ix, CP_cut_path))
            continue

    print('\nDone with periodicity analysis for {:s}.\n\n'.format(field_name))

    if DEBiL_write:
        # Write DEBiL "input list" of HAT-IDs and periods.
        write_path = '../data/DEBiL_heads/'+field+'_DSP'+str(DSP_lim)+'.txt'
        if not os.path.exists(write_path):
            f_id = open(write_path, 'wb+')
            data = np.array([out.index, out['PERIOD']])
            np.savetxt(f_id, data.T, fmt=['%15s', '%.6f'])
            f_id.close()


if __name__ == '__main__':
    '''See docstring at top of program'''
    parser = argparse.ArgumentParser()
    parser.add_argument('args', nargs='+')
    ao = parser.parse_args()
    assert len(ao.args) == 2, 'two args allowed: HAT field_id and DSP lim'

    # Params needed to identify LCs of interest for periodicity analysis
    Ntra_min = 3
    field_id = ao.args[0] # keep as string
    DSP_lim = int(ao.args[1])

    scp_blsanalsum_output(field_id=field_id)
    parse_blsanalsum_for_BLS_peaks(max_DSP=DSP_lim, Ntra_min=Ntra_min)
    out = download_parsed_LCs(DSP_lim=DSP_lim, field_id=field_id)
    periodicity_analysis(out, DSP_lim=DSP_lim, field_id=field_id)
