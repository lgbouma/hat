'''
Process visual inspection results, parsing them to csv files with path
information and basic facts from periodicity analysis.

Usage (from within hat/src/):
> python 04_summarize_visual_inspection.py remakecsvs

Args:
    remakecsvs: 0 or 1 for whether you want to regenerate csv files to be
    saved in ../data/visual_inspection_results/. This takes a while because it
    needs to read thousands of pickle files.
'''

import os
import itertools
import numpy as np
import pandas as pd
import datetime
import argparse
from shutil import copyfile
try:
    import cPickle as pickle
except:
    import pickle


def get_DEB_list():
    '''
    After looking at and flagging DEBs, get a list of them.
    Also get properties of the stars -- their RMS and SDSS r magnitudes.
    '''

    with open('04_DATA_processed_fields.txt') as f:
        fieldstrs = f.readlines()

    fieldstrs = [f.strip('\n') for f in fieldstrs]

    debpaths = []
    deblengths, debperiods, debrmss, debrmags = [], [], [], []
    weirdpaths = []
    allpaths = []
    for fieldstr in fieldstrs:
        #FIXME: be smarter about where checkplot data is stored?
        cppath = '/media/luke/LGB_tess_data/hat/CPs_cut/G'+fieldstr+'_20'
        #cppath = '../data/CPs_cut/G'+fieldstr+'_20'
        pklpaths = [cppath+'/'+f for f in os.listdir(cppath) if f.endswith('.pkl')]
        allpaths.append(pklpaths)

        for ix, pklpath in enumerate(pklpaths):
            if ix % 100 == 0:
                print('{:d} of {:d}'.format(ix, len(pklpaths)))

            thislc = pickle.load(open(pklpath, 'rb'))

            if thislc['varinfo']['vartags'] == 'detached EB':
                debpaths.append(pklpath)
                deblengths.append(len(thislc['magseries']['times']))
                debperiods.append(thislc['varinfo']['varperiod'])

                #follow Hartman et al. (2010) Eq 2 to calculate unbiased rms.
                meanmag = np.mean(thislc['magseries']['mags'])
                mags = thislc['magseries']['mags']
                N_p = 544
                N = thislc['objectinfo']['ndet']
                if N > N_p:
                    rms = np.sqrt(np.sum( (mags - meanmag)**2) / (N - N_p))
                else:
                    rms = -99
                debrmss.append(rms)

                #transform from 2MASS to Sloan r following Hartman+ 2010 Eq 1.
                if thislc['objectinfo']['sdssr']:
                    rmag = thislc['objectinfo']['sdssr']

                elif thislc['objectinfo']['jmag'] and \
                        thislc['objectinfo']['hmag'] and \
                        thislc['objectinfo']['kmag']:

                    rmag = 0.6975 + 2.9782 * thislc['objectinfo']['jmag'] \
                            - 0.8809 * thislc['objectinfo']['hmag'] \
                            - 1.1230 * thislc['objectinfo']['kmag']

                else:
                    rmag = -99

                debrmags.append(rmag)


            elif thislc['varinfo']['vartags'] == 'weird variability':
                weirdpaths.append(pklpath)


        print('{:s} retrieved.'.format(fieldstr))

    allpaths = list(itertools.chain.from_iterable(allpaths))

    return debpaths, weirdpaths, allpaths, deblengths, debperiods, \
            debrmss, debrmags



def print_list_summaries(debpaths, weirdpaths, allpaths):

    pathd = {'debs':debpaths, 'weirds':weirdpaths, 'all':allpaths}
    for k in pathd.keys():
        print('{:s}: length {:d}, fraction {:.3g}'.format(\
                k, len(pathd[k]), len(pathd[k])/len(pathd['all'])))


def _write_checkpoint_csvs(debpaths, weirdpaths, allpaths, deblengths,
        debperiods, debrmss, debrmags):
    '''
    Given results of `get_DEB_list`, turn them into saveable csv files to not
    have to call the slow "read through all of these pickle files" every time.

    This will save a "most up-to-date" version, and also a dated catalog with
    today's date (at most one per day).
    '''

    writepath = '../data/visual_inspection_results/'

    debhatids = [path.split('/')[-1][10:-4] for path in debpaths]
    weirdhatids = [path.split('/')[-1][10:-4] for path in weirdpaths]

    debinfo = pd.DataFrame(data={'period':debperiods,
             'rms_no_correction':debrmss,
             'rmag':debrmags,
             'length': deblengths,
             'pklpath': debpaths},
             index=debhatids)

    weirdinfo = pd.DataFrame(data={'pklpath': weirdpaths},
             index=weirdhatids)

    today = str(datetime.date.today())
    debinfo.to_csv(writepath+'deb_info.csv', index=True)
    debinfo.to_csv(writepath+'/archive/deb_info-'+today+'.csv', index=True)

    weirdinfo.to_csv(writepath+'weird_info.csv', index=True)
    weirdinfo.to_csv(writepath+'/archive/weird_info-'+today+'.csv', index=True)

    print('\nWrote detached EB info and weird object info to csv files...\n')


def _read_checkpoint_csvs(csvtype):
    '''
    Get pandas dataframe for checkpoint csv files. Detached EB dataframe comes
    with period, (non-corrected for EB variability) RMS, r mag, number of
    points in LC, and Checkplot paths.
    "Weird" dataframe comes just with Checkplot paths.

    Arg:
        csvtype (str): 'deb' or 'weird'
    '''

    assert csvtype == 'deb' or csvtype == 'weird', \
        'Only have done detached eclipsing binary and weird variability flags'

    if csvtype == 'deb':
        fname = 'deb_info.csv'
    elif csvtype == 'weird':
        fname = 'weird_info.csv'

    readpath = '../data/visual_inspection_results/'

    df = pd.read_csv(readpath+fname)
    df = df.set_index(['Unnamed: 0'])
    df.index.rename('hatid', inplace=True)

    return df


def write_pkl_dirs():
    '''
    Make ../data/debpkls and ../data/weirdpkls. These are directories of the
    pickle files of things flagged "weird" and "detached EB". They also serve
    as input to view all the LCs with the checkplotserver.
    '''

    debs = _read_checkpoint_csvs('deb')
    weirds = _read_checkpoint_csvs('weird')

    print('Writing checkplot pickles flagged as weird to ../data/weirdpkls/')
    for wp in weirds['pklpath']:
        copyfile(wp, '../data/weirdpkls/'+wp.split('/')[-1])

    print('Writing checkplot pickles flagged as DEB to ../data/debpkls/')
    for debp in debs['pklpath']:
        copyfile(debp, '../data/debpkls/'+debp.split('/')[-1])

    print('Done with write_pkl_dirs')


def main(remakecsvs):
    '''See docstring at top of program'''

    if remakecsvs:
        debpaths, weirdpaths, allpaths, deblengths, debperiods, \
                debrmss, debrmags = get_DEB_list()
        print_list_summaries(debpaths, weirdpaths, allpaths)

        _write_checkpoint_csvs(debpaths, weirdpaths, allpaths, deblengths,
                debperiods, debrmss, debrmags)

        write_pkl_dirs()


    else:
        print('Didn\'t make any new csvs. Why did you call this?')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('args', nargs='+')
    ao = parser.parse_args()
    assert len(ao.args) == 1, 'one arg allowed: whether you want to remake csvs'
    assert ao.args[0] == '0' or ao.args[0] == '1'
    remakecsvs = bool(int(ao.args[0]))

    main(remakecsvs)

