import os
import itertools
import numpy as np
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
        #FIXME: be smarter about where checkplot data is stored
        cppath = '/media/luke/LGB_tess_data/hat/CPs_cut/G'+fieldstr+'_20'
        #cppath = '../data/CPs_cut/G'+fieldstr+'_20'
        pklpaths = [cppath+'/'+f for f in os.listdir(cppath) if f.endswith('.pkl')]
        allpaths.append(pklpaths)

        for ix, pklpath in enumerate(pklpaths):
            if ix % 100 == 0:
                print('{:d} of {:d}', ix, len(pklpaths))

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

def get_DEB_ids(debpaths):
    '''
    Given debpaths, parse it into a list of dictionaries where each dictionary
    gets: hatid (str), pickle path (str), original LC path (for sigma clipping
    reversion).
    '''



if __name__ == '__main__':

    debpaths, weirdpaths, allpaths, deblengths, debperiods, \
            debrmss, debrmags = get_DEB_list()
    print_list_summaries(debpaths, weirdpaths, allpaths)

