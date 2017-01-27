import os
import itertools
try:
    import cPickle as pickle
except:
    import pickle


def get_DEB_list():
    '''
    After looking at and flagging DEBs, get a list of them.
    '''

    with open('04_DATA_processed_fields.txt') as f:
        fieldstrs = f.readlines()

    fieldstrs = [f.strip('\n') for f in fieldstrs]

    debpaths = []
    weirdpaths = []
    allpaths = []
    for fieldstr in fieldstrs:
        cppath = '../data/CPs_cut/G'+fieldstr+'_20'
        pklpaths = [cppath+'/'+f for f in os.listdir(cppath) if f.endswith('.pkl')]
        allpaths.append(pklpaths)

        for pklpath in pklpaths:
            thislc = pickle.load(open(pklpath, 'rb'))
            if thislc['varinfo']['vartags'] == 'detached EB':
                debpaths.append(pklpath)
            elif thislc['varinfo']['vartags'] == 'weird variability':
                weirdpaths.append(pklpath)

        print('{:s} retrieved.'.format(fieldstr))

    allpaths = list(itertools.chain.from_iterable(allpaths))

    return debpaths, weirdpaths, allpaths


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

    debpaths, weirdpaths, allpaths = get_DEB_list()
    print_list_summaries(debpaths, weirdpaths, allpaths)

