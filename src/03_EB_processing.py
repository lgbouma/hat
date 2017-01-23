from astrobase import hatlc, checkplot, periodbase, plotbase
import numpy as np, matplotlib.pyplot as plt
import gzip
import pickle
import batman
from scipy.interpolate import interp1d

####################
# HELPER FUNCTIONS #
####################

def _load_sqlitecurve(hatid):
    '''
    For a given HAT LC that is known to be a detached eclipsing binary, load in
    its raw lightcurve (sqlitecurve format), filter it, and normalize it.
    '''

    detachedEBpath = '../data/detachedEBs/sqlcs/'
    ptail = '-V0-DR0-hatlc.sqlite.gz'
    lcd, msg = hatlc.read_and_filter_sqlitecurve(detachedEBpath+hatid+ptail)
    normlcd = hatlc.normalize_lcdict(lcd)
    # Select recommended EPD aperture with 'G' flag.
    # Get normalized unphased LC times, mags, errs
    ap = next(iter(lcd['lcbestaperture']['ap']))
    times = normlcd['rjd'][normlcd['aiq_'+ap]=='G']
    mags = normlcd['aep_'+ap][normlcd['aiq_'+ap]=='G']
    errs = normlcd['aie_'+ap][normlcd['aiq_'+ap]=='G']

    return normlcd, times, mags, errs


def _load_cp_picklefile(hatid):
    '''
    For the HAT LC corresponding to above sqlitecurve, read in the results from
    periodicty analysis of 02_get_EB_checkplots.py.
    '''
    #TODO: fix varperiod/varepoch to actually be read for general case

    ## Loading steps to get checkplot dictionary. N.b. this is _already normalized_ mags
    detachedEBpath = '../data/detachedEBs/checkplots/checkplot-'
    ptail = '.pkl.gz'

    cpd = checkplot._read_checkplot_picklefile(detachedEBpath+hatid+ptail)
    ## checkplot dictionary should have varperiod and varepoch
    cpd['varinfo']['varperiod'] = 2.101282
    cpd['varinfo']['varepoch'] = 53956.94080
    # Results from EB periodicity analysis:
    varperiod = cpd['varinfo']['varperiod']
    varepoch = cpd['varinfo']['varepoch']

    return cpd, varperiod, varepoch


def _initialize_batman_model(hatid, varepoch, varperiod, times, mags):
    '''
    Injects a periodic planet to passed magnitudes, returning injected
    magnitudes (a numpy array).
    '''
    #TODO: allow general planet parameters. Make more like a realistic CBP.
    #TODO: measure actual exposure times? This might matter...

    # Initialize batman.
    params = batman.TransitParams()
    params.t0 = varepoch - 1.3            #time of inferior conjunction: not a multiple of main transit
    params.per = varperiod*3.64159        #orbital period: ~3.2x that of EB (we want to find it)
    params.rp = 0.02**(1/2.)              #planet radius: a 2% dip (units: stellar radii)
    params.a = 9.                         #semi-major axis (in units of stellar radii)
    params.inc = 89.                      #orbital inclination (in degrees)
    params.ecc = 0.                       #eccentricity
    params.w = 90.                        #longitude of periastron (in degrees)
    params.u = [0.1, 0.3]                 #limb darkening coefficients
    params.limb_dark = "quadratic"        #limb darkening model

    exp_time = 0.00388                    #rough exposure time (days). ~about 5.5minutes

    m = batman.TransitModel(params,       #initialize model
                            times,
                            supersample_factor = 7,
                            exp_time = exp_time)

    relflux = m.light_curve(params)       #calculates light curve
    magdiff = -5/2. * np.log10(relflux)   #difference from median mag

    # Inject planet, make plots
    injmags = mags+magdiff

    plt.close('all')
    plt.scatter(times, magdiff)
    plt.savefig('../data/detachedEBs/plots/0_transitmodel-'+hatid+'.png', dpi=150)
    plt.close('all')
    plt.scatter(times, mags)
    plt.savefig('../data/detachedEBs/plots/1_noinj-'+hatid+'.png', dpi=150)
    plt.close('all')
    plt.scatter(times, injmags)
    plt.savefig('../data/detachedEBs/plots/2_injunphased-'+hatid+'.png', dpi=150)
    plt.close('all')

    return injmags


def _periodicity_analysis(hatid, times, mags, errs, normlcd,
        varperiod=None, isresidual=False):
    '''
    Given times, magnitudes, and errors, run Phase Dispersion Minimization and
    Box Least Squares, then write the phase-folded LCs, periodograms, and
    header data to a checkplot pickle.
    Typically, in 03_EB_processing.py, mags will be (data+injected planet), or
    (data + injected planet - EB model).

    Args:
        hatid (str): HAT-XXX-XXXXXXX
        times (np.array): times of measurement (typically RJD)
        mags (np.array): measured magnitudes at observing times
        errs (np.array): measured error on mags
        normlcd: lightcurve dictionary from early reading / astrobase parsing
    Keyword Args:
        varperiod (float, default None): EB period, used to set upper and lower
        bounds of frequency search.
        isresidual (bool, default False): if `mags` is from residuals, use
        this. It will do a denser frequency search, at >~2x the EB period.
    Returns:
        nothing
    '''
    #TODO: fix args, add docs

    if not varperiod:
        #injection has happened
        smallest_p = 0.5
        cpwritepre = '../data/detachedEBs/injs/checkplot-inj-'
    elif isinstance(varperiod,float) and isresidual:
        smallest_p = varperiod*2. + 0.1
        cpwritepre = '../data/detachedEBs/injres/checkplot-injresiduals-'
    else:
        raise ValueError('smallest_p never properly assigned')

    biggest_p = min((times[-1] - times[0])/2.01, 100.)

    print('\nStellingwerf...\n')
    spdmp = periodbase.stellingwerf_pdm(times,mags,errs,
                        autofreq=True,
                        startp=smallest_p,
                        endp=biggest_p,
                        normalize=False,
                        stepsize=1.0e-5,
                        phasebinsize=0.03,
                        mindetperbin=9,
                        nbestpeaks=5,
                        periodepsilon=0.1, # 0.1days
                        sigclip=None, # no sigma clipping
                        nworkers=None)

    varinfo = {'objectisvar':True,
              'vartags':None,
              'varisperiodic':True,
              'varperiod':spdmp['bestperiod'],
              'varepoch':None}

    print('\nBLS...\n')
    blsp = periodbase.bls_parallel_pfind(times, mags, errs,
                                        startp=smallest_p,
                                        endp=biggest_p,
                                        stepsize=1.0e-5,
                                        mintransitduration=0.01,
                                        maxtransitduration=0.6,
                                        nphasebins=200,
                                        autofreq=False,
                                        nbestpeaks=5,
                                        periodepsilon=0.1,
                                        nworkers=None,
                                        sigclip=None)

    lspinfolist = [spdmp, blsp]
    cp = checkplot.checkplot_pickle(lspinfolist, times, mags, errs,
                                   nperiodstouse=3,
                                   objectinfo=normlcd['objectinfo'],
                                   varinfo=varinfo,
                                   findercmap='gray_r',
                                   normto='globalmedian',
                                   normmingap=4.,
                                   outfile=cpwritepre+hatid+'.pkl.gz',
                                   sigclip=[10.,-3.],
                                   varepoch='min',
                                   phasewrap=True,
                                   phasesort=True,
                                   phasebin=0.002,
                                   plotxlim=[-0.6,0.6],
                                   plotdpi=150,
                                   returndict=False,
                                   pickleprotocol=3,
                                   greenhighlight=False,
                                   xgridlines=[-0.5,0.,0.5])

    #since we don't have easy checkplot_pickle to png written yet, use twolsp_checkplot_png
    #TODO: remove any png writing once happy with output
    if not varperiod:
        injpath = '../data/detachedEBs/plots/3_checkplot-inj-'
    elif isinstance(varperiod,float) and isresidual:
        injpath = '../data/detachedEBs/plots/5_checkplot-injresiduals-'
    else:
        raise ValueError('injpath never properly assigned')

    checkplot.twolsp_checkplot_png(spdmp, blsp, times, mags, errs,
                                   objectinfo=normlcd['objectinfo'],
                                   findercmap='gray_r',
                                   normto='globalmedian',
                                   normmingap=4.,
                                   outfile=injpath+hatid+'.png',
                                   sigclip=[10.,-3.],
                                   varepoch='min',
                                   phasewrap=True,
                                   phasesort=True,
                                   phasebin=0.002,
                                   plotxlim=[-0.6,0.6],
                                   plotdpi=150)
    plt.close('all')



def _subtract_median_EB_signal(hatid, varperiod, varepoch):
    '''
    Once injected planet periodograms have been run, get the median EB
    signal and subtract it to get residuals.

    Args:
        hatid (str): HAT-XXX-XXXXXXX
        varperiod (float): best EB period from PDM in periodicity analysis
        varepoch (float): epoch of folded LC
    Returns:
        residualmags (np array): data + injection - median model
    '''
    #TODO Propagate errors properly.

    #confirm EB period.
    injEBpath = '../data/detachedEBs/injs/checkplot-inj-'
    ptail = '.pkl.gz'
    cpd = checkplot._read_checkplot_picklefile(injEBpath+hatid+ptail)

    #unphased LC times, mags, errs
    times = cpd['magseries']['times']
    mags = cpd['magseries']['mags']
    errs = cpd['magseries']['errs']
    #phased data, then binned data (blue points)
    phase, phasedmags = cpd['pdm'][0]['phase'], cpd['pdm'][0]['phasedmags']
    binphase, binphasedmags = cpd['pdm'][0]['binphase'], cpd['pdm'][0]['binphasedmags']

    twophasetimes = varepoch + (binphase*varperiod) #times over two phases (stored in binphase originally)

    #now extend twophasetimes to beginning and end of baseline
    modeltimes = np.array(twophasetimes)
    modelmags = binphasedmags

    exitf,upf,downf = False,True,True
    diff = modeltimes[-1] - modeltimes[0]

    while not exitf:

        while upf:
            onetimeup = np.array(modeltimes[-1] - modeltimes[0] + twophasetimes)[1:]
            onemagup = np.array(binphasedmags[1:])
            #oneerrup = np.array() #obnoxiously, these aren't kept thru the phase curve-making. should be needed to propagate uncertainties

            modeltimes = np.append(modeltimes, onetimeup)
            modelmags = np.append(modelmags, onemagup)
            if max(modeltimes) > max(times):
                upf = False

        k = 1
        while downf:
                #careful of lengths
                onetimedown = np.array(-k*diff + twophasetimes)[:-1]
                onemagdown = np.array(binphasedmags[:-1])
                #oneerrup = np.array() #not kept thru phase curve making 

                modeltimes = np.insert(modeltimes, 0, onetimedown)
                modelmags = np.insert(modelmags, 0, onemagdown)

                k+=1

                if min(modeltimes) < min(times):
                    downf = False
                    exitf = True

    f = interp1d(modeltimes, modelmags, kind='linear')

    interpmodelmags = f(times)

    residualmags = mags - interpmodelmags # n.b. these mags are the injected ones

    #save figure showing model and residuals
    plt.close('all')
    f, axs = plt.subplots(3, sharex=True, figsize=(15,8*1.5))
    axs[0].scatter(times, mags, lw=0, c='black', s=10)
    axs[0].set(ylabel='DEB data (w/ planet injected)')
    axs[1].scatter(times, interpmodelmags, lw=0, c='black', s=10)
    axs[1].set(ylabel='median model\n(medians from binned phase-folded DEB LC)')
    axs[2].scatter(times, residualmags, lw=0, c='black', s=10)
    axs[2].set(ylabel='data - model')
    f.tight_layout()
    f.savefig('../data/detachedEBs/plots/4_median_EB_modelLC_and_residuals.png', dpi=150)
    plt.close('all')
    print('\nMade figure of median EB model and residuals...\n')

    return times, residualmags, errs


def main():
    #TODO: fix all save paths to somewhere in data (for figures being written)
    #TODO: propagate model uncertainty to residual magnitudes

    hatid = 'HAT-199-0000452'

    normlcd, times, mags, errs = _load_sqlitecurve(hatid)

    cpd, varperiod, varepoch = _load_cp_picklefile(hatid)

    injmags = _initialize_batman_model(hatid, varepoch, varperiod, times, mags)

    _periodicity_analysis(hatid, times, injmags, errs, normlcd)

    times, residualmags, errs = _subtract_median_EB_signal(hatid, varperiod, varepoch)

    _periodicity_analysis(hatid, times, residualmags, errs, normlcd,
            varperiod=varperiod, isresidual=True)

    print('done.')

if __name__ == '__main__':
    main()
