
Note: Unless specifies otherwise, each one of these source files is self contained.


(1) Programs:
=============
debil.c - The DEBiL model fitter*
debil_noEcc.c - The DEBiL model fitter, assuming circular orbits (faster)*
debil_noDoublePeriod.c - The DEBiL model fitter, without period doubling*
debil_noEcc_noDoublePeriod - The DEBiL model fitter, assuming circular orbits (faster), without period doubling*
periodS2.c - Finds the period of a given light curve
periodS2list.c - Finds the periods of a given list of light curves

* These programs all use the same source code, but with different switch settings


(2) DEBiL utilities:
====================
utils/splitPeriod.c - Assesses the EB light curve and creates appropriate alternate periods (double, half, or both).
utils/periodFilter.c - Filters out all but the strongly periodic light curves in a given list
utils/makeScatter.c - Adds statistical information to period database (obsolete)
utils/scatter2hist.c - Makes a histogram of the scatter information
utils/show2.c - An ASCII display of a given light curve (2-column format)
utils/show3.c - An ASCII display of a given light curve (3-column format)
utils/lc3sim.c   - Simulates a light curve with given parameters and noise
utils/addColor.c - Combines the Sumi database (with color info) with the DEBiL database
utils/makeEffTable.c - Precalculates a table for addProbability.c & addProbability_OGLE.c
utils/makeEffTable.res - A 1000x1000 table made from makeEffTable.c
utils/addProbability.c - Adds the (inverse) probability of seeing each binary eclipse
utils/debil2hist.c - Make histograms of the combined database
utils/filterDebil.c - filter DEBiL where a given column is above/below some threshold
utils/filterDebil2.c - filter DEBiL where a one column is above/below another column
utils/filterRoche.c - A more physical measure of detached systems. Filters if none, one or both lobes are filled
utils/filterSimilarDips.c - filters out light curves with similar dips, to remove double-true period cases 
utils/showDebil.c - makes *.fit and *.data for plotting the DEBiL fit (31-param format)
utils/plotAll.c - creates a SM script (plotAll.sm) for scatter plotting everything vs. everything
utils/histAll.c - creates a SM script (histAll.sm) for histogramming everything


(3) OGLE-versions:
======================
OGLE/periodS2.c - Finds the period of a given light curve (OGLE format)
OGLE/makeScatter.c - Adds statistical information to period database (obsolete)
OGLE/debil.c    - The DEBiL model fitter
OGLE/showDebil.c - makes *.fit and *.data for plotting the DEBiL fit (45-param format)
OGLE/show.c - Makes a series of ASCII displays for a list of light curves (OGLE format)
OGLE/show1.c - An ASCII display of a given light curve (OGLE format)
OGLE/showConf.c - Same as show.c, but skipping low confidence periods
OGLE/addProbability_OGLE.c - Adds the (inverse) probability of seeing each binary eclipse
OGLE/filterDebil.c - filter DEBiL (/w color) where a given column is above/below some threshold
OGLE/filterDebil2.c - filter DEBiL (/w color) where a one column is above/below another column
OGLE/filterDetached.c - filter DEBiL (/w color) where (R1+R2) is above/below some threshold
OGLE/filterRoche.c - A more physical measure of detached systems. Filters if none, one or both lobes are filled
OGLE/filterSimilarDips.c - filters out light curves with similar dips, to remove double-true period cases 
OGLE/fitCatalog.c - creates a SM script (fitCatalog.sm) displaying all the fits
OGLE/fitCatalogLarge.c - like fitCatalog.c only larger plots and thicker lines
OGLE/fitCatalogNoResidual.c - like fitCatalog.c only larger and without plots of the residuals
OGLE/makeOgle.c - Convert from "<time> <mag> <err>" to OGLE format
OGLE/makeOgleErr.c - Convert from "<time> <mag>" to OGLE format
