/*********************************************************************
 *
 * This utility converts the scatter plot for the light curve periods into histograms
 *
 * compile:  gcc scatter2hist.c -o scatter2hist -Wall -lm -O4
 * run:   scatter2hist ../periodFinder/results1/scatter.txt 100 15.0 >! hist1.txt
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h> // MAXFLOAT

int main(const int argc, const char **argv)
{
  FILE *finScatter ;
  double period, conf, period2, conf2, outAvr, outStdDiv, minConf, maxPeriod ;
  int numRid, val, numHist ;
  unsigned long *histPeriod, *histConf, *histPeriod2, *histConf2 ;
  unsigned long *histRid, *histAvr, *histStdDiv ;

  if ((argc < 3) || (argc > 5))
    {
      printf ("%s <scatter file> <num hist> [min confidence] [max period]\n", argv[0]) ;
      return (1) ;
    }

  numHist = atoi (argv[2]) ;

  histPeriod = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histConf = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histPeriod2 = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histConf2 = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histRid = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histAvr = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histStdDiv = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;

  if (!histPeriod || !histConf || !histPeriod2 || !histConf2 || !histRid || !histAvr || !histStdDiv)
    {
      if (histPeriod) free (histPeriod) ;
      if (histConf) free (histConf) ;
      if (histPeriod2) free (histPeriod2) ;
      if (histConf2) free (histConf2) ;
      if (histRid) free (histRid) ;
      if (histAvr) free (histAvr) ;
      if (histStdDiv) free (histStdDiv) ;
    }

  if (argc >= 4)
    minConf = atof (argv[3]) ;
  else
    minConf = 0.0 ;

  if (argc == 5)
    maxPeriod = atof (argv[4]) ;
  else
    maxPeriod = MAXFLOAT ;
  
  finScatter = fopen (argv[1], "rt") ;
  if (!finScatter)
    {
      printf ("ERROR: couldn't open input file '%s'\n", argv[1]) ;
      free (histPeriod) ;
      free (histConf) ;
      free (histPeriod2) ;
      free (histConf2) ;
      free (histRid) ;
      free (histAvr) ;
      free (histStdDiv) ;
      return (2) ;
    }

  while (fscanf (finScatter, "%*s %lf %lf %lf %lf %d %lf %lf\n", &period, &conf, &period2, &conf2, &numRid, &outAvr, &outStdDiv) == 7)
    if ((conf >= minConf) && (period <= maxPeriod)) 
      {
	val = (int)(log10(period * 100) * numHist * 0.2) ;      // transform [0.01, 1000] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histPeriod[val]++ ;

	//----
	val = (int)(conf * numHist / 20) ;                      // transform [0, 20] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histConf[val]++ ;

	//----

	val = (int)(log10(period2 * 100) * numHist * 0.2) ;      // transform [0.01, 1000] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histPeriod2[val]++ ;

	//----
	val = (int)(conf2 * numHist / 20) ;                      // transform [0, 20] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histConf2[val]++ ;

	//----
	val = numRid ;
	if ((val >= 0) && (val < numHist))
	  histRid[val]++ ;

	//----
	val = (int)(outAvr * numHist / 25) ;                    // transform [0, 25] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histAvr[val]++ ;

	//----
	val = (int)(log10(outStdDiv * 10000) * numHist * 0.2) ;   // transform [0.0001, 10] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histStdDiv[val]++ ;
      }

  for (val = 0 ; val < numHist ; val++)
    printf ("%f %lu  %f %lu  %f %lu  %f %lu  %d %lu  %f %lu  %f %lu\n",
	    pow(10, val / (numHist * 0.2)) / 100, histPeriod[val],
	    (float)val * 20.0 / numHist, histConf[val],
	    pow(10, val / (numHist * 0.2)) / 100, histPeriod2[val],
	    (float)val * 20.0 / numHist, histConf2[val],
	    val, histRid[val],
	    (float)val * 25.0 / numHist, histAvr[val],
	    pow(10, val / (numHist * 0.2)) / 10000, histStdDiv[val]) ;

  fclose (finScatter) ;
  free (histPeriod) ;
  free (histConf) ;
  free (histPeriod2) ;
  free (histConf2) ;
  free (histRid) ;
  free (histAvr) ;
  free (histStdDiv) ;
  return (0) ;
}
