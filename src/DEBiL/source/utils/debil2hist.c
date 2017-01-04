/*********************************************************************
 *
 * This utility converts the scatter plots from DEBiL (with color) into histograms
 *
 * compile:  gcc debil2hist.c -o debil2hist -Wall -lm
 * 
 * run:   debil2hist ag.chi2.out 100 >! hist_chi2.txt
 *        debil2hist ag.chi3.out 100 > ! hist_chi3.txt
 *        debil2hist ag.chi4.out 100 > ! hist_chi4.txt
 *        debil2hist ag.chiAll.out 100 > ! hist_chiAll.txt
 *
 **********************************************************************/

// The minimal eccentricity needed for the argument of perihelion to be robust
#define MIN_OMEGA_ECC 0.1 

//#define MIN_NUMBER 1  // to prevent log(0) errors (must be a positive integer)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h> // MAXFLOAT

int main(const int argc, const char **argv)
{
  FILE *finScatter ;
  double period, ecc, r1, r2, b1, b2, sini, time0, omega, chi2, chiScore, sigmaMin2Amp ;
  double sigmaMax1Amp, sigmaMaxDiff, plateauWaviness, midScatterScore, meanDensity, maxDensity ;
  double maxChi2, maxSigmaMin2Amp, I1, I2, I12, color, invProbab ;
  int val, numHist ;
  unsigned long *histPeriod, *histEcc, *histR1,*histR2, *histSini, *histSiniZoom, *histTime0, *histOmega ;
  unsigned long *histChi2, *histChiScore, *histSigmaMin2Amp, *histSigmaMax1Amp, *histSigmaMaxDiff ;
  unsigned long *histPlateauWaviness, *histMidScatterScore, *histMeanDensity, *histMaxDensity ;
  unsigned long *histI1, *histI2, *histI12, *histColor, *histQ, *histDiffI, *histR1pR2 ;
  double *histPeriodCorrected ;

  if ((argc < 3) || (argc > 5))
    {
      printf ("Usage: %s <DEBiL /w color> <num hist> [max chi2] [max sigmaMin2Amp]\n", argv[0]) ;
      return (1) ;
    }
  finScatter = fopen (argv[1], "rt") ;
  if (!finScatter)
    {
      printf ("ERROR: couldn't open input file '%s'\n", argv[1]) ;
      return (2) ;
    }

  numHist = atoi (argv[2]) ;

  histPeriod = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histPeriodCorrected = (double*)calloc(numHist, sizeof(double)) ;
  histEcc = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histR1 = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histR2 = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histSini = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histSiniZoom = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histTime0 = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histOmega = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histChi2 = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histChiScore = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histSigmaMin2Amp = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histSigmaMax1Amp = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histSigmaMaxDiff = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histPlateauWaviness = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ; 
  histMidScatterScore = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histMeanDensity = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histMaxDensity = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ; 
  histI1 = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histI2 = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histI12 = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histColor = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histQ = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ;
  histDiffI = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ; 
  histR1pR2 = (unsigned long*)calloc(numHist, sizeof(unsigned long)) ; 

  if (!histPeriod || !histPeriodCorrected || !histEcc || !histR1 || !histR2 || !histSini || !histSiniZoom || !histTime0 ||
      !histOmega || !histChi2 || !histChiScore || !histSigmaMin2Amp || !histSigmaMax1Amp || !histSigmaMaxDiff ||
      !histPlateauWaviness || !histMidScatterScore || !histMeanDensity || !histMaxDensity || !histI1 || !histI2 ||
      !histI12 || !histColor || !histQ || !histDiffI || !histR1pR2)
    {
      printf ("ERROR: not enough memory\n") ;
      if (histPeriod) free (histPeriod) ;
      if (histPeriodCorrected) free (histPeriodCorrected) ;
      if (histEcc) free (histEcc) ;
      if (histR1) free (histR1) ;
      if (histR2) free (histR2) ;
      if (histSini) free (histSini) ;
      if (histSiniZoom) free (histSiniZoom) ;
      if (histTime0) free (histTime0) ;
      if (histOmega) free (histOmega) ;
      if (histChi2) free (histChi2) ;
      if (histChiScore) free (histChiScore) ;
      if (histSigmaMin2Amp) free (histSigmaMin2Amp) ;
      if (histSigmaMax1Amp) free (histSigmaMax1Amp) ;
      if (histSigmaMaxDiff) free (histSigmaMaxDiff) ;
      if (histPlateauWaviness) free (histPlateauWaviness) ;
      if (histMidScatterScore) free (histMidScatterScore) ;
      if (histMeanDensity) free (histMeanDensity) ;
      if (histMaxDensity) free (histMaxDensity) ;
      if (histI1) free (histI1) ;
      if (histI2) free (histI2) ;
      if (histI12) free (histI12) ;
      if (histColor) free (histColor) ;
      if (histQ) free (histQ) ;
      if (histDiffI) free (histDiffI) ;
      if (histR1pR2) free (histR1pR2) ;
      fclose (finScatter) ;
      return (3) ;
    }

  if (argc >= 4)
    maxChi2 = atof (argv[3]) ;
  else
    maxChi2 = MAXFLOAT ;

  if (argc == 5)
    maxSigmaMin2Amp = atof (argv[4]) ;
  else
    maxSigmaMin2Amp = MAXFLOAT ;


  while (23 == fscanf (finScatter, 
		       "%*d %*d %lf %lf %*f %lf %*f %lf %*f %lf %*f %lf %*f %lf %*f %lf %*f %lf %*f %*d %*d %lf %*f %*f %*f %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %*f %lf %*f %*f %*f %*f %*s %*s %lf\n", 
		       &period, &ecc, &r1, &r2, &b1, &b2, &sini, &time0,
		       &omega, &chi2, &chiScore, &sigmaMin2Amp, &sigmaMax1Amp, &sigmaMaxDiff,
		       &plateauWaviness, &midScatterScore, &meanDensity, &maxDensity,
		       &I1, &I2, &I12, &color, &invProbab))
    if ((chi2 <= maxChi2) && (sigmaMin2Amp <= maxSigmaMin2Amp)) 
      {
	val = (int)(log10(period * 10.0) * numHist / 4.0) ;      // transform [0.1, 1000] --> [0 , numHist]

	if ((val >= 0) && (val < numHist))
         {
	  histPeriod[val]++ ;
          histPeriodCorrected[val] += invProbab ;
         }

	//----

	val = (int)(ecc * numHist) ;                     // transform [0, 1] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histEcc[val]++ ;

	//----

	val = (int)(r1 * numHist) ;                      // transform [0, 1] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histR1[val]++ ;

	//----

	val = (int)(r2 * numHist) ;                     // transform [0.01, 1000] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histR2[val]++ ;

	//----

	val = (int)(((asin(sini) * 180.0 / M_PI) - 70.0) * numHist / 20.000001) ;  // transform [70, 90.000001] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histSini[val]++ ;

	//----

	val = (int)(((asin(sini) * 180.0 / M_PI) - 88.0) * numHist / 2.000001) ;  // transform [88, 90.000001] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histSiniZoom[val]++ ;

	//----
	val = (int)(time0 * numHist) ;                   // transform [0, 1] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histTime0[val]++ ;

	//----

	val = (int)(omega * numHist / 360.0) ;      // transform [0, 360] --> [0 , numHist]
	if ((val >= 0) && (val < numHist) && (ecc > MIN_OMEGA_ECC))
	  histOmega[val]++ ;

	//----

	val = (int)(chi2 * numHist / 10.0) ;    // transform [0, 10] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histChi2[val]++ ;

	//----

	val = (int)((chiScore - 0.2) * numHist) ;    // transform [0.2, 1.2] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histChiScore[val]++ ;

	//----

	val = (int)(sigmaMin2Amp * numHist / 20.0) ;  // transform [0, 20] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histSigmaMin2Amp[val]++ ;

	//----

	val = (int)((sigmaMax1Amp + 5.0) * numHist / 10.0) ;  // transform [-5, 5] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histSigmaMax1Amp[val]++ ;

	//----

	val = (int)(sigmaMaxDiff * numHist / 10.0) ;  // transform [0, 10] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histSigmaMaxDiff[val]++ ;

	//----

	val = (int)(plateauWaviness * numHist) ;  // transform [0, 1] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histPlateauWaviness[val]++ ;

	//----

	val = (int)(midScatterScore * numHist) ;   // transform [0, 1] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histMidScatterScore[val]++ ;

	//----

	val = (int)(log10(meanDensity * 1.0e6) * numHist / 9.0) ;      // transform [1e-6, 1000] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histMeanDensity[val]++ ;

	//----

	val = (int)(log10(maxDensity * 1.0e6) * numHist / 9.0) ;      // transform [1e-6, 1000] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histMaxDensity[val]++ ;

	//----

	val = (int)((I1 - 5.0) * numHist / 25.0) ;      // transform [5, 30] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histI1[val]++ ;

	//----

	val = (int)((I2 - 5.0) * numHist / 25.0) ;      // transform [5, 30] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histI2[val]++ ;

	//----

	val = (int)((I12 - 5.0) * numHist / 25.0) ;     // transform [5, 30] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histI12[val]++ ;

	//----

	val = (int)((color + 3.0) * numHist / 6.0) ;    // transform [-3, 3] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histColor[val]++ ;

	//----

	val = (int)((r2 / r1) * numHist) ;              // transform [0, 1] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histQ[val]++ ;

	//----

	val = (int)((3.0 + b2 - b1) * numHist / 10.0) ; // transform [-3, 7] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histDiffI[val]++ ;

	//----

	val = (int)((r1 + r2) * numHist) ;              // transform [0, 1] --> [0 , numHist]
	if ((val >= 0) && (val < numHist))
	  histR1pR2[val]++ ;
      }

#ifdef MIN_NUMBER
  for (val = 0 ; val < numHist ; val++)
    {
      if (histPeriod[val] == 0)  histPeriod[val] = MIN_NUMBER ;
      if (histPeriodCorrected[val] == 0.0)  histPeriod[val] = MIN_NUMBER ;
      if (histEcc[val] == 0)     histEcc[val] = MIN_NUMBER ;
      if (histR1[val] == 0)      histR1[val] = MIN_NUMBER ;
      if (histR2[val] == 0)      histR2[val] = MIN_NUMBER ;
      if (histSini[val] == 0)    histSini[val] = MIN_NUMBER ;
      if (histSiniZoom[val] == 0)         histSiniZoom[val] = MIN_NUMBER ;
      if (histTime0[val] == 0)   histTime0[val] = MIN_NUMBER ;
      if (histOmega[val] == 0)   histOmega[val] = MIN_NUMBER ;
      if (histChi2[val] == 0)    histChi2[val] = MIN_NUMBER ;
      if (histChiScore[val] == 0)         histChiScore[val] = MIN_NUMBER ;
      if (histSigmaMin2Amp[val] == 0)     histSigmaMin2Amp[val] = MIN_NUMBER ;
      if (histSigmaMax1Amp[val] == 0)     histSigmaMax1Amp[val] = MIN_NUMBER ;
      if (histSigmaMaxDiff[val] == 0)     histSigmaMaxDiff[val] = MIN_NUMBER ;
      if (histPlateauWaviness[val] == 0)  histPlateauWaviness[val] = MIN_NUMBER ;
      if (histMidScatterScore[val] == 0)  histMidScatterScore[val] = MIN_NUMBER ;
      if (histMeanDensity[val] == 0)      histMeanDensity[val] = MIN_NUMBER ;
      if (histMaxDensity[val] == 0)       histMaxDensity[val] = MIN_NUMBER ;
      if (histI1[val] == 0)      histI1[val] = MIN_NUMBER ;
      if (histI2[val] == 0)      histI2[val] = MIN_NUMBER ;
      if (histI12[val] == 0)     histI12[val] = MIN_NUMBER ;
      if (histColor[val] == 0)   histColor[val] = MIN_NUMBER ;
      if (histQ[val] == 0)       histQ[val] = MIN_NUMBER ;
      if (histDiffI[val] == 0)   histDiffI[val] = MIN_NUMBER ;
      if (histR1pR2[val] == 0)   histR1pR2[val] = MIN_NUMBER ;
    }
#endif
  
  for (val = 0 ; val < numHist ; val++)
    printf ("%f %lu %f  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu  %f %lu\n",
	    ((float)val * 4.0 / numHist) - 1.0, histPeriod[val], histPeriodCorrected[val],
	    ((float)val / numHist), histEcc[val],
	    ((float)val / numHist), histR1[val],
	    ((float)val / numHist), histR2[val],
	    ((float)val * 20.000001 / numHist) + 70.0, histSini[val],
	    ((float)val * 2.000001 / numHist) + 88.0, histSiniZoom[val],
	    ((float)val / numHist), histTime0[val],
	    ((float)val * 360.0 / numHist), histOmega[val],
	    ((float)val * 10.0 / numHist), histChi2[val],
	    ((float)val / numHist) + 0.2, histChiScore[val],
	    ((float)val * 20.0 / numHist), histSigmaMin2Amp[val],
	    ((float)val * 10.0 / numHist) - 5.0, histSigmaMax1Amp[val],
	    ((float)val * 10.0 / numHist), histSigmaMaxDiff[val],
	    ((float)val / numHist), histPlateauWaviness[val],
	    ((float)val / numHist), histMidScatterScore[val],
	    ((float)val * 9.0 / numHist) - 6.0, histMeanDensity[val],
	    ((float)val * 9.0 / numHist) - 6.0, histMaxDensity[val],
	    ((float)val * 25.0 / numHist) + 5.0, histI1[val],
	    ((float)val * 25.0 / numHist) + 5.0, histI2[val],
	    ((float)val * 25.0 / numHist) + 5.0, histI12[val],
	    ((float)val * 6.0 / numHist) - 3.0, histColor[val],
	    ((float)val / numHist), histQ[val],
	    ((float)val * 10.0 / numHist) - 3.0, histDiffI[val], 
	    ((float)val / numHist), histR1pR2[val]) ;


  free (histPeriod) ;
  free (histPeriodCorrected) ;
  free (histEcc) ;
  free (histR1) ;
  free (histR2) ;
  free (histSini) ;
  free (histSiniZoom) ;
  free (histTime0) ;
  free (histOmega) ;
  free (histChi2) ;
  free (histChiScore) ;
  free (histSigmaMin2Amp) ;
  free (histSigmaMax1Amp) ;
  free (histSigmaMaxDiff) ;
  free (histPlateauWaviness) ;
  free (histMidScatterScore) ;
  free (histMeanDensity) ;
  free (histMaxDensity) ;
  free (histI1) ;
  free (histI2) ;
  free (histI12) ;
  free (histColor) ;
  free (histQ) ;
  free (histDiffI) ;
  free (histR1pR2) ;
  fclose (finScatter) ;
  return (0) ;
}
