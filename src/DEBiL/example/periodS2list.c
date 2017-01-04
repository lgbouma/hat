/**************************************************************************
 * this program finds the period of a given 3-column light curve table
 * using the minimum variance with a linear best fit, in each bin of the 
 * phase histogram. Uses the Schwarzenberg-Czerny (1989) algorithm.
 * Provides the second best period as well.
 *
 * compile:  gcc periodS2list.c -o periodS2list -Wall -O4 -lm
 * run:      periodS2 /data/clusters/OGLE/Bulge/BUL_SC4/bul_sc4_7143.dat out.txt err.txt --> (p=0.623) 
 *           periodS2 /data/clusters/OGLE/Bulge/BUL_SC11/bul_sc11_1177.dat out.txt err.txt --> problematic (p=0.9629 or p=26.9)
 *
 **************************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>  // memset
#include <stdlib.h>


#define HIST_SIZE 8
// identify outliers: must be larger than 1.0 or else, in theory, all the points could be outliers
#define STDDIV_LIMIT 5.0
#define MIN_PERIOD 0.1         // starting point for scan
// resolution of scan: must be smaller than 1.0 or else periods will be completely missed
#define SUBSAMPLE 0.05  // resolution of "fine tune" scan
#define FINE_TUNE 0.0005
#define MAX_DOUBLE_CHECK_MULTIPLE 19
#define MAX_PERIOD_DIFF_MULTIPLE 9

#define MIN_NUM_IN_BIN 3
#define EPSILON 10E-10

#define ERROR_SCORE 1000.0  // assumed to be bigger than any valid score (>> 1)

#define OGLE2_PERIOD_FILTER
//#define DEBUG

unsigned long  histN[HIST_SIZE] ;
double histA[HIST_SIZE], histB[HIST_SIZE], histC[HIST_SIZE], histD[HIST_SIZE], histE[HIST_SIZE] ;

//inline
double sqr (double x)
{
  return (x * x) ;
}

// Returns the modulo 1 value of x
double mod1 (double x)
{
  return (x - floor(x)) ;
}

//-------------------------------


// returns the reduced size of the array
unsigned int normalize (unsigned int size, double *mag, double *time,
			double *outAvr, double *outStdDiv)
{
  unsigned int i, doPrune ;
  double Avr, StdDiv ;

  // step 1: calc the average and standard deviation
  Avr = mag[0] ;
  for (i = 1 ; i < size ; i++)
    Avr += mag[i] ;
  Avr /= size ;

  StdDiv = sqr(mag[0] - Avr) ;
  for (i = 1 ; i < size ; i++)
    StdDiv += sqr(mag[i] - Avr) ;
  StdDiv = sqrt(StdDiv / size) ;


  // step 2: prune outlier
#ifdef STDDIV_LIMIT
  do
    {
      doPrune = 0 ;
      i = 0 ;

      while (i < size)
	{
	  if (fabs(mag[i] - Avr) > (StdDiv * STDDIV_LIMIT))
	    {
#ifdef DEBUG
	      printf ("outlier: %f (sigma = %f %f) \n", mag[i], fabs(mag[i] - Avr) / StdDiv, StdDiv) ;
#endif
	      size-- ;
	      time[i] = time[size] ;
	      mag[i] = mag[size] ;
	      doPrune = 1 ;
	    }
	  else i++ ;
	}
#endif

      Avr = mag[0] ;
      for (i = 1 ; i < size ; i++)
	Avr += mag[i] ;
      Avr /= size ;

      StdDiv = sqr(mag[0] - Avr) ;
      for (i = 1 ; i < size ; i++)
	StdDiv += sqr(mag[i] - Avr) ;
      StdDiv = sqrt(StdDiv / size) ;
    }
  while (doPrune) ;


  // step 3: normalize to [-1,1]
  for (i = 0 ; i < size ; i++)
    mag[i] = (mag[i] - Avr) / StdDiv ;

  *outAvr = Avr ;
  *outStdDiv = StdDiv ;
  return (size) ;
}



double testPeriod (unsigned long size, double *time, double *mag, double period)
{
  int i, index ;
  double sum = 0.0, X, Y, tmp, s1 = 0.0, L2 = 0.0 ;
  unsigned long N ;

#ifdef OGLE2_PERIOD_FILTER
  //---------------- filtering out problematic periods ----------------------
  // filtering out bands around 1 siderial day and its rational multiples (harmonics)
  // these values are specificly tuned for the window function of OGLE II bulge

  if (fabs(period - 0.999 ) < 0.01  ) return (ERROR_SCORE) ; //   1
  if (fabs(period - 1.997 ) < 0.012 ) return (ERROR_SCORE) ; //   2
  if (fabs(period - 2.995 ) < 0.016 ) return (ERROR_SCORE) ; //   3
  if (fabs(period - 3.993 ) < 0.02  ) return (ERROR_SCORE) ; //   4
  if (fabs(period - 4.989 ) < 0.02  ) return (ERROR_SCORE) ; //   5
  if (fabs(period - 5.985 ) < 0.015 ) return (ERROR_SCORE) ; //   6
  if (fabs(period - 0.4991) < 0.002 ) return (ERROR_SCORE) ; //  1/2
  if (fabs(period - 1.4965) < 0.004 ) return (ERROR_SCORE) ; //  3/2
  if (fabs(period - 2.493 ) < 0.003 ) return (ERROR_SCORE) ; //  5/2
  if (fabs(period - 0.3327) < 0.0008) return (ERROR_SCORE) ; //  1/3
  if (fabs(period - 0.6651) < 0.0008) return (ERROR_SCORE) ; //  2/3
  if (fabs(period - 1.3299) < 0.0012) return (ERROR_SCORE) ; //  4/3
  if (fabs(period - 0.2494) < 0.0005) return (ERROR_SCORE) ; //  1/4
#endif

  //----------------- initialize ----------------------------
  
  memset (histN, 0, HIST_SIZE * sizeof (unsigned long)) ;
  memset (histA, 0, HIST_SIZE * sizeof (double)) ;
  memset (histB, 0, HIST_SIZE * sizeof (double)) ;
  memset (histC, 0, HIST_SIZE * sizeof (double)) ;
  memset (histD, 0, HIST_SIZE * sizeof (double)) ;
  memset (histE, 0, HIST_SIZE * sizeof (double)) ;
  
  for (i = 0 ; i < size ; i++)
    {
      X = mod1(time[i] / period) ;
      index = (unsigned int)(HIST_SIZE * X) ;
      Y = mag[i] ;
      
      histN[index]++ ;
      histA[index] += X ;
      histB[index] += X * X ;
      histC[index] += X * Y ;
      histD[index] += Y ;
      histE[index] += Y * Y ;
    }
  
  //--------------------- scatter --------------------------------------
  
  for (i = 0 ; i < HIST_SIZE ; i++)
    {
      N = histN[i] ;
      if (N < MIN_NUM_IN_BIN) return (ERROR_SCORE) ;
      Y = histD[i] ;
      X = Y * Y / N ;

      tmp = sqr(histA[i]) - (histB[i] * N) ;
      if (fabs(tmp) < EPSILON) return (ERROR_SCORE) ;  // will happen if all bin points are vertical	
      sum += Y ;
      s1 += X ;
      L2 += (sqr((histC[i] * N) - (Y * histA[i])) / (N * tmp)) + histE[i] - X ;
    }

  s1 -= sum * sum / size ;  // no need to divide s1 by a constant  (HIST_SIZE-1)
  L2 /= size - HIST_SIZE ;

  if ((s1 <= 0.0) || (L2 <= 0.0)) return (ERROR_SCORE) ;
  return (log(L2 / s1)) ;     // -ln(theta_AOV)
}



void printError (const char *errFilename, const char *filename, const char *errStr)
{
  FILE *ferr = fopen (errFilename, "at") ;

  if (!ferr)
    {
      printf ("ERROR: couldn't open the output error file ('%s')\n", errFilename) ;
      return ;
    }

  fprintf (ferr, "%s: %s\n", filename, errStr) ;
  fclose (ferr) ;
} 



int main(int argc, char **argv)
{
  unsigned long size, prevSize, numRid ;
  unsigned int a, b, aBest, bBest, isMultiple ;
  unsigned long N ;
  double T, dT, period, score, bestScore, bestPeriod ;
  double *time, *mag ;
  double sum, sumSqr, stddiv, tmpTime, tmpMag, tmpErr, outAvr, outStdDiv ;
  char str[32], filename[256] ;
  FILE *finList, *fin, *fout ;


  if (argc != 4)
    {
      printf ("usage: %s <input LC list> <output results filename> <output error filename>\n", argv[0]) ;
      return (1) ;
    }

  if (!strcmp(argv[1], argv[2]) ||
      !strcmp(argv[1], argv[3]) ||
      !strcmp(argv[2], argv[3]))
    {
      printf ("ERROR: all the input/output filenames must be different\n") ;
      return (2) ;
    }

  if (!(finList = fopen (argv[1], "rt")))
    {
      printError (argv[3], argv[1], "Can't open LC list") ;
      return (3) ;
    }

  if (!(fout = fopen (argv[2], "at")))
    {
      printError (argv[3], argv[2], "Can't open output file") ;
      return (4) ;
    }

  //-----------------

  while (1 == fscanf(finList, "%s\n", filename))
    {
      printf ("%s\n", filename) ;
   
      if (!(fin = fopen (filename, "rt")))
	{
	  printError (argv[3], filename, "Can't open input file") ;
	  continue ;
	}

      size = 0 ;
      while (1 == fscanf(fin, "%*f %*f %lf\n", &tmpErr))
	if (tmpErr > 0)  // sign of invalid magnitude
	  size++ ;

      if (size < (HIST_SIZE * MIN_NUM_IN_BIN))
	{
	  sprintf (str, "Not enough data points (%lu)\n", size) ;
	  printError (argv[3], filename, str) ;
	  fclose (fin) ;
	  continue ;
	}

      time = (double*)malloc(size * sizeof(double)) ;
      mag = (double*)malloc(size * sizeof(double)) ;

      if (!time || !mag)
	{
	  printError (argv[3], filename, "Not enough memory") ;  
	  if (time) free (time) ;
	  if (mag) free (mag) ;
	  fclose (fin) ;
	  continue ;
	}


      rewind(fin) ;
      prevSize = size ;
      size = 0 ;
      while (3 == fscanf(fin, "%lf %lf %lf\n", &tmpTime, &tmpMag, &tmpErr))
	if ((tmpErr > 0) && (size < prevSize))  // sign of invalid magnitude
	  {
	    time[size] = tmpTime ;
	    mag[size] = tmpMag ;
	    size++ ;
	  }

      if (prevSize != size)
	{
	  sprintf (str, "Unequal read sizes (%lu  %lu)\n", prevSize, size) ;
	  printError (argv[3], filename, str) ;
	  free (mag) ;
	  free (time) ;
	  fclose (fin) ;		 
	  continue ;
	}

      T = time[size-1] - time[0] ;
      numRid = size ;
      size = normalize (size, mag, time, &outAvr, &outStdDiv) ;  // jumbles up the order
      numRid -=size ;

      N = 0 ;
      sum = 0.0 ;
      sumSqr = 0.0 ;
      bestScore = ERROR_SCORE ;
      bestPeriod = 0.0 ;

      // ****** full scan ********
      for (period = MIN_PERIOD ; period <= T ; period += (SUBSAMPLE * period * period / T))
	{
	  score = testPeriod(size, time, mag, period) ;

	  if (score != ERROR_SCORE)
	    {
	      N++ ;
	      sum += score ;
	      sumSqr += (score * score) ;
	    }

	  if (score < bestScore)
	    {
	      bestScore = score ;
	      bestPeriod = period ;

#ifdef DEBUG
	      printf ("1: %f %f\n", period, score) ;
#endif
	    }
	}


      if (bestScore == ERROR_SCORE)
	{
	  printError (argv[3], filename, "No valid periods") ;
	  fclose (fin) ;
	  free (mag) ;
	  free (time) ;
	  continue ;
	}

      // ****** fine tune ********
      // note that in principal, the fine-tune could make the period deviate from [MIN_PERIOD, T]
      dT = bestPeriod * bestPeriod / T ;  // period step - have a lot of overlap on purpose
      for (period = bestPeriod - dT ; (period < bestPeriod + dT) && (period <= T) ; period += (FINE_TUNE * dT))
	{
	  score = testPeriod (size, time, mag, period) ;
	  
	  if (score < bestScore)
	    {
	      bestScore = score ;
	      bestPeriod = period ;
	    }
	}

      aBest = 1 ;
      bBest = 1 ;
      isMultiple = 0 ;

      // ****** double-check period multiples ******
      // (should be very rare, but costs almost nothing)
      for (a = 1 ; a <= MAX_DOUBLE_CHECK_MULTIPLE ; a++)
	for (b = 1 ; b <= MAX_DOUBLE_CHECK_MULTIPLE ; b++)
	  if (a != b)
	    {
	      period = (bestPeriod * a) / b ;

	      if ((period > MIN_PERIOD) && (period <= T))
		{
		  score = testPeriod (size, time, mag, period) ;

		  if (score < bestScore)
		    {
		      isMultiple = 1 ;
		      aBest = a ;
		      bBest = b ;
		      bestScore = score ;
		    }
		}
	    }

      sum /= N ;    // average
      sumSqr /= N ; // average square
      stddiv = sqrt(sumSqr - sqr(sum)) ;

#ifdef DEBUG
      printf ("avr= %f   stddiv= %f\n", sum, stddiv) ;
#endif

      if (isMultiple)
	{
	  bestPeriod = (bestPeriod * aBest) / bBest ;
	  sprintf (str, "found a better multiple (%u/%u)", aBest, bBest) ;
	  printError (argv[3], filename, str) ;

#ifdef DEBUG
	  printf ("Warning: found a better multiple (%u/%u)\n", aBest, bBest) ;
#endif
	}


      //----------------------------------------------------------------------------------------

      fprintf (fout, "%s %.12f %f %lu %f %f\n", filename, bestPeriod, (sum - bestScore) / stddiv, numRid, outAvr, outStdDiv) ;
      fflush (fout) ;

      free (time) ;
      free (mag) ;
      fclose (fin) ;
    }

  fclose (fout) ;
  fclose (finList) ;
  return (0) ;
}


/* todo:
   1. find best number of bins (r)
   2. or allow for variable size of bins  r = (int)sqrt(size)
   3. fold periodagaram ?
   4. post - fine tune in a bisection search
   5. add uncertainty info into the fit
*/
