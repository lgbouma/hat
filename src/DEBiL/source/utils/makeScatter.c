/*********************************************************************
 *
 * This utility makes the all-inclusive scatter plot data for the light curves
 *
 * compile:  gcc makeScatter.c -o makeScatter -Wall -lm 
 * run: makeScatter allBig.txt >! scatter.txt
 *
 **********************************************************************/

#include <stdio.h>
#include <math.h>

#define STDDIV_LIMIT 3.0
#define MAX_SIZE 500

//#define DEBUG



inline double sqr (double x)
{
  return (x * x) ;
}


// returns the reduced size of the array
unsigned long measure (unsigned long size, double *mag, double *time, double *outAvr, double *outStdDiv)
{
  unsigned long N, i ;
  double Sum, Vsum ;
  
  // step 1: first pass
  N = 1 ;
  Sum = mag[0] ;
  Vsum = 0.0 ;
  for (i = 1 ; i < size ; i++)
    {
      N++ ;
      Sum += mag[i] ;
      Vsum += (sqr((N * mag[i]) - Sum) / (N * (N-1))) ;
    }
  
  Sum /= N ;   // Sum --> Average
  Vsum = sqrt(Vsum / N) ;   // Vsum --> Std Div
  
#ifdef DEBUG
  printf (" avr = %f   std div = %f\n", Sum, Vsum) ;
#endif
  
  // step 2: prune outlier
#ifdef STDDIV_LIMIT
  i = 0 ;
  while (i < size)
    {
      if (fabs(mag[i] - Sum) > (Vsum * STDDIV_LIMIT))
	{
#ifdef DEBUG
	  printf ("outlier: %f  \n", mag[i]) ;
#endif
	  size-- ;
	  time[i] = time[size] ;
	  mag[i] = mag[size] ;
	}
      else i++ ;
    }
#endif
  
  // step 3: second pass
  N = 1 ;
  Sum = mag[0] ;
  Vsum = 0.0 ;
  for (i = 1 ; i < size ; i++)
    {
      N++ ;
      Sum += mag[i] ;
      Vsum += (sqr((N * mag[i]) - Sum) / (N * (N-1))) ;
    }

  //--------------------------

  *outAvr = Sum / N ;             // Sum --> Average
  *outStdDiv = sqrt(Vsum / N) ;   // Vsum --> Std Div
  
  return (size) ;
}


int main(const int argc, const char **argv)
{
  FILE *finList, *finLC ;
  double period, conf ;
  char filename[512] ;
  double time[MAX_SIZE], mag[MAX_SIZE] ;
  double tmpTime, tmpMag, outAvr, outStdDiv ;
  unsigned long numRid, size ;

  if (argc != 2)
    {
      printf ("%s <period list file>\n", argv[0]) ;
      return (1) ;
    }

  finList = fopen (argv[1], "rt") ;
  if (!finList)
    {
      printf ("ERROR: couldn't open input file '%s'\n", argv[1]) ;
      return (2) ;
    }

  while (fscanf (finList, "%s period= %lf confidence= %lf\n", filename, &period, &conf) == 3)
    {
      if (conf < 20.0) continue ;

      finLC = fopen (filename, "rt") ;
      if (!finLC)
	{
	  printf ("ERROR: couldn't open file '%s'\n", filename) ;
	  getc(stdin) ;
	  continue ;
	}

      //------------

      size = 0 ;
      while (2 == fscanf(finLC, "%lf %lf %*f\n", &tmpTime, &tmpMag))
	{
	  if (size == MAX_SIZE)
	    {
	      //  printf ("ERROR: there are more than %u data points\n", MAX_SIZE-1) ;
	      break ;
	    }
	  
	  if (tmpMag > 0)  // sign of invalid magnitude
	    {
	      time[size] = tmpTime ;
	      mag[size] = tmpMag ;
	      size++ ;
	    }
	}

      fclose (finLC) ;

      numRid = size ;
      size = measure (size, mag, time, &outAvr, &outStdDiv) ;
      numRid -=size ;

      // paramaters to play with:  period, conf, numRid, outAvr, outStdDiv

      printf ("%s %f %f %lu %f %f\n", filename, period, conf, numRid, outAvr, outStdDiv) ;
    }
  
  fclose (finList) ;
  return (0) ;
}
