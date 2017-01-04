/*********************************************************************
 *
 * This utility displays a list of light curves (LCs) in ASCII form, 
 * wrapped around (modulo) a given period
 *
 * compile:  gcc show.c -o show -Wall -lm 
 * run: show ../periodFinder/scripts/big.txt 0 big.comment.txt
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>   // for atol
#include <math.h>
#include <strings.h>


#define SCREEN_X 75    // for printing: < 78
#define SCREEN_Y 24    // for printing: 24
#define SCREEN_X_REPEAT 1   // (recommend: 1 or 2)
// note that if you use 2, you should lower SCREEN_X to half
// and/or set the shell window: Option-> Terminal...-> End-of-line Wrapping=off

#define SCREEN_Y_MARGIN 0.4   // [-sigma, sigma] --> [SCREEN_Y_MARGIN, 1-SCREEN_Y_MARGIN]
#define STDDIV_LIMIT 3.0
//#define DEBUG

//#define THINNING 19

unsigned char screenArr[SCREEN_X][SCREEN_Y] ;


inline double sqr (double x)
{
  return (x * x) ;
}

// Returns the modulo 1 value of x
double mod1 (double x)
{
  return (x - floor(x)) ;
}

//--------------------


// returns the reduced size of the array
// gaussian-normalized around [-1, 1]
unsigned long normalize (unsigned long size, double *mag, double *time)
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
  if (N <= 1) printf ("ERROR 1\n") ;
  
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
  if (N <= 1) printf ("ERROR 2\n") ;
  Sum /= N ;                // Sum --> Average
  Vsum = sqrt(Vsum / N) ;   // Vsum --> Std Div
  
  // step 4: normalize to [-1,1]
  for (i = 0 ; i < size ; i++)
    mag[i] = (mag[i] - Sum) / Vsum ;
  
  return (size) ;
}


void printPeriod (unsigned long size, double *time, double *mag, double period)
{
  int i, j, k, posX, posY ;
  double X, Y ;
    
  memset (screenArr, 0, SCREEN_X * SCREEN_Y * sizeof (unsigned char)) ;
  
  for (i = 0 ; i < size ; i++)
    {
      X = mod1(time[i] / period) ;
      posX = (int)(SCREEN_X * X) ;
      Y = (mag[i] * (0.5 - SCREEN_Y_MARGIN)) + 0.5 ;   // transforms:  [-1, 1] --> [SCREEN_Y_MARGIN, 1-SCREEN_Y_MARGIN]
      posY = (int)(SCREEN_Y * Y) ;
      
      if ((posX >= 0) && (posX < SCREEN_X) && (posY >= 0) && (posY < SCREEN_Y) && (screenArr[posX][posY] < 10))
	screenArr[posX][posY]++ ;
    }
  
  //--------------------- scatter --------------------------------------
  
  for (i = 0 ; i < SCREEN_Y ; i++)
    {
      for (k = 0 ; k < SCREEN_X_REPEAT ; k++)
	{
	  for (j = 0 ; j < SCREEN_X ; j++)
	    {
	      if (screenArr[j][i] == 0) printf (" ") ;
	      else if (screenArr[j][i] == 10) printf ("#") ;
	      else printf ("%c", '0' + screenArr[j][i]) ;
	    }

	  printf ("|") ;
	}
      
      printf ("\n") ;
    }
}


int main(const int argc, const char **argv)
{
  FILE *finList, *finLC, *fout = 0 ;
  double period ;
  char c, strUser[512], strFilename[512] ;
  double *time, *mag ;
  double tmpTime, tmpMag ;
  unsigned long numRid, size, counter = 0, offset = 0 ;

  if ((argc < 2) || (argc > 4))
    {
      printf ("%s <period list file> [offset] [output file]\n", argv[0]) ;
      return (1) ;
    }

  finList = fopen (argv[1], "rt") ;
  if (!finList)
    {
      printf ("ERROR: couldn't open input file '%s'\n", argv[1]) ;
      return (2) ;
    }

  if (argc >= 3)
    offset = atol (argv[2]) ;

  if (argc >= 4)
    {
      if (strcmp(argv[1], argv[3]) == 0)
	{
	  printf ("ERROR: the input and output files are the same ('%s')\n", argv[3]) ;
	  fclose (finList) ;
	  return (3) ;
	}

      fout = fopen (argv[3], "at") ;
      if (!fout)
	{
	  printf ("ERROR: couldn't open output file '%s'\n", argv[3]) ;
	  fclose (finList) ;
	  return (4) ;
	}
    }

  while (fscanf (finList, "%s %lf", strFilename, &period) == 2)
    {
      while ((c = getc(finList)) != '\n') 
	printf ("%c", c) ;  // rids remainind fields

      counter++ ;
      if (counter < offset) continue ;

#ifdef THINNING
      if (counter % THINNING) continue ; 
#endif

      finLC = fopen (strFilename, "rt") ;
      if (!finLC)
	{
	  printf ("ERROR: couldn't open file '%s'\n", strFilename) ;
	  getc(stdin) ;
	  continue ;
	}

      printf ("\n\n #%lu %s  period= %f \n\n", counter, strFilename, period) ;

      //------------

      size = 0 ;
      while (1 == fscanf(finLC, "%*f %*f %*f %lf %*f %*d\n", &tmpMag))
	if (tmpMag > 0)  // sign of invalid magnitude
	  size++ ;

      time = (double*)malloc (size *sizeof(double)) ;
      mag = (double*)malloc (size *sizeof(double)) ;
      
      if (!time || !mag)
	{
	  printf ("ERROR: not enought memory (size = %lu)\n", size) ;
	  if (time) free (time) ;
	  if (mag) free (mag) ;
	  fclose (finLC) ;
	  return (3) ;
	}
      
      size = 0 ;
      rewind(finLC) ;
      while (2 == fscanf(finLC, "%lf %*f %*f %lf %*f %*d\n", &tmpTime, &tmpMag)) 
	if (tmpMag > 0)  // sign of invalid magnitude
	  {
	    time[size] = tmpTime ;
	    mag[size] = tmpMag ;
	    size++ ;
	  }

      fclose (finLC) ;

      numRid = size ;
      size = normalize (size, mag, time) ;
      numRid -=size ;

      printPeriod (size, time, mag, period) ;

      free (mag) ;
      free (time) ;

      gets (strUser) ;

      if (strUser[0] == 'q') break ;

      if (fout)
	{
	  fprintf (fout, "%s %f", strFilename, period) ;

	  if (strUser[0]) fprintf (fout, "  user: %s", strUser) ;

	  fprintf (fout, "\n") ;
	}
      else if (strUser[0])
	{
	  printf ("Warning: your comment ('%s') is ignored since there is no output file\n", strUser) ;
	  getc (stdin) ;
	}
    }
  
  fclose (finList) ;
  if (fout) fclose (fout) ;		 
  return (0) ;
}
