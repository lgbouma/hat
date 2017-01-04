/**********************************************************
 *
 * This program simulates light curves, for full 31-parameter DEBiL files
 *
 * compile:  gcc showLC.c -o showLC -Wall -lm -O4
 *
 * run:      showLC ag.periods ag.sm
 *           
 *           sm
 *           : inp ag.sm
 *           : quit
 *
 *           ps2pdf ag.ps
 *           acroread ag.pdf &
 *
 **********************************************************/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>   // for strcmp() and strlen()


#define NUM_BOXES_X 4  // for SM file
#define NUM_BOXES_Y 6

#define WRITE_CURVE_FIT 500   // The number of fitted points to write

#define EPSILON 1.0e-9  // to avoid division by zero


// Returns the modulo 1 value of x
double mod1 (double x)
{
  return (x - floor(x)) ;
}

//-------------------------------------

int writeCurveFit (float *time, float *amp, float *errMag, int size, 
		   double period, char *localFilename)
{
  FILE *fout ;
  char foutName[256] ;
  int i ;

  sprintf (foutName, "%s.P%.3f.data", localFilename, period) ;
  fout = fopen (foutName, "wt") ;
  if (!fout)
    return (1) ;

  for (i = 0 ; i < size ; i++)
    {
      fprintf (fout, "%f %f\n", time[i], -2.5 * log10(amp[i])) ;
    }
  fclose (fout) ;

  return (0) ;
}


void printSM (FILE *foutSM, char *filename, double period)
{
  static int i = 0, j = NUM_BOXES_Y ;

  i++ ;

  if (i > NUM_BOXES_X)
    {
      i = 1 ;
      j-- ;
    }

  if (j <= 0)
    {
      j = NUM_BOXES_Y ;
      fprintf (foutSM, "page\n\n") ;
    }

  
  fprintf (foutSM, "window %d %d %d %d\n", NUM_BOXES_X, NUM_BOXES_Y, i, j) ;
  
  fprintf (foutSM, "data %s.P%.3f.data\n", filename, period) ;
  fprintf (foutSM, "read {x 1 y 2}\n") ;
  fprintf (foutSM, "limits 0 1 y\n") ;
  fprintf (foutSM, "limits 0 1 $fy2 $fy1\n") ; // invert axis direction
  // add margin:	   limits 0 1 $($fy2+0.2) $($fy1-0.1)
  
  fprintf (foutSM, "points x y\n") ;  
  fprintf (foutSM, "ticksize 0 0 0 0\n") ;
  fprintf (foutSM, "xlabel %s\n", filename) ;
  fprintf (foutSM, "ylabel P = %.3f days\n", period) ;
  fprintf (foutSM, "box\n\n") ;
  //------------------------------
  
  fprintf (foutSM, "set x = 0\n") ;
  fprintf (foutSM, "set y = 0\n") ;
}



int main (int argc, char **argv)
{
  FILE *fin, *finLC, *foutSM ;
  int size, i, res ;
  double period ;
  char filename[256], *localFilename ;
  double tmpTime, tmpMag, tmpErr ;
  float *time, *amp, *errMag ;

  if ((argc != 3) && (argc != 5))
    {
      printf ("usage: %s <LC list> <output SM file>\n", argv[0]) ;
      return (1) ;
    }

  if (!strcmp(argv[1], argv[2]))
    {
      printf ("ERROR: Input and output filenames must be different.\n") ;
      return (2) ;
    }

  if (!(fin = fopen (argv[1], "rt")))
    {
      printf ("ERROR: Couldn't open the DEBiL database ('%s').\n", argv[1]) ;
      return (3) ;
    }

  if (!(foutSM = fopen (argv[2], "wt")))
    {
      printf ("ERROR: Couldn't open the output SM file ('%s').\n", argv[2]) ;
      return (4) ;
    }

  //-------------------------------------

  fprintf (foutSM, "device postportfile %s.ps\n\n", argv[2]) ;
  fprintf (foutSM, "expand 0.4\n") ;
  fprintf (foutSM, "ptype 1 1\n\n") ;
  fprintf (foutSM, "ctype black\n") ;


  while (2 == (res = fscanf (fin, "%s %lf\n", filename, &period)))
    {
      printf ("%s\n", filename) ;

      finLC = fopen (filename, "rt") ;
 
      if (!finLC)
	{
	  printf ("Warning: Couldn't open LC file ('%s').\n", filename) ;
	  continue ;
	}

      size = 0 ;
      while (1 == fscanf(finLC, "%*f %*f %lf\n", &tmpErr))
	if (tmpErr > 0.0)  // Sign of an invalid magnitude
	  size++ ;
       
      time = (float*)malloc (size * sizeof(float)) ;
      amp = (float*)malloc (size * sizeof(float)) ;
      errMag = (float*)malloc (size * sizeof(float)) ;
      
      if ((!time) || (!amp) || (!errMag))
	{
	  if (time) free(time) ;
	  if (amp) free(amp) ;
	  if (errMag) free(errMag) ;

	  printf ("Warning: Not enough memory for loading '%s'.", filename) ;
	  fclose(finLC) ;
	  continue ;
	}

      //-------------------------

      rewind(finLC) ;

      i = 0 ;
      while ((i < size) && (3 == fscanf(finLC, "%lf %lf %lf\n", &tmpTime, &tmpMag, &tmpErr)))
	if (tmpErr > 0.0)  // Sign of an invalid magnitude
	  {
	    time[i] = (float)mod1(tmpTime / period) ;
	    amp[i] = (float)pow(10.0, -0.4 * tmpMag) ;  // Convert magnitude (logarithmic) to amplitude (linear)
	    errMag[i] = (float)tmpErr ;
	    i++ ;
	  }

      if (i != size)
	{
	  printf ("Warning: Number mismatch.\n") ;
	  free(time) ;
	  free(amp) ;
	  free(errMag) ;
	  fclose (finLC) ;
	  continue ;
	}

      //--------------------------

      i = strlen(filename) ;
      for (i-- ; (i >= 0) && (filename[i] != '/') ; i--) ;  // goes to slash or -1
      localFilename =  &filename[i+1] ;
      
      i = writeCurveFit (time, amp, errMag, size, period, localFilename) ;
       
      if (i)
	printf ("Warning: Couldn't write a *.data or *.fit file for '%s' (P=%.f3).", localFilename, period) ;
      else
	printSM (foutSM, localFilename, period) ;

      free(time) ;
      free(amp) ;
      free(errMag) ;
      fclose (finLC) ;
    }
   
  fprintf (foutSM, "hardcopy\n") ;

  if ((res > 0) || !feof(fin))
    printf ("\nWarning: Premature end of input file.\n") ;

  printf ("\nCreated files:\n") ;
  printf ("  *.data - light curve data\n") ;
  printf ("  %s - SM file (a SuperMongo script)\n", argv[2]) ;
  printf ("\nTo run the SM file, enter: 'sm', then 'inp %s', then 'quit'.\n", argv[2]) ;
  printf ("This will create a postscript file ('%s.ps') of the model fits catalogue.\n\n", argv[2]) ;


  fclose (fin) ;
  fclose (foutSM) ;
  return (0) ;
}
