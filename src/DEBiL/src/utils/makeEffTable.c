/*****************************************************
 *
 * Written by Jonathan Devor  (Dec. 2004)
 *
 * This program calculates (numerically/monte-carlo) the geometric
 * selection effect of eccentric orbits. The final result
 * is stored in a 2-dimentional matrix ((r1+r2)/a ; e) and the
 * value in the array is the probablity that this orbit will
 * be seen as an ecclipsing binary.
 * Note that if the stars touch, the orbit will always be
 * eclipsing. This happens when:  (r1+r2)/a > (1-e)
 *
 * The final result is a (4-way) interpolation of the matrix
 *
 * O( numIter * numEcc )
 *
 * compile:  gcc makeEffTable.c -o makeEffTable -lm -O3 -Wall
 *
 * run:      makeEffTable 100 100 10000 makeEffTable.res &
 *        or
 *           makeEffTable 1000 1000 10000000 makeEffTable.res &
 *
 *           eff_OGLE makeEffTable.res agc.out agcp.out
 *
 ****************************************************/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define INIT_STEP_SIZE 0.1

#define SEED 1977

#define MIN_PROBABILITY 1.0e-8


double sqr (double x)
{
  return (x * x) ;
}


// the distance between two point in a 2-body eccentric orbit at a spesific moment.
// given their eccentrisity (e), angle of perihelion (sin omega, cos_omega),
// inclination angle (sin_i) and mean anomoly (E)
double dist2 (double e, double se, double sin_omega, double cos_omega, double sin_i, double E)
{
  double cosE = cos(E) ;

  return (sqr(1.0 - (e * cosE)) - sqr((((cosE - e) * sin_omega) + (se * sin(E) * cos_omega)) * sin_i)) ;
}


// find minimal projected distance between two points in a 2-body eccentric orbit
double findMinDist2 (double e, double omega, double i)
{
  double x, x1 = 0.0, x2, x3, y, y1, y2, y3, subStep = 0.5 ;
  double se = sqrt (1.0 - (e*e)) ;
  double sin_omega = sin(omega) ;
  double cos_omega = cos(omega) ;
  double sin_i = sin(i) ;
  int loop ;

  y1 = 2.0 ; // should be larger than 1.0

  for (x = 2.0 * M_PI ; x > 0.0 ; x -= INIT_STEP_SIZE)
    {
      y = dist2 (e, se, sin_omega, cos_omega, sin_i, x) ;

      if (y < y1)
	{
	  x1 = x ;
	  y1 = y ;
	}
    }

  for (loop = 0 ; loop < 3 ; loop++, subStep *= 0.1)
    {
      x2 = x1 + (subStep * INIT_STEP_SIZE) ;
      x3 = x1 - (subStep * INIT_STEP_SIZE) ;
      y2 = dist2 (e, se, sin_omega, cos_omega, sin_i, x2) ;
      y3 = dist2 (e, se, sin_omega, cos_omega, sin_i, x3) ;

      x = (x1 * (y2-y3)) + (x2 * (y3-y1)) + (x3 * (y1-y2)) ;

      if (x == 0.0)
	return (y1) ;

      x1 = 0.5 * ((x1 * x1 * (y2-y3)) + (x2 * x2 * (y3-y1)) + (x3 * x3 * (y1-y2))) / x ;
      y1 = dist2 (e, se, sin_omega, cos_omega, sin_i, x1) ;
    }

  return (y1) ;
}



void makeArr (unsigned long numRa, unsigned long numEcc, unsigned long numIter, unsigned long *arr)
{
  int r_a_count, ecc_count ;
  double x, omega, i, e ;
  unsigned long iter ;

  // e = 0 ;
  for (r_a_count = 1 ; r_a_count <= numRa ; r_a_count++)
    arr[r_a_count] = (unsigned long)(((double)numIter * r_a_count / numRa) + 0.5) ;

  // 0 < e < 1
  srand (SEED) ;
  for (ecc_count = 1 ; ecc_count <= numEcc ; ecc_count++)
    {
      e = ((double)ecc_count) / numEcc ;

      printf ("%f\n", e) ;

      for (iter = 0 ; iter < numIter ; iter++)
	{
	  omega = 2.0 * M_PI * rand() / (((double)RAND_MAX) + 1.0) ;
	  i = acos(rand() / ((double)RAND_MAX)) ;

	  x = findMinDist2 (e, omega, i) ;

	  if (x > 0.0)
	    x = sqrt(x) ;
	  else
	    x = 0.0 ;

	  for (r_a_count = (int)ceil(x * numRa) ; r_a_count <= numRa ; r_a_count++)
	    arr[(ecc_count * (numRa+1)) + r_a_count]++ ;
	}
    }
}



int main(int argc, char **argv)
{
  FILE *fout ;
  long numIter, numEcc, numRa, ecc_count, r_a_count ;
  unsigned long *arr ;


  if (argc != 5)
    {
      printf ("Usage: %s <num ecc. columns> <num (R/a) rows> <num iter> <output table file>\n", argv[0]) ;
      return (1) ;
    }

  numEcc = atol (argv[1]) ;
  numRa = atol (argv[2]) ;
  numIter = atol (argv[3]) ;

  if ((numEcc < 2) || (numRa < 2) || (numIter < 2))
    {
      printf ("ERROR: At least one of the input values is too small (<2)\n") ;
      return (2) ;
    }

  arr = (unsigned long*)calloc ((numRa+1) * (numEcc+1), sizeof(unsigned long)) ;

  if (!arr)
    {
      printf ("ERROR: not enough memory\n") ;
      return (3) ;
    }

  fout = fopen (argv[4], "wt") ;

  if (!fout)
    {
      printf ("ERROR: couldn't open output file '%s' for writing\n", argv[4]) ;
      free (arr) ;
      return (4) ;
    }

  //-----------------------------------------

  makeArr (numRa, numEcc, numIter, arr) ;


  fprintf (fout, "%lu  %lu\n", numEcc, numRa) ;

  for (ecc_count = 0 ; ecc_count < numEcc ; ecc_count++)
    {
      for (r_a_count = 0 ; r_a_count < numRa ; r_a_count++)
	fprintf (fout, "%f ", (double)arr[(ecc_count * (numRa+1)) + r_a_count] / numIter) ;
      
      fprintf (fout, "\n") ;
    } 

 
  free(arr) ;
  fclose(fout) ;
  return (0) ;
}


