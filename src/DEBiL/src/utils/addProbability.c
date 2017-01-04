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
 * compile:  gcc addProbability.c -o addProbability -lm -Wall -O3
 *
 * run: 
 *       makeEffTable 1000 1000 10000000 makeEffTable.res &
 *
 *       eff makeEffTable.res agc.out agcp.out
 *
 ****************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>  // for log10

#define MIN_PROBABILITY 1.0e-8

#define NUM_HIST 100


// calculates the 4-way interpolation (using all 4 neighboring values) of the given matrix
// returns the probability of an orbit with a given inclination (i) and eccentrisity (e)
double avr2d (double r, double e, double *arr, unsigned long numRa, unsigned long numEcc)
{
  double z00, z10, z01, z11, zmid ;
  int x = (int)(r * numRa) ;
  int y = (int)(e * numEcc) ;

  if (r == 1.0)   // contact limit
    return (1.0) ;

  if ((r > 1.0) || (e >= 1.0))
    {
      printf ("error 1  (r=%f ; e=%f)\n", r, e) ;
      return (1.0) ;
    }

  if (r < 0.0)
    {
      printf ("error 2\n") ;
      return (0.0) ;
    }

  if (e < 0.0)
    {
      printf ("error 3\n") ;
      return (r) ;
    }

  //-------------

  r = (r * numRa) - x ;  // normalize [0,1)
  e = (e * numEcc) - y ;

  z00 = arr[x + (y * (numRa+1))] ;
  z10 = arr[(x+1) + (y * (numRa+1))] ;
  z01 = arr[x + ((y+1) * (numRa+1))] ;
  z11 = arr[(x+1) + ((y+1) * (numRa+1))] ;

  zmid = 0.25 * (z00 + z10 + z01 + z11) ;

  if (r > e)
    {
      if ((r + e) < 1.0)
	return ((r * (z10-z00)) + (e * (zmid+zmid-z00-z10)) + z00) ;
      else
	return ((r * (z10+z11-zmid-zmid)) + (e * (z11-z10)) + (zmid+zmid-z11)) ;


    }
  else // r <= e
    {
      if ((r + e) < 1.0)
	return ((r * (zmid+zmid-z00-z01)) + (e * (z01-z00)) + z00) ;
      else
	return ((r * (z11-z01)) + (e * (z01+z11-zmid-zmid)) + (zmid+zmid-z11)) ;
    }
}



int main(int argc, char **argv)
{
  double p ;
  FILE *finTable, *fin, *fout, *fhist ;
  double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31 ;
  int i19, i20, val ;
  char str[256] ;
  unsigned long numEcc, numRa, ecc_count, r_a_count ;
  double *arr, *histPeriod ;

  if ((argc != 4) && (argc != 5))
    {
      printf ("Usage: %s <Eff table> <input DEBiL 31-param> <output DEBiL 31-param /w weight> [period histogram]\n", argv[0]) ;
      return (1) ;
    }

 
  finTable = fopen (argv[1], "rt") ;

  if (!finTable)
    {
      printf ("ERROR: couldn't open input table file '%s' for reading\n", argv[1]) ;
      return (2) ;
    }

  fin = fopen (argv[2], "rt") ;

  if (!fin)
    {
      printf ("ERROR: couldn't open input DEBiL file '%s' for reading\n", argv[2]) ;
      fclose (finTable) ;
      return (3) ;
    }

  fout = fopen (argv[3], "wt") ;

  if (!fout)
    {
      printf ("ERROR: couldn't open output file '%s' for writing\n", argv[3]) ;
      fclose (finTable) ;
      fclose (fin) ;
      return (4) ;
    }

  if (argc == 5)
    {
      fhist = fopen (argv[4], "wt") ;

      if (!fhist)
	{
	  printf ("ERROR: couldn't open histogram file '%s' for writing\n", argv[4]) ;
	  fclose (finTable) ;
	  fclose (fin) ;
	  fclose (fout) ;
	  return (5) ;
	}

      histPeriod = (double*)calloc(NUM_HIST, sizeof(double)) ;

      if (!histPeriod)
	{
	  printf ("ERROR: not enough memory for histPeriod\n") ;
	  fclose (finTable) ;
	  fclose (fhist) ;
	  fclose (fin) ;
	  fclose (fout) ;
	  return (6) ;
	}
    }
  else 
    {
      fhist = 0 ;
      histPeriod = 0 ; 
    }

  fscanf (finTable, "%lu  %lu\n", &numEcc, &numRa) ;

  arr = (double*)calloc ((numRa+1) * (numEcc+1), sizeof(double)) ;

  if (!arr)
    {
      printf ("ERROR: not enough memory for array\n") ;
      free (histPeriod) ;
      fclose (finTable) ;
      fclose (fhist) ;
      fclose (fin) ;
      fclose (fout) ;
      return (7) ;
    }

  //-----------------------------------------

  printf ("reading table...\n") ;
  
  for (ecc_count = 0 ; ecc_count < numEcc ; ecc_count++)
    {
      for (r_a_count = 0 ; r_a_count < numRa ; r_a_count++)
	fscanf (finTable, "%lf ", &arr[(ecc_count * (numRa+1)) + r_a_count]) ;
   
      fscanf (finTable, "\n") ;

      arr[(ecc_count * (numRa+1)) + numRa] = 1.0 ;
    }

  for (r_a_count = 0 ; r_a_count < numRa ; r_a_count++)
    arr[(numEcc * (numRa+1)) + r_a_count] = 1.0 ;

  //-----------

  /*
  for (ecc_count = 0 ; ecc_count <= numEcc ; ecc_count++)
    {
      for (r_a_count = 0 ; r_a_count <= numRa ; r_a_count++)
	printf ("%f ", arr[(ecc_count * (numRa+1)) + r_a_count]) ;
      
      printf ("\n") ;
    }
 */
   
  printf ("writing...\n") ;
  
  while (31 == fscanf (fin, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		       str, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		       &x17, &x18, &i19, &i20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, &x31))
    {
      p = avr2d(x5+x7 , x3, arr, numRa, numEcc) ; // (r1+r2, e, arr, numRa, numEcc)
      
      fprintf (fout, "%s %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f %f\n",
	       str, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
	       x17, x18, i19, i20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, 1.0 / p) ;
      
      if (histPeriod && fhist)
	{
	  val = (int)(log10(x2 * 10.0) * NUM_HIST / 4.0) ;      // transform [0.1, 1000] --> [0 , numHist]
	  if ((val >= 0) && (val < NUM_HIST))
	    {
	      if (p < MIN_PROBABILITY)
		p = MIN_PROBABILITY ;
	      
	      histPeriod[val] += (1.0 / p) ;
	    }
	}
    }
  
  if (histPeriod && fhist)
    {
      printf ("histogram...\n") ;
      
      for (val = 0 ; val < NUM_HIST ; val++)
	fprintf (fhist, "%f %f\n", ((float)val * 4.0 / NUM_HIST) - 1.0, histPeriod[val]) ;
    }

  printf ("done.\n") ;

  free (arr) ;
  fclose (fin) ;
  fclose (fout) ;
  fclose (finTable) ;
  if (histPeriod) free (histPeriod) ;
  if (fhist) fclose (fhist) ;
  return (0) ;
}


