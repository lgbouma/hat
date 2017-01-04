/************************************************************
 *
 * This utility filters the DEBiL (31-param) database
 * if the stellar radius of one of the binary conponent is larger/smaller
 * than the Roche limit (i.e semi-detached or contact)
 *
 * Compile:  gcc filterRoche.c -o filterRoche -Wall -lm
 *
 * Run:   filterRoche ag.chi4.out 0 > ag.roche.out
 *
 ************************************************************/

/* Column description:

1.  File name (with full path if not in local directory)
2.  Period
3.  Orbital eccentricity
4.  Absolute error in (3)
5.  Radius of large star (in units of semi-major axis)
6.  Absolute error in (5)
7.  Radius of small star (in units of semi-major axis)
8.  Absolute error in (7)
9.  Brightness of large star (I-band magnitudes)
10. Absolute error in (9)
11. Brightness of small star (I-band magnitudes)
12. Absolute error in (11)
13. Sine of inclination - sin(i)
14. Absolute error in (13)
15. Epoch of perihelion (modulo the period, in units of the period)
16. Absolute error in (15)
17. Argument of perihelion (in degrees)
18. Absolute error in (17)
19. Number of used data points (not including outliers)
20. Number of outliers
21. Reduced chi2 of model
22. Reduced chi2 relative to a constant, set to the average value
23. Reduced chi2 relative to a second order spline
24. Reduced chi2 relative to a sinusoidal best fit (min of mode 1 & 2)
25. Significance of the secondary dip depth (in sigma)
26. Significance of the hump height (at midpoint between dips in sigma)
27. Significance of the hump difference between the two humps
28. Waviness (see manual)
29. Scatter score (see manual)
30. Mean density (grams per cm^3 ; see manual)
31. Max density (grams per cm^3 ; see manual)
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double rRoche (double q)
{
 double logq, Ra ;

 logq = log10(q) ;

 if (q > 0.1)
   Ra = 0.37771 + (0.20247 * logq) + (0.01838 * logq * logq) + (0.02275 * logq * logq * logq) ;
 else
   Ra = 0.3771 + (0.2131 * logq)  - (0.008 * logq * logq) + (0.0066 * logq * logq * logq) ;

  return (Ra) ;
}


int main(int argc, char **argv)
{
  FILE *fin ;
  double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31 ;
  int x19, x20 ;
  char x1[256] ;
  int numFilled, countFilled ;

  if (argc != 3)
    {
      printf (" Usage:  %s  <DEBiL database (31-param)> <num filled lobes (0,1,2)>\n", argv[0]) ;
      return (1) ;
    }

  numFilled = atoi (argv[2]) ;

  if ((numFilled < 0) || (numFilled > 2))
  {
    printf ("ERROR: Invalid number of filled lobes (%d)\n", numFilled) ;
    return (2) ;
  }

  fin = fopen (argv[1], "rt") ;

  if (!fin)
    {
      printf ("ERROR: couldn't open file '%s'\n", argv[1]) ;
      return (3) ;
    }  
  
  
  while (31 == fscanf (fin, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		       x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		       &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, &x31))
    {
      countFilled = 0 ;

      if (x5 >= rRoche(pow (x5 / x7, 1.534)))  // early type main sequence:  R ~ M^0.652
	countFilled++ ;

      if (x7 >= rRoche(pow (x7 / x5, 1.534)))
	countFilled++ ;
  
      if (countFilled == numFilled)
	printf ("%s %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f\n",
		x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
		x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31) ;
    }

fclose (fin) ;
return (0) ;
}
