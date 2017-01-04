/************************************************************
 *
 * This utility takes two DEBiL (31-param) databases and picks the best fit of each light curve;
 * It's specifically designed for collating ecc / non-ecc debil fits.
 * Note that we assume that the two DBs contain the same light curves in the same order.
 * 
 *
 * Compile:  gcc collate2lc.c -o collate2lc -Wall -lm
 *
 * Run:      collate2lc listGood5.0.split.dbl listGood5.0.ecc.dbl > listGood5.0.collate.dbl
 *           collate2lc listGood5.0.split.dbl listGood5.0.ecc.dbl 2.0 > listGood5.0.collate.f2.dbl
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
#include <string.h>
#include <stdlib.h>

int main(int argc, char **argv)
{
  FILE *finx, *finy ;
  double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31 ;
  int x19, x20 ;
  double y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16 ;
  double y17, y18, y21, y22, y23, y24, y25, y26, y27, y28, y29, y30, y31 ;
  int y19, y20 ;
  char x1[256], y1[256] ;
  float pref1Factor ;

  if ((argc != 3) && (argc != 4))
    {
      printf (" Usage:  %s  <DEBiL database #1 (31-param)>  <DEBiL database #2 (31-param)>  [factor for preffering #1]\n", argv[0]) ;
      return (1) ;
    }

  finx = fopen (argv[1], "rt") ;
  if (!finx)
    {
      printf ("ERROR: couldn't open file #1: '%s'\n", argv[1]) ;
      return (2) ;
    }  

  finy = fopen (argv[2], "rt") ;
  if (!finy)
    {
      printf ("ERROR: couldn't open file #2: '%s'\n", argv[2]) ;
      fclose (finx) ;
      return (3) ;
    }  

  if (argc > 3)
    pref1Factor = atof (argv[3]) ;
  else
    pref1Factor = 1.0 ;

    
  while ((31 == fscanf (finx, 
			"%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		       x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		       &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, &x31)) &&
	 (31 == fscanf (finy, 
			"%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		    y1, &y2, &y3, &y4, &y5, &y6, &y7, &y8, &y9, &y10, &y11, &y12, &y13, &y14, &y15, &y16, 
		    &y17, &y18, &y19, &y20, &y21, &y22, &y23, &y24, &y25, &y26, &y27, &y28, &y29, &y30, &y31)) &&
	 (strcmp (x1, y1) == 0))
    {

      if (x21 < (pref1Factor * y21))
	printf ("%s %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f\n",
		x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
		x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31) ;
      else
	printf ("%s %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f\n",
		y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, 
		y17, y18, y19, y20, y21, y22, y23, y24, y25, y26, y27, y28, y29, y30, y31) ;
    }

  if (31 == fscanf (finx,
		    "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		     x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		     &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, &x31))
    printf ("file #1 ('%s') didn't reach its end.\n", argv[1]) ;


  if (31 == fscanf (finy,
		    "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		    y1, &y2, &y3, &y4, &y5, &y6, &y7, &y8, &y9, &y10, &y11, &y12, &y13, &y14, &y15, &y16, 
		    &y17, &y18, &y19, &y20, &y21, &y22, &y23, &y24, &y25, &y26, &y27, &y28, &y29, &y30, &y31))
    printf ("file #2 ('%s') didn't reach its end.\n", argv[2]) ;
    
  fclose (finx) ;
  fclose (finy) ;
  return (0) ;
}
