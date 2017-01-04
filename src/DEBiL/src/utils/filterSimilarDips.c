/************************************************************
 *
 * This utility filters out the DEBiL (31-param) database
 * light curves which could be erroneous double-true-period. This happens
 * when one of the dips is not detected, and causes to identical dips to form.
 * This filter takes these light curves out, acknowledging that some of them
 * may be correct (e.g. same stellar types, orbiting in a circular orbit)
 *
 * Compile:  gcc filterSimilarDips.c -o filterSimilarDips -Wall -lm
 *
 * Run:   filterSimilarDips ag.chi4.out -0.5 > ag.noSimilarDips.out
 *     or
 *        filterSimilarDips ag.chi4.out > ag.noSimilarDips.out
 *        //filterSimilarDips ag.chi4.out -0.5 | grep == | sort | uniq -c   (verbose)
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
#include <stdlib.h>  // for atof()
#include <math.h>


double sqr (double x)
{
  return (x * x) ;
}


int main(int argc, char **argv)
{
  FILE *fin ;
  double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31 ;
  int x19, x20 ;
  char x1[256] ;
  double logJp, logJs, epsilon_logJpJs, thresholdMidHump = 0.0 ;
  
  if ((argc != 2) && (argc != 3))
    {
      printf (" Usage:  %s  <DEBiL database (31-param)>  [max mid hump]\n", argv[0]) ;
      return (1) ;
    }

  fin = fopen (argv[1], "rt") ;

  if (!fin)
    {
      printf ("ERROR: couldn't open file '%s'\n", argv[1]) ;
      return (3) ;
    }  
  
  if (argc == 3)
    thresholdMidHump = atof (argv[2]) ;


 while (31 == fscanf (fin, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		       x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		      &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, &x31))
   {
     
     logJp = (-0.4 * x9) - (2.0 * log10(x5)) ; 
     logJs = (-0.4 * x11) - (2.0 * log10(x7)) ;
     epsilon_logJpJs = sqrt((0.16 * (sqr(x10) + sqr(x12))) + (0.754446788 * (sqr(x6/x5) + sqr(x8/x7)))) ; // 0.754446788 = sqr(2 / ln(10))
     
     // printf ("== %d  %d %d %d\n", (x26 > thresholdMidHump), (x3 > x4) , (x26 > 0.0) , (fabs(logJp - logJs) > epsilon_logJpJs)) ;
     
     if (((argc < 3) || (x26 > thresholdMidHump)) && ((x3 > x4) || (x26 > 0.0) || (fabs(logJp - logJs) > epsilon_logJpJs)))
       printf ("%s %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f\n",
	       x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
	       x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31) ;
   }
 
 fclose (fin) ;
 return (0) ;
}
