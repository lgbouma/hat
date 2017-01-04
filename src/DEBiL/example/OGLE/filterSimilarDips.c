/************************************************************
 *
 * This utility filters out the DEBiL (/w color) database
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

1.  OGLE II bulge field number
2.  Light curve number (unique to field)
3.  Period
4.  Orbital eccentricity
5.  Absolute error in (4)
6.  Radius of large star (in units of semi-major axis)
7.  Absolute error in (6)
8.  Radius of small star (in units of semi-major axis)
9.  Absolute error in (8)
10. Brightness of large star (I-band magnitudes)
11. Absolute error in (10)
12. Brightness of small star (I-band magnitudes)
13. Absolute error in (12)
14. Sine of inclination - sin(i)
15. Absolute error in (14)
16. Epoch of perihelion (modulo the period, in units of the period)
17. Absolute error in (16)
18. Argument of perihelion (in degrees)
19. Absolute error in (18)
20. Number of used data points (not including outliers)
21. Number of outliers
22. Chi2 of model
23. Chi2 relative to the average value (should be worse than the model)
24. Chi2 relative to a second order spline
25. Chi2 relative to a sinusoidal best fit
26. Fitness score (position of (22) between (23) and (24))
27. Significance of the secondary dip depth (in sigma)
28. Significance of the hump height (at midpoint between dips in sigma)
29. Significance of the hump difference between the two humps
30. Waviness (see manual)
31. Scatter score (see manual)
32. Mean density (grams per cm^3 ; see manual)
33. Max density (grams per cm^3 ; see manual)
34. Extinction corrected brightness of large star (I-band magnitudes)
35. Extinction corrected brightness of small star (I-band magnitudes))
36. Extinction corrected total binary brightness (I-band magnitudes)
37. Extinction corrected I-band magnitude, from Udalski's catalog - similar to (34)
38. Extinction corrected color (V-I), from Udalski's catalog
39. Uncorrected I-band magnitude
40. Uncorrected color (V-I)
41. Absolute error in V-band magnitude
42. Absolute error in I-band magnitude
43. Binary right-ascension (RA)
44. Binary declination (Dec)
45. Inverse of probablity
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
  double x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x19, x22, x23, x24, x25, x26, x27, x28, x29, x30 ;
  double x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x45 ;
  int x1, x2, x20, x21 ;
  char x43[256], x44[256] ;
  double logJp, logJs, epsilon_logJpJs, thresholdMidHump = 0.0 ;
  
  if ((argc != 2) && (argc != 3))
    {
      printf (" Usage:  %s  <DEBiL database with color>  [max mid hump]\n", argv[0]) ;
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


  while (45 == fscanf (fin, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %s %lf\n",
		       &x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		       &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, 
		       &x31, &x32, &x33, &x34, &x35, &x36, &x37, &x38, &x39, &x40, &x41, &x42, x43, x44, &x45))
    {
     
      logJp = (-0.4 * x10) - (2.0 * log10(x6)) ; 
      logJs = (-0.4 * x12) - (2.0 * log10(x8)) ;
      epsilon_logJpJs = sqrt((0.16 * (sqr(x11) + sqr(x13))) + (0.754446788 * (sqr(x7/x6) + sqr(x9/x8)))) ; // 0.754446788 = sqr(2 / ln(10))

      // printf ("== %d  %d %d %d\n", (x28 > thresholdMidHump), (x4 > x5) , (x28 > 0.0) , (fabs(logJp - logJs) > epsilon_logJpJs)) ;

      if (((argc < 3) || (x28 > thresholdMidHump)) && ((x4 > x5) || (x28 > 0.0) || (fabs(logJp - logJs) > epsilon_logJpJs)))
	printf ("%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %f\n",
		x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
		x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, 
		x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45) ;
    }

fclose (fin) ;
return (0) ;
}
