/************************************************************
 *
 * This utility filters the DEBiL (/w color) database
 * above or below a given threshold for some column
 *
 * Compile:  gcc filterDebil.c -o filterDebil -Wall
 *
 * Run:      grep -v 1000.0 ag.chi4.out >! ag.chi4.noNaN.out
 *           filterDebil ag.chi4.noNaN.out 2 less 10.0 >! ag.filter.out
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
#include <stdlib.h>

int main(int argc, char **argv)
{
  FILE *fin ;
  double x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x19, x22, x23, x24, x25, x26, x27, x28, x29, x30 ;
  double x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x45 ;
  int x1, x2, x20, x21 ;
  char x43[256], x44[256] ;
  int numColumn, doLess ;
  double threshold, val ;

  if (argc != 5)
    {
      printf (" Usage:  %s  <DEBiL database with color> <column num> <more / less> <threshold>\n", argv[0]) ;
      return (1) ;
    }
 
  numColumn = atoi (argv[2]) ;

  if ((numColumn <= 0) || (numColumn > 42))
    {
      printf ("ERROR: Column number is out of range (%d)\n", numColumn) ;
      return (2) ;
    }

  if (argv[3][0] == 'm')
    doLess = 0 ;
  else if (argv[3][0] == 'l')
    doLess = 1 ; 
  else
    {
      printf ("ERROR: Use 'more' or 'less'  (not '%s')\n", argv[3]) ;
      return (3) ;
    }

  threshold = atof (argv[4]) ;

  fin = fopen (argv[1], "rt") ;

  if (!fin)
    {
      printf ("ERROR: couldn't open file '%s'\n", argv[1]) ;
      return (4) ;
    }  
  
  
  while (45 == fscanf (fin, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %s %lf\n",
		       &x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		       &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, 
		       &x31, &x32, &x33, &x34, &x35, &x36, &x37, &x38, &x39, &x40, &x41, &x42, x43, x44, &x45))
    {

      switch (numColumn)
	{
	case 1: val = x1 ; break ;
	case 2: val = x2 ; break ;
	case 3: val = x3 ; break ;
	case 4: val = x4 ; break ;
	case 5: val = x5 ; break ;
	case 6: val = x6 ; break ;
	case 7: val = x7 ; break ;
	case 8: val = x8 ; break ;
	case 9: val = x9 ; break ;
	case 10: val = x10 ; break ;
	case 11: val = x11 ; break ;
	case 12: val = x12 ; break ;
	case 13: val = x13 ; break ;
	case 14: val = x14 ; break ;
	case 15: val = x15 ; break ;
	case 16: val = x16 ; break ;
	case 17: val = x17 ; break ;
	case 18: val = x18 ; break ;
	case 19: val = x19 ; break ;
	case 20: val = x20 ; break ;
	case 21: val = x21 ; break ;
	case 22: val = x22 ; break ;
	case 23: val = x23 ; break ;
	case 24: val = x24 ; break ;
	case 25: val = x25 ; break ;
	case 26: val = x26 ; break ;
	case 27: val = x27 ; break ;
	case 28: val = x28 ; break ;
	case 29: val = x29 ; break ;
	case 30: val = x30 ; break ;
	case 31: val = x31 ; break ;
	case 32: val = x32 ; break ;
	case 33: val = x33 ; break ;
	case 34: val = x34 ; break ;
	case 35: val = x35 ; break ;
	case 36: val = x36 ; break ;
	case 37: val = x37 ; break ;
	case 38: val = x38 ; break ;
	case 39: val = x39 ; break ;
	case 40: val = x40 ; break ;
	case 41: val = x41 ; break ;
	case 42: val = x42 ; break ;
        case 45: val = x45 ; break ;
	}
  
      if ((doLess & (val < threshold)) ||
	  (!doLess & (val > threshold)))
	printf ("%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %f\n",
		x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
		x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, 
		x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45) ;
    }

fclose (fin) ;
return (0) ;
}
