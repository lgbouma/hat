/************************************************************
 *
 * This utility filters the DEBiL (/w color) database
 * where one column must be more or less another
 *
 * Compile:  gcc filterDebil2.c -o filterDebil2 -Wall
 *
 * Run:      grep -v 1000.0 ag.chiAll.out >! ag.chiAll.noNaN.out
 *           filterDebil2 ag.chiAll.noNaN.out 22 less 25 >! ag.notSin.out
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
  int numColumn1, numColumn2, doLess ;
  double val1, val2 ;

  if (argc != 5)
    {
      printf (" Usage:  %s  <DEBiL database with color> <column num> <more / less> <column num>\n", argv[0]) ;
      return (1) ;
    }
 
  numColumn1 = atoi (argv[2]) ;
  numColumn2 = atoi (argv[4]) ;

  if ((numColumn1 <= 0) || (numColumn2 <= 0) || (numColumn1 > 42) || (numColumn2 > 42))
    {
      printf ("ERROR: Column number is out of range (%d %d)\n", numColumn1, numColumn2) ;
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

      switch (numColumn1)
	{
	case 1: val1 = x1 ; break ;
	case 2: val1 = x2 ; break ;
	case 3: val1 = x3 ; break ;
	case 4: val1 = x4 ; break ;
	case 5: val1 = x5 ; break ;
	case 6: val1 = x6 ; break ;
	case 7: val1 = x7 ; break ;
	case 8: val1 = x8 ; break ;
	case 9: val1 = x9 ; break ;
	case 10: val1 = x10 ; break ;
	case 11: val1 = x11 ; break ;
	case 12: val1 = x12 ; break ;
	case 13: val1 = x13 ; break ;
	case 14: val1 = x14 ; break ;
	case 15: val1 = x15 ; break ;
	case 16: val1 = x16 ; break ;
	case 17: val1 = x17 ; break ;
	case 18: val1 = x18 ; break ;
	case 19: val1 = x19 ; break ;
	case 20: val1 = x20 ; break ;
	case 21: val1 = x21 ; break ;
	case 22: val1 = x22 ; break ;
	case 23: val1 = x23 ; break ;
	case 24: val1 = x24 ; break ;
	case 25: val1 = x25 ; break ;
	case 26: val1 = x26 ; break ;
	case 27: val1 = x27 ; break ;
	case 28: val1 = x28 ; break ;
	case 29: val1 = x29 ; break ;
	case 30: val1 = x30 ; break ;
	case 31: val1 = x31 ; break ;
	case 32: val1 = x32 ; break ;
	case 33: val1 = x33 ; break ;
	case 34: val1 = x34 ; break ;
	case 35: val1 = x35 ; break ;
	case 36: val1 = x36 ; break ;
	case 37: val1 = x37 ; break ;
	case 38: val1 = x38 ; break ;
	case 39: val1 = x39 ; break ;
	case 40: val1 = x40 ; break ;
	case 41: val1 = x41 ; break ;
	case 42: val1 = x42 ; break ;
	case 45: val1 = x45 ; break ;
	}
      
      switch (numColumn2)
	{
	case 1: val2 = x1 ; break ;
	case 2: val2 = x2 ; break ;
	case 3: val2 = x3 ; break ;
	case 4: val2 = x4 ; break ;
	case 5: val2 = x5 ; break ;
	case 6: val2 = x6 ; break ;
	case 7: val2 = x7 ; break ;
	case 8: val2 = x8 ; break ;
	case 9: val2 = x9 ; break ;
	case 10: val2 = x10 ; break ;
	case 11: val2 = x11 ; break ;
	case 12: val2 = x12 ; break ;
	case 13: val2 = x13 ; break ;
	case 14: val2 = x14 ; break ;
	case 15: val2 = x15 ; break ;
	case 16: val2 = x16 ; break ;
	case 17: val2 = x17 ; break ;
	case 18: val2 = x18 ; break ;
	case 19: val2 = x19 ; break ;
	case 20: val2 = x20 ; break ;
	case 21: val2 = x21 ; break ;
	case 22: val2 = x22 ; break ;
	case 23: val2 = x23 ; break ;
	case 24: val2 = x24 ; break ;
	case 25: val2 = x25 ; break ;
	case 26: val2 = x26 ; break ;
	case 27: val2 = x27 ; break ;
	case 28: val2 = x28 ; break ;
	case 29: val2 = x29 ; break ;
	case 30: val2 = x30 ; break ;
	case 31: val2 = x31 ; break ;
	case 32: val2 = x32 ; break ;
	case 33: val2 = x33 ; break ;
	case 34: val2 = x34 ; break ;
	case 35: val2 = x35 ; break ;
	case 36: val2 = x36 ; break ;
	case 37: val2 = x37 ; break ;
	case 38: val2 = x38 ; break ;
	case 39: val2 = x39 ; break ;
	case 40: val2 = x40 ; break ;
	case 41: val2 = x41 ; break ;
	case 42: val2 = x42 ; break ;
	case 45: val2 = x45 ; break ;
	}

  
      if ((doLess & (val1 < val2)) ||
	  (!doLess & (val1 > val2)))
	printf ("%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %s %f\n",
		x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
		x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, 
		x31, x32, x33, x34, x35, x36, x37, x38, x39, x40, x41, x42, x43, x44, x45) ;
    }

fclose (fin) ;
return (0) ;
}
