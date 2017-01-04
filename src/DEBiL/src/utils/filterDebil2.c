/************************************************************
 *
 * This utility filters the DEBiL (31-param) database
 * where one column must be more or less another
 *
 * Compile:  gcc filterDebil2.c -o filterDebil2 -Wall
 *
 * Run:      grep -v 1000.0 ag.chiAll.out >! ag.chiAll.noNaN.out
 *           filterDebil2 ag.chiAll.noNaN.out 22 less 25 >! ag.notSin.out
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

int main(int argc, char **argv)
{
  FILE *fin ;
  double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31 ;
  int x19, x20 ;
  char x1[256] ;
  int numColumn1, numColumn2, doLess ;
  double val1, val2 ;

  if (argc != 5)
    {
      printf (" Usage:  %s  <DEBiL database (31-param)> <column num> <more / less> <column num>\n", argv[0]) ;
      return (1) ;
    }
 
  numColumn1 = atoi (argv[2]) ;
  numColumn2 = atoi (argv[4]) ;

  if ((numColumn1 <= 1) || (numColumn2 <= 1) || (numColumn1 > 31) || (numColumn2 > 31))
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
  
  
  while (31 == fscanf (fin, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		       x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		       &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, &x31))
    {

      switch (numColumn1)
	{
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
	}
      
      switch (numColumn2)
	{
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
	}

      if ((doLess & (val1 < val2)) || (!doLess & (val1 > val2)))
	printf ("%s %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f\n",
		x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
		x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31) ;
    }

fclose (fin) ;
return (0) ;
}
