/************************************************************
 *
 * This utility filters the DEBiL (31-param) database
 * above or below a given threshold for some column
 *
 * Compile:  gcc filterDebil.c -o filterDebil -Wall -lm
 *
 * Run:      grep -v 1000.0 ag.chi4.out >! ag.chi4.noNaN.out
 *           filterDebil ag.chi4.noNaN.out 2 less 10.0 >! ag.filter.out
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

int main(int argc, char **argv)
{
  FILE *fin ;
  double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31 ;
  int x19, x20 ;
  char x1[256] ;
  int numColumn, doComp ;
  double threshold, val ;

  if (argc != 5)
    {
      printf (" Usage:  %s  <DEBiL database (31-param)> <column num> <more / equal/ less> <threshold>\n", argv[0]) ;
      printf ("\nSpecial column num = 101: combined magnitude\n") ;
      return (1) ;
    }
 
  numColumn = atoi (argv[2]) ;

  if ((numColumn <= 1) || ((numColumn > 31) && (numColumn != 101)))
    {
      printf ("ERROR: Column number is out of range (%d)\n", numColumn) ;
      return (2) ;
    }

  if ((argv[3][0] == 'm') || (argv[3][0] == 'M'))
    doComp = 0 ;
  else if ((argv[3][0] == 'l') || (argv[3][0] == 'L'))
    doComp = 1 ; 
  else if ((argv[3][0] == 'e') || (argv[3][0] == 'E'))
    doComp = 2 ; 
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
  
  
  while (31 == fscanf (fin, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		       x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		       &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, &x31))
    {
      switch (numColumn)
	{
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
	case 101: val = -2.5 * log10(pow(10.0, -0.4 * x9) + pow(10.0, -0.4 * x11)) ;
	}
  
      if (((doComp == 0) && (val > threshold)) || 
	  ((doComp == 1) && (val < threshold)) ||
	  ((doComp == 2) && (val == threshold)))
	printf ("%s %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f\n",
		x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, 
		x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31) ;
    }

fclose (fin) ;
return (0) ;
}
