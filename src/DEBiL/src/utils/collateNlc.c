/************************************************************
 *
 * This utility takes any number of DEBiL (31-param) databases and picks the best fit of each light curve;
 * It's specifically designed for collating ecc / non-ecc debil fits.
 * Note that we assume that the two DBs contain the same light curves in the same order.
 * 
 * Compile:  gcc collateNlc.c -o collateNlc -Wall -lm
 *
 * Run:      
 *
 * foreach file (*.dbl)
 * echo $file
 * cat $file >>! concat.dbl
 * end
 *
 * sort concat.dbl > concat.sort.dbl
 *
 * collateNlc concat.sort.dbl > collate.dbl
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


#define EPSILON 0.1


// return: 1= if the light curves are different
//         0= if they are the same
int isDifferent (char *name1, double period1, char *name2, double period2)
{
  return (((period1 - period2) > (period2 * EPSILON)) || ((period2 - period1) > (period1 * EPSILON)) || (strcmp (name1, name2) != 0)) ;
}


int main(int argc, char **argv)
{
  FILE *fin ;
  double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31 ;
  int x19, x20 ;
  double y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16 ;
  double y17, y18, y21, y22, y23, y24, y25, y26, y27, y28, y29, y30, y31 ;
  int y19, y20 ;
  char x1[256], y1[256] ;
  double x21_f, y21_f ;
  float pref1Factor ;

  if ((argc != 2) && (argc != 3))
    {
      printf (" Usage:  %s  <sorted DEBiL database (31-param)>  [factor for preffering (e=0)]\n\n", argv[0]) ;
      return (1) ;
    }

  fin = fopen (argv[1], "rt") ;
  if (!fin)
    {
      printf ("ERROR: couldn't open file ('%s')\n", argv[1]) ;
      return (2) ;
    }  

 if (argc > 2)
    pref1Factor = atof (argv[2]) ;
  else
    pref1Factor = 1.0 ;

  y1[0] = 0 ; // makes it an empty sting
    
  while (31 == fscanf (fin, 
		       "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		       x1, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		       &x17, &x18, &x19, &x20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, &x31))
    {
      if (x3 != 0.0)
	x21_f = x21 * pref1Factor ;
      else
	x21_f = x21 ;


      if ((y1[0] != 0) && isDifferent (x1, x2, y1, y2))
	printf ("%s %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f\n",
		y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, 
		y17, y18, y19, y20, y21, y22, y23, y24, y25, y26, y27, y28, y29, y30, y31) ;
      
      if (isDifferent (x1, x2, y1, y2) || (x21_f < y21_f))
	{
	  strcpy (y1, x1) ;
	  y2 = x2 ;
	  y3 = x3 ;
	  y4 = x4 ;
	  y5 = x5 ;
	  y6 = x6 ;
	  y7 = x7 ;
	  y8 = x8 ;
	  y9 = x9 ;
	  y10 = x10 ;
	  y11 = x11 ;
	  y12 = x12 ;
	  y13 = x13 ;
	  y14 = x14 ;
	  y15 = x15 ;
	  y16 = x16 ;
	  y17 = x17 ;
	  y18 = x18 ;
	  y19 = x19 ;
	  y20 = x20 ;
	  y21 = x21 ;
	  y22 = x22 ;
	  y23 = x23 ;
	  y24 = x24 ;
	  y25 = x25 ;
	  y26 = x26 ;
	  y27 = x27 ;
	  y28 = x28 ;
	  y29 = x29 ;
	  y30 = x30 ;
	  y31 = x31 ;
	  y21_f = x21_f ;
	}
    }

  if (y1[0] != 0)
    printf ("%s %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f\n",
	    y1, y2, y3, y4, y5, y6, y7, y8, y9, y10, y11, y12, y13, y14, y15, y16, 
	    y17, y18, y19, y20, y21, y22, y23, y24, y25, y26, y27, y28, y29, y30, y31) ;
  
  fclose (fin) ;
  return (0) ;
}
