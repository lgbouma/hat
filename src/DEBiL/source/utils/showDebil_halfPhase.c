/**********************************************************
 *
 * This program simulates light curves, for full 31-parameter DEBiL files
 *
 * compile:  gcc showDebil_halfPhase.c -o showDebil_halfPhase -Wall -lm -O4
 *
 * run:      showDebil_halfPhase ag.dbl ag.sm
 *           
 *           sm
 *           : inp ag.sm
 *           : quit
 *
 *           ps2pdf ag.ps
 *           acroread ag.pdf &
 *
 **********************************************************/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>   // for strcmp() and strlen()


#define NUM_BOXES_X 1  // for SM file
#define NUM_BOXES_Y 6

#define WRITE_CURVE_FIT 10000   // The number of fitted points to write

#define INTEGRATION_STEP 0.002  // acuracy of light curve simulation

//Eccentric anomaly calculation parameter:
#define ECC_ANOMALY_MAX_ERROR 1.0e-4  // Convergence threshold (average error is less than a third of this) 

#define EPSILON 1.0e-9  // to avoid division by zero

// General bisection implementation parameter:
// only needed for pathological cases (e.g. overflows or huge roundoff errors)
#define BISECTION_SAFETY_MAX_ITERATIONS 1000 

#define DEFAULT_LIMB_QUAD_A 0.3542 // 0.282   [Claret 2003 // Claret 1998]
#define DEFAULT_LIMB_QUAD_B 0.1939 // 0.311


double sqr (double x)
{
  return (x * x) ;
}

double cube (double x)
{
  return (x * x * x) ;
}


// Returns the modulo 1 value of x
double mod1 (double x)
{
  return (x - floor(x)) ;
}

//-------------------------------------

// Uses a combined Newton-Raphson + bisection method for finding the root (E) of:
//  M = E - (e * sin(E))     [Kepler's equation]
// Given: M = mean anomaly ; e = eccentricity of orbit (0 <= e < 1)
// Returns: E = Eccentric Anomaly
double eccentricAnomaly (double M, double e)
{
  int counter = BISECTION_SAFETY_MAX_ITERATIONS ;
  double Emin = M - 1.0, Emax = M + 1.0, E = M ;
  double f, dfm, dE ;

  do
    {
      f = M + (e * sin(E)) - E ;  // May be optimized by setting cos=sqrt(1-sin^2)
      dfm = 1.0 - (e * cos(E)) ;  // Minus differential of f (always non-negative)
      dE = dfm * ECC_ANOMALY_MAX_ERROR ;

      if (f > dE)
	{
	  Emin = E ;
	  dE = Emax - Emin ;

	  if ((dfm * dE) > f)
	    E += (f / dfm) ;
	  else
	    E = 0.5 * (Emax + Emin) ;
	}
      else if (f < (-dE))
	{
	  Emax = E ;
	  dE = Emax - Emin ;

	  if ((dfm * dE) > -f)
	    E += (f / dfm) ;
	  else
	    E = 0.5 * (Emax + Emin) ;
	}
      else return (E) ;
    }
  while ((dE > ECC_ANOMALY_MAX_ERROR) && ((--counter) > 0)) ;

  return (E) ;  // should almost never exits here
}


//========================== Limb Darkening ====================================

// Quadratic limb darkening coefficients for a solar-like star in I filter (Claret 1998)
// Note that a "square root" law give a better fit (especially for I filter), but
// but the error in assuming a solar star, dominates by far, and the total CPU requirements
// would almost double (doing a linear fit, on the other hand, won't save much CPU).

static double limbConstA = 0.0 ;
static double limbConstB = 0.0 ;
static double limbConstC = 0.0 ;
static double limbConstD = 0.0 ;
static double limbConstE = 0.0 ;
static double limbConstF = 0.0 ;
static double limbConstG = 0.0 ;
static double limbConstH = 0.0 ;


// sets the quadratic limb darkening parameters to arbitrary given values.
void setLimbDarkening (double quadA, double quadB)
{
  limbConstA = 1.0 - quadA - quadB ;
  limbConstB = quadA + quadB + quadB ;
  limbConstC = -quadB ;
  limbConstD = M_PI * (1.0 - quadA - quadB - quadB) ;
  limbConstE = 0.5 * M_PI * quadB ;
  limbConstF = 2.0 * M_PI * (quadA + quadB + quadB) / 3.0 ;
  limbConstG = M_PI * (1.0 - quadA - quadB) ;
  limbConstH = M_PI * (6.0 - quadA - quadA - quadB) / 6.0 ;
}


// cosAngle2 = the square of the cosine of the angle between the zenith of
// the star to the given point on the surface (relative to the star's center)
// returns the solar limb darkening coefficient
// note that it is normalized so that at the star's zenith (cosAngle2 = 1), it returns 1
double limbDarkening (double cosAngle2)
{
  return (limbConstA + (limbConstB * sqrt(cosAngle2)) + (limbConstC * cosAngle2)) ;
}


// Analytically integrates the limb darkened disk, from r=0 to r=rBack
// rBack = the full radius of the star's disk
double integrateWholeDisk (double rBack)
{
  return (limbConstH * rBack * rBack) ;
}


// Analytically integrates the limb darkened disk, from 0 to the given finish
// finish = where to end the integration
// rBack = the radius of the back star's disk (the one seen as a crescent)
double integrateDiskFromStart (double finish, double rBack)
{
  const double x = sqr(finish / rBack) ;

  return (rBack * rBack * (((limbConstD + (limbConstE * x)) * x) + (limbConstF * (1.0 - cube(sqrt(1.0 - x)))))) ;
}


// Analytically integrates the limb darkened disk, from the given start to rBack
// start = from where to begin the integration
// rBack = the radius of the back star's disk (the one seen as a crescent)
double integrateDiskToFinish (double start, double rBack)
{
  const double x = 1.0 - sqr(start / rBack) ;

  return (rBack * rBack * (((limbConstG - (limbConstE * x)) * x) + (limbConstF * cube(sqrt(x))))) ;
}



// Makes a simple trapezoid numerical integration with varying step size
// (smaller, the closer an asimptote) for finding the brightness of
// a partially hidden star. This method works better than Simpson's.
// start = where to begin the integration
// finish = where to end the integration
// rFront = the radius of the front star's disk (the one seen as a disk)
// rBack = the radius of the back star's disk (the one seen as a crescent)
// D = the distance between the centers of the two stars
// dr0 = the infinitesimal integration step scale
// Assumes:  D + rFront >= finish >= start >= |D - rFront| >= 0
// [num iterations] = 4 / dr0   ;  [error] = 0.1 * rFront^2 * dr0^2
double integrateCrescent (double start, double finish, double rFront, double rBack,
			  double D, double dr0)
{
  const double rBack_m2 = 1.0 / (rBack * rBack) ;
  const double rFront2mD2 = (rFront * rFront) - (D * D) ;
  const double middle = 0.5 * (start + finish) ;
  const double D2 = 2.0 * D ;
  const double limitUp = ((D + rFront < rBack) ? D + rFront : rBack) ;
  const double limitDown = fabs(D - rFront) ;
  double dr, r, r2, rStep, sum = 0.0 ;

  if (D < EPSILON)  // do avoid division by zero
    return (0.0) ;

  dr0 *= sqrt(rFront) ;  // since dr0 is unitless and we want dr to have
  // units of "distance" (1 = sum semi-major axese)

  // 1. integrate from middle up
  dr = dr0 * sqrt(limitUp - middle) ;
  for (rStep = middle + dr ; rStep < finish ; rStep += dr)
    {
      r = rStep - (0.5 * dr) ;
      r2 = r * r ;
      sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening(1.0 - (rBack_m2 * r2)) * dr) ;

      dr = dr0 * sqrt(limitUp - r) ;
    }

  // 2. add sliver at upper edge
  r = 0.5 * (finish + rStep - dr) ;
  r2 = r * r ;
  dr += (finish - rStep) ;
  sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening(1.0 - (rBack_m2 * r2)) * dr) ;

  // 3. integrate from middle down
  dr = dr0 * sqrt(middle - limitDown) ;
  for (rStep = middle - dr ; rStep > start ; rStep -= dr)
    {
      r = rStep + (0.5 * dr) ;
      r2 = r * r ;
      sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening(1.0 - (rBack_m2 * r2)) * dr) ;

      dr = dr0 * sqrt(r - limitDown) ;
    }

  // 4. add sliver at bottom edge
  r = 0.5 * (start + rStep + dr) ;
  r2 = r * r ;
  dr += (rStep - start) ;
  sum += (r * acos((rFront2mD2 - r2) / (D2 * r)) * limbDarkening(1.0 - (rBack_m2 * r2)) * dr) ;

  return (2.0 * sum) ;
}

// =============== calculate the observed radiation flux from the binary system ===============

// rFront = the radius of the star in front
// rBack = the radius of the star in back
// D = the distance between the centers of the two stars
// dr = the infinitesimal integration step
// returns the relative brightness of the crescent of the bigger star,
// when it is partially eclipsed by the smaller star
double integrateBackOverlap (double rFront, double rBack, double D, double dr)
{
  if (rFront < D)
    {
      if (rBack > D + rFront)
	return (integrateDiskFromStart (D - rFront, rBack) +
		integrateCrescent (D - rFront, D + rFront, rFront, rBack, D, dr) +
		integrateDiskToFinish (D + rFront, rBack)) ;                  // E
      
      return (integrateDiskFromStart (D - rFront, rBack) +
	      integrateCrescent (D - rFront, rBack, rFront, rBack, D, dr)) ;  // B
    }

  if (rFront > D + rBack)
    return (0.0) ;                                                            // D
  
  if (rBack > D + rFront)
    return (integrateCrescent (rFront - D, rFront + D, rFront, rBack, D, dr) +
	    integrateDiskToFinish (rFront + D, rBack)) ;                      // F
  
  return (integrateCrescent (rFront - D, rBack, rFront, rBack, D, dr)) ;      // C
}


/******************************************************************
 * returns the total expected flux from both stars of the binary
 * ignores limb darkening, gravity darkening, etc.
 * time = scan parameter 0 <= time < 1
 * p[DIM] = parameter vector
 * dr = the infinitesimal integration step
 * doB12max = 1:recalculate the light curve "plateau" brightness ; 0:don't
 *
 * NOTE: [definition of 1 unit distance]
 *  the physical distance between the stars at perihelion (minimal) is (1-e) units
 *  thus the physical distance between the stars at aphelion (maximal) is (1+e) units
 *  in other words, 1 unit is the sum of the distances of the two semi-major axes
 *
 ******************************************************************/
double flux (double time, double time0, double sin_i,
	     double r1, double r2, const double B1, const double B2,
	     double e, double omega, double dr)
{
  const double meanAnomaly = 2.0 * M_PI * (time - time0) ;
  double D2, D, E ;

  E = eccentricAnomaly(meanAnomaly, e) ;

  D2 = sqr(1.0 - (e * cos(E))) - sqr((((cos(E) - e) * sin(omega)) + (sqrt(1.0 - (e * e)) * sin(E) * cos(omega))) * sin_i) ;

  D = sqrt (D2) ;

  if (D >= r1 + r2)   // side by side
    return ((B1 * integrateWholeDisk(r1)) + (B2 * integrateWholeDisk(r2))) ;
  
  if (sin(E + omega) > 0.0)  // is star 2 in front?
    return ((B2 * integrateWholeDisk(r2)) + (B1 * integrateBackOverlap(r2, r1, D, dr))) ;
  
  // star 1 is in front
  return ((B1 * integrateWholeDisk(r1)) + (B2 * integrateBackOverlap(r1, r2, D, dr))) ;
}

//-------------------------------------------------------


int writeCurveFit (float *time, float *amp, float *errMag, int size, 
		   double period, double time0, double sin_i,
		   double r1, double r2, const double B1, const double B2,
		   double e, double omega, char *localFilename)
{
  FILE *fout ;
  char foutName[256] ;
  int i ;
  double F ;

  sprintf (foutName, "%s.P%.3f.data", localFilename, period) ;
  fout = fopen (foutName, "wt") ;
  if (!fout)
    return (1) ;

  for (i = 0 ; i < size ; i++)
    {
      F = flux (time[i], time0, sin_i, r1, r2, B1, B2, e, omega, INTEGRATION_STEP) ;

      fprintf (fout, "%f %f %f %f\n", fmod(time[i], 0.5), -2.5 * log10(amp[i]), 
	       errMag[i], 2.5 * log10(F / amp[i])) ;
    }
  fclose (fout) ;

  //---------------------

  sprintf (foutName, "%s.P%.3f.fit", localFilename, period) ;
  fout = fopen (foutName, "wt") ;
  if (!fout)
    return (2) ;

  for (i = 0 ; i < WRITE_CURVE_FIT ; i++)
    {
      F = flux ((float)i / WRITE_CURVE_FIT, time0, sin_i, r1, r2, B1, B2, e, omega, INTEGRATION_STEP) ;
      fprintf (fout, "%f %f\n", fmod((float)i / WRITE_CURVE_FIT, 0.5), -2.5 * log10(F)) ;
    }

  fclose (fout) ;
  return (0) ;
}

//-------------------------------------------------------


void printSM (FILE *foutSM, char *filename, double period)
{
  static int i = 0, j = NUM_BOXES_Y ;

  i++ ;

  if (i > NUM_BOXES_X)
    {
      i = 1 ;
      j -= 2 ;
    }

  if (j <= 0)
    {
      j = NUM_BOXES_Y ;
      fprintf (foutSM, "page\n\n") ;
    }

  
  fprintf (foutSM, "window %d %d %d %d\n", NUM_BOXES_X, NUM_BOXES_Y, i, j) ;
  
  fprintf (foutSM, "data %s.P%.3f.fit\n", filename, period) ;
  fprintf (foutSM, "read {xf 1 yf 2}\n") ;
  fprintf (foutSM, "data %s.P%.3f.data\n", filename, period) ;
  fprintf (foutSM, "read {x 1 y 2 r 4}\n") ;
  fprintf (foutSM, "limits 0 0.5 y\n") ;
  fprintf (foutSM, "limits 0 0.5 $fy2 $fy1\n") ; // invert axis direction
  // add margin:	   limits 0 1 $($fy2+0.2) $($fy1-0.1)
  
  fprintf (foutSM, "ctype blue\n") ;
  fprintf (foutSM, "points xf yf\n") ;
  fprintf (foutSM, "ctype black\n") ;
  fprintf (foutSM, "points x y\n") ;
  
  fprintf (foutSM, "ticksize 0 0 0 0\n") ;
  fprintf (foutSM, "xlabel %s\n", filename) ;
  fprintf (foutSM, "ylabel P = %.3f days\n", period) ;
  fprintf (foutSM, "box\n\n") ;
  //------------------------------
  
  fprintf (foutSM, "window %d %d %d %d\n", NUM_BOXES_X, NUM_BOXES_Y, i, j-1) ;
  fprintf (foutSM, "limits 0 0.5 0.3 -0.3\n") ; // invert axis direction
  fprintf (foutSM, "points x r\n") ;
  
  fprintf (foutSM, "ticksize 0.05 0.2 0.05 0.1\n") ;
  fprintf (foutSM, "xlabel %s\n", filename) ;
  fprintf (foutSM, "ylabel Residuals\n") ;
  fprintf (foutSM, "box\n\n") ;

  fprintf (foutSM, "set xf = 0\n") ;
  fprintf (foutSM, "set yf = 0\n") ;
  fprintf (foutSM, "set x = 0\n") ;
  fprintf (foutSM, "set y = 0\n") ;
  fprintf (foutSM, "set r = 0\n") ;
}



int main (int argc, char **argv)
{
  FILE *fin, *finLC, *foutSM ;
  double quadA = DEFAULT_LIMB_QUAD_A, quadB = DEFAULT_LIMB_QUAD_B ;
  int size, i, res ;
  double period, e, r1, r2, B1, B2, sin_i, time0, omega ;
  char filename[256], *localFilename ;
  double tmpTime, tmpMag, tmpErr ;
  float *time, *amp, *errMag ;

  if ((argc != 3) && (argc != 5))
    {
      printf ("usage: %s <DEBiL database (31-param)> <output SM file> [<limb darkening A> <limb darkening B>]\n", argv[0]) ;
      return (1) ;
    }

  if (!strcmp(argv[1], argv[2]))
    {
      printf ("ERROR: Input and output filenames must be different.\n") ;
      return (2) ;
    }

  if (!(fin = fopen (argv[1], "rt")))
    {
      printf ("ERROR: Couldn't open the DEBiL database ('%s').\n", argv[1]) ;
      return (3) ;
    }

  if (!(foutSM = fopen (argv[2], "wt")))
    {
      printf ("ERROR: Couldn't open the output SM file ('%s').\n", argv[2]) ;
      return (4) ;
    }

  if (argc > 3)
    {
      quadA = atof(argv[3]) ;
      quadB = atof(argv[4]) ;
    }

  setLimbDarkening (quadA, quadB) ;

  printf ("Setting limb darkening coefficients:  A=%f   B=%f\n", quadA, quadB) ;

  //-------------------------------------

  fprintf (foutSM, "device postportfile %s.ps\n\n", argv[2]) ;
  fprintf (foutSM, "expand 0.4\n") ;
  fprintf (foutSM, "ptype 1 1\n\n") ;


  while (10 == (res = fscanf (fin, "%s %lf %lf %*f %lf %*f %lf %*f %lf %*f %lf %*f %lf %*f %lf %*f %lf %*f %*d %*d %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f %*f\n", filename, &period, &e, &r1, &r2, &B1, &B2, &sin_i, &time0, &omega)))
    {
      B1 = pow(10.0, -0.4 * B1) / integrateWholeDisk(r1) ;
      B2 = pow(10.0, -0.4 * B2) / integrateWholeDisk(r2) ;
      omega = omega * M_PI / 180.0 ;
      
      printf ("%s\n", filename) ;

      finLC = fopen (filename, "rt") ;
 
      if (!finLC)
	{
	  printf ("Warning: Couldn't open LC file ('%s').\n", filename) ;
	  continue ;
	}

      size = 0 ;
      while (1 == fscanf(finLC, "%*f %*f %lf\n", &tmpErr))
	if (tmpErr > 0.0)  // Sign of an invalid magnitude
	  size++ ;
       
      time = (float*)malloc (size * sizeof(float)) ;
      amp = (float*)malloc (size * sizeof(float)) ;
      errMag = (float*)malloc (size * sizeof(float)) ;
       
      if ((!time) || (!amp) || (!errMag))
	{
	  if (time) free(time) ;
	  if (amp) free(amp) ;
	  if (errMag) free(errMag) ;

	  printf ("Warning: Not enough memory for loading '%s'.", filename) ;
	  fclose(finLC) ;
	  continue ;
	}

      //-------------------------

      rewind(finLC) ;

      i = 0 ;
      while ((i < size) && (3 == fscanf(finLC, "%lf %lf %lf\n", &tmpTime, &tmpMag, &tmpErr)))
	if (tmpErr > 0.0)  // Sign of an invalid magnitude
	  {
	    time[i] = (float)mod1(tmpTime / period) ;
	    amp[i] = (float)pow(10.0, -0.4 * tmpMag) ;  // Convert magnitude (logarithmic) to amplitude (linear)
	    errMag[i] = (float)tmpErr ;
	    i++ ;
	  }

      if (i != size)
	{
	  printf ("Warning: Number mismatch.\n") ;
	  free(time) ;
	  free(amp) ;
	  free(errMag) ;
	  fclose (finLC) ;
	  continue ;
	}

      //--------------------------

      i = strlen(filename) ;
      for (i-- ; (i >= 0) && (filename[i] != '/') ; i--) ;  // goes to slash or -1
      localFilename =  &filename[i+1] ;
      
      i = writeCurveFit (time, amp, errMag, size, period, time0, sin_i, r1, r2, B1, B2, e, omega, localFilename) ;
       
      if (i)
	printf ("Warning: Couldn't write a *.data or *.fit file for '%s' (P=%.f3).", localFilename, period) ;
      else
	printSM (foutSM, localFilename, period) ;

      free(time) ;
      free(amp) ;
      free(errMag) ;
      fclose (finLC) ;
    }
   
  fprintf (foutSM, "hardcopy\n") ;

  if ((res > 0) || !feof(fin))
    printf ("\nWarning: Premature end of input file.\n") ;

  printf ("\nCreated files:\n") ;
  printf ("  *.fit - fitted light curve model\n") ;
  printf ("  *.data - light curve data and model residuals\n") ;
  printf ("  %s - SM file (a SuperMongo script)\n", argv[2]) ;
  printf ("\nTo run the SM file, enter: 'sm', then 'inp %s', then 'quit'.\n", argv[2]) ;
  printf ("This will create a postscript file ('%s.ps') of the model fits catalogue.\n\n", argv[2]) ;


  fclose (fin) ;
  fclose (foutSM) ;
  return (0) ;
}
