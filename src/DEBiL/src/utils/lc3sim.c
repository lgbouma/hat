/******************************************
 *
 * This program simulates light curves
 *
 * compile:  gcc lc3sim.c -o lc3sim -Wall -lm
 * run:      lc3sim lc3sim.txt 1.0 0.03 0.4 0.3 11.0 12.0 0.95 0.1 0.6
 *
 *  <output LC filename> <period> <e> <r1> <r2> <B1> <B2> <sin(i)> <time0> <omega>
 *****************************************/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_DATA_POINTS 1500

// time between samples in units of period (shouldn't be near a "small value" rational number)
#define SAMPLING_PERIOD 0.005413  // 0.3343

#define MAX_SAMPLING_ERROR 0.0  // fractional error in sample time

#define MAX_FLUX_ERROR 0.01  // fractional error in flux (must be positive)

//Eccentric anomaly calculation parameter:
#define ECC_ANOMALY_MAX_ERROR 1.0e-4  // Convergence threshold (average error is less than a third of this) 

#define INTEGRATION_STEP 1e-4

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

  if (r1 + r2 > 1.0 - e)
    {
      printf ("ERROR f1 %f %f %f", r1, r2, e) ;
      return (0.0) ;
    }

  E = eccentricAnomaly(meanAnomaly, e) ;

  D2 = sqr(1.0 - (e * cos(E))) - sqr((((cos(E) - e) * sin(omega)) + (sqrt(1.0 - (e * e)) * sin(E) * cos (omega))) * sin_i) ;

  if (D2 < 0.0)
    {
      printf ("ERROR f2  (%f)", D2) ;
      return (0.0) ;
    }

  D = sqrt (D2) ;

  if (D >= r1 + r2)   // side by side
    return ((B1 * integrateWholeDisk(r1)) + (B2 * integrateWholeDisk(r2))) ;
  
  if (sin(E + omega) > 0.0)  // is star 2 in front?
    return ((B2 * integrateWholeDisk(r2)) + (B1 * integrateBackOverlap(r2, r1, D, dr))) ;
  
  // star 1 is in front
  return ((B1 * integrateWholeDisk(r1)) + (B2 * integrateBackOverlap(r1, r2, D, dr))) ;
}

//-------------------------------------------------------


int main (int argc, char **argv)
{
  int i ;
  FILE *foutLC ;
  double time, F, period, time0, sin_i, r1, r2, B1, B2, e, omega ;
  double quadA = DEFAULT_LIMB_QUAD_A, quadB = DEFAULT_LIMB_QUAD_B ;

  /*--------------
  argc = 13 ;
  argv[1] = "lc3sim.txt" ;
  argv[2] = "1.0" ;  // period
  argv[3] = "0.08" ; // e
  argv[4] = "0.22" ;  // r1
  argv[5] = "0.02" ;  // r2
  argv[6] = "12.5" ;  // B1 (mag)
  argv[7] = "15.0" ;  // B2 (mag)
  argv[8] = "0.998" ; // sin_i
  argv[9] = "0.5" ;  // time0
  argv[10] = "50.0" ; // omega
  argv[11] = "0.3542" ; // quadratic limb darkening param a
  argv[12] = "0.1939" ; // quadratic limb darkening param b
  ---------------------------*/

  if ((argc != 11) && (argc != 13))
    {
      printf ("usage: %s <output LC filename> <period> <e> <r1> <r2> <B1> <B2> <sin(i)> <time0> <omega> [<limb a> <limb b>]\n", argv[0]) ;
      return (1) ;
    }

  if (!(foutLC = fopen (argv[1], "wt")))
    {
      printf ("couldn't open the output LC file") ;
      return (2) ;
    }

  if (argc == 13) // this block must come before integrateWholeDisk()
    {
      quadA = atof(argv[11]) ;
      quadB = atof(argv[12]) ;
    }

  setLimbDarkening (quadA, quadB) ;

  period = atof(argv[2]) ;
  e = atof(argv[3]) ;
  r1 = atof(argv[4]) ;
  r2 = atof(argv[5]) ;
  B1 = pow(10.0, -0.4 * atof(argv[6])) / integrateWholeDisk(r1) ;
  B2 = pow(10.0, -0.4 * atof(argv[7])) / integrateWholeDisk(r2) ;
  sin_i = atof(argv[8]) ;
  time0 = atof(argv[9]) ;
  omega = atof(argv[10]) * M_PI / 180.0 ;

  if (argc == 13)
    {
      quadA = atof(argv[11]) ;
      quadB = atof(argv[12]) ;
    }

  setLimbDarkening (quadA, quadB) ;

  for (i = 1 ; i <= NUM_DATA_POINTS ; i++)
    {
      time = i * SAMPLING_PERIOD  ;
      time += ((2.0 * MAX_SAMPLING_ERROR * (float)rand() / RAND_MAX) - MAX_SAMPLING_ERROR) * time ;

      F = flux (time, time0, sin_i, r1, r2, B1, B2, e, omega, INTEGRATION_STEP) ;
      F += ((2.0 * MAX_FLUX_ERROR * (float)rand() / RAND_MAX) - MAX_FLUX_ERROR) * F ;

      // (2.5 / (sqrt(3) * ln(10))) converts stddiv to error magnitude
      //  fprintf (foutLC, " %f  %e  0.0  %f  %f  0\n", time * period, F, -2.5 * log10(F), 0.62685009 * MAX_FLUX_ERROR) ;  // OGLE format
      fprintf (foutLC, " %f  %f  %f\n", time * period, -2.5 * log10(F), 0.62685009 * MAX_FLUX_ERROR) ;
    }

  fclose (foutLC) ;
  return (0) ;
}
