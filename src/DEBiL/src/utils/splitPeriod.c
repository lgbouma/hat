/******************************************************************************************
 * This program is based on DEBiL
 * It is meant to check if the found light curve period is an eclipsing binary, and whether a 
 * half/double-period possibility should be added for the DEBiL list. After running splitPeriod,
 * you should run DEBiL WITHOUT period doubling.
 *
 * Compile:  gcc splitPeriod.c -o splitPeriod -Wall -lm -O4
 * Run: splitPeriod periodList.out periodList.sp.txt  periodList.sp.err
 *
 *
 * ~/work/source/periodS2 LC_list.txt periodList.txt periodList.err
 * ~/work/source/periodFilter periodList.txt 4 1000 > periodList.out
 * splitPeriod periodList.out periodList.sp.txt >! periodList.sp.err
 * ~/work/source/debil  periodList.sp.txt  periodList.sp.dbl  periodList.sp.dbl_err
 *******************************************************************************************/

// Should be much smaller than (KERNEL_HALF_WIDTH / size), especially when the fit is crude
#define SEARCH_STEP_SIZE  0.0001   

// Smoothing kernel parameter:
#define NUM_POINTS_HALF_KERNEL 8  // Must be smaller than a half a dip width, but at least 2

// Outlier identification parameter:
#define OUTLIER_VARIANCE_LIMIT 8.0  // "sigma clipping" 


// Error thresholds:
#define MIN_STDDIV_BELOW_MEDIAN 2.5    // threshold for doubling the period (CAUTION)
#define MAX_STDDIV_ABOVE_MEDIAN 2.5    // midway hump
#define MIN_DIP_SIG_DIFF 2.0           // when the dips are significantly different

// Minimum number of data points in the dips and in the plateau
#define MIN_NUM_DATA 5  // Must be at least a third of DIM

#define EPSILON 1.0e-9  // To avoid division by zero

// A useful constant for magnitude error conversion (2.5 / ln(10))
#define MAGNITUDE_ERROR_FACTOR 1.085736205

// ==================== Includes ===================================

#include <math.h>
#include <stdio.h>  // For getc()
#include <stdlib.h> // For malloc() + qsort()
#include <values.h> // For MAXFLOAT
#include <string.h> // For memmove() + strlen() + strcmp()


// ================ General utility functions ======================


// Returns the square of a given x
double sqr (double x)
{
  return (x * x) ;
}

// Returns the cube of a given x
double cube (double x)
{
  return (x * x * x) ;
}


// Returns the modulo 1 value of x
double mod1 (double x)
{
  return (x - floor(x)) ;
}


// Swaps two given values
void swap (double *x1, double *x2)
{
  double tmp = *x1 ;
  *x1 = *x2 ;
  *x2 = tmp ;
}


// =================== Removing outliers ==========================

// Performs one iteration of statistical summations
void updateVars (float x, float y,
		 double *A1, double *A2, double *A3, double *A4,
		 double *B0, double *B1, double *B2, double *C0)
{
  const double x2 = x * x ;

  (*A1) += x ;
  (*A2) += x2 ;
  (*A3) += x * x2 ;
  (*A4) += x2 * x2 ;
  (*B0) += y ;
  (*B1) += y * x ;
  (*B2) += y * x2 ;
  (*C0) += y * y ;
}


// Calculates a second order regression (parabola:  a*x^2 + b*x + c) fit to given data (time, amp)
// around a given instant (sampTime) going up and down by a certain number of data points
// to minimize round off errors, time[] is centered around sampTime and the results are corrected
// for the shift at the end.
// Assumes that the data are sorted in time
// Output: pointers to: a b c  and the variance around the fit (pVariance, if not null)
void regressionOrder2 (float sampTime, float *time, float *amp, int size, double *pa, double *pb, double *pc, double *pVariance)
{
  int i, indexDown = 0, indexUp = size - 1 ;
  int A0 = NUM_POINTS_HALF_KERNEL + NUM_POINTS_HALF_KERNEL ;
  double A1 = 0.0, A2 = 0.0, A3 = 0.0, A4 = 0.0 ;
  double B0 = 0.0, B1 = 0.0, B2 = 0.0, C0 = 0.0 ;
  double denom, a, b, c ;

  // Step 1: find the closes index above it
  // Step 1.1: run the bisection algorithm
  while ((indexDown + 1) < indexUp)
    {
      i = (indexDown + indexUp) / 2 ;

      if (time[i] > sampTime)
	indexUp = i ;
      else
	indexDown = i ;
    }

  // Step 1.2: take care of the ends
  if (time[indexDown] > sampTime)
    indexUp = 0 ;
  else
    if (time[indexUp] <= sampTime)
      indexDown = size - 1 ;


  // Step 2: do statistics
  if (indexUp < NUM_POINTS_HALF_KERNEL)
    {
      for (i = indexUp + size - NUM_POINTS_HALF_KERNEL ; i < size ; i++)
	updateVars (time[i] - 1.0 - sampTime, amp[i], &A1, &A2, &A3, &A4, &B0, &B1, &B2, &C0) ;

      for (i = 0 ; i < indexUp + NUM_POINTS_HALF_KERNEL ; i++)
	updateVars (time[i] - sampTime, amp[i], &A1, &A2, &A3, &A4, &B0, &B1, &B2, &C0) ;
    }
  else if (indexDown >= (size - NUM_POINTS_HALF_KERNEL))
    {
      for (i = indexDown + 1 - NUM_POINTS_HALF_KERNEL ; i < size ; i++)
	updateVars (time[i] - sampTime, amp[i], &A1, &A2, &A3, &A4, &B0, &B1, &B2, &C0) ;

      for (i = 0 ; i <= (indexDown + NUM_POINTS_HALF_KERNEL - size) ; i++)
	updateVars (time[i] + 1.0 - sampTime, amp[i], &A1, &A2, &A3, &A4, &B0, &B1, &B2, &C0) ;
    }
  else  // Normal case
    {
      for (i = indexUp - NUM_POINTS_HALF_KERNEL ; i < (indexUp + NUM_POINTS_HALF_KERNEL) ; i++)
	updateVars (time[i] - sampTime, amp[i], &A1, &A2, &A3, &A4, &B0, &B1, &B2, &C0) ;
    }

  denom = (A2 * ((A0*A4) + (2.0*A1*A3) - (A2*A2))) - (A0*A3*A3) - (A1*A1*A4) ;  // Denominator

  if (denom == 0.0)
    denom = EPSILON ;    // Prevents division by zero

  *pa = ((A0 * ((A2*B2) - (A3*B1))) + (A1 * ((A3*B0) - (A1*B2))) + (A2 * ((A1*B1) - (A2*B0)))) / denom ;
  *pb = ((A0 * ((A4*B1) - (A3*B2))) + (A1 * ((A2*B2) - (A4*B0))) + (A2 * ((A3*B0) - (A2*B1)))) / denom ;
  *pc = ((A1 * ((A3*B2) - (A4*B1))) + (A2 * ((A4*B0) - (A2*B2))) + (A3 * ((A2*B1) - (A3*B0)))) / denom ;

  if (pVariance)
    {
      a = *pa ;
      b = *pb ;
      c = *pc ;
      *pVariance = ((A4*a*a) + (2.0*a*b*A3) + (A2 * ((b*b) + (2.0*a*c))) + (2.0*b*c*A1) + (A0*c*c) - (2.0 * ((B2*a) + (B1*b) + (B0*c))) + C0) / (A0-3) ;
    }

  // Moves back the "zero point":  a(x-x0)^2 + b(x-x0) + c  -->  ax^2 + bx + c
  *pc += (((*pa) * sampTime * sampTime) - ((*pb) * sampTime)) ;
  *pb -= (2.0 * (*pa) * sampTime) ;
}


// Interpolates a 2nd order regression spline (parabola) and returns the amplitude at
// an arbitrary given value (sampTime)
// sampTime - moment to interpolate around (+/- kernelHalfWidth)
// time - input: array of time stamps (may be rearranged)
// amp - input: matching array of amplitudes (may be rearranged)
// size - size of input arrays
// Returns the 2nd order regression at the given point (sampTime)
double interpolateSmooth (double sampTime, float *time, float *amp, int size)
{
  double a, b, c ;

  regressionOrder2 (sampTime, time, amp, size, &a, &b, &c, 0) ;

  return ((a * sampTime * sampTime) + (b * sampTime) + c) ;
}



// Interpolates a 2nd order regression spline (fit to a parabola) and returns the
// time of a given amplitude (goal). It will return the closer (to sampTime) of the
// two solutions, as long as it's within SEARCH_STEP_SIZE
// sampTime - time to fit around
// time - input: array of time stamps (may be rearranged)
// amp - input: matching array of amplitudes (may be rearranged)
// size - size of input arrays
// Returns the nearest x-axis value of the fitted parabola for a given y-axis value
//  if it's out of range (sampTime +/- SEARCH_STEP_SIZE), returns sampTime
double findGoalSmooth (double sampTime, double goal, float *time, float *amp, int size)
{
  double a, b, c, D, x1, x2, x ;

  regressionOrder2 (sampTime, time, amp, size, &a, &b, &c, 0) ;

  D = (b*b) - (4.0 * a * (c - goal)) ;   // The determinant

  if ((D < 0.0) || (a == 0.0))
    return (sampTime) ;

  x1 = (-b + sqrt(D)) / (2.0 * a) ;
  x2 = (-b - sqrt(D)) / (2.0 * a) ;

  if (fabs(x1 - sampTime) < fabs(x2 - sampTime))
    x = x1 ;
  else
    x = x2 ;

  if (fabs(sampTime - x) > SEARCH_STEP_SIZE) // Should be less than half SEARCH_STEP_SIZE
    return (sampTime) ;

  return (mod1(x)) ;
}


// Removes the outliers from the input data array
// by comparing the data against a 2nd order regression spline.
// the data is assumed to be sorted by time (modulo period = phase)
// even when points are removes, the sorted order is maintained
// note that it will keep "lonely points" with no other data points near it
// time - input: array of time stamps
// amp - input: matching array of amplitudes
// size - size of input arrays
// *psize = outputs the new size of the array (will be reduced if outliers are found)
// Returns the chi2 around the spline (MAXFLOAT in case of error)
double ridOutliers (float *time, float *amp, float *err, int *psize)
{
  int index, doRid, numVariance ;
  double x, sumVariance, variance, varianceFit, sumChi2 ;
  double a, b, c ;

  do
    {
      doRid = 0 ;
      index = 0 ;
      numVariance = 0 ;
      sumVariance = 0.0 ;
      sumChi2 = 0.0 ;

      while (index < (*psize))
	{
	  x = time[index] ;

	  regressionOrder2 (x, time, amp, *psize, &a, &b, &c, &varianceFit) ;

	  variance = sqr((a * x * x) + (b * x) + c - amp[index]) ;

	  if (variance > (varianceFit * OUTLIER_VARIANCE_LIMIT))
	    {
	      doRid = 1 ;
	      (*psize)-- ;

	      if (index < (*psize))
		{
		  memmove (&time[index], &time[index+1], ((*psize) - index) * sizeof(float)) ;
		  memmove (&amp[index], &amp[index+1], ((*psize) - index) * sizeof(float)) ;
		  memmove (&err[index], &err[index+1], ((*psize) - index) * sizeof(float)) ;
		}
	    }
	  else
	    {
	      sumVariance += variance ;
	      sumChi2 += variance / sqr(err[index]) ;
	      numVariance++ ;
	      index++ ;
	    }
	}
    }
  while (doRid) ;

  if (numVariance == 0)
    return (MAXFLOAT) ;

  return (sumChi2 / numVariance) ;
}


//======================= Make a "first guess" for the initial parameters ===============================

// Data structure of a single light curve observation
struct LCdata
{
  float primary ;
  float secondary ;
  float tertiary ;
};


// A comparator for a pair-o-float
// Returns: -1 if x.primary < y.primary
//           0 if x.primary = y.primary
//           1 if x.primary > y.primary
int comparLCdata (const void *x, const void *y)
{
  if (((struct LCdata*)x)->primary < ((struct LCdata*)y)->primary)
    return (-1) ;
  else
    return (((struct LCdata*)x)->primary > ((struct LCdata*)y)->primary) ;
}


// Returns 0 = success ; 1 = failure (malloc)
int sortSamples (float *time, float *amp, float *err, int size)
{
  int i ;
  struct LCdata *arr = (struct LCdata*)malloc (size * sizeof(struct LCdata)) ;

  if (!arr)
    return (1) ;

  for (i = 0 ; i < size ; i++)
    {
      arr[i].primary = time[i] ;
      arr[i].secondary = amp[i] ;
      arr[i].tertiary = err[i] ;
    }

  qsort (arr, size, sizeof(struct LCdata), comparLCdata) ;

  for (i = 0 ; i < size ; i++)
    {
      time[i] = arr[i].primary ;
      amp[i] = arr[i].secondary ;
      err[i] = arr[i].tertiary ;
    }

  free (arr) ;
  return (0) ;
}


// A comparator function for sorting an array of floating point values
int comparFloat (const void *x, const void *y)
{
  if ((*((float*)x)) < (*((float*)y)))
    return (-1) ;
  else
    return ((*((float*)x)) > (*((float*)y))) ;
}



// Given a floating point array
// Returns the median amplitude  (or MAXFLOAT in case of error - malloc failure)
double getMedian (float *data, int size)
{
  double median ;
  float *arr = (float*)malloc (size * sizeof(float)) ;

  if (!arr)
    return (MAXFLOAT) ;

  memcpy (arr, data, size * sizeof(float)) ;

  qsort (arr, size, sizeof(float), comparFloat) ;

  if (size & 1)  // is odd
    median = arr[(size-1)/2] ;
  else           // is even
    median = 0.5 * (arr[size/2] + arr[(size/2)-1]) ;

  free (arr) ;
  return (median) ;
}



//----------------------- Find dips ------------------------------------

// Finds the two deepest minima of the light curve
// The primary minimum is simply at the smallest value
// The secondary minimum is the place that has the largest difference
// between the value itself and the maxima between it and the primary minimum
// Note: if a secondary minimum can't be found (very unusual), will set *pmin2Amp = MAXFLOAT
// May return an "edge minimum"- near (*pmin1Time + 0.5)
// *pmax1Amp >= *pmax2Amp >= *pmin1Amp ; *pmax1Amp >= *pmin2Amp >= *pmin1Amp
// ideally:  *pmax2Amp > *pmin2Amp  but there is no guarantee (can be used as a test)
void findDips (float *time, float *amp, int size,
	       double *pmin1Time, double *pmin2Time,
	       double *pmin1Amp, double *pmin2Amp)
{
  double maxAmp, min1TimePlusHalf ;
  double diff, sampTime, sampAmp ;

  // Step 1: search for the primary dip
  *pmin1Time = 0.0 ;
  *pmin1Amp = interpolateSmooth (0.0, time, amp, size) ;

  for (sampTime = SEARCH_STEP_SIZE ; sampTime < 1.0 ; sampTime += SEARCH_STEP_SIZE)
    {
      sampAmp = interpolateSmooth (sampTime, time, amp, size) ;

      if (*pmin1Amp > sampAmp)
	{
	  *pmin1Amp = sampAmp ;
	  *pmin1Time = sampTime ;
	}
    }

  // Step 2: search for the secondary dip min1TimePlusHalf = min1Time + 0.5 ;
  *pmin2Amp = MAXFLOAT ;
  min1TimePlusHalf = *pmin1Time + 0.5 ;
  maxAmp = *pmin1Amp ;
  diff = 0.0 ;

  if (min1TimePlusHalf < 1.0)
    {
      // 2.1: look up
      for (sampTime = *pmin1Time + SEARCH_STEP_SIZE ; sampTime <= min1TimePlusHalf ; sampTime += SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (maxAmp < sampAmp)
	    maxAmp = sampAmp ;

	  if (diff < maxAmp - sampAmp)
	    {
	      diff = maxAmp - sampAmp ;
	      *pmin2Time = sampTime ;
	      *pmin2Amp = sampAmp ;
	    }
	}

      // 2.2: look down
      maxAmp = *pmin1Amp ;

      for (sampTime = *pmin1Time - SEARCH_STEP_SIZE ; sampTime >= 0.0 ; sampTime -= SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (maxAmp < sampAmp)
	    maxAmp = sampAmp ;

	  if (diff < maxAmp - sampAmp)
	    {
	      diff = maxAmp - sampAmp ;
	      *pmin2Time = sampTime ;
	      *pmin2Amp = sampAmp ;
	    }
	}

      for (sampTime = 1.0 - SEARCH_STEP_SIZE ; sampTime >= min1TimePlusHalf ; sampTime -= SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (maxAmp < sampAmp)
	    maxAmp = sampAmp ;

	  if (diff < maxAmp - sampAmp)
	    {
	      diff = maxAmp - sampAmp ;
	      *pmin2Time = sampTime ;
	      *pmin2Amp = sampAmp ;
	    }
	}
    }
  else  // min1TimePlusHalf >= 1.0
    {
      min1TimePlusHalf -= 1.0 ;  // modulo 1.0

      // 2.1: look up
      for (sampTime = *pmin1Time + SEARCH_STEP_SIZE ; sampTime < 1.0 ; sampTime += SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (maxAmp < sampAmp)
	    maxAmp = sampAmp ;

	  if (diff < (maxAmp - sampAmp))
	    {
	      diff = maxAmp - sampAmp ;
	      *pmin2Time = sampTime ;
	      *pmin2Amp = sampAmp ;
	    }
	}

      for (sampTime = 0.0 ; sampTime <= min1TimePlusHalf ; sampTime += SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (maxAmp < sampAmp)
	    maxAmp = sampAmp ;

	  if (diff < maxAmp - sampAmp)
	    {
	      diff = maxAmp - sampAmp ;
	      *pmin2Time = sampTime ;
	      *pmin2Amp = sampAmp ;
	    }
	}

      // 2.2: look down
      maxAmp = *pmin1Amp ;

      for (sampTime = *pmin1Time - SEARCH_STEP_SIZE ; sampTime >= min1TimePlusHalf ; sampTime -= SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (maxAmp < sampAmp)
	    maxAmp = sampAmp ;

	  if (diff < (maxAmp - sampAmp))
	    {
	      diff = maxAmp - sampAmp ;
	      *pmin2Time = sampTime ;
	      *pmin2Amp = sampAmp ;
	    }
	}
    }
}


// Finds the values (should be about max) between the dips (pmax1Amp, pmax2Amp)
// given the average of the two dip times (min1Time, min2Time)
// Note that we can't use a simple max, since it may pick the interpolation
// overshoot at the edges of the dips. This approach is far more robust.
// Sets:  *pmax1Amp >= *pmax2Amp
void findMidMaxs (float *time, float *amp, int size, double avrTime, double *pmax1Amp, double *pmax2Amp)
{
  *pmax1Amp = interpolateSmooth(avrTime, time, amp, size) ;
  *pmax2Amp = interpolateSmooth(mod1(avrTime + 0.5), time, amp, size) ;

  if (*pmax1Amp < *pmax2Amp)
    swap(pmax1Amp, pmax2Amp) ;
}



// Assumed that goalAmp is larger than the amplitude at *pminTime
// Returns the half-width of the given dip, to the point it crosses a given goal amplitude (goalAmp)
// Returns 0.0 upon fatal error (couldn't reach goalAmp)
double findHalfWidth (float *time, float *amp, int size, double goalAmp, double *pminTime)
{
  double minTimePlusHalf = *pminTime + 0.5 ;
  double sampAmp, sampTime ;
  double limitUp = MAXFLOAT, limitDown = MAXFLOAT ;
  double limitUpFineTune, limitDownFineTune ;

  // Step 1: find dip edges
  if (minTimePlusHalf < 1.0)
    {
      // 1.1: look up
      for (sampTime = *pminTime + SEARCH_STEP_SIZE ; (limitUp == MAXFLOAT) && (sampTime <= minTimePlusHalf) ; sampTime += SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (sampAmp > goalAmp)
	    limitUp = sampTime ;
	}

      // 1.2: look down
      for (sampTime = *pminTime - SEARCH_STEP_SIZE ; (limitDown == MAXFLOAT) && (sampTime >= 0.0) ; sampTime -= SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (sampAmp > goalAmp)
	    limitDown = sampTime ;
	}

      for (sampTime = 1.0 - SEARCH_STEP_SIZE ; (limitDown == MAXFLOAT) && (sampTime >= minTimePlusHalf) ; sampTime -= SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (sampAmp > goalAmp)
	    limitDown = sampTime ;
	}
    }
  else  // minTimePlusHalf >= 1.0
    {
      minTimePlusHalf -= 1.0 ;  // modulo 1.0

      // 1.1: look up
      for (sampTime = *pminTime + SEARCH_STEP_SIZE ; (limitUp == MAXFLOAT) && (sampTime < 1.0) ; sampTime += SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (sampAmp > goalAmp)
	    limitUp = sampTime ;
	}

      for (sampTime = 0.0 ; (limitUp == MAXFLOAT) && (sampTime <= minTimePlusHalf) ; sampTime += SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (sampAmp > goalAmp)
	    limitUp = sampTime ;
	}

      // 1.2: look down
      for (sampTime = *pminTime - SEARCH_STEP_SIZE ; (limitDown == MAXFLOAT) && (sampTime >= minTimePlusHalf) ; sampTime -= SEARCH_STEP_SIZE)
	{
	  sampAmp = interpolateSmooth(sampTime, time, amp, size) ;

	  if (sampAmp > goalAmp)
	    limitDown = sampTime ;
	}
    }

  if ((limitUp == MAXFLOAT) || (limitDown == MAXFLOAT))
    return (0.0) ;  // Fatal error - couldn't reach goalAmp


  // Step 2: correct for the search overshoot:  (error <= 0.5)
  limitUp -= (0.5 * SEARCH_STEP_SIZE) ;
  if (limitUp < 0.0) limitUp += 1.0 ;

  limitDown += (0.5 * SEARCH_STEP_SIZE) ;
  if (limitDown >= 1.0) limitDown -= 1.0 ;

  // Step 3: fine tune output results:
  limitUpFineTune = findGoalSmooth (limitUp, goalAmp, time, amp, size) ;
  limitDownFineTune = findGoalSmooth (limitDown, goalAmp, time, amp, size) ;

  // Step 4: update minimum:  (especially important for square dips)
  if (limitDownFineTune < limitUpFineTune)
    {
      *pminTime = 0.5 * (limitUpFineTune + limitDownFineTune) ;

      return (0.5 * (limitUpFineTune - limitDownFineTune)) ;
    }

  *pminTime = 0.5 * (1.0 + limitUpFineTune + limitDownFineTune) ;

  if (*pminTime > 1.0)
    *pminTime -= 1.0 ;

  return (0.5 * (1.0 + limitUpFineTune - limitDownFineTune)) ;
}


// =============== Take dips into account =========================


// Returns 1 if t is in the given dip, otherwise returns 0
int isInDip (double t, double minTime, double halfWidth)
{
  double diff = fabs(t - minTime) ;

  return ((diff < halfWidth) || ((1.0 - diff) < halfWidth)) ;
}


// Makes sure that there is a sufficient amount of data to contain a given dip
// minTime = time of dip's minimum
// halfWidth = the dip's half width
// returns: 1 = OK, plenty of data for dip  ;  0 = bad, not enough data
int isEnoughDipData (float *time, int size, double minTime, double halfWidth)
{
  int i, dipDataNum = 0 ;

  for (i = 0 ; (dipDataNum < MIN_NUM_DATA) && (i < size) ; i++)
    if (isInDip (time[i], minTime, halfWidth))
      dipDataNum++ ;

  return (dipDataNum == MIN_NUM_DATA) ;
}


// This may only correct by a small amount. But it's very important to get the
// depth of the small dip as accurately as possible
// pPlateauStdDiv = (output pointer) the standard deviation of the plateau amplitudes
// pPlateauWaviness = (output pointer) a statistical measure for how wavy the plateau is.
//  The range is [-1,1], where 1 is very wavy and 0 is not wavy at all (negative values are unphysical).
// note that this is better than stddiv around a spline because a spline,
// has the nasty effect overshooting around the dips even worse, the spline stddiv is strongly
// determined by NUM_POINTS_HALF_KERNEL and might trace the deviations too closely (or not close enough)
// Returns: 0 = success ; 1 = malloc failure ; 2 = goal couldn't be reached ; 3 = not enough data in plateau ; 4 = zero variance
int findPlateauMedian (float *time, float *amp, int size, double min1Time, double min2Time,
		       double *pMedianAmp, double *pPlateauStdDiv, double *pPlateauWaviness)
{
  int i, reducedSize = 0 ;
  double min1HalfWidth, min2HalfWidth ;
  double sumVarianceAmp, sumNeighborVariance ;
  float *reducedAmp = (float*)malloc (size * sizeof(float)) ;

  if (!reducedAmp)
    return (1) ;

  min1HalfWidth = findHalfWidth(time, amp, size, *pMedianAmp, &min1Time) ;
  min2HalfWidth = findHalfWidth(time, amp, size, *pMedianAmp, &min2Time) ;

  if ((min1HalfWidth == 0.0) || (min2HalfWidth == 0.0))
    {
      free (reducedAmp) ;
      return (2) ;
    }

  for (i = 0 ; i < size ; i++)
    if (!isInDip (time[i], min1Time, min1HalfWidth) && !isInDip(time[i], min2Time, min2HalfWidth))
      {
	reducedAmp[reducedSize] = amp[i] ;
	reducedSize++ ;
      }

  // Check the number of data points in the plateau
  if (reducedSize < MIN_NUM_DATA)
    {
      free (reducedAmp) ;
      return (3) ;
    }

  *pMedianAmp = getMedian(reducedAmp, size) ;

  // Note that I'm calculating around the median, instead of the average, so to make it more robust
  // doing this will make the variance (or sumVarianceAmp) larger than it would have been otherwise.
  // but usually the average will be very close to the median, so the difference will be small.
  // Also note that I'm calculating the variance in this manner to reduce numerical errors, which
  // could become a serious problem.
  // Waviness measure - from the ratio of the mean square difference between neighbor to the variance
  sumVarianceAmp = sqr(reducedAmp[0] - (*pMedianAmp)) ;
  sumNeighborVariance = sqr(reducedAmp[size-1] - reducedAmp[0]) ;
  for (i = 1 ; i < reducedSize ; i++)
    {
      sumVarianceAmp += sqr(reducedAmp[i] - (*pMedianAmp)) ;
      sumNeighborVariance += sqr(reducedAmp[i-1] - reducedAmp[i]) ;
    }

  if (sumVarianceAmp <= 0.0)
    return (4) ;  // This should never happen

  *pPlateauStdDiv = sqrt(sumVarianceAmp / reducedSize) ;
  *pPlateauWaviness = 1.0 - (0.5 * sumNeighborVariance / sumVarianceAmp) ;

  free (reducedAmp) ;
  return (0) ;
}

//-------------------------------------------------------


// returns the sum of:
//           1  - period OK
//           2  - double period OK
//           4  - half period OK
int NumEclipses (double period, FILE *finLC, float *time, float *amp, float *err, int size)
{
  int i, res = 7 ;
  double tmpTime, tmpMag, tmpErr ;
  double min1Time, min2Time, min1Amp, min2Amp, max1Amp, max2Amp ;
  double medianAmp, min1HalfWidth, min2HalfWidth ;
  double goalAmp, plateauStdDiv, plateauWaviness ;
  double sigmaMin1Amp, sigmaMin2Amp, sigmaMax1Amp ; // sigmaMaxDiff ;

  rewind(finLC) ;
  i = 0 ;

  while ((i < size) && (3 == fscanf(finLC, "%lf %lf %lf\n", &tmpTime, &tmpMag, &tmpErr)))
    if (tmpErr > 0.0)  // Sign of an invalid magnitude
      {
	time[i] = (float)mod1(tmpTime / period) ;
	amp[i] = (float)pow(10.0, -0.4 * tmpMag) ;  // Convert magnitude (logarithmic) to amplitude (linear)
	err[i] = (float)(amp[i] * tmpErr / MAGNITUDE_ERROR_FACTOR) ;  // Convert to absolute error
	i++ ;
      }

  if (i != size) 
    {
      printf ("ERROR: Number mismatch\n") ;
      return (0) ;
    }

  if (sortSamples(time, amp, err, size))
    {
      printf ("ERROR: sortSamples() malloc failure\n") ;
      return (0) ;
    }

  ridOutliers(time, amp, err, &size) ;

  if (size < (3 * MIN_NUM_DATA))  // Need MIN_NUM_DATA for both the dips and the plateau
    {
     printf ("Warning: Not enough non-outlier data points\n") ;
     return (0) ;
    }

  //-------------------------------------------------

  // Find median:
  medianAmp = getMedian(amp, size) ;
  if (medianAmp == MAXFLOAT)
    {
      printf ("ERROR: getMedian() malloc failure\n") ;
      return (0) ;
    }

  //--------------------------------------------------
  // Get the location of the dips:

  // Coarse search (secondary may be an "edge minimum")
  findDips(time, amp, size, &min1Time, &min2Time, &min1Amp, &min2Amp) ;

  // Fine-tune dips locations and finds their widths
  goalAmp = 0.5 * (medianAmp + min1Amp) ;  // Half width half min (coarse)
  (void)findHalfWidth(time, amp, size, goalAmp, &min1Time) ;
  min1Amp = interpolateSmooth (min1Time, time, amp, size) ;
  goalAmp = 0.5 * (medianAmp + min2Amp) ;  // Half width half min (coarse)
  (void)findHalfWidth(time, amp, size, goalAmp, &min2Time) ;
  min2Amp = interpolateSmooth (min2Time, time, amp, size) ;

  // find "maxs" midway between the dips
  findMidMaxs(time, amp, size, 0.5 * (min1Time + min2Time), &max1Amp, &max2Amp) ;

  // fine-tune median:
  i = findPlateauMedian(time, amp, size, min1Time, min2Time,
			&medianAmp, &plateauStdDiv, &plateauWaviness) ;
  if (i)
    {
      printf ("Error: #%d in findPlateauMedian()\n", i) ;
      return (0) ;
    }

  goalAmp = 0.5 * (medianAmp + min1Amp) ;  // half width half min (fine tune)
  min1HalfWidth = findHalfWidth(time, amp, size, goalAmp, &min1Time) ;
  min1Amp = interpolateSmooth (min1Time, time, amp, size) ;
  goalAmp = 0.5 * (medianAmp + min2Amp) ;  // half width half min (fine tune)
  min2HalfWidth = findHalfWidth(time, amp, size, goalAmp, &min2Time) ;
  min2Amp = interpolateSmooth (min2Time, time, amp, size) ;

  sigmaMin1Amp = (medianAmp - min1Amp) / plateauStdDiv ;
  sigmaMin2Amp = (medianAmp - min2Amp) / plateauStdDiv ;
  sigmaMax1Amp = (max1Amp - medianAmp) / plateauStdDiv ;
  //  sigmaMaxDiff = (max1Amp - max2Amp) / plateauStdDiv ;

    printf ("== %f %f %f %f\n", min1Time, sigmaMin1Amp, min2Time, sigmaMin2Amp) ;

  //------------- Tests: ---------------------
  // Make sure that the dips are okay
  // I don't use err[] here, in case the data is inherently volatile

  if ((min1HalfWidth == 0.0) || (min2HalfWidth == 0.0))
    {
      printf ("ERROR: Couldn't reach goal amplitude\n") ;
      return (0) ;
    }

  if (sigmaMin1Amp < MIN_STDDIV_BELOW_MEDIAN)
    {
      printf ("Warning: Primary dip is too small\n") ;
      return (0) ;
    }

  if (sigmaMax1Amp > MAX_STDDIV_ABOVE_MEDIAN)
    {
      printf ("Warning: Large hump at mid-plateau\n") ;
      res &= 6 ;
    }

  if (sigmaMin2Amp < MIN_STDDIV_BELOW_MEDIAN)
    {
      printf ("Warning: No secondary dip\n") ;
      res &= 3 ;
    }
  else if ((sigmaMin1Amp - sigmaMin2Amp) > MIN_DIP_SIG_DIFF)  // A "good" light curve
    res &= 1 ;


  if ((fabs(min1Time - min2Time) < (min1HalfWidth + min2HalfWidth)) ||
      (fabs(min1Time - min2Time) > (1.0 - min1HalfWidth - min2HalfWidth)))
    {
      printf ("Warning: Dips are overlapping\n") ;
      res &= 6 ;
    }

  if ((min1HalfWidth >= 0.25) || (min2HalfWidth >= 0.25))
    {
      printf ("Warning: Out of range half width\n") ;
      res &= 6 ;
    }

  if (!isEnoughDipData(time, size, min1Time, min1HalfWidth) ||
      !isEnoughDipData(time, size, min2Time, min2HalfWidth))
    {
      printf ("Warning: Not enough data to constrain one of the dips\n") ;
      res &= 6 ;
    }

  if ((min1Amp + min2Amp) < medianAmp)
    {
      printf ("Warning: The dips are too deep\n") ;
      res &= 4 ;
    }
  else if ((min1Amp + min1Amp) < medianAmp)
    {
      printf ("Warning: The primry dip is too deep for doubling\n") ;
      res &= 5 ;
    }
 

  return (res) ;
}


//-------------------------------------------------------

int main (int argc, char **argv)
{
  int size, res1, res2 ;
  float *time, *amp, *err ;
  double period, tmpErr ;
  char filename[1000] ;
  FILE *finLC, *fout, *finList ;

  if (argc != 3)
    {
      printf ("Usage: %s <input LC list filename> <output LC list filename>\n", argv[0]) ;
      return (1) ;
    }

  if (!strcmp(argv[1], argv[2]))
    {
      printf ("ERROR: all the input/output filenames must be different\n") ;
      return (2) ;
    }

  if (!(finList = fopen (argv[1], "rt")))
    {
      printf ("ERROR: couldn't open the input LC file ('%s')\n", argv[1]) ;
      return (3) ;
    }

  if (!(fout = fopen (argv[2], "wt")))
    {
      printf ("ERROR: couldn't open the output data file ('%s')\n", argv[2]) ;
      fclose (finList) ;
      return (4) ;
    }

  //-----------------------------------------------------------------------

  while (2 == fscanf(finList, "%s %lf\n", filename, &period))
    {
      printf ("%s\n", filename) ;

      if (period <= 0.0)
	{
	  printf ("Warning: Invalid period\n") ;
	  continue ;
	}

      if (!(finLC = fopen (filename, "rt")))
	{
	  printf ("Warning: Couldn't open the input LC file\n") ;
	  continue ;
	}

      //----------------------------------------------------------------

      size = 0 ;
      while (1 == fscanf(finLC, "%*f %*f %lf\n", &tmpErr))
	if (tmpErr > 0.0)  // Sign of an invalid magnitude
	  size++ ;

      if (size < (3 * MIN_NUM_DATA))  // Need MIN_NUM_DATA for both the dips and the plateau
	{
	  printf ("Warning: Not enough data points\n") ;
	  fclose (finLC) ;
	  continue ;
	}

      time = (float*)malloc (size * sizeof(float)) ;
      amp = (float*)malloc (size * sizeof(float)) ;
      err = (float*)malloc (size * sizeof(float)) ;

      if ((!time) || (!amp) || (!err))
	{
	  if (time) free(time) ;
	  if (amp) free(amp) ;
	  if (err) free(err) ;

	  printf ("ERROR: Not enough memory\n") ;
	  fclose(finLC) ;
	  continue ;
	}

      //-------------------------

      // returns the sum of:
      //           1  - period OK
      //           2  - double period OK
      //           4  - half period OK
      res1 = NumEclipses (period, finLC, time, amp, err, size) ;
      printf ("norm: res = %d\n", res1) ;
      if (res1 & 1)
	fprintf (fout, "%s  %.12f\n", filename, period) ;

      if (res1 & 2)
	{
	  res2 = NumEclipses (2.0 * period, finLC, time, amp, err, size) ;
	  printf ("double: res = %d\n", res2) ;
	  
	  if (res2 & 1)
	    fprintf (fout, "%s  %.12f\n", filename, 2.0 * period) ;
	}

      if (res1 & 4)
	{
	  res2 = NumEclipses (0.5 * period, finLC, time, amp, err, size) ;
	  printf ("half: res = %d\n", res2) ;

	  if (res2 & 1)
	    fprintf (fout, "%s  %.12f\n", filename, 0.5 * period) ;
	}
      

      fclose (finLC) ;
      free (time) ;
      free (amp) ;
      free (err) ;
    }
  
  fclose (finList) ;
  fclose (fout) ;
 
  return (0) ;
}
