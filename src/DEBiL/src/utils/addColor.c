/************************************************************
 *
 * This utility adds color information from the Sumi (2004)
 * extinction map to the DEBiL database
 *
 * Compile:  gcc addColor.c -o addColor -Wall -lm
 *
 * Run:      addColor ag.out agc.out /home/jdevor/work/debil/sumy-color/
 *           sort -n agc.out >! ag.dbl
 *           rm agc.out
 *
 ************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <values.h>


#define MAX_CHI2 4.0

#define ERROR_VALUE -1000.0
#define ERROR_STRING "???"

#define MIN_VALID_MAG_VALUE 9.99   // Larger than 9.999
#define MAX_VALID_MAG_VALUE 50.0   // Smaller than 99.999

#define MAX_VALID_MAG_DIFF 9.99    // Smaller than 9.999
#define MAX_VALID_MAG_ERR 9.99     // Smaller than 9.999

#define FILE_STRING_LEN 1000
#define SMALL_STRING_LEN 100
#define NUM_DIRECTORY_LEVELS 5
#define NAME_OFFSET 6


// Combines the flux of two given magnitudes
double combineMag (double m1, double m2)
{
  return (-2.5 * log10 (pow(10.0, -0.4 * m1) + pow(10.0, -0.4 * m2))) ;
}


// Returns the variable catalogue index of the light curve
// Returns (-1) upon an error
int parseFilename (char *str, int *numGroup)
{
  int i, j, len, countSlash = 0 ;
  char numStr[SMALL_STRING_LEN] ;

  len = strlen (str) ;
  for (i = 1 ; (i < len) && (countSlash < NUM_DIRECTORY_LEVELS) ; i++)
    if (str[i] == '/')
      countSlash++ ;

  i += NAME_OFFSET ;

  if (i >= len)
    return (-1) ;

  for (j = 0 ; (i < len) && (j < SMALL_STRING_LEN) && (str[i] >= '0') && (str[i] <= '9') ; i++, j++)
    numStr[j] = str[i] ;

  i++ ;

  if ((j >= SMALL_STRING_LEN) || (i >= len))
    return (-1) ;

  numStr[j] = 0 ;
  (*numGroup) = atoi(numStr) ;

  for (j = 0 ; (i < len) && (j < SMALL_STRING_LEN) && (str[i] >= '0') && (str[i] <= '9') ; i++, j++)
    numStr[j] = str[i] ;

  if ((j >= SMALL_STRING_LEN) || (i >= len))
    return (-1) ;

  numStr[j] = 0 ;
  return (atoi(numStr)) ;
}



int main(int argc, char **argv)
{
  FILE *fin, *fout, *fcat = 0 ;
  double x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 ;
  double x17, x18, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31 ;
  double chiScore, catI, catVI, catIc, catVIc, catIe, catVe, AI, x9c, x11c, combc ;
  int num, numGroup, prevNumGroup = -1, i19, i20, index, count = 0 ;
  char str[FILE_STRING_LEN], catFilename[FILE_STRING_LEN] ;
  char catRA[SMALL_STRING_LEN], catDec[SMALL_STRING_LEN] ;

  if (argc != 4)
    {
      printf (" Usage:  %s  <input DEBiL database> <output DEBiL with color> <Sumy color info path>\n", argv[0]) ;
      printf (" [note: path must have '/' in the end]\n") ;
      return (1) ;
    }

  fin = fopen (argv[1], "rt") ;

  if (!fin)
    {
      printf ("ERROR: Couldn't open the input file ('%s')\n", argv[1]) ;
      return (3) ;
    }

  fout = fopen (argv[2], "wt") ;

  if (!fout)
    {
      printf ("ERROR: Couldn't open the output file ('%s')\n", argv[2]) ;
      fclose (fin) ;
      return (4) ;
    }
  
  //=========================================================================

  while (31 == fscanf (fin, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
		       str, &x2, &x3, &x4, &x5, &x6, &x7, &x8, &x9, &x10, &x11, &x12, &x13, &x14, &x15, &x16, 
		       &x17, &x18, &i19, &i20, &x21, &x22, &x23, &x24, &x25, &x26, &x27, &x28, &x29, &x30, &x31))
      {
	count++ ;

	num = parseFilename (str, &numGroup) ;

	if (num == -1)
	  {
	    printf ("ERROR: Invalid filename ('%s')\n", str) ;
	    break ;
	  }

	if (numGroup != prevNumGroup)
	  {
	    prevNumGroup = numGroup ;

	    if (fcat)
	      fclose (fcat) ;

	    sprintf (catFilename, "%ssc%d.vi", argv[3], numGroup) ;
	  
	    fcat = fopen (catFilename, "rt") ;

	    index = -1 ;
	    printf ("opening new group\n") ;

	    if (!fcat)
	      {
		printf ("ERROR: Couldn't open catalogue ('%s')\n", catFilename) ;
		break ;
	      } 

	    fscanf (fcat, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n") ;
	  }

	//---------------

	if (index >= num)
	  {
	    rewind (fcat) ;
	    fscanf (fcat, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n") ;
	    printf ("rewind\n") ;
	  }

	while ((9 == fscanf (fcat, "%d %*f %*f %*f %*f %*d %*f %*d %s %s %*f %*f %*f %lf %lf %*d %*d %lf %*d %*d %lf %lf %lf\n", &index, catRA, catDec, &catVI, &catI, &catVe, &catIe, &catIc, &catVIc)) && (index < num)) ;

	if (num != index)
	  {
	    printf ("Warning: Couldn't find the given index in the catalogue (%d in '%s')\n", num, catFilename) ;
	
	    x9c = ERROR_VALUE ;
	    x11c = ERROR_VALUE ;

	    catVI = ERROR_VALUE ;
	    catI = ERROR_VALUE ;
	    combc = ERROR_VALUE ;
	    catIc = ERROR_VALUE ;
	    catVIc = ERROR_VALUE ;
	    catVe = ERROR_VALUE ;
	    catIe = ERROR_VALUE ;

	    sprintf (catRA, ERROR_STRING) ;
	    sprintf (catDec, ERROR_STRING) ;
	  }
	else
	  {
	    printf ("%d %d %d %s\n", count, numGroup, num, str) ;


	    // Standardize the error values
	    if ((catI >= MIN_VALID_MAG_VALUE) && (catI <= MAX_VALID_MAG_VALUE) &&
		(catIc >= MIN_VALID_MAG_VALUE) && (catIc <= MAX_VALID_MAG_VALUE))
	      {
		AI = catI - catIc ;  //  Extinction coefficient (absorption in I-band)
		x9c = x9 - AI ;
		x11c = x11 - AI ;
		combc = combineMag(x9, x11) - AI ;
	      }
	    else
	      {
		x9c = ERROR_VALUE ;
		x11c = ERROR_VALUE ;
		combc = ERROR_VALUE ;
	      }

	    if ((catVI < MIN_VALID_MAG_VALUE) || (catVI > MAX_VALID_MAG_VALUE))
	      catVI = ERROR_VALUE ;

	    if ((catIc < MIN_VALID_MAG_VALUE) || (catIc > MAX_VALID_MAG_VALUE))
	      catIc = ERROR_VALUE ;
	
	    if (fabs(catVIc) > MAX_VALID_MAG_DIFF)
	      catVIc = ERROR_VALUE ;

	    if (fabs(catVe) > MAX_VALID_MAG_ERR)
	      catVe = ERROR_VALUE ;

	    if (fabs(catIe) > MAX_VALID_MAG_ERR)
	      catIe = ERROR_VALUE ;
	  }

	//---------------------

	if (x22 <= x23)
	  chiScore = ERROR_VALUE ;
	else
	  chiScore = (x21 - x22) / (x23 - x22) ;

	//---------------------

	fprintf (fout, "%02d %05d  %.12f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f %f  %f %f %f %f %f %f %f %f %f %s %s\n",
		 numGroup, num, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12,
		 x13, x14, x15, x16, x17, x18, i19, i20, x21, x22, x23, x24, chiScore,
		 x25, x26, x27, x28, x29, x30, x31, x9c, x11c, combc, catIc, catVIc,
		 catI, catVI, catVe, catIe, catRA, catDec) ;
      }
		       
  if (fcat)
    fclose (fcat) ;

  fclose (fin) ;
  fclose (fout) ;
  return (0) ;
}
