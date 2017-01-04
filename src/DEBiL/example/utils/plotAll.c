
/******************************************************************************
 *
 * This program creates a Super Mongo (SM) script for ploting everything
 * versus everything in the full DEBiL database
 *
 * compile: gcc plotAll.c -o plotAll -Wall
 *
 * run:
 *
 * > plotAll ag.filter90.out
 * > sm
 *  : inp plotAll.sm
 *  : quit
 *
 * > ps2pdf plotAll.ps 
 *
 *******************************************************************************/


#include <stdio.h>
#include <stdlib.h>  // for system()

#define NUM_BOXES_X 2
#define NUM_BOXES_Y 3

#define MAX_PARAMS 48           // including added parameters
#define DATABASE_NUM_PARAM 42   // not including RA and Dec
#define NUM_SWAP_PARAM 7        // num paramters that need axis swap
#define NUM_PARAMS 19           // size of subset of parameters used in plots

int paramSwapArr[NUM_SWAP_PARAM] = {10, 12, 34, 35, 36, 37, 39} ;

int paramArr[NUM_PARAMS] = {45, 38, 4, 6, 8, 22, 26, 28, 29, 30, 31, 34, 35, 36, 43, 44, 46, 47, 48} ;


char paramStr[MAX_PARAMS][64] = {
  "Field number",           // 1
  "Light curve number",     // 2
  "Period",                 // 3
  "Eccentricity",           // 4
  "Error in eccentricity",  // 5
  "r1",                     // 6
  "Error in r1",            // 7
  "r2",                     // 8
  "Error in r2",            // 9
  "I1",                     // 10
  "Error in I1",            // 11 
  "I2",                     // 12
  "Error in I2",            // 13
  "sin(i)",                 // 14
  "Error in sin(i)",        // 15
  "Epoch of perihelion",    // 16
  "Error in epoch of perihelion",     // 17
  "Argument of perihelion [degrees]", // 18
  "Error in argument of perihelion",  // 19
  "Num data points",        // 20
  "Num outliers",           // 21
  "Chi2 of model",          // 22
  "Chi2 of median",         // 23
  "Chi2 of spline",         // 24
  "Chi2 of sinusoidal",     // 25
  "Fitness score",          // 26
  "Secondary dip depth",    // 27
  "Hump height",            // 28
  "Hump difference",        // 29
  "Waviness",               // 30
  "Scatter score",          // 31
  "Mean density [g/cm^3]",  // 32
  "Max density [g/cm^3]",   // 33
  "Corrected I1",           // 34
  "Corrected I2",           // 35
  "Corrected I1 & I2",      // 36
  "I  (Udalski's catalog)", // 37
  "Color (V-I)",            // 38
  "Uncorrected I",          // 39
  "Uncorrected color (V-I)",// 40
  "Error in V",             // 41
  "Error in I",             // 42
  "R2 / R1",                // 43 - added
  "I2 - I1",                // 44
  "log(period)",            // 45
  "log(mean density)",      // 46
  "log(max density)",       // 47
  "r1 + r2"                 // 48
};


int doSwap (paramNum)
{
  int i ;
  
  for (i = 0 ; i < NUM_SWAP_PARAM ; i++)
    if (paramNum == paramSwapArr[i])
      return (1) ;

  return (0) ;
}



// return: 0 = finished ; 1 = continue
int doPage (FILE *fout)
{
  static int paramX = 0, paramY = 0 ;
  int i, j ;
  
  for (j = NUM_BOXES_Y ; j > 0  ; j--)
    for (i = 1 ; i <= NUM_BOXES_X ; i++)
      {
	paramY++ ;

	if (paramY == NUM_PARAMS)
	  {
	    paramX++ ;
	    paramY = paramX + 1 ;
	  }

	if (paramY == NUM_PARAMS)
	  return(0) ;
       
	fprintf (fout, "window %d %d %d %d\n", NUM_BOXES_X, NUM_BOXES_Y, i, j) ;
	fprintf (fout, "limits x%d x%d\n", paramArr[paramX], paramArr[paramY]) ;

	if (doSwap (paramArr[paramX]))
	  fprintf (fout, "limits $fx2 $fx1 $fy1 $fy2\n") ; // invert axis direction
	
	if (doSwap (paramArr[paramY]))
	  fprintf (fout, "limits $fx1 $fx2 $fy2 $fy1\n") ; // invert axis direction

	fprintf (fout, "points x%d x%d\n", paramArr[paramX], paramArr[paramY]) ;
	fprintf (fout, "xlabel %s\n", paramStr[paramArr[paramX]-1]) ;
	fprintf (fout, "ylabel %s\n", paramStr[paramArr[paramY]-1]) ;
	fprintf (fout, "box\n\n") ;
      }
  
  return(1) ;
}



int main(int argc, char **argv)
{
  int i, j ;
  FILE *fout ;

  if (argc != 2)
    {
      printf ("Usage %s  <DEBiL database /w color>\n", argv[0]) ;
      return (1) ;
    }

  fout = fopen ("plotAll.sm", "wt") ;

  if (!fout)
    {
      printf ("ERROR: couldn't create output file 'plotAll.sm'\n\n") ;
      return (2) ;
    }

  printf ("writing SM scripts file into 'plotAll.sm'\n\n") ;

  fprintf (fout, "device postportfile %s.ps\n\n", argv[1]) ;
  fprintf (fout, "expand 0.8\n") ;
  fprintf (fout, "ptype 1 1\n\n") ;
  fprintf (fout, "ctype black\n") ;
  fprintf (fout, "data %s\n", argv[1]) ;
  
  for (i = 0 ; i < DATABASE_NUM_PARAM ;)
    {
      fprintf (fout, "read {") ;
      
      for (j = 0 ; (j < 10) && (i < DATABASE_NUM_PARAM) ; j++, i++)
	fprintf (fout," x%d %d", i+1, i+1) ;

      fprintf (fout,"}\n") ;
    }

  // added:
  fprintf (fout, "set x43 = x8 / x6\n") ;
  fprintf (fout, "set x44 = x12 - x10\n") ;
  fprintf (fout, "set x45 = LG(x3)\n") ;
  fprintf (fout, "set x46 = LG(x32)\n") ;
  fprintf (fout, "set x47 = LG(x33)\n") ;
  fprintf (fout, "set x48 = x6 + x8\n") ;
  fprintf (fout, "\n\n") ;

  while (doPage (fout)) 
    fprintf (fout, "page\n") ;

  fprintf (fout, "hardcopy\n") ;
  printf ("done.\n") ;

  fclose (fout) ;
  return (0) ;
}
