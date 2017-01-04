
/******************************************************************************
 *
 * This program creates a Super Mongo (SM) script for ploting all the
 * histograms relevant to the full DEBiL database
 *
 * compile: gcc histAll.c -o histAll -Wall
 *
 * run:
 *
 * > histAll hist_chi3.txt 
 * > sm
 *  : inp histAll.sm
 *  : quit
 *
 * > ps2pdf histAll.ps 
 *
 *******************************************************************************/


#include <stdio.h>
#include <stdlib.h>  // for system()

#define NUM_BOXES_X 2
#define NUM_BOXES_Y 3

#define NUM_PARAMS 25


char paramStr[NUM_PARAMS][64] = {
  "Log(period)",            // 1
  "Log(periodCorrected)",   // 2
  "Eccentricity",           // 3
  "r1",                     // 4
  "r2",                     // 5
  "sin(i)",                 // 6
  "sin(i)",                 // 7 - zoomed
  "Epoch of perihelion",    // 8
  "Argument of perihelion [degrees]", // 9
  "Chi2 of model",          // 10
  "Fitness score of model", // 11
  "Secondary dip depth",    // 12
  "Hump height",            // 13
  "Hump difference",        // 14
  "Waviness",               // 15
  "Scatter score",          // 16
  "Log(mean density [g/cm^3])",  // 17
  "Log(max density [g/cm^3])",   // 18
  "Corrected I1",           // 19
  "Corrected I2",           // 20
  "Corrected I1 & I2",      // 21
  "Color (V-I)",            // 22
  "r2 / r1",                // 23
  "I2 - I1",                // 24
  "r1 + r2"                 // 25
};


// return: 0 = finished ; 1 = continue
int doPage (FILE *fout)
{
  static int paramNum = 0 ;
  int i, j ;
  
  for (j = NUM_BOXES_Y ; j > 0  ; j--)
    for (i = 1 ; i <= NUM_BOXES_X ; i++)
      {
	paramNum++ ;
       
	if (paramNum > NUM_PARAMS) 
	  return (0) ;

	fprintf (fout, "window %d %d %d %d\n", NUM_BOXES_X, NUM_BOXES_Y, i, j) ;
	fprintf (fout, "limits x%d n%d\n", paramNum, paramNum) ;
	fprintf (fout, "histogram x%d n%d\n", paramNum, paramNum) ;
	fprintf (fout, "xlabel %s\n", paramStr[paramNum-1]) ;
	fprintf (fout, "ylabel Log(number)\n") ;
	fprintf (fout, "box\n\n") ;
      }
  
  return (1) ;
}



int main(int argc, char **argv)
{
  int i, j ;
  FILE *fout ;

  if (argc != 2)
    {
      printf ("Usage %s  <DEBiL histogram*>\n", argv[0]) ;
      printf ("  [*The output of running debit2hist]\n") ;
      return (1) ;
    }

  fout = fopen ("histAll.sm", "wt") ;

  if (!fout)
    {
      printf ("ERROR: couldn't create output file 'histAll.sm'\n\n") ;
      return (2) ;
    }

  printf ("writing SM scripts file into 'histAll.sm'\n\n") ;

  fprintf (fout, "device postportfile %s.ps\n\n", argv[1]) ;
  fprintf (fout, "expand 0.4\n") ;
  fprintf (fout, "ptype 1 1\n\n") ;
  fprintf (fout, "ctype black\n") ;
  fprintf (fout, "data %s\n", argv[1]) ;
  
  fprintf (fout, "read { x1 1 n1 2 }\n") ;
  fprintf (fout, "read { x2 1 n2 3 }\n") ;

  for (i = 2 ; i < NUM_PARAMS ;)
    {
      fprintf (fout, "read {") ;
      
      for (j = 0 ; (j < 5) && (i < NUM_PARAMS) ; j++, i++)
	fprintf (fout," x%d %d n%d %d ", i + 1, 2 * i, i + 1, (2 * i) + 1) ;

      fprintf (fout,"}\n") ;
    }

  fprintf (fout, "\n") ;

  while (doPage (fout)) 
    fprintf (fout, "page\n") ;

  fprintf (fout, "hardcopy\n") ;
  printf ("done.\n") ;

  fclose (fout) ;
  return (0) ;
}
