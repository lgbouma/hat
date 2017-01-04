/************************************************
 *
 * Converts format: <time> <mag>
 * Into OGLE format
 *
 * Compile:   gcc makeOgleErr.c -o makeOgleErr -Wall
 *
 * Run:  makeOgleErr lightcurve.txt > lightcurve.dat
 *
 ***********************************************/

#include <stdio.h>
#include <stdlib.h>

int main (int argc, char **argv)
{
  FILE *fin ;
  float time, mag, err ;

  if (argc != 3)
    {
      printf ("usage: %s <input file> <magnitude error>\n", argv[0]) ;
      return (1) ;
    }

  err = atof(argv[2]) ;

  if (err <= 0.0)
    {
      printf ("ERROR: magnitude error must be positive\n") ;
      return (2) ;
    }
  
  fin = fopen (argv[1], "rt") ;

  if (!fin) 
    {
      printf ("ERROR: couldn't open input file '%s'\n", argv[1]) ;
      return (3) ;
    }

  while (2 == fscanf(fin, "%f %f\n", &time, &mag))
    printf ("%f  0.0 0.0 %f %f 0\n", time, mag, err) ;

  fclose (fin) ;
  return (0) ;
}
