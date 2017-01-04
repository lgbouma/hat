/************************************************
 *
 * Converts format: <time> <mag> <err>
 * Into OGLE format
 *
 * Compile:   gcc makeOgle.c -o makeOgle -Wall
 *
 * Run:  makeOgle lightcurve.txt > lightcurve.dat
 *
 ***********************************************/

#include <stdio.h>

int main (int argc, char **argv)
{
  FILE *fin ;
  float time, mag, err ;

  if (argc != 2)
    {
      printf ("usage: %s <input file>\n", argv[0]) ;
      return (1) ;
    }

  fin = fopen (argv[1], "rt") ;

  if (!fin) 
    {
      printf ("ERROR: couldn't open input file '%s'\n", argv[1]) ;
      return (2) ;
    }

  while (3 == fscanf(fin, "%f %f %f\n", &time, &mag, &err))
    printf ("%f  0.0 0.0 %f %f 0\n", time, mag, err) ;

  fclose (fin) ;
  return (0) ;
}
