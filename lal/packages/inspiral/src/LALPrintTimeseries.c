# include <stdio.h>
  void LALPrintTimeseries (int n, double *signal, double delta, double t0) ;
  void LALPrintTimeseries (int n, double *signal, double delta, double t0) 
{
  int i=0;
  FILE *outfile;
  outfile=fopen("chirp.dat","w");

  do 
     fprintf (outfile, "%e %e\n", i*delta+t0, *(signal+i));

  while (n-++i); 

  printf("&\n");
}
