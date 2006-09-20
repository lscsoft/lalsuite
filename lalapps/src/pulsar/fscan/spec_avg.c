#include <lal/SFTfileIO.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/*---------- empty initializers ---------- */
LALStatus empty_status;
SFTConstraints empty_constraints;
/*---------- Global variables ----------*/

INT4 lalDebugLevel = 1;

int main(int argc, char **argv)
{
  FILE *fp,*fp2;
  LALStatus status = empty_status;

  SFTCatalog *catalog = NULL;
  SFTVector *sft_vect = NULL;
  INT4 i,j,k;
  UINT4 numBins, nSFT;
  double f_min,f_max; 
  SFTConstraints constraints=empty_constraints;
  LIGOTimeGPS startTime, endTime; 
  double avg;
  char outfile[128],outfile2[128];
  if (argc!=7)
  {
   fprintf(stderr, "startGPS endGPS Detector F_min f_max where?\n"); 
   return(0);
  }  
  
  startTime.gpsSeconds = atoi(argv[1]);
  startTime.gpsNanoSeconds = 0;
  constraints.startTime = &startTime; 
  
  endTime.gpsSeconds = atoi(argv[2]);
  endTime.gpsNanoSeconds = 0;
  constraints.endTime = &endTime;
  constraints.detector = argv[3];
  f_min = atof(argv[4]);
  f_max = atof(argv[5]);
  LALSFTdataFind ( &status, &catalog,argv[6], &constraints );
  LALLoadSFTs ( &status, &sft_vect, catalog, f_min,f_max);
  LALDestroySFTCatalog( &status, &catalog);
  
   numBins = sft_vect->data->data->length;
   nSFT = sft_vect->length;

  
   
   fprintf(stderr, "nSFT = %d\tnumBins = %d\tf0 = %f\n", nSFT, numBins,sft_vect->data->f0);
  sprintf(outfile, "spec_%.2f_%.2f_%s_%d_%d", f_min,f_max,constraints.detector,startTime.gpsSeconds,endTime.gpsSeconds);
  sprintf(outfile2, "spec_%.2f_%.2f_%s_%d_%d_timestamps", f_min,f_max,constraints.detector,startTime.gpsSeconds, endTime.gpsSeconds);
 
  fp = fopen(outfile, "w");
  fp2 = fopen(outfile2, "w");

  for (j=0;j<nSFT;j++)
  { 
 //     fprintf(stderr, "sft %d  out of %d\n", j, nSFT);
	fprintf(fp2, "%d.\t%d\n", j, sft_vect->data[j].epoch.gpsSeconds);

    for ( i=0; i < (numBins-2); i+=180)
    {
      avg = 0.0;
      if (i+180>numBins) {printf("Error\n");return(2);}
      for (k=0;k<180;k++)
        avg += sqrt(sft_vect->data[j].data->data[i+k].re*sft_vect->data[j].data->data[i+k].re + 
                       sft_vect->data[j].data->data[i+k].im*sft_vect->data[j].data->data[i+k].im);
      fprintf(fp,"%e\t",avg/180.0); 
    }
    fprintf(fp,"\n");
  }
  
 LALDestroySFTVector (&status, &sft_vect );
 fclose(fp);
 fclose(fp2);

return(0);

} //main
