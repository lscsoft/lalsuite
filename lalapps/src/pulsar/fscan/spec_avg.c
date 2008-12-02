/*
*  Copyright (C) 2007 Gregory Mendell
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <lalapps.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/UserInput.h>
#include <lal/SFTfileIO.h>
#include <lal/NormalizeSFTRngMed.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

RCSID("$Id");

int main(int argc, char **argv)
{
  FILE *fp  = NULL;
  FILE *fp2 = NULL;
  FILE *fp3 = NULL;
  LALStatus status = blank_status;

  SFTCatalog *catalog = NULL;
  SFTVector *sft_vect = NULL;
  INT4 i,j,k;
  UINT4 numBins, nSFT;
  SFTConstraints constraints=empty_SFTConstraints;
  LIGOTimeGPS startTime, endTime; 
  double avg;
  REAL4 *timeavg;
  REAL8 f;
  CHAR outbase[128],outfile[128],outfile2[128],outfile3[128];

  BOOLEAN help = 0;
  CHAR *SFTpatt = NULL;
  CHAR *IFO = NULL;
  INT4 startGPS = 0;
  INT4 endGPS = 0;
  REAL8 f_min = 0.0;
  REAL8 f_max = 0.0;
  INT4 blocksRngMed = 101;
  CHAR *outputBname = NULL;

  lalDebugLevel = 0;
  LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);

  LAL_CALL(LALRegisterBOOLUserVar  (&status, "help",         'h', UVAR_HELP,     "Print this help message",     &help        ), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "SFTs",         'I', UVAR_REQUIRED, "SFT location/pattern",        &SFTpatt     ), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "IFO",          'I', UVAR_REQUIRED, "Detector",                    &IFO         ), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "startGPS",     's', UVAR_REQUIRED, "Starting GPS time",           &startGPS    ), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "endGPS",       'e', UVAR_REQUIRED, "Ending GPS time",             &endGPS      ), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "fMin",         'f', UVAR_REQUIRED, "Minimum frequency",           &f_min       ), &status);
  LAL_CALL(LALRegisterREALUserVar  (&status, "fMax",         'F', UVAR_REQUIRED, "Maximum frequency",           &f_max       ), &status);
  LAL_CALL(LALRegisterINTUserVar   (&status, "blocksRngMed", 'w', UVAR_OPTIONAL, "Running Median window size",  &blocksRngMed), &status);
  LAL_CALL(LALRegisterSTRINGUserVar(&status, "outputBname",  'o', UVAR_OPTIONAL, "Base name of output files",   &outputBname ), &status);

  LAL_CALL(LALUserVarReadAllInput(&status, argc, argv), &status);
  if (help)
    return(0);

  startTime.gpsSeconds = startGPS;
  startTime.gpsNanoSeconds = 0;
  constraints.startTime = &startTime; 
  
  endTime.gpsSeconds = endGPS;
  endTime.gpsNanoSeconds = 0;
  constraints.endTime = &endTime;
  constraints.detector = IFO;
  LALSFTdataFind ( &status, &catalog,SFTpatt, &constraints );
  LALLoadSFTs ( &status, &sft_vect, catalog, f_min,f_max);
  LALDestroySFTCatalog( &status, &catalog);

  numBins = sft_vect->data->data->length;
  nSFT = sft_vect->length;

  fprintf(stderr, "nSFT = %d\tnumBins = %d\tf0 = %f\n", nSFT, numBins,sft_vect->data->f0);
  if (LALUserVarWasSet(&outputBname))
    strcpy(outbase, outputBname);
  else
    sprintf(outbase, "spec_%.2f_%.2f_%s_%d_%d", f_min,f_max,constraints.detector,startTime.gpsSeconds,endTime.gpsSeconds);
  sprintf(outfile,  "%s", outbase);
  sprintf(outfile2, "%s_timestamps", outbase);
  sprintf(outfile3, "%s_timeaverage", outbase);
 
  fp = fopen(outfile, "w");
  fp2 = fopen(outfile2, "w");
  fp3 = fopen(outfile3, "w");

  for (j=0;j<nSFT;j++)
  { 
    /* fprintf(stderr, "sft %d  out of %d\n", j, nSFT); */
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
 
 /* Find time average of normalized SFTs */
 LALNormalizeSFTVect(&status, sft_vect, blocksRngMed);   
 timeavg = (REAL4 *)LALMalloc(numBins*sizeof(REAL4));
 for (j=0;j<nSFT;j++)
 { 
    for ( i=0; i < numBins; i++)
    {
        if (j == 0) {
          timeavg[i] = sft_vect->data[j].data->data[i].re*sft_vect->data[j].data->data[i].re + 
                       sft_vect->data[j].data->data[i].im*sft_vect->data[j].data->data[i].im;
        } else {
          timeavg[i] += sft_vect->data[j].data->data[i].re*sft_vect->data[j].data->data[i].re + 
                       sft_vect->data[j].data->data[i].im*sft_vect->data[j].data->data[i].im;
        }
    }

 }
 for ( i=0; i < numBins; i++)
 {
      f = sft_vect->data->f0 + ((REAL4)i)*sft_vect->data->deltaF;
      fprintf(fp3,"%16.8f %g\n",f,timeavg[i]/((REAL4)nSFT));
 } 

 LALDestroySFTVector (&status, &sft_vect );
 LALFree(timeavg);

 fclose(fp);
 fclose(fp2);
 fclose(fp3);

 LAL_CALL(LALDestroyUserVars(&status), &status);

 return(0);

}
/* END main */
