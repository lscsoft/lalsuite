/*
*  Copyright (C) 2007  Rejean Dupuis
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

/**
 * \file
 * \ingroup pulsarApps
 * \author R.J. Dupuis
 * \brief This program injects a signal in the Bk's
 */

/****************************************************************************
R.J. Dupuis 18 April, 2004
This program injects a signal in the Bk's... 

reads in (from inputs) S2_h95_inj.txt,and pulsar parameter file 
*************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/SkyCoordinates.h> 
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 1; } else ((void)0)


int main(int argc, char **argv)
{
  static LALStatus      status;  
  FILE *fpout, *fpin, *fp, *psrfp;
  char infile[256], indir[256], outdir[256], inname[256], outname[256], dect[16];
  char pulsarname[256];
  REAL8 t, cosiota, cos2iota, phase, psi, h0, RA, DEC;
  int iRA=0, iDEC=0;
  COMPLEX16 h, B;
  REAL8 val;
  char txt[16], psrinput[256];
  
  LIGOTimeGPS tgps;
  LALDetector 		detector;
  LALSource 		pulsar;
  LALDetAndSource 	detAndSource;
  LALDetAMResponse 	computedResponse;  
  LALGPSandAcc		pGPSandAcc;

  
  
  
  if(argv[1]==NULL||argv[2]==NULL||argv[3]==NULL||argv[4]==NULL) 
  {
    printf("Error! infile indir outdir detector\n");
    return(0);
  }  
  
  sprintf(infile, "%s", argv[1]);
  sprintf(indir, "%s", argv[2]);
  sprintf(outdir, "%s", argv[3]);
  sprintf(dect, "%s", argv[4]);
        
  if (!strcmp(dect,"GEO")) detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF]; 
  else if (!strcmp(dect,"L1"))  detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if (!strcmp(dect,"H1") || !strcmp(dect,"H2"))  detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  else 
  {
    printf("non valid dectector!\n");
    return(1);
  }
  
  
  fp = fopen(infile, "r");
  
  while (5==fscanf(fp, "%s %lf %lf %lf %lf", &pulsarname[0], &h0, &cosiota, &phase, &psi))
  {
     printf("%s\t%e\t%f\t%f\t%f\n", pulsarname, h0, cosiota, phase, psi);
    /************** BEGIN READING PULSAR PARAMETERS ******************/     
    sprintf(psrinput,"inputs/%s", pulsarname);
    psrfp=fopen(psrinput,"r"); 
    while (2==fscanf(psrfp,"%s %lf", &txt[0], &val))
  {
    if( !strcmp(txt,"ra") || !strcmp(txt,"RA")) {
      RA = val;
      iRA = 1;
    }  
    
    else if( !strcmp(txt,"dec") || !strcmp(txt,"DEC")) {
      DEC = val;
      iDEC = 1;
    }
  }  
  if ((iRA+iDEC) != 2) {
    fprintf(stderr, "missing RA DEC");
    return(3);
  }
  
  fclose(psrfp);
  
    /************** END READING PULSAR PARAMETERS ******************/   
  
  pulsar.equatorialCoords.longitude = RA;	      
  pulsar.equatorialCoords.latitude = DEC;	      
  pulsar.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;
  pulsar.orientation = psi;			      
  strcpy(pulsar.name, "123400");		      

  detAndSource.pSource = &pulsar;
  detAndSource.pDetector = &detector;   
  
  sprintf(inname, "%s/outfine.%s_%s.S2", indir, pulsarname, dect);
  fpin = fopen(inname, "r");
  
  sprintf(outname, "%s/outfine.%s_%s.S2", outdir, pulsarname, dect);
  fpout = fopen(outname, "w");
  
  cos2iota = 1.0 + cosiota*cosiota;
  
  while (3==fscanf(fpin, "%lf %lf %lf", &t, &B.re, &B.im))
  {
    tgps.gpsSeconds = (INT4)floor(t);
    tgps.gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));
    pGPSandAcc.gps =   tgps; 
    pGPSandAcc.accuracy = 1.0; 
    LALComputeDetAMResponse(&status, &computedResponse,&detAndSource, &pGPSandAcc);
    h.re = 0.25*computedResponse.plus*h0*cos2iota*cos(phase*2.0)+ 0.5*h0*computedResponse.cross*cosiota*sin(phase*2.0);
    h.im = 0.25*computedResponse.plus*h0*cos2iota*sin(phase*2.0) -0.5*h0*computedResponse.cross*cosiota*cos(phase*2.0);
    fprintf(fpout,"%f\t%e\t%e\n",  (double)t,  h.re+B.re ,  h.im+B.im);
  }
   fclose(fpin);
   fclose(fpout);
}

  fclose(fp);
  return(0);
}
