/*----------------------------------------------------------------------- 
 * 
 * File Name: Coherent.c
 *
 * Author: Seader, S. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */


/* USAGE:  ./a.out -HH1.dat -hH2.dat -LL.dat -GG.dat -VV.dat -TT.dat */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <getopt.h>

#include <lal/LALRCSID.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/DetectorSite.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/LIGOMetadataTables.h>
#include "CoherentInspiral.h"
/* CHECK: above: include <lal/CoherentInspiral.h>
 */

NRCSID (COHERENTINSPIRALFILTERC, "$Id$");

static void 
ParseOptions (int argc, char *argv[]);

extern char    *optarg;
extern int      optind;
static FILE *fp[4];
static FILE *fp2[4];

/* Parse Options - Hard code these and provide error checking */

char H1filename[256];
char H2filename[256];
char L1filename[256];
char GEOfilename[256];
char VIRGOfilename[256];
char TAMAfilename[256];

INT4 siteID[6] = {0,0,1,3,2,4}; /*H1 H2 L G V T*/
INT4 caseID[6] = {0,0,0,0,0,0};

static REAL4 cartesianInnerProduct(REAL4 x[3], REAL4 y[3])
{
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

void
LALCoherentInspiralFilterSegment (
    LALStatus                             *status,
    MultiInspiralTable                   **eventList,
    CoherentInspiralFilterInput           *input,
    CoherentInspiralFilterParams          *params
    )
{
  INT4                                k,l,p,q,w;
  INT4                                dsites[4] = {0,0,0,0};
  INT4                                nlines = 5;
  INT4                                nlines2 = 2; 
  INT4                                sampleRate; /* say, 1024Hz */
  INT4                                numPoints = 0;
  INT4                                slidePoints1 = 0;
  INT4                                slidePoints2 = 0;
  INT4                                slidePoints3 = 0;
  UINT4                               numdet = 0;  
  REAL4                               i,j,m,n;
  REAL4                               s[3],s2[3],s3[3];
  REAL4                               distance,distance2,distance3;
  REAL4                               deltaT;
  REAL4                               timeDelay[3] = {0,0,0};
  REAL4                               cohSNRLocal = 0;
  REAL4                               nHatVect[3] = {0,0,0};
  REAL4                              *cohSNR = NULL;
  LALDetectorPair                     detectors,detectors2,detectors3;
  COMPLEX8TimeSeries                 *zData;
  CoherentInspiralZVector            *multiZData;


  /*Now the structures that will store the coefficient data if necessary */
 
  REAL8                               detcoeff[4][nlines2][4]; 
 
  char  namearray[4][256] = {"0","0","0","0"};
  char  namearray2[4][256] = {"0","0","0","0"};
  char *argv[4];
  int   argc;

  ParseOptions (argc, argv);

  /* Now initialize the declared arrays to zero */
  /* if full coherent snr vector is required, set it to zero */
  numPoints = params->numPoints;
  
  if ( params->cohSNRVec )
    memset( params->cohSNRVec->data->data, 0, numPoints * sizeof( REAL4 ) );

  for (p=0;p<4;p++)
    {
    for (l=0;l<nlines2;l++)
      {
        for (k=0;k<4;k++)
	  {
	    detcoeff[p][l][k] = 0;
	  }
      }
    }

  for (l=0;l<7;l++)
    {
      if(caseID[l] == 1)
	numdet++;
    }


  fprintf(stdout, "You have specified %2d  detector(s).\n",numdet);
  fprintf(stdout, "The caseID is: %d %d %d %d %d %d (H1,H2,L1,GEO,VIRGO,TAMA) \n",caseID[0],caseID[1],caseID[2],caseID[3],caseID[4],caseID[5]);

  if (numdet > 4)
    {
    fprintf(stdout, "Too many detectors specified - exiting. \n");
    exit( 1 );
    }
  if (numdet == 0)
    {
      fprintf(stdout, "You must specify data filename(s) for 1 to 4 detectors - exiting. \n");
      exit( 1 );
    }

  /* Now read in all the filenames and store them in arrays */
  /* This will keep the files in order: H1 H2 L1 GEO VIRGO TAMA */

  if(caseID[0])
    {
    strcpy(namearray[0],H1filename);
    strcpy(namearray2[0],"H1coeff.txt");
    }

  if (caseID[0] && caseID[1])
    {
    strcpy(namearray[1],H2filename);
    strcpy(namearray2[1],"H2coeff.txt");
    }
  else if(caseID[1])
    {
    strcpy(namearray[0],H2filename);
    strcpy(namearray2[0],"H2coeff.txt");
    }

  if (caseID[0] && caseID[1] && caseID[2])
    {
    strcpy(namearray[2],L1filename);
    strcpy(namearray2[2],"Lcoeff.txt");
    }
  else if((caseID[0] || caseID[1]) && caseID[2])
    {
    strcpy(namearray[1],L1filename);
    strcpy(namearray2[1],"Lcoeff.txt");
    }
  else if(caseID[2])
    {
    strcpy(namearray[0],L1filename);
    strcpy(namearray2[0],"Lcoeff.txt");
    }

  if (caseID[0] && caseID[1] && caseID[2] && caseID[3])
    {
    strcpy(namearray[3],GEOfilename);
    strcpy(namearray2[3],"Gcoeff.txt");
    }
  else if(((caseID[0] && caseID[1]) || (caseID[0] && caseID[2]) || (caseID[1] && caseID[2])) && caseID[3])
    {
    strcpy(namearray[2],GEOfilename);
    strcpy(namearray2[2],"Gcoeff.txt");
    }
  else if((caseID[0] || caseID[1] || caseID[2]) && caseID[3])
    {
    strcpy(namearray[1],GEOfilename);
    strcpy(namearray2[1],"Gcoeff.txt");
    }
  else if(caseID[3])
    {
    strcpy(namearray[0],GEOfilename);
    strcpy(namearray2[0],"Gcoeff.txt");
    }

  if (caseID[0] && caseID[1] && caseID[2] && caseID[3])
    {
    }
  else if(((caseID[0] && caseID[1] && caseID[2]) || (caseID[0] && caseID[1] && caseID[3]) || (caseID[0] && caseID[2] && caseID[3]) || (caseID[1] && caseID[2] && caseID[3])) && caseID[4])
    {
    strcpy(namearray[3],VIRGOfilename);
    strcpy(namearray2[3],"Vcoeff.txt");
    }
  else if(((caseID[0] && caseID[1]) || (caseID[0] && caseID[2]) || (caseID[0] && caseID[3]) || (caseID[1] && caseID[2]) || (caseID[1] && caseID[3]) || (caseID[2] && caseID[3])) && caseID[4])
    {
    strcpy(namearray[2],VIRGOfilename);
    strcpy(namearray2[2],"Vcoeff.txt");
    }
  else if((caseID[0] || caseID[1] || caseID[2] || caseID[3]) && caseID[4])
    {
    strcpy(namearray[1],VIRGOfilename);
    strcpy(namearray2[1],"Vcoeff.txt");
    }
  else if(caseID[4])
    {
    strcpy(namearray[0],VIRGOfilename);
    strcpy(namearray2[0],"Vcoeff.txt");
    }

  if (caseID[0] && caseID[1] && caseID[2] && caseID[3] && caseID[4])
    {
    }
  else if (caseID[0] && caseID[1] && caseID[2] && caseID[3])
    {
    }
  else if(((caseID[0] && caseID[1] && caseID[2]) || (caseID[0] && caseID[1] && caseID[3]) || (caseID[0] && caseID[2] && caseID[3]) || (caseID[1] && caseID[2] && caseID[3])|| (caseID[0] && caseID[1] && caseID[4]) || (caseID[0] && caseID[2] && caseID[4]) || (caseID[0] && caseID[3] && caseID[4]) || (caseID[1] && caseID[2] && caseID[4]) || (caseID[1] && caseID[3] && caseID[4]) || (caseID[2] && caseID[3] && caseID[4])) && caseID[5])
    {
    strcpy(namearray[3],TAMAfilename);
    strcpy(namearray2[3],"Tcoeff.txt");
    }
  else if(((caseID[0] && caseID[1]) || (caseID[0] && caseID[2]) || (caseID[0] && caseID[3]) || (caseID[1] && caseID[2]) || (caseID[1] && caseID[3]) || (caseID[2] && caseID[3]) || (caseID[0] && caseID[4]) || (caseID[1] && caseID[4]) || (caseID[2] && caseID[4]) || (caseID[3] && caseID[4])) && caseID[5])
    {
    strcpy(namearray[2],TAMAfilename);
    strcpy(namearray2[2],"Tcoeff.txt");
    }
  else if((caseID[0] || caseID[1] || caseID[2] || caseID[3] || caseID[4]) && caseID[5])
    {
    strcpy(namearray[1],TAMAfilename);
    strcpy(namearray2[1],"Tcoeff.txt");
    }
  else if(caseID[5])
    {
    strcpy(namearray[0],TAMAfilename);
    strcpy(namearray2[0],"Tcoeff.txt");
    }

  /* Now that the order is sorted, the data will be read in */
  /* First, the files will be tested for length consistency */

  for(l=0;l<numdet;l++)
    {
      fp[l] = fopen(namearray[l], "r");

    }

  /* read in the z-data for multiple detectors */
  multiZData->length = input->multiZData->length;
  multiZData->zData = input->multiZData->zData;

  for ( l=0 ; l < multiZData->length ; l++) 
    {
      zData[l] = input->multiZData->zData[l]; 
    }


  /* read in the coefficients for each specified detector.  This should be moved to the case statements below where the coherent snr is computed.  It should be put in each of the last three cases: 3b,4a,4b.  This saves some unnecessary reading in of data since the coefficients are only used in the last 3 cases. Note that for networks with H1 and H2, there is some redundancy because they have the same coefficients.  currently it reads in a separate file for each of them. */

 for (l=0;l<numdet;l++)
      {
      fp2[l] = fopen(namearray2[l], "r");
      if(!fp2[l])
	{
	  fprintf(stdout,"The file %s containing the coefficients could not be found - exiting...\n",namearray2[l]);
	  exit(1);
	}
      for (k=0;k<nlines2;k++)
	  {
	  fscanf(fp2[l],"%f, %f, %f, %f",&i,&j,&m,&n);
	  detcoeff[l][k][0] = i;
	  detcoeff[l][k][1] = j;
	  detcoeff[l][k][2] = m;
	  detcoeff[l][k][3] = n; 
	  }
      fclose(fp2[l]);
      } 
   

  fprintf(stdout,"%8s %8s %8s %8s \n",namearray[0],namearray[1],namearray[2],namearray[3]);
  fprintf(stdout,"%8s %8s %8s %8s \n",namearray2[0],namearray2[1],namearray2[2],namearray2[3]);

  /* Now construct the appropriate coherent SNR */


  deltaT = 1/sampleRate;

  switch (numdet)
    {
    case 1:
      fprintf(stdout,"Network of 1 detector - no coherent search will  be performed. \n");
      /* CHECK: still need to check for threshold crossers 
	 and save then in events table */
      for (k=0;k<numPoints;k++)
	{
	  cohSNR[k] = sqrt(pow(zData[0].data->data[k].re,2) + pow(zData[0].data->data[k].im,2));
        }
      break;
    case 2:
      if(caseID[0] && caseID[1])
	{
	  fprintf(stdout,"Network of 2 detectors - includes H1 and H2 \n");
	  for (k=0;k<numPoints;k++)
	    {
             cohSNR[k] = (1 / sqrt(2)) * sqrt( pow(zData[0].data->data[k].re + zData[1].data->data[k].re,2) + pow( zData[0].data->data[k].im + zData[1].data->data[k].im,2));
	    } 
	}
      else 
	{
	  /*Here, the time delay looping must start */
	  fprintf(stdout,"Network of 2 detectors \n");

	  p=0;
	  for (l=0;l<7;l++)
	    {
	      if(caseID[l])
		{
		  dsites[p] = siteID[l];
		  p++;
		}
	    }

	  fprintf(stdout,"%d %d \n",dsites[0],dsites[1]);

	  /*Now compute the distance between the sites */

          detectors.detectorOne = lalCachedDetectors[dsites[0]];
          detectors.detectorTwo = lalCachedDetectors[dsites[1]];

          for (l=0;l<3;l++)
            {
             s[l] = (REAL4) ( detectors.detectorOne.location[l] - detectors.detectorTwo.location[l]);
            }
         /* Now calculate the distance in meters */
         distance = sqrt( cartesianInnerProduct(s,s) );

	 timeDelay[0] = distance/LAL_C_SI;
	 slidePoints1 = ceilf( fabsf(timeDelay[0])/deltaT );
	 fprintf(stdout,"%d \n",numPoints);

	 for(l=0;l<numPoints;l++)
	   {
	     for (n = l-slidePoints1;n < l+slidePoints1;n++)
	       {
		 if(n >= 0 && n < numPoints)
		   {
		     cohSNRLocal = sqrt(pow(zData[0].data->data[l].re,2) + pow(zData[1].data->data[(INT4) n].re,2) + pow(zData[0].data->data[l].im,2) + pow(zData[1].data->data[(INT4) n].im,2));
		   }
		 if(cohSNRLocal > cohSNR[l])
		   {
		   cohSNR[l] = cohSNRLocal;
		   }
	       }
	   }
	}
      break;
    case 3:
      if(caseID[0] && caseID[1])
	{
         fprintf(stdout,"Network of 3 detectors - includes H1 and H2 \n");
	 /*The SNR for the H1-H2 pair is first computed, then the looping is done with one time delay */

          p=0;
	  for (l=0;l<7;l++)
	    {
	      if(caseID[l])
		{
		  dsites[p] = siteID[l];
		  p++;
		}
	    }

	  fprintf(stdout,"%d %d %d \n",dsites[0],dsites[1],dsites[2]);

	  /*Now compute the distance between the sites */

          detectors.detectorOne = lalCachedDetectors[dsites[0]];
          detectors.detectorTwo = lalCachedDetectors[dsites[2]];

          for (l=0;l<3;l++)
            {
             s[l] = (REAL4) ( detectors.detectorOne.location[l] - detectors.detectorTwo.location[l]);
            }
         /* Now calculate the distance in meters */
         distance = sqrt( cartesianInnerProduct(s,s) );

	 timeDelay[0] = distance/LAL_C_SI;
	 slidePoints1 = ceilf( fabsf(timeDelay[0])/deltaT );
	 fprintf(stdout,"%d \n",slidePoints1);

           for(l=0;l<numPoints;l++)
	     {
	       for (n = l-slidePoints1;n < l+slidePoints1;n++)
	         {
		   if(n >= 0 && n < numPoints)
		     {
		       cohSNRLocal = sqrt( 0.5 * (pow(zData[0].data->data[l].re + zData[1].data->data[l].re,2) + pow(zData[0].data->data[l].im + zData[1].data->data[l].im,2)) + pow(zData[2].data->data[(INT4) n].re,2) + pow(zData[2].data->data[(INT4) n].im,2));


		     }
		   if(cohSNRLocal > cohSNR[l])
		     {
		     cohSNR[l] = cohSNRLocal;
		     }
	         }
	     }
	
	}
      else
	{
	  /* Now the last 3 cases will involve the looping over the coefficients */
          p=0;
	  for (l=0;l<7;l++)
	    {
	      if(caseID[l])
		{
		  dsites[p] = siteID[l];
		  p++;
		}
	    }
	  fprintf(stdout,"%d %d %d \n",dsites[0],dsites[1],dsites[2]);
	  for (l=0;l < nlines2; l++)
	    {
	      /*Now compute the distance between the sites */
	      
	      detectors.detectorOne = lalCachedDetectors[dsites[0]];
	      detectors.detectorTwo = lalCachedDetectors[dsites[1]];
	      detectors2.detectorOne = lalCachedDetectors[dsites[1]];
	      detectors2.detectorTwo = lalCachedDetectors[dsites[2]];
	      
	      for (k=0;k<3;k++)
		{
		  s[k] = (REAL4) ( detectors.detectorOne.location[k] - detectors.detectorTwo.location[k]);
		  s2[k] = (REAL4) ( detectors2.detectorOne.location[k] - detectors2.detectorTwo.location[k]);
		}
	      
	      
	      /* Now calculate the distance in meters */
	      nHatVect[0] = cos(detcoeff[0][l][0] * LAL_PI_180)*sin(detcoeff[0][l][1] * LAL_PI_180);
	      nHatVect[1] = sin(detcoeff[0][l][0] * LAL_PI_180)*sin(detcoeff[0][l][1] * LAL_PI_180);
	      nHatVect[2] = cos(detcoeff[0][l][1] * LAL_PI_180);
	      
	      distance = cartesianInnerProduct(s,nHatVect);
	      distance2 = cartesianInnerProduct(s2,nHatVect);
	      timeDelay[0] = distance/LAL_C_SI;
	      timeDelay[1] = distance2/LAL_C_SI;
	      slidePoints1 = ceilf( fabsf(timeDelay[0])/deltaT );
	      slidePoints2 = ceilf( fabsf(timeDelay[1])/deltaT );
	      
	      fprintf(stdout,"%d %d\n",slidePoints1,slidePoints2);
	      
	      for(p=0;p<numPoints;p++)
		{
		  for (n = p-slidePoints1;n < p+slidePoints1;n++)
		    {
		      if(n >= 0 && n < numPoints)
			{
			  
			  for (q = n-slidePoints2; q < n+slidePoints2;q++)
			    {
			      if (q >= 0 && q < numPoints)
				{
				  cohSNRLocal = sqrt( (pow(detcoeff[0][l][2],2) + pow(detcoeff[0][l][3],2))*(pow(zData[0].data->data[p].re,2) + pow(zData[0].data->data[p].im,2)) + (pow(detcoeff[1][l][2],2) + pow(detcoeff[1][l][3],2))*(pow(zData[1].data->data[(INT4) n].re,2) + pow(zData[1].data->data[(INT4) n].im,2)) + (pow(detcoeff[2][l][2],2) + pow(detcoeff[2][l][3],2))*(pow(zData[2].data->data[(INT4) q].re,2) + pow(zData[2].data->data[(INT4) q].im,2)) + 2*(detcoeff[0][l][2]*detcoeff[1][l][2] + detcoeff[0][l][3]*detcoeff[1][l][3])*(zData[0].data->data[p].re * zData[1].data->data[(INT4) n].re + zData[0].data->data[p].im * zData[1].data->data[(INT4) n].im) + 2*(detcoeff[0][l][2]*detcoeff[2][l][2] + detcoeff[0][l][3]*detcoeff[2][l][3])*(zData[0].data->data[p].re*zData[2].data->data[(INT4) n].re + zData[0].data->data[p].im * zData[2].data->data[(INT4) n].im) + 2*(detcoeff[1][l][2]*detcoeff[2][l][2] + detcoeff[1][l][3]*detcoeff[2][l][3])*(zData[1].data->data[(INT4) n].re * zData[2].data->data[(INT4) q].re + zData[1].data->data[(INT4) n].im *zData[2].data->data[(INT4) q].im));

				  fprintf(stdout,"%f \n", cohSNRLocal);
				}
			      if(cohSNRLocal > cohSNR[p])
				{
				  cohSNR[p] = cohSNRLocal;
				}
			      
			    }
			  
			}
		      
		    }
		}
	      
	    }/* outer loop end */
	  
	} /* else statement end */
      break;
    case 4:
      if(caseID[0] && caseID[1])
	{
	  fprintf(stdout,"Network of 4 detectors - includes H1 and H2 \n");

         p=0;
	  for (l=0;l<7;l++)
	    {
	      if(caseID[l])
		{
		  dsites[p] = siteID[l];
		  p++;
		}
	    }
	  fprintf(stdout,"%d %d %d %d\n",dsites[0],dsites[1],dsites[2],dsites[3]);

	  /*start search looping */

	  for (l=0;l < nlines2; l++) 
	    {
	      /*Now compute the distance between the sites */
	      
	      detectors.detectorOne = lalCachedDetectors[dsites[0]];
	      detectors.detectorTwo = lalCachedDetectors[dsites[1]];
	      detectors2.detectorOne = lalCachedDetectors[dsites[1]];
	      detectors2.detectorTwo = lalCachedDetectors[dsites[2]];
	      detectors3.detectorOne = lalCachedDetectors[dsites[2]];
	      detectors3.detectorTwo = lalCachedDetectors[dsites[3]];
	      
	      for (k=0;k<3;k++)
		{
		  s[k] = (REAL4) ( detectors.detectorOne.location[k] - detectors.detectorTwo.location[k]);
		  s2[k] = (REAL4) ( detectors2.detectorOne.location[k] - detectors2.detectorTwo.location[k]);
		  s3[k] = (REAL4) ( detectors3.detectorOne.location[k] - detectors3.detectorTwo.location[k]);
		}
	      
	      /* Now calculate the distance in meters */
	      nHatVect[0] = cos(detcoeff[0][l][0] * LAL_PI_180)*sin(detcoeff[0][l][1] * LAL_PI_180);
	      nHatVect[1] = sin(detcoeff[0][l][0] * LAL_PI_180)*sin(detcoeff[0][l][1] * LAL_PI_180);
	      nHatVect[2] = cos(detcoeff[0][l][1] * LAL_PI_180);
	      
	      distance = cartesianInnerProduct(s,nHatVect);
	      distance2 = cartesianInnerProduct(s2,nHatVect);
	      distance3 = cartesianInnerProduct(s3,nHatVect);
	      timeDelay[0] = distance/LAL_C_SI; /*check:this should be zero (H1-H2 delay)*/
	      timeDelay[1] = distance2/LAL_C_SI;
	      timeDelay[2] = distance3/LAL_C_SI;
	      slidePoints1 = ceilf( fabsf(timeDelay[0])/deltaT ); /*This should be zero since the delay=0*/
	      slidePoints2 = ceilf( fabsf(timeDelay[1])/deltaT );
	      slidePoints3 = ceilf( fabsf(timeDelay[2])/deltaT );
	      
	      fprintf(stdout,"%d %d %d\n",slidePoints1,slidePoints2,slidePoints3);
	      
	      for(p=0;p<numPoints;p++)
		{
		  for (n = p-slidePoints2;n < p+slidePoints2;n++)
		    {
		      if(n >= 0 && n < numPoints)
			{
			  
			  for (q = n-slidePoints3; q < n+slidePoints3;q++)
			    {
			      if (q >= 0 && q < numPoints)
				{
				  cohSNRLocal = sqrt( 0.5*(pow(detcoeff[0][l][2],2) + pow(detcoeff[0][l][3],2))*(pow(zData[0].data->data[p].re,2) + pow(zData[0].data->data[p].im,2) + pow(zData[1].data->data[p].re,2) + pow(zData[1].data->data[p].im,2) + 2*(zData[0].data->data[p].re*zData[1].data->data[p].re + zData[0].data->data[p].im*zData[1].data->data[p].im)) + (pow(detcoeff[2][l][2],2) + pow(detcoeff[2][l][3],2))*(pow(zData[2].data->data[(INT4) n].re,2) + pow(zData[2].data->data[(INT4) n].im,2)) + (pow(detcoeff[3][l][2],2) + pow(detcoeff[3][l][3],2))*(pow(zData[3].data->data[(INT4) q].re,2) + pow(zData[3].data->data[(INT4) q].im,2)) + sqrt(2)*(detcoeff[0][l][2]*detcoeff[2][l][2] + detcoeff[0][l][3]*detcoeff[2][l][3])*(zData[0].data->data[p].re*zData[2].data->data[(INT4) n].re + zData[0].data->data[p].im*zData[2].data->data[(INT4) n].im + zData[1].data->data[p].re*zData[2].data->data[(INT4) n].re + zData[1].data->data[p].im*zData[2].data->data[(INT4) n].im) + sqrt(2)*(detcoeff[0][l][2]*detcoeff[3][l][2] + detcoeff[0][l][3]*detcoeff[3][l][3])*(zData[0].data->data[p].re*zData[3].data->data[(INT4) q].re + zData[0].data->data[p].im*zData[3].data->data[(INT4) q].im + zData[1].data->data[p].re*zData[3].data->data[(INT4) q].re + zData[1].data->data[p].im*zData[3].data->data[(INT4) q].im) + 2*(detcoeff[2][l][2]*detcoeff[3][l][2] + detcoeff[2][l][3]*detcoeff[3][l][3])*(zData[2].data->data[(INT4) n].re*zData[3].data->data[(INT4) q].re + zData[2].data->data[(INT4) n].im*zData[3].data->data[(INT4) q].im));
				  
				  fprintf(stdout,"%f \n", cohSNRLocal);
				}
			      if(cohSNRLocal > cohSNR[p])
				{
				  cohSNR[p] = cohSNRLocal;
				}
			      
			    }
			  
			}
		      
		    }
		}
	      
	      
	    } /* end for statement prior to computing distances*/
	  
	  
	  
	} /* end outer if statement in case4*/
      else
	{
	  fprintf(stdout,"Network of four detctors \n");
	  /* there will be one extra loop over the last case since there are 3 nonzero time delays*/
	  
	  p=0;
	  for (l=0;l<7;l++)
	    {
	      if(caseID[l])
		{
		  dsites[p] = siteID[l];
		  p++;
		}
	    }
	  fprintf(stdout,"%d %d %d %d\n",dsites[0],dsites[1],dsites[2],dsites[3]);
	  
	  /*start search looping */
	  
	  for (l=0;l < nlines2; l++)
	    {
	      /*Now compute the distance between the sites */
	      
	      detectors.detectorOne = lalCachedDetectors[dsites[0]];
	      detectors.detectorTwo = lalCachedDetectors[dsites[1]];
	      detectors2.detectorOne = lalCachedDetectors[dsites[1]];
	      detectors2.detectorTwo = lalCachedDetectors[dsites[2]];
	      detectors3.detectorOne = lalCachedDetectors[dsites[2]];
	      detectors3.detectorTwo = lalCachedDetectors[dsites[3]];
	      
	      for (k=0;k<3;k++)
		{
                  s[k] = (REAL4) ( detectors.detectorOne.location[k] - detectors.detectorTwo.location[k]);
	          s2[k] = (REAL4) ( detectors2.detectorOne.location[k] - detectors2.detectorTwo.location[k]);
		  s3[k] = (REAL4) ( detectors3.detectorOne.location[k] - detectors3.detectorTwo.location[k]);
		}
	      
	      /* Now calculate the distance in meters */
	      nHatVect[0] = cos(detcoeff[0][l][0] * LAL_PI_180)*sin(detcoeff[0][l][1] * LAL_PI_180);
	      nHatVect[1] = sin(detcoeff[0][l][0] * LAL_PI_180)*sin(detcoeff[0][l][1] * LAL_PI_180);
	      nHatVect[2] = cos(detcoeff[0][l][1] * LAL_PI_180);
	      
	      distance = cartesianInnerProduct(s,nHatVect);
	      distance2 = cartesianInnerProduct(s2,nHatVect);
	      distance3 = cartesianInnerProduct(s3,nHatVect);
	      timeDelay[0] = distance/LAL_C_SI; 
	      timeDelay[1] = distance2/LAL_C_SI;
	      timeDelay[2] = distance3/LAL_C_SI;
	      slidePoints1 = ceilf( fabsf(timeDelay[0])/deltaT ); 
	      slidePoints2 = ceilf( fabsf(timeDelay[1])/deltaT );
	      slidePoints3 = ceilf( fabsf(timeDelay[2])/deltaT );
	      
	      fprintf(stdout,"%d %d %d\n",slidePoints1,slidePoints2,slidePoints3);
	      
	      for(p=0;p<numPoints;p++)
		{
		  for (n = p-slidePoints1;n < p+slidePoints1;n++)
		    {
		      if(n >= 0 && n < numPoints)
			{
			  for (q = n-slidePoints2; q < n+slidePoints2;q++)
			    {
			      if (q >= 0 && q < numPoints)
				{
				  for(w = q-slidePoints3; w < q+slidePoints3; w++)
				    {
				      if(w >=0 && w < numPoints)
					{
					  cohSNRLocal = sqrt( (pow(detcoeff[0][l][2],2) + pow(detcoeff[0][l][3],2))*(pow(zData[0].data->data[p].re,2) + pow(zData[0].data->data[p].im,2)) + (pow(detcoeff[1][l][2],2) + pow(detcoeff[1][l][3],2))*(pow(zData[1].data->data[(INT4) n].re,2) + pow(zData[1].data->data[(INT4) n].im,2)) + (pow(detcoeff[2][l][2],2) + pow(detcoeff[2][l][3],2))*(pow(zData[2].data->data[(INT4) q].re,2) + pow(zData[2].data->data[(INT4) q].im,2)) + (pow(detcoeff[3][l][2],2) + pow(detcoeff[3][l][3],2))*(pow(zData[3].data->data[(INT4) w].re,2) + pow(zData[3].data->data[(INT4) w].im,2)) + 2*(detcoeff[0][l][2]*detcoeff[1][l][2] + detcoeff[0][l][3]*detcoeff[1][l][3])*(zData[0].data->data[p].re * zData[1].data->data[(INT4) n].re + zData[0].data->data[p].im * zData[1].data->data[(INT4) n].im) + 2*(detcoeff[0][l][2]*detcoeff[2][l][2] + detcoeff[0][l][3]*detcoeff[2][l][3])*(zData[0].data->data[p].re*zData[2].data->data[(INT4) n].re + zData[0].data->data[p].im * zData[2].data->data[(INT4) n].im) + 2*(detcoeff[0][l][2]*detcoeff[3][l][2] + detcoeff[0][l][3]*detcoeff[3][l][3])*(zData[0].data->data[p].re*zData[3].data->data[(INT4) w].re + zData[0].data->data[p].im*zData[3].data->data[(INT4) w].im) + 2*(detcoeff[1][l][2]*detcoeff[2][l][2] + detcoeff[1][l][3]*detcoeff[2][l][3])*(zData[1].data->data[(INT4) n].re * zData[2].data->data[(INT4) q].re + zData[1].data->data[(INT4) n].im *zData[2].data->data[(INT4) q].im) + 2*(detcoeff[1][l][2]*detcoeff[3][l][2] + detcoeff[1][l][3]*detcoeff[3][l][3])*(zData[1].data->data[(INT4) n].re*zData[3].data->data[(INT4) w].re + zData[1].data->data[(INT4) n].im*zData[3].data->data[(INT4) w].im) + 2*(detcoeff[2][l][2]*detcoeff[3][l][2] + detcoeff[2][l][3]*detcoeff[3][l][3])*(zData[2].data->data[(INT4) q].re*zData[3].data->data[(INT4) w].re + zData[2].data->data[(INT4) q].im * zData[3].data->data[(INT4) w].im));
					  
					  
					  fprintf(stdout,"%f \n", cohSNRLocal);
					}
				      if(cohSNRLocal > cohSNR[p])
					{
					  cohSNR[p] = cohSNRLocal;
					}
				    }
				}
			    }
			}
		    }
		}
	      
	    } /*end the outermost for statement prior to computing distances*/
	  
	} /*end else statement */
      break; 
    } /* end case statement */
  
  
  
    
  
  fprintf(stdout,"%6.2f %6.2f %6.2f %6.2f %6.2f \n",cohSNR[0],cohSNR[1],cohSNR[2],cohSNR[3],cohSNR[4]);
  
  exit(0);
  
} /* end main */


  /* may want to add a verbose option flag to print out extra info here and there */
static void
ParseOptions (
	      int         argc, 
	      char       *argv[]
	      )
{
  while (1)
    {
      int c = -1;
      
      c = getopt (argc, argv, "H:""h:""L:""G:""V:""T:");
      if (c == -1)
	{
	  break;
	}
      
      switch (c)
	{
	case 'H': 
	  strcpy(H1filename,optarg);
	  caseID[0] = 1;
	  break;
	case 'h':
	  strcpy(H2filename,optarg);
	  caseID[1] = 1;
	  break;
	case 'L':
	  strcpy(L1filename,optarg);
	  caseID[2] = 1;
	  break;
	case 'G':
	  strcpy(GEOfilename,optarg);
	  caseID[3] = 1;  
	  break;
	case 'V':
	  strcpy(VIRGOfilename,optarg);
	  caseID[4] = 1;	   
	  break;
	case 'T':
	  strcpy(TAMAfilename,optarg);
	  caseID[5] = 1;      
	  break;
        default:
        fprintf( stderr, "unknown error while parsing options\n" );
        exit( 1 );

	}
    }
  return;
}
