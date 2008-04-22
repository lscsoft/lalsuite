/*
 *  Copyright (C) 2007 Badri Krishnan  
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



#include "./radiometer.h"
/* lalapps includes */
#include <lalapps.h>
#include <lal/DopplerScan.h>
#include <gsl/gsl_permutation.h>

RCSID( "$Id$");



/* globals, constants and defaults */


extern int lalDebugLevel;


#define EARTHEPHEMERIS "/Users/artax/work/opt/lscsoft/lal/share/lal/earth05-09.dat"
#define SUNEPHEMERIS "/Users/artax/work/opt/lscsoft/lal/share/lal/sun05-09.dat"

#define F0 100
#define FBAND 1

#define BLOCKSRNGMED 51
#define MAXFILENAMELENGTH 512 /* maximum # of characters  of a filename */

#define DIROUT "./out"   /* output directory */
#define BASENAMEOUT "radio"    /* prefix file output */

#define SKYFILE "./skypatchfile"      
#define SKYREGION "allsky" 

#define TRUE (1==1)
#define FALSE (1==0)

void CreateSFTPairsIndicesFrom2SFTvectors(LALStatus                *status,
					 INT4VectorSequence        *out,
					 const SFTVector           *in1,
					 const SFTVector           *in2,
					  SFTPairParams             *par);




/** create pairs of sfts */
void CreateSFTIndexPairs(LALStatus                *status,
			 INT4VectorSequence       *out,
			 MultiSFTVector           *inputSFTs,
			 SFTPairParams            *par)
{  


  UINT4 numifo;
  UINT4 i;
  INITSTATUS (status, "CreateSFTIndexPairs", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (out, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (inputSFTs, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);

  numifo = inputSFTs->length;
  /* for now require exactly 2 ifos -- to be relaxed very soon in the future */
  if ( numifo != 2) {
    ABORT ( status, RADIOMETER_EBAD, RADIOMETER_MSGEBAD );
  }

  TRY( CreateSFTPairsIndicesFrom2SFTvectors(status->statusPtr, out, inputSFTs->data[0], 
					    inputSFTs->data[1], par), status);


  for (i=0; i<numifo; i++) {
	printf("ifo %i is %s\n",i+1,inputSFTs->data[i]->data->name);
  }


  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);
} /* CreateSFTIndexPairs */


void CreateSFTPairsIndicesFrom2SFTvectors(LALStatus                *status,
					 INT4VectorSequence        *out,
					 const SFTVector           *in1,
					 const SFTVector           *in2,
					 SFTPairParams             *par)
{
  
  UINT4 i, j, numsft1, numsft2, numPairs, numsftMin;
  INT4 thisPair;
  COMPLEX8FrequencySeries  *thisSFT1, *thisSFT2;	       

  INITSTATUS (status, "CreateSFTPairsIndicesFrom2SFTvectors", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (out, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (in1, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (in2, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (par, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);

  /* number of SFTs from the two ifos */
  numsft1 = in1->length;
  numsft2 = in2->length;
  numsftMin = (numsft1 < numsft2)? numsft1 : numsft2;

  out->length = numsft1*numsft2;

  numPairs = out->length; /* initialization */

  /* increase length of sftpair vector */
  out->length += numsftMin;
  out->data = LALRealloc(out->data, out->length * out->vectorLength * sizeof(out->data[0]));

  thisPair = 0;

  
  /* go over all sfts in first vector */
  for (i=0; i<numsft1; i++) {

    thisSFT1 = in1->data + i;
    
    /* go over all sfts in second vector and check if it should be paired with thisSFT1 */
    /* this can be made more efficient if necessary */
    for (j=0; j<numsft2; j++) {

      LIGOTimeGPS t1, t2;
      REAL8 timeDiff;

      thisSFT2 = in2->data + j;

      /* calculate time difference */      
      t1 = thisSFT1->epoch;
      t2 = thisSFT2->epoch;
      timeDiff = XLALDeltaFloatGPS( &t1, &t2);

      /* decide whether to add this pair or not */
      if ( fabs(timeDiff) < par->lag ) {

	numPairs++;

	if ( numPairs < out->length) {
	  /* there is enough memory 
	  thisPair++;
	  *thisPair = i;
	  thisPair++;
	  *thisPair = j;*/
	out->data[thisPair++] = i;
	out->data[thisPair++] = j;
/*printf("pair indices %i %i\n", out->data[thisPair-2], out->data[thisPair-1]);*/

	} 
	else {
	  /* there not enough memory -- allocate memory and add the pair */
	  out->length += numsftMin;
	  out->data = LALRealloc(out->data, out->length * out->vectorLength * sizeof(out->data[0]));

	 /* thisPair++;
	  *thisPair = i;
	  thisPair++;
	  *thisPair = j;*/
	out->data[thisPair++] = i;
	out->data[thisPair++] = j;

	}

      } /* if ( numPairs < out->length)  */
    } /* end loop over second sft set */
  } /* end loop over first sft set */ 



  /* realloc memory -- reduce to exact number of pairs*/
  out->length = numPairs;
  out->data = LALRealloc(out->data, out->length * out->vectorLength * sizeof(out->data[0]));


  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

} /* CreateSFTPairsIndicesFrom2SFTvectors */




/** create pairs of sfts */
void CreateSFTPairs(LALStatus                *status,
		    SFTPairVec               *out,
		    MultiSFTVector           *inputSFTs,
		    MultiDetectorStateSeries *mdetStates,
		    SFTPairParams            *par)
{  

  UINT4 numifo;

  INITSTATUS (status, "CreateSFTPairs", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (out, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (inputSFTs, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (mdetStates, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (inputSFTs->length == mdetStates->length, status, RADIOMETER_EBAD, RADIOMETER_MSGEBAD);

  numifo = inputSFTs->length;
  /* for now require exactly 2 ifos -- to be relaxed very soon in the future */
  if ( numifo != 2) {
    ABORT ( status, RADIOMETER_EBAD, RADIOMETER_MSGEBAD );
  }

  TRY( CreateSFTPairsFrom2SFTvectors(status->statusPtr, out, inputSFTs->data[0], 
				     inputSFTs->data[1], mdetStates->data[0], 
				     mdetStates->data[1], par), status);

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);
}





/** create pairs of sfts from a pair of sft vectors*/
void CreateSFTPairsFrom2SFTvectors(LALStatus                 *status,
				   SFTPairVec                *out,
				   const SFTVector           *in1,
				   const SFTVector           *in2,
				   const DetectorStateSeries *det1,
				   const DetectorStateSeries *det2,
				   SFTPairParams             *par)
{
  
  UINT4 i, j, numsft1, numsft2, numPairs, numsftMin, counter=0;

  COMPLEX8FrequencySeries  *thisSFT1, *thisSFT2;	       
  DetectorState *thisDetState1, *thisDetState2;
  SingleSFTpair *thisPair = NULL;

  INITSTATUS (status, "CreateSFTPairsFrom2SFTvectors", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (out, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (in1, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (in2, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (det1, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (det2, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (par, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);

  /* number of SFTs from the two ifos */
  numsft1 = in1->length;
  numsft2 = in2->length;
  numsftMin = (numsft1 < numsft2)? numsft1 : numsft2;

  numPairs = out->length; /* initialization */

  /* increase length of sftpair vector */
  out->length += numsftMin;
  out->data = LALRealloc(out->data, out->length * sizeof(out->data[0]));

  /* need to allocate memory for singleSFTpair */
  thisPair = (SingleSFTpair *) LALCalloc(1, sizeof(out->data[0]));

printf("first detector start time %i\n", det1->data->tGPS.gpsSeconds);
printf("second detector start time %i\n",det2->data->tGPS.gpsSeconds);

  /* go over all sfts in first vector */
  for (i=0; i<numsft1; i++) {

    thisSFT1 = in1->data + i;
    thisDetState1 = det1->data + i;
    
    /* go over all sfts in second vector and check if it should be paired with thisSFT1 */
    /* this can be made more efficient if necessary */
    for (j=0; j<numsft2; j++) {

      LIGOTimeGPS t1, t2;
      REAL8 timeDiff;

      thisSFT2 = in2->data + j;
      thisDetState2 = det2->data + i;

      /* calculate time difference */      
      t1 = thisSFT1->epoch;
      t2 = thisSFT2->epoch;
      timeDiff = XLALDeltaFloatGPS( &t1, &t2);

      /* decide whether to add this pair or not */
      if ( fabs(timeDiff) < par->lag ) {
	numPairs++;


	if ( numPairs < out->length) {
	  /* there is enough memory */
/*	  thisPair++;*/




	  TRY( FillSFTPair( status->statusPtr, thisPair, thisSFT1, thisSFT2, thisDetState1, thisDetState2), status);
	out->data[counter] = *thisPair;
	counter++;
/*printf("counter %i\n",counter);*/
	} 
	else {
	  /* there not enough memory -- allocate memory and add the pair */

	  out->length += numsftMin;
	  out->data = LALRealloc(out->data, out->length * sizeof(out->data[0]));

/*	  thisPair++;*/

	  TRY( FillSFTPair( status->statusPtr, thisPair, thisSFT1, thisSFT2, thisDetState1, thisDetState2), status);
	out->data[counter] = *thisPair;
	counter++;
	}

      }
    } /* end loop over second sft set */
  } /* end loop over first sft set */ 


  /* realloc memory */
  out->length = numPairs;
  out->data = LALRealloc(out->data, out->length * sizeof(out->data[0]));

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}





/** little helper function for filling up sft pair */
void FillSFTPair(LALStatus                 *status,
		 SingleSFTpair             *out,
		 COMPLEX8FrequencySeries   *sft1, 
		 COMPLEX8FrequencySeries   *sft2, 
		 DetectorState             *det1,
		 DetectorState             *det2)
{  
  INITSTATUS (status, "FillSFTPairs", rcsid);
  ATTATCHSTATUSPTR (status);
    

  out->sft1 = sft1;
  out->sft2 = sft2;

  out->vel1[0] = det1->vDetector[0];
  out->vel1[1] = det1->vDetector[1];
  out->vel1[2] = det1->vDetector[2];
  
  out->vel2[0] = det2->vDetector[0];
  out->vel2[1] = det2->vDetector[1];
  out->vel2[2] = det2->vDetector[2];
  
  out->pos1[0] = det1->rDetector[0];
  out->pos1[1] = det1->rDetector[1];
  out->pos1[2] = det1->rDetector[2];
  
  out->pos2[0] = det2->rDetector[0];
  out->pos2[1] = det2->rDetector[1];
  out->pos2[2] = det2->rDetector[2];

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}


/** Correlate a single pair of SFT at a parameter space point*/
void CorrelateSingleSFTPair(LALStatus                *status,
			    COMPLEX16                *out,
			    COMPLEX8FrequencySeries  *sft1,
			    COMPLEX8FrequencySeries  *sft2,
			    REAL8                    *freq1,
			    REAL8                    *freq2)
{  
  UINT8 bin1, bin2;
  REAL8 deltaF;
  REAL8 re1, re2, im1, im2;

  INITSTATUS (status, "CorrelateSingleSFTPair", rcsid);
  ATTATCHSTATUSPTR (status);

  /* assume both sfts have the same freq. resolution */
  deltaF = sft1->deltaF;

  bin1 = (UINT8)( (*freq1 - sft1->f0) / deltaF);
  bin2 = (UINT8)( (*freq2 - sft2->f0)/ deltaF);

/*printf("freq1 %f freq2 %f\n", *freq1, *freq2);
printf("bin1 %i bin2 %i\n", (INT4)bin1, (INT4)bin2);*/

  re1 = sft1->data->data[bin1].re;
  im1 = sft1->data->data[bin1].im;
  re2 = sft2->data->data[bin2].re;
  im2 = sft2->data->data[bin2].im;

  out->re = deltaF * deltaF * (re1*re2 + im1*im2);
  out->im = deltaF * deltaF * (re1*im2 - re2*im1);

/*printf("real1, real2 %f %f, im1, im2, %f, %f\n",re1, re2, im1,im2);*/
/*printf("real %f imaginary %f\n",out->re, out->im); */

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}




/** Correlate a single pair of SFT at a parameter space point*/
void GetSignalFrequencyInSFT(LALStatus                *status,
			     REAL8                    *out,
			     COMPLEX8FrequencySeries  *sft1,
			     PulsarDopplerParams      *dopp,
			     REAL8Vector              *vel)
{  
  UINT4 k;
  REAL8 alpha, delta;
  REAL8 vDotn, fhat, factor, timeDiff;
  
  INITSTATUS (status, "GetSignalFrequencyInSFT", rcsid);
  ATTATCHSTATUSPTR (status);

  alpha = dopp->Alpha;
  delta = dopp->Delta;

  vDotn = sin(delta) * cos(alpha) * vel->data[0]
    + sin(delta) * sin(alpha) * vel->data[1]
    + cos(delta) * vel->data[2];

  /* now calculate the intrinsic signal frequency in the SFT */
  /* fhat = f_0 + f_1(t-t0) + f_2(t-t0)^2/2 + ... */

  /* this is the sft reference time  - the pulsar reference time */
  timeDiff = XLALDeltaFloatGPS( &(sft1->epoch), &(dopp->refTime));

  fhat = dopp->fkdot[0]; /* initialization */
  factor = 1.0;
  for (k = 1;  k < PULSAR_MAX_SPINS; k++) {
    factor *= timeDiff / k;  
    fhat += dopp->fkdot[k] * factor;
  }

  *out = fhat * (1 + vDotn);  

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);
}



/** Get signal phase at a given epoch */
void GetSignalPhaseInSFT(LALStatus            *status,
			 REAL8                *out,
			 LIGOTimeGPS          *epoch,
			 PulsarDopplerParams  *dopp,
			 REAL8Vector          *pos)
{  
  UINT4 k;
  REAL8 alpha, delta;
  REAL8 rDotn, fhat, phihat, factor, timeDiff;
  
  INITSTATUS (status, "GetSignalPhaseInSFT", rcsid);
  ATTATCHSTATUSPTR (status);

  alpha = dopp->Alpha;
  delta = dopp->Delta;

  rDotn = sin(delta) * cos(alpha) * pos->data[0]
    + sin(delta) * sin(alpha) * pos->data[1]
    + cos(delta) * pos->data[2];

  /* now calculate the intrinsic signal frequency in the SFT */
  /* fhat = f_0 + f_1(t-t0) + f_2(t-t0)^2/2 + ... */
  timeDiff = XLALDeltaFloatGPS( epoch, &(dopp->refTime));

  fhat = dopp->fkdot[0]; /* initialization */
  phihat = 0.0;

  factor = 1.0;
  for (k = 1;  k < PULSAR_MAX_SPINS; k++) {
    factor *= timeDiff / k;  

    fhat += dopp->fkdot[k] * factor;
    phihat += dopp->fkdot[k-1] * factor;
  }

  *out = LAL_TWOPI * ( phihat + fhat * rDotn);


  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);
}

