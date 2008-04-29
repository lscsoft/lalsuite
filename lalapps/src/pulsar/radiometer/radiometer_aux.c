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


/** create pairs of sfts */
/*
void CreateSFTIndexPairs(LALStatus                *status,
			 INT4VectorSequence       *out,
			 SFTVector           *inputSFTs,
			 SFTPairParams            *par)
{  


  UINT4 numpairs;
  UINT4 i;
  INITSTATUS (status, "CreateSFTIndexPairs", rcsid);
  ATTATCHSTATUSPTR (status);


  numpairs = inputSFTs->length; 

  TRY( CreateSFTPairsIndicesFrom2SFTvectors(status->statusPtr, out, inputSFTs, 
					    par), status);


  DETATCHSTATUSPTR (status);
	
   normal exit 	
  RETURN (status);
} CreateSFTIndexPairs 
*/

void CreateSFTPairsIndicesFrom2SFTvectors(LALStatus                *status,
					 INT4VectorSequence        **out,
					 SFTVector           *in,
					 SFTPairParams             *par)
{
  
  UINT4 i, j, numsft, numPairs;
  INT4 thisPair;
  INT4Vector *List1, *List2;
  INT4VectorSequence *ret;
  COMPLEX8FrequencySeries  *thisSFT1, *thisSFT2;	       


  INITSTATUS (status, "CreateSFTPairsIndicesFrom2SFTvectors", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (in, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (par, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);

  /* number of SFTs */
  numsft = in->length;
  
  List1 = XLALCreateINT4Vector(numsft*numsft);
  List2 = XLALCreateINT4Vector(numsft*numsft);

  numPairs = 0; /* initialization */

  /* increase length of sftpair vector */

  thisPair = 0;

  
  /* go over all sfts */
  for (i=0; i<numsft; i++) {

    thisSFT1 = in->data + i;
    
    /* go over all sfts again and check if it should be paired with thisSFT1 */
    /* this can be made more efficient if necessary */
    for (j=i; j<numsft; j++) {

      LIGOTimeGPS t1, t2;
      REAL8 timeDiff;

      thisSFT2 = in->data + j;

      /* calculate time difference */      
      t1 = thisSFT1->epoch;
      t2 = thisSFT2->epoch;
      timeDiff = XLALGPSDiff( &t1, &t2);

      /* decide whether to add this pair or not */
      if ( fabs(timeDiff) < par->lag ) {

	numPairs++;

	List1->data[thisPair] = i;
	List2->data[thisPair++] = j;


      } /* if ( numPairs < out->length)  */
    } /* end loop over second sft set */
  } /* end loop over first sft set */ 



  /* initialise pair list vector*/
  ret = (INT4VectorSequence *) LALCalloc(1, sizeof(INT4VectorSequence));
  ret->length = 2;
  ret->vectorLength = numPairs;
  ret->data = LALCalloc(ret->length * ret->vectorLength, sizeof(INT4));

  for (i=0; i < numPairs; i++) {

  ret->data[i] = List1->data[i];
  ret->data[i+numPairs] = List2->data[i];
  }

  (*out) = ret;

  XLALDestroyINT4Vector(List1);
  XLALDestroyINT4Vector(List2);

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

} /* CreateSFTPairsIndicesFrom2SFTvectors */



/*
\** create pairs of sfts *\/
void CreateSFTPairs(LALStatus                *status,
		    SFTPairVec               *out,
		    SFTVector           *inputSFTs,
		    PSDVector	     *inputPSDs,
		    DetectorStateSeries *mdetStates,
		    SFTPairParams            *par)
{  

  UINT4 numifo;

  INITSTATUS (status, "CreateSFTPairs", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (out, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (inputSFTs, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (inputPSDs, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (mdetStates, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (inputSFTs->length == mdetStates->length, status, RADIOMETER_EBAD, RADIOMETER_MSGEBAD);

  numifo = inputSFTs->length;
  \* for now require exactly 2 ifos -- to be relaxed very soon in the future *\/
  if ( numifo != 2) {
    ABORT ( status, RADIOMETER_EBAD, RADIOMETER_MSGEBAD );
  }

  TRY( CreateSFTPairsFrom2SFTvectors(status->statusPtr, out, inputSFTs, 
				     inputPSDs, 
				     mdetStates, 
				     par), status);

  DETATCHSTATUSPTR (status);
	
  \* normal exit *\/	
  RETURN (status);
}





\** create pairs of sfts from a pair of sft vectors*\/
void CreateSFTPairsFrom2SFTvectors(LALStatus                 *status,
				   SFTPairVec                *out,
				   const SFTVector           *in1,
				   const PSDVector	     *psdin1,
				   const DetectorStateSeries *det1,
				   SFTPairParams             *par)
{
  
  UINT4 i, j, numsft1, numsft2, numPairs, numsftMin, counter=0;

  COMPLEX8FrequencySeries  *thisSFT1, *thisSFT2;
  REAL8FrequencySeries 	*thisPSD1, *thisPSD2;	       
  DetectorState *thisDetState1, *thisDetState2;
  SingleSFTpair *thisPair = NULL;

  INITSTATUS (status, "CreateSFTPairsFrom2SFTvectors", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (out, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (in1, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (in2, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (psdin1, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (psdin2, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);

  ASSERT (det1, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (det2, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (par, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);

  \* number of SFTs from the two ifos *\/
  numsft1 = in1->length;
  numsft2 = in2->length;
  numsftMin = (numsft1 < numsft2)? numsft1 : numsft2;

  numPairs = out->length; \* initialization *\/

  \* increase length of sftpair vector *\/
  out->length += numsftMin;
  out->data = LALRealloc(out->data, out->length * sizeof(out->data[0]));

  \* need to allocate memory for singleSFTpair *\/
  thisPair = (SingleSFTpair *) LALCalloc(1, sizeof(out->data[0]));

printf("first detector start time %i\n", det1->data[0].tGPS.gpsSeconds);
printf("second detector start time %i\n",det2->data->tGPS.gpsSeconds);
printf("cf sft epoch %i\n",in1->data[0].epoch.gpsSeconds);

  \* go over all sfts in first vector *\/
  for (i=0; i<numsft1; i++) {
    thisSFT1 = in1->data + i;
    thisPSD1 = psdin1->data + i;
    thisDetState1 = det1->data + i;

    
    \* go over all sfts in second vector and check if it should be paired with thisSFT1 *\/
    for (j=0; j<numsft2; j++) {

      LIGOTimeGPS t1, t2;
      REAL8 timeDiff;

      thisSFT2 = in2->data + j;
      thisPSD2 = psdin2->data + j;
      thisDetState2 = det2->data + j;

      \* calculate time difference *\/      
      t1 = thisSFT1->epoch;
      t2 = thisSFT2->epoch;
      timeDiff = XLALGPSDiff( &t1, &t2);

      \* decide whether to add this pair or not *\/
      if ( fabs(timeDiff) < par->lag ) {
	numPairs++;


	if ( numPairs < out->length) {
	  \* there is enough memory *\/

	  TRY( FillSFTPair( status->statusPtr, thisPair, thisSFT1, thisSFT2, thisPSD1, thisPSD2, thisDetState1, thisDetState2), status);
	out->data[counter] = *thisPair;
	counter++;
	} 
	else {
	  \* there not enough memory -- allocate memory and add the pair *\/

	  out->length += numsftMin;
	  out->data = LALRealloc(out->data, out->length * sizeof(out->data[0]));


	  TRY( FillSFTPair( status->statusPtr, thisPair, thisSFT1, thisSFT2, thisPSD1, thisPSD2, thisDetState1, thisDetState2), status);
	out->data[counter] = *thisPair;
	counter++;
	}

      }
    } \* end loop over second sft set *\/
  } \* end loop over first sft set *\/ 


  \* realloc memory *\/
  out->length = numPairs;
  out->data = LALRealloc(out->data, out->length * sizeof(out->data[0]));

  DETATCHSTATUSPTR (status);
	
  \* normal exit *\/	
  RETURN (status);

}





\* little helper function for filling up sft pair *\/
void FillSFTPair(LALStatus                 *status,
		 SingleSFTpair             *out,
		 COMPLEX8FrequencySeries   *sft1, 
		 COMPLEX8FrequencySeries   *sft2, 
		 REAL8FrequencySeries	   *psd1,
       		 REAL8FrequencySeries      *psd2,
		 DetectorState             *det1,
		 DetectorState             *det2)
{  
  INITSTATUS (status, "FillSFTPairs", rcsid);
  ATTATCHSTATUSPTR (status);
    

  out->sft1 = sft1;
  out->sft2 = sft2;

  out->psd1 = psd1;
  out->psd2 = psd2;

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
   normal exit 	
  RETURN (status);

}
*/

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
  bin1 = (UINT8)( (*freq1 - sft1->f0)*(sft1->data->length) / (sft1->data->length * deltaF));
  bin2 = (UINT8)( (*freq2 - sft2->f0)*(sft2->data->length)/ (sft2->data->length * deltaF));

/*printf("f01 %f, bin location %f, total bins %i\n", sft1->f0, *freq1 - sft1->f0, sft1->data->length);*/

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
  timeDiff = XLALGPSDiff( &(sft1->epoch), &(dopp->refTime));

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
void GetSignalPhaseInSFT(LALStatus               *status,
			 REAL8                   *out,
			 COMPLEX8FrequencySeries *sft1,
			 PulsarDopplerParams     *dopp,
			 REAL8Vector             *pos)
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

  rDotn /= LAL_C_SI;	/* is this right? */ 

  /* now calculate the phase of the SFT */
  /* phi(t) = phi_0 + 2pi(f_0 t + 0.5 f_1 t^2) + 2pi (f_0 + f_1 t) r.n/c */
  /* this is an approximation... need to change in the future? */

  /* this is the sft reference time  - the pulsar reference time */
  timeDiff = XLALGPSDiff( &(sft1->epoch), &(dopp->refTime));

  fhat = dopp->fkdot[0]; /* initialization */
  phihat = 0.0;

  factor = 1.0;
  for (k = 1;  k < PULSAR_MAX_SPINS; k++) {
    factor *= timeDiff / k;  
    fhat += dopp->fkdot[k] * factor;
    phihat += dopp->fkdot[k-1] * factor;

  }

    factor *= timeDiff / k;
    phihat += dopp->fkdot[k-1] * factor;


  *out = LAL_TWOPI * ( phihat + fhat * rDotn);


  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);
}

/** Calculate pair weights (G_alpha) **/
void CalculateUalpha(LALStatus *status,
		      COMPLEX16	*out,
		      REAL8	*Aplus,
		      REAL8	*Across,
		      REAL8	*phiI,
		      REAL8	*phiJ,
		      REAL8	*FplusI,
		      REAL8	*FplusJ,
		      REAL8	*FcrossI,
		      REAL8	*FcrossJ,
		      REAL8 	*freq1,
		      REAL8	*freq2,
		      REAL8FrequencySeries *psd1,
		      REAL8FrequencySeries *psd2)
{
  REAL8 deltaPhi;
  UINT8 bin1, bin2;
  REAL8 deltaF;
  REAL8 re, im, sigmasq;
    
  INITSTATUS (status, "CalculateUalpha", rcsid);
  ATTATCHSTATUSPTR (status);

  deltaF = psd1->deltaF;
  deltaPhi = *phiI - *phiJ;

  re = 0.25 * cos(deltaPhi) * ((*FplusI * (*FplusJ) * (*Aplus) * (*Aplus)) 
			    + (*FcrossI * (*FcrossJ) * (*Across) * (*Across)) );
  im = 0.25 * sin(-deltaPhi) * ((*FplusI * (*FcrossJ) - (*FcrossI) * (*FplusJ)) 
				 * (*Aplus) * (*Across) );

  bin1 = (UINT8)( (*freq1 - psd1->f0)*(psd1->data->length) / (psd1->data->length * deltaF));
  bin2 = (UINT8)( (*freq2 - psd2->f0)*(psd2->data->length)/ (psd2->data->length * deltaF));

  sigmasq = 0.25 * deltaF * deltaF * psd1->data->data[bin1] * psd2->data->data[bin2];

  out->re = re/sigmasq;
  out->im = -im/sigmasq;


  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}

void CalculateWeights(LALStatus       *status,
		      REAL8Vector     *out,
		      COMPLEX16Vector *yalpha,
		      COMPLEX16Vector *ualpha)
{
  INT4 i;

  INITSTATUS (status, "CalculateWeights", rcsid);
  ATTATCHSTATUSPTR (status);

  ASSERT (out, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
  ASSERT (yalpha, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
/*  ASSERT (ualpha, status, RADIOMETER_ENULL, RADIOMETER_MSGENULL);
*/

  for (i=0; i < (INT4)out->length; i++) {

  out->data[i] = 2.0 * (yalpha->data[i].re * ualpha->data[i].re + yalpha->data[i].im * ualpha->data[i].im);

  }

  DETATCHSTATUSPTR (status);
	
  /* normal exit */	
  RETURN (status);

}


