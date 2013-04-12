/*
 *  Copyright (C) 2007 Badri Krishnan
 *  Copyright (C) 2008 Christine Chung, Badri Krishnan and John Whelan
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

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/DopplerScan.h>
#include <lal/PulsarCrossCorr.h>
#include <gsl/gsl_permutation.h>

#define SQUARE(x) (x*x)

/* this function is currently unused (20/07/09) */
void LALCreateSFTPairsIndicesFrom2SFTvectors(LALStatus          *status,
					     INT4VectorSequence **out,
					     SFTListElement     *in,
					     REAL8	        lag,
					     INT4		listLength,
					     DetChoice 		detChoice,
					     BOOLEAN	 	autoCorrelate)
{

  INT4 i, j, numPairs, initj;
  INT4 thisPair;
  INT4 sameDet = 0;
  INT4Vector *List1, *List2;
  INT4VectorSequence *ret;
  COMPLEX8FrequencySeries  *thisSFT1, *thisSFT2;
  LIGOTimeGPS t1, t2;
  REAL8 timeDiff;
  SFTListElement *sfttmp;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (in, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);

  *out = NULL;


  List1 = XLALCreateINT4Vector(listLength);
  List2 = XLALCreateINT4Vector(listLength);

  numPairs = 0; /* initialization */

  thisPair = 0;
  sfttmp = in;

  /*just need to check the first sft against the others */
  thisSFT1 = &(sfttmp->sft);
  i = 0;

  /*if we're autocorrelating, then j starts from 0, otherwise it starts from 1*/
  if (autoCorrelate) {
    initj = 0;
  } else {
    initj = 1;
  }

  for (j=initj; j<listLength; j++) {

      if (j > 0) {
	sfttmp = (SFTListElement *)sfttmp->nextSFT;
      }

      thisSFT2 = &(sfttmp->sft);
      /* calculate time difference */
      t1 = thisSFT1->epoch;
      t2 = thisSFT2->epoch;
      timeDiff = XLALGPSDiff( &t1, &t2);

      /*strcmp returns 0 if strings are equal, >0 if strings are different*/
      sameDet = strcmp(thisSFT1->name, thisSFT2->name);

      /* if they are different, set sameDet to 1 so that it will match if
	 detChoice == DIFFERENT */
      if (sameDet != 0) { sameDet = 1; }

      /* however, if detChoice = ALL, then we want sameDet to match it */
      if (detChoice == ALL) { sameDet = detChoice; }

      /* decide whether to add this pair or not */
      if ((sameDet == (INT4)detChoice) && (fabs(timeDiff) <= lag) ) {
	numPairs++;

	List1->data[thisPair] = i;
	List2->data[thisPair++] = j;
      } /* if ( numPairs < out->length)  */

      thisSFT2 = NULL;
    } /* end loop over second sft set */

    thisSFT1 = NULL;
 /* } end loop over first sft set */


  /* initialise pair list vector, if we have pairs*/
  if (numPairs > 0) {
    ret = XLALCreateINT4VectorSequence(2, numPairs);

    for (i=0; i < numPairs; i++) {

      ret->data[i] = List1->data[i];
      ret->data[i+numPairs] = List2->data[i];
    }
  }
  else {
    ret = NULL;
  }

  (*out) = ret;

  XLALDestroyINT4Vector(List1);
  XLALDestroyINT4Vector(List2);
  sfttmp = NULL;

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);

} /* CreateSFTPairsIndicesFrom2SFTvectors */


/** Correlate a single pair of SFT at a parameter space point. This function calculates Y_alpha
 *  according to Eqn 4.1 in Dhurandar et al 2008, where
 *  Y_alpha = (xI* xJ)/Delta T^2
 *  sft1 and sft2 have been normalised by LALNormalizeSFT in pulsar_crosscorr.c, so they are actually
 *  sft1 = xI/sqrt(psd1), sft2 = xJ/sqrt(psd2)
 *  Therefore, when calculating the output, we need to have
 *  out = sft1*sqrt(psd1)*sft2*sqrt(psd2)/Delta T^2 */
void LALCorrelateSingleSFTPair(LALStatus                *status,
			       COMPLEX16                *out,
			       COMPLEX8FrequencySeries  *sft1,
			       COMPLEX8FrequencySeries  *sft2,
			       REAL8FrequencySeries     *psd1,
			       REAL8FrequencySeries     *psd2,
			       REAL8                    freq1,
			       REAL8                    freq2)
{
  INT4 bin1, bin2;
  REAL8 deltaF;
  REAL8 re1, re2, im1, im2;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  ASSERT (sft1, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT (sft2, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT (psd1, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT (psd2, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);

  /* assume both sfts have the same freq. resolution */
  deltaF = sft1->deltaF;

  /* check that frequencies are in the right range */
  if ((freq1 < sft1->f0) || (freq1 > (sft1->f0 + deltaF*sft1->data->length)) ) {
	ABORT(status, PULSARCROSSCORR_EVAL, PULSARCROSSCORR_MSGEVAL);
  }
  /* check that frequencies are in the right range */
  if ((freq2 < sft2->f0) || (freq2 > (sft2->f0 + deltaF*sft2->data->length)) ) {
	ABORT(status, PULSARCROSSCORR_EVAL, PULSARCROSSCORR_MSGEVAL);
  }



  bin1 = (INT4)ceil( ((freq1 - sft1->f0) / (deltaF)) - 0.5);
  bin2 = (INT4)ceil( ((freq2 - sft2->f0)/ (deltaF)) - 0.5);

  re1 = crealf(sft1->data->data[bin1]);
  im1 = cimagf(sft1->data->data[bin1]);
  re2 = crealf(sft2->data->data[bin2]);
  im2 = cimagf(sft2->data->data[bin2]);

  out->re = (deltaF * deltaF * sqrt(psd1->data->data[bin1] * psd2->data->data[bin2])) * (re1*re2 + im1*im2);
  out->im = (deltaF * deltaF * sqrt(psd1->data->data[bin1] * psd2->data->data[bin2])) * (re1*im2 - re2*im1);

/*
printf("psd1 psd2 %1.15g %1.15g\n", psd1->data->data[bin1], psd2->data->data[bin2]);
printf("sft1.re sft1.im %1.15g %1.15g\n", re1*sqrt(psd1->data->data[bin1]), im1*sqrt(psd1->data->data[bin1]));
printf("sft2.re sft2.im %1.15g %1.15g\n\n", re2*sqrt(psd2->data->data[bin2]), im2*sqrt(psd2->data->data[bin2]));
*/
  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);

}




/** Calculate the frequency of the SFT at a given epoch */
/* This is according to Eqns 2.11 and 2.12 of Dhurandhar et al 2008 */
void LALGetSignalFrequencyInSFT(LALStatus                *status,
				REAL8                    *out,
				LIGOTimeGPS		 *epoch,
				PulsarDopplerParams      *dopp,
				REAL8Vector              *vel_c)
{
  UINT4 k;
  REAL8 alpha, delta;
  REAL8 vDotn_c, fhat, factor, timeDiff;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT(epoch, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT(dopp, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT(vel_c, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);

  alpha = dopp->Alpha;
  delta = dopp->Delta;


  vDotn_c = cos(delta) * cos(alpha) * vel_c->data[0]
    + cos(delta) * sin(alpha)  * vel_c->data[1]
    + sin(delta) * vel_c->data[2];

  /* now calculate the intrinsic signal frequency in the SFT */
  /* fhat = f_0 + f_1(t-t0) + f_2(t-t0)^2/2 + ... */

  /* this is the sft reference time  - the pulsar reference time */
  timeDiff = XLALGPSDiff( (epoch), &(dopp->refTime));
  fhat = dopp->fkdot[0]; /* initialization */
  factor = 1.0;
  for (k = 1;  k < PULSAR_MAX_SPINS; k++) {
    factor *= timeDiff / k;
    fhat += dopp->fkdot[k] * factor;
  }
  *out = fhat * (1 + vDotn_c);

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}



/** Get signal phase at a given epoch */
/* This is according to Eqn 2.13 in Dhurandar et al 2008 */
/* and includes the 2nd order term that is left out in the paper */
void LALGetSignalPhaseInSFT(LALStatus               *status,
			    REAL8                   *out,
			    LIGOTimeGPS 	    *epoch,
			    PulsarDopplerParams     *dopp,
			    REAL8Vector             *r_c)
{
  UINT4 k;
  REAL8 alpha, delta;
  REAL8 rDotn_c, phihat, factor, timeDiff;
  REAL8 epoch_plus_rdotn;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT(epoch, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT(dopp, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT(r_c, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);

  alpha = dopp->Alpha;
  delta = dopp->Delta;

  rDotn_c = cos(delta) * cos(alpha) * r_c->data[0]
    + cos(delta) * sin(alpha) * r_c->data[1]
    + sin(delta) * r_c->data[2];


  /* now calculate the phase of the SFT */
  /* phi(t) = phi_0 + 2pi(f_0 t + 0.5 f_1 t^2) + 2pi (f_0 + f_1 t) r.n/c */

  /* this is the sft reference time  - the pulsar reference time */
  /* we need to convert epoch to REAL8 first before adding rdotn
   * because if LIGOTimeGPS only has INT4 accuracy. converting rdotn into LIGOTimeGPS
   * will introduce rounding errors*/ 
  epoch_plus_rdotn = XLALGPSGetREAL8(epoch) + rDotn_c;
  timeDiff = epoch_plus_rdotn - XLALGPSGetREAL8(&(dopp->refTime));

  phihat = 0.0;

  factor = 1.0;

 for (k = 1;  k <= PULSAR_MAX_SPINS; k++) {
    factor *= timeDiff / k;
    phihat += dopp->fkdot[k-1] * factor;
  }

  *out = LAL_TWOPI * ( phihat );
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/* This function calculates sigma_alpha^2 according to Eqn 4.13 in Dhurandar
 * et al 2008, where
 * sigma_alpha^2 = Sn1*Sn2/(4 DeltaT^2)
 * psd1 and psd2 as returned by LALNormalizeSFT are
 * psd1 = DeltaT*Sn1/2, psd2 = DeltaT*Sn2/2
 * so when calculating the output, we need to compute
 * out = psd1*psd2/DeltaT^4
 * where the factor of DeltaT^2/4 is absorbed in psd1 and psd2*/
void LALCalculateSigmaAlphaSq(LALStatus            *status,
			      REAL8                *out,
			      REAL8                freq1,
			      REAL8                freq2,
			      REAL8FrequencySeries *psd1,
			      REAL8FrequencySeries *psd2)
{
  INT8 bin1, bin2;
  REAL8 deltaF;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT(psd1, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT(psd2, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);

  deltaF = psd1->deltaF;

  bin1 = (INT8)ceil( ((freq1 - psd1->f0) / (deltaF)) - 0.5);
  bin2 = (INT8)ceil( ((freq2 - psd2->f0)/ (deltaF)) - 0.5);
  *out = SQUARE(deltaF)*SQUARE(deltaF) * psd1->data->data[bin1] * psd2->data->data[bin2];
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

}

/** Calculate pair weights (U_alpha) for an average over Psi and cos(iota) **/
/* G_IJ is calculated according to Eqn 3.17 in Dhurandhar et al 2008
 * and Ualpha is calculated according to Eqn 4.10 of the same paper */
void LALCalculateAveUalpha(LALStatus *status,
			COMPLEX16 *out,
			REAL8     phiI,
			REAL8     phiJ,
			REAL8 	  freqI,
			REAL8 	  freqJ,
			REAL8 	  deltaF,
			CrossCorrBeamFn beamfnsI,
			CrossCorrBeamFn beamfnsJ,
			REAL8     sigmasq)
{
  REAL8 deltaPhi;
  REAL8 re, im;
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  deltaPhi = phiI - phiJ + LAL_PI*(freqI - freqJ)/deltaF;
  /*calculate G_IJ. In this case, we have <G_IJ> = 0.1*(exp^(-i delta phi)) * (aIaJ + bIbJ)*/
  re = 0.1 * cos(deltaPhi) * ((beamfnsI.a * beamfnsJ.a) + (beamfnsI.b * beamfnsJ.b));
  im = - 0.1 * sin(deltaPhi) * ((beamfnsI.a * beamfnsJ.a) + (beamfnsI.b * beamfnsJ.b));

  /*calculate Ualpha*/
  out->re = re/(sigmasq);
  out->im = -im/(sigmasq);

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);


}
/** Calculate pair weights (U_alpha) for the general case **/
/* G_IJ is calculated according to Eqn 3.10 of Dhurandar et al 2008.
 * deltaPhi is calculated from Eqn 3.11, and Ualpha is calculated according to
 * Eqn 4.10 of the same paper.
 * The Fplus, Fcross terms are calculated according to Eqns 2.9a and 2.9b */
void LALCalculateUalpha(LALStatus *status,
			COMPLEX16 *out,
			CrossCorrAmps amplitudes,
			REAL8     phiI,
			REAL8     phiJ,
			REAL8 	  freqI,
			REAL8 	  freqJ,
			REAL8 	  deltaF,
			CrossCorrBeamFn beamfnsI,
			CrossCorrBeamFn beamfnsJ,
			REAL8     sigmasq,
			REAL8     *psi,
			COMPLEX16 *gplus,
			COMPLEX16 *gcross)
{
  REAL8 deltaPhi;
  REAL8 re, im;
  REAL8 FplusI, FplusJ, FcrossI, FcrossJ;
  REAL8 FpIFpJ, FcIFcJ, FpIFcJ, FcIFpJ;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  deltaPhi = phiI - phiJ + LAL_PI*(freqI - freqJ)/deltaF;

 
  /*if not averaging over psi, calculate F+, Fx exactly*/
  if (psi) {
    FplusI = (beamfnsI.a * cos(2.0*(*psi))) + (beamfnsI.b * sin(2.0*(*psi)));
    FplusJ = (beamfnsJ.a * cos(2.0*(*psi))) + (beamfnsJ.b * sin(2.0*(*psi)));
    FcrossI = (beamfnsI.b * cos(2.0 * (*psi))) - (beamfnsI.a * sin(2.0 * (*psi)));
    FcrossJ = (beamfnsJ.b * cos(2.0 * (*psi))) - (beamfnsJ.a * sin(2.0 * (*psi)));;

    /*calculate G_IJ*/
    re = 0.25 * ( (cos(deltaPhi)* (FplusI*FplusJ * amplitudes.Aplussq + FcrossI*FcrossJ * amplitudes.Acrosssq) )
		- (sin(deltaPhi)*((FplusI*FcrossJ - FcrossI*FplusJ) * amplitudes.AplusAcross)) );


    im = 0.25 * ( -(cos(deltaPhi)*((FplusI*FcrossJ - FcrossI*FplusJ)*amplitudes.AplusAcross))
	          - (sin(deltaPhi) * (FplusI*FplusJ*amplitudes.Aplussq + FcrossI*FcrossJ * amplitudes.Acrosssq)) );

  /*calculate estimators*/
  gplus->re = 0.25*cos(deltaPhi)*FplusI*FplusJ;
  gplus->im = 0.25*(-sin(deltaPhi))*FplusI*FplusJ;

  gcross->re = 0.25*cos(deltaPhi)*FcrossI*FcrossJ;
  gcross->im = 0.25*(-sin(deltaPhi))*FcrossI*FcrossJ;


  }
  else {
    FpIFpJ = 0.5*((beamfnsI.a * beamfnsJ.a) + (beamfnsI.b * beamfnsJ.b));
    FcIFcJ = 0.5*((beamfnsI.a * beamfnsJ.a) + (beamfnsI.b * beamfnsJ.b));
    FpIFcJ = 0.5*((beamfnsI.a * beamfnsJ.b) - (beamfnsI.b * beamfnsJ.a));
    FcIFpJ = 0.5*((beamfnsJ.a * beamfnsI.b) - (beamfnsJ.b * beamfnsI.a));

    /*calculate G_IJ*/
    re = 0.25 * ( (cos(deltaPhi)*
		   ((FpIFpJ * amplitudes.Aplussq) + (FcIFcJ * amplitudes.Acrosssq)) )
  		  -(sin(deltaPhi)*((FpIFcJ - FcIFpJ) * amplitudes.AplusAcross)) );


    im = 0.25 * ( -(cos(deltaPhi) * ((FpIFcJ - FcIFpJ)*amplitudes.AplusAcross))
	       - (sin(deltaPhi) *
		   ((FpIFpJ * amplitudes.Aplussq) + (FcIFcJ * amplitudes.Acrosssq))) );

  }


  /*calculate Ualpha*/
  out->re = re/(sigmasq);
  out->im = -im/(sigmasq);



  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);

}


/* Calculate rho, the cross correlation power,
 * according to Eqn 4.4 in Dhurandhar et al 2008 */
void LALCalculateCrossCorrPower(LALStatus       *status,
				REAL8	        *out,
				COMPLEX16Vector *yalpha,
				COMPLEX16Vector *ualpha)
{
  INT4 i;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (out, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT (yalpha, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT (ualpha, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);

  *out = 0;

  for (i=0; i < (INT4)yalpha->length; i++) {

  *out += 2.0 * ((yalpha->data[i].re * ualpha->data[i].re) - (yalpha->data[i].im * ualpha->data[i].im));

  }

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);

}

/* Calculate the variance of the cross correlation power.
 * The variance is given by sigma^2 in Eqn 4.6 of Dhurandhar et al 2008 */
void LALNormaliseCrossCorrPower(LALStatus        *status,
				REAL8		 *out,
				COMPLEX16Vector  *ualpha,
				REAL8Vector      *sigmaAlphasq)
{
  INT4 i;
  REAL8 variance = 0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (out, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT (ualpha, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT (sigmaAlphasq, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);


  for (i=0; i < (INT4)ualpha->length; i++) {
	variance += (SQUARE(ualpha->data[i].re) + SQUARE(ualpha->data[i].im)) * sigmaAlphasq->data[i];

  }

  variance *= 2.0;
  *out = variance;
  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);


}

/* Calculate the estimators for Aplus, Asq.
 * Eqn 6.4 of Dhurandhar et al 2008 */
void LALCalculateEstimators(LALStatus    *status,
				REAL8 *aplussq1,
				REAL8 *aplussq2,
				REAL8 *acrossq1,
				REAL8 *acrossq2,
				COMPLEX16Vector  *yalpha,
				COMPLEX16Vector  *gplus,
				COMPLEX16Vector  *gcross,
				REAL8Vector      *sigmaAlphasq)
{
  INT4 i;
  REAL8 ap1 = 0, ap2 = 0, ac1 = 0, ac2 = 0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  ASSERT (aplussq1, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT (yalpha, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);
  ASSERT (sigmaAlphasq, status, PULSARCROSSCORR_ENULL, PULSARCROSSCORR_MSGENULL);


  for (i=0; i < (INT4)yalpha->length; i++) {
	ap1 += 2.0*(SQUARE(gplus->data[i].re) + SQUARE(gplus->data[i].im))/sigmaAlphasq->data[i];
	ac1 += 2.0*(SQUARE(gcross->data[i].re) + SQUARE(gcross->data[i].im))/sigmaAlphasq->data[i];
	ap2 += 2.0*((yalpha->data[i].re * gplus->data[i].re) + (yalpha->data[i].im * gplus->data[i].im))/sigmaAlphasq->data[i];
	ac2 += 2.0*((yalpha->data[i].re * gcross->data[i].re) + (yalpha->data[i].im * gcross->data[i].im))/sigmaAlphasq->data[i];
  }

  *aplussq1 = ap1;
  *aplussq2 = ap2;
  *acrossq1 = ac1;
  *acrossq2 = ac2;
  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);


}

