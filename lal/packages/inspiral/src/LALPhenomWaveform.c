/*
*  Copyright (C) 2008 P. Ajith, Badri Krishnan, Lucia Santamaria
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


#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>
#include <lal/SeqFactories.h>
#include <lal/NRWaveInject.h>
#include <lal/RealFFT.h>
#include <lal/TimeFreqFFT.h>
#include <lal/SeqFactories.h>
#include <lal/VectorOps.h>
#include <lal/BBHPhenomCoeffs.h>
#include <lal/AVFactories.h>

#include <math.h>

typedef struct tagBBHPhenomParams{
  REAL8 fMerger;
  REAL8 fRing;
  REAL8 fCut;
  REAL8 sigma;
  REAL8 psi0;
  REAL8 psi1;
  REAL8 psi2;
  REAL8 psi3;
  REAL8 psi4;
  REAL8 psi5;
  REAL8 psi6;
  REAL8 psi7;
}
BBHPhenomParams;


static void XLALComputePhenomParams( BBHPhenomParams *phenParams,
			      InspiralTemplate *params);

static REAL8 XLALLorentzianFn ( REAL8 freq,
			 REAL8 fRing,
			 REAL8 sigma);


static void XLALBBHPhenWaveFD ( BBHPhenomParams  *params,
			 InspiralTemplate *insp_template,
			 REAL4Vector *signalvec);


static void XLALComputeInstantFreq( REAL4Vector *Freq,
			     REAL4Vector *hp,
			     REAL4Vector *hc,
			     REAL8 dt);

static REAL4Vector *XLALCutAtFreq( REAL4Vector     *h,
				       REAL4Vector     *freq,
				       REAL8           cutFreq,
				       REAL8           deltaT);

NRCSID (LALPHENOMWAVEFORMC, "$Id: LALPhenomWaveform.c,v 1.13 2009/04/01 00:54:59 ajith Exp $");


/*********************************************************************/
/*** Top level function to generate the phenomenological waveform ****/
/*********************************************************************/
void LALBBHPhenWaveFreqDom ( LALStatus        *status,
			     REAL4Vector      *signalvec,
			     InspiralTemplate *params)
{

  BBHPhenomParams phenParams;

  /* check inputs */
  INITSTATUS (status, "LALBBHPhenWaveFreqDom", LALPHENOMWAVEFORMC);
  ATTATCHSTATUSPTR(status);
  ASSERT (signalvec,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (signalvec->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (signalvec->length>2,  status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);

    /* compute the phenomenological parameters */
    XLALComputePhenomParams(&phenParams, params);

    /* generate the phenomenological waveform in frequency domain */
    XLALBBHPhenWaveFD (&phenParams, params, signalvec);

/*     params->fFinal = params->fCut; */

    DETATCHSTATUSPTR(status);
    RETURN(status);
}


static void XLALComputePhenomParams( BBHPhenomParams  *phenParams,
			      InspiralTemplate *params)

{

  REAL8 totalMass, piM, eta, fMerg_a, fMerg_b, fMerg_c, fRing_a, fRing_b;
  REAL8 fRing_c, sigma_a, sigma_b, sigma_c, fCut_a, fCut_b, fCut_c;
  REAL8 psi0_a, psi0_b, psi0_c, psi2_a, psi2_b, psi2_c, psi3_a, psi3_b, psi3_c;
  REAL8 psi4_a, psi4_b, psi4_c, psi6_a, psi6_b, psi6_c, psi7_a, psi7_b, psi7_c;

  /* calculate the total mass and symmetric mass ratio */

  if (params) {

    totalMass = params->mass1+params->mass2;
    eta = params->mass1*params->mass2/pow(totalMass,2.);
    piM = totalMass*LAL_PI*LAL_MTSUN_SI;
  }
  else {
    return;
  }

  fMerg_a = BBHPHENOMCOEFFSH_FMERG_A;
  fMerg_b = BBHPHENOMCOEFFSH_FMERG_B;
  fMerg_c = BBHPHENOMCOEFFSH_FMERG_C;

  fRing_a = BBHPHENOMCOEFFSH_FRING_A;
  fRing_b = BBHPHENOMCOEFFSH_FRING_B;
  fRing_c = BBHPHENOMCOEFFSH_FRING_C;

  sigma_a = BBHPHENOMCOEFFSH_SIGMA_A;
  sigma_b = BBHPHENOMCOEFFSH_SIGMA_B;
  sigma_c = BBHPHENOMCOEFFSH_SIGMA_C;

  fCut_a = BBHPHENOMCOEFFSH_FCUT_A;
  fCut_b = BBHPHENOMCOEFFSH_FCUT_B;
  fCut_c = BBHPHENOMCOEFFSH_FCUT_C;

  psi0_a = BBHPHENOMCOEFFSH_PSI0_X;
  psi0_b = BBHPHENOMCOEFFSH_PSI0_Y;
  psi0_c = BBHPHENOMCOEFFSH_PSI0_Z;

  psi2_a = BBHPHENOMCOEFFSH_PSI2_X;
  psi2_b = BBHPHENOMCOEFFSH_PSI2_Y;
  psi2_c = BBHPHENOMCOEFFSH_PSI2_Z;

  psi3_a = BBHPHENOMCOEFFSH_PSI3_X;
  psi3_b = BBHPHENOMCOEFFSH_PSI3_Y;
  psi3_c = BBHPHENOMCOEFFSH_PSI3_Z;

  psi4_a = BBHPHENOMCOEFFSH_PSI4_X;
  psi4_b = BBHPHENOMCOEFFSH_PSI4_Y;
  psi4_c = BBHPHENOMCOEFFSH_PSI4_Z;

  psi6_a = BBHPHENOMCOEFFSH_PSI6_X;
  psi6_b = BBHPHENOMCOEFFSH_PSI6_Y;
  psi6_c = BBHPHENOMCOEFFSH_PSI6_Z;

  psi7_a = BBHPHENOMCOEFFSH_PSI7_X;
  psi7_b = BBHPHENOMCOEFFSH_PSI7_Y;
  psi7_c = BBHPHENOMCOEFFSH_PSI7_Z;

  /* Evaluate the polynomials. See Eq. (4.18) of P. Ajith et al
   * arXiv:0710.2335 [gr-qc] */
  if (phenParams) {

    phenParams->fCut  = (fCut_a*eta*eta  + fCut_b*eta  + fCut_c)/piM;
    phenParams->fMerger  = (fMerg_a*eta*eta  + fMerg_b*eta  + fMerg_c)/piM;
    phenParams->fRing  = (fRing_a*eta*eta + fRing_b*eta + fRing_c)/piM;
    phenParams->sigma = (sigma_a*eta*eta + sigma_b*eta + sigma_c)/piM;

    phenParams->psi0 = (psi0_a*eta*eta + psi0_b*eta + psi0_c)/(eta*pow(piM, 5./3.));
    phenParams->psi1 = 0.;
    phenParams->psi2 = (psi2_a*eta*eta + psi2_b*eta + psi2_c)/(eta*pow(piM, 3./3.));
    phenParams->psi3 = (psi3_a*eta*eta + psi3_b*eta + psi3_c)/(eta*pow(piM, 2./3.));
    phenParams->psi4 = (psi4_a*eta*eta + psi4_b*eta + psi4_c)/(eta*pow(piM, 1./3.));
    phenParams->psi5 = 0.;
    phenParams->psi6 = (psi6_a*eta*eta + psi6_b*eta + psi6_c)/(eta*pow(piM, -1./3.));
    phenParams->psi7 = (psi7_a*eta*eta + psi7_b*eta + psi7_c)/(eta*pow(piM, -2./3.));
  }

  return;

}


static void XLALBBHPhenWaveFD ( BBHPhenomParams  *params,
				InspiralTemplate *insp_template,
				REAL4Vector *signalvec) {

    REAL8 df, shft, phi, amp0, ampEff, psiEff, fMerg, fNorm;
    REAL8 f, fRing, sigma, totalMass, eta;
    INT4 i, j, n;

    /* freq resolution and the low-freq bin */
    df = insp_template->tSampling/signalvec->length;
    n = signalvec->length;

    /* If we want to pad with zeroes in the beginning then the instant of
    * coalescence will be the chirp time + the duration for which padding
    * is needed. Thus, in the equation below nStartPad occurs with a +ve sign.
    * This code doesn't support non-zero start-time. i.e. params->startTime
    * should be necessarily zero.*/
    shft = 2.*LAL_PI * ((REAL4)signalvec->length/insp_template->tSampling +
            insp_template->nStartPad/insp_template->tSampling + insp_template->startTime);
    phi  = insp_template->startPhase;

    /* phenomenological  parameters*/
    fMerg = params->fMerger;
    fRing = params->fRing;
    sigma = params->sigma;
    totalMass = insp_template->mass1 + insp_template->mass2;
    eta = insp_template->mass1 * insp_template->mass2 / pow(totalMass, 2.);

    /* Now compute the amplitude.  NOTE the params->distance is assumed to
     * me in meters. This is, in principle, inconsistent with the LAL
     * documentation (inspiral package). But this seems to be the convention
     * employed in the injection codes */
    amp0 = pow(LAL_MTSUN_SI*totalMass, 5./6.)*pow(fMerg,-7./6.)/pow(LAL_PI,2./3.);
    amp0 *= pow(5.*eta/24., 1./2.)/(insp_template->distance/LAL_C_SI);

    /* fill the zero and Nyquist frequency with zeros */
    *(signalvec->data+0) = 0.;
    *(signalvec->data+n/2) = 0.;

    /* now generate the waveform at all frequency bins */
    for (i=1; i<n/2; i++) {

        /* this is the index of the imaginary part */
        j = n-i;

        /* fourier frequency corresponding to this bin */
      	f = i * df;
        fNorm = f/fMerg;

    	/* compute the amplitude */
        if ((f < insp_template->fLower) || (f > params->fCut)) {
            ampEff = 0.;
        }
        else if (f <= fMerg) {
            ampEff = amp0*pow(fNorm, -7./6.);
        }
        else if ((f > fMerg) & (f <= fRing)) {
            ampEff = amp0*pow(fNorm, -2./3.);
        }
        else if (f > fRing) {
            ampEff = XLALLorentzianFn ( f, fRing, sigma);
            ampEff *= amp0*LAL_PI_2*pow(fRing/fMerg,-2./3.)*sigma;
        }

        /* now compute the phase */
       	psiEff = shft*f + phi
                    + params->psi0*pow(f,-5./3.)
                    + params->psi1*pow(f,-4./3.)
                    + params->psi2*pow(f,-3./3.)
                    + params->psi3*pow(f,-2./3.)
                    + params->psi4*pow(f,-1./3.)
                    + params->psi5*pow(f,0.)
                    + params->psi6*pow(f,1./3.)
                    + params->psi7*pow(f,2./3.);

       	/* generate the waveform */
       	*(signalvec->data+i) = (REAL4) (ampEff * cos(psiEff));     /* real */
        *(signalvec->data+j) = (REAL4) (ampEff * sin(psiEff));    /* imag */
    }

}



static REAL8 XLALLorentzianFn ( REAL8 freq,
			 REAL8 fRing,
			 REAL8 sigma)
{
  REAL8 out;

  out = sigma / (2 * LAL_PI * ((freq - fRing)*(freq - fRing)
			       + sigma*sigma / 4.0));

  return(out);
}


void LALBBHPhenWaveFreqDomTemplates( LALStatus        *status,
				     REAL4Vector      *signalvec1,
				     REAL4Vector      *signalvec2,
				     InspiralTemplate *params)
{

  INITSTATUS(status, "LALBBHPhenWaveFreqDomTemplates", LALPHENOMWAVEFORMC);
  ATTATCHSTATUSPTR(status);

  ASSERT(signalvec1, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signalvec2, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signalvec1->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signalvec2->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  /* Initially the waveforms are empty */
  memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));
  memset(signalvec2->data, 0, signalvec2->length * sizeof(REAL4));

  /* generate one waveform with startPhase specified by the user */
  LALBBHPhenWaveFreqDom(status->statusPtr, signalvec1, params);
  CHECKSTATUSPTR(status);

  /* generate another waveform orthogonal to it */
  params->startPhase += LAL_PI_2;
  LALBBHPhenWaveFreqDom(status->statusPtr, signalvec2, params);
  CHECKSTATUSPTR(status);

  DETATCHSTATUSPTR(status);
  RETURN (status);

}


void LALBBHPhenWaveTimeDom ( LALStatus        *status,
			     REAL4Vector      *signalvec,
			     InspiralTemplate *insp_template)
{

  REAL8 fLower;
  REAL8 fCut, fRes, f, totalMass, softWin;
  REAL8 fLowerOrig, eta, tau0;
  REAL8 winFLo, winFHi, sigLo, sigHi;
  REAL4 windowLength;
  INT4 i, k, n;
  /* INT4 kMin, kMax; */
  /* REAL8 dfM;*/
  REAL4Vector *signalFD1 = NULL;
  REAL4FFTPlan *revPlan = NULL;
  /* FILE *filePtr; */
  BBHPhenomParams phenParams;

  /* check inputs */
  INITSTATUS (status, "LALBBHPhenWaveTimeDom", LALPHENOMWAVEFORMC);
  ATTATCHSTATUSPTR(status);
  ASSERT (signalvec,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (signalvec->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (insp_template, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (signalvec->length>2,  status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
  ASSERT (insp_template->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (insp_template->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (insp_template->fLower > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (insp_template->tSampling > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


  /* compute the phenomenological parameters */
  XLALComputePhenomParams(&phenParams, insp_template);

  totalMass = insp_template->mass1 + insp_template->mass2;
  eta = insp_template->mass1*insp_template->mass2/pow(totalMass,2.);

  /* we will generate the waveform from a frequency which is lower than the
   * fLower chosen. Also the cutoff frequency is higher than the fCut. We
   * will later apply a window function, and truncate the time-domain waveform
   * below an instantaneous frequency  fLower */
  fLowerOrig = insp_template->fLower;    /* this is the low-freq set by the user */

    /* Find an optimum value for fLower (using the definition of Newtonian chirp time)
     * such that the waveform has a minimum length of tau0. This is necessary to avoid
     * FFT artifacts */
    tau0 = 32.;
    fLower = pow((tau0*256.*eta*pow(totalMass*LAL_MTSUN_SI,5./3.)/5.),-3./8.)/LAL_PI; /* Actually, this seems to be a better ansatz
                                                                                        as this makes use of the definition of Newtonian
                                                                                        chirp time in order to find an fLower that gives
                                                                                        a minimum length of tau0 (which can be fixed) for
                                                                                        the waveform - Ajith, 18 Nov 08 */
  /* Better ansatz (Lucia Oct 08) */
    /*   fLower = 2.E-3/(totalMass*LAL_MTSUN_SI); */
    /* fLower = 18. - 3.*totalMass/25.; */
  fCut = (1.025)*phenParams.fCut;

  /* make sure that these frequencies are not too out of range */
  if (fLower > fLowerOrig) fLower = fLowerOrig;
  if (fCut > insp_template->tSampling/2.-100.) fCut = insp_template->tSampling/2.-100.;

  /* generate waveforms over this frequency range */
  insp_template->fLower = fLower;
  phenParams.fCut = insp_template->tSampling/2.;

  /* make sure that fLower is not too low */
   if (insp_template->fLower < 0.5) insp_template->fLower = 0.5;

  /* generate the phenomenological waveform in frequency domain */
  n = signalvec->length;
  signalFD1 = XLALCreateREAL4Vector(n);
  XLALBBHPhenWaveFD (&phenParams, insp_template, signalFD1);

  /* apply the softening window function */
  fRes = insp_template->tSampling/n;

  /********************************* DEBUG ********************************/
  /*      filePtr = fopen("FreqDomPhenWave.txt","a");
	  for (i = 1; i < n/2; i++) {
	  fprintf(filePtr,"%e\t%e\t%e\n", i*fRes, signalFD1->data[i],signalFD1->data[n-i]);
	  }
	  fclose(filePtr); */
  /************************************************************************/

  winFLo = (fLowerOrig + fLower)/2.;
  winFHi = (fCut + phenParams.fCut)/2.;
  sigLo = 4.;
  sigHi = 4.;

  signalFD1->data[0] *= 0.;
  for (k = 1; k <= n/2; k++) {
    f = k*fRes;
    softWin = (1+tanh((4.*(f-winFLo)/sigLo)))*(1-tanh(4.*(f-fCut)/sigHi))/4.;
    signalFD1->data[k] *= softWin;
    signalFD1->data[n-k] *= softWin;
    }


  /********************************* DEBUG ********************************/
  /*       filePtr = fopen("FreqDomPhenWave_Wind.txt","a");
	   for (i = 1; i < n/2; i++) {
	   fprintf(filePtr,"%e\t%e\t%e\n", i*fRes, signalFD1->data[i],signalFD1->data[n-i]);
	   }
	   fclose(filePtr); */
  /************************************************************************/

  /* Inverse Fourier transform */
  LALCreateReverseREAL4FFTPlan(status->statusPtr, &revPlan, n, 0 );
  LALREAL4VectorFFT(status->statusPtr, signalvec, signalFD1, revPlan);
  XLALDestroyREAL4Vector(signalFD1);
  LALDestroyREAL4FFTPlan(status->statusPtr, &revPlan);

  /* FFT normalisation. The LAL implementation of the FFT omits the factor 1/n.
   * Also we change the sign of the waveform so that the initialPhase = 0 and pi/2
   * will match to the 'plus' and 'cross' polarisations of the hybrid waveforms,
     * respectively*/
  for (i = 0; i < n; i++) {
    signalvec->data[i] *= -insp_template->tSampling/n;
  }

  /********************************* DEBUG ********************************/
  /*filePtr = fopen("TimeDomPhenWave.txt","a");
    for (i = 0; i < n; i++) {
    fprintf(filePtr,"%e\t%e\n", i/insp_template->tSampling, signalvec->data[i]);
    }
    fclose(filePtr);*/
  /************************************************************************/

  /* apply a linearly increasing/decresing window at the beginning and at the end
   * of the waveform in order to avoid edge effects. This could be made fancier */
   windowLength = 10.*totalMass * LAL_MTSUN_SI*insp_template->tSampling;
   for (i=0; i< windowLength; i++){
         signalvec->data[n-i-1] *= i/windowLength;
   }
   windowLength = 1000.*totalMass * LAL_MTSUN_SI*insp_template->tSampling;
   for (i=0; i< windowLength; i++){
        signalvec->data[i] *= i/windowLength;
   }

  /********************************* DEBUG ********************************/
  /* CHAR fileName[1000];
     sprintf(fileName, "TimeDomPhenWave_Wind_Phi0%4.3f.txt",insp_template->startPhase);
     filePtr = fopen(fileName,"a");
     for (i = 0; i < n; i++) {
     fprintf(filePtr,"%e\t%e\n", i/insp_template->tSampling, signalvec->data[i]);
     }
     fclose(filePtr); */
  /************************************************************************/

    /* reassign the original value of fLower */
  insp_template->fLower = fLowerOrig;
  insp_template->fFinal = phenParams.fCut;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


void LALBBHPhenWaveTimeDomTemplates( LALStatus        *status,
				     REAL4Vector      *signalvec1,
				     REAL4Vector      *signalvec2,
				     InspiralTemplate *params)
{

    UINT4 n, i, peakAmpIdx;
    REAL4Vector *f=NULL, *a=NULL;
    REAL8 peakAmp, dt;

  INITSTATUS(status, "LALBBHPhenWaveTimeDomTemplates", LALPHENOMWAVEFORMC);
  ATTATCHSTATUSPTR(status);

  ASSERT(signalvec1, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signalvec2, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signalvec1->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(signalvec2->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  /* Initially the waveforms are empty */
  memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));
  memset(signalvec2->data, 0, signalvec2->length * sizeof(REAL4));

  /* generate one waveform with startPhase specified by the user */
  LALBBHPhenWaveTimeDom(status->statusPtr, signalvec1, params);
  CHECKSTATUSPTR(status);

  /* generate another waveform orthogonal to it */
  params->startPhase += LAL_PI_2;
  LALBBHPhenWaveTimeDom(status->statusPtr, signalvec2, params);
  CHECKSTATUSPTR(status);

    /* compute the instantaneous frequency */
    dt = 1./params->tSampling;
    n = signalvec1->length;
    f = XLALCreateREAL4Vector(n);
    a = XLALCreateREAL4Vector(n);
    if (f) XLALComputeInstantFreq(f, signalvec1, signalvec2, dt);
    peakAmp = 0.;
    peakAmpIdx = 0;

    /* find the peak amplitude of the waveform */
    for (i=0; i<n; i++){

        a->data[i] = sqrt(pow(signalvec1->data[i],2.) + pow(signalvec2->data[i],2.));

        /* find the peak amplitude*/
        if (a->data[i] > peakAmp) {
            peakAmp = a->data[i];
            peakAmpIdx = i;
        }
    }

    /* if the instantaneous amplitude is less than 1/1000 of the peak amplitude, set the
     * instantaneous freq to be zero. This frequency estimation can very well be corrput due
     * to the very low amplitude of the signal, and is dominated by noise arising from the
     * edge effects */
    for (i=0; i< f->length; i++) {
        if (a->data[i] < 1.0e-3*peakAmp) {
            f->data[i] = 0.0;
        }
    }

    /* cut the waveform at the low freq given by the user */
    signalvec1 = XLALCutAtFreq( signalvec1, f, params->fLower, dt);
    signalvec2 = XLALCutAtFreq( signalvec2, f, params->fLower, dt);

    /* store some paramteters for record keeping */
    params->vFinal = 0.5;           /* this parameter has realy no meaning here*/
    params->tC = peakAmpIdx*dt;     /* time of coalescence. defined as the time
                                       corresponding to the peak amplitude*/

    /* free the memory allocated to f and a */
    XLALDestroyREAL4Vector(f);
    XLALDestroyREAL4Vector(a);

  DETATCHSTATUSPTR(status);
  RETURN (status);

}


void LALBBHPhenTimeDomEngine( LALStatus        *status,
			      REAL4Vector      *signalvec1,
			      REAL4Vector      *signalvec2,
			      REAL4Vector      *h,
			      REAL4Vector      *a,
			      REAL4Vector      *f,
			      REAL8Vector      *phiOut,
			      InspiralTemplate *params)
{

    INT4 i, j, k, n;
    REAL8 dt, cosI;
    REAL8Vector *phi=NULL;

    INITSTATUS(status, "LALBBHPhenTimeDomEngine", LALPHENOMWAVEFORMC);
    ATTATCHSTATUSPTR(status);

    /* check inputs */
    ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    ASSERT (params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    ASSERT (params->fLower > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    ASSERT (params->tSampling > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

    ASSERT(signalvec1, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(signalvec2, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(signalvec1->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    ASSERT(signalvec2->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

    dt = 1./params->tSampling;

    cosI = cos( params->inclination );

    /* generate two orthogonal waveforms */
    LALBBHPhenWaveTimeDomTemplates(status->statusPtr, signalvec1, signalvec2, params);
    CHECKSTATUSPTR(status);

    /* compute the instantaneous frequency */
    if (f) XLALComputeInstantFreq(f, signalvec1, signalvec2, dt);

    /* allocate memory for the temporary phase vector */
    n = signalvec1->length;
    if (phiOut) phi = XLALCreateREAL8Vector(n);

    /* compute the amplitude, phase and frequency.  Fill the polarisation vector h
     * in the prescribed format */
    for (i=0; i<n; i++){
        j = 2*i;
        k = j+1;

        if (phiOut) {
            phi->data[i] = -atan2(signalvec2->data[i], signalvec1->data[i]);
        }

        /* fill the amplitude vector, if required. Currently we assume that both
         * polarisations are of equal ampliude, which is defined as [hp^2+hc^2]^0.5 */
        if (a) {
            a->data[j] = sqrt(pow(signalvec1->data[i],2.) + pow(signalvec2->data[i],2.));
            a->data[k] = sqrt(pow(signalvec1->data[i],2.) + pow(signalvec2->data[i],2.));
        }

        /* fill in the h vector, if required */
        if (h) {
            REAL8 z1 = 0.25 * (( 1. + cosI ) * ( 1. + cosI ) + ( 1. - cosI ) * ( 1. - cosI ));
            REAL8 z2 = 0.25 * (( 1. - cosI ) * ( 1. - cosI ) - ( 1. + cosI ) * ( 1. + cosI ));

            h->data[j] = z1 * signalvec1->data[i];
            h->data[k] = z2 * signalvec2->data[i];
        }
     }

    /* unwrap the phase */
    if (phiOut) {
        LALUnwrapREAL8Angle (status->statusPtr, phiOut, phi);
        XLALDestroyREAL8Vector(phi);
    }

    DETATCHSTATUSPTR(status);
    RETURN (status);
}




void LALBBHPhenWaveTimeDomForInjection (LALStatus        *status,
					CoherentGW       *waveform,
					InspiralTemplate *params,
					PPNParamStruc    *ppnParams) {

  REAL4Vector *a=NULL;      /* amplitude  data */
  REAL4Vector *h=NULL;      /* polarization data */
  REAL4Vector *ff=NULL;     /* frequency data */
  REAL8Vector *phi=NULL;    /* phase data */

  REAL4Vector *hp=NULL;     /* signal with initial phase 0 - temporary */
  REAL4Vector *hc=NULL;     /* signal with initial phase pi/2 - temporary */

  UINT4 count, i;
  REAL8 s, phiC;            /* phase at coalescence */
  CHAR message[256];
  InspiralInit paramsInit;
  CreateVectorSequenceIn in;

  INITSTATUS(status, "LALBBHPhenWaveTimeDomForInjection", LALPHENOMWAVEFORMC);
  ATTATCHSTATUSPTR(status);


  /* Make sure parameter and waveform structures exist. */
  ASSERT( params, status, LALINSPIRALH_ENULL,  LALINSPIRALH_MSGENULL );
  ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  /* check inputs */
  ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (params->fLower > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  ASSERT (params->tSampling > 0.0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

  /* Make sure waveform fields don't exist. */
  ASSERT( !( waveform->a ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->h ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->f ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->phi ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->shift ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );

  params->ampOrder = 0;
  sprintf(message, "WARNING: Amp Order has been reset to %d", params->ampOrder);
  LALInfo(status, message);

  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);

  count = paramsInit.nbins;

  if (paramsInit.nbins==0) {
      DETATCHSTATUSPTR(status);
      RETURN (status);
  }

  /* Now we can allocate memory and vector for coherentGW structure*/
  LALSCreateVector(status->statusPtr, &ff, paramsInit.nbins);
  CHECKSTATUSPTR(status);
  LALSCreateVector(status->statusPtr, &a, 2*paramsInit.nbins);
  CHECKSTATUSPTR(status);
  LALDCreateVector(status->statusPtr, &phi, paramsInit.nbins);
  CHECKSTATUSPTR(status);
  LALSCreateVector(status->statusPtr, &hp, paramsInit.nbins);
  CHECKSTATUSPTR(status);
  LALSCreateVector(status->statusPtr, &hc, paramsInit.nbins);
  CHECKSTATUSPTR(status);

  /* By default the waveform is empty */
  memset(ff->data, 0, paramsInit.nbins * sizeof(REAL4));
  memset(a->data, 0, 2 * paramsInit.nbins * sizeof(REAL4));
  memset(phi->data, 0, paramsInit.nbins * sizeof(REAL8));

  if( params->approximant == IMRPhenomA ) {
      LALSCreateVector(status->statusPtr, &h, 2*paramsInit.nbins);
      CHECKSTATUSPTR(status);
      memset(h->data,  0, 2*paramsInit.nbins * sizeof(REAL4));
  }

  /* generate two orthogonal waveforms */
  params->startPhase = ppnParams->phi;
  LALBBHPhenTimeDomEngine(status->statusPtr, hp, hc, h, a, ff, phi, params);

  BEGINFAIL(status) {
    LALSDestroyVector(status->statusPtr, &ff);
    CHECKSTATUSPTR(status);
    LALSDestroyVector(status->statusPtr, &a);
    CHECKSTATUSPTR(status);
    LALDDestroyVector(status->statusPtr, &phi);
    CHECKSTATUSPTR(status);
    if( params->approximant == IMRPhenomA ){
      LALSDestroyVector(status->statusPtr, &h);
      CHECKSTATUSPTR(status);
    }
  }
  ENDFAIL(status);

  /* Check an empty waveform hasn't been returned */
  for (i = 0; i < phi->length; i++) {
    if (phi->data[i] != 0.0) break;
    if (i == phi->length - 1){
      LALSDestroyVector(status->statusPtr, &ff);
      CHECKSTATUSPTR(status);
      LALSDestroyVector(status->statusPtr, &a);
      CHECKSTATUSPTR(status);
      LALDDestroyVector(status->statusPtr, &phi);
      CHECKSTATUSPTR(status);
      if( params->approximant == IMRPhenomA){
	LALSDestroyVector(status->statusPtr, &h);
	CHECKSTATUSPTR(status);
      }

      DETATCHSTATUSPTR( status );
      RETURN( status );
    }
  }

  /* print some messages */
  sprintf(message, "fFinal = %f", params->fFinal);
  LALInfo(status, message);

  s = 0.5 * phi->data[count - 1];
  sprintf(message, "cycles = %f", s/3.14159);
  LALInfo(status, message);

  /* CHECK THIS CHECK THIS $$$$$$$$$$$$$$$$$$$$$$$$$$$*/
  sprintf(message, "final coalescence phase with respet to actual data =%f ",
          (ff->data[count-1]-ff->data[count-2])/2/3.14159);
  LALInfo(status, message);

  if ( (s/LAL_PI) < 2 ){
      sprintf(message, "The waveform has only %f cycles; we don't keep waveform with less than 2 cycles.",
	      (double) s/ (double)LAL_PI );
      LALWarning(status, message);
  }
  else {

      phiC =  phi->data[count-1] ;

      for (i=0; i<count;i++) {
	    phi->data[i] =  -phiC + phi->data[i] + ppnParams->phi;
      }

      /* Allocate the waveform structures. */
      if ( ( waveform->a = (REAL4TimeVectorSeries *)
                  LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
          ABORT( status, LALINSPIRALH_EMEM,
                  LALINSPIRALH_MSGEMEM );
      }
      memset( waveform->a, 0, sizeof(REAL4TimeVectorSeries) );
      if ( ( waveform->f = (REAL4TimeSeries *)
                  LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
          LALFree( waveform->a ); waveform->a = NULL;
          ABORT( status, LALINSPIRALH_EMEM,
                  LALINSPIRALH_MSGEMEM );
      }
      memset( waveform->f, 0, sizeof(REAL4TimeSeries) );
      if ( ( waveform->phi = (REAL8TimeSeries *)
                  LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
          LALFree( waveform->a ); waveform->a = NULL;
          LALFree( waveform->f ); waveform->f = NULL;
          ABORT( status, LALINSPIRALH_EMEM,
                  LALINSPIRALH_MSGEMEM );
      }
      memset( waveform->phi, 0, sizeof(REAL8TimeSeries) );

      in.length = (UINT4)count;
      in.vectorLength = 2;
      LALSCreateVectorSequence( status->statusPtr,
				&( waveform->a->data ), &in );
      CHECKSTATUSPTR(status);
      LALSCreateVector( status->statusPtr,
			&( waveform->f->data ), count);
      CHECKSTATUSPTR(status);
      LALDCreateVector( status->statusPtr,
			&( waveform->phi->data ), count );
      CHECKSTATUSPTR(status);

      /* copy the frequency, amplitude and phase data to the waveform structure */
      memcpy(waveform->f->data->data , ff->data, count*(sizeof(REAL4)));
      memcpy(waveform->a->data->data , a->data, 2*count*(sizeof(REAL4)));
      memcpy(waveform->phi->data->data ,phi->data, count*(sizeof(REAL8)));

      /* also set other parameters in the waveform structure */
      waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT = 1./params->tSampling;
      waveform->a->sampleUnits = lalStrainUnit;
      waveform->f->sampleUnits = lalHertzUnit;
      waveform->phi->sampleUnits = lalDimensionlessUnit;
      waveform->position = ppnParams->position;
      waveform->psi = ppnParams->psi;

      /* assign names */
      LALSnprintf( waveform->a->name,
	  	LALNameLength, "Phenom inspiral amplitudes");
      LALSnprintf( waveform->f->name,
		  LALNameLength, "Phenom inspiral frequency");
      LALSnprintf( waveform->phi->name,
	  	LALNameLength, "Phenom inspiral phase");

      /* fill some output */
      ppnParams->tc     = (double)(count-1) / params->tSampling ;
      ppnParams->length = count;
      ppnParams->dfdt   = ((REAL4)(waveform->f->data->data[count-1]
				   - waveform->f->data->data[count-2]))* ppnParams->deltaT;
      ppnParams->fStop  = params->fFinal;
      ppnParams->termCode        = GENERATEPPNINSPIRALH_EFSTOP;
      ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;
      ppnParams->fStart   = ppnParams->fStartIn;

      if( params->approximant == IMRPhenomA ){
	if ( ( waveform->h = (REAL4TimeVectorSeries *)
	       LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ){
	  ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
	}
	memset( waveform->h, 0, sizeof(REAL4TimeVectorSeries) );
	LALSCreateVectorSequence( status->statusPtr,&( waveform->h->data ), &in );
	CHECKSTATUSPTR(status);
	memcpy(waveform->h->data->data , h->data, 2*count*(sizeof(REAL4)));
	waveform->h->deltaT = 1./params->tSampling;
	waveform->h->sampleUnits = lalStrainUnit;
	LALSnprintf( waveform->h->name,
		     LALNameLength, "Phenom inspiral polarizations");
	LALSDestroyVector(status->statusPtr, &h);
	CHECKSTATUSPTR(status);
      }
  } /* end phase condition*/

  /* free memory */
  LALSDestroyVector(status->statusPtr, &ff);
  CHECKSTATUSPTR(status);
  LALSDestroyVector(status->statusPtr, &a);
  CHECKSTATUSPTR(status);
  LALDDestroyVector(status->statusPtr, &phi);
  CHECKSTATUSPTR(status);

  LALSDestroyVector(status->statusPtr, &hc);
  CHECKSTATUSPTR(status);

  LALSDestroyVector(status->statusPtr, &hp);
  CHECKSTATUSPTR(status);

/*   LALSDestroyVector(status->statusPtr, &h); */
/*   CHECKSTATUSPTR(status); */

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



static void XLALComputeInstantFreq( REAL4Vector *Freq,
			     REAL4Vector *hp,
			     REAL4Vector *hc,
			     REAL8 dt)
{
    REAL4Vector *hpDot = NULL, *hcDot = NULL;
    UINT4 k, len;

    len = hp->length;

    hpDot= XLALCreateREAL4Vector(len);
    hcDot= XLALCreateREAL4Vector(len);

    /* Construct the dot vectors (2nd order differencing) */
    hpDot->data[0] = 0.0;
    hpDot->data[len-1] = 0.0;
    hcDot->data[0] = 0.0;
    hcDot->data[len-1] = 0.0;
    for( k = 1; k < len-1; k++) {
        hpDot->data[k] = 1./(2.*dt)*(hp->data[k+1]-hp->data[k-1]);
        hcDot->data[k] = 1./(2.*dt)*(hc->data[k+1]-hc->data[k-1]);
    }

    /* Compute frequency using the fact that  */
    /*h(t) = A(t) e^(i Phi) = Re(h) + i Im(h) */
    for( k = 0; k < len; k++) {
        Freq->data[k] = -hcDot->data[k] * hp->data[k] +hpDot->data[k] * hc->data[k];
        Freq->data[k] /= LAL_TWOPI;
        Freq->data[k] /= (pow(hp->data[k],2.) + pow(hc->data[k], 2.));
    }

    /* free the memory allocated for the derivative vectors */
    XLALDestroyREAL4Vector(hpDot);
    XLALDestroyREAL4Vector(hcDot);

    return;

}

static REAL4Vector *XLALCutAtFreq( REAL4Vector     *h,
				   REAL4Vector     *freq,
				   REAL8           cutFreq,
				   REAL8           deltaT)
{
  REAL8 dt;
  UINT4 k, k0, kMid, len;
  REAL4 currentFreq;

  len = freq->length;
  dt = deltaT;

  /* Since the boundaries of this freq vector are likely to have   */
  /* FFT crap, let's scan the freq values starting from the middle */
  kMid = len/2;
  currentFreq = freq->data[kMid];
  k = kMid;

  /* freq is an increasing function of time */
  /* If we are above the cutFreq we move to the left; else to the right */
  if (currentFreq > cutFreq && k > 0)
    {
      while(currentFreq > cutFreq)
	{
	  currentFreq = freq->data[k];
	  k--;
	}
      k0 = k;
    }
  else
    {
      while(currentFreq < cutFreq && k < len)
	{
	  currentFreq = freq->data[k];
	  k++;
	}
      k0 = k;
    }

  /* Allocate memory for the frequency series */

  for(k = 0; k < k0; k++)
    {
      h->data[k] = 0.0;
    }

  return h;

}
