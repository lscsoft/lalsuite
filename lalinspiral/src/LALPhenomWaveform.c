/*
*  Copyright (C) 2008-2009 P. Ajith, Badri Krishnan, Lucia Santamaria, Nickolas Fotopoulos
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
#include <lal/TimeSeries.h>
#include <lal/SeqFactories.h>
#include <lal/VectorOps.h>
#include <lal/BBHPhenomCoeffs.h>
#include <lal/AVFactories.h>

#include <math.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

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
  REAL8 psi8;
}
BBHPhenomParams;


static void XLALComputePhenomParams(
    BBHPhenomParams *phenParams,
    InspiralTemplate *params);

static void XLALComputePhenomParams2(
    BBHPhenomParams *phenParams,
    InspiralTemplate *params);

static REAL8 XLALLorentzianFn (
    REAL8 freq,
    REAL8 fRing,
    REAL8 sigma);

static void XLALBBHPhenWaveFD (
    BBHPhenomParams *params,
    InspiralTemplate *insp_template,
    REAL4Vector *signalvec);

static void XLALBBHPhenWaveFD2 (
    BBHPhenomParams *params,
    InspiralTemplate *insp_template,
    REAL4Vector *signalvec);

static void XLALComputeInstantFreq(
    REAL4Vector *Freq,
    REAL4Vector *hp,
    REAL4Vector *hc,
    REAL8 dt);

static UNUSED REAL4Vector *XLALCutAtFreq(
    REAL4Vector *h,
    REAL4Vector *freq,
    REAL8 cutFreq,
    REAL8 deltaT);

/*
 *
 * Wrapper functions to be called from LALInspiralWave
 *
 */


/* A = non-spinning binaries. Ref. http://arxiv.org/pdf/0710.2335 */
int XLALBBHPhenWaveAFreqDom (
    REAL4Vector      *signalvec,  /**< output array */
    InspiralTemplate *params      /**< inspiral parameters */
    ) {
  BBHPhenomParams phenParams;

  /* check inputs */
  if (!signalvec || !(signalvec->data) || !params) {
    XLALPrintError(LALINSPIRALH_MSGENULL);
    XLAL_ERROR(XLAL_EFAULT);
  }
  if ((signalvec->length < 2)
      || (params->mass1 <= 0) || (params->mass2 <= 0)
      || (params->spin1[0] != 0) || (params->spin1[1] != 0) || (params->spin1[2] != 0)
      || (params->spin2[0] != 0) || (params->spin2[1] != 0) || (params->spin2[2] != 0)) {
    XLALPrintError(LALINSPIRALH_MSGECHOICE);
    XLAL_ERROR(XLAL_EDOM);
  }

  /* compute the phenomenological parameters */
  XLALComputePhenomParams(&phenParams, params);

  /* generate the phenomenological waveform in frequency domain */
  XLALBBHPhenWaveFD (&phenParams, params, signalvec);
  return XLAL_SUCCESS;
}

/* B = aligned-spin binaries. Ref. http://arxiv.org/pdf/0909.2867 */
int XLALBBHPhenWaveBFreqDom (
    REAL4Vector      *signalvec,  /**< output array */
    InspiralTemplate *params      /**< inspiral parameters */
    ) {
  BBHPhenomParams phenParams;

  /* check inputs */
  if (!signalvec || !(signalvec->data) || !params) {
    XLALPrintError(LALINSPIRALH_MSGENULL);
    XLAL_ERROR(XLAL_EFAULT);
  }
  if ((signalvec->length < 2)
      || (params->mass1 <= 0) || (params->mass2 <= 0)
      || (params->spin1[0] != 0) || (params->spin1[1] != 0)
      || (params->spin2[0] != 0) || (params->spin2[1] != 0)) {
    XLALPrintError(LALINSPIRALH_MSGECHOICE);
    XLAL_ERROR(XLAL_EDOM);
  }

  /* compute the phenomenological parameters */
  XLALComputePhenomParams2(&phenParams, params);

  /* generate the phenomenological waveform in frequency domain */
  XLALBBHPhenWaveFD2 (&phenParams, params, signalvec);

  return XLAL_SUCCESS;
}

/* A = non-spinning binaries. Ref. http://arxiv.org/pdf/0710.2335 */
int XLALBBHPhenWaveAFreqDomTemplates(
                                     REAL4Vector      *signalvec1,
                                     REAL4Vector      *signalvec2,
                                     InspiralTemplate *params) {
  if (!signalvec1 || !signalvec2 ||
      !(signalvec1->data) || !(signalvec2->data)) {
    XLALPrintError(LALINSPIRALH_MSGENULL);
    XLAL_ERROR(LALINSPIRALH_ENULL);
  }

  /* Initially the waveforms are empty */
  memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));
  memset(signalvec2->data, 0, signalvec2->length * sizeof(REAL4));

  /* generate one waveform with startPhase specified by the user */
  if (XLALBBHPhenWaveAFreqDom(signalvec1, params))
    XLAL_ERROR(XLAL_EFUNC);

  /* generate another waveform orthogonal to it */
  params->startPhase += LAL_PI_2;
  if (XLALBBHPhenWaveAFreqDom(signalvec2, params))
    XLAL_ERROR(XLAL_EFUNC);

  return 0;
}

/* B = aligned-spin binaries. Ref. http://arxiv.org/pdf/0909.2867 */
int XLALBBHPhenWaveBFreqDomTemplates(
                                     REAL4Vector      *signalvec1,
                                     REAL4Vector      *signalvec2,
                                     InspiralTemplate *params) {
  if (!signalvec1 || !signalvec2
      || !(signalvec1->data) || !(signalvec2->data)
      || !params) {
    XLALPrintError(LALINSPIRALH_MSGENULL);
    XLAL_ERROR(LALINSPIRALH_ENULL);
  }

  /* Initially the waveforms are empty */
  memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));
  memset(signalvec2->data, 0, signalvec2->length * sizeof(REAL4));

  /* generate one waveform with startPhase specified by the user */
  if (XLALBBHPhenWaveBFreqDom(signalvec1, params))
    XLAL_ERROR(XLAL_EFUNC);

  /* generate another waveform orthogonal to it */
  params->startPhase += LAL_PI_2;
  if (XLALBBHPhenWaveBFreqDom(signalvec2, params))
    XLAL_ERROR(XLAL_EFUNC);

  return 0;
}

int XLALBBHPhenWaveTimeDom(
    REAL4Vector      *signalvec1,
    InspiralTemplate *insp_template) {
  UINT4 count;

  /* check inputs */
  if (!signalvec1 || !(signalvec1->data)) {
    XLALPrintError(LALINSPIRALH_MSGENULL);
    XLAL_ERROR(LALINSPIRALH_ENULL);
  }
  if (signalvec1->length <= 2) {
    XLALPrintError(LALINSPIRALH_MSGECHOICE);
    XLAL_ERROR(XLAL_EDOM);
  }
  /* Initially the waveforms are empty */
  memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));

  if (XLALBBHPhenTimeDomEngine(signalvec1, NULL, NULL, NULL, NULL, NULL, &count, insp_template))
    XLAL_ERROR(XLAL_EFUNC);

  return 0;
}

int XLALBBHPhenWaveTimeDomTemplates(
    REAL4Vector      *signalvec1,
    REAL4Vector      *signalvec2,
    InspiralTemplate *insp_template) {
  UINT4 count;

  /* check inputs */
  if (!signalvec1 || !(signalvec1->data)
      || !signalvec2 || !(signalvec2->data)
      || !insp_template) {
    XLALPrintError(LALINSPIRALH_MSGENULL);
    XLAL_ERROR(XLAL_EFAULT);
  }
  if ((signalvec1->length <= 2) || (signalvec2->length <= 2)) {
    XLALPrintError(LALINSPIRALH_MSGECHOICE);
    XLAL_ERROR(XLAL_EDOM);
  }

  /* Initially the waveforms are empty */
  memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));
  memset(signalvec2->data, 0, signalvec2->length * sizeof(REAL4));

  return XLALBBHPhenTimeDomEngine(signalvec1, signalvec2, NULL, NULL, NULL, NULL, &count, insp_template);
}

int XLALBBHPhenTimeDomEngine(
    REAL4Vector      *signalvec1,  /**< optional output waveform with phi_c = 0 */
    REAL4Vector      *signalvec2,  /**< optional output waveform with phi_c = pi/2 */
    REAL4Vector      *h,  /**< optional output waveforms, alternating h+ and hx components */
    REAL4Vector      *aVec,   /**< optional output inst. amplitude, alternating A+ and Ax components, assuming that h+ & hx have equal ...*/
    REAL4Vector      *freqVec,  /**< optional output instant. freq */
    REAL8Vector      *phiVec,  /**< optional output phase evolution */
    UINT4            *countback,  /**< output number of non-zero samples */
    InspiralTemplate *params) /**< UNDOCUMENTED */ {
  REAL8 dt, cosI, fLower, peakAmp, fCut, fRes, f, totalMass, softWin, z1, z2;
  REAL8 fLowerOrig, eta, tau0, winFLo, sigLo, sigHi, tF0, expectedAmplRatio;
  REAL8 phaseShift, sig1, sig2, startPhaseOrig, phiC;
  REAL4 windowLength;
  UINT4 i, j, k, l, n, peakAmpIdx, sigLength, iLower;
  REAL4Vector *signalFD1 = NULL, *signalFD2 = NULL, *signalTD1 = NULL, *signalTD2 = NULL, *a=NULL, *fVec=NULL;
  REAL8Vector *phi=NULL;
  REAL4FFTPlan *revPlan = NULL;
  BBHPhenomParams phenParams;

  /* check inputs */
  if (!params) {
    XLALPrintError(LALINSPIRALH_MSGENULL);
    XLAL_ERROR(XLAL_EFAULT);
  }
  if ((params->nStartPad < 0) || (params->nEndPad < 0)
      || (params->fLower <= 0) || (params->tSampling <= 0)
      || (params->mass1 < 0) || (params->mass2 < 0)
      || (params->spin1[0] != 0) || (params->spin1[1] != 0)
      || (params->spin2[0] != 0) || (params->spin2[1] != 0)) {
    XLALPrintError(LALINSPIRALH_MSGECHOICE);
    XLAL_ERROR(XLAL_EDOM);
  }

  dt = 1./params->tSampling;
  cosI = cos( params->inclination );

  /* compute the phenomenological parameters */
  switch (params->approximant) {
      /* non-spinning binaries. Ref. http://arxiv.org/pdf/0710.2335 */
    case IMRPhenomA:
      if ((params->spin1[2] != 0) || (params->spin2[2] != 0)) {
        XLALPrintError(LALINSPIRALH_MSGECHOICE);
        XLAL_ERROR(XLAL_EDOM);
      }
      XLALComputePhenomParams(&phenParams, params);
      break;
      /* aligned-spin binaries. Ref. http://arxiv.org/abs/0909.2867 */
    case IMRPhenomB:
      XLALComputePhenomParams2(&phenParams, params);
      break;
    default:
      XLALPrintError(LALINSPIRALH_MSGESWITCH);
      XLAL_ERROR(XLAL_ENAME);
  }

  totalMass = params->mass1 + params->mass2;
  eta = params->mass1*params->mass2/pow(totalMass,2.);

  /* we will generate the waveform from a frequency which is lower than the
   * fLower chosen. Also the cutoff frequency is higher than the fCut. We
   * will later apply a window function, and truncate the time-domain waveform
   * below an instantaneous frequency  fLower */
  fLowerOrig = params->fLower;    /* this is the low-freq set by the user */
  startPhaseOrig = params->startPhase; /* this is the coalescence phase set by the user*/

  /* Find an optimum value for fLower (using the definition of Newtonian chirp time)
   * such that the waveform has a minimum length of tau0. This is necessary to avoid
   * FFT artifacts */
  tau0 = 64.;
  fLower = pow((tau0*256.*eta*pow(totalMass*LAL_MTSUN_SI,5./3.)/5.),-3./8.)/LAL_PI;
  fCut = (1.025)*phenParams.fCut;

  /* make sure that these frequencies are not too out of range */
  if (fLower > fLowerOrig) fLower = fLowerOrig;
  if (fLower < 0.5) fLower = 0.5;
  if (fCut > params->tSampling/2.-100.) fCut = params->tSampling/2.-100.;

  /* generate waveforms over this frequency range */
  params->fLower = fLower;
  phenParams.fCut = params->tSampling/2.;

  /* allocate memory for the freq and domain tempaltes - templates are generated 
   * using a lower value of fLower than prescribed by the user */
  tau0 = 5.*totalMass*LAL_MTSUN_SI/(256.*eta*pow(LAL_PI*totalMass*LAL_MTSUN_SI*fLower, 8./3.));
  n = pow(2., ceil(log2(tau0*params->tSampling)));
  signalFD1 = XLALCreateREAL4Vector(n);
  signalFD2 = XLALCreateREAL4Vector(n);
  signalTD1 = XLALCreateREAL4Vector(n);
  signalTD2 = XLALCreateREAL4Vector(n);
  if (!signalFD1 || !signalFD2 || !signalTD1 || !signalTD2) {
    if (signalFD1) XLALDestroyREAL4Vector(signalFD1);
    if (signalFD2) XLALDestroyREAL4Vector(signalFD2);
    if (signalTD1) XLALDestroyREAL4Vector(signalTD1);
    if (signalTD2) XLALDestroyREAL4Vector(signalTD2);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* generate the phenomenological waveform in frequency domain */
  switch (params->approximant) {
      /* non-spinning binaries. Ref. http://arxiv.org/pdf/0710.2335 */
    case IMRPhenomA:
      XLALBBHPhenWaveFD (&phenParams, params, signalFD1);
      params->startPhase += LAL_PI_2; 
      XLALBBHPhenWaveFD (&phenParams, params, signalFD2);
      break;
      /* aligned-spin binaries. Ref. http://arxiv.org/abs/0909.2867 */
    case IMRPhenomB:
      XLALBBHPhenWaveFD2 (&phenParams, params, signalFD1);
      params->startPhase += LAL_PI_2; 
      XLALBBHPhenWaveFD2 (&phenParams, params, signalFD2);
      break;
    default:
      XLALDestroyREAL4Vector(signalFD1);
      XLALDestroyREAL4Vector(signalFD2);
      XLALDestroyREAL4Vector(signalTD1);
      XLALDestroyREAL4Vector(signalTD2);
      XLALPrintError(LALINSPIRALH_MSGESWITCH);
      XLAL_ERROR(XLAL_ENAME);
  }

  /* apply the softening window function */
  fRes = params->tSampling/n;

  winFLo = (fLowerOrig + fLower)/2.;
  // UNUSED!!: REAL8 winFHi = (fCut + phenParams.fCut)/2.;
  sigLo = 4.;
  sigHi = 4.;

  signalFD1->data[0] = 0.;
  signalFD2->data[0] = 0.;
  for (k = 1; k <= n/2; k++) {
    f = k*fRes;
    softWin = (1+tanh((4.*(f-winFLo)/sigLo)))*(1-tanh(4.*(f-fCut)/sigHi))/4.;
    signalFD1->data[k] *= softWin;
    signalFD1->data[n-k] *= softWin;
    signalFD2->data[k] *= softWin;
    signalFD2->data[n-k] *= softWin;
  }

  /* Inverse Fourier transform */
  revPlan = XLALCreateReverseREAL4FFTPlan(n, 0);
  if (!revPlan) {
    XLALDestroyREAL4Vector(signalFD1);
    XLALDestroyREAL4Vector(signalFD2);
    XLALDestroyREAL4Vector(signalTD1);
    XLALDestroyREAL4Vector(signalTD2);
    XLALPrintError(LALINSPIRALH_MSGEMEM);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  XLALREAL4VectorFFT(signalTD1, signalFD1, revPlan);
  XLALREAL4VectorFFT(signalTD2, signalFD2, revPlan);
  XLALDestroyREAL4Vector(signalFD1);
  XLALDestroyREAL4Vector(signalFD2);
  XLALDestroyREAL4FFTPlan(revPlan);

  /* FFT normalisation. The LAL implementation of the FFT omits the factor 1/n.*/
  for (i = 0; i < n; i++) {
    signalTD1->data[i] *= params->tSampling/n;
    signalTD2->data[i] *= params->tSampling/n;
  }

  /* apply a linearly decresing window at the end
   * of the waveform in order to avoid edge effects. */
  windowLength = 10.*totalMass * LAL_MTSUN_SI*params->tSampling;
  for (i=0; i< windowLength; i++){
    signalTD1->data[n-i-1] *= i/windowLength;
    signalTD2->data[n-i-1] *= i/windowLength;
  }

  /* compute the instantaneous frequency */
  a = XLALCreateREAL4Vector(n);
  fVec = XLALCreateREAL4Vector(n);
  phi = XLALCreateREAL8Vector(n);

  XLALComputeInstantFreq(fVec, signalTD1, signalTD2, dt);
  peakAmp = 0.;
  peakAmpIdx = 0;

  /* compute the amplitude and phase. find the peak amplitude */
  for (i=0; i<n; i++){
    a->data[i] = sqrt(pow(signalTD1->data[i],2.) + pow(signalTD2->data[i],2.));
    phi->data[i] = -atan2(signalTD2->data[i], signalTD1->data[i]);

    /* find the peak amplitude*/
    if (a->data[i] > peakAmp) {
      peakAmp = a->data[i];
      peakAmpIdx = i;
    }
  }

  /* If the instantaneous amplitude a(t) is too low, noise will corrupt the
   * estimation of the instantaneous frequency f(t). Hence, if the ratio of
   * a(t) with the peak amplitude is lower than a threshold, we set f(t) = 0
   * for those time bins. Choosing of such a threshold is a tricky business!
   * The current value is chosen using the following argument: a(t) is
   * proportional to the square of the speed v(t) of the binary (the PN
   * parameter). At the merger, v(t) ~ 1. At the start of the waveform
   * v(t) = (pi M f_lower)^{1/3}. Thus, the expected ratio at the start of
   * the waveform is a_0/a_peak = (pi M f_lower)^{-2/3}. Divide by 2 as a
   * safety margin */
  expectedAmplRatio = pow(LAL_PI*LAL_MTSUN_SI*totalMass*fLowerOrig, 2./3.)/2.;

  /* Also find the index of the fVec corresponding to fLowerOrig */
  iLower = 0;
  for (i=0; i< fVec->length; i++) {
    if (a->data[i] < expectedAmplRatio*peakAmp) {
      fVec->data[i] = 0.0;
    }
    if ((iLower == 0) && (fVec->data[i] >= fLowerOrig)) iLower = i;
  }

  /* unwrap the phase */
  XLALREAL8VectorUnwrapAngle(phi, phi);

  /* reassign the original values of fLower, fCut and startPhase */
  params->fLower = fLowerOrig;
  params->fFinal = fCut;
  params->startPhase = startPhaseOrig;

  /* apply a phase-shift to the waveforms such that the phase at coalescence is
   * equal to params->startPhase */
  phaseShift = phi->data[peakAmpIdx] + params->startPhase;
  phiC = phi->data[peakAmpIdx] - params->startPhase;
  for (i=iLower; i < fVec->length; i++) {
    sig1 = signalTD1->data[i]*cos(phaseShift) - signalTD2->data[i]*sin(phaseShift);
    sig2 = signalTD1->data[i]*sin(phaseShift) + signalTD2->data[i]*cos(phaseShift);
    signalTD1->data[i] = sig1;
    signalTD2->data[i] = sig2;
    phi->data[i] -= phiC;
  }

  /* copy the frequency components above fLower to the signal vector to be returned */
  j = params->nStartPad;
  tF0 = iLower*dt;           /* time bin corresponding to f(t) = fLower */

  sigLength = 0;
  if (signalvec1) sigLength = signalvec1->length;
  else if (signalvec2) sigLength = signalvec2->length;
  else if (freqVec) sigLength = freqVec->length;
  else if (phiVec) sigLength = phiVec->length;
  else if (aVec) sigLength = aVec->length / 2;
  else if (h) sigLength = h->length / 2;

  /* inclination-weights on two polarizations */
  z1 = -0.5*(1. + cosI*cosI);
  z2 = -cosI;

  for (i=iLower; i < fVec->length; i++) {
    if (j<sigLength) {
      k = 2*j;
      l = k+1;

      if (signalvec1) signalvec1->data[j] = signalTD1->data[i]; /* signal with phi_c = 0 */
      if (signalvec2) signalvec2->data[j] = signalTD2->data[i]; /* signal with phic_c = pi/2 */
      if (freqVec) freqVec->data[j] = fVec->data[i]; /* instant. freq */
      if (phiVec) phiVec->data[j] = phi->data[i];  /* phase evolution */
      if (aVec) {
        aVec->data[k] = a->data[i]; /* inst. amplitude, assuming that h+ & hx have equal ...*/ 
        aVec->data[l] = a->data[i]; /* ... amplitude. valid in the absence of precession */
      }
      if (h) {
        h->data[k] = z1 * signalTD1->data[i]; /* h+ ~ -[ 1+cos^2(iota) ]/2 cos (phi) */
        h->data[l] = z2 * signalTD2->data[i]; /* hx ~ -cos(iota) sin (phi) */
      }
      j++;
    }
  }
  *countback = j;

  /* free the memory */
  XLALDestroyREAL4Vector(a);
  XLALDestroyREAL4Vector(fVec);
  XLALDestroyREAL8Vector(phi);
  XLALDestroyREAL4Vector(signalTD1);
  XLALDestroyREAL4Vector(signalTD2);

  /* store some parameters for record keeping */
  params->vFinal = pow(LAL_PI*LAL_MTSUN_SI*totalMass*params->fFinal, 1./3.);
  params->tC = peakAmpIdx*dt-tF0;   /* time of coalescence. defined as the time
                                     corresponding to the peak amplitude*/

  return 0;
}

int XLALBBHPhenWaveTimeDomForInjection (
                                        CoherentGW       *waveform,  /**< allocated, but completely zeroed CoherentGW structure; this function allocates sub-structures that you must free */
                                        InspiralTemplate *params,	/**< UNDOCUMENTED */
                                        PPNParamStruc    *ppnParams) /**< UNDOCUMENTED */ {
  REAL4Vector *a=NULL;      /* amplitude  data */
  REAL4Vector *h=NULL;      /* polarization data */
  REAL4Vector *ff=NULL;     /* frequency data */
  REAL8Vector *phi=NULL;    /* phase data */

  UINT4 count, i;
  REAL8 s, phiC;            /* phase at coalescence */
  CHAR message[256];
  LIGOTimeGPS zero_time = {0, 0};
  InspiralInit paramsInit;
  int returnval = XLAL_SUCCESS;

  /* check inputs */
  if (!params || !waveform) {
    XLALPrintError(LALINSPIRALH_MSGENULL);
    XLAL_ERROR(XLAL_EFAULT);
  }
  if ((params->nStartPad < 0) || (params->nEndPad < 0)
      || (params->fLower <= 0) || (params->tSampling <= 0)) {
    XLALPrintError(LALINSPIRALH_MSGECHOICE);
    XLAL_ERROR(XLAL_EDOM);
  }

  /* Make sure waveform fields don't exist */
  if (waveform->a || waveform->h || waveform->f
      || waveform->phi || waveform->shift) {
    XLALPrintError(LALINSPIRALH_MSGENULL);
    XLAL_ERROR(XLAL_EFAULT);
  }

  if (params->ampOrder) {
    params->ampOrder = (LALPNOrder) 0;
    snprintf(message, 256, "WARNING: Amp Order has been reset to %d", params->ampOrder);
    XLALPrintInfo(message);
  }

  /* Compute some parameters */
  if (XLALInspiralInit(params, &paramsInit))
    XLAL_ERROR(XLAL_EFUNC);
  if (paramsInit.nbins == 0) return XLAL_SUCCESS;

  /* allocate temporary structures */
  h = XLALCreateREAL4Vector(2 * paramsInit.nbins);
  a = XLALCreateREAL4Vector(2 * paramsInit.nbins);
  ff = XLALCreateREAL4Vector(paramsInit.nbins);
  phi = XLALCreateREAL8Vector(paramsInit.nbins);
  memset(h->data,  0, 2 * paramsInit.nbins * sizeof(REAL4));
  memset(a->data, 0, 2 * paramsInit.nbins * sizeof(REAL4));
  memset(ff->data, 0, paramsInit.nbins * sizeof(REAL4));
  memset(phi->data, 0, paramsInit.nbins * sizeof(REAL8));
  if (!h || !a || !ff || !phi) {
    XLALError(__func__, __FILE__, __LINE__, XLAL_EFAULT);
    returnval = XLAL_FAILURE;
    goto done;
  }

  /* generate two orthogonal waveforms */
  params->startPhase = ppnParams->phi;
  if (XLALBBHPhenTimeDomEngine(NULL, NULL, h, a, ff, phi, &count, params)) {
    XLALError(__func__, __FILE__, __LINE__, XLAL_EFUNC);
    returnval = XLAL_FAILURE;
    goto done;
  }

  /* Check an empty waveform hasn't been returned */
  for (i = 0; i < phi->length; i++) {
    if (phi->data[i] != 0.0) break;
  }
  if (i == phi->length) {
    XLALError(__func__, __FILE__, __LINE__, XLAL_ERANGE);
    returnval = XLAL_FAILURE;
    goto done;
  }

  /* print some messages */
  snprintf(message, 256, "fFinal = %f", params->fFinal);
  XLALPrintInfo(message);
  s = 0.5 * (phi->data[count-1] - phi->data[0]);
  snprintf(message, 256, "cycles = %f", s/3.14159);
  XLALPrintInfo(message);
  if (s < LAL_TWOPI){
    snprintf(message, 256, "The waveform has only %f cycles; we don't "
             "keep waveform with less than 2 cycles.", (double) s/ (double)LAL_PI );
    XLALPrintWarning(message);
    goto done;
  }

  /* compute phase */
  phiC =  phi->data[count-1] ;
  for (i=0; i<count;i++) {
    phi->data[i] =  -phiC + phi->data[i] + ppnParams->phi;
  }

  /* Allocate the CoherentGW structures. */
  waveform->h = (REAL4TimeVectorSeries *) XLALMalloc(sizeof(REAL4TimeVectorSeries));
  waveform->a = (REAL4TimeVectorSeries *) XLALMalloc(sizeof(REAL4TimeVectorSeries));
  if (!(waveform->h) || !(waveform->a)) {
    XLALError(__func__, __FILE__, __LINE__, XLAL_ENOMEM);
    returnval = XLAL_FAILURE;
    goto done;
  }
  memset(waveform->h, 0, sizeof(REAL4TimeVectorSeries));
  memset(waveform->a, 0, sizeof(REAL4TimeVectorSeries));
  waveform->h->data = XLALCreateREAL4VectorSequence(count, 2);
  waveform->a->data = XLALCreateREAL4VectorSequence(count, 2);
  waveform->f = XLALCreateREAL4TimeSeries("Phenom inspiral frequency", &zero_time, 0, 1. / params->tSampling, &lalHertzUnit, count);
  waveform->phi = XLALCreateREAL8TimeSeries("Phenom inspiral phase", &zero_time, 0, 1 / params->tSampling, &lalDimensionlessUnit, count);
  if (!(waveform->h->data) || !(waveform->a->data) || !(waveform->f) || !(waveform->phi)) {
    XLALError(__func__, __FILE__, __LINE__, XLAL_ENOMEM);
    returnval = XLAL_FAILURE;
    goto done;
  }

  /* copy the frequency, amplitude and phase data to the waveform structure */
  /*
   * Truncation (count < paramsInit.nbins in general) is necessary
   * because tC is defined at the end of the vector in some codes and not
   * others.
   */
  memcpy(waveform->h->data->data, h->data, 2*count*(sizeof(REAL4)));
  memcpy(waveform->a->data->data, a->data, 2*count*(sizeof(REAL4)));
  memcpy(waveform->f->data->data, ff->data, count*(sizeof(REAL4)));
  memcpy(waveform->phi->data->data, phi->data, count*(sizeof(REAL8)));

  /* also set other parameters in the waveform structure */
  waveform->h->sampleUnits = waveform->a->sampleUnits = lalStrainUnit;
  waveform->h->deltaT = waveform->a->deltaT = 1./params->tSampling;
  waveform->position = ppnParams->position;
  waveform->psi = ppnParams->psi;
  snprintf(waveform->h->name, LALNameLength, "Phenom inspiral polarizations");
  snprintf(waveform->a->name, LALNameLength, "Phenom inspiral amplitudes");
  /* epoch is set elsewhere */

  /* fill some output */
  ppnParams->tc     = (double)(count-1) / params->tSampling;
  ppnParams->length = count;
  ppnParams->dfdt   = ((REAL4)(waveform->f->data->data[count-1]
                               - waveform->f->data->data[count-2]))* ppnParams->deltaT;
  ppnParams->fStop  = params->fFinal;
  ppnParams->termCode        = GENERATEPPNINSPIRALH_EFSTOP;
  ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;
  ppnParams->fStart   = ppnParams->fStartIn;

done:

  if (ff) XLALDestroyREAL4Vector(ff);
  if (a) XLALDestroyREAL4Vector(a);
  if (phi) XLALDestroyREAL8Vector(phi);
  if (h) XLALDestroyREAL4Vector(h);
  if (returnval != XLAL_SUCCESS) {
    if (waveform->h && waveform->h->data) XLALDestroyREAL4VectorSequence(waveform->h->data);
    if (waveform->h) XLALFree(waveform->h);
    if (waveform->a && waveform->a->data) XLALDestroyREAL4VectorSequence(waveform->a->data);
    if (waveform->a) XLALFree(waveform->a);
    if (waveform->f) XLALDestroyREAL4TimeSeries(waveform->f);
    if (waveform->phi) XLALDestroyREAL8TimeSeries(waveform->phi);
  }
  return returnval;
}


/*
 *
 * Core numerical routines; the Phenom prescription is in FD
 *
 */

/*********************************************************************/
/* Compute phenomenological parameters for non-spinning binaries     */
/* Ref. Eq.(4.18) of http://arxiv.org/pdf/0710.2335 and              */
/* Table I of http://arxiv.org/pdf/0712.0343                         */
/*********************************************************************/
static void XLALComputePhenomParams(
    BBHPhenomParams  *phenParams,
    InspiralTemplate *params) {
  REAL8 totalMass, piM, eta, fMerg_a, fMerg_b, fMerg_c, fRing_a, fRing_b, etap2;
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
    etap2 = eta*eta;
    phenParams->fCut  = (fCut_a*etap2  + fCut_b*eta  + fCut_c)/piM;
    phenParams->fMerger  = (fMerg_a*etap2  + fMerg_b*eta  + fMerg_c)/piM;
    phenParams->fRing  = (fRing_a*etap2 + fRing_b*eta + fRing_c)/piM;
    phenParams->sigma = (sigma_a*etap2 + sigma_b*eta + sigma_c)/piM;

    phenParams->psi0 = (psi0_a*etap2 + psi0_b*eta + psi0_c)/(eta*pow(piM, 5./3.));
    phenParams->psi1 = 0.;
    phenParams->psi2 = (psi2_a*etap2 + psi2_b*eta + psi2_c)/(eta*pow(piM, 3./3.));
    phenParams->psi3 = (psi3_a*etap2 + psi3_b*eta + psi3_c)/(eta*pow(piM, 2./3.));
    phenParams->psi4 = (psi4_a*etap2 + psi4_b*eta + psi4_c)/(eta*pow(piM, 1./3.));
    phenParams->psi5 = 0.;
    phenParams->psi6 = (psi6_a*etap2 + psi6_b*eta + psi6_c)/(eta*pow(piM, -1./3.));
    phenParams->psi7 = (psi7_a*etap2 + psi7_b*eta + psi7_c)/(eta*pow(piM, -2./3.));
  }

  return;

}

/*********************************************************************/
/* Compute phenomenological parameters for aligned-spin binaries     */
/* Ref. Eq.(2) and Table I of http://arxiv.org/pdf/0909.2867         */
/*********************************************************************/
static void XLALComputePhenomParams2(
  BBHPhenomParams  *phenParams,
  InspiralTemplate *params) {

  REAL8 totalMass, piM, eta, chi, delta;
  REAL8 etap2, chip2, etap3, etap2chi, etachip2, etachi; 

  if (!params || !phenParams) return;

  /* calculate the total mass and symmetric mass ratio */
  totalMass = params->mass1+params->mass2;
  eta = params->mass1*params->mass2/pow(totalMass,2.);
  piM = totalMass*LAL_PI*LAL_MTSUN_SI;
  delta = sqrt(1.-4.*eta); /* asymmetry parameter */

  /* spin parameter used for the search */
  chi = 0.5*(params->spin1[2]*(1.+delta) + params->spin2[2]*(1.-delta));

  /* spinning phenomenological waveforms */
  etap2 = eta*eta;
  chip2 = chi*chi;
  etap3 = etap2*eta;
  etap2chi = etap2*chi;
  etachip2 = eta*chip2;	
  etachi = eta*chi;

  phenParams->psi0 = 3./(128.*eta);

  phenParams->psi2 = 3715./756. +
  -9.2091e+02*eta + 4.9213e+02*etachi + 1.3503e+02*etachip2 +
  6.7419e+03*etap2 + -1.0534e+03*etap2chi +
  -1.3397e+04*etap3 ;

  phenParams->psi3 = -16.*LAL_PI + 113.*chi/3. +
  1.7022e+04*eta + -9.5659e+03*etachi + -2.1821e+03*etachip2 +
  -1.2137e+05*etap2 + 2.0752e+04*etap2chi +
  2.3859e+05*etap3 ;

  phenParams->psi4 = 15293365./508032. - 405.*chip2/8. +
  -1.2544e+05*eta + 7.5066e+04*etachi + 1.3382e+04*etachip2 +
  8.7354e+05*etap2 + -1.6573e+05*etap2chi +
  -1.6936e+06*etap3 ;

  phenParams->psi6 = -8.8977e+05*eta + 6.3102e+05*etachi + 5.0676e+04*etachip2 +
  5.9808e+06*etap2 + -1.4148e+06*etap2chi +
  -1.1280e+07*etap3 ;

  phenParams->psi7 = 8.6960e+05*eta + -6.7098e+05*etachi + -3.0082e+04*etachip2 +
  -5.8379e+06*etap2 + 1.5145e+06*etap2chi +
  1.0891e+07*etap3 ;

  phenParams->psi8 = -3.6600e+05*eta + 3.0670e+05*etachi + 6.3176e+02*etachip2 +
  2.4265e+06*etap2 + -7.2180e+05*etap2chi + 
  -4.5524e+06*etap3;

  phenParams->fMerger =  1. - 4.4547*pow(1.-chi,0.217) + 3.521*pow(1.-chi,0.26) +
  6.4365e-01*eta + 8.2696e-01*etachi + -2.7063e-01*etachip2 +
  -5.8218e-02*etap2 + -3.9346e+00*etap2chi +
  -7.0916e+00*etap3 ;

  phenParams->fRing = (1. - 0.63*pow(1.-chi,0.3))/2. +
  1.4690e-01*eta + -1.2281e-01*etachi + -2.6091e-02*etachip2 +
  -2.4900e-02*etap2 + 1.7013e-01*etap2chi +
  2.3252e+00*etap3 ;

  phenParams->sigma = (1. - 0.63*pow(1.-chi,0.3))*pow(1.-chi,0.45)/4. +
  -4.0979e-01*eta + -3.5226e-02*etachi + 1.0082e-01*etachip2 +
  1.8286e+00*etap2 + -2.0169e-02*etap2chi +
  -2.8698e+00*etap3 ;

  phenParams->fCut = 3.2361e-01 + 4.8935e-02*chi + 1.3463e-02*chip2 +
  -1.3313e-01*eta + -8.1719e-02*etachi + 1.4512e-01*etachip2 +
  -2.7140e-01*etap2 + 1.2788e-01*etap2chi +
  4.9220e+00*etap3 ;

  phenParams->fCut   /= piM;
  phenParams->fMerger/= piM;
  phenParams->fRing  /= piM;
  phenParams->sigma  /= piM;

  phenParams->psi1    = 0.;
  phenParams->psi5    = 0.;
}

/*********************************************************************/
/* Compute frequency-domain waveforms for non-spinning binaries      */
/* Ref. Eq.(4.13) and (4.16) of http://arxiv.org/pdf/0710.2335       */
/*********************************************************************/
static void XLALBBHPhenWaveFD(
    BBHPhenomParams  *params,
    InspiralTemplate *insp_template,
    REAL4Vector *signalvec) {

  REAL8 df, shft, phi, amp0, ampEff=0, psiEff, fMerg, fNorm;
  REAL8 f, fRing, sigma, totalMass, eta;
  INT4 i, j, n, nby2;

  /* freq resolution and the low-freq bin */
  df = insp_template->tSampling/signalvec->length;
  n = signalvec->length;

  shft = LAL_TWOPI*(insp_template->nStartPad/insp_template->tSampling + insp_template->startTime);
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

  nby2 = n/2; 

  /* now generate the waveform at all frequency bins */
  for (i=1; i<nby2; i++) {
    /* this is the index of the imaginary part */
    j = n-i;

    /* fourier frequency corresponding to this bin */
    f = i * df;
    fNorm = f/fMerg;

    /* compute the amplitude */
    if ((f < insp_template->fLower) || (f > params->fCut))
        ampEff = 0.;
    else if (f <= fMerg)
        ampEff = amp0*pow(fNorm, -7./6.);
    else if ((f > fMerg) & (f <= fRing))
        ampEff = amp0*pow(fNorm, -2./3.);
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
    *(signalvec->data+i) = (REAL4) (ampEff * cos(psiEff));  /* real */
    *(signalvec->data+j) = (REAL4) (ampEff * sin(psiEff));  /* imag */
  }

}


/*********************************************************************/
/* Compute frequency-domain waveforms for aligned-spin binaries      */
/* Ref. Eq.(1) of http://arxiv.org/pdf/0909.2867                     */
/*********************************************************************/
static void XLALBBHPhenWaveFD2(
    BBHPhenomParams  *params,
    InspiralTemplate *insp_template,
    REAL4Vector *signalvec) {

  REAL8 df, shft, phi, amp0, ampEff=0, psiEff, fMerg, fNorm;
  REAL8 f, fRing, sigma, totalMass, eta;
  INT4 i, j, n, nby2;
  REAL8 v, alpha2, alpha3, w1, vMerg, v2, v3, v4, v5, v6, v7, v8;
  REAL8 epsilon_1, epsilon_2, w2, vRing, chi, mergPower, delta;

  /* freq resolution and the low-freq bin */
  df = insp_template->tSampling/signalvec->length;
  n = signalvec->length;

  shft = LAL_TWOPI*(insp_template->nStartPad/insp_template->tSampling + insp_template->startTime);

  /* the negative sign is introduced such that the definition of the
   * polarisations (phi = 0 for plus, phi = pi/2 for cross) is consistent
   * with the IMRPhenomA waveforms */
  phi  = -insp_template->startPhase;

  /* phenomenological  parameters*/
  fMerg = params->fMerger;
  fRing = params->fRing;
  sigma = params->sigma;

  /* physical parameters*/
  totalMass = insp_template->mass1 + insp_template->mass2;
  eta = insp_template->eta = insp_template->mass1 * insp_template->mass2 / pow(totalMass, 2.);
  delta = sqrt(1.-4.*eta);
  chi = 0.5*(insp_template->spin1[2]*(1.+delta) + insp_template->spin2[2]*(1.-delta));

  /* Now compute the amplitude.  NOTE the params->distance is assumed to
   * me in meters. This is, in principle, inconsistent with the LAL
   * documentation (inspiral package). But this seems to be the convention
   * employed in the injection codes */
  amp0 = pow(LAL_MTSUN_SI*totalMass, 5./6.)*pow(fMerg,-7./6.)/pow(LAL_PI,2./3.);
  amp0 *= pow(5.*eta/24., 1./2.)/(insp_template->distance/LAL_C_SI);

  /* fill the zero and Nyquist frequency with zeros */
  *(signalvec->data+0) = 0.;
  *(signalvec->data+n/2) = 0.;

  /***********************************************************************/
  /* these are the parameters required for the "new" phenomenological IMR
   * waveforms*/
  /***********************************************************************/

  /* PN corrections to the frequency domain amplitude of the (2,2) mode */
  alpha2   = -323./224. + 451.*insp_template->eta/168.;
  alpha3   = (27./8. - 11.*insp_template->eta/6.)*chi;

  /* leading order power law of the merger amplitude */
  mergPower = -2./3.;

  /* spin-dependant corrections to the merger amplitude */
  epsilon_1 =  1.4547*chi - 1.8897;
  epsilon_2 = -1.8153*chi + 1.6557;

  /* normalisation constant of the inspiral amplitude */
  vMerg = pow(LAL_PI*totalMass*LAL_MTSUN_SI*fMerg, 1./3.);
  vRing = pow(LAL_PI*totalMass*LAL_MTSUN_SI*fRing, 1./3.);

  w1 = 1. + alpha2*pow(vMerg,2.) + alpha3*pow(vMerg,3.);
  w1 = w1/(1. + epsilon_1*vMerg + epsilon_2*vMerg*vMerg);
  w2 = w1*(LAL_PI*sigma/2.)*pow(fRing/fMerg, mergPower)*(1. + epsilon_1*vRing
          + epsilon_2*vRing*vRing);

  /***********************************************************************/
  /* now generate the waveform at all frequency bins */
  /***********************************************************************/
  nby2 = n/2; 
  for (i=1; i<nby2; i++) {
    /* this is the index of the imaginary part */
    j = n-i;

    /* fourier frequency corresponding to this bin */
    f = i * df;
    fNorm = f/fMerg;

    /* PN expansion parameter */
    v = pow(LAL_PI*totalMass*LAL_MTSUN_SI*f, 1./3.);

    v2 = v*v; v3 = v2*v; v4 = v2*v2; v5 = v4*v; v6 = v3*v3; v7 = v6*v, v8 = v7*v;

    /* compute the amplitude */
    if ((f < insp_template->fLower) || (f > params->fCut))
      ampEff = 0.;
    else if (f <= fMerg)
      ampEff = pow(fNorm, -7./6.)*(1. + alpha2*v2 + alpha3*v3);
    else if ((f > fMerg) & (f <= fRing))
      ampEff = w1*pow(fNorm, mergPower)*(1. + epsilon_1*v + epsilon_2*v2);
    else if (f > fRing)
      ampEff = w2*XLALLorentzianFn ( f, fRing, sigma);

    /* now compute the phase */
    psiEff =  shft*f + phi
      + 3./(128.*eta*v5)*(1 + params->psi2*v2
      + params->psi3*v3 + params->psi4*v4
      + params->psi5*v5 + params->psi6*v6
      + params->psi7*v7 + params->psi8*v8);

    /* generate the waveform */
    *(signalvec->data+i) = (REAL4) (amp0 * ampEff * cos(psiEff));  /* real */
    *(signalvec->data+j) = (REAL4) (-amp0 * ampEff * sin(psiEff));  /* imag */
  }
}


/*
 * Private support functions
 */

static void XLALComputeInstantFreq(
    REAL4Vector *Freq,
    REAL4Vector *hp,
    REAL4Vector *hc,
    REAL8 dt) {
  REAL8Vector *hpDot = NULL, *hcDot = NULL;
  UINT4 k, len;
  REAL8 fOfT;

  len = hp->length;

  hpDot= XLALCreateREAL8Vector(len);
  hcDot= XLALCreateREAL8Vector(len);

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
  /*h(t) = A(t) e^(-i Phi) = Re(h) - i Im(h) = h_+ - i h_x */
  for( k = 0; k < len; k++) {
    fOfT = -hcDot->data[k] * hp->data[k] +hpDot->data[k] * hc->data[k];
    fOfT /= LAL_TWOPI;
    fOfT /= (pow(hp->data[k],2.) + pow(hc->data[k], 2.));
    Freq->data[k] = (REAL4) fOfT;
  }

  /* free the memory allocated for the derivative vectors */
  XLALDestroyREAL8Vector(hpDot);
  XLALDestroyREAL8Vector(hcDot);

  return;

}

static REAL8 XLALLorentzianFn(
    REAL8 freq,
    REAL8 fRing,
    REAL8 sigma) {
  return sigma / (2 * LAL_PI * ((freq - fRing)*(freq - fRing)
    + sigma*sigma / 4.0));
}

static REAL4Vector *XLALCutAtFreq(
    REAL4Vector     *h,
    REAL4Vector     *freq,
    REAL8           cutFreq,
    REAL8           UNUSED deltaT) {
  UINT4 k, k0, kMid, len;
  REAL4 currentFreq;
  len = freq->length;

  /* Since the boundaries of this freq vector are likely to have   */
  /* FFT crap, let's scan the freq values starting from the middle */
  kMid = len/2;
  currentFreq = freq->data[kMid];
  k = kMid;

  /* freq is an increasing function of time */
  /* If we are above the cutFreq we move to the left; else to the right */
  if (currentFreq > cutFreq && k > 0) {
    while (currentFreq > cutFreq)
      currentFreq = freq->data[k--];
    k0 = k;
  }
  else {
    while (currentFreq < cutFreq && k < len)
      currentFreq = freq->data[k++];
    k0 = k;
  }

  /* zero the waveform below the cutoff frequency */
  for (k = 0; k < k0; k++)
      h->data[k] = 0.0;

  return h;
}


/*
 * Deprecated LAL wrapper functions
 */

void LALBBHPhenWaveFreqDom(
    LALStatus        *status,
    REAL4Vector      *signalvec,
    InspiralTemplate *params) {
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALBBHPhenWaveFreqDom");
  switch (params->approximant) {
    case IMRPhenomA:
      if (XLALBBHPhenWaveAFreqDom(signalvec, params)) ABORTXLAL(status);
      break;
    case IMRPhenomB:
      if (XLALBBHPhenWaveBFreqDom(signalvec, params)) ABORTXLAL(status);
      break;
    default:
      ABORT(status, LALINSPIRALH_ESWITCH, LALINSPIRALH_MSGESWITCH);
  }
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void LALBBHPhenWaveFreqDomTemplates(
    LALStatus        *status,
    REAL4Vector      *signalvec1,
    REAL4Vector      *signalvec2,
    InspiralTemplate *params) {
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALBBHPhenWaveFreqDomTemplates");
  switch (params->approximant) {
    case IMRPhenomA:
      if (XLALBBHPhenWaveAFreqDomTemplates(signalvec1, signalvec2, params)) ABORTXLAL(status);
      break;
    case IMRPhenomB:
      if (XLALBBHPhenWaveAFreqDomTemplates(signalvec1, signalvec2, params)) ABORTXLAL(status);
      break;
    default:
      ABORT(status, LALINSPIRALH_ESWITCH, LALINSPIRALH_MSGESWITCH);
  }
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void LALBBHPhenWaveTimeDom(
    LALStatus        *status,
    REAL4Vector      *signalvec1,
    InspiralTemplate *insp_template) {
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALBBHPhenWaveTimeDom");
  if (XLALBBHPhenWaveTimeDom(signalvec1, insp_template)) ABORTXLAL(status);
  CHECKSTATUSPTR(status);
  DETATCHSTATUSPTR(status);
  RETURN (status);
}

void LALBBHPhenWaveTimeDomTemplates(
    LALStatus        *status,
    REAL4Vector      *signalvec1,
    REAL4Vector      *signalvec2,
    InspiralTemplate *insp_template) {
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALBBHPhenWaveTimeDomTemplates");
  if (XLALBBHPhenWaveTimeDomTemplates(signalvec1, signalvec2, insp_template))
    ABORTXLAL(status);
  CHECKSTATUSPTR(status);
  DETATCHSTATUSPTR(status);
  RETURN (status);
}

void LALBBHPhenTimeDomEngine(
    LALStatus        *status,
    REAL4Vector      *signalvec1,
    REAL4Vector      *signalvec2,
    REAL4Vector      *h,
    REAL4Vector      *aVec,
    REAL4Vector      *freqVec,
    REAL8Vector      *phiVec,
    UINT4            *countback,
    InspiralTemplate *params) {
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALBBHPhenTimeDomEngine");
  if (XLALBBHPhenTimeDomEngine(signalvec1, signalvec2, h, aVec,
    freqVec, phiVec, countback, params)) ABORTXLAL(status);
  DETATCHSTATUSPTR(status);
  RETURN (status);
}

void LALBBHPhenWaveTimeDomForInjection(
    LALStatus        *status,
    CoherentGW       *waveform,
    InspiralTemplate *params,
    PPNParamStruc    *ppnParams) {
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);
  XLAL_PRINT_DEPRECATION_WARNING("XLALBBHPhenWaveTimeDomForInjection");
  if (XLALBBHPhenWaveTimeDomForInjection(waveform, params, ppnParams))
    ABORTXLAL(status);
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
