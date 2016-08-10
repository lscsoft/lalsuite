/*
 *  LALInferenceTemplate.c: Bayesian Followup, template calls to LAL
 *  template functions. Temporary GeneratePPN
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever,
 *  Marc van der Sluys, John Veitch, Will M. Farr and Salvatore Vitale
 *
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

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/Date.h>
#include <lal/VectorOps.h>
#include <lal/TimeFreqFFT.h>
#include <lal/GenerateInspiral.h>
#include <lal/GenerateInspRing.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALInference.h>
#include <lal/XLALError.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALInferenceMultibanding.h>

/* LIB imports*/
#include <lal/LALInferenceBurstRoutines.h>


#define PROGRAM_NAME "LALInferenceTemplate.c"
#define CVS_ID_STRING "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define CVS_NAME_STRING "$Name$"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* Max amplitude orders found in LALSimulation (not accessible from outside of LALSim) */
#define MAX_NONPRECESSING_AMP_PN_ORDER 6
#define MAX_PRECESSING_AMP_PN_ORDER 3

#define Pi_p2 9.8696044010893586188344909998761511
#define Pi_p2by3 2.1450293971110256000774441009412356
#define log4 1.3862943611198906188344642429163531

static void q2masses(double mc, double q, double *m1, double *m2);

/* list of testing GR parameters to be passed to the waveform */

const char list_extra_parameters[34][16] = {"dchi0","dchi1","dchi2","dchi3","dchi4","dchi5","dchi5l","dchi6","dchi6l","dchi7","aPPE","alphaPPE","bPPE","betaPPE","betaStep","fStep","dxi1","dxi2","dxi3","dxi4","dxi5","dxi6","dalpha1","dalpha2","dalpha3","dalpha4","dalpha5","dbeta1","dbeta2","dbeta3","dsigma1","dsigma2","dsigma3","dsigma4"};

const UINT4 N_extra_params = 34;


static int InterpolateWaveform(REAL8Vector *freqs, COMPLEX16FrequencySeries *src, COMPLEX16FrequencySeries *dest);
static int InterpolateWaveform(REAL8Vector *freqs, COMPLEX16FrequencySeries *src, COMPLEX16FrequencySeries *dest)
{
  REAL8 deltaF = dest->deltaF;
  UINT4 j=ceil(freqs->data[0] / deltaF);
  COMPLEX16 *d=dest->data->data;
  memset(d, 0, sizeof(*(d))*j);
  
  /* Loop over reduced frequency set */
  for(UINT4 i=0;i<freqs->length-1;i++)
  {  
    double startpsi = carg(src->data->data[i]);
    double startamp = cabs(src->data->data[i]);
    double endpsi = carg(src->data->data[i+1]);
    double endamp = cabs(src->data->data[i+1]);

    double startf=freqs->data[i];
    double endf=freqs->data[i+1];

    double df = endf - startf; /* Big freq step */
    
    /* linear interpolation setup */
    double dpsi = (endpsi - startpsi);

    /* Catch case where phase wraps around */
    /* NOTE: If changing this check that waveforms are not corrupted
     * at high frequencies when dpsi/df can go slightly -ve without
     * the phase wrapping around (e.g. TF2 1.4-1.4 srate=4096)
     */
    if (dpsi/df<-LAL_PI ) {dpsi+=LAL_TWOPI;}
    
    double dpsidf = dpsi/df;
    double dampdf = (endamp - startamp)/df;

    double damp = dampdf *deltaF;
    
    const double dim = sin(dpsidf*deltaF);
    const double dre = 2.0*sin(dpsidf*deltaF*0.5)*sin(dpsidf*deltaF*0.5);

    /* Loop variables */
    double newRe,newIm,f,re,im,a;
    for(f=j*deltaF,
	re = cos(startpsi), im = sin(startpsi),
        a = startamp;
        
        f<endf;
        
        j++, f+=deltaF,
        newRe = re - dre*re-dim*im,
        newIm = im + re*dim-dre*im,
        re=newRe, im = newIm,
        a += damp )
    {
      d[j] = a * (re + I*im);
    }
  }
  memset(&(d[j]), 0, sizeof(d[j])*(dest->data->length - j));
  return 0;
}


void LALInferenceTemplateNullFreqdomain(LALInferenceModel *model)
/**********************************************/
/* returns a frequency-domain 'null' template */
/* (all zeroes, implying no signal present).  */
/**********************************************/
{
  UINT4 i;
  if ((model->freqhPlus==NULL) || (model->freqhCross==NULL)) {
    XLALPrintError(" ERROR in templateNullFreqdomain(): encountered unallocated 'freqModelhPlus/-Cross'.\n");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
  for (i=0; i<model->freqhPlus->data->length; ++i){
    model->freqhPlus->data->data[i] = 0.0;
    model->freqhCross->data->data[i] = 0.0;
  }
  model->domain = LAL_SIM_DOMAIN_FREQUENCY;
  return;
}



void LALInferenceTemplateNullTimedomain(LALInferenceModel *model)
/*********************************************/
/* returns a time-domain 'null' template     */
/* (all zeroes, implying no signal present). */
/*********************************************/
{
  UINT4 i;
  if ((model->timehPlus==NULL) || (model->timehCross==NULL)) {
    XLALPrintError(" ERROR in templateNullTimedomain(): encountered unallocated 'timeModelhPlus/-Cross'.\n");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
  for (i=0; i<model->timehPlus->data->length; ++i){
    model->timehPlus->data->data[i]  = 0.0;
    model->timehCross->data->data[i] = 0.0;
  }
  model->domain = LAL_SIM_DOMAIN_TIME;
  return;
}



/* ============ LAL template wrapper function: ========== */


static void mc2masses(double mc, double eta, double *m1, double *m2);

static void mc2masses(double mc, double eta, double *m1, double *m2)
/*  Compute individual companion masses (m1, m2)   */
/*  for given chirp mass (m_c) & mass ratio (eta)  */
/*  (note: m1 >= m2).                              */
{
  double root = sqrt(0.25-eta);
  double fraction = (0.5+root) / (0.5-root);
  *m2 = mc * (pow(1+fraction,0.2) / pow(fraction,0.6));
  *m1 = mc * (pow(1+1.0/fraction,0.2) / pow(1.0/fraction,0.6));
  return;
}

static void q2masses(double mc, double q, double *m1, double *m2)
/*  Compute individual companion masses (m1, m2)   */
/*  for given chirp mass (m_c) & asymmetric mass   */
/*  ratio (q).  note: q = m2/m1, where m1 >= m2    */
{
  *m1 = mc * pow(q, -3.0/5.0) * pow(q+1, 1.0/5.0);
  *m2 = (*m1) * q;
  return;
}

void LALInferenceROQWrapperForXLALSimInspiralChooseFDWaveformSequence(LALInferenceModel *model){
/*************************************************************************************************************************/
  Approximant approximant = (Approximant) 0;
  INT4 order=-1;
  INT4 amporder;

  int ret=0;
  INT4 errnum=0;

  model->roq->hptildeLinear=NULL, model->roq->hctildeLinear=NULL;
  model->roq->hptildeQuadratic=NULL, model->roq->hctildeQuadratic=NULL;
  REAL8 mc;
  REAL8 phi0, m1, m2, distance, inclination;

  REAL8 *m1_p,*m2_p;

  if (LALInferenceCheckVariable(model->params, "LAL_APPROXIMANT"))
    approximant = *(Approximant*) LALInferenceGetVariable(model->params, "LAL_APPROXIMANT");
  else {
    XLALPrintError(" ERROR in templateLALGenerateInspiral(): (INT4) \"LAL_APPROXIMANT\" parameter not provided!\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }

  if (LALInferenceCheckVariable(model->params, "LAL_PNORDER"))
    order = *(INT4*) LALInferenceGetVariable(model->params, "LAL_PNORDER");
  else {
    XLALPrintError(" ERROR in templateLALGenerateInspiral(): (INT4) \"LAL_PNORDER\" parameter not provided!\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }

  /* Explicitly set the default amplitude order if one is not specified.
   *   This serves two purposes:
   *     1) The default behavior of the code won't change unexpectedly due to changes in LALSimulation.
   *     2) We need to know the amplitude order in order to set the starting frequency of the waveform properly. */
  if (LALInferenceCheckVariable(model->params, "LAL_AMPORDER"))
    amporder = *(INT4*) LALInferenceGetVariable(model->params, "LAL_AMPORDER");
  else
    amporder = -1;
  REAL8 f_ref = 100.0;
  if (LALInferenceCheckVariable(model->params, "f_ref")) f_ref = *(REAL8 *)LALInferenceGetVariable(model->params, "f_ref");

  REAL8 fTemp = f_ref;
  if(LALInferenceCheckVariable(model->params,"chirpmass"))
    {
      mc  = *(REAL8*) LALInferenceGetVariable(model->params, "chirpmass");
      if (LALInferenceCheckVariable(model->params,"q")) {
	REAL8 q = *(REAL8 *)LALInferenceGetVariable(model->params,"q");
	q2masses(mc, q, &m1, &m2);
      } else {
	REAL8 eta = *(REAL8*) LALInferenceGetVariable(model->params, "eta");
	mc2masses(mc, eta, &m1, &m2);
      }
    }
  else if((m1_p=(REAL8 *)LALInferenceGetVariable(model->params, "mass1")) && (m2_p=(REAL8 *)LALInferenceGetVariable(model->params, "mass2")))
    {
      m1=*m1_p;
      m2=*m2_p;
    }
  else
    {
      fprintf(stderr,"No mass parameters found!");
      exit(0);
    }

  distance = exp(LALInferenceGetREAL8Variable(model->params, "logdistance"))* LAL_PC_SI * 1.0e6;        /* distance (1 Mpc) in units of metres */

  phi0 = LALInferenceGetREAL8Variable(model->params, "phase"); /* START phase as per lalsimulation convention, radians*/
  /* Zenith angle between J and N in radians. Also known as inclination angle when spins are aligned */
  REAL8 thetaJN = acos(LALInferenceGetREAL8Variable(model->params, "costheta_jn"));     /* zenith angle between J and N in radians */

  /* ==== SPINS ==== */
  /* We will default to spinless signal and then add in the spin components if required */
  /* If there are non-aligned spins, we must convert between the System Frame coordinates
   * and the cartestian coordinates */

  /* The cartesian spin coordinates (default 0), as passed to LALSimulation */
  REAL8 spin1x = 0.0;
  REAL8 spin1y = 0.0;
  REAL8 spin1z = 0.0;
  REAL8 spin2x = 0.0;
  REAL8 spin2y = 0.0;
  REAL8 spin2z = 0.0;

  /* System frame coordinates as used for jump proposals */
  REAL8 a_spin1 = 0.0;  /* Magnitude of spin1 */
  REAL8 a_spin2 = 0.0;  /* Magnitude of spin2 */
  REAL8 phiJL  = 0.0;  /* azimuthal angle of L_N on its cone about J radians */
  REAL8 tilt1   = 0.0;  /* zenith angle between S1 and LNhat in radians */
  REAL8 tilt2   = 0.0;  /* zenith angle between S2 and LNhat in radians */
  REAL8 phi12   = 0.0;  /* difference in azimuthal angle btwn S1, S2 in radians */

  /* Now check if we have spin amplitudes */
  if(LALInferenceCheckVariable(model->params, "a_spin1"))    a_spin1   = *(REAL8*) LALInferenceGetVariable(model->params, "a_spin1");
  if(LALInferenceCheckVariable(model->params, "a_spin2"))    a_spin2   = *(REAL8*) LALInferenceGetVariable(model->params, "a_spin2");

  /* Check if we have spin angles too */
  if(LALInferenceCheckVariable(model->params, "phi_jl"))
      phiJL = LALInferenceGetREAL8Variable(model->params, "phi_jl");
  if(LALInferenceCheckVariable(model->params, "tilt_spin1"))
      tilt1 = LALInferenceGetREAL8Variable(model->params, "tilt_spin1");
  if(LALInferenceCheckVariable(model->params, "tilt_spin2"))
      tilt2 = LALInferenceGetREAL8Variable(model->params, "tilt_spin2");
  if(LALInferenceCheckVariable(model->params, "phi12"))
      phi12 = LALInferenceGetREAL8Variable(model->params, "phi12");

  /* If we have tilt angles zero, then the spins are aligned and we just set the z component */
  /* However, if the waveform supports precession then we still need to get the right coordinate components */
  SpinSupport spin_support=XLALSimInspiralGetSpinSupportFromApproximant(approximant);
  if(tilt1==0.0 && tilt2==0.0 && (spin_support==LAL_SIM_INSPIRAL_SPINLESS || spin_support==LAL_SIM_INSPIRAL_ALIGNEDSPIN))
  {
      spin1z=a_spin1;
      spin2z=a_spin2;
      inclination = thetaJN; /* Inclination angle is just thetaJN */
  }
  else
  {   /* Template is not aligned-spin only. */
      /* Set all the other spin components according to the angles we received above */
      /* The transformation function doesn't know fLow, so f_ref==0 isn't interpretted as a request to use the starting frequency for reference. */
      XLAL_TRY(ret=XLALSimInspiralTransformPrecessingNewInitialConditions(
                    &inclination, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z,
                    thetaJN, phiJL, tilt1, tilt2, phi12, a_spin1, a_spin2, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, fTemp), errnum);
      if (ret == XLAL_FAILURE)
      {
        XLALPrintError(" ERROR in XLALSimInspiralTransformPrecessingNewInitialConditions(): error converting angles. errnum=%d\n",errnum );
        return;
      }
  }

  /* ==== TIDAL PARAMETERS ==== */
  REAL8 lambda1 = 0.;
  if(LALInferenceCheckVariable(model->params, "lambda1")) lambda1 = *(REAL8*) LALInferenceGetVariable(model->params, "lambda1");
  REAL8 lambda2 = 0.;
  if(LALInferenceCheckVariable(model->params, "lambda2")) lambda2 = *(REAL8*) LALInferenceGetVariable(model->params, "lambda2");
  REAL8 lambdaT = 0.;
  REAL8 dLambdaT = 0.;
  REAL8 sym_mass_ratio_eta = 0.;
  if(LALInferenceCheckVariable(model->params, "lambdaT")&&LALInferenceCheckVariable(model->params, "dLambdaT")){
    lambdaT = *(REAL8*) LALInferenceGetVariable(model->params, "lambdaT");
    dLambdaT = *(REAL8*) LALInferenceGetVariable(model->params, "dLambdaT");
    sym_mass_ratio_eta = m1*m2/((m1+m2)*(m1+m2));
    LALInferenceLambdaTsEta2Lambdas(lambdaT,dLambdaT,sym_mass_ratio_eta,&lambda1,&lambda2);
  }


  /* Only use GR templates */
  LALSimInspiralTestGRParam *nonGRparams = NULL;

  /* Fill in the extra parameters for testing GR, if necessary */
  for (UINT4 k=0; k<N_extra_params; k++)
  {
    if(LALInferenceCheckVariable(model->params,list_extra_parameters[k]))
    {
      XLALSimInspiralAddTestGRParam(&nonGRparams,list_extra_parameters[k],*(REAL8 *)LALInferenceGetVariable(model->params,list_extra_parameters[k]));
    }
  }
  /* Fill in PPE params if they are available */
  char PPEparam[64]="";
  const char *PPEnames[]={"aPPE","alphaPPE","bPPE","betaPPE",NULL};
  for(UINT4 idx=0;PPEnames[idx];idx++)
  {
    for(UINT4 ppeidx=0;;ppeidx++)
    {
      sprintf(PPEparam, "%s%d",PPEnames[idx],ppeidx);
      if(LALInferenceCheckVariable(model->params,PPEparam))
        XLALSimInspiralAddTestGRParam(&nonGRparams,PPEparam,LALInferenceGetREAL8Variable(model->params,PPEparam));
      else
        break;
    }
  }

  /* ==== Call the waveform generator ==== */
  XLAL_TRY(ret=XLALSimInspiralChooseFDWaveformSequence (&(model->roq->hptildeLinear), &(model->roq->hctildeLinear), phi0, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI,
                spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, f_ref, distance, inclination, lambda1, lambda2,
                model->waveFlags, nonGRparams, amporder, order, approximant, (model->roq->frequencyNodesLinear)), errnum);

  XLAL_TRY(ret=XLALSimInspiralChooseFDWaveformSequence (&(model->roq->hptildeQuadratic), &(model->roq->hctildeQuadratic), phi0, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI,
                  spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, f_ref, distance, inclination, lambda1, lambda2,
                  model->waveFlags, nonGRparams, amporder, order, approximant, (model->roq->frequencyNodesQuadratic)), errnum);
    /* Destroy the nonGr params */
    XLALSimInspiralDestroyTestGRParam(nonGRparams);

    REAL8 instant = model->freqhPlus->epoch.gpsSeconds + 1e-9*model->freqhPlus->epoch.gpsNanoSeconds;
    LALInferenceSetVariable(model->params, "time", &instant);

        return;
}

void LALInferenceTemplateSineGaussian(LALInferenceModel *model)
/*****************************************************/
/* Sine-Gaussian (burst) template.                   */
/* Signal is (by now?) linearly polarised,           */
/* i.e., the cross-waveform remains zero.            */
/* * * * * * * * * * * * * * * * * * * * * * * * * * */
/* The (plus-) waveform is:                          */
/*   a * exp(-((t-mu)/sigma)^2) * sin(2*pi*f*t-phi)  */
/* Note that by setting f=0, phi=pi/2 you also get   */
/* a `pure' Gaussian template.                       */
/*                                                   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * ************************************/
/* Required (`model->params') parameters are:                                    */
/*   - "time"       (the "mu" parameter of the Gaussian part; REAL8, GPS sec.)          */
/*   - "sigma"      (width, the "sigma" parameter of the Gaussian part; REAL8, seconds) */
/*   - "frequency"  (frequency of the sine part; REAL8, Hertz)                          */
/*   - "phase"      (phase (at above "mu"); REAL8, radians)                             */
/*   - "amplitude"  (amplitude, REAL8)                                                  */
/****************************************************************************************/
{
  double endtime  = *(REAL8*) LALInferenceGetVariable(model->params, "time");       /* time parameter ("mu"), GPS sec.  */
  double sigma = *(REAL8*) LALInferenceGetVariable(model->params, "sigma");      /* width parameter, seconds         */
  double f     = *(REAL8*) LALInferenceGetVariable(model->params, "frequency");  /* frequency, Hz                    */
  double phi   = *(REAL8*) LALInferenceGetVariable(model->params, "phase");      /* phase, rad                       */
  double a     = *(REAL8*) LALInferenceGetVariable(model->params, "amplitude");  /* amplitude                        */
  double t, tsigma, twopif = LAL_TWOPI*f;
  double epochGPS = XLALGPSGetREAL8(&(model->timehPlus->epoch));
  unsigned long i;
  if (sigma <= 0.0) {
    fprintf(stderr, " ERROR in templateSineGaussian(): zero or negative \"sigma\" parameter (sigma=%e).\n", sigma);
    exit(1);
  }
  if (f < 0.0)
    fprintf(stderr, " WARNING in templateSineGaussian(): negative \"frequency\" parameter (f=%e).\n", f);
  if (a < 0.0)
    fprintf(stderr, " WARNING in templateSineGaussian(): negative \"amplitude\" parameter (a=%e).\n", a);
  for (i=0; i<model->timehPlus->data->length; ++i){
    t = ((double)i)*model->deltaT + (epochGPS-endtime);  /* t-mu         */
    tsigma = t/sigma;                                             /* (t-mu)/sigma */
    if (fabs(tsigma) < 5.0)   /*  (only do computations within a 10 sigma range)  */
      model->timehPlus->data->data[i] = a * exp(-0.5*tsigma*tsigma) * sin(twopif*t+phi);
    else
      model->timehPlus->data->data[i] = 0.0;
    model->timehCross->data->data[i] = 0.0;
  }
  model->domain = LAL_SIM_DOMAIN_TIME;
  return;
}

void LALInferenceTemplateDampedSinusoid(LALInferenceModel *model)
/*****************************************************/
/* Damped Sinusoid (burst) template.                 */
/* Signal is linearly polarized,                     */
/* i.e., cross term is zero.                         */
/* * * * * * * * * * * * * * * * * * * * * * * * * * */
/* The (plus-) waveform is an exponentially decaying */
/* sine wave:                                        */
/*   a * exp((t-time)/tau) * sin(2*pi*f*(t-time))    */
/* where "time" is the time parameter denoting the   */
/* instant where the signal starts.                  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * **************************/
/* Required (`model->params') parameters are:                          */
/*   - "time"       (the instant at which the signal starts; REAL8, GPS sec.) */
/*   - "tau"        (width parameter; REAL8, seconds)                         */
/*   - "frequency"  (frequency of the sine part; REAL8, Hertz)                */
/*   - "amplitude"  (amplitude, REAL8)                                        */
/******************************************************************************/
{
  double endtime  = *(REAL8*) LALInferenceGetVariable(model->params, "time");       /* time parameter ("mu"), GPS sec.  */
  double tau   = *(REAL8*) LALInferenceGetVariable(model->params, "tau");        /* width parameter, seconds         */
  double f     = *(REAL8*) LALInferenceGetVariable(model->params, "frequency");  /* frequency, Hz                    */
  double a     = *(REAL8*) LALInferenceGetVariable(model->params, "amplitude");  /* amplitude                        */
  double t, ttau, twopif = LAL_TWOPI*f;
  double epochGPS = XLALGPSGetREAL8(&(model->timehPlus->epoch));
  unsigned long i;
  if (tau <= 0.0) {
    fprintf(stderr, " ERROR in templateDampedSinusoid(): zero or negative \"tau\" parameter (tau=%e).\n", tau);
    exit(1);
  }
  if (f < 0.0)
    fprintf(stderr, " WARNING in templateDampedSinusoid(): negative \"frequency\" parameter (f=%e).\n", f);
  for (i=0; i<model->timehPlus->data->length; ++i){
    t = ((double)i)*model->deltaT + (epochGPS-endtime);  /* t-mu       */
    if ((t>0.0) && ((ttau=t/tau) < 10.0)) /*  (only do computations within a 10 tau range)  */
      model->timehPlus->data->data[i] = a * exp(-ttau) * sin(twopif*t);
    else
      model->timehPlus->data->data[i] = 0.0;
    model->timehCross->data->data[i] = 0.0;
  }
  model->domain = LAL_SIM_DOMAIN_TIME;
  return;
}



void LALInferenceTemplateSinc(LALInferenceModel *model)
/*****************************************************/
/* Sinc function (burst) template.                   */
/* Signal is linearly polarized,                     */
/* i.e., cross term is zero.                         */
/* * * * * * * * * * * * * * * * * * * * * * * * * * */
/* The (plus-) waveform is a sinc function of given  */
/* frequency:                                        */
/*   a * sinc(2*pi*f*(t-time))                       */
/*   = a * sin(2*pi*f*(t-time)) / (2*pi*f*(t-time))  */
/* where "time" is the time parameter denoting the   */
/* signal's central peak location.                   */
/* * * * * * * * * * * * * * * * * * * * * * * * * * *************************/
/* Required (`model->params') parameters are:                         */
/*   - "time"       (the instant at which the signal peaks; REAL8, GPS sec.) */
/*   - "frequency"  (frequency of the sine part; REAL8, Hertz)               */
/*   - "amplitude"  (amplitude, REAL8)                                       */
/*****************************************************************************/
{
  double endtime  = *(REAL8*) LALInferenceGetVariable(model->params, "time");       /* time parameter ("mu"), GPS sec.  */
  double f     = *(REAL8*) LALInferenceGetVariable(model->params, "frequency");  /* frequency, Hz                    */
  double a     = *(REAL8*) LALInferenceGetVariable(model->params, "amplitude");  /* amplitude                        */
  double t, sinArg, sinc, twopif = LAL_TWOPI*f;
  double epochGPS = XLALGPSGetREAL8(&(model->timehPlus->epoch));
  unsigned long i;
  if (f < 0.0)
    fprintf(stderr, " WARNING in templateSinc(): negative \"frequency\" parameter (f=%e).\n", f);
  for (i=0; i<model->timehPlus->data->length; ++i){
    t = ((double)i)*model->deltaT + (epochGPS-endtime);  /* t-mu       */
    sinArg = twopif*t;
    sinc = (sinArg==0.0) ? 1.0 : sin(sinArg)/sinArg;
    model->timehPlus->data->data[i] = a * sinc;
    model->timehCross->data->data[i] = 0.0;
  }
  model->domain = LAL_SIM_DOMAIN_TIME;
  return;
}


void LALInferenceTemplateASinOmegaT(LALInferenceModel *model)
/************************************************************/
/* Trivial h(t)=A*sin(Omega*t) template						*/
/*  Required (`model->params') parameters are:       */
/*   - "A"       (dimensionless amplitude, REAL8)			*/
/*   - "Omega"   (frequency; REAL8, radians/sec)            */
/************************************************************/
{
  double A		= *(REAL8*) LALInferenceGetVariable(model->params, "A");				/* dim-less	   */
  double Omega	= *(REAL8*) LALInferenceGetVariable(model->params, "Omega");			/* rad/sec     */
  double t;
  double epochGPS = XLALGPSGetREAL8(&(model->timehPlus->epoch));

  unsigned long i;
  for (i=0; i<model->timehPlus->data->length; ++i){
    t = ((double)i)*model->deltaT + (epochGPS);  /* t-mu       */
    model->timehPlus->data->data[i] = A * sin(Omega*t);
    model->timehCross->data->data[i] = 0.0;
  }
  model->domain = LAL_SIM_DOMAIN_TIME;
  return;
}

void LALInferenceTemplateXLALSimInspiralChooseWaveform(LALInferenceModel *model)
/*************************************************************************************************************************/
/* Wrapper for LALSimulation waveforms:						                                                             */
/* XLALSimInspiralChooseFDWaveform() and XLALSimInspiralChooseTDWaveform().                                              */
/*                                                                                                                       */
/*  model->params parameters are:										                                                 */
/*  - "name" description; type OPTIONAL (default value)										                             */
/*										                                                                                 */
/*   MODEL PARAMETERS										                                                             */
/*   - "LAL_APPROXIMANT     Approximant;        Approximant                                                              */
/*   - "LAL_PNORDER"        Phase PN order;     INT4                                                                     */
/*   - "LAL_AMPORDER"       Amplitude PN order; INT4 OPTIONAL (-1)                                                       */
/*   - "spinO"              Spin order;         LALSimInspiralSpinOrder OPTIONAL (LAL_SIM_INSPIRAL_SPIN_ORDER_DEFAULT)   */
/*   - "tideO"              Tidal order;        LALSimInspiralTidalOrder OPTIONAL (LAL_SIM_INSPIRAL_TIDAL_ORDER_DEFAULT) */
/*   - "f_ref"               frequency at which the (frequency dependent) parameters are defined; REAL8 OPTIONAL (0.0)    */
/*   - "fLow"               lower frequency bound; REAL8 OPTIONAL (model->fLow)                                          */
/*                                                                                                                       */
/*   MASS PARAMETERS; either:                                                                                            */
/*      - "mass1"           mass of object 1 in solar mass; REAL8								                         */
/*      - "mass2"		        mass of object 1 in solar mass; REAL8								                     */
/*      OR                                                                                                               */
/*      - "chirpmass"       chirpmass in solar mass; REAL8                                                               */
/*      - "q"  asymmetric mass ratio m2/m1, 0<q<1; REAL8                                      */
/*      OR                                                                                                               */
/*      - "chirpmass"       chirpmass in solar mass; REAL8                                                               */
/*      - "eta"             symmetric mass ratio (m1*m2)/(m1+m2)^2; REAL8                                                */
/*                                                                                                                       */
/*   ORIENTATION AND SPIN PARAMETERS                                                                                     */
/*   - "phi0"               reference phase as per LALSimulation convention; REAL8                                       */
/*   - "logdistance"           distance in Mpc                                                                              */
/*   - "costheta_jn");      cos of zenith angle between J and N in radians;            REAL8                                    */
/*   - "phi_jl");        azimuthal angle of L_N on its cone about J radians; REAL8                                    */
/*   - "tilt_spin1");    zenith angle between S1 and LNhat in radians;       REAL8                                    */
/*   - "tilt_spin2");    zenith angle between S2 and LNhat in radians;       REAL8                                    */
/*   - "phi12");         difference in azimuthal angle between S1, S2 in radians;   REAL8                             */
/*   - "a_spin1"            magnitude of spin 1 in general configuration, -1<a_spin1<1; REAL8 OPTIONAL (0.0)              */
/*   - "a_spin2"            magnitude of spin 2 in general configuration, -1<a_spin1<1; REAL8 OPTIONAL (0.0)              */
/*                                                                                                                       */
/*   OTHER PARAMETERS                                                                                                    */
/*   - "lambda1"            tidal parameter of object 1; REAL8  OPTIONAL (0.0)                                           */
/*   - "lambda2"            tidal parameter of object 1; REAL8  OPTIONAL (0.0)                                           */
/*                                                                                                                       */
/*   - "time"               used as an OUTPUT only; REAL8								                                 */
/*                                                                                                                       */
/*                                                                                                                       */
/*   model needs to also contain:                                                                                        */
/*   - model->fLow Unless  - "fLow" OPTIONAL                                                                             */
/*   - model->deltaT                                                                                                     */
/*   - if model->domain == LAL_SIM_DOMAIN_FREQUENCY                                                                      */
/*      - model->deltaF                                                                                                  */
/*      - model->freqhCross                                                                                              */
/*      - model->freqhPlus                                                                                               */
/*   - else                                                                                                              */
/*      - model->timehPlus                                                                                               */
/*      - model->timehCross                                                                                              */
/*************************************************************************************************************************/
{

  Approximant approximant = (Approximant) 0;
  INT4 order=-1;
  INT4 amporder;

  static int sizeWarning = 0;
  int ret=0;
  INT4 errnum=0;

  REAL8TimeSeries *hplus=NULL;  /**< +-polarization waveform [returned] */
  REAL8TimeSeries *hcross=NULL; /**< x-polarization waveform [returned] */
  COMPLEX16FrequencySeries *hptilde=NULL, *hctilde=NULL;
  REAL8 mc;
  REAL8 phi0, deltaT, m1, m2, f_low, f_start, distance, inclination;

  REAL8 *m1_p,*m2_p;
  REAL8 deltaF, f_max;

  /* Sampling rate for time domain models */
  deltaT = model->deltaT;

  if (LALInferenceCheckVariable(model->params, "LAL_APPROXIMANT"))
    approximant = *(Approximant*) LALInferenceGetVariable(model->params, "LAL_APPROXIMANT");
  else {
    XLALPrintError(" ERROR in templateLALGenerateInspiral(): (INT4) \"LAL_APPROXIMANT\" parameter not provided!\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }

  if (LALInferenceCheckVariable(model->params, "LAL_PNORDER"))
    order = *(INT4*) LALInferenceGetVariable(model->params, "LAL_PNORDER");
  else {
    XLALPrintError(" ERROR in templateLALGenerateInspiral(): (INT4) \"LAL_PNORDER\" parameter not provided!\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }

  /* Explicitly set the default amplitude order if one is not specified.
   *   This serves two purposes:
   *     1) The default behavior of the code won't change unexpectedly due to changes in LALSimulation.
   *     2) We need to know the amplitude order in order to set the starting frequency of the waveform properly. */
  if (LALInferenceCheckVariable(model->params, "LAL_AMPORDER"))
    amporder = *(INT4*) LALInferenceGetVariable(model->params, "LAL_AMPORDER");
  else
    amporder = -1;

  REAL8 f_ref = 100.0;
  if (LALInferenceCheckVariable(model->params, "f_ref")) f_ref = *(REAL8 *)LALInferenceGetVariable(model->params, "f_ref");

  REAL8 fTemp = f_ref;

  if(LALInferenceCheckVariable(model->params,"chirpmass"))
    {
      mc  = *(REAL8*) LALInferenceGetVariable(model->params, "chirpmass");
      if (LALInferenceCheckVariable(model->params,"q")) {
	REAL8 q = *(REAL8 *)LALInferenceGetVariable(model->params,"q");
	q2masses(mc, q, &m1, &m2);
      } else {
	REAL8 eta = *(REAL8*) LALInferenceGetVariable(model->params, "eta");
	mc2masses(mc, eta, &m1, &m2);
      }
    }
  else if((m1_p=(REAL8 *)LALInferenceGetVariable(model->params, "mass1")) && (m2_p=(REAL8 *)LALInferenceGetVariable(model->params, "mass2")))
    {
      m1=*m1_p;
      m2=*m2_p;
    }
  else
    {
      fprintf(stderr,"No mass parameters found!");
      exit(0);
    }

  distance	= exp(LALInferenceGetREAL8Variable(model->params, "logdistance"))* LAL_PC_SI * 1.0e6;        /* distance (1 Mpc) in units of metres */

  phi0		= LALInferenceGetREAL8Variable(model->params, "phase"); /* START phase as per lalsimulation convention, radians*/

  /* Zenith angle between J and N in radians. Also known as inclination angle when spins are aligned */
  REAL8 thetaJN = acos(LALInferenceGetREAL8Variable(model->params, "costheta_jn"));     /* zenith angle between J and N in radians */

  /* Check if fLow is a model parameter, otherwise use data structure definition */
  if(LALInferenceCheckVariable(model->params, "flow"))
    f_low = *(REAL8*) LALInferenceGetVariable(model->params, "flow");
  else
    f_low = model->fLow;

  f_start = XLALSimInspiralfLow2fStart(f_low, amporder, approximant);
  f_max = 0.0; /* for freq domain waveforms this will stop at ISCO. Previously found using model->fHigh causes NaNs in waveform (see redmine issue #750)*/

  /* ==== SPINS ==== */
  /* We will default to spinless signal and then add in the spin components if required */
  /* If there are non-aligned spins, we must convert between the System Frame coordinates
   * and the cartestian coordinates */

  /* The cartesian spin coordinates (default 0), as passed to LALSimulation */
  REAL8 spin1x = 0.0;
  REAL8 spin1y = 0.0;
  REAL8 spin1z = 0.0;
  REAL8 spin2x = 0.0;
  REAL8 spin2y = 0.0;
  REAL8 spin2z = 0.0;

  /* System frame coordinates as used for jump proposals */
  REAL8 a_spin1 = 0.0;  /* Magnitude of spin1 */
  REAL8 a_spin2 = 0.0;  /* Magnitude of spin2 */
  REAL8 phiJL  = 0.0;  /* azimuthal angle of L_N on its cone about J radians */
  REAL8 tilt1   = 0.0;  /* zenith angle between S1 and LNhat in radians */
  REAL8 tilt2   = 0.0;  /* zenith angle between S2 and LNhat in radians */
  REAL8 phi12   = 0.0;  /* difference in azimuthal angle btwn S1, S2 in radians */

  /* Now check if we have spin amplitudes */
  if(LALInferenceCheckVariable(model->params, "a_spin1"))    a_spin1   = *(REAL8*) LALInferenceGetVariable(model->params, "a_spin1");
  if(LALInferenceCheckVariable(model->params, "a_spin2"))    a_spin2   = *(REAL8*) LALInferenceGetVariable(model->params, "a_spin2");

  /* Check if we have spin angles too */
  if(LALInferenceCheckVariable(model->params, "phi_jl"))
      phiJL = LALInferenceGetREAL8Variable(model->params, "phi_jl");
  if(LALInferenceCheckVariable(model->params, "tilt_spin1"))
      tilt1 = LALInferenceGetREAL8Variable(model->params, "tilt_spin1");
  if(LALInferenceCheckVariable(model->params, "tilt_spin2"))
      tilt2 = LALInferenceGetREAL8Variable(model->params, "tilt_spin2");
  if(LALInferenceCheckVariable(model->params, "phi12"))
      phi12 = LALInferenceGetREAL8Variable(model->params, "phi12");

  /* If we have tilt angles zero, then the spins are aligned and we just set the z component */
  /* However, if the waveform supports precession then we still need to get the right coordinate components */
  SpinSupport spin_support=XLALSimInspiralGetSpinSupportFromApproximant(approximant);
  if(tilt1==0.0 && tilt2==0.0 && (spin_support==LAL_SIM_INSPIRAL_SPINLESS || spin_support==LAL_SIM_INSPIRAL_ALIGNEDSPIN))
  {
      spin1z=a_spin1;
      spin2z=a_spin2;
      inclination = thetaJN; /* Inclination angle is just thetaJN */
  }
  else
  {   /* Template is not aligned-spin only. */
      /* Set all the other spin components according to the angles we received above */
      /* The transformation function doesn't know fLow, so f_ref==0 isn't interpretted as a request to use the starting frequency for reference. */
      if(fTemp==0.0)
        fTemp = f_start;

      XLAL_TRY(ret=XLALSimInspiralTransformPrecessingNewInitialConditions(
                    &inclination, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z,
                    thetaJN, phiJL, tilt1, tilt2, phi12, a_spin1, a_spin2, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, fTemp), errnum);
      if (ret == XLAL_FAILURE)
      {
        XLALPrintError(" ERROR in XLALSimInspiralTransformPrecessingNewInitialConditions(): error converting angles. errnum=%d: %s\n",errnum, XLALErrorString(errnum) );
        return;
      }
  }


  /* ==== TIDAL PARAMETERS ==== */
  REAL8 lambda1 = 0.;
  if(LALInferenceCheckVariable(model->params, "lambda1")) lambda1 = *(REAL8*) LALInferenceGetVariable(model->params, "lambda1");
  REAL8 lambda2 = 0.;
  if(LALInferenceCheckVariable(model->params, "lambda2")) lambda2 = *(REAL8*) LALInferenceGetVariable(model->params, "lambda2");
  REAL8 lambdaT = 0.;
  REAL8 dLambdaT = 0.;
  REAL8 sym_mass_ratio_eta = 0.;
  if(LALInferenceCheckVariable(model->params, "lambdaT")&&LALInferenceCheckVariable(model->params, "dLambdaT")){
    lambdaT = *(REAL8*) LALInferenceGetVariable(model->params, "lambdaT");
    dLambdaT = *(REAL8*) LALInferenceGetVariable(model->params, "dLambdaT");
    sym_mass_ratio_eta = m1*m2/((m1+m2)*(m1+m2));
    LALInferenceLambdaTsEta2Lambdas(lambdaT,dLambdaT,sym_mass_ratio_eta,&lambda1,&lambda2);
  }


  /* Only use GR templates */
  LALSimInspiralTestGRParam *nonGRparams = NULL;

  /* Fill in the extra parameters for testing GR, if necessary */
  for (UINT4 k=0; k<N_extra_params; k++)
  {
      if(LALInferenceCheckVariable(model->params,list_extra_parameters[k]))
      {
          XLALSimInspiralAddTestGRParam(&nonGRparams,list_extra_parameters[k],*(REAL8 *)LALInferenceGetVariable(model->params,list_extra_parameters[k]));
      }
  }
  /* Fill in PPE params if they are available */
  char PPEparam[64]="";
  const char *PPEnames[]={"aPPE","alphaPPE","bPPE","betaPPE",NULL};
  for(UINT4 idx=0;PPEnames[idx];idx++)
  {
    for(UINT4 ppeidx=0;;ppeidx++)
    {
      sprintf(PPEparam, "%s%d",PPEnames[idx],ppeidx);
      if(LALInferenceCheckVariable(model->params,PPEparam))
        XLALSimInspiralAddTestGRParam(&nonGRparams,PPEparam,LALInferenceGetREAL8Variable(model->params,PPEparam));
      else
        break;
    }
  }



  /* ==== Call the waveform generator ==== */
  if(model->domain == LAL_SIM_DOMAIN_FREQUENCY) {
    deltaF = model->deltaF;

	XLAL_TRY(ret=XLALSimInspiralChooseFDWaveformFromCache(&hptilde, &hctilde, phi0,
            deltaF, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z,
            spin2x, spin2y, spin2z, f_start, f_max, f_ref, distance, inclination,lambda1, lambda2, model->waveFlags, nonGRparams, amporder, order,
            approximant,model->waveformCache, NULL), errnum);


    /* if the waveform failed to generate, fill the buffer with zeros
     * so that the previous waveform is not left there
     */
    if(ret!=XLAL_SUCCESS){
      memset(model->freqhPlus->data->data,0,sizeof(model->freqhPlus->data->data[0])*model->freqhPlus->data->length);
      memset(model->freqhCross->data->data,0,sizeof(model->freqhCross->data->data[0])*model->freqhCross->data->length);
      if ( hptilde ) XLALDestroyCOMPLEX16FrequencySeries(hptilde);
      if ( hctilde ) XLALDestroyCOMPLEX16FrequencySeries(hctilde);
      XLALSimInspiralDestroyTestGRParam(nonGRparams);
      errnum&=~XLAL_EFUNC; /* Mask out the internal function failure bit */
      switch(errnum)
      {
        case XLAL_EDOM:
          /* The waveform was called outside its domain. Return an empty vector but not an error */
          XLAL_ERROR_VOID(XLAL_EUSR0);
        default:
          /* Another error occurred that we can't handle. Propogate upward */
          XLALSetErrno(errnum);
          XLAL_ERROR_VOID(errnum,"%s: Template generation failed in XLALSimInspiralChooseFDWaveformFromCache:\n\
XLALSimInspiralChooseFDWaveformFromCache(&hptilde, &hctilde, \
%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, \
model->waveFlags(%d,%d,%d,%d,numreldata),nonGRparams,%d,%d,%d,model->waveformCache)\n",__func__,
            phi0, deltaF, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
            f_start, f_max, f_ref, distance, inclination, lambda1, lambda2, (int) XLALSimInspiralGetSpinOrder(model->waveFlags),
            (int) XLALSimInspiralGetTidalOrder(model->waveFlags),(int) XLALSimInspiralGetFrameAxis(model->waveFlags),
            (int) XLALSimInspiralGetModesChoice(model->waveFlags),amporder, order, approximant);
      }
    }

    if (hptilde==NULL || hptilde->data==NULL || hptilde->data->data==NULL ) {
	  XLALPrintError(" ERROR in %s: encountered unallocated 'hptilde'.\n",__func__);
	  XLAL_ERROR_VOID(XLAL_EFAULT);
	}
	if (hctilde==NULL || hctilde->data==NULL || hctilde->data->data==NULL ) {
	  XLALPrintError(" ERROR in %s: encountered unallocated 'hctilde'.\n",__func__);
	  XLAL_ERROR_VOID(XLAL_EFAULT);
	}

    INT4 rem=0;
    UINT4 size=hptilde->data->length;
    if(size>model->freqhPlus->data->length) size=model->freqhPlus->data->length;
    memcpy(model->freqhPlus->data->data,hptilde->data->data,sizeof(hptilde->data->data[0])*size);
    if( (rem=(model->freqhPlus->data->length - size)) > 0)
        memset(&(model->freqhPlus->data->data[size]),0, rem*sizeof(hptilde->data->data[0]) );

    size=hctilde->data->length;
    if(size>model->freqhCross->data->length) size=model->freqhCross->data->length;
    memcpy(model->freqhCross->data->data,hctilde->data->data,sizeof(hctilde->data->data[0])*size);
    if( (rem=(model->freqhCross->data->length - size)) > 0)
        memset(&(model->freqhCross->data->data[size]),0, rem*sizeof(hctilde->data->data[0]) );

    REAL8 instant = model->freqhPlus->epoch.gpsSeconds + 1e-9*model->freqhPlus->epoch.gpsNanoSeconds;
    LALInferenceSetVariable(model->params, "time", &instant);

  } else {

    XLAL_TRY(ret=XLALSimInspiralChooseTDWaveformFromCache(&hplus, &hcross, phi0, deltaT,
            m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z,
            spin2x, spin2y, spin2z, f_start, f_ref, distance,
            inclination, lambda1, lambda2, model->waveFlags, nonGRparams,
            amporder, order, approximant,model->waveformCache), errnum);
    /* if the waveform failed to generate, fill the buffer with zeros
     * so that the previous waveform is not left there
     */
    if(ret!=XLAL_SUCCESS){
      memset(model->freqhPlus->data->data,0,sizeof(model->freqhPlus->data->data[0])*model->freqhPlus->data->length);
      memset(model->freqhCross->data->data,0,sizeof(model->freqhCross->data->data[0])*model->freqhCross->data->length);
      if ( hptilde ) XLALDestroyCOMPLEX16FrequencySeries(hptilde);
      if ( hctilde ) XLALDestroyCOMPLEX16FrequencySeries(hctilde);
      if ( hplus) XLALDestroyREAL8TimeSeries(hplus);
      if ( hcross ) XLALDestroyREAL8TimeSeries(hcross);
      XLALSimInspiralDestroyTestGRParam(nonGRparams);
      errnum&=~XLAL_EFUNC; /* Mask out the internal function failure bit */
      switch(errnum)
      {
        case XLAL_EDOM:
          /* The waveform was called outside its domain. Return an empty vector but not an error */
          XLAL_ERROR_VOID(XLAL_EUSR0);
        default:
          /* Another error occurred that we can't handle. Propogate upward */
          XLALSetErrno(errnum);
          XLAL_ERROR_VOID(errnum,"%s: Template generation failed in XLALSimInspiralChooseTDWaveformFromCache\n\
XLALSimInspiralChooseTDWaveformFromCache(&hplus, &hcross, \
%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, \
model->waveFlags(%d,%d,%d,%d,numreldata),nonGRparams,%d,%d,%d,model->waveformCache)\n",__func__,
            phi0, deltaT, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
            f_start, f_ref, distance, inclination, lambda1, lambda2, (int) XLALSimInspiralGetSpinOrder(model->waveFlags),
            (int) XLALSimInspiralGetTidalOrder(model->waveFlags),(int) XLALSimInspiralGetFrameAxis(model->waveFlags),
            (int) XLALSimInspiralGetModesChoice(model->waveFlags),amporder, order, approximant);
      }
    }

    /* The following complicated mess is a result of the following considerations:

       1) The discrete time samples of the template and the timeModel
       buffers will not, in general line up.

       2) The likelihood function will timeshift the template in the
       frequency domain to align it properly with the desired tc in
       each detector (these are different because the detectors
       receive the signal at different times).  Because this
       timeshifting is done in the frequency domain, the effective
       time-domain template is periodic.  We want to avoid the
       possibility of non-zero template samples wrapping around from
       the start/end of the buffer, since real templates are not
       periodic!

       3) If the template apporaches the ends of the timeModel buffer,
       then it should be tapered in the same way as the timeData
       (currently 0.4 seconds, hard-coded! Tukey window; see
       LALInferenceReadData.c, near line 233) so that template and
       signal in the data match.  However, as an optimization, we
       perform only one tapering and FFT-ing in the likelihood
       function; subsequent timeshifts for the different detectors
       will cause the tapered regions of the template and data to
       become mis-aligned.

       The algorthim we use is the following:

       1) Inject the template to align with the nearest sample in the
       timeModel buffer to the desired geocent_end time.

       2) Check whether either the start or the end of the template
       overlaps the tapered region, plus a safety buffer corresponding
       to a conservative estimate of the largest geocenter <-->
       detector timeshift.

         a) If there is no overlap at the start or end of the buffer,
         we're done.

	 b) If there is an overlap, issue one warning per process
	 (which can be disabled by setting the LAL debug level) about
	 a too-short segment length, and return.
*/

    size_t waveLength = hplus->data->length;
    size_t bufLength = model->timehPlus->data->length;

    /* 2*Rearth/(c*deltaT)---2 is safety factor---is the maximum time
       shift for any earth-based detector. */
    size_t maxShift = (size_t)lround(4.255e-2/deltaT);

    /* Taper 0.4 seconds at start and end (hard-coded! in
       LALInferenceReadData.c, around line 233). */
    size_t taperLength = (size_t)lround(0.4/deltaT);

    /* Within unsafeLength of ends of buffer, possible danger of
       wrapping and/or tapering interactions. */
    size_t unsafeLength = taperLength + maxShift;

    REAL8 desiredTc = *(REAL8 *)LALInferenceGetVariable(model->params, "time");
    REAL8 tStart = XLALGPSGetREAL8(&(model->timehPlus->epoch));
    REAL8 tEnd = tStart + deltaT * model->timehPlus->data->length;

    if (desiredTc < tStart || desiredTc > tEnd) {
      XLALDestroyREAL8TimeSeries(hplus);
      XLALDestroyREAL8TimeSeries(hcross);

      XLAL_PRINT_ERROR("desired tc (%.4f) outside data buffer\n", desiredTc);
      XLAL_ERROR_VOID(XLAL_EDOM);
    }

    /* The nearest sample in model buffer to the desired tc. */
    size_t tcSample = (size_t)lround((desiredTc - XLALGPSGetREAL8(&(model->timehPlus->epoch)))/deltaT);

    /* The acutal coalescence time that corresponds to the buffer
       sample on which the waveform's tC lands. */
    REAL8 injTc = XLALGPSGetREAL8(&(model->timehPlus->epoch)) + tcSample*deltaT;

    /* The sample at which the waveform reaches tc. */
    size_t waveTcSample = (size_t)lround(-XLALGPSGetREAL8(&(hplus->epoch))/deltaT);

    /* 1 + (number of samples post-tc in waveform) */
    size_t wavePostTc = waveLength - waveTcSample;

    size_t bufStartIndex = (tcSample >= waveTcSample ? tcSample - waveTcSample : 0);
    size_t bufEndIndex = (wavePostTc + tcSample <= bufLength ? wavePostTc + tcSample : bufLength);
    size_t bufWaveLength = bufEndIndex - bufStartIndex;
    size_t waveStartIndex = (tcSample >= waveTcSample ? 0 : waveTcSample - tcSample);

    if (bufStartIndex < unsafeLength || (bufLength - bufEndIndex) <= unsafeLength) {
      /* The waveform could be timeshifted into a region where it will
	 be tapered improperly, or even wrap around from the periodic
	 timeshift.  Issue warning. */
      if (!sizeWarning) {
	fprintf(stderr, "WARNING: Generated template is too long to guarantee that it will not\n");
	fprintf(stderr, "WARNING:  (a) lie in a tapered region of the time-domain buffer\n");
	fprintf(stderr, "WARNING:  (b) wrap periodically when timeshifted in likelihood computation\n");
	fprintf(stderr, "WARNING: Either of these may cause differences between the template and the\n");
	fprintf(stderr, "WARNING: correct GW waveform in each detector.\n");
	fprintf(stderr, "WARNING: Parameter estimation will continue, but you should consider\n");
	fprintf(stderr, "WARNING: increasing the data segment length (using the --seglen) option.\n");
	sizeWarning = 1;
      }
    }

    /* Clear model buffers */
    memset(model->timehPlus->data->data, 0, sizeof(REAL8)*model->timehPlus->data->length);
    memset(model->timehCross->data->data, 0, sizeof(REAL8)*model->timehCross->data->length);

    /* Inject */
    memcpy(model->timehPlus->data->data + bufStartIndex,
	   hplus->data->data + waveStartIndex,
	   bufWaveLength*sizeof(REAL8));
    memcpy(model->timehCross->data->data + bufStartIndex,
	   hcross->data->data + waveStartIndex,
	   bufWaveLength*sizeof(REAL8));

    LALInferenceSetVariable(model->params, "time", &injTc);
  }

  if ( hplus ) XLALDestroyREAL8TimeSeries(hplus);
  if ( hcross ) XLALDestroyREAL8TimeSeries(hcross);
  if ( hptilde ) XLALDestroyCOMPLEX16FrequencySeries(hptilde);
  if ( hctilde ) XLALDestroyCOMPLEX16FrequencySeries(hctilde);

  return;
}


void LALInferenceTemplateXLALSimInspiralChooseWaveformPhaseInterpolated(LALInferenceModel *model)
/*************************************************************************************************************************/
/* Wrapper for LALSimulation waveforms:						                                                             */
/* XLALSimInspiralChooseFDWaveform() and XLALSimInspiralChooseTDWaveform().                                              */
/*                                                                                                                       */
/*  model->params parameters are:										                                                 */
/*  - "name" description; type OPTIONAL (default value)										                             */
/*										                                                                                 */
/*   MODEL PARAMETERS										                                                             */
/*   - "LAL_APPROXIMANT     Approximant;        Approximant                                                              */
/*   - "LAL_PNORDER"        Phase PN order;     INT4                                                                     */
/*   - "LAL_AMPORDER"       Amplitude PN order; INT4 OPTIONAL (-1)                                                       */
/*   - "spinO"              Spin order;         LALSimInspiralSpinOrder OPTIONAL (LAL_SIM_INSPIRAL_SPIN_ORDER_DEFAULT)   */
/*   - "tideO"              Tidal order;        LALSimInspiralTidalOrder OPTIONAL (LAL_SIM_INSPIRAL_TIDAL_ORDER_DEFAULT) */
/*   - "f_ref"               frequency at which the (frequency dependent) parameters are defined; REAL8 OPTIONAL (0.0)    */
/*   - "fLow"               lower frequency bound; REAL8 OPTIONAL (model->fLow)                                          */
/*                                                                                                                       */
/*   MASS PARAMETERS; either:                                                                                            */
/*      - "mass1"           mass of object 1 in solar mass; REAL8								                         */
/*      - "mass2"		        mass of object 1 in solar mass; REAL8								                     */
/*      OR                                                                                                               */
/*      - "chirpmass"       chirpmass in solar mass; REAL8                                                               */
/*      - "q"  asymmetric mass ratio m2/m1, 0<q<1; REAL8                                      */
/*      OR                                                                                                               */
/*      - "chirpmass"       chirpmass in solar mass; REAL8                                                               */
/*      - "eta"             symmetric mass ratio (m1*m2)/(m1+m2)^2; REAL8                                                */
/*                                                                                                                       */
/*   ORIENTATION AND SPIN PARAMETERS                                                                                     */
/*   - "phi0"               reference phase as per LALSimulation convention; REAL8                                       */
/*   - "logdistance"           distance in Mpc                                                                              */
/*   - "costheta_jn");      cos of zenith angle between J and N in radians;            REAL8                                    */
/*   - "phi_jl");        azimuthal angle of L_N on its cone about J radians; REAL8                                    */
/*   - "tilt_spin1");    zenith angle between S1 and LNhat in radians;       REAL8                                    */
/*   - "tilt_spin2");    zenith angle between S2 and LNhat in radians;       REAL8                                    */
/*   - "phi12");         difference in azimuthal angle between S1, S2 in radians;   REAL8                             */
/*   - "a_spin1"            magnitude of spin 1 in general configuration, -1<a_spin1<1; REAL8 OPTIONAL (0.0)              */
/*   - "a_spin2"            magnitude of spin 2 in general configuration, -1<a_spin1<1; REAL8 OPTIONAL (0.0)              */
/*                                                                                                                       */
/*   OTHER PARAMETERS                                                                                                    */
/*   - "lambda1"            tidal parameter of object 1; REAL8  OPTIONAL (0.0)                                           */
/*   - "lambda2"            tidal parameter of object 1; REAL8  OPTIONAL (0.0)                                           */
/*                                                                                                                       */
/*   - "time"               used as an OUTPUT only; REAL8								                                 */
/*                                                                                                                       */
/*                                                                                                                       */
/*   model needs to also contain:                                                                                        */
/*   - model->fLow Unless  - "fLow" OPTIONAL                                                                             */
/*   - model->deltaT                                                                                                     */
/*   - if model->domain == LAL_SIM_DOMAIN_FREQUENCY                                                                      */
/*      - model->deltaF                                                                                                  */
/*      - model->freqhCross                                                                                              */
/*      - model->freqhPlus                                                                                               */
/*   - else                                                                                                              */
/*      - model->timehPlus                                                                                               */
/*      - model->timehCross                                                                                              */
/*************************************************************************************************************************/
{
    Approximant approximant = (Approximant) 0;
    INT4 order=-1;
    INT4 amporder;

    static int sizeWarning = 0;
    int ret=0;
    INT4 errnum=0;
    REAL8TimeSeries *hplus=NULL;  /**< +-polarization waveform [returned] */
    REAL8TimeSeries *hcross=NULL; /**< x-polarization waveform [returned] */
    COMPLEX16FrequencySeries *hptilde=NULL, *hctilde=NULL;
    REAL8 mc;
    REAL8 phi0, deltaT, m1, m2, f_low, f_start, distance, inclination;

    REAL8 *m1_p,*m2_p;
    REAL8 f_max;

    /* Sampling rate for time domain models */
    deltaT = model->deltaT;

    if (LALInferenceCheckVariable(model->params, "LAL_APPROXIMANT"))
        approximant = *(Approximant*) LALInferenceGetVariable(model->params, "LAL_APPROXIMANT");
    else {
        XLALPrintError(" ERROR in templateLALGenerateInspiral(): (INT4) \"LAL_APPROXIMANT\" parameter not provided!\n");
        XLAL_ERROR_VOID(XLAL_EDATA);
    }

    if (LALInferenceCheckVariable(model->params, "LAL_PNORDER"))
        order = *(INT4*) LALInferenceGetVariable(model->params, "LAL_PNORDER");
    else {
        XLALPrintError(" ERROR in templateLALGenerateInspiral(): (INT4) \"LAL_PNORDER\" parameter not provided!\n");
        XLAL_ERROR_VOID(XLAL_EDATA);
    }

    /* Explicitly set the default amplitude order if one is not specified.
     *   This serves two purposes:
     *     1) The default behavior of the code won't change unexpectedly due to changes in LALSimulation.
     *     2) We need to know the amplitude order in order to set the starting frequency of the waveform properly. */

    if (LALInferenceCheckVariable(model->params, "LAL_AMPORDER"))
        amporder = *(INT4*) LALInferenceGetVariable(model->params, "LAL_AMPORDER");
    else
        amporder = -1;

    REAL8 f_ref = 100.0;
    if (LALInferenceCheckVariable(model->params, "f_ref")) f_ref = *(REAL8 *)LALInferenceGetVariable(model->params, "f_ref");

    REAL8 fTemp = f_ref;

    if(LALInferenceCheckVariable(model->params,"chirpmass"))
    {
        mc  = *(REAL8*) LALInferenceGetVariable(model->params, "chirpmass");
        if (LALInferenceCheckVariable(model->params,"q")) {
            REAL8 q = *(REAL8 *)LALInferenceGetVariable(model->params,"q");
            q2masses(mc, q, &m1, &m2);
        } else {
            REAL8 eta = *(REAL8*) LALInferenceGetVariable(model->params, "eta");
            mc2masses(mc, eta, &m1, &m2);
        }
    }
    else if((m1_p=(REAL8 *)LALInferenceGetVariable(model->params, "mass1")) && (m2_p=(REAL8 *)LALInferenceGetVariable(model->params, "mass2")))
    {
        m1=*m1_p;
        m2=*m2_p;
    }
    else
    {
        fprintf(stderr,"No mass parameters found!");
        exit(0);
    }

    distance	= exp(LALInferenceGetREAL8Variable(model->params, "logdistance"))* LAL_PC_SI * 1.0e6;        /* distance (1 Mpc) in units of metres */

    phi0		= LALInferenceGetREAL8Variable(model->params, "phase"); /* START phase as per lalsimulation convention, radians*/

    /* Zenith angle between J and N in radians. Also known as inclination angle when spins are aligned */
    REAL8 thetaJN = acos(LALInferenceGetREAL8Variable(model->params, "costheta_jn"));     /* zenith angle between J and N in radians */

    /* Check if fLow is a model parameter, otherwise use data structure definition */
    if(LALInferenceCheckVariable(model->params, "flow"))
        f_low = *(REAL8*) LALInferenceGetVariable(model->params, "flow");
    else
        f_low = model->fLow;

    f_start = XLALSimInspiralfLow2fStart(f_low, amporder, approximant);
    f_max = 0.0; /* for freq domain waveforms this will stop at ISCO. Previously found using model->fHigh causes NaNs in waveform (see redmine issue #750)*/

    /* ==== SPINS ==== */
    /* We will default to spinless signal and then add in the spin components if required */
    /* If there are non-aligned spins, we must convert between the System Frame coordinates
     * and the cartestian coordinates */

    /* The cartesian spin coordinates (default 0), as passed to LALSimulation */
    REAL8 spin1x = 0.0;
    REAL8 spin1y = 0.0;
    REAL8 spin1z = 0.0;
    REAL8 spin2x = 0.0;
    REAL8 spin2y = 0.0;
    REAL8 spin2z = 0.0;

    /* System frame coordinates as used for jump proposals */
    REAL8 a_spin1 = 0.0;  /* Magnitude of spin1 */
    REAL8 a_spin2 = 0.0;  /* Magnitude of spin2 */
    REAL8 phiJL  = 0.0;  /* azimuthal angle of L_N on its cone about J radians */
    REAL8 tilt1   = 0.0;  /* zenith angle between S1 and LNhat in radians */
    REAL8 tilt2   = 0.0;  /* zenith angle between S2 and LNhat in radians */
    REAL8 phi12   = 0.0;  /* difference in azimuthal angle btwn S1, S2 in radians */

    /* Now check if we have spin amplitudes */
    if(LALInferenceCheckVariable(model->params, "a_spin1"))    a_spin1   = *(REAL8*) LALInferenceGetVariable(model->params, "a_spin1");
    if(LALInferenceCheckVariable(model->params, "a_spin2"))    a_spin2   = *(REAL8*) LALInferenceGetVariable(model->params, "a_spin2");

    /* Check if we have spin angles too */
    if(LALInferenceCheckVariable(model->params, "phi_jl"))
        phiJL = LALInferenceGetREAL8Variable(model->params, "phi_jl");
    if(LALInferenceCheckVariable(model->params, "tilt_spin1"))
        tilt1 = LALInferenceGetREAL8Variable(model->params, "tilt_spin1");
    if(LALInferenceCheckVariable(model->params, "tilt_spin2"))
        tilt2 = LALInferenceGetREAL8Variable(model->params, "tilt_spin2");
    if(LALInferenceCheckVariable(model->params, "phi12"))
        phi12 = LALInferenceGetREAL8Variable(model->params, "phi12");

    /* If we have tilt angles zero, then the spins are aligned and we just set the z component */
    /* However, if the waveform supports precession then we still need to get the right coordinate components */
    SpinSupport spin_support=XLALSimInspiralGetSpinSupportFromApproximant(approximant);
    if(tilt1==0.0 && tilt2==0.0 && (spin_support==LAL_SIM_INSPIRAL_SPINLESS || spin_support==LAL_SIM_INSPIRAL_ALIGNEDSPIN))
    {
        spin1z=a_spin1;
        spin2z=a_spin2;
        inclination = thetaJN; /* Inclination angle is just thetaJN */
    }
    else
    {   /* Template is not aligned-spin only. */
        /* Set all the other spin components according to the angles we received above */
        /* The transformation function doesn't know fLow, so f_ref==0 isn't interpretted as a request to use the starting frequency for reference. */
        if(fTemp==0.0)
            fTemp = f_start;

        XLAL_TRY(ret=XLALSimInspiralTransformPrecessingNewInitialConditions(
                                                                            &inclination, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z,
                                                                            thetaJN, phiJL, tilt1, tilt2, phi12, a_spin1, a_spin2, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, fTemp), errnum);
        if (ret == XLAL_FAILURE)
        {
            XLALPrintError(" ERROR in XLALSimInspiralTransformPrecessingNewInitialConditions(): error converting angles. errnum=%d: %s\n",errnum, XLALErrorString(errnum) );
            return;
        }
    }

    /* ==== TIDAL PARAMETERS ==== */
    REAL8 lambda1 = 0.;
    if(LALInferenceCheckVariable(model->params, "lambda1")) lambda1 = *(REAL8*) LALInferenceGetVariable(model->params, "lambda1");
    REAL8 lambda2 = 0.;
    if(LALInferenceCheckVariable(model->params, "lambda2")) lambda2 = *(REAL8*) LALInferenceGetVariable(model->params, "lambda2");
    REAL8 lambdaT = 0.;
    REAL8 dLambdaT = 0.;
    REAL8 sym_mass_ratio_eta = 0.;
    if(LALInferenceCheckVariable(model->params, "lambdaT")&&LALInferenceCheckVariable(model->params, "dLambdaT")){
        lambdaT = *(REAL8*) LALInferenceGetVariable(model->params, "lambdaT");
        dLambdaT = *(REAL8*) LALInferenceGetVariable(model->params, "dLambdaT");
        sym_mass_ratio_eta = m1*m2/((m1+m2)*(m1+m2));
        LALInferenceLambdaTsEta2Lambdas(lambdaT,dLambdaT,sym_mass_ratio_eta,&lambda1,&lambda2);
    }

    /* Only use GR templates */
    LALSimInspiralTestGRParam *nonGRparams = NULL;
    /* Fill in the extra parameters for testing GR, if necessary */
    for (UINT4 k=0; k<N_extra_params; k++)
    {
        if(LALInferenceCheckVariable(model->params,list_extra_parameters[k]))
        {
            XLALSimInspiralAddTestGRParam(&nonGRparams,list_extra_parameters[k],*(REAL8 *)LALInferenceGetVariable(model->params,list_extra_parameters[k]));
        }
    }
    /* Fill in PPE params if they are available */
    char PPEparam[64]="";
    const char *PPEnames[]= {"aPPE","alphaPPE","bPPE","betaPPE",NULL};
    for(UINT4 idx=0; PPEnames[idx]; idx++)
    {
        for(UINT4 ppeidx=0;; ppeidx++)
        {
            sprintf(PPEparam, "%s%d",PPEnames[idx],ppeidx);
            if(LALInferenceCheckVariable(model->params,PPEparam))
                XLALSimInspiralAddTestGRParam(&nonGRparams,PPEparam,LALInferenceGetREAL8Variable(model->params,PPEparam));
            else
                break;
        }
    }
    INT4 Nbands=-1; /* Use optimum number of bands */
    double mc_min=1.0/pow(2,0.2); /* For min 1.0-1.0 waveform */

    /* Vector of frequencies at which to compute FD template */
    static REAL8Sequence *frequencies = NULL;

    /* ==== Call the waveform generator ==== */
    if(model->domain == LAL_SIM_DOMAIN_FREQUENCY) {
        if(!frequencies) frequencies = LALInferenceMultibandFrequencies(Nbands,f_start,0.5/deltaT, model->deltaF, mc_min);


        XLAL_TRY(ret=XLALSimInspiralChooseFDWaveformFromCache(&hptilde, &hctilde, phi0,
                                                              0.0, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z,
                                                              spin2x, spin2y, spin2z, f_start, f_max, f_ref, distance, inclination,lambda1, lambda2, model->waveFlags, nonGRparams, amporder, order,
                                                              approximant,model->waveformCache, frequencies), errnum);

        /* if the waveform failed to generate, fill the buffer with zeros
         * so that the previous waveform is not left there
         */
        if(ret!=XLAL_SUCCESS){
            memset(model->freqhPlus->data->data,0,sizeof(model->freqhPlus->data->data[0])*model->freqhPlus->data->length);
            memset(model->freqhCross->data->data,0,sizeof(model->freqhCross->data->data[0])*model->freqhCross->data->length);
            if ( hptilde ) XLALDestroyCOMPLEX16FrequencySeries(hptilde);
            if ( hctilde ) XLALDestroyCOMPLEX16FrequencySeries(hctilde);
            XLALSimInspiralDestroyTestGRParam(nonGRparams);
            errnum&=~XLAL_EFUNC; /* Mask out the internal function failure bit */
            switch(errnum)
            {
                case XLAL_EDOM:
                    /* The waveform was called outside its domain. Return an empty vector but not an error */
                    XLAL_ERROR_VOID(XLAL_EUSR0);
                default:
                    /* Another error occurred that we can't handle. Propogate upward */
                    XLALSetErrno(errnum);
                    XLAL_ERROR_VOID(errnum,"%s: Template generation failed in XLALSimInspiralChooseFDWaveformFromCache:\n\
XLALSimInspiralChooseFDWaveformFromCache(&hptilde, &hctilde, \
%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, \
model->waveFlags(%d,%d,%d,%d,numreldata),nonGRparams,%d,%d,%d,model->waveformCache)\n",__func__,
		      phi0, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z,
		      f_start, f_max, f_ref, distance, inclination, lambda1, lambda2, (int) XLALSimInspiralGetSpinOrder(model->waveFlags),
		      (int) XLALSimInspiralGetTidalOrder(model->waveFlags),(int) XLALSimInspiralGetFrameAxis(model->waveFlags),
		      (int) XLALSimInspiralGetModesChoice(model->waveFlags),amporder, order, approximant);
            }
        }

        if (hptilde==NULL || hptilde->data==NULL || hptilde->data->data==NULL ) {
            XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'hptilde'.\n");
            XLAL_ERROR_VOID(XLAL_EFAULT);
        }
        if (hctilde==NULL || hctilde->data==NULL || hctilde->data->data==NULL ) {
            XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'hctilde'.\n");
            XLAL_ERROR_VOID(XLAL_EFAULT);
        }
        
        
        InterpolateWaveform(frequencies, hptilde, model->freqhPlus);
        InterpolateWaveform(frequencies, hctilde, model->freqhCross);

        /* Destroy the nonGr params */
        XLALSimInspiralDestroyTestGRParam(nonGRparams);

        REAL8 instant = model->freqhPlus->epoch.gpsSeconds + 1e-9*model->freqhPlus->epoch.gpsNanoSeconds;
        LALInferenceSetVariable(model->params, "time", &instant);


    } else {

        XLAL_TRY(ret=XLALSimInspiralChooseTDWaveformFromCache(&hplus, &hcross, phi0, deltaT,
                                                              m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z,
                                                              spin2x, spin2y, spin2z, f_start, f_ref, distance,
                                                              inclination, lambda1, lambda2, model->waveFlags, nonGRparams,
                                                              amporder, order, approximant,model->waveformCache), errnum);
        /* if the waveform failed to generate, fill the buffer with zeros
         * so that the previous waveform is not left there
         */
        if(ret!=XLAL_SUCCESS){
            memset(model->freqhPlus->data->data,0,sizeof(model->freqhPlus->data->data[0])*model->freqhPlus->data->length);
            memset(model->freqhCross->data->data,0,sizeof(model->freqhCross->data->data[0])*model->freqhCross->data->length);
            if ( hptilde ) XLALDestroyCOMPLEX16FrequencySeries(hptilde);
            if ( hctilde ) XLALDestroyCOMPLEX16FrequencySeries(hctilde);
            XLALSimInspiralDestroyTestGRParam(nonGRparams);
            errnum&=~XLAL_EFUNC; /* Mask out the internal function failure bit */
            switch(errnum)
            {
                case XLAL_EDOM:
                    /* The waveform was called outside its domain. Return an empty vector but not an error */
                    XLAL_ERROR_VOID(XLAL_EUSR0);
                default:
                    /* Another error occurred that we can't handle. Propogate upward */
                    XLALSetErrno(errnum);
                    XLAL_ERROR_VOID(errnum,"%s: Template generation failed in XLALSimInspiralChooseTDWaveformFromCache",__func__);
            }
        }

        /* The following complicated mess is a result of the following considerations:

         1) The discrete time samples of the template and the timeModel
         buffers will not, in general line up.

         2) The likelihood function will timeshift the template in the
         frequency domain to align it properly with the desired tc in
         each detector (these are different because the detectors
         receive the signal at different times).  Because this
         timeshifting is done in the frequency domain, the effective
         time-domain template is periodic.  We want to avoid the
         possibility of non-zero template samples wrapping around from
         the start/end of the buffer, since real templates are not
         periodic!

         3) If the template apporaches the ends of the timeModel buffer,
         then it should be tapered in the same way as the timeData
         (currently 0.4 seconds, hard-coded! Tukey window; see
         LALInferenceReadData.c, near line 233) so that template and
         signal in the data match.  However, as an optimization, we
         perform only one tapering and FFT-ing in the likelihood
         function; subsequent timeshifts for the different detectors
         will cause the tapered regions of the template and data to
         become mis-aligned.

         The algorthim we use is the following:

         1) Inject the template to align with the nearest sample in the
         timeModel buffer to the desired geocent_end time.

         2) Check whether either the start or the end of the template
         overlaps the tapered region, plus a safety buffer corresponding
         to a conservative estimate of the largest geocenter <-->
         detector timeshift.

         a) If there is no overlap at the start or end of the buffer,
         we're done.

         b) If there is an overlap, issue one warning per process
         (which can be disabled by setting the LAL debug level) about
         a too-short segment length, and return.
         */

        size_t waveLength = hplus->data->length;
        size_t bufLength = model->timehPlus->data->length;

        /* 2*Rearth/(c*deltaT)---2 is safety factor---is the maximum time
         shift for any earth-based detector. */
        size_t maxShift = (size_t)lround(4.255e-2/deltaT);

        /* Taper 0.4 seconds at start and end (hard-coded! in
         LALInferenceReadData.c, around line 233). */
        size_t taperLength = (size_t)lround(0.4/deltaT);

        /* Within unsafeLength of ends of buffer, possible danger of
         wrapping and/or tapering interactions. */
        size_t unsafeLength = taperLength + maxShift;

        REAL8 desiredTc = *(REAL8 *)LALInferenceGetVariable(model->params, "time");
        REAL8 tStart = XLALGPSGetREAL8(&(model->timehPlus->epoch));
        REAL8 tEnd = tStart + deltaT * model->timehPlus->data->length;

        if (desiredTc < tStart || desiredTc > tEnd) {
            XLALDestroyREAL8TimeSeries(hplus);
            XLALDestroyREAL8TimeSeries(hcross);

            XLAL_PRINT_ERROR("desired tc (%.4f) outside data buffer\n", desiredTc);
            XLAL_ERROR_VOID(XLAL_EDOM);
        }

        /* The nearest sample in model buffer to the desired tc. */
        size_t tcSample = (size_t)lround((desiredTc - XLALGPSGetREAL8(&(model->timehPlus->epoch)))/deltaT);

        /* The acutal coalescence time that corresponds to the buffer
         sample on which the waveform's tC lands. */
        REAL8 injTc = XLALGPSGetREAL8(&(model->timehPlus->epoch)) + tcSample*deltaT;

        /* The sample at which the waveform reaches tc. */
        size_t waveTcSample = (size_t)lround(-XLALGPSGetREAL8(&(hplus->epoch))/deltaT);

        /* 1 + (number of samples post-tc in waveform) */
        size_t wavePostTc = waveLength - waveTcSample;

        size_t bufStartIndex = (tcSample >= waveTcSample ? tcSample - waveTcSample : 0);
        size_t bufEndIndex = (wavePostTc + tcSample <= bufLength ? wavePostTc + tcSample : bufLength);
        size_t bufWaveLength = bufEndIndex - bufStartIndex;
        size_t waveStartIndex = (tcSample >= waveTcSample ? 0 : waveTcSample - tcSample);

        if (bufStartIndex < unsafeLength || (bufLength - bufEndIndex) <= unsafeLength) {
            /* The waveform could be timeshifted into a region where it will
             be tapered improperly, or even wrap around from the periodic
             timeshift.  Issue warning. */
            if (!sizeWarning) {
                fprintf(stderr, "WARNING: Generated template is too long to guarantee that it will not\n");
                fprintf(stderr, "WARNING:  (a) lie in a tapered region of the time-domain buffer\n");
                fprintf(stderr, "WARNING:  (b) wrap periodically when timeshifted in likelihood computation\n");
                fprintf(stderr, "WARNING: Either of these may cause differences between the template and the\n");
                fprintf(stderr, "WARNING: correct GW waveform in each detector.\n");
                fprintf(stderr, "WARNING: Parameter estimation will continue, but you should consider\n");
                fprintf(stderr, "WARNING: increasing the data segment length (using the --seglen) option.\n");
                sizeWarning = 1;
            }
        }

        /* Clear model buffers */
        memset(model->timehPlus->data->data, 0, sizeof(REAL8)*model->timehPlus->data->length);
        memset(model->timehCross->data->data, 0, sizeof(REAL8)*model->timehCross->data->length);

        /* Inject */
        memcpy(model->timehPlus->data->data + bufStartIndex,
               hplus->data->data + waveStartIndex,
               bufWaveLength*sizeof(REAL8));
        memcpy(model->timehCross->data->data + bufStartIndex,
               hcross->data->data + waveStartIndex,
               bufWaveLength*sizeof(REAL8));

        LALInferenceSetVariable(model->params, "time", &injTc);
    }

    if ( hplus ) XLALDestroyREAL8TimeSeries(hplus);
    if ( hcross ) XLALDestroyREAL8TimeSeries(hcross);
    if ( hptilde ) XLALDestroyCOMPLEX16FrequencySeries(hptilde);
    if ( hctilde ) XLALDestroyCOMPLEX16FrequencySeries(hctilde);

    return;
}

void LALInferenceTemplateXLALSimBurstChooseWaveform(LALInferenceModel *model)
/*************************************************************************************************************************/
/* Wrapper for LALSimulation waveforms:						                                                             */
/* XLALSimBurstChooseFDWaveform() and XLALSimBurstChooseTDWaveform().                                              */
/*                                                                                                                       */
/*  model->params parameters are:										                                         */
/*  - "name" description; type OPTIONAL (default value)										                             */
/*   "LAL_APPROXIMANT" burst approximant, BurstApproximant */
/*	"frequency" central frequency, REAL8                                                                          */
/*   "Q" quality, REAL8 (optional, depending on the WF)                                              */
/*   "duration" duration, REAL8 (optional, depending on the WF)                                      */
/*   "polar_angle" ellipticity polar angle, REAL8 */
/*   "polar_eccentricity" ellipticity ellipse eccentricity, REAL8 (optional)                                     */
/*                                                                                                                      */
/*************************************************************************************************************************/
{

  BurstApproximant approximant = (BurstApproximant) 0;
  unsigned long	i;
  static int sizeWarning = 0;
  int ret=0;
  INT4 errnum=0;
  REAL8 instant;


  REAL8TimeSeries *hplus=NULL;  /**< +-polarization waveform [returned] */
  REAL8TimeSeries *hcross=NULL; /**< x-polarization waveform [returned] */
  COMPLEX16FrequencySeries *hptilde=NULL, *hctilde=NULL;
  REAL8 deltaT,deltaF,
  freq=0.0,
  quality=0.0,
  duration=0.0,
  f_low, f_max,
  hrss=1.0,
  polar_ecc=1.0,polar_angle=LAL_PI/2.;
  LALSimBurstExtraParam *extraParams = NULL;

  if (LALInferenceCheckVariable(model->params, "LAL_APPROXIMANT"))
    approximant = *(BurstApproximant*) LALInferenceGetVariable(model->params, "LAL_APPROXIMANT");
  else {
    XLALPrintError(" ERROR in LALInferenceTemplateXLALSimBurstChooseWaveform(): (INT4) \"LAL_APPROXIMANT\" parameter not provided!\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }

  if(LALInferenceCheckVariable(model->params,"frequency"))
    {
      freq=*(REAL8*) LALInferenceGetVariable(model->params, "frequency");
  }

  if(LALInferenceCheckVariable(model->params,"quality"))
    {
      quality=*(REAL8*) LALInferenceGetVariable(model->params, "quality");
  }
  if(LALInferenceCheckVariable(model->params,"duration"))
    {
      duration=*(REAL8*) LALInferenceGetVariable(model->params, "duration");
  }
  /*
   * not needed since template is calculated at hrss=1 and scaled back in the likelihood
    if(LALInferenceCheckVariable(model->params,"hrss"))
      {
        hrss=*(REAL8*) LALInferenceGetVariable(model->params, "hrss");
    } else if (LALInferenceCheckVariable(model->params,"loghrss"))
        hrss=exp(*(REAL8*) LALInferenceGetVariable(model->params, "hrss"));

  if(LALInferenceCheckVariable(model->params,"alpha"))
    {
      alpha=*(REAL8*) LALInferenceGetVariable(model->params, "alpha");
      if (!extraParams) extraParams=XLALSimBurstCreateExtraParam("alpha",alpha);
      else XLALSimBurstAddExtraParam(&extraParams,"alpha",alpha);
      polar_angle=alpha;
  }
  */
  /* If someone wants to use old parametrization, allow for */
  if(LALInferenceCheckVariable(model->params,"polar_angle"))
    polar_angle=*(REAL8*) LALInferenceGetVariable(model->params, "polar_angle");
  if(LALInferenceCheckVariable(model->params,"polar_eccentricity"))
    polar_ecc=*(REAL8*) LALInferenceGetVariable(model->params, "polar_eccentricity");

  f_low=0.0; // These are not really used for burst signals.
  f_max = 0.0; /* for freq domain waveforms this will stop at Nyquist of lower, if the WF allows.*/


  if (model->timehCross==NULL) {
    XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'timeData'.\n");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
  deltaT = model->timehCross->deltaT;

  if(model->domain == LAL_SIM_DOMAIN_FREQUENCY) {
    if (model->freqhCross==NULL) {
      XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'freqhCross'.\n");
      XLAL_ERROR_VOID(XLAL_EFAULT);
    }

    deltaF = model->deltaF;

	/*Create BurstExtra params here and set depending on approx or let chooseFD do that*/


  XLAL_TRY(ret=XLALSimBurstChooseFDWaveformFromCache(&hptilde, &hctilde, deltaF,deltaT,freq,quality,duration,f_low,f_max,hrss,polar_angle,polar_ecc,extraParams,approximant,model->burstWaveformCache), errnum);
  //XLAL_TRY(ret=XLALSimBurstChooseFDWaveform(&hptilde, &hctilde, deltaF,deltaT,freq,quality,duration,f_low,f_max,hrss,polar_angle,polar_ecc,extraParams,approximant), errnum);
  if (ret == XLAL_FAILURE)
      {
        XLALPrintError(" ERROR in LALInferenceTemplateXLALSimBurstChooseWaveform(). errnum=%d: %s\n",errnum ,XLALErrorString(errnum));
        return;
      }
	if (hptilde==NULL || hptilde->data==NULL || hptilde->data->data==NULL ) {
	  XLALPrintError(" ERROR in LALInferenceTemplateXLALSimBurstChooseWaveform: encountered unallocated 'hptilde'.\n");
	  XLAL_ERROR_VOID(XLAL_EFAULT);
	}
	if (hctilde==NULL || hctilde->data==NULL || hctilde->data->data==NULL ) {
	  XLALPrintError(" ERROR in LALInferenceTemplateXLALSimBurstChooseWaveform: encountered unallocated 'hctilde'.\n");
	  XLAL_ERROR_VOID(XLAL_EFAULT);
	}

	COMPLEX16 *dataPtr = hptilde->data->data;

    for (i=0; i<model->freqhCross->data->length; ++i) {
      dataPtr = hptilde->data->data;
      if(i < hptilde->data->length){
        model->freqhPlus->data->data[i] = dataPtr[i];
      }else{
        model->freqhPlus->data->data[i] = 0.0;
      }
    }
    for (i=0; i<model->freqhCross->data->length; ++i) {
      dataPtr = hctilde->data->data;
      if(i < hctilde->data->length){
        model->freqhCross->data->data[i] = dataPtr[i];
      }else{
        model->freqhCross->data->data[i] = 0.0;
      }
    }

    /* Destroy the extra params */
    XLALSimBurstDestroyExtraParam(extraParams);

    instant= (model->timehCross->epoch.gpsSeconds + 1e-9*model->timehCross->epoch.gpsNanoSeconds);
    LALInferenceSetVariable(model->params, "time", &instant);
  }
 else {
    /*Time domain WF*/

   // XLAL_TRY(ret=XLALSimBurstChooseTDWaveform(&hplus, &hcross,deltaT,freq,quality,duration,f_low,f_max,hrss,polar_angle,polar_ecc,extraParams,approximant), errnum);
    XLAL_TRY(ret=XLALSimBurstChooseTDWaveformFromCache(&hplus, &hcross,deltaT,freq,quality,duration,f_low,f_max,hrss,polar_angle,polar_ecc,extraParams,approximant,model->burstWaveformCache), errnum);
    XLALSimBurstDestroyExtraParam(extraParams);
    if (ret == XLAL_FAILURE || hplus == NULL || hcross == NULL)
      {
        XLALPrintError(" ERROR in XLALSimBurstChooseTDWaveform(): error generating waveform. errnum=%d: %s\n",errnum ,XLALErrorString(errnum));
        for (i=0; i<model->timehCross->data->length; i++){
          model->timehPlus->data->data[i] = 0.0;
          model->timehCross->data->data[i] = 0.0;
        }
        return;
      }

    /* The following complicated mess is a result of the following considerations:

       1) The discrete time samples of the template and the timeModel
       buffers will not, in general line up.

       2) The likelihood function will timeshift the template in the
       frequency domain to align it properly with the desired tc in
       each detector (these are different because the detectors
       receive the signal at different times).  Because this
       timeshifting is done in the frequency domain, the effective
       time-domain template is periodic.  We want to avoid the
       possibility of non-zero template samples wrapping around from
       the start/end of the buffer, since real templates are not
       periodic!

       3) If the template apporaches the ends of the timeModel buffer,
       then it should be tapered in the same way as the timeData
       (currently 0.4 seconds, hard-coded! Tukey window; see
       LALInferenceReadData.c, near line 233) so that template and
       signal in the data match.  However, as an optimization, we
       perform only one tapering and FFT-ing in the likelihood
       function; subsequent timeshifts for the different detectors
       will cause the tapered regions of the template and data to
       become mis-aligned.

       The algorthim we use is the following:

       1) Inject the template to align with the nearest sample in the
       timeModel buffer to the desired geocent_end time.

       2) Check whether either the start or the end of the template
       overlaps the tapered region, plus a safety buffer corresponding
       to a conservative estimate of the largest geocenter <-->
       detector timeshift.

         a) If there is no overlap at the start or end of the buffer,
         we're done.

	 b) If there is an overlap, issue one warning per process
	 (which can be disabled by setting the LAL debug level) about
	 a too-short segment length, and return.
*/

    size_t waveLength = hplus->data->length;
    size_t bufLength = model->timehCross->data->length;

    /* 2*Rearth/(c*deltaT)---2 is safety factor---is the maximum time
       shift for any earth-based detector. */
    size_t maxShift = (size_t)lround(4.255e-2/hplus->deltaT);

    /* Taper 0.4 seconds at start and end (hard-coded! in
       LALInferenceReadData.c, around line 233). */
    REAL8 pad=model->padding;
    size_t taperLength = (size_t)lround(pad/hplus->deltaT);

    /* Within unsafeLength of ends of buffer, possible danger of
       wrapping and/or tapering interactions. */
    size_t unsafeLength = taperLength + maxShift;

    REAL8 desiredTc = *(REAL8 *)LALInferenceGetVariable(model->params, "time");
    REAL8 tStart = XLALGPSGetREAL8(&(model->timehCross->epoch));
    REAL8 tEnd = tStart + model->deltaT * model->timehCross->data->length;

    if (desiredTc < tStart || desiredTc > tEnd) {
      XLALDestroyREAL8TimeSeries(hplus);
      XLALDestroyREAL8TimeSeries(hcross);

      XLAL_PRINT_ERROR("desired tc (%.4f) outside data buffer\n", desiredTc);
      XLAL_ERROR_VOID(XLAL_EDOM);
    }

    /* The nearest sample in model buffer to the desired tc. */
    size_t tcSample = (size_t)lround((desiredTc - XLALGPSGetREAL8(&(model->timehPlus->epoch)))/model->deltaT);

    /* The acutal coalescence time that corresponds to the buffer
       sample on which the waveform's tC lands. */
    REAL8 injTc = XLALGPSGetREAL8(&(model->timehPlus->epoch)) + tcSample*model->deltaT;

    /* The sample at which the waveform reaches tc. */
    size_t waveTcSample = (size_t)lround(-XLALGPSGetREAL8(&(hplus->epoch))/hplus->deltaT);

    /* 1 + (number of samples post-tc in waveform) */
    size_t wavePostTc = waveLength - waveTcSample;

    size_t bufStartIndex = (tcSample >= waveTcSample ? tcSample - waveTcSample : 0);
    size_t bufEndIndex = (wavePostTc + tcSample <= bufLength ? wavePostTc + tcSample : bufLength);
    size_t bufWaveLength = bufEndIndex - bufStartIndex;
    size_t waveStartIndex = (tcSample >= waveTcSample ? 0 : waveTcSample - tcSample);

    if (bufStartIndex < unsafeLength || (bufLength - bufEndIndex) <= unsafeLength) {
      /* The waveform could be timeshifted into a region where it will
	 be tapered improperly, or even wrap around from the periodic
	 timeshift.  Issue warning. */
      if (!sizeWarning) {
	fprintf(stderr, "WARNING: Generated template is too long to guarantee that it will not\n");
	fprintf(stderr, "WARNING:  (a) lie in a tapered region of the time-domain buffer\n");
	fprintf(stderr, "WARNING:  (b) wrap periodically when timeshifted in likelihood computation\n");
	fprintf(stderr, "WARNING: Either of these may cause differences between the template and the\n");
	fprintf(stderr, "WARNING: correct GW waveform in each detector.\n");
	fprintf(stderr, "WARNING: Parameter estimation will continue, but you should consider\n");
	fprintf(stderr, "WARNING: increasing the data segment length (using the --seglen) option.\n");
	sizeWarning = 1;

      }
    }

    /* Clear IFOdata buffers */
    memset(model->timehPlus->data->data, 0, sizeof(REAL8)*model->timehPlus->data->length);
    memset(model->timehCross->data->data, 0, sizeof(REAL8)*model->timehCross->data->length);

    /* Inject */
    memcpy(model->timehPlus->data->data + bufStartIndex,
	   hplus->data->data + waveStartIndex,
	   bufWaveLength*sizeof(REAL8));
    memcpy(model->timehCross->data->data + bufStartIndex,
	   hcross->data->data + waveStartIndex,
	   bufWaveLength*sizeof(REAL8));

    LALInferenceSetVariable(model->params, "time", &injTc);
  }


  if ( hplus ) XLALDestroyREAL8TimeSeries(hplus);
  if ( hcross ) XLALDestroyREAL8TimeSeries(hcross);
  if ( hptilde ) XLALDestroyCOMPLEX16FrequencySeries(hptilde);
  if ( hctilde ) XLALDestroyCOMPLEX16FrequencySeries(hctilde);

  return;
}


void LALInferenceTemplateXLALSimBurstSineGaussianF(LALInferenceModel *model)
/*************************************************************************************************************************/
/* Wrapper for LALSimulation waveforms:						                                                             */
/* XLALInferenceBurstChooseFDWaveform() and XLALInferenceBurstChooseTDWaveform().                                              */
/*                                                                                                                       */
/*  model->params parameters are:										                                         */
/*  - "name" description; type OPTIONAL (default value)										                             */
/*   "LAL_APPROXIMANT" burst approximant, BurstApproximant */
/*	"frequency" central frequency, REAL8                                                                          */
/*   "Q" quality, REAL8 (optional, depending on the WF)                                              */
/*   "duration" duration, REAL8 (optional, depending on the WF)                                      */
/*   "polar_angle" ellipticity polar angle, REAL8 */
/*   "polar_eccentricity" ellipticity ellipse eccentricity, REAL8 (optional)                                     */
/*                                                                                                                      */
/*************************************************************************************************************************/
{

  unsigned long	i;
  int ret=0;
  INT4 errnum=0;
  REAL8 instant;


  REAL8TimeSeries *hplus=NULL;  /**< +-polarization waveform [returned] */
  REAL8TimeSeries *hcross=NULL; /**< x-polarization waveform [returned] */
  COMPLEX16FrequencySeries *hptilde=NULL, *hctilde=NULL;
  REAL8 deltaT,deltaF,
  freq=0.0,
  quality=0.0,tau=0.0,
  hrss=1.0;

  freq=*(REAL8*) LALInferenceGetVariable(model->params, "frequency");
  quality=*(REAL8*) LALInferenceGetVariable(model->params, "quality");
  if(LALInferenceCheckVariable(model->params,"duration"))
  {
    tau=*(REAL8*) LALInferenceGetVariable(model->params, "duration");
    quality=tau*freq*LAL_SQRT2*LAL_PI;
  }
  REAL8 polar_angle=*(REAL8*) LALInferenceGetVariable(model->params, "polar_angle");
  /* If someone wants to use old parametrization, allow for */
  REAL8 polar_ecc=*(REAL8*) LALInferenceGetVariable(model->params, "polar_eccentricity");

  if (model->timehCross==NULL) {
    XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'timeData'.\n");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
  deltaT = model->timehCross->deltaT;

  if (model->freqhCross==NULL) {
      XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'freqhCross'.\n");
      XLAL_ERROR_VOID(XLAL_EFAULT);
  }

  deltaF = model->deltaF;
  XLAL_TRY(ret=XLALInferenceBurstSineGaussianFFast(&hptilde, &hctilde, quality,freq,hrss, polar_ecc, polar_angle,deltaF,deltaT), errnum);
  if (ret == XLAL_FAILURE)
      {
        XLALPrintError(" ERROR in LALInferenceTemplateXLALSimBurstChooseWaveform(). errnum=%d: %s\n",errnum,XLALErrorString(errnum) );
        return;
      }
	if (hptilde==NULL || hptilde->data==NULL || hptilde->data->data==NULL ) {
	  XLALPrintError(" ERROR in LALInferenceTemplateXLALInferenceBurstChooseWaveform: encountered unallocated 'hptilde'.\n");
	  XLAL_ERROR_VOID(XLAL_EFAULT);
	}
	if (hctilde==NULL || hctilde->data==NULL || hctilde->data->data==NULL ) {
	  XLALPrintError(" ERROR in LALInferenceTemplateXLALInferenceBurstChooseWaveform: encountered unallocated 'hctilde'.\n");
	  XLAL_ERROR_VOID(XLAL_EFAULT);
	}
  size_t lower =(size_t) ( hctilde->f0/hctilde->deltaF);
  size_t upper= (size_t) ( hctilde->data->length + lower);
  COMPLEX16 *dataPtrP = hptilde->data->data;
  COMPLEX16 *dataPtrC = hctilde->data->data;
  for (i=0; i<lower; ++i) {
    model->freqhPlus->data->data[i] = 0.0;
    model->freqhCross->data->data[i] = 0.0;
  }

  if (upper>model->freqhPlus->data->length)
    upper=model->freqhPlus->data->length;
  else{
    for (i=upper; i<model->freqhPlus->data->length; ++i) {
      model->freqhPlus->data->data[i] = 0.0;
      model->freqhCross->data->data[i] = 0.0;
    }
  }
  for (i=lower; i<upper; ++i) {
    model->freqhPlus->data->data[i] = dataPtrP[i-lower];
    model->freqhCross->data->data[i] = dataPtrC[i-lower];
  }

  instant= (model->timehCross->epoch.gpsSeconds + 1e-9*model->timehCross->epoch.gpsNanoSeconds);
  LALInferenceSetVariable(model->params, "time", &instant);


  if ( hplus ) XLALDestroyREAL8TimeSeries(hplus);
  if ( hcross ) XLALDestroyREAL8TimeSeries(hcross);
  if ( hptilde ) XLALDestroyCOMPLEX16FrequencySeries(hptilde);
  if ( hctilde ) XLALDestroyCOMPLEX16FrequencySeries(hctilde);

  return;
}


void LALInferenceDumptemplateFreqDomain(LALInferenceVariables *currentParams,
                                        LALInferenceModel *model,
                                        const char *filename)
/* de-bugging function writing (frequency-domain) template to a CSV file */
/* File contains real & imaginary parts of plus & cross components.      */
/* Template amplitude is scaled to 1Mpc distance.                        */
{
  FILE *outfile=NULL;
  REAL8 f;
  UINT4 i;

  REAL8 deltaF = model->deltaF;

  LALInferenceCopyVariables(currentParams, model->params);
  model->templt(model);
  if (model->domain == LAL_SIM_DOMAIN_TIME)
    LALInferenceExecuteFT(model);

  outfile = fopen(filename, "w");
  /*fprintf(outfile, "f PSD dataRe dataIm signalPlusRe signalPlusIm signalCrossRe signalCrossIm\n");*/
  fprintf(outfile, "\"f\",\"signalPlusRe\",\"signalPlusIm\",\"signalCrossRe\",\"signalCrossIm\"\n");
  for (i=0; i<model->freqhPlus->data->length; ++i){
    f = ((double) i) * deltaF;
    fprintf(outfile, "%f,%e,%e,%e,%e\n",
            f,
            creal(model->freqhPlus->data->data[i]),
            cimag(model->freqhPlus->data->data[i]),
            creal(model->freqhCross->data->data[i]),
            cimag(model->freqhCross->data->data[i]));
  }
  fclose(outfile);
  fprintf(stdout, " wrote (frequency-domain) template to CSV file \"%s\".\n", filename);
}


void LALInferenceDumptemplateTimeDomain(LALInferenceVariables *currentParams,
                                        LALInferenceModel *model,
					                    const char *filename)
/* de-bugging function writing (frequency-domain) template to a CSV file */
/* File contains real & imaginary parts of plus & cross components.      */
/* Template amplitude is scaled to 1Mpc distance.                        */
{
  FILE *outfile=NULL;
  REAL8 deltaT, t, epoch; // deltaF - set but not used
  UINT4 i;

  LALInferenceCopyVariables(currentParams, model->params);
  model->templt(model);
  if (model->domain == LAL_SIM_DOMAIN_FREQUENCY)
    LALInferenceExecuteInvFT(model);

  outfile = fopen(filename, "w");
  fprintf(outfile, "\"t\",\"signalPlus\",\"signalCross\"\n");
  deltaT = model->deltaT;
  epoch = XLALGPSGetREAL8(&model->timehPlus->epoch);
  for (i=0; i<model->timehPlus->data->length; ++i){
    t =  epoch + ((double) i) * deltaT;
    fprintf(outfile, "%f,%e,%e\n",
            t,
            model->timehPlus->data->data[i],
            model->timehCross->data->data[i]);
  }
  fclose(outfile);
  fprintf(stdout, " wrote (time-domain) template to CSV file \"%s\".\n", filename);
}
