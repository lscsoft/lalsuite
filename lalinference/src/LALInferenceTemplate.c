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
#include <lal/Sequence.h>
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
#include <lal/LALSimRingdown.h>
#include <lal/LALInferenceTemplate.h>
#include <lal/LALSimBurstWaveformCache.h>
#include <lal/LALSimBurst.h>
#include <lal/NRWaveInject.h>
#include <gsl/gsl_interp.h>
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

static void q2eta(double q, double *eta);
static void q2masses(double mc, double q, double *m1, double *m2);


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

static void q2eta(double q, double *eta)
/* Compute symmetric mass ratio (eta) for a given  */
/* asymmetric mass ratio (q).                       */
/* (note: q = m2/m1, where m1 >= m2)               */
{
  *eta = q/pow(1+q,2.0);
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

void LALInferenceTemplateROQ(LALInferenceModel *model)
{
double m1,m2,mc,eta,q,iota=0;
/* Prefer m1 and m2 if available (i.e. for injection template) */
  if(LALInferenceCheckVariable(model->params,"mass1")&&LALInferenceCheckVariable(model->params,"mass2"))
  {
    m1=*(REAL8 *)LALInferenceGetVariable(model->params,"mass1");
    m2=*(REAL8 *)LALInferenceGetVariable(model->params,"mass2");
    eta=m1*m2/((m1+m2)*(m1+m2));
    mc=pow(eta , 0.6)*(m1+m2);
  }
  else
  {
    if (LALInferenceCheckVariable(model->params,"q")) {
      q = *(REAL8 *)LALInferenceGetVariable(model->params,"q");
      q2eta(q, &eta);
    }
    else
      eta = *(REAL8*) LALInferenceGetVariable(model->params, "eta");
    mc       = *(REAL8*) LALInferenceGetVariable(model->params, "chirpmass");
    mc2masses(mc, eta, &m1, &m2);
  }

  iota = *(REAL8*) LALInferenceGetVariable(model->params, "theta_jn");     /* zenith angle between J and N in radians */

  double cosiota = cos(iota);
  double plusCoef  = 0.5 * (1.0 + cosiota*cosiota);
  double crossCoef = cosiota;
  /* external: SI; internal: solar masses */
  const REAL8 m = m1 + m2;
  const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
  const REAL8 etap2 = eta * eta;
  const REAL8 etap3 = etap2 * eta;
  const REAL8 piM = LAL_PI * m_sec;
  const REAL8 mSevenBySix = -7./6.;
  const REAL8 phic = *(REAL8 *)LALInferenceGetVariable(model->params,"phase");
  const REAL8 r = 1e6*LAL_PC_SI;
  REAL8 logv0 = log(1.); //the standard tf2 definition is log(v0), but I've changed it to reflect Scott's convention
  REAL8 shft, amp0;//, f_max;
  REAL8 psiNewt, psi2, psi3, psi4, psi5, psi6, psi6L, psi7, psi3S, psi4S, psi5S;
  REAL8 eta_fac = -113. + 76. * eta;
  REAL8 chi=0; //NOTE: chi isn't used here yet, so we just set it to zero
  gsl_complex h_i;

  /* extrinsic parameters */
  amp0 = -pow(m_sec, 5./6.) * sqrt(5.*eta / 24.) / (Pi_p2by3 * r / LAL_C_SI);
  shft = 0;//LAL_TWOPI * (tStart.gpsSeconds + 1e-9 * tStart.gpsNanoSeconds);

  /* spin terms in the amplitude and phase (in terms of the reduced
   * spin parameter */
  psi3S = 113.*chi/3.;
  psi4S = 63845.*(-81. + 4.*eta)*chi*chi/(8. * eta_fac * eta_fac);
  psi5S = -565.*(-146597. + 135856.*eta + 17136.*etap2)*chi/(2268.*eta_fac);

  /* coefficients of the phase at PN orders from 0 to 3.5PN */
  psiNewt = 3./(128.*eta);
  psi2 = 3715./756. + 55.*eta/9.;
  psi3 = psi3S - 16.*LAL_PI;
  psi4 = 15293365./508032. + 27145.*eta/504. + 3085.*eta*eta/72. + psi4S;
  psi5 = (38645.*LAL_PI/756. - 65.*LAL_PI*eta/9. + psi5S);
  psi6 = 11583231236531./4694215680. - (640.*Pi_p2)/3. - (6848.*LAL_GAMMA)/21.
           + (-5162.983708047263 + 2255.*Pi_p2/12.)*eta
           + (76055.*etap2)/1728. - (127825.*etap3)/1296.;
  psi6L = -6848./21.;
  psi7 = (77096675.*LAL_PI)/254016. + (378515.*LAL_PI*eta)/1512.
           - (74045.*LAL_PI*eta*eta)/756.;

  for (unsigned int i = 0; i < model->roq->frequencyNodes->size; i++) {
    /* fourier frequency corresponding to this bin */
    const REAL8 f = gsl_vector_get(model->roq->frequencyNodes, i);
    const REAL8 v3 = piM*f;

    /* PN expansion parameter */
    REAL8 v, v2, v4, v5, v6, v7, logv, Psi, amp;
    v = cbrt(v3);
    v2 = v*v; v4 = v3*v; v5 = v4*v; v6 = v3*v3; v7 = v6*v;
    logv = log(v);

    /* compute the phase and amplitude */
    Psi = psiNewt / v5 * (1.
         + psi2 * v2 + psi3 * v3 + psi4 * v4
         + psi5 * v5 * (1. + 3. * (logv - logv0))
         + (psi6 + psi6L * (log4 + logv)) * v6 + psi7 * v7);

    amp = amp0 * pow(f, mSevenBySix);

    GSL_SET_COMPLEX(&h_i, amp * cos(Psi + shft * f - 2.*phic - LAL_PI_4), - amp * sin(Psi + shft * f - 2.*phic - LAL_PI_4));

    gsl_vector_complex_set(model->roq->hplus, i, gsl_complex_mul_real(h_i,plusCoef));
    gsl_vector_complex_set(model->roq->hcross, i, gsl_complex_mul_real(gsl_complex_mul_imag(h_i,-1.0),crossCoef));

  }
	return;
}
void LALInferenceTemplateROQ_amp_squared(LALInferenceModel *model)
{

double m1,m2,mc,eta,q;
/* Prefer m1 and m2 if available (i.e. for injection template) */
 if(LALInferenceCheckVariable(model->params,"mass1")&&LALInferenceCheckVariable(model->params,"mass2"))
    {
      m1=*(REAL8 *)LALInferenceGetVariable(model->params,"mass1");
      m2=*(REAL8 *)LALInferenceGetVariable(model->params,"mass2");
      eta=m1*m2/((m1+m2)*(m1+m2));
      mc=pow(eta , 0.6)*(m1+m2);
    }
  else
    {
      if (LALInferenceCheckVariable(model->params,"q")) {
        q = *(REAL8 *)LALInferenceGetVariable(model->params,"q");
        q2eta(q, &eta);
      }
      else
        eta = *(REAL8*) LALInferenceGetVariable(model->params, "eta");
        mc       = *(REAL8*) LALInferenceGetVariable(model->params, "chirpmass");
      mc2masses(mc, eta, &m1, &m2);
    } 
    /* external: SI; internal: solar masses */
    const REAL8 m = m1 + m2;
    const REAL8 m_sec = m * LAL_MTSUN_SI;  /* total mass in seconds */
    const REAL8 r = 1e6*LAL_PC_SI;
    double amp_squared;

    amp_squared = pow( pow(m_sec, 5./6.) * sqrt(5.*eta / 24.) / (Pi_p2by3 * r / LAL_C_SI), 2. );
  
    *(model->roq->amp_squared) = amp_squared;

}


REAL8 fLow2fStart(REAL8 fLow, INT4 ampOrder, INT4 approximant)
/*  Compute the minimum frequency for waveform generation */
/*  using amplitude orders above Newtonian.  The waveform */
/*  generator turns on all orders at the orbital          */
/*  associated with fMin, so information from higher      */
/*  orders is not included at fLow unless fMin is         */
/*  sufficiently low.                                     */
{
  if (ampOrder == -1) {
      if (approximant == SpinTaylorT2 || approximant == SpinTaylorT4)
          ampOrder = MAX_PRECESSING_AMP_PN_ORDER;
      else
          ampOrder = MAX_NONPRECESSING_AMP_PN_ORDER;
  }

    REAL8 fStart;
    fStart = fLow * 2./(ampOrder+2);
    return fStart;
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
    static int sizeWarning = 0;

    REAL8TimeSeries *hplus=NULL;  /**< +-polarization waveform */
    REAL8TimeSeries *hcross=NULL; /**< x-polarization waveform */
    //    REAL8TimeSeries       *signalvecREAL8=NULL;
    REAL8 Q, centre_frequency,hrss,eccentricity,polar_angle;
    UINT8 windowshift=(UINT8) ceil(model->padding/model->deltaT);
    UINT4 i=0;

    Q = *(REAL8*) LALInferenceGetVariable(model->params, "quality");
    centre_frequency=  *(REAL8*) LALInferenceGetVariable(model->params, "frequency");  
    /*if(LALInferenceCheckVariable(model->params,"hrss"))
    hrss=*(REAL8*) LALInferenceGetVariable(model->params, "hrss"); 
    else if(LALInferenceCheckVariable(model->params,"loghrss"))
    {
    loghrss=*(REAL8*) LALInferenceGetVariable(model->params, "loghrss"); 
    hrss=exp(loghrss);
    }
    else {fprintf(stderr,"ERROR (In LALInferenceTemplate): modelParams does not contain hrss or loghrss. Exiting...\n"); exit(1);}
    */

    /*Always calculate the template at fixed hrss of 1. That will avoid recalculation of the template unless Q, f, polar_angle, polar_eccentricity are varied */
    hrss=1.0;

    polar_angle=*(REAL8*) LALInferenceGetVariable(model->params, "polar_angle"); 
    eccentricity=*(REAL8*) LALInferenceGetVariable(model->params, "eccentricity"); 

    XLALSimBurstDampedSinusoid(&hplus,&hcross, Q, centre_frequency,hrss,eccentricity,polar_angle,model->deltaT);
 
    REAL8 instant = model->freqhPlus->epoch.gpsSeconds + 1e-9*model->freqhPlus->epoch.gpsNanoSeconds;
    LALInferenceSetVariable(model->params, "time", &instant);

    /* write template (time axis) location in "->modelParams" so that     */
    /* template corresponds to stored parameter values                    */
    /* and other functions may time-shift template to where they want it: */
    
    instant=instant+(INT8)windowshift*model->deltaT; //leave enough room for the tuckey windowing of the data.
    LALInferenceSetVariable(model->params, "time", &instant);
    
    if(hplus->data && hcross->data){
      if(hplus->data->length+2*windowshift<=model->timehPlus->data->length){ //check whether the model->timehCross->data vector is long enough to store the waveform produced
        for (i=0; i<model->timehCross->data->length; i++){
          if(i>=((unsigned long int)(hplus->data->length) + windowshift)  || i<windowshift || isnan(hplus->data->data[i-(INT8)windowshift]) || isnan(hcross->data->data[i-(INT8)windowshift])){
            model->timehPlus->data->data[i] = 0;
            model->timehCross->data->data[i] = 0;
          }else{
            model->timehPlus->data->data[i] = hplus->data->data[i-(INT8)windowshift];
            model->timehCross->data->data[i] = hcross->data->data[i-(INT8)windowshift];
          }
        }
      }else{
        if (!sizeWarning) {
          sizeWarning = 1;
          fprintf(stderr, "WARNING: hplus->data->length = %d is longer than model->timehPlus->data->length = %d minus windowshift = %d.\n", hplus->data->length, model->timehCross->data->length, (int) windowshift);
          if(hplus->data->length + (int) windowshift > model->timehCross->data->length)
            fprintf(stderr, "The waveform template used will be missing its first %d points. Consider increasing the segment length (--seglen). (in %s, line %d)\n",hplus->data->length - model->timehCross->data->length + (int) windowshift , __FILE__, __LINE__);
          else
            fprintf(stderr, "The waveform template used will have its first %d points tapered. Consider increasing the segment length (--seglen). (in %s, line %d)\n",hplus->data->length - model->timehCross->data->length + 2*(int)windowshift , __FILE__, __LINE__);
        }
        for (i=0; i<model->timehCross->data->length; i++){
          if((INT8)i>=(INT8)model->timehCross->data->length-(INT8)windowshift || (INT8)i+(INT8)hplus->data->length-(INT8)model->timehCross->data->length+(INT8)windowshift < 0 || isnan(hplus->data->data[(INT8)i+(INT8)hplus->data->length-(INT8)model->timehCross->data->length+(INT8)windowshift]) || isnan(hcross->data->data[(INT8)i+(INT8)hcross->data->length-(INT8)model->timehCross->data->length+(INT8)windowshift]) ){
            model->timehPlus->data->data[i] = 0.0;
            model->timehCross->data->data[i] = 0.0;
          }else{
      //        printf("Im here i=%d",(INT4) i);         

            model->timehPlus->data->data[i] = hplus->data->data[(INT8)i+(INT8)hplus->data->length-(INT8)model->timehCross->data->length+(INT8)windowshift];
            model->timehCross->data->data[i] = hcross->data->data[(INT8)i+(INT8)hcross->data->length-(INT8)model->timehCross->data->length+(INT8)windowshift];
          }
        }
        instant-= ((INT8)hplus->data->length-(INT8)model->timehCross->data->length+2*(INT8)windowshift)*model->timehCross->deltaT;
        LALInferenceSetVariable(model->params, "time", &instant);
      }
    }else{
      for (i=0; i<model->timehCross->data->length; i++){
        model->timehPlus->data->data[i] = 0;
        model->timehCross->data->data[i] = 0;
      }
      fprintf( stderr, " ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): no generated waveform.\n");
    }

REAL8 max=0.0;
 for (i=0; i<model->timehCross->data->length; i++){
        max=max>model->timehPlus->data->data[i] ?max:model->timehPlus->data->data[i];
        max=max>model->timehCross->data->data[i] ?max:model->timehCross->data->data[i];
      }
      if (max==0.0)
      printf("zeromax\n");

                if ( hplus ) XLALDestroyREAL8TimeSeries(hplus);
                if ( hcross ) XLALDestroyREAL8TimeSeries(hcross);
 
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
/*   - "fRef"               frequency at which the (frequency dependent) parameters are defined; REAL8 OPTIONAL (0.0)    */
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
/*   - "distance"           distance in Mpc                                                                              */
/*   - "theta_jn");      zenith angle between J and N in radians;            REAL8                                    */
/*   - "phi_jl");        azimuthal angle of L_N on its cone about J radians; REAL8                                    */
/*   - "tilt_spin1");    zenith angle between S1 and LNhat in radians;       REAL8                                    */
/*   - "tilt_spin2");    zenith angle between S2 and LNhat in radians;       REAL8                                    */
/*   - "phi12");         difference in azimuthal angle between S1, S2 in radians;   REAL8                             */
/*   - "a_spin1"            magnitude of spin 1 in general configuration, 0<a_spin1<1; REAL8 OPTIONAL (0.0)              */
/*   - "a_spin2"            magnitude of spin 2 in general configuration, 0<a_spin1<1; REAL8 OPTIONAL (0.0)              */
/*   - "spin1"              magnitude of spin 1 in aligned configuration, -1<spin1<1;  REAL8 OPTIONAL (0.0)              */
/*   - "spin2"              magnitude of spin 2 in aligned configuration, -1<spin1<1;  REAL8 OPTIONAL (0.0)              */
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

  unsigned long	i;
  static int sizeWarning = 0;
  int ret=0;
  INT4 errnum=0;
  
  REAL8TimeSeries *hplus=NULL;  /**< +-polarization waveform [returned] */
  REAL8TimeSeries *hcross=NULL; /**< x-polarization waveform [returned] */
  COMPLEX16FrequencySeries *hptilde=NULL, *hctilde=NULL;
  REAL8 mc;
  REAL8 phi0, deltaT, m1, m2, spin1x=0.0, spin1y=0.0, spin1z=0.0, spin2x=0.0, spin2y=0.0, spin2z=0.0, f_low, f_start, distance, inclination;
  
  REAL8 *m1_p,*m2_p;
  REAL8 deltaF, f_max;
  
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
    
  REAL8 fRef = 100.0;
  if (LALInferenceCheckVariable(model->params, "fRef")) fRef = *(REAL8 *)LALInferenceGetVariable(model->params, "fRef");

  REAL8 fTemp = fRef;

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

  distance	= LALInferenceGetREAL8Variable(model->params,"distance")* LAL_PC_SI * 1.0e6;        /* distance (1 Mpc) in units of metres */

  /* Default to spinless signals if spin amplitude are not present in model->params */
  REAL8 a_spin1		= 0.0;
  if(LALInferenceCheckVariable(model->params, "a_spin1"))    a_spin1   = *(REAL8*) LALInferenceGetVariable(model->params, "a_spin1");
  REAL8 a_spin2    = 0.0;
  if(LALInferenceCheckVariable(model->params, "a_spin2"))    a_spin2   = *(REAL8*) LALInferenceGetVariable(model->params, "a_spin2");

  phi0		= *(REAL8*) LALInferenceGetVariable(model->params, "phase"); /* START phase as per lalsimulation convention*/

  /* Check if fLow is a model parameter, otherwise use data structure definition */
  if(LALInferenceCheckVariable(model->params, "flow"))
    f_low = *(REAL8*) LALInferenceGetVariable(model->params, "flow");
  else
    f_low = model->fLow /** 0.9 */;

  f_start = fLow2fStart(f_low, amporder, approximant);
  f_max = 0.0; /* for freq domain waveforms this will stop at ISCO. Previously found using model->fHigh causes NaNs in waveform (see redmine issue #750)*/

  int aligned_spins=0;
  /* We first check if we are deadling with a spin-aligned only template, for which we use "spin1" and "spin2" names */
  if(LALInferenceCheckVariable(model->params, "spin1")){
    spin1z= *(REAL8*) LALInferenceGetVariable(model->params, "spin1");
    spin1x=0.0;
    spin1y=0.0;
    aligned_spins+=1;
  }
  if(LALInferenceCheckVariable(model->params, "spin2")){
    spin2z= *(REAL8*) LALInferenceGetVariable(model->params, "spin2");
    spin2x=0.0;
    spin2y=0.0;
    aligned_spins+=1;
  }

  /* Set inclination to something sensible now, because for spin aligned we won't enter in the blocks below, where theta_JN 
   *  is set for no-spin or precessing spins.  */
  inclination       = *(REAL8*) LALInferenceGetVariable(model->params, "theta_jn");           /* inclination in radian */
  
  if (aligned_spins==0){
    /* Template is not spin-aligned only.
    * Set the all other spins variables (that is an overkill. We can just check if there are spins in the first place, and if there aren't don't bother calculating them (they'll be zero) ) */
    REAL8 phiJL=0.0;
    REAL8 tilt1=0.0;
    REAL8 tilt2=0.0;
    REAL8 phi12=0.0;

    /* IMPORTANT NOTE: default to spin aligned case (i.e. tilt1=tilt2=0) if no angles are provided for the spins.
    * If you change this, must also change LALInferenceInitCBC.c
    */
    REAL8 thetaJN = *(REAL8*) LALInferenceGetVariable(model->params, "theta_jn");     /* zenith angle between J and N in radians */
    if(LALInferenceCheckVariable(model->params, "phi_jl"))
      phiJL = *(REAL8*) LALInferenceGetVariable(model->params, "phi_jl");     /* azimuthal angle of L_N on its cone about J radians */
    if(LALInferenceCheckVariable(model->params, "tilt_spin1"))
      tilt1 = *(REAL8*) LALInferenceGetVariable(model->params, "tilt_spin1");     /* zenith angle between S1 and LNhat in radians */
    if(LALInferenceCheckVariable(model->params, "tilt_spin2"))
      tilt2 = *(REAL8*) LALInferenceGetVariable(model->params, "tilt_spin2");     /* zenith angle between S2 and LNhat in radians */
    if(LALInferenceCheckVariable(model->params, "phi12"))
      phi12 = *(REAL8*) LALInferenceGetVariable(model->params, "phi12");      /* difference in azimuthal angle btwn S1, S2 in radians */

    /* The transformation function doesn't know fLow, so fRef==0 isn't interpretted as a request to use the starting frequency for reference. */
    if(fTemp==0.0)
      fTemp = f_start;

    XLAL_TRY(ret=XLALSimInspiralTransformPrecessingInitialConditions(
    &inclination, &spin1x, &spin1y, &spin1z, &spin2x, &spin2y, &spin2z,
    thetaJN, phiJL, tilt1, tilt2, phi12, a_spin1, a_spin2, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, fTemp), errnum);
    if (ret == XLAL_FAILURE)
    {
      XLALPrintError(" ERROR in XLALSimInspiralTransformPrecessingInitialConditions(): error converting angles. errnum=%d\n",errnum );
      return;
    }

  }
	
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

  LALSimInspiralTestGRParam *nonGRparams = NULL;
  
  deltaT = model->deltaT;
  
  if(model->domain == LAL_SIM_DOMAIN_FREQUENCY) {
    deltaF = model->deltaF;
    
	XLAL_TRY(ret=XLALSimInspiralChooseFDWaveformFromCache(&hptilde, &hctilde, phi0,
            deltaF, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z,
            spin2x, spin2y, spin2z, f_start, f_max, fRef, distance, inclination,lambda1, lambda2, model->waveFlags, nonGRparams, amporder, order,
            approximant,model->waveformCache), errnum);

	if (hptilde==NULL || hptilde->data==NULL || hptilde->data->data==NULL ) {
	  XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'hptilde'.\n");
	  XLAL_ERROR_VOID(XLAL_EFAULT);
	}
	if (hctilde==NULL || hctilde->data==NULL || hctilde->data->data==NULL ) {
	  XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'hctilde'.\n");
	  XLAL_ERROR_VOID(XLAL_EFAULT);
	}
      
	COMPLEX16 *dataPtr = hptilde->data->data;

    for (i=0; i<model->freqhPlus->data->length; ++i) {
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
    
    
    /* Destroy the nonGr params */
    XLALSimInspiralDestroyTestGRParam(nonGRparams);
    
    REAL8 instant = model->freqhPlus->epoch.gpsSeconds + 1e-9*model->freqhPlus->epoch.gpsNanoSeconds;
    LALInferenceSetVariable(model->params, "time", &instant);
    
  } else {

    XLAL_TRY(ret=XLALSimInspiralChooseTDWaveformFromCache(&hplus, &hcross, phi0, deltaT,
            m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z,
            spin2x, spin2y, spin2z, f_start, fRef, distance,
            inclination, lambda1, lambda2, model->waveFlags, nonGRparams,
            amporder, order, approximant,model->waveformCache), errnum);
    XLALSimInspiralDestroyTestGRParam(nonGRparams);
    if (ret == XLAL_FAILURE || hplus == NULL || hcross == NULL)
      {
	XLALPrintError(" ERROR in XLALSimInspiralChooseWaveformFromCache(): error generating waveform. errnum=%d\n",errnum );
	for (i=0; i<model->timehPlus->data->length; i++){
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
    size_t bufLength = model->timehPlus->data->length;

    /* 2*Rearth/(c*deltaT)---2 is safety factor---is the maximum time
       shift for any earth-based detector. */
    size_t maxShift = (size_t)lround(4.255e-2/deltaT); 

    /* Taper pad seconds at start and end */
    REAL8 pad=model->padding;
    size_t taperLength = (size_t)lround(pad/deltaT);

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
    REAL8 injTc = XLALGPSGetREAL8(&(model->timehPlus->epoch)) +(tcSample)*deltaT;
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
/*   "alpha"  ellipticity, REAL8 (optional, depending on the WF)                                     */
/*   "phase" phase, REAL8 (optional)                                                                 */
/*   "polar_angle" ellipticity polar angle, REAL8 (optional, together with polar_eccentricity may replace alpha)*/
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
  REAL8 phi0, deltaT,deltaF, 
  freq=0.0,
  quality=0.0,
  duration=0.0,
  f_low, f_max,
  hrss=1.0,
  polar_ecc=1.0,polar_angle=LAL_PI/2.,alpha=LAL_PI/2.; 
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
  */
  if(LALInferenceCheckVariable(model->params,"alpha"))
    {
      alpha=*(REAL8*) LALInferenceGetVariable(model->params, "alpha");
      if (!extraParams) extraParams=XLALSimBurstCreateExtraParam("alpha",alpha);
      else XLALSimBurstAddExtraParam(&extraParams,"alpha",alpha);
      polar_angle=alpha;
  } 
  
  if(LALInferenceCheckVariable(model->params,"phase"))
    {
      phi0=*(REAL8*) LALInferenceGetVariable(model->params, "phase");
      if (!extraParams) extraParams=XLALSimBurstCreateExtraParam("phase",phi0);
      else XLALSimBurstAddExtraParam(&extraParams,"phase",phi0);
  }

  /* If someone wants to use old parametrization, allow for */
  if(LALInferenceCheckVariable(model->params,"polar_angle"))
    polar_angle=*(REAL8*) LALInferenceGetVariable(model->params, "polar_angle");
  if(LALInferenceCheckVariable(model->params,"polar_eccentricity"))
    polar_ecc=*(REAL8*) LALInferenceGetVariable(model->params, "polar_eccentricity");
    
    
  /* Check if fLow is a model parameter, otherwise use data structure definition */
  if(LALInferenceCheckVariable(model->params, "fLow"))
    f_low = *(REAL8*) LALInferenceGetVariable(model->params, "fLow");
  else
    f_low = model->fLow /** 0.9 */;

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
        XLALPrintError(" ERROR in LALInferenceTemplateXLALSimBurstChooseWaveform(). errnum=%d\n",errnum );
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

    XLAL_TRY(ret=XLALSimBurstChooseTDWaveformFromCache(&hplus, &hcross,deltaT,freq,quality,duration,f_low,f_max,hrss,polar_angle,polar_ecc,extraParams,approximant,model->burstWaveformCache), errnum);
    XLALSimBurstDestroyExtraParam(extraParams);
    if (ret == XLAL_FAILURE || hplus == NULL || hcross == NULL)
      {
        XLALPrintError(" ERROR in XLALSimBurstChooseTDWaveform(): error generating waveform. errnum=%d\n",errnum );
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

void LALInferenceTemplatePrincipalComp(LALInferenceModel *model)
{

    UINT4 i=0, lower=120;
    REAL8 hrss = 1.0e-4;

    /* Principle Component Coefficients */
    REAL8 betas[model->pcs->nPCs];

    if (model->pcs->nPCs>=1)
        betas[0] = *(REAL8*) LALInferenceGetVariable(model->params, "beta1");
    if (model->pcs->nPCs>=2)
        betas[1] = *(REAL8*) LALInferenceGetVariable(model->params, "beta2");
    if (model->pcs->nPCs>=3)
        betas[2] = *(REAL8*) LALInferenceGetVariable(model->params, "beta3");
    if (model->pcs->nPCs>=4)
        betas[3] = *(REAL8*) LALInferenceGetVariable(model->params, "beta4");
    if (model->pcs->nPCs>=5)
        betas[4] = *(REAL8*) LALInferenceGetVariable(model->params, "beta5");
    if (model->pcs->nPCs>=6)
        betas[5] = *(REAL8*) LALInferenceGetVariable(model->params, "beta6");
    if (model->pcs->nPCs>=7)
        betas[6] = *(REAL8*) LALInferenceGetVariable(model->params, "beta7");
    if (model->pcs->nPCs>=8)
        betas[7] = *(REAL8*) LALInferenceGetVariable(model->params, "beta8");
    if (model->pcs->nPCs>=9)
        betas[8] = *(REAL8*) LALInferenceGetVariable(model->params, "beta9");
    if (model->pcs->nPCs>=10)
        betas[9] = *(REAL8*) LALInferenceGetVariable(model->params, "beta10");
    if (model->pcs->nPCs>=11)
        betas[10] = *(REAL8*) LALInferenceGetVariable(model->params, "beta11");
    if (model->pcs->nPCs>=12)
        betas[11] = *(REAL8*) LALInferenceGetVariable(model->params, "beta12");
    if (model->pcs->nPCs>=13)
        betas[12] = *(REAL8*) LALInferenceGetVariable(model->params, "beta13");
    if (model->pcs->nPCs>=14)
        betas[13] = *(REAL8*) LALInferenceGetVariable(model->params, "beta14");
    if (model->pcs->nPCs>=15)
        betas[14] = *(REAL8*) LALInferenceGetVariable(model->params, "beta15");
    if (model->pcs->nPCs>=16)
        betas[15] = *(REAL8*) LALInferenceGetVariable(model->params, "beta16");
    if (model->pcs->nPCs>=17)
        betas[16] = *(REAL8*) LALInferenceGetVariable(model->params, "beta17");
    if (model->pcs->nPCs>=18)
        betas[17] = *(REAL8*) LALInferenceGetVariable(model->params, "beta18");
    if (model->pcs->nPCs>=19)
        betas[18] = *(REAL8*) LALInferenceGetVariable(model->params, "beta19");
    if (model->pcs->nPCs==20)
        betas[19] = *(REAL8*) LALInferenceGetVariable(model->params, "beta20");

    /* Build Template */

    /* Set to zero below flow (30Hz) */
    for(i = 0; i < lower; i++) {
      model->freqhPlus->data->data[i] = 0.0;
    }

    for(i=lower; i<model->freqhPlus->data->length; i++){

        model->freqhPlus->data->data[i] =
            betas[0] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 0)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 0)))
            + betas[1] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 1)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 1)))
            + betas[2] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 2)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 2)))
            + betas[3] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 3)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 3)))
            + betas[4] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 4)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 4)))
            + betas[5] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 5)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 5)))
            + betas[6] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 6)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 6)))
            + betas[7] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 7)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 7)))
            + betas[8] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 8)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 8)))
            + betas[9] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 9)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 9)))
            + betas[10] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 10)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 10)))
            + betas[11] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 11)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 11)))
            + betas[12] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 12)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 12)))
            + betas[13] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 13)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 13)))
            + betas[14] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 14)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 14)))
            + betas[15] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 15)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 15)))
            + betas[16] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 16)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 16)))
            + betas[17] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 17)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 17)))
            + betas[18] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 18)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 18)))
            + betas[19] * (GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 19)) + I*GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, i, 19)));

        /* XXX: hrss is only really applied in the likelihood function */
        model->freqhPlus->data->data[i] *= hrss;

        /* FIXME: do something about the cross polarisation!!! */
        model->freqhCross->data->data[i] = 0.0;

    }

    model->domain = LAL_SIM_DOMAIN_FREQUENCY;

    return;


}

void LALInferenceTemplatePrincipalCompSpec(LALInferenceModel *model)
{
  //fprintf(stderr, "Beginning template function\n");
  INT4 lower = 1, start = 0;//120, i = 0;
  INT4 times_N = model->spechPlus->times_N;
  INT4 freq_N = model->spechPlus->freq_N;
  REAL8 hrss = 1.0;

  /* Principle Component Coefficients */
  REAL8 betas[model->pcs->nPCs];

  if (model->pcs->nPCs>=1)
    betas[0] = *(REAL8*) LALInferenceGetVariable(model->params, "beta1");
  if (model->pcs->nPCs>=2)
    betas[1] = *(REAL8*) LALInferenceGetVariable(model->params, "beta2");
  if (model->pcs->nPCs>=3)
    betas[2] = *(REAL8*) LALInferenceGetVariable(model->params, "beta3");
  if (model->pcs->nPCs>=4)
    betas[3] = *(REAL8*) LALInferenceGetVariable(model->params, "beta4");
  if (model->pcs->nPCs>=5)
    betas[4] = *(REAL8*) LALInferenceGetVariable(model->params, "beta5");
  if (model->pcs->nPCs>=6)
    betas[5] = *(REAL8*) LALInferenceGetVariable(model->params, "beta6");
  if (model->pcs->nPCs>=7)
    betas[6] = *(REAL8*) LALInferenceGetVariable(model->params, "beta7");
  if (model->pcs->nPCs>=8)
    betas[7] = *(REAL8*) LALInferenceGetVariable(model->params, "beta8");
  if (model->pcs->nPCs>=9)
    betas[8] = *(REAL8*) LALInferenceGetVariable(model->params, "beta9");
  if (model->pcs->nPCs>=10)
    betas[9] = *(REAL8*) LALInferenceGetVariable(model->params, "beta10");
  if (model->pcs->nPCs>=11)
    betas[10] = *(REAL8*) LALInferenceGetVariable(model->params, "beta11");
  if (model->pcs->nPCs>=12)
    betas[11] = *(REAL8*) LALInferenceGetVariable(model->params, "beta12");
  if (model->pcs->nPCs>=13)
    betas[12] = *(REAL8*) LALInferenceGetVariable(model->params, "beta13");
  if (model->pcs->nPCs>=14)
    betas[13] = *(REAL8*) LALInferenceGetVariable(model->params, "beta14");
  if (model->pcs->nPCs>=15)
    betas[14] = *(REAL8*) LALInferenceGetVariable(model->params, "beta15");
  if (model->pcs->nPCs>=16)
    betas[15] = *(REAL8*) LALInferenceGetVariable(model->params, "beta16");
  if (model->pcs->nPCs>=17)
    betas[16] = *(REAL8*) LALInferenceGetVariable(model->params, "beta17");
  if (model->pcs->nPCs>=18)
    betas[17] = *(REAL8*) LALInferenceGetVariable(model->params, "beta18");
  if (model->pcs->nPCs>=19)
    betas[18] = *(REAL8*) LALInferenceGetVariable(model->params, "beta19");
  if (model->pcs->nPCs==20)
    betas[19] = *(REAL8*) LALInferenceGetVariable(model->params, "beta20");

  //fprintf(stderr, "betas [0] = %e , [1] = %e , [2] = %e. [3] = %e , [4] = %e , [5] = %e , [6] = %e , [7] = %e , [8] = %e\n", betas[0], betas[1], betas[2], betas[3], betas[4], betas[5], betas[6], betas[7], betas[8]);

  /* Build Template */
  /* Set to zero below flow */
  for (int t = 0; t < times_N; t++) {
    for (int f = 0; f < lower; f++) {
      model->spechPlus->mag_data[f][t] = 0.0;
      model->spechCross->mag_data[f][t] = 0.0;
    }
  }
  //fprintf(stderr, "Beginning template 2nd double for loop\n");
  for (int t = 0; t < times_N; t++) {
    start = t * freq_N;
    for (int f = lower; f < freq_N; f++) {
      model->spechPlus->mag_data[f][t] = 
	betas[0] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 0)) + betas[1] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 1))
	+ betas[2] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 2)) + betas[3] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 3))
	+ betas[4] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 4)) + betas[5] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 5))
	+ betas[6] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 6)) + betas[7] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 7))
	+ betas[8] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 8)) + betas[9] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 9))
	+ betas[10] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 10)) + betas[11] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 11))
	+ betas[12] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 12)) + betas[13] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 13))
	+ betas[14] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 14)) + betas[15] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 15))
	+ betas[16] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 16)) + betas[17] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 17))
	+ betas[18] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 18)) + betas[19] * GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 19));
      
	/*fprintf(stderr, "betas[0] = %e , gsl_thing = %e\n", betas[0], GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 0)));
	  fprintf(stderr, "betas[2] = %e , gsl_thing = %e\n", betas[2], GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, start + f, 2)));*/
      //fprintf(stderr, "spec entry = %e\n", model->spechPlus->mag_data[f][t]);

      //fprintf(stderr, "mag_data[%d][%d] = %e\n", f, t, model->spechPlus->mag_data[f][t]);
      /* FIXME: do something about the cross polarization!!! */
      model->spechCross->mag_data[f][t] = 0;

      /* XXX: hrss is only really applied in the likelihood function */
      model->spechPlus->mag_data[f][t] *= hrss;
      model->spechCross->mag_data[f][t] *= hrss;
    }
  }
  /*fprintf(stderr, "[0][0]real = %e , [0][0]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 0)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 0, 0)));
  fprintf(stderr, "[1][0]real = %e , [1][0]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 1, 0)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 1, 0)));
  fprintf(stderr, "[2][0]real = %e , [2][0]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 2, 0)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 2, 0)));
  fprintf(stderr, "[3][0]real = %e , [3][0]imag = %e\n", GSL_REAL(gsl_matrix_complex_get(model->pcs->pcs_plus, 3, 0)), GSL_IMAG(gsl_matrix_complex_get(model->pcs->pcs_plus, 3, 0)));*/

  // Rescaling to have hrss ~ 1
//  double S[65]; // Average PSD.  Essentially doing a Welch
//  double x = 0;
//  for (int f = 0; f < 65; f++) {
//    x = 0;
//    for (int t = 0; t < 936; t++) {
//      x += model->spechPlus->mag_data[f][t] * model->spechPlus->mag_data[f][t];
//    }
//    x /= 936.0;
//    S[f] = x;
    //fprintf(stderr, "S[f], x = %e\n", x);
//  }

//  for (int f = 0; f < 65; f++) {
//    hrss += S[f];
//  }
//  hrss *= 32.0;
//  hrss = sqrt(hrss);
//  double factor = 1.0 / hrss;

  //fprintf(stderr, "hrss = %e\n", hrss);
  //fprintf(stderr, "factor = %e\n", factor);

//  for (int f = 0; f < 65; f++) {
//    for (int t = 0; t < 936; t++) {
//      model->spechPlus->mag_data[f][t] = factor * model->spechPlus->mag_data[f][t];
//    }
//  }

  //model->domain = LAL_SIM_DOMAIN_SPECTROGRAM;
  model->domain = LAL_SIM_DOMAIN_FREQUENCY;
  //fprintf(stderr, "Finishing template function\n"); 
    return;

 }


