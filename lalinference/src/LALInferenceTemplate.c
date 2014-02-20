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

extern int newswitch; //temporay global variable to use the new LALSTPN
static void q2masses(double mc, double q, double *m1, double *m2);


void LALInferenceTemplateNullFreqdomain(LALInferenceIFOData *IFOdata)
/**********************************************/
/* returns a frequency-domain 'null' template */
/* (all zeroes, implying no signal present).  */
/**********************************************/
{
  UINT4 i;
  if ((IFOdata->freqModelhPlus==NULL) || (IFOdata->freqModelhCross==NULL)) {
    XLALPrintError(" ERROR in templateNullFreqdomain(): encountered unallocated 'freqModelhPlus/-Cross'.\n");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
  for (i=0; i<IFOdata->freqModelhPlus->data->length; ++i){
    IFOdata->freqModelhPlus->data->data[i] = 0.0;
    IFOdata->freqModelhCross->data->data[i] = 0.0;
  }
  IFOdata->modelDomain = LAL_SIM_DOMAIN_FREQUENCY;
  return;
}



void LALInferenceTemplateNullTimedomain(LALInferenceIFOData *IFOdata)
/*********************************************/
/* returns a time-domain 'null' template     */
/* (all zeroes, implying no signal present). */
/*********************************************/
{
  UINT4 i;
  if ((IFOdata->timeModelhPlus==NULL) || (IFOdata->timeModelhCross==NULL)) {
    XLALPrintError(" ERROR in templateNullTimedomain(): encountered unallocated 'timeModelhPlus/-Cross'.\n");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
  for (i=0; i<IFOdata->timeModelhPlus->data->length; ++i){
    IFOdata->timeModelhPlus->data->data[i]  = 0.0;
    IFOdata->timeModelhCross->data->data[i] = 0.0;
  }
  IFOdata->modelDomain = LAL_SIM_DOMAIN_TIME;
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




void LALInferenceTemplateSineGaussian(LALInferenceIFOData *IFOdata)
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
/* Required (`IFOdata->modelParams') parameters are:                                    */
/*   - "time"       (the "mu" parameter of the Gaussian part; REAL8, GPS sec.)          */
/*   - "sigma"      (width, the "sigma" parameter of the Gaussian part; REAL8, seconds) */
/*   - "frequency"  (frequency of the sine part; REAL8, Hertz)                          */
/*   - "phase"      (phase (at above "mu"); REAL8, radians)                             */
/*   - "amplitude"  (amplitude, REAL8)                                                  */
/****************************************************************************************/
{
  double endtime  = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "time");       /* time parameter ("mu"), GPS sec.  */
  double sigma = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "sigma");      /* width parameter, seconds         */
  double f     = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "frequency");  /* frequency, Hz                    */
  double phi   = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phase");      /* phase, rad                       */
  double a     = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "amplitude");  /* amplitude                        */
  double t, tsigma, twopif = LAL_TWOPI*f;
  double epochGPS = XLALGPSGetREAL8(&(IFOdata->timeData->epoch));
  unsigned long i;
  if (sigma <= 0.0) {
    fprintf(stderr, " ERROR in templateSineGaussian(): zero or negative \"sigma\" parameter (sigma=%e).\n", sigma);
    exit(1);
  }
  if (f < 0.0)
    fprintf(stderr, " WARNING in templateSineGaussian(): negative \"frequency\" parameter (f=%e).\n", f);
  if (a < 0.0)
    fprintf(stderr, " WARNING in templateSineGaussian(): negative \"amplitude\" parameter (a=%e).\n", a);
  for (i=0; i<IFOdata->timeData->data->length; ++i){
    t = ((double)i)*IFOdata->timeData->deltaT + (epochGPS-endtime);  /* t-mu         */
    tsigma = t/sigma;                                             /* (t-mu)/sigma */
    if (fabs(tsigma) < 5.0)   /*  (only do computations within a 10 sigma range)  */
      IFOdata->timeModelhPlus->data->data[i] = a * exp(-0.5*tsigma*tsigma) * sin(twopif*t+phi);
    else 
      IFOdata->timeModelhPlus->data->data[i] = 0.0;
    IFOdata->timeModelhCross->data->data[i] = 0.0;
  }
  IFOdata->modelDomain = LAL_SIM_DOMAIN_TIME;
  return;
}



void LALInferenceTemplateDampedSinusoid(LALInferenceIFOData *IFOdata)
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
/* Required (`IFOdata->modelParams') parameters are:                          */
/*   - "time"       (the instant at which the signal starts; REAL8, GPS sec.) */
/*   - "tau"        (width parameter; REAL8, seconds)                         */
/*   - "frequency"  (frequency of the sine part; REAL8, Hertz)                */
/*   - "amplitude"  (amplitude, REAL8)                                        */
/******************************************************************************/
{
  double endtime  = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "time");       /* time parameter ("mu"), GPS sec.  */
  double tau   = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "tau");        /* width parameter, seconds         */
  double f     = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "frequency");  /* frequency, Hz                    */
  double a     = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "amplitude");  /* amplitude                        */
  double t, ttau, twopif = LAL_TWOPI*f;
  double epochGPS = XLALGPSGetREAL8(&(IFOdata->timeData->epoch));
  unsigned long i;
  if (tau <= 0.0) {
    fprintf(stderr, " ERROR in templateDampedSinusoid(): zero or negative \"tau\" parameter (tau=%e).\n", tau);
    exit(1);
  }
  if (f < 0.0)
    fprintf(stderr, " WARNING in templateDampedSinusoid(): negative \"frequency\" parameter (f=%e).\n", f);
  for (i=0; i<IFOdata->timeData->data->length; ++i){
    t = ((double)i)*IFOdata->timeData->deltaT + (epochGPS-endtime);  /* t-mu       */
    if ((t>0.0) && ((ttau=t/tau) < 10.0)) /*  (only do computations within a 10 tau range)  */
      IFOdata->timeModelhPlus->data->data[i] = a * exp(-ttau) * sin(twopif*t);
    else 
      IFOdata->timeModelhPlus->data->data[i] = 0.0;
    IFOdata->timeModelhCross->data->data[i] = 0.0;
  }
  IFOdata->modelDomain = LAL_SIM_DOMAIN_TIME;
  return;
}



void LALInferenceTemplateSinc(LALInferenceIFOData *IFOdata)
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
/* Required (`IFOdata->modelParams') parameters are:                         */
/*   - "time"       (the instant at which the signal peaks; REAL8, GPS sec.) */
/*   - "frequency"  (frequency of the sine part; REAL8, Hertz)               */
/*   - "amplitude"  (amplitude, REAL8)                                       */
/*****************************************************************************/
{
  double endtime  = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "time");       /* time parameter ("mu"), GPS sec.  */
  double f     = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "frequency");  /* frequency, Hz                    */
  double a     = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "amplitude");  /* amplitude                        */
  double t, sinArg, sinc, twopif = LAL_TWOPI*f;
  double epochGPS = XLALGPSGetREAL8(&(IFOdata->timeData->epoch));
  unsigned long i;
  if (f < 0.0)
    fprintf(stderr, " WARNING in templateSinc(): negative \"frequency\" parameter (f=%e).\n", f);
  for (i=0; i<IFOdata->timeData->data->length; ++i){
    t = ((double)i)*IFOdata->timeData->deltaT + (epochGPS-endtime);  /* t-mu       */
    sinArg = twopif*t;
    sinc = (sinArg==0.0) ? 1.0 : sin(sinArg)/sinArg;    
    IFOdata->timeModelhPlus->data->data[i] = a * sinc;
    IFOdata->timeModelhCross->data->data[i] = 0.0;
  }
  IFOdata->modelDomain = LAL_SIM_DOMAIN_TIME;
  return;
}


void LALInferenceTemplateASinOmegaT(LALInferenceIFOData *IFOdata)
/************************************************************/
/* Trivial h(t)=A*sin(Omega*t) template						*/
/*  Required (`IFOdata->modelParams') parameters are:       */
/*   - "A"       (dimensionless amplitude, REAL8)			*/
/*   - "Omega"   (frequency; REAL8, radians/sec)            */
/************************************************************/
{
  double A		= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "A");				/* dim-less	   */
  double Omega	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "Omega");			/* rad/sec     */
  double t;
  double epochGPS = XLALGPSGetREAL8(&(IFOdata->timeData->epoch));	

  unsigned long i;
  for (i=0; i<IFOdata->timeData->data->length; ++i){
    t = ((double)i)*IFOdata->timeData->deltaT + (epochGPS);  /* t-mu       */   
    IFOdata->timeModelhPlus->data->data[i] = A * sin(Omega*t);
    IFOdata->timeModelhCross->data->data[i] = 0.0;
  }
  IFOdata->modelDomain = LAL_SIM_DOMAIN_TIME;
  return;
}

void LALInferenceTemplateXLALSimInspiralChooseWaveform(LALInferenceIFOData *IFOdata)
/*************************************************************************************************************************/
/* Wrapper for LALSimulation waveforms:						                                                             */
/* XLALSimInspiralChooseFDWaveform() and XLALSimInspiralChooseTDWaveform().                                              */
/*                                                                                                                       */
/*  IFOdata->modelParams parameters are:										                                         */
/*  - "name" description; type OPTIONAL (default value)										                             */
/*										                                                                                 */
/*   MODEL PARAMETERS										                                                             */
/*   - "LAL_APPROXIMANT"	  Approximant;        Approximant                                                            */
/*   - "LAL_PNORDER"        Phase PN order;     INT4                                                                     */
/*   - "LAL_AMPORDER"       Amplitude PN order; INT4 OPTIONAL (-1)                                                       */
/*   - "LALINFERENCE_FRAME" reference frame;    LALInferenceFrame OPTIONAL (LALINFERENCE_FRAME_RADIATION)                */
/*   - "spinO"              Spin order;         LALSimInspiralSpinOrder OPTIONAL (LAL_SIM_INSPIRAL_SPIN_ORDER_DEFAULT)   */
/*   - "tideO"              Tidal order;        LALSimInspiralTidalOrder OPTIONAL (LAL_SIM_INSPIRAL_TIDAL_ORDER_DEFAULT) */
/*   - "fRef"               frequency at which the (frequency dependent) parameters are defined; REAL8 OPTIONAL (0.0)    */
/*   - "fLow"               lower frequency bound; REAL8 OPTIONAL (IFOdata->fLow)                                        */
/*                                                                                                                       */
/*   MASS PARAMETERS; either:                                                                                            */
/*      - "mass1"           mass of object 1 in solar mass; REAL8								                         */
/*      - "mass2"		        mass of object 1 in solar mass; REAL8								                     */
/*      OR                                                                                                               */
/*      - "chirpmass"       chirpmass in solar mass; REAL8                                                               */
/*      - "asym_massratio"  asymmetric mass ration m2/m1, 0<asym_massratio<1; REAL8                                      */
/*      OR                                                                                                               */
/*      - "chirpmass"       chirpmass in solar mass; REAL8                                                               */
/*      - "massratio"       symmetric mass ratio (m1*m2)/(m1+m2)^2; REAL8                                                */
/*                                                                                                                       */
/*   ORIENTATION AND SPIN PARAMETERS                                                                                     */
/*   - "phi0"               reference phase as per LALSimulation convention; REAL8                                       */
/*   - "distance"           distance in Mpc                                                                              */
/*   - if LALINFERENCE_FRAME == LALINFERENCE_FRAME_SYSTEM  (default)                                                     */
/*      - "theta_JN");      zenith angle between J and N in radians;            REAL8                                    */
/*      - "phi_JL");        azimuthal angle of L_N on its cone about J radians; REAL8                                    */
/*      - "tilt_spin1");    zenith angle between S1 and LNhat in radians;       REAL8                                    */
/*      - "tilt_spin2");    zenith angle between S2 and LNhat in radians;       REAL8                                    */
/*      - "phi12");         difference in azimuthal angle between S1, S2 in radians;   REAL8                             */
/*   - else if LALINFERENCE_FRAME == LALINFERENCE_FRAME_RADIATION                                                        */
/*      - "inclination"	    inclination angle L.N in radians;                            REAL8                           */
/*      - "theta_spin1"     polar angle of spin 1, default to the spin aligned case;     REAL8 OPTIONAL (inclination)    */
/*      - "phi_spin1"       azimuthal angle of spin 1, default to the spin aligned case; REAL8  OPTIONAL (0.0)           */
/*      - "theta_spin2"     polar angle of spin 2, default to the spin aligned case;     REAL8 OPTIONAL (inclination)    */
/*      - "phi_spin2"       azimuthal angle of spin 1, default to the spin aligned case; REAL8  OPTIONAL (0.0)           */
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
/*   IFOdata needs to also contain:                                                                                      */
/*   - IFOdata->fLow Unless  - "fLow" OPTIONAL                                                                           */
/*   - IFOdata->timeData                                                                                                 */
/*      - IFOdata->timeData->deltaT                                                                                      */
/*   - if IFOdata->modelDomain == LAL_SIM_DOMAIN_FREQUENCY                                                               */
/*      - IFOdata->freqData                                                                                              */
/*          - IFOdata->freqData->deltaF                                                                                  */
/*      - IFOdata->freqModelhCross                                                                                       */
/*      - IFOdata->freqModelhPlus                                                                                        */
/*   - else                                                                                                              */
/*      - IFOdata->timeModelhPlus                                                                                        */
/*      - IFOdata->timeModelhCross                                                                                       */
/*************************************************************************************************************************/
{

  Approximant approximant = (Approximant) 0;
  INT4 order=-1;
  INT4 amporder;
  LALInferenceFrame frame=LALINFERENCE_FRAME_SYSTEM;

  unsigned long	i;
  static int sizeWarning = 0;
  int ret=0;
  INT4 errnum=0;
  REAL8 instant;
  
  
  REAL8TimeSeries *hplus=NULL;  /**< +-polarization waveform [returned] */
  REAL8TimeSeries *hcross=NULL; /**< x-polarization waveform [returned] */
  COMPLEX16FrequencySeries *hptilde=NULL, *hctilde=NULL;
  REAL8 mc;
  REAL8 phi0, deltaT, m1, m2, spin1x=0.0, spin1y=0.0, spin1z=0.0, spin2x=0.0, spin2y=0.0, spin2z=0.0, f_low, f_start, distance, inclination;
  
  REAL8 *m1_p,*m2_p;
  REAL8 deltaF, f_max;
  
  if (LALInferenceCheckVariable(IFOdata->modelParams, "LAL_APPROXIMANT"))
    approximant = *(Approximant*) LALInferenceGetVariable(IFOdata->modelParams, "LAL_APPROXIMANT");
  else {
    XLALPrintError(" ERROR in templateLALGenerateInspiral(): (INT4) \"LAL_APPROXIMANT\" parameter not provided!\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }
	
  if (LALInferenceCheckVariable(IFOdata->modelParams, "LAL_PNORDER"))
    order = *(INT4*) LALInferenceGetVariable(IFOdata->modelParams, "LAL_PNORDER");
  else {
    XLALPrintError(" ERROR in templateLALGenerateInspiral(): (INT4) \"LAL_PNORDER\" parameter not provided!\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }

  /* Explicitly set the default amplitude order if one is not specified.
   *   This serves two purposes:
   *     1) The default behavior of the code won't change unexpectedly due to changes in LALSimulation.
   *     2) We need to know the amplitude order in order to set the starting frequency of the waveform properly. */
  if (LALInferenceCheckVariable(IFOdata->modelParams, "LAL_AMPORDER"))
    amporder = *(INT4*) LALInferenceGetVariable(IFOdata->modelParams, "LAL_AMPORDER");
  else
    amporder = -1;

  if (LALInferenceCheckVariable(IFOdata->modelParams, "LALINFERENCE_FRAME"))
    frame = *(LALInferenceFrame*) LALInferenceGetVariable(IFOdata->modelParams, "LALINFERENCE_FRAME");

  REAL8 fRef = 100.0;
  if (LALInferenceCheckVariable(IFOdata->modelParams, "fRef")) fRef = *(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams, "fRef");

  REAL8 fTemp = fRef;

  if(LALInferenceCheckVariable(IFOdata->modelParams,"chirpmass"))
    {
      mc  = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "chirpmass");
      if (LALInferenceCheckVariable(IFOdata->modelParams,"asym_massratio")) {
	REAL8 q = *(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"asym_massratio");
	q2masses(mc, q, &m1, &m2);
      } else {
	REAL8 eta = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "massratio");
	mc2masses(mc, eta, &m1, &m2);
      }
    }
  else if((m1_p=(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams, "mass1")) && (m2_p=(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams, "mass2")))
    {
      m1=*m1_p;
      m2=*m2_p;
    }
  else
    {
      fprintf(stderr,"No mass parameters found!");
      exit(0);
    }

  distance	= LALInferenceGetREAL8Variable(IFOdata->modelParams,"distance")* LAL_PC_SI * 1.0e6;        /* distance (1 Mpc) in units of metres */

  /* Default to spinless signals if spin amplitude are not present in modelParams */
  REAL8 a_spin1		= 0.0;
  if(LALInferenceCheckVariable(IFOdata->modelParams, "a_spin1"))    a_spin1   = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "a_spin1");
  REAL8 a_spin2    = 0.0;
  if(LALInferenceCheckVariable(IFOdata->modelParams, "a_spin2"))    a_spin2   = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "a_spin2");

  phi0		= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phase"); /* START phase as per lalsimulation convention*/

  /* Check if fLow is a model parameter, otherwise use data structure definition */
  if(LALInferenceCheckVariable(IFOdata->modelParams, "fLow"))
    f_low = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "fLow");
  else
    f_low = IFOdata->fLow /** 0.9 */;

  f_start = fLow2fStart(f_low, amporder, approximant);
  f_max = 0.0; /* for freq domain waveforms this will stop at ISCO. Previously found using IFOdata->fHigh causes NaNs in waveform (see redmine issue #750)*/

  int aligned_spins=0;
  /* We first check if we are deadling with a spin-aligned only template, for which we use "spin1" and "spin2" names */
  if(LALInferenceCheckVariable(IFOdata->modelParams, "spin1")){
    spin1z= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "spin1");
    spin1x=0.0;
    spin1y=0.0;
    aligned_spins+=1;
  }
  if(LALInferenceCheckVariable(IFOdata->modelParams, "spin2")){
    spin2z= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "spin2");
    spin2x=0.0;
    spin2y=0.0;
    aligned_spins+=1;
  }
    
  if (aligned_spins==0){
    /* Template is not spin-aligned only.
     * Set the all other spins variables (that is an overkill. We can just check if there are spins in the first place, and if there aren't don't bother calculating them (they'll be zero) ) */
    REAL8 phiJL=0.0;
    REAL8 tilt1=0.0;
    REAL8 tilt2=0.0;
    REAL8 phi12=0.0;
    if(frame==LALINFERENCE_FRAME_SYSTEM) {
       /* IMPORTANT NOTE: default to spin aligned case (i.e. tilt1=tilt2=0) if no angles are provided for the spins.
       * If you change this, must also change LALInferenceInitCBC.c
       */
       REAL8 thetaJN = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "theta_JN");     /* zenith angle between J and N in radians */
       if(LALInferenceCheckVariable(IFOdata->modelParams, "phi_JL"))
           phiJL = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi_JL");     /* azimuthal angle of L_N on its cone about J radians */
       if(LALInferenceCheckVariable(IFOdata->modelParams, "tilt_spin1"))
           tilt1 = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "tilt_spin1");     /* zenith angle between S1 and LNhat in radians */
       if(LALInferenceCheckVariable(IFOdata->modelParams, "tilt_spin2"))
           tilt2 = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "tilt_spin2");     /* zenith angle between S2 and LNhat in radians */
       if(LALInferenceCheckVariable(IFOdata->modelParams, "phi12"))
           phi12 = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi12");      /* difference in azimuthal angle btwn S1, S2 in radians */

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
    else if(frame==LALINFERENCE_FRAME_RADIATION){
      /* IMPORTANT NOTE: default to spin aligned case (i.e. theta1=theta2=iota) if no angles are provided for the spins.
       * If you change this, must also change LALInferenceInitCBC.c
       */
      inclination	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "inclination");	    /* inclination in radian */
      REAL8 theta_spin1	= inclination;
      if(LALInferenceCheckVariable(IFOdata->modelParams, "theta_spin1"))	theta_spin1	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "theta_spin1");
      REAL8 phi_spin1		= 0.0;
      if(LALInferenceCheckVariable(IFOdata->modelParams, "phi_spin1"))	phi_spin1	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi_spin1");
      REAL8 theta_spin2	= inclination;
      if(LALInferenceCheckVariable(IFOdata->modelParams, "theta_spin2"))	theta_spin2	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "theta_spin2");
      REAL8 phi_spin2		= 0.0;
      if(LALInferenceCheckVariable(IFOdata->modelParams, "phi_spin2"))	phi_spin2	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi_spin2");

      /* These are the components of spins in the frame N || z */
      spin1x = (a_spin1 * sin(theta_spin1) * cos(phi_spin1));
      spin1y = (a_spin1 * sin(theta_spin1) * sin(phi_spin1));
      spin1z = (a_spin1 * cos(theta_spin1));

      spin2x = (a_spin2 * sin(theta_spin2) * cos(phi_spin2));
      spin2y = (a_spin2 * sin(theta_spin2) * sin(phi_spin2));
      spin2z = (a_spin2 * cos(theta_spin2));

    }
    else {
        XLALPrintError("Error: unknown frame %i\n",frame);
        XLAL_ERROR_VOID(XLAL_EFAULT);
    }
  }
	
  REAL8 lambda1 = 0.;
  if(LALInferenceCheckVariable(IFOdata->modelParams, "lambda1")) lambda1 = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "lambda1");
  REAL8 lambda2 = 0.;
  if(LALInferenceCheckVariable(IFOdata->modelParams, "lambda2")) lambda2 = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "lambda2");
  REAL8 lambdaT = 0.;
  REAL8 dLambdaT = 0.;
  REAL8 sym_mass_ratio_eta = 0.;
  if(LALInferenceCheckVariable(IFOdata->modelParams, "lambdaT")&&LALInferenceCheckVariable(IFOdata->modelParams, "dLambdaT")){
    lambdaT = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "lambdaT");
    dLambdaT = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "dLambdaT");
    sym_mass_ratio_eta = m1*m2/((m1+m2)*(m1+m2));
    LALInferenceLambdaTsEta2Lambdas(lambdaT,dLambdaT,sym_mass_ratio_eta,&lambda1,&lambda2);
  }

  LALSimInspiralTestGRParam *nonGRparams = NULL;
  
  if (IFOdata->timeData==NULL) {
    XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'timeData'.\n");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
  deltaT = IFOdata->timeData->deltaT;
  
  
  if(IFOdata->modelDomain == LAL_SIM_DOMAIN_FREQUENCY) {
    if (IFOdata->freqData==NULL) {
      XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'freqData'.\n");
      XLAL_ERROR_VOID(XLAL_EFAULT);
    }

    deltaF = IFOdata->freqData->deltaF;
    
	XLAL_TRY(ret=XLALSimInspiralChooseFDWaveformFromCache(&hptilde, &hctilde, phi0,
            deltaF, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z,
            spin2x, spin2y, spin2z, f_start, f_max, distance, inclination,lambda1, lambda2, IFOdata->waveFlags, nonGRparams, amporder, order,
            approximant,IFOdata->waveformCache), errnum);

	if (hptilde==NULL || hptilde->data==NULL || hptilde->data->data==NULL ) {
	  XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'hptilde'.\n");
	  XLAL_ERROR_VOID(XLAL_EFAULT);
	}
	if (hctilde==NULL || hctilde->data==NULL || hctilde->data->data==NULL ) {
	  XLALPrintError(" ERROR in LALInferenceTemplateXLALSimInspiralChooseWaveform(): encountered unallocated 'hctilde'.\n");
	  XLAL_ERROR_VOID(XLAL_EFAULT);
	}
      
	COMPLEX16 *dataPtr = hptilde->data->data;

    for (i=0; i<IFOdata->freqModelhPlus->data->length; ++i) {
      dataPtr = hptilde->data->data;
      if(i < hptilde->data->length){
        IFOdata->freqModelhPlus->data->data[i] = dataPtr[i];
      }else{
        IFOdata->freqModelhPlus->data->data[i] = 0.0;
      }
    }
    for (i=0; i<IFOdata->freqModelhCross->data->length; ++i) {
      dataPtr = hctilde->data->data;
      if(i < hctilde->data->length){
        IFOdata->freqModelhCross->data->data[i] = dataPtr[i];
      }else{
        IFOdata->freqModelhCross->data->data[i] = 0.0;
      }
    }
    
    
    /* Destroy the nonGr params */
    XLALSimInspiralDestroyTestGRParam(nonGRparams);
    
    instant= (IFOdata->timeData->epoch.gpsSeconds + 1e-9*IFOdata->timeData->epoch.gpsNanoSeconds);
    LALInferenceSetVariable(IFOdata->modelParams, "time", &instant);
    
  } else {

    XLAL_TRY(ret=XLALSimInspiralChooseTDWaveformFromCache(&hplus, &hcross, phi0, deltaT,
            m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z,
            spin2x, spin2y, spin2z, f_start, fRef, distance,
            inclination, lambda1, lambda2, IFOdata->waveFlags, nonGRparams,
            amporder, order, approximant,IFOdata->waveformCache), errnum);
    XLALSimInspiralDestroyTestGRParam(nonGRparams);
    if (ret == XLAL_FAILURE || hplus == NULL || hcross == NULL)
      {
	XLALPrintError(" ERROR in XLALSimInspiralChooseWaveformFromCache(): error generating waveform. errnum=%d\n",errnum );
	for (i=0; i<IFOdata->timeData->data->length; i++){
	  IFOdata->timeModelhPlus->data->data[i] = 0.0;
	  IFOdata->timeModelhCross->data->data[i] = 0.0;
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
    size_t bufLength = IFOdata->timeData->data->length;

    /* 2*Rearth/(c*deltaT)---2 is safety factor---is the maximum time
       shift for any earth-based detector. */
    size_t maxShift = (size_t)lround(4.255e-2/hplus->deltaT); 

    /* Taper 0.4 seconds at start and end (hard-coded! in
       LALInferenceReadData.c, around line 233). */
    size_t taperLength = (size_t)lround(0.4/hplus->deltaT); 

    /* Within unsafeLength of ends of buffer, possible danger of
       wrapping and/or tapering interactions. */
    size_t unsafeLength = taperLength + maxShift;

    REAL8 desiredTc = *(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams, "time");
    REAL8 tStart = XLALGPSGetREAL8(&(IFOdata->timeModelhPlus->epoch));
    REAL8 tEnd = tStart + IFOdata->timeModelhPlus->deltaT * IFOdata->timeModelhPlus->data->length;

    if (desiredTc < tStart || desiredTc > tEnd) {
      XLALDestroyREAL8TimeSeries(hplus);
      XLALDestroyREAL8TimeSeries(hcross);

      XLAL_PRINT_ERROR("desired tc (%.4f) outside data buffer\n", desiredTc);
      XLAL_ERROR_VOID(XLAL_EDOM);
    }

    /* The nearest sample in model buffer to the desired tc. */
    size_t tcSample = (size_t)lround((desiredTc - XLALGPSGetREAL8(&(IFOdata->timeModelhPlus->epoch)))/IFOdata->timeModelhPlus->deltaT);

    /* The acutal coalescence time that corresponds to the buffer
       sample on which the waveform's tC lands. */
    REAL8 injTc = XLALGPSGetREAL8(&(IFOdata->timeModelhPlus->epoch)) + tcSample*IFOdata->timeModelhPlus->deltaT;

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
    memset(IFOdata->timeModelhPlus->data->data, 0, sizeof(REAL8)*IFOdata->timeModelhPlus->data->length);
    memset(IFOdata->timeModelhCross->data->data, 0, sizeof(REAL8)*IFOdata->timeModelhCross->data->length);
    
    /* Inject */
    memcpy(IFOdata->timeModelhPlus->data->data + bufStartIndex,
	   hplus->data->data + waveStartIndex,
	   bufWaveLength*sizeof(REAL8));
    memcpy(IFOdata->timeModelhCross->data->data + bufStartIndex,
	   hcross->data->data + waveStartIndex,
	   bufWaveLength*sizeof(REAL8));

    LALInferenceSetVariable(IFOdata->modelParams, "time", &injTc);
  }
  if ( hplus ) XLALDestroyREAL8TimeSeries(hplus);
  if ( hcross ) XLALDestroyREAL8TimeSeries(hcross);
  if ( hptilde ) XLALDestroyCOMPLEX16FrequencySeries(hptilde);
  if ( hctilde ) XLALDestroyCOMPLEX16FrequencySeries(hctilde);
  
  return;
}




void LALInferenceDumptemplateFreqDomain(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
					LALInferenceTemplateFunction templt, const char *filename)
/* de-bugging function writing (frequency-domain) template to a CSV file */
/* File contains real & imaginary parts of plus & cross components.      */
/* Template amplitude is scaled to 1Mpc distance.                        */
{
  FILE *outfile=NULL; 
  LALInferenceIFOData *dataPtr;
  double deltaT, deltaF, f;
  UINT4 i;

  LALInferenceCopyVariables(currentParams, data->modelParams);
  dataPtr = data;
  while (dataPtr != NULL) { /* this loop actually does nothing (yet) here. */
    templt(data);
    if (data->modelDomain == LAL_SIM_DOMAIN_TIME)
      LALInferenceExecuteFT(data);

    outfile = fopen(filename, "w");
    /*fprintf(outfile, "f PSD dataRe dataIm signalPlusRe signalPlusIm signalCrossRe signalCrossIm\n");*/
    fprintf(outfile, "\"f\",\"PSD\",\"signalPlusRe\",\"signalPlusIm\",\"signalCrossRe\",\"signalCrossIm\"\n");
    deltaT = dataPtr->timeData->deltaT;
    deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT);
    for (i=0; i<data->freqModelhPlus->data->length; ++i){
      f = ((double) i) * deltaF;
      fprintf(outfile, "%f,%e,%e,%e,%e,%e\n",
              f, data->oneSidedNoisePowerSpectrum->data->data[i],
              /*data->freqData->data->data[i].re, data->freqData->data->data[i].im,*/
              creal(data->freqModelhPlus->data->data[i]),
              cimag(data->freqModelhPlus->data->data[i]),
              creal(data->freqModelhCross->data->data[i]),
              cimag(data->freqModelhCross->data->data[i]));
    }
    fclose(outfile);
    dataPtr = NULL;
  }
  fprintf(stdout, " wrote (frequency-domain) template to CSV file \"%s\".\n", filename);
}


void LALInferenceDumptemplateTimeDomain(LALInferenceVariables *currentParams, LALInferenceIFOData * data, 
					LALInferenceTemplateFunction templt, const char *filename)
/* de-bugging function writing (frequency-domain) template to a CSV file */
/* File contains real & imaginary parts of plus & cross components.      */
/* Template amplitude is scaled to 1Mpc distance.                        */
{
  FILE *outfile=NULL; 
  LALInferenceIFOData *dataPtr;
  double deltaT, t, epoch; // deltaF - set but not used
  UINT4 i;

  LALInferenceCopyVariables(currentParams, data->modelParams);
  dataPtr = data;
  while (dataPtr != NULL) { /* this loop actually does nothing (yet) here. */
    templt(data);
    if (data->modelDomain == LAL_SIM_DOMAIN_FREQUENCY)
      LALInferenceExecuteInvFT(data);

    outfile = fopen(filename, "w");
    fprintf(outfile, "\"t\",\"signalPlus\",\"signalCross\"\n");
    deltaT = dataPtr->timeData->deltaT;
    //deltaF = 1.0 / (((double)dataPtr->timeData->data->length) * deltaT); - set but not used
    epoch = XLALGPSGetREAL8(&data->timeData->epoch);
    for (i=0; i<data->timeModelhPlus->data->length; ++i){
      t =  epoch + ((double) i) * deltaT;
      fprintf(outfile, "%f,%e,%e\n",
              t,
              data->timeModelhPlus->data->data[i],
              data->timeModelhCross->data->data[i]);
    }
    fclose(outfile);
    dataPtr = NULL;
  }
  fprintf(stdout, " wrote (time-domain) template to CSV file \"%s\".\n", filename);
}


