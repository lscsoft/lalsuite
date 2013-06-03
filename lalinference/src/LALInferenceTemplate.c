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

#define LAL_USE_OLD_COMPLEX_STRUCTS
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
static void destroyCoherentGW( CoherentGW *waveform );
static void q2eta(double q, double *eta);
static void q2masses(double mc, double q, double *m1, double *m2);
static REAL8 fLow2fStart(REAL8 fLow, INT4 ampOrder);


void LALInferenceTemplateStatPhase(LALInferenceIFOData *IFOdata)
/*************************************************************/
/* returns the (analytic) frequency-domain template.         */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* 2.5PN stationary-phase approximation                      */
/* following  Tanaka/Tagoshi (2000), Phys.Rev.D 62(8):082001 */
/* or Christensen/Meyer (2001), Phys.Rev.D 64(2):022001.     */
/* By supplying the optional IFOdata->modelParams "PNOrder"  */
/* parameter, one may request a 2.0PN (instead of 2.5PN)     */
/* template.                                                 */
/* Signal's amplitude corresponds to a luminosity distance   */
/* of 1 Mpc; re-scaling will need to be taken care of e.g.   */
/* in the calling likelihood function.                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *********************************/
/* Required (`IFOdata->modelParams') parameters are:                                         */
/*   - "chirpmass"      (REAL8,units of solar masses)                                        */
/*   - "massratio"      (symmetric mass ratio:  0 < eta <= 0.25, REAL8) <or asym_massratio>  */
/*   - "asym_massratio" (asymmetric mass ratio:  0 < q <= 1.0, REAL8)   <or massratio>       */
/*   - "phase"          (coalescence phase; REAL8, radians)                                  */
/*   - "time"           (coalescence time; REAL8, GPS seconds)                               */
/*   - "inclination"    (inclination angle, REAL8, radians)                                  */
/*********************************************************************************************/
{
  double mc   = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "chirpmass");
  double phi  = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phase");
  double iota = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "inclination");
  double tc   = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "time");
  
  double eta; 
  double fLow;
  if (LALInferenceCheckVariable(IFOdata->modelParams,"asym_massratio")) {
    double q = *(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"asym_massratio");
    q2eta(q, &eta);
  }
  else
    eta = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "massratio");
 
  double PNOrder = 2.5;  /* (default) */
  double fraction = (0.5+sqrt(0.25-eta)) / (0.5-sqrt(0.25-eta));
  double mt = mc * ((pow(1.0+fraction,0.2) / pow(fraction,0.6))
                    + (pow(1.0+1.0/fraction,0.2) / pow(1.0/fraction,0.6))) *  LAL_MSUN_SI;  /* (total mass, kg (!)) */
  double log_q   = log(mt) + log(LAL_PI) + log(LAL_G_SI) - 3.0*log((double) LAL_C_SI);
  double log_eta = log(eta);
  double a[5];
  double ampliConst, plusCoef, crossCoef;
  double dataStart, NDeltaT, phaseArg;
  double f, f01, f02, f04, f06, f07, f10, Psi, twopitc;
  double plusRe, plusIm, crossRe, crossIm;
  UINT4 i, lower, upper;

  if (IFOdata->timeData==NULL){
    XLALPrintError(" ERROR in templateStatPhase(): encountered unallocated 'timeData'.\n");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
  if ((IFOdata->freqModelhPlus==NULL) || (IFOdata->freqModelhCross==NULL)) {
    XLALPrintError(" ERROR in templateStatPhase(): encountered unallocated 'freqModelhPlus/-Cross'.\n");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
  if (LALInferenceCheckVariable(IFOdata->modelParams, "PNOrder"))
    PNOrder = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "PNOrder");
  if ((PNOrder!=2.5) && (PNOrder!=2.0)) {
    XLALPrintError(" ERROR in templateStatPhase(): only PN orders 2.0 or 2.5 allowed.");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
  ampliConst  = 0.5*log(5.0) + (5.0/6.0)*log(LAL_G_SI) - log(2.0) - 0.5*log(6.0) - (2.0/3.0)*log(LAL_PI) - 1.5*log((double)LAL_C_SI);
  ampliConst  = exp(ampliConst + 0.5*log_eta + (5.0/6.0)*log(mt) - (log(LAL_PC_SI)+log(1.0e+6)));
  /* leaving out the following term makes freqDomain template scaling match that of "XLALREAL8TimeFreqFFT()" output: */
  /* ampliConst /= IFOdata->timeData->deltaT; */
  plusCoef  = (-0.5*(1.0+pow(cos(iota),2.0)));
  crossCoef = cos(iota);//was   crossCoef = (-1.0*cos(iota));, change iota to -iota+Pi to match HW injection definitions.
  dataStart = XLALGPSGetREAL8(&(IFOdata->timeData->epoch));
  twopitc = LAL_TWOPI * (tc - dataStart);
  a[0] =  exp(log(3.0/128.0) - (5.0/3.0)*log_q - log_eta);
  a[1] =  exp(log(3715.0/84.0+55.0*eta) - log(384.0) - log_eta - log_q);
  a[2] = -exp(log(48.0*LAL_PI/128.0) - (2.0/3.0)*log_q - log_eta);
  a[3] =  exp(log(3.0/128.0) - log_eta - (1.0/3.0)*log_q
              + log(15293365.0/508032.0 + ((27145.0/504.0) + (3085.0/72.0)*eta)*eta));
  a[4] =  exp(log(LAL_PI/128.0)-log_eta+log(38645.0/252.0+5.0*eta));
  NDeltaT = ((double) IFOdata->timeData->data->length) * IFOdata->timeData->deltaT;

  /* Check if fLow is a model parameter, otherwise use data structure definition */
  if(LALInferenceCheckVariable(IFOdata->modelParams, "fLow"))
    fLow = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "fLow");
  else
    fLow = IFOdata->fLow; 
  lower = ceil(fLow * NDeltaT);
  upper = floor(IFOdata->fHigh * NDeltaT);
  /* loop over frequency bins: */
  for (i=0; i<IFOdata->freqModelhPlus->data->length; ++i){
    if ((i > upper) || (i < lower)) /* (no computations outside freq. range) */
      plusRe = plusIm = crossRe = crossIm = 0.0;
    else {
      f   = ((double)i) / NDeltaT;
      f01 = pow(f, -1.0/6.0);       /* = f^-1/6  */
      f02 = f01*f01;                /* = f^-2/6  */
      f04 = f02*f02;                /* = f^-4/6  */
      f06 = f04*f02;                /* = f^-6/6  */
      f07 = f06*f01;                /* = f^-7/6  */
      f10 = f06*f04;                /* = f^-10/6 */
      Psi = a[0]*f10 + a[1]*f06 + a[2]*f04 + a[3]*f02;
      if (PNOrder > 2.0)  /*  5th coefficient ignored for 2.0 PN order  */
        Psi += a[4]*log(f); 
      phaseArg = Psi + twopitc*f + phi;
      plusRe  =  ampliConst * f07 * cos(phaseArg);
      plusIm  =  ampliConst * f07 * (-sin(phaseArg));
      crossRe =  -1.0*plusIm * crossCoef;
      crossIm =  plusRe * crossCoef;
      plusRe  *= plusCoef;
      plusIm  *= plusCoef;
    }
    /* copy f'domain waveform over to IFOdata: */
    IFOdata->freqModelhPlus->data->data[i].real_FIXME  = plusRe;
    IFOdata->freqModelhPlus->data->data[i].imag_FIXME  = plusIm;
    IFOdata->freqModelhCross->data->data[i].real_FIXME = crossRe;
    IFOdata->freqModelhCross->data->data[i].imag_FIXME = crossIm;
  }
  IFOdata->modelDomain = LAL_SIM_DOMAIN_FREQUENCY;
  return;
}



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
    IFOdata->freqModelhPlus->data->data[i].real_FIXME  = 0.0;
    IFOdata->freqModelhPlus->data->data[i].imag_FIXME  = 0.0;
    IFOdata->freqModelhCross->data->data[i].real_FIXME = 0.0;
    IFOdata->freqModelhCross->data->data[i].imag_FIXME = 0.0;
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

static REAL8 fLow2fStart(REAL8 fLow, INT4 ampOrder)
/*  Compute the minimum frequency for waveform generation */
/*  using amplitude orders above Newtonian.  The waveform */
/*  generator turns on all orders at the orbital          */
/*  associated with fMin, so information from higher      */
/*  orders is not included at fLow unless fMin is         */
/*  sufficiently low.                                     */
{
    REAL8 fStart;
    fStart = fLow * 2./(ampOrder+2);
    return fStart;
}

void LALInferenceTemplatePSTRD(LALInferenceIFOData *IFOdata)

/** Template function for PhenSpinTaylorRingDown waveforms. 
    THIS HAS NOT BEEN TESTED! */
{
  static LALStatus status;
  memset(&status,0,sizeof(LALStatus));
  InspiralTemplate template;
  memset(&template,0,sizeof(InspiralTemplate));
  UINT4 idx=0;
	
  /* spin variables still need to be initialised */
  double a_spin1=0.		;
  double theta_spin1=0.	;
  double phi_spin1=0.	;
	
  double a_spin2=0.	;
  double theta_spin2=0.	;
  double phi_spin2=0.	;
	
  /* spin variables still need to be initialised */	
	
  /* spin variables still need to be initialised */
  if (LALInferenceCheckVariable(IFOdata->modelParams, "a_spin1")){		
    a_spin1 = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "a_spin1");
  }
	
  if (LALInferenceCheckVariable(IFOdata->modelParams, "theta_spin1")){
    theta_spin1	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "theta_spin1");
  }
	
  if (LALInferenceCheckVariable(IFOdata->modelParams, "phi_spin1")){
    phi_spin1= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi_spin1");
  }
	
  if (LALInferenceCheckVariable(IFOdata->modelParams, "a_spin2")){		
    a_spin2 = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "a_spin2");
  }
	
  if (LALInferenceCheckVariable(IFOdata->modelParams, "theta_spin2")){
    theta_spin2	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "theta_spin2");
  }
	
  if (LALInferenceCheckVariable(IFOdata->modelParams, "phi_spin2")){
    phi_spin2= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi_spin2");
  }
	
  //double distance = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams,"logdistance");
  //template.distance = exp(distance)*LAL_PC_SI*1.e6;  

  /* spin variables still need to be initialised */
	
  //double mc       = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "chirpmass");
  double logmc = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "logmc");
  double mc = exp(logmc);
  double phi      = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phase");       /* here: startPhase !! */
  double iota     = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "inclination");

  double eta;	
  if (LALInferenceCheckVariable(IFOdata->modelParams,"asym_massratio")) {
    double q = *(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"asym_massratio");
    q2eta(q, &eta);
  }
  else
    eta = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "massratio");
	
  REAL8 mtot=mc/pow(eta,3./5.);	
	
  /* fill the template structure */
  //double distance = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams,"logdistance");
  template.spin1[0]=a_spin1*sin(theta_spin1)*cos(phi_spin1);
  template.spin1[1]=a_spin1*sin(theta_spin1)*sin(phi_spin1);
  template.spin1[2]=a_spin1*cos(theta_spin1); 
  template.spin2[0]=a_spin2*sin(theta_spin2)*cos(phi_spin2);
  template.spin2[1]=a_spin2*sin(theta_spin2)*sin(phi_spin2);
  template.spin2[2]=a_spin2*cos(theta_spin2);
  template.totalMass = mtot;
  template.eta = eta;
  template.massChoice = totalMassAndEta;
  template.tSampling = 1./IFOdata->timeData->deltaT;
  template.fCutoff = 0.5/IFOdata->timeData->deltaT-1.0;
  template.nStartPad = 0;
  template.nEndPad =0;
  template.startPhase = phi;
  template.startTime = 0.0;
  template.ieta = 1;
  template.inclination=iota;
  //template.distance = exp(distance)*LAL_PC_SI*1.e6;
  template.distance = LAL_PC_SI*1.e6;
  int order = *(INT4*) LALInferenceGetVariable(IFOdata->modelParams, "LAL_PNORDER");
  template.order= (LALPNOrder) order; //check order is set correctly
  if (LALInferenceCheckVariable(IFOdata->modelParams, "LAL_APPROXIMANT")){
    template.approximant = *(Approximant*) LALInferenceGetVariable(IFOdata->modelParams, "LAL_APPROXIMANT");
    if(template.approximant!=PhenSpinTaylorRD) {
      XLALPrintError("Error, LALInferenceTemplatePSTRD can only use PhenSpinTaylorRD approximant!");
      XLAL_ERROR_VOID(XLAL_EDATA);
    }
  }
	
  /* Check if fLow is a model parameter, otherwise use data structure definition */
  if(LALInferenceCheckVariable(IFOdata->modelParams, "fLow"))
    template.fLower = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "fLow");
  else
    template.fLower = IFOdata->fLow;	

  template.next = NULL;
  template.fine = NULL;
  int UNUSED errnum;
  XLAL_TRY(LALInspiralParameterCalc(&status,&template),errnum);
	
  REAL4Vector *hPlus = XLALCreateREAL4Vector(IFOdata->timeModelhPlus->data->length);
  REAL4Vector *hCross = XLALCreateREAL4Vector(IFOdata->timeModelhCross->data->length);
	
  XLAL_TRY(LALPSpinInspiralRDTemplates(&status,hPlus,hCross,&template),errnum);

  //REAL4 WinNorm = sqrt(IFOdata->window->sumofsquares/IFOdata->window->data->length);
  for(idx=0;idx<hPlus->length;idx++) IFOdata->timeModelhPlus->data->data[idx]= (REAL8)hPlus->data[idx];
  for(idx=0;idx<hCross->length;idx++) IFOdata->timeModelhCross->data->data[idx]= (REAL8)hCross->data[idx];
  //for(idx=0;idx<hPlus->length;idx++) IFOdata->timeModelhPlus->data->data[idx]*=IFOdata->window->data->data[idx]/WinNorm;
  //for(idx=0;idx<hCross->length;idx++) IFOdata->timeModelhCross->data->data[idx]*=IFOdata->window->data->data[idx]/WinNorm;

  XLALDestroyREAL4Vector(hPlus);
  XLALDestroyREAL4Vector(hCross);

  //executeFT(LALIFOData *IFOdata); //for phenspin we need to transform each of the states separately so i think you can do it with this function, but can you check just incase

  //XLALREAL8TimeFreqFFT(IFOdata->freqModelhPlus, IFOdata->timeModelhPlus, IFOdata->timeToFreqFFTPlan);
  //XLALREAL8TimeFreqFFT(IFOdata->freqModelhCross, IFOdata->timeModelhCross, IFOdata->timeToFreqFFTPlan);
  //for(idx=0;idx<hPlus->length;idx++) fprintf(stderr,"%12.6e\t %12.6ei\n",IFOdata->freqModelhCross->data->data[idx].re, IFOdata->freqModelhCross->data->data[idx].im);	
  //IFOdata->modelDomain = LAL_SIM_DOMAIN_FREQUENCY;

  /*	for(idx=0;idx<IFOdata->timeModelhPlus->data->data[idx];idx++){
	IFOdata->freqModelhPlus->data->data[idx].re*=IFOdata->timeData->deltaT;
	IFOdata->freqModelhPlus->data->data[idx].im*=IFOdata->timeData->deltaT;
	IFOdata->freqModelhCross->data->data[idx].re*=IFOdata->timeData->deltaT;
	IFOdata->freqModelhCross->data->data[idx].im*=IFOdata->timeData->deltaT;
	}
  */		
  double tc       = IFOdata->epoch.gpsSeconds + 1.e-9*IFOdata->epoch.gpsNanoSeconds + template.tC;
  LALInferenceSetVariable(IFOdata->modelParams, "time", &tc);

	
  return;
}


void LALInferenceTemplateLAL(LALInferenceIFOData *IFOdata)
/*************************************************************************************************/
/* Wrapper function to call LAL functions for waveform generation.                               */
/* Will always return frequency-domain templates (numerically FT'ed                              */
/* in case the LAL function returns time-domain).                                                */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`IFOdata->modelParams') parameters are:                                             */
/*   - "chirpmass"        (REAL8,units of solar masses)                                          */
/*   - "massratio"        (symmetric mass ratio:  0 < eta <= 0.25, REAL8) <or asym_massratio>    */
/*   - "asym_massratio"   (asymmetric mass ratio:  0 < q <= 1.0, REAL8)   <or massratio>         */
/*   - "phase"            (here: 'startPhase', not coalescence phase; REAL8, radians)            */
/*   - "time"             (coalescence time, or equivalent/analog/similar; REAL8, GPS sec.)      */
/*   - "inclination"      (inclination angle, REAL8, radians)                                    */
/*   - "LAL_APPROXIMANT"  (INT4 value corresponding to `enum approximant' definition             */
/*                         in `LALInspiral.h'.                                                   */
/*                         Templates that (seem to) work by now are:                             */
/*                         TaylorF2, TaylorT1, TaylorT2, TaylorT3, BCV, IMRPhenomA, EOB, EOBNR)  */
/*   - "LAL_PNORDER"      (INT4 value corresponding to `enum LALPNOrder' definition              */
/*                         in `LALInspiral.h'.)                                                  */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* 'problematic' templates are:                                                                  */
/*  - Taylor F1 :  returns with error ("Maximum iterations exceeded")                            */
/*  - Taylor T2 :  fails at low mass (error: "Attempting to write beyond the end of vector")     */
/*  - Taylor T4 :  funny scaling, ~ 16 orders of magnitude too large                             */
/*  - EOBNR     :  fails for low masses (e.g.: mc=3, eta=0.24)                                   */
/*  - BCV       :  amplitude is "arbitrary" (as stated in documentation)                         */
/*                                                                                               */
/*************************************************************************************************/
{
  static LALStatus status;
  memset(&status,0,sizeof(status));

  static InspiralTemplate params;
  static REAL4Vector *LALSignal=NULL;
  UINT4 n;
  unsigned long i,j, jmax=0;
  double pj, pmax, pleft, pright;
  double m1,m2,mc,eta,q;
  /* Prefer m1 and m2 if available (i.e. for injection template) */
  if(LALInferenceCheckVariable(IFOdata->modelParams,"mass1")&&LALInferenceCheckVariable(IFOdata->modelParams,"mass2"))
    {
      m1=*(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"mass1");
      m2=*(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"mass2");
      eta=m1*m2/((m1+m2)*(m1+m2));
      mc=pow(eta , 0.6)*(m1+m2);
    }
  else
    {
      if (LALInferenceCheckVariable(IFOdata->modelParams,"asym_massratio")) {
        q = *(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"asym_massratio");
        q2eta(q, &eta);
      }
      else
        eta = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "massratio");
      mc       = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "chirpmass");
      mc2masses(mc, eta, &m1, &m2);
    }
  double phi      = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phase");       /* here: startPhase !! */
  double tc       = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "time");
  double iota     = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "inclination");
  double spin1    = 0.0;
  double spin2    = 0.0;
  /* Just two spins specified - assume they are the z-components */
  if (LALInferenceCheckVariable(IFOdata->modelParams, "spin1")){
    spin1 =  *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "spin1");
    params.spin1[2]=spin1;
  }
  if (LALInferenceCheckVariable(IFOdata->modelParams, "spin2")) {
    spin2 =  *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "spin2");
    params.spin2[2]=spin2;
  }
  
  if(LALInferenceCheckVariable(IFOdata->modelParams,"a_spin1")&&LALInferenceCheckVariable(IFOdata->modelParams,"theta_spin1")&&LALInferenceCheckVariable(IFOdata->modelParams,"phi_spin1"))
    {
      REAL8 theta=*(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"theta_spin1");
      REAL8 phi_spin=*(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"phi_spin1");
      REAL8 a=*(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"a_spin1");
      params.spin1[0]=a*sin(theta)*cos(phi_spin);
      params.spin1[1]=a*sin(theta)*sin(phi_spin);
      params.spin1[2]=a*cos(theta);
    }
  if(LALInferenceCheckVariable(IFOdata->modelParams,"a_spin2")&&LALInferenceCheckVariable(IFOdata->modelParams,"theta_spin2")&&LALInferenceCheckVariable(IFOdata->modelParams,"phi_spin2"))
    {
      REAL8 theta=*(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"theta_spin2");
      REAL8 phi_spin=*(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"phi_spin2");
      REAL8 a=*(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"a_spin2");
      params.spin2[0]=a*sin(theta)*cos(phi_spin);
      params.spin2[1]=a*sin(theta)*sin(phi_spin);
      params.spin2[2]=a*cos(theta);
    }


  int approximant=0, order=0;
  int FDomain;    /* (denotes domain of the _LAL_ template!) */
  double chirptime, deltaT;
  double plusCoef  = -0.5 * (1.0 + pow(cos(iota),2.0));
  double crossCoef = cos(iota);//was   crossCoef = (-1.0*cos(iota));, change iota to -iota+Pi to match HW injection definitions.
  double instant;
  int forceTimeLocation;
  double twopit, f, deltaF, re, im, templateReal, templateImag;
  double fLow;

  if (LALInferenceCheckVariable(IFOdata->modelParams, "LAL_APPROXIMANT"))
    approximant = *(INT4*) LALInferenceGetVariable(IFOdata->modelParams, "LAL_APPROXIMANT");
  else {
    XLALPrintError(" ERROR in templateLAL(): (INT4) \"LAL_APPROXIMANT\" parameter not provided!\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }

  if (LALInferenceCheckVariable(IFOdata->modelParams, "LAL_PNORDER"))
    order = *(INT4*) LALInferenceGetVariable(IFOdata->modelParams, "LAL_PNORDER");
  else {
    XLALPrintError(" ERROR in templateLAL(): (INT4) \"LAL_PNORDER\" parameter not provided!\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }

  /*fprintf(stdout, " templateLAL() - approximant = %d,  PN order = %d\n", approximant, order);*/

  /* little consistency check (otherwise no output without warning): */
  if (((approximant==EOBNR) || (approximant==EOB)) && (order!=LAL_PNORDER_PSEUDO_FOUR)) {
    XLALPrintError(" ERROR in templateLAL(): \"EOB\" and \"EOBNR\" templates require \"LAL_PNORDER_PSEUDO_FOUR\" PN order!\n");  
    XLAL_ERROR_VOID(XLAL_EDATA);
  }
    
  if (IFOdata->timeData==NULL) {
    XLALPrintError(" ERROR in templateLAL(): encountered unallocated 'timeData'.\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }
  if ((IFOdata->freqModelhPlus==NULL) || (IFOdata->freqModelhCross==NULL)) {
    XLALPrintError(" ERROR in templateLAL(): encountered unallocated 'freqModelhPlus/-Cross'.\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }
  deltaT = IFOdata->timeData->deltaT;

  /* Check if fLow is a model parameter, otherwise use data structure definition */
  if(LALInferenceCheckVariable(IFOdata->modelParams, "fLow"))
    fLow = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "fLow");
  else
    fLow = IFOdata->fLow;	
  
  params.OmegaS      = 0.0;     /* (?) */
  params.Theta       = 0.0;     /* (?) */
  /* params.Zeta2    = 0.0; */  /* (?) */
  params.ieta        = 1; 
  params.nStartPad   = 0;
  params.nEndPad     = 0;
  params.massChoice  = m1Andm2;
  params.approximant = (Approximant) approximant;  /*  TaylorT1, ...   */
  params.order       = (LALPNOrder) order;        /*  Newtonian, ...  */
  params.fLower      = fLow * 0.9;
  params.fCutoff     = (IFOdata->freqData->data->length-1) * IFOdata->freqData->deltaF;  /* (Nyquist freq.) */
  params.tSampling   = 1.0 / deltaT;
  params.startTime   = 0.0;

  /* actual inspiral parameters: */
  params.mass1       = m1;
  params.mass2       = m2;
  params.startPhase  = phi;
  if ((params.approximant == EOB) 
      || (params.approximant == EOBNR)
      || (params.approximant == TaylorT3)
      || (params.approximant == IMRPhenomA)
      || (params.approximant == TaylorF2RedSpin))
    params.distance  = LAL_PC_SI * 1.0e6;        /* distance (1 Mpc) in units of metres */
  else if ((params.approximant == TaylorT1)
           || (params.approximant == TaylorT2)
           || (params.approximant == PadeT1)
           || (params.approximant == TaylorF1)
           || (params.approximant == TaylorF2)
           || (params.approximant == PadeF1)
           || (params.approximant == BCV))
    params.distance  = 1.0;                                          /* distance in Mpc */
  else                                                     
    params.distance  = LAL_PC_SI * 1.0e6 / ((double) LAL_C_SI);  /* distance in seconds */

  /* ensure proper "fCutoff" setting: */
  if (params.fCutoff >= 0.5*params.tSampling)
    params.fCutoff = 0.5*params.tSampling - 0.5*IFOdata->freqData->deltaF;
  if (! (params.tSampling > 2.0*params.fCutoff)){
    fprintf(stderr," ERROR in templateLAL(): 'LALInspiralSetup()' (called within 'LALInspiralWavelength()')\n");
    fprintf(stderr,"                         requires (tSampling > 2 x fCutoff) !!\n");
    fprintf(stderr," (settings are:  tSampling = %f s,  fCutoff = %f Hz)  \n", params.tSampling, params.fCutoff);
    exit(1);
  }

  /* ensure compatible sampling rate: */
  if ((params.approximant == EOBNR)
      && (fmod(log((double)params.tSampling)/log(2.0),1.0) != 0.0)) {
    fprintf(stderr, " ERROR in templateLAL(): \"EOBNR\" templates require power-of-two sampling rates!\n");
    fprintf(stderr, "                         (params.tSampling = %f Hz)\n", params.tSampling);
    exit(1);
  }

  /* compute other elements of `params', check out the `.tC' value, */
  /* shift the start time to match the coalescence time,            */
  /* and eventually re-do parameter calculations:                   */
  /* Reset errno. */

  LALInspiralParameterCalc(&status, &params);
  chirptime = params.tC;
  if ((params.approximant != TaylorF2) && (params.approximant != TaylorF2RedSpin) && (params.approximant != BCV)) {
    params.startTime = (tc - XLALGPSGetREAL8(&IFOdata->timeData->epoch)) - chirptime;
    LALInspiralParameterCalc(&status, &params); /* (re-calculation necessary? probably not...) */
  }

  if (params.approximant == TaylorF2 || params.approximant == TaylorF2RedSpin) {	
    expnCoeffs ak;
    expnFunc expnFunction;
    memset(&ak,0,sizeof(expnCoeffs));
    /* Calculate the time of ISCO (v = 6^(-1/2) ) */
    LALInspiralSetup(&status,&ak,&params);
    LALInspiralChooseModel(&status,&expnFunction,&ak,&params);
    chirptime=ak.tn;
  }

  /* compute "params.signalAmplitude" slot: */
  LALInspiralRestrictedAmplitude(&status, &params);

  /* figure out inspiral length & set `n': */
  /* LALInspiralWaveLength(&status, &n, params); */
  n = IFOdata->timeData->data->length;

  /* domain of LAL template as returned by LAL function: */
  FDomain = ((params.approximant == TaylorF1)
             || (params.approximant == TaylorF2)
             || (params.approximant == TaylorF2RedSpin)
             || (params.approximant == PadeF1)
             || (params.approximant == BCV));
  if (FDomain && (n % 2 != 0)){
    fprintf(stderr, " ERROR in templateLAL(): frequency-domain LAL waveforms require even number of samples!\n");
    fprintf(stderr, "                         (N = IFOdata->timeData->data->length = %d)\n", n);
    exit(1);
  }

  /* allocate (temporary) waveform vector: */
  LALCreateVector(&status, &LALSignal, n);
  
  /*--  ACTUAL WAVEFORM COMPUTATION:  --*/
  /* Catch any previous errors */
  if (status.statusCode != 0) {
    fprintf(stderr, " ERROR in templateLAL(): encountered non-zero status code.\n");
    fprintf(stderr, " Template parameters:\n");
    LALInferencePrintVariables(IFOdata->modelParams);
    fprintf(stderr, " LAL Status:\n");
    REPORTSTATUS(&status);
    exit(1);
  }

  /* Check for the Integration overflow error and work around it */
  if(XLALGetBaseErrno() == XLAL_EMAXITER )
    {
      XLALPrintError("Template generation failed!");
      if(LALSignal)     LALDestroyVector(&status, &LALSignal);                                                                                                                  
      XLAL_ERROR_VOID(XLAL_FAILURE);
    }

  memset(LALSignal->data,0,LALSignal->length*sizeof(LALSignal->data[0]));

  
  if (params.fCutoff >= 0.5*params.tSampling)
    params.fCutoff = 0.5*params.tSampling - 0.5*IFOdata->freqData->deltaF; //should not be needed.
  
  // lal_errhandler = LAL_ERR_RTRN;
  // REPORTSTATUS(&status); 
  LALInspiralWave(&status, LALSignal, &params);
  // REPORTSTATUS(&status); 
  // lal_errhandler = LAL_ERR_DFLT; 
  if (status.statusCode != 0) {
    fprintf(stderr, "\n ERROR in templateLAL(): \"LALInspiralWave()\" call returned with non-zero status.\n");
    fprintf(stderr, " Template parameters:\n");
    LALInferencePrintVariables(IFOdata->modelParams);
    fprintf(stderr, " LAL Status:\n");
    REPORTSTATUS(&status);
    exit(1);
  }

  if (! FDomain) {   /*  (LAL function returns TIME-DOMAIN template)       */
    IFOdata->modelDomain = LAL_SIM_DOMAIN_TIME;

    /* copy over, normalise: */
    for (i=0; i<n; ++i) {
      IFOdata->timeModelhPlus->data->data[i]  = LALSignal->data[i];
      IFOdata->timeModelhCross->data->data[i] = 0.0;  /* (no cross waveform) */
    }
    LALDestroyVector(&status, &LALSignal);
    /* apply window & execute FT of plus component: */
    if (IFOdata->window==NULL) {
      XLALPrintError(" ERROR in templateLAL(): ran into uninitialized 'IFOdata->window'.\n");
      XLAL_ERROR_VOID(XLAL_EFAULT);
    }
    XLALDDVectorMultiply(IFOdata->timeModelhPlus->data, IFOdata->timeModelhPlus->data, IFOdata->window->data);
    if (IFOdata->timeToFreqFFTPlan==NULL) {
      XLALPrintError(" ERROR in templateLAL(): ran into uninitialized 'IFOdata->timeToFreqFFTPlan'.\n");
      XLAL_ERROR_VOID(XLAL_EFAULT);
    }
    XLALREAL8TimeFreqFFT(IFOdata->freqModelhPlus, IFOdata->timeModelhPlus, IFOdata->timeToFreqFFTPlan);
    /* Normalise by RMS of window (same as injections and data) */
    REAL8 WinNorm=sqrt(IFOdata->window->sumofsquares/IFOdata->window->data->length);
    for(i=0;i<IFOdata->freqModelhPlus->data->length;i++) {
      IFOdata->freqModelhPlus->data->data[i].real_FIXME/=WinNorm;
      IFOdata->freqModelhPlus->data->data[i].imag_FIXME/=WinNorm;
    }
  }  
  else
    {             /*  (LAL function returns FREQUENCY-DOMAIN template)  */
      IFOdata->modelDomain = LAL_SIM_DOMAIN_FREQUENCY;

      /* copy over: */
      IFOdata->freqModelhPlus->data->data[0].real_FIXME = ((REAL8) LALSignal->data[0]);
      IFOdata->freqModelhPlus->data->data[0].imag_FIXME = 0.0;
      for (i=1; i<IFOdata->freqModelhPlus->data->length-1; ++i) {
	IFOdata->freqModelhPlus->data->data[i].real_FIXME = ((REAL8) LALSignal->data[i]);
	IFOdata->freqModelhPlus->data->data[i].imag_FIXME = ((REAL8) LALSignal->data[n-i]);
      }
      IFOdata->freqModelhPlus->data->data[IFOdata->freqModelhPlus->data->length-1].real_FIXME = LALSignal->data[IFOdata->freqModelhPlus->data->length-1];
      IFOdata->freqModelhPlus->data->data[IFOdata->freqModelhPlus->data->length-1].imag_FIXME = 0.0;
      LALDestroyVector(&status, &LALSignal);
      /* nomalise (apply same scaling as in XLALREAL8TimeFreqFFT()") : */
      for (i=0; i<IFOdata->freqModelhPlus->data->length; ++i) {
	IFOdata->freqModelhPlus->data->data[i].real_FIXME *= ((REAL8) n) * deltaT;
	IFOdata->freqModelhPlus->data->data[i].imag_FIXME *= ((REAL8) n) * deltaT;
      }
      if(LALInferenceCheckVariable(IFOdata->modelParams, "ppealpha") && LALInferenceCheckVariable(IFOdata->modelParams, "ppeuppera") &&
	 LALInferenceCheckVariable(IFOdata->modelParams, "ppelowera") && LALInferenceCheckVariable(IFOdata->modelParams, "ppebeta") &&
	 LALInferenceCheckVariable(IFOdata->modelParams, "ppeupperb") && LALInferenceCheckVariable(IFOdata->modelParams, "ppelowerb")){
      
	REAL8 alpha, A, a, beta, B, b, ppE_amp, ppE_phase, cos_ppE_phase, sin_ppE_phase;
	alpha =  *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "ppealpha");
	A     =  *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "ppeuppera");
	a     =  *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "ppelowera");
	beta  =  *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "ppebeta");
	B     =  *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "ppeupperb");
	b     =  *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "ppelowerb");
      
	for (i=0; i<IFOdata->freqModelhPlus->data->length; ++i) {
	  ppE_amp = 1.0+alpha*pow(4.0*eta,A)*pow(LAL_PI*mc*(fLow*0.9 + ((REAL8) i)*IFOdata->freqData->deltaF),a);
	  ppE_phase = beta*pow(4.0*eta,B)*pow(LAL_PI*mc*(fLow*0.9 + ((REAL8) i)*IFOdata->freqData->deltaF),b);
	  cos_ppE_phase = cos(ppE_phase);
	  sin_ppE_phase = sin(ppE_phase);
      
	  IFOdata->freqModelhPlus->data->data[i].real_FIXME = (ppE_amp)*(creal(IFOdata->freqModelhPlus->data->data[i])*cos_ppE_phase-cimag(IFOdata->freqModelhPlus->data->data[i])*sin_ppE_phase);
	  IFOdata->freqModelhPlus->data->data[i].imag_FIXME = (ppE_amp)*(creal(IFOdata->freqModelhPlus->data->data[i])*sin_ppE_phase+cimag(IFOdata->freqModelhPlus->data->data[i])*cos_ppE_phase);
	}
      }
    }

  /* (now frequency-domain plus-waveform has been computed, either directly or via FFT)   */

  /*  cross waveform is "i x plus" :  */
  for (i=1; i<IFOdata->freqModelhCross->data->length-1; ++i) {
    IFOdata->freqModelhCross->data->data[i].real_FIXME = -cimag(IFOdata->freqModelhPlus->data->data[i]);
    IFOdata->freqModelhCross->data->data[i].imag_FIXME = creal(IFOdata->freqModelhPlus->data->data[i]);
    // consider inclination angle's effect:
    IFOdata->freqModelhPlus->data->data[i].real_FIXME  *= plusCoef;
    IFOdata->freqModelhPlus->data->data[i].imag_FIXME  *= plusCoef;
    IFOdata->freqModelhCross->data->data[i].real_FIXME *= crossCoef;
    IFOdata->freqModelhCross->data->data[i].imag_FIXME *= crossCoef;
  }

  /*
   * NOTE: the dirty trick here is to assume the LAL waveform to constitute
   *       the cosine chirp and then derive the corresponding sine chirp 
   *       as the orthogonal ("i x cosinechirp") waveform.
   *       In general they should not necessarily be only related 
   *       by a mere phase shift though...
   */

  /* Now...template is not (necessarily) located at specified coalescence time  */
  /* and/or we don't know even where it actually is located...                  */
  /* Figure out time location corresponding to template just computed:          */

  /* default: assume template to have correctly considered              */
  /* the supplied "params.tC" value                                     */
  /* (Roughly OK for "TaylorF1" (?), "TaylorT1", "TaylorT3", "EOB",     */
  /* "EOBNR", and "PadeT1".                                             */
  /* May still by off by tens/hundreds of milliseconds.):               */
  instant = tc;

  /* Signal simply evolved from start of template on,         */
  /* for approximately "chirptime" seconds:                   */
  if ((params.approximant == TaylorT2) 
      || (params.approximant == TaylorF2)
      || (params.approximant == TaylorF2RedSpin))
    instant = XLALGPSGetREAL8(&IFOdata->timeData->epoch) + chirptime;

  /* Coalescence happens at very end of signal template:      */
  else if (params.approximant == BCV) 
    instant = XLALGPSGetREAL8(&IFOdata->timeData->epoch) + ((double) IFOdata->timeData->data->length) * deltaT;

  /* No idea where signal lies; brute-force search for amplitude peak: */
  /* (this is time-comsuming and should be avoided where possible!!)   */
  else  if (params.approximant == IMRPhenomA) {
    /* Inv-FT back to time domain: */
    /* (admittedly, this extra FT is time-consuming not elegant...  */
    /* but might be ok given that once generated, templates may be  */
    /* re-used at different timeshifts/skylocations/etc.)           */
    LALInferenceExecuteInvFT(IFOdata);
    /* find amplitude peak & two neighbouring bins: */
    pmax = 0.0;
    for (j=0; j<IFOdata->timeModelhPlus->data->length; ++j) {
      pj = IFOdata->timeModelhPlus->data->data[j] * IFOdata->timeModelhPlus->data->data[j]
	+ IFOdata->timeModelhCross->data->data[j] * IFOdata->timeModelhCross->data->data[j];
      if (pj > pmax){
        pmax = pj;
        jmax = j;
      }
    }
    j = (jmax>0) ? jmax-1 : IFOdata->timeModelhPlus->data->length-1;
    pleft = sqrt(IFOdata->timeModelhPlus->data->data[j] * IFOdata->timeModelhPlus->data->data[j]
                 + IFOdata->timeModelhCross->data->data[j] * IFOdata->timeModelhCross->data->data[j]);
    j = (jmax<IFOdata->timeModelhPlus->data->length-1) ? jmax+1 : 0;
    pright = sqrt(IFOdata->timeModelhPlus->data->data[j] * IFOdata->timeModelhPlus->data->data[j]
                  + IFOdata->timeModelhCross->data->data[j] * IFOdata->timeModelhCross->data->data[j]);
    pmax = sqrt(pmax);
    /* do some ad-hoc corrections to ensure actually having a peak: */
    if (!((pleft<pmax) || (pright<pmax)))
      pleft = pright = pmax - 1.0;
    else if (!(pleft<pmax)) pleft = 0.5*(pmax+pright);
    else if (!(pright<pmax)) pright = 0.5*(pmax+pleft);
    /*  do a quadratic interpolation                        */
    /*  to determine peak location to sub-deltaT accuracy:  */
    instant = (pleft-pright) / (2.0*pleft-4.0*pmax+2.0*pright);
    instant = (XLALGPSGetREAL8(&IFOdata->timeData->epoch) + jmax*deltaT) + instant*deltaT;
    /* fprintf(stdout, " interpolated location: %.8f GPS sec.\n", instant); */
  }

  /* now either time-shift template or just store the time value: */
  /* (time-shifting should not be necessary in general,           */
  /* but may be neat to have for de-bugging etc.)                 */
  forceTimeLocation = 0;  /* default: zero! */
  if (instant != tc) {
    if (forceTimeLocation) { /* time-shift the frequency-domain template: */
      twopit = LAL_TWOPI * (tc - instant);
      deltaF = 1.0 / (((double)IFOdata->timeData->data->length) * deltaT);
      for (i=1; i<IFOdata->freqModelhPlus->data->length; ++i){
        f = ((double) i) * deltaF;
        /* real & imag parts of  exp(-2*pi*i*f*deltaT): */
        re = cos(twopit * f);
        im = - sin(twopit * f);
        templateReal = creal(IFOdata->freqModelhPlus->data->data[i]);
        templateImag = cimag(IFOdata->freqModelhPlus->data->data[i]);
        IFOdata->freqModelhPlus->data->data[i].real_FIXME = templateReal*re - templateImag*im;
        IFOdata->freqModelhPlus->data->data[i].imag_FIXME = templateReal*im + templateImag*re;
        templateReal = creal(IFOdata->freqModelhCross->data->data[i]);
        templateImag = cimag(IFOdata->freqModelhCross->data->data[i]);
        IFOdata->freqModelhCross->data->data[i].real_FIXME = templateReal*re - templateImag*im;
        IFOdata->freqModelhCross->data->data[i].imag_FIXME = templateReal*im + templateImag*re;
      }
    }
    else {
      /* write template (time axis) location in "->modelParams" so that     */
      /* template corresponds to stored parameter values                    */
      /* and other functions may time-shift template to where they want it: */
      LALInferenceSetVariable(IFOdata->modelParams, "time", &instant);
    }
  }

  IFOdata->modelDomain = LAL_SIM_DOMAIN_FREQUENCY;
  return;
}



void LALInferenceTemplate3525TD(LALInferenceIFOData *IFOdata)
/*****************************************************************/
/* 3.5PN phase / 2.5PN amplitude time-domain inspiral templates  */
/* following                                                     */
/*   Blanchet et al. 2001   gr-qc/0104084                        */
/*   Blanchet at al. 2002   PRD 65(6):061501    gr-qc/0105099    */
/*   Blanchet at al. 2005   PRD 71(12):129902                    */
/*   Arun et al. 2004       CQG 21(15):3771                      */
/*   Arun et al. 2004       CQG 22(14):3115                      */
/*   Blanchet et al. 2004   PRL 93(9):091101                     */
/* This is basically the implementation that was also used in    */
/* the "Roever/Meyer/Guidi/Vicere/Christensen (2007)" paper      */
/* (CQG 24(19):S607).                                            */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Formula numbers (x.xx) refer to the 2001 Blanchet paper,      */
/* numbers (xx) refer to the more recent 2002 paper.             */
/* Numbers referring to Arun et al (2004) are explicitly marked. */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *********************************************/
/* Required (`IFOdata->modelParams') parameters are:                                                         */
/*   - "chirpmass"        (REAL8, chirp mass, in units of solar masses)                                      */
/*   - "massratio"        (REAL8, symmetric mass ratio:  0 < eta <= 0.25, dimensionless) <or asym_massratio> */
/*   - "asym_massratio"   (REAL8, asymmetric mass ratio:  0 < q <= 1.0, dimensionless)   <or massratio>      */
/*   - "phase"            (REAL8, coalescence phase, radians)                                                */
/*   - "time"             (REAL8, coalescence time, GPS seconds)                                             */
/*   - "inclination"      (REAL8, inclination angle, radians)                                                */
/*************************************************************************************************************/
{
  double mc    = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "chirpmass");      /* chirp mass m_c, solar masses           */
  double tc    = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "time");           /* coalescence time, GPS sec.             */
  double phase = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phase");          /* coalescence phase, rad                 */
  double eta;                                                                               /* mass ratio eta, dimensionless          */
  if (LALInferenceCheckVariable(IFOdata->modelParams,"asym_massratio")) {
    double q = *(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"asym_massratio");  /* asymmetric mass ratio q, dimensionless */
    q2eta(q, &eta);
  }
  else
    eta = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "massratio");
  double m1, m2;
  mc2masses(mc, eta, &m1, &m2);                   /* (in units of Msun) */
  double mt         = m1 + m2;
  double dmm        = (m2-m1)/mt;                 /*  = (delta m) / mt  (dimensionless) */
  double log_mt     = log(mt) + log(LAL_MSUN_SI); /* (in Kg) */
  double log_eta    = log(eta);
  double eta2       = eta * eta;
  double eta3       = eta2 * eta;
  double log_mu     = log_eta + log_mt;
  double log_omega0 = log(4.0*LAL_PI);
  double log_tau0   = 0.0;  /* = log(1.0) */
  double t, phi, psi;
  double taucoef = 3.0*log((double) LAL_C_SI)-log(5.0)-log(LAL_G_SI) + log_eta - log_mt; /*  (4.17) or (11) */
  double log_tau, tau18, tau28, tau38, tau48, tau58, tau68, tau78;
  double ci  =  cos(*(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "inclination"));
  double ci2 = ci*ci,     ci4 = ci2*ci2,   ci6 = ci4*ci2;
  double si2 = (1.0-ci2), si  = sqrt(si2), si4 = si2*si2, si5 = si4*si;
  double h_plus, h_cross;
  double Hp00, Hp05, Hp10, Hp15, Hp20, Hp25;
  double Hc00, Hc05, Hc10, Hc15, Hc20, Hc25;
  double plus10a  = (1.0/6.0)*((19.0+9.0*ci2-2.0*ci4)-eta*(19.0-11.0*ci2-6.0*ci4));   /* (6.4) */
  double plus10b  = (4.0/3.0)*si2*(1.0+ci2)*(1.0-3.0*eta);
  double plus15a  = ((57.0 + 60.0*ci2-ci4) - 2.0*eta*(49.0-12.0*ci2-ci4));            /* (6.5) */
  double plus15b  = 13.5*((73.0+40.0*ci2-9.0*ci4) - 2.0*eta*(25.0-8.0*ci2-9.0*ci4));
  double plus15c  = 312.5*(1.0-2.0*eta)*si2*(1.0+ci2);
  double plus20a  = (1.0/120.0)*((22.0+396.0*ci2+145.0*ci4-5.0*ci6)                   /* (6.6) */
				 + (5.0/3.0)*eta*(706.0-216.0*ci2-251.0*ci4+15.0*ci6)
				 -5.0*eta2*(98.0-108.0*ci2+7.0*ci4+5.0*ci6));
  double plus20b  = (2.0/15.0)*si2*((59.0+35.0*ci2-8.0*ci4)
				    -(5.0/3.0)*eta*(131.0+59.0*ci2-24.0*ci4)
				    +5.0*eta2*(21.0-3.0*ci2-8.0*ci4));
  double plus20c  = 2.025*(1.0-5.0*eta+5.0*eta2)*si4*(1.0+ci2);
  double plus20d  = (11.0+7.0*ci2+10.0*(5.0+ci2)*LAL_LN2);
  double plus20e  = 27.0*(7.0-10.0*log(1.5));
  double plus25a  = si*dmm*((1771.0/5120.0)-(1667.0/5120.0)*ci2+(217.0/9216.0)*ci4-(1.0/9216.0)*ci6
                            +eta*((681.0/256.0)+(13.0/768.0)*ci2-(35.0/768.0)*ci4+(1.0/2304.0)*ci6)
                            +eta2*(-(3451.0/9216.0)+(673.0/3072.0)*ci2-(5.0/9216.0)*ci4-(1.0/3072.0)*ci6)); /* Arun (5.9) */
  double plus25b  = LAL_PI*((19.0/3.0)+3.0*ci2-(2.0/3.0)*ci4
                            +eta*(-(16.0/3.0)+(14.0/3.0)*ci2+2.0*ci4));
  double plus25c  = si*dmm*((3537.0/1024.0)-(22977.0/5120.0)*ci2-(15309.0/5120.0)*ci4+(729.0/5120.0)*ci6
                            +eta*(-(23829.0/1280.0)+(5529.0/1280.0)*ci2+(7749.0/1280.0)*ci4-(729.0/1280.0)*ci6)
	                    +eta2*((29127.0/5120.0)-(27267.0/5120.0)*ci2-(1647.0/5120.0)*ci4+(2187.0/5120.0)*ci6));
  double plus25d  = (-(16.0/3.0)*LAL_PI*(1.0+ci2)*si2*(1.0-3.0*eta));
  double plus25e  = si*dmm*(-(108125.0/9216.0)+(40625.0/9216.0)*ci2+(83125.0/9216.0)*ci4-(15625.0/9216.0)*ci6
                            +eta*((8125.0/265.0)-(40625.0/2304.0)*ci2-(48125.0/2304.0)*ci4+(15625.0/2304.0)*ci6)
                            +eta2*(-(119375.0/9216.0)+(40625.0/3072.0)*ci2+(44375.0/9216.0)*ci4-(15625.0/3072.0)*ci6));
  double plus25f  = dmm*((117649.0/46080.0)*si5*(1.0+ci2)*(1.0-4.0*eta+3.0*eta2));
  double plus25g  = (-1.8+2.8*ci2+1.4*ci4+eta*(19.2-1.6*ci2-5.6*ci4));
  double plus25h  = si2*(1.0+ci2)*(11.2 - 32.0*LAL_LN2/3.0 - eta*(1193.0/30.0 - 32.0*LAL_LN2));

  double cross10a = (ci/3.0)*((17.0-4.0*ci2)-eta*(13.0-12.0*ci2));                    /* (6.9) */
  double cross10b = (8.0/3.0)*(1.0-3.0*eta)*ci*si2;
  double cross15a = ((63.0-5.0*ci2)-2.0*eta*(23.0-5.0*ci2));                          /* (6.10) */
  double cross15b = 13.5*((67.0-15.0*ci2)-2.0*eta*(19.0-15.0*ci2));
  double cross15c = 312.5*(1.0-2.0*eta)*si2;
  double cross20a = (ci/60.0)*((68.0+226.0*ci2-15.0*ci4)+(5.0/3.0)*eta*(572.0-490.0*ci2+45.0*ci4)
			       -5.0*eta2*(56.0-70.0*ci2+15.0*ci4));                              /* (6.11) */
  double cross20b = (4.0/15.0)*ci*si2*((55.0-12.0*ci2)-(5.0/3.0)*eta*(119.0-36.0*ci2)
				       +5.0*eta2*(17.0-12.0*ci2));
  double cross20c = 4.05*(1.0-5.0*eta+5.0*eta2)*ci*si4;
  double cross20d = 3.0+10*LAL_LN2;
  double cross20e = 9.0*(7.0-10.0*log(1.5));
  double cross25a = 1.2*si2*ci*eta;                                                   /* Arun (5.10) */
  double cross25b = ci*(2.0-4.4*ci2+eta*(-30.8+18.8*ci2));
  double cross25c = ci*si2*((-112.0/5.0 + (64.0/3.0)*LAL_LN2)+eta*(1193.0/15.0 - 64.0*LAL_LN2));
  double cross25d = si*ci*dmm*(-(913.0/7680.0)+(1891.0/11520.0)*ci2-(7.0/4608.0)*ci4
                               +eta*((1165.0/384.0)-(235.0/576.0)*ci2+(7.0/1152.0)*ci4)
                               +eta2*(-(1301.0/4608.0)+(301.0/2304.0)*ci2-(7.0/1536.0)*ci4));
  double cross25e = LAL_PI*ci*((34.0/3.0)-(8.0/3.0)*ci2-eta*((20.0/3.0)-8.0*ci2));
  double cross25f = si*ci*dmm*((12501.0/2560.0)-(12069.0/1260.0)*ci2+(1701.0/2560.0)*ci4
                               +eta*(-(19581.0/640.0)+(7821.0/320.0)*ci2-(1701.0/640.0)*ci4)
                               +eta2*((18903.0/2560.0)-(11403.0/1280.0)*ci2+(5103.0/2560.0)*ci4));
  double cross25g = si2*ci*(-((32.0/3.0)*LAL_PI)*(1.0-3.0*eta));
  double cross25h = dmm*si*ci*(-(101875.0/4608.0)+(6875.0/256.0)*ci2-(21875.0/4608.0)*ci4
                               +eta*((66875.0/1152.0)-(44375.0/576.0)*ci2+(21875.0/1152.0)*ci4)
                               +eta2*(-(100625.0/4608.0)+(83125.0/2304.0)*ci2-(21875.0/1536.0)*ci4));
  double cross25i = dmm*si5*ci*((117649.0/23040.0)*(1.0-4.0*eta+3.0*eta2));
  double sin1psi, sin2psi, sin3psi, sin4psi, sin5psi, sin6psi, sin7psi;
  double cos1psi, cos2psi, cos3psi, cos4psi, cos5psi, cos6psi, cos7psi;
  double constfactor = exp(LAL_LN2+log(LAL_G_SI)-2.0*log((double)LAL_C_SI) + log_mu - log(LAL_PC_SI*1.0e6));  
  double x, sqrtx, oldx=0.0;                                                          /* (6.01); distance is 1 Mpc here. */
  double omega, omegacoef=exp(3.0*log((double) LAL_C_SI) - log(LAL_G_SI) - log_mt);   /* = (c^3)/(G*mt) */
  double EulerGamma = 0.57721566490153286; /* Euler constant */
  double xi     = -9871.0/9240.0;          /* Blanchet et al (2004): PRL 93(9):091101 */
  double kappa  = 0.0;                     /* (ibid.)                                 */
  double zeta   = -7.0/33.0;               /* (ibid.)                                 */
  double theta  = xi + 2.0*kappa + zeta;    
  double lambda = -(1987.0/3080);           
  double PI2    = LAL_PI * LAL_PI;
  double xcoef1 =    (743.0/4032.0)   +    (11.0/48.0)    *eta;                       /* (12) */
  double xcoef2 =  (19583.0/254016.0) + (24401.0/193536.0)*eta + (31.0/288.0)*eta2;
  double xcoef3 = -(11891.0/53760.0)  +   (109.0/1920.0)  *eta;
  double xcoef4 = (-10052469856691.0/6008596070400.0 + PI2/6.0 + (107.0/420.0)*EulerGamma)
    + (15335597827.0/3901685760.0 - (451.0/3072.0)*PI2 - (77.0/72.0)*lambda + (11.0/24.0)*theta) *eta 
    - (15211.0/442368.0)*eta2 + (25565.0/331776.0)*eta3;
  double xcoef5 = -(113868647.0/433520640.0)*LAL_PI - (31821.0/143360.0)*LAL_PI*eta + (294941.0/3870720.0)*LAL_PI*eta2;
  double log256 = 8.0 * LAL_LN2;
  double phicoef1 =  (3715.0/8064.0)  +  (55.0/96.0) *eta;                            /* (13) */
  double phicoef2 =  (9275495.0/14450688.0) + (284875.0/258048.0)*eta + (1855.0/2048.0)*eta2;
  double phicoef3 = -(38645.0/172032.0)*LAL_PI + (65.0/2048.0)*LAL_PI*eta;
  double phicoef4 = (831032450749357.0/57682522275840.0 - (53.0/40.0)*PI2 - (107.0/56.0)*EulerGamma)
    + (-123292747421.0/4161798144.0 + (2255.0/2048.0)*PI2 + (385.0/48.0)*lambda - (55.0/16.0)*theta) * eta 
    + (154565.0/1835008.0)*eta2 - (1179625/1769472)*eta3;
  double phicoef5 =  (188516689.0/173408256.0)*LAL_PI  +  (488825.0/516096.0)*LAL_PI*eta - (141769.0/516096.0)*LAL_PI*eta2;
  double x_isco = 1.0/6.0; /* pow( (pi * f_isco)/omegacoef , 2.0/3.0); */
  int terminate=0;
  UINT4 i;
  double epochGPS = XLALGPSGetREAL8(&(IFOdata->timeData->epoch));

  /* fill `timeModelhPlus' & `timeModelhCross' with time-domain template: */
  for (i=0; i<IFOdata->timeData->data->length; ++i){
    /* determine time left until coalescence, "(t_c-t)" in (4.17)/(11): */
    t = (tc - epochGPS) - ((double)i)*IFOdata->timeData->deltaT; 
    if ((t>0.0) && (!terminate)) {  /*  (before t_c and before frequency reaches its maximum) */
      /*  determine `dimensionless time variable' tau: */
      log_tau = taucoef + log(t);                                                /*  (4.17), (11) */
      tau18   = exp(0.125 * log_tau);   /* = tau ^ (1/8) */
      tau28   = exp(0.25  * log_tau);   /* = tau ^ (2/8) */
      tau38   = exp(0.375 * log_tau);   /* = tau ^ (3/8) */
      tau48   = exp(0.5   * log_tau);   /* = tau ^ (4/8) */
      tau58   = exp(0.625 * log_tau);   /* = tau ^ (5/8) */
      tau68   = exp(0.75  * log_tau);   /* = tau ^ (6/8) */
      tau78   = exp(0.875 * log_tau);   /* = tau ^ (7/8) */
      /* determine (dimensionless) `frequency' x: */
      x = (0.25/tau28) * (1.0 + xcoef1/tau28 - (LAL_PI/5.0)/tau38
                          + xcoef2/tau48 + xcoef3/tau58
                          + (xcoef4-(107.0/3360.0)*(log_tau-log256))/tau68 
                          + xcoef5/tau78);                                        /*  (12)  */
      if ((x > x_isco) || (x < oldx)){  /* (frequency decreases  ==>  signal is terminated) */
        h_plus = h_cross = 0.0; 
        terminate = 1;
      }
      else {                    /*  (frequency still increasing  ==>  keep on computing...) */
        oldx    = x;
        sqrtx   = sqrt(x);
        /* derive angular frequency omega: (omega/pi gives frequency in Hz) */
        omega   = omegacoef*x*sqrtx;   /*  = ((c^3)/(G*mt)) * x^(3/2)                (4.13) */
        /* determine phase phi: */
	phi     = phase - (1.0/eta) * 
	  (tau58 + phicoef1*tau38 - (0.75*LAL_PI)*tau28
	   + phicoef2*tau18 + phicoef3*(log_tau-log_tau0)
	   + (phicoef4 + (107.0/448.0)*(log_tau-log256))/tau18
	   + phicoef5/tau28);                                             /*  (13)    */
        /* derive `basic phase' psi: */
        /* psi     = phi - 2.0*x*sqrtx * (log(omega)-log_omega0); */              /*  (6.12)  */
	psi     = phi - 2.0*x*sqrtx * (log(omega)-log_omega0) * (1.0-(eta/2.0)*x); /* Arun et al. (5.6) */
	sin1psi = sin(psi);      cos1psi = cos(psi);
	sin2psi = sin(2.0*psi);  cos2psi = cos(2.0*psi);
	sin3psi = sin(3.0*psi);  cos3psi = cos(3.0*psi);
	sin4psi = sin(4.0*psi);  cos4psi = cos(4.0*psi);
	sin5psi = sin(5.0*psi);  cos5psi = cos(5.0*psi);
	sin6psi = sin(6.0*psi);  cos6psi = cos(6.0*psi);
	sin7psi = sin(7.0*psi);  cos7psi = cos(7.0*psi);
        /* determine PN plus- & cross-terms: */
	Hp00    = -(1.0+ci2)*cos2psi - (si2/96.0)*(17.0+ci2);                     /*  (6.02), Arun et al (5.7a) */
	Hp05    = -(si/8.0)*dmm * ((5.0+ci2)*cos1psi - 9.0*(1.0+ci2)*cos3psi);    /*  (6.03)  */
	Hp10    = plus10a*cos2psi - plus10b*cos4psi;                              /*  (6.04)  */
	Hp15    = (si/192.0)*dmm * (plus15a*cos1psi - plus15b*cos3psi + plus15c*cos5psi) 
	  - LAL_TWOPI*(1.0+ci2)*cos2psi;                          /*  (6.05)  */
	Hp20    = plus20a*cos2psi + plus20b*cos4psi - plus20c*cos6psi
	  +si/40.0*dmm*(plus20d*sin1psi-(5.0*LAL_PI)*(5.0+ci2)*cos1psi 
			-plus20e*(1.0+ci2)*sin3psi+(135.0*LAL_PI)*(1.0+ci2)*cos3psi);   /*  (6.06)  */
        Hp25    = cos1psi*plus25a + cos2psi*plus25b + cos3psi*plus25c
	  + cos4psi*plus25d + cos5psi*plus25e + cos7psi*plus25f
	  + sin2psi*plus25g + sin4psi*plus25h;                            /*  Arun & al. (5.09) */
	Hc00    = -2.0*ci*sin2psi;                                                /*  (6.07)  */
	Hc05    = -0.75*si*ci*dmm*(sin1psi-3.0*sin3psi);                          /*  (6.08)  */
	Hc10    = cross10a*sin2psi - cross10b*sin4psi;                            /*  (6.09)  */
	Hc15    = ((si*ci)/96.0)*dmm * 
	  (cross15a*sin1psi - cross15b*sin3psi + cross15c*sin5psi)
	  -(4.0*LAL_PI)*ci*sin2psi;                                       /*  (6.10)  */
	Hc20    = cross20a*sin2psi + cross20b*sin4psi - cross20c*sin6psi
	  -0.15*si*ci*dmm*(cross20d*cos1psi+(5.0*LAL_PI)*sin1psi
			   -cross20e*cos3psi - (45.0*LAL_PI)*sin3psi);                     /*  (6.11)  */
        Hc25    = cross25a + cos2psi*cross25b + cos4psi*cross25c
	  + sin1psi*cross25d + sin2psi*cross25e + sin3psi*cross25f
	  + sin4psi*cross25g + sin5psi*cross25h + sin7psi*cross25i;       /*  Arun & al. (5.10) */
        /* and finally - the actual signal: */
	h_plus  = h_cross = constfactor * x;
	h_plus  *= Hp00 + sqrtx*(Hp05 + sqrtx*(Hp10 + sqrtx*(Hp15 + sqrtx*(Hp20 + sqrtx*Hp25))));
	h_cross *= Hc00 + sqrtx*(Hc05 + sqrtx*(Hc10 + sqrtx*(Hc15 + sqrtx*(Hc20 + sqrtx*Hc25))));/* (6.01) */
      }
    }
    else h_plus = h_cross = 0.0;  /*  (after t_c or after termination) */
    IFOdata->timeModelhPlus->data->data[i]  = h_plus;
    IFOdata->timeModelhCross->data->data[i] = h_cross;
  }
  IFOdata->modelDomain = LAL_SIM_DOMAIN_TIME;
  return;
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


void LALInferenceTemplateLALGenerateInspiral(LALInferenceIFOData *IFOdata)
/********************************************************************************************/
/* LALGenerateInspiral wrapper.																*/
/*  Required (`IFOdata->modelParams') parameters are:										*/
/*   - "m1"				(mass of object 1; REAL8, solar mass)								*/
/*   - "m2"				(mass of object 1; REAL8, solar mass)								*/
/*   - "inclination"	(inclination angle; REAL8, radians)                                 */
/*   - "coa_phase"      (phase angle; REAL8, radians)                                       */
/*   - "spin1x"			(x component of the spin of object 1; REAL8) (if SpinTaylor approx)	*/
/*   - "spin1y"			(y component of the spin of object 1; REAL8) (if SpinTaylor approx)	*/
/*   - "spin1z"			(z component of the spin of object 1; REAL8) (if SpinTaylor approx)	*/
/*   - "spin2x"			(x component of the spin of object 2; REAL8) (if SpinTaylor approx)	*/
/*   - "spin2y"			(y component of the spin of object 2; REAL8) (if SpinTaylor approx)	*/
/*   - "spin2z"			(z component of the spin of object 2; REAL8) (if SpinTaylor approx)	*/
/*	 - "shift0"			(shift offset; REAL8, radians)			                            */
/*   - "time"			(coalescence time, or equivalent/analog/similar; REAL8, GPS sec.)	*/
/*	 - "PNorder"		(Phase PN order; REAL8)												*/
/********************************************************************************************/
{
	
  static LALStatus    status;
  CoherentGW          waveform;
  SimInspiralTable    injParams;
  PPNParamStruc       ppnParams;
  Approximant			approximant=(Approximant)0;
  LALPNOrder			order=(LALPNOrder)0;
  CHAR				approximant_order[LIGOMETA_WAVEFORM_MAX];
  unsigned long				i;
  int					forceTimeLocation;
  static int sizeWarning = 0;
  
  REAL8 a1,a2,phi,shift;
  REAL8 m1,m2,mc,eta;
  REAL8 chirplength;
	
  REAL8 padding=0.4; // hard coded value found in LALInferenceReadData(). Padding (in seconds) for the tuckey window.
  UINT8 windowshift=(UINT8) ceil(padding/IFOdata->timeData->deltaT);
  
  memset( &status, 0, sizeof(LALStatus) );
  memset( &waveform, 0, sizeof(CoherentGW) );
  memset( &injParams, 0, sizeof(SimInspiralTable) );
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
	
  newswitch = 0;
  if (LALInferenceCheckVariable(IFOdata->modelParams, "newswitch"))
    newswitch = *(INT4*) LALInferenceGetVariable(IFOdata->modelParams, "newswitch"); //temporay global variable to use the new LALSTPN
  
  IFOdata->modelDomain = LAL_SIM_DOMAIN_TIME;
	
  if (LALInferenceCheckVariable(IFOdata->modelParams, "LAL_APPROXIMANT"))
    approximant = *(Approximant*) LALInferenceGetVariable(IFOdata->modelParams, "LAL_APPROXIMANT");
  else {
    XLALPrintError(" ERROR in templateLALGenerateInspiral(): (INT4) \"LAL_APPROXIMANT\" parameter not provided!\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }
	
  if (LALInferenceCheckVariable(IFOdata->modelParams, "LAL_PNORDER"))
    order = *(LALPNOrder*) LALInferenceGetVariable(IFOdata->modelParams, "LAL_PNORDER");
  else {
    XLALPrintError(" ERROR in templateLALGenerateInspiral(): (INT4) \"LAL_PNORDER\" parameter not provided!\n");
    XLAL_ERROR_VOID(XLAL_EDATA);
  }
	
  XLALInspiralGetApproximantString( approximant_order, LIGOMETA_WAVEFORM_MAX, (Approximant) approximant, (LALPNOrder)  order);
  //LALSnprintf(injParams.waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),approximant_order);
  snprintf(injParams.waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"%s",approximant_order);

  if (approximant == SpinQuadTaylor){
    injParams.qmParameter1 = 1.;
    injParams.qmParameter2 = 1.;
  }
  //	if (approximant == SpinTaylorT3){
  //    snprintf(injParams.waveform,LIGOMETA_WAVEFORM_MAX*sizeof(CHAR),"SpinTaylorFramelessthreePointFivePN");
  //  }

  /* Prefer m1 and m2 if available */
  REAL8 *m1p=NULL; 
  REAL8 *m2p=NULL;
  if( (m1p=(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"mass1")) && (m2p=(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"mass2")))
    {
      m1=*m1p; m2=*m2p;
      eta=(m1*m2)/((m1+m2)*(m1+m2));
      mc=(m1+m2)*pow(eta,0.6);
    }
  else{
    mc  = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "chirpmass");
    if (LALInferenceCheckVariable(IFOdata->modelParams,"asym_massratio")) {
      REAL8 q = *(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams,"asym_massratio");
      q2eta(q, &eta);
    }
    else
      eta = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "massratio");
    mc2masses(mc, eta, &m1, &m2);
  }
	
  injParams.mass1			= m1;				/* stellar mass */
  injParams.mass2			= m2;			    /* stellar mass */
  injParams.eta			= eta;
  injParams.mchirp		= mc;

  injParams.inclination	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "inclination");	    /* inclination in radian */
  injParams.coa_phase		= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phase");
	
	
	
  REAL8 a_spin1		= 0.0;
  if(LALInferenceCheckVariable(IFOdata->modelParams, "a_spin1"))		a_spin1		= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "a_spin1");
  REAL8 theta_spin1	= injParams.inclination; //default to spin aligned case if no angles are provided for the spins. 
  if(LALInferenceCheckVariable(IFOdata->modelParams, "theta_spin1"))	theta_spin1	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "theta_spin1");
  REAL8 phi_spin1		= 0.0;
  if(LALInferenceCheckVariable(IFOdata->modelParams, "phi_spin1"))	phi_spin1	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi_spin1");
	
  REAL8 a_spin2		= 0.0;
  if(LALInferenceCheckVariable(IFOdata->modelParams, "a_spin2"))		a_spin2		= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "a_spin2");
  REAL8 theta_spin2	= injParams.inclination; //default to spin aligned case if no angles are provided for the spins.
  if(LALInferenceCheckVariable(IFOdata->modelParams, "theta_spin2"))	theta_spin2	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "theta_spin2");
  REAL8 phi_spin2		= 0.0;
  if(LALInferenceCheckVariable(IFOdata->modelParams, "phi_spin2"))	phi_spin2	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi_spin2");
	
  injParams.spin1x = (a_spin1 * sin(theta_spin1) * cos(phi_spin1));
  injParams.spin1y = (a_spin1 * sin(theta_spin1) * sin(phi_spin1));
  injParams.spin1z = (a_spin1 * cos(theta_spin1));
	
  injParams.spin2x = (a_spin2 * sin(theta_spin2) * cos(phi_spin2));
  injParams.spin2y = (a_spin2 * sin(theta_spin2) * sin(phi_spin2));
  injParams.spin2z = (a_spin2 * cos(theta_spin2));
	
  injParams.distance	= 1.;																	/* distance set at 1 Mpc */
	
  if (IFOdata->timeData==NULL) {
    XLALPrintError(" ERROR in templateLALGenerateInspiral(): encountered unallocated 'timeData'.\n");
    XLAL_ERROR_VOID(XLAL_EFAULT);
  }
	
  ppnParams.deltaT = IFOdata->timeData->deltaT;
  double deltaT = IFOdata->timeData->deltaT;

  //injParams.f_final = IFOdata->fHigh; //(IFOdata->freqData->data->length-1) * IFOdata->freqData->deltaF;  /* (Nyquist freq.) */
  /* Check if fLow is a model parameter, otherwise use data structure definition */
  if(LALInferenceCheckVariable(IFOdata->modelParams, "fLow"))
    injParams.f_lower = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "fLow");
  else
    injParams.f_lower = IFOdata->fLow;	
  //ppnParams.fStartIn = IFOdata->fLow;
  //ppnParams.lengthIn = 0;
  //ppnParams.ppn      = NULL;


  REAL8 desired_tc		= *(REAL8 *)LALInferenceGetVariable(IFOdata->modelParams, "time");   			/* time at coalescence */

  if(desired_tc < (IFOdata->timeData->epoch.gpsSeconds + 1e-9*IFOdata->timeData->epoch.gpsNanoSeconds)){
    fprintf(stderr, "ERROR: Desired tc %f is before start of segment %f (in %s, line %d)\n",desired_tc,(IFOdata->timeData->epoch.gpsSeconds + 1e-9*IFOdata->timeData->epoch.gpsNanoSeconds), __FILE__, __LINE__);
    exit(1);
  }
	
  REAL8 instant;
  INT4 errnum;
  
  //lal_errhandler = LAL_ERR_RTRN;
  //REPORTSTATUS(&status);
  /* LAL_CALL( LALGenerateInspiral( &status, &waveform, &injParams, &ppnParams ),&status); */
  XLAL_TRY( LALGenerateInspiral( &status, &waveform, &injParams, &ppnParams ), errnum);
  //REPORTSTATUS(&status);
	
  if ( status.statusCode )
    {
      fprintf( stderr, " ERROR in templateLALGenerateInspiral(): error generating waveform. errnum=%d\n",errnum );
      REPORTSTATUS(&status);
      for (i=0; i<IFOdata->timeData->data->length; i++){
			
	IFOdata->timeModelhPlus->data->data[i] = 0.0;
	IFOdata->timeModelhPlus->data->data[i] = 0.0;
      }
      destroyCoherentGW( &waveform );	
      return;
    }
	
  instant= (IFOdata->timeData->epoch.gpsSeconds + 1e-9*IFOdata->timeData->epoch.gpsNanoSeconds)+ppnParams.tc;
	
  /* now either time-shift template or just store the time value: */
  /* (time-shifting should not be necessary in general,           */
  /* but may be neat to have for de-bugging etc.)                 */
  forceTimeLocation = 0;  /* default: zero! */
  if (forceTimeLocation) { /* time-shift the time-domain template: */
			
    chirplength=ppnParams.tc;	/*The waveform duration up to tc */
    REAL8 timeShift = desired_tc - (chirplength + IFOdata->timeData->epoch.gpsSeconds + 1e-9*IFOdata->timeData->epoch.gpsNanoSeconds);   /* This is the difference between the desired start time and the actual start time */
    INT4 integerLeftShift = ceil(-timeShift/deltaT);
    REAL8 fractionalRightShift = (deltaT*integerLeftShift+timeShift)/deltaT;
			
    for (i=0; i<IFOdata->timeData->data->length; i++){		
      if(deltaT*i>desired_tc || (i+integerLeftShift+1)>=(waveform.phi->data->length - 1) || ((long)i+integerLeftShift)<0){	//set waveform to zero after desired tc, or if need to go past end of input
	IFOdata->timeModelhPlus->data->data[i] = 0;
	IFOdata->timeModelhCross->data->data[i] = 0;		
      }
      /* Shifting waveform to account for timeShift: */
      else{
	a1  = (1.0-fractionalRightShift)*waveform.a->data->data[2*(i+integerLeftShift)]+fractionalRightShift*waveform.a->data->data[2*(i+integerLeftShift)+2];
	a2  = (1.0-fractionalRightShift)*waveform.a->data->data[2*(i+integerLeftShift)+1]+fractionalRightShift*waveform.a->data->data[2*(i+integerLeftShift)+3];
	phi     = (1.0-fractionalRightShift)*waveform.phi->data->data[i+integerLeftShift]+fractionalRightShift*waveform.phi->data->data[i+integerLeftShift+1];
	shift   = (1.0-fractionalRightShift)*waveform.shift->data->data[i+integerLeftShift]+fractionalRightShift*waveform.shift->data->data[i+integerLeftShift+1];
					
	IFOdata->timeModelhPlus->data->data[i] = a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi);
	IFOdata->timeModelhCross->data->data[i] = a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi);
      }
    }
		
  }
  else {
    /* write template (time axis) location in "->modelParams" so that     */
    /* template corresponds to stored parameter values                    */
    /* and other functions may time-shift template to where they want it: */
      
    instant=instant+(INT8)windowshift*IFOdata->timeData->deltaT; //leave enough room for the tuckey windowing of the data.
    LALInferenceSetVariable(IFOdata->modelParams, "time", &instant);
			
			
    if(waveform.a && waveform.phi){
      if(waveform.phi->data->length+2*windowshift<=IFOdata->timeData->data->length){ //check whether the IFOdata->timeData->data vector is long enough to store the waveform produced
	for (i=0; i<IFOdata->timeData->data->length; i++){
	  if(i>=(waveform.phi->data->length + windowshift) || i<windowshift){
	    IFOdata->timeModelhPlus->data->data[i] = 0;
	    IFOdata->timeModelhCross->data->data[i] = 0;		
	  }else{
	    a1		= waveform.a->data->data[2*((INT8)i-(INT8)windowshift)];
	    a2		= waveform.a->data->data[2*((INT8)i-(INT8)windowshift)+1];
	    phi     = waveform.phi->data->data[(INT8)i-(INT8)windowshift];
	    if (waveform.shift) shift   = waveform.shift->data->data[(INT8)i-(INT8)windowshift];
	    else shift = 0.0;
					
	    IFOdata->timeModelhPlus->data->data[i] = a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi);
	    IFOdata->timeModelhCross->data->data[i]= a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi);
	  }
	}
      }else{
	if (!sizeWarning) {
	  sizeWarning = 1;
	  fprintf(stderr, "WARNING: waveform.phi->data->length = %d is longer than IFOdata->timeData->data->length = %d minus windowshift = %d.\n", waveform.phi->data->length, IFOdata->timeData->data->length,(int) windowshift);
	  if(waveform.phi->data->length + (int) windowshift > IFOdata->timeData->data->length)
	    fprintf(stderr, "The waveform template used will be missing its first %d points. Consider increasing the segment length (--seglen). (in %s, line %d)\n",waveform.phi->data->length - IFOdata->timeData->data->length + (int) windowshift , __FILE__, __LINE__);
	  else
	    fprintf(stderr, "The waveform template used will have its first %d points tapered. Consider increasing the segment length (--seglen). (in %s, line %d)\n",waveform.phi->data->length - IFOdata->timeData->data->length + 2*(int)windowshift , __FILE__, __LINE__);
	}
	for (i=0; i<IFOdata->timeData->data->length; i++){
	  if((INT8)i>=(INT8)IFOdata->timeData->data->length-(INT8)windowshift || (INT8)i+(INT8)waveform.phi->data->length-(INT8)IFOdata->timeData->data->length+(INT8)windowshift < 0){
	    IFOdata->timeModelhPlus->data->data[i] = 0.0;
	    IFOdata->timeModelhCross->data->data[i] = 0.0;
	  }else{
	    a1		= waveform.a->data->data[2*((INT8)i+(INT8)waveform.phi->data->length-(INT8)IFOdata->timeData->data->length+(INT8)windowshift)];
	    a2		= waveform.a->data->data[2*((INT8)i+(INT8)waveform.phi->data->length-(INT8)IFOdata->timeData->data->length+(INT8)windowshift)+1];
	    phi     = waveform.phi->data->data[(INT8)i+(INT8)waveform.phi->data->length-(INT8)IFOdata->timeData->data->length+(INT8)windowshift];
	    if (waveform.shift) shift   = waveform.shift->data->data[(INT8)i+(INT8)waveform.phi->data->length-(INT8)IFOdata->timeData->data->length+(INT8)windowshift];
	    else shift = 0.0;
              
	    IFOdata->timeModelhPlus->data->data[i] = a1*cos(shift)*cos(phi) - a2*sin(shift)*sin(phi);
	    IFOdata->timeModelhCross->data->data[i]= a1*sin(shift)*cos(phi) + a2*cos(shift)*sin(phi);
	  }
	}
	instant-= ((INT8)waveform.phi->data->length-(INT8)IFOdata->timeData->data->length+2*(INT8)windowshift)*IFOdata->timeData->deltaT;
	LALInferenceSetVariable(IFOdata->modelParams, "time", &instant);
      }
    }else if(waveform.h){
      if(waveform.h->data->length+2*windowshift<=IFOdata->timeData->data->length){ //check whether the IFOdata->timeData->data vector is long enough to store the waveform produced
	for (i=0; i<IFOdata->timeData->data->length; i++){
	  if(i>=((unsigned long int)(waveform.h->data->length) + windowshift -1 )  || i<windowshift || isnan(waveform.h->data->data[2*(i-(INT8)windowshift)]) || isnan(waveform.h->data->data[2*(i-(INT8)windowshift)+1]) ){
	    IFOdata->timeModelhPlus->data->data[i] = 0;
	    IFOdata->timeModelhCross->data->data[i] = 0;		
	  }else{
	    IFOdata->timeModelhPlus->data->data[i] = waveform.h->data->data[2*(i-(INT8)windowshift)];
	    IFOdata->timeModelhCross->data->data[i] = waveform.h->data->data[2*(i-(INT8)windowshift)+1];
	  }
	}
      }else{
	if (!sizeWarning) {
	  sizeWarning = 1;
	  fprintf(stderr, "WARNING: waveform.h->data->length = %d is longer than IFOdata->timeData->data->length = %d minus windowshift = %d.\n", waveform.h->data->length, IFOdata->timeData->data->length, (int) windowshift);
	  if(waveform.h->data->length + (int) windowshift > IFOdata->timeData->data->length)
	    fprintf(stderr, "The waveform template used will be missing its first %d points. Consider increasing the segment length (--seglen). (in %s, line %d)\n",waveform.h->data->length - IFOdata->timeData->data->length + (int) windowshift , __FILE__, __LINE__);
	  else
	    fprintf(stderr, "The waveform template used will have its first %d points tapered. Consider increasing the segment length (--seglen). (in %s, line %d)\n",waveform.h->data->length - IFOdata->timeData->data->length + 2*(int)windowshift , __FILE__, __LINE__);
	}
	for (i=0; i<IFOdata->timeData->data->length; i++){
	  if((INT8)i>=(INT8)IFOdata->timeData->data->length-(INT8)windowshift || (INT8)i+(INT8)waveform.h->data->length-(INT8)IFOdata->timeData->data->length+(INT8)windowshift < 0 || isnan(waveform.h->data->data[2*((INT8)i+(INT8)waveform.h->data->length-(INT8)IFOdata->timeData->data->length+(INT8)windowshift)]) || isnan(waveform.h->data->data[2*((INT8)i+(INT8)waveform.h->data->length-(INT8)IFOdata->timeData->data->length+(INT8)windowshift)]+1) ){
	    IFOdata->timeModelhPlus->data->data[i] = 0.0;
	    IFOdata->timeModelhCross->data->data[i] = 0.0;
	  }else{                
	    IFOdata->timeModelhPlus->data->data[i] = waveform.h->data->data[2*((INT8)i+(INT8)waveform.h->data->length-(INT8)IFOdata->timeData->data->length+(INT8)windowshift)];
	    IFOdata->timeModelhCross->data->data[i] = waveform.h->data->data[2*((INT8)i+(INT8)waveform.h->data->length-(INT8)IFOdata->timeData->data->length+(INT8)windowshift)+1];
	  }
	}
	instant-= ((INT8)waveform.h->data->length-(INT8)IFOdata->timeData->data->length+2*(INT8)windowshift)*IFOdata->timeData->deltaT;
	LALInferenceSetVariable(IFOdata->modelParams, "time", &instant);
      }
    }else{
      for (i=0; i<IFOdata->timeData->data->length; i++){
	IFOdata->timeModelhPlus->data->data[i] = 0;
	IFOdata->timeModelhCross->data->data[i] = 0;
      }
      fprintf( stderr, " ERROR in templateLALGenerateInspiral(): no generated waveform.\n");
    }
  }
	
  if(LALInferenceCheckVariable(IFOdata->modelParams, "INFERENCE_TAPER")){
		
    if(*(LALInferenceApplyTaper*)LALInferenceGetVariable(IFOdata->modelParams, "INFERENCE_TAPER")<5 && *(LALInferenceApplyTaper*)LALInferenceGetVariable(IFOdata->modelParams, "INFERENCE_TAPER")>0){
			
      LALSimInspiralApplyTaper bookends = *(LALSimInspiralApplyTaper*) LALInferenceGetVariable(IFOdata->modelParams, "INFERENCE_TAPER");
			
      REAL4Vector *tempVec = NULL;
      tempVec = (REAL4Vector *)XLALCreateREAL4Vector(IFOdata->timeData->data->length);
			
      for (i=0; i<IFOdata->timeData->data->length; i++){
	tempVec->data[i]=(REAL4) IFOdata->timeModelhPlus->data->data[i];
      }
      XLALSimInspiralREAL4WaveTaper(tempVec,bookends);
      for (i=0; i<IFOdata->timeData->data->length; i++){
	IFOdata->timeModelhPlus->data->data[i]=(REAL8) tempVec->data[i];
      }
			
      for (i=0; i<IFOdata->timeData->data->length; i++){
	tempVec->data[i]=(REAL4) IFOdata->timeModelhCross->data->data[i];
      }
      XLALSimInspiralREAL4WaveTaper(tempVec,bookends);
      for (i=0; i<IFOdata->timeData->data->length; i++){
	IFOdata->timeModelhCross->data->data[i]=(REAL8) tempVec->data[i];
      }
      XLALDestroyREAL4Vector(tempVec);
			
    }
  }
	

	
  destroyCoherentGW( &waveform );	
	
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
/*   - if LALINFERENCE_FRAME == LALINFERENCE_FRAME_RADIATION (default)                                                   */
/*      - "inclination"	    inclination angle L.N in radians;                            REAL8                           */
/*      - "theta_spin1"     polar angle of spin 1, default to the spin aligned case;     REAL8 OPTIONAL (inclination)    */
/*      - "phi_spin1"       azimuthal angle of spin 1, default to the spin aligned case; REAL8  OPTIONAL (0.0)           */
/*      - "theta_spin2"     polar angle of spin 2, default to the spin aligned case;     REAL8 OPTIONAL (inclination)    */
/*      - "phi_spin2"       azimuthal angle of spin 1, default to the spin aligned case; REAL8  OPTIONAL (0.0)           */
/*   - else if LALINFERENCE_FRAME == LALINFERENCE_FRAME_SYSTEM                                                           */
/*      - "theta_JN");      zenith angle between J and N in radians;            REAL8                                    */
/*      - "phi_JL");        azimuthal angle of L_N on its cone about J radians; REAL8                                    */
/*      - "tilt_spin1");    zenith angle between S1 and LNhat in radians;       REAL8                                    */
/*      - "tilt_spin2");    zenith angle between S2 and LNhat in radians;       REAL8                                    */
/*      - "phi12");         difference in azimuthal angle between S1, S2 in radians;   REAL8                             */
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
  LALInferenceFrame frame=LALINFERENCE_FRAME_RADIATION;

  unsigned long	i;
  static int sizeWarning = 0;
  int ret=0;
  INT4 errnum=0;
  REAL8 instant;
  
  
  REAL8TimeSeries *hplus=NULL;  /**< +-polarization waveform [returned] */
  REAL8TimeSeries *hcross=NULL; /**< x-polarization waveform [returned] */
  COMPLEX16FrequencySeries *hptilde=NULL, *hctilde=NULL;
  
  REAL8 mc;
  REAL8 phi0, deltaT, m1, m2, spin1x, spin1y, spin1z, spin2x, spin2y, spin2z, f_low, f_start, distance, inclination;
  
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
  else if (approximant == SpinTaylorT2 || approximant == SpinTaylorT4)
    amporder = MAX_PRECESSING_AMP_PN_ORDER;
  else
    amporder = MAX_NONPRECESSING_AMP_PN_ORDER;

  if (LALInferenceCheckVariable(IFOdata->modelParams, "LALINFERENCE_FRAME"))
    frame = *(LALInferenceFrame*) LALInferenceGetVariable(IFOdata->modelParams, "LALINFERENCE_FRAME");

  REAL8 fRef = 0.0;
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

  f_start = fLow2fStart(f_low, amporder);
  f_max = 0.0; /* for freq domain waveforms this will stop at ISCO. Previously found using IFOdata->fHigh causes NaNs in waveform (see redmine issue #750)*/
  
  if(frame==LALINFERENCE_FRAME_RADIATION){
    inclination	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "inclination");	    /* inclination in radian */
    REAL8 theta_spin1	= inclination; //default to spin aligned case if no angles are provided for the spins. 
    if(LALInferenceCheckVariable(IFOdata->modelParams, "theta_spin1"))	theta_spin1	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "theta_spin1");
    REAL8 phi_spin1		= 0.0;
    if(LALInferenceCheckVariable(IFOdata->modelParams, "phi_spin1"))	phi_spin1	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi_spin1");
    
    REAL8 theta_spin2	= inclination; //default to spin aligned case if no angles are provided for the spins.
    if(LALInferenceCheckVariable(IFOdata->modelParams, "theta_spin2"))	theta_spin2	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "theta_spin2");
    REAL8 phi_spin2		= 0.0;
    if(LALInferenceCheckVariable(IFOdata->modelParams, "phi_spin2"))	phi_spin2	= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi_spin2");
    
    spin1x = (a_spin1 * sin(theta_spin1) * cos(phi_spin1));
    spin1y = (a_spin1 * sin(theta_spin1) * sin(phi_spin1));
    spin1z = (a_spin1 * cos(theta_spin1));
    
    spin2x = (a_spin2 * sin(theta_spin2) * cos(phi_spin2));
    spin2y = (a_spin2 * sin(theta_spin2) * sin(phi_spin2));
    spin2z = (a_spin2 * cos(theta_spin2));

  } else if(frame==LALINFERENCE_FRAME_SYSTEM) {
    REAL8 thetaJN = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "theta_JN");     /* zenith angle between J and N in radians */
    REAL8 phiJL = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi_JL");     /* azimuthal angle of L_N on its cone about J radians */
    REAL8 tilt1 = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "tilt_spin1");     /* zenith angle between S1 and LNhat in radians */
    REAL8 tilt2 = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "tilt_spin2");     /* zenith angle between S2 and LNhat in radians */
    REAL8 phi12 = *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "phi12");      /* difference in azimuthal angle btwn S1, S2 in radians */
 
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

  if(LALInferenceCheckVariable(IFOdata->modelParams, "spin1"))		spin1z		= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "spin1");
  if(LALInferenceCheckVariable(IFOdata->modelParams, "spin2"))		spin2z		= *(REAL8*) LALInferenceGetVariable(IFOdata->modelParams, "spin2");
  
  distance	= LAL_PC_SI * 1.0e6;        /* distance (1 Mpc) in units of metres */
	
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
  LALSimInspiralWaveformFlags *waveFlags = XLALSimInspiralCreateWaveformFlags();
  if(LALInferenceCheckVariable(IFOdata->modelParams, "spinO")) XLALSimInspiralSetSpinOrder(waveFlags, *(LALSimInspiralSpinOrder*) LALInferenceGetVariable(IFOdata->modelParams, "spinO"));
  if(LALInferenceCheckVariable(IFOdata->modelParams, "tideO")) XLALSimInspiralSetTidalOrder(waveFlags, *(LALSimInspiralTidalOrder*) LALInferenceGetVariable(IFOdata->modelParams, "tideO"));
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
    
	XLAL_TRY(ret=XLALSimInspiralChooseFDWaveform(&hptilde, &hctilde, phi0,
            deltaF, m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z,
            spin2x, spin2y, spin2z, f_start, f_max, distance, inclination,
            lambda1, lambda2, waveFlags, nonGRparams, amporder, order,
            approximant), errnum);

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
        IFOdata->freqModelhPlus->data->data[i].real_FIXME = 0.0;
        IFOdata->freqModelhPlus->data->data[i].imag_FIXME = 0.0;
      }
    }
    for (i=0; i<IFOdata->freqModelhCross->data->length; ++i) {
      dataPtr = hctilde->data->data;
      if(i < hctilde->data->length){
        IFOdata->freqModelhCross->data->data[i] = dataPtr[i];
      }else{
        IFOdata->freqModelhCross->data->data[i].real_FIXME = 0.0;
        IFOdata->freqModelhCross->data->data[i].imag_FIXME = 0.0;
      }
    }
    
    
    /* Destroy the WF flags and the nonGr params */
    XLALSimInspiralDestroyWaveformFlags(waveFlags);
    XLALSimInspiralDestroyTestGRParam(nonGRparams);
    
    instant= (IFOdata->timeData->epoch.gpsSeconds + 1e-9*IFOdata->timeData->epoch.gpsNanoSeconds);
    LALInferenceSetVariable(IFOdata->modelParams, "time", &instant);
    
  } else {

    XLAL_TRY(ret=XLALSimInspiralChooseTDWaveform(&hplus, &hcross, phi0, deltaT,
            m1*LAL_MSUN_SI, m2*LAL_MSUN_SI, spin1x, spin1y, spin1z,
            spin2x, spin2y, spin2z, f_start, fRef, distance,
            inclination, lambda1, lambda2, waveFlags, nonGRparams,
            amporder, order, approximant), errnum);
    XLALSimInspiralDestroyWaveformFlags(waveFlags);
    XLALSimInspiralDestroyTestGRParam(nonGRparams);
    if (ret == XLAL_FAILURE || hplus == NULL || hcross == NULL)
      {
	XLALPrintError(" ERROR in XLALSimInspiralChooseWaveform(): error generating waveform. errnum=%d\n",errnum );
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
    size_t maxShift = (size_t)lround(4.255e-4/hplus->deltaT); 

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



static void destroyCoherentGW( CoherentGW *waveform )
{
  if ( waveform->h )
    {
      XLALDestroyREAL4VectorSequence( waveform->h->data );
      LALFree( waveform->h );
    }
  if ( waveform->a )
    {
      XLALDestroyREAL4VectorSequence( waveform->a->data );
      LALFree( waveform->a );
    }
  if ( waveform->phi )
    {
      XLALDestroyREAL8Vector( waveform->phi->data );
      LALFree( waveform->phi );
    }
  if ( waveform->f )
    {
      XLALDestroyREAL4Vector( waveform->f->data );
      LALFree( waveform->f );
    }
  if ( waveform->shift )
    {
      XLALDestroyREAL4Vector( waveform->shift->data );
      LALFree( waveform->shift );
    }
	
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


