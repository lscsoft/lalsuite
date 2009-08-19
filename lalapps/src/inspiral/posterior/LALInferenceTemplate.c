/* 
 *  LALInferenceTemplate.c:  Bayesian Followup, template calls to LAL template functions. Temporary GeneratePPN
 *
 *  Copyright (C) 2009 Ilya Mandel, Vivien Raymond, Christian Roever, Marc van der Sluys and John Veitch
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
#include <lal/Date.h>
#include <lal/VectorOps.h>
#include <lal/TimeFreqFFT.h>
#include "LALInference.h"


void LALTemplateGeneratePPN(LALIFOData *IFOdata){

	static LALStatus stat;								/* status structure */	

	IFOdata->modelDomain = timeDomain;
	
	
	REAL4 m1=*(REAL4 *)getVariable(IFOdata->modelParams,"m1");			/* binary masses */	
	REAL4 m2=*(REAL4 *)getVariable(IFOdata->modelParams,"m2");
	
	REAL4 dist=1.0;//*(REAL4 *)getVariable(IFOdata->modelParams,"dist");      /* binary distance SET AS FIDUCIAL  - 1 Mpc*/
	REAL4 inc=*(REAL4 *)getVariable(IFOdata->modelParams,"inc");		/* inclination and coalescence phase */
	REAL4 phii=*(REAL4 *)getVariable(IFOdata->modelParams,"phii");
	
	REAL8 desired_tc = *(REAL8 *)getVariable(IFOdata->modelParams,"time");  
	
	REAL4 f_min = IFOdata->fLow; 
	REAL4 f_max = IFOdata->fHigh;			/* start and stop frequencies */
	
	
	REAL8 deltaT = IFOdata->timeData->deltaT;				/* waveform-generation data-sampling interval */
	INT4 order = 4;										/* PN order */
	
	/* Other variables. */
	UINT4 i;                      /* index */
	PPNParamStruc params;         /* input parameters */
	CoherentGW waveform;          /* output waveform */

	/* Make sure that values won't crash the system or anything. */
//	CHECKVAL( order, -1, 5 );
//	CHECKVAL( dt, LAL_REAL4_MIN, LAL_REAL4_MAX );
//	CHECKVAL( deltat, 0.0, LAL_REAL4_MAX );
	
	
	/*******************************************************************
	 * INPUT SETUP                                                     *
	 *******************************************************************/
	
	/* Fixed parameters. Set them when injecting....*/
	params.position.latitude = 0.0;//*(REAL4 *)getVariable(ifo->modelParams,"latitude");
	params.position.longitude = 0.0;//*(REAL4 *)getVariable(ifo->modelParams,"longitude");
	params.position.system = COORDINATESYSTEM_EQUATORIAL;
	params.psi = 0.0;
	params.lengthIn = 0;
	params.epoch.gpsSeconds = IFOdata->timeData->epoch.gpsSeconds;
	params.epoch.gpsNanoSeconds = IFOdata->timeData->epoch.gpsNanoSeconds;

	/* Variable parameters. */

	params.deltaT = deltaT;
	params.mTot = m1 + m2;
	params.eta = m1*m2/( params.mTot*params.mTot );
	params.inc = inc;
	params.phi = phii;//0.0;
	params.d = dist*LAL_PC_SI*1.0e3;
	params.fStartIn = f_min;
	params.fStopIn = f_max;

	/* PPN parameter. */
	params.ppn = NULL;
	LALSCreateVector( &stat, &(params.ppn), order + 1 );
	params.ppn->data[0] = 1.0;
	if ( order > 0 )
		params.ppn->data[1] = 0.0;
	for ( i = 2; i <= (UINT4)( order ); i++ )
		params.ppn->data[i] = 1.0;
	/* Output parameters. */
	memset( &waveform, 0, sizeof(CoherentGW) );
	
	
	/*******************************************************************
	 * OUTPUT GENERATION                                               *
	 *******************************************************************/

	/* Generate waveform. */
	LALGeneratePPNInspiral( &stat, &waveform, &params );
	
	REAL8 chirplength=params.tc;	/*The waveform duration up to tc */
	printf("desired_tc %10.10f chirplength %g epoch %10.10f\n", desired_tc, chirplength, IFOdata->timeData->epoch.gpsSeconds + 1e-9*IFOdata->timeData->epoch.gpsNanoSeconds);
	
	
	/* This is the difference between the desired start time and the actual start time */
	REAL8 timeShift = desired_tc - (chirplength + IFOdata->timeData->epoch.gpsSeconds + 1e-9*IFOdata->timeData->epoch.gpsNanoSeconds);
	
	printf("Timeshift %g\n", timeShift);
	
	if(desired_tc < (IFOdata->timeData->epoch.gpsSeconds + 1e-9*IFOdata->timeData->epoch.gpsNanoSeconds)){
		fprintf(stderr, "ERROR: Desired tc is before start of segment\n");
		exit(1);
	}
	if(timeShift > 0){ //If we rightshift, we should window first
		//if(!IFOData->window)
		//	IFOdata[i].window=XLALCreateTukeyREAL8Window(seglen,(REAL8)2.0*padding*SampleRate/(REAL8)seglen);
		//XLALDDVectorMultiply(waveform.a->data->data,waveform.a->data->data,IFOdata[i].window->data);
		//fprintf(stderr, "ERROR: Desired tc is greater than generated tc; can't right-shift waveform\n");
		//exit(1);
	}
	
	/* Check if sampling interval was too large. */
	if ( params.dfdt > 2.0 ) {
		printf(
				 "Waveform sampling interval is too large:\n"
				 "\tmaximum df*dt = %f", params.dfdt );
		//WARNING( message );
	}
	
	
	/* Shifting waveform to account for timeShift: */
			
	REAL8 p, ap, ac;
	UINT4 integerLeftShift = ceil(-timeShift/deltaT);
	REAL8 fractionalRightShift = (deltaT*integerLeftShift+timeShift)/deltaT;
		
	//printf("deltaT %g, iLS %d, fRS %g\n", deltaT, integerLeftShift, fractionalRightShift);
	//printf("t %d, a %d, phi %d\n", IFOdata->timeData->data->length, waveform.a->data->length, waveform.phi->data->length);
	
	UINT4 length = IFOdata->timeData->data->length;//waveform.a->data->length-1;  //Why are waveform.a and waveform.phi of same length, yet .a contains both plus and cross?
	REAL8 *phiData = waveform.phi->data->data;
	REAL4 *aData = waveform.a->data->data;
	
	for(i=0; i<length; i++){
		if(deltaT*i>desired_tc || i+integerLeftShift+1>=waveform.phi->data->length - 1
			|| i+integerLeftShift<0){	//set waveform to zero after desired tc, or if need to go past end of input
			IFOdata->timeModelhPlus->data->data[i] = 0;
			IFOdata->timeModelhCross->data->data[i] = 0;
		}
		else{
			p = (1.0-fractionalRightShift)*phiData[i+integerLeftShift] + fractionalRightShift*phiData[i+integerLeftShift+1];
			ap = (1.0-fractionalRightShift)*aData[2*(i+integerLeftShift)] + fractionalRightShift*aData[2*(i+integerLeftShift)+2];
			ac = (1.0-fractionalRightShift)*aData[2*(i+integerLeftShift)+1] + fractionalRightShift*aData[2*(i+integerLeftShift)+3];
			IFOdata->timeModelhPlus->data->data[i] = ap*cos(p);
			IFOdata->timeModelhCross->data->data[i] = ap*sin(p);
		}
	}
/*			
			REAL8 dx = deltat/dt;
			REAL8 xMax = waveform.a->data->length - 1;
			REAL8 *phiData = waveform.phi->data->data;
			//REAL4 *fData = waveform.f->data->data;
			REAL4 *aData = waveform.a->data->data;
			for ( ; x < xMax; x += dx, t += deltat ) {
				UINT4 j = floor( x );
				if(j < IFOdata->timeData->data->length ){
					REAL8 frac = x - j;
					REAL8 p = frac*phiData[j+1] + ( 1.0 - frac )*phiData[j];
					//REAL8 f = frac*fData[j+1] + ( 1.0 - frac )*fData[j];
					REAL8 ap = frac*aData[2*j+2] + ( 1.0 - frac )*aData[2*j];
					REAL8 ac = frac*aData[2*j+3] + ( 1.0 - frac )*aData[2*j+1];
					IFOdata->timeModelhPlus->data->data[j] = ap*cos( p );
					IFOdata->timeModelhCross->data->data[j] = ac*sin( p );
				}
			}
*/
	//INT4 k = 0;
	//for(k=0 ; k < IFOdata->timeData->data->length; k++ ){
//		fprintf(stdout,"%d\t%13.6e\t%13.6e\n",k,IFOdata->timeModelhPlus->data->data[k],IFOdata->timeModelhCross->data->data[k]);
//	    }
	
	
	/*******************************************************************
	 * CLEANUP                                                         *
	 *******************************************************************/
	
	LALSDestroyVector( &stat, &(params.ppn) );
	LALSDestroyVectorSequence( &stat, &(waveform.a->data) );
	LALSDestroyVector( &stat, &(waveform.f->data) );
	LALDDestroyVector( &stat, &(waveform.phi->data) );
	LALFree( waveform.a );
	LALFree( waveform.f );
	LALFree( waveform.phi );
	
	
	
	LALCheckMemoryLeaks();
	
	
//	INFO( GENERATEPPNINSPIRALTESTC_MSGENORM );
//	return GENERATEPPNINSPIRALTESTC_ENORM;
	
}




void templateStatPhase(LALIFOData *IFOdata)
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
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **********/
/* Required (`IFOdata->modelParams') parameters are:                  */
/*   - "chirpmass"    (REAL8,units of solar masses)                   */
/*   - "massratio"    (symmetric mass ratio:  0 < eta <= 0.25, REAL8) */
/*   - "phase"        (coalescence phase; REAL8, radians)             */
/*   - "time"         (coalescence time; REAL8, GPS seconds)          */
/*   - "inclination"  (inclination angle, REAL8, radians)             */
/**********************************************************************/
{
  double mc   = *(REAL8*) getVariable(IFOdata->modelParams, "chirpmass");
  double eta  = *(REAL8*) getVariable(IFOdata->modelParams, "massratio");
  double phi  = *(REAL8*) getVariable(IFOdata->modelParams, "phase");
  double iota = *(REAL8*) getVariable(IFOdata->modelParams, "inclination");
  double tc   = *(REAL8*) getVariable(IFOdata->modelParams, "time");
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
  int i, lower, upper;
 
  if (IFOdata->timeData==NULL)
    die(" ERROR in templateStatPhase(): encountered unallocated 'timeData'.\n");
  if ((IFOdata->freqModelhPlus==NULL) || (IFOdata->freqModelhCross==NULL)) 
    die(" ERROR in templateStatPhase(): encountered unallocated 'freqModelhPlus/-Cross'.\n");
  if (checkVariable(IFOdata->modelParams, "PNOrder"))
    PNOrder = *(REAL8*) getVariable(IFOdata->modelParams, "PNOrder");
  if ((PNOrder!=2.5) && (PNOrder!=2.0)) die(" ERROR in templateStatPhase(): only PN orders 2.0 or 2.5 allowed.");
  ampliConst  = 0.5*log(5.0) + (5.0/6.0)*log(LAL_G_SI) - log(2.0) - 0.5*log(6.0) - (2.0/3.0)*log(LAL_PI) - 1.5*log((double)LAL_C_SI);
  ampliConst  = exp(ampliConst + 0.5*log_eta + (5.0/6.0)*log(mt) - (log(LAL_PC_SI)+log(1.0e+6)));
  /* leaving out the following term makes freqDomain template scaling match that of "XLALREAL8TimeFreqFFT()" output: */
  /* ampliConst /= IFOdata->timeData->deltaT; */
  plusCoef  = (-0.5*(1.0+pow(cos(iota),2.0)));
  crossCoef = (-1.0*cos(iota));
  dataStart = XLALGPSGetREAL8(&(IFOdata->timeData->epoch));
  twopitc = LAL_TWOPI * (tc - dataStart);
  a[0] =  exp(log(3.0/128.0) - (5.0/3.0)*log_q - log_eta);
  a[1] =  exp(log(3715.0/84.0+55.0*eta) - log(384.0) - log_eta - log_q);
  a[2] = -exp(log(48.0*LAL_PI/128.0) - (2.0/3.0)*log_q - log_eta);
  a[3] =  exp(log(3.0/128.0) - log_eta - (1.0/3.0)*log_q
              + log(15293365.0/508032.0 + ((27145.0/504.0) + (3085.0/72.0)*eta)*eta));
  a[4] =  exp(log(LAL_PI/128.0)-log_eta+log(38645.0/252.0+5.0*eta));
  NDeltaT = ((double) IFOdata->timeData->data->length) * IFOdata->timeData->deltaT;
  lower = ceil(IFOdata->fLow * NDeltaT);
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
    IFOdata->freqModelhPlus->data->data[i].re  = plusRe;
    IFOdata->freqModelhPlus->data->data[i].im  = plusIm;
    IFOdata->freqModelhCross->data->data[i].re = crossRe;
    IFOdata->freqModelhCross->data->data[i].im = crossIm;
  }
  IFOdata->modelDomain = frequencyDomain;
  return;
}



void templateNullFreqdomain(LALIFOData *IFOdata)
/**********************************************/
/* returns a frequency-domain 'null' template */
/* (all zeroes, implying no signal present).  */
/**********************************************/
{
  int i;
  if ((IFOdata->freqModelhPlus==NULL) || (IFOdata->freqModelhCross==NULL)) 
    die(" ERROR in templateNullFreqdomain(): encountered unallocated 'freqModelhPlus/-Cross'.\n");
  for (i=0; i<IFOdata->freqModelhPlus->data->length; ++i){
    IFOdata->freqModelhPlus->data->data[i].re  = 0.0;
    IFOdata->freqModelhPlus->data->data[i].im  = 0.0;
    IFOdata->freqModelhCross->data->data[i].re = 0.0;
    IFOdata->freqModelhCross->data->data[i].im = 0.0;
  }
  IFOdata->modelDomain = frequencyDomain;
  return;
}



void templateNullTimedomain(LALIFOData *IFOdata)
/*********************************************/
/* returns a time-domain 'null' template     */
/* (all zeroes, implying no signal present). */
/*********************************************/
{
  int i;
  if ((IFOdata->timeModelhPlus==NULL) || (IFOdata->timeModelhCross==NULL)) 
    die(" ERROR in templateNullTimedomain(): encountered unallocated 'timeModelhPlus/-Cross'.\n");
  for (i=0; i<IFOdata->timeModelhPlus->data->length; ++i){
    IFOdata->timeModelhPlus->data->data[i]  = 0.0;
    IFOdata->timeModelhCross->data->data[i] = 0.0;
  }
  IFOdata->modelDomain = timeDomain;
  return;
}



/* ============ LAL template wrapper function: ========== */



void mc2masses(double mc, double eta, double *m1, double *m2)
/*  Compute individual companion masses (m1, m2)   */
/*  for given chirp mass (m_c) & mass ratio (eta)  */
/*  (note: m2 >= m1).                              */
{
  double root = sqrt(0.25-eta);
  double fraction = (0.5+root) / (0.5-root);
  *m1 = mc * (pow(1+fraction,0.2) / pow(fraction,0.6));
  *m2 = mc * (pow(1+1.0/fraction,0.2) / pow(1.0/fraction,0.6));
}



void templateLAL(LALIFOData *IFOdata)
/*************************************************************************************************/
/* Wrapper function to call LAL functions for waveform generation.                               */
/* Will always return frequency-domain templates (numerically FT'ed                              */
/* in case the LAL function returns time-domain).                                                */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`IFOdata->modelParams') parameters are:                                             */
/*   - "chirpmass"        (REAL8,units of solar masses)                                          */
/*   - "massratio"        (symmetric mass ratio:  0 < eta <= 0.25, REAL8)                        */
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
  static InspiralTemplate params;
  static REAL4Vector *LALSignal=NULL;
  UINT4 n;
  long i; 
  long j, jmax;
  double pj, pmax, pleft, pright;

  double mc       = *(REAL8*) getVariable(IFOdata->modelParams, "chirpmass");
  double eta      = *(REAL8*) getVariable(IFOdata->modelParams, "massratio");
  double phi      = *(REAL8*) getVariable(IFOdata->modelParams, "phase");       /* here: startPhase !! */
  double tc       = *(REAL8*) getVariable(IFOdata->modelParams, "time");
  double iota     = *(REAL8*) getVariable(IFOdata->modelParams, "inclination");
  int approximant, order;
  int FDomain;    /* (denotes domain of the _LAL_ template!) */
  double m1, m2, chirptime, deltaT;
  double plusCoef  = -0.5 * (1.0 + pow(cos(iota),2.0));
  double crossCoef = -1.0 * cos(iota);
  double instant;
  int forceTimeLocation;
  double twopit, f, deltaF, re, im, templateReal, templateImag;

  if (checkVariable(IFOdata->modelParams, "LAL_APPROXIMANT"))
    approximant = *(INT4*) getVariable(IFOdata->modelParams, "LAL_APPROXIMANT");
  else die(" ERROR in templateLAL(): (INT4) \"LAL_APPROXIMANT\" parameter not provided!\n");

  if (checkVariable(IFOdata->modelParams, "LAL_PNORDER"))
    order = *(INT4*) getVariable(IFOdata->modelParams, "LAL_PNORDER");
  else die(" ERROR in templateLAL(): (INT4) \"LAL_PNORDER\" parameter not provided!\n");

  /*fprintf(stdout, " templateLAL() - approximant = %d,  PN order = %d\n", approximant, order);*/

  /* little consistency check (otherwise no output without warning): */
  if (((approximant==EOBNR) || (approximant==EOB)) && (order!=LAL_PNORDER_PSEUDO_FOUR))
    die(" ERROR in templateLAL(): \"EOB\" and \"EOBNR\" templates require \"LAL_PNORDER_PSEUDO_FOUR\" PN order!\n");  

  if (IFOdata->timeData==NULL) 
    die(" ERROR in templateLAL(): encountered unallocated 'timeData'.\n");
  if ((IFOdata->freqModelhPlus==NULL) || (IFOdata->freqModelhCross==NULL)) 
    die(" ERROR in templateLAL(): encountered unallocated 'freqModelhPlus/-Cross'.\n");
  deltaT = IFOdata->timeData->deltaT;

  mc2masses(mc, eta, &m1, &m2);
  params.OmegaS      = 0.0;     /* (?) */
  params.Theta       = 0.0;     /* (?) */
  /* params.Zeta2    = 0.0; */  /* (?) */
  params.ieta        = 1; 
  params.nStartPad   = 0;
  params.nEndPad     = 0;
  params.massChoice  = m1Andm2;
  params.approximant = approximant;  /*  TaylorT1, ...   */
  params.order       = order;        /*  Newtonian, ...  */
  params.fLower      = IFOdata->fLow * 0.9;
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
      || (params.approximant == IMRPhenomA))
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

  LALInspiralParameterCalc(&status, &params);
  chirptime = params.tC;
  if ((params.approximant != TaylorF2) && (params.approximant != BCV)) {
    params.startTime = (tc - XLALGPSGetREAL8(&IFOdata->timeData->epoch)) - chirptime;
    LALInspiralParameterCalc(&status, &params); /* (re-calculation necessary? probably not...) */
  }

  /* compute "params.signalAmplitude" slot: */
  LALInspiralRestrictedAmplitude(&status, &params);

  /* figure out inspiral length & set `n': */
  /* LALInspiralWaveLength(&status, &n, params); */
  n = IFOdata->timeData->data->length;

  /* domain of LAL template as returned by LAL function: */
  FDomain = ((params.approximant == TaylorF1)
             || (params.approximant == TaylorF2)
             || (params.approximant == PadeF1)
             || (params.approximant == BCV));
  if (FDomain && (n % 2 != 0)){
    fprintf(stderr, " ERROR in templateLAL(): frequency-domain LAL waveforms require even number of samples!\n");
    fprintf(stderr, "                         (N = IFOdata->timeData->data->length = %d)\n", n);
    exit(1);
  }

  /* allocate (temporary) waveform vector: */
  LALCreateVector(&status, &LALSignal, n);
  for (i=0; i<n; ++i) LALSignal->data[i] = 0.0;


  /*--  ACTUAL WAVEFORM COMPUTATION:  --*/
  if (status.statusCode != 0) {
    fprintf(stderr, " ERROR in templateLAL(): encountered non-zero status code.\n");
    fprintf(stderr, " Template parameters:\n");
    printVariables(IFOdata->modelParams);
    fprintf(stderr, " LAL Status:\n");
    REPORTSTATUS(&status);
    exit(1);
  }
  /* REPORTSTATUS(&status); */
  LALInspiralWave(&status, LALSignal, &params);
  /* REPORTSTATUS(&status); */
  if (status.statusCode != 0) {
    fprintf(stderr, "\n ERROR in templateLAL(): \"LALInspiralWave()\" call returned with non-zero status.\n");
    fprintf(stderr, " Template parameters:\n");
    printVariables(IFOdata->modelParams);
    fprintf(stderr, " LAL Status:\n");
    REPORTSTATUS(&status);
    exit(1);
  }

  if (! FDomain) {   /*  (LAL function returns TIME-DOMAIN template)       */
    /* copy over, normalise: */
    for (i=0; i<n; ++i) {
      IFOdata->timeModelhPlus->data->data[i]  = LALSignal->data[i];
      IFOdata->timeModelhCross->data->data[i] = 0.0;  /* (no cross waveform) */
    }
    LALDestroyVector(&status, &LALSignal);
    /* apply window & execute FT of plus component: */
    if (IFOdata->window==NULL) 
      die(" ERROR in templateLAL(): ran into uninitialized 'IFOdata->window'.\n");
    XLALDDVectorMultiply(IFOdata->timeModelhPlus->data, IFOdata->timeModelhPlus->data, IFOdata->window->data);
    if (IFOdata->timeToFreqFFTPlan==NULL)
      die(" ERROR in templateLAL(): ran into uninitialized 'IFOdata->timeToFreqFFTPlan'.\n");
    XLALREAL8TimeFreqFFT(IFOdata->freqModelhPlus, IFOdata->timeModelhPlus, IFOdata->timeToFreqFFTPlan);
  }

  else {             /*  (LAL function returns FREQUENCY-DOMAIN template)  */
    /* copy over: */
    IFOdata->freqModelhPlus->data->data[0].re = ((REAL8) LALSignal->data[0]);
    IFOdata->freqModelhPlus->data->data[0].im = 0.0;
    for (i=1; i<IFOdata->freqModelhPlus->data->length-1; ++i) {
      IFOdata->freqModelhPlus->data->data[i].re = ((REAL8) LALSignal->data[i]);
      IFOdata->freqModelhPlus->data->data[i].im = ((REAL8) LALSignal->data[n-i]);
    }
    IFOdata->freqModelhPlus->data->data[IFOdata->freqModelhPlus->data->length-1].re = LALSignal->data[IFOdata->freqModelhPlus->data->length-1];
    IFOdata->freqModelhPlus->data->data[IFOdata->freqModelhPlus->data->length-1].im = 0.0;
    LALDestroyVector(&status, &LALSignal);
    /* nomalise (apply same scaling as in XLALREAL8TimeFreqFFT()") : */
    for (i=0; i<IFOdata->freqModelhPlus->data->length; ++i) {
      IFOdata->freqModelhPlus->data->data[i].re *= ((REAL8) n) * deltaT;
      IFOdata->freqModelhPlus->data->data[i].im *= ((REAL8) n) * deltaT;
    }
  }

  /* (now frequency-domain plus-waveform has been computed, either directly or via FFT)   */

  /*  cross waveform is "i x plus" :  */
  for (i=1; i<IFOdata->freqModelhCross->data->length-1; ++i) {
    IFOdata->freqModelhCross->data->data[i].re = -IFOdata->freqModelhPlus->data->data[i].im;
    IFOdata->freqModelhCross->data->data[i].im = IFOdata->freqModelhPlus->data->data[i].re;
    // consider inclination angle's effect:
    IFOdata->freqModelhPlus->data->data[i].re  *= plusCoef;
    IFOdata->freqModelhPlus->data->data[i].im  *= plusCoef;
    IFOdata->freqModelhCross->data->data[i].re *= crossCoef;
    IFOdata->freqModelhCross->data->data[i].im *= crossCoef;
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
      || (params.approximant == TaylorF2))
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
    executeInvFT(IFOdata);
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
        templateReal = IFOdata->freqModelhPlus->data->data[i].re;
        templateImag = IFOdata->freqModelhPlus->data->data[i].im;
        IFOdata->freqModelhPlus->data->data[i].re = templateReal*re - templateImag*im;
        IFOdata->freqModelhPlus->data->data[i].im = templateReal*im + templateImag*re;
        templateReal = IFOdata->freqModelhCross->data->data[i].re;
        templateImag = IFOdata->freqModelhCross->data->data[i].im;
        IFOdata->freqModelhCross->data->data[i].re = templateReal*re - templateImag*im;
        IFOdata->freqModelhCross->data->data[i].im = templateReal*im + templateImag*re;
      }
    }
    else {
      /* write template (time axis) location in "->modelParams" so that     */
      /* template corresponds to stored parameter values                    */
      /* and other functions may time-shift template to where they want it: */
      setVariable(IFOdata->modelParams, "time", &instant);
    }
  }

  IFOdata->modelDomain = frequencyDomain;
  return;
}



void template3525TD(LALIFOData *IFOdata)
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
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *****************************/
/* Required (`IFOdata->modelParams') parameters are:                                         */
/*   - "chirpmass"        (REAL8,units of solar masses)                                      */
/*   - "massratio"        (symmetric mass ratio:  0 < eta <= 0.25, REAL8)                    */
/*   - "phase"            (here: 'startPhase', not coalescence phase; REAL8, radians)        */
/*   - "time"             (coalescence time, or equivalent/analog/similar; REAL8, GPS sec.)  */
/*   - "inclination"      (inclination angle, REAL8, radians)                                */
/*********************************************************************************************/
{
  double mc    = *(REAL8*) getVariable(IFOdata->modelParams, "chirpmass");  /* chirp mass m_c, solar masses  */
  double eta   = *(REAL8*) getVariable(IFOdata->modelParams, "massratio");  /* mass ratio eta, dimensionless */
  double tc    = *(REAL8*) getVariable(IFOdata->modelParams, "time");       /* coalescence time, GPS sec.    */
  double phase = *(REAL8*) getVariable(IFOdata->modelParams, "phase");      /* coalescence phase, rad        */
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
  double ci  =  cos(*(REAL8*) getVariable(IFOdata->modelParams, "inclination"));
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
                                                                                      /* (6.01); distance is 1 Mpc here. */
  double x, sqrtx, oldx=0.0;
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
  int i, terminate=0;
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
  IFOdata->modelDomain = timeDomain;
}
