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



/* ============ Christian's attempt of a LAL template wrapper function: ========== */


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
/*********************************************************************************************/
/* Wrapper function to call LAL functions for waveform generation.                           */
/* Will always return frequency-domain templates (numerically FT'ed                          */
/* in case the LAL function returns time-domain).                                            */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* Required (`IFOdata->modelParams') parameters are:                                         */
/*   - "chirpmass"        (REAL8,units of solar masses)                                      */
/*   - "massratio"        (symmetric mass ratio:  0 < eta <= 0.25, REAL8)                    */
/*   - "phase"            (here: 'startPhase', not coalescence phase; REAL8, radians)        */
/*   - "time"             (coalescence time, or equivalent/analog/similar; REAL8, GPS sec.)  */
/*   - "inclination"      (inclination angle, REAL8, radians)                                */
/*   - "LAL_APPROXIMANT"  (INT4 value corresponding to `enum approximant' definition         */
/*                         in `LALInspiral.h'. Allowed values by now are:                    */
/*                         TaylorT2, TaylorT3, EOBNR)                                        */
/*   - "LAL_PNORDER"      (INT4 value corresponding to `enum LALPNOrder' definition          */
/*                         in `LALInspiral.h'.)                                              */
/*********************************************************************************************/
{
  static LALStatus status;
  static InspiralTemplate params;
  static REAL4Vector *LALSignal=NULL;
  UINT4 n;
  long i;

  double mc       = *(REAL8*) getVariable(IFOdata->modelParams, "chirpmass");
  double eta      = *(REAL8*) getVariable(IFOdata->modelParams, "massratio");
  double phi      = *(REAL8*) getVariable(IFOdata->modelParams, "phase");       /* here: startPhase !! */
  double tc       = *(REAL8*) getVariable(IFOdata->modelParams, "time");
  double iota     = *(REAL8*) getVariable(IFOdata->modelParams, "inclination");
  int approximant = *(INT4*) getVariable(IFOdata->modelParams, "LAL_APPROXIMANT");
  int order       = *(INT4*) getVariable(IFOdata->modelParams, "LAL_PNORDER");
  int FDomain;    /* (denotes domain of the _LAL_ template!) */
  double m1, m2, chirptime;
  double plusCoef  = -0.5*(1.0+pow(cos(iota),2.0));
  double crossCoef = -1.0*cos(iota);

  if (IFOdata->timeData==NULL) 
    die(" ERROR in templateLALwrap(): encountered unallocated 'timeData'.\n");
  if ((IFOdata->freqModelhPlus==NULL) || (IFOdata->freqModelhCross==NULL)) 
    die(" ERROR in templateLALwrap(): encountered unallocated 'freqModelhPlus/-Cross'.\n");

  mc2masses(mc, eta, &m1, &m2);
  params.OmegaS      = 0.0;     /* (?) */
  params.Theta       = 0.0;     /* (?) */
  /* params.Zeta2    = 0.0; */  /* (?) */
  params.ieta        = 1; 
  params.nStartPad   = 0;
  params.nEndPad     = 0;
  params.massChoice  = m1Andm2;
  params.approximant = approximant;  /*  TaylorT1, ...              */
  params.order       = order;        /*  0=Newtonian, ..., 7=3.5PN  */
  params.fLower      = IFOdata->fLow * 0.9;
  params.fCutoff     = (IFOdata->freqData->data->length-1) * IFOdata->freqData->deltaF;  /* (Nyquist freq.) */
  params.tSampling   = 1.0 / IFOdata->timeData->deltaT;
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
    fprintf(stderr," WARNING: 'LALInspiralSetup()' (called within 'LALInspiralWavelength()')\n");
    fprintf(stderr,"          requires (tSampling > 2 x fCutoff) !!\n");
    fprintf(stderr," (settings are:  tSampling = %f s,  fCutoff = %f Hz)  \n", params.tSampling, params.fCutoff);
    exit(1);
  }

  /* ensure compatible sampling rate: */
  if ((params.approximant == EOBNR)
      && (fmod(log((double)params.tSampling)/log(2.0),1.0) != 0.0)) {
    fprintf(stderr, " ERROR: \"EOBNR\" templates require power-of-two sampling rates!\n");
    fprintf(stderr, "        (params.tSampling = %f Hz)\n", params.tSampling);
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
    fprintf(stderr, " ERROR: frequency-domain LAL waveforms require even number of samples!\n");
    fprintf(stderr, "        (n = %d)\n", n);
    exit(1);
  }

  /* allocate (temporary) waveform vector: */
  LALCreateVector(&status, &LALSignal, n);
  for (i=0; i<n; ++i) LALSignal->data[i] = 0.0;


  /*--  ACTUAL WAVEFORM COMPUTATION:  --*/
  /* REPORTSTATUS(&status); */
  LALInspiralWave(&status, LALSignal, &params);
  /* REPORTSTATUS(&status); */


  if (! FDomain) {   /*  (LAL function returns time-domain template)       */
    /* copy over, normalise: */
    for (i=0; i<n; ++i) {
      IFOdata->timeModelhPlus->data->data[i]  = LALSignal->data[i] * ((REAL8) n);
      IFOdata->timeModelhCross->data->data[i] = 0.0;  /* (no cross waveform) */
    }
    LALDestroyVector(&status, &LALSignal);
    /* apply window & execute FT of plus component: */
    if (IFOdata->window==NULL) 
      die(" ERROR in templateLALwrap(): ran into uninitialized 'IFOdata->window'.\n");
    XLALDDVectorMultiply(IFOdata->timeModelhPlus->data, IFOdata->timeModelhPlus->data, IFOdata->window->data);
    if (IFOdata->timeToFreqFFTPlan==NULL)
      die(" ERROR in templateLALwrap(): ran into uninitialized 'IFOdata->timeToFreqFFTPlan'.\n");
    XLALREAL8TimeFreqFFT(IFOdata->freqModelhPlus, IFOdata->timeModelhPlus, IFOdata->timeToFreqFFTPlan);
  }

  else {             /*  (LAL function returns frequency-domain template)  */
    /* copy over, normalise: */
    IFOdata->freqModelhPlus->data->data[0].re = ((REAL8) LALSignal->data[0]) * ((REAL8) n);
    IFOdata->freqModelhPlus->data->data[0].im = 0.0;
    for (i=1; i<IFOdata->freqModelhPlus->data->length-1; ++i) {
      IFOdata->freqModelhPlus->data->data[i].re = ((REAL8) LALSignal->data[i]) * ((REAL8) n);
      IFOdata->freqModelhPlus->data->data[i].im = ((REAL8) LALSignal->data[n-i]) * ((REAL8) n);
    }
    IFOdata->freqModelhPlus->data->data[IFOdata->freqModelhPlus->data->length-1].re = LALSignal->data[IFOdata->freqModelhPlus->data->length-1] * ((double) n);
    IFOdata->freqModelhPlus->data->data[IFOdata->freqModelhPlus->data->length-1].im = 0.0;
    LALDestroyVector(&status, &LALSignal);
  }

  /* (now frequency-domain plus-waveform has been computed, either directly or via FFT) */

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

  IFOdata->modelDomain = frequencyDomain;
  return;
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
                    + (pow(1.0+1.0/fraction,0.2) / pow(1.0/fraction,0.6)));  /* (total mass) */
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
  ampliConst  = 0.5*log(5.0)+(5.0/6.0)*log(LAL_G_SI)-log(2.0)-0.5*log(6.0)-(2.0/3.0)*log(LAL_PI)-1.5*log((double)LAL_C_SI);
  ampliConst  = exp(ampliConst+0.5*log_eta+(5.0/6.0)*log(mt)-(log(LAL_PC_SI)+log(1e6)));
  ampliConst /= IFOdata->timeData->deltaT;
  plusCoef  = ampliConst * (-0.5*(1.0+pow(cos(iota),2.0)));
  crossCoef = ampliConst * (-1.0*cos(iota));
  dataStart = XLALGPSGetREAL8(&(IFOdata->timeData->epoch));
  twopitc = LAL_TWOPI * (tc - dataStart);
  a[0] =  exp(log(3.0/128.0) - (5.0/3.0)*log_q - log_eta);
  a[1] =  exp(log(3715.0/84.0+55.0*eta) - log(1.0/384.0) - log_eta - log_q);
  a[2] = -exp(log(48.0*LAL_PI/128.0) - (2.0/3.0)*log_q - log_eta);
  a[3] =  exp(log(3.0/128.0) - log_eta - (1.0/3.0)*log_q
              + log(15293365.0/508032.0+(27145.0/504.0)*eta+(3085.0/72.0)*exp(2.0*log_eta)));
  a[4] =  exp(log(LAL_PI/128.0)-log_eta+log(38645.0/252.0+5.0*eta));
 
  NDeltaT = IFOdata->timeData->data->length * IFOdata->timeData->deltaT;
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
      plusRe  =  f07 * cos(phaseArg);
      plusIm  =  f07 * sin(phaseArg);
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
