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


void templateLALwrap2(LALIFOData *IFOdata)
{
//  static LALStatus status;
  static InspiralTemplate params;
/*  static REAL4Vector *LALSignal=NULL;
  UINT4 n;
  long i;
  double mc   = *(REAL8*) getVariable(IFOdata->modelParams, "chirpmass");
  double eta  = *(REAL8*) getVariable(IFOdata->modelParams, "massratio");
  double phi  = *(REAL8*) getVariable(IFOdata->modelParams, "phase");
  double iota = *(REAL8*) getVariable(IFOdata->modelParams, "inclination");
  double tc   = *(REAL8*) getVariable(IFOdata->modelParams, "time");
  double chirptime;*/

  params.OmegaS      = 0.0;     /* (?) */
  params.Theta       = 0.0;     /* (?) */
  /* params.Zeta2    = 0.0; */  /* (?) */
  params.ieta        = 1; 
  params.nStartPad   = 0;
  params.nEndPad     = 0;
  params.massChoice  = m1Andm2;
	params.approximant = TaylorT1;//approximant;  /*  TaylorT1, ...              */
	params.order       = 0;//order;        /*  0=Newtonian, ..., 7=3.5PN  */
	params.fLower      = 40.0;//DF->minF * 0.9;
	params.fCutoff     = 500.0;//	(DF->FTSize-1)*DF->FTDeltaF;  /* (Nyquist freq.) */
	params.tSampling   = 1.0/100.0;//DF->dataDeltaT;
  params.startTime   = 0.0;

  
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
/* parameter, one may request a 2.0PN (insead of 2.5PN)      */
/* template.                                                 */
/* Signal's amplitude corresponds to a luminosity distance   */
/* of 1 Mpc; re-scaling will need to be taken care of e.g.   */
/* in the calling likelihood function.                       */
/*************************************************************/
{
  double mc   = *(REAL8*) getVariable(IFOdata->modelParams, "chirpmass");
  double eta  = *(REAL8*) getVariable(IFOdata->modelParams, "massratio");
  double phi  = *(REAL8*) getVariable(IFOdata->modelParams, "phase");
  double iota = *(REAL8*) getVariable(IFOdata->modelParams, "inclination");
  double tc   = *(REAL8*) getVariable(IFOdata->modelParams, "time");
  double PNOrder = 2.5;  /* (default) */
  double fraction = (0.5+sqrt(0.25-eta)) / (0.5-sqrt(0.25-eta));
  double mt = mc * ((pow(1.0+fraction,0.2) / pow(fraction,0.6)) /* (total mass) */
                    + (pow(1.0+1.0/fraction,0.2) / pow(1.0/fraction,0.6)));
  double log_q   = log(mt) + log(LAL_PI) + log(LAL_G_SI) - 3.0*log((double) LAL_C_SI);
  double log_eta = log(eta);
  double a[5];
  double ampliConst, plusCoef, crossCoef;
  double dataStart, NDeltaT, phaseArg;
  double f, f01, f02, f04, f06, f07, f10, Psi, twopitc;
  double plusRe, plusIm, crossRe, crossIm;
  int i, lower, upper;
 
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
  for (i=0; i<IFOdata->timeData->data->length; ++i){
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

/*    if(!IFOdata->freqModelhPlus)
      IFOdata->freqModelhPlus=(COMPLEX16FrequencySeries *)XLALCreateCOMPLEX16FrequencySeries("freqData",&(IFOdata->timeData->epoch),0.0,
                IFOdata->freqData->deltaF,&lalDimensionlessUnit,IFOdata->freqData->data->length);
*/
//These aren't properly allocated yet.
//    IFOdata->freqModelhPlus->data->data[i].re  = plusRe;
//    IFOdata->freqModelhPlus->data->data[i].im  = plusIm;
//    IFOdata->freqModelhCross->data->data[i].re = crossRe;
//    IFOdata->freqModelhCross->data->data[i].im = crossIm;
  }
  IFOdata->modelDomain = frequencyDomain;
  return;
}
