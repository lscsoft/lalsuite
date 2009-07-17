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


void LALTemplateWrapper(LALIFOData *IFOdata){

	static LALStatus stat;								/* status structure */	

	IFOdata->modelDomain = timeDomain;
	
	
	REAL4 m1=*(REAL4 *)getVariable(IFOdata->modelParams,"m1");			/* binary masses */	
	REAL4 m2=*(REAL4 *)getVariable(IFOdata->modelParams,"m2");
	
	REAL4 dist=1.0;//*(REAL4 *)getVariable(IFOdata->modelParams,"dist");      /* binary distance SET AS FIDUCIAL */
	REAL4 inc=*(REAL4 *)getVariable(IFOdata->modelParams,"inc");		/* inclination and coalescence phase */
	REAL4 phii=*(REAL4 *)getVariable(IFOdata->modelParams,"phii");

	
	REAL4 f_min = IFOdata->fLow; 
	REAL4 f_max = IFOdata->fHigh;			/* start and stop frequencies */
	
	
	REAL8 dt = IFOdata->timeData->deltaT;					/* sampling interval */
	REAL8 deltat = IFOdata->timeData->deltaT;				/* wave sampling interval */
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
	params.epoch.gpsSeconds = (INT4)(0);
	params.epoch.gpsNanoSeconds = (INT4)(0);

	/* Variable parameters. */

	params.deltaT = dt;
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
	
	/* Print termination information. */
	//snprintf( message, MSGLENGTH, "%d: %s", params.termCode,
	//		 params.termDescription );
	//INFO( message );
	
	/* Print coalescence phase.
	 snprintf( message, MSGLENGTH,
	 "Waveform ends %.3f cycles before coalescence",
	 -waveform.phi->data->data[waveform.phi->data->length-1]
	 / (REAL4)( LAL_TWOPI ) ); */
	/*{
		INT4 code = sprintf( message,
							"Waveform ends %.3f cycles before coalescence",
							-waveform.phi->data->data[waveform.phi->data->length
													  -1]
							/ (REAL4)( LAL_TWOPI ) );
		if ( code >= MSGLENGTH || code < 0 ) {
			ERROR( GENERATEPPNINSPIRALTESTC_EPRINT,
				  GENERATEPPNINSPIRALTESTC_MSGEPRINT, 0 );
			return GENERATEPPNINSPIRALTESTC_EPRINT;
		}
	}
	INFO( message );*/
	
	/* Check if sampling interval was too large. */
	if ( params.dfdt > 2.0 ) {
		printf(
				 "Waveform sampling interval is too large:\n"
				 "\tmaximum df*dt = %f", params.dfdt );
		//WARNING( message );
	}
	
	/* Renormalize phase. */
	//phii -= waveform.phi->data->data[0];
	//for ( i = 0; i < waveform.phi->data->length; i++ )
	//	waveform.phi->data->data[i] += phii;
	
	/* Write output. */
//		if ( ( fp = fopen( outfile, "w" ) ) == NULL ) {
//			ERROR( GENERATEPPNINSPIRALTESTC_EFILE,
//				  GENERATEPPNINSPIRALTESTC_MSGEFILE, outfile );
//			return GENERATEPPNINSPIRALTESTC_EFILE;
//		}
		
		/* Amplitude and phase functions: */
//		if ( deltat == 0.0 ) {
//			REAL8 t = 0.0; /* time */
//			for ( i = 0; i < waveform.a->data->length; i++, t += dt )
//				fprintf( fp, "%f %.3f %10.3e %10.3e %10.3e\n", t,
//						waveform.phi->data->data[i],
//						waveform.f->data->data[i],
//						waveform.a->data->data[2*i],
//						waveform.a->data->data[2*i+1] );
//		}
		
		/* Waveform: */
			REAL8 t = 0.0;
			REAL8 x = 0.0;
			REAL8 dx = deltat/dt;
			REAL8 xMax = waveform.a->data->length - 1;
			REAL8 *phiData = waveform.phi->data->data;
			//REAL4 *fData = waveform.f->data->data;
			REAL4 *aData = waveform.a->data->data;
			for ( ; x < xMax; x += dx, t += deltat ) {
				UINT4 j = floor( x );
				REAL8 frac = x - j;
				REAL8 p = frac*phiData[j+1] + ( 1.0 - frac )*phiData[j];
				//REAL8 f = frac*fData[j+1] + ( 1.0 - frac )*fData[j];
				REAL8 ap = frac*aData[2*j+2] + ( 1.0 - frac )*aData[2*j];
				REAL8 ac = frac*aData[2*j+3] + ( 1.0 - frac )*aData[2*j+1];
				
				if(j < IFOdata->timeData->data->length ){
					//fprintf(stdout,"%13.6e\t%13.6e\n",ap*cos( p ),ac*sin( p ));
					IFOdata->timeModelhPlus->data->data[j] = ap*cos( p );
					IFOdata->timeModelhCross->data->data[j] = ac*sin( p );
				}
			}
		
		
	
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
  //static LALStatus status;
  //static InspiralTemplate params;
  //static REAL4Vector *LALSignal=NULL;
  //UINT4 n;
  
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
/*************************************************************/
{
  double mc   = *(REAL8*) getVariable(IFOdata->modelParams, "chirpmass");
  double eta  = *(REAL8*) getVariable(IFOdata->modelParams, "massratio");
  double phi  = *(REAL8*) getVariable(IFOdata->modelParams, "phase");
  double iota = *(REAL8*) getVariable(IFOdata->modelParams, "inclination");
  double tc   = *(REAL8*) getVariable(IFOdata->modelParams, "time");
  double PNOrder=2.5;
  double root = sqrt(0.25-eta);
  double fraction = (0.5+root) / (0.5-root);
  double inversefraction = (0.5-root) / (0.5+root);
  double mt = mc * ((pow(1+fraction,0.2) / pow(fraction,0.6)) /* (total mass) */
                    + (pow(1+inversefraction,0.2) / pow(inversefraction,0.6)));
  double log_q   = log(mt) + log(LAL_PI) + log(LAL_G_SI) - 3.0*log((double) LAL_C_SI);
  double log_eta = log(eta);
  double a[5];
  long i;
  double f, f01, f02, f04, f06, f07, f10, Psi, twopitc;
  double ampliConst;
  //double sineCoef, cosineCoef;
  double plusCoef, crossCoef;
  double NDeltaT, phaseArg;
  double plusRe, plusIm, crossRe, crossIm;
  double dataStart;
  int lower, upper;
 
  if (checkVariable(IFOdata->modelParams, "PNOrder"))
    PNOrder = *(REAL8*) getVariable(IFOdata->modelParams, "PNOrder");
  if ((PNOrder!=2.5) && (PNOrder!=2.0)) die(" ERROR in templateStatPhase(): only PN orders 2.0 & 2.5 allowed.");
  ampliConst  = 0.5*log(5.0)+(5.0/6.0)*log(LAL_G_SI)-log(2.0)-0.5*log(6.0)-(2.0/3.0)*log(LAL_PI)-1.5*log((double)LAL_C_SI);
  ampliConst  = exp(ampliConst+0.5*log_eta+(5.0/6.0)*log(mt)-(log(LAL_PC_SI)+log(1e6)));
  ampliConst /= IFOdata->timeData->deltaT;
  //cosineCoef = Fplus  * (-0.5*(1.0+pow(cos(iota),2.0)));
  //sineCoef   = Fcross * (-1.0*cos(iota));
  plusCoef  = ampliConst * (-0.5*(1.0+pow(cos(iota),2.0)));
  crossCoef = ampliConst * (-1.0*cos(iota));
  //twopitc = 2.0 * pi * (vectorGetValue(parameter,"time") - DF->dataStart);
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
  for (i=0; i<IFOdata->timeData->data->length; ++i){
    if ((i > upper) || (i < lower)) /* (no computations outside freq. range) */
      plusRe = plusIm = crossRe = crossIm = 0.0;
    else {
      f    = ((double)i) / NDeltaT;
      f01  = pow(f, -1.0/6.0);             /* = f^-1/6  */
      f02  = f01*f01;                      /* = f^-2/6  */
      f04  = f02*f02;                      /* = f^-4/6  */
      f06  = f04*f02;                      /* = f^-6/6  */
      f07  = f06*f01;                      /* = f^-7/6  */
      f10  = f06*f04;                      /* = f^-10/6 */
      Psi = a[0]*f10 + a[1]*f06 + a[2]*f04 + a[3]*f02;
      if (PNOrder>2.0) /*  5th coefficient ignored for 2.0 PN order  */
        Psi += a[4]*log(f); 
      phaseArg = Psi + twopitc*f + phi;
      plusRe  =  f07 * cos(phaseArg);
      plusIm  =  f07 * sin(phaseArg);
      crossRe =  -1.0*plusIm * crossCoef;
      crossIm =  plusRe * crossCoef;
      plusRe  *= plusCoef;
      plusIm  *= plusCoef;
      /* copy over to IFOdata: */
      IFOdata->freqModelhPlus->data->data[i].re  = plusRe;
      IFOdata->freqModelhPlus->data->data[i].im  = plusIm;
      IFOdata->freqModelhCross->data->data[i].re = crossRe;
      IFOdata->freqModelhCross->data->data[i].im = crossIm;
    }
  }
  IFOdata->modelDomain = frequencyDomain;
  return;
}
