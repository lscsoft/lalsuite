#include <lal/LALSTPNWaveform2.h>
#include <lal/LALSTPNWaveformErrors.h>

typedef struct LALSTPNstructparams {
	REAL8 eta, m1m2, m2m1;
	REAL8 wdotnew, wdotorb[9], wspin15, wspin20;
	REAL8 LNhdot15, LNhdot20;
	REAL8 S1dot15, S2dot15, Sdot20;
	REAL8 epnorb[9];
} LALSTPNparams;

static void XLALSTPNAdaptiveSetParams(LALSTPNparams *mparams,InspiralTemplate *params,InspiralInit *paramsInit) {
	mparams->m2m1 = params->mass2 / params->mass1;
  mparams->m1m2 = params->mass1 / params->mass2;

	mparams->eta = (params->mass1 * params->mass2) / pow(params->mass1 + params->mass2,2.0);

  mparams->wdotnew = (96.0/5.0) * params->eta;

  for(int j = LAL_PNORDER_NEWTONIAN; j <= 8; j++) {
    mparams->wdotorb[j] = mparams->epnorb[j] = 0.0;
  }

  mparams->epnorb[0] = 1.0;

	for(unsigned int j = LAL_PNORDER_NEWTONIAN; j <= params->order; j++) {
		mparams->wdotorb[j] = paramsInit->ak.ST[j];
	}

	/* note: in the original code epnorb[2] is set even for PNORDER_NEWTONIAN or PNORDER_HALF */
	if (params->order >= LAL_PNORDER_ONE) {
		mparams->epnorb[2]  = -(1.0/12.0) * (9.0 + params->eta);
	}

	for(int j = params->order + 1; j <= 8; j++) {
  	mparams->wdotorb[j] = 0;
	}

	if (params->order >= LAL_PNORDER_ONE_POINT_FIVE) {
    mparams->wspin15 	= -(1.0/12.0);
    mparams->LNhdot15 = 0.5;
    mparams->S1dot15 	= (4.0 + 3.0 * mparams->m2m1) / 2.0 ;
    mparams->S2dot15 	= (4.0 + 3.0 * mparams->m1m2) / 2.0 ;
	} else {
		mparams->wspin15 	= 0.0;
    mparams->LNhdot15 = 0.0;
    mparams->S1dot15 	= 0.0;
    mparams->S2dot15 	= 0.0;
	}

	if (params->order >= LAL_PNORDER_TWO) {
    mparams->wspin20 	= -(1.0/48.0) / params->eta;	/* TO DO: these will give infinity in test mass case */
    mparams->LNhdot20 = -1.5 / params->eta;
    mparams->Sdot20 	= 0.5;

		mparams->epnorb[4]  = (1.0/24.0) * (-81.0 + 57.0*params->eta - params->eta*params->eta);
	} else {
    mparams->wspin20 	= 0.0;
    mparams->LNhdot20 = 0.0;
    mparams->Sdot20 	= 0.0;
	}

  if (params->order >= LAL_PNORDER_THREE) {
    mparams->epnorb[6]  = ( - (675.0/64.0)
		      									+ ( (209323.0/4032.0) -(205.0/96.0)*LAL_PI*LAL_PI - (110.0/9.0) * (-1987.0/3080.0) ) * params->eta
		                        - (155.0/96.0)*(params->eta*params->eta) - (35.0/5184.0)*(params->eta*params->eta*params->eta) );
	}

  if (params->order == LAL_PNORDER_THREE) {
    mparams->wdotorb[(int)(LAL_PNORDER_THREE+1)] = paramsInit->ak.ST[(int)(LAL_PNORDER_THREE+1)];
	}

  if (params->order == LAL_PNORDER_THREE_POINT_FIVE) {
    mparams->wdotorb[8] = paramsInit->ak.ST[8];
	}

	return;
}

#define UNUSED(expr) do { (void)(expr); } while (0)

static int XLALSTPNAdaptiveTest(double t,const double values[],double dvalues[],void *mparams) {
	REAL8 omega, v, test;
	LALSTPNparams *params = (LALSTPNparams*)mparams;

	UNUSED(t);

	omega = values[1];
	v = pow(omega,oneby3);

  test = -0.5 * params->eta * ( (2.0/3.0) * (1.0/v) * params->epnorb[0]	+
                                params->epnorb[1] +
																(4.0/3.0) * v * (params->epnorb[2] +
																	(5.0/4.0) * v * (params->epnorb[3] +
																		(6.0/5.0) * v * (params->epnorb[4] +
																			(7.0/6.0) * v * (params->epnorb[5] +
																				(8.0/7.0) * v * (params->epnorb[6] +
																					(9.0/8.0) * v * (params->epnorb[7] +
																						(10.0/9.0)* v * params->epnorb[8]
																													)
																												)
																											)
																										)
																									)
																								)
															);

    test = -0.5 * params->eta * (
				 (2.0/3.0) * (1.0/v) * params->epnorb[0]
				 + params->epnorb[1]
				 + (4.0/3.0) * v * (params->epnorb[2]
				 + (5.0/4.0) * v * (params->epnorb[3]
				 + (6.0/5.0) * v * (params->epnorb[4]
				 + (7.0/6.0) * v * (params->epnorb[5]
				 + (8.0/7.0) * v * (params->epnorb[6]
				 + (9.0/8.0) * v * (params->epnorb[7]
				 + (10.0/9.0)* v *  params->epnorb[8] )))))) );

	if (params->wspin15 != 0.0) {
		REAL8 LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z;
		REAL8 v2, dotLNS1, dotLNS2;

		LNhx = values[2]; LNhy  = values[3]; LNhz = values[4] ;
		S1x  = values[5]; S1y   = values[6]; S1z  = values[7] ;
		S2x  = values[8]; S2y   = values[9]; S2z  = values[10];

		v2 = v * v;

	  dotLNS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
	  dotLNS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);

		test += -0.5 * params->eta * (5.0/3.0) * v2 * ( (8.0/3.0 + 2.0*params->m2m1)*dotLNS1 + (8.0/3.0 + 2.0*params->m1m2)*dotLNS2 );

	  if (params->wspin20 != 0.0) {
			REAL8 dotS1S2;

			dotS1S2 = (S1x*S2x + S1y*S2y + S1z*S2z);

	    test += -(v*v2)  * (dotS1S2 - 3.0 * dotLNS1 * dotLNS2);
		}
	}

	if (test > 0.0) {							 									   /* energy test fails! */
		return LALSTPN_TEST_ENERGY;
	} else if (dvalues[1] < 0.0) { 									   /* omegadot < 0! */
		return LALSTPN_TEST_OMEGADOT;
	} else if isnan(omega) {													 /* omega is nan */
		return LALSTPN_TEST_OMEGANAN;
	} else {
		return GSL_SUCCESS;
	}
}

static int XLALSTPNAdaptiveDerivatives(double t,const double values[],double dvalues[],void *mparams) {
	/* coordinates and derivatives */
  REAL8  omega,  LNhx,  LNhy,  LNhz,  S1x,  S1y,  S1z,  S2x,  S2y,  S2z;
  REAL8 ds, domega, dLNhx, dLNhy, dLNhz, dS1x, dS1y, dS1z, dS2x, dS2y, dS2z;

	/* auxiliary variables */
	REAL8 v, v2, v3, v4, v7, v11;
	REAL8 dotLNS1, dotLNS2, dotS1S2;
	REAL8 omega2, tmpx, tmpy, tmpz;
  REAL8 LNmag, crossx, crossy, crossz;

  LALSTPNparams *params = (LALSTPNparams*)mparams;

	UNUSED(t);

	/* copy variables */
	// UNUSED!!: REAL8 s    = values[0];
        omega = values[1];
	LNhx = values[2]; LNhy  = values[3]; LNhz = values[4] ;
	S1x  = values[5]; S1y   = values[6]; S1z  = values[7] ;
	S2x  = values[8]; S2y   = values[9]; S2z  = values[10];

	if (omega <= 0.0) {
		return LALSTPN_DERIVATIVE_OMEGANONPOS;
	}

  v = pow(omega,oneby3);
  v2  = v * v; v3 = v2 * v;	v4 = v3 * v; v7 = v4 * v3; v11 = v7 * v4;

  dotLNS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
  dotLNS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);
  dotS1S2 = (S1x*S2x  + S1y*S2y  + S1z*S2z );

	/* domega */
  domega = params->wdotorb[0] + v * (params->wdotorb[1] +
 																	 	 v * (params->wdotorb[2] +
																		    	v * (params->wdotorb[3] +
																				     	 v * (params->wdotorb[4] +
																						      	v * (params->wdotorb[5] +
																											 	 v * (params->wdotorb[6] + params->wdotorb[7] * log(omega) +
																											      	v * params->wdotorb[8]
																											     	 )
																												)
																								 	 )
																							)
																			 	 )
																  	);

	domega += params->wspin15 * omega * ( LNhx * (113.0 * S1x + 113.0 * S2x + 75.0 * params->m2m1 * S1x + 75.0 * params->m1m2 * S2x) +
																	      LNhy * (113.0 * S1y + 113.0 * S2y + 75.0 * params->m2m1 * S1y + 75.0 * params->m1m2 * S2y) +
																	      LNhz * (113.0 * S1z + 113.0 * S2z + 75.0 * params->m2m1 * S1z + 75.0 * params->m1m2 * S2z) );

  domega += params->wspin20 * v4 * (247.0 * dotS1S2 - 721.0 * dotLNS1 * dotLNS2);

  domega *= params->wdotnew * v11; /* Newtonian term */

  /* dLN */
  omega2 = omega * omega;

  tmpx = (4.0 + 3.0*params->m2m1) * S1x + (4.0 + 3.0*params->m1m2) * S2x;
  tmpy = (4.0 + 3.0*params->m2m1) * S1y + (4.0 + 3.0*params->m1m2) * S2y;
  tmpz = (4.0 + 3.0*params->m2m1) * S1z + (4.0 + 3.0*params->m1m2) * S2z;

  dLNhx = params->LNhdot15 * omega2 * (-tmpz*LNhy + tmpy*LNhz);
  dLNhy = params->LNhdot15 * omega2 * (-tmpx*LNhz + tmpz*LNhx);
  dLNhz = params->LNhdot15 * omega2 * (-tmpy*LNhx + tmpx*LNhy);

  tmpx = dotLNS2 * S1x + dotLNS1 * S2x;
  tmpy = dotLNS2 * S1y + dotLNS1 * S2y;
  tmpz = dotLNS2 * S1z + dotLNS1 * S2z;

  dLNhx += params->LNhdot20 * v7 * (-tmpz*LNhy + tmpy*LNhz);
  dLNhy += params->LNhdot20 * v7 * (-tmpx*LNhz + tmpz*LNhx);
  dLNhz += params->LNhdot20 * v7 * (-tmpy*LNhx + tmpx*LNhy);

	/* dS1 */
	LNmag = params->eta / v;

	crossx = (LNhy*S1z - LNhz*S1y);
	crossy = (LNhz*S1x - LNhx*S1z);
	crossz = (LNhx*S1y - LNhy*S1x);

	dS1x = params->S1dot15 * omega2 * LNmag * crossx;
	dS1y = params->S1dot15 * omega2 * LNmag * crossy;
	dS1z = params->S1dot15 * omega2 * LNmag * crossz;

	tmpx = S1z*S2y - S1y*S2z;
	tmpy = S1x*S2z - S1z*S2x;
	tmpz = S1y*S2x - S1x*S2y;

	dS1x += params->Sdot20 * omega2 * (tmpx - 3.0 * dotLNS2 * crossx);
	dS1y += params->Sdot20 * omega2 * (tmpy - 3.0 * dotLNS2 * crossy);
	dS1z += params->Sdot20 * omega2 * (tmpz - 3.0 * dotLNS2 * crossz);

	/* dS2 */
  crossx = (LNhy*S2z - LNhz*S2y);
  crossy = (LNhz*S2x - LNhx*S2z);
  crossz = (LNhx*S2y - LNhy*S2x);

  dS2x = params->S2dot15 * omega2 * LNmag * crossx;
  dS2y = params->S2dot15 * omega2 * LNmag * crossy;
  dS2z = params->S2dot15 * omega2 * LNmag * crossz;

  dS2x += params->Sdot20 * omega2 * (-tmpx - 3.0 * dotLNS1 * crossx);
  dS2y += params->Sdot20 * omega2 * (-tmpy - 3.0 * dotLNS1 * crossy);
  dS2z += params->Sdot20 * omega2 * (-tmpz - 3.0 * dotLNS1 * crossz);

	/* dphi */
  if(LNhx*LNhx + LNhy*LNhy > 0.0) {
    ds = omega - (-LNhz * (LNhy*dLNhx - LNhx*dLNhy) / (LNhx*LNhx + LNhy*LNhy)); /* term in parens is alphadotcosi */
  } else {
    ds = omega;
  }

	dvalues[0] = ds   ; dvalues[1] = domega;
	dvalues[2] = dLNhx; dvalues[3] = dLNhy ; dvalues[4] = dLNhz;
	dvalues[5] = dS1x ; dvalues[6] = dS1y  ; dvalues[7] = dS1z ;
	dvalues[8] = dS2x ; dvalues[9] = dS2y  ; dvalues[10]= dS2z ;

	return GSL_SUCCESS;
}


void
LALSTPNAdaptiveWaveformEngine( LALStatus *status,
                							 REAL4Vector *signalvec1,REAL4Vector *signalvec2,
                							 REAL4Vector *a,REAL4Vector *ff,REAL8Vector *phi,REAL4Vector *shift,
                							 UINT4 *countback,
                							 InspiralTemplate *params,InspiralInit *paramsInit
                						 )
{
	/* PN parameters */
  LALSTPNparams mparams;

	/* needed for integration */
  ark4GSLIntegrator *integrator;
	unsigned int len;
	int intreturn;
	REAL8 yinit[11];
  REAL8Array *yout;

  /* other computed values */
  REAL8 unitHz, dt, m, lengths, norm;

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

 	/* Make sure parameter and waveform structures exist. */
  ASSERT(params,     status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT(paramsInit, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

  /* units */
  unitHz = params->totalMass * LAL_MTSUN_SI * (REAL8)LAL_PI;
  dt = 1.0/params->tSampling;	  /* tSampling is in Hz, so dt is in seconds */
  m = params->totalMass * LAL_MTSUN_SI;

  /* length estimation (Newtonian); since integration is adaptive, we could use a better estimate */
  lengths = (5.0/256.0) * pow(LAL_PI,-8.0/3.0) * pow(params->chirpMass * LAL_MTSUN_SI * params->fLower,-5.0/3.0) / params->fLower;

  /* setup coefficients for PN equations */
	XLALSTPNAdaptiveSetParams(&mparams,params,paramsInit);

  /* initialize the coordinates */
	yinit[0] = 0.0;                     							/* vphi */
	yinit[1] = params->fLower * unitHz; 							/* omega (really pi M f) */

	yinit[2] = sin(params->inclination);							/* LNh(x,y,z) */
	yinit[3] = 0.0;
	yinit[4] = cos(params->inclination);

	norm = pow(params->mass1/params->totalMass,2.0);
	yinit[5] = norm * params->spin1[0];								/* S1(x,y,z) */
	yinit[6] = norm * params->spin1[1];
	yinit[7] = norm * params->spin1[2];

	norm = pow(params->mass2/params->totalMass,2.0);	/* S2(x,y,z) */
	yinit[8] = norm * params->spin2[0];
	yinit[9] = norm * params->spin2[1];
	yinit[10]= norm * params->spin2[2];

  xlalErrno = 0;

	/* allocate the integrator */
	integrator = XLALAdaptiveRungeKutta4Init(11,XLALSTPNAdaptiveDerivatives,XLALSTPNAdaptiveTest,1.0e-6,1.0e-6);
  if (!integrator) {
		fprintf(stderr,"LALSTPNWaveform2: Cannot allocate integrator.\n");
    if (XLALClearErrno() == XLAL_ENOMEM)
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
    else
      ABORTXLAL(status);
  }

	/* stop the integration only when the test is true */
	integrator->stopontestonly = 1;

	/* run the integration; note: time is measured in units of total mass */
	len = XLALAdaptiveRungeKutta4(integrator,(void *)&mparams,yinit,0.0,lengths/m,dt/m,&yout);

	intreturn = integrator->returncode;
	XLALAdaptiveRungeKutta4Free(integrator);

	if (!len) {
    if (XLALClearErrno() == XLAL_ENOMEM) {
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
    } else {
			fprintf(stderr,"LALSTPNWaveform2: integration failed with errorcode %d.\n",intreturn);
			ABORTXLAL(status);
		}
	}

	/* report on abnormal termination (TO DO: throw some kind of LAL error?) */
	if (intreturn != 0 && intreturn != LALSTPN_TEST_ENERGY && intreturn != LALSTPN_TEST_OMEGADOT) {
		fprintf(stderr,"LALSTPNWaveform2 WARNING: integration terminated with code %d.\n",intreturn);
    fprintf(stderr,"                          Waveform parameters were m1 = %e, m2 = %e, s1 = (%e,%e,%e), s2 = (%e,%e,%e), inc = %e.",
     							 params->mass1, params->mass2,
     							 params->spin1[0], params->spin1[1], params->spin1[2],
     							 params->spin2[0], params->spin2[1], params->spin2[2],
     							 params->inclination);
	}

	/* check that we're not above Nyquist */
	if (yinit[1]/unitHz > 0.5 * params->tSampling) {
		fprintf(stderr,"LALSTPNWaveform2 WARNING: final frequency above Nyquist.\n");
	}

	/* if we have enough space, compute the waveform components; otherwise abort */
  if ((signalvec1 && len >= signalvec1->length) || (ff && len >= ff->length)) {
		if (signalvec1) {
			fprintf(stderr,"LALSTPNWaveform2: no space to write in signalvec1: %d vs. %d\n",len,signalvec1->length);
		} else if (ff) {
			fprintf(stderr,"LALSTPNWaveform2: no space to write in ff: %d vs. %d\n",len,ff->length);
		} else {
			fprintf(stderr,"LALSTPNWaveform2: no space to write anywhere!\n");
		}
		ABORT(status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
  } else {
		/* set up some aliases for the returned arrays; note vector 0 is time */

  //  REAL8 *thet = yout->data;
		REAL8 *vphi = &yout->data[1*len]; REAL8 *omega = &yout->data[2*len];
		REAL8 *LNhx = &yout->data[3*len]; REAL8 *LNhy  = &yout->data[4*len];	REAL8 *LNhz  = &yout->data[5*len];

		/* these are not needed for the waveforms:
		REAL8 *S1x  = &yout->data[6*len]; REAL8 *S1y   = &yout->data[7*len];  REAL8 *S1z   = &yout->data[8*len];
		REAL8 *S2x  = &yout->data[9*len]; REAL8 *S2y   = &yout->data[10*len]; REAL8 *S2z   = &yout->data[11*len];	*/

		*countback = len;

		if (signalvec1) { /* return polarizations */
			REAL8 v=0, amp=0, alpha=0, alpha0 = atan2(LNhy[0],LNhx[0]);

			for(unsigned int i=0;i<len;i++) {
				v = pow(omega[i],oneby3);
				amp = params->signalAmplitude * (v*v);
				if(LNhx[i]*LNhx[i] + LNhy[i]*LNhy[i] > 0.0) {
          alpha = atan2(LNhy[i],LNhx[i]); alpha0 = alpha;
        } else {
          alpha = alpha0;
        }

				signalvec1->data[i]   = (REAL4)(-0.5 * amp * cos(2*vphi[i]) * cos(2*alpha) * (1.0 + LNhz[i]*LNhz[i]) \
				                                     + amp * sin(2*vphi[i]) * sin(2*alpha) * LNhz[i]);

				if (signalvec2) {
					signalvec2->data[i] = (REAL4)(-0.5 * amp * cos(2*vphi[i]) * sin(2*alpha) * (1.0 + LNhz[i]*LNhz[i]) \
																				     - amp * sin(2*vphi[i]) * cos(2*alpha) * LNhz[i]);
				}
			}

			params->fFinal = pow(v,3.0)/(LAL_PI*m);
			if (!signalvec2) params->tC = yout->data[len-1];	/* TO DO: why only in this case? */
		} else if (a) {	/* return coherentGW components */
			REAL8 apcommon, f2a, alpha, alpha0 = atan2(LNhy[0],LNhx[0]);

			/* (minus) amplitude for distance in m; should be (1e6 * LAL_PC_SI * params->distance) for distance in Mpc */
			apcommon = -4.0 * params->mu * LAL_MRSUN_SI/(params->distance);

			for(unsigned int i=0;i<len;i++) {
				f2a = pow(omega[i],twoby3);
				if(LNhx[i]*LNhx[i] + LNhy[i]*LNhy[i] > 0.0) {
          alpha = atan2(LNhy[i],LNhx[i]); alpha0 = alpha;
        } else {
          alpha = alpha0;
        }

			  ff   ->data[i]     = (REAL4)(omega[i]/unitHz);
			  a    ->data[2*i]   = (REAL4)(apcommon * f2a * 0.5 * (1 + LNhz[i]*LNhz[i]));
			  a    ->data[2*i+1] = (REAL4)(apcommon * f2a * LNhz[i]);
			  phi  ->data[i]     = (REAL8)(2.0 * vphi[i]);
			  shift->data[i]     = (REAL4)(2.0 * alpha);
			}

			params->fFinal = ff->data[len-1];
		}
	}

	if (yout) XLALDestroyREAL8Array(yout);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
