#include <lal/LALSTPNWaveformFrameless.h>
#include <lal/LALSTPNWaveformErrors.h>
#include <lal/LALAdaptiveRungeKutta4.h>
#include <stdio.h>

#define UNUSED(expr) do { (void)(expr); } while (0)

typedef struct LALSTPNstructparams {
	REAL8 eta, m1m2, m2m1;
	REAL8 wdotnew, wdotorb[9], wspin15, wspin20;
	REAL8 LNhdot15, LNhdot20;
	REAL8 S1dot15, S2dot15, Sdot20;
	REAL8 epnorb[9];
} LALSTPNparams;

static int XLALSTPNAdaptiveSetParams(LALSTPNparams *mparams, 
    InspiralTemplate *params, InspiralInit *paramsInit) 
{
  // Check the relevant pointers
  if( !params || !mparams || !paramsInit )
    XLAL_ERROR(XLAL_EFAULT);

  /* Check the PN order is sane */
  /* LAL_PNORDER_NUM_ORDER is 1 more than the maximum PN order implemented */
  if( params->order > LAL_PNORDER_NUM_ORDER )
    XLAL_ERROR(XLAL_EINVAL);

  REAL8 m1, m2;
  m1 = params->mass1;
  m2 = params->mass2;
  mparams->m2m1 = m2 / m1;
  mparams->m1m2 = m1 / m2;

  mparams->eta = (m1 * m2) / (m1 + m2) / (m1 + m2);

  mparams->wdotnew = (96.0/5.0) * params->eta;

  for(int j = LAL_PNORDER_NEWTONIAN; j <= 8; j++) 
  {
    mparams->wdotorb[j] = mparams->epnorb[j] = 0.0;
  }

  mparams->epnorb[0] = 1.0;

  for(unsigned int j = LAL_PNORDER_NEWTONIAN; j <= params->order; j++) 
  {
    mparams->wdotorb[j] = paramsInit->ak.ST[j];
  }

  /* note: in the original code epnorb[2] is set even for */
  /* PNORDER_NEWTONIAN or PNORDER_HALF */
  if (params->order >= LAL_PNORDER_ONE) 
  {
    mparams->epnorb[2]  = -(1.0/12.0) * (9.0 + params->eta);
  }

  for(int j = params->order + 1; j <= 8; j++) 
  {
    mparams->wdotorb[j] = 0;
  }

  if (params->order >= LAL_PNORDER_ONE_POINT_FIVE) 
  {
    mparams->wspin15 	= -(1.0/12.0);
    mparams->LNhdot15 	= 0.5;
    mparams->S1dot15 	= (4.0 + 3.0 * mparams->m2m1) / 2.0 ;
    mparams->S2dot15 	= (4.0 + 3.0 * mparams->m1m2) / 2.0 ;
  } 
  else 
  {
    mparams->wspin15 	= 0.0;
    mparams->LNhdot15 	= 0.0;
    mparams->S1dot15 	= 0.0;
    mparams->S2dot15 	= 0.0;
  }

  if (params->order >= LAL_PNORDER_TWO) 
  {
    mparams->wspin20 	= -(1.0/48.0) / params->eta;
    mparams->LNhdot20 	= -1.5 / params->eta;
    mparams->Sdot20 	= 0.5;
    mparams->epnorb[4]  = (1.0/24.0) 
        * (-81.0 + 57.0*params->eta - params->eta*params->eta);
  } 
  else 
  {
    mparams->wspin20 	= 0.0;
    mparams->LNhdot20 	= 0.0;
    mparams->Sdot20 	= 0.0;
  }

  if (params->order >= LAL_PNORDER_THREE) 
  {
    mparams->epnorb[6]  = ( - (675.0/64.0) + ( (209323.0/4032.0) 
        -(205.0/96.0)*LAL_PI*LAL_PI - (110.0/9.0) * (-1987.0/3080.0) ) * params
        ->eta - (155.0/96.0)*(params->eta*params->eta) 
        - (35.0/5184.0)*(params->eta*params->eta*params->eta) );
  }

  if (params->order == LAL_PNORDER_THREE) 
  {
    mparams->wdotorb[(int)(LAL_PNORDER_THREE+1)] 
        = paramsInit->ak.ST[(int)(LAL_PNORDER_THREE+1)];
  }

  if (params->order == LAL_PNORDER_THREE_POINT_FIVE) 
  {
    mparams->wdotorb[8] = paramsInit->ak.ST[8];
  }

  return XLAL_SUCCESS;
}

static int XLALSTPNAdaptiveTest(double t, const double values[],
    double dvalues[],void *mparams) 
{
  REAL8 omega, v, test;
  LALSTPNparams *params = (LALSTPNparams*)mparams;

  UNUSED(t);

  omega = values[1];
  v = pow(omega,(1./3.));

  test = -0.5 * params->eta * ( (2.0/3.0) * (1.0/v) * params->epnorb[0]
      + params->epnorb[1] + (4.0/3.0) * v * (params->epnorb[2]
      + (5.0/4.0) * v * (params->epnorb[3] + (6.0/5.0) * v * (params->epnorb[4]
      + (7.0/6.0) * v * (params->epnorb[5] + (8.0/7.0) * v * (params->epnorb[6]
      + (9.0/8.0) * v * (params->epnorb[7] 
      + (10.0/9.0)* v *  params->epnorb[8] )))))) );

  if (params->wspin15 != 0.0) 
  {
    REAL8 LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z;
    REAL8 v2, dotLNS1, dotLNS2;

    LNhx = values[2]; LNhy  = values[3]; LNhz = values[4] ;
    S1x  = values[5]; S1y   = values[6]; S1z  = values[7] ;
    S2x  = values[8]; S2y   = values[9]; S2z  = values[10];

    v2 = v * v;

    dotLNS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
    dotLNS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);

    test += -0.5 * params->eta * (5.0/3.0) * v2 
        * ( (8.0/3.0 + 2.0*params->m2m1)*dotLNS1 
        + (8.0/3.0 + 2.0*params->m1m2)*dotLNS2 );

    if (params->wspin20 != 0.0) 
    {
      REAL8 dotS1S2;

      dotS1S2 = (S1x*S2x + S1y*S2y + S1z*S2z);
      test += -(v*v2)  * (dotS1S2 - 3.0 * dotLNS1 * dotLNS2);
    }
  }

  if (test > 0.0) /* energy test fails! */
  {
    return LALSTPN_TEST_ENERGY;
  } 
  else if (dvalues[1] < 0.0) /* omegadot < 0! */
  {
    return LALSTPN_TEST_OMEGADOT;
  } 
  else if isnan(omega) /* omega is nan! */
  {
    return LALSTPN_TEST_OMEGANAN;
  } 
  else /* Step is valid, continue the integration */
  {
    return GSL_SUCCESS;
  }
}

static int XLALSTPNFramelessAdaptiveDerivatives(double t, 
    const double values[], double dvalues[], void *mparams) 
{
  /* coordinates and derivatives */
  REAL8 omega, LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z, E1x, E1y, E1z;
  REAL8 ds, domega, dLNhx, dLNhy, dLNhz, dS1x, dS1y, dS1z, dS2x, dS2y, dS2z; 
  REAL8 dE1x, dE1y, dE1z;

  /* auxiliary variables */
  REAL8 v, v2, v3, v4, v7, v11;
  REAL8 dotLNS1, dotLNS2, dotS1S2;
  REAL8 omega2, tmpx, tmpy, tmpz, tmp1;
  REAL8 LNmag, crossx, crossy, crossz;

  LALSTPNparams *params = (LALSTPNparams*)mparams;

  UNUSED(t);

  /* copy variables */
  // UNUSED!!: s    = values[0];
  omega = values[1];
  LNhx = values[2];  LNhy  = values[3];  LNhz = values[4] ;
  S1x  = values[5];  S1y   = values[6];  S1z  = values[7] ;
  S2x  = values[8];  S2y   = values[9];  S2z  = values[10];
  E1x  = values[11]; E1y   = values[12]; E1z  = values[13];

  if (omega <= 0.0) 
  {
    return LALSTPN_DERIVATIVE_OMEGANONPOS;
  }

  v  = pow(omega,(1./3.));
  v2 = v * v; v3 = v2 * v; v4 = v3 * v; v7 = v4 * v3; v11 = v7 * v4;

  dotLNS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
  dotLNS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);
  dotS1S2 = (S1x*S2x  + S1y*S2y  + S1z*S2z );

  /* domega - non-spinning PN corrections... */
  domega = params->wdotorb[0] + v * (params->wdotorb[1] + 
      v * (params->wdotorb[2] + v * (params->wdotorb[3] + 
      v * (params->wdotorb[4] + v * (params->wdotorb[5] + 
      v * (params->wdotorb[6] + params->wdotorb[7] * log(omega) + 
      v * params->wdotorb[8] ))))));

  /* ...now add in spin-dependent corrections... */
  domega += params->wspin15 * omega * ( LNhx * (113.0 * S1x + 113.0 * S2x 
      + 75.0 * params->m2m1 * S1x + 75.0 * params->m1m2 * S2x) 
      + LNhy * (113.0 * S1y + 113.0 * S2y + 75.0 * params->m2m1 * S1y 
      + 75.0 * params->m1m2 * S2y) + LNhz * (113.0 * S1z + 113.0 * S2z 
      + 75.0 * params->m2m1 * S1z + 75.0 * params->m1m2 * S2z) );

  domega += params->wspin20 * v4 *(247.0 * dotS1S2 - 721.0 * dotLNS1 * dotLNS2);

  domega *= params->wdotnew * v11; /* ... now multiply by Newtonian term */

  /* dLN */
  omega2 = omega * omega;

  tmpx = params->LNhdot15 * omega2 * ((4.0 + 3.0*params->m2m1) * S1x 
      + (4.0 + 3.0*params->m1m2) * S2x);
  tmpy = params->LNhdot15 * omega2 * ((4.0 + 3.0*params->m2m1) * S1y 
      + (4.0 + 3.0*params->m1m2) * S2y);
  tmpz = params->LNhdot15 * omega2 * ((4.0 + 3.0*params->m2m1) * S1z 
      + (4.0 + 3.0*params->m1m2) * S2z);

  tmpx += params->LNhdot20 * v7 * (dotLNS2 * S1x + dotLNS1 * S2x);
  tmpy += params->LNhdot20 * v7 * (dotLNS2 * S1y + dotLNS1 * S2y);
  tmpz += params->LNhdot20 * v7 * (dotLNS2 * S1z + dotLNS1 * S2z);

  dLNhx = (-tmpz*LNhy + tmpy*LNhz);
  dLNhy = (-tmpx*LNhz + tmpz*LNhx);
  dLNhz = (-tmpy*LNhx + tmpx*LNhy);

  /* dE1; reuses tmpxyz above, which is PBCV's Omega_L */

  tmp1 = tmpx * LNhx + tmpy * LNhy + tmpz * LNhz;
  tmpx -= tmp1 * LNhx;
  tmpy -= tmp1 * LNhy;
  tmpz -= tmp1 * LNhz; /* Now tmpxyz is PBCV's Omega_e */

  dE1x = (-tmpz*E1y + tmpy*E1z);
  dE1y = (-tmpx*E1z + tmpz*E1x);
  dE1z = (-tmpy*E1x + tmpx*E1y);

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
  /* No coordinate singularity; instead, frame is evolving */
  ds = omega;

  dvalues[0]  = ds   ; dvalues[1]  = domega;
  dvalues[2]  = dLNhx; dvalues[3]  = dLNhy ; dvalues[4]  = dLNhz;
  dvalues[5]  = dS1x ; dvalues[6]  = dS1y  ; dvalues[7]  = dS1z ;
  dvalues[8]  = dS2x ; dvalues[9]  = dS2y  ; dvalues[10] = dS2z ;
  dvalues[11] = dE1x ; dvalues[12] = dE1y  ; dvalues[13] = dE1z ;

  return GSL_SUCCESS;
}

void LALSTPNFramelessWaveform(LALStatus *status, REAL4Vector *signalvec, 
    InspiralTemplate *params)
{
  XLALPrintDeprecationWarning("LALSTPNFramelessWaveform",
      "XLALSTPNFramelessWaveform");
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if( XLALSTPNFramelessWaveform(signalvec, params) )
      ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

int XLALSTPNFramelessWaveform(REAL4Vector *signalvec, InspiralTemplate *params)
{
  // Check the relevant pointers
  if( !signalvec || !signalvec->data || !params )
    XLAL_ERROR(XLAL_EFAULT);

  // Check the parameters are sane
  if( params->nStartPad < 0 || params->nEndPad < 0 || params->fLower <= 0
      || params->tSampling <= 0 || params->totalMass <= 0.)
    XLAL_ERROR(XLAL_EINVAL);

  UINT4 count;
  InspiralInit paramsInit;

  if( XLALInspiralSetup(&(paramsInit.ak), params) )
    XLAL_ERROR(XLAL_EFUNC);
  if( XLALInspiralChooseModel(&(paramsInit.func), &(paramsInit.ak), params) )
    XLAL_ERROR(XLAL_EFUNC);

  memset(signalvec->data, 0, signalvec->length * sizeof( REAL4 ));   
  /* Call the engine function */
  if( XLALSTPNFramelessAdaptiveWaveformEngine(signalvec, NULL, &count, params, 
      &paramsInit) )
    XLAL_ERROR(XLAL_EFUNC);

  return XLAL_SUCCESS;    
}

void LALSTPNFramelessWaveformTemplates(LALStatus *status, 
    REAL4Vector *signalvec1, REAL4Vector *signalvec2, InspiralTemplate *params)
{
  XLALPrintDeprecationWarning("LALSTPNFramelessWaveformTemplates",
      "XLALSTPNFramelessWaveformTemplates");
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if( XLALSTPNFramelessWaveformTemplates(signalvec1, signalvec2, params) )
      ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}



int XLALSTPNFramelessWaveformTemplates(REAL4Vector *signalvec1, 
    REAL4Vector *signalvec2, InspiralTemplate *params)
{
  // Check the relevant pointers
  if( !signalvec1 || !signalvec1->data || !signalvec2 
      || !signalvec2->data || !params )
    XLAL_ERROR(XLAL_EFAULT);

  // Check the parameters are sane
  if( params->nStartPad < 0 || params->nEndPad < 0 || params->fLower <= 0
      || params->tSampling <= 0 || params->totalMass <= 0.)
    XLAL_ERROR(XLAL_EINVAL);


  UINT4 count;
  InspiralInit paramsInit;

  if( XLALInspiralInit(params, &paramsInit) )
    XLAL_ERROR(XLAL_EFUNC);

  memset(signalvec1->data, 0, signalvec1->length * sizeof( REAL4 ));   
  memset(signalvec2->data, 0, signalvec2->length * sizeof( REAL4 ));   
  /* Call the engine function */
  if( XLALSTPNFramelessAdaptiveWaveformEngine(signalvec1, signalvec2, 
      &count, params, &paramsInit) )
    XLAL_ERROR(XLAL_EFUNC);

  return XLAL_SUCCESS;
}

void LALSTPNFramelessWaveformForInjection(LALStatus *status, 
    CoherentGW *waveform, InspiralTemplate *params, PPNParamStruc *ppnParams)
{
  XLALPrintDeprecationWarning("LALSTPNFramelessWaveformForInjection",
      "XLALSTPNFramelessWaveformForInjection");
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if( XLALSTPNFramelessWaveformForInjection(waveform, params, ppnParams) )
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

int XLALSTPNFramelessWaveformForInjection(CoherentGW *waveform, 
    InspiralTemplate *params, PPNParamStruc *ppnParams)
{
  // Check the relevant pointers exist
  if( !waveform || !params )
    XLAL_ERROR(XLAL_EFAULT);

  // CoherentGW waveform->h shouldn't exist yet
  if( waveform->h )
    XLAL_ERROR(XLAL_EFAULT);

  UINT4 count, i;

  REAL4Vector *hplus  = NULL;
  REAL4Vector *hcross = NULL;

  InspiralInit paramsInit;

  /* Compute some parameters*/

  /* The start phase is passed via coa_phase in the SimInspiralTable, */
  /* then ppnParams->phi due to convention in GenerateInspiral.c */
  params->startPhase = ppnParams->phi; 
  if( XLALInspiralInit(params, &paramsInit) )
    XLAL_ERROR(XLAL_EFUNC);
  if (paramsInit.nbins==0)
    {
      XLALPrintWarning("Warning! Waveform of zero length requested in %s. Returning empty waveform.\n ", __func__);
      return XLAL_SUCCESS;
    }

  /* Now we can allocate memory and vector for coherentGW structure*/
  hplus  = XLALCreateREAL4Vector(paramsInit.nbins);
  hcross = XLALCreateREAL4Vector(paramsInit.nbins);

  /* By default the waveform is empty */
  memset(hplus->data,  0, paramsInit.nbins * sizeof(REAL4));
  memset(hcross->data, 0, paramsInit.nbins * sizeof(REAL4));


  /* Call the engine function */
  if( XLALSTPNFramelessAdaptiveWaveformEngine(hplus, hcross, &count, 
      params, &paramsInit) )
    XLAL_ERROR(XLAL_EFUNC);

  /* Check an empty waveform hasn't been returned
     Is this something that we should actually check for?
     If so, uncomment and switch to checking hplus
  for (i = 0; i < a->length; i++)
  {
    if (a->data[i] != 0.0) break;
    if (i == a->length - 1)
    {
      XLALDestroyREAL4Vector(ff);
      XLALDestroyREAL4Vector(a);
      XLALDestroyREAL4Vector(phi);
      XLALDestroyREAL4Vector(shift);
    }
  } */


  /* Allocate the waveform structures. */
  if ( ( waveform->h = (REAL4TimeVectorSeries *) 
      XLALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) 
  {
    XLAL_ERROR(XLAL_ENOMEM);
  }
  memset( waveform->h, 0, sizeof(REAL4TimeVectorSeries) );

  // UNUSED!! CreateVectorSequenceIn in;
  // in.length = (UINT4)count;
  // in.vectorLength = 2;
  // '2' is the number of vectors - 1 for each of h+ and hx
  waveform->h->data = XLALCreateREAL4VectorSequence( (UINT4) count, 2);

  /* We alternate hplus and hcross. LALInjectStrainGW actually wants */
  /* a continuous block of hplus followed by a continuous block of hcross */
  /* But, LALFindChirpInjectSignals reorders h assuming it is passed */
  /* into the function as alternating hplus and hcross */
  for( i = 0; i < count; i++)
  {
    waveform->h->data->data[2*i]   = hplus->data[i];
    waveform->h->data->data[2*i+1] = hcross->data[i];
  }

  waveform->h->deltaT = 1./params->tSampling;
  waveform->h->sampleUnits = lalStrainUnit;

  waveform->position = ppnParams->position;
  waveform->psi = ppnParams->psi;


  snprintf(waveform->h->name, LALNameLength, "STPN inspiral strain");


  /* --- fill some output ---*/
  ppnParams->tc     = (double)(count-1) / params->tSampling ;
  ppnParams->length = count;
  /*ppnParams->dfdt = ((REAL4)(waveform->f->data->data[count-1]
	- waveform->f->data->data[count-2])) * ppnParams->deltaT; */
  ppnParams->fStop  = params->fFinal;
  ppnParams->termCode        = GENERATEPPNINSPIRALH_EFSTOP;
  ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;

  ppnParams->fStart = ppnParams->fStartIn;

  /* --- free memory --- */
  XLALDestroyREAL4Vector(hplus);
  XLALDestroyREAL4Vector(hcross);

  return XLAL_SUCCESS;
}

void LALSTPNFramelessAdaptiveWaveformEngine( LALStatus *status, 
    REAL4Vector *signalvec1, REAL4Vector *signalvec2, UINT4 *countback, 
    InspiralTemplate *params, InspiralInit *paramsInit)
{
  XLALPrintDeprecationWarning("LALSTPNFramelessAdaptiveWaveformEngine",
      "XLALSTPNFramelessAdaptiveWaveformEngine");

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if( XLALSTPNFramelessAdaptiveWaveformEngine(signalvec1, signalvec2, 
      countback, params, paramsInit) )
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


int XLALSTPNFramelessAdaptiveWaveformEngine(REAL4Vector *signalvec1, 
    REAL4Vector *signalvec2, UINT4 *countback, InspiralTemplate *params, 
    InspiralInit *paramsInit)
{
  /* PN parameters */
  LALSTPNparams mparams;

  /* needed for integration */
  ark4GSLIntegrator *integrator;
  unsigned int len;
  int intreturn;
  REAL8 yinit[14];
  REAL8Array *yout;

  /* other computed values */
  REAL8 unitHz, dt, m, lengths, norm;
  REAL8 E2x, E2y;
  REAL8 hpluscos, hplussin, hcrosscos, hcrosssin;

  /* units */
  unitHz = params->totalMass * LAL_MTSUN_SI * (REAL8)LAL_PI;
  dt = 1.0/params->tSampling;	  /* tSampling is in Hz, so dt is in seconds */
  m = params->totalMass * LAL_MTSUN_SI;

  /* length estimation (Newtonian) */ 
  /* since integration is adaptive, we could use a better estimate */
  lengths = (5.0/256.0) * pow(LAL_PI,-8.0/3.0) * pow(params->chirpMass 
      * LAL_MTSUN_SI * params->fLower,-5.0/3.0) / params->fLower;

  /* setup coefficients for PN equations */
  XLALSTPNAdaptiveSetParams(&mparams, params, paramsInit);

  /* initialize the coordinates */
  yinit[0] = params->startPhase / 2.0; /* initial vphi sets start phase */
  yinit[1] = params->fLower * unitHz;  /* omega (really pi M f) */

  yinit[2] = sin(params->inclination); /* LNh(x,y,z) */
  yinit[3] = 0.0;
  yinit[4] = cos(params->inclination);

  norm = pow(params->mass1/params->totalMass,2.0);
  yinit[5] = norm * params->spin1[0];  /* S1(x,y,z) */
  yinit[6] = norm * params->spin1[1];
  yinit[7] = norm * params->spin1[2];

  norm = pow(params->mass2/params->totalMass,2.0);
  yinit[8] = norm * params->spin2[0];  /* S2(x,y,z) */
  yinit[9] = norm * params->spin2[1];
  yinit[10]= norm * params->spin2[2];

  yinit[11] = cos(params->inclination); /* E1(x,y,z) */
  yinit[12] = 0.0;
  yinit[13] = -1.0 * sin(params->inclination);

  xlalErrno = 0;

  /* allocate the integrator */
  integrator = XLALAdaptiveRungeKutta4Init(14, 
      XLALSTPNFramelessAdaptiveDerivatives, XLALSTPNAdaptiveTest, 
      1.0e-6, 1.0e-6);
  if (!integrator) 
  {
    XLALPrintError("LALSTPNFramelessAdaptiveWaveformEngine: Cannot allocate integrator.\n");
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* stop the integration only when the test is true */
  integrator->stopontestonly = 1;

  /* run the integration; note: time is measured in units of total mass */
  len = XLALAdaptiveRungeKutta4(integrator, (void *) &mparams, yinit, 0.0, 
      lengths/m, dt/m, &yout);

  intreturn = integrator->returncode;
  XLALAdaptiveRungeKutta4Free(integrator);

  if (!len) 
  {
    XLALPrintError("LALSTPNFramelessAdaptiveWaveformEngine: integration failed with errorcode %d.\n",intreturn);
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* report on abnormal termination (TO DO: throw some kind of LAL error?) */
  if (intreturn != 0 && intreturn != LALSTPN_TEST_ENERGY 
      && intreturn != LALSTPN_TEST_OMEGADOT) 
  {
    XLALPrintWarning("LALSTPNFramelessAdaptiveWaveformEngine WARNING: integration terminated with code %d.\n",intreturn);
    XLALPrintWarning("Waveform parameters were m1 = %e, m2 = %e, s1 = (%e,%e,%e), s2 = (%e,%e,%e), inc = %e.\n", 
        params->mass1, params->mass2, params->spin1[0], params->spin1[1], 
        params->spin1[2], params->spin2[0], params->spin2[1], 
        params->spin2[2], params->inclination);
  }

  /* check that we're not above Nyquist. Need to find a way to pass 
   * this information along without creating GB worth of warnings. 
   * a --no-warning flag?
   * if (yinit[1]/unitHz > 0.5 * params->tSampling) {
   *   fprintf(stderr,"LALSTPNWaveform2 WARNING: final frequency above Nyquist.\n");
   *  }
   */
  /* if we have enough space, compute the waveform components; otherwise abort*/
  if (signalvec1 && len >= signalvec1->length) 
  {
    XLALPrintError("LALSTPNFramelessAdaptiveWaveformEngine: no space to write in signalvec1: %d vs. %d\n", len, signalvec1->length);
    XLAL_ERROR(XLAL_EBADLEN);
  }
  if (signalvec2 && len >= signalvec2->length) 
  {
    XLALPrintError("LALSTPNFramelessAdaptiveWaveformEngine: no space to write in signalvec2: %d vs. %d\n", len, signalvec2->length);
    XLAL_ERROR(XLAL_EBADLEN);
  }

  /* set up some aliases for the returned arrays; note vector 0 is time */
  REAL8 *vphi = &yout->data[1*len]; REAL8 *omega = &yout->data[2*len];
  REAL8 *LNhx = &yout->data[3*len]; REAL8 *LNhy  = &yout->data[4*len]; 
  REAL8 *LNhz = &yout->data[5*len];

  /*-- these are not needed for the waveforms, but uncomment for debugging: --*/
  /*REAL8 *S1x  = &yout->data[6*len];  REAL8 *S1y   = &yout->data[7*len];  
  REAL8 *S1z  = &yout->data[8*len];  REAL8 *S2x   = &yout->data[9*len]; 
  REAL8 *S2y  = &yout->data[10*len]; REAL8 *S2z   = &yout->data[11*len];
  REAL8 *that = &yout->data[0];*/
  /*--------------------------------------------------------------------------*/

  REAL8 *E1x = &yout->data[12*len]; REAL8 *E1y = &yout->data[13*len]; 
  REAL8 *E1z = &yout->data[14*len];

  *countback = len;

  REAL8 amp, f2a;
  amp = -4.0 * params->mu * LAL_MRSUN_SI / (params->distance);
  /*---- Uncomment the next 2 lines for debugging ----*/
  /*FILE *outFile=NULL;
  outFile = fopen("frameless-dynamics.dat", "w");*/
  /*--------------------------------------------------*/
  for(unsigned int i=0;i<len;i++) 
  {
    f2a = pow(omega[i],2.0/3.0);

    E2x = LNhy[i]*E1z[i] - LNhz[i]*E1y[i]; /* E2 = LNhat x E1 */
    E2y = LNhz[i]*E1x[i] - LNhx[i]*E1z[i];

    /*---- Uncomment the next line for debugging ----*/
    /*    E2z = LNhx[i]*E1y[i] - LNhy[i]*E1x[i];
          fprintf(outFile,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
		that[i] * m, vphi[i], omega[i], LNhx[i], LNhy[i], LNhz[i], 
		S1x[i], S1y[i], S1z[i], S2x[i], S2y[i], S2z[i], 
		E1x[i], E1y[i], E1z[i], E2x, E2y, E2z);*/
    /* ----------------------------------------------*/

    if (signalvec1) 
    {
      /* Polarization tensors projected into viewer's frame */
      hpluscos  = 0.5 * (E1x[i]*E1x[i] - E1y[i]*E1y[i] - E2x*E2x + E2y*E2y);
      hplussin  = E1x[i]*E2x - E1y[i]*E2y;

      signalvec1->data[i] = (REAL4) ( amp * f2a * \
          ( hpluscos * cos(2*vphi[i]) + hplussin * sin(2*vphi[i]) ) );
    }

    if (signalvec2) 
    {
      /* Polarization tensors projected into viewer's frame */
      hcrosscos = E1x[i]*E1y[i] - E2x*E2y;
      hcrosssin = E1y[i]*E2x + E1x[i]*E2y;

      signalvec2->data[i] = (REAL4) ( amp * f2a * \
          ( hcrosscos * cos(2*vphi[i]) + hcrosssin * sin(2*vphi[i]) ) );
    }
  }
  /*---- Uncomment next line for debugging ----*/
  /*fclose(outFile);*/
  /*-------------------------------------------*/
      
  params->fFinal = (REAL4)(omega[len-1]/unitHz);
  params->tC = yout->data[len-1] * m; /* m converts from t/m to t in seconds */

  if (yout) XLALDestroyREAL8Array(yout);

  return XLAL_SUCCESS;
}

