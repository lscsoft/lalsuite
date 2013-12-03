/*
*  Copyright (C) 2007 Stas Babak, David Churches, B.S. Sathyaprakash, Craig Robinson , Thomas Cokelaer
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

/**
 * \author Sathyaprakash, B. S.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief These modules generate a time-domain chirp waveform of type #TaylorT3.
 *
 * ### Prototypes ###
 *
 * <tt>LALInspiralWave3()</tt>
 * <ul>
 * <li> \c output: Output containing the inspiral waveform.</li>
 * <li> \c params: Input containing binary chirp parameters.</li>
 * </ul>
 *
 * <tt>LALInspiralWave3Templates()</tt>
 * <ul>
 * <li> \c output1: Output containing the 0-phase inspiral waveform.</li>
 * <li> \c output2: Output containing the \f$\pi/2\f$-phase inspiral waveform.</li>
 * <li> \c params: Input containing binary chirp parameters.</li>
 * </ul>
 *
 * ### Description ###
 *
 * LALInspiralWave3() generates #TaylorT3 approximant which
 * corresponds to the case wherein
 * the phase of the waveform is given as an explicit function of time
 * as in \eqref{eq_InspiralWavePhase3}.
 *
 * LALInspiralWave3Templates() simultaneously generates
 * two inspiral waveforms and the two differ in
 * phase by \f$\pi/2\f$.
 *
 * ### Algorithm ###
 *
 *
 * ### Uses ###
 *
 * \code
 * LALInspiralParameterCalc()
 * LALInspiralChooseModel()
 * LALInspiralSetup()
 * LALInspiralPhasing3 (via expnFunc)()
 * LALInspiralFrequency3 (via expnFunc)()
 * \endcode
 *
 */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>


typedef struct
{
	REAL8 (*func)(REAL8 tC, expnCoeffs *ak);
	expnCoeffs ak;
}
ChirptimeFromFreqIn;

static REAL8 XLALInspiralFrequency3Wrapper(REAL8 tC, void *pars);

static int
XLALInspiralWave3Engine(
                REAL4Vector      *output1,
                REAL4Vector      *output2,
                REAL4Vector      *h,
                REAL4Vector      *a,
                REAL4Vector      *ff,
                REAL8Vector      *phi,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                );

void
LALInspiralWave3 (
   LALStatus        *status,
   REAL4Vector      *output,
   InspiralTemplate *params
   )
{
   XLAL_PRINT_DEPRECATION_WARNING("XLALInspiralWave3");
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   if( XLALInspiralWave3(output, params) )
      ABORTXLAL(status);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

int
XLALInspiralWave3 (
   REAL4Vector      *output,
   InspiralTemplate *params
   )
{
  INT4 count;
  InspiralInit paramsInit;

  if (output == NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (output->data == NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (params == NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (params->nStartPad < 0)
    XLAL_ERROR(XLAL_EDOM);
  if (params->fLower <= 0)
    XLAL_ERROR(XLAL_EDOM);
  if (params->tSampling <= 0)
    XLAL_ERROR(XLAL_EDOM);

  if (XLALInspiralParameterCalc(params))
    XLAL_ERROR(XLAL_EFUNC);
  if (XLALInspiralSetup(&(paramsInit.ak), params))
    XLAL_ERROR(XLAL_EFUNC);
  if (XLALInspiralChooseModel(&(paramsInit.func), &(paramsInit.ak), params))
    XLAL_ERROR(XLAL_EFUNC);

  if (params->totalMass <= 0.)
    XLAL_ERROR(XLAL_EDOM);
  if (params->eta < 0.)
    XLAL_ERROR(XLAL_EDOM);

  memset( output->data, 0, output->length * sizeof(REAL4) );

  /* Call the engine function */
  count = XLALInspiralWave3Engine(output, NULL, NULL,
			NULL, NULL, NULL, params, &paramsInit);
  if (count < 0)
    XLAL_ERROR(XLAL_EFUNC);

  return XLAL_SUCCESS;
}

static REAL8 XLALInspiralFrequency3Wrapper(REAL8 tC, void *pars)
{
  ChirptimeFromFreqIn *in;
  REAL8 freq, f;

  in = (ChirptimeFromFreqIn *) pars;
  freq = in->func(tC, &(in->ak));
  if (XLAL_IS_REAL8_FAIL_NAN(freq))
    XLAL_ERROR_REAL8(XLAL_EFUNC);
  f = freq - in->ak.f0;

  /*
  fprintf(stderr, "Here freq=%e f=%e tc=%e f0=%e\n", freq, *f, tC, in->ak.f0);
   */

  return f;
}

void
LALInspiralWave3Templates (
   LALStatus        *status,
   REAL4Vector      *output1,
   REAL4Vector      *output2,
   InspiralTemplate *params
   )
{
   XLAL_PRINT_DEPRECATION_WARNING("XLALInspiralWave3Templates");
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   if( XLALInspiralWave3Templates(output1, output2, params) )
      ABORTXLAL(status);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

int
XLALInspiralWave3Templates (
   REAL4Vector      *output1,
   REAL4Vector      *output2,
   InspiralTemplate *params
   )
{
  INT4 count;

  InspiralInit paramsInit;

  if (output1 == NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (output2 == NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (output1->data == NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (output2->data == NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (params == NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (params->nStartPad < 0)
    XLAL_ERROR(XLAL_EDOM);
  if (params->fLower <= 0)
    XLAL_ERROR(XLAL_EDOM);
  if (params->tSampling <= 0)
    XLAL_ERROR(XLAL_EDOM);

  if (XLALInspiralParameterCalc(params))
    XLAL_ERROR(XLAL_EFUNC);
  if (XLALInspiralSetup(&(paramsInit.ak), params))
    XLAL_ERROR(XLAL_EFUNC);
  if (XLALInspiralChooseModel(&(paramsInit.func), &(paramsInit.ak), params))
    XLAL_ERROR(XLAL_EFUNC);

  if (params->totalMass <= 0.)
    XLAL_ERROR(XLAL_EDOM);
  if (params->eta < 0.)
    XLAL_ERROR(XLAL_EDOM);

  /* Initialise the waveforms to zero */
  memset(output1->data, 0, output1->length * sizeof(REAL4));
  memset(output2->data, 0, output2->length * sizeof(REAL4));

  /* Call the engine function */
  count = XLALInspiralWave3Engine(output1, output2, NULL,
			    NULL, NULL, NULL, params, &paramsInit);
  if (count < 0)
    XLAL_ERROR(XLAL_EFUNC);

  return XLAL_SUCCESS;
}

void
LALInspiralWave3ForInjection (
   LALStatus        *status,
   CoherentGW       *waveform,
   InspiralTemplate *params,
   PPNParamStruc  *ppnParams
   )
{
   XLAL_PRINT_DEPRECATION_WARNING("XLALInspiralWave3ForInjection");
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   if( XLALInspiralWave3ForInjection(waveform, params, ppnParams) )
      ABORTXLAL(status);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

int
XLALInspiralWave3ForInjection (
   CoherentGW       *waveform,
   InspiralTemplate *params,
   PPNParamStruc  *ppnParams
   )
{
  INT4 count;
  UINT4 i;
  REAL4Vector *h=NULL;
  REAL4Vector *a=NULL;
  REAL4Vector *ff=NULL ;
  REAL8Vector *phiv=NULL;

  REAL8 phiC;/* phase at coalescence */
  CHAR message[256];


  InspiralInit paramsInit;

 /** -- -- */

  /* Make sure parameter and waveform structures exist. */
  if (params == NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (waveform == NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (waveform->h != NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (waveform->a != NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (waveform->f != NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (waveform->phi != NULL)
    XLAL_ERROR(XLAL_EFAULT);
  if (waveform->shift != NULL)
    XLAL_ERROR(XLAL_EFAULT);

  params->ampOrder = (LALPNOrder) 0;
  sprintf(message, "WARNING: Amp Order has been reset to %d", params->ampOrder);
  XLALPrintInfo(message);
  /* Compute some parameters*/
  if (XLALInspiralInit(params, &paramsInit))
    XLAL_ERROR(XLAL_EFUNC);

  if (paramsInit.nbins==0)
  {
    /* FIXME: is this the correct thing to return? */
    return XLAL_SUCCESS;
  }

  /* Now we can allocate memory and vector for coherentGW structure*/
  ff = XLALCreateREAL4Vector(paramsInit.nbins);
  if (ff == NULL)
    XLAL_ERROR(XLAL_ENOMEM);
  a = XLALCreateREAL4Vector(2*paramsInit.nbins);
  if (a == NULL)
    XLAL_ERROR(XLAL_ENOMEM);
  phiv = XLALCreateREAL8Vector(paramsInit.nbins);
  if (phiv == NULL)
    XLAL_ERROR(XLAL_ENOMEM);

 /* By default the waveform is empty */

  memset(ff->data, 0, ff->length * sizeof(REAL4));
  memset(a->data,  0, a->length * sizeof(REAL4));
  memset(phiv->data, 0, phiv->length * sizeof(REAL8));

  if( params->ampOrder )
  {
    h = XLALCreateREAL4Vector(2*paramsInit.nbins);
    if (h == NULL)
      XLAL_ERROR(XLAL_ENOMEM);
    memset(h->data,  0, h->length * sizeof(REAL4));
  }

  /* Call the engine function */
  count = XLALInspiralWave3Engine(NULL, NULL, h, a, ff, phiv, params, &paramsInit);
  if (count < 0)
  {
    XLALDestroyREAL4Vector(ff);
    XLALDestroyREAL4Vector(a);
    XLALDestroyREAL8Vector(phiv);
    if( h )
    {
      XLALDestroyREAL4Vector(h);
    }
    XLAL_ERROR(XLAL_EFUNC);
  }

  /* Check an empty waveform hasn't been returned */
  for (i = 0; i < phiv->length; i++)
  {
    if (phiv->data[i] != 0.0) break;
    /* If the waveform returned is empty, return now */
    if (i == phiv->length - 1)
    {
      XLALDestroyREAL4Vector(ff);
      XLALDestroyREAL4Vector(a);
      XLALDestroyREAL8Vector(phiv);
      if( h )
      {
        XLALDestroyREAL4Vector(h);
      }

      /* FIXME: is this the correct thing to return? */
      return XLAL_SUCCESS;
    }
  }

  /*  if ( (phase/2./LAL_PI) < 2. ){
    sprintf(message, "The waveform has only %lf cycles; we don't keep waveform with less than 2 cycles.",
	       (double)phase/2./(double)LAL_PI );
    LALWarning(status, message);


  }
  else*/
  {

    /*wrap the phase vector*/
    phiC =  phiv->data[count-1] ;
    for (i=0; i<(UINT4)count;i++)
    {
      phiv->data[i] =  phiv->data[i] -phiC + ppnParams->phi;
    }
    /* Allocate the waveform structures. */
    waveform->a = (REAL4TimeVectorSeries *) XLALMalloc( sizeof(REAL4TimeVectorSeries) );
    if ( waveform->a == NULL )
      XLAL_ERROR(XLAL_ENOMEM);
    memset( waveform->a, 0, sizeof(REAL4TimeVectorSeries) );

    waveform->f = (REAL4TimeSeries *) LALMalloc( sizeof(REAL4TimeSeries) );
    if ( waveform->f == NULL )
    {
      XLALFree( waveform->a );
      waveform->a = NULL;
      XLAL_ERROR(XLAL_ENOMEM);
    }
    memset( waveform->f, 0, sizeof(REAL4TimeSeries) );

    waveform->phi = (REAL8TimeSeries *) LALMalloc( sizeof(REAL8TimeSeries) );
    if ( waveform->phi == NULL )
    {
      XLALFree( waveform->a );
      waveform->a = NULL;
      XLALFree( waveform->f );
      waveform->f = NULL;
      XLAL_ERROR(XLAL_ENOMEM);
    }
    memset( waveform->phi, 0, sizeof(REAL8TimeSeries) );

    waveform->a->data = XLALCreateREAL4VectorSequence(count, 2);
    if (waveform->a->data == NULL)
      XLAL_ERROR(XLAL_ENOMEM);
    waveform->f->data = XLALCreateREAL4Vector(count);
    if (waveform->f->data == NULL)
      XLAL_ERROR(XLAL_ENOMEM);
    waveform->phi->data = XLALCreateREAL8Vector(count);
    if (waveform->phi->data == NULL)
      XLAL_ERROR(XLAL_ENOMEM);

    memcpy(waveform->f->data->data , ff->data, count*(sizeof(REAL4)));
    memcpy(waveform->a->data->data , a->data, 2*count*(sizeof(REAL4)));
    memcpy(waveform->phi->data->data ,phiv->data, count*(sizeof(REAL8)));

    waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT = 1./params->tSampling;

    waveform->a->sampleUnits = lalStrainUnit;
    waveform->f->sampleUnits = lalHertzUnit;
    waveform->phi->sampleUnits = lalDimensionlessUnit;
    waveform->position = ppnParams->position;
    waveform->psi = ppnParams->psi;

    snprintf( waveform->a->name, LALNameLength, "T3 inspiral amplitudes" );
    snprintf( waveform->f->name, LALNameLength, "T3 inspiral frequency" );
    snprintf( waveform->phi->name, LALNameLength, "T3 inspiral phase" );

    /* --- fill some output ---*/

    ppnParams->tc     = (double)(count-1) / params->tSampling ;
    ppnParams->length = count;
    ppnParams->dfdt   = ((REAL4)(waveform->f->data->data[count-1] - waveform->f->data->data[count-2])) * ppnParams->deltaT;
    ppnParams->fStop  = params->fFinal;
    ppnParams->termCode        = GENERATEPPNINSPIRALH_EFSTOP;
    ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;

    ppnParams->fStart   = ppnParams->fStartIn;

    if( params->ampOrder )
    {
      waveform->h = (REAL4TimeVectorSeries *) XLALMalloc( sizeof(REAL4TimeVectorSeries) );
      if ( waveform->h == NULL )
        XLAL_ERROR(XLAL_ENOMEM);
      memset( waveform->h, 0, sizeof(REAL4TimeVectorSeries) );

      waveform->h->data = XLALCreateREAL4VectorSequence(count, 2);
      if ( waveform->h->data == NULL )
        XLAL_ERROR(XLAL_ENOMEM);
      memcpy(waveform->h->data->data , h->data, 2*count*(sizeof(REAL4)));
      waveform->h->deltaT = 1./params->tSampling;
      waveform->h->sampleUnits = lalStrainUnit;
      snprintf( waveform->h->name, LALNameLength, "T3 inspiral polarizations" );
      XLALDestroyREAL4Vector(h);
      h = NULL;
    }
  }    /*end of coherentGW storage */


  /* --- free memory --- */
  XLALDestroyREAL4Vector(ff);
  XLALDestroyREAL4Vector(a);
  XLALDestroyREAL8Vector(phiv);

  return XLAL_SUCCESS;
}


/* Engine function used to generate the waveforms */
static int
XLALInspiralWave3Engine(
                REAL4Vector      *output1,
                REAL4Vector      *output2,
                REAL4Vector      *h,
                REAL4Vector      *a,
                REAL4Vector      *ff,
                REAL8Vector      *phiv,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                )

{
  INT4 i, startShift, count;
  REAL8 dt, fu, eta, tc, totalMass, t, td, c1, phi0, phi1, phi;
  REAL8 v, v2, f, fHigh, tmax, fOld, phase, omega;
  REAL8 xmin,xmax,xacc;
  REAL8 (*frequencyFunction)(REAL8, void *);
  expnFunc func;
  expnCoeffs ak;
  ChirptimeFromFreqIn timeIn;
  void *pars;

  REAL8 temp, tempMax=0, tempMin = 0;

  /* Only used in injection case */
  REAL8 unitHz = 0.;
  REAL8 cosI = 0.;/* cosine of system inclination */
  REAL8 apFac = 0., acFac = 0.;/* extra factor in plus and cross amplitudes */


  ak   = paramsInit->ak;
  func = paramsInit->func;

  if (output2 || a)
      params->nStartPad = 0; /* must be zero for templates and injections */

  eta = params->eta;   /* Symmetric mass ratio  */
  /* Only in injection case.. */
  unitHz = params->totalMass*LAL_MTSUN_SI*(REAL8)LAL_PI;
  cosI   = cos( params->inclination );
  apFac  = acFac = -2.0 * eta *params->totalMass * LAL_MRSUN_SI/params->distance;
  apFac *= 1.0 + cosI*cosI;
  acFac *= 2.0 * cosI;

  dt = 1.0/(params->tSampling);    /* sampling rate  */
  fu = params->fCutoff;            /* upper frequency cutoff  */
  phi = params->startPhase;        /* initial phase  */
  startShift = params->nStartPad;  /* number of zeros at the start of the wave  */


  tc=params->tC;       /* Instant of coalescence of the compact objects */
  totalMass = (params->totalMass)*LAL_MTSUN_SI; /* mass of the system in seconds */


/*
   If flso is less than the user inputted upper frequency cutoff fu,
   then the waveforn is truncated at f=flso.
*/

  if (fu)
     fHigh = (fu < ak.flso) ? fu : ak.flso;
  else
     fHigh = ak.flso;

/*
   Check that the highest frequency is less than half the sampling frequency -
   the Nyquist theorum
*/

  if (fHigh >= 0.5/dt)
  {
    XLALPrintError("fHigh must be less than Nyquist frequency\n");
    XLAL_ERROR(XLAL_EDOM);
  }

  if (fHigh <= params->fLower)
  {
    XLALPrintError("fHigh must be larger than fLower\n");
    XLAL_ERROR(XLAL_EDOM);
  }

/* Here's the part which calculates the waveform */

  c1 = eta/(5.*totalMass);

  i = startShift;

  /*
   * In Jan 2003 we realized that the tC determined as a sum of chirp times is
   * not quite the tC that should enter the definition of Theta in the expression
   * for the frequency as a function of time (see DIS3, 2000). This is because
   * chirp times are obtained by inverting t(f). Rather tC should be obtained by
   * solving the equation f0 - f(tC) = 0. This is what is implemented below.
   */

  timeIn.func = func.frequency3;
  timeIn.ak = ak;
  frequencyFunction = &XLALInspiralFrequency3Wrapper;
  xmin = c1*params->tC/2.;
  xmax = c1*params->tC*2.;
  xacc = 1.e-6;
  pars = (void*) &timeIn;
  /* tc is the instant of coalescence */

  xmax = c1*params->tC*3 + 5.; /* we add 5 so that if tC is small then xmax
                                  is always greater than a given value (here 5)*/

  /* for x in [xmin, xmax], we search the value which gives the max frequency.
   and keep the corresponding rootIn.xmin. */



  for (tc = c1*params->tC/1000.; tc < xmax; tc+=c1*params->tC/1000.){
    temp = XLALInspiralFrequency3Wrapper(tc , pars);
    if (XLAL_IS_REAL8_FAIL_NAN(temp))
      XLAL_ERROR(XLAL_EFUNC);
    if (temp > tempMax) {
      xmin = tc;
      tempMax = temp;
    }
    if (temp < tempMin) {
      tempMin = temp;
    }
  }

  /* if we have found a value positive then everything should be fine in the
     BissectionFindRoot function */
  if (tempMax > 0  &&  tempMin < 0){
    tc = XLALDBisectionFindRoot (frequencyFunction, xmin, xmax, xacc, pars);
    if (XLAL_IS_REAL8_FAIL_NAN(tc))
      XLAL_ERROR(XLAL_EFUNC);
  }
  else if (a)
  {
    /* Otherwise we return an empty waveform for injection */
    return 0;
  }
  else
  {
    /* Or abort if not injection */
    XLALPrintError("Can't find good bracket for BisectionFindRoot\n");
    XLAL_ERROR(XLAL_EFAILED);
  }

  tc /= c1;

  tc += params->startTime;       /* Add user given startTime to instant of
                                     coalescence of the compact objects */

  t=0.0;
  td = c1*(tc-t);
  phase = func.phasing3(td, &ak);
  if (XLAL_IS_REAL8_FAIL_NAN(phase))
    XLAL_ERROR(XLAL_EFUNC);
  f = func.frequency3(td, &ak);
  if (XLAL_IS_REAL8_FAIL_NAN(f))
    XLAL_ERROR(XLAL_EFUNC);
  phi0=-phase+phi;
  phi1=phi0+LAL_PI_2;

  count = 0;
  tmax = tc - dt;
  fOld = 0.0;

/* We stop if any of the following conditions fail */

  while (f < fHigh && t < tmax && f > fOld)
  {
    /* Check we don't write past the end of the vector */
    if ((output1 && ((UINT4)i >= output1->length)) || (ff && ((UINT4)count >= ff->length)))
    {
      XLALPrintError("Attempting to write beyond the end of vector\n");
      XLAL_ERROR(XLAL_EBADLEN);
    }

    fOld = f;
    v = pow(f*LAL_PI*totalMass, (1./3.));
    v2 = v*v;
    if (output1)
    {


      /*
      output1->data[i]   = LALInspiralHPlusPolarization( phase+phi0, v, params );
      */
      output1->data[i]   = (REAL4) (apFac * v2) * cos(phase+phi0);
      if (output2)
      /*
        output2->data[i] = LALInspiralHCrossPolarization( phase+phi1, v, params );
      */
        output2->data[i] = (REAL4) (apFac * v2) * cos(phase+phi1);
    }
    else if (a)
    {
      int ice, ico;
      ice = 2*count;
      ico = ice + 1;
      omega = v*v*v;

      ff->data[count]     = (REAL4)(omega/unitHz);
      a->data[ice]        = (REAL4)(apFac * v2 );
      a->data[ico]        = (REAL4)(acFac * v2 );
      phiv->data[count]   = (REAL8)(phase );

      if (h)
      {
        h->data[ice] = LALInspiralHPlusPolarization( phase, v, params );
        h->data[ico] = LALInspiralHCrossPolarization( phase, v, params );
      }
    }
    ++i;
    ++count;
    t=count*dt;
    td = c1*(tc-t);
    phase = func.phasing3(td, &ak);
    if (XLAL_IS_REAL8_FAIL_NAN(phase))
      XLAL_ERROR(XLAL_EFUNC);
    f = func.frequency3(td, &ak);
    if (XLAL_IS_REAL8_FAIL_NAN(f))
      XLAL_ERROR(XLAL_EFUNC);
  }
  params->fFinal = fOld;
  if (output1 && !output2) params->tC = t;
/*
  fprintf(stderr, "%e %e\n", f, fHigh);
*/

  return count;
}
