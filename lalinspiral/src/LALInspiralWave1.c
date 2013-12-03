/*
*  Copyright (C) 2007 Stas Babak, David Churches, Alexander Dietz, B.S. Sathyaprakash, Craig Robinson , Thomas Cokelaer
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
 * \brief The code \c LALInspiralWave1() generates an time-domain inspiral waveform corresponding to the
 * ::Approximant #TaylorT1 and #PadeT1 as outlined in the documentation for the function \c LALInspiralWave().
 *
 * ### Prototypes ###
 *
 * <tt>LALInspiralWave1()</tt>
 * <ul>
 * <li> \c signalvec: Output containing the inspiral waveform.</li>
 * <li> \c params: Input containing binary chirp parameters.</li>
 * </ul>
 *
 * <tt>LALInspiralWave1Templates()</tt>
 * <ul>
 * <li> \c signalvec1: Output containing the 0-phase inspiral waveform.</li>
 * <li> \c signalvec2: Output containing the \f$\pi/2\f$-phase inspiral waveform.</li>
 * <li> \c params: Input containing binary chirp parameters.</li>
 * </ul>
 *
 * ### Description ###
 *
 * LALInspiralWave1() is called if the user has specified the
 * \c enum ::Approximant to be
 * either #TaylorT1 or #PadeT1.
 * LALInspiralWave1Templates() is exactly the same as LALInspiralWave1(), except that
 * it generates two templates one for which the starting phase is
 * <tt>params.startPhase</tt> and the other for which the phase is
 * <tt>params.startPhase + \f$\pi/2\f$</tt>.
 *
 * ### Algorithm ###
 *
 * This code uses a fourth-order Runge-Kutta algorithm to solve the ODEs
 * in \eqref{eq_ode2}.
 *
 * ### Uses ###
 *
 * \code
 * LALInspiralSetup()
 * LALInspiralChooseModel()
 * LALInspiralVelocity()
 * LALInspiralPhasing1()
 * LALInspiralDerivatives()
 * LALRungeKutta4()
 * \endcode
 *
 * ### Notes ###
 *
 */

/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>

static int
XLALInspiralWave1Engine(
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
LALInspiralWave1(
   LALStatus        *status,
   REAL4Vector      *output,
   InspiralTemplate *params
   )
{
   XLAL_PRINT_DEPRECATION_WARNING("XLALInspiralWave1");
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   if( XLALInspiralWave1(output, params) )
      ABORTXLAL(status);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

int
XLALInspiralWave1(
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

   /*Call the engine function*/
   count = XLALInspiralWave1Engine(output, NULL, NULL, NULL, NULL, NULL, params, &paramsInit);
   if (count < 0)
      XLAL_ERROR(XLAL_EFUNC);

   return XLAL_SUCCESS;
}



/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/

void
LALInspiralWave1Templates(
   LALStatus        *status,
   REAL4Vector      *output1,
   REAL4Vector      *output2,
   InspiralTemplate *params
   )
{
   XLAL_PRINT_DEPRECATION_WARNING("XLALInspiralWave1Templates");
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   if( XLALInspiralWave1Templates(output1, output2, params) )
      ABORTXLAL(status);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

int
XLALInspiralWave1Templates(
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
   count = XLALInspiralWave1Engine(output1, output2, NULL, NULL, NULL, NULL, params, &paramsInit);
   if (count < 0)
      XLAL_ERROR(XLAL_EFUNC);

   return XLAL_SUCCESS;
}

/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms for injection packages T.Cokelaer sept 2003
*/

void
LALInspiralWave1ForInjection(
  LALStatus        *status,
  CoherentGW       *waveform,
  InspiralTemplate *params,
  PPNParamStruc  *ppnParams
  )
{
  XLAL_PRINT_DEPRECATION_WARNING("XLALInspiralWave1ForInjection");
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if( XLALInspiralWave1ForInjection(waveform, params, ppnParams) )
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

int
XLALInspiralWave1ForInjection(
  CoherentGW       *waveform,
  InspiralTemplate *params,
  PPNParamStruc  *ppnParams
  )
{
  INT4        count, i;
  REAL8       p, phiC;

  REAL4Vector *a   = NULL;      /* pointers to generated amplitude  data */
  REAL4Vector *h   = NULL;      /* pointers to generated polarizations */
  REAL4Vector *ff  = NULL;      /* pointers to generated  frequency data */
  REAL8Vector *phi = NULL;      /* pointer to generated phase data */

  CHAR message[256];

  InspiralInit paramsInit;


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

  if (paramsInit.nbins == 0)
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
  phi = XLALCreateREAL8Vector(paramsInit.nbins);
  if (phi == NULL)
    XLAL_ERROR(XLAL_ENOMEM);

  /* By default the waveform is empty */
  memset(ff->data, 0, paramsInit.nbins * sizeof(REAL4));
  memset(a->data, 0, 2 * paramsInit.nbins * sizeof(REAL4));
  memset(phi->data, 0, paramsInit.nbins * sizeof(REAL8));

  if( params->ampOrder )
  {
    h = XLALCreateREAL4Vector(2*paramsInit.nbins);
    if (h == NULL)
      XLAL_ERROR(XLAL_ENOMEM);
    memset(h->data,  0, h->length * sizeof(REAL4));
  }

  /* Call the engine function */
  count = XLALInspiralWave1Engine(NULL, NULL, h, a, ff,
			     phi, params, &paramsInit);
  if (count < 0)
  {
    XLALDestroyREAL4Vector(ff);
    XLALDestroyREAL4Vector(a);
    XLALDestroyREAL8Vector(phi);
    if( h )
    {
      XLALDestroyREAL4Vector(h);
    }
    XLAL_ERROR(XLAL_EFUNC);
  }

  p = phi->data[count-1];

  params->fFinal = ff->data[count-1];
  sprintf(message, "cycles = %f", fabs(p)/(double)LAL_TWOPI);
  XLALPrintInfo(message);

  if ( (INT4)(p/LAL_TWOPI) < 2 ){
    sprintf(message, "The waveform has only %f cycles; we don't keep waveform with less than 2 cycles.",
	       fabs(p)/(double)LAL_TWOPI );
    XLALPrintWarning(message);
  }
  else
  {
    /*wrap the phase vector*/
    phiC =  phi->data[count-1] ;
    for (i=0; i<count;i++)
    {
      phi->data[i] =  phi->data[i] -phiC + ppnParams->phi;
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
    memcpy(waveform->phi->data->data ,phi->data, count*(sizeof(REAL8)));

    waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT = ppnParams->deltaT;

    waveform->a->sampleUnits = lalStrainUnit;
    waveform->f->sampleUnits = lalHertzUnit;
    waveform->phi->sampleUnits = lalDimensionlessUnit;
    waveform->position = ppnParams->position;
    waveform->psi = ppnParams->psi;

    snprintf( waveform->a->name, LALNameLength,   "T1 inspiral amplitude" );
    snprintf( waveform->f->name, LALNameLength,   "T1 inspiral frequency" );
    snprintf( waveform->phi->name, LALNameLength, "T1 inspiral phase" );

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
  XLALDestroyREAL8Vector(phi);

  return XLAL_SUCCESS;
}

/*
 *  Engine function for use by other LALInspiralWave1* functions
 *  Craig Robinson April 2005
 */

int
XLALInspiralWave1Engine(
		REAL4Vector      *signalvec1,
		REAL4Vector      *signalvec2,
		REAL4Vector      *h,
		REAL4Vector      *a,
		REAL4Vector      *ff,
		REAL8Vector      *phi,
		InspiralTemplate *params,
		InspiralInit     *paramsInit
)
{
   INT4 n=2, count;
   REAL8 omega;
   REAL8 amp, m, dt, t, v, p, h1, h2, f, fu, fHigh, piM;
   REAL8Vector dummy, values, dvalues, valuesNew, yt, dym, dyt;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   rk4GSLIntegrator *integrator;
   void *funcParams;
   expnCoeffs ak;
   expnFunc func;

   REAL8 mTot = 0;
   REAL8 unitHz = 0;
   REAL8 f2a = 0;
   REAL8 mu = 0;
   REAL8 cosI = 0;/* cosine of system inclination */
   REAL8 etab = 0;
   REAL8 fFac = 0; /* SI normalization for f and t */
   REAL8 f2aFac = 0;/* factor multiplying f in amplitude function */
   REAL8 apFac = 0, acFac = 0;/* extra factor in plus and cross amplitudes */

   ak   = paramsInit->ak;
   func = paramsInit->func;

   if (params == NULL)
      XLAL_ERROR(XLAL_EFAULT);
   if (params->nStartPad < 0)
      XLAL_ERROR(XLAL_EDOM);
   if (params->nEndPad < 0)
      XLAL_ERROR(XLAL_EDOM);
   if (params->fLower <= 0)
      XLAL_ERROR(XLAL_EDOM);
   if (params->tSampling <= 0)
      XLAL_ERROR(XLAL_EDOM);

   values.length = dvalues.length = valuesNew.length =
   yt.length = dym.length = dyt.length = n;
   dummy.length = n * 6;
   if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * n * 6))) {
      XLAL_ERROR(XLAL_ENOMEM);
   }

   values.data = &dummy.data[0];
   dvalues.data = &dummy.data[n];
   valuesNew.data = &dummy.data[2*n];
   yt.data = &dummy.data[3*n];
   dym.data = &dummy.data[4*n];
   dyt.data = &dummy.data[5*n];

   m = ak.totalmass;
   dt = 1./params->tSampling;

   if (a || h)
   {
      mTot   = params->mass1 + params->mass2;
      etab   = params->mass1 * params->mass2;
      etab  /= mTot;
      etab  /= mTot;
      unitHz = mTot *LAL_MTSUN_SI*(REAL8)LAL_PI;
      cosI   = cos( params->inclination );
      mu     = etab * mTot;
      fFac   = 1.0 / ( 4.0*LAL_TWOPI*LAL_MTSUN_SI*mTot );
      f2aFac = LAL_PI*LAL_MTSUN_SI*mTot*fFac;
      apFac  = acFac = -2.0 * mu * LAL_MRSUN_SI/params->distance;
      apFac *= 1.0 + cosI*cosI;
      acFac *= 2.0*cosI;
      params->nStartPad = 0;
   }

   if (params->totalMass <= 0)
      XLAL_ERROR(XLAL_EDOM);

   t = 0.0;
   in1.t = t;
   in1.t0=ak.t0;
   in1.v0 = ak.v0;
   in1.vlso = ak.vlso;
   in1.totalmass = ak.totalmass;
   in1.dEnergy = func.dEnergy;
   in1.flux = func.flux;
   in1.coeffs = &ak;

   in2.v0 = ak.v0;
   in2.phi0 = params->startPhase;
   in2.dEnergy = func.dEnergy;
   in2.flux = func.flux;
   in2.coeffs = &ak;

   in3.totalmass = ak.totalmass;
   in3.dEnergy = func.dEnergy;
   in3.flux = func.flux;
   in3.coeffs = &ak;
   funcParams = (void *) &in3;

   v = XLALInspiralVelocity(&in1);
   if (XLAL_IS_REAL8_FAIL_NAN(v))
      XLAL_ERROR(XLAL_EFUNC);

   piM = LAL_PI * m;
   f = (v*v*v)/piM;

   fu = params->fCutoff;
   if (fu)
      fHigh = (fu < ak.flso) ? fu : ak.flso;
   else
      fHigh = ak.flso;
   f = (v*v*v)/(LAL_PI*m);

/*
    Check that the highest frequency is less than half
    the sampling frequency - the Nyquist theorem
*/

  if (fHigh >= 0.5/dt)
     XLAL_ERROR(XLAL_EDOM);
  if (fHigh <= params->fLower)
     XLAL_ERROR(XLAL_EDOM);


   p = XLALInspiralPhasing1(v, &in2);
   if (XLAL_IS_REAL8_FAIL_NAN(p))
      XLAL_ERROR(XLAL_EFUNC);

   *(values.data) = v;
   *(values.data+1) = p;

   in4.function = LALInspiralDerivatives;
   in4.x = t;
   in4.y = &values;
   in4.h = dt;
   in4.n = n;
   in4.yt = &yt;
   in4.dym = &dym;
   in4.dyt = &dyt;

   /* Initialize GSL integrator */
   if (!(integrator = XLALRungeKutta4Init(n, &in4)))
     XLAL_ERROR(XLAL_EFUNC);

   count = 0;
   if (signalvec2) {
   params->nStartPad = 0;} /* for template generation, that value must be zero*/

   else if (signalvec1) {
     count = params->nStartPad;
   }

   t = 0.0;
   do {
      /* Free up memory and abort if writing beyond the end of vector*/
      if ((signalvec1 && (UINT4)count >= signalvec1->length) || (ff && (UINT4)count >= ff->length))
      {
          XLALRungeKutta4Free( integrator );
          XLALFree(dummy.data);
          XLAL_ERROR(XLAL_EBADLEN);
      }

      /* Non-injection case */
      if (signalvec1)
      {
         amp = params->signalAmplitude * v*v;
         h1 = amp * cos(p);
         *(signalvec1->data + count) = (REAL4) h1;
	 if (signalvec2)
	 {
            h2 = amp * cos(p+LAL_PI_2);
            *(signalvec2->data + count) = (REAL4) h2;
	 }
      }
      /* Injection case */
      else if (a)
      {
	  int ice, ico;
          ice = 2*count;
          ico = ice + 1;
          omega = v*v*v;

          ff->data[count]       = (REAL4)(omega/unitHz);
          f2a                   = pow (f2aFac * omega, 2./3.);
          a->data[ice]      = (REAL4)(4.*apFac * f2a);
          a->data[ico]    = (REAL4)(4.*acFac * f2a);
          phi->data[count]      = (REAL8)(p);

          if (h)
	  {
	    h->data[ice] = LALInspiralHPlusPolarization( p, v, params );
	    h->data[ico] = LALInspiralHCrossPolarization( p, v, params );
          }

      }

      LALInspiralDerivatives(&values, &dvalues, funcParams);

      in4.dydx = &dvalues;
      in4.x=t;

      if(XLALRungeKutta4(&valuesNew, integrator, funcParams))
        XLAL_ERROR(XLAL_EFUNC);

      *(values.data) = v = *(valuesNew.data);
      *(values.data+1) = p = *(valuesNew.data+1);

      t = (++count-params->nStartPad) * dt;
      f = v*v*v/piM;

   } while (t < ak.tn && f<fHigh);

   params->vFinal = p;
   params->fFinal = f;
   params->tC = t;

   XLALRungeKutta4Free( integrator );
   XLALFree(dummy.data);

   return count;
}
