/*
*  Copyright (C) 2007 Stas Babak, David Churches, Jolien Creighton, B.S. Sathyaprakash, Craig Robinson , Thomas Cokelaer, Drew Keppel
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
 * \brief These modules generate a time-domain chirp waveform of type #TaylorT2.
 *
 * \heading{Prototypes}
 *
 * <tt>LALInspiralWave2()</tt>
 * <ul>
 * <li> \c output: Output containing the inspiral waveform.</li>
 * <li> \c params: Input containing binary chirp parameters.</li>
 * </ul>
 *
 * <tt>LALInspiralWave2Templates()</tt>
 * <ul>
 * <li> \c output1: Output containing the 0-phase inspiral waveform.</li>
 * <li> \c output2: Output containing the \f$\pi/2\f$-phase inspiral waveform.</li>
 * <li> \c params: Input containing binary chirp parameters.</li>
 * </ul>
 *
 * \heading{Description}
 *
 * LALInspiralWave2() generates #TaylorT2 approximant wherein
 * the phase of the waveform is given as an implicit function of time
 * as in Equation.\eqref{eq_InspiralWavePhase2}. A template is required
 * to be sampled at equal intervals of time. Thus, first of the equations
 * in Equation.\eqref{eq_InspiralWavePhase2} is solved for \f$v\f$ at equally
 * spaced values of the time steps
 * \f$t_k\f$ and the resulting value of \f$v_k\f$ is used in the second equation to
 * obtain the phase \f$\phi_k\f$.
 *
 * LALInspiralWave2Templates() is exactly the same as LALInspiralWave2()
 * except that it generates two waveforms that differ in phase by \f$\pi/2.\f$
 *
 * \heading{Uses}
 * \code
 * LALInspiralParameterCalc()
 * LALDBisectionFindRoot()
 * LALInspiralPhasing2()
 * \endcode
 *
 */

#include <lal/LALStdlib.h>
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>

static int
XLALInspiralWave2Engine(
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
LALInspiralWave2 (
   LALStatus        *status,
   REAL4Vector      *output,
   InspiralTemplate *params
   )
{
  XLALPrintDeprecationWarning("LALInspiralWave2", "XLALInspiralWave2");
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if( XLALInspiralWave2(output, params) )
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

int
XLALInspiralWave2(
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
  count = XLALInspiralWave2Engine(output, NULL, NULL, NULL,
			NULL, NULL, params, &paramsInit);
  if (count < 0)
    XLAL_ERROR(XLAL_EFUNC);

  return XLAL_SUCCESS;
}

void
LALInspiralWave2Templates (
  LALStatus        *status,
  REAL4Vector      *output1,
  REAL4Vector      *output2,
  InspiralTemplate *params
  )
{
  XLALPrintDeprecationWarning("LALInspiralWave2Templates", "XLALInspiralWave2Templates");
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if( XLALInspiralWave2Templates(output1, output2, params) )
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

int
XLALInspiralWave2Templates(
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
  count = XLALInspiralWave2Engine(output1, output2, NULL, NULL,
			   NULL, NULL, params, &paramsInit);

  if (count < 0)
    XLAL_ERROR(XLAL_EFUNC);

  return XLAL_SUCCESS;
}

void
LALInspiralWave2ForInjection (
  LALStatus        *status,
  CoherentGW       *waveform,
  InspiralTemplate *params,
  PPNParamStruc  *ppnParams
  )
{
  XLALPrintDeprecationWarning("LALInspiralWave2ForInjection", "XLALInspiralWave2ForInjection");
  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if( XLALInspiralWave2ForInjection(waveform, params, ppnParams) )
    ABORTXLAL(status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

int
XLALInspiralWave2ForInjection(
  CoherentGW *waveform,
  InspiralTemplate *params,
  PPNParamStruc  *ppnParams
  )
{
  INT4 count;
  UINT4 i;
  REAL4Vector *a   = NULL;/* pointers to generated amplitude  data */
  REAL4Vector *h   = NULL;/* pointers to generated polarizations */
  REAL4Vector *ff  = NULL;/* pointers to generated  frequency data */
  REAL8Vector *phi = NULL;/* pointer to generated phase data */

  REAL8 phiC;/* phase at coalescence */
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
  count = XLALInspiralWave2Engine(NULL, NULL, h, a, ff,
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

  if ( fabs(phi->data[count-1]/2.)/LAL_PI < 2. ){
        sprintf(message, "The waveform has only %f cycles; we don't keep waveform with less than 2 cycles.",
	       (double)(fabs(phi->data[count-1]/2.)/LAL_PI) );
    XLALPrintWarning(message);
  }
  else
    {
    /*wrap the phase vector*/
    phiC =  phi->data[count-1] ;
    for (i=0; i<(UINT4)count;i++)
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

    snprintf( waveform->a->name, LALNameLength,   "T2 inspiral amplitude" );
    snprintf( waveform->f->name, LALNameLength,   "T2 inspiral frequency" );
    snprintf( waveform->phi->name, LALNameLength, "T2 inspiral phase" );


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

/* 'Engine' function upon which all the other functions invoke
    Craig Robinson 04/05 */



static int
XLALInspiralWave2Engine(
                REAL4Vector      *output1,
                REAL4Vector      *output2,
                REAL4Vector      *h,
                REAL4Vector      *a,
                REAL4Vector      *ff,
                REAL8Vector      *phi,
                InspiralTemplate *params,
		InspiralInit     *paramsInit
                )
{

  REAL8 amp, dt, fs, fu, fHigh, phase0, phase1, tC;
  REAL8 phase, v, totalMass, fLso, freq, fOld;
  INT4 i, startShift, count;
  REAL8 xmin,xmax,xacc;
  REAL8 (*timing2)(REAL8, void *);
  InspiralToffInput toffIn;
  void *funcParams;
  expnCoeffs ak;
  expnFunc func;

  /* Variables only used in injection case.. */
  REAL8 omega;
  REAL8 unitHz = 0;
  REAL8 f2a = 0;
  REAL8 mu = 0;
  REAL8 mTot = 0;
  REAL8 cosI = 0;/* cosine of system inclination */
  REAL8 etab =0;
  REAL8 fFac = 0; /* SI normalization for f and t */
  REAL8 f2aFac = 0;/* factor multiplying f in amplitude function */
  REAL8 apFac = 0, acFac = 0;/* extra factor in plus and cross amplitudes */


  ak   = paramsInit->ak;
  func = paramsInit->func;

  if (output2 || a)
         params->nStartPad = 0;   /* that value must be zero for template generation */
  dt = 1.0/(params->tSampling);   /* sampling interval */
  fs = params->fLower;            /* lower frequency cutoff */
  fu = params->fCutoff;           /* upper frequency cutoff */
  startShift = params->nStartPad; /* number of bins to pad at the beginning */
  phase0 = params->startPhase;    /* initial phasea */
  phase1 = phase0 + LAL_PI_2;

  timing2 = func.timing2; /* function to solve for v, given t:*/

  if (a || h)           /* Only used in injection case */
  {
    mTot   =  params->mass1 + params->mass2;
    etab   =  params->mass1 * params->mass2;
    etab  /= mTot;
    etab  /= mTot;
    unitHz = (mTot) *LAL_MTSUN_SI*(REAL8)LAL_PI;
    cosI   = cos( params->inclination );
    mu     = etab * mTot;
    fFac   = 1.0 / ( 4.0*LAL_TWOPI*LAL_MTSUN_SI*mTot );
    f2aFac = LAL_PI*LAL_MTSUN_SI*mTot*fFac;
    apFac  = acFac = -2.0 * mu * LAL_MRSUN_SI/params->distance;
    apFac *= 1.0 + cosI*cosI;
    acFac *= 2.0 * cosI;
  }

  if (params->totalMass <= 0)
    XLAL_ERROR(XLAL_EDOM);
  if (params->eta < 0 || params->eta > 0.25)
    XLAL_ERROR(XLAL_EDOM);

  totalMass = params->totalMass*LAL_MTSUN_SI; /* convert from solar masses to seconds */

  toffIn.tN = ak.tvaN;
  toffIn.t2 = ak.tva2;
  toffIn.t3 = ak.tva3;
  toffIn.t4 = ak.tva4;
  toffIn.t5 = ak.tva5;
  toffIn.t6 = ak.tva6;
  toffIn.t7 = ak.tva7;
  toffIn.tl6 = ak.tvl6;
  toffIn.piM = ak.totalmass * LAL_PI;

  /* Determine the total chirp-time tC: the total chirp time is
     timing2(v0;tC,t) with t=tc=0*/

  toffIn.t = 0.;
  toffIn.tc = 0.;
  funcParams = (void *) &toffIn;
  tC = func.timing2(fs, funcParams);
  if (XLAL_IS_REAL8_FAIL_NAN(tC))
    XLAL_ERROR(XLAL_EFUNC);
  /* Reset chirp time in toffIn structure */
  toffIn.tc = -tC;

  /* Determine the initial phase: it is phasing2(v0) with ak.phiC=0 */
  v = pow(fs * LAL_PI * totalMass, (1./3.));
  ak.phiC = 0.0;
  phase = func.phasing2(v, &ak);
  if (XLAL_IS_REAL8_FAIL_NAN(phase))
    XLAL_ERROR(XLAL_EFUNC);
  ak.phiC = -phase;

  /*
     If flso is less than the user inputted upper frequency cutoff fu,
  */

  fLso = ak.fn;
  if (fu)
    fHigh = (fu < fLso) ? fu : fLso;
  else
    fHigh = fLso;

  /* Is the sampling rate large enough? */

  if (fHigh >= 0.5/dt)
    XLAL_ERROR(XLAL_EDOM);
  if (fHigh <= params->fLower)
    XLAL_ERROR(XLAL_EDOM);

  xmax = 1.2*fu;
  xacc = 1.0e-8;
  xmin = 0.999999*fs;

  i = startShift;

  /* Now cast the input structure to argument 4 of BisectionFindRoot so that it
     of type void * rather than InspiralToffInput  */

  funcParams = (void *) &toffIn;

  toffIn.t = 0.0;
  freq = fs;
  count=0;
  do
    {
    /*
    Check we're not writing past the end of the vector
    */
    if ((output1 && ((UINT4)i >= output1->length)) || (ff && ((UINT4)count >= ff->length)))
    {
      XLAL_ERROR(XLAL_EBADLEN);
    }

    fOld = freq;
    v = pow(freq*toffIn.piM, (1./3.));
    phase = func.phasing2(v, &ak); /* phase at given v */
    if (XLAL_IS_REAL8_FAIL_NAN(phase))
      XLAL_ERROR(XLAL_EFUNC);

    amp = params->signalAmplitude*v*v;

    if (output1)
    {
      output1->data[i]=(REAL4)(amp*cos(phase+phase0));
      if (output2)
        output2->data[i]=(REAL4)(amp*cos(phase+phase1));
    }
    else
    {
      int ice, ico;
      ice = 2*count;
      ico = ice + 1;
      omega = v*v*v;

      ff->data[count]= (REAL4)(omega/unitHz);
      f2a = pow (f2aFac * omega, 2./3.);
      a->data[ice]          = (REAL4)(4.*apFac * f2a);
      a->data[ico]        = (REAL4)(4.*acFac * f2a);
      phi->data[count]          = (REAL8)(phase);

      if(h)
      {
        h->data[ice] = LALInspiralHPlusPolarization( phase, v, params );
        h->data[ico] = LALInspiralHCrossPolarization( phase, v, params );
      }
    }
    i++;
    ++count;
    toffIn.t=count*dt;
    /*
       Determine the frequency at the current time by solving timing2(v;tC,t)=0
    */
    xmin = 0.8*freq;
    freq = XLALDBisectionFindRoot(timing2, xmin, xmax, xacc, funcParams);
    if (XLAL_IS_REAL8_FAIL_NAN(freq))
      XLAL_ERROR(XLAL_EFUNC);
    } while (freq < fHigh && freq > fOld && toffIn.t < -tC);
  params->fFinal = fOld;
  if (output1 && !(output2))   params->tC = toffIn.t;

  return count;
}
