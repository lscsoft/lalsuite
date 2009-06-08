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

/*  <lalVerbatim file="LALInspiralWave1CV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralWave1.c} and \texttt{LALInspiralWave1Templates.c}}

The code \texttt{LALInspiralWave1} generates an time-domain inspiral waveform corresponding to the
\texttt{approximant} \texttt{TaylorT1} and \texttt{PadeT1} as outlined in the
documentation for the function \texttt{LALInspiralWave}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralWave1CP}
\index{\verb&LALInspiralWave1()&}
\begin{itemize}
\item {\tt signalvec:} Output containing the inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\input{LALInspiralWave1TemplatesCP}
\index{\verb&LALInspiralWave1Templates()&}
\begin{itemize}
\item {\tt signalvec1:} Output containing the 0-phase inspiral waveform.
\item {\tt signalvec2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\subsubsection*{Description}

\texttt{LALInspiralWave1} is called if the user has specified the
\texttt{enum} \texttt{approximant} to be
either \texttt{TaylorT1} or \texttt{PadeT1}.
{\tt LALInspiralWave1Templates} is exactly the same as \texttt{LALInspiralWave1,} except that
it generates two templates one for which the starting phase is
\texttt{params.startPhase} and the other for which the phase is
\texttt{params.startPhase + $\pi/2$}.


\subsubsection*{Algorithm}
This code uses a fourth-order Runge-Kutta algorithm to solve the ODEs
in Equation (\ref{eq:ode2}).

\subsubsection*{Uses}

\texttt{LALInspiralSetup}\\
\texttt{LALInspiralChooseModel}\\
\texttt{LALInspiralVelocity}\\
\texttt{LALInspiralPhasing1}\\
\texttt{LALInspiralDerivatives}\\
\texttt{LALRungeKutta4}.


\subsubsection*{Notes}

\vfill{\footnotesize\input{LALInspiralWave1CV}}

</lalLaTeX>  */

/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>

static void
LALInspiralWave1Engine(
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   REAL4Vector      *h,
   REAL4Vector      *a,
   REAL4Vector      *ff,
   REAL8Vector      *phi,
   INT4             *countback,
   InspiralTemplate *params
   );


NRCSID (LALINSPIRALWAVE1C, "$Id$");

/*  <lalVerbatim file="LALInspiralWave1CP"> */
void
LALInspiralWave1(
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

   INT4 count;

   INITSTATUS(status, "LALInspiralWave1", LALINSPIRALWAVE1C);
   ATTATCHSTATUSPTR(status);

   ASSERT(signalvec, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


   /* Initially the waveform is empty*/
   memset(signalvec->data, 0, signalvec->length*sizeof(REAL4));

   /*Call the engine function*/
   LALInspiralWave1Engine(status->statusPtr, signalvec, NULL, NULL, NULL, NULL, NULL, &count, params);
   CHECKSTATUSPTR(status);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}



/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/

NRCSID (LALINSPIRALWAVE1TEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave1TemplatesCP"> */
void
LALInspiralWave1Templates(
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

   INT4 count;

   INITSTATUS(status, "LALInspiralWave1Templates", LALINSPIRALWAVE1TEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signalvec1, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec2, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec1->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec2->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

   /* Initially the waveforms are empty */
   memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));
   memset(signalvec2->data, 0, signalvec2->length * sizeof(REAL4));

   /* Call the engine function */
   LALInspiralWave1Engine(status->statusPtr, signalvec1, signalvec2, NULL, NULL, NULL, NULL, &count, params);
   CHECKSTATUSPTR(status);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}

/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms for injection packages T.Cokelaer sept 2003
*/

NRCSID (LALINSPIRALWAVE1FORINJECTIONC, "$Id$");

/*  <lalVerbatim file="LALInspiralWave1ForInjectionCP"> */
void
LALInspiralWave1ForInjection(
			     LALStatus        *status,
			     CoherentGW       *waveform,
			     InspiralTemplate *params,
			     PPNParamStruc  *ppnParams
			     )
{ /* </lalVerbatim>  */

  INT4        count, i;
  REAL8       p, phiC;

  REAL4Vector *a   = NULL;      /* pointers to generated amplitude  data */
  REAL4Vector *h   = NULL;      /* pointers to generated polarizations */
  REAL4Vector *ff  = NULL;      /* pointers to generated  frequency data */
  REAL8Vector *phi = NULL;      /* pointer to generated phase data */


  CreateVectorSequenceIn in;

  CHAR message[256];

  InspiralInit paramsInit;


  INITSTATUS(status, "LALInspiralWave1ForInjection", LALINSPIRALWAVE1TEMPLATESC);
  ATTATCHSTATUSPTR(status);

  /* Make sure parameter and waveform structures exist. */
  ASSERT( params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT( !( waveform->a ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->h ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->f ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->phi ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );

  params->ampOrder = 0;
  sprintf(message, "WARNING: Amp Order has been reset to %d", params->ampOrder);
  LALInfo(status, message);
  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);

  if (paramsInit.nbins == 0){
      DETATCHSTATUSPTR(status);
      RETURN (status);
  }

  /* Now we can allocate memory and vector for coherentGW structure*/
  LALSCreateVector(status->statusPtr, &ff, paramsInit.nbins);
  CHECKSTATUSPTR(status);
  LALSCreateVector(status->statusPtr, &a, 2*paramsInit.nbins);
  CHECKSTATUSPTR(status);
  LALDCreateVector(status->statusPtr, &phi, paramsInit.nbins);
  CHECKSTATUSPTR(status);

  /* By default the waveform is empty */
  memset(ff->data, 0, paramsInit.nbins * sizeof(REAL4));
  memset(a->data, 0, 2 * paramsInit.nbins * sizeof(REAL4));
  memset(phi->data, 0, paramsInit.nbins * sizeof(REAL8));


  if( params->ampOrder )
  {
    LALSCreateVector(status->statusPtr, &h, 2*paramsInit.nbins);
    CHECKSTATUSPTR(status);
    memset(h->data, 0, 2 * paramsInit.nbins * sizeof(REAL4));
  }

  count = 0;

  /* Call the engine function */
  LALInspiralWave1Engine(status->statusPtr, NULL, NULL, h, a, ff,
		             phi, &count, params);
  BEGINFAIL( status )
  {
     LALSDestroyVector(status->statusPtr, &ff);
     CHECKSTATUSPTR(status);
     LALSDestroyVector(status->statusPtr, &a);
     CHECKSTATUSPTR(status);
     LALDDestroyVector(status->statusPtr, &phi);
     CHECKSTATUSPTR(status);
     if( params->ampOrder )
     {
       LALSDestroyVector(status->statusPtr, &h);
       CHECKSTATUSPTR(status);
     }
  }
  ENDFAIL( status );

  p = phi->data[count-1];

  params->fFinal = ff->data[count-1];
  sprintf(message, "cycles = %f", fabs(p)/(double)LAL_TWOPI);
  LALInfo(status, message);

  if ( (INT4)(p/LAL_TWOPI) < 2 ){
    sprintf(message, "The waveform has only %f cycles; we don't keep waveform with less than 2 cycles.",
	       fabs(p)/(double)LAL_TWOPI );
    XLALPrintError(message);
    LALWarning(status, message);
  }
  else
  {

      /*wrap the phase vector*/
      phiC =  phi->data[count-1] ;
      for (i = 0; i < count; i++)
	{
	  phi->data[i] =  phi->data[i] - phiC + ppnParams->phi;
	}

      /* Allocate the waveform structures. */
      if ( ( waveform->a = (REAL4TimeVectorSeries *)
	     LALCalloc(1, sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }
      if ( ( waveform->f = (REAL4TimeSeries *)
	     LALCalloc(1, sizeof(REAL4TimeSeries) ) ) == NULL ) {
	LALFree( waveform->a ); waveform->a = NULL;
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }
      if ( ( waveform->phi = (REAL8TimeSeries *)
	     LALCalloc(1, sizeof(REAL8TimeSeries) ) ) == NULL ) {
	LALFree( waveform->a ); waveform->a = NULL;
	LALFree( waveform->f ); waveform->f = NULL;
	ABORT( status, LALINSPIRALH_EMEM,
	       LALINSPIRALH_MSGEMEM );
      }

      in.length = (UINT4)(count);
      in.vectorLength = 2;

      LALSCreateVectorSequence( status->statusPtr, &( waveform->a->data ), &in );
      CHECKSTATUSPTR(status);

      LALSCreateVector( status->statusPtr, &( waveform->f->data ), count);
      CHECKSTATUSPTR(status);

      LALDCreateVector( status->statusPtr, &( waveform->phi->data ), count );
      CHECKSTATUSPTR(status);


      memcpy(waveform->f->data->data , ff->data, count*(sizeof(REAL4)));
      memcpy(waveform->a->data->data , a->data, 2*count*(sizeof(REAL4)));
      memcpy(waveform->phi->data->data ,phi->data, count*(sizeof(REAL8)));

      waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT
	= ppnParams->deltaT;

      waveform->a->sampleUnits    = lalStrainUnit;
      waveform->f->sampleUnits    = lalHertzUnit;
      waveform->phi->sampleUnits  = lalDimensionlessUnit;
      waveform->position = ppnParams->position;
      waveform->psi = ppnParams->psi;

      LALSnprintf( waveform->a->name, LALNameLength,   "T1 inspiral amplitude" );
      LALSnprintf( waveform->f->name, LALNameLength,   "T1 inspiral frequency" );
      LALSnprintf( waveform->phi->name, LALNameLength, "T1 inspiral phase" );

      /* --- fill some output ---*/
      ppnParams->tc     = (double)(count-1) / params->tSampling ;
      ppnParams->length = count;
      ppnParams->dfdt   = ((REAL4)(waveform->f->data->data[count-1]
			- waveform->f->data->data[count-2])) * ppnParams->deltaT;
      ppnParams->fStop  = params->fFinal;
      ppnParams->termCode        = GENERATEPPNINSPIRALH_EFSTOP;
      ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;

      ppnParams->fStart   = ppnParams->fStartIn;

      if ( params->ampOrder )
      {
        if ( ( waveform->h = (REAL4TimeVectorSeries *)
	       LALCalloc(1, sizeof(REAL4TimeVectorSeries) ) ) == NULL )
        {
	   ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
        }
        LALSCreateVectorSequence( status->statusPtr, &( waveform->h->data ), &in );
        CHECKSTATUSPTR(status);
        memcpy(waveform->h->data->data , h->data, 2*count*(sizeof(REAL4)));
        waveform->h->deltaT = ppnParams->deltaT;
        waveform->h->sampleUnits    = lalStrainUnit;
        LALSnprintf( waveform->h->name, LALNameLength,   "T1 inspiral polarizations" );
        LALSDestroyVector(status->statusPtr, &h);
        CHECKSTATUSPTR(status);
      }
  }

  /* --- free memory --- */
  LALSDestroyVector(status->statusPtr, &ff);
  CHECKSTATUSPTR(status);
  LALSDestroyVector(status->statusPtr, &a);
  CHECKSTATUSPTR(status);
  LALDDestroyVector(status->statusPtr, &phi);
  CHECKSTATUSPTR(status);

  DETATCHSTATUSPTR(status);
  RETURN (status);
}

/*
 *  Engine function for use by other LALInspiralWave1* functions
 *  Craig Robinson April 2005
 */

NRCSID (LALINSPIRALWAVE1ENGINEC, "$Id$");

void
LALInspiralWave1Engine(
		LALStatus        *status,
		REAL4Vector      *signalvec1,
		REAL4Vector      *signalvec2,
		REAL4Vector      *h,
		REAL4Vector      *a,
		REAL4Vector      *ff,
		REAL8Vector      *phi,
		INT4             *countback,
		InspiralTemplate *params)
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


   INITSTATUS(status, "LALInspiralWave1Engine", LALINSPIRALWAVE1ENGINEC);
   ATTATCHSTATUSPTR(status);

   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   values.length = dvalues.length = valuesNew.length =
   yt.length = dym.length = dyt.length = n;
   dummy.length = n * 6;
   if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * n * 6))) {
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
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

   ASSERT(ak.totalmass > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

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

   LALInspiralVelocity(status->statusPtr, &v, &in1);
   CHECKSTATUSPTR(status);

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

   ASSERT(fHigh < 0.5/dt, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(fHigh > params->fLower, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


   LALInspiralPhasing1(status->statusPtr, &p, v, &in2);
   CHECKSTATUSPTR(status);


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

   xlalErrno = 0;
   /* Initialize GSL integrator */
   if (!(integrator = XLALRungeKutta4Init(n, &in4)))
   {
     INT4 errNum = XLALClearErrno();
     LALFree(dummy.data);

     if (errNum == XLAL_ENOMEM)
       ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
     else
       ABORTXLAL( status );
   }

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
          LALFree(dummy.data);
          ABORT(status, LALINSPIRALH_EVECTOR, LALINSPIRALH_MSGEVECTOR);
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
      CHECKSTATUSPTR(status);

      in4.dydx = &dvalues;
      in4.x=t;

      LALRungeKutta4(status->statusPtr, &valuesNew, integrator, funcParams);
      CHECKSTATUSPTR(status);

      *(values.data) = v = *(valuesNew.data);
      *(values.data+1) = p = *(valuesNew.data+1);

      t = (++count-params->nStartPad) * dt;
      f = v*v*v/piM;

   } while (t < ak.tn && f<fHigh);

   params->vFinal = p;
   params->fFinal = f;
   params->tC = t;

   *countback = count;

   XLALRungeKutta4Free( integrator );
   LALFree(dummy.data);

   DETATCHSTATUSPTR(status);
   RETURN (status);

}
