/*
*  Copyright (C) 2007 Duncan Brown, David McKechan, Thomas Cokelaer
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

/*  <lalVerbatim file="LALInspiralAmplitudeCorrectedWaveCV">
Author: Cokelaer T, McKechan D
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALInspiralAmplitudeCorrectedWave.c} and \texttt{LALInspiralAmplitudeCorrectedWaveTemplates.c}}

The code \texttt{LALInspiralAmplitudeCorrectedWave} generates an time-domain inspiral waveform corresponding to the
\texttt{approximant} \texttt{TaylorT1} and \texttt{PadeT1} as outlined in the
documentation for the function \texttt{LALInspiralWave}.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALInspiralAmplitudeCorrectedWaveCP}
\index{\verb&LALInspiralAmplitudeCorrectedWave()&}
\begin{itemize}
\item {\tt signalvec:} Output containing the inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\input{LALInspiralAmplitudeCorrectedWaveTemplatesCP}
\index{\verb&LALInspiralAmplitudeCorrectedWaveTemplates()&}
\begin{itemize}
\item {\tt signalvec1:} Output containing the 0-phase inspiral waveform.
\item {\tt signalvec2:} Output containing the $\pi/2$-phase inspiral waveform.
\item {\tt params:} Input containing binary chirp parameters.
\end{itemize}

\subsubsection*{Description}

\texttt{LALInspiralAmplitudeCorrectedWave} is called if the user has specified the
\texttt{enum} \texttt{approximant} to be
either \texttt{TaylorT1} or \texttt{PadeT1}.
{\tt LALInspiralAmplitudeCorrectedWaveTemplates} is exactly the same as \texttt{LALInspiralAmplitudeCorrectedWave,} except that
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

\vfill{\footnotesize\input{LALInspiralAmplitudeCorrectedWaveCV}}

</lalLaTeX>  */

/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/GenerateInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/DetResponse.h>
#include <lal/LIGOMetadataUtils.h>
static void
LALInspiralAmplitudeCorrectedWaveEngine(
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   REAL4Vector      *a,
   REAL4Vector      *ff,
   REAL8Vector      *phi,
   INT4             *countback,
   InspiralTemplate *params
   );


NRCSID (LALINSPIRALAMPLITUDECORRECTEDWAVEC, "$Id$");

/*  <lalVerbatim file="LALInspiralAmplitudeCorrectedWaveCP"> */
void
LALInspiralAmplitudeCorrectedWave(
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

   INT4 count;

   INITSTATUS(status, "LALInspiralAmplitudeCorrectedWave",LALINSPIRALAMPLITUDECORRECTEDWAVEC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signalvec, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

   /* Initially the waveform is empty*/
   memset(signalvec->data, 0, signalvec->length*sizeof(REAL4));

   /*Call the engine function*/
   LALInspiralAmplitudeCorrectedWaveEngine(status->statusPtr, signalvec, NULL, NULL, NULL, NULL, &count, params);
   CHECKSTATUSPTR(status);

   DETATCHSTATUSPTR(status);
   RETURN (status);
}



NRCSID (LALINSPIRALAMPLITUDECORRECTEDWAVETEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALInspiralAmplitudeCorrectedWaveTemplatesCP"> */
void
LALInspiralAmplitudeCorrectedWaveTemplates(
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   InspiralTemplate *params
   )
 { /* </lalVerbatim>  */

   INT4 count;

   INITSTATUS(status, "LALInspiralAmplitudeCorrectedWaveTemplates",LALINSPIRALAMPLITUDECORRECTEDWAVETEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signalvec1, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec2, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec1->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec2->data, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

   /* Initially the waveforms are empty */
   memset(signalvec1->data, 0, signalvec1->length * sizeof(REAL4));
   memset(signalvec2->data, 0, signalvec2->length * sizeof(REAL4));

   /* Call the engine function */
   LALInspiralAmplitudeCorrectedWaveEngine(status->statusPtr, signalvec1, signalvec2, NULL, NULL, NULL, &count, params);
   CHECKSTATUSPTR(status);

   DETATCHSTATUSPTR(status);
   RETURN (status);
}

/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms for injection packages T.Cokelaer sept 2003
*/

NRCSID (LALINSPIRALAMPLITUDECORRECTEDWAVEFORINJECTIONC, "$Id$");

/*  <lalVerbatim file="LALInspiralAmplitudeCorrectedWaveForInjectionCP"> */
void
LALInspiralAmplitudeCorrectedWaveForInjection(
			     LALStatus        *status,
			     CoherentGW       *waveform,
			     InspiralTemplate *params,
			     PPNParamStruc  *ppnParams
			     )
{ /* </lalVerbatim>  */

  INT4        count, i;
  REAL8       p, phiC;

  REAL4Vector a;           /* pointers to generated amplitude  data */
  REAL4Vector ff;          /* pointers to generated  frequency data */
  REAL8Vector phi;         /* generated phase data */

  CreateVectorSequenceIn in;

  CHAR message[256];

  InspiralInit paramsInit;

  INITSTATUS(status,"LALInspiralAmplitudeCorrectedWaveForInjection",
  LALINSPIRALAMPLITUDECORRECTEDWAVETEMPLATESC);
  ATTATCHSTATUSPTR(status);

  /* Make sure parameter and waveform structures exist. */
  ASSERT( params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT( !( waveform->a ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->f ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->phi ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );

  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);

  if (paramsInit.nbins == 0){
      DETATCHSTATUSPTR(status);
      RETURN (status);
  }

  /* Now we can allocate memory and vector for coherentGW structure*/

  ff.length  = paramsInit.nbins;
  a.length   = 2* paramsInit.nbins;
  phi.length = paramsInit.nbins;

  ff.data = (REAL4 *) LALCalloc(paramsInit.nbins, sizeof(REAL4));
  a.data  = (REAL4 *) LALCalloc(2 * paramsInit.nbins, sizeof(REAL4));
  phi.data= (REAL8 *) LALCalloc(paramsInit.nbins, sizeof(REAL8));

  /* Check momory allocation is okay */
  if (!(ff.data) || !(a.data) || !(phi.data))
  {
    if (ff.data)  LALFree(ff.data);
    if (a.data)   LALFree(a.data);
    if (phi.data) LALFree(phi.data);


    ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
  }

  count = 0;

  /* Call the engine function */
  LALInspiralAmplitudeCorrectedWaveEngine(status->statusPtr, NULL, NULL, &a, &ff, &phi, &count, params);
  BEGINFAIL( status )
  {
    LALFree(ff.data);
    LALFree(a.data);
    LALFree(phi.data);
  }
  ENDFAIL( status );

  p = phi.data[count-1];

  params->fFinal = ff.data[count-1];
  sprintf(message, "cycles = %f", p/(double)LAL_TWOPI);
  LALInfo(status, message);

  if ( (INT4)(p/LAL_TWOPI) < 2 ){
    sprintf(message, "The waveform has only %f cycles; we don't keep waveform with less than 2 cycles.",
	       p/(double)LAL_TWOPI );
    XLALPrintError(message);
    LALWarning(status, message);
  }

      /*wrap the phase vector*/
      phiC =  phi.data[count-1] ;
      for (i = 0; i < count; i++)
	{
	  phi.data[i] =  phi.data[i] - phiC + ppnParams->phi;
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

      LALSCreateVectorSequence( status->statusPtr, &( waveform->h->data ), &in );
      CHECKSTATUSPTR(status);

      LALSCreateVectorSequence( status->statusPtr, &( waveform->a->data ), &in );
      CHECKSTATUSPTR(status);

      LALSCreateVector( status->statusPtr, &( waveform->f->data ), count);
      CHECKSTATUSPTR(status);

      LALDCreateVector( status->statusPtr, &( waveform->phi->data ), count );
      CHECKSTATUSPTR(status);

      memcpy(waveform->f->data->data , ff.data, count*(sizeof(REAL4)));
      memcpy(waveform->h->data->data , a.data, 2*count*(sizeof(REAL4)));
      memcpy(waveform->a->data->data , a.data, 2*count*(sizeof(REAL4)));
      memcpy(waveform->phi->data->data ,phi.data, count*(sizeof(REAL8)));

      waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT = waveform->h->deltaT
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
				   - waveform->f->data->data[count-2]))
	* ppnParams->deltaT;
      ppnParams->fStop  = params->fFinal;
      ppnParams->termCode        = GENERATEPPNINSPIRALH_EFSTOP;
      ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;

      ppnParams->fStart   = ppnParams->fStartIn;

  /* --- free memory --- */
  LALFree(ff.data);
  LALFree(a.data);
  LALFree(phi.data);

  DETATCHSTATUSPTR(status);
  RETURN (status);
}

/*
 *  Engine function for use by other LALInspiralAmplitudeCorrectedWave* functions
 *  Craig Robinson April 2005
 */

NRCSID (LALINSPIRALAMPLITUDECORRECTEDWAVEENGINEC, "$Id$");

void
LALInspiralAmplitudeCorrectedWaveEngine(
		LALStatus        *status,
		REAL4Vector      *signalvec1,
		REAL4Vector      *signalvec2,
		REAL4Vector      *a,
		REAL4Vector      *ff,
		REAL8Vector      *phi,
		INT4             *countback,
		InspiralTemplate *params)
{
   PPNParamStruc ppnParams;
   CoherentGW 	waveform;
   INT4 i, count;
   REAL8 dt;
   REAL8 mTot = 0;
   REAL8 unitHz = 0;
   REAL8 mu = 0;
   REAL8 cosI = 0;/* cosine of system inclination */
   REAL8 etab = 0;
   REAL8 fFac = 0; /* SI normalization for f and t */
   REAL8 f2aFac = 0;/* factor multiplying f in amplitude function */
   REAL8 apFac = 0, acFac = 0;/* extra factor in plus and cross amplitudes */

   REAL4 hPlus, hCross;
   double fPlus, fCross;

  /* LALDetector det;
   InterferometerNumber ifoNumber = LAL_UNKNOWN_IFO;
   REAL4 longitude,latitude,polarization,gmst;
*/
   INITSTATUS(status, "LALInspiralAmplitudeCorrectedWaveEngine", LALINSPIRALAMPLITUDECORRECTEDWAVEENGINEC);
   ATTATCHSTATUSPTR(status);

   ASSERT (params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);


   /* First, we compute the fplus and fcross*/
/*   memset( &det, 0, sizeof(LALDetector));
   ifoNumber = XLALIFONumber("H1");
   XLALReturnDetector(&det, ifoNumber);
   longitude = 0.4;
   latitude = 0.4;
   gmst = 0.4;
   XLALComputeDetAMResponse(&fPlus, &fCross, det.response, longitude, latitude, params->polarisationAngle, gmst);
  */
  /* For  overhead injection, we just need to compute these 2 expresions for fplus and fcross */
  fPlus = cos(2.*params->polarisationAngle);
  fCross = sin(2.*params->polarisationAngle);

   dt = 1./params->tSampling;

   /* some values to be used to compute h(t)*/
   {
      mTot   = params->mass1 + params->mass2;
      etab   = params->mass1 * params->mass2;
      etab  /= mTot;
      etab  /= mTot;
      unitHz = mTot *LAL_MTSUN_SI*(REAL8)LAL_PI;
      /*cosI and the following parameters are probably useless. apFac and acFac
       * are computed within GeneratePPNAmpCor*/
      cosI   = cos( params->inclination );
      mu     = etab * mTot;
      fFac   = 1.0 / ( 4.0*LAL_TWOPI*LAL_MTSUN_SI*mTot );
      f2aFac = LAL_PI*LAL_MTSUN_SI*mTot*fFac;
      apFac  = acFac = -2.0 * mu * LAL_MRSUN_SI/params->distance;
      apFac *= 1.0 + cosI*cosI;
      acFac *= 2.0*cosI;
      params->nStartPad = 0;
   }

   /* this parameters are not used in GeneratePPNAmp*/
   ppnParams.position.latitude = ppnParams.position.longitude = 0.;
   ppnParams.position.system = COORDINATESYSTEM_EQUATORIAL;
   ppnParams.psi = 0.0;
   ppnParams.lengthIn = 0.0;

   /* The following fields are used by GeneratePPNAmpCor function */
   /* Variable Parameters */
   ppnParams.mTot = mTot;
   ppnParams.eta = etab;
   ppnParams.d = LAL_PC_SI*1.0e3;
   ppnParams.inc = params->inclination;
   /* GeneratePPN does not set the starting phase.*/
   ppnParams.phi = params->startPhase;
   ppnParams.fStartIn = params->fLower;
   /* set to zero so that GeneratePPNAmp recomputes the flso*/
   ppnParams.fStopIn = 0.;
   ppnParams.ppn = NULL;
   ppnParams.ampOrder = params->ampOrder;
   /* set ther PN order of the flux */
   LALSCreateVector( status->statusPtr, &(ppnParams.ppn), params->order + 1 );
   ppnParams.ppn->data[0] = 1.0;
   if ( params->order > 0 )
     ppnParams.ppn->data[1] = 0.0;
   for ( i = 2; i <= (INT4)( params->order ); i++ )
     ppnParams.ppn->data[i] = 1.0;
	ppnParams.fStopIn = 0;
	ppnParams.deltaT = dt;


   count = 0;
   if (signalvec2) {
   params->nStartPad = 0;
   } /* for template genera  memset( &waveform, 0, sizeof(CoherentGW) );
tion, that value must be zero*/
   else if (signalvec1) {
     count = params->nStartPad;
   }

   memset( &waveform, 0, sizeof(CoherentGW) );
   LALGeneratePPNAmpCorInspiral( status->statusPtr, &waveform, &ppnParams );





   count = 0;
  for (i=0;i<(INT4)waveform.h->data->length; i++)
  {
	   /* Non-injection case */
      if (signalvec1)
      {
         /* For amplitude corrected waveforms, we do not want only h+ or only hx but F+h+ + FxHx*/
         hPlus  = (REAL4) waveform.h->data->data[2*i];
         hCross = (REAL4) waveform.h->data->data[2*i+1];
         *(signalvec1->data + count) = fPlus * hPlus + fCross * hCross;
         /* todo: add an Abort if signalvec2<>0*/
	 if (signalvec2)
	 {
         *(signalvec2->data + count) = (REAL4) waveform.h->data->data[2*i+1];
	 }
      }

      /* Injection case */
      else if (a)
      {
          #if 0
          omega = v*v*v;
          ff->data[count]       = (REAL4)(omega/unitHz);
          f2a                   = pow (f2aFac * omega, 2./3.);
          /* waveform->a does not exist. Only waveform->h is populated within
          GeneratePPNAmpCorr function.
          */
          a->data[2*count]      = (REAL4) waveform.a->data->data[count];
          a->data[2*count+1]    = (REAL4) waveform.a->data->data[count];
          phi->data[count]      = (REAL8) waveform.phi->data->data[count];
          #endif
      }

 		count++;
   }


   params->tC = count*dt;

   /* The final frequency needs to take into account the amplitude corrected order */
   params->fFinal = waveform.f->data->data[count-1] * (params->ampOrder+2.)/2.;


   *countback = count;

	/* destroy the waveform signal */
   if ( waveform.a ){
     LALSDestroyVectorSequence( status->statusPtr, &(waveform.a->data) );
     CHECKSTATUSPTR( status );
     LALFree( waveform.a );
   }
   if ( waveform.shift )
    {
      LALSDestroyVector( status->statusPtr, &(waveform.shift->data) );
      CHECKSTATUSPTR( status );
      LALFree( waveform.shift );
    }

    if ( waveform.f ){
    LALSDestroyVector( status->statusPtr, &(waveform.f->data) );
    CHECKSTATUSPTR( status );
      LALFree( waveform.f );
    }

    if ( waveform.phi ){
    LALDDestroyVector( status->statusPtr, &(waveform.phi->data) );
    CHECKSTATUSPTR( status );
      LALFree( waveform.phi );
    }

   if ( waveform.h ){
    LALSDestroyVectorSequence( status->statusPtr, &(waveform.h->data) );
    CHECKSTATUSPTR( status );
    LALFree( waveform.h );
    }


   DETATCHSTATUSPTR(status);
   RETURN (status);

}
