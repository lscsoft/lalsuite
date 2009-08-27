/*
*  Copyright (C) 2008 B.S. Sathyaprakash, Bala R. Iyer
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

/*  <lalVerbatim file="LALTaylorT4WaveformCV">
Author: Sathyaprakash, B. S., Cokelaer T.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTaylorT4Waveform.c} and
\texttt{LALTaylorT4WaveformTemplates.c}}

\vfill{\footnotesize\input{LALTaylorT4WaveformCV}}

</lalLaTeX>  */
/*-------------------------------------------------------------------*/
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/FindRoot.h>
#include <lal/Units.h>

void LALTaylorT4Derivatives4PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);

void LALTaylorT4Derivatives5PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);

void LALTaylorT4Derivatives6PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);

void LALTaylorT4Derivatives7PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);

void LALTaylorT4Waveform (
  LALStatus        *status,
  REAL4Vector      *signalvec,
  InspiralTemplate *params
);

void LALTaylorT4WaveformEngine (
  LALStatus        *status,
  REAL4Vector      *signalvec,
  REAL4Vector      *a,
  REAL4Vector      *ff,
  REAL8Vector      *phi,
  UINT4            *countback,
  InspiralTemplate *params,
  InspiralInit     *paramsInit
);

NRCSID (LALTAYLORT4WAVEFORMC,
"$Id$");

void LALTaylorT4Derivatives4PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   InspiralDerivativesIn *ak;
   REAL8 v, v2, v3, v4, v5, v6, v7, v9, eta, eta2, eta3, pisq, fourpi;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   v = values->data[0];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
   v2 = v*v;
   v3 = v2*v;
   v4 = v3*v;
   v5 = v4*v;
   v6 = v5*v;
   v7 = v6*v;
   v9 = v6*v3;
   pisq = LAL_PI*LAL_PI;
   fourpi = 4.*LAL_PI;

   /*--------------------------------------------------------------------*/
   /* The Expression for dv/dt in the code is simplied after             */
   /* evaluating all the fractions so as to speed of the code            */
   /* The PN expansion from which it results is given by: dv/dt =        */
   /* (32*v^9*\[Eta]*(1 + 4*Pi*v^3 + v^2*(-743/336 - (11*\[Eta])/4)      */
   /* + v^5*((-4159*Pi)/672 - (189*Pi*\[Eta])/8)                         */
   /* + v^4*(34103/18144 + (13661*\[Eta])/2016 + (59*\[Eta]^2)/18)       */
   /* + v^7*((-4415*Pi)/4032 + (358675*Pi*\[Eta])/6048                   */
   /* + (91495*Pi*\[Eta]^2)/1512) + v^6*(16447322263/139708800           */
   /* - (1712*EulerGamma)/105 + (16*Pi^2)/3 - (56198689*\[Eta])/217728   */
   /* + (451*Pi^2*\[Eta])/48 + (541*\[Eta]^2)/896 - (5605*\[Eta]^3)/2592 */
   /* - (856*Log[16])/105 - (1712*Log[v])/105)))/5                       */
   /*--------------------------------------------------------------------*/

   dvalues->data[0] = 6.4*v9*eta/ak->totalmass*(1.
     + v2*(-2.2113095238095237 - 2.75*eta) + v3*fourpi
     + v4*(1.8795745149911816 + 6.77628968*eta  + 3.27777778*eta2));

   dvalues->data[1] =  v3/ak->totalmass;
}

void LALTaylorT4Derivatives5PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   InspiralDerivativesIn *ak;
   REAL8 v, v2, v3, v4, v5, v6, v7, v9, eta, eta2, eta3, pisq, fourpi;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   v = values->data[0];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
   v2 = v*v;
   v3 = v2*v;
   v4 = v3*v;
   v5 = v4*v;
   v6 = v5*v;
   v7 = v6*v;
   v9 = v6*v3;
   pisq = LAL_PI*LAL_PI;
   fourpi = 4.*LAL_PI;

   /*--------------------------------------------------------------------*/
   /* The Expression for dv/dt in the code is simplied after             */
   /* evaluating all the fractions so as to speed of the code            */
   /* The PN expansion from which it results is given by: dv/dt =        */
   /* (32*v^9*\[Eta]*(1 + 4*Pi*v^3 + v^2*(-743/336 - (11*\[Eta])/4)      */
   /* + v^5*((-4159*Pi)/672 - (189*Pi*\[Eta])/8)                         */
   /* + v^4*(34103/18144 + (13661*\[Eta])/2016 + (59*\[Eta]^2)/18)       */
   /* + v^7*((-4415*Pi)/4032 + (358675*Pi*\[Eta])/6048                   */
   /* + (91495*Pi*\[Eta]^2)/1512) + v^6*(16447322263/139708800           */
   /* - (1712*EulerGamma)/105 + (16*Pi^2)/3 - (56198689*\[Eta])/217728   */
   /* + (451*Pi^2*\[Eta])/48 + (541*\[Eta]^2)/896 - (5605*\[Eta]^3)/2592 */
   /* - (856*Log[16])/105 - (1712*Log[v])/105)))/5                       */
   /*--------------------------------------------------------------------*/

   dvalues->data[0] = 6.4*v9*eta/ak->totalmass*(1.
     + v2*(-2.2113095238095237 - 2.75*eta) + v3*fourpi
     + v4*(1.8795745149911816 + 6.77628968*eta  + 3.27777778*eta2)
     + v5*(-6.1889881*LAL_PI - 23.625*LAL_PI*eta));

   dvalues->data[1] =  v3/ak->totalmass;
}

void LALTaylorT4Derivatives6PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   InspiralDerivativesIn *ak;
   REAL8 v, v2, v3, v4, v5, v6, v7, v9, eta, eta2, eta3, pisq, fourpi;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   v = values->data[0];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
   v2 = v*v;
   v3 = v2*v;
   v4 = v3*v;
   v5 = v4*v;
   v6 = v5*v;
   v7 = v6*v;
   v9 = v6*v3;
   pisq = LAL_PI*LAL_PI;
   fourpi = 4.*LAL_PI;

   /*--------------------------------------------------------------------*/
   /* The Expression for dv/dt in the code is simplied after             */
   /* evaluating all the fractions so as to speed of the code            */
   /* The PN expansion from which it results is given by: dv/dt =        */
   /* (32*v^9*\[Eta]*(1 + 4*Pi*v^3 + v^2*(-743/336 - (11*\[Eta])/4)      */
   /* + v^5*((-4159*Pi)/672 - (189*Pi*\[Eta])/8)                         */
   /* + v^4*(34103/18144 + (13661*\[Eta])/2016 + (59*\[Eta]^2)/18)       */
   /* + v^7*((-4415*Pi)/4032 + (358675*Pi*\[Eta])/6048                   */
   /* + (91495*Pi*\[Eta]^2)/1512) + v^6*(16447322263/139708800           */
   /* - (1712*EulerGamma)/105 + (16*Pi^2)/3 - (56198689*\[Eta])/217728   */
   /* + (451*Pi^2*\[Eta])/48 + (541*\[Eta]^2)/896 - (5605*\[Eta]^3)/2592 */
   /* - (856*Log[16])/105 - (1712*Log[v])/105)))/5                       */
   /*--------------------------------------------------------------------*/

   dvalues->data[0] = 6.4*v9*eta/ak->totalmass*(1.
     + v2*(-2.2113095238095237 - 2.75*eta) + v3*fourpi
     + v4*(1.8795745149911816 + 6.77628968*eta  + 3.27777778*eta2)
     + v5*(-6.1889881*LAL_PI - 23.625*LAL_PI*eta)
     + v6*(117.72574285227559 - 16.3047619*0.577215664901
	     + 5.3333333*pisq - 258.114202*eta
	     + 9.39583333*pisq*eta + 0.603794643*eta2
	     - 2.16242284*eta3 - 8.15238095*log(16*v2)));

   dvalues->data[1] =  v3/ak->totalmass;
}

void LALTaylorT4Derivatives7PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   InspiralDerivativesIn *ak;
   REAL8 v, v2, v3, v4, v5, v6, v7, v9, eta, eta2, eta3, pisq, fourpi;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   v = values->data[0];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
   v2 = v*v;
   v3 = v2*v;
   v4 = v3*v;
   v5 = v4*v;
   v6 = v5*v;
   v7 = v6*v;
   v9 = v6*v3;
   pisq = LAL_PI*LAL_PI;
   fourpi = 4.*LAL_PI;

   /*--------------------------------------------------------------------*/
   /* The Expression for dv/dt in the code is simplied after             */
   /* evaluating all the fractions so as to speed of the code            */
   /* The PN expansion from which it results is given by: dv/dt =        */
   /* (32*v^9*\[Eta]*(1 + 4*Pi*v^3 + v^2*(-743/336 - (11*\[Eta])/4)      */
   /* + v^5*((-4159*Pi)/672 - (189*Pi*\[Eta])/8)                         */
   /* + v^4*(34103/18144 + (13661*\[Eta])/2016 + (59*\[Eta]^2)/18)       */
   /* + v^7*((-4415*Pi)/4032 + (358675*Pi*\[Eta])/6048                   */
   /* + (91495*Pi*\[Eta]^2)/1512) + v^6*(16447322263/139708800           */
   /* - (1712*EulerGamma)/105 + (16*Pi^2)/3 - (56198689*\[Eta])/217728   */
   /* + (451*Pi^2*\[Eta])/48 + (541*\[Eta]^2)/896 - (5605*\[Eta]^3)/2592 */
   /* - (856*Log[16])/105 - (1712*Log[v])/105)))/5                       */
   /*--------------------------------------------------------------------*/

   dvalues->data[0] = 6.4*v9*eta/ak->totalmass*(1.
     + v2*(-2.2113095238095237 - 2.75*eta)
     + v3*fourpi
     + v4*(1.8795745149911816 + 6.77628968*eta
	     + 3.27777778*eta2)
     + v5*(-6.1889881*LAL_PI - 23.625*LAL_PI*eta)
     + v6*(117.72574285227559 - 16.3047619*0.577215664901
	     + 5.3333333*pisq - 258.114202*eta
	     + 9.39583333*pisq*eta + 0.603794643*eta2
	     - 2.16242284*eta3 - 8.15238095*log(16*v2))
     + v7*(-1.09499008*LAL_PI + 59.3047288*LAL_PI*eta
	     + 60.5125661*LAL_PI*eta2));

   dvalues->data[1] =  v3/ak->totalmass;
}

/*  <lalVerbatim file="LALTaylorT4WaveformCP"> */
void LALTaylorT4Waveform (
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{ /* </lalVerbatim> */

   InspiralInit paramsInit;
   UINT4        count;
   INITSTATUS(status, "LALTaylorT4Waveform", LALTAYLORT4WAVEFORMC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signalvec,  status,
	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec->data,  status,
   	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(params,  status,
   	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(params->nStartPad >= 0, status,
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->nEndPad >= 0, status,
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->fLower > 0, status,
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->tSampling > 0, status,
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->totalMass > 0., status,
   	LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   LALInspiralSetup (status->statusPtr, &(paramsInit.ak), params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &(paramsInit.func),
					 &(paramsInit.ak), params);
   CHECKSTATUSPTR(status);

   memset(signalvec->data, 0, signalvec->length * sizeof( REAL4 ));

   /* Call the engine function */
   LALTaylorT4WaveformEngine(status->statusPtr, signalvec, 
             NULL, NULL, NULL, &count, params, &paramsInit);
   CHECKSTATUSPTR( status );

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

void
LALTaylorT4WaveformForInjection(
                             LALStatus        *status,
                             CoherentGW       *waveform,
                             InspiralTemplate *params,
                             PPNParamStruc  *ppnParams
                             )
{
  UINT4        count, i;
  REAL8       phiC;

  REAL4Vector *a   = NULL;      /* pointers to generated amplitude  data */
  REAL4Vector *ff  = NULL;      /* pointers to generated  frequency data */
  REAL8Vector *phi = NULL;      /* pointer to generated phase data */

  InspiralInit paramsInit;


  INITSTATUS(status, "LALInspirallTaylorT4WaveformForInjection", LALTAYLORT4WAVEFORMC);
  ATTATCHSTATUSPTR(status);

  /* Make sure parameter and waveform structures exist. */
  ASSERT( params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  ASSERT( !( waveform->a ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->h ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->f ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->phi ), status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );

  params->ampOrder = 0;
  XLALPrintInfo( "WARNING: Amp Order has been reset to %d", params->ampOrder);

  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);

  if (paramsInit.nbins == 0)
  {
      XLALPrintError( "Error: Estimated length of injection is zero.\n" );
      ABORT( status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE );
  }

  /* Now we can allocate memory and vector for coherentGW structure*/
  ff  = XLALCreateREAL4Vector( paramsInit.nbins );
  a   = XLALCreateREAL4Vector( 2*paramsInit.nbins );
  phi = XLALCreateREAL8Vector( paramsInit.nbins );
  if ( !ff || !a || !phi )
  {
    if ( ff )
      XLALDestroyREAL4Vector( ff );
    if ( a )
      XLALDestroyREAL4Vector( a );
    if ( phi )
      XLALDestroyREAL8Vector( phi );
    ABORTXLAL( status );
  }

  /* Call the engine function */  
  LALTaylorT4WaveformEngine(status->statusPtr, NULL, a, ff,
                             phi, &count, params, &paramsInit);
  BEGINFAIL( status )
  {
    XLALDestroyREAL4Vector( ff );
    XLALDestroyREAL4Vector( a );
    XLALDestroyREAL8Vector( phi );
  }
  ENDFAIL( status );

  /*wrap the phase vector*/
  phiC =  phi->data[count-1] ;
  for (i = 0; i < count; i++)
  {
    phi->data[i] =  phi->data[i] - phiC + ppnParams->phi;
  }

  /* Allocate the waveform structures. */
  if ( ( waveform->a = (REAL4TimeVectorSeries *)
         LALCalloc(1, sizeof(REAL4TimeVectorSeries) ) ) == NULL ) 
  {
    XLALDestroyREAL4Vector( ff );
    XLALDestroyREAL4Vector( a );
    XLALDestroyREAL8Vector( phi );
    ABORT( status, LALINSPIRALH_EMEM,
           LALINSPIRALH_MSGEMEM );
  }
  if ( ( waveform->f = (REAL4TimeSeries *)
         LALCalloc(1, sizeof(REAL4TimeSeries) ) ) == NULL ) 
  {
    LALFree( waveform->a ); waveform->a = NULL;
    XLALDestroyREAL4Vector( ff );
    XLALDestroyREAL4Vector( a );
    XLALDestroyREAL8Vector( phi );
    ABORT( status, LALINSPIRALH_EMEM,
           LALINSPIRALH_MSGEMEM );
  }
  if ( ( waveform->phi = (REAL8TimeSeries *)
         LALCalloc(1, sizeof(REAL8TimeSeries) ) ) == NULL ) 
  {
    LALFree( waveform->a ); waveform->a = NULL;
    LALFree( waveform->f ); waveform->f = NULL;
    XLALDestroyREAL4Vector( ff );
    XLALDestroyREAL4Vector( a );
    XLALDestroyREAL8Vector( phi );
    ABORT( status, LALINSPIRALH_EMEM,
               LALINSPIRALH_MSGEMEM );
  }

  waveform->a->data = XLALCreateREAL4VectorSequence( (UINT4)count, 2 );
  waveform->f->data = XLALCreateREAL4Vector( count );
  waveform->phi->data = XLALCreateREAL8Vector( count );

  if ( !waveform->a->data || !waveform->f->data || !waveform->phi->data )
  {
    if ( waveform->a->data )
      XLALDestroyREAL4VectorSequence( waveform->a->data );
    if ( waveform->f->data )
      XLALDestroyREAL4Vector( waveform->f->data );
    if ( waveform->phi->data )
      XLALDestroyREAL8Vector( waveform->phi->data );
    LALFree( waveform->a ); waveform->a = NULL;
    LALFree( waveform->f ); waveform->f = NULL;
    XLALDestroyREAL4Vector( ff );
    XLALDestroyREAL4Vector( a );
    XLALDestroyREAL8Vector( phi );
    ABORTXLAL( status );
  }

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

  snprintf( waveform->a->name, LALNameLength,   "T4 inspiral amplitude" );
  snprintf( waveform->f->name, LALNameLength,   "T4 inspiral frequency" );
  snprintf( waveform->phi->name, LALNameLength, "T4 inspiral phase" );

  /* --- fill some output ---*/
  ppnParams->tc     = (double)(count-1) / params->tSampling ;
  ppnParams->length = count;
  ppnParams->dfdt   = ((REAL4)(waveform->f->data->data[count-1]
                    - waveform->f->data->data[count-2])) * ppnParams->deltaT;
  ppnParams->fStop  = params->fFinal;
  ppnParams->termCode        = GENERATEPPNINSPIRALH_EFSTOP;
  ppnParams->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;

  ppnParams->fStart   = ppnParams->fStartIn;


  /* --- free memory --- */
  XLALDestroyREAL4Vector( ff );
  XLALDestroyREAL4Vector( a );
  XLALDestroyREAL8Vector( phi );

  DETATCHSTATUSPTR(status);
  RETURN (status);
}

/*---------------------------------------------------------*/
void
LALTaylorT4WaveformEngine (
                LALStatus        *status,
                REAL4Vector      *signalvec,
                REAL4Vector      *a,
                REAL4Vector      *ff,
                REAL8Vector      *phiVec,
                UINT4            *countback,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                )
{

   void                  *funcParams;
   UINT4                 length=0, count, ndx;
   INT4                  nn=2;
   REAL8                 h, omega, omegaMax, t, dt, m, eta, phi, v;
   REAL8Vector           dummy, values, dvalues, newvalues, yt, dym, dyt;
   rk4GSLIntegrator      *integrator = NULL;
   InspiralDerivativesIn in2;
   rk4In                 in4;
   expnCoeffs            ak;
   expnFunc              func;
   /*DFindRootIn           rootIn;*/
   CHAR message[256];

   /* Variables for injection */
   REAL8 mTot = 0;
   REAL8 unitHz = 0;
   REAL8 f2a = 0;
   REAL8 mu = 0;
   REAL8 cosI = 0;/* cosine of system inclination */
   REAL8 etab = 0;
   REAL8 fFac = 0; /* SI normalization for f and t */
   REAL8 f2aFac = 0;/* factor multiplying f in amplitude function */
   REAL8 apFac = 0, acFac = 0;/* extra factor in plus and cross amplitudes */


   INITSTATUS(status, "LALTaylorT4WaveformEngine", LALTAYLORT4WAVEFORMC);
   ATTATCHSTATUSPTR(status);

   ASSERT( signalvec || ( ff && a && phiVec ), status, 
            LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );

/* Allocate all the memory required to dummy and then point the various
   arrays to dummy - this makes it easier to handle memory failures */

   dummy.length = nn * 6;

   if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * nn * 6))) {
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
   }

   values.length    = nn;
   dvalues.length   = nn;
   newvalues.length = nn;
   yt.length        = nn;
   dym.length       = nn;
   dyt.length       = nn;

   values.data    = &dummy.data[0];
   dvalues.data   = &dummy.data[nn];
   newvalues.data = &dummy.data[2*nn];
   yt.data        = &dummy.data[3*nn];
   dym.data       = &dummy.data[4*nn];
   dyt.data       = &dummy.data[5*nn];


   if ( a )
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

   /* Set dt to sampling interval specified by user */

   dt = 1./params->tSampling;
   ak   = paramsInit->ak;
   func = paramsInit->func;
   if ( signalvec )
   {
     length = signalvec->length;
   }
   else
   {
     length = ff->length;
   }
   eta = ak.eta;
   m = ak.totalmass;

   /* Begin initial conditions */
   /* Given omega compute v    */
   omega = LAL_PI*params->fLower;
   v = pow(m*omega,1./3.);
   /* End of initial conditions */

   values.data[0] = v;
   values.data[1] = phi = params->startPhase;

   /* fprintf(stdout, "Initial v=%e, phi=%e\n", values.data[0], values.data[1]); */
   t = 0.0;

   /* Initialize the GSL integrator */
   switch (params->order)
   {
     case LAL_PNORDER_TWO:
       in4.function = LALTaylorT4Derivatives4PN;
       break;
     case LAL_PNORDER_TWO_POINT_FIVE:
       in4.function = LALTaylorT4Derivatives5PN;
       break;
     case LAL_PNORDER_THREE:
       in4.function = LALTaylorT4Derivatives6PN;
       break;
     case LAL_PNORDER_THREE_POINT_FIVE:
       in4.function = LALTaylorT4Derivatives7PN;
       break;
     default:
       snprintf(message, 256, "There are no T4 waveforms at order %d\n", params->order);
       LALError( status, message );
       LALFree(dummy.data);
       ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
   }

   in4.y = &values;
   in4.h = dt;
   in4.n = nn;
   in4.yt = &yt;
   in4.dym = &dym;
   in4.dyt = &dyt;
   in4.x = 0.;

   if (!(integrator = XLALRungeKutta4Init(nn, &in4)))
   {
     LALFree(dummy.data);
     ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
   }

   in2.totalmass = ak.totalmass;
   in2.dEnergy = func.dEnergy;
   in2.flux = func.flux;
   in2.coeffs = &ak;
   funcParams = (void *) &in2;
   /* Calculate the initial value of omega */
   in4.function(&values, &dvalues, funcParams);

   /* Begin integration loop here */

   omegaMax = 1./(pow(6.,1.5)*m);
   t = 0.0;
   ndx = 0;
   count = 0;

   /* fprintf(stdout, "fMin=%e, fMax=%e, dv/dt=%e, dphi/dt=%e\n", */
   /* omega/LAL_PI, omegaMax/LAL_PI,                              */
   /* dvalues.data[0], dvalues.data[1]);                          */

   while (omega<omegaMax)
   {
      if (count >= length)
      {
        XLALRungeKutta4Free( integrator );
	LALFree(dummy.data);
	ABORT(status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
      }

      h = 4 * m * eta * v*v * cos(2.*phi);
      if ( signalvec )
      {
        signalvec->data[ndx] = h;
      }
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
          phiVec->data[count]      = (REAL8)(phi);
      }
      /* fprintf(stdout, "%e %e %e\n", t, h, omega/(2.*m*LAL_PI)); */

      /* Integrate one step forward */
      in4.dydx = &dvalues;
      in4.x = t;
      LALRungeKutta4(status->statusPtr, &newvalues, integrator, funcParams);
      BEGINFAIL( status )
      {
        XLALRungeKutta4Free( integrator );
        LALFree(dummy.data);
      }
      ENDFAIL( status );

      /* Update the values of the dynamical variables */
      v = values.data[0] = newvalues.data[0];
      phi  = values.data[1] = newvalues.data[1];

      /* Compute the derivaties at the new location */
      in4.function(&values, &dvalues, funcParams);
      omega = dvalues.data[1];

      t = (++count-params->nStartPad) * dt;
      ndx++;
   }

   /*----------------------------------------------------------------------*/
   /* Record the final cutoff frequency of BD Waveforms for record keeping */
   /* ---------------------------------------------------------------------*/
   *countback = count;
   v = pow(m*omega, 1./3.);
   params->vFinal = v;
   params->tC = t;
   params->fFinal = omega/LAL_PI;
   /* fprintf(stdout, "Final velocity=%e, time=%e, frequency=%e\n", v, t, params->fFinal); */

   XLALRungeKutta4Free( integrator );
   LALFree(dummy.data);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
