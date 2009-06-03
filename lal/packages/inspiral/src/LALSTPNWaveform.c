/*
*  Copyright (C) 2007 Drew Keppel, Gareth Jones, Craig Robinson , Thomas Cokelaer, Michele Vallisneri
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

/*<lalVerbatim file="LALSTPNWaveformCV">
  Author: Vallisneri, M.  Cokelaer, T.
  $Id$
  </lalVerbatim>  */



/*  <lalLaTeX>

\subsection{Module \texttt{LALSTPNWaveform.c}}
DOCUMENTATION IN PROGRESS

Module to generate STPN (spinning binaries) waveforms in agreement with
the injecttion  package (return a CoherentGW structure).

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALSTPNWaveformForInjectionCP}
\index{\verb&LALSTPNWaveformForInjection&}


\subsubsection*{Description}

\subsubsection*{Algorithm}


\subsubsection*{Uses}
\begin{verbatim}
   LALSTPNderivatives
\end{verbatim}

\subsubsection*{Notes}


\vfill{\footnotesize\input{LALSTPNWaveformCV}}

</lalLaTeX>  */



#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>


NRCSID (LALSTPNWAVEFORMC, "$Id$");

/* private structure with PN parameters*/

void LALSTPNderivatives(REAL8Vector *values, REAL8Vector *dvalues, void *mparams);
typedef struct LALSTPNstructparams {
	REAL8 eta;
	REAL8 m1m2;
	REAL8 m2m1;
	REAL8 wdotnew;
	REAL8 wdotorb[9];
	REAL8 wspin15;
        REAL8 wspin20;
	REAL8 LNhdot15;
	REAL8 LNhdot20;
	REAL8 S1dot15;
	REAL8 S2dot15;
        REAL8 Sdot20;
	REAL8 epnorb[9]; /*added for stopping test */
} LALSTPNparams;

/* my function to set STPN derivatives*/





void LALSTPNderivatives(REAL8Vector *values, REAL8Vector *dvalues, void *mparams) {

    REAL8 s;
    REAL8 omega;
    REAL8 LNhx,LNhy,LNhz,LNmag;
    REAL8 S1x,S1y,S1z;
    REAL8 S2x,S2y,S2z;
    REAL8 alphadotcosi;
    REAL8 tmpx,tmpy,tmpz;
    REAL8 dotLNS1,dotLNS2;
    REAL8 crossx,crossy,crossz;
    REAL8 ds,domega,domeganew;
    REAL8 dLNhx,dLNhy,dLNhz;
    REAL8 dS1x,dS1y,dS1z;
    REAL8 dS2x,dS2y,dS2z;
    REAL8 omega2;
    REAL8 test; /*Michele-041208: added test, see below*/
    REAL8 dotS1S2; /*Michele-041208: addede since it is used twice, see below*/
    LALSTPNparams *params = (LALSTPNparams*)mparams;

    REAL8 v, v2, v3, v4, v7, v11;

    /* --- computation start here --- */
    s     = values->data[0];
    omega = values->data[1];

    LNhx  = values->data[2];
    LNhy  = values->data[3];
    LNhz  = values->data[4];

    S1x  = values->data[5];
    S1y  = values->data[6];
    S1z  = values->data[7];

    S2x  = values->data[8];
    S2y  = values->data[9];
    S2z  = values->data[10];

    /* domega: could be optimized by computing omega^(1/3) and then multiplying it*/
    /* Thomas comments: in order to improve the speed I introduced the variable v and
     * replace the previous* function by a new one taken into account that variable.
     * the idea is to avoid most of the pow function which are quite slow*/


    v = pow(omega, 1.0/3.0);
    v2 = v * v;
    v3 = v2 * v;
    v4 = v3 * v;     /* effectively used farther */
    v7 = v4 * v3;    /* effectively used farther */
    v11 = v7 * v4;   /* effectively used farther */

    domeganew = params->wdotnew * v11;

    if (omega < 0.0){
       fprintf(stderr, "WARNING: Omega has become -ve, this should lead to nan's \n");
    }

    domega =
	params->wdotorb[0]
	+ v * (params->wdotorb[1]
	+ v * ( params->wdotorb[2]
	+ v * ( params->wdotorb[3]
	+ v * (	params->wdotorb[4]
	+ v * (	params->wdotorb[5]
	+ v * ( params->wdotorb[6] +  params->wdotorb[7] *  log(omega)
	+ v * ( params->wdotorb[8] ) ) ) ) ) ) ) ;

  /* Michele-041208: this evaluates the part of the test coming from the orbital energy */

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

    domega += params->wspin15 * omega *
	( LNhx * (113.0 * S1x + 113.0 * S2x + 75.0 * params->m2m1 * S1x + 75.0 * params->m1m2 * S2x) +
	  LNhy * (113.0 * S1y + 113.0 * S2y + 75.0 * params->m2m1 * S1y + 75.0 * params->m1m2 * S2y) +
	  LNhz * (113.0 * S1z + 113.0 * S2z + 75.0 * params->m2m1 * S1z + 75.0 * params->m1m2 * S2z) );

    dotLNS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
    dotLNS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);

    /* Michele-041208: added dotS1S2 since it's used twice, also in test */
    dotS1S2 = (S1x*S2x + S1y*S2y + S1z*S2z);

    domega += params->wspin20 * v4 *
	( 247.0 * dotS1S2 -
	  721.0 * dotLNS1 * dotLNS2 );

    domega *= domeganew;

    /* Michele-041208: this evaluates the part of the test coming from the spin energy terms */
    /* If more terms come in, they need to be added by hand */

    if(params->wspin15 != 0.0)
      test += -0.5 * params->eta *
	(5.0/3.0) * v2 * ( (8.0/3.0 + 2.0*params->m2m1)*dotLNS1 + (8.0/3.0 + 2.0*params->m1m2)*dotLNS2 );
    /* Michele-041208: here we multiply and divide by eta to keep parallelism
       with the formula above; hopefully the compiler will simplify*/

    if(params->wspin20 != 0.0)
          test += -v3  * (dotS1S2 - 3.0 * dotLNS1 * dotLNS2);

    /* dLN, 1.5PN*/

    omega2 = omega * omega;

    tmpx = (4.0 + 3.0*params->m2m1) * S1x + (4.0 + 3.0*params->m1m2) * S2x;
    tmpy = (4.0 + 3.0*params->m2m1) * S1y + (4.0 + 3.0*params->m1m2) * S2y;
    tmpz = (4.0 + 3.0*params->m2m1) * S1z + (4.0 + 3.0*params->m1m2) * S2z;

    dLNhx = params->LNhdot15 * omega2 * (-tmpz*LNhy + tmpy*LNhz);
    dLNhy = params->LNhdot15 * omega2 * (-tmpx*LNhz + tmpz*LNhx);
    dLNhz = params->LNhdot15 * omega2 * (-tmpy*LNhx + tmpx*LNhy);

    /* dLN, 2PN*/

    tmpx = dotLNS2 * S1x + dotLNS1 * S2x;
    tmpy = dotLNS2 * S1y + dotLNS1 * S2y;
    tmpz = dotLNS2 * S1z + dotLNS1 * S2z;

    dLNhx += params->LNhdot20 * v7 * (-tmpz*LNhy + tmpy*LNhz);
    dLNhy += params->LNhdot20 * v7 * (-tmpx*LNhz + tmpz*LNhx);
    dLNhz += params->LNhdot20 * v7 * (-tmpy*LNhx + tmpx*LNhy);

    /* dS1, 1.5PN*/

    LNmag = params->eta / v ;

    crossx = (LNhy*S1z - LNhz*S1y);
    crossy = (LNhz*S1x - LNhx*S1z);
    crossz = (LNhx*S1y - LNhy*S1x);

    dS1x = params->S1dot15 * omega2 * LNmag * crossx;
    dS1y = params->S1dot15 * omega2 * LNmag * crossy;
    dS1z = params->S1dot15 * omega2 * LNmag * crossz;

    /* dS1, 2PN*/

    tmpx = S1z*S2y - S1y*S2z;
    tmpy = S1x*S2z - S1z*S2x;
    tmpz = S1y*S2x - S1x*S2y;

    dS1x += params->Sdot20 * omega2 *
	(tmpx - 3.0 * dotLNS2 * crossx);
    dS1y += params->Sdot20 * omega2 *
	(tmpy - 3.0 * dotLNS2 * crossy);
    dS1z += params->Sdot20 * omega2 *
	(tmpz - 3.0 * dotLNS2 * crossz);

    /* dS2, 1.5PN*/

    crossx = (-LNhz*S2y + LNhy*S2z);
    crossy = (-LNhx*S2z + LNhz*S2x);
    crossz = (-LNhy*S2x + LNhx*S2y);

    dS2x = params->S2dot15 * omega2 * LNmag * crossx;
    dS2y = params->S2dot15 * omega2 * LNmag * crossy;
    dS2z = params->S2dot15 * omega2 * LNmag * crossz;

    /* dS2, 2PN*/

    dS2x += params->Sdot20 * omega2 * (-tmpx - 3.0 * dotLNS1 * crossx);
    dS2y += params->Sdot20 * omega2 * (-tmpy - 3.0 * dotLNS1 * crossy);
    dS2z += params->Sdot20 * omega2 * (-tmpz - 3.0 * dotLNS1 * crossz);

    /* -? the original equations had one additional spin term*/
    /* Michele-041208: but I think it was an experimental term introduced by Yi and Alessandra */

    /* dphi*/

    alphadotcosi = -LNhz * (LNhy*dLNhx - LNhx*dLNhy) / (LNhx*LNhx + LNhy*LNhy);
    ds = omega - alphadotcosi;

    /* copy back into dvalues structure*/

    dvalues->data[0] = ds;
    dvalues->data[1] = domega;

    dvalues->data[2] = dLNhx;
    dvalues->data[3] = dLNhy;
    dvalues->data[4] = dLNhz;

    dvalues->data[5] = dS1x;
    dvalues->data[6] = dS1y;
    dvalues->data[7] = dS1z;

    dvalues->data[8] = dS2x;
    dvalues->data[9] = dS2y;
    dvalues->data[10]= dS2z;

    dvalues->data[11] = 0;

    /* Michele-041208: not a derivative, but pass it back here anyway */

    values->data[11] = test;
}






/*  <lalVerbatim file="LALEOBWaveformCP"> */
void
LALSTPNWaveform (
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{ /* </lalVerbatim> */

   UINT4 count;
   InspiralInit paramsInit;
   INITSTATUS(status, "LALSTPNWaveform", LALSTPNWAVEFORMC);
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
   LALInspiralChooseModel(status->statusPtr, &(paramsInit.func), &(paramsInit.ak), params);
   CHECKSTATUSPTR(status);

   memset(signalvec->data, 0, signalvec->length * sizeof( REAL4 ));
   /* Call the engine function */
   LALSTPNWaveformEngine(status->statusPtr, signalvec, NULL,NULL, NULL, NULL, NULL, &count, params, &paramsInit);
   CHECKSTATUSPTR( status );

   DETATCHSTATUSPTR(status);
   RETURN(status);
}



/*  <lalVerbatim file="LALEOBWaveformCP"> */
void
LALSTPNWaveformTemplates (
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   InspiralTemplate *params
   )
{ /* </lalVerbatim> */

   UINT4 count;

   InspiralInit paramsInit;

   INITSTATUS(status, "LALSTPNWaveform", LALSTPNWAVEFORMC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signalvec1,  status,
	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec1->data,  status,
   	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec2,  status,
	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signalvec2->data,  status,
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

   /*LALInspiralSetup (status->statusPtr, &(paramsInit.ak), params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &(paramsInit.func), &(paramsInit.ak), params);
   CHECKSTATUSPTR(status);
   */
   LALInspiralInit(status->statusPtr, params, &paramsInit);


   memset(signalvec1->data, 0, signalvec1->length * sizeof( REAL4 ));
   memset(signalvec2->data, 0, signalvec2->length * sizeof( REAL4 ));

   /* Call the engine function */
   LALSTPNWaveformEngine(status->statusPtr, signalvec1, signalvec2, NULL, NULL, NULL, NULL, &count, params, &paramsInit);
   CHECKSTATUSPTR( status );

   DETATCHSTATUSPTR(status);
   RETURN(status);
}







NRCSID (LALSTPNWAVEFORMFORINJECTIONC,
"$Id$");

/*  <lalVerbatim file="LALSTPNWaveformForInjectionCP"> */
void
LALSTPNWaveformForInjection (
			     LALStatus        *status,
			     CoherentGW       *waveform,
			     InspiralTemplate *params,
			     PPNParamStruc    *ppnParams
			    )
{
  /* </lalVerbatim> */
 UINT4 count, i;

  REAL4Vector *a=NULL;/* pointers to generated amplitude  data */
  REAL4Vector *ff=NULL ;/* pointers to generated  frequency data */
  REAL8Vector *phi=NULL;/* pointer to generated phase data */
  REAL4Vector *shift=NULL;/* pointer to generated phase data */

  REAL8 phiC;/* phase at coalescence */
  InspiralInit paramsInit;

  CreateVectorSequenceIn in;

  INITSTATUS(status, "LALSRPNWaveformForInjection", LALSTPNWAVEFORMFORINJECTIONC);
  ATTATCHSTATUSPTR(status);


  /* Make sure parameter and waveform structures exist. */
  ASSERT( params, status, LALINSPIRALH_ENULL,  LALINSPIRALH_MSGENULL );
  ASSERT(waveform, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  /* Make sure waveform fields don't exist. */
  ASSERT( !( waveform->a ), status,
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->f ), status,
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->phi ), status,
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );
  ASSERT( !( waveform->shift ), status,
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );

  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);
  if (paramsInit.nbins==0)
    {
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
  LALSCreateVector(status->statusPtr, &shift, paramsInit.nbins);
  CHECKSTATUSPTR(status);


  /* By default the waveform is empty */
  memset(ff->data, 0, paramsInit.nbins * sizeof(REAL4));
  memset(a->data, 0, 2 * paramsInit.nbins * sizeof(REAL4));
  memset(phi->data, 0, paramsInit.nbins * sizeof(REAL8));
  memset(shift->data, 0, paramsInit.nbins * sizeof(REAL4));


  /* Call the engine function */
  LALSTPNWaveformEngine(status->statusPtr, NULL, NULL, a, ff, phi, shift,&count, params, &paramsInit);

  BEGINFAIL( status )
  {
     LALSDestroyVector(status->statusPtr, &ff);
     CHECKSTATUSPTR(status);
     LALSDestroyVector(status->statusPtr, &a);
     CHECKSTATUSPTR(status);
     LALDDestroyVector(status->statusPtr, &phi);
     CHECKSTATUSPTR(status);
     LALSDestroyVector(status->statusPtr, &shift);
     CHECKSTATUSPTR(status);
  }
  ENDFAIL( status );

  /* Check an empty waveform hasn't been returned */
  for (i = 0; i < phi->length; i++)
  {
    if (phi->data[i] != 0.0) break;
    if (i == phi->length - 1)
    {
      LALSDestroyVector(status->statusPtr, &ff);
      CHECKSTATUSPTR(status);
      LALSDestroyVector(status->statusPtr, &a);
      CHECKSTATUSPTR(status);
      LALDDestroyVector(status->statusPtr, &phi);
      CHECKSTATUSPTR(status);
      LALSDestroyVector(status->statusPtr, &shift);
      CHECKSTATUSPTR(status);

      DETATCHSTATUSPTR( status );
      RETURN( status );
    }
  }

  {
    phiC =  phi->data[count-1] ;


    for (i=0; i<count;i++)
      {
	phi->data[i] =  -phiC + phi->data[i] +ppnParams->phi;
      }

    /* Allocate the waveform structures. */
    if ( ( waveform->a = (REAL4TimeVectorSeries *)
	   LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
      ABORT( status, LALINSPIRALH_EMEM,
	     LALINSPIRALH_MSGEMEM );
    }
    memset( waveform->a, 0, sizeof(REAL4TimeVectorSeries) );

    if ( ( waveform->f = (REAL4TimeSeries *)
	   LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
      LALFree( waveform->a ); waveform->a = NULL;
      ABORT( status, LALINSPIRALH_EMEM,
	     LALINSPIRALH_MSGEMEM );
    }
    memset( waveform->f, 0, sizeof(REAL4TimeSeries) );

    if ( ( waveform->phi = (REAL8TimeSeries *)
	   LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
      LALFree( waveform->a ); waveform->a = NULL;
      LALFree( waveform->f ); waveform->f = NULL;
      ABORT( status, LALINSPIRALH_EMEM,
	     LALINSPIRALH_MSGEMEM );
    }
    memset( waveform->phi, 0, sizeof(REAL8TimeSeries) );

    if ( ( waveform->shift = (REAL4TimeSeries *)
	   LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
      LALFree( waveform->a ); waveform->a = NULL;
      LALFree( waveform->f ); waveform->f = NULL;
      LALFree( waveform->phi ); waveform->phi = NULL;
      ABORT( status, LALINSPIRALH_EMEM,
	     LALINSPIRALH_MSGEMEM );
    }
    memset( waveform->shift, 0, sizeof(REAL4TimeSeries) );



    in.length = (UINT4)count;
    in.vectorLength = 2;
    LALSCreateVectorSequence( status->statusPtr,
			      &( waveform->a->data ), &in );
    CHECKSTATUSPTR(status);
    LALSCreateVector( status->statusPtr,
		      &( waveform->f->data ), count);
    CHECKSTATUSPTR(status);
    LALDCreateVector( status->statusPtr,
		      &( waveform->phi->data ), count );
    CHECKSTATUSPTR(status);

    LALSCreateVector( status->statusPtr,
		      &( waveform->shift->data ), count );
    CHECKSTATUSPTR(status);




    memcpy(waveform->f->data->data , ff->data, count*(sizeof(REAL4)));
    memcpy(waveform->a->data->data , a->data, 2*count*(sizeof(REAL4)));
    memcpy(waveform->phi->data->data ,phi->data, count*(sizeof(REAL8)));
    memcpy(waveform->shift->data->data ,shift->data, count*(sizeof(REAL4)));

    waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT = waveform->shift->deltaT
      = 1./params->tSampling;

    waveform->a->sampleUnits = lalStrainUnit;
    waveform->f->sampleUnits = lalHertzUnit;
    waveform->phi->sampleUnits = lalDimensionlessUnit;
    waveform->shift->sampleUnits 	= lalDimensionlessUnit;

    waveform->position = ppnParams->position;
    waveform->psi = ppnParams->psi;


    LALSnprintf( waveform->a->name, 	LALNameLength, "STPN inspiral amplitudes" );
    LALSnprintf( waveform->f->name, 	LALNameLength, "STPN inspiral frequency" );
    LALSnprintf( waveform->phi->name, 	LALNameLength, "STPN inspiral phase" );
    LALSnprintf( waveform->shift->name, LALNameLength, "STPN inspiral polshift" );


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
  } /* end phase condition*/

  /* --- free memory --- */


  LALSDestroyVector(status->statusPtr, &ff);
  CHECKSTATUSPTR(status);
  LALSDestroyVector(status->statusPtr, &a);
  CHECKSTATUSPTR(status);
  LALDDestroyVector(status->statusPtr, &phi);
  CHECKSTATUSPTR(status);
  LALSDestroyVector(status->statusPtr, &shift);
  CHECKSTATUSPTR(status);


  DETATCHSTATUSPTR(status);
  RETURN(status);
}







/*  <lalVerbatim file="LALSTPNWaveformInjectionCP"> */
void
LALSTPNWaveformEngine (
                LALStatus        *status,
                REAL4Vector      *signalvec1,
                REAL4Vector      *signalvec2,
                REAL4Vector      *a,
                REAL4Vector      *ff,
                REAL8Vector      *phi,
		REAL4Vector      *shift,
                UINT4            *countback,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                )
  /* </lalVerbatim> */
{
  /* declare model parameters*/

  LALSTPNparams STPNparameters;
  LALSTPNparams *mparams;

  /* declare code parameters and variables*/
  INT4 		nn = 11+1;              /* number of dynamical variables*/
                                        /* Michele-041208: added one for "test" */
  UINT4 	count;                  /* integration steps performed*/
  INT4 		length;                 /* memory allocation structure*/
  INT4 		j;                      /* counter*/

  rk4In 	in4;                    /* used to setup the Runge-Kutta integration*/
  rk4GSLIntegrator *integrator;

  expnCoeffs 	ak;
  expnFunc 	func;


  REAL8Vector 	dummy, values, dvalues, newvalues, yt, dym, dyt;

  REAL8 	lengths;
  REAL4 v=0,h1,h2, amp=0;
  REAL8 	m;
  REAL8 	t;                   /* time (units of total mass)*/
  REAL8 	unitHz;
  REAL8  	dt;
  REAL8 LNhztol = 1.0e-8;

  /* declare initial values of dynamical variables*/
  REAL8 initphi;
  REAL8 initLNhx, initLNhy, initLNhz;
  REAL8 initS1x, initS1y, initS1z;
  REAL8 initS2x, initS2y, initS2z;

  /* declare dynamical variables*/
  REAL8 vphi, omega, LNhx, LNhy, LNhz, S1x, S1y, S1z, S2x, S2y, S2z;
  REAL8 test=-2;
  REAL8 alpha, omegadot;
  REAL8 f2a, apcommon;

  INITSTATUS(status, "LALSTPNWaveform", LALSTPNWAVEFORMC);
  ATTATCHSTATUSPTR(status);

  /* set parameters from InspiralTemplate structure*/

  /* I'm getting parameters from the "InspiralTemplate" structure: this is how it looks.
     typedef struct tagInspiralTemplate
     {
     REAL8 alpha;
     REAL8 alpha1;
     REAL8 alpha2;
     Approximant approximant;
     REAL8 beta;
     REAL8 chirpMass;
     REAL8 distance; ... USE AS DISTANCE (meters)
     REAL8 eccentricity;
     REAL8 eta; ... SET TO symmetric mass ratio
     REAL8 fCutoff;
     REAL8 fendBCV;
     REAL8 fFinal;
     REAL8 fLower; ... USE AS INITIAL FREQUENCY (Hz)
     INT4  ieta;
     REAL8 inclination; ... NOT USED
     INT4  level;
     REAL4 minMatch;
     REAL8 mass1; ... USE AS MASS1
     REAL8 mass2; ... USE AS MASS2
     InputMasses massChoice;
     REAL8 mu; ... SET TO reduced mass
     INT4  number;
     INT4  nStartPad;
     INT4  nEndPad;
     REAL8 OmegaS;
     REAL8 orbitTheta0; ... USE TO SET initLNh; CANNOT BE ZERO!
     REAL8 orbitPhi0; ... USE TO SET initLNh
     Order order; ... USE AS common PN order
     REAL8 psi0;
     REAL8 psi3;
     REAL8 rFinal;
     REAL8 rInitial;
     REAL8 rLightRing;
     REAL8 signalAmplitude; ... NOT USED
     INT4Vector *segmentIdVec;
     REAL8 spin1[3]; ... USE AS spin1 (not normalized) in units of mass1^2
     REAL8 spin2[3]; ... USE AS spin2 (not normalized) in units of mass2^2
     REAL8 sourceTheta; ... this is the position of the source (colatitude)
     REAL8 sourcePhi;   ... this is the position of the source (azimuth)
     REAL8 startPhase; ... NOT USED
     REAL8 startTime; ... NOT USED
     REAL8 t0;
     REAL8 t2;
     REAL8 t3;
     REAL8 t4;
     REAL8 t5;
     REAL8 tC;
     REAL8 Theta;
     REAL8 totalMass; ... SET TO mass1 + mass2
     REAL8 tSampling; ... USE AS SAMPLING AND INTEGRATION TIME
     REAL8 vFinal;
     REAL8 vInitial;
     REAL8 Zeta2;
     struct tagInspiralTemplate *next;
     struct tagInspiralTemplate *fine;
     } InspiralTemplate; */


  /* Michele-041208: in my code I used to compute some members of the
     "params" structure (such as totalMass, eta) from more basic ones (such as mass1, mass2).
     Am I to assume that all the parameter fields have already been set when "params" is passed?

     Thomas 10/12/2004 : yes ny using LALInspiralParameterCalc function called in LALInspiralInit */


 /* Make sure parameter and waveform structures exist. */
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);


  /*===================*/

 /* Compute some parameters*/

  if (paramsInit->nbins==0)
    {
      DETATCHSTATUSPTR(status);
      RETURN (status);
    }
  func = paramsInit->func;
  ak   = paramsInit->ak;


  mparams = &STPNparameters;

  /* -? I hope the units in the following are correct*/
  /* MV-050829: they are correct if the distance is given in meters;
     use (1e6 * LAL_PC_SI * params->distance) if distance is given
     in Mpc; notice the minus here */

  apcommon = -4.0 * params->mu * LAL_MRSUN_SI/(params->distance);

  /* set units*/

  m = params->totalMass * LAL_MTSUN_SI;
  unitHz = params->totalMass * LAL_MTSUN_SI * (REAL8)LAL_PI;

  /* -? dt is set from params; but Thomas' code also contained instructions to*/
  /*    set it as -1. * eta / ( params->tSampling * 5.0*LAL_MTSUN_SI*params->totalMass );*/

  /*    tSampling is in Hz, so dt is in seconds*/

  dt = 1.0/params->tSampling;

  /* -? if the integration timestep is too large, it could be set to a fraction of the total mass*/
  /* dt = 10e-3 * m;*/

  /* -- length in seconds from Newtonian formula; chirpm in seconds*/
  lengths = (5.0/256.0) * pow(LAL_PI,-8.0/3.0)
	  * pow(params->chirpMass * LAL_MTSUN_SI * params->fLower,-5.0/3.0) / params->fLower;

  length = ceil(log10(lengths/dt)/log10(2.0));
  length = pow(2,length);
  length *= 2;

  /* set initial values of dynamical variables*/

  initphi = 0.0; /* -? see code at the end; initial phase is disabled for the moment*/

  /* note that Theta0 cannot be 0.0!*/
  /*initLNhx = sin(params->orbitTheta0)*cos(params->orbitPhi0);
  initLNhy = sin(params->orbitTheta0)*sin(params->orbitPhi0);
  initLNhz = cos(params->orbitTheta0);*/

  initLNhx = sin(params->inclination);
  initLNhy = 0.;
  initLNhz = cos(params->inclination);

  initS1x = params->spin1[0] * (params->mass1 * params->mass1) / (params->totalMass * params->totalMass);
  initS1y = params->spin1[1] * (params->mass1 * params->mass1) / (params->totalMass * params->totalMass);
  initS1z = params->spin1[2] * (params->mass1 * params->mass1) / (params->totalMass * params->totalMass);

  initS2x = params->spin2[0] * (params->mass2 * params->mass2) / (params->totalMass * params->totalMass);
  initS2y = params->spin2[1] * (params->mass2 * params->mass2) / (params->totalMass * params->totalMass);
  initS2z = params->spin2[2] * (params->mass2 * params->mass2) / (params->totalMass * params->totalMass);


  dummy.length = nn * 6;

  values.length = dvalues.length = newvalues.length =
    yt.length = dym.length = dyt.length = nn;

  if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * nn * 6))) {
    ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
  }

  values.data 		= &dummy.data[0];
  dvalues.data 		= &dummy.data[nn];
  newvalues.data 	= &dummy.data[2*nn];
  yt.data 		= &dummy.data[3*nn];
  dym.data 		= &dummy.data[4*nn];
  dyt.data 		= &dummy.data[5*nn];

    /* setup coefficients for PN equations*/

  mparams->m2m1 = params->mass2 / params->mass1;
  mparams->m1m2 = params->mass1 / params->mass2;
  /* params->eta might have been set up before but just for safety, we
   * recompute it here below.*/
  params->eta =  (params->mass1 * params->mass2)/
    (params->mass1 + params->mass2)/(params->mass1 + params->mass2);
  mparams->eta = params->eta;
  mparams->wdotnew = (96.0/5.0) * params->eta;

 /* Michele-041208: I trust (but did not check), that the same numbering
     that I used for the wdotorb[] coefficients is shared by the ak
     terms, especially given that wdotorb[6] and wdotorb[7] are both
     3PN terms for me

     Thomas-10Dec04 Indeed wdotorb[6] and [7] are part of 3PN and I took that
     into account in the LALInspiralSetup function. */

  /* Michele-041208: added epnorb initialization; is newtonian = 0? */

  for (j = newtonian; j <= 8; j++) {
    mparams->wdotorb[j] = ak.ST[j];
    mparams->epnorb[j] = 0.0;
  }

  mparams->epnorb[0] = 1.0;

  /* Michele-041208: modified for epnorb initialization; perhaps the logic is not the best */


  switch (params->order){
  case     newtonian:
  case     oneHalfPN:
  case     onePN:
    for (j = params->order + 1; j <= 8; j++){
      mparams->wdotorb[j] = 0;
    }
    mparams->wspin15 	= 0.0;
    mparams->LNhdot15 	= 0.0;
    mparams->S1dot15 	= 0.0;
    mparams->S2dot15 	= 0.0;
    mparams->wspin20 	= 0.0;
    mparams->LNhdot20 	= 0.0;
    mparams->Sdot20 	= 0.0;
    mparams->epnorb[2]  = -(1.0/12.0) * (9.0 + params->eta); /* 1PN*/

    break;
  case     onePointFivePN:
    for (j = params->order + 1; j <= 8; j++){
      mparams->wdotorb[j] = 0;
    }
    mparams->wspin15 	= -(1.0/12.0);
    mparams->LNhdot15 	= 0.5;
    mparams->S1dot15 	= (4.0 + 3.0 * mparams->m2m1) / 2.0 ;
    mparams->S2dot15 	= (4.0 + 3.0 * mparams->m1m2) / 2.0 ;
    mparams->wspin20 	= 0.0;
    mparams->LNhdot20 	= 0.0;
    mparams->Sdot20 	= 0.0;
    mparams->epnorb[2]  = -(1.0/12.0) * (9.0 + params->eta); /* 1PN*/

    break;
  case     twoPN:
  case     twoPointFivePN:
  case     threePN:
  case     threePointFivePN:
    for (j = params->order +1; j <= 8; j++){
      mparams->wdotorb[j] = 0;
    }
    mparams->wspin15 	= -(1.0/12.0);
    mparams->wspin20 	= -(1.0/48.0) / params->eta ;
    mparams->LNhdot20 	= -1.5 / params->eta;
    mparams->Sdot20 	= 0.5;
    mparams->LNhdot15 	= 0.5;
    mparams->S1dot15 	= (4.0 + 3.0 * mparams->m2m1) / 2.0 ;
    mparams->S2dot15 	= (4.0 + 3.0 * mparams->m1m2) / 2.0 ;
    mparams->epnorb[2]  = -(1.0/12.0) * (9.0 + params->eta); /* 1PN*/
    mparams->epnorb[4]  = (1.0/24.0) * (-81.0 + 57.0*params->eta - params->eta*params->eta); /* 2PN*/
    if(params->order == threePN || params->order == threePointFivePN)
      mparams->epnorb[6]  = ( -(675.0/64.0) /* 3PN*/
			      +( (209323.0/4032.0) -(205.0/96.0)*LAL_PI*LAL_PI
				 -(110.0/9.0) * (-1987.0/3080.0) ) * params->eta
			      -(155.0/96.0)*params->eta*params->eta - (35.0/5184.0)*params->eta*params->eta*params->eta );
    break;
  }

  if (params->order == threePN)
    mparams->wdotorb[(int)(threePN+1)] = ak.ST[(int)(threePN+1)];
  if (params->order == threePointFivePN)
    mparams->wdotorb[8] = ak.ST[8];



  /* setup initial conditions for dynamical variables*/

  vphi = initphi;
  omega = params->fLower * unitHz;

  LNhx = initLNhx;
  LNhy = initLNhy;
  LNhz = initLNhz;

  S1x = initS1x;
  S1y = initS1y;
  S1z = initS1z;

  S2x = initS2x;
  S2y = initS2y;
  S2z = initS2z;

  /* copy everything in the "values" structure*/

  values.data[0] = vphi;
  values.data[1] = omega;

  values.data[2] = LNhx;
  values.data[3] = LNhy;
  values.data[4] = LNhz;

  values.data[5] = S1x;
  values.data[6] = S1y;
  values.data[7] = S1z;

  values.data[8] = S2x;
  values.data[9] = S2y;
  values.data[10]= S2z;

  /* setup LALRungeKutta4: parameters are
     - "newvalues" (some kind of array allocated above with "dummy")
     - "in4":
     TestFunction *function -> is set to the derivative function
     REAL8 x (time)         -> is set to the time over m
     REAL8Vector *y         -> is a pointer set to "&values": current variables
     REAL8Vector *dydx      -> is a pointer set to "&dvalues": current derivatives
     REAL8Vector *yt        -> is a pointer set to "&yt"; probably a workbuffer
     REAL8Vector *dym       -> is a pointer set to "&dym"; probably a workbuffer
     REAL8Vector *dyt       -> is a pointer set to "&dyt"; probably a workbuffer
     REAL8 h                -> is the timestep, set to dt/m
     INT4 n
     - "funcParams": void pointer to some parameters needed by derivatives */

  in4.function 	= LALSTPNderivatives;
  in4.y 	= &values;
  in4.dydx 	= &dvalues;
  in4.h 	= dt/m;
  in4.n 	= nn;
  in4.yt 	= &yt;
  in4.dym 	= &dym;
  in4.dyt 	= &dyt;

  xlalErrno = 0;
  /* Initialize GSL integrator */
  if (!(integrator = XLALRungeKutta4Init(nn, &in4)))
  {
    INT4 errNum = XLALClearErrno();
    LALFree(dummy.data);

    if (errNum == XLAL_ENOMEM)
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
    else
      ABORTXLAL( status );
  }

  /* main integration loop*/

  t = 0.0;
  count = 0;

  /* -- for the first step let's assume omegadot is positive;*/
  /*    otherwise I should call the derivative function here*/

  LALSTPNderivatives(&values,&dvalues,(void*)mparams);
  omegadot = dvalues.data[1];

  /* -? the original BCV target model had one additional stopping test*/
  /*    that I should probably reinstate*/
  /* Michele-041208: you do a do-while loop, while I had a while; we're trusting that the first
     step does not exit already, which is probably acceptable */

  /* Thomas-10Dev04: I do not remember why I did .... but not for nothing for sure.
   It might be an empty bin at beginning or a shift im time of one bin in the coherent
   code...*/

  if (a || signalvec2)
    params->nStartPad = 0; /* must be zero for templates and injection */

  do {

    /*    fprintf(stderr,"%d %d %15.12lf %15.12lf %15.12lf %15.12lf %15.12lf\n", count, length, omegadot, LNhz*LNhz, LNhztol, omega, unitHz);*/
    if ((signalvec1 && count >= signalvec1->length) || (ff && count >= ff->length))
      {
        XLALRungeKutta4Free( integrator );
	LALFree(dummy.data);
	ABORT(status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
      }

      /* From SimulateCoherentGW.h:
	 <<We therefore write the waveforms in terms of two polarization amplitudes $A_1(t)$
	 and $A_2(t)$, a single phase function $\phi(t)$, and a polarization shift $\Phi(t)$:
	 h_+(t) = A_1(t)\cos\Phi(t)\cos\phi(t) - A_2(t)\sin\Phi(t)\sin\phi(t)
	 h_\times(t) = A_1(t)\sin\Phi(t)\cos\phi(t) + A_2(t)\cos\Phi(t)\sin\phi(t)>>

	 Note: I can only do it if BCV2's Theta = 0. Hopefully that degree of freedom can
	 be recovered from the orientation of the initial angular momentum. */

      /* now setting the wave from the dynamical variables*/

      alpha = atan2(LNhy, LNhx);

      /* I don't really need i, because I can use the explicit formulae below*/
      /* i = acos(LNhz);*/

      /* -? if omega is in units of total masses, then the following conversion is correct*/
      /*    was previously f2a = pow (f2aFac * omega, 2./3.)*/
      f2a = pow(omega, 2./3.);

      /* the ff assignment is fine, since unitHz is pi M Msun [s]*/
      /* a->data[2*count]          = (REAL4)(apcommon * f2a * (-0.25) * (3.0 + cos(2*i)));*/
      /* a->data[2*count+1]        = (REAL4)(apcommon * f2a * (-cos(i)));                 */

      if (signalvec1)
	{
	  v = pow(omega, oneby3);
	  amp = params->signalAmplitude*v*v;


	  h1 = -0.5 * amp * cos(2.*vphi) * cos(2.*alpha) * (1. + LNhz*LNhz) + amp * sin(2.*vphi) * sin(2.*alpha) * LNhz;
	  *(signalvec1->data + count) = (REAL4) h1;

	  if (signalvec2)
	    {
	      h2 = -0.5 * amp * cos(2.*vphi) * sin(2.*alpha) * (1. + LNhz*LNhz) - amp * sin(2.*vphi) * cos(2.*alpha) * LNhz;
	      *(signalvec2->data + count) = (REAL4) h2;
	    }

	}
      else if (a)
	{
	  /* MV-050829: note that the two amplitudes here differ by a
             minus w.r.t. the definitions commented above. But these are the correct ones. */

	  ff->data[count]= (REAL4)(omega/unitHz);
	  a->data[2*count]          = (REAL4)(apcommon * f2a * 0.5 * (1 + LNhz*LNhz));
	  a->data[2*count+1]        = (REAL4)(apcommon * f2a * LNhz);
	  phi->data[count]          = (REAL8)(2.0 * vphi);
	  shift->data[count]        = (REAL4)(2.0 * alpha);
	}

      /* Debugging: it can be occasionally useful to store dynamical variables
	 in the waveform output structure. Keep this here.

	 ff->data[count] = omega; phi->data[count] = vphi;
	 a->data[2*count] = LNhx; a->data[2*count+1] = LNhy; shift->data[count] = LNhz; */

      /* advancing time and calling stepper; time is in units of total mass*/

      in4.x = t/m;

      /* This integrator should be replaced, eventually,*/
      /* by a variable-timestep algorithm*/

      LALRungeKutta4(status->statusPtr, &newvalues, integrator, (void*)mparams);
      CHECKSTATUSPTR(status);

      /* updating values of dynamical variables*/

      vphi  = values.data[0] = newvalues.data[0];
      omega = values.data[1] = newvalues.data[1];

      LNhx  = values.data[2] = newvalues.data[2];
      LNhy  = values.data[3] = newvalues.data[3];
      LNhz  = values.data[4] = newvalues.data[4];

      S1x   = values.data[5] = newvalues.data[5];
      S1y   = values.data[6] = newvalues.data[6];
      S1z   = values.data[7] = newvalues.data[7];

      S2x   = values.data[8] = newvalues.data[8];
      S2y   = values.data[9] = newvalues.data[9];
      S2z   = values.data[10]= newvalues.data[10];

      LALSTPNderivatives(&values,&dvalues,(void*)mparams);
      omegadot = dvalues.data[1];

      test = values.data[11];

      t = (++count - params->nStartPad) * dt;

  }
 /* Test that omega/unitHz < NYQUIST */
 while(test < 0.0 && omegadot > 0 && LNhz*LNhz < 1.0 - LNhztol && omega/unitHz < params->tSampling/2. && !(isnan(omega))) ;

 /* if code stopped since evolving quantities became nan write an error message */
 if (isnan(omega)){
    fprintf(stderr,
    "WARNING: evolving quantities have become nan. "
     "m1: %e, "
     "m2: %e, "
     "spin1x: %e, "
     "spin1y: %e, "
     "spin1z: %e, "
     "spin2x: %e, "
     "spin2y: %e, "
     "spin2z: %e, "
     "inclination: %e\n ",
      params->mass1, params->mass2,
      params->spin1[0], params->spin1[1], params->spin1[2],
      params->spin2[0], params->spin2[1], params->spin2[2],
      params->inclination);
 }
 /* if code stopped due to co-ord singularity write an error message */
 else if (!(LNhz*LNhz < 1.0 - LNhztol)){
     fprintf( stderr,
     "WARNING: Injection terminated, co-ord singularity. "
     "m1: %e, "
     "m2: %e, "
     "spin1x: %e, "
     "spin1y: %e, "
     "spin1z: %e, "
     "spin2x: %e, "
     "spin2y: %e, "
     "spin2z: %e, "
     "inclination: %e\n ",
      params->mass1, params->mass2,
      params->spin1[0], params->spin1[1], params->spin1[2],
      params->spin2[0], params->spin2[1], params->spin2[2],
      params->inclination);
 }

  /* Michele-041208: modified, added test */
  /* Michele-041208: to check, is omega really in Hz? */
  /* Michele-041208: the original Mathematica code, stopped the
     inspiral at the ISCO for order = newtonian. This is not currently implemented here. */


  /* -? the EOB version saves some final values in params; I'm doing only fFinal*/

  if (signalvec1 || signalvec2){
    params->fFinal = pow(v,3.)/(LAL_PI*m);
  }
  else if (a){
    params->fFinal =ff->data[count-1];
  }

  /* -? this looks like the final phase is being subtracted from the phase history*/
  /*    this operations negates the use of initphi, which is set to 0.0 anyway*/
  /* -? I will comment this out to compare with my Mathematica code*/

     if (signalvec1 && !signalvec2) params->tC = t;
   *countback = count;


  XLALRungeKutta4Free( integrator );
  LALFree(dummy.data);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

