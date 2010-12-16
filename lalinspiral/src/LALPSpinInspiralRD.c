/*
*  Copyright (C) 2010 Riccardo Sturani, based on LALEOBWaveform.c by
*  Stas Babak, David Churches, Duncan Brown, David Chin, Jolien Creighton,
*  B.S. Sathyaprakash, Craig Robinson , Thomas Cokelaer, Evan Ochsner
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

/****  <lalVerbatim file="LALPSpinInspiralRDCV">
 * $Id$
 **** </lalVerbatim> */

/****  <lalLaTeX>
 *
 * \subsection{Module \texttt{LALPSpinInspiralRD.c},
 * \texttt{LALPSpinInspiralTemplates} and \texttt{LALPSpinInspiralForInjection}} 
 * \label{ss:LALPSpinInspiralRD.c}
 *
 * Module to generate generic spinning binaries waveforms complete with ring-down 
 * 
 * \subsubsection*{Prototypes}
 * \vspace{0.1in}
 * \input{LALPSpinInspiralRDCP}
 * \idx{\verb&LALPSpinInspiralRD()&}
 * \begin{description}
 * \item {\tt signalvec:} Output containing the inspiral waveform.
 * \item {\tt params:} Input containing binary chirp parameters.
 * \end{description}
 *
 * \input{LALPSpinInspiralRDTemplatesCP}
 * \idx{\verb&LALPSpinInspiralRDTemplates()&}
 * \begin{description}
 * \item {\tt signalvec1:} Output containing the $+$ inspiral waveform.
 * \item {\tt signalvec2:} Output containing the $\times$ inspiral waveform.
 * \item {\tt params:} Input containing binary chirp parameters.
 * \end{description}
 *
 * \input{LALPSpinInspiralRDInjectionCP}
 * \idx{\verb&LALPSpinInspiralRDInjection()&}
 * \begin{description}
 * \item {\tt signalvec:} Output containing the inspiral waveform.
 * \item {\tt params:} Input containing binary chirp parameters.
 * \end{description}
 *
 * \subsubsection*{Description}
 * This codes provide complete waveforms for generically spinning binary systems.
 * In order to construct the waveforms three phases are joined together: 
 * an initial inspiral phase, a phenomenological phase encompassing the description
 * of the merger and the ring-down of the final black hole.
 * During the inspiral phase the system is evolved according to the standard
 * PN formulas, valid up to 3.5PN for the orbital motion, 
 * to 2.5PN level for spin-orbital momentum and to 2PN for spin-spin contributions
 * to the orbital phase.
 * Then a phenomenological phase is added during which the frequency of the 
 * waveform has a pole-like behaviour. The stitching is performed in order to 
 * ensure continuity of the phase, the frequency and its first and second 
 * derivatives. Finally a ring-down phase is attached.
 * 
 * \subsubsection*{Algorithm}
 *
 * \subsubsection*{Uses}
 * \begin{verbatim}
 * LALPSpinInspiralRDderivatives
 * LALInspiralSetup
 * LALInspiralChooseModel
 * LALRungeKutta4
 * \end{verbatim}
 *
 * \subsubsection*{Notes}
 * 
 * \vfill{\footnotesize\input{LALPSpinInspiralRDCV}}
 * 
 **** </lalLaTeX>  */

/** \defgroup psird Complete phenomenological spin-inspiral waveforms
 * 
 * This code provides complete waveforms for generically spinning binary 
 * systems.
 *
 */

#include <lal/Units.h>
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/NRWaveInject.h>
#include <lal/RealFFT.h>

NRCSID (LALPSPININSPIRALRDC, "$Id$");

typedef struct LALPSpinInspiralRDstructparams {
  REAL8 eta;                  ///< symmetric mass ratio
  REAL8 dm;                   ///< \f$m_1-m_2\f$
  REAL8 m1m2;                 ///< \f$m_1/m_2\f$
  REAL8 m2m1;                 ///< \f$m_2/m_1\f$
  REAL8 m2m;
  REAL8 m1m;
  REAL8 wdotorb[8];           ///< Coefficients of the analytic PN expansion of \f$\dot\omega_{orb}\f$
  REAL8 wdotorblog;           ///< Log coefficient of the PN expansion of of \f$\dot\omega_{orb}\f$
  REAL8 wdotspin15S1LNh;
  REAL8 wdotspin15S2LNh;
  REAL8 wdotspin20S1S2;
  REAL8 wdotspin20S1S1;      ///< Coeff. of the \f$s_1s_1\f$ cntrb. to \f$\dot\omega\f$
  REAL8 wdotspin20S1S2LNh;
  REAL8 wdotspin25S1LNh;
  REAL8 wdotspin25S2LNh;     ///< Coeff. of the \f$s_2\cdot \hat L_N\f$ cntrb. to \f$\dot\omega\f$
  REAL8 S1dot15;
  REAL8 S2dot15;
  REAL8 Sdot20;
  REAL8 S1dot25;
  REAL8 S2dot25;
  REAL8 LNhdot15;
  REAL8 LNhdot20;
  REAL8 epnorb[4];           ///< Coefficients of the PN expansion of the energy
  REAL8 epnspin15S1dotLNh;   ///< Coeff. of the \f$S_1\cdot L\f$ term in energy
  REAL8 epnspin15S2dotLNh;   ///< Coeff. of the \f$S_2\cdot L\f$ term in energy
  REAL8 epnspin20S1S2;       ///< Coeff. of the \f$S_1\cdot S_2\f$ term in energy
  REAL8 epnspin20S1S2dotLNh; ///< Coeff. of the \f$S_{1,2}\cdot L\f$ term in energy
  REAL8 epnspin20S1S1;       ///< Coeff. of the \f$S_1\cdot S_1\f$ term in energy
  REAL8 epnspin20S1S1dotLNh; 
  REAL8 epnspin20S2S2;       ///< Coeff. of the \f$S_2\cdot S_2\f$ term in energy
  REAL8 epnspin20S2S2dotLNh;
  REAL8 epnspin25S1dotLNh;
  REAL8 epnspin25S2dotLNh;

} LALPSpinInspiralRDparams;

/* 
 *
 * function to set derivatives: values and mparams input, dvalues output
 *
 */

/**
 * \ingroup psird
 * \brief Module to compute detivative of dynamical variables 
 */
void LALPSpinInspiralRDderivatives(REAL8Vector *values, REAL8Vector *dvalues, void *mparams ) {

  REAL8 Phi;                    // half of the main GW phase, this is \f$Phi\f$ of eq.3.11 of arXiv:0810.5336
  REAL8 omega;                  // time-derivative of the orbital phase
  REAL8 omega2;                 // omega squared
  REAL8 LNhx,LNhy,LNhz;         // orbital angolar momentum unit vector
  REAL8 S1x,S1y,S1z;            // dimension-less spin variable S/M^2
  REAL8 S2x,S2y,S2z;
  REAL8 alphadotcosi;           // alpha is the right ascension of L, i(iota) the angle between L and J 
  REAL8 LNhS1,LNhS2;            // scalar products
  REAL8 domega;                 // derivative of omega
  REAL8 dLNhx,dLNhy,dLNhz;      // derivatives of \f$\hat L_N\f$ components
  REAL8 dS1x,dS1y,dS1z;         // derivative of \f$S_i\f$
  REAL8 dS2x,dS2y,dS2z;
  REAL8 S1S2,S1S1,S2S2;         // Scalar products
  REAL8 energy;
  
  LALPSpinInspiralRDparams *params = (LALPSpinInspiralRDparams*)mparams;

  REAL8 v, v2, v3, v4, v5, v6, v7, v11;                                        // support variables
  REAL8 tmpx,tmpy,tmpz,cross1x,cross1y,cross1z,cross2x,cross2y,cross2z,LNhxy;  // support variables

  /* --- computation start here --- */
  Phi   = values->data[0];
  omega = values->data[1];

  if (omega < 0.0){
    fprintf(stderr, "WARNING: Omega has become -ve, this should lead to nan's \n");
  }

  LNhx  = values->data[2];
  LNhy  = values->data[3];
  LNhz  = values->data[4];
  
  S1x   = values->data[5];
  S1y   = values->data[6];
  S1z   = values->data[7];
  
  S2x   = values->data[8];
  S2y   = values->data[9];
  S2z   = values->data[10];

  v = pow(omega, 1.0/3.0);
  v2 = v * v;
  v3 = v2 * v;
  v4 = v2 * v2;
  v5 = v3 * v2;
  v6 = v4 * v2;
  v7 = v5 * v2;
  v11= v7 * v4;

  // Omega derivative without spin effects up to 3.5 PN
  // this formula does not include the 1.5PN shift mentioned in arXiv:0810.5336, five lines below (3.11)
  domega =
    params->wdotorb[0]
    + v * (params->wdotorb[1]
	   + v * ( params->wdotorb[2]
		   + v * ( params->wdotorb[3]
			   + v * ( params->wdotorb[4]
				   + v * ( params->wdotorb[5]
					   + v * ( params->wdotorb[6] +  params->wdotorblog *  log(omega)
						   + v * params->wdotorb[7] ) ) ) ) ) );

  // energy=-eta/2 * v^2 * [1-(9+\eta)/12 v^2 +...] up to 3PN without spin effects
  energy = ( 1. + v2 * ( params->epnorb[1] +
			 v2 * ( params->epnorb[2] + 
				v2 * params->epnorb[3] ) ) );

  // Adding spin effects to omega and energy
  // L dot S1,2
  LNhS1 = (LNhx*S1x + LNhy*S1y + LNhz*S1z);
  LNhS2 = (LNhx*S2x + LNhy*S2y + LNhz*S2z);

  // wdotspin15SiLNh = -1/12 (...)
  domega+= v3 * ( params->wdotspin15S1LNh * LNhS1 + params->wdotspin15S2LNh * LNhS2 );      // see e.g. Buonanno et al. gr-qc/0211087
  energy+= v3 * ( params->epnspin15S1dotLNh * LNhS1 + params->epnspin15S2dotLNh * LNhS2 );  // see e.g. Blanchet et al. gr-qc/0605140

  // wdotspin20S1S1 = -1/48 eta
  S1S1= (S1x*S1x + S1y*S1y + S1z*S1z);
  S2S2= (S2x*S2x + S2y*S2y + S2z*S2z);
  S1S2 = (S1x*S2x + S1y*S2y + S1z*S2z);
  domega+= params->wdotspin20S1S2 * v4 * ( 247.0 * S1S2 - 721.0 * LNhS1 * LNhS2 );                // see e.g. Buonanno et al. arXiv:0810.5336
  domega+= params->wdotspin20S1S1 * v4 * (719.*( LNhS1*LNhS1 + LNhS2*LNhS2 ) - 233.*(S1S1+S2S2)); // see Racine et al. arXiv:0812.4413

  energy+= v4 * ( params->epnspin20S1S2*S1S2 + params->epnspin20S1S2dotLNh * LNhS1 * LNhS2);      // see e.g. Buonanno et al. as above
  energy+= v4 * ( params->epnspin20S1S1*S1S1 + params->epnspin20S2S2*S2S2 + params->epnspin20S1S1dotLNh * LNhS1*LNhS1 + params->epnspin20S2S2 * LNhS2*LNhS2 );      // see Racine et al. as above
  
  // wdotspin25SiLNh = see below
  domega+= v5 * ( params->wdotspin25S1LNh * LNhS1 + params->wdotspin25S2LNh * LNhS2);       //see (8.3) of Blanchet et al.  
  energy+= v5 * ( params->epnspin25S1dotLNh * LNhS1 + params->epnspin25S2dotLNh * LNhS2 );  //see (7.9) of Blanchet et al.

  // Setting the right pre-factor 
  omega2 = omega * omega;
  domega *= 96./5. * params->eta * v5 * omega2;

  energy *= params->epnorb[0] * v2;

  /*Derivative of the angular momentum and spins*/

  /* dS1, 1.5PN*/
  /* S1dot15= (4+3m2/m1)/2 * eta*/

  cross1x = (LNhy*S1z - LNhz*S1y);
  cross1y = (LNhz*S1x - LNhx*S1z);
  cross1z = (LNhx*S1y - LNhy*S1x);

  cross2x = ( LNhy*S2z - LNhz*S2y );
  cross2y = ( LNhz*S2x - LNhx*S2z );
  cross2z = ( LNhx*S2y - LNhy*S2x ); 
  
  dS1x = params->S1dot15 * v5 * cross1x;
  dS1y = params->S1dot15 * v5 * cross1y;
  dS1z = params->S1dot15 * v5 * cross1z;
  
  /* dS1, 2PN*/
  /* Sdot20= 0.5*/
  tmpx = S1z*S2y - S1y*S2z;
  tmpy = S1x*S2z - S1z*S2x;
  tmpz = S1y*S2x - S1x*S2y;

  // S1S2 contribution
  dS1x += params->Sdot20 * omega2 * (tmpx - 3. * LNhS2 * cross1x);
  dS1y += params->Sdot20 * omega2 * (tmpy - 3. * LNhS2 * cross1y);
  dS1z += params->Sdot20 * omega2 * (tmpz - 3. * LNhS2 * cross1z);
  // S1S1 contribution
  dS1x -= 3. * params->Sdot20 * omega2 * LNhS1 * cross1x * (1. + params->m2m1) * params->m2m;
  dS1y -= 3. * params->Sdot20 * omega2 * LNhS1 * cross1y * (1. + params->m2m1) * params->m2m;
  dS1z -= 3. * params->Sdot20 * omega2 * LNhS1 * cross1z * (1. + params->m2m1) * params->m2m;

  // dS1, 2.5PN, eq. 7.8 of Blanchet et al. gr-qc/0605140
  // S1dot25= 9/8-eta/2.+eta+mparams->eta*29./24.+mparams->m1m2*(-9./8.+5./4.*mparams->eta)
  dS1x += params->S1dot25 * v7 * cross1x;
  dS1y += params->S1dot25 * v7 * cross1y;
  dS1z += params->S1dot25 * v7 * cross1z;
  
  /* dS2, 1.5PN*/
  /*Attenzione all'errore qui!*/
  dS2x = params->S2dot15 * v5 * cross2x;
  dS2y = params->S2dot15 * v5 * cross2y;
  dS2z = params->S2dot15 * v5 * cross2z;

  /* dS2, 2PN*/
  dS2x += params->Sdot20 * omega2 * (-tmpx - 3.0 * LNhS1 * cross2x);
  dS2y += params->Sdot20 * omega2 * (-tmpy - 3.0 * LNhS1 * cross2y);
  dS2z += params->Sdot20 * omega2 * (-tmpz - 3.0 * LNhS1 * cross2z);
  // S2S2 contribution
  dS2x -= 3. * params->Sdot20 * omega2 * LNhS2 * cross2x * (1. + params->m1m2) * params->m1m;
  dS2y -= 3. * params->Sdot20 * omega2 * LNhS2 * cross2y * (1. + params->m1m2) * params->m1m;
  dS2z -= 3. * params->Sdot20 * omega2 * LNhS2 * cross2z * (1. + params->m1m2) *params->m1m;

  // dS2, 2.5PN, eq. 7.8 of Blanchet et al. gr-qc/0605140
  dS2x += params->S2dot25 * v7 * cross2x;
  dS2y += params->S2dot25 * v7 * cross2y;
  dS2z += params->S2dot25 * v7 * cross2z;

  dLNhx = -( dS1x + dS2x ) * v / params->eta;
  dLNhy = -( dS1y + dS2y ) * v / params->eta;
  dLNhz = -( dS1z + dS2z ) * v / params->eta;

  /* dphi*/
  LNhxy= LNhx*LNhx + LNhy*LNhy;

  if (LNhxy > 0.0) alphadotcosi = LNhz * (LNhx*dLNhy - LNhy*dLNhx) / LNhxy;
  else alphadotcosi = 0.;

  /* dvalues->data[0] is the phase derivative*/
  /* omega is the derivative of the orbital phase omega \neq dvalues->data[0]*/

  dvalues->data[0] = omega - alphadotcosi;
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

  dvalues->data[11]= 0.;

  // energy value is stored
  values->data[11] = energy;
  
} /* end of LALPSpinInspiralRDderivatives*/


/**
 * \ingroup psird
 * \brief Main module to produce waveforms 
 */
void LALPSpinInspiralRD (
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{

   UINT4 count;
   InspiralInit paramsInit;
   INITSTATUS(status, "LALPSpinInspiralRD", LALPSPININSPIRALRDC);
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
   LALPSpinInspiralRDEngine(status->statusPtr, signalvec, NULL,NULL, NULL, NULL, NULL, &count, params, &paramsInit);
   CHECKSTATUSPTR( status );

   DETATCHSTATUSPTR(status);
   RETURN(status);
}


NRCSID (LALPSPININSPIRALRDTEMPLATESC,"$Id$");

/**
 * \ingroup psird
 * \brief Module to produce waveform templates 
 */
void LALPSpinInspiralRDTemplates (
   LALStatus        *status,
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   InspiralTemplate *params
   )
{
   UINT4 count;

   InspiralInit paramsInit;

   INITSTATUS(status, "LALPSpinInspiralRDTemplates", LALPSPININSPIRALRDTEMPLATESC);
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

   LALInspiralSetup (status->statusPtr, &(paramsInit.ak), params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &(paramsInit.func), &(paramsInit.ak), params);
   CHECKSTATUSPTR(status);

   memset(signalvec1->data, 0, signalvec1->length * sizeof( REAL4 ));
   memset(signalvec2->data, 0, signalvec2->length * sizeof( REAL4 ));

   LALPSpinInspiralRDEngine(status->statusPtr, signalvec1, signalvec2, NULL, NULL, NULL, NULL, &count, params, &paramsInit);
   CHECKSTATUSPTR( status );

   DETATCHSTATUSPTR(status);
   RETURN(status);
}


NRCSID (LALPSPININSPIRALRDINJECTIONC,"$Id$");

/**
 * \ingroup psird
 * \brief Module to produce injection waveforms  
 */
void LALPSpinInspiralRDForInjection (
			     LALStatus        *status,
			     CoherentGW       *waveform,
			     InspiralTemplate *params,
			     PPNParamStruc    *ppnParams
			    )
{
  UINT4 count,i;

  REAL4Vector *hh=NULL;/* pointers to generated amplitude  data */
  REAL4Vector *ff=NULL ;/* pointers to generated  frequency data */
  REAL8Vector *phi=NULL;/* pointer to generated phase data */
  REAL4Vector *alpha=NULL;/* pointer to generated phase data */

  InspiralInit paramsInit;
  UINT4 nbins;


  INITSTATUS(status, "LALPSpinInspiralRDInjection", LALPSPININSPIRALRDINJECTIONC);
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
  ASSERT( !( waveform->h ), status,
  	LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL );

  /* Compute some parameters*/
  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);
  if (paramsInit.nbins==0)
    {
      DETATCHSTATUSPTR(status);
      RETURN (status);
    }

  nbins=2*paramsInit.nbins;

  /* Now we can allocate memory and vector for coherentGW structure*/
  ff   = XLALCreateREAL4Vector(   nbins);
  hh   = XLALCreateREAL4Vector( 2*nbins);
  phi = XLALCreateREAL8Vector(  nbins);
  alpha = XLALCreateREAL4Vector( nbins);

  /* By default the waveform is empty */
  memset(ff->data, 0, nbins * sizeof(REAL4));
  memset(hh->data, 0, 2 * nbins * sizeof(REAL4));
  memset(phi->data, 0, nbins * sizeof(REAL8));
  memset(alpha->data, 0, nbins * sizeof(REAL4));

  /* Call the engine function */

  LALPSpinInspiralRDEngine(status->statusPtr, NULL, NULL, hh, ff, phi, alpha,&count, params, &paramsInit);

  BEGINFAIL( status )
  {
     XLALDestroyREAL4Vector(ff);
     XLALDestroyREAL4Vector(hh);
     XLALDestroyREAL8Vector(phi);
     XLALDestroyREAL4Vector(alpha);
  }
  ENDFAIL( status );

  /* Check an empty waveform hasn't been returned */
  for (i = 0; i < phi->length; i++)
  {
    if (phi->data[i] != 0.0) break;
    if (i == phi->length - 1)
    {
      XLALDestroyREAL4Vector(ff);
      XLALDestroyREAL4Vector(hh);
      XLALDestroyREAL8Vector(phi);
      XLALDestroyREAL4Vector(alpha);
    }
  }

  /* Allocate the waveform structures. */
  if ( ( waveform->h = (REAL4TimeVectorSeries *)
	 LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
    ABORT( status, LALINSPIRALH_EMEM,
	   LALINSPIRALH_MSGEMEM );
  }
  memset( waveform->h, 0, sizeof(REAL4TimeVectorSeries) );

  if ( ( waveform->f = (REAL4TimeSeries *)
	 LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
    LALFree( waveform->h ); waveform->a = NULL;
    ABORT( status, LALINSPIRALH_EMEM,
	   LALINSPIRALH_MSGEMEM );
  }
  memset( waveform->f, 0, sizeof(REAL4TimeSeries) );

  if ( ( waveform->phi = (REAL8TimeSeries *)
	 LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
    LALFree( waveform->h ); waveform->h = NULL;
    LALFree( waveform->f ); waveform->f = NULL;
    ABORT( status, LALINSPIRALH_EMEM,
	   LALINSPIRALH_MSGEMEM );
  }
  memset( waveform->phi, 0, sizeof(REAL8TimeSeries) );

  if ( ( waveform->shift = (REAL4TimeSeries *)
	 LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
    LALFree( waveform->h ); waveform->h = NULL;
    LALFree( waveform->f ); waveform->f = NULL;
    LALFree( waveform->phi ); waveform->phi = NULL;
    ABORT( status, LALINSPIRALH_EMEM,
	   LALINSPIRALH_MSGEMEM );
  }
  memset( waveform->shift, 0, sizeof(REAL4TimeSeries) );

  waveform->h->data=XLALCreateREAL4VectorSequence(count,2);
  waveform->f->data=XLALCreateREAL4Vector(count);
  waveform->phi->data=XLALCreateREAL8Vector(count);
  waveform->shift->data=XLALCreateREAL4Vector(count);

  memcpy(waveform->f->data->data ,     ff->data ,    count*(sizeof(REAL4)));
  memcpy(waveform->h->data->data ,     hh->data ,  2*count*(sizeof(REAL4)));
  memcpy(waveform->phi->data->data ,   phi->data ,   count*(sizeof(REAL8)));
  memcpy(waveform->shift->data->data , alpha->data , count*(sizeof(REAL4)));

  waveform->h->deltaT = waveform->f->deltaT = waveform->phi->deltaT = waveform->shift->deltaT
    = 1./params->tSampling;

  waveform->h->sampleUnits = lalStrainUnit;
  waveform->f->sampleUnits = lalHertzUnit;
  waveform->phi->sampleUnits = lalDimensionlessUnit;
  waveform->shift->sampleUnits 	= lalDimensionlessUnit;

  waveform->position = ppnParams->position;
  waveform->psi = ppnParams->psi;

  snprintf( waveform->h->name,     LALNameLength, "PSpinInspiralRD amplitudes" );
  snprintf( waveform->f->name,     LALNameLength, "PSpinInspiralRD frequency" );
  snprintf( waveform->phi->name,   LALNameLength, "PSpinInspiralRD phase" );
  snprintf( waveform->shift->name, LALNameLength, "PSpinInspiralRD alpha" );

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

  XLALDestroyREAL4Vector(ff);
  XLALDestroyREAL4Vector(hh);
  XLALDestroyREAL8Vector(phi);
  XLALDestroyREAL4Vector(alpha);

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* End LALPSpinInspiralRDForInjection */

void LALPSpinInspiralRDFreqDom (
				LALStatus        *status,
				REAL4Vector      *signalvec,
				InspiralTemplate *params
				)
{

  REAL4Vector *tsignalvec = NULL;
  REAL4Vector *fsignalvec = NULL;
  REAL4FFTPlan *forwPlan = NULL;

  InspiralInit paramsInit;

  REAL4 mod,ph;

  UINT4 count,nbins,sub;
  UINT4 i,j,iperiod;

  INITSTATUS(status, "LALPSpinInspiralRDFReqDom", LALPSPININSPIRALRDC);
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

  LALInspiralInit(status->statusPtr, params, &paramsInit);
  CHECKSTATUSPTR(status);

  if (paramsInit.nbins==0)
    {
      DETATCHSTATUSPTR(status);
      RETURN (status);
    }

  nbins=paramsInit.nbins;

  if (nbins<signalvec->length) nbins=signalvec->length;

  tsignalvec   = XLALCreateREAL4Vector( nbins);
  fsignalvec   = XLALCreateREAL4Vector( nbins);

  memset(signalvec->data, 0, signalvec->length * sizeof( REAL4 ));
  memset(tsignalvec->data, 0, nbins * sizeof(REAL4));
  memset(fsignalvec->data, 0, nbins * sizeof(REAL4));

  /* Call the engine function */

  LALPSpinInspiralRDEngine(status->statusPtr, tsignalvec, NULL,NULL, NULL, NULL, NULL, &count, params, &paramsInit);
  CHECKSTATUSPTR( status );

  iperiod=(INT4)(1./params->tSampling/params->fLower);
  for (i=0;i<iperiod;i++)
    tsignalvec->data[i]*=exp(7.*(((REAL4)(i))/((REAL4)(iperiod))-1.));

  REAL8 norm;
  norm=0.;
  for (i=0;i<tsignalvec->length;i++) {
    norm+=tsignalvec->data[i]*tsignalvec->data[i];
  }

  forwPlan = XLALCreateForwardREAL4FFTPlan(nbins, 0);
  if (forwPlan == NULL) ABORTXLAL(status);

  XLALREAL4VectorFFT(fsignalvec, tsignalvec, forwPlan);
  XLALDestroyREAL4Vector(tsignalvec);
  XLALDestroyREAL4FFTPlan(forwPlan);  
 
  sub=nbins/signalvec->length;

  mod=0.;
  ph=0.;
  norm=0.;
  j=1;
  for (i=1;i<nbins/2;i++) {
    mod+=sqrt(fsignalvec->data[i]*fsignalvec->data[i]+fsignalvec->data[nbins-i]*fsignalvec->data[nbins-i]);
    ph+=atan2(fsignalvec->data[nbins-i],fsignalvec->data[i]);
    if (i%sub==0) {
      mod/=(REAL4)(sub);
      ph/=(REAL4)(sub);
      signalvec->data[j]=mod*cos(ph);
      signalvec->data[nbins-j]=mod*sin(ph);
      norm+=2.*mod*mod;
      mod=0.;
      ph=0.;
      j++;
    }
  }
  signalvec->data[0]=0.;
  if (i%sub==0)
    signalvec->data[1]=(mod+fsignalvec->data[1])/((REAL4)(sub));
  else {
    signalvec->data[1]=fsignalvec->data[1];
  }
  norm+=signalvec->data[0]*signalvec->data[0]+signalvec->data[1]*signalvec->data[1];

  XLALDestroyREAL4Vector(fsignalvec);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


/*
 *
 * Main function
 *
 */

NRCSID (LALPSPININSPIRALRDENGINEC,"$Id$");

/**
 * \ingroup psird
 * \brief Module actually computing PSIRD waveforms  
 */
void LALPSpinInspiralRDEngine (
                LALStatus        *status,
                REAL4Vector      *signalvec1,
                REAL4Vector      *signalvec2,
                REAL4Vector      *hh,
                REAL4Vector      *ff,
                REAL8Vector      *phi,
                REAL4Vector      *shift,
                UINT4            *countback,
                InspiralTemplate *params,
		InspiralInit     *paramsInit
                )

{

  /* declare model parameters*/
  LALPSpinInspiralRDparams PSIRDparameters;
  
  /* declare code parameters and variables*/
  INT4		nn = 11+1;              // number of dynamical variables, the extra one is the energy
  UINT4 	count,apcount,apcountM,write;    // integration steps performed 
  UINT4		estlength;              // maximum signal vector length
  UINT4 	length;                 // signal vector length
  UINT4		i,j,k,l;                // counters          
  UINT4         subsampling=1;          // multiplies the rate           

  rk4In 	in4;                    // used to setup the Runge-Kutta integration
  rk4GSLIntegrator *integrator;

  expnCoeffs 	ak;                     //Coefficients in a generic PN expansion (E, flux...)

  REAL8Vector 	dummy, values, dvalues, newvalues, yt, dym, dyt;
  // support variables

  REAL8 	lengths;                //length in seconds
  REAL4         v=0.;
  REAL4         v2=0.;
  REAL8 	m;                      // Total mass in SI units
  REAL8 	t;                      // time (units of total mass)
  REAL8 	unitHz;
  REAL8  	dt;
  REAL8         LNhztol = 1.0e-8;

  /* declare initial values of dynamical variables */
  REAL8 initPhi,initomega,initomwrite,initv;
                      
  REAL8 initLNh[3],initS1[3],initS2[3];
  REAL8 iJ[3],iJmod,initJh[3];
  REAL8 iS1[3],iS2[3];
  REAL8 ry[3][3],rz[3][3];
  REAL8 thetaJ,phiJ;

  REAL8 ci,si,ci2,si2,c2i,s2i,c2i2,s2i2,c4i2,s4i2,c5i2,s5i2,c6i2,s6i2;
  REAL8 amp22,amp20,amp33,amp22ini;
  /* Useful variables to define GW multipolar modes*/

  REAL4Vector  *h2P2, *h2M2, *h2P1, *h2M1, *h20, *h3P3, *h3M3;
  REAL4Vector  *sig1, *sig2;
  REAL4Vector  *hap, *fap, *shift22;
  REAL8Vector  *phap;

  /* support variables*/

  /* declare dynamical variables*/
  REAL8 Phi, omega, LNhx, LNhy, LNhz, LNhxy, S1x, S1y, S1z, S2x, S2y, S2z;
  REAL8 S1dotL,S2dotL,S1dotS1,S2dotS2,S1dotS2;
  REAL8 dLNhx,dLNhy,dLNhz=0.;
  REAL8 Psi=0.;
  REAL8 omegaold=0.;
  REAL8 omegadot=0.;
  REAL8 omegadotold;
  REAL8 omegaddot;
  REAL8 alpha=0.;
  REAL8 iota=0.;
  REAL8 alphadot=0.;
  REAL8 iotadot=0.;
  REAL8 iotadot0,iotadot1,alphadot0,alphadot1;
  REAL8 iotaddot,iotadotold,alphaddot,alphadotold;

  REAL4 mchi1,mchi2;
  REAL4 cosa1,cosa2,cosg;
  REAL4 LNmag;

  REAL8 energy=0.;
  REAL8 enold;
  REAL8 v2old;

  LALPSpinInspiralRDparams *mparams;

  //  CHAR message[256];

  /* used in the phen-phase attach procedure*/
  REAL8           tAs,om1,om0;
  REAL8           t0,Psi0,alpha0,iota0;

  INT4 errcode,errcode2;
  REAL8 finalMass, finalSpin;
  
  /* used in the ring-down attach procedure*/
  REAL8           omegamatch,omegaRD=0.;
  COMPLEX8Vector  *modefreqs;
  INT4            xlalSt2P2,xlalSt2M2=0;
  UINT4           nmodes;
  COMPLEX16       MultSphHarm2P2;     // Spin-weighted spherical harmonics 2,2
  COMPLEX16       MultSphHarm2M2;     // Spin-weighted spherical harmonics 2,-2
  COMPLEX16       MultSphHarm2P1;     // Spin-weighted spherical harmonics 2,1
  COMPLEX16       MultSphHarm2M1;     // Spin-weighted spherical harmonics 2,-1
  COMPLEX16       MultSphHarm20;      // Spin-weighted spherical harmonics 2,0
  COMPLEX16       MultSphHarm3P3;     // Spin-weighted spherical harmonics 3,3
  COMPLEX16       MultSphHarm3M3;     // Spin-weighted spherical harmonics 3,-3
  REAL4           x0, x1, x2, x3;

  REAL4 inc;
  REAL4 fracRD;
  REAL4 sqrtOneMinus4Eta;

  /* switch to keep track of matching of the linear frequency growth phase*/
  INT4 rett=0;

  const double omM0   =  0.0595;
  const double omMz1p2 = -5.02e-3;
  const double omM12   = -4.29e-4;
  const double omMsq   =  5.78e-3;
  const double omMz12  =  2.66e-3;
  const double omMzsq  = -9.27e-3;

  const double frac0   =  0.57;
  const double frac1p2 = 1.42e-2;
  const double frac12  = -3.71e-3;
  const double fracsq  = -1.201e-2;
  const double fracz12 = -2.447e-2;
  const double fraczsq = -1.930e-2;


  INITSTATUS(status, "LALPSpinInspiralRDEngine", LALPSPININSPIRALRDENGINEC);
  ATTATCHSTATUSPTR(status);

  /* set parameters from InspiralTemplate structure*/

  /* Make sure parameter and waveform structures exist. */
  ASSERT(params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
  /*===================*/
  
  /* Compute some parameters*/
  
  if (paramsInit->nbins==0)
    {
      DETATCHSTATUSPTR(status);
      RETURN (status);
    }
  ak   = paramsInit->ak;

  mparams = &PSIRDparameters;

  m = params->totalMass * LAL_MTSUN_SI;
  unitHz = params->totalMass * LAL_MTSUN_SI * (REAL8)LAL_PI;
  /*    tSampling is in Hz, so dt is in seconds*/
  dt = 1.0/params->tSampling;

  /* -- length in seconds from Newtonian formula; */
  lengths = (5.0/256.0)/LAL_PI * pow(LAL_PI * params->chirpMass * LAL_MTSUN_SI * params->fLower,-5.0/3.0) / params->fLower;

  estlength = ceil(log10(lengths/dt)/log10(2.0));
  estlength = pow(2,estlength);
  estlength *= 2;

  /* The analytical formula for omega_match is in the do-while loop. However
     omegamatch can be controlled by fCutoff by un-commenting the following
     line and commenting the definition of omegamatch in the loop.*/
  //omegamatch = params->fCutoff *unitHz;
  if (params->eta>=0.25) sqrtOneMinus4Eta=0.;
  else                   sqrtOneMinus4Eta=sqrt(1.-4.*params->eta);
  omegamatch = omM0 +6.05e-3*sqrtOneMinus4Eta;

  while ((omegamatch * 16./unitHz) >  (REAL4)(subsampling)*params->tSampling ) subsampling*=2; 
  dt/= (REAL4)(subsampling);

  initPhi = params->startPhase;

  /* set initial values of dynamical variables*/
  initomega = params->fLower * unitHz;
  initomwrite=initomega;
  if (initomega>0.9*omegamatch) {
    initomega=0.9*omegamatch;
  } 
  initv = pow( initomega, oneby3 );  

  /* Here we use the following convention:
     the coordinates of the spin vectors params->spin1,2 are assumed to 
     be fed to this engine macro either in the frame set by the orbital angual 
     momentum (params->directionChoice= OrbitalL or TotalJ) or int the frame 
     set by the viewing direction (params->directionChoice=View).
     The gw spherical modes are coputed in wither of the three frames 
     specified by the params->directionChoice variable.
     The spin magnitude are nrmalized to the individual mass^2, i.e. 
     they are dimension-less.
     The modulus of the initial angular momentum is fixed by m1,m2 and 
     initial frequency, the inclination is the
     angle between the view direction and the orbital angular momentum.
     The polarization angle is not used here, it enters the pattern
     functions along with the angles marking the sky position of the
     source. */

  LNmag= params->eta * params->totalMass * params->totalMass / initv;

  mchi1=sqrt( params->spin1[0]*params->spin1[0] + params->spin1[1]*params->spin1[1] + params->spin1[2]*params->spin1[2] );
  mchi2=sqrt( params->spin2[0]*params->spin2[0] + params->spin2[1]*params->spin2[1] + params->spin2[2]*params->spin2[2] );

  /* Cosinus of the angle between the spin vectors */
  cosg  = (params->spin1[0]*params->spin2[0] + params->spin1[1]*params->spin2[1] + params->spin1[2]*params->spin2[2])/ mchi1/mchi2;

  switch (params->axisChoice) {
    
  case OrbitalL:
    cosa1 = params->spin1[2]/mchi1;
    cosa2 = params->spin2[2]/mchi2;
    initLNh[0]=0.;
    initLNh[1]=0.;
    initLNh[2]=1.;
    for (i=0;i<3;i++) {
      initS1[i]=params->spin1[i]*params->mass1*params->mass1;
      initS2[i]=params->spin2[i]*params->mass2*params->mass2;
    }
    inc=params->inclination; 
    break;
    
  case TotalJ:
    cosa1 = params->spin1[2]/mchi1;
    cosa2 = params->spin2[2]/mchi2;
    for (j=0;j<3;j++) {
      iS1[j] = params->spin1[j] * params->mass1 * params->mass1;
      iS2[j] = params->spin2[j] * params->mass2 * params->mass2;
      iJ[j] = iS1[j] + iS2[j];
    }
    iJ[2] += LNmag;
    iJmod = sqrt ( iJ[0]*iJ[0] + iJ[1]*iJ[1] + iJ[2]*iJ[2] );
    for (j=0;j<3;j++) {
      initJh[j] = iJ[j]/iJmod;
      initLNh[j]=0.;
      initS1[j]=0.;
      initS2[j]=0.;
    }
    if (initJh[0]==0.) phiJ=0.;
    else phiJ=atan2(initJh[1],initJh[0]);
    thetaJ=acos(initJh[2]);    
    rz[0][0]=cos(phiJ)  ; rz[0][1]=sin(phiJ); rz[0][2]=0.;
    rz[1][0]=-sin(phiJ) ; rz[1][1]=cos(phiJ); rz[1][2]=0.;
    rz[2][0]=0.         ; rz[2][1]=0.       ; rz[2][2]=1.;
    ry[0][0]=cos(thetaJ); ry[0][1]=0        ; ry[0][2]=-sin(thetaJ);
    ry[1][0]=0.         ; ry[1][1]=1.       ; ry[1][2]=0.;
    ry[2][0]=sin(thetaJ); ry[2][1]=0.       ; ry[2][2]=cos(thetaJ);
    for (j=0;j<3;j++) {
      for (k=0;k<3;k++) {
	initLNh[j] += ry[j][k] * rz[k][2];
	for (l=0;l<3;l++) {
	  initS1[j] += ry[j][k] * rz[k][l] * iS1[l];
	  initS2[j] += ry[j][k] * rz[k][l] * iS2[l];
	}
      } 
    }
    inc=params->inclination;
    break;
  default:
    //case View:
    cosa1 = (params->spin1[2]*cos(params->inclination)+params->spin1[0]*sin(params->inclination))/mchi1; 
    cosa2 = (params->spin1[2]*cos(params->inclination)+params->spin2[0]*sin(params->inclination))/mchi2; 
    for (i=0;i<3;i++) {
      initS1[i]=params->spin1[i]*params->mass1*params->mass1;
      initS2[i]=params->spin2[i]*params->mass2*params->mass2;
    }
    initLNh[0]=sin(params->inclination);
    initLNh[1]=0.;
    initLNh[2]=cos(params->inclination);
    inc=0.;
    break;
  }

  /*All the PN formulas used in this code assume that the spin variables 
    are the physical ones divided by totalmasss^2, here we introduce the 
    correct normalization, changing the input one.*/
  for (j=0;j<3;j++) {
    initS1[j] /= params->totalMass*params->totalMass;
    initS2[j] /= params->totalMass*params->totalMass;
  }

  if (signalvec1) {
    length=signalvec1->length; 
  }
  else {
    if (ff) length = ff->length;
    else length=0;
  }

  dummy.length = nn * 6;
  
  values.length = dvalues.length = newvalues.length = yt.length = dym.length = dyt.length = nn;
  
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
  mparams->m2m1 = params->mass2/params->mass1;
  mparams->m1m2 = params->mass1/params->mass2;
  mparams->m1m  = params->mass1/params->totalMass;
  mparams->m2m  = params->mass2/params->totalMass;
  mparams->dm   = (params->mass1-params->mass2)/params->totalMass;

  /* params->eta might have been set up before but just for safety, we
   * recompute it here below.*/
  params->eta =  (params->mass1 * params->mass2)/
    (params->mass1 + params->mass2)/(params->mass1 + params->mass2);
  mparams->eta = params->eta;

  for (j = LAL_PNORDER_NEWTONIAN; j <= params->order; j++) {
    mparams->wdotorb[j] = ak.ST[j];
  }
  mparams->wdotorblog=0.;
  for (j = params->order +1; j < 8; j++){
    mparams->wdotorb[j] = 0.;
  }
  if ((params->order)>=6) {
    mparams->wdotorblog=ak.ST[7];
    if ((params->order)==7) mparams->wdotorb[7]=ak.ST[8];
  }

  mparams->epnorb[0]=ak.ETaN;

  switch (params->order){
  case     LAL_PNORDER_NEWTONIAN:
  case     LAL_PNORDER_HALF:
    break;
  case     LAL_PNORDER_ONE:
    mparams->epnorb[1]     = ak.ETa1;
    break;
  case     LAL_PNORDER_ONE_POINT_FIVE:

    mparams->epnorb[1]     = ak.ETa1;
    mparams->epnspin15S1dotLNh = 8./3.+2.*mparams->m2m1;
    mparams->epnspin15S2dotLNh = 8./3.+2.*mparams->m1m2;

    mparams->wdotspin15S1LNh  = -( 113.0 + 75.0 * mparams->m2m1 ) / 12.0;
    mparams->wdotspin15S2LNh  = -( 113.0 + 75.0 * mparams->m1m2 ) / 12.0;

    mparams->LNhdot15      = 0.5;

    mparams->S1dot15 	   = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
    mparams->S2dot15 	   = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;

    mparams->LNhdot15 	= 0.5;
    mparams->S1dot15 	= (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
    mparams->S2dot15 	= (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;

    break;
  case     LAL_PNORDER_TWO:

    mparams->epnorb[1]         = ak.ETa1;
    mparams->epnspin15S1dotLNh = 8./3.+2.*mparams->m2m1;
    mparams->epnspin15S2dotLNh = 8./3.+2.*mparams->m1m2;
    mparams->epnorb[2]         = ak.ETa2;

    mparams->wdotspin15S1LNh   = -( 113.0 + 75.0 * mparams->m2m1 ) / 12.0;
    mparams->wdotspin15S2LNh   = -( 113.0 + 75.0 * mparams->m1m2 ) / 12.0;
    mparams->wdotspin20S1S2    = -(1.0/48.0) / mparams->eta ;
    mparams->wdotspin20S1S1    = 1./96.;

    mparams->LNhdot15 	       = 0.5;
    mparams->LNhdot20          = -1.5 / mparams->eta;

    mparams->S1dot15 	       = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
    mparams->S2dot15 	       = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;
    mparams->Sdot20 	       = 0.5 ;
    break;
  case     LAL_PNORDER_TWO_POINT_FIVE:

    mparams->epnorb[1]         = ak.ETa1;
    mparams->epnspin15S1dotLNh = 8./3.+2.*mparams->m2m1;
    mparams->epnspin15S2dotLNh = 8./3.+2.*mparams->m1m2;
    mparams->epnorb[2]         = ak.ETa2;
    mparams->epnspin25S1dotLNh = 8. - 31./9.*mparams->eta + (3.-10./3.*mparams->eta)*mparams->m2m1;
    mparams->epnspin25S2dotLNh = 8. - 31./9.*mparams->eta + (3.-10./3.*mparams->eta)*mparams->m1m2;

    mparams->wdotspin15S1LNh   = -( 113.0 + 75.0 * mparams->m2m1 ) / 12.0;
    mparams->wdotspin15S2LNh   = -( 113.0 + 75.0 * mparams->m1m2 ) / 12.0;
    mparams->wdotspin20S1S2    = -(1.0/48.0) / mparams->eta ;
    mparams->wdotspin20S1S1    = 1./96.;
    mparams->wdotspin25S1LNh   = -31319./1008. + 1159./24.* mparams->eta + (-809./84.+281./8.*mparams->eta)*mparams->m2m1;
    mparams->wdotspin25S2LNh   = -31319./1008. + 1159./24.* mparams->eta + (-809./84.+281./8.*mparams->eta)*mparams->m1m2;

    mparams->LNhdot15 	       = 0.5;
    mparams->LNhdot20 	       = -1.5 / mparams->eta;

    mparams->S1dot15 	       = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
    mparams->S2dot15 	       = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;
    mparams->Sdot20            = 0.5;
    mparams->S1dot25 	       = 0.5625 + 1.25*mparams->eta - mparams->eta*mparams->eta/24. + mparams->dm*(-0.5625+0.625*mparams->eta);
    mparams->S2dot25 	       = 0.5625 + 1.25*mparams->eta - mparams->eta*mparams->eta/24. - mparams->dm*(-0.5625+0.625*mparams->eta);
    break;
  case     LAL_PNORDER_THREE:
  case     LAL_PNORDER_THREE_POINT_FIVE:

    mparams->epnorb[1]           = ak.ETa1;
    mparams->epnspin15S1dotLNh   = 8./3.+2.*mparams->m2m1;
    mparams->epnspin15S2dotLNh   = 8./3.+2.*mparams->m1m2;
    mparams->epnorb[2]           = ak.ETa2;
    mparams->epnspin20S1S2       = 1./mparams->eta;
    mparams->epnspin20S1S2dotLNh = -3./mparams->eta;
    mparams->epnspin20S1S1       = (1.+mparams->m2m1)*(1.+mparams->m2m1)/2.;
    mparams->epnspin20S2S2       = (1.+mparams->m1m2)*(1.+mparams->m1m2)/2.;
    mparams->epnspin20S1S1dotLNh = -3.*(1.+mparams->m2m1)*(1.+mparams->m2m1)/2.;
    mparams->epnspin20S2S2dotLNh = -3.*(1.+mparams->m1m2)*(1.+mparams->m1m2)/2.;
    mparams->epnspin25S1dotLNh   = 8. - 31./9.*mparams->eta + (3.-10./3.*mparams->eta)*mparams->m2m1;
    mparams->epnspin25S2dotLNh   = 8. - 31./9.*mparams->eta + (3.-10./3.*mparams->eta)*mparams->m1m2;
    mparams->epnorb[3]           = ak.ETa3;

    mparams->wdotspin15S1LNh     = -( 113.0 + 75.0 * mparams->m2m1 ) / 12.0;
    mparams->wdotspin15S2LNh     = -( 113.0 + 75.0 * mparams->m1m2 ) / 12.0;
    mparams->wdotspin20S1S2      = -(1.0/48.0) / mparams->eta ;
    mparams->wdotspin20S1S1      = 1./96.;
    mparams->wdotspin25S1LNh   = -31319./1008. + 1159./24.* mparams->eta + (-809./84.+281./8.*mparams->eta)*mparams->m2m1;
    mparams->wdotspin25S2LNh   = -31319./1008. + 1159./24.* mparams->eta + (-809./84.+281./8.*mparams->eta)*mparams->m1m2;
    mparams->S1dot15 	         = (4.0 + 3.0 * mparams->m2m1) / 2.0 * mparams->eta;
    mparams->S2dot15 	         = (4.0 + 3.0 * mparams->m1m2) / 2.0 * mparams->eta;
    mparams->Sdot20              = 0.5;

    mparams->S1dot25 	       = 0.5625 + 1.25*mparams->eta - mparams->eta*mparams->eta/24. + mparams->dm*(-0.5625+0.625*mparams->eta);
    mparams->S2dot25 	       = 0.5625 + 1.25*mparams->eta - mparams->eta*mparams->eta/24. - mparams->dm*(-0.5625+0.625*mparams->eta);
    break;
  case     LAL_PNORDER_PSEUDO_FOUR:
    fprintf(stderr,"*** LALPhenSpinInspiralRD ERROR: PhenSpin approximant not available at pseudo4PN order\n");
    break;
  case     LAL_PNORDER_NUM_ORDER:
    fprintf(stderr,"*** LALPhenSpinInspiralRD ERROR: NUM_ORDER not a valid PN order\n");
  }

  /* setup initial conditions for dynamical variables*/

  Phi  = initPhi;
  omega = initomega;

  LNhx = initLNh[0];
  LNhy = initLNh[1];
  LNhz = initLNh[2];

  S1x = initS1[0];
  S1y = initS1[1];
  S1z = initS1[2];

  S2x = initS2[0];
  S2y = initS2[1];
  S2z = initS2[2];

  /* copy everything in the "values" structure*/

  values.data[0] = Phi;
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

  in4.function 	= LALPSpinInspiralRDderivatives;
  in4.y 	= &values;
  in4.dydx 	= &dvalues;
  in4.h 	= dt/m;
  in4.n 	= nn;
  in4.yt 	= &yt;
  in4.dym 	= &dym;
  in4.dyt 	= &dyt;

  /* Allocate memory for temporary arrays */
  h2P2  = XLALCreateREAL4Vector ( length*2 );
  h2M2  = XLALCreateREAL4Vector ( length*2 );
  h2P1  = XLALCreateREAL4Vector ( length*2 );
  h2M1  = XLALCreateREAL4Vector ( length*2 );
  h20   = XLALCreateREAL4Vector ( length*2 );
  h3P3  = XLALCreateREAL4Vector ( length*2 );
  h3M3  = XLALCreateREAL4Vector ( length*2 );
  sig1 = XLALCreateREAL4Vector ( length);
  sig2 = XLALCreateREAL4Vector ( length);
  hap  = XLALCreateREAL4Vector ( length*2);
  fap  = XLALCreateREAL4Vector ( length );
  phap = XLALCreateREAL8Vector ( length );
  shift22 = XLALCreateREAL4Vector ( length );
  
  if ( !h2P2 || !h2M2 || !h2P1 || !h2M1 || !h20 || !sig1 || !sig2 || !fap || !phap || !shift22 || !hap || !h3P3 || !h3M3 )
    {
      if ( h2P2 ) XLALDestroyREAL4Vector( h2P2 );
      if ( h2M2 ) XLALDestroyREAL4Vector( h2M2 );
      if ( h2P1 ) XLALDestroyREAL4Vector( h2P1 );
      if ( h2M2 ) XLALDestroyREAL4Vector( h2M1 );
      if ( h20 ) XLALDestroyREAL4Vector( h20 );
      if ( h3P3 ) XLALDestroyREAL4Vector( h3P3 );
      if ( h3M3 ) XLALDestroyREAL4Vector( h3M3 );
      if ( sig1 ) XLALDestroyREAL4Vector( sig1 );
      if ( sig2 ) XLALDestroyREAL4Vector( sig2 );
      if ( fap ) XLALDestroyREAL4Vector( fap );
      if ( hap ) XLALDestroyREAL4Vector( hap );
      if ( phap ) XLALDestroyREAL8Vector( phap );
      if ( shift22 ) XLALDestroyREAL4Vector( shift22 );
      LALFree( dummy.data );
      ABORT( status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM );
    }
  
  memset(h2P2->data, 0, h2P2->length * sizeof( REAL4 ));
  memset(h2M2->data, 0, h2M2->length * sizeof( REAL4 ));
  memset(h2P1->data, 0, h2P1->length * sizeof( REAL4 ));
  memset(h2M1->data, 0, h2P1->length * sizeof( REAL4 ));
  memset(h20->data, 0, h20->length * sizeof( REAL4 ));
  memset(h3P3->data, 0, h3P3->length * sizeof( REAL4 ));
  memset(h3M3->data, 0, h3M3->length * sizeof( REAL4 ));
  memset(sig1->data, 0, sig1->length * sizeof( REAL4 ));
  memset(sig2->data, 0, sig2->length * sizeof( REAL4 ));
  memset(hap->data, 0, hap->length * sizeof( REAL4 ));
  memset(fap->data, 0, fap->length * sizeof( REAL4 ));
  memset(phap->data, 0, phap->length * sizeof( REAL8 ));
  memset(shift22->data, 0, shift22->length * sizeof( REAL4 ));
  
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
  write = 0;

  alpha=atan2(LNhy,LNhx);

  LALPSpinInspiralRDderivatives(&values,&dvalues,(void*)mparams);

  /* Injection: hh,ff; template: signalvec1,2*/

  if ( hh || signalvec2)
    params->nStartPad = 0; /* must be zero for templates and injection */  

  /* The number of Ring Down modes is hard-coded here, it cannot exceed 3*/
  nmodes = 2;

  /* For RD, check that the 220 QNM freq. is less than the Nyquist freq. */
  /* Get QNM frequencies */
  modefreqs = XLALCreateCOMPLEX8Vector( nmodes );

  /* Call XLALFinalMassSpin() to get mass and spin of the final black hole */
  errcode = XLALPSpinFinalMassSpin(&finalMass, &finalSpin, params, energy, initLNh);
  if ( errcode != XLAL_SUCCESS )
  {
    XLALDestroyCOMPLEX8Vector( modefreqs );
    ABORTXLAL( status );
  }

  xlalSt2P2 = XLALPSpinGenerateQNMFreq( modefreqs, params, 2, 2, nmodes, finalMass, finalSpin );
  if ( xlalSt2P2 != XLAL_SUCCESS )
    {
      XLALDestroyCOMPLEX8Vector( modefreqs );
      ABORTXLAL( status );
    }  

  omegaRD=modefreqs->data[0].re * unitHz / LAL_PI / 2.;    
  /* If Nyquist freq. <  220 QNM freq., print warning message, but go on*/
  /* Note that we cancelled a factor of 2 occuring on both sides */
  if ( params->tSampling < modefreqs->data[0].re / LAL_PI ) {
    fprintf(stdout, "LALPhenSpin WARNING : Estimated ringdown freq larger than Nyquist freq.\n" );
    /* XLALDestroyCOMPLEX8Vector( modefreqs ); 
    LALError(status->statusPtr, message);
    ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE); */
  }

  v2=pow(omega,2./3.);
  v=sqrt(v2);


  params->ampOrder=0;

  if (params->distance > 0.) amp22ini= -2.0 * params->mu * LAL_MRSUN_SI/(params->distance) * sqrt( 16.*LAL_PI/5.);
  else amp22ini  = 2. * sqrt( LAL_PI / 5.0) * params->signalAmplitude;

  do {

    if ( write >= length) {
      XLALRungeKutta4Free( integrator );
      LALFree(dummy.data);
      ABORT(status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
    }

    amp22 = amp22ini * v2;
    amp20 = amp22 * sqrt(1.5);
    amp33 = -amp22/4. * sqrt(5./42.);

    omegaold    = omega;
    omegadotold = omegadot;
    enold       = energy;

    Psi = Phi - 2. * omega * log(omega);

    ci=(LNhz);
    c2i=ci*ci;
    s2i=1.-ci*ci;
    si=sqrt(s2i);
    c2i2 = (1. + ci)/2.;
    s2i2 = (1. - ci)/2.;
    ci2  = sqrt(c2i2);
    si2  = sqrt(s2i2);
    c4i2 = c2i2*c2i2;
    s4i2 = s2i2*s2i2;
    c5i2 = c4i2*ci2;
    s5i2 = s4i2*si2;
    c6i2 = c4i2*c2i2;
    s6i2 = s4i2*s2i2;

    if (omega>initomwrite) {

      if (count%subsampling==0) {

	// amp22= -2.0 * params->mu * LAL_MRSUN_SI/(params->distance) * sqrt( 16.*LAL_PI/5.)*v2;
	// amp20 = amp22*sqrt(3/2)
	// Y22 \pm Y2-2= sqrt(5/PI)    ((1+cos^2 t)/4, (cos t)/2) 
	// Y21 \pm Y2-1= sqrt(5/PI)    ((sin t)/2, (sin 2t)/4)
	// Y20         = sqrt(15/2 PI) (sin^2 t)/4

	h2P2->data[2*write]   = (REAL4)(amp22 * ( 1./(1.+v2/42.*(107.-55.*mparams->eta)) * ( cos(2.*(Psi+alpha))*c4i2 + cos(2.*(Psi-alpha))*s4i2 ) + v * mparams->dm/3. * si * ( cos(Psi-2.*alpha)*s2i2 + cos(Psi+2.*alpha)*c2i2 ) ) );
	h2P2->data[2*write+1] = (REAL4)(amp22 * ( 1./(1.+v2/42.*(107.-55.*mparams->eta)) * ( sin(2.*(Psi+alpha))*c4i2 - sin(2.*(Psi-alpha))*s4i2 ) + v * mparams->dm/3. * si * (-sin(Psi-2.*alpha)*s2i2 + sin(Psi+2.*alpha)*c2i2 ) ) );

	h2M2->data[2*write]   = (REAL4)(amp22 * ( 1./(1.+v2/42.*(107.-55.*mparams->eta)) * ( cos(2.*(Psi+alpha))*c4i2 + cos(2.*(Psi-alpha))*s4i2 ) - v * mparams->dm/3. * si * ( cos(Psi-2.*alpha)*s2i2 + cos(Psi+2.*alpha)*c2i2 ) ) );
	h2M2->data[2*write+1] = (REAL4)(amp22 * ( 1./(1.+v2/42.*(107.-55.*mparams->eta)) * ( sin(2.*(Psi+alpha))*c4i2 - sin(2.*(Psi-alpha))*s4i2 ) - v * mparams->dm/3. * si * (-sin(Psi-2.*alpha)*s2i2 + sin(Psi+2.*alpha)*c2i2 ) ) );
	
	h2P1->data[2*write]   = (REAL4)(amp22 * (si / (1.+v2/84.*(107.-55.*mparams->eta)) * (-cos(2.*Psi-alpha)*s2i2 + cos(2.*Psi+alpha)*c2i2 ) + v * mparams->dm/3. * ( -cos(Psi+alpha)*(ci+c2i)/2. - cos(Psi-alpha)*s2i2*(1.+2.*ci) ) ) );
	h2P1->data[2*write+1] = (REAL4)(amp22 * (si / (1.+v2/84.*(107.-55.*mparams->eta)) * ( sin(2.*Psi-alpha)*s2i2 + sin(2.*Psi+alpha)*c2i2 ) + v * mparams->dm/3. * (-sin(Psi+alpha)*(ci+c2i)/2. + sin(Psi-alpha)*s2i2*(1.+2.*ci) ) ) ); 

	h2M1->data[2*write]   = (REAL4)(amp22 * (si / (1.+v2/84.*(107.-55.*mparams->eta)) * (-cos(2.*Psi-alpha)*s2i2 + cos(2.*Psi+alpha)*c2i2 ) - v * mparams->dm/3. * ( -cos(Psi+alpha)*(ci+c2i)/2. - cos(Psi-alpha)*s2i2*(1.+2.*ci) ) ) );
	h2M1->data[2*write+1] = (REAL4)(amp22 * (si / (1.+v2/84.*(107.-55.*mparams->eta)) * ( sin(2.*Psi-alpha)*s2i2 + sin(2.*Psi+alpha)*c2i2 ) -  v * mparams->dm/3. * (-sin(Psi+alpha)*(ci+c2i)/2. + sin(Psi-alpha)*s2i2*(1.+2.*ci) ) ) ); 

	h20->data[2*write]   = (REAL4)(amp20 * ( s2i * cos( 2.* Psi) ) );
	h20->data[2*write+1] = (REAL4)(amp20 * ( - v * mparams->dm* 8/9. * si * sin(Psi) ) );
	
	h3P3->data[2*write]   = (REAL4)(amp33 * ( v * mparams->dm * (-9.*cos(3.*(Psi-alpha))*s6i2 - cos(Psi-3.*alpha)*s4i2*c2i2 + cos(Psi+3.*alpha)*s2i2*c4i2 + 9.*cos(3.*(Psi+alpha))*c6i2 ) +  v2 * 8.*(1.-3.*mparams->eta) * (-cos(2.*Psi+3.*alpha)*s5i2*ci2 + cos(2.*Psi-3.*alpha)*c5i2*si2 ) ) );
	h3P3->data[2*write+1] = (REAL4)(amp33 * ( v * mparams->dm * ( 9.*sin(3.*(Psi-alpha))*s6i2 + sin(Psi-3.*alpha)*s4i2*c2i2 + sin(Psi+3.*alpha)*s2i2*c4i2 + 9.*sin(3.*(Psi+alpha))*c6i2 ) + v2 * 8.*(1.-3.*mparams->eta) * (-sin(2.*Psi+3.*alpha)*s5i2*ci2 + sin(2.*Psi-3.*alpha)*c5i2*si2 ) ) );

	h3M3->data[2*write]   = (REAL4)(amp33 * (-v * mparams->dm * (-9.*cos(3.*(Psi-alpha))*s6i2 - cos(Psi-3.*alpha)*s4i2*c2i2 + cos(Psi+3.*alpha)*s2i2*c4i2 + 9.*cos(3.*(Psi+alpha))*c6i2 ) +  v2 * 8.*(1.-3.*mparams->eta) * (-cos(2.*Psi+3.*alpha)*s5i2*ci2 + cos(2.*Psi-3.*alpha)*c5i2*si2 ) ) );
	h3M3->data[2*write+1] = (REAL4)(amp33 * (-v * mparams->dm * ( 9.*sin(3.*(Psi-alpha))*s6i2 + sin(Psi-3.*alpha)*s4i2*c2i2 + sin(Psi+3.*alpha)*s2i2*c4i2 + 9.*sin(3.*(Psi+alpha))*c6i2 ) + v2 * 8.*(1.-3.*mparams->eta) * (-sin(2.*Psi+3.*alpha)*s5i2*ci2 + sin(2.*Psi-3.*alpha)*c5i2*si2 ) ) );

	fap->data[write] =  (REAL4)( omega );
	
	phap->data[write] =  (REAL8)( Psi );
	
	shift22->data[write] = alpha;
	
	write++;
	
      }

    }

    in4.x = t/m;

    LALRungeKutta4(status->statusPtr, &newvalues, integrator, (void*)mparams);
    CHECKSTATUSPTR(status);

    /* updating values of dynamical variables*/

    Phi  = values.data[0] = newvalues.data[0];

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

    alpha = atan2(LNhy,LNhx);
    
    S1dotL=(S1x*LNhx+S1y*LNhy+S1z*LNhz)*params->totalMass*params->totalMass/params->mass1/params->mass1;
    S2dotL=(S2x*LNhx+S2y*LNhy+S2z*LNhz)*params->totalMass*params->totalMass/params->mass2/params->mass2;
    S1dotS1=(S1x*S1x+S1y*S1y+S1z*S1z)*params->totalMass*params->totalMass*params->totalMass*params->totalMass/params->mass1/params->mass1/params->mass1/params->mass1;
    S2dotS2=(S2x*S2x+S2y*S2y+S2z*S2z)*params->totalMass*params->totalMass*params->totalMass*params->totalMass/params->mass2/params->mass2/params->mass2/params->mass2;
    S1dotS2=(S1x*S2x+S1y*S2y+S1z*S2z)*params->totalMass*params->totalMass*params->totalMass*params->totalMass/params->mass1/params->mass1/params->mass2/params->mass2;

    LALPSpinInspiralRDderivatives(&values,&dvalues,(void*)mparams);

    dLNhx=dvalues.data[2];
    dLNhy=dvalues.data[3];
    dLNhz=dvalues.data[4];

    LNhxy=sqrt(LNhx*LNhx+LNhy*LNhy);
    alphadotold=alphadot;
    iotadotold=iotadot;
    if (LNhxy>0.) {
      iotadot  = - dLNhz / LNhxy;
      alphadot = (LNhx*dLNhy - LNhy*dLNhx) / LNhxy;
      iotaddot =(iotadot-iotadotold)/dt*m;
      alphaddot=(alphadot-alphadotold)/dt*m;
    }
    else {
      iotadot=0.;
      alphadot=0.;
      iotaddot=0.;
      alphaddot=0.;
    }

    energy = values.data[11];
    v2=pow(omega,2./3.);
    v=sqrt(v2);

    omegadot = dvalues.data[1];
    omegaddot= (omegadot-omegadotold)/dt*m;

    t = (++count - params->nStartPad) * dt;

    //adjourn ommatch

    omegamatch= omM0 + 6.05e-3*sqrtOneMinus4Eta + omMz1p2*(S1dotL+S2dotL) + omM12*(S1dotS2-S1dotL*S2dotL) + omMsq*(S1dotS1+S2dotS2-S1dotL*S1dotL-S2dotL*S2dotL) + omMz12*(S1dotL*S2dotL) +omMzsq*(S1dotL*S1dotL+S2dotL*S2dotL);
  
  }

 /* Test that omega/unitHz < NYQUIST */
  
  while ( (energy <= 0.99*enold) && (omega >= 0.99*omegaold) && (omega/unitHz < params->tSampling) && (!(isnan(omega))) && (omega < omegamatch) );

 /* if code stopped since evolving quantities became nan write an error message */
  if ( omega/unitHz >= params->tSampling ) {
    fprintf(stderr,
            "** LALPhenSpinInspiralRD ERROR **: f=%11.4e  rate=%11.4e "
	    "m1: %e, "
	    "m2: %e, "
	    "spin1x: %e, "
	    "spin1y: %e, "
	    "spin1z: %e, "
	    "spin2x: %e, "
	    "spin2y: %e, "
	    "spin2z: %e, "
	    "inclination: %e\n ",
	    omega/unitHz,params->tSampling,
	    params->mass1, params->mass2,
	    params->spin1[0], params->spin1[1], params->spin1[2],
	    params->spin2[0], params->spin2[1], params->spin2[2],
	    params->inclination);
    ABORTXLAL( status );
  }
  else if (isnan(omega)){
    fprintf(stderr,
	    "** LALPhenSpinInspiralRD ERROR **: omega has become nan. "
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
    ABORTXLAL( status );
  }
  /* if code stopped due to co-ord singularity write an error message */
  else if ((LNhx*LNhx + LNhy*LNhy + LNhz*LNhz) < (1.0 - LNhztol)){
    fprintf( stderr,
	     "** LALPhenSpinInspiralRD ERROR **: waveform terminated, coord singularity. "
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
    ABORTXLAL(status);
  }
  else if ( ( energy > enold) && (omega<omegamatch) ) {
    fprintf(stderr, "** LALPhenSpinInspiralRD ERROR **: waveform terminated, energy increases %11.3e > %11.3e  int. step=%d  omega=%11.3e\n",energy,enold,count,omega);
    ABORTXLAL(status);
  }
  else if ( omega < omegaold) {
    fprintf( stderr,
	     "** LALPhenSpinInspiralRD ERROR **: waveform terminated, omega decrease "
	     "om: %e,  omold: %e"
	     "m1: %e, "
	     "m2: %e, "
	     "spin1x: %e, "
	     "spin1y: %e, "
	     "spin1z: %e, "
	     "spin2x: %e, "
	     "spin2y: %e, "
	     "spin2z: %e, "
	     "inclination: %e\n ",
	     omega,omegaold,
	     params->mass1, params->mass2,
	     params->spin1[0], params->spin1[1], params->spin1[2],
	     params->spin2[0], params->spin2[1], params->spin2[2],
	     params->inclination);
  }
  else if ( omega >= omegamatch) {
    rett=1;
  }

  t0=t-dt;
  tAs=t0+2.*omegadot/omegaddot*m;
  om1=omegadot*tAs/m*(1.-t0/tAs)*(1.-t0/tAs);
  om0=omega-om1/(1.-t0/tAs);
  
  iotadot1=iotaddot*tAs/m*(1.-t0/tAs)*(1.-t0/tAs);
  iotadot0=iotadot-iotadot1/(1.-t0/tAs);
  
  alphadot1=alphaddot*tAs/m*(1.-t0/tAs)*(1.-t0/tAs);
  alphadot0=alphadot-alphadot1/(1.-t0/tAs);

  if ((tAs<0)||(om1<0.)) {
    fprintf(stderr,"** LALPSpinInspiralRD ERROR **: Could not attach phen part\n");
    rett=0;
    ABORTXLAL(status);
  }

  dt=1./params->tSampling;
  Psi0=Psi+tAs*om1*log(1.-t0/tAs)/m;
  alpha0=alpha+tAs*alphadot1*log(1.-t0/tAs)/m;
  iota0=acos(LNhz)+iotadot1*tAs/m*log(1.-t0/tAs);

  /* Get QNM frequencies */
  modefreqs = XLALCreateCOMPLEX8Vector( nmodes );

  errcode = XLALPSpinFinalMassSpin(&finalMass, &finalSpin, params, energy, initLNh);
  if ( errcode != XLAL_SUCCESS )
    {
      ABORTXLAL( status );
    }

  xlalSt2P2 = XLALPSpinGenerateQNMFreq( modefreqs, params, 2, 2, nmodes, finalMass, finalSpin );
  if ( xlalSt2P2 != XLAL_SUCCESS )
    {
      XLALDestroyCOMPLEX8Vector( modefreqs );
      ABORTXLAL( status );
    }
  
  omegaRD=modefreqs->data[0].re * unitHz / LAL_PI /2.;

  fracRD = frac0 + frac1p2*(S1dotL+S2dotL) + frac12*(S1dotS2-S1dotL*S2dotL) + fracsq*(S1dotS1+S2dotS2-S1dotL*S1dotL-S2dotL*S2dotL) + fracz12*(S1dotL*S2dotL) + fraczsq*(S1dotL*S1dotL+S2dotL*S2dotL);

  /* Now the phenomenological part is added */
  if (rett==1) {

    count=write;

    do {
      
      if ( count >= length) {
	fprintf(stderr,"** LALPhenSpinInspiralRD ERROR**: phen part exceeds array length");
	ABORT(status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
      }

      omegaold=omega;
      omega = om1/(1.-t/tAs)+om0;
      fap->data[count] =  omega;
      Psi = Psi0-tAs*om1*log(1.-t/tAs)/m+om0*(t-t0)/m;      
      iota=iotadot0*(t-t0)/m-iotadot1*tAs/m*log(1.-t/tAs)+iota0;
      alpha=alpha0+alphadot0*(t-t0)/m-alphadot1*tAs/m*log(1.-t/tAs);
      v2old=v2;
      v2=pow(omega,2./3.);
      v=sqrt(v2);
      amp22*=v2/v2old;
      amp20=amp22*sqrt(1.5);
      amp33=-amp22/4.*sqrt(5./42.);

      ci=cos(iota);
      c2i2=(1.+ci)/2.;
      s2i2=1.-c2i2;
      ci2=sqrt(c2i2);
      si=sqrt(1.-ci*ci);
      s2i=1.-ci*ci;
      c4i2=c2i2*c2i2;
      s4i2=s2i2*s2i2;
      c5i2=c4i2*ci2;
      s5i2=s4i2*si2;
      c6i2=c4i2*c2i2;
      s6i2=s4i2*s2i2;

      h2P2->data[2*count]   = (REAL4)(amp22 * ( 1./(1.+v2/42.*(107.-55.*mparams->eta)) * ( cos(2.*(Psi+alpha))*c4i2 + cos(2.*(Psi-alpha))*s4i2 ) + v * mparams->dm/3. * si * ( cos(Psi-2.*alpha)*s2i2 + cos(Psi+2.*alpha)*c2i2 ) ) );
      h2P2->data[2*count+1] = (REAL4)(amp22 * ( 1./(1.+v2/42.*(107.-55.*mparams->eta)) * ( sin(2.*(Psi+alpha))*c4i2 - sin(2.*(Psi-alpha))*s4i2 ) + v * mparams->dm/3. * si * (-sin(Psi-2.*alpha)*s2i2 + sin(Psi+2.*alpha)*c2i2 ) ) );
      
      h2M2->data[2*count]   = (REAL4)(amp22 * ( 1./(1.+v2/42.*(107.-55.*mparams->eta)) * ( cos(2.*(Psi+alpha))*c4i2 + cos(2.*(Psi-alpha))*s4i2 ) - v * mparams->dm/3. * si * ( cos(Psi-2.*alpha)*s2i2 + cos(Psi+2.*alpha)*c2i2 ) ) );
      h2M2->data[2*count+1] = (REAL4)(amp22 * ( 1./(1.+v2/42.*(107.-55.*mparams->eta)) * ( sin(2.*(Psi+alpha))*c4i2 - sin(2.*(Psi-alpha))*s4i2 ) - v * mparams->dm/3. * si * (-sin(Psi-2.*alpha)*s2i2 + sin(Psi+2.*alpha)*c2i2 ) ) );
      
      h2P1->data[2*count]   = (REAL4)(amp22 * (si / (1.+v2/84.*(107.-55.*mparams->eta)) * (-cos(2.*Psi-alpha)*s2i2 + cos(2.*Psi+alpha)*c2i2 ) + v * mparams->dm/3. * ( -cos(Psi+alpha)*(ci+c2i)/2. - cos(Psi-alpha)*s2i2*(1.+2.*ci) ) ) );
      h2P1->data[2*count+1] = (REAL4)(amp22 * (si / (1.+v2/84.*(107.-55.*mparams->eta)) * ( sin(2.*Psi-alpha)*s2i2 + sin(2.*Psi+alpha)*c2i2 ) + v * mparams->dm/3. * (-sin(Psi+alpha)*(ci+c2i)/2. + sin(Psi-alpha)*s2i2*(1.+2.*ci) ) ) ); 
      
      h2M1->data[2*count]   = (REAL4)(amp22 * (si / (1.+v2/84.*(107.-55.*mparams->eta)) * (-cos(2.*Psi-alpha)*s2i2 + cos(2.*Psi+alpha)*c2i2 ) - v * mparams->dm/3. * ( -cos(Psi+alpha)*(ci+c2i)/2. - cos(Psi-alpha)*s2i2*(1.+2.*ci) ) ) );
      h2M1->data[2*count+1] = (REAL4)(amp22 * (si / (1.+v2/84.*(107.-55.*mparams->eta)) * ( sin(2.*Psi-alpha)*s2i2 + sin(2.*Psi+alpha)*c2i2 ) - v * mparams->dm/3. * (-sin(Psi+alpha)*(ci+c2i)/2. + sin(Psi-alpha)*s2i2*(1.+2.*ci) ) ) ); 
      
      h20->data[2*count]   = (REAL4)(amp20 * ( s2i * cos( 2.* Psi) ) );
      h20->data[2*count+1] = (REAL4)(amp20 * ( - v * mparams->dm* 8/9. * si * sin(Psi) ) );
      
      h3P3->data[2*count]   = (REAL4)(amp33 * ( v * mparams->dm * (-9.*cos(3.*(Psi-alpha))*s6i2 - cos(Psi-3.*alpha)*s4i2*c2i2 + cos(Psi+3.*alpha)*s2i2*c4i2 + 9.*cos(3.*(Psi+alpha))*c6i2 ) +  v2 * 8.*(1.-3.*mparams->eta) * (-cos(2.*Psi+3.*alpha)*s5i2*ci2 + cos(2.*Psi-3.*alpha)*c5i2*si2 ) ) );
      h3P3->data[2*count+1] = (REAL4)(amp33 * ( v * mparams->dm * ( 9.*sin(3.*(Psi-alpha))*s6i2 + sin(Psi-3.*alpha)*s4i2*c2i2 + sin(Psi+3.*alpha)*s2i2*c4i2 + 9.*sin(3.*(Psi+alpha))*c6i2 ) + v2 * 8.*(1.-3.*mparams->eta) * (-sin(2.*Psi+3.*alpha)*s5i2*ci2 + sin(2.*Psi-3.*alpha)*c5i2*si2 ) ) );

      h3M3->data[2*count]   = (REAL4)(amp33 * (-v * mparams->dm * (-9.*cos(3.*(Psi-alpha))*s6i2 - cos(Psi-3.*alpha)*s4i2*c2i2 + cos(Psi+3.*alpha)*s2i2*c4i2 + 9.*cos(3.*(Psi+alpha))*c6i2 ) +  v2 * 8.*(1.-3.*mparams->eta) * (-cos(2.*Psi+3.*alpha)*s5i2*ci2 + cos(2.*Psi-3.*alpha)*c5i2*si2 ) ) );
      h3M3->data[2*count+1] = (REAL4)(amp33 * (-v * mparams->dm * ( 9.*sin(3.*(Psi-alpha))*s6i2 + sin(Psi-3.*alpha)*s4i2*c2i2 + sin(Psi+3.*alpha)*s2i2*c4i2 + 9.*sin(3.*(Psi+alpha))*c6i2 ) + v2 * 8.*(1.-3.*mparams->eta) * (-sin(2.*Psi+3.*alpha)*s5i2*ci2 + sin(2.*Psi-3.*alpha)*c5i2*si2 ) ) );
      
      fap->data[count] =  (REAL4)( omega );      
      
      phap->data[count] =  (REAL8)( Psi );
      
      shift22->data[count] = alpha;
      
      count++;
      t+=dt;    

      /*aggiungere la formula per frac*/


    } while ((omega < fracRD*omegaRD));
    
  }
  else {
    fprintf(stderr,"** LALPhenSpinInspiralRD ERROR **: No phen part added\n");
    fprintf(stderr,"** m (%11.4e  %11.4e)  f0 %11.4e\n",params->mass1,params->mass2,params->fLower);
    fprintf(stderr,"** S1 (%8.4f  %8.4f  %8.4f)\n",initS1[0],initS1[1],initS1[2]);
    fprintf(stderr,"** S2 (%8.4f  %8.4f  %8.4f)\n",initS2[0],initS2[1],initS2[2]);
    fprintf(stderr,"** om (%11.4e  %11.4e)  E (%11.4e, %11.4e) count %d\n",omega,omegaold,energy,enold,count);

    ABORTXLAL(status);

  }

  *countback=count;

  /*Now ringdown is inserted*/

  XLALRungeKutta4Free( integrator );
  LALFree(dummy.data);
  
  /*--------------------------------------------------------------
   * Attach the ringdown waveform to the end of inspiral
   -------------------------------------------------------------*/  

  apcount= count;
  apcountM= count;
  xlalSt2P2 = XLALPSpinInspiralAttachRingdownWave( h2P2 , params , &apcount, nmodes , 2 , 2 , finalMass, finalSpin);
  xlalSt2M2 = XLALPSpinInspiralAttachRingdownWave( h2M2 , params , &apcountM, nmodes , 2 ,-2 , finalMass, finalSpin);
  if ((xlalSt2P2 != XLAL_SUCCESS)||(xlalSt2M2!= XLAL_SUCCESS)){
    XLALDestroyREAL4Vector( h2P2 );
    XLALDestroyREAL4Vector( h2M2 );
    XLALDestroyREAL4Vector( h2P1 );
    XLALDestroyREAL4Vector( h2M1 );
    XLALDestroyREAL4Vector( h3P3 );
    XLALDestroyREAL4Vector( h3M3 );
    XLALDestroyREAL4Vector( fap );
    XLALDestroyREAL8Vector( phap );
    XLALDestroyREAL4Vector( shift22 );
    ABORTXLAL( status );
  }
  else {
    for (i=2*apcount;i<2*length;i++) h2P2->data[i]=0.;
    if ( apcount > *countback ) *countback= apcount;
    for (i=2*apcountM;i<2*length;i++) h2M2->data[i]=0.;
    if ( apcountM > *countback ) *countback= apcountM;
  }

  apcount= count;
  errcode = XLALPSpinInspiralAttachRingdownWave( h2P1 , params , &apcount, nmodes , 2 , 1 , finalMass, finalSpin );
  if (errcode != XLAL_SUCCESS ) {
    XLALDestroyREAL4Vector( h2P1 );
  }
  else {
    for (i=2*apcount;i<2*length;i++) h2P1->data[i]=0.;
    if (apcount > *countback) *countback = apcount;
  }
  apcount= count;
  errcode = XLALPSpinInspiralAttachRingdownWave( h2M1 , params , &apcount, nmodes , 2 ,-1 , finalMass, finalSpin );
  if (errcode != XLAL_SUCCESS) {
    XLALDestroyREAL4Vector( h2M1 );
  }
  else {
    for (i=2*apcount;i<2*length;i++) h2M1->data[i]=0.;
    if (apcount > *countback) *countback = apcount;
  }

  apcount=count;
  errcode = XLALPSpinInspiralAttachRingdownWave( h20 , params , &apcount, nmodes , 2 , 0 , finalMass, finalSpin);
  if (errcode != XLAL_SUCCESS ) {
    XLALDestroyREAL4Vector( h20 );
  }
  else {
    for (i=2*apcount;i<2*length;i++) h20->data[i]=0.;
    if (apcount > *countback) *countback= apcount;
  }

  apcount= count;
  errcode = XLALPSpinInspiralAttachRingdownWave( h3P3 , params , &apcount, nmodes , 3 , 3 , finalMass, finalSpin );
  if (errcode != XLAL_SUCCESS ) XLALDestroyREAL4Vector( h3P3 );
  else {
    for (i=2*apcount;i<2*length;i++) h3P3->data[i]=0.;
    if (apcount > *countback ) *countback= apcount;
  }
  apcount= count;
  errcode = XLALPSpinInspiralAttachRingdownWave( h3M3 , params , &apcount, nmodes , 3 ,-3 , finalMass, finalSpin );
  if (errcode != XLAL_SUCCESS ) XLALDestroyREAL4Vector( h3M3 );
  else {
    for (i=2*apcount;i<2*length;i++) h3M3->data[i]=0.;
    if (apcount > *countback ) *countback= apcount;
  }

  /*-------------------------------------------------------------------
   * Compute the spherical harmonics required for constructing (h+,hx).
   -------------------------------------------------------------------*/

  /* The angles theta for the spherical harmonics has been set according to 
     the input inclination parameter and the axisChoice*/

  /* -----------------------------------------------------------------
   * Attaching the (2,\pm 2), (2,\pm 1), (2,0), (3,\pm 3), Spherical Harmonics
   *----------------------------------------*/

   xlalSt2P2 = XLALSphHarm( &MultSphHarm2P2, 2, 2, inc , 0. );
   xlalSt2M2 = XLALSphHarm( &MultSphHarm2M2, 2,-2, inc , 0. );
   if ((xlalSt2P2 != XLAL_SUCCESS )||(xlalSt2M2 != XLAL_SUCCESS))
   {
     XLALDestroyREAL4Vector( h2P2 );
     XLALDestroyREAL4Vector( h2M2 );
     XLALDestroyREAL4Vector( h2P1 );
     XLALDestroyREAL4Vector( h2M1 );
     XLALDestroyREAL4Vector( h20 );
     XLALDestroyREAL4Vector( h3P3 );
     XLALDestroyREAL4Vector( h3M3 );
     XLALDestroyREAL4Vector( fap );
     XLALDestroyREAL8Vector( phap );
     XLALDestroyREAL4Vector( shift22 );
     fprintf(stderr,"** LALPSpinInspiralRD ERROR **: impossible to create Y22\n");
     ABORTXLAL( status );
   }

   for ( i = 0; i < length; i++) {
     fap->data[i] /= unitHz;
     x0 = h2P2->data[2*i];
     x1 = h2P2->data[2*i+1];
     x2 = h2M2->data[2*i];
     x3 = h2M2->data[2*i+1];
     sig1->data[i] = x0*MultSphHarm2P2.re - x1*MultSphHarm2P2.im + x2*MultSphHarm2M2.re - x3*MultSphHarm2M2.im;
     sig2->data[i] = x0*MultSphHarm2P2.im + x1*MultSphHarm2P2.re + x2*MultSphHarm2M2.im - x3*MultSphHarm2M2.re;
   }

   errcode  = XLALSphHarm( &MultSphHarm2P1, 2, 1, inc , 0. );
   errcode2 = XLALSphHarm( &MultSphHarm2M1, 2, -1, inc , 0. );
   if ((errcode != XLAL_SUCCESS )||(errcode2 != XLAL_SUCCESS))
     {
       XLALDestroyREAL4Vector( h2P1 );
       XLALDestroyREAL4Vector( h2M1 );
       fprintf(stdout,"** WARNING **: impossible to create Y21\n");
     }
   else {
     for ( i = 0; i < length; i++)
       {
	 x0 = h2P1->data[2*i];
	 x1 = h2P1->data[2*i+1];
	 x2 = h2M1->data[2*i];
	 x3 = h2M1->data[2*i+1];
	 sig1->data[i]+= x0*MultSphHarm2P1.re - x1*MultSphHarm2P1.im + x2*MultSphHarm2M1.re - x3*MultSphHarm2M1.im;
	 sig2->data[i]+= x0*MultSphHarm2P1.im + x1*MultSphHarm2P1.re + x2*MultSphHarm2M1.im - x3*MultSphHarm2M1.re;
       }
   }

   errcode = XLALSphHarm( &MultSphHarm20, 2, 0, inc , 0. );
   if (errcode != XLAL_SUCCESS )
     {
       XLALDestroyREAL4Vector( h20 );
       fprintf(stdout,"** WARNING **: impossible to create Y20\n");
     }
   else {
     for ( i = 0; i < length; i++)
       {
	 x0 = h20->data[2*i];
	 x1 = h20->data[2*i+1];
	 sig1->data[i]+= x1*MultSphHarm20.re - x1*MultSphHarm20.im;
	 sig2->data[i]+= x1*MultSphHarm20.im + x1*MultSphHarm20.re;
       }
   }

   errcode = XLALSphHarm( &MultSphHarm3P3, 3, 3, inc , 0. );
   errcode2 = XLALSphHarm( &MultSphHarm3M3, 3, -3, inc , 0. );
   if ((errcode != XLAL_SUCCESS )||(errcode2 != XLAL_SUCCESS)) {
     XLALDestroyREAL4Vector( h3P3 );
     XLALDestroyREAL4Vector( h3M3 );
     fprintf(stdout,"** WARNING **: impossible to create Y33,Y3-3\n");
   }   
   else {
     for ( i = 0; i < length; i++)
       {
	 x0 = h3P3->data[2*i];
	 x1 = h3P3->data[2*i+1];
	 x2 = h3M3->data[2*i];
	 x3 = h3M3->data[2*i+1];
	 sig1->data[i]+= x0*MultSphHarm3P3.re - x1*MultSphHarm3P3.im + x2*MultSphHarm3M3.re - x3*MultSphHarm3M3.im;
	 sig2->data[i]+= x0*MultSphHarm3P3.im + x1*MultSphHarm3P3.re + x2*MultSphHarm3M3.im + x3*MultSphHarm3M3.re;
       }
   }

   params->fFinal = params->tSampling/2.;

   /*------------------------------------------------------
    * If required by the user copy other data sets to the
    * relevant arrays
    ------------------------------------------------------*/
   if (hh)
   {
     for(i = 0; i < length; i++)
     {
       j=2*i;
       k=2*i+1;
       hap->data[j] = sig1->data[i];
       hap->data[k] = sig2->data[i];
     }
   }

   if (signalvec1) memcpy(signalvec1->data , sig1->data,    length * (sizeof(REAL4)));
   if (signalvec2) memcpy(signalvec2->data , sig2->data,    length * (sizeof(REAL4)));
   if (hh)         memcpy(hh->data         , hap->data,   2*length * (sizeof(REAL4))); 
   if (ff)         memcpy(ff->data         , fap->data,     length * (sizeof(REAL4)));
   if (phi)        memcpy(phi->data        , phap->data,    length * (sizeof(REAL8)));
   if (shift)      memcpy(shift->data      , shift22->data, length * (sizeof(REAL4)));

   /* Clean up */
   XLALDestroyREAL4Vector ( h2P2 );
   XLALDestroyREAL4Vector ( h2M2 );
   XLALDestroyREAL4Vector ( h2P1 );
   XLALDestroyREAL4Vector ( h2M1 );
   XLALDestroyREAL4Vector ( h20 );
   XLALDestroyREAL4Vector ( h3P3 );
   XLALDestroyREAL4Vector ( h3M3 );
   XLALDestroyREAL4Vector ( shift22 );
   XLALDestroyREAL4Vector ( fap );
   XLALDestroyREAL8Vector ( phap );
   XLALDestroyREAL4Vector ( sig1 );
   XLALDestroyREAL4Vector ( sig2 );
   XLALDestroyREAL4Vector ( hap );

   DETATCHSTATUSPTR(status);
   RETURN(status);

  /*End*/
   
}
