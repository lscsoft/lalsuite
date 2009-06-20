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

/*  <lalVerbatim file="LALTaylorEtWaveformCV">
Author: Sathyaprakash, B. S., Cokelaer T.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTaylorEtWaveform.c} and
\texttt{LALTaylorEtWaveformTemplates.c}}

\vfill{\footnotesize\input{LALTaylorEtWaveformCV}}

</lalLaTeX>  */
/*-------------------------------------------------------------------*/
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/FindRoot.h>
#include <lal/LALConstants.h>

typedef struct tagzetaInitIn {
  REAL8 eta, omega, e0;
} zetaInitIn;

static void LALzetaInit4PN(
  LALStatus *status,
  REAL8     *x,
  REAL8      zeta,
  void      *params);

static void LALzetaInit5PN(
  LALStatus *status,
  REAL8     *x,
  REAL8      zeta,
  void      *params);

static void LALzetaInit6PN(
  LALStatus *status,
  REAL8     *x,
  REAL8      zeta,
  void      *params);

static void LALzetaInit7PN(
  LALStatus *status,
  REAL8     *x,
  REAL8      zeta,
  void      *params);


void LALTaylorEtDerivatives4PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);
void LALTaylorEtDerivatives5PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);
void LALTaylorEtDerivatives6PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);
void LALTaylorEtDerivatives7PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);

void LALTaylorEtWaveform (
  LALStatus        *status,
  REAL4Vector      *signalvec,
  InspiralTemplate *params
);

void
LALTaylorEtWaveformEngine (
  LALStatus        *status,
  REAL4Vector      *signalvec,
  InspiralTemplate *params,
  InspiralInit     *paramsInit
);

NRCSID (LALTAYLORETWAVEFORMC,
"$Id$");

static void LALzetaInit4PN(
   LALStatus *status,
   REAL8     *x,
   REAL8      zeta,
   void      *params)
{
   zetaInitIn *in;
   REAL8 zeta2, zeta3, zeta32, eta, eta2, eta3, pisq;

   status = NULL;
   in = (zetaInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   eta3 = eta2*eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta32 = pow(zeta, 1.5);
   pisq = LAL_PI*LAL_PI;

   *x = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2 );
   *x -= in->omega;
}

static void LALzetaInit5PN(
   LALStatus *status,
   REAL8     *x,
   REAL8      zeta,
   void      *params)
{
   zetaInitIn *in;
   REAL8 zeta2, zeta3, zeta32, eta, eta2, eta3, pisq;

   status = NULL;
   in = (zetaInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   eta3 = eta2*eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta32 = pow(zeta, 1.5);
   pisq = LAL_PI*LAL_PI;

   *x = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2 );
   *x -= in->omega;
}

static void LALzetaInit6PN(
   LALStatus *status,
   REAL8     *x,
   REAL8      zeta,
   void      *params)
{
   zetaInitIn *in;
   REAL8 zeta2, zeta3, zeta32, eta, eta2, eta3, pisq;

   status = NULL;
   in = (zetaInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   eta3 = eta2*eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta32 = pow(zeta, 1.5);
   pisq = LAL_PI*LAL_PI;

   *x = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2
     + (41445./1024. - (309715./3072. - 205./64.*pisq) * eta
     + 1215./1024.*eta2 + 45./1024*eta3) * zeta3);
   *x -= in->omega;
}

static void LALzetaInit7PN(
   LALStatus *status,
   REAL8     *x,
   REAL8      zeta,
   void      *params)
{
   zetaInitIn *in;
   REAL8 zeta2, zeta3, zeta32, eta, eta2, eta3, pisq;

   status = NULL;
   in = (zetaInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   eta3 = eta2*eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta32 = pow(zeta, 1.5);
   pisq = LAL_PI*LAL_PI;

   *x = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2
     + (41445./1024. - (309715./3072. - 205./64.*pisq) * eta
     + 1215./1024.*eta2 + 45./1024*eta3) * zeta3);
   *x -= in->omega;
}

void LALTaylorEtDerivatives4PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   InspiralDerivativesIn *ak;
   REAL8 zeta, zeta2, zeta3, zeta5, zeta32, zeta52, zeta72, eta, eta2, eta3, pisq, fourpi;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   zeta = values->data[1];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta5 = zeta2*zeta3;
   zeta32 = pow(zeta, 1.5);
   zeta52 = pow(zeta, 2.5);
   zeta72 = pow(zeta, 3.5);
   pisq = LAL_PI*LAL_PI;
   fourpi = 4.*LAL_PI;

   dvalues->data[0] = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2);
   /*
   dvalues->data[1] = 64./5.* eta * zeta5 * (1. + (13./336. - 2.5*eta)*zeta
     + fourpi*zeta32 + (117857./18144. - 12017./2016.*eta + 2.5*eta2) * zeta2);
   */
   dvalues->data[1] = 64./5.* eta * zeta5 * (1. + (13./336. - 2.5*eta)*zeta
     + fourpi*zeta32 + (117857./18144. - 12017./2016.*eta + 2.5*eta2) * zeta2);
}

void LALTaylorEtDerivatives5PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   InspiralDerivativesIn *ak;
   REAL8 zeta, zeta2, zeta3, zeta5, zeta32, zeta52, zeta72, eta, eta2, eta3, pisq, fourpi;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   zeta = values->data[1];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta5 = zeta2*zeta3;
   zeta32 = pow(zeta, 1.5);
   zeta52 = pow(zeta, 2.5);
   zeta72 = pow(zeta, 3.5);
   pisq = LAL_PI*LAL_PI;
   fourpi = 4.*LAL_PI;

   dvalues->data[0] = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2 );
   /*
   dvalues->data[1] = 64./5.* eta * zeta5 * (1. + (13./336. - 2.5*eta)*zeta
     + fourpi*zeta32 + (117857./18144. - 12017./2016.*eta + 2.5*eta2) * zeta2);
   */
   dvalues->data[1] = 64./5.* eta * zeta5 * (1. + (13./336. - 2.5*eta)*zeta
     + fourpi*zeta32 + (117857./18144. - 12017./2016.*eta + 2.5*eta2) * zeta2
     + (4913.*LAL_PI/672. - 177.*LAL_PI*eta/8.) * zeta52 );
}

void LALTaylorEtDerivatives6PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   InspiralDerivativesIn *ak;
   REAL8 zeta, zeta2, zeta3, zeta5, zeta32, zeta52, zeta72, eta, eta2, eta3, pisq, fourpi;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   zeta = values->data[1];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta5 = zeta2*zeta3;
   zeta32 = pow(zeta, 1.5);
   zeta52 = pow(zeta, 2.5);
   zeta72 = pow(zeta, 3.5);
   pisq = LAL_PI*LAL_PI;
   fourpi = 4.*LAL_PI;

   dvalues->data[0] = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2
     + (41445./1024. - (309715./3072. - 205./64.*pisq) * eta
     + 1215./1024.*eta2 + 45./1024*eta3) * zeta3);
   /*
   dvalues->data[1] = 64./5.* eta * zeta5 * (1. + (13./336. - 2.5*eta)*zeta
     + fourpi*zeta32 + (117857./18144. - 12017./2016.*eta + 2.5*eta2) * zeta2);
   */
   dvalues->data[1] = 64./5.* eta * zeta5 * (1. + (13./336. - 2.5*eta)*zeta
     + fourpi*zeta32 + (117857./18144. - 12017./2016.*eta + 2.5*eta2) * zeta2
     + (4913.*LAL_PI/672. - 177.*LAL_PI*eta/8.) * zeta52
     + (-85./64.*eta3 + 488849./16128.*eta2 + 369.*pisq*eta/32.
	     - 24861497.*eta/72576. - 856.*log(16.*zeta)/105. + 16.*pisq/3.
	     - 1712.*LAL_GAMMA/105. + 37999588601./279417600.) * zeta3 );
}

void LALTaylorEtDerivatives7PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   InspiralDerivativesIn *ak;
   REAL8 zeta, zeta2, zeta3, zeta5, zeta32, zeta52, zeta72, eta, eta2, eta3, pisq, fourpi;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   zeta = values->data[1];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta5 = zeta2*zeta3;
   zeta32 = pow(zeta, 1.5);
   zeta52 = pow(zeta, 2.5);
   zeta72 = pow(zeta, 3.5);
   pisq = LAL_PI*LAL_PI;
   fourpi = 4.*LAL_PI;

   dvalues->data[0] = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2
     + (41445./1024. - (309715./3072. - 205./64.*pisq) * eta
     + 1215./1024.*eta2 + 45./1024*eta3) * zeta3);
   /*
   dvalues->data[1] = 64./5.* eta * zeta5 * (1. + (13./336. - 2.5*eta)*zeta
     + fourpi*zeta32 + (117857./18144. - 12017./2016.*eta + 2.5*eta2) * zeta2);
   */
   dvalues->data[1] = 64./5.* eta * zeta5 * (1. + (13./336. - 2.5*eta)*zeta
     + fourpi*zeta32 + (117857./18144. - 12017./2016.*eta + 2.5*eta2) * zeta2
     + (4913.*LAL_PI/672. - 177.*LAL_PI*eta/8.) * zeta52
     + (-85./64.*eta3 + 488849./16128.*eta2 + 369.*pisq*eta/32.
	     - 24861497.*eta/72576. - 856.*log(16.*zeta)/105. + 16.*pisq/3.
	     - 1712.*LAL_GAMMA/105. + 37999588601./279417600.) * zeta3
     + (613373.*LAL_PI*eta2/12096. - 3207739.*LAL_PI*eta/48384.
	     + 129817.*LAL_PI/2304.) * zeta72);
}


/*  <lalVerbatim file="LALTaylorEtWaveformCP"> */
void LALTaylorEtWaveform (
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{ /* </lalVerbatim> */

   InspiralInit paramsInit;
   INITSTATUS(status, "LALTaylorEtWaveform", LALTAYLORETWAVEFORMC);
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
   LALTaylorEtWaveformEngine(status->statusPtr, signalvec, params, &paramsInit);
   CHECKSTATUSPTR( status );

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/*---------------------------------------------------------*/
void
LALTaylorEtWaveformEngine (
                LALStatus        *status,
                REAL4Vector      *signalvec,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                )
{

   void                  *funcParams;
   UINT4                 length=0, count, ndx;
   INT4                  nn=2;
   REAL8                 h, omega, omegaMax, t, dt, m, eta, phi, zeta, v;
   REAL8Vector           dummy, values, dvalues, newvalues, yt, dym, dyt;
   rk4GSLIntegrator      *integrator = NULL;
   InspiralDerivativesIn in2;
   zetaInitIn              in3;
   rk4In                 in4;
   expnCoeffs            ak;
   expnFunc              func;
   DFindRootIn           rootIn;
   CHAR			 message[256];


   INITSTATUS(status, "LALTaylorEtWaveformEngine", LALTAYLORETWAVEFORMC);
   ATTATCHSTATUSPTR(status);

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


   /* Set dt to sampling interval specified by user */

   dt = 1./params->tSampling;
   ak   = paramsInit->ak;
   func = paramsInit->func;
   length = signalvec->length;
   eta = ak.eta;
   m = ak.totalmass;

   /* Begin initial conditions */
   /* Given omega compute zeta by solving Eq.(18b) of TG08 */
   omega = LAL_PI * params->fLower * m;
   v = pow(omega,1./3.);
   rootIn.xacc = 1.0e-16;
   rootIn.xmax = 1.;
   rootIn.xmin = v*v/2.;

   /* Initialize the GSL integrator */
   switch (params->order)
   {
	case LAL_PNORDER_TWO:
	   rootIn.function = LALzetaInit4PN;
	   break;
	case LAL_PNORDER_TWO_POINT_FIVE:
	   rootIn.function = LALzetaInit5PN;
	   break;
	case LAL_PNORDER_THREE:
	   rootIn.function = LALzetaInit6PN;
	   break;
	case LAL_PNORDER_THREE_POINT_FIVE:
	   rootIn.function = LALzetaInit7PN;
	   break;
	default:
	   snprintf(message, 256, "There are no Et waveforms at order %d\n", params->order);
	   LALError( status, message );
	   LALFree(dummy.data);
	   ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
	   break;
   }
   in3.eta = ak.eta;
   in3.omega = omega;
   funcParams = (void *) &in3;


   LALDBisectionFindRoot(status->statusPtr, &zeta, &rootIn, funcParams);
   CHECKSTATUSPTR(status);

   /* End of initial conditions */

   values.data[0] = phi = params->startPhase;
   values.data[1] = zeta;

   /* fprintf(stdout, "Initial phi=%e, zeta=%e\n", values.data[0], values.data[1]); */
   t = 0.0;

   /* Initialize the GSL integrator */
   switch (params->order)
   {
	case LAL_PNORDER_TWO:
	   in4.function = LALTaylorEtDerivatives4PN;
	   break;
	case LAL_PNORDER_TWO_POINT_FIVE:
	   in4.function = LALTaylorEtDerivatives5PN;
	   break;
	case LAL_PNORDER_THREE:
	   in4.function = LALTaylorEtDerivatives6PN;
	   break;
case LAL_PNORDER_THREE_POINT_FIVE:
	   in4.function = LALTaylorEtDerivatives7PN;
	   break;
	default:
	   snprintf(message, 256, "There are no Et waveforms at order %d\n", params->order);
	   LALError( status, message );
	   LALFree(dummy.data);
	   ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
	   break;
   }
   in4.y = &values;
   in4.h = dt/m;
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

   omegaMax = 1./pow(6.,1.5);
   t = 0.0;
   ndx = 0;
   count = 0;

   /* fprintf(stdout, "fMin=%e, fMax=%e, f=%e, dzeta/dt=%e\n",  */
   /* omega/(m*LAL_PI), omegaMax/(m*LAL_PI),                  */
   /* dvalues.data[0]/(m*LAL_PI), dvalues.data[1]);           */

   while (omega<omegaMax)
   {
      if (count > length)
      {
        XLALRungeKutta4Free( integrator );
	LALFree(dummy.data);
	ABORT(status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
      }

      h = 4 * m * eta * zeta * sin(2.*phi)/1.e14;
      signalvec->data[ndx] = h;
      /* fprintf(stdout, "%e %e %e\n", t, h, omega/(m*LAL_PI)); */

      /* Integrate one step forward */
      in4.dydx = &dvalues;
      in4.x = t/m;
      LALRungeKutta4(status->statusPtr, &newvalues, integrator, funcParams);
      BEGINFAIL( status )
      {
        XLALRungeKutta4Free( integrator );
        LALFree(dummy.data);
      }
      ENDFAIL( status );

      /* Update the values of the dynamical variables */
      phi = values.data[0] = newvalues.data[0];
      zeta  = values.data[1] = newvalues.data[1];
      /* fprintf(stdout, "phi=%e, zeta=%e, omega=%e, dzeta/dt=%e\n",  */
      /* phi, zeta, omega/(m*LAL_PI), dvalues.data[1]);             */

      /* Compute the derivaties at the new location */
      in4.function(&values, &dvalues, funcParams);
      omega = dvalues.data[0];

      t = (++count-params->nStartPad) * dt;
      ndx++;
   }

   /*----------------------------------------------------------------------*/
   /* Record the final cutoff frequency of BD Waveforms for record keeping */
   /* ---------------------------------------------------------------------*/
   v = pow(omega, 1./3.);
   params->vFinal = v;
   params->tC = t;
   params->fFinal = omega/(LAL_PI*m);
   /* fprintf(stdout, "Final velocity=%e, time=%e, frequency=%e\n", v, t, params->fFinal); */

   XLALRungeKutta4Free( integrator );
   LALFree(dummy.data);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
/*---------------------------------------------------------*/
