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

/*  <lalVerbatim file="LALTaylorNWaveformCV">
Author: Sathyaprakash, B. S., Cokelaer T.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALTaylorNWaveform.c} and
\texttt{LALTaylorNWaveformTemplates.c}}

\vfill{\footnotesize\input{LALTaylorNWaveformCV}}

</lalLaTeX>  */
/*-------------------------------------------------------------------*/
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/FindRoot.h>

typedef struct tagxiInitIn {
  REAL8 eta, omega, e0;
} xiInitIn;

static void LALxiInit4PN(
  LALStatus *status,
  REAL8     *x,
  REAL8      xi,
  void      *params);


static void LALxiInit5PN(
  LALStatus *status,
  REAL8     *x,
  REAL8      xi,
  void      *params);


static void LALxiInit6PN(
  LALStatus *status,
  REAL8     *x,
  REAL8      xi,
  void      *params);


static void LALxiInit7PN(
  LALStatus *status,
  REAL8     *x,
  REAL8      xi,
  void      *params);


void LALTaylorNDerivatives4PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);
void LALTaylorNDerivatives5PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);
void LALTaylorNDerivatives6PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);
void LALTaylorNDerivatives7PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);

void
LALTaylorNWaveformEngine (
  LALStatus        *status,
  REAL4Vector      *signalvec,
  InspiralTemplate *params,
  InspiralInit     *paramsInit
);

NRCSID (LALTAYLORETWAVEFORMC,
"$Id$");

static void LALxiInit4PN(
   LALStatus *status,
   REAL8     *x,
   REAL8      xi,
   void      *params)
{

   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   xiInitIn *in;
   REAL8 xi2, xi43, xi23, eta, eta2, pisq;

   status = NULL;
   in = (xiInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   xi2 = xi*xi;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   pisq = LAL_PI*LAL_PI;

   *x = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43 );
   *x -= in->omega;
}

static void LALxiInit5PN(
   LALStatus *status,
   REAL8     *x,
   REAL8      xi,
   void      *params)
{

   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   xiInitIn *in;
   REAL8 xi2, xi43, xi23, eta, eta2, pisq;

   status = NULL;
   in = (xiInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   xi2 = xi*xi;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   pisq = LAL_PI*LAL_PI;

   *x = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43 );
   *x -= in->omega;
}

static void LALxiInit6PN(
   LALStatus *status,
   REAL8     *x,
   REAL8      xi,
   void      *params)
{

   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   xiInitIn *in;
   REAL8 xi2, xi43, xi23, eta, eta2, pisq;

   status = NULL;
   in = (xiInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   xi2 = xi*xi;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   pisq = LAL_PI*LAL_PI;

   *x = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43
        + (315./2. + (-817./4. + 123./32.*pisq) * eta + 7*eta2) * xi2);
   *x -= in->omega;
}

static void LALxiInit7PN(
   LALStatus *status,
   REAL8     *x,
   REAL8      xi,
   void      *params)
{

   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   xiInitIn *in;
   REAL8 xi2, xi43, xi23, eta, eta2, pisq;

   status = NULL;
   in = (xiInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   xi2 = xi*xi;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   pisq = LAL_PI*LAL_PI;

   *x = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43
        + (315./2. + (-817./4. + 123./32.*pisq) * eta + 7*eta2) * xi2);
   *x -= in->omega;
}


void LALTaylorNDerivatives4PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   InspiralDerivativesIn *ak;
   REAL8 xi, xi2, xi23, xi43, xi53, xi73, xi113, eta, eta2, pisq, gamma;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   xi = values->data[1];
   eta2 = eta * eta;
   gamma = 0.577215664901532;
   xi2 = xi*xi;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   xi53 = pow(xi, 5./3.);
   xi73 = pow(xi, 7./3.);
   xi113 = pow(xi, 11./3.);
   pisq = LAL_PI*LAL_PI;

   dvalues->data[0] = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43);

   dvalues->data[1] = 96./5.*eta*xi113*(1. + (1273./336. - 11./4.*eta) * xi23
		   + 4.*LAL_PI*xi + (438887./18144. - 49507./2016.*eta + 59./18.*eta2)*xi43);
}

void LALTaylorNDerivatives5PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   InspiralDerivativesIn *ak;
   REAL8 xi, xi2, xi23, xi43, xi53, xi73, xi113, eta, eta2, pisq, gamma;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   xi = values->data[1];
   eta2 = eta * eta;
   gamma = 0.577215664901532;
   xi2 = xi*xi;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   xi53 = pow(xi, 5./3.);
   xi73 = pow(xi, 7./3.);
   xi113 = pow(xi, 11./3.);
   pisq = LAL_PI*LAL_PI;


   dvalues->data[0] = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43);
   dvalues->data[1] = 96./5.*eta*xi113*(1. + (1273./336. - 11./4.*eta) * xi23
		   + 4.*LAL_PI*xi + (438887./18144. - 49507./2016.*eta + 59./18.*eta2)*xi43
		   + (20033./672. - 189./8.*eta) * LAL_PI * xi53);
}

void LALTaylorNDerivatives6PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   InspiralDerivativesIn *ak;
   REAL8 xi, xi2, xi13, xi23, xi43, xi53, xi73, xi113, eta, eta2, eta3, pisq, gamma;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   xi = values->data[1];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
   gamma = 0.577215664901532;
   xi2 = xi*xi;
   xi13 = pow(xi, 1./3.);
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   xi53 = pow(xi, 5./3.);
   xi73 = pow(xi, 7./3.);
   xi113 = pow(xi, 11./3.);
   pisq = LAL_PI*LAL_PI;

   dvalues->data[0] = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43
		   + (315./2. + (-817./4. + 123./32.*pisq) * eta + 7*eta2) * xi2);

   dvalues->data[1] = 96./5.*eta*xi113*(1. + (1273./336. - 11./4.*eta) * xi23
		   + 4.*LAL_PI*xi + (438887./18144. - 49507./2016.*eta + 59./18.*eta2)*xi43
		   + (20033./672. - 189./8.*eta) * LAL_PI * xi53 + (38047038863./139708800.
			   + (16./3. + 287./24. * eta)*pisq - 16554367./31104.*eta
			   + 617285./8064.*eta2 - 5605./2592.*eta3 - 1712./105.*(gamma + log(4*xi13)))*xi2);
}

void LALTaylorNDerivatives7PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   InspiralDerivativesIn *ak;
   REAL8 xi, xi2, xi13, xi23, xi43, xi53, xi73, xi113, eta, eta2, eta3, pisq, gamma;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   xi = values->data[1];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
   gamma = 0.577215664901532;
   xi2 = xi*xi;
   xi13 = pow(xi, 1./3.);
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   xi53 = pow(xi, 5./3.);
   xi73 = pow(xi, 7./3.);
   xi113 = pow(xi, 11./3.);
   pisq = LAL_PI*LAL_PI;

   dvalues->data[0] = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43
		   + (315./2. + (-817./4. + 123./32.*pisq) * eta + 7*eta2) * xi2);

   dvalues->data[1] = 96./5.*eta*xi113*(1. + (1273./336. - 11./4.*eta) * xi23
		   + 4.*LAL_PI*xi + (438887./18144. - 49507./2016.*eta + 59./18.*eta2)*xi43
		   + (20033./672. - 189./8.*eta) * LAL_PI * xi53 + (38047038863./139708800.
			   + (16./3. + 287./24. * eta)*pisq - 16554367./31104.*eta
			   + 617285./8064.*eta2 - 5605./2592.*eta3 - 1712./105.*(gamma + log(4*xi13)))*xi2
		   + (971011./4032. - 1608185./6048.*eta + 91495./1512.*eta2)*LAL_PI*xi73);
}


/*  <lalVerbatim file="LALTaylorNWaveformCP"> */
void LALTaylorNWaveform (
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{ /* </lalVerbatim> */

   InspiralInit paramsInit;
   INITSTATUS(status, "LALTaylorNWaveform", LALTAYLORETWAVEFORMC);
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
   LALTaylorNWaveformEngine(status->statusPtr, signalvec, params, &paramsInit);
   CHECKSTATUSPTR( status );

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/*---------------------------------------------------------*/
void
LALTaylorNWaveformEngine (
                LALStatus        *status,
                REAL4Vector      *signalvec,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                )
{

   void                  *funcParams;
   UINT4                 length=0, count, ndx;
   INT4                  nn=2;
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   REAL8                 h, omega, omegaMax, t, dt, m, eta, phi, xi, v;
   REAL8Vector           dummy, values, dvalues, newvalues, yt, dym, dyt;
   rk4GSLIntegrator      *integrator = NULL;
   InspiralDerivativesIn in2;
   xiInitIn              in3;
   rk4In                 in4;
   expnCoeffs            ak;
   expnFunc              func;
   DFindRootIn           rootIn;
   CHAR			 message[256];


   INITSTATUS(status, "LALTaylorNWaveformEngine", LALTAYLORETWAVEFORMC);
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
   /* Given omega compute xi by solving Eq.(18b) of TG08 */
   omega = LAL_PI * params->fLower * m;
   rootIn.xacc = 1.0e-16;
   rootIn.xmax = 2.*omega;
   rootIn.xmin = omega/2.;
   switch (params->order)
   {
	case LAL_PNORDER_TWO:
	   rootIn.function = LALxiInit4PN;
	   break;
	case LAL_PNORDER_TWO_POINT_FIVE:
	   rootIn.function = LALxiInit5PN;
	   break;
	case LAL_PNORDER_THREE:
	   rootIn.function = LALxiInit6PN;
	   break;
	case LAL_PNORDER_THREE_POINT_FIVE:
	   rootIn.function = LALxiInit7PN;
	   break;
	default:
	   LALSnprintf(message, 256, "There are no Et waveforms at order %d\n", params->order);
	   LALError( status, message );
	   LALFree(dummy.data);
	   ABORT( status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);
	   break;
   }
   in3.eta = ak.eta;
   in3.omega = omega;
   funcParams = (void *) &in3;


   LALDBisectionFindRoot(status->statusPtr, &xi, &rootIn, funcParams);
   CHECKSTATUSPTR(status);

   /* End of initial conditions */

   values.data[0] = phi = params->startPhase;
   values.data[1] = xi;

   /* fprintf(stdout, "Initial phi=%e, xi=%e\n", values.data[0], values.data[1]); */
   t = 0.0;

   /* Initialize the GSL integrator */
   switch (params->order)
   {
	case LAL_PNORDER_TWO:
	   in4.function = LALTaylorNDerivatives4PN;
	   break;
	case LAL_PNORDER_TWO_POINT_FIVE:
	   in4.function = LALTaylorNDerivatives5PN;
	   break;
	case LAL_PNORDER_THREE:
	   in4.function = LALTaylorNDerivatives6PN;
	   break;
	case LAL_PNORDER_THREE_POINT_FIVE:
	   in4.function = LALTaylorNDerivatives7PN;
	   break;
	default:
	   LALSnprintf(message, 256, "There are no Et waveforms at order %d\n", params->order);
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

   /* fprintf(stdout, "fMin=%e, fMax=%e, f=%e, dxi/dt=%e\n",  */
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

      h = 4 * m * eta * xi * sin(2.*phi)/1.e14;
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
      xi  = values.data[1] = newvalues.data[1];
      /* fprintf(stdout, "phi=%e, xi=%e, omega=%e, dxi/dt=%e\n",  */
      /* phi, xi, omega/(m*LAL_PI), dvalues.data[1]);             */

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
