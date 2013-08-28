/*
*  Copyright (C) 2008 B.S. Sathyaprakash, Bala R. Iyer, Drew Keppel
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
 * \author Sathyaprakash, B. S., Cokelaer T.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief NONE
 *
 */
/*-------------------------------------------------------------------*/
#include <lal/LALInspiral.h>
#include <lal/SeqFactories.h>
#include <lal/FindRoot.h>
#include <lal/LALConstants.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

typedef struct tagxiInitIn {
  REAL8 eta, omega, e0;
} xiInitIn;

static REAL8 XLALxiInit4PN(
  REAL8      xi,
  void      *params);


static REAL8 XLALxiInit5PN(
  REAL8      xi,
  void      *params);


static REAL8 XLALxiInit6PN(
  REAL8      xi,
  void      *params);


static REAL8 XLALxiInit7PN(
  REAL8      xi,
  void      *params);


void XLALTaylorNDerivatives4PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);

void XLALTaylorNDerivatives5PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);

void XLALTaylorNDerivatives6PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);

void XLALTaylorNDerivatives7PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);

int
XLALTaylorNWaveformEngine (
  REAL4Vector      *output1,
  REAL4Vector      *output2,
  InspiralTemplate *params,
  InspiralInit     *paramsInit
);

static REAL8 XLALxiInit4PN(
   REAL8      xi,
   void      *params)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   xiInitIn *in;
   REAL8 xi43, xi23, eta, x;

   in = (xiInitIn *) params;
   eta = in->eta;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);

   x = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43 );
   x -= in->omega;

   return x;
}

static REAL8 XLALxiInit5PN(
   REAL8      xi,
   void      *params)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   xiInitIn *in;
   REAL8 xi43, xi23, eta, x;

   in = (xiInitIn *) params;
   eta = in->eta;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);

   x = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43 );
   x -= in->omega;

   return x;
}

static REAL8 XLALxiInit6PN(
   REAL8      xi,
   void      *params)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   xiInitIn *in;
   REAL8 xi2, xi43, xi23, eta, eta2, pisq, x;

   in = (xiInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   xi2 = xi*xi;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   pisq = LAL_PI*LAL_PI;

   x = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43
        + (315./2. + (-817./4. + 123./32.*pisq) * eta + 7*eta2) * xi2);
   x -= in->omega;

   return x;
}

static REAL8 XLALxiInit7PN(
   REAL8      xi,
   void      *params)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   xiInitIn *in;
   REAL8 xi2, xi43, xi23, eta, eta2, pisq, x;

   in = (xiInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   xi2 = xi*xi;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   pisq = LAL_PI*LAL_PI;

   x = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43
        + (315./2. + (-817./4. + 123./32.*pisq) * eta + 7*eta2) * xi2);
   x -= in->omega;

   return x;
}


void XLALTaylorNDerivatives4PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   InspiralDerivativesIn *ak;
   REAL8 xi, xi23, xi43, xi113, eta, eta2;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   xi = values->data[1];
   eta2 = eta * eta;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   xi113 = pow(xi, 11./3.);

   dvalues->data[0] = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43);

   dvalues->data[1] = 96./5.*eta*xi113*(1. + (1273./336. - 11./4.*eta) * xi23
		   + 4.*LAL_PI*xi + (438887./18144. - 49507./2016.*eta + 59./18.*eta2)*xi43);
}

void XLALTaylorNDerivatives5PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   InspiralDerivativesIn *ak;
   REAL8 xi, xi23, xi43, xi53, xi113, eta, eta2;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   xi = values->data[1];
   eta2 = eta * eta;
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   xi53 = pow(xi, 5./3.);
   xi113 = pow(xi, 11./3.);

   dvalues->data[0] = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43);
   dvalues->data[1] = 96./5.*eta*xi113*(1. + (1273./336. - 11./4.*eta) * xi23
		   + 4.*LAL_PI*xi + (438887./18144. - 49507./2016.*eta + 59./18.*eta2)*xi43
		   + (20033./672. - 189./8.*eta) * LAL_PI * xi53);
}

void XLALTaylorNDerivatives6PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   InspiralDerivativesIn *ak;
   REAL8 xi, xi2, xi13, xi23, xi43, xi53, xi113, eta, eta2, eta3, pisq;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   xi = values->data[1];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
   xi2 = xi*xi;
   xi13 = pow(xi, 1./3.);
   xi23 = pow(xi, 2./3.);
   xi43 = pow(xi, 4./3.);
   xi53 = pow(xi, 5./3.);
   xi113 = pow(xi, 11./3.);
   pisq = LAL_PI*LAL_PI;

   dvalues->data[0] = xi * ( 1. + 3. * xi23 + (39./2. - 7*eta) * xi43
		   + (315./2. + (-817./4. + 123./32.*pisq) * eta + 7*eta2) * xi2);

   dvalues->data[1] = 96./5.*eta*xi113*(1. + (1273./336. - 11./4.*eta) * xi23
		   + 4.*LAL_PI*xi + (438887./18144. - 49507./2016.*eta + 59./18.*eta2)*xi43
		   + (20033./672. - 189./8.*eta) * LAL_PI * xi53 + (38047038863./139708800.
			   + (16./3. + 287./24. * eta)*pisq - 16554367./31104.*eta
			   + 617285./8064.*eta2 - 5605./2592.*eta3 - 1712./105.*(LAL_GAMMA + log(4*xi13)))*xi2);
}

void XLALTaylorNDerivatives7PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   InspiralDerivativesIn *ak;
   REAL8 xi, xi2, xi13, xi23, xi43, xi53, xi73, xi113, eta, eta2, eta3, pisq;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   xi = values->data[1];
   eta2 = eta * eta;
   eta3 = eta2 * eta;
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
			   + 617285./8064.*eta2 - 5605./2592.*eta3 - 1712./105.*(LAL_GAMMA + log(4*xi13)))*xi2
		   + (971011./4032. - 1608185./6048.*eta + 91495./1512.*eta2)*LAL_PI*xi73);
}

void LALTaylorNWaveform (
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{
   XLALPrintDeprecationWarning("LALTaylorNWaveform", 
         "XLALTaylorNWaveform");
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   if (XLALTaylorNWaveform(signalvec, params))
      ABORTXLAL(status);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

int XLALTaylorNWaveform (
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

   if (XLALInspiralInit(params, &paramsInit))
      XLAL_ERROR(XLAL_EFUNC);

   if (params->totalMass <= 0.)
      XLAL_ERROR(XLAL_EDOM);
   if (params->eta < 0.)
      XLAL_ERROR(XLAL_EDOM);

   memset( output->data, 0, output->length * sizeof(REAL4) );

   /* Call the engine function */
   count = XLALTaylorNWaveformEngine(output, NULL, params, &paramsInit);
   if (count < 0)
      XLAL_ERROR(XLAL_EFUNC);

   return XLAL_SUCCESS;
}

/*---------------------------------------------------------*/
int
XLALTaylorNWaveformEngine (
                REAL4Vector      *output1,
                REAL4Vector      *output2,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                )
{

   void                  *funcParams;
   UINT4                 length=0, count, ndx;
   INT4                  nn=2;
   /* The variable xi=GMn/c^3, where omega = n (1 + k) */
   REAL8                 h, omega, omegaMax, t, dt, m, eta, phi, xi, v, xmax, xmin, xacc;
   REAL8Vector           dummy, values, dvalues, newvalues, yt, dym, dyt;
   rk4GSLIntegrator      *integrator = NULL;
   InspiralDerivativesIn in2;
   xiInitIn              in3;
   rk4In                 in4;
   expnCoeffs            ak;
   expnFunc              func;
   REAL8 (*rootfunction)(REAL8, void *);
   CHAR			 message[256];


/* Allocate all the memory required to dummy and then point the various
   arrays to dummy - this makes it easier to handle memory failures */

   dummy.length = nn * 6;

   if (!(dummy.data = (REAL8 * ) XLALMalloc(sizeof(REAL8) * nn * 6)))
      XLAL_ERROR(XLAL_ENOMEM);

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
   length = output1->length;
   eta = ak.eta;
   m = ak.totalmass;

   /* Begin initial conditions */
   /* Given omega compute xi by solving Eq.(18b) of TG08 */
   omega = LAL_PI * params->fLower * m;
   xacc = 1.0e-16;
   xmax = 2.*omega;
   xmin = omega/2.;
   switch (params->order)
   {
	case LAL_PNORDER_TWO:
	   rootfunction = &XLALxiInit4PN;
	   break;
	case LAL_PNORDER_TWO_POINT_FIVE:
	   rootfunction = &XLALxiInit5PN;
	   break;
	case LAL_PNORDER_THREE:
	   rootfunction = &XLALxiInit6PN;
	   break;
	case LAL_PNORDER_THREE_POINT_FIVE:
	   rootfunction = &XLALxiInit7PN;
	   break;
	default:
	   snprintf(message, 256, "There are no Et waveforms at order %d\n", params->order);
	   XLALPrintError( message );
	   XLALFree(dummy.data);
	   XLAL_ERROR(XLAL_EINVAL);
	   break;
   }
   in3.eta = ak.eta;
   in3.omega = omega;
   funcParams = (void *) &in3;


   xi = XLALDBisectionFindRoot(rootfunction, xmin, xmax, xacc, funcParams);
   if (XLAL_IS_REAL8_FAIL_NAN(xi))
      XLAL_ERROR(XLAL_EFUNC);

   /* End of initial conditions */

   values.data[0] = phi = params->startPhase;
   values.data[1] = xi;

   /* fprintf(stdout, "Initial phi=%e, xi=%e\n", values.data[0], values.data[1]); */
   t = 0.0;

   /* Initialize the GSL integrator */
   switch (params->order)
   {
	case LAL_PNORDER_TWO:
	   in4.function = XLALTaylorNDerivatives4PN;
	   break;
	case LAL_PNORDER_TWO_POINT_FIVE:
	   in4.function = XLALTaylorNDerivatives5PN;
	   break;
	case LAL_PNORDER_THREE:
	   in4.function = XLALTaylorNDerivatives6PN;
	   break;
	case LAL_PNORDER_THREE_POINT_FIVE:
	   in4.function = XLALTaylorNDerivatives7PN;
	   break;
	default:
	   snprintf(message, 256, "There are no Et waveforms at order %d\n", params->order);
	   XLALPrintError( message );
	   XLALFree(dummy.data);
	   XLAL_ERROR(XLAL_EINVAL);
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
      XLALFree(dummy.data);
      XLAL_ERROR(XLAL_ENOMEM);
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
         XLALFree(dummy.data);
         XLAL_ERROR(XLAL_EBADLEN);
      }

      h = 4 * m * eta * xi * sin(2.*phi)/1.e14;
      output1->data[ndx] = h;
      if (output2)
      {
        h = 4 * m * eta * xi * cos(2.*phi)/1.e14;
        output2->data[ndx] = h;
      }
      /* fprintf(stdout, "%e %e %e\n", t, h, omega/(m*LAL_PI)); */

      /* Integrate one step forward */
      in4.dydx = &dvalues;
      in4.x = t/m;
      if (XLALRungeKutta4(&newvalues, integrator, funcParams))
      {
         XLALRungeKutta4Free( integrator );
         XLALFree(dummy.data);
         XLAL_ERROR(XLAL_EFUNC);
      }

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
   XLALFree(dummy.data);

   return count;
}
/*---------------------------------------------------------*/
