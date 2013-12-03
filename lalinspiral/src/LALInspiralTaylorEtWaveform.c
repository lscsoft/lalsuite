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

/**
 * \author Sathyaprakash, B. S., Cokelaer T.
 * \file
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

typedef struct tagzetaInitIn {
  REAL8 eta, omega, e0;
} zetaInitIn;

static REAL8 XLALzetaInit4PN(
  REAL8      zeta,
  void      *params);

static REAL8 XLALzetaInit5PN(
  REAL8      zeta,
  void      *params);

static REAL8 XLALzetaInit6PN(
  REAL8      zeta,
  void      *params);

static REAL8 XLALzetaInit7PN(
  REAL8      zeta,
  void      *params);


void XLALTaylorEtDerivatives4PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);
void XLALTaylorEtDerivatives5PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);
void XLALTaylorEtDerivatives6PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);
void XLALTaylorEtDerivatives7PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
);

void LALTaylorEtWaveformEngine (
  LALStatus        *status,
  REAL4Vector      *signalvec,
  InspiralTemplate *params,
  InspiralInit     *paramsInit
);

int XLALTaylorEtWaveformEngine (
  REAL4Vector      *signalvec1,
  REAL4Vector      *signalvec2,
  InspiralTemplate *params,
  InspiralInit     *paramsInit
);

static REAL8 XLALzetaInit4PN(
   REAL8      zeta,
   void      *params)
{
   if( !params )
      XLAL_ERROR_REAL8(XLAL_EFAULT);

   zetaInitIn *in;
   REAL8 x, zeta2, zeta32, eta, eta2;

   in = (zetaInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   zeta2 = zeta*zeta;
   zeta32 = pow(zeta, 1.5);

   x = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2 );
   x -= in->omega;

   return x;
}

static REAL8 XLALzetaInit5PN(
   REAL8      zeta,
   void      *params)
{
   if( !params )
      XLAL_ERROR(XLAL_EFAULT);

   zetaInitIn *in;
   REAL8 x, zeta2, zeta32, eta, eta2;

   in = (zetaInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   zeta2 = zeta*zeta;
   zeta32 = pow(zeta, 1.5);

   x = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2 );
   x -= in->omega;

   return x;
}

static REAL8 XLALzetaInit6PN(
   REAL8      zeta,
   void      *params)
{
   if( !params )
      XLAL_ERROR(XLAL_EFAULT);

   zetaInitIn *in;
   REAL8 x, zeta2, zeta3, zeta32, eta, eta2, eta3, pisq;

   in = (zetaInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   eta3 = eta2*eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta32 = pow(zeta, 1.5);
   pisq = LAL_PI*LAL_PI;

   x = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2
     + (41445./1024. - (309715./3072. - 205./64.*pisq) * eta
     + 1215./1024.*eta2 + 45./1024*eta3) * zeta3);
   x -= in->omega;

   return x;
}

static REAL8 XLALzetaInit7PN(
   REAL8      zeta,
   void      *params)
{
   if( !params )
      XLAL_ERROR(XLAL_EFAULT);

   zetaInitIn *in;
   REAL8 x, zeta2, zeta3, zeta32, eta, eta2, eta3, pisq;

   in = (zetaInitIn *) params;
   eta = in->eta;
   eta2 = eta*eta;
   eta3 = eta2*eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta32 = pow(zeta, 1.5);
   pisq = LAL_PI*LAL_PI;

   x = zeta32 * ( 1. + 0.125 * (9.+eta)*zeta
     + (891./128. - 201./64.*eta + 11/128.*eta2) * zeta2
     + (41445./1024. - (309715./3072. - 205./64.*pisq) * eta
     + 1215./1024.*eta2 + 45./1024*eta3) * zeta3);
   x -= in->omega;

   return x;
}

void XLALTaylorEtDerivatives4PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /*if( !values || !dvalues || !funcParams )
      XLAL_ERROR_NULL(XLAL_EFAULT);*/

   InspiralDerivativesIn *ak;
   REAL8 zeta, zeta2, zeta3, zeta5, zeta32, eta, eta2, fourpi;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   zeta = values->data[1];
   eta2 = eta * eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta5 = zeta2*zeta3;
   zeta32 = pow(zeta, 1.5);
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

void XLALTaylorEtDerivatives5PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /*if( !values || !dvalues || !funcParams )
      XLAL_ERROR_NULL(XLAL_EFAULT);*/

   InspiralDerivativesIn *ak;
   REAL8 zeta, zeta2, zeta3, zeta5, zeta32, zeta52, eta, eta2, fourpi;

   ak = (InspiralDerivativesIn *) funcParams;
   eta = ak->coeffs->eta;

   zeta = values->data[1];
   eta2 = eta * eta;
   zeta2 = zeta*zeta;
   zeta3 = zeta2*zeta;
   zeta5 = zeta2*zeta3;
   zeta32 = pow(zeta, 1.5);
   zeta52 = pow(zeta, 2.5);
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

void XLALTaylorEtDerivatives6PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /*if( !values || !dvalues || !funcParams )
      XLAL_ERROR_NULL(XLAL_EFAULT);*/

   InspiralDerivativesIn *ak;
   REAL8 zeta, zeta2, zeta3, zeta5, zeta32, zeta52, eta, eta2, eta3, pisq, fourpi;

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

void XLALTaylorEtDerivatives7PN(
  REAL8Vector *values,
  REAL8Vector *dvalues,
  void        *funcParams
)
{
   /*if( !values || !dvalues || !funcParams )
      XLAL_ERROR_NULL(XLAL_EFAULT);*/

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


void LALTaylorEtWaveform (
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{
   XLAL_PRINT_DEPRECATION_WARNING("XLALTaylorEtWaveform");
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   if( XLALTaylorEtWaveform(signalvec, params) )
      ABORTXLAL(status);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

int XLALTaylorEtWaveform (
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{

   InspiralInit paramsInit;

   /* Check the relevant pointers */
   if( !signalvec || !(signalvec->data) || !params )
      XLAL_ERROR(XLAL_EFAULT);

   /* Check the parameters are sane */
   if( params->nStartPad < 0 || params->nEndPad < 0 || params->fLower <= 0 
         || params->tSampling <= 0 || params->totalMass <= 0. )
      XLAL_ERROR(XLAL_EINVAL);

   if( XLALInspiralSetup(&(paramsInit.ak), params) )
      XLAL_ERROR(XLAL_EFUNC);

   if( XLALInspiralChooseModel(&(paramsInit.func), &(paramsInit.ak), params) )
      XLAL_ERROR(XLAL_EFUNC);

   memset(signalvec->data, 0, signalvec->length * sizeof( REAL4 ));

   /* Call the engine function */
   if( XLALTaylorEtWaveformEngine(signalvec, NULL, params, &paramsInit) )
      XLAL_ERROR(XLAL_EFUNC);

   return XLAL_SUCCESS;
}

int XLALTaylorEtWaveformTemplates (
   REAL4Vector      *signalvec1,
   REAL4Vector      *signalvec2,
   InspiralTemplate *params
   )
{

   InspiralInit paramsInit;

   /* Check the relevant pointers */
   if( !signalvec1 || !(signalvec1->data) || !signalvec2 || !(signalvec2->data) || !params )
      XLAL_ERROR(XLAL_EFAULT);

   /* Check the parameters are sane */
   if( params->nStartPad < 0 || params->nEndPad < 0 || params->fLower <= 0 
         || params->tSampling <= 0 || params->totalMass <= 0. )
      XLAL_ERROR(XLAL_EINVAL);

   if( XLALInspiralSetup(&(paramsInit.ak), params) )
      XLAL_ERROR(XLAL_EFUNC);

   if( XLALInspiralChooseModel(&(paramsInit.func), &(paramsInit.ak), params) )
      XLAL_ERROR(XLAL_EFUNC);

   memset(signalvec1->data, 0, signalvec1->length * sizeof( REAL4 ));
   memset(signalvec2->data, 0, signalvec2->length * sizeof( REAL4 ));

   /* Call the engine function */
   if( XLALTaylorEtWaveformEngine(signalvec1, signalvec2, params, &paramsInit) )
      XLAL_ERROR(XLAL_EFUNC);

   return XLAL_SUCCESS;
}

/*---------------- Main engine function --------------*/
/* LAL wrapper of XLAL engine function */
void LALTaylorEtWaveformEngine (
                LALStatus        *status,
                REAL4Vector      *signalvec,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                )
{
   XLAL_PRINT_DEPRECATION_WARNING("XLALTaylorEtWaveformEngine");
   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   if( XLALTaylorEtWaveformEngine(signalvec, NULL, params, paramsInit) )
      ABORTXLAL(status);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/* Actual engine */
int XLALTaylorEtWaveformEngine (
                REAL4Vector      *signalvec1,
                REAL4Vector      *signalvec2,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                )
{

   void                  *funcParams;
   UINT4                 length=0, count, ndx;
   INT4                  nn=2;
   REAL8                 xmin, xmax, xacc;
   REAL8                 h, omega, omegaMax, t, dt, m, eta, phi, zeta, v;
   REAL8Vector           dummy, values, dvalues, newvalues, yt, dym, dyt;
   rk4GSLIntegrator      *integrator = NULL;
   InspiralDerivativesIn in2;
   zetaInitIn            in3;
   rk4In                 in4;
   expnCoeffs            ak;
   expnFunc              func;
   REAL8 (*rootFunction)(REAL8, void *);


/* Allocate all the memory required to dummy and then point the various
   arrays to dummy - this makes it easier to handle memory failures */

   dummy.length = nn * 6;

   if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * nn * 6)))
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
   length = signalvec1->length;
   eta = ak.eta;
   m = ak.totalmass;

   /* Begin initial conditions */
   /* Given omega compute zeta by solving Eq.(18b) of TG08 */
   omega = LAL_PI * params->fLower * m;
   v = pow(omega,1./3.);

   xmin = v*v/2.;
   xmax = 1.;
   xacc = 1.0e-16;

   /* Initialize the GSL integrator */
   switch (params->order)
   {
	case LAL_PNORDER_TWO:
	   rootFunction = &XLALzetaInit4PN;
	   break;
	case LAL_PNORDER_TWO_POINT_FIVE:
	   rootFunction = &XLALzetaInit5PN;
	   break;
	case LAL_PNORDER_THREE:
	   rootFunction = &XLALzetaInit6PN;
	   break;
	case LAL_PNORDER_THREE_POINT_FIVE:
	   rootFunction = &XLALzetaInit7PN;
	   break;
	default:
           XLALPrintError("XLAL Error: %s - There are no Et waveforms at order %d\n", __func__, params->order);
           XLAL_ERROR(XLAL_EINVAL);
	   break;
   }
   in3.eta = ak.eta;
   in3.omega = omega;
   funcParams = (void *) &in3;


   zeta = XLALDBisectionFindRoot(rootFunction, xmin, xmax, xacc, funcParams);
   if (XLAL_IS_REAL8_FAIL_NAN(zeta))
      XLAL_ERROR(XLAL_EFUNC);
   /* End of initial conditions */

   values.data[0] = phi = params->startPhase;
   values.data[1] = zeta;

   /* fprintf(stdout, "Initial phi=%e, zeta=%e\n", values.data[0], values.data[1]); */
   t = 0.0;

   /* Initialize the GSL integrator */
   switch (params->order)
   {
	case LAL_PNORDER_TWO:
	   in4.function = XLALTaylorEtDerivatives4PN;
	   break;
	case LAL_PNORDER_TWO_POINT_FIVE:
	   in4.function = XLALTaylorEtDerivatives5PN;
	   break;
	case LAL_PNORDER_THREE:
	   in4.function = XLALTaylorEtDerivatives6PN;
	   break;
case LAL_PNORDER_THREE_POINT_FIVE:
	   in4.function = XLALTaylorEtDerivatives7PN;
	   break;
	default:
           XLALPrintError("XLAL Error: %s - There are no Et waveforms at order %d\n", __func__, params->order);
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

   /* Initialize the integrator */
   if (!(integrator = XLALRungeKutta4Init(nn, &in4)))
      XLAL_ERROR(XLAL_EFUNC);

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
	XLALFree(dummy.data);
        XLAL_ERROR(XLAL_EBADLEN);
      }

      h = 4 * m * eta * zeta * sin(2.*phi)/1.e14;
      signalvec1->data[ndx] = h;
      if (signalvec2)
         signalvec2->data[ndx] = 4 * m * eta * zeta * cos(2.*phi)/1.e14;
      /* fprintf(stdout, "%e %e %e\n", t, h, omega/(m*LAL_PI)); */

      /* Integrate one step forward */
      in4.dydx = &dvalues;
      in4.x = t/m;
      if( XLALRungeKutta4(&newvalues, integrator, funcParams) )
      {
         XLALRungeKutta4Free( integrator );
         XLALFree(dummy.data);
         XLAL_ERROR(XLAL_EFUNC);
      }

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
   XLALFree(dummy.data);

   return XLAL_SUCCESS;
}
/*---------------------------------------------------------*/
