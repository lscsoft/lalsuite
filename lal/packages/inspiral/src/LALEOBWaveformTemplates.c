/*  <lalVerbatim file="LALEOBWaveformTemplatesCV">
Author: Sathyaprakash, B. S.
$Id$
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{LALEOBWaveformTemplates.c}}

Module to generate two effective-one-body waveforms that differ in
phase by $\pi/2$ (see documentation on \texttt{LALEOBWaveform}
on how these waveforms are generated).

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALEOBWaveformTemplatesCP}
\index{\verb&LALEOBWaveformTemplates()&}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALEOBWaveformTemplatesCV}}

</lalLaTeX>  */
#include <lal/LALInspiral.h>
#include <lal/FindRoot.h>

typedef struct tagrOfOmegaIn {
   REAL8 eta, omega;
} rOfOmegaIn;

void LALHCapDerivatives(REAL8Vector *values, REAL8Vector *dvalues, void *funcParams);
void LALprInit(REAL8 *pr, REAL8 r, InspiralDerivativesIn *ak);
void LALpphiInit(REAL8 *phase, REAL8 r, REAL8 eta);
void LALlightRingRadius(LALStatus *status, REAL8 *x, REAL8 r, void *params);
void LALrOfOmega (LALStatus *status, REAL8 *x, REAL8 r, void *params);

NRCSID (LALEOBWAVEFORMTEMPLATESC, "$Id$");

/*  <lalVerbatim file="LALEOBWaveformTemplatesCP"> */

void LALEOBWaveformTemplates (LALStatus *status,
                     REAL4Vector *signal1,
                     REAL4Vector *signal2,
                     InspiralTemplate *params) 
{ /* </lalVerbatim> */

   INT4 count, nn=4;
   REAL8 amp, eta, m, rn, r, rOld, s, p, q, dt, t, h1, h2, v, omega, f;
   REAL8Vector dummy, values, dvalues, newvalues, yt, dym, dyt;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   void *funcParams;
   expnCoeffs ak;
   expnFunc func;
   rOfOmegaIn rofomegain;
   DFindRootIn rootIn;
   
   INITSTATUS(status, "LALEOBWaveformTemplates", LALEOBWAVEFORMTEMPLATESC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signal1,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signal2,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signal1->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signal2->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->totalMass > 0., status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   ASSERT(ak.totalmass/LAL_MTSUN_SI > 0.4, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(ak.totalmass/LAL_MTSUN_SI < 100, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

/* Allocate all the memory required to dummy and then point the various
   arrays to dummy - this makes it easier to handle memory failures */

   dummy.length = nn * 6;

   values.length = dvalues.length = newvalues.length =
   yt.length = dym.length = dyt.length = nn;

   if (!(dummy.data = (REAL8 * ) LALMalloc(sizeof(REAL8) * nn * 6))) {
      ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
   }

   values.data = &dummy.data[0];
   dvalues.data = &dummy.data[nn];
   newvalues.data = &dummy.data[2*nn];
   yt.data = &dummy.data[3*nn];
   dym.data = &dummy.data[4*nn];
   dyt.data = &dummy.data[5*nn];

   dt = 1./params->tSampling;
   eta = ak.eta;
   m = ak.totalmass;

   t = 0.0;
   in1.t = t;
   in1.t0 = ak.t0;
   in1.v0 = ak.v0;
   in1.vlso = ak.vlso;
   in1.totalmass = ak.totalmass;
   in1.dEnergy = func.dEnergy;
   in1.flux = func.flux;
   in1.coeffs = &ak;

   LALInspiralVelocity(status->statusPtr, &v, &in1);
   CHECKSTATUSPTR(status);

   omega = pow(v,3.);
   f = omega/(LAL_PI*m);

   in2.v0 = ak.v0;
   in2.phi0 = params->startPhase;
   in2.dEnergy = func.dEnergy;
   in2.flux = func.flux;
   in2.coeffs = &ak;
   LALInspiralPhasing1(status->statusPtr, &s, v, &in2);
   CHECKSTATUSPTR(status);
   s = s/2.;

   rofomegain.eta = eta;
   rofomegain.omega = omega;

   rootIn.xacc = 1.0e-8;
   rootIn.function = LALlightRingRadius;
   rootIn.xmax = 2.;
   rootIn.xmin = 4.;
   funcParams = (void *) &rofomegain;

/*-------------------------------------------------------------------
Userful for debugging: Make sure a solution for r exists.
--------------------------------------------------------
   for (r=0.01; r<10.; r+=.01) {
      LALlightRingRadius(status->statusPtr, &x, r, funcParams);
      CHECKSTATUSPTR(status);
      printf("%e %e\n", r, x);
   }
   printf("&\n");
   for (r=1; r<100.; r+=.1) {
      LALrOfOmega(status->statusPtr, &x, r, funcParams);
      CHECKSTATUSPTR(status);
      printf("%e %e\n", r, x);
   }
-------------------------------------------------------------------*/

   LALDBisectionFindRoot(status->statusPtr, &rn, &rootIn, funcParams);
   CHECKSTATUSPTR(status);

   rootIn.function = LALrOfOmega;
   rootIn.xmax = 100.;
   rootIn.xmin = 6.;
   LALDBisectionFindRoot(status->statusPtr, &r, &rootIn, funcParams);
   CHECKSTATUSPTR(status);

   params->rInitial = r;
   params->vInitial = v;
   params->rLightRing = rn;

/* 
   LALInspiralPhasing1(v) gives the GW phase (= twice the orbital phase).
   The ODEs we solve give the orbital phase. Therefore, set the
   initial phase to be half the GW pahse.
*/
   in3.totalmass = ak.totalmass;
   in3.dEnergy = func.dEnergy;
   in3.flux = func.flux;
   in3.coeffs = &ak;
   funcParams = (void *) &in3;

   LALprInit(&p, r, &in3);
   LALpphiInit(&q, r, eta);

   values.data[0] = r;
   values.data[1] = s;
   values.data[2] = p;
   values.data[3] = q;

   
   in4.function = LALHCapDerivatives;
   in4.y = &values;
   in4.h = dt/m;
   in4.n = nn;
   in4.yt = &yt;
   in4.dym = &dym;
   in4.dyt = &dyt;

   t = 0.0;
   count = 0;
   while (count < params->nStartPad) 
   {
      *(signal1->data + count) = *(signal2->data + count) = 0.;
      count++;
   }

   t = 0.0;
   rOld = r+0.1;
   while (r>=rn && r<rOld) {
      ASSERT(count< (INT4)signal1->length, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
      rOld = r;
      v = pow(omega, oneby3);
      amp = params->signalAmplitude *v*v;
      h1 = amp * cos(2.*s);
      h2 = amp * cos(2.*s + LAL_PI_2);
      *(signal1->data + count) = (REAL4) h1;
      *(signal2->data + count) = (REAL4) h2;
      LALHCapDerivatives(&values, &dvalues, funcParams);
      CHECKSTATUSPTR(status);
      omega = dvalues.data[1];
      in4.dydx = &dvalues;
      in4.x = t/m;
      LALRungeKutta4(status->statusPtr, &newvalues, &in4, funcParams);
      CHECKSTATUSPTR(status);

      r = values.data[0] = newvalues.data[0];
      s = values.data[1] = newvalues.data[1];
      p = values.data[2] = newvalues.data[2];
      q = values.data[3] = newvalues.data[3];

      t = (++count-params->nStartPad) * dt;
/*----------------------------------------------------------
      printf("%e %e %e %e %e %e %e\n", t, r, v, s, p, q, h);
      if (v>ak->vlso) printf("TLSO=%e\n", t);
      printf("&\n");
----------------------------------------------------------*/
   }  

/*----------------------------------------------------------------- 
Record the final cutoff frequency of BD Waveforms for record keeping 
-----------------------------------------------------------------*/
   params->rFinal = rOld;
   params->fFinal = pow(v,3.)/(LAL_PI*m);
   params->vFinal = v;
   while (count < (INT4)signal1->length) 
   {
      *(signal1->data + count) = *(signal2->data + count) = 0.;
      count++;
   }

   DETATCHSTATUSPTR(status);
   RETURN(status);
   LALFree(dummy.data);
}

