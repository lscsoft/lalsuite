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

NRCSID (LALEOBWAVEFORMC, "$Id$");

void LALEOBWaveform (LALStatus *status,
                     REAL4Vector *signal,
                     InspiralTemplate *params) 
{ 
   INT4 count, i, nn=4;
   REAL8 eta, m, rn, r, s, p, q, dt, t, h, v, omega, f, x;
   REAL8Vector values, dvalues, newvalues, yt, dym, dyt;
   TofVIn in1;
   InspiralPhaseIn in2;
   InspiralDerivativesIn in3;
   rk4In in4;
   void *funcParams;
   expnCoeffs ak;
   expnFunc func;
   rOfOmegaIn rofomegain;
   DFindRootIn rootIn;
   
   INITSTATUS(status, "LALEOBWaveform", LALEOBWAVEFORMC);
   ATTATCHSTATUSPTR(status);

   ASSERT(signal,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT(params->nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->fLower > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->totalMass > 0.4, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(params->totalMass < 100, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   ASSERT(ak.totalmass/LAL_MTSUN_SI > 0.4, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT(ak.totalmass/LAL_MTSUN_SI < 100, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   values.data = (REAL8 *)LALMalloc(sizeof(REAL8)*nn);
   ASSERT (values.data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   dvalues.data = (REAL8 *)LALMalloc(sizeof(REAL8)*nn);
   ASSERT (dvalues.data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   newvalues.data = (REAL8 *)LALMalloc(sizeof(REAL8)*nn);
   ASSERT (newvalues.data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

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
   LALInspiralPhase(status->statusPtr, &s, v, &in2);
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
   LALInspiralPhase(v) gives the GW phase (= twice the orbital phase).
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

   dym.data = (REAL8 *)LALMalloc(sizeof(REAL8)*nn);
   ASSERT (dym.data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   dyt.data = (REAL8 *)LALMalloc(sizeof(REAL8)*nn);
   ASSERT (dyt.data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   yt.data = (REAL8 *)LALMalloc(sizeof(REAL8)*nn);
   ASSERT (yt.data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   
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
      *(signal->data + count++) = 0.;

   t = 0.0;
   while (r>=rn) {
      ASSERT(count< signal->length, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
      v = pow(omega, oneby3);
      h = params->signalAmplitude * v*v * cos(2.*s);
      *(signal->data + count) = (REAL4) h;
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
   params->rFinal = r;
   params->vFinal = v;
   f = pow(v,3.)/(LAL_PI*m);
   while (count < signal->length) *(signal->data + count++) = 0.;

   LALFree(yt.data);
   LALFree(dyt.data);
   LALFree(dym.data);
   LALFree(values.data);
   LALFree(dvalues.data);
   LALFree(newvalues.data);
}

/*----------------------------------------------------------------------*/

void LALprInit(
   REAL8 *pr, 
   REAL8 r, 
   InspiralDerivativesIn *ak) 
{

/* 
This combines Eq. (4.13), (4.14) (4.16) of BD2 to get the initial value of pr 
*/

   REAL8 eta, z, omega, jadiab, djadiabdr, H0cap, v, FDIS; 
   REAL8 A, r2, r3, cr;

   eta = ak->coeffs->eta;
   r2 = r*r;
   r3 = r2*r;

   A = 1. - 2./r + 2.*eta/r3;
   z = r3 * A * A/(r3 - 3.*r2 + 5.*eta);
   omega = sqrt((1.-3*eta/r2)/(r3 * (1. + 2*eta*(sqrt(z)-1.))));
   jadiab = sqrt (r2 * (r2 - 3.*eta)/(r3 - 3.*r2 + 5.*eta));
   djadiabdr = (r3*r3 - 6.*r3*r2 + 3.*eta*r2*r2 +20.*eta*r3-30.*eta*eta*r)/
               (pow(r3 - 3.*r2 + 5.*eta, 2.)*2.*jadiab);
   H0cap = sqrt(1. + 2.*eta*(-1. + sqrt(z)))/eta;
   cr = A*A/((1.-6.*eta/r2) * eta * H0cap * sqrt(z));
   v = pow(omega,oneby3);
   FDIS = -ak->flux(v, ak->coeffs)/(eta * v*v*v); 
   *pr = FDIS/(djadiabdr*cr);
}

/*----------------------------------------------------------------------*/
void LALpphiInit(
   REAL8 *phase, 
   REAL8 r, 
   REAL8 eta)
{
   REAL8 r2, r3;
   r2 = r*r;
   r3 = r2*r;
   *phase = pow(r2 * (r2 - 3.*eta) / (r3 - 3.* r2 + 5.*eta), 0.5);
}

/*----------------------------------------------------------------------*/
NRCSID (LALROFOMEGAC, "$Id$");
void LALrOfOmega (
   LALStatus *status, 
   REAL8 *x, 
   REAL8 r, 
   void *params) 
{
   REAL8 r2, r3, A, z, eta, omega;
   rOfOmegaIn *rofomegain;

   INITSTATUS(status, "LALrOfOmega", LALROFOMEGAC);
   rofomegain = (rOfOmegaIn *) params;
   eta = rofomegain->eta;
   omega = rofomegain->omega;
   r2 = r*r;
   r3 = r2*r;

   A = 1. - 2./r + 2.*eta/r3;
   z = r3 * A * A/(r3 - 3.*r2 + 5.*eta);
   *x = omega - sqrt((1.-3*eta/r2)/(r3 * (1. + 2*eta*(sqrt(z)-1.))));
   RETURN(status);
}

/*----------------------------------------------------------------------*/
NRCSID (LALLIGHTRINGRADIUSC, "$Id$");
void LALlightRingRadius(
   LALStatus *status, 
   REAL8 *x, 
   REAL8 r, 
   void *params) 
{
   REAL8 eta;
   rOfOmegaIn *rofomegain;

   INITSTATUS(status, "LALlightRingRadius", LALLIGHTRINGRADIUSC);
   rofomegain = (rOfOmegaIn *) params;
   eta = rofomegain->eta;
   *x = pow(r,3.) - 3.*r*r + 5.* eta;
   RETURN(status);
}


void LALHCapDerivatives(
   REAL8Vector *values, 
   REAL8Vector *dvalues, 
   void *funcParams) 
{
   REAL8 r, s, p, q, r2, r3, p2, q2, A, B, dA, dB, hcap, Hcap, etahH;
   REAL8 omega, v, eta;
   InspiralDerivativesIn *ak;

   ak = (InspiralDerivativesIn *) funcParams;

   eta = ak->coeffs->eta;

   r = values->data[0];
   s = values->data[1];
   p = values->data[2];
   q = values->data[3];

   r2 = r*r;
   r3 = r2*r;
   p2 = p*p;
   q2 = q*q;
   A = 1. - 2./r + 2.*eta/r3;
   B = (1. - 6.*eta/r2)/A;
   dA = 2./r2 - 6.*eta/pow(r,4.);
   dB = (-dA * B + 12.*eta/r3)/A;
   hcap = pow (A*(1. + p2/B + q*q/r2), 0.5);
   Hcap = pow (1. + 2.*eta*(hcap - 1.), 0.5) / eta;
   etahH = eta*hcap*Hcap;

   dvalues->data[0] = A*A*p/((1. - 6.*eta/r2) * etahH);
   dvalues->data[1] = omega = A * q / (r2 * etahH);

   v = pow(omega,oneby3);

   dvalues->data[2] = -0.5 * (dA * hcap * hcap/A - p2 * A * dB/(B*B) - 2. * A * q2/r3) / etahH;
   dvalues->data[3] = -ak->flux(v, ak->coeffs)/(eta * v*v*v); 
}
