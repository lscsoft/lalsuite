#include <lal/LALInspiral.h>

void LALUSPA (
   LALStatus status, 
   REAL4Vector *signal,
   InspiralTemplate *params
) 
{
   REAL8 t, pimmc, f, v, df, shft, phi, amp0, amp, psif, psi, sign;
   INT4 n0, nn, nn1, i;
   void *funcParams;
   DIntegrateIn intinp;
   TofVIntegrandIn psiIn;
   expnCoeffs ak;
   expnFunc func;

   INITSTATUS (status, "LALUSPA", LALUSPAC);
   ATTATCHSTATUSPTR(status);
   ASSERT (signal,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   psiIn.dEnergy = func.dEnergy;
   psiIn.flux = func.flux;
   psiIn.coeffs = ak;

   df = params->tSampling/signal->length;
   pimmc = pi * ak.totalmass;
   t = 0.0;
   v = InspiralVelocity(t);
   f = pow(v,3.)/pimmc;

   shft = twopi * (params->nStartPad/params->tSampling + params->startTime);
   phi =  params->startPhase + piby4;
   amp0 = 0.5 * params->signalAmplitude * ak.totalmass * pow(pi/3., 0.5) * 
          params->tSampling; /* * (2. / signal->length); */

   n0 = 2 * ceil (f/df);
   nn = 2 * ceil (ak.fn/df);
   nn1 = 2 * ceil (params->fCutoff/df);
   if (nn1>nn) nn = nn1;

/* 
   All frequency components below f0 and above fn are set to zero  
*/
   i = 0; 
   while (i<n0) *(signal->data + i++) = 0.;
   i = nn; 
   while (i<signal->length) *(signal->data + i++) = 0.;

/* 
   Compute the standard stationary phase approximation. 
*/
   funcParams = (void *) &psiIn;
   intinp.function = LALPsiOfT;
   intinp.type = ClosedInterval;
   for (i=n0; i<nn; i+=2) {
      f = (i/2) * df;
      ak.vf = v = pow(pimmc * f, oneby3);
      if (v==ak.v0) {
         psif = 0.;

      } else {
         if (ak.v0 > v) {
            intinp.xmin = ak.v;
            intinp.xmax = ak.v0;
	    sign = -1.0;
         } else {
            intinp.xmin = ak.v0;
            intinp.xmax = ak.v;
            sign = 1.0;
         }
         LALDRombergIntegrate (status->statusPtr, &psif, &intinp, funcParams);
         psif *= sign;
         CHECKSTATUSPTR(status);
      }
      psi = shft * f - phi + psif;
/* 
      dEnergybyFlux computes 1/(dv/dt) while what we need is 1/(dF/dt):
      dF/dt=(dF/dv)(dv/dt)=-dEnergybyFlux/(dF/dv)=-dEnergybyFlux 3v^2/(pi m^2)
*/
      amp = amp0 * pow(-dEnergybyFlux(v), 0.5) * v;
      h1 = *(signal->data+i) = (REAL4) amp * cos(psi);
      h2 = *(signal->data+i+1) = (REAL4) amp * sin(psi);
/*
      printf ("%e %e \n", v, psif);
      printf ("%e %e %e %e %e\n", f, pow(h1,2.)+pow(h2,2.), h2, psi, psif);
   printf ("&\n");
*/
   }
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

void LALPsiOfT(
   LALStatus *stauts, 
   REAL8 *psioft, 
   REAL8 v, 
   void *param
) {
   REAL8 vf, dE, F;
   TofVIntegrandIn *param

   par = (TofVIntegrandIn *) param;

   dE = par->dEnergy(v, par->coeffs);
   F = par->flux(v, par->coeffs);
   vf = par->coeffs.vf;
   &psioft = 2. * (v*v*v - vf*vf*vf) * dE/F;
}
