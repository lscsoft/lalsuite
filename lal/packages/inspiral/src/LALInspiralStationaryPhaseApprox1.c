/**** <lalVerbatim file="LALInspiralStationaryPhaseApprox1CV">
 * Author: B.S. Sathyaprakash
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{LALInspiralStationaryPhaseApprox1.c}}
 * This module computes the usual stationary phase approximation to the
 * Fourier transform of a chirp waveform.
 * %% A one-line description of the function(s) defined in this module.
 *
 * \subsubsection*{Prototypes}
 * \input{LALInspiralStationaryPhaseApprox1CP}
 * \idx{LALInspiralStationaryPhaseApprox1()}
 *
 * \subsubsection*{Description}
 *
 * %% A description of the data analysis task performed by this function;
 * %% this is the main place to document the module.
 *
 * \subsubsection*{Algorithm}
 *
 * %% A description of the method used to perform the calculation.
 *
 * \subsubsection*{Uses}
 *
 * %% List of any external functions called by this function.
 * \begin{verbatim}
 * None
 * \end{verbatim}
 * \subsubsection*{Notes}
 *
 * %% Any relevant notes.
 *
 * \vfill{\footnotesize\input{LALInspiralStationaryPhaseApprox1CV}}
 *
 **** </lalLaTeX> */

#include <lal/LALInspiral.h>
#include <lal/Integrate.h>

/* a local function to compute the phase of the Fourier transform */
void LALPsiOfT(
   LALStatus *stauts, 
   REAL8 *psioft, 
   REAL8 v, 
   void *param
);

NRCSID (LALINSPIRALSTATIONARYPHASEAPPROX1C, "$Id$");

/* This is the main function to compute the stationary phase approximation */

/*  <lalVerbatim file="LALInspiralStationaryPhaseApprox1CP"> */
void LALInspiralStationaryPhaseApprox1 (
   LALStatus *status, 
   REAL4Vector *signal,
   InspiralTemplate *params) 
{ /* </lalVerbatim>  */
   REAL8 t, pimmc, f0, fn, f, v, df, shft, phi, amp0, amp, psif, psi, sign;
   INT4 n, i, nby2;
   void *funcParams;
   DIntegrateIn intinp;
   TofVIntegrandIn psiIn;
   TofVIn tofvin;
   expnCoeffs ak;
   expnFunc func;

   INITSTATUS (status, "LALInspiralStationaryPhaseApprox1", LALINSPIRALSTATIONARYPHASEAPPROX1C);
   ATTATCHSTATUSPTR(status);
   ASSERT (signal,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signal->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

/* Set up the coefficients in post-Newtonian expansion, vlso, etc. */
   LALInspiralSetup (status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
/* Set up the functions required for the chosen signal approximation sheme */
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

/* Cast the struct required for the phase function */
   psiIn.dEnergy = func.dEnergy;
   psiIn.flux = func.flux;
   psiIn.coeffs = &ak;

/* Cast the struct required for computing initial conditions */
   n = signal->length;
   nby2 = n/2;
   df = params->tSampling/signal->length;
   pimmc = LAL_PI * ak.totalmass;
   t = 0.0;
   tofvin.t = t;
   tofvin.v0 = ak.v0;
   tofvin.t0 = ak.t0;
   tofvin.vlso = ak.vlsoT2;
   tofvin.totalmass = ak.totalmass;
   tofvin.dEnergy = func.dEnergy;
   tofvin.flux = func.flux;
   tofvin.coeffs = &ak;

/* Compute the initial velocity frequency  */
   LALInspiralVelocity(status->statusPtr, &v, &tofvin);
   f0 = pow(v,3.L)/pimmc;
   fn = (params->fCutoff < ak.fn) ? params->fCutoff : ak.fn; 

   /* If we want to pad with zeroes in the beginning then the instant of
    * coalescence will be the chirp time + the duration for which padding
    * is needed. Thus, in the equation below nStartPad occurs with a +ve sign.
    */
   shft = LAL_PI*2.L * (params->nStartPad/params->tSampling + params->startTime);
   phi =  params->startPhase - LAL_PI/4.L;
   amp0 = params->signalAmplitude * ak.totalmass * pow(LAL_PI/12.L, 0.5L) * 
          params->tSampling / (REAL8) signal->length; 

   signal->data[0] = 0.;
   signal->data[nby2] = 0.;
/* 
   Compute the standard stationary phase approximation. 
*/
   funcParams = (void *) &psiIn;
   intinp.function = LALPsiOfT;
   intinp.type = ClosedInterval;
   for (i=1; i<nby2; i++) 
   {
      f = i * df;
      if (f<f0 || f>fn)
      {
/* 
   All frequency components below f0 and above fn are set to zero  
*/
         signal->data[i] = 0.;
         signal->data[n-i] = 0.;
      }
      else
      {
         ak.vf = v = pow(pimmc * f, oneby3);
         if (v==ak.v0) 
         {
            psif = 0.;

         } 
	 else 
         {
            if (ak.v0 > ak.vf) 
	    {
               intinp.xmin = ak.vf;
               intinp.xmax = ak.v0;
	       sign = -1.0;
            } 
	    else 
	    {
               intinp.xmin = ak.v0;
               intinp.xmax = ak.vf;
               sign = 1.0;
            }
	    LALDRombergIntegrate (status->statusPtr, &psif, &intinp, funcParams);
	    CHECKSTATUSPTR(status);
	    psif *= sign;
	 }
	 psi = shft * f + phi + psif;
/* 
      dEnergybyFlux computes 1/(dv/dt) while what we need is 1/(dF/dt):
      dF/dt=(dF/dv)(dv/dt)=-dEnergybyFlux/(dF/dv)=-dEnergybyFlux 3v^2/(LAL_PI m^2)
*/
	 amp = amp0 * pow(-func.dEnergy(v,&ak)/func.flux(v,&ak),0.5) * v;
	 signal->data[i] = (REAL4) (amp * cos(psi));
	 signal->data[n-i] = -(REAL4) (amp * sin(psi));
      }
   }
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

void LALPsiOfT(
   LALStatus *status, 
   REAL8 *psioft, 
   REAL8 v, 
   void *param
) {
   REAL8 vf, dE, F;
   TofVIntegrandIn *par;

   status = NULL;
   par = (TofVIntegrandIn *) param;

   /* The integrand below has an overall -ve sign as compared to
    * Eq. (3.5) of DIS3; this is because as oppsed to Eq. (3.5) of
    * DIS3 here we integrate from v0 to vf instead of from vf to v0.
    */
   dE = par->dEnergy(v, par->coeffs);
   F = par->flux(v, par->coeffs);
   vf = par->coeffs->vf;
   *psioft = 2. * (v*v*v - vf*vf*vf) * dE/F;
}
