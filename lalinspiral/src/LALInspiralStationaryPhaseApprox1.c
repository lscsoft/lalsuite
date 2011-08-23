/*
*  Copyright (C) 2007 B.S. Sathyaprakash, Thomas Cokelaer
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
 * \author B.S. Sathyaprakash
\file
\ingroup LALInspiral_h

\brief This module computes the stationary phase approximation to the
Fourier transform of a chirp waveform by integrating Eq.\eqref{eq_InspiralFourierPhase}.

 * \heading{Prototypes}
 *
 * <tt>LALInspiralStationaryPhaseApprox1()</tt>
 * <ul>
 * <li> \c signalvec: Output containing the inspiral waveform.
 * </li><li> \c params: Input containing binary chirp parameters.
 * </li></ul>
 *
 * \heading{Description}
 *
 * This module generates the Fourier domain waveform that is analogous of
 * the time-domain approximant <tt>TaylorT1.</tt> Instead of re-expanding the
 * the energy and flux functions they are kept in tact and the integral
 * in Eq.\eqref{eq_InspiralFourierPhase} is solved numerically.
 * The code returns the Fourier transform packed in the same way as fftw
 * would for the Fourier transform of a real vector.  For a signal vector
 * of length <tt>n=signalvec->length</tt> (\c n even):
 * <ul>
 * <li> <tt>signalvec->data[0]</tt> is the \e real 0th frequency component of the Fourier transform.
 * </li><li> <tt>signalvec->data[n/2]</tt> is the \e real Nyquist frequency component of the Fourier transform.
 * </li><li> <tt>signalvec->data[k]</tt> and <tt>signalvec->data[n-k],</tt> for <tt>k=1,..., n/2-1,</tt> are
 * the real and imaginary parts of the Fourier transform at a frequency \f$k\Delta f=k/T,\f$ \f$T\f$ being
 * the duration of the signal and \f$\Delta f=1/T\f$ is the frequency resolution.
 * </li></ul>
 *
 * \heading{Algorithm}
 *
 * The lal code \c LALDRombergIntegrate() is used to solve the
 * integral in Eq.\eqref{eq_InspiralFourierPhase}.
 * The reference points are chosen so that on inverse Fourier transforming
 * the time-domain waveform will
 * <ul>
 * <li> be padded with zeroes in the first <tt>params->nStartPad</tt> bins,
 * </li><li> begin with a phase shift of <tt>params->nStartPhase</tt> radians,
 * </li><li> have an amplitude of \f$n v^2.\f$
 * </li></ul>
 *
 * \heading{Uses}
 *
 * \code
   LALInspiralSetup()
   LALInspiralChooseModel()
   LALDRombergIntegrate()
 * \endcode
 * \heading{Notes}
 *
 * If it is required to compare the output of this module with a time domain
 * signal one should use an inverse Fourier transform routine that packs data
 * in the same way as fftw. Moreover, one should divide the resulting inverse
 * Fourier transform by a factor \f$n/2\f$ to be consistent with the
 * amplitude used in time-domain signal models.
 *
*/

#include <lal/LALInspiral.h>
#include <lal/Integrate.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* a local function to compute the phase of the Fourier transform */
void
LALPsiOfT (
   LALStatus *stauts,
   REAL8     *psioft,
   REAL8      v,
   void      *param
   );

NRCSID (LALINSPIRALSTATIONARYPHASEAPPROX1C, "$Id$");

/* This is the main function to compute the stationary phase approximation */


void
LALInspiralStationaryPhaseApprox1 (
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{
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
   ASSERT (signalvec,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
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
   n = signalvec->length;
   nby2 = n/2;
   df = params->tSampling/signalvec->length;
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
   shft = LAL_PI*2.L * (ak.tn + params->nStartPad/params->tSampling + params->startTime);
   phi =  params->startPhase + LAL_PI/4.L;
   amp0 = params->signalAmplitude * ak.totalmass * pow(LAL_PI/12.L, 0.5L) * df;

   signalvec->data[0] = 0.;
   signalvec->data[nby2] = 0.;
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
         signalvec->data[i] = 0.;
         signalvec->data[n-i] = 0.;
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
	 signalvec->data[i] = (REAL4) (amp * cos(psi));
	 signalvec->data[n-i] = -(REAL4) (amp * sin(psi));
      }
   }
   params->fFinal = fn;
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

void
LALPsiOfT(
   LALStatus UNUSED *status,
   REAL8     *psioft,
   REAL8      v,
   void      *param
   )
{
   REAL8 vf, dE, F;
   TofVIntegrandIn *par;

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
