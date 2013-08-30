/*
*  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Craig Robinson, Thomas Cokelaer, Drew Keppel
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
 * \author Sathyaprakash, B. S.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief This module generates the inspiral waveform from a binary consisting of
 * two spinning compact stars.
 *
 * ### Prototypes ###
 *
 * <tt>LALInspiralSpinningBHBinary()</tt>
 * <ul>
 * <li> \c signalvec: Output containing the spin modulated inspiral waveform.</li>
 * <li> \c in: Input containing binary chirp parameters.</li>
 * </ul>
 *
 * ### Description ###
 *
 * Using the formalism described in Apostolatos
 * et al [\ref ACST94] and Blanchet et al. [\ref BDIWW1995] and formulas
 * summarized in Sec.\ \ref sec_smirches this module computes
 * the spin-modulated chirps from a pair of compact stars in orbit around
 * each other.
 *
 * ### Algorithm ###
 *
 * This code uses a fourth-order Runge-Kutta algorithm to solve the nine
 * first-order, coupled, ordinary differential equations in \eqref{eqn_precession1}
 * \eqref{eqn_precession2} and \eqref{eqn_precession3}. The solution is then used
 * in \eqref{eqn_waveform} (and following equations) to get the waveform  emitted
 * by a spinning black hole binary.
 *
 * ### Uses ###
 *
 * \code
 * LALInspiralSetup()
 * LALInspiralChooseModel()
 * LALInspiralVelocity()
 * LALInspiralPhasing3()
 * LALRungeKutta4()
 * \endcode
 *
 */

/*
   Interface routine needed to generate time-domain T- or a P-approximant
   waveforms by solving the ODEs using a 4th order Runge-Kutta; April 5, 00.
*/
#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>


/* computes the polarisation angle psi */
static void LALInspiralPolarisationAngle( REAL8 *psi, REAL8 psiOld, REAL8 *NCap, REAL8 *L);
/* compute beam factors Fplus and Fcross */
static void LALInspiralBeamFactors(REAL8 *Fplus, REAL8 *Fcross, REAL8 Fp0, REAL8 Fp1, REAL8 Fc0, REAL8 Fc1, REAL8 psi);
/* compute polarisation phase phi */
static void LALInspiralPolarisationPhase(REAL8 *phi, REAL8 phiOld, REAL8 NCapDotLByL, REAL8 Fplus, REAL8 Fcross);
/* computes modulated amplitude A */
static void LALInspiralModulatedAmplitude(REAL8 *amp, REAL8 v, REAL8 amp0, REAL8 NCapDotLByL, REAL8 Fplus, REAL8 Fcross);
/* computes the derivatives of vectors L, S1 and S2 */
void LALACSTDerivatives (REAL8Vector *values, REAL8Vector *dvalues, void *funcParams);

/* computes carrier phase Phi */
/* func.phasing3(LALStatus *status, REAL8 *Phi, REAL8 Theta, REAL8 *ak); */

/* computes precession correction to phase dThi */

/* Routine to generate inspiral waveforms from binaries consisting of spinning objects */

void
LALInspiralSpinModulatedWave(
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *in
   )
{


	UINT4 nDEDim=9/*, Dim=3*/;
        UINT4 j, count;
	/* polarisation and antenna pattern */
	REAL8 psi, psiOld, Fplus, Fcross, Fp0, Fp1, Fc0, Fc1, magL, amp, amp0, phi, phiOld, /*Amplitude,*/ Phi/*, dPhi*/;
	REAL8 v, t, tMax, dt, f, fOld, fn, phase, phi0, MCube, Theta, etaBy5M, NCapDotL, NCapDotLByL;
	/* source direction angles */
	InspiralACSTParams acstPars;
	REAL8Vector dummy, values, dvalues, newvalues, yt, dym, dyt;
	void *funcParams;
	rk4In rk4in;
        rk4GSLIntegrator *integrator;
	expnFunc func;
	expnCoeffs ak;


	INITSTATUS(status);
	ATTATCHSTATUSPTR(status);
	ASSERT(signalvec,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(signalvec->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
	ASSERT(in,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);

       	LALInspiralSetup (status->statusPtr, &ak, in);
	CHECKSTATUSPTR(status);
	LALInspiralChooseModel(status->statusPtr, &func, &ak, in);
	CHECKSTATUSPTR(status);

	/* Allocate space for all the vectors in one go ... */

	dummy.length = nDEDim * 6;
	values.length = dvalues.length = newvalues.length = yt.length = dym.length = dyt.length = nDEDim;

	if (!(dummy.data = (REAL8 *) LALMalloc(sizeof(REAL8) * dummy.length)))
	{
		ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
	}

	/* and distribute as necessary */
	values.data = &dummy.data[0];             /* values of L, S1, S2 at current time */
	dvalues.data = &dummy.data[nDEDim];       /* values of dL/dt, dS1/dt, dS2/dt at current time */
	newvalues.data = &dummy.data[2*nDEDim];   /* new values of L, S1, S2 after time increment by dt*/
	yt.data = &dummy.data[3*nDEDim];          /* Vector space needed by LALRungeKutta4 */
	dym.data = &dummy.data[4*nDEDim];         /* Vector space needed by LALRungeKutta4 */
	dyt.data = &dummy.data[5*nDEDim];         /* Vector space needed by LALRungeKutta4 */

	v = pow(LAL_PI * in->totalMass * LAL_MTSUN_SI * in->fLower, 1.L/3.L);
	MCube = pow(LAL_MTSUN_SI * in->totalMass, 3.L);
	/* Fill in the structure needed by the routine that computes the derivatives */
	/* constant spins of the two bodies */
	acstPars.magS1 = in->spin1[0]*in->spin1[0] + in->spin1[1]*in->spin1[1] + in->spin1[2]*in->spin1[2];
	acstPars.magS2 = in->spin2[0]*in->spin2[0] + in->spin2[1]*in->spin2[1] + in->spin2[2]*in->spin2[2];
	/* Direction to the source in the solar system barycenter */
	acstPars.NCap[0] = sin(in->sourceTheta)*cos(in->sourcePhi);
	acstPars.NCap[1] = sin(in->sourceTheta)*sin(in->sourcePhi);
	acstPars.NCap[2] = cos(in->sourceTheta);
	/* total mass in seconds */
	acstPars.M = LAL_MTSUN_SI * in->totalMass;
	/* Combination of masses that appears in the evolution of L, S1 and S2 */
	acstPars.fourM1Plus = (4.L*in->mass1 + 3.L*in->mass2)/(2.*in->mass1*MCube);
	acstPars.fourM2Plus = (4.L*in->mass2 + 3.L*in->mass1)/(2.*in->mass2*MCube);
	acstPars.oneBy2Mcube = 0.5L/MCube;
	acstPars.threeBy2Mcube = 1.5L/MCube;
	acstPars.thirtytwoBy5etc = (32.L/5.L) * pow(in->eta,2.) * acstPars.M;

	/* Constant amplitude of GW */
	amp0 = 2.L*in->eta * acstPars.M/in->distance;
	/* Initial magnitude of the angular momentum, determined by the initial frequency/velocity */
	magL = in->eta * (acstPars.M * acstPars.M)/v;
	/* Initial angular momentum and spin vectors */
	/* Angular momentum */
	values.data[0] = magL * sin(in->orbitTheta0) * cos(in->orbitPhi0);
	values.data[1] = magL * sin(in->orbitTheta0) * sin(in->orbitPhi0);
	values.data[2] = magL * cos(in->orbitTheta0);
	/* Spin of primary */
	values.data[3] = in->spin1[0];
	values.data[4] = in->spin1[1];
	values.data[5] = in->spin1[2];
	/* Spin of secondarfy */
	values.data[6] = in->spin2[0];
	values.data[7] = in->spin2[1];
	values.data[8] = in->spin2[2];
	/* Structure needed by RungeKutta4 for the evolution of the differential equations */
	rk4in.function = LALACSTDerivatives;
	rk4in.y = &values;
	rk4in.h = 1.L/in->tSampling;
	rk4in.n = nDEDim;
	rk4in.yt = &yt;
	rk4in.dym = &dym;
	rk4in.dyt = &dyt;

        xlalErrno = 0;
        /* Initialize GSL integrator */
        if (!(integrator = XLALRungeKutta4Init(nDEDim, &rk4in)))
        {
          INT4 errNum = XLALClearErrno();
          LALFree(dummy.data);

          if (errNum == XLAL_ENOMEM)
            ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
          else
            ABORTXLAL( status );
        }

	/* Pad the first nStartPad elements of the signal array with zero */
	count = 0;
	while ((INT4)count < in->nStartPad)
	{
		signalvec->data[count] = 0.L;
		count++;
	}

	/* Get the initial conditions before the evolution */
	t = 0.L;                                                 /* initial time */
	dt = 1.L/in->tSampling;                                  /* Sampling interval */
	etaBy5M = in->eta/(5.L*acstPars.M);                      /* Constant that appears in ... */
	Theta = etaBy5M * (in->tC - t);                          /* ... Theta */
	tMax = in->tC - dt;                                      /* Maximum duration of the signal */
	fn = (ak.flso > in->fCutoff) ? ak.flso : in->fCutoff;    /* Frequency cutoff, smaller of user given f or flso */

	Phi = func.phasing3(Theta, &ak);                         /* Carrier phase at the initial time */
	if (XLAL_IS_REAL8_FAIL_NAN(Phi))
		ABORTXLAL(status);
	f = func.frequency3(Theta, &ak);                         /* Carrier Freqeuncy at the initial time */
	if (XLAL_IS_REAL8_FAIL_NAN(f))
		ABORTXLAL(status);

	/* Constants that appear in the antenna pattern */
	Fp0 = (1.L + cos(in->sourceTheta)*cos(in->sourceTheta)) * cos(2.L*in->sourcePhi);
	Fp1 = cos(in->sourceTheta) * sin(2.L*in->sourcePhi);
	Fc0 = Fp0;
	Fc1 = -Fp1;

	psiOld = 0.L;
	phiOld = 0.L;
	magL = sqrt(values.data[0]*values.data[0] + values.data[1]*values.data[1] + values.data[2]*values.data[2]);
	NCapDotL = acstPars.NCap[0]*values.data[0] + acstPars.NCap[1]*values.data[1] + acstPars.NCap[2]*values.data[2];
	NCapDotLByL = NCapDotL/magL;

	/* Initial values of */
	/* the polarization angle */
	LALInspiralPolarisationAngle(&psi, psiOld, &acstPars.NCap[0], &values.data[0]);
	/* beam pattern factors */
	LALInspiralBeamFactors(&Fplus, &Fcross, Fp0, Fp1, Fc0, Fc1, psi);
	/* the polarization phase */
	LALInspiralPolarisationPhase(&phi, phiOld, NCapDotLByL, Fplus, Fcross);
	/* modulated amplitude */
	LALInspiralModulatedAmplitude(&amp, v, amp0, NCapDotLByL, Fplus, Fcross);

	phase = Phi + phi;
	phi0 = -phase + in->startPhase;

	fOld = f - 0.1;
	while (f < fn && t < tMax && f>fOld)
	{
		ASSERT((INT4)count < (INT4)signalvec->length, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
		/* Subtract the constant initial phase (chosen to have a phase of in->startPhase) */
		signalvec->data[count] = amp*cos(phase+phi0);

		/*
		double s1;
		s1 = pow(pow(values.data[3],2.0) + pow(values.data[4], 2.0) +pow(values.data[5], 2.0), 0.5);
		printf("%e %e %e %e %e\n", t, values.data[3], values.data[4], values.data[5], s1);
		printf("%e %e %e %e %e %e %e %e\n", t, signalvec->data[count], Phi, phi, psi, NCapDotL/magL, Fcross, Fplus);
		*/

		/* Record the old values of frequency, polarisation angle and phase */
		fOld = f;
		psiOld = psi;
		phiOld = phi;

		count++;
		t = (count-in->nStartPad) * dt;

		/* The velocity parameter */
		acstPars.v = v;
		/* Cast the input structure into a void struct as required by the RungeKutta4 routine */
		funcParams = (void *) &acstPars;
		/* Compute the derivatives at the current time. */
		LALACSTDerivatives(&values, &dvalues, funcParams);
		/* Supply the integration routine with derivatives and current time*/
		rk4in.dydx = &dvalues;
		rk4in.x = t;
		/* Compute the carrier phase of the signal at the current time */
		Theta = etaBy5M * (in->tC - t);
		Phi = func.phasing3(Theta, &ak);
		if (XLAL_IS_REAL8_FAIL_NAN(Phi))
			ABORTXLAL(status);
		/* Compute the post-Newtonian frequency of the signal and 'velocity' at the current time */
		f = func.frequency3(Theta, &ak);
		if (XLAL_IS_REAL8_FAIL_NAN(f))
			ABORTXLAL(status);
		v = pow(LAL_PI * in->totalMass * LAL_MTSUN_SI * f, 1.L/3.L);
		/* Integrate the equations one step forward */
		LALRungeKutta4(status->statusPtr, &newvalues, integrator, funcParams);
		CHECKSTATUSPTR(status);

		/* re-record the new values of the variables */
		for (j=0; j<nDEDim; j++) values.data[j] = newvalues.data[j];

		/* Compute the magnitude of the angular-momentum and its component along line-of-sight. */
		magL = sqrt(values.data[0]*values.data[0] + values.data[1]*values.data[1] + values.data[2]*values.data[2]);
		NCapDotL = acstPars.NCap[0]*values.data[0]+acstPars.NCap[1]*values.data[1]+acstPars.NCap[2]*values.data[2];
		NCapDotLByL = NCapDotL/magL;

		/* Compute the polarisation angle, beam factors, polarisation phase and modulated amplitude */
		LALInspiralPolarisationAngle(&psi, psiOld, &acstPars.NCap[0], &values.data[0]);
		LALInspiralBeamFactors(&Fplus, &Fcross, Fp0, Fp1, Fc0, Fc1, psi);
		LALInspiralPolarisationPhase(&phi, phiOld, NCapDotLByL, Fplus, Fcross);
		LALInspiralModulatedAmplitude(&amp, v, amp0, NCapDotLByL, Fplus, Fcross);

		/* The new phase of the signal is ... */
		phase = Phi + phi;
	}
	while (count < signalvec->length)
	{
		signalvec->data[count] = 0.L;
		count++;
	}
        XLALRungeKutta4Free( integrator );
	LALFree(dummy.data);
	DETATCHSTATUSPTR(status);
	RETURN(status);

}


static void
LALInspiralPolarisationAngle(
		REAL8 *psi,
		REAL8 psiOld,
		REAL8 *NCap,
		REAL8 *L)
{
	REAL8 NCapDotL, NCapDotZ, LCapDotZ, NCapDotLCrossZ;

	/* page 6278, Eq. (20) of ACST */
	NCapDotL = NCap[0]*L[0] + NCap[1]*L[1] + NCap[2]*L[2];
	NCapDotZ = NCap[2];
	LCapDotZ = L[2];
	NCapDotLCrossZ = NCap[0]*L[1] - NCap[1]*L[0];
	if (NCapDotLCrossZ)
		*psi = atan((LCapDotZ - NCapDotL*NCapDotZ)/NCapDotLCrossZ);
	else
		*psi = LAL_PI_2;

	/* If you require your polarisation angle to be continuous, then uncomment the line below */
	if (psiOld) while (fabs(*psi-psiOld)>LAL_PI_2) *psi = (psiOld > *psi) ? *psi+LAL_PI : *psi-LAL_PI;
}

static void
LALInspiralBeamFactors(
		REAL8 *Fplus,
	        REAL8 *Fcross,
		REAL8 Fp0,
		REAL8 Fp1,
		REAL8 Fc0,
		REAL8 Fc1,
	        REAL8 psi)
{


	REAL8 cosTwoPsi, sinTwoPsi;
	/* page 6276, Eqs. (4a) and (4b) of ACST */

	cosTwoPsi = cos(2.L*psi);
	sinTwoPsi = sin(2.L*psi);

	*Fplus = Fp0 * cosTwoPsi + Fp1 * sinTwoPsi;
	*Fcross = Fc0 * sinTwoPsi + Fc1 * cosTwoPsi;
}

static void
LALInspiralPolarisationPhase(
		REAL8 *phi,
		REAL8 phiOld,
		REAL8 NCapDotLByL,
		REAL8 Fplus,
		REAL8 Fcross
		)
{
	/* page 6278, Eqs. (19b) of ACST */
	*phi=atan(2.L*NCapDotLByL*Fcross/((1.L + NCapDotLByL*NCapDotLByL)*Fplus));

	/* If you require your polarisation phase to be continuous, then uncomment the line below */
	if (phiOld) while (fabs(*phi-phiOld)>LAL_PI_2) *phi = (phiOld>*phi) ? *phi+LAL_PI : *phi-LAL_PI;
}

static void
LALInspiralModulatedAmplitude(
		REAL8 *amp,
		REAL8 v,
		REAL8 amp0,
		REAL8 NCapDotL,
		REAL8 Fplus,
		REAL8 Fcross
		)
{
	REAL8 NCapDotLSq;
	/* page 6278, Eqs. (19a) of ACST */
	NCapDotLSq = NCapDotL * NCapDotL;
	*amp = amp0 * v*v * sqrt ( pow(1.L+NCapDotLSq,2.L) * Fplus*Fplus + 4.L*NCapDotLSq*Fcross*Fcross);
}

/* computes the derivatives of the angular mom precession eqns for rkdumb */
void
LALACSTDerivatives
(
 REAL8Vector *values,
 REAL8Vector *dvalues,
 void *funcParams
 )
{


	/* derivatives of vectors L,S1,S2 */
	/* page 6277, Eqs. (11a)-(11c) of ACST
	 * Note that ACST has a mis-print in Eq. (11b) */
	/* loop variables */
	enum { Dim=3 };
	UINT4 i, j, k, p, q;
	/* magnitudes of S1, S2, L etc. */
	REAL8 /*Theta,*/ v, magL, S1DotL, S2DotL;
        REAL8 fourM1Plus, fourM2Plus, oneBy2Mcube, Lsq, dL0, c2, v6;
        REAL8 L[Dim], S1[Dim], S2[Dim], S1CrossL[Dim], S2CrossL[Dim], S1CrossS2[Dim];
	InspiralACSTParams *ACSTIn;


	ACSTIn = (InspiralACSTParams *) funcParams;
	/* extract 'velocity' and masses */
	v = ACSTIn->v;
	fourM1Plus = ACSTIn->fourM1Plus;
	fourM2Plus = ACSTIn->fourM2Plus;
	oneBy2Mcube = ACSTIn->oneBy2Mcube;
	magL = sqrt(values->data[0]*values->data[0] + values->data[1]*values->data[1] + values->data[2]*values->data[2]);

	/* extract vectors, angular momentum and spins */
	for (i=0; i<Dim; i++)
	{
		L[i]=values->data[i];
		S1[i]=values->data[i+Dim];
		S2[i]=values->data[i+2*Dim];
	}

	S1DotL = S2DotL = 0.L;
	for (i=0; i<Dim; i++)
	{
		j = (i+1) % Dim;
		k = (i+2) % Dim;
		S1DotL += S1[i] * L[i];
		S2DotL += S2[i] * L[i];
		S1CrossL[i] = S1[j] * L[k] - S1[k] * L[j];
		S2CrossL[i] = S2[j] * L[k] - S2[k] * L[j];
		S1CrossS2[i] = S1[j] * S2[k] - S1[k] * S2[j];
	}


	Lsq = magL*magL;
	dL0 = ACSTIn->thirtytwoBy5etc/magL;
	c2 = ACSTIn->threeBy2Mcube/Lsq;
	v6 = pow(v,6.L);

	/*
	printf("%e %e %e\n", L[0], L[1], L[2]);
	printf("%e %e %e\n", S1[0], S1[1], S1[2]);
	printf("%e %e %e\n", S2[0], S2[1], S2[2]);
	printf("%e %e %e\n", S1CrossS2[0], S1CrossS2[1], S1CrossS2[2]);
	printf("%e %e %e\n", S1CrossL[0], S1CrossL[1], S1CrossL[2]);
	printf("%e %e %e %e %e %e\n", S1DotL, S2DotL, Lsq, dL0, c2, v6);
	*/

	for(i=0;i<Dim;i++)
	{
		p = i+Dim;
		q = i+2*Dim;
		/* compute the derivatives */
		dvalues->data[i] = ((fourM1Plus - c2 * S2DotL) * S1CrossL[i]
	              +  (fourM2Plus - c2 * S1DotL) * S2CrossL[i] - dL0 * L[i]*v) * v6;;

		dvalues->data[p] = ((c2 * S2DotL - fourM1Plus) * S1CrossL[i] - oneBy2Mcube * S1CrossS2[i]) * v6;

		dvalues->data[q] = ((c2 * S1DotL - fourM2Plus) * S2CrossL[i] + oneBy2Mcube * S1CrossS2[i]) * v6;
	}

}



/* Routine to inject inspiral waveforms from binaries consisting of spinning objects */


/* NOT DONE FOR THE MOMENT should remove the polarisation effects which are already
   taken into account in inject package */
void
LALInspiralSpinModulatedWaveForInjection(
					 LALStatus        *status,
					 CoherentGW *waveform,
					 InspiralTemplate *params,
					 PPNParamStruc  *ppnParams
					 )
{


  UINT4 nDEDim=9/*, Dim=3*/;
  UINT4 i,j, count, length = 0;
  /* polarisation and antenna pattern */
  REAL8 omega,psi, psiOld, Fplus, Fcross, Fp0, Fp1, Fc0, Fc1, magL, amp, amp0, phi, phiOld, /*Amplitude,*/ Phi/*, dPhi*/;
  REAL8 v, t, tMax, dt, f, fOld, fn, phase, MCube, Theta, etaBy5M, NCapDotL, NCapDotLByL;
  /* source direction angles */
  InspiralACSTParams acstPars;
  REAL8Vector dummy, values, dvalues, newvalues, yt, dym, dyt;
  void *funcParams;
    rk4In rk4in;
    rk4GSLIntegrator *integrator;
    expnFunc func;
    expnCoeffs ak;
    REAL4Vector *a=NULL;
    REAL4Vector *ff=NULL ;
    REAL8Vector *phiv=NULL;
    REAL8 unitHz,f2a, mu, mTot, cosI, etab, fFac,  f2aFac, apFac, acFac, phiC;
    CreateVectorSequenceIn in;

    mTot   =  params->mass1 + params->mass2;
    etab   =  params->mass1 * params->mass2;
    etab  /= mTot;
    etab  /= mTot;
    unitHz = (mTot) *LAL_MTSUN_SI*(REAL8)LAL_PI;
    cosI   = cos( params->inclination );
    mu     = etab * mTot;
    fFac   = 1.0 / ( 4.0*LAL_TWOPI*LAL_MTSUN_SI*mTot );
    dt     = -1. * etab / ( params->tSampling * 5.0*LAL_MTSUN_SI*mTot );
    f2aFac = LAL_PI*LAL_MTSUN_SI*mTot*fFac;
    apFac  = acFac = -2.0 * mu * LAL_MRSUN_SI/params->distance;
    apFac *= 1.0 + cosI*cosI;
    acFac *= 2.0*cosI;


    INITSTATUS(status);
    ATTATCHSTATUSPTR(status);
    ASSERT(params,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
    /* Make sure waveform fields don't exist. */
    ASSERT( !( waveform->a ), status, LALINSPIRALH_ENULL,
	    LALINSPIRALH_MSGENULL );
    ASSERT( !( waveform->f ), status, LALINSPIRALH_ENULL,
	    LALINSPIRALH_MSGENULL );
    ASSERT( !( waveform->phi ), status, LALINSPIRALH_ENULL,
	    LALINSPIRALH_MSGENULL );
    ASSERT( !( waveform->shift ), status, LALINSPIRALH_ENULL,
	    LALINSPIRALH_MSGENULL );


    LALInspiralSetup (status->statusPtr, &ak, params);
    CHECKSTATUSPTR(status);
    LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
    CHECKSTATUSPTR(status);
    LALInspiralWaveLength(status->statusPtr, &length, *params);
    CHECKSTATUSPTR(status);
    LALInspiralParameterCalc(status->statusPtr, params);
    CHECKSTATUSPTR(status);

    /* memset(&ff,0,sizeof(REAL4TimeSeries));*/
    LALSCreateVector(status->statusPtr, &ff, length);
    CHECKSTATUSPTR(status);

    LALSCreateVector(status->statusPtr, &a, 2*length);
    CHECKSTATUSPTR(status);

    LALDCreateVector(status->statusPtr, &phiv, length);
    CHECKSTATUSPTR(status);


    /* Allocate space for all the vectors in one go ... */

    dummy.length = nDEDim * 6;
    values.length = dvalues.length = newvalues.length = yt.length = dym.length = dyt.length = nDEDim;

    if (!(dummy.data = (REAL8 *) LALMalloc(sizeof(REAL8) * dummy.length)))
      {
	ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
      }

    /* and distribute as necessary */
    values.data = &dummy.data[0];             /* values of L, S1, S2 at current time */
    dvalues.data = &dummy.data[nDEDim];       /* values of dL/dt, dS1/dt, dS2/dt at current time */
    newvalues.data = &dummy.data[2*nDEDim];   /* new values of L, S1, S2 after time increment by dt*/
    yt.data = &dummy.data[3*nDEDim];          /* Vector space needed by LALRungeKutta4 */
    dym.data = &dummy.data[4*nDEDim];         /* Vector space needed by LALRungeKutta4 */
    dyt.data = &dummy.data[5*nDEDim];         /* Vector space needed by LALRungeKutta4 */

    v = pow(LAL_PI * params->totalMass * LAL_MTSUN_SI * params->fLower, 1.L/3.L);
    MCube = pow(LAL_MTSUN_SI * params->totalMass, 3.L);
    /* Fill in the structure needed by the routine that computes the derivatives */
    /* constant spins of the two bodies */
    acstPars.magS1 = params->spin1[0]*params->spin1[0] + params->spin1[1]*params->spin1[1] + params->spin1[2]*params->spin1[2];
    acstPars.magS2 = params->spin2[0]*params->spin2[0] + params->spin2[1]*params->spin2[1] + params->spin2[2]*params->spin2[2];
    /* Direction to the source in the solar system barycenter */
    acstPars.NCap[0] = sin(params->sourceTheta)*cos(params->sourcePhi);
    acstPars.NCap[1] = sin(params->sourceTheta)*sin(params->sourcePhi);
    acstPars.NCap[2] = cos(params->sourceTheta);
    /* total mass in seconds */
    acstPars.M = LAL_MTSUN_SI * params->totalMass;
    /* Combination of masses that appears in the evolution of L, S1 and S2 */
    acstPars.fourM1Plus = (4.L*params->mass1 + 3.L*params->mass2)/(2.*params->mass1*MCube);
    acstPars.fourM2Plus = (4.L*params->mass2 + 3.L*params->mass1)/(2.*params->mass2*MCube);
    acstPars.oneBy2Mcube = 0.5L/MCube;
    acstPars.threeBy2Mcube = 1.5L/MCube;
    acstPars.thirtytwoBy5etc = (32.L/5.L) * pow(params->eta,2.) * acstPars.M;

    /* Constant amplitude of GW */
    amp0 = 2.L*params->eta * acstPars.M/params->distance;
    /* Initial magnitude of the angular momentum, determined by the initial frequency/velocity */
    magL = params->eta * (acstPars.M * acstPars.M)/v;
    /* Initial angular momentum and spin vectors */
    /* Angular momentum */
    values.data[0] = magL * sin(params->orbitTheta0) * cos(params->orbitPhi0);
    values.data[1] = magL * sin(params->orbitTheta0) * sin(params->orbitPhi0);
    values.data[2] = magL * cos(params->orbitTheta0);
    /* Spin of primary */
    values.data[3] = params->spin1[0];
    values.data[4] = params->spin1[1];
    values.data[5] = params->spin1[2];
    /* Spin of secondarfy */
    values.data[6] = params->spin2[0];
    values.data[7] = params->spin2[1];
    values.data[8] = params->spin2[2];
    /* Structure needed by RungeKutta4 for the evolution of the differential equations */
    rk4in.function = LALACSTDerivatives;
    rk4in.y = &values;
    rk4in.h = 1.L/params->tSampling;
    rk4in.n = nDEDim;
    rk4in.yt = &yt;
    rk4in.dym = &dym;
    rk4in.dyt = &dyt;

    xlalErrno = 0;
    /* Initialize GSL integrator */
    if (!(integrator = XLALRungeKutta4Init(nDEDim, &rk4in)))
    {
      INT4 errNum = XLALClearErrno();
      LALFree(dummy.data);

      if (errNum == XLAL_ENOMEM)
        ABORT(status, LALINSPIRALH_EMEM, LALINSPIRALH_MSGEMEM);
      else
        ABORTXLAL( status );
    }

    /* Pad the first nStartPad elements of the signal array with zero */
    count = 0;


    /* Get the initial conditions before the evolution */
    t = 0.L;                                                 /* initial time */
    dt = 1.L/params->tSampling;                                  /* Sampling interval */
    etaBy5M = params->eta/(5.L*acstPars.M);                      /* Constant that appears in ... */
    Theta = etaBy5M * (params->tC - t);                          /* ... Theta */
    tMax = params->tC - dt;                                      /* Maximum duration of the signal */
    fn = (ak.flso > params->fCutoff) ? ak.flso : params->fCutoff;    /* Frequency cutoff, smaller of user given f or flso */

    Phi = func.phasing3(Theta, &ak);                             /* Carrier phase at the initial time */
    if (XLAL_IS_REAL8_FAIL_NAN(Phi))
        ABORTXLAL(status);
    f = func.frequency3(Theta, &ak);                             /* Carrier Freqeuncy at the initial time */
    if (XLAL_IS_REAL8_FAIL_NAN(f))
        ABORTXLAL(status);

    /* Constants that appear in the antenna pattern */
    Fp0 = (1.L + cos(params->sourceTheta)*cos(params->sourceTheta)) * cos(2.L*params->sourcePhi);
    Fp1 = cos(params->sourceTheta) * sin(2.L*params->sourcePhi);
    Fc0 = Fp0;
    Fc1 = -Fp1;

    psiOld = 0.L;
    phiOld = 0.L;
    magL = sqrt(values.data[0]*values.data[0] + values.data[1]*values.data[1] + values.data[2]*values.data[2]);
    NCapDotL = acstPars.NCap[0]*values.data[0] + acstPars.NCap[1]*values.data[1] + acstPars.NCap[2]*values.data[2];
    NCapDotLByL = NCapDotL/magL;

    /* Initial values of */
    /* the polarization angle */
    LALInspiralPolarisationAngle(&psi, psiOld, &acstPars.NCap[0], &values.data[0]);
    /* beam pattern factors */
    LALInspiralBeamFactors(&Fplus, &Fcross, Fp0, Fp1, Fc0, Fc1, psi);
    /* the polarization phase */
    LALInspiralPolarisationPhase(&phi, phiOld, NCapDotLByL, Fplus, Fcross);
    /* modulated amplitude */
    LALInspiralModulatedAmplitude(&amp, v, amp0, NCapDotLByL, Fplus, Fcross);


    phase = Phi + phi;

    fOld = f - 0.1;
    while (f < fn && t < tMax && f>fOld)
      {
	/*ASSERT((INT4)count < (INT4)signalvec->length, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);*/
	/* Subtract the constant initial phase (chosen to have a phase of in->startPhase) */



	omega = v*v*v;
	ff->data[count]= (REAL4)(omega/unitHz);
	f2a = pow (f2aFac * omega, 2./3.);
	a->data[2*count]          = (REAL4)(4.*apFac * f2a);
	a->data[2*count+1]        = (REAL4)(4.*acFac * f2a);
	phiv->data[count]          = (REAL8)(phase);



	/*
	  double s1;
	  s1 = pow(pow(values.data[3],2.0) + pow(values.data[4], 2.0) +pow(values.data[5], 2.0), 0.5);
	  printf("%e %e %e %e %e\n", t, values.data[3], values.data[4], values.data[5], s1);
	  printf("%e %e %e %e %e %e %e %e\n", t, signalvec->data[count], Phi, phi, psi, NCapDotL/magL, Fcross, Fplus);
	*/

	/* Record the old values of frequency, polarisation angle and phase */
	fOld = f;
	psiOld = psi;
	phiOld = phi;

	count++;
	t = (count-params->nStartPad) * dt;

	/* The velocity parameter */
	acstPars.v = v;
	/* Cast the input structure into a void struct as required by the RungeKutta4 routine */
	funcParams = (void *) &acstPars;
	/* Compute the derivatives at the current time. */
	LALACSTDerivatives(&values, &dvalues, funcParams);
	/* Supply the integration routine with derivatives and current time*/
	rk4in.dydx = &dvalues;
	rk4in.x = t;
	/* Compute the carrier phase of the signal at the current time */
	Theta = etaBy5M * (params->tC - t);
	Phi = func.phasing3(Theta, &ak);
	if (XLAL_IS_REAL8_FAIL_NAN(Phi))
		ABORTXLAL(status);
	/* Compute the post-Newtonian frequency of the signal and 'velocity' at the current time */
	f = func.frequency3(Theta, &ak);
	if (XLAL_IS_REAL8_FAIL_NAN(Phi))
		ABORTXLAL(status);
	v = pow(LAL_PI * params->totalMass * LAL_MTSUN_SI * f, 1.L/3.L);
	/* Integrate the equations one step forward */
	LALRungeKutta4(status->statusPtr, &newvalues, integrator, funcParams);
	CHECKSTATUSPTR(status);

	/* re-record the new values of the variables */
	for (j=0; j<nDEDim; j++) values.data[j] = newvalues.data[j];

	/* Compute the magnitude of the angular-momentum and its component along line-of-sight. */
	magL = sqrt(values.data[0]*values.data[0] + values.data[1]*values.data[1] + values.data[2]*values.data[2]);
	NCapDotL = acstPars.NCap[0]*values.data[0]+acstPars.NCap[1]*values.data[1]+acstPars.NCap[2]*values.data[2];
	NCapDotLByL = NCapDotL/magL;

	/* Compute the polarisation angle, beam factors, polarisation phase and modulated amplitude */
	LALInspiralPolarisationAngle(&psi, psiOld, &acstPars.NCap[0], &values.data[0]);
	LALInspiralBeamFactors(&Fplus, &Fcross, Fp0, Fp1, Fc0, Fc1, psi);
	LALInspiralPolarisationPhase(&phi, phiOld, NCapDotLByL, Fplus, Fcross);
	LALInspiralModulatedAmplitude(&amp, v, amp0, NCapDotLByL, Fplus, Fcross);

	/* The new phase of the signal is ... */
	phase = Phi + phi;
      }


    params->fFinal = f;

    /*wrap the phase vector*/
    phiC =  phiv->data[count-1] ;
    for (i=0; i<count;i++)
      {
	phiv->data[i] =  phiv->data[i] -phiC;
      }
    /* Allocate the waveform structures. */
    if ( ( waveform->a = (REAL4TimeVectorSeries *)
	   LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
      ABORT( status, LALINSPIRALH_EMEM,
	     LALINSPIRALH_MSGEMEM );
    }
    memset( waveform->a, 0, sizeof(REAL4TimeVectorSeries) );
    if ( ( waveform->f = (REAL4TimeSeries *)
	   LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
      LALFree( waveform->a ); waveform->a = NULL;
      ABORT( status, LALINSPIRALH_EMEM,
	     LALINSPIRALH_MSGEMEM );
    }
    memset( waveform->f, 0, sizeof(REAL4TimeSeries) );
    if ( ( waveform->phi = (REAL8TimeSeries *)
	   LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
      LALFree( waveform->a ); waveform->a = NULL;
      LALFree( waveform->f ); waveform->f = NULL;
      ABORT( status, LALINSPIRALH_EMEM,
	     LALINSPIRALH_MSGEMEM );
    }
    memset( waveform->phi, 0, sizeof(REAL8TimeSeries) );



    in.length = (UINT4)count;
    in.vectorLength = 2;
    LALSCreateVectorSequence( status->statusPtr,
			      &( waveform->a->data ), &in );
    CHECKSTATUSPTR(status);
    LALSCreateVector( status->statusPtr,
		      &( waveform->f->data ), count);
    CHECKSTATUSPTR(status);
    LALDCreateVector( status->statusPtr,
		      &( waveform->phi->data ), count );
    CHECKSTATUSPTR(status);




    memcpy(waveform->f->data->data , ff->data, count*(sizeof(REAL4)));
    memcpy(waveform->a->data->data , a->data, 2*count*(sizeof(REAL4)));
    memcpy(waveform->phi->data->data ,phiv->data, count*(sizeof(REAL8)));



    dt = -1. * etab / ( params->tSampling * 5.0*LAL_MTSUN_SI*mTot );

    waveform->a->deltaT = waveform->f->deltaT = waveform->phi->deltaT
      = 1./params->tSampling;

    waveform->a->sampleUnits = lalStrainUnit;
    waveform->f->sampleUnits = lalHertzUnit;
    waveform->phi->sampleUnits = lalDimensionlessUnit;
    waveform->position = ppnParams->position;
    waveform->psi = ppnParams->psi;


    snprintf( waveform->a->name, LALNameLength, "Spinning inspiral amplitudes" );
    snprintf( waveform->f->name, LALNameLength, "Spinning inspiral frequency" );
    snprintf( waveform->phi->name, LALNameLength, "Spinning inspiral phase" );


    params->tC = count / params->tSampling ;
    params->tSampling = (REAL4)(waveform->f->data->data[count-1]
				-
				waveform->f->data->data[count-2]);
    /* - (waveform->f->data[count-2]);*/
    params->tSampling /= (REAL4)dt;
    params->nStartPad = count;

    LALFree(a->data);
    LALFree(ff->data);
    LALFree(phiv->data);

    XLALRungeKutta4Free( integrator );
    LALFree(dummy.data);
    DETATCHSTATUSPTR(status);
    RETURN(status);

  }

