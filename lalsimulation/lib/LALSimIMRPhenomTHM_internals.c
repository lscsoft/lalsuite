/*
 * Copyright (C) 2020 Hector Estelles
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
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */


/**
 * \author Hector Estelles
 */

#include <math.h>

#include "LALSimIMRPhenomXUtilities.h"
#include "LALSimIMRPhenomTHM_fits.c"

/* LAL Header Files */
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/Units.h>

/* GSL Header Files */
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_roots.h>

/* ************************************* */
/* ******* STRUCTS INITIALIZATION ****** */
/* ************************************* */

/* Function for populating the IMRPhenomTWaveformStruct.
   This will compute and store the needed parameters for evaluating the fits and phase/amplitude related quantities. */

int IMRPhenomTSetWaveformVariables(
	IMRPhenomTWaveformStruct *wf,  /* Waveform struct */
	const REAL8 m1_SI,			   /* Mass of companion 1 (kg) */
  	const REAL8 m2_SI,             /* Mass of companion 2 (kg) */
  	const REAL8 chi1L_In,          /* Dimensionless aligned spin of companion 1 */
  	const REAL8 chi2L_In,          /* Dimensionless aligned spin of companion 2 */
  	const REAL8 distance,          /* Luminosity distance (m) */
  	const REAL8 deltaT,            /* sampling interval (s) */
  	const REAL8 fmin,              /* starting GW frequency (Hz) */
  	const REAL8 fRef,              /* reference GW frequency (Hz) */
  	const REAL8 phiRef,             /* reference orbital phase (rad) */
  	LALDict *lalParams				/* Dictionary with LAL parameters */
)
{
	/* Rescale to mass in solar masses */
	REAL8 m1_In      = m1_SI / LAL_MSUN_SI; // Mass 1 in solar masses
	REAL8 m2_In      = m2_SI / LAL_MSUN_SI; // Mass 2 in solar masses

	REAL8 m1, m2, chi1L, chi2L;

	/* Check if m1 >= m2, if not then swap masses/spins */
	if(m1_In >= m2_In)
	{
		chi1L = chi1L_In;
		chi2L = chi2L_In;
		m1    = m1_In;
		m2    = m2_In;
	}
	else
	{
		XLAL_PRINT_WARNING("Warning: m1 < m2, swapping the masses and spins.\n");
		chi1L = chi2L_In;
		chi2L = chi1L_In;
		m1    = m2_In;
		m2    = m1_In;
	}

	// Symmetric mass ratio
	REAL8 delta = fabs((m1 - m2) / (m1+m2));
	REAL8 eta   = fabs(0.25 * (1.0 - delta*delta) ); // Symmetric mass ratio, use fabs to prevent negative sign due to roundoff
	REAL8 q   = ((m1 > m2) ? (m1 / m2) : (m2 / m1));

	/* If eta > 0.25, e.g. roundoff, then set to 0.25 */
	if(eta > 0.25)
	{
		eta = 0.25;
	}
	if(eta == 0.25) q = 1.;

	/* Masses definitions. Note that m1 and m2 are the component masses in solar masses. */
	wf->m1_SI     = m1 * LAL_MSUN_SI;						// Mass 1 (larger) in SI units
	wf->m2_SI     = m2 * LAL_MSUN_SI;						// Mass 2 (smaller) in SI units
	wf->q         = q;										// mass ratio q>1
	wf->eta       = eta;									// Symmetric mass ratio
	wf->Mtot_SI   = wf->m1_SI + wf->m2_SI;					// Total mass in SI units
	wf->Mtot      = m1 + m2;								// Total mass in solar masses
	wf->m1        = m1 / (wf->Mtot);						// Mass 1 (larger), dimensionless (i.e. m1 \in [0,1])
	wf->m2        = m2 / (wf->Mtot);						// Mass 2 (smaller), dimensionless
	wf->M_sec     = wf->Mtot * LAL_MTSUN_SI;				// Conversion factor Hz to geometric frequency, total mass in seconds
	wf->delta     = delta;									// PN symmetry coefficient
	wf->dist_sec  = distance * (LAL_MTSUN_SI/LAL_MRSUN_SI);	// distance in seconds
	wf->phiRef  = phiRef;									// Reference orbital phase

	/* Spins */
	wf->chi1L     = chi1L;
	wf->chi2L     = chi2L;

	/* Spin parameterisations for calling the calibrated fits*/
	wf->Shat     = ((wf->m1)*(wf->m1)*chi1L + (wf->m2)*(wf->m2)*chi2L)/((wf->m1)*(wf->m1) + (wf->m2)*(wf->m2));
	wf->dchi      = chi1L - chi2L;

	/* Final Mass and Spin (we employ the XLAL functions of PhenomX) */
	wf->Mfinal    = XLALSimIMRPhenomXFinalMass2017(wf->eta,wf->chi1L,wf->chi2L);
	wf->afinal    = XLALSimIMRPhenomXFinalSpin2017(wf->eta,wf->chi1L,wf->chi2L);

	wf->afinal_prec = XLALSimIMRPhenomXFinalSpin2017(wf->eta,wf->chi1L,wf->chi2L);

	/* Distance*/
	wf->distance    = distance;

	/* Impose fRef=fmin if fRef=0 */
	wf->fRef = fRef;
	if(fRef==0.0){wf->fRef = fmin;}
	wf->fmin = fmin;

	/*Dimensionless minimum and reference frequency */
	wf->MfRef = wf->M_sec*wf->fRef;
	wf->Mfmin = wf->M_sec*wf->fmin;

	wf->deltaT = deltaT;
	wf->dtM = deltaT/wf->M_sec; // Dimensionless sampling rate

	wf->ampfac = (wf->M_sec/wf->dist_sec); // Amplitude conversion factor from NR to physical units

	wf->inspVersion = XLALSimInspiralWaveformParamsLookupPhenomTHMInspiralVersion(lalParams); // This will select if (2,2) inspiral phase/frequency are reconstructed with 3 or 4 regions

	return XLAL_SUCCESS;
}

/* Function for populating the IMRPhenomTPhase22Struct.
   This will solve the linear systems for obtaining ansatz coefficients from 22 frequency collocation points.
   PN Coefficients of TaylorT3 approximant employed in the inspiral ansatz are defined here.
   Minimum and reference time, as well as waveform length, will be computed here. */

int IMRPhenomTSetPhase22Coefficients(IMRPhenomTPhase22Struct *pPhase, IMRPhenomTWaveformStruct *wf){

	/* ***** Initialize parameters and variables ****** */

	REAL8 eta = wf->eta;     // Symmetric mass ratio
	REAL8 chi1L = wf->chi1L; // Dimensionless aligned spin of companion 1
	REAL8 chi2L = wf->chi2L; // Dimensionless aligned spin of companion 2
	REAL8 S = wf->Shat;      // Effective dimensionless spin parameter
	REAL8 dchi = wf->dchi;   // Dimensionless spin difference chi1L - chi2L
	REAL8 delta = wf->delta; // Asymmetry parameter
	pPhase->dtM = wf->dtM; // Dimensionless sampling rate
	
	/* Calibrated value of TaylorT3 t0 parameter for matching frequency at TaylorT3 theta=0.33.
	If 4 regions are selected for reconstruction, first inspiral region will be employ TaylorT3 with this value of t0.
	With this, it is ensured that the early inspiral is described by the proper TaylorT3 PN description.
	If 3 regions are selected, this will be still employed for computing the value of TaylorT3 at theta=0.33 for setting the first collocation point.*/
	REAL8 tt0 = IMRPhenomT_Inspiral_TaylorT3_t0(eta, S, dchi, delta);
	pPhase->tt0 = tt0;

	/* Ringdown and damping frequency of final BH for the 22 mode (defined on LALSimIMRPhenomTHM_fits.c) */
	REAL8 fRING     = evaluate_QNMfit_fring22(wf->afinal) / (wf->Mfinal); // 22 mode ringdown frequency
	REAL8 fDAMP     = evaluate_QNMfit_fdamp22(wf->afinal) / (wf->Mfinal); // 22 mode damping frequency of ground state (n=1) QNM

	/* Angular ringdown and damping frequencies (omega =2*pi*f) */
	REAL8 alpha1RD  = 2*LAL_PI*fDAMP;
	REAL8 omegaRING = 2*LAL_PI*fRING;

	pPhase->omegaRING = omegaRING;
	pPhase->alpha1RD  = alpha1RD;

	pPhase->omegaRING_prec = 2*LAL_PI*evaluate_QNMfit_fring22(wf->afinal_prec) / (wf->Mfinal);
	pPhase->alpha1RD_prec  = 2*LAL_PI*evaluate_QNMfit_fdamp22(wf->afinal_prec) / (wf->Mfinal);

	pPhase->EulerRDslope = 2*LAL_PI*(evaluate_QNMfit_fring22(wf->afinal_prec) / (wf->Mfinal) - evaluate_QNMfit_fring21(wf->afinal_prec) / (wf->Mfinal)); // FIXME: 
	if(wf->afinal<0)
	{
		pPhase->EulerRDslope = -pPhase->EulerRDslope;
	}

	/* Dimensionless minimum and reference frequencies */
	pPhase->MfRef = wf->MfRef;
	pPhase->Mfmin = wf->Mfmin;

	/* GSL objects for solving system of equations via LU decomposition */
	gsl_vector *b, *x;
	gsl_matrix *A;
	gsl_permutation *p;
	int s; /* Sign of permutation */

	/* ********************************** */
	/* *** RINGDOWN COEFFICIENTS ******** */
	/* ********************************** */

	/* For phase and frequency ringdown ansatz, coefficients are obtained from ringdown and damping frequencies, frequency at peak amplitude
	   and phenomenological fits of two free coefficients. No linear system solving is needed.
	   Coefficient c1 is defined in Eq. [9] of Damour&Nagar 2014 (10.1103/PhysRevD.90.024054, https://arxiv.org/abs/1406.0401).
	   Coefficients c2 and c3 are calibrated to NR/Teukolsky data.
	   Coefficient c4 is set to zero.
	   Explained also in Sec II.C, eq. 26 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */

	pPhase->omegaPeak = IMRPhenomT_PeakFrequency_22(eta, S, dchi, delta); // 22 Frequency at the peak 22 amplitude time (t=0 by definition)

	pPhase->c3 = IMRPhenomT_RD_Freq_D3_22(eta, S, dchi, delta);
	pPhase->c2 = IMRPhenomT_RD_Freq_D2_22(eta, S, dchi, delta);
	pPhase->c4 = 0.0;
	pPhase->c1 = (1 + pPhase->c3 + pPhase->c4)*(pPhase->omegaRING - pPhase->omegaPeak)/pPhase->c2/(pPhase->c3 + 2*pPhase->c4);
	pPhase->c1_prec = (1 + pPhase->c3 + pPhase->c4)*(pPhase->omegaRING_prec - pPhase->omegaPeak)/pPhase->c2/(pPhase->c3 + 2*pPhase->c4);

	/* ********************************** */
	/* *** INSPIRAL COEFFICIENTS ******** */
	/* ********************************** */

	/* PN coefficients of the TaylorT3 approximant employed in the inspiral ansatz are defined here.
	   Unknown coefficients for the higher order extension in the inspiral ansatz are obtained through fitted collocation points
	   solving a linear system of equations. */


	/*PN coefficients corresponding to the TaylorT3 implementation of the model,
	  as defined in Appendix A1 of Estelles et al 2020 (https://arxiv.org/pdf/2004.08302.pdf) */

	pPhase->omega1PN = 0.27641369047619047 + (11*eta)/32.;

	pPhase->omega1halfPN = (-19*(chi1L + chi2L)*eta)/80. + (-113*(chi2L*(-1 + delta) - chi1L*(1 + delta)) - 96*LAL_PI)/320.;

	pPhase->omega2PN = (1855099 + 1714608*chi2L*chi2L*(-1 + delta) - 1714608*chi1L*chi1L*(1 + delta))/1.4450688e7 + ((56975 + 61236*chi1L*chi1L - 119448*chi1L*chi2L + 61236*chi2L*chi2L)*eta)/258048. + \
   (371*eta*eta)/2048.;

    pPhase->omega2halfPN = (-17*(chi1L + chi2L)*eta*eta)/128. + (-146597*(chi2L*(-1 + delta) - chi1L*(1 + delta)) - 46374*LAL_PI)/129024. + (eta*(-2*(chi1L*(1213 - 63*delta) + chi2L*(1213 + 63*delta)) + \
     117*LAL_PI))/2304.;

    pPhase->omega3PN = -2.499258364444952 - (16928263*chi1L*chi1L)/1.376256e8 - (16928263*chi2L*chi2L)/1.376256e8 - (16928263*chi1L*chi1L*delta)/1.376256e8 + (16928263*chi2L*chi2L*delta)/1.376256e8 + \
   ((-2318475 + 18767224*chi1L*chi1L - 54663952*chi1L*chi2L + 18767224*chi2L*chi2L)*eta*eta)/1.376256e8 + (235925*eta*eta*eta)/1.769472e6 + (107*LAL_GAMMA)/280. - (6127*chi1L*LAL_PI)/12800. - \
   (6127*chi2L*LAL_PI)/12800. - (6127*chi1L*delta*LAL_PI)/12800. + (6127*chi2L*delta*LAL_PI)/12800. + (53*LAL_PI*LAL_PI)/200. + \
   (eta*(632550449425 + 35200873512*chi1L*chi1L - 28527282000*chi1L*chi2L + 9605339856*chi1L*chi1L*delta - 1512*chi2L*chi2L*(-23281001 + 6352738*delta) + 34172264448*(chi1L + chi2L)*LAL_PI - \
        22912243200*LAL_PI*LAL_PI))/1.040449536e11 + (107*log(2))/280.;

    pPhase->omega3halfPN = (-12029*(chi1L + chi2L)*eta*eta*eta)/92160. + (eta*eta*(507654*chi1L*chi2L*chi2L - 838782*chi2L*chi2L*chi2L + chi2L*(-840149 + 507654*chi1L*chi1L - 870576*delta) +
        chi1L*(-840149 - 838782*chi1L*chi1L + 870576*delta) + 1701228*LAL_PI))/1.548288e7 +
   (eta*(218532006*chi1L*chi2L*chi2L*(-1 + delta) - 1134*chi2L*chi2L*chi2L*(-206917 + 71931*delta) - chi2L*(1496368361 - 429508815*delta + 218532006*chi1L*chi1L*(1 + delta)) +
        chi1L*(-1496368361 - 429508815*delta + 1134*chi1L*chi1L*(206917 + 71931*delta)) - 144*(488825 + 923076*chi1L*chi1L - 1782648*chi1L*chi2L + 923076*chi2L*chi2L)*LAL_PI))/1.8579456e8 +
   (-6579635551*chi2L*(-1 + delta) + 535759434*chi2L*chi2L*chi2L*(-1 + delta) - chi1L*(-6579635551 + 535759434*chi1L*chi1L)*(1 + delta) +
      (-565550067 - 465230304*chi2L*chi2L*(-1 + delta) + 465230304*chi1L*chi1L*(1 + delta))*LAL_PI)/1.30056192e9;


   /* Solve inspiral coefficients */

   /* Inspiral reconstruction is based on the TaylorT3 approximant at 3.5PN with spin-orbit and spin-spin contributions, whose coefficients are defined above.
   TaylorT3 is parameterised in terms of a 'PN time parameter' theta=pow(eta*(tt0 - tEarly)/5,-1./8), which depends on an integration constant tt0 (coming from the TaylorT2 derivation) which has to be fixed
   in order to impose a determined frequency at a determined time. TaylorT3 will diverge when t=tt0, so this parameter can be understood as the merger time prediction of TaylorT3, which in general will not be
   the actual merger time. In order to overcome the divergence issue in cases where this occurs too early, we originally decided to set tt0 to the peak time of the 22. This produces a small frequency shift that
   over many cycles can reduce the accuracy of the model for long waveforms. To overcome this, an early collocation point at theta=0.33 is set to the actual value of TaylorT3 (through a tt0 fit) to impose
   the low frequency regime to follow the correct description. Along with this, the early inspiral region with theta<0.33 was selected to be modeled with pure TaylorT3 with the tt0 fit value.
   However, several tests suggest that actually imposing the theta=0.33 value is enough and two regions are not needed. We let as an option in the model to select the inspiral reconstruction with or without this split.
   Depending on the value of the LAL parameter 'PhenomTHMInspiralVersion' reconstruction of the inspiral phase can be splitted into two different regions. For default value (0) this will not be done.
   In default reconstruction, TaylorT3 with t0=0 will be employed, plus 6 additional higher order terms to be solved with the value of 5 fitted collocation points at fixed theta positions {0.45, 0.55, 0.65, 0.75, 0.82}
   + and early collocation point at theta=0.33 whose value will be the value of TaylorT3 with t0 computed from a fit.
   In non-default reconstruction, region with theta<0.33 will be modelled by TaylorT3 with the tt0 fit value and without extra coefficients.*/

   /* Initialize higher order extension coefficients to 0. In this way, when calling the ansatz IMRPhenomTInspiralOmegaAnsatz22 it will return PN TaylorT3 value. */
	pPhase->omegaInspC1 = 0.0;
	pPhase->omegaInspC2 = 0.0;
	pPhase->omegaInspC3 = 0.0;
	pPhase->omegaInspC4 = 0.0;
	pPhase->omegaInspC5 = 0.0;
	pPhase->omegaInspC6 = 0.0;

   /* Set collocation points */

    REAL8 thetapoints[6] = {0.33,0.45,0.55,0.65,0.75,0.82};  // Collocation point times as defined in Eq. 11 of PhenomTHM paper https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */

    REAL8 omegainsppoints[6];

    /* Boundary time between early and late inspiral regions (only needed if non-default 4 region reconstruction is selected). Boundary is fixed at theta=0.33 for all cases. Inverting the definition of theta, the boundary time is obtained.
    Each inspiral region (if non-default reconstruction is selected) has a different parameterization. For the early inspiral region, theta is defined employing the calibrated quantity t0 (tt0 in the code),
    while in the late inspiral region t0 is set to 0 for avoiding singularities. */
    REAL8 tEarly = -5.0/(eta*pow(thetapoints[0],8));
    pPhase->tEarly = tEarly;
    REAL8 thetaini = pow(eta*(tt0 - tEarly)/5,-1./8); // This thetaini is the reparameterization of theta=0.33 for the early inspiral region theta definition.

    /* initialize collocation point values for solving inspiral coefficient system */
    omegainsppoints[0] = IMRPhenomTTaylorT3(thetaini, pPhase);
    omegainsppoints[1] = IMRPhenomT_Inspiral_Freq_CP1_22(eta, S, dchi, delta);
    omegainsppoints[2] = IMRPhenomT_Inspiral_Freq_CP2_22(eta, S, dchi, delta);
    omegainsppoints[3] = IMRPhenomT_Inspiral_Freq_CP3_22(eta, S, dchi, delta);
    omegainsppoints[4] = IMRPhenomT_Inspiral_Freq_CP4_22(eta, S, dchi, delta);
    omegainsppoints[5] = IMRPhenomT_Inspiral_Freq_CP5_22(eta, S, dchi, delta);

	/*Set linear system, which is rank 6 */
	p = gsl_permutation_alloc(6);
	b = gsl_vector_alloc(6);
	x = gsl_vector_alloc(6);
	A = gsl_matrix_alloc(6,6);

	/*Set A matrix and b vector*/

	REAL8 theta, theta8, theta9, theta10, theta11, theta12, theta13, T3offset; // Initialize theta powers the diagonal A matrix
	/* theta is a weighted dimensionless time parameter defined in Eq. 315 of Blanchet 2014 (https://arxiv.org/abs/1310.1528).
	   Notice however that here is defined at the (-1/8) power, in order to represent the Post-Newtonian order (1/c). */

	/* Set up inspiral coefficient system.
	This system of equations is explained in eq. 10 of THM paper https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */
	for (UINT4 idx=0; idx<6; idx++)
	{
		/* Needed powers of theta */
		theta = thetapoints[idx];
		theta8 = pow(theta,8);
		theta9 = theta*theta8;
		theta10 = theta*theta9;
		theta11 = theta*theta10;
		theta12 = theta*theta11;
		theta13 = theta*theta12;

		T3offset = IMRPhenomTTaylorT3(theta, pPhase); // TaylorT3 value at the specific collocation point time */

		gsl_vector_set(b,idx,(4/theta/theta/theta)*(omegainsppoints[idx] - T3offset));

		gsl_matrix_set(A,idx,0,theta8);
		gsl_matrix_set(A,idx,1,theta9);
		gsl_matrix_set(A,idx,2,theta10);
		gsl_matrix_set(A,idx,3,theta11);
		gsl_matrix_set(A,idx,4,theta12);
		gsl_matrix_set(A,idx,5,theta13);

	}

	/* We now solve the system A x = b via an LU decomposition */
	gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,x);

	/* Set inspiral phenomenological coefficients from solution to A x = b */
	pPhase->omegaInspC1 = gsl_vector_get(x,0);
	pPhase->omegaInspC2 = gsl_vector_get(x,1);
	pPhase->omegaInspC3 = gsl_vector_get(x,2);
	pPhase->omegaInspC4 = gsl_vector_get(x,3);
	pPhase->omegaInspC5 = gsl_vector_get(x,4);
	pPhase->omegaInspC6 = gsl_vector_get(x,5);

	/* Boundary time between late inspiral and merger regions.
	This is selected to correspond to theta=0.81, earlier than the last collocation point at theta=0.82, because in this way
	the derivative at the boundary with the merger region, needed for the merger reconstruction, is less forced. */
    REAL8 tCut = -5.0/(eta*pow(0.81,8));
    pPhase->tCut22 = tCut;
	REAL8 omegaCut =  IMRPhenomTInspiralOmegaAnsatz22(0.81, pPhase);// 22 frequency value at tCut

	/* Tidy up in preparation for next GSL solve ... */
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_matrix_free(A);
	gsl_permutation_free(p);

	/*Set linear system for solving merger coefficients*/
	p = gsl_permutation_alloc(3);
	b = gsl_vector_alloc(3);
	x = gsl_vector_alloc(3);
	A = gsl_matrix_alloc(3,3);

	/* ********************************** */
	/* *** MERGER COEFFICIENTS ********** */
	/* ********************************** */

	/* In order to obtain the value of the 3 free coefficients of the ansatz defined in IMRPhenomTMergerOmegaAnsatz22, we need to solve a linear system
      where the ansatz is imposed to match a collocation point value at theta(t)=0.95 and to satisfy continuity with the inspiral ansatz and differentiability with the ringdown ansatz.
      Notice that the ansatz definition IMRPhenomTMergerOmegaAnsatz22 differs from eq [27-26] of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf)
      in that here 2 of the coefficients are analytically solved for imposing continuity with the ringdown ansatz and differentiability with the inspiral ansatz.

      Reminder: the ansatz described in IMRPhenomTMergerOmegaAnsatz22 corresponds to the rescaled frequency \bar{\omega}=1 - (\omega / \omega_ring) (eq. 27 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf),
      so the system is solved for the rescaled frequency. For the analytical phase of IMRPhenomTMergerPhaseAnsatz22, which is the quantity employed
      in the waveform construction, the factors are already included to produce the correct phase. */


	REAL8 tMerger22 = -5.0/(eta*pow(0.95,8)); // Collocation point time of the merger region. This is placed at a fixed position theta=0.95.
	REAL8 omegaMergerCP = 1. - IMRPhenomT_Merger_Freq_CP1_22(eta, S, dchi, delta)/pPhase->omegaRING; // Collocation point.
	REAL8 omegaCutBar = 1. - omegaCut/pPhase->omegaRING; // Boundary frequency between inspiral and merger ringdown, rescaled.

	/* Now we need the derivative values at the boundaries.
	   Being the inspiral and ringdown ansatz differentiable analytical functions, the analytical derivative could be computed.
	   However, for code saving, it is enough to compute a numerical derivative at the boundary. Since the expressions are differentiable,
	   the numerical derivative is clean from any noise, and using a first order finite difference scheme is enough to obtain the derivative
	   at sufficient precission. */

	REAL8 theta2 = pow(-eta*tCut/5.,-1./8);
	REAL8 theta1 = pow(-eta*(tCut-0.0000001)/5,-1./8);
	REAL8 domegaCut = -(IMRPhenomTInspiralOmegaAnsatz22(theta2, pPhase) - IMRPhenomTInspiralOmegaAnsatz22(theta1, pPhase))/(0.0000001)/pPhase->omegaRING; // Derivative of rescale frequency at the inspiral boundary
	REAL8 domegaPeak = -(IMRPhenomTRDOmegaAnsatz22(0.0000001, pPhase) - IMRPhenomTRDOmegaAnsatz22(0., pPhase))/(0.0000001)/pPhase->omegaRING; // Derivative of rescale frequency at the ringdown boundary
	pPhase->domegaPeak = domegaPeak;

	/* Now we set the linear system for solving merger ansatz coefficients of IMRPhenomTMergerOmegaAnsatz22,
	as explained in eq. 30 of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf)
	Reminder: Here two coefficients are already analytically solved, so the system contains only three equations for the three remaining coefficients. */

	/* A_{0,i} and b_{0}. Here we impose continuity with the inspiral region.
	 The value of the solution vector coefficient b_{0} is the rescaled frequency at the inspiral boundary, minus the terms of the ansatz which are already fixed analytically.
	 The value of the first row of basis matrix are the arcsinh powers of the ansatz evaluated at the boundary time tCut. */

	gsl_complex phi = gsl_complex_rect(pPhase->alpha1RD*tCut,0);
	REAL8 ascut = GSL_REAL(gsl_complex_arcsinh(phi));
	REAL8 ascut2 = ascut*ascut;
	REAL8 ascut3 = ascut*ascut2;
	REAL8 ascut4 = ascut*ascut3;

	gsl_vector_set(b,0,omegaCutBar - (1. - pPhase->omegaPeak/pPhase->omegaRING) - (pPhase->domegaPeak/pPhase->alpha1RD)*ascut);

	gsl_matrix_set(A,0,0,ascut2);
	gsl_matrix_set(A,0,1,ascut3);
	gsl_matrix_set(A,0,2,ascut4);

	/* A_{1,i} and b_{1}. Here we impose the collocation point at theta=0.95.
	 The value of the solution vector coefficient b_{1} is the rescaled value of the collocation point, minus the terms of the ansatz which are already fixed analytically.
	 The value of the second row of basis matrix are the arcsinh powers of the ansatz evaluated at the collocation point time tMerger22 (determined from theta(t)=0.95). */

	phi = gsl_complex_rect(pPhase->alpha1RD*tMerger22,0);
	REAL8 as025cut = GSL_REAL(gsl_complex_arcsinh(phi));
	REAL8 as025cut2 = as025cut*as025cut;
	REAL8 as025cut3 = as025cut*as025cut2;
	REAL8 as025cut4 = as025cut*as025cut3;

	gsl_vector_set(b,1,omegaMergerCP - (1. - pPhase->omegaPeak/pPhase->omegaRING) - (pPhase->domegaPeak/pPhase->alpha1RD)*as025cut);

	gsl_matrix_set(A,1,0,as025cut2);
	gsl_matrix_set(A,1,1,as025cut3);
	gsl_matrix_set(A,1,2,as025cut4);

	/* A_{2,i} and b_{2}. Here we impose differentiability with the ringdown region.
	 The value of the solution vector coefficient b_{2} is the rescaled frequency derivative at the ringdown boundary, minus the terms of the ansatz derivative which are already fixed analytically.
	 The value of the third row of basis matrix are the arcsinh powers of the ansatz evaluated at the boundary time tCut. */

	REAL8 dencut = sqrt(1.0 + tCut*tCut*pPhase->alpha1RD*pPhase->alpha1RD); // Factor that appears from the derivative of the ansatz

	gsl_matrix_set(A,2,0,2.0*pPhase->alpha1RD*ascut/dencut);
	gsl_matrix_set(A,2,1,3.0*pPhase->alpha1RD*ascut2/dencut);
	gsl_matrix_set(A,2,2,4.0*pPhase->alpha1RD*ascut3/dencut);

	gsl_vector_set(b,2,domegaCut - pPhase->domegaPeak/dencut);

	/* We now solve the system A x = b via an LU decomposition */
	gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,x);

	/* Set merger phenomenological coefficients from solution to A x = b */
	pPhase->omegaMergerC1 = gsl_vector_get(x,0);
	pPhase->omegaMergerC2 = gsl_vector_get(x,1);
	pPhase->omegaMergerC3 = gsl_vector_get(x,2);

	// Free the gsl objects employed in solving the coefficient system
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_matrix_free(A);
	gsl_permutation_free(p);

	/* *********************************************** */
	/* *** PHASE CONTINUITY BETWEEN REGIONS ********** */
	/* *********************************************** */

	/* Phase of the different regions correspond to the analytical integration of the frequency ansatz. 
	Here we adapt the corresponding offsets for each phase ansatz in order to obtain continuity between regions */

	/* First we set the offsets to zero in order to call the phase ansatzs without the offsets */
	pPhase->phOffInsp = 0.;
	pPhase->phOffMerger = 0.;
	pPhase->phOffRD = 0.;

	/* phOffInsp contains the difference between the early and late inspiral regions at the early inspiral boundary defined by tEarly (theta=0.33).
	This is only needed if reconstruction is non-default (4 regions) but it is harmless if reconstruction employs 3 regions */
	REAL8 thetabarini = pow(eta*(tt0 - tEarly),-1./8);
	REAL8 thetabarini2 = pow(-eta*(tEarly),-1./8);
	pPhase->phOffInsp = IMRPhenomTInspiralPhaseTaylorT3(thetabarini, wf, pPhase) - IMRPhenomTInspiralPhaseAnsatz22(tEarly, thetabarini2, wf, pPhase);

	REAL8 thetabarCut = pow(-eta*tCut,-1./8);
	REAL8 phMECOinsp = IMRPhenomTInspiralPhaseAnsatz22(tCut, thetabarCut, wf, pPhase); //Value of inspiral phase at inspiral-merger boundary
	REAL8 phMECOmerger = IMRPhenomTMergerPhaseAnsatz22(tCut, pPhase); //Value of merger phase at merger-ringdown boundary

	pPhase->phOffMerger = phMECOinsp - phMECOmerger; // Needed offset for merger ansatz in order to match inspiral phase value at the boundary
	pPhase->phOffRD = IMRPhenomTMergerPhaseAnsatz22(0, pPhase); //Neded offset for ringdown ansatz in order to match merger phase value at the boundary 


    /* *********************************************** */
	/* *** EPOCH AND WAVEFORM LENGHT ***************** */
	/* *********************************************** */


	/* Waveform lenght is determined by the time spent from the starting frequency to the peak amplitude time of the 22 mode, with the addition of 500M after that time for 
	having a sufficient ringdown window for sane Fourier tranforms */

	/* Starting time is determined as the time at which starting frequency occurs, and it is determined by root finding employing the frequency function */


	/* First we set the gsl objects needed for root finding the minimum time */
    const gsl_root_fsolver_type *solver_type;
    gsl_root_fsolver *solver;
    gsl_function F;

    /* Populate an auxiliary struct needed for gsl root finder */
    struct FindRootParams frparams;
    frparams.f0 = pPhase->Mfmin;
    frparams.pPhase = pPhase;
    frparams.wf = wf;

    /* Populate the gsl function needed for the root finder */
    F.function = &GetTimeOfFreq;
    F.params = &frparams;

    /* Allocate a brent solver and set it to use F */
    /* Left boundary of the solver is set at a sufficient early time such that above one solar mass, it exists a solution. */
    solver_type = gsl_root_fsolver_brent;
    gsl_set_error_handler_off();

    int status;
    int status_solver = GSL_CONTINUE;
    solver = gsl_root_fsolver_alloc(solver_type);
    status = gsl_root_fsolver_set(solver, &F, -1000000000.0,tEnd);

    double r = 0.0;
    for (UINT4 i = 1; i <= 1000 && status_solver == GSL_CONTINUE; ++i) {
        /* iterate one step of the solver */
        status_solver = gsl_root_fsolver_iterate(solver);
        if (status_solver != GSL_SUCCESS)
            break;

        /* get the solver's current best solution and bounds */
        r = gsl_root_fsolver_root(solver);
        double x_lo = gsl_root_fsolver_x_lower(solver);
        double x_hi = gsl_root_fsolver_x_upper(solver);

        /* Check to see if the solution is within 0.0001 */
        status_solver = gsl_root_test_interval(x_lo, x_hi, 0, 0.0001);
        }
    XLAL_CHECK(GSL_SUCCESS == status, XLAL_EFUNC, "Error: Root finder unable to solve minimum time. Minimum frequency may be too high for these parameters. Try reducing fmin below the peak frequency: %.8f Hz. \n", pPhase->omegaPeak/(LAL_TWOPI*wf->M_sec));

    pPhase->tmin = r; // Minimum time is selected as the solution of the above root finding operation
    gsl_root_fsolver_free(solver); // Free the gsl solver

    /* Now we repeat the same procedure for determining the time at which the specified reference frequency occurs. This is needed to set the reference phase */

    /* First we check if fmin and fref coincide, to save the computation */
    if(wf->fRef == wf->fmin)
    {
    	pPhase->tRef = pPhase->tmin;
    }
    /* If not, we repeat the same operation to obtain the reference time */
    else
    {
    	frparams.f0 = pPhase->MfRef;
    	frparams.pPhase = pPhase;
    	frparams.wf = wf;

    	F.function = &GetTimeOfFreq;
    	F.params = &frparams;

    	/* Allocate a brent solver and set it to use F */
    	solver_type = gsl_root_fsolver_brent;
    	solver = gsl_root_fsolver_alloc(solver_type);
    	status = gsl_root_fsolver_set(solver, &F, -1000000000.0,tEnd);
    	status_solver = GSL_CONTINUE;

    	for (UINT4 i = 1; i <= 1000 && status_solver == GSL_CONTINUE; ++i) {
        	/* iterate one step of the solver */
        	status_solver = gsl_root_fsolver_iterate(solver);
        	if (status_solver != GSL_SUCCESS)
            	break;

        	/* get the solver's current best solution and bounds */
        	r = gsl_root_fsolver_root(solver);
        	double x_lo = gsl_root_fsolver_x_lower(solver);
        	double x_hi = gsl_root_fsolver_x_upper(solver);

        	/* Check to see if the solution is within 0.001 */
        	status_solver = gsl_root_test_interval(x_lo, x_hi, 0, 0.0001);
        	}
        	XLAL_CHECK(GSL_SUCCESS == status, XLAL_EFUNC, "Error: Root finder unable to solve reference time. Reference frequency may be too high for these parameters. Try reducing fRef below the peak frequency: %.8f Hz. \n", pPhase->omegaPeak/(LAL_TWOPI*wf->M_sec));

    	pPhase->tRef = r;
    	gsl_root_fsolver_free(solver);
    }

    pPhase->tminSec = wf->M_sec*pPhase->tmin;

    /* Required waveform length. We select maximum time of 500M since it is enough to contain the non-negligible ringdown content of the modes up to m=5.
    Length of early and late inspiral regions are computed since in the case of non-default reconstruction (4 regions) this is needed to compute frequencies
    in both regimes with the two different TaylorT3 implementations. If default reconstruction, it is harmless. */

    pPhase->wflength = floor((tEnd - pPhase->tmin)/wf->dtM);
    if(tEarly<=pPhase->tmin && tCut>pPhase->tmin)
    {
    	pPhase->wflength_insp_late = floor((tCut - pPhase->tmin + wf->dtM)/wf->dtM);
    	pPhase->wflength_insp_early = 0;
    }
    else if(tCut<=pPhase->tmin)
    {
    	pPhase->wflength_insp_late = 0;
    	pPhase->wflength_insp_early = 0;
    }
    else{
    	pPhase->wflength_insp_late = floor((tCut - tEarly + wf->dtM)/wf->dtM);
    	pPhase->wflength_insp_early = floor((tEarly - pPhase->tmin + wf->dtM)/wf->dtM);
    }
    


	return XLAL_SUCCESS;
}

/* ---------------------------------------------------------------------------------------------------- */

/* **************************************************** */
/* ******* HIGHER MODES STRUCTURE INITIALIZATION ****** */
/* **************************************************** */

/* Function for populating the IMRPhenomTHMAmpStruct.
   This will solve the linear systems for obtaining ansatz coefficients from lm amplitude collocation points.
   PN Coefficients of the mode amplitude are defined and stored for each mode.
   Inspiral PN amplitude is a complex quantity, with its argument contributing to the mode's phase. Since  amplitude collocation points are taken from the absolute value,
   some preprocessing of the fits is needed, in particular to obtain the real projection once known the argument angle. */

int IMRPhenomTSetHMAmplitudeCoefficients(int l, int m, IMRPhenomTHMAmpStruct *pAmp, IMRPhenomTPhase22Struct *pPhase, IMRPhenomTWaveformStruct *wf)
{
	REAL8 eta   = wf->eta;   // Symmetric mass ratio
	REAL8 chi1 = wf->chi1L; // Dimensionless aligned spin of companion 1
	REAL8 chi2 = wf->chi2L; // Dimensionless aligned spin of companion 2
	REAL8 S     = wf->Shat;  // Dimensionless effective spin parameters employed for fit evaluation
	REAL8 dchi  = wf->dchi;  // Dimensionless spin difference chi1 - chi2
	REAL8 delta = wf->delta; // Mass asymmetry parameter
	REAL8 tCut  = tCUT_Amp;  // tCUT_Amp = -150

	pAmp->fac0 = 2*eta*sqrt(16*LAL_PI/5); // Precomputed amplitude global factor (see eq. 12-14 of PhenomTHM paper https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf)

	/* GSL objects for solving system of equations via LU decomposition */
	gsl_vector *b, *x;
	gsl_matrix *A;
	gsl_permutation *p;
	int s;

	/* For the selected mode, all needed fit values, PN coefficients and ringdown coefficients to fully determine the amplitude ansatz are computed. */

	REAL8 ampInspCP[3]; // Inspiral collocation point values.
	REAL8 ampMergerCP1, ampPeak;    // Merger collocation point value and peak amplitude
	REAL8 ampRDC3, fDAMP, fDAMPn2, fDAMP_prec, fDAMPn2_prec, coshc3, tanhc3; // Ringdown ansatz c3 free coefficient, damping frequencies and needed quantities to compute ringdown ansatz coefficients.
	/* Ringdown ansatz coefficients c_{1,2,3,4} as defined in eq. [6-8] of Damour&Nagar 2014 (https://arxiv.org/pdf/1406.0401.pdf), also explained in eq.26 of THM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf)
	 and alpha_1, alpha_2, alpha_21 as defined in Sec IV of the same reference are stored in the amplitude struct since they are directly passed to the ringdown ansatz.
	 tshift, also passed to the amplitude struct, corresponds to the peak amplitude time of the (l,m) mode (0 for the 22 by construction). */
	gsl_complex phi; // Argument for the gsl hyperbolic functions. It is a real quantity but gsl functions expect a complex number, so it is constructed as (phi,0).

	/* The PN amplitude coefficients of the inspiral ansatz are defined in Appendix A of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf).
	   They are constructed from:
	   - 3PN expressions from Blanchet et al. 2008 (Class.Quant.Grav.25:165003,2008, https://arxiv.org/abs/0802.1249)
	   - 2PN spin corrections from Buonanno et al. 2013 (Phys. Rev D87, 044009, 2013, https://arxiv.org/abs/1209.6349)
	   - 1.5PN contributions from Arun et al. 2008 (Phys.Rev.D79:104023,2009; Erratum-ibid.D84:049901,2011, https://arxiv.org/abs/0810.5336) */

	if (l==2 && m==2)
	{
		/* Needed fits for collocation points and ringdown coefficient */
		ampInspCP[0] = IMRPhenomT_Inspiral_Amp_CP1_22(eta, S, dchi, delta);
		ampInspCP[1] = IMRPhenomT_Inspiral_Amp_CP2_22(eta, S, dchi, delta);
		ampInspCP[2] = IMRPhenomT_Inspiral_Amp_CP3_22(eta, S, dchi, delta);
		ampMergerCP1 = IMRPhenomT_Merger_Amp_CP1_22(eta, S, dchi, delta);
		ampPeak = IMRPhenomT_PeakAmp_22(eta, S, dchi, delta);
		ampRDC3 = IMRPhenomT_RD_Amp_C3_22(eta, S);

		fDAMP     = evaluate_QNMfit_fdamp22(wf->afinal) / (wf->Mfinal); //damping frequency of 122 QNM
		fDAMPn2   = evaluate_QNMfit_fdamp22n2(wf->afinal) / (wf->Mfinal); //damping frequency of 222 QNM

		fDAMP_prec = evaluate_QNMfit_fdamp22(wf->afinal_prec) / (wf->Mfinal);
		fDAMPn2_prec   = evaluate_QNMfit_fdamp22n2(wf->afinal_prec) / (wf->Mfinal);

		pAmp->tshift = 0.0; //Peak time, by construction 0 for l=2, m=2

		/* PN Amplitude coefficients, defined in eq.A1 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf (before global rotation applied, see eq.13 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf) */

		pAmp->ampN = 1.0;
		pAmp->amp0halfPNreal = 0.0;
		pAmp->amp0halfPNimag = 0.0;

		pAmp->amp1PNreal = -2.5476190476190474 + (55*eta)/42.;
		pAmp->amp1PNimag = 0.0;

		pAmp->amp1halfPNreal = (-2*chi1)/3. - (2*chi2)/3. - (2*chi1*delta)/(3.*((1 - delta)/2. + (1 + delta)/2.)) + (2*chi2*delta)/(3.*((1 - delta)/2. + (1 + delta)/2.)) + (2*chi1*eta)/3. + (2*chi2*eta)/3. + 2*LAL_PI;
    	pAmp->amp1halfPNimag = 0.0;

    	pAmp->amp2PNreal = -1.437169312169312 + pow(chi1,2)/2. + pow(chi2,2)/2. + (pow(chi1,2)*delta)/2. - (pow(chi2,2)*delta)/2. - (1069*eta)/216. - pow(chi1,2)*eta + 2*chi1*chi2*eta - pow(chi2,2)*eta + (2047*pow(eta,2))/1512.;
    	pAmp->amp2PNimag = 0.0;

    	pAmp->amp2halfPNreal = - (107*LAL_PI)/21. + (34*eta*LAL_PI)/21.;
    	pAmp->amp2halfPNimag =  -24.*eta;

    	pAmp->amp3PNreal = 41.78634662956092 - (278185*eta)/33264. - (20261*pow(eta,2))/2772. + (114635*pow(eta,3))/99792. - (856*LAL_GAMMA)/105. + (2*pow(LAL_PI,2))/3. + (41*eta*pow(LAL_PI,2))/96.;
    	pAmp->amp3PNimag = (428./105)*LAL_PI;

    	pAmp->amp3halfPNreal = (-2173*LAL_PI)/756. - (2495*eta*LAL_PI)/378. + (40*pow(eta,2)*LAL_PI)/27.;
    	pAmp->amp3halfPNimag = (14333*eta)/162. - (4066*pow(eta,2))/945.;
 
    	pAmp->amplog = -428/105.;
	}

	else if (l==2 && m==1)
	{
		/* Needed fits for collocation points and ringdown coefficient */
		ampInspCP[0] = IMRPhenomT_Inspiral_Amp_CP1_21(eta, S, dchi, delta);
		ampInspCP[1] = IMRPhenomT_Inspiral_Amp_CP2_21(eta, S, dchi, delta);
		ampInspCP[2] = IMRPhenomT_Inspiral_Amp_CP3_21(eta, S, dchi, delta);
		ampMergerCP1 = IMRPhenomT_Merger_Amp_CP1_21(eta, S, dchi, delta);
		ampPeak = IMRPhenomT_PeakAmp_21(eta, S, dchi, delta);
		ampRDC3 = IMRPhenomT_RD_Amp_C3_21(eta, S, dchi);

		fDAMP     = evaluate_QNMfit_fdamp21(wf->afinal) / (wf->Mfinal); //damping frequency of 121 QNM
		fDAMPn2   = evaluate_QNMfit_fdamp21n2(wf->afinal) / (wf->Mfinal); //damping frequency of 221 QNM

		fDAMP_prec = evaluate_QNMfit_fdamp21(wf->afinal_prec) / (wf->Mfinal);
		fDAMPn2_prec   = evaluate_QNMfit_fdamp21n2(wf->afinal_prec) / (wf->Mfinal);

		pAmp->tshift = IMRPhenomT_tshift_21(eta, S, dchi); //Peak time of l=2, m=1 mode

		/* PN Amplitude coefficients, defined in eq.A2 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf (before global rotation applied, see eq.13 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf) */

		pAmp->ampN = 0.0;
		pAmp->amp0halfPNreal = delta/3.;
		pAmp->amp0halfPNimag = 0.0;

		pAmp->amp1PNreal = -chi1/4. + chi2/4. - (chi1*delta)/(4.*((1 - delta)/2. + (1 + delta)/2.)) - (chi2*delta)/(4.*((1 - delta)/2. + (1 + delta)/2.));
		pAmp->amp1PNimag = 0.0;

		pAmp->amp1halfPNreal = -17*delta/84. + (5*delta*eta)/21.;
    	pAmp->amp1halfPNimag = 0.0;

    	pAmp->amp2PNimag = -delta/6. - (2*delta*log(2))/3.;
    	pAmp->amp2PNreal = (79*chi1)/84. - (79*chi2)/84. + (79*chi1*delta)/84. + (79*chi2*delta)/84. - (43*chi1)/(42.*((1 - delta)/2. + (1 + delta)/2.)) + (43*chi2)/(42.*((1 - delta)/2. + (1 + delta)/2.)) - (43*chi1*delta)/(42.*((1 - delta)/2. + (1 + delta)/2.)) - (43*chi2*delta)/(42.*((1 - delta)/2. + (1 + delta)/2.)) - (139*chi1*eta)/84. + (139*chi2*eta)/84. - (139*chi1*delta*eta)/84. - (139*chi2*delta*eta)/84. + (86*chi1*eta)/(21.*((1 - delta)/2. + (1 + delta)/2.)) - (86*chi2*eta)/(21.*((1 - delta)/2. + (1 + delta)/2.)) + (43*chi1*delta*eta)/(21.*((1 - delta)/2. + (1 + delta)/2.)) + (43*chi2*delta*eta)/(21.*((1 - delta)/2. + (1 + delta)/2.)) + (delta*LAL_PI)/3.;

    	pAmp->amp2halfPNimag = 0.0;
    	pAmp->amp2halfPNreal = (-43*delta)/378. - (509*delta*eta)/378. + (79*delta*pow(eta,2))/504.;

    	pAmp->amp3PNimag = -((-17*delta)/168. + (353*delta*eta)/84. - (17*delta*log(2))/42. + (delta*eta*log(2))/7.);
    	pAmp->amp3PNreal = (-17*delta*LAL_PI)/84. + (delta*eta*LAL_PI)/14.;

    	pAmp->amp3halfPNreal = 0.0;
    	pAmp->amp3halfPNimag = 0.0;

    	pAmp->amplog = 0.0;
	}

	else if (l==3 && m==3)
	{
		/* Needed fits for collocation points and ringdown coefficient */
		ampInspCP[0] = IMRPhenomT_Inspiral_Amp_CP1_33(eta, S, dchi, delta);
		ampInspCP[1] = IMRPhenomT_Inspiral_Amp_CP2_33(eta, S, dchi, delta);
		ampInspCP[2] = IMRPhenomT_Inspiral_Amp_CP3_33(eta, S, dchi, delta);
		ampMergerCP1 = IMRPhenomT_Merger_Amp_CP1_33(eta, S, dchi, delta);
		ampPeak = IMRPhenomT_PeakAmp_33(eta, S, dchi, delta);
		ampRDC3 = IMRPhenomT_RD_Amp_C3_33(eta, S);

		fDAMP     = evaluate_QNMfit_fdamp33(wf->afinal) / (wf->Mfinal); //damping frequency of 133 QNM
		fDAMPn2   = evaluate_QNMfit_fdamp33n2(wf->afinal) / (wf->Mfinal); //damping frequency of 233 QNM

		fDAMP_prec = evaluate_QNMfit_fdamp33(wf->afinal_prec) / (wf->Mfinal);
		fDAMPn2_prec   = evaluate_QNMfit_fdamp33n2(wf->afinal_prec) / (wf->Mfinal);

		pAmp->tshift = IMRPhenomT_tshift_33(eta, S);

		/* PN Amplitude coefficients, defined in eq.A3 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf (before global rotation applied, see eq.13 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf) */

		pAmp->ampN = 0.0;
		pAmp->amp0halfPNimag = 0.0;
		pAmp->amp0halfPNreal = 0.7763237542601484*delta;

		pAmp->amp1PNreal = 0.0;
		pAmp->amp1PNimag = 0.0;

		pAmp->amp1halfPNimag = 0.0;
    	pAmp->amp1halfPNreal = -3.1052950170405937*delta + 1.5526475085202969*delta*eta;

    	pAmp->amp2PNimag = -1.371926598204461*delta;
    	pAmp->amp2PNreal = -(-0.5822428156951114*chi1 + 0.5822428156951114*chi2 - 7.316679009572791*delta - 0.5822428156951114*chi1*delta - 0.5822428156951114*chi2*delta + (1.3585665699552598*chi1)/((1 - delta)/2. + (1 + delta)/2.) - (1.3585665699552598*chi2)/((1 - delta)/2. + (1 + delta)/2.) + (1.3585665699552598*chi1*delta)/((1 - delta)/2. + (1 + delta)/2.) + (1.3585665699552598*chi2*delta)/((1 - delta)/2. + (1 + delta)/2.) + 1.7467284470853341*chi1*eta - 1.7467284470853341*chi2*eta + 1.7467284470853341*chi1*delta*eta + 1.7467284470853341*chi2*delta*eta - (5.434266279821039*chi1*eta)/((1 - delta)/2. + (1 + delta)/2.) + (5.434266279821039*chi2*eta)/((1 - delta)/2. + (1 + delta)/2.) - (2.7171331399105196*chi1*delta*eta)/((1 - delta)/2. + (1 + delta)/2.) - (2.7171331399105196*chi2*delta*eta)/((1 - delta)/2. + (1 + delta)/2.));

    	pAmp->amp2halfPNimag = 0.0;
    	pAmp->amp2halfPNreal = -(-0.08680711070363478*delta + 8.647776123213047*delta*eta - 2.0866641516022777*delta*pow(eta,2));

    	pAmp->amp3PNreal = 0.0;
    	pAmp->amp3PNimag = 0.0;

    	pAmp->amp3halfPNreal = 0.0;
    	pAmp->amp3halfPNimag = 0.0;

    	pAmp->amplog = 0.0;
	}

	else if (l==4 && m==4)
	{
		/* Needed fits for collocation points and ringdown coefficient */
		ampInspCP[0] = IMRPhenomT_Inspiral_Amp_CP1_44(eta, S, dchi, delta);
		ampInspCP[1] = IMRPhenomT_Inspiral_Amp_CP2_44(eta, S, dchi, delta);
		ampInspCP[2] = IMRPhenomT_Inspiral_Amp_CP3_44(eta, S, dchi, delta);
		ampMergerCP1 = IMRPhenomT_Merger_Amp_CP1_44(eta, S, dchi, delta);
		ampPeak = IMRPhenomT_PeakAmp_44(eta, S, dchi, delta);
		ampRDC3 = IMRPhenomT_RD_Amp_C3_44(eta, S);

		fDAMP     = evaluate_QNMfit_fdamp44(wf->afinal) / (wf->Mfinal); //damping frequency of 144 QNM
		fDAMPn2   = evaluate_QNMfit_fdamp44n2(wf->afinal) / (wf->Mfinal); //damping frequency of 244 QNM

		fDAMP_prec = evaluate_QNMfit_fdamp44(wf->afinal_prec) / (wf->Mfinal);
		fDAMPn2_prec   = evaluate_QNMfit_fdamp44n2(wf->afinal_prec) / (wf->Mfinal);

		pAmp->tshift = IMRPhenomT_tshift_44(eta, S); //Peak time of l=4, m=4 mode

		/* PN Amplitude coefficients, defined in eq.A4 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf (before global rotation applied, see eq.13 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf) */

		pAmp->ampN = 0.0;
		pAmp->amp0halfPNreal = 0.0;
		pAmp->amp0halfPNimag = 0.0;

		pAmp->amp1PNreal = 0.751248226425348*(1 - 3*eta);
		pAmp->amp1PNimag = 0.0;

		pAmp->amp1halfPNreal = 0.0;
    	pAmp->amp1halfPNimag = 0.0;

    	pAmp->amp2PNreal = -4.049910893365739 + 14.489984730901032*eta - 5.9758381647470875*pow(eta,2);
    	pAmp->amp2PNimag = 0.0;

    	pAmp->amp2halfPNreal = 0.751248226425348*(4*LAL_PI - 12*eta*LAL_PI);
    	pAmp->amp2halfPNimag = 0.751248226425348*(-2.854822555520438 + 13.189467666561313*eta);

    	pAmp->amp3PNreal = -(-8*pow(0.7142857142857143,0.5)*(5.338016983016983 - (1088119*eta)/28600. + (146879*pow(eta,2))/2340. - (226097*pow(eta,3))/17160.))/9.;
    	pAmp->amp3PNimag = 0.0;

    	pAmp->amp3halfPNreal = 0.0;
    	pAmp->amp3halfPNimag = 0.0;

    	pAmp->amplog = 0.0;
	}

	else if (l==5 && m==5)
	{
		/* Needed fits for collocation points and ringdown coefficient */
		ampInspCP[0] = IMRPhenomT_Inspiral_Amp_CP1_55(eta, S, dchi, delta);
		ampInspCP[1] = IMRPhenomT_Inspiral_Amp_CP2_55(eta, S, dchi, delta);
		ampInspCP[2] = IMRPhenomT_Inspiral_Amp_CP3_55(eta, S, dchi, delta);
		ampMergerCP1 = IMRPhenomT_Merger_Amp_CP1_55(eta, S, dchi, delta);
		ampPeak = IMRPhenomT_PeakAmp_55(eta, S, dchi, delta);
		ampRDC3 = IMRPhenomT_RD_Amp_C3_55(eta, S, dchi);

		fDAMP     = evaluate_QNMfit_fdamp55(wf->afinal) / (wf->Mfinal); //damping frequency of 155 QNM
		fDAMPn2   = evaluate_QNMfit_fdamp55n2(wf->afinal) / (wf->Mfinal); //damping frequency of 255 QNM

		fDAMP_prec = evaluate_QNMfit_fdamp55(wf->afinal_prec) / (wf->Mfinal);
		fDAMPn2_prec   = evaluate_QNMfit_fdamp55n2(wf->afinal_prec) / (wf->Mfinal);

		pAmp->tshift = IMRPhenomT_tshift_55(eta, S); //Peak time of l=5, m=5 mode

		/* PN Amplitude coefficients, defined in eq.A5 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf (before global rotation applied, see eq.13 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf) */

		pAmp->ampN = 0.0;
		pAmp->amp0halfPNreal = 0.0;
		pAmp->amp0halfPNimag = 0.0;

		pAmp->amp1PNreal = 0.0;
		pAmp->amp1PNimag = 0.0;

		pAmp->amp1halfPNimag = 0.0;
    	pAmp->amp1halfPNreal = 0.8013768943966973*delta*(1 - 2*eta);

    	pAmp->amp2PNreal = 0.0;
    	pAmp->amp2PNimag = 0.0;

    	pAmp->amp2halfPNimag = 0.0;
    	pAmp->amp2halfPNreal = 0.8013768943966973*delta*(-6.743589743589744 + (688*eta)/39. - (256*pow(eta,2))/39.);

    	pAmp->amp3PNimag = -3.0177162096765713*delta + 12.454250695829877*delta*eta;
    	pAmp->amp3PNreal = 12.58799882096634*delta - 25.175997641932675*delta*eta;

    	pAmp->amp3halfPNreal = 0.0;
    	pAmp->amp3halfPNimag = 0.0;

    	pAmp->amplog = 0.0;
	}

	else
	{
		XLAL_ERROR(XLAL_EFUNC, "Mode not implemented in PhenomTHM. Modes available: [2,2], [2,1], [3,3], [4,4], [5,5].");
	}

	/***********************************************************/
	/************** RINGDOWN ANSATZ COEFFICIENTS ***************/
	/***********************************************************/

	/* Ringdown ansatz coefficients as defined in in eq. [6-8] of Damour&Nagar 2014 (https://arxiv.org/pdf/1406.0401.pdf).
	See also eq.26c-e of THM paper: https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf
	Essentially, c3 is the only calibrated coefficient. c2 accounts for the effect of the first overtone in the early ringdown, and c1 and c4 fix the peak amplitude and null derivative at peak */
	
	pAmp->alpha1RD = 2*LAL_PI*fDAMP;
	pAmp->alpha2RD = 2*LAL_PI*fDAMPn2;
	pAmp->alpha21RD = 0.5*(pAmp->alpha2RD - pAmp->alpha1RD); //Coefficient c_2 of ringdown amplitude ansatz as defined in equation 7 of Damour&Nagar 2014 (https://arxiv.org/pdf/1406.0401.pdf) and eq.26d of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf

	pAmp->alpha1RD_prec = 2*LAL_PI*fDAMP_prec;
	pAmp->alpha2RD_prec = 2*LAL_PI*fDAMPn2_prec;
	pAmp->alpha21RD_prec = 0.5*(pAmp->alpha2RD_prec - pAmp->alpha1RD_prec);

	pAmp->c3 = ampRDC3;
	pAmp->c2 = 0.5*(pAmp->alpha2RD - pAmp->alpha1RD);
	pAmp->c2_prec = 0.5*(pAmp->alpha2RD_prec - pAmp->alpha1RD_prec);

	phi = gsl_complex_rect(pAmp->c3,0); // Needed complex parameter for gsl hyperbolic functions
	coshc3 =  GSL_REAL(gsl_complex_cosh(phi));
	tanhc3 =  GSL_REAL(gsl_complex_tanh(phi));

	/* This condition ensures that the second derivative of the amplitude at the mode peak is always zero or negative, not producing then a second peak in the ringdown */
	if(fabs(pAmp->c2) > fabs(0.5*pAmp->alpha1RD/tanhc3))
	{
		pAmp->c2 = -0.5*pAmp->alpha1RD/tanhc3;
	}
	if(fabs(pAmp->c2_prec) > fabs(0.5*pAmp->alpha1RD_prec/tanhc3))
	{
		pAmp->c2_prec = -0.5*pAmp->alpha1RD_prec/tanhc3;
	}
	

	pAmp->c1 = ampPeak*pAmp->alpha1RD*coshc3*coshc3/pAmp->c2;
	pAmp->c1_prec = ampPeak*pAmp->alpha1RD_prec*coshc3*coshc3/pAmp->c2_prec;
	pAmp->c4 = ampPeak - pAmp->c1*tanhc3;
	pAmp->c4_prec = ampPeak - pAmp->c1_prec*tanhc3;

	/***********************************************************/
	/************** INSPIRAL COEFFICIENTS SOLUTION *************/
	/***********************************************************/

	/* In order to obtain the value of the 3 unknown extra coefficients of the inspiral amplitude ansatz, we need to solve a linear system
      where the ansatz is imposed to match collocation point values at the corresponding collocation point times.
      See equation 15 of THM paper: https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */

	/* Initialise the extra coefficients to zero, so calls to the amplitude ansatz function returns pure PN */
	pAmp->inspC1 = 0.0;
	pAmp->inspC2 = 0.0;
	pAmp->inspC3 = 0.0;

	REAL8 tinsppoints[3]     = {-2000., -250., -150.0}; // Collocation point times as defined in Eq. 16 of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf)

	/* We allocate a rank three linear system to solve the coefficients */
	p = gsl_permutation_alloc(3);
	b = gsl_vector_alloc(3);
	x = gsl_vector_alloc(3);
	A = gsl_matrix_alloc(3,3);

	REAL8 theta, omega, xx, x4, x4half, x5; // Needed powers of PN parameter x=v^2=(\omega_orb)^(2/3)=(0.5\omega_22)^(2/3)
	REAL8 ampoffset; // Known PN part of the amplitude
	REAL8 bi; // CP value - known PN amplitude vector

	/* In this loop over collocation points, the components of the solution vector b and the basis matrix A are established.
	See equation 15 of THM paper: https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */
	for (UINT4 idx=0; idx<3; idx++)
	{
		theta = pow(-eta*tinsppoints[idx]/5,-1./8);
		omega = IMRPhenomTomega22(tinsppoints[idx],theta, wf, pPhase); // Twice the orbital frequency at the collocation point time
		xx = pow(0.5*omega,2./3); // PN expansion parameter of the amplitude
		x4 = xx*xx*xx*xx; // Needed powers
		x4half = x4*sqrt(xx);
		x5 = x4*xx;

		ampoffset = creal(IMRPhenomTInspiralAmpAnsatzHM(xx, pAmp)); // Real part of the known PN contribution
		bi = (1./pAmp->fac0/xx)*(ampInspCP[idx] - ampoffset); // Solution vector: collocation point value minus the know PN part of the ansatz, factored by the amplitude factor to not include it in each basis function (the powers of x)

		gsl_vector_set(b,idx,bi); // Set b vector
		
		/*Set basis matrix elements, Basis functions are the higher order powers of x that we add to the PN ansatz */
		gsl_matrix_set(A,idx,0,x4);
		gsl_matrix_set(A,idx,1,x4half);
		gsl_matrix_set(A,idx,2,x5);
	}

	/* We now solve the system A x = b via an LU decomposition */
	gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,x);

	/* Set the extra pseudo-PN coefficients with the solutions of the system */
	pAmp->inspC1 = gsl_vector_get(x,0);
	pAmp->inspC2 = gsl_vector_get(x,1);
	pAmp->inspC3 = gsl_vector_get(x,2);

	/* Free the gsl solver */
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_matrix_free(A);
	gsl_permutation_free(p);

	/***********************************************************/
	/************** MERGER COEFFICIENTS SOLUTION ***************/
	/***********************************************************/


	/* We need to solve for 4 unknown coefficients of the ansatz defined in IMRPhenomTMergerAmpAnsatzHM.
	   Three of them are obtained by imposing continuity and differentiability at the boundary with the inspiral
	   region (i.e in t = tCUT_Amp) and continuity with the ringdown region (i.e imposing amp(t_peak)=AmpPeak).
	   Differentiability at the ringdown boundary is satisfied by default since both the merger and the ringdown
	   ansatz describe a peak (derivative zero) at t_peak.
	   4th coefficient is obtained by imposing the ansatz to match a collocation point at t=-25M.
	   See equation 31 of THM paper: https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */

	/* Reallocate gsl linear system, this time rank 4*/
	p = gsl_permutation_alloc(4);
	b = gsl_vector_alloc(4);
	x = gsl_vector_alloc(4);
	A = gsl_matrix_alloc(4,4);

	/* Set A_{0,i} and b_{0}: here we impose continuity with the inspiral region, essentially equating the value of the merger ansatz at the inspiral boundary time tCut
	with the value of the inspiral amplitude at that time. */

	theta = pow(-eta*tCut/5,-1./8);
	xx = pow(0.5*IMRPhenomTomega22(tCut,theta, wf, pPhase),2./3); // PN expansion parameter at tCut
	REAL8 ampinsp = copysign(1.0,creal(IMRPhenomTInspiralAmpAnsatzHM(xx, pAmp)))*cabs(IMRPhenomTInspiralAmpAnsatzHM(xx, pAmp)); // Value of the absolute inspiral amplitude, carrying the sign
	gsl_vector_set(b,0,ampinsp); // Set solution vector: Continuity with the inspiral region

	/* Here we compute the needed hyperbolic secant functions for the merger ansatz basis matrix */
	/* Time parameterisation of merger ansatz is in tau=t-tshift, so peak occurs at tau=0 */ 
	phi = gsl_complex_rect(pAmp->alpha1RD*(tCut-pAmp->tshift),0);
	gsl_complex phi2 = gsl_complex_rect(2*pAmp->alpha1RD*(tCut-pAmp->tshift),0);
	REAL8 sech1 = GSL_REAL(gsl_complex_sech(phi));
	REAL8 sech2 = GSL_REAL(gsl_complex_sech(phi2));

	/* Set the first row of the basis matrix. Just the functions that multiply each unknown coefficient of the ansatz */
	gsl_matrix_set(A,0,0,1.0);
	gsl_matrix_set(A,0,1,sech1);
	gsl_matrix_set(A,0,2,pow(sech2,1./7));
	gsl_matrix_set(A,0,3,(tCut-pAmp->tshift)*(tCut-pAmp->tshift));


	/* Set A_{1,i} and b_{1}: here we impose a collocation point value, essentially equating the value of the merger ansatz at the collocation point time
	with the value of the collocation point. */

	gsl_vector_set(b,1,ampMergerCP1); // Imposing collocation point value

	/* Here we compute the needed hyperbolic secant functions for the merger ansatz basis matrix */
	phi = gsl_complex_rect(pAmp->alpha1RD*(tcpMerger-pAmp->tshift),0);
	phi2 = gsl_complex_rect(2*pAmp->alpha1RD*(tcpMerger-pAmp->tshift),0);
	sech1 = GSL_REAL(gsl_complex_sech(phi));
	sech2 = GSL_REAL(gsl_complex_sech(phi2));

	/* Set the second row of the basis matrix. Just the functions that multiply each unknown coefficient of the ansatz */
	gsl_matrix_set(A,1,0,1.0);
	gsl_matrix_set(A,1,1,sech1);
	gsl_matrix_set(A,1,2,pow(sech2,1./7));
	gsl_matrix_set(A,1,3,(tcpMerger-pAmp->tshift)*(tcpMerger-pAmp->tshift));

	/* Set A_{2,i} and b_{2}: here we impose the peak amplitude value, essentially equating the value of the merger ansatz at the peak time tshift with
	with the value of the peak amplitude. */

	gsl_vector_set(b,2,ampPeak); // Imposing peak amplitude, that guarantees continuity with ringdown region.

	/* Set the second row of the basis matrix. Just the functions that multiply each unknown coefficient of the ansatz, once evaluated in the peak time */
	gsl_matrix_set(A,2,0,1.0);
	gsl_matrix_set(A,2,1,1.0); // sech(tau=0)=1
	gsl_matrix_set(A,2,2,1.0); // sech(tau=0)=1
	gsl_matrix_set(A,2,3,0.0); // tau*tau = 0 in tau=0

	/* Set A_{3,i} and b_{3}: here we impose the differentiability at the inspiral merger boundary, essentially by equating the value of the merger ansatz derivative at the
	boundary time tCut with the value of the inspiral amplitude derivative at that time. */

	/* First we compute the numerical derivatives with inspiral region for imposing differentiability at boundary.
	For this, first we compute the values of theta at two differentially close points in the boundary, from that we compute twice the orbital frequency at those points
	and then the value of the PN expansion parameter x at those points. Derivative is the amplitude difference between these two points weighted by the differential step.
	We carry the sign of the inspiral amplitude derivative. */
	REAL8 theta2 = pow(-eta*tCut/5.,-1./8);
	REAL8 theta1 = pow(-eta*(tCut-0.000001)/5,-1./8);
	REAL8 omega2 = IMRPhenomTomega22(tCut, theta2, wf, pPhase); 
	REAL8 omega1 = IMRPhenomTomega22(tCut-0.000001, theta1, wf, pPhase);
	REAL8 x1 = pow(0.5*omega1,2./3);
	REAL8 x2 = pow(0.5*omega2,2./3);
	REAL8 dampMECO = copysign(1.0,creal(IMRPhenomTInspiralAmpAnsatzHM(x2, pAmp)))*(cabs(IMRPhenomTInspiralAmpAnsatzHM(x2, pAmp)) - cabs(IMRPhenomTInspiralAmpAnsatzHM(x1, pAmp)))/0.000001; // Value of inspiral derivative at boundary time.

	gsl_vector_set(b,3,dampMECO); // We set this value to the solution vector

	/* Here we compute the needed hyperbolic  functions for the merger ansatz derivative basis matrix */
	phi = gsl_complex_rect(pAmp->alpha1RD*(tCut-pAmp->tshift),0);
	phi2 = gsl_complex_rect(2*pAmp->alpha1RD*(tCut-pAmp->tshift),0);
	sech1 = GSL_REAL(gsl_complex_sech(phi));
	sech2 = GSL_REAL(gsl_complex_sech(phi2));
	REAL8 tanh = GSL_REAL(gsl_complex_tanh(phi));
	REAL8 sinh = GSL_REAL(gsl_complex_sinh(phi2));

	/* Basis functions of the analytical time derivative of the merger ansatz */
	REAL8 aux1 = -pAmp->alpha1RD*sech1*tanh;
	REAL8 aux2 = (-2./7)*pAmp->alpha1RD*sinh*pow(sech2,8./7);
	REAL8 aux3 = 2*(tCut-pAmp->tshift);

	/*We set the value of the basis matrix with the previous elements */
	gsl_matrix_set(A,3,0,0.0);
	gsl_matrix_set(A,3,1,aux1);
	gsl_matrix_set(A,3,2,aux2);
	gsl_matrix_set(A,3,3,aux3);

	/* Once we have set up the system, we now solve the system A x = b via an LU decomposition */
	gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,x);

	/* Initialize the unkown merger ansatz coefficients with the solution of the linear system */
	pAmp->mergerC1 = gsl_vector_get(x,0);
	pAmp->mergerC2 = gsl_vector_get(x,1);
	pAmp->mergerC3 = gsl_vector_get(x,2);
	pAmp->mergerC4 = gsl_vector_get(x,3);

	/* Deallocate the gsl linear system objects */
	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_matrix_free(A);
	gsl_permutation_free(p);


	/************************************************************/
	/*** COMPLEX AMPLITUDE PHASING AT INSPIRAL MERGER BOUNDARY **/
	/************************************************************/

	/* Inspiral amplitude is a complex quantity, and then its argument contributes to the phase and frequency
	of the modes. In order to not have a phase/frequency discontinuity between inspiral and merger regions,
	we need to store the value of this phase contribution at the boundary, so it can be added
	to the merger phase/frequency ansatz later. See discussion on Sec. IID of THM paper: https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf*/


	pAmp->omegaCutPNAMP = -(ComplexAmpOrientation(x2, pAmp) - ComplexAmpOrientation(x1, pAmp))/0.000001; // Derivative of the phase contribution at boundary. Needed for frequency continuity
	pAmp->phiCutPNAMP = atan2(cimag(IMRPhenomTInspiralAmpAnsatzHM(x2, pAmp)),creal(IMRPhenomTInspiralAmpAnsatzHM(x2, pAmp))); // Phase contribution at boundary. Needed for phase continuity

	/* We need to compute an extra Pi factor for adding to the merger phase if the inspiral amplitude real part is negative */
	if(copysign(1.0,creal(IMRPhenomTInspiralAmpAnsatzHM(x2, pAmp)))==-1.0)
	{
		pAmp->phiCutPNAMP += LAL_PI;
	}

	return XLAL_SUCCESS;
}

/* Function for populating the IMRPhenomTHMPhaseStruct.
   This will solve the linear systems for obtaining lm phase ansatz coefficients from lm frequency collocation points.
   Inspiral is not present here since for (l,m)!=(2,2) the inspiral phase is rescaled from the (2,2) phase.
   For each mode, ringdown quantities are defined and stored and the merger linear system of free coefficients is solved.
   Essentially, the code is the same than for the 22 merger-ringdown, but with the right quantities for each mode. */

int IMRPhenomTSetHMPhaseCoefficients(int l, int m, IMRPhenomTHMPhaseStruct *pPhaseHM, IMRPhenomTPhase22Struct *pPhase, IMRPhenomTHMAmpStruct *pAmplm, IMRPhenomTWaveformStruct *wf)
{

	REAL8 eta = wf->eta;  // Symmetric mass ratio
	REAL8 S = wf->Shat;   // Dimensionless effective spin parameters employed for fit evaluation
	REAL8 dchi = wf->dchi; // Dimensionless spin difference chi1L - chi2L
	REAL8 delta = wf->delta; // PN Asymmetry parameter
	REAL8 tCut = tCUT_Freq;  // Inspiral ending time (t=-150M)

	/*Store twice orbital frequency at t=-150, for higher mode reconstruction*/
	REAL8 thetaCut = pow(-eta*tCut/5,-1./8); 
	REAL8 omegaCut =  IMRPhenomTomega22(tCut, thetaCut, wf, pPhase);

	REAL8 omegaCutPNAMP = pAmplm->omegaCutPNAMP; // Value of the phase derivative coming from complex inspiral amplitude at the boundary time.

	pPhaseHM->emm = m; // Magnetic number of the mode

	/* GSL objects for solving system of equations via LU decomposition */
	gsl_vector *b, *x;
	gsl_matrix *A;
	gsl_permutation *p;
	int s;

	REAL8 fRING, fDAMP, fDAMPn2; // Ringdown and damping frequencies of each mode
	REAL8 omegaCutBar, omegaMergerCP; // Rescaled frequencies at the inspiral boundary time tCut and collocation point time.
	REAL8 domegaPeak; // Rescaled frequency derivative at the ringdown boundary

	REAL8 theta2 = pow(-eta*tCut/5.,-1./8);
	REAL8 theta1 = pow(-eta*(tCut-0.0000001)/5,-1./8);
	REAL8 domegaCut = (IMRPhenomTomega22(tCut, theta2, wf, pPhase) - IMRPhenomTomega22(tCut-0.0000001, theta1, wf, pPhase))/(0.0000001); // Insoiral frequency derivative at the boundary time.

	/* Set quantities for each mode.
	Essentially, we compute the needed ringdown and damping frequencies of each mode for the ringdown ansatz:
	   -Coefficient c1 is defined in Eq. [9] of Damour&Nagar 2014 (10.1103/PhysRevD.90.024054, https://arxiv.org/abs/1406.0401).
	   -Coefficients c2 and c3 are calibrated to NR/Teukolsky data.
	   -Coefficient c4 is set to zero. 
	We also store the value of the rescaled frequencies and frequency derivatives at the boundaries, and the collocation point value. */


	if (l==2 && m==2)
	{
		/* Ringdown and damping frequencies*/
		fRING   = evaluate_QNMfit_fring22(wf->afinal) / (wf->Mfinal);
		fDAMP   = evaluate_QNMfit_fdamp22(wf->afinal) / (wf->Mfinal);
		fDAMPn2   = evaluate_QNMfit_fdamp22n2(wf->afinal) / (wf->Mfinal);

		pPhaseHM->omegaRING = 2*LAL_PI*fRING;
		pPhaseHM->alpha1RD = 2*LAL_PI*fDAMP;
		pPhaseHM->alpha2RD = 2*LAL_PI*fDAMPn2;
		pPhaseHM->alpha21RD = 0.5*(pPhaseHM->alpha2RD - pPhaseHM->alpha1RD);

		pPhaseHM->omegaRING_prec = 2*LAL_PI*evaluate_QNMfit_fring22(wf->afinal_prec) / (wf->Mfinal);
		pPhaseHM->alpha1RD_prec  = 2*LAL_PI*evaluate_QNMfit_fdamp22(wf->afinal_prec) / (wf->Mfinal);


		omegaCutBar = 1. - (omegaCut + omegaCutPNAMP)/pPhaseHM->omegaRING; // Value of the rescaled inspiral frequency at the boundary time. Notice that it includes the contribution from the complex amplitude.
		omegaMergerCP = 1. - IMRPhenomT_Merger_Freq_CP1_22(eta, S, dchi, delta)/pPhaseHM->omegaRING; // Value of the rescaled collocation point.
		pPhaseHM->omegaPeak = IMRPhenomT_PeakFrequency_22(eta, S, dchi, delta); // Value of the frequency at the merger-ringdown boundary.

		domegaCut = -domegaCut/pPhaseHM->omegaRING; //Rescaled frequency derivative at the inspiral boundary
		domegaPeak = -(IMRPhenomTRDOmegaAnsatzHM(0.0000001, pPhaseHM) - IMRPhenomTRDOmegaAnsatzHM(0., pPhaseHM))/(0.0000001)/pPhaseHM->omegaRING; //Rescaled numerical frequency derivative at the ringdown boundary
		pPhaseHM->domegaPeak = domegaPeak;

		/* Ringdown ansatz coefficients, as defined Damour&Nagar 2014 (10.1103/PhysRevD.90.024054, https://arxiv.org/abs/1406.0401), with the difference that we set c4=0 and we free c2 to acomodate a fit. */
		pPhaseHM->c3 = IMRPhenomT_RD_Freq_D3_22(eta, S, dchi, delta);
		pPhaseHM->c2 = IMRPhenomT_RD_Freq_D2_22(eta, S, dchi, delta);
		pPhaseHM->c4 = 0.0;
		pPhaseHM->c1 = (1 + pPhaseHM->c3 + pPhaseHM->c4)*(pPhaseHM->omegaRING - pPhaseHM->omegaPeak)/pPhaseHM->c2/(pPhaseHM->c3 + 2*pPhaseHM->c4);
		pPhaseHM->c1_prec = (1 + pPhaseHM->c3 + pPhaseHM->c4)*(pPhaseHM->omegaRING_prec - pPhaseHM->omegaPeak)/pPhaseHM->c2/(pPhaseHM->c3 + 2*pPhaseHM->c4);

	}

	else if (l==2 && m==1)
	{
		/* Ringdown and damping frequencies*/
		fRING   = evaluate_QNMfit_fring21(wf->afinal) / (wf->Mfinal);
		fDAMP   = evaluate_QNMfit_fdamp21(wf->afinal) / (wf->Mfinal);
		fDAMPn2   = evaluate_QNMfit_fdamp21n2(wf->afinal) / (wf->Mfinal);

		pPhaseHM->omegaRING = 2*LAL_PI*fRING;
		pPhaseHM->alpha1RD = 2*LAL_PI*fDAMP;
		pPhaseHM->alpha2RD = 2*LAL_PI*fDAMPn2;
		pPhaseHM->alpha21RD = 0.5*(pPhaseHM->alpha2RD - pPhaseHM->alpha1RD);

		pPhaseHM->omegaRING_prec = 2*LAL_PI*evaluate_QNMfit_fring21(wf->afinal_prec) / (wf->Mfinal);
		pPhaseHM->alpha1RD_prec  = 2*LAL_PI*evaluate_QNMfit_fdamp21(wf->afinal_prec) / (wf->Mfinal);

		omegaCut = 0.5*omegaCut; // Frequency at the inspiral boundary coming from the 2,2 rescaling
		omegaCutBar = 1. - (omegaCut + omegaCutPNAMP)/pPhaseHM->omegaRING; // Value of the rescaled inspiral frequency at the boundary time. Notice that it includes the contribution from the complex amplitude.
		omegaMergerCP = 1. - IMRPhenomT_Merger_Freq_CP1_21(eta, S, dchi, delta)/pPhaseHM->omegaRING; // Value of the rescaled collocation point.
		pPhaseHM->omegaPeak = IMRPhenomT_PeakFrequency_21(eta, S, dchi, delta); // Value of the frequency at the merger-ringdown boundary.

		/* Ringdown ansatz coefficients, as defined Damour&Nagar 2014 (10.1103/PhysRevD.90.024054, https://arxiv.org/abs/1406.0401), with the difference that we set c4=0 and we free c2 to acomodate a fit. */
		pPhaseHM->c3 = IMRPhenomT_RD_Freq_D3_21(eta, S, dchi, delta);
		pPhaseHM->c2 = IMRPhenomT_RD_Freq_D2_21(eta, S, dchi, delta);
		pPhaseHM->c4 = 0.0;
		pPhaseHM->c1 = (1 + pPhaseHM->c3 + pPhaseHM->c4)*(pPhaseHM->omegaRING - pPhaseHM->omegaPeak)/pPhaseHM->c2/(pPhaseHM->c3 + 2*pPhaseHM->c4);
		pPhaseHM->c1_prec = (1 + pPhaseHM->c3 + pPhaseHM->c4)*(pPhaseHM->omegaRING_prec - pPhaseHM->omegaPeak)/pPhaseHM->c2/(pPhaseHM->c3 + 2*pPhaseHM->c4);

		domegaCut = -0.5*domegaCut/pPhaseHM->omegaRING; //Rescaled frequency derivative at the inspiral boundary
		domegaPeak = -(IMRPhenomTRDOmegaAnsatzHM(0.00001, pPhaseHM) - IMRPhenomTRDOmegaAnsatzHM(0., pPhaseHM))/(0.00001)/pPhaseHM->omegaRING; //Rescaled numerical frequency derivative at the ringdown boundary
		pPhaseHM->domegaPeak = domegaPeak;

	}

	else if (l==3 && m==3)
	{
		/* Ringdown and damping frequencies*/
		fRING   = evaluate_QNMfit_fring33(wf->afinal) / (wf->Mfinal);
		fDAMP   = evaluate_QNMfit_fdamp33(wf->afinal) / (wf->Mfinal);
		fDAMPn2   = evaluate_QNMfit_fdamp33n2(wf->afinal) / (wf->Mfinal);

		pPhaseHM->omegaRING = 2*LAL_PI*fRING;
		pPhaseHM->alpha1RD = 2*LAL_PI*fDAMP;
		pPhaseHM->alpha2RD = 2*LAL_PI*fDAMPn2;
		pPhaseHM->alpha21RD = 0.5*(pPhaseHM->alpha2RD - pPhaseHM->alpha1RD);

		pPhaseHM->omegaRING_prec = 2*LAL_PI*evaluate_QNMfit_fring33(wf->afinal_prec) / (wf->Mfinal);
		pPhaseHM->alpha1RD_prec  = 2*LAL_PI*evaluate_QNMfit_fdamp33(wf->afinal_prec) / (wf->Mfinal);

		omegaCut = 1.5*omegaCut; // Frequency at the inspiral boundary coming from the 2,2 rescaling
		omegaCutBar = 1. - (omegaCut + omegaCutPNAMP)/pPhaseHM->omegaRING; // Value of the rescaled inspiral frequency at the boundary time. Notice that it includes the contribution from the complex amplitude.
		omegaMergerCP = 1. - IMRPhenomT_Merger_Freq_CP1_33(eta, S, dchi, delta)/pPhaseHM->omegaRING; // Value of the rescaled collocation point.
		pPhaseHM->omegaPeak = IMRPhenomT_PeakFrequency_33(eta, S, dchi); // Value of the frequency at the merger-ringdown boundary.

		/* Ringdown ansatz coefficients, as defined Damour&Nagar 2014 (10.1103/PhysRevD.90.024054, https://arxiv.org/abs/1406.0401), with the difference that we set c4=0 and we free c2 to acomodate a fit. */
		pPhaseHM->c3 = IMRPhenomT_RD_Freq_D3_33(eta, S, dchi, delta);
		pPhaseHM->c2 = IMRPhenomT_RD_Freq_D2_33(eta, S, dchi, delta);
		pPhaseHM->c4 = 0.0;
		pPhaseHM->c1 = (1 + pPhaseHM->c3 + pPhaseHM->c4)*(pPhaseHM->omegaRING - pPhaseHM->omegaPeak)/pPhaseHM->c2/(pPhaseHM->c3 + 2*pPhaseHM->c4);
		pPhaseHM->c1_prec = (1 + pPhaseHM->c3 + pPhaseHM->c4)*(pPhaseHM->omegaRING_prec - pPhaseHM->omegaPeak)/pPhaseHM->c2/(pPhaseHM->c3 + 2*pPhaseHM->c4);

		domegaCut = -1.5*domegaCut/pPhaseHM->omegaRING; //Rescaled frequency derivative at the inspiral boundary
		domegaPeak = -(IMRPhenomTRDOmegaAnsatzHM(0.00001, pPhaseHM) - IMRPhenomTRDOmegaAnsatzHM(0., pPhaseHM))/(0.00001)/pPhaseHM->omegaRING; //Rescaled numerical frequency derivative at the ringdown boundary
		pPhaseHM->domegaPeak = domegaPeak;

	}

	else if (l==4 && m==4)
	{
		/* Ringdown and damping frequencies*/
		fRING   = evaluate_QNMfit_fring44(wf->afinal) / (wf->Mfinal);
		fDAMP   = evaluate_QNMfit_fdamp44(wf->afinal) / (wf->Mfinal);
		fDAMPn2   = evaluate_QNMfit_fdamp44n2(wf->afinal) / (wf->Mfinal);

		pPhaseHM->omegaRING = 2*LAL_PI*fRING;
		pPhaseHM->alpha1RD = 2*LAL_PI*fDAMP;
		pPhaseHM->alpha2RD = 2*LAL_PI*fDAMPn2;
		pPhaseHM->alpha21RD = 0.5*(pPhaseHM->alpha2RD - pPhaseHM->alpha1RD);

		pPhaseHM->omegaRING_prec = 2*LAL_PI*evaluate_QNMfit_fring44(wf->afinal_prec) / (wf->Mfinal);
		pPhaseHM->alpha1RD_prec  = 2*LAL_PI*evaluate_QNMfit_fdamp44(wf->afinal_prec) / (wf->Mfinal);

		omegaCut = 2*omegaCut; // Frequency at the inspiral boundary coming from the 2,2 rescaling
		omegaCutBar = 1. - (omegaCut + omegaCutPNAMP)/pPhaseHM->omegaRING; // Value of the rescaled inspiral frequency at the boundary time. Notice that it includes the contribution from the complex amplitude.
		omegaMergerCP = 1. - IMRPhenomT_Merger_Freq_CP1_44(eta, S, dchi, delta)/pPhaseHM->omegaRING; // Value of the rescaled collocation point.
		pPhaseHM->omegaPeak = IMRPhenomT_PeakFrequency_44(eta, S, dchi, delta); // Value of the frequency at the merger-ringdown boundary.

		/* Ringdown ansatz coefficients, as defined Damour&Nagar 2014 (10.1103/PhysRevD.90.024054, https://arxiv.org/abs/1406.0401), with the difference that we set c4=0 and we free c2 to acomodate a fit. */
		pPhaseHM->c3 = IMRPhenomT_RD_Freq_D3_44(eta, S, dchi, delta);
		pPhaseHM->c2 = IMRPhenomT_RD_Freq_D2_44(eta, S, dchi, delta);
		pPhaseHM->c4 = 0.0;
		pPhaseHM->c1 = (1 + pPhaseHM->c3 + pPhaseHM->c4)*(pPhaseHM->omegaRING - pPhaseHM->omegaPeak)/pPhaseHM->c2/(pPhaseHM->c3 + 2*pPhaseHM->c4);
		pPhaseHM->c1_prec = (1 + pPhaseHM->c3 + pPhaseHM->c4)*(pPhaseHM->omegaRING_prec - pPhaseHM->omegaPeak)/pPhaseHM->c2/(pPhaseHM->c3 + 2*pPhaseHM->c4);

		domegaCut = -2.0*domegaCut/pPhaseHM->omegaRING; //Rescaled frequency derivative at the inspiral boundary
		domegaPeak = -(IMRPhenomTRDOmegaAnsatzHM(0.0000001, pPhaseHM) - IMRPhenomTRDOmegaAnsatzHM(0., pPhaseHM))/(0.0000001)/pPhaseHM->omegaRING; //Rescaled numerical frequency derivative at the ringdown boundary
		pPhaseHM->domegaPeak = domegaPeak;

	}

	else if (l==5 && m==5)
	{
		/* Ringdown and damping frequencies*/
		fRING   = evaluate_QNMfit_fring55(wf->afinal) / (wf->Mfinal);
		fDAMP   = evaluate_QNMfit_fdamp55(wf->afinal) / (wf->Mfinal);
		fDAMPn2   = evaluate_QNMfit_fdamp55n2(wf->afinal) / (wf->Mfinal);

		pPhaseHM->omegaRING = 2*LAL_PI*fRING;
		pPhaseHM->alpha1RD = 2*LAL_PI*fDAMP;
		pPhaseHM->alpha2RD = 2*LAL_PI*fDAMPn2;
		pPhaseHM->alpha21RD = 0.5*(pPhaseHM->alpha2RD - pPhaseHM->alpha1RD);

		pPhaseHM->omegaRING_prec = 2*LAL_PI*evaluate_QNMfit_fring55(wf->afinal_prec) / (wf->Mfinal);
		pPhaseHM->alpha1RD_prec  = 2*LAL_PI*evaluate_QNMfit_fdamp55(wf->afinal_prec) / (wf->Mfinal);

		omegaCut = 2.5*omegaCut; // Frequency at the inspiral boundary coming from the 2,2 rescaling
		omegaCutBar = 1. - (omegaCut + omegaCutPNAMP)/pPhaseHM->omegaRING; // Value of the rescaled inspiral frequency at the boundary time. Notice that it includes the contribution from the complex amplitude.
		omegaMergerCP = 1. - IMRPhenomT_Merger_Freq_CP1_55(eta, S, dchi, delta)/pPhaseHM->omegaRING; // Value of the rescaled collocation point.
		pPhaseHM->omegaPeak = IMRPhenomT_PeakFrequency_55(eta, S, dchi, delta); // Value of the frequency at the merger-ringdown boundary.

		/* Ringdown ansatz coefficients, as defined Damour&Nagar 2014 (10.1103/PhysRevD.90.024054, https://arxiv.org/abs/1406.0401), with the difference that we set c4=0 and we free c2 to acomodate a fit. */
		pPhaseHM->c3 = IMRPhenomT_RD_Freq_D3_55(eta, S, dchi, delta);
		pPhaseHM->c2 = IMRPhenomT_RD_Freq_D2_55(eta, S, dchi, delta);
		pPhaseHM->c4 = 0.0;
		pPhaseHM->c1 = (1 + pPhaseHM->c3 + pPhaseHM->c4)*(pPhaseHM->omegaRING - pPhaseHM->omegaPeak)/pPhaseHM->c2/(pPhaseHM->c3 + 2*pPhaseHM->c4);
		pPhaseHM->c1_prec = (1 + pPhaseHM->c3 + pPhaseHM->c4)*(pPhaseHM->omegaRING_prec - pPhaseHM->omegaPeak)/pPhaseHM->c2/(pPhaseHM->c3 + 2*pPhaseHM->c4);

		domegaCut = -2.5*domegaCut/pPhaseHM->omegaRING; //Rescaled frequency derivative at the inspiral boundary
		domegaPeak = -(IMRPhenomTRDOmegaAnsatzHM(0.0000001, pPhaseHM) - IMRPhenomTRDOmegaAnsatzHM(0., pPhaseHM))/(0.0000001)/pPhaseHM->omegaRING; //Rescaled numerical frequency derivative at the ringdown boundary
		pPhaseHM->domegaPeak = domegaPeak;

	}

	else
	{
		XLAL_ERROR(XLAL_EFUNC, "Mode not implemented in PhenomTHM. Modes available: [2,2], [2,1], [3,3], [4,4], [5,5].");
	}

	/* In order to obtain the value of the 3 free coefficients of the ansatz defined in IMRPhenomTMergerOmegaAnsatzHM, we need to solve a linear system
      where the ansatz is imposed to match a collocation point value at t=-25M and to satisfy continuity with the inspiral ansatz and differentiability with the ringdown ansatz.
      Notice that the ansatz definition IMRPhenomTMergerOmegaAnsatzHM differs from eq [27-28] of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf)
      in that here 2 of the coefficients are analytically solved for imposing continuity with the ringdown ansatz and differentiability with the inspiral ansatz.

      Reminder: the ansatz described in IMRPhenomTMergerOmegaAnsatzHM corresponds to the rescaled frequency \bar{\omega}=1 - (\omega / \omega_ring) (eq. 27 of https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf),
      so the system is solved for the rescaled frequency. For the analytical phase of IMRPhenomTMergerPhaseAnsatzHM, which is the quantity employed
      in the waveform construction, the factors are already included to produce the correct phase. */

	p = gsl_permutation_alloc(3);
	b = gsl_vector_alloc(3);
	x = gsl_vector_alloc(3);
	A = gsl_matrix_alloc(3,3);

	/* This is essentially the same procedure as in IMRPhenomTSetPhase22Coefficients */

	/* Set linear system */

	/* A_{0,i}, imposing continuity at the inspiral boundary. */

	gsl_complex phi = gsl_complex_rect(pPhaseHM->alpha1RD*tCut,0);
	REAL8 ascut = GSL_REAL(gsl_complex_arcsinh(phi));
	REAL8 ascut2 = ascut*ascut;
	REAL8 ascut3 = ascut*ascut2;
	REAL8 ascut4 = ascut*ascut3;

	gsl_vector_set(b,0,omegaCutBar - (1. - pPhaseHM->omegaPeak/pPhaseHM->omegaRING) - (pPhaseHM->domegaPeak/pPhaseHM->alpha1RD)*ascut);

	gsl_matrix_set(A,0,0,ascut2);
	gsl_matrix_set(A,0,1,ascut3);
	gsl_matrix_set(A,0,2,ascut4);

	/* A_{1,i}, imposing collocation point. */

	phi = gsl_complex_rect(pPhaseHM->alpha1RD*tcpMerger,0);
	REAL8 as025cut = GSL_REAL(gsl_complex_arcsinh(phi));
	REAL8 as025cut2 = as025cut*as025cut;
	REAL8 as025cut3 = as025cut*as025cut2;
	REAL8 as025cut4 = as025cut*as025cut3;

	gsl_vector_set(b,1,omegaMergerCP - (1. - pPhaseHM->omegaPeak/pPhaseHM->omegaRING) - (pPhaseHM->domegaPeak/pPhaseHM->alpha1RD)*as025cut);

	gsl_matrix_set(A,1,0,as025cut2);
	gsl_matrix_set(A,1,1,as025cut3);
	gsl_matrix_set(A,1,2,as025cut4);

	/* A_{2,i}, imposing differentiability at the ringdown boundary. */

	REAL8 dencut = sqrt(1.0 + tCut*tCut*pPhaseHM->alpha1RD*pPhaseHM->alpha1RD);

	gsl_matrix_set(A,2,0,2.0*pPhaseHM->alpha1RD*ascut/dencut);
	gsl_matrix_set(A,2,1,3.0*pPhaseHM->alpha1RD*ascut2/dencut);
	gsl_matrix_set(A,2,2,4.0*pPhaseHM->alpha1RD*ascut3/dencut);

	gsl_vector_set(b,2,domegaCut - pPhaseHM->domegaPeak/dencut);

	/* We now solve the system A x = b via an LU decomposition */
	gsl_linalg_LU_decomp(A,p,&s);
	gsl_linalg_LU_solve(A,p,b,x);

	/* Set merger phenomenological coefficients from solution to A x = b */

	pPhaseHM->omegaMergerC1 = gsl_vector_get(x,0);
	pPhaseHM->omegaMergerC2 = gsl_vector_get(x,1);
	pPhaseHM->omegaMergerC3 = gsl_vector_get(x,2);

	gsl_vector_free(b);
	gsl_vector_free(x);
	gsl_matrix_free(A);
	gsl_permutation_free(p);

	/* Solve phase offsets */

	pPhaseHM->phOffMerger = 0.;
	pPhaseHM->phOffRD = 0.;

	REAL8 thetabarCut = pow(-eta*tCut,-1./8);
	REAL8 phMECOinsp = (m/2.)*IMRPhenomTPhase22(tCut, thetabarCut, wf, pPhase);
	REAL8 phMECOmerger = IMRPhenomTMergerPhaseAnsatzHM(tCut, pPhaseHM);

	pPhaseHM->phOffMerger = phMECOinsp - phMECOmerger;

	pPhaseHM->phOffRD = IMRPhenomTMergerPhaseAnsatzHM(0, pPhaseHM);

	return XLAL_SUCCESS;
}

/* ---------------------------------------------------------------------------------------------------- */

/* ************************************* */
/* ******* ANSATZ DEFINITION ****** */
/* ************************************* */

/* ******* 22 ANSATZ ****** */
/* Ansatz for inspiral-merger-ringdown regions of frequency and phase of the (2,2) mode.
   Amplitude ansatz for the (2,2) is included in the general (l,m) amplitude ansatz. */

/* ******* FREQUENCY ANSATZ ****** */

/* 22 inspiral frequency ansatz is contructed with a higher order extension of TaylorT3, as described in equation 9 of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf) */
double IMRPhenomTTaylorT3(REAL8 theta, IMRPhenomTPhase22Struct *pPhase){

   REAL8 theta2 = theta*theta;
   REAL8 theta3 = theta2*theta;
   REAL8 theta4 = theta2*theta2;
   REAL8 theta5 = theta3*theta2;
   REAL8 theta6 = theta3*theta3;
   REAL8 theta7 = theta4*theta3;

   REAL8 fac = theta3/8; 
   REAL8 logterm = (107*log(theta))/280.;

   REAL8 out = (1 + pPhase->omega1PN*theta2 + pPhase->omega1halfPN*theta3 + pPhase->omega2PN*theta4 + pPhase->omega2halfPN*theta5 + pPhase->omega3PN*theta6 + logterm*theta6 + pPhase->omega3halfPN*theta7);

   return 2*fac*out;
}

double IMRPhenomTInspiralOmegaAnsatz22(REAL8 theta, IMRPhenomTPhase22Struct *pPhase){

   REAL8 theta8 = pow(theta,8);
   REAL8 theta9 = theta8*theta;
   REAL8 theta10 = theta9*theta;
   REAL8 theta11 = theta10*theta;
   REAL8 theta12 = theta11*theta;
   REAL8 theta13 = theta12*theta;

   REAL8 fac = theta*theta*theta/8; 

   REAL8 taylort3 = IMRPhenomTTaylorT3(theta, pPhase);

   REAL8 out =  pPhase->omegaInspC1*theta8 + pPhase->omegaInspC2*theta9 + pPhase->omegaInspC3*theta10 + pPhase->omegaInspC4*theta11 + pPhase->omegaInspC5*theta12 + pPhase->omegaInspC6*theta13;

   return taylort3 + 2*fac*out;
}

/* 22 merger frequency ansatz is a phenomenological expression described in equation [27-28] of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf)
   The ansatz described here corresponds to the rescaled frequency \bar{\omega}=1 - (\omega / \omega_ring) of equation 27. */
double IMRPhenomTMergerOmegaAnsatz22(REAL8 t, IMRPhenomTPhase22Struct *pPhase){

	gsl_complex phi = gsl_complex_rect(pPhase->alpha1RD*t,0);
   REAL8 x = GSL_REAL(gsl_complex_arcsinh(phi));
   REAL8 x2 = x*x;
   REAL8 x3 = x2*x;
   REAL8 x4 = x2*x2;

   REAL8 out = 1. - pPhase->omegaPeak/pPhase->omegaRING + (pPhase->domegaPeak/pPhase->alpha1RD)*x + pPhase->omegaMergerC1*x2 + pPhase->omegaMergerC2*x3 + pPhase->omegaMergerC3*x4;

   return out;
}

/* 22 ringdown frequency ansatz is the analytical derivate of the expression in IMRPhenomTRDPhaseAnsatz22, described in equation 24 of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf) */
double IMRPhenomTRDOmegaAnsatz22(REAL8 t, IMRPhenomTPhase22Struct *pPhase){

   REAL8 c3 = pPhase->c3;
   REAL8 c4 = pPhase->c4;
   REAL8 c2 = pPhase->c2;
   REAL8 c1 = pPhase->c1;


   REAL8 expC = exp(-1*c2*t);
   REAL8 expC2 = expC*expC;



   REAL8 num = c1*(-2*c2*c4*expC2 - c2*c3*expC);
   REAL8 den = 1 + c4*expC2 + c3*expC;

   return (num/den + pPhase->omegaRING);
}

/* ******** PHASE ANSATZ ********* */

double IMRPhenomTInspiralPhaseTaylorT3(REAL8 thetabar, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase){

	REAL8 eta = wf->eta;

	REAL8 aux = (pow(5,-0.625)*pow(eta,-1)*pow(thetabar,-5)*(-168 - 280*pPhase->omega1PN*pow(5,0.25)*pow(thetabar,2) - 420*pPhase->omega1halfPN*pow(5,0.375)*pow(thetabar,3) - 
       840*pPhase->omega2PN*pow(5,0.5)*pow(thetabar,4) + 840*pPhase->omega2halfPN*log(thetabar)*pow(5,0.625)*pow(thetabar,5) - 321*pow(5,0.75)*pow(thetabar,6) + 
       840*pPhase->omega3PN*pow(5,0.75)*pow(thetabar,6) + 321*log(thetabar*pow(5,0.125))*pow(5,0.75)*pow(thetabar,6) + 420*pPhase->omega3halfPN*pow(5,0.875)*pow(thetabar,7)))/84.;

    return aux;
}

/* 22 inspiral phase is obtained from analytical integral of the frequency ansatz, which containts TaylorT3 + the extra calibrated terms.
   Then, it depends on the omega{N}PN coefficients of TaylorT3 and on the solution to the extra coefficients omegaInspCi. */
double IMRPhenomTInspiralPhaseAnsatz22(REAL8 t, REAL8 thetabar, IMRPhenomTWaveformStruct *wf, IMRPhenomTPhase22Struct *pPhase){

	REAL8 eta = wf->eta;

	REAL8 aux = -(pow(5,-0.625)*pow(eta,-2)*pow(t,-1)*pow(thetabar,-7)*(3*(-107 + 280*pPhase->omega3PN)*pow(5,0.75) + 321*log(thetabar*pow(5,0.125))*pow(5,0.75) + 420*pPhase->omega3halfPN*thetabar*pow(5,0.875) + 
        56*(25*pPhase->omegaInspC1 + 3*eta*t)*pow(thetabar,2) + 1050*pPhase->omegaInspC2*pow(5,0.125)*pow(thetabar,3) + 280*(3*pPhase->omegaInspC3 + eta*pPhase->omega1PN*t)*pow(5,0.25)*pow(thetabar,4) + 
        140*(5*pPhase->omegaInspC4 + 3*eta*pPhase->omega1halfPN*t)*pow(5,0.375)*pow(thetabar,5) + 120*(5*pPhase->omegaInspC5 + 7*eta*pPhase->omega2PN*t)*pow(5,0.5)*pow(thetabar,6) + 
        525*pPhase->omegaInspC6*pow(5,0.625)*pow(thetabar,7) + 105*eta*pPhase->omega2halfPN*t*log(-t)*pow(5,0.625)*pow(thetabar,7)))/84.;

    return aux + pPhase->phOffInsp;
}

/* 22 merger phase ansatz is an analytical integral of the 22 frequency ansatz of IMRPhenomTMergerOmegaAnsatz22 */
double IMRPhenomTMergerPhaseAnsatz22(REAL8 t, IMRPhenomTPhase22Struct *pPhase){

	gsl_complex phi = gsl_complex_rect(pPhase->alpha1RD*t,0);
   REAL8 x = GSL_REAL(gsl_complex_arcsinh(phi));

   REAL8 alpha1RD = pPhase->alpha1RD;
   REAL8 omegaPeak = pPhase->omegaPeak;
   REAL8 domegaPeak = pPhase->domegaPeak;
   REAL8 omegaRING = pPhase->omegaRING;

   REAL8 cc = pPhase->omegaMergerC1;
   REAL8 dd = pPhase->omegaMergerC2;
   REAL8 ee = pPhase->omegaMergerC3;

   REAL8 aux = omegaRING*t - omegaRING*(2*cc*t + 24*ee*t + 6*dd*t*x + domegaPeak*t*x*pow(alpha1RD,-1) + t*(1 - omegaPeak*pow(omegaRING,-1)) + cc*t*pow(x,2) + 12*ee*t*pow(x,2) + dd*t*pow(x,3) +
      ee*t*pow(x,4) - domegaPeak*pow(alpha1RD,-2)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5) - 6*dd*pow(alpha1RD,-1)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5) -
      2*cc*x*pow(alpha1RD,-1)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5) - 24*ee*x*pow(alpha1RD,-1)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5) -
      3*dd*pow(alpha1RD,-1)*pow(x,2)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5) - 4*ee*pow(alpha1RD,-1)*pow(x,3)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5));

   return (aux + pPhase->phOffMerger);
}

/* 22 ringdown phase ansatz is taken from equation 5 in Damour&Nagar 2014 (https://arxiv.org/pdf/1406.0401.pdf).
   By construction is the analytical integral of the ansatz in IMRPhenomTRDOmegaAnsatz22. */
double IMRPhenomTRDPhaseAnsatz22(REAL8 t, IMRPhenomTPhase22Struct *pPhase){

   REAL8 c3 = pPhase->c3;
   REAL8 c4 = pPhase->c4;
   REAL8 c2 = pPhase->c2;
   REAL8 c1 = pPhase->c1_prec;


   REAL8 expC = exp(-c2*t);
   REAL8 expC2 = expC*expC;

   REAL8 num = 1 + c3*expC + c4*expC2;
   REAL8 den = 1 + c3 + c4;
   REAL8 aux = log(num/den);

   return (c1*aux + pPhase->omegaRING_prec*t + pPhase->phOffRD);
}

/* ---------------------------------------------------------------------------------------------------- */

/* ******************* HIGHER MODES ANSATZ ******************* */

/* Ansatz for inspiral-merger-ringdown regions of amplitude of (l,m) modes,
   and ansatz for merger-ringdown regions of frequency and phase of (l,m) modes.
   (inspiral region for frequency and phase is rescaled from the 22 result). */


/* ***** FREQUENCY ANSATZ ********* */
/* Same ansatz presented in IMRPhenomTMergerOmegaAnsatz22 and IMRPhenomTRDOmegaAnsatz22,
   but with ringdown quantities and coefficient fits for the particular (l,m) mode,
   passed through the IMRPhenomTHMPhaseStruct */

double IMRPhenomTMergerOmegaAnsatzHM(REAL8 t, IMRPhenomTHMPhaseStruct *pPhase){

	gsl_complex phi = gsl_complex_rect(pPhase->alpha1RD*t,0);
   REAL8 x = GSL_REAL(gsl_complex_arcsinh(phi));
   REAL8 x2 = x*x;
   REAL8 x3 = x2*x;
   REAL8 x4 = x2*x2;

   REAL8 out = 1. - pPhase->omegaPeak/pPhase->omegaRING + (pPhase->domegaPeak/pPhase->alpha1RD)*x + pPhase->omegaMergerC1*x2 + pPhase->omegaMergerC2*x3 + pPhase->omegaMergerC3*x4;

   return out;
}

double IMRPhenomTRDOmegaAnsatzHM(REAL8 t, IMRPhenomTHMPhaseStruct *pPhase){

   REAL8 c3 = pPhase->c3;
   REAL8 c4 = pPhase->c4;
   REAL8 c2 = pPhase->c2;
   REAL8 c1 = pPhase->c1;


   REAL8 expC = exp(-c2*t);
   REAL8 expC2 = expC*expC;



   REAL8 num = c1*(-2*c2*c4*expC2 - c2*c3*expC);
   REAL8 den = 1 + c4*expC2 + c3*expC;

   return (num/den + pPhase->omegaRING);
}

/* ***** PHASE ANSATZ ********* */
/* Same ansatz presented in IMRPhenomTMergerPhaseAnsatz22 and IMRPhenomTRDPhaseAnsatz22,
   but with ringdown quantities and coefficient fits for the particular (l,m) mode,
   passed through the IMRPhenomTHMPhaseStruct */

double IMRPhenomTMergerPhaseAnsatzHM(REAL8 t, IMRPhenomTHMPhaseStruct *pPhase){

	gsl_complex phi = gsl_complex_rect(pPhase->alpha1RD*t,0);
   REAL8 x = GSL_REAL(gsl_complex_arcsinh(phi));

   REAL8 alpha1RD = pPhase->alpha1RD;
   REAL8 omegaPeak = pPhase->omegaPeak;
   REAL8 domegaPeak = pPhase->domegaPeak;
   REAL8 omegaRING = pPhase->omegaRING;

   REAL8 cc = pPhase->omegaMergerC1;
   REAL8 dd = pPhase->omegaMergerC2;
   REAL8 ee = pPhase->omegaMergerC3;

   REAL8 aux = omegaRING*t - omegaRING*(2*cc*t + 24*ee*t + 6*dd*t*x + domegaPeak*t*x*pow(alpha1RD,-1) + t*(1 - omegaPeak*pow(omegaRING,-1)) + cc*t*pow(x,2) + 12*ee*t*pow(x,2) + dd*t*pow(x,3) +
      ee*t*pow(x,4) - domegaPeak*pow(alpha1RD,-2)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5) - 6*dd*pow(alpha1RD,-1)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5) -
      2*cc*x*pow(alpha1RD,-1)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5) - 24*ee*x*pow(alpha1RD,-1)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5) -
      3*dd*pow(alpha1RD,-1)*pow(x,2)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5) - 4*ee*pow(alpha1RD,-1)*pow(x,3)*pow(1 + pow(alpha1RD,2)*pow(t,2),0.5));

   return (aux + pPhase->phOffMerger);
}

double IMRPhenomTRDPhaseAnsatzHM(REAL8 t, IMRPhenomTHMPhaseStruct *pPhase){

   REAL8 c3 = pPhase->c3;
   REAL8 c4 = pPhase->c4;
   REAL8 c2 = pPhase->c2;
   REAL8 c1 = pPhase->c1_prec;


   REAL8 expC = exp(-1*c2*t);
   REAL8 expC2 = expC*expC;

   REAL8 num = 1 + c3*expC + c4*expC2;
   REAL8 den = 1 + c3 + c4;
   REAL8 aux = log(num/den);

   return (c1*aux + pPhase->omegaRING_prec*t + pPhase->phOffRD);
}

/* *********** AMPLITUDE ANSATZ ************** */

/* Ansatz for the inspiral amplitude of a general (l,m) mode. PN amplitude coefficients for each mode are defined in IMRPhenomTSetHMAmplitudeCoefficients
   and are passed through the IMRPhenomTHMAmpStruct defined for the partucular mode. Expression is extended to higher order through the addition of three extra terms
   whose coefficients inspCi are solved in IMRPhenomTSetHMAmplitudeCoefficients.
   General form of the expression is described in eq [14] of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf).*/
COMPLEX16 IMRPhenomTInspiralAmpAnsatzHM(REAL8 x, IMRPhenomTHMAmpStruct *pAmp)
{
	REAL8 fac = pAmp->fac0*x;

	REAL8 xhalf = sqrt(x);
	REAL8 x1half = x*xhalf;
	REAL8 x2 = x*x;
	REAL8 x2half = x2*xhalf;
	REAL8 x3 = x2*x;
	REAL8 x3half = x3*xhalf;
	REAL8 x4 = x2*x2;
	REAL8 x4half = x4*xhalf;
	REAL8 x5 = x3*x2;

	REAL8 ampreal = pAmp->ampN + pAmp->amp0halfPNreal*xhalf + pAmp->amp1PNreal*x + pAmp->amp1halfPNreal*x1half + pAmp->amp2PNreal*x2  + pAmp->amp2halfPNreal*x2half + pAmp->amp3PNreal*x3 + pAmp->amp3halfPNreal*x3half + pAmp->amplog*log(16*x)*x3;
	REAL8 ampimag = pAmp->amp0halfPNimag*xhalf + pAmp->amp1PNimag*x + pAmp->amp1halfPNimag*x1half + pAmp->amp2PNimag*x2  + pAmp->amp2halfPNimag*x2half + pAmp->amp3PNimag*x3 + pAmp->amp3halfPNimag*x3half;
	COMPLEX16 amp = crect(ampreal + pAmp->inspC1*x4 + pAmp->inspC2*x4half + pAmp->inspC3*x5, ampimag);
	
	return fac*amp;
}

/* Wrapper to compute directly the phase contribution of the complex inspiral amplitude at a specific value of the PN expansion parameter x.
See equation 18 of PhenomTHM paper: https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf */
double ComplexAmpOrientation(REAL8 xref, IMRPhenomTHMAmpStruct *pAmp)
{
	REAL8 xhalf = sqrt(xref);
	REAL8 x1half = xref*xhalf;
	REAL8 x2 = xref*xref;
	REAL8 x2half = x2*xhalf;
	REAL8 x3 = x2*xref;
	REAL8 x3half = x3*xhalf;
	REAL8 x4 = x2*x2;
	REAL8 x4half = x4*xhalf;
	REAL8 x5 = x3*x2;

	REAL8 ampreal = pAmp->ampN + pAmp->amp0halfPNreal*xhalf + pAmp->amp1PNreal*xref + pAmp->amp1halfPNreal*x1half + pAmp->amp2PNreal*x2  + pAmp->amp2halfPNreal*x2half + pAmp->amp3PNreal*x3 + pAmp->amp3halfPNreal*x3half + pAmp->amplog*log(16*xref)*x3;
	REAL8 ampimag = pAmp->amp0halfPNimag*xhalf + pAmp->amp1PNimag*xref + pAmp->amp1halfPNimag*x1half + pAmp->amp2PNimag*x2  + pAmp->amp2halfPNimag*x2half + pAmp->amp3PNimag*x3 + pAmp->amp3halfPNimag*x3half;
	
	return atan2(ampimag,ampreal + pAmp->inspC1*x4 + pAmp->inspC2*x4half + pAmp->inspC3*x5);
}

/* Ansatz for the merger amplitude of a general (l,m) mode. A phenomenological expression in terms of hyperbolic functions is employed,
   as described in eq 29 of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf) */
double IMRPhenomTMergerAmpAnsatzHM(REAL8 t, IMRPhenomTHMAmpStruct *pAmp)
{
	double tpeak = pAmp->tshift;

	gsl_complex phi = gsl_complex_rect(pAmp->alpha1RD*(t-tpeak),0);
	gsl_complex phi2 = gsl_complex_rect(2*pAmp->alpha1RD*(t-tpeak),0);

	REAL8 sech1 = GSL_REAL(gsl_complex_sech(phi));
	REAL8 sech2 = GSL_REAL(gsl_complex_sech(phi2)); 

	REAL8 aux = pAmp->mergerC1 + pAmp->mergerC2*sech1 + pAmp->mergerC3*pow(sech2,1./7) + pAmp->mergerC4*(t-tpeak)*(t-tpeak);
	return aux;
}

/* Ansatz for the ringdown amplitude of a general (l,m) mode. Equivalent to eq 4 of Damour&Nagar 2014 (https://arxiv.org/pdf/1406.0401.pdf).
   Also described in eq 25 of PhenomTHM paper (https://dcc.ligo.org/DocDB/0172/P2000524/001/PhenomTHM_SH-3.pdf). */
double IMRPhenomTRDAmpAnsatzHM(REAL8 t, IMRPhenomTHMAmpStruct *pAmp)
{
	double tpeak = pAmp->tshift;
	REAL8 c3 = pAmp->c3;
  REAL8 c4 = pAmp->c4_prec;
  REAL8 c2 = pAmp->c2_prec;
  REAL8 c1 = pAmp->c1_prec;

  gsl_complex phi = gsl_complex_rect(c2*(t-tpeak) + c3,0);
  REAL8 tanh = GSL_REAL(gsl_complex_tanh(phi));

  REAL8 expAlpha = exp(-pAmp->alpha1RD_prec*(t-tpeak));

  return expAlpha*(c1*tanh + c4);
}

/* ---------------------------------------------------------------------------------------------------- */

/* ****************************************** */
/* ******* PIECEWISE EXPRESSIONS ************ */
/* ****************************************** */

/* ***** Functions for obtaining 22 phase and frequency for the whole time series **** */
/* These are needed as separate functions from the general (l,m) piecewise functions because
   22 phase and frequency need to be precomputed before each mode structure is invoked. */

double IMRPhenomTomega22(
  REAL8 t,
  REAL8 theta,
  IMRPhenomTWaveformStruct *pWF,
  IMRPhenomTPhase22Struct *pPhase
)
{
  REAL8 w;

  	if(t < pPhase->tEarly && pWF->inspVersion!=0) // If non-default recounstruction (4 regions) computes early inspiral region with pure TaylorT3
  	{
  		w = IMRPhenomTTaylorT3(theta, pPhase);
  	}
  	else if(t < pPhase->tCut22 - pPhase->dtM) // For times earlier than the inspiral-merger boundary, computes inspiral frequency
    {
      	w = IMRPhenomTInspiralOmegaAnsatz22(theta, pPhase);
    }
    else if(t > 0) // For times later than the 22 peak amplitude time, computes ringdown frequencies
    {
        w = IMRPhenomTRDOmegaAnsatz22(t, pPhase);
    }
    else // Remaining thing is to compute merger frequency
    {
        w = IMRPhenomTMergerOmegaAnsatz22(t, pPhase);
        w = pPhase->omegaRING*(1. - w); // This is done for correcting the rescaling in the merger, as explained in IMRPhenomTMergerOmegaAnsatz22
    }

  return w;
}

double IMRPhenomTPhase22(
  REAL8 t,
  REAL8 thetabar,
  IMRPhenomTWaveformStruct *pWF,
  IMRPhenomTPhase22Struct *pPhase
)
{
  REAL8 ph;

  if(t < pPhase->tEarly && pWF->inspVersion!=0) // If non-default recounstruction (4 regions) computes early inspiral region with pure TaylorT3
  	{
  		ph = IMRPhenomTInspiralPhaseTaylorT3(thetabar, pWF, pPhase);
  	}
  else if(t < pPhase->tCut22 - pPhase->dtM) // For times earlier than the inspiral-merger boundary, computes inspiral phase
        {
          ph = IMRPhenomTInspiralPhaseAnsatz22(t, thetabar, pWF, pPhase);
        }
        else if(t > 0)
        {
          ph = IMRPhenomTRDPhaseAnsatz22(t, pPhase); // For times later than the 22 peak amplitude time, computes ringdown phase
        }
        else
        {
          ph = IMRPhenomTMergerPhaseAnsatz22(t, pPhase); // Remaining thing is to compute merger phase
        }

  return ph;
}

/* ***** Functions for obtaining (l,m) phase and amplitude for the whole time series.
Explanation of the piece-wise function is the same as above, but now boundary times are different and there are only 3 regions for each mode. **** */

double IMRPhenomTHMPhase(
  	REAL8 t,
  	REAL8 phiInsp,
  	IMRPhenomTHMPhaseStruct *pPhaseHM,
  	UNUSED IMRPhenomTHMAmpStruct *pAmpHM
)	
{
  REAL8 ph;

  if(t < tCUT_Freq)
        {
        	ph = (pPhaseHM->emm/2.)*phiInsp;
        }
        else if(t > 0)
        {
        	ph = IMRPhenomTRDPhaseAnsatzHM(t, pPhaseHM) - pAmpHM->phiCutPNAMP;
        }
        else
        {
        	ph = IMRPhenomTMergerPhaseAnsatzHM(t, pPhaseHM) - pAmpHM->phiCutPNAMP;
        }

  return ph;
}

COMPLEX16 IMRPhenomTHMAmp(
  REAL8 t,
  UNUSED REAL8 x,
  IMRPhenomTHMAmpStruct *pAmp
)
{
    COMPLEX16 amp;

    if(t < tCUT_Amp)
        {
          amp = IMRPhenomTInspiralAmpAnsatzHM(x, pAmp);
        }
        else if(t > pAmp->tshift)
        {
          amp = IMRPhenomTRDAmpAnsatzHM(t, pAmp);
        }
        else
        {
          amp = IMRPhenomTMergerAmpAnsatzHM(t, pAmp);
        }

    return amp;
}

/* ---------------------------------------------------------------------------------------------------- */

/******** Wrapper function for GSL root finder *****/

double GetTimeOfFreq(double t, void *params)
{
	struct FindRootParams *p;

	p = (struct FindRootParams *)params;

	REAL8 theta;

	if(t < p->pPhase->tEarly)
	{
		theta = pow(p->wf->eta*(p->pPhase->tt0-t)/5,-1./8);
	}
	else
	{
		theta = pow(p->wf->eta*(-t)/5,-1./8);
	}

	return(LAL_TWOPI*p->f0 - IMRPhenomTomega22(t, theta, p->wf, p->pPhase));
}

REAL8 GetEulerSlope(REAL8 af, REAL8 mf)
{
	REAL8 EulerRDslope = 2*LAL_PI*(evaluate_QNMfit_fring22(af) / (mf) - evaluate_QNMfit_fring21(af) / (mf)); // FIXME: 
	if(af<0)
	{
		EulerRDslope = -EulerRDslope;
	}

	return EulerRDslope;
}
