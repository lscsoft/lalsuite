/**
 * \author Yi Pan, Craig Robinson, Prayush Kumar, Stas Babak, Andrea Taracchini

 *
 * \brief Module to compute the ring-down waveform as linear combination
 * of quasi-normal-modes decaying waveforms, which can be attached to
 * the inspiral part of the compat binary coalescing waveform.
 * The method is describe in Sec. II C of Pan et al. PRD 84, 124052 (2011),
 * specifically Eqs. 30 - 32.
 * Eqs. 30 and 31 are written in explicity linear equation systems in
 * DCC document T1100433.
 * This method is currently used for EOBNRv2 and SEOBNRv1 models. The only difference
 * between the two models in ring-down waveform is the pseudo-QNM introduced
 * in the latter (see Taracchini et al. PRD 86, 024011 (2012) for more details).
 */
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_blas.h>

//#include "LALSimIMREOBNRv2.h"
//#include "LALSimBlackHoleRingdown.h"
#include "LALSimBlackHoleRingdownPrec.h"
#include "LALSimIMREOBNQCCorrection.c"
//#include "LALSimIMREOBHybridRingdown.c"

#ifndef _LALSIMIMREOBHYBRIDRINGDOWNPREC_C
#define _LALSIMIMREOBHYBRIDRINGDOWNPREC_C

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Generates the ringdown wave associated with the given real
 * and imaginary parts of the inspiral waveform. The parameters of
 * the ringdown, such as amplitude and phase offsets, are determined
 * by solving the linear equations defined in the DCC document T1100433.
 * In the linear equations Ax=y,
 * A is a 16-by-16 matrix depending on QNM (complex) frequencies,
 * x is a 16-d vector of the 8 unknown complex QNM amplitudes,
 * y is a 16-d vector depending on inspiral-plunge waveforms and their derivatives near merger.
 */
static INT4 XLALSimIMREOBHybridRingdownWave(
                                            REAL8Vector          *rdwave1,   /**<< OUTPUT, Real part of ringdown waveform */
                                            REAL8Vector          *rdwave2,   /**<< OUTPUT, Imag part of ringdown waveform */
                                            const REAL8           dt,        /**<< Sampling interval */
                                            const REAL8           mass1,     /**<< First component mass (in Solar masses) */
                                            const REAL8           mass2,     /**<< Second component mass (in Solar masses) */
                                            REAL8VectorSequence  *inspwave1, /**<< Values and derivs of real part inspiral waveform */
                                            REAL8VectorSequence  *inspwave2, /**<< Values and derivs of imag part inspiral waveform */
                                            COMPLEX16Vector      *modefreqs, /**<< Complex freqs of ringdown (scaled by total mass) */
                                            REAL8Vector          *matchrange /**<< Times which determine the comb of ringdown attachment */
)
{

    /* XLAL error handling */
    INT4 errcode = XLAL_SUCCESS;

    /* For checking GSL return codes */
    INT4 gslStatus;
    INT4 debugSB = 0;

    UINT4 i, j, k, nmodes = 8;

    /* Sampling rate from input */
    REAL8 t1, t2, t3, t4, t5, rt;
    gsl_matrix *coef;
    gsl_vector *hderivs;
    gsl_vector *x;
    gsl_permutation *p;
    REAL8Vector *modeamps;
    int s;
    REAL8 tj;
    REAL8 m;

    /* mass in geometric units */
    m  = (mass1 + mass2) * LAL_MTSUN_SI;
    t5 = (matchrange->data[0] - matchrange->data[1]) * m;
    rt = -t5 / 5.;

    t4 = t5 + rt;
    t3 = t4 + rt;
    t2 = t3 + rt;
    t1 = t2 + rt;

    if ( inspwave1->length != 3 || inspwave2->length != 3 ||
        modefreqs->length != nmodes )
    {
        XLAL_ERROR( XLAL_EBADLEN );
    }

    /* Solving the linear system for QNMs amplitude coefficients using gsl routine */
    /* Initiate matrices and supporting variables */
    XLAL_CALLGSL( coef = (gsl_matrix *) gsl_matrix_alloc(2 * nmodes, 2 * nmodes) );
    XLAL_CALLGSL( hderivs = (gsl_vector *) gsl_vector_alloc(2 * nmodes) );
    XLAL_CALLGSL( x = (gsl_vector *) gsl_vector_alloc(2 * nmodes) );
    XLAL_CALLGSL( p = (gsl_permutation *) gsl_permutation_alloc(2 * nmodes) );

    /* Check all matrices and variables were allocated */
    if ( !coef || !hderivs || !x || !p )
    {
        if (coef)    gsl_matrix_free(coef);
        if (hderivs) gsl_vector_free(hderivs);
        if (x)       gsl_vector_free(x);
        if (p)       gsl_permutation_free(p);

        XLAL_ERROR( XLAL_ENOMEM );
    }

    /* Define the linear system Ax=y */
    /* Matrix A (2*n by 2*n) has block symmetry. Define half of A here as "coef" */
    /* The half of A defined here corresponds to matrices M1 and -M2 in the DCC document T1100433 */
    /* Define y here as "hderivs" */
    for (i = 0; i < nmodes; ++i)
    {
        gsl_matrix_set(coef, 0, i, 1);
        gsl_matrix_set(coef, 1, i, - cimag(modefreqs->data[i]));
        gsl_matrix_set(coef, 2, i, exp(-cimag(modefreqs->data[i])*t1) * cos(creal(modefreqs->data[i])*t1));
        gsl_matrix_set(coef, 3, i, exp(-cimag(modefreqs->data[i])*t2) * cos(creal(modefreqs->data[i])*t2));
        gsl_matrix_set(coef, 4, i, exp(-cimag(modefreqs->data[i])*t3) * cos(creal(modefreqs->data[i])*t3));
        gsl_matrix_set(coef, 5, i, exp(-cimag(modefreqs->data[i])*t4) * cos(creal(modefreqs->data[i])*t4));
        gsl_matrix_set(coef, 6, i, exp(-cimag(modefreqs->data[i])*t5) * cos(creal(modefreqs->data[i])*t5));
        gsl_matrix_set(coef, 7, i, exp(-cimag(modefreqs->data[i])*t5) *
                       (-cimag(modefreqs->data[i]) * cos(creal(modefreqs->data[i])*t5)
                        -creal(modefreqs->data[i]) * sin(creal(modefreqs->data[i])*t5)));
        gsl_matrix_set(coef, 8, i, 0);
        gsl_matrix_set(coef, 9, i, - creal(modefreqs->data[i]));
        gsl_matrix_set(coef, 10, i, -exp(-cimag(modefreqs->data[i])*t1) * sin(creal(modefreqs->data[i])*t1));
        gsl_matrix_set(coef, 11, i, -exp(-cimag(modefreqs->data[i])*t2) * sin(creal(modefreqs->data[i])*t2));
        gsl_matrix_set(coef, 12, i, -exp(-cimag(modefreqs->data[i])*t3) * sin(creal(modefreqs->data[i])*t3));
        gsl_matrix_set(coef, 13, i, -exp(-cimag(modefreqs->data[i])*t4) * sin(creal(modefreqs->data[i])*t4));
        gsl_matrix_set(coef, 14, i, -exp(-cimag(modefreqs->data[i])*t5) * sin(creal(modefreqs->data[i])*t5));
        gsl_matrix_set(coef, 15, i, exp(-cimag(modefreqs->data[i])*t5) *
                       ( cimag(modefreqs->data[i]) * sin(creal(modefreqs->data[i])*t5)
                        -creal(modefreqs->data[i]) * cos(creal(modefreqs->data[i])*t5)));
    }
    for (i = 0; i < 2; ++i)
    {
        gsl_vector_set(hderivs, i, inspwave1->data[(i + 1) * inspwave1->vectorLength - 1]);
        gsl_vector_set(hderivs, i + nmodes, inspwave2->data[(i + 1) * inspwave2->vectorLength - 1]);
        gsl_vector_set(hderivs, i + 6, inspwave1->data[i * inspwave1->vectorLength]);
        gsl_vector_set(hderivs, i + 6 + nmodes, inspwave2->data[i * inspwave2->vectorLength]);
    }
    gsl_vector_set(hderivs, 2, inspwave1->data[4]);
    gsl_vector_set(hderivs, 2 + nmodes, inspwave2->data[4]);
    gsl_vector_set(hderivs, 3, inspwave1->data[3]);
    gsl_vector_set(hderivs, 3 + nmodes, inspwave2->data[3]);
    gsl_vector_set(hderivs, 4, inspwave1->data[2]);
    gsl_vector_set(hderivs, 4 + nmodes, inspwave2->data[2]);
    gsl_vector_set(hderivs, 5, inspwave1->data[1]);
    gsl_vector_set(hderivs, 5 + nmodes, inspwave2->data[1]);

    /* Complete the definition for the rest half of A */
    for (i = 0; i < nmodes; ++i)
    {
        for (k = 0; k < nmodes; ++k)
        {
            gsl_matrix_set(coef, i, k + nmodes, - gsl_matrix_get(coef, i + nmodes, k));
            gsl_matrix_set(coef, i + nmodes, k + nmodes, gsl_matrix_get(coef, i, k));
        }
    }

#if 0
    /* print ringdown-matching linear system: coefficient matrix and RHS vector */
    XLAL_PRINT_INFO("\nRingdown matching matrix:\n");
    for (i = 0; i < 16; ++i)
    {
        for (j = 0; j < 16; ++j)
        {
            XLAL_PRINT_INFO("%.12e ",gsl_matrix_get(coef,i,j));
        }
        XLAL_PRINT_INFO("\n");
    }
    XLAL_PRINT_INFO("RHS:  ");
    for (i = 0; i < 16; ++i)
    {
        XLAL_PRINT_INFO("%.12e   ",gsl_vector_get(hderivs,i));
    }
    XLAL_PRINT_INFO("\n");
#endif

    /* Call gsl LU decomposition to solve the linear system */
    XLAL_CALLGSL( gslStatus = gsl_linalg_LU_decomp(coef, p, &s) );
    if ( gslStatus == GSL_SUCCESS )
    {
        XLAL_CALLGSL( gslStatus = gsl_linalg_LU_solve(coef, p, hderivs, x) );
    }
    if ( gslStatus != GSL_SUCCESS )
    {
        gsl_matrix_free(coef);
        gsl_vector_free(hderivs);
        gsl_vector_free(x);
        gsl_permutation_free(p);
        XLAL_ERROR( XLAL_EFUNC );
    }

    /* Putting solution to an XLAL vector */
    modeamps = XLALCreateREAL8Vector(2 * nmodes);

    if ( !modeamps )
    {
        gsl_matrix_free(coef);
        gsl_vector_free(hderivs);
        gsl_vector_free(x);
        gsl_permutation_free(p);
        XLAL_ERROR( XLAL_ENOMEM );
    }

    for (i = 0; i < nmodes; ++i)
    {
        modeamps->data[i] = gsl_vector_get(x, i);
        modeamps->data[i + nmodes] = gsl_vector_get(x, i + nmodes);
    }

    /* Free all gsl linear algebra objects */
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);

    /* Build ring-down waveforms */

    REAL8 timeOffset = fmod( matchrange->data[1], dt/m) * dt;

    if (debugSB){
       /// Print the solution: 
       for (i=0; i<nmodes; i++){
           printf("RD info: QNM: (re) %.16e, (im) %.16e, Amp %16e \n", creal(modefreqs->data[i]), 
                    cimag(modefreqs->data[i]), modeamps->data[i]);
       }

    }

    for (j = 0; j < rdwave1->length; ++j)
    {
        tj = j * dt - timeOffset;
        rdwave1->data[j] = 0;
        rdwave2->data[j] = 0;
        for (i = 0; i < nmodes; ++i)
        {
            rdwave1->data[j] += exp(- tj * cimag(modefreqs->data[i]))
            * ( modeamps->data[i] * cos(tj * creal(modefreqs->data[i]))
               +   modeamps->data[i + nmodes] * sin(tj * creal(modefreqs->data[i])) );
            rdwave2->data[j] += exp(- tj * cimag(modefreqs->data[i]))
            * (- modeamps->data[i] * sin(tj * creal(modefreqs->data[i]))
               +   modeamps->data[i + nmodes] * cos(tj * creal(modefreqs->data[i])) );
        }
    }

    XLALDestroyREAL8Vector(modeamps);
    return errcode;
}

/**
 * Function which calculates the value of the waveform, plus its
 * first and second derivatives, for the points which will be required
 * in the hybrid comb attachment of the ringdown.
 */
static INT4 XLALGenerateHybridWaveDerivatives (
                                               REAL8Vector	*rwave,      /**<< OUTPUT, values of the waveform at comb points */
                                               REAL8Vector	*dwave,      /**<< OUTPUT, 1st deriv of the waveform at comb points */
                                               REAL8Vector	*ddwave,     /**<< OUTPUT, 2nd deriv of the waveform at comb points */
                                               REAL8Vector	*timeVec,    /**<< Vector containing the time */
                                               REAL8Vector	*wave,       /**<< Last part of inspiral waveform */
                                               REAL8Vector	*matchrange, /**<< Times which determine the size of the comb */
                                               REAL8           dt,          /**<< Sample time step */
                                               REAL8           mass1,       /**<< First component mass (in Solar masses) */
                                               REAL8           mass2        /**<< Second component mass (in Solar masses) */
)
{

    /* XLAL error handling */
    INT4 errcode = XLAL_SUCCESS;

    /* For checking GSL return codes */
    INT4 gslStatus;

    UINT4 j;
    UINT4 vecLength;
    REAL8 m;
    double *y;
    double ry, dy, dy2;
    double rt;
    double *tlist;
    gsl_interp_accel *acc;
    gsl_spline *spline;

    /* Total mass in geometric units */
    m  = (mass1 + mass2) * LAL_MTSUN_SI;

    tlist = (double *) LALMalloc(6 * sizeof(double));
    rt = (matchrange->data[1] - matchrange->data[0]) / 5.;
    tlist[0] = matchrange->data[0];
    tlist[1] = tlist[0] + rt;
    tlist[2] = tlist[1] + rt;
    tlist[3] = tlist[2] + rt;
    tlist[4] = tlist[3] + rt;
    tlist[5] = matchrange->data[1];

    /* Set the length of the interpolation vectors */
    vecLength = ( m * matchrange->data[2] / dt ) + 1;

    /* Getting interpolation and derivatives of the waveform using gsl spline routine */
    /* Initiate arrays and supporting variables for gsl */
    y = (double *) LALMalloc(vecLength * sizeof(double));

    if ( !y )
    {
        XLAL_ERROR( XLAL_ENOMEM );
    }

    for (j = 0; j < vecLength; ++j)
    {
        y[j] = wave->data[j];
    }


    XLAL_CALLGSL( acc = (gsl_interp_accel*) gsl_interp_accel_alloc() );
    XLAL_CALLGSL( spline = (gsl_spline*) gsl_spline_alloc(gsl_interp_cspline, vecLength) );
    if ( !acc || !spline )
    {
        if ( acc )    gsl_interp_accel_free(acc);
        if ( spline ) gsl_spline_free(spline);
        LALFree( y );
        XLAL_ERROR( XLAL_ENOMEM );
    }

    /* Gall gsl spline interpolation */
    gslStatus = gsl_spline_init(spline, timeVec->data, y, vecLength);
    if ( gslStatus != GSL_SUCCESS )
    {
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
        LALFree( y );
        XLAL_ERROR( XLAL_EFUNC );
    }

    /* Getting first and second order time derivatives from gsl interpolations */
    for (j = 0; j < 6; ++j)
    {
        gslStatus = gsl_spline_eval_e( spline, tlist[j], acc, &ry );
        if ( gslStatus == GSL_SUCCESS )
        {
            gslStatus = gsl_spline_eval_deriv_e(spline, tlist[j], acc, &dy );
            gslStatus = gsl_spline_eval_deriv2_e(spline, tlist[j], acc, &dy2 );
        }
        if (gslStatus != GSL_SUCCESS )
        {
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            LALFree( y );
            XLAL_ERROR( XLAL_EFUNC );
        }
        rwave->data[j]  = (REAL8)(ry);
        dwave->data[j]  = (REAL8)(dy/m);
        ddwave->data[j] = (REAL8)(dy2/m/m);

    }

    /* Free gsl variables */
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    LALFree( tlist );
    LALFree(y);

    return errcode;
}

/**
 * The main workhorse function for performing the ringdown attachment for EOB
 * models EOBNRv2 and SEOBNRv1. This is the function which gets called by the
 * code generating the full IMR waveform once generation of the inspiral part
 * has been completed.
 * The ringdown is attached using the hybrid comb matching detailed in
 * The method is describe in Sec. II C of Pan et al. PRD 84, 124052 (2011),
 * specifically Eqs. 30 - 32.. Further details of the
 * implementation of the found in the DCC document T1100433.
 * In SEOBNRv1, the last physical overtone is replace by a pseudoQNM. See
 * Taracchini et al. PRD 86, 024011 (2012) for details.
 * STEP 1) Get mass and spin of the final black hole and the complex ringdown frequencies
 * STEP 2) Based on least-damped-mode decay time, allocate memory for rigndown waveform
 * STEP 3) Get values and derivatives of inspiral waveforms at matching comb points
 * STEP 4) Solve QNM coefficients and generate ringdown waveforms
 * STEP 5) Stitch inspiral and ringdown waveoforms
 * SEOBNRv2 prescriptions can be found in  https://dcc.ligo.org/T1400476
 * SEONBRv3 prescriptions can be found in XXX
 */
static INT4
XLALSimIMREOBHybridAttachRingdownPrec(
				  REAL8Vector * signal1,	/**<< OUTPUT, Real of inspiral waveform to which we attach ringdown */
				  REAL8Vector * signal2,	/**<< OUTPUT, Imag of inspiral waveform to which we attach ringdown */
				  const INT4 l,	/**<< Current mode l */
				  const INT4 m,	/**<< Current mode m */
				  const REAL8 dt,	/**<< Sample time step (in seconds) */
				  const REAL8 mass1,	/**<< First component mass (in Solar masses) */
				  const REAL8 mass2,	/**<< Second component mass (in Solar masses) */
				  const REAL8 spin1x,	/**<<The spin of the first object; only needed for spin waveforms */
				  const REAL8 spin1y,	/**<<The spin of the first object; only needed for spin waveforms */
				  const REAL8 spin1z,	/**<<The spin of the first object; only needed for spin waveforms */
				  const REAL8 spin2x,	/**<<The spin of the second object; only needed for spin waveforms */
				  const REAL8 spin2y,	/**<<The spin of the second object; only needed for spin waveforms */
				  const REAL8 spin2z,	/**<<The spin of the second object; only needed for spin waveforms */
				  REAL8Vector * timeVec,	/**<< Vector containing the time values */
				  REAL8Vector * matchrange,	/**<< Time values chosen as points for performing comb matching */
				  Approximant approximant,	/**<<The waveform approximant being used */
                  const REAL8 JLN           /**<< cosine of the angle between J and LN at the light ring */
)
{
	INT4		debugout = 0;

	COMPLEX16Vector *modefreqs;
	//COMPLEX16 freq7sav;
	UINT4		Nrdwave;
	UINT4		j;

	UINT4		nmodes;
	REAL8Vector    *rdwave1;
	REAL8Vector    *rdwave2;
	REAL8Vector    *rinspwave;
	REAL8Vector    *dinspwave;
	REAL8Vector    *ddinspwave;
	REAL8VectorSequence *inspwaves1;
	REAL8VectorSequence *inspwaves2;
	REAL8		eta     , a, chi, NRPeakOmega22;	/* To generate pQNM
								 * frequency */
	REAL8		sh      , kk, kt1, kt2;	/* To generate pQNM frequency */
	REAL8		mTot;	/* In geometric units */
	REAL8		spin1    [3] = {spin1x, spin1y, spin1z};
	REAL8		spin2    [3] = {spin2x, spin2y, spin2z};
	REAL8		finalMass, finalSpin;
	REAL8		chi1    , chi2, theta1, theta2;
    INT4 mHere = 0;

    mTot = (mass1 + mass2) * LAL_MTSUN_SI;
	eta = mass1 * mass2 / ((mass1 + mass2) * (mass1 + mass2));

	/*
	 * STEP 1) Get mass and spin of the final black hole and the complex
	 * ringdown frequencies
	 */

	/* Create memory for the QNM frequencies */
	nmodes = 8;
	modefreqs = XLALCreateCOMPLEX16Vector(nmodes);
	if (!modefreqs) {
		XLAL_ERROR(XLAL_ENOMEM);
	}
	if (XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs, mass1, mass2, spin1, spin2, l, m, nmodes, approximant) == XLAL_FAILURE) {
		XLALDestroyCOMPLEX16Vector(modefreqs);
		XLAL_ERROR(XLAL_EFUNC);
	}
	if (approximant == SEOBNRv3) {

        if (JLN > 0){
            mHere = (int)fabs((REAL8) m);
        }else{
            mHere = -(int)fabs((REAL8) m);
        }

		if (XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs, mass1, mass2, spin1, spin2, l, (int)fabs((REAL8) m), nmodes, approximant) == XLAL_FAILURE) {
			XLALDestroyCOMPLEX16Vector(modefreqs);
			XLAL_ERROR(XLAL_EFUNC);
		}
		// FIXME
        if (m < 0 && JLN >= 0) {
			for (j = 0; j < nmodes; j++) {
				modefreqs->data[j] = conjl(-1.0 * modefreqs->data[j]);
			}
		}
        if (m > 0 && JLN < 0) {
            for (j = 0; j < nmodes; j++) {
                modefreqs->data[j] = conjl(-1.0 * modefreqs->data[j]);
            }
        }
        if (m == 0) {
			for (j = 5; j < nmodes; j++) {
				modefreqs->data[j] = conjl(-1.0 * modefreqs->data[j - 5]);
			}
		}
	}

	/*
	 * Call XLALSimIMREOBFinalMassSpinPrec() to get mass and spin of the
	 * final black hole
	 */
	if (XLALSimIMREOBFinalMassSpinPrec(&finalMass, &finalSpin, mass1, mass2, spin1, spin2, approximant) == XLAL_FAILURE) {
		XLAL_ERROR(XLAL_EFUNC);
	}
	if (debugout) {
		XLAL_PRINT_INFO("(Mf, af) = (%3.10f, %3.10f)\n", finalMass, finalSpin);
	}
	if (approximant == SEOBNRv1) {
		/* Replace the last QNM with pQNM */
		/* We assume aligned/antialigned spins here */
		a = (spin1[2] + spin2[2]) / 2. * (1.0 - 2.0 * eta) + (spin1[2] - spin2[2]) / 2. * (mass1 - mass2) / (mass1 + mass2);
		NRPeakOmega22 = GetNRSpinPeakOmega(l, m, eta, a) / mTot;
		/*
		 * XLAL_PRINT_INFO("a and NRomega in QNM freq: %.16e %.16e %.16e %.16e
		 * %.16e\n",spin1[2],spin2[2],
		 * mTot/LAL_MTSUN_SI,a,NRPeakOmega22*mTot);
		 */
		modefreqs->data[7] = (NRPeakOmega22 / finalMass + creal(modefreqs->data[0])) / 2.;
		modefreqs->data[7] += I * 10. / 3. * cimag(modefreqs->data[0]);
	}
	if (approximant == SEOBNRv2)
		//See pages 6 to 12 of the dcc document T1400476 - v3 for expressions
		//	in this block.
		{
			/* Replace the last two QNMs with pQNMs */
			/* We assume aligned/antialigned spins here */
			/*
			 * Definitions of a, chi and NRPeakOmega22, where the
			 * last one is an approximation of \phi'[tmatch] in
			 * T1400476-v3.
			 */
			a = (spin1[2] + spin2[2]) / 2. * (1.0 - 2.0 * eta) + (spin1[2] - spin2[2]) / 2. * (mass1 - mass2) / (mass1 + mass2);
			NRPeakOmega22 = GetNRSpinPeakOmegav2(l, m, eta, a) / mTot;

			/* Define chi */
			chi = (spin1[2] + spin2[2]) / 2. + (spin1[2] - spin2[2]) / 2. * ((mass1 - mass2) / (mass1 + mass2)) / (1. - 2. * eta);

			/*
			 * For extreme chi (>= 0.8), there are scale factors
			 * in both complex pseudo-QNM frequencies. kk, kt1,
			 * kt2 describe those factors.
			 */
			//Below definitions of kk, kt1 and kt2 are likely obsolete
				kk = kt1 = kt2 = 1.;
			if (chi >= 0.8) {
				kk = 0.7 + 0.3 * exp(100. * (eta - 0.25));
				kt1 = 0.5 * sqrt(1. + 800.0 * eta * eta / 3.0) - 0.125;
				kt2 = 0.5 * pow(1. + 0.5 * eta * sqrt(eta) / 0.0225, 2. / 3.) - 0.2;
			}
			//Above definitions of kk, kt1 and kt2 are likely obsolete
			/*
			 * XLAL_PRINT_INFO("a, chi and NRomega in QNM freq: %.16e
			 * %.16e %.16e %.16e %.16e %.16e\n",
			 * spin1[2],spin2[2],mTot/LAL_MTSUN_SI,a,chi,NRPeakOme
			 * ga22*mTot);
			 */
				sh = 0.;
			//freq7sav = modefreqs->data[7];

			/*
			 * Cases 1, 2 and 3 in T1400476-v3. Note that the
			 * difference between the chi1=chi2=0 case and the
			 * chi<0.7 cases is only in Dtcomb, which is not
			 * specified or used in this file.
			 */
			modefreqs->data[7] = (2. / 3. * NRPeakOmega22 / finalMass) + (1. / 3. * creal(modefreqs->data[0]));
			modefreqs->data[7] += I * 3.5 / 0.9 * cimag(modefreqs->data[0]);
			modefreqs->data[6] = (3. / 4. * NRPeakOmega22 / finalMass) + (1. / 4. * creal(modefreqs->data[0]));
			modefreqs->data[6] += I * 3.5 * cimag(modefreqs->data[0]);

			/*
			 * For extreme chi (>= 0.8), the ringdown attachment
			 * should end slightly before the value given passed
			 * to this function. sh gives exactly how long
			 * before.
			 */
			if (chi >= 0.7 && chi < 0.8) {
				sh = -9. * (eta - 0.25);
			}
			if ((eta > 30. / 31. / 31. && eta <= 10. / 121. && chi >= 0.8) || (eta <= 30. / 31. / 31. && chi >= 0.8 && chi < 0.9)) {
				//This is case 4 in T1400476 - v3
					sh = -9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
				kk = 0.7 + 0.3 * exp(100. * (eta - 0.25));
				kt1 = 0.5 * sqrt(1. + 800.0 * eta * eta / 3.0) - 0.125;
				kt2 = 0.5 * pow(1. + 0.5 * eta * sqrt(eta) / 0.0225, 2. / 3.) - 0.2;
				modefreqs->data[4] = 0.4 * (1. + kk) * creal(modefreqs->data[6])
					+ I * cimag(modefreqs->data[6]) / (2.5 * kt2 * exp(-(eta - 0.005) / 0.03));
				modefreqs->data[5] = 0.4 * (1. + kk) * creal(modefreqs->data[7])
					+ I * cimag(modefreqs->data[7]) / (1.5 * kt1 * exp(-(eta - 0.005) / 0.03));
				modefreqs->data[6] = kk * creal(modefreqs->data[6]) + I * cimag(modefreqs->data[6]) / kt2;
				modefreqs->data[7] = kk * creal(modefreqs->data[7]) + I * cimag(modefreqs->data[7]) / kt1;
			}
			if (eta < 30. / 31. / 31. && chi >= 0.9) {
				//This is case 5 in T1400476 - v3
					sh = 0.55 - 9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
				kk = 0.7 + 0.3 * exp(100. * (eta - 0.25));
				kt1 = 0.5 * sqrt(1. + 800.0 * eta * eta / 3.0) - 0.125;
				kt2 = 0.5 * pow(1. + 0.5 * eta * sqrt(eta) / 0.0225, 2. / 3.) - 0.2;
				modefreqs->data[4] = 1.1 * 0.4 * (1. + kk) * creal(modefreqs->data[6])
					+ I * cimag(modefreqs->data[6]) / (1.05 * 2.5 * kt2 * exp(-(eta - 0.005) / 0.03));
				modefreqs->data[5] = 0.4 * (1. + kk) * creal(modefreqs->data[7])
					+ I * cimag(modefreqs->data[7]) / (1.05 * 1.5 * kt1 * exp(-(eta - 0.005) / 0.03));
				modefreqs->data[6] = kk * creal(modefreqs->data[6]) + I * cimag(modefreqs->data[6]) / kt2;
				modefreqs->data[7] = kk * creal(modefreqs->data[7]) + I * cimag(modefreqs->data[7]) / kt1;
			}
			if (eta > 10. / 121. && chi >= 0.8) {
				//This is case 6 in T1400476 - v3
					sh = 1. - 9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
				kk = 0.7 + 0.3 * exp(100. * (eta - 0.25));
				kt1 = 0.45 * sqrt(1. + 200.0 * eta * eta / 3.0) - 0.125;
				kt2 = 0.5 * pow(1. + 0.5 * eta * sqrt(eta) / 0.0225, 2. / 3.) - 0.2;
				modefreqs->data[6] = kk * creal(modefreqs->data[6]) + I * cimag(modefreqs->data[6]) / 0.95 / kt2;
				modefreqs->data[7] = kk * creal(modefreqs->data[7]) + I * cimag(modefreqs->data[7]) / kt1;
			}
			//The last line of T1400476 - v3
				matchrange->data[0] -= sh;
			matchrange->data[1] -= sh;
			/*
			 * modefreqs->data[7] = 0.38068371/mTot +
			 * I/1.4677128/mTot; modefreqs->data[6] =
			 * 0.37007703/mTot + I/1.3359367/mTot;
			 * modefreqs->data[5] = 0.36980703/mTot +
			 * I/1.7791212/mTot; modefreqs->data[4] =
			 * 0.3595034/mTot + I/2.6989764/mTot; XLAL_PRINT_INFO("sh =
			 * %f\n",sh); XLAL_PRINT_INFO("PeakOmega = %f, mTot =
			 * %f\n",NRPeakOmega22,mTot); XLAL_PRINT_INFO("w0 = %f, t0 =
			 * %f\n",creal(modefreqs->data[0])*mTot,
			 * 1./cimag(modefreqs->data[0])/mTot); XLAL_PRINT_INFO("w1 =
			 * %f, t1 = %f\n",creal(modefreqs->data[6])*mTot,
			 * 1./cimag(modefreqs->data[6])/mTot); XLAL_PRINT_INFO("w2 =
			 * %f, t2 = %f\n",creal(modefreqs->data[7])*mTot,
			 * 1./cimag(modefreqs->data[7])/mTot); XLAL_PRINT_INFO("w3 =
			 * %f, t3 = %f\n",creal(modefreqs->data[4])*mTot,
			 * 1./cimag(modefreqs->data[4])/mTot); XLAL_PRINT_INFO("w4 =
			 * %f, t4 = %f\n",creal(modefreqs->data[5])*mTot,
			 * 1./cimag(modefreqs->data[5])/mTot);
			 */
		}
	if (approximant == SEOBNRv3) {
        REAL8 kappa_thr = 0.175;
        //REAL8 eJL_thr = 7.5e-3;
        REAL8 eJL_thr = 5.0e-2;
        //REAL8 eJL_thr = 1.0e-2;
		chi1 = sqrt(spin1[0] * spin1[0] + spin1[1] * spin1[1] + spin1[2] * spin1[2]);
		chi2 = sqrt(spin2[0] * spin2[0] + spin2[1] * spin2[1] + spin2[2] * spin2[2]);
		if (chi1 < 1.0e-15) {
			theta1 = 0.;
		} else {
			theta1 = acos(spin1[2] / chi1);
		}
		if (chi2 < 1.0e-15) {
			theta2 = 0.;
		} else {
			theta2 = acos(spin2[2] / chi2);
		}
		chi1 = chi1 * cos(theta1);
		chi2 = chi2 * cos(theta2);

        /*Compute wf freq at matching point*/
        double *y;
        double ry, dy;
        gsl_interp_accel *acc;
        gsl_spline     *spline;
        //UINT4 vecLength = ((mass1 + mass2) * LAL_MTSUN_SI * matchrange->data[2] / dt) + 1;
        //y = (double *)LALMalloc(vecLength * sizeof(double));
        y = (double *)LALMalloc(timeVec->length * sizeof(double));
        double hRe, dhRe, hIm, dhIm;

        if (!y) {
            XLAL_ERROR(XLAL_ENOMEM);
        }
        //FILE *out1 = fopen( "Andrea1.dat","w");
        for (j = 0; j < timeVec->length; ++j) {
            y[j] = signal1->data[j];
            //fprintf(out1, "%.16e %.16e\n", timeVec->data[j], y[j]);
            //XLAL_PRINT_INFO("%.16e %.16e\n", timeVec->data[j], y[j]);
        }
        //fclose(out1);
        //exit(0);

        XLAL_CALLGSL(acc = (gsl_interp_accel *) gsl_interp_accel_alloc());
        //XLAL_CALLGSL(spline = (gsl_spline *) gsl_spline_alloc(gsl_interp_cspline, vecLength));
        XLAL_CALLGSL(spline = (gsl_spline *) gsl_spline_alloc(gsl_interp_cspline, timeVec->length));
        if (!acc || !spline) {
            if (acc)
                gsl_interp_accel_free(acc);
            if (spline)
                gsl_spline_free(spline);
            LALFree(y);
            XLAL_ERROR(XLAL_ENOMEM);
        }

        //XLAL_PRINT_INFO("here1....  last point %f, match_t %f \n", timeVec->data[timeVec->length-1], matchrange->data[1]);
        //INT4 gslStatus = gsl_spline_init(spline, timeVec->data, y, vecLength);
        INT4 gslStatus = gsl_spline_init(spline, timeVec->data, y, timeVec->length);
        if (gslStatus != GSL_SUCCESS) {
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            LALFree(y);
            XLAL_ERROR(XLAL_EFUNC);
        }
        gslStatus = gsl_spline_eval_e(spline, matchrange->data[1], acc, &ry);
        //int ind_stas = (int) matchrange->data[1]*(((mass1 + mass2) * LAL_MTSUN_SI / dt));
        //XLAL_PRINT_INFO("here2.... %f,  %f,  %f, %f \n", ry, timeVec->data[ind_stas], matchrange->data[1], y[ind_stas]);
        if (gslStatus == GSL_SUCCESS) {
            gslStatus = gsl_spline_eval_deriv_e(spline, matchrange->data[1], acc, &dy);
            //XLAL_PRINT_INFO("here2.... ");
        }
        if (gslStatus != GSL_SUCCESS) {
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            LALFree(y);
            XLAL_ERROR(XLAL_EFUNC);
        }
        hRe = (REAL8) (ry);
        dhRe = (REAL8) (dy);

        //XLAL_PRINT_INFO("here hRe = %f, dhRe = %f \n", hRe, dhRe);

        //FILE *out2 = fopen( "Andrea2.dat","w");
        //for (j = 0; j < vecLength; ++j) {
        for (j = 0; j < timeVec->length; ++j) {
            y[j] = signal2->data[j];
            //fprintf(out2, "%.16e %.16e\n", timeVec->data[j], y[j]);
        }
        //fclose(out2);
        //XLAL_CALLGSL(acc = (gsl_interp_accel *) gsl_interp_accel_alloc());
        //XLAL_CALLGSL(spline = (gsl_spline *) gsl_spline_alloc(gsl_interp_cspline, vecLength));
        //XLAL_CALLGSL(spline = (gsl_spline *) gsl_spline_alloc(gsl_interp_cspline, signal2->length));
        if (!acc || !spline) {
            if (acc)
                gsl_interp_accel_free(acc);
            if (spline)
                gsl_spline_free(spline);
            LALFree(y);
            XLAL_ERROR(XLAL_ENOMEM);
        }
        //gslStatus = gsl_spline_init(spline, timeVec->data, y, vecLength);
        gslStatus = gsl_spline_init(spline, timeVec->data, y, timeVec->length);
        if (gslStatus != GSL_SUCCESS) {
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            LALFree(y);
            XLAL_ERROR(XLAL_EFUNC);
        }
        gslStatus = gsl_spline_eval_e(spline, matchrange->data[1], acc, &ry);
        if (gslStatus == GSL_SUCCESS) {
            gslStatus = gsl_spline_eval_deriv_e(spline, matchrange->data[1], acc, &dy);
        }
        if (gslStatus != GSL_SUCCESS) {
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            LALFree(y);
            XLAL_ERROR(XLAL_EFUNC);
        }
        hIm = (REAL8) (ry);
        dhIm = (REAL8) (dy);

        //XLAL_PRINT_INFO("Stas, check hRe = %f, dhRe = %f, hIm = %f, dhIm = %f \n", hRe, dhRe, hIm, dhIm);
        double hNorm2 = hRe*hRe + hIm*hIm;
        double omegaWavePeak = creal(modefreqs->data[0]);
        if (hNorm2 != 0.0){
            omegaWavePeak = ((-dhRe*hIm + dhIm*hRe)/hNorm2) / mTot;
        }
        else{
            if(debugout) {XLAL_PRINT_INFO("hNorm=0.0 for (l,m)=(%d,%d)\n",l,m);}
//            XLAL_ERROR(XLAL_EFAILED);
        }


		a = (chi1 + chi2) / 2. * (1.0 - 2.0 * eta) + (chi1 - chi2) / 2. * (mass1 - mass2) / (mass1 + mass2);
        NRPeakOmega22 = fabs(omegaWavePeak);
//        NRPeakOmega22 = GetNRSpinPeakOmegav2(l, m, eta, a) / mTot;
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
        LALFree(y);
        // FIXME
        //NRPeakOmega22 = 0.3;
		//NRPeakOmega22 = GetNRSpinPeakOmegav2(l, m, eta, a) / mTot;
//        NRPeakOmega22 = omegaWavePeak/mTot;
//        XLAL_PRINT_INFO("(hRe, hIm, dhRe, dhIm)=(%.16e, %.16e, %.16e, %.16e)\n", hRe, hIm, dhRe, dhIm);
//        XLAL_PRINT_INFO("hNorm %.16e, tmatch %.16e\n",sqrt(hNorm2), matchrange->data[1]);
//        XLAL_PRINT_INFO("NRPeakOmega22 %.16e, omegaWavePeak %.16e\n",NRPeakOmega22,omegaWavePeak);
		chi = (chi1 + chi2) / 2. + (chi1 - chi2) / 2. * ((mass1 - mass2) / (mass1 + mass2)) / (1. - 2. * eta);

		sh = 0.;
		if (chi >= 0.7 && chi < 0.8) {
			sh = -9. * (eta - 0.25);
		}
		if ((eta > 30. / 31. / 31. && eta <= 10. / 121. && chi >= 0.8) || (eta <= 30. / 31. / 31. && chi >= 0.8 && chi < 0.9)) {
			//This is case 4 in T1400476 - v3
				sh = -9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
		}
		if (eta < 30. / 31. / 31. && chi >= 0.9) {
			//This is case 5 in T1400476 - v3
				sh = 0.55 - 9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
		}
		if (eta > 10. / 121. && chi >= 0.8) {
			//This is case 6 in T1400476 - v3
				sh = 1. - 9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
		}
		if ((int)fabs((REAL8) m) == 2 && l == 2) {
			/* if (m == 2 && l ==2){ */

			/*
			 * For extreme chi (>= 0.8), there are scale factors
			 * in both complex pseudo-QNM frequencies. kk, kt1,
			 * kt2 describe those factors.
			 */
			/*
			 * XLAL_PRINT_INFO("Stas: a, chi and NRomega in QNM freq:
			 * %.16e %.16e %.16e %.16e %.16e %.16e\n",
			 * spin1[2],spin2[2],mTot/LAL_MTSUN_SI,a,chi,NRPeakOme
			 * ga22*mTot);
			 */
			kk = kt1 = kt2 = 1.;
			if (chi >= 0.8) {
				kk = 0.7 + 0.3 * exp(100. * (eta - 0.25));
				kt1 = -0.125 + sqrt(1. + 200. * pow(eta, 2) / 3.) / 2.;
				kt2 = -0.2 + pow(1. + 200. * pow(eta, 3. / 2.) / 9., 2. / 3.) / 2.;
			}
			//Computing pQNMs
            // FIXME
			//if (m < 0) {
			//if (m > 0) {
            //printf("Stas: J.L = %f, J.L*m = %f  Re(mode) = %f \n", JLN, JLN*m, creal(modefreqs->data[0]));
            if (JLN*m < 0.0){
				modefreqs->data[6] = (-3. / 4. * NRPeakOmega22 / finalMass) + (1. / 4. * creal(modefreqs->data[0]));
				modefreqs->data[7] = (-2. / 3. * NRPeakOmega22 / finalMass) + (1. / 3. * creal(modefreqs->data[0]));
			} else {
				modefreqs->data[6] = (3. / 4. * NRPeakOmega22 / finalMass) + (1. / 4. * creal(modefreqs->data[0]));
				modefreqs->data[7] = (2. / 3. * NRPeakOmega22 / finalMass) + (1. / 3. * creal(modefreqs->data[0]));
                //XLAL_PRINT_INFO("Stas , here m = %d modes: 0 = %4.10f, 6 = %4.10f, 7 = %4.10f \n", m, creal(modefreqs->data[0])*mTot, creal(modefreqs->data[6])*mTot, creal(modefreqs->data[7])*mTot);
			}

			modefreqs->data[7] += I * 3.5 / 0.9 * cimag(modefreqs->data[0]);
			modefreqs->data[6] += I * 3.5 * cimag(modefreqs->data[0]);

			if ((eta > 30. / 31. / 31. && eta <= 10. / 121. && chi >= 0.8) || (eta <= 30. / 31. / 31. && chi >= 0.8 && chi < 0.9)) {
				//This is case 4 in T1400476 - v3
                if (debugout)
                    XLAL_PRINT_INFO("Stas, this is case4 for pQNM eta = %f, chi = %f \n", eta, chi);
					sh = -9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
				kk = 0.7 + 0.3 * exp(100. * (eta - 0.25));
				kt1 = 0.5 * sqrt(1. + 800.0 * eta * eta / 3.0) - 0.125;
				kt2 = 0.5 * pow(1. + 0.5 * eta * sqrt(eta) / 0.0225, 2. / 3.) - 0.2;

				modefreqs->data[4] = 0.4 * (1. + kk) * creal(modefreqs->data[6])
					+ I * cimag(modefreqs->data[6]) / (2.5 * kt2 * exp(-(eta - 0.005) / 0.03));
				modefreqs->data[5] = 0.4 * (1. + kk) * creal(modefreqs->data[7])
					+ I * cimag(modefreqs->data[7]) / (1.5 * kt1 * exp(-(eta - 0.005) / 0.03));
				modefreqs->data[6] = kk * creal(modefreqs->data[6]) + I * cimag(modefreqs->data[6]) / kt2;
				modefreqs->data[7] = kk * creal(modefreqs->data[7]) + I * cimag(modefreqs->data[7]) / kt1;
			}
			if (eta < 30. / 31. / 31. && chi >= 0.9) {
				//This is case 5 in T1400476 - v3
                if (debugout)
                    XLAL_PRINT_INFO("Stas, this is case5 for pQNM eta = %f, chi = %f \n", eta, chi);
					sh = 0.55 - 9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
				kk = 0.7 + 0.3 * exp(100. * (eta - 0.25));
				kt1 = 0.5 * sqrt(1. + 800.0 * eta * eta / 3.0) - 0.125;
				kt2 = 0.5 * pow(1. + 0.5 * eta * sqrt(eta) / 0.0225, 2. / 3.) - 0.2;
				modefreqs->data[4] = 1.1 * 0.4 * (1. + kk) * creal(modefreqs->data[6])
					+ I * cimag(modefreqs->data[6]) / (1.05 * 2.5 * kt2 * exp(-(eta - 0.005) / 0.03));
				modefreqs->data[5] = 0.4 * (1. + kk) * creal(modefreqs->data[7])
					+ I * cimag(modefreqs->data[7]) / (1.05 * 1.5 * kt1 * exp(-(eta - 0.005) / 0.03));
				modefreqs->data[6] = kk * creal(modefreqs->data[6]) + I * cimag(modefreqs->data[6]) / kt2;
				modefreqs->data[7] = kk * creal(modefreqs->data[7]) + I * cimag(modefreqs->data[7]) / kt1;
			}
			if (eta > 10. / 121. && chi >= 0.8) {
				//This is case 6 in T1400476 - v3
                if (debugout)
                    XLAL_PRINT_INFO("Stas, this is case6 for pQNM eta = %f, chi = %f \n", eta, chi);
					sh = 1. - 9. * (eta - 0.25) * (1. + 2. * exp(-(chi - 0.85) * (chi - 0.85) / 0.05 / 0.05)) * (1. + 1. / (1. + exp((eta - 0.01) / 0.001)));
				kk = 0.7 + 0.3 * exp(100. * (eta - 0.25));
				kt1 = 0.45 * sqrt(1. + 200.0 * eta * eta / 3.0) - 0.125;
				kt2 = 0.5 * pow(1. + 0.5 * eta * sqrt(eta) / 0.0225, 2. / 3.) - 0.2;
				modefreqs->data[6] = kk * creal(modefreqs->data[6]) + I * cimag(modefreqs->data[6]) / 0.95 / kt2;
				modefreqs->data[7] = kk * creal(modefreqs->data[7]) + I * cimag(modefreqs->data[7]) / kt1;
			}
            if (debugout)
                XLAL_PRINT_INFO("Stas, default (if nothing specified above) for pQNM eta = %f, chi = %f chi1 = %f \n", eta, chi,  chi1);

            // FIXME It is not compatible with v2!!!!
//            if (chi1 >= 0.96){//  && eta < 10.0/11./11.){
//            //    modefreqs->data[7] += I * 3.5 / 0.9 * cimag(modefreqs->data[0]) - I * cimag(modefreqs->data[7]);
//                if (debugout){
//                    XLAL_PRINT_INFO("Stas, the decay time will be modified for pQNM  chi1 = %f, from %f to %f \n", chi1, 1.0/cimag(modefreqs->data[7])/mTot,   1.0/(5.0 * cimag(modefreqs->data[0]))/mTot );
//                }
//                  modefreqs->data[7] += I * 5.0 * cimag(modefreqs->data[0]) - I * cimag(modefreqs->data[7]);
//            }
			//The last line of T1400476 - v3
//				matchrange->data[0] -= sh;
//			matchrange->data[1] -= sh;

			/*
			 * This is a v1 of pQNM in RD attachment if (m<0){
			 * modefreqs->data[7] = (-NRPeakOmega22/finalMass +
			 * creal(modefreqs->data[0])) / 2.;       } else {
			 * modefreqs->data[7] = (NRPeakOmega22/finalMass +
			 * creal(modefreqs->data[0])) / 2.; }
			 * modefreqs->data[7] += I * 10./3. *
			 * cimag(modefreqs->data[0]);
			 */
//            if(m==-2 || m==-2) {
//            modefreqs->data[6] = 0.495264/mTot + I/mTot/3.716525;
//            modefreqs->data[7] = 0.514602/mTot + I/mTot/3.344873;
//            }
//            else if (m==-2){
//                modefreqs->data[6] = -0.495264/mTot + I/mTot/3.716525;
//                modefreqs->data[7] = -0.514602/mTot + I/mTot/3.344873;
//            }
            
//            XLAL_PRINT_INFO("finalSpin = %f\n",finalSpin);
//            XLAL_PRINT_INFO("finalMass = %f\n",finalMass);
//            XLAL_PRINT_INFO("PeakOmega = %f\n",NRPeakOmega22*mTot);
//            XLAL_PRINT_INFO("w0 = %f, t0 = %f\n",creal(modefreqs->data[0])*mTot, 1./cimag(modefreqs->data[0])/mTot);
//            XLAL_PRINT_INFO("w1 = %f, t1 = %f\n",creal(modefreqs->data[1])*mTot, 1./cimag(modefreqs->data[1])/mTot);
//            XLAL_PRINT_INFO("w2 = %f, t2 = %f\n",creal(modefreqs->data[2])*mTot, 1./cimag(modefreqs->data[2])/mTot);
//            XLAL_PRINT_INFO("w3 = %f, t3 = %f\n",creal(modefreqs->data[3])*mTot, 1./cimag(modefreqs->data[3])/mTot);
//            XLAL_PRINT_INFO("w4 = %f, t4 = %f\n",creal(modefreqs->data[4])*mTot, 1./cimag(modefreqs->data[4])/mTot);
//            XLAL_PRINT_INFO("w5 = %f, t5 = %f\n",creal(modefreqs->data[5])*mTot, 1./cimag(modefreqs->data[5])/mTot);
//            XLAL_PRINT_INFO("w6 = %f, t6 = %f\n",creal(modefreqs->data[6])*mTot, 1./cimag(modefreqs->data[6])/mTot);
//            XLAL_PRINT_INFO("w7 = %f, t7 = %f\n",creal(modefreqs->data[7])*mTot, 1./cimag(modefreqs->data[7])/mTot);

    

             // FIXME
             // {{{
            if (debugout) {
                XLAL_PRINT_INFO("Stas: JLN*eta = %f, \n", JLN*eta);
                XLAL_PRINT_INFO("NRPeakOmega22 = %3.10f,  %3.10f,  %3.10f,  %3.10f\n", NRPeakOmega22 * mTot /finalMass, NRPeakOmega22 / finalMass,  finalMass,  mTot);
                for (j = 0; j < nmodes; j++) {
                    XLAL_PRINT_INFO("QNM frequencies: %d %d %d %3.10f %3.10f\n", l, m, j, creal(modefreqs->data[j]) * mTot, 1. / cimag(modefreqs->data[j]) / mTot);
                }
            }

             /*if (JLN*m < 0.0){
				modefreqs->data[5] = (-3. / 4. * NRPeakOmega22 / finalMass) + (1. / 4. * creal(modefreqs->data[0]));
				modefreqs->data[4] = (-2. / 3. * NRPeakOmega22 / finalMass) + (1. / 3. * creal(modefreqs->data[0]));
			 } else {
				modefreqs->data[5] = (3. / 4. * NRPeakOmega22 / finalMass) + (1. / 4. * creal(modefreqs->data[0]));
				modefreqs->data[4] = (2. / 3. * NRPeakOmega22 / finalMass) + (1. / 3. * creal(modefreqs->data[0]));
                //XLAL_PRINT_INFO("Stas , here m = %d modes: 0 = %4.10f, 6 = %4.10f, 7 = %4.10f \n", m, creal(modefreqs->data[0])*mTot, creal(modefreqs->data[6])*mTot, creal(modefreqs->data[7])*mTot);
			 }
			 modefreqs->data[5] += I * 3.5 / 0.9 * cimag(modefreqs->data[0]);
			 modefreqs->data[4] += I * 3.5 * cimag(modefreqs->data[0]);*/

	         COMPLEX16Vector *modefreqs_xtr;
             modefreqs_xtr = XLALCreateCOMPLEX16Vector(nmodes);
             if (!modefreqs_xtr) {
                 XLAL_ERROR(XLAL_ENOMEM);
             }
             //modefreqs->data[4] =    ;
             //modefreqs->data[5] =  0.0;
             //modefreqs->data[6] =   0.0;
             //modefreqs->data[7] =  0.0;
             //XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, 1, nmodes, approximant);
             //if (m == -2){
             //    XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, -1, nmodes, approximant);
             //}
             //modefreqs->data[4] = modefreqs_xtr->data[0];
             //modefreqs->data[5] = modefreqs_xtr->data[1];

/*
            XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, -2, nmodes, approximant);
             if ((JLN > 0.0 && JLN < 0.98) || (eta*JLN < eJL_thr && eta*JLN>0.0)){
                 modefreqs->data[5] =  modefreqs_xtr->data[0];
                 if (JLN < kappa_thr || eta*JLN < eJL_thr){
                     modefreqs->data[6] =  modefreqs_xtr->data[1];
                 }
                 if (m == -2){
                     modefreqs->data[5] =  conjl(-1.0 * modefreqs_xtr->data[0]);
                     if (JLN < kappa_thr || eta*JLN < eJL_thr){
                         modefreqs->data[6] =  conjl(-1.0 *  modefreqs_xtr->data[1]);
                     }
                 }
             }
 */
            if ((JLN > 0.0 && JLN < 0.98) || (eta*JLN < eJL_thr && eta*JLN>0.0)){
                                 XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, -2, nmodes, approximant);
                                 modefreqs->data[5] = modefreqs_xtr->data[0];
                                 if (m == -2){
                                         modefreqs->data[5] =  conjl(-1.0 * modefreqs_xtr->data[0]);
                                     }
                             }
 
             if ((JLN < 0.0 && JLN > -0.98) || (eta*JLN > -eJL_thr && eta*JLN<0.0) ){
                 spin1[0] *= -1;
                 spin1[1] *= -1;
                 spin1[2] *= -1;
                 spin2[0] *= -1;
                 spin2[1] *= -1;
                 spin2[2] *= -1;
                 XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l,  -mHere, nmodes, approximant);
                 spin1[0] *= -1;
                 spin1[1] *= -1;
                 spin1[2] *= -1;
                 spin2[0] *= -1;
                 spin2[1] *= -1;
                 spin2[2] *= -1;
                 modefreqs->data[5] =  modefreqs_xtr->data[0];
//                 if (JLN > -1.0*kappa_thr || fabs(eta*JLN) < eJL_thr){
//                     modefreqs->data[6] =  modefreqs_xtr->data[1];
//                 }
                 if (m == -2){
                     modefreqs->data[5] =  conjl(-1.0 * modefreqs->data[5]);
//                     if (JLN > -1.0*kappa_thr || fabs(eta*JLN) < eJL_thr){
//                        modefreqs->data[6] =  conjl(-1.0 * modefreqs_xtr->data[1]);
 //                    }
                 }
             }
             XLALDestroyCOMPLEX16Vector(modefreqs_xtr);
	if (debugout) {
            XLAL_PRINT_INFO("l,m = %d %d\n",l,m);
            XLAL_PRINT_INFO("finalSpin = %f\n",finalSpin);
            XLAL_PRINT_INFO("finalMass = %f\n",finalMass);
            XLAL_PRINT_INFO("PeakOmega = %f\n",NRPeakOmega22*mTot);
            XLAL_PRINT_INFO("w0 = %f, t0 = %f\n",creal(modefreqs->data[0])*mTot, 1./cimag(modefreqs->data[0])/mTot);
            XLAL_PRINT_INFO("w1 = %f, t1 = %f\n",creal(modefreqs->data[1])*mTot, 1./cimag(modefreqs->data[1])/mTot);
            XLAL_PRINT_INFO("w2 = %f, t2 = %f\n",creal(modefreqs->data[2])*mTot, 1./cimag(modefreqs->data[2])/mTot);
            XLAL_PRINT_INFO("w3 = %f, t3 = %f\n",creal(modefreqs->data[3])*mTot, 1./cimag(modefreqs->data[3])/mTot);
            XLAL_PRINT_INFO("w4 = %f, t4 = %f\n",creal(modefreqs->data[4])*mTot, 1./cimag(modefreqs->data[4])/mTot);
            XLAL_PRINT_INFO("w5 = %f, t5 = %f\n",creal(modefreqs->data[5])*mTot, 1./cimag(modefreqs->data[5])/mTot);
            XLAL_PRINT_INFO("w6 = %f, t6 = %f\n",creal(modefreqs->data[6])*mTot, 1./cimag(modefreqs->data[6])/mTot);
            XLAL_PRINT_INFO("w7 = %f, t7 = %f\n",creal(modefreqs->data[7])*mTot, 1./cimag(modefreqs->data[7])/mTot);
            XLAL_PRINT_INFO("matchrange %.16e,%.16e,%.16e\n", matchrange->data[0] * mTot / dt, matchrange->data[1] * mTot / dt, matchrange->data[2] * mTot / dt - 2);
    }
             /*XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, 1, nmodes, approximant);



             modefreqs->data[6] =  modefreqs_xtr->data[0];
             modefreqs->data[7] =  modefreqs_xtr->data[1];
             //modefreqs->data[6] =  conjl(-1.0 * modefreqs->data[0]);
             //modefreqs->data[7] =  conjl(-1.0 *  modefreqs->data[1]);
             //modefreqs->data[6] =  conjl(-1.0 * modefreqs_xtr->data[0]);
             //modefreqs->data[7] =  conjl(-1.0 *  modefreqs_xtr->data[1]);
             if (m <0) {
             //    modefreqs->data[6] =  modefreqs_xtr->data[0];
             //    modefreqs->data[7] =  modefreqs_xtr->data[1];
                 modefreqs->data[6] =  conjl(-1.0 * modefreqs_xtr->data[0]);
                 modefreqs->data[7] =  conjl(-1.0 *  modefreqs_xtr->data[1]);
             }
             XLALDestroyCOMPLEX16Vector(modefreqs_xtr);*/



             // }}}

		} //end if m= 2, l = 2
        if ((int)fabs((REAL8) m) == 1 && l == 2) {

            if (debugout)
                XLAL_PRINT_INFO("Stas, default (if nothing specified above) for pQNM eta = %f, chi = %f chi1 = %f \n", eta, chi,  chi1);
            if (debugout) {
                XLAL_PRINT_INFO("Stas: JLN*eta = %f, \n", JLN*eta);
                XLAL_PRINT_INFO("NRPeakOmega22 = %3.10f,  %3.10f,  %3.10f,  %3.10f\n", NRPeakOmega22 * mTot /finalMass, NRPeakOmega22 / finalMass,  finalMass,  mTot);
                for (j = 0; j < nmodes; j++) {
                    XLAL_PRINT_INFO("QNM frequencies: %d %d %d %3.10f %3.10f\n", l, m, j, creal(modefreqs->data[j]) * mTot, 1. / cimag(modefreqs->data[j]) / mTot);
                }
            }
            
            COMPLEX16Vector *modefreqs_xtr;
            modefreqs_xtr = XLALCreateCOMPLEX16Vector(nmodes);
            if (!modefreqs_xtr) {
                XLAL_ERROR(XLAL_ENOMEM);
            }
            
            if ( 1==1 || (JLN < 0.0 && eta <= 0.1)){
                spin1[0] *= -1;
                spin1[1] *= -1;
                spin1[2] *= -1;
                spin2[0] *= -1;
                spin2[1] *= -1;
                spin2[2] *= -1;
                XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l,  -mHere, nmodes, approximant);
                spin1[0] *= -1;
                spin1[1] *= -1;
                spin1[2] *= -1;
                spin2[0] *= -1;
                spin2[1] *= -1;
                spin2[2] *= -1;
                modefreqs->data[7] =  modefreqs_xtr->data[0];
                modefreqs->data[6] = NRPeakOmega22 + I/mTot/((1./cimag(modefreqs->data[0])/mTot)/2.);

                if (m == -1){
                    modefreqs->data[6] =  conjl(-1.0 * modefreqs->data[6]);
                    modefreqs->data[7] =  conjl(-1.0 * modefreqs->data[7]);
                }
            }
 //           modefreqs->data[7] = 0.5*(NRPeakOmega22 + modefreqs->data[0]) + I/mTot/((1./cimag(modefreqs->data[0])/mTot)/3.);
            XLALDestroyCOMPLEX16Vector(modefreqs_xtr);
	if (debugout) {
            XLAL_PRINT_INFO("l,m = %d %d\n",l,m);
            XLAL_PRINT_INFO("finalSpin = %f\n",finalSpin);
            XLAL_PRINT_INFO("finalMass = %f\n",finalMass);
            XLAL_PRINT_INFO("PeakOmega = %f\n",NRPeakOmega22*mTot);
            XLAL_PRINT_INFO("w0 = %f, t0 = %f\n",creal(modefreqs->data[0])*mTot, 1./cimag(modefreqs->data[0])/mTot);
            XLAL_PRINT_INFO("w1 = %f, t1 = %f\n",creal(modefreqs->data[1])*mTot, 1./cimag(modefreqs->data[1])/mTot);
            XLAL_PRINT_INFO("w2 = %f, t2 = %f\n",creal(modefreqs->data[2])*mTot, 1./cimag(modefreqs->data[2])/mTot);
            XLAL_PRINT_INFO("w3 = %f, t3 = %f\n",creal(modefreqs->data[3])*mTot, 1./cimag(modefreqs->data[3])/mTot);
            XLAL_PRINT_INFO("w4 = %f, t4 = %f\n",creal(modefreqs->data[4])*mTot, 1./cimag(modefreqs->data[4])/mTot);
            XLAL_PRINT_INFO("w5 = %f, t5 = %f\n",creal(modefreqs->data[5])*mTot, 1./cimag(modefreqs->data[5])/mTot);
            XLAL_PRINT_INFO("w6 = %f, t6 = %f\n",creal(modefreqs->data[6])*mTot, 1./cimag(modefreqs->data[6])/mTot);
            XLAL_PRINT_INFO("w7 = %f, t7 = %f\n",creal(modefreqs->data[7])*mTot, 1./cimag(modefreqs->data[7])/mTot);
            XLAL_PRINT_INFO("matchrange %.16e,%.16e,%.16e\n", matchrange->data[0] * mTot / dt, matchrange->data[1] * mTot / dt, matchrange->data[2] * mTot / dt - 2);
        }
        }
            // FIXME
            if (1==0 &&(m==1 || m==-1)){
                /*NRPeakOmega22 *= 0.5;
                if (JLN*m < 0.0){
                    modefreqs->data[5] = (-3. / 4. * NRPeakOmega22 / finalMass) + (1. / 4. * creal(modefreqs->data[0]));
                    modefreqs->data[4] = (-2. / 3. * NRPeakOmega22 / finalMass) + (1. / 3. * creal(modefreqs->data[0]));
                 } else {
                    modefreqs->data[5] = (3. / 4. * NRPeakOmega22 / finalMass) + (1. / 4. * creal(modefreqs->data[0]));
                    modefreqs->data[4] = (2. / 3. * NRPeakOmega22 / finalMass) + (1. / 3. * creal(modefreqs->data[0]));
                    //XLAL_PRINT_INFO("Stas , here m = %d modes: 0 = %4.10f, 6 = %4.10f, 7 = %4.10f \n", m, creal(modefreqs->data[0])*mTot, creal(modefreqs->data[6])*mTot, creal(modefreqs->data[7])*mTot);
                 }
                 modefreqs->data[5] += I * 3.5 / 0.9 * cimag(modefreqs->data[0]);
                 modefreqs->data[4] += I * 3.5 * cimag(modefreqs->data[0]);*/

                 COMPLEX16Vector *modefreqs_xtr;
                 modefreqs_xtr = XLALCreateCOMPLEX16Vector(nmodes);
                 if (!modefreqs_xtr) {
                     XLAL_ERROR(XLAL_ENOMEM);
                 }
                 //XLAL_PRINT_INFO("Stas, generating QNM freqs... ");
                 if (JLN > 0.0){
                     XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, -1, nmodes, approximant);
                     //XLAL_PRINT_INFO("done\n");
                     modefreqs->data[5] =  modefreqs_xtr->data[0];
                     //modefreqs->data[6] =  modefreqs_xtr->data[1];
                     //modefreqs->data[7] =  modefreqs_xtr->data[2];
                     if (JLN < kappa_thr || eta*JLN < eJL_thr){
                         modefreqs->data[6] =  modefreqs_xtr->data[1];
                         modefreqs->data[7] =  modefreqs_xtr->data[2];
                     }
                     if (m == -1){
                         //XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, -1, nmodes, approximant);
                         modefreqs->data[5] =  conjl(-1.0 * modefreqs_xtr->data[0]);
                         //modefreqs->data[6] =  conjl(-1.0 * modefreqs_xtr->data[1]);
                         //modefreqs->data[7] =  conjl(-1.0 * modefreqs_xtr->data[2]);
                         if (JLN < kappa_thr || eta*JLN < eJL_thr){
                             modefreqs->data[6] =  conjl(-1.0 *  modefreqs_xtr->data[1]);
                             modefreqs->data[7] =  conjl(-1.0 * modefreqs_xtr->data[2]);
                         }
                     }
                 }
                 //XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, -2, nmodes, approximant);
                 //if (JLN > 0.0){
                 //    modefreqs->data[6] =  modefreqs_xtr->data[0];
                     //modefreqs->data[5] =  modefreqs_xtr->data[1];
                //     if (m == -1){
                         //XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, 2, nmodes, approximant);
                //         modefreqs->data[6] =  conjl(-1.0 * modefreqs_xtr->data[0]);
                         //modefreqs->data[5] =  conjl(-1.0 *  modefreqs_xtr->data[1]);
                //     }
                // }
                 if (JLN < 0.0){
                     //modefreqs->data[4] =  conjl(-1.0 * modefreqs->data[4]);
                     //modefreqs->data[5] =  conjl(-1.0 *  modefreqs->data[5]);
                     XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, 1, nmodes, approximant);
                     //XLAL_PRINT_INFO("done\n");
                     modefreqs->data[5] =  modefreqs_xtr->data[0];
                     //modefreqs->data[6] =  modefreqs_xtr->data[1];
                     //modefreqs->data[7] =  modefreqs_xtr->data[2];
                     if (JLN > -1.0*kappa_thr || fabs(eta*JLN) < eJL_thr){
                         modefreqs->data[6] =  modefreqs_xtr->data[1];
                         modefreqs->data[7] =  modefreqs_xtr->data[2];
                     }
                     if (m == -1){
                         //XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, -1, nmodes, approximant);
                         modefreqs->data[5] =  conjl(-1.0 * modefreqs_xtr->data[0]);
                         //modefreqs->data[6] =  conjl(-1.0 * modefreqs_xtr->data[1]);
                         //modefreqs->data[7] =  conjl(-1.0 * modefreqs_xtr->data[2]);
                         if (JLN > -1.0*kappa_thr || fabs(eta*JLN) < eJL_thr){
                             modefreqs->data[6] =  conjl(-1.0 *  modefreqs_xtr->data[1]);
                             modefreqs->data[7] =  conjl(-1.0 * modefreqs_xtr->data[2]);
                         }
                     }
                     //XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, 2, nmodes, approximant);
                     //modefreqs->data[6] =  modefreqs_xtr->data[0];
                     //modefreqs->data[5] =  modefreqs_xtr->data[1];
                     //if (m == -1){
                         //XLALSimIMREOBGenerateQNMFreqV2Prec(modefreqs_xtr, mass1, mass2, spin1, spin2, l, 2, nmodes, approximant);
                     //    modefreqs->data[6] =  conjl(-1.0 * modefreqs_xtr->data[0]);
                         //modefreqs->data[5] =  conjl(-1.0 *  modefreqs_xtr->data[1]);
                     //}
                 }


                 XLALDestroyCOMPLEX16Vector(modefreqs_xtr);
            }
            } //v3

        // FIXME !!!!
        //for (j =0; j<nmodes; j++){
        //if (m > 0) {
		//	for (j = 0; j < nmodes; j++) {
		//		modefreqs->data[j] = conjl(-1.0 * modefreqs->data[j]);
		//	}
		//}


				// Move ringdown comb boundaries to sampling points to avoid numerical artifacts.
//				matchrange->data[0] -= fmod(matchrange->data[0], dt / mTot);
//		matchrange->data[1] -= fmod(matchrange->data[1], dt / mTot);
		if (debugout) {
			XLAL_PRINT_INFO("NRPeakOmega22 = %3.10f,  %3.10f,  %3.10f,  %3.10f\n", NRPeakOmega22 * mTot /finalMass, NRPeakOmega22 / finalMass,  finalMass,  mTot);
			for (j = 0; j < nmodes; j++) {
				XLAL_PRINT_INFO("QNM frequencies: %d %d %d %3.10f %3.10f\n", l, m, j, creal(modefreqs->data[j]) * mTot, 1. / cimag(modefreqs->data[j]) / mTot);
			}
		}
		/*
		 * Ringdown signal length: 10 times the decay time of the n=0
		 * mode
		 */
		Nrdwave = (INT4) (EOB_RD_EFOLDS / cimag(modefreqs->data[0]) / dt);
        //XLAL_PRINT_INFO("STas, EOB_RD_EFOLDS = %f, freq = %f \n", EOB_RD_EFOLDS, cimag(modefreqs->data[0])*dt);

		/*
		 * Check the value of attpos, to prevent memory access
		 * problems later
		 */
		if (matchrange->data[0] * mTot / dt < 5 || matchrange->data[1] * mTot / dt > matchrange->data[2] * mTot / dt - 2) {
			XLALPrintError("More inspiral points needed for ringdown matching.\n");
			XLAL_PRINT_INFO("%.16e,%.16e,%.16e\n", matchrange->data[0] * mTot / dt, matchrange->data[1] * mTot / dt, matchrange->data[2] * mTot / dt - 2);
			XLALDestroyCOMPLEX16Vector(modefreqs);
			XLAL_ERROR(XLAL_EFAILED);
		}
		/*
		 * STEP 2) Based on least-damped-mode decay time, allocate
		 * memory for rigndown waveform
		 */

		/*
		 * Create memory for the ring-down and full waveforms, and
		 * derivatives of inspirals
		 */

		rdwave1 = XLALCreateREAL8Vector(Nrdwave);
		rdwave2 = XLALCreateREAL8Vector(Nrdwave);
		rinspwave = XLALCreateREAL8Vector(6);
		dinspwave = XLALCreateREAL8Vector(6);
		ddinspwave = XLALCreateREAL8Vector(6);
		inspwaves1 = XLALCreateREAL8VectorSequence(3, 6);
		inspwaves2 = XLALCreateREAL8VectorSequence(3, 6);

		/* Check memory was allocated */
		if (!rdwave1 || !rdwave2 || !rinspwave || !dinspwave
		    || !ddinspwave || !inspwaves1 || !inspwaves2) {
			XLALDestroyCOMPLEX16Vector(modefreqs);
			if (rdwave1)
				XLALDestroyREAL8Vector(rdwave1);
			if (rdwave2)
				XLALDestroyREAL8Vector(rdwave2);
			if (rinspwave)
				XLALDestroyREAL8Vector(rinspwave);
			if (dinspwave)
				XLALDestroyREAL8Vector(dinspwave);
			if (ddinspwave)
				XLALDestroyREAL8Vector(ddinspwave);
			if (inspwaves1)
				XLALDestroyREAL8VectorSequence(inspwaves1);
			if (inspwaves2)
				XLALDestroyREAL8VectorSequence(inspwaves2);
			XLAL_ERROR(XLAL_ENOMEM);
		}
		memset(rdwave1->data, 0, rdwave1->length * sizeof(REAL8));
		memset(rdwave2->data, 0, rdwave2->length * sizeof(REAL8));

		/*
		 * STEP 3) Get values and derivatives of inspiral waveforms
		 * at matching comb points
		 */

		/* Generate derivatives of the last part of inspiral waves */
		/* Get derivatives of signal1 */
		if (XLALGenerateHybridWaveDerivatives(rinspwave, dinspwave, ddinspwave, timeVec, signal1,
			    matchrange, dt, mass1, mass2) == XLAL_FAILURE) {
			XLALDestroyCOMPLEX16Vector(modefreqs);
			XLALDestroyREAL8Vector(rdwave1);
			XLALDestroyREAL8Vector(rdwave2);
			XLALDestroyREAL8Vector(rinspwave);
			XLALDestroyREAL8Vector(dinspwave);
			XLALDestroyREAL8Vector(ddinspwave);
			XLALDestroyREAL8VectorSequence(inspwaves1);
			XLALDestroyREAL8VectorSequence(inspwaves2);
			XLAL_ERROR(XLAL_EFUNC);
		}
		for (j = 0; j < 6; j++) {
			inspwaves1->data[j] = rinspwave->data[j];
			inspwaves1->data[j + 6] = dinspwave->data[j];
			inspwaves1->data[j + 12] = ddinspwave->data[j];
		}

		/* Get derivatives of signal2 */
		if (XLALGenerateHybridWaveDerivatives(rinspwave, dinspwave, ddinspwave, timeVec, signal2,
			    matchrange, dt, mass1, mass2) == XLAL_FAILURE) {
			XLALDestroyCOMPLEX16Vector(modefreqs);
			XLALDestroyREAL8Vector(rdwave1);
			XLALDestroyREAL8Vector(rdwave2);
			XLALDestroyREAL8Vector(rinspwave);
			XLALDestroyREAL8Vector(dinspwave);
			XLALDestroyREAL8Vector(ddinspwave);
			XLALDestroyREAL8VectorSequence(inspwaves1);
			XLALDestroyREAL8VectorSequence(inspwaves2);
			XLAL_ERROR(XLAL_EFUNC);
		}
		for (j = 0; j < 6; j++) {
			inspwaves2->data[j] = rinspwave->data[j];
			inspwaves2->data[j + 6] = dinspwave->data[j];
			inspwaves2->data[j + 12] = ddinspwave->data[j];
		}


		/*
		 * STEP 4) Solve QNM coefficients and generate ringdown
		 * waveforms
		 */

		/* Generate ring-down waveforms */
		if (XLALSimIMREOBHybridRingdownWave(rdwave1, rdwave2, dt, mass1, mass2, inspwaves1, inspwaves2,
				   modefreqs, matchrange) == XLAL_FAILURE) {
			XLALDestroyCOMPLEX16Vector(modefreqs);
			XLALDestroyREAL8Vector(rdwave1);
			XLALDestroyREAL8Vector(rdwave2);
			XLALDestroyREAL8Vector(rinspwave);
			XLALDestroyREAL8Vector(dinspwave);
			XLALDestroyREAL8Vector(ddinspwave);
			XLALDestroyREAL8VectorSequence(inspwaves1);
			XLALDestroyREAL8VectorSequence(inspwaves2);
			XLAL_ERROR(XLAL_EFUNC);
		}

        //XLAL_PRINT_INFO("Stas we're down with RD solving, attaching it now....");
		/*
		 * STEP 5) Stitch inspiral and ringdown waveoforms
		 */

		/*
		 * Generate full waveforms, by stitching inspiral and
		 * ring-down waveforms
		 */
		//UINT4 attachIdx = matchrange->data[1] * mTot / dt;
		UINT4		attachIdx = round(matchrange->data[1] * mTot / dt);
		for (j = 1; j < Nrdwave; ++j) {
			signal1->data[j + attachIdx] = rdwave1->data[j];
			signal2->data[j + attachIdx] = rdwave2->data[j];
		}
        //XLAL_PRINT_INFO(" Stas, done 1 length = %d, Nrdwave = %d, attachIdx = %d \n", signal1->length, Nrdwave, attachIdx);

		memset(signal1->data + Nrdwave + attachIdx, 0, (signal1->length - Nrdwave - attachIdx) * sizeof(REAL8));
		memset(signal2->data + Nrdwave + attachIdx, 0, (signal2->length - Nrdwave - attachIdx) * sizeof(REAL8));

        //XLAL_PRINT_INFO(" Stas, done 2 \n");

		/* Free memory */
		XLALDestroyCOMPLEX16Vector(modefreqs);
		XLALDestroyREAL8Vector(rdwave1);
		XLALDestroyREAL8Vector(rdwave2);
		XLALDestroyREAL8Vector(rinspwave);
		XLALDestroyREAL8Vector(dinspwave);
		XLALDestroyREAL8Vector(ddinspwave);
		XLALDestroyREAL8VectorSequence(inspwaves1);
		XLALDestroyREAL8VectorSequence(inspwaves2);

        //XLAL_PRINT_INFO(" Stas, done 3 \n");
//    matchrange->data[0] += fmod(matchrange->data[0], dt / mTot);
//    matchrange->data[1] += fmod(matchrange->data[1], dt / mTot);

		return XLAL_SUCCESS;
	}
#endif				/* _LALSIMIMREOBHYBRIDRINGDOWNPREC_C */
