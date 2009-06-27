/*
 * Copyright (C) 2009 Reinhard Prix
 * Copyright (C) 2007 Chris Messenger
 * Copyright (C) 2006 John T. Whelan, Badri Krishnan
 * Copyright (C) 2005, 2006, 2007 Reinhard Prix
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

/** \author R. Prix, J. T. Whelan
 * \ingroup pulsarCoherent
 * \file
 * \brief
 * Functions to calculate the so-called F-statistic for a given point in parameter-space,
 * following the equations in \ref JKS98.
 *
 * This code is partly a descendant of an earlier implementation found in
 * LALDemod.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens, Bruce Allen
 * ComputSky.[ch] by Jolien Creighton, Reinhard Prix, Steve Berukoff
 * LALComputeAM.[ch] by Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens
 *
 * NOTE: this file contains specialized versions of the central CFSv2 routines, that are aimed at GPU optimization.
 * At the moment this only means they are internally using only single precision, but still aggree to within 1%
 * for Tobs ~ 1day and fmax ~ 1kHz.
 *
 */

/*---------- INCLUDES ----------*/
#define __USE_ISOC99 1
#include <math.h>

#include <lal/AVFactories.h>
#include <lal/ComputeFstat.h>

#include "ComputeFstatGPU.h"

NRCSID( COMPUTEFSTATC, "$Id$");

/*---------- local DEFINES ----------*/
#define LD_SMALL4       ((REAL4)2.0e-4)		/**< "small" number for REAL4*/

#define TWOPI_FLOAT     6.28318530717958f  	/**< single-precision 2*pi */
#define OOTWOPI_FLOAT   (1.0f / TWOPI_FLOAT)	/**< single-precision 1 / (2pi) */

/*----- Macros ----- */
#define SQ(x) ( (x) * (x) )
#define REM(x) ( (x) - (INT4)(x) )

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
static const REAL4 inv_fact[PULSAR_MAX_SPINS] = { 1.0, 1.0, (1.0/2.0), (1.0/6.0), (1.0/24.0), (1.0/120.0), (1.0/720.0) };

/* global sin-cos lookup table */
#define LUT_RES         	64      /* resolution of lookup-table */
#define OO_LUT_RES		(1.0f / LUT_RES )

static REAL4 sinVal[LUT_RES+1], cosVal[LUT_RES+1];

/* empty initializers  */
static const LALStatus empty_LALStatus;
static const AMCoeffs empty_AMCoeffs;

const SSBtimes empty_SSBtimes;
const MultiSSBtimes empty_MultiSSBtimes;
const AntennaPatternMatrix empty_AntennaPatternMatrix;
const MultiAMCoeffs empty_MultiAMCoeffs;
const Fcomponents empty_Fcomponents;
const ComputeFBuffer empty_ComputeFBuffer;
const PulsarSpins_REAL4 empty_PulsarSpins_REAL4;
const ComputeFBuffer_REAL4 empty_ComputeFBuffer_REAL4;
const Fcomponents_REAL4 empty_Fcomponents_REAL4;

/*---------- internal prototypes ----------*/
int finite(double x);
void sin_cos_2PI_LUT_REAL4 (REAL4 *sin2pix, REAL4 *cos2pix, REAL4 x);
void sin_cos_LUT_REAL4 (REAL4 *sinx, REAL4 *cosx, REAL4 x);
void init_sin_cos_LUT_REAL4 (void);

/*==================== FUNCTION DEFINITIONS ====================*/

/** REAL4 and GPU-ready version of ComputeFStatFreqBand().
    Function to compute a vector of Fstatistic values for a number of frequency bins.
    This function is simply a wrapper for ComputeFstat() which is called repeatedly for
    every frequency value.  The output, i.e. fstatVector must be properly allocated
    before this function is called.  The values of the start frequency, the step size
    in the frequency and the number of frequency values for which the Fstatistic is
    to be calculated are read from fstatVector.  The other parameters are not checked and
    they must be correctly set outside this function.
*/
int
XLALComputeFStatFreqBand (   REAL4FrequencySeries *fstatVector, 	/**< [out] Vector of Fstat values */
                             const PulsarDopplerParams *doppler,	/**< parameter-space point to compute F for */
                             const MultiSFTVector *multiSFTs, 		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
                             const MultiNoiseWeights *multiWeights,	/**< noise-weights of all SFTs */
                             const MultiDetectorStateSeries *multiDetStates,/**< 'trajectories' of the different IFOs */
                             const ComputeFParams *params		/**< addition computational params */
                             )
{
  static const char *fn = "XLALComputeFStatFreqBand()";

  UINT4 numIFOs, numBins, k;
  REAL8 deltaF;
  REAL4 Fstat;
  PulsarDopplerParams thisPoint;
  ComputeFBuffer_REAL4 cfBuffer = empty_ComputeFBuffer_REAL4;

  /* check input consistency */
  if ( !doppler || !multiSFTs || !multiWeights || !multiDetStates || !params ) {
    XLALPrintError ("%s: illegal NULL input pointer.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  numIFOs = multiSFTs->length;
  if ( multiDetStates->length != numIFOs ) {
    XLALPrintError ("%s: number of IFOs inconsisten between multiSFTs (%d) and multiDetStates (%d).\n", fn, numIFOs, multiDetStates->length );
    XLAL_ERROR (fn, XLAL_EINVAL );
  }

  if ( !fstatVector || !fstatVector->data || !fstatVector->data->data || !fstatVector->data->length ) {
    XLALPrintError ("%s: illegal NULL or empty output pointer 'fstatVector'.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* a check that the f0 values from thisPoint and fstatVector are
     at least close to each other -- this is only meant to catch
     stupid errors but not subtle ones */
  if ( fabs(fstatVector->f0 - doppler->fkdot[0]) >= fstatVector->deltaF ) {
    XLALPrintError ("%s: fstatVector->f0 = %f differs from doppler->fkdot[0] = %f by more than deltaF = %f\n",
                    fn, fstatVector->f0, doppler->fkdot[0], fstatVector->deltaF );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* copy values from 'doppler' to local variable 'thisPoint' */
  thisPoint = *doppler;

  numBins = fstatVector->data->length;
  deltaF = fstatVector->deltaF;

  /* loop over frequency values and fill up values in fstatVector */
  for ( k = 0; k < numBins; k++)
    {
      if ( XLALDriverFstatGPU ( &Fstat, &thisPoint, multiSFTs, multiWeights, multiDetStates, params->Dterms, &cfBuffer ) != XLAL_SUCCESS ) {
        XLALPrintError ("%s: XLALDriverFstatGPU() failed with errno = %d.\n", fn, xlalErrno );
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      fstatVector->data->data[k] = Fstat;

      thisPoint.fkdot[0] += deltaF;
    }

  XLALEmptyComputeFBuffer_REAL4 ( &cfBuffer );

  return XLAL_SUCCESS;

} /* XLALComputeFStatFreqBand() */





/** Host-bound 'driver' function for the central F-stat computation
 *  of a single F-stat value for one parameter-space point.
 *
 * This is a GPU-adapted replacement for ComputeFstat(), and implements a 'wrapper'
 * around the core F-stat routines that can be executed as kernels on a GPU device.
 *
 */
int
XLALDriverFstatGPU ( REAL4 *Fstat,	                 		/**< [out] Fstatistic value 'F' */
                     const PulsarDopplerParams *doppler, 		/**< parameter-space point to compute F for */
                     const MultiSFTVector *multiSFTs,    		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
                     const MultiNoiseWeights *multiWeights,		/**< noise-weights of all SFTs */
                     const MultiDetectorStateSeries *multiDetStates,	/**< 'trajectories' of the different IFOs */
                     UINT4 Dterms,       				/**< number of terms to use in Dirichlet kernel interpolation */
                     ComputeFBuffer_REAL4 *cfBuffer       		/**< CF-internal buffering structure */
                     )
{
  static const char *fn = "XLALDriverFstatGPU()";

  MultiSSBtimes_REAL4 *multiSSB4 = NULL;
  MultiAMCoeffs *multiAMcoef = NULL;

  /* check input consistency */
  if ( !Fstat || !doppler || !multiSFTs || !multiDetStates || !cfBuffer ) {
    XLALPrintError("%s: Illegal NULL pointer input.\n", fn );
    XLAL_ERROR (fn, XLAL_EINVAL );
  }
  UINT4 numIFOs = multiSFTs->length;
  if ( multiDetStates->length != numIFOs ) {
    XLALPrintError("%s: inconsistent number of IFOs in SFTs (%d) and detector-states (%d).\n", fn, numIFOs, multiDetStates->length);
    XLAL_ERROR (fn, XLAL_EINVAL );
  }
  if ( multiWeights && multiWeights->length != numIFOs ) {
    XLALPrintError("%s: inconsistent number of IFOs in SFTs (%d) and noise-weights (%d).\n", fn, numIFOs, multiWeights->length);
    XLAL_ERROR (fn, XLAL_EINVAL );
  }

  /* make sure sin/cos lookup-tables are initialized */
  init_sin_cos_LUT_REAL4();

  /* check if for this skyposition and data, the SSB+AMcoef were already buffered */
  if ( cfBuffer->Alpha == doppler->Alpha && cfBuffer->Delta == doppler->Delta && multiDetStates == cfBuffer->multiDetStates )
    { /* yes ==> reuse */
      multiSSB4   = cfBuffer->multiSSB;
      multiAMcoef = cfBuffer->multiAMcoef;
    }
  else
    {
      /* compute new SSB timings */
      if ( (multiSSB4 = XLALGetMultiSSBtimes_REAL4 ( multiDetStates, doppler->Alpha, doppler->Delta, doppler->refTime)) == NULL ) {
        XLALPrintError ( "%s: XLALGetMultiSSBtimes_REAL4() failed. xlalErrno = %d.\n", fn, xlalErrno );
	XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      /* compute new AM-coefficients */
      SkyPosition skypos;
      LALStatus status = empty_LALStatus;
      skypos.system =   COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = doppler->Alpha;
      skypos.latitude  = doppler->Delta;

      LALGetMultiAMCoeffs ( &status, &multiAMcoef, multiDetStates, skypos );
      if ( status.statusCode ) {
        XLALPrintError ("%s: LALGetMultiAMCoeffs() failed with statusCode=%d, '%s'\n", fn, status.statusCode, status.statusDescription );
        XLAL_ERROR ( fn, XLAL_EFAILED );
      }

      /* apply noise-weights to Antenna-patterns and compute A,B,C */
      if ( XLALWeighMultiAMCoeffs ( multiAMcoef, multiWeights ) != XLAL_SUCCESS ) {
	XLALPrintError("%s: XLALWeighMultiAMCoeffs() failed with error = %d\n", fn, xlalErrno );
	XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      /* store these in buffer */
      cfBuffer->Alpha = doppler->Alpha;
      cfBuffer->Delta = doppler->Delta;
      cfBuffer->multiDetStates = multiDetStates ;

      XLALDestroyMultiAMCoeffs ( cfBuffer->multiAMcoef );
      XLALDestroyMultiSSBtimes_REAL4 ( cfBuffer->multiSSB );

      cfBuffer->multiSSB = multiSSB4;
      cfBuffer->multiAMcoef = multiAMcoef;

    } /* if we could NOT reuse previously buffered quantites */


  /* ---------- prepare REAL4 version of PulsarSpins to be passed into core functions */
  PulsarSpins_REAL4 fkdot4 = empty_PulsarSpins_REAL4;
  UINT4 s;
  fkdot4.FreqMain = (INT4)doppler->fkdot[0];
  /* fkdot[0] now only carries the *remainder* of Freq wrt FreqMain */
  fkdot4.fkdot[0] = (REAL4)( doppler->fkdot[0] - (REAL8)fkdot4.FreqMain );
  /* the remaining spins are simply down-cast to REAL4 */
  for ( s=1; s < PULSAR_MAX_SPINS; s ++ )
    fkdot4.fkdot[s] = (REAL4) doppler->fkdot[s];


  /* call the core function to compute one multi-IFO F-statistic */

  /* >>>>> this function could run on the GPU device <<<<< */
  XLALCoreFstatGPU (Fstat, &fkdot4, multiSFTs, multiSSB4, multiAMcoef, Dterms );

  if ( xlalErrno ) {
    XLALPrintError ("%s: XLALCoreFstatGPU() failed with errno = %d.\n", xlalErrno );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;


} /* XLALDriverFstatGPU() */


/** This function computes a multi-IFO F-statistic value for given frequency (+fkdots),
 * antenna-pattern functions, SSB-timings and data ("SFTs").
 *
 * It is *only* using single-precision quantities. The aim is that this function can easily
 * be turned into code to be run on a GPU device.
 *
 */
void
XLALCoreFstatGPU (REAL4 *Fstat,				/**< [out] multi-IFO F-statistic value 'F' */
                  PulsarSpins_REAL4 *fkdot4,		/**< [in] frequency and derivatives in REAL4-safe format */
                  const MultiSFTVector *multiSFTs,	/**< [in] input multi-IFO data ("SFTs") */
                  MultiSSBtimes_REAL4 *multiSSB4,	/**< [in] multi-IFO SSB-timings in REAL4-safe format */
                  MultiAMCoeffs *multiAMcoef,		/**< [in] multi-IFO antenna-pattern functions */
                  UINT4 Dterms				/**< [in] parameter: number of Dterms to use in Dirichlet kernel */
                  )
{
  static const char *fn = "XLALCoreFstatGPU()";

  UINT4 numIFOs, X;
  REAL4 Fa_re, Fa_im, Fb_re, Fb_im;

#ifndef LAL_NDEBUG
  /* check input consistency */
  if ( !Fstat || !fkdot4 || !multiSFTs || !multiSSB4 || !multiAMcoef ) {
    XLALPrintError ("%s: illegal NULL input.\n", fn);
    XLAL_ERROR_VOID ( fn, XLAL_EINVAL );
  }
  if ( multiSFTs->length == 0 || multiSSB4->length == 0 || multiAMcoef->length == 0 ||
       !multiSFTs->data || !multiSSB4->data || !multiAMcoef->data )
    {
      XLALPrintError ("%s: invalid empty input.\n", fn);
      XLAL_ERROR_VOID ( fn, XLAL_EINVAL );
    }
#endif

  numIFOs = multiSFTs->length;

#ifndef LAL_NDEBUG
  if ( multiSSB4->length != numIFOs || multiAMcoef->length != numIFOs ) {
    XLALPrintError ("%s: inconsistent number of IFOs between multiSFTs, multiSSB4 and multiAMcoef.\n", fn);
    XLAL_ERROR_VOID ( fn, XLAL_EINVAL );
  }
#endif

  /* ----- loop over detectors and compute all detector-specific quantities ----- */
  Fa_re = Fa_im = Fb_re = Fb_im = 0;

  for ( X=0; X < numIFOs; X ++)
    {
      Fcomponents_REAL4 FcX;

      XLALComputeFaFb_REAL4 (&FcX, multiSFTs->data[X], fkdot4, multiSSB4->data[X], multiAMcoef->data[X], Dterms);

#ifndef LAL_NDEBUG
      if ( xlalErrno ) {
        XLALPrintError ("%s: XALComputeFaFb_REAL4() failed\n", fn );
        XLAL_ERROR_VOID ( fn, XLAL_EFUNC );
      }
      if ( !finite(FcX.Fa.re) || !finite(FcX.Fa.im) || !finite(FcX.Fb.re) || !finite(FcX.Fb.im) ) {
	XLALPrintError("%s: XLALComputeFaFb_REAL4() returned non-finite: Fa_X=(%f,%f), Fb_X=(%f,%f) for X=%d\n",
                       fn, FcX.Fa.re, FcX.Fa.im, FcX.Fb.re, FcX.Fb.im, X );
	XLAL_ERROR_VOID ( fn, XLAL_EFPINVAL );
      }
#endif

      /* Fa = sum_X Fa_X */
      Fa_re += FcX.Fa.re;
      Fa_im += FcX.Fa.im;

      /* Fb = sum_X Fb_X */
      Fb_re += FcX.Fb.re;
      Fb_im += FcX.Fb.im;

    } /* for  X < numIFOs */

  /* ----- compute final Fstatistic-value ----- */

  /* convenient shortcuts */
  REAL4 Ad = multiAMcoef->Mmunu.Ad;
  REAL4 Bd = multiAMcoef->Mmunu.Bd;
  REAL4 Cd = multiAMcoef->Mmunu.Cd;
  REAL4 Dd_inv = 1.0f / multiAMcoef->Mmunu.Dd;

  /* NOTE: the data MUST be normalized by the DOUBLE-SIDED PSD (using LALNormalizeMultiSFTVect),
   * therefore there is a factor of 2 difference with respect to the equations in JKS, which
   * where based on the single-sided PSD.
   */
  (*Fstat) = Dd_inv * ( Bd * (SQ(Fa_re) + SQ(Fa_im) )
                      + Ad * ( SQ(Fb_re) + SQ(Fb_im) )
                      - 2.0f * Cd *( Fa_re * Fb_re + Fa_im * Fb_im )
                      );

  return;

} /* XLALCoreFstatGPU() */



/** Revamped version of LALDemod() (based on TestLALDemod() in CFS).
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 *
 * Note: this is a single-precision version aimed for GPU parallelization.
 *
 */
void
XLALComputeFaFb_REAL4 ( Fcomponents_REAL4 *FaFb,		/**< [out] single-IFO Fa/Fb for this parameter-space point */
                        const SFTVector *sfts,			/**< [in] single-IFO input data ("SFTs") */
                        const PulsarSpins_REAL4 *fkdot4,	/**< [in] frequency (and derivatives) in REAL4-safe format */
                        const SSBtimes_REAL4 *tSSB,		/**< [in] single-IFO SSB-timings in REAL4-safe format */
                        const AMCoeffs *amcoe,			/**< [in] single-IFO antenna-pattern coefficients */
                        UINT4 Dterms)				/**< [in] number of Dterms to use in Dirichlet kernel */
{
  static const char *fn = "XLALComputeFaFb_REAL4()";

  UINT4 alpha;                 	/* loop index over SFTs */
  UINT4 spdnOrder;		/* maximal spindown-orders */
  UINT4 numSFTs;		/* number of SFTs (M in the Notes) */
  COMPLEX8 Fa, Fb;
  REAL4 Tsft; 			/* length of SFTs in seconds */
  INT4 freqIndex0;		/* index of first frequency-bin in SFTs */
  INT4 freqIndex1;		/* index of last frequency-bin in SFTs */

  REAL4 *a_al, *b_al;		/* pointer to alpha-arrays over a and b */

  REAL4 *DeltaT_int_al, *DeltaT_rem_al, *TdotM1_al;	/* pointer to alpha-arrays of SSB-timings */
  SFTtype *SFT_al;		/* SFT alpha  */

  REAL4 norm = OOTWOPI_FLOAT;

  /* ----- check validity of input */
#ifndef LAL_NDEBUG
  if ( !FaFb ) {
    XLALPrintError ("%s: Output-pointer is NULL !\n", fn);
    XLAL_ERROR_VOID ( fn, XLAL_EINVAL);
  }

  if ( !sfts || !sfts->data ) {
    XLALPrintError ("%s: Input SFTs are NULL!\n", fn);
    XLAL_ERROR_VOID ( fn, XLAL_EINVAL);
  }

  if ( !fkdot4 || !tSSB || !tSSB->DeltaT_int || !tSSB->DeltaT_rem || !tSSB->TdotM1 || !amcoe || !amcoe->a || !amcoe->b )
    {
      XLALPrintError ("%s: Illegal NULL in input !\n", fn);
      XLAL_ERROR_VOID ( fn, XLAL_EINVAL);
    }
#endif

  /* ----- prepare convenience variables */
  numSFTs = sfts->length;
  Tsft = (REAL4)( 1.0 / sfts->data[0].deltaF );

  REAL4 dFreq = sfts->data[0].deltaF;
  freqIndex0 = (UINT4) ( sfts->data[0].f0 / dFreq + 0.5); /* lowest freqency-index */
  freqIndex1 = freqIndex0 + sfts->data[0].data->length;

  REAL4 f0 = fkdot4->FreqMain;
  REAL4 df = fkdot4->fkdot[0];

  /* find highest non-zero spindown-entry */
  for ( spdnOrder = PULSAR_MAX_SPINS - 1;  spdnOrder > 0 ; spdnOrder --  )
    if ( fkdot4->fkdot[spdnOrder] )
      break;

  Fa.re = 0.0f;
  Fa.im = 0.0f;
  Fb.re = 0.0f;
  Fb.im = 0.0f;

  /* convenient shortcuts, pointers to beginning of alpha-arrays */
  a_al = amcoe->a->data;
  b_al = amcoe->b->data;
  DeltaT_int_al = tSSB->DeltaT_int->data;
  DeltaT_rem_al = tSSB->DeltaT_rem->data;
  TdotM1_al = tSSB->TdotM1->data;

  SFT_al = sfts->data;

  /* Loop over all SFTs  */
  for ( alpha = 0; alpha < numSFTs; alpha++ )
    {
      REAL4 a_alpha, b_alpha;

      INT4 kstar;		/* central frequency-bin k* = round(xhat_alpha) */
      INT4 k0, k1;

      COMPLEX8 *Xalpha = SFT_al->data->data; /* pointer to current SFT-data */
      COMPLEX8 *Xalpha_l; 	/* pointer to frequency-bin k in current SFT */
      REAL4 s_alpha, c_alpha;	/* sin(2pi kappa_alpha) and (cos(2pi kappa_alpha)-1) */
      REAL4 realQ, imagQ;	/* Re and Im of Q = e^{-i 2 pi lambda_alpha} */
      REAL4 realXP, imagXP;	/* Re/Im of sum_k X_ak * P_ak */
      REAL4 realQXP, imagQXP;	/* Re/Im of Q_alpha R_alpha */

      REAL4 lambda_alpha;
      REAL4 kappa_max, kappa_star;

      /* ----- calculate kappa_max and lambda_alpha */
      {
	UINT4 s; 		/* loop-index over spindown-order */

	REAL4 Tas;		/* temporary variable to calculate (DeltaT_alpha)^s */
        REAL4 T0, dT, deltaT;
        REAL4 Dphi_alpha_int, Dphi_alpha_rem;
        REAL4 phi_alpha_rem;
        REAL4 TdotM1;		/* defined as Tdot_al - 1 */

        /* 1st oder: s = 0 */
        TdotM1 = *TdotM1_al;
        T0 = *DeltaT_int_al;
        dT = *DeltaT_rem_al;
        deltaT = T0 + dT;

        /* phi_alpha = f * Tas; */
        phi_alpha_rem = f0 * dT + df * T0 + df * dT;
        phi_alpha_rem = REM ( phi_alpha_rem );
        Dphi_alpha_int = f0;
        Dphi_alpha_rem = df * (1.0f + TdotM1) + f0 * TdotM1;

        /* higher-order spindowns */
        Tas = deltaT;
	for (s=1; s <= spdnOrder; s++)
	  {
	    REAL4 fsdot = fkdot4->fkdot[s];
	    Dphi_alpha_rem += fsdot * Tas * inv_fact[s]; 	/* here: Tas = DT^s */
	    Tas *= deltaT;					/* now:  Tas = DT^(s+1) */
	    phi_alpha_rem += fsdot * Tas * inv_fact[s+1];
	  } /* for s <= spdnOrder */

	/* Step 3: apply global factor of Tsft to complete Dphi_alpha */
        Dphi_alpha_int *= Tsft;
	Dphi_alpha_rem *= Tsft;

        REAL4 tmp = REM( 0.5f * Dphi_alpha_int ) + REM ( 0.5f * Dphi_alpha_rem );
	lambda_alpha = phi_alpha_rem - tmp;

	/* real- and imaginary part of e^{-i 2 pi lambda_alpha } */
	sin_cos_2PI_LUT_REAL4 ( &imagQ, &realQ, - lambda_alpha );

        kstar = (INT4)Dphi_alpha_int + (INT4)Dphi_alpha_rem;
	kappa_star = REM(Dphi_alpha_int) + REM(Dphi_alpha_rem);
	kappa_max = kappa_star + 1.0f * Dterms - 1.0f;

	/* ----- check that required frequency-bins are found in the SFTs ----- */
	k0 = kstar - Dterms + 1;
	k1 = k0 + 2 * Dterms - 1;
	if ( (k0 < freqIndex0) || (k1 > freqIndex1) )
	  {
	    XLALPrintError ("%s: Required frequency-bins [%d, %d] not covered by SFT-interval [%d, %d]\n\n",
                            fn, k0, k1, freqIndex0, freqIndex1 );
	    XLAL_ERROR_VOID( fn, XLAL_EDOM);
	  }

      } /* compute kappa_star, lambda_alpha */

      /* NOTE: sin[ 2pi (Dphi_alpha - k) ] = sin [ 2pi Dphi_alpha ], therefore
       * the trig-functions need to be calculated only once!
       * We choose the value sin[ 2pi(Dphi_alpha - kstar) ] because it is the
       * closest to zero and will pose no numerical difficulties !
       */
      sin_cos_2PI_LUT_REAL4 ( &s_alpha, &c_alpha, kappa_star );
      c_alpha -= 1.0f;

      /* ---------- calculate the (truncated to Dterms) sum over k ---------- */

      /* ---------- ATTENTION: this the "hot-loop", which will be
       * executed many millions of times, so anything in here
       * has a HUGE impact on the whole performance of the code.
       *
       * DON'T touch *anything* in here unless you really know
       * what you're doing !!
       *------------------------------------------------------------
       */

      Xalpha_l = Xalpha + k0 - freqIndex0;  /* first frequency-bin in sum */

      realXP = 0;
      imagXP = 0;

      /* if no danger of denominator -> 0 */
      if ( ( kappa_star > LD_SMALL4 ) && (kappa_star < 1.0 - LD_SMALL4) )
	{
	  /* improved hotloop algorithm by Fekete Akos:
	   * take out repeated divisions into a single common denominator,
	   * plus use extra cleverness to compute the nominator efficiently...
	   */
	  REAL4 Sn = (*Xalpha_l).re;
	  REAL4 Tn = (*Xalpha_l).im;
	  REAL4 pn = kappa_max;
	  REAL4 qn = pn;
	  REAL4 U_alpha, V_alpha;

	  /* recursion with 2*Dterms steps */
	  UINT4 l;
	  for ( l = 1; l < 2*Dterms; l ++ )
	    {
	      Xalpha_l ++;

	      pn = pn - 1.0f; 			/* p_(n+1) */
	      Sn = pn * Sn + qn * (*Xalpha_l).re;	/* S_(n+1) */
	      Tn = pn * Tn + qn * (*Xalpha_l).im;	/* T_(n+1) */
	      qn *= pn;				/* q_(n+1) */
	    } /* for l <= 2*Dterms */

	  U_alpha = Sn / qn;
	  V_alpha = Tn / qn;

#ifndef LAL_NDEBUG
	  if ( !finite(U_alpha) || !finite(V_alpha) || !finite(pn) || !finite(qn) || !finite(Sn) || !finite(Tn) ) {
	    XLAL_ERROR_VOID (fn, XLAL_EFPINVAL );
	  }
#endif

	  realXP = s_alpha * U_alpha - c_alpha * V_alpha;
	  imagXP = c_alpha * U_alpha + s_alpha * V_alpha;

	} /* if |remainder| > LD_SMALL4 */
      else
	{ /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
	  UINT4 ind0;
  	  if ( kappa_star <= LD_SMALL4 ) ind0 = Dterms - 1;
  	  else ind0 = Dterms;
	  realXP = TWOPI_FLOAT * Xalpha_l[ind0].re;
	  imagXP = TWOPI_FLOAT * Xalpha_l[ind0].im;
	} /* if |remainder| <= LD_SMALL4 */

      realQXP = realQ * realXP - imagQ * imagXP;
      imagQXP = realQ * imagXP + imagQ * realXP;

      /* we're done: ==> combine these into Fa and Fb */
      a_alpha = (*a_al);
      b_alpha = (*b_al);

      Fa.re += a_alpha * realQXP;
      Fa.im += a_alpha * imagQXP;

      Fb.re += b_alpha * realQXP;
      Fb.im += b_alpha * imagQXP;


      /* advance pointers over alpha */
      a_al ++;
      b_al ++;

      DeltaT_int_al ++;
      DeltaT_rem_al ++;
      TdotM1_al ++;

      SFT_al ++;

    } /* for alpha < numSFTs */

  /* return result */
  FaFb->Fa.re = norm * Fa.re;
  FaFb->Fa.im = norm * Fa.im;
  FaFb->Fb.re = norm * Fb.re;
  FaFb->Fb.im = norm * Fb.im;

  return;

} /* XLALComputeFaFb_REAL4() */



/** Destroy a MultiSSBtimes_REAL4 structure.
 *  Note, this is "NULL-robust" in the sense that it will not crash
 *  on NULL-entries anywhere in this struct, so it can be used
 *  for failure-cleanup even on incomplete structs
 */
void
XLALDestroyMultiSSBtimes_REAL4 ( MultiSSBtimes_REAL4 *multiSSB )
{
  UINT4 X;

  if ( ! multiSSB )
    return;

  if ( multiSSB->data )
    {
      for ( X=0; X < multiSSB->length; X ++ )
	{
          XLALDestroySSBtimes_REAL4 ( multiSSB->data[X] );
	} /* for X < numDetectors */
      LALFree ( multiSSB->data );
    }

  LALFree ( multiSSB );

  return;

} /* XLALDestroyMultiSSBtimes_REAL4() */


/** Destroy a SSBtimes_REAL4 structure.
 *  Note, this is "NULL-robust" in the sense that it will not crash
 *  on NULL-entries anywhere in this struct, so it can be used
 *  for failure-cleanup even on incomplete structs
 */
void
XLALDestroySSBtimes_REAL4 ( SSBtimes_REAL4 *tSSB )
{
  if ( ! tSSB )
    return;

  if ( tSSB->DeltaT_int )
    XLALDestroyREAL4Vector ( tSSB->DeltaT_int );

  if ( tSSB->DeltaT_rem )
    XLALDestroyREAL4Vector ( tSSB->DeltaT_rem );

  if ( tSSB->TdotM1 )
    XLALDestroyREAL4Vector ( tSSB->TdotM1 );

  LALFree ( tSSB );

  return;

} /* XLALDestroySSBtimes_REAL4() */


/** Multi-IFO version of LALGetSSBtimes_REAL4().
 * Get all SSB-timings for all input detector-series in REAL4 representation.
 *
 */
MultiSSBtimes_REAL4 *
XLALGetMultiSSBtimes_REAL4 ( const MultiDetectorStateSeries *multiDetStates, 	/**< [in] detector-states at timestamps t_i */
                             REAL8 Alpha, REAL8 Delta,				/**< source sky-location, in equatorial coordinates  */
                             LIGOTimeGPS refTime
                             )
{
  static const char *fn = "XLALGetMultiSSBtimes_REAL4()";

  UINT4 X, numDetectors;
  MultiSSBtimes_REAL4 *ret;

  /* check input */
  if ( !multiDetStates || multiDetStates->length==0 ) {
    XLALPrintError ("%s: illegal NULL or empty input 'multiDetStates'.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  numDetectors = multiDetStates->length;

  if ( ( ret = XLALCalloc( 1, sizeof( *ret ) )) == NULL ) {
    XLALPrintError ("%s: XLALCalloc( 1, %d ) failed.\n", fn, sizeof( *ret ) );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  ret->length = numDetectors;
  if ( ( ret->data = XLALCalloc ( numDetectors, sizeof ( *ret->data ) )) == NULL ) {
    XLALFree ( ret );
    XLALPrintError ("%s: XLALCalloc( %d, %d ) failed.\n", fn, numDetectors, sizeof( *ret->data ) );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  for ( X=0; X < numDetectors; X ++ )
    {
      if ( (ret->data[X] = XLALGetSSBtimes_REAL4 (multiDetStates->data[X], Alpha, Delta, refTime)) == NULL ) {
        XLALPrintError ("%s: XLALGetSSBtimes_REAL4() failed. xlalErrno = %d\n", fn, xlalErrno );
        goto failed;
      }

    } /* for X < numDet */

  goto success;

 failed:
  /* free all memory allocated so far */
  XLALDestroyMultiSSBtimes_REAL4 ( ret );
  XLAL_ERROR_NULL (fn, XLAL_EFAILED );

 success:
  return ret;

} /* LALGetMultiSSBtimes() */



/** XLAL REAL4-version of LALGetSSBtimes()
 *
 */
SSBtimes_REAL4 *
XLALGetSSBtimes_REAL4 ( const DetectorStateSeries *DetectorStates,	/**< [in] detector-states at timestamps t_i */
                        REAL8 Alpha, REAL8 Delta,			/**< source sky-location, in equatorial coordinates  */
                        LIGOTimeGPS refTime
                        )
{
  static const char *fn = "XLALGetSSBtimes_REAL4()";

  UINT4 numSteps, i;
  REAL8 refTimeREAL8;
  SSBtimes_REAL4 *ret;

  /* check input consistency */
  if ( !DetectorStates || DetectorStates->length==0 ) {
    XLALPrintError ("%s: illegal NULL or empty input 'DetectorStates'.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  numSteps = DetectorStates->length;		/* number of timestamps */

  /* convenience variables */
  refTimeREAL8 = XLALGPSGetREAL8( &refTime );

  /* allocate return container */
  if ( ( ret = XLALMalloc( sizeof(*ret))) == NULL ) {
    XLALPrintError ("%s: XLALMalloc(%d) failed.\n", fn, sizeof(*ret) );
    goto failed;
  }
  if ( (ret->DeltaT_int  = XLALCreateREAL4Vector ( numSteps )) == NULL ||
       (ret->DeltaT_rem  = XLALCreateREAL4Vector ( numSteps )) == NULL ||
       (ret->TdotM1      = XLALCreateREAL4Vector ( numSteps )) == NULL )
    {
      XLALPrintError ("%s: XLALCreateREAL4Vector ( %d ) failed.\n", fn, numSteps );
      goto failed;
    }

  /* loop over all SFTs and compute timing info */
  for (i=0; i < numSteps; i++ )
    {
      BarycenterInput baryinput = empty_BarycenterInput;
      EmissionTime emit;
      DetectorState *state = &(DetectorStates->data[i]);
      LALStatus status = empty_LALStatus;
      REAL8 deltaT;
      REAL4 deltaT_int;

      baryinput.tgps = state->tGPS;
      baryinput.site = DetectorStates->detector;
      baryinput.site.location[0] /= LAL_C_SI;
      baryinput.site.location[1] /= LAL_C_SI;
      baryinput.site.location[2] /= LAL_C_SI;

      baryinput.alpha = Alpha;
      baryinput.delta = Delta;
      baryinput.dInv = 0;

      LALBarycenter(&status, &emit, &baryinput, &(state->earthState) );
      if ( status.statusCode ) {
        XLALPrintError ("%s: LALBarycenter() failed with status = %d, '%s'\n", fn, status.statusCode, status.statusDescription );
        goto failed;
      }

      deltaT = XLALGPSGetREAL8 ( &emit.te ) - refTimeREAL8;
      deltaT_int = (INT4)deltaT;

      ret->DeltaT_int->data[i] = deltaT_int;
      ret->DeltaT_rem->data[i]  = (REAL4)( deltaT - (REAL8)deltaT_int );
      ret->TdotM1->data[i] = (REAL4) ( emit.tDot - 1.0 );

    } /* for i < numSteps */

  /* finally: store the reference-time used into the output-structure */
  ret->refTime = refTime;

  goto success;

 failed:
  XLALDestroySSBtimes_REAL4 ( ret );
  XLAL_ERROR_NULL ( fn, XLAL_EFAILED );

 success:
  return ret;

} /* XLALGetSSBtimes_REAL4() */



/** Destruction of a ComputeFBuffer_REAL4 *contents*,
 * i.e. the multiSSB and multiAMcoeff, while the
 * buffer-container is not freed (which is why it's passed
 * by value and not by reference...) */
void
XLALEmptyComputeFBuffer_REAL4 ( ComputeFBuffer_REAL4 *cfb)
{
  if ( !cfb )
    return;

  XLALDestroyMultiSSBtimes_REAL4 ( cfb->multiSSB );
  cfb->multiSSB = NULL;

  XLALDestroyMultiAMCoeffs ( cfb->multiAMcoef );
  cfb->multiAMcoef = NULL;

  return;

} /* XLALDestroyComputeFBuffer_REAL4() */



/* ---------- pure REAL4 version of sin/cos lookup tables */

/** REAL4 version of sin_cos_LUT()
 *
 * Calculate sin(x) and cos(x) to roughly 1e-7 precision using
 * a lookup-table and Taylor-expansion.
 *
 * NOTE: this function will fail for arguments larger than
 * |x| > INT4_MAX = 2147483647 ~ 2e9 !!!
 *
 * return = 0: OK, nonzero=ERROR
 */
void
sin_cos_LUT_REAL4 (REAL4 *sinx, REAL4 *cosx, REAL4 x)
{
  sin_cos_2PI_LUT_REAL4 ( sinx, cosx, x * OOTWOPI_FLOAT );
  return;
}

/* initialize the global sin/cos lookup table */
void
init_sin_cos_LUT_REAL4 (void)
{
  UINT4 k;
  static int firstCall = 1;

  if ( !firstCall )
    return;

  for (k=0; k <= LUT_RES; k++)
    {
      sinVal[k] = sin( LAL_TWOPI * k / LUT_RES );
      cosVal[k] = cos( LAL_TWOPI * k / LUT_RES );
    }

  firstCall = 0;

  return;

} /* init_sin_cos_LUT_REAL4() */


/* REAL4 version of sin_cos_2PI_LUT() */
void
sin_cos_2PI_LUT_REAL4 (REAL4 *sin2pix, REAL4 *cos2pix, REAL4 x)
{
  REAL4 xt;
  INT4 i0;
  REAL4 d, d2;
  REAL4 ts, tc;

  /* we only need the fractional part of 'x', which is number of cylces,
   * this was previously done using
   *   xt = x - (INT4)x;
   * which is numerically unsafe for x > LAL_INT4_MAX ~ 2e9
   * for saftey we therefore rather use modf(), even if that
   * will be somewhat slower...
   */
  xt = x - (INT8)x;	/* xt in (-1, 1) */

  if ( xt < 0.0 )
    xt += 1.0f;			/* xt in [0, 1 ) */

#ifndef LAL_NDEBUG
  if ( xt < 0.0 || xt > 1.0 )
    {
      XLALPrintError("\nFailed numerics in sin_cos_2PI_LUT_REAL4(): xt = %f not in [0,1)\n\n", xt );
      XLAL_ERROR_VOID ( "sin_cos_2PI_LUT_REAL4()", XLAL_EFPINEXCT );
    }
#endif

  i0 = (INT4)( xt * LUT_RES + 0.5f );	/* i0 in [0, LUT_RES ] */
  d = d2 = LAL_TWOPI * (xt - OO_LUT_RES * i0);
  d2 *= 0.5f * d;

  ts = sinVal[i0];
  tc = cosVal[i0];

  /* use Taylor-expansions for sin/cos around LUT-points */
  (*sin2pix) = ts + d * tc - d2 * ts;
  (*cos2pix) = tc - d * ts - d2 * tc;

  return;

} /* sin_cos_2PI_LUT_REAL4() */
