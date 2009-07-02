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

#include "ComputeFstatREAL4.h"
#include "../HierarchicalSearch.h"

NRCSID( COMPUTEFSTATC, "$Id$");

#ifdef EAH_OPTIMIZATION
#include "LocalOptimizationFlags.h"
#endif

/*----- SIN/COS Lookup table code inserted here -----*/
/* In particular this defines
  - a macro SINCOS_TRIM_X(y,x) which trims the value x to interval [0..2)
  - a function void local_sin_cos_2PI_LUT_init(void) that inits the lookup tables
  - a function local_sin_cos_2PI_LUT_trimmed(*sin,*cos,x) that uses the
    lookup tables to evaluate sin and cos values of 2*Pi*x
  - macros to be used in "hotloop" variants to interleave hotloop & sincos calculations
*/
#include "sincos.ci"

/*---------- local DEFINES ----------*/
#define TRUE  (1==1)
#define FALSE (1==0)

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
const PulsarSpinsREAL4 empty_PulsarSpinsREAL4;
const ComputeFBufferREAL4 empty_ComputeFBufferREAL4;
const ComputeFBufferREAL4V empty_ComputeFBufferREAL4V;
const FcomponentsREAL4 empty_FcomponentsREAL4;


/*---------- internal prototypes ----------*/
/* don't use finite() from win_lib.cpp here for performance reasons */
#ifdef _MSC_VER
#define finite _finite
#else
int finite(double x);
#endif
void init_sin_cos_LUT_REAL4(void);


/*==================== FUNCTION DEFINITIONS ====================*/

/** REAL4 and GPU-ready version of ComputeFStatFreqBand(), extended to loop over segments as well.
 *
 * Computes a vector of Fstatistic values for a number of frequency bins, for each segment
 */
int
XLALComputeFStatFreqBandVector (   REAL4FrequencySeriesVector *fstatBandV, 		/**< [out] Vector of Fstat frequency-bands */
                                   const PulsarDopplerParams *doppler,			/**< parameter-space point to compute F for */
                                   const MultiSFTVectorSequence *multiSFTsV, 		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
                                   const MultiNoiseWeightsSequence *multiWeightsV,	/**< noise-weights of all SFTs */
                                   const MultiDetectorStateSeriesSequence *multiDetStatesV,/**< 'trajectories' of the different IFOs */
                                   UINT4 Dterms,					/**< number of Dirichlet kernel Dterms to use */
                                   ComputeFBufferREAL4V *cfvBuffer			/**< buffer quantities that don't need to be recomputed */
                                   )
{
  static const char *fn = "XLALComputeFStatFreqBandVector()";

  UINT4 numBins, k;
  UINT4 numSegments, n;
  REAL8 f0, deltaF;
  REAL4 Fstat;
  REAL8 Freq0, Freq;

  /* check input consistency */
  if ( !doppler || !multiSFTsV || !multiWeightsV || !multiDetStatesV || !cfvBuffer ) {
    XLALPrintError ("%s: illegal NULL input pointer.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  if ( !fstatBandV || !fstatBandV->length ) {
    XLALPrintError ("%s: illegal NULL or empty output pointer 'fstatBandV'.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  numSegments = fstatBandV->length;

  if ( (multiSFTsV->length != numSegments) || (multiDetStatesV->length != numSegments ) ) {
    XLALPrintError ("%s: inconsistent number of segments between fstatBandV (%d), multiSFTsV(%d) and multiDetStatesV (%d)\n",
                    fn, numSegments, multiSFTsV->length, multiDetStatesV->length );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( multiWeightsV->length != numSegments ) {
    XLALPrintError ("%s: inconsistent number of segments between fstatBandV (%d) and multiWeightsV (%d)\n",
                    fn, numSegments, multiWeightsV->length );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  numBins = fstatBandV->data[0].data->length;
  f0      = fstatBandV->data[0].f0;
  deltaF  = fstatBandV->data[0].deltaF;

  /* a check that the f0 values from thisPoint and fstatVector are
   * at least close to each other -- this is only meant to catch
   * stupid errors but not subtle ones */
  if ( fabs(f0 - doppler->fkdot[0]) >= deltaF ) {
    XLALPrintError ("%s: fstatVector->f0 = %f differs from doppler->fkdot[0] = %f by more than deltaF = %g\n", fn, f0, doppler->fkdot[0], deltaF );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* ---------- prepare REAL4 version of PulsarSpins to be passed into core functions */
  PulsarSpinsREAL4 fkdot4 = empty_PulsarSpinsREAL4;
  UINT4 s;
  Freq = Freq0 = doppler->fkdot[0];
  fkdot4.FreqMain = (INT4)Freq0;
  /* fkdot[0] now only carries the *remainder* of Freq wrt FreqMain */
  fkdot4.fkdot[0] = (REAL4)( Freq0 - (REAL8)fkdot4.FreqMain );
  /* the remaining spins are simply down-cast to REAL4 */
  for ( s=1; s < PULSAR_MAX_SPINS; s ++ )
    fkdot4.fkdot[s] = (REAL4) doppler->fkdot[s];

  /* make sure sin/cos lookup-tables are initialized */
  init_sin_cos_LUT_REAL4();

  /* ---------- Buffering quantities that don't need to be recomputed ---------- */

  /* check if for this skyposition and data, the SSB+AMcoefs were already buffered */
  if ( cfvBuffer->Alpha != doppler->Alpha ||
       cfvBuffer->Delta != doppler->Delta ||
       cfvBuffer->multiDetStatesV != multiDetStatesV ||
       cfvBuffer->numSegments != numSegments
       )
    { /* no ==> compute and buffer */

      SkyPosition skypos;
      skypos.system =   COORDINATESYSTEM_EQUATORIAL;
      skypos.longitude = doppler->Alpha;
      skypos.latitude  = doppler->Delta;

      /* prepare buffer */
      XLALEmptyComputeFBufferREAL4V ( cfvBuffer );

      /* store new precomputed quantities in buffer */
      cfvBuffer->Alpha = doppler->Alpha;
      cfvBuffer->Delta = doppler->Delta;
      cfvBuffer->multiDetStatesV = multiDetStatesV ;
      cfvBuffer->numSegments = numSegments;

      if ( (cfvBuffer->multiSSB4V = XLALCalloc ( numSegments, sizeof(*cfvBuffer->multiSSB4V) )) == NULL ) {
        XLALPrintError ("%s: XLALCalloc ( %d, %d) failed.\n", fn, numSegments, sizeof(*cfvBuffer->multiSSB4V) );
        XLAL_ERROR ( fn, XLAL_ENOMEM );
      }
      if ( (cfvBuffer->multiAMcoefV = XLALCalloc ( numSegments, sizeof(*cfvBuffer->multiAMcoefV) )) == NULL ) {
        XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
        XLALPrintError ("%s: XLALCalloc ( %d, %d) failed.\n", fn, numSegments, sizeof(*cfvBuffer->multiAMcoefV) );
        XLAL_ERROR ( fn, XLAL_ENOMEM );
      }

      for ( n=0; n < numSegments; n ++ )
        {
          /* compute new SSB timings over all segments */
          if ( (cfvBuffer->multiSSB4V[n] = XLALGetMultiSSBtimesREAL4 ( multiDetStatesV->data[n], doppler->Alpha, doppler->Delta, doppler->refTime)) == NULL ) {
            XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
            XLALPrintError ( "%s: XLALGetMultiSSBtimesREAL4() failed. xlalErrno = %d.\n", fn, xlalErrno );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }

          LALStatus status = empty_LALStatus;
          LALGetMultiAMCoeffs ( &status, &(cfvBuffer->multiAMcoefV[n]), multiDetStatesV->data[n], skypos );
          if ( status.statusCode ) {
            XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
            XLALPrintError ("%s: LALGetMultiAMCoeffs() failed with statusCode=%d, '%s'\n", fn, status.statusCode, status.statusDescription );
            XLAL_ERROR ( fn, XLAL_EFAILED );
          }

          /* apply noise-weights to Antenna-patterns and compute A,B,C */
          if ( XLALWeighMultiAMCoeffs ( cfvBuffer->multiAMcoefV[n], multiWeightsV->data[n] ) != XLAL_SUCCESS ) {
            XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
            XLALPrintError("%s: XLALWeighMultiAMCoeffs() failed with error = %d\n", fn, xlalErrno );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }

        } /* for n < numSegments */

    } /* if we could NOT reuse previously buffered quantites */


  /* loop over all segments and compute FstatVector over frequencies for each */
  for ( n = 0; n < numSegments; n ++ )
    {

      Freq = Freq0;	/* reset frequency to start value for next segment */

      /* loop over frequency values and fill up values in fstatVector */
      for ( k = 0; k < numBins; k++)
        {
          /* REAL4-split frequency value */
          fkdot4.FreqMain = (INT4)Freq;
          fkdot4.fkdot[0] = (REAL4)( Freq - (REAL8)fkdot4.FreqMain );

          /* call the core function to compute one multi-IFO F-statistic */

          /* >>>>> this function could run on the GPU device <<<<< */
          XLALCoreFstatREAL4 ( &Fstat, &fkdot4, multiSFTsV->data[n], cfvBuffer->multiSSB4V[n], cfvBuffer->multiAMcoefV[n], Dterms );

          if ( xlalErrno ) {
            XLALEmptyComputeFBufferREAL4V ( cfvBuffer );
            XLALPrintError ("%s: XLALCoreFstatREAL4() failed with errno = %d in loop n=%d, k=%d.\n", fn, xlalErrno, n, k );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }

          fstatBandV->data[n].data->data[k] = Fstat;

          /* increase frequency by deltaF */
          Freq += deltaF;

        } /* for k < numBins */

    } /* for n < numSegments */

  return XLAL_SUCCESS;

} /* XLALComputeFStatFreqBandVector() */



/** Host-bound 'driver' function for the central F-stat computation
 *  of a single F-stat value for one parameter-space point.
 *
 * This is a GPU-adapted replacement for ComputeFstat(), and implements a 'wrapper'
 * around the core F-stat routines that can be executed as kernels on a GPU device.
 *
 */
int
XLALDriverFstatREAL4 ( REAL4 *Fstat,	                 		/**< [out] Fstatistic value 'F' */
                     const PulsarDopplerParams *doppler, 		/**< parameter-space point to compute F for */
                     const MultiSFTVector *multiSFTs,    		/**< normalized (by DOUBLE-sided Sn!) data-SFTs of all IFOs */
                     const MultiNoiseWeights *multiWeights,		/**< noise-weights of all SFTs */
                     const MultiDetectorStateSeries *multiDetStates,	/**< 'trajectories' of the different IFOs */
                     UINT4 Dterms,       				/**< number of terms to use in Dirichlet kernel interpolation */
                     ComputeFBufferREAL4 *cfBuffer       		/**< CF-internal buffering structure */
                     )
{
  static const char *fn = "XLALDriverFstatREAL4()";

  MultiSSBtimesREAL4 *multiSSB4 = NULL;
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
      if ( (multiSSB4 = XLALGetMultiSSBtimesREAL4 ( multiDetStates, doppler->Alpha, doppler->Delta, doppler->refTime)) == NULL ) {
        XLALPrintError ( "%s: XLALGetMultiSSBtimesREAL4() failed. xlalErrno = %d.\n", fn, xlalErrno );
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
      XLALEmptyComputeFBufferREAL4 ( cfBuffer );

      cfBuffer->Alpha = doppler->Alpha;
      cfBuffer->Delta = doppler->Delta;
      cfBuffer->multiDetStates = multiDetStates ;

      cfBuffer->multiSSB = multiSSB4;
      cfBuffer->multiAMcoef = multiAMcoef;

    } /* if we could NOT reuse previously buffered quantites */


  /* ---------- prepare REAL4 version of PulsarSpins to be passed into core functions */
  PulsarSpinsREAL4 fkdot4 = empty_PulsarSpinsREAL4;
  UINT4 s;
  fkdot4.FreqMain = (INT4)doppler->fkdot[0];
  /* fkdot[0] now only carries the *remainder* of Freq wrt FreqMain */
  fkdot4.fkdot[0] = (REAL4)( doppler->fkdot[0] - (REAL8)fkdot4.FreqMain );
  /* the remaining spins are simply down-cast to REAL4 */
  for ( s=1; s < PULSAR_MAX_SPINS; s ++ )
    fkdot4.fkdot[s] = (REAL4) doppler->fkdot[s];


  /* call the core function to compute one multi-IFO F-statistic */

  /* >>>>> this function could run on the GPU device <<<<< */
  XLALCoreFstatREAL4 (Fstat, &fkdot4, multiSFTs, multiSSB4, multiAMcoef, Dterms );

  if ( xlalErrno ) {
    XLALPrintError ("%s: XLALCoreFstatREAL4() failed with errno = %d.\n", fn, xlalErrno );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;


} /* XLALDriverFstatREAL4() */


/** This function computes a multi-IFO F-statistic value for given frequency (+fkdots),
 * antenna-pattern functions, SSB-timings and data ("SFTs").
 *
 * It is *only* using single-precision quantities. The aim is that this function can easily
 * be turned into code to be run on a GPU device.
 *
 */
void
XLALCoreFstatREAL4 (REAL4 *Fstat,				/**< [out] multi-IFO F-statistic value 'F' */
                  PulsarSpinsREAL4 *fkdot4,		/**< [in] frequency and derivatives in REAL4-safe format */
                  const MultiSFTVector *multiSFTs,	/**< [in] input multi-IFO data ("SFTs") */
                  MultiSSBtimesREAL4 *multiSSB4,	/**< [in] multi-IFO SSB-timings in REAL4-safe format */
                  MultiAMCoeffs *multiAMcoef,		/**< [in] multi-IFO antenna-pattern functions */
                  UINT4 Dterms				/**< [in] parameter: number of Dterms to use in Dirichlet kernel */
                  )
{
  static const char *fn = "XLALCoreFstatREAL4()";

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
      FcomponentsREAL4 FcX;

      XLALComputeFaFbREAL4 (&FcX, multiSFTs->data[X], fkdot4, multiSSB4->data[X], multiAMcoef->data[X], Dterms);

#ifndef LAL_NDEBUG
      if ( xlalErrno ) {
        XLALPrintError ("%s: XALComputeFaFbREAL4() failed\n", fn );
        XLAL_ERROR_VOID ( fn, XLAL_EFUNC );
      }
      if ( !finite(FcX.Fa.re) || !finite(FcX.Fa.im) || !finite(FcX.Fb.re) || !finite(FcX.Fb.im) ) {
	XLALPrintError("%s: XLALComputeFaFbREAL4() returned non-finite: Fa_X=(%f,%f), Fb_X=(%f,%f) for X=%d\n",
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

} /* XLALCoreFstatREAL4() */



/** Revamped version of LALDemod() (based on TestLALDemod() in CFS).
 * Compute JKS's Fa and Fb, which are ingredients for calculating the F-statistic.
 *
 * Note: this is a single-precision version aimed for GPU parallelization.
 *
 */
void
XLALComputeFaFbREAL4 ( FcomponentsREAL4 *FaFb,		/**< [out] single-IFO Fa/Fb for this parameter-space point */
                       const SFTVector *sfts,		/**< [in] single-IFO input data ("SFTs") */
                       const PulsarSpinsREAL4 *fkdot4,	/**< [in] frequency (and derivatives) in REAL4-safe format */
                       const SSBtimesREAL4 *tSSB,	/**< [in] single-IFO SSB-timings in REAL4-safe format */
                       const AMCoeffs *amcoe,		/**< [in] single-IFO antenna-pattern coefficients */
                       UINT4 Dterms)			/**< [in] number of Dterms to use in Dirichlet kernel */
{
  static const char *fn = "XLALComputeFaFbREAL4()";

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
  REAL4 tau = 1.0f / df;
  REAL4 Freq = f0 + df;

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
        REAL4 T0rem = fmodf ( T0, tau );
        phi_alpha_rem = f0 * dT;
        phi_alpha_rem += T0rem * df;
        phi_alpha_rem += df * dT;
        Dphi_alpha_int = f0;
        Dphi_alpha_rem = df + Freq * TdotM1;

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
	/* the sin/cos LUT calculations are woven into the hotloop for performance reasons
	   sin_cos_2PI_LUT_REAL4 ( &imagQ, &realQ, - lambda_alpha );
	*/

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
      /* the sin/cos LUT calculations are woven into the hotloop for performance reasons
	 sin_cos_2PI_LUT_REAL4 ( &s_alpha, &c_alpha, kappa_star );
	 c_alpha -= 1.0f;
      */

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
#ifdef __GNUC__
	/** somehow the branch prediction of gcc-4.1.2 terribly failes
	    with the current case distinction in the hot-loop,
	    having a severe impact on runtime of the E@H Linux App.
	    So let's allow to give gcc a hint which path has a higher probablility */
	if (__builtin_expect((kappa_star > LD_SMALL4) && (kappa_star < 1.0 - LD_SMALL4), (0==0)))
#else
      if ( ( kappa_star > LD_SMALL4 ) && (kappa_star < 1.0 - LD_SMALL4) )
#endif
#ifdef AUTOVECT_HOTLOOP
#include "hotloop_autovect.ci"
#elif __ALTIVEC__
#include "hotloop_altivec.ci"
#elif __SSE__
#ifdef _MSC_VER
#include "hotloop_sse_msc.ci"
#else
#include "hotloop_precalc.ci"
#endif /* MSC_VER */
#else
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

	      pn = pn - 1.0f; 				/* p_(n+1) */
	      Sn = pn * Sn + qn * (*Xalpha_l).re;	/* S_(n+1) */
	      Tn = pn * Tn + qn * (*Xalpha_l).im;	/* T_(n+1) */
	      qn *= pn;					/* q_(n+1) */
	    } /* for l <= 2*Dterms */

          REAL4 qn_inv = 1.0f / qn;

	  U_alpha = Sn * qn_inv;
	  V_alpha = Tn * qn_inv;

#ifndef LAL_NDEBUG
	  if ( !finite(U_alpha) || !finite(V_alpha) || !finite(pn) || !finite(qn) || !finite(Sn) || !finite(Tn) ) {
	    XLAL_ERROR_VOID (fn, XLAL_EFPINVAL );
	  }
#endif

	  {
	    REAL8 kstar8 = kappa_star;
	    SINCOS_2PI_TRIMMED(&s_alpha, &c_alpha, kstar8);
	  }
	  c_alpha -= 1.0f;
	
 	  realXP = s_alpha * U_alpha - c_alpha * V_alpha;
 	  imagXP = c_alpha * U_alpha + s_alpha * V_alpha;
	} /* if |remainder| > LD_SMALL4 */

        {
	  REAL8 _lambda_alpha = -lambda_alpha;
	  SINCOS_TRIM_X (_lambda_alpha,_lambda_alpha);
	  SINCOS_2PI_TRIMMED( &imagQ, &realQ, _lambda_alpha );
        }


#endif /* hotloop variant */

      else

	{ /* otherwise: lim_{rem->0}P_alpha,k  = 2pi delta_{k,kstar} */
	  UINT4 ind0;

	  /* realQ/imagQ are calculated in the hotloop in the other case; in this case we have to do it too */
	  REAL8 _lambda_alpha = -lambda_alpha;
	  SINCOS_TRIM_X (_lambda_alpha,_lambda_alpha);
	  SINCOS_2PI_TRIMMED( &imagQ, &realQ, _lambda_alpha );

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

} /* XLALComputeFaFbREAL4() */



/** Destroy a MultiSSBtimesREAL4 structure.
 *  Note, this is "NULL-robust" in the sense that it will not crash
 *  on NULL-entries anywhere in this struct, so it can be used
 *  for failure-cleanup even on incomplete structs
 */
void
XLALDestroyMultiSSBtimesREAL4 ( MultiSSBtimesREAL4 *multiSSB )
{
  UINT4 X;

  if ( ! multiSSB )
    return;

  if ( multiSSB->data )
    {
      for ( X=0; X < multiSSB->length; X ++ )
	{
          XLALDestroySSBtimesREAL4 ( multiSSB->data[X] );
	} /* for X < numDetectors */
      LALFree ( multiSSB->data );
    }

  LALFree ( multiSSB );

  return;

} /* XLALDestroyMultiSSBtimesREAL4() */


/** Destroy a SSBtimesREAL4 structure.
 *  Note, this is "NULL-robust" in the sense that it will not crash
 *  on NULL-entries anywhere in this struct, so it can be used
 *  for failure-cleanup even on incomplete structs
 */
void
XLALDestroySSBtimesREAL4 ( SSBtimesREAL4 *tSSB )
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

} /* XLALDestroySSBtimesREAL4() */


/** Multi-IFO version of LALGetSSBtimesREAL4().
 * Get all SSB-timings for all input detector-series in REAL4 representation.
 *
 */
MultiSSBtimesREAL4 *
XLALGetMultiSSBtimesREAL4 ( const MultiDetectorStateSeries *multiDetStates, 	/**< [in] detector-states at timestamps t_i */
                            REAL8 Alpha, REAL8 Delta,				/**< source sky-location, in equatorial coordinates  */
                            LIGOTimeGPS refTime
                            )
{
  static const char *fn = "XLALGetMultiSSBtimesREAL4()";

  UINT4 X, numDetectors;
  MultiSSBtimesREAL4 *ret;

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
      if ( (ret->data[X] = XLALGetSSBtimesREAL4 (multiDetStates->data[X], Alpha, Delta, refTime)) == NULL ) {
        XLALPrintError ("%s: XLALGetSSBtimesREAL4() failed. xlalErrno = %d\n", fn, xlalErrno );
        goto failed;
      }

    } /* for X < numDet */

  goto success;

 failed:
  /* free all memory allocated so far */
  XLALDestroyMultiSSBtimesREAL4 ( ret );
  XLAL_ERROR_NULL (fn, XLAL_EFAILED );

 success:
  return ret;

} /* LALGetMultiSSBtimes() */



/** XLAL REAL4-version of LALGetSSBtimes()
 *
 */
SSBtimesREAL4 *
XLALGetSSBtimesREAL4 ( const DetectorStateSeries *DetectorStates,	/**< [in] detector-states at timestamps t_i */
                       REAL8 Alpha, REAL8 Delta,			/**< source sky-location, in equatorial coordinates  */
                       LIGOTimeGPS refTime
                       )
{
  static const char *fn = "XLALGetSSBtimesREAL4()";

  UINT4 numSteps, i;
  REAL8 refTimeREAL8;
  SSBtimesREAL4 *ret;

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
  XLALDestroySSBtimesREAL4 ( ret );
  XLAL_ERROR_NULL ( fn, XLAL_EFAILED );

 success:
  return ret;

} /* XLALGetSSBtimesREAL4() */



/** Destruction of a ComputeFBufferREAL4 *contents*,
 * i.e. the multiSSB and multiAMcoeff, while the
 * buffer-container is not freed (which is why it's passed
 * by value and not by reference...) */
void
XLALEmptyComputeFBufferREAL4 ( ComputeFBufferREAL4 *cfb)
{
  if ( !cfb )
    return;

  XLALDestroyMultiSSBtimesREAL4 ( cfb->multiSSB );
  cfb->multiSSB = NULL;

  XLALDestroyMultiAMCoeffs ( cfb->multiAMcoef );
  cfb->multiAMcoef = NULL;

  return;

} /* XLALDestroyComputeFBufferREAL4() */



/** Destruction of a ComputeFBufferREAL4V *contents*,
 * i.e. the arrays of multiSSB and multiAMcoeff, while the
 * buffer-container is not freed (which is why it's passed
 * by value and not by reference...) */
void
XLALEmptyComputeFBufferREAL4V ( ComputeFBufferREAL4V *cfbv )
{
  if ( !cfbv )
    return;

  UINT4 numSegments = cfbv->numSegments;
  UINT4 i;
  for (i=0; i < numSegments; i ++ )
    {
      XLALDestroyMultiSSBtimesREAL4 ( cfbv->multiSSB4V[i] );
      XLALDestroyMultiAMCoeffs ( cfbv->multiAMcoefV[i] );
    }

  XLALFree ( cfbv->multiSSB4V );
  cfbv->multiSSB4V = NULL;

  XLALFree ( cfbv->multiAMcoefV );
  cfbv->multiAMcoefV = NULL;

  return;

} /* XLALDestroyComputeFBufferREAL4V() */



/* ---------- pure REAL4 version of sin/cos lookup tables */

void init_sin_cos_LUT_REAL4(void) {
  static int firstcall = TRUE;
  if (firstcall) {
    /* init sin/cos lookup tables */
    local_sin_cos_2PI_LUT_init();
    
    /* write out optimization settings */

#ifdef EAH_OPTIMIZATION
    
    fprintf(stderr,
	    "\n$Revision$ REV:%s, OPT:%d, "
	    "SCVAR:%d, SCTRIM:%d, "
	    "HOTVAR:%d, HOTDIV:%d, "
	    "HGHPRE:%d, HGHBAT:%d\n",
	    EAH_OPTIMIZATION_REVISION, EAH_OPTIMIZATION,
	    EAH_SINCOS_VARIANT,  EAH_SINCOS_ROUND,
	    EAH_HOTLOOP_VARIANT, EAH_HOTLOOP_DIVS,
	    EAH_HOUGH_PREFETCH,  EAH_HOUGH_BATCHSIZE_LOG2);
    
#endif /* def EAH_OPTIMIZATION */

    firstcall = FALSE;
  }
}
