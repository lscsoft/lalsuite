/*
 * Copyright (C) 2009 Reinhard Prix
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

/*********************************************************************************/
/**
 * \author R. Prix
 * \file
 * \brief
 * Some helper functions useful for "transient CWs", mostly applying transient window
 * functions.
 *
 */

#ifndef _TRANSIENTCW_UTILS_H
#define _TRANSIENTCW_UTILS_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

#include <lal/ProbabilityDensity.h>

/* ---------- System includes ---------- */
/* gsl-includes */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/LogPrintf.h>
#include <lal/SFTutils.h>
#include <lal/PulsarDataTypes.h>
#include <lal/ComputeFstat.h>

/* ---------- exported API defines ---------- */

#define DAY24 (24 * 3600)	/* standard 24h day = 86400 seconds ==> this is what's used in the definition of 'tauDays' */

#define TRANSIENT_EXP_EFOLDING	3.0      /**< e-folding parameter for exponential window, after which we truncate
                                          * the window for efficiency. 3 e-foldings means we lose only
                                          * about e^(-2x3) ~1e-8 of signal power! */

/* ---------- exported API types ---------- */

/** Struct defining a range of transient windows */
typedef struct tagtransientWindowRange_t
{
  transientWindowType_t type;	/**< window-type: none, rectangular, exponential, .... */
  UINT4 t0;			/**< earliest GPS start-time 't0' in seconds */
  UINT4 t0Band;			/**< range of start-times 't0' to search, in seconds */
  UINT4 dt0;			/**< stepsize to search t0-range with, in seconds */
  UINT4 tau;			/**< shortest transient timescale tau in seconds */
  UINT4 tauBand;		/**< range of transient timescales tau to search, in seconds */
  UINT4 dtau;			/**< stepsize to search tau-range with, in seconds */
} transientWindowRange_t;

/**
 * Struct holding a transient-window "F-statistic map" over start-time and timescale {t0, tau}.
 * This contains a 2D matrix F_mn, with m = index over start-times t0, and n = index over timescales tau,
 * in steps of dt0 in [t0, t0+t0Band], and dtau in [tau, tau+tauBand] as defined in transientWindowRange.
 *
 */
typedef struct tagtransientFstatMap_t {
  gsl_matrix *F_mn;			/**< "payload" F-map: F_mn for t0_m = t0 + m*dt0, and tau_n = tau + n*dtau */
  REAL8 maxF;				/**< maximal F-value obtained over transientWindowRange */
  UINT4 t0_ML;				/**< maximum-likelihood estimator for start-time t0 of  max{2F} over transientWindowRange (in GPS seconds) */
  UINT4 tau_ML;				/**< maximum-likelihood estimator for duration Tcoh of max{2F} over the transientWindowRange (in seconds) */
} transientFstatMap_t;


/** Struct holding a transient CW candidate */
typedef struct tagtransientCandidate_t {
  PulsarDopplerParams doppler;		/**< Doppler params of this 'candidate' */
  transientWindowRange_t windowRange;	/**< type and parameters specifying the transient window range in {t0, tau} covered */
  transientFstatMap_t *FstatMap;	/**< F-statistic over transient-window range {t0, tau} AND ML-estimators { Fmax, t0_Fmax, tau_Fmax } */
  REAL8 logBstat;			/**< log of Bayes-factor, marginalized over transientWindowRange */
  REAL8 t0_MP;				/**< maximum-posterior estimate for t0 */
  REAL8 tau_MP;				/**< maximum-posterior estimate for tau */
} transientCandidate_t;

/* empty struct initializers */
extern const transientCandidate_t empty_transientCandidate;
extern const transientWindow_t empty_transientWindow;
extern const transientWindowRange_t empty_transientWindowRange;
extern const transientFstatMap_t empty_transientFstatMap;

/* ---------- exported API prototypes ---------- */
int XLALParseTransientWindowName ( const char *windowName );

int XLALGetTransientWindowTimespan ( UINT4 *t0, UINT4 *t1, transientWindow_t transientWindow );

int XLALApplyTransientWindow ( REAL4TimeSeries *series, transientWindow_t TransientWindowParams );

int XLALApplyTransientWindow2NoiseWeights ( MultiNoiseWeights *multiNoiseWeights,
                                            const MultiLIGOTimeGPSVector *multiTS,
                                            transientWindow_t TransientWindowParams );

int write_transientCandidate_to_fp ( FILE *fp, const transientCandidate_t *thisTransCand );


transientFstatMap_t *XLALComputeTransientFstatMap ( const MultiFstatAtomVector *multiFstatAtoms,
                                                    transientWindowRange_t windowRange,
                                                    BOOLEAN useFReg );

REAL8 XLALComputeTransientBstat ( transientWindowRange_t windowRange, const transientFstatMap_t *FstatMap );
pdf1D_t *XLALComputeTransientPosterior_t0  ( transientWindowRange_t windowRange, const transientFstatMap_t *FstatMap );
pdf1D_t *XLALComputeTransientPosterior_tau ( transientWindowRange_t windowRange, const transientFstatMap_t *FstatMap );


void XLALDestroyTransientFstatMap ( transientFstatMap_t *FstatMap );
void XLALDestroyTransientCandidate ( transientCandidate_t *cand );

/* these functions operate on the module-local lookup-table for negative-exponentials,
 * which will dynamically be generated on first use of XLALFastNegExp(), and can
 * be destroyed at any time using XLALDestroyExpLUT()
 */
REAL8 XLALFastNegExp ( REAL8 mx );
void XLALDestroyExpLUT( void );

/* ---------- Fstat-atoms related functions ----------*/
int write_MultiFstatAtoms_to_fp ( FILE *fp, const MultiFstatAtomVector *multiAtoms );
CHAR* XLALPulsarDopplerParams2String ( const PulsarDopplerParams *par );

FstatAtomVector *XLALmergeMultiFstatAtomsBinned ( const MultiFstatAtomVector *multiAtoms, UINT4 deltaT );


/* ---------- INLINE function definitions ---------- */

/**
 * Function to compute the value of a rectangular transient-window at a given timestamp.
 * This is the central function defining the rectangular window properties.
 */
static inline REAL8
XLALGetRectangularTransientWindowValue ( UINT4 timestamp,	/**< timestamp for which to compute window-value */
                                         UINT4 t0, 		/**< start-time of rectangular window */
                                         UINT4 t1		/**< end-time of rectangular window */
                                         )
{
  if ( timestamp < t0 || timestamp > t1 )
    return 0.0;
  else
    return 1.0;

} /* XLALGetRectangularTransientWindowValue() */

/**
 * Function to compute the value of an exponential transient-window at a given timestamp.
 *
 * This is the central function defining the exponential window properties.
 */
static inline REAL8
XLALGetExponentialTransientWindowValue ( UINT4 timestamp,	/**< timestamp for which to compute window-value */
                                         UINT4 t0, 		/**< start-time of exponential window */
                                         UINT4 t1, 		/**< end-time of exponential window */
                                         UINT4 tau		/**< characteristic time of the exponential window */
                                         )
{
  REAL8 ret;

  if ( timestamp < t0 || timestamp > t1 )
    ret = 0.0;
  else
    {
      REAL8 x = 1.0 * ( timestamp - t0 ) / tau;
      ret = XLALFastNegExp ( x );	// computes e^(-x)
    }

  return ret;

} /* XLALGetExponentialTransientWindowValue() */

/**
 * Function to compute the value of a given transient-window function at a given timestamp.
 *
 * This is a simple wrapper to the actual window-defining functions
 */
static inline REAL8
XLALGetTransientWindowValue ( UINT4 timestamp,	/**< timestamp for which to compute window-value */
                              UINT4 t0, 	/**< start-time of window */
                              UINT4 t1, 	/**< end-time of window */
                              UINT4 tau,	/**< characteristic time of window */
                              transientWindowType_t type /**< window type */
                              )
{
  REAL8 val;

  switch ( type )
    {
    case TRANSIENT_NONE:
      val = 1.0;
      break;

    case TRANSIENT_RECTANGULAR:
      val = XLALGetRectangularTransientWindowValue ( timestamp, t0, t1 );
      break;

    case TRANSIENT_EXPONENTIAL:
      val = XLALGetExponentialTransientWindowValue ( timestamp, t0, t1, tau );
      break;

    default:
      XLALPrintError ("invalid transient window type %d not in [%d, %d].\n",
                      type, TRANSIENT_NONE, TRANSIENT_LAST -1 );
      return -1;	/* cop out because we're in an inline function */

    } /* switch window-type */

  /* return result */
  return val;

} /* XLALGetTransientWindowValue() */




#ifdef  __cplusplus
}
#endif
/* C++ protection. */

#endif  /* Double-include protection. */
