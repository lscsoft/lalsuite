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
/** \author R. Prix
 * \file
 * \brief
 * Some helper functions useful for "transient CWs", mostly applying transient window
 * functions.
 *
 *********************************************************************************/

#ifndef _TRANSIENTCW_UTILS_H
#define _TRANSIENTCW_UTILS_H

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

#include "config.h"

/* ---------- System includes ---------- */

/* gsl-includes */
#define GSL_RANGE_CHECK_OFF 1
#define GSL_C99_INLINE 1
#define HAVE_INLINE 1

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* LAL-includes */
#include <lal/AVFactories.h>
#include <lal/LogPrintf.h>
#include <lal/SFTutils.h>

#include <lal/ComputeFstat.h>

/* ---------- exported API defines ---------- */

#define DAY24 (24 * 3600)	/* standard 24h day = 86400 seconds ==> this is what's used in the definition of 'tauDays' */

/** Struct to define parameters of a 'transient window' to be applied to obtain transient signals */
typedef enum {
  TRANSIENT_NONE = 0,		/**< Note: in this case the window-parameters will be ignored, and treated as rect={data},
                                 * i.e. a simple rectangular window covering all the data => this should always reproduce the
                                 * standard F-statistic computation.
                                 */

  TRANSIENT_RECTANGULAR = 1,	/**< standard rectangular window covering [t0, t0+tau] */

  TRANSIENT_EXPONENTIAL,	/**< exponentially decaying window e^{-t0/tau} starting at t0.
                                 * Note: we'll truncate this at some small (eg 3x) e-folding TRANSIENT_EXP_EFOLDING
                                 */
  TRANSIENT_LAST
} transientWindowType_t;

#define TRANSIENT_EXP_EFOLDING	3.0      /**< e-folding parameter for exponential window, after which we truncate
                                          * the window for efficiency. 3 e-foldings means we lose only
                                          * about e^(-2x3) ~1e-8 of signal power! */

/* ---------- exported API types ---------- */

/** Struct defining one transient window instance */
typedef struct
{
  transientWindowType_t type;	/**< window-type: none, rectangular, exponential, .... */
  UINT4 t0;			/**< GPS start-time 't0' */
  UINT4 tau;			/**< transient timescale tau in seconds */
} transientWindow_t;

/** Struct defining a range of transient windows */
typedef struct
{
  transientWindowType_t type;	/**< window-type: none, rectangular, exponential, .... */
  UINT4 t0;			/**< earliest GPS start-time 't0' in seconds */
  UINT4 t0Band;			/**< range of start-times 't0' to search, in seconds */
  UINT4 dt0;			/**< stepsize to search t0-range with, in seconds */
  UINT4 tau;			/**< shortest transient timescale tau in seconds */
  UINT4 tauBand;		/**< range of transient timescales tau to search, in seconds */
  UINT4 dtau;			/**< stepsize to search tau-range with, in seconds */
  gsl_matrix *exp_buffer;	/**< buffer matrix storing exp(-(t0-ti)/tau) values for this window-range (optional) */
} transientWindowRange_t;

/** Struct holding a transient CW candidate */
typedef struct {
  PulsarDopplerParams doppler;		/**< Doppler params of this 'candidate' */
  REAL8 twoFtotal;			/**< 2F obtained in the full search over all SFTs */
  REAL8 maxTwoF;			/**< maximal 2F value obtained over transientWindowRange */
  UINT4 t0offs_maxF;			/**< start-time offset from transient_t0 of max{2F} over transientWindowRange (in GPS seconds)*/
  UINT4 tau_maxF;			/**< duration Tcoh where max{2F} occurred over the transientWindowRange (in seconds) */
  REAL8 logBstat;			/**< log of Bayes-factor, marginalized over transientWindowRange */
} TransientCandidate_t;

/* empty struct initializers */
extern const TransientCandidate_t empty_TransientCandidate;
extern const transientWindow_t empty_transientWindow;
extern const transientWindowRange_t empty_transientWindowRange;

/* ---------- exported API prototypes ---------- */
int XLALGetTransientWindowTimespan ( UINT4 *t0, UINT4 *t1, transientWindow_t transientWindow );

int XLALApplyTransientWindow ( REAL4TimeSeries *series, transientWindow_t TransientWindowParams );

int XLALApplyTransientWindow2NoiseWeights ( MultiNoiseWeights *multiNoiseWeights,
                                            const MultiLIGOTimeGPSVector *multiTS,
                                            transientWindow_t TransientWindowParams );

int write_TransientCandidate_to_fp ( FILE *fp, const TransientCandidate_t *thisTransCand );

int XLALFillExpWindowBuffer ( transientWindowRange_t *windowRange );
int XLALComputeTransientBstat ( TransientCandidate_t *transientCand, const MultiFstatAtomVector *multiFstatAtoms, transientWindowRange_t windowRange, BOOLEAN useFReg );


/* ---------- Fstat-atoms related functions ----------*/
int write_MultiFstatAtoms_to_fp ( FILE *fp, const MultiFstatAtomVector *multiAtoms );
CHAR* XLALPulsarDopplerParams2String ( const PulsarDopplerParams *par );

FstatAtomVector *XLALmergeMultiFstatAtomsBinned ( const MultiFstatAtomVector *multiAtoms, UINT4 deltaT );


/* ---------- INLINE function definitions ---------- */

/** Function to compute the value of a rectangular transient-window at a given timestamp.
 *
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

/** Function to compute the value of an exponential transient-window at a given timestamp.
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
  if ( timestamp < t0 || timestamp > t1 )
    return 0.0;
  else
    return exp ( - 1.0 * ( timestamp - t0 ) / tau );

} /* XLALGetExponentialTransientWindowValue() */

/** Function to compute the value of a given transient-window function at a given timestamp.
 *
 * This is a simple wrapper to the actual window-defining functions
 */
static inline REAL8
XLALGetTransientWindowValue ( UINT4 timestamp,	/**< timestamp for which to compute window-value */
                              UINT4 t0, 	/**< start-time of window */
                              UINT4 t1, 	/**< end-time of window */
                              UINT4 tau,	/**< characteristic time of window */
                              transientWindowType_t type	/**< window type */
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
