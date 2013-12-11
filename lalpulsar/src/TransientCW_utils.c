/*
 * Copyright (C) 2010 Reinhard Prix, Stefanos Giampanis
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
 * \author R. Prix, S. Giampanis
 * \file
 * \brief
 * Some helper functions useful for "transient CWs", mostly applying transient window
 * functions.
 *
 */
#include "config.h"

/* System includes */
#include <math.h>

/* LAL-includes */
#include <lal/XLALError.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/LogPrintf.h>
#include <lal/LALString.h>

#include <lal/ProbabilityDensity.h>
#include <lal/TransientCW_utils.h>


/* ----- MACRO definitions ---------- */
#define SQ(x) ((x)*(x))
#define LAL_INT4_MAX 2147483647

/* ---------- internal prototypes ---------- */
/* empty struct initializers */
const transientCandidate_t empty_transientCandidate;
const transientWindow_t empty_transientWindow;
const transientWindowRange_t empty_transientWindowRange;
const transientFstatMap_t empty_transientFstatMap;

/* ----- module-local fast lookup-table handling of negative exponentials ----- */

/**
 * Lookup-table for negative exponentials e^(-x)
 * Holds an array 'data' of 'length' for values e^(-x) for x in the range [0, xmax]
 */
#define EXPLUT_XMAX 	20.0	// LUT down to e^(-20) = 2.0612e-09
#define EXPLUT_LENGTH 	2000	// number of LUT values to pre-compute
static gsl_vector *expLUT = NULL; 	/**< module-global lookup-table for negative exponentials e^(-x) */
#define EXPLUT_DXINV  ((EXPLUT_LENGTH)/(EXPLUT_XMAX))	// 1/dx with dx = xmax/length

static int XLALCreateExpLUT ( void );	/* only ever used internally, destructor is in exported API */

static const char *transientWindowNames[TRANSIENT_LAST] =
  {
    [TRANSIENT_NONE]	 	= "none",
    [TRANSIENT_RECTANGULAR]	= "rect",
    [TRANSIENT_EXPONENTIAL]	= "exp"
  };

/* ==================== function definitions ==================== */
/// Parse a transient window name string into the corresponding transientWindowType
int
XLALParseTransientWindowName ( const char *windowName )
{
  XLAL_CHECK ( windowName != NULL, XLAL_EINVAL );

  // convert input window-name into lower-case first
  char windowNameLC [ strlen(windowName) + 1 ];
  strcpy ( windowNameLC, windowName );
  XLALStringToLowerCase ( windowNameLC );

  int winType = -1;
  for ( UINT4 j=0; j < TRANSIENT_LAST; j ++ )
    {
      if ( !strcmp ( windowNameLC, transientWindowNames[j] ) ) {
        winType = j;
        break;
      }
    } // j < TRANSIENT_LAST

  if ( winType == -1 )
    {
      XLALPrintError ("Invalid transient Window-name '%s', allowed are (case-insensitive): [%s", windowName, transientWindowNames[0] );
      for ( UINT4 j = 1; j < TRANSIENT_LAST; j ++ ) {
        XLALPrintError (", %s", transientWindowNames[j] );
      }
      XLALPrintError ("]\n");
      XLAL_ERROR ( XLAL_EINVAL );
    } // if windowName not valid

  return winType;

} // XLALParseTransientWindowName()

/**
 * Helper-function to determine the total timespan of
 * a transient CW window, ie. the earliest and latest timestamps
 * of non-zero window function.
 */
int
XLALGetTransientWindowTimespan ( UINT4 *t0,				/**< [out] window start-time */
                                 UINT4 *t1,				/**< [out] window end-time */
                                 transientWindow_t transientWindow	/**< [in] window-parameters */
                                 )
{
  /* check input consistency */
  if ( !t0 || !t1 ) {
    XLALPrintError ("%s: invalid NULL input 't0=%p', 't1=%p'\n", __func__, t0, t1 );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  UINT4 win_t0 = transientWindow.t0;
  UINT4 win_tau = transientWindow.tau;

  switch ( transientWindow.type )
    {
    case TRANSIENT_NONE:
      (*t0) = 0;
      (*t1) = LAL_INT4_MAX;
      break;
    case TRANSIENT_EXPONENTIAL:
      (*t0) = win_t0;
      /* for given tau, what Tcoh does should the exponential window cover?
       * for speed reasons we want to truncate Tcoh = tau * TRANSIENT_EXP_EFOLDING
       * with the e-folding factor chosen such that the window-value
       * is practically negligible after that, where it will be set to 0
       */
      (*t1) = (UINT4)( win_t0 + TRANSIENT_EXP_EFOLDING * win_tau + 0.5 );
      break;
    case TRANSIENT_RECTANGULAR:
      (*t0) = win_t0;
      (*t1) = win_t0 + win_tau;
      break;
    default:
      XLALPrintError ("invalid transient window type %d not in [%d, %d].\n",
                      transientWindow.type, TRANSIENT_NONE, TRANSIENT_LAST -1 );
      XLAL_ERROR ( XLAL_EINVAL );

    } /* switch window-type */

  return XLAL_SUCCESS;

} /* XLALGetTransientWindowTimespan() */


/**
 * apply a "transient CW window" described by TransientWindowParams to the given
 * timeseries
 */
int
XLALApplyTransientWindow ( REAL4TimeSeries *series,		/**< input timeseries to apply window to */
                           transientWindow_t transientWindow	/**< transient-CW window to apply */
                           )
{
  /* check input consistency */
  if ( !series || !series->data ){
    XLALPrintError ("%s: Illegal NULL in input timeseries!\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* special time-saving break-condition: do nothing for window=none */
  if ( transientWindow.type == TRANSIENT_NONE )
    return XLAL_SUCCESS;

  /* deal with non-trivial windows */
  UINT4 ts_t0 = series->epoch.gpsSeconds;
  UINT4 ts_length = series->data->length;
  REAL8 ts_dt = series->deltaT;

  UINT4 t0, t1;
  if ( XLALGetTransientWindowTimespan ( &t0, &t1, transientWindow ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALGetTransientWindowTimespan() failed.\n", __func__ );
    XLAL_ERROR_REAL8 ( XLAL_EFUNC );
  }

  UINT4 i;
  switch ( transientWindow.type )
    {
    case TRANSIENT_RECTANGULAR:
      for ( i = 0; i < ts_length; i ++ )
        {
          UINT4 ti = (UINT4) ( ts_t0 + i * ts_dt + 0.5 );	// integer round: floor(x+0.5)
          if ( ti < t0 || ti > t1 ) { // outside rectangular window: set to zero
            series->data->data[i] = 0;
          } // otherwise do nothing
        } /* for i < length */
      break;

    case TRANSIENT_EXPONENTIAL:
      for ( i = 0; i < ts_length; i ++ )
        {
          UINT4 ti = (UINT4) ( ts_t0 + i * ts_dt + 0.5 );
          REAL8 win = XLALGetExponentialTransientWindowValue ( ti, t0, t1, transientWindow.tau );
          series->data->data[i] *= win;
        } /* for i < length */
      break;

    default:
      XLALPrintError ("%s: invalid transient window type %d not in [%d, %d].\n",
                      __func__, transientWindow.type, TRANSIENT_NONE, TRANSIENT_LAST -1 );
      XLAL_ERROR ( XLAL_EINVAL );
      break;

    } /* switch (window.type) */

  return XLAL_SUCCESS;

} /* XLALApplyTransientWindow() */


/**
 * apply transient window to give multi noise-weights, associated with given
 * multi timestamps
 */
int
XLALApplyTransientWindow2NoiseWeights ( MultiNoiseWeights *multiNoiseWeights,	/**< [in/out] noise weights to apply transient window to */
                                        const MultiLIGOTimeGPSVector *multiTS,	/**< [in] associated timestamps of noise-weights */
                                        transientWindow_t transientWindow	/**< [in] transient window parameters */
                                        )
{
  UINT4 numIFOs, X;
  UINT4 numTS, i;

  /* check input consistency */
  if ( !multiNoiseWeights || multiNoiseWeights->length == 0 ) {
    XLALPrintError ("%s: empty or NULL input 'multiNoiseWeights'.\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( !multiTS || multiTS->length == 0 ) {
    XLALPrintError ("%s: empty or NULL input 'multiTS'.\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  numIFOs = multiNoiseWeights->length;
  if ( multiTS->length != numIFOs ) {
    XLALPrintError ("%s: inconsistent numIFOs between 'multiNoiseWeights' (%d) and 'multiTS' (%d).\n", __func__, numIFOs, multiTS->length );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* special time-saving break-condition: do nothing for window=none */
  if ( transientWindow.type == TRANSIENT_NONE )
    return XLAL_SUCCESS;

  /* deal with non-trivial windows */
  UINT4 t0, t1;
  if ( XLALGetTransientWindowTimespan ( &t0, &t1, transientWindow ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALGetTransientWindowTimespan() failed.\n", __func__ );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* loop over all detectors X */
  for ( X = 0; X < numIFOs; X ++ )
    {
      numTS = multiNoiseWeights->data[X]->length;

      if ( multiTS->data[X]->length != numTS ) {
        XLALPrintError ("%s: inconsistent number of timesteps 'multiNoiseWeights[%d]' (%d) and 'multiTS[%d]' (%d).\n", __func__, X, numTS, X, multiTS->data[X]->length );
        XLAL_ERROR ( XLAL_EINVAL );
      }

      switch ( transientWindow.type )
        {
        case TRANSIENT_RECTANGULAR:
          for ( i=0; i < numTS; i ++ )
            {
              UINT4 ti = multiTS->data[X]->data[i].gpsSeconds;
              REAL8 win = XLALGetRectangularTransientWindowValue ( ti, t0, t1 );
              multiNoiseWeights->data[X]->data[i] *= win;
            } /* for i < length */
          break;

        case TRANSIENT_EXPONENTIAL:
          for ( i=0; i < numTS; i ++ )
            {
              UINT4 ti = multiTS->data[X]->data[i].gpsSeconds;
              REAL8 win = XLALGetExponentialTransientWindowValue ( ti, t0, t1, transientWindow.tau );
              multiNoiseWeights->data[X]->data[i] *= win;
            } /* for i < length */
          break;

        default:
          XLALPrintError ("%s: invalid transient window type %d not in [%d, %d].\n",
                          __func__, transientWindow.type, TRANSIENT_NONE, TRANSIENT_LAST -1 );
          XLAL_ERROR ( XLAL_EINVAL );
          break;

        } /* switch (window.type) */

    } /* for X < numIFOs */

  return XLAL_SUCCESS;

} /* XLALApplyTransientWindow2NoiseWeights() */


/**
 * Turn pulsar doppler-params into a single string that can be used for filenames
 * The format is
 * tRefNNNNNN_RAXXXXX_DECXXXXXX_FreqXXXXX[_f1dotXXXXX][_f2dotXXXXx][_f3dotXXXXX]
 */
CHAR*
XLALPulsarDopplerParams2String ( const PulsarDopplerParams *par )
{
#define MAXLEN 1024
  CHAR buf[MAXLEN];
  CHAR *ret = NULL;
  int len;
  UINT4 i;

  if ( !par )
    {
      LogPrintf(LOG_CRITICAL, "%s: NULL params input.\n", __func__ );
      XLAL_ERROR_NULL( XLAL_EDOM);
    }

  len = snprintf ( buf, MAXLEN, "tRef%09d_RA%.9g_DEC%.9g_Freq%.15g",
		      par->refTime.gpsSeconds,
		      par->Alpha,
		      par->Delta,
		      par->fkdot[0] );
  if ( len >= MAXLEN )
    {
      LogPrintf(LOG_CRITICAL, "%s: filename-size (%d) exceeded maximal length (%d): '%s'!\n", __func__, len, MAXLEN, buf );
      XLAL_ERROR_NULL( XLAL_EDOM);
    }

  for ( i = 1; i < PULSAR_MAX_SPINS; i++)
    {
      if ( par->fkdot[i] )
	{
	  CHAR buf1[MAXLEN];
	  len = snprintf ( buf1, MAXLEN, "%s_f%ddot%.7g", buf, i, par->fkdot[i] );
	  if ( len >= MAXLEN )
	    {
	      LogPrintf(LOG_CRITICAL, "%s: filename-size (%d) exceeded maximal length (%d): '%s'!\n", __func__, len, MAXLEN, buf1 );
	      XLAL_ERROR_NULL( XLAL_EDOM);
	    }
	  strcpy ( buf, buf1 );
	}
    }

  if ( par->orbit )
    {
      LogPrintf(LOG_NORMAL, "%s: orbital params not supported in Doppler-filenames yet\n", __func__ );
    }

  len = strlen(buf) + 1;
  if ( (ret = LALMalloc ( len )) == NULL )
    {
      LogPrintf(LOG_CRITICAL, "%s: failed to LALMalloc(%d)!\n", __func__, len );
      XLAL_ERROR_NULL( XLAL_ENOMEM);
    }

  strcpy ( ret, buf );

  return ret;
} /* PulsarDopplerParams2String() */


/**
 * Compute transient-CW Bayes-factor B_SG = P(x|HypS)/P(x|HypG)  (where HypG = Gaussian noise hypothesis),
 * marginalized over start-time and timescale of transient CW signal, using given type and parameters
 * of transient window range.
 *
 * Note: this function is a C-implemention, partially based-on/inspired-by Stefanos Giampanis'
 * original matlab implementation of this search function.
 *
 * Note2: if window->type == none, uses a single rectangular window covering all the data.
 */
REAL8
XLALComputeTransientBstat ( transientWindowRange_t windowRange,		/**< [in] type and parameters specifying transient window range */
                            const transientFstatMap_t *FstatMap		/**< [in] pre-computed transient-Fstat map F_mn over {t0, tau} ranges */
                            )
{
  /* ----- check input consistency */
  if ( !FstatMap || !FstatMap->F_mn ) {
    XLALPrintError ("%s: invalid NULL input 'FstatMap' or 'FstatMap->F_mn'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( windowRange.type >= TRANSIENT_LAST ) {
    XLALPrintError ("%s: unknown window-type (%d) passes as input. Allowed are [0,%d].\n", __func__, windowRange.type, TRANSIENT_LAST-1);
    XLAL_ERROR ( XLAL_EINVAL );
  }

  /* ----- step through F_mn array subtract maxF and sum e^{F_mn - maxF}*/
  /*
   * The maximum-likelihood Fmax is globally subtracted from F_mn, and stored separatedly in the struct, because in most
   * expressions it is numerically more robust to compute e^(F_mn - Fmax), which at worst can underflow, while
   * e^F_mn can overflow (for F>~700). The constant offset e^Fmax is irrelevant for posteriors (normalization constant), or
   * can be handled separately, eg by computing log(B) = Fmax + log(sum e^(Fmn-Fmax)) for the Bayes-factor.
   */
  UINT4 N_t0Range  = FstatMap->F_mn->size1;
  UINT4 N_tauRange = FstatMap->F_mn->size2;
  UINT4 m, n;
  REAL8 sum_eB = 0;
  for ( m=0; m < N_t0Range; m ++ )
    {
      for ( n=0; n < N_tauRange; n ++ )
        {
          REAL8 DeltaF = FstatMap->maxF - gsl_matrix_get ( FstatMap->F_mn, m, n );	// always >= 0, exactly ==0 at {m,n}_max

          //sum_eB += exp ( - DeltaF );
          sum_eB += XLALFastNegExp ( DeltaF );

        } /* for n < N_tauRange */

    } /* for m < N_t0Range */

  /* combine this to final log(Bstat) result with proper normalization (assuming rhohMax=1) : */

  REAL8 logBhat = FstatMap->maxF + log ( sum_eB );	// unnormalized Bhat

  REAL8 normBh = 70.0 / ( N_t0Range * N_tauRange );	// normalization factor assuming rhohMax=1

  /* final normalized Bayes factor, assuming rhohMax=1 */
  /* NOTE: correct this for different rhohMax by adding "- 4 * log(rhohMax)" to logB*/
  REAL8 logBstat = log ( normBh ) +  logBhat;	/* - 4.0 * log ( rhohMax ) */

  // printf ( "\n\nlogBhat = %g, normBh = %g, log(normBh) = %g\nN_t0Range = %d, N_tauRange=%d\n\n", logBhat, normBh, log(normBh), N_t0Range, N_tauRange );

  /* free mem */
  XLALDestroyExpLUT();

  /* ----- return ----- */
  return logBstat;

} /* XLALComputeTransientBstat() */

/**
 * Compute transient-CW posterior (normalized) on start-time t0, using given type and parameters
 * of transient window range.
 *
 * NOTE: the returned pdf has a number of sample-points N_t0Range given by the size
 * of the input matrix  FstatMap (namely N_t0Range = t0Band / dt0)
 *
 */
pdf1D_t *
XLALComputeTransientPosterior_t0 ( transientWindowRange_t windowRange,		/**< [in] type and parameters specifying transient window range */
                                   const transientFstatMap_t *FstatMap		/**< [in] pre-computed transient-Fstat map F_mn over {t0, tau} ranges */
                                   )
{
  /* ----- check input consistency */
  if ( !FstatMap || !FstatMap->F_mn ) {
    XLALPrintError ("%s: invalid NULL input 'FstatMap' or 'FstatMap->F_mn'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  if ( windowRange.type >= TRANSIENT_LAST ) {
    XLALPrintError ("%s: unknown window-type (%d) passes as input. Allowed are [0,%d].\n", __func__, windowRange.type, TRANSIENT_LAST-1);
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* ----- step through F_mn array subtract maxF and sum e^{F_mn - maxF}*/
  /*
   * It is numerically more robust to marginalize over e^(F_mn - Fmax), which at worst can underflow, while
   * e^F_mn can overflow (for F>~700). The constant offset e^Fmax is irrelevant for posteriors (normalization constant).
   */
  UINT4 N_t0Range  = FstatMap->F_mn->size1;
  UINT4 N_tauRange = FstatMap->F_mn->size2;

  REAL8 t0 = windowRange.t0;
  REAL8 t1 = t0 + windowRange.t0Band;

  pdf1D_t *ret;

  /* ----- handle special cases: 1) point-like support, 2) uniform pdf-value over 1 bin ----- */
  if ( N_t0Range == 1 && (windowRange.t0Band == 0) )
    {
      if ( (ret = XLALCreateSingularPDF1D ( t0 )) == NULL ) {
        XLALPrintError ("%s: failed to create singular pdf for t0 = %g\n", __func__, t0 );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }
      return ret;
    } /* if singular pdf in t0 */
  if ( (N_t0Range == 1) && (windowRange.t0Band > 0) )
    {
      if ( (ret = XLALCreateUniformPDF1D ( t0, t1 )) == NULL ) {
        XLALPrintError ( "%s: failed to created unform pdf over [%g, %g]\n", __func__, t0, t1 );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }
      return ret;
    } /* if uniform pdf over small band t0Band */

  /* ----- general N>1 point pdf case ----- */
  if ( ( ret = XLALCreateDiscretePDF1D ( t0, t1, N_t0Range )) == NULL ) {
    XLALPrintError ("%s: XLALCreateDiscretePDF1D() failed with xlalErrno = %d\n", __func__, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  UINT4 m, n;
  for ( m=0; m < N_t0Range; m ++ )	/* loop over start-times t0 */
    {
      REAL8 sum_eF = 0;
      for ( n=0; n < N_tauRange; n ++ )	/* loop over timescales tau */
        {
          REAL8 DeltaF = FstatMap->maxF - gsl_matrix_get ( FstatMap->F_mn, m, n );	// always >= 0, exactly ==0 at {m,n}_max

          //sum_eB += exp ( - DeltaF );
          sum_eF += XLALFastNegExp ( DeltaF );

        } /* for n < N_tauRange */

      ret->probDens->data[m] = sum_eF;

    } /* for m < N_t0Range */

  /* free mem */
  XLALDestroyExpLUT();

  /* normalize this PDF */
  if ( XLALNormalizePDF1D ( ret ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: failed to normalize posterior pdf ..\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* ----- return ----- */
  return ret;

} /* XLALComputeTransientPosterior_t0() */

/**
 * Compute transient-CW posterior (normalized) on timescale tau, using given type and parameters
 * of transient window range.
 *
 * NOTE: the returned pdf has a number of sample-points N_tauRange given by the size
 * of the input matrix  FstatMap (namely N_tauRange = tauBand / dtau)
 *
 */
pdf1D_t *
XLALComputeTransientPosterior_tau ( transientWindowRange_t windowRange,		/**< [in] type and parameters specifying transient window range */
                                    const transientFstatMap_t *FstatMap		/**< [in] pre-computed transient-Fstat map F_mn over {t0, tau} ranges */
                                    )
{
  /* ----- check input consistency */
  if ( !FstatMap || !FstatMap->F_mn ) {
    XLALPrintError ("%s: invalid NULL input 'FstatMap' or 'FstatMap->F_mn'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  if ( windowRange.type >= TRANSIENT_LAST ) {
    XLALPrintError ("%s: unknown window-type (%d) passes as input. Allowed are [0,%d].\n", __func__, windowRange.type, TRANSIENT_LAST-1);
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* ----- step through F_mn array subtract maxF and sum e^{F_mn - maxF}*/
  /*
   * It is numerically more robust to marginalize over e^(F_mn - Fmax), which at worst can underflow, while
   * e^F_mn can overflow (for F>~700). The constant offset e^Fmax is irrelevant for posteriors (normalization constant).
   */
  UINT4 N_t0Range  = FstatMap->F_mn->size1;
  UINT4 N_tauRange = FstatMap->F_mn->size2;

  REAL8 tau0 = windowRange.tau;
  REAL8 tau1 = tau0 + windowRange.tauBand;

  pdf1D_t *ret;

  /* ----- handle special cases: 1) point-like support, 2) uniform pdf-value over 1 bin ----- */
  if ( N_tauRange == 1 && (windowRange.tauBand == 0) )
    {
      if ( (ret = XLALCreateSingularPDF1D ( tau0 )) == NULL ) {
        XLALPrintError ("%s: failed to create singular pdf for tau0 = %g\n", __func__, tau0 );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }
      return ret;
    } /* if singular pdf in tau */
  if ( (N_tauRange == 1) && (windowRange.tauBand > 0) )
    {
      if ( (ret = XLALCreateUniformPDF1D ( tau0, tau1 )) == NULL ) {
        XLALPrintError ( "%s: failed to created unform pdf over [%g, %g]\n", __func__, tau0, tau1 );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }
      return ret;
    } /* if uniform pdf over small band tauBand */

  /* ----- general N>1 point pdf case ----- */
  if ( ( ret = XLALCreateDiscretePDF1D ( tau0, tau1, N_tauRange )) == NULL ) {
    XLALPrintError ("%s: XLALCreateDiscretePDF1D() failed with xlalErrno = %d\n", __func__, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  UINT4 m, n;
  for ( n=0; n < N_tauRange; n ++ )	/* loop over timescales tau */
    {
      REAL8 sum_eF = 0;
      for ( m=0; m < N_t0Range; m ++ )	/* loop over start-times t0 */
        {
          REAL8 DeltaF = FstatMap->maxF - gsl_matrix_get ( FstatMap->F_mn, m, n );	// always >= 0, exactly ==0 at {m,n}_max

          //sum_eB += exp ( - DeltaF );
          sum_eF += XLALFastNegExp ( DeltaF );

        } /* for m < N_t0Range */

      ret->probDens->data[n] = sum_eF;

    } /* for n < N_tauRange */

  /* free mem */
  XLALDestroyExpLUT();

  /* normalize this PDF */
  if ( XLALNormalizePDF1D ( ret ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: failed to normalize posterior pdf ..\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* ----- return ----- */
  return ret;

} /* XLALComputeTransientPosterior_tau() */



/**
 * Function to compute transient-window "F-statistic map" over start-time and timescale {t0, tau}.
 * Returns a 2D matrix F_mn, with m = index over start-times t0, and n = index over timescales tau,
 * in steps of dt0 in [t0, t0+t0Band], and dtau in [tau, tau+tauBand] as defined in transientWindowRange
 *
 * Note: if window->type == none, we compute a single rectangular window covering all the data.
 *
 * Note2: if the experimental switch useFReg is true, returns FReg=F - log(D) instead of F. This option is of
 * little practical interest, except for demonstrating that marginalizing (1/D)e^F is *less* sensitive
 * than marginalizing e^F (see transient methods-paper [in prepartion])
 *
 */
transientFstatMap_t *
XLALComputeTransientFstatMap ( const MultiFstatAtomVector *multiFstatAtoms, 	/**< [in] multi-IFO F-statistic atoms */
                               transientWindowRange_t windowRange,		/**< [in] type and parameters specifying transient window range to search */
                               BOOLEAN useFReg					/**< [in] experimental switch: compute FReg = F - log(D) instead of F */
                               )
{
  /* check input consistency */
  if ( !multiFstatAtoms || !multiFstatAtoms->data || !multiFstatAtoms->data[0]) {
    XLALPrintError ("%s: invalid NULL input.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  if ( windowRange.type >= TRANSIENT_LAST ) {
    XLALPrintError ("%s: unknown window-type (%d) passes as input. Allowed are [0,%d].\n", __func__, windowRange.type, TRANSIENT_LAST-1);
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* ----- pepare return container ----- */
  transientFstatMap_t *ret;
  if ( (ret = XLALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    XLALPrintError ("%s: XLALCalloc(1,%s) failed.\n", sizeof(*ret) );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  /* ----- first combine all multi-atoms into a single atoms-vector with *unique* timestamps */
  FstatAtomVector *atoms;
  UINT4 TAtom = multiFstatAtoms->data[0]->TAtom;
  UINT4 TAtomHalf = TAtom/2;	/* integer division */

  if ( (atoms = XLALmergeMultiFstatAtomsBinned ( multiFstatAtoms, TAtom )) == NULL ) {
    XLALPrintError ("%s: XLALmergeMultiFstatAtomsSorted() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  UINT4 numAtoms = atoms->length;
  /* actual data spans [t0_data, t0_data + numAtoms * TAtom] in steps of TAtom */
  UINT4 t0_data = atoms->data[0].timestamp;
  UINT4 t1_data = atoms->data[numAtoms-1].timestamp + TAtom;

  /* ----- special treatment of window_type = none ==> replace by rectangular window spanning all the data */
  if ( windowRange.type == TRANSIENT_NONE )
    {
      windowRange.type = TRANSIENT_RECTANGULAR;
      windowRange.t0 = t0_data;
      windowRange.t0Band = 0;
      windowRange.dt0 = TAtom;	/* irrelevant */
      windowRange.tau = numAtoms * TAtom;
      windowRange.tauBand = 0;
      windowRange.dtau = TAtom;	/* irrelevant */
    }

  /* NOTE: indices {i,j} enumerate *actual* atoms and their timestamps t_i, while the
   * indices {m,n} enumerate the full grid of values in [t0_min, t0_max]x[Tcoh_min, Tcoh_max] in
   * steps of deltaT. This allows us to deal with gaps in the data in a transparent way.
   *
   * NOTE2: we operate on the 'binned' atoms returned from XLALmergeMultiFstatAtomsBinned(),
   * which means we can safely assume all atoms to be lined up perfectly on a 'deltaT' binned grid.
   *
   * The mapping used will therefore be {i,j} -> {m,n}:
   *   m = offs_i  / deltaT		= start-time offset from t0_min measured in deltaT
   *   n = Tcoh_ij / deltaT		= duration Tcoh_ij measured in deltaT,
   *
   * where
   *   offs_i  = t_i - t0_min
   *   Tcoh_ij = t_j - t_i + deltaT
   *
   */

  /* We allocate a matrix  {m x n} = t0Range * TcohRange elements
   * covering the full timerange the transient window-range [t0,t0+t0Band]x[tau,tau+tauBand]
   */
  UINT4 N_t0Range  = (UINT4) floor ( windowRange.t0Band / windowRange.dt0 ) + 1;
  UINT4 N_tauRange = (UINT4) floor ( windowRange.tauBand / windowRange.dtau ) + 1;

  if ( ( ret->F_mn = gsl_matrix_calloc ( N_t0Range, N_tauRange )) == NULL ) {
    XLALPrintError ("%s: failed ret->F_mn = gsl_matrix_calloc ( %d, %d )\n", __func__, N_tauRange, N_t0Range );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  transientWindow_t win_mn;
  win_mn.type = windowRange.type;
  ret->maxF = 0;	// keep track of loudest F-stat point
  UINT4 m, n;
  /* ----- OUTER loop over start-times [t0,t0+t0Band] ---------- */
  for ( m = 0; m < N_t0Range; m ++ ) /* m enumerates 'binned' t0 start-time indices  */
    {
      /* compute Fstat-atom index i_t0 in [0, numAtoms) */
      win_mn.t0 = windowRange.t0 + m * windowRange.dt0;
      INT4 i_tmp = ( win_mn.t0 - t0_data + TAtomHalf ) / TAtom;	// integer round: floor(x+0.5)
      if ( i_tmp < 0 ) i_tmp = 0;
      UINT4 i_t0 = (UINT4)i_tmp;
      if ( i_t0 >= numAtoms ) i_t0 = numAtoms - 1;

      /* ----- INNER loop over timescale-parameter tau ---------- */
      REAL8 Ad=0, Bd=0, Cd=0, Fa_re=0, Fa_im=0, Fb_re=0, Fb_im=0;
      UINT4 i_t1_last = i_t0;

      for ( n = 0; n < N_tauRange; n ++ )
        {
          /* translate n into an atoms end-index for this search interval [t0, t0+Tcoh],
           * giving the index range of atoms to sum over
           */
          win_mn.tau = windowRange.tau + n * windowRange.dtau;

          /* get end-time t1 of this transient-window search */
          UINT4 t0, t1;
          if ( XLALGetTransientWindowTimespan ( &t0, &t1, win_mn ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: XLALGetTransientWindowTimespan() failed.\n", __func__ );
            XLAL_ERROR_NULL ( XLAL_EFUNC );
          }

          /* compute window end-time Fstat-atom index i_t1 in [0, numAtoms) */
          i_tmp = ( t1 - t0_data + TAtomHalf ) / TAtom  - 1;	// integer round: floor(x+0.5)
          if ( i_tmp < 0 ) i_tmp = 0;
          UINT4 i_t1 = (UINT4)i_tmp;
          if ( i_t1 >= numAtoms ) i_t1 = numAtoms - 1;

          /* protection against degenerate 1-atom case: (this implies D=0 and therefore F->inf) */
          if ( i_t1 == i_t0 ) {
            XLALPrintError ("%s: encountered a single-atom Fstat-calculation. This is degenerate and cannot be computed!\n", __func__ );
            XLALPrintError ("Window-values m=%d (t0=%d=t0_data + %d), n=%d (tau=%d) ==> t1_data - t0 = %d\n",
                            m, win_mn.t0, i_t0 * TAtom, n, win_mn.tau, t1_data - win_mn.t0 );
            XLALPrintError ("The most likely cause is that your t0-range covered all of your data: t0 must stay away *at least* 2*TAtom from the end of the data!\n");
            XLAL_ERROR_NULL ( XLAL_EDOM );
          }

          /* now we have two valid atoms-indices [i_t0, i_t1] spanning our Fstat-window to sum over,
           * using weights according to the window-type
           */
          switch ( windowRange.type )
            {
            case TRANSIENT_RECTANGULAR:
#if 0
              /* 'vanilla' unoptimized method, for sanity checks with 'optimized' method */
              Ad=0; Bd=0; Cd=0; Fa_re=0; Fa_im=0; Fb_re=0; Fb_im=0;
              for ( UINT4 i = i_t0; i <= i_t1; i ++ )
#else
              /* special optimiziation in the rectangular-window case: just add on to previous tau values
               * ie re-use the sum over [i_t0, i_t1_last] from the pevious tau-loop iteration
               */
              for ( UINT4 i = i_t1_last; i <= i_t1; i ++ )
#endif
                {
                  FstatAtom *thisAtom_i = &atoms->data[i];

                  /* now add on top of previous values, summed from [i_t0, i_t1_last] */
                  Ad += thisAtom_i->a2_alpha;
                  Bd += thisAtom_i->b2_alpha;
                  Cd += thisAtom_i->ab_alpha;

                  Fa_re += crealf(thisAtom_i->Fa_alpha);
                  Fa_im += cimagf(thisAtom_i->Fa_alpha);

                  Fb_re += crealf(thisAtom_i->Fb_alpha);
                  Fb_im += cimagf(thisAtom_i->Fb_alpha);

                } /* for i = i_t1_last : i_t1 */

              i_t1_last = i_t1 + 1;		/* keep track of up to where we summed for the next iteration */

              break;

            case TRANSIENT_EXPONENTIAL:
              /* reset all values */
              Ad=0; Bd=0; Cd=0; Fa_re=0; Fa_im=0; Fb_re=0; Fb_im=0;

              for ( UINT4 i = i_t0; i <= i_t1; i ++ )
                {
                  FstatAtom *thisAtom_i = &atoms->data[i];
                  UINT4 t_i = thisAtom_i->timestamp;

                  REAL8 win_i;
                  win_i = XLALGetExponentialTransientWindowValue ( t_i, t0, t1, win_mn.tau );

                  REAL8 win2_i = win_i * win_i;

                  Ad += thisAtom_i->a2_alpha * win2_i;
                  Bd += thisAtom_i->b2_alpha * win2_i;
                  Cd += thisAtom_i->ab_alpha * win2_i;

                  Fa_re += crealf(thisAtom_i->Fa_alpha) * win_i;
                  Fa_im += cimagf(thisAtom_i->Fa_alpha) * win_i;

                  Fb_re += crealf(thisAtom_i->Fb_alpha) * win_i;
                  Fb_im += cimagf(thisAtom_i->Fb_alpha) * win_i;

                } /* for i in [i_t0, i_t1] */
              break;

            default:
              XLALPrintError ("%s: invalid transient window type %d not in [%d, %d].\n",
                              __func__, windowRange.type, TRANSIENT_NONE, TRANSIENT_LAST -1 );
              XLAL_ERROR_NULL ( XLAL_EINVAL );
              break;

            } /* switch window.type */


          /* generic F-stat calculation from A,B,C, Fa, Fb */
          REAL8 DdInv = 1.0 / ( Ad * Bd - Cd * Cd );
          REAL8 F = DdInv * ( Bd * (SQ(Fa_re) + SQ(Fa_im) ) + Ad * ( SQ(Fb_re) + SQ(Fb_im) )
                              - 2.0 * Cd *( Fa_re * Fb_re + Fa_im * Fb_im )
                              );

          /* keep track of loudest F-stat value encountered over the m x n matrix */
          if ( F > ret->maxF )
            {
              ret->maxF = F;
              ret->t0_ML  = win_mn.t0;	/* start-time t0 corresponding to Fmax */
              ret->tau_ML = win_mn.tau;	/* timescale tau corresponding to Fmax */
            }

          /* if requested: use 'regularized' F-stat: log ( 1/D * e^F ) = F + log(1/D) */
          if ( useFReg )
            F += log( DdInv );

          /* and store this in Fstat-matrix as element {m,n} */
          gsl_matrix_set ( ret->F_mn, m, n, F );

        } /* for n in n[tau] : n[tau+tauBand] */

    } /* for m in m[t0] : m[t0+t0Band] */

  /* free internal mem */
  XLALDestroyFstatAtomVector ( atoms );

  /* return end product: F-stat map */
  return ret;

} /* XLALComputeTransientFstatMap() */




/**
 * Combine N Fstat-atoms vectors into a single 'canonical' binned and ordered atoms-vector.
 * The function pre-sums all atoms on a regular 'grid' of timestep bins deltaT covering the full data-span.
 * Atoms with timestamps falling into the bin i : [t_i, t_{i+1} ) are pre-summed and returned as atoms[i],
 * where t_i = t_0 + i * deltaT.
 *
 * Note: this pre-binning is equivalent to using a rectangular transient window on the deltaT timescale,
 * which is OK even with a different transient window, provided deltaT << transient-window timescale!
 *
 * Bins containing no atoms are returned with all values set to zero.
 */
FstatAtomVector *
XLALmergeMultiFstatAtomsBinned ( const MultiFstatAtomVector *multiAtoms, UINT4 deltaT )
{
  if ( !multiAtoms || !multiAtoms->length || !multiAtoms->data[0] || (deltaT==0) ) {
    XLALPrintError ("%s: invalid NULL input or deltaT=0.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  UINT4 numDet = multiAtoms->length;
  UINT4 X;
  UINT4 TAtom = multiAtoms->data[0]->TAtom;

  /* check consistency of time-step lengths between different IFOs */
  for ( X=0; X < numDet; X ++ ) {
    if ( multiAtoms->data[X]->TAtom != TAtom ) {
      XLALPrintError ("%s: Invalid input, atoms baseline TAtom=%d must be identical for all multiFstatAtomVectors (IFO=%d: TAtom=%d)\n",
                      __func__, TAtom, X, multiAtoms->data[X]->TAtom );
      XLAL_ERROR_NULL ( XLAL_EINVAL );
    }
  } /* for X < numDet */

  /* get earliest and latest atoms timestamps across all input detectors */
  UINT4 tMin = LAL_INT4_MAX - 1;
  UINT4 tMax = 0;
  for ( X=0; X < numDet; X ++ )
    {
      UINT4 numAtomsX = multiAtoms->data[X]->length;

      if ( multiAtoms->data[X]->data[0].timestamp < tMin )
        tMin = multiAtoms->data[X]->data[0].timestamp;

      if ( multiAtoms->data[X]->data[numAtomsX-1].timestamp > tMax )
        tMax = multiAtoms->data[X]->data[numAtomsX-1].timestamp;

    } /* for X < numDet */


  /* prepare 'canonical' binned atoms output vector */
  UINT4 NBinnedAtoms = (UINT4)floor( 1.0 * (tMax - tMin) / deltaT ) + 1; /* round up this way to make sure tMax is always included in the last bin */

  FstatAtomVector *atomsOut;
  if ( (atomsOut = XLALCreateFstatAtomVector ( NBinnedAtoms )) == NULL ) {	/* NOTE: these atoms are pre-initialized to zero already! */
    XLALPrintError ("%s: failed to XLALCreateFstatAtomVector ( %d )\n", __func__, NBinnedAtoms );
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  atomsOut->TAtom = deltaT;	/* output atoms-vector has new atoms baseline 'deltaT' */

  /* Step through all input atoms, and sum them together into output bins */
  for ( X=0; X < numDet; X ++ )
    {
      UINT4 i;
      UINT4 numAtomsX = multiAtoms->data[X]->length;
      for ( i=0; i < numAtomsX; i ++ )
        {
          FstatAtom *atom_X_i = &multiAtoms->data[X]->data[i];
          UINT4 t_X_i = atom_X_i -> timestamp;

          /* determine target bin-index j such that t_i in [ t_j, t_{j+1} )  */
          UINT4 j = (UINT4) floor ( 1.0 * ( t_X_i - tMin ) / deltaT );

          /* add atoms i to target atoms j */
          FstatAtom *destAtom = &atomsOut->data[j];
          destAtom->timestamp = tMin + i * deltaT;	/* set binned output atoms timestamp */

          destAtom->a2_alpha += atom_X_i->a2_alpha;
          destAtom->b2_alpha += atom_X_i->b2_alpha;
          destAtom->ab_alpha += atom_X_i->ab_alpha;
          destAtom->Fa_alpha += atom_X_i->Fa_alpha;
          destAtom->Fb_alpha += atom_X_i->Fb_alpha;

        } /* for i < numAtomsX */
    } /* for X < numDet */

  return atomsOut;

} /* XLALmergeMultiFstatAtomsBinned() */

/**
 * Write one line for given transient CW candidate into output file.
 *
 * NOTE: if input thisCand == NULL, we write a header comment-line explaining the fields
 *
 */
int
write_transientCandidate_to_fp ( FILE *fp, const transientCandidate_t *thisCand )
{
  /* sanity checks */
  if ( !fp ) {
    XLALPrintError ( "%s: invalid NULL filepointer input.\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }


  if ( thisCand == NULL )	/* write header-line comment */
    {
      fprintf (fp, "%%%% Freq[Hz]            Alpha[rad]          Delta[rad]          fkdot[1]  fkdot[2]  fkdot[3]    t0_ML[d]      tau_ML[d]    maxTwoF        logBstat      t0_MP[d]      tau_MP[d]\n");
    }
  else
    {
      if ( !thisCand->FstatMap ) {
        XLALPrintError ("%s: incomplete: transientCand->FstatMap == NULL!\n", __func__ );
        XLAL_ERROR ( XLAL_EINVAL );
      }
      UINT4 t0 = thisCand->windowRange.t0;
      REAL8 t0_d_ML = 1.0 * (thisCand->FstatMap->t0_ML - t0) / DAY24;
      REAL8 tau_d_ML= 1.0 *  thisCand->FstatMap->tau_ML / DAY24;
      REAL8 maxTwoF = 2.0 *  thisCand->FstatMap->maxF;
      REAL8 t0_d_MP = 1.0 * ( thisCand->t0_MP - t0 ) / DAY24;
      REAL8 tau_d_MP= 1.0 * thisCand->tau_MP / DAY24;

      fprintf (fp, "  %- 18.16f %- 19.16f %- 19.16f %- 9.6g %- 9.5g %- 9.5g    %-8.5f      %-8.5f    %- 11.8g    %- 11.8g    %-8.5f      %8.5f\n",
               thisCand->doppler.fkdot[0], thisCand->doppler.Alpha, thisCand->doppler.Delta,
               thisCand->doppler.fkdot[1], thisCand->doppler.fkdot[2], thisCand->doppler.fkdot[3],
               t0_d_ML, tau_d_ML, maxTwoF,
               thisCand->logBstat,
               t0_d_MP, tau_d_MP
               );
    }

  return XLAL_SUCCESS;

} /* write_TransCandidate_to_fp() */


/**
 * Write multi-IFO F-stat atoms 'multiAtoms' into output stream 'fstat'.
 */
int
write_MultiFstatAtoms_to_fp ( FILE *fp, const MultiFstatAtomVector *multiAtoms )
{
  UINT4 X, alpha;

  if ( !fp || !multiAtoms ) {
    XLALPrintError ( "%s: invalid NULL input.\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  fprintf ( fp, "%%%% GPS[s]     a^2(t_i)   b^2(t_i)  ab(t_i)            Fa(t_i)                  Fb(t_i)\n");

  for ( X=0; X < multiAtoms->length; X++ )
    {
      FstatAtomVector *thisAtomVector = multiAtoms->data[X];
      for ( alpha=0; alpha < thisAtomVector->length; alpha ++ )
	{
          FstatAtom *thisAtom = &thisAtomVector->data[alpha];
	  fprintf ( fp, "%d   % f  % f  %f    % f  % f     % f  % f\n",
		    thisAtom->timestamp,
		    thisAtom->a2_alpha,
		    thisAtom->b2_alpha,
		    thisAtom->ab_alpha,
		    crealf(thisAtom->Fa_alpha), cimagf(thisAtom->Fa_alpha),
		    crealf(thisAtom->Fb_alpha), cimagf(thisAtom->Fb_alpha)
		    );
	} /* for alpha < numSFTs */
    } /* for X < numDet */

  return XLAL_SUCCESS;

} /* write_MultiFstatAtoms_to_fp() */


/**
 * Generate an exponential lookup-table expLUT for e^(-x)
 * over the interval x in [0, xmax], using 'length' points.
 */
int
XLALCreateExpLUT ( void )
{
  /* create empty output LUT */
  gsl_vector *ret;
  if ( ( ret = gsl_vector_alloc ( EXPLUT_LENGTH + 1)) == NULL ) {
    XLALPrintError ("%s: failed to gsl_vector_alloc (%s)\n", __func__, EXPLUT_LENGTH +1 );
    XLAL_ERROR ( XLAL_ENOMEM );
  }

  /* fill output LUT */
  REAL8 dx = EXPLUT_XMAX / EXPLUT_LENGTH;
  UINT4 i;
  for ( i=0; i <= EXPLUT_LENGTH; i ++ )
    {
      REAL8 xi = i * dx;

      gsl_vector_set ( ret, i, exp( - xi ) );

    } /* for i < length() */

  /* 'return' this by setting the global vector */
  expLUT = ret;

  return XLAL_SUCCESS;

} /* XLALCreateExpLUT() */

/**
 * Destructor function for expLUT_t lookup table
 */
void
XLALDestroyExpLUT ( void )
{
  if ( !expLUT )
    return;

  gsl_vector_free ( expLUT );

  expLUT = NULL;

  return;

} /* XLALDestroyExpLUT() */

/**
 * Fast exponential function e^-x using lookup-table (LUT).
 * We need to compute exp(-x) for x >= 0, typically in a B-stat
 * integral of the form int e^-x dx: this means that small values e^(-x)
 * will not contribute much to the integral and are less important than
 * values close to 1. Therefore we pre-compute a LUT of e^(-x) for x in [0, xmax],
 * in Npoints points, and set e^(-x) = 0 for x < xmax.
 *
 * NOTE: if module-global expLUT=NULL, we create it here
 * NOTE: if argument is negative, we use math-lib exp(-x) instead of LUT
 */
REAL8
XLALFastNegExp ( REAL8 mx )
{
  if ( mx > EXPLUT_XMAX )	/* for values smaller than e^(-xmax) we truncate to 0 */
    return 0.0;

  if ( mx < 0 )
    return exp ( - mx  );

  /* if lookup table doesn't exist yet: generate it now */
  if ( !expLUT && ( XLALCreateExpLUT() != XLAL_SUCCESS) ) {
    XLAL_ERROR_REAL8 ( XLAL_EFUNC );
  }

  /* find index of closest point xp in LUT to xm */
  UINT4 i0 = (UINT4) ( mx * EXPLUT_DXINV + 0.5 );

  return gsl_vector_get ( expLUT, i0 );

} /* XLALFastNegExp() */

/**
 * Standard destructor for transientFstatMap_t
 * Fully NULL-robust as usual.
 */
void
XLALDestroyTransientFstatMap ( transientFstatMap_t *FstatMap )
{
  if ( !FstatMap )
    return;

  if ( FstatMap->F_mn )
    gsl_matrix_free ( FstatMap->F_mn );

  XLALFree ( FstatMap );

  return;

} /* XLALDestroyTransientFstatMap() */

/**
 * Standard destructor for transientCandidate_t
 * Fully NULL-robust as usual.
 */
void
XLALDestroyTransientCandidate ( transientCandidate_t *cand )
{
  if ( !cand )
    return;

  if ( cand->FstatMap )
    XLALDestroyTransientFstatMap ( cand->FstatMap );

  XLALFree ( cand );

  return;

} /* XLALDestroyTransientCandidate() */
