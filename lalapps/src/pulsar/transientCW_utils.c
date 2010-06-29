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
#include "config.h"

#include "transientCW_utils.h"

/* System includes */
#include <math.h>

/* LAL-includes */
#include <lal/XLALError.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/LogPrintf.h>

/* ----- MACRO definitions ---------- */
#define SQ(x) ((x)*(x))

/* ---------- internal prototypes ---------- */
REAL4Vector *XLALGetTransientWindowVals ( const LIGOTimeGPSVector *tGPS, const transientWindow_t *TransientWindowParams );

/* empty struct initializers */
const TransientCandidate_t empty_TransientCandidate;


/* ==================== function definitions ==================== */

/** apply a "transient CW window" described by TransientWindowParams to the given
 * timeseries
 */
int
XLALApplyTransientWindow ( REAL4TimeSeries *series, transientWindow_t TransientWindowParams )
{
  const CHAR *fn = "XLALApplyTransientWindow()";
  UINT4 i;

  REAL8 ts_t0, ts_dt, ts_T;	/* start-time, stepsize and duration of input timeseries */
  REAL8 ti;
  INT4 i0, i1;			/* time-series index corresonding to start-time (and end-time) of transient window */

  if ( !series || !series->data ){
    XLALPrintError ("%s: Illegal NULL in input timeseries!\n", fn );
    return XLAL_EINVAL;
  }
  ts_t0 = XLALGPSGetREAL8 ( &series->epoch );
  ts_dt = series->deltaT;
  ts_T  = ts_dt * series->data->length;

  i0 = ( TransientWindowParams.t0 - ts_t0 ) / ts_dt;
  if ( i0 < 0 ) i0 = 0;

  switch ( TransientWindowParams.type )
    {
    case TRANSIENT_NONE:
      return XLAL_SUCCESS;	/* nothing to be done here */
      break;

    case TRANSIENT_RECTANGULAR:	/* standard 'rectangular window */
      for ( i = 0; i < (UINT4)i0; i ++ ) {
	series->data->data[i] = 0;
      }
      i1 = (TransientWindowParams.t0 + TransientWindowParams.tau - ts_t0 ) / ts_dt + 1;
      if ( i1 < 0 ) i1 = 0;
      if ( (UINT4)i1 >= series->data->length ) i1 = series->data->length - 1;

      for ( i = i1; i < series->data->length; i ++) {
	series->data->data[i] = 0;
      }
      break;

    case TRANSIENT_EXPONENTIAL:
      for ( i = 0; i < (UINT4)i0; i ++ ) {
	series->data->data[i] = 0;
      }
      ti = 0;
      for ( i=i0; i < series->data->length; i ++)
	{
	  REAL8 fact = exp( - ti / TransientWindowParams.tau );
	  ti += ts_dt;
	  series->data->data[i] *= fact;
	}
      break;

    default:
      XLALPrintError("Illegal transient-signal window type specified '%d'\n", TransientWindowParams.type );
      return XLAL_EINVAL;
      break;
    }

  return XLAL_SUCCESS;

} /* XLALApplyTransientWindow() */


/** Compute values of given transient window function at GPS times 'tGPS'.
 */
REAL4Vector *
XLALGetTransientWindowVals ( const LIGOTimeGPSVector *tGPS, const transientWindow_t *TransientWindowParams )
{
  static const char *fn = "XLALGetTransientWindowVals()";
  UINT4 i, numTS;
  REAL4Vector *ret;
  REAL4 t0, t1, tau;
  REAL8 ti;

  if ( !tGPS || tGPS->length == 0 || !TransientWindowParams ) {
    XLALPrintError ("%s: invalid NULL input.\n", fn );
    XLAL_ERROR_NULL( fn, XLAL_EINVAL );
  }
  numTS = tGPS->length;

  if ( ( ret = XLALCreateREAL4Vector ( numTS )) == NULL ) {
    XLALPrintError ( "%s: XLALCreateREAL4Vector(%d) failed.\n", fn, numTS );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  t0 = TransientWindowParams->t0;
  tau = TransientWindowParams->tau;
  t1 = t0 + tau;

  switch ( TransientWindowParams->type )
    {
    case TRANSIENT_NONE:
      for ( i = 0; i < numTS; i ++ )
        ret->data[i] = 1;
      break;

    case TRANSIENT_RECTANGULAR:		/* standard 'rectangular window */
      for ( i = 0; i < numTS; i ++ )
        {
          ti = XLALGPSGetREAL8 ( &tGPS->data[i] );
          if ( ( ti >= t0 ) && (ti <= t1) )
            ret->data[i] = 1;
          else
            ret->data[i] = 0;
        } /* for i < numTS */
      break;

    case TRANSIENT_EXPONENTIAL:
      for ( i=0; i < numTS; i ++ )
        {
          ti = XLALGPSGetREAL8 ( &tGPS->data[i] );
          if ( ti >= t0 )
            ret->data[i] = exp( (t0 - ti) / tau );
          else
            ret->data[i] = 0;
        } /* for i < numTS */
      break;

    default:
      XLALPrintError ("%s: invalid transient window type %d not in [%d, %d].\n",
                      fn, TransientWindowParams->type, TRANSIENT_NONE, TRANSIENT_LAST -1 );
      XLALDestroyREAL4Vector ( ret );
      XLAL_ERROR_NULL( fn, XLAL_EINVAL );

      break;

    } /* switch transient type */

  return ret;

} /* XLALGetTransientWindowVals() */


/** apply transient window to give multi noise-weights, associated with given
 * multi timestamps
 */
int
XLALApplyTransientWindow2NoiseWeights ( MultiNoiseWeights *multiNoiseWeights,	/**< [in/out] noise weights to apply transient window to */
                                        const MultiLIGOTimeGPSVector *multiTS,	/**< [in] associated timestamps of noise-weights */
                                        transientWindow_t TransientWindowParams	/**< [in] transient window parameters */
                                        )
{
  static const char *fn = "XLALApplyTransientWindow2NoiseWeights()";

  UINT4 numIFOs, X;
  UINT4 numTS, i;
  REAL4Vector *win;

  /* check input consistency */
  if ( !multiNoiseWeights || multiNoiseWeights->length == 0 ) {
    XLALPrintError ("%s: empty or NULL input 'multiNoiseWeights'.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( !multiTS || multiTS->length == 0 ) {
    XLALPrintError ("%s: empty or NULL input 'multiTS'.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  numIFOs = multiNoiseWeights->length;
  if ( multiTS->length != numIFOs ) {
    XLALPrintError ("%s: inconsistent numIFOs between 'multiNoiseWeights' (%d) and 'multiTS' (%d).\n", fn, numIFOs, multiTS->length );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  for ( X = 0; X < numIFOs; X ++ )
    {
      numTS = multiNoiseWeights->data[X]->length;

      if ( multiTS->data[X]->length != numTS ) {
        XLALPrintError ("%s: inconsistent number of timesteps 'multiNoiseWeights[%d]' (%d) and 'multiTS[%d]' (%d).\n", fn, X, numTS, X, multiTS->data[X]->length );
        XLAL_ERROR ( fn, XLAL_EINVAL );
      }

      if ( ( win = XLALGetTransientWindowVals ( multiTS->data[X], &TransientWindowParams )) == NULL ) {
        XLALPrintError ("%s: XLALGetTransientWindowVals() failed. xlalErrno = %d.\n", fn, xlalErrno );
        XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      for ( i=0; i < numTS; i ++ )
        {
          multiNoiseWeights->data[X]->data[i] *= win->data[i];
        } /* for i < numTS */

      XLALDestroyREAL4Vector ( win );

    } /* for X < numIFOs */

  return XLAL_SUCCESS;

} /* XLALApplyTransientWindow2NoiseWeights() */




int
XLALoutputMultiFstatAtoms ( FILE *fp, MultiFstatAtomVector *multiAtoms )
{
  const char *fn = __func__;
  UINT4 X, alpha;

  if ( !fp || !multiAtoms )
    XLAL_ERROR (fn, XLAL_EINVAL );

  fprintf ( fp, "%% GPS[s]              a2(t_i)     b2(t_i)            Fa(t_i)                 Fb(t_i)\n");

  for ( X=0; X < multiAtoms->length; X++ )
    {
      FstatAtomVector *thisAtomVector = multiAtoms->data[X];
      for ( alpha=0; alpha < thisAtomVector->length; alpha ++ )
	{
          FstatAtom *thisAtom = &thisAtomVector->data[alpha];
	  fprintf ( fp, "%d   % f  % f     % f  % f     % f  % f\n",
		    thisAtom->timestamp,
		    thisAtom->a2_alpha,
		    thisAtom->b2_alpha,
		    thisAtom->Fa_alpha.re, thisAtom->Fa_alpha.im,
		    thisAtom->Fb_alpha.re, thisAtom->Fb_alpha.im
		    );
	} /* for alpha < numSFTs */
    } /* for X < numDet */

  return XLAL_SUCCESS;
} /* XLALoutputMultiFstatAtoms() */

/** Turn pulsar doppler-params into a single string that can be used for filenames
 * The format is
 * tRefNNNNNN_RAXXXXX_DECXXXXXX_FreqXXXXX[_f1dotXXXXX][_f2dotXXXXx][_f3dotXXXXX]
 */
CHAR*
XLALPulsarDopplerParams2String ( const PulsarDopplerParams *par )
{
  const CHAR *fn = "XLALPulsarDopplerParams2String()";
#define MAXLEN 1024
  CHAR buf[MAXLEN];
  CHAR *ret = NULL;
  int len;
  UINT4 i;

  if ( !par )
    {
      LogPrintf(LOG_CRITICAL, "%s: NULL params input.\n", fn );
      XLAL_ERROR_NULL( fn, XLAL_EDOM);
    }

  len = snprintf ( buf, MAXLEN, "tRef%09d_RA%.9g_DEC%.9g_Freq%.15g",
		      par->refTime.gpsSeconds,
		      par->Alpha,
		      par->Delta,
		      par->fkdot[0] );
  if ( len >= MAXLEN )
    {
      LogPrintf(LOG_CRITICAL, "%s: filename-size (%d) exceeded maximal length (%d): '%s'!\n", fn, len, MAXLEN, buf );
      XLAL_ERROR_NULL( fn, XLAL_EDOM);
    }

  for ( i = 1; i < PULSAR_MAX_SPINS; i++)
    {
      if ( par->fkdot[i] )
	{
	  CHAR buf1[MAXLEN];
	  len = snprintf ( buf1, MAXLEN, "%s_f%ddot%.7g", buf, i, par->fkdot[i] );
	  if ( len >= MAXLEN )
	    {
	      LogPrintf(LOG_CRITICAL, "%s: filename-size (%d) exceeded maximal length (%d): '%s'!\n", fn, len, MAXLEN, buf1 );
	      XLAL_ERROR_NULL( fn, XLAL_EDOM);
	    }
	  strcpy ( buf, buf1 );
	}
    }

  if ( par->orbit )
    {
      LogPrintf(LOG_NORMAL, "%s: orbital params not supported in Doppler-filenames yet\n", fn );
    }

  len = strlen(buf) + 1;
  if ( (ret = LALMalloc ( len )) == NULL )
    {
      LogPrintf(LOG_CRITICAL, "%s: failed to LALMalloc(%d)!\n", fn, len );
      XLAL_ERROR_NULL( fn, XLAL_ENOMEM);
    }

  strcpy ( ret, buf );

  return ret;
} /* PulsarDopplerParams2String() */


/** Function to compute marginalized B-statistic over start-time and duration
 * of transient CW signal, using given type and parameters of transient window range.
 */
REAL8
XLALComputeTransientBstat ( const MultiFstatAtomVector *multiFstatAtoms,	/**< [in] multi-IFO F-statistic atoms */
                            transientWindowRange_t windowRange )	/**< [in] type and parameters specifying transient window range to search */
{
  const char *fn = __func__;

  REAL8 tau0 = windowRange.tau_min; // smallest search duration-window in seconds
  REAL8 tau1 = windowRange.tau_max; // longest search duration-window in seconds

  REAL8 t0 = windowRange.t0_min;	// earliest GPS start-time
  REAL8 t1 = windowRange.t0_max;	// latest GPS start-time


  /* check input consistency */
  if ( !multiFstatAtoms ) {
    XLALPrintError ("%s: invalid NULL input.\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EINVAL );
  }

  if ( t1 < t0 ) {
    XLALPrintError ("%s: invalid input arguments t0 (=%d), t1 (=%d): must t1>t0 \n", fn, t0, t1 );
    XLAL_ERROR_REAL8 ( fn, XLAL_EDOM );
  }

  if ( tau1 < tau0 ) {
    XLALPrintError ("%s: invalid input arguments tau0 (=%d), tau1 (=%d): must tau1>tau0 \n", fn, tau0, tau1 );
    XLAL_ERROR_REAL8 ( fn, XLAL_EDOM );
  }

  if ( windowRange.type != TRANSIENT_RECTANGULAR ) {
    XLALPrintError ("%s: Sorry, only rectangular window implemented right now!\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EDOM );
  }


  UINT4 t_i;        // t index (t-summation)
  REAL8 tau_i;      // tau index (tau-summation)
  REAL8 logBAYES = 0;  // return value of function
  return logBAYES;

} /* XLALComputeTransientBstat() */

/** Write one line for given transient CW candidate into output file.
 * If input candidate == NULL, write a header comment-line explaining fields
 */
int
write_TransientCandidate_to_fp ( FILE *fp, const TransientCandidate_t *thisCand )
{
  if ( !fp )
    return -1;

  if ( thisCand == NULL )	/* write header-line comment */
    fprintf (fp, "\n\n%%%%        fkdot[0]         Alpha[rad]         Delta[rad]  fkdot[1] fkdot[2] fkdot[3]     2F_full     t0_max    tau_max       2F_max     logBstat\n");
  else
    fprintf (fp, "%18.16g %18.16g %18.16g %8.6g %8.5g %8.5g  %11.9g  %09d  %09d  %11.9g  %11.9g\n",
             thisCand->doppler.fkdot[0], thisCand->doppler.Alpha, thisCand->doppler.Delta,
             thisCand->doppler.fkdot[1], thisCand->doppler.fkdot[2], thisCand->doppler.fkdot[3],
             thisCand->fullFstat,
             thisCand->maxt0, thisCand->maxtau, thisCand->maxFstat,
             thisCand->logBstat
             );

  return XLAL_SUCCESS;

} /* write_TransCandidate_to_fp() */

/** Combine N Fstat-atoms vectors into a single one, with timestamps ordered by increasing GPS time,
 * Atoms with identical timestamp are immediately merged into one, so the final timestamps list only
 * contains unique (and ordered) entries.
 */
FstatAtomVector *
XLALmergeMultiFstatAtomsSorted ( const MultiFstatAtomVector *multiAtoms )
{
  const char *fn = __func__;

  if ( !multiAtoms ) {
    XLALPrintError ("%s: invalid NULL input.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  UINT4 numDet = multiAtoms->length;
  UINT4 X;
  UINT4 maxNumAtoms = 0;	/* upper limit on total number of atoms: sum over all detectors (assumes no coincident timestamps) */
  for ( X=0; X < numDet; X ++ )
    maxNumAtoms += multiAtoms->data[X]->length;

  /* first allocate an atoms-vector with maxNumAtoms length, then truncate at the end when we're done*/
  FstatAtomVector *atoms;
  if ( (atoms = XLALCreateFstatAtomVector ( maxNumAtoms )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCreateFstatAtomVector ( %d )\n", fn, maxNumAtoms );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  /* simply combine all atoms-vector by concatentation first */
  

  return NULL;

} /* XLALmergeMultiFstatAtoms() */
