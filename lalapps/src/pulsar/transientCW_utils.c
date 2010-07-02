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
/** \author R. Prix, S. Giampanis
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
int compareAtoms(const void *in1, const void *in2);


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
 *
 * Note: this function is a C-implemention, partially based-on/inspired-by Stefanos Giampanis'
 * original matlab implementation of this search function.
 */
int
XLALComputeTransientBstat ( TransientCandidate_t *cand, 		/**< [out] transient candidate info */
                            const MultiFstatAtomVector *multiFstatAtoms,/**< [in] multi-IFO F-statistic atoms */
                            transientWindowRange_t windowRange )	/**< [in] type and parameters specifying transient window range to search */
{
  const char *fn = __func__;

  UINT4 tau_min = windowRange.tau_min; // smallest search duration-window in seconds
  UINT4 tau_max = windowRange.tau_max; // longest search duration-window in seconds

  UINT4 t0_min = windowRange.t0_min;	// earliest GPS start-time
  UINT4 t0_max = windowRange.t0_max;	// latest GPS start-time

  /* initialize empty return, in case sth goes wrong */
  TransientCandidate_t ret = empty_TransientCandidate;
  (*cand) = ret;

  /* check input consistency */
  if ( !multiFstatAtoms ) {
    XLALPrintError ("%s: invalid NULL input.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  if ( t0_max < t0_min ) {
    XLALPrintError ("%s: invalid input arguments t0 (=%d), t1 (=%d): must t1>t0 \n", fn, t0_min, t0_max );
    XLAL_ERROR ( fn, XLAL_EDOM );
  }

  if ( tau_max < tau_min ) {
    XLALPrintError ("%s: invalid input arguments tau0 (=%d), tau1 (=%d): must tau1>tau0 \n", fn, tau_min, tau_max );
    XLAL_ERROR ( fn, XLAL_EDOM );
  }

  if ( windowRange.type != TRANSIENT_RECTANGULAR ) {
    XLALPrintError ("%s: Sorry, only rectangular window implemented right now!\n", fn );
    XLAL_ERROR ( fn, XLAL_EDOM );
  }

  /* combine all multi-atoms into a single atoms-vector with *unique* timestamps */
  FstatAtomVector *atoms;
  if ( (atoms = XLALmergeMultiFstatAtomsSorted ( multiFstatAtoms )) == NULL ) {
    XLALPrintError ("%s: XLALmergeMultiFstatAtomsSorted() failed with code %d\n", fn, xlalErrno );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }
  UINT4 numAtoms = atoms->length;
  UINT4 deltaT = atoms->deltaT;

  /* find indices corresponding to t0_min, t0_max */
  UINT4 i_min=0;
  while ( (i_min < numAtoms) && (atoms->data[i_min].timestamp < t0_min) )
    i_min++;
  if ( i_min == atoms->length ) {
    XLALPrintError ( "%s: earliest start-time %d later than latest atoms timestamp %d\n", fn, t0_min, atoms->data[i_min-1].timestamp );
    XLAL_ERROR (fn, XLAL_EDOM );
  }

  UINT4 i_max = i_min;
  while ( (i_max < numAtoms) && (atoms->data[i_max].timestamp < t0_max) )
    i_max ++;
  if ( i_max > 0 ) i_max --;	/* we want last timestamp that still satifies t_max <= t0_max */

  /* It is often numerically impossible to compute e^F and sum these values, because of range-overflow
   * instead we first determine max{F_ij}, then compute the logB = log ( e^Fmax * sum_{ij} 1/Dij * e^{Fij - Fmax} )
   * which is logB = Fmax + log( sum_{ij} e^FReg_ij ), where FReg_ij = -log(Dij) + Fij - Fmax.
   * This avoids numerical problems.
   *
   * As we don't know Fmax before having computed the full matrix F_ij, we keep a list of
   * 'regularized' F-stats FReg_ij over the field of {t0, tau} values.
   * As we don't know exactly the size of that matrix (because of possible gaps in the data),
   *  we use the maximal (conservative) estimate of that size t0Range* tauRange elements
   */
  UINT4 t0Range = ( t0_max - t0_min ) / deltaT + 1;
  UINT4 tauRange = (tau_max - tau_min) / deltaT + 1;
  UINT4 maxNumSummands = t0Range * tauRange;
  REAL8 *regFList;	/* will be initialized to zeros ! */
  if ( ( regFList = XLALCalloc ( maxNumSummands, sizeof(REAL8) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( %d, sizeof(REAL8)\n", fn, maxNumSummands );
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }
  UINT4 counter = 0;

  /* ----- OUTER loop over start-times t0 ---------- */
  ret.maxFstat = 0;	// keep track of loudest 2F-value over {t0, tau} space
  UINT4 i;
  REAL8 norm = 1.0 / SQ(LAL_TWOPI);

  for ( i = i_min; i <= i_max; i ++ )
    {
      UINT4 t0_i = atoms->data[i].timestamp;

      /* ----- INNER loop over durations tau ---------- */
      REAL8 Ad=0, Bd=0, Cd=0, Fa_re=0, Fa_im=0, Fb_re=0, Fb_im=0;

      UINT4 j = i;
      while ( j < numAtoms )
        {
          FstatAtom *thisAtom = &atoms->data[j];
          UINT4 t_j = thisAtom->timestamp;
          UINT4 tau_j = t_j - t0_i + deltaT;

          if ( tau_j > tau_max )
            break;

          Ad += thisAtom->a2_alpha;
          Bd += thisAtom->b2_alpha;
          Cd += thisAtom->ab_alpha;

          Fa_re += thisAtom->Fa_alpha.re;
          Fa_im += thisAtom->Fa_alpha.im;

          Fb_re += thisAtom->Fb_alpha.re;
          Fb_im += thisAtom->Fb_alpha.im;

          if ( tau_j >= tau_min )
            {
              REAL8 Dd, twoF;
              Dd = Ad * Bd - Cd * Cd;

              twoF = norm * (2.0 / Dd) * ( Bd * (SQ(Fa_re) + SQ(Fa_im) ) + Ad * ( SQ(Fb_re) + SQ(Fb_im) )
                                           - 2.0 * Cd *( Fa_re * Fb_re + Fa_im * Fb_im )
                                           );

              if ( twoF > ret.maxFstat )
                {
                  ret.maxFstat = twoF;
                  ret.t0_maxF  = t0_i;
                  ret.tau_maxF = tau_j;
                }

              /* compute 'regularized' F-stat: log ( 1/D * e^F ) = -logD + F */
              regFList[counter] = - log( Dd ) + 0.5 * twoF;
              counter ++;

              if ( counter > maxNumSummands ) {
                XLALPrintError ("%s: something went badly wrong, or repr can't count! ... numSummands=%d > maxNumSummand=%d\n", fn, counter, maxNumSummands );
                XLAL_ERROR ( fn, XLAL_EBADLEN );
              }

            } 	/* if inside [tau_min, tau_max] => compute Fstat */

          j ++ ;

        } /* j < numAtoms */

    } /* for i in [i_min, i_max] */

  UINT4 numSummands = counter;
  /* now step through list of FReg_ij, subtract maxFstat and sum e^{FReg - Fmax}*/
  REAL8 sum_eB = 0;
  for ( i=0; i < numSummands; i ++ )
    sum_eB += exp ( regFList[i] - 0.5 * ret.maxFstat );

  /* combine this to final log(Bstat) result: */
  ret.logBstat = 0.5 * ret.maxFstat + 2.0 * log ( deltaT ) + log ( sum_eB );

  /* free mem */
  XLALDestroyFstatAtomVector ( atoms );
  XLALFree ( regFList );

  /* return */
  (*cand) = ret;
  return XLAL_SUCCESS;

} /* XLALComputeTransientBstat() */


/** Combine N Fstat-atoms vectors into a single one, with timestamps ordered by increasing GPS time,
 * Atoms with identical timestamp are immediately merged into one, so the final timestamps list only
 * contains unique (and ordered) entries.
 */
FstatAtomVector *
XLALmergeMultiFstatAtomsSorted ( const MultiFstatAtomVector *multiAtoms )
{
  const char *fn = __func__;

  if ( !multiAtoms || !multiAtoms->length || !multiAtoms->data[0] ) {
    XLALPrintError ("%s: invalid NULL input.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  UINT4 numDet = multiAtoms->length;
  UINT4 X;
  UINT4 maxNumAtoms = 0;	/* upper limit on total number of atoms: sum over all detectors (assumes no coincident timestamps) */
  UINT4 deltaT = multiAtoms->data[0]->deltaT;

  /* check consistency of time-step lengths between different IFOs */
  for ( X=0; X < numDet; X ++ ) {
    if ( multiAtoms->data[X]->deltaT != deltaT ) {
      XLALPrintError ("%s: Invalid input, timestep-length deltaT=%d must be identical for all multiFstatAtomVectors (IFO=%d: deltaT=%d)\n",
                      fn, deltaT, X, multiAtoms->data[X]->deltaT );
      XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
    }
  } /* for X < numDet */

  for ( X=0; X < numDet; X ++ )
    maxNumAtoms += multiAtoms->data[X]->length;

  /* first allocate an atoms-vector with maxNumAtoms length, then truncate as needed at the end */
  FstatAtomVector *atomsIn;
  if ( (atomsIn = XLALCreateFstatAtomVector ( maxNumAtoms )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCreateFstatAtomVector ( %d )\n", fn, maxNumAtoms );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  /* simply combine all atoms-vector by concatentation first */
  UINT4 offset = 0;
  for ( X=0; X < numDet; X ++ )
    {
      FstatAtomVector *thisVect = multiAtoms->data[X];
      memcpy ( atomsIn->data + offset, thisVect->data, thisVect->length * sizeof(*thisVect->data) );
      offset += thisVect->length;
    } /* for X < numDet */

  /* now sort by increasing GPS time */
  qsort( atomsIn->data, maxNumAtoms, sizeof(*atomsIn->data), compareAtoms );

  /* finally: step through and 'merge' equal-timestamp atoms */
  FstatAtomVector *atomsOut;
  if ( (atomsOut = XLALCreateFstatAtomVector ( maxNumAtoms )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCreateFstatAtomVector ( %d )\n", fn, maxNumAtoms );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  atomsOut->deltaT = deltaT;

  FstatAtom *destAtom = &atomsOut->data[0];
  /* handle first atom by hand */
  (*destAtom) = atomsIn->data[0];
  UINT4 counter=1;
  /* and step through the rest of them */
  UINT4 i;
  for ( i=1; i < maxNumAtoms; i ++ )
    {
      FstatAtom *srcAtom = &atomsIn->data[i];

      /* still same timestamp? merge entries */
      if ( srcAtom->timestamp == destAtom->timestamp )
        {
          destAtom->a2_alpha += srcAtom->a2_alpha;
          destAtom->b2_alpha += srcAtom->b2_alpha;
          destAtom->ab_alpha += srcAtom->ab_alpha;
          destAtom->Fa_alpha.re += srcAtom->Fa_alpha.re;
          destAtom->Fa_alpha.im += srcAtom->Fa_alpha.im;
          destAtom->Fb_alpha.re += srcAtom->Fb_alpha.re;
          destAtom->Fb_alpha.im += srcAtom->Fb_alpha.im;
        } /* if same timestamp */
      else
        {
          counter ++;
          destAtom ++;
          (*destAtom) = (*srcAtom);
        } /* add new timestamp atom */

    } /* for i < maxNumAtoms */

  /* free concat vector */
  XLALDestroyFstatAtomVector ( atomsIn );

  /* now resize output Vector to that actually needed */
  UINT4 newsize = counter * sizeof(*atomsOut->data);
  atomsOut->length = counter;
  if ( (atomsOut->data = XLALRealloc ( atomsOut->data, newsize )) == NULL ) {
    XLALPrintError ("%s: failed to XLALRealloc() atomsOut to new size %d\n", fn, newsize );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
  }

  return atomsOut;

} /* XLALmergeMultiFstatAtoms() */

/* comparison function for atoms: sort by GPS time */
int
compareAtoms(const void *in1, const void *in2)
{
  const FstatAtom *atom1 = (const FstatAtom*)in1;
  const FstatAtom *atom2 = (const FstatAtom*)in2;

  if ( atom1->timestamp < atom2->timestamp )
    return -1;
  else if ( atom1->timestamp == atom2->timestamp )
    return 0;
  else
    return 1;

} /* compareAtoms() */

/** Write one line for given transient CW candidate into output file.
 * Note: if input candidate == NULL, write a header comment-line explaining fields
 */
int
write_TransientCandidate_to_fp ( FILE *fp, const TransientCandidate_t *thisCand )
{
  const char *fn = __func__;

  if ( !fp ) {
    XLALPrintError ( "%s: invalid NULL filepointer input.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  if ( thisCand == NULL )	/* write header-line comment */
    fprintf (fp, "\n\n%%%%        fkdot[0]         Alpha[rad]         Delta[rad]  fkdot[1] fkdot[2] fkdot[3]     2F_full     t0_Fmax   tau_Fmax      2F_max     logBstat\n");
  else
    fprintf (fp, "%18.16g %18.16g %18.16g %8.6g %8.5g %8.5g  %11.9g  %9d  %9d  %11.9g  %11.9g\n",
             thisCand->doppler.fkdot[0], thisCand->doppler.Alpha, thisCand->doppler.Delta,
             thisCand->doppler.fkdot[1], thisCand->doppler.fkdot[2], thisCand->doppler.fkdot[3],
             thisCand->fullFstat,
             thisCand->t0_maxF, thisCand->tau_maxF, thisCand->maxFstat,
             thisCand->logBstat
             );

  return XLAL_SUCCESS;

} /* write_TransCandidate_to_fp() */

/** Write multi-IFO F-stat atoms 'multiAtoms' into output stream 'fstat'.
 */
int
write_MultiFstatAtoms_to_fp ( FILE *fp, const MultiFstatAtomVector *multiAtoms )
{
  const char *fn = __func__;
  UINT4 X, alpha;

  if ( !fp || !multiAtoms ) {
    XLALPrintError ( "%s: invalid NULL input.\n", fn );
    XLAL_ERROR (fn, XLAL_EINVAL );
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
		    thisAtom->Fa_alpha.re, thisAtom->Fa_alpha.im,
		    thisAtom->Fb_alpha.re, thisAtom->Fb_alpha.im
		    );
	} /* for alpha < numSFTs */
    } /* for X < numDet */

  return XLAL_SUCCESS;

} /* write_MultiFstatAtoms_to_fp() */

