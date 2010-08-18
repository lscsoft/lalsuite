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
#define LAL_INT4_MAX 2147483647

/* macro to map indices {m,n} over {t0, tau} space into 1-dimensional array index, with tau-index n faster varying */
#define IND_MN(m,n) ( (m) * N_tauRange + (n) )

/* ---------- internal prototypes ---------- */
int compareAtoms(const void *in1, const void *in2);

/* empty struct initializers */
const TransientCandidate_t empty_TransientCandidate;
const transientWindow_t empty_transientWindow;
const transientWindowRange_t empty_transientWindowRange;

/* ==================== function definitions ==================== */

/** Helper-function to determine the total timespan of
 * a transient CW window, ie. the earliest and latest timestamps
 * of non-zero window function.
 */
int
XLALGetTransientWindowTimespan ( UINT4 *t0,				/**< [out] window start-time */
                                 UINT4 *t1,				/**< [out] window end-time */
                                 transientWindow_t transientWindow	/**< [in] window-parameters */
                                 )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !t0 || !t1 ) {
    XLALPrintError ("%s: invalid NULL input 't0=%p', 't1=%p'\n", fn, t0, t1 );
    XLAL_ERROR ( fn, XLAL_EINVAL );
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
      XLAL_ERROR_REAL8 ( fn, XLAL_EINVAL );

    } /* switch window-type */

  return XLAL_SUCCESS;

} /* XLALGetTransientWindowTimespan() */


/** apply a "transient CW window" described by TransientWindowParams to the given
 * timeseries
 */
int
XLALApplyTransientWindow ( REAL4TimeSeries *series,		/**< input timeseries to apply window to */
                           transientWindow_t transientWindow	/**< transient-CW window to apply */
                           )
{
  const CHAR *fn = __func__;

  /* check input consistency */
  if ( !series || !series->data ){
    XLALPrintError ("%s: Illegal NULL in input timeseries!\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* special time-saving break-condition: do nothing for window=none */
  if ( transientWindow.type == TRANSIENT_NONE )
    return XLAL_SUCCESS;

  /* deal with non-trivial windows */
  REAL8 ts_t0 = XLALGPSGetREAL8 ( &series->epoch );
  REAL8 ts_dt = series->deltaT;
  UINT4 ts_length = series->data->length;

  UINT4 t0, t1;
  if ( XLALGetTransientWindowTimespan ( &t0, &t1, transientWindow ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALGetTransientWindowTimespan() failed.\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EFUNC );
  }

  UINT4 i;
  switch ( transientWindow.type )
    {
    case TRANSIENT_RECTANGULAR:
      for ( i = 0; i < ts_length; i ++ )
        {
          UINT4 ti = (UINT4) round ( ts_t0 + i * ts_dt );
          REAL8 win = XLALGetRectangularTransientWindowValue ( ti, t0, t1 );
          series->data->data[i] *= win;
        } /* for i < length */
      break;

    case TRANSIENT_EXPONENTIAL:
      for ( i = 0; i < ts_length; i ++ )
        {
          UINT4 ti = (UINT4) round ( ts_t0 + i * ts_dt );
          REAL8 win = XLALGetExponentialTransientWindowValue ( ti, t0, t1, transientWindow.tau );
          series->data->data[i] *= win;
        } /* for i < length */
      break;

    default:
      XLALPrintError ("%s: invalid transient window type %d not in [%d, %d].\n",
                      fn, transientWindow.type, TRANSIENT_NONE, TRANSIENT_LAST -1 );
      XLAL_ERROR ( fn, XLAL_EINVAL );
      break;

    } /* switch (window.type) */

  return XLAL_SUCCESS;

} /* XLALApplyTransientWindow() */


/** apply transient window to give multi noise-weights, associated with given
 * multi timestamps
 */
int
XLALApplyTransientWindow2NoiseWeights ( MultiNoiseWeights *multiNoiseWeights,	/**< [in/out] noise weights to apply transient window to */
                                        const MultiLIGOTimeGPSVector *multiTS,	/**< [in] associated timestamps of noise-weights */
                                        transientWindow_t transientWindow	/**< [in] transient window parameters */
                                        )
{
  static const char *fn = __func__;

  UINT4 numIFOs, X;
  UINT4 numTS, i;

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

  /* special time-saving break-condition: do nothing for window=none */
  if ( transientWindow.type == TRANSIENT_NONE )
    return XLAL_SUCCESS;

  /* deal with non-trivial windows */
  UINT4 t0, t1;
  if ( XLALGetTransientWindowTimespan ( &t0, &t1, transientWindow ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALGetTransientWindowTimespan() failed.\n", fn );
    XLAL_ERROR_REAL8 ( fn, XLAL_EFUNC );
  }

  /* loop over all detectors X */
  for ( X = 0; X < numIFOs; X ++ )
    {
      numTS = multiNoiseWeights->data[X]->length;

      if ( multiTS->data[X]->length != numTS ) {
        XLALPrintError ("%s: inconsistent number of timesteps 'multiNoiseWeights[%d]' (%d) and 'multiTS[%d]' (%d).\n", fn, X, numTS, X, multiTS->data[X]->length );
        XLAL_ERROR ( fn, XLAL_EINVAL );
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
                          fn, transientWindow.type, TRANSIENT_NONE, TRANSIENT_LAST -1 );
          XLAL_ERROR ( fn, XLAL_EINVAL );
          break;

        } /* switch (window.type) */

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
  const CHAR *fn = __func__;
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
 *
 * Note2: if window->type == none, we compute a single rectangular window covering all the data.
 */
int
XLALComputeTransientBstat ( TransientCandidate_t *cand, 		/**< [out] transient candidate info */
                            const MultiFstatAtomVector *multiFstatAtoms,/**< [in] multi-IFO F-statistic atoms */
                            transientWindowRange_t windowRange )	/**< [in] type and parameters specifying transient window range to search */
{
  const char *fn = __func__;

  /* initialize empty return, in case sth goes wrong */
  TransientCandidate_t ret = empty_TransientCandidate;
  (*cand) = ret;

  /* check input consistency */
  if ( !multiFstatAtoms || !multiFstatAtoms->data || !multiFstatAtoms->data[0]) {
    XLALPrintError ("%s: invalid NULL input.\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  if ( windowRange.type >= TRANSIENT_LAST ) {
    XLALPrintError ("%s: unknown window-type (%d) passes as input. Allowed are [0,%d].\n", fn, windowRange.type, TRANSIENT_LAST-1);
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* combine all multi-atoms into a single atoms-vector with *unique* timestamps */
  FstatAtomVector *atoms;
  UINT4 TAtom = multiFstatAtoms->data[0]->TAtom;
  if ( (atoms = XLALmergeMultiFstatAtomsBinned ( multiFstatAtoms, TAtom )) == NULL ) {
    XLALPrintError ("%s: XLALmergeMultiFstatAtomsSorted() failed with code %d\n", fn, xlalErrno );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }
  UINT4 numAtoms = atoms->length;
  /* actual data spans [t0_data, t0_data + numAtoms * TAtom] in steps of TAtom */
  UINT4 t0_data = atoms->data[0].timestamp;
  UINT4 t1_data = atoms->data[numAtoms-1].timestamp + TAtom;

  /* special treatment of window_type = none ==> replace by rectangular window spanning all the data */
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

  /* It is often numerically impossible to compute e^F and sum these values, because of range-overflow
   * instead we first determine max{F_mn}, then compute the logB = log ( e^Fmax * sum_{mn} 1/D_mn * e^{Fmn - Fmax} )
   * which is logB = Fmax + log( sum_{mn} e^FReg_mn ), where FReg_mn = -log(D_mn) + Fmn - Fmax.
   * This avoids numerical problems.
   *
   * As we don't know Fmax before having computed the full matrix F_mn, we keep the full array of
   * 'regularized' F-stats FReg_mn over the field of {t0, tau} values in steps of dt0 x dtau.
   *
   * NOTE2: indices {i,j} enumerate *actual* atoms and their timestamps t_i, while the
   * indices {m,n} enumerate the full grid of values in [t0_min, t0_max]x[Tcoh_min, Tcoh_max] in
   * steps of deltaT. This allows us to deal with gaps in the data in a transparent way.
   *
   * NOTE3: we operate on the 'binned' atoms returned from XLALmergeMultiFstatAtomsBinned(),
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

  REAL8 *FReg_mn;	/* 2D matrix {m x n} of FReg values, will be initialized to zeros ! */
  UINT4 N_grid2D   = N_t0Range * N_tauRange;
  if ( ( FReg_mn = XLALCalloc ( N_grid2D, sizeof(REAL8) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( %d, sizeof(REAL8)\n", fn, N_grid2D );
    XLAL_ERROR ( fn, XLAL_ENOMEM );
  }

  ret.maxFstat = 0;	// keep track of loudest 2F-value over t0Band x tauBand space
  UINT4 m, n;

  transientWindow_t window;
  window.type = windowRange.type;

  /* ----- OUTER loop over start-times [t0,t0+t0Band] ---------- */
  for ( m = 0; m < N_t0Range; m ++ ) /* m enumerates 'binned' t0 start-time indices  */
    {
      /* compute Fstat-atom index i_t0 in [0, numAtoms) */
      window.t0 = windowRange.t0 + m * windowRange.dt0;
      INT4 i_tmp = (INT4)round ( 1.0 * ( window.t0 - t0_data ) / TAtom );
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
          window.tau = windowRange.tau + n * windowRange.dtau;

          /* get end-time t1 of this transient-window search */
          UINT4 t0, t1;
          if ( XLALGetTransientWindowTimespan ( &t0, &t1, window ) != XLAL_SUCCESS ) {
            XLALPrintError ("%s: XLALGetTransientWindowTimespan() failed.\n", fn );
            XLAL_ERROR ( fn, XLAL_EFUNC );
          }

          /* compute window end-time Fstat-atom index i_t1 in [0, numAtoms) */
          i_tmp = (INT4) round ( 1.0 * ( t1 - t0_data ) / TAtom ) - 1;
          if ( i_tmp < 0 ) i_tmp = 0;
          UINT4 i_t1 = (UINT4)i_tmp;
          if ( i_t1 >= numAtoms ) i_t1 = numAtoms - 1;

          /* protection against degenerate 1-atom case: (this implies D=0 and therefore F->inf) */
          if ( i_t1 == i_t0 ) {
            XLALPrintError ("%s: encountered a single-atom Fstat-calculation. This is degenerate and cannot be computed!\n", fn );
            XLALPrintError ("Window-values m=%d (t0=%d=t0_data + %d), n=%d (tau=%d) ==> t1_data - t0 = %d\n",
                            m, window.t0, i_t0 * TAtom, n, window.tau, t1_data - window.t0 );
            XLALPrintError ("The most likely cause is that your t0-range covered all of your data: t0 must stay away *at least* 2*TAtom from the end of the data!\n");
            XLAL_ERROR ( fn, XLAL_EDOM );
          }

          REAL8 Dd, twoF;

          /* now we have two valid atoms-indices [i_t0, i_t1] spanning our Fstat-window to sum over,
           * using weights according to the window-type
           */
          switch ( windowRange.type )
            {
            case TRANSIENT_RECTANGULAR:
              /* special optimiziation in the rectangular-window case: just add on to previous tau values
               * ie re-use the sum over [i_t0, i_t1_last] from the pevious tau-loop iteration
               */
              for ( UINT4 i = i_t1_last; i <= i_t1; i ++ )
                {
                  FstatAtom *thisAtom_i = &atoms->data[i];

                  /* now add on top of previous values, summed from [i_t0, i_t1_last] */
                  Ad += thisAtom_i->a2_alpha;
                  Bd += thisAtom_i->b2_alpha;
                  Cd += thisAtom_i->ab_alpha;

                  Fa_re += thisAtom_i->Fa_alpha.re;
                  Fa_im += thisAtom_i->Fa_alpha.im;

                  Fb_re += thisAtom_i->Fb_alpha.re;
                  Fb_im += thisAtom_i->Fb_alpha.im;

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

                  REAL8 win_i = XLALGetExponentialTransientWindowValue ( t_i, t0, t1, window.tau );
                  REAL8 win2_i = win_i * win_i;

                  Ad += thisAtom_i->a2_alpha * win2_i;
                  Bd += thisAtom_i->b2_alpha * win2_i;
                  Cd += thisAtom_i->ab_alpha * win2_i;

                  Fa_re += thisAtom_i->Fa_alpha.re * win_i;
                  Fa_im += thisAtom_i->Fa_alpha.im * win_i;

                  Fb_re += thisAtom_i->Fb_alpha.re * win_i;
                  Fb_im += thisAtom_i->Fb_alpha.im * win_i;

                } /* for i in [i_t0, i_t1] */
              break;

            default:
              XLALPrintError ("%s: invalid transient window type %d not in [%d, %d].\n",
                              fn, windowRange.type, TRANSIENT_NONE, TRANSIENT_LAST -1 );
              XLAL_ERROR ( fn, XLAL_EINVAL );
              break;

            } /* switch window.type */


          /* generic F-stat calculation from A,B,C, Fa, Fb */
          Dd = Ad * Bd - Cd * Cd;

          twoF = (2.0 / Dd) * ( Bd * (SQ(Fa_re) + SQ(Fa_im) ) + Ad * ( SQ(Fb_re) + SQ(Fb_im) )
                                - 2.0 * Cd *( Fa_re * Fb_re + Fa_im * Fb_im )
                                );

          /* keep track of loudest F-stat value encountered over the m x n matrix */
          if ( twoF > ret.maxFstat )
            {
              ret.maxFstat = twoF;
              ret.t0offs_maxF  = window.t0 - windowRange.t0;	/* start-time offset from earliest t0 in window-range*/
              ret.tau_maxF = window.tau;
            }

          /* compute 'regularized' F-stat: log ( 1/D * e^F ) = -logD + F */
          REAL8 FReg = - log( Dd ) + 0.5 * twoF;
          /* and store this in Fstat-matrix as element {m,n} */
          FReg_mn[IND_MN(m,n)] = FReg;

        } /* for n in n[tau] : n[tau+tauBand] */

    } /* for m in m[t0] : m[t0+t0Band] */

  /* now step through FReg_mn array subtract maxFstat and sum e^{FReg - Fmax}*/
  REAL8 sum_eB = 0;

  for ( m=0; m < N_t0Range; m ++ )
    {
      for ( n=0; n < N_tauRange; n ++ )
        {
          REAL8 FReg = FReg_mn [ IND_MN(m,n) ];

          sum_eB += exp ( FReg - 0.5 * ret.maxFstat );

        } /* for n < N_tauRange */

    } /* for m < N_t0Range */

  /* combine this to final log(Bstat) result with proper normalization (assuming hmaxhat=1) : */

  REAL8 logBhat = 0.5 * ret.maxFstat + log ( windowRange.dt0 * windowRange.dtau ) + log ( sum_eB );	/* unnormalized Bhat */
  /* final normalized Bayes factor, assuming hmaxhat=1 */
  /* NOTE: correct for different hmaxhat by adding "- 4 * log(hmaxhat)" to this */

  REAL8 normBh = 70.0 / ( N_t0Range * N_tauRange * TAtom * TAtom );
  ret.logBstat = log ( normBh ) +  logBhat;	/* - 4.0 * log ( hmaxhat ) */

  /* free mem */
  XLALDestroyFstatAtomVector ( atoms );
  XLALFree ( FReg_mn );

  /* return */
  (*cand) = ret;
  return XLAL_SUCCESS;

} /* XLALComputeTransientBstat() */


/** Combine N Fstat-atoms vectors into a single 'canonical' binned and ordered atoms-vector.
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
  const char *fn = __func__;

  if ( !multiAtoms || !multiAtoms->length || !multiAtoms->data[0] || (deltaT==0) ) {
    XLALPrintError ("%s: invalid NULL input or deltaT=0.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }

  UINT4 numDet = multiAtoms->length;
  UINT4 X;
  UINT4 TAtom = multiAtoms->data[0]->TAtom;

  /* check consistency of time-step lengths between different IFOs */
  for ( X=0; X < numDet; X ++ ) {
    if ( multiAtoms->data[X]->TAtom != TAtom ) {
      XLALPrintError ("%s: Invalid input, atoms baseline TAtom=%d must be identical for all multiFstatAtomVectors (IFO=%d: TAtom=%d)\n",
                      fn, TAtom, X, multiAtoms->data[X]->TAtom );
      XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
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
    XLALPrintError ("%s: failed to XLALCreateFstatAtomVector ( %d )\n", fn, NBinnedAtoms );
    XLAL_ERROR_NULL ( fn, XLAL_ENOMEM );
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
          destAtom->Fa_alpha.re += atom_X_i->Fa_alpha.re;
          destAtom->Fa_alpha.im += atom_X_i->Fa_alpha.im;
          destAtom->Fb_alpha.re += atom_X_i->Fb_alpha.re;
          destAtom->Fb_alpha.im += atom_X_i->Fb_alpha.im;

        } /* for i < numAtomsX */
    } /* for X < numDet */

  return atomsOut;

} /* XLALmergeMultiFstatAtomsBinned() */

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
    fprintf (fp, "\n%%%%        fkdot[0]         Alpha[rad]         Delta[rad]  fkdot[1] fkdot[2] fkdot[3]   twoFtotal  t0offs_maxF[d] tau_maxF[d]      maxFstat       logBstat\n");
  else
    fprintf (fp, "%18.16g %18.16g %18.16g %8.6g %8.5g %8.5g  %11.9g        %7.5f      %7.5f   %11.9g    %11.9g\n",
             thisCand->doppler.fkdot[0], thisCand->doppler.Alpha, thisCand->doppler.Delta,
             thisCand->doppler.fkdot[1], thisCand->doppler.fkdot[2], thisCand->doppler.fkdot[3],
             thisCand->twoFtotal,
             1.0 * thisCand->t0offs_maxF / DAY24, 1.0 * thisCand->tau_maxF / DAY24, thisCand->maxFstat,
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

