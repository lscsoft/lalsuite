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
XLALoutputMultiFstatAtoms ( FILE *fp, MultiFstatAtoms *multiAtoms )
{
  const char *fn = "XLALoutputMultiFstatAtoms()";
  UINT4 X, alpha;

  if ( !fp || !multiAtoms )
    XLAL_ERROR (fn, XLAL_EINVAL );

  fprintf ( fp, "%% GPS                a(t_i)     b(t_i)            Fa(t_i)                 Fb(t_i)\n");

  for ( X=0; X < multiAtoms->length; X++ )
    {
      FstatAtoms *thisAtom = multiAtoms->data[X];
      for ( alpha=0; alpha < multiAtoms->data[X]->length; alpha ++ )
	{
	  fprintf ( fp, "%f   % f  % f     % f  % f     % f  % f\n",
		    XLALGPSGetREAL8( &thisAtom->timestamps[alpha] ),
		    thisAtom->a_alpha[alpha],
		    thisAtom->b_alpha[alpha],
		    thisAtom->Fa_alpha[alpha].re, thisAtom->Fa_alpha[alpha].im,
		    thisAtom->Fb_alpha[alpha].re, thisAtom->Fb_alpha[alpha].im
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
XLALComputeTransientBstat ( const MultiFstatAtoms *multiFstatAtoms,	/**< [in] multi-IFO F-statistic atoms */
                            transientWindowRange_t windowRange )	/**< [in] type and parameters specifying transient window range to search */
{
  const char *fn = __func__;

  REAL8 tau0 = windowRange.tau_min; // smallest search duration-window in seconds
  REAL8 tau1 = windowRange.tau_max; // longest search duration-window in seconds

  REAL8 t0 = windowRange.t0_min;	// earliest GPS start-time
  REAL8 t1 = windowRange.t0_max;	// latest GPS start-time

  // check input argument consistency
  if ( t1 < t0 )
    {
      XLALPrintError ("%s: invalid input arguments t0 (=%d), t1 (=%d): must t1>t0 \n", fn, t0, t1 );
      XLAL_ERROR_REAL8 ( fn, XLAL_EDOM );
    }

  if ( tau1 < tau0 )
    {
      XLALPrintError ("%s: invalid input arguments tau0 (=%d), tau1 (=%d): must tau1>tau0 \n", fn, tau0, tau1 );
      XLAL_ERROR_REAL8 ( fn, XLAL_EDOM );
    }


  UINT4 t_i;        // t index (t-summation)
  REAL8 tau_i;      // tau index (tau-summation)
  REAL8 logBAYES = 0;  // return value of function

  // Define and initiate F-Stat variables
  REAL8 A = 0;
  REAL8 B = 0;
  REAL8 C = 0;
  REAL8 D = 0;
  REAL8 Ds = 0;
  COMPLEX8 FA = {0,0};
  COMPLEX8 FB = {0,0};
  REAL8 F = 0;

  for (tau_i=tau0; tau_i <= tau1; tau_i+=1800) // FIXME use SFT duration
    {

      for (t_i=t0; t_i<=t1; t_i+=1800)
        {
          A = 0;
          B = 0;
          C = 0;
          D = 0;
          FA.re = 0;
          FA.im = 0;
          FB.re = 0;
          FB.im = 0;
          for ( UINT4 X=0; X < multiFstatAtoms->length; X++ ) // Loop over detectors
            {
              // Per IFO data ("atoms")
              FstatAtoms *atoms_X = multiFstatAtoms->data[X];
              UINT4 Natoms = atoms_X->length;
              LIGOTimeGPS *t = atoms_X->timestamps;
              REAL8 *a = atoms_X->a_alpha;
              REAL8 *b = atoms_X->b_alpha;
              COMPLEX8 *Fa = atoms_X->Fa_alpha;
              COMPLEX8 *Fb = atoms_X->Fb_alpha;

              UINT4 j = 0;
              while ( (j < Natoms) && (t[j].gpsSeconds < t_i) )
                j ++;


              //for (j=t_i; j<t_i+tau_i; j++) t0-tau summation
              while ( (j < Natoms) && ( t[j].gpsSeconds <= t_i + tau_i) )
                {
                  A += SQ(a[j]);
                  B += SQ(b[j]);
                  C += (a[j]*b[j]);

                  FA.re += Fa[j].re;
                  FA.im += Fa[j].im;
                  FB.re += Fb[j].re;
                  FB.im += Fb[j].im;

                  j++;
                }


            } // for X < numDet
          //printf("j=%d\n",j);
          D = A*B - SQ(C);
          F += (1.0/D)*( B*(SQ(FA.re)+SQ(FA.im)) +
                         A*(SQ(FB.re)+SQ(FB.im)) -
                         2*C*(FA.re*FB.re+FA.im*FB.im) ); /* Compute summed F-Statistic */
          Ds += D;
        } // for t0 in t0Range
    } // for tau in tauRange
  logBAYES = F - log(Ds);
  printf("t_i=%d tau_i=%f \n",t_i-1800,(tau_i-1800)/(3600*24));
  printf("A=%f \n",A);
  printf("B=%f \n",B);
  printf("C=%f \n",C);
  printf("D=%f \n",D);
  printf("sum(D)=%f \n",Ds);
  printf("sum(2F)=%f \n",2*F);
  printf("log(Bayes)=%f \n",logBAYES);

  return logBAYES;

} /* XLALComputeTransientBstat() */
