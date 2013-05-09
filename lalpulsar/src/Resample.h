/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

#ifndef _RESAMPLE_H
#define _RESAMPLE_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
   \defgroup Resample_h Header Resample.h
   \ingroup pkg_pulsarCommon
   \author Creighton, T. D.
   \brief Provides routines for resampling time series according to a new canonical time coordinate.

\heading{Synopsis}
\code
#include <lal/Resample.h>
\endcode

One of the crucial problems in searching for
constant-frequency astrophysical signals is removing the effects of
Doppler modulation due to the Earth's motion.  This is normally
accomplished by constructing a canonical time coordinate \f$\tau\f$ of an
inertial frame (i.e.\ the <em>barycentred time</em>), and
decimating/resampling the data at fixed intervals in \f$\tau\f$.  The
reconstructed \f$\tau\f$ depends on the direction to the source relative
to the Earth's motion; in addition, slow intrinsic parameterized
modulations in the source frequency can also be corrected by this
coordinate transformation.

Most of the routines in this module assume that \f$\tau\f$ can be
piecewise expanded as a Taylor series in \f$t\f$.  That is, one defines a
set of fitting \e regions \f$T_i=[t_{\mathrm{bound}(i-1)},
t_{\mathrm{bound}(i)}]\f$, and a set of fitting \e points
\f$t_{(i)}\in T_i\f$.  In each region one then writes:
\anchor eq_tau
\f{equation}{
\label{eq_tau}
\tau(t) = \sum_{k=0} \frac{1}{k!}c_{k(i)}(t-t_{(i)})^k \; .
\f}
Since one is normally interested in tracking the difference
\f$\tau(t)-t\f$, one can also write the expansion as:
\anchor eq_delta-tau
\f{equation}{
\label{eq_delta-tau}
\tau(t)-t = \sum_{k=0} a_{k(i)}(t-t_{(i)})^k \; ,
\f}
where
\anchor eq_a_c
\f{eqnarray}{
a_{0(i)} & = & c_{0(i)}-t_{(i)}           \; , \nonumber\\
a_{1(i)} & = & c_{1(i)}-1                 \; , \nonumber\\
a_{k(i)} & = & c_{k(i)}/k! \; , \; k\geq2 \; . \nonumber
\label{eq_a_c}
\f}
These are the polynomial coefficients normally assumed in the modules
under this header.

The procedure for resampling according to \f$\tau\f$ is normally combined
with \e decimating the time series.  That is, one takes a time
series sampled at constant intervals \f$\Delta t\f$ in \f$t\f$, and samples it
at constant intervals \f$d\Delta t\f$ in \f$\tau\f$, where the
<em>decimation factor</em> \f$d\f$ is normally taken to be an integer
\f$\geq1\f$.  When \f$\tau\f$ and \f$t\f$ are drifting out of phase relatively
slowly, this means that most of the time every \f$d^\mathrm{th}\f$ sample
in the original time series becomes the next sample in the decimated
time series.  However, when \f$\tau\f$ and \f$t\f$ drift out of synch by an
amount \f$\pm\Delta t\f$, one can force the decimated time series to track
\f$\tau\f$ (rather than \f$t\f$) by sampling the \f$d\pm1^\mathrm{th}\f$ next
datum (rather than the \f$d^\mathrm{th}\f$).  If the drift is sufficiently
rapid or \f$d\f$ is sufficiently large, one may be forced to choose the
point \f$d\pm2\f$, \f$d\pm3\f$, etc.; the size of this adjustment is called
the correction \e shift.  The number of (resampled) time intervals
between one correction point and the next is called the correction
\e interval.

Unless otherwise specified, all time variables and parameters in the
functions under this header can be assumed to measure the detector
time coordinate \f$t\f$.  Canonical times are specified by giving the
difference \f$\tau-t\f$.

\heading{Caveat emptor:} The inclusion of this header and its
associated modules into LAL is provisional at this time.  The routines
and the test code appear to work, but a later standalone code,
operating on much larger datasets, appeared to encounter a memory
leak.  I have not yet determined whether this leak was in the
standalone code or in these LAL routines.

*/
/*@{*/

/** \name Error Codes */
/*@{*/
#define RESAMPLEH_ENUL    1
#define RESAMPLEH_EOUT    2
#define RESAMPLEH_EMEM    3
#define RESAMPLEH_EDTPOS  4
#define RESAMPLEH_ELENGTH 5
#define RESAMPLEH_ETIME   6

#define RESAMPLEH_MSGENUL    "Unexpected null pointer in arguments"
#define RESAMPLEH_MSGEOUT    "Output handle points to a non-null pointer"
#define RESAMPLEH_MSGEMEM    "Memory allocation error"
#define RESAMPLEH_MSGELENGTH "Vector lengths in polyco structure don't argree"
#define RESAMPLEH_MSGEDTPOS  "Sampling interval is not positive"
#define RESAMPLEH_MSGETIME   "Requested output time span extends beyond range of validity of input"
/*@}*/


/**
 * The rules for taking a time series \f$t\f$, sampled at constant intervals \f$\Delta t\f$, and resampling it at
 * constant intervals \f$d\Delta t\f$ in the canonical time coordinate \f$\tau\f$.
 */
typedef struct tagResampleRules {
  LIGOTimeGPS start; /**< The initial time for which the rules apply */
  LIGOTimeGPS stop;  /**< The final time for which the rules apply */
  INT4 length;       /**< The number of correction points, i.e.\ points where the resampling interval
                      * is adjusted from \f$d\Delta t\f$ to \f$(d\pm n)\Delta t\f$ */
  INT4 *interval;    /**< An array giving the number of resampled time intervals between correction points */
  INT2 *shift;       /**< An array giving the size of the correction shift (i.e.\ the number \f$n\f$ above) at each correction point */
  INT4 decimate;     /**< decimation factor \f$d\f$ */
  REAL8 deltaT;      /**< The sampling interval before decimation, n seconds */
  REAL8 startDiff;   /**< The difference \f$\tau-t\f$ at the time #start, in seconds */
  REAL8 stopDiff;    /**< The difference \f$\tau-t\f$ at the time #stop, in seconds */
} ResampleRules;


/** Parameters of the piecewise polynomial fit of \f$\tau-t\f$ as a function of \f$t\f$, see
 *  Eq.\eqref{eq_delta-tau} for notation.
 */
typedef struct tagPolycoStruc {
  REAL4 ra;  			/**< Right ascension angle of the source, in \e radians in the range \f$[0,2\pi)\f$ */
  REAL4 dec; 			/**< Declination angle of the source, in \e radians in the range \f$[-\pi/2,\pi/2]\f$ */
  REAL4Vector *spindown; 	/**< A vector \f$\vec\lambda=(\lambda_0,\ldots,\lambda_{n-1})\f$ of parameters
                                 * describing a slow intrisic frequency drift \f$f=f(t)\f$ of the source:
                                 * \f$\lambda_k=f^{-1}d^{k+1}f/dt^{k+1}\f$ at the time given by #start. */

  LIGOTimeGPS start;  		/**< The initial time over which the polynomial fit applies */

  REAL4Sequence *tBound; 	/**< The sequence of times \f$t_{\mathrm{bound}(i)}\f$ defining the endpoints
                                 * of the fitting regions, given in seconds after the time #start;
                                 * The first fitting region \f$i=0\f$ runs from #start to
                                 * #start + \f$t_{\mathrm{bound}(0)}\f$, the next from there to
                                 * #start + \f$t_{\mathrm{bound}(1)}\f$, and so on. */

  REAL4Sequence *t0;     	/**< The sequence of times \f$t_{(i)}\f$ in each fitting region at which the polynomial
                                 * fits are computed, given in seconds after the time #start. */

  REAL4VectorSequence *polyco;	/**< A sequence of vectors \f$\vec a_{(i)}=(a_{0(i)},a_{1(i)},\ldots)\f$ giving the
                                 * coefficients of the polynomial fit at each time \f$t_{(i)}\f$;
                                 * each element \f$a_{k(i)}\f$ has units of \f$\mathrm{s}^{1-k}\f$. */
} PolycoStruc;


/** Extra parameters required to construct a ResampleRules object from a PolycoStruc object.
 */
typedef struct tagResampleParamStruc{
  LIGOTimeGPS start;    /**< The initial time for which the resample rules will apply */
  LIGOTimeGPS stop;     /**< The final time for which the resample rules will apply */
  REAL8       deltaT;   /**< The sampling interval before decimation, in seconds */
  INT4        decimate; /**< The decimation factor */
} ResampleParamStruc;


/* Function prototypes. */
void
LALCreateResampleRules( LALStatus          *status,
			ResampleRules      **rules,
			PolycoStruc        *polyco,
			ResampleParamStruc *params );




void
LALDestroyResampleRules( LALStatus     *status,
			 ResampleRules **rules );




void
LALApplyResampleRules( LALStatus       *status,
		       REAL4TimeSeries *output,
		       REAL4TimeSeries *input,
		       ResampleRules   *rules );

/* Patrick: I don't know if you want ApplyResampleRules() to be a
   lower-level routine that operates on vectors, or whether time
   series are okay. */




void
LALPolycoToTimingDifference( LALStatus       *status,
			     REAL4TimeSeries *difference,
			     PolycoStruc     *polyco );




void
LALRulesToTimingDifference( LALStatus       *status,
			    REAL4TimeSeries *difference,
			    ResampleRules   *rules );



/*@}*/

#ifdef __cplusplus
}
#endif

#endif /* _RESAMPLE_H */
