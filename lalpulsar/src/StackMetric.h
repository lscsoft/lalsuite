/*
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix, Teviet Creighton
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

#ifndef _STACKMETRIC_H
#define _STACKMETRIC_H

#include <lal/LALStdlib.h>
#include <lal/PulsarTimes.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup StackMetric_h Header StackMetric.h
 * \author Creighton, T. D.
 * \date 2000 -- 2003
 * \ingroup pkg_pulsarMetric
 * \brief Provides routines to compute parameter-space metrics for coherent or stacked pulsar searches.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/StackMetric.h>
 * \endcode
 *
 * This module covers routines that determine the metric
 * coefficients for the mismatch function (ambiguity function) on the
 * parameter space for a pulsar search.  The assumed search method is
 * stacking one or more Fourier power spectra, after some suitable
 * demodulation.
 *
 * The method for determining the parameter metric is discussed in detail
 * in Sec.\ II of \ref Brady_P2000; we present the key results here in
 * brief.  We assume that a model waveform in our search is described by
 * an overall frequency scale \f$f_0\f$, and by some modulation about that
 * frequency described by ``shape'' parameters
 * \f$\vec{\lambda}=(\lambda^1,\ldots,\lambda^n)\f$, such that the
 * parameterized phase of the waveform is \f$\phi[t;\mathbf{\lambda}] = 2\pi
 * f_0\tau[t;\vec{\lambda}]\f$.  Here \f$\mathbf{\lambda} = (\lambda^0,\vec{\lambda})
 * = (f_0,\lambda^1,\ldots,\lambda^n)\f$ represents the total set of
 * parameters we must search over, and \f$\tau[t;\vec{\lambda}]\f$ is a
 * canonical time coordinate describing the shape of the waveform.
 *
 * A (local) maximum in detected power \f$P\f$ occurs if a signal is filtered
 * (or demodulated in time and the Fourier spectrum is sampled) using a
 * phase model that matches the true phase of the signal.  If the
 * parameters \f$\mathbf{\lambda}\f$ do not match the true parameters of the
 * signal, then the detected power will be degraded.  The fractional
 * power loss \f$\Delta P/P\f$ thus has a (local) minimum of 0 for matched
 * parameters and increases for mismatched parameters; it can be thought
 * of as describing a distance between the two (nearby) parameter sets.
 * The \em metric of this distance measure is simply the set of
 * quadratic coefficients of the Taylor expansion of \f$\Delta P/P\f$ about
 * its minimum.
 *
 * Clearly the power will degrade rapidly with variation in some
 * parameter \f$\lambda^\alpha\f$ if the phase function \f$\phi\f$ depends
 * strongly on that parameter.  It turns out that if the detected power
 * is computed from a coherent power spectrum of a time interval \f$\Delta t\f$,
 * then the metric components are given simply by the covariances of
 * the phase derivatives \f$\partial\phi/\partial\lambda^\alpha\f$ over the
 * time interval:
 * \f{equation}{
 * \label{eq_gab_phi}
 * g_{\alpha\beta}(\mathbf{\lambda}) =
 * \left\langle
 * \frac{\partial\phi[t;\mathbf{\lambda}]}{\partial\lambda^\alpha}
 * \frac{\partial\phi[t;\mathbf{\lambda}]}{\partial\lambda^\beta}
 * \right\rangle
 * -
 * \left\langle
 * \frac{\partial\phi[t;\mathbf{\lambda}]}{\partial\lambda^\alpha}
 * \right\rangle
 * \left\langle
 * \frac{\partial\phi[t;\mathbf{\lambda}]}{\partial\lambda^\beta}
 * \right\rangle \; ,
 * \f}
 * where \f$\langle\ldots\rangle\f$ denotes a time average over the interval
 * \f$\Delta t\f$, and \f$\alpha\f$ and \f$\beta\f$ are indecies running from 0 to
 * \f$n\f$.  The partial derivatives are evaluated at the point \f$\mathbf{\lambda}\f$
 * in parameter space.  If instead the detected power is computed from
 * the sum of several power spectra computed from separate time intervals
 * (of the same length), then the overall metric is the \em average of
 * the metrics from each time interval.
 *
 * When power spectra are computed using fast Fourier transforms, the
 * entire frequency band from DC to Nyquist is computed at once; one then
 * scans all frequencies for significant peaks.  In this case one is
 * concerned with how the peak power (maximized over frequency) is
 * reduced by mismatch in the remaining ``shape'' parameters
 * \f$\vec{\lambda}\f$.  This is given by the the \em projected metric
 * \f$\gamma_{ij}(\vec{\lambda})\f$, where \f$i\f$ and \f$j\f$ run from 1 to \f$n\f$:
 * \f{equation}{
 * \label{eq_gij_gab}
 * \gamma_{ij}(\vec{\lambda}) = \left[g_{ij}-\frac{g_{0i}g_{0j}}{g_{00}}
 * \right]_{\lambda^0=f_\mathrm{max}} \; .
 * \f}
 * Here \f$f_\mathrm{max}\f$ is the highest-frequency signal expected to be
 * present, which ensures that, for lower-frequency signals,
 * \f$\gamma_{ij}\f$ will \em overestimate the detection scheme's
 * sensitivity to the ``shape'' parameters.
 */
/*@{*/

/** \name Error conditions */
/*@{*/
#define STACKMETRICH_ENUL 1
#define STACKMETRICH_EBAD 2

#define STACKMETRICH_MSGENUL "Null pointer"
#define STACKMETRICH_MSGEBAD "Bad parameter values"
/*@}*/

/**
 * \brief This structure stores and passes parameters for computing a parameter-space metric.
 * It points to the canonical time function used
 * to compute the metric and to the parameters required by this function.
 * In addition, this structure must indicate the timespan over which the
 * timing differences accumulate, and whether this accumulation is
 * coherent or divided into stacks which are summed in power.
 */
typedef struct tagMetricParamStruc{
  void (*dtCanon)(LALStatus *, REAL8Vector *, REAL8Vector *, PulsarTimesParamStruc * ); /**< The function to compute the canonical
											  * time coordinate and its derivatives. */
  PulsarTimesParamStruc *constants; 	/**< The constant parameters used by *dt(). */
  REAL8 start;  			/**< Start time of search, measured relative to constants->epoch. */
  REAL8 deltaT;				/**< Length of each stack, in s. */
  UINT4 n;				/**< Number of stacks. */
  BOOLEAN errors;			/**< Whether to estimate errors in metric. */
} MetricParamStruc;


void
LALCoherentMetric( LALStatus *,
		   REAL8Vector *metric,
		   REAL8Vector *lambda,
		   MetricParamStruc *params );

void
LALStackMetric( LALStatus *,
		REAL8Vector *metric,
		REAL8Vector *lambda,
		MetricParamStruc *params );

void
LALProjectMetric( LALStatus *, REAL8Vector *metric, BOOLEAN errors );

/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _STACKMETRIC_H */
