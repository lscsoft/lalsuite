/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton, John Whelan
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

#ifndef _INJECT_H
#define _INJECT_H

#include <lal/LALStdlib.h>
#include <lal/Random.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup Inject_h
 * \author Creighton, T. D.
 *
 * \brief Provides routines to inject a signal into detector output.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/Inject.h>
 * \endcode
 *
 * This header provides simple routines to inject a signal, stored as a
 * floating-point time series, into an integer time series that
 * represents the ADC output of a detector channel.
 *
 * The basic concept at work here is that of \e dithering.  That is,
 * to add a real signal \f$x(t_k)\f$ to the integer output \f$n(t_k)\f$ of an
 * ADC, we cannot simply compute \f$x+n\f$ at each time \f$t_k\f$ and round to
 * the nearest integer.  To see this, consider injecting a sinusoid with
 * an amplitude less than half the ADC resolution, \f$|x|<0.5\f$.  Then
 * adding \f$x+n\f$ and rounding will always give back \f$n\f$, and the signal
 * injection will have no effect on the output.
 *
 * Instead, what we would like to do is to add \f$x\f$ to the ADC
 * \e input before rounding occurs, and then round to an integer.
 * That is, the ADC input was actually \f$n+d\f$, where \f$d\in[-0.5,0.5)\f$; we
 * want to compute the rounded value of \f$n+d+x\f$, which may occasionally
 * be different from \f$n\f$.  In principle, a Fourier transforms can detect
 * a sinusoidal signal with an output much smaller than the ADC
 * resolution: by integrating enough data, one can eventually detect a
 * statistically significant phase correlation between the occasional
 * increments and decrements in \f$n\f$.
 *
 * Of course given the output time series \f$n(t_k)\f$ we can only guess at
 * the input series \f$n(t_k)+d(t_k)\f$ that produced it.  The simplest guess
 * is to assume that each \f$d(t_k)\f$ is an independent random variable with
 * a flat distribution over the range \f$[-0.5,0.5)\f$.  This is a reasonable
 * guess to make if the root mean square variation between succesive
 * output values \f$n\f$ is a few or more ADC counts; i.e. if the
 * dimensionless power spectral density \f$\sqrt{fS(f)}\f$ has a value of a
 * few or more around the sampling freqeuncy \f$f\f$.  This is almost always
 * true of any detector designed to work at or near its noise limit: the
 * input to the ADC will first be whitened so that \f$\sqrt{fS(f)}\f$ is
 * nearly flat, and then amplified so that \f$\sqrt{fS(f)}\f$ is on the order
 * of several (or more) ADC counts.
 *
 * In the routines covered by this header we will take it for granted
 * that the above is a reasonable approximation, and will not check for
 * it.  We will further assume that the signal to be injected has already
 * been subjected to the same whitening and amplification, so that the
 * units of \f$x(t_k)\f$ are normalized ADC counts (although it is still a
 * real number, not an integer).
 *
 * The dithering routines should be used whenever one is injecting a
 * signal into a time series representing raw digitized data.  In some
 * data storage specifications, ADC output is not stored as an integer,
 * but as a floating-point number representing an integer value.  Such
 * data must be cast to integers before being passed to the digitizing
 * routines.
 *
 * This header also provides even simpler routines for injecting a signal
 * into floating-point data, without dithering.  These should only be
 * used when the data is genuinely continuous in character.  This can
 * include data derived by applying floating-point operations on ADC
 * channels (e.g.\ digital filters, linear combinations of channels,
 * etc.), but not data that simply represents ADC output in
 * floating-point format.  The assumption here is that the numerical
 * post-processing of the ADC data completely masks any statistical
 * signiatures of the digitization.
 *
 * @{
 * \defgroup InjectVector_c Module InjectVector.c
 * \defgroup InjectTimeSeries_c Module InjectTimeSeries.c
 * @}
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define INJECTH_ENUL  1	/**< Unexpected null pointer in arguments. */
#define INJECTH_EBAD  2 /**< A sampling interval is (effectively) zero */
#define INJECTH_EUNIT 3	/**< Input or output is not in units of ADC counts */
/*@}*/
/*@}*/

#define INJECTH_MSGENUL  "Unexpected null pointer in arguments"
#define INJECTH_MSGEBAD  "A sampling interval is (effectively) zero"
#define INJECTH_MSGEUNIT "Input or output is not in units of ADC counts"

/* Function prototypes. */

void
LALSI2InjectVector( LALStatus    *,
		    INT2Vector   *output,
		    REAL4Vector  *signalvec,
		    RandomParams *params );

void
LALSSInjectVector( LALStatus    *,
		   REAL4Vector  *output,
		   REAL4Vector  *signalvec );

void
LALSI2InjectTimeSeries( LALStatus       *,
			INT2TimeSeries  *output,
			REAL4TimeSeries *signalvec,
			RandomParams    *params );

void
LALSSInjectTimeSeries( LALStatus       *,
		       REAL4TimeSeries *output,
		       REAL4TimeSeries *signalvec );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _INJECT_H */
