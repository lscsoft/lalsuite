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

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \defgroup Inject_h Header Inject.h
 * \ingroup lal_inject
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
 * This header also provides simple routines for injecting a signal
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
LALSSInjectTimeSeries( LALStatus       *,
		       REAL4TimeSeries *output,
		       REAL4TimeSeries *signalvec );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _INJECT_H */
