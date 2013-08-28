/*
*  Copyright (C) 2007 David Churches, Duncan Brown, Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
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

/**
 * \author Sathyaprakash, B. S.
 * \file
 * \ingroup LALInspiral_h
 *
 * \brief Module to calculate the number of data points (to the nearest power of 2)
 * needed to store a waveform.
 *
 * \heading{Prototypes}
 *
 * <tt>LALInspiralWaveLength()</tt>
 * <ul>
 * <li> \c length: output, number of bins to the nearest power of 2
 * greater than the minimum length required to store a wave of parameters as in \c params.</li>
 * <li> \c params: input, parameters of the binary system.</li>
 * </ul>
 *
 * \heading{Description}
 *
 * This module first calls LALInspiralChooseModel(), which gives the length of the
 * waveform in seconds. That function returns an estimated waveform length. However, the
 * length might not be appropriate in some extreme cases (large masses and large lower
 * cut-off frequency). It is especially true in the EOB case. Therefore, we introduce
 * two constants namely LALINSPIRAL_LENGTHOVERESTIMATION (in percentage) which
 * overestimate the length of the waveform and LALINSPIRAL_MINIMALWAVELENGTH which is
 * the minimal waveform length in seconds. Multiplying this by the sampling rate
 * <tt>params.tSampling</tt> gives the minimum number of samples needed to hold the waveform.
 * To this are added the number of bins of leading and trailing zeroes requested by the user in
 * <tt>params.nStartPad</tt> and <tt>params.nEndPad.</tt> The resulting number is rounded to
 * an upward power of 2 and returned in \c length.
 *
 * \heading{Algorithm}
 *
 * \heading{Uses}
 * This function calls:
 * \code
 * LALInspiralSetup()
 * LALInspiralChooseModel()
 * \endcode
 *
 */

#include <lal/LALInspiral.h>
#include <lal/LALStdlib.h>

#define  LALINSPIRAL_LENGTHOVERESTIMATION  0.1       /* 10 % */
#define  LALINSPIRAL_MINIMALWAVELENGTH     0.03125   /* 64 bins with a 2048Hz sampling*/

void
LALInspiralWaveLength(
   LALStatus        *status,
   UINT4            *length,
   InspiralTemplate  params
   )
{

   INT4 ndx;
   REAL8 x;
   CHAR message[256];
   expnCoeffs ak;
   expnFunc func;

   INITSTATUS(status);
   ATTATCHSTATUSPTR(status);

   ASSERT (params.nStartPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params.nEndPad >= 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);
   ASSERT (params.tSampling > 0, status, LALINSPIRALH_ESIZE, LALINSPIRALH_MSGESIZE);

   LALInspiralSetup (status->statusPtr, &ak, &params);
   CHECKSTATUSPTR(status);

   LALInspiralChooseModel(status->statusPtr, &func, &ak, &params);
   /* not checkstatus before but we check if everything is fine,
      then we return the length otherwise length is zero*/
   if (status->statusPtr->statusCode == 0){
	   /*we add a minimal value and 10 % of overestimation */
      	   x	= (1.+ LALINSPIRAL_LENGTHOVERESTIMATION) * (ak.tn + LALINSPIRAL_MINIMALWAVELENGTH ) * params.tSampling
	     + params.nStartPad + params.nEndPad;
	   ndx 	= ceil(log10(x)/log10(2.));
      	   *length = pow(2, ndx);

	DETATCHSTATUSPTR(status);
	RETURN(status);
   }
   else {
	sprintf(message,
	       	"size is zero for the following waveform: totalMass = %f, fLower = %f, approximant = %d @ %fPN"
	       , params.mass1 + params.mass2, params.fLower, params.approximant, params.order/2.);
       	LALWarning(status, message);
       	*length = 0;

       DETATCHSTATUSPTR(status);
       RETURN(status);
     }
}
