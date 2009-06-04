
/*
 *  Copyright (C) 2007 Anand Sengupta
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

/*----------------------------------------------------------------------------
 *
 * File Name: BandPassInspiralTemplate.c
 *
 * Author: Anand Sengupta
 *
 * $Id$
 *
 *---------------------------------------------------------------------------*/

#if 0
<lalVerbatim file="BandPassInspiralTemplateCV">
Author: Anand Sengupta
$Id$
</lalVerbatim>
#endif


#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALDatatypes.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>

#include <lal/BandPassTimeSeries.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>

int XLALBandPassInspiralTemplate(
        REAL4Sequence  *sequence,
        REAL4          fLow,
        REAL4          fHigh,
        REAL4          fSampling
        )
{

    static const char *func = "XLALBandPassInspiralTemplate";

    REAL4TimeSeries   *signal = NULL;
    PassBandParamStruc signalHighpassParam;
    PassBandParamStruc signalLowpassParam;

    /* Check that the input makes sense */
    if ( !sequence || sequence->length <= 0 || !sequence->data ||
            fLow > fHigh || fLow <= 0.0 || fHigh >= fSampling/2.0 )
          XLAL_ERROR( func, XLAL_EFAULT );

    /* Initialize signal */
    signal = (REAL4TimeSeries *)
            LALCalloc(1, sizeof(REAL4TimeSeries));

    signal->deltaT = 1.0/fSampling;
    signal->data   = sequence;

    /* Now first high pass the time series beyond fSeismic */
    signalHighpassParam.nMax = 10;
    signalHighpassParam.f1 = -1.0;
    signalHighpassParam.a1 = -1.0;
    signalHighpassParam.f2 = fLow;
    signalHighpassParam.a2 = 0.98;

    /* Call the Butterworth routine and check its success */
    if ( XLALButterworthREAL4TimeSeries(
                signal, &signalHighpassParam ) != XLAL_SUCCESS )
    {
        if (signal) LALFree( signal );
        XLAL_ERROR(func, XLAL_EFUNC);
    }

    /* Now first low pass the time series below fFinal */
    signalLowpassParam.nMax = 10;
    signalLowpassParam.f1 =  fHigh;
    signalLowpassParam.a1 =  0.02;
    signalLowpassParam.f2 = -1.0;
    signalLowpassParam.a2 = -1.0;

    /* Call the Butterworth routine and check its success */
    if ( XLALButterworthREAL4TimeSeries(
                signal, &signalLowpassParam ) != XLAL_SUCCESS)
    {
        if (signal) LALFree( signal );
        XLAL_ERROR(func, XLAL_EFUNC);
    }

    /* That's it - we are done and can now free the signal */
    if (signal) LALFree( signal );

    return XLAL_SUCCESS;

}

