/*
 * $Id$
 *
 * Copyright (C) 2007  Kipp Cannon and Flanagan, E
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#ifndef _EPSEARCH_H
#define _EPSEARCH_H


#include <lal/LALDatatypes.h>
#include <lal/LALRCSID.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/Random.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Window.h>


#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif


NRCSID(EPSEARCHH, "$Id$");


/*
 * liblal.so can't resolve symbols from liblalsupport.so, so to call
 * diagnostics dump functions from lal, they have to be passed in as
 * pointers.  The prototypes here must match those in LIGOLwXMLArray.h or
 * memory problems will occur.  Since this is only used for diagnostics,
 * not production runs, there isn't any need to do this more robustly.
 * Note, also that the LIGOLwXMLStream structure is defined in a header
 * from liblalsupport, so here it has to be refered to as a void *.
 */


struct XLALEPSearchDiagnostics {
	void *LIGOLwXMLStream;
	int (*XLALWriteLIGOLwXMLArrayREAL4FrequencySeries)(void *, const char *, const REAL4FrequencySeries *);
	int (*XLALWriteLIGOLwXMLArrayREAL4TimeSeries)(void *, const char *, const REAL4TimeSeries *);
	int (*XLALWriteLIGOLwXMLArrayCOMPLEX8FrequencySeries)(void *, const char *, const COMPLEX8FrequencySeries *);
};


typedef struct tagEPSearchParams {
	struct XLALEPSearchDiagnostics *diagnostics;
	REAL4Window *window;
	REAL8 confidence_threshold;
	AvgSpecMethod method;
	/* time-frequency plane parameters */
	int useOverWhitening;	
	REAL8 flow;
	REAL8 bandwidth;
	/* t.f. plane tiling parameters */
	REAL8 fractional_stride;
	REAL8 maxTileBandwidth;
	REAL8 maxTileDuration;
} EPSearchParams;


SnglBurstTable *XLALEPSearch(
	const REAL4TimeSeries  *tseries,
	EPSearchParams   *params
);


int XLALEPConditionData(
	REAL4TimeSeries  *series,
	REAL8             flow,
	REAL8             resampledeltaT,
	INT4              corruption
);


COMPLEX8FrequencySeries *XLALWindowedREAL4ForwardFFT(
	const REAL4TimeSeries *tseries,
	const REAL4Window *window,
	const REAL4FFTPlan *plan
);


COMPLEX16FrequencySeries *XLALWindowedREAL8ForwardFFT(
	const REAL8TimeSeries *tseries,
	const REAL8Window *window,
	const REAL8FFTPlan *plan
);


REAL4Sequence *XLALREAL4AddWhiteNoise(
	REAL4Sequence *sequence,
	REAL4 rms,
	RandomParams *params
);


COMPLEX8Sequence *XLALCOMPLEX8AddWhiteNoise(
	COMPLEX8Sequence *sequence,
	REAL8 rms,
	RandomParams *params
);


#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif
