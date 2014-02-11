/*
 * Copyright (C) 2014 Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */


#ifndef _TIMESERIESINTERP_H_
#define _TIMESERIESINTERP_H_


#include <lal/LALDatatypes.h>


#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/**
 * Opaque LALREAL8TimeSeriesInterp structure.
 */

typedef struct tagLALREAL8TimeSeriesInterp LALREAL8TimeSeriesInterp;


LALREAL8TimeSeriesInterp *XLALREAL8TimeSeriesInterpCreate(const REAL8TimeSeries *, int);
void XLALREAL8TimeSeriesInterpDestroy(LALREAL8TimeSeriesInterp *);
REAL8 XLALREAL8TimeSeriesInterpEval(LALREAL8TimeSeriesInterp *, const LIGOTimeGPS *);


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif	/* _TIMESERIESINTERP_H_ */
