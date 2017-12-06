/*
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
#include <lal/LIGOMetadataTables.h>
#include <lal/Random.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Window.h>


#ifdef  __cplusplus   /* C++ protection. */
extern "C" {
#endif

/**
 * \defgroup EPSearch_h Header EPSearch.h
 * \ingroup pkg_burstsearch
 *
 * \brief UNDOCUMENTED
 */
  /*@{*/


/**
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
	int (*XLALWriteLIGOLwXMLArrayREAL8FrequencySeries)(void *, const char *, const REAL8FrequencySeries *);
	int (*XLALWriteLIGOLwXMLArrayREAL8TimeSeries)(void *, const char *, const REAL8TimeSeries *);
	int (*XLALWriteLIGOLwXMLArrayCOMPLEX16FrequencySeries)(void *, const char *, const COMPLEX16FrequencySeries *);
};


SnglBurst *XLALEPSearch(
	struct XLALEPSearchDiagnostics *diagnostics,
	const REAL8TimeSeries  *tseries,
	REAL8Window *window,
	REAL8 flow,
	REAL8 bandwidth,
	REAL8 confidence_threshold,
	/* t.f. plane tiling parameters */
	REAL8 fractional_stride,
	REAL8 maxTileBandwidth,
	REAL8 maxTileDuration
);

  /*@}*/

#ifdef  __cplusplus
}
#endif  /* C++ protection. */

#endif
