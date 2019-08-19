/*
 * Copyright (C) 2007 Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */


#ifndef _LIGOLWXMLARRAY_H
#define _LIGOLWXMLARRAY_H


#include <lal/LIGOLwXML.h>
#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


int XLALWriteLIGOLwXMLArrayREAL4TimeSeries(
	LIGOLwXMLStream *xml,
	const char *comment,
	const REAL4TimeSeries *series
);


int XLALWriteLIGOLwXMLArrayREAL8TimeSeries(
	LIGOLwXMLStream *xml,
	const char *comment,
	const REAL8TimeSeries *series
);


int XLALWriteLIGOLwXMLArrayREAL4FrequencySeries(
	LIGOLwXMLStream *xml,
	const char *comment,
	const REAL4FrequencySeries *series
);


int XLALWriteLIGOLwXMLArrayREAL8FrequencySeries(
	LIGOLwXMLStream *xml,
	const char *comment,
	const REAL8FrequencySeries *series
);


int XLALWriteLIGOLwXMLArrayCOMPLEX8FrequencySeries(
	LIGOLwXMLStream *xml,
	const char *comment,
	const COMPLEX8FrequencySeries *series
);


int XLALWriteLIGOLwXMLArrayCOMPLEX16FrequencySeries(
	LIGOLwXMLStream *xml,
	const char *comment,
	const COMPLEX16FrequencySeries *series
);


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LIGOLIGOLWXMLARRAY_H */
