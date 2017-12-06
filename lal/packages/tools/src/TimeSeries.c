/*
 *
 * Copyright (C) 2007  Kipp Cannon
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


#include <complex.h>
#include <math.h>
#include <string.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/XLALError.h>


REAL4TimeSeries *XLALConvertREAL8TimeSeriesToREAL4(
	const REAL8TimeSeries *series
)
{
	REAL4TimeSeries *new = XLALCreateREAL4TimeSeries(series->name, &series->epoch, series->f0, series->deltaT, &series->sampleUnits, series->data->length);
	unsigned i;

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	for(i = 0; i < new->data->length; i++)
		new->data->data[i] = series->data->data[i];

	return new;
}


REAL8TimeSeries *XLALConvertREAL4TimeSeriesToREAL8(
	const REAL4TimeSeries *series
)
{
	REAL8TimeSeries *new = XLALCreateREAL8TimeSeries(series->name, &series->epoch, series->f0, series->deltaT, &series->sampleUnits, series->data->length);
	unsigned i;

	if(!new)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	for(i = 0; i < new->data->length; i++)
		new->data->data[i] = series->data->data[i];

	return new;
}


#define DATATYPE REAL4
#include "TimeSeries_source.c"
#undef DATATYPE

#define DATATYPE REAL8
#include "TimeSeries_source.c"
#undef DATATYPE

#define DATATYPE COMPLEX8
#include "TimeSeries_source.c"
#undef DATATYPE

#define DATATYPE COMPLEX16
#include "TimeSeries_source.c"
#undef DATATYPE

#define DATATYPE INT2
#include "TimeSeries_source.c"
#undef DATATYPE

#define DATATYPE UINT2
#include "TimeSeries_source.c"
#undef DATATYPE

#define DATATYPE INT4
#include "TimeSeries_source.c"
#undef DATATYPE

#define DATATYPE UINT4
#include "TimeSeries_source.c"
#undef DATATYPE

#define DATATYPE INT8
#include "TimeSeries_source.c"
#undef DATATYPE

#define DATATYPE UINT8
#include "TimeSeries_source.c"
#undef DATATYPE
