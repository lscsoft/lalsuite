/*
 * Author: Kipp Cannon <kipp@gravity.phys.uwm.edu>
 * Revision: $Id$
 */

#ifndef _TIMESERIES_H
#define _TIMESERIES_H

#include <stddef.h>
#include <lal/LALDatatypes.h>

define(`DATATYPE',COMPLEX8)
include(TimeSeriesH.m4)

define(`DATATYPE',COMPLEX16)
include(TimeSeriesH.m4)

define(`DATATYPE',REAL4)
include(TimeSeriesH.m4)

define(`DATATYPE',REAL8)
include(TimeSeriesH.m4)

define(`DATATYPE',INT2)
include(TimeSeriesH.m4)

define(`DATATYPE',INT4)
include(TimeSeriesH.m4)

define(`DATATYPE',INT8)
include(TimeSeriesH.m4)

define(`DATATYPE',UINT2)
include(TimeSeriesH.m4)

define(`DATATYPE',UINT4)
include(TimeSeriesH.m4)

define(`DATATYPE',UINT8)
include(TimeSeriesH.m4)

#endif  /* _TIMESERIES_H */
