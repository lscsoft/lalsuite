/*
 * Author: Kipp Cannon <kipp@gravity.phys.uwm.edu>
 * Revision: $Id$
 */

#include <string.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALErrno.h>
#include <lal/LALError.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>

NRCSID(TIMESERIESC, "$Id$");

define(`DATATYPE',REAL4)
include(TimeSeriesC.m4)

define(`DATATYPE',REAL8)
include(TimeSeriesC.m4)

define(`DATATYPE',COMPLEX8)
include(TimeSeriesC.m4)

define(`DATATYPE',COMPLEX16)
include(TimeSeriesC.m4)

define(`DATATYPE',INT2)
include(TimeSeriesC.m4)

define(`DATATYPE',UINT2)
include(TimeSeriesC.m4)

define(`DATATYPE',INT4)
include(TimeSeriesC.m4)

define(`DATATYPE',UINT4)
include(TimeSeriesC.m4)

define(`DATATYPE',INT8)
include(TimeSeriesC.m4)

define(`DATATYPE',UINT8)
include(TimeSeriesC.m4)
