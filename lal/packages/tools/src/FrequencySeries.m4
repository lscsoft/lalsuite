/*
 * Author: Kipp Cannon <kipp@gravity.phys.uwm.edu>
 * Revision: $Id$
 */

#include <string.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALErrno.h>
#include <lal/LALError.h>
#include <lal/Sequence.h>
#include <lal/Units.h>
#include <lal/XLALError.h>

NRCSID(FREQUENCYSERIESC, "$Id$");

define(`DATATYPE',REAL4)
include(FrequencySeriesC.m4)

define(`DATATYPE',REAL8)
include(FrequencySeriesC.m4)

define(`DATATYPE',COMPLEX8)
include(FrequencySeriesC.m4)

define(`DATATYPE',COMPLEX16)
include(FrequencySeriesC.m4)

define(`DATATYPE',INT2)
include(FrequencySeriesC.m4)

define(`DATATYPE',UINT2)
include(FrequencySeriesC.m4)

define(`DATATYPE',INT4)
include(FrequencySeriesC.m4)

define(`DATATYPE',UINT4)
include(FrequencySeriesC.m4)

define(`DATATYPE',INT8)
include(FrequencySeriesC.m4)

define(`DATATYPE',UINT8)
include(FrequencySeriesC.m4)
