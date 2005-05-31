/*
 * Author: Kipp Cannon <kipp@gravity.phys.uwm.edu>
 * Revision: $Id$
 */

#ifndef _FREQUENCYSERIES_H
#define _FREQUENCYSERIES_H

#include <stddef.h>
#include <lal/LALDatatypes.h>

define(`DATATYPE',COMPLEX8)
include(FrequencySeriesH.m4)

define(`DATATYPE',COMPLEX16)
include(FrequencySeriesH.m4)

define(`DATATYPE',REAL4)
include(FrequencySeriesH.m4)

define(`DATATYPE',REAL8)
include(FrequencySeriesH.m4)

define(`DATATYPE',INT2)
include(FrequencySeriesH.m4)

define(`DATATYPE',INT4)
include(FrequencySeriesH.m4)

define(`DATATYPE',INT8)
include(FrequencySeriesH.m4)

define(`DATATYPE',UINT2)
include(FrequencySeriesH.m4)

define(`DATATYPE',UINT4)
include(FrequencySeriesH.m4)

define(`DATATYPE',UINT8)
include(FrequencySeriesH.m4)

#endif  /* _FREQUENCYSERIES_H */
