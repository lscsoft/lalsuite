/* <lalVerbatim file="SequenceCV">
Author: Kipp Cannon <kipp@gravity.phys.uwm.edu>
Revision: $Id$
</lalVerbatim>
 */

#include <string.h>
#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStatusMacros.h>
#include <lal/LALErrno.h>
#include <lal/LALError.h>
#include <lal/Sequence.h>

NRCSID(SEQUENCEC, "$Id$");

define(DATATYPE,COMPLEX8)
include(SequenceC.m4)

define(DATATYPE,COMPLEX16)
include(SequenceC.m4)

define(DATATYPE,REAL4)
include(SequenceC.m4)

define(DATATYPE,REAL8)
include(SequenceC.m4)

define(DATATYPE,INT2)
include(SequenceC.m4)

define(DATATYPE,INT4)
include(SequenceC.m4)

define(DATATYPE,INT8)
include(SequenceC.m4)

define(DATATYPE,UINT2)
include(SequenceC.m4)

define(DATATYPE,UINT4)
include(SequenceC.m4)

define(DATATYPE,UINT8)
include(SequenceC.m4)
