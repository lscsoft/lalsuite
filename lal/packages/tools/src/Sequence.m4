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

NRCSID(SEQUENCEC, "$Id$");

/*
 * Shift the bytes in the buffer buff, whose length is length, count bytes to
 * higher addresses.  If the magnitude of count is greater than or equal to
 * that of length, then nothing is done.
 */

static void memshift(void *buff, size_t length, int count)
{
	if(count >= 0) {
		if(length > (size_t) count)
			memmove((char *) buff + count, buff, length - count);
	} else {
		if(length > (size_t) -count)
			memmove(buff, (char *) buff - count, length + count);
	}
}


/*
 * Begin library functions...
 */

define(`DATATYPE',REAL4)
include(SequenceC.m4)

define(`DATATYPE',REAL8)
include(SequenceC.m4)

define(`DATATYPE',COMPLEX8)
include(SequenceC.m4)

define(`DATATYPE',COMPLEX16)
include(SequenceC.m4)

define(`DATATYPE',INT2)
include(SequenceC.m4)

define(`DATATYPE',UINT2)
include(SequenceC.m4)

define(`DATATYPE',INT4)
include(SequenceC.m4)

define(`DATATYPE',UINT4)
include(SequenceC.m4)

define(`DATATYPE',INT8)
include(SequenceC.m4)

define(`DATATYPE',UINT8)
include(SequenceC.m4)
