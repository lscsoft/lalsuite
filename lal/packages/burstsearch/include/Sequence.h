/* <lalVerbatim file="SequenceHV">
Author: Kipp Cannon <kipp@gravity.phys.uwm.edu>
$Id$
</lalVerbatim>
 */

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <lal/LALDatatypes.h>

void XLALDestroyCOMPLEX8Sequence(
	COMPLEX8Sequence *sequence
);

void LALDestroyCOMPLEX8Sequence(
	LALStatus *status,
	COMPLEX8Sequence *sequence
);

void XLALDestroyREAL4Sequence(
	REAL4Sequence *sequence
);

void LALDestroyREAL4Sequence(
	LALStatus *status,
	REAL4Sequence *sequence
);

COMPLEX8Sequence *XLALCreateCOMPLEX8Sequence(
	size_t length
);

void LALCreateCOMPLEX8Sequence(
	LALStatus *status,
	COMPLEX8Sequence **output,
	size_t length
);

REAL4Sequence *XLALCreateREAL4Sequence(
	size_t length
);

void LALCreateREAL4Sequence(
	LALStatus *status,
	REAL4Sequence **output,
	size_t length
);

REAL4Sequence *XLALCutREAL4Sequence(
	REAL4Sequence *sequence,
	size_t first,
	size_t length
);

void LALCutREAL4Sequence(
	LALStatus *status,
	REAL4Sequence **output,
	REAL4Sequence *input,
	size_t first,
	size_t length
);

REAL4Sequence *XLALShrinkREAL4Sequence(
	REAL4Sequence *sequence,
	size_t first,
	size_t length
);

void LALShrinkREAL4Sequence(
	LALStatus *status,
	REAL4Sequence **sequence,
	size_t first,
	size_t length
);

#endif  /* _SEQUENCE_H */
