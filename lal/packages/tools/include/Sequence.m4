/*
 * Author: Kipp Cannon <kipp@gravity.phys.uwm.edu>
 * Revision: $Id$
 */

#ifndef _SEQUENCE_H
#define _SEQUENCE_H

#include <stddef.h>
#include <lal/LALDatatypes.h>

define(`DATATYPE',COMPLEX8)
include(SequenceH.m4)

define(`DATATYPE',COMPLEX16)
include(SequenceH.m4)

define(`DATATYPE',REAL4)
include(SequenceH.m4)

define(`DATATYPE',REAL8)
include(SequenceH.m4)

define(`DATATYPE',INT2)
include(SequenceH.m4)

define(`DATATYPE',INT4)
include(SequenceH.m4)

define(`DATATYPE',INT8)
include(SequenceH.m4)

define(`DATATYPE',UINT2)
include(SequenceH.m4)

define(`DATATYPE',UINT4)
include(SequenceH.m4)

define(`DATATYPE',UINT8)
include(SequenceH.m4)

#endif  /* _SEQUENCE_H */
