/* -*- C -*- */

#include <math.h>
#include "VectorIndexRange.h"

NRCSID( VECTORINDEXRANGEC, "$Id$" );

define(`TYPECODE',`CHAR')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`I2')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`I4')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`I8')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`U2')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`U4')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`U8')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`S')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`D')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`C')
include(`VectorIndexRangeBaseC.m4')

define(`TYPECODE',`Z')
include(`VectorIndexRangeBaseC.m4')



