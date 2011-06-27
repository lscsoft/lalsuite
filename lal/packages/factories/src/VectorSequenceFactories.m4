/**
\file
\ingroup SeqFactories_h
\brief Create/destroy \<datatype\>VectorSequence objects.

\heading{Description}

The \c CreateVectorSequence family of functions create a
\<datatype\>VectorSequence of the appropriate dimensions.

The \c DestroyVectorSequence family of functions return the storage
allocated by the \c CreateVectorSequence functions to the system.

\heading{Algorithm}

\heading{Uses}
\code
LALMalloc()
LALFree()
\endcode

\heading{Notes}

*/

#include "LALStdlib.h"
#include "SeqFactories.h"


NRCSID( VECTORSEQUENCEFACTORIESC, "$Id$" );


define(`TYPECODE',`Z')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`C')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`D')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`S')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`I2')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`I4')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`I8')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`U2')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`U4')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`U8')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`CHAR')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')

define(`TYPECODE',`')
include(`CreateVectorSequence.m4')
include(`DestroyVectorSequence.m4')
