/**

\file
\ingroup SeqFactories_h

\brief Create/destroy \<datatype\>ArraySequence objects.

\heading{Description}

The \c CreateArraySequence family of functions create a
\<datatype\>ArraySequence of the appropriate dimensions.

The \c DestroyArraySequence family of functions return the storage
allocated by the \c CreateArraySequence functions to the system.

\heading{Algorithm}

\heading{Uses}
\code
LALMalloc()
LALFree()
\endcode

\heading{Notes}

*/

#include "LALStdlib.h"
#include "AVFactories.h"
#include "SeqFactories.h"


NRCSID( ARRAYSEQUENCEFACTORIESC, "$Id$" );


define(`TYPECODE',`Z')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`C')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`D')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`S')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`I2')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`I4')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`I8')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`U2')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`U4')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`U8')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')

define(`TYPECODE',`')
include(`CreateArraySequence.m4')
include(`DestroyArraySequence.m4')
