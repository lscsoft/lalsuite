/**

\file
\ingroup AVFactories_h

\brief Create/destroy \<datatype\>Array objects.

\heading{Description}

The \c CreateArray family of functions create a
\<datatype\>Array of the appropriate dimensions.

The \c DestroyArray family of functions return the storage allocated by
the \c CreateArray functions to the system.

\heading{Algorithm}

\heading{Uses}
\code
LALMalloc()
LALFree()
\endcode

\heading{Notes}

*/

#include <string.h>
#include "LALStdlib.h"
#include "AVFactories.h"


NRCSID( ARRAYFACTORIESC, "$Id$" );


define(`TYPECODE',`Z')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`C')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`D')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`S')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I2')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I4')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`I8')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U2')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U4')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`U8')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')

define(`TYPECODE',`')
include(`CreateArray.m4')
include(`ResizeArray.m4')
include(`DestroyArray.m4')
