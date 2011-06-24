/*-----------------------------------------------------------------------

File Name: VectorFactories.c

-------------------------------------------------------------------------*/

/**
\file
\ingroup AVFactories_h

\brief Create/destroy \<datatype\>%Vector objects.

\heading{Description}

The \c CreateVector family of functions create a
\<datatype\>%Vector of the appropriate dimensions.

The \c ResizeVector family of functions changes the amount of
storage allocated by the \c CreateVector functions.

The \c DestroyVector family of functions return the storage allocated by
the \c CreateVector functions to the system.

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


NRCSID( VECTORFACTORIESC, "$Id$" );


define(`TYPECODE',`Z')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`C')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`D')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`S')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I2')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I4')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`I8')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U2')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U4')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`U8')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`CHAR')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')

define(`TYPECODE',`')
include(`CreateVector.m4')
include(`ResizeVector.m4')
include(`DestroyVector.m4')
