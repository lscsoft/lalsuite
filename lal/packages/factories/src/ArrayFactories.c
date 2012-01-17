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

#define TYPECODE Z
#define TYPE COMPLEX16
#include "CreateArray_source.c"
#include "ResizeArray_source.c"
#include "DestroyArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE C
#define TYPE COMPLEX8
#include "CreateArray_source.c"
#include "ResizeArray_source.c"
#include "DestroyArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE D
#define TYPE REAL8
#include "CreateArray_source.c"
#include "ResizeArray_source.c"
#include "DestroyArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "CreateArray_source.c"
#include "ResizeArray_source.c"
#include "DestroyArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I2
#define TYPE INT2
#include "CreateArray_source.c"
#include "ResizeArray_source.c"
#include "DestroyArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I4
#define TYPE INT4
#include "CreateArray_source.c"
#include "ResizeArray_source.c"
#include "DestroyArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I8
#define TYPE INT8
#include "CreateArray_source.c"
#include "ResizeArray_source.c"
#include "DestroyArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U2
#define TYPE UINT2
#include "CreateArray_source.c"
#include "ResizeArray_source.c"
#include "DestroyArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U4
#define TYPE UINT4
#include "CreateArray_source.c"
#include "ResizeArray_source.c"
#include "DestroyArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U8
#define TYPE UINT8
#include "CreateArray_source.c"
#include "ResizeArray_source.c"
#include "DestroyArray_source.c"
#undef TYPECODE
#undef TYPE

#define TYPE REAL4
#include "CreateArray_source.c"
#include "ResizeArray_source.c"
#include "DestroyArray_source.c"
#undef TYPE
