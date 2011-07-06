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

#define TYPECODE Z
#define TYPE COMPLEX16
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE C
#define TYPE COMPLEX8
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE D
#define TYPE REAL8
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I2
#define TYPE INT2
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I4
#define TYPE INT4
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I8
#define TYPE INT8
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U2
#define TYPE UINT2
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U4
#define TYPE UINT4
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U8
#define TYPE UINT8
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE CHAR
#define TYPE CHAR
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPECODE
#undef TYPE

#define TYPE REAL4
#include "CreateVector_source.c"
#include "DestroyVector_source.c"
#include "ResizeVector_source.c"
#undef TYPE
