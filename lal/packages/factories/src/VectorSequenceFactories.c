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

#define TYPECODE Z
#define TYPE COMPLEX16
#include "CreateVectorSequence_source.c"
#include "DestroyVectorSequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE C
#define TYPE COMPLEX8
#include "CreateVectorSequence_source.c"
#include "DestroyVectorSequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE D
#define TYPE REAL8
#include "CreateVectorSequence_source.c"
#include "DestroyVectorSequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "CreateVectorSequence_source.c"
#include "DestroyVectorSequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I2
#define TYPE INT2
#include "CreateVectorSequence_source.c"
#include "DestroyVectorSequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I4
#define TYPE INT4
#include "CreateVectorSequence_source.c"
#include "DestroyVectorSequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I8
#define TYPE INT8
#include "CreateVectorSequence_source.c"
#include "DestroyVectorSequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U2
#define TYPE UINT2
#include "CreateVectorSequence_source.c"
#include "DestroyVectorSequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U4
#define TYPE UINT4
#include "CreateVectorSequence_source.c"
#include "DestroyVectorSequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U8
#define TYPE UINT8
#include "CreateVectorSequence_source.c"
#include "DestroyVectorSequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPE REAL4
#include "CreateVectorSequence_source.c"
#include "DestroyVectorSequence_source.c"
#undef TYPE
