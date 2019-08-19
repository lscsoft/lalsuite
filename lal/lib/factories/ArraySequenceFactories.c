#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>

#define TYPECODE Z
#define TYPE COMPLEX16
#include "CreateArraySequence_source.c"
#include "DestroyArraySequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE C
#define TYPE COMPLEX8
#include "CreateArraySequence_source.c"
#include "DestroyArraySequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE D
#define TYPE REAL8
#include "CreateArraySequence_source.c"
#include "DestroyArraySequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "CreateArraySequence_source.c"
#include "DestroyArraySequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I2
#define TYPE INT2
#include "CreateArraySequence_source.c"
#include "DestroyArraySequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I4
#define TYPE INT4
#include "CreateArraySequence_source.c"
#include "DestroyArraySequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I8
#define TYPE INT8
#include "CreateArraySequence_source.c"
#include "DestroyArraySequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U2
#define TYPE UINT2
#include "CreateArraySequence_source.c"
#include "DestroyArraySequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U4
#define TYPE UINT4
#include "CreateArraySequence_source.c"
#include "DestroyArraySequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U8
#define TYPE UINT8
#include "CreateArraySequence_source.c"
#include "DestroyArraySequence_source.c"
#undef TYPECODE
#undef TYPE

#define TYPE REAL4
#include "CreateArraySequence_source.c"
#include "DestroyArraySequence_source.c"
#undef TYPE
