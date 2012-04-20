#include <math.h>
#include <lal/LALMoment.h>

#define TYPECODE D
#define TYPE REAL8
#include "LALMoment_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "LALMoment_source.c"
#undef TYPECODE
#undef TYPE
