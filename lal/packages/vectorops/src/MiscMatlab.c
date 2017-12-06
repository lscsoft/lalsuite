#include <math.h>
#include <lal/Matrix.h>

#define TYPECODE D
#define TYPE REAL8
#include "MiscMatlab_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "MiscMatlab_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I2
#define TYPE INT2
#include "MiscMatlab_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I4
#define TYPE INT4
#include "MiscMatlab_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I8
#define TYPE INT8
#include "MiscMatlab_source.c"
#undef TYPECODE
#undef TYPE
