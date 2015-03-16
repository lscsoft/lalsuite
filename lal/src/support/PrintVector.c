#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALDatatypes.h>
#include <lal/PrintVector.h>

#define TYPECODE Z
#define TYPE COMPLEX16
#define FMT "%i %g %g\n"
#define ARG creal(vector->data[i]),cimag(vector->data[i])
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG

#define TYPECODE C
#define TYPE COMPLEX8
#define FMT "%i %g %g\n"
#define ARG crealf(vector->data[i]),cimagf(vector->data[i])
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG

#define TYPECODE D
#define TYPE REAL8
#define FMT "%i %g\n"
#define ARG vector->data[i]
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG

#define TYPECODE S
#define TYPE REAL4
#define FMT "%i %g\n"
#define ARG vector->data[i]
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG

#define TYPECODE I2
#define TYPE INT2
#define FMT "%i %i\n"
#define ARG vector->data[i]
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG

#define TYPECODE I4
#define TYPE INT4
#define FMT "%i %i\n"
#define ARG vector->data[i]
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG

/* Note that LALI8PrintVector does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */
#define TYPECODE I8
#define TYPE INT8
#define FMT "%i %0.0f\n"
#define ARG (REAL8)vector->data[i]
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG

#define TYPECODE U2
#define TYPE UINT2
#define FMT "%i %i\n"
#define ARG vector->data[i]
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG

#define TYPECODE U4
#define TYPE UINT4
#define FMT "%i %i\n"
#define ARG vector->data[i]
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG

/* Note that LALU8PrintVector does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */
#define TYPECODE U8
#define TYPE UINT8
#define FMT "%i %0.0f\n"
#define ARG (REAL8)vector->data[i]
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG

#define TYPECODE CHAR
#define TYPE CHAR
#define FMT "%i %c\n"
#define ARG vector->data[i]
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG

#define TYPECODE
#define TYPE REAL4
#define FMT "%i %f\n"
#define ARG vector->data[i]
#include "PrintVector_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG
