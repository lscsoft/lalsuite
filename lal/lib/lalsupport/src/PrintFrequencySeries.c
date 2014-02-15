#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/PrintFTSeries.h>

#define TYPECODE Z
#define TYPE COMPLEX16
#define FMT "%e\t%le\t%le\n"
#define HEADER "# Freq (Hz)\tRe(Value)\tIm(Value)\n"
#define ARG creal(*data),cimag(*data)
#include "PrintFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE C
#define TYPE COMPLEX8
#define FMT "%e\t%e\t%e\n"
#define HEADER "# Freq (Hz)\tRe(Value)\tIm(Value)\n"
#define ARG crealf(*data),cimagf(*data)
#include "PrintFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE D
#define TYPE REAL8
#define FMT "%e\t%le\n"
#define HEADER "# Freq (Hz)\tValue\n"
#define ARG *data
#include "PrintFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE S
#define TYPE REAL4
#define FMT "%e\t%e\n"
#define HEADER "# Freq (Hz)\tValue\n"
#define ARG *data
#include "PrintFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE I2
#define TYPE INT2
#define FMT "%g\t%i\n"
#define HEADER "# Freq (Hz)\tValue\n"
#define ARG *data
#include "PrintFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE I4
#define TYPE INT4
#define FMT "%g\t%i\n"
#define HEADER "# Freq (Hz)\tValue\n"
#define ARG *data
#include "PrintFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

/* Note that LALI8PrintFrequencySeries does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */
#define TYPECODE I8
#define TYPE INT8
#define FMT "%g\t%0.0f\n"
#define HEADER "# Freq (Hz)\tValue\n"
#define ARG (REAL8)*data
#include "PrintFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE U2
#define TYPE UINT2
#define FMT "%g\t%i\n"
#define HEADER "# Freq (Hz)\tValue\n"
#define ARG *data
#include "PrintFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE U4
#define TYPE UINT4
#define FMT "%g\t%i\n"
#define HEADER "# Freq (Hz)\tValue\n"
#define ARG *data
#include "PrintFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

/* Note that LALU8PrintFrequencySeries does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */
#define TYPECODE U8
#define TYPE UINT8
#define FMT "%g\t%0.0f\n"
#define HEADER "# Freq (Hz)\tValue\n"
#define ARG (REAL8)*data
#include "PrintFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE
#define TYPE REAL4
#define FMT "%g\t%f\n"
#define HEADER "# Freq (Hz)\tValue\n"
#define ARG *data
#include "PrintFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG
