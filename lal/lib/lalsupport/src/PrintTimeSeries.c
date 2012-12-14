#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALDatatypes.h>
#include <lal/PrintFTSeries.h>

void LALCHARCreateVector( LALStatus *, CHARVector **, UINT4 );
void LALCHARDestroyVector( LALStatus *, CHARVector ** );
void LALUnitAsString( LALStatus *status, CHARVector *output,
                      const LALUnit *input );

enum { LALUnitTextSize = sizeof("10^-32768 m^-32768/32767 kg^-32768/32767 "
				"s^-32768/32767 A^-32768/32767 " 
				"K^-32768/32767 strain^-32768/32767 "
				"count^-32768/32767") };

#define TYPECODE Z
#define TYPE COMPLEX16
#define FMT "%e\t%le\t%le\n"
#define HEADER "# Seconds since epoch\tRe(Value)\tIm(Value)\n"
#define ARG creal(*data),cimag(*data)
#include "PrintTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE C
#define TYPE COMPLEX8
#define FMT "%e\t%e\t%e\n"
#define HEADER "# Seconds since epoch\tRe(Value)\tIm(Value)\n"
#define ARG crealf(*data),cimagf(*data)
#include "PrintTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE D
#define TYPE REAL8
#define FMT "%e\t%le\n"
#define HEADER "# Seconds since epoch\tValue\n"
#define ARG *data
#include "PrintTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE S
#define TYPE REAL4
#define FMT "%e\t%e\n"
#define HEADER "# Seconds since epoch\tValue\n"
#define ARG *data
#include "PrintTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE I2
#define TYPE INT2
#define FMT "%g\t%i\n"
#define HEADER "# Seconds since epoch\tValue\n"
#define ARG *data
#include "PrintTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE I4
#define TYPE INT4
#define FMT "%g\t%i\n"
#define HEADER "# Seconds since epoch\tValue\n"
#define ARG *data
#include "PrintTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

/* Note that LALI8PrintTimeSeries does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */
#define TYPECODE I8
#define TYPE INT8
#define FMT "%g\t%0.0f\n"
#define HEADER "# Seconds since epoch\tValue\n"
#define ARG (REAL8)*data
#include "PrintTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE U2
#define TYPE UINT2
#define FMT "%g\t%i\n"
#define HEADER "# Seconds since epoch\tValue\n"
#define ARG *data
#include "PrintTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE U4
#define TYPE UINT4
#define FMT "%g\t%i\n"
#define HEADER "# Seconds since epoch\tValue\n"
#define ARG *data
#include "PrintTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

/* Note that LALU8PrintTimeSeries does a typecast to REAL8 and is thus
 * inaccurate for numbers >~ 1e15 
 */
#define TYPECODE U8
#define TYPE UINT8
#define FMT "%g\t%0.0f\n"
#define HEADER "# Seconds since epoch\tValue\n"
#define ARG (REAL8)*data
#include "PrintTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG

#define TYPECODE
#define TYPE REAL4
#define FMT "%g\t%f\n"
#define HEADER "# Seconds since epoch\tValue\n"
#define ARG *data
#include "PrintTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef HEADER
#undef ARG
