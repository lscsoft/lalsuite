/**
\author Whelan, J. T.
\file

\heading{Module \ref PrintFrequencySeries.c}
\latexonly\label{ss_PrintFrequencySeries_c}\endlatexonly

Print a \f$\langle\mbox{datatype}\rangle\f$FrequencySeries object into a
file.  For use in non-production and test code only.

\heading{Prototypes}














\heading{Description}

Each member of this family of functions prints the elements of
\f$\langle\mbox{datatype}\rangle\f$\c FrequencySeries into a file.
Note: the file name is specified using a character string.  This
function is for debugging use only: its arguments do not conform to
LAL standards so it should not be used in any real analysis codes.

\heading{Algorithm}

\heading{Uses}

\code
LALFopen()
LALFclose()
LALCHARCreateVector()
LALCHARDestroyVector()
LALUnitAsString()
\endcode

\heading{Notes}

This function's arguments do not conform to the LAL spec.  For this
reason it should only be used for debugging purposes in test
functions, not in any production code.

Additionally, since printf cannot handle INT8 as integers, the
functions <tt>LALI8PrintFrequencySeries()</tt> and
<tt>LALU8PrintFrequencySeries()</tt> use a typecast to REAL8 and are
thus only valid for numbers between around \f$-10^{15}\f$ and \f$10^{15}\f$.

The first four lines of the file are a header containing:
<ol>
<li> the name of the series</li>
<li> heterodyning information, if any</li>
<li> the starting epoch, relative to the GPS reference epoch (1980 January 6)</li>
<li> the units expressed in terms of the basic SI units</li>
<li> column labels</li>
</ol>
after which come the data, one per line.

The output format is two or three tab-separated columns: the first
column is the frequency in hertz corresponding to the row in question;
for real and integer frequency series, the second column is the value
of the series; for complex frequency series, the second column is the
real part and the third the imaginary part of the value.

Note that the frequency given is the physical frequency.  In the case
of a heterodyned frequency series, this is the heterodyning frequency
plus the frequency offset.  A frequency series of length \f$[N]\f$ is
assumed to be packed so that the 0th element corresponds to zero
frequency offset, elements 1 through \f$[N/2]\f$ to positive frequency
offsets (in ascending order), and elements \f$N-[N/2]\f$ to \f$N-1\f$ to
negative frequency offsets (also in ascending order, so that the
frequency corresponding to the \f$N-1\f$st element is just below that of
the 0th element).  If \f$N\f$ is even, the element in position \f$N/2\f$ is
assumed to correspond both the maximum poitive and negative frequency
offset.



*/


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


NRCSID( PRINTFREQUENCYSERIESC, "$Id$" );


#define TYPECODE Z
#define TYPE COMPLEX16
#define FMT "%e\t%le\t%le\n"
#define HEADER "# Freq (Hz)\tRe(Value)\tIm(Value)\n"
#define ARG data->re,data->im
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
#define ARG data->re,data->im
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
