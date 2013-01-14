#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <string.h>
#include <lal/LALDatatypes.h>
#include <lal/ReadFTSeries.h>

void LALCHARCreateVector( LALStatus *, CHARVector **, UINT4 );
void LALCHARDestroyVector( LALStatus *, CHARVector ** );
void LALDCreateVector( LALStatus *, REAL8Vector **, UINT4 );
void LALDDestroyVector( LALStatus *, REAL8Vector ** );
void LALParseUnitString ( LALStatus *status,
	       	          LALUnit *output,
		          const CHARVector *input );
static const LALUnit lalDimensionlessUnit 
	= {  0, { 0, 0, 0, 0, 0, 0, 0}, { 0, 0, 0, 0, 0, 0, 0} };

/* Change the first instance of the target to '\0'; returns 0 on success,
   1 on failure */
static INT2 changeCharToNull (CHAR *string, CHAR target, CHAR *offEndOfString)
{
  CHAR *charPtr;

  for ( charPtr=string; charPtr<offEndOfString; ++charPtr )
  {
    if ( *charPtr == target )
    {
      *charPtr = '\0';
      return 0;
    }
  }
  return 1;
}

#define TYPECODE Z
#define TYPE COMPLEX16
#define BASETYPE REAL8
#define FMT "%lf\t%lf\t%lf\n"
#define ARG &(data.array[0]),&(data.array[1])
#define NARGS 2
#include "ReadFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG
#undef NARGS

#define TYPECODE C
#define TYPE COMPLEX8
#define BASETYPE REAL4
#define FMT "%lf\t%f\t%f\n"
#define ARG &(data.array[0]),&(data.array[1])
#define NARGS 2
#include "ReadFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG
#undef NARGS

#define TYPECODE D
#define TYPE REAL8
#define FMT "%lf\t%lf\n"
#define ARG data.array
#define NARGS 1
#include "ReadFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG
#undef NARGS

#define TYPECODE S
#define TYPE REAL4
#define FMT "%lf\t%f\n"
#define ARG data.array
#define NARGS 1
#include "ReadFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG
#undef NARGS

#define TYPECODE
#define TYPE REAL4
#define FMT "%lf\t%f\n"
#define ARG data.array
#define NARGS 1
#include "ReadFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG
#undef NARGS
