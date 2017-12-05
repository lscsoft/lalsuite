#include <complex.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <string.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/ReadFTSeries.h>

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
#include "ReadTimeSeries_source.c"
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
#include "ReadTimeSeries_source.c"
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
#include "ReadTimeSeries_source.c"
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
#include "ReadTimeSeries_source.c"
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
#include "ReadTimeSeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG
#undef NARGS
