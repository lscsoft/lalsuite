/**
\author Torres, C. V.
\file

\heading{Module \ref ReadFrequencySeries.c}
\latexonly\label{ss_ReadFrequencySeries_c}\endlatexonly

\heading{Prototypes}








\heading{Description}

Each member of this family of functions reads from a file the output
of the corresponding \c PrintFrequencySeries routine.

\heading{Algorithm}

\heading{Uses}

\code
LALOpenDataFile()
LALParseUnitString()
LALCHARCreateVector()
LALCHARDestroyVector()
LALDCreateVector()
LALDDestroyVector()
\endcode

\heading{Notes}

These functions perform I/O operations, which are not a part of LAL
proper. They should only be used for debugging purposes in test
functions, not in any production code.



*/

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


NRCSID( READFREQUENCYSERIESC, "$Id$" );


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
#define FMT "%lf\t%lf\t%lf\n"
#define ARG &(data.re),&(data.im)
#define NARGS 2
#include "ReadFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG
#undef NARGS

#define TYPECODE C
#define TYPE COMPLEX8
#define FMT "%lf\t%f\t%f\n"
#define ARG &(data.re),&(data.im)
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
#define ARG &data
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
#define ARG &data
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
#define ARG &data
#define NARGS 1
#include "ReadFrequencySeries_source.c"
#undef TYPECODE
#undef TYPE
#undef FMT
#undef ARG
#undef NARGS
