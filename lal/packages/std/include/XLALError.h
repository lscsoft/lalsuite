#ifndef XLALERROR_H
#define XLALERROR_H

#include <lal/LALAtomicDatatypes.h>

NRCSID( XLALERRORH, "$Id$" );

#ifdef __cplusplus
extern "C" {
#pragma }
#endif


/*
 *
 * Use these functions to print arbitrary messages as errors or warnings.
 *
 */

/* prints an error message if error printing is enabled by lalDebugLevel */
int XLALPrintError( const char *fmt, ... );

/* prints a warning message if warning printing is enabled by lalDebugLevel */
int XLALPrintWarning( const char *fmt, ... );


/* silence gcc warnings about certain (possibly) unused symbols */
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*
 *
 * The LAL specification requires that XLAL functions that return a
 * floating-point type (either REAL4 or REAL8) should return a particular
 * value to indicate an error.  These values are given by the macros
 * XLAL_REAL4_FAIL_NAN and XLAL_REAL8_FAIL_NAN (they are Not a Number
 * or NaN values).  To implement these we choose hexadecimal representations
 * and then provide static functions that return the equivalent REAL4 or
 * REAL8 values.  The macros then invoke these functions.  This is done
 * so that the compiler can easily inline the functions (or eliminate them
 * if they are not used).  Conversion from the hexadecimal representation
 * to the floating-point representation is done using a union.
 *
 * The LAL specification also requires that there be two macros,
 * XLAL_IS_REAL4_FAIL_NAN(val) and XLAL_IS_REAL8_FAIL_NAN(val) that will
 * test if val is one of these XLAL-specific fail NaNs.  Again these macros
 * invoke static functions that return the result of the comparison.  The
 * comparison itself is done with the hexadecimal representation.
 *
 */

/* hexadecimal representation of the REAL4 and REAL8 NaN failure bit pattern */
#define XLAL_REAL4_FAIL_NAN_INT 0x7fc001a1
#define XLAL_REAL8_FAIL_NAN_INT LAL_INT8_C(0x7ff80000000001a1)

/* the floating point values themselves are returned by static functions that
 * can be easily inlined by the compiler; similarly, the routines to test
 * if a value is the LAL failure NaN can also be inlined */

/* returns the value of the XLAL REAL4 failure NaN */
static REAL4 UNUSED XLALREAL4FailNaN( void )
{
  const union { INT4 i; REAL4 x; } val = { XLAL_REAL4_FAIL_NAN_INT } ;
  return val.x;
}

/* returns the value of the XLAL REAL8 failure NaN */
static REAL8 UNUSED XLALREAL8FailNaN( void )
{
  const union { INT8 i; REAL8 x; } val = { XLAL_REAL8_FAIL_NAN_INT } ;
  return val.x;
}

/* tests if a value is an XLAL REAL4 failure NaN */
static int UNUSED XLALIsREAL4FailNaN( REAL4 val )
{
  return ( *(INT4 *)(&val) == XLAL_REAL4_FAIL_NAN_INT );
}

/* tests if a value is an XLAL REAL8 failure NaN */
static int UNUSED XLALIsREAL8FailNaN( REAL8 val )
{
  return ( *(INT8 *)(&val) == XLAL_REAL8_FAIL_NAN_INT );
}
#undef UNUSED

/* here are the macro constants for the fail NaNs */
#define XLAL_REAL4_FAIL_NAN ( XLALREAL4FailNaN() )
#define XLAL_REAL8_FAIL_NAN ( XLALREAL8FailNaN() )

/* here are the macros to test for fail NaNs */
#define XLAL_IS_REAL4_FAIL_NAN(val) XLALIsREAL4FailNaN(val)
#define XLAL_IS_REAL8_FAIL_NAN(val) XLALIsREAL8FailNaN(val)



/*
 *
 * The LAL specification requires particular return code and error values.
 * These are implemented here as enumeration constants.
 *
 */

/* XLAL error numbers and return values */
enum {
  /* these are result codes, not error numbers */
  XLAL_SUCCESS =  0, /* Success */
  XLAL_FAILURE = -1, /* Failure */

  /* these are standard error numbers */
  XLAL_EIO     =  5,  /* I/O error */
  XLAL_ENOMEM  = 12,  /* Memory allocation error */
  XLAL_EFAULT  = 14,  /* Invalid pointer */
  XLAL_EINVAL  = 22,  /* Invalid argument */
  XLAL_EDOM    = 33,  /* Input domain error */
  XLAL_ERANGE  = 34,  /* Output range error */

  /* extended error numbers start at 128 ... should be beyond normal errnos */

  /* these are common errors for XLAL functions */
  XLAL_EFAILED = 128, /* Generic failure */
  XLAL_EBADLEN = 129, /* Inconsistent or invalid vector length */

  /* specific mathematical and numerical errors start at 256 */

  /* IEEE floating point errors */
  XLAL_EFPINVAL  = 256, /* Invalid floating point operation, eg sqrt(-1), 0/0 */
  XLAL_EFPDIV0   = 257, /* Division by zero floating point error */
  XLAL_EFPOVRFLW = 258, /* Floating point overflow error */
  XLAL_EFPUNDFLW = 259, /* Floating point underflow error */
  XLAL_EFPINEXCT = 260, /* Floating point inexact error */

  /* numerical algorithm errors */
  XLAL_EMAXITER  = 261, /* Exceeded maximum number of iterations */
  XLAL_EDIVERGE  = 262, /* Series is diverging */
  XLAL_ESING     = 263, /* Apparent singularity detected */
  XLAL_ETOL      = 264, /* Failed to reach specified tolerance */
  XLAL_ELOSS     = 265, /* Loss of accuracy */

  /* failure from within a function call: "or" error number with this */
  XLAL_EFUNC     = 1024 /* Internal function call failed */
};

/*
 *
 * These functions provide message associated with an error code and print
 * an error message associated with the error code.  The macro XLAL_PERROR
 * fills in the current file and line information and uses the current
 * value of xlalErrno as the error number.
 *
 */

/* returns the error message associated with an error number */
const char * XLALErrorString( int errnum );

/* prints an error message for a particular error code in a standard format */
void XLALPerror( const char *func, const char *file, int line, int errnum );

/* prints an error message for the current value of xlalErrno */
#define XLAL_PERROR( func ) XLALPerror( func, __FILE__, __LINE__, xlalErrno )


/*
 *
 * Here is the XLAL error handler type and the routines that set it.
 * Also provide is the default error handler.
 *
 */

/* the XLAL error handler type */
typedef void XLALErrorHandlerType( const char *func, const char *file, int line, int errnum ); 

/* the default XLAL error handler */
void XLALDefaultErrorHandler( const char *func, const char *file, int line, int errnum );

/* sets the error handler to a new handler and returns the old handler */
XLALErrorHandlerType * XLALSetErrorHandler( XLALErrorHandlerType *newHandler );

/* sets the error handler to the default handler and returns the old handler */
XLALErrorHandlerType * XLALSetDefaultErrorHandler( void );


/*
 *
 * Here are the routines that set or clear the XLAL error number.
 *
 */

/* sets the XLAL error number to errnum */
void XLALSetErrno( int errnum );

/* clears the XLAL error number */
void XLALClearErrno( void );

/*
 *
 * The LAL specifiation requires that the XLAL error number be a modifiable
 * lvalue.  Similarly, the function pointer to the XLAL error handler is
 * a modifiable lvalue.  These are implemented as macros that dereference
 * pointers to the current value (in the current thread).  The pointer is
 * returned by the functions XLALGetErrnoPtr and XLALGetErrorHandlerPtr.
 * Here these functions and macros are defined.
 *
 */

/* function to return pointer to the XLAL error number */
int * XLALGetErrnoPtr( void );

/* function to return pointer to the XLAL error handler function pointer */
XLALErrorHandlerType ** XLALGetErrorHandlerPtr( void );

/* these are the modifiable lvalues for xlalErrno and XLALErrorHandler */
#define xlalErrno ( * XLALGetErrnoPtr() )
#define XLALErrorHandler ( * XLALGetErrorHandlerPtr() )


/*
 *
 * Here are the routines and macros that are used to report errors when
 * an XLAL function fails.  They (i) set the XLAL error number and (ii)
 * invoke the XLAL error handler.  The macros also (iii) return the
 * appropriate failure codes.  The macros should be used to report all
 * failures.
 *
 */

/* routine to set the XLAL error number and invoke the XLAL error handler...
 * it is used by the error macros below */
void XLALError( const char *func, const char *file, int line, int errnum );

/* invoke the XLALError function and return with code val... is should not
 * really be used itself, but forms the basis for the macros below */
#define XLAL_ERROR_VAL( func, errnum, val ) \
        do { \
          XLALError( func, __FILE__, __LINE__, errnum ); \
          return val; \
        } while (0)

/* how to invoke a failure from a XLAL routine returning an integer */
#define XLAL_ERROR( func, errnum ) \
    XLAL_ERROR_VAL( func, errnum, XLAL_FAILURE )

/* how to invoke a failure from a XLAL routine returning a pointer */
#define XLAL_ERROR_NULL( func, errnum ) \
    XLAL_ERROR_VAL( func, errnum, NULL )

/* how to invoke a failure from a XLAL routine returning void */
#define XLAL_ERROR_VOID( func, errnum ) \
    XLAL_ERROR_VAL( func, errnum, /* void */ )

/* how to invoke a failure from a XLAL routine returning a REAL4 */
#define XLAL_ERROR_REAL4( func, errnum ) \
    XLAL_ERROR_VAL( func, errnum, XLAL_REAL4_FAIL_NAN )

/* how to invoke a failure from a XLAL routine returning a REAL8 */
#define XLAL_ERROR_REAL8( func, errnum ) \
    XLAL_ERROR_VAL( func, errnum, XLAL_REAL8_FAIL_NAN )

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* XLALERROR_H */
