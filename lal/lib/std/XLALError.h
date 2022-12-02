/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

#ifndef XLALERROR_H
#define XLALERROR_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALString.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif


/**
 * \defgroup XLALError_h Header XLALError.h
 * \ingroup lal_std
 * \author Creighton, J. D. E.
 * \date 2005
 * \brief This header covers routines to provide the XLAL interface error
 * handling.
 *
 * ### XLAL Errors ###
 *
 * When an XLAL routine fails, the routine should set the <tt>xlalErrno</tt> to
 * an appropriate error number and return with the appropriate error code.  The
 * return value depends on the return type of the XLAL function.  Furthermore,
 * the XLAL error handler should be invoked.
 *
 * Whenever possible (i.e., always), standard XLAL error macros should be used
 * when generating an error.  These macros (i) invoke the current error handler,
 * (ii) set the error code to the specified value, and (iii) return with the
 * correct return value.  In addition, these macros may take an optional
 * printf-like format string (along with additional parameters for this format
 * string) to provide additional information about the nature of the failure.
 * The error macros that should be used are:
 *
 * <tt> #XLAL_ERROR(errnum [, fmt [, ...]])</tt> for XLAL routines returning an
 * integer type.
 *
 * <tt> #XLAL_ERROR_VOID(errnum [, fmt [, ...]])</tt> for XLAL routines with no
 * return value.
 *
 * <tt> #XLAL_ERROR_NULL(errnum [, fmt [, ...]])</tt> for XLAL routines returning
 * a pointer.
 *
 * <tt> #XLAL_ERROR_REAL4(errnum [, fmt [, ...]])</tt> for XLAL routines
 * returning a <tt>REAL4</tt> floating-point value.
 *
 * <tt> #XLAL_ERROR_REAL8(errnum [, fmt [, ...]])</tt> for XLAL routines
 * returning a <tt>REAL8</tt> floating-point value.
 *
 * Assert-like error checking can be performed with <tt>#XLAL_CHECK</tt>-style
 * macros.  Unlike <tt>assert()</tt> statements, <tt>#XLAL_CHECK</tt> macros
 * do <i>not</i> get removed when the code is compiled with <tt>-DNDEBUG</tt>.
 *
 * Additional error, warning, and informational messages can be generated using
 * the routines <tt>XLALPrintError()</tt>, <tt>XLALPrintWarning()</tt> and
 * <tt>XLALPrintInfo()</tt>.  These routines (which work just like
 * <tt>printf()</tt>) print or suppress the message depending on the value of
 * <tt>lalDebugLevel</tt>.  To print error/warning/info messages with a
 * standard format, use the macros
 * <tt>#XLAL_PRINT_ERROR(fmt [, ...])</tt>
 * <tt>#XLAL_PRINT_WARNING(fmt [, ...])</tt>
 * <tt>#XLAL_PRINT_INFO(fmt [, ...])</tt>
 *
 * On rare occations, you may be prepared for an XLAL routine to fail, and may
 * want to handle the failure immediately.  In these circumstances, the XLAL
 * error handler needs to be disabled before the routine is called so that the
 * failure can be caught.  The <tt>#XLAL_TRY(statement,errnum)</tt> macro is
 * designed to be used in these situations.  Here is an example:
 * \code
 * REAL8 XLALLogFactorial(INT4 n)
 * {
 * 	REAL8 y;
 * 	int errnum;
 * 	XLAL_TRY(y = XLALGammaFunction(n + 1), errnum);
 * 	if (XLAL_IS_REAL8_FAIL_NAN(y))
 * 		switch (errnum) {
 * 			case XLAL_ERANGE:
 * 				y  = n * (log(n) - 1);
 * 				y += 0.5 * log(2.0 * LAL_PI * n);
 * 				return y;
 * 			default:
 * 				XLALSetErrno(errnum);
 * 				XLAL_ERROR_REAL8(XLAL_EFUNC);
 * 		}
 * 	return log(y);
 * }
 * \endcode
 *
 * ### XLAL Function Return Codes ###
 *
 * XLAL functions that return an integer-type will return <tt>#XLAL_FAILURE</tt>
 * on failure.  XLAL functions that return a pointer will return <tt>NULL</tt>
 * on failure.
 *
 * The LAL specification requires that XLAL functions that return a
 * floating-point type (either ::REAL4 or ::REAL8) should return
 * a particular value to indicate an error.  These values are given by the
 * macros <tt>#XLAL_REAL4_FAIL_NAN</tt> and <tt>#XLAL_REAL8_FAIL_NAN</tt> (they
 * are Not a Number or NaN values).  To implement these we choose hexadecimal
 * representations and then provide static functions that return the equivalent
 * ::REAL4 or ::REAL8 values.  The macros then invoke these
 * functions.  This is done so that the compiler can easily inline the
 * functions (or eliminate them if they are not used).  Conversion from the
 * hexadecimal representation to the floating-point representation is done
 * using a union.
 *
 * The LAL specification also requires that there be two macros,
 * <tt>#XLAL_IS_REAL4_FAIL_NAN(val)</tt> and
 * <tt>#XLAL_IS_REAL8_FAIL_NAN(val)</tt> that will test if val is one of these
 * XLAL-specific fail NaNs.  Again these macros invoke static functions that
 * return the result of the comparison.  The cmparison itself is done with the
 * hexadecimal representation.
 *
 * ### XLAL Error Codes ###
 *
 * The LAL specification requires particular return code and error values.
 * These are implemented here as enumeration constants in the
 * ::XLALErrorValue enumeration.
 */
/** @{ */


#ifndef SWIG    /* exclude from SWIG interface */


/*
 *
 * Use these functions to print arbitrary messages as errors or warnings.
 *
 */

/** Prints an error message if error printing is enabled by lalDebugLevel. */
int XLALPrintError(const char *fmt, ...) _LAL_GCC_PRINTF_FORMAT_(1,2);

/** Prints a warning message if warning printing is enabled by lalDebugLevel. */
int XLALPrintWarning(const char *fmt, ...) _LAL_GCC_PRINTF_FORMAT_(1,2);

/** Prints an info message if info printing is enabled by lalDebugLevel. */
int XLALPrintInfo(const char *fmt, ...) _LAL_GCC_PRINTF_FORMAT_(1,2);

/** Prints an error message if error printing is enabled by lalDebugLevel. */
int XLALVPrintError(const char *fmt, va_list ap) _LAL_GCC_VPRINTF_FORMAT_(1);

/** Prints a warning message if warning printing is enabled by lalDebugLevel. */
int XLALVPrintWarning(const char *fmt, va_list ap) _LAL_GCC_VPRINTF_FORMAT_(1);

/** Prints an info message if info printing is enabled by lalDebugLevel. */
int XLALVPrintInfo(const char *fmt, va_list ap) _LAL_GCC_VPRINTF_FORMAT_(1);


/*
 *
 * Miscelaneous routines to print information with standard formatting.
 *
 */

/**
 * Print an error message with standard XLAL formatting (if error messages
 * are enabled by lalDebugLevel).
 */
void XLALPrintErrorMessage(const char *func, const char *file, int line,
                           const char *fmt, ...) _LAL_GCC_PRINTF_FORMAT_(4,5);

/**
 * Print an warning message with standard XLAL formatting (if warning messages
 * are enabled by lalDebugLevel).
 */
void XLALPrintWarningMessage(const char *func, const char *file, int line,
                             const char *fmt, ...) _LAL_GCC_PRINTF_FORMAT_(4,5);

/**
 * Print an info message with standard XLAL formatting (if info messages
 * are enabled by lalDebugLevel).
 */
void XLALPrintInfoMessage(const char *func, const char *file, int line,
                          const char *fmt, ...) _LAL_GCC_PRINTF_FORMAT_(4,5);

/**
 * Print an error message with standard XLAL formatting (if error messages
 * are enabled by lalDebugLevel).
 */
void XLALVPrintErrorMessage(const char *func, const char *file, int line,
                            const char *fmt, va_list ap) _LAL_GCC_VPRINTF_FORMAT_(4);

/**
 * Print an warning message with standard XLAL formatting (if warning messages
 * are enabled by lalDebugLevel).
 */
void XLALVPrintWarningMessage(const char *func, const char *file, int line,
                              const char *fmt, va_list ap) _LAL_GCC_VPRINTF_FORMAT_(4);

/**
 * Print an error message with standard XLAL formatting (if error messages
 * are enabled by lalDebugLevel).
 */
void XLALVPrintInfoMessage(const char *func, const char *file, int line,
                           const char *fmt, va_list ap) _LAL_GCC_VPRINTF_FORMAT_(4);

/** Prints a progress bar at the "info" verbosity level. */
int XLALPrintProgressBar(double);

/** Prints a deprecation warning at the "warning" verbosity level. */
#define XLAL_PRINT_DEPRECATION_WARNING(replacement) \
  do { \
    static int _xlal_print_deprecation_warning_ = 1; \
    if (_xlal_print_deprecation_warning_) { \
      XLALPrintWarning( \
        "\nDEPRECATION WARNING: program has invoked obsolete function %s(). " \
        "Please see %s() for information about a replacement.\n", \
        __func__, replacement); \
      _xlal_print_deprecation_warning_ = 0; \
    } \
  } while(0)


/*
 *
 * Macros that will print error/warning/info messages with a standard format.
 *
 */

/**
 * \brief Macro that will print an error message with a standard format.
 *
 * Prototype: <b>XLAL_PRINT_ERROR(fmt [, ...])</b>
 *
 * \b Parameters:<ul>
 * <li> \b fmt A printf-like format string.
 * <li> \b ... (Optional) Arguments to the format string.
 * </ul>
 */
#define XLAL_PRINT_ERROR(...) \
	XLALPrintErrorMessage(__func__, __FILE__, __LINE__, __VA_ARGS__)

/**
 * \brief Macro that will print a warning message with a standard format.
 *
 * Prototype: <b>XLAL_PRINT_WARNING(fmt [, ...])</b>
 *
 * \b Parameters:<ul>
 * <li> \b fmt A printf-like format string.
 * <li> \b ... (Optional) Arguments to the format string.
 * </ul>
 */
#define XLAL_PRINT_WARNING(...) \
	XLALPrintWarningMessage(__func__, __FILE__, __LINE__, __VA_ARGS__)

/**
 * \brief Macro that will print an info message with a standard format.
 *
 * Prototype: <b>XLAL_PRINT_INFO(fmt [, ...])</b>
 *
 * \b Parameters:<ul>
 * <li> \b fmt A printf-like format string.
 * <li> \b ... (Optional) Arguments to the format string.
 * </ul>
 */
#define XLAL_PRINT_INFO(...) \
	XLALPrintInfoMessage(__func__, __FILE__, __LINE__, __VA_ARGS__)


/*
 *
 * The LAL specification requires that XLAL functions that return a
 * floating-point type (either <tt>REAL4</tt> or <tt>REAL8</tt>) should return
 * a particular value to indicate an error.  These values are given by the
 * macros <tt>XLAL_REAL4_FAIL_NAN</tt> and <tt>XLAL_REAL8_FAIL_NAN</tt> (they
 * are Not a Number or NaN values).  To implement these we choose hexadecimal
 * representations and then provide static functions that return the equivalent
 * <tt>REAL4</tt> or <tt>REAL8</tt> values.  The macros then invoke these
 * functions.  This is done so that the compiler can easily inline the
 * functions (or eliminate them if they are not used).  Conversion from the
 * hexadecimal representation to the floating-point representation is done
 * using a union.
 *
 * The LAL specification also requires that there be two macros,
 * <tt>XLAL_IS_REAL4_FAIL_NAN(val)</tt> and
 * <tt>XLAL_IS_REAL8_FAIL_NAN(val)</tt> that will
 * test if val is one of these XLAL-specific fail NaNs.  Again these macros
 * invoke static functions that return the result of the comparison.  The
 * comparison itself is done with the hexadecimal representation.
 *
 */

/* Hexadecimal representation of the <tt>REAL4</tt> and <tt>REAL8</tt> NaN
 * failure bit pattern. */
#define XLAL_REAL4_FAIL_NAN_INT 0x7fc001a1 /**< Hexadecimal representation of <tt>REAL4</tt> NaN failure bit pattern */
#define XLAL_REAL8_FAIL_NAN_INT LAL_INT8_C(0x7ff80000000001a1) /**< Hexadecimal representation of <tt>REAL8</tt> NaN failure bit pattern */

/*
 * The floating point values themselves are returned by static functions that
 * can be easily inlined by the compiler; similarly, the routines to test if a
 * value is the LAL failure NaN can also be inlined.
 */

/** Returns the value of the XLAL <tt>REAL4</tt> failure NaN. */
static _LAL_INLINE_ REAL4 XLALREAL4FailNaN(void);
static _LAL_INLINE_ REAL4 XLALREAL4FailNaN(void)
{
    volatile const union {
        INT4 i;
        REAL4 x;
    } val = {
    XLAL_REAL4_FAIL_NAN_INT};
    return val.x;
}

/** Returns the value of the XLAL <tt>REAL8</tt> failure NaN. */
static _LAL_INLINE_ REAL8 XLALREAL8FailNaN(void);
static _LAL_INLINE_ REAL8 XLALREAL8FailNaN(void)
{
    volatile const union {
        INT8 i;
        REAL8 x;
    } val = {
    XLAL_REAL8_FAIL_NAN_INT};
    return val.x;
}

/** Tests if a value is an XLAL <tt>REAL4</tt> failure NaN. */
static _LAL_INLINE_ int XLALIsREAL4FailNaN(REAL4 val);
static _LAL_INLINE_ int XLALIsREAL4FailNaN(REAL4 val)
{
    volatile const union {
        INT4 i;
        unsigned char s[4];
    } a = {
    XLAL_REAL4_FAIL_NAN_INT};
    volatile union {
        REAL4 x;
        unsigned char s[4];
    } b;
    size_t n;
    b.x = val;
    for (n = 0; n < sizeof(val); ++n)
        if (a.s[n] != b.s[n])
            return 0;
    return 1;
}

/** Tests if a value is an XLAL <tt>REAL8</tt> failure NaN. */
static _LAL_INLINE_ int XLALIsREAL8FailNaN(REAL8 val);
static _LAL_INLINE_ int XLALIsREAL8FailNaN(REAL8 val)
{
    volatile const union {
        INT8 i;
        unsigned char s[8];
    } a = {
    XLAL_REAL8_FAIL_NAN_INT};
    volatile union {
        REAL8 x;
        unsigned char s[8];
    } b;
    size_t n;
    b.x = val;
    for (n = 0; n < sizeof(val); ++n)
        if (a.s[n] != b.s[n])
            return 0;
    return 1;
}

/* Here are the macro constants for the fail NaNs. */
#define XLAL_REAL4_FAIL_NAN ( XLALREAL4FailNaN() ) /**< Floating-point value of the XLAL <tt>REAL4</tt> failure NaN. */
#define XLAL_REAL8_FAIL_NAN ( XLALREAL8FailNaN() ) /**< Floating-point value of the XLAL <tt>REAL8</tt> failure NaN. */

/* Here are the macros to test for fail NaNs. */
#define XLAL_IS_REAL4_FAIL_NAN(val) XLALIsREAL4FailNaN(val) /**< Tests if <tt>val</tt> is a XLAL <tt>REAL4</tt> failure NaN. */
#define XLAL_IS_REAL8_FAIL_NAN(val) XLALIsREAL8FailNaN(val) /**< Tests if <tt>val</tt> is a XLAL <tt>REAL8</tt> failure NaN. */


#endif /* SWIG */


/** XLAL error numbers and return values. */
enum XLALErrorValue {
    XLAL_SUCCESS = 0,      /**< Success return value (not an error number) */
    XLAL_FAILURE = -1,     /**< Failure return value (not an error number) */

    /* these are standard error numbers */
    XLAL_ENOENT = 2,        /**< No such file or directory */
    XLAL_EIO = 5,           /**< I/O error */
    XLAL_ENOMEM = 12,       /**< Memory allocation error */
    XLAL_EFAULT = 14,       /**< Invalid pointer */
    XLAL_EINVAL = 22,       /**< Invalid argument */
    XLAL_EDOM = 33,         /**< Input domain error */
    XLAL_ERANGE = 34,       /**< Output range error */
    XLAL_ENOSYS = 38,       /**< Function not implemented */

    /* extended error numbers start at 128 ...
     * should be beyond normal errnos */

    /* these are common errors for XLAL functions */
    XLAL_EFAILED = 128,     /**< Generic failure */
    XLAL_EBADLEN = 129,     /**< Inconsistent or invalid length */
    XLAL_ESIZE = 130,       /**< Wrong size */
    XLAL_EDIMS = 131,       /**< Wrong dimensions */
    XLAL_ETYPE = 132,       /**< Wrong or unknown type */
    XLAL_ETIME = 133,       /**< Invalid time */
    XLAL_EFREQ = 134,       /**< Invalid freqency */
    XLAL_EUNIT = 135,       /**< Invalid units */
    XLAL_ENAME = 136,       /**< Wrong name */
    XLAL_EDATA = 137,       /**< Invalid data */

    /* user-defined errors */
    XLAL_EUSR0 = 200,       /**< User-defined error 0 */
    XLAL_EUSR1 = 201,       /**< User-defined error 1 */
    XLAL_EUSR2 = 202,       /**< User-defined error 2 */
    XLAL_EUSR3 = 203,       /**< User-defined error 3 */
    XLAL_EUSR4 = 204,       /**< User-defined error 4 */
    XLAL_EUSR5 = 205,       /**< User-defined error 5 */
    XLAL_EUSR6 = 206,       /**< User-defined error 6 */
    XLAL_EUSR7 = 207,       /**< User-defined error 7 */
    XLAL_EUSR8 = 208,       /**< User-defined error 8 */
    XLAL_EUSR9 = 209,       /**< User-defined error 9 */

    /* external or internal errors */
    XLAL_ESYS = 254,        /**< System error */
    XLAL_EERR = 255,        /**< Internal error */

    /* specific mathematical and numerical errors start at 256 */

    /* IEEE floating point errors */
    XLAL_EFPINVAL = 256,      /**< IEEE Invalid floating point operation, eg sqrt(-1), 0/0 */
    XLAL_EFPDIV0 = 257,       /**< IEEE Division by zero floating point error */
    XLAL_EFPOVRFLW = 258,     /**< IEEE Floating point overflow error */
    XLAL_EFPUNDFLW = 259,     /**< IEEE Floating point underflow error */
    XLAL_EFPINEXCT = 260,     /**< IEEE Floating point inexact error */

    /* numerical algorithm errors */
    XLAL_EMAXITER = 261,      /**< Exceeded maximum number of iterations */
    XLAL_EDIVERGE = 262,      /**< Series is diverging */
    XLAL_ESING = 263,         /**< Apparent singularity detected */
    XLAL_ETOL = 264,          /**< Failed to reach specified tolerance */
    XLAL_ELOSS = 265,         /**< Loss of accuracy */

    /* failure from within a function call: "or" error number with this */
    XLAL_EFUNC = 1024         /**< Internal function call failed bit: "or" this with existing error number */
};


#ifndef SWIG    /* exclude from SWIG interface */


/*
 *
 * These functions provide message associated with an error code and print
 * an error message associated with the error code.  The macro XLAL_PERROR
 * fills in the current file and line information and uses the current
 * value of xlalErrno as the error number.
 *
 */

/** Returns the error message associated with an error number. */
const char *XLALErrorString(int errnum);

/** Prints an error message for a particular error code in a standard format. */
void XLALPerror(const char *func, const char *file, int line, int errnum);

/** Prints an error message for the current value of <tt>xlalErrno</tt>. */
#define XLAL_PERROR( ) XLALPerror(__func__, __FILE__, __LINE__, xlalErrno)


/*
 *
 * Here is the XLAL error handler type and the routines that set it.
 * Also provide is the default error handler.
 *
 */

/** The XLAL error handler type. */
typedef void XLALErrorHandlerType(const char *func, const char *file,
                                  int line, int errnum);

/** The default XLAL error handler. */
void XLALDefaultErrorHandler(const char *func, const char *file, int line,
                             int errnum);
/** A silent XLAL error handler. */
void XLALSilentErrorHandler(const char *func, const char *file, int line,
                            int errnum);

/* Other useful XLAL error handlers. */
/** The XLAL error handler that raises SIGABRT. */
void XLALAbortErrorHandler(const char *func, const char *file, int line,
                           int errnum);
/** The XLAL error handler that calls exit. */
void XLALExitErrorHandler(const char *func, const char *file, int line,
                          int errnum);
/** The XLAL error handler that prints a function call backtrace then raises SIGABRT. */
void XLALBacktraceErrorHandler(const char *func, const char *file,
                               int line, int errnum);

/** Function to return pointer to the XLAL error handler function pointer. */
XLALErrorHandlerType **XLALGetErrorHandlerPtr(void);

/** Sets the error handler to a new handler and returns the old handler. */
XLALErrorHandlerType *XLALSetErrorHandler(XLALErrorHandlerType *
                                          newHandler);

/** Sets the error handler to the default handler and returns the old handler. */
XLALErrorHandlerType *XLALSetDefaultErrorHandler(void);
/** Sets the error handler to a silent handler and returns the old handler. */
XLALErrorHandlerType *XLALSetSilentErrorHandler(void);


#endif /* SWIG */


/*
 *
 * Here are the routines that set or clear the XLAL error number.
 *
 */

#ifdef SWIG     /* SWIG interface directives */
SWIGLAL(DISABLE_EXCEPTIONS(XLALSetErrno, XLALGetBaseErrno, XLALClearErrno));
#endif /* SWIG */

/** Sets the XLAL error number to errnum, returns the new value. */
int XLALSetErrno(int errnum);

/** Gets the XLAL base error number ignoring the internal-function-failed flag. */
int XLALGetBaseErrno(void);

/** Clears the XLAL error number, returns the old value. */
int XLALClearErrno(void);


#ifndef SWIG    /* exclude from SWIG interface */


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

/** Function to return pointer to the XLAL error number. */
int *XLALGetErrnoPtr(void);

/* these are the modifiable lvalues for xlalErrno and XLALErrorHandler */
#define xlalErrno ( * XLALGetErrnoPtr() ) /**< Modifiable lvalue containing the XLAL error number */
#define XLALErrorHandler ( * XLALGetErrorHandlerPtr() ) /**< Modifiable lvalue containing the XLAL error handler */


/**
 * A macro to (i) disable the XLAL error handling and preserve the
 * current value of xlalErrno (ii) perform a statement that involves an
 * XLAL function call and (iii) restore the XLAL error handler and value of
 * xlalErrno while setting variable errnum to the xlalErrno set by the
 * statement.
 */
#define XLAL_TRY( statement, errnum ) \
	do { \
		XLALErrorHandlerType *xlalSaveErrorHandler; \
		int xlalSaveErrno; \
		xlalSaveErrorHandler = XLALSetSilentErrorHandler(); \
		xlalSaveErrno = xlalErrno; \
		XLALClearErrno(); \
		statement ; \
		errnum = xlalErrno; \
		xlalErrno = xlalSaveErrno; \
		XLALSetErrorHandler(xlalSaveErrorHandler); \
	} while (0)

/**
 * Performs the same actions as XLAL_TRY(), but additionally silences
 * any error/warning/etc. messages being printed while statement is
 * executed, regardless of the value of #lalDebugLevel.
 */
#define XLAL_TRY_SILENT( statement, errnum ) \
	do { \
		int xlalSaveDebugLevel = lalDebugLevel; \
		XLALClobberDebugLevel(xlalSaveDebugLevel & ~(LALERRORBIT | LALWARNINGBIT | LALINFOBIT | LALTRACEBIT)); \
		XLAL_TRY(statement, errnum); \
		XLALClobberDebugLevel(xlalSaveDebugLevel); \
	} while (0)

/*
 *
 * Here are the routines and macros that are used to report errors when
 * an XLAL function fails.  They (i) set the XLAL error number and (ii)
 * invoke the XLAL error handler.  The macros also (iii) return the
 * appropriate failure codes.  The macros should be used to report all
 * failures.
 *
 */

/**
 * Routine to set the XLAL error number and invoke the XLAL error handler.
 * It is used by the error macros.
 */
void XLALError(const char *func,
                          /**< name of function where the error occurs */
               const char *file,
                          /**< source file name (use the __FILE__ macro) */
               int line,  /**< source line number (use the __LINE__ macro) */
               int errnum /**< error code */
    );

/** \cond DONT_DOXYGEN */
/*
 * Helper macros for internal use only:
 * To allow for a possibly empty error message, these macros use
 *   snprintf(buf, sizeof(buf), "X" __VA_ARGS__)
 * to print any error message preceded by "X" (to silence -Wformat-zero-length)
 * to a buffer 'buf', then print the error message with XLAL_PRINT_ERROR() only
 * if 'buf' contains any characters after the "X". This construct allows for
 *   XLAL_ERROR(XLAL_EFUNC);
 *   XLAL_ERROR(XLAL_EFUNC, "%i < 0", n);
 * It does not allow a non-literal format string, e.g.
 *   const char *fmt = "%i < 0"; XLAL_ERROR(XLAL_EFUNC, fmt, n);
 * but since this is a security risk, -Wformat does not allow this anyway.
 */

#define _XLAL_ERROR_IMPL_(statement, errnum, ...) \
	do { \
		char _XLAL_ERROR_IMPL_buf_[1024]; \
		XLALStringPrint(_XLAL_ERROR_IMPL_buf_, sizeof(_XLAL_ERROR_IMPL_buf_), "X" __VA_ARGS__); \
		if (_XLAL_ERROR_IMPL_buf_[1] != 0) { \
			XLAL_PRINT_ERROR("%s", &_XLAL_ERROR_IMPL_buf_[1]); \
		} \
		XLALError(__func__, __FILE__, __LINE__, errnum); \
		statement; \
	} while (0)

#define _XLAL_CHECK_IMPL_(statement, assertion, errnum, ...) \
	do { \
		if (!(assertion)) { \
			char _XLAL_CHECK_IMPL_buf_[1024]; \
			XLALStringPrint(_XLAL_CHECK_IMPL_buf_, sizeof(_XLAL_CHECK_IMPL_buf_), "X" __VA_ARGS__); \
			if (_XLAL_CHECK_IMPL_buf_[1] != 0) { \
				XLAL_PRINT_ERROR("%s", &_XLAL_CHECK_IMPL_buf_[1]); \
			} else { \
				XLAL_PRINT_ERROR("Check failed: %s", #assertion); \
			} \
			XLALError(__func__, __FILE__, __LINE__, errnum); \
			statement; \
		} \
	} while (0)

/** \endcond */

/**
 * \brief Macro to invoke the <tt>XLALError()</tt> function and return
 * with code val (it should not really be used itself, but forms the basis for
 * other macros).
 *
 * Prototype: <b>XLAL_ERROR_VAL(val, errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b val The value to return.
 * <li> \b errnum The XLAL error number to set.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_ERROR_VAL(val, ...) _XLAL_ERROR_IMPL_(return val, __VA_ARGS__)

/**
 * Macro to invoke a failure from a XLAL routine returning an integer.
 *
 * Prototype: <b>XLAL_ERROR(errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b errnum The XLAL error number to set.
 * <li> \b fmt     (Optional) Format string for additional error information.
 * <li> \b ...     (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_ERROR(...) _XLAL_ERROR_IMPL_(return (int)XLAL_FAILURE, __VA_ARGS__)

/**
 * Macro to invoke a failure from a XLAL routine returning a pointer.
 *
 * Prototype: <b>XLAL_ERROR_NULL(errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b errnum The XLAL error number to set.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_ERROR_NULL(...) _XLAL_ERROR_IMPL_(return NULL, __VA_ARGS__)

/**
 * \brief Macro to invoke a failure from a XLAL routine returning void.
 *
 * Prototype: <b>XLAL_ERROR_VOID(errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b errnum The XLAL error number to set.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_ERROR_VOID(...) _XLAL_ERROR_IMPL_(return, __VA_ARGS__)

/**
 * \brief Macro to invoke a failure from a XLAL routine returning a <tt>REAL4</tt>.
 *
 * Prototype: <b>XLAL_ERROR_REAL4(errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b errnum The XLAL error number to set.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_ERROR_REAL4(...) _XLAL_ERROR_IMPL_(return XLAL_REAL4_FAIL_NAN, __VA_ARGS__)

/**
 * \brief Macro to invoke a failure from a XLAL routine returning a <tt>REAL8</tt>.
 *
 * Prototype <b>XLAL_ERROR_REAL8(errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b errnum The XLAL error number to set.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_ERROR_REAL8(...) _XLAL_ERROR_IMPL_(return XLAL_REAL8_FAIL_NAN, __VA_ARGS__)

/**
 * \brief Macro to invoke a failure from a C <tt>main()</tt> routine.
 *
 * Prototype <b>XLAL_ERROR_MAIN(errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b errnum The XLAL error number to set.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_ERROR_MAIN(...) _XLAL_ERROR_IMPL_(return EXIT_FAILURE, __VA_ARGS__)

/**
 * \brief Macro to invoke a failure by jumping to a <tt>XLAL_FAIL</tt> label.
 *
 * Prototype <b>XLAL_ERROR_FAIL(errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b errnum The XLAL error number to set.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_ERROR_FAIL(...) _XLAL_ERROR_IMPL_(goto XLAL_FAIL, __VA_ARGS__)

/**
 * \brief Macro to test an assertion; if it is not true, invoke the
 * <tt>XLALError()</tt> function and return with code val (it should not really
 * be used itself, but forms the basis for other macros).
 *
 * Prototype: <b>XLAL_CHECK_VAL(val, assertion, errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b val The value to return.
 * <li> \b assertion The assertion to test.
 * <li> \b errnum The XLAL error number to set if the assertion is false.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_CHECK_VAL(val, assertion, ...) _XLAL_CHECK_IMPL_(return val, assertion, __VA_ARGS__)

/**
 * \brief Macro to test an assertion and invoke a failure if it is not true
 * in a function that returns an integer.
 *
 * Prototype: <b>XLAL_CHECK(assertion, errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b assertion The assertion to test.
 * <li> \b errnum The XLAL error number to set if the assertion is false.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_CHECK(assertion, ...) _XLAL_CHECK_IMPL_(return (int)XLAL_FAILURE, assertion, __VA_ARGS__)

/**
 * \brief Macro to test an assertion and invoke a failure if it is not true
 * in a function that returns a pointer.
 *
 * Prototype: <b>XLAL_CHECK_NULL(assertion, errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b assertion The assertion to test.
 * <li> \b errnum The XLAL error number to set if the assertion is false.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_CHECK_NULL(assertion, ...) _XLAL_CHECK_IMPL_(return NULL, assertion, __VA_ARGS__)

/**
 * \brief Macro to test an assertion and invoke a failure if it is not true
 * in a function that returns void.
 *
 * Prototype: <b>XLAL_CHECK_VOID(assertion, errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b assertion The assertion to test.
 * <li> \b errnum The XLAL error number to set if the assertion is false.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_CHECK_VOID(assertion, ...) _XLAL_CHECK_IMPL_(return, assertion, __VA_ARGS__)

/**
 * \brief Macro to test an assertion and invoke a failure if it is not true
 * in a function that returns a <tt>REAL4</tt>.
 *
 * Prototype: <b>XLAL_CHECK_REAL4(assertion, errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b assertion The assertion to test.
 * <li> \b errnum The XLAL error number to set if the assertion is false.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_CHECK_REAL4(assertion, ...) _XLAL_CHECK_IMPL_(return XLAL_REAL4_FAIL_NAN, assertion, __VA_ARGS__)

/**
 * \brief Macro to test an assertion and invoke a failure if it is not true
 * in a function that returns a <tt>REAL8</tt>.
 *
 * Prototype: <b>XLAL_CHECK_REAL8(assertion, errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b assertion The assertion to test.
 * <li> \b errnum The XLAL error number to set if the assertion is false.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_CHECK_REAL8(assertion, ...) _XLAL_CHECK_IMPL_(return XLAL_REAL8_FAIL_NAN, assertion, __VA_ARGS__)

/**
 * \brief Macro to test an assertion and invoke a failure if it is not true
 * in a C <tt>main()</tt> routine.
 *
 * Prototype: <b>XLAL_CHECK_MAIN(assertion, errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b assertion The assertion to test.
 * <li> \b errnum The XLAL error number to set if the assertion is false.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_CHECK_MAIN(assertion, ...) _XLAL_CHECK_IMPL_(return EXIT_FAILURE, assertion, __VA_ARGS__)

/**
 * \brief Macro to test an assertion and invoke a failure if it is not true
 * by jumping to a <tt>XLAL_FAIL</tt> label.
 *
 * Prototype: <b>XLAL_CHECK_FAIL(assertion, errnum [, fmt [, ...]])</b>
 *
 * \b Parameters:<ul>
 * <li> \b assertion The assertion to test.
 * <li> \b errnum The XLAL error number to set if the assertion is false.
 * <li> \b fmt (Optional) Format string for additional error information.
 * <li> \b ... (Optional) Additional arguments for printf-like format.
 * </ul>
 */
#define XLAL_CHECK_FAIL(assertion, ...) _XLAL_CHECK_IMPL_(goto XLAL_FAIL, assertion, __VA_ARGS__)

/**
 * \brief Macro to test an assertion and invoke a failure if it is not true
 * by calling <tt>lalAbortHook()</tt>.
 *
 * Prototype: <b>XLAL_CHECK_ABORT(assertion)</b>
 *
 * \b Parameters:<ul>
 * <li> \b assertion The assertion to test.
 * </ul>
 */
#define XLAL_CHECK_ABORT(assertion) \
	do { \
		if (!(assertion)) { \
			XLAL_PRINT_ERROR("Check failed: %s", #assertion); \
			lalAbortHook("XLAL_CHECK_ABORT() failed"); \
		} \
	} while (0)

/**
 * \brief Macro to test an assertion and invoke a failure if it is not true
 * by calling <tt>exit(1)</tt>.
 *
 * Prototype: <b>XLAL_CHECK_EXIT(assertion)</b>
 *
 * \b Parameters:<ul>
 * <li> \b assertion The assertion to test.
 * </ul>
 */
#define XLAL_CHECK_EXIT(assertion) \
	do { \
		if (!(assertion)) { \
			XLAL_PRINT_ERROR("Check failed: %s", #assertion); \
			exit(1); \
		} \
	} while (0)


#endif /* SWIG */


/** @} */

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* XLALERROR_H */
