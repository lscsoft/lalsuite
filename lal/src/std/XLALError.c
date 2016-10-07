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
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* - NOTE: API is doxygen-documented in header file XLALError.h - */

#include <config.h>

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#ifdef HAVE_EXECINFO_H
#include <execinfo.h>
#define BACKTRACE_LEVELMAX 0100
#endif

#include <lal/XLALError.h>
#include <lal/LALStdlib.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/*
 *
 * Routines to print generic error messages and warning messages.
 *
 */

/* Prints an error message if error printing is enabled by lalDebugLevel. */
int XLALVPrintError(const char *fmt, va_list ap)
{
    return (lalDebugLevel & LALERROR) ? vfprintf(stderr, fmt, ap) : 0;
}

/* Prints a warning message if warning printing is enabled by lalDebugLevel. */
int XLALVPrintWarning(const char *fmt, va_list ap)
{
    return (lalDebugLevel & LALWARNING) ? vfprintf(stderr, fmt, ap) : 0;
}

/* Prints an info message if info printing is enabled by lalDebugLevel. */
int XLALVPrintInfo(const char *fmt, va_list ap)
{
    return (lalDebugLevel & LALINFO) ? vfprintf(stderr, fmt, ap) : 0;
}

/* Prints an error message if error printing is enabled by lalDebugLevel. */
int XLALPrintError(const char *fmt, ...)
{
    int n = 0;
    va_list ap;
    va_start(ap, fmt);
    n = XLALVPrintError(fmt, ap);
    va_end(ap);
    return n;
}

/* Prints a warning message if warning printing is enabled by lalDebugLevel. */
int XLALPrintWarning(const char *fmt, ...)
{
    int n = 0;
    va_list ap;
    va_start(ap, fmt);
    n = XLALVPrintWarning(fmt, ap);
    va_end(ap);
    return n;
}

/* Prints an info message if info printing is enabled by lalDebugLevel. */
int XLALPrintInfo(const char *fmt, ...)
{
    int n = 0;
    va_list ap;
    va_start(ap, fmt);
    n = XLALVPrintInfo(fmt, ap);
    va_end(ap);
    return n;
}

/*
 * Prints a standard-formatted error message
 * (if error printing is enabled by lalDebugLevel).
 */
void XLALVPrintErrorMessage(const char *func, const char *file, int line,
                            const char *fmt, va_list ap)
{
    XLALPrintError("XLAL Error");
    if (func && *func)
        XLALPrintError(" - %s", func);
    if (file && *file)
        XLALPrintError(" (%s:%d)", file, line);
    XLALPrintError(": ");
    XLALVPrintError(fmt, ap);
    XLALPrintError("\n");
    return;
}

void XLALVPrintWarningMessage(const char *func, const char *file, int line,
                              const char *fmt, va_list ap)
{
    XLALPrintWarning("XLAL Warning");
    if (func && *func)
        XLALPrintWarning(" - %s", func);
    if (file && *file)
        XLALPrintWarning(" (%s:%d)", file, line);
    XLALPrintWarning(": ");
    XLALVPrintWarning(fmt, ap);
    XLALPrintWarning("\n");
    return;
}

void XLALVPrintInfoMessage(const char *func, const char *file, int line,
                           const char *fmt, va_list ap)
{
    XLALPrintInfo("XLAL Info");
    if (func && *func)
        XLALPrintInfo(" - %s", func);
    if (file && *file)
        XLALPrintInfo(" (%s:%d)", file, line);
    XLALPrintInfo(": ");
    XLALVPrintInfo(fmt, ap);
    XLALPrintInfo("\n");
    return;
}

void XLALPrintErrorMessage(const char *func, const char *file, int line,
                           const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    XLALVPrintErrorMessage(func, file, line, fmt, ap);
    va_end(ap);
    return;
}

void XLALPrintWarningMessage(const char *func, const char *file, int line,
                             const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    XLALVPrintWarningMessage(func, file, line, fmt, ap);
    va_end(ap);
    return;
}

void XLALPrintInfoMessage(const char *func, const char *file, int line,
                          const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    XLALVPrintInfoMessage(func, file, line, fmt, ap);
    va_end(ap);
    return;
}



/*
 * Prints a progress bar at the "info" verbosity level.
 */
int XLALPrintProgressBar(double fraction)
{
    static const char mrk[] =
        "+++++++++++++++++++++++++++++++++++++++++++++++++)";
    static const char spc[] =
        "-------------------------------------------------)";
    int l = XLAL_NUM_ELEM(mrk) - 1;
    int offset =
        floor((fraction < 0.0 ? 0.0 : fraction >
               1.0 ? 1.0 : fraction) * l + 0.5);

    return XLALPrintInfo("[%s%s %.1f%%", mrk + l - offset, spc + offset,
                         100.0 * fraction);
}

/*
 *
 * Implementation of xlalErrno and XLALErrorHandler.
 * If code must be POSIX thread safe then the code is somewhat more complicated.
 *
 */

#ifndef LAL_PTHREAD_LOCK        /* non-pthread-safe code */

/* XLAL error number is just a global variable */
int xlalErrnoGlobal = 0;

/* XLALGetErrnoPtr just returns the address of the global variable */
int *XLALGetErrnoPtr(void)
{
    return &xlalErrnoGlobal;
}

/* XLAL error handler is just a global variable */
XLALErrorHandlerType *xlalErrorHandlerGlobal = NULL;

/* XLALGetErrorHandlerPtr just returns the address of the global variable */
XLALErrorHandlerType **XLALGetErrorHandlerPtr(void)
{
    return &xlalErrorHandlerGlobal;
}

#else /* pthread safe code */

/* Note: malloc and free are used here rather than LALMalloc and LALFree...
 * this is so that if a user checks for memory leaks within a thread before
 * rejoining to the main thread (which shouldn't be done) then at least
 * these routines won't report any leaks.  */

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <pthread.h>

pthread_key_t xlalErrnoKey;
pthread_once_t xlalErrnoKeyOnce = PTHREAD_ONCE_INIT;
pthread_key_t xlalErrorHandlerKey;
pthread_once_t xlalErrorHandlerKeyOnce = PTHREAD_ONCE_INIT;

/* routine to free the XLAL error number pointer */
static void XLALDestroyErrnoPtr(void *xlalErrnoPtr)
{
    free(xlalErrnoPtr);
    return;
}

/* routine to free the XLAL error handler pointer */
static void XLALDestroyErrorHandlerPtr(void *xlalErrorHandlerPtr)
{
    free(xlalErrorHandlerPtr);
    return;
}

/* routine to create the XLAL error number key */
static void XLALCreateErrnoKey(void)
{
    pthread_key_create(&xlalErrnoKey, XLALDestroyErrnoPtr);
    return;
}

/* routine to create the XLAL error handler key */
static void XLALCreateErrorHandlerKey(void)
{
    pthread_key_create(&xlalErrorHandlerKey, XLALDestroyErrorHandlerPtr);
    return;
}

/* return the pointer to the XLAL error number in this thread */
int *XLALGetErrnoPtr(void)
{
    int *xlalErrnoPtr;

    /* create key on the first call only */
    pthread_once(&xlalErrnoKeyOnce, XLALCreateErrnoKey);

    /* get the pointer to the XLAL error number in this thread */
    xlalErrnoPtr = pthread_getspecific(xlalErrnoKey);
    if (!xlalErrnoPtr) {        /* haven't allocated pointer yet... do it now */
        xlalErrnoPtr = malloc(sizeof(*xlalErrnoPtr));
        if (!xlalErrnoPtr)
            lalAbortHook
                ("could not set xlal error number: malloc failed\n");
        *xlalErrnoPtr = 0;      /* raises segv if memory allocation fails */
        /* now set the value of the pointer in this thread in the key */
        if (pthread_setspecific(xlalErrnoKey, xlalErrnoPtr))
            lalAbortHook
                ("could not set xlal error number: pthread_setspecific failed\n");
    }
    return xlalErrnoPtr;
}

/* return the pointer to the XLAL error handler in this thread */
XLALErrorHandlerType **XLALGetErrorHandlerPtr(void)
{
    XLALErrorHandlerType **xlalErrorHandlerPtr;

    /* create key on the first call only */
    pthread_once(&xlalErrorHandlerKeyOnce, XLALCreateErrorHandlerKey);

    /* get the pointer to the XLAL error handler in this thread */
    xlalErrorHandlerPtr = pthread_getspecific(xlalErrorHandlerKey);
    if (!xlalErrorHandlerPtr) { /* haven't allocated pointer yet... do it now */
        xlalErrorHandlerPtr = malloc(sizeof(*xlalErrorHandlerPtr));
        if (!xlalErrorHandlerPtr)
            lalAbortHook
                ("could not set xlal error handler: malloc failed\n");
        *xlalErrorHandlerPtr = NULL;    /* raises segv if memory allocation fails */
        /* now set the value of the pointer in this thread in the key */
        if (pthread_setspecific(xlalErrorHandlerKey, xlalErrorHandlerPtr))
            lalAbortHook
                ("could not set xlal error handler: pthread_setspecific failed\n");
    }
    return xlalErrorHandlerPtr;
}

#endif /* end of pthread-safe code */


/*
 *
 * Here are the routines to set the error number or error handler.
 *
 */


/* Set the XLAL error number to errnum. */
int XLALSetErrno(int errnum)
{
    if (errnum == 0) {
        xlalErrno = 0;
        return xlalErrno;
    }

    /*
     * if this is an error indicating an internal error then set the bit
     * that indicates this; otherwise, xlalErrno should presumably be zero
     */
    if (errnum & XLAL_EFUNC) {
        xlalErrno |= XLAL_EFUNC;        /* make sure XLAL_EFUNC bit is set */
        return xlalErrno;
    }

    /*
     * if xlalErrno is not zero, probably forgot to deal with previous
     * error
     */
    if (xlalErrno)
        XLAL_PRINT_WARNING("Ignoring previous error (xlalErrno=%d) %s\n",
                           xlalErrno, XLALErrorString(xlalErrno));
    xlalErrno = errnum;
    return xlalErrno;
}


/* Gets the basic error number ignoring the internal-function-failed flag. */
int XLALGetBaseErrno(void)
{
    return xlalErrno & ~XLAL_EFUNC;
}


/* Clears the XLAL error number. */
int XLALClearErrno(void)
{
    int olderrno = xlalErrno;
    xlalErrno = 0;
    return olderrno;
}


/* Set the XLAL error handler to newHandler; return the old handler. */
XLALErrorHandlerType *XLALSetErrorHandler(XLALErrorHandlerType *
                                          newHandler)
{
    XLALErrorHandlerType *oldHandler;
    oldHandler = XLALErrorHandler;
    XLALErrorHandler = newHandler;
    return oldHandler;
}


/* Set the XLAL error handler to the default handler; return the old handler.  */
XLALErrorHandlerType *XLALSetDefaultErrorHandler(void)
{
    XLALErrorHandlerType *oldHandler;
    oldHandler = XLALErrorHandler;
    XLALErrorHandler = XLALDefaultErrorHandler;
    return oldHandler;
}

/* Set the XLAL error handler to a silent handler; return the old handler. */
XLALErrorHandlerType *XLALSetSilentErrorHandler(void)
{
    XLALErrorHandlerType *oldHandler;
    oldHandler = XLALErrorHandler;
    XLALErrorHandler = XLALSilentErrorHandler;
    return oldHandler;
}


/*
 *
 * Routines to give the error message associated with a given error number.
 *
 */


/* Return the error message associated with an error number or return value. */
const char *XLALErrorString(int code)
{

    if (code <= 0) {    /* this is a return code, not an error number */
        if (code == 0)
            return "Success";
        else if (code == -1)
            return "Failure";
        else
            return "Unknown return code";
    }

    /* check to see if an internal function call has failed, but the error
     * number was not "or"ed against the mask XLAL_EFUNC */
    if (code == XLAL_EFUNC)
        return "Internal function call failed";

    /* use this to report error strings... deals with possible mask for
     * errors arising from internal function calls */
# define XLAL_ERROR_STRING(s) \
    ( ( code & XLAL_EFUNC ) ? "Internal function call failed: " s : (const char *) s )
    switch (code & ~XLAL_EFUNC) {
        /* these are standard error numbers */
    case XLAL_ENOENT:
        return XLAL_ERROR_STRING("No such file or directory");
    case XLAL_EIO:
        return XLAL_ERROR_STRING("I/O error");
    case XLAL_ENOMEM:
        return XLAL_ERROR_STRING("Memory allocation error");
    case XLAL_EFAULT:
        return XLAL_ERROR_STRING("Invalid pointer");
    case XLAL_EINVAL:
        return XLAL_ERROR_STRING("Invalid argument");
    case XLAL_EDOM:
        return XLAL_ERROR_STRING("Input domain error");
    case XLAL_ERANGE:
        return XLAL_ERROR_STRING("Output range error");
    case XLAL_ENOSYS:
        return XLAL_ERROR_STRING("Function not implemented");

        /* extended error numbers start at 128 ...
         * should be beyond normal errnos */

        /* these are common errors for XLAL functions */
    case XLAL_EFAILED:
        return XLAL_ERROR_STRING("Generic failure");
    case XLAL_EBADLEN:
        return XLAL_ERROR_STRING("Inconsistent or invalid vector length");
    case XLAL_ESIZE:
        return XLAL_ERROR_STRING("Wrong size");
    case XLAL_EDIMS:
        return XLAL_ERROR_STRING("Wrong dimensions");
    case XLAL_ETYPE:
        return XLAL_ERROR_STRING("Wrong or unknown type");
    case XLAL_ETIME:
        return XLAL_ERROR_STRING("Invalid time");
    case XLAL_EFREQ:
        return XLAL_ERROR_STRING("Invalid freqency");
    case XLAL_EUNIT:
        return XLAL_ERROR_STRING("Invalid units");
    case XLAL_ENAME:
        return XLAL_ERROR_STRING("Wrong name");
    case XLAL_EDATA:
        return XLAL_ERROR_STRING("Invalid data");

        /* user-defined errors */
    case XLAL_EUSR0:
        return XLAL_ERROR_STRING("User-defined error 0");
    case XLAL_EUSR1:
        return XLAL_ERROR_STRING("User-defined error 1");
    case XLAL_EUSR2:
        return XLAL_ERROR_STRING("User-defined error 2");
    case XLAL_EUSR3:
        return XLAL_ERROR_STRING("User-defined error 3");
    case XLAL_EUSR4:
        return XLAL_ERROR_STRING("User-defined error 4");
    case XLAL_EUSR5:
        return XLAL_ERROR_STRING("User-defined error 5");
    case XLAL_EUSR6:
        return XLAL_ERROR_STRING("User-defined error 6");
    case XLAL_EUSR7:
        return XLAL_ERROR_STRING("User-defined error 7");
    case XLAL_EUSR8:
        return XLAL_ERROR_STRING("User-defined error 8");
    case XLAL_EUSR9:
        return XLAL_ERROR_STRING("User-defined error 9");

        /* external or internal errors */
    case XLAL_ESYS:
        return XLAL_ERROR_STRING("System error");
    case XLAL_EERR:
        return XLAL_ERROR_STRING("Internal error");

        /* specific mathematical and numerical errors start at 256 */

        /* IEEE floating point errors */
    case XLAL_EFPINVAL:
        return
            XLAL_ERROR_STRING
            ("Invalid floating point operation, eg sqrt(-1), 0/0");
    case XLAL_EFPDIV0:
        return XLAL_ERROR_STRING("Division by zero floating point error");
    case XLAL_EFPOVRFLW:
        return XLAL_ERROR_STRING("Floating point overflow error");
    case XLAL_EFPUNDFLW:
        return XLAL_ERROR_STRING("Floating point underflow error");
    case XLAL_EFPINEXCT:
        return XLAL_ERROR_STRING("Floating point inexact error");

        /* numerical algorithm errors */
    case XLAL_EMAXITER:
        return XLAL_ERROR_STRING("Exceeded maximum number of iterations");
    case XLAL_EDIVERGE:
        return XLAL_ERROR_STRING("Series is diverging");
    case XLAL_ESING:
        return XLAL_ERROR_STRING("Apparent singularity detected");
    case XLAL_ETOL:
        return XLAL_ERROR_STRING("Failed to reach specified tolerance");
    case XLAL_ELOSS:
        return XLAL_ERROR_STRING("Loss of accuracy");

        /* unrecognized error number */
    default:
        return "Unknown error";
    }
# undef XLAL_ERROR_STRING
    return NULL;        /* impossible to get here */
}

/* Print an error message associated with an error number or return code. */
void XLALPerror(const char *func, const char *file, int line, int code)
{
    if (code > 0)
        XLALPrintError("XLAL Error");
    else
        XLALPrintError("XLAL Result");
    if (func && *func)
        XLALPrintError(" - %s", func);
    if (file && *file)
        XLALPrintError(" (%s:%d)", file, line);
    XLALPrintError(": %s\n", XLALErrorString(code));
    return;
}


/*
 *
 * Here is the default error handler.
 *
 */

/* Default XLAL error handler */
void XLALDefaultErrorHandler(const char *func, const char *file, int line,
                             int errnum)
{
    XLALPerror(func, file, line, errnum);
    return;
}

/* Silent XLAL error handler */
void XLALSilentErrorHandler(const char UNUSED * func,
                            const char UNUSED * file, int UNUSED line,
                            int UNUSED errnum)
{
    return;
}


/*
 *
 * Routine to set the error number and invoke the current error handler.
 *
 */
void XLALError(const char *func, const char *file, int line, int errnum)
{
    XLALSetErrno(errnum);
    if (!XLALErrorHandler)
        XLALErrorHandler = XLALDefaultErrorHandler;
    XLALErrorHandler(func, file, line, xlalErrno);
    return;
}


/*
 *
 * Useful standard error handlers
 *
 */

/* XLAL error handler to abort on error. */
void XLALAbortErrorHandler(const char *func, const char *file, int line,
                           int errnum)
{
    XLALPerror(func, file, line, errnum);
    abort();
}

/* XLAL error handler to exit on error. */
void XLALExitErrorHandler(const char *func, const char *file, int line,
                          int errnum)
{
    XLALPerror(func, file, line, errnum);
    exit(1);
}

/* XLAL error handler to abort on error and print a backtrace (if possible). */
void XLALBacktraceErrorHandler(const char *func, const char *file,
                               int line, int errnum)
{
    XLALPerror(func, file, line, errnum);
#if defined(HAVE_BACKTRACE) && defined(BACKTRACE_LEVELMAX)
    void *callstack[BACKTRACE_LEVELMAX];
    size_t frames = backtrace(callstack, BACKTRACE_LEVELMAX);
    fprintf(stderr, "backtrace:\n");
    backtrace_symbols_fd(callstack, frames, fileno(stderr));
#endif
    abort();
}
