/*
*  Copyright (C) 2007 Jolien Creighton
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

/* ---------- SEE LALError.dox for doxygen documentation ---------- */

#ifndef _LALERROR_H
#define _LALERROR_H

#include <lal/LALDebugLevel.h>

#include <lal/LALDatatypes.h>

#include <lal/XLALError.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

#ifndef SWIG    /* exclude from SWIG interface */

extern int (*lalRaiseHook) (int, const char *, ...);
extern void (*lalAbortHook) (const char *, ...);

/** \addtogroup LALError_h */ /*@{ */

int LALPrintError(const char *fmt, ...);

int LALRaise(int sig, const char *fmt, ...);

void LALAbort(const char *fmt, ...);

int LALError(LALStatus * status, const char *statement);

int LALWarning(LALStatus * status, const char *warning);

int LALInfo(LALStatus * status, const char *info);

int LALTrace(LALStatus * status, int exitflg);

/*@}*/

int
LALInitStatus(LALStatus * status, const char *function,
              const char *id, const char *file, const int line);

int LALPrepareReturn(LALStatus * status, const char *file, const int line);

int
LALAttatchStatusPtr(LALStatus * status, const char *file, const int line);

int
LALDetatchStatusPtr(LALStatus * status, const char *file, const int line);

int
LALPrepareAbort(LALStatus * status, const INT4 code, const char *mesg,
                const char *file, const int line);

int
LALPrepareAssertFail(LALStatus * status, const INT4 code,
                     const char *mesg, const char *statement,
                     const char *file, const int line);

int
LALCheckStatusPtr(LALStatus * status, const char *statement,
                  const char *file, const int line);

#ifdef NOLALMACROS

void FREESTATUSPTR(LALStatus * status);

void REPORTSTATUS(LALStatus * status);

#endif

#ifdef NDEBUG

#define LALError( statusptr, statement ) 0
#define LALWarning( statusptr, warning ) 0
#define LALInfo( statusptr, info )       0
#define LALTrace( statusptr, exitflg )   0

#else

#ifndef NOLALMACROS

#define LALError( statusptr, statement )                                    \
  ( lalDebugLevel & LALERROR ?                                                 \
    LALPrintError( "Error[%d] %d: function %s, file %s, line %d, %s\n"       \
        "        %s %s\n", (statusptr)->level, (statusptr)->statusCode,      \
        (statusptr)->function, (statusptr)->file, (statusptr)->line,        \
        (statusptr)->Id, (statement),                    \
        (statusptr)->statusDescription ) : 0 )

#define LALWarning( statusptr, warning )                                    \
  ( lalDebugLevel & LALWARNING ?                                               \
    LALPrintError( "Warning[%d]: function %s, file %s, line %d, %s\n"        \
        "        %s\n", (statusptr)->level, (statusptr)->function,          \
        (statusptr)->file, (statusptr)->line, (statusptr)->Id, (warning) )  \
      : 0 )

#define LALInfo( statusptr, info )                                          \
  ( lalDebugLevel & LALINFO ?                                                  \
    LALPrintError( "Info[%d]: function %s, file %s, line %d, %s\n"          \
        "        %s\n", (statusptr)->level, (statusptr)->function,          \
        (statusptr)->file, (statusptr)->line, (statusptr)->Id, (info) )     \
      : 0 )

#define LALTrace( statusptr, exitflg ) \
  ( lalDebugLevel & LALTRACE ? \
    LALPrintError( "%s[%d]: function %s, file %s, line %d, %s\n",      \
        (exitflg) ? "Leave" : "Enter", (statusptr)->level, \
        (statusptr)->function, (statusptr)->file, (statusptr)->line, \
        (statusptr)->Id )     \
      : 0 )

#endif /* NOLALMACROS */

#endif /* NDEBUG */

#endif /* SWIG */

/*
 * Error codes and corresponding error messages.
 */

#define LAL_FAIL_ERR	XLAL_EFAILED
#define LAL_FAIL_MSG	"operation failed"
#define LAL_NULL_ERR	XLAL_EFAULT
#define LAL_NULL_MSG	"unexpected NULL pointer"
#define LAL_NNULL_ERR	XLAL_EFAULT
#define LAL_NNULL_MSG	"unexpected non-NULL pointer"
#define LAL_NOMEM_ERR	XLAL_ENOMEM
#define LAL_NOMEM_MSG	"out of memory"
#define LAL_RANGE_ERR	XLAL_ERANGE
#define LAL_RANGE_MSG	"parameter out of range"
#define LAL_BADPARM_ERR XLAL_EINVAL
#define LAL_BADPARM_MSG "invalid parameter value"

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALERROR_H */
