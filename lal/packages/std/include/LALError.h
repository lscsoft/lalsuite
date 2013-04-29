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

#include <lal/LALDatatypes.h>

#include <lal/XLALError.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

/* lalDebugLevel bit field values: */
enum {
    LALERRORBIT   = 0001,
    LALWARNINGBIT = 0002,
    LALINFOBIT    = 0004,
    LALTRACEBIT   = 0010,
    LALMEMDBGBIT  = 0020,
    LALMEMPADBIT  = 0040,
    LALMEMTRKBIT  = 0100,
    LALMEMINFOBIT = 0200
};

/* composite lalDebugLevels: */
enum {
    LALNDEBUG   = 0,
    LALERROR    = LALERRORBIT,
    LALWARNING  = LALWARNINGBIT,
    LALINFO     = LALINFOBIT,
    LALTRACE    = LALTRACEBIT,
    LALMSGLVL1  = LALERRORBIT,
    LALMSGLVL2  = LALERRORBIT | LALWARNINGBIT,
    LALMSGLVL3  = LALERRORBIT | LALWARNINGBIT | LALINFOBIT,
    LALMEMDBG   = LALMEMDBGBIT | LALMEMPADBIT | LALMEMTRKBIT,
    LALMEMTRACE = LALTRACEBIT | LALMEMDBG | LALMEMINFOBIT,
    LALALLDBG   = ~LALNDEBUG
};

#ifndef SWIG    /* exclude from SWIG interface */

extern int (*lalRaiseHook) (int, const char *, ...);
extern void (*lalAbortHook) (const char *, ...);

/** \addtogroup LALError_h *//*@{ */

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

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALERROR_H */
