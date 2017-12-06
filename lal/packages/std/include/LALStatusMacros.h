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

/* ---------- SEE LALStatusMacros.dox for doxygen documentation ---------- */

#ifndef _LALSTATUSMACROS_H
#define _LALSTATUSMACROS_H

#include <lal/LALConfig.h>
#ifdef NDEBUG
#ifndef LAL_NDEBUG
#define LAL_NDEBUG
#endif
#endif

#include <stdlib.h>
#include <string.h>

#include <lal/LALDebugLevel.h>
#include <lal/LALMalloc.h>
#include <lal/LALDatatypes.h>
#include <lal/LALError.h>

#ifdef  __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

extern const int lalNoDebug;

#define LAL_EXLAL     16384
#define LAL_MSGEXLAL  "Failure in an XLAL routine"
#define ABORTXLAL(sp) ABORT(sp,LAL_EXLAL,LAL_MSGEXLAL)

#ifndef NOLALMACROS

#define INITSTATUS( statusptr )                                               \
  do { if ( (statusptr) )                                                     \
  {                                                                           \
    INT4 level_ = (statusptr)->level ;                                        \
    INT4 statp_ = (statusptr)->statusPtr ? 1 : 0 ;                            \
    memset( (statusptr), 0, sizeof( LALStatus ) ); /* possible memory leak */ \
    (statusptr)->level    = level_ > 0 ? level_ : 1 ;                         \
    (statusptr)->Id       = "$Id$";                                           \
    (statusptr)->function = __func__;                                         \
    SETSTATUSFILELINE( statusptr );                                           \
    (void) LALTrace( statusptr, 0 );                                          \
    if ( statp_ )                                                             \
      ABORT( statusptr, -2, "INITSTATUS: non-null status pointer" );          \
    else if ( xlalErrno )                                                     \
      ABORT( statusptr, -16, "INITSTATUS: non-zero xlalErrno" );              \
  }                                                                           \
  else                                                                        \
    lalAbortHook( "Abort: function %s, file %s, line %d, %s\n"                \
                  "       Null status pointer passed to function\n",          \
                  __func__, __FILE__, __LINE__, "$Id$" );                     \
  } while ( 0 )

#define RETURN( statusptr )                                                   \
  do { if ( 1 )                                                               \
  {                                                                           \
    SETSTATUSFILELINE( statusptr );                                           \
    if ( (statusptr)->statusCode )                                            \
      (void) LALError( statusptr, "RETURN:" );                                \
    (void) LALTrace( statusptr, 1 );                                          \
    if ( xlalErrno )                                                          \
      ABORT( statusptr, -32, "RETURN: untrapped XLAL error" );                \
    return;                                                                   \
  }                                                                           \
  } while ( 0 )

#define ATTATCHSTATUSPTR(statusptr)                                           \
  do { if ( !(statusptr)->statusPtr )                                         \
  {                                                                           \
    (statusptr)->statusPtr = (LALStatus *)LALCalloc( 1, sizeof( LALStatus ) );\
    if ( !(statusptr)->statusPtr )                                            \
      ABORT( statusptr, -4, "ATTATCHSTATUSPTR: memory allocation error" );    \
    (statusptr)->statusPtr->level = (statusptr)->level + 1;                   \
  }                                                                           \
  else                                                                        \
    ABORT( statusptr, -2, "ATTATCHSTATUSPTR: non-null status pointer" );      \
  } while ( 0 )

#define DETATCHSTATUSPTR( statusptr )                                         \
  do { if ( (statusptr)->statusPtr )                                          \
  {                                                                           \
    FREESTATUSPTR( statusptr );                                               \
    (statusptr)->statusCode = 0;                                              \
    (statusptr)->statusDescription = NULL;                                    \
  }                                                                           \
  else                                                                        \
    ABORT( statusptr, -8, "DETATCHSTATUSPTR: null status pointer" );          \
  } while ( 0 )

#define ABORT( statusptr, code, mesg )                                        \
  do { if ( 1 )                                                               \
  {                                                                           \
    if ( (statusptr)->statusPtr ) FREESTATUSPTR( statusptr );                   \
    SETSTATUS( statusptr, code, mesg );                                       \
    if ( code )                                                               \
      (void) LALError( statusptr, "ABORT:" );                                 \
    (void) LALTrace( statusptr, 1 );                                          \
    return;                                                                   \
  }                                                                           \
  } while ( 0 )

#ifdef LAL_NDEBUG
#define ASSERT( assertion, statusptr, code, mesg )
#else
#define ASSERT( assertion, statusptr, code, mesg )                            \
  do { if ( !(assertion) )                                                    \
  {                                                                           \
    if ( (statusptr)->statusPtr )                                               \
      FREESTATUSPTR( statusptr );                                             \
    SETSTATUS( statusptr, code, mesg );                                       \
    (void) LALError( statusptr, "Assertion \"" #assertion "\" failed:" );     \
    (void) LALTrace( statusptr, 1 );                                          \
    return;                                                                   \
  }                                                                           \
  } while ( 0 )
#endif

#define TRY( func, statusptr )                                                \
  do { if ( (func), (statusptr)->statusPtr->statusCode )                      \
  {                                                                           \
    SETSTATUS( statusptr, -1, "Recursive error" );                            \
    (void) LALError( statusptr, "Function call \"" #func "\" failed:" );      \
    (void) LALTrace( statusptr, 1 );                                          \
    return;                                                                   \
  }                                                                           \
  } while ( 0 )

#define CHECKSTATUSPTR( statusptr )                                           \
  do { if ( (statusptr)->statusPtr->statusCode )                              \
  {                                                                           \
    SETSTATUS( statusptr, -1, "Recursive error" );                            \
    (void) LALError( statusptr, "CHECKSTATUSPTR:" );                          \
    (void) LALTrace( statusptr, 1 );                                          \
    return;                                                                   \
  }                                                                           \
  } while ( 0 )

#define FREESTATUSPTR( statusptr )                                            \
  do                                                                          \
  {                                                                           \
    LALStatus *next_ = (statusptr)->statusPtr->statusPtr;                     \
    LALFree( (statusptr)->statusPtr );                                        \
    (statusptr)->statusPtr = next_;                                           \
  }                                                                           \
  while ( (statusptr)->statusPtr )

#define REPORTSTATUS( statusptr )                                             \
  do                                                                          \
  {                                                                           \
    LALStatus *ptr_;                                                          \
    for ( ptr_ = (statusptr); ptr_; ptr_ = ptr_->statusPtr )                  \
    {                                                                         \
      LALPrintError( "\nLevel %i: %s\n", ptr_->level, ptr_->Id );             \
      if ( ptr_->statusCode )                                                 \
        LALPrintError( "\tStatus code %i: %s\n", ptr_->statusCode,            \
                       ptr_->statusDescription );                             \
      else                                                                    \
        LALPrintError( "\tStatus code 0: Nominal\n" );                        \
      LALPrintError( "\tfunction %s, file %s, line %i\n",                     \
                     ptr_->function, ptr_->file, ptr_->line );                \
    }                                                                         \
  } while ( 0 )

#else /* NOLALMACROS */

#define INITSTATUS( statusptr ) \
  do { if ( LALInitStatus( statusptr, __func__, "$Id$", __FILE__, __LINE__ ) ) return; } while ( 0 )

#define RETURN( statusptr ) \
  do { if ( LALPrepareReturn( statusptr, __FILE__, __LINE__ ), 1 ) return; } while ( 0 )

#define ATTATCHSTATUSPTR( statusptr ) \
  do { if ( LALAttatchStatusPtr( statusptr, __FILE__, __LINE__ ) ) return; } while ( 0 )

#define DETATCHSTATUSPTR( statusptr ) \
  do { if ( LALDetatchStatusPtr( statusptr, __FILE__, __LINE__ ) ) return; } while ( 0 )

#define ABORT( statusptr, code, mesg ) \
  do { if ( LALPrepareAbort( statusptr, code, mesg, __FILE__, __LINE__ ), 1 ) return; } while ( 0 )

#ifdef LAL_NDEBUG
#define ASSERT( assertion, statusptr, code, mesg )
#else
#define ASSERT( assertion, statusptr, code, mesg )                            \
  do { if ( !(assertion) )                                                    \
  {                                                                           \
    LALPrepareAssertFail( statusptr, code, mesg,                              \
                          "Assertion \"" #assertion "\" failed:",             \
                          __FILE__, __LINE__ );                               \
    return;                                                                   \
  }                                                                           \
  } while ( 0 )
#endif

#define TRY( func, statusptr )                                                \
  do                                                                          \
  {                                                                           \
    (func);                                                                   \
    if ( LALCheckStatusPtr( statusptr, "Function call \"" #func "\" failed:", \
                            __FILE__, __LINE__ ) )                            \
      return;                                                                 \
  }                                                                           \
  while ( 0 )

#define CHECKSTATUSPTR( statusptr )                                           \
  do { if ( LALCheckStatusPtr( statusptr, "CHECKSTATUSPTR:", __FILE__, __LINE__ ) ) return; } while ( 0 )

#endif /* NOLALMACROS */

/* these just have to be macros... */

#define BEGINFAIL( statusptr )                                                \
do {                                                                          \
  if ( !(statusptr) )                                                         \
    ABORT( statusptr, -8, "BEGINFAIL: null status pointer" );                 \
  if ( !( (statusptr)->statusPtr ) )                                          \
    ABORT( statusptr, -8, "BEGINFAIL: null status pointer pointer" );         \
  if ( (statusptr)->statusPtr->statusCode ) {                                 \
    LALStatus *statusPtrSave_ = (statusptr)->statusPtr;                       \
    (statusptr)->statusPtr = NULL;                                            \
    ATTATCHSTATUSPTR( statusptr );                                            \
    do

#define ENDFAIL( statusptr )                                                  \
    while ( 0 );                                                              \
    DETATCHSTATUSPTR( statusptr );                                            \
    (statusptr)->statusPtr = statusPtrSave_;                                  \
    SETSTATUS( statusptr, -1, "Recursive error" );                            \
    (void) LALError( statusptr, "ENDFAIL:" );                                 \
    (void) LALTrace( statusptr, 1 );                                          \
    return;                                                                   \
  }                                                                           \
} while ( 0 )

#define SETSTATUSFILELINE( statusptr ) \
  ( ( void ) ( (statusptr)->file = __FILE__, (statusptr)->line = __LINE__ ) )

#define SETSTATUS( statusptr, code, mesg )                                    \
  ( SETSTATUSFILELINE( statusptr ),                                           \
    (statusptr)->statusDescription = (mesg),                                  \
    (statusptr)->statusCode = (code) )


#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /* _LALSTATUSMACROS_H */
