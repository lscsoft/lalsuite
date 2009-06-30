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

/**** <lalVerbatim file="ExampleCV">
 * Author: Al A. Lal
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \subsection{Module \texttt{Example.c}}
 *
 * %% A one-line description of the function(s) defined in this module.
 *
 * \subsubsection*{Prototypes}
 * \input{ExampleCP}
 * %% \idx{LALExample()}
 *
 * \subsubsection*{Description}
 *
 * %% A description of the data analysis task performed by this function;
 * %% this is the main place to document the module.
 *
 * \subsubsection*{Algorithm}
 *
 * %% A description of the method used to perform the calculation.
 *
 * \subsubsection*{Uses}
 *
 * %% List of any external functions called by this function.
 * \begin{verbatim}
 * \end{verbatim}
 *
 * \subsubsection*{Notes}
 *
 * %% Any relevant notes.
 *
 * \vfill{\footnotesize\input{ExampleCV}}
 *
 **** </lalLaTeX> */

/**
 ** INCLUDE STANDARD LIBRARY HEADERS
 ** note LALStdLib.h already includes stdio.h and stdarg.h
 **
 ** for example:
 **
 **   #include <string.h>
 **/

/** INCLUDE ANY LAL HEADERS **/
#include <lal/LALStdlib.h>
#include <lal/Example.h>

/** DEFINE RCS ID STRING **/
NRCSID( EXAMPLEC, "$Id$" );

/**
 ** DEFINE LOCAL MACROS AND VARIABLES
 ** (these should be used sparingly)
 **
 ** macros:
 **
 **   #define STRSZ 64
 **
 **
 ** local variables must be both static and constant, for example:
 **
 **   static const INT4 strsz = STRSZ;
 **
 **
 ** enums are also allowed:
 **
 **   enum { StrSz = STRSZ };
 **
 **/

/**
 ** DECLARE LOCAL FUNCTIONS
 ** (they may be defined here or later)
 **
 ** these functions must all be static so that they do not have
 ** external linkage, for example:
 **
 **   static void
 **   LALSubroutine(
 **       LALStatus  *status,
 **       CHAR       *string,
 **       const CHAR *message,
 **       UINT4       strlen
 **       )
 **   {
 **     INITSTATUS( status, "LALSubroutine", EXAMPLEC );
 **     ASSERT( strlen > 0, status, EXAMPLEH_EIEIO, EXAMPLEH_MSGEIEIO );
 **     ASSERT( message, status, EXAMPLEH_ENULL, EXAMPLEH_MSGENULL );
 **     ASSERT( string, status, EXAMPLEH_ENULL, EXAMPLEH_MSGENULL );
 **     memcpy( string, message, strlen );
 **     string[strlen] = 0;
 **     RETURN( status );
 **   }
 **
 ** local functions need not satisfy all of the LAL convensions for
 ** functions
 **
 **/

/** DEFINE GLOBAL FUNCTIONS **/

/* <lalVerbatim file="ExampleCP"> */
void
LALExample(
    LALStatus     *status,
    ExampleOutput *output,
    ExampleInput  *input,
    ExampleParams *params
    )
/* </lalVerbatim> */
{
  /**
   ** variable declaration
   **
   ** variables should normally be declared one-per-line, for example:
   **
   **   CHAR *str;
   **   INT4  strlen;
   **
   ** no function (not even a local one) should use static variables:
   **
   **   static REAL4 x;  !!! BAD !!!
   **
   **/

  /** initialize status **/
  INITSTATUS( status, "LALExample", EXAMPLEC );
  ATTATCHSTATUSPTR( status ); /** only if subroutines are called **/

  /**
   ** check validity of arguments, for example:
   **
   **   ASSERT( output, status, EXAMPLEH_ENULL, EXAMPLEH_MSGENULL );
   **
   ** all assert statements are removed when debugging is disabled, so
   ** you should never do anything like:
   **
   **   ASSERT( strlen = strsize - 1 > 0, status, EXAMPLEH_EMOTE,
   **       EXAMPLEH_MSGEMOTE );  !!! BAD !!!
   **
   ** instead do:
   **
   **   strlen = strsize - 1;
   **   ASSERT( strlen > 0, status, EXAMPLEH_EMOTE, EXAMPLEH_MSGEMOTE );
   **
   **/

  /**
   ** bulk of code
   **
   ** here is how to call a subroutine:
   **
   **   LALSubroutine( status->statusPtr, str, "hello", strlen );
   **   CHECKSTATUSPTR( status );
   **
   ** (note: you must have attatched [sic] the status pointer)
   **
   **
   ** always check return status of memory allocation, for example:
   **
   **   str = LALMalloc( strsz );
   **   if ( ! str )
   **   {
   **     ABORT( status, EXAMPLEH_EALOC, EXAMPLEH_MSGEALOC );
   **   }
   **
   ** (note: always use ABORT rather than ASSERT to do this)
   **
   **
   ** if a subroutine might fail, but there is memory that needs to be
   ** deallocated before returning, use BEGINFAIL/ENDFAIL rather than
   ** CHECKSTATUSPTR as follows:
   **
   **   LALSubroutine( status->statusPtr, str, "world", strlen );
   **   BEGINFAIL( status )
   **   {
   **     LALFree( str );
   **   }
   **   ENDFAIL( status );
   **
   **
   ** never use I/O:
   **
   **   puts( str ); !!! BAD !!!
   **
   **
   ** make sure that all memory allocated in the function is deallocated
   ** before exiting the function (unless the function is specifically
   ** designed to create an object):
   **
   **   LALFree( str );
   **/

  DETATCHSTATUSPTR( status ); /** only if status pointer was attatched [sic] **/
  RETURN(status);
}
