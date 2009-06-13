/*
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix
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

/************************************ <lalVerbatim file="UserInputTestCV">
Author: Reinhard Prix
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Program \texttt{UserInputTest.c}}
\label{s:UserInputTest.c}

Tests the routines in \verb@UserInput.h@.

\subsubsection*{Usage}
\begin{verbatim}
UserInputTest
\end{verbatim}

\subsubsection*{Description}

Do some standard-tests for the config-file reading routines. No
extensive error-condition checking is done here, we only check if the
basic functionality works.

\subsubsection*{Exit codes}

\input{UserInputErrors}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{UserInputTestCV}}

</lalLaTeX> */

#include <lal/UserInput.h>

NRCSID (USERINPUTTESTC, "$Id$");

/* Error codes and messages */

/************** <lalErrTable file="UserInputErrors"> */
#define USERINPUTTESTC_ENORM 		0

#define USERINPUTTESTC_MSGENORM 	"Normal exit"

/******************************************** </lalErrTable> */


/* Default parameters. */

INT4 lalDebugLevel=3;


/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    LALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                   "        %s %s\n", (code), *argv, __FILE__,       \
              __LINE__, USERINPUTTESTC, statement ? statement :  \
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    LALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              USERINPUTTESTC, (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( USERINPUTTESTC_ESUB, USERINPUTTESTC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return USERINPUTTESTC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

int main(int argc, char *argv[])
{

  /* EMPTY STUB */
/*
  printf ("User said: %d arguments, which are %s %s %s ...\n", argc, argv[0], argv[1], argv[2]);
*/
  argc = 0; argv = NULL;


  /* sorry, something still should be written here */
  LALCheckMemoryLeaks();

  return (USERINPUTTESTC_ENORM);
}


