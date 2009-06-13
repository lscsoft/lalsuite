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

/************************************ <lalVerbatim file="LALHelloCV">
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{LALHello.c}}
\label{ss:LALHello.c}

Function to print ``hello, LSC!''

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALHelloCP}
\idx{LALHello()}

\subsubsection*{Description}

The routine \verb+LALHello()+ prints the message ``hello, LSC!'' to the file
specified by \verb+fileName+, or to \verb+stdout+ if this is \verb+NULL+.
(The actual printing is handled by a subroutine \verb+LALPrintMessage()+,
which is static to the file \verb+LALHello.c+.)

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALFopen()
LALFclose()
\end{verbatim}

\subsubsection*{Notes}
This is a sample function only.

\vfill{\footnotesize\input{LALHelloCV}}

</lalLaTeX> */


#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALHello.h>

/* <lalVerbatim file="LALHelloNRCSID"> */
NRCSID( LALHELLOC, "$Id$" );
/* </lalVerbatim> */


/* don't want this function to be visible outside this file */
static void
LALPrintMessage( LALStatus *status, const CHAR *message, const CHAR *fileName )
{
  INT4  error;
  INT4  numChar;
  FILE *filePtr;

  INITSTATUS( status, "LALPrintMessage", LALHELLOC );
  ASSERT( message, status, LALHELLOH_ENMSG, LALHELLOH_MSGENMSG );

  if ( fileName )
  {
    filePtr = LALFopen( fileName, "w" );
  }
  else
  {
    filePtr = stdout;
  }

  if ( !filePtr )
  {
    ABORT( status, LALHELLOH_EOPEN, LALHELLOH_MSGEOPEN );
  }

  numChar = fprintf( filePtr, message );

  if ( numChar < 0 )
  {
    ABORT( status, LALHELLOH_EWRITE, LALHELLOH_MSGEWRITE );
  }

  if ( fileName )
  {
    error = LALFclose( filePtr );
    if ( error )
    {
      ABORT( status, LALHELLOH_ECLOSE, LALHELLOH_MSGECLOSE );
    }
  }
  else
  {
    error = fflush( stdout );
    if ( error )
    {
      ABORT( status, LALHELLOH_EFLUSH, LALHELLOH_MSGEFLUSH );
    }
  }

  RETURN( status );
}


/* <lalVerbatim file="LALHelloCP"> */
void
LALHello( LALStatus *status, const CHAR *fileName )
{ /* </lalVerbatim> */
  INITSTATUS( status, "LALHello", LALHELLOC );
  ATTATCHSTATUSPTR( status );

  LALPrintMessage( status->statusPtr, "hello, LSC!\n", fileName );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
