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
\index{\verb&LALHello()&}

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
#include "LALStdlib.h"
#include "LALHello.h"

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
