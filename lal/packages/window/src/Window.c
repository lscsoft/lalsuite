/****************************************************** <lalVerbatim file=WindowCV>
Author: Bruce Allen
$Id$
****************************************************** </lalVerbatim>*/
/****************************************************** <lalLaTeX>

\subsection{Module \texttt{Window.c}}
\label{ss:Window.c}

Creates vector structure containing a window (also called
a taper, lag window, or apodization function).  The choices
currently available are:
\begin{itemize}
\item Rectangular
\item Hann
\item Welch
\item Bartlett
\item Parzen
\item Papoulis
\item Hamming
\end{itemize}
It should be straighforward to add additional window functions if
they are desired.
\subsubsection*{Prototypes}
\input{WindowCP}
\idx{LALWindow()}
Note that the \texttt{paramters} argument handles both input and output.

\subsubsection*{Description}
This function creates a time-domain window function in a vector of
specified length.  Note that this function was not written to be
particularly efficient.  If you need a window lots of times, calculate
it once then save it, please.

The window functions are defined for $j=0,\cdots,N-1$ by the following
formulae.  Note that $N$ is the vector.  In these formulae, let
$x=2 \pi j/N$, and $y=|2j/N-1|$,
\begin{eqnarray*}
{\rm Rectangular:\ } w_j &=& 1 \\
{\rm Hann:\ } w_j &=& {1 \over 2} ( 1 - \cos  x  ) \\
{\rm Welch:\ } w_j &=& 1 -  y^2 \\
{\rm Bartlett:\ } w_j &=& 1 -  y \\
{\rm Parzen:\ } w_j &=&  1 - 6 y^2 + 6 y^3  {\rm\ if\ } y\le 1/2\\
                    &=&  2 (1-y)^3 {\rm\ if\ } y>1/2\\
{\rm Papoulis:\ } w_j &=& {1 \over \pi} \sin (\pi  y  ) + ( 1 -  y  ) \cos (\pi  y  )\\
{\rm Hamming:\ } w_j &=& 1-0.46 (1 + \cos x ) \\
\end{eqnarray*}
These window functions are shown in Fig.~\ref{f:window} for $N=1024$.
\begin{figure}
\noindent\includegraphics[angle=-90,width=.9\linewidth]{windowFig}
\caption{\label{f:window} Examples of the window functions for length 1024}
\end{figure}

****************************************************** </lalLaTeX> */

#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/Window.h>

NRCSID (WINDOW, "$Id$");

static const char *WindowTypeNames[] = WINDOWNAMELIST;

/* <lalVerbatim file="WindowCP"> */
void 
LALWindow( LALStatus       *status, 
           REAL4Vector     *vector, 
           LALWindowParams *parameters ) 
     /* </lalVerbatim> */
{
  UINT4 i;
  UINT4 length;
  INT4 windowtype;
  REAL8 wss;    /* window summed and squared */
  REAL8 win;
  REAL8 x,y,z;

  /* Initialize status structure   */
  INITSTATUS(status,"LALWindow Function",WINDOW);

  /* Check that parameter block is there. */ 
  ASSERT(parameters!=NULL,status,WINDOWH_ENULLPARAM,WINDOWH_MSGENULLPARAM);

  /* check that the vector is not null */
  ASSERT(vector!=NULL,status,WINDOWH_ENULLHANDLE,WINDOWH_MSGENULLHANDLE);

  /* Check that window length is reasonable. */ 
  length=parameters->length;
  ASSERT(length>0,status,WINDOWH_EELENGTH,WINDOWH_MSGEELENGTH);

  /* Make sure that window is of a known type */
  windowtype=parameters->type;
  ASSERT(windowtype>=Rectangular && windowtype<NumberWindowTypes,status,
         WINDOWH_ETYPEUNKNOWN,WINDOWH_MSGETYPEUNKNOWN);

  /* vector is apparently already allocated.  Check length, data area */
  ASSERT(vector->length==length,status,
         WINDOWH_EWRONGLENGTH,WINDOWH_MSGEWRONGLENGTH);
  ASSERT(vector->data!=NULL,status,WINDOWH_ENULLDATA,WINDOWH_MSGENULLDATA);

  wss=0.0;
  for (i=0;i<length;i++)
  {
    x=(2.0*LAL_PI*i)/length;
    y=fabs(2.0*i/length-1.0);

    switch (windowtype)
    {
    /* rectangular (no) window */
    case Rectangular:
      win=1.0;
      break;

    /* Hann window */
    case Hann:
      win=0.5*(1.0-cos(x));
      break;

    /* Welch window */
    case Welch:
      win=1.0-y*y;
      break;

    /* Bartlett window */
    case Bartlett:
      win=1.0-y;
      break;

    /* Parzen window */
    case Parzen:
      z=1.0-y;
      if (y<=0.5)
        win=1.0-6.0*y*y*z;
      else
        win=2.0*z*z*z;
      break;

    /* Papoulis window */
    case Papoulis:
      win=1.0/LAL_PI*sin(LAL_PI*y)+(1.0-y)*cos(LAL_PI*y);
      break;

    case Hamming:
      win=1.0-0.46*(1.0+cos(x));
      break;

    /* Default case -- this will NEVER happen -- it is trapped above! */
    default:
      ABORT(status,WINDOWH_ETYPEUNKNOWN,WINDOWH_MSGETYPEUNKNOWN);
      break;
    }
    wss+=win*win;
    vector->data[i]=(REAL4)win;
  }
  parameters->sumofsquares=wss;
  parameters->windowname=WindowTypeNames[parameters->type];

  RETURN(status);
}

void LALCreateREAL4Window (
    LALStatus    *status,
    REAL4Window **output,
    UINT4         length,
    WindowType    type
    )
{
  LALWindowParams       wpars;

  INITSTATUS( status, "LALCreateREAL4Window", WINDOW );
  ATTATCHSTATUSPTR( status );

  ASSERT( output, status, WINDOWH_ENULL, WINDOWH_MSGENULL );
  ASSERT( ! *output, status, WINDOWH_ENNUL, WINDOWH_MSGENNUL );

  /* allocate the storage for the window vector */
  *output = (REAL4Window *) LALCalloc( 1, sizeof(REAL4Window) );
  LALCreateVector( status->statusPtr, &((*output)->data), length );
  CHECKSTATUSPTR( status );

  /* compute the window */
  (*output)->type = wpars.type = type;
  wpars.length = (INT4) length;
  LALWindow( status->statusPtr, (*output)->data, &wpars );
  CHECKSTATUSPTR( status );

  /* copy the output params to the structure */
  (*output)->sumofsquares = wpars.sumofsquares;
  strncpy( (*output)->windowname, wpars.windowname, 
      LALNameLength * sizeof(CHAR) );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

void LALDestroyREAL4Window (
    LALStatus     *status,
    REAL4Window  **output
    )
{
  INITSTATUS( status, "LALCreateREAL4Window", WINDOW );
  ATTATCHSTATUSPTR( status );

  ASSERT( output, status, WINDOWH_ENULL, WINDOWH_MSGENULL );
  ASSERT( *output, status, WINDOWH_ENULL, WINDOWH_MSGENULL );

  /* destroy the window */
  LALDestroyVector( status->statusPtr, (*output)->data );
  CHECKSTATUSPTR( status );
  LALFree( *output );
  *output = NULL;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
