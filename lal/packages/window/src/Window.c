/************************************************* <lalVerbatim file=WindowCV>
Authors: Allen, B. and Brown, D. A.
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
\item Kaiser
\end{itemize}
It should be straighforward to add additional window functions if
they are desired.
\subsubsection*{Prototypes}
\input{WindowCP}
\idx{LALWindow()}
Note that the the function \verb|LALWindow()| is depricated and will soon be
deleted. Windows should be created and destroyed by calles to the 
\verb|LALCreateREAL4Window()| and \verb|LALDestroyREAL4Window()| functions.

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
{\rm Kaiser:\ } w_j &=& I_0\left( \beta\sqrt{1 - (y - 1)^2} \right)/I_0(\beta) \\
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

static REAL8 BesselI0( REAL8 x )
{
  /*
   * 
   * Calculates the value of the 0th order, modified Bessel function of the 
   * first kind using a power series expansion. (See "Handbook of Math.
   * Functions," eds. Abramowitz and Stegun, 9.6.12 for details.) NOTE: the
   * accuracy of the expansion is chosen to be 2e-9. Stolen from Philip.
   *
   */


  REAL8 ds = 1.0;
  REAL8 d  = 0.0;
  REAL8 s  = 1.0;

  do 
  { 
    d  += 2.0; 
    ds *= x*x/(d*d); 
    s  += ds; 
  } 
  while (ds > 2e-9); 

  return s; 
} 

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
  REAL8 beta, betaI0;

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

  /* check that if the case of a Kaiser window, beta is positive */
  if ( windowtype == Kaiser )
  {
    ASSERT(parameters->beta >= 0,status, WINDOWH_EBETA,WINDOWH_MSGEBETA);
    beta = parameters->beta;
    betaI0 = BesselI0( beta );
  }

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

    case Kaiser:
      {
        REAL8 kai = (i - (length-1.0)/2.0)*2.0/((REAL8)length-1.0);
        win = BesselI0( beta * sqrt(1.0 - kai*kai) )/betaI0;
      }
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

/* <lalVerbatim file="WindowCP"> */
void LALCreateREAL4Window (
    LALStatus    *status,
    REAL4Window **output,
    LALWindowParams *params
/* </lalVerbatim> */
    )
{
  INITSTATUS( status, "LALCreateREAL4Window", WINDOW );
  ATTATCHSTATUSPTR( status );

  ASSERT( output, status, WINDOWH_ENULL, WINDOWH_MSGENULL );
  ASSERT( ! *output, status, WINDOWH_ENNUL, WINDOWH_MSGENNUL );
  ASSERT( params, status, WINDOWH_ENULL, WINDOWH_MSGENULL );

  /* allocate the storage for the window vector */
  *output = (REAL4Window *) LALCalloc( 1, sizeof(REAL4Window) );
  LALCreateVector( status->statusPtr, &((*output)->data), params->length );
  CHECKSTATUSPTR( status );

  /* compute the window */
  (*output)->type = params->type;
  LALWindow( status->statusPtr, (*output)->data, params );
  CHECKSTATUSPTR( status );

  /* copy the output params to the structure */
  (*output)->sumofsquares = params->sumofsquares;
  (*output)->beta = params->beta;
  strncpy( (*output)->windowname, params->windowname, 
      LALNameLength * sizeof(CHAR) );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="WindowCP"> */
void LALDestroyREAL4Window (
    LALStatus     *status,
    REAL4Window  **output
    )
/* </lalVerbatim> */
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
