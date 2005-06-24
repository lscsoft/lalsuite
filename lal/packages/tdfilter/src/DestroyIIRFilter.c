/***************************** <lalVerbatim file="DestroyIIRFilterCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{DestroyIIRFilter.c}}
\label{ss:DestroyIIRFilter.c}

Destroys IIR filter objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DestroyIIRFilterCP}
\idx{LALDestroyREAL4IIRFilter()}
\idx{LALDestroyREAL8IIRFilter()}

\subsubsection*{Description}

These functions destroy an object \verb@**input@ of type
\texttt{REAL4IIRFilter} or \texttt{REAL8IIRFilter}, and set
\verb@*input@ to \verb@NULL@.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
void LALFree()
void LALSDestroyVector()
void LALDDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{DestroyIIRFilterCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/IIRFilter.h>

NRCSID(DESTROYIIRFILTERC,"$Id$");

void XLALDestroyREAL4IIRFilter( REAL4IIRFilter *filter )
{
  if ( filter )
  {
    XLALDestroyREAL4Vector( filter->directCoef );
    XLALDestroyREAL4Vector( filter->recursCoef );
    XLALDestroyREAL4Vector( filter->history );
    LALFree( filter );
  }
  return;
}

void XLALDestroyREAL8IIRFilter( REAL8IIRFilter *filter )
{
  if ( filter )
  {
    XLALDestroyREAL8Vector( filter->directCoef );
    XLALDestroyREAL8Vector( filter->recursCoef );
    XLALDestroyREAL8Vector( filter->history );
    LALFree( filter );
  }
  return;
}

/* <lalVerbatim file="DestroyIIRFilterCP"> */
void
LALDestroyREAL4IIRFilter( LALStatus      *stat,
			  REAL4IIRFilter **input )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALDestroyREAL4IIRFilter",DESTROYIIRFILTERC);

  /* Make sure handle is non-null, and points to a non-null pointer.
     (The routine LALSDestroyVector will check that the data fields are
     non-null.) */
  ASSERT(input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(*input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  /* Free the filter, then point the handle to NULL. */
  XLALDestroyREAL4IIRFilter( *input );
  *input=NULL;

  /* Normal exit */
  RETURN(stat);
}


/* <lalVerbatim file="DestroyIIRFilterCP"> */
void
LALDestroyREAL8IIRFilter( LALStatus      *stat,
			  REAL8IIRFilter **input )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALDestroyREAL8IIRFilter",DESTROYIIRFILTERC);

  /* Make sure handle is non-null, and points to a non-null pointer.
     (The routine LALDDestroyVector will check that the data fields are
     non-null.) */
  ASSERT(input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(*input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  /* Free the filter, then point the handle to NULL. */
  XLALDestroyREAL8IIRFilter( *input );

  /* Normal exit */
  RETURN(stat);
}
