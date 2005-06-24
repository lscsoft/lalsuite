/************************************ <lalVerbatim file="IIRFilterCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{IIRFilter.c}}
\label{ss:IIRFilter.c}

Computes an instant-by-instant IIR filter response.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{IIRFilterCP}
\idx{LALIIRFilterREAL4()}
\idx{LALIIRFilterREAL8()}
\idx{LALSIIRFilter()}
\idx{LALDIIRFilter()}

\subsubsection*{Description}

These functions pass a time-domain datum to an object \verb@*filter@
of type \verb@REAL4IIRFilter@ or \verb@REAL8IIRFilter@, and return the
filter response.  This is done using the auxiliary data series
formalism described in the header \verb@IIRFilter.h@.

There are two pairs of routines in this module.  The functions
\verb@IIRFilterReal4()@ and \verb@LALIIRFilterREAL8()@ conform to the LAL
standard, with status handling and error trapping; the input datum is
passed in as \verb@input@ and the response is returned in
\verb@*output@.  The functions \verb@LALSIIRFilter()@ and
\verb@LALDIIRFilter()@ are non-standard lightweight routines, which may
be more suitable for multiple callings from the inner loops of
programs; they have no status handling or error trapping.  The input
datum is passed in by the variable \verb@x@, and the response is
returned through the function's return statement.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{IIRFilterCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/IIRFilter.h>

NRCSID(IIRFILTERC,"$Id$");


/*
 *
 * WARNING: THIS FUNCTION IS OBSOLETE
 *
 */
/* <lalVerbatim file="IIRFilterCP"> */
void
LALIIRFilterREAL4( LALStatus      *stat,
		   REAL4          *output,
		   REAL4          input,
		   REAL4IIRFilter *filter )
{ /* </lalVerbatim> */
  INT4 j;      /* Index for filter coefficients. */
  INT4 jmax;   /* Number of filter coefficients. */
  REAL4 *coef; /* Values of filter coefficients. */
  REAL4 *hist; /* Values of filter history. */

  INITSTATUS(stat,"LALIIRFilterREAL4",IIRFILTERC);

  /* Check all the passed parameters for null pointers. */
  ASSERT(output,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->history,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->history->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  /* Compute the auxiliary datum. */
  jmax=filter->recursCoef->length;
  coef=filter->recursCoef->data+1;
  hist=filter->history->data;
  for(j=1;j<jmax;j++)
    input+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Compute the filter output. */
  jmax=filter->directCoef->length;
  coef=filter->directCoef->data;
  *output=(*(coef++))*input;
  for(j=1;j<jmax;j++)
    *output+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Updata the filter history. */
  jmax=filter->history->length-1;
  hist+=jmax;
  for(j=jmax;j>0;j--,hist--)
    *hist=hist[-1];
  *hist=input;

  /* Normal exit */
  RETURN(stat);
}


/* <lalVerbatim file="IIRFilterCP"> */
void
LALIIRFilterREAL8( LALStatus      *stat,
		   REAL8          *output,
		   REAL8          input,
		   REAL8IIRFilter *filter )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALIIRFilterREAL8",IIRFILTERC);

  /* Check all the passed parameters for null pointers. */
  ASSERT(output,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->history,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->history->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);

  *output=XLALIIRFilterREAL8(input,filter);

  /* Normal exit */
  RETURN(stat);
}


/*
 *
 * WARNING: THIS FUNCTION IS OBSOLETE
 *
 */
/* <lalVerbatim file="IIRFilterCP"> */
REAL4
LALSIIRFilter( REAL4 x, REAL4IIRFilter *filter )
{ /* </lalVerbatim> */
  INT4 j;      /* Index for filter coefficients. */
  INT4 jmax;   /* Number of filter coefficients. */
  REAL4 *coef; /* Values of filter coefficients. */
  REAL4 *hist; /* Values of filter history. */
  REAL4 w;     /* Auxiliary datum. */
  REAL4 y;     /* Output datum. */

  /* Compute the auxiliary datum. */
  jmax=filter->recursCoef->length;
  coef=filter->recursCoef->data+1;
  hist=filter->history->data;
  w=x;
  for(j=1;j<jmax;j++)
    w+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Compute the filter output. */
  jmax=filter->directCoef->length;
  coef=filter->directCoef->data;
  y=(*(coef++))*w;
  for(j=1;j<jmax;j++)
    y+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Updata the filter history. */
  jmax=filter->history->length-1;
  hist+=jmax;
  for(j=jmax;j>0;j--,hist--)
    *hist=hist[-1];
  *hist=w;

  /* Normal exit */
  return y;
}


/* <lalVerbatim file="IIRFilterCP"> */
REAL8 XLALIIRFilterREAL8( REAL8 x, REAL8IIRFilter *filter )
{ /* </lalVerbatim> */
  INT4 j;      /* Index for filter coefficients. */
  INT4 jmax;   /* Number of filter coefficients. */
  REAL8 *coef; /* Values of filter coefficients. */
  REAL8 *hist; /* Values of filter history. */
  REAL8 w;     /* Auxiliary datum. */
  REAL8 y;     /* Output datum. */

  /* Compute the auxiliary datum. */
  jmax=filter->recursCoef->length;
  coef=filter->recursCoef->data+1;
  hist=filter->history->data;
  w=x;
  for(j=1;j<jmax;j++)
    w+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Compute the filter output. */
  jmax=filter->directCoef->length;
  coef=filter->directCoef->data;
  y=(*(coef++))*w;
  for(j=1;j<jmax;j++)
    y+=(*(coef++))*(*(hist++));
  hist-=(jmax-1);

  /* Updata the filter history. */
  jmax=filter->history->length-1;
  hist+=jmax;
  for(j=jmax;j>0;j--,hist--)
    *hist=hist[-1];
  *hist=w;

  /* Normal exit */
  return y;
}

/* <lalVerbatim file="IIRFilterCP"> */
REAL4 XLALIIRFilterREAL4( REAL4 x, REAL8IIRFilter *filter )
{ /* </lalVerbatim> */
  return (REAL4)XLALIIRFilterREAL8((REAL8)x,filter);
}
