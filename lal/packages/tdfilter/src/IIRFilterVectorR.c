/***************************** <lalVerbatim file="IIRFilterVectorRCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{IIRFilterVectorR.c}}
\label{ss:IIRFilterVectorR.c}

Applies a time-reversed IIR filter to a data stream.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{IIRFilterVectorRCP}
\idx{LALIIRFilterREAL4VectorR()}
\idx{LALIIRFilterREAL8VectorR()}
\idx{LALDIIRFilterREAL4VectorR()}

\subsubsection*{Description}

These functions apply a generic time-domain filter \verb@*filter@ to a
time series \verb@*vector@, as with the routines
\verb@LALIIRFilterREAL4Vector()@, \verb@LALIIRFilterREAL8Vector()@,
and \verb@LALDIIRFilterREAL4Vector()@, but do so in a time-reversed
manner.  By successively applying normal and time-reversed IIR filters
to the same data, one squares the magnitude of the frequency response
while canceling the phase shift.  This can be significant when one
wishes to preserve phase correlations across wide frequency bands.

\subsubsection*{Algorithm}

Because these filter routines are inherently acausal, the
\verb@filter->history@ vector is meaningless and unnecessary.  These
routines neither use nor modify this data array.  They effectively
treat the ``future'' data as zero.

(An alternative implementation would be to assume that the filter
``history'', when invoked by these routines, stores the \emph{future}
values of the auxiliary sequence.  This would allow a large vector to
be broken into chunks and time-reverse filtered, yielding the same
result as if the whole vector had been time-reverse filtered.  I can
switch to this implementation if there is any demand for it.)

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{IIRFilterVectorRCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/IIRFilter.h>

NRCSID(IIRFILTERVECTORRC,"$Id$");


int XLALIIRFilterReverseREAL4Vector( REAL4Vector *vector, REAL8IIRFilter *filter )
{
  static const char *func = "XLALIIRFilterReverseREAL4Vector";
  INT4 j;            /* Index for filter coeficients. */
  INT4 length;       /* Length of vector. */
  REAL4 *data;       /* Vector data. */
  REAL8 w, datum;    /* Current auxiliary and output values. */
  INT4 directOrder;  /* Number of direct filter coefficients. */
  INT4 recursOrder;  /* Number of recursive filter coefficients. */
  INT4 numHist;      /* The number of auxiliary data kept. */
  REAL8 *directCoef; /* Direct filter coefficients. */
  REAL8 *recursCoef; /* Recursive filter coefficients. */
  REAL8 *temp=NULL;  /* Temporary storage for auxiliary sequence. */

  /* Make sure all the structures have been initialized. */
  if ( ! vector || ! filter )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! vector->data )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! filter->directCoef || ! filter->recursCoef || ! filter->history
      || !  filter->directCoef->data || ! filter->recursCoef->data
      || !  filter->history->data )
    XLAL_ERROR( func, XLAL_EINVAL );

  length=vector->length;
  data=vector->data+length-1;
  directOrder=filter->directCoef->length;
  recursOrder=filter->recursCoef->length;
  directCoef=filter->directCoef->data;
  recursCoef=filter->recursCoef->data;
  numHist=filter->history->length+1;
  temp = LALMalloc( numHist*sizeof(*temp) );
  if ( ! temp )
    XLAL_ERROR( func, XLAL_ENOMEM );
  memset(temp,0,(numHist-1)*sizeof(*temp));

  /* Run through the vector. */
  while(length--){

    /* Compute the auxiliary variable. */
    for(j=numHist-1;j>=recursOrder;j--)
      temp[j]=temp[j-1];
    w=*data;
    for(;j;j--)
      w+=recursCoef[j]*(temp[j]=temp[j-1]);

    /* Compute filter output. */
    datum=*directCoef*(*temp=w);
    for(j=1;j<directOrder;j++)
      datum+=directCoef[j]*temp[j];
    *(data--)=datum;
  }
  LALFree(temp);

  /* Normal exit */
  return 0;
}


int XLALIIRFilterReverseREAL8Vector( REAL8Vector *vector, REAL8IIRFilter *filter )
{
  static const char *func = "XLALIIRFilterReverseREAL8Vector";
  INT4 i;            /* Loop counter for data vector. */
  INT4 j;            /* Index for filter coeficients. */
  INT4 length;       /* Length of vector. */
  REAL8 *data;       /* Vector data. */
  REAL8 datum;       /* Temporary working variable. */
  INT4 directOrder;  /* Number of direct filter coefficients. */
  INT4 recursOrder;  /* Number of recursive filter coefficients. */
  REAL8 *directCoef; /* Direct filter coefficients. */
  REAL8 *recursCoef; /* Recursive filter coefficients. */

  /* Make sure all the structures have been initialized. */
  if ( ! vector || ! filter )
    XLAL_ERROR( func, XLAL_EFAULT );
  if ( ! vector->data )
    XLAL_ERROR( func, XLAL_EINVAL );
  if ( ! filter->directCoef || ! filter->recursCoef || ! filter->history
      || !  filter->directCoef->data || ! filter->recursCoef->data
      || !  filter->history->data )
    XLAL_ERROR( func, XLAL_EINVAL );

  length=vector->length;
  data=vector->data+length-1;
  directOrder=filter->directCoef->length;
  recursOrder=filter->recursCoef->length;
  directCoef=filter->directCoef->data;
  recursCoef=filter->recursCoef->data;

  /* Perform the auxilliary piece of the filter. */
  for(i=0;(i<recursOrder)&&(i<length);i++,data--){
    datum=*data;
    for(j=1;j<=i;j++)
      datum+=data[j]*recursCoef[j];
    *data=datum;
  }
  for(;i<length;i++,data--){
    datum=*data;
    for(j=1;j<recursOrder;j++)
      datum+=data[j]*recursCoef[j];
    *data=datum;
  }
  data++;

  /* Perform the direct piece of the filter. */
  for(;i>directOrder;i--,data++){
    datum=*data*directCoef[0];
    for(j=1;j<directOrder;j++)
      datum+=data[j]*directCoef[j];
    *data=datum;
  }
  for(;i>0;i--,data++){
    datum=*data*directCoef[0];
    for(j=1;j<i;j++)
      datum+=data[j]*directCoef[j];
    *data=datum;
  }

  /* Normal exit */
  return 0;
}



/*
 *
 * WARNING: THIS FUNCTION IS OBSOLETE!
 *
 */

/* <lalVerbatim file="IIRFilterVectorRCP"> */
void
LALIIRFilterREAL4VectorR( LALStatus      *stat,
			  REAL4Vector    *vector,
			  REAL4IIRFilter *filter )
{ /* </lalVerbatim> */
  INT4 i;            /* Loop counter for data vector. */
  INT4 j;            /* Index for filter coeficients. */
  INT4 length;       /* Length of vector. */
  REAL4 *data;       /* Vector data. */
  REAL8 datum;       /* Temporary working variable. */
  INT4 directOrder;  /* Number of direct filter coefficients. */
  INT4 recursOrder;  /* Number of recursive filter coefficients. */
  REAL4 *directCoef; /* Direct filter coefficients. */
  REAL4 *recursCoef; /* Recursive filter coefficients. */

  INITSTATUS(stat,"LALIIRFilterREAL4VectorR",IIRFILTERVECTORRC);

  /* Make sure all the structures have been initialized. */
  ASSERT(vector,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(vector->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTERH_ENUL,
      IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTERH_ENUL,
      IIRFILTERH_MSGENUL);
  length=vector->length;
  data=vector->data+length-1;
  directOrder=filter->directCoef->length;
  recursOrder=filter->recursCoef->length;
  directCoef=filter->directCoef->data;
  recursCoef=filter->recursCoef->data;

  /* Perform the auxilliary piece of the filter. */
  for(i=0;(i<recursOrder)&&(i<length);i++,data--){
    datum=*data;
    for(j=1;j<=i;j++)
      datum+=data[j]*recursCoef[j];
    *data=datum;
  }
  for(;i<length;i++,data--){
    datum=*data;
    for(j=1;j<recursOrder;j++)
      datum+=data[j]*recursCoef[j];
    *data=datum;
  }
  data++;

  /* Perform the direct piece of the filter. */
  for(;i>directOrder;i--,data++){
    datum=*data*directCoef[0];
    for(j=1;j<directOrder;j++)
      datum+=data[j]*directCoef[j];
    *data=datum;
  }
  for(;i>0;i--,data++){
    datum=*data*directCoef[0];
    for(j=1;j<i;j++)
      datum+=data[j]*directCoef[j];
    *data=datum;
  }

  /* Normal exit */
  RETURN(stat);
}


/* <lalVerbatim file="IIRFilterVectorRCP"> */
void
LALIIRFilterREAL8VectorR( LALStatus      *stat,
			  REAL8Vector    *vector,
			  REAL8IIRFilter *filter )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALIIRFilterREAL8VectorR",IIRFILTERVECTORRC);

  /* Make sure all the structures have been initialized. */
  ASSERT(vector,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(vector->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);

  if ( XLALIIRFilterReverseREAL8Vector( vector, filter ) < 0 )
  {
    XLALClearErrno();
    ABORT(stat,IIRFILTERH_EMEM,IIRFILTERH_MSGEMEM);
  }

  /* Normal exit */
  RETURN(stat);
}


/* <lalVerbatim file="IIRFilterVectorRCP"> */
void
LALDIIRFilterREAL4VectorR( LALStatus      *stat,
			   REAL4Vector    *vector,
			   REAL8IIRFilter *filter )
{ /* </lalVerbatim> */
  INITSTATUS(stat,"LALDIIRFilterREAL4VectorR",IIRFILTERVECTORRC);

  /* Make sure all the structures have been initialized. */
  ASSERT(vector,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(vector->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTERH_ENUL,
	 IIRFILTERH_MSGENUL);

  if ( XLALIIRFilterReverseREAL4Vector( vector, filter ) < 0 )
  {
    XLALClearErrno();
    ABORT(stat,IIRFILTERH_EMEM,IIRFILTERH_MSGEMEM);
  }

  /* Normal exit */
  RETURN(stat);
}
