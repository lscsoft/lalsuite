/***************************** <lalVerbatim file="IIRFilterVectorRCV">
Author: Creighton, T. D.
$Id$
****************************** </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{IIRFilterVectorR.c}}
\label{ss:IIRFilterVectorR.c}

Applies a time-reversed IIR filter to a data stream.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{IIRFilterVectorRCP}

\subsubsection*{Description}

These functions apply a generic time-domain filter \verb@*filter@ to a
time series \verb@*vector@, as with the routines
\verb@LALIIRFilterREAL4Vector()@ and \verb@LALIIRFilterREAL8Vector()@, but
do so in a time-reversed manner.  By successively applying normal and
time-reversed IIR filters to the same data, one squares the magnitude
of the frequency response while canceling the phase shift.  This can
be significant when one wishes to preserve phase correlations across
wide frequency bands.

\subsubsection*{Algorithm}

Because these filter routines are inherently acausal, the
\verb@filter->history@ vector is meaningless and unnecessary.  These
routines neither use nor modify this data array.

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{IIRFilterVectorRCV}}

</lalLaTeX> */

#include "LALStdlib.h"
#include "IIRFilter.h"

NRCSID(IIRFILTERVECTORRC,"$Id$");


/* <lalVerbatim file="IIRFilterVectorRCP"> */
void LALIIRFilterREAL4VectorR(LALStatus         *stat,
			   REAL4Vector    *vector,
			   REAL4IIRFilter *filter)
{ /* </lalVerbatim> */
  INT4 i;            /* Loop counter for data vector. */
  INT4 j;            /* Index for filter coeficients. */
  INT4 length;       /* Length of vector. */
  REAL4 *data;       /* Vector data. */
  INT4 directOrder;  /* Number of direct filter coefficients. */
  INT4 recursOrder;  /* Number of recursive filter coefficients. */
  REAL4 *directCoef; /* Direct filter coefficients. */
  REAL4 *recursCoef; /* Recursive filter coefficients. */

  INITSTATUS(stat,"LALIIRFilterREAL4VectorR",IIRFILTERVECTORRC);

  /* Make sure all the structures have been initialized. */
  ASSERT(vector,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(vector->data,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTER_ENUL,
	 IIRFILTER_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTER_ENUL,
	 IIRFILTER_MSGENUL);
  length=vector->length;
  data=vector->data+length-1;
  directOrder=filter->directCoef->length;
  recursOrder=filter->recursCoef->length;
  directCoef=filter->directCoef->data;
  recursCoef=filter->recursCoef->data;

  /* Perform the auxilliary piece of the filter. */
  for(i=0;i<recursOrder;i++,data--)
    for(j=1;j<=i;j++)
      *data+=data[j]*recursCoef[j];
  for(;i<length;i++,data--)
    for(j=1;j<recursOrder;j++)
      *data+=data[j]*recursCoef[j];
  data++;

  /* Perform the direct piece of the filter. */
  for(;i>directOrder;i--,data++){
    *data*=directCoef[0];
    for(j=1;j<directOrder;j++)
      *data+=data[j]*directCoef[j];
  }
  for(;i>0;i--,data++){
    *data*=directCoef[0];
    for(j=1;j<i;j++)
      *data+=data[j]*directCoef[j];
  }

  /* Normal exit */
  RETURN(stat);
}


/* <lalVerbatim file="IIRFilterVectorRCP"> */
void LALIIRFilterREAL8VectorR(LALStatus         *stat,
			   REAL8Vector    *vector,
			   REAL8IIRFilter *filter)
{ /* </lalVerbatim> */
  INT4 i;            /* Loop counter for data vector. */
  INT4 j;            /* Index for filter coeficients. */
  INT4 length;       /* Length of vector. */
  REAL8 *data;       /* Vector data. */
  INT4 directOrder;  /* Number of direct filter coefficients. */
  INT4 recursOrder;  /* Number of recursive filter coefficients. */
  REAL8 *directCoef; /* Direct filter coefficients. */
  REAL8 *recursCoef; /* Recursive filter coefficients. */

  INITSTATUS(stat,"LALIIRFilterREAL8VectorR",IIRFILTERVECTORRC);

  /* Make sure all the structures have been initialized. */
  ASSERT(vector,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(vector->data,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTER_ENUL,
	 IIRFILTER_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTER_ENUL,
	 IIRFILTER_MSGENUL);
  length=vector->length;
  data=vector->data+length-1;
  directOrder=filter->directCoef->length;
  recursOrder=filter->recursCoef->length;
  directCoef=filter->directCoef->data;
  recursCoef=filter->recursCoef->data;

  /* Perform the auxilliary piece of the filter. */
  for(i=0;i<recursOrder;i++,data--)
    for(j=1;j<=i;j++)
      *data+=data[j]*recursCoef[j];
  for(;i<length;i++,data--)
    for(j=1;j<recursOrder;j++)
      *data+=data[j]*recursCoef[j];
  data++;

  /* Perform the direct piece of the filter. */
  for(;i>directOrder;i--,data++){
    *data*=directCoef[0];
    for(j=1;j<directOrder;j++)
      *data+=data[j]*directCoef[j];
  }
  for(;i>0;i--,data++){
    *data*=directCoef[0];
    for(j=1;j<i;j++)
      *data+=data[j]*directCoef[j];
  }

  /* Normal exit */
  RETURN(stat);
}
