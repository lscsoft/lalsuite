/****************************** <lalVerbatim file="IIRFilterVectorCV">
Author: Creighton, T. D.
$Id$
******************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{IIRFilterVector.c}}
\label{ss:IIRFilterVector.c}

Applies an IIR filter to a data stream.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{IIRFilterVectorCP}
\index{\texttt{LALIIRFilterREAL4Vector()}}
\index{\texttt{LALIIRFilterREAL8Vector()}}

\subsubsection*{Description}

These functions apply a generic time-domain filter given by an object
\verb@*filter@ of type \verb@REAL4IIRFilter@ or \verb@REAL8IIRFilter@
to a list \verb@*vector@ of data representing a time series.  This is
done in place using the auxiliary data series formalism described in
\verb@IIRFilter.h@, so as to accomodate potentially large data series.
To filter a piece of a larger dataset, the calling routine may pass a
vector structure whose data pointer and length fields specify a subset
of a larger vector.

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()
LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{IIRFilterVectorCV}}

</lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/IIRFilter.h>

NRCSID(IIRFILTERVECTORC,"$Id$");


/* <lalVerbatim file="IIRFilterVectorCP"> */
void
LALIIRFilterREAL4Vector( LALStatus      *stat,
			 REAL4Vector    *vector,
			 REAL4IIRFilter *filter )
{ /* </lalVerbatim> */
  INT4 i;            /* Loop counter for data vector. */
  INT4 j;            /* Index for filter coeficients. */
  INT4 k;            /* Index for filter history. */
  INT4 length;       /* Length of vector. */
  REAL4 *data;       /* Vector data. */
  INT4 directOrder;  /* Number of direct filter coefficients. */
  INT4 recursOrder;  /* Number of recursive filter coefficients. */
  INT4 numHist;      /* The number of history data. */
  REAL4 *directCoef; /* Direct filter coefficients. */
  REAL4 *recursCoef; /* Recursive filter coefficients. */
  REAL4 *history;    /* Filter history. */
  REAL4 *temp=NULL;  /* Temporary storage for the filter history. */

  INITSTATUS(stat,"LALIIRFilterREAL4Vector",IIRFILTERVECTORC);

  /* Make sure all the structures have been initialized. */
  ASSERT(vector,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(vector->data,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->history,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTER_ENUL,
	 IIRFILTER_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTER_ENUL,
	 IIRFILTER_MSGENUL);
  ASSERT(filter->history->data,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  length=vector->length;
  data=vector->data;
  directOrder=filter->directCoef->length;
  recursOrder=filter->recursCoef->length;
  numHist=filter->history->length;
  directCoef=filter->directCoef->data;
  recursCoef=filter->recursCoef->data;
  history=filter->history->data;
  temp=(REAL4 *)LALMalloc(numHist*sizeof(REAL4));
  if ( !temp ) {
    ABORT(stat,IIRFILTER_EMEM,IIRFILTER_MSGEMEM);
  }

  /* Compute the auxiliary data series. */
  for(i=0;(i<recursOrder)&&(i<length);i++,data++){
    for(j=1;j<=i;j++)
      *data+=data[-j]*recursCoef[j];
    for(k=0;j<recursOrder;j++,k++)
      *data+=history[k]*recursCoef[j];
  }
  for(;i<length;i++,data++){
    for(j=1;j<recursOrder;j++)
      *data+=data[-j]*recursCoef[j];
  }
  data--;

  /* Store the last few auxiliary data to the temporary history. */
  for(k=numHist-1;k>=length;k--)
    temp[k]=history[k-length];
  for(;k>=0;k--)
    temp[k]=data[-k];

  /* Compute the output data series. */
  for(;i>directOrder;i--,data--){
    *data*=directCoef[0];
    for(j=1;j<directOrder;j++)
      *data+=data[-j]*directCoef[j];
  }
  for(;i>0;i--,data--){
    *data*=directCoef[0];
    for(j=1;j<i;j++)
      *data+=data[-j]*directCoef[j];
    for(k=0;j<directOrder;j++,k++)
      *data+=history[k]*directCoef[j];
  }

  /* Update the filter history from the temporary history. */
  for(k=0;k<numHist;k++)
    history[k]=temp[k];
  LALFree(temp);

  /* Normal exit */
  RETURN(stat);
}


/* <lalVerbatim file="IIRFilterVectorCP"> */
void
LALIIRFilterREAL8Vector( LALStatus      *stat,
			 REAL8Vector    *vector,
			 REAL8IIRFilter *filter )
{ /* </lalVerbatim> */
  INT4 i;            /* Loop counter for data vector. */
  INT4 j;            /* Index for filter coeficients. */
  INT4 k;            /* Index for filter history. */
  INT4 length;       /* Length of vector. */
  REAL8 *data;       /* Vector data. */
  INT4 directOrder;  /* Number of direct filter coefficients. */
  INT4 recursOrder;  /* Number of recursive filter coefficients. */
  INT4 numHist;      /* The number of history data. */
  REAL8 *directCoef; /* Direct filter coefficients. */
  REAL8 *recursCoef; /* Recursive filter coefficients. */
  REAL8 *history;    /* Filter history. */
  REAL8 *temp=NULL;  /* Temporary storage for the filter history. */

  INITSTATUS(stat,"LALIIRFilterREAL8Vector",IIRFILTERVECTORC);

  /* Make sure all the structures have been initialized. */
  ASSERT(vector,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(vector->data,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->directCoef,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->recursCoef,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->history,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(filter->directCoef->data,stat,IIRFILTER_ENUL,
	 IIRFILTER_MSGENUL);
  ASSERT(filter->recursCoef->data,stat,IIRFILTER_ENUL,
	 IIRFILTER_MSGENUL);
  ASSERT(filter->history->data,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  length=vector->length;
  data=vector->data;
  directOrder=filter->directCoef->length;
  recursOrder=filter->recursCoef->length;
  numHist=filter->history->length;
  directCoef=filter->directCoef->data;
  recursCoef=filter->recursCoef->data;
  history=filter->history->data;
  temp=(REAL8 *)LALMalloc(numHist*sizeof(REAL8));
  if ( !temp ) {
    ABORT(stat,IIRFILTER_EMEM,IIRFILTER_MSGEMEM);
  }

  /* Compute the auxiliary data series. */
  for(i=0;(i<recursOrder)&&(i<length);i++,data++){
    for(j=1;j<=i;j++)
      *data+=data[-j]*recursCoef[j];
    for(k=0;j<recursOrder;j++,k++)
      *data+=history[k]*recursCoef[j];
  }
  for(;i<length;i++,data++){
    for(j=1;j<recursOrder;j++)
      *data+=data[-j]*recursCoef[j];
  }
  data--;

  /* Store the last few auxiliary data to the temporary history. */
  for(k=numHist-1;k>=length;k--)
    temp[k]=history[k-length];
  for(;k>=0;k--)
    temp[k]=data[-k];

  /* Compute the output data series. */
  for(;i>directOrder;i--,data--){
    *data*=directCoef[0];
    for(j=1;j<directOrder;j++)
      *data+=data[-j]*directCoef[j];
  }
  for(;i>0;i--,data--){
    *data*=directCoef[0];
    for(j=1;j<i;j++)
      *data+=data[-j]*directCoef[j];
    for(k=0;j<directOrder;j++,k++)
      *data+=history[k]*directCoef[j];
  }

  /* Update the filter history from the temporary history. */
  for(k=0;k<numHist;k++)
    history[k]=temp[k];
  LALFree(temp);

  /* Normal exit */
  RETURN(stat);
}
