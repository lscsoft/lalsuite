/************************************ <lalVerbatim file="LALWienerCV">
Author: Dwyer, Steven J.
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{LALWiener.c}}
\label{ss:LALWiener.c}

This package passes input data through a wiener filter

\subsubsection*{Prototypes}
\input{LALWienerCP}
\index{\texttt{LALUnFormattedWiener()}}
\index{\texttt{LALFormattedWiener()}}

\subsubsection*{Description}
This package contains two functions that utilize wiener filtering to remove noise.  The first function accepts as input a template REAL4Vector and a source REAL4Vector containing the unfiltered data.  The second function accepts a COMPLEX8 template that was fourier transformed, a source REAL4Vector containing the unfiltered data, and an INT2 of the length of the untransformed template plus the length of the source data. 


\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALWienerCV}}

*********************************************************** </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/LALWiener.h>


NRCSID( LALWIENERC, "$Id$" );

/* <lalVerbatim file="LALWienerCP"> */
void LALUnFormattedWiener
( 
     LALStatus *stat, 
     WienerOutput *output, 
     WienerUnFormattedInput *input 
)
{ /* </lalVerbatim> */

  /*  variable delcarations  */
  INT2			length;		/*  length of complex vector passed to LALFormattedWiener */ 
  INT2			count1, count2/*, count3*/; /* iterators */
  REAL4Vector		*sz ;		/* zero padded version of WienerUnformattedInput.s */
  WienerFormattedInput	hstl;		/* finalproduct sentto LALFormattedWiener */
  RealFFTPlan		*plan;		/* plan used to create stil */

  INITSTATUS( stat, "LALUnFormattedWiener", LALWIENERC );

  ATTATCHSTATUSPTR( stat );
  CHECKSTATUSPTR(stat);

  /*  Input must be defined  */
  ASSERT ( input, stat, LALWIENERH_ENULL, LALWIENERH_MSGENULL);
  CHECKSTATUSPTR(stat);

  /*  Output must be defined  */
  ASSERT ( output, stat, LALWIENERH_ENULL, LALWIENERH_MSGENULL);
  CHECKSTATUSPTR(stat);

  /*  data must be defined  */
  ASSERT( (input->h->data),stat,LALWIENERH_EDATA,LALWIENERH_MSGEDATA);
  CHECKSTATUSPTR(stat);

  /*  data must be defined  */
  ASSERT( (input->s->data),stat,LALWIENERH_EDATA,LALWIENERH_MSGEDATA);
  CHECKSTATUSPTR(stat);

  /*  length of h must be >= length of s  */
  ASSERT ( input->h->length >= input->s->length, stat, LALWIENERH_ESIZE, LALWIENERH_MSGESIZE);
  CHECKSTATUSPTR(stat);

  /*  Variable Initialize  */
  plan	  = NULL;
  sz	  = NULL;
  hstl.h  = NULL;
  hstl.st = NULL;

  length = (input->s->length + input->h->length);

  /*  Allocate space for data  */
  LALSCreateVector(stat->statusPtr, &sz, length);
  CHECKSTATUSPTR(stat);

  LALSCreateVector(stat->statusPtr, &hstl.h, input->h->length);
  CHECKSTATUSPTR(stat);

  LALCCreateVector(stat->statusPtr, &hstl.st,(length/2)+1);
  CHECKSTATUSPTR(stat);

  /* zero padding */
  for (count1 = 0; count1 <(INT4)(input->h->length)   ; count1++)
      hstl.h->data[count1] = input->h->data[count1];

  for (count1 = 0; count1 <(INT4)(input->s->length)   ; count1++)
      sz->data[count1] = input->s->data[count1];

  for (count2 = count1; count2 < (length); count2++)
      sz->data[count2] = 0;
 
  LALCreateForwardRealFFTPlan(stat->statusPtr, &plan, (length), 0);
  CHECKSTATUSPTR(stat);

  LALForwardRealFFT(stat->statusPtr, hstl.st, sz, plan); 
  CHECKSTATUSPTR(stat);

  hstl.TotalLength = &length;

  LALFormattedWiener(stat->statusPtr, output, &hstl);
  CHECKSTATUSPTR(stat);

  /* freeing up memory */
  LALSDestroyVector(stat->statusPtr, &sz);
  CHECKSTATUSPTR(stat);

  LALSDestroyVector(stat->statusPtr, &hstl.h);
  CHECKSTATUSPTR(stat);

  LALCDestroyVector(stat->statusPtr, &hstl.st);
  CHECKSTATUSPTR(stat);

  LALDestroyRealFFTPlan(stat->statusPtr, &plan);
  CHECKSTATUSPTR(stat);

  DETATCHSTATUSPTR( stat );  
}

/* <lalVerbatim file="LALWienerCP"> */
void LALFormattedWiener
( 
     LALStatus *stat,
     WienerOutput *output,
     WienerFormattedInput *input
)
{ /* </lalVerbatim> */

  /*  Variable Declarations  */
  /*INT2 delta;*/			/* difference in length between hstl.st & hstl.h */ 
  INT2 count1, count2;	        /* iterators */
  REAL4 N;			/* this is a real4 version of input->TotalLength */
  REAL4Vector *hz;		/* zero padded version of hstl.h */
  REAL4Vector *largeq;		/* a superset of the output q */
  COMPLEX8Vector *htil;		/* transformed version of h */
  COMPLEX8Vector *qprime;	/* q before it's transformed into end result*/
  COMPLEX8Vector *conjs;	/* complex conjucate of hstl.st */
  RealFFTPlan *plan;		/* plan used to create htil */
  RealFFTPlan *plan2;		/* plan used to create largeq */
  
  INITSTATUS( stat, "LALFormattedWiener", LALWIENERC );
  ATTATCHSTATUSPTR( stat );
  N = *(input->TotalLength);
  
  /*  Input must be defined  */
  ASSERT ( input, stat, LALWIENERH_ENULL, LALWIENERH_MSGENULL);
  CHECKSTATUSPTR(stat);

  /*  Output must be defined  */
  ASSERT ( output, stat, LALWIENERH_ENULL, LALWIENERH_MSGENULL);
  CHECKSTATUSPTR(stat);

  /*  data must be defined  */
  ASSERT( (input->st->data),stat,LALWIENERH_EDATA,LALWIENERH_MSGEDATA);
  CHECKSTATUSPTR(stat);

  /*  data must be defined  */
  ASSERT( (input->h->data),stat,LALWIENERH_EDATA,LALWIENERH_MSGEDATA);
  CHECKSTATUSPTR(stat);

  /*  Initialize Variables  */
  plan	 = NULL;
  plan2  = NULL;
  hz	 = NULL;
  largeq = NULL;
  htil	 = NULL;
  qprime = NULL;
  conjs  = NULL;

  /*  Allocate space for data  */
  LALSCreateVector(stat->statusPtr, &hz, *(input->TotalLength));
  CHECKSTATUSPTR(stat);
  LALSCreateVector(stat->statusPtr, &largeq, *(input->TotalLength));
  CHECKSTATUSPTR(stat);
  LALCCreateVector(stat->statusPtr, &htil, input->st->length);
  CHECKSTATUSPTR(stat);
  LALCCreateVector(stat->statusPtr, &qprime, input->st->length);
  CHECKSTATUSPTR(stat);
  LALCCreateVector(stat->statusPtr, &conjs, input->st->length);
  CHECKSTATUSPTR(stat);


  /* zero padding the source data */

  for (count1 = 0; count1 < (INT4)(input->h->length); count1++)
  {
      hz->data[count1] = input->h->data[count1];
  }

  for (count2 = count1; count2 < *(input->TotalLength); count2++)
  { 
      hz->data[count2] = 0;
      
  }


  /* taking the fourier transform of the source data */

  LALCreateForwardRealFFTPlan(stat->statusPtr, &plan, *(input->TotalLength), 0);  CHECKSTATUSPTR(stat);
  CHECKSTATUSPTR(stat);
  LALForwardRealFFT(stat->statusPtr, htil, hz, plan); 
  CHECKSTATUSPTR(stat);

  /* taking the conplex conjucate of the transformed template */

  for (count1 = 0; count1 < (INT4)input->st->length; count1++)  
  {    
     conjs->data[count1].re = input->st->data[count1].re;
     conjs->data[count1].im = -(input->st->data[count1].im);
  }


  /* multitplying that result by the transformed source data */
  for (count1 = 0; count1 < (INT4)input->st->length; count1++)  
  {  
     qprime->data[count1].re = (conjs->data[count1].re * htil->data[count1].re)-(conjs->data[count1].im * htil->data[count1].im);
     qprime->data[count1].im = (conjs->data[count1].re * htil->data[count1].im)+(conjs->data[count1].im * htil->data[count1].re); 
  }

  /* taking the inverse fourier transform of the last result */

  LALCreateReverseRealFFTPlan(stat->statusPtr, &plan2, *(input->TotalLength), 0);
  CHECKSTATUSPTR(stat);

  LALReverseRealFFT(stat->statusPtr, largeq, qprime, plan2);  
  CHECKSTATUSPTR(stat);

  /* separating the data into two vectors: one containg the filtered data
  and the other containing the remainder. */

 




  for (count1 = 0; count1 < (INT4)input->h->length; count1++)
      output->q->data[count1] = (1/N)*largeq->data[count1];



  for (count2 = 0; count2 < (INT4)(*(input->TotalLength)-input->h->length); count2++)
  {
      output->z->data[count2] = (1/N)*largeq->data[count1];
      count1++;
      
  }

  /*  Clean up local variables  */
  LALSDestroyVector(stat->statusPtr, &hz);
  CHECKSTATUSPTR(stat);

  LALSDestroyVector(stat->statusPtr, &largeq);
  CHECKSTATUSPTR(stat);

  LALCDestroyVector(stat->statusPtr, &htil); 
  CHECKSTATUSPTR(stat);

  LALCDestroyVector(stat->statusPtr, &qprime);
  CHECKSTATUSPTR(stat);

  LALCDestroyVector(stat->statusPtr, &conjs);
  CHECKSTATUSPTR(stat);

  LALDestroyRealFFTPlan(stat->statusPtr, &plan);
  CHECKSTATUSPTR(stat);

  LALDestroyRealFFTPlan(stat->statusPtr, &plan2);
  CHECKSTATUSPTR(stat); 

  DETATCHSTATUSPTR( stat );

  RETURN( stat );


}





















