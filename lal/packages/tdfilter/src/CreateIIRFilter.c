/****************************** <lalVerbatim file="CreateIIRFilterCV">
Author: Creighton, T. D.
$Id$
******************************* </lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{CreateIIRFilter.c}}
\label{ss:CreateIIRFilter.c}

Creates IIR filter objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CreateIIRFilterCP}
\index{\verb&CreateREAL4IIRFilter()&}
\index{\verb&CreateREAL8IIRFilter()&}

\subsubsection*{Description}

These functions create an object \verb@**output@ of type
\verb@<datatype>IIRFilter@, where \verb@<datatype>@ is \verb@REAL4@ or
\verb@REAL8@.  The filter coefficients are computed from the zeroes,
poles, and gain of an input object \verb@*input@ of type
\verb@COMPLEX8ZPGFilter@ or \verb@COMPLEX16ZPGFilter@, respectively.
The ZPG filter should express the factored transfer function in the
$z=\exp(2\pi if)$ plane.  Initially the output handle must be a valid
handle (\verb@output@$\neq$\verb@NULL@) but should not point to an
existing object (\verb@*output@=\verb@NULL@)

\subsubsection*{Algorithm}

An IIR filter is a real time-domain filter, which imposes certain
constraints on the zeros, poles, and gain of the filter transfer
function.  The function \verb@Create<datatype>IIRFilter()@ deals with
the constraints either by aborting if they are not met, or by
adjusting the filter response so that they are met.  In the latter
case, warning messages will be issued if the external parameter
\verb@debuglevel@ is 1 or more.  The specific constraints, and how
they are dealt with, are as follows:

First, the filter must be \emph{causal}; that is, the output at any
time can be generated entirely from the input at previous times.  In
practice this means that the number of (finite) poles in the $z$ plane
must equal or exceed the number of (finite) zeros.  If this is not the
case, \verb@Create<datatype>IIRFilter()@ will add additional poles at
$z=0$.  Effectively this is just adding a delay to the filter response
in order to make it causal.

Second, the filter should be \emph{stable}, which means that all poles
should be located on or within the circle $|z|=1$.  This is not
enforced by \verb@Create<datatype>IIRFilter()@, which can be used to
make unstable filters; however, warnings will be issued if
\verb@debuglevel@ is 1 or more.  (In some sense the first condition is
a special case of this one, since a transfer function with more zeros
than poles actually has corresponding poles at infinity.)

Third, the filter must be \emph{physically realizable}; that is, the
transfer function should expand to a rational function of $z$ with
real coefficients.  Necessary and sufficient conditions for this are
that the gain be real, and that all zeros and poles either be real or
come in complex conjugate pairs.  The routine
\verb@Create<datatype>IIRFilter()@ enforces this by using only the
real part of the gain, and only the real or positive-imaginary zeros
and poles; it assumes that the latter are paired with
negative-imaginary conjugates.  The routine will abort if this
assumption results in a change in the given number of zeros or poles,
but will otherwise simply modify the filter response.  This allows
\verb@debuglevel@=0 runs to proceed without lengthy and usually
unnecessary error trapping; when \verb@debuglevel@ is 1 or more, the
routine checks to make sure that each nonreal zero or pole does in
fact have a complex-conjugate partner.

\subsubsection*{Uses}
\begin{verbatim}
debuglevel
LALMalloc()
SCreateVector()
DCreateVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{CreateIIRFilterCV}}

</lalLaTeX> */

#include "LALStdlib.h"
#include "LALConstants.h"
#include "AVFactories.h"
#include <math.h>
#include "IIRFilter.h"

NRCSID(CREATEIIRFILTERC,"$Id$");

extern INT4 debuglevel;

/* <lalVerbatim file="CreateIIRFilterCP"> */
void CreateREAL4IIRFilter(Status            *stat,
			  REAL4IIRFilter    **output,
			  COMPLEX8ZPGFilter *input)
{ /* </lalVerbatim> */
  INT4 i;          /* Index counter for zeros and poles. */
  INT4 numZeros;   /* The number of zeros. */
  INT4 numPoles;   /* The number of poles. */
  INT4 numDirect;  /* The number of direct filter coefficients. */
  INT4 numRecurs;  /* The number of recursive filter coefficients. */
  INT4 num;        /* An extra counter for error checking. */
  COMPLEX8 *zeros; /* The zeros of the transfer function. */
  COMPLEX8 *poles; /* The poles of the transfer function. */
  REAL4 *direct;   /* The direct filter coefficients. */
  REAL4 *recurs;   /* The recursive filter coefficients. */
  REAL4 *history;  /* The filter history. */

  INITSTATUS(stat,"CreateREAL4IIRFilter",CREATEIIRFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure all the input structures have been initialized. */
  ASSERT(input,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(input->zeros,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(input->zeros->data,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(input->poles,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(input->poles->data,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  numZeros=input->zeros->length;
  numPoles=input->poles->length;
  zeros=input->zeros->data;
  poles=input->poles->data;

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(!*output,stat,IIRFILTER_EOUT,IIRFILTER_MSGEOUT);

  /* Check that zeros are appropriately paired.  Also, keep track of
     the number of zeros at z=0, since these will reduce the number of
     direct coefficients required. */
  numDirect=1;
  for(i=0,num=0;i<numZeros;i++)
    if(zeros[i].im==0.0){
      num+=1;
      if(zeros[i].re==0.0)
	numDirect-=1;
    }
    else if(zeros[i].im>0.0){
      num+=2;
      if(debuglevel>0){
	/* Check to see that another zero is an actual conjugate.
           This is not a foolproof test, as it can be fooled by
           multiple zeros at the same location. */
	BOOLEAN ok=0;
	INT4 j=0;
	for(;j<numZeros;j++)
	  if((zeros[j].re==zeros[i].re)&&(zeros[j].im==-zeros[i].im))
	    ok=1;
	if(!ok){
	  LALPrintError("Warning: createIIRFilter: zero number %i"
			" has no complex conjugate pair.\n",i);
	  if(debuglevel>1)
	    LALPrintError("         Zero location:      %.8e +"
			  " i*%.8e\n",zeros[i].re,zeros[i].im);
	  if(debuglevel>2){
	    /* Find the nearest approximation to a conjugate, to check
               for possible typos or roundoff errors. */
            REAL4 x = zeros[i].re - zeros[0].re;
            REAL4 y = zeros[i].im + zeros[0].im;
            REAL4 sep = x*x + y*y;
            INT4 k = 0;
	    for(j=1;j<numZeros;j++){
	      x=zeros[i].re-zeros[j].re;
	      y=zeros[i].im+zeros[j].im;
	      if(sep>x*x+y*y){
		sep=x*x+y*y;
		k=j;
	      }
	    }
	    LALPrintError("         Nearest conjugate:  %.8e +"
			  " i*%.8e\n",zeros[k].re,zeros[k].im);
	  }
	}
      }
    }
  ASSERT(num==numZeros,stat,IIRFILTER_EPAIR,IIRFILTER_MSGEPAIR);

  /* Check that poles are appropriately paired.  Also, keep track of
     the number of poles at z=0, since these will reduce the number of
     recursive coefficients required. */
  numRecurs=1;
  for(i=0,num=0;i<numPoles;i++)
    if(poles[i].im==0.0){
      num+=1;
      if(poles[i].re==0.0)
	numRecurs-=1;
    }
    else if(poles[i].im>0.0){
      num+=2;
      if(debuglevel>0){
	/* Check to see that another zero is an actual conjugate.
           This is not a foolproof test, as it can be fooled by
           multiple poles at the same location. */
	BOOLEAN ok=0;
	INT4 j=0;
	for(;j<numPoles;j++)
	  if((poles[j].re==poles[i].re)&&(poles[j].im==-poles[i].im))
	    ok=1;
	if(!ok){
	  LALPrintError("Warning: createIIRFilter: pole number %i"
			" has no complex conjugate pair.\n",i);
	  if(debuglevel>1)
	    LALPrintError("         Pole location:      %.8e +"
			  " i*%.8e\n",poles[i].re,poles[i].im);
	  if(debuglevel>2){
	    /* Find the nearest approximation to a conjugate, to check
               for possible typos or roundoff errors. */
	    REAL4 x=poles[i].re-poles[0].re;
	    REAL4 y=poles[i].im+poles[0].im;
	    REAL4 sep=x*x + y*y;
	    INT4 k = 0;
	    for(j=1;j<numPoles;j++){
	      x=poles[i].re-poles[j].re;
	      y=poles[i].im+poles[j].im;
	      if(sep>x*x+y*y){
		sep=x*x+y*y;
		k=j;
	      }
	    }
	    LALPrintError("         Nearest conjugate:  %.8e +"
			  " i*%.8e\n",poles[k].re,poles[k].im);
	  }
	}
      }
    }
  ASSERT(num==numPoles,stat,IIRFILTER_EPAIR,IIRFILTER_MSGEPAIR);

  if(debuglevel>0){
    /* Issue a warning if the gain is nonreal. */
    if(input->gain.im!=0.0){
      LALPrintError("Warning: createIIRFilter: gain is non-real.\n");
      if(debuglevel>1)
	LALPrintError("         Value: %.8e + i*%.8e\n",
		      input->gain.re,input->gain.im);
    }
    /* Issue a warning if there are any ``removeable'' poles. */
    for(i=0;i<numPoles;i++){
      INT4 j=0;
      for(;j<numZeros;j++)
	if((poles[i].re==zeros[j].re)&&(poles[i].im==zeros[j].im)){
	  LALPrintError("Warning: createIIRFilter: pole number %i"
			" overlaps with zero number %i\n",i,j);
	  if(debuglevel>1)
	    LALPrintError("         Location: %.8e + i*%.8e\n",
			  poles[i].re,poles[i].im);
	}
    }
    /* Issue a warning if extra factors of 1/z will be applied. */
    if(numPoles<numZeros)
      LALPrintError("Warning: createIIRFilter: %i more zeros than"
		    " poles\n",numZeros-numPoles);
    /* Issue a warning if any poles are outside |z|=1. */
    for(i=0;i<numPoles;i++)
      if(poles[i].re*poles[i].re+poles[i].im*poles[i].im>1.0){
	LALPrintError("Warning: createIIRFilter: pole number %i lies"
		      " outside |z|=1\n",i);
	if(debuglevel>1)
	  LALPrintError("         Location: %.8e + i*%.8e\n",
			poles[i].re,poles[i].im);
      }
  }

  /* Everything seems okay, so initialize the filter. */
  *output=(REAL4IIRFilter *)LALMalloc(sizeof(REAL4IIRFilter));
  ASSERT(*output,stat,IIRFILTER_EMEM,IIRFILTER_MSGEMEM);
  memset(*output,0,sizeof(REAL4IIRFilter));
  if(numPoles>=numZeros)
    num=numZeros;
  else
    num=numPoles;
  numDirect+=num;
  numRecurs+=num;

  TRY(SCreateVector(stat->statusPtr,&((*output)->directCoef),numDirect),
      stat);
  TRY(SCreateVector(stat->statusPtr,&((*output)->recursCoef),numRecurs),
      stat);
  direct=(*output)->directCoef->data;
  recurs=(*output)->recursCoef->data;

  /* Expand the denominator as a polynomial in z. */
  *recurs=-1.0;
  for(i=1;i<numRecurs;i++)
    recurs[i]=0.0;
  for(i=0;i<numPoles;i++){
    INT4 j=numRecurs-1;
    REAL4 x=poles[i].re;
    REAL4 y=poles[i].im;
    if(y==0.0)
      for(;j>0;j--)
	recurs[j]-=x*recurs[j-1];
    else if(y>0.0){
      for(;j>1;j--)
	recurs[j]-=2.0*x*recurs[j-1]-(x*x+y*y)*recurs[j-2];
      recurs[j]-=2.0*x*recurs[j-1];
    }
  }

  /* Expand the numerator as a polynomial in z. */
  for(i=0;i<numDirect;i++)
    direct[i]=0.0;
  direct[num-numZeros]=input->gain.re;
  for(i=0;i<numZeros;i++){
    INT4 j=numDirect-1;
    REAL4 x=zeros[i].re;
    REAL4 y=zeros[i].im;
    if(y==0.0)
      for(;j>num-numZeros;j--)
	direct[j]-=x*direct[j-1];
    else if(y>0.0){
      for(;j>num-numZeros+1;j--)
	direct[j]-=2.0*x*direct[j-1]-(x*x+y*y)*direct[j-2];
      direct[j]-=2.0*x*direct[j-1];
    }
  }

  /* Initialize filter history. */
  if(numDirect>=numRecurs)
    num=numDirect-1;
  else
    num=numRecurs-1;
  TRY(SCreateVector(stat->statusPtr,&((*output)->history),num),stat);
  history=(*output)->history->data;
  for(i=0;i<num;i++)
    history[i]=0.0;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="CreateIIRFilterCP"> */
void CreateREAL8IIRFilter(Status             *stat,
			  REAL8IIRFilter     **output,
			  COMPLEX16ZPGFilter *input)
{ /* </lalVerbatim> */
  INT4 i;           /* Index counter for zeros and poles. */
  INT4 numZeros;    /* The number of zeros. */
  INT4 numPoles;    /* The number of poles. */
  INT4 numDirect;   /* The number of direct filter coefficients. */
  INT4 numRecurs;   /* The number of recursive filter coefficients. */
  INT4 num;         /* An extra counter for error checking. */
  COMPLEX16 *zeros; /* The zeros of the transfer function. */
  COMPLEX16 *poles; /* The poles of the transfer function. */
  REAL8 *direct;    /* The direct filter coefficients. */
  REAL8 *recurs;    /* The recursive filter coefficients. */
  REAL8 *history;   /* The filter history. */

  INITSTATUS(stat,"CreateREAL8IIRFilter",CREATEIIRFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure all the input structures have been initialized. */
  ASSERT(input,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(input->zeros,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(input->zeros->data,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(input->poles,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(input->poles->data,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  numZeros=input->zeros->length;
  numPoles=input->poles->length;
  zeros=input->zeros->data;
  poles=input->poles->data;

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,IIRFILTER_ENUL,IIRFILTER_MSGENUL);
  ASSERT(!*output,stat,IIRFILTER_EOUT,IIRFILTER_MSGEOUT);

  /* Check that zeros are appropriately paired.  Also, keep track of
     the number of zeros at z=0, since these will reduce the number of
     direct coefficients required. */
  numDirect=1;
  for(i=0,num=0;i<numZeros;i++)
    if(zeros[i].im==0.0){
      num+=1;
      if(zeros[i].re==0.0)
	numDirect-=1;
    }
    else if(zeros[i].im>0.0){
      num+=2;
      if(debuglevel>0){
	/* Check to see that another zero is an actual conjugate. */
	BOOLEAN ok=0;
	INT4 j=0;
	for(;j<numZeros;j++)
	  if((zeros[j].re==zeros[i].re)&&(zeros[j].im==-zeros[i].im))
	    ok=1;
	if(!ok){
	  LALPrintError("Warning: createIIRFilter: zero number %i"
			" has no complex conjugate pair.\n",i);
	  if(debuglevel>1)
	    LALPrintError("         Zero location:      %.16e +"
			  " i*%.16e\n",zeros[i].re,zeros[i].im);
	  if(debuglevel>2){
	    /* Find the nearest approximation to a conjugate, to check
               for possible typos or roundoff errors. */
	    REAL8 x=zeros[i].re-zeros[0].re;
	    REAL8 y=zeros[i].im+zeros[0].im;
	    REAL8 sep=x*x + y*y;
	    INT4 k = 0;
	    for(j=1;j<numZeros;j++){
	      x=zeros[i].re-zeros[j].re;
	      y=zeros[i].im+zeros[j].im;
	      if(sep>x*x+y*y){
		sep=x*x+y*y;
		k=j;
	      }
	    }
	    LALPrintError("         Nearest conjugate:  %.16e +"
			  " i*%.16e\n",zeros[k].re,zeros[k].im);
	  }
	}
      }
    }
  ASSERT(num==numZeros,stat,IIRFILTER_EPAIR,IIRFILTER_MSGEPAIR);

  /* Check that poles are appropriately paired.  Also, keep track of
     the number of poles at z=0, since these will reduce the number of
     recursive coefficients required. */
  numRecurs=1;
  for(i=0,num=0;i<numPoles;i++)
    if(poles[i].im==0.0){
      num+=1;
      if(poles[i].re==0.0)
	numRecurs-=1;
    }
    else if(poles[i].im>0.0){
      num+=2;
      if(debuglevel>0){
	/* Check to see that another zero is an actual conjugate. */
	BOOLEAN ok=0;
	INT4 j=0;
	for(;j<numPoles;j++)
	  if((poles[j].re==poles[i].re)&&(poles[j].im==-poles[i].im))
	    ok=1;
	if(!ok){
	  LALPrintError("Warning: createIIRFilter: pole number %i"
			" has no complex conjugate pair.\n",i);
	  if(debuglevel>1)
	    LALPrintError("         Pole location:      %.16e +"
			  " i*%.16e\n",poles[i].re,poles[i].im);
	  if(debuglevel>2){
	    /* Find the nearest approximation to a conjugate, to check
               for possible typos or roundoff errors. */
	    REAL8 x=poles[i].re-poles[0].re;
	    REAL8 y=poles[i].im+poles[0].im;
	    REAL8 sep=x*x + y*y;
	    INT4 k = 0;
	    for(j=1;j<numPoles;j++){
	      x=poles[i].re-poles[j].re;
	      y=poles[i].im+poles[j].im;
	      if(sep>x*x+y*y){
		sep=x*x+y*y;
		k=j;
	      }
	    }
	    LALPrintError("         Nearest conjugate:  %.16e +"
			  " i*%.16e\n",poles[k].re,poles[k].im);
	  }
	}
      }
    }
  ASSERT(num==numPoles,stat,IIRFILTER_EPAIR,IIRFILTER_MSGEPAIR);

  if(debuglevel>0){
    /* Issue a warning if the gain is nonreal. */
    if(input->gain.im!=0.0){
      LALPrintError("Warning: createIIRFilter: gain is non-real.\n");
      if(debuglevel>1)
	LALPrintError("         Value: %.16e + i*%.16e\n",
		      input->gain.re,input->gain.im);
    }
    /* Issue a warning if there are any ``removeable'' poles. */
    for(i=0;i<numPoles;i++){
      INT4 j=0;
      for(;j<numZeros;j++)
	if((poles[i].re==zeros[j].re)&&(poles[i].im==zeros[j].im)){
	  LALPrintError("Warning: createIIRFilter: pole number %i"
			" overlaps with zero number %i\n",i,j);
	  if(debuglevel>1)
	    LALPrintError("         Location: %.16e + i*%.16e\n",
			  poles[i].re,poles[i].im);
	}
    }
    /* Issue a warning if extra factors of 1/z will be applied. */
    if(numPoles<numZeros)
      LALPrintError("Warning: createIIRFilter: %i more zeros than"
		    " poles\n",numZeros-numPoles);
    /* Issue a warning if any poles are outside |z|=1. */
    for(i=0;i<numPoles;i++)
      if(poles[i].re*poles[i].re+poles[i].im*poles[i].im>1.0){
	LALPrintError("Warning: createIIRFilter: pole number %i lies"
		      " outside |z|=1\n",i);
	if(debuglevel>1)
	  LALPrintError("         Location: %.16e + i*%.16e\n",
			poles[i].re,poles[i].im);
      }
  }

  /* Everything seems okay, so initialize the filter. */
  *output=(REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter));
  ASSERT(*output,stat,IIRFILTER_EMEM,IIRFILTER_MSGEMEM);
  memset(*output,0,sizeof(REAL8IIRFilter));
  if(numPoles>=numZeros)
    num=numZeros;
  else
    num=numPoles;
  numDirect+=num;
  numRecurs+=num;

  TRY(DCreateVector(stat->statusPtr,&((*output)->directCoef),numDirect),
      stat);
  TRY(DCreateVector(stat->statusPtr,&((*output)->recursCoef),numRecurs),
      stat);
  direct=(*output)->directCoef->data;
  recurs=(*output)->recursCoef->data;
  for(i=0;i<numDirect;i++)
    direct[i]=0.0;
  direct[num-numZeros]=input->gain.re;
  *recurs=-1.0;
  for(i=1;i<numRecurs;i++)
    recurs[i]=0.0;

  /* Expand the numerator as a polynomial in z. */
  for(i=0;i<numZeros;i++){
    INT4 j=num-numZeros+1;
    REAL8 x=zeros[i].re;
    REAL8 y=zeros[i].im;
    if(y==0.0)
      for(;j<numDirect;j++)
	direct[j]-=x*direct[j-1];
    else if(y>0.0){
      direct[j]-=2.0*x*direct[j-1];
      for(j++;j<numDirect;j++)
	direct[j]-=x*direct[j-1]-(x*x+y*y)*direct[j-2];
    }
  }

  /* Expand the denominator as a polynomial in z. */
  for(i=0;i<numPoles;i++){
    INT4 j=1;
    REAL8 x=poles[i].re;
    REAL8 y=poles[i].im;
    if(y==0.0)
      for(;j<numRecurs;j++)
	recurs[j]-=x*recurs[j-1];
    else if(y>0.0){
      recurs[j]-=2.0*x*recurs[j-1];
      for(j++;j<numRecurs;j++)
	recurs[j]-=x*recurs[j-1]-(x*x+y*y)*recurs[j-2];
    }
  }

  /* Initialize filter history. */
  if(numDirect>=numRecurs)
    num=numDirect-1;
  else
    num=numRecurs-1;
  TRY(DCreateVector(stat->statusPtr,&((*output)->history),num),stat);
  history=(*output)->history->data;
  for(i=0;i<num;i++)
    history[i]=0.0;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}
