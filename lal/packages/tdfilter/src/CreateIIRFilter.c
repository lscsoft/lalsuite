/****************************** <lalVerbatim file="CreateIIRFilterCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{CreateIIRFilter.c}}
\label{ss:CreateIIRFilter.c}

Creates IIR filter objects.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CreateIIRFilterCP}
\idx{LALCreateREAL4IIRFilter()}
\idx{LALCreateREAL8IIRFilter()}

\subsubsection*{Description}

These functions create an object \verb@**output@ of type
\verb@<datatype>IIRFilter@, where \verb@<datatype>@ is \verb@REAL4@ or
\verb@REAL8@.  The filter coefficients are computed from the zeroes,
poles, and gain of an input object \verb@*input@ of type
\verb@COMPLEX8ZPGFilter@ or \verb@COMPLEX16ZPGFilter@, respectively;
the sampling time interval is taken directly from
\verb@input->deltaT@.  The ZPG filter should express the factored
transfer function in the $z=\exp(2\pi if)$ plane.  Initially the
output handle must be a valid handle (\verb@output@$\neq$\verb@NULL@)
but should not point to an existing object
(\verb@*output@=\verb@NULL@)

\subsubsection*{Algorithm}

An IIR filter is a real time-domain filter, which imposes certain
constraints on the zeros, poles, and gain of the filter transfer
function.  The function \verb@Create<datatype>IIRFilter()@ deals with
the constraints either by aborting if they are not met, or by
adjusting the filter response so that they are met.  In the latter
case, warning messages will be issued if the external parameter
\verb@lalDebugLevel@ is set to allow such messages.  The specific
constraints, and how they are dealt with, are as follows:

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
make unstable filters; however, warnings will be issued.  (In some
sense the first condition is a special case of this one, since a
transfer function with more zeros than poles actually has
corresponding poles at infinity.)

Third, the filter must be \emph{physically realizable}; that is, the
transfer function should expand to a rational function of $z$ with
real coefficients.  Necessary and sufficient conditions for this are
that the gain be real, and that all zeros and poles either be real or
come in complex conjugate pairs.  The routine
\verb@Create<datatype>IIRFilter()@ enforces this by using only the
real part of the gain, and only the real or positive-imaginary zeros
and poles; it assumes that the latter are paired with
negative-imaginary conjugates.  The routine will abort if this
assumption results in a change in the given number of zeros or poles.
If \verb@lalDebugLevel@ is set to allow warnings, the routine will
actually check to see that each pair of nonreal poles or zeros are in
fact complex conjugates, and will issue a warning if an unmatched pair
is detected; however, the algorithm will then simply proceed as if the
negative-imaginary points were relocated to the ``correct'' positions.

The code associated with the warning messages is potentially rather
cumbersome for production algorithms; therefore, the value of
\verb@lalDebugLevel@ is tested before performing any other tests
associated with warning messages.  Furthermore, this code block is
surrounded with compiler directives to exclude the code entirely if
the module is compiled with the \verb@NDEBUG@ flag set.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALWarning()                    LALPrintError()
LALMalloc()                     LALFree()
LALSCreateVector()              LALSDestroyVector()
LALDCreateVector()              LALDDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{CreateIIRFilterCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <math.h>
#include <lal/IIRFilter.h>

NRCSID(CREATEIIRFILTERC,"$Id$");

extern INT4 lalDebugLevel;

/* <lalVerbatim file="CreateIIRFilterCP"> */
void LALCreateREAL4IIRFilter( LALStatus         *stat,
			      REAL4IIRFilter    **output,
			      COMPLEX8ZPGFilter *input )
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

  INITSTATUS(stat,"LALCreateREAL4IIRFilter",CREATEIIRFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure all the input structures have been initialized. */
  ASSERT(input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->zeros,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->zeros->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->poles,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->poles->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  numZeros=input->zeros->length;
  numPoles=input->poles->length;
  zeros=input->zeros->data;
  poles=input->poles->data;

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(!*output,stat,IIRFILTERH_EOUT,IIRFILTERH_MSGEOUT);

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

#ifndef NDEBUG
      if(lalDebugLevel&LALWARNING){
	/* Check to see that another zero is an actual conjugate.
           This is not a foolproof test, as it can be fooled by
           multiple zeros at the same location. */
	INT4 j=0;
	INT4 k=0;
	REAL4 x = zeros[i].re - zeros[0].re;
	REAL4 y = zeros[i].im + zeros[0].im;
	REAL4 sep = x*x + y*y;
	for(j=1;j<numZeros;j++){
	  x=zeros[i].re-zeros[j].re;
	  y=zeros[i].im+zeros[j].im;
	  if(sep>x*x+y*y){
	    sep=x*x+y*y;
	    k=j;
	  }
	}
	if(sep>LAL_REAL4_EPS){
	  LALWarning(stat,"Complex zero has no conjugate pair:");
	  LALPrintError("\tUnmatched zero z_%i = %.8e + i*%.8e\n",i,
			zeros[i].re,zeros[i].im);
	  LALPrintError("\tNearest pair   z_%i = %.8e + i*%.8e\n",k,
			zeros[k].re,zeros[k].im);
	}
      }
#endif
    }
  ASSERT(num==numZeros,stat,IIRFILTERH_EPAIR,IIRFILTERH_MSGEPAIR);

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

#ifndef NDEBUG
      if(lalDebugLevel&LALWARNING){
	/* Check to see that another pole is an actual conjugate.
           This is not a foolproof test, as it can be fooled by
           multiple poles at the same location. */
	INT4 j=0;
	INT4 k=0;
	REAL4 x = poles[i].re - poles[0].re;
	REAL4 y = poles[i].im + poles[0].im;
	REAL4 sep = x*x + y*y;
	for(j=1;j<numPoles;j++){
	  x=poles[i].re-poles[j].re;
	  y=poles[i].im+poles[j].im;
	  if(sep>x*x+y*y){
	    sep=x*x+y*y;
	    k=j;
	  }
	}
	if(sep>LAL_REAL4_EPS){
	  LALWarning(stat,"Complex pole has no conjugate pair:");
	  LALPrintError("\tUnmatched pole p_%i = %.8e + i*%.8e\n",i,
			poles[i].re,poles[i].im);
	  LALPrintError("\tNearest pair   p_%i = %.8e + i*%.8e\n",k,
			poles[k].re,poles[k].im);
	}
      }
#endif
    }
  ASSERT(num==numPoles,stat,IIRFILTERH_EPAIR,IIRFILTERH_MSGEPAIR);

#ifndef NDEBUG
  if(lalDebugLevel&LALWARNING){
    /* Issue a warning if the gain is nonreal. */
    if(fabs(input->gain.im)>fabs(LAL_REAL4_EPS*input->gain.re)){
      LALWarning(stat,"Gain is non-real:");
      LALPrintError("\tg = %.8e + i*%.8e\n", input->gain.re,
		    input->gain.im);
    }
    /* Issue a warning if there are any ``removeable'' poles. */
    for(i=0;i<numPoles;i++){
      INT4 j=0;
      for(;j<numZeros;j++)
	if((poles[i].re==zeros[j].re)&&(poles[i].im==zeros[j].im)){
	  LALWarning(stat,"Removeable pole:");
	  LALPrintError("\tp_%i = z_%i = %.8e + i*%.8e\n",i,j,
			poles[i].re,poles[i].im);
	}
    }
    /* Issue a warning if extra factors of 1/z will be applied. */
    if(numPoles<numZeros){
      LALWarning(stat,"Filter has more zeros than poles:");
      LALPrintError("\t%i poles added at complex origin\n",
		    numZeros-numPoles);
    }
    /* Issue a warning if any poles are outside |z|=1. */
    for(i=0;i<numPoles;i++){
      REAL4 zAbs=poles[i].re*poles[i].re+poles[i].im*poles[i].im;
      if(zAbs>1.0){
	LALWarning(stat,"Filter has pole outside of unit circle:");
	LALPrintError("\tp_%i = %.8e + i*%.8e, |p_%i| = %.8e\n",i,
		      poles[i].re,poles[i].im,i,zAbs);
      }
    }
  }
#endif

  /* Everything seems okay, so initialize the filter. */
  *output=(REAL4IIRFilter *)LALMalloc(sizeof(REAL4IIRFilter));
  if ( !(*output) ) {
    ABORT(stat,IIRFILTERH_EMEM,IIRFILTERH_MSGEMEM);
  }
  memset(*output,0,sizeof(REAL4IIRFilter));
  num = (numPoles>=numZeros) ? numZeros : numPoles;
  numDirect+=num;
  numRecurs+=num;

  (*output)->deltaT=input->deltaT;
  LALSCreateVector(stat->statusPtr,&((*output)->directCoef),
		   numDirect);
  BEGINFAIL(stat) {
    LALFree(*output);
    *output=NULL;
  } ENDFAIL(stat);
  LALSCreateVector(stat->statusPtr,&((*output)->recursCoef),
		   numRecurs);
  BEGINFAIL(stat) {
    TRY(LALSDestroyVector(stat->statusPtr,&((*output)->directCoef)),
	stat);
    LALFree(*output);
    *output=NULL;
  } ENDFAIL(stat);
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
  LALSCreateVector(stat->statusPtr,&((*output)->history),num);
  BEGINFAIL(stat) {
    TRY(LALSDestroyVector(stat->statusPtr,&((*output)->directCoef)),
	stat);
    TRY(LALSDestroyVector(stat->statusPtr,&((*output)->recursCoef)),
	stat);
    LALFree(*output);
    *output=NULL;
  } ENDFAIL(stat);
  history=(*output)->history->data;
  for(i=0;i<num;i++)
    history[i]=0.0;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}


/* <lalVerbatim file="CreateIIRFilterCP"> */
void LALCreateREAL8IIRFilter( LALStatus          *stat,
			      REAL8IIRFilter     **output,
			      COMPLEX16ZPGFilter *input )
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

  INITSTATUS(stat,"LALCreateREAL8IIRFilter",CREATEIIRFILTERC);
  ATTATCHSTATUSPTR(stat);

  /* Make sure all the input structures have been initialized. */
  ASSERT(input,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->zeros,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->zeros->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->poles,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(input->poles->data,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  numZeros=input->zeros->length;
  numPoles=input->poles->length;
  zeros=input->zeros->data;
  poles=input->poles->data;

  /* Make sure that the output handle exists, but points to a null
     pointer. */
  ASSERT(output,stat,IIRFILTERH_ENUL,IIRFILTERH_MSGENUL);
  ASSERT(!*output,stat,IIRFILTERH_EOUT,IIRFILTERH_MSGEOUT);

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

#ifndef NDEBUG
      if(lalDebugLevel&LALWARNING){
	/* Check to see that another zero is an actual conjugate.
           This is not a foolproof test, as it can be fooled by
           multiple zeros at the same location. */
	INT4 j=0;
	INT4 k=0;
	REAL8 x = zeros[i].re - zeros[0].re;
	REAL8 y = zeros[i].im + zeros[0].im;
	REAL8 sep = x*x + y*y;
	for(j=1;j<numZeros;j++){
	  x=zeros[i].re-zeros[j].re;
	  y=zeros[i].im+zeros[j].im;
	  if(sep>x*x+y*y){
	    sep=x*x+y*y;
	    k=j;
	  }
	}
	if(sep>LAL_REAL8_EPS){
	  LALWarning(stat,"Complex zero has no conjugate pair:");
	  LALPrintError("\tUnmatched zero z_%i = %.8e + i*%.8e\n",i,
			zeros[i].re,zeros[i].im);
	  LALPrintError("\tNearest pair   z_%i = %.8e + i*%.8e\n",k,
			zeros[k].re,zeros[k].im);
	}
      }
#endif
    }
  ASSERT(num==numZeros,stat,IIRFILTERH_EPAIR,IIRFILTERH_MSGEPAIR);

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

#ifndef NDEBUG
      if(lalDebugLevel&LALWARNING){
	/* Check to see that another pole is an actual conjugate.
           This is not a foolproof test, as it can be fooled by
           multiple poles at the same location. */
	INT4 j=0;
	INT4 k=0;
	REAL8 x = poles[i].re - poles[0].re;
	REAL8 y = poles[i].im + poles[0].im;
	REAL8 sep = x*x + y*y;
	for(j=1;j<numPoles;j++){
	  x=poles[i].re-poles[j].re;
	  y=poles[i].im+poles[j].im;
	  if(sep>x*x+y*y){
	    sep=x*x+y*y;
	    k=j;
	  }
	}
	if(sep>LAL_REAL8_EPS){
	  LALWarning(stat,"Complex pole has no conjugate pair:");
	  LALPrintError("\tUnmatched pole p_%i = %.8e + i*%.8e\n",i,
			poles[i].re,poles[i].im);
	  LALPrintError("\tNearest pair   p_%i = %.8e + i*%.8e\n",k,
			poles[k].re,poles[k].im);
	}
      }
#endif
    }
  ASSERT(num==numPoles,stat,IIRFILTERH_EPAIR,IIRFILTERH_MSGEPAIR);

#ifndef NDEBUG
  if(lalDebugLevel&LALWARNING){
    /* Issue a warning if the gain is nonreal. */
    if(fabs(input->gain.im)>fabs(LAL_REAL8_EPS*input->gain.re)){
      LALWarning(stat,"Gain is non-real:");
      LALPrintError("\tg = %.8e + i*%.8e\n", input->gain.re,
		    input->gain.im);
    }
    /* Issue a warning if there are any ``removeable'' poles. */
    for(i=0;i<numPoles;i++){
      INT4 j=0;
      for(;j<numZeros;j++)
	if((poles[i].re==zeros[j].re)&&(poles[i].im==zeros[j].im)){
	  LALWarning(stat,"Removeable pole:");
	  LALPrintError("\tp_%i = z_%i = %.8e + i*%.8e\n",i,j,
			poles[i].re,poles[i].im);
	}
    }
    /* Issue a warning if extra factors of 1/z will be applied. */
    if(numPoles<numZeros){
      LALWarning(stat,"Filter has more zeros than poles:");
      LALPrintError("\t%i poles added at complex origin\n",
		    numZeros-numPoles);
    }
    /* Issue a warning if any poles are outside |z|=1. */
    for(i=0;i<numPoles;i++){
      REAL8 zAbs=poles[i].re*poles[i].re+poles[i].im*poles[i].im;
      if(zAbs>1.0){
	LALWarning(stat,"Filter has pole outside of unit circle:");
	LALPrintError("\tp_%i = %.8e + i*%.8e, |p_%i| = %.8e\n",i,
		      poles[i].re,poles[i].im,i,zAbs);
      }
    }
  }
#endif

  /* Everything seems okay, so initialize the filter. */
  *output=(REAL8IIRFilter *)LALMalloc(sizeof(REAL8IIRFilter));
  if ( !(*output) ) {
    ABORT(stat,IIRFILTERH_EMEM,IIRFILTERH_MSGEMEM);
  }
  memset(*output,0,sizeof(REAL8IIRFilter));
  num = (numPoles>=numZeros) ? numZeros : numPoles;
  numDirect+=num;
  numRecurs+=num;

  (*output)->deltaT=input->deltaT;
  LALDCreateVector(stat->statusPtr,&((*output)->directCoef),
		   numDirect);
  BEGINFAIL(stat) {
    LALFree(*output);
    *output=NULL;
  } ENDFAIL(stat);
  LALDCreateVector(stat->statusPtr,&((*output)->recursCoef),
		   numRecurs);
  BEGINFAIL(stat) {
    TRY(LALDDestroyVector(stat->statusPtr,&((*output)->directCoef)),
	stat);
    LALFree(*output);
    *output=NULL;
  } ENDFAIL(stat);
  direct=(*output)->directCoef->data;
  recurs=(*output)->recursCoef->data;

  /* Expand the denominator as a polynomial in z. */
  *recurs=-1.0;
  for(i=1;i<numRecurs;i++)
    recurs[i]=0.0;
  for(i=0;i<numPoles;i++){
    INT4 j=numRecurs-1;
    REAL8 x=poles[i].re;
    REAL8 y=poles[i].im;
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
    REAL8 x=zeros[i].re;
    REAL8 y=zeros[i].im;
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
  LALDCreateVector(stat->statusPtr,&((*output)->history),num);
  BEGINFAIL(stat) {
    TRY(LALDDestroyVector(stat->statusPtr,&((*output)->directCoef)),
	stat);
    TRY(LALDDestroyVector(stat->statusPtr,&((*output)->recursCoef)),
	stat);
    LALFree(*output);
    *output=NULL;
  } ENDFAIL(stat);
  history=(*output)->history->data;
  for(i=0;i<num;i++)
    history[i]=0.0;

  /* Normal exit */
  DETATCHSTATUSPTR(stat);
  RETURN(stat);
}
