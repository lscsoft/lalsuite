/************************************ <lalVerbatim file="LALDemodCV">
Author: Berukoff, S.J.,  Papa, M.A., Allen, B. $Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>
\subsection{Module \texttt{LALDemod.c}}\label{ss:LALDemod.c}
Computes a demodulated Fourier transform (DeFT) given a set of input short Fourier transforms (SFT).

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALDemodD}
\index{\texttt{LALDemod()}}

\subsubsection*{Description}
As described in sections \ref{demodulation_intro} and \ref{s:LALDemod.h}, the
demodulation routine takes a set of SFTs and combines them into a longer time
baseline, higher frequency resolution, transform.  If the constraints on the time baseline of the
SFTs with respect to the phase evolution of the SFT signal are met, then the total
power of the signal is concentrated, within a few percent, in one or two
frequency bins.  The normalization adopted here is the following: suppose a
signal that appears as monochromatic during the SFT time baseline has power
(here defined as $|SFT|^2$) $\cal{A}$, at the "instantaneous" signal's
frequency bin\footnote{This frequency bin in general is different for
different SFTs - the signal's instantaneous frequency changes in time.}.  The
power observed after demodulation, at the signal's intrinsic frequency,
$|DeFT|^2$, will be $M\cal{A}$, with $M$ being the number of SFTs that have
been combined together in the DeFT. These normalization conditions are valid
within a factor of $\sim 2$ due to frequency discretization effects.

The {\bf parameter structure} contains the necessary information about the
frequency band of the DeFT (\verb@if0Max@ and \verb@if0Min@), the frequency
band of the input SFTs (\verb@ifMin@), how many SFTs have to be combined to
make a DeFT (\verb@mCohSFT@) and the source template parameters to be used in
the demodulation (\verb@*spinDwn and @\verb@*skyConst@).  The parameter
structure contains two other fields other than the ones mentioned above:
\verb@mObsCoh@ and \verb@iCoh@.  These are superfluous if one wants to produce
only a single DeFT, as the former is unity and the latter, zero.  However, if
one wants to produce more than one, these may be useful parameters to manage
input and output data (see below).  \verb@mObsCoh@ is the total number of
DeFTs that one wants to produce and \verb@iCoh@ labels the order number of the
one currently under production.  Finally, \verb@amcoe@ contains the values of the amplitude modulation functions $a$ and $b$.

The {\bf input} is: \verb@**input@, an array of structures of type \verb@FFT@,
containing not only the data itself but quantities relevant to the search.  It
is composed of a LAL data structure \verb@COMPLEX8FrequencySeries *@ and a
\verb@ParameterSet@ structure; the latter allows one to store specific
information about an FFT.  Note further that the
\verb@COMPLEX8FrequencySeries@ structure particularly contains
\begin{description}
\item[\texttt{input[k]->fft->data->data[i]}]	${\rm{i^{th}}}$ complex data point of the ${\rm{k^{th}}}$ SFT.
\item[\texttt{input[k]->fft->data->length}]		Length of data set.
\end{description}
These are perhaps the most oft-referenced parts of these structures.  

The {\bf output} is a pointer to the DeFT, '\verb@xHat@',
that has been produced.  This is a structure of type \verb@DeFTPeriodogram **@.  In general \verb@xHat@ may hold
more than a single DeFT.  By knowing  \verb@mObsCoh@ and \verb@iCoh@ the
demodulation routine will fill up the appropriate part of \verb@xHat@. 

\subsubsection*{Algorithm}

The routine implements the analytical result of eq. \ref{DeFT_algo}.  It thus
uses a nested-loop structure, which computes the frequency components of the
DeFT, and then their square moduli, with each loop, placing them in an array (\verb@xHat@) for future use.
The two outer most loops loop over two parameters $l$ and $\beta$ ($0<l<N-1$
and $0<\beta <M-1$) which comprise the DeFT frequency index $b=\beta+Ml$
(following the notation of section \ref{s:LALDemod.h}).  The next loop is over
$\alpha$, which identifies the SFTs to be used in order to produce the DeFT.
The value of  $k^*$ is then computed using the second of Eq. \ref{DeFT_defs},
and thus the summation over $k$ of Eq.\ref{DeFT_algo} is carried out, with a
loop over NTERMS.  In this loop the product $\tilde{x}_{\alpha k}P_{\alpha k}$
is calculated.  Once this loop completes, 
$e^{iy_\alpha}$ is computed, the summation over $\alpha$ performed and,
finally,  the code yields the DeFT $\hat{x}_b$.  It can be seen that the code
closely follows the analytical development of the formalism. 

Finally, note that in order to avoid repeated
trigonometric function computations, a look-up-table (LUT) for sine and cosine is 
constructed at the beginning of the routine. 

\subsubsection*{Uses}
\begin{verbatim}
None
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALDemodCV}}

</lalLaTeX> */

/* loop protection */
#ifndef LALDEMOD_C
#define LALDEMOD_C
#endif

#include "LALDemod.h"


NRCSID(LALDEMODC, "$Id$");

/* <lalVerbatim file="LALDemodD"> */

void LALDemod(LALStatus *stat, 
	DeFTPeriodogram *xHat, 
	FFT **input, 
	DemodPar *params) 
{ 
/* </lalVerbatim> */

  INT4 l,beta,alpha;  /* loop indices */
  INT4 xHatIndex=-1; /* indices for DeFFT vector */
  INT4 alpha1, alpha2;	/* indices labelling SFT being worked on */
  REAL8	*xSum=NULL, *ySum=NULL;	/* temp variables for computation of fs*as and fs*bs */
  INT4 s=0;		/* local variable for spinDwn calcs. */
  REAL8	xTemp;	/* temp variable for phase model */
  REAL8	deltaF;	/* width of SFT band */
  INT4	k, k1;	/* defining the sum over which input->fft->data->data*P is calculated */
  REAL8 	*skyConst;	/* vector of sky constants data */
  REAL8 	*spinDwn;	/* vector of spinDwn parameters (maybe a structure? */
  INT4	spOrder;	/* maximum spinDwn order */
  REAL8	x;		/* local variable for holding x */
  REAL8	realTemp, imagTemp;	/* temp variables used in computation of input->fft->data->data*P */
  REAL8	realP, imagP;	/* real and imaginary parts of P, see CVS */
  INT4	nDeltaF;	/* number of frequency bins per SFT band */
  INT4	sftIndex;	/* more temp variables */
  REAL8	y;		/* local variable for holding y */
  REAL8 realQ, imagQ;
  INT4 	iCoh;

  REAL8 *sinVal,*cosVal; /*LUT values computed by the routine do_trig_lut*/
  INT4  res;           /*resolution of the argument of the trig functions: 2pi/res.*/

  REAL8 sqrtMCohSFT, oneOverSqrtMCohSFT, oneOverMCohSFT;
  INT4 tempMCohSFT, tempif0Max, tempif0Min, tempifMin;
  INT4 *tempInt1;
  INT4 index;
  COMPLEX16Vector *tempXhata = NULL;	  
  COMPLEX16Vector *tempXhatb = NULL;
  REAL8Vector *fA = NULL;
  REAL8Vector *fB = NULL;
  REAL8Vector *fAB = NULL;

  INITSTATUS(stat, "LALDemod", LALDEMODC);
  ATTATCHSTATUSPTR(stat);
  
  /* Make sure there are SFTs to work on  */
  ASSERT((*input)!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);          	/* the input structure pointer */
  ASSERT((*input)->fft!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);     	/* the fft (REAL8Fr.Ser.) sub-structure */
  ASSERT((*input)->fft->data!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);	/* the data (REAL8Vec.) sub-strucure */
  ASSERT((*input)->fft->data->data!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL); /* the data itself (REAL8) */
  
  /* Make sure parameter pointers are pointing somewhere */
  ASSERT(params!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);               /* the entire parameter structure */
  ASSERT(params->skyConst!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);     /* the sky constants, a and b */
  ASSERT(params->spinDwn!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);      /* the spindown constants */
  
  /* Make sure parameters are reasonable */
  ASSERT(params->if0Min>=0, stat, LALDEMODH_ENNEG, LALDEMODH_MSGENNEG);
  ASSERT(params->if0Max>=0, stat, LALDEMODH_ENNEG, LALDEMODH_MSGENNEG);
  ASSERT(params->if0Min<params->if0Max, stat, LALDEMODH_EORDR, LALDEMODH_MSGEORDR);
  ASSERT(params->mCohSFT>0, stat, LALDEMODH_ENNEG, LALDEMODH_MSGENNEG);
  ASSERT(params->mObsCoh>0, stat, LALDEMODH_ENNEG, LALDEMODH_MSGENNEG);
  ASSERT(params->iCoh>=0, stat, LALDEMODH_ENNEG, LALDEMODH_MSGENNEG);
  ASSERT(params->spinDwnOrder>=0, stat, LALDEMODH_ENNEG, LALDEMODH_MSGENNEG);
  
  /* Make sure input structure parameters are reasonable */
  ASSERT((*input)->fft->data->length>0, stat, LALDEMODH_ENNEG, LALDEMODH_MSGENNEG);
  ASSERT((*input)->fft->deltaF>0, stat, LALDEMODH_ENNEG, LALDEMODH_MSGENNEG);
  
  /* variable redefinitions for code readability */
  spOrder=params->spinDwnOrder;
  spinDwn=params->spinDwn;
  skyConst=params->skyConst;
  deltaF=(*input)->fft->deltaF;
  nDeltaF=(*input)->fft->data->length;
  iCoh=params->iCoh;
  sqrtMCohSFT=sqrt((REAL8)params->mCohSFT);
  oneOverSqrtMCohSFT=1.0/sqrtMCohSFT;
  oneOverMCohSFT=1.0/(REAL8)params->mCohSFT;
  tempMCohSFT=params->mCohSFT;
  tempif0Max=params->if0Max;
  tempif0Min=params->if0Min;
  tempifMin=params->ifMin;

  /* CHECK ALLOCATION, SO WE HAVE SOMETHING TO WRITE TO */
  ASSERT(xHat!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);
  ASSERT(xHat->fft!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);
  ASSERT(xHat->fft->data!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);
 
  
  /* Allocate memory for the temporary xHats */
  LALZCreateVector(stat->statusPtr, &(tempXhata), xHat->fft->data->length);
  LALZCreateVector(stat->statusPtr, &(tempXhatb), xHat->fft->data->length);

  /* Allocate memory for the individual statistics */
  LALDCreateVector(stat->statusPtr, &(fA), xHat->fft->data->length);
  LALDCreateVector(stat->statusPtr, &(fB), xHat->fft->data->length);
  LALDCreateVector(stat->statusPtr, &(fAB), xHat->fft->data->length);
  /* Allocate memory for LUTs. These are stored in "sinVal" and in "cosVal". 
     The integer "res" defines the resolution of the argument of the trig functions: 2pi/res. 
     The size of the LUT is res+1 in order to protect against rounding errors in the index 
     calculation.*/

  /* res=10*(params->mCohSFT); */
  /* This size LUT gives errors ~ 10^-7 with a three-term Taylor series */
   res=64;
   sinVal=(REAL8 *)LALMalloc((res+1)*sizeof(REAL8));
   cosVal=(REAL8 *)LALMalloc((res+1)*sizeof(REAL8)); 
   for (k=0; k<=res; k++){
     sinVal[k]=sin((LAL_TWOPI*k)/res);
     cosVal[k]=cos((LAL_TWOPI*k)/res);
   }

  /* setting some values */
  alpha1=iCoh*tempMCohSFT;
  alpha2=alpha1+tempMCohSFT;
  /*
    tempXhat=(COMPLEX16 *)LALMalloc((tempif0Max-tempif0Min)*tempMCohSFT*sizeof(COMPLEX16));
  */
  /* this loop computes the values of the phase model */
  xSum=(REAL8 *)LALMalloc(tempMCohSFT*sizeof(REAL8));
  ySum=(REAL8 *)LALMalloc(tempMCohSFT*sizeof(REAL8));
  tempInt1=(INT4 *)LALMalloc(tempMCohSFT*sizeof(INT4));
  for(alpha=alpha1;alpha<alpha2;alpha++){
    tempInt1[alpha-alpha1]=2*alpha*(spOrder+1)+1;
    xSum[alpha-alpha1]=0.0;
    ySum[alpha-alpha1]=0.0;
    for(s=0; s<spOrder;s++) {
      xSum[alpha-alpha1] += spinDwn[s] * skyConst[tempInt1[alpha-alpha1]+2+2*s]; 	
      ySum[alpha-alpha1] += spinDwn[s] * skyConst[tempInt1[alpha-alpha1]+1+2*s];
    }
  }

   /* the loops begin */
  for(l=tempif0Min;l<tempif0Max+1;l++)
    {
      for(beta=0;beta<tempMCohSFT;beta++)
	{
	  COMPLEX16 tXa;
	  COMPLEX16 tXb;
	  
	  xHatIndex++;
  	  tXa.re = tempXhata->data[xHatIndex].re = 0.0; 
  	  tXa.im = tempXhata->data[xHatIndex].im = 0.0; 
  	  tXb.re = tempXhatb->data[xHatIndex].re = 0.0; 
  	  tXb.im = tempXhatb->data[xHatIndex].im = 0.0;	   
	  
	  for(alpha=alpha1;alpha<alpha2;alpha++)
	    {
	      REAL8 tsin, tcos, tempFreq;
	      COMPLEX8 *tempPtr=input[alpha]->fft->data->data;
	      REAL4 a = params->amcoe->a->data[alpha];
	      REAL4 b = params->amcoe->b->data[alpha];


	      /* This looks complicated - and it is. Look at the documentation, specifically, 
		 look at eq. \ref{DeFT_algo}. We have tried to use the same names for the 
		 variables in this routine and in the equations of the documentation. */
	      
	      xTemp=((REAL8)beta*oneOverMCohSFT+(REAL8)l)*
		skyConst[tempInt1[alpha-alpha1]]*deltaF+
		xSum[alpha-alpha1];
	      
	      realTemp=0.0;
	      imagTemp=0.0;
	      
	      /* find correct index into LUT -- pick closest point */
	      tempFreq=xTemp-(INT4)xTemp;
	      index=(INT4)(tempFreq*res+0.5);
	      /*     tsin=sinVal[index];
		     tcos=cosVal[index]-1.0;
	      */
	      
	      {
		/* Taylor series approximant for better accuracy with small LUT:
		   sin(x+d)=sin(x)+d cos(x) -1/2 d^2 sin(x) -1/6 d^3 cos(x)
		   cos(x+d)=cos(x)-d sin(x) -1/2 d^2 cos(x) +1/6 d^3 sin(x)
		*/
		REAL8 d=LAL_TWOPI*(tempFreq-(REAL8)index/(REAL8)res);
		REAL8 d2=0.5*d*d;
		/* REAL8 d3=0.3333333333333333333333333333333*d2*d; */
		REAL8 ts=sinVal[index];
		REAL8 tc=cosVal[index];
		
		tsin=ts+d*tc-d2*ts /* -d3*tc */;
		tcos=tc-d*ts-d2*tc /* +d3*ts */-1.0;
	     }

	     /* Should yield same result as the LUT defined above*/
	     /* tsin=sin(LAL_TWOPI*tempFreq+SMALL); */
	     /* tcos=cos(LAL_TWOPI*tempFreq)-1.0; */
	     
	     tempFreq=LAL_TWOPI*(tempFreq+NTERM_COH_DIV_TWO-1);
	     k1=(INT4)xTemp-NTERM_COH_DIV_TWO+1;
	     for(k=0;k<2*NTERM_COH_DIV_TWO;k++)
	     {
		 COMPLEX8 tempPtr2;

		 
		 x=tempFreq-LAL_TWOPI*(REAL8)k;
		 /* Note that this is a problem if x is small and negative!!! */
		 realP=(tsin)/(x+SMALL);
		 imagP=tcos/(x+SMALL);
		 		 
		 sftIndex=k1+k-tempifMin;

		 /* these four lines compute P*xtilde */
		 tempPtr2=tempPtr[sftIndex];
		 realTemp += tempPtr2.re*realP;
		 realTemp -= tempPtr2.im*imagP;
		 imagTemp += tempPtr2.re*imagP;
		 imagTemp += tempPtr2.im*realP;
	     }

	     y=-LAL_TWOPI*(((REAL8)beta*oneOverMCohSFT+(REAL8)l) * 
			   skyConst[tempInt1[alpha-alpha1]-1]*deltaF+
			   ySum[alpha-alpha1]);

	     realQ = cos(y);
	     imagQ = sin(y);

	     /* implementation of amplitude demodulation */
	     {
	       REAL8 temp1 = realTemp*realQ-imagTemp*imagQ;
	       REAL8 temp2 = realTemp*imagQ+imagTemp*realQ;
	       tXa.re += a*temp1;
	       tXa.im += a*temp2;
	       tXb.re += b*temp1;
	       tXb.im += b*temp2;
	     }
	    
	    }
	  
	  /* Now, compute normalised periodograms for each statistic */
	  /* |F_a|^2, |F_b|^2, Re[(F_a)(F_b}^*] */
	  
	  fA->data[xHatIndex] = oneOverSqrtMCohSFT*(tXa.re*tXa.re+tXa.im*tXa.im);
	  fB->data[xHatIndex] = oneOverSqrtMCohSFT*(tXb.re*tXb.re+tXb.im*tXb.im);
	  fAB->data[xHatIndex] = oneOverSqrtMCohSFT*(tXa.re*tXb.re+tXa.im*tXb.im);
	  
	  /* Compute statistic for this DeFT bin */
	  /* F = 1/D * (B*|F_a|^2+A*|F_b|^2-2C*Re(F_a * F_b^*)) */
	  xHat->fft->data->data[xHatIndex] = (1.0/params->amcoe->D) * 
	    ((params->amcoe->B)*fA->data[xHatIndex] + 
	    (params->amcoe->A)*fB->data[xHatIndex] -
	    2*(params->amcoe->C)*fAB->data[xHatIndex]);
	}
    }
  
  /* Clean up */
  
  LALZDestroyVector(stat->statusPtr, &tempXhata);
  LALZDestroyVector(stat->statusPtr, &tempXhatb);

  LALDDestroyVector(stat->statusPtr, &fA);
  LALDDestroyVector(stat->statusPtr, &fB);
  LALDDestroyVector(stat->statusPtr, &fAB);

  LALFree(tempInt1);
  LALFree(xSum);
  LALFree(ySum);
  
  LALFree(sinVal);
  LALFree(cosVal);

 DETATCHSTATUSPTR(stat);
 RETURN(stat);
 
}


