/************************************ <lalVerbatim file="LALDemodCV">
Author: Berukoff, S.J.,  Papa, M.A., $Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>
\subsection{Module \texttt{LALDemod.c}}\label{ss:LALDemod.c}
Computes a demodulated Fourier transform (DeFT) given a set of input short Fourier transforms (SFT).

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALDemodD}
\idx{LALDemod()}

\subsubsection*{Description}
As described in sections \ref{demodulation_intro} and \ref{s:LALDemod.h}, the
demodulation routine takes a set of SFTs and combines them into a longer time
baseline, higher frequency resolution, transform.  If the constraints
described in the previous sections are met, regarding the time baseline of the
SFTs with respect to the phase evolution of the SFT signal, then the total
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
one currently under production.

The {\bf input} is: \verb@**input@, an array of structures of type \verb@FFT@,
containing not only the data itself but quantities relevant to the search.  It
is composed of a LAL data structure \verb@COMPLEX16FrequencySeries *@ and a
\verb@ParameterSet@ structure; the latter allows one to store specific
information about an FFT.  Note further that the
\verb@COMPLEX16FrequencySeries@ structure particularly contains
\begin{description}
\item[\texttt{input[k]->fft->data->data[i]}]	${\rm{i^{th}}}$ complex data point of the ${\rm{k^{th}}}$ SFT.
\item[\texttt{input[k]->fft->data->length}]		Length of data set.
\end{description}
These are perhaps the most oft-referenced parts of these structures.  

The {\bf output} is a pointer to the first datum of the DeFT, '\verb@xHat@',
that has been produced.  This too is a structure of type \verb@FFT **@, like
the data structure which holds the SFT data.  In general \verb@xHat@ may hold
more than a single DeFT.  By knowing  \verb@mObsCoh@ and \verb@iCoh@ the
demodulation routine will fill up the appropriate part of \verb@xHat@. 

\subsubsection*{Algorithm}

The routine implements the analytical result of eq. \ref{DeFT_algo}.  It thus
uses a nested-loop structure, which computes the frequency components of the
DeFT with each loop, placing them in an array (\verb@xHat@) for future use.
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

#include <lal/LALDemod.h>

NRCSID(LALDEMODC, "$Id$");

/* <lalVerbatim file="LALDemodD"> */
void LALDemod(LALStatus *stat, 
	FFT **xHat, 
	FFT **input, 
	DemodPar *params) 
{ 
/* </lalVerbatim> */

INT4	l,beta,alpha;	/* loop indices */
INT4	xHatIndex=-1;  /* indices for DeFFT vector */
INT4	alpha1, alpha2;	/* indices labelling SFT being worked on */
REAL8	xSum, ySum;	/* temp variables for computation of fs*as and fs*bs */
INT4	s=0;		/* local variable for spinDwn calcs. */
REAL8	xTemp;	/* temp variable for phase model */
REAL8 	deltaF;	/* width of SFT band */
INT4	k, k1, k2;	/* defining the sum over which input->fft->data->data*P is calculated */
REAL8 	*skyConst;	/* vector of sky constants data */
REAL8 	*spinDwn;	/* vector of spinDwn parameters (maybe a structure? */
INT4	spOrder;	/* maximum spinDwn order */
REAL8	x;		/* local variable for holding x */
REAL8	realTemp, imagTemp;	/* temp variables used in computation of input->fft->data->data*P */
REAL8	realP, imagP;	/* real and imaginary parts of P, see CVS */
INT4	nDeltaF;	/* number of frequency bins per SFT band */
INT4	sftIndex;	/* more temp variables */
REAL8	y;		/* local variable for holding y */
REAL8	realQ, imagQ;  /* real and imaginary parts of Q, see CVS */
INT4 	iCoh;

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

/* CHECK ALLOCATION, SO WE HAVE SOMETHING TO WRITE TO */
ASSERT(xHat[iCoh]!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);
ASSERT(xHat[iCoh]->fft!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);
ASSERT(xHat[iCoh]->fft->data!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);
ASSERT(xHat[iCoh]->fft->data->data!=NULL, stat, LALDEMODH_ENULL, LALDEMODH_MSGENULL);

/* setting some values */

alpha1=iCoh*params->mCohSFT;
alpha2=alpha1+params->mCohSFT;
/* xHatIndex=params->xHatIndex; */

/* the loops begin */
for(l=params->if0Min;l<params->if0Max;l++)
{
  for(beta=0;beta<params->mCohSFT;beta++)
  {
		xHatIndex++;
		xHat[iCoh]->fft->data->data[xHatIndex].re=0.0;	
		xHat[iCoh]->fft->data->data[xHatIndex].im=0.0;
		
    for(alpha=alpha1;alpha<alpha2;alpha++)
    {
      xSum=0.0;
      ySum=0.0;
      
	/* this loop computes the values of the phase model */
	for(s=0; s<spOrder;s++) {
		xSum += *(spinDwn+s) * (skyConst[2*alpha*(spOrder+1)+2*s+3]); 	
		ySum += *(spinDwn+s) * (skyConst[2*alpha*(spOrder+1)+2*s+2]);
		}
      xTemp=((REAL8)beta/(REAL8)params->mCohSFT+(REAL8)l)*skyConst[2*alpha*(spOrder+1)+1]*deltaF+xSum;
      realTemp=0.0;
      imagTemp=0.0;
      k1=(INT4)xTemp-NTERM_COH_DIV_TWO+1;
      k2=(INT4)xTemp+NTERM_COH_DIV_TWO;
      for(k=k1;k<=k2;k++)
      {
		x=LAL_TWOPI*((REAL8)xTemp-(REAL8)k);
		realP=1.0;
		imagP=0.0;
		if(fabs(x)>SMALL)
		{
	  		realP=sin(x)/x;
	  		imagP=(cos(x)-1.0)/x;
		}
	
		sftIndex=k-params->ifMin;
		
		/* these four lines compute P*xtilde */
		realTemp+= input[alpha]->fft->data->data[sftIndex].re*realP;
		realTemp-= input[alpha]->fft->data->data[sftIndex].im*imagP;
		imagTemp+= input[alpha]->fft->data->data[sftIndex].re*imagP;
		imagTemp+= input[alpha]->fft->data->data[sftIndex].im*realP;
      }

y=-LAL_TWOPI*(((REAL8)beta/(REAL8)params->mCohSFT+(REAL8)l) * skyConst[2*alpha*(spOrder+1)]*deltaF+ySum);
      realQ=cos(y);
      imagQ=sin(y);
	
	/* these four lines compute the DeFT */
      xHat[iCoh]->fft->data->data[xHatIndex].re += realTemp*realQ;
      xHat[iCoh]->fft->data->data[xHatIndex].re -= imagTemp*imagQ;
      xHat[iCoh]->fft->data->data[xHatIndex].im += realTemp*imagQ;
	xHat[iCoh]->fft->data->data[xHatIndex].im += imagTemp*realQ;
	
    }
	
	/* Normalization */
	xHat[iCoh]->fft->data->data[xHatIndex].re /= sqrt(params->mCohSFT);
	xHat[iCoh]->fft->data->data[xHatIndex].im /= sqrt(params->mCohSFT);
  }
}


DETATCHSTATUSPTR(stat);
RETURN(stat);

}
