/************************************ <lalVerbatim file="LALDemodCV">
Author: Berukoff, S.J.,  Papa, M.A., Allen, B., Siemens, X. $Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>
\subsection{Module \texttt{LALDemod.c}}\label{ss:LALDemod.c}
Computes a demodulated Fourier transform (DeFT) given a set of input short Fourier transforms (SFT).

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{LALDemodCP}
\idx{LALDemod()}

\subsubsection*{Description}

This routine computes the $F$ statistic for a set of templates that are defined by: one sky position, a set of spin-down parameters and a band of possible signal frequencies. The $F$ statistic is described in JKS, Phys Rev D 58, 063001 (1998). Here, it has been adapted to a single emission frequency model.

The {\bf parameter structure} defines the search frequency band 
(\verb@f0@ and \verb@imax@), the search frequency resolution (\verb@df@) the first frequency of the input SFTs (\verb@ifmin@), how many SFTs have to be combined (\verb@SFTno@) and template parameters (\verb@*spinDwnOrder@, \verb@*spinDwn@ and @\verb@*skyConst@). \verb@amcoe@ contains the values of the amplitude modulation functions $a$ and $b$.  \verb@Dterms@ represents the numbers of terms to be summed to compute the Dirichlet kernel on each side of the instantaneous frequency.

The {\bf input} is: \verb@**input@, an array of structures of type \verb@FFT@. This data type will soon disappear as it is just a complex8frequencyseries.

The {\bf output} is a pointer a structure of type \verb+LALFstat+ containing an array of the values of $\mathcal{F}$. 
In addition, if \verb+DemodPar->returnFaFb == TRUE+, the values of $F_a$ and $F_b$ will be returned in addition.
(Memory has to be allocated correctly beforehand!)

\subsubsection*{Algorithm}

The routine implements the analytical result of eq. \ref{DeFT_algo}.  It thus
uses a nested-loop structure, which computes $F$ for all the template frequencies.

The outer most loop is over the search frequencies.   The next loop is over
$\alpha$, which identifies the SFTs.
The value of  $k^*$ is then computed using the second of Eq. \ref{DeFT_defs},
and thus the summation over $k$ of Eq.\ref{DeFT_algo} is carried out, with a
loop over Dterms.  In this loop the product $\tilde{x}_{\alpha k}P_{\alpha k}$
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

#include <lal/LALStatusMacros.h>
#include <lal/LALDemod.h>
NRCSID( LALDEMODC, "$Id$" );

/* <lalVerbatim file="LALDemodCP"> */
void LALDemod(LALStatus *status, LALFstat *Fstat, FFT **input, DemodPar *params) 
/* </lalVerbatim> */
{ 

  INT4 alpha,i;                 /* loop indices */
  REAL8	*xSum=NULL, *ySum=NULL;	/* temp variables for computation of fs*as and fs*bs */
  INT4 s;		        /* local variable for spinDwn calcs. */
  REAL8	xTemp;	                /* temp variable for phase model */
  REAL8	deltaF;	                /* width of SFT band */
  INT4	k, k1;	                /* defining the sum over which is calculated */
  REAL8 *skyConst;	        /* vector of sky constants data */
  REAL8 *spinDwn;	        /* vector of spinDwn parameters (maybe a structure? */
  INT4	spOrder;	        /* maximum spinDwn order */
  REAL8	x;		        /* local variable for holding x */
  REAL8	realXP, imagXP; 	/* temp variables used in computation of */
  REAL8	realP, imagP;	        /* real and imaginary parts of P, see CVS */
  INT4	nDeltaF;	        /* number of frequency bins per SFT band */
  INT4	sftIndex;	        /* more temp variables */
  REAL8	y;		        /* local variable for holding y */
  REAL8 realQ, imagQ;
  REAL8 *sinVal,*cosVal;        /*LUT values computed by the routine do_trig_lut*/
  INT4  res;                    /*resolution of the argument of the trig functions: 2pi/res.*/
  INT4 *tempInt1;
  INT4 myindex;
  REAL8 FaSq;
  REAL8 FbSq;
  REAL8 FaFb;
  COMPLEX16 Fa, Fb;

  REAL8 f;

  REAL8 A=params->amcoe->A,B=params->amcoe->B,C=params->amcoe->C,D=params->amcoe->D;
  INT4 M=params->SFTno;

  INITSTATUS( status, "LALDemod", LALDEMODC );

  /* catch some obvious programming errors */
  ASSERT ( (Fstat != NULL)&&(Fstat->F != NULL), status, LALDEMODH_ENULL, LALDEMODH_MSGENULL );
  if (params->returnFaFb)
    {
      ASSERT ( (Fstat->Fa != NULL)&&(Fstat->Fb != NULL), status, LALDEMODH_ENULL, LALDEMODH_MSGENULL );
    }

  /* variable redefinitions for code readability */
  spOrder=params->spinDwnOrder;
  spinDwn=params->spinDwn;
  skyConst=params->skyConst;
  deltaF=(*input)->fft->deltaF;
  nDeltaF=(*input)->fft->data->length;

  /* res=10*(params->mCohSFT); */
  /* This size LUT gives errors ~ 10^-7 with a three-term Taylor series */
   res=64;
   sinVal=(REAL8 *)LALMalloc((res+1)*sizeof(REAL8));
   cosVal=(REAL8 *)LALMalloc((res+1)*sizeof(REAL8)); 
   for (k=0; k<=res; k++){
     sinVal[k]=sin((LAL_TWOPI*k)/res);
     cosVal[k]=cos((LAL_TWOPI*k)/res);
   }

  /* this loop computes the values of the phase model */
  xSum=(REAL8 *)LALMalloc(params->SFTno*sizeof(REAL8));
  ySum=(REAL8 *)LALMalloc(params->SFTno*sizeof(REAL8));
  tempInt1=(INT4 *)LALMalloc(params->SFTno*sizeof(INT4));
  for(alpha=0;alpha<params->SFTno;alpha++){
    tempInt1[alpha]=2*alpha*(spOrder+1)+1;
    xSum[alpha]=0.0;
    ySum[alpha]=0.0;
    for(s=0; s<spOrder;s++) {
      xSum[alpha] += spinDwn[s] * skyConst[tempInt1[alpha]+2+2*s]; 	
      ySum[alpha] += spinDwn[s] * skyConst[tempInt1[alpha]+1+2*s];
    }
  }


  /* Loop over frequencies to be demodulated */
  for(i=0 ; i< params->imax  ; i++ )
  {
    Fa.re =0.0;
    Fa.im =0.0;
    Fb.re =0.0;
    Fb.im =0.0;

    f=params->f0+i*params->df;

    /* Loop over SFTs that contribute to F-stat for a given frequency */
    for(alpha=0;alpha<params->SFTno;alpha++)
      {
	REAL8 tsin, tcos, tempFreq;
	COMPLEX8 *Xalpha=input[alpha]->fft->data->data;
	REAL4 a = params->amcoe->a->data[alpha];
	REAL4 b = params->amcoe->b->data[alpha];

	xTemp=f*skyConst[tempInt1[alpha]]+xSum[alpha];

	realXP=0.0;
	imagXP=0.0;
	      
	/* find correct index into LUT -- pick closest point */
	tempFreq=xTemp-(INT4)xTemp;
	myindex=(INT4)(tempFreq*res+0.5);
	      
	{
	  REAL8 d=LAL_TWOPI*(tempFreq-(REAL8)myindex/(REAL8)res);
	  REAL8 d2=0.5*d*d;
	  REAL8 ts=sinVal[myindex];
	  REAL8 tc=cosVal[myindex];
		
	  tsin=ts+d*tc-d2*ts;
	  tcos=tc-d*ts-d2*tc-1.0;
	}
		     
	tempFreq=LAL_TWOPI*(tempFreq+params->Dterms-1);
	k1=(INT4)xTemp-params->Dterms+1;
	/* Loop over terms in dirichlet Kernel */
	for(k=0;k<2*params->Dterms;k++)
	  {
	    COMPLEX8 Xalpha_k;
	    x=tempFreq-LAL_TWOPI*(REAL8)k;
	    realP=tsin/x;
	    imagP=tcos/x;

	    /* If x is small we need correct x->0 limit of Dirichlet kernel */
	    if(fabs(x) < SMALL) 
	      {
		realP=1.0;
		imagP=0.0;
	      }	 
 
	    sftIndex=k1+k-params->ifmin;

	    /* these four lines compute P*xtilde */
	    Xalpha_k=Xalpha[sftIndex];
	    realXP += Xalpha_k.re*realP;
	    realXP -= Xalpha_k.im*imagP;
	    imagXP += Xalpha_k.re*imagP;
	    imagXP += Xalpha_k.im*realP;
	  }
      
	y=-LAL_TWOPI*(f*skyConst[tempInt1[alpha]-1]+ySum[alpha]);

	realQ = cos(y);
	imagQ = sin(y);

	/* implementation of amplitude demodulation */
	{
	  REAL8 realQXP = realXP*realQ-imagXP*imagQ;
	  REAL8 imagQXP = realXP*imagQ+imagXP*realQ;
	  Fa.re += a*realQXP;
	  Fa.im += a*imagQXP;
	  Fb.re += b*realQXP;
	  Fb.im += b*imagQXP;
	}
    }      

    FaSq = Fa.re*Fa.re+Fa.im*Fa.im;
    FbSq = Fb.re*Fb.re+Fb.im*Fb.im;
    FaFb = Fa.re*Fb.re+Fa.im*Fb.im;
	  	  	
    Fstat->F[i] = (4.0/(M*D))*(B*FaSq + A*FbSq - 2.0*C*FaFb);
    if (params->returnFaFb)
      {
	Fstat->Fa[i] = Fa;
	Fstat->Fb[i] = Fb;
      }


  }
  /* Clean up */
  LALFree(tempInt1);
  LALFree(xSum);
  LALFree(ySum);
  
  LALFree(sinVal);
  LALFree(cosVal);

  RETURN( status );

} /* LALDemod() */

