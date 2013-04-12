/*
*  Copyright (C) 2007 Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens, Bruce Allen
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStatusMacros.h>
#include <lal/LALDemod.h>


/**
\author Berukoff, S.J.,  Papa, M.A., Allen, B., Siemens, X.

\brief Computes a demodulated Fourier transform (DeFT) given a set of input short Fourier transforms (SFT).

\heading{Description}

This routine computes the \f$F\f$ statistic for a set of templates that are defined by: one sky
position, a set of spin-down parameters and a band of possible signal frequencies. The \f$F\f$
statistic is described in JKS, Phys Rev D 58, 063001 (1998). Here, it has been adapted to a single
emission frequency model.

The <tt>parameter structure</tt> defines the search frequency band (\c f0 and \c imax), the search
frequency resolution (\c df) the first frequency of the input SFTs (\c ifmin), how many SFTs have to
be combined (\c SFTno) and template parameters (<tt>*spinDwnOrder</tt>, <tt>*spinDwn</tt> and
<tt>*skyConst</tt>). \c amcoe contains the values of the amplitude modulation functions \f$a\f$ and
\f$b\f$.  \c Dterms represents the numbers of terms to be summed to compute the Dirichlet kernel on
each side of the instantaneous frequency.

The \c input is: <tt>**input</tt>, an array of structures of type \c FFT. This data type will soon
disappear as it is just a COMPLEX8FrequencySeries.

The \c output is a pointer a structure of type \c LALFstat containing an array of the values of
\f$\mathcal{F}\f$. In addition, if <tt>DemodPar-\>returnFaFb == TRUE</tt>, the values of \f$F_a\f$
and \f$F_b\f$ will be returned in addition. (Memory has to be allocated correctly beforehand!)

\heading{Algorithm}

The routine implements the analytical result of Eq.\eqref{DeFT_algo}.  It thus
uses a nested-loop structure, which computes \f$F\f$ for all the template frequencies.

The outer most loop is over the search frequencies.   The next loop is over
\f$\alpha\f$, which identifies the SFTs.
The value of  \f$k^*\f$ is then computed using the second of Eq.\eqref{DeFT_defs},
and thus the summation over \f$k\f$ of Eq.\eqref{DeFT_algo} is carried out, with a
loop over Dterms.  In this loop the product \f$\tilde{x}_{\alpha k}P_{\alpha k}\f$
is calculated.  Once this loop completes,
\f$e^{iy_\alpha}\f$ is computed, the summation over \f$\alpha\f$ performed and,
finally,  the code yields the DeFT \f$\hat{x}_b\f$.  It can be seen that the code
closely follows the analytical development of the formalism.

Finally, note that in order to avoid repeated
trigonometric function computations, a look-up-table (LUT) for sine and cosine is
constructed at the beginning of the routine.
*/
void LALDemod(LALStatus *status, LALFstat *Fstat, FFT **input, DemodPar *params)

{

  INT4 alpha,i;                 /* loop indices */
  REAL8	*xSum=NULL, *ySum=NULL;	/* temp variables for computation of fs*as and fs*bs */
  INT4 s;		        /* local variable for spinDwn calcs. */
  REAL8	xTemp;	                /* temp variable for phase model */
  INT4	k, k1;	                /* defining the sum over which is calculated */
  REAL8 *skyConst;	        /* vector of sky constants data */
  REAL8 *spinDwn;	        /* vector of spinDwn parameters (maybe a structure? */
  INT4	spOrder;	        /* maximum spinDwn order */
  REAL8	x;		        /* local variable for holding x */
  REAL8	realXP, imagXP; 	/* temp variables used in computation of */
  REAL8	realP, imagP;	        /* real and imaginary parts of P, see CVS */
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

  INITSTATUS(status);

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
	    realXP += crealf(Xalpha_k)*realP;
	    realXP -= Xalpha_k.im*imagP;
	    imagXP += crealf(Xalpha_k)*imagP;
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
