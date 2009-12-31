/*
*  Copyright (C) 2007 Jolien Creighton
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

/*-----------------------------------------------------------------------
 *
 * File Name: RefInterference.c
 *
 * Author: Sintes, A. M.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 *  RefInterference.c
 *
 * SYNOPSIS
 *
 *
 * DESCRIPTION
 * Generates the reference interference signal.
 *
 * DIAGNOSTICS
 *
 * CALLS
 *
 * NOTES
 *
 *
 *
 *-----------------------------------------------------------------------
 */


/************************************ <lalVerbatim file="RefInterferenceCV">
Author: Sintes, A. M.
$Id$
************************************* </lalVerbatim> */


/* <lalLaTeX>

\subsection{Module \texttt{RefInterference.c}}
\label{ss:RefInterference.c}
Generates a reference interference signal.


\subsubsection*{Prototypes}
\vspace{0.1in}
\input{RefInterferenceD}
\idx{LALRefInterference()}

\subsubsection*{Description}
Given the complex vector  \verb@*in1@ of length $n/2+1$, containing
the Fourier transform of the data $\tilde x(\nu)$,
\begin{description}
\item[\texttt{in1->length}] The number of elements in
            \texttt{in1->data} $=n/2+1$.
\item[\texttt{in1->data}]   The data $\tilde x(\nu)$,
\end{description}
 and given another vector \verb@*par@ containing the  information related
to the harmonics from which we want to construct the reference signal
(i.e., indices, and  initial and final frequency bin locations),
\begin{description}
\item[\texttt{par->length}] The number of elements in \texttt{par->data}.
       This is equal to three times the number of harmonics that will be
        used to construct the reference signal $M(t)$.
\item[\texttt{par->data}]    $\{ k,\nu_{ik}, \nu_{fk} \} $,
      e.g.,  $\{$3, 9868, 9894, 5, 16449, 16487, 9, 29607, 29675$\ldots\}$,
\end{description}
it  generates the time domain
$M(t)$ reference interference signal, \verb@*out@. This is a complex
vector of length $n$.
\begin{description}
\item[\texttt{out->length}] The number of elements in
            \texttt{out->data} $=n$.
\item[\texttt{out->data}]   $M(t)$ complex data.
\end{description}
$M(t)$ corresponds to a nearly monochromatic function near the
frequency  $f_0$, that is implicit in the information
given in \verb@*par@.

\subsubsection*{Algorithm}
Described before.

\subsubsection*{Uses}
\begin{verbatim}
LALSCreateVector()
LALCCreateVector()
LALZCreateVector()
LALCreateReverseComplexFFTPlan()
LALCOMPLEX8VectorFFT()
LALSDestroyVector()
LALCDestroyVector()
LALZDestroyVector()
LALDestroyComplexFFTPlan()
\end{verbatim}

\subsubsection*{Notes}

The harmonics selected to construct the reference signal
should not be (if possible) buried with other strong lines,
 such as violin modes. Choose the strongest harmonics, those
that clearly stand over the noise level.

\vfill{\footnotesize\input{RefInterferenceCV}}

</lalLaTeX> */


#include <lal/CLR.h>


NRCSID (REFINTERFERENCEC, "$Id$");

/* <lalVerbatim file="RefInterferenceD"> */
void LALRefInterference (LALStatus    *status,
		   COMPLEX8Vector     *out,   /*  M(t), size n */
		   COMPLEX8Vector     *in1,   /*  x(f), size n/2+1 */
		   INT4Vector         *par)
{ /* </lalVerbatim> */

  INT4    i,j;
  INT4    n;
  INT4    l;
  INT4    k,binini,binfin;

/* increased precision */

  REAL4      px;
  REAL8      mod2,ampB,phB,inv2k,invk,invN;
  REAL8      /*cumsum,*/phaseI,phaseII,diffph;
  REAL8      dr,di,br,bi;

  COMPLEX16  lambda,lambdan;
  REAL8      lambdad;
/*-----------------------------------*/
  INT8       countPI;

  COMPLEX8   *x;
  COMPLEX8   *m;

  INT4    *harmo;

  COMPLEX8Vector  *zf  = NULL;
  COMPLEX8Vector  *zt  = NULL;
  COMPLEX16Vector  *b1t = NULL;
  COMPLEX16Vector  *bt  = NULL;

  REAL4Vector   *invarb = NULL;
  COMPLEX8Vector  *snum = NULL;
  REAL4Vector     *sden = NULL;

  ComplexFFTPlan  *pinv = NULL;
/* --------------------------------------------- */

  INITSTATUS (status, "LALRefInterference", REFINTERFERENCEC);
  ATTATCHSTATUSPTR (status);


  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (in1 , status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (par, status, CLRH_ENULL, CLRH_MSGENULL);

  /*   Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (in1->data, status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (par->data, status, CLRH_ENULL, CLRH_MSGENULL);

  /*   Make sure that the size is correct:  */
  ASSERT (out->length > 0, status, CLRH_ESIZE, CLRH_MSGESIZE);
  ASSERT (par->length > 0, status, CLRH_ESIZE, CLRH_MSGESIZE);

  /*   Make sure that the parameter length is a multiple of 3:  */
  ASSERT (par->length%3 == 0, status, CLRH_ESIZE, CLRH_MSGESIZE);

  /*   Make sure that the lengths are correct (size mismatch): */
  ASSERT (in1->length == (out->length)/2+1, status, CLRH_ESZMM, CLRH_MSGESZMM);
  /* -------------------------------------------   */

  n = out->length;
  invN = 1.0/n;
  m = out->data;
  x = in1->data;
  l = (par->length)/3;
  harmo = par->data;
  /* -------------------------------------------   */

  /* Create Vectors and fft plan */
  TRY(LALSCreateVector(status->statusPtr, &sden  ,n), status);
  TRY(LALSCreateVector(status->statusPtr, &invarb,n), status);

  TRY(LALCCreateVector(status->statusPtr, &zf, n), status);
  TRY(LALCCreateVector(status->statusPtr, &zt, n), status);

  TRY(LALZCreateVector(status->statusPtr, &bt, n), status);
  TRY(LALZCreateVector(status->statusPtr, &b1t,  n), status);

  TRY(LALCCreateVector(status->statusPtr, &snum, n), status);

  pinv = XLALCreateReverseCOMPLEX8FFTPlan(n, 0);
  if (pinv == NULL)
    ABORTXLAL(status);
  /* -------------------------------------------   */

  /* -------------------------------------------   */
  /* ------ For the 1st harmonic considered ----   */
  /* -------------------------------------------   */

  k = *harmo;
  ++harmo;
  binini =  *harmo;
  ++harmo;
  binfin =  *harmo;

  /*   Make sure that the frequency interval is correct:  */
  ASSERT (binini < binfin, status, CLRH_EINT, CLRH_MSGEINT);
  ASSERT (0< binini, status, CLRH_EINT, CLRH_MSGEINT);
  ASSERT (binfin < (n/2)+1, status, CLRH_EINT, CLRH_MSGEINT);

  /* Calculate px */

  px= 0.0;

  for (i=10; i > 0; --i) {
    px += x[binini-i].re*x[binini-i].re + x[binini-i].im*x[binini-i].im;
  }
  for (i=1; i < 11; ++i) {
    px += x[binfin+i].re*x[binfin+i].re + x[binfin+i].im*x[binfin+i].im;
  }

  px = px *(binfin-binini);        /* proportional to the width */

  /* Build z_k(nu) */
  for (i=0; i< binini; ++i) {
    zf->data[i].re = 0.0;
    zf->data[i].im = 0.0;
  }
  for (i=binini; i< binfin+1; ++i) {
    zf->data[i].re = x[i].re;
    zf->data[i].im = x[i].im;
  }
  for (i=binfin+1; i< n; ++i) {
    zf->data[i].re = 0.0;
    zf->data[i].im = 0.0;
  }

  /* Calculate z_k(t) by performing FFTs */
  if (XLALCOMPLEX8VectorFFT(zt, zf, pinv) != 0)
    ABORTXLAL(status);

  /* calculate invarb, B_k(1)(t) and initializes sden and snum */

  px  = (k*k)/px;
  inv2k= 0.5/k;
  invk= 2.0*inv2k;

  /* cumsum = 0.0; */
  countPI = 0;

  dr = zt->data[0].re * invN;
  di = zt->data[0].im * invN;

  mod2 = dr*dr+di*di;

  if( mod2< LAL_REAL4_MIN)
    {  phaseI= 0.0; /* to avoid NaN */ }
  else
    {  phaseI = atan2(di,dr);  }

  /* calculation of B_k(1)(t) */
  ampB = pow(mod2,inv2k);
  phB  = phaseI*invk;

  b1t->data[0].re = ampB * cos(phB);
  b1t->data[0].im = ampB * sin(phB);

  invarb->data[0] = px*mod2;
  sden->data[0]   = invarb->data[0];

  snum->data[0].re = b1t->data[0].re * invarb->data[0];
  snum->data[0].im = b1t->data[0].im * invarb->data[0];


  for (i=1; i< n; ++i) {
    dr = zt->data[i].re * invN;
    di = zt->data[i].im * invN;
    /* calculation of modulus^2 and phase, and unwrap phase */
    mod2 = dr*dr+di*di;

    if( mod2< LAL_REAL4_MIN)
      {  phaseII= 0.0; /* to avoid NaN */ }
    else
      {  phaseII = atan2(di,dr);  }

    diffph = phaseII - phaseI;
    phaseI = phaseII;
    /* cumsum += LAL_TWOPI*( (diffph < -LAL_PI) - (diffph >  LAL_PI) ); */
    countPI += (diffph < -LAL_PI) - (diffph >  LAL_PI);

    /* calculation of B_k(1)(t) */

    /*  phB = (phaseII + cumsum)*invk; */
    phB = (phaseII + LAL_TWOPI*( countPI%k ) )*invk;

    ampB = pow(mod2,inv2k);

    b1t->data[i].re = ampB * cos(phB);
    b1t->data[i].im = ampB * sin(phB);

    invarb->data[i] = px*mod2;
    sden->data[i] = invarb->data[i];

    snum->data[i].re = b1t->data[i].re * invarb->data[i];
    snum->data[i].im = b1t->data[i].im * invarb->data[i];
  }

  /* ----------------------------------------------   */
  /* ------ For the other harmonics considered ----   */
  /* ----------------------------------------------   */

  for(j=1; j< l; ++j) {
    ++harmo;
    k = *harmo;
    ++harmo;
    binini =  *harmo;
    ++harmo;
    binfin =  *harmo;

    /*   Make sure that the frequency interval is correct:  */
    ASSERT (binini < binfin, status, CLRH_EINT, CLRH_MSGEINT);
    ASSERT (0< binini, status, CLRH_EINT, CLRH_MSGEINT);
    ASSERT (binfin < (n/2)+1, status, CLRH_EINT, CLRH_MSGEINT);

    /* Calculate px */

    px= 0.0;

    for (i=10; i > 0; --i)
      {
	px += x[binini-i].re*x[binini-i].re + x[binini-i].im*x[binini-i].im;
      }
    for (i=1; i < 11; ++i)
      {
	px += x[binfin+i].re*x[binfin+i].re + x[binfin+i].im*x[binfin+i].im;
      }

    px = px *(binfin-binini);        /* proportional to the width */

    /* Build z_k(nu) */
    for (i=0; i< binini; ++i) {
      zf->data[i].re = 0.0;
      zf->data[i].im = 0.0;
    }
    for (i=binini; i< binfin+1; ++i) {
      zf->data[i].re = x[i].re;
      zf->data[i].im = x[i].im;
    }
    for (i=binfin+1; i< n; ++i) {
      zf->data[i].re = 0.0;
      zf->data[i].im = 0.0;
    }

    /* Calculate z_k(t) by performing FFTs */
    if (XLALCOMPLEX8VectorFFT(zt, zf, pinv) != 0)
      ABORTXLAL(status);

    /* calculate invarb, B_k(t)  */

    px  = (k*k)/px;
    inv2k= 0.5/k;
    invk= 2.0*inv2k;

    /* cumsum = 0.0; */
    countPI = 0;
    dr = zt->data[0].re * invN;
    di = zt->data[0].im * invN;

    mod2 = dr*dr+di*di;

    if( mod2< LAL_REAL4_MIN)
      {  phaseI= 0.0; /* to avoid NaN */ }
    else
      {  phaseI = atan2(di,dr);  }

    /* calculation of B_k(t) */
    ampB = pow(mod2,inv2k);
    phB = phaseI*invk;

    bt->data[0].re = ampB * cos(phB);
    bt->data[0].im = ampB * sin(phB);

    /* initialize calculation of Lambda_k */
    br =  bt->data[0].re;
    bi =  bt->data[0].im;
    dr =  b1t->data[0].re;
    di =  b1t->data[0].im;

    lambdan.re =  br*dr + bi*di;
    lambdan.im = -dr*bi + br*di;
    lambdad    =  br*br + bi*bi;

    /* inverse of the variance of beta */
    invarb->data[0] = px*mod2;

    for (i=1; i< n; ++i) {
      dr = zt->data[i].re * invN;
      di = zt->data[i].im * invN;
      /* calculation of modulus^2 and phase, and unwrap phase */
      mod2 = dr*dr+di*di;

      if( mod2< LAL_REAL4_MIN)
	{  phaseII= 0.0; /* to avoid NaN */ }
      else
	{  phaseII = atan2(di,dr);  }

      diffph = phaseII - phaseI;
      phaseI = phaseII;
      /*  cumsum += LAL_TWOPI*( (diffph < -LAL_PI) - (diffph> LAL_PI) ); */
      countPI += (diffph < -LAL_PI) - (diffph> LAL_PI);

      /* calculation of B_k(t) */
      /*  phB = (phaseII + cumsum)*invk; */
      phB = (phaseII + LAL_TWOPI*( countPI%k ) )*invk;

      ampB = pow(mod2,inv2k);

      bt->data[i].re = ampB * cos(phB);
      bt->data[i].im = ampB * sin(phB);

      /* for the calculation of Lambda_k */
      br =  bt->data[i].re;
      bi =  bt->data[i].im;
      dr =  b1t->data[i].re;
      di =  b1t->data[i].im;

      lambdan.re +=  br*dr + bi*di;
      lambdan.im += -dr*bi + br*di;
      lambdad    +=  br*br + bi*bi;

      /* inverse of the variance of beta */
      invarb->data[i] = px*mod2;
    }

    /* update  sden and snum */
    lambda.re = lambdan.re/lambdad;
    lambda.im = lambdan.im/lambdad;

    for(i=0; i< n; ++i) {
      br =  bt->data[i].re;
      bi =  bt->data[i].im;
      snum->data[i].re += (br*lambda.re - bi*lambda.im) * invarb->data[i];
      snum->data[i].im += (br*lambda.im + bi*lambda.re) * invarb->data[i];
      sden->data[i]    += invarb->data[i];
    }

  } /* the harmonics */

  /* -------------------------------------------   */
  /*   calculation of M(t)                         */
  /* -------------------------------------------   */

  for(i=0; i< n; ++i){
    m[i].re = snum->data[i].re / sden->data[i] ;
    m[i].im = snum->data[i].im / sden->data[i] ;
  }
  /* -------------------------------------------   */


  /* Destroy Vectors and fft plan */
  TRY(LALSDestroyVector(status->statusPtr, &sden), status);
  TRY(LALSDestroyVector(status->statusPtr, &invarb), status);

  TRY(LALCDestroyVector(status->statusPtr, &zf), status);
  TRY(LALCDestroyVector(status->statusPtr, &zt), status);
  TRY(LALZDestroyVector(status->statusPtr, &bt), status);
  TRY(LALZDestroyVector(status->statusPtr, &b1t), status);
  TRY(LALCDestroyVector(status->statusPtr, &snum), status);

  XLALDestroyCOMPLEX8FFTPlan(pinv);
  /* -------------------------------------------   */



  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}
