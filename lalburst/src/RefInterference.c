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

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/CLR.h>


/**
\author Sintes, A. M.

\brief Generates a reference interference signal.

\heading{Description}
Given the complex vector  <tt>*in1</tt> of length \f$n/2+1\f$, containing
the Fourier transform of the data \f$\tilde x(\nu)\f$,
<dl>
<dt><tt>in1->length</tt></dt><dd> The number of elements in
            <tt>in1->data</tt> \f$=n/2+1\f$.</dd>
<dt><tt>in1->data</tt></dt><dd>   The data \f$\tilde x(\nu)\f$,</dd>
</dl>
 and given another vector <tt>*par</tt> containing the  information related
to the harmonics from which we want to construct the reference signal
(i.e., indices, and  initial and final frequency bin locations),
<dl>
<dt><tt>par->length</tt></dt><dd> The number of elements in <tt>par->data</tt>.
       This is equal to three times the number of harmonics that will be
        used to construct the reference signal \f$M(t)\f$.</dd>
<dt><tt>par->data</tt></dt><dd>    \f$\{ k,\nu_{ik}, \nu_{fk} \} \f$,
      e.g.,  \f$\{3, 9868, 9894, 5, 16449, 16487, 9, 29607, 29675\ldots\}\f$,</dd>
</dl>
it  generates the time domain
\f$M(t)\f$ reference interference signal, <tt>*out</tt>. This is a complex
vector of length \f$n\f$.
<dl>
<dt><tt>out->length</tt></dt><dd> The number of elements in
            <tt>out->data</tt> \f$=n\f$.</dd>
<dt><tt>out->data</tt></dt><dd>   \f$M(t)\f$ complex data.</dd>
</dl>
\f$M(t)\f$ corresponds to a nearly monochromatic function near the
frequency  \f$f_0\f$, that is implicit in the information
given in <tt>*par</tt>.

\heading{Algorithm}
Described before.

\heading{Notes}

The harmonics selected to construct the reference signal
should not be (if possible) buried with other strong lines,
 such as violin modes. Choose the strongest harmonics, those
that clearly stand over the noise level.

*/
void LALRefInterference (LALStatus    *status,/**< LAL status pointer */
		   COMPLEX8Vector     *out,   /**<  M(t), size n */
		   COMPLEX8Vector     *in1,   /**<  x(f), size n/2+1 */
                   INT4Vector         *par)   /**< UNDOCUMENTED */
{

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

  INITSTATUS(status);
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
    px += crealf(x[binini-i])*crealf(x[binini-i]) + cimagf(x[binini-i])*cimagf(x[binini-i]);
  }
  for (i=1; i < 11; ++i) {
    px += crealf(x[binfin+i])*crealf(x[binfin+i]) + cimagf(x[binfin+i])*cimagf(x[binfin+i]);
  }

  px = px *(binfin-binini);        /* proportional to the width */

  /* Build z_k(nu) */
  for (i=0; i< binini; ++i) {
    zf->data[i].realf_FIXME = 0.0;
    zf->data[i].imagf_FIXME = 0.0;
  }
  for (i=binini; i< binfin+1; ++i) {
    zf->data[i].realf_FIXME = crealf(x[i]);
    zf->data[i].imagf_FIXME = cimagf(x[i]);
  }
  for (i=binfin+1; i< n; ++i) {
    zf->data[i].realf_FIXME = 0.0;
    zf->data[i].imagf_FIXME = 0.0;
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

  dr = crealf(zt->data[0]) * invN;
  di = cimagf(zt->data[0]) * invN;

  mod2 = dr*dr+di*di;

  if( mod2< LAL_REAL4_MIN)
    {  phaseI= 0.0; /* to avoid NaN */ }
  else
    {  phaseI = atan2(di,dr);  }

  /* calculation of B_k(1)(t) */
  ampB = pow(mod2,inv2k);
  phB  = phaseI*invk;

  b1t->data[0].real_FIXME = ampB * cos(phB);
  b1t->data[0].im = ampB * sin(phB);

  invarb->data[0] = px*mod2;
  sden->data[0]   = invarb->data[0];

  snum->data[0].realf_FIXME = creal(b1t->data[0]) * invarb->data[0];
  snum->data[0].imagf_FIXME = b1t->data[0].im * invarb->data[0];


  for (i=1; i< n; ++i) {
    dr = crealf(zt->data[i]) * invN;
    di = cimagf(zt->data[i]) * invN;
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

    b1t->data[i].real_FIXME = ampB * cos(phB);
    b1t->data[i].im = ampB * sin(phB);

    invarb->data[i] = px*mod2;
    sden->data[i] = invarb->data[i];

    snum->data[i].realf_FIXME = creal(b1t->data[i]) * invarb->data[i];
    snum->data[i].imagf_FIXME = b1t->data[i].im * invarb->data[i];
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
	px += crealf(x[binini-i])*crealf(x[binini-i]) + cimagf(x[binini-i])*cimagf(x[binini-i]);
      }
    for (i=1; i < 11; ++i)
      {
	px += crealf(x[binfin+i])*crealf(x[binfin+i]) + cimagf(x[binfin+i])*cimagf(x[binfin+i]);
      }

    px = px *(binfin-binini);        /* proportional to the width */

    /* Build z_k(nu) */
    for (i=0; i< binini; ++i) {
      zf->data[i].realf_FIXME = 0.0;
      zf->data[i].imagf_FIXME = 0.0;
    }
    for (i=binini; i< binfin+1; ++i) {
      zf->data[i].realf_FIXME = crealf(x[i]);
      zf->data[i].imagf_FIXME = cimagf(x[i]);
    }
    for (i=binfin+1; i< n; ++i) {
      zf->data[i].realf_FIXME = 0.0;
      zf->data[i].imagf_FIXME = 0.0;
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
    dr = crealf(zt->data[0]) * invN;
    di = cimagf(zt->data[0]) * invN;

    mod2 = dr*dr+di*di;

    if( mod2< LAL_REAL4_MIN)
      {  phaseI= 0.0; /* to avoid NaN */ }
    else
      {  phaseI = atan2(di,dr);  }

    /* calculation of B_k(t) */
    ampB = pow(mod2,inv2k);
    phB = phaseI*invk;

    bt->data[0].real_FIXME = ampB * cos(phB);
    bt->data[0].im = ampB * sin(phB);

    /* initialize calculation of Lambda_k */
    br =  creal(bt->data[0]);
    bi =  bt->data[0].im;
    dr =  creal(b1t->data[0]);
    di =  b1t->data[0].im;

    lambdan.real_FIXME =  br*dr + bi*di;
    lambdan.im = -dr*bi + br*di;
    lambdad    =  br*br + bi*bi;

    /* inverse of the variance of beta */
    invarb->data[0] = px*mod2;

    for (i=1; i< n; ++i) {
      dr = crealf(zt->data[i]) * invN;
      di = cimagf(zt->data[i]) * invN;
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

      bt->data[i].real_FIXME = ampB * cos(phB);
      bt->data[i].im = ampB * sin(phB);

      /* for the calculation of Lambda_k */
      br =  creal(bt->data[i]);
      bi =  bt->data[i].im;
      dr =  creal(b1t->data[i]);
      di =  b1t->data[i].im;

      lambdan.real_FIXME +=  br*dr + bi*di;
      lambdan.im += -dr*bi + br*di;
      lambdad    +=  br*br + bi*bi;

      /* inverse of the variance of beta */
      invarb->data[i] = px*mod2;
    }

    /* update  sden and snum */
    lambda.real_FIXME = creal(lambdan)/lambdad;
    lambda.im = lambdan.im/lambdad;

    for(i=0; i< n; ++i) {
      br =  creal(bt->data[i]);
      bi =  bt->data[i].im;
      snum->data[i].realf_FIXME += (br*creal(lambda) - bi*lambda.im) * invarb->data[i];
      snum->data[i].imagf_FIXME += (br*lambda.im + bi*creal(lambda)) * invarb->data[i];
      sden->data[i]    += invarb->data[i];
    }

  } /* the harmonics */

  /* -------------------------------------------   */
  /*   calculation of M(t)                         */
  /* -------------------------------------------   */

  for(i=0; i< n; ++i){
    m[i].realf_FIXME = crealf(snum->data[i]) / sden->data[i] ;
    m[i].imagf_FIXME = cimagf(snum->data[i]) / sden->data[i] ;
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
