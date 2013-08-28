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
 * File Name: CleanAll.c
 *
 * Author: Sintes, A. M.
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 *  CleanAll.c
 *
 * SYNOPSIS
 *
 *
 * DESCRIPTION
 *   Given the data and a reference signal, removes all the harmonics
 *
 * DIAGNOSTICS
 *
 * CALLS
 *
 * NOTES
 *
 *
 *-----------------------------------------------------------------------
 */

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/CLR.h>

/**
 * \author Sintes, A. M.
 *
 * \brief Gets data cleaned from line harmonic interference given  a time domain reference signal.
 *
 * \heading{Description}
 * This routine cleans data in the time domain from line harmonic interference
 * (from the first harmonic up to the Nyquist frequency). The inputs are:
 *
 * <tt>*in1</tt> the time domain data of type  \c REAL4TVectorCLR,
 * containing also the interference fundamental frequency \f$f_0\f$ and the
 * sampling spacing. This information is needed in order to obtain
 * the total  number of harmonics contained in the data.
 * <dl>
 * <dt><tt>in1->length</tt></dt><dd> The number of elements in <tt>in1->data</tt> \f$=n\f$.</dd>
 * <dt><tt>in1->data</tt></dt><dd>   The (real) time domain data,  \f$x(t)\f$.</dd>
 * <dt><tt>in1->deltaT</tt></dt><dd> The sample spacing in seconds.</dd>
 * <dt><tt>in1->fLine</tt></dt><dd>  The interference fundamental frequency \f$f_0\f$
 * (in Hz), e.g., 60 Hz.</dd>
 * </dl>
 *
 * <tt>*in2</tt> the time domain reference signal (a complex vector).
 * <dl>
 * <dt><tt>in2->length</tt></dt><dd> The number of elements in
 * <tt>in2->data</tt> \f$=n\f$.</dd>
 * <dt><tt>in2->data</tt></dt><dd>    The \f$M(t)\f$ complex data.</dd>
 * </dl>
 *
 * The output <tt>*out</tt> is a real vector containing the clean data.
 * <dl>
 * <dt><tt>out->length</tt></dt><dd> The number of elements in
 * <tt>out->data</tt> \f$=n\f$.</dd>
 * <dt><tt>out->data</tt></dt><dd>    The clean (real) time domain data.</dd>
 * </dl>
 *
 * \heading{Algorithm}
 * It takes the reference signal \f$M(t)\f$ and, for all possible harmonics
 * \f$j\f$
 * (\f$j=1,\ldots,\f$<tt>floor(1.0/fabs( 2.02* in1->deltaT * in1->fLine))</tt> ),
 * from the fundamental frequency up to the Nyquist frequency,
 * constructs \f$M(t)^j\f$,  performs a least-squares fit, i.e.,
 * minimizes the power \f$\vert x(t) -\rho_j M(t)^j\vert^2\f$ with
 * respect to \f$\rho_j\f$, and  subtracts \f$\rho_j M(t)^j\f$ from the
 * original data, \f$x(t)\f$.
 */
void LALCleanAll (LALStatus     *status,/**< LAL status pointer */
               REAL4Vector      *out,  /**< clean data */
               COMPLEX8Vector   *in2,  /**< M(t), ref. interference */
               REAL4TVectorCLR  *in1)  /**< x(t), data + information */
{

  INT4    n;
  INT4    i,j;
  INT4    maxH;


  REAL4      *x;
  REAL4      *xc;
  COMPLEX8   *m;

  COMPLEX16  rho,rhon;
  REAL8      rhod;
  REAL8      mr,mi;
  REAL8      amj,phj;

  COMPLEX16Vector *mj    = NULL;
  REAL8Vector     *ampM2 = NULL;
  REAL8Vector     *phaM  = NULL;

/* --------------------------------------------- */

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (out, status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (in1, status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (in2, status, CLRH_ENULL, CLRH_MSGENULL);

  /*   Make sure the data pointers are not NULL: */
  ASSERT (out->data, status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (in1->data, status, CLRH_ENULL, CLRH_MSGENULL);
  ASSERT (in2->data, status, CLRH_ENULL, CLRH_MSGENULL);

  /*   Make sure that the size is correct:  */
  ASSERT (out->length > 0, status, CLRH_ESIZE, CLRH_MSGESIZE);

  /*   Make sure that the lengths are correct (size mismatch): */
  ASSERT (in1->length == out->length, status, CLRH_ESZMM, CLRH_MSGESZMM);
  ASSERT (in2->length == out->length, status, CLRH_ESZMM, CLRH_MSGESZMM);

  /*   Make sure the arguments are not pointing to the same memory location */
  /*   ASSERT (out != in2, status, CLRH_ESAME, CLRH_MSGESAME); */

  /*   Make sure F_Nyquist > fLine  */
  ASSERT (fabs(in1->deltaT * in1->fLine) < 0.495,
                status, CLRH_EFREQ, CLRH_MSGEFREQ);

  /* margin 0.495 and not 0.5 in order to avoid errors */


  /* -------------------------------------------   */

  n  = out->length;
  m  = in2->data;
  x  = in1->data;
  xc = out->data;

  /* -------------------------------------------   */
  /*     Removing the fundamental harmonic         */
  /* -------------------------------------------   */

  rhon.real_FIXME = 0.0;
  rhon.imag_FIXME = 0.0;
  rhod    = 0.0;

  for (i=0; i<n; ++i) {
    mr = crealf(m[i]);
    mi = cimagf(m[i]);

    rhon.real_FIXME +=  x[i]*mr;
    rhon.imag_FIXME += -x[i]*mi;
    rhod    +=  mr*mr + mi*mi;
  }

  rho.real_FIXME = creal(rhon) / rhod;
  rho.imag_FIXME = cimag(rhon) / rhod;

  for (i=0; i<n; ++i) {
    xc[i] = x[i] - 2.0*(creal(rho) * crealf(m[i]) - cimag(rho) * cimagf(m[i]) );
  }


  /* -------------------------------------------   */
  /*     Removing the other harmonics              */
  /* -------------------------------------------   */

  maxH= floor(1.0/fabs( 2.02* in1->deltaT * in1->fLine));

  if (maxH > 1) /* in case there are more harmonics */
  {

    /* Create Vectors amplitude and phase m(t) */
    TRY(LALDCreateVector(status->statusPtr, &ampM2, n), status);
    TRY(LALDCreateVector(status->statusPtr, &phaM,  n), status);
    /* reserve space for m^j(t) */
    TRY(LALZCreateVector(status->statusPtr, &mj  ,n), status);

    for (i=0; i<n; ++i) {
      /* calculation of amplitude^2 and phase of m(t)  */
      mr = crealf(m[i]);
      mi = cimagf(m[i]);

      ampM2->data[i] = mr*mr+ mi*mi;

      if( ampM2->data[i] < LAL_REAL4_MIN)
        {  phaM->data[i] = 0.0; /* to avoid NaN */ }
      else
        {  phaM->data[i] = atan2(mi,mr);  }
    }

    for(j=2; j<= maxH; ++j) {
      rhon.real_FIXME = 0.0;
      rhon.imag_FIXME = 0.0;
      rhod    = 0.0;

      for (i=0; i<n; ++i) {
	/* calculation of m^j(t) for a fixed t */
	amj = pow( ampM2->data[i], 0.5*j);
	phj = j * phaM->data[i];
	mj->data[i].real_FIXME = amj * cos(phj);
	mr =  creal(mj->data[i]);
	mj->data[i].imag_FIXME = amj * sin(phj);
	mi =  cimag(mj->data[i]);

	rhon.real_FIXME +=  xc[i]*mr;
	rhon.imag_FIXME += -xc[i]*mi;
	rhod    +=  mr*mr + mi*mi;
      }

      rho.real_FIXME = creal(rhon) / rhod;
      rho.imag_FIXME = cimag(rhon) / rhod;

     for (i=0; i<n; ++i) {
       xc[i] += - 2.0*(creal(rho) * creal(mj->data[i]) - cimag(rho) * cimag(mj->data[i]) );
     }

    } /* closes for all harmonics */

    /* Destroy Vectors  */
    TRY(LALDDestroyVector(status->statusPtr, &ampM2), status);
    TRY(LALDDestroyVector(status->statusPtr, &phaM), status);
    TRY(LALZDestroyVector(status->statusPtr, &mj), status);

  } /* closes if */

  /* -------------------------------------------   */

  DETATCHSTATUSPTR (status);

  /* normal exit */
  RETURN (status);
}
