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
 * Revision: $Id:
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


/************************************ <lalVerbatim file="CleanAllCV">
Author: Sintes, A. M.
$Id$
************************************* </lalVerbatim> */


/* <lalLaTeX>

\subsection{Module \texttt{CleanAll.c}}
\label{ss:CleanAll.c}
Gets data cleaned from line harmonic interference given  a time domain
reference signal.


\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CleanAllD}
\idx{LALCleanAll()}

\subsubsection*{Description}
This routine cleans data in the time domain from line harmonic interference
(from the first harmonic up to the Nyquist frequency). The inputs are:

\verb@*in1@ the time domain data of type  \verb@REAL4TVectorCLR@,
containing also the interference fundamental frequency $f_0$ and the
sampling spacing. This information is needed in order to obtain
the total  number of harmonics contained in the data.
\begin{description}
\item[\texttt{in1->length}] The number of elements in \texttt{in1->data} $=n$.
\item[\texttt{in1->data}]   The (real) time domain data,  $x(t)$.
\item[\texttt{in1->deltaT}] The sample spacing in seconds.
\item[\texttt{in1->fLine}]  The interference fundamental frequency $f_0$
       (in Hz), e.g., 60 Hz.
\end{description}

\verb@*in2@ the time domain reference signal (a complex vector).
\begin{description}
\item[\texttt{in2->length}] The number of elements in
            \texttt{in2->data} $=n$.
\item[\texttt{in2->data}]    The $M(t)$ complex data.
\end{description}

The output \verb@*out@ is a real vector containing the clean data.
\begin{description}
\item[\texttt{out->length}] The number of elements in
            \texttt{out->data} $=n$.
\item[\texttt{out->data}]    The clean (real) time domain data.
\end{description}

\subsubsection*{Algorithm}
It takes the reference signal $M(t)$ and, for all possible harmonics
$j$
($j=1,\ldots,$\texttt{floor(1.0/fabs( 2.02* in1->deltaT * in1->fLine))} ),
from the fundamental frequency up to the Nyquist frequency,
constructs $M(t)^j$,  performs a least-squares fit, i.e.,
minimizes the power $\vert x(t) -\rho_j M(t)^j\vert^2$ with
respect to $\rho_j$, and  subtracts $\rho_j M(t)^j$ from the
original data, $x(t)$.

\subsubsection*{Uses}
\begin{verbatim}
LALDCreateVector()
LALZCreateVector()
LALDDestroyVector()
LALZDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{CleanAllCV}}

</lalLaTeX> */



#include <lal/CLR.h>

NRCSID (CLEANALLC, "$Id$");

/* <lalVerbatim file="CleanAllD"> */
void LALCleanAll (LALStatus     *status,
               REAL4Vector      *out,  /* clean data */
               COMPLEX8Vector   *in2,  /* M(t), ref. interference */
               REAL4TVectorCLR  *in1)  /* x(t), data + information */
{ /* </lalVerbatim> */

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

  INITSTATUS (status, "LALCleanAll", CLEANALLC);
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

  rhon.re = 0.0;
  rhon.im = 0.0;
  rhod    = 0.0;

  for (i=0; i<n; ++i) {
    mr = m[i].re;
    mi = m[i].im;

    rhon.re +=  x[i]*mr;
    rhon.im += -x[i]*mi;
    rhod    +=  mr*mr + mi*mi;
  }

  rho.re = rhon.re / rhod;
  rho.im = rhon.im / rhod;

  for (i=0; i<n; ++i) {
    xc[i] = x[i] - 2.0*(rho.re * m[i].re - rho.im * m[i].im );
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
      mr = m[i].re;
      mi = m[i].im;

      ampM2->data[i] = mr*mr+ mi*mi;

      if( ampM2->data[i] < LAL_REAL4_MIN)
        {  phaM->data[i] = 0.0; /* to avoid NaN */ }
      else
        {  phaM->data[i] = atan2(mi,mr);  }
    }

    for(j=2; j<= maxH; ++j) {
      rhon.re = 0.0;
      rhon.im = 0.0;
      rhod    = 0.0;

      for (i=0; i<n; ++i) {
	/* calculation of m^j(t) for a fixed t */
	amj = pow( ampM2->data[i], 0.5*j);
	phj = j * phaM->data[i];
	mj->data[i].re = amj * cos(phj);
	mr =  mj->data[i].re;
	mj->data[i].im = amj * sin(phj);
	mi =  mj->data[i].im;

	rhon.re +=  xc[i]*mr;
	rhon.im += -xc[i]*mi;
	rhod    +=  mr*mr + mi*mi;
      }

      rho.re = rhon.re / rhod;
      rho.im = rhon.im / rhod;

     for (i=0; i<n; ++i) {
       xc[i] += - 2.0*(rho.re * mj->data[i].re - rho.im * mj->data[i].im );
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

