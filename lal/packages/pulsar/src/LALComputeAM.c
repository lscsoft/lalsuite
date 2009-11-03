/*
*  Copyright (C) 2007 Jolien Creighton, Maria Alessandra Papa, Reinhard Prix, Steve Berukoff, Xavier Siemens
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

/********************************* <lalVerbatim file="LALComputeAMCV">
Author:  Berukoff, S.J.  $Id$
********************************** </lalVerbatim> */

/* <lalLaTeX>
\subsection{Module \texttt{LALComputeAM.c}}\label{ss:LALComputeAM.c}
Computes quantities for amplitude demodulation.

\subsection*{Prototypes}
\vspace{0.1in}
\input{LALComputeAMD}
\index{\texttt{LALComputeAM()}}

\subsubsection*{Description}

This routine computes the quantities $a(t)$ and $b(t)$ as defined in
Jaranowski, Krolak, and Schutz (gr-qc/9804014), hereafter JKS.  These
functions quantify the dependence of the detector output on the
beam-pattern functions $F_{+}$ and $F_{\times}$; in fact, $a(t)$ and
$b(t)$ {\it{are}} the beam-pattern functions, without the dependence
on polarization angle and detector arm angle.  Since the
\verb@LALDemod()@ suite is an attempt to compute an optimal statistic,
it is necessary to include these quantiti es in the computation.
Otherwise, the motion of the Earth as it revolves about its axis will
smear the signal into several neighboring bins centered about the
search frequency, consequently losing valuable SNR.

\subsubsection*{Algorithm}

The routine is really simple.  From JKS,
\begin{eqnarray}
F_{+} = \sin
\zeta [ a(t) \cos 2 \psi + b(t) \sin 2 \psi ] \\
F_{\times} = \sin
\zeta [ b(t) \cos 2 \psi - a(t) \sin 2 \psi ]
\end{eqnarray}
We use the routine \verb@LALComputeDetAMResponse()@ to calculate
$F_{+}$ and $F_{\times}$ for a given polarization angle, and then
extract $a(t)$ and $b(t)$, once for each timestamp $t$.  Additionally,
computation of the optimal statistic requires that we compute inner
products of these two quantities for later use.  See \ref{Deft algo}
for more on the optimal statistic.

\subsubsection*{Uses}
\begin{verbatim}
LALBarycenter()
LALBarycenterEarth()
LALComputeDetAMResponse()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALComputeAMCV}}

</lalLaTeX> */

/* loop protection */
#ifndef LALCOMPUTEAM_C
#define LALCOMPUTEAM_C
#endif

#include <lal/LALComputeAM.h>

#define SQ(x) (x) * (x)

NRCSID (LALCOMPUTEAMC, "$Id LALComputeAM.c $");

/* <lalVerbatim file="LALComputeAMD"> */
void LALComputeAM (LALStatus          *status,
		   AMCoeffs           *coe,
		   LIGOTimeGPS        *ts,
		   AMCoeffsParams     *params)
{  /* </lalVerbatim> */

  REAL4 zeta;                  /* sine of angle between detector arms        */
  INT4 i;                      /* temporary loop index                       */
  LALDetAMResponse response;   /* output of LALComputeDetAMResponse          */
  LALGPSandAcc timeAndAcc;     /* parameter structure to LALComputeAMResponse*/

  REAL4 sumA2=0.0;
  REAL4 sumB2=0.0;
  REAL4 sumAB=0.0;             /* variables to store scalar products         */
  INT4 length=coe->a->length;  /* length of input time series                */

  REAL4 cos2psi;
  REAL4 sin2psi;               /* temp variables                             */

  INITSTATUS(status, "LALComputeAM", LALCOMPUTEAMC);
  ATTATCHSTATUSPTR(status);

  /* Must put an ASSERT checking that ts vec and coe vec are same length */

  /* Compute the angle between detector arms, then the reciprocal */
  {
    LALFrDetector det = params->das->pDetector->frDetector;
    zeta = 1.0/(sin(det.xArmAzimuthRadians - det.yArmAzimuthRadians));
    if(params->das->pDetector->type == LALDETECTORTYPE_CYLBAR) zeta=1.0;
  }

  cos2psi = cos(2.0 * params->polAngle);
  sin2psi = sin(2.0 * params->polAngle);

  /* Note the length is the same for the b vector */
  for(i=0; i<length; ++i)
    {
      REAL4 *a = coe->a->data;
      REAL4 *b = coe->b->data;

      timeAndAcc.gps=ts[i];
      /* Compute F_plus, F_cross */
      LALComputeDetAMResponse(status->statusPtr, &response, params->das, &timeAndAcc);

      /*  Compute a, b from JKS eq 10,11
       *  a = zeta * (F_plus*cos(2\psi)-F_cross*sin(2\psi))
       *  b = zeta * (F_cross*cos(2\psi)+Fplus*sin(2\psi))
       */
      a[i] = zeta * (response.plus*cos2psi-response.cross*sin2psi);
      b[i] = zeta * (response.cross*cos2psi+response.plus*sin2psi);

      /* Compute scalar products */
      sumA2 += SQ(a[i]);                       /*  A  */
      sumB2 += SQ(b[i]);                       /*  B  */
      sumAB += (a[i]) * (b[i]);                /*  C  */
    }

  {
    /* Normalization factor */
    REAL8 L = 2.0/(REAL8)length;

    /* Assign output values and normalise */
    coe->A = L*sumA2;
    coe->B = L*sumB2;
    coe->C = L*sumAB;
    coe->D = ( coe->A * coe->B - SQ(coe->C) );
    /* protection against case when AB=C^2 */
    if(coe->D == 0) coe->D=1.0e-9;
  }

  /* Normal exit */

  DETATCHSTATUSPTR(status);
  RETURN(status);
}











