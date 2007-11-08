/*
*  Copyright (C) 2007 Yi Pan, Duncan Brown
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

#if 0
<lalVerbatim file="LALInspiralComputePTFMetricCV">
Author: Yi Pan, Duncan Brown
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{LALInspiralComputePTFMetric.c}}

Module to compute the components of the metric which is used to describe
distances on Physical Template Family signal manifold.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{XLALInspiralComputePTFIntrinsicMetricCP}
\idx{XLALInspiralComputePTFIntrinsicMetric()}
\begin{itemize}
   \item \texttt{metric,} Output, the metric at the lattice point defined by \texttt{params}
   \item \texttt{psd,} Input, the power spectral density of the data
   \item \texttt{params,} Input, the parameters where metric must be computed
   in the computation of the metric.
\end{itemize}

\input{XLALInspiralComputePTFFullMetricCP}
\idx{XLALInspiralComputePTFFullMetric()}
\begin{itemize}
   \item \texttt{metric,} Output, the metric at the lattice point defined by \texttt{params}
   \item \texttt{psd,} Input, the power spectral density of the data
   \item \texttt{params,} Input, the parameters where metric must be computed
   in the computation of the metric.
\end{itemize}

\subsubsection*{Description}
We calculate the components of the metric using the procedure outlined 
by Yi.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\begin{verbatim}
LALMalloc
LALFree
\end{verbatim}

\subsubsection*{Notes}
 
\vfill{\footnotesize\input{LALInspiralComputePTFMetricCV}}

</lalLaTeX>
#endif

#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>

/* <lalVerbatim file="XLALInspiralComputePTFIntrinsicMetricCP">  */
INT4 XLALInspiralComputePTFIntrinsticMetric (
    InspiralMetric             *metric,
    REAL8FrequencySeries       *psd,
    InspiralTemplate           *params
    )
/* </lalVerbatim> */
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  static const char* func = "XLALInspiralComputePTFIntrinsticMetric";

  /* number of points in a time-domain segment */
  UINT4 N = 2 * (psd->data->length - 1); 
  UINT4 i, j, k;

  /* some useful numbers */
  REAL4 sqrtoftwo      = sqrt(2.0);
  REAL4 onebysqrtoftwo = 1.0 / sqrtoftwo;
  REAL4 onebysqrtofsix = 1.0 / sqrt(6.0);

  /* get a local copy of the intrinstic parameters */
  REAL8 chirpMass = params->chirpMass;
  REAL8 eta = params->eta;
  REAL8 chi = params->chi;
  REAL8 kappa = params->kappa;

  /* get a local copy of the extrinstic parameters */
  REAL8 Theta = params->sourceTheta;
  REAL8 Phi = params->sourcePhi;
  REAL8 Psi = params->polarisationAngle;
  REAL8 phi0 = params->startPhase;
  REAL8 t0 = params->startTime;

  /* bounds on the power spectrum integration for the moments */
  REAL8 fLower = params->fLower;
  REAL8 fCutoff = params->fCutoff;
  REAL8 deltaT = params->tSampling;

  /* local pointers to the Q and Qtilde data */
  REAL4Vector Q[5];
  COMPLEX8Vector Qtilde[5];

  /* should really check that deltaT from the template agrees with N */
  /* and deltaF from the power spectrum                              */

  /* these variables, and the memory allocation following them should be     */
  /* moved to a differnent function at some point in the future (preferably  */
  /* before this function is ever called multiple times inside a loop)       */
  REAL4VectorSequence          *PTFQ;
  COMPLEX8VectorSequence       *PTFQtilde;
  REAL4Vector                  *PTFphi;
  REAL4Vector                  *PTFomega_2_3;
  REAL4VectorSequence          *PTFe1;
  REAL4VectorSequence          *PTFe2;
  RealFFTPlan                  *fwdPlan;
  
  PTFQ = XLALCreateVectorSequence( 5, N );
  PTFQtilde = XLALCreateCOMPLEX8VectorSequence( 5, N / 2 + 1 );
  PTFphi = XLALCreateVector( N );
  PTFomega_2_3 = XLALCreateVector( N );
  PTFe1 = XLALCreateVectorSequence( 3, N );
  PTFe2 = XLALCreateVectorSequence( 3, N );
  fwdPlan = XLALCreateForwardREAL4FFTPlan( N, 0 );

  /* call the PTF waveform code */
  errcode = XLALFindChirpPTFWaveform( PTFphi, PTFomega_2_3, PTFe1, PTFe2,
      params, deltaT);

  if ( errcode != XLAL_SUCCESS ) XLAL_ERROR( func, errcode );
  
  /* point the dummy variables Q and Qtilde to the actual output structures */
  for ( i = 0; i < 5; ++i )
  {
    Q[i].length      = N;
    Qtilde[i].length = N / 2 + 1;
    Q[i].data        = PTFQ->data + (i * N);
    Qtilde[i].data   = PTFQtilde->data + (i * (N / 2 + 1)) ;
  }  

  /* evaluate the Q^I factors from the dynamical variables */
  for ( j = 0; j < N; ++j )
  {
    REAL4 omega_2_3 = PTFomega_2_3->data[j];
    REAL4 phi       = PTFphi->data[j];
    REAL4 e1x       = PTFe1->data[j];
    REAL4 e1y       = PTFe1->data[N + j];
    REAL4 e1z       = PTFe1->data[2 * N + j]; 
    REAL4 e2x       = PTFe2->data[j];
    REAL4 e2y       = PTFe2->data[N + j];
    REAL4 e2z       = PTFe2->data[2 * N + j];
    
    Q[0].data[j] = omega_2_3 * onebysqrtoftwo * ( cos(2 * phi) * ( e1x * e1x +
          e2y * e2y - e2x * e2x - e1y * e1y ) + 2 * sin(2 * phi) *
        ( e1x * e2x - e1y * e2y ));
    Q[1].data[j] = omega_2_3 * sqrtoftwo * ( cos(2 * phi) * ( e1x * e1y - 
          e2x * e2y ) + sin(2 * phi) * ( e1x * e2y + e1y * e2x ));
    Q[2].data[j] = omega_2_3 * sqrtoftwo * ( cos(2 * phi) * ( e1x * e1z - 
          e2x * e2z ) + sin(2 * phi) * ( e1x * e2z + e1z * e2x ));
    Q[3].data[j] = omega_2_3 * sqrtoftwo * ( cos(2 * phi) * ( e1y * e1z - 
          e2y * e2z ) + sin(2 * phi) * ( e1y * e2z + e1z * e2y ));
    Q[4].data[j] = omega_2_3 * onebysqrtofsix * ( cos(2 * phi) * 
        ( 2 * e2z * e2z - 2 * e1z * e1z + e1x * e1x + e1y * e1y - 
          e2x * e2x - e2y * e2y ) + 2 * sin(2 * phi) * ( e1x * e2x +
            e1y * e2y - 2 * e1z * e2z ));                              
  }
  
  /* Fourier transform the Q's into the Qtilde's */
  for ( i = 0; i < 5; ++i )
  {
    XLALREAL4ForwardFFT( &Qtilde[i], &Q[i], fwdPlan);
  }
  
  /* now compute the metric... */
  for ( i = 0; i < 10; ++i )
  {
    metric->Gamma[i] = (REAL8) i;
    params->Gamma[i] = metric->Gamma[i];
  }

  /* this memory deallocation code should be moved to a separate function */
  XLALDestroyVectorSequence( PTFQ );
  XLALDestroyCOMPLEX8VectorSequence( PTFQtilde );
  XLALDestroyVector( PTFphi );
  XLALDestroyVector( PTFomega_2_3 );
  XLALDestroyVectorSequence( PTFe1 );
  XLALDestroyVectorSequence( PTFe2 );

  /* normal exit */
  return XLAL_SUCCESS;
}

/* <lalVerbatim file="XLALInspiralComputePTFFullMetricCP">  */
INT4 XLALInspiralComputePTFFullMetric (
    InspiralMetric             *metric,
    REAL8FrequencySeries       *psd,
    InspiralTemplate           *params
    )
/* </lalVerbatim> */
{

}

