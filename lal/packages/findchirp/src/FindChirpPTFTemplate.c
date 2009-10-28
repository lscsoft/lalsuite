/*
*  Copyright (C) 2007 Diego Fazi, Duncan Brown
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
 * File Name: FindChirpPTFTemplate.c
 *
 * Author: Brown, D. A., and Fazi, D.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpPTFTemplateCV">
Author: Brown, D. A., and Fazi, D.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpPTFTemplate.c}}
\label{ss:FindChirpPTFTemplate.c}

Provides functions to create physical template family templates in a
form that can be used by the \texttt{FindChirpPTFFilter()} function.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpPTFTemplateCP}
\idx{LALFindChirpPTFTemplate()}

The function \texttt{LALFindChirpPTFTemplate()} creates a physical template
family template as described by the algorithm below.

\subsubsection*{Algorithm}

Blah.

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpPTFTemplateCV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpPTF.h>
#include <lal/MatrixUtils.h>

NRCSID(FINDCHIRPPTFTEMPLATEC, "$Id$");

/* <lalVerbatim file="FindChirpPTFTemplateCP"> */
void
LALFindChirpPTFTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *InspTmplt,
    FindChirpTmpltParams       *params
    )
/* </lalVerbatim> */
{
  UINT4 errcode;
  /* local variables */
  UINT4 i, N;
  REAL4 phi, omega_2_3, e1x, e1y, e1z, e2x, e2y, e2z, sqrtoftwo,
        onebysqrtoftwo, onebysqrtofsix;
  REAL4Vector Q[5];
  COMPLEX8Vector Qtilde[5];

  sqrtoftwo      = sqrt(2.0);
  onebysqrtoftwo = 1.0 / sqrtoftwo;
  onebysqrtofsix = 1.0 / sqrt(6.0);
  N = params->PTFphi->length;

  INITSTATUS( status, "LALFindChirpPTFTemplate", FINDCHIRPPTFTEMPLATEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */

  /* check that the output structures exist */
  ASSERT( fcTmplt, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( fcTmplt->PTFQtilde, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( fcTmplt->PTFQtilde->length == 5, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( fcTmplt->PTFQtilde->data, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the parameter structure exists */
  ASSERT( params, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->PTFQ, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->PTFQ->length == 5, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->PTFQ->data, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  ASSERT( params->fwdPlan, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the timestep is positive */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );

  /* check that the input exists */
  ASSERT( InspTmplt, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  if ( params->approximant != FindChirpPTF )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }
  LALInfo( status, "Generating template using FindChirpPTF" );

  /* copy the template parameters to the finchirp template structure */
  memcpy( &(fcTmplt->tmplt), InspTmplt, sizeof(InspiralTemplate) );
  fcTmplt->tmplt.approximant = params->approximant;

  /* XXX delete this line if the low frequency cutoff XXX */
  /* XXX should be read from the template bank        XXX */
  InspTmplt->fLower = fcTmplt->tmplt.fLower = params->fLow;

  /* Zero out the Q and Qtilde vectors */
  memset( params->PTFQ->data, 0, 5 * N * sizeof(REAL4) );
  memset( fcTmplt->PTFQtilde->data, 0, 5 * (N /2 + 1) * sizeof(COMPLEX8) );

 /* Point the dummy variables Q and Qtilde to the actual output structures */
  for ( i = 0; i < 5; ++i )
  {
    Q[i].length      = N;
    Qtilde[i].length = N / 2 + 1;
    Q[i].data        = params->PTFQ->data + (i * N);
    Qtilde[i].data   = fcTmplt->PTFQtilde->data + (i * (N / 2 + 1)) ;
  }


  /* call the waveform generation function */

  errcode = XLALFindChirpPTFWaveform( params->PTFphi, params->PTFomega_2_3,
                                      params->PTFe1, params->PTFe2, InspTmplt,
                                      params->deltaT);
  if ( errcode != XLAL_SUCCESS )
  {
    ABORT( status, FINDCHIRPH_EPTFW, FINDCHIRPH_MSGEPTFW );
  }


  /* evaluate the Q^I factors from the dynamical variables */
  for( i = 0; i < N; ++i)
  {
    omega_2_3 = params->PTFomega_2_3->data[i];
    phi       = params->PTFphi->data[i];
    e1x       = params->PTFe1->data[i];
    e1y       = params->PTFe1->data[N + i];
    e1z       = params->PTFe1->data[2 * N + i];
    e2x       = params->PTFe2->data[i];
    e2y       = params->PTFe2->data[N + i];
    e2z       = params->PTFe2->data[2 * N + i];

    Q[0].data[i] = omega_2_3 * onebysqrtoftwo * ( cos(2 * phi) * ( e1x * e1x +
          e2y * e2y - e2x * e2x - e1y * e1y ) + 2 * sin(2 * phi) *
        ( e1x * e2x - e1y * e2y ));
    Q[1].data[i] = omega_2_3 * sqrtoftwo * ( cos(2 * phi) * ( e1x * e1y -
          e2x * e2y ) + sin(2 * phi) * ( e1x * e2y + e1y * e2x ));
    Q[2].data[i] = omega_2_3 * sqrtoftwo * ( cos(2 * phi) * ( e1x * e1z -
          e2x * e2z ) + sin(2 * phi) * ( e1x * e2z + e1z * e2x ));
    Q[3].data[i] = omega_2_3 * sqrtoftwo * ( cos(2 * phi) * ( e1y * e1z -
          e2y * e2z ) + sin(2 * phi) * ( e1y * e2z + e1z * e2y ));
    Q[4].data[i] = omega_2_3 * onebysqrtofsix * ( cos(2 * phi) *
        ( 2 * e2z * e2z - 2 * e1z * e1z + e1x * e1x + e1y * e1y -
          e2x * e2x - e2y * e2y ) + 2 * sin(2 * phi) * ( e1x * e2x +
            e1y * e2y - 2 * e1z * e2z ));
  }


  /* Fourier transform the Q's into the Qtilde's */
  for ( i = 0; i < 5; ++i )
  {
    LALForwardRealFFT( status->statusPtr, &Qtilde[i], &Q[i],
        params->fwdPlan);
  }

  /* XXX set this to be the correct values XXX */
  fcTmplt->tmplt.tC = InspTmplt->tC; /* length of template in seconds */
  fcTmplt->tmplt.fFinal = InspTmplt->fFinal; /* upper freq of template in Hz */

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpPTFTemplateCP"> */
void
LALFindChirpPTFNormalize(
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    FindChirpSegment           *fcSeg,
    FindChirpDataParams        *params
    )
/* </lalVerbatim> */
{
  UINT4         i, j, k, kmin, len, kmax;
  REAL4         f_min, deltaT, deltaF, fFinal;
  REAL4        *det         = NULL;
  REAL4        *PTFB        = NULL;
  COMPLEX8     *wtilde      = NULL;
  COMPLEX8     *PTFQtilde   = NULL;

  /* wtilde contains 1/S_n(f) (up to dynamic range issues) */
  wtilde    = params->wtildeVec->data;
  PTFQtilde = fcTmplt->PTFQtilde->data;
  PTFB      = fcTmplt->PTFB->data;
  len       = params->wtildeVec->length;
  deltaT    = (REAL4) fcSeg->deltaT;
  deltaF    = 1.0 / ( deltaT * 2 * ( (REAL4)len - 1) );
  f_min     = (REAL4) fcTmplt->tmplt.fLower;
  kmin      = f_min / deltaF > 1 ?  f_min / deltaF : 1;
  fFinal    = (REAL4) fcTmplt->tmplt.fFinal;
  kmax      = fFinal / deltaF < (len - 1) ? fFinal / deltaF : (len - 1);

  INITSTATUS( status, "LALFindChirpPTFNormalize", FINDCHIRPPTFTEMPLATEC );
  ATTATCHSTATUSPTR( status );

  /* check the required input exists */
  ASSERT( fcTmplt, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( fcSeg, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  ASSERT( params, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  ASSERT( params->wtildeVec, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->wtildeVec->data, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the parameter structure is set to a time domain approximant */
  if ( params->approximant != FindChirpPTF )
  {
    ABORT( status, FINDCHIRPH_EMAPX, FINDCHIRPH_MSGEMAPX );
  }

  /* output is B_{IJ}^{-1} */

  /*
   *
   * compute PTF normalization matrix
   *
   */

  /* Zero out the element sof matrix B and Binverse*/
  memset( fcTmplt->PTFB->data, 0, 25 * sizeof(REAL4) );
  memset( fcTmplt->PTFBinverse->data, 0, 25 * sizeof(REAL4) );

  /* Compute B_ij from Qtilde_i and Qtilde_j */
  for( i = 0; i < 5; ++i )
  {
    for ( j = 0; j < i + 1; ++j )
    {
      for ( k = kmin; k < kmax ; ++k )
      {
        PTFB[5 * i + j] += (PTFQtilde[k + i * len].re *
                            PTFQtilde[k + j * len].re +
                            PTFQtilde[k + i * len].im *
                            PTFQtilde[k + j * len].im )
                            * wtilde[k].re ;
      }
      PTFB[5 * i + j] *= 4.0 * deltaF ;
      /* Use the symmetry of B */
      PTFB[5 * j + i] = PTFB[5 * i + j];
    }
  }

  /* Invert B and store the output in Binverse */
  LALSMatrixInverse ( status->statusPtr,
      det, fcTmplt->PTFB, fcTmplt->PTFBinverse );
  CHECKSTATUSPTR( status );

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
