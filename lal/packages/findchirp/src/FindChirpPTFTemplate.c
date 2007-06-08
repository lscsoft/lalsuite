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
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>

NRCSID(FINDCHIRPPTFTEMPLATEC, "$Id$");

/* <lalVerbatim file="FindChirpPTFTemplateCP"> */
void
LALFindChirpPTFTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpTmpltParams       *params
    )
/* </lalVerbatim> */
{
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
  ASSERT( tmplt, status, 
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  if ( params->approximant != FindChirpPTF )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }
  LALInfo( status, "Generating template using FindChirpPTF" );

  fprintf( stderr, 
      "LALFindChirpPTFTemplate() called with "
      "m1 = %e, m2 = %e, chi = %e, kappa = %e\n", 
      tmplt->mass1, tmplt->mass2, tmplt->chi, tmplt->kappa );

  /* copy the template parameters to the finchirp template structure */
  memcpy( &(fcTmplt->tmplt), tmplt, sizeof(InspiralTemplate) );
  fcTmplt->tmplt.approximant = params->approximant;
  fcTmplt->tmplt.tC = 25.0; /* length of template in seconds */
  fcTmplt->tmplt.fFinal = 1000.0; /* upper freq of template in Hz */

  /* local variables */
  UINT4 i, j, len;
  REAL4 phi, omega_2_3, e1x, e1y, e1z, e2x, e2y, e2z, sqrtoftwo, 
        onebysqrtoftwo, onebysqrtofsix;
  REAL4Vector Q[5];
  COMPLEX8Vector Qtilde[5];
  
  sqrtoftwo      = sqrt(2.0);
  onebysqrtoftwo = 1.0 / sqrtoftwo;
  onebysqrtofsix = 1.0 / sqrt(6.0);
  len = params->PTFphi->length;
  
  /* Point the dummy variables Q and Qtilde to the actual output structures */
  for (i=0; i<5; i++)
  {
    Q[i].length      = len;
    Qtilde[i].length = len / 2 + 1;
    Q[i].data        = params->PTFQ->data + i * len;
    Qtilde[i].data   = fcTmplt->PTFQtilde->data + i * (len / 2 + 1) ;
  }  
  
  /* evaluate the Q^I factors from the dynamical variables */
  for( i=0; i < len; i++)
  {
    omega_2_3 = params->PTFomega_2_3->data[i];
    phi       = params->PTFphi->data[i];
    e1x       = params->PTFe1->data[i];
    e1y       = params->PTFe1->data[len + i];
    e1z       = params->PTFe1->data[2 * len + i]; 
    e2x       = params->PTFe2->data[i];
    e2y       = params->PTFe2->data[len + i];
    e2z       = params->PTFe2->data[2 * len + i];
    
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
  for ( i=0; i<5; i++ )
  {
    LALForwardRealFFT( status->statusPtr, &Qtilde[i], &Q[i],
                       params->fwdPlan);
  }
  
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
  UINT4         i, j, k, kmin, len;
  REAL4         deltaT, fmin;
  REAL8         deltaF;
  REAL4        *det         = NULL;
  REAL4        *PTFB        = NULL; 
  COMPLEX8     *wtilde      = NULL;
  COMPLEX8     *PTFQtilde   = NULL;

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

  /* wtilde contains 1/S_n(f) (up to dynamic range issues) */
  wtilde    = params->wtildeVec->data;
  PTFQtilde = fcTmplt->PTFQtilde->data;
  PTFB      = fcTmplt->PTFB->data;
  len       = params->wtildeVec->length; 
  deltaT    = fcSeg->deltaT;
  deltaF    = 1.0 / ( (REAL4) deltaT * (REAL4) len);
  fmin      = fcTmplt->tmplt.fLower;
  kmin      = fmin / deltaF > 1 ?  fmin / deltaF : 1;
  
  
  /* output is B_{IJ}^{-1} */

  fprintf( stderr, 
      "LALFindChirpPTFNormalize() called with wtilde at %p, PTFQtilde at %p\n",
      wtilde, PTFQtilde );

  /*
   *
   * compute PTF normalization matrix
   *
   */

  /* Compute B_ij from Qtilde_i and Qtilde_j */
  for( i=0; i<5; i++ )
  {
    for ( j=0; j<i+1; j++ )
    {  
        for ( k=kmin; k<len; ++k )
        {  
          PTFB[5 * i + j] += (PTFQtilde[k + i * len].re * 
                              PTFQtilde[k + j * len].re +
                              PTFQtilde[k + i * len].im * 
                              PTFQtilde[k + j * len].im ) / 
                              wtilde[k].re ;
        }    
      PTFB[5 * i + j] *= 4.0 * deltaF ;
      /* Use the symmetry of B */
      PTFB[5 * i + j] = PTFB[5 * j + i];
    }    
  }  

  /* Invert B and store the outut in Binverse */
  LALSMatrixInverse ( status, det, fcTmplt->PTFB, fcTmplt->PTFBinverse );

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}




