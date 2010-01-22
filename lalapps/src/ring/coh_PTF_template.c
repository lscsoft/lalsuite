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

#include "coh_PTF.h"

#include "lalapps.h"
#include "errutil.h"

RCSID( "$Id$" );

NRCSID(FINDCHIRPPTFTEMPLATEC, "$Id: FindChirpPTFTemplate.c,v 1.7 2008/06/26 19:05:07 dfazi Exp $");

#define sanity_check( condition ) \
  ( condition ? 0 : ( fputs( #condition " not satisfied\n", stderr ), error( "sanity check failed\n" ) ) ) 

/* <lalVerbatim file="FindChirpPTFTemplateCP"> */
void
cohPTFTemplate (
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *InspTmplt,
    FindChirpTmpltParams       *params
    )
/* </lalVerbatim> */
{
/*  LALStatus status;
  LAL_CALL (LALFindChirpPTFTemplate( &status,fcTmplt,InspTmplt,params), &status);*/
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

  /*
   *
   * check that the arguments are reasonable
   *
   */

  /* check that the output structures exist */
  sanity_check( fcTmplt );
  sanity_check( fcTmplt->PTFQtilde );
  sanity_check( fcTmplt->PTFQtilde->length == 5 );
  sanity_check( fcTmplt->PTFQtilde->data );

  /* check that the parameter structure exists */
  sanity_check( params );
  sanity_check( params->PTFQ );
  sanity_check( params->PTFQ->length == 5 );
  sanity_check( params->PTFQ->data );

  sanity_check( params->fwdPlan );

  /* check that the timestep is positive */
  sanity_check( params->deltaT > 0 );

  /* check that the input exists */
  sanity_check( InspTmplt );

  N = params->PTFphi->length;

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  sanity_check( InspTmplt->approximant == FindChirpPTF );
  sanity_check( InspTmplt->fLower );

  /* copy the template parameters to the finchirp template structure */
  memcpy( &(fcTmplt->tmplt), InspTmplt, sizeof(InspiralTemplate) );

  /* XXX delete this line if the low frequency cutoff XXX */
  /* XXX should be read from the template bank        XXX */
  fcTmplt->tmplt.fLower = params->fLow = InspTmplt->fLower;

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
  sanity_check( errcode == XLAL_SUCCESS );

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

  /* This block outputs a file containing all 5 Q time series (and time) */
  /* Need to make a option to print this if required*/
  /*FILE *outfile;
  outfile = fopen("Q_timeseries.dat","w");
  for ( i = 0; i < N; ++i)
  {
    fprintf (outfile,"%f %f %f %f %f %f \n",params->deltaT*i,Q[0].data[i],Q[1].data[i],Q[2].data[i],Q[3].data[i],Q[4].data[i]);
  }
  fclose(outfile);*/
  

  /* Fourier transform the Q's into the Qtilde's */
  for ( i = 0; i < 5; ++i )
  {
    XLALREAL4ForwardFFT( &Qtilde[i], &Q[i],
        params->fwdPlan);
  }

  /* This block outputs a file containing all 5 Q frequency series (and time) */
  /* It prints only the magnitude of the frequency */
  /* Need to make a option to print this if required and use deltaF not T*/
  /*FILE *outfile2;
  outfile2 = fopen("Q_freqseries.dat","w");
  for ( i = 0; i < N/2 + 1; ++i)
  {
    fprintf (outfile2,"%f %f %f %f %f %f %f %f %f %f %f\n",1./(N*params->deltaT)*i,Qtilde[0].data[i].re,Qtilde[0].data[i].im,Qtilde[1].data[i].re,Qtilde[1].data[i].im,Qtilde[2].data[i].re,Qtilde[2].data[i].im,Qtilde[3].data[i].re,Qtilde[3].data[i].im,Qtilde[4].data[i].re,Qtilde[4].data[i].im);
  }
  fclose(outfile2);*/


  /* XXX set this to be the correct values XXX */
  fcTmplt->tmplt.tC = InspTmplt->tC; /* length of template in seconds */
  fcTmplt->tmplt.fFinal = InspTmplt->fFinal; /* upper freq of template in Hz */

}


/* <lalVerbatim file="FindChirpPTFTemplateCP"> */
void
cohPTFNormalize(
    FindChirpTemplate          *fcTmplt,
    REAL4FrequencySeries       *invspec,
    REAL8Array                 *PTFM,
    REAL8Array                 *PTFN,
    COMPLEX8VectorSequence     *PTFqVec,
    COMPLEX8FrequencySeries    *sgmnt,
    COMPLEX8FFTPlan            *invPlan,
    UINT4                      spinTemplate
    )
/* Note that I use a different notation than Diego. In our notation M is
 * the same as Diego's B matrix (Q_0 | Q_0). We also have
 * A = ( Q_0 | s ) and B = ( Q_{\pi/2} | s ) Diego's A matrix is equal to
 * sqrt(A**2 + B**2) in this notation. */

/* </lalVerbatim> */
{
  UINT4         i, j, k, kmin, len, kmax,numPoints,vecLength;
  REAL8         f_min, deltaF,deltaT, fFinal, r, s, x, y, length;
  COMPLEX8     *PTFq, *qtilde, *inputData;
  COMPLEX8Vector *qtildeVec,qVec;
  REAL4        *det         = NULL;
  COMPLEX8     *PTFQtilde   = NULL;
  

  /* check the required input exists */
  sanity_check( fcTmplt );
  sanity_check( PTFM );
  sanity_check( PTFM->data );
  sanity_check( PTFqVec );
  sanity_check( PTFqVec->data );
  sanity_check( invspec );
  sanity_check( invspec->data );
  sanity_check( invspec->deltaF );
  sanity_check( sgmnt );
  sanity_check( sgmnt->data );
  sanity_check( invPlan );

  /* wtilde contains 1/S_n(f) (up to dynamic range issues) */
  numPoints   = PTFqVec->vectorLength;
  PTFQtilde = fcTmplt->PTFQtilde->data;
  len       = invspec->data->length;
  deltaF    = invspec->deltaF;
  deltaT    = 1.0 / ( deltaF * (REAL4) numPoints);
  /* This is explicit as I want f_min of template lower than f_min of filter*/
  f_min     = (REAL4) fcTmplt->tmplt.fLower;
  f_min     = 40.;
  kmin      = f_min / deltaF > 1 ?  f_min / deltaF : 1;
  fFinal    = (REAL4) fcTmplt->tmplt.fFinal;
  fFinal    = 200.;
  kmax      = fFinal / deltaF < (len - 1) ? fFinal / deltaF : (len - 1);
  qVec.length = numPoints;
  qtildeVec    = XLALCreateCOMPLEX8Vector( numPoints );
  qtilde = qtildeVec->data;

  /* Data params */
  inputData   = sgmnt->data->data;
  length      = sgmnt->data->length;

  /* Check that these input values are sane */
  sanity_check ( deltaT > 0 );
  sanity_check ( deltaF > 0 );
  sanity_check ( fcTmplt->tmplt.tC > 0 );
  /*Segment, response function and PTFQtilde must be the same length */
  sanity_check ( len == length ) ; 
  sanity_check ( fcTmplt->PTFQtilde->vectorLength == len);

  /* check that the parameter structure is set to a time domain approximant */
  sanity_check ( fcTmplt->tmplt.approximant == FindChirpPTF );

  /* Set parameters to determine spin/nonspin */
  vecLength = 2;
  if (spinTemplate == 1)
  {
    vecLength = 5;
  } 
  else
  {
    for ( k = kmin; k < kmax ; ++k )
    {
      PTFQtilde[k].im = PTFQtilde[k+len].re;
      PTFQtilde[k + len].im = -PTFQtilde[k].re;
    }
  }

  /*
   *
   * compute PTF normalization matrix
   *
   */

  /* Compute M_ij from Qtilde_i and Qtilde_j */
  for( i = 0; i < vecLength; ++i )
  {
    for ( j = 0; j < i + 1; ++j )
    {
      for ( k = kmin; k < kmax ; ++k )
      {
        PTFM->data[5 * i + j] += (PTFQtilde[k + i * len].re *
                            PTFQtilde[k + j * len].re +
                            PTFQtilde[k + i * len].im *
                            PTFQtilde[k + j * len].im )
                            * invspec->data->data[k] ;
      }
      PTFM->data[5 * i + j] *= 4.0 * deltaF ;
      /* Use the symmetry of M */
      PTFM->data[5 * j + i] = PTFM->data[5 * i + j];
    }
  }
  
  if (PTFN)
  {
    /* Compute N_ij */
    for( i = 0; i < vecLength; ++i )
    {
      for ( j = 0; j < i + 1; ++j )
      {
        for ( k = kmin; k < kmax ; ++k )
        {
          PTFN->data[5 * i + j] += (PTFQtilde[k + i * len].re *
                              PTFQtilde[k + j * len].im -
                              PTFQtilde[k + i * len].im *
                              PTFQtilde[k + j * len].re )
                              * invspec->data->data[k] ;
        }
        PTFN->data[5 * i + j] *= 4.0 * deltaF ;
        /* Use the anti-symmetry of M */
        PTFN->data[5 * j + i] = -PTFN->data[5 * i + j];
      }
    }
  }


  /* A routine to print out the M^IJ information */
 
  FILE *outfile;
  /*
  outfile = fopen("M_array.dat","w");
  for ( i = 0; i < 5; ++i )
  {
    for ( j = 0; j < 5; ++j )
    {
      fprintf(outfile,"%f ", PTFM->data[5*i + j]);
    }
  fprintf(outfile,"\n");
  }
  fclose(outfile);
  */
  

  for ( i = 0; i < vecLength; ++i )
  {

    /* compute qtilde using data and Qtilde */

    memset( qtildeVec->data, 0,
        qtildeVec->length * sizeof(COMPLEX8) );

    /* qtilde positive frequency, not DC or nyquist */
    for ( k = kmin; k < kmax ; ++k )
    {
      r = inputData[k].re;
      s = inputData[k].im;
      x = PTFQtilde[i * (numPoints / 2 + 1) + k].re;
      y = 0 - PTFQtilde[i * (numPoints / 2 + 1) + k].im; /* cplx conj */

      qtilde[k].re = 4. * (r*x - s*y)*deltaF;
      qtilde[k].im = 4. * (r*y + s*x)*deltaF;
    }

    qVec.data = PTFqVec->data + (i * numPoints);

    /* inverse fft to get q */
    XLALCOMPLEX8VectorFFT( &qVec, qtildeVec, invPlan );
  }

  /*Routines for printing the A and B time series */
  
  /*
  outfile = fopen("A_timeseries.dat","w");
  for ( i = 0; i < numPoints; ++i)
  {
    fprintf (outfile,"%f %f %f %f %f %f\n",deltaT*i,PTFqVec->data[i].re,PTFqVec->data[i+numPoints].re,PTFqVec->data[i+2*numPoints].re,PTFqVec->data[i+3*numPoints].re,PTFqVec->data[i+4*numPoints].re);
  }
  fclose(outfile);

  outfile = fopen("B_timeseries.dat","w");
  for ( i = 0; i < numPoints; ++i)
  {
    fprintf (outfile,"%f %f %f %f %f %f\n",deltaT*i,PTFqVec->data[i].im,PTFqVec->data[i+numPoints].im,PTFqVec->data[i+2*numPoints].im,PTFqVec->data[i+3*numPoints].im,PTFqVec->data[i+4*numPoints].im);
  }
  fclose(outfile);
  */
  XLALDestroyCOMPLEX8Vector( qtildeVec );
}
