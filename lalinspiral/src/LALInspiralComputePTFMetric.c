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

/**
 * \author Yi Pan, Duncan Brown
 * \file
 * \ingroup LALInspiralBank_h
 *
 * \brief Module to compute the components of the metric which is used to describe
 * distances on Physical Template Family signal manifold.
 *
 * ### Prototypes ###
 *
 * <tt>XLALInspiralComputePTFIntrinsicMetric()</tt>:
 * <ul>
 * <li> <tt>metric,</tt> Output, the metric at the lattice point defined by \c params
 * </li><li> <tt>psd,</tt> Input, the power spectral density of the data
 * </li><li> <tt>params,</tt> Input, the parameters where metric must be computed
 * in the computation of the metric.</li>
 * </ul>
 *
 * <tt>XLALInspiralComputePTFFullMetric()</tt>:
 * <ul>
 * <li> <tt>metric,</tt> Output, the metric at the lattice point defined by \c params
 * </li><li> <tt>psd,</tt> Input, the power spectral density of the data
 * </li><li> <tt>params,</tt> Input, the parameters where metric must be computed
 * in the computation of the metric.</li>
 * </ul>
 *
 * <tt>XLALInspiralComputePTFWaveform()</tt>:
 * <ul>
 * <li> <tt>ptfwave,</tt> Output, the waveform at the lattice point defined
 * by \c params
 * </li><li> <tt>params,</tt> Input, the parameters where metric must be computed
 * in the computation of the metric.</li>
 * </ul>
 *
 * <tt>XLALInspiralComputePTFWDeriv()</tt>:
 * <ul>
 * <li> <tt>Wderiv,</tt> Output, the time derivative of waveform at the lattice
 * point defined by \c params
 * </li><li> <tt>psd,</tt>  Input, the power spectral density of the data
 * </li><li> <tt>params,</tt> Input, the parameters where metric must be computed
 * in the computation of the metric
 * </li><li> <tt>paramid,</tt> Input, id of the parameter to take derivative on
 * </li><li> <tt>initdelta,</tt> Input, initial difference in parameters
 * </li><li> <tt>tolerance,</tt> Input, stop iteration when difference between two
 * bisections is smaller than tolerance.</li>
 * </ul>
 *
 * <tt>XLALInspiralComputePTFQDeriv()</tt>:
 * <ul>
 * <li> <tt>Qderiv,</tt> Output, the time derivative of Q at the lattice point
 * defined by \c params
 * </li><li> <tt>params,</tt> Input, the parameters where metric must be computed
 * in the computation of the metric.</li>
 * </ul>
 *
 * ### Description ###
 *
 * We calculate the components of the metric using the procedure outlined
 * by Yi.
 *
 * ### Algorithm ###
 *
 *
 * ### Uses ###
 *
 * \code
 * LALMalloc
 * LALFree
 * \endcode
 *
 * ### Notes ###
 *
 */

#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <stdio.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>
#include <lal/MatrixUtils.h>
#include <lal/FindChirpPTF.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


INT4 XLALInspiralComputePTFIntrinsicMetric (
    InspiralMetric			   *metric,
    REAL8Vector				   *fullmetric,
    REAL8FrequencySeries                   *psd,
    InspiralTemplate                       *params
    )

{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;

  /* number of points in a time-domain segment */
  UINT4 N = 2 * (psd->data->length - 1);
  UINT4 i, j, k;

  /* some useful numbers */
  REAL4 sqrtoftwo      = sqrt(2.0);
  REAL4 onebysqrtoftwo = 1.0 / sqrtoftwo;
  REAL4 onebysqrtofsix = 1.0 / sqrt(6.0);

  /* get a local copy of the extrinstic parameters */
  REAL8 Theta = params->sourceTheta;
  REAL8 Phi = params->sourcePhi;
  REAL8 Psi = params->polarisationAngle;
  REAL8 phi0 = params->startPhase;

  /* bounds on the power spectrum integration for the moments */
  REAL8 deltaT = params->tSampling;
  REAL8 deltaF = 1.0 / ((REAL8)N * deltaT);

  /* local pointers to the Q and Qtilde data */
  REAL8				P[5];
  REAL8				PdTh[5];
  REAL8				PdPh[5];
  REAL8				PdPs[5];
  REAL8Vector		        Q[5];
  COMPLEX16Vector	        Qtilde[5];

  FILE *derivsfile;
  /* should really check that deltaT from the template agrees with N */
  /* and deltaF from the power spectrum                              */

  /* these variables, and the memory allocation following them should be     */
  /* moved to a differnent function at some point in the future (preferably  */
  /* before this function is ever called multiple times inside a loop)       */
  REAL8VectorSequence			*PTFQ;
  COMPLEX16VectorSequence		*PTFQtilde;
  REAL4Vector				*PTFphi;
  REAL4Vector				*PTFomega_2_3;
  REAL4VectorSequence			*PTFe1;
  REAL4VectorSequence			*PTFe2;
  REAL8FFTPlan				*fwdPlan;
  REAL8FFTPlan				*revPlan;
  COMPLEX16Vector			*fsig0;
  COMPLEX16Vector			*fsig1;
  COMPLEX16Vector			*fsig;
  COMPLEX16Vector			*intrinsicderiv;
  COMPLEX16VectorSequence	        *derivs;
  REAL8		                         initdelta = 0.001;
  REAL8		                         tolerance = 0.001;
  REAL8Vector	                 	*invpsd;
  REAL8Vector	                	*tempnorm;
  COMPLEX16Vector               	*tempnormtilde;
  REAL8                                  wavenorm;
  REAL8Array                            *IntInt_matrix;
  REAL8Array                            *ExtExt_matrix;
  REAL8Array                            *ExtExt_inverse;
  REAL8Array                            *IntExt_matrix;
  REAL8Array                            *IntExt_transp;
  REAL8Array                            *matrix_5_x_4;
  REAL8Array                            *matrix_4_x_4;
  REAL8                                 *det;
  LALStatus                              status;

  PTFQ           = XLALCreateREAL8VectorSequence( 5, N );
  PTFQtilde      = XLALCreateCOMPLEX16VectorSequence( 5, N / 2 + 1 );
  PTFphi         = XLALCreateVector( N );
  PTFomega_2_3   = XLALCreateVector( N );
  PTFe1          = XLALCreateVectorSequence( 3, N );
  PTFe2          = XLALCreateVectorSequence( 3, N );
  fwdPlan        = XLALCreateForwardREAL8FFTPlan( N, 0 );
  revPlan        = XLALCreateReverseREAL8FFTPlan( N, 0 );
  IntInt_matrix  = XLALCreateREAL8ArrayL ( 2, 4, 4 );
  ExtExt_matrix  = XLALCreateREAL8ArrayL ( 2, 5, 5 );
  ExtExt_inverse = XLALCreateREAL8ArrayL ( 2, 5, 5 );
  IntExt_matrix  = XLALCreateREAL8ArrayL ( 2, 4, 5 );
  IntExt_transp  = XLALCreateREAL8ArrayL ( 2, 5, 4 );
  matrix_5_x_4   = XLALCreateREAL8ArrayL ( 2, 5, 4 );
  matrix_4_x_4   = XLALCreateREAL8ArrayL ( 2, 4, 4 );

  det            = NULL;

  memset( &status, 0, sizeof(LALStatus) );

  /* call the PTF waveform code */
  errcode = XLALFindChirpPTFWaveform( PTFphi, PTFomega_2_3, PTFe1, PTFe2,
      params, deltaT);

  if ( errcode != XLAL_SUCCESS ) XLAL_ERROR( errcode );

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
    REAL8 omega_2_3 = PTFomega_2_3->data[j];
    REAL8 phi       = PTFphi->data[j];
    REAL8 e1x       = PTFe1->data[j];
    REAL8 e1y       = PTFe1->data[N + j];
    REAL8 e1z       = PTFe1->data[2 * N + j];
    REAL8 e2x       = PTFe2->data[j];
    REAL8 e2y       = PTFe2->data[N + j];
    REAL8 e2z       = PTFe2->data[2 * N + j];

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
    XLALREAL8ForwardFFT( &Qtilde[i], &Q[i], fwdPlan);
  }

  /* evaluate the P_I factors from the extrinsic parameters */
  P[0] = -(1 + cos(Theta) * cos(Theta)) * cos(2 * Phi) * cos(Psi)
    + 2 * cos(Theta) * sin(2 * Phi) * sin(Psi);
  P[1] = - 2 * cos(Theta) * cos(2 * Phi) * sin(Psi)
    -(1 + cos(Theta) * cos(Theta)) * sin(2 * Phi) * cos(Psi);
  P[2] = 2 * sin(Theta) * (-sin(Phi) * sin(Psi) + cos(Theta) * cos(Phi) * cos(Psi));
  P[3] = 2 * sin(Theta) * ( cos(Phi) * sin(Psi) + cos(Theta) * sin(Phi) * cos(Psi));
  P[4] = 3 * sin(Theta) * sin(Theta) * cos(Psi);

  /* P_I derivative over Theta */
  PdTh[0] = -(1 - sin(2 * Theta)) * cos(2 * Phi) * cos(Psi)
    - 2 * sin(Theta) * sin(2 * Phi) * sin(Psi);
  PdTh[1] =  2 * sin(Theta) * cos(2 * Phi) * sin(Psi)
    -(1 - sin(2 * Theta)) * sin(2 * Phi) * cos(Psi);
  PdTh[2] = 2 * cos(Theta) * (-sin(Phi) * sin(Psi) - sin(Theta) * cos(Phi) * cos(Psi));
  PdTh[3] = 2 * cos(Theta) * ( cos(Phi) * sin(Psi) - sin(Theta) * sin(Phi) * cos(Psi));
  PdTh[4] = 3 * sin(2 * Theta) * cos(Psi);

  /* P_I derivative over Phi */
  PdPh[0] = - 2 * (1 + cos(Theta) * cos(Theta)) * sin(2 * Phi) * cos(Psi)
    + 4 * cos(Theta) * cos(2 * Phi) * sin(Psi);
  PdPh[1] =   4 * cos(Theta) * sin(2 * Phi) * sin(Psi)
    - 2 * (1 + cos(Theta) * cos(Theta)) * cos(2 * Phi) * cos(Psi);
  PdPh[2] = 2 * sin(Theta) * (-cos(Phi) * sin(Psi) - cos(Theta) * sin(Phi) * cos(Psi));
  PdPh[3] = 2 * sin(Theta) * (-sin(Phi) * sin(Psi) + cos(Theta) * cos(Phi) * cos(Psi));
  PdPh[4] = 0;

  /* P_I derivative over Psi */
  PdPs[0] = (1 + cos(Theta) * cos(Theta)) * cos(2 * Phi) * sin(Psi)
    + 2 * cos(Theta) * sin(2 * Phi) * cos(Psi);
  PdPs[1] = - 2 * cos(Theta) * cos(2 * Phi) * cos(Psi)
    + (1 + cos(Theta) * cos(Theta)) * sin(2 * Phi) * sin(Psi);
  PdPs[2] = 2 * sin(Theta) * (-sin(Phi) * cos(Psi) - cos(Theta) * cos(Phi) * sin(Psi));
  PdPs[3] = 2 * sin(Theta) * ( cos(Phi) * cos(Psi) - cos(Theta) * sin(Phi) * sin(Psi));
  PdPs[4] = - 3 * sin(Theta) * sin(Theta) * sin(Psi);

  /* Prepare frequency domain signals */
  fsig0 = XLALCreateCOMPLEX16Vector(N / 2 + 1);
  fsig1 = XLALCreateCOMPLEX16Vector(N / 2 + 1);
  fsig  = XLALCreateCOMPLEX16Vector(N / 2 + 1);

  for (k = 0; k < N / 2 + 1; ++k)
  {
    fsig0->data[k] = 0.0;
    fsig1->data[k] = 0.0;
    for (i = 0; i < 5; ++i)
    {
      fsig0->data[k] += (((REAL8) P[i]) * Qtilde[i].data[k]);
    }
    fsig1->data[k] = crect( cimag(fsig0->data[k]), - creal(fsig0->data[k]) );
  }
  fsig1->data[0] = fsig0->data[0];
  fsig1->data[N/2] = fsig0->data[N / 2];
  for (k = 0; k < N / 2 + 1; ++k)
  {
    fsig->data[k] = crect( cos(phi0) * creal(fsig0->data[k]) + sin(phi0) * creal(fsig1->data[k]), cos(phi0) * cimag(fsig0->data[k]) + sin(phi0) * cimag(fsig1->data[k]) );
  }

  /* Waveform derivatives */
  intrinsicderiv = XLALCreateCOMPLEX16Vector(N / 2 + 1);
  derivs = XLALCreateCOMPLEX16VectorSequence(9, N / 2 + 1);

  /* Compute extrinsic derivatives */
  for (k = 0; k < N / 2 + 1; ++k)
  {
    /* t0 */
    derivs->data[k] = crect( - LAL_TWOPI * deltaF * k * cimag(fsig->data[k]), LAL_TWOPI * deltaF * k * creal(fsig->data[k]) );
    /* phi0 */
    derivs->data[N/2+1+k] = crect( - sin(phi0) * creal(fsig0->data[k]) + cos(phi0) * creal(fsig1->data[k]), - sin(phi0) * cimag(fsig0->data[k]) + cos(phi0) * cimag(fsig1->data[k]) );
    /* Theta */
    derivs->data[N+2+k] = 0.0;
    for (i = 0; i < 5; ++i)
    {
      derivs->data[N+2+k] += (((REAL8) PdTh[i]) * Qtilde[i].data[k]);
    }
    /* Phi */
    derivs->data[3*N/2+3+k] = 0.0;
    for (i = 0; i < 5; ++i)
    {
      derivs->data[3*N/2+3+k] += (((REAL8) PdPh[i]) * Qtilde[i].data[k]);
    }
    /* Psi */
    derivs->data[2*N+4+k] = 0.0;
    for (i = 0; i < 5; ++i)
    {
      derivs->data[2*N+4+k] += (((REAL8) PdPs[i]) * Qtilde[i].data[k]);
    }
  }

  /* Compute intrinsic derivatives */
  for (i = 5; i < 9; ++i)
  {
    errcode = XLALInspiralComputePTFWDeriv(intrinsicderiv, psd, params, i - 4, initdelta, tolerance);
    if ( errcode != XLAL_SUCCESS )
    {
      fprintf( stderr, "XLALInspiralComputePTFDeriv failed\n" );
      exit( 1 );
    }
    for (k = 0; k < N / 2 + 1; ++k)
    {
      derivs->data[i*(N/2+1)+k] = crect( cos(phi0) * creal(intrinsicderiv->data[k]) + sin(phi0) * cimag(intrinsicderiv->data[k]), cos(phi0) * cimag(intrinsicderiv->data[k]) - sin(phi0) * creal(intrinsicderiv->data[k]) );
    }
    derivs->data[i*(N/2+1)] = crect( (cos(phi0) + sin(phi0)) * creal(intrinsicderiv->data[0]), (cos(phi0) + sin(phi0)) * cimag(intrinsicderiv->data[0]) );
    derivs->data[i*(N/2+1)+N/2] = crect( (cos(phi0) + sin(phi0)) * creal(intrinsicderiv->data[N / 2]), (cos(phi0) + sin(phi0)) * cimag(intrinsicderiv->data[N / 2]) );
  }


  derivsfile = fopen( "myderivs.dat", "w" );

  if( derivsfile != NULL )
  {
    for (j = 0; j < N / 2 + 1; ++j)
    {
      fprintf( derivsfile, "%f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f   %f\n",
          creal(derivs->data[j]), cimag(derivs->data[j]),
          creal(derivs->data[(N / 2 + 1) + j]), cimag(derivs->data[(N / 2 + 1) + j]),
          creal(derivs->data[2 * (N / 2 + 1) + j]), cimag(derivs->data[2 * (N / 2 + 1) + j]),
          creal(derivs->data[3 * (N / 2 + 1) + j]), cimag(derivs->data[3 * (N / 2 + 1) + j]),
          creal(derivs->data[4 * (N / 2 + 1) + j]), cimag(derivs->data[4 * (N / 2 + 1) + j]),
          creal(derivs->data[5 * (N / 2 + 1) + j]), cimag(derivs->data[5 * (N / 2 + 1) + j]),
          creal(derivs->data[6 * (N / 2 + 1) + j]), cimag(derivs->data[6 * (N / 2 + 1) + j]),
          creal(derivs->data[7 * (N / 2 + 1) + j]), cimag(derivs->data[7 * (N / 2 + 1) + j]),
          creal(derivs->data[8 * (N / 2 + 1) + j]), cimag(derivs->data[8 * (N / 2 + 1) + j]));
    }
  }
  fclose(derivsfile);

  /* Compute full metric */
  invpsd		= XLALCreateREAL8Vector(N / 2 + 1);
  tempnorm		= XLALCreateREAL8Vector(N);
  tempnormtilde = XLALCreateCOMPLEX16Vector(N / 2 + 1);

  for (k = 0; k < N / 2 + 1; ++k)
  {
    if (psd->data->data[k] == 0.)
      invpsd->data[k] = 0.;
    else
      invpsd->data[k] = 1.0 / psd->data->data[k];
  }

  for (i = 0; i < 9; ++i)
  {
    for (j = 0; j < i + 1; ++j)
    {
      for (k = 0; k < N / 2 + 1; ++k)
      {
        tempnormtilde->data[k] = crect( ( creal(derivs->data[i * (N/2+1) + k]) * creal(derivs->data[j * (N/2+1) + k]) + cimag(derivs->data[i * (N/2+1) + k]) * cimag(derivs->data[j * (N/2+1) + k])) * invpsd->data[k], ( cimag(derivs->data[i * (N/2+1) + k]) * creal(derivs->data[j * (N/2+1) + k]) - creal(derivs->data[i * (N/2+1) + k]) * cimag(derivs->data[j * (N/2+1) + k])) * invpsd->data[k] );
      }
      /* Inverse Fourier of tempnorm */
      XLALREAL8ReverseFFT(tempnorm, tempnormtilde, revPlan);
      fullmetric->data[i * (i + 1) / 2 + j] = 4.0 * tempnorm->data[0];
    }
  }

  wavenorm = fullmetric->data[2];
  for (i = 0; i < 45; ++i)
  {
    fullmetric->data[i] /= wavenorm;
  }

  /* Calculating the extrinsic-extrinsic block of the metric */
  for (i = 0; i < 5; ++i)
  {
    for (j = 0; j < i + 1; ++j)
    {
      ExtExt_matrix->data[ 5 * i + j] = ExtExt_matrix->data[ 5 * j + i] =
        fullmetric->data[i * (i + 1) / 2 + j];
    }
  }

  /* Calculating the inverse of ext-ext */
  LALDMatrixInverse( &status, det, ExtExt_matrix, ExtExt_inverse );
  if ( status.statusCode )
  {
    REPORTSTATUS( &status );
    exit( 1 );
  }

  /* Calculating the intrinsic-intrinsic block of the metric */
  for (i = 5; i < 9; ++i)
  {
    for (j = 5; j < i + 1; ++j)
    {
      IntInt_matrix->data[ 4 * (i-5) + j - 5] = IntInt_matrix->data[ 4 * (j-5) + i - 5]=
        fullmetric->data[i * (i + 1) / 2 + j];
    }
  }

  /* Calculating the mixed intrinsic-extrinsic block of the metric */
  for (i = 5; i < 9; ++i)
  {
    for (j = 0; j < 5; ++j)
    {
      IntExt_matrix->data[ 5 * (i-5) + j] = fullmetric->data[i * (i + 1) / 2 + j];
    }
  }

  /* Calculating the transpose of int-ext */
  LALDMatrixTranspose ( &status, IntExt_transp, IntExt_matrix );
  if ( status.statusCode )
  {
    REPORTSTATUS( &status );
    exit( 1 );
  }

  for (i = 0; i < 5; ++i)
  {
    for (j=0; j < 4; ++j)
    {
      IntExt_transp->data[ 4 * i + j] = IntExt_matrix->data[ j * 5 + i];
    }
  }

  /* Compute the projected metric from eq.(72) of PBCV*/

  LALDMatrixMultiply ( &status, matrix_5_x_4, ExtExt_inverse, IntExt_transp );
  if ( status.statusCode )
  {
    REPORTSTATUS( &status );
    exit( 1 );
  }

  LALDMatrixMultiply ( &status, matrix_4_x_4, IntExt_matrix, matrix_5_x_4 );
  if ( status.statusCode )
  {
    REPORTSTATUS( &status );
    exit( 1 );
  }

  for ( i = 0; i < 4; ++i )
  {
    for ( j = 0; j < i + 1; ++j )
    {
      params->Gamma[i * (i + 1) / 2 + j] = metric->Gamma[i * (i + 1) / 2 + j] =
        IntInt_matrix->data[i * 4 + j] - matrix_4_x_4->data[i * 4 + j];
    }
  }

  /* this memory deallocation code should be moved to a separate function */
  XLALDestroyCOMPLEX16Vector(fsig0);
  XLALDestroyCOMPLEX16Vector(fsig1);
  XLALDestroyCOMPLEX16Vector(fsig);
  XLALDestroyCOMPLEX16Vector(intrinsicderiv);
  XLALDestroyCOMPLEX16VectorSequence(derivs);
  XLALDestroyCOMPLEX16Vector(tempnormtilde);
  XLALDestroyREAL8Vector(tempnorm);
  XLALDestroyREAL8VectorSequence( PTFQ );
  XLALDestroyCOMPLEX16VectorSequence( PTFQtilde );
  XLALDestroyVector( PTFphi );
  XLALDestroyVector( PTFomega_2_3 );
  XLALDestroyVectorSequence( PTFe1 );
  XLALDestroyVectorSequence( PTFe2 );
  XLALDestroyREAL8FFTPlan( fwdPlan );
  XLALDestroyREAL8FFTPlan( revPlan );
  XLALDestroyREAL8Array( IntInt_matrix );
  XLALDestroyREAL8Array( ExtExt_matrix );
  XLALDestroyREAL8Array( IntExt_matrix );
  XLALDestroyREAL8Array( IntExt_transp );
  XLALDestroyREAL8Array( ExtExt_inverse );
  XLALDestroyREAL8Array( matrix_5_x_4 );
  XLALDestroyREAL8Array( matrix_4_x_4 );

  /* normal exit */
  return errcode;
}


INT4 XLALInspiralComputePTFFullMetric (
    InspiralMetric             UNUSED *metric,
    REAL8FrequencySeries       UNUSED *psd,
    InspiralTemplate           UNUSED *params
    )

{
  return XLAL_SUCCESS;
}



INT4 XLALInspiralComputePTFWaveform (
    REAL8Vector				   *ptfwave,
    InspiralTemplate           *params
    )

{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;

  /* number of points in a time-domain segment */
  UINT4 N = ptfwave->length;
  UINT4 i, j;

  /* some useful numbers */
  REAL8 sqrtoftwo      = sqrt(2.0);
  REAL8 onebysqrtoftwo = 1.0 / sqrtoftwo;
  REAL8 onebysqrtofsix = 1.0 / sqrt(6.0);

  /* get a local copy of the extrinstic parameters */
  REAL8 Theta = params->sourceTheta;
  REAL8 Phi = params->sourcePhi;
  REAL8 Psi = params->polarisationAngle;
  REAL8 phi0 = params->startPhase;

  /* bounds on the power spectrum integration for the moments */
  REAL8 deltaT = params->tSampling;

  /* P-factors */
  REAL8 P[5];
  /* local pointers to the Q and Qtilde data */
  REAL8Vector Q[5];
  REAL8Vector Qp[5];

  /* should really check that deltaT from the template agrees with N */
  /* and deltaF from the power spectrum                              */

  /* these variables, and the memory allocation following them should be     */
  /* moved to a differnent function at some point in the future (preferably  */
  /* before this function is ever called multiple times inside a loop)       */
  REAL8VectorSequence          *PTFQ;
  REAL8VectorSequence          *PTFQp;
  REAL4Vector                  *PTFphi;
  REAL4Vector                  *PTFomega_2_3;
  REAL4VectorSequence          *PTFe1;
  REAL4VectorSequence          *PTFe2;
  REAL8FFTPlan                  *fwdPlan;

  PTFQ = XLALCreateREAL8VectorSequence( 5, N );
  PTFQp = XLALCreateREAL8VectorSequence( 5, N );
  PTFphi = XLALCreateREAL4Vector( N );
  PTFomega_2_3 = XLALCreateREAL4Vector( N );
  PTFe1 = XLALCreateREAL4VectorSequence( 3, N );
  PTFe2 = XLALCreateREAL4VectorSequence( 3, N );
  fwdPlan = XLALCreateForwardREAL8FFTPlan( N, 0 );

  /* call the PTF waveform code */
  errcode = XLALFindChirpPTFWaveform( PTFphi, PTFomega_2_3, PTFe1, PTFe2,
      params, deltaT);

  if ( errcode != XLAL_SUCCESS ) XLAL_ERROR( errcode );

  /* point the dummy variables Q and Qp to the actual output structures */
  for ( i = 0; i < 5; ++i )
  {
    Q[i].length      = N;
    Qp[i].length     = N;
    Q[i].data        = PTFQ->data + (i * N);
    Qp[i].data       = PTFQp->data + (i * N);
  }

  /* evaluate the Q^I factors from the dynamical variables */
  for ( j = 0; j < N; ++j )
  {
    REAL8 omega_2_3 = PTFomega_2_3->data[j];
    REAL8 phi       = PTFphi->data[j];
    REAL8 e1x       = PTFe1->data[j];
    REAL8 e1y       = PTFe1->data[N + j];
    REAL8 e1z       = PTFe1->data[2 * N + j];
    REAL8 e2x       = PTFe2->data[j];
    REAL8 e2y       = PTFe2->data[N + j];
    REAL8 e2z       = PTFe2->data[2 * N + j];

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

	Qp[0].data[j] = omega_2_3 * onebysqrtoftwo * ( - sin(2 * phi) * ( e1x * e1x +
          e2y * e2y - e2x * e2x - e1y * e1y ) + 2 * cos(2 * phi) *
        ( e1x * e2x - e1y * e2y ));
	Qp[1].data[j] = omega_2_3 * sqrtoftwo * ( - sin(2 * phi) * ( e1x * e1y -
          e2x * e2y ) + cos(2 * phi) * ( e1x * e2y + e1y * e2x ));
	Qp[2].data[j] = omega_2_3 * sqrtoftwo * ( - sin(2 * phi) * ( e1x * e1z -
          e2x * e2z ) + cos(2 * phi) * ( e1x * e2z + e1z * e2x ));
	Qp[3].data[j] = omega_2_3 * sqrtoftwo * ( - sin(2 * phi) * ( e1y * e1z -
          e2y * e2z ) + cos(2 * phi) * ( e1y * e2z + e1z * e2y ));
	Qp[4].data[j] = omega_2_3 * onebysqrtofsix * ( - sin(2 * phi) *
        ( 2 * e2z * e2z - 2 * e1z * e1z + e1x * e1x + e1y * e1y -
          e2x * e2x - e2y * e2y ) + 2 * cos(2 * phi) * ( e1x * e2x +
            e1y * e2y - 2 * e1z * e2z ));
  }

  /* evaluate the P_I factors from the extrinsic parameters */
  P[0] = -(1 + cos(Theta) * cos(Theta)) * cos(2 * Phi) * cos(Psi)
		 + 2 * cos(Theta) * sin(2 * Phi) * sin(Psi);
  P[1] = - 2 * cos(Theta) * cos(2 * Phi) * sin(Psi)
		 -(1 + cos(Theta) * cos(Theta)) * sin(2 * Phi) * cos(Psi);
  P[2] = 2 * sin(Theta) * (-sin(Phi) * sin(Psi) + cos(Theta) * cos(Phi) * cos(Psi));
  P[3] = 2 * sin(Theta) * ( cos(Phi) * sin(Psi) + cos(Theta) * sin(Phi) * cos(Psi));
  P[4] = 3 * sin(Theta) * sin(Theta) * cos(Psi);

  /* Build the time domain PTF waveform */
  /* fprintf( stdout, "initial phase: %f\n", phi0); */
  for ( j = 0; j < N; ++j )
  {
    ptfwave->data[j] = P[0] * ( Q[0].data[j] * cos(phi0) + Qp[0].data[j] * sin(phi0) )
					 + P[1] * ( Q[1].data[j] * cos(phi0) + Qp[1].data[j] * sin(phi0) )
					 + P[2] * ( Q[2].data[j] * cos(phi0) + Qp[2].data[j] * sin(phi0) )
					 + P[3] * ( Q[3].data[j] * cos(phi0) + Qp[3].data[j] * sin(phi0) )
					 + P[4] * ( Q[4].data[j] * cos(phi0) + Qp[4].data[j] * sin(phi0) );
  }

  /* this memory deallocation code should be moved to a separate function */
  XLALDestroyREAL8VectorSequence( PTFQ );
  XLALDestroyREAL8VectorSequence( PTFQp );
  XLALDestroyREAL4Vector( PTFphi );
  XLALDestroyREAL4Vector( PTFomega_2_3 );
  XLALDestroyREAL4VectorSequence( PTFe1 );
  XLALDestroyREAL4VectorSequence( PTFe2 );
  XLALDestroyREAL8FFTPlan( fwdPlan );

  /* normal exit */
  return XLAL_SUCCESS;

}


INT4 XLALInspiralComputePTFWDeriv (
	COMPLEX16Vector			   *Wderiv,
	REAL8FrequencySeries       *psd,
    InspiralTemplate           *params,
	INT4					   paramid,
	REAL8					   initdelta,
	REAL8					   tolerance
    )

{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;

  /* number of points in a time-domain segment */
  UINT4 N = (Wderiv->length - 1) * 2;
  UINT4 j, k;

  /* get a local copy of the intrinstic parameters */
  REAL8 totalMass = params->totalMass;
  REAL8 eta = params->eta;
  REAL8 chi = params->chi;
  REAL8 kappa = params->kappa;

  /* bounds on the power spectrum integration for the moments */
  REAL8 deltaT = params->tSampling;

  /* deviated template parameters */
  InspiralTemplate pparams = *params;
  InspiralTemplate mparams = *params;

  /* Waveform length recorder */
  UINT4 clen, plen, mlen;
  INT4 dlen;

  /* local pointers to the waveforms data */
  /* this memory deallocation code should be moved to a separate function */
  REAL8Vector				*cwave;
  REAL8Vector				*pwavenew;
  REAL8Vector				*mwavenew;
  REAL8Vector				*pwaveold;
  REAL8Vector				*mwaveold;
  REAL8Vector				*pwaveold2;
  REAL8Vector				*mwaveold2;
  REAL8Vector				*wavederivnew;
  REAL8Vector				*wavederivold;
  REAL8Vector				*wavederivdiff;
  REAL8Vector				*wavederiv_2;
  REAL8Vector				*wavederivdiff_2;
  COMPLEX16Vector			*derivtilde;
  COMPLEX16Vector			*derivdifftilde;
  COMPLEX16Vector			*derivpowertilde;
  COMPLEX16Vector			*derivdiffpowertilde;
  REAL8FFTPlan				*fwdPlan;
  REAL8FFTPlan				*revPlan;
  UINT4 iter = 1;
  UINT4 maxiter = 20;
  REAL8 reldelta = initdelta;
  REAL8 absdelta = 0;
  REAL8 relderivdiff = 1.0;
  REAL8 invpsd;
  REAL8 powerderiv;
  REAL8 powerderivdiff;

  LALStatus status;
  memset( &status, 0, sizeof(LALStatus) );
  pparams.massChoice = totalMassAndEta;
  mparams.massChoice = totalMassAndEta;

  cwave = XLALCreateREAL8Vector( N );
  pwavenew = XLALCreateREAL8Vector( N );
  mwavenew = XLALCreateREAL8Vector( N );
  pwaveold = XLALCreateREAL8Vector( N );
  mwaveold = XLALCreateREAL8Vector( N );
  pwaveold2 = XLALCreateREAL8Vector( N );
  mwaveold2 = XLALCreateREAL8Vector( N );
  wavederivnew = XLALCreateREAL8Vector( N );
  wavederivold = XLALCreateREAL8Vector( N );
  wavederivdiff = XLALCreateREAL8Vector( N );
  wavederiv_2 = XLALCreateREAL8Vector( N );
  wavederivdiff_2 = XLALCreateREAL8Vector( N );
  derivtilde = XLALCreateCOMPLEX16Vector( N / 2 + 1 );
  derivdifftilde = XLALCreateCOMPLEX16Vector( N / 2 + 1 );
  derivpowertilde = XLALCreateCOMPLEX16Vector( N / 2 + 1 );
  derivdiffpowertilde = XLALCreateCOMPLEX16Vector( N / 2 + 1 );
  fwdPlan = XLALCreateForwardREAL8FFTPlan( N, 0 );
  revPlan = XLALCreateReverseREAL8FFTPlan( N, 0 );

  /* main loop */

  /* Generate central waveform */
  errcode = XLALInspiralComputePTFWaveform(cwave, params);
  if ( errcode != XLAL_SUCCESS )
  {
	fprintf( stderr, "XLALInspiralComputePTFWaveform failed\n" );
	exit( 1 );
  }
  clen = params->tC / deltaT + 0.5;

  /* Finite differencing */
  while (relderivdiff > tolerance && iter < maxiter + 1)
  {
	switch (paramid)
	{
	  case 1:
		absdelta = totalMass * reldelta;
		pparams.totalMass = totalMass + absdelta;
		mparams.totalMass = totalMass - absdelta;
		break;
	  case 2:
		absdelta = params->eta * reldelta;
		pparams.eta		  = eta		  + absdelta;
		mparams.eta		  = eta		  - absdelta;
		break;
	  case 3:
		absdelta = reldelta;
		pparams.chi		  = chi		  + absdelta;
		mparams.chi		  = chi		  - absdelta;
		break;
	  case 4:
		absdelta = reldelta;
		pparams.kappa	  = kappa	  + absdelta;
		mparams.kappa	  = kappa	  - absdelta;
		break;
	  default:
	    fprintf( stdout, "XLALInspiralComputePTFWDeriv failed\n" );
		break;
	}
	/* Update masses according to totalMass and eta */
	LALInspiralParameterCalc(&status, &pparams);
	if ( status.statusCode )
	{
      REPORTSTATUS( &status );
      exit( 1 );
	}
	LALInspiralParameterCalc(&status, &mparams);
	if ( status.statusCode )
	{
      REPORTSTATUS( &status );
      exit( 1 );
	}

	/* Generate new waveforms */
	errcode = XLALInspiralComputePTFWaveform(pwavenew, &pparams);
	errcode = XLALInspiralComputePTFWaveform(mwavenew, &mparams);
	if ( errcode != XLAL_SUCCESS )
	{
	  fprintf( stderr, "XLALInspiralComputePTFWaveform failed\n" );
      exit( 1 );
	}
	plen = pparams.tC / deltaT + 0.5;
	mlen = mparams.tC / deltaT + 0.5;

	/* Shift deviated waveforms to start the same time as the center waveform */
	dlen = plen - clen;
	if (dlen > 0)
	{
	  memmove(pwavenew->data + dlen, pwavenew->data, (N - dlen) * sizeof(REAL8));
	  memset(pwavenew->data, 0, dlen * sizeof(REAL8));
	}
	else if (dlen < 0)
	{
	  memmove(pwavenew->data, pwavenew->data - dlen, (N + dlen) * sizeof(REAL8));
	  memset(pwavenew->data + (N + dlen), 0, (- dlen) * sizeof(REAL8));
	}
	dlen = mlen - clen;
	if (dlen > 0)
	{
	  memmove(mwavenew->data + dlen, mwavenew->data, (N - dlen) * sizeof(REAL8));
	  memset(mwavenew->data, 0, dlen * sizeof(REAL8));
	}
	else if (dlen < 0)
	{
	  memmove(mwavenew->data, mwavenew->data - dlen, (N + dlen) * sizeof(REAL8));
	  memset(mwavenew->data + (N + dlen), 0, (- dlen) * sizeof(REAL8));
	}

	/* Calculate waveform derivative using finite differencing formula */
	if (iter == 1)
	{
	  for (j = 0; j < N; ++j)
	  {
	    wavederivnew->data[j] =
			(pwavenew->data[j]/2 - mwavenew->data[j]/2) / absdelta;
		wavederivdiff->data[j] = wavederivnew->data[j];
		wavederivold->data[j] = wavederivnew->data[j];
		pwaveold2->data[j] = pwavenew->data[j];
		mwaveold2->data[j] = mwavenew->data[j];
		pwaveold->data[j] = pwavenew->data[j];
		mwaveold->data[j] = mwavenew->data[j];
	  }
	}
	else /* if (iter == 2) */
	{
	  for (j = 0; j < N; ++j)
	  {
	    wavederivnew->data[j] =
			(- pwaveold->data[j] / 12 + 2 * pwavenew->data[j] / 3
			 - 2 * mwavenew->data[j] / 3 + mwaveold->data[j] / 12) / absdelta;
		wavederivdiff->data[j] = wavederivnew->data[j] - wavederivold->data[j];
		wavederivold->data[j] = wavederivnew->data[j];
		pwaveold2->data[j] = pwaveold->data[j];
		mwaveold2->data[j] = mwaveold->data[j];
		pwaveold->data[j] = pwavenew->data[j];
		mwaveold->data[j] = mwavenew->data[j];
	  }
	}
	/* else
	{
	  for (j = 0; j < N; ++j)
	  {
	    wavederivnew->data[j] =
			(pwaveold2->data[j] / 60 - pwaveold->data[j] * 3 / 20 + pwavenew->data[j] * 3 / 4
			 - mwavenew->data[j] * 3 / 4 + mwaveold->data[j] * 3 / 20 - mwaveold2->data[j] / 60) / absdelta;
		wavederivdiff->data[j] = wavederivnew->data[j] - wavederivold->data[j];
		wavederivold->data[j] = wavederivnew->data[j];
		pwaveold2->data[j] = pwaveold->data[j];
		mwaveold2->data[j] = mwaveold->data[j];
		pwaveold->data[j] = pwavenew->data[j];
		mwaveold->data[j] = mwavenew->data[j];
	  }
	} */

	/* Fourier transform of derivative and derivative difference of waveforms */
	XLALREAL8ForwardFFT( derivtilde, wavederivold, fwdPlan);
    XLALREAL8ForwardFFT( derivdifftilde, wavederivdiff, fwdPlan);

	for (k = 0; k < N/2 + 1; ++k)
	{
	  if (psd->data->data[k] == 0.)
	    invpsd = 0.;
	  else
		invpsd = 1.0 / psd->data->data[k];
	  derivpowertilde->data[k] = crect( (creal(derivtilde->data[k]) * creal(derivtilde->data[k]) +cimag(derivtilde->data[k]) * cimag(derivtilde->data[k])) * invpsd, 0 );
	  derivdiffpowertilde->data[k] = crect( (creal(derivdifftilde->data[k]) * creal(derivdifftilde->data[k]) +cimag(derivdifftilde->data[k]) * cimag(derivdifftilde->data[k])) * invpsd, 0 );
	}

	/* Inverse Fourier of derivdifftilde */
    XLALREAL8ReverseFFT( wavederiv_2, derivpowertilde, revPlan);
    XLALREAL8ReverseFFT( wavederivdiff_2, derivdiffpowertilde, revPlan);
	/* I will take the power for now */
	powerderiv = 4.0 * wavederiv_2->data[0];
	powerderivdiff = 4.0 * wavederivdiff_2->data[0];

	relderivdiff = powerderivdiff / powerderiv;

	if (relderivdiff < tolerance)
	{
	  for (k = 0; k < N / 2 + 1; ++k)
	  {
	    Wderiv->data[k] = derivtilde->data[k];
	  }
	  break;
	}

	reldelta *= 0.5;
	iter += 1;
  }

  /* this memory deallocation code should be moved to a separate function */
  XLALDestroyREAL8Vector( cwave );
  XLALDestroyREAL8Vector( pwavenew );
  XLALDestroyREAL8Vector( mwavenew );
  XLALDestroyREAL8Vector( pwaveold );
  XLALDestroyREAL8Vector( mwaveold );
  XLALDestroyREAL8Vector( pwaveold2 );
  XLALDestroyREAL8Vector( mwaveold2 );
  XLALDestroyREAL8Vector( wavederivold );
  XLALDestroyREAL8Vector( wavederivnew );
  XLALDestroyREAL8Vector( wavederivdiff );
  XLALDestroyREAL8Vector( wavederiv_2 );
  XLALDestroyREAL8Vector( wavederivdiff_2 );
  XLALDestroyCOMPLEX16Vector( derivtilde );
  XLALDestroyCOMPLEX16Vector( derivdifftilde );
  XLALDestroyCOMPLEX16Vector( derivpowertilde );
  XLALDestroyCOMPLEX16Vector( derivdiffpowertilde );
  XLALDestroyREAL8FFTPlan( fwdPlan );
  XLALDestroyREAL8FFTPlan( revPlan );

  /* normal exit */
  return XLAL_SUCCESS;
}


INT4 XLALInspiralComputePTFQDeriv (
    REAL8VectorSequence        UNUSED *Qderiv,
    InspiralTemplate           UNUSED *params
    )

{
  return XLAL_SUCCESS;
}
