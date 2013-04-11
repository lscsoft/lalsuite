/*
*  Copyright (C) 2007 Ian Harry, Diego Fazi, Duncan Brown
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


#define LAL_USE_OLD_COMPLEX_STRUCTS
#include "coh_PTF.h"

void coh_PTF_template (
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *InspTmplt,
    FindChirpTmpltParams       *params
    )
{
  UINT4 i;
  LALStatus status = blank_status;
  switch ( params->approximant )
  {
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorT4:
    case GeneratePPN:
    case PadeT1:
    case EOB:
    case EOBNR:
    case IMRPhenomB:
      LALFindChirpTDTemplate( &status,fcTmplt,InspTmplt,params );
      break;
    case FindChirpSP:
      LALFindChirpSPTemplate( &status,fcTmplt,InspTmplt,params );
      for (i=0 ; i < params->xfacVec->length ; i++ )
      {
        fcTmplt->data->data[i].re = fcTmplt->data->data[i].re * params->PTFphi->data[i];
        fcTmplt->data->data[i].im = fcTmplt->data->data[i].im * params->PTFphi->data[i];
      }
      break;
    case FindChirpPTF:
      coh_PTF_template_PTF(fcTmplt,InspTmplt,params);
      break;
    default:
      fprintf(stderr,"Waveform approximant not recognized at template generation\n");
      exit(1);
  }
}

void
coh_PTF_template_PTF (
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *InspTmplt,
    FindChirpTmpltParams       *params
    )
{

  // This function generates Q_{1-5} as a function of time and then FFTs
  // them before returning Q(f) in fcTmplt->PTFQtilde
  // It would be nice to use the stuff in LAL to do most of this ....

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
  sanity_check( fcTmplt->PTFQ );
  sanity_check( fcTmplt->PTFQ->length == 5 );
  sanity_check( fcTmplt->PTFQ->data );

  /* check that the parameter structure exists */
  sanity_check( params );
  sanity_check( params->fwdPlan );

  /* check that the timestep is positive */
  sanity_check( params->deltaT > 0 );

  /* check that the input exists */
  sanity_check( InspTmplt );

  N = params->PTFphi->length;

  /* set the parameter structure */
  /* to the correct waveform approximant       */
  InspTmplt->approximant = FindChirpPTF;
  sanity_check( InspTmplt->fLower );

  /* copy the template parameters to the finchirp template structure */
  memcpy( &(fcTmplt->tmplt), InspTmplt, sizeof(InspiralTemplate) );

  /* XXX delete this line if the low frequency cutoff XXX */
  /* XXX should be read from the template bank        XXX */
  fcTmplt->tmplt.fLower = params->fLow = InspTmplt->fLower;

  /* Zero out the Q and Qtilde vectors */
  memset( fcTmplt->PTFQ->data, 0, 5 * N * sizeof(REAL4) );
  memset( fcTmplt->PTFQtilde->data, 0, 5 * (N /2 + 1) * sizeof(COMPLEX8) );

 /* Point the dummy variables Q and Qtilde to the actual output structures */
  for ( i = 0; i < 5; ++i )
  {
    Q[i].length      = N;
    Qtilde[i].length = N / 2 + 1;
    Q[i].data        = fcTmplt->PTFQ->data + (i * N);
    Qtilde[i].data   = fcTmplt->PTFQtilde->data + (i * (N / 2 + 1)) ;
  }

  /* call the waveform generation function */
  errcode = XLALFindChirpPTFWaveform( params->PTFphi, params->PTFomega_2_3,
                                      params->PTFe1, params->PTFe2, InspTmplt,
                                      params->deltaT);
  sanity_check( errcode == XLAL_SUCCESS );

  if (InspTmplt->tC > params->maxTempLength )
  {
    fprintf(stderr,"Template generated is longer than max. Template must not ");
    fprintf(stderr,"be longer than this as it causes wrapping issues in the ");
    fprintf(stderr,"FFT. Template length is %lf \n",InspTmplt->tC);
    exit(1);
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
    XLALREAL4ForwardFFT( &Qtilde[i], &Q[i],
        params->fwdPlan);
  }

  /* FIXME: Q and Qtilde should be outputtable. Probably easier outside of
     this function as template number is available. */

  fcTmplt->tmplt.tC = InspTmplt->tC; /* length of template in seconds */
  fcTmplt->tmplt.fFinal = InspTmplt->fFinal; /* upper freq of template in Hz */
}


void coh_PTF_normalize(
    struct coh_PTF_params      *params,
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
{
  // This function calculates the various filters used to calculate SNR
  // It calculates (Q_0 | Q_0), (Q | s) and if necessary it will calculate
  // (Q_0 | Q_{\pi/2}) as well (for the spin checker code).


  UINT4         i, j, k, kmin, len, kmax,numPoints,vecLength;
  REAL8         f_min, deltaF,deltaT, fFinal, r, s, x, y, length;
  COMPLEX8     *qtilde, *inputData;
  COMPLEX8Vector *qtildeVec,qVec;
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
  f_min     = params->lowFilterFrequency;
  kmin      = f_min / deltaF > 1 ?  f_min / deltaF : 1;
  fFinal    = params->highFilterFrequency;
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

//  sanity_check ( fcTmplt->tmplt.approximant == FindChirpPTF );

  /* Set parameters to determine spin/nonspin */
  /* For non-spin we only need one filter. For PTF all 5 are needed */
  vecLength = 1;
  if (spinTemplate == 1)
  {
    vecLength = 5;
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

  /* FIXME: We would like to be able to print off A,B and M. 
     like above, this may be better in the main function */
  
  XLALDestroyCOMPLEX8Vector( qtildeVec );
}

void coh_PTF_template_overlaps(
    struct coh_PTF_params      *params,
    FindChirpTemplate          *fcTmplt1,
    FindChirpTemplate          *fcTmplt2,
    REAL4FrequencySeries       *invspec,
    UINT4                      spinBank,
    REAL8Array                 *PTFM
    )
{
  // This function calculates the real part of the overlap between two templates

  UINT4         i, j, k, kmin, kmax, len,vecLen;
  REAL8         f_min, deltaF, fFinal;
  COMPLEX8     *PTFQtilde1   = NULL;
  COMPLEX8     *PTFQtilde2   = NULL;


  PTFQtilde1 = fcTmplt1->PTFQtilde->data;
  PTFQtilde2 = fcTmplt2->PTFQtilde->data;

  len       = invspec->data->length;
  deltaF    = invspec->deltaF;
  f_min     = params->lowFilterFrequency;
  kmin      = f_min / deltaF > 1 ?  f_min / deltaF : 1;
  fFinal    = params->highFilterFrequency;
  kmax      = fFinal / deltaF < (len - 1) ? fFinal / deltaF : (len - 1);

  vecLen = 5;
  if (! spinBank )
  {
    vecLen = 1;
  }

  for( i = 0; i < vecLen; ++i )
  {
    for ( j = 0; j < vecLen; ++j )
    {
      for ( k = kmin; k < kmax ; ++k )
      {
        PTFM->data[vecLen * i + j] += (PTFQtilde1[k + i * len].re *
                            PTFQtilde2[k + j * len].re +
                            PTFQtilde1[k + i * len].im *
                            PTFQtilde2[k + j * len].im )
                            * invspec->data->data[k] ;
      
      }
      PTFM->data[vecLen * i + j] *= 4.0 * deltaF ;
    }
  }

}

void coh_PTF_complex_template_overlaps(
    struct coh_PTF_params      *params,
    FindChirpTemplate          *fcTmplt1,
    FindChirpTemplate          *fcTmplt2,
    REAL4FrequencySeries       *invspec,
    UINT4                      spinBank,
    COMPLEX8Array                 *PTFM
    )
{

  // This function calculates the complex overlap between two templates
  UINT4         i, j, k, kmin, kmax, len,vecLen;
  REAL8         f_min, deltaF, fFinal;
  COMPLEX8     *PTFQtilde1   = NULL;
  COMPLEX8     *PTFQtilde2   = NULL;


  PTFQtilde1 = fcTmplt1->PTFQtilde->data;
  PTFQtilde2 = fcTmplt2->PTFQtilde->data;

  len       = invspec->data->length;
  deltaF    = invspec->deltaF;
  f_min     = params->lowFilterFrequency;
  kmin      = f_min / deltaF > 1 ?  f_min / deltaF : 1;
  fFinal    = params->highFilterFrequency;
  kmax      = fFinal / deltaF < (len - 1) ? fFinal / deltaF : (len - 1);

  vecLen = 5;
  if (! spinBank )
  {
    vecLen = 1;
  }

  for( i = 0; i < vecLen; ++i )
  {
    for ( j = 0; j < vecLen; ++j )
    {
      for ( k = kmin; k < kmax ; ++k )
      {
        PTFM->data[vecLen * i + j].re += (PTFQtilde1[k + i * len].re *
                            PTFQtilde2[k + j * len].re +
                            PTFQtilde1[k + i * len].im *
                            PTFQtilde2[k + j * len].im )
                            * invspec->data->data[k] ;
        PTFM->data[vecLen * i + j].im += (-PTFQtilde1[k + i * len].re *
                            PTFQtilde2[k + j * len].im +
                            PTFQtilde1[k + i * len].im *
                            PTFQtilde2[k + j * len].re )
                            * invspec->data->data[k] ;
      }
      PTFM->data[vecLen * i + j].re *= 4.0 * deltaF ;
      PTFM->data[vecLen * i + j].im *= 4.0 * deltaF ;
    }
  }

}

void coh_PTF_bank_filters(
    struct coh_PTF_params      *params,
    FindChirpTemplate          *fcTmplt,
    UINT4                      spinBank,
    COMPLEX8FrequencySeries    *sgmnt,
    COMPLEX8FFTPlan            *invBankPlan,
    COMPLEX8VectorSequence     *PTFqVec,
    COMPLEX8VectorSequence     *PTFBankqVec,
    REAL8                      f_min,
    REAL8                      fFinal)
{
  /* This function calculates (Q|s) for the bank veto. It only returns the 
   * middle half of the time series with some buffer to allow for time shifts */

  /* FIXME: Can this function be merged with normalize?? */

  UINT4          i, j, k, kmin, len, kmax,vecLen;
  REAL8          deltaF, r, s, x, y;
  COMPLEX8       *inputData,*qtilde;
  COMPLEX8Vector *qtildeVec,qVec;
  COMPLEX8       *PTFQtilde   = NULL;

  len       = sgmnt->data->length;
  PTFQtilde = fcTmplt->PTFQtilde->data;
  deltaF    = sgmnt->deltaF;
/*  deltaT    = 1.0 / ( deltaF * (REAL4) numPoints); */

  /* F_min and F_max are used to do the chisquared limited filters */
  if (! f_min)
  {
    f_min     = params->lowFilterFrequency;
  }
  kmin      = f_min / deltaF > 1 ?  f_min / deltaF : 1;
  if (! fFinal)
  {
    fFinal    = params->highFilterFrequency;
  }
  kmax      = fFinal / deltaF < (len - 1) ? fFinal / deltaF : (len - 1);
  qVec.length = params->numTimePoints;
  qtildeVec    = XLALCreateCOMPLEX8Vector( params->numTimePoints );
  qtilde = qtildeVec->data;

  if (! spinBank )
  {
    vecLen = 1;
  }
  else
    vecLen = 5;

  /* Data params */
  inputData   = sgmnt->data->data;

  for ( i = 0; i < vecLen; ++i )
  {
    memset( qtildeVec->data, 0,
        qtildeVec->length * sizeof(COMPLEX8) );
    /* qtilde positive frequency, not DC or nyquist */
    for ( k = kmin; k < kmax ; ++k )
    {
      r = inputData[k].re;
      s = inputData[k].im;
      x = PTFQtilde[i * (params->numFreqPoints) + k].re;
      y = 0 - PTFQtilde[i * (params->numFreqPoints) + k].im; /* cplx conj */

      qtilde[k].re = 4. * (r*x - s*y)*deltaF;
      qtilde[k].im = 4. * (r*y + s*x)*deltaF;
    }

    qVec.data = PTFqVec->data + (i * params->numTimePoints);

    /* inverse fft to get q */
    XLALCOMPLEX8VectorFFT( &qVec, qtildeVec, invBankPlan );
  }
  XLALDestroyCOMPLEX8Vector( qtildeVec );

  for ( i = 0; i < vecLen ; i++ )
  {
    for ( j = params->analStartPointBuf; j < params->analEndPointBuf; ++j )
    {
      PTFBankqVec->data[i*(params->numAnalPointsBuf) + \
                        (j-params->analStartPointBuf)]\
               = PTFqVec->data[i*params->numTimePoints + j];
    }
  }
}

void coh_PTF_auto_veto_overlaps(
    struct coh_PTF_params      *params,
    FindChirpTemplate          *fcTmplt,
    struct bankComplexTemplateOverlaps *autoTempOverlaps,
    REAL4FrequencySeries       *invspec,
    COMPLEX8FFTPlan            *invBankPlan,
    UINT4                      spinBank,
    UINT4                      numAutoPoints, 
    UINT4                      timeStepPoints,
    UINT4                      ifoNumber )
{
  // This function calculate (Q | Q(delta_t) ) at various different points
  // for the auto veto

  UINT4          i, j, k, kmin, len, kmax,vecLen,numPoints;
  REAL8          f_min, deltaF, fFinal, r, s, x, y;
  COMPLEX8       *qtilde;
  COMPLEX8Vector *qtildeVec,*qVec;
  COMPLEX8       *PTFQtilde   = NULL;

  len         = params->numFreqPoints;
  PTFQtilde = fcTmplt->PTFQtilde->data;
  numPoints   = params->numTimePoints;
  deltaF      = invspec->deltaF;
//  deltaT    = 1.0 / ( deltaF * (REAL4) len);

  f_min     = params->lowFilterFrequency;
  kmin      = f_min / deltaF > 1 ?  f_min / deltaF : 1;
  fFinal    = params->highFilterFrequency;
  kmax      = fFinal / deltaF < (len - 1) ? fFinal / deltaF : (len - 1);
  qVec = XLALCreateCOMPLEX8Vector( numPoints );
  qtildeVec    = XLALCreateCOMPLEX8Vector( numPoints );
  qtilde = qtildeVec->data;

  if (! spinBank )
    vecLen = 1;
  else
    vecLen = 5;

  for ( i = 0; i < vecLen; ++i )
  {
    for ( j=0; j < vecLen; ++j )
    {
      memset( qtildeVec->data, 0,
          qtildeVec->length * sizeof(COMPLEX8) );
      /* qtilde positive frequency, not DC or nyquist */
      for ( k = kmin; k < kmax ; ++k )
      {
        r = PTFQtilde[i * (len ) + k].re;
        s = PTFQtilde[i * (len ) + k].im;
        x = PTFQtilde[j * (len ) + k].re;
        y = 0 - PTFQtilde[j * (len ) + k].im; /* cplx conj */

        qtilde[k].re = 4. * (r*x - s*y)*deltaF * invspec->data->data[k];
        qtilde[k].im = 4. * (r*y + s*x)*deltaF * invspec->data->data[k];
      }

      /* inverse fft to get q */
      XLALCOMPLEX8VectorFFT( qVec, qtildeVec, invBankPlan );
      for ( k = 1; k < numAutoPoints+1 ; k++ )
      {
        autoTempOverlaps[k-1].timeStepDelay = - k * timeStepPoints;
        autoTempOverlaps[k-1].PTFM[ifoNumber]->data[i*vecLen + j] = qVec->data[numPoints - k * timeStepPoints];
      }
    }
  }

  /* FIXME: We should be able to print off the autocorrelation timeseries */

  XLALDestroyCOMPLEX8Vector( qtildeVec);
  XLALDestroyCOMPLEX8Vector( qVec );
}
