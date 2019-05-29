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
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A., and Fazi, D.
 * \ingroup FindChirpPTF_h
 * \file
 *
 * \brief Provides functions to create physical template family templates in a
 * form that can be used by the <tt>FindChirpPTFFilter()</tt> function.
 *
 * ### Prototypes ###
 *
 * The function <tt>LALFindChirpPTFTemplate()</tt> creates a physical template
 * family template as described by the algorithm below.
 *
 * ### Algorithm ###
 *
 * Blah.
 *
 */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpPTF.h>
#include <lal/MatrixUtils.h>

void
LALFindChirpPTFTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *InspTmplt,
    FindChirpTmpltParams       *params
    )

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

  INITSTATUS(status);
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
  ASSERT( params, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( fcTmplt->PTFQ, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( fcTmplt->PTFQ->length == 5, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( fcTmplt->PTFQ->data, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the parameter structure exists */

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
