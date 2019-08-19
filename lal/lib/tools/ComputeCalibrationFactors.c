/*
*  Copyright (C) 2007 Jolien Creighton, Xavier Siemens
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

#include <complex.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/Calibration.h>

/* Independently computes the calibration factors \alpha(t)*\beta(t) and \alpha (t) from */


void LALComputeCalibrationFactors(
    LALStatus               *status,
    CalFactors              *output,
    UpdateFactorsParams    *params
    )
{
  const REAL8 tiny = LAL_REAL8_MIN;
  COMPLEX16 alpha;
  COMPLEX16 alphabeta;
  COMPLEX16 beta;
  COMPLEX16 DARM_CTRL;
  COMPLEX16 AS_Q;
  COMPLEX16 EXC;
  COMPLEX16 G0,D0,W0;
  REAL8 f0;
  REAL8 a;
  REAL8 b;
  REAL8 c;
  REAL8 s;
  UINT4 i;

  INITSTATUS(status);

  ASSERT( output, status, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params, status, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );

  /* sanity check on params */
  ASSERT( params->darmCtrl, status, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->asQ, status, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->exc, status, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->darmCtrl->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->asQ->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->exc->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( (int)params->darmCtrl->data->length > 0, status,
      CALIBRATIONH_ESIZE, CALIBRATIONH_MSGESIZE );
  ASSERT( (int)params->asQ->data->length > 0, status,
      CALIBRATIONH_ESIZE, CALIBRATIONH_MSGESIZE );
  ASSERT( (int)params->exc->data->length > 0, status,
      CALIBRATIONH_ESIZE, CALIBRATIONH_MSGESIZE );
  ASSERT( params->darmCtrl->data->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->asQ->data->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->exc->data->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );

  f0 = params->lineFrequency;

  G0 = params->openloop;
  D0 = params->digital;
  W0 = params->whitener;


  /* ------------------------------------- match filtering ----------------------------------- */

  /* filter DARM_CTRL against sinusoids at f0 */
  DARM_CTRL = 0;
  a = cos( 2.0 * LAL_PI * f0 * params->darmCtrl->deltaT );
  b = sin( 2.0 * LAL_PI * f0 * params->darmCtrl->deltaT );
  c = 1;
  s = 0;
  for ( i = 0; i < params->darmCtrl->data->length; ++i )
  {
    REAL8 tmp = a * c - b * s;
    DARM_CTRL += c * params->darmCtrl->data->data[i];
    DARM_CTRL += I * s * params->darmCtrl->data->data[i];
    s = b * c + a * s;
    c = tmp;
  }

  DARM_CTRL = -conj(DARM_CTRL) * params->darmCtrl->deltaT;  /* The conjugation is needed for correct FT convention */

  /* De-whiten DARM_CTRL */
  if(W0 != 0.0)
    {
     DARM_CTRL *= W0;
    }

  output->darm=DARM_CTRL;

  /* filter AS_Q against sinusoids at f0 */
  AS_Q = 0;
  a = cos( 2.0 * LAL_PI * f0 * params->asQ->deltaT );
  b = sin( 2.0 * LAL_PI * f0 * params->asQ->deltaT );
  c = 1;
  s = 0;
  for ( i = 0; i < params->asQ->data->length; ++i )
  {
    REAL8 tmp = a * c - b * s;
    AS_Q += c * params->asQ->data->data[i];
    AS_Q += I * s * params->asQ->data->data[i];
    s = b * c + a * s;
    c = tmp;
  }

  AS_Q = conj(AS_Q) * params->asQ->deltaT;

  output->asq=AS_Q;

  /* filter EXC against sinusoids at f0 */
  EXC = 0;
  a = cos( 2.0 * LAL_PI * f0 * params->exc->deltaT );
  b = sin( 2.0 * LAL_PI * f0 * params->exc->deltaT );
  c = 1;
  s = 0;
  for ( i = 0; i < params->exc->data->length; ++i )
  {
    REAL8 tmp = a * c - b * s;
    EXC += c * params->exc->data->data[i];
    EXC += I * s * params->exc->data->data[i];
    s = b * c + a * s;
    c = tmp;
  }

  EXC = conj(EXC) * params->exc->deltaT;

  output->exc=EXC;

  if (( fabs( creal(EXC) ) < tiny && fabs( cimag(EXC) ) < tiny )) /* check on DARM_CTRL too?? */
  {
    output->alphabeta=0.0;
    output->alpha=0.0;
    RETURN(status);
  }


  /* ------------------------------ compute alpha*beta ------------------------------------- */

  {
    COMPLEX16 RD;

    RD = EXC / DARM_CTRL;
    alphabeta = (RD - 1.0) / G0;
  }

  /* ------------------------------- compute alpha ----------------------------------------- */
  {
    COMPLEX16 RQ;

    RQ = AS_Q / DARM_CTRL;
    alpha = -RQ * D0 / G0;
  }
  /* ------------------------------- compute beta ----------------------------------------- */

  {
    COMPLEX16 DmLoQ;

    DmLoQ = (DARM_CTRL - EXC) / AS_Q;
    beta = DmLoQ / D0;
  }
  /* ------------------ Done, now put alpha, beta  and alpha*beta in the output structure --------- */

  output->alphabeta=alphabeta;
  output->alpha=alpha;
  output->beta=beta;

  RETURN( status );
}
