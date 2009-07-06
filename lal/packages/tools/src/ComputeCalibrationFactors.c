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

/* <lalLaTeX>
\subsection{Module \texttt{ComputeCalibrationFactors.c}}
Independently computes the calibration factors \alpha(t)*\beta(t) and \alpha (t)
from .

\subsubsection*{Description}


\subsubsection*{Algorithm}


\subsubsection*{Uses}
\begin{verbatim}
None
\end{verbatim}

\subsubsection*{Notes}

</lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALConstants.h>
#include <lal/Calibration.h>

NRCSID( UPDATEFACTORSC, "$Id$" );
RCSID("$Id$");

static COMPLEX16 *cmul( COMPLEX16 *pc, COMPLEX16 *pa, COMPLEX16 *pb )
{
  COMPLEX16 a = *pa;
  COMPLEX16 b = *pb;
  COMPLEX16 c;
  c.re = a.re * b.re - a.im * b.im;
  c.im = a.re * b.im + a.im * b.re;
  *pc = c;
  return pc;
}

static COMPLEX16 *cdiv( COMPLEX16 *pc, COMPLEX16 *pa, COMPLEX16 *pb )
{
  COMPLEX16 a = *pa;
  COMPLEX16 b = *pb;
  COMPLEX16 c;
  REAL8 rat;
  REAL8 den;
  if ( fabs( b.re ) > fabs( b.im ) )
  {
    rat = b.im / b.re;
    den = b.re + rat * b.im;
    c.re = ( a.re + rat * a.im ) / den;
    c.im = ( a.im - rat * a.re ) / den;
  }
  else
  {
    rat = b.re / b.im;
    den = b.im + rat * b.re;
    c.re = ( a.re * rat + a.im ) / den;
    c.im = ( a.im * rat - a.re ) / den;
  }
  *pc = c;
  return pc;
}

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

  INITSTATUS( status, "LALComputeCalibrationFactors", UPDATEFACTORSC );

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
  DARM_CTRL.re = DARM_CTRL.im = 0;
  a = cos( 2.0 * LAL_PI * f0 * params->darmCtrl->deltaT );
  b = sin( 2.0 * LAL_PI * f0 * params->darmCtrl->deltaT );
  c = 1;
  s = 0;
  for ( i = 0; i < params->darmCtrl->data->length; ++i )
  {
    REAL8 tmp = a * c - b * s;
    DARM_CTRL.re += c * params->darmCtrl->data->data[i];
    DARM_CTRL.im += s * params->darmCtrl->data->data[i];
    s = b * c + a * s;
    c = tmp;
  }

  DARM_CTRL.re *= params->darmCtrl->deltaT;
  DARM_CTRL.im *= -params->darmCtrl->deltaT;  /* The negative sign is needed for correct FT convention */

  /* De-whiten DARM_CTRL */
  if(W0.re != 0 && W0.im != 0)
    {
     cmul( &DARM_CTRL, &DARM_CTRL, &W0); 
    }

  output->darm=DARM_CTRL;

  /* filter AS_Q against sinusoids at f0 */
  AS_Q.re = AS_Q.im = 0;
  a = cos( 2.0 * LAL_PI * f0 * params->asQ->deltaT );
  b = sin( 2.0 * LAL_PI * f0 * params->asQ->deltaT );
  c = 1;
  s = 0;
  for ( i = 0; i < params->asQ->data->length; ++i )
  {
    REAL8 tmp = a * c - b * s;
    AS_Q.re += c * params->asQ->data->data[i];
    AS_Q.im += s * params->asQ->data->data[i];
    s = b * c + a * s;
    c = tmp;
  }
  
  AS_Q.re *= params->asQ->deltaT;
  AS_Q.im *= -params->asQ->deltaT;

  output->asq=AS_Q;
 
  /* filter EXC against sinusoids at f0 */
  EXC.re = EXC.im = 0;
  a = cos( 2.0 * LAL_PI * f0 * params->exc->deltaT );
  b = sin( 2.0 * LAL_PI * f0 * params->exc->deltaT );
  c = 1;
  s = 0;
  for ( i = 0; i < params->exc->data->length; ++i )
  {
    REAL8 tmp = a * c - b * s;
    EXC.re += c * params->exc->data->data[i];
    EXC.im += s * params->exc->data->data[i];
    s = b * c + a * s;
    c = tmp;
  }

  EXC.re *= params->exc->deltaT;
  EXC.im *= -params->exc->deltaT;

  output->exc=EXC;

  if (( fabs( EXC.re ) < tiny && fabs( EXC.im ) < tiny )) /* check on DARM_CTRL too?? */
  {
    output->alphabeta.re=0.0;
    output->alphabeta.im=0.0;
    output->alpha.re=0.0;
    output->alpha.im=0.0;
    RETURN(status);
  }


  /* ------------------------------ compute alpha*beta ------------------------------------- */

  {
    COMPLEX16 RD, RDminus1;

    cdiv( &RD, &EXC, &DARM_CTRL);
    RDminus1.re=RD.re-1.0;
    RDminus1.im=RD.im;
    cdiv( &alphabeta, &RDminus1, &G0 );
  }

  /* ------------------------------- compute alpha ----------------------------------------- */
  {
    COMPLEX16 RQ;

    cdiv( &RQ, &AS_Q, &DARM_CTRL );
    cmul( &alpha, &RQ, &D0);
    cdiv( &alpha, &alpha, &G0 );
    alpha.re *= -1.0;
    alpha.im *= -1.0;
  }
  /* ------------------------------- compute beta ----------------------------------------- */

  {
    COMPLEX16 DmL,DmLoQ;

    DmL.re = DARM_CTRL.re - EXC.re;
    DmL.im = DARM_CTRL.im - EXC.im;
    cdiv( &DmLoQ, &DmL , &AS_Q );
    cdiv( &beta, &DmLoQ, &D0);
  }
  /* ------------------ Done, now put alpha, beta  and alpha*beta in the output structure --------- */

  output->alphabeta=alphabeta;
  output->alpha=alpha;
  output->beta=beta;

  RETURN( status );
}
