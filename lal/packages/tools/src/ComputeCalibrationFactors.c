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
#include <lal/LALConstants.h>
#include <lal/Calibration.h>

NRCSID( UPDATEFACTORSC, "$Id$" );

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
  COMPLEX16 DARM_CTRL;
  COMPLEX16 AS_Q;
  COMPLEX16 EXC;
  COMPLEX16 RQ;
  COMPLEX16 RD,OneMinusRD;
  COMPLEX16 C0;
  COMPLEX16 A0;
  COMPLEX16 R0;
  COMPLEX16 H0;
  COMPLEX16 H;
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
  
  R0 = params->responseFactor;
  C0 = params->sensingFactor;
  A0 = params->actuationFactor;

  cmul(&H0,&R0,&C0);
  H0.re -= 1.0;        /* Open loop gain */


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

  if ( fabs( EXC.re ) < tiny && fabs( EXC.im ) < tiny )
  {
    output->alphabeta.re=0.0;
    output->alphabeta.im=0.0;
    output->alpha.re=0.0;
    output->alpha.im=0.0;
    RETURN(status);
  }

  /* Adjust excitation for the fact that it does not see the entire actuation function */

 
  DARM_CTRL.re -= EXC.re;
  DARM_CTRL.im -= EXC.im;
  

  EXC.re *= params->mu; 
  EXC.im *= params->mu; 
  
  /* ------------------------------ compute alpha*beta ------------------------------------- */

  cdiv( &RD, &DARM_CTRL, &EXC );

  OneMinusRD.re=1.0-RD.re;
  OneMinusRD.im=-RD.im;
  
  if ( fabs( OneMinusRD.re ) < tiny && fabs( OneMinusRD.im ) < tiny )
  {
    ABORT( status, CALIBRATIONH_EZERO, CALIBRATIONH_MSGEZERO );
  }

  cdiv( &alphabeta, &RD, &OneMinusRD );
  cdiv( &alphabeta, &alphabeta, &H0 );

  /* ------------------------------- compute alpha ----------------------------------------- */

  cdiv( &RQ, &AS_Q, &EXC );
  cdiv( &alpha, &RQ, &OneMinusRD);
  cdiv( &alpha, &alpha, &A0 );
  cdiv( &alpha, &alpha, &C0 );

  /* ------------------ Done, now put alpha and alpha*beta in the output structure --------- */

  output->alphabeta=alphabeta;
  output->alpha=alpha;

  RETURN( status );
}
