/********************************* <lalVerbatim file="ComputeTransferCV">
Author: Patrick Brady, Jolien Creighton
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{ComputeTransfer.c}}
\label{ss:ComputeTransfer.c}

Computes the transfer function from zero-pole-gain representation.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ComputeTransferCP}
\idx{LALComputeTransfer()}
\idx{LALUpdateCalibration()}
\idx{LALResponseConvert()}

\subsubsection*{Description}

A transfer function can either be specified as a list of coefficients or a 
list of poles and zeros. The function \verb@LALComputeTransfer()@ computes the 
frequency representation of the transfer function \verb+calrec->transfer+
described by the zeroes,  poles and gain in \verb@*calrec@.   The memory for 
the frequency series should be allocated before calling this routine which uses
\verb+calrec->transfer->deltaF+ and \verb+calrec->transfer->data->npoints+.

The routine \texttt{LALUpdateCalibration()} updates the response function
and the sensing function from some reference functions to the current
functions using information about the calibration lines.  The two calibration
lines yield two constants (as a slowly-varying function of time) that are
used as coefficients to the reference response and sensing functions to
compute the current response and sensing functions.  These coefficients are
stored in time series in the parameter structure, along with the current
epoch and duration for which the calibration functions are to be computed.  If 
the duration is zero, the calibration factors are the first ones at or after 
the given epoch.  If the duration is non-zero, then the calibration factors are
an average of all calibrations between epoch and epoch + duration.

The routine \texttt{LALResponseConvert()} takes a given frequency series
and converts it to a new frequency series by performing the following steps:
(i) the original frequency series is interpolated (using linear interpolation
of the real and imaginary parts independently) to the frequencies required
in the output frequency series;  (ii) if the output frequency series has units
that are the inverse of those of the input frequency series, the data is
inverted;  (iii) the data is scaled by an appropriate power of ten computed
from the input units and the output units.  For example you can convert from 
strain per count to counts per atto-strain.

\subsubsection*{Algorithm}

The transfer function is deduced from the poles and zeros as follows:
\begin{equation}
T(f) = c_{\mathrm{gain}} 
{\prod_i \textrm{zero}(f,z_i)}{ \prod_{i} \textrm{pole}(f,p_i)}
\end{equation}
where 
\begin{equation}
\textrm{zero}(f,z) = \left\{ \begin{array}{ll}
(i f / z) + 1 & \textrm{ when } z\neq 0 \\
i f & \textrm{ when } z = 0
\end{array}
\right.
\end{equation}
and 
\begin{equation}
\textrm{pole}(f,p) = \left\{ \begin{array}{ll}
\frac{1}{(i f / p) + 1} & \textrm{ when } p \neq 0 \\
\frac{1}{i f} & \textrm{ when } p = 0
\end{array}
\right.
\end{equation}
For both the transfer function and the pole-zero notation the units for
frequency is Hz rather than rad/s (angular frequency).  In particular, poles
and zeros are specified by their location in frequency space

To update the response function from one epoch to another, two functions are
needed.  These are the sensing function $C(f)$ and the response function
$R(f)$, which are related by
\begin{equation}
  R(f) = \frac{1+H(f)}{C(f)}
\end{equation}
where $H(f)$ is the open-loop gain function.  If the sensing function and
the open-loop gain function are known at some reference time ($C_0(f)$ and
$H_0(f)$) then the sensing function and open-loop gain function can be
calculated at a later time.  They are $C(f)=\alpha C_0(f)$ and
$H(f)=\alpha\beta H_0(f)$ where $\alpha$ and $\beta$ are slowly varying
coefficients that account for overall changes in the gains of the sensing
function and the open-loop gain.  The coefficients $\alpha$ and $\alpha\beta$
can be determined, as slowly-varying functions of time, by monitoring the
two injected calibration lines.  Thus, an updated sensing function and response
function can be computed from reference sensing function and response function,
$C_0(f)$ and $R_0(f)$ via:
\begin{equation}
  C(f) = \alpha C_0(f)
\end{equation}
and
\begin{equation}
  R(f) = \frac{1+\alpha\beta[C_0(f)R_0(f)-1]}{\alpha C_0(f)}
\end{equation}
where $\alpha$ and $\beta$ are those values of the coefficients that are
appropriate for the particular epoch.

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}
The DC component of \verb+calrec->transfer+ is always filled with $1 + i 0$.  
In most cases,  this should be irrelevant for gravitational wave data analysis,  
but care should be taken if DC is relevant when this function is used.  

\vfill{\footnotesize\input{ComputeTransferCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Calibration.h>
#include <lal/Units.h>
#include <lal/Date.h>

#define CAL_S2START 729273613
#define CAL_S2END 734367613

NRCSID( COMPUTETRANSFERC, "$Id$" );

static void product(COMPLEX8 *c,COMPLEX8 *a, COMPLEX8 *b) {
  
  c->re = a->re * b->re - a->im * b->im;
  c->im = a->re * b->im + a->im * b->re;
  
  return;
}

static void ratio(COMPLEX8 *c,COMPLEX8 *a, COMPLEX8 *b) {
  REAL4 norm;

  norm = b->re * b->re + b->im * b->im;
  
  c->re = (a->re * b->re + a->im * b->im)/norm;
  c->im = (- a->re * b->im + a->im * b->re)/norm;
  
  return;
}

/* <lalVerbatim file="ComputeTransferCP"> */
void
LALComputeTransfer( LALStatus                 *stat,
                    CalibrationRecord         *calrec
                    )
/* </lalVerbatim> */
{ 
  UINT4         i, j;                    /* indexes               */
  UINT4         jmin;                    /* used to handle DC     */
  COMPLEX8      dummy, dummyC, factor;   /* dummy complex numbers */
  REAL4         f,df;                    /* freq and interval     */
  REAL8         norm;

  INITSTATUS( stat, "LALComputeTransfer", COMPUTETRANSFERC );
  ATTATCHSTATUSPTR (stat);

  /* Make sure parameter structures and their fields exist. */
  ASSERT( calrec, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( calrec->zeros, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( calrec->poles, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );  
  ASSERT( calrec->transfer, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );  

  /* initialise everything */ 
  dummyC.re = factor.re = 1.0;
  dummyC.im = factor.im = 0.0;
  df = calrec->transfer->deltaF;
  jmin = 0;

  /* compute the normalization constant */
  norm = calrec->gain;
  
  for ( i=0 ; i < calrec->poles->length ; i++)
    if ( calrec->poles->data[i] != 0 )
    {
      norm *= calrec->poles->data[i];
    }

  for ( i=0 ; i < calrec->zeros->length ; i++)
    if ( calrec->zeros->data[i] != 0 )
    {
      norm /= calrec->zeros->data[i];
    }

  /* Handle DC if necessary */
  if ( calrec->transfer->f0 == 0.0 )
  {
    calrec->transfer->data->data[0].re = 1.0;
    calrec->transfer->data->data[0].im = 0.0;
    jmin = 1;
  }

  /* loop over frequency in the output */
  for ( j=jmin ; j<calrec->transfer->data->length ; j++)
  {
    /* the frequency */
    f = calrec->transfer->f0 + (REAL4) j * df;
    dummyC.re=1.0;
    dummyC.im=0.0;

    /* loop over zeroes */
    for (i = 0 ; i < calrec->zeros->length ; i++)
    {
      factor.re = (calrec->zeros->data[i]);
      factor.im = f;
      product( &dummy, &dummyC, &factor);
      dummyC.re=dummy.re;
      dummyC.im=dummy.im;
    }

    /* loop over poles */
    for (i = 0 ; i < calrec->poles->length ; i++)
    {
      factor.re = (calrec->poles->data[i]);
      factor.im = f;
      ratio( &dummy, &dummyC, &factor);
      dummyC.re=dummy.re;
      dummyC.im=dummy.im;
    }

    /* fill the frequency series */
    calrec->transfer->data->data[j].re = norm * dummyC.re;
    calrec->transfer->data->data[j].im = norm * dummyC.im;
  }

 
  /* we're out of here */
  DETATCHSTATUSPTR (stat);
  RETURN( stat );
}





#define cini COMPLEX8 tmpa, tmpb, tmpc; REAL4 tmpx, tmpy

#define cmul( a, b ) \
( tmpa = (a), tmpb = (b), \
  tmpc.re = tmpa.re * tmpb.re - tmpa.im * tmpb.im, \
  tmpc.im = tmpa.re * tmpb.im + tmpa.im * tmpb.re, \
  tmpc )

#define cdiv( a, b ) \
( tmpa = (a), tmpb = (b), \
  fabs( tmpb.re ) >= fabs( tmpb.im ) ? \
    ( tmpx = tmpb.im / tmpb.re, \
      tmpy = tmpb.re + tmpx * tmpb.im, \
      tmpc.re = ( tmpa.re + tmpx * tmpa.im ) / tmpy, \
      tmpc.im = ( tmpa.im - tmpx * tmpa.re ) / tmpy, \
      tmpc ) : \
    ( tmpx = tmpb.re / tmpb.im, \
      tmpy = tmpb.im + tmpx * tmpb.re, \
      tmpc.re = ( tmpa.re * tmpx + tmpa.im ) / tmpy, \
      tmpc.im = ( tmpa.im * tmpx - tmpa.re ) / tmpy, \
      tmpc ) )

#define cinv( b ) \
( tmpb = (b), \
  fabs( tmpb.re ) >= fabs( tmpb.im ) ? \
    ( tmpx = tmpb.im / tmpb.re, \
      tmpy = tmpb.re + tmpx * tmpb.im, \
      tmpc.re = 1 / tmpy, \
      tmpc.im = -tmpx / tmpy, \
      tmpc ) : \
    ( tmpx = tmpb.re / tmpb.im, \
      tmpy = tmpb.im + tmpx * tmpb.re, \
      tmpc.re = tmpx / tmpy, \
      tmpc.im = -1 / tmpy, \
      tmpc ) )


/* <lalVerbatim file="ComputeTransferCP"> */
void
LALUpdateCalibration(
    LALStatus               *status,
    CalibrationFunctions    *output,
    CalibrationFunctions    *input,
    CalibrationUpdateParams *params
    )
{ /* </lalVerbatim> */
  const REAL4 tiny = 1e-6;
  cini;
  COMPLEX8Vector *save;
  COMPLEX8 *R;
  COMPLEX8 *C;
  COMPLEX8 *R0;
  COMPLEX8 *C0;
  COMPLEX8 a;
  COMPLEX8 ab;
  REAL8 epoch;
  REAL8 first_cal;
  REAL8 duration;
  REAL4 dt;
  REAL4 i_r4;
  UINT4 n;
  UINT4 i;
  UINT4 length = 0;
  CHAR  warnMsg[512];

  INITSTATUS( status, "LALUpdateCalibration", COMPUTETRANSFERC );
  ATTATCHSTATUSPTR( status );

  /* check input */
  ASSERT( input, status, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( input->sensingFunction, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( input->responseFunction, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( input->sensingFunction->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( input->responseFunction->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( input->sensingFunction->data->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( input->responseFunction->data->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  n = input->sensingFunction->data->length;
  ASSERT( (int)n > 0, status, CALIBRATIONH_ESIZE, CALIBRATIONH_MSGESIZE );
  ASSERT( input->responseFunction->data->length == n, status,
      CALIBRATIONH_ESZMM, CALIBRATIONH_MSGESZMM );

  /* check output */
  ASSERT( output, status, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( output->sensingFunction, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( output->responseFunction, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( output->sensingFunction->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( output->responseFunction->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( output->sensingFunction->data->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( output->responseFunction->data->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( output->sensingFunction->data->length == n, status,
      CALIBRATIONH_ESZMM, CALIBRATIONH_MSGESZMM );
  ASSERT( output->responseFunction->data->length == n, status,
      CALIBRATIONH_ESZMM, CALIBRATIONH_MSGESZMM );

  /* check params */
  ASSERT( params, status, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->sensingFactor, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->openLoopFactor, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->sensingFactor->deltaT, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->sensingFactor->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->openLoopFactor->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->sensingFactor->data->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->openLoopFactor->data->data, status,
      CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( params->openLoopFactor->data->length ==
      params->sensingFactor->data->length, status,
      CALIBRATIONH_ESZMM, CALIBRATIONH_MSGESZMM );

  R0 = input->responseFunction->data->data;
  C0 = input->sensingFunction->data->data;
  R = output->responseFunction->data->data;
  C = output->sensingFunction->data->data;

  save = output->responseFunction->data;
  output->responseFunction = input->responseFunction;
  output->responseFunction->data = save;
  output->responseFunction->epoch = params->epoch;

  save = output->sensingFunction->data;
  output->sensingFunction = input->sensingFunction;
  output->sensingFunction->data = save;
  output->sensingFunction->epoch = params->epoch;

  /* locate correct values of a and ab */
  TRY( LALGPStoFloat( status->statusPtr, &epoch, &(params->epoch)), 
      status );
  TRY( LALGPStoFloat( status->statusPtr, &first_cal,
    &(params->sensingFactor->epoch)), status );
  TRY( LALGPStoFloat( status->statusPtr, &duration, &(params->duration)), 
      status );

  dt = epoch - first_cal;
  
  /* find the first point at or before the requested time */
  if ( (i_r4 = floor( dt / params->sensingFactor->deltaT ) ) < 0 )
  {
    ABORT( status, CALIBRATIONH_ETIME, CALIBRATIONH_MSGETIME );
  }
  else
  {
    i = (UINT4) i_r4;
  }
 
  /* compute the sum of the calibration factors */
  a.re = a.im = ab.re = ab.im = 0;
  length = 0;
  do
  {
    COMPLEX8 this_a;
    COMPLEX8 this_ab;
    
    if ( i > params->sensingFactor->data->length - 1 )
    {
      ABORT( status, CALIBRATIONH_ETIME, CALIBRATIONH_MSGETIME );
    }

    this_a = params->sensingFactor->data->data[i];
    this_ab = params->openLoopFactor->data->data[i];

    /* JC: I CHANGED THE LOGIC HERE TO WHAT I THOUGHT IT SHOULD BE! */
    if ( ( fabs( this_a.re ) < tiny && fabs( this_a.im ) < tiny ) ||
         ( fabs( this_ab.re ) < tiny && fabs( this_ab.im ) < tiny ) )
    {
      /* this is a hack for the broken S2 calibration frame data */
      if ( (params->epoch.gpsSeconds >= CAL_S2START) && 
          (params->epoch.gpsSeconds < CAL_S2END ) )
      {
        /* if the zero is during S2 print a warning... */
        LALSnprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
            "Zero calibration factors found during S2 at GPS %10.9f",
            first_cal + (REAL8) i * params->sensingFactor->deltaT );
        LALWarning( status, warnMsg );
      }
      else
      {
        /* ...or abort if we are outside S2 */
        LALSnprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
            "Zero calibration factor found at GPS %10.9f",
            first_cal + (REAL8) i * params->sensingFactor->deltaT );
        LALWarning( status, warnMsg );
        ABORT( status, CALIBRATIONH_EZERO, CALIBRATIONH_MSGEZERO );
      }
    }
    else
    {
      /* increment the count of factors if we are adding a non-zero value */
      ++length;
    }

    /* add this value to the sum */
    a.re += this_a.re;
    a.im += this_a.im;
    ab.re += this_ab.re;
    ab.im += this_ab.im;

    /* increment the calibration factor index */
    ++i;
  }
  while ( (first_cal + (REAL8) i * params->sensingFactor->deltaT) < 
      (epoch + duration) );

  /* if all the calibration factors are zero the abort */
  if ( ! length || 
      (fabs( a.re ) < tiny && fabs( a.im ) < tiny) ||
      (fabs( ab.re ) < tiny && fabs( ab.im ) < tiny) )
  {
    LALSnprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
        "Got %d calibration samples\nalpha and/or beta are zero:\n"
        "a.re = %e\ta.im = %e\nab.re = %e\tab.im = %e",
        length, a.re, a.im, ab.re, ab.im );
    LALWarning( status, warnMsg );
    ABORT( status, CALIBRATIONH_EZERO, CALIBRATIONH_MSGEZERO );
  }
  
  /* compute the mean of the calibration factors from the sum */
  a.re /= length;
  a.im /= length;
  ab.re /= length;
  ab.im /= length;

  /* return the used values of alpha and alphabeta */
  params->alpha = a;
  params->alphabeta = ab;
  LALSnprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
      "Got %d calibration samples\n"
      "a.re = %e\ta.im = %e\nab.re = %e\tab.im = %e",
      length, a.re, a.im, ab.re, ab.im );
  LALInfo( status, warnMsg );

  for ( i = 0; i < n; ++i )
  {
    COMPLEX8 tmp;
    /* compute the reference open loop function */
    tmp = cmul( C0[i], R0[i] );
    tmp.re -= 1;
    /* update the open loop function */
    tmp = cmul( ab, tmp );
    /* update the sensing function */
    C[i] = cmul( a, C0[i] );
    /* compute the updated response function */
    tmp.re += 1;
    R[i] = cdiv( tmp, C[i] );
  }
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/* <lalVerbatim file="ComputeTransferCP"> */
void
LALResponseConvert(
    LALStatus               *status,
    COMPLEX8FrequencySeries *output,
    COMPLEX8FrequencySeries *input
    )
{ /* </lalVerbatim> */
  COMPLEX8 tmpb, tmpc;
  REAL4 tmpx, tmpy;
  LALUnit unitOne;
  LALUnit unitTwo;
  UINT4 i;
  INT4 inv;
  INT4 fac;
  INT4 bad;

  INITSTATUS( status, "LALResponseConvert", COMPUTETRANSFERC );
  ATTATCHSTATUSPTR( status );

  output->epoch = input->epoch;

  /*
   * Interpolate to requested frequencies.
   * Just do linear interpolation of real and imag components.
   */
  for ( i = 0; i < output->data->length; ++i )
  {
    REAL4 x;
    UINT4 j;
    x = i * output->deltaF / input->deltaF;
    j = floor( x );
    if ( j > input->data->length - 2 )
      j = input->data->length - 2;
    x -= j;
    output->data->data[i].re = input->data->data[j].re 
      + x * ( input->data->data[j+1].re - input->data->data[j].re );
    output->data->data[i].im = input->data->data[j].im 
      + x * ( input->data->data[j+1].im - input->data->data[j].im );
  }


  /*
   * Use output units to figure out:
   *   1. Whether we want strain/ct or ct/strain
   *   2. Overall (power of ten) factor to apply.
   */

  /* determine if units need to be inverted or not (or if they are bad) */
  LALUnitNormalize( status->statusPtr, &unitOne, &output->sampleUnits );
  CHECKSTATUSPTR( status );
  LALUnitNormalize( status->statusPtr, &unitTwo, &input->sampleUnits );
  CHECKSTATUSPTR( status );
  bad = 0;
  inv = -1;
  for ( i = 0; i < LALNumUnits; ++i )
  {
    if ( unitOne.unitDenominatorMinusOne[i] != unitTwo.unitDenominatorMinusOne[i] )
    {
      bad = 1;
      break;
    }
    if ( unitOne.unitNumerator[i] == unitTwo.unitNumerator[i] )
    {
      if ( unitOne.unitNumerator[i] ) /* if this unit exists */
      {
        inv = 0; /* we don't need to invert */
        if ( inv == 1 ) /* error: some units need to be inverted, others not */
        {
          bad = 1;
          break;
        }
      }
    }
    else
    {
      if ( unitOne.unitNumerator[i] == -unitTwo.unitNumerator[i] )
      {
        /* this unit needs to be inverted */
        inv = 1;
      }
      else /* error: output units not equal to input or inverse of input */
      {
        bad = 1;
        break;
      }
    }
  }
  if ( bad ) /* units were bad: abort */
  {
    ABORT( status, CALIBRATIONH_EUNIT, CALIBRATIONH_MSGEUNIT );
  }

  /* determine if there is a scale factor that needs to be applied */
  fac = unitOne.powerOfTen - ( inv ? -unitTwo.powerOfTen : unitTwo.powerOfTen );

  /* perform conversion(s) */

  if ( inv ) /* invert data */
  {
    for ( i = 0; i < output->data->length; ++i )
    {
      output->data->data[i] = cinv( output->data->data[i] );
    }
  }

  if ( fac ) /* scale data */
  {
    REAL4 scale = pow( 10.0, -fac );
    for ( i = 0; i < output->data->length; ++i )
    {
      output->data->data[i].re *= scale;
      output->data->data[i].im *= scale;
    }
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
