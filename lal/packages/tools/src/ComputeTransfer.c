/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Patrick Brady, Saikat Ray-Majumder, Stephen Fairhurst
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
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Calibration.h>
#include <lal/Units.h>
#include <lal/Date.h>

#define CAL_S2START 729273613
#define CAL_S2END 734367613

/**
 * \defgroup ComputeTransfer_c Module ComputeTransfer.c
 * \ingroup Calibration_h
 * \author Patrick Brady, Jolien Creighton
 *
 * \brief Computes the transfer function from zero-pole-gain representation.
 *
 * A transfer function can either be specified as a list of coefficients or a
 * list of poles and zeros. The function LALComputeTransfer() computes the
 * frequency representation of the transfer function <tt>calrec->transfer</tt>
 * described by the zeroes,  poles and gain in <tt>*calrec</tt>.   The memory for
 * the frequency series should be allocated before calling this routine which uses
 * <tt>calrec->transfer->deltaF</tt> and <tt>calrec->transfer->data->npoints</tt>.
 *
 * The routine LALUpdateCalibration() updates the response function
 * and the sensing function from some reference functions to the current
 * functions using information about the calibration lines.  The two calibration
 * lines yield two constants (as a slowly-varying function of time) that are
 * used as coefficients to the reference response and sensing functions to
 * compute the current response and sensing functions.  These coefficients are
 * stored in time series in the parameter structure, along with the current
 * epoch and duration for which the calibration functions are to be computed.  If
 * the duration is zero, the calibration factors are the first ones at or after
 * the given epoch.  If the duration is non-zero, then the calibration factors are
 * an average of all calibrations between epoch and epoch + duration.
 *
 * The routine LALResponseConvert() takes a given frequency series
 * and converts it to a new frequency series by performing the following steps:
 * (i) the original frequency series is interpolated (using linear interpolation
 * of the real and imaginary parts independently) to the frequencies required
 * in the output frequency series;  (ii) if the output frequency series has units
 * that are the inverse of those of the input frequency series, the data is
 * inverted;  (iii) the data is scaled by an appropriate power of ten computed
 * from the input units and the output units.  For example you can convert from
 * strain per count to counts per atto-strain.
 *
 * \heading{Algorithm}
 *
 * The transfer function is deduced from the poles and zeros as follows:
 * \f{equation}{
 * T(f) = c_{\mathrm{gain}}
 * {\prod_i \textrm{zero}(f,z_i)}{ \prod_{i} \textrm{pole}(f,p_i)}
 * \f}
 * where
 * \f{equation}{
 * \textrm{zero}(f,z) = \left\{ \begin{array}{ll}
 * (i f / z) + 1 & \textrm{ when } z\neq 0 \\
 * i f & \textrm{ when } z = 0
 * \end{array}
 * \right.
 * \f}
 * and
 * \f{equation}{
 * \textrm{pole}(f,p) = \left\{ \begin{array}{ll}
 * \frac{1}{(i f / p) + 1} & \textrm{ when } p \neq 0 \\
 * \frac{1}{i f} & \textrm{ when } p = 0
 * \end{array}
 * \right.
 * \f}
 * For both the transfer function and the pole-zero notation the units for
 * frequency is Hz rather than rad/s (angular frequency).  In particular, poles
 * and zeros are specified by their location in frequency space
 *
 * To update the response function from one epoch to another, two functions are
 * needed.  These are the sensing function \f$C(f)\f$ and the response function
 * \f$R(f)\f$, which are related by
 * \f{equation}{
 * R(f) = \frac{1+H(f)}{C(f)}
 * \f}
 * where \f$H(f)\f$ is the open-loop gain function.  If the sensing function and
 * the open-loop gain function are known at some reference time (\f$C_0(f)\f$ and
 * \f$H_0(f)\f$) then the sensing function and open-loop gain function can be
 * calculated at a later time.  They are \f$C(f)=\alpha C_0(f)\f$ and
 * \f$H(f)=\alpha\beta H_0(f)\f$ where \f$\alpha\f$ and \f$\beta\f$ are slowly varying
 * coefficients that account for overall changes in the gains of the sensing
 * function and the open-loop gain.  The coefficients \f$\alpha\f$ and \f$\alpha\beta\f$
 * can be determined, as slowly-varying functions of time, by monitoring the
 * two injected calibration lines.  Thus, an updated sensing function and response
 * function can be computed from reference sensing function and response function,
 * \f$C_0(f)\f$ and \f$R_0(f)\f$ via:
 * \f{equation}{
 * C(f) = \alpha C_0(f)
 * \f}
 * and
 * \f{equation}{
 * R(f) = \frac{1+\alpha\beta[C_0(f)R_0(f)-1]}{\alpha C_0(f)}
 * \f}
 * where \f$\alpha\f$ and \f$\beta\f$ are those values of the coefficients that are
 * appropriate for the particular epoch.
 *
 * \heading{Uses}
 *
 * \heading{Notes}
 * The DC component of <tt>calrec->transfer</tt> is always filled with \f$1 + i 0\f$.
 * In most cases,  this should be irrelevant for gravitational wave data analysis,
 * but care should be taken if DC is relevant when this function is used.
 *
 */
/*@{*/

/** UNDOCUMENTED */
void
LALComputeTransfer( LALStatus                 *stat,
                    CalibrationRecord         *calrec
                    )

{
  UINT4         i, j;                    /* indexes               */
  UINT4         jmin;                    /* used to handle DC     */
  REAL4         f,df;                    /* freq and interval     */
  REAL8         norm;

  INITSTATUS(stat);
  ATTATCHSTATUSPTR (stat);

  /* Make sure parameter structures and their fields exist. */
  ASSERT( calrec, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( calrec->zeros, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( calrec->poles, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );
  ASSERT( calrec->transfer, stat, CALIBRATIONH_ENULL, CALIBRATIONH_MSGENULL );

  /* initialise everything */
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
    calrec->transfer->data->data[0] = 1.0;
    jmin = 1;
  }

  /* loop over frequency in the output */
  for ( j=jmin ; j<calrec->transfer->data->length ; j++)
  {
    COMPLEX8 T = 1.0;
    /* the frequency */
    f = calrec->transfer->f0 + (REAL4) j * df;

    /* loop over zeroes */
    for (i = 0 ; i < calrec->zeros->length ; i++)
      T *= calrec->zeros->data[i] + I * f;

    /* loop over poles */
    for (i = 0 ; i < calrec->poles->length ; i++)
      T /= calrec->zeros->data[i] + I * f;

    /* fill the frequency series */
    calrec->transfer->data->data[j] = norm * T;
  }


  /* we're out of here */
  DETATCHSTATUSPTR (stat);
  RETURN( stat );
}


/** UNDOCUMENTED */
void
LALUpdateCalibration(
    LALStatus               *status,
    CalibrationFunctions    *output,
    CalibrationFunctions    *input,
    CalibrationUpdateParams *params
    )
{
  const REAL4 tiny = 1e-6;
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

  INITSTATUS(status);
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
  epoch = XLALGPSGetREAL8(&(params->epoch));
  first_cal = XLALGPSGetREAL8(&(params->sensingFactor->epoch));
  duration = XLALGPSGetREAL8(&(params->duration));

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
  a = ab = 0;
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
    if ( ( fabs( creal(this_a) ) < tiny && fabs( cimag(this_a) ) < tiny ) ||
         ( fabs( creal(this_ab) ) < tiny && fabs( cimag(this_ab) ) < tiny ) )
    {
      /* this is a hack for the broken S2 calibration frame data */
      if ( (params->epoch.gpsSeconds >= CAL_S2START) &&
          (params->epoch.gpsSeconds < CAL_S2END ) )
      {
        /* if the zero is during S2 print a warning... */
        snprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
            "Zero calibration factors found during S2 at GPS %10.9f",
            first_cal + (REAL8) i * params->sensingFactor->deltaT );
        LALWarning( status, warnMsg );
      }
      else
      {
        /* ...or abort if we are outside S2 */
        snprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
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
    a += this_a;
    ab += this_ab;

    /* increment the calibration factor index */
    ++i;
  }
  while ( (first_cal + (REAL8) i * params->sensingFactor->deltaT) <
      (epoch + duration) );

  /* if all the calibration factors are zero the abort */
  if ( ! length ||
      (fabs( creal(a) ) < tiny && fabs( cimag(a) ) < tiny) ||
      (fabs( creal(ab) ) < tiny && fabs( cimag(ab) ) < tiny) )
  {
    snprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
        "Got %d calibration samples\nalpha and/or beta are zero:\n"
        "Re a = %e\tIm a = %e\nRe ab = %e\tIm ab = %e",
        length, creal(a), cimag(a), creal(ab), cimag(ab) );
    LALWarning( status, warnMsg );
    ABORT( status, CALIBRATIONH_EZERO, CALIBRATIONH_MSGEZERO );
  }

  /* compute the mean of the calibration factors from the sum */
  a /= length;
  ab /= length;

  /* return the used values of alpha and alphabeta */
  params->alpha = a;
  params->alphabeta = ab;
  snprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
      "Got %d calibration samples\n"
      "Re a = %e\tIm a = %e\nRe ab = %e\tIm ab = %e",
      length, creal(a), cimag(a), creal(ab), cimag(ab) );
  LALInfo( status, warnMsg );

  for ( i = 0; i < n; ++i )
  {
    C[i] = a * C0[i];
    R[i] = (ab * (C0[i] * R0[i] - 1.0) + 1.0) / C[i];
  }
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/** UNDOCUMENTED */
void
LALResponseConvert(
    LALStatus               *status,
    COMPLEX8FrequencySeries *output,
    COMPLEX8FrequencySeries *input
    )
{
  LALUnit unitOne;
  LALUnit unitTwo;
  UINT4 i;
  INT4 inv;
  INT4 fac;
  INT4 bad;

  INITSTATUS(status);
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
    output->data->data[i] = input->data->data[j]
      + x * ( input->data->data[j+1] - input->data->data[j] );
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
      output->data->data[i] = 1.0 / output->data->data[i];
    }
  }

  if ( fac ) /* scale data */
  {
    REAL4 scale = pow( 10.0, -fac );
    for ( i = 0; i < output->data->length; ++i )
    {
      output->data->data[i] *= scale;
    }
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/** UNDOCUMENTED */
INT4
XLALResponseConvert(
    COMPLEX8FrequencySeries *output,
    COMPLEX8FrequencySeries *input
    )
{
  LALUnit unitOne;
  LALUnit unitTwo;
  UINT4 i;
  INT4 inv;
  INT4 fac;
  INT4 bad;

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
    output->data->data[i] = input->data->data[j]
      + x * ( input->data->data[j+1] - input->data->data[j] );
  }


  /*
   * Use output units to figure out:
   *   1. Whether we want strain/ct or ct/strain
   *   2. Overall (power of ten) factor to apply.
   */

  /* determine if units need to be inverted or not (or if they are bad) */
  XLALUnitNormalize( &output->sampleUnits );
  XLALUnitNormalize( &input->sampleUnits );
  unitOne = output->sampleUnits;
  unitTwo = input->sampleUnits;

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
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* determine if there is a scale factor that needs to be applied */
  fac = unitOne.powerOfTen - ( inv ? -unitTwo.powerOfTen : unitTwo.powerOfTen );

  /* perform conversion(s) */

  if ( inv ) /* invert data */
  {
    for ( i = 0; i < output->data->length; ++i )
    {
      output->data->data[i] = 1.0 / output->data->data[i];
    }
  }

  if ( fac ) /* scale data */
  {
    REAL4 scale = pow( 10.0, -fac );
    for ( i = 0; i < output->data->length; ++i )
    {
      output->data->data[i] *= scale;
    }
  }

  return 0;
}

/*@}*/
