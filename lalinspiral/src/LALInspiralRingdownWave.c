/*
*  Copyright (C) 2008 Yi Pan, B.S. Sathyaprakash (minor modificaitons)
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
 * \author Yi Pan
 * \file
 *
 * \brief Module to compute the ring-down waveform as linear combination
 * of quasi-normal-modes decaying waveforms, which can be attached to
 * the inspiral part of the compat binary coalescing waveform.
 *
 * \heading{Prototypes}
 *
 * <tt>XLALXLALInspiralRingdownWave()</tt>
 * <ul>
 * <li> <tt>rdwave1,</tt> Output, the real part of the ring-down waveform
 * </li><li> <tt>rdwave2,</tt> Output, the imaginary part of the ring-down waveform
 * </li><li> <tt>params,</tt> Input, the parameters where ring-down waveforms are computed
 * </li><li> <tt>inspwave1,</tt> Input, the real part of the ring-down waveform
 * </li><li> <tt>inspwave2,</tt> Input, the real part of the ring-down waveform
 * </li><li> <tt>modefreqs,</tt> Input, the frequencies of the quasi-normal-modes
 * </li><li> <tt>nmode,</tt> Input, the number of quasi-normal-modes to be combined.</li>
 * </ul>
 *
 * <tt>XLALGenerateWaveDerivatives()</tt>
 * <ul>
 * <li> <tt>dwave,</tt> Output, time derivative of the input waveform
 * </li><li> <tt>ddwave,</tt> Output, two time derivative of the input waveform
 * </li><li> <tt>wave,</tt> Input, waveform to be differentiated in time
 * </li><li> <tt>params,</tt> Input, the parameters of the input waveform.</li>
 * </ul>
 *
 * <tt>XLALGenerateQNMFreq()</tt>
 * <ul>
 * <li> <tt>ptfwave,</tt> Output, the frequencies of the quasi-normal-modes
 * </li><li> <tt>params,</tt> Input, the parameters of the binary system
 * </li><li> <tt>l,</tt> Input, the l of the modes
 * </li><li> <tt>m,</tt> Input, the m of the modes
 * </li><li> <tt>nmodes,</tt> Input, the number of overtones considered.</li>
 * </ul>
 *
 * <tt>XLALFinalMassSpin()</tt>
 * <ul>
 * <li> <tt>finalMass,</tt> Output, the mass of the final Kerr black hole
 * </li><li> <tt>finalSpin,</tt>  Input, the spin of the final Kerr balck hole
 * </li><li> <tt>params,</tt> Input, the parameters of the binary system.</li>
 * </ul>
 *
 * \heading{Description}
 * Generating ring-down waveforms.
 *
 * \heading{Algorithm}
 *
 * \heading{Uses}
 *
 * \code
 * LALMalloc
 * LALFree
 * \endcode
 *
 * \heading{Notes}
 *
 */

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>
#include <lal/LALConstants.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


static INT4 XLALFinalMassSpin(
	REAL8		 *finalMass,
	REAL8		 *finalSpin,
	InspiralTemplate *params
	);


INT4 XLALInspiralHybridRingdownWave (
	REAL4Vector			*rdwave1,
	REAL4Vector			*rdwave2,
	InspiralTemplate		*params,
	REAL4VectorSequence		*inspwave1,
	REAL4VectorSequence		*inspwave2,
	COMPLEX8Vector			*modefreqs,
	REAL8Vector			*matchrange
	)
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;

  /* For checking GSL return codes */
  INT4 gslStatus;

  UINT4 i, j, k, nmodes = 8;

  /* Sampling rate from input */
  REAL8 dt, t1, t2, t3, t4, t5, rt;
  gsl_matrix *coef;
  gsl_vector *hderivs;
  gsl_vector *x;
  gsl_permutation *p;
  REAL8Vector *modeamps;
  int s;
  REAL8 tj;
  REAL8 m;

  dt = 1.0 / params -> tSampling;
  m  = (params->mass1 + params->mass2) * LAL_MTSUN_SI;
  t5 = (matchrange->data[0] - matchrange->data[1]) * m;
  rt = -t5 / 5.;

  t4 = t5 + rt;
  t3 = t4 + rt;
  t2 = t3 + rt;
  t1 = t2 + rt;
  
  if ( inspwave1->length != 3 || inspwave2->length != 3 ||
		modefreqs->length != nmodes )
  {
    XLAL_ERROR( XLAL_EBADLEN );
  }

  /* Solving the linear system for QNMs amplitude coefficients using gsl routine */
  /* Initiate matrices and supporting variables */
  XLAL_CALLGSL( coef = (gsl_matrix *) gsl_matrix_alloc(2 * nmodes, 2 * nmodes) );
  XLAL_CALLGSL( hderivs = (gsl_vector *) gsl_vector_alloc(2 * nmodes) );
  XLAL_CALLGSL( x = (gsl_vector *) gsl_vector_alloc(2 * nmodes) );
  XLAL_CALLGSL( p = (gsl_permutation *) gsl_permutation_alloc(2 * nmodes) );

  /* Check all matrices and variables were allocated */
  if ( !coef || !hderivs || !x || !p )
  {
    if (coef)    gsl_matrix_free(coef);
    if (hderivs) gsl_vector_free(hderivs);
    if (x)       gsl_vector_free(x);
    if (p)       gsl_permutation_free(p);

    XLAL_ERROR( XLAL_ENOMEM );
  }

  /* Define the linear system Ax=y */
  /* Matrix A (2*n by 2*n) has block symmetry. Define half of A here as "coef" */
  /* Define y here as "hderivs" */
  for (i = 0; i < nmodes; ++i)
  {
	gsl_matrix_set(coef, 0, i, 1);
	gsl_matrix_set(coef, 1, i, - cimagf(modefreqs->data[i]));
	gsl_matrix_set(coef, 2, i, exp(-cimagf(modefreqs->data[i])*t1) * cos(crealf(modefreqs->data[i])*t1));
	gsl_matrix_set(coef, 3, i, exp(-cimagf(modefreqs->data[i])*t2) * cos(crealf(modefreqs->data[i])*t2));
	gsl_matrix_set(coef, 4, i, exp(-cimagf(modefreqs->data[i])*t3) * cos(crealf(modefreqs->data[i])*t3));
	gsl_matrix_set(coef, 5, i, exp(-cimagf(modefreqs->data[i])*t4) * cos(crealf(modefreqs->data[i])*t4));
	gsl_matrix_set(coef, 6, i, exp(-cimagf(modefreqs->data[i])*t5) * cos(crealf(modefreqs->data[i])*t5));
	gsl_matrix_set(coef, 7, i, exp(-cimagf(modefreqs->data[i])*t5) * 
				      (-cimagf(modefreqs->data[i]) * cos(crealf(modefreqs->data[i])*t5)
				       -crealf(modefreqs->data[i]) * sin(crealf(modefreqs->data[i])*t5)));
	gsl_matrix_set(coef, 8, i, 0);
	gsl_matrix_set(coef, 9, i, - crealf(modefreqs->data[i]));
	gsl_matrix_set(coef, 10, i, -exp(-cimagf(modefreqs->data[i])*t1) * sin(crealf(modefreqs->data[i])*t1));
	gsl_matrix_set(coef, 11, i, -exp(-cimagf(modefreqs->data[i])*t2) * sin(crealf(modefreqs->data[i])*t2));
	gsl_matrix_set(coef, 12, i, -exp(-cimagf(modefreqs->data[i])*t3) * sin(crealf(modefreqs->data[i])*t3));
	gsl_matrix_set(coef, 13, i, -exp(-cimagf(modefreqs->data[i])*t4) * sin(crealf(modefreqs->data[i])*t4));
	gsl_matrix_set(coef, 14, i, -exp(-cimagf(modefreqs->data[i])*t5) * sin(crealf(modefreqs->data[i])*t5));
	gsl_matrix_set(coef, 15, i, exp(-cimagf(modefreqs->data[i])*t5) * 
				      ( cimagf(modefreqs->data[i]) * sin(crealf(modefreqs->data[i])*t5)
				       -crealf(modefreqs->data[i]) * cos(crealf(modefreqs->data[i])*t5)));
  }
  for (i = 0; i < 2; ++i)
  {
	gsl_vector_set(hderivs, i, inspwave1->data[(i + 1) * inspwave1->vectorLength - 1]);
	gsl_vector_set(hderivs, i + nmodes, inspwave2->data[(i + 1) * inspwave2->vectorLength - 1]);
	gsl_vector_set(hderivs, i + 6, inspwave1->data[i * inspwave1->vectorLength]);
	gsl_vector_set(hderivs, i + 6 + nmodes, inspwave2->data[i * inspwave2->vectorLength]);
  }
  gsl_vector_set(hderivs, 2, inspwave1->data[4]);
  gsl_vector_set(hderivs, 2 + nmodes, inspwave2->data[4]);
  gsl_vector_set(hderivs, 3, inspwave1->data[3]);
  gsl_vector_set(hderivs, 3 + nmodes, inspwave2->data[3]);
  gsl_vector_set(hderivs, 4, inspwave1->data[2]);
  gsl_vector_set(hderivs, 4 + nmodes, inspwave2->data[2]);
  gsl_vector_set(hderivs, 5, inspwave1->data[1]);
  gsl_vector_set(hderivs, 5 + nmodes, inspwave2->data[1]);
  
  /* Complete the definition for the rest half of A */
  for (i = 0; i < nmodes; ++i)
  {
	for (k = 0; k < nmodes; ++k)
	{
	  gsl_matrix_set(coef, i, k + nmodes, - gsl_matrix_get(coef, i + nmodes, k));
	  gsl_matrix_set(coef, i + nmodes, k + nmodes, gsl_matrix_get(coef, i, k));
	}
  }

#if 0
  /* print ringdown-matching linear system: coefficient matrix and RHS vector */
  printf("\nRingdown matching matrix:\n");
  for (i = 0; i < 16; ++i)
  {
    for (j = 0; j < 16; ++j)
    {
      printf("%.12e ",gsl_matrix_get(coef,i,j));
    }
    printf("\n");
  }
  printf("RHS:  ");
  for (i = 0; i < 16; ++i)
  {
    printf("%.12e   ",gsl_vector_get(hderivs,i));
  }
  printf("\n");
 #endif

  /* Call gsl LU decomposition to solve the linear system */
  XLAL_CALLGSL( gslStatus = gsl_linalg_LU_decomp(coef, p, &s) );
  if ( gslStatus == GSL_SUCCESS )
  {
    XLAL_CALLGSL( gslStatus = gsl_linalg_LU_solve(coef, p, hderivs, x) );
  }
  if ( gslStatus != GSL_SUCCESS )
  {
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Putting solution to an XLAL vector */
  modeamps = XLALCreateREAL8Vector(2 * nmodes);

  if ( !modeamps )
  {
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    XLAL_ERROR( XLAL_ENOMEM );
  }

  for (i = 0; i < nmodes; ++i)
  {
	modeamps->data[i] = gsl_vector_get(x, i);
	modeamps->data[i + nmodes] = gsl_vector_get(x, i + nmodes);
  }

  /* Free all gsl linear algebra objects */
  gsl_matrix_free(coef);
  gsl_vector_free(hderivs);
  gsl_vector_free(x);
  gsl_permutation_free(p);

  /* Build ring-down waveforms */

  REAL8 timeOffset = fmod( matchrange->data[1], dt/m) * dt;

  for (j = 0; j < rdwave1->length; ++j)
  {
	tj = j * dt - timeOffset;
	rdwave1->data[j] = 0;
	rdwave2->data[j] = 0;
	for (i = 0; i < nmodes; ++i)
	{
	  rdwave1->data[j] += exp(- tj * cimagf(modefreqs->data[i]))
			* ( modeamps->data[i] * cos(tj * crealf(modefreqs->data[i]))
			+   modeamps->data[i + nmodes] * sin(tj * crealf(modefreqs->data[i])) );
	  rdwave2->data[j] += exp(- tj * cimagf(modefreqs->data[i]))
			* (- modeamps->data[i] * sin(tj * crealf(modefreqs->data[i]))
			+   modeamps->data[i + nmodes] * cos(tj * crealf(modefreqs->data[i])) );
	}
  }

  XLALDestroyREAL8Vector(modeamps);
  return errcode;
}

INT4 XLALInspiralRingdownWave (
	REAL4Vector			*rdwave1,
	REAL4Vector			*rdwave2,
	InspiralTemplate		*params,
	REAL4VectorSequence		*inspwave1,
	REAL4VectorSequence		*inspwave2,
	COMPLEX8Vector			*modefreqs,
	UINT4				nmodes
	)

{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;

  /* For checking GSL return codes */
  INT4 gslStatus;

  UINT4 i, j, k;

  /* Sampling rate from input */
  REAL8 dt;
  gsl_matrix *coef;
  gsl_vector *hderivs;
  gsl_vector *x;
  gsl_permutation *p;
  REAL8Vector *modeamps;
  int s;
  REAL8 tj;

  dt = 1.0 / params -> tSampling;

  if ( inspwave1->length != nmodes || inspwave2->length != nmodes ||
		modefreqs->length != nmodes )
  {
    XLAL_ERROR( XLAL_EBADLEN );
  }

  /* Solving the linear system for QNMs amplitude coefficients using gsl routine */
  /* Initiate matrices and supporting variables */
  XLAL_CALLGSL( coef = (gsl_matrix *) gsl_matrix_alloc(2 * nmodes, 2 * nmodes) );
  XLAL_CALLGSL( hderivs = (gsl_vector *) gsl_vector_alloc(2 * nmodes) );
  XLAL_CALLGSL( x = (gsl_vector *) gsl_vector_alloc(2 * nmodes) );
  XLAL_CALLGSL( p = (gsl_permutation *) gsl_permutation_alloc(2 * nmodes) );

  /* Check all matrices and variables were allocated */
  if ( !coef || !hderivs || !x || !p )
  {
    if (coef)    gsl_matrix_free(coef);
    if (hderivs) gsl_vector_free(hderivs);
    if (x)       gsl_vector_free(x);
    if (p)       gsl_permutation_free(p);

    XLAL_ERROR( XLAL_ENOMEM );
  }

  /* Define the linear system Ax=y */
  /* Matrix A (2*n by 2*n) has block symmetry. Define half of A here as "coef" */
  /* Define y here as "hderivs" */
  for (i = 0; i < nmodes; ++i)
  {
	gsl_matrix_set(coef, 0, i, 1);
	gsl_matrix_set(coef, 1, i, - cimagf(modefreqs->data[i]));
	gsl_matrix_set(coef, 2, i, cimagf(modefreqs->data[i]) * cimagf(modefreqs->data[i])
			- crealf(modefreqs->data[i]) * crealf(modefreqs->data[i]));
	gsl_matrix_set(coef, 3, i, 0);
	gsl_matrix_set(coef, 4, i, - crealf(modefreqs->data[i]));
	gsl_matrix_set(coef, 5, i,  2 * crealf(modefreqs->data[i]) * cimagf(modefreqs->data[i]));

	gsl_vector_set(hderivs, i, inspwave1->data[(i + 1) * inspwave1->vectorLength - 1]);
	gsl_vector_set(hderivs, i + nmodes, inspwave2->data[(i + 1) * inspwave2->vectorLength - 1]);
  }
  /* Complete the definition for the rest half of A */
  for (i = 0; i < nmodes; ++i)
  {
	for (k = 0; k < nmodes; ++k)
	{
	  gsl_matrix_set(coef, i, k + nmodes, - gsl_matrix_get(coef, i + nmodes, k));
	  gsl_matrix_set(coef, i + nmodes, k + nmodes, gsl_matrix_get(coef, i, k));
	}
  }

  /* Call gsl LU decomposition to solve the linear system */
  XLAL_CALLGSL( gslStatus = gsl_linalg_LU_decomp(coef, p, &s) );
  if ( gslStatus == GSL_SUCCESS )
  {
    XLAL_CALLGSL( gslStatus = gsl_linalg_LU_solve(coef, p, hderivs, x) );
  }

  if ( gslStatus != GSL_SUCCESS )
  {
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Putting solution to an XLAL vector */
  modeamps = XLALCreateREAL8Vector(2 * nmodes);

  if ( !modeamps )
  {
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    XLAL_ERROR( XLAL_ENOMEM );
  }

  for (i = 0; i < nmodes; ++i)
  {
	modeamps->data[i] = gsl_vector_get(x, i);
	modeamps->data[i + nmodes] = gsl_vector_get(x, i + nmodes);
  }

  /* Free all gsl linear algebra objects */
  gsl_matrix_free(coef);
  gsl_vector_free(hderivs);
  gsl_vector_free(x);
  gsl_permutation_free(p);

  /* Build ring-down waveforms */
  for (j = 0; j < rdwave1->length; ++j)
  {
	tj = j * dt;
	rdwave1->data[j] = 0;
	rdwave2->data[j] = 0;
	for (i = 0; i < nmodes; ++i)
	{
	  rdwave1->data[j] += exp(- tj * cimagf(modefreqs->data[i]))
			* ( modeamps->data[i] * cos(tj * crealf(modefreqs->data[i]))
			+   modeamps->data[i + nmodes] * sin(tj * crealf(modefreqs->data[i])) );
	  rdwave2->data[j] += exp(- tj * cimagf(modefreqs->data[i]))
			* (- modeamps->data[i] * sin(tj * crealf(modefreqs->data[i]))
			+   modeamps->data[i + nmodes] * cos(tj * crealf(modefreqs->data[i])) );
	}
  }

  XLALDestroyREAL8Vector(modeamps);
  return errcode;
}

INT4 XLALGenerateHybridWaveDerivatives (
	REAL4Vector				*rwave,
	REAL4Vector				*dwave,
	REAL4Vector				*ddwave,
        REAL8Vector				*timeVec,
	REAL4Vector				*wave,
	REAL8Vector				*matchrange,
	InspiralTemplate			*params
	)
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;

  /* For checking GSL return codes */
  INT4 gslStatus;

  UINT4 j;
  UINT4 vecLength;
  REAL8 dt;
  REAL8 m;
  double *y;
  double ry, dy, dy2;
  double rt;
  double *tlist;
  gsl_interp_accel *acc;
  gsl_spline *spline;

  /* Sampling rate from input */
  dt = 1.0 / params -> tSampling;
  m  = (params->mass1 + params->mass2) * LAL_MTSUN_SI;

  tlist = (double *) LALMalloc(6 * sizeof(double));
  rt = (matchrange->data[1] - matchrange->data[0]) / 5.;
  tlist[0] = matchrange->data[0];
  tlist[1] = tlist[0] + rt;
  tlist[2] = tlist[1] + rt;
  tlist[3] = tlist[2] + rt;
  tlist[4] = tlist[3] + rt;
  tlist[5] = matchrange->data[1];

  /* Set the length of the interpolation vectors */
  vecLength = ( m * matchrange->data[2] / dt ) + 1;

  /* Getting interpolation and derivatives of the waveform using gsl spline routine */
  /* Initiate arrays and supporting variables for gsl */
  y = (double *) LALMalloc(vecLength * sizeof(double));

  if ( !y )
  {
    XLAL_ERROR( XLAL_ENOMEM );
  }

  for (j = 0; j < vecLength; ++j)
  {
	y[j] = wave->data[j];
  }


  XLAL_CALLGSL( acc = (gsl_interp_accel*) gsl_interp_accel_alloc() );
  XLAL_CALLGSL( spline = (gsl_spline*) gsl_spline_alloc(gsl_interp_cspline, vecLength) );
  if ( !acc || !spline )
  {
    if ( acc )    gsl_interp_accel_free(acc);
    if ( spline ) gsl_spline_free(spline);
    LALFree( y );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  /* Gall gsl spline interpolation */
  gslStatus = gsl_spline_init(spline, timeVec->data, y, vecLength);
  if ( gslStatus != GSL_SUCCESS )
  { 
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    LALFree( y );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Getting first and second order time derivatives from gsl interpolations */
  for (j = 0; j < 6; ++j)
  {
    gslStatus = gsl_spline_eval_e( spline, tlist[j], acc, &ry );
    if ( gslStatus == GSL_SUCCESS )
    {
      gslStatus = gsl_spline_eval_deriv_e(spline, tlist[j], acc, &dy );
      gslStatus = gsl_spline_eval_deriv2_e(spline, tlist[j], acc, &dy2 );
    }
    if (gslStatus != GSL_SUCCESS )
    {
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
      LALFree( y );
      XLAL_ERROR( XLAL_EFUNC );
    }
    rwave->data[j]  = (REAL4)(ry);
    dwave->data[j]  = (REAL4)(dy/m);
    ddwave->data[j] = (REAL4)(dy2/m/m);

  }
  
  /* Free gsl variables */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  LALFree( tlist );
  LALFree(y);

  return errcode;
}

INT4 XLALGenerateWaveDerivatives (
	REAL4Vector			*dwave,
	REAL4Vector			*ddwave,
	REAL4Vector			*wave,
	InspiralTemplate		*params
	)

{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;

  /* For checking GSL return codes */
  INT4 gslStatus;

  UINT4 j;
  REAL8 dt;
  double *x, *y;
  double dy, dy2;
  gsl_interp_accel *acc;
  gsl_spline *spline;

  /* Sampling rate from input */
  dt = 1.0 / params -> tSampling;

  /* Getting interpolation and derivatives of the waveform using gsl spline routine */
  /* Initiate arrays and supporting variables for gsl */
  x = (double *) LALMalloc(wave->length * sizeof(double));
  y = (double *) LALMalloc(wave->length * sizeof(double));

  if ( !x || !y )
  {
    if ( x ) LALFree (x);
    if ( y ) LALFree (y);
    XLAL_ERROR( XLAL_ENOMEM );
  }

  for (j = 0; j < wave->length; ++j)
  {
	x[j] = j;
	y[j] = wave->data[j];
  }

  XLAL_CALLGSL( acc = (gsl_interp_accel*) gsl_interp_accel_alloc() );
  XLAL_CALLGSL( spline = (gsl_spline*) gsl_spline_alloc(gsl_interp_cspline, wave->length) );
  if ( !acc || !spline )
  {
    if ( acc )    gsl_interp_accel_free(acc);
    if ( spline ) gsl_spline_free(spline);
    LALFree( x );
    LALFree( y );
    XLAL_ERROR( XLAL_ENOMEM );
  }

  /* Gall gsl spline interpolation */
  XLAL_CALLGSL( gslStatus = gsl_spline_init(spline, x, y, wave->length) );
  if ( gslStatus != GSL_SUCCESS )
  { 
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    LALFree( x );
    LALFree( y );
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Getting first and second order time derivatives from gsl interpolations */
  for (j = 0; j < wave->length; ++j)
  {
    XLAL_CALLGSL(gslStatus = gsl_spline_eval_deriv_e( spline, j, acc, &dy ) );
    if ( gslStatus == GSL_SUCCESS )
    {
      XLAL_CALLGSL(gslStatus = gsl_spline_eval_deriv2_e(spline, j, acc, &dy2 ) );
    }
    if (gslStatus != GSL_SUCCESS )
    {
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
      LALFree( x );		
      LALFree( y );
      XLAL_ERROR( XLAL_EFUNC );
    }
    dwave->data[j]  = (REAL4)(dy / dt);
    ddwave->data[j] = (REAL4)(dy2 / dt / dt);

  }

  /* Free gsl variables */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  LALFree(x);
  LALFree(y);

  return errcode;
}

INT4 XLALGenerateQNMFreq(
	COMPLEX8Vector		*modefreqs,
	InspiralTemplate	*params,
	UINT4			UNUSED l,
	UINT4			UNUSED m,
	UINT4			nmodes
	)

{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  UINT4 i;
  REAL8 totalMass, finalMass, finalSpin;
  /* Fitting coefficients for QNM frequencies from PRD73, 064030 */
  REAL4 BCWre[8][3] = { {1.5251, -1.1568,  0.1292}, {1.3673, -1.0260,  0.1628}, { 1.3223, -1.0257,  0.1860} };
  REAL4 BCWim[8][3] = { {0.7000,  1.4187, -0.4990}, {0.1000,  0.5436, -0.4731}, {-0.1000,  0.4206, -0.4256} };

  /* Get a local copy of the intrinstic parameters */
  totalMass = params->totalMass;
  finalMass = 0;
  finalSpin = 0;

  /* Call XLALFinalMassSpin() to get mass and spin of the final black hole */
  errcode = XLALFinalMassSpin(&finalMass, &finalSpin, params);
  if ( errcode != XLAL_SUCCESS )
  {
	  XLAL_ERROR( XLAL_EFUNC );
  }

  /* QNM frequencies from the fitting given in PRD73, 064030 */
  for (i = 0; i < nmodes; ++i)
  {
	modefreqs->data[i].realf_FIXME = BCWre[i][0] + BCWre[i][1] * pow(1.- finalSpin, BCWre[i][2]);
	modefreqs->data[i].imagf_FIXME = crealf(modefreqs->data[i]) / 2
			     / (BCWim[i][0] + BCWim[i][1] * pow(1.- finalSpin, BCWim[i][2]));
	modefreqs->data[i].realf_FIXME *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
	modefreqs->data[i].imagf_FIXME *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
  }
  return errcode;
}


/**
 * As with the above function, this generates the quasinormal mode frequencies for a black
 * hole ringdown. However, this function is more general than the other function, which
 * only works for the (2,2) mode, and only the first three overtones.
 */
INT4 XLALGenerateQNMFreqV2(
        COMPLEX8Vector          *modefreqs,
        InspiralTemplate        *params,
        UINT4                    l,
        UINT4                    m,
        UINT4                   nmodes
        )
{

  /* Data for interpolating the quasinormal mode frequencies is taken from */
  /* The webpage of Emanuele Berti, http://www.phy.olemiss.edu/~berti/qnms.html */

  static const double afinallist[50] = {0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97, 0.975, 0.98, 0.982, 0.984, 0.986, 0.988, 0.99, 0.992, 0.994, 0.995, 0.996, 0.997, 0.9975, 0.998, 0.9985, 0.999, 0.9992, 0.9994, 0.9995, 0.9996, 0.9997, 0.9998, 0.9999};

  /* 2, 2 mode */

  const double reomegaqnm22[8][50] = {{0.373672, 0.380146, 0.387018, 0.394333, 0.402145, 0.410518, 0.419527, 0.429264, 0.439842, 0.451402, 0.464123, 0.478235, 0.494045, 0.511969, 0.5326, 0.556817, 0.586017, 0.59958, 0.614539, 0.631206, 0.650018, 0.671614, 0.696995, 0.727875, 0.74632, 0.767674, 0.793208, 0.808235, 0.825429, 0.8331, 0.841343, 0.850272, 0.860046, 0.870893, 0.883162, 0.897446, 0.905664, 0.914902, 0.925581, 0.931689, 0.938524, 0.946385, 0.955854, 0.960358, 0.965514, 0.968438, 0.97169, 0.975404, 0.979841, 0.985674}, 
{0.346711, 0.354101, 0.36191, 0.370183, 0.378976, 0.388353, 0.39839, 0.409183, 0.420847, 0.433527, 0.447407, 0.462728, 0.479807, 0.499079, 0.521161, 0.546964, 0.577922, 0.592247, 0.608001, 0.625499, 0.645173, 0.667658, 0.693938, 0.725708, 0.744582, 0.766349, 0.792272, 0.807482, 0.824852, 0.832591, 0.840901, 0.849896, 0.859735, 0.870645, 0.882977, 0.897322, 0.905568, 0.914834, 0.925538, 0.931657, 0.938502, 0.946371, 0.955847, 0.960353, 0.96551, 0.968435, 0.971688, 0.975402, 0.97984, 0.985673}, 
{0.301053, 0.309873, 0.319153, 0.328939, 0.339285, 0.350255, 0.361927, 0.374396, 0.387779, 0.402225, 0.417925, 0.43513, 0.454179, 0.475545, 0.499906, 0.528277, 0.56224, 0.577926, 0.595148, 0.614222, 0.635576, 0.659827, 0.687923, 0.721496, 0.741239, 0.763831, 0.790516, 0.806079, 0.823779, 0.831643, 0.840075, 0.849189, 0.859144, 0.870168, 0.882612, 0.897067, 0.905368, 0.914687, 0.925444, 0.931587, 0.938454, 0.946343, 0.955833, 0.960343, 0.965503, 0.96843, 0.971683, 0.975399, 0.979837, 0.985671}, 
{0.251505, 0.261348, 0.271749, 0.282749, 0.294397, 0.306753, 0.319886, 0.333884, 0.348856, 0.364938, 0.382309, 0.40121, 0.421977, 0.445099, 0.471336, 0.501925, 0.538956, 0.556285, 0.575464, 0.59683, 0.620788, 0.647869, 0.678894, 0.715333, 0.736426, 0.760279, 0.788101, 0.804173, 0.822338, 0.830374, 0.83897, 0.848242, 0.858347, 0.869515, 0.8821, 0.896696, 0.905068, 0.91446, 0.925289, 0.931468, 0.938371, 0.946292, 0.955811, 0.960329, 0.965496, 0.968424, 0.971679, 0.975395, 0.979834, 0.985669}, 
{0.207515, 0.217594, 0.228433, 0.24006, 0.252511, 0.265825, 0.280057, 0.295269, 0.311546, 0.328992, 0.347741, 0.367967, 0.389902, 0.413872, 0.440385, 0.470392, 0.506263, 0.523592, 0.543861, 0.568205, 0.597077, 0.629836, 0.666171, 0.707148, 0.730203, 0.755806, 0.785154, 0.801892, 0.82065, 0.828901, 0.837698, 0.847157, 0.857436, 0.868766, 0.881504, 0.896248, 0.904696, 0.914167, 0.925079, 0.931302, 0.938248, 0.946214, 0.955773, 0.960306, 0.965484, 0.968417, 0.971675, 0.975393, 0.979831, 0.985667}, 
{0.169299, 0.179085, 0.19014, 0.202441, 0.215947, 0.230623, 0.246457, 0.263464, 0.281691, 0.301218, 0.322157, 0.344655, 0.368896, 0.39509, 0.42345, 0.454055, 0.486283, 0.499037, 0.511101, 0.522095, 0.532461, 0.536917, 0.527548, 0.521465, 0.518438, 0.515417, 0.51241, 0.510911, 0.509415, 0.508818, 0.508221, 0.507624, 0.507028, 0.506432, 0.505837, 0.505242, 0.504944, 0.504647, 0.50435, 0.504202, 0.504054, 0.503904, 0.50376, 0.503692, 0.503638, 0.503632, 0.503585, 0.503474, 0.503674, 0.502985}, 
{0.133252, 0.142539, 0.154841, 0.169536, 0.185971, 0.203721, 0.222585, 0.242517, 0.263571, 0.285865, 0.309569, 0.334899, 0.362133, 0.391629, 0.423854, 0.459429, 0.499165, 0.516428, 0.534543, 0.553525, 0.573596, 0.602227, 0.649295, 0.697231, 0.722841, 0.750611, 0.781797, 0.799333, 0.818797, 0.827303, 0.836336, 0.846012, 0.856487, 0.867993, 0.880887, 0.895775, 0.904293, 0.913837, 0.92483, 0.931097, 0.938091, 0.946107, 0.955717, 0.960269, 0.965464, 0.968404, 0.971669, 0.975386, 0.979829, 0.985667}, 
{0.0928224, 0.102996, 0.122468, 0.144233, 0.166117, 0.187881, 0.209753, 0.232042, 0.255048, 0.279057, 0.304354, 0.331244, 0.360084, 0.391319, 0.425537, 0.463567, 0.506652, 0.52572, 0.546107, 0.568081, 0.592057, 0.618799, 0.650217, 0.691189, 0.716737, 0.745623, 0.778341, 0.796657, 0.816858, 0.825638, 0.834931, 0.844846, 0.855539, 0.867237, 0.880294, 0.895318, 0.903898, 0.913504, 0.924565, 0.930871, 0.93791, 0.945976, 0.955643, 0.960217, 0.965434, 0.968355, 0.971639, 0.975378, 0.979825, 0.985667}};

  const double imomegaqnm22[8][50] = {{0.0889623, 0.0888489, 0.0887057, 0.0885283, 0.0883112, 0.0880477, 0.0877293, 0.0873453, 0.086882, 0.0863212, 0.0856388, 0.0848021, 0.0837652, 0.0824618, 0.0807929, 0.0786025, 0.0756296, 0.0741258, 0.072378, 0.0703215, 0.0678642, 0.0648692, 0.0611186, 0.0562313, 0.053149, 0.0494336, 0.0447904, 0.0419586, 0.0386302, 0.0371155, 0.0354676, 0.033659, 0.0316516, 0.0293904, 0.0267908, 0.0237095, 0.0219107, 0.0198661, 0.0174737, 0.0160919, 0.014534, 0.0127274, 0.0105306, 0.00947799, 0.00826689, 0.00757704, 0.0068074, 0.00592534, 0.00486726, 0.00346901}, 
{0.273915, 0.273232, 0.272452, 0.271562, 0.270546, 0.269383, 0.268049, 0.266512, 0.264733, 0.262661, 0.260225, 0.257331, 0.253847, 0.24958, 0.244238, 0.237357, 0.228149, 0.223524, 0.218165, 0.211876, 0.204376, 0.195252, 0.183847, 0.169019, 0.159687, 0.148458, 0.134453, 0.125927, 0.115917, 0.111365, 0.106415, 0.100985, 0.0949599, 0.0881754, 0.080377, 0.0711342, 0.0657383, 0.0596046, 0.0524264, 0.04828, 0.043605, 0.0381837, 0.0315916, 0.0284334, 0.0247999, 0.0227304, 0.0204216, 0.0177757, 0.0146015, 0.0104061}, 
{0.478277, 0.475955, 0.473463, 0.470778, 0.467872, 0.464713, 0.46126, 0.457463, 0.453259, 0.448569, 0.443287, 0.437269, 0.430315, 0.422131, 0.412262, 0.399957, 0.383895, 0.375919, 0.366715, 0.35594, 0.34311, 0.327518, 0.308058, 0.282827, 0.266998, 0.248008, 0.22441, 0.210086, 0.193307, 0.18569, 0.177413, 0.16834, 0.158283, 0.146966, 0.133965, 0.118563, 0.109572, 0.0993515, 0.0873891, 0.0804781, 0.0726848, 0.0636465, 0.0526555, 0.0473899, 0.0413327, 0.037883, 0.0340348, 0.0296251, 0.0243356, 0.0242817}, 
{0.705148, 0.699946, 0.694481, 0.688721, 0.682628, 0.676155, 0.669243, 0.661823, 0.653808, 0.645091, 0.635536, 0.624969, 0.613154, 0.599764, 0.584301, 0.56591, 0.542888, 0.531648, 0.518699, 0.50348, 0.485224, 0.462864, 0.434823, 0.398473, 0.375741, 0.348572, 0.314974, 0.294663, 0.270945, 0.2602, 0.248542, 0.235779, 0.221649, 0.205769, 0.187547, 0.16598, 0.153396, 0.139094, 0.122354, 0.112682, 0.101773, 0.0891181, 0.0737264, 0.0663517, 0.0578682, 0.053037, 0.0476481, 0.0414736, 0.0340683, 0.0312184}, 
{0.946845, 0.937987, 0.928725, 0.919026, 0.908845, 0.898118, 0.886764, 0.874681, 0.861738, 0.847776, 0.832595, 0.815947, 0.797527, 0.776957, 0.7538, 0.727605, 0.697962, 0.684902, 0.670601, 0.65375, 0.632042, 0.603295, 0.56587, 0.517102, 0.486783, 0.450771, 0.406552, 0.379969, 0.349048, 0.33508, 0.319949, 0.303412, 0.285135, 0.264631, 0.241142, 0.213382, 0.197201, 0.178818, 0.157307, 0.144878, 0.130859, 0.114594, 0.0948041, 0.0853202, 0.0744092, 0.0681953, 0.0612642, 0.0533231, 0.0535335, 0.0381542}, 
{1.19561, 1.18221, 1.16825, 1.15377, 1.13878, 1.1232, 1.10691, 1.08976, 1.07153, 1.05194, 1.03067, 1.0073, 0.981265, 0.951825, 0.91795, 0.878201, 0.830803, 0.809534, 0.787531, 0.766572, 0.750871, 0.748723, 0.741744, 0.732374, 0.728044, 0.723765, 0.719565, 0.717492, 0.715438, 0.714621, 0.713807, 0.712996, 0.712188, 0.711383, 0.710581, 0.709781, 0.709382, 0.708984, 0.708587, 0.708388, 0.70819, 0.707992, 0.707796, 0.70772, 0.70762, 0.707592, 0.707598, 0.7075, 0.707413, 0.707316}, 
{1.44791, 1.42768, 1.40711, 1.38682, 1.36683, 1.34689, 1.32662, 1.30561, 1.28347, 1.25976, 1.234, 1.20558, 1.17376, 1.13753, 1.09541, 1.04518, 0.983084, 0.953463, 0.919986, 0.881204, 0.833927, 0.770337, 0.707988, 0.641106, 0.601525, 0.555418, 0.499613, 0.466355, 0.427874, 0.41055, 0.391822, 0.371395, 0.348864, 0.323638, 0.2948, 0.260789, 0.240992, 0.218519, 0.192238, 0.177057, 0.159933, 0.140065, 0.115884, 0.104293, 0.0909554, 0.0833586, 0.0748845, 0.0770313, 0.073006, 0.0589629}, 
{1.70384, 1.66882, 1.63893, 1.61423, 1.59206, 1.57064, 1.54892, 1.52619, 1.5019, 1.47555, 1.44657, 1.41428, 1.37781, 1.33598, 1.28712, 1.22872, 1.15671, 1.12263, 1.08453, 1.04141, 0.991742, 0.933178, 0.861906, 0.773444, 0.721837, 0.663505, 0.594647, 0.554166, 0.507666, 0.486825, 0.464348, 0.439888, 0.412969, 0.382898, 0.348601, 0.308246, 0.284798, 0.25821, 0.227147, 0.209213, 0.188989, 0.165524, 0.136961, 0.123266, 0.107504, 0.113694, 0.102135, 0.08889, 0.0827456, 0.0728406}};

  /* 2, 1 mode */

  const double reomegaqnm21[8][50] = {{0.373672, 0.376931, 0.380432, 0.384197, 0.388248, 0.392615, 0.39733, 0.402436, 0.407979, 0.41402, 0.420632, 0.427909, 0.435968, 0.444968, 0.455121, 0.466727, 0.480231, 0.486308, 0.492859, 0.499965, 0.507729, 0.516291, 0.525845, 0.536673, 0.542693, 0.549213, 0.556329, 0.560146, 0.564155, 0.565814, 0.567505, 0.569227, 0.570976, 0.572749, 0.574535, 0.576322, 0.577208, 0.578084, 0.578948, 0.579374, 0.579795, 0.580212, 0.580623, 0.580784, 0.580942, 0.581018, 0.581093, 0.581171, 0.581283, 0.581826}, 
{0.346711, 0.350448, 0.354478, 0.358823, 0.363506, 0.368558, 0.374014, 0.379916, 0.386313, 0.393268, 0.400856, 0.409172, 0.418337, 0.428509, 0.439897, 0.452793, 0.467612, 0.474208, 0.481261, 0.488834, 0.497004, 0.505862, 0.515506, 0.526011, 0.531572, 0.537263, 0.542874, 0.545479, 0.547723, 0.548438, 0.548984, 0.54929, 0.549241, 0.548642, 0.54714, 0.544185, 0.542606, 0.542325, 0.541626, 0.540962, 0.540617, 0.540156, 0.539726, 0.539553, 0.539377, 0.539288, 0.539198, 0.539103, 0.539007, 0.539097}, 
{0.301053, 0.30554, 0.31041, 0.315684, 0.321389, 0.327554, 0.334218, 0.341421, 0.349217, 0.357668, 0.366852, 0.376864, 0.387825, 0.39989, 0.413261, 0.428206, 0.445088, 0.452488, 0.46031, 0.468593, 0.477369, 0.48665, 0.49639, 0.506348, 0.511184, 0.515607, 0.51903, 0.519992, 0.5201, 0.51983, 0.519361, 0.51871, 0.517943, 0.517204, 0.516728, 0.516696, 0.516098, 0.513376, 0.509813, 0.508083, 0.506412, 0.504789, 0.503194, 0.502559, 0.501923, 0.501605, 0.501286, 0.500966, 0.500947, 0.500626}, 
{0.251505, 0.256534, 0.262058, 0.268101, 0.274689, 0.281852, 0.289628, 0.298059, 0.307197, 0.317104, 0.327855, 0.339542, 0.352278, 0.366207, 0.38151, 0.398418, 0.417227, 0.425361, 0.433878, 0.442797, 0.452127, 0.461855, 0.471918, 0.482195, 0.48739, 0.49271, 0.498495, 0.501813, 0.5056, 0.507251, 0.508939, 0.510585, 0.512031, 0.513028, 0.513289, 0.512758, 0.512367, 0.512091, 0.511921, 0.51135, 0.5096, 0.50715, 0.504694, 0.503744, 0.502809, 0.502344, 0.50188, 0.501414, 0.501243, 0.500775}, 
{0.207515, 0.212679, 0.21848, 0.224939, 0.232083, 0.23994, 0.248543, 0.257928, 0.268139, 0.279228, 0.291258, 0.304303, 0.318456, 0.33383, 0.350563, 0.368831, 0.388855, 0.397422, 0.406346, 0.41567, 0.425462, 0.435857, 0.447149, 0.45999, 0.467349, 0.475546, 0.484594, 0.489367, 0.494295, 0.496339, 0.498459, 0.500692, 0.50308, 0.505635, 0.508214, 0.510197, 0.510572, 0.510387, 0.509864, 0.509635, 0.509456, 0.50879, 0.506257, 0.504974, 0.503704, 0.503082, 0.502466, 0.501855, 0.501243, 0.500922}, 
{0.169299, 0.174306, 0.180192, 0.186977, 0.194671, 0.203282, 0.212811, 0.223263, 0.234649, 0.246985, 0.260297, 0.274622, 0.290009, 0.306518, 0.324226, 0.343234, 0.363727, 0.372432, 0.381518, 0.391111, 0.401435, 0.412868, 0.425995, 0.441563, 0.45052, 0.460463, 0.471717, 0.477974, 0.484665, 0.487438, 0.490248, 0.493085, 0.49596, 0.498923, 0.502084, 0.505527, 0.507222, 0.508526, 0.508875, 0.508627, 0.508251, 0.507926, 0.507181, 0.506113, 0.504624, 0.503839, 0.50306, 0.502293, 0.501826, 0.501211}, 
{0.133252, 0.137863, 0.144021, 0.151672, 0.160684, 0.170902, 0.182186, 0.194429, 0.207556, 0.221528, 0.236327, 0.251955, 0.268421, 0.285744, 0.30395, 0.323112, 0.343522, 0.352232, 0.361443, 0.371382, 0.382361, 0.394771, 0.409089, 0.426092, 0.436065, 0.4474, 0.460415, 0.467675, 0.475563, 0.478922, 0.482402, 0.485993, 0.489665, 0.493377, 0.497113, 0.500995, 0.503086, 0.505273, 0.507186, 0.507677, 0.507628, 0.507177, 0.506714, 0.506383, 0.505391, 0.504576, 0.503665, 0.502738, 0.502119, 0.501502}, 
{0.0928224, 0.0967608, 0.105267, 0.116748, 0.129817, 0.143719, 0.158124, 0.17292, 0.188094, 0.203673, 0.219697, 0.236201, 0.253208, 0.270727, 0.288784, 0.307531, 0.327628, 0.336416, 0.345907, 0.35633, 0.367923, 0.380968, 0.395965, 0.413866, 0.424377, 0.436366, 0.450442, 0.458506, 0.467375, 0.471171, 0.475126, 0.479254, 0.483568, 0.488052, 0.492639, 0.497232, 0.499557, 0.501993, 0.504598, 0.505838, 0.506697, 0.506742, 0.506156, 0.505924, 0.505541, 0.505047, 0.504219, 0.503188, 0.502416, 0.501648}};

  const double imomegaqnm21[8][50] = {{0.0889623, 0.0888968, 0.0887983, 0.0886636, 0.0884885, 0.0882679, 0.0879952, 0.0876618, 0.0872571, 0.086767, 0.086173, 0.0854501, 0.0845642, 0.0834665, 0.0820852, 0.0803094, 0.077955, 0.0767847, 0.0754394, 0.0738744, 0.072027, 0.0698043, 0.067061, 0.0635492, 0.0613721, 0.058791, 0.0556436, 0.0537759, 0.0516427, 0.0506978, 0.0496908, 0.0486136, 0.0474565, 0.0462084, 0.0448564, 0.0433879, 0.0426068, 0.0417939, 0.0409501, 0.0405173, 0.0400777, 0.0396319, 0.0391802, 0.0389974, 0.0388111, 0.038714, 0.0386091, 0.0384831, 0.0382843, 0.0376257}, 
{0.273915, 0.273543, 0.273058, 0.272452, 0.271711, 0.270819, 0.269754, 0.268489, 0.266992, 0.265218, 0.263108, 0.260587, 0.257547, 0.253838, 0.249238, 0.243404, 0.235766, 0.231999, 0.227686, 0.222686, 0.2168, 0.209731, 0.201006, 0.189803, 0.18282, 0.174484, 0.164205, 0.158029, 0.150891, 0.1477, 0.14428, 0.140606, 0.136657, 0.132444, 0.128112, 0.124511, 0.123925, 0.123321, 0.12167, 0.121162, 0.120732, 0.12014, 0.119631, 0.119429, 0.119227, 0.119127, 0.119027, 0.118922, 0.118789, 0.118488}, 
{0.478277, 0.477047, 0.475592, 0.473893, 0.471929, 0.469669, 0.467076, 0.464104, 0.460691, 0.456761, 0.45221, 0.446903, 0.440653, 0.433196, 0.424144, 0.412892, 0.398434, 0.391383, 0.383355, 0.374094, 0.363227, 0.350188, 0.334044, 0.313093, 0.299808, 0.283583, 0.262719, 0.249472, 0.233119, 0.225306, 0.216483, 0.206365, 0.194563, 0.180538, 0.16345, 0.141483, 0.12731, 0.111596, 0.0954876, 0.0868316, 0.0774193, 0.066857, 0.0544355, 0.0486318, 0.0420655, 0.0383765, 0.0343031, 0.029688, 0.0342671, 0.0312969}, 
{0.705148, 0.702423, 0.699317, 0.695803, 0.691851, 0.687416, 0.682445, 0.676867, 0.670592, 0.663504, 0.655451, 0.646227, 0.635554, 0.623036, 0.60809, 0.589806, 0.566651, 0.555455, 0.54276, 0.528156, 0.511038, 0.490465, 0.464847, 0.431189, 0.409563, 0.382936, 0.348824, 0.32773, 0.302947, 0.291759, 0.279694, 0.26659, 0.252176, 0.235931, 0.216835, 0.192946, 0.178198, 0.160681, 0.138951, 0.125667, 0.110911, 0.095128, 0.0772103, 0.0689279, 0.0595841, 0.0543422, 0.048559, 0.0420118, 0.0443244, 0.0383948}, 
{0.946845, 0.942236, 0.937057, 0.931283, 0.924873, 0.917775, 0.909918, 0.901208, 0.891524, 0.880707, 0.868545, 0.854753, 0.838939, 0.820547, 0.798756, 0.772277, 0.738922, 0.722836, 0.704612, 0.683652, 0.659084, 0.629562, 0.592913, 0.54542, 0.515651, 0.480063, 0.436051, 0.409334, 0.377864, 0.363486, 0.347803, 0.330563, 0.311454, 0.290066, 0.265802, 0.237412, 0.220683, 0.201065, 0.17693, 0.162373, 0.145371, 0.124789, 0.100358, 0.0893955, 0.0771778, 0.0703585, 0.0628494, 0.054359, 0.0443244, 0.0454973}, 
{1.19561, 1.18868, 1.18098, 1.17251, 1.16324, 1.15314, 1.14213, 1.13009, 1.11686, 1.10224, 1.08592, 1.06753, 1.04652, 1.02214, 0.993243, 0.958068, 0.913601, 0.892095, 0.867699, 0.839644, 0.806856, 0.767797, 0.720146, 0.659944, 0.622848, 0.578732, 0.524237, 0.491317, 0.452927, 0.435542, 0.416661, 0.395943, 0.372912, 0.346898, 0.316949, 0.28167, 0.261342, 0.238436, 0.211307, 0.195088, 0.176156, 0.153337, 0.124273, 0.110319, 0.0949385, 0.0864639, 0.0771851, 0.0667293, 0.0644759, 0.059715}, 
{1.44791, 1.43752, 1.42619, 1.41412, 1.4014, 1.38801, 1.37384, 1.35872, 1.34237, 1.32445, 1.30456, 1.28213, 1.25643, 1.22645, 1.19065, 1.14668, 1.09053, 1.06323, 1.03225, 0.996712, 0.955487, 0.906959, 0.848541, 0.775494, 0.730792, 0.678016, 0.613349, 0.574379, 0.528885, 0.508284, 0.485946, 0.461517, 0.434493, 0.404108, 0.369135, 0.327547, 0.303278, 0.275875, 0.244209, 0.225953, 0.204982, 0.179727, 0.147552, 0.131555, 0.113035, 0.102763, 0.09161, 0.0791341, 0.0745748, 0.0739536}, 
{1.70384, 1.68568, 1.6674, 1.65022, 1.63402, 1.61822, 1.60219, 1.58534, 1.56716, 1.54712, 1.52462, 1.49894, 1.46917, 1.43397, 1.39142, 1.33846, 1.27014, 1.23683, 1.19912, 1.15609, 1.10649, 1.04841, 0.97878, 0.892317, 0.839828, 0.778168, 0.702974, 0.657884, 0.605411, 0.581667, 0.555906, 0.527711, 0.496517, 0.46151, 0.421398, 0.373867, 0.346042, 0.314367, 0.277467, 0.25637, 0.232725, 0.204918, 0.169568, 0.151985, 0.131247, 0.119312, 0.106197, 0.0916018, 0.0847009, 0.0810854}};

  /* 3, 3 mode */
  const double reomegaqnm33[8][50] = {{0.599443, 0.609823, 0.620796, 0.632425, 0.644787, 0.657972, 0.672086, 0.68726, 0.70365, 0.721455, 0.740921, 0.762369, 0.786223, 0.813057, 0.843687, 0.879323, 0.921885, 0.941521, 0.963088, 0.987016, 1.01391, 1.04464, 1.08058, 1.1241, 1.14998, 1.17986, 1.21547, 1.23637, 1.26023, 1.27086, 1.28227, 1.29462, 1.30812, 1.32308, 1.33999, 1.35965, 1.37094, 1.38363, 1.39829, 1.40666, 1.41603, 1.42679, 1.43975, 1.44591, 1.45295, 1.45695, 1.46139, 1.46646, 1.47252, 1.48047}, 
{0.582644, 0.593642, 0.60525, 0.617535, 0.630573, 0.644453, 0.659285, 0.675198, 0.692352, 0.710943, 0.731221, 0.753507, 0.778225, 0.805952, 0.837504, 0.874094, 0.917645, 0.937687, 0.959667, 0.984015, 1.01133, 1.0425, 1.07889, 1.12285, 1.14897, 1.17907, 1.21491, 1.23591, 1.25988, 1.27055, 1.28201, 1.2944, 1.30794, 1.32294, 1.33989, 1.35958, 1.37089, 1.38359, 1.39826, 1.40664, 1.41601, 1.42678, 1.43974, 1.4459, 1.45295, 1.45695, 1.46139, 1.46646, 1.47251, 1.48047}, 
{0.551685, 0.563804, 0.576567, 0.590041, 0.604302, 0.619445, 0.635578, 0.652833, 0.671371, 0.691391, 0.713145, 0.736956, 0.763249, 0.792606, 0.825846, 0.864187, 0.90956, 0.930355, 0.953102, 0.978234, 1.00636, 1.03835, 1.07558, 1.12042, 1.14697, 1.17752, 1.21379, 1.23502, 1.2592, 1.26995, 1.28149, 1.29396, 1.30758, 1.32266, 1.33967, 1.35943, 1.37078, 1.38351, 1.3982, 1.4066, 1.41598, 1.42676, 1.43973, 1.44589, 1.45294, 1.45694, 1.46139, 1.46646, 1.47251, 1.48047}, 
{0.511962, 0.525402, 0.539534, 0.554428, 0.570163, 0.586836, 0.604556, 0.62346, 0.643711, 0.665511, 0.689115, 0.71485, 0.743145, 0.774586, 0.809998, 0.850607, 0.898362, 0.920146, 0.943908, 0.970086, 0.999286, 1.0324, 1.0708, 1.11686, 1.14404, 1.17523, 1.21215, 1.23369, 1.25818, 1.26906, 1.28072, 1.2933, 1.30704, 1.32223, 1.33936, 1.35922, 1.37061, 1.38339, 1.39812, 1.40653, 1.41593, 1.42673, 1.43971, 1.44588, 1.45294, 1.45694, 1.46138, 1.46645, 1.47251, 1.48047}, 
{0.470174, 0.484725, 0.500042, 0.516195, 0.533269, 0.551359, 0.570583, 0.591077, 0.61301, 0.636588, 0.662069, 0.689786, 0.720173, 0.753822, 0.791564, 0.834636, 0.885005, 0.907885, 0.932779, 0.960131, 0.990559, 1.02497, 1.06475, 1.11229, 1.14026, 1.17226, 1.21, 1.23195, 1.25686, 1.26789, 1.2797, 1.29245, 1.30633, 1.32167, 1.33893, 1.35893, 1.37038, 1.38322, 1.39801, 1.40645, 1.41587, 1.42669, 1.43969, 1.44587, 1.45293, 1.45693, 1.46137, 1.46645, 1.47251, 1.48047}, 
{0.431386, 0.446612, 0.462696, 0.479713, 0.49775, 0.516906, 0.537303, 0.559082, 0.582416, 0.607517, 0.634649, 0.66415, 0.696462, 0.732182, 0.77215, 0.81761, 0.870543, 0.894505, 0.920521, 0.949042, 0.980704, 1.01643, 1.05767, 1.10684, 1.1357, 1.16864, 1.20737, 1.22983, 1.25522, 1.26646, 1.27846, 1.29139, 1.30546, 1.32098, 1.33842, 1.35857, 1.3701, 1.38301, 1.39787, 1.40634, 1.41579, 1.42664, 1.43966, 1.44585, 1.45291, 1.45692, 1.46137, 1.46644, 1.4725, 1.48047}, 
{0.39766, 0.413152, 0.429604, 0.447094, 0.465712, 0.485563, 0.506772, 0.529486, 0.553883, 0.580183, 0.608655, 0.639647, 0.673609, 0.711146, 0.753105, 0.800733, 0.856013, 0.880965, 0.908002, 0.937583, 0.970359, 1.00729, 1.04987, 1.10064, 1.13044, 1.16442, 1.20427, 1.22731, 1.2533, 1.26476, 1.27699, 1.29015, 1.30444, 1.32016, 1.3378, 1.35815, 1.36977, 1.38276, 1.3977, 1.40621, 1.4157, 1.42658, 1.43963, 1.44582, 1.4529, 1.45691, 1.46135, 1.46643, 1.4725, 1.48046}, 
{0.368992, 0.384475, 0.401019, 0.418708, 0.437635, 0.45791, 0.479661, 0.50304, 0.528232, 0.555461, 0.585008, 0.617226, 0.652577, 0.691676, 0.735379, 0.784937, 0.842317, 0.868149, 0.896082, 0.926573, 0.960272, 0.998162, 1.04182, 1.09394, 1.1246, 1.15963, 1.20071, 1.22441, 1.25107, 1.2628, 1.2753, 1.28872, 1.30326, 1.31923, 1.3371, 1.35766, 1.36938, 1.38248, 1.39751, 1.40606, 1.41559, 1.4265, 1.43959, 1.44579, 1.45288, 1.45689, 1.46134, 1.46643, 1.4725, 1.48046}};

  const double imomegaqnm33[8][50] = {{0.092703, 0.0925869, 0.0924305, 0.0922281, 0.0919726, 0.0916556, 0.0912666, 0.0907928, 0.0902179, 0.0895213, 0.0886763, 0.087647, 0.0863849, 0.0848213, 0.0828557, 0.0803327, 0.0769953, 0.075339, 0.0734361, 0.0712234, 0.0686106, 0.0654629, 0.0615646, 0.056537, 0.0533881, 0.0496087, 0.0449046, 0.0420439, 0.0386881, 0.0371631, 0.0355053, 0.0336873, 0.0316714, 0.0294027, 0.0267968, 0.023711, 0.0219107, 0.0198652, 0.0174725, 0.0160908, 0.014533, 0.0127266, 0.0105299, 0.00947732, 0.00826624, 0.0075764, 0.00680679, 0.00592479, 0.0048668, 0.00346854}, 
{0.281298, 0.280807, 0.280191, 0.279434, 0.278515, 0.277407, 0.27608, 0.274494, 0.272602, 0.27034, 0.26763, 0.264363, 0.260394, 0.255518, 0.249435, 0.241682, 0.231489, 0.226451, 0.220675, 0.213971, 0.206071, 0.19657, 0.184822, 0.169692, 0.160224, 0.148867, 0.134739, 0.12615, 0.116076, 0.111499, 0.106524, 0.101069, 0.0950193, 0.088212, 0.0803933, 0.0711348, 0.0657335, 0.0595967, 0.0524182, 0.0482728, 0.0435993, 0.03818, 0.0315898, 0.028432, 0.0247987, 0.0227292, 0.0204204, 0.0177744, 0.0146004, 0.0104057}, 
{0.479093, 0.477786, 0.476267, 0.474506, 0.472466, 0.470105, 0.467369, 0.464195, 0.4605, 0.456181, 0.451103, 0.44509, 0.437899, 0.42919, 0.418467, 0.404962, 0.387407, 0.37879, 0.36895, 0.357573, 0.344212, 0.328195, 0.308447, 0.283082, 0.267235, 0.248247, 0.224647, 0.210309, 0.193501, 0.185866, 0.177568, 0.168469, 0.158382, 0.147032, 0.133998, 0.118564, 0.10956, 0.0993313, 0.0873659, 0.0804562, 0.0726665, 0.0636339, 0.0526499, 0.0473868, 0.0413313, 0.0378821, 0.034034, 0.029624, 0.024334, 0.0173428}, 
{0.690337, 0.68755, 0.684459, 0.68102, 0.677183, 0.672887, 0.668056, 0.662599, 0.656401, 0.649315, 0.641153, 0.631664, 0.620511, 0.607221, 0.591101, 0.571089, 0.545428, 0.532943, 0.518759, 0.502439, 0.48336, 0.460589, 0.432625, 0.396831, 0.37452, 0.347819, 0.314676, 0.294558, 0.270986, 0.260282, 0.248651, 0.235902, 0.221769, 0.20587, 0.187614, 0.166001, 0.153393, 0.13907, 0.122317, 0.112642, 0.101735, 0.0890889, 0.0737103, 0.0663418, 0.057864, 0.053035, 0.0476477, 0.0414737, 0.0340677, 0.0242799}, 
{0.915649, 0.910755, 0.905456, 0.899695, 0.893406, 0.88651, 0.878908, 0.870479, 0.861073, 0.850501, 0.838516, 0.824795, 0.808901, 0.790225, 0.767883, 0.740519, 0.7059, 0.689212, 0.670355, 0.648776, 0.623682, 0.593881, 0.557455, 0.511016, 0.482143, 0.44764, 0.404868, 0.378929, 0.348555, 0.334768, 0.319791, 0.303378, 0.285189, 0.264731, 0.241246, 0.213447, 0.197233, 0.178815, 0.157271, 0.144831, 0.130807, 0.114545, 0.0947715, 0.0852972, 0.0743968, 0.068188, 0.0612614, 0.0533233, 0.0438014, 0.031217}, 
{1.15215, 1.14481, 1.13695, 1.12848, 1.11934, 1.10942, 1.09859, 1.08672, 1.07361, 1.05902, 1.04265, 1.02409, 1.00282, 0.978074, 0.948772, 0.913265, 0.868853, 0.847622, 0.823753, 0.796583, 0.765161, 0.728048, 0.682917, 0.625639, 0.590118, 0.547732, 0.49525, 0.463448, 0.42623, 0.409344, 0.391005, 0.370913, 0.348652, 0.323623, 0.294897, 0.260904, 0.241081, 0.218564, 0.192229, 0.177023, 0.159881, 0.140004, 0.115834, 0.104253, 0.09093, 0.0833412, 0.0748751, 0.0770227, 0.0632687, 0.0589657}, 
{1.39591, 1.38607, 1.37555, 1.36428, 1.35215, 1.33904, 1.32481, 1.30928, 1.29221, 1.27333, 1.25226, 1.22852, 1.20145, 1.17016, 1.13335, 1.08907, 1.03413, 1.00803, 0.978819, 0.945725, 0.907651, 0.862934, 0.808869, 0.7406, 0.698384, 0.648072, 0.585828, 0.54813, 0.50403, 0.484027, 0.462309, 0.43852, 0.412171, 0.382553, 0.348572, 0.308374, 0.284937, 0.258319, 0.227191, 0.209218, 0.188957, 0.165464, 0.136898, 0.123211, 0.107464, 0.0984948, 0.102103, 0.0888725, 0.0730023, 0.0659028}, 
{1.64384, 1.63157, 1.61846, 1.60441, 1.58931, 1.57303, 1.55539, 1.53617, 1.5151, 1.49185, 1.46597, 1.43689, 1.40384, 1.36575, 1.32111, 1.2676, 1.20156, 1.17032, 1.13545, 1.09609, 1.05101, 0.998352, 0.935078, 0.85569, 0.806777, 0.748562, 0.676569, 0.632966, 0.58196, 0.558828, 0.533715, 0.506213, 0.475757, 0.441533, 0.40228, 0.355861, 0.328804, 0.29808, 0.262156, 0.241415, 0.218035, 0.190926, 0.157962, 0.142169, 0.123998, 0.113649, 0.115718, 0.100723, 0.0827359, 0.0728399}};


  /* 4, 4 mode */
  const double reomegaqnm44[8][50] = {{0.809178, 0.823517, 0.83866, 0.854693, 0.871718, 0.889853, 0.909242, 0.930054, 0.9525, 0.976839, 1.0034, 1.03259, 1.06498, 1.10131, 1.14265, 1.19057, 1.24755, 1.27374, 1.30245, 1.33422, 1.36984, 1.41042, 1.45773, 1.51478, 1.54862, 1.58759, 1.6339, 1.66102, 1.69194, 1.7057, 1.72045, 1.73641, 1.75384, 1.77314, 1.79492, 1.82022, 1.83474, 1.85104, 1.86985, 1.8806, 1.8926, 1.9064, 1.92299, 1.93088, 1.93989, 1.945, 1.95068, 1.95716, 1.9649, 1.97507}, 
{0.796632, 0.811434, 0.827055, 0.843581, 0.861114, 0.879773, 0.899703, 0.921074, 0.944098, 0.969034, 0.996209, 1.02604, 1.05909, 1.09611, 1.13816, 1.1868, 1.24453, 1.27103, 1.30005, 1.33214, 1.36806, 1.40896, 1.45658, 1.51396, 1.54795, 1.58707, 1.63353, 1.66073, 1.69172, 1.7055, 1.72028, 1.73626, 1.75372, 1.77304, 1.79485, 1.82017, 1.8347, 1.85101, 1.86983, 1.88058, 1.89259, 1.90639, 1.92299, 1.93087, 1.93989, 1.945, 1.95068, 1.95716, 1.9649, 1.97507}, 
{0.77271, 0.788392, 0.804918, 0.822375, 0.840868, 0.860518, 0.88147, 0.903897, 0.928012, 0.954076, 0.982419, 1.01346, 1.04777, 1.08608, 1.12947, 1.17951, 1.23867, 1.26576, 1.29538, 1.32806, 1.36459, 1.4061, 1.45434, 1.51232, 1.54662, 1.58604, 1.6328, 1.66014, 1.69127, 1.70511, 1.71994, 1.73598, 1.75348, 1.77286, 1.79471, 1.82008, 1.83463, 1.85096, 1.8698, 1.88055, 1.89257, 1.90638, 1.92298, 1.93087, 1.93989, 1.945, 1.95068, 1.95716, 1.9649, 1.97507}, 
{0.739837, 0.75669, 0.774421, 0.793122, 0.812899, 0.833875, 0.856197, 0.880043, 0.905628, 0.933216, 0.96314, 0.995828, 1.03184, 1.07193, 1.11717, 1.16913, 1.2303, 1.25821, 1.28866, 1.32219, 1.35958, 1.40195, 1.45106, 1.50993, 1.54467, 1.58453, 1.63172, 1.65927, 1.69061, 1.70453, 1.71944, 1.73555, 1.75313, 1.77258, 1.79451, 1.81994, 1.83452, 1.85088, 1.86975, 1.88051, 1.89254, 1.90636, 1.92297, 1.93086, 1.93988, 1.94499, 1.95068, 1.95716, 1.9649, 1.97507}, 
{0.701516, 0.71962, 0.738651, 0.758702, 0.779882, 0.802317, 0.82616, 0.851591, 0.878829, 0.908143, 0.939873, 0.974449, 1.01244, 1.05461, 1.10202, 1.15627, 1.21984, 1.24875, 1.28021, 1.31477, 1.35322, 1.39666, 1.44687, 1.50684, 1.54214, 1.58256, 1.63031, 1.65813, 1.68974, 1.70376, 1.71878, 1.73499, 1.75267, 1.77222, 1.79424, 1.81975, 1.83438, 1.85077, 1.86967, 1.88046, 1.8925, 1.90633, 1.92296, 1.93085, 1.93987, 1.94499, 1.95067, 1.95716, 1.9649, 1.97507}, 
{0.661572, 0.680777, 0.700971, 0.722252, 0.744732, 0.768542, 0.793839, 0.820807, 0.849672, 0.880708, 0.914263, 0.950774, 0.990817, 1.03516, 1.0849, 1.14162, 1.20782, 1.23782, 1.27041, 1.30613, 1.34576, 1.39042, 1.44188, 1.50314, 1.53909, 1.58017, 1.62858, 1.65673, 1.68866, 1.70281, 1.71796, 1.7343, 1.7521, 1.77176, 1.7939, 1.81952, 1.83419, 1.85064, 1.86958, 1.88039, 1.89245, 1.9063, 1.92294, 1.93083, 1.93986, 1.94498, 1.95067, 1.95715, 1.9649, 1.97507}, 
{0.623109, 0.643109, 0.664176, 0.686412, 0.709933, 0.734875, 0.761399, 0.789696, 0.819999, 0.85259, 0.887824, 0.92615, 0.968153, 1.01462, 1.06665, 1.12585, 1.19475, 1.22589, 1.25966, 1.29659, 1.33748, 1.38345, 1.43625, 1.49892, 1.53559, 1.57741, 1.62657, 1.6551, 1.6874, 1.7017, 1.717, 1.73348, 1.75142, 1.77123, 1.79349, 1.81924, 1.83398, 1.85048, 1.86948, 1.8803, 1.89239, 1.90626, 1.92291, 1.93082, 1.93985, 1.94497, 1.95066, 1.95715, 1.9649, 1.97506},
{0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72}};

  const double imomegaqnm44[8][50] = {{0.094164, 0.0940768, 0.0939478, 0.0937705, 0.0935374, 0.0932393, 0.092865, 0.0924006, 0.0918287, 0.0911274, 0.0902679, 0.0892124, 0.0879094, 0.086287, 0.0842401, 0.0816078, 0.0781255, 0.0763992, 0.0744185, 0.0721192, 0.0694103, 0.066156, 0.0621396, 0.056982, 0.0537633, 0.0499111, 0.0451315, 0.0422323, 0.0388376, 0.037297, 0.0356236, 0.0337901, 0.0317587, 0.0294747, 0.0268536, 0.0237529, 0.0219453, 0.0198926, 0.0174928, 0.0161075, 0.0145463, 0.0127365, 0.0105364, 0.00948248, 0.00827008, 0.00757959, 0.00680933, 0.00592669, 0.00486806, 0.00346918}, 
{0.284334, 0.28399, 0.283519, 0.282901, 0.282113, 0.281129, 0.279913, 0.278426, 0.276614, 0.274412, 0.271733, 0.268465, 0.264452, 0.25948, 0.253233, 0.24523, 0.234678, 0.229458, 0.223476, 0.216539, 0.208375, 0.198576, 0.186494, 0.170992, 0.161324, 0.149756, 0.135408, 0.126707, 0.116519, 0.111896, 0.106875, 0.101374, 0.0952787, 0.0884259, 0.0805622, 0.0712595, 0.0658364, 0.0596782, 0.0524786, 0.0483228, 0.043639, 0.0382095, 0.0316092, 0.0284475, 0.0248103, 0.0227388, 0.020428, 0.0177801, 0.0146042, 0.0104075}, 
{0.479908, 0.479056, 0.477986, 0.476666, 0.47506, 0.47312, 0.470791, 0.468003, 0.464669, 0.460679, 0.45589, 0.450113, 0.44309, 0.434465, 0.423712, 0.410032, 0.392108, 0.383276, 0.373175, 0.361487, 0.347759, 0.331313, 0.311071, 0.28514, 0.268985, 0.249669, 0.225726, 0.21121, 0.194221, 0.186512, 0.17814, 0.168968, 0.158807, 0.147383, 0.134275, 0.118768, 0.109729, 0.099465, 0.0874652, 0.0805385, 0.072732, 0.0636828, 0.0526822, 0.0474125, 0.0413505, 0.037898, 0.0340467, 0.0296335, 0.0243403, 0.0173459}, 
{0.683924, 0.682151, 0.680068, 0.677629, 0.674782, 0.671462, 0.667589, 0.663066, 0.657772, 0.651551, 0.644203, 0.635463, 0.624972, 0.612229, 0.596503, 0.576677, 0.550915, 0.538287, 0.523888, 0.507273, 0.48781, 0.464555, 0.436002, 0.399509, 0.376808, 0.349693, 0.316109, 0.295761, 0.271954, 0.261154, 0.249426, 0.236578, 0.222347, 0.206349, 0.187993, 0.166281, 0.153625, 0.139254, 0.122453, 0.112755, 0.101826, 0.0891563, 0.0737552, 0.0663776, 0.0578907, 0.0530572, 0.0476654, 0.0414869, 0.0340764, 0.0242843}, 
{0.898239, 0.895043, 0.891442, 0.887378, 0.88278, 0.877566, 0.871633, 0.864855, 0.857075, 0.848091, 0.837647, 0.8254, 0.81089, 0.79347, 0.772205, 0.745659, 0.711484, 0.694831, 0.675907, 0.654142, 0.628729, 0.598457, 0.5614, 0.514171, 0.484847, 0.449863, 0.406581, 0.380377, 0.349729, 0.335831, 0.32074, 0.304211, 0.285904, 0.265327, 0.24172, 0.213799, 0.197524, 0.179045, 0.157442, 0.144973, 0.13092, 0.11463, 0.0948285, 0.0853429, 0.0744311, 0.0682165, 0.0612841, 0.0533403, 0.0438126, 0.038161}, 
{1.12298, 1.1179, 1.11232, 1.10615, 1.09933, 1.09174, 1.08325, 1.07371, 1.06294, 1.05067, 1.0366, 1.0203, 1.00121, 0.97854, 0.951144, 0.917273, 0.874063, 0.853135, 0.829435, 0.802272, 0.770663, 0.733138, 0.687353, 0.629184, 0.593147, 0.550213, 0.497162, 0.465071, 0.427558, 0.410551, 0.392089, 0.371872, 0.349482, 0.32432, 0.295457, 0.261323, 0.241428, 0.218839, 0.192433, 0.177192, 0.160015, 0.140105, 0.115902, 0.104308, 0.0909715, 0.0833759, 0.0749029, 0.0651937, 0.0535487, 0.0450994}, 
{1.35669, 1.34943, 1.34155, 1.33296, 1.32356, 1.31324, 1.30183, 1.28915, 1.27497, 1.259, 1.24085, 1.22004, 1.19589, 1.16745, 1.13339, 1.09162, 1.03877, 1.01331, 0.984584, 0.951766, 0.913705, 0.868675, 0.81392, 0.744589, 0.701738, 0.650764, 0.587868, 0.549856, 0.505448, 0.485323, 0.46348, 0.439564, 0.413085, 0.383331, 0.349205, 0.308853, 0.285337, 0.258637, 0.227427, 0.209413, 0.189112, 0.16558, 0.136976, 0.123274, 0.107512, 0.0985354, 0.0885218, 0.0770472, 0.0632849, 0.0659145}, 
{0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28}};

  /* 5, 5 mode */
  const double reomegaqnm55[8][50] = {{1.0123, 1.03056, 1.04985, 1.07027, 1.09194, 1.11502, 1.13968, 1.16614, 1.19467, 1.22558, 1.25928, 1.29631, 1.33734, 1.38332, 1.43555, 1.496, 1.56774, 1.60067, 1.63671, 1.67655, 1.72115, 1.77188, 1.83092, 1.90197, 1.94403, 1.99239, 2.04977, 2.08333, 2.12154, 2.13852, 2.15673, 2.17641, 2.1979, 2.22168, 2.24849, 2.27961, 2.29746, 2.31749, 2.34058, 2.35377, 2.3685, 2.38542, 2.40575, 2.41541, 2.42645, 2.43271, 2.43967, 2.4476, 2.45707, 2.46951}, 
{1.00222, 1.02086, 1.04054, 1.06135, 1.08343, 1.10693, 1.13203, 1.15894, 1.18793, 1.21932, 1.25352, 1.29105, 1.33262, 1.37915, 1.43196, 1.493, 1.56534, 1.59852, 1.63481, 1.67491, 1.71976, 1.77074, 1.83004, 1.90134, 1.94352, 1.992, 2.0495, 2.08311, 2.12137, 2.13838, 2.15661, 2.17631, 2.19781, 2.22161, 2.24844, 2.27957, 2.29743, 2.31747, 2.34057, 2.35376, 2.36849, 2.38541, 2.40575, 2.41541, 2.42645, 2.43271, 2.43967, 2.4476, 2.45707, 2.46951}, 
{0.982696, 1.00206, 1.02247, 1.04405, 1.06692, 1.09123, 1.11716, 1.14494, 1.17482, 1.20713, 1.24229, 1.28082, 1.32342, 1.37101, 1.42493, 1.48711, 1.56064, 1.5943, 1.63108, 1.67167, 1.71702, 1.7685, 1.82829, 1.90008, 1.9425, 1.99122, 2.04895, 2.08267, 2.12104, 2.13809, 2.15636, 2.1761, 2.19764, 2.22147, 2.24834, 2.27951, 2.29738, 2.31743, 2.34054, 2.35374, 2.36848, 2.3854, 2.40575, 2.41541, 2.42645, 2.43271, 2.43966, 2.4476, 2.45707, 2.46951}, 
{0.955004, 0.975376, 0.996827, 1.01947, 1.04344, 1.06888, 1.09599, 1.12497, 1.1561, 1.18971, 1.22621, 1.26614, 1.31018, 1.35928, 1.41476, 1.47858, 1.5538, 1.58816, 1.62564, 1.66694, 1.713, 1.76521, 1.82572, 1.89823, 1.941, 1.99007, 2.04813, 2.08202, 2.12055, 2.13766, 2.15599, 2.17578, 2.19738, 2.22127, 2.24819, 2.2794, 2.2973, 2.31737, 2.34051, 2.35371, 2.36846, 2.38539, 2.40574, 2.4154, 2.42645, 2.43271, 2.43966, 2.4476, 2.45707, 2.46951}, 
{0.921082, 0.94264, 0.965314, 0.989218, 1.01449, 1.04128, 1.06978, 1.1002, 1.13284, 1.168, 1.20613, 1.24774, 1.29355, 1.34449, 1.4019, 1.46773, 1.54506, 1.58029, 1.61866, 1.66085, 1.70782, 1.76095, 1.82238, 1.89581, 1.93904, 1.98856, 2.04706, 2.08116, 2.1199, 2.13709, 2.1555, 2.17537, 2.19704, 2.221, 2.24799, 2.27927, 2.29719, 2.31729, 2.34045, 2.35367, 2.36843, 2.38537, 2.40573, 2.41539, 2.42644, 2.4327, 2.43966, 2.4476, 2.45707, 2.46951}, 
{0.883336, 0.90611, 0.930046, 0.955261, 0.981893, 1.0101, 1.04008, 1.07205, 1.10629, 1.14315, 1.18304, 1.22651, 1.27427, 1.32727, 1.38685, 1.45497, 1.53471, 1.57094, 1.61033, 1.65357, 1.70161, 1.75581, 1.81834, 1.89287, 1.93665, 1.98671, 2.04575, 2.08011, 2.1191, 2.13638, 2.15489, 2.17486, 2.19662, 2.22067, 2.24774, 2.2791, 2.29706, 2.31719, 2.34039, 2.35362, 2.36839, 2.38534, 2.40571, 2.41538, 2.42643, 2.4327, 2.43966, 2.4476, 2.45707, 2.46951}, 
{0.844248, 0.868117, 0.893205, 0.919633, 0.947544, 0.977101, 1.0085, 1.04198, 1.07781, 1.11634, 1.15802, 1.20338, 1.25315, 1.30829, 1.37015, 1.44071, 1.52307, 1.56039, 1.6009, 1.64529, 1.69451, 1.74992, 1.81368, 1.88946, 1.93386, 1.98455, 2.0442, 2.07886, 2.11815, 2.13555, 2.15417, 2.17425, 2.19612, 2.22028, 2.24745, 2.2789, 2.29691, 2.31708, 2.34031, 2.35356, 2.36835, 2.38531, 2.4057, 2.41537, 2.42643, 2.43269, 2.43965, 2.44759, 2.45707, 2.46951}, 
{0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9}};

  const double imomegaqnm55[8][50] = {{0.0948705, 0.094803, 0.0946932, 0.0945348, 0.0943198, 0.0940389, 0.0936804, 0.0932301, 0.09267, 0.0919774, 0.0911228, 0.0900671, 0.0887575, 0.0871199, 0.0850466, 0.0823729, 0.0788289, 0.0770706, 0.0750526, 0.0727103, 0.0699513, 0.0666388, 0.0625545, 0.0573174, 0.054054, 0.0501531, 0.0453206, 0.0423931, 0.038969, 0.0374162, 0.0357305, 0.0338844, 0.0318402, 0.0295431, 0.0269089, 0.0237947, 0.0219803, 0.0199207, 0.017514, 0.0161252, 0.0145605, 0.0127472, 0.0105435, 0.00948819, 0.00827437, 0.00758316, 0.00681219, 0.00592884, 0.00486949, 0.0173495}, 
{0.285817, 0.28556, 0.285175, 0.284643, 0.28394, 0.283038, 0.281902, 0.280489, 0.278746, 0.276603, 0.273974, 0.27074, 0.266744, 0.261763, 0.255474, 0.247385, 0.236685, 0.231383, 0.225304, 0.218251, 0.20995, 0.199989, 0.187715, 0.171984, 0.162185, 0.150475, 0.135971, 0.127186, 0.116912, 0.112252, 0.107194, 0.101656, 0.0955223, 0.0886307, 0.0807275, 0.0713845, 0.0659411, 0.0597624, 0.0525421, 0.0483758, 0.0436815, 0.0382415, 0.0316306, 0.0284646, 0.0248231, 0.0227495, 0.0204366, 0.0177865, 0.0146085, 0.0104097}, 
{0.480328, 0.479717, 0.478889, 0.477812, 0.476446, 0.474745, 0.472651, 0.470091, 0.466978, 0.463197, 0.4586, 0.452993, 0.446112, 0.437588, 0.426883, 0.413177, 0.395124, 0.386201, 0.375984, 0.364147, 0.350231, 0.333554, 0.313026, 0.286746, 0.270386, 0.250845, 0.226651, 0.212, 0.194868, 0.1871, 0.178668, 0.169434, 0.15921, 0.147722, 0.134549, 0.118976, 0.109903, 0.0996047, 0.0875705, 0.0806266, 0.0728027, 0.063736, 0.0527177, 0.047441, 0.0413719, 0.0379158, 0.034061, 0.0296442, 0.0243475, 0.0173495}, 
{0.680557, 0.679315, 0.677766, 0.675862, 0.673549, 0.67076, 0.667416, 0.663415, 0.658633, 0.65291, 0.64604, 0.637751, 0.627674, 0.615292, 0.599856, 0.580217, 0.554498, 0.541831, 0.527355, 0.510615, 0.490971, 0.467468, 0.438585, 0.401662, 0.378702, 0.351292, 0.317377, 0.296847, 0.272846, 0.261966, 0.250155, 0.237224, 0.222906, 0.20682, 0.188374, 0.16657, 0.153867, 0.139448, 0.1226, 0.112878, 0.101924, 0.0892306, 0.0738049, 0.0664175, 0.0579207, 0.0530822, 0.0476854, 0.0415019, 0.0340865, 0.0242893}, 
{0.888198, 0.88596, 0.88332, 0.880219, 0.876586, 0.872336, 0.867364, 0.86154, 0.854705, 0.846653, 0.83712, 0.825754, 0.812083, 0.795442, 0.774871, 0.748895, 0.715106, 0.698537, 0.679644, 0.657847, 0.632324, 0.601851, 0.564477, 0.516789, 0.487173, 0.451846, 0.408167, 0.381741, 0.350856, 0.336857, 0.321664, 0.30503, 0.286614, 0.265926, 0.242205, 0.214167, 0.197833, 0.179293, 0.15763, 0.14513, 0.131046, 0.114725, 0.0948923, 0.085394, 0.0744695, 0.0682486, 0.0613098, 0.0533596, 0.0438255, 0.0312291}, 
{1.10418, 1.10054, 1.09641, 1.0917, 1.08634, 1.08021, 1.07319, 1.06513, 1.05582, 1.04501, 1.03238, 1.0175, 0.999787, 0.978439, 0.952276, 0.919504, 0.877184, 0.856528, 0.833035, 0.806, 0.774421, 0.736806, 0.690779, 0.632177, 0.595837, 0.552533, 0.499038, 0.466693, 0.428905, 0.411781, 0.393199, 0.372856, 0.350338, 0.325043, 0.296043, 0.261768, 0.241802, 0.219141, 0.192661, 0.177383, 0.160169, 0.140221, 0.11598, 0.104371, 0.0910184, 0.083415, 0.0749342, 0.0652173, 0.0535645, 0.0451087}, 
{1.3285, 1.32309, 1.31709, 1.3104, 1.30292, 1.29452, 1.28507, 1.27436, 1.26216, 1.24818, 1.23202, 1.2132, 1.19101, 1.1645, 1.13228, 1.09223, 1.04089, 1.01595, 0.987663, 0.955194, 0.917365, 0.872418, 0.817554, 0.747869, 0.704729, 0.653377, 0.590005, 0.551715, 0.507, 0.486743, 0.464764, 0.440707, 0.41408, 0.384173, 0.34989, 0.309375, 0.285775, 0.25899, 0.227694, 0.209637, 0.189292, 0.165716, 0.137067, 0.123347, 0.107567, 0.0985815, 0.0885587, 0.077075, 0.0633035, 0.0520485}, 
{0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28}};

  /* To actually calculate the ringdown, we will use the following pointers to point to the correct mode */
  const double (*reomegaqnm)[50] = NULL;
  const double (*imomegaqnm)[50] = NULL;

  REAL8 totalMass, finalMass, finalSpin;

  /* Stuff for interpolating the data */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;

  UINT4 i;

  /* Check we do not require more overtones than we have */
  if ( nmodes > 8 )
  {
    XLALPrintError( "Requesting more overtones than we have data to generate!\n");
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* Choose the appropriate data bases on the user requested l an m */
  switch ( l )
  {
    case 2:
      if ( m == 2 )
      {
        reomegaqnm = reomegaqnm22;
        imomegaqnm = imomegaqnm22;
      }
      else if ( m == 1 )
      {
        reomegaqnm = reomegaqnm21;
        imomegaqnm = imomegaqnm21;
      }
      else
      {
        XLALPrintError( "Unsupported combination of l, m (%d, %d)\n", l, m );
        XLAL_ERROR( XLAL_EINVAL );
      }
      break;
    case 3:
      if ( l == 3 )
      {
        reomegaqnm = reomegaqnm33;
        imomegaqnm = imomegaqnm33;
      }
      else
      {
        XLALPrintError( "Unsupported combination of l, m (%d, %d)\n", l, m );
        XLAL_ERROR( XLAL_EINVAL ); 
      }
      break;
    case 4:
      if ( l == 4 )
      {
        reomegaqnm = reomegaqnm44;
        imomegaqnm = imomegaqnm44;
      }
      else
      {
        XLALPrintError( "Unsupported combination of l, m (%d, %d)\n", l, m );
        XLAL_ERROR( XLAL_EINVAL );
      }
      break;
    case 5:
      if ( l == 5 )
      {
        reomegaqnm = reomegaqnm55;
        imomegaqnm = imomegaqnm55;
      }
      else
      {
        XLALPrintError( "Unsupported combination of l, m (%d, %d)\n", l, m );
        XLAL_ERROR( XLAL_EINVAL );
      }
      break;
    default:
      XLALPrintError( "Unsupported combination of l, m (%d, %d)\n", l, m );
      XLAL_ERROR( XLAL_EINVAL );
      break;
  }

  spline = gsl_spline_alloc( gsl_interp_cspline, 50 );
  acc    = gsl_interp_accel_alloc();

  totalMass = params->totalMass;

 /* Call XLALFinalMassSpin() to get mass and spin of the final black hole */
  if ( XLALFinalMassSpin(&finalMass, &finalSpin, params) == XLAL_FAILURE )
  {
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Now get the QNM frequencies from interpolating the above data */
  for ( i = 0; i < nmodes; i++ )
  {
    gsl_spline_init( spline, afinallist, reomegaqnm[i], 50 );
    gsl_interp_accel_reset( acc );
    
    modefreqs->data[i].realf_FIXME = gsl_spline_eval( spline, finalSpin, acc );

    gsl_spline_init( spline, afinallist, imomegaqnm[i], 50 );
    gsl_interp_accel_reset( acc );

    modefreqs->data[i].imagf_FIXME = gsl_spline_eval( spline, finalSpin, acc );

    /* Scale by the appropriate mass factors */
    modefreqs->data[i].realf_FIXME *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
    modefreqs->data[i].imagf_FIXME *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
  }

  /* Free memory and exit */
  gsl_spline_free( spline );
  gsl_interp_accel_free( acc );

  return XLAL_SUCCESS;
}

static INT4 XLALFinalMassSpin(
	REAL8		 *finalMass,
	REAL8		 *finalSpin,
	InspiralTemplate *params
	)

{
  static const REAL8 root9ovr8minus1 = -0.057190958417936644;
  static const REAL8 root12          = 3.4641016151377544;

  REAL8 eta, eta2, eta3;

  /* get a local copy of the intrinsic parameters */
  eta = params->eta;
  eta2 = eta * eta;

  if ( params->approximant == EOBNRv2 || params->approximant == EOBNRv2HM )
  {
    eta3 = eta2 * eta;
    /* Final mass and spin given by a fitting in Pan et al, arXiv:1106.1021v1 [gr-qc] */
    *finalMass = 1. + root9ovr8minus1 * eta - 0.4333 * eta2 - 0.4392 * eta3;
    *finalSpin = root12 * eta - 3.871 * eta2 + 4.028 * eta3;
  }
  else
  {
    /* Final mass and spin given by a fitting in PRD76, 104049 */
    *finalMass = 1 - 0.057191 * eta - 0.498 * eta2;
    *finalSpin = 3.464102 * eta - 2.9 * eta2;
  }

  return XLAL_SUCCESS;
}

INT4 XLALInspiralHybridAttachRingdownWave (
      REAL4Vector 	*signal1,
      REAL4Vector 	*signal2,
      INT4              l,
      INT4              m,
      REAL8Vector       *timeVec,
      REAL8Vector	*matchrange,
      InspiralTemplate 	*params)
{
      COMPLEX8Vector *modefreqs;
      UINT4 Nrdwave;
      UINT4 j;
      INT4 errcode;

      UINT4 nmodes;
      REAL4Vector		*rdwave1;
      REAL4Vector		*rdwave2;
      REAL4Vector		*rinspwave;
      REAL4Vector		*dinspwave;
      REAL4Vector		*ddinspwave;
      REAL4VectorSequence	*inspwaves1;
      REAL4VectorSequence	*inspwaves2;
      REAL8 dt;
      REAL8 mTot; /* In geometric units */

   /* Attaching position set by omega_match */
   /* Omega_match is given by Eq.(37) of PRD 76, 104049 (2007) */
   /* -0.01 because the current EOBNR 4PN setting can't reach omega_match */

      dt = 1./params->tSampling;
      // UNUSED!!: REAL8 tmatch = matchrange->data[1];

      // UNUSED!!: REAL8 c1 = 1./(LAL_PI*LAL_MTSUN_SI*params->totalMass);
      mTot  = params->totalMass * LAL_MTSUN_SI;

      /* Create memory for the QNM frequencies */
      nmodes = 8;
      modefreqs = XLALCreateCOMPLEX8Vector( nmodes );
      if ( !modefreqs )
      {
        XLAL_ERROR( XLAL_ENOMEM );
      }

      errcode = XLALGenerateQNMFreqV2( modefreqs, params, l, m, nmodes );
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( XLAL_EFUNC );
      }

      /* Ringdown signal length: 10 times the decay time of the n=0 mode */
      Nrdwave = (INT4) (10 / cimagf(modefreqs->data[0]) / dt);

      /* Check the value of attpos, to prevent memory access problems later */
      if ( matchrange->data[0] * mTot / dt < 5 || matchrange->data[1]*mTot/dt > matchrange->data[2] *mTot/dt - 2 )
      {
        XLALPrintError( "More inpiral points needed for ringdown matching.\n" );
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( XLAL_EFAILED );
      }

      /* Create memory for the ring-down and full waveforms, and derivatives of inspirals */

      rdwave1 = XLALCreateREAL4Vector( Nrdwave );
      rdwave2 = XLALCreateREAL4Vector( Nrdwave );
      rinspwave = XLALCreateREAL4Vector( 6 );
      dinspwave = XLALCreateREAL4Vector( 6 );
      ddinspwave = XLALCreateREAL4Vector( 6 );
      inspwaves1 = XLALCreateREAL4VectorSequence( 3, 6 );
      inspwaves2 = XLALCreateREAL4VectorSequence( 3, 6 );

      /* Check memory was allocated */
      if ( !rdwave1 || !rdwave2 || !rinspwave || !dinspwave 
	   || !ddinspwave || !inspwaves1 || !inspwaves2 )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        if (rdwave1)    XLALDestroyREAL4Vector( rdwave1 );
        if (rdwave2)    XLALDestroyREAL4Vector( rdwave2 );
        if (rinspwave)  XLALDestroyREAL4Vector( rinspwave );
        if (dinspwave)  XLALDestroyREAL4Vector( dinspwave );
        if (ddinspwave) XLALDestroyREAL4Vector( ddinspwave );
        if (inspwaves1) XLALDestroyREAL4VectorSequence( inspwaves1 );
        if (inspwaves2) XLALDestroyREAL4VectorSequence( inspwaves2 );
        XLAL_ERROR( XLAL_ENOMEM );
      }

      memset( rdwave1->data, 0, rdwave1->length * sizeof( REAL4 ) );
      memset( rdwave2->data, 0, rdwave2->length * sizeof( REAL4 ) );

      /* Generate derivatives of the last part of inspiral waves */
      /* Get derivatives of signal1 */
      errcode = XLALGenerateHybridWaveDerivatives( rinspwave, dinspwave, ddinspwave, timeVec, signal1, 
									matchrange, params );
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLALDestroyREAL4Vector( rdwave1 );
        XLALDestroyREAL4Vector( rdwave2 );
        XLALDestroyREAL4Vector( rinspwave );
        XLALDestroyREAL4Vector( dinspwave );
        XLALDestroyREAL4Vector( ddinspwave );
        XLALDestroyREAL4VectorSequence( inspwaves1 );
        XLALDestroyREAL4VectorSequence( inspwaves2 );
        XLAL_ERROR( XLAL_EFUNC );
      }
      for (j = 0; j < 6; j++)
      {
	    inspwaves1->data[j] = rinspwave->data[j];
	    inspwaves1->data[j + 6] = dinspwave->data[j];
	    inspwaves1->data[j + 12] = ddinspwave->data[j];
      }

      /* Get derivatives of signal2 */
      errcode = XLALGenerateHybridWaveDerivatives( rinspwave, dinspwave, ddinspwave, timeVec, signal2, 
									matchrange, params );
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLALDestroyREAL4Vector( rdwave1 );
        XLALDestroyREAL4Vector( rdwave2 );
        XLALDestroyREAL4Vector( rinspwave );
        XLALDestroyREAL4Vector( dinspwave );
        XLALDestroyREAL4Vector( ddinspwave );
        XLALDestroyREAL4VectorSequence( inspwaves1 );
        XLALDestroyREAL4VectorSequence( inspwaves2 );
        XLAL_ERROR( XLAL_EFUNC );
      }
      for (j = 0; j < 6; j++)
      {
	    inspwaves2->data[j] = rinspwave->data[j];
	    inspwaves2->data[j + 6] = dinspwave->data[j];
	    inspwaves2->data[j + 12] = ddinspwave->data[j];
      }


      /* Generate ring-down waveforms */
      errcode = XLALInspiralHybridRingdownWave( rdwave1, rdwave2, params, inspwaves1, inspwaves2,
								      modefreqs, matchrange );
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLALDestroyREAL4Vector( rdwave1 );
        XLALDestroyREAL4Vector( rdwave2 );
        XLALDestroyREAL4Vector( rinspwave );
        XLALDestroyREAL4Vector( dinspwave );
        XLALDestroyREAL4Vector( ddinspwave );
        XLALDestroyREAL4VectorSequence( inspwaves1 );
        XLALDestroyREAL4VectorSequence( inspwaves2 );
        XLAL_ERROR( XLAL_EFUNC );
      }

      /* Generate full waveforms, by stitching inspiral and ring-down waveforms */
      UINT4 attachIdx = matchrange->data[1] * mTot / dt;
      for (j = 1; j < Nrdwave; ++j)
      {
	    signal1->data[j + attachIdx] = rdwave1->data[j];
	    signal2->data[j + attachIdx] = rdwave2->data[j];
      }

      memset( signal1->data+Nrdwave+attachIdx, 0, (signal1->length - Nrdwave - attachIdx)*sizeof(REAL4) );
      memset( signal2->data+Nrdwave+attachIdx, 0, (signal2->length - Nrdwave - attachIdx)*sizeof(REAL4) );

      /* Free memory */
      XLALDestroyCOMPLEX8Vector( modefreqs );
      XLALDestroyREAL4Vector( rdwave1 );
      XLALDestroyREAL4Vector( rdwave2 );
      XLALDestroyREAL4Vector( rinspwave );
      XLALDestroyREAL4Vector( dinspwave );
      XLALDestroyREAL4Vector( ddinspwave );
      XLALDestroyREAL4VectorSequence( inspwaves1 );
      XLALDestroyREAL4VectorSequence( inspwaves2 );

      return errcode;
}

INT4 XLALInspiralAttachRingdownWave (
      REAL4Vector 	*Omega,
      REAL4Vector 	*signal1,
      REAL4Vector 	*signal2,
      InspiralTemplate 	*params)
{
      COMPLEX8Vector *modefreqs;
      UINT4 Nrdwave, Npatch;
      UINT4 attpos = 0;
      UINT4 j;
      INT4 errcode;

      UINT4 nmodes;
      INT4 l, m;
      REAL4Vector		*rdwave1;
      REAL4Vector		*rdwave2;
      REAL4Vector		*inspwave;
      REAL4Vector		*dinspwave;
      REAL4Vector		*ddinspwave;
      REAL4VectorSequence	*inspwaves1;
      REAL4VectorSequence	*inspwaves2;
      REAL8 omegamatch, dt;

   /* Attaching position set by omega_match */
   /* Omega_match is given by Eq.(37) of PRD 76, 104049 (2007) */
   /* -0.01 because the current EOBNR 4PN setting can't reach omega_match */

      dt = 1./params->tSampling;
      omegamatch = -0.01 + 0.133 + 0.183 * params->eta + 0.161 * params->eta * params->eta;

      for (j = 0; j < Omega->length; ++j)
      {

        if(Omega->data[j] > omegamatch)
	    {
	      attpos = j - 1;
	      break;
	    }

      }

      /* Create memory for the QNM frequencies */
      nmodes = 3;
      l = 2;
      m = 2;
      modefreqs = XLALCreateCOMPLEX8Vector( nmodes );
      if ( !modefreqs )
      {
        XLAL_ERROR( XLAL_ENOMEM );
      }

      errcode = XLALGenerateQNMFreq( modefreqs, params, l, m, nmodes );
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( XLAL_EFUNC );
      }

      /* Ringdown signal length: 10 times the decay time of the n=0 mode */
      Nrdwave = (INT4) (10 / cimagf(modefreqs->data[0]) / dt);
      /* Patch length, centered around the matching point "attpos" */
      Npatch = 11;

      /* Check the value of attpos, to prevent memory access problems later */
      if ( attpos < ( Npatch + 1 ) / 2 || attpos + (Npatch - 1) / 2 >= signal1->length )
      {
        XLALPrintError( "Value of attpos inconsistent with given value of Npatch.\n" );
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( XLAL_EFAILED );
      }

      /* Create memory for the ring-down and full waveforms, and derivatives of inspirals */

      rdwave1 = XLALCreateREAL4Vector( Nrdwave );
      rdwave2 = XLALCreateREAL4Vector( Nrdwave );
      inspwave = XLALCreateREAL4Vector( Npatch );
      dinspwave = XLALCreateREAL4Vector( Npatch );
      ddinspwave = XLALCreateREAL4Vector( Npatch );
      inspwaves1 = XLALCreateREAL4VectorSequence( 3, (Npatch + 1) / 2 );
      inspwaves2 = XLALCreateREAL4VectorSequence( 3, (Npatch + 1) / 2 );

      /* Check memory was allocated */
      if ( !rdwave1 || !rdwave2 || !inspwave || !dinspwave || !ddinspwave
           || !inspwaves1 || !inspwaves2 )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        if (rdwave1)    XLALDestroyREAL4Vector( rdwave1 );
        if (rdwave2)    XLALDestroyREAL4Vector( rdwave2 );
        if (inspwave)   XLALDestroyREAL4Vector( inspwave );
        if (dinspwave)  XLALDestroyREAL4Vector( dinspwave );
        if (ddinspwave) XLALDestroyREAL4Vector( ddinspwave );
        if (inspwaves1) XLALDestroyREAL4VectorSequence( inspwaves1 );
        if (inspwaves2) XLALDestroyREAL4VectorSequence( inspwaves2 );
        XLAL_ERROR( XLAL_ENOMEM );
      }

      /* Generate derivatives of the last part of inspiral waves */
      /* Take the last part of signal1 */
      for (j = 0; j < Npatch; j++)
      {
	    inspwave->data[j] = signal1->data[attpos - (Npatch + 1) / 2 + j];
      }
      /* Get derivatives of signal1 */
      errcode = XLALGenerateWaveDerivatives( dinspwave, ddinspwave, inspwave, params );
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLALDestroyREAL4Vector( rdwave1 );
        XLALDestroyREAL4Vector( rdwave2 );
        XLALDestroyREAL4Vector( inspwave );
        XLALDestroyREAL4Vector( dinspwave );
        XLALDestroyREAL4Vector( ddinspwave );
        XLALDestroyREAL4VectorSequence( inspwaves1 );
        XLALDestroyREAL4VectorSequence( inspwaves2 );
        XLAL_ERROR( XLAL_EFUNC );
      }
      for (j = 0; j < (Npatch + 1) / 2; j++)
      {
	    inspwaves1->data[j] = inspwave->data[j];
	    inspwaves1->data[j + (Npatch + 1) / 2] = dinspwave->data[j];
	    inspwaves1->data[j + 2 * (Npatch + 1) / 2] = ddinspwave->data[j];
      }

      /* Take the last part of signal2 */
      for (j = 0; j < Npatch; j++)
      {
	    inspwave->data[j] = signal2->data[attpos - (Npatch + 1) / 2 + j];
      }
      /* Get derivatives of signal2 */
      errcode = XLALGenerateWaveDerivatives( dinspwave, ddinspwave, inspwave, params );
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLALDestroyREAL4Vector( rdwave1 );
        XLALDestroyREAL4Vector( rdwave2 );
        XLALDestroyREAL4Vector( inspwave );
        XLALDestroyREAL4Vector( dinspwave );
        XLALDestroyREAL4Vector( ddinspwave );
        XLALDestroyREAL4VectorSequence( inspwaves1 );
        XLALDestroyREAL4VectorSequence( inspwaves2 );
        XLAL_ERROR( XLAL_EFUNC );
      }
      for (j = 0; j < (Npatch + 1) / 2; j++)
      {
	    inspwaves2->data[j] = inspwave->data[j];
	    inspwaves2->data[j + (Npatch + 1) / 2] = dinspwave->data[j];
	    inspwaves2->data[j + 2 * (Npatch + 1) / 2] = ddinspwave->data[j];
      }


      /* Generate ring-down waveforms */
      errcode = XLALInspiralRingdownWave( rdwave1, rdwave2, params, inspwaves1, inspwaves2,
								      modefreqs, nmodes );
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLALDestroyREAL4Vector( rdwave1 );
        XLALDestroyREAL4Vector( rdwave2 );
        XLALDestroyREAL4Vector( inspwave );
        XLALDestroyREAL4Vector( dinspwave );
        XLALDestroyREAL4Vector( ddinspwave );
        XLALDestroyREAL4VectorSequence( inspwaves1 );
        XLALDestroyREAL4VectorSequence( inspwaves2 );
        XLAL_ERROR( XLAL_EFUNC );
      }
      /* Generate full waveforms, by stitching inspiral and ring-down waveforms */
      for (j = 1; j < Nrdwave; ++j)
      {
	    signal1->data[j + attpos - 1] = rdwave1->data[j];
	    signal2->data[j + attpos - 1] = rdwave2->data[j];
      }

      /* Free memory */
      XLALDestroyCOMPLEX8Vector( modefreqs );
      XLALDestroyREAL4Vector( rdwave1 );
      XLALDestroyREAL4Vector( rdwave2 );
      XLALDestroyREAL4Vector( inspwave );
      XLALDestroyREAL4Vector( dinspwave );
      XLALDestroyREAL4Vector( ddinspwave );
      XLALDestroyREAL4VectorSequence( inspwaves1 );
      XLALDestroyREAL4VectorSequence( inspwaves2 );

      return errcode;
}
