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

\author Yi Pan
\file

\brief Module to compute the ring-down waveform as linear combination
of quasi-normal-modes decaying waveforms, which can be attached to
the inspiral part of the compat binary coalescing waveform.

\heading{Prototypes}


<tt>XLALXLALInspiralRingdownWave()</tt>
<ul>
   <li> <tt>rdwave1,</tt> Output, the real part of the ring-down waveform
   </li><li> <tt>rdwave2,</tt> Output, the imaginary part of the ring-down waveform
   </li><li> <tt>params,</tt> Input, the parameters where ring-down waveforms are computed
   </li><li> <tt>inspwave1,</tt> Input, the real part of the ring-down waveform
   </li><li> <tt>inspwave2,</tt> Input, the real part of the ring-down waveform
   </li><li> <tt>modefreqs,</tt> Input, the frequencies of the quasi-normal-modes
   </li><li> <tt>nmode,</tt> Input, the number of quasi-normal-modes to be combined.</li>
</ul>


<tt>XLALGenerateWaveDerivatives()</tt>
<ul>
   <li> <tt>dwave,</tt> Output, time derivative of the input waveform
   </li><li> <tt>ddwave,</tt> Output, two time derivative of the input waveform
   </li><li> <tt>wave,</tt> Input, waveform to be differentiated in time
   </li><li> <tt>params,</tt> Input, the parameters of the input waveform.</li>
</ul>


<tt>XLALGenerateQNMFreq()</tt>
<ul>
   <li> <tt>ptfwave,</tt> Output, the frequencies of the quasi-normal-modes
   </li><li> <tt>params,</tt> Input, the parameters of the binary system
   </li><li> <tt>l,</tt> Input, the l of the modes
   </li><li> <tt>m,</tt> Input, the m of the modes
   </li><li> <tt>nmodes,</tt> Input, the number of overtones considered.</li>
</ul>


<tt>XLALFinalMassSpin()</tt>
<ul>
   <li> <tt>finalMass,</tt> Output, the mass of the final Kerr black hole
   </li><li> <tt>finalSpin,</tt>  Input, the spin of the final Kerr balck hole
   </li><li> <tt>params,</tt> Input, the parameters of the binary system.</li>
</ul>

\heading{Description}
Generating ring-down waveforms.

\heading{Algorithm}

\heading{Uses}

\code
LALMalloc
LALFree
\endcode

\heading{Notes}

*/

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


/* <lalVerbatim file="XLALInspiralHybridRingdownWaveCP">  */
INT4 XLALInspiralHybridRingdownWave (
	REAL4Vector			*rdwave1,
	REAL4Vector			*rdwave2,
	InspiralTemplate		*params,
	REAL4VectorSequence		*inspwave1,
	REAL4VectorSequence		*inspwave2,
	COMPLEX8Vector			*modefreqs,
	REAL8Vector			*matchrange
	)
/* </lalVerbatim> */
{

  static const char *func = "XLALInspiralHybridRingdownWave";

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
    XLAL_ERROR( func, XLAL_EBADLEN );
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

    XLAL_ERROR( func, XLAL_ENOMEM );
  }

  /* Define the linear system Ax=y */
  /* Matrix A (2*n by 2*n) has block symmetry. Define half of A here as "coef" */
  /* Define y here as "hderivs" */
  for (i = 0; i < nmodes; ++i)
  {
	gsl_matrix_set(coef, 0, i, 1);
	gsl_matrix_set(coef, 1, i, - modefreqs->data[i].im);
	gsl_matrix_set(coef, 2, i, exp(-modefreqs->data[i].im*t1) * cos(modefreqs->data[i].re*t1));
	gsl_matrix_set(coef, 3, i, exp(-modefreqs->data[i].im*t2) * cos(modefreqs->data[i].re*t2));
	gsl_matrix_set(coef, 4, i, exp(-modefreqs->data[i].im*t3) * cos(modefreqs->data[i].re*t3));
	gsl_matrix_set(coef, 5, i, exp(-modefreqs->data[i].im*t4) * cos(modefreqs->data[i].re*t4));
	gsl_matrix_set(coef, 6, i, exp(-modefreqs->data[i].im*t5) * cos(modefreqs->data[i].re*t5));
	gsl_matrix_set(coef, 7, i, exp(-modefreqs->data[i].im*t5) * 
				      (-modefreqs->data[i].im * cos(modefreqs->data[i].re*t5)
				       -modefreqs->data[i].re * sin(modefreqs->data[i].re*t5)));
	gsl_matrix_set(coef, 8, i, 0);
	gsl_matrix_set(coef, 9, i, - modefreqs->data[i].re);
	gsl_matrix_set(coef, 10, i, -exp(-modefreqs->data[i].im*t1) * sin(modefreqs->data[i].re*t1));
	gsl_matrix_set(coef, 11, i, -exp(-modefreqs->data[i].im*t2) * sin(modefreqs->data[i].re*t2));
	gsl_matrix_set(coef, 12, i, -exp(-modefreqs->data[i].im*t3) * sin(modefreqs->data[i].re*t3));
	gsl_matrix_set(coef, 13, i, -exp(-modefreqs->data[i].im*t4) * sin(modefreqs->data[i].re*t4));
	gsl_matrix_set(coef, 14, i, -exp(-modefreqs->data[i].im*t5) * sin(modefreqs->data[i].re*t5));
	gsl_matrix_set(coef, 15, i, exp(-modefreqs->data[i].im*t5) * 
				      ( modefreqs->data[i].im * sin(modefreqs->data[i].re*t5)
				       -modefreqs->data[i].re * cos(modefreqs->data[i].re*t5)));
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
  gslStatus = gsl_linalg_LU_decomp(coef, p, &s);
  if ( gslStatus == GSL_SUCCESS )
  {
    gslStatus = gsl_linalg_LU_solve(coef, p, hderivs, x);
  }
  if ( gslStatus != GSL_SUCCESS )
  {
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    XLAL_ERROR( func, XLAL_EFUNC );
  }

  /* Putting solution to an XLAL vector */
  modeamps = XLALCreateREAL8Vector(2 * nmodes);

  if ( !modeamps )
  {
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    XLAL_ERROR( func, XLAL_ENOMEM );
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
	  rdwave1->data[j] += exp(- tj * modefreqs->data[i].im)
			* ( modeamps->data[i] * cos(tj * modefreqs->data[i].re)
			+   modeamps->data[i + nmodes] * sin(tj * modefreqs->data[i].re) );
	  rdwave2->data[j] += exp(- tj * modefreqs->data[i].im)
			* (- modeamps->data[i] * sin(tj * modefreqs->data[i].re)
			+   modeamps->data[i + nmodes] * cos(tj * modefreqs->data[i].re) );
	}
  }

  XLALDestroyREAL8Vector(modeamps);
  return errcode;
}

/* <lalVerbatim file="XLALInspiralRingdownWaveCP">  */
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

  static const char *func = "XLALInspiralRingdownWave";

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
    XLAL_ERROR( func, XLAL_EBADLEN );
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

    XLAL_ERROR( func, XLAL_ENOMEM );
  }

  /* Define the linear system Ax=y */
  /* Matrix A (2*n by 2*n) has block symmetry. Define half of A here as "coef" */
  /* Define y here as "hderivs" */
  for (i = 0; i < nmodes; ++i)
  {
	gsl_matrix_set(coef, 0, i, 1);
	gsl_matrix_set(coef, 1, i, - modefreqs->data[i].im);
	gsl_matrix_set(coef, 2, i, modefreqs->data[i].im * modefreqs->data[i].im
			- modefreqs->data[i].re * modefreqs->data[i].re);
	gsl_matrix_set(coef, 3, i, 0);
	gsl_matrix_set(coef, 4, i, - modefreqs->data[i].re);
	gsl_matrix_set(coef, 5, i,  2 * modefreqs->data[i].re * modefreqs->data[i].im);

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
    XLAL_ERROR( func, XLAL_EFUNC );
  }

  /* Putting solution to an XLAL vector */
  modeamps = XLALCreateREAL8Vector(2 * nmodes);

  if ( !modeamps )
  {
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    XLAL_ERROR( func, XLAL_ENOMEM );
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
	  rdwave1->data[j] += exp(- tj * modefreqs->data[i].im)
			* ( modeamps->data[i] * cos(tj * modefreqs->data[i].re)
			+   modeamps->data[i + nmodes] * sin(tj * modefreqs->data[i].re) );
	  rdwave2->data[j] += exp(- tj * modefreqs->data[i].im)
			* (- modeamps->data[i] * sin(tj * modefreqs->data[i].re)
			+   modeamps->data[i + nmodes] * cos(tj * modefreqs->data[i].re) );
	}
  }

  XLALDestroyREAL8Vector(modeamps);
  return errcode;
}

/* <lalVerbatim file="XLALGenerateHybridWaveDerivatives">  */
INT4 XLALGenerateHybridWaveDerivatives (
	REAL4Vector				*rwave,
	REAL4Vector				*dwave,
	REAL4Vector				*ddwave,
        REAL8Vector				*time,
	REAL4Vector				*wave,
	REAL8Vector				*matchrange,
	InspiralTemplate			*params
	)
/* </lalVerbatim> */
{
  static const char *func = "XLALGenerateHybridWaveDerivatives";

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
    XLAL_ERROR( func, XLAL_ENOMEM );
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
    XLAL_ERROR( func, XLAL_ENOMEM );
  }

  /* Gall gsl spline interpolation */
  gslStatus = gsl_spline_init(spline, time->data, y, vecLength);
  if ( gslStatus != GSL_SUCCESS )
  { 
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    LALFree( y );
    XLAL_ERROR( func, XLAL_EFUNC );
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
      XLAL_ERROR( func, XLAL_EFUNC );
    }
    rwave->data[j]  = (REAL4)(ry);
    dwave->data[j]  = (REAL4)(dy/m);
    ddwave->data[j] = (REAL4)(dy2/m/m);

  }
  
  /* Free gsl variables */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  LALFree(y);

  return errcode;
}

/* <lalVerbatim file="XLALGenerateWaveDerivatives">  */
INT4 XLALGenerateWaveDerivatives (
	REAL4Vector			*dwave,
	REAL4Vector			*ddwave,
	REAL4Vector			*wave,
	InspiralTemplate		*params
	)

{
  static const char *func = "XLALGenerateWaveDerivatives";

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
    XLAL_ERROR( func, XLAL_ENOMEM );
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
    XLAL_ERROR( func, XLAL_ENOMEM );
  }

  /* Gall gsl spline interpolation */
  XLAL_CALLGSL( gslStatus = gsl_spline_init(spline, x, y, wave->length) );
  if ( gslStatus != GSL_SUCCESS )
  { 
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    LALFree( x );
    LALFree( y );
    XLAL_ERROR( func, XLAL_EFUNC );
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
      XLAL_ERROR( func, XLAL_EFUNC );
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

/* <lalVerbatim file="XLALGenerateQNMFreqCP">  */
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
	  XLAL_ERROR( __func__, XLAL_EFUNC );
  }

  /* QNM frequencies from the fitting given in PRD73, 064030 */
  for (i = 0; i < nmodes; ++i)
  {
	modefreqs->data[i].re = BCWre[i][0] + BCWre[i][1] * pow(1.- finalSpin, BCWre[i][2]);
	modefreqs->data[i].im = modefreqs->data[i].re / 2
			     / (BCWim[i][0] + BCWim[i][1] * pow(1.- finalSpin, BCWim[i][2]));
	modefreqs->data[i].re *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
	modefreqs->data[i].im *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
  }
  return errcode;
}


/**
 * As with the above function, this generates the quasinormal mode freuqencies for a black
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

  static const double afinallist[15] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.998};

  /* 2, 2 mode */

  const double reomegaqnm22[8][15] = {{0.37367168, 0.38701754, 0.40214532, 0.41952668, 0.43984192, 0.46412303, 0.49404478, 0.51196924, 
                                      0.53260024, 0.55681739, 0.58601697, 0.62263533, 0.67161427, 0.74632, 0.93852357},
                                      {0.346711, 0.36190968, 0.37897636, 0.39839033, 0.42084668, 0.44740704, 0.47980667, 0.49907908, 
                                      0.52116077, 0.54696353, 0.5779224, 0.61650844, 0.66765755, 0.74458247, 0.93850185},
                                      {0.30105345, 0.31915347, 0.33928514, 0.36192723, 0.38777925, 0.41792504, 0.45417909, 0.4755449, 
                                       0.49990626, 0.52827695, 0.56223982, 0.60442974, 0.65982668, 0.74123888, 0.9384539},
                                      {0.25150496, 0.27174872, 0.29439701, 0.31988611, 0.34885586, 0.3823089, 0.42197661, 0.44509933, 
                                       0.4713363, 0.50192525, 0.53895598, 0.58585025, 0.64786903, 0.73642637, 0.93837084},
                                      {0.20751458, 0.22843281, 0.25251051, 0.28005651, 0.31154608, 0.34774125, 0.3899022, 0.41387206, 
                                       0.44038499, 0.47039243, 0.50626299, 0.55546523, 0.62983606, 0.73020283, 0.93824837},
                                      {0.1692994, 0.1901402, 0.21594701, 0.24645722, 0.281691, 0.32215673, 0.36889572, 0.39509023, 
                                       0.42345037, 0.45405474, 0.48628343, 0.5167307, 0.53691672, 0.51843771, 0.50405383},
                                      {0.13325234, 0.15484078, 0.18597128, 0.22258463, 0.26357081, 0.3095685, 0.36213312, 0.39162879, 
                                       0.42385357, 0.4594292, 0.49916491, 0.5439254, 0.60222748, 0.72284134, 0.93691856},
                                      {0.092822398, 0.12246772, 0.16611733, 0.20975295, 0.25504803, 0.30435388, 0.36008433, 0.39131906, 
                                       0.42553732, 0.46356655, 0.50665162, 0.55687406, 0.61879884, 0.71673712, 0.93708621}};

  const double imomegaqnm22[8][15] = {{0.088962316, 0.088705699, 0.088311166, 0.087729272, 0.086881962, 0.085638835, 0.083765202, 0.08246182,
                                      0.080792873, 0.078602544, 0.075629552, 0.071393306, 0.064869236, 0.053149008, 0.014533966},
                                      {0.27391488, 0.27245211, 0.2705456, 0.26804857, 0.26473345, 0.26022455, 0.25384686, 0.24957997, 
                                      0.24423832, 0.23735737, 0.22814894, 0.21515199, 0.19525207, 0.1596868, 0.043604956},
                                      {0.47827698, 0.47346308, 0.46787212, 0.46125969, 0.4532595, 0.44328686, 0.4303152, 0.42213116, 
                                      0.41226154, 0.39995695, 0.38389521, 0.36154981, 0.32751838, 0.26699789, 0.072684788},
                                      {0.7051482, 0.6944809, 0.68262839, 0.66924293, 0.65380811, 0.63553649, 0.61315357, 0.59976369, 
                                      0.58430072, 0.56591037, 0.54288812, 0.51141591, 0.46286425, 0.37574114, 0.10177254},
                                      {0.94684489, 0.92872504, 0.90884467, 0.88676436, 0.86173803, 0.83259475, 0.79752655, 0.77695669, 
                                      0.75379965, 0.72760503, 0.69796239, 0.66262211, 0.60329505, 0.48678255, 0.130859},
                                      {1.1956081, 1.168248, 1.1387787, 1.106915, 1.0715254, 1.030673, 0.98126515, 0.9518251, 0.91794953, 
                                      0.87820149, 0.83080266, 0.7767326, 0.74872346, 0.72804394, 0.70819007},
                                      {1.4479106, 1.4071149, 1.3668301, 1.3266175, 1.2834733, 1.233995, 1.1737642, 1.1375257, 1.0954075, 
                                      1.0451815, 0.98308361, 0.90139025, 0.77033698, 0.60152481, 0.39231616},
                                      {1.7038413, 1.6389345, 1.5920622, 1.5489188, 1.5019035, 1.446567, 1.377806, 1.3359779, 1.2871194, 
                                      1.2287195, 1.1567076, 1.0636775, 0.93317832, 0.72183665, 0.3341267}};

  /* 2, 1 mode */

  const double reomegaqnm21[8][15] = {{0.3736716844, 0.380432255, 0.3882478199, 0.3973303975, 0.4079791087, 0.4206323931, 0.4359684716, 0.4449684385, 
                                      0.4551214855, 0.4667274647, 0.480230706, 0.4963372414, 0.5162914042, 0.5426930312, 0.5797954879},
                                      {0.3467109969, 0.3544782503, 0.3635061125, 0.3740143534, 0.38631329, 0.4008560371, 0.4183370293, 0.4285085622, 
                                      0.4398971095, 0.4527929231, 0.4676119321, 0.4849776663, 0.505861689, 0.5315724913, 0.5406171253},
                                      {0.3010534545, 0.3104096213, 0.3213886443, 0.3342175299, 0.3492168074, 0.3668517051, 0.3878251641, 0.399890218, 
                                      0.4132608711, 0.4282058471, 0.4450882227, 0.4643918906, 0.4866503513, 0.5111835709, 0.5076280007},
                                      {0.2515049622, 0.2620584381, 0.2746885407, 0.2896276416, 0.3071971794, 0.3278553229, 0.3522784211, 0.3662071434, 
                                      0.3815095774, 0.3984179086, 0.4172271468, 0.4382861751, 0.4618546661, 0.4873898485, 0.5066965377},
                                      {0.2075145798, 0.2184795252, 0.2320829538, 0.2485429494, 0.2681392837, 0.2912579006, 0.318456088, 0.3338297059, 
                                      0.3505633537, 0.3688308697, 0.3888549057, 0.4109548422, 0.4358571062, 0.4673486496, 0.5015154925},
                                      {0.1692994034, 0.1801918011, 0.1946714946, 0.2128105953, 0.2346488155, 0.2602967174, 0.2900089739, 0.3065184561, 
                                      0.324226136, 0.3432341487, 0.3637271461, 0.3862403279, 0.4128682374, 0.4505195661, 0.4960487709},
                                      {0.1332523443, 0.144021248, 0.1606836045, 0.1821862747, 0.2075564398, 0.2363272864, 0.2684211272, 0.2857440639, 
                                      0.3039496637, 0.323111668, 0.3435224899, 0.3663039717, 0.3947709493, 0.4360646256, 0.4923249706},
                                      {0.09282239791, 0.1052666743, 0.1298169054, 0.1581240633, 0.1880936731, 0.219696674, 0.253207736, 0.2707274255, 
                                      0.2887837658, 0.3075307815, 0.3276276643, 0.3509876159, 0.3809678468, 0.4243771567, 0.4900375545}};

  const double imomegaqnm21[8][15] = {{0.08896231569, 0.08879830091, 0.08848852054, 0.08799516823, 0.08725714474, 0.08617299428, 0.08456418202, 
                                      0.08346645962, 0.08208522616, 0.08030936156, 0.07795498003, 0.07468776326, 0.06980434566, 0.06137210553, 0.04007773558},
                                      {0.2739148753, 0.2730584933, 0.2717113664, 0.2697535596, 0.2669920281, 0.2631084667, 0.2575469251, 0.2538377582, 
                                      0.2492378562, 0.2434038652, 0.2357661024, 0.2252822875, 0.2097307961, 0.1828201364, 0.1207320089},
                                      {0.4782769831, 0.4755919695, 0.4719287462, 0.4670759855, 0.4606910631, 0.452210076, 0.440653184, 0.4331961436, 
                                      0.4241435126, 0.412891901, 0.398433546, 0.3788988822, 0.3501876236, 0.2998076248, 0.2049817544},
                                      {0.7051482024, 0.6993167803, 0.6918507693, 0.6824446898, 0.6705918616, 0.6554506444, 0.6355543685, 0.6230364122, 
                                      0.6080901386, 0.5898058886, 0.56665055, 0.5357285396, 0.4904654652, 0.4095631339, 0.2327245081},
                                      {0.9468448909, 0.93705744, 0.9248731358, 0.9099175682, 0.8915244957, 0.8685449148, 0.8389385625, 0.8205471685, 
                                      0.7987559615, 0.7722765031, 0.738921947, 0.6945204084, 0.6295619034, 0.5156508662, 0.3170232474},
                                      {1.195608054, 1.180977518, 1.163244724, 1.142130482, 1.116862815, 1.085922353, 1.046522009, 1.022136197, 0.9932426816, 
                                      0.9580677259, 0.9136012774, 0.8541870404, 0.7677973187, 0.6228480793, 0.4316186791},
                                      {1.447910632, 1.426190421, 1.401396021, 1.373844971, 1.342365373, 1.304556664, 1.256434162, 1.22644742, 1.190653231, 
                                      1.146675625, 1.090529258, 1.015111168, 0.9069587978, 0.7307917488, 0.5179222621},
                                      {1.703841327, 1.66739584, 1.634024328, 1.602185642, 1.567162999, 1.524615916, 1.469166037, 1.433973808, 1.391423129, 
                                      1.338463823, 1.270135027, 1.178335412, 1.048413135, 0.8398281191, 0.5756154785}};

  /* 3, 3 mode */
  const double reomegaqnm33[8][15] = {{0.5994432884, 0.6207956124, 0.6447869769, 0.6720860938, 0.7036502752, 0.7409210977, 0.7862226687, 0.8130568917, 
                                      0.8436868622, 0.8793227375, 0.9218847991, 0.9747239629, 1.044637093, 1.149982884, 1.416011054},
                                      {0.582643803, 0.6052503876, 0.6305726013, 0.659285235, 0.6923520524, 0.731221479, 0.7782251065, 0.8059518959, 
                                      0.8375040485, 0.8740938201, 0.917644637, 0.9715122733, 1.042500234, 1.148967877, 1.415934546},
                                      {0.5516849008, 0.5765672483, 0.6043024574, 0.6355778444, 0.6713709165, 0.7131452055, 0.7632493457, 0.7926061616, 
                                      0.8258461407, 0.864186844, 0.9095602818, 0.9653371401, 1.03834596, 1.146969719, 1.415794397},
                                      {0.5119619111, 0.5395338943, 0.5701634339, 0.6045563477, 0.6437111989, 0.6891154314, 0.7431451362, 0.7745858568, 
                                      0.8099978544, 0.8506074215, 0.8983618867, 0.9566622235, 1.032398165, 1.144044702, 1.415794397},
                                      {0.4701740058, 0.5000419614, 0.5332687341, 0.5705829142, 0.6130100012, 0.6620691136, 0.7201731215, 0.7538215306, 
                                      0.7915637172, 0.8346360247, 0.8850046365, 0.9461143218, 1.024965274, 1.140262781, 1.415587572},
                                      {0.4313864786, 0.462696355, 0.4977495717, 0.5373030628, 0.5824161957, 0.6346489161, 0.6964620562, 0.7321822672, 
                                      0.7721502155, 0.8176102461, 0.8705426742, 0.9344338412, 1.016433028, 1.135699691, 1.41516014},
                                      {0.3976595242, 0.42960427, 0.4657117304, 0.5067718402, 0.5538833351, 0.6086553826, 0.6736092568, 0.711146066, 
                                      0.7531045776, 0.8007332133, 0.8560126947, 0.9224394758, 1.007287679, 1.1304398, 1.41516014},
                                      {0.3689922759, 0.4010191846, 0.437635084, 0.4796608711, 0.5282321113, 0.5850077967, 0.6525769469, 0.6916760588, 
                                      0.7353790205, 0.7849369255, 0.8423170891, 0.9109727593, 0.9981623474, 1.124602562, 1.415160139}};

  const double imomegaqnm33[8][15] = {{0.09270304794, 0.09243049006, 0.09197261689, 0.09126659479, 0.09021787993, 0.08867625749, 0.08638485474, 0.08482127648, 
                                      0.08285568204, 0.08033274097, 0.07699529971, 0.07237329939, 0.06546290444, 0.05338811775, 0.04359929065},
                                      {0.2812981134, 0.2801911651, 0.2785149377, 0.2760799967, 0.2726018653, 0.2676295675, 0.2603936922, 0.2555181302, 
                                      0.2494354597, 0.2416818903, 0.2314894674, 0.2174534353, 0.1965703496, 0.1602243192, 0.1017354847},
                                      {0.479092751, 0.4762674395, 0.472465889, 0.4673691578, 0.4604997228, 0.4511034113, 0.437898944, 0.4291902821, 
                                      0.418467013, 0.4049624308, 0.3874072795, 0.3634776939, 0.3281945301, 0.2672354045, 0.1598806417},
                                      {0.690337096, 0.6844586019, 0.6771829398, 0.6680562408, 0.6564011379, 0.6411526103, 0.620511334, 0.6072206907, 
                                      0.5911010792, 0.5710893065, 0.5454280092, 0.5108990231, 0.4605887925, 0.3745195582, 0.1598806418},
                                      {0.9156493925, 0.9054557014, 0.893406404, 0.878907582, 0.8610734918, 0.8385157223, 0.8089008967, 0.7902253313, 
                                      0.7678830094, 0.7405185187, 0.7059000151, 0.6599481309, 0.5938806138, 0.4821425913, 0.2180346072},
                                      {1.152151362, 1.136946593, 1.119340112, 1.098594431, 1.073606158, 1.042646297, 1.002819632, 0.9780736342, 
                                      0.948772119, 0.9132651079, 0.8688534887, 0.810631559, 0.7280475867, 0.5901180618, 0.305277036},
                                      {1.395912243, 1.375551676, 1.352146664, 1.324811907, 1.292214086, 1.252259659, 1.201449635, 1.170159734, 
                                      1.133351471, 1.089066145, 1.034126577, 0.9628160519, 0.8629344761, 0.6983837287, 0.305277036},
                                      {1.643844528, 1.61845552, 1.589314895, 1.555388359, 1.515100994, 1.465968346, 1.403837974, 1.365751704, 1.321105272, 
                                      1.267604615, 1.20155758, 1.116399914, 0.9983518377, 0.806776984, 0.3052770359}};


  /* 4, 4 mode */
  const double reomegaqnm44[8][15] = {{0.8091783775, 0.8386602441, 0.8717179971, 0.909241654, 0.9525001933, 1.003396037, 1.06498138, 1.101314963, 1.142653897, 
                                      1.190569266, 1.247547008, 1.317912819, 1.41041612, 1.548623269, 1.892594773},
                                      {0.796631532, 0.8270552357, 0.8611135727, 0.8997025841, 0.9440978347, 0.9962087118, 1.059094229, 1.096111622, 1.138157398, 
                                      1.186802528, 1.244532042, 1.315669283, 1.4089583, 1.54795097, 1.892543953},
                                      {0.7727095326, 0.804917568, 0.8408681858, 0.8814697803, 0.9280117764, 0.9824187938, 1.047765773, 1.086080611, 1.129470308, 
                                      1.179506551, 1.238673858, 1.311292766, 1.406099635, 1.54662288, 1.892503315},
                                      {0.73983673, 0.7744212381, 0.8128988806, 0.8561971755, 0.9056276286, 0.9631401563, 1.03183795, 1.071929456, 1.117168066, 
                                      1.16912842, 1.230296848, 1.304993508, 1.401949633, 1.544671097, 1.892503315},
                                      {0.7015155093, 0.7386509975, 0.7798815979, 0.8261602268, 0.8788286478, 0.9398729033, 1.012438934, 1.05460587, 1.10202274, 
                                      1.156270628, 1.219841781, 1.297061094, 1.396662185, 1.542141688, 1.892320149},
                                      {0.6615724994, 0.7009708715, 0.7447317957, 0.7938388505, 0.8496715121, 0.9142626743, 0.9908166022, 1.035164713, 1.084900932, 
                                      1.14161724, 1.207817215, 1.287837432, 1.390424934, 1.539092775, 1.892238519},
                                      {0.6231088747, 0.664176248, 0.7099329205, 0.7613987368, 0.8199986402, 0.8878236757, 0.9681534025, 1.014622099, 1.066652892, 
                                      1.125853947, 1.194746541, 1.277686809, 1.383447464, 1.535592747, 1.891932278},
                                      {0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72}};

  const double imomegaqnm44[8][15] = {{0.09416396099, 0.09394777714, 0.0935373891, 0.09286495746, 0.09182873056, 0.0902678738, 0.08790941931, 0.08628700793, 
                                      0.08424014071, 0.08160780777, 0.0781254881, 0.07331354964, 0.06615598078, 0.05376329412, 0.04363899457},
                                      {0.2843343494, 0.2835188445, 0.2821130849, 0.2799133375, 0.2766142018, 0.2717332243, 0.2644520637, 0.2594796575, 
                                      0.253232947, 0.2452297324, 0.2346778877, 0.220141078, 0.1985760368, 0.1613239058, 0.1018255988},
                                      {0.4799081751, 0.4779855396, 0.4750599634, 0.470791113, 0.4646692749, 0.4558900375, 0.4430904403, 0.4344645829, 
                                      0.4237121852, 0.410031924, 0.3921078361, 0.3675538828, 0.3313130547, 0.2689849723, 0.1309199583},
                                      {0.683924319, 0.6800677823, 0.6747822157, 0.6675889375, 0.6577718751, 0.6442025355, 0.6249719514, 0.6122289479, 
                                      0.5965033728, 0.5766770684, 0.5509149841, 0.5158911006, 0.4645548837, 0.3768081545, 0.1309199583},
                                      {0.8982389718, 0.8914421903, 0.8827803342, 0.8716332875, 0.8570746544, 0.8376469755, 0.8108899014, 0.7934704356, 
                                      0.7722047073, 0.7456591662, 0.711483886, 0.6654232199, 0.5984574414, 0.4848474408, 0.2182097941},
                                      {1.122976754, 1.112315977, 1.099328208, 1.083248751, 1.062936067, 1.036596799, 1.0012099, 0.9785401456, 0.9511440936, 
                                      0.9172726107, 0.8740632363, 0.8163395551, 0.7331382111, 0.5931467615, 0.247309098},
                                      {1.356686268, 1.3415477, 1.323564658, 1.301825067, 1.274966924, 1.240852217, 1.195885839, 1.167453117, 1.13338524, 
                                      1.091615012, 1.038766645, 0.9687480012, 0.8686752793, 0.7017384806, 0.3346155071},
                                      {0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28}};

  /* 5, 5 mode */
  const double reomegaqnm55[8][15] = {{1.012295312, 1.049851616, 1.091940772, 1.139682353, 1.194669326, 1.259284164, 1.337338848, 1.383315031, 1.435553557, 
                                      1.496002653, 1.567738903, 1.656104376, 1.771879595, 1.944029193, 2.36850041},
                                      {1.002221028, 1.04053517, 1.083428403, 1.132025534, 1.187925439, 1.253517091, 1.332619192, 1.379147454, 1.431957628, 
                                      1.49299789, 1.565343843, 1.654334717, 1.770744093, 1.943518167, 2.368493022},
                                      {0.9826957608, 1.022473707, 1.066917524, 1.1171627, 1.174820361, 1.242293291, 1.323415056, 1.371009549, 1.424925441, 
                                      1.487111413, 1.560641928, 1.65085167, 1.768501947, 1.942504562, 2.36847825},
                                      {0.9550040061, 0.996827188, 1.043436943, 1.095985312, 1.156102765, 1.22621453, 1.310179483, 1.359280485, 1.414763553, 
                                      1.478579454, 1.553802746, 1.645763809, 1.765209202, 1.941004962, 2.368456103},
                                      {0.9210818439, 0.9653140523, 1.01448695, 1.069775208, 1.132836271, 1.206127549, 1.29354586, 1.344489472, 1.401899686, 
                                      1.467732064, 1.545064101, 1.639224429, 1.760945894, 1.939043498, 2.368426595},
                                      {0.8833357607, 0.9300460355, 0.9818925575, 1.04007822, 1.106294845, 1.183042074, 1.274268603, 1.327267527, 1.386845151, 
                                      1.454965586, 1.534713631, 1.63142103, 1.755811623, 1.936651014, 2.368389744},
                                      {0.8442481668, 0.8932048876, 0.9475435534, 1.008500191, 1.077807901, 1.158018504, 1.253147249, 1.30828747, 1.370149231, 
                                      1.440709836, 1.523067123, 1.622562923, 1.749920255, 1.933864051, 2.368345567},
                                      {0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9}};

  const double imomegaqnm55[8][15] = {{0.09487051608, 0.09469320694, 0.0943197987, 0.09368038152, 0.09267001286, 0.09112281116, 0.08875752731, 0.08711991414, 
                                      0.08504662409, 0.08237290235, 0.0788288767, 0.07392697215, 0.06663876869, 0.05405397505, 0.0145604752},
                                      {0.2858173818, 0.2851754898, 0.2839402351, 0.2819018594, 0.2787455963, 0.2739740359, 0.2667440213, 0.2617629027, 
                                      0.2554743347, 0.2473845535, 0.2366847744, 0.2219137328, 0.1999889311, 0.1621853591, 0.04368148906},
                                      {0.480328456, 0.4788887352, 0.4764460478, 0.4726506621, 0.466978052, 0.4585997397, 0.4461119952, 0.4375878166, 
                                      0.426883346, 0.4131771317, 0.3951237351, 0.3702922687, 0.3335540436, 0.2703864339, 0.07280269303},
                                      {0.6805569086, 0.6777655047, 0.6735486008, 0.6674159912, 0.6586329957, 0.6460400166, 0.6276736425, 0.6152919829, 
                                      0.5998556529, 0.5802171896, 0.5544975676, 0.5193021255, 0.4674679048, 0.3787016615, 0.1019242128},
                                      {0.8881975924, 0.8833204554, 0.8765864548, 0.8673635371, 0.8547053659, 0.8371197493, 0.8120826774, 0.7954424475, 
                                      0.7748705464, 0.7488949816, 0.7151062528, 0.6691524835, 0.6018509125, 0.4871726397, 0.131046173},
                                      {1.104182464, 1.09640866, 1.086337843, 1.073193872, 1.055816042, 1.032376375, 0.9997868778, 0.9784387696, 0.9522763791, 
                                      0.9195040232, 0.8771841992, 0.820014025, 0.7368063238, 0.5958372163, 0.1601686969},
                                      {1.328496983, 1.317085574, 1.302916522, 1.285067238, 1.262159349, 1.232024641, 1.191005457, 1.16449526, 1.132276893, 
                                      1.092230885, 1.040893034, 0.9720143596, 0.8724179408, 0.7047288414, 0.1892919055},
                                      {0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28}};

  /* To actually calculate the ringdown, we will use the following pointers to point to the correct mode */
  const double (*reomegaqnm)[15] = NULL;
  const double (*imomegaqnm)[15] = NULL;

  REAL8 totalMass, finalMass, finalSpin;

  /* Stuff for interpolating the data */
  gsl_spline    *spline = NULL;
  gsl_interp_accel *acc = NULL;

  UINT4 i;

  /* Check we do not require more overtones than we have */
  if ( nmodes > 8 )
  {
    XLALPrintError( "Requesting more overtones than we have data to generate!\n");
    XLAL_ERROR( __func__, XLAL_EINVAL );
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
        XLAL_ERROR( __func__, XLAL_EINVAL );
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
        XLAL_ERROR( __func__, XLAL_EINVAL ); 
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
        XLAL_ERROR( __func__, XLAL_EINVAL );
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
        XLAL_ERROR( __func__, XLAL_EINVAL );
      }
      break;
    default:
      XLALPrintError( "Unsupported combination of l, m (%d, %d)\n", l, m );
      XLAL_ERROR( __func__, XLAL_EINVAL );
      break;
  }

  spline = gsl_spline_alloc( gsl_interp_cspline, 15 );
  acc    = gsl_interp_accel_alloc();

  totalMass = params->totalMass;

 /* Call XLALFinalMassSpin() to get mass and spin of the final black hole */
  if ( XLALFinalMassSpin(&finalMass, &finalSpin, params) == XLAL_FAILURE )
  {
    XLAL_ERROR( __func__, XLAL_EFUNC );
  }

  /* Now get the QNM frequencies from interpolating the above data */
  for ( i = 0; i < nmodes; i++ )
  {
    gsl_spline_init( spline, afinallist, reomegaqnm[i], 15 );
    gsl_interp_accel_reset( acc );
    
    modefreqs->data[i].re = gsl_spline_eval( spline, finalSpin, acc );

    gsl_spline_init( spline, afinallist, imomegaqnm[i], 15 );
    gsl_interp_accel_reset( acc );

    modefreqs->data[i].im = gsl_spline_eval( spline, finalSpin, acc );

    /* Scale by the appropriate mass factors */
    printf( "%d, %d frequency = %e + i %e\n", l, m, modefreqs->data[i].re, modefreqs->data[i].im );
    modefreqs->data[i].re *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
    modefreqs->data[i].im *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
  }

  /* Free memory and exit */
  gsl_spline_free( spline );
  gsl_interp_accel_free( acc );

  return XLAL_SUCCESS;
}


INT4 XLALFinalMassSpin(
	REAL8		 *finalMass,
	REAL8		 *finalSpin,
	InspiralTemplate *params
	)

{
  static const REAL8 root9ovr8minus1 = -0.057190958417936644;
  static const REAL8 root12          = 3.4641016151377544;

  REAL8 eta, eta2, eta3;

  /* get a local copy of the intrinstic parameters */
  eta = params->eta;
  eta2 = eta * eta;

  if ( params->approximant == EOBNRv2 )
  {
    eta3 = eta2 * eta;
    /* Final mass and spin given by a fitting in Pan et al, in preparation */
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

/* <lalVerbatim file="XLALInspiralHybridAttachRingdownWaveCP">  */
INT4 XLALInspiralHybridAttachRingdownWave (
      REAL4Vector 	*signal1,
      REAL4Vector 	*signal2,
      INT4              l,
      INT4              m,
      REAL8Vector       *time,
      REAL8Vector	*matchrange,
      InspiralTemplate 	*params)
{

      static const char *func = "XLALInspiralAttachRingdownWave";

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
      REAL8 dt, c1;
      REAL8 tmatch;
      REAL8 mTot; /* In geometric units */

   /* Attaching position set by omega_match */
   /* Omega_match is given by Eq.(37) of PRD 76, 104049 (2007) */
   /* -0.01 because the current EOBNR 4PN setting can't reach omega_match */

      dt = 1./params->tSampling;
      tmatch = matchrange->data[1];

      c1 = 1./(LAL_PI*LAL_MTSUN_SI*params->totalMass);
      mTot  = params->totalMass * LAL_MTSUN_SI;

      /* Create memory for the QNM frequencies */
      nmodes = 8;
      modefreqs = XLALCreateCOMPLEX8Vector( nmodes );
      if ( !modefreqs )
      {
        XLAL_ERROR( func, XLAL_ENOMEM );
      }

      errcode = XLALGenerateQNMFreqV2( modefreqs, params, l, m, nmodes );
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( func, XLAL_EFUNC );
      }

      /* Ringdown signal length: 10 times the decay time of the n=0 mode */
      Nrdwave = (INT4) (10 / modefreqs->data[0].im / dt);

      /* Check the value of attpos, to prevent memory access problems later */
      if ( matchrange->data[0] * mTot / dt < 5 || matchrange->data[1]*mTot/dt > matchrange->data[2] *mTot/dt - 2 )
      {
        XLALPrintError( "More inpiral points needed for ringdown matching.\n" );
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( func, XLAL_EFAILED );
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
        XLAL_ERROR( func, XLAL_ENOMEM );
      }

      /* Generate derivatives of the last part of inspiral waves */
      /* Get derivatives of signal1 */
      errcode = XLALGenerateHybridWaveDerivatives( rinspwave, dinspwave, ddinspwave, time, signal1, 
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
        XLAL_ERROR( func, XLAL_EFUNC );
      }
      for (j = 0; j < 6; j++)
      {
	    inspwaves1->data[j] = rinspwave->data[j];
	    inspwaves1->data[j + 6] = dinspwave->data[j];
	    inspwaves1->data[j + 12] = ddinspwave->data[j];
      }

      /* Get derivatives of signal2 */
      errcode = XLALGenerateHybridWaveDerivatives( rinspwave, dinspwave, ddinspwave, time, signal2, 
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
        XLAL_ERROR( func, XLAL_EFUNC );
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
        XLAL_ERROR( func, XLAL_EFUNC );
      }

      /* Generate full waveforms, by stitching inspiral and ring-down waveforms */
      UINT4 attachIdx = matchrange->data[1] * mTot / dt;
      for (j = 1; j < Nrdwave; ++j)
      {
	    signal1->data[j + attachIdx] = rdwave1->data[j];
	    signal2->data[j + attachIdx] = rdwave2->data[j];
      }

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

/* <lalVerbatim file="XLALInspiralAttachRingdownWaveCP">  */
INT4 XLALInspiralAttachRingdownWave (
      REAL4Vector 	*Omega,
      REAL4Vector 	*signal1,
      REAL4Vector 	*signal2,
      InspiralTemplate 	*params)
{

      static const char *func = "XLALInspiralAttachRingdownWave";

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
        XLAL_ERROR( func, XLAL_ENOMEM );
      }

      errcode = XLALGenerateQNMFreq( modefreqs, params, l, m, nmodes );
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( func, XLAL_EFUNC );
      }

      /* Ringdown signal length: 10 times the decay time of the n=0 mode */
      Nrdwave = (INT4) (10 / modefreqs->data[0].im / dt);
      /* Patch length, centered around the matching point "attpos" */
      Npatch = 11;

      /* Check the value of attpos, to prevent memory access problems later */
      if ( attpos < ( Npatch + 1 ) / 2 || attpos + (Npatch - 1) / 2 >= signal1->length )
      {
        XLALPrintError( "Value of attpos inconsistent with given value of Npatch.\n" );
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( func, XLAL_EFAILED );
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
        XLAL_ERROR( func, XLAL_ENOMEM );
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
        XLAL_ERROR( func, XLAL_EFUNC );
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
        XLAL_ERROR( func, XLAL_EFUNC );
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
        XLAL_ERROR( func, XLAL_EFUNC );
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
