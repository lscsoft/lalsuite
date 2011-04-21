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
	UINT4Vector			*matchrange
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

  dt = 1.0 / params -> tSampling;
  t5 = ((int)matchrange->data[0] - (int)matchrange->data[1]) * dt;
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

  /* print ringdown-matching linear system: coefficient matrix and RHS vector */
#if 0  
  printf("matching matrix:\n");
  for (i = 0; i < 16; ++i)
  {
    for (j = 0; j < 16; ++j)
    {
      printf("%8.2f ",gsl_matrix_get(coef,i,j));
    }
    printf("\n");
  }
  printf("RHS:  ");
  for (i = 0; i < 16; ++i)
  {
    printf("%e   ",gsl_vector_get(hderivs,i));
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
	REAL4Vector				*wave,
	UINT4Vector				*matchrange,
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
  double *x, *y;
  double ry, dy, dy2;
  UINT4 rt;
  UINT4 *tlist;
  gsl_interp_accel *acc;
  gsl_spline *spline;

  /* Sampling rate from input */
  dt = 1.0 / params -> tSampling;

  tlist = (UINT4 *) LALMalloc(6 * sizeof(UINT4));
  rt = (matchrange->data[1] - matchrange->data[0]) / 5.;
  tlist[0] = matchrange->data[0];
  tlist[1] = tlist[0] + rt;
  tlist[2] = tlist[1] + rt;
  tlist[3] = tlist[2] + rt;
  tlist[4] = tlist[3] + rt;
  tlist[5] = matchrange->data[1];

  /* Set the length of the interpolation vectors */
  vecLength = matchrange->data[2] + 1;

  /* Getting interpolation and derivatives of the waveform using gsl spline routine */
  /* Initiate arrays and supporting variables for gsl */
  x = (double *) LALMalloc(vecLength * sizeof(double));
  y = (double *) LALMalloc(vecLength * sizeof(double));

  if ( !x || !y )
  {
    if ( x ) LALFree (x);
    if ( y ) LALFree (y);
    XLAL_ERROR( func, XLAL_ENOMEM );
  }

  for (j = 0; j < vecLength; ++j)
  {
	x[j] = j;
	y[j] = wave->data[j];
  }


  XLAL_CALLGSL( acc = (gsl_interp_accel*) gsl_interp_accel_alloc() );
  XLAL_CALLGSL( spline = (gsl_spline*) gsl_spline_alloc(gsl_interp_cspline, vecLength) );
  if ( !acc || !spline )
  {
    if ( acc )    gsl_interp_accel_free(acc);
    if ( spline ) gsl_spline_free(spline);
    LALFree( x );
    LALFree( y );
    XLAL_ERROR( func, XLAL_ENOMEM );
  }

  /* Gall gsl spline interpolation */
  XLAL_CALLGSL( gslStatus = gsl_spline_init(spline, x, y, vecLength) );
  if ( gslStatus != GSL_SUCCESS )
  { 
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    LALFree( x );
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
      LALFree( x );		
      LALFree( y );
      XLAL_ERROR( func, XLAL_EFUNC );
    }
    rwave->data[j]  = (REAL4)(ry);
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
        UINT4                   UNUSED l,
        UINT4                   UNUSED m,
        UINT4                   nmodes
        )
{

  /* Data for interpolating the quasinormal mode frequencies is taken from */
  /* The webpage of Emanuele Berti, http://www.phy.olemiss.edu/~berti/qnms.html */

  static const double afinallist[15] = {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.998};

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
    gsl_spline_init( spline, afinallist, reomegaqnm22[i], 15 );
    gsl_interp_accel_reset( acc );
    
    modefreqs->data[i].re = gsl_spline_eval( spline, finalSpin, acc );

    gsl_spline_init( spline, afinallist, imomegaqnm22[i], 15 );
    gsl_interp_accel_reset( acc );

    modefreqs->data[i].im = gsl_spline_eval( spline, finalSpin, acc );

    /* Scale by the appropriate mass factors */
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
      UINT4Vector	*matchrange,
      InspiralTemplate 	*params)
{

      static const char *func = "XLALInspiralAttachRingdownWave";

      COMPLEX8Vector *modefreqs;
      UINT4 Nrdwave;
      UINT4 j;
      INT4 errcode;

      UINT4 nmodes;
      INT4 l, m;
      REAL4Vector		*rdwave1;
      REAL4Vector		*rdwave2;
      REAL4Vector		*rinspwave;
      REAL4Vector		*dinspwave;
      REAL4Vector		*ddinspwave;
      REAL4VectorSequence	*inspwaves1;
      REAL4VectorSequence	*inspwaves2;
      REAL8 dt, c1;
      UINT4 tmatch;

   /* Attaching position set by omega_match */
   /* Omega_match is given by Eq.(37) of PRD 76, 104049 (2007) */
   /* -0.01 because the current EOBNR 4PN setting can't reach omega_match */

      dt = 1./params->tSampling;
      tmatch = matchrange->data[1];

      c1 = 1./(LAL_PI*LAL_MTSUN_SI*params->totalMass);

      /* Create memory for the QNM frequencies */
      nmodes = 8;
      l = 2;
      m = 2;
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
      if ( matchrange->data[0] < 5 || matchrange->data[1] > matchrange->data[2] - 2 )
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
      errcode = XLALGenerateHybridWaveDerivatives( rinspwave, dinspwave, ddinspwave, signal1, 
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
      errcode = XLALGenerateHybridWaveDerivatives( rinspwave, dinspwave, ddinspwave, signal2, 
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
      for (j = 1; j < Nrdwave; ++j)
      {
	    signal1->data[j + matchrange->data[1] - 1] = rdwave1->data[j];
	    signal2->data[j + matchrange->data[1] - 1] = rdwave2->data[j];
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
