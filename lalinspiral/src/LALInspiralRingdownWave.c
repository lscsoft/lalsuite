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
  REAL8 dt, t1, t2, t3, rt;
  gsl_matrix *coef;
  gsl_vector *hderivs;
  gsl_vector *x;
  gsl_permutation *p;
  REAL8Vector *modeamps;
  int s;
  REAL8 tj;

  dt = 1.0 / params -> tSampling;
  printf("dt: %f\n",dt);
  printf("tSampling: %f\n",params -> tSampling);
  printf("matchrange: %d   %d   %d\n",matchrange->data[0],matchrange->data[1],matchrange->data[2]);
  printf("range diff: %d \n",matchrange->data[0] - matchrange->data[1]);
  t3 = ((int)matchrange->data[0] - (int)matchrange->data[1]) * dt;
  rt = -t3 / 3.;
  t2 = t3 + rt;
  t1 = t2 + rt;
  
  printf("time list: %f   %f   %f \n",t1,t2,t3);
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

  printf("time list: %f   %f   %f \n",t1,t2,t3);
  /* Define the linear system Ax=y */
  /* Matrix A (2*n by 2*n) has block symmetry. Define half of A here as "coef" */
  /* Define y here as "hderivs" */
  for (i = 0; i < nmodes; ++i)
  {
	gsl_matrix_set(coef, 0, i, 1);
	gsl_matrix_set(coef, 1, i, - modefreqs->data[i].im);
	gsl_matrix_set(coef, 2, i, modefreqs->data[i].im * modefreqs->data[i].im
			- modefreqs->data[i].re * modefreqs->data[i].re);
	gsl_matrix_set(coef, 3, i, exp(-modefreqs->data[i].im*t1) * cos(modefreqs->data[i].re*t1));
	gsl_matrix_set(coef, 4, i, exp(-modefreqs->data[i].im*t2) * cos(modefreqs->data[i].re*t2));
	gsl_matrix_set(coef, 5, i, exp(-modefreqs->data[i].im*t3) * cos(modefreqs->data[i].re*t3));
	gsl_matrix_set(coef, 6, i, exp(-modefreqs->data[i].im*t3) * 
				      (-modefreqs->data[i].im * cos(modefreqs->data[i].re*t3)
				       -modefreqs->data[i].re * sin(modefreqs->data[i].re*t3)));
	gsl_matrix_set(coef, 7, i, exp( -modefreqs->data[i].im*t3) * 
				      (( modefreqs->data[i].im*modefreqs->data[i].im
					-modefreqs->data[i].re*modefreqs->data[i].re) 
					*cos(modefreqs->data[i].re*t3)
					+2.0*modefreqs->data[i].re*modefreqs->data[i].im 
					*sin(modefreqs->data[i].re*t3)));
	gsl_matrix_set(coef, 8, i, 0);
	gsl_matrix_set(coef, 9, i, - modefreqs->data[i].re);
	gsl_matrix_set(coef, 10, i,  2 * modefreqs->data[i].re * modefreqs->data[i].im);
	gsl_matrix_set(coef, 11, i, -exp(-modefreqs->data[i].im*t1) * sin(modefreqs->data[i].re*t1));
	gsl_matrix_set(coef, 12, i, -exp(-modefreqs->data[i].im*t2) * sin(modefreqs->data[i].re*t2));
	gsl_matrix_set(coef, 13, i, -exp(-modefreqs->data[i].im*t3) * sin(modefreqs->data[i].re*t3));
	gsl_matrix_set(coef, 14, i, exp(-modefreqs->data[i].im*t3) * 
				      ( modefreqs->data[i].im * sin(modefreqs->data[i].re*t3)
				       -modefreqs->data[i].re * cos(modefreqs->data[i].re*t3)));
	gsl_matrix_set(coef, 15, i, exp( -modefreqs->data[i].im*t3) * 
				      (( modefreqs->data[i].re*modefreqs->data[i].re
					-modefreqs->data[i].im*modefreqs->data[i].im) 
					*sin(modefreqs->data[i].re*t3)
					+2.0*modefreqs->data[i].re*modefreqs->data[i].im 
					*cos(modefreqs->data[i].re*t3)));
  }
  for (i = 0; i < 3; ++i)
  {
	gsl_vector_set(hderivs, i, inspwave1->data[(i + 1) * inspwave1->vectorLength - 1]);
	gsl_vector_set(hderivs, i + nmodes, inspwave2->data[(i + 1) * inspwave2->vectorLength - 1]);
	gsl_vector_set(hderivs, i + 5, inspwave1->data[i * inspwave1->vectorLength]);
	gsl_vector_set(hderivs, i + 5 + nmodes, inspwave2->data[i * inspwave2->vectorLength]);
  }
  gsl_vector_set(hderivs, 3, inspwave1->data[2]);
  gsl_vector_set(hderivs, 3 + nmodes, inspwave2->data[2]);
  gsl_vector_set(hderivs, 4, inspwave1->data[1]);
  gsl_vector_set(hderivs, 4 + nmodes, inspwave2->data[1]);
  
  /* Complete the definition for the rest half of A */
  for (i = 0; i < nmodes; ++i)
  {
	for (k = 0; k < nmodes; ++k)
	{
	  gsl_matrix_set(coef, i, k + nmodes, - gsl_matrix_get(coef, i + nmodes, k));
	  gsl_matrix_set(coef, i + nmodes, k + nmodes, gsl_matrix_get(coef, i, k));
	}
  }

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
 
  /* Call gsl LU decomposition to solve the linear system */
  gslStatus = gsl_linalg_LU_decomp(coef, p, &s);
  if ( gslStatus == GSL_SUCCESS )
  {
    gslStatus = gsl_linalg_LU_solve(coef, p, hderivs, x);
  }
  printf("coefs:  ");
  for (i = 0; i < 16; ++i)
  {
    printf("%8.2e   ",gsl_vector_get(x,i));
  }
  printf("\n");
  printf("Eqn1: %e \n", gsl_vector_get(x,0) + gsl_vector_get(x,1) + gsl_vector_get(x,2) + gsl_vector_get(x,3) + gsl_vector_get(x,4) + gsl_vector_get(x,5) + gsl_vector_get(x,6) + gsl_vector_get(x,7) );
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
  REAL8 dt;
  double *x, *y;
  double ry, dy, dy2;
  UINT4 rt;
  UINT4 *tlist;
  gsl_interp_accel *acc;
  gsl_spline *spline;

  /* Sampling rate from input */
  dt = 1.0 / params -> tSampling;

  tlist = (UINT4 *) LALMalloc(4 * sizeof(UINT4));
  rt = (matchrange->data[1] - matchrange->data[0]) / 3;
  tlist[0] = matchrange->data[0];
  tlist[1] = tlist[0] + rt;
  tlist[2] = tlist[1] + rt;
  tlist[3] = matchrange->data[1];
  printf("tlist in deriv: %d, %d, %d, %d \n",tlist[0], tlist[1], tlist[2], tlist[3]);

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

  printf( "%u %u %u %u \n", tlist[0], tlist[1], tlist[2], tlist[3] );
  printf( "%e %e %e %e \n", y[tlist[0]], y[tlist[1]], y[tlist[2]], y[tlist[3]] );

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
  for (j = 0; j < 4; ++j)
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
  printf("In Deriv: \n");
  printf("rwave: %e, %e, %e, %e \n", rwave->data[0],rwave->data[1],rwave->data[2], rwave->data[3]);
  printf("dwave: %e, %e, %e, %e \n", dwave->data[0],dwave->data[1],dwave->data[2], dwave->data[3]);
  printf("ddwave: %e, %e, %e, %e \n", ddwave->data[0], ddwave->data[1], ddwave->data[2], ddwave->data[3]); 
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
  REAL4 BCWre[8][3] = { {1.5251, -1.1568,  0.1292}, {1.3673, -1.0260,  0.1628}, { 1.3223, -1.0257,  0.1860}, 
			{ 1.25, -1.0257,  0.1860},{ 1.2, -1.0257,  0.1860},{ 1.15, -1.0257,  0.1860}, 
			{ 1.1, -1.0257,  0.1860},{ 1.05, -1.0257,  0.1860} };
  REAL4 BCWim[8][3] = { {0.7000,  1.4187, -0.4990}, {0.1000,  0.5436, -0.4731}, {-0.1000,  0.4206, -0.4256},
			{-0.1000,  0.4206, -0.4156},{-0.1000,  0.4206, -0.4056},{-0.1000,  0.4206, -0.3956}, 
			{-0.1000,  0.4206, -0.3856}, {-0.1000,  0.4206, -0.3756} };

  /* Get a local copy of the intrinstic parameters */
  totalMass = params->totalMass;
  finalMass = 0;
  finalSpin = 0;

  /* Call XLALFinalMassSpin() to get mass and spin of the final black hole */
  errcode = XLALFinalMassSpin(&finalMass, &finalSpin, params);
  if ( errcode != XLAL_SUCCESS )
  {
	  XLAL_ERROR( "XLALGenerateQNMFreq", XLAL_EFUNC );
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


INT4 XLALFinalMassSpin(
	REAL8		 *finalMass,
	REAL8		 *finalSpin,
	InspiralTemplate *params
	)

{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  REAL8 eta, eta2;

  /* get a local copy of the intrinstic parameters */
  eta = params->eta;
  eta2 = eta * eta;

  /* Final mass and spin given by a fitting in PRD76, 104049 */
  *finalMass = 1 - 0.057191 * eta - 0.498 * eta2;
  *finalSpin = 3.464102 * eta - 2.9 * eta2;

  return errcode;
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

      printf("length of signals: %d %d \n", signal1->length, signal2->length);
      for ( j = 0; j < 160; ++j)
      {
        printf("%e      %e \n",signal1->data[j],signal2->data[j]);
      }

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

      errcode = XLALGenerateQNMFreq( modefreqs, params, l, m, nmodes );
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( func, XLAL_EFUNC );
      }

      /*
      fprintf(stderr, "f0=%e, f1=%e f2=%e\n",
		      (modefreqs->data[0].re)/LAL_TWOPI,
		      (modefreqs->data[1].re)/LAL_TWOPI,
		      (modefreqs->data[2].re)/LAL_TWOPI);
      */

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
      rinspwave = XLALCreateREAL4Vector( 4 );
      dinspwave = XLALCreateREAL4Vector( 4 );
      ddinspwave = XLALCreateREAL4Vector( 4 );
      inspwaves1 = XLALCreateREAL4VectorSequence( 3, 4 );
      inspwaves2 = XLALCreateREAL4VectorSequence( 3, 4 );

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
      for (j = 0; j < 4; j++)
      {
	    inspwaves1->data[j] = rinspwave->data[j];
	    inspwaves1->data[j + 4] = dinspwave->data[j];
	    inspwaves1->data[j + 8] = ddinspwave->data[j];
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
      for (j = 0; j < 4; j++)
      {
	    inspwaves2->data[j] = rinspwave->data[j];
	    inspwaves2->data[j + 4] = dinspwave->data[j];
	    inspwaves2->data[j + 8] = ddinspwave->data[j];
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

      /*
      fprintf(stderr, "f0=%e, f1=%e f2=%e\n",
		      (modefreqs->data[0].re)/LAL_TWOPI,
		      (modefreqs->data[1].re)/LAL_TWOPI,
		      (modefreqs->data[2].re)/LAL_TWOPI);
      */

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
