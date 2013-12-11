/*
*  Copyright (C) 2010 R. Sturani, adapted from LALInspiralRingdownWave.c
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
 * \file
 *
 * \brief Module to compute the ring-down waveform as linear combination
 * of quasi-normal-modes decaying waveforms, which can be attached to
 * the phenomenological spin Taylor waveform.
 *
 * ### Prototypes ###
 *
 * <tt>XLALXLALPSpinInspiralRingdownWave()</tt>
 * <ul>
 * <li> <tt>rdwave,</tt> Output, the ring-down waveform
 * </li><li> <tt>params,</tt> Input, the parameters where ring-down waveforms are computed
 * </li><li> <tt>inspwave,</tt> Input, the inspiral waveform with given multiple
 * </li><li> <tt>modefreqs,</tt> Input, the frequencies of the quasi-normal-modes
 * </li><li> <tt>nmodes,</tt> Input, the number of quasi-normal-modes to be combined.
 * </li></ul>
 *
 * <tt>XLALGenerateWaveDerivative()</tt>
 * <ul>
 * <li> <tt>dwave,</tt> Output, time derivative of the input waveform
 * </li><li> <tt>wave,</tt> Input, waveform to be differentiated in time
 * </li><li> <tt>params,</tt> Input, the parameters of the input waveform.
 * </li></ul>
 *
 * <tt>XLALPSpinGenerateQNMFreq()</tt>
 * <ul>
 * <li> <tt>ptfwave,</tt> Output, the frequencies of the quasi-normal-modes
 * </li><li> <tt>params,</tt> Input, the parameters of the binary system
 * </li><li> <tt>l,</tt> Input, the l of the modes
 * </li><li> <tt>m,</tt> Input, the m of the modes
 * </li><li> <tt>nmodes,</tt> Input, the number of overtones considered.
 * </li></ul>
 *
 * <tt>XLALPSpinFinalMassSpin()</tt>
 * <ul>
 * <li> <tt>finalMass,</tt> Output, the mass of the final Kerr black hole
 * </li><li> <tt>finalSpin,</tt>  Input, the spin of the final Kerr balck hole
 * </li><li> <tt>params,</tt> Input, the parameters of the binary system.
 * </li><li> <tt>energy,</tt> Input, the binding energy at the time final time.
 * </li></ul>
 *
 * <tt>XLALPSpinInspiralAttachRingdownWave()</tt>
 * <ul>
 * <li> <tt>sigl,</tt> Output, the waveform filled with ring-down phase
 * </li> <li> <tt>params,</tt> Input, inspiral parameters
 * </li> <li> <tt>attpos,</tt> Input, position of the start of the ring-down
 * </li> <li> <tt>nmodes,</tt> Input, number of ring-down modes
 * </li> <li> <tt>l,</tt> Input, spherical harmonic l-number of the ring-down mode
 * </li> <li> <tt>m,</tt> Input, spherical harmonic m-number of the ring-down mode
 * </li> <li> <tt>finalMass,</tt> Input, estimated final mass of the black hole
 * </li> <li> <tt>finalSpin,</tt> Input, estimated final spin of the black hole</li>
 * </ul>
 *
 * ### Description ###
 *
 * This module generate ring-down waveforms.
 *
 * ### Algorithm ###
 *
 *
 * ### Uses ###
 *
 * \code
 * LALMalloc
 * LALFree
 * \endcode
 *
 * ### Notes ###
 *
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
#include <lal/XLALGSL.h>

#define XLAL_BEGINGSL \
        { \
          gsl_error_handler_t *saveGSLErrorHandler_; \
          XLALGSL_PTHREAD_MUTEX_LOCK; \
          saveGSLErrorHandler_ = gsl_set_error_handler_off();

#define XLAL_ENDGSL \
          gsl_set_error_handler( saveGSLErrorHandler_ ); \
          XLALGSL_PTHREAD_MUTEX_UNLOCK; \
        }

INT4 XLALPSpinInspiralRingdownWave (
	REAL8Vector		*rdwave,
	InspiralTemplate	*params,
	REAL8Vector	        *matchinspwave,
	COMPLEX8Vector		*modefreqs,
	UINT4			nmodes
	)

{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;

  /* Needed to check GSL return codes */
  INT4 gslStatus;

  UINT4 i, j;

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

  if ( modefreqs->length != nmodes )
  {
    XLAL_ERROR( XLAL_EBADLEN );
  }

  /* Solving the linear system for QNMs amplitude coefficients using gsl routine */
  /* Initialize matrices and supporting variables */

  XLAL_BEGINGSL;
  coef = (gsl_matrix *) gsl_matrix_alloc(2 * nmodes, 2 * nmodes);
  hderivs = (gsl_vector *) gsl_vector_alloc(2 * nmodes);
  x = (gsl_vector *) gsl_vector_alloc(2 * nmodes);
  p = (gsl_permutation *) gsl_permutation_alloc(2 * nmodes);

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
  /* Matrix A (2*nmodes by 2*nmodes) has block symmetry. Define half of A here as "coef" */
  /* Define y here as "hderivs" */

  j=0;
  while (j<nmodes) {
    switch (j) {
      case 0:
        for (i = 0; i < nmodes; i++) {
          gsl_matrix_set(coef, 2*j, i, 1.);
          gsl_matrix_set(coef, 2*j, i+nmodes, 0.);
          gsl_matrix_set(coef, 2*j+1, i, -cimagf(modefreqs->data[i]));
          gsl_matrix_set(coef, 2*j+1, i+nmodes, crealf(modefreqs->data[i]));
        }
        break;
      case 1:
        for (i = 0; i < nmodes; i++) {
          gsl_matrix_set(coef, 2*j, i, cimagf(modefreqs->data[i])*cimagf(modefreqs->data[i])-crealf(modefreqs->data[i])*crealf(modefreqs->data[i]));
          gsl_matrix_set(coef, 2*j, i+nmodes, -2.*cimagf(modefreqs->data[i])*crealf(modefreqs->data[i]));
          gsl_matrix_set(coef, 2*j+1, i, -cimagf(modefreqs->data[i])*cimagf(modefreqs->data[i])*cimagf(modefreqs->data[i])+3.*cimagf(modefreqs->data[i])*crealf(modefreqs->data[i])*crealf(modefreqs->data[i]));
          gsl_matrix_set(coef, 2*j+1, i+nmodes, -crealf(modefreqs->data[i])*crealf(modefreqs->data[i])*crealf(modefreqs->data[i])+3.*crealf(modefreqs->data[i])*cimagf(modefreqs->data[i])*cimagf(modefreqs->data[i]));
        }
        break;
      case 2:
        for (i = 0; i < nmodes; i++) {
          gsl_matrix_set(coef, 2*j, i, pow(cimagf(modefreqs->data[i]),4.)+pow(crealf(modefreqs->data[i]),4.)-6.*pow(crealf(modefreqs->data[i])*cimagf(modefreqs->data[i]),2.));
          gsl_matrix_set(coef, 2*j, i+nmodes, -4.*pow(cimagf(modefreqs->data[i]),3.)*crealf(modefreqs->data[i])+4.*pow(crealf(modefreqs->data[i]),3.)*cimagf(modefreqs->data[i]));
          gsl_matrix_set(coef, 2*j+1, i, -pow(cimagf(modefreqs->data[i]),5.)+10.*pow(cimagf(modefreqs->data[i]),3.)*pow(crealf(modefreqs->data[i]),2.)-5.*cimagf(modefreqs->data[i])*pow(crealf(modefreqs->data[i]),4.));
          gsl_matrix_set(coef, 2*j+1, i+nmodes, 5.*pow(cimagf(modefreqs->data[i]),4.)*crealf(modefreqs->data[i])-10.*pow(cimagf(modefreqs->data[i]),2.)*pow(crealf(modefreqs->data[i]),3.)+pow(crealf(modefreqs->data[i]),5.));
        }
        break;
      default:
        XLALPrintError("*** LALPSpinInspiralRingDown ERROR ***: nmode must be <=2, %d selected\n",nmodes);
        gsl_matrix_free(coef);
        gsl_vector_free(hderivs);
        gsl_vector_free(x);
        gsl_permutation_free(p);
        XLAL_ERROR( XLAL_EDOM );
    }
    gsl_vector_set(hderivs, 2*j, matchinspwave->data[2*j]);
    gsl_vector_set(hderivs, 2*j+1, matchinspwave->data[2*j+1]);
    j++;
  }

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
    XLAL_ERROR( XLAL_EFUNC );
  }

  /* Putting solution to an XLAL vector */
  modeamps = XLALCreateREAL8Vector(2*nmodes);

  if ( !modeamps )
  {
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    XLAL_ERROR( XLAL_ENOMEM );
  }

  for (i = 0; i < 2*nmodes; i++) {
    modeamps->data[i] = gsl_vector_get(x, i);
  }

  /* Free all gsl linear algebra objects */
  gsl_matrix_free(coef);
  gsl_vector_free(hderivs);
  gsl_vector_free(x);
  gsl_permutation_free(p);
  XLAL_ENDGSL;

  /* Build ring-down waveforms */
  UINT4 Nrdwave=rdwave->length;
  for (j = 0; j < Nrdwave; j++) {
    tj = j * dt;
    rdwave->data[j] = 0.;
    for (i = 0; i < nmodes; i++) {
      rdwave->data[j] += exp(- tj * cimagf(modefreqs->data[i]))
	* ( modeamps->data[i] * cos(tj * crealf(modefreqs->data[i]))
	    +   modeamps->data[i + nmodes] * sin(tj * crealf(modefreqs->data[i])) );
    }
  }

  XLALDestroyREAL8Vector(modeamps);
  return errcode;
} /*End of XLALPSpinInspiralRingdownWave */


INT4 XLALGenerateWaveDerivative (
	REAL8Vector		*dwave,
	REAL8Vector	        *wave,
	REAL8                    dt
    )
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;

  /* For checking GSL return codes */
  INT4 gslStatus;

  UINT4 j;
  double *x, *y;
  double dy;
  gsl_interp_accel *acc;
  gsl_spline *spline;

  if (wave->length!=dwave->length)
    XLAL_ERROR( XLAL_EFUNC );

  /* Getting interpolation and derivatives of the waveform using gsl spline routine */
  /* Initialize arrays and supporting variables for gsl */

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
    if (gslStatus != GSL_SUCCESS )
    {
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
      LALFree( x );
      LALFree( y );
      XLAL_ERROR( XLAL_EFUNC );
    }
    dwave->data[j]  = (REAL8)(dy / dt);

  }


  /* Free gsl variables */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  LALFree(x);
  LALFree(y);

  return errcode;
}

INT4 XLALPSpinGenerateQNMFreq(
	COMPLEX8Vector		*modefreqs,
	InspiralTemplate	*params,
	UINT4			l,
	INT4			m,
	UINT4			nmodes,
	REAL8                   finalMass,
	REAL8                   finalSpin
	)

{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  UINT4 i;
  REAL8 totalMass;
  /* Fitting coefficients for QNM frequencies from PRD73, 064030, gr-qc/0512160, tables VIII and IX */
  REAL4 BCW22re[3][3]  = { {1.5251, -1.1568,  0.1292}, {1.3673, -1.0260,  0.1628}, { 1.3223, -1.0257,  0.1860} };
  REAL4 BCW22im[3][3]  = { {0.7000,  1.4187, -0.4990}, {0.1000,  0.5436, -0.4731}, {-0.1000,  0.4206, -0.4256} };

  /*REAL4 BCW2m2re[3][3] = { {0.2938,  0.0782,  1.3546}, {0.2528,  0.0921,  1.3344}, { 0.1873,  0.1117,  1.3322} };
    REAL4 BCW2m2im[3][3] = { {1.6700,  0.4192,  1.4700}, {0.4550,  0.1729,  1.3617}, { 0.1850,  0.1266,  1.3661} };*/

  REAL4 BCW21re[3][3]  = { {0.60000, -0.2339, 0.4175}, {0.5800, -0.2416, 0.4708}, { 0.5660, -0.2740, 0.4960} };
  REAL4 BCW21im[3][3]  = { {-0.30000, 2.3561, -0.2277}, {-0.3300, 0.9501, -0.2072}, { -0.1000, 0.4173, -0.2774} };

  /*REAL4 BCW2m1re[3][3] = { {0.3441, 0.0293, 2.0010}, {0.3165, 0.0301, 2.3415}, {0.2696, 0.0315, 2.7755} };
    REAL4 BCW2m1im[3][3] = { {2.0000, 0.1078, 5.0069}, {0.6100, 0.0276, 13.1683}, {0.2900, 0.0276, 6.4715} };*/

  REAL4 BCW20re[3][3]  = { {0.4437, -0.0739,  0.3350}, {0.4185, -0.0768,  0.4355}, { 0.3734, -0.0794,  0.6306} };
  REAL4 BCW20im[3][3]  = { {4.0000,  -1.9550, 0.1420}, {1.2500,  -0.6359, 0.1614}, {0.5600,  -0.2589, -0.3034} };

  REAL4 BCW33re[3][3]  = { {1.8596, -1.3043, 0.1818}, {1.8566, -1.2818, 0.1934}, {1.8004, -1.2558, 0.2133} };
  REAL4 BCW33im[3][3]  = { {0.9000, 2.3430, -0.4810}, {0.2274, 0.8173, -0.4731}, {0.0400, 0.5445, -0.4539} };

  /*REAL4 BCW3m3re[3][3] = { {0.4673, 0.1296, 1.3255}, {0.4413, 0.1387, 1.3178}, {0.3933, 0.1555, 1.3037} };
    REAL4 BCW3m3im[3][3] = { {2.5500, 0.6576, 1.3378}, {0.7900, 0.2381, 1.3706}, {0.4070, 0.1637, 1.3819} };*/

  REAL4 BCW32re[3][3]  = { {1.1481, -0.5552, 0.3002}, {1.1226, -0.5471, 0.3264}, {1.0989, -0.5550, 0.3569} };
  REAL4 BCW32im[3][3]  = { {0.8313, 2.3773, -0.3655}, {0.2300, 0.8025, -0.3684}, {0.1000, 0.4804, -0.3784}};

  /*REAL4 BCW3m2re[3][3] = { {0.5158, 0.8195, 1.408}, {0.4413, 0.1378, 1.3178}, {0.4567, 0.09300, 1.4469} };
    REAL4 BCW3m2im[3][3] = { {2.9000, 0.3365, 2.3050}, {0.9000, 0.1295, 1.6142}, {0.4900, 0.0848, 1.9737} };*/

  REAL4 BCW31re[3][3]  = { {0.8345, -0.2405, 0.4095}, {0.8105, -0.2342, 0.4660}, {0.7684, -0.2252, 0.5805} };
  REAL4 BCW31im[3][3]  = { {23.8450, -20.724, 0.03837}, {8.8530, -7.8506, 0.03418}, {2.1800, -1.6273, 0.1163} };

  /*REAL4 BCW3m1re[3][3] = { {0.5751, 0.02508, 3.1360}, {0.5584, 0.02514, 3.4154}, {0.5271, 0.02561, 3.8011} };
    REAL4 BCW3m1im[3][3] = { {3.0464, 0.1162, -0.2812}, {1.2000, -0.1928, 0.1037}, {1.0000, -0.4424, 0.02467} };*/

  REAL4 BCW30re[3][3]  = { {0.6873, -0.09282, 0.3479}, {0.6687, -0.09155, 0.4021}, {0.6343, -0.08915, 0.5117} };
  REAL4 BCW30im[3][3]  = { {6.7841, -3.6112, 0.09480}, {2.0075, -0.9930, 0.1197}, {0.9000, -0.3409, 0.2679} };

  REAL4 BCW44re[3][3]  = { {2.3, -1.5056, 0.2244}, {2.3, -1.5173, 0.2271}, {2.3, -1.5397, 0.2321} };
  REAL4 BCW44im[3][3]  = { {1.1929, 3.1191, -0.4825}, {0.3, 1.1034, -0.4703}, {0.11, 0.6997, -0.4607} };

  /*REAL4 BCW4m4re[3][3]  = { {0.6256, 0.18, 1.3218}, {0.6061, 0.1869, 1.3168}, {0.5686, 0.2003, 1.3068} };
    REAL4 BCW4m4im[3][3]  = { {3.4, 0.8696, 1.4074}, {1.08, 0.3095, 1.3279}, {0.5980, 0.2015, 1.3765} };*/

  REAL4 BCW43re[3][3] = { {1.6869, -0.8862, 0.2822}, {1.6722, -0.8843, 0.2923}, {1.6526, -0.8888, 0.3081} };
  REAL4 BCW43im[3][3] = { {1.4812, 2.8096, -0.4271}, {0.4451, 0.9569, -0.425}, {0.22, 0.5904, -0.4236} };

  /*REAL4 BCW4m3re[3][3] = { {0.6728, 0.1338, 1.3413}, {0.6562, 0.1377, 1.3456}, {0.6244, 0.1454, 1.3513} };
    REAL4 BCW4m3im[3][3] = { {3.7, 0.5829, 1.6681}, {1.18, 0.2111, 1.4129}, {0.66, 0.1385, 1.3742} };*/

  REAL4 BCW42re[3][3]  = { {1.2702, -0.4685, 0.3835}, {1.2462, -0.4580, 0.4139}, {1.2025, -0.4401, 0.4769} };
  REAL4 BCW42im[3][3]  = { {-3.6, 7.7749, -0.1491}, {-1.5, 2.8601, -0.1392}, {-1.5, 2.2784, -0.1124}};

  /*REAL4 BCW4m2re[3][3] = { {0.7294, 0.07842, 1.5646}, {0.7154, 0.07979, 1.5852}, {0.6885, 0.08259, 1.6136} };
    REAL4 BCW4m2im[3][3] = { {4., 0.2777, 2.0647}, {1.32, 0.08694, 4.3255}, {0.75, 0.05803, 3.7971} };*/

  REAL4 BCW41re[3][3]  = { {1.0507, -0.2478, 0.4348}, {1.0337, -0.2439, 0.4695}, {1.0019, -0.2374, 0.5397} };
  REAL4 BCW41im[3][3]  = { {14., -9.8240, 0.09047}, {4.2, -2.8399, 0.1081}, {2.2, -1.4195, 0.1372} };

  /*REAL4 BCW4m1re[3][3] = { {0.7908, 0.02024, 5.4628}, {0.7785, 0.02005, 5.8547}, {0.7549, 0.01985, 6.5272} };
    REAL4 BCW4m1im[3][3] = { {4.6, -0.4038, 0.4629}, {1.6, -0.2323, 0.2306}, {1.6, -0.8136, 0.03163} };*/

  REAL4 BCW40re[3][3]  = { {0.9175, -0.1144, 0.3511}, {0.9028, -0.1127, 0.3843}, {0.8751, -0.1096, 0.4516} };
  REAL4 BCW40im[3][3]  = { {7.0, -2.7934, 0.1708}, {2.2, -0.8308, 0.2023}, {1.2, -0.4159, 0.2687} };


  /* Get a local copy of the intrinstic parameters */
  totalMass = params->totalMass;

  /* QNM frequencies from the fitting given in PRD73, 064030 */

  if ((l==2)&&(abs(m)==2)) {
    for (i = 0; i < nmodes; i++)
      {
	modefreqs->data[i] = crectf( BCW22re[i][0] + BCW22re[i][1] * pow(1.- finalSpin, BCW22re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW22im[i][0] + BCW22im[i][1] * pow(1.- finalSpin, BCW22im[i][2])) );
	modefreqs->data[i] *= ((REAL4) 1./ finalMass / (totalMass * LAL_MTSUN_SI));
      }
  }
  else {
    if ((l==2)&&(m==0)) {
      for (i = 0; i < nmodes; i++)
	{
	  modefreqs->data[i] = crectf( BCW20re[i][0] + BCW20re[i][1] * pow(1.- finalSpin, BCW20re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW20im[i][0] + BCW20im[i][1] * pow(1.- finalSpin, BCW20im[i][2])) );
	  modefreqs->data[i] /= ((REAL4) finalMass * totalMass * LAL_MTSUN_SI);
	}
    }
    else {
      if ((l==2)&&(abs(m)==1)) {
	for (i = 0; i < nmodes; i++) {
	  modefreqs->data[i] = crectf( BCW21re[i][0] + BCW21re[i][1] * pow(1.- finalSpin, BCW21re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW21im[i][0] + BCW21im[i][1] * pow(1.- finalSpin, BCW21im[i][2])) );
	  modefreqs->data[i] /= ((REAL4) finalMass * totalMass * LAL_MTSUN_SI);
	}
      }
      else {
	if ((l==3)&&(abs(m)==3)) {
	  for (i = 0; i < nmodes; i++) {
	    modefreqs->data[i] = crectf( BCW33re[i][0] + BCW33re[i][1] * pow(1.- finalSpin, BCW33re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW33im[i][0] + BCW33im[i][1] * pow(1.- finalSpin, BCW33im[i][2])) );
	    modefreqs->data[i] /= ((REAL4) finalMass * totalMass * LAL_MTSUN_SI);
	  }
	}
	else
	  if ((l==3)&&(abs(m)==2)) {
	    for (i = 0; i < nmodes; i++) {
	      modefreqs->data[i] = crectf( BCW32re[i][0] + BCW32re[i][1] * pow(1.- finalSpin, BCW32re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW32im[i][0] + BCW32im[i][1] * pow(1.- finalSpin, BCW32im[i][2])) );
	      modefreqs->data[i] /= ((REAL4) finalMass * totalMass * LAL_MTSUN_SI);
	    }
	  }
	  else {
	    if ((l==3)&&(abs(m)==1)) {
	      for (i = 0; i < nmodes; i++) {
		modefreqs->data[i] = crectf( BCW31re[i][0] + BCW31re[i][1] * pow(1.- finalSpin, BCW31re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW31im[i][0] + BCW31im[i][1] * pow(1.- finalSpin, BCW31im[i][2])) );
		modefreqs->data[i] /= ((REAL4) finalMass * totalMass * LAL_MTSUN_SI);
	      }
	    }
	    else {
	      if ((l==3)&&(m==0)) {
		for (i = 0; i < nmodes; i++) {
		  modefreqs->data[i] = crectf( BCW30re[i][0] + BCW30re[i][1] * pow(1.- finalSpin, BCW30re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW30im[i][0] + BCW30im[i][1] * pow(1.- finalSpin, BCW30im[i][2])) );
		  modefreqs->data[i] /= ((REAL4) finalMass * totalMass * LAL_MTSUN_SI);
		}
	      }
	      else {
		if ((l==4)&&(abs(m)==4)) {
		  for (i = 0; i < nmodes; i++) {
		    modefreqs->data[i] = crectf( BCW44re[i][0] + BCW44re[i][1] * pow(1.- finalSpin, BCW44re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW44im[i][0] + BCW44im[i][1] * pow(1.- finalSpin, BCW44im[i][2])) );
		    modefreqs->data[i] /= ((REAL4) finalMass * totalMass * LAL_MTSUN_SI);
		  }
		}
		else {
		  if ((l==4)&&(abs(m)==3)) {
		    for (i = 0; i < nmodes; i++) {
		      modefreqs->data[i] = crectf( BCW43re[i][0] + BCW43re[i][1] * pow(1.- finalSpin, BCW43re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW43im[i][0] + BCW43im[i][1] * pow(1.- finalSpin, BCW43im[i][2])) );
		      modefreqs->data[i] /= ((REAL4) finalMass * totalMass * LAL_MTSUN_SI);
		    }
		  }
		  else {
		    if ((l==4)&&(abs(m)==2)) {
		      for (i = 0; i < nmodes; i++) {
			modefreqs->data[i] = crectf( BCW42re[i][0] + BCW42re[i][1] * pow(1.- finalSpin, BCW42re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW42im[i][0] + BCW42im[i][1] * pow(1.- finalSpin, BCW42im[i][2])) );
			modefreqs->data[i] /= ((REAL4) finalMass * totalMass * LAL_MTSUN_SI);
		      }
		    }
		    else {
		      if ((l==4)&&(abs(m)==1)) {
			for (i = 0; i < nmodes; i++) {
			  modefreqs->data[i] = crectf( BCW41re[i][0] + BCW41re[i][1] * pow(1.- finalSpin, BCW41re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW41im[i][0] + BCW41im[i][1] * pow(1.- finalSpin, BCW41im[i][2])) );
			  modefreqs->data[i] /= ((REAL4) finalMass * totalMass * LAL_MTSUN_SI);
			}
		      }
		      else {
			if ((l==4)&&(m==0)) {
			  for (i = 0; i < nmodes; i++) {
			    modefreqs->data[i] = crectf( BCW40re[i][0] + BCW40re[i][1] * pow(1.- finalSpin, BCW40re[i][2]), crealf(modefreqs->data[i]) / 2. / (BCW40im[i][0] + BCW40im[i][1] * pow(1.- finalSpin, BCW40im[i][2])) );
			    modefreqs->data[i] /= ((REAL4) finalMass * totalMass * LAL_MTSUN_SI);
			  }
			}
			else {
			  fprintf(stderr,"*** LALPSpinInspiralRingdownWave ERROR: Ringdown modes for l=%d m=%d not availbale\n",l,m);
			  XLAL_ERROR( XLAL_EDOM );
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
      }
    }
  }

  return errcode;
}

INT4 XLALPSpinFinalMassSpin(
	REAL8		 *finalMass,
	REAL8		 *finalSpin,
	InspiralTemplate *params,
	REAL8            energy,
	REAL8            *LNhvec
	)
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  REAL8 qq,ll,eta;

  /* See eq.(6) in arXiv:0904.2577 */
  REAL8 ma1,ma2,a12,a12l;
  REAL8 cosa1=0.;
  REAL8 cosa2=0.;
  REAL8 cosa12=0.;

  REAL8 t0=-2.9;
  REAL8 t3=2.6;
  REAL8 s4=-0.123;
  REAL8 s5=0.45;
  REAL8 t2=16.*(0.6865-t3/64.-sqrt(3.)/2.);

  /* get a local copy of the intrinstic parameters */
  qq=params->mass2/params->mass1;
  eta = params->eta;
  /* done */
  ma1=sqrt( params->spin1[0]*params->spin1[0] + params->spin1[1]*params->spin1[1] + params->spin1[2]*params->spin1[2] );
  ma2=sqrt( params->spin2[0]*params->spin2[0] + params->spin2[1]*params->spin2[1] + params->spin2[2]*params->spin2[2] );

  if (ma1>0.) cosa1 = (params->spin1[0]*LNhvec[0]+params->spin1[1]*LNhvec[1]+params->spin1[2]*LNhvec[2])/ma1;
  else cosa1=0.;
  if (ma2>0.) cosa2 = (params->spin2[0]*LNhvec[0]+params->spin2[1]*LNhvec[1]+params->spin2[2]*LNhvec[2])/ma2;
  else cosa2=0.;
  if ((ma1>0.)&&(ma2>0.)) {
    cosa12  = (params->spin1[0]*params->spin2[0] + params->spin1[1]*params->spin2[1] + params->spin1[2]*params->spin2[2])/ma1/ma2;
  }
  else cosa12=0.;

  a12  = ma1*ma1 + ma2*ma2*qq*qq*qq*qq + 2.*ma1*ma2*qq*qq*cosa12 ;
  a12l = ma1*cosa1 + ma2*cosa2*qq*qq ;
  ll = 2.*sqrt(3.)+ t2*eta + t3*eta*eta + s4*a12/(1.+qq*qq)/(1.+qq*qq) + (s5*eta+t0+2.)/(1.+qq*qq)*a12l;

  /* Estimate final mass by adding the negative binding energy to the rest mass*/
  *finalMass = 1. + energy;
  if (*finalMass < 0.) {
    fprintf(stderr,"*** LALPSpinInspiralRingdownWave ERROR: Estimated final mass <0 : %12.6f\n ",*finalMass);
    fprintf(stderr,"***                                    Final mass set to initial mass\n");
    XLAL_ERROR( XLAL_ERANGE);
    *finalMass = 1.;
  }

  *finalSpin = sqrt( a12 + 2.*ll*qq*a12l + ll*ll*qq*qq)/(1.+qq)/(1.+qq);
  if ((*finalSpin > 1.)||(*finalSpin < 0.)) {
    if ((*finalSpin>=1.)&&(*finalSpin<1.01)) {
	fprintf(stderr,"*** LALPSpinInspiralRingdownWave WARNING: Estimated final Spin slightly >1 : %11.3e\n ",*finalSpin);
	fprintf(stderr,"      (m1=%8.3f  m2=%8.3f s1=(%8.3f,%8.3f,%8.3f) s2=(%8.3f,%8.3f,%8.3f) ) final spin set to 1 and code goes on\n",params->mass1,params->mass2,params->spin1[0],params->spin1[1],params->spin1[2],params->spin2[0],params->spin2[1],params->spin2[2]);
	*finalSpin = .99999;
      }
    else {
      fprintf(stderr,"*** LALPSpinInspiralRingdownWave ERROR: Unphysical estimation of final Spin : %11.3e\n ",*finalSpin);
     fprintf(stderr,"      (m1=%8.3f  m2=%8.3f s1=(%8.3f,%8.3f,%8.3f) s2=(%8.3f,%8.3f,%8.3f) )\n",params->mass1,params->mass2,params->spin1[0],params->spin1[1],params->spin1[2],params->spin2[0],params->spin2[1],params->spin2[2]); 
     fprintf(stderr,"***                                    Code aborts\n");
      *finalSpin = 0.;
      XLAL_ERROR( XLAL_ERANGE);
    }
  }

  /*For reference these are the formula used in the EOBNR construction*/
  //*finalMass = 1. - 0.057191 * eta - 0.498 * eta*eta;
  //*finalSpin = 3.464102 * eta - 2.9 * eta*eta;

  return errcode;
}

INT4 XLALPSpinInspiralAttachRingdownWave (
      REAL8Vector 	*sigl,
      InspiralTemplate 	*params,
      UINT4             *attpos,
      UINT4              nmodes,
      UINT4              l,
      INT4               m,
      REAL8              finalMass,
      REAL8              finalSpin
    )
{
      const UINT4 Npatch=40;
      const UINT4 offsetAttch = 2;

      COMPLEX8Vector *modefreqs;
      UINT4 Nrdwave;

      UINT4 i=0;
      UINT4 j=0;
      UINT4 k=0;
      UINT4 atpos;
      INT4 errcode;

      REAL8Vector	*rdwave;
      REAL8Vector	*inspwave,*dinspwave;
      REAL8Vector	*matchinspwave;
      REAL8 dt;

      dt = 1./params->tSampling;
      atpos=(*attpos);

      /* Create memory for the QNM frequencies */
      modefreqs = XLALCreateCOMPLEX8Vector( nmodes );
      if ( !modefreqs )
      {
        XLAL_ERROR( XLAL_ENOMEM );
      }
      errcode = XLALPSpinGenerateQNMFreq( modefreqs, params, l, m, nmodes, finalMass, finalSpin);
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( XLAL_EFUNC );
      }

      /* Ringdown signal length: 10 times the decay time of the n=0 mode */
      Nrdwave = (INT4) (10. / cimagf(modefreqs->data[0]) / dt);
      /* Patch length, centered around the matching point "attpos" */

      (*attpos)+=Nrdwave;

      /* Check the value of attpos, to prevent memory access problems later */
      if ( atpos < Npatch || atpos + Npatch >= sigl->length )
      {
        XLALPrintError( "Value of attpos inconsistent with given value of Npatch: atpos=%d  Npatch=%d, sign->length=%d, m1=%11.5f  m2=%11.5f  s1z=%8.3f  s2z=%8.3f  fL=%11.3e\n",atpos,Npatch,sigl->length,params->mass1,params->mass2,params->spin1[2],params->spin2[2],params->fLower);
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( XLAL_EFAILED );
      }

      /* Create memory for the ring-down and full waveforms, derivatives of inspirals 
	 and waveforms and its derivative values at the attach point */

      rdwave = XLALCreateREAL8Vector( Nrdwave );
      inspwave = XLALCreateREAL8Vector( Npatch );
      dinspwave = XLALCreateREAL8Vector( Npatch );
      matchinspwave = XLALCreateREAL8Vector( 2*nmodes );

      /* Check memory was allocated */
      if ( !rdwave || !inspwave || !dinspwave || !matchinspwave )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        if (rdwave)         XLALDestroyREAL8Vector( rdwave );
        if (inspwave)       XLALDestroyREAL8Vector( inspwave );
        if (dinspwave)      XLALDestroyREAL8Vector( dinspwave );
        if (matchinspwave) XLALDestroyREAL8Vector( matchinspwave );
        XLAL_ERROR( XLAL_ENOMEM );
      }

      /* Generate derivatives of the last part of inspiral waves */
      /* Take the last part of sigl1 */

      for (i=0; i<2; i++) {
	/* i=0(1) for real(imaginary) part */
	for (j = 0; j < Npatch; j++) {
	  inspwave->data[j]    = sigl->data[2*(atpos - Npatch + j)+i];
	}

	for (k=0;k<2*nmodes;k++) {
	  matchinspwave->data[k] = inspwave->data[Npatch-1-offsetAttch];
	  if ((k+1)<2*nmodes) {
	    errcode = XLALGenerateWaveDerivative( dinspwave, inspwave, dt);
	    if ( (errcode != XLAL_SUCCESS) ) {
	      XLALDestroyCOMPLEX8Vector( modefreqs );
	      XLALDestroyREAL8Vector( rdwave );
	      XLALDestroyREAL8Vector( inspwave );
	      XLALDestroyREAL8Vector( dinspwave );
	      XLALDestroyREAL8Vector( matchinspwave );
	      XLAL_ERROR( XLAL_EFUNC );
	    }
	    for (j=0; j<Npatch; j++) {
	      inspwave->data[j]=dinspwave->data[j];
	    }
	  }
	}

	errcode = XLALPSpinInspiralRingdownWave( rdwave, params, matchinspwave, modefreqs, nmodes );

	if ( errcode != XLAL_SUCCESS ) {
	  XLALDestroyCOMPLEX8Vector( modefreqs );
	  XLALDestroyREAL8Vector( rdwave );
	  XLALDestroyREAL8Vector( inspwave );
	  XLALDestroyREAL8Vector( dinspwave );
	  XLALDestroyREAL8Vector( matchinspwave );
	  XLAL_ERROR( XLAL_EFUNC );
	}
	/* Generate full waveforms, by stitching inspiral and ring-down waveforms */

	for (j = 0; j < Nrdwave; j++) {
	  sigl->data[2*j + 2*(atpos - 1 - offsetAttch) + i ] = rdwave->data[j];
	}

      }

      /* Free memory */
      XLALDestroyCOMPLEX8Vector( modefreqs );
      XLALDestroyREAL8Vector( rdwave );
      XLALDestroyREAL8Vector( inspwave );
      XLALDestroyREAL8Vector( dinspwave );
      XLALDestroyREAL8Vector( matchinspwave );

      return errcode;
}
