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


/**** <lalVerbatim file="LALPSpinInspiralRingdownWaveCV">
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Module \texttt{LALPSpinInspiralRingdownWave.c}}
 * 
 * Module to compute the ring-down waveform as linear combination
 * of quasi-normal-modes decaying waveforms, which can be attached to
 * the phenomenological spin Taylor waveform.
 * 
 * \subsubsection*{Prototypes}
 * \vspace{0.1in}
 * \input{XLALXLALPSpinInspiralRingdownWaveCP}
 * \idx{XLALXLALPSpinInspiralRingdownWave()}
 * \begin{itemize}
 * \item \texttt{rdwave,} Output, the ring-down waveform
 * \item \texttt{params,} Input, the parameters where ring-down waveforms are computed
 * \item \texttt{inspwave,} Input, the inspiral waveform with given multiple
 * \item \texttt{modefreqs,} Input, the frequencies of the quasi-normal-modes
 * \item \texttt{nmodes,} Input, the number of quasi-normal-modes to be combined.
 * \end{itemize}
 *
 * \input{XLALPSpinGenerateWaveDerivativeCP}
 * \idx{XLALPSpinGenerateWaveDerivative()}
 * \begin{itemize}
 * \item \texttt{dwave,} Output, time derivative of the input waveform
 * \item \texttt{wave,} Input, waveform to be differentiated in time
 * \item \texttt{params,} Input, the parameters of the input waveform.
 * \end{itemize}
 *
 * \input{XLALPSpinGenerateQNMFreqCP}
 * \idx{XLALPSpinGenerateQNMFreq()}
 * \begin{itemize}
 * \item \texttt{ptfwave,} Output, the frequencies of the quasi-normal-modes
 * \item \texttt{params,} Input, the parameters of the binary system
 * \item \texttt{l,} Input, the l of the modes
 * \item \texttt{m,} Input, the m of the modes
 * \item \texttt{nmodes,} Input, the number of overtones considered.
 * \end{itemize}
 *
 * \input{XLALPSpinFinalMassSpinCP}
 * \idx{XLALPSpinFinalMassSpin()}
 * \begin{itemize}
 * \item \texttt{finalMass,} Output, the mass of the final Kerr black hole
 * \item \texttt{finalSpin,}  Input, the spin of the final Kerr balck hole
 * \item \texttt{params,} Input, the parameters of the binary system.
 * \item \texttt{energy,} Input, the binding energy at the time final time.
 * \end{itemize}

 * \subsubsection*{Description}
 * This module generate ring-down waveforms.
 *
 * \subsubsection*{Algorithm}
 *
 * \subsubsection*{Uses}
 * \begin{verbatim}
 * LALMalloc
 * LALFree
 * \end{verbatim}
 *
 * \subsubsection*{Notes}
 *
 * \vfill{\footnotesize\input{LALPSpinInspiralRingdownWaveCV}}
 *
 **** </lalLaTeX> */

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

/* <lalVerbatim file="XLALInspiralRingdownWaveCP">  */
INT4 XLALPSpinInspiralRingdownWave (
	REAL4Vector		*rdwave,
	InspiralTemplate	*params,
	REAL4Vector	        *matchinspwave,
	COMPLEX8Vector		*modefreqs,
	UINT4			nmodes
	)
/* </lalVerbatim> */
{

  static const char *func = "XLALPSpinInspiralRingdownWave";

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
      XLAL_ERROR( func, XLAL_EBADLEN );
    }

  /* Solving the linear system for QNMs amplitude coefficients using gsl routine */
  /* Initialize matrices and supporting variables */
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
  /* Matrix A (2*nmodes by 2*nmodes) has block symmetry. Define half of A here as "coef" */
  /* Define y here as "hderivs" */

  j=0;
  while (j<nmodes) {
    if (j==0) {
      for (i = 0; i < nmodes; i++) {
	gsl_matrix_set(coef, j, i, 1.);
	gsl_matrix_set(coef, j, i+nmodes, 0.);
	gsl_matrix_set(coef, j+1, i, -modefreqs->data[i].im);
	gsl_matrix_set(coef, j+1, i+nmodes, modefreqs->data[i].re);
      }
    }
    else {
      if (j==1) {
	for (i = 0; i < nmodes; i++) {
	  gsl_matrix_set(coef, 2*j, i, modefreqs->data[i].im*modefreqs->data[i].im-modefreqs->data[i].re*modefreqs->data[i].re);
	  gsl_matrix_set(coef, 2*j, i+nmodes, -2.*modefreqs->data[i].im*modefreqs->data[i].re);
	  gsl_matrix_set(coef, 2*j+1, i, -modefreqs->data[i].im*modefreqs->data[i].im*modefreqs->data[i].im+3.*modefreqs->data[i].im*modefreqs->data[i].re*modefreqs->data[i].re);
	  gsl_matrix_set(coef, 2*j+1, i+nmodes, -modefreqs->data[i].re*modefreqs->data[i].re*modefreqs->data[i].re+3*modefreqs->data[i].re*modefreqs->data[i].im*modefreqs->data[i].im);
	}
      }
      else {
	fprintf(stderr,"*** ERROR: nmode must be <=2, %d selected\n",nmodes);
	XLAL_ERROR( func, XLAL_EDOM );
      }
    }
    gsl_vector_set(hderivs, 2*j, matchinspwave->data[2*j]);
    gsl_vector_set(hderivs, 2*j+1, matchinspwave->data[2*j+1]);
    j++;
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
  modeamps = XLALCreateREAL8Vector(2*nmodes);

  if ( !modeamps )
  {
    gsl_matrix_free(coef);
    gsl_vector_free(hderivs);
    gsl_vector_free(x);
    gsl_permutation_free(p);
    XLAL_ERROR( func, XLAL_ENOMEM );
  }

  for (i = 0; i < 2*nmodes; ++i) {
    modeamps->data[i] = gsl_vector_get(x, i);
  }

  /* Free all gsl linear algebra objects */
  gsl_matrix_free(coef);
  gsl_vector_free(hderivs);
  gsl_vector_free(x);
  gsl_permutation_free(p);

  /* Build ring-down waveforms */
  UINT4 Nrdwave=rdwave->length;
  for (j = 0; j < Nrdwave; j++)
    {
      tj = j * dt;
      rdwave->data[j] = 0.;
      for (i = 0; i < nmodes; i++)
	{
	  rdwave->data[j] += exp(- tj * modefreqs->data[i].im)
	    * ( modeamps->data[i] * cos(tj * modefreqs->data[i].re)
		+   modeamps->data[i + nmodes] * sin(tj * modefreqs->data[i].re) );

	}

    }

  XLALDestroyREAL8Vector(modeamps);
  return errcode;
} /*End of XLALPSpinInspiralRingdownWave */

/* <lalVerbatim file="XLALPSpinGenerateWaveDerivative">  */
INT4 XLALPSpinGenerateWaveDerivative (
	REAL4Vector		*dwave,
	REAL4Vector	        *wave,
	InspiralTemplate	*params
	)
/* </lalVerbatim> */
{
  static const char *func = "XLALPSpinGenerateWaveDerivative";

  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;

  /* For checking GSL return codes */
  INT4 gslStatus;

  UINT4 j;
  REAL8 dt;
  double *x, *y;
  double dy;
  gsl_interp_accel *acc;
  gsl_spline *spline;

  /* Sampling rate from input */
  dt = 1.0 / params -> tSampling;

  /* Getting interpolation and derivatives of the waveform using gsl spline routine */
  /* Initialize arrays and supporting variables for gsl */

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
    if (gslStatus != GSL_SUCCESS )
    {
      gsl_spline_free(spline);
      gsl_interp_accel_free(acc);
      LALFree( x );
      LALFree( y );
      XLAL_ERROR( func, XLAL_EFUNC );
    }
    dwave->data[j]  = (REAL4)(dy / dt);

  }


  /* Free gsl variables */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  LALFree(x);
  LALFree(y);

  return errcode;
}

/* <lalVerbatim file="XLALPSpinGenerateQNMFreqCP">  */
INT4 XLALPSpinGenerateQNMFreq(
	COMPLEX8Vector		*modefreqs,
	InspiralTemplate	*params,
	UINT4			l,
	INT4			m,
	UINT4			nmodes,
	REAL8                   finalMass,
	REAL8                   finalSpin
	)

/* </lalVerbatim> */
{

  static const char *func = "XLALPSpinGenerateQNMFreq";

  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  UINT4 i;
  REAL8 totalMass;
  /* Fitting coefficients for QNM frequencies from PRD73, 064030, gr-qc/0512160 */
  REAL4 BCW22re[3][3]  = { {1.5251, -1.1568,  0.1292}, {1.3673, -1.0260,  0.1628}, { 1.3223, -1.0257,  0.1860} };
  REAL4 BCW22im[3][3]  = { {0.7000,  1.4187, -0.4990}, {0.1000,  0.5436, -0.4731}, {-0.1000,  0.4206, -0.4256} };

  REAL4 BCW2m2re[3][3] = { {0.2938,  0.0782,  1.3546}, {0.2528,  0.0921,  1.3344}, { 0.1873,  0.1117,  1.3322} };
  REAL4 BCW2m2im[3][3] = { {1.6700,  0.4192,  1.4700}, {0.4550,  0.1729,  1.3617}, { 0.1850,  0.1266,  1.3661} };

  REAL4 BCW21re[3][3]  = { {0.60000, -0.2339, 0.4175}, {0.5800, -0.2416, 0.4708}, { 0.5660, -0.2740, 0.4960} };
  REAL4 BCW21im[3][3]  = { {-0.30000, 2.3561, -0.2277}, {-0.3300, 0.9501, -0.2072}, { -0.1000, 0.4173, -0.2774} };

  REAL4 BCW2m1re[3][3] = { {0.3441, 0.0293, 2.0010}, {0.3165, 0.0301, 2.3415}, {0.2696, 0.0315, 2.7755} };
  REAL4 BCW2m1im[3][3] = { {2.0000, 0.1078, 5.0069}, {0.6100, 0.0276, 13.1683}, {0.2900, 0.0276, 6.4715} };

  REAL4 BCW20re[3][3]  = { {0.4437, -0.0739,  0.3350}, {0.4185, -0.0768,  0.4355}, { 0.3734, -0.0794,  0.6306} };
  REAL4 BCW20im[3][3]  = { {4.0000,  -1.9550, 0.1420}, {1.2500,  -0.6359, 0.1614}, {0.5600,  -0.2589, -0.3034} };

  REAL4 BCW33re[3][3]  = { {1.8596, -1.3043, 0.1818}, {1.8566, -1.2818, 0.1934}, {1.8004, -1.2558, 0.2133} };
  REAL4 BCW33im[3][3]  = { {0.9000, 2.3430, -0.4810}, {0.2274, 0.8173, -0.4731}, {0.0400, 0.5445, -0.4539} };

  REAL4 BCW3m3re[3][3] = { {0.4673, 0.1296, 1.3255}, {0.4413, 0.1387, 1.3178}, {0.3933, 0.1555, 1.3037} };
  REAL4 BCW3m3im[3][3] = { {2.5500, 0.6576, 1.3378}, {0.7900, 0.2381, 1.3706}, {0.4070, 0.1637, 1.3819} };

  /*
  REAL4 BCW32re[3][3]  = { {1.1481, -0.5552, 0.3002}, {1.1226, -0.5471, 0.3264}, {1.0989, -0.5550, 0.3569} };
  REAL4 BCW32im[3][3]  = { {0.8313, 2.3773, -0.3655}, {0.2300, 0.8025, -0.3684}, {0.1000, 0.4804, -0.3784}};

  REAL4 BCW3m2re[3][3] = { {0.5158, 0.8195, 1.408}, {0.4413, 0.1378, 1.3178}, {0.4567, 0.09300, 1.4469} };
  REAL4 BCW3m2im[3][3] = { {2.9000, 0.3365, 2.3050}, {0.9000, 0.1295, 1.6142}, {0.4900, 0.0848, 1.9737} };

  */

  /* Get a local copy of the intrinstic parameters */
  totalMass = params->totalMass;

  /* QNM frequencies from the fitting given in PRD73, 064030 */
  /* Other modes than l=2, m=\pm 2, 0 should be added*/

  if ((l==2)&&(m==2)) {
    for (i = 0; i < nmodes; i++)
      {
	modefreqs->data[i].re = BCW22re[i][0] + BCW22re[i][1] * pow(1.- finalSpin, BCW22re[i][2]);
	modefreqs->data[i].im = modefreqs->data[i].re / 2.
	  / (BCW22im[i][0] + BCW22im[i][1] * pow(1.- finalSpin, BCW22im[i][2]));
	modefreqs->data[i].re *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
	modefreqs->data[i].im *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
      }
  }
  else {
    if ((l==2)&&(m==-2)) {
      for (i = 0; i < nmodes; i++)
	{
	  modefreqs->data[i].re = BCW2m2re[i][0] + BCW2m2re[i][1] * pow(1.- finalSpin, BCW2m2re[i][2]);
	  modefreqs->data[i].im = modefreqs->data[i].re / 2.
	    / (BCW2m2im[i][0] + BCW2m2im[i][1] * pow(1.- finalSpin, BCW2m2im[i][2]));
	  modefreqs->data[i].re *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
	  modefreqs->data[i].im *= 1./ finalMass / (totalMass * LAL_MTSUN_SI);
	}
    }
    else {
      if ((l==2)&&(m==0)) {
	for (i = 0; i < nmodes; i++)
	  {
	    modefreqs->data[i].re = BCW20re[i][0] + BCW20re[i][1] * pow(1.- finalSpin, BCW20re[i][2]);
	    modefreqs->data[i].im = modefreqs->data[i].re / 2.
	      / (BCW20im[i][0] + BCW20im[i][1] * pow(1.- finalSpin, BCW20im[i][2]));
	    modefreqs->data[i].re /= finalMass * totalMass * LAL_MTSUN_SI;
	    modefreqs->data[i].im /= finalMass * totalMass * LAL_MTSUN_SI;
	  }
      }
      else {
	if ((l==2)&&(m==1)) {
	  for (i = 0; i < nmodes; i++) {
	    modefreqs->data[i].re = BCW21re[i][0] + BCW21re[i][1] * pow(1.- finalSpin, BCW21re[i][2]);
	    modefreqs->data[i].im = modefreqs->data[i].re / 2.
	      / (BCW21im[i][0] + BCW21im[i][1] * pow(1.- finalSpin, BCW21im[i][2]));
	    modefreqs->data[i].re /= finalMass * totalMass * LAL_MTSUN_SI;
	    modefreqs->data[i].im /= finalMass * totalMass * LAL_MTSUN_SI;
	  }
	}
	else {
	  if ((l==2)&&(m==-1)) {
	    for (i = 0; i < nmodes; i++) {
	      modefreqs->data[i].re = BCW2m1re[i][0] + BCW2m1re[i][1] * pow(1.- finalSpin, BCW2m1re[i][2]);
	      modefreqs->data[i].im = modefreqs->data[i].re / 2.
		/ (BCW2m1im[i][0] + BCW2m1im[i][1] * pow(1.- finalSpin, BCW2m1im[i][2]));
	      modefreqs->data[i].re /= finalMass * totalMass * LAL_MTSUN_SI;
	      modefreqs->data[i].im /= finalMass * totalMass * LAL_MTSUN_SI;
	    }
	  }
	  else {
	    if ((l==3)&&(m==3)) {
	      for (i = 0; i < nmodes; i++) {
		modefreqs->data[i].re = BCW33re[i][0] + BCW33re[i][1] * pow(1.- finalSpin, BCW33re[i][2]);
		modefreqs->data[i].im = modefreqs->data[i].re / 2.
		  / (BCW33im[i][0] + BCW33im[i][1] * pow(1.- finalSpin, BCW33im[i][2]));
		modefreqs->data[i].re /= finalMass * totalMass * LAL_MTSUN_SI;
		modefreqs->data[i].im /= finalMass * totalMass * LAL_MTSUN_SI;
	      }
	    }
	    else {
	      if ((l==3)&&(m==-3)) {
		for (i = 0; i < nmodes; i++) {
		  modefreqs->data[i].re = BCW3m3re[i][0] + BCW3m3re[i][1] * pow(1.- finalSpin, BCW3m3re[i][2]);
		  modefreqs->data[i].im = modefreqs->data[i].re / 2.
		    / (BCW3m3im[i][0] + BCW3m3im[i][1] * pow(1.- finalSpin, BCW3m3im[i][2]));
		  modefreqs->data[i].re /= finalMass * totalMass * LAL_MTSUN_SI;
		  modefreqs->data[i].im /= finalMass * totalMass * LAL_MTSUN_SI;
		}
	      }
	      else {
		fprintf(stderr,"*** LALPSpinInspiralRingdownWave ERROR: Ringdown modes for l=%d m=%d not availbale\n",l,m);
		XLAL_ERROR( func , XLAL_EDOM );
	      }
	    }
	  }
	}
      }
    }
  }

  return errcode;
}

/* <lalVerbatim file="XLALPSpinFinalMassSpinCP">  */
INT4 XLALPSpinFinalMassSpin(
	REAL8		 *finalMass,
	REAL8		 *finalSpin,
	InspiralTemplate *params,
	REAL8            energy,
	REAL8            *LNhvec
	)
/* </lalVerbatim> */
{

  static const char *func = "XLALPSpinFinalMassSpin";

  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  REAL8 qq,ll,eta, eta2;

  /* See eq.(6) in arXiv:0904.2577 */
  REAL8 ma1,ma2,a12,a12l;
  REAL8 cosa1=0.;
  REAL8 cosa2=0.;
  REAL8 cosa12=0.;
  REAL8 unitHz;

  REAL8 t0=-2.9;
  REAL8 t3=2.6;
  REAL8 s4=-0.123;
  REAL8 s5=0.45;
  REAL8 t2=16.*(0.6865-t3/64.-sqrt(3.)/2.);

  /* get a local copy of the intrinstic parameters */
  qq=params->mass2/params->mass1;
  unitHz = params->totalMass * LAL_MTSUN_SI * (REAL8)LAL_PI;
  eta = params->eta;
  eta2 = eta * eta;
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
  ll = 2.*sqrt(3.)+ t2*params->eta + t3*params->eta*params->eta + s4*a12/(1.+qq*qq)/(1.+qq*qq) + (s5*params->eta+t0+2.)/(1.+qq*qq)*a12l;

  /* Estimate final mass by adding the negative binding energy to the rest mass*/
  *finalMass = 1. + energy;
  if (*finalMass < 0.) {
    fprintf(stderr,"*** LALPSpinInspiralRingdownWave ERROR: Estimated final mass <0 : %12.6f\n ",*finalMass);
    fprintf(stderr,"***                                    Final mass set to initial mass\n");
    XLAL_ERROR( func, XLAL_ERANGE);
    *finalMass = 1.;
  }

  *finalSpin = sqrt( a12 + 2.*ll*qq*a12l + ll*ll*qq*qq)/(1.+qq)/(1.+qq);
  if ((*finalSpin > 1.)||(*finalSpin < 0.)) {
    if ((*finalSpin>=1.)&&(*finalSpin<1.005)) {
	fprintf(stderr,"*** LALPSpinInspiralRingdownWave WARNING: Estimated final Spin slightly >1 : %11.3e\n ",*finalSpin);
	fprintf(stderr,"      (m1=%8.3f  m2=%8.3f) final spin set to one and code goes on\n",params->mass1,params->mass2);
	*finalSpin = .99999;
      }
    else {
      fprintf(stderr,"*** LALPSpinInspiralRingdownWave ERROR: Estimated final Spin >1 : %11.3e\n ",*finalSpin);
      fprintf(stderr,"***                                    Final Spin set to 0\n");
      *finalSpin = 0.;
      XLAL_ERROR( func, XLAL_ERANGE);
    }
  }

  /*For reference these are the formula used in the EOBNR construction*/
  //*finalMass = 1. - 0.057191 * eta - 0.498 * eta2;
  //*finalSpin = 3.464102 * eta - 2.9 * eta2;

  return errcode;
}

/* <lalVerbatim file="XLALPSpinInspiralAttachRingdownWaveCP">  */
INT4 XLALPSpinInspiralAttachRingdownWave (
      REAL4Vector 	*sigl,
      InspiralTemplate 	*params,
      UINT4              *attpos,
      UINT4              nmodes,
      UINT4              l,
      INT4               m,
      REAL8              finalMass,
      REAL8              finalSpin
    )
{

      static const char *func = "XLALInspiralAttachRingdownWave";

      COMPLEX8Vector *modefreqs;
      UINT4 Nrdwave, Npatch;
      UINT4 i=0;
      UINT4 j=0;
      UINT4 k=0;
      UINT4 atpos;
      INT4 errcode;

      REAL4Vector	*rdwave;
      REAL4Vector	*inspwave,*dinspwave;
      REAL4Vector	*matchinspwave;
      REAL8 dt;

      dt = 1./params->tSampling;
      atpos=(*attpos);

      /* Create memory for the QNM frequencies */
      modefreqs = XLALCreateCOMPLEX8Vector( nmodes );
      if ( !modefreqs )
      {
        XLAL_ERROR( func, XLAL_ENOMEM );
      }
      errcode = XLALPSpinGenerateQNMFreq( modefreqs, params, l, m, nmodes, finalMass, finalSpin);
      if ( errcode != XLAL_SUCCESS )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( func, XLAL_EFUNC );
      }

      /* Ringdown signal length: 10 times the decay time of the n=0 mode */
      Nrdwave = (INT4) (10. / modefreqs->data[0].im / dt);
      /* Patch length, centered around the matching point "attpos" */

      (*attpos)+=Nrdwave;
      Npatch = 6;

      /* Check the value of attpos, to prevent memory access problems later */
      if ( atpos < Npatch || atpos + Npatch >= sigl->length )
      {
        XLALPrintError( "Value of attpos inconsistent with given value of Npatch: atpos=%d  Npatch=%d, sign->length=%d, m1=%11.5f  m2=%11.5f  s1z=%8.3f  s2z=%8.3f  fL=%11.3e\n",atpos,Npatch,sigl->length,params->mass1,params->mass2,params->spin1[2],params->spin2[2],params->fLower);
        XLALDestroyCOMPLEX8Vector( modefreqs );
        XLAL_ERROR( func, XLAL_EFAILED );
      }

      /* Create memory for the ring-down and full waveforms, and eventual derivatives of inspirals */

      rdwave = XLALCreateREAL4Vector( Nrdwave );
      inspwave = XLALCreateREAL4Vector( Npatch );
      dinspwave = XLALCreateREAL4Vector( Npatch );
      matchinspwave = XLALCreateREAL4Vector( 2*nmodes );

      /* Check memory was allocated */
      if ( !rdwave || !inspwave || !dinspwave || !matchinspwave )
      {
        XLALDestroyCOMPLEX8Vector( modefreqs );
        if (rdwave)         XLALDestroyREAL4Vector( rdwave );
        if (inspwave)       XLALDestroyREAL4Vector( inspwave );
        if (dinspwave)      XLALDestroyREAL4Vector( dinspwave );
        if (matchinspwave) XLALDestroyREAL4Vector( matchinspwave );
        XLAL_ERROR( func, XLAL_ENOMEM );
      }

      /* Generate derivatives of the last part of inspiral waves */
      /* Take the last part of signal1 */

      for (i=0; i<2; i++) {

	for (j = 0; j < Npatch; j++) {
	  inspwave->data[j]    = sigl->data[2*(atpos - Npatch + j)+i];	
	}

	for (k=0;k<2*nmodes;k++) {
	  matchinspwave->data[k] = inspwave->data[Npatch-1];
	  if ((k+1)<2*nmodes) {
	    errcode = XLALPSpinGenerateWaveDerivative( dinspwave, inspwave, params );
	    if ( (errcode != XLAL_SUCCESS) ) {
	      XLALDestroyCOMPLEX8Vector( modefreqs );
	      XLALDestroyREAL4Vector( rdwave );
	      XLALDestroyREAL4Vector( inspwave );
	      XLALDestroyREAL4Vector( dinspwave );
	      XLALDestroyREAL4Vector( matchinspwave );
	      XLAL_ERROR( func, XLAL_EFUNC );
	    }
	    for (j=0; j<Npatch; j++) {
	      inspwave->data[j]=dinspwave->data[j];
	    }
	  }
	}

	errcode = XLALPSpinInspiralRingdownWave( rdwave, params, matchinspwave, modefreqs, nmodes );

	if ( errcode != XLAL_SUCCESS ) {
	  XLALDestroyCOMPLEX8Vector( modefreqs );
	  XLALDestroyREAL4Vector( rdwave );
	  XLALDestroyREAL4Vector( inspwave );
	  XLALDestroyREAL4Vector( dinspwave );
	  XLALDestroyREAL4Vector( matchinspwave );
	  XLAL_ERROR( func, XLAL_EFUNC );
	}
	/* Generate full waveforms, by stitching inspiral and ring-down waveforms */

	for (j = 0; j < Nrdwave; j++) {
	  sigl->data[2*j + 2*atpos - 2 + i] = rdwave->data[j];
	}

      }

      /* Free memory */
      XLALDestroyCOMPLEX8Vector( modefreqs );
      XLALDestroyREAL4Vector( rdwave );
      XLALDestroyREAL4Vector( inspwave );
      XLALDestroyREAL4Vector( dinspwave );
      XLALDestroyREAL4Vector( matchinspwave );

      return errcode;
}
