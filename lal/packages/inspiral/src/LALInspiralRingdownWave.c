/*
*  Copyright (C) 2008 Yi Pan
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

#if 0
<lalVerbatim file="LALInspiralRingdownWaveCV">
Author: Yi Pan
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{LALInspiralRingdownWave.c}}

Module to compute the ring-down waveform as linear combination 
of quasi-normal-modes decaying waveforms, which can be attached to
the inspiral part of the compat binary coalescing waveform.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{XLALXLALInspiralRingdownWaveCP}
\idx{XLALXLALInspiralRingdownWave()}
\begin{itemize}
   \item \texttt{rdwave1,} Output, the real part of the ring-down waveform
   \item \texttt{rdwave2,} Output, the imaginary part of the ring-down waveform
   \item \texttt{params,} Input, the parameters where ring-down waveforms are computed
   \item \texttt{inspwave1,} Input, the real part of the ring-down waveform
   \item \texttt{inspwave2,} Input, the real part of the ring-down waveform
   \item \texttt{modefreqs,} Input, the frequencies of the quasi-normal-modes
   \item \texttt{nmode,} Input, the number of quasi-normal-modes to be combined.
\end{itemize}

\input{XLALGenerateWaveDerivativesCP}
\idx{XLALGenerateWaveDerivatives()}
\begin{itemize}
   \item \texttt{dwave,} Output, time derivative of the input waveform
   \item \texttt{ddwave,} Output, two time derivative of the input waveform
   \item \texttt{wave,} Input, waveform to be differentiated in time
   \item \texttt{params,} Input, the parameters of the input waveform.
\end{itemize}

\input{XLALGenerateQNMFreqCP}
\idx{XLALGenerateQNMFreq()}
\begin{itemize}
   \item \texttt{ptfwave,} Output, the frequencies of the quasi-normal-modes
   \item \texttt{params,} Input, the parameters of the binary system
   \item \texttt{l,} Input, the l of the modes
   \item \texttt{m,} Input, the m of the modes
   \item \texttt{nmodes,} Input, the number of overtones considered.
\end{itemize}

\input{XLALFinalMassSpinCP}
\idx{XLALFinalMassSpin()}
\begin{itemize}
   \item \texttt{finalMass,} Output, the mass of the final Kerr black hole
   \item \texttt{finalSpin,}  Input, the spin of the final Kerr balck hole
   \item \texttt{params,} Input, the parameters of the binary system.
\end{itemize}

\subsubsection*{Description}
Generating ring-down waveforms.

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\begin{verbatim}
LALMalloc
LALFree
\end{verbatim}

\subsubsection*{Notes}
 
\vfill{\footnotesize\input{LALInspiralRingdownWaveCV}}

</lalLaTeX>
#endif

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
INT4 XLALInspiralRingdownWave (
	REAL4Vector				*rdwave1,
	REAL4Vector				*rdwave2,
	InspiralTemplate		*params,
	REAL4VectorSequence		*inspwave1,
	REAL4VectorSequence		*inspwave2,
	COMPLEX8Vector			*modefreqs,
	UINT4					nmodes
	)
/* </lalVerbatim> */
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  static const char* func = "XLALInspiralRingdownWave";
  
  INT4 i, j, k;
  
  /* Sampling rate from input */
  REAL8 dt = 1.0 / params -> tSampling;
  
  if ( inspwave1->length != nmodes || inspwave2->length != nmodes || 
		modefreqs->length != nmodes )
  {
	errcode = 0;
	fprintf( stdout, "Wrong number of waveform derivatives, or frequencies of QNMs.\n");
	exit(1);
  }
  
  /* Solving the linear system for QNMs amplitude coefficients using gsl routine */
  /* Initiate matrices and supporting variables */
  gsl_matrix *coef	  = gsl_matrix_alloc(2 * nmodes, 2 * nmodes);
  gsl_vector *hderivs = gsl_vector_alloc(2 * nmodes);
  gsl_vector *x		  = gsl_vector_alloc(2 * nmodes);
  gsl_permutation *p  = gsl_permutation_alloc(2 * nmodes);
  int s;
  
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
  gsl_linalg_LU_decomp(coef, p, &s);
  gsl_linalg_LU_solve(coef, p, hderivs, x);
  
  /* Putting solution to an XLAL vector */
  REAL8Vector *modeamps = XLALCreateREAL8Vector(2 * nmodes);
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
  REAL8 tj;
  for (j = 0; j < rdwave1->length; ++j)
  {
	tj = j * dt;
	rdwave1->data[j] = 0;
	rdwave2->data[j] = 0;
	for (i = 0; i < nmodes; ++i)
	{
	  rdwave1->data[j] += exp(- tj * modefreqs->data[i].im) 
						* (  modeamps->data[i]			* cos(tj * modefreqs->data[i].re)
						   + modeamps->data[i + nmodes] * sin(tj * modefreqs->data[i].re) );
	  rdwave2->data[j] += exp(- tj * modefreqs->data[i].im) 
						* (- modeamps->data[i]			* sin(tj * modefreqs->data[i].re)
						   + modeamps->data[i + nmodes] * cos(tj * modefreqs->data[i].re) );
	}
  }
  
  return errcode;
}

/* <lalVerbatim file="XLALGenerateWaveDerivatives">  */
INT4 XLALGenerateWaveDerivatives (
	REAL4Vector				*dwave,
	REAL4Vector				*ddwave,
	REAL4Vector				*wave,
	InspiralTemplate		*params
	)
/* </lalVerbatim> */
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  static const char* func = "XLALGenerateQNMFreq";
  
  INT4 j;
  
  /* Sampling rate from input */
  REAL8 dt = 1.0 / params -> tSampling;
  
  /* Getting interpolation and derivatives of the waveform using gsl spline routine */
  /* Initiate arrays and supporting variables for gsl */
  double *x, *y;
  x = (double *) malloc(wave->length * sizeof(double));
  y = (double *) malloc(wave->length * sizeof(double));
  for (j = 0; j < wave->length; ++j)
  {
	x[j] = j;
	y[j] = wave->data[j];
  }
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, wave->length);
  
  /* Gall gsl spline interpolation */
  gsl_spline_init(spline, x, y, wave->length);
  
  /* Getting first and second order time derivatives from gsl interpolations */
  for (j = 0; j < wave->length; ++j)
  {
	dwave->data[j] = gsl_spline_eval_deriv(spline, j, acc) / dt;
	ddwave->data[j] = gsl_spline_eval_deriv2(spline, j, acc) / dt / dt;
  }
  
  /* Free gsl variables */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  free(x);
  free(y);
  
  return errcode;
}

/* <lalVerbatim file="XLALGenerateQNMFreqCP">  */
INT4 XLALGenerateQNMFreq(
	COMPLEX8Vector			*modefreqs,
	InspiralTemplate		*params,
	UINT4					l,
	UINT4					m,
	UINT4					nmodes
	)
/* </lalVerbatim> */
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  static const char* func = "XLALGenerateQNMFreq";
    
  INT4 i;
  
  /* Get a local copy of the intrinstic parameters */
  REAL8 totalMass = params->totalMass;
  REAL8 eta = params->eta;
  REAL8 finalMass = 0;
  REAL8 finalSpin = 0;
  
  /* Call XLALFinalMassSpin() to get mass and spin of the final black hole */
  errcode = XLALFinalMassSpin(&finalMass, &finalSpin, params);
  if ( errcode != XLAL_SUCCESS )
  {
	  fprintf( stderr, "XLALFinalMassSpin failed\n" );
      exit( 1 );
  }
  
  /* Fitting coefficients for QNM frequencies from PRD73, 064030 */
  REAL4 BCWre[3][3] = { {1.5251, -1.1568, 0.1292},
						{1.3673, -1.0260, 0.1628},
						{1.3223, -1.0257, 0.186} };
  REAL4 BCWim[3][3] = { { 0.7, 1.4187, -0.4990},
						{ 0.1, 0.5436, -0.4731},
						{-0.1, 0.4206, -0.4256} };
						  
  /* QNM frequencies from the fitting given in PRD73, 064030 */
  for (i = 0; i < nmodes; ++i)
  {
	modefreqs->data[i].re = BCWre[i][0] + BCWre[i][1] * pow(1 - finalSpin, BCWre[i][2]);
	modefreqs->data[i].im = modefreqs->data[i].re / 2
							/ (BCWim[i][0] + BCWim[i][1] * pow(1 - finalSpin, BCWim[i][2]));
	modefreqs->data[i].re *= 1 / finalMass / (totalMass * LAL_MTSUN_SI);
	modefreqs->data[i].im *= 1 / finalMass / (totalMass * LAL_MTSUN_SI);
  }
  return errcode;
}

/* <lalVerbatim file="XLALFinalMassSpinCP">  */
INT4 XLALFinalMassSpin(
	REAL8					*finalMass,
	REAL8					*finalSpin,
	InspiralTemplate		*params
	)
/* </lalVerbatim> */
{
  /* XLAL error handling */
  INT4 errcode = XLAL_SUCCESS;
  static const char* func = "XLALFinalMassSpin";
  
  /* get a local copy of the intrinstic parameters */
  REAL8 eta = params->eta;
  REAL8 eta2 = eta * eta;
  
  /* Final mass and spin given by a fitting in PRD76, 104049 */
  *finalMass = 1 - 0.057191 * eta - 0.498 * eta2;
  *finalSpin = 3.464102 * eta - 2.9 * eta2;
  
  return errcode;
}
