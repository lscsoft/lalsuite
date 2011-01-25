/*

  This code relates to Infinite Impulse Response filters that correspond to an
  inspiral waveform.  The idea is that a sum a set of delayed first order IIR
  filters with one feedback coefficient (a1) and one feedforward (b0)
  coefficient that will approximate the correlation of the input data and the
  inspiral waveform. 

  I.E the total impulse response is approximately a time reversed inspiral
  waveform.

  To generate the IIR set of a1's, b0's and delays, you need to provide a
  amplitude and phase time series.

  Created by Shaun Hooper 2010-05-28

*/

#include <lal/LALInspiral.h>
#include <math.h>
#include <complex.h>

int XLALInspiralGenerateIIRSet(REAL8Vector *amp, REAL8Vector *phase, double epsilon, double alpha, double beta, COMPLEX16Vector **a1, COMPLEX16Vector **b0, INT4Vector **delay)
{  
	int j, jmax, jstep, k;
	int nfilters = 0;
	double phase_ddot, phase_dot;

	/* FIXME: Add error checking for lengths of amp and phase */
	if (amp->length != phase->length) 
	         XLAL_ERROR(__func__, XLAL_EINVAL);

	jmax = amp->length;

	/* FIXME: Error check that *b0, *a1, *delay ARE null pointers */
	/*if (!*a1 || !*b0 || !*delay)
	         XLAL_ERROR(__func__, XLAL_EFAULT);*/

	phase_ddot = (phase->data[2] -2.0 * phase->data[1] + phase->data[0]) / ( 2.0 * LAL_PI); // Second derivative of the phase at the first data point

	jstep = (int) floor(sqrt(2.0 * epsilon / phase_ddot) + 0.5);
	j = jstep;

	while (j < jmax - 2 && jstep != 0) {
		/* FIXME: Check that the flooring of k really is correct */
		nfilters++;
		k = (int ) floor((double ) j - alpha * ((double ) jstep) + 0.5) - 1;
		phase_dot = (-phase->data[k+2] + 8 * (phase->data[k+1] - phase->data[k-1]) + phase->data[k-2]) / 12.0; // Five-point stencil first derivative of phase

		/* FIXME: Should think about being smarter about allocating memory for these (linked list??) */
		*a1 = XLALResizeCOMPLEX16Vector(*a1, nfilters);
		*b0 = XLALResizeCOMPLEX16Vector(*b0, nfilters);
		*delay = XLALResizeINT4Vector(*delay, nfilters);

		/* Record a1, b0 and delay */
		(*a1)->data[nfilters-1] = XLALCOMPLEX16Polar((double) exp(-beta / ((double) jstep)), -phase_dot);
		(*b0)->data[nfilters-1] = XLALCOMPLEX16Polar(amp->data[k], phase->data[k] + phase_dot * ((double) (j-(k+1))) );
		(*delay)->data[nfilters-1] = jmax - j;

		/* Calculate the next data point step */
		phase_ddot = (phase->data[j+1] - 2.0 * phase->data[j] + phase->data[j-1]) / (2.0 * LAL_PI);
		jstep = (int) floor(sqrt(2.0 * epsilon / phase_ddot)+0.5);
		j += jstep;

	}

	return 0;
}


int XLALInspiralIIRSetResponse(int N, COMPLEX16Vector *a1, COMPLEX16Vector *b0, INT4Vector *delay, COMPLEX16Vector **response)
{  
	int f, j, nfilters;
	COMPLEX16 a1f, b0f;
	int delayf;

	/* FIX ME: Check if a1, b0, delay actually exist */
	nfilters = a1->length;
	*response = XLALResizeCOMPLEX16Vector(*response, N);

	for (j = 0; j < N; j++)
		(*response)->data[j] = XLALCOMPLEX16Rect(0.0, 0.0);

	for (f = 0; f < nfilters; f++)
		{
			a1f = a1->data[f];
			b0f = b0->data[f];
			delayf = delay->data[f];

			for (j = delayf; j < N; j++)
				(*response)->data[j] = XLALCOMPLEX16Add((*response)->data[j], XLALCOMPLEX16Mul(b0f, XLALCOMPLEX16PowReal(a1f, (double ) (j - delayf))));

		}

	return 0;
}

int XLALInspiralGenerateIIRSetFourierTransform(int j, int jmax, COMPLEX16 a1, COMPLEX16 b0, INT4 delay, COMPLEX16 *hfcos, COMPLEX16 *hfsin)
{
	double loga1, arga1, pf;
	COMPLEX16 scl, ft, ftconj;

	/* FIXME: Check if a1, b0, delay exist */

	loga1 = XLALCOMPLEX16LogAbs(a1);
	arga1 = XLALCOMPLEX16Arg(a1);
	pf = 2.0 * LAL_PI * ((double ) j) / ((double ) jmax);
	scl = XLALCOMPLEX16Polar(0.5, - pf * ((double ) (jmax - delay)));

	ft = XLALCOMPLEX16Div(b0, XLALCOMPLEX16Rect(-loga1, -arga1 - pf));
	ftconj = XLALCOMPLEX16Div(XLALCOMPLEX16Conjugate(b0), XLALCOMPLEX16Rect(-loga1, arga1 - pf));

	*hfcos = XLALCOMPLEX16Mul(scl, XLALCOMPLEX16Add(ft, ftconj));
	*hfsin = XLALCOMPLEX16Mul(scl, XLALCOMPLEX16Sub(ft, ftconj));

	return 0;
}

int XLALInspiralCalculateIIRSetInnerProduct(COMPLEX16Vector *a1, COMPLEX16Vector *b0, INT4Vector *delay, REAL8Vector *psd, double *ip)
{
	UINT4 k, j;
	COMPLEX16 hA;
	COMPLEX16 hfcos = XLALCOMPLEX16Rect(0.0, 0.0);
	COMPLEX16 hfsin = XLALCOMPLEX16Rect(0.0, 0.0);

	*ip = 0.0; 

	for (j = 0; j < psd->length; j++)
		{
			hA = XLALCOMPLEX16Rect(0.0, 0.0);
			for (k = 0; k < a1->length; k++)
				{
					XLALInspiralGenerateIIRSetFourierTransform(j, 2 * psd->length, a1->data[k], b0->data[k], delay->data[k], &hfcos, &hfsin);
					hA = XLALCOMPLEX16Add(hA, hfcos);
				}
			*ip += XLALCOMPLEX16Abs2(hA) / (psd->data[j] * ((double ) psd->length)) * 1.0; 
		}

	return 0;
}
