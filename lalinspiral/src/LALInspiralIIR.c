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

static REAL8 clogabs(COMPLEX16 z)
{
  REAL8 xabs = fabs(creal(z));
  REAL8 yabs = fabs(cimag(z));
  REAL8 max, u;
  if (xabs >= yabs)
  {
    max = xabs;
    u = yabs / xabs;
  }
  else
  {
    max = yabs;
    u = xabs / yabs;
  }
  return log(max) + 0.5 * log1p(u * u);
}

int XLALInspiralGenerateIIRSet(REAL8Vector *amp, REAL8Vector *phase, double epsilon, double alpha, double beta, double padding, COMPLEX16Vector **a1, COMPLEX16Vector **b0, INT4Vector **delay)
{  
	int j = amp->length-1, jstep, k;
	int nfilters = 0, decimationFactor = 1;
	double phase_ddot, phase_dot;

	/* FIXME: Add error checking for lengths of amp and phase */
	if (amp->length != phase->length) 
	         XLAL_ERROR(XLAL_EINVAL);

	*a1 = XLALCreateCOMPLEX16Vector(0);
	*b0 = XLALCreateCOMPLEX16Vector(0);
	*delay = XLALCreateINT4Vector(0);

	while (j >= 2 ) {
		//int prej = j;
		/* Reset j so that the delay will be an integar number of decimated rate */
		//j = amp->length-1 - (int) floor((amp->length-1-j)/(double) decimationFactor + 0.5)*decimationFactor;

		/* Get error term */
		phase_ddot = (phase->data[j-2] - 2.0 * phase->data[j-1] + phase->data[j]) / (2.0 * LAL_PI);
		jstep = (int) floor(sqrt(2.0 * epsilon / phase_ddot)+0.5);
		k = (int ) floor((double ) j - alpha * (double ) jstep + 0.5);

		if (k <= 2) break;
		nfilters++;

		if (k > (int) amp->length-3) {
			phase_dot = (11.0/6.0*phase->data[k] - 3.0*phase->data[k-1] + 1.5*phase->data[k-2] - 1.0/3.0*phase->data[k-3]);
		}
		else {
			phase_dot = (-phase->data[k+2] + 8 * (phase->data[k+1] - phase->data[k-1]) + phase->data[k-2]) / 12.0; // Five-point stencil first derivative of phase
		}
		//fprintf(stderr, "%3.0d, %6.0d, %3.0d, %11.2f, %11.8f\n",nfilters, amp->length-1-j, decimationFactor, ((double) (amp->length-1-j))/((double) decimationFactor), phase_dot/(2.0*LAL_PI)*2048.0);
		decimationFactor = ((int ) pow(2.0,-ceil(log(2.0*padding*phase_dot/(2.0*LAL_PI))/log(2.0))));
		if (decimationFactor < 1 ) decimationFactor = 1;

		//fprintf(stderr, "filter = %d, prej = %d, j = %d, k=%d, jstep = %d, decimation rate = %d, nFreq = %e, phase[k] = %e\n", nfilters, prej, j, k, jstep, decimationFactor, phase_dot/(2.0*LAL_PI), phase->data[k]);
		/* FIXME: Should think about being smarter about allocating memory for these (linked list??) */
		*a1 = XLALResizeCOMPLEX16Vector(*a1, nfilters);
		*b0 = XLALResizeCOMPLEX16Vector(*b0, nfilters);
		*delay = XLALResizeINT4Vector(*delay, nfilters);

		/* Record a1, b0 and delay */
		(*a1)->data[nfilters-1] = cpolar((double) exp(-beta / ((double) jstep)), -phase_dot);
		(*b0)->data[nfilters-1] = cpolar(amp->data[k], phase->data[k] + phase_dot * ((double) (j - k)) );
		(*delay)->data[nfilters-1] = amp->length - 1 - j;



		/* Calculate the next data point step */
		j -= jstep;

	}

	return 0;
}


int XLALInspiralIIRSetResponse(COMPLEX16Vector *a1, COMPLEX16Vector *b0, INT4Vector *delay, COMPLEX16Vector *response)
{
	int length, numFilters;
	complex double y;
	complex double *a1_last;
	complex double *a1f = (complex double *) a1->data;
	complex double *b0f = (complex double *) b0->data;
	int *delayf = delay->data;

	if(a1->length != b0->length || a1->length != delay->length)
		XLAL_ERROR(XLAL_EBADLEN);

	memset(response->data, 0, sizeof(complex double) * response->length);

	numFilters = a1->length;
	for (a1_last = a1f + numFilters; a1f < a1_last; a1f++)
		{
			y = *b0f / *a1f;
			length = (int) (logf((1e-13)/cabs(*b0f)))/(logf(cabs(*a1f)));// + *delayf;
			int maxlength = response->length - *delayf;
			if (length > maxlength)
				length = maxlength;

			complex double *response_data = (complex double *) &response->data[*delayf];
			complex double *response_data_last;

			for (response_data_last = response_data + length; response_data < response_data_last; response_data++ )
				{
					y *= *a1f;
					*response_data += y;
				}
			b0f++;
			delayf++;
		}
	return 0;
}

int XLALInspiralGenerateIIRSetFourierTransform(int j, int jmax, COMPLEX16 a1, COMPLEX16 b0, INT4 delay, COMPLEX16 *hfcos, COMPLEX16 *hfsin)
{
	double loga1, arga1, pf;
	COMPLEX16 scl, ft, ftconj;

	/* FIXME: Check if a1, b0, delay exist */

	loga1 = clogabs(a1);
	arga1 = carg(a1);
	pf = 2.0 * LAL_PI * ((double ) j) / ((double ) jmax);
	scl = cpolar(0.5, - pf * ((double ) (jmax - delay)));

	ft = b0 / crect(-loga1, -arga1 - pf);
	ftconj = conj(b0) / crect(-loga1, arga1 - pf);

	*hfcos = scl * (ft + ftconj);
	*hfsin = scl * (ft - ftconj);

	return 0;
}

int XLALInspiralCalculateIIRSetInnerProduct(COMPLEX16Vector *a1, COMPLEX16Vector *b0, INT4Vector *delay, REAL8Vector *psd, double *ip)
{
	UINT4 k, j;
	COMPLEX16 hA;
	COMPLEX16 hfcos = 0.0;
	COMPLEX16 hfsin = 0.0;

	*ip = 0.0; 

	for (j = 0; j < psd->length; j++)
		{
			hA = 0.0;
			for (k = 0; k < a1->length; k++)
				{
					XLALInspiralGenerateIIRSetFourierTransform(j, 2 * psd->length, a1->data[k], b0->data[k], delay->data[k], &hfcos, &hfsin);
					hA = hA + hfcos;
				}
			*ip += cabs(hA) * cabs(hA) / (psd->data[j] * ((double ) psd->length)) * 1.0;
		}

	return 0;
}
