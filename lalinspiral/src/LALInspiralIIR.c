/*
  This is a code that takes as its arugment a complex time series, or an
  amplitude and phase time series and produces a set of IIR filters that when
  combined, approximated recreate the original complex time series.
  
  Created by Shaun Hooper 2010-05-28

*/

#include <lal/LALInspiral.h>
#include <math.h>

int XLALCreateInspiralIIRFilters(REAL8Vector *amp, REAL8Vector *phase, double epsilon, double alpha, double beta, COMPLEX16Vector **a1, COMPLEX16Vector **b0, INT4Vector **delay)
{  
	UINT4 j, jmax, jstep, k;
	UINT4 nfilters = 0;
	REAL8 phase_ddot, phase_dot;
	UINT4 i;
	COMPLEX16Vector *resp;

	/* FIXME: Add error checking for lengths of amp and phase */
	jmax = amp->length;

	/* FIXME: Error check that *b0, *a1, *delay ARE null pointers */

	phase_ddot = (phase->data[2] -2.0 * phase->data[1] + phase->data[0]) / ( 2.0 * LAL_PI); // Second derivative of the phase at the first data point
	jstep = (UINT4) floor(sqrt(2.0 * epsilon / phase_ddot) + 0.5);
	j = jstep;

	resp = XLALCreateCOMPLEX16Vector(jmax);
	for (i=0; i<jmax; i++)
		{
			resp->data[j].re = 0.0;
			resp->data[j].im = 0.0;
		}

	while (j < jmax - 2 && jstep != 0) {
		/* FIXME: Check that the flooring of k really is correct */
		nfilters++;
		k = (UINT4 ) floor((REAL8 ) j - alpha * ((REAL8 ) jstep) + 0.5) - 1;
		phase_dot = (-phase->data[k+2] + 8 * (phase->data[k+1] - phase->data[k-1]) + phase->data[k-2]) / 12.0; // Five-point stencil first derivative of phase

		/* FIXME: Should think about being smarter about allocating memory for these (linked list??) */
		*a1 = XLALResizeCOMPLEX16Vector(*a1, nfilters);
		*b0 = XLALResizeCOMPLEX16Vector(*b0, nfilters);
		*delay = XLALResizeINT4Vector(*delay, nfilters);

		/* Record a1, b0 and delay */

		(*a1)->data[nfilters-1] = XLALCOMPLEX16Polar((REAL8) exp(-beta / ((double) jstep)), -phase_dot);		
		(*b0)->data[nfilters-1] = XLALCOMPLEX16Polar(amp->data[k], phase->data[k] + phase_dot * ((REAL8) (j-(k+1))) );
		(*delay)->data[nfilters-1] = jmax - j;
		
		/* I was thinking about putting the resonse function in here */

		/* Calculate the next data point step */
		phase_ddot = (phase->data[j+1] - 2.0 * phase->data[j] + phase->data[j-1]) / (2.0 * LAL_PI);
		jstep = (UINT4) floor(sqrt(2.0 * epsilon / phase_ddot)+0.5);
		j += jstep;

	}
	
	return 0;
}


int XLALCreateIIRResponseSeries(UINT4 N, COMPLEX16Vector *a1, COMPLEX16Vector *b0, INT4Vector *delay, COMPLEX16Vector **response)
{  
	UINT4 f, j, nfilters;
	COMPLEX16 a1f, b0f;
	UINT4 delayf;

	/* FIX ME: Check if a1, b0, delay actually exist */
	nfilters = a1->length;
	*response = XLALResizeCOMPLEX16Vector(*response, N);
	fprintf(stderr, "About to start response vector setup\n");
	for (j = 0; j < N; j++)
		(*response)->data[j] = XLALCOMPLEX16Rect(0.0, 0.0);
	fprintf(stderr, "About to looping filters\n");
	for (f = 0; f < nfilters; f++)
		{
			a1f = a1->data[f];
			b0f = b0->data[f];
			delayf = delay->data[f];

			for (j = delayf; j < N; j++)
				(*response)->data[j] = XLALCOMPLEX16Add((*response)->data[j], XLALCOMPLEX16Mul(b0f, XLALCOMPLEX16PowReal(a1f, (REAL8 ) (j - delayf))));

		}

	return 0;
}

int XLALIIRFreqSeries(UINT4 j, UINT4 jmax, COMPLEX16 a1, COMPLEX16 b0, INT4 delay, COMPLEX16 *hfcos, COMPLEX16 *hfsin)
{
	REAL8 loga1, arga1, pf;
	COMPLEX16 scl, ft, ftconj;

	/* FIXME: Check if a1, b0, delay exist */

	loga1 = XLALCOMPLEX16LogAbs(a1);
	arga1 = XLALCOMPLEX16Arg(a1);
	pf = 2.0 * LAL_PI * ((REAL8 ) j) / ((REAL8 ) jmax);
	scl = XLALCOMPLEX16Polar(0.5, - pf * ((REAL8 ) (jmax - delay)));

	ft = XLALCOMPLEX16Div(b0, XLALCOMPLEX16Rect(-loga1, -arga1 - pf));
	ftconj = XLALCOMPLEX16Div(XLALCOMPLEX16Conjugate(b0), XLALCOMPLEX16Rect(-loga1, arga1 - pf));

	*hfcos = XLALCOMPLEX16Mul(scl, XLALCOMPLEX16Add(ft, ftconj));
	*hfsin = XLALCOMPLEX16Mul(scl, XLALCOMPLEX16Sub(ft, ftconj));

	return 0;
}

int XLALIIRInnerProduct(COMPLEX16Vector *a1, COMPLEX16Vector *b0, INT4Vector *delay, REAL8Vector *psd, REAL8 *ip)
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
					XLALIIRFreqSeries(j, 2 * psd->length, a1->data[k], b0->data[k], delay->data[k], &hfcos, &hfsin);
					hA = XLALCOMPLEX16Add(hA, hfcos);
				}
			*ip += XLALCOMPLEX16Abs2(hA) / (psd->data[j] * ((REAL8 ) psd->length)) * 1.0; 
		}

	return 0;
}
