/*----------------------------------------------------------------------- 
 * 
 * File Name:  CrossCorr.c
 * 
 * Author: Steve Drasco
 * 
 * Revision: $Id$ 
 * 
 *----------------------------------------------------------------------- 
 * 
 * NAME 
 * CrossCorr.c
 * 
 * SYNOPSIS 
 * void LALCrossCorr (LALStatus *status, REAL4 *out, CCIn *in);
 *
 * typedef struct tagCCIn {
 *      REAL4TimeSeries	*g1;
 *      REAL4TimeSeries	*g2;
 *      REAL4Vector	*QmaxTilde;
 *      INT2		plan;
 * } CCIn;
 * 
 *
 * DESCRIPTION 
 * Calculates the cross-correlation statistic.
 *
 * DIAGNOSTICS 
 * 
 * CALLS
 *  InitStatus()
 *  LALSCreateVector()
 *  LALSDestroyVector()
 *  LALCCreateVector()
 *  LALCDestroyVector()
 *  LALMeasureFwdRealFFTPlan()
 *  EstimateFrdRealFFTPlan()
 *  LALFwdRealFFT()
 *  LALDestroyRealFFTPlan()
 * 
 * NOTES
 *
 *-----------------------------------------------------------------------
 */


#include <math.h>
#include "LALStdlib.h"
#include "AVFactories.h"
#include "RealFFT.h"
#include "CrossCorr.h"

NRCSID (CROSSCORRC, "$Id$");

extern INT4 LALDebugLevel;

void 
LALCrossCorr ( LALStatus *status,
            REAL4  *out,
            CCIn   *in     )
{
	INT4		i, N;
	REAL4Vector	*g1Bar=NULL;
	REAL4Vector	*g2Bar=NULL;
	RealFFTPlan	*g1Plan=NULL;
	RealFFTPlan	*g2Plan=NULL;
        COMPLEX8Vector	*g1BarTilde=NULL;
	COMPLEX8Vector	*g2BarTilde=NULL;

	/* initialize status structure */
	INITSTATUS(status, "LALCrossCorr", CROSSCORRC);
	ATTATCHSTATUSPTR (status);

	/* check address of input structure */
	ASSERT(in != NULL, status, CROSSCORR_EIN, CROSSCORR_MSGEIN);

        /* check address of output */
        ASSERT(out != NULL, status, CROSSCORR_EOUT, CROSSCORR_MSGEOUT);

	/* check addresses of detector outputs */
	ASSERT(in->g1->data != NULL, status, CROSSCORR_ENULL, CROSSCORR_MSGENULL);
	ASSERT(in->g2->data != NULL, status, CROSSCORR_ENULL, CROSSCORR_MSGENULL);

	/* store length of time series data */
	N = in->g1->data->length;

	/* compare lengths of time series vectors */
	ASSERT(in->g2->data->length == N, status, CROSSCORR_ESIZE1, CROSSCORR_MSGESIZE1);

	/* check length of kernel vector */
	ASSERT(in->QmaxTilde->length == N, status, CROSSCORR_ESIZE2, CROSSCORR_MSGESIZE2);

        /* allocate padded vectors and FFT vectors */
	LALSCreateVector(status->statusPtr,&g1Bar,2*N-1);
	LALSCreateVector(status->statusPtr,&g2Bar,2*N-1);
	LALCCreateVector(status->statusPtr,&g1BarTilde,N);
	LALCCreateVector(status->statusPtr,&g2BarTilde,N);

	/* fill */
        for (i=0; i<N; i++) {
		g1Bar->data[i] = in->g1->data->data[i];
		g2Bar->data[i] = in->g2->data->data[i];
        }

        /* pad */
	for (i=N; i<2*N-1; i++) {
		g1Bar->data[i] = g2Bar->data[i] = 0.0;
	}

	/* FFT plans */
	if (in->plan != 1) {
		LALMeasureFwdRealFFTPlan(status->statusPtr,&g1Plan,2*N-1);
		LALMeasureFwdRealFFTPlan(status->statusPtr,&g2Plan,2*N-1);
	} else {
                LALEstimateFwdRealFFTPlan(status->statusPtr,&g1Plan,2*N-1);
                LALEstimateFwdRealFFTPlan(status->statusPtr,&g2Plan,2*N-1);
	}

	/* FFT */
	LALFwdRealFFT(status->statusPtr,g1BarTilde,g1Bar,g1Plan);
	LALFwdRealFFT(status->statusPtr,g2BarTilde,g2Bar,g2Plan);

	/* initialize output with zero index part */
	*out = in->QmaxTilde->data[0] * 0.5 * ( g1BarTilde->data[0].re * g2BarTilde->data[0].re
                                              + g1BarTilde->data[0].im * g2BarTilde->data[0].im  );

        /* sum over indices bigger than zero */
        for (i=1; i < N; i++){  
		*out += in->QmaxTilde->data[i] * ( g1BarTilde->data[i].re * g2BarTilde->data[i].re 
                                                 + g1BarTilde->data[i].im * g2BarTilde->data[i].im  );
	} 

        /* normalize */
        *out /= (REAL4) N-0.5;

	/* free the memory */
	LALSDestroyVector(status->statusPtr,&g1Bar);
	LALSDestroyVector(status->statusPtr,&g2Bar);
	LALCDestroyVector(status->statusPtr,&g1BarTilde);
	LALCDestroyVector(status->statusPtr,&g2BarTilde);
	LALDestroyRealFFTPlan(status->statusPtr,&g1Plan);
	LALDestroyRealFFTPlan(status->statusPtr,&g2Plan);
	

	/* normal exit */
	DETATCHSTATUSPTR (status);
	RETURN (status);

} /* LALCrossCorr() */
