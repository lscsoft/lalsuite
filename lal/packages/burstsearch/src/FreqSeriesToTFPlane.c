/******** <lalVerbatim file="FreqSeriesToTFPlaneCV"> ********
Author: Flanagan, E
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID(FREQSERIESTOTFPLANEC, "$Id$");

#include <math.h>
#include <lal/Date.h>
#include <lal/LALConstants.h>
#include <lal/LALErrno.h>
#include <lal/RealFFT.h>
#include <lal/SeqFactories.h>
#include <lal/TFTransform.h>

/******** <lalVerbatim file="FreqSeriesToTFPlaneCP"> ********/
int XLALFreqSeriesToTFPlane(
	COMPLEX8TimeFrequencyPlane *tfp,
	const COMPLEX8FrequencySeries *freqSeries,
	UINT4 windowShift,
	REAL4 *normalisation,
	const REAL4FrequencySeries *psd
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALFreqSeriesToTFPlane";
	REAL4Vector *filter = NULL;
	REAL4Vector *snr = NULL;
	COMPLEX8Vector *tmp = NULL;
	COMPLEX8Vector *fcorr = NULL;
	INT4 i;
	INT4 j;
	INT4 numpts = 0;
	INT4 nt = tfp->params->timeBins;
	INT4 nf = tfp->params->freqBins;
	INT4 ntotal;

	INT4 flow1;
	INT4 fseglength;
	INT4 fcut;
	INT4 fwindow;

	REAL4 delF = 0;
	REAL4 delT = 0;
	REAL4 dt = 0;
	REAL4 twopiOverNumpts = 0;
	INT4 filterlen = 0;

	RealFFTPlan *pfwd = NULL;
	RealFFTPlan *prev = NULL;

	/* check input parameters */
	if ((nt <= 0) || (nf <= 0) || (freqSeries->deltaF <= 0.0) || (tfp->params->deltaT <= 0.0) || (freqSeries->f0 < 0.0))
		XLAL_ERROR(func, XLAL_EDOM);

	/* Lowest freq of time freq plane >= lowest freq of freq series */
	if (tfp->params->flow < freqSeries->f0)
		XLAL_ERROR(func, XLAL_EDATA);

	/* 
	 * delF is the frequency resoltion of the time-frequency plane.  
	 * delT is the time resolution of the plane.
	 * Note that the frequency resolution of the frequency series does
	 * not have to be the same as delF.
	 */
	delT = tfp->params->deltaT;
	delF = tfp->params->deltaF;

	/* number of bins of the frequency series corresponding to delF */
	fseglength = 0.5 + delF / freqSeries->deltaF;

	/* low frequency cutoff:  in terms of number of bins of frequency series */
	flow1 = (tfp->params->flow - freqSeries->f0) / freqSeries->deltaF;

	/* compute total number of data points to be used to construct TF plane */
	ntotal = nf * fseglength;

	/* make sure have enough data points in freq series */
	if (ntotal + flow1 > (INT4) freqSeries->data->length)
		XLAL_ERROR(func, XLAL_EDATA);

	/* create temporary vectors for filter and correlations */
	numpts = 2 * (freqSeries->data->length - 1);
	twopiOverNumpts = 2.0 * LAL_PI / (float) numpts;
	filter = XLALCreateREAL4Vector(numpts);
	tmp = XLALCreateCOMPLEX8Vector(freqSeries->data->length);
	fcorr = XLALCreateCOMPLEX8Vector(freqSeries->data->length);
	snr = XLALCreateREAL4Vector(numpts);
	if (!filter || !tmp || !fcorr || !snr)
		XLAL_ERROR(func, XLAL_EFUNC);

	/* sampling rate of time series which gave freqSeries */
	dt = 1.0 / (((REAL4) numpts) * freqSeries->deltaF);

	/* set the epoch of the TF plane */
	tfp->epoch = freqSeries->epoch;

	/* number of points from peak of filter to first zero */
	filterlen = numpts / fseglength;

	/* Create FFTW plans for forward and reverse REAL4 FFTs */
	pfwd = XLALCreateForwardREAL4FFTPlan(numpts, 0);
	prev = XLALCreateReverseREAL4FFTPlan(numpts, 0);
	if (!pfwd || !prev)
		XLAL_ERROR(func, XLAL_EFUNC);

	/* 
	 * PRB - test code to add a delta function to the segment and
	 * confirm the output of the code
	 */
	/*
	   memset(filter->data, 0, numpts * sizeof(REAL4));
	   filter->data[numpts/2] = 1.0;
	   LALForwardRealFFT( status->statusPtr, freqSeries->data, filter, pfwd);
	   CHECKSTATUSPTR (status);
	 */

	/* PRB - generate the time domain filter function.  This filter should
	 * depends on the frequency flow1 and the fseglength.  It should
	 * be non-zero for some amount that depends on fseglength and
	 * the total bandwidth of the input frequency series.  
	 */
	memset(filter->data, 0, numpts * sizeof(REAL4));
	filter->data[0] = twopiOverNumpts * fseglength / (LAL_PI * dt);
	for (j = 1; j < filterlen; j++)
		filter->data[j] = filter->data[numpts - j] = (sin(twopiOverNumpts * j * (flow1 + fseglength)) - sin(twopiOverNumpts * j * flow1)) / (LAL_PI * 1.0 * ((float) j) * dt);

	/* PRB - Fourier transform the filter into the frequency domain */
	if (XLALREAL4ForwardFFT(tmp, filter, pfwd))
		XLAL_ERROR(func, XLAL_EFUNC);

	/* set the frequency window(bins); for each channel contributions
	 * will be zeroed out from fwindow before and after
	 */
	fwindow = 100;
	fcut = flow1 - fwindow;

	/* loop over different basis vectors in the frequency domain */
	for (i = 0; i < nf; i++) {
		normalisation[i] = 0.0;

		/* All values below (i * df - fwindow) to zero */
		for (j = 0; j < (fcut + (i * fseglength)); j++) {
			REAL4 reFilter = 0.0;
			REAL4 imFilter = 0.0;
			REAL4 reData = freqSeries->data->data[j].re;
			REAL4 imData = freqSeries->data->data[j].im;

			fcorr->data[j].re = reFilter * reData + imFilter * imData;
			fcorr->data[j].im = reFilter * imData - imFilter * reData;

			normalisation[i] += (reFilter * reFilter + imFilter * imFilter);
		}

		/* PRB - Multiply the filter by the data.  Don't forget complex
		 * conjugate and any other relevant information */
		for (j = (fcut + (i * fseglength)); j < (fcut + (2 * fwindow) + ((i + 1) * fseglength)); j++) {
			REAL4 reFilter = tmp->data[j - i * fseglength].re;
			REAL4 imFilter = tmp->data[j - i * fseglength].im;
			REAL4 reData = freqSeries->data->data[j].re;
			REAL4 imData = freqSeries->data->data[j].im;

			if(psd) {
				reFilter /= sqrt(psd->data->data[j]);
				imFilter /= sqrt(psd->data->data[j]);
			}

			fcorr->data[j].re = reFilter * reData + imFilter * imData;
			fcorr->data[j].im = reFilter * imData - imFilter * reData;

			normalisation[i] += reFilter * reFilter + imFilter * imFilter;
		}

		for (j = (fcut + (2 * fwindow) + ((i + 1) * fseglength)); (unsigned) j < freqSeries->data->length; j++) {
			REAL4 reFilter = 0.0;
			REAL4 imFilter = 0.0;
			REAL4 reData = freqSeries->data->data[j].re;
			REAL4 imData = freqSeries->data->data[j].im;

			fcorr->data[j].re = reFilter * reData + imFilter * imData;
			fcorr->data[j].im = reFilter * imData - imFilter * reData;

			normalisation[i] += (reFilter * reFilter + imFilter * imFilter);
		}


		normalisation[i] = sqrt(4.0 * normalisation[i]);

		/* clean up Nyquist */
		fcorr->data[freqSeries->data->length - 1].re = 0.0;
		fcorr->data[freqSeries->data->length - 1].im = 0.0;

		/* PRB - Inverse transform the product so that we get a time series at
		 * the full sample rate.   Make sure to check that the sample
		 * rate agrees with what you had before. */
		if (XLALREAL4ReverseFFT(snr, fcorr, prev))
			XLAL_ERROR(func, XLAL_EFUNC);

		/* 
		 * PRB & SRM - to output the time series when the delta function
		 * is inserted instead of other data 
		 */
		/*
		   {
		   char fname[248];
		   sprintf(fname, "z%02d.dat", i);
		   fp = fopen(fname,"a");
		   for(j=0; j<numpts ; j++)
		   {
		   fprintf(fp,"%i %f\n",j,snr->data[j]/norm);
		   }
		   fclose(fp);
		   }
		 */

		/* PRB - copy the data back into the time series.  In this
		 * process,  one only takes the samples corresponding to the
		 * uncorupted times.   This means that the data should be longer
		 * than the 1/dfmin.   This can mean a change compared to what
		 * has been done before.   We'll have to check that carefully.
		 * Notice that the originally complex plane is now real.  I still
		 * don't understand why the difference arises with the Eanns's
		 * original implementation.   We'll have to look at that.  */
		/* Copy the result into appropriate spot in output structure */
		for (j = 0; j < nt; j++) {
			tfp->data[j * nf + i].re = snr->data[(j * (INT4) (delT / dt)) + windowShift];
			tfp->data[j * nf + i].im = 0.0;
		}
	}

	/* Get rid of all temporary memory */
	XLALDestroyREAL4FFTPlan(pfwd);
	XLALDestroyREAL4FFTPlan(prev);
	XLALDestroyREAL4Vector(filter);
	XLALDestroyCOMPLEX8Vector(tmp);
	XLALDestroyCOMPLEX8Vector(fcorr);
	XLALDestroyREAL4Vector(snr);

	/* normal exit */
	return (0);
}
