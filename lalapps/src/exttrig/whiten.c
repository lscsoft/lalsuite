#include <math.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALComplex.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/TimeFreqFFT.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeSeries.h>
#include <lal/Window.h>
#include "whiten.h"

REAL8Sequence *XLALMoveREAL8Sequence( REAL8Sequence *destination, const REAL8Sequence *source, int first )
{
	static const char *func = "XLALMoveREAL8Sequence";
	if ( ! destination || ! source )
		XLAL_ERROR_NULL( func, XLAL_EFAULT );
	if ( ! destination->data || ! source->data )
		XLAL_ERROR_NULL( func, XLAL_EINVAL );
	if ( first < -(int)destination->length || first > (int)source->length )
		memset( destination->data, 0, destination->length * sizeof(*destination->data) );
	else {
		REAL8 *src = source->data;
		UINT4  len = source->length;
		REAL8 *dst = destination->data;
		UINT4  cpy = destination->length;
		if ( first < 0 ) {
			memset( dst, 0, (-first) * sizeof(*dst) );
			cpy -= -first;
			dst += -first;
			first = 0;
		}
		src += first;
		len -= first;
		if ( cpy > len ) {
			memset( dst + len, 0, (cpy - len) * sizeof(*dst) );
			cpy = len;
		}
		memcpy( dst, src, cpy * sizeof(*dst) );
	}
	return destination;
}


REAL8TimeSeries *XLALMoveREAL8TimeSeries( REAL8TimeSeries *destination, const REAL8TimeSeries *source, int first )
{
	static const char *func = "XLALMoveREAL8TimeSeries";
	REAL8Sequence *sequence;
	if ( ! destination || ! source )
		XLAL_ERROR_NULL( func, XLAL_EFAULT );
	sequence = destination->data;
	if ( NULL == XLALMoveREAL8Sequence( sequence, source->data, first ) )
		XLAL_ERROR_NULL( func, XLAL_EFUNC );
	*destination = *source;
	destination->data = sequence;
	XLALGPSAdd( &destination->epoch, first * destination->deltaT );
	return destination;
}


/* IN PLACE whitening */
static int datacond1( REAL8TimeSeries *series, REAL8Window *window, REAL8FrequencySeries *psd1, REAL8FrequencySeries *psd2, COMPLEX16FrequencySeries *condresp, COMPLEX16FrequencySeries *work, REAL8FFTPlan *fwdplan, REAL8FFTPlan *revplan, REAL8 fmin, REAL8 flow, REAL8 fhigh, REAL8 fmax )
{

	UINT4 seglen = series->data->length;
	REAL8 duration = series->deltaT * seglen;
	UINT4 kmin  = fmin  * duration;
	UINT4 klow  = flow  * duration;
	UINT4 khigh = fhigh * duration;
	UINT4 kmax  = fmax  * duration;
	UINT4 k;

	if ( kmax > seglen/2 )
		kmax = seglen/2;
	if ( khigh > seglen/2 )
		khigh = seglen/2;
	if ( klow > seglen/2 || kmin > seglen/2 )
		XLAL_ERROR( "datacond1", XLAL_EFAILED ) /* error */ ;

	if ( window ) {
		UINT4 j;
		for ( j = 0; j < seglen; ++j )
			series->data->data[j] *= window->data->data[j];
	}

	XLALREAL8TimeFreqFFT( work, series, fwdplan );

	for ( k = 0; k < kmin; ++k )
		LAL_SET_REAL( &work->data->data[k], 0.0 );
	if ( psd1 && psd2 )
		for ( ; k < kmax; ++k )
			work->data->data[k] = XLALCOMPLEX16DivReal( work->data->data[k], sqrt(0.5*(psd1->data->data[k]+psd2->data->data[k])) );
	else if ( psd1 )
		for ( ; k < kmax; ++k )
			work->data->data[k] = XLALCOMPLEX16DivReal( work->data->data[k], sqrt(psd1->data->data[k]) );
	else if ( psd2 )
		for ( ; k < kmax; ++k )
			work->data->data[k] = XLALCOMPLEX16DivReal( work->data->data[k], sqrt(psd2->data->data[k]) );
	/* else do nothing: just band-pass filter */
	for ( k = kmax; k <= seglen/2; ++k )
		LAL_SET_REAL( &work->data->data[k], 0.0 );

#if 0
	/* frequency domain band-pass filtering */
	/* window from min to low and from max to high */
	for ( k = kmin; k < klow; ++k )
		work->data->data[k] = XLALCOMPLEX16MulReal( work->data->data[k], cos((LAL_PI*(klow - k))/(2.0*(klow - kmin))) );
	for ( k = khigh; k < kmax; ++k )
		work->data->data[k] = XLALCOMPLEX16MulReal( work->data->data[k], cos((LAL_PI*(k - khigh))/(2.0*(kmax - khigh))) );
	
	/* phase correction */
	if ( condresp )
		for ( k = kmin; k < kmax; ++k )
			work->data->data[k] = XLALCOMPLEX16Mul(work->data->data[k],XLALCOMPLEX16Exp(XLALCOMPLEX16MulReal(LAL_COMPLEX16_I,XLALCOMPLEX16Arg(response->data->data[k]))));
#endif
	for ( k = 0; k < work->data->length; ++k )
		work->data->data[k] = XLALCOMPLEX16Mul( work->data->data[k], condresp->data->data[k] );


	XLALREAL8FreqTimeFFT( series, work, revplan );
	return 0;
}


static int datacond2( REAL8TimeSeries *whitened, REAL8TimeSeries *unwhitened, REAL8Window *window, UINT4 stride, COMPLEX16FrequencySeries *response, REAL8FFTPlan *fwdplan, REAL8FFTPlan *revplan, REAL8 fmin, REAL8 flow, REAL8 fhigh, REAL8 fmax )
{
	COMPLEX16FrequencySeries *thisfft;
	REAL8FrequencySeries     *prevpsd;
	REAL8FrequencySeries     *thispsd;
	REAL8FrequencySeries     *nextpsd;
	REAL8TimeSeries          *thisseg;
	REAL8TimeSeries          *nextseg;

	LIGOTimeGPS epoch = unwhitened->epoch;
	UINT4 reclen = unwhitened->data->length;
	UINT4 seglen = window->data->length;
	UINT4 numseg = 1 + (reclen - seglen)/stride;
	REAL8 deltaT = unwhitened->deltaT;
	REAL8 deltaF = 1.0 / (seglen * deltaT);

	UINT4 seg;
	UINT4 j;

	memset( whitened->data->data, 0, whitened->data->length * sizeof(*whitened->data->data) );

	thisfft = XLALCreateCOMPLEX16FrequencySeries( "FFT", &epoch, 0.0, deltaF, &lalDimensionlessUnit, seglen/2 + 1 );

	prevpsd = XLALCreateREAL8FrequencySeries( "PSD A", &epoch, 0.0, deltaF, &lalDimensionlessUnit, seglen/2 + 1 );
	thispsd = XLALCreateREAL8FrequencySeries( "PSD B", &epoch, 0.0, deltaF, &lalDimensionlessUnit, seglen/2 + 1 );
	nextpsd = XLALCreateREAL8FrequencySeries( "PSD C", &epoch, 0.0, deltaF, &lalDimensionlessUnit, seglen/2 + 1 );

	thisseg = XLALCutREAL8TimeSeries( unwhitened, 0, seglen );
	nextseg = XLALCutREAL8TimeSeries( unwhitened, stride, seglen );

	XLALREAL8ModifiedPeriodogram( thispsd, thisseg, window, fwdplan );

	for ( seg = 0; seg < numseg; ++seg ) {
		REAL8FrequencySeries *temppsd;
		REAL8TimeSeries      *tempseg;
		if ( seg == numseg - 1 ) {
			datacond1( thisseg, window, prevpsd, NULL, response, thisfft, fwdplan, revplan, fmin, flow, fhigh, fmax );
		}
		else {
			XLALMoveREAL8TimeSeries( nextseg, unwhitened, (seg + 1) * stride );
			XLALREAL8ModifiedPeriodogram( nextpsd, nextseg, window, fwdplan );
			datacond1( thisseg, window, seg ? prevpsd : NULL, nextpsd, response, thisfft, fwdplan, revplan, fmin, flow, fhigh, fmax );
		}

		for ( j = (seglen - stride)/2; j < (seglen + stride)/2; ++j )
			whitened->data->data[j + seg*stride] = thisseg->data->data[j];

		/* roll the psd pointers */
		temppsd = prevpsd;
		prevpsd = thispsd;
		thispsd = nextpsd;
		nextpsd = temppsd;
		
		/* swap the segment pointers */
		tempseg = thisseg;
		thisseg = nextseg;
		nextseg = tempseg;
	}

	XLALDestroyREAL8TimeSeries( nextseg );
	XLALDestroyREAL8TimeSeries( thisseg );
	XLALDestroyREAL8FrequencySeries( nextpsd );
	XLALDestroyREAL8FrequencySeries( thispsd );
	XLALDestroyREAL8FrequencySeries( prevpsd );
	XLALDestroyCOMPLEX16FrequencySeries( thisfft );

	return 0;
}


COMPLEX16FrequencySeries * condition_response( COMPLEX8FrequencySeries *response, REAL8 flow, REAL8 fhigh, REAL8 fmin, REAL8 fmax, UINT4 seglen, REAL8 deltaT )
{
	COMPLEX16FrequencySeries *condresp;
	LIGOTimeGPS epoch = { 0 , 0 };  /* irrelevant */
	REAL8 duration = seglen * deltaT;
	UINT4 kmin  = fmin  * duration;
	UINT4 klow  = flow  * duration;
	UINT4 khigh = fhigh * duration;
	UINT4 kmax  = fmax  * duration;
	UINT4 k;

	condresp = XLALCreateCOMPLEX16FrequencySeries( "CONDRESP", &epoch, 0.0, 1.0/duration, &lalDimensionlessUnit, seglen/2 + 1 );
	if ( kmax > condresp->data->length )
		kmax = condresp->data->length;
	if ( khigh > condresp->data->length )
		khigh = condresp->data->length;

	/* frequency domain band-pass filtering */
	/* window from min to low and from max to high */
	for ( k = 0; k < kmin; ++k )
		condresp->data->data[k] = LAL_COMPLEX16_ZERO;
	for ( ; k < klow; ++k )
		condresp->data->data[k] = XLALCOMPLEX16MulReal( LAL_COMPLEX16_ONE, cos((LAL_PI*(klow - k))/(2.0*(klow - kmin))) );
	for ( ; k < khigh; ++k )
		condresp->data->data[k] = LAL_COMPLEX16_ONE;
	for ( ; k < kmax; ++k )
		condresp->data->data[k] = XLALCOMPLEX16MulReal( LAL_COMPLEX16_ONE, cos((LAL_PI*(k - khigh))/(2.0*(kmax - khigh))) );
	for ( ; k < condresp->data->length; ++k )
		condresp->data->data[k] = LAL_COMPLEX16_ZERO;

	/* if response present, extract phases */
	if ( response ) {
		/* TODO: check lengths */
		for ( k = kmin; k < kmax; ++k ) {
			COMPLEX16 phasefac;
			COMPLEX16 phase;
			COMPLEX16 r16;
			COMPLEX8  r8;
			r8 = response->data->data[k];
		       	LAL_SET_COMPLEX( &r16, LAL_REAL(r8), LAL_IMAG(r8) );

			phase = XLALCOMPLEX16MulReal(LAL_COMPLEX16_I,XLALCOMPLEX16Arg(r16));
			phasefac = XLALCOMPLEX16Exp(phase);

			condresp->data->data[k] = XLALCOMPLEX16Mul(condresp->data->data[k],phasefac);
		}
	}
	return condresp;
}


REAL8TimeSeries * XLALLeonorWhitenREAL8TimeSeries( REAL8TimeSeries *unwhitened, REAL8FFTPlan *fwdplan, REAL8FFTPlan *revplan, UINT4 seglen, COMPLEX8FrequencySeries *response )
{
	const REAL8 beta  = 0.75;
	const REAL8 flow  = 40.0;
	const REAL8 fhigh = 2000.0;
	const REAL8 fmin  = 30.0;
	const REAL8 fmax  = 2200.0;
	UINT4 taperlen = floor( 0.5 + beta * seglen );
	UINT4 stride = seglen - taperlen;

	COMPLEX16FrequencySeries *condresp;
	REAL8TimeSeries *whitened;
	REAL8Window     *window;

	window = XLALCreateTukeyREAL8Window( seglen, beta );
	whitened = XLALCreateREAL8TimeSeries( unwhitened->name, &unwhitened->epoch, unwhitened->f0, unwhitened->deltaT, &unwhitened->sampleUnits, unwhitened->data->length );

	condresp = condition_response( response, flow, fhigh, fmin, fmax, seglen, unwhitened->deltaT );
	datacond2( whitened, unwhitened, window, stride, condresp, fwdplan, revplan, fmin, flow, fhigh, fmax );

	XLALDestroyCOMPLEX16FrequencySeries(condresp);
	XLALDestroyREAL8Window( window );

	/* first and last taperlen/2 points are zero ... remove them */
	/* after that, the first and last stride points are based on
	 * a single psd rather than two ... remove them too */
	XLALResizeREAL8TimeSeries( whitened, taperlen/2 + stride, whitened->data->length - taperlen - 2*stride );

	return whitened;
}
