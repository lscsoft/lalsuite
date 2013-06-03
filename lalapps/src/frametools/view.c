/*
*  Copyright (C) 2007 Jolien Creighton
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

/* vim: set noet ts=4 sw=4: */
#include <complex.h>
#include <getopt.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lal/LALFrameIO.h>
#include <lal/LALCache.h>
#include <lal/FrameStream.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/Audio.h>
#include <lal/Window.h>
#include <lal/RealFFT.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataBurstUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/FindChirp.h>
#include <lal/GenerateBurst.h>



/* output type enumeration */
enum { ASCII_OUTPUT, WAVE_OUTPUT, AU_OUTPUT, PSD_OUTPUT };

/* input data enumeration */
enum { READ_INPUT, ZERO_INPUT, IMPULSE_INPUT };

LIGOTimeGPS start_;
REAL8 minfreq_;
REAL8 maxfreq_;
REAL8 srate_;
REAL8 duration_ = 16.0;
const char *channel_ = NULL;
const char *datafile_ = NULL;
const char *calibfile_ = NULL;
const char *inspinjfile_ = NULL;
const char *burstinjfile_ = NULL;
const char *outfile_ = NULL;
UINT4 numave_ = 0;
REAL8 dynrange_ = 1;
int texactflg_ = 0;
int cachefileflg_ = -1;
int calibflg_ = 0;
int outtype_ = ASCII_OUTPUT;
int intype_ = READ_INPUT;
int verboseflg_ = 0;
int debugflg_ = 0;

extern char *optarg;
extern int optind;

int verbose( const char *fmt, ... );
int usage( const char *program );
int parseopts( int argc, char **argv );
int dbg_tsdump( REAL4TimeSeries *series, const char *fname );
int dbg_fsdump( COMPLEX8FrequencySeries *series, const char *fname );
int dbg_specdump( REAL4FrequencySeries *series, const char *fname );
REAL4TimeSeries * setdata( int intype, const char *channel, LIGOTimeGPS *start, REAL8 duration, REAL8 samplerate );
REAL4TimeSeries * getdata( const char *path, int cachefileflg, const char *channel, LIGOTimeGPS *start, REAL8 duration );
// COMPLEX8FrequencySeries * getresp( const char *calibfile, const char *channel, LIGOTimeGPS *epoch, REAL4 deltaF, UINT4 seglen, REAL8 dynrange );
int inspinj( REAL4TimeSeries *series, const char *inspinjfile, const char *calfile );
int burstinj( REAL4TimeSeries *series, const char *burstinjfile, const char *calfile );
int resample( REAL4TimeSeries *series, REAL8 srate );
int filter( REAL4TimeSeries *series, REAL8 minfreq, REAL8 maxfreq );
int calibrate( REAL4TimeSeries *tseries, const char *calfile, REAL8 f_min );
REAL4FrequencySeries * powerspec( REAL4TimeSeries *series, REAL8 segdur, LIGOTimeGPS *epoch, const char *calibfile, int intype );
int output_psd( const char *outfile, REAL4FrequencySeries *series );
int output( const char *outfile, int outtype, REAL4TimeSeries *series );


int main( int argc, char *argv[] )
{
	REAL4TimeSeries *series;
	REAL8 segdur = 0.0;
	REAL8 discard;
	UINT4 ndiscard;

	XLALSetErrorHandler( XLALAbortErrorHandler );

	parseopts( argc, argv );

	discard = duration_;
	if ( outtype_ == PSD_OUTPUT ) {
		/* power spectrum: acutal duration required is more */
		segdur     = duration_;
		duration_  = ( duration_ * ( numave_ + 1 ) ) / 2;
	}
	duration_ += discard; /* we'll discard data from beginning and end */

	/* if we want result to start at exactly the right time, need to start
	 * reading so much earlier */
	if ( texactflg_ )
		XLALGPSAdd( &start_, -0.5*discard );

	/* get the data */
	if ( intype_ == READ_INPUT )
		series = getdata( datafile_, cachefileflg_, channel_, &start_, duration_ );
	else
		series = setdata( intype_, channel_, &start_, duration_, srate_ ? srate_ : 16384.0 );

	dbg_tsdump( series, "tseries0.dat" );

	/* inject inspirals and bursts */
	inspinj( series, inspinjfile_, calibfile_ );
	burstinj( series, burstinjfile_, calibfile_ );

	dbg_tsdump( series, "tseries1.dat" );

	/* condition the data */
	resample( series, srate_ );
	filter( series, minfreq_, maxfreq_ );
	if ( calibflg_ && outtype_ != PSD_OUTPUT ) /* calibrate elsewhere if doing a psd */
		calibrate( series, calibfile_, minfreq_ );

	dbg_tsdump( series, "tseries2.dat" );

	/* discarding extra data from beginning and end */
	ndiscard = 0.5 * discard / series->deltaT;
	XLALResizeREAL4TimeSeries( series, ndiscard, series->data->length - 2*ndiscard);

	/* output results */
	if ( outtype_ == PSD_OUTPUT ) { /* power spectrum */
		REAL4FrequencySeries *spectrum;
		if ( calibflg_ )
			spectrum = powerspec( series, segdur, &start_, calibfile_, intype_ );
		else
			spectrum = powerspec( series, segdur, &start_, NULL, intype_ );
		output_psd( outfile_, spectrum );
	} else
		output( outfile_, outtype_, series );

	return 0;
}


int parseopts( int argc, char **argv )
{
	struct option long_options[] = {
		{ "help", no_argument, 0, 'h' },
		{ "power-spectrum", required_argument, 0, 'P' },
		{ "gps-start-time", required_argument, 0, 't' },
		{ "exact-gps-start-time", required_argument, 0, 'T' },
		{ "duration", required_argument, 0, 'd' },
		{ "min-frequency", required_argument, 0, 'm' },
		{ "max-frequency", required_argument, 0, 'M' },
		{ "channel-name", required_argument, 0, 'c' },
		{ "calibration-file", required_argument, 0, 'C' },
		{ "frame-files", required_argument, 0, 'f' },
		{ "cache-file", required_argument, 0, 'F' },
		{ "output-file", required_argument, 0, 'o' },
		{ "output-type", required_argument, 0, 'O' },
		{ "sample-rate", required_argument, 0, 's' },
		{ "inject-inspiral", required_argument, 0, 'I' },
		{ "inject-burst", required_argument, 0, 'B' },
		{ "unit-impulse", no_argument, 0, '1' },
		{ "zero-data", no_argument, 0, '0' },
		{ "no-calibration", no_argument, &calibflg_, 0 },
		{ "verbose", no_argument, &verboseflg_, 1 },
		{ "debug", no_argument, &debugflg_, 1 },
		{ 0, 0, 0, 0 } };
	char args[] = "01P:t:T:d:m:M:c:C:f:F:o:O:s:I:B:";

	while ( 1 ) {
		int option_index = 0;
		int c;
		c = getopt_long_only( argc, argv, args, long_options, &option_index );
		if ( c == -1 )
			break;
		switch ( c ) {
			case 0:
				if ( long_options[option_index].flag )
					break; /* option set flag: nothing else to do */
				fprintf( stderr, "error parsing option %s with argument %s\n", long_options[option_index].name, optarg );
				exit( 1 );
			case '0':
				intype_ = ZERO_INPUT;
				break;
			case '1':
				intype_ = IMPULSE_INPUT;
				break;
			case 'h':
				usage( argv[0] );
				exit( 0 );
			case 'P':
				numave_ = atoi( optarg );
				outtype_ = PSD_OUTPUT;
				break;
			case 't':
				XLALStrToGPS( &start_, optarg, NULL );
				texactflg_ = 0;
				break;
			case 'T':
				XLALStrToGPS( &start_, optarg, NULL );
				texactflg_ = 1;
				break;
			case 'c':
				channel_ = optarg;
				break;
			case 'C':
				calibfile_ = optarg;
				calibflg_ = 1;
				break;
			case 'f':
				cachefileflg_ = 0;
				datafile_ = optarg;
				break;
			case 'F':
				cachefileflg_ = 1;
				datafile_ = optarg;
				break;
			case 'I':
				inspinjfile_ = optarg;
				break;
			case 'B':
				burstinjfile_ = optarg;
				break;
			case 'd':
				duration_ = atof( optarg );
				break;
			case 'o':
				outfile_ = optarg;
				break;
			case 'O':
				if ( strstr( optarg, "ASC" ) || strstr( optarg, "asc" ) )
					outtype_ = ASCII_OUTPUT;
				else if ( strstr( optarg, "WAV" ) || strstr( optarg, "wav" ) )
					outtype_ = WAVE_OUTPUT;
				else if ( strstr( optarg, "AU" ) || strstr( optarg, "au" ) )
					outtype_ = AU_OUTPUT;
				else {
					fprintf( stderr, "error: unrecognized output type\n" );
					exit( 1 );
				}
				break;
			case 's':
				srate_ = atof( optarg );
				break;
			case 'm':
				minfreq_ = atof( optarg );
				break;
			case 'M':
				maxfreq_ = atof( optarg );
				break;
			case '?':
			default:
				fprintf( stderr, "unknown error while parsing options\n" );
				exit( 1 );
		}
	}

	if ( optind < argc ) {
		fprintf( stderr, "extraneous command line arguments:\n" );
		while ( optind < argc )
			fprintf( stderr, "%s\n", argv[optind++] );
		exit( 1 );
	}

	return 0;
}


REAL4TimeSeries * setdata( int intype, const char *channel, LIGOTimeGPS *start, REAL8 duration, REAL8 samplerate )
{
	REAL4TimeSeries *series;
	UINT4 length;
	UINT4 i;

	verbose( "generating %g seconds of %s data starting at time %d.%09d\n",
			duration, channel, start->gpsSeconds, start->gpsNanoSeconds );

	length = duration * samplerate;
	series = XLALCreateREAL4TimeSeries( channel, start, 0, 1.0/samplerate, &lalADCCountUnit, length );
	for ( i = 0; i < series->data->length; ++i )
		series->data->data[i] = 0;
	if ( intype == IMPULSE_INPUT )
		series->data->data[series->data->length/2] = 1.0;
	return series;
}

REAL4TimeSeries * getdata( const char *path, int cachefileflg, const char *channel, LIGOTimeGPS *start, REAL8 duration )
{
	LALCache *cache   = NULL;
	FrStream *stream = NULL;
	int mode = LAL_FR_VERBOSE_MODE;
	int tstype;
	REAL4TimeSeries *series;

	verbose( "reading %g seconds of %s data starting at time %d.%09d\n",
			duration, channel, start->gpsSeconds, start->gpsNanoSeconds );

	/* open the data cache and use it to get a frame stream */
	if ( cachefileflg ) {
		verbose( "importing cache file %s\n", path );
		cache = XLALCacheImport( path );
	} else {
		char pathcpy[FILENAME_MAX];
		char *basename;
		char *dirname = NULL;
		strncpy( pathcpy, path, sizeof(pathcpy) - 1 );
		basename = strrchr( pathcpy, '/' );
		if ( basename ) {
			*basename++ = 0;
			dirname = pathcpy;
		} else {
			basename = pathcpy;
			dirname = NULL;
		}
		verbose( "using data file(s) %s\n", path );
		cache = XLALCacheGlob( dirname, basename );
	}
	stream = XLALFrCacheOpen( cache );
	XLALDestroyCache( cache );

	/* set the mode of the frame stream */
	XLALFrSetMode( stream, mode );

	tstype = XLALFrGetTimeSeriesType( channel, stream );
	if ( tstype == LAL_S_TYPE_CODE ) {
		series = XLALFrReadREAL4TimeSeries( stream, channel, start, duration, 0 );
	} else if ( tstype == LAL_D_TYPE_CODE ) { /* assume strain data */
		REAL8TimeSeries *dblseries;
		UINT4 i;
		dynrange_ = 1e20;
		dblseries = XLALFrReadREAL8TimeSeries( stream, channel, start, duration, 0 );
		/* TODO: this shouldn't be hard-coded! */
		XLALHighPassREAL8TimeSeries( dblseries, 40.0, 0.9, 8 );
		series = XLALCreateREAL4TimeSeries( dblseries->name, &dblseries->epoch, dblseries->f0, dblseries->deltaT, &dblseries->sampleUnits, dblseries->data->length );
		for ( i = 0; i < series->data->length; ++i )
			series->data->data[i] = dynrange_ * dblseries->data->data[i];
	} else {
		return NULL;
	}

	/* close the stream */
	XLALFrClose( stream );

	return series;
}

#if 0
COMPLEX8FrequencySeries * getresp( const char *calibfile, const char *channel, LIGOTimeGPS *epoch, REAL4 deltaF, UINT4 seglen, REAL8 dynrange )
{
	if ( ! calibfile )
		return NULL;

        XLAL_ERROR_NULL(XLAL_EERR, "Calibration frames no longer supported");
}
#endif


int inspinj( REAL4TimeSeries *series, const char *inspinjfile, const char *calfile )
{
	static LALStatus status;
	SimInspiralTable *injections = NULL;
	COMPLEX8FrequencySeries *response;
	REAL4TimeSeries keep;
	REAL8 deltaF;
	int ninj;
	int tbeg;
	int tend;

	if ( ! inspinjfile )
		return 0;
	if ( ! calfile ) {
		fprintf( stderr, "warning: cannot perform injections without calibration\n" );
		return 0;
	}

	tbeg = series->epoch.gpsSeconds;
	tend = tbeg + ceil( series->deltaT * series->data->length );
	ninj = SimInspiralTableFromLIGOLw( &injections, inspinjfile, tbeg, tend );

	verbose( "injecting %d inspirals listed in file %s between times %d and %d\n", ninj, inspinjfile, tbeg, tend );

	if ( ninj < 0 ) {
		fprintf( stderr, "error: could not read file %s\n", inspinjfile );
		exit( 1 );
	} else if ( ninj == 0 ) {
		fprintf( stderr, "warning: no relevant injections found in %s\n", inspinjfile );
		return 0;
	}

	deltaF = 1.0 / ( series->deltaT * series->data->length );
	// response = getresp( calfile, series->name, &series->epoch, deltaF, series->data->length, 1.0 );
	response = NULL;
	keep = *series;
	LALFindChirpInjectSignals( &status, series, injections, response );
	*series = keep;
	if ( status.statusCode ) {
		fprintf( stderr, "error: signal injection failed\n" );
		exit( 1 );
	}

	while ( injections ) {
		SimInspiralTable *thisInjection = injections;
		injections = injections->next;
		LALFree( thisInjection );
	}
	XLALDestroyCOMPLEX8FrequencySeries( response );

	return 0;
}


int burstinj( REAL4TimeSeries *series, const char *burstinjfile, const char *calfile )
{
	TimeSlide *time_slide;
	SimBurst *sim_burst;
	REAL8TimeSeries *injections;
	COMPLEX8FrequencySeries *response;
	REAL8 deltaF;
	LIGOTimeGPS tbeg = series->epoch;
	LIGOTimeGPS tend = series->epoch;
	unsigned i;

	if ( ! burstinjfile )
		return 0;
	if ( ! calfile ) {
		fprintf( stderr, "warning: cannot perform injections without calibration\n" );
		return 0;
	}

	XLALGPSAdd(&tend, series->deltaT * series->data->length);

	time_slide = XLALTimeSlideTableFromLIGOLw( burstinjfile );
	sim_burst = XLALSimBurstTableFromLIGOLw( burstinjfile, &tbeg, &tend );
	if ( !sim_burst || !time_slide ) {
		fprintf( stderr, "error: could not read file %s\n", burstinjfile );
		exit( 1 );
	}

	verbose( "injecting bursts listed in file %s between times %d and %d\n", burstinjfile, tbeg, tend );

	deltaF = 1.0 / ( series->deltaT * series->data->length );
	// response = getresp( calfile, series->name, &series->epoch, deltaF, series->data->length, 1.0 );
	response = NULL;

	injections = XLALCreateREAL8TimeSeries(series->name, &series->epoch, series->f0, series->deltaT, &series->sampleUnits, series->data->length);
	/* FIXME:  new injection code requires double precision respose */
	//if(XLALBurstInjectSignals( injections, sim_burst, time_slide, /*response*/ NULL )) {
	//	fprintf( stderr, "error: signal injection failed\n" );
	//	exit( 1 );
	//}
	for(i = 0; i < series->data->length; i++)
		series->data->data[i] += injections->data->data[i];
	XLALDestroyREAL8TimeSeries(injections);

	XLALDestroyTimeSlideTable(time_slide);
	XLALDestroySimBurstTable(sim_burst);
	XLALDestroyCOMPLEX8FrequencySeries( response );

	return 0;
}


int resample( REAL4TimeSeries *series, REAL8 srate )
{
	if ( srate > 0 ) {
		verbose( "resampling to rate %g Hz\n", srate );
		XLALResampleREAL4TimeSeries( series, 1.0/srate );
	}
	return 0;
}

int filter( REAL4TimeSeries *series, REAL8 minfreq, REAL8 maxfreq )
{
	const INT4  filtorder = 8;
	const REAL8 amplitude = 0.9; /* 10% attenuation at specified frequency */

	if ( minfreq > 0 ) {
		verbose( "high-pass filtering at frequency %g Hz\n", minfreq );
		XLALHighPassREAL4TimeSeries( series, minfreq, amplitude, filtorder );
	}
	if ( maxfreq > 0 ) {
		verbose( "low-pass filtering at frequency %g Hz\n", maxfreq );
		XLALLowPassREAL4TimeSeries( series, maxfreq, amplitude, filtorder );
	}

	return 0;
}

int calibrate( REAL4TimeSeries *tseries, const char *calfile, REAL8 f_min )
{
	COMPLEX8FrequencySeries *response;
	COMPLEX8FrequencySeries *fseries;
	REAL4FFTPlan *pfwd, *prev;
	REAL8 deltaF;
	UINT4 kmin, k;

	if ( ! calfile )
		return 0;

	verbose( "calibrating data\n" );

	pfwd = XLALCreateForwardREAL4FFTPlan( tseries->data->length, 0 );
	prev = XLALCreateReverseREAL4FFTPlan( tseries->data->length, 0 );

	fseries = XLALCreateCOMPLEX8FrequencySeries( tseries->name, &tseries->epoch, 0, 1.0/tseries->deltaT, &lalDimensionlessUnit, tseries->data->length/2 + 1 );
	XLALREAL4TimeFreqFFT( fseries, tseries, pfwd );

	deltaF = 1.0 / ( tseries->deltaT * tseries->data->length );
	dynrange_ = 1e20;
	// response = getresp( calfile, tseries->name, &tseries->epoch, deltaF, tseries->data->length, dynrange_ );
	response = NULL;

	if ( f_min < 30.0 ) {
		fprintf( stderr, "warning: setting minimum frequency to 30 Hz for calibration\n" );
		f_min = 30.0;
	}
	kmin = f_min / fseries->deltaF;
	for ( k = 0; k < kmin; ++k )
		fseries->data->data[k] = 0.0;
	for ( k = kmin; k < fseries->data->length; ++k )
		fseries->data->data[k] = fseries->data->data[k] * response->data->data[k];
	XLALREAL4FreqTimeFFT( tseries, fseries, prev );

	XLALDestroyCOMPLEX8FrequencySeries( fseries );
	XLALDestroyCOMPLEX8FrequencySeries( response );
	XLALDestroyREAL4FFTPlan( prev );
	XLALDestroyREAL4FFTPlan( pfwd );
	return 0;
}

REAL4FrequencySeries * powerspec( REAL4TimeSeries *series, REAL8 segdur, LIGOTimeGPS *epoch, const char *calibfile, int intype )
{
	REAL4FrequencySeries *spectrum;
	UINT4 seglen;
	UINT4 stride;
	seglen = segdur / series->deltaT;
	stride = seglen / 2;


	spectrum = XLALCreateREAL4FrequencySeries( "spectrum", epoch, 0.0, 1.0/segdur, &lalDimensionlessUnit, seglen/2 + 1 );
	if ( intype != IMPULSE_INPUT ) {
		REAL4FFTPlan *plan;
		REAL4Window *window;
		verbose( "computing power spectrum\n" );
		plan = XLALCreateForwardREAL4FFTPlan( seglen, 0 );
		window = XLALCreateHannREAL4Window( seglen );
		XLALREAL4AverageSpectrumWelch( spectrum, series, seglen, stride, window, plan );
		XLALDestroyREAL4Window( window );
		XLALDestroyREAL4FFTPlan( plan );
	} else { /* impulse input: set spectrum to unity */
		UINT4 k;
		verbose( "setting spectrum to be unity\n" );
		spectrum->data->data[0] = 0.0;
		for ( k = 1; k < spectrum->data->length; ++k )
			spectrum->data->data[k] = 1.0;
	}

	if ( ! calibfile )
		return spectrum;
	else {
		COMPLEX8FrequencySeries *response;
		UINT4 k;
		dynrange_ = 1e20;
		// response = getresp( calibfile, series->name, &series->epoch, spectrum->deltaF, seglen, dynrange_ );
		response = NULL;
		for ( k = 1; k < spectrum->data->length; ++k ) {
			REAL4 re = crealf(response->data->data[k]);
			REAL4 im = cimagf(response->data->data[k]);
			spectrum->data->data[k] *= re*re + im*im;
		}
		XLALDestroyCOMPLEX8FrequencySeries( response );
	}

	return spectrum;
}

int output_psd( const char *outfile, REAL4FrequencySeries *series )
{
	UINT4 i;
	FILE *fp;
	fp = outfile ? fopen( outfile, "w" ) : stdout;

	verbose( "output psd of %s data beginning at GPS time %d.%09d to file %s\n",
			series->name, series->epoch.gpsSeconds,
			series->epoch.gpsNanoSeconds, outfile );

	/* note: omit DC */
	for ( i = 1; i < series->data->length; ++i )
		fprintf( fp, "%e\t%e\n", i * series->deltaF, sqrt( series->data->data[i] ) / dynrange_ );

	if ( fp != stdout )
		fclose( fp );
	return 0;
}

int output( const char *outfile, int outtype, REAL4TimeSeries *series )
{
	UINT4 i;
	FILE *fp;
	fp = outfile ? fopen( outfile, "w" ) : stdout;

	verbose( "output %s data beginning at GPS time %d.%09d to file %s\n",
			series->name, series->epoch.gpsSeconds,
			series->epoch.gpsNanoSeconds, outfile );

	switch ( outtype ) {
		case ASCII_OUTPUT:
			for ( i = 0; i < series->data->length; ++i ) {
				REAL8 val;
				val  = series->data->data[i];
				val /= dynrange_;
				fprintf( fp, "%.9f\t%e\n", i * series->deltaT, val );
			}
			break;
		case WAVE_OUTPUT:
			XLALAudioWAVRecordREAL4TimeSeries( fp, series );
			break;
		case AU_OUTPUT:
			XLALAudioAURecordREAL4TimeSeries( fp, series );
			break;
		default:
			fprintf( stderr, "error: invalid output type\n" );
			exit( 1 );
	}

	if ( fp != stdout )
		fclose( fp );
	return 0;
}

int usage( const char * program )
{
	fprintf( stderr, "USAGE\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t%s [options]\n", program );

	fprintf( stderr, "\n" );
	fprintf( stderr, "OPTIONS\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-c CHANNEL\n\t--channel-name=CHANNEL\n" );
	fprintf( stderr, "\t\tread channel CHANNEL\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-C CALFILE\n\t--calibration-file=CALFILE\n" );
	fprintf( stderr, "\t\tuse calibration file CALFILE\n" );
	fprintf( stderr, "\t\tthis means that data will be calibrated\n" );
	fprintf( stderr, "\t\tunless the --no-calibration option is used\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-d DURATION\n\t--duration=DURATION\n" );
	fprintf( stderr, "\t\tread DURATION seconds of data\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-f FILES\n\t--frame-files=FILES\n" );
	fprintf( stderr, "\t\tread data from frame files FILES\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-F CACHE\n\t--cache-file=CACHE\n" );
	fprintf( stderr, "\t\tread data from files in frame cache file CACHE\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-h\n\t--help\n\t\tprint this message\n\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-m MINFREQ\n\t--min-frequency=MINFREQ\n" );
	fprintf( stderr, "\t\thighpass data at MINFREQ hertz\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-M MAXFREQ\n\t--max-frequency=MAXNFREQ\n" );
	fprintf( stderr, "\t\tlowpass data at MAXFREQ hertz\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-o OUTFILE\n\t--output-file=OUTFILE\n" );
	fprintf( stderr, "\t\toutput to file OUTFILE\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-O OTYPE\n\t--output-type=OTYPE\n" );
	fprintf( stderr, "\t\toutput data type OTYPE [ \"ASCII\" | \"AU\" | \"WAVE\" ]\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-P NUMAVG\n\t--power-spectrum=NUMAVE\n" );
	fprintf( stderr, "\t\tcompute power spectrum with NUMAVE averages\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-s SRATE\n\t--sample-rate=SRATE\n" );
	fprintf( stderr, "\t\tresample to sample rate SRATE hertz\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-t GPSTIME\n\t--gps-start-time=GPSTIME\n" );
	fprintf( stderr, "\t\tstart reading data at time GPSTIME\n" );
	fprintf( stderr, "\t\tnote: output results for data a little later\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t-T GPSTIME\n\t--exact-gps-start-time=GPSTIME\n" );
	fprintf( stderr, "\t\tdata output for EXACTLY time GPSTIME\n" );
	fprintf( stderr, "\t\tnote: data must exist before this time\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t--no-calibration\n" );
	fprintf( stderr, "\t\tdo not apply response function to data\n" );
	fprintf( stderr, "\t\t(it may still be required for injections)\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t--verbose\n" );
	fprintf( stderr, "\t\tprint messages describing actions that are taken\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\t--debug\n" );
	fprintf( stderr, "\t\tdump intermediate products to data files\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "EXAMPLES\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\tRead some data, condition it, and write it as an audio file\n" );
	fprintf( stderr, "\n\t\t%s --verbose --channel L1:LSC-DARM_ERR --frame-files=L-RDS_R_L3-795259680-256.gwf --min-frequency=80 --max-frequency=500 --sample-rate=1024 --output-file=data.au --output-type=au --gps-start-time=795259680\n", program );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\tImpulse response of the response function\n" );
	fprintf( stderr, "\n\t\t%s --verbose --channel L1:LSC-DARM_ERR --calibration-file=L-L1_CAL_S4_V4-793128493-2801040.gwf --min-frequency=30 --max-frequency=500 --unit-impulse --output-file=resp.dat --gps-start-time=795259680", program );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\tPower spectrum of the response function\n" );
	fprintf( stderr, "\n\t\t%s --verbose --channel L1:LSC-DARM_ERR --calibration-file=L-L1_CAL_S4_V4-793128493-2801040.gwf --unit-impulse --output-file=rpsd.dat --gps-start-time=795259680 --power-spectrum=1", program );

	fprintf( stderr, "\n" );
	fprintf( stderr, "\tMake a strain sensitivity curve with 32 averages\n" );
	fprintf( stderr, "\n\t\t%s --channel=L1:LSC-DARM_ERR --calibration-file=L-L1_CAL_S4_V4-793128493-2801040.gwf --duration=16 --frame-file=\"L-RDS_R_L3-*.gwf\" --output-file=spec.dat --gps-start-time=795259680 --power-spectrum=32\n", program );
	fprintf( stderr, "\n\t(note quotes around wildcard in --frame-file option)\n" );
	fprintf( stderr, "\n" );
	fprintf( stderr, "\tMake an audio file of a chirp injected into strain data\n"
			"\tresampled to 1024 Hz and high-pass filtered at 70 Hz\n" );
	fprintf( stderr, "\n\t\t%s -c L1:LSC-DARM_ERR -C L-L1_CAL_S4_V4-793128493-2801040.gwf -d 16 -f L-RDS_R_L3-795259680-256.gwf -m 70 -o chirph.au -s 1024 -t 795259680 -I inj.xml -O au\n", program );
	fprintf( stderr, "\n" );
	fprintf( stderr, "\tAs above but in terms of ADC counts rather than strain\n" );
	fprintf( stderr, "\n\t\t%s -c L1:LSC-DARM_ERR -C L-L1_CAL_S4_V4-793128493-2801040.gwf -d 16 -f L-RDS_R_L3-795259680-256.gwf -m 70 -o chirpv.au -s 1024 -t 795259680 -I inj.xml -O au --no-cal\n", program );
	fprintf( stderr, "\n" );
	fprintf( stderr, "\tAs above but just with the injection (no data!)\n" );
	fprintf( stderr, "\n\t\t%s -c L1:LSC-DARM_ERR -C L-L1_CAL_S4_V4-793128493-2801040.gwf -d 16 -m 70 -o chirph.au -s 1024 -t 795259680 -I inj.xml -O au --no-cal -0\n", program );

	return 0;
}


int dbg_tsdump( REAL4TimeSeries *series, const char *fname )
{
	if ( debugflg_ ) {
		UINT4 j;
		FILE *fp;
		fp = fopen( fname, "w" );
		for ( j = 0; j < series->data->length; ++j )
			fprintf( fp, "%.9f\t%e\n", j * series->deltaT, series->data->data[j] );
		fclose( fp );
	}
	return 0;
}

int dbg_fsdump( COMPLEX8FrequencySeries *series, const char *fname )
{
	if ( debugflg_ ) {
		UINT4 k;
		FILE *fp;
		fp = fopen( fname, "w" );
		for ( k = 1; k < series->data->length; ++k )
			fprintf( fp, "%e\t%e\t%e\n", k * series->deltaF,
					cabsf(series->data->data[k]), cargf(series->data->data[k]) );
		fclose( fp );
	}
	return 0;
}

int dbg_specdump( REAL4FrequencySeries *series, const char *fname )
{
	if ( debugflg_ ) {
		UINT4 k;
		FILE *fp;
		fp = fopen( fname, "w" );
		for ( k = 1; k < series->data->length; ++k )
			fprintf( fp, "%e\t%e\n", k * series->deltaF, series->data->data[k] );
		fclose( fp );
	}
	return 0;
}

int verbose( const char *fmt, ... )
{
	if ( verboseflg_ ) {
		va_list ap;
		va_start( ap, fmt );
		fprintf( stderr, "verbose: " );
		vfprintf( stderr, fmt, ap );
		va_end( ap );
	}
	return 0;
}
