/* vim: set noet ts=4 sw=4: */
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lal/FrameCache.h>
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

enum { ASCII_OUTPUT, WAVE_OUTPUT, AU_OUTPUT };

LIGOTimeGPS start;
REAL8 minfreq, maxfreq;
REAL8 srate;
REAL8 duration = 16.0;
INT4  outtype;
const char *channel = "DARM_ERR";
const char *path    = NULL;
FILE       *outfile;
int cachefile = -1;
UINT4 numave = 0;

extern char *optarg;
extern int optind;

int usage( const char *program );

int parseopts( int argc, char **argv )
{
	struct option long_options[] = {
		{ "help", no_argument, 0, 'h' },
		{ "power-spectrum", required_argument, 0, 'P' },
		{ "gps-start-time", required_argument, 0, 't' },
		{ "duration", required_argument, 0, 'd' },
		{ "min-frequency", required_argument, 0, 'm' },
		{ "max-frequency", required_argument, 0, 'M' },
		{ "channel-name", required_argument, 0, 'c' },
		{ "frame-files", required_argument, 0, 'f' },
		{ "cache-file", required_argument, 0, 'F' },
		{ "output-file", required_argument, 0, 'o' },
		{ "output-type", required_argument, 0, 'O' },
		{ "sample-rate", required_argument, 0, 's' },
		{ 0, 0, 0, 0 } };
	char args[] = "P:t:d:m:M:c:f:F:o:O:s:";

	/* init defaults */
	outfile = stdout;

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
			case 'h':
				usage( argv[0] );
				exit( 0 );
			case 'P':
				numave = atoi( optarg );
				break;
			case 't':
				XLALStrToGPS( &start, optarg, NULL );
				break;
			case 'T':
				start.gpsNanoSeconds = atol( optarg );
				break;
			case 'c':
				channel = optarg;
				break;
			case 'f':
				cachefile = 0;
				path = optarg;
				break;
			case 'F':
				cachefile = 1;
				path = optarg;
				break;
			case 'd':
				duration = atof( optarg );
				break;
			case 'o':
				outfile = fopen( optarg, "w" );
				break;
			case 'O':
				if ( strstr( optarg, "ASC" ) || strstr( optarg, "asc" ) )
					outtype = ASCII_OUTPUT;
				else if ( strstr( optarg, "WAV" ) || strstr( optarg, "wav" ) )
					outtype = WAVE_OUTPUT;
				else if ( strstr( optarg, "AU" ) || strstr( optarg, "au" ) )
					outtype = AU_OUTPUT;
				else {
					fprintf( stderr, "error: unrecognized output type\n" );
					exit( 1 );
				}
				break;
			case 's':
				srate = atof( optarg );
				break;
			case 'm':
				minfreq = atof( optarg );
				break;
			case 'M':
				maxfreq = atof( optarg );
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



REAL4TimeSeries * getdata(
		const char *path,
		int cachefile,
		const char *channel,
		LIGOTimeGPS *start,
		REAL8 duration )
{
	FrCache *cache   = NULL;
	FrStream *stream = NULL;
	int mode = LAL_FR_VERBOSE_MODE;
	REAL4TimeSeries *series;

	/* open the data cache and use it to get a frame stream */
	if ( cachefile )
		cache = XLALFrImportCache( path );
	else {
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
		cache = XLALFrGenerateCache( dirname, basename );
	}
	stream = XLALFrCacheOpen( cache );
	XLALFrDestroyCache( cache );

	/* set the mode of the frame stream */
	XLALFrSetMode( stream, mode );

	series = XLALFrReadREAL4TimeSeries( stream, channel, start, duration, 0 );

	/* close the stream */
	XLALFrClose( stream );

	return series;
}

int resample( REAL4TimeSeries *series, REAL8 srate )
{
	if ( srate > 0 )
		XLALResampleREAL4TimeSeries( series, 1.0/srate );
	return 0;
}

int filter( REAL4TimeSeries *series, REAL8 minfreq, REAL8 maxfreq )
{
	const INT4  filtorder = 8;
	const REAL8 amplitude = 0.9; /* 10% attenuation at specified frequency */
	UINT4 sec;

	if ( minfreq > 0 )
		XLALHighPassREAL4TimeSeries( series, minfreq, amplitude, filtorder );
	if ( maxfreq > 0 )
		XLALLowPassREAL4TimeSeries( series, minfreq, amplitude, filtorder );

	/* number of points in a second */
	sec = ceil(1.0/series->deltaT);

	/* resize the array by discarding first and last second (do this always) */
	XLALResizeREAL4TimeSeries( series, sec, series->data->length - 2*sec);

	return 0;
}

REAL4FrequencySeries * powerspec( REAL4TimeSeries *series, REAL8 segdur, LIGOTimeGPS *epoch )
{
	REAL4FrequencySeries *spectrum;
	REAL4FFTPlan *plan;
	REAL4Window *window;
	UINT4 seglen;
	UINT4 stride;
	seglen = segdur / series->deltaT;
	stride = seglen / 2;
	
	spectrum = XLALCreateREAL4FrequencySeries( "untitled", epoch, 0.0, 1.0/series->deltaT, &lalDimensionlessUnit, seglen/2 + 1 );
	plan = XLALCreateForwardREAL4FFTPlan( seglen, 0 );
	window = XLALCreateHannREAL4Window( seglen );
	XLALREAL4AverageSpectrumWelch( spectrum, series, seglen, stride, window, plan );
	XLALDestroyREAL4Window( window );
	XLALDestroyREAL4FFTPlan( plan );

	return spectrum;
}

int output_psd( FILE *outfile, REAL4FrequencySeries *series )
{
	UINT4 i;
	/* note: omit DC */
	for ( i = 1; i < series->data->length; ++i )
		fprintf( outfile, "%e\t%e\n", i * series->deltaF, sqrt( series->data->data[i] ) );
	return 0;
}

int output( FILE *outfile, int outtype, REAL4TimeSeries *series )
{
	UINT4 i;
	int nfreq;

	switch ( outtype ) {
		case ASCII_OUTPUT:
			for ( i = 0; i < series->data->length; ++i ) 
				fprintf( outfile, "%20.9f\t%e\n", i * series->deltaT, series->data->data[i] );
			break;
		case WAVE_OUTPUT:
			nfreq = 1.0/series->deltaT;
			XLALAudioWAVRecordREAL4TimeSeries( outfile, series );
			break;
		case AU_OUTPUT:
			nfreq = 1.0/series->deltaT;
			XLALAudioAURecordREAL4TimeSeries( outfile, series );
			break;
		default:
			fprintf( stderr, "error: invalid output type\n" );
			exit( 1 );
	}
	return 0;
}


int main( int argc, char *argv[] )
{
	REAL4TimeSeries *series;
	REAL8 segdur = 0.0;

	lalDebugLevel = 1;
	XLALSetErrorHandler( XLALAbortErrorHandler );

	parseopts( argc, argv );
	if ( numave ) { /* power spectrum */
		segdur = duration;
		duration = ( duration * ( numave + 1 ) ) / 2;
	}
	duration += 2; /* because we'll discard two seconds */
	series = getdata( path, cachefile, channel, &start, duration );
	resample( series, srate );
	filter( series, minfreq, maxfreq );
	if ( numave ) { /* power spectrum */
		REAL4FrequencySeries *spectrum;
		spectrum = powerspec( series, segdur, &start );
		output_psd( outfile, spectrum );
	} else
		output( outfile, outtype, series );
	return 0;
}


int usage( const char * program )
{
	fprintf( stderr, "usage: %s [options]\n", program );

	fprintf( stderr, "\t-c CHANNEL\n\t--channel-name=CHANNEL\n" );
	fprintf( stderr, "\t\tread channel CHANNEL\n" );

	fprintf( stderr, "\t-d DURATION\n\t--duration=DURATION\n" );
	fprintf( stderr, "\t\tread DURATION seconds of data\n" );

	fprintf( stderr, "\t-f FILES\n\t--frame-files=FILES\n" );
	fprintf( stderr, "\t\tread data from frame files FILES\n" );

	fprintf( stderr, "\t-F CACHE\n\t--cache-file=CACHE\n" );
	fprintf( stderr, "\t\tread data from files in frame cache file CACHE\n" );

	fprintf( stderr, "\t-h\n\t--help\n\t\tprint this message\n\n" );

	fprintf( stderr, "\t-m MINFREQ\n\t--min-frequency=MINFREQ\n" );
	fprintf( stderr, "\t\thighpass data at MINFREQ hertz\n" );

	fprintf( stderr, "\t-M MAXFREQ\n\t--max-frequency=MAXNFREQ\n" );
	fprintf( stderr, "\t\tlowpass data at MAXFREQ hertz\n" );

	fprintf( stderr, "\t-o OUTFILE\n\t--output-file=OUTFILE\n" );
	fprintf( stderr, "\t\toutput to file OUTFILE\n" );

	fprintf( stderr, "\t-O OTYPE\n\t--output-type=OTYPE\n" );
	fprintf( stderr, "\t\toutput data type OTYPE [ \"ASCII\" | \"AU\" | \"WAVE\" ]\n" );

	fprintf( stderr, "\t-P NUMAVG\n\t--power-spectrum=NUMAVE\n" );
	fprintf( stderr, "\t\tcompute power spectrum with NUMAVE averages\n" );

	fprintf( stderr, "\t-s SRATE\n\t--sample-rate=SRATE\n" );
	fprintf( stderr, "\t\tresample to sample rate SRATE hertz\n" );

	fprintf( stderr, "\t-t GPSTIME\n\t--gps-start-time=GPSTIME\n" );
	fprintf( stderr, "\t\tstart reading data at time GPSTIME\n" );
	fprintf( stderr, "\t\tnote: actual start time is one second later!\n" );

	return 0;
}
