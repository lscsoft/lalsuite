#include <math.h>
#include <getopt.h>
#include <lal/LALStdio.h>
#include <lal/Units.h>
#include <lal/LALFrameIO.h>
#include <lal/LALCalibration.h>
#include <lal/FrequencySeries.h>
#include <lal/VectorOps.h>
#include <lal/PrintFTSeries.h>

extern int lalDebugLevel;
int main( int argc, char *argv[] )
{
	COMPLEX16FrequencySeries *response;
	REAL8FrequencySeries     *responseabs;
	REAL8FrequencySeries     *responsearg;
	LALCalData *caldata;
	const char *channel = NULL;
	const char *calfile = NULL;
	LIGOTimeGPS epoch = {0, 0};
	REAL8 duration = 0.0;
	REAL8 deltaF = 0.0;
	REAL8 fmax = 0.0;
	char *program = argv[0];
	struct option long_options[] = {
		{ "help",                    no_argument, 0, 'h' },
		{ "gps-start-time",    required_argument, 0, 't' },
		{ "duration",          required_argument, 0, 'd' },
		{ "channel-name",      required_argument, 0, 'c' },
		{ "calibration-file",  required_argument, 0, 'C' },
		{ "frequency-step",    required_argument, 0, 'f' },
		{ "maximum-frequency", required_argument, 0, 'F' },
		{ 0, 0, 0, 0 }
	};
	char args[] = "ht:d:c:C:f:F:";

	lalDebugLevel = 7;

	XLALSetErrorHandler( XLALAbortErrorHandler );

	/* parse options */
	while ( 1 ) {
		int option_index = 0;
		int c;
		c = getopt_long_only( argc, argv, args, long_options, &option_index );
		if ( c == -1 )
			break;
		switch ( c ) {
			case 0: /* if option set a flag, nothing else to do */
				if ( long_options[option_index].flag )
					break;
				fprintf( stderr, "error parsing option %s with argument %s\n", long_options[option_index].name, optarg );
				exit(1);
			case 'h':
				fprintf( stderr, "usage: %s options\n", program );
				fprintf( stderr, "options:\n" );
				fprintf( stderr, "--help                     print this message\n" );
				fprintf( stderr, "--gps-start-time=time      GPS epoch of response function [required]\n" );
				fprintf( stderr, "--duration=dt              duration (seconds) for averaging factors [optional]\n" );
				fprintf( stderr, "--channel-name=chan        readout channel (e.g., \"H1:LSC-DARM_ERR\") [required]\n" );
				fprintf( stderr, "--calibration-file=fname   calibration frame file name [required]\n" );
				fprintf( stderr, "--frequency-step=df        frequency resolution (Hz) of response [optional]\n" );
				fprintf( stderr, "--maximum-frequency=fmax   maximum frequency (Hz) of response [optional]\n" );
				exit(0);
			case 't':
				epoch.gpsSeconds = atol(optarg);
				break;
			case 'd':
				duration = atof(optarg);
				break;
			case 'c':
				channel = optarg;
				break;
			case 'C':
				calfile = optarg;
				break;
			case 'f':
				deltaF = atof(optarg);
				break;
			case 'F':
				fmax = atof(optarg);
				break;
			case '?':
				fprintf(stderr, "unknown error while parsing options\n"); 
				exit(1);
			default:
				fprintf(stderr, "unknown error while parsing options\n"); 
				exit(1);
		}
	}
	if ( optind < argc ) {
		fprintf(stderr, "extraneous command line arguments:\n"); 
		while ( optind < argc )
			fprintf(stderr, "%s\n", argv[optind++]); 
		exit(1);
	}


	/* check that required values have been set */
	if ( ! epoch.gpsSeconds ) {
		fprintf(stderr, "use --gps-start-time to calibration epoch\n"); 
		exit(1);
	}
	if ( ! channel ) {
		fprintf(stderr, "use --channel-name to specify readout channel\n"); 
		exit(1);
	}
	if ( ! calfile ) {
		fprintf(stderr, "use --calibration-file to specify calibration frame file name\n"); 
		exit(1);
	}


	caldata = XLALFrGetCalData(&epoch, channel, calfile);
	if ( ! caldata ) {
		fprintf(stderr, "error: could not get calibration data\n");
		exit(1);
	}

	deltaF   = deltaF > 0.0 ? deltaF : caldata->responseReference->deltaF;
	fmax     = fmax > 0.0 ? fmax : caldata->responseReference->deltaF * caldata->responseReference->data->length;
	response = XLALCreateCOMPLEX16Response( &epoch, duration, deltaF, (UINT4)floor(fmax/deltaF + 0.5), caldata );

	responseabs = XLALCreateREAL8FrequencySeries( "response_abs", &response->epoch, response->f0, response->deltaF, &response->sampleUnits, response->data->length );
	responsearg = XLALCreateREAL8FrequencySeries( "response_arg", &response->epoch, response->f0, response->deltaF, &lalDimensionlessUnit, response->data->length );
	XLALCOMPLEX16VectorAbs( responseabs->data, response->data );
	XLALCOMPLEX16VectorArg( responsearg->data, response->data );

	LALZPrintFrequencySeries( response, "response.dat" );
	LALDPrintFrequencySeries( responseabs, "response_abs.dat" );
	LALDPrintFrequencySeries( responsearg, "response_pvarg.dat" );

	XLALREAL8VectorUnwrapAngle( responsearg->data, responsearg->data );
	LALDPrintFrequencySeries( responsearg, "response_arg.dat" );

	XLALDestroyREAL8FrequencySeries( responsearg );
	XLALDestroyREAL8FrequencySeries( responseabs );
	XLALDestroyCOMPLEX16FrequencySeries( response );
	XLALDestroyCalData( caldata );

	LALCheckMemoryLeaks();
	return 0;
}
