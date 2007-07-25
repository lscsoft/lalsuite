#include <math.h>
#include <stdio.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/TimeSeries.h>
#include <lal/FrameStream.h>
#include <lal/PrintFTSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/LALCalibration.h>
#include <lal/LALFrameIO.h>
#include <lal/LALDetectors.h>
#include <lal/SkyCoordinates.h>
#include "whiten.h"
#include "crosscorr.h"
#include "signal.h"
#include "inject.h"

#define CHANNEL_TYPE "LSC-STRAIN"
/* #define CHANNEL_TYPE "LSC-DARM_ERR" */
#define IFO_1 "H1"
#define IFO_2 "H2"
#define CHANNEL_1 IFO_1 ":" CHANNEL_TYPE
#define CHANNEL_2 IFO_2 ":" CHANNEL_TYPE

COMPLEX8FrequencySeries * getresp( const char *fname, const char *channel, LIGOTimeGPS *epoch, REAL8 deltaF, UINT4 seglen );
int inject( REAL8TimeSeries *raw_series_1, REAL8TimeSeries *raw_series_2, COMPLEX8FrequencySeries *response_1, COMPLEX8FrequencySeries *response_2, LIGOTimeGPS *injepoch );

int verbose = 1;
int debug = 0;

int lalDebugLevel = 7;
int main( void )
{
	const REAL8 recdur = 180.0;
	/* const REAL8 recdur = 16.0; */
	const REAL8 segdur = 1.0;
	const REAL8 intdur = 0.1;
	LIGOTimeGPS      epoch = {  854378484, 000000000 };
	LIGOTimeGPS      injepoch = {  854378486, 142857143 };
        FrStream        *stream_1;
        FrStream        *stream_2;
        REAL8TimeSeries *raw_series_1;
        REAL8TimeSeries *raw_series_2;
        REAL8TimeSeries *conditioned_series_1;
        REAL8TimeSeries *conditioned_series_2;
	REAL8FFTPlan    *fwdplan;
	REAL8FFTPlan    *revplan;
	COMPLEX8FrequencySeries *response_1 = NULL;
	COMPLEX8FrequencySeries *response_2 = NULL;
	REAL8 ccmax;
	UINT4 intlen;
	UINT4 seglen;

        XLALSetErrorHandler( XLALAbortErrorHandler );

	/* read in the two data streams */
        stream_1 = XLALFrOpen( NULL, "H-H1_RDS_*.gwf" );
        stream_2 = XLALFrOpen( NULL, "H-H2_RDS_*.gwf" );
	verbose ? fprintf( stderr, "reading raw_series_1\n" ) : 0;
	raw_series_1 = XLALFrInputREAL8TimeSeries( stream_1, CHANNEL_1, &epoch, recdur, 0 );
	raw_series_1->sampleUnits = lalStrainUnit;
	verbose ? fprintf( stderr, "reading raw_series_2\n" ) : 0;
	raw_series_2 = XLALFrInputREAL8TimeSeries( stream_2, CHANNEL_2, &epoch, recdur, 0 );
	raw_series_2->sampleUnits = lalStrainUnit;
	XLALFrClose( stream_2 );
	XLALFrClose( stream_1 );
	if ( debug ) {
		verbose ? fprintf( stderr, "writing raw_series_1\n" ) : 0;
		LALDPrintTimeSeries( raw_series_1, "raw_series_1.dat" );
		verbose ? fprintf( stderr, "writing raw_series_2\n" ) : 0;
		LALDPrintTimeSeries( raw_series_2, "raw_series_2.dat" );
	}


	seglen = floor( 0.5 + segdur / raw_series_1->deltaT );
	intlen = floor( 0.5 + intdur / raw_series_1->deltaT );

	if ( 0 != strcmp( CHANNEL_TYPE, "LSC-STRAIN" ) ) {
		response_1 = getresp( "H-H1_CAL_S5_V2-815268733-999999999.gwf", CHANNEL_1, &epoch, 1.0/segdur, seglen );
		response_2 = getresp( "H-H2_CAL_S5_V2-815100013-999999999.gwf", CHANNEL_2, &epoch, 1.0/segdur, seglen );
	}

	/* inject a sine gaussian */
	if ( 1 )
		inject( raw_series_1, raw_series_2, response_1, response_2, &injepoch );

	/* condition the two data streams */
	fwdplan = XLALCreateForwardREAL8FFTPlan( seglen, 0 );
	revplan = XLALCreateReverseREAL8FFTPlan( seglen, 0 );

	verbose ? fprintf( stderr, "conditioning raw_series_1\n" ) : 0;
	conditioned_series_1 = XLALLeonorWhitenREAL8TimeSeries( raw_series_1, fwdplan, revplan, seglen, response_1 );
	verbose ? fprintf( stderr, "conditioning raw_series_2\n" ) : 0;
	conditioned_series_2 = XLALLeonorWhitenREAL8TimeSeries( raw_series_2, fwdplan, revplan, seglen, response_2 );
	if ( debug ) {
		verbose ? fprintf( stderr, "writing conditioned_series_1\n" ) : 0;
		LALDPrintTimeSeries( conditioned_series_1, "conditioned_series_1.dat" );
		verbose ? fprintf( stderr, "writing conditioned_series_2\n" ) : 0;
		LALDPrintTimeSeries( conditioned_series_2, "conditioned_series_2.dat" );
	}

	XLALDestroyCOMPLEX8FrequencySeries( response_2 );
	XLALDestroyCOMPLEX8FrequencySeries( response_1 );
	XLALDestroyREAL8FFTPlan( revplan );
	XLALDestroyREAL8FFTPlan( fwdplan );
	XLALDestroyREAL8TimeSeries( raw_series_2 );
	XLALDestroyREAL8TimeSeries( raw_series_1 );

	ccmax = XLALPearsonMaxCrossCorrelationREAL8Vector( conditioned_series_1->data, conditioned_series_2->data, intlen );
	fprintf( stdout, "ccmax = %g\n", ccmax );

	XLALDestroyREAL8TimeSeries( conditioned_series_2 );
	XLALDestroyREAL8TimeSeries( conditioned_series_1 );

	LALCheckMemoryLeaks();

	return 0;
}

/* get calibration */
COMPLEX8FrequencySeries * getresp( const char *fname, const char *channel, LIGOTimeGPS *epoch, REAL8 deltaF, UINT4 seglen ) 
{
	COMPLEX8FrequencySeries *response;
	LALCalData *caldata;
	verbose ? fprintf( stderr, "getting calibration\n" ) : 0;
	response = XLALCreateCOMPLEX8FrequencySeries( "response", epoch, 0.0, deltaF, &lalDimensionlessUnit, seglen/2 + 1 );
	caldata = XLALFrGetCalData( epoch, channel, fname );
	/* XLALUpdateResponse( response, recdur, caldata ); */
	XLALUpdateResponse( response, 0.0, caldata );
	XLALDestroyCalData( caldata );
	return response;
}

/* make a sine gaussian */
int inject( REAL8TimeSeries *raw_series_1, REAL8TimeSeries *raw_series_2, COMPLEX8FrequencySeries *response_1, COMPLEX8FrequencySeries *response_2, LIGOTimeGPS *injepoch )
{
	REAL8TimeSeries *hplus;
	REAL8TimeSeries *hcross;
	REAL8TimeSeries *h_1;
	REAL8TimeSeries *h_2;
	REAL8TimeSeries *injection_1;
	REAL8TimeSeries *injection_2;
	SkyPosition position = { 0.0, 0.0, COORDINATESYSTEM_EQUATORIAL };
	LALDetector detector_1 = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
	LALDetector detector_2 = lalCachedDetectors[LAL_LHO_2K_DETECTOR];
	LIGOTimeGPS epoch;

	verbose ? fprintf( stderr, "injecting simulated signal\n" ) : 0;
	epoch = raw_series_1->epoch;

	XLALMakeSineGaussian( &hplus, &hcross, injepoch, raw_series_1->deltaT, 8.9, 150.0, 3e-21, 0.0, 0.0 );
	h_1 = XLALSignalDetectorStrainREAL8TimeSeries( hplus, hcross, &position, 0.0, &detector_1 );
	injection_1 = XLALInjectionREAL8TimeSeries( h_1, &epoch, raw_series_1->deltaT, raw_series_1->data->length, response_1 );
	h_2 = XLALSignalDetectorStrainREAL8TimeSeries( hplus, hcross, &position, 0.0, &detector_2 );
	injection_2 = XLALInjectionREAL8TimeSeries( h_2, &epoch, raw_series_1->deltaT, raw_series_1->data->length, response_2 );
	if ( debug ) {
		verbose ? fprintf( stderr, "writing hplus\n" ) : 0;
		LALDPrintTimeSeries( hplus, "hplus.dat" );
		verbose ? fprintf( stderr, "writing hcross\n" ) : 0;
		LALDPrintTimeSeries( hcross, "hcross.dat" );
		verbose ? fprintf( stderr, "writing h_1\n" ) : 0;
		LALDPrintTimeSeries( h_1, "h_1.dat" );
		verbose ? fprintf( stderr, "writing injection_1\n" ) : 0;
		LALDPrintTimeSeries( injection_1, "injection_1.dat" );
		verbose ? fprintf( stderr, "writing h_2\n" ) : 0;
		LALDPrintTimeSeries( h_2, "h_2.dat" );
		verbose ? fprintf( stderr, "writing injection_2\n" ) : 0;
		LALDPrintTimeSeries( injection_2, "injection_2.dat" );
	}
	XLALAddREAL8TimeSeries( raw_series_1, injection_1 );
	XLALAddREAL8TimeSeries( raw_series_2, injection_2 );
	XLALDestroyREAL8TimeSeries( injection_2 );
	XLALDestroyREAL8TimeSeries( injection_1 );
	XLALDestroyREAL8TimeSeries( h_2 );
	XLALDestroyREAL8TimeSeries( h_1 );
	XLALDestroyREAL8TimeSeries( hcross );
	XLALDestroyREAL8TimeSeries( hplus );
	return 0;
}
