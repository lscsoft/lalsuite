/* ****************************************************************
File: LALIndependentTestDetReponce.c

Authors: Greg Mendell, Malik Rakhmanov Started: May 20 2003.

Code: Independently tests of LALComputeDetAMResponse using
model for F_+ and F_x given in Jaranowski, Krolak, and Schutz gr-qc/9804014.

Original Author: Patrick Brian Cameron.  Started: SURF Project Summer 2002.
Original code: generated test continuous graviational wave signals using
the model given in Jaranowski, Krolak, and Schutz gr-qc/9804014.

*******************************************************************/

/* Changes:
05/14/03 gam: Change vertexLongitudeDegrees and vertexLatitudeDegrees to vertexLongitudeRadians vertexLatitudeRadians when using FrDetector structure.
05/15/03 gam: xArmAzimuthRadians and yArmAzimuthRadians are now measured clockwise from North; gamma is the angle to the bisector of the arms, measured
              counterclockwise from East.
05/20/03 gam: Make code compatible with LAL instead of LALapps.
05/20/03 gam: Last argument in LALGPStoLMST1 changed to type LALMSTUnitsAndAcc.
05/20/03 gam: Add Brian Cameron's LALResponseFunc.c into function GenerateResponseFuncUsingLAL(int argc, char *argv[]);
05/20/03 gam: Last argument of LALComputeDetAMResponse changed to type LALGPSandAcc.
05/23/04 Dave Chin globally changed "time" to "timevec" to fix warnings about shadowing global declarations.
05/30/03 gam: Change GenerateResponseFuncUsingLAL to work with params in the config file.
09/26/03 gam: Clean up GenerateResponseFuncUsingLAL, removed unnecessary code.
09/26/03 gam: Need to set LALTimeIntervalAndNSample time_info.accuracy = LALLEAPSEC_STRICT.
09/29/03 gam: Return data from GenerateResponseFuncUsingLAL
09/30/03 gam: Change name of LALComputeResponseFunc to more descriptive name GenerateResponseFuncNotUsingLAL.
09/30/03 gam: Clean up GenerateResponseFuncNotUsingLAL, removed unnecessary code.
09/30/03 gam: Always call GenerateResponseFuncUsingLAL
09/30/03 gam: Removed all other unnecessary code.
09/30/03 gam: Add code to test difference between independent and LAL values for F_+ and F_x.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/* #include <lalapps.h> */ /* 05/20/03 gam */
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALConfig.h>
#include <lal/AVFactories.h>
#include <lal/DetectorSite.h>
#include <lal/DetResponse.h>
#include <lal/Date.h>
#include <lal/Random.h>

#define usgfmt \
     "Usage: %s [options]\n" \
     " -h                  Prints this message\n" \
     " -c configFile       Use this Configure File\n" \
     " -d lalDebugLevel    Set debug level\n" \
     "\n" \
     "NOTE: This program can only support filenames up to 25 characters\n" \
     "      in length in the config file (this is the string length so includes\n" \
     "      the directories in the name as well), but the configFile entered at\n" \
     "      the command line can be as long as you like.\n" \
     "      To change this you can edit the parameter valueSize in the source\n" \
     "      code, but make sure the none of the comments in the config file start\n" \
     "      before the new value, or extend past the line size of 80 characters.\n"

#define usage( program ) fprintf( stderr, usgfmt, program )

NRCSID( LALINDEPENDENTTESTDETRESPONSEC, "$Id$" );

void GenerateResponseFuncNotUsingLAL(LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet, REAL8 phiStart, REAL8Vector *timevec, REAL8Vector *fPlus, REAL8Vector *fCross);

void PrintLALDetector(LALDetector *detector);

void GenerateResponseFuncUsingLAL(LALStatus *status, LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet, REAL8 sampleRate, LIGOTimeGPS *gps, LALDetAMResponseSeries  *am_response_series_ptr);

extern char *optarg;
REAL8 omegaEarth = LAL_TWOPI/LAL_DAYSID_SI;
REAL8 omegaEarthSun = LAL_TWOPI/LAL_YRSID_SI;

INT4 lalDebugLevel = 0;

int main( int argc, char *argv[] )
{
  /* Parsing Command Line, and General Use */
  const CHAR *program = argv[0];
  static LALStatus status;
  INT4 opt;
  UINT4 i;

  /* For reading Config File */
  UINT4 lineSize = 80, valueSize  = 25;
  CHARVector *lineString = NULL, *valueString = NULL;
  INT4 checkFileName = 0;
  /* const CHAR *configFileName = NULL; */ /* 09/23/03 gam */
  CHAR *configFileName = NULL;
  FILE *configFile = NULL;

  /* Flags */ /* 09/30/03 gam */
  INT2 outputIndependentFPlusFCross;
  INT2 outputLALFPlusFCross;
  INT2 outputFPlusFCrossDiffs;

  /* For computing differences */ /* 09/30/03 gam */
  REAL4 maxFPlusDiff = 0.0;
  REAL4 maxFCrossDiff = 0.0;
  REAL4 tmpFPlusDiff = 0.0;
  REAL4 tmpFCrossDiff = 0.0;
  REAL4 tolerance = 0.0;
  /* REAL4 tmpFPlusDenom = 0.0; */ /* Only used if testing fractional difference. Not currently used. */
  /* REAL4 tmpFCrossDenom = 0.0; */ /* Only used if testing fractional difference. Not currently used. */
  /* REAL4 TINY = 1.e-4; */ /* Only used if testing fractional difference. Not currently used. */

  /* For File Output */
  FILE *outFileFPlus = NULL, *outFileFCross = NULL;
  FILE *outFileLALPlus = NULL, *outFileLALCross = NULL;
  FILE *outFilePlusDiff = NULL, *outFileCrossDiff = NULL;

  /* Used LAL structures */
  LIGOTimeGPS gpsTime;
  LALSource pulsar;
  LALDetector detector;
  LALPlaceAndGPS place_and_gps;
  LALDetAMResponseSeries am_response_series; /* 09/30/03 gam */
  REAL4TimeSeries plus_series, cross_series, scalar_series; /* 09/30/03 gam */

  /* first num is longitude, second is latitude (radians) 05/14/03 gam */
  LALFrDetector testFrDetector;

  /* Vectors For Results */
  REAL8Vector *fPlus=NULL, *fCross=NULL;

  /* Other Stuff */
  REAL8Vector *timevec = NULL;
  REAL8 sampleRate = 0, duration = 0;
  UINT4 lgthDataSet = 0;
  REAL8 phiStart = 0.0;

  /* LALMSTUnitsAndAcc uandacc = { MST_RAD, LALLEAPSEC_STRICT}; */ /* 05/20/03 gam */
  LALMSTUnitsAndAcc uandacc; /* 05/20/03 gam */

  uandacc.units = MST_RAD;
  uandacc.accuracy = LALLEAPSEC_STRICT;

  /* parse options */
  if( argc == 1)
    {
      usage( program );
      return 0;
    }
  /* while ( 0 < ( opt = getopt( argc, argv, "hc:d:" ) ) ) */ /* Added case L */ /* 05/20/03 gam */
  while ( 0 < ( opt = getopt( argc, argv, "hc:Ld:" ) ) )
    {
      switch ( opt )
        {
        case 'h':
          usage( program );
          return 0;
        case 'c':
	  configFileName = optarg;
          break;
        case 'L':
	  break;
	case 'd':
	  lalDebugLevel = atoi( optarg );
	  break;
	default:
          usage( program );
          return 1;
        }
    }

  /* set_debug_level( "NDEBUG" ); */ /* 05/20/03 gam */

  /*
   * Read in Configuration File
   */
  configFile = fopen( configFileName, "r" );
  if(configFile == NULL){
    fprintf( stderr, "Error Opening ConfigFile\n" );
    usage( program );
    return 0;
  }

  LALCHARCreateVector(&status, &lineString, lineSize);
  LALCHARCreateVector(&status, &valueString, valueSize);

  /* Throw Away First Line Comment */
  fgets(lineString->data, lineSize, configFile);

  /* Test Name */
  strcpy(pulsar.name, "TEST");
  pulsar.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;

  /* R.A. of the Source */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  pulsar.equatorialCoords.longitude = LAL_PI_180*(atof( valueString->data ));

  /* Declination of the Source */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  pulsar.equatorialCoords.latitude = LAL_PI_180*(atof( valueString->data ));

  /* Polarization Angle */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  pulsar.orientation = LAL_PI_180*atof( valueString->data );

  /* Start Time */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  /* gpsTime.gpsSeconds = atoi( valueString->data ); */
  gpsTime.gpsSeconds = atol( valueString->data ); /* 09/30/03 gam */
  gpsTime.gpsNanoSeconds = 0;

  /* Sample Rate */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  sampleRate = atof( valueString->data );

  /* Duration */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  duration = atof( valueString->data );

  /* lgthDataSet = (UINT4)(duration * sampleRate) + 1; */ /* 09/30/03 gam */
  lgthDataSet = (UINT4)(duration * sampleRate);

  /* Detector Site */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  if(valueString->data[0] == 'H')
    detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  if(valueString->data[0] == 'L')
    detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  if(valueString->data[0] == 'D')
    {
      fprintf( stderr,  "Enter Properties of test detector\n" );
      /* fprintf( stderr, "Detector Longitude (deg): " );
      scanf( "%lf", &testFrDetector.vertexLongitudeDegrees ); */ /* 05/14/03 gam */
      fprintf( stderr, "Detector Longitude (radians): " );
      scanf( "%lf", &testFrDetector.vertexLongitudeRadians );

      /* fprintf( stderr, "Detector Latitude (deg): " );
      scanf( "%lf" , &testFrDetector.vertexLatitudeDegrees ); */ /* 05/14/03 gam */
      fprintf( stderr, "Detector Latitude (radians): " );
      scanf( "%lf" , &testFrDetector.vertexLatitudeRadians );

      /* fprintf( stderr, "Detector xArmAzimuth (deg): " ); */ /* 05/14/03 gam */
      fprintf( stderr, "Detector xArmAzimuth (radians): " );
      scanf( "%f", &testFrDetector.xArmAzimuthRadians );

      /* testFrDetector.xArmAzimuthRadians *= LAL_PI_180; */ /* 05/14/03 gam */
      /* testFrDetector.yArmAzimuthRadians = testFrDetector.xArmAzimuthRadians + LAL_PI_2; */ /* 05/15/03 gam */
      testFrDetector.yArmAzimuthRadians = testFrDetector.xArmAzimuthRadians - LAL_PI_2;
      testFrDetector.vertexElevation = 0.0;
      testFrDetector.xArmAltitudeRadians = 0.0;
      testFrDetector.yArmAltitudeRadians = 0.0;
      LALCreateDetector( &status, &detector, &testFrDetector, LALDETECTORTYPE_IFODIFF);
    }

  /* tolerance */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  tolerance = atof( valueString->data );

  /* Open Various Output Files */

  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  outputIndependentFPlusFCross = atoi( valueString->data );

  if (outputIndependentFPlusFCross) {
    /*  F_Plus Filename */
    fgets(lineString->data, lineSize, configFile);
    strncpy(valueString->data, lineString->data, valueSize);
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    outFileFPlus = fopen(valueString->data, "w");

    /* F_Cross Filename */
    fgets(lineString->data, lineSize, configFile);
    strncpy(valueString->data, lineString->data, valueSize);
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    outFileFCross = fopen(valueString->data, "w");
 }

  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  outputLALFPlusFCross = atoi( valueString->data );

  if (outputLALFPlusFCross) {
    /*  F_Plus LAL Filename */
    fgets(lineString->data, lineSize, configFile);
    strncpy(valueString->data, lineString->data, valueSize);
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    outFileLALPlus = fopen(valueString->data, "w");

    /* F_Cross LAL Filename */
    fgets(lineString->data, lineSize, configFile);
    strncpy(valueString->data, lineString->data, valueSize);
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    outFileLALCross = fopen(valueString->data, "w");
 }

  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  outputFPlusFCrossDiffs = atoi( valueString->data );

  if (outputFPlusFCrossDiffs) {

    /*  F_Plus LAL Filename */
    fgets(lineString->data, lineSize, configFile);
    strncpy(valueString->data, lineString->data, valueSize);
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    outFilePlusDiff = fopen(valueString->data, "w");

    /* F_Cross LAL Filename */
    fgets(lineString->data, lineSize, configFile);
    strncpy(valueString->data, lineString->data, valueSize);
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    outFileCrossDiff = fopen(valueString->data, "w");

 }
 fclose(configFile);

  /* Echo Parameters */
  if( lalDebugLevel > 1)
    {
      fprintf( stderr,  "\n" );
      fprintf( stderr, "Source Info\n" );
      fprintf( stderr, "RA: %e\n", pulsar.equatorialCoords.longitude );
      fprintf( stderr, "Dec: %e\n", pulsar.equatorialCoords.latitude );
      fprintf( stderr, "Polarization: %f\n", pulsar.orientation );
      fprintf( stderr, "\n" );

      fprintf( stderr, "Date Info\n" );
      fprintf( stderr, "For GPS Seconds, NanoSeconds: %10d, %9d\n", gpsTime.gpsSeconds, gpsTime.gpsNanoSeconds );
      fprintf( stderr, "\n" );

      fprintf( stderr, "Sampling Rate: %f\n", sampleRate );
      fprintf( stderr, "Length of Data Set: %f\n", duration );
      fprintf( stderr, "\n");

      PrintLALDetector( &detector );

      fprintf( stderr, "Ang. Vel. of Earth = %e\n", omegaEarth );
      fprintf( stderr, "Ang. Vel. of Earth about the Sun = %e\n", omegaEarthSun );
      fprintf( stderr, "\n");
    }

  /*
   * Start Computing Signal
   */
  place_and_gps.p_detector = &detector;
  place_and_gps.p_gps = &gpsTime;

  /* Compute Local Sidereal Time at the start of the data set*/
  /* LALGPStoLMST1(&status, &phiStart,  &place_and_gps, MST_RAD); */ /* 05/20/03 gam */
  LALGPStoLMST1(&status, &phiStart,  &place_and_gps, &uandacc);

  if(lalDebugLevel > 1) fprintf( stderr, "LMST (radians)= %f\n", phiStart );

  /* Intialize Arrays */
  if(lalDebugLevel > 0) fprintf( stderr, "Creating Vectors.\n" );

  LALDCreateVector(&status, &fPlus, lgthDataSet);
  LALDCreateVector(&status, &fCross, lgthDataSet);
  LALDCreateVector(&status, &timevec, lgthDataSet);

  plus_series.data = NULL;
  cross_series.data = NULL;
  scalar_series.data = NULL;
  am_response_series.pPlus   = &(plus_series);
  am_response_series.pCross  = &(cross_series);
  am_response_series.pScalar = &(scalar_series);

  LALSCreateVector(&status, &(am_response_series.pPlus->data), lgthDataSet);
  LALSCreateVector(&status, &(am_response_series.pCross->data), lgthDataSet);
  LALSCreateVector(&status, &(am_response_series.pScalar->data), lgthDataSet);

  /* make time array */
  for(i = 0;i<lgthDataSet;i++)
    {
      timevec->data[i] = i*(1.0/sampleRate);
    }

  if(lalDebugLevel > 0) {
    fprintf(stdout, "Computing Independent Response Functions.\n" );
    fflush(stdout);
  }

  GenerateResponseFuncNotUsingLAL(&pulsar, &detector, lgthDataSet, phiStart, timevec, fPlus, fCross );

  if(lalDebugLevel > 0) {
    fprintf(stdout, "Computing Response Functions Using LAL.\n" );
    fflush(stdout);
  }

  GenerateResponseFuncUsingLAL(&status, &pulsar, &detector, lgthDataSet, sampleRate, &gpsTime, &am_response_series);

  /* Write out files */
  if (outputIndependentFPlusFCross) {
    if(lalDebugLevel > 0) fprintf( stderr, "Writing Independent F_+ and F_x to File.\n" );
    for (i = 0;i<lgthDataSet; i++) {
        fprintf(outFileFPlus, "%.15f\n", fPlus->data[i]);
        fprintf(outFileFCross, "%.15f\n", fCross->data[i]);
    }
    fclose(outFileFPlus);
    fclose(outFileFCross);
  }
  if (outputLALFPlusFCross) {
    if(lalDebugLevel > 0) fprintf( stderr, "Writing LAL F_+ and F_x to File.\n" );
    for (i = 0; i < lgthDataSet; ++i) {
        fprintf(outFileLALPlus, "%.15f\n", am_response_series.pPlus->data->data[i]);
        fprintf(outFileLALCross, "%.15f\n", am_response_series.pCross->data->data[i]);
    }
    fclose(outFileLALPlus);
    fclose(outFileLALCross);

  }

  /* Compute differences; write out differences if requested */ /* 09/30/03 gam */
  if(lalDebugLevel > 0) fprintf( stderr, "Computing differences LAL F_+ - Independent F_+ and LAL F_x - Independent F_x.\n" );
  if (outputFPlusFCrossDiffs && lalDebugLevel > 0) fprintf( stderr, "Writing differences LAL F_+ - Independent F_+ and LAL F_x - Independent F_x to File.\n" );
  maxFPlusDiff = 0.0;  /* initialize */
  maxFCrossDiff = 0.0; /* initialize */
  for (i = 0; i < lgthDataSet; ++i) {
    tmpFPlusDiff = am_response_series.pPlus->data->data[i] - (REAL4)fPlus->data[i];
    tmpFCrossDiff = am_response_series.pCross->data->data[i] - (REAL4)fCross->data[i];
    /* tmpFPlusDenom and tmpFCrossDenom used for testing fractional difference. Not currently used. */
    /* tmpFPlusDenom = fabs((REAL4)fPlus->data[i]); */
    /* tmpFCrossDenom = fabs((REAL4)fCross->data[i]); */
    if (outputFPlusFCrossDiffs) {
      fprintf(outFilePlusDiff, "%.15f\n",  tmpFPlusDiff);
      fprintf(outFileCrossDiff, "%.15f\n", tmpFCrossDiff);
    }
    /* Code for testing absolute difference */
    tmpFPlusDiff = fabs(tmpFPlusDiff);
    if (tmpFPlusDiff > maxFPlusDiff) maxFPlusDiff = tmpFPlusDiff;
    tmpFCrossDiff = fabs(tmpFCrossDiff);
    if (tmpFCrossDiff > maxFCrossDiff) maxFCrossDiff = tmpFCrossDiff;
    /* Code for testing fractional difference. Not currently used. */
    /* if (tmpFPlusDenom > TINY) {
      tmpFPlusDiff = fabs(tmpFPlusDiff)/tmpFPlusDenom;
      if (tmpFPlusDiff > maxFPlusDiff) maxFPlusDiff = tmpFPlusDiff;
    }
    if (tmpFCrossDenom >  TINY) {
      tmpFCrossDiff = fabs(tmpFCrossDiff)/tmpFCrossDenom;
      if (tmpFCrossDiff > maxFCrossDiff) maxFCrossDiff = tmpFCrossDiff;
    } */
  }
  if (outputFPlusFCrossDiffs) {
    fclose(outFilePlusDiff);
    fclose(outFileCrossDiff);
  }
  if(lalDebugLevel > 0) {
    fprintf(stdout, "\nmaxFPlusDiff, tolerance = %.6e, %.6e\n",maxFPlusDiff, tolerance );
    fprintf(stdout, "\nmaxFCrossDiff, tolerance = %.6e, %.6e\n",maxFCrossDiff, tolerance);
    fflush(stdout);
  }

  /* Free memory, and then test difference */
  LALSDestroyVector(&status, &(am_response_series.pPlus->data));
  LALSDestroyVector(&status, &(am_response_series.pCross->data));
  LALSDestroyVector(&status, &(am_response_series.pScalar->data));

  LALCHARDestroyVector( &status, &lineString );
  LALCHARDestroyVector( &status, &valueString );
  LALDDestroyVector( &status, &fPlus );
  LALDDestroyVector( &status, &fCross );
  LALDDestroyVector( &status, &timevec );

  /* Test differences. End with error if difference exceed tolerance. */
  if (maxFPlusDiff > tolerance || maxFCrossDiff > tolerance) {
    if (maxFPlusDiff > tolerance) {
      fprintf(stderr, "\nLAL F_+ - Independent F_+ Difference Test Failed!\n");
      fflush(stderr);
    }
    if (maxFCrossDiff > tolerance) {
      fprintf(stderr, "\nLAL F_x - Independent F_x Difference Test Failed!\n");
      fflush(stderr);
    }
    exit(1);
  }

  if(lalDebugLevel > 0) {
    fprintf(stdout, "\nALL LALIndependentTestDetResponse Test Passed!\n");
    fflush(stdout);
  }

  return 0;
}

void GenerateResponseFuncNotUsingLAL(LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet,
			    REAL8 phiStart, REAL8Vector *timevec,  REAL8Vector *fPlus, REAL8Vector *fCross)
{
  REAL8 aCoeff, bCoeff;
  /* REAL8 gammaAngle = detector->frDetector.xArmAzimuthRadians + LAL_PI_4; */ /* 05/15/03 gam */
  REAL8 gammaAngle = LAL_PI_2 + LAL_PI_4 - detector->frDetector.xArmAzimuthRadians;
  /* REAL8 detLatitude = LAL_PI_180*detector->frDetector.vertexLatitudeDegrees; */ /* 05/14/03 gam */
  REAL8 detLatitude = detector->frDetector.vertexLatitudeRadians;
  INT4 i;

  /* Compute Fplus and Fcross from JKS */ /* 09/30/03 gam aCoeff and bCoeff no longer vectors */
  for(i = 0;i<lgthDataSet; i++) {
      aCoeff = sin(2.0*gammaAngle)*(3.0 - cos(2.0*detLatitude))*(3.0 - cos(2.0*pulsar->equatorialCoords.latitude))
	*cos(2.0*(pulsar->equatorialCoords.longitude - phiStart - omegaEarth*timevec->data[i]))/16.0;
      aCoeff -= 0.25*cos(2.0*gammaAngle)*sin(detLatitude)*(3 - cos(2.0*pulsar->equatorialCoords.latitude))
	*sin(2.0*(pulsar->equatorialCoords.longitude - phiStart - omegaEarth*timevec->data[i]));
      aCoeff += 0.25*sin(2.0*gammaAngle)*sin(2.0*detLatitude)*sin(2.0*pulsar->equatorialCoords.latitude)
	*cos(pulsar->equatorialCoords.longitude - phiStart - omegaEarth*timevec->data[i]);
      aCoeff -= 0.5*cos(2.0*gammaAngle)*cos(detLatitude)*sin(2.0*pulsar->equatorialCoords.latitude)
	*sin(pulsar->equatorialCoords.longitude - phiStart - omegaEarth*timevec->data[i]);
      aCoeff += 0.75*sin(2.0*gammaAngle)*cos(detLatitude)*cos(detLatitude)
	*cos(pulsar->equatorialCoords.latitude)*cos(pulsar->equatorialCoords.latitude);

      bCoeff = cos(2*gammaAngle)*sin(detLatitude)*sin(pulsar->equatorialCoords.latitude)
	*cos(2.0*(pulsar->equatorialCoords.longitude-phiStart-omegaEarth*timevec->data[i]));
      bCoeff += 0.25*sin(2.0*gammaAngle)*(3.0-cos(2.0*detLatitude))*sin(pulsar->equatorialCoords.latitude)
	*sin(2.0*(pulsar->equatorialCoords.longitude-phiStart-omegaEarth*timevec->data[i]));
      bCoeff += cos(2*gammaAngle)*cos(detLatitude)*cos(pulsar->equatorialCoords.latitude)
	*cos(pulsar->equatorialCoords.longitude-phiStart-omegaEarth*timevec->data[i]);
      bCoeff += 0.5*sin(2.0*gammaAngle)*sin(2.0*detLatitude)*cos(pulsar->equatorialCoords.latitude)
	*sin(pulsar->equatorialCoords.longitude-phiStart-omegaEarth*timevec->data[i]);
      /* Note: Assumes Interferometer Arms are Perpendicular */
      fPlus->data[i] = aCoeff*cos(2.0*pulsar->orientation) + bCoeff*sin(2.0*pulsar->orientation);
      fCross->data[i] = bCoeff*cos(2.0*pulsar->orientation) - aCoeff*sin(2.0*pulsar->orientation);
  }
  return;
}

void PrintLALDetector(LALDetector *detector)
{
  /* REAL8 gamma = detector->frDetector.xArmAzimuthRadians + LAL_PI_4; */ /* 05/15/03 gam */
  REAL8 gamma = LAL_PI_2 + LAL_PI_4 - detector->frDetector.xArmAzimuthRadians;

  fprintf( stderr,  "Location  = { { %.6e, %.6e, %.6e },\n", detector->location[0],
	   detector->location[1], detector->location[2] );
  fprintf( stderr, " \"%s\",\n", detector->frDetector.name);
  /* fprintf( stderr, "Location =   %.6g, %.6g, %.6g,\n",
	   detector->frDetector.vertexLongitudeDegrees,
	   detector->frDetector.vertexLatitudeDegrees,
	   detector->frDetector.vertexElevation); */ /* 05/14/03 gam */
  fprintf( stderr, "Location =   %.6g, %.6g, %.6g,\n",
	   detector->frDetector.vertexLongitudeRadians,
	   detector->frDetector.vertexLatitudeRadians,
	   detector->frDetector.vertexElevation);
  fprintf( stderr, "Angle to Arm Bisector %f\n", gamma );
  fprintf( stderr, "\n" );

  return;
}

void GenerateResponseFuncUsingLAL(LALStatus *status, LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet, REAL8 sampleRate, LIGOTimeGPS *gps, LALDetAMResponseSeries  *am_response_series_ptr)
{
  LALDetAndSource   det_and_pulsar;
  LALDetAMResponse  am_response;
  LALPlaceAndGPS    place_and_gps; REAL8 lmsthours;

  LALTimeIntervalAndNSample time_info;

  /* LALMSTUnitsAndAcc uandacc = { MST_RAD, LALLEAPSEC_STRICT}; */ /* 05/20/03 gam */
  LALMSTUnitsAndAcc uandacc; /* 05/20/03 gam */
  LALGPSandAcc gps_and_acc; /* 05/20/03 gam */

  INITSTATUS (status, "GenerateResponseFuncUsingLAL", LALINDEPENDENTTESTDETRESPONSEC);
  ATTATCHSTATUSPTR(status);

  uandacc.units = MST_RAD;
  uandacc.accuracy = LALLEAPSEC_STRICT;

  if (lalDebugLevel > 0)  {
     fprintf(stdout,"status->statusCode = %i \n", status->statusCode);
     fprintf(stdout,"status->statusPtr->statusCode = %i \n", status->statusPtr->statusCode);
     printf("Source Info\n");
     printf("RA: %e\n", pulsar->equatorialCoords.longitude);
     printf("Dec: %e\n", pulsar->equatorialCoords.latitude);
     printf("Psi: %e\n", pulsar->orientation);
     printf("Det #1 location: (%7.4e, %7.4e, %7.4e)\n",
	   detector->location[0], detector->location[1],
           detector->location[2]);
     fprintf(stderr, "GPS = {%10d, %9d}\n", gps->gpsSeconds, gps->gpsNanoSeconds);
  }

  /* Add detector and source into structures */
  place_and_gps.p_detector = detector;
  place_and_gps.p_gps = gps;
  det_and_pulsar.pDetector = detector;
  det_and_pulsar.pSource   = pulsar;

  /* LALGPStoLMST1(&status, &lmsthours, &place_and_gps, MST_RAD); */ /* 05/20/03 gam */
  LALGPStoLMST1(status->statusPtr, &lmsthours, &place_and_gps, &uandacc);
  CHECKSTATUSPTR (status);

  if (lalDebugLevel > 0)  {
    fprintf(stdout, "LMST = %7e \n", lmsthours);
  }

  /* Compute instantaneous AM response */
  gps_and_acc.gps.gpsSeconds     = gps->gpsSeconds;
  gps_and_acc.gps.gpsNanoSeconds = gps->gpsNanoSeconds;
  gps_and_acc.accuracy = LALLEAPSEC_STRICT; /* 05/20/03 gam */

  /* LALComputeDetAMResponse(&status, &am_response, &det_and_pulsar, &gps); */ /* 05/20/03 gam */
  LALComputeDetAMResponse(status->statusPtr, &am_response, &det_and_pulsar, &gps_and_acc);
  CHECKSTATUSPTR (status);

   if (lalDebugLevel > 0) {
   }

  time_info.epoch.gpsSeconds     = gps->gpsSeconds;
  time_info.epoch.gpsNanoSeconds = gps->gpsNanoSeconds;
  time_info.deltaT               = 1.0/sampleRate;
  time_info.nSample              = lgthDataSet;
  time_info.accuracy = LALLEAPSEC_STRICT; /* 09/26/03 gam */

   if (lalDebugLevel > 0) {
      printf("Start computing AM response vectors\n");
      /* fprintf(stdout, "am_response_series_ptr->pPlus->data->length = %d\n", am_response_series_ptr->pPlus->data->length);
      fprintf(stdout, "am_response_series_ptr->pCross->data->length = %d\n", am_response_series_ptr->pCross->data->length);
      fprintf(stdout, "am_response_series_ptr->pScalar->data->length = %d\n", am_response_series_ptr->pScalar->data->length); */
      fprintf(stdout, "timeinfo = %d %d %g %d\n",time_info.epoch.gpsSeconds,time_info.epoch.gpsNanoSeconds,time_info.deltaT,time_info.nSample);
      fflush(stdout);
   }

  LALComputeDetAMResponseSeries(status->statusPtr, am_response_series_ptr, &det_and_pulsar, &time_info);
  CHECKSTATUSPTR (status);

  if (lalDebugLevel > 0) {
      printf("Done computing AM response vectors\n");
      /* printf("am_response_series_ptr->pPlus->data->length = %d\n", am_response_series_ptr->pPlus->data->length);
      printf("am_response_series_ptr->pCross->data->length = %d\n", am_response_series_ptr->pCross->data->length);
      printf("am_response_series_ptr->pScalar->data->length = %d\n", am_response_series_ptr->pScalar->data->length); */
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
