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
     " -L                  Generate the response using LAL only and return \n" \
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

NRCSID( GENERATECWSIGNALC, "$Id$" );

static void
LALComputeResponseFunc(LALStatus *status, LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet,
		       REAL8 phiStart, REAL8Vector *time, REAL8Vector *fPlus, REAL8Vector *fCross);
static void
LALComputePolarFuncs(LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet, REAL8 phiStart, REAL8Vector *time,
		     REAL8Vector *freqParams, REAL8Vector *hPlus, REAL8Vector *hCross);
static void
PrintLALDetector(LALDetector *detector);

int GenerateResponseFuncUsingLAL(int argc, char *argv[]); /* 05/20/03 gam */

extern char *optarg;
REAL8 omegaEarth = LAL_TWOPI/LAL_DAYSID_SI;
REAL8 omegaEarthSun = LAL_TWOPI/LAL_YRSID_SI;
REAL8 phiSun;
REAL8 plusSignalAmplitude;
REAL8 crossSignalAmplitude;


int main( int argc, char *argv[] )
{
  /* Parsing Command Line, and General Use */
  const CHAR *program = argv[0];
  static LALStatus status;
  INT4 opt;
  UINT4 i;
  INT4 lalDebugLevel = 0;

  /* For reading Config File */
  UINT4 lineSize = 80, valueSize  = 25;
  CHARVector *lineString = NULL, *valueString = NULL;
  INT4 checkFileName = 0;
  const CHAR *configFileName = NULL;
  FILE *configFile = NULL;

  /* For Printing and File Output */
  FILE *outFileFPlus = NULL, *outFileFCross = NULL;
  FILE *outFileHPlus = NULL, *outFileHCross = NULL;
  FILE *outFileNoise = NULL, *outFileSignal = NULL;
  FILE *weightFile = NULL;

  /* Used LAL structures */
  LIGOTimeGPS gpsTime;
  LALSource pulsar;
  LALDetector detector;
  LALPlaceAndGPS place_and_gps;

  /* first num is longitude, second is latitude (radians) 05/14/03 gam */
  LALFrDetector testFrDetector;

  /* Vectors For Results */
  static REAL4Vector  *normalDeviates = NULL;
  static RandomParams *randpar = NULL;
  REAL8Vector *fPlus=NULL, *fCross=NULL;
  REAL8Vector *hPlus=NULL, *hCross=NULL;
  REAL8Vector *cwSignal=NULL;
  REAL4Vector *weightVector = NULL;

  /* Other Stuff */
  REAL8Vector *time = NULL;
  REAL8Vector *frequencyParams = NULL;
  REAL8 sampleRate = 0, duration = 0;
  UINT4 numFreqParams = 6;
  UINT4 lgthDataSet = 0;
  REAL8 phiStart = 0.0;
  REAL8 randSigma = 0.0;

 LALMSTUnitsAndAcc uandacc = { MST_RAD, LALLEAPSEC_STRICT}; /* 05/20/03 gam */

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
          GenerateResponseFuncUsingLAL(argc, argv); /* Added case L */ /* 05/20/03 gam */
          return 0;
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
  if(lalDebugLevel > 2)
    REPORTSTATUS (&status);
  LALCHARCreateVector(&status, &valueString, valueSize);
  if(lalDebugLevel > 2)
    REPORTSTATUS (&status);

  /* Throw Away First Line, Formatting Comment */
  fgets(lineString->data, lineSize, configFile);


  /* Start by Reading Pulsar Coordinates */
  strcpy(pulsar.name, "INPUT PULSAR");
  pulsar.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;

  /* R.A. of the Source */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  pulsar.equatorialCoords.longitude = LAL_PI_180*(atof( valueString->data ));

  /* Declination of the Source */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  pulsar.equatorialCoords.latitude = LAL_PI_180*(atof( valueString->data ));

  /* Pulsar Polarization Angle */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  pulsar.orientation = LAL_PI_180*atof( valueString->data );

  /* Parse Date */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  gpsTime.gpsSeconds = atoi( valueString->data );
  gpsTime.gpsNanoSeconds = 0;

  /* Read in Sample Rate */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  sampleRate = atof( valueString->data );

  /* Length of the Data Set */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  duration = atof( valueString->data );

  lgthDataSet = (UINT4)(duration * sampleRate) + 1;

  /* Read in Frequency Paramters */
  LALDCreateVector(&status, &frequencyParams, numFreqParams);
  if(lalDebugLevel > 2)
    REPORTSTATUS (&status);

  for(i=0;i<numFreqParams;i++)
    {
      fgets(lineString->data, lineSize, configFile);
      strncpy(valueString->data, lineString->data, valueSize);
      frequencyParams->data[i] = atof( valueString->data );
    }


  /* Read in the Detector Site */
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

  /* Options For Random Number Creation */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  randSigma = atof( valueString->data );

  /*
   * Multiply Earth's spin or orbit by some constant
   * Will allow us to simplify by letting the earth stand still
   * or speed the earth up to facilitate smaller data sets.
   */

  /* Earth's Rotation Axis */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  omegaEarth *= atof( valueString->data );

  /* Rotation about the Sun */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  omegaEarthSun *= atof( valueString->data );

  /* Read in Phase of the earth's orbit */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  phiSun = atof( valueString->data );

  /* Read in Signal Amplitudes */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  plusSignalAmplitude = atof( valueString->data );

  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  crossSignalAmplitude = atof( valueString->data );

  /* Open Various Output Files */

  /* Weighting File */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  if(valueString->data[0] != '0'){
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    weightFile = fopen(valueString->data, "r");
  }

  /*  F_Plus Filename */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  if(valueString->data[0] != '0'){
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    outFileFPlus = fopen(valueString->data, "w");
  }

  /* F_Cross Filename */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  if(valueString->data[0] != '0'){
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

  /* H_Plus Filename */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  if(valueString->data[0] != '0'){
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    outFileHPlus = fopen(valueString->data, "w");
  }

  /* H_Cross Filename */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  if(valueString->data[0] != '0'){
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    outFileHCross = fopen(valueString->data, "w");
  }

  /* Noise Filename */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  if(valueString->data[0] != '0'){
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    outFileNoise = fopen(valueString->data, "w");
  }

  /* Signal Filename */
  fgets(lineString->data, lineSize, configFile);
  strncpy(valueString->data, lineString->data, valueSize);
  if(valueString->data[0] != '0'){
    i=0; checkFileName = 0;
    while( !checkFileName ) {
      if(valueString->data[i] == ' ') {
	valueString->data[i] = '\0';
	checkFileName = 1;
      }
      i++;
    }
    outFileSignal = fopen(valueString->data, "w");
  }


  /* Echo Parameters */
  if( lalDebugLevel > 1)
    {
      fprintf( stderr,  "\n" );
      fprintf( stderr, "Source Info\n" );
      fprintf( stderr, "RA: %e\n", pulsar.equatorialCoords.longitude );
      fprintf( stderr, "Dec: %e\n", pulsar.equatorialCoords.latitude );
      fprintf( stderr, "Polarization: %f\n", pulsar.orientation );
      fprintf( stderr, "\n" );

      fprintf( stderr, "Freq. Params\n" );
      fprintf( stderr, "%e %e %e %e %e %e\n", frequencyParams->data[0], frequencyParams->data[1],
	       frequencyParams->data[2], frequencyParams->data[3],frequencyParams->data[4]
	       ,frequencyParams->data[5] );
      fprintf( stderr, "\n" );

      fprintf( stderr, "Date Info\n" );
      fprintf( stderr, "For GPS: %10d, %9d\n", gpsTime.gpsSeconds, gpsTime.gpsNanoSeconds );
      fprintf( stderr, "\n" );

      fprintf( stderr, "Sampling Rate: %f\n", sampleRate );
      fprintf( stderr, "Length of Data Set: %f\n", duration );
      fprintf( stderr, "\n");

      PrintLALDetector( &detector );

      fprintf( stderr, "Ang. Vel. of Earth = %e\n", omegaEarth );
      fprintf( stderr, "Ang. Vel. of Earth about the Sun = %e\n", omegaEarthSun );
      fprintf( stderr, "\n");

      fprintf( stderr, "Plus Signal Amplitude %f\n", plusSignalAmplitude );
      fprintf( stderr, "Cross Signal Amplitude %f\n", crossSignalAmplitude );
      fprintf( stderr, "\n" );

      fprintf( stderr, "Standard Deviation of Noise %f\n", randSigma );
      fprintf( stderr, "\n" );
    }

  /*
   * Start Computing Signal
   */
  place_and_gps.p_detector = &detector;
  place_and_gps.p_gps = &gpsTime;

  /* Compute Local Sidereal Time at the start of the data set*/
  /* LALGPStoLMST1(&status, &phiStart,  &place_and_gps, MST_RAD); */ /* 05/20/03 gam */
  LALGPStoLMST1(&status, &phiStart,  &place_and_gps, &uandacc);

  if(lalDebugLevel > 1)
    fprintf( stderr, "LMST (radians)= %f\n", phiStart );
  

  /* Intialize Arrays */
  if(lalDebugLevel > 0)
    fprintf( stderr, "Creating Vectors.\n" );

  LALDCreateVector(&status, &fPlus, lgthDataSet);
  if(lalDebugLevel > 2)
    REPORTSTATUS (&status);
  LALDCreateVector(&status, &fCross, lgthDataSet);
  if(lalDebugLevel > 2)
    REPORTSTATUS (&status);
  LALDCreateVector(&status, &hPlus, lgthDataSet);
  if(lalDebugLevel > 2)
    REPORTSTATUS (&status);
  LALDCreateVector(&status, &hCross, lgthDataSet);
  if(lalDebugLevel > 2)
    REPORTSTATUS (&status);
  LALDCreateVector(&status, &cwSignal, lgthDataSet);
  if(lalDebugLevel > 2)
    REPORTSTATUS (&status);
  LALCreateVector (&status, &normalDeviates, lgthDataSet);
  if(lalDebugLevel > 2)
    REPORTSTATUS (&status);
  LALDCreateVector(&status, &time, lgthDataSet);
  if(lalDebugLevel > 2)
    REPORTSTATUS (&status);
  LALCreateVector(&status, &weightVector, lgthDataSet);
  if(lalDebugLevel > 2)
    REPORTSTATUS(&status);

  /* Read in Weighting File */
  if( weightFile != NULL ) {
    for(i = 0; i<=lgthDataSet; i++) {
      fgets(lineString->data, lineSize, weightFile);
      strncpy(valueString->data, lineString->data, valueSize);
      weightVector->data[i] = atof( valueString->data );
    }
  } else {
    if( lalDebugLevel > 0 )
      fprintf( stderr, "Weighting Function not in use.\n" );
    for(i = 0; i<=lgthDataSet; i++)
      weightVector->data[i] = 1.0;
  }


  /* make time array */
  for(i = 0;i<lgthDataSet;i++)
    {
      time->data[i] = i*(1.0/sampleRate);
    }
  /* Create Random Number Distribution */
  if(lalDebugLevel > 0)
    fprintf( stderr, "Generating Random Numbers.\n" );

  LALCreateRandomParams (&status, &randpar, 0);
  LALNormalDeviates (&status, normalDeviates, randpar);
  for( i=0; i<lgthDataSet ; i++ )
    normalDeviates->data[i] *= randSigma;

  if(lalDebugLevel > 0)
    fprintf( stderr, "Computing Response Functions.\n" );

  LALComputeResponseFunc( &status, &pulsar, &detector, lgthDataSet, phiStart, time, fPlus, fCross );

  if(lalDebugLevel > 0)
    fprintf( stderr, "Computing Polarization Functions.\n" );

  LALComputePolarFuncs( &pulsar, &detector,  lgthDataSet, phiStart, time, frequencyParams, hPlus, hCross );
  
  /*
   * Combine Signals with Noise
   */
  if(lalDebugLevel > 0)
    fprintf( stderr, "Creating (possibly) Weighted Signal With Noise.\n" );

  for(i=0;i<lgthDataSet;i++) 
    cwSignal->data[i] = weightVector->data[i]*(fPlus->data[i]*hPlus->data[i] + fCross->data[i]*hCross->data[i] 
					       + normalDeviates->data[i]);

  if(lalDebugLevel > 0)
    fprintf( stderr, "Writing Results to File.\n" );
  /*
   * Write out Response Functions
   */
  if(outFileFPlus != NULL)
    {
      for(i = 0;i<lgthDataSet; i++)
	fprintf(outFileFPlus, "%.15f\n", fPlus->data[i] );
    }
  if(outFileFCross != NULL)
    {
      for(i = 0;i<lgthDataSet; i++)
	fprintf(outFileFCross, "%.15f\n", fCross->data[i] );
    }

  /*
   * Write out Polarization Functions
   */
  if(outFileHPlus != NULL)
    {
      for(i = 0;i<lgthDataSet; i++)
	fprintf(outFileHPlus, "%.15f\n", hPlus->data[i] );
    }
  if(outFileHCross != NULL)
    {
      for(i = 0;i<lgthDataSet; i++)
	fprintf(outFileHCross, "%.15f\n", hCross->data[i] );
    }
  
  /*
   * Write out Noise
   */
  if(outFileNoise != NULL)
    {
      for(i = 0;i<lgthDataSet; i++)
	fprintf(outFileNoise, "%.15f\n", normalDeviates->data[i] );
    }

  /*
   * Write out Signal
   */
  if(outFileSignal != NULL)
    {
      for(i = 0;i<lgthDataSet; i++)
	fprintf(outFileSignal, "%.15f\n", cwSignal->data[i] );
    }

  LALCHARDestroyVector( &status, &lineString ); 
  LALCHARDestroyVector( &status, &valueString ); 
  LALDDestroyVector( &status, &fPlus );
  LALDDestroyVector( &status, &fCross );
  LALDDestroyVector( &status, &hPlus );
  LALDDestroyVector( &status, &hCross );
  LALDDestroyVector( &status, &cwSignal );
  LALDestroyVector ( &status, &normalDeviates );
  LALDestroyVector( &status, &weightVector );
  LALDestroyRandomParams (&status, &randpar );
  
  LALCheckMemoryLeaks();
  return 0;
}

void LALComputeResponseFunc(LALStatus *status, LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet, 
			    REAL8 phiStart, REAL8Vector *time,  REAL8Vector *fPlus, REAL8Vector *fCross)
{
  REAL8Vector *aCoeff = NULL, *bCoeff = NULL;
  /* REAL8 gammaAngle = detector->frDetector.xArmAzimuthRadians + LAL_PI_4; */ /* 05/15/03 gam */
  REAL8 gammaAngle = LAL_PI_2 + LAL_PI_4 - detector->frDetector.xArmAzimuthRadians;
  /* REAL8 detLatitude = LAL_PI_180*detector->frDetector.vertexLatitudeDegrees; */ /* 05/14/03 gam */
  REAL8 detLatitude = detector->frDetector.vertexLatitudeRadians;
  INT4 i;

  INITSTATUS (status, "LALComputeResponseFunc", GENERATECWSIGNALC);
  ATTATCHSTATUSPTR(status);

  /* a(t) and b(t) coefficients from JKS */
  LALDCreateVector( status->statusPtr, &aCoeff, lgthDataSet ); 
  LALDCreateVector( status->statusPtr, &bCoeff, lgthDataSet );

  for(i = 0;i<lgthDataSet; i++)
    {
      aCoeff->data[i] = sin(2.0*gammaAngle)*(3.0 - cos(2.0*detLatitude))*(3.0 - cos(2.0*pulsar->equatorialCoords.latitude))
	*cos(2.0*(pulsar->equatorialCoords.longitude - phiStart - omegaEarth*time->data[i]))/16.0; 
      aCoeff->data[i] -= 0.25*cos(2.0*gammaAngle)*sin(detLatitude)*(3 - cos(2.0*pulsar->equatorialCoords.latitude))
	*sin(2.0*(pulsar->equatorialCoords.longitude - phiStart - omegaEarth*time->data[i]));
      aCoeff->data[i] += 0.25*sin(2.0*gammaAngle)*sin(2.0*detLatitude)*sin(2.0*pulsar->equatorialCoords.latitude)
	*cos(pulsar->equatorialCoords.longitude - phiStart - omegaEarth*time->data[i]); 
      aCoeff->data[i] -= 0.5*cos(2.0*gammaAngle)*cos(detLatitude)*sin(2.0*pulsar->equatorialCoords.latitude)
	*sin(pulsar->equatorialCoords.longitude - phiStart - omegaEarth*time->data[i]);
      aCoeff->data[i] += 0.75*sin(2.0*gammaAngle)*cos(detLatitude)*cos(detLatitude)
	*cos(pulsar->equatorialCoords.latitude)*cos(pulsar->equatorialCoords.latitude);
      
      bCoeff->data[i] = cos(2*gammaAngle)*sin(detLatitude)*sin(pulsar->equatorialCoords.latitude)
	*cos(2.0*(pulsar->equatorialCoords.longitude-phiStart-omegaEarth*time->data[i]));
      bCoeff->data[i] += 0.25*sin(2.0*gammaAngle)*(3.0-cos(2.0*detLatitude))*sin(pulsar->equatorialCoords.latitude)
	*sin(2.0*(pulsar->equatorialCoords.longitude-phiStart-omegaEarth*time->data[i]));
      bCoeff->data[i] += cos(2*gammaAngle)*cos(detLatitude)*cos(pulsar->equatorialCoords.latitude)
	*cos(pulsar->equatorialCoords.longitude-phiStart-omegaEarth*time->data[i]);
      bCoeff->data[i] += 0.5*sin(2.0*gammaAngle)*sin(2.0*detLatitude)*cos(pulsar->equatorialCoords.latitude)
	*sin(pulsar->equatorialCoords.longitude-phiStart-omegaEarth*time->data[i]);
      
    }

  /* Compute Fplus and Fcross from JKS */
  
  /* Note: Assumes Interferometer Arms are Perpendicular */
  for(i = 0;i<lgthDataSet;i++)
    {
      fPlus->data[i] = aCoeff->data[i]*cos(2.0*pulsar->orientation) + bCoeff->data[i]*sin(2.0*pulsar->orientation);
      fCross->data[i] = bCoeff->data[i]*cos(2.0*pulsar->orientation) - aCoeff->data[i]*sin(2.0*pulsar->orientation);
    }
  
  LALDDestroyVector( status, &aCoeff );
  LALDDestroyVector( status, &bCoeff );
  
  DETATCHSTATUSPTR(status);
  LALCheckMemoryLeaks();
  
}


/*
 *  The strain of the detector is given by JKS eq. (20), and the polarization
 *  magnitudes by JKS (21) & (22)
 */

void LALComputePolarFuncs( LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet, REAL8 phiStart, REAL8Vector *time, 
			   REAL8Vector *freqParams, REAL8Vector *hPlus, REAL8Vector *hCross)
{
  REAL8 wavePhase = 0.0;
  /* REAL8 detLatitude = LAL_PI_180*detector->frDetector.vertexLatitudeDegrees; */ /* 05/14/03 gam */
  REAL8 detLatitude = detector->frDetector.vertexLatitudeRadians;
  REAL8 epsilon = 23.5*LAL_PI_180;
  INT4 i;

  for(i = 0; i<lgthDataSet; i++)
    {
      wavePhase = LAL_AU_SI*(cos(pulsar->equatorialCoords.longitude)*cos(pulsar->equatorialCoords.latitude)
			     *cos(phiSun + omegaEarthSun*time->data[i]) + (cos(epsilon)*sin(pulsar->equatorialCoords.longitude)
									 *cos(pulsar->equatorialCoords.latitude) 
									 + sin(epsilon)*sin(pulsar->equatorialCoords.latitude))
			     *sin(phiSun + omegaEarthSun*time->data[i]));
      
      wavePhase += LAL_REARTH_SI*(sin(detLatitude)*sin(pulsar->equatorialCoords.latitude)
				  + cos(detLatitude)*cos(pulsar->equatorialCoords.latitude)
				  *cos(pulsar->equatorialCoords.longitude - phiStart - omegaEarth*time->data[i]));

      wavePhase *= (LAL_TWOPI/LAL_C_SI)*(freqParams->data[0] + freqParams->data[1]*time->data[i] 
					 + freqParams->data[2]*time->data[i]*time->data[i]/2.0 
					 + freqParams->data[3]*time->data[i]*time->data[i]*time->data[i]/6.0
					 + freqParams->data[4]*pow(time->data[i], 4)/24.0
					 + freqParams->data[5]*pow(time->data[i], 5)/120.0);

      wavePhase += LAL_TWOPI*(freqParams->data[0]*time->data[i] + (freqParams->data[1]*time->data[i]*time->data[i]/2.0) 
			      + (freqParams->data[2]*time->data[i]*time->data[i]*time->data[i]/6.0)
			      + (freqParams->data[3]*pow(time->data[i], 4)/24.0)
			      + (freqParams->data[4]*pow(time->data[i], 5)/120.0)
			      + (freqParams->data[5]*pow(time->data[i], 6)/720.0) );		     

      hPlus->data[i] = plusSignalAmplitude*cos(wavePhase); 
      hCross->data[i] = crossSignalAmplitude*sin(wavePhase); 
      
      wavePhase = 0.0;
    }
  
  LALCheckMemoryLeaks();
  
}

static void
PrintLALDetector(LALDetector *detector)
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

/* This function is based on the code from Brian Cameron's LALResponseFunc. c */ /* 05/20/03 gam */
int GenerateResponseFuncUsingLAL(int argc, char *argv[])
{
  static LALStatus  status;
  LALSource         pulsar;
  LALDetector       detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
  LIGOTimeGPS       gps;
  LALDetAndSource   det_and_pulsar;
  LALDetAMResponse  am_response;
  LALPlaceAndGPS    place_and_gps; REAL8 lmsthours;
  LALDate date;
  LALDetAMResponseSeries    am_response_series = {NULL,NULL,NULL};
  REAL4TimeSeries           plus_series, cross_series, scalar_series;
  LALTimeIntervalAndNSample time_info;
  LALLeapSecAccuracy leapAccuracy = LALLEAPSEC_STRICT;
  CHARVector          *timestamp = NULL;
  UINT4 i;

  LALMSTUnitsAndAcc uandacc = { MST_RAD, LALLEAPSEC_STRICT}; /* 05/20/03 gam */
  LALGPSandAcc gps_and_acc; /* 05/20/03 gam */

  FILE *outFilePlus, *outFileCross;
  outFilePlus = fopen("fPlusLAL.out", "w");
  outFileCross = fopen("fCrossLAL.out", "w");

  if (argc == 2)
    lalDebugLevel = atoi(argv[1]);

  /*
   * Set Siurce as Crab Pulsar
   */
  strcpy(pulsar.name, "CRAB PULSAR");
  pulsar.equatorialCoords.longitude = 83.63*LAL_PI_180;  /* RA */
  pulsar.equatorialCoords.latitude  = 22.014*LAL_PI_180;  /* Dec */
  pulsar.equatorialCoords.system    = COORDINATESYSTEM_EQUATORIAL;
  pulsar.orientation                = 22.5*LAL_PI_180;  /* orientation */

  /*
   * Print Source Info
   */
  printf("Source Info\n");
  printf("RA: %e\n", pulsar.equatorialCoords.longitude);
  printf("Dec: %e\n", pulsar.equatorialCoords.latitude);
  printf("Psi: %e\n", pulsar.orientation);

  if (lalDebugLevel > 0)
    printf("Det #1 location: (%7.4e, %7.4e, %7.4e)\n",
	   detector.location[0], detector.location[1],
           detector.location[2]);

  /*
   * Set Date
   */
  LALCHARCreateVector(&status, &timestamp, (UINT4)64);

  date.unixDate.tm_year = 100;
  date.unixDate.tm_mon  =  0;
  date.unixDate.tm_mday =  1;
  date.unixDate.tm_hour =  0;
  date.unixDate.tm_min  =  0;
  date.unixDate.tm_sec  =  0;
  date.residualNanoSeconds = 0;

  LALUTCtoGPS (&status, &gps, &date, &leapAccuracy);

  LALDateString(&status, timestamp, &date);

  gps.gpsSeconds = 709398013;
  gps.gpsNanoSeconds = 0;

  fprintf(stderr, "For: %s\n", timestamp->data);
  fprintf(stderr, "GPS = {%10d, %9d}\n", gps.gpsSeconds,
          gps.gpsNanoSeconds);

  /*
   * Stuff Detector and Source into structure
   */
  place_and_gps.p_detector = &detector;
  place_and_gps.p_gps   = &gps;
  /* LALGPStoLMST1(&status, &lmsthours, &place_and_gps, MST_RAD); */ /* 05/20/03 gam */
  LALGPStoLMST1(&status, &lmsthours, &place_and_gps, &uandacc);

  fprintf(stderr, "LMST = %7e \n", lmsthours);

  /*
   * Stuff Detector and Source into structure
   */
  det_and_pulsar.pDetector = &detector;
  det_and_pulsar.pSource   = &pulsar;

  /*
   * Compute instantaneous AM response
   */
  gps_and_acc.gps = gps; /* 05/20/03 gam */
  gps_and_acc.accuracy = LALLEAPSEC_STRICT; /* 05/20/03 gam */

  /* LALComputeDetAMResponse(&status, &am_response, &det_and_pulsar, &gps); */ /* 05/20/03 gam */
  LALComputeDetAMResponse(&status, &am_response, &det_and_pulsar, &gps_and_acc);

  if (status.statusCode && lalDebugLevel > 0)
    {
      /* fprintf(stderr,
              "LALTestDetResponse0: error in LALComputeDetAMResponse, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C); */ /* 05/20/03 gam */
      fprintf(stderr, "Error in LALComputeDetAMResponse\n");
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  /*
   * Compute a time series AM response
   */
  if (lalDebugLevel > 0)
    {
      printf("Starting vector test\n");
    }

  plus_series.data = NULL;
  cross_series.data = NULL;
  scalar_series.data = NULL;

  am_response_series.pPlus   = &(plus_series);
  am_response_series.pCross  = &(cross_series);
  am_response_series.pScalar = &(scalar_series);

  LALSCreateVector(&status, &(am_response_series.pPlus->data), 1);
  LALSCreateVector(&status, &(am_response_series.pCross->data), 1);
  LALSCreateVector(&status, &(am_response_series.pScalar->data), 1);

  if (lalDebugLevel > 0)
    {
      printf("am_response_series.pPlus->data->length = %d\n",
             am_response_series.pPlus->data->length);
      printf("am_response_series.pCros->data->length = %d\n",
             am_response_series.pCross->data->length);
      printf("am_response_series.pScalar->data->length = %d\n",
             am_response_series.pScalar->data->length);
    }

  time_info.epoch.gpsSeconds     = gps.gpsSeconds;
  time_info.epoch.gpsNanoSeconds = gps.gpsNanoSeconds;
  time_info.deltaT               = 1;
  time_info.nSample              = 86401;

  LALComputeDetAMResponseSeries(&status,
                                &am_response_series,
                                &det_and_pulsar,
                                &time_info);

  if (status.statusCode && lalDebugLevel > 0)
    {
      /* fprintf(stderr,
              "LALTestDetResponse0: error in LALComputeDetAMResponseSeries, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSE0C); */ /* 05/20/03 gam */
      fprintf(stderr, "Error in LALComputeDetAMResponse\n");
      return status.statusCode;
    }

  if (lalDebugLevel > 0)
    {
      printf("Done computing AM response vectors\n");

      printf("am_response_series.pPlus->data->length = %d\n",
             am_response_series.pPlus->data->length);
      printf("am_response_series.pCross->data->length = %d\n",
             am_response_series.pCross->data->length);
      printf("am_response_series.pScalar->data->length = %d\n",
             am_response_series.pScalar->data->length);
    }


  for (i = 0; i < time_info.nSample; ++i)
    {
      fprintf(outFilePlus, "%f\n", am_response_series.pPlus->data->data[i]);
      fprintf(outFileCross, "%f\n", am_response_series.pCross->data->data[i]);
    }
  LALCHARDestroyVector(&status, &timestamp);
  LALSDestroyVector(&status, &(am_response_series.pPlus->data));
  LALSDestroyVector(&status, &(am_response_series.pCross->data));
  LALSDestroyVector(&status, &(am_response_series.pScalar->data));

  LALCheckMemoryLeaks();

  return 0;
}

