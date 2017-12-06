/*
*  Copyright (C) 2007 David Chin, Gregory Mendell
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

/*****************************************************************
File: LALIndependentTestDetReponce.c

Authors: Greg Mendell, Malik Rakhmanov Started: May 20 2003.

Code: Independently tests of LALComputeDetAMResponse using
model for F_+ and F_x given in Jaranowski, Krolak, and Schutz gr-qc/9804014.

Original Author: Patrick Brian Cameron.  Started: SURF Project Summer 2002.
Original code: generated test of continuous graviational wave signals using
the model given in Jaranowski, Krolak, and Schutz gr-qc/9804014.

*******************************************************************/

/* Changes:
05/14/03 gam: Change vertexLongitudeDegrees and vertexLatitudeDegrees to vertexLongitudeRadians vertexLatitudeRadians when using FrDetector structure.
05/15/03 gam: xArmAzimuthRadians and yArmAzimuthRadians are now measured clockwise from North; lal_gamma is the angle to the bisector of the arms, measured
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
10/13/03 gam: Fix bug - when some output files requested but not others, no output occurs.
10/13/03 gam: Use independent code from Jolien Creighton to convert GPS to Sidereal time.
10/13/03 gam: Include GEO detector as an option.
10/13/03 gam: Use constants independent of LAL.
10/14/04 gam: Update definition of lal_gamma when angle between arms, zeta, != pi/2.
10/14/04 gam: Change input RA, DEC and orientation angle (polarization angle) in config file to be in radians.
10/14/04 gam: Use independent detector geometry values when doing independent calculation.
10/15/04 gam: Fix bug M_PI not defined when configuring lal with --enable-gcc-flags.
              WARNING: LHO AND LLO VALUES WERE TAKEN FROM OTHER LIGO SOURCES; GEO VALUES ARE NOT INDEPENDENT BUT TAKEN FROM LAL
*/

#include <config.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
/* #include <lalapps.h> */ /* 05/20/03 gam */
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
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
     "      the command line can be as long as you like.\n"

#define usage( program ) fprintf( stdout, usgfmt, program )

/* 10/13/04 gam; independent values for pi/2, pi, and 2pi */
#define LALIND_PI_2  1.57079633
#define LALIND_PI    3.14159265
#define LALIND_TWOPI 6.28318531
/* 10/13/04 gam; allowed timing difference between LMST from LAL and independent code */
#define LALIND_TIMETOLERANCE 32

/* void GenerateResponseFuncNotUsingLAL(LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet, REAL8 phiStart, REAL8Vector *timevec, REAL8Vector *fPlus, REAL8Vector *fCross); */ /* 10/14/04 gam */
void GenerateResponseFuncNotUsingLAL(LALSource *pulsar, REAL8 inputXArmAzimuthRadians, REAL8 inputVertexLatitudeRadians, INT4 lgthDataSet, REAL8 phiStart, REAL8Vector *timevec,  REAL8Vector *fPlus, REAL8Vector *fCross);

void PrintLALDetector(const LALDetector *detector);

void GenerateResponseFuncUsingLAL(LALStatus *status, LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet, REAL8 sampleRate, LIGOTimeGPS *gps, LALDetAMResponseSeries  *am_response_series_ptr);

/* 10/13/03 gam: Use independent code from Jolien Creighton to convert GPS to Sidereal time. */
REAL8 greenwich_mean_sidereal_time( INT4 gpssec, INT4 gpsnan, INT4 taiutc );

extern char *optarg;
/* REAL8 omegaEarth = LAL_TWOPI/LAL_DAYSID_SI; */
/* REAL8 omegaEarthSun = LAL_TWOPI/LAL_YRSID_SI; */
REAL8 omegaEarth = LALIND_TWOPI/(23.0*3600.0 + 56.0*60.0 + 4.091);  /* Number of ordinary seconds in one sidereal day */
REAL8 omegaEarthSun = LALIND_TWOPI/(365.0*86400.0 + 6.0*3600.0 + 9.0*60.0 + 10.0); /* Number of ordinary seconds in one sidereal year */
REAL8 zeta = ((REAL8)LALIND_PI_2);               /* 10/13/04 gam; angle between the detector arms */
REAL8 sinZeta = 1.0;                             /* 10/13/04 gam */
REAL8 zetaGEO = 94.33*((REAL8)LALIND_PI)/180.0;  /* 10/13/04 gam; value from JKS paper */

int main( int argc, char *argv[] )
{
  /* Parsing Command Line, and General Use */
  const CHAR *program = argv[0];
  static LALStatus status;
  INT4 opt;
  UINT4 i;
  char *s;
  int rc;

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

  /* 10/14/04 gam; for independent values for detector geometry: */
  REAL8 indXArmAzimuthRadians = 0.0;
  REAL8 indVertexLongitudeRadians = 0.0;
  REAL8 indVertexLatitudeRadians = 0.0;

  /* Vectors For Results */
  REAL8Vector *fPlus=NULL, *fCross=NULL;

  /* Other Stuff */
  REAL8Vector *timevec = NULL;
  REAL8 sampleRate = 0, duration = 0;
  UINT4 lgthDataSet = 0;
  REAL8 phiStart = 0.0;
  REAL8 phiStartLAL = 0.0;  /* 10/13/04 gam */

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
    return 1;
  }

  LALCHARCreateVector(&status, &lineString, lineSize);
  LALCHARCreateVector(&status, &valueString, valueSize);

  /* Throw Away First Line Comment */
  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }

  /* Test Name */
  strcpy(pulsar.name, "TEST");
  pulsar.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;

  /* R.A. of the Source */
  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }
  strncpy(valueString->data, lineString->data, valueSize);
  /* pulsar.equatorialCoords.longitude = LAL_PI_180*(atof( valueString->data )); */
  pulsar.equatorialCoords.longitude = atof( valueString->data ); /* 10/14/04 gam; input needs to be in radians already */

  /* Declination of the Source */
  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }
  strncpy(valueString->data, lineString->data, valueSize);
  /* pulsar.equatorialCoords.latitude = LAL_PI_180*(atof( valueString->data )); */
  pulsar.equatorialCoords.latitude = atof( valueString->data ); /* 10/14/04 gam; input needs to be in radians already */

  /* Polarization Angle */
  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }
  strncpy(valueString->data, lineString->data, valueSize);
  /* pulsar.orientation = LAL_PI_180*atof( valueString->data ); */
  pulsar.orientation = atof( valueString->data ); /* 10/14/04 gam; input needs to be in radians already */

  /* Start Time */
  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }
  strncpy(valueString->data, lineString->data, valueSize);
  /* gpsTime.gpsSeconds = atoi( valueString->data ); */
  gpsTime.gpsSeconds = atol( valueString->data ); /* 09/30/03 gam */
  gpsTime.gpsNanoSeconds = 0;

  /* Sample Rate */
  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }
  strncpy(valueString->data, lineString->data, valueSize);
  sampleRate = atof( valueString->data );

  /* Duration */
  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }
  strncpy(valueString->data, lineString->data, valueSize);
  duration = atof( valueString->data );

  /* lgthDataSet = (UINT4)(duration * sampleRate) + 1; */ /* 09/30/03 gam */
  lgthDataSet = (UINT4)(duration * sampleRate);

  /* Detector Site */
  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }
  strncpy(valueString->data, lineString->data, valueSize);
  if(valueString->data[0] == 'H') {
    detector = lalCachedDetectors[LALDetectorIndexLHODIFF];
    /* 10/14/04 gam; fill in independent values for detector geometry: */
    indXArmAzimuthRadians = 5.6548781395;
    indVertexLongitudeRadians = -2.08405709267;
    indVertexLatitudeRadians = 0.810795009136;
  }
  if(valueString->data[0] == 'L') {
    detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
    /* 10/14/04 gam; fill in independent values for detector geometry: */
    indXArmAzimuthRadians = 4.40317821503;
    indVertexLongitudeRadians = -1.5843089819;
    indVertexLatitudeRadians = 0.533423006535;
  }
  if(valueString->data[0] == 'G') {
    detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
    /* 10/14/04 gam; fill in values for detector geometry: WARNING, THESE ARE NOT INDEPENDENT BUT TAKEN FROM LAL */
    indXArmAzimuthRadians = 1.1936010122;
    indVertexLongitudeRadians = 0.17116780435;
    indVertexLatitudeRadians = 0.91184982752;
    zeta = zetaGEO;         /* 10/13/04 gam; use value from above from JKS paper */
    sinZeta = sin(zetaGEO); /* 10/13/04 gam; use value from above from JKS paper */
  }
  if(valueString->data[0] == 'D') {
      fprintf( stdout,  "Enter Properties of test detector\n" );
      /* fprintf( stderr, "Detector Longitude (deg): " );
      rc = scanf( "%lf", &testFrDetector.vertexLongitudeDegrees ); */ /* 05/14/03 gam */
      fprintf( stdout, "Detector Longitude (radians): " );
      rc = scanf( "%lf", &testFrDetector.vertexLongitudeRadians );
      if (rc != 1)
      {
        fprintf(stderr, "Error: Unable to read input\n");
        exit(1);
      }

      /* fprintf( stderr, "Detector Latitude (deg): " );
      rc = scanf( "%lf" , &testFrDetector.vertexLatitudeDegrees ); */ /* 05/14/03 gam */
      /* if (rc != 1)
      {
        fprintf(stderr, "Error: Unable to read input\n");
        exit(1);
      } */
      fprintf( stdout, "Detector Latitude (radians): " );
      rc = scanf( "%lf" , &testFrDetector.vertexLatitudeRadians );
      if (rc != 1)
      {
        fprintf(stderr, "Error: Unable to read input\n");
        exit(1);
      }

      /* fprintf( stderr, "Detector xArmAzimuth (deg): " ); */ /* 05/14/03 gam */
      fprintf( stdout, "Detector xArmAzimuth (radians): " );
      rc = scanf( "%f", &testFrDetector.xArmAzimuthRadians );
      if (rc != 1)
      {
        fprintf(stderr, "Error: Unable to read input\n");
        exit(1);
      }

      /* testFrDetector.xArmAzimuthRadians *= LAL_PI_180; */ /* 05/14/03 gam */
      /* testFrDetector.yArmAzimuthRadians = testFrDetector.xArmAzimuthRadians + LAL_PI_2; */ /* 05/15/03 gam */
      /* testFrDetector.yArmAzimuthRadians = testFrDetector.xArmAzimuthRadians - LAL_PI_2; */ /* 10/13/04 gam */
      testFrDetector.yArmAzimuthRadians = testFrDetector.xArmAzimuthRadians - zeta;
      testFrDetector.vertexElevation = 0.0;
      testFrDetector.xArmAltitudeRadians = 0.0;
      testFrDetector.yArmAltitudeRadians = 0.0;
      LALCreateDetector( &status, &detector, &testFrDetector, LALDETECTORTYPE_IFODIFF);
      /* 10/14/04 gam; fill in independent values for detector geometry: */
      indXArmAzimuthRadians = testFrDetector.xArmAzimuthRadians;
      indVertexLongitudeRadians = testFrDetector.vertexLongitudeRadians;
      indVertexLatitudeRadians = testFrDetector.vertexLatitudeRadians;
  }

  /* tolerance */
  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }
  strncpy(valueString->data, lineString->data, valueSize);
  tolerance = atof( valueString->data );

  /* Open Various Output Files */

  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }
  strncpy(valueString->data, lineString->data, valueSize);
  outputIndependentFPlusFCross = atoi( valueString->data );

  if (outputIndependentFPlusFCross) {
    /*  F_Plus Filename */
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
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
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
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
 } else {
    /* 10/13/03 gam; still need to get filenames */
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
 }

  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }
  strncpy(valueString->data, lineString->data, valueSize);
  outputLALFPlusFCross = atoi( valueString->data );

  if (outputLALFPlusFCross) {
    /*  F_Plus LAL Filename */
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
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
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
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
 } else {
    /* 10/13/03 gam; still need to get filenames */
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
 }

  s = fgets(lineString->data, lineSize, configFile);
  if (s == NULL)
  {
    fprintf(stderr, "Error: Unable to read ConfigFile\n");
    exit(1);
  }
  strncpy(valueString->data, lineString->data, valueSize);
  outputFPlusFCrossDiffs = atoi( valueString->data );

  if (outputFPlusFCrossDiffs) {

    /*  F_Plus LAL Filename */
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
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
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
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
 } else {
    /* 10/13/03 gam; still need to get filenames */
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
    s = fgets(lineString->data, lineSize, configFile);
    if (s == NULL)
    {
      fprintf(stderr, "Error: Unable to read ConfigFile\n");
      exit(1);
    }
 }

 fclose(configFile);

  /* Echo Parameters */
  if( lalDebugLevel > 0)
    {
      /* 10/14/04 gam; reformat: */
      fprintf( stdout,  "\n" );
      fprintf( stdout, "Source Info:\n" );
      fprintf( stdout, "RA (radians): %e\n", pulsar.equatorialCoords.longitude );
      fprintf( stdout, "Dec (radians): %e\n", pulsar.equatorialCoords.latitude );
      fprintf( stdout, "Polarization (radians): %f\n", pulsar.orientation );
      fprintf( stdout, "For GPS Seconds, NanoSeconds: %10d, %10d\n", gpsTime.gpsSeconds, gpsTime.gpsNanoSeconds );
      fprintf( stdout, "Sampling Rate (Hz): %f\n", sampleRate );
      fprintf( stdout, "Length of Data Set (seconds): %f\n", duration );
      fprintf( stdout, "Ang. Vel. of Earth (radians/s) = %e\n", omegaEarth );
      fprintf( stdout, "Ang. Vel. of Earth about the Sun (radians/s) = %e\n", omegaEarthSun );
      fprintf( stdout, "Independent detector vertex longitude: %.11g\n", indVertexLongitudeRadians ); /* 10/14/04 gam */
      fprintf( stdout, "Independent detector vertex latitude: %.11g\n", indVertexLatitudeRadians  );  /* 10/14/04 gam */
      fprintf( stdout, "Independent detector xArm Azimuth: %.11g\n", indXArmAzimuthRadians );         /* 10/14/04 gam */
      PrintLALDetector( &detector );
    }

  /*
   * Start Computing Signal
   */
  place_and_gps.p_detector = &detector;
  place_and_gps.p_gps = &gpsTime;

  /* Compute Local Sidereal Time at the start of the data set*/

  phiStartLAL = XLALGreenwichMeanSiderealTime(place_and_gps.p_gps) + atan2(place_and_gps.p_detector->location[1], place_and_gps.p_detector->location[0]);
  phiStartLAL = fmod(phiStartLAL,LAL_TWOPI);  /* put into interval 0 to 2pi */
  if(lalDebugLevel > 0) {
     fprintf(stdout, "Local Mean Sidereal Time from LAL (radians) = %f\n", phiStartLAL);
     fflush(stdout);
  }

  /****************************************************************************/
  /* !!! 10/13/03 gam; get Local Mean Sidereal Time from independent code !!! */
  /* Use routine from Jolien Creighton. the last argument is integer TAI - UTC which is the number
     of leap seconds (the difference between atomic time and UTC); assume 32 is accurate enough. */
  /* Also note sidereal time is alway the angle the vernal equinoix is west of the local
     meridian; i.e., between 0 and 2pi in radians. However, the detector longitude is measured
     + when east of Greenwich and - when west of Greenwich. */ /* 10/14/04 gam; update with ind detector geometry. */
  phiStart = greenwich_mean_sidereal_time( (INT4)gpsTime.gpsSeconds, (INT4)gpsTime.gpsNanoSeconds, 32 );
  if (indVertexLongitudeRadians > 0.0) {
     phiStart = phiStart + indVertexLongitudeRadians; /* convert GMST to LMST for eastern longitude */
  } else {
     phiStart = phiStart + ((REAL8)LALIND_TWOPI) + indVertexLongitudeRadians; /* convert GMST to LMST for western logitude */
  }
  phiStart = fmod(phiStart,LALIND_TWOPI);  /* put into interval 0 to 2pi */
  if(lalDebugLevel > 0) {
    fprintf( stdout, "Local Mean Sidereal Time from independent code (radians) = %f\n", phiStart );
    fflush(stdout);
  }

  if ( fabs(phiStart - phiStartLAL)*3600.0*24.0/((REAL8)LALIND_TWOPI) > ((REAL8)LALIND_TIMETOLERANCE)) {
      fprintf(stderr, "\nERROR: Local Mean Sidereal Time from LAL and independent code differ by more than %.6g seconds!\n\n",((REAL8)LALIND_TIMETOLERANCE));
      fflush(stderr);
      exit(1);
  }

  /* Intialize Arrays */
  if(lalDebugLevel > 1) fprintf( stdout, "Creating Vectors.\n" );

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

  if(lalDebugLevel > 1) {
    fprintf(stdout, "Computing Independent Response Functions.\n" );
    fflush(stdout);
  }

  /* GenerateResponseFuncNotUsingLAL(&pulsar, &detector, lgthDataSet, phiStart, timevec, fPlus, fCross ); */ /* 10/14/04 gam */
  GenerateResponseFuncNotUsingLAL(&pulsar, indXArmAzimuthRadians, indVertexLatitudeRadians, lgthDataSet, phiStart, timevec, fPlus, fCross );

  if(lalDebugLevel > 1) {
    fprintf(stdout, "Computing Response Functions Using LAL.\n" );
    fflush(stdout);
  }

  GenerateResponseFuncUsingLAL(&status, &pulsar, &detector, lgthDataSet, sampleRate, &gpsTime, &am_response_series);

  /* Write out files */
  if (outputIndependentFPlusFCross) {
    if(lalDebugLevel > 1) fprintf( stdout, "Writing Independent F_+ and F_x to File.\n" );
    for (i = 0;i<lgthDataSet; i++) {
        fprintf(outFileFPlus, "%.15f\n", fPlus->data[i]);
        fprintf(outFileFCross, "%.15f\n", fCross->data[i]);
    }
    fclose(outFileFPlus);
    fclose(outFileFCross);
  }
  if (outputLALFPlusFCross) {
    if(lalDebugLevel > 1) fprintf( stdout, "Writing LAL F_+ and F_x to File.\n" );
    for (i = 0; i < lgthDataSet; ++i) {
        fprintf(outFileLALPlus, "%.15f\n", am_response_series.pPlus->data->data[i]);
        fprintf(outFileLALCross, "%.15f\n", am_response_series.pCross->data->data[i]);
    }
    fclose(outFileLALPlus);
    fclose(outFileLALCross);

  }

  /* Compute differences; write out differences if requested */ /* 09/30/03 gam */
  if(lalDebugLevel > 1) fprintf( stdout, "Computing differences LAL F_+ - Independent F_+ and LAL F_x - Independent F_x.\n" );
  if (outputFPlusFCrossDiffs && lalDebugLevel > 1) fprintf( stdout, "Writing differences LAL F_+ - Independent F_+ and LAL F_x - Independent F_x to File.\n" );
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
      fprintf(stderr, "\nERROR: LAL F_+ - Independent F_+ Difference Test Failed!\n\n");
      fflush(stderr);
    }
    if (maxFCrossDiff > tolerance) {
      fprintf(stderr, "\nERROR: LAL F_x - Independent F_x Difference Test Failed!\n\n");
      fflush(stderr);
    }
    exit(1);
  }

  if(lalDebugLevel > 0) {
    fprintf(stdout, "\nALL LALIndependentTestDetResponse Tests Passed!\n");
    fflush(stdout);
  }

  return 0;
}

/* void GenerateResponseFuncNotUsingLAL(LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet,
			    REAL8 phiStart, REAL8Vector *timevec,  REAL8Vector *fPlus, REAL8Vector *fCross) */
/* 10/14/04 gam; make inputXArmAzimuthRadians and inputVertexLatitudeRadians inputs */
void GenerateResponseFuncNotUsingLAL(LALSource *pulsar, REAL8 inputXArmAzimuthRadians, REAL8 inputVertexLatitudeRadians,
     INT4 lgthDataSet, REAL8 phiStart, REAL8Vector *timevec,  REAL8Vector *fPlus, REAL8Vector *fCross)
{
  REAL8 aCoeff, bCoeff;
  /* REAL8 gammaAngle = detector->frDetector.xArmAzimuthRadians + LAL_PI_4; */ /* 05/15/03 gam */
  /* REAL8 gammaAngle = LAL_PI_2 + LAL_PI_4 - detector->frDetector.xArmAzimuthRadians; */ /* 10/14/04 gam */
  /* REAL8 gammaAngle = ((REAL8)LALIND_PI_2) + zeta/2.0 - detector->frDetector.xArmAzimuthRadians; */ /* 10/14/04 gam */
  REAL8 gammaAngle = ((REAL8)LALIND_PI_2) + zeta/2.0 - inputXArmAzimuthRadians;

  /* REAL8 detLatitude = LAL_PI_180*detector->frDetector.vertexLatitudeDegrees; */ /* 05/14/03 gam */
  /* REAL8 detLatitude = detector->frDetector.vertexLatitudeRadians; */ /* 10/14/04 gam */
  REAL8 detLatitude = inputVertexLatitudeRadians;
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
      fPlus->data[i] = sinZeta*( aCoeff*cos(2.0*pulsar->orientation) + bCoeff*sin(2.0*pulsar->orientation) );  /* 10/13/04 gam; add sinZeta */
      fCross->data[i] = sinZeta*( bCoeff*cos(2.0*pulsar->orientation) - aCoeff*sin(2.0*pulsar->orientation) ); /* 10/13/04 gam; add sinZeta */
  }
  return;
}

void PrintLALDetector(const LALDetector *detector)
{
  /* REAL8 lal_gamma = detector->frDetector.xArmAzimuthRadians + LAL_PI_4; */ /* 05/15/03 gam */
  /* REAL8 lal_gamma = LAL_PI_2 + LAL_PI_4 - detector->frDetector.xArmAzimuthRadians; */ /* 10/14/04 gam */
  REAL8 lal_gamma = ((REAL8)LAL_PI_2) + zeta/2.0 - detector->frDetector.xArmAzimuthRadians;
  REAL8 computedZeta = fabs(detector->frDetector.xArmAzimuthRadians - detector->frDetector.yArmAzimuthRadians); /* 10/14/04 gam */

  /* 10/14/04 gam; add in rest of structure and reformat: */

  if ( computedZeta > ((REAL8)LAL_PI) ) {
       computedZeta = ((REAL8)LAL_TWOPI) - computedZeta;
  }

  fprintf( stdout,"Detector name: %s,\n", detector->frDetector.name);
  fprintf( stdout,"Detector struct location array  = { %.6e, %.6e, %.6e },\n", detector->location[0],
           detector->location[1], detector->location[2] );

  /* fprintf( stdout, "Location =   %.6g, %.6g, %.6g,\n",
	   detector->frDetector.vertexLongitudeDegrees,
	   detector->frDetector.vertexLatitudeDegrees,
	   detector->frDetector.vertexElevation); */ /* 05/14/03 gam */
  fprintf( stdout, "frDetector struct vertex long. (radians), lat. (radians), elev. =   %.11g, %.11g, %.11g\n",
	   detector->frDetector.vertexLongitudeRadians,
	   detector->frDetector.vertexLatitudeRadians,
	   detector->frDetector.vertexElevation);

  fprintf( stdout, "frDetector struct xarm azimuth, yarm azimuth, xarm altitude (clockwise from north), yarm altitude (all in radians) =  %.11g, %.11g, %.11g, %.11g\n",
           detector->frDetector.xArmAzimuthRadians,
           detector->frDetector.yArmAzimuthRadians,
           detector->frDetector.xArmAltitudeRadians,
           detector->frDetector.yArmAltitudeRadians);

  fprintf( stdout, "Detector angle between arms and angle to arm bisector counter-clockwise from east = %.11g, %.11g\n", computedZeta, lal_gamma);
  fprintf( stdout, "\n" );
  return;
}

void GenerateResponseFuncUsingLAL(LALStatus *status, LALSource *pulsar, LALDetector *detector, INT4 lgthDataSet, REAL8 sampleRate, LIGOTimeGPS *gps, LALDetAMResponseSeries  *am_response_series_ptr)
{
  LALDetAndSource   det_and_pulsar;
  LALDetAMResponse  am_response;
  LALPlaceAndGPS    place_and_gps; REAL8 lmsthours;

  LALTimeIntervalAndNSample time_info;

  INITSTATUS(status);
  ATTATCHSTATUSPTR(status);

  if (lalDebugLevel > 1)  {
     /* fprintf(stdout,"status->statusCode = %i \n", status->statusCode);
     fprintf(stdout,"status->statusPtr->statusCode = %i \n", status->statusPtr->statusCode); */
     printf("Source Info\n");
     printf("RA: %e\n", pulsar->equatorialCoords.longitude);
     printf("Dec: %e\n", pulsar->equatorialCoords.latitude);
     printf("Psi: %e\n", pulsar->orientation);
     printf("Detector location: (%7.4e, %7.4e, %7.4e)\n",
	   detector->location[0], detector->location[1],
           detector->location[2]);
     fprintf(stdout, "GPS = {%10d, %9d}\n", gps->gpsSeconds, gps->gpsNanoSeconds);
  }

  /* Add detector and source into structures */
  place_and_gps.p_detector = detector;
  place_and_gps.p_gps = gps;
  det_and_pulsar.pDetector = detector;
  det_and_pulsar.pSource   = pulsar;

  lmsthours = XLALGreenwichMeanSiderealTime(place_and_gps.p_gps) + atan2(place_and_gps.p_detector->location[1], place_and_gps.p_detector->location[0]);

  if (lalDebugLevel > 1)  {
    fprintf(stdout, "In GenerateResponseFuncUsingLAL LMST = %7e \n", lmsthours); /* 10/13/04 gam */
  }

  /* LALComputeDetAMResponse(&status, &am_response, &det_and_pulsar, &gps); */ /* 05/20/03 gam */
  LALComputeDetAMResponse(status->statusPtr, &am_response, &det_and_pulsar, gps);
  CHECKSTATUSPTR (status);

  time_info.epoch.gpsSeconds     = gps->gpsSeconds;
  time_info.epoch.gpsNanoSeconds = gps->gpsNanoSeconds;
  time_info.deltaT               = 1.0/sampleRate;
  time_info.nSample              = lgthDataSet;

   if (lalDebugLevel > 1) {
      printf("Start computing AM response vectors\n");
      /* fprintf(stdout, "am_response_series_ptr->pPlus->data->length = %d\n", am_response_series_ptr->pPlus->data->length);
      fprintf(stdout, "am_response_series_ptr->pCross->data->length = %d\n", am_response_series_ptr->pCross->data->length);
      fprintf(stdout, "am_response_series_ptr->pScalar->data->length = %d\n", am_response_series_ptr->pScalar->data->length); */
      fprintf(stdout, "timeinfo = %d %d %g %d\n",time_info.epoch.gpsSeconds,time_info.epoch.gpsNanoSeconds,time_info.deltaT,time_info.nSample);
      fflush(stdout);
   }

  LALComputeDetAMResponseSeries(status->statusPtr, am_response_series_ptr, &det_and_pulsar, &time_info);
  CHECKSTATUSPTR (status);

  if (lalDebugLevel > 1) {
      printf("Done computing AM response vectors\n");
      /* printf("am_response_series_ptr->pPlus->data->length = %d\n", am_response_series_ptr->pPlus->data->length);
      printf("am_response_series_ptr->pCross->data->length = %d\n", am_response_series_ptr->pCross->data->length);
      printf("am_response_series_ptr->pScalar->data->length = %d\n", am_response_series_ptr->pScalar->data->length); */
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


/* 10/13/03 gam: Use independent code from Jolien Creighton to convert GPS to Sidereal time. */
/* Notes from Jolien:
Here is my GPS -> Sidereal conversion routine, which is entirely
independent of the code in LAL.  (I didn't write the code in LAL,
and I didn't even look at it while writing this routine.)  There
are some references in comments in the code.

For some reason, I represent sidereal time in the range [0,2pi).
You can change a couple of lines at the end if you want it in a
different range.

The code was not written to be greatly efficient.  Instead it is
a rather plodding and straightforward implementation for testing
purposes.

I have checked the routine with the USNO observatory sidereal time
calculator and I find that it agrees to around ~1ms or less, I think...
(you should double check).

You need to give the routine the integer TAI - UTC which is the number
of leap seconds (the difference between atomic time and UTC).  The
value is currently 32 I think (and has been for the past year or two).
This is, of course, non-algorithmic so you just have to look it up
(and wait for new leap seconds to be announced).

If you want accuracy better than about 1s of sidereal time, you probably
need to do a very careful monitoring of the Earth's motion, so this
routine will not be sufficient.  It does give you the "correct" sidereal
time, but the exact position of quasars will drift around by up to a
fraction of a second before a leap second is added.  So unless you are
monitoring this motion, I'm not sure how meaningful it is to require
too much more accuracy in sidereal time.

Cheers
Jolien
*/
REAL8 greenwich_mean_sidereal_time( INT4 gpssec, INT4 gpsnan, INT4 taiutc )
{
   /* cf. S. Aoki et al., A&A 105, 359 (1982) eqs. 13 & 19 */
   /* also cf. http://aa.usno.navy.mil */
   /* Note: 00h UT 01 Jan 2000 has JD=2451544.5 and GPS=630720013 */
   const REAL8 JD_12h_01_Jan_2000     = 2451545.0;
   const REAL8 JD_00h_01_Jan_2000     = 2451544.5;
   const REAL8 GPS_00h_01_Jan_2000    = 630720013;
   const REAL8 TAIUTC_00h_01_Jan_2000 = 32; /* leap seconds: TAI - UTC */


   REAL8 t;
   REAL8 dpU;
   REAL8 TpU;
   REAL8 gmst;

   /* compute number of seconds since 00h UT 01 Jan 2000 */
   t  = gpssec - GPS_00h_01_Jan_2000;
   t += 1e-9 * gpsnan;
   t += taiutc - TAIUTC_00h_01_Jan_2000;

   /* compute number of days since 12h UT 01 Jan 2000 */
   dpU  = floor( t / ( 24.0 * 3600.0 ) ); /* full days since 0h UT 01
Jan 2000 */
   dpU += JD_00h_01_Jan_2000 - JD_12h_01_Jan_2000; /* i.e., -0.5 */

   /* compute number of centuries since 12h UT 31 Dec 1899 */
   TpU = dpU / 36525.0;

   /* compute the gmst at 0h of the current day */
   gmst = 24110.54841
     + TpU * ( 8640184.812866
         + TpU * ( 0.093104
           - TpU * 6.2e-6 ) ); /* seconds */

   /* add the sidereal time since the start of the day */
   t = fmod( t, 24.0 * 3600.0 ); /* seconds since start of day */
   gmst += t * 1.002737909350795; /* corrections omitted */

   /* convert to fractions of a day and to radians */
   gmst = fmod( gmst / ( 24.0 * 3600.0 ), 1.0 ); /* fraction of day */
   /* gmst *= 2.0 * M_PI; */ /* radians */ /* 10/15/04 gam */
   gmst *= ((REAL8)LALIND_TWOPI); /* radians */

   return gmst;
}
