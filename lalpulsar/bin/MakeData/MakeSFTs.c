/*
*  Copyright (C) 2007 Gregory Mendell
*  Copyright (C) 2010,2011,2016 Bernd Machenschalk
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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/

/**
 * \file
 * \ingroup lalpulsar_bin_SFTTools
 * \author Gregory Mendell, Xavier Siemens, Bruce Allen, Bernd Machenschalk
 * \brief generate SFTs
 */

/*********************************************************************************/
/*                                                                               */
/* File: MakeSFTs.c                                                              */
/* Purpose: generate SFTs                                                        */
/* Origin: first written by Xavier Siemens, UWM - May 2005,                      */
/*         based on make_sfts.c by Bruce Allen                                   */
/*********************************************************************************/

/* REVISIONS: */
/* 11/02/05 gam; To save memory, change lalDebugLevel to 0; note when this is set to 3 that 3x the memory is used! */
/* 11/02/05 gam; To save memory, do not save the window function in memory; just recompute this for each SFT. */
/* 11/02/05 gam; To save memory, do not hold dataSingle.data and vtilde in memory at the same time. */
/* 11/03/05 gam; Add TRACKMEMUSE preprocessor flag and function for tracking memory usage; copied from make_sfts.c */
/* 11/19/05 gam; Reorganize code so that function appear in order called */
/* 11/19/05 gam; Add command line option to use single rather than double precision; will use LALDButterworthREAL4TimeSeries in single precision case. */
/* 11/19/05 gam; Rename vtilde fftDataDouble; add in fftDatasingle; rename fplan fftPlanDouble; add in fftPlanSingle */
/* 11/29/05 gam; Add PRINTEXAMPLEDATA to print the first NUMTOPRINT, middle NUMTOPRINT, and last NUMTOPRINT input/ouput data at various stages */
/* 12/27/05 gam; Add option make-gps-dirs, -D <num>, to make directory based on this many GPS digits. */
/* 12/28/05 gam; Add option misc-desc, -X <string> giving misc. part of the SFT description field in the filename */
/* 12/28/05 gam; Make FMIN = 48.0 and DF = 2000.0 global variables, add start-freq -F and band -B options to enter these */
/* 12/28/05 gam; Add in window-type options; 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window */
/* 12/28/05 gam; Add option --overlap-fraction -P (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 0.0) */
/* 12/28/05 gam; Add sft-version, -v option to select output SFT version (1 = default is version 1 SFTs; 2 = version 2 SFTs */
/* 12/28/05 gam; Add comment-field, -c option, for comment for version 2 SFTs */
/* 01/05/06 gam; Add in version 2 normalization; add function to print example version 2 SFT data; add memory checking */
/* 01/09/06 gam; Add make-tmp-file, -Z option; write SFT to .*.tmp file, then move to final file name. */
/* 01/10/07 gam; Add -u --frame-struct-type option; specified the input data type in the frames (default ADC_REAL4) */
/* 01/14/07 gam; Add -i --ifo option to specify the ifo independent of the channel name which can begin with H0, L0, or G0. */
/* 06/26/07 gam; Write all command line arguments to commentField of version 2 SFTs, based on lalapps/src/calibration/ComputeStrainDriver.c */
/* 06/26/07 gam; Use finite to check that data does not contains a non-FINITE (+/- Inf, NaN) values, based on sftlib/SFTvalidate.c */
/* 10/05/12 gam; Add to version 2 normalization one over the root mean square of the window function (defined here as winFncRMS) as per RedMine LALSuite CW Bug #560*/
/* 24/07/14 eag; Change default SFT output to version 2 per RedMine LALSuite CW patch #1518 */

#include "config.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <glob.h>
#include <errno.h>
#include <stdarg.h>

#include <lal/LALDatatypes.h>
#include <lal/LALgetopt.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/AVFactories.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/Window.h>
#include <lal/Calibration.h>
#include <lal/LALConstants.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/SFTfileIO.h>
#include <lal/LALVCSInfo.h>
#include <lal/LALPulsarVCSInfo.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

/***************************************************************************/

/* STRUCTURES */
struct {
  REAL8 HPf;              /* High pass filtering frequency */
  INT4 T;                 /* SFT duration */
  char *stringT;          /* 12/27/05 gam; string with SFT duration */
  INT4 GPSStart;
  INT4 GPSEnd;
  char *commentField;      /* 12/28/05 gam; string comment for version 2 SFTs */
  char *FrCacheFile;       /* Frame cache file */
  char *ChannelName;
  char *SFTpath;           /* path to SFT file location */
  INT4 windowOption;       /* 12/28/05 gam; window options; 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window */
  REAL8 windowR;
  REAL8 overlapFraction;   /* 12/28/05 gam; overlap fraction (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 1.0). */
} CLA;


/***************************************************************************/

/* GLOBAL VARIABLES */
extern REAL8 winFncRMS;
extern REAL8TimeSeries dataDouble;
extern REAL4TimeSeries dataSingle;

REAL8 FMIN = 48.0; /* default start frequency */
REAL8 DF = 2000.0; /* default band */

static LALStatus status;
LALCache *framecache;         /* frame reading variables */
LALFrStream *framestream = NULL;

INT4 SegmentDuration;
LIGOTimeGPS gpsepoch;

REAL8FFTPlan *fftPlanDouble;           /* fft plan and data container, double precision case */
COMPLEX16Vector *fftDataDouble = NULL;

CHAR allargs[16384]; /* 06/26/07 gam; copy all command line args into commentField, based on lalapps/src/calibration/ComputeStrainDriver.c */
/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* Reads the command line */
int ReadCommandLine( int argc, char *argv[] );

/* Windows data */
int WindowData( REAL8 r );
int WindowDataTukey2( void );
int WindowDataHann( void );

/* prototypes */
void getSFTDescField( CHAR *sftDescField, CHAR *numSFTs, CHAR *ifo, CHAR *stringT, CHAR *typeMisc );
void mkSFTFilename( CHAR *sftFilename, CHAR *site, CHAR *numSFTs, CHAR *ifo, CHAR *stringT, CHAR *typeMisc, CHAR *gpstime );
void mvFilenames( CHAR *filename1, CHAR *filename2 );


/* -------------------- function definitions -------------------- */


/* 12/27/05 gam; function to create sftDescField for output directory or filename. */
void getSFTDescField( CHAR *sftDescField, CHAR *numSFTs, CHAR *ifo, CHAR *stringT, CHAR *typeMisc )
{
  strcpy( sftDescField, numSFTs );
  strcat( sftDescField, "_" );
  strcat( sftDescField, ifo );
  strcat( sftDescField, "_" );
  strcat( sftDescField, stringT );
  strcat( sftDescField, "SFT" );
  if ( typeMisc != NULL ) {
    strcat( sftDescField, "_" );
    strcat( sftDescField, typeMisc );
  }
}

/* 12/27/05 gam; make SFT file name according to LIGO T040164-01 specification */
void mkSFTFilename( CHAR *sftFilename, CHAR *site, CHAR *numSFTs, CHAR *ifo, CHAR *stringT, CHAR *typeMisc, CHAR *gpstime )
{
  CHAR sftDescField[256];
  strcpy( sftFilename, site );
  strcat( sftFilename, "-" );
  getSFTDescField( sftDescField, numSFTs, ifo, stringT, typeMisc );
  strcat( sftFilename, sftDescField );
  strcat( sftFilename, "-" );
  strcat( sftFilename, gpstime );
  strcat( sftFilename, "-" );
  strcat( sftFilename, stringT );
  strcat( sftFilename, ".sft" );
}

/* 01/09/06 gam; move filename1 to filename2 */
void mvFilenames( CHAR *filename1, CHAR *filename2 )
{
  CHAR mvFilenamesCommand[512 + 4];
  sprintf( mvFilenamesCommand, "mv %s %s", filename1, filename2 );
  if ( system( mvFilenamesCommand ) ) {
    XLALPrintError( "system() returned non-zero status\n" );
  }
}



/************************************* MAIN PROGRAM *************************************/

int main( int argc, char *argv[] )
{

  if ( ReadCommandLine( argc, argv ) ) {
    return 1;
  }
  SegmentDuration = CLA.GPSEnd - CLA.GPSStart ;

  /* create Frame cache, open frame stream and delete frame cache */
  framecache = XLALCacheImport( CLA.FrCacheFile );
  LALFrCacheOpen( &status, &framestream, framecache );
  TESTSTATUS( &status );
  XLALDestroyCache( framecache );

  if ( SegmentDuration < CLA.T ) {
    fprintf( stderr, "Cannot fit an SFT of duration %d between %d and %d\n",
             CLA.T, CLA.GPSStart, CLA.GPSEnd );
    return 0;;
  }

  gpsepoch.gpsSeconds = CLA.GPSStart;
  gpsepoch.gpsNanoSeconds = 0;

  /* Allocates space for data */
  {
    static FrChanIn chanin;

    chanin.name  = CLA.ChannelName;

    /* These calls just return deltaT for the channel */
    chanin.type  = ProcDataChannel;
    /* Get channel time step size by calling LALFrGetREAL8TimeSeries */
    LALFrSeek( &status, &gpsepoch, framestream );
    TESTSTATUS( &status );
    LALFrGetREAL8TimeSeries( &status, &dataDouble, &chanin, framestream );
    TESTSTATUS( &status );
    dataSingle.deltaT = dataDouble.deltaT;

    LALDCreateVector( &status, &dataDouble.data, ( UINT4 )( CLA.T / dataDouble.deltaT + 0.5 ) );
    TESTSTATUS( &status );

    fftPlanDouble = XLALCreateForwardREAL8FFTPlan( dataDouble.data->length, 0 );
    XLAL_CHECK( fftPlanDouble != NULL, XLAL_EFUNC );

  }

  while ( gpsepoch.gpsSeconds + CLA.T <= CLA.GPSEnd ) {

    /* Reads T seconds of data */
    {
      static FrChanIn chanin;
      chanin.name  = CLA.ChannelName;
      LALFrSeek( &status, &gpsepoch, framestream );
      TESTSTATUS( &status );

      chanin.type  = ProcDataChannel;
      LALFrGetREAL8TimeSeries( &status, &dataDouble, &chanin, framestream );
      TESTSTATUS( &status );
    }

    /* High-pass data with Butterworth filter */
    {
      PassBandParamStruc filterpar;
      char tmpname[] = "Butterworth High Pass";

      filterpar.name  = tmpname;
      filterpar.nMax  = 10;
      filterpar.f2    = CLA.HPf;
      filterpar.a2    = 0.5;
      filterpar.f1    = -1.0;
      filterpar.a1    = -1.0;

      if ( CLA.HPf > 0.0 ) {

        /* High pass the time series */
        LALButterworthREAL8TimeSeries( &status, &dataDouble, &filterpar );
        TESTSTATUS( &status );

      }

    }

    /* Window data; 12/28/05 gam; add options */
    if ( CLA.windowOption == 1 ) {
      if ( WindowData( CLA.windowR ) ) {
        return 5;  /* CLA.windowOption==1 is the default */
      }
    } else if ( CLA.windowOption == 2 ) {
      if ( WindowDataTukey2( ) ) {
        return 5;
      }
    } else if ( CLA.windowOption == 3 ) {
      if ( WindowDataHann( ) ) {
        return 5;
      }
    } else {
      /* Continue with no windowing; parsing of command line args makes sure options are one of the above or 0 for now windowing. */
    }

    /* create an SFT */
    /* 11/02/05 gam; allocate container for SFT data */
    LALZCreateVector( &status, &fftDataDouble, dataDouble.data->length / 2 + 1 );
    TESTSTATUS( &status );

    /* compute sft */
    XLAL_CHECK( XLALREAL8ForwardFFT( fftDataDouble, dataDouble.data, fftPlanDouble ) == XLAL_SUCCESS, XLAL_EFUNC );

    /* write out sft */
    {
      char sftname[256];
      char sftFilename[256];
      char sftnameFinal[256]; /* 01/09/06 gam */
      char numSFTs[2]; /* 12/27/05 gam */
      char site[2];    /* 12/27/05 gam */
      char ifo[3];     /* 12/27/05 gam; allow 3rd charactor for null termination */
      int firstbin = ( INT4 )( FMIN * CLA.T + 0.5 ), k;
      char gpstime[11]; /* 12/27/05 gam; allow for 10 digit GPS times and null termination */
      SFTtype *oneSFT = NULL;
      INT4 nBins = ( INT4 )( DF * CLA.T + 0.5 );
      REAL8 doubleDeltaT = 0.0; /* 01/05/06 gam */

      /* 12/27/05 gam; set up the number of SFTs, site, and ifo as null terminated strings */
      numSFTs[0] = '1';
      numSFTs[1] = '\0'; /* null terminate */
      strncpy( site, CLA.ChannelName, 1 );
      site[1] = '\0'; /* null terminate */
      strncpy( ifo, CLA.ChannelName, 2 );
      ifo[2] = '\0'; /* null terminate */
      sprintf( gpstime, "%09d", gpsepoch.gpsSeconds );

      strcpy( sftname, CLA.SFTpath );

      strcat( sftname, "/" );
      mkSFTFilename( sftFilename, site, numSFTs, ifo, CLA.stringT, NULL, gpstime );
      /* 01/09/06 gam; sftname will be temporary; will move to sftnameFinal. */
      /* set up sftnameFinal with usual SFT name */
      strcpy( sftnameFinal, sftname );
      strcat( sftnameFinal, sftFilename );
      /* sftname begins with . and ends in .tmp */
      strcat( sftname, "." );
      strcat( sftname, sftFilename );
      strcat( sftname, ".tmp" );

      /* make container to store the SFT data */
      XLAL_CHECK( ( oneSFT = XLALCreateSFT( ( ( UINT4 )nBins ) ) ) != NULL, XLAL_EFUNC );

      /* copy the data to oneSFT */
      strcpy( oneSFT->name, ifo );
      oneSFT->epoch.gpsSeconds = gpsepoch.gpsSeconds;
      oneSFT->epoch.gpsNanoSeconds = gpsepoch.gpsNanoSeconds;
      oneSFT->f0 = FMIN;
      oneSFT->deltaF = 1.0 / ( ( REAL8 )CLA.T );
      oneSFT->data->length = nBins;

      doubleDeltaT = ( REAL8 )( dataDouble.deltaT / winFncRMS ); /* include 1 over window function RMS */
      for ( k = 0; k < nBins; k++ ) {
        oneSFT->data->data[k] = crectf( doubleDeltaT * creal( fftDataDouble->data[k + firstbin] ), doubleDeltaT * cimag( fftDataDouble->data[k + firstbin] ) );
        /* 06/26/07 gam; use finite to check that data does not contains a non-FINITE (+/- Inf, NaN) values */
        if ( !isfinite( crealf( oneSFT->data->data[k] ) ) || !isfinite( cimagf( oneSFT->data->data[k] ) ) ) {
          fprintf( stderr, "Infinite or NaN data at freq bin %d.\n", k );
          return 7;
        }
      }
      LALZDestroyVector( &status, &fftDataDouble );
      TESTSTATUS( &status );

      /* write the SFT */
      XLAL_CHECK( XLALWriteSFT2NamedFile(oneSFT, sftname, "unknown" /* FIXME */, 0, CLA.commentField) == XLAL_SUCCESS, XLAL_EFUNC );

      /* 01/09/06 gam; sftname is temporary; move to sftnameFinal. */
      mvFilenames( sftname, sftnameFinal );

      XLALDestroySFT( oneSFT );
    }

    gpsepoch.gpsSeconds = gpsepoch.gpsSeconds + ( INT4 )( ( 1.0 - CLA.overlapFraction ) * ( ( REAL8 )CLA.T ) );
    gpsepoch.gpsNanoSeconds = 0;
  }

  LALFrClose( &status, &framestream );
  TESTSTATUS( &status );

  LALDDestroyVector( &status, &dataDouble.data );
  TESTSTATUS( &status );
  XLALDestroyREAL8FFTPlan( fftPlanDouble );

  LALCheckMemoryLeaks();

  return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/

/*** FUNCTIONS ***/

/*******************************************************************************/
int ReadCommandLine( int argc, char *argv[] )
{
  INT4 errflg = 0;
  INT4 i;              /* 06/26/07 gam */
  struct LALoption long_options[] = {
    {"high-pass-freq",       required_argument, NULL,          'f'},
    {"sft-duration",         required_argument, NULL,          't'},
    {"sft-write-path",       required_argument, NULL,          'p'},
    {"frame-cache",          required_argument, NULL,          'C'},
    {"channel-name",         required_argument, NULL,          'N'},
    {"gps-start-time",       required_argument, NULL,          's'},
    {"gps-end-time",         required_argument, NULL,          'e'},
    /* {"sft-version",          required_argument, NULL,          'v'}, */
    {"comment-field",        required_argument, NULL,          'c'},
    {"start-freq",           required_argument, NULL,          'F'},
    {"band",                 required_argument, NULL,          'B'},
    /* {"make-gps-dirs",        required_argument, NULL,          'D'}, */
    /* {"make-tmp-file",        required_argument, NULL,          'Z'}, */
    /* {"misc-desc",            required_argument, NULL,          'X'}, */
    /* {"frame-struct-type",    required_argument, NULL,          'u'}, */
    /* {"ifo",                  required_argument, NULL,          'i'}, */
    {"window-type",          required_argument, NULL,          'w'},
    {"window-radius",        required_argument, NULL,          'r'},
    {"overlap-fraction",     required_argument, NULL,          'P'},
    /* {"ht-data",              no_argument,       NULL,          'H'}, */
    /* {"use-single",           no_argument,       NULL,          'S'}, */
    {"help",                 no_argument,       NULL,          'h'},
    {"version",              no_argument,       NULL,          'V'},
    {0, 0, 0, 0}
  };
  char args[] = "hHZSf:t:C:N:i:s:e:v:c:F:B:D:X:u:w:r:P:p:ab:";

  /* Initialize default values */
  CLA.HPf = -1.0;
  CLA.T = 0.0;
  CLA.stringT = NULL; /* 12/27/05 gam */
  CLA.FrCacheFile = NULL;
  CLA.GPSStart = 0;
  CLA.GPSEnd = 0;
  CLA.ChannelName = NULL;
  CLA.SFTpath = NULL;
  CLA.commentField = NULL; /* 12/28/05 gam; comment for version 2 SFT header. */
  CLA.windowOption = 1; /* 12/28/05 gam; window options; 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window */
  CLA.windowR = 0.001;
  CLA.overlapFraction = 0.0; /* 12/28/05 gam; overlap fraction (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 0.0). */

  strcat( allargs, "\nMakeSFTs " );
  strcat( allargs, lalVCSIdentInfo.vcsId );
  strcat( allargs, lalVCSIdentInfo.vcsStatus );
  strcat( allargs, "\nMakeSFTs " );
  strcat( allargs, lalPulsarVCSIdentInfo.vcsId );
  strcat( allargs, lalPulsarVCSIdentInfo.vcsStatus );
  strcat( allargs, "\nMakeSFTs command line args: " ); /* 06/26/07 gam; copy all command line args into commentField */
  for ( i = 0; i < argc; i++ ) {
    /* if the argument describes an additioanl comment, don't record it
       in the command-line as well, as it will be recoeded in the comment later.
       if the option is a single argument and does not include the comment itself,
       skip the next argument as well, as this should be the comment then. */
    if ( ( strstr( argv[i], "-c" ) == argv[i] ) ) {
      strcat( allargs, "-c ... " );
      if ( strcmp( argv[i], "-c" ) == 0 ) {
        i++;
      }
      continue;
    }
    if ( ( strstr( argv[i], "--comment-field" ) == argv[i] ) ) {
      strcat( allargs, "--comment-field ... " );
      if ( strcmp( argv[i], "--comment-field" ) == 0 ) {
        i++;
      }
      continue;
    }
    strcat( allargs, argv[i] );
    strcat( allargs, " " );
  }
  strcat( allargs, "\n" );
  CLA.commentField = allargs;

  /* Scan through list of command line arguments */
  while ( 1 ) {
    int option_index = 0; /* LALgetopt_long stores long option here */
    int c;

    c = LALgetopt_long_only( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) { /* end of options */
      break;
    }

    switch ( c ) {
    case 'f':
      /* high pass frequency */
      CLA.HPf = atof( LALoptarg );
      break;
    case 't':
      /* SFT time */
      CLA.stringT = LALoptarg;  /* 12/27/05 gam; keep pointer to string that gives the SFT duration */
      CLA.T = atoi( LALoptarg );
      break;
    case 'C':
      /* name of frame cache file */
      CLA.FrCacheFile = LALoptarg;
      break;
    case 's':
      /* GPS start */
      CLA.GPSStart = atof( LALoptarg );
      break;
    case 'e':
      /* GPS end */
      CLA.GPSEnd = atof( LALoptarg );
      break;
    case 'F':
      /* 12/28/05 gam; start frequency */
      FMIN = ( REAL8 )atof( LALoptarg );
      break;
    case 'B':
      /* 12/28/05 gam; band */
      DF = ( REAL8 )atof( LALoptarg );
      break;
    case 'c':
      /* 12/28/05 gam; comment for version 2 SFTs */
      strcat( CLA.commentField, "MakeSFTs additional comment: " ); /* 06/26/07 gam; copy all command line args into commentField */
      strcat( CLA.commentField, LALoptarg );
      strcat( CLA.commentField, "\n" );
      break;
    case 'w':
      /* 12/28/05 gam; window options; 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window */
      CLA.windowOption = atoi( LALoptarg );
      break;
    case 'r':
      /* defulat 0.001 */
      CLA.windowR = ( REAL8 )atof( LALoptarg );
      break;
    case 'P':
      /* 12/28/05 gam; overlap fraction (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 1.0). */
      CLA.overlapFraction = ( REAL8 )atof( LALoptarg );
      break;
    case 'N':
      CLA.ChannelName = LALoptarg;
      break;
    case 'p':
      CLA.SFTpath = LALoptarg;
      break;
    case 'h':
      /* print usage/help message */
      fprintf( stdout, "Arguments are:\n" );
      fprintf( stdout, "\thigh-pass-freq (-f)\tFLOAT\t High pass filtering frequency in Hz.\n" );
      fprintf( stdout, "\tsft-duration (-t)\tFLOAT\t SFT duration in seconds.\n" );
      fprintf( stdout, "\tsft-write-path (-p)\tFLOAT\t Location of output SFTs.\n" );
      fprintf( stdout, "\tframe-cache (-C)\tSTRING\t Path to frame cache file (including the filename).\n" );
      fprintf( stdout, "\tgps-start-time (-s)\tINT\t GPS start time of segment.\n" );
      fprintf( stdout, "\tgps-end-time (-e)\tINT\t GPS end time of segment.\n" );
      fprintf( stdout, "\tchannel-name (-N)\tSTRING\t Name of channel to read within a frame.\n" );
      fprintf( stdout, "\tifo (-i)\t\tSTRING\t (optional) Name of IFO, i.e., H1, H2, L1, or G1; use if channel name begins with H0, L0, or G0; default: use first two characters from channel name.\n" );
      fprintf( stdout, "\tcomment-field (-c)\tSTRING\t (optional) Comment for version 2 SFT header.\n" );
      fprintf( stdout, "\tstart-freq (-F) \tFLOAT\t (optional) Start frequency of the SFTs (default is 48 Hz).\n" );
      fprintf( stdout, "\tband (-B)       \tFLOAT\t (optional) Frequency band of the SFTs (default is 2000 Hz).\n" );
      fprintf( stdout, "\tmake-gps-dirs (-D)\tINT\t (optional) Make directories for output SFTs based on this many digits of the GPS time.\n" );
      fprintf( stdout, "\tmake-tmp-file (-Z)\tINT\t (optional) Write SFT to .*.tmp file, then move to final filename.\n" );
      fprintf( stdout, "\twindow-type (-w)\tINT\t (optional) 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window\n" );
      fprintf( stdout, "\twindow-radius (-r)\tFLOAT\t (optional) default = 0.001\n" );
      fprintf( stdout, "\toverlap-fraction (-P)\tFLOAT\t (optional) Overlap fraction (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 0.0).\n" );
      fprintf( stdout, "\tversion (-V)\t\tFLAG\t Print LAL & LALPulsar version and exit.\n" );
      fprintf( stdout, "\thelp (-h)\t\tFLAG\t This message.\n" );
      exit( 0 );
    case 'V':
      /* print version */
      fprintf( stdout, "MakeSFTs %s %s\n", lalVCSIdentInfo.vcsId, lalVCSIdentInfo.vcsStatus );
      fprintf( stdout, "MakeSFTs %s %s\n", lalPulsarVCSIdentInfo.vcsId, lalPulsarVCSIdentInfo.vcsStatus );
      exit( 0 );
    default:
      /* unrecognized option */
      errflg++;
      if ( ( c >= 48 ) && ( c < 128 ) ) {
        fprintf( stderr, "Unrecognized option '%c'\n", c );
      } else {
        fprintf( stderr, "Unrecognized option %d\n", c );
      }
      exit( 1 );
      break;
    }
  }

  if ( CLA.HPf < 0 ) {
    fprintf( stderr, "No high pass filtering frequency specified.\n" );
    fprintf( stderr, "If you don't want to high pass filter set the frequency to 0.\n" );
    fprintf( stderr, "Try %s -h \n", argv[0] );
    return 1;
  }
  if ( CLA.T == 0.0 ) {
    fprintf( stderr, "No SFT duration specified.\n" );
    fprintf( stderr, "Try %s -h \n", argv[0] );
    return 1;
  }
  if ( CLA.GPSStart == 0 ) {
    fprintf( stderr, "No GPS start time specified.\n" );
    fprintf( stderr, "Try %s -h \n", argv[0] );
    return 1;
  }
  if ( CLA.GPSEnd == 0 ) {
    fprintf( stderr, "No GPS end time specified.\n" );
    fprintf( stderr, "Try %s -h \n", argv[0] );
    return 1;
  }
  if ( FMIN < 0.0 ) {
    fprintf( stderr, "Illegal start-freq option given.\n" );
    fprintf( stderr, "Try %s -h \n", argv[0] );
    return 1;
  }
  if ( DF < 0.0 ) {
    fprintf( stderr, "Illegal band option given.\n" );
    fprintf( stderr, "Try %s -h \n", argv[0] );
    return 1;
  }
  if ( ( CLA.windowOption < 0 ) || ( CLA.windowOption > 3 ) ) {
    fprintf( stderr, "Illegal window-type given.\n" );
    fprintf( stderr, "Try %s -h \n", argv[0] );
    return 1;
  }
  if ( ( CLA.overlapFraction < 0.0 ) || ( CLA.overlapFraction >= 1.0 ) ) {
    fprintf( stderr, "Illegal overlap-fraction given.\n" );
    fprintf( stderr, "Try %s -h \n", argv[0] );
    return 1;
  }
  if ( CLA.FrCacheFile == NULL ) {
    fprintf( stderr, "No frame cache file specified.\n" );
    fprintf( stderr, "Try %s -h \n", argv[0] );
    return 1;
  }
  if ( CLA.ChannelName == NULL ) {
    fprintf( stderr, "No data channel name specified.\n" );
    fprintf( stderr, "Try %s -h \n", argv[0] );
    return 1;
  }
  if ( CLA.SFTpath == NULL ) {
    fprintf( stderr, "No output path specified for SFTs.\n" );
    fprintf( stderr, "Try %s -h \n", argv[0] );
    return 1;
  }

  return errflg;
}
/*******************************************************************************/
