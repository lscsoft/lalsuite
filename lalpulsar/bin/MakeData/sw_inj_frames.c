/*
 * Copyright (C) 2010 Erin Macdonald, Matthew Pitkin
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

/* Code to create frames and add them to existing data files
Input format as $ ./sw_inj_frames framefile duration epoch
example$ ./lalpulsar_sw_inj_frames -p /Users/erinmacdonald/lsc/analyses/test_par_files -g /Users/erinmacdonald/lsc/analyses/frames -c /Users/erinmacdonald/lsc/analyses/CWINJframes -I H1 -e /Users/erinmacdonald/lsc/src/lscsoft/lalsuite/lalpulsar/test -y 09-11
*/

/* 2/3/11 v. 2 (which is not entirely accurate, but I'm starting from here): added a log file, version number, and separate channel for i(t) -- E. Macdonald */
/* 3/3/11 v. 3: placed log comments throughout */
/* 7/3/11 v. 4: output file to dump ascii of failed frames */
/* 21/3/11 v. 5: fixed closing .par file L520 */
/* 9/4/11 v.6: including proper motion calculations */
/* 12/5/11 v.7: fix memory leak */
/* 31/5/11 v.8: Matt Pitkin's fixes for the memory leak*/
/* 23/6/11 v.9; Changed naming of log files to append time, sloppy but works */
/* 29/6/11 v.10; XLAL functions memory leaks fixed */
/* 11/7/11 v.11; Outputs surpressed and better log output */
/* 30/1/12 v.12: Set signal generation barycenter delay lookup table step size */
/* 27/06/13 v.13: Allow 2nd and 3rd frequency derivatives and make compatible with new frame functions */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <errno.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <signal.h>
#include <math.h>
#include <dirent.h>
#include <sys/stat.h>

#include <lal/Units.h>
#include <lal/LALFrStream.h>
#include <lal/LALFrameIO.h>
#include <lal/LALCache.h>
#include <lal/TimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/UserInput.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALDatatypes.h>
#include <lal/DetectorSite.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/RngMedBias.h>
#include <lal/LALRunningMedian.h>
#include <lal/BinaryPulsarTiming.h>
#include <lal/LogPrintf.h>
#include <lal/LALString.h>
#include <lal/LALPulsarVCSInfo.h>

#define STRINGLENGTH 256              /* the length of general string */

/* ---------- local type definitions ---------- */
typedef struct {
  /*    CHAR *out_chan;  Channel output ie. H1_LDAS_C02_L2_CWINJ*/
  /*    CHAR *in_chan;  Channel input from .gwf ie. H1:LDAS-STRAIN*/
  REAL8 srate; /* sample rate -- default 16384*/
  /*  REAL8 duration; duration (sec) Now read from file name */
  /*  REAL8 start; epoch (GPSSeconds) Now read from file name */
  CHAR *inputdir; /* directory for .par files*/
  CHAR *gwfdir; /*directory for .gwf files*/
  CHAR *ephemEarth;             /**< Earth ephemeris file to use */
  CHAR *ephemSun;               /**< Sun ephemeris file to use */
  CHAR *IFO; /*detector */
  CHAR *logDir; /*directory for the log files */
  CHAR *outputdir;
} UserInput_t;

UserInput_t uvar_struct;

/* ---------- local function prototypes ---------- */
int InitUserVars( UserInput_t *uvar, int argc, char **argv );  /*Initiates user variables*/

int XLALFrameFileName( char *fname, size_t size, const char *chname, const LIGOTimeGPS *epoch, double duration );

/* ---------- Function definitions ---------- */

int main( int argc, char **argv )
{
  const char *fn = __func__;

  UserInput_t *uvar = &uvar_struct; /*structure for user-defined input variables*/

  /* read all user input */
  if ( InitUserVars( uvar, argc, argv ) != XLAL_SUCCESS ) {
    XLAL_ERROR( XLAL_EFUNC );
  }

  char version[256];
  sprintf( version, "v13" ); /*manually change */

  /*Get current time for log */
  time_t result;
  result = time( NULL );
  struct tm *btm = localtime( &result );

  /*Start the log file */
  FILE *logfile;
  char log_file[256];

  char logfiledir[6];
  strncpy( logfiledir, ( strrchr( uvar->gwfdir, '/' ) + 1 ), 5 );
  logfiledir[5] = '\0';

  sprintf( log_file, "%s/%s/injections-%s.log", uvar->outputdir, uvar->logDir, logfiledir );
  /*fprintf(stderr, "Your log file is here: %s\n", log_file);*/

  /*Writing to .log file*/
  if ( ( logfile = fopen( log_file, "a" ) ) == NULL ) {
    fprintf( stderr, "Error opening .log file! \n" );
    return 0;
  } else {
    fprintf( logfile, "\nsw_inj_frames.c version number: %s\nCurrent time: %s \n", version, asctime( btm ) );
    fprintf( logfile, "User inputs:\n Sample rate: %f\n Pulsar file directory: %s\n Raw frame file directory: %s\n Output directory: %s\n \n Name of directory for log (w/in output dir): %s\n IFO: %s\n", uvar->srate, uvar->inputdir, uvar->gwfdir, uvar->outputdir, uvar->logDir, uvar->IFO );
    fclose( logfile );
  }
  /* </log file> */

  struct dirent **parnamelist;
  struct dirent **gwfnamelist;
  char pulin[256];

  LALStatus XLAL_INIT_DECL( status );

  /*init ephemeris-data */
  EphemerisData *edat;
  XLAL_CHECK( ( edat = XLALInitBarycenter( uvar->ephemEarth, uvar->ephemSun ) ) != NULL, XLAL_EFUNC );

  /*init detector info */
  LALDetector *site;
  if ( ( site = XLALGetSiteInfo( uvar -> IFO ) ) == NULL ) {
    XLALPrintError( "%s: Failed to get site-info for detector '%s'\n", fn, uvar->IFO );
    XLAL_ERROR( XLAL_EFUNC );
  }


  /*  if ( (params.pulsar.spindown = XLALCreateREAL8Vector(1))== NULL ) { */
  /*     XLALPrintError("Out of memory"); */
  /*     XLAL_ERROR ( fn, XLAL_EFUNC ); */
  /*   } */

  REAL8 srate; /*sample rate defaulted to 16384 */
  srate = uvar->srate;

  /*creates the .gwf channels */
  CHAR in_chan[256];
  sprintf( in_chan, "%s:LDAS-STRAIN", uvar->IFO );
  /*fprintf( stderr, "\nin_chan = %s\n", in_chan);*/

  CHAR out_chan[256];
  sprintf( out_chan, "%s:CWINJ_TOT", uvar->IFO );
  /*fprintf( stderr, "\nout_chan = %s\n", out_chan);*/

  CHAR inj_chan[256];
  sprintf( inj_chan, "%s:CWINJ_SIG", uvar->IFO );

  /*Writing to .log file*/
  if ( ( logfile = fopen( log_file, "a" ) ) == NULL ) {
    fprintf( stderr, "Error opening .log file! \n" );
    return 0;
  } else {
    fprintf( logfile, "Channels:\n Reading from raw: %s\n Injection+noise: %s\n Pure injection: %s\n", in_chan, out_chan, inj_chan );
    fclose( logfile );
  }
  /* </log file> */

  /* Major loop to run through gwf files */
  int m;
  UINT4 i;

  char out_file[256]; /*need a different way to do this */
  FILE *inject; /*for .par file from uvar->inputdir*/

  CHAR gwf_dir[256]; /* for use with XLALFrStreamOpen */
  sprintf( gwf_dir, "%s/.", uvar->gwfdir );

  LIGOTimeGPS epoch;
  char strepoch[10];

  UINT4 ndata;

  /** Put the pulsar files in the log (so as not to loop every time) **/

  int n = scandir( uvar->inputdir, &parnamelist, 0, alphasort );
  if ( n < 0 ) {
    XLALPrintError( "%s: scandir() failed for directory '%s'\n", fn, uvar->inputdir );
    XLAL_ERROR( XLAL_EIO );
  }

  UINT4 numParFiles = ( UINT4 )n;
  UINT4 h = 0;

  char parname[256];
  PulsarParameters *pulparams[numParFiles];

  /*clear . and .. in directory */
  free( parnamelist[0] );
  free( parnamelist[1] );

  for ( h = 2; h < numParFiles; h++ ) {
    if ( strstr( parnamelist[h]->d_name, ".par" ) == NULL ) {

      free( parnamelist[h] );

      continue;
    } else {
      /*Writing to .log file*/
      if ( ( logfile = fopen( log_file, "a" ) ) == NULL ) {
        fprintf( stderr, "Error opening .log file! \n" );
        return 0;
      } else {
        sprintf( parname, "%s/%s", uvar->inputdir, parnamelist[h]->d_name );
        char t[ 100 ];
        struct stat b;
        if ( !stat( parname, &b ) ) {
          strftime( t, 100, "%d/%m/%Y %H:%M:%S", localtime( &b.st_mtime ) );
          fprintf( logfile, "%s/%s Last Modified: %s\n", uvar->inputdir, parnamelist[h]->d_name, t );
          fclose( logfile );
        }
      }
      /* </log file> */

      /* Save the pulsar parameter files to memory to be called later*/
      /*sprintf(pulin,"%s/%s", uvar->inputdir, parnamelist[h]->d_name);*/
      if ( ( inject = fopen( parname, "r" ) ) == NULL ) {
        fprintf( stderr, "Error opening file: %s\n", parname );
        XLAL_ERROR( XLAL_EIO );
      }

      pulparams[h] = XLALReadTEMPOParFile( parname );

      fclose( inject );
      /*fprintf(stderr,"You have now read in %s, saved in index %i", parname, h);*/

    } /*checking if file ends in .par*/

    free( parnamelist[h] );

  }/*looping through files in uvar->inputdir */

  free( parnamelist );

  /** Done logging .par files **/

  if ( ( m = scandir( uvar->gwfdir, &gwfnamelist, 0, alphasort ) ) == -1 ) {
    XLALPrintError( "%s: scandir('%s',...) failed.\n", fn, uvar->gwfdir );
    XLAL_ERROR( XLAL_EIO );
  }

  free( gwfnamelist[0] );
  free( gwfnamelist[1] );

  UINT4 k = 0;
  for ( k = 2; k < ( UINT4 )m; k++ ) {
    if ( strstr( gwfnamelist[k]->d_name, ".gwf" ) == NULL ) {

      free( gwfnamelist[k] );

      continue; /*make sure it's a .gwf file */
    } else {
      LALFrStream *gwffile = NULL;

      if ( ( gwffile = XLALFrStreamOpen( uvar->gwfdir, gwfnamelist[k]->d_name ) ) == NULL ) {
        /*XLAL_ERROR ( fn, XLAL_EFUNC ); -- don't want to abort, but save elsewhere test that it's an acceptable file */

        /*Writing to failed file*/
        FILE *frames;
        char framefilename[256];
        sprintf( framefilename, "%s/%s/failed_frames.txt", uvar->outputdir, uvar->logDir );
        if ( ( frames = fopen( framefilename, "a" ) ) == NULL ) {
          fprintf( stderr, "Error opening file %s! \n", framefilename );
          return 0;
        } else {
          fprintf( frames, "%s\n", gwfnamelist[k]->d_name );
          fclose( frames );
        }
        /* </failed file> */
        continue; /*don't exit program if .gwf file fails, continue through*/
      } else {
        /*Writing to .log file*/
        if ( ( logfile = fopen( log_file, "a" ) ) == NULL ) {
          fprintf( stderr, "Error opening .log file! \n" );
          return 0;
        } else {
          fprintf( logfile, "Read file %s/%s\n", uvar->gwfdir, gwfnamelist[k]->d_name );
          fclose( logfile );
        }
        /* </log file> */

        REAL8TimeSeries *gwfseries = NULL;
        REAL8TimeSeries *series = NULL;

        /** extract epoch and duration from gwf file name **/

        strncpy( strepoch, strchr( gwfnamelist[k]->d_name, '9' ), 9 ); /*All epochs in S6 begin with 9... potential problem in future */
        XLAL_LAST_ELEM( strepoch ) = '\0'; /*Null terminate the string*/
        epoch.gpsSeconds = atoi( strepoch ); /* convert to integer from string */
        epoch.gpsNanoSeconds = 0; /* no nanosecond precision */
        /*      fprintf(stderr, "epoch = %i\n", epoch.gpsSeconds);*/

        char strdur[4];
        strncpy( strdur, ( strrchr( gwfnamelist[k]->d_name, '-' ) + 1 ), 3 ); /* duration is last number in frame file */
        XLAL_LAST_ELEM( strdur ) = '\0';
        /* assigns duration from .gwf frame */
        ndata = atoi( strdur );

        /*Writing to .log file*/
        if ( ( logfile = fopen( log_file, "a" ) ) == NULL ) {
          fprintf( stderr, "Error opening .log file! \n" );
          return 0;
        } else {
          fprintf( logfile, "Using epoch: %i and duration: %i\n", epoch.gpsSeconds, ndata );
          fclose( logfile );
        }
        /* </log file> */

        /* acquire time series from frame file */
        if ( ( gwfseries = XLALCreateREAL8TimeSeries( in_chan, &epoch, 0., 1. / srate, &lalSecondUnit, ( int )( ndata * srate ) ) ) == NULL ) {
          XLAL_ERROR( XLAL_EFUNC );
        }

        if ( XLALFrStreamGetREAL8TimeSeries( gwfseries, gwffile ) != XLAL_SUCCESS ) {
          /*Writing to failed file*/
          FILE *frames;
          char framefilename[256];
          sprintf( framefilename, "%s/%s/failed_frames.txt", uvar->outputdir, uvar->logDir );
          if ( ( frames = fopen( framefilename, "a" ) ) == NULL ) {
            fprintf( stderr, "Error opening file %s! \n", framefilename );
            return 0;
          } else {
            fprintf( frames, "%s\n", gwfnamelist[k]->d_name );
            fclose( frames );
          }

          /* </failed file> */
          XLALDestroyREAL8TimeSeries( gwfseries );
          continue; /*don't exit program if .gwf file fails, continue through*/
        }

        int detflags;

        /* define output .gwf file */
        detflags = XLALFrameFileName( out_file, sizeof( out_file ), out_chan, &epoch, ndata );
        //sprintf(out_file, "%c-%s-%d-%d.gwf", uvar->IFO[0], out_chan, epoch.gpsSeconds, ndata);

        /*Writing to .log file*/
        if ( ( logfile = fopen( log_file, "a" ) ) == NULL ) {
          fprintf( stderr, "Error opening .log file! \n" );
          return 0;
        } else {
          fprintf( logfile, "Writing to %s/%s\n", uvar->outputdir, out_file );
          fclose( logfile );
        }
        /* </log file> */

        /* read in and test generated frame with XLAL function*/

        LALFrameH *outFrame   = NULL;        /* frame data structure */

        if ( ( outFrame = XLALFrameNew( &epoch, ( REAL8 )ndata, "CW_INJ", 0, 0, detflags ) ) == NULL ) {
          LogPrintf( LOG_CRITICAL, "%s : XLALFrameNew() failed with error = %d.\n", fn, xlalErrno );
          XLAL_ERROR( XLAL_EFAILED );
        }

        series = NULL;
        if ( ( series = XLALCreateREAL8TimeSeries( out_chan, &epoch, 0., 1. / srate, &lalSecondUnit, ( int )( ndata * srate ) ) ) == NULL ) {
          XLAL_ERROR( XLAL_EFUNC );
        }

        FILE *filetest;
        CHAR fffile[256];
        snprintf( fffile, STRINGLENGTH, "%s/%s", uvar->outputdir, out_file );

        if ( ( filetest = fopen( fffile, "w+" ) ) == NULL ) {
          fprintf( stderr, "fail" );
        }

        fclose( filetest );

        /*create series to be the sum, general series to add to the noise */
        REAL8TimeSeries *total_inject = NULL;
        if ( ( total_inject = XLALCreateREAL8TimeSeries( inj_chan, &epoch, 0., 1. / srate, &lalSecondUnit, ( UINT4 )( ndata * srate ) ) ) == NULL ) {
          XLAL_ERROR( XLAL_EFUNC );
        }

        /* need to set all values in total_inject to zero so not pre-allocated */
        UINT4 counter = 0;
        for ( counter = 0; counter < total_inject->data->length; counter++ ) {
          total_inject->data->data[counter] = 0;
        }

        /*Read in all pulsar files from inputdir*/
        n = scandir( uvar->inputdir, &parnamelist, 0, alphasort );
        if ( n < 0 ) {
          XLALPrintError( "%s: scandir() failed for directory '%s'\n", fn, uvar->inputdir );
          XLAL_ERROR( XLAL_EIO );
        }
        numParFiles = ( UINT4 )n;
        h = 0;

        for ( h = 0; h < numParFiles; h++ ) {
          if ( strstr( parnamelist[h]->d_name, ".par" ) == NULL ) {
            /*fprintf(stderr, "Not using file %s\n",parnamelist[h]->d_name);*/

            free( parnamelist[h] );

            continue;
          } else {
            /*sprintf(pulin, "%s/%s", uvar->inputdir, parnamelist[h]->d_name);*/
            /*fprintf(stderr, "This is your pulsar file:%s\n", pulin);*/
            /*if (( inject = fopen ( pulin, "r" )) == NULL ){*/
            /*fprintf(stderr, "Error opening file: %s\n", pulin);*/
            /*XLAL_ERROR ( fn, XLAL_EIO );*/
            /*}*/
            PulsarSignalParams XLAL_INIT_DECL( params );
            /* set signal generation barycenter delay look-up table step size */
            params.dtDelayBy2 = 10.; /* generate table every 10 seconds */

            /*BinaryPulsarParams pulparams; read from the .par file */
            if ( ( params.pulsar.spindown = XLALCreateREAL8Vector( 3 ) ) == NULL ) {
              XLALPrintError( "Out of memory" );
              XLAL_ERROR( XLAL_EFUNC );
            }

            const REAL8Vector *fs = PulsarGetREAL8VectorParam( pulparams[h], "F" );
            fprintf( stderr, "Your f0 is %f\n", fs->data[0] );

            /*Convert location with proper motion */
            INT4 dtpos = 0;
            if ( PulsarCheckParam( pulparams[h], "POSEPOCH" ) ) {
              dtpos = epoch.gpsSeconds - ( INT4 )PulsarGetREAL8Param( pulparams[h], "POSEPOCH" );  /*calculate time since posepoch */
            } else {
              dtpos = epoch.gpsSeconds - ( INT4 )PulsarGetREAL8Param( pulparams[h], "PEPOCH" );  /*calculate time since pepoch */
            }

            REAL8 ra = 0., dec = 0.;
            if ( PulsarCheckParam( pulparams[h], "RAJ" ) ) {
              ra = PulsarGetREAL8Param( pulparams[h], "RAJ" );
            } else if ( PulsarCheckParam( pulparams[h], "RA" ) ) {
              ra = PulsarGetREAL8Param( pulparams[h], "RA" );
            } else {
              XLALPrintError( "No right ascension found" );
              XLAL_ERROR( XLAL_EFUNC );
            }
            if ( PulsarCheckParam( pulparams[h], "DECJ" ) ) {
              dec = PulsarGetREAL8Param( pulparams[h], "DECJ" );
            } else if ( PulsarCheckParam( pulparams[h], "DEC" ) ) {
              dec = PulsarGetREAL8Param( pulparams[h], "DEC" );
            } else {
              XLALPrintError( "No declination found" );
              XLAL_ERROR( XLAL_EFUNC );
            }

            params.pulsar.position.latitude = dec + ( REAL8 )dtpos * PulsarGetREAL8ParamOrZero( pulparams[h], "PMDEC" );
            params.pulsar.position.longitude = ra + ( REAL8 )dtpos * PulsarGetREAL8ParamOrZero( pulparams[h], "PMRA" ) / cos( params.pulsar.position.latitude );
            params.pulsar.position.system = COORDINATESYSTEM_EQUATORIAL;

            params.pulsar.f0 = 2.*fs->data[0];
            for ( UINT4 kk = 0; kk < 3; kk++ ) {
              // set frequency derivaties (up to third derivative)
              if ( fs->length > kk ) {
                params.pulsar.spindown->data[kk] = 2.*fs->data[kk + 1];
              } else {
                params.pulsar.spindown->data[kk] = 0.;
              }
            }
            if ( ( XLALGPSSetREAL8( &( params.pulsar.refTime ), PulsarGetREAL8Param( pulparams[h], "PEPOCH" ) ) ) == NULL ) {
              XLAL_ERROR( XLAL_EFUNC );
            }
            params.pulsar.psi = PulsarGetREAL8ParamOrZero( pulparams[h], "PSI" );
            params.pulsar.phi0 = PulsarGetREAL8ParamOrZero( pulparams[h], "PHI0" );
            /*Conversion from h0 and cosi to plus and cross */
            REAL8 cosiota = PulsarGetREAL8ParamOrZero( pulparams[h], "COSIOTA" );
            REAL8 h0 = PulsarGetREAL8ParamOrZero( pulparams[h], "H0" );
            params.pulsar.aPlus = 0.5 * h0 * ( 1. + cosiota * cosiota );
            params.pulsar.aCross = h0 * cosiota;

            /*if binary */
            if ( PulsarCheckParam( pulparams[h], "BINARY" ) ) {
              /*fprintf(stderr, "You have a binary! %s", pulin);*/
              params.orbit.asini = PulsarGetREAL8ParamOrZero( pulparams[h], "A1" );
              params.orbit.period = PulsarGetREAL8ParamOrZero( pulparams[h], "PB" ) * 86400;

              /*Taking into account ELL1 model option */
              if ( strstr( PulsarGetStringParam( pulparams[h], "BINARY" ), "ELL1" ) != NULL ) {
                REAL8 w, e, eps1, eps2;
                eps1 = PulsarGetREAL8ParamOrZero( pulparams[h], "EPS1" );
                eps2 = PulsarGetREAL8ParamOrZero( pulparams[h], "EPS2" );
                w = atan2( eps1, eps2 );
                e = sqrt( eps1 * eps1 + eps2 * eps2 );
                params.orbit.argp = w;
                params.orbit.ecc = e;
              } else {
                params.orbit.argp = PulsarGetREAL8ParamOrZero( pulparams[h], "OM" );
                params.orbit.ecc = PulsarGetREAL8ParamOrZero( pulparams[h], "ECC" );
              }
              if ( strstr( PulsarGetStringParam( pulparams[h], "BINARY" ), "ELL1" ) != NULL ) {
                REAL8 fe, uasc, Dt, T0val = 0.;
                fe = sqrt( ( 1.0 - params.orbit.ecc ) / ( 1.0 + params.orbit.ecc ) );
                uasc = 2.0 * atan( fe * tan( params.orbit.argp / 2.0 ) );
                Dt = ( params.orbit.period / LAL_TWOPI ) * ( uasc - params.orbit.ecc * sin( uasc ) );
                T0val = PulsarGetREAL8ParamOrZero( pulparams[h], "TASC" ) + Dt;
                PulsarAddParam( pulparams[h], "T0", &T0val, PULSARTYPE_REAL8_t );
              }

              params.orbit.tp.gpsSeconds = ( UINT4 )floor( PulsarGetREAL8ParamOrZero( pulparams[h], "T0" ) );
              params.orbit.tp.gpsNanoSeconds = ( UINT4 )floor( ( PulsarGetREAL8ParamOrZero( pulparams[h], "T0" ) - params.orbit.tp.gpsSeconds ) * 1e9 );
            } /* if is a BINARY */
            else {
              XLAL_INIT_MEM( params.orbit ); /*if not a binary pulsar*/
            }
            params.site = site;
            params.ephemerides = edat;
            params.startTimeGPS = epoch;
            params.duration = ndata; /*maybe want to take out conversion and keep this uvar->duration*/
            params.samplingRate = srate; /*same as above*/
            params.fHeterodyne = 0.;

            REAL4TimeSeries *TSeries = NULL;
            /*fprintf(stderr, "status%i\n", status.statusCode);*/
            LALGeneratePulsarSignal( &status, &TSeries, &params );
            if ( status.statusCode ) {
              if ( ( logfile = fopen( log_file, "a" ) ) == NULL ) {
                fprintf( stderr, "Error opening .log file! \n" );
                return 0;
              } else {
                fprintf( logfile, "LAL Routine failed at pulsar %s\n", pulin );
                XLAL_ERROR( XLAL_EFAILED );
                fclose( logfile );
              }
              fprintf( stderr, "LAL Routine failed\n" );
              XLAL_ERROR( XLAL_EFAILED );
            }
            /* add that timeseries to a common one */
            for ( i = 0; i < TSeries->data->length; i++ ) {
              total_inject->data->data[i] += TSeries->data->data[i];
            }
            /*    if ((total_inject = XLALFrameAddREAL8TimeSeriesProcData(frfile, TSeries)) == NULL)*/
            /*fprintf(stderr, "your pulsar is %s\n", pulin);*/
            XLALDestroyREAL4TimeSeries( TSeries );
            XLALDestroyREAL8Vector( params.pulsar.spindown );
            /*fclose( inject ); close .par file */
          } /* for strstr .par  */

          free( parnamelist[h] );

        } /*for k<num par files */

        free( parnamelist );
        /*add to series*/
        for ( i = 0; i < series->data->length; i++ ) {
          series->data->data[i] = gwfseries->data->data[i] + total_inject->data->data[i];
          /*fprintf(stderr, "%le\n", series->data->data[i]);*/
        }

        if ( XLALFrameAddREAL8TimeSeriesProcData( outFrame, series ) ) {
          LogPrintf( LOG_CRITICAL, "%s : XLALFrameAddINT4TimeSeries() failed with error = %d.\n", fn, xlalErrno );
          XLAL_ERROR( XLAL_EFAILED );
        }

        if ( XLALFrameAddREAL8TimeSeriesProcData( outFrame, total_inject ) ) {
          LogPrintf( LOG_CRITICAL, "%s : XLALFrameAddREAL8TimeSeriesProcData() failed with eerror = %d.\n", fn, xlalErrno );
          XLAL_ERROR( XLAL_EFAILED );
        }


        /* write frame structure to file (opens, writes, and closes file) - last argument is compression level */
        if ( XLALFrameWrite( outFrame, fffile ) ) {
          LogPrintf( LOG_CRITICAL, "%s : XLALFrameWrite() failed with error = %d.\n", fn, xlalErrno );
          XLAL_ERROR( XLAL_EFAILED );
        }

        /*free(fffile);*/

        XLALFrameFree( outFrame );

        XLALDestroyREAL8TimeSeries( gwfseries );
        XLALDestroyREAL8TimeSeries( total_inject );
        XLALDestroyREAL8TimeSeries( series );

        /* Free up the memory */

        /*  XLALDestroyREAL8Vector(injsig);*/

        /*  XLALDestroyREAL8TimeSeries( gwfseries );*/

        /*  XLALDestroyREAL8TimeSeries( series );*/

      } /*ends loop for creating CWINJ .gwf file */


      free( gwfnamelist[k] );

      /*Writing to .log file*/
      if ( ( logfile = fopen( log_file, "a" ) ) == NULL ) {
        fprintf( stderr, "Error opening .log file! \n" );
        return 0;
      } else {
        fprintf( logfile, "%s created successfully!\n", out_file );
        fclose( logfile );
      }
      /* </log file> */
      /*fprintf(stderr, "you created %s\n", out_file);*/

      XLALFrStreamClose( gwffile );

    } /*ends loop through all .gwf files in .gwf directory*/
  }

  free( gwfnamelist );
  XLALDestroyEphemerisData( edat );

  LALFree( site );

  fprintf( stderr, "done\n" );

  return 0;

} /* main() */


/**
 * Function to register and read all user input
 */
int
InitUserVars( UserInput_t *uvar,       /**< [out] UserInput structure to be filled */
              int argc,               /**< [in] number of argv element */
              char **argv             /**< [in] array of input arguments */
            )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !uvar ) {
    XLALPrintError( "%s: invalid NULL input 'uvar'\n", fn );
    XLAL_ERROR( XLAL_EINVAL );
  }
  if ( !argv ) {
    XLALPrintError( "%s: invalid NULL input 'argv'\n", fn );
    XLAL_ERROR( XLAL_EINVAL );
  }

  /* some defaults */
  uvar->srate = 16384;
  uvar->ephemEarth = XLALStringDuplicate( "earth00-40-DE405.dat.gz" );
  uvar->ephemSun = XLALStringDuplicate( "sun00-40-DE405.dat.gz" );

  /* Register User Variables*/
  /*    XLALRegisterUvarMember(out_chan,   STRING, 'o', OPTIONAL, "Output channel i.e. (IFO)_LDAS_C02_L2_CWINJ");*/
  /*    XLALRegisterUvarMember(in_chan,        STRING, 'i', OPTIONAL, "Input channel from .gwf file, i.e. (IFO):LDAS-STRAIN");*/
  XLALRegisterUvarMember( srate,            REAL8, 'r', OPTIONAL, "user defined sample rate, default = 16384" );
  /*  XLALRegisterUvarMember(duration,       REAL8, 'd', OPTIONAL, "duration of frame (sec)"); */
  /*  XLALRegisterUvarMember(start,            REAL8, 's', OPTIONAL, "epoch in GPS Seconds"); */
  XLALRegisterUvarMember( inputdir,       STRING, 'p', OPTIONAL, "directory for .par files" );
  XLALRegisterUvarMember( gwfdir,     STRING, 'g', OPTIONAL, "directory for .gwf files" );
  XLALRegisterUvarMember( outputdir,  STRING, 'c', OPTIONAL, "directory for CWINJ files" );
  XLALRegisterUvarMember( ephemEarth,   STRING, 0,  OPTIONAL,     "Earth ephemeris file to use" );
  XLALRegisterUvarMember( ephemSun,     STRING, 0,  OPTIONAL,     "Sun ephemeris file to use" );
  XLALRegisterUvarMember( IFO,       STRING, 'I', REQUIRED, "Detector: 'G1', 'L1', 'H1', 'H2', 'V1'..." );
  XLALRegisterUvarMember( logDir, STRING, 'L', OPTIONAL, "Directory to put .log file" );

  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv, lalPulsarVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( should_exit ) {
    exit( 1 );
  }

  return XLAL_SUCCESS;

} /* InitUserVars() */

/* compare to chars - taken from LALFrameIO.c */
static int charcmp( const void *c1, const void *c2 )
{
  char a = *( const char * )c1;
  char b = *( const char * )c2;
  return ( a > b ) - ( a < b );
}

/* function to set the Frame file name from the input data - taken from LALFrameIO.c */
int XLALFrameFileName( char *fname, size_t size, const char *chname, const LIGOTimeGPS *epoch, double duration )
{
  char site[LAL_NUM_DETECTORS + 1] = "";
  char *desc;
  const char *cs;
  char *s;
  int detflgs = 0;
  int t0;
  int dt;

  /* parse chname to get identified sites and detectors */
  /* strip out detectors from "XmYn...:"-style prefix */
  for ( cs = chname; *cs; cs += 2 ) {
    int d;
    /* when you get to a colon, you're done! */
    if ( *cs == ':' ) {
      break;
    }
    /* see if this is an unexpected format */
    if ( strlen( cs ) <= 2 || !isupper( cs[0] ) || !isdigit( cs[1] ) ) {
      /* parse error so reset detflgs and site */
      detflgs = 0;
      site[0] = 0;
      break;
    }
    /* try to find this detector */
    for ( d = 0; d < LAL_NUM_DETECTORS; ++d )
      if ( 0 == strncmp( cs, lalCachedDetectors[d].frDetector.prefix, 2 ) ) {
        /* found it: put it in sites and detflgs */
        detflgs |= 1 << 2 * d;
        strncat( site, cs, 1 );
      }
  }

  /* sort and uniqify sites */
  qsort( site, strlen( site ), 1, charcmp );
  cs = s = site;
  for ( cs = s = site; *s; ++cs )
    if ( *s != *cs ) {
      *++s = *cs;
    }

  /* description is a modified version of chname */
  /* replace invalid description char with '_' */
  desc = XLALStringDuplicate( chname );
  for ( s = desc; *s; ++s )
    if ( !isalnum( *s ) ) {
      *s = '_';
    }

  /* determine start time field and duration field */
  t0 = epoch->gpsSeconds;
  dt = ( int )ceil( XLALGPSGetREAL8( epoch ) + duration ) - t0;

  /* now format the file name */
  snprintf( fname, size, "%s-%s-%d-%d.gwf", *site ? site : "X", desc, t0,
            dt );

  LALFree( desc );
  return detflgs;
}
