/*
 *  create_pulsar_signal_frame.c: GW pulsar injection code
 *
 *  Copyright (C) 2012 Erin Macdonald, Matthew Pitkin
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

#include <stdio.h>
#include <unistd.h>
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
#include <getopt.h>

/*LAL Functions */
#include <lalapps.h>
#include <lal/Units.h>
#include <lal/FrameStream.h>
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

#define STRINGLENGTH 256              /* the length of general string */

/* create a detector at the GEO centre based on Virgo's arms */
//LALFrDetector FrGEOcentre = {
//  "C1", // the detector name
//  "C1", // detector prefix
//  0.0, // the vertex longitude (rads)
//  LAL_PI_2, // the vertex latitude (rads)
//  -LAL_BWGS84_SI, // the vertex elevation (m) - set to Earth centre
//  0.0, // the x-arm altitude (rads)
//  0.33916285222, // the x-arm azimuth (rads) - same as Virgo
//  0.0, // y-arm altitude (rads)
//  5.05155183261, // y-arm azimuth (rads) - same as Virgo
//  1500.0, // x-arm midpoint (m)
//  1500.0}; // y-arm midpoint (m)

#define USAGE \
"Usage: %s [options]\n\n"\
" --help, -h        display this message\n"\
" --detector, -i    the detector for which to create a frame (e.g. H1)\n"\
" --channel, -c     the channel name into which the signals will be added\n"\
" --epoch, -e       (int) the start epoch of the created frame\n"\
" --geocentre, g    set the detector to the geocentre\n"\
" --duration, -d    (int) the duration (in seconds) of the frame\n"\
" --pulsar-dir, -p  the directory containing pulsar (.par) files of signals\n\
                   to be added to the frame\n"\
" --output-dir, -o  the directory in which to output the frame\n"\
" --output-str, -s  a string to start the frame file name\n"\
" --ephem-dir, -m   the directory containing ephemeris files\n"\
" --ephem-type, -y  the ephemeris file type to use (e.g. DE405 [default])\n"\
" --dbg-lvl, -l     (int) the code debug level\n"\
"\n"

/*function to read ephemeris files*/
EphemerisData * InitEphemeris (const CHAR *ephemType, const CHAR *ephemDir );

LALStatus empty_LALStatus;

typedef struct tagInputParams{
  CHAR *ephemDir; /* ephemeris directory */
  CHAR *ephemType; /* ephemeric year */

  CHAR *pulsarDir; /* directory containing par files */

  CHAR *det; /* detector */

  CHAR *channel; /* channel name */

  CHAR *outDir; /* output directory */
  CHAR *outStr; /* string for output file name */

  UINT4 frDur; /* duration of a frame */
  UINT4 epoch; /* start epoch of a frame */

  UINT4 geocentre; /* a flag to set the detector to the geocentre */
} InputParams;

const InputParams empty_InputParams;

void ReadInput(InputParams *inputParams, int argc, char *argv[]);

/*--------------main function---------------*/
int main(int argc, char **argv){
  const CHAR *fn = __func__;

  InputParams inputs = empty_InputParams;

  REAL8 srate = 16384.0; /*sample rate defaulted to 16384 */

  /* read in command line input args */
  ReadInput( &inputs, argc, argv );

  LALStatus status = empty_LALStatus;

  EphemerisData *edat;
  if ( (edat = InitEphemeris ( inputs.ephemType, inputs.ephemDir)) == NULL ){
    XLALPrintError ( "%s: Failed to init ephemeris data\n", fn );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /*init detector info */
  LALDetector *site;
  if ( ( site = XLALGetSiteInfo ( inputs.det )) == NULL ){
    XLALPrintError("%s: Failed to get site-info for detector '%s'\n", fn,
                   inputs.det );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  if( inputs.geocentre ){ /* set site to the geocentre */
    site->location[0] = 0.0;
    site->location[1] = 0.0;
    site->location[2] = 0.0;
  }

  struct dirent **pulsars;
  INT4 n=scandir(inputs.pulsarDir, &pulsars, 0, alphasort);
  if ( n < 0){
    XLALPrintError("scandir failed\n");
    XLAL_ERROR(XLAL_EIO);
  }

  UINT4 numpulsars = (UINT4)n;
  UINT4 h=0;

  CHAR parname[256];
  BinaryPulsarParams pulparams[numpulsars];

  for(h=2; h<numpulsars; h++){
    if(strstr(pulsars[h]->d_name,".par") == NULL){
      free(pulsars[h]);
      continue;
    }
    else{
      sprintf(parname,"%s/%s", inputs.pulsarDir, pulsars[h]->d_name);
      fprintf(stderr, "%s\n", parname);
      FILE *inject;

      if (( inject = fopen ( parname, "r" )) == NULL ){
        fprintf(stderr,"Error opening file: %s\n", parname);
        XLAL_ERROR ( XLAL_EIO );
      }

      XLALReadTEMPOParFile( &pulparams[h], parname );
      fclose( inject );
    }
  }
  LIGOTimeGPS epoch;

  UINT4 ndata;

  epoch.gpsSeconds = inputs.epoch;
  epoch.gpsNanoSeconds = 0;

  ndata = inputs.frDur;

  REAL8TimeSeries *series=NULL;

  CHAR out_file[256];
  sprintf(out_file, "%s-%s-%d-%d.gwf", inputs.det, inputs.outStr,
          epoch.gpsSeconds, ndata );

  struct FrameH *outFrame = NULL;

  if ((outFrame = XLALFrameNew( &epoch, (REAL8)ndata, inputs.channel, 1, 0,
       0 )) == NULL) {
    LogPrintf(LOG_CRITICAL, "%s : XLALFrameNew() filed with error = %d.\n", fn, xlalErrno);
    XLAL_ERROR( XLAL_EFAILED);
  }

  if ((series = XLALCreateREAL8TimeSeries( inputs.channel, &epoch, 0.,
    1./srate,&lalSecondUnit, (int)(ndata*srate) )) == NULL) {
    XLAL_ERROR( XLAL_EFUNC );
  }

  UINT4 counter=0;
  for (counter = 0; counter < series->data->length; counter++)
    series->data->data[counter] = 0;

  /*** Read Pulsar Data ***/
  for (h=0; h < numpulsars; h++){
    if(strstr(pulsars[h]->d_name,".par")==NULL){
      free(pulsars[h]);
      continue;
    }
    else{
      PulsarSignalParams params = empty_PulsarSignalParams;

      /* set signal generation barycenter delay look-up table step size */
      params.dtDelayBy2 = 10.; /* generate table every 10 seconds */

      if (( params.pulsar.spindown = XLALCreateREAL8Vector(1)) == NULL ){
        XLALPrintError("Out of memory");
        XLAL_ERROR ( XLAL_EFUNC );
      }

      if (pulparams[h].posepoch == 0.0 )
        pulparams[h].posepoch = pulparams[h].pepoch;
      INT4 dtpos = epoch.gpsSeconds - (INT4)pulparams[h].posepoch;

      params.pulsar.position.latitude = pulparams[h].dec + dtpos * pulparams[h].pmdec;
      params.pulsar.position.longitude = pulparams[h].ra +dtpos * pulparams[h].pmra / cos(params.pulsar.position.latitude);
      params.pulsar.position.system = COORDINATESYSTEM_EQUATORIAL;

      params.pulsar.f0 = 2.*pulparams[h].f0;
      params.pulsar.spindown->data[0] = 2.*pulparams[h].f1;
      if (( XLALGPSSetREAL8(&(params.pulsar.refTime), pulparams[h].pepoch) ) == NULL )
        XLAL_ERROR ( XLAL_EFUNC );
      params.pulsar.psi = pulparams[h].psi;
      params.pulsar.phi0 = pulparams[h].phi0;
      params.pulsar.aPlus = 0.5 * pulparams[h].h0 * (1. + pulparams[h].cosiota * pulparams[h].cosiota );
      params.pulsar.aCross = pulparams[h].h0 * pulparams[h].cosiota;

      /*Add binary later if needed!*/

      params.site = site;
      params.ephemerides = edat;
      params.startTimeGPS = epoch;
      params.duration = ndata;
      params.samplingRate = srate;
      params.fHeterodyne = 0.;

      REAL4TimeSeries *TSeries = NULL;

      LALGeneratePulsarSignal( &status, &TSeries, &params );

      if (status.statusCode){
        fprintf(stderr, "LAL Routine failed!\n");
        XLAL_ERROR (XLAL_EFAILED);
      }
      UINT4 i;
      for (i=0; i < TSeries->data->length; i++)
        series->data->data[i] += TSeries->data->data[i];

      XLALDestroyREAL4TimeSeries(TSeries);
      XLALDestroyREAL8Vector(params.pulsar.spindown);
    }
  }

  if (XLALFrameAddREAL8TimeSeriesProcData(outFrame,series)){
      LogPrintf(LOG_CRITICAL, "%s : XLALFrameAddREAL8TimeSeries() failed with error = %d.\n",fn,xlalErrno);
      XLAL_ERROR(XLAL_EFAILED);
  }

  CHAR OUTFILE[256];
  sprintf(OUTFILE, "%s/%s", inputs.outDir, out_file);

  if (  XLALFrameWrite(outFrame, OUTFILE,1)){
    LogPrintf(LOG_CRITICAL, "%s : XLALFrameWrite() failed with error = %d.\n", fn, xlalErrno);
    XLAL_ERROR( XLAL_EFAILED );
  }

  FrameFree(outFrame);
  XLALDestroyREAL8TimeSeries( series );

  return 0;
}

EphemerisData *InitEphemeris (const CHAR *ephemType, const CHAR *ephemDir){
  const CHAR *fn = __func__;
#define FNAME_LENGTH 1024
  CHAR EphemEarth[FNAME_LENGTH];  /* filename of earth-ephemeris data */
  CHAR EphemSun[FNAME_LENGTH];    /* filename of sun-ephemeris data */

  /* check input consistency */
  if ( !ephemType ) {
    XLALPrintError ("%s: invalid NULL input for 'ephemType'\n", fn );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  snprintf(EphemEarth, FNAME_LENGTH, "%s/earth00-19-%s.dat.gz", ephemDir, ephemType);
  snprintf(EphemSun, FNAME_LENGTH, "%s/sun00-19-%s.dat.gz", ephemDir, ephemType);

  EphemEarth[FNAME_LENGTH-1]=0;
  EphemSun[FNAME_LENGTH-1]=0;

  EphemerisData *edat;
  if ( (edat = XLALInitBarycenter ( EphemEarth, EphemSun)) == NULL ) {
    XLALPrintError ("%s: XLALInitBarycenter() failed.\n", fn );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  /* return ephemeris */
  return edat;
} /* InitEphemeris() */

void ReadInput(InputParams *inputParams, int argc, char *argv[]){
  struct option long_options[] =
  {
    { "help",                     no_argument,        0, 'h' },
    { "detector",                 required_argument,  0, 'i' },
    { "channel",                  required_argument,  0, 'c' },
    { "epoch",                    required_argument,  0, 'e' },
    { "geocentre",                no_argument,        0, 'g' },
    { "duration",                 required_argument,  0, 'd' },
    { "pulsar-dir",               required_argument,  0, 'p' },
    { "output-dir",               required_argument,  0, 'o' },
    { "output-str",               required_argument,  0, 's' },
    { "ephem-dir",                required_argument,  0, 'm' },
    { "ephem-type",               required_argument,  0, 'y' },
    { "dbg-lvl",                  required_argument,  0, 'l' },
    { 0, 0, 0, 0 }
  };

  const CHAR *fn = __func__;
  CHAR args[] = "hi:c:e:gd:p:o:s:m:y:l:";
  CHAR *program = argv[0];

  /* default the debug level to 1 */

  /* default to no set the detector to the geocentre */
  inputParams->geocentre = 0;

  /* default ephemeris to use DE405 */
  inputParams->ephemType = XLALStringDuplicate( "DE405" );

  /* get input arguments */
  while(1){
    INT4 option_index = 0;
    INT4 c;

    c = getopt_long( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch(c){
      case 0: /* if option set a flag, nothing else to do */
        if ( long_options[option_index].flag )
          break;
        else
          fprintf(stderr, "Error parsing option %s with argument %s\n",
            long_options[option_index].name, optarg );
      case 'h': /* help message */
        fprintf(stderr, USAGE, program);
        exit(0);
      case 'l': /* debug level */
        break;
      case 'i': /* interferometer/detector */
        inputParams->det = XLALStringDuplicate( optarg );
        break;
      case 'c': /* channel name */
        inputParams->channel = XLALStringDuplicate( optarg );
        break;
      case 'e': /* frame epoch */
        inputParams->epoch = atoi(optarg);
        break;
      case 'g': /* geocentre flag */
        inputParams->geocentre = 1;
        break;
      case 'd': /* frame duration */
        inputParams->frDur = atoi(optarg);
        break;
      case 'p': /* pulsar par file directory */
        inputParams->pulsarDir = XLALStringDuplicate( optarg );
        break;
      case 'o': /* output directory */
        inputParams->outDir = XLALStringDuplicate( optarg );
        break;
      case 's': /* output name string */
        inputParams->outStr = XLALStringDuplicate( optarg );
        break;
      case 'm': /* ephemeris file directory */
        inputParams->ephemDir = XLALStringDuplicate( optarg );
        break;
      case 'y': /* ephemeris file year */
        inputParams->ephemType = XLALStringDuplicate( optarg );
        break;
      case '?':
        fprintf(stderr, "unknown error while parsing options\n" );
      default:
        fprintf(stderr, "unknown error while parsing options\n" );
    }
  }

  if( inputParams->epoch == 0 || inputParams->frDur == 0 ){
    XLALPrintError ("%s: Frame epoch or duration are 0!\n", fn );
    XLAL_ERROR_VOID ( XLAL_EFUNC );
  }
}
