/*
 * Copyright (C) 2010 Erin Macdonald
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

/**
 * \file
 * \ingroup pulsarApps
 * \author Erin Macdonald
 * \brief
 * Code to create frames and add them to existing data files
 */

/* Code to create frames and add them to existing data files
Input format as $ ./sw_inj_frames framefile duration epoch
example$ ./lalapps_sw_inj_frames -p /Users/erinmacdonald/lsc/analyses/test_par_files -g /Users/erinmacdonald/lsc/analyses/frames -c /Users/erinmacdonald/lsc/analyses/CWINJframes -I H1 -e /Users/erinmacdonald/lsc/src/lscsoft/lalsuite/lalpulsar/test -y 09-11
*/

/* 2/3/11 v. 2 (which is not entirely accurate, but I'm starting from here): added a log file, version number, and separate channel for i(t) -- E. Macdonald */
/* 3/3/11 v. 3: placed log comments throughout */
/* 7/3/11 v. 4: output file to dump ascii of failed frames */

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

/*LAL Functions */
#include <lalapps.h>
#include <lal/Units.h>
#include <lal/FrameStream.h>
#include <lal/LALFrameIO.h>
#include <lal/FrameCache.h>
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
/*#include <lal/PulsarDataTypes.h>*/

#define STRINGLENGTH 256              /* the length of general string */

/* ---------- local type definitions ---------- */
typedef struct{
  BOOLEAN help; /* Trigger output of help string*/
/*    CHAR *out_chan;  Channel output ie. H1_LDAS_C02_L2_CWINJ*/
/*    CHAR *in_chan;  Channel input from .gwf ie. H1:LDAS-STRAIN*/
  REAL8 srate; /* sample rate -- default 16384*/
  /*  REAL8 duration; duration (sec) Now read from file name */
  /*  REAL8 start; epoch (GPSSeconds) Now read from file name */
  CHAR *inputdir; /* directory for .par files*/
  CHAR *gwfdir; /*directory for .gwf files*/
  CHAR *outputdir; /*directory for CWINJ files */
  CHAR *ephemDir; /*directory for ephemeris files*/
  CHAR *ephemYear; /*ephemeris years*/
  CHAR *IFO; /*detector */
  CHAR *logDir; /*directory for the log files */
} UserInput_t;
  
UserInput_t uvar_struct;    

/* ---------- local function prototypes ---------- */
EphemerisData * XLALInitEphemeris (const CHAR *ephemYear, const CHAR *ephemDir ); /*function to read ephemeris files*/

int InitUserVars ( UserInput_t *uvar, int argc, char **argv ); /*Initiates user variables*/

/* empty initialiser */
LALStatus empty_LALStatus;


/* ---------- Function definitions ---------- */

int main(int argc, char **argv)
{
  const char *fn = __func__;
  
  UserInput_t *uvar = &uvar_struct; /*structure for user-defined input variables*/
 
  /* read all user input */
  if ( InitUserVars ( uvar, argc, argv ) != XLAL_SUCCESS ) {
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }
 
  if ( uvar->help )
    return 0;  

  char version[256];
  sprintf(version, "v4"); /*manually change */

  /*Get current time for log */
  time_t result;
  result = time(NULL);
  struct tm* btm = localtime(&result);

  /*Start the log file */  
  FILE *logfile;
  char log_file[256];

  sprintf(log_file, "%s/%s/test.log", uvar->outputdir, uvar->logDir); /*How to name log files separately?*/
  fprintf(stderr, "Your log file is here: %s\n", log_file);

  /*Writing to .log file*/
  if (( logfile = fopen( log_file, "a" )) == NULL ) {
    fprintf (stderr, "Error opening .log file! \n" );
    return 0;
  }
  else {
    fprintf (logfile, "\nsw_inj_frames.c version number: %s\nCurrent time: %s \n", version, asctime(btm));
    fprintf (logfile, "User inputs:\n Sample rate: %f\n Pulsar file directory: %s\n Raw frame file directory: %s\n Output directory: %s\n Years for ephemeris: %s\n Directory for ephemeris: %s\n Name of directory for log (w/in output dir): %s\n IFO: %s\n", uvar->srate, uvar->inputdir, uvar->gwfdir, uvar->outputdir, uvar->ephemYear, uvar->ephemDir, uvar->logDir, uvar->IFO);
    fclose(logfile);
  }
  /* </log file> */

  struct dirent **parnamelist;
  struct dirent **gwfnamelist;
  char pulin[256];
  
  lalDebugLevel = 1;	/* debug level for this code */
  LALStatus status = empty_LALStatus;
  
  PulsarSignalParams params = empty_PulsarSignalParams; /*pulsar parameter structure*/
  BinaryPulsarParams pulparams; /* read from the .par file */
  
  /*  sprintf(lalframefile, "%s", argv[1]);  User defined frame file -- need to simplify */
  /*  sprintf(gwfframefile, "%s", argv[2]);  Frame file to be read in*/
  
  /*  srate = atoi(argv[3]);  User defined sample rate (16384)*/
  
  /*  ndata = atoi(argv[4]);  length of data set */
  
  /*  epoch.gpsSeconds = atoi(argv[5]);  User defined gps epoch */
  /*  epoch.gpsNanoSeconds = 0;*/
  
  /*init ephemeris-data */
  EphemerisData *edat;
  if ( (edat = XLALInitEphemeris ( uvar->ephemYear, uvar->ephemDir )) == NULL ) {
    XLALPrintError ( "%s: Failed to init ephemeris data for year-span '%s'\n", fn, uvar->ephemYear );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }
  
  /*init detector info */
  LALDetector *site;
  if ( ( site = XLALGetSiteInfo (uvar -> IFO )) == NULL ){
    XLALPrintError("%s: Failed to get site-info for detector '%s'\n", fn, uvar->IFO );
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }
  
  
  if ( (params.pulsar.spindown = XLALCreateREAL8Vector(1))== NULL ) {
    XLALPrintError("Out of memory");
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }
  
  /*Error Checks*/
  /*    if (uvar->in_chan == 0){*/
  /*        XLALPrintError ("\nNeed an input channel!\n");*/
  /*        return 1;*/
  /*    }*/
  
  /*  UINT4 ndata; -- done in loop below from file name */
  /* ndata = uvar->duration; */

  REAL8 srate; /*sample rate defaulted to 16384 */
  srate = uvar->srate;
  
  /*creates the .gwf channels */
  CHAR in_chan[256];
  sprintf( in_chan, "%s:LDAS-STRAIN", uvar->IFO );
  /*fprintf( stderr, "\nin_chan = %s\n", in_chan);*/

  CHAR out_chan[256];
  sprintf( out_chan, "%s_LDAS_C02_L2_CWINJ_TOT", uvar->IFO );
  /*fprintf( stderr, "\nout_chan = %s\n", out_chan);*/

  CHAR inj_chan[256];
  sprintf( inj_chan, "%s_LDAS_C02_L2_CWINJ_SIG", uvar->IFO );

  /*Writing to .log file*/
  if (( logfile = fopen( log_file, "a" )) == NULL ) {
    fprintf (stderr, "Error opening .log file! \n" );
    return 0;
  }
  else {
    fprintf (logfile, "Channels:\n Reading from raw: %s\n Injection+noise: %s\n Pure injection: %s\n", in_chan, out_chan, inj_chan);
    fclose(logfile);
  }
  /* </log file> */


  /*  extract .gwf file name from inputs*/
  /*  pos = strchr(uvar->out_chan, '_');*/
  /*  ipos = pos-(uvar->out_chan);*/
  /*  strncpy(detname, uvar->out_chan, ipos);*/
  /*  fprintf(stderr, "%s\n", detname);*/
  
  /*strcpy(channame, uvar->out_chan[ipos+1]);
    
    strncpy(injpos, uvar->out_chan, 14);*/
  
  /* Major loop to run through gwf files */
  int m;
  UINT4 i;
  
  /*  REAL8TimeSeries *gwfseries=NULL;*/
  /*REAL8TimeSeries *series=NULL;*/
  
  char out_file[256]; /*need a different way to do this */
  FILE *inject; /*for .par file from uvar->inputdir*/

  CHAR gwf_dir[256]; /* for use with XLALFrOpen */
  sprintf( gwf_dir, "%s/.", uvar->gwfdir );

  LIGOTimeGPS epoch;
  char strepoch[10];

  UINT4 ndata;
  size_t filength;

  /* Put the pulsar files in the log (so as not to loop every time) */
  /*Read in all pulsar files from inputdir*/
  int n = scandir(uvar->inputdir, &parnamelist, 0, alphasort);
  if ( n < 0 ) {
    XLALPrintError ("%s: scandir() failed for directory '%s'\n", fn, uvar->inputdir );
    XLAL_ERROR ( fn, XLAL_EIO );
  }
  UINT4 numParFiles = (UINT4)n;
  UINT4 h=0;
  /*h starts at 2 to skip over . and .. found in scandir, they are the first two entries */
  for (h=2; h < numParFiles; h++) {
    if(strstr(parnamelist[h]->d_name, ".par") == NULL){
      fprintf(stderr, "Not using file %s\n",parnamelist[h]->d_name);
      continue;
    }
    else{
      /*Writing to .log file*/
      if (( logfile = fopen( log_file, "a" )) == NULL ) {
	fprintf (stderr, "Error opening .log file! \n" );
	return 0;
      }
      else {
	char parname[256];
	sprintf(parname, "%s/%s", uvar->inputdir, parnamelist[h]->d_name);
	char t[ 100 ];
	struct stat b;
	if (!stat(parname, &b)) {
	  strftime(t, 100, "%d/%m/%Y %H:%M:%S", localtime( &b.st_mtime));
	  fprintf (logfile, "%s/%s Last Modified: %s\n", uvar->inputdir, parnamelist[h]->d_name, t);
	  fclose(logfile);
	}
      }
      /* </log file> */
    } /*checking if file ends in .par*/
  }/*looping through files in uvar->inputdir */

  if ( (m = scandir(uvar->gwfdir, &gwfnamelist, 0, alphasort)) == -1) {
    XLALPrintError ("%s: scandir('%s',...) failed.\n", fn, uvar->gwfdir );
    XLAL_ERROR ( fn, XLAL_EIO );
  }
  
  UINT4 k=0;
  for (k=2; k < (UINT4)m; k++){
    if (strstr(gwfnamelist[k]->d_name, ".gwf") == NULL){
      continue; /*make sure it's a .gwf file */
    }
    else{
      FrStream *gwffile=NULL; 
      
      if (( gwffile = XLALFrOpen( uvar->gwfdir, gwfnamelist[k]->d_name)) == NULL ) {
	/*XLAL_ERROR ( fn, XLAL_EFUNC ); -- don't want to abort, but save elsewhere test that it's an acceptable file */
	
	/*Writing to failed file*/
	FILE *frames;
	char framefilename[256];
	sprintf(framefilename, "%s/%s/failed_frames.txt", uvar->outputdir, uvar->logDir);
        if (( frames = fopen( framefilename, "a" )) == NULL ) {
          fprintf (stderr, "Error opening file %s! \n", framefilename );
          return 0;
        }
        else {
          fprintf (frames, "%s\n", gwfnamelist[k]->d_name);
	  fclose(logfile);
        }
        /* </failed file> */
	continue; /*don't exit program if .gwf file fails, continue through*/
      }
      else {
	/*Writing to .log file*/
	if (( logfile = fopen( log_file, "a" )) == NULL ) {
	  fprintf (stderr, "Error opening .log file! \n" );
	  return 0;
	}
	else {
	  fprintf (logfile, "Read file %s/%s\n", uvar->gwfdir, gwfnamelist[k]->d_name);
	  fclose(logfile);
	}
	/* </log file> */

	REAL8TimeSeries *gwfseries=NULL;
	REAL8TimeSeries *series=NULL;

	/** extract epoch and duration from gwf file name **/
	/*LIGOTimeGPS epoch;*/
	/*char strepoch[10];*/
	strncpy(strepoch, strchr(gwfnamelist[k]->d_name, '9'), 9 ); /*All epochs in S6 begin with 9... potential problem in future */
	strepoch[sizeof(strepoch)-1] = '\0'; /*Null terminate the string*/
	epoch.gpsSeconds = atoi(strepoch);  /* convert to integer from string */
	epoch.gpsNanoSeconds = 0; /* no nanosecond precision */
	/*	fprintf(stderr, "epoch = %i\n", epoch.gpsSeconds);*/
	
	/*UINT4 ndata;*/
	/*size_t filength;*/
	filength = strlen(gwfnamelist[k]->d_name); 
	char strdur[4];
	strncpy(strdur, (strrchr(gwfnamelist[k]->d_name, '-')+1), 3); /* duration is last number in frame file */
	strdur[sizeof(strdur)-1] = '\0';
	/*	fprintf(stderr, "Duration = %s\n", strdur);*/
	/* assigns duration from .gwf frame */
	ndata = atoi(strdur);

	/*Writing to .log file*/
	if (( logfile = fopen( log_file, "a" )) == NULL ) {
	  fprintf (stderr, "Error opening .log file! \n" );
	  return 0;
	}
	else {
	  fprintf (logfile, "Using epoch: %i and duration: %i\n", epoch.gpsSeconds, ndata);
	  fclose(logfile);
	}
	/* </log file> */
	
	/* acquire time series from frame file */
	if ( (gwfseries = XLALCreateREAL8TimeSeries( in_chan, &epoch, 0., 1./srate, &lalSecondUnit, (int)(ndata*srate) )) == NULL ) {
	  XLAL_ERROR ( fn, XLAL_EFUNC );
	}
	
	if ( XLALFrGetREAL8TimeSeries( gwfseries, gwffile ) != XLAL_SUCCESS ) {
	  XLAL_ERROR ( fn, XLAL_EFUNC );
	}
	
	/* define output .gwf file */
	sprintf(out_file, "%c-%s-%d-%d.gwf", uvar->IFO[0], out_chan, epoch.gpsSeconds, ndata);
	
	/*Writing to .log file*/
	if (( logfile = fopen( log_file, "a" )) == NULL ) {
	  fprintf (stderr, "Error opening .log file! \n" );
	  return 0;
	}
	else {
	  fprintf (logfile, "Writing to %s/%s\n", uvar->outputdir, out_file);
	  fclose(logfile);
	}
	/* </log file> */
	
	/* read in and test generated frame with XLAL function*/
	
	struct FrameH *outFrame   = NULL;        /* frame data structure */
	
	if ((outFrame = XLALFrameNew(&epoch,(REAL8)ndata,"CW_INJ",1,0,0)) == NULL) {
	  LogPrintf(LOG_CRITICAL, "%s : XLALFrameNew() failed with error = %d.\n",fn,xlalErrno);
	  XLAL_ERROR(fn,XLAL_EFAILED);
	}
	
	series=NULL;
	if ((series = XLALCreateREAL8TimeSeries( out_chan, &epoch, 0., 1./srate, &lalSecondUnit, (int)(ndata*srate) )) == NULL){
	  XLAL_ERROR ( fn, XLAL_EFUNC );
	}
	CHAR fffile[256];
	snprintf(fffile,STRINGLENGTH,"%s/%s",uvar->outputdir,out_file);
	if (fopen(fffile , "w+") == NULL){
	  fprintf(stderr, "fail");
	}
	
	/*if (( frfile = XLALFrOpen( uvar->outputdir, out_file )) == NULL){
	  XLAL_ERROR ( fn, XLAL_EFUNC );
	  }*/
	/*create series to be the sum, general series to add to the noise */
	REAL8TimeSeries *total_inject=NULL;
	if (( total_inject = XLALCreateREAL8TimeSeries(inj_chan, &epoch, 0., 1./srate, &lalSecondUnit, (UINT4)(ndata*srate)) ) == NULL) {
	  XLAL_ERROR ( fn, XLAL_EFUNC );
	}
	
	/* need to set all values in total_inject to zero so not pre-allocated */
	UINT4 counter=0;
	for (counter = 0; counter < total_inject->data->length; counter++)
	  total_inject->data->data[counter] = 0;
	
	/*Read in all pulsar files from inputdir*/
	n = scandir(uvar->inputdir, &parnamelist, 0, alphasort);
	if ( n < 0 ) {
	  XLALPrintError ("%s: scandir() failed for directory '%s'\n", fn, uvar->inputdir );
	  XLAL_ERROR ( fn, XLAL_EIO );
	}
	numParFiles=(UINT4)n;
	h=0;
	/*h starts at 2 to skip over . and .. found in scandir, they are the first two entries */
	for (h=2; h < numParFiles; h++) {
	  if(strstr(parnamelist[h]->d_name, ".par") == NULL){
	    fprintf(stderr, "Not using file %s\n",parnamelist[h]->d_name);
	    continue;
	  }
	  else{
	    sprintf(pulin, "%s/%s", uvar->inputdir, parnamelist[h]->d_name);
	    /*fprintf(stderr, "This is your pulsar file:%s\n", pulin);*/
	    if (( inject = fopen ( pulin, "r" )) == NULL ){
	      fprintf(stderr, "Error opening file: %s\n", pulin);
	      XLAL_ERROR ( fn, XLAL_EIO );
	    }
	    
	    /*read in parameters from .par file */
	    XLALReadTEMPOParFile( &pulparams, pulin);
	    params.pulsar.position.longitude = pulparams.ra;
	    params.pulsar.position.latitude = pulparams.dec;
	    params.pulsar.position.system = COORDINATESYSTEM_EQUATORIAL;
	    params.pulsar.f0 = 2.*pulparams.f0;
	    params.pulsar.spindown->data[0] = 2.*pulparams.f1; /*spindown is REAL8Vector ?? */
	    if (( XLALGPSSetREAL8(&(params.pulsar.refTime),pulparams.pepoch) ) == NULL ){
	      XLAL_ERROR ( fn, XLAL_EFUNC );
	    }
	    params.pulsar.psi = pulparams.psi;
	    params.pulsar.phi0 = pulparams.phi0;
	    /*Conversion from h0 and cosi to plus and cross */
	    params.pulsar.aPlus = 0.5 * pulparams.h0 * (1. + pulparams.cosiota * pulparams.cosiota ); 
	    params.pulsar.aCross = pulparams.h0 * pulparams.cosiota;
	    
	    /*if binary */
	    if (pulparams.model != NULL) {
	      params.orbit->asini = pulparams.x;
	      params.orbit->period = pulparams.Pb*86400;
	      
	      /*Taking into account ELL1 model option */
	      if (strstr(pulparams.model,"ELL1") != NULL) {
		REAL8 w,e,eps1,eps2;
		eps1 = pulparams.eps1;
		eps2 = pulparams.eps2;
		w = atan2(eps1,eps2);
		e = sqrt(eps1*eps1+eps2*eps2);
		params.orbit->argp = w;
		params.orbit->ecc = e;
	      }
	      else {
		params.orbit->argp = pulparams.w0;
		params.orbit->ecc = pulparams.e;
	      }
	      if (strstr(pulparams.model,"ELL1") != NULL) {
		REAL8 fe, uasc,Dt;
		fe = sqrt((1.0-params.orbit->ecc)/(1.0+params.orbit->ecc));
		uasc = 2.0*atan(fe*tan(params.orbit->argp/2.0));
		Dt = (params.orbit->period/LAL_TWOPI)*(uasc-params.orbit->ecc*sin(uasc));
		pulparams.T0 = pulparams.Tasc + Dt;
	      }
	      
	      params.orbit->tp.gpsSeconds = (UINT4)floor(pulparams.T0);
	      params.orbit->tp.gpsNanoSeconds = (UINT4)floor((pulparams.T0 - params.orbit->tp.gpsSeconds)*1e9);
	    } /* if pulparams.model is BINARY */
	    else
	      params.orbit = NULL; /*if pulparams.model = NULL -- isolated pulsar*/
            
	    params.site = site;
	    params.ephemerides = edat;
	    params.startTimeGPS = epoch;
	    params.duration = ndata; /*maybe want to take out conversion and keep this uvar->duration*/
	    params.samplingRate = srate; /*same as above*/
	    params.fHeterodyne = 0.; /*not sure what this is, will check */
	    
	    REAL4TimeSeries *TSeries = NULL;
	    LALGeneratePulsarSignal (&status, &TSeries, &params);
	    if ( status.statusCode ) {
	      fprintf( stderr, "LAL Routine failed\n" );
	      XLAL_ERROR ( fn, XLAL_EFAILED );
	    }
	    /* add that timeseries to a common one */
	    for (i=0; i < TSeries->data->length; i++){
	      total_inject->data->data[i] += TSeries->data->data[i];
	    }
	    /*	  if ((total_inject = XLALFrameAddREAL8TimeSeriesProcData(frfile, TSeries)) == NULL)*/
	    XLALDestroyREAL4TimeSeries (TSeries);
	  } /* for strstr .par  */
	} /*for k<num par files */
	
	/*add to series*/
	for (i=0; i < series->data->length; i++){
	  series->data->data[i] = gwfseries->data->data[i] + total_inject->data->data[i];
	  /*fprintf(stderr, "%le\n", series->data->data[i]);*/
	}
	
	if (XLALFrameAddREAL8TimeSeriesProcData(outFrame,series)) {
	  LogPrintf(LOG_CRITICAL, "%s : XLALFrameAddINT4TimeSeries() failed with error = %d.\n",fn,xlalErrno);
	  XLAL_ERROR(fn,XLAL_EFAILED);
	}
	
	if ( XLALFrameAddREAL8TimeSeriesProcData(outFrame, total_inject)) {
	  LogPrintf(LOG_CRITICAL, "%s : XLALFrameAddREAL8TimeSeriesProcData() failed with eerror = %d.\n",fn, xlalErrno );
	  XLAL_ERROR(fn, XLAL_EFAILED);
	}
	
	/*fprintf(stderr,"%e\n", series->data->data[250]);*/
	
	/* write frame structure to file (opens, writes, and closes file) - last argument is compression level */
	if (XLALFrameWrite(outFrame,fffile,1)) {
	  LogPrintf(LOG_CRITICAL, "%s : XLALFrameWrite() failed with error = %d.\n",fn,xlalErrno);
	  XLAL_ERROR(fn,XLAL_EFAILED);
	}
	
	if (( XLALFrWriteREAL8TimeSeries( series, 1 ) ) != XLAL_SUCCESS ){
	  XLAL_ERROR( fn, XLAL_EFUNC );
	}
	
	/****Test for Matlab****/
	FILE *outtest;
        if (( outtest = fopen( "/Users/erinmacdonald/lsc/analyses/sw_injections/test.txt", "w" )) == NULL)
          fprintf(stderr, "Error opening file\n");
        else{
          for (i = 0; i < series->data->length; i++){
            fprintf(outtest,"%e\n",series->data->data[i] );
          }
          fclose(outtest);
        }
	FILE *injtest;
	if (( injtest = fopen( "/Users/erinmacdonald/lsc/analyses/sw_injections/injtest.txt","w" )) == NULL)
	  fprintf(stderr, "Error opening injection file\n");
	else{
	  for (i = 0; i < total_inject->data->length; i++){
	    fprintf(injtest,"%e\n",total_inject->data->data[i] );
	  }
	  fclose(injtest);
	}
        /****End test****/
	
	XLALDestroyREAL8TimeSeries( gwfseries);
	XLALDestroyREAL8TimeSeries( total_inject );
	XLALDestroyREAL8TimeSeries( series );
	
	/* Free up the memory */
	
	/*  XLALDestroyREAL8Vector(injsig);*/
	
	/*  XLALDestroyREAL8TimeSeries( gwfseries );*/
	
	/*  XLALDestroyREAL8TimeSeries( series );*/

      } /*ends loop for creating CWINJ .gwf file */
      /*Writing to .log file*/
      if (( logfile = fopen( log_file, "a" )) == NULL ) {
	fprintf (stderr, "Error opening .log file! \n" );
	return 0;
      }
      else {
	fprintf (logfile, "%s created successfully!\n", out_file);
	fclose(logfile);
      }
      /* </log file> */

    } /*ends loop through all .gwf files in .gwf directory*/
  }
  return 1;
    
} /* main() */


/** Function to register and read all user input
 */
int
InitUserVars ( UserInput_t *uvar,	/**< [out] UserInput structure to be filled */
	       int argc,		/**< [in] number of argv element */
	       char **argv		/**< [in] array of input arguments */
	       )
{
  const char *fn = __func__;

  /* check input consistency */
  if ( !uvar ) {
    XLALPrintError ("%s: invalid NULL input 'uvar'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }
  if ( !argv ) {
    XLALPrintError ("%s: invalid NULL input 'argv'\n", fn );
    XLAL_ERROR ( fn, XLAL_EINVAL );
  }

  /* some defaults */
#define EPHEM_YEARS "09-11"  
  uvar->ephemYear = XLALCalloc(1, strlen(EPHEM_YEARS)+1);
  strcpy (uvar->ephemYear, EPHEM_YEARS);
  
  uvar->srate=16384;

  /* Register User Variables*/
  XLALregBOOLUserStruct( help,            'h', UVAR_HELP, "Print this message");
  /*    XLALregSTRINGUserStruct(out_chan,   'o', UVAR_OPTIONAL, "Output channel i.e. (IFO)_LDAS_C02_L2_CWINJ");*/
  /*    XLALregSTRINGUserStruct(in_chan,        'i', UVAR_OPTIONAL, "Input channel from .gwf file, i.e. (IFO):LDAS-STRAIN");*/
  XLALregREALUserStruct(srate,            'r', UVAR_OPTIONAL, "user defined sample rate, default = 16384");
  /*  XLALregREALUserStruct(duration,       'd', UVAR_OPTIONAL, "duration of frame (sec)"); */
  /*  XLALregREALUserStruct(start,            's', UVAR_OPTIONAL, "epoch in GPS Seconds"); */
  XLALregSTRINGUserStruct(inputdir,       'p', UVAR_OPTIONAL, "directory for .par files");
  XLALregSTRINGUserStruct(gwfdir,     'g', UVAR_OPTIONAL,"directory for .gwf files");
  XLALregSTRINGUserStruct(outputdir,  'c', UVAR_OPTIONAL, "directory for CWINJ files");
  XLALregSTRINGUserStruct(ephemYear,      'y', UVAR_OPTIONAL,"Year(or range of years) for ephemeris files to be used");
  XLALregSTRINGUserStruct(ephemDir,   'e', UVAR_OPTIONAL,"Directory to find ephemeris files");
  XLALregSTRINGUserStruct( IFO,       'I', UVAR_REQUIRED, "Detector: 'G1', 'L1', 'H1', 'H2', 'V1'...");
  XLALregSTRINGUserStruct( logDir, 'L', UVAR_OPTIONAL, "Directory to put .log file");
  
  if (XLALUserVarReadAllInput (argc, argv ) != XLAL_SUCCESS) {
    XLALPrintError ("%s: XLALUserVarReadAllInput() failed with errno=%d\n", fn, xlalErrno);
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

} /* InitUserVars() */

EphemerisData *
XLALInitEphemeris (const CHAR *ephemYear, const CHAR *ephemDir)
{
  const char *fn = __func__;
#define FNAME_LENGTH 1024
  CHAR EphemEarth[FNAME_LENGTH];	/* filename of earth-ephemeris data */
  CHAR EphemSun[FNAME_LENGTH];	/* filename of sun-ephemeris data */
  
  /* check input consistency */
  if ( !ephemYear ) {
    XLALPrintError ("%s: invalid NULL input for 'ephemYear'\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EINVAL );
  }
  
  snprintf(EphemEarth, FNAME_LENGTH, "%s/earth%s.dat", ephemDir, ephemYear);
  snprintf(EphemSun, FNAME_LENGTH, "%s/sun%s.dat", ephemDir, ephemYear);  

  EphemEarth[FNAME_LENGTH-1]=0;
  EphemSun[FNAME_LENGTH-1]=0;
  
  EphemerisData *edat;
  if ( (edat = XLALInitBarycenter ( EphemEarth, EphemSun)) == NULL ) {
    XLALPrintError ("%s: XLALInitBarycenter() failed.\n", fn );
    XLAL_ERROR_NULL ( fn, XLAL_EFUNC );
  }
  
  /* return ephemeris */
  return edat;
  
} /* XLALInitEphemeris() */
