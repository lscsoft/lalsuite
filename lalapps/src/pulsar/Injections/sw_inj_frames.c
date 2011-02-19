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
example$ ./lalapps_sw_inj_frames -d 128 -s 945625472 -p /Users/erinmacdonald/lsc/analyses/test_par_files -g /Users/erinmacdonald/lsc/analyses/frames -c /Users/erinmacdonald/lsc/analyses/CWINJframes -I H1 -y 09-11
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
/*#include <lal/PulsarDataTypes.h>*/

/* ---------- local type definitions ---------- */
typedef struct{
  BOOLEAN help; /* Trigger output of help string*/
/*    CHAR *out_chan;  Channel output ie. H1_LDAS_C02_L2_CWINJ*/
/*    CHAR *in_chan;  Channel input from .gwf ie. H1:LDAS-STRAIN*/
  REAL8 srate; /* sample rate */
  REAL8 duration; /*duration (sec)*/
  REAL8 start; /*epoch (GPSSeconds)*/
  CHAR *inputdir; /* directory for .par files*/
  CHAR *gwfdir; /*directory for .gwf files*/
  CHAR *outputdir; /*directory for CWINJ files */
  CHAR *ephemDir; /*directory for ephemeris files*/
  CHAR *ephemYear; /*ephemeris years*/
  CHAR *IFO; /*detector */
} UserInput_t;
  
UserInput_t uvar_struct;    

/* ---------- local function prototypes ---------- */
EphemerisData * XLALInitEphemeris (const CHAR *ephemYear );
int InitUserVars ( UserInput_t *uvar, int argc, char **argv );


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
  
  /*char lalframefile[256];
    char gwfframefile[256];
    
    double srate;
    double ndata; *//* data taken every 1 sec */
  
  
  /*  REAL8Vector *Tstamp=NULL;*/
  /*  FILE *outtest;*/
  
  FrStream *frfile=NULL;
  FrStream *gwffile=NULL;
  
  struct dirent **namelist;
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
  
  /*   define the input directory for the mfd_files*/
  /*  sprintf(inputdir, "/home/emacdonald/par_files/generated/%s/mfd_files", argv[6]);*/
  /*  fprintf(stderr, "%s\n", inputdir);*/
  
  LIGOTimeGPS epoch;
  epoch.gpsSeconds = uvar->start;
  epoch.gpsNanoSeconds = 0;
  
  /*init ephemeris-data */
  EphemerisData *edat;
  if ( (edat = XLALInitEphemeris ( uvar->ephemYear )) == NULL ) {
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
  
  UINT4 ndata;
  REAL8 srate;
  
  ndata = uvar->duration;
  srate = uvar->srate;
  
  /*creates the .gwf channels */
  CHAR in_chan[256];
  sprintf( in_chan, "%s:LDAS-STRAIN", uvar->IFO );
  /*fprintf( stderr, "\nin_chan = %s\n", in_chan);*/

  CHAR out_chan[256];
  sprintf( out_chan, "%s_LDAS_C02_L2_CWINJ", uvar->IFO );
  /*fprintf( stderr, "\nout_chan = %s\n", out_chan);*/

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
  
  REAL8TimeSeries *gwfseries=NULL;
  REAL8TimeSeries *series=NULL;
  
  char out_file[256]; /*need a different way to do this */
  FILE *inject; /*for .par file from uvar->inputdir*/

  CHAR gwf_dir[256]; /* for use with XLALFrOpen */
  sprintf( gwf_dir, "%s/.", uvar->gwfdir );

  if ( (m = scandir(uvar->gwfdir, &namelist, 0, alphasort)) == -1) {
    XLALPrintError ("%s: scandir('%s',...) failed.\n", fn, uvar->gwfdir );
    XLAL_ERROR ( fn, XLAL_EIO );
  }
  UINT4 k=0;
  for (k=2; k < (UINT4)m; k++){
    fprintf(stderr, "\ngwf_dir = %s \nfile name = %s\n", uvar->gwfdir, namelist[k]->d_name );
    if (( gwffile = XLALFrOpen( uvar->gwfdir, namelist[k]->d_name)) == NULL ) {
      XLAL_ERROR ( fn, XLAL_EFUNC );
    }
    else {
      
      /****HERE WE WANT TO EXTRACT DURATION, EPOCH TO CREATE SERIES FROM .gwf FILE NAMES ****/
      
      if ( (gwfseries = XLALCreateREAL8TimeSeries( in_chan, &epoch, 0., 1./srate, &lalSecondUnit, (int)(ndata*srate) )) == NULL ) {
	XLAL_ERROR ( fn, XLAL_EFUNC );
      }

      if ( XLALFrGetREAL8TimeSeries( gwfseries, gwffile ) != XLAL_SUCCESS ) {
	XLAL_ERROR ( fn, XLAL_EFUNC );
      }
      
   
      /* define output .gwf file */
      sprintf(out_file, "%c-%s-%d-%d.gwf", uvar->IFO[0], out_chan, epoch.gpsSeconds, ndata);
      fprintf(stderr, "out_file = %s\n", out_file);

      /* read in and test generated frame with XLAL function*/

      /**** This is where we need to change to Chris's code for writing frames ********/
      /*This is wher the error is coming in... need to have a file defined and written already */
      if (( frfile = XLALFrOpen( uvar->outputdir, out_file )) == NULL){
	XLAL_ERROR ( fn, XLAL_EFUNC);
	}
      else {
      series=NULL;
	if ((series = XLALCreateREAL8TimeSeries( out_chan, &epoch, 0., 1./srate, &lalSecondUnit, (int)(ndata*srate) )) == NULL){
	  XLAL_ERROR ( fn, XLAL_EFUNC );
	}
	CHAR *fffile;
	sprintf(fffile, "%s/%s", uvar->outputdir, out_file);
	if (fopen(fffile , 'w+') == NULL){
	  fprintf(stderr, "fail");
	}

	if (( frfile = XLALFrOpen( uvar->outputdir, out_file )) == NULL){
	  XLAL_ERROR ( fn, XLAL_EFUNC );
	}
	/*create series to be the sum, general series to add to the noise */
	REAL8TimeSeries *total_inject=NULL;
	if (( total_inject = XLALCreateREAL8TimeSeries(out_chan, &epoch, 0., 1./srate, &lalSecondUnit, (UINT4)(ndata*srate)) ) == NULL) {
	  XLAL_ERROR ( fn, XLAL_EFUNC );
	}
	
	/*  Tstamp = XLALCreateREAL8Vector((UINT4)(ndata*srate));*/
	/*fprintf(stderr, "length = %d\n", injsig->length);*/
	
	/*Read in all pulsar files from inputdir*/
	int n = scandir(uvar->inputdir, &namelist, 0, alphasort);
	if ( n < 0 ) {
	  XLALPrintError ("%s: scandir() failed for directory '%s'\n", fn, uvar->inputdir );
	  XLAL_ERROR ( fn, XLAL_EIO );
	}
	UINT4 numParFiles = (UINT4)n;

	/*k starts at 2 to skip over . and .. found in scandir, they are the first two entries */
	for (k=2; k < numParFiles; k++) {
	  sprintf(pulin, "%s/%s", uvar->inputdir, namelist[k]->d_name);
	  fprintf(stderr, "This is your pulsar file:%s\n", pulin);
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
	  } /* if pulparams.model */
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
	} /* for k < numParFiles */
      
	/*  UINT4 i;*/
	/*  for (k=2; k < n; k++){*/
	/*fprintf( stderr, "%s\n",namelist[k]->d_name );*/
	/*    sprintf(fin, "%s/%s",inputdir, namelist[k]->d_name);*/
	/*    fprintf(stderr, "%s\n", fin);*/
	/*    if (( inject = fopen( fin, "r" )) ==NULL)*/
	/*      fprintf(stderr, "Error opening file\n" );*/
	/*    else{*/
	/*      double temp2=0.;*/
	
	/*      j=0;*/
	
	/*      While(!feof(inject)){*/
	/*        fscanf(inject,"%lf %lf", &Tstamp->data[j], &temp2);*/
	/*	injsig->data[j] += temp2;*/
	/*	j++;*/
	/*      }*/
	
	/*      for(i=0; i<3; i++)*/
	/*	fprintf(stderr, "%g\n", injsig->data[i]);*/
	/*      fclose(inject);*/
	/*    }*/
	/*  }*/
	/* junk */
	
	/*add to series*/
	for (i=0; i < series->data->length; i++){
	  series->data->data[i] = gwfseries->data->data[i] + total_inject->data->data[i]; /*read in makefakedata file -- how to add?*/
	  /*fprintf(stderr, "%le\n", series->data->data[i]);*/
	}
	fprintf(stderr,"%e\n", series->data->data[1]);
	
	if (( XLALFrWriteREAL8TimeSeries( series, 1 ) ) != XLAL_SUCCESS ){
	  XLAL_ERROR( fn, XLAL_EFUNC );
	}
	XLALDestroyREAL8TimeSeries( gwfseries);
	XLALDestroyREAL8TimeSeries( total_inject );
	XLALDestroyREAL8TimeSeries( series );
	/****Test for Matlab****/
	/*  if (( outtest = fopen( "/home/erinmacdonald/lsc/analyses/sw_injections/test.txt", "w+" )) == NULL)*/
	/*    fprintf(stderr, "Error opening file\n");*/
	/*  else{*/
	/*    for (i = 0; i < series->data->length; i++){*/
	/*      fprintf( outtest,"%le\n", series->data->data[i] );*/
	/*    }*/
	/*    fclose(outtest);*/
	/*  }*/
	/****End test****/
	
      }
	/* Free up the memory */
	
	/*  XLALDestroyREAL8Vector(injsig);*/
	
	/*  XLALDestroyREAL8TimeSeries( gwfseries );*/
	
	/*  XLALDestroyREAL8TimeSeries( series );*/
    } /*ends loop for creating CWINJ .gwf file */
  } /*ends loop through all .gwf files in .gwf directory*/

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
  XLALregBOOLUserStruct( help,            'h', UVAR_HELP, "Print this meessage");
  /*    XLALregSTRINGUserStruct(out_chan,   'o', UVAR_OPTIONAL, "Output channel i.e. (IFO)_LDAS_C02_L2_CWINJ");*/
  /*    XLALregSTRINGUserStruct(in_chan,        'i', UVAR_OPTIONAL, "Input channel from .gwf file, i.e. (IFO):LDAS-STRAIN");*/
  XLALregREALUserStruct(srate,            'r', UVAR_OPTIONAL, "user defined sample rate, default = 16384");
  XLALregREALUserStruct(duration,       'd', UVAR_OPTIONAL, "duration of frame (sec)");
  XLALregREALUserStruct(start,            's', UVAR_OPTIONAL, "epoch in GPS Seconds");
  XLALregSTRINGUserStruct(inputdir,       'p', UVAR_OPTIONAL, "directory for .par files");
  XLALregSTRINGUserStruct(gwfdir,     'g', UVAR_OPTIONAL,"directory for .gwf files");
  XLALregSTRINGUserStruct(outputdir,  'c', UVAR_OPTIONAL, "directory for CWINJ files");
  XLALregSTRINGUserStruct(ephemYear,      'y', UVAR_OPTIONAL,"Year(or range of years) for ephemeris files to be used");
  XLALregSTRINGUserStruct(ephemDir,   'e', UVAR_OPTIONAL,"Directory to find ephemeris files");
  XLALregSTRINGUserStruct( IFO,       'I', UVAR_REQUIRED, "Detector: 'G1', 'L1', 'H1', 'H2', 'V1'...");
  
  if (XLALUserVarReadAllInput (argc, argv ) != XLAL_SUCCESS) {
    XLALPrintError ("%s: XLALUserVarReadAllInput() failed with errno=%d\n", fn, xlalErrno);
    XLAL_ERROR ( fn, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

} /* InitUserVars() */


EphemerisData *
XLALInitEphemeris (const CHAR *ephemYear )
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
  
  snprintf(EphemEarth, FNAME_LENGTH, "earth%s.dat", ephemYear);
  snprintf(EphemSun, FNAME_LENGTH, "sun%s.dat",  ephemYear);  

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
