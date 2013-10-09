/*  Copyright (C) 2013 Chris Messenger
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

/** \author C.Messenger
 * \ingroup pulsarApps
 * \file
 * \brief
 * This code is designed to convert a binary timeseries file into a frame file.
 *
 */

/***********************************************************************************************/
/* includes */
#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <lal/TimeSeries.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/SFTutils.h>
#include <lal/SFTfileIO.h>
#include <lal/ComplexFFT.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lalapps.h>
#include <lal/BandPassTimeSeries.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/LALConstants.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>

#include <lal/TimeSeries.h>
#include <lal/GeneratePulsarSignal.h>
#include <lal/TransientCW_utils.h>
#include <lal/LALString.h>
#include <lal/StringVector.h>

#include <gsl/gsl_rng.h>           /* for random number generation */ 
#include <gsl/gsl_randist.h>       /* for random number generation */ 

/** A structure that stores user input variables 
 */
typedef struct { 
  BOOLEAN help;		            /**< trigger output of help string */
  CHAR *outLabel;                   /**< 'misc' entry in SFT-filenames or 'description' entry of frame filenames */
  CHAR *outputdir;                  /**< the output directory */
  CHAR *cachefile;                  /**< the name of the input cache file */
  REAL8 fsamp;                      /**< the sampling frequency of the data */
  INT4 tsft;                        /**< the length of the SFTs */
  REAL8 freq;                       /**< the starting frequency */
  REAL8 freqband;                   /**< the band width */
  REAL8 highpassf;                  /**< the high pass filter frequency */
  BOOLEAN outSingleSFT;	            /**< use to output a single concatenated SFT */
  REAL8 amp_inj;                /**< if set we inject a fake signal with this fractional amplitude */
  INT4 seed;
  REAL8 f_inj;
  REAL8 asini_inj;
  REAL8 tasc_inj;
  REAL8 P_inj;
  REAL8 phi_inj;
} UserInput_t;

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
int XLALReadUserVars(int argc,char *argv[],UserInput_t *uvar);
int compare_function(const void *a,const void *b);
int XLALInitgslrand(gsl_rng **gslrnd,INT8 seed);
int XLALCopySFT (SFTtype *dest, const SFTtype *src);
int XLALAppendSFT2Vector (SFTVector *vect,const SFTtype *sft);
			  
/***********************************************************************************************/
/* empty initializers */
UserInput_t empty_UserInput;
static const LALUnit empty_LALUnit;

/** The main function of binary2sft.c
 *
 */
int main( int argc, char *argv[] )  {

  UserInput_t uvar = empty_UserInput;           /* user input variables */
  INT4 i,j,k;                                     /* counter */
  INT4 sftcnt = 0;                              /* a counter for the number of SFTs */

  /**********************************************************************************/
  /* register and read all user-variables */
  if (XLALReadUserVars(argc,argv,&uvar)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadUserVars() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : read in uservars\n",__func__);

  /**********************************************************************************/
  /* read in the cache file */
  FILE *cachefp = NULL;
  if ((cachefp = fopen(uvar.cachefile,"r")) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to open binary input file %s\n",__func__,uvar.cachefile);
    return 1;
  }
  i = 0;
  while (fscanf(cachefp,"%*s %*d %*d")!=EOF) i++;
  INT4 Nfiles = i;
  fclose(cachefp);
  LogPrintf(LOG_DEBUG,"%s : counted %d files listed in the cache file.\n",__func__,Nfiles);

  /* allocate memory */
  char **filenames = LALCalloc(Nfiles,sizeof(char*));
  LIGOTimeGPSVector fileStart;
  fileStart.data = LALCalloc(Nfiles,sizeof(LIGOTimeGPS));
  for (i=0;i<Nfiles;i++) filenames[i] = LALCalloc(512,sizeof(char));
  
  if ((cachefp = fopen(uvar.cachefile,"r")) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to open binary input file %s\n",__func__,uvar.cachefile);
    return 1;
  }

  for (i=0;i<Nfiles;i++) {
    fscanf(cachefp,"%s %d %d",filenames[i],&(fileStart.data[i].gpsSeconds),&(fileStart.data[i].gpsNanoSeconds));
  }
  fclose(cachefp);
  
  /* initialise results vector */
  SFTVector *SFTvect = LALCalloc(1,sizeof(SFTVector));

  /* initialise the random number generator */
  gsl_rng * r;
  if (XLALInitgslrand(&r,uvar.seed)) {
    LogPrintf(LOG_CRITICAL,"%s: XLALinitgslrand() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EFAULT);
  }  

  /**********************************************************************************/
  /* loop over the input files */
  for (j=0;j<Nfiles;j++) {

    LogPrintf(LOG_DEBUG,"%s : working on file %s\n",__func__,filenames[j]);
    
    /**********************************************************************************/
    /* open the binary file and count the number of elements */
    FILE *binfp = NULL;                              /* file pointer for input file */
    REAL4 dummy;                                  /* dummy variable */
    if ((binfp = fopen(filenames[j],"r")) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to open binary input file %s\n",__func__,filenames[j]);
      return 1;
    }
    i = 0;
    while (fread(&dummy,sizeof(REAL4),1,binfp)) i++;
    INT4 N = i;
    fclose(binfp);
    LogPrintf(LOG_DEBUG,"%s : counted %d samples in the file.\n",__func__,N);

    /* NOTE: a timeseries of length N*dT has no timestep at N*dT !! (convention) */
    REAL8 dt = 1.0 / uvar.fsamp;
    LIGOTimeGPS startTimeGPS;                     /* the start time of the observation */
    startTimeGPS.gpsSeconds = fileStart.data[j].gpsSeconds;
    startTimeGPS.gpsNanoSeconds = fileStart.data[j].gpsNanoSeconds;
    LIGOTimeGPS endTimeGPS = startTimeGPS;
    XLALGPSAdd(&endTimeGPS,N*dt);
    REAL4TimeSeries *Tseries = XLALCreateREAL4TimeSeries ( "X1", &(startTimeGPS), 0, dt, &empty_LALUnit, N);
    REAL8 sum = 0;

    /* read in the data to the timeseries */
    if ((binfp = fopen(filenames[j],"r")) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to open binary input file %s\n",__func__,filenames[j]);
      return 1;
    }
    i = 0;
    while (fread(&dummy,sizeof(REAL4),1,binfp)) {
      Tseries->data->data[i] = (REAL4)dummy;
      sum += Tseries->data->data[i];
      /* fprintf(stdout,"%f\n",Tseries->data->data[i]); */
      i++;
    }
    REAL8 rate = sum/N;    /* the rate per bin */
    /* fprintf(stdout,"sft %d : counts/sec = %f\n",j,sum/(Tseries->data->length*dt)); */
    fclose(binfp);
  
    /**********************************************************************************/
    /* make timestamps */
    LIGOTimeGPSVector timestamps;
    INT4 Nsft = (INT4)floor(N*dt/uvar.tsft);        /* compute the number of SFTs */
    timestamps.length = Nsft;
    
    /* only continue if we have any data in this SFT */
    if (Nsft>0) {

      timestamps.data = LALCalloc(Nsft,sizeof(LIGOTimeGPS));
      for (i=0;i<Nsft;i++) {
	timestamps.data[i] = startTimeGPS;
	XLALGPSAdd(&(timestamps.data[i]),i*uvar.tsft);
      }
      
      /* if we are injecting a signal */
      if (uvar.amp_inj>0) {
	
	REAL8 nu = uvar.f_inj;
	REAL8 a = uvar.asini_inj;
	REAL8 tasc = uvar.tasc_inj;
	REAL8 omega = LAL_TWOPI/uvar.P_inj;
	REAL8 phi0 = uvar.phi_inj;

	/* loop over time */
	for (k=0;k<(INT4)Tseries->data->length;k++) {
	  	  
	  REAL8 t = k*dt + XLALGPSGetREAL8(&startTimeGPS);   
	  REAL8 phase = phi0 + LAL_TWOPI*nu*(t - a*sin(omega*(t-tasc)));
	  REAL8 lambda = rate + uvar.amp_inj*rate*sin(phase);
	  INT4 temp = gsl_ran_poisson (r, lambda);
	  Tseries->data->data[k] = (REAL4)temp;
	  
	}

      }

      /**********************************************************************************/
      /* define the high pass filter params and high pass filter the data */
      PassBandParamStruc filterpar;
      char tmpname[] = "Butterworth High Pass";
      filterpar.name  = tmpname;
      filterpar.nMax  = 10;
      filterpar.f2    = uvar.highpassf;
      filterpar.a2    = 0.5;
      filterpar.f1    = -1.0;
      filterpar.a1    = -1.0;
      XLALButterworthREAL4TimeSeries(Tseries, &filterpar);
      
      /**********************************************************************************/
      /* compute SFTs from timeseries */
      SFTParams sftParams = empty_SFTParams;
      sftParams.Tsft = uvar.tsft;
      sftParams.timestamps = &timestamps;
      sftParams.noiseSFTs = NULL;       // not used here any more!
      sftParams.window = NULL;
      SFTVector *sftVect;
      XLAL_CHECK ( (sftVect = XLALSignalToSFTs (Tseries, &sftParams)) != NULL, XLAL_EFUNC );
      
      /*  for (k=0;k<sftVect->data[0].data->length;k++) { */
      /*       fprintf(stdout,"%f %f\n",crealf(sftVect->data[0].data->data[k]),cimagf(sftVect->data[0].data->data[k])); */
      /*     } */
      /*     exit(0); */
      
      /**********************************************************************************/
      /* extract effective band from this, if neccessary (ie if faster-sampled output SFTs) */
      SFTVector *tempSFTvect = XLALExtractBandFromSFTVector ( sftVect, uvar.freq, uvar.freqband );
      XLALDestroySFTVector ( sftVect ); 
      
      /**********************************************************************************/
      /* append these SFTs to the full list of SFTs */
      for (k=0;k<Nsft;k++) XLALAppendSFT2Vector(SFTvect,&(tempSFTvect->data[k]));  
      sftcnt += Nsft;
      
      /**********************************************************************************/
      
      XLALFree(timestamps.data);
      XLALDestroySFTVector (tempSFTvect);

    }
    
    /* free memory inside the loop */
    XLALDestroyREAL4TimeSeries (Tseries);
    
  }  /* end loop over input files */

 /*  for (i=0;i<sftcnt;i++) { */
/*     fprintf(stdout,"%d %d\n",SFTvect->data[i].epoch.gpsSeconds,SFTvect->data[i].epoch.gpsNanoSeconds); */
/*   } */

  /**********************************************************************************/
  /* generate comment string */
  char *VCSInfoString = XLALGetVersionString(0);
  XLAL_CHECK ( VCSInfoString != NULL, XLAL_EFUNC, "XLALGetVersionString(0) failed.\n" );
  CHAR *logstr;
  size_t len; 
  XLAL_CHECK ( (logstr = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );
  char *comment = XLALCalloc ( 1, len = strlen ( logstr ) + strlen(VCSInfoString) + 512 );
  XLAL_CHECK ( comment != NULL, XLAL_ENOMEM, "XLALCalloc(1,%d) failed.\n", len );
  sprintf ( comment, "Generated by:\n%s\n%s\n", logstr, VCSInfoString );
  
  /**********************************************************************************/
  /* either write whole SFT-vector to single concatenated file */
  if ( uvar.outSingleSFT ) {
    XLAL_CHECK ( XLALWriteSFTVector2File( SFTvect, uvar.outputdir, comment, uvar.outLabel ) == XLAL_SUCCESS, XLAL_EFUNC );
  } else {	/* or as individual SFT-files */
    XLAL_CHECK ( XLALWriteSFTVector2Dir( SFTvect, uvar.outputdir, comment, uvar.outLabel ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  
  /**********************************************************************************/
  /* free memory */
  XLALDestroySFTVector ( SFTvect );
  XLALFree ( logstr );
  XLALFree ( comment );
 
  LALCheckMemoryLeaks();

  return 0;

}

/*******************************************************************************/

/** Read in input user arguments
 *
 */
int XLALReadUserVars(int argc,            /**< [in] the command line argument counter */ 
		     char *argv[],        /**< [in] the command line arguments */
		     UserInput_t *uvar    /**< [out] the user input structure */
		     )
{


  /* initialise user variables */
  uvar->outLabel = NULL; 
  uvar->outputdir = NULL;
  uvar->cachefile = NULL;
  uvar->fsamp = 2048;
  uvar->tsft = 100;
  uvar->freq = 550;
  uvar->freqband = 1;
  uvar->highpassf = 40;
  uvar->outSingleSFT = 0;
  uvar->amp_inj = 0;
  uvar->f_inj = 550.2;
  uvar->asini_inj = 2.0;
  uvar->tasc_inj = 800000000;
  uvar->P_inj = 68212;
  uvar->phi_inj = 1.0;
  uvar->seed = 0;

  /* ---------- register all user-variables ---------- */
  XLALregBOOLUserStruct(help, 		        'h', UVAR_HELP,     "Print this message");
  XLALregSTRINGUserStruct(outLabel, 	        'n', UVAR_REQUIRED, "'misc' entry in SFT-filenames or 'description' entry of frame filenames"); 
  XLALregSTRINGUserStruct(outputdir, 	        'o', UVAR_REQUIRED, "The output directory name"); 
  XLALregSTRINGUserStruct(cachefile, 	        'i', UVAR_REQUIRED, "The input binary file name"); 
  XLALregREALUserStruct(freq,                   'f', UVAR_OPTIONAL, "The starting frequency (Hz)");
  XLALregREALUserStruct(freqband,   	        'b', UVAR_OPTIONAL, "The frequency band (Hz)");
  XLALregINTUserStruct(tsft,                    't', UVAR_OPTIONAL, "The length of SFTs (sec)");
  XLALregREALUserStruct(fsamp,           	's', UVAR_OPTIONAL, "The sampling frequency");
  XLALregREALUserStruct(highpassf,           	'p', UVAR_OPTIONAL, "The high pass filter frequency");
  XLALregBOOLUserStruct(outSingleSFT,           'S', UVAR_OPTIONAL, "Write a single concatenated SFT file instead of individual files" );
  XLALregREALUserStruct(amp_inj,          	'I', UVAR_OPTIONAL, "Fractional amplitude of injected signal");
  XLALregREALUserStruct(f_inj,            	'A', UVAR_OPTIONAL, "frequency of injected signal");
  XLALregREALUserStruct(asini_inj,          	'B', UVAR_OPTIONAL, "projected semi-major axis of injected signal");
  XLALregREALUserStruct(tasc_inj,          	'C', UVAR_OPTIONAL, "time of ascension of injected signal");
  XLALregREALUserStruct(P_inj,           	'D', UVAR_OPTIONAL, "orbital period of injected signal");
  XLALregREALUserStruct(phi_inj,          	'E', UVAR_OPTIONAL, "initial phase of injected signal");
  XLALregINTUserStruct(seed,                    'r', UVAR_OPTIONAL, "The random seed");

  /* do ALL cmdline and cfgfile handling */
  if (XLALUserVarReadAllInput(argc, argv)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALUserVarReadAllInput() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EINVAL);
  }

  /* if help was requested, we're done here */
  if (uvar->help) exit(0);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/** this function initialises the gsl random number generation
 *
 * If the input seed is zero then a random seed is drawn from
 * /dev/urandom.
 *
 */
int XLALInitgslrand(gsl_rng **gslrnd,     /**< [out] the gsl random number generator */
		    INT8 seed             /**< [in] the random number generator seed */
		    )
{
  FILE *devrandom = NULL;      /* pointer to the /dev/urandom file */
  
  /* if the seed is 0 then we draw a random seed from /dev/urandom */
  if (seed == 0) {
    
    /* open /dev/urandom */
    if ((devrandom=fopen("/dev/urandom","r")) == NULL)  {
      LogPrintf(LOG_CRITICAL,"%s: Error, unable to open device /dev/random\n",__func__);
      XLAL_ERROR(XLAL_EINVAL);
    }
    
    /* read a random seed */
    if (fread((void*)&seed,sizeof(INT8),1,devrandom) != 1) {
      LogPrintf(LOG_CRITICAL,"%s: Error, unable to read /dev/random\n",__func__);
      XLAL_ERROR(XLAL_EINVAL);
    }
    fclose(devrandom);
    
  }
  
  /* setup gsl random number generation */
  *gslrnd = gsl_rng_alloc(gsl_rng_taus2);
  gsl_rng_set(*gslrnd,seed);
 
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/** Append the given SFTtype to the SFT-vector (no SFT-specific checks are done!) */
int XLALAppendSFT2Vector (SFTVector *vect,		/**< destinatino SFTVector to append to */
			  const SFTtype *sft            /**< the SFT to append */
			  )	
{
  UINT4 oldlen = vect->length;

  if ( (vect->data = LALRealloc ( vect->data, (oldlen + 1)*sizeof( *vect->data ) )) == NULL ) {
     LogPrintf(LOG_CRITICAL,"%s: Error, unable to allocate memory\n",__func__);
     XLAL_ERROR(XLAL_EINVAL);
  }
  memset ( &(vect->data[oldlen]), 0, sizeof( vect->data[0] ) );
  vect->length ++;

  XLALCopySFT(&vect->data[oldlen], sft );
  
  return XLAL_SUCCESS;
  
} /* XLALAppendSFT2Vector() */

/** Copy an entire SFT-type into another.
 * We require the destination-SFT to have a NULL data-entry, as the
 * corresponding data-vector will be allocated here and copied into
 *
 * Note: the source-SFT is allowed to have a NULL data-entry,
 * in which case only the header is copied.
 */
int XLALCopySFT (SFTtype *dest, 	/**< [out] copied SFT (needs to be allocated already) */
		 const SFTtype *src	/**< input-SFT to be copied */
		 )
{

  /* copy complete head (including data-pointer, but this will be separately alloc'ed and copied in the next step) */
  memcpy ( dest, src, sizeof(*dest) );

  /* copy data (if there's any )*/
  if ( src->data )
    {
      UINT4 numBins = src->data->length;
      if ( (dest->data = XLALCreateCOMPLEX8Vector ( numBins )) == NULL ) {
	LogPrintf(LOG_CRITICAL,"%s: Error, unable to allocate memory\n",__func__);
	XLAL_ERROR(XLAL_EINVAL);
      }
      memcpy (dest->data->data, src->data->data, numBins * sizeof (src->data->data[0]));
    }

   return XLAL_SUCCESS;

} /* XLALCopySFT() */


