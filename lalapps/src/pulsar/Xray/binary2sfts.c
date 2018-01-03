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
 * \ingroup lalapps_pulsar_Xray
 * \file
 * \brief
 * This code is designed to convert a binary timeseries file into a frame file.
 *
 */

/***********************************************************************************************/
/* includes */
#include "config.h"
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

#include "SemiCoherent.h"

/** A structure that stores user input variables
 */
typedef struct {
  CHAR *outLabel;                   /**< 'misc' entry in SFT-filenames or 'description' entry of frame filenames */
  CHAR *outputdir;                  /**< the output directory */
  CHAR *cachefile;                  /**< the name of the input cache file */
  REAL8 tsamp;                      /**< the sampling time of the data */
  INT4 tsft;                        /**< the length of the SFTs */
  REAL8 freq;                       /**< the starting frequency */
  REAL8 freqband;                   /**< the band width */
  REAL8 highpassf;                  /**< the high pass filter frequency */
  BOOLEAN outSingleSFT;	            /**< use to output a single concatenated SFT */
  REAL8 amp_inj;                    /**< if set we inject a fake signal with this fractional amplitude */
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

/***********************************************************************************************/
/* empty initializers */
UserInput_t empty_UserInput;

/** The main function of binary2sft.c
 *
 */
int main( int argc, char *argv[] )  {

  UserInput_t uvar = empty_UserInput;           /* user input variables */
  INT4 i,j;                                     /* counter */
  SFTVector *SFTvect = NULL;
  char *noisestr = XLALCalloc(1,sizeof(char));

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
    fscanf(cachefp,"%s %d %d %*d",filenames[i],&(fileStart.data[i].gpsSeconds),&(fileStart.data[i].gpsNanoSeconds));
  }
  fclose(cachefp);

  /* initialise the random number generator */
  gsl_rng * r;
  if (XLALInitgslrand(&r,uvar.seed)) {
    LogPrintf(LOG_CRITICAL,"%s: XLALinitgslrand() failed with error = %d\n",__func__,xlalErrno);
    XLAL_ERROR(XLAL_EFAULT);
  }

  /* setup the binaryToSFT parameters */
  BinaryToSFTparams par;
  par.tsft = uvar.tsft;
  par.freq = uvar.freq;
  par.freqband = uvar.freqband;
  par.tsamp = uvar.tsamp;
  par.highpassf = uvar.highpassf;
  par.amp_inj = uvar.amp_inj;
  par.f_inj = uvar.f_inj;
  par.asini_inj = uvar.asini_inj;
  XLALGPSSetREAL8(&(par.tasc_inj),uvar.tasc_inj);
  par.tref = fileStart.data[0];
  par.P_inj = uvar.P_inj;
  par.phi_inj = uvar.phi_inj;
  par.r = r;

  /**********************************************************************************/
  /* loop over the input files */
  long int ntot = 0;
  for (j=0;j<Nfiles;j++) {

    UINT4 k = 0;
    INT8Vector *np = NULL;
    REAL8Vector *R = NULL;
    par.tstart = fileStart.data[j];
    REAL8 norm1 = par.tsamp/par.tsft;
    REAL8 norm2 = 1.0/(par.tsamp*par.tsft);
    UINT4 oldlen;
    if (SFTvect==NULL) oldlen = 0;
    else oldlen = SFTvect->length;
    LogPrintf(LOG_DEBUG,"%s : working on file %s\n",__func__,filenames[j]);

    if (XLALBinaryToSFTVector(&SFTvect,filenames[j],&(fileStart.data[j]),&par,&np,&R)) {
      LogPrintf(LOG_CRITICAL,"%s : failed to convert binary input file %s to sfts\n",__func__,filenames[j]);
      return 1;
    }
    if ((np!=NULL) && (R!=NULL)) {
      for (k=0;k<np->length;k++) {
        ntot += np->data[k];
        char temp[64];
        sprintf(temp,"%d %e %e\n",SFTvect->data[oldlen+k].epoch.gpsSeconds,(REAL8)np->data[k]*norm1,R->data[k]*norm2);
        noisestr = (char *)XLALRealloc(noisestr,sizeof(char)*(1+strlen(noisestr)+strlen(temp)));
        strcat(noisestr,temp);
      }
      XLALDestroyINT8Vector(np);
      XLALDestroyREAL8Vector(R);
    }

  }  /* end loop over input files */

  /**********************************************************************************/
  /* create a noise string */


  /**********************************************************************************/
  /* generate comment string */
  char *VCSInfoString = XLALGetVersionString(0);
  XLAL_CHECK ( VCSInfoString != NULL, XLAL_EFUNC, "XLALGetVersionString(0) failed.\n" );
  CHAR *logstr;
  size_t len;
  XLAL_CHECK ( (logstr = XLALUserVarGetLog ( UVAR_LOGFMT_CMDLINE )) != NULL, XLAL_EFUNC );
  char *comment = XLALCalloc ( 1, len = strlen ( logstr ) + strlen(VCSInfoString) + strlen(noisestr) + 512 );
  XLAL_CHECK ( comment != NULL, XLAL_ENOMEM, "XLALCalloc(1,%zd) failed.\n", len );
  sprintf ( comment, "Generated by:\n%s\n%s\nTotal number of photons = %ld\n%s\n", logstr, VCSInfoString, ntot, noisestr );

  /**********************************************************************************/
  /* either write whole SFT-vector to single concatenated file */
  if ( uvar.outSingleSFT ) {
    XLAL_CHECK ( XLALWriteSFTVector2File( SFTvect, uvar.outputdir, comment, uvar.outLabel ) == XLAL_SUCCESS, XLAL_EFUNC );
  } else {	/* or as individual SFT-files */
    XLAL_CHECK ( XLALWriteSFTVector2Dir( SFTvect, uvar.outputdir, comment, uvar.outLabel ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  /**********************************************************************************/
  /* free memory */
  XLALDestroySFTVector(SFTvect);
  XLALFree(logstr);
  XLALFree(comment);
  XLALFree(noisestr);

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
  uvar->tsamp = 1.0/8192;
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
  XLALRegisterUvarMember(outLabel, 	        STRING, 'n', REQUIRED, "'misc' entry in SFT-filenames or 'description' entry of frame filenames");
  XLALRegisterUvarMember(outputdir, 	        STRING, 'o', REQUIRED, "The output directory name");
  XLALRegisterUvarMember(cachefile, 	        STRING, 'i', REQUIRED, "The input binary file name");
  XLALRegisterUvarMember(freq,                   REAL8, 'f', OPTIONAL, "The starting frequency (Hz)");
  XLALRegisterUvarMember(freqband,   	        REAL8, 'b', OPTIONAL, "The frequency band (Hz)");
  XLALRegisterUvarMember(tsft,                    INT4, 't', OPTIONAL, "The length of SFTs (sec)");
  XLALRegisterUvarMember(tsamp,           	REAL8, 's', OPTIONAL, "The sampling time (sec)");
  XLALRegisterUvarMember(highpassf,           	REAL8, 'p', OPTIONAL, "The high pass filter frequency");
  XLALRegisterUvarMember(outSingleSFT,           BOOLEAN, 'S', OPTIONAL, "Write a single concatenated SFT file instead of individual files" );
  XLALRegisterUvarMember(amp_inj,          	REAL8, 'I', OPTIONAL, "Fractional amplitude of injected signal");
  XLALRegisterUvarMember(f_inj,            	REAL8, 'A', OPTIONAL, "frequency of injected signal");
  XLALRegisterUvarMember(asini_inj,          	REAL8, 'B', OPTIONAL, "projected semi-major axis of injected signal");
  XLALRegisterUvarMember(tasc_inj,          	REAL8, 'C', OPTIONAL, "time of ascension of injected signal");
  XLALRegisterUvarMember(P_inj,           	REAL8, 'D', OPTIONAL, "orbital period of injected signal");
  XLALRegisterUvarMember(phi_inj,          	REAL8, 'E', OPTIONAL, "initial phase of injected signal");
  XLALRegisterUvarMember(seed,                    INT4, 'r', OPTIONAL, "The random seed");

  /* do ALL cmdline and cfgfile handling */
  BOOLEAN should_exit = 0;
  if (XLALUserVarReadAllInput(&should_exit, argc, argv)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALUserVarReadAllInput failed with error = %d\n",__func__,xlalErrno);
    return XLAL_EFAULT;
  }
  if (should_exit) exit(1);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}
