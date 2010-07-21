/*
 * Copyright (C) 2010 Chris Messenger
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
 * \ingroup pulsarCoherent
 * \file
 * \brief
 * This code is designed to compute statistical quantities from a specified channel in a given 
 * frame file.
 *
 * It reads in the timeseries from the frame file and returns the quantiles of the timeseries 
 * distribution as well as the min, max, mean and standard deviation values.
 *
 */

/***********************************************************************************************/
/* includes */
#include <math.h>
#include <glob.h>
#include <lal/TimeSeries.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/LALFrameIO.h>
#include <lal/SFTutils.h>
#include <lal/FrameStream.h>
#include <lalappsfrutils.h>
#include <lalapps.h>

/* remove me */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/Date.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALFrameIO.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/LALDatatypes.h>
#include <lal/FrameCache.h>

/***********************************************************************************************/
/* some macros */

/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define NQUANTILE 10                       /* the number of quantiles to compute */
#define STRINGLENGTH 1024

/***********************************************************************************************/
/* define internal structures */

/** A structure that stores statistical information gathered from a timeseries
 */
typedef struct { 
  CHAR filename[1024];              /**< the frame filename */
  INT4 epoch;                       /**< the GPS start time */
  REAL8 duration;                   /**< the duration */
  REAL8 dt;                         /**< the sampling time */
  REAL8 min;                        /**< the minimum sample in the timeseries */
  REAL8 max;                        /**< the maximum sample in the timeseries */
  REAL8 median;                     /**< the median of the samples in the timeseries */
  REAL8 mean;                       /**< the mean of the samples in the timeseries */
  REAL8 var;                        /**< the variance of the samples in the timeseries */
  REAL8 quant[NQUANTILE+1];         /**< the quantiles of the samples in the timeseries */
} Stats;

/** A structure that stores statistical information gathered from a timeseries
 */
typedef struct { 
  Stats *data;                      /**< a vector of stats info */
  INT4 length;                      /**< the length of the vector */
} StatsVector;

/** A structure that sores user input variables 
 */
typedef struct { 
  BOOLEAN help;		            /**< trigger output of help string */
  CHAR *inputfile;                  /**< input frame file name pattern */
  CHAR *channel;                    /**< the frame channel to read */
  CHAR *outputfile;                 /**< name of output file */
  BOOLEAN version;	            /**< output version-info */
} UserInput_t;

/***********************************************************************************************/
/* global variables */
extern int vrbflg;	 	/**< defined in lalapps.c */
RCSID( "$Id$");		/* FIXME: use git-ID instead to set 'rcsid' */

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
void ReadUserVars(LALStatus *status,int argc,char *argv[],UserInput_t *uvar);
int XLALReadFrameDir(glob_t *filelist, CHAR *inputfile);
int XLALReadFrameINT4TimeSeries(INT4TimeSeries **ts,CHAR *filename,CHAR *channel);
int XLALComputeINT4TimeSeriesStats(Stats *stats,INT4TimeSeries *ts,CHAR *filename);
static int compareINT4(const void *p1, const void *p2);
int XLALOutputStats(StatsVector *stats, CHAR *outputfile);

/***********************************************************************************************/
/* empty initializers */
UserInput_t empty_UserInput;

/** The main function of framestats.c
 *
 * Here we read in a channel from a single frame file and return statistical quantities 
 * computed from the timeseries.  
 *
 */
int main( int argc, char *argv[] )  {

  static const char *fn = __func__;             /* store function name for log output */
  LALStatus status = blank_status;              /* empty LAL status structure */
  UserInput_t uvar = empty_UserInput;           /* user input variables */
  glob_t filelist;                              /* stores the matching frame file names */
  StatsVector stats;                            /* a structure containing statistical info */
  INT4 i;                                       /* counter */

  lalDebugLevel = 1;
  vrbflg = 1;	                        /* verbose error-messages */

  /* turn off default GSL error handler */
  gsl_set_error_handler_off();

  /* setup LAL debug level */
  LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LogSetLevel(lalDebugLevel);

  /* register and read all user-variables */
  LAL_CALL (ReadUserVars(&status,argc,argv,&uvar), &status);
  LogPrintf(LOG_DEBUG,"%s : read in uservars\n",fn);
 
  /**********************************************************************************/
  /* READ FILES */
  /**********************************************************************************/

  /* get a list of frame file names */
  if (XLALReadFrameDir(&filelist,uvar.inputfile)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadFrameDir() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }

  /* allocate memory for the stats output */
  stats.length = (INT4)filelist.gl_pathc;
  if ((stats.data = (Stats *)XLALCalloc(stats.length,sizeof(Stats))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to allocate memory for output stats info.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /**********************************************************************************/
  /* FOR EACH FILE READ DATA AND COMPUTE STATS */
  /**********************************************************************************/

  /* for each frame file */
  for (i=0;i<(INT4)filelist.gl_pathc;i++) {

    INT4TimeSeries *ts = NULL;              /* a timeseries */
   
    /* read frame into a timeseries structure */
    if (XLALReadFrameINT4TimeSeries(&ts,filelist.gl_pathv[i],uvar.channel)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALReadFrameINT4TimeSeries() failed with error = %d\n",fn,xlalErrno);
      return 1;
    }
    
    /* compute timeseries statistics */   
    if (XLALComputeINT4TimeSeriesStats(&(stats.data[i]),ts,filelist.gl_pathv[i])) {
      LogPrintf(LOG_CRITICAL,"%s : XLALComputeINT4TimeSeriesStats() failed with error = %d\n",fn,xlalErrno);
      return 1;
    }
      
    /* free the timeseries */
    XLALDestroyINT4TimeSeries(ts);

  } /* end the loop over files */
    
  /**********************************************************************************/
  /* OUTPUT RESULTS */
  /**********************************************************************************/
  
  if (XLALOutputStats(&stats,uvar.outputfile)) { 
    LogPrintf(LOG_CRITICAL,"%s : XLALOutputStats() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }

  /**********************************************************************************/
  /* CLEAN UP */
  /**********************************************************************************/
  
  /* free filelist */
  globfree(&filelist);
  
  /* free stats */
  XLALFree(stats.data);

  /* Free config-Variables and userInput stuff */
  LAL_CALL (LALDestroyUserVars (&status), &status);

  /* did we forget anything ? */
  LALCheckMemoryLeaks();
  LogPrintf(LOG_DEBUG,"%s : successfully checked memory leaks.\n",fn);

  LogPrintf(LOG_DEBUG,"%s : successfully completed.\n",fn);
  return 0;
  
}

/** Read in input user arguments
 */
void ReadUserVars(LALStatus *status,int argc,char *argv[],UserInput_t *uvar)
{
  
  CHAR *version_string;

  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* initialise user variables */
  uvar->inputfile = NULL; 
  uvar->outputfile = NULL; 
  uvar->channel = NULL;

  /* ---------- register all user-variables ---------- */
  LALregBOOLUserStruct(status, 	help, 		'h', UVAR_HELP,     "Print this message");
  LALregSTRINGUserStruct(status,inputfile, 	'i', UVAR_REQUIRED, "The input frame name pattern"); 
  LALregSTRINGUserStruct(status,outputfile, 	'o', UVAR_REQUIRED, "The output statistics file"); 
  LALregSTRINGUserStruct(status,channel, 	'c', UVAR_REQUIRED, "The frame channel to be read"); 
  LALregBOOLUserStruct(status,	version,        'V', UVAR_SPECIAL,  "Output code version");

  /* do ALL cmdline and cfgfile handling */
  LAL_CALL (LALUserVarReadAllInput(status->statusPtr, argc, argv), status->statusPtr);
 
  /* if help was requested, we're done here */
  if (uvar->help) exit(0);

  if ((version_string = XLALGetVersionString(0)) == NULL) {
    XLALPrintError("XLALGetVersionString(0) failed.\n");
    exit(1);
  }
  
  if (uvar->version) {
    printf("%s\n",version_string);
    exit(0);
  }
  XLALFree(version_string);

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

/** Read in list of frame files from input directory
 */
int XLALReadFrameDir(glob_t *filelist,         /**< [out] a structure containing a list of all input frame channels */
		     CHAR *inputfile           /**< [in] the input frame file pattern */
		     )
{
  
  const CHAR *fn = __func__;      /* store function name for log output */
  INT4 i;                         /* counter */

  /* check input arguments */
  if (filelist == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input file list structure == NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }  
  if (inputfile == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input file string == NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }  
  printf("%s\n",inputfile);
  /* search for matching filenames */
  if (glob(inputfile,0,NULL,filelist)) {
    LogPrintf(LOG_CRITICAL,"%s : glob() failed to return a filelist.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }

  if (filelist->gl_pathc>0) {
    for (i=0;i<(INT4)filelist->gl_pathc;i++) LogPrintf(LOG_DEBUG,"%s : found file %s\n",fn,filelist->gl_pathv[i]);
  }
  else {
    LogPrintf(LOG_CRITICAL,"%s : could not find any frame files matching %s.\n",fn,inputfile);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
 
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}

/** this function combines the files listed in the combination plan into a single timeseries 
 */
int XLALReadFrameINT4TimeSeries(INT4TimeSeries **ts,           /**< [out] the timeseries */
				CHAR *filename,                /**< [in] the input frame file name */
				CHAR *channel                  /**< [in] the channel to be read */
				)
{

  const CHAR *fn = __func__;      /* store function name for log output */
  FrStream *fs = NULL;
  LIGOTimeGPS epoch;
  REAL8 duration;

  /* check input arguments */
  if ((*ts) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output INT4TimeSeries structure != NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EFAULT);
  } 
  if (filename == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input filename string == NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EFAULT);
  } 
  if (channel == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input frame channel string == NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }  
  LogPrintf(LOG_DEBUG,"%s : checked input\n",fn);
  
  /* if we can open the frame file */
  if ((fs = XLALFrOpen(NULL,filename)) != NULL) {
   
    LogPrintf(LOG_DEBUG,"%s: opened frame file %s.\n",fn,filename);
    
    /* define start and duration */
    XLALGPSSetREAL8(&epoch,(REAL8)fs->flist->t0);
    duration = fs->flist->dt;
    
    /* seek to the start of the frame */
    if (XLALFrSeek(fs,&epoch)) {
      LogPrintf(LOG_CRITICAL,"%s: unable to seek to start of frame file %s.\n",fn,filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    
    /* read in timeseries from this file - final arg is limit on length of timeseries (0 = no limit) */
    if (((*ts) = XLALFrReadINT4TimeSeries(fs,channel,&epoch,duration,0)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to read channel %s from frame file %s.\n",fn,filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    LogPrintf(LOG_DEBUG,"%s: reading channel %s\n",fn,channel);
  
    /* close the frame file */
    XLALFrClose(fs);
    
  }
  /* otherwise allocate space for an empty timeseries and extract frame params from the filename */
  else {
    LogPrintf(LOG_DEBUG,"%s: unable to open frame file %s.\n",fn,filename);
    xlalErrno = 0;
    INT4 i;

    /* extract epoch and duration from filename */
    {
      REAL8 newepoch;
      CHAR *c1 = strrchr(filename,'-');
      INT4 length = strlen(c1) - 5;
      CHAR *tempdur = (CHAR *)XLALCalloc(length+1,sizeof(CHAR));
      CHAR *tempepoch = (CHAR *)XLALCalloc(10,sizeof(CHAR));
      strncpy(tempdur,c1+1,length);
      strncpy(tempepoch,c1-9,9);
      newepoch = atof(tempepoch);
      XLALGPSSetREAL8(&epoch,newepoch);
      duration = atof(tempdur);
      XLALFree(tempepoch);
      XLALFree(tempdur);
    }

    if (((*ts) = XLALCreateINT4TimeSeries("BAD_FILE",&epoch,0,1,&lalDimensionlessUnit,NQUANTILE)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to allocate empty timeseries\n.",fn);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
    for (i=0;i<(INT4)(*ts)->data->length;i++) (*ts)->data->data[i] = 0;

  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}

/** this function computes a series of statistics based on the input timeseries
 */
int XLALComputeINT4TimeSeriesStats(Stats *stats,             /**< [out] the timeseries statistics */
				   INT4TimeSeries *ts,       /**< [in] the input timeseries */ 
				   CHAR *filename            /**< [in] the input filename */
				   )
{  

  const CHAR *fn = __func__;      /* store function name for log output */
  INT8 N;                         /* number of samples */
  INT8 i;                         /* counter */
  REAL8 sum = 0.0;                /* initialise sum */
  REAL8 sqsum = 0.0;              /* initialise sum of squares  */

  /* check input arguments */
  if (stats == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output statistics structure == NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EFAULT);
  } 
  if (ts == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input timeseries == NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EFAULT);
  } 
  LogPrintf(LOG_DEBUG,"%s : checked input\n",fn);
  
  N = ts->data->length;
  if (N < 1) {
    LogPrintf(LOG_CRITICAL,"%s: length of timeseries < 1.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* record the epoch and duration */
  snprintf(stats->filename,STRINGLENGTH,"%s",filename);
  stats->epoch = GPS2REAL8(ts->epoch);
  stats->duration = (REAL8)(ts->deltaT*ts->data->length);
  stats->dt = ts->deltaT;

  /* compute mean */
  for (i=0;i<N;i++) sum += ts->data->data[i];
  stats->mean = sum/(REAL8)N;
  LogPrintf(LOG_DEBUG,"%s: computed mean as %f\n",fn,stats->mean);
 
  /* compute variance */
  for (i=0;i<N;i++) sqsum += (ts->data->data[i]-stats->mean)*(ts->data->data[i]-stats->mean);
  stats->var = sqsum/(REAL8)N;
  LogPrintf(LOG_DEBUG,"%s: computed variance as %f\n",fn,stats->var);

  /* compute min and max */
  stats->min = ts->data->data[0];
  stats->max = ts->data->data[0];
  for (i=0;i<N;i++) {
    if (ts->data->data[i]>stats->max) stats->max = ts->data->data[i];
    if (ts->data->data[i]<stats->min) stats->min = ts->data->data[i];
  }
  LogPrintf(LOG_DEBUG,"%s: computed min/max as %f/%f\n",fn,stats->min,stats->max);
  
  /* sort timeseries */
  qsort(ts->data->data, N, sizeof(INT4), compareINT4);
  
  /* compute the median */
  stats->median = ts->data->data[N/2];
  LogPrintf(LOG_DEBUG,"%s: computed median as %f\n",fn,stats->median);
  
  /* compute quantiles at 10% intervals */
  for (i=0;i<=NQUANTILE;i++) {

    INT8 idx = (INT8)(N*i/NQUANTILE);
    if (idx>N-1) idx = N - 1;

    stats->quant[i] = ts->data->data[idx];
    LogPrintf(LOG_DEBUG,"%s: computed quantile %d as %f\n",fn,i,stats->quant[i]);
  }
 
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}

/** comparison function for use with qsort
 */
static int compareINT4(const void *a, const void *b)
{
  
  const INT4 *ia = (const INT4 *)a; 
  const INT4 *ib = (const INT4 *)b;
  return *ia  - *ib; 
  
  /* integer comparison: returns negative if b > a 
	and positive if a > b */
   
}
/** function to output the stats results to file
 */
int XLALOutputStats(StatsVector *stats,      /**< [in] the output stats results */
		    CHAR *outputfile         /**< [in] the output filename */
		    ) 
{
  
  const CHAR *fn = __func__;      /* store function name for log output */
  FILE *fp = NULL;                /* file pointer */
  INT4 i,j;                        /* counters */

  /* open output file for writing */
  if ((fp = fopen(outputfile,"w")) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to open output file %s.\n",fn,outputfile);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* loop over input frame files */
  fprintf(fp,"## filename start duration dt mean median min max var quantile1 quantile2 ... quantileN\n");
  for (i=0;i<stats->length;i++) {

    fprintf(fp,"%s %d %f %f %f %f %f %f %f ",stats->data[i].filename,stats->data[i].epoch,stats->data[i].duration,stats->data[i].dt,
	    stats->data[i].mean,stats->data[i].median,stats->data[i].min,stats->data[i].max,stats->data[i].var);
    for (j=0;j<NQUANTILE;j++) fprintf(fp,"%f ",stats->data[i].quant[j]);
    fprintf(fp,"\n");

  }

  /* close the file */
  fclose(fp);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}
