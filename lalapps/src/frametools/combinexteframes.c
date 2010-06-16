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
 * This code is designed to combine (R)XTE frames from simultaneous observations.
 *
 * A particular observation using (R)XTE can have multiple timeseries generated 
 * simultaneaously each with a different "data mode".  This code combines the data
 * into a single timeseries using user-defined choices regarding the type of combination.
 *
 */

/***********************************************************************************************/
/* includes */
#include <math.h>
#include <dirent.h>
#include <lal/TimeSeries.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/LALFrameIO.h>
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
/* some global constants */
#define MAX_DT 0.000976562            /* the maximum sample time we will use (equiv max freq 512 Hz) */
#define MIN_DT 0.000244140625         /* the minimum sample time we will use (equiv max freq 2048 Hz) */
#define APIDLENGTH 5                  /* the length of APID HEX strings */
#define STRINGLENGTH 256              /* the length of general string */
#define LONGSTRINGLENGTH 1024         /* the length of general string */
#define NAPID 6                       /* the number of valid APIDs we can currently read */
#define MINFRAMELENGTH 1              /* the minimum duration of output frame file in seconds */
#define NPCU 5                        /* the number of PCUs on XTE */
#define PCUCOUNTTHRESH 50             /* threshold on the counts per PCU that determine operational status */

/***********************************************************************************************/
/* define internal structures */

/** A structure that sores user input variables 
 */
typedef struct { 
  UINT4 length;	                    /**< the number of files */
  CHAR dir[STRINGLENGTH];           /**< the directory containing the files */
  struct dirent **list;             /**< a file list structure */
} FileList;

/** A structure that sores user input variables 
 */
typedef struct { 
  UINT4 length;	                    /**< the number of intervals */
  LIGOTimeGPS *start;               /**< the GPS start time of an interval */
  LIGOTimeGPS *end;                 /**< the GPS end time of an interval */
  UINT2 *pcucount;                  /**< the number of operational PCUs */
} GoodPCUIntervals;

/** A structure that sores user input variables 
 */
typedef struct { 
  BOOLEAN help;		            /**< trigger output of help string */
  CHAR *inputdir;                   /**< directory containing input frame files */
  CHAR *outputdir;                  /**< name of output directory */
  REAL8 deltat;                     /**< the desired sampling time */
  BOOLEAN version;	            /**< output version-info */
} UserInput_t;

/***********************************************************************************************/
/* global variables */
extern int vrbflg;	 	/**< defined in lalapps.c */

/* fill in the HEX APID names */
const char *APID[6] = {"FS37","FS3b","FS3f","FS4f","XENO","FS46"};

/* the PCU channel names for FS46 fits files */
const char *PCUchannel[5] = {"XeCntPcu0","XeCntPcu1","XeCntPcu2","XeCntPcu3","XeCntPcu4"};

RCSID( "$Id$");		/* FIXME: use git-ID instead to set 'rcsid' */

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
void ReadUserVars(LALStatus *status,int argc,char *argv[],UserInput_t *uvar, CHAR *clargs);
int XLALReadFrameDir(FileList **framefiles, CHAR *inputdir);
int XLALFrameFilter(const struct dirent *x);
int XLALReadGoodPCUInterval(GoodPCUIntervals **pcu,FileList *framefiles);

/***********************************************************************************************/
/* empty initializers */
UserInput_t empty_UserInput;

/** The main function of xtefitstoframe.c
 *
 * Here we read in a single XTE FITS file containing PCA data, generate a timeseries,
 * barycenter the data if requested, and output as a frame file. 
 *
 */
int main( int argc, char *argv[] )  {

  static const char *fn = __func__;             /* store function name for log output */
  LALStatus status = blank_status;              /* empty LAL status structure */
  UserInput_t uvar = empty_UserInput;           /* user input variables */
  FileList *framefiles = NULL;                  /* list of input frame file names */
  GoodPCUIntervals *pcu = NULL;                 /* the operational pcu count */ 
  CHAR clargs[LONGSTRINGLENGTH];                /* store the command line args */
  INT4 i;                                       /* counter */

  lalDebugLevel = 1;
  vrbflg = 1;	                        /* verbose error-messages */

  /* turn off default GSL error handler */
  gsl_set_error_handler_off();

  /* setup LAL debug level */
  LAL_CALL (LALGetDebugLevel(&status, argc, argv, 'v'), &status);
  LogSetLevel(lalDebugLevel);

  /* register and read all user-variables */
  LAL_CALL (ReadUserVars(&status,argc,argv,&uvar,clargs), &status);
  LogPrintf(LOG_DEBUG,"%s : read in uservars\n",fn);
 
  /**********************************************************************************/
  /* READ FILES */
  /**********************************************************************************/

  /* get a list of frame file names */
  if (XLALReadFrameDir(&framefiles,uvar.inputdir)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadFrameDir() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }

  /* get good pcu time intervals based on FS46 files photon counts */
  if (XLALReadGoodPCUInterval(&pcu,framefiles)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadGoodPCUInterval() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }
  for (i=0;i<(INT4)pcu->length;i++) LogPrintf(LOG_DEBUG,"%s : pcu interals -> %d %d %d\n",fn,pcu->start[i].gpsSeconds,pcu->end[i].gpsSeconds,pcu->pcucount[i]);

  /* read in data from the input FITS file */
  /* if (XLALReadFITSFile(&fitsdata,uvar.inputfile)) { */
/*     LogPrintf(LOG_CRITICAL,"%s : XLALReadFitsFile() failed with error = %d\n",fn,xlalErrno); */
/*     return 1; */
/*   }    */
/*   LogPrintf(LOG_DEBUG,"%s : Read in fits data from file %s\n",fn,fitsdata->header->filename);   */
  
  /**********************************************************************************/
  /* CLEAN UP */
  /**********************************************************************************/

  /* free filelist */
  for (i=0;i<(INT4)framefiles->length;i++) free(framefiles->list[i]);
  free(framefiles->list);
  XLALFree(framefiles);
 
  /* free PCU count structure */
  XLALFree(pcu->start);
  XLALFree(pcu->end);
  XLALFree(pcu->pcucount);
  XLALFree(pcu);

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
void ReadUserVars(LALStatus *status,int argc,char *argv[],UserInput_t *uvar,CHAR *clargs)
{
  
  CHAR *version_string;
  INT4 i;

  INITSTATUS( status, "initUserVars", rcsid );
  ATTATCHSTATUSPTR (status);

  /* initialise user variables */
  uvar->inputdir = NULL; 
  uvar->outputdir = NULL; 
  uvar->deltat = MIN_DT;

  /* ---------- register all user-variables ---------- */
  LALregBOOLUserStruct(status, 	help, 		'h', UVAR_HELP,     "Print this message");
  LALregSTRINGUserStruct(status,inputdir, 	'i', UVAR_REQUIRED, "The input frame directory"); 
  LALregSTRINGUserStruct(status,outputdir, 	'o', UVAR_REQUIRED, "The output frame file directory"); 
  LALregREALUserStruct(status,  deltat,         't', UVAR_OPTIONAL, "The output sampling time (in seconds)");
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

  /* put clargs into string */
  strcpy(clargs,"");
  for (i=0;i<argc;i++) {
    strcat(clargs,argv[i]);
    strcat(clargs," ");
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);

}

/** Read in list of frame files from input directory
 */
int XLALReadFrameDir(FileList **framefiles,    /**< [out] a structure containing a list of all input frame files */
		     CHAR *inputdir                 /**< [in] the input frame directory */
		     )
{
  
  const CHAR *fn = __func__;      /* store function name for log output */
  INT4 N;                         /* the number of returned files from scandir */
    
  /* check input arguments */
  if ((*framefiles) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input frame file list structure != NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }  
  if (inputdir == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input directory string == NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }  
         
  /* allocate memory for filelist structure */
  if (((*framefiles) = (FileList *)LALCalloc(1,sizeof(FileList))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for filelist structure.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* get the mode1 FITS file (apID 55) */
  N = scandir(inputdir, &((*framefiles)->list), XLALFrameFilter, alphasort);
  if (N>0) {
    INT4 i;
    for (i=0;i<N;i++) printf("%s\n",(*framefiles)->list[i]->d_name);

  }
  else if (N == -1) {
    LogPrintf(LOG_CRITICAL,"%s: scandir() failed.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
  else if (N == 0) {
    LogPrintf(LOG_CRITICAL,"%s: could not find any frame files in director %s.\n",fn,inputdir);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
   
  /* add inputdir to output structure */
  strncpy((*framefiles)->dir,inputdir,STRINGLENGTH);
  (*framefiles)->length = N;

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}


/** Filter function for reading in frame files
 */
int XLALFrameFilter(const struct dirent *x)
{
  
  /* if not *.gwf then return 0 */
  if (strstr(x->d_name,".gwf")!=NULL) return 1;
  else return 0;

}

/** Read in pcu counts from FS46 files and generate a vector of PCU numbers 
 */
int XLALReadGoodPCUInterval(GoodPCUIntervals **pcu,   /**< [out] the PCU interval information */
			    FileList *framefiles      /**< [in] the framefile list */
			    )
{
  
  const CHAR *fn = __func__;      /* store function name for log output */
  INT4 i;         
  INT4 count = 0;

  /* check input arguments */
  if (framefiles == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input frame file list structure == NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }  
  if ((*pcu) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output pcu structure != NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }  

  /* allocate memory for the output structure */
  if (((*pcu) = (GoodPCUIntervals *)LALCalloc(1,sizeof(GoodPCUIntervals))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to allocate memory for GoodPCUIntervals structure.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*pcu)->start = (LIGOTimeGPS *)LALCalloc(framefiles->length,sizeof(LIGOTimeGPS))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to allocate memory for pcu start data.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*pcu)->end = (LIGOTimeGPS *)LALCalloc(framefiles->length,sizeof(LIGOTimeGPS))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to allocate memory for pcu end data.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*pcu)->pcucount = (UINT2 *)LALCalloc(framefiles->length,sizeof(UINT2))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to allocate memory for pcu counts data.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  (*pcu)->length = framefiles->length;

  /* loop over each frame file */
  for (i=0;i<(INT4)framefiles->length;i++) {
    
    /* check if it is an FS46 file using the filename */
    if (strstr(framefiles->list[i]->d_name,"_FS46-") != NULL) {
    
      FrStream *fs = NULL;
      LIGOTimeGPS epoch;
      LIGOTimeGPS end;
      REAL8 duration;
      INT4 j;
      UINT2 pcucount = 0;
      LogPrintf(LOG_DEBUG,"%s: opening frame file %s.\n",fn,framefiles->list[i]->d_name);

      /* open the frame file */
      if ((fs = XLALFrOpen(framefiles->dir,framefiles->list[i]->d_name)) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s: unable to open FS46 frame file %s.\n",fn,framefiles->list[i]->d_name);
	XLAL_ERROR(fn,XLAL_EINVAL);
      }
      LogPrintf(LOG_DEBUG,"%s: opened frame file %s.\n",fn,framefiles->list[i]->d_name);

      /* extract start and duration of frame from stream */
      XLALGPSSetREAL8(&epoch,(REAL8)fs->flist->t0);
      duration = (REAL8)fs->flist->dt;
      LogPrintf(LOG_DEBUG,"%s: frame start time = %d.\n",fn,epoch.gpsSeconds);
      LogPrintf(LOG_DEBUG,"%s: frame duration = %f.\n",fn,duration);
    
      /* read in data for each PCU */ 
      for (j=0;j<NPCU;j++) {
	
	INT4TimeSeries *ts = NULL;
	INT4 k;
	REAL8 mean = 0;
	
	/* seek to the start of the frame */
	if (XLALFrSeek(fs,&epoch)) {
	  LogPrintf(LOG_CRITICAL,"%s: unable to seek to start of frame file %s.\n",fn,framefiles->list[i]->d_name);
	  XLAL_ERROR(fn,XLAL_EINVAL);
	}

	/* read in timeseries from this file - final arg is limit on length of timeseries (0 = no limit) */
	if ((ts = XLALFrReadINT4TimeSeries(fs,PCUchannel[j],&epoch,duration,0)) == NULL) {
	  LogPrintf(LOG_CRITICAL,"%s: unable to read channel %s from frame file %s.\n",fn,PCUchannel[j],framefiles->list[i]->d_name);
	  XLAL_ERROR(fn,XLAL_EINVAL);
	}
	LogPrintf(LOG_DEBUG,"%s: reading channel %s from %s.\n",fn,PCUchannel[j],framefiles->list[i]->d_name);
	printf("timeseries N = %d dt = %6.12f T = %f\n",ts->data->length,ts->deltaT,ts->data->length*ts->deltaT);
	
	/* determine whether each PCU was working or not */
	/* from ED MORGANS email ----------------------------------------------------------------------------------------------*/
	/* All 5 PCUs are all still functional to some extent. Any given observation may have  between 1 and 5 PCU turned on. */
	/* The average number of PCUs has been declining more or less linearly over the last 13 years. */
	/* The way I  determine the number of PCUs is to look apid 70. The first 5 columns contain the count rate per PCU per .125 s. */
	/* So I rebin to a longer time (e.g. 128s) count the number of non zero columns (of the 5). Actually even off detectors occasionally */
	/* get a count, count the number of columns with cts > a few.  */
	for (k=0;k<(INT4)ts->data->length;k++) mean += ts->data->data[k];
	mean /= (REAL8)ts->data->length;
	if (mean > PCUCOUNTTHRESH) pcucount++;

	/* free the PCU timeseries data */
	XLALDestroyINT4TimeSeries(ts);

      }

      /* fill in output structure */
      (*pcu)->start[count].gpsSeconds = epoch.gpsSeconds;
      (*pcu)->start[count].gpsNanoSeconds = epoch.gpsNanoSeconds;
      XLALGPSSetREAL8(&end,(REAL8)fs->flist->t0 + duration);
      (*pcu)->end[count].gpsSeconds = end.gpsSeconds;
      (*pcu)->end[count].gpsNanoSeconds = end.gpsNanoSeconds;
      (*pcu)->pcucount[count] = (UINT2)pcucount;
  
      printf("pcucount = %d\n",pcucount);
      /* close the frame file */
      XLALFrClose(fs);

      /* increment file count */
      count++;
      
    }

  }
  
  /* check if any FS46 files were found */
  if (count == 0) {
    LogPrintf(LOG_CRITICAL,"%s: could not find any FS46 (PCU) files.  Exiting.\n",fn);
    exit(0);
  }

  (*pcu)->length = count;

  /* resize the output vectors */
  if (((*pcu)->start = (LIGOTimeGPS *)LALRealloc((*pcu)->start,(*pcu)->length*sizeof(LIGOTimeGPS))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to resize memory for pcu start data.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*pcu)->end = (LIGOTimeGPS *)LALCalloc((*pcu)->length,sizeof(LIGOTimeGPS))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to resize memory for pcu end data.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*pcu)->pcucount = (UINT2 *)LALCalloc((*pcu)->length,sizeof(UINT2))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to resize memory for pcu counts data.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  (*pcu)->length = framefiles->length;
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}


