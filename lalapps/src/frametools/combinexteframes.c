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
/* some macros */

/** convert GPS-time to REAL8 */
#define GPS2REAL8(gps) (1.0 * (gps).gpsSeconds + 1.e-9 * (gps).gpsNanoSeconds )

#define MYMAX(x,y) ( (x) > (y) ? (x) : (y) )
#define MYMIN(x,y) ( (x) < (y) ? (x) : (y) )

/* Useful macros */
#define SECNAN_TO_I8TIME( sec, nan ) \
  ((INT8)1000000000*(INT8)(sec)+(INT8)(nan))
/* Dangerous!!!: */
#define EPOCH_TO_I8TIME( epoch ) \
  SECNAN_TO_I8TIME( (epoch).gpsSeconds, (epoch).gpsNanoSeconds )
#define SET_EPOCH( pepoch, i8time ) \
  do { INT8 t=(i8time); LIGOTimeGPS *pe=(pepoch); \
    pe->gpsSeconds=t/(INT8)1000000000; pe->gpsNanoSeconds=t%(INT8)1000000000; \
  } while( 0 )

/***********************************************************************************************/
/* define internal structures */

/** A structure that contains information regarding a specific frame file  
 */
typedef struct { 
  CHAR filename[STRINGLENGTH];      /**< the name of the file */
  LIGOTimeGPS epoch;                /**< file start time */
  LIGOTimeGPS end;                  /**< file end time */
  REAL8 duration;                   /**< file duration */
  REAL8 dt;                         /**< the sampling time */
  INT4 lld;                         /**< the lld flag */
  INT4 energy[2];                   /**< the energy channel range */
  CHAR *history;                    /**< contains the history field of the frame file */
} FrameFile;

/** A structure that stores information about a collection of Frame files
 */
typedef struct { 
  UINT4 length;	                    /**< the number of files */
  CHAR dir[STRINGLENGTH];           /**< the directory containing the files */
  FrameFile *file;                  /**< a pointert to FrameFile structures */
} FrameFileList;

/** A structure that sores user input variables 
 */
typedef struct { 
  UINT4 length;	                    /**< the number of intervals */
  LIGOTimeGPS *start;               /**< the GPS start time of an interval */
  LIGOTimeGPS *end;                 /**< the GPS end time of an interval */
  UINT2 *pcucount;                  /**< the number of operational PCUs */
} GoodPCUIntervals;

/** A structure that contains information regarding how to combine coincident frames
 */
typedef struct { 
  FrameFileList filelist;           /**< a list of files */
  REAL8 dt;                         /**< the common sampling time */
  LIGOTimeGPS epoch;                /**< combination start time */
  LIGOTimeGPS end;                  /**< combination end time */
  REAL8 duration;                   /**< combination duration */
} FrameCombinationPlan;

/** A vector of structures containing information regarding how to combine coincident frames
 */
typedef struct { 
  FrameCombinationPlan *data;       /**< a vector of plans */
  INT4 length;                      /**< the vector length */
} FrameCombinationPlanVector;

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

const char xtechannelname[16] = "X1:PHOTONCOUNTS";

/* the PCU channel names for FS46 fits files */
const char *PCUchannel[5] = {"XeCntPcu0","XeCntPcu1","XeCntPcu2","XeCntPcu3","XeCntPcu4"};

RCSID( "$Id$");		/* FIXME: use git-ID instead to set 'rcsid' */

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
void ReadUserVars(LALStatus *status,int argc,char *argv[],UserInput_t *uvar, CHAR *clargs);
int XLALReadFrameDir(FrameFileList **framefiles, CHAR *inputdir);
int XLALFrameFilter(struct dirent *x);
int XLALReadGoodPCUInterval(GoodPCUIntervals **pcu,FrameFileList *framefiles);
int XLALFindFramesInInterval(FrameFileList **subframefiles, FrameFileList *framefiles,LIGOTimeGPS intstart,LIGOTimeGPS intend);
int XLALCreateCombinationPlan(FrameCombinationPlanVector *plans,FrameFileList *framefiles);
static int compareGPS(const void *p1, const void *p2);
int XLALCombinationPlanToINT4TimeSeries(INT4TimeSeries **ts,FrameCombinationPlan *plan);
int XLALINT4TimeSeriesToFrame(CHAR *outdir,INT4TimeSeries *ts,FrameCombinationPlan *plan,CHAR *clargs);
int XLALReadFrameHistory(CHAR **history_string, FrFile *file);

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
  FrameFileList *framefiles = NULL;             /* list of input frame file names */
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

  /**********************************************************************************/
  /* FIND TIME INTERVALS CONTAINING PCU DATA */
  /**********************************************************************************/

  /* get good pcu time intervals based on FS46 files photon counts */
  if (XLALReadGoodPCUInterval(&pcu,framefiles)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadGoodPCUInterval() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }
  for (i=0;i<(INT4)pcu->length;i++) LogPrintf(LOG_DEBUG,"%s : pcu intervals -> %d %d %d\n",fn,pcu->start[i].gpsSeconds,pcu->end[i].gpsSeconds,pcu->pcucount[i]);

  /**********************************************************************************/
  /* FOR EACH INTERVAL FIND DATA AND COMBINE APPROPRIATELY */
  /**********************************************************************************/

  /* for each time interval for which we have PCU data we now find any coincident data stretches */
  for (i=0;i<(INT4)pcu->length;i++) {

    FrameFileList *subframefiles = NULL;                  /* list of frame file names within given interval */
    FrameCombinationPlanVector plans;                     /* defines how to combine frames */
    INT4 j,k;                                             /* counter */

    /* find any frames coincident with this data stretch */
    if (XLALFindFramesInInterval(&subframefiles,framefiles,pcu->start[i],pcu->end[i])) {
      LogPrintf(LOG_CRITICAL,"%s : XLALFindFramesInInterval() failed with error = %d\n",fn,xlalErrno);
      return 1;
    }

    /* if any coincident frames were found */
    if (subframefiles->length) {
      
      /* of these frames find out which ones to use and how to combine them */
      if (XLALCreateCombinationPlan(&plans,subframefiles)) {
	LogPrintf(LOG_CRITICAL,"%s : XLALCreateCombinationPlan() failed with error = %d\n",fn,xlalErrno);
	return 1;
      }

      /**********************************************************************************/
      /* FOR EACH SEGMENT OF DATA WITHIN THE INTERVAL */
      /**********************************************************************************/
      
      /* take each plan and combine the corresponding files appropriately into a timeseries */
      for (j=0;j<plans.length;j++) {
	
	INT4TimeSeries *ts = NULL;

	if (XLALCombinationPlanToINT4TimeSeries(&ts,&(plans.data[j]))) {
	  LogPrintf(LOG_CRITICAL,"%s : XLALCombinationPlanToINT4TimeSeries() failed with error = %d\n",fn,xlalErrno);
	  return 1;
	}

	/**********************************************************************************/
	/* NOW GENERATE AN OUTPUT FRAME */
	/**********************************************************************************/
	
	if (XLALINT4TimeSeriesToFrame(uvar.outputdir,ts,&(plans.data[j]),clargs)) {
	  LogPrintf(LOG_CRITICAL,"%s : XLALINT4TimeSeriesToFrame() failed with error = %d\n",fn,xlalErrno);
	  return 1;
	}

	/* free the timeseries */
	XLALDestroyINT4TimeSeries(ts);

      }
    
      /* free plan memory */
      if (plans.data) {
	for (j=0;j<plans.length;j++) {
	  for (k=0;k<(INT4)plans.data[j].filelist.length;k++) XLALFree(plans.data[j].filelist.file[k].history);
	  XLALFree(plans.data[j].filelist.file);
	}
	XLALFree(plans.data);
      }

    }

    /* free mem */
    for (j=0;j<(INT4)subframefiles->length;j++) XLALFree(subframefiles->file[j].history);
    XLALFree(subframefiles->file);
    XLALFree(subframefiles);

  }

  /**********************************************************************************/
  /* CLEAN UP */
  /**********************************************************************************/

  /* free filelist */
  for (i=0;i<(INT4)framefiles->length;i++) XLALFree(framefiles->file[i].history);
  XLALFree(framefiles->file);
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
int XLALReadFrameDir(FrameFileList **framefiles,    /**< [out] a structure containing a list of all input frame files */
		     CHAR *inputdir                 /**< [in] the input frame directory */
		     )
{
  
  const CHAR *fn = __func__;      /* store function name for log output */
  INT4 N;                         /* the number of returned files from scandir */
  INT4 i;                         /* counter */
  struct dirent **list;           /* temporary structure for storing filenames */

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
  if (((*framefiles) = (FrameFileList *)LALCalloc(1,sizeof(FrameFileList))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for filelist structure.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* get the mode1 FITS file (apID 55) */
  N = scandir(inputdir, &list, XLALFrameFilter, alphasort);
  if (N>0) {
    for (i=0;i<N;i++) printf("%s\n",list[i]->d_name);
  }
  else if (N == -1) {
    LogPrintf(LOG_CRITICAL,"%s: scandir() failed.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
  else if (N == 0) {
    LogPrintf(LOG_CRITICAL,"%s: could not find any frame files in director %s.\n",fn,inputdir);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* allocate memory for filelist metadata */
  if (((*framefiles)->file = (FrameFile *)LALCalloc(N,sizeof(FrameFile))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FrameFile structures.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
 
  /* open each file and read in some header information */
  for (i=0;i<N;i++) {

    FrStream *fs = NULL;
    INT4TimeSeries ts;
    LogPrintf(LOG_DEBUG,"%s: opening frame file %s.\n",fn,list[i]->d_name);
    
    /* fill in some timeseries params */
    snprintf(ts.name,LALNameLength,"%s0",xtechannelname);
    
    /* open the frame file */
    if ((fs = XLALFrOpen(inputdir,list[i]->d_name)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to open FS46 frame file %s.\n",fn,list[i]->d_name);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    LogPrintf(LOG_DEBUG,"%s: opened frame file %s.\n",fn,list[i]->d_name);
    
    /* try to get history information */
    if (XLALReadFrameHistory(&((*framefiles)->file[i].history),fs->file)) {
      LogPrintf(LOG_CRITICAL,"%s: XLALReadFrameHistory() unable to read history from frame file %s.\n",fn,(*framefiles)->file[i].filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    LogPrintf(LOG_DEBUG,"%s: read history field from file %s.\n",fn,list[i]->d_name);
  
    /* extract start and duration of frame from stream */
    strncpy((*framefiles)->file[i].filename,list[i]->d_name,STRINGLENGTH);
    XLALGPSSetREAL8(&((*framefiles)->file[i].epoch),(REAL8)fs->flist->t0);
    XLALGPSSetREAL8(&((*framefiles)->file[i].end),(REAL8)(fs->flist->t0 + fs->flist->dt));
    (*framefiles)->file[i].duration = (REAL8)fs->flist->dt;
    
    /* read in timeseries from this file - final arg is limit on length of timeseries (0 = no limit) */
    if (XLALFrGetINT4TimeSeriesMetadata(&ts,fs)) {
      LogPrintf(LOG_CRITICAL,"%s: unable to read from frame file %s.\n",fn,(*framefiles)->file[i].filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    (*framefiles)->file[i].dt = ts.deltaT;

    /* extract energy and LLD information from the filename */
    /* the lld flag is immediately after the 5th "_" and the min energy after the 6th and the max energy after the 7th */
     { 
       CHAR *c1;
       CHAR lldtemp[1], entemp[4];
       INT4 j;
       c1 = strstr(list[i]->d_name,"-");
       for (j=0;j<5;j++) {
	 char *c2 = strstr(c1+1,"_");
	 c1 = c2;	
       }
       strncpy(lldtemp,c1+1,1);    
       (*framefiles)->file[i].lld = atoi(lldtemp);
       c1 = strstr(c1+1,"_");     
       strncpy(entemp,c1+1,3);    
       (*framefiles)->file[i].energy[0] = atoi(entemp);    
       c1 = strstr(c1+1,"_");    
       strncpy(entemp,c1+1,3);    
       (*framefiles)->file[i].energy[1] = atoi(entemp);
       
    }

    LogPrintf(LOG_DEBUG,"%s: frame filename = %s.\n",fn,(*framefiles)->file[i].filename);
    LogPrintf(LOG_DEBUG,"%s: frame start time = %d.\n",fn,(*framefiles)->file[i].epoch.gpsSeconds);
    LogPrintf(LOG_DEBUG,"%s: frame end time = %d.\n",fn,(*framefiles)->file[i].end.gpsSeconds);
    LogPrintf(LOG_DEBUG,"%s: frame duration = %f.\n",fn,(*framefiles)->file[i].duration);
    LogPrintf(LOG_DEBUG,"%s: frame lld flag = %d.\n",fn,(*framefiles)->file[i].lld);
    LogPrintf(LOG_DEBUG,"%s: frame min energy = %d.\n",fn,(*framefiles)->file[i].energy[0]);
    LogPrintf(LOG_DEBUG,"%s: frame max energy = %d.\n",fn,(*framefiles)->file[i].energy[1]);
    LogPrintf(LOG_DEBUG,"%s: frame dt = %6.12f.\n",fn,(*framefiles)->file[i].dt);

    
    /* close the frame file */
    XLALFrClose(fs);

  }

  /* add inputdir and length to output structure */
  strncpy((*framefiles)->dir,inputdir,STRINGLENGTH);
  (*framefiles)->length = N;

  /* free original filelist */
  for (i=0;i<N;i++) free(list[i]);
  free(list);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}


/** Filter function for reading in frame files
 */
int XLALFrameFilter(struct dirent *x)
{
  
  /* if not *.gwf then return 0 */
  if (strstr(x->d_name,".gwf")!=NULL) return 1;
  else return 0;

}

/** Read in pcu counts from FS46 files and generate a vector of PCU numbers 
 */
int XLALReadGoodPCUInterval(GoodPCUIntervals **pcu,   /**< [out] the PCU interval information */
			    FrameFileList *framefiles      /**< [in] the framefile list */
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
    if (strstr(framefiles->file[i].filename,"_FS46_") != NULL) {
    
      FrStream *fs = NULL;
      LIGOTimeGPS epoch;
      LIGOTimeGPS end;
      REAL8 duration;
      INT4 j;
      UINT2 pcucount = 0;
      LogPrintf(LOG_DEBUG,"%s: opening frame file %s.\n",fn,framefiles->file[i].filename);

      /* open the frame file */
      if ((fs = XLALFrOpen(framefiles->dir,framefiles->file[i].filename)) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s: unable to open FS46 frame file %s.\n",fn,framefiles->file[i].filename);
	XLAL_ERROR(fn,XLAL_EINVAL);
      }
      LogPrintf(LOG_DEBUG,"%s: opened frame file %s.\n",fn,framefiles->file[i].filename);

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
	CHAR channelname[STRINGLENGTH];

	/* fill in some timeseries params */
	snprintf(channelname,STRINGLENGTH,"%s%d",xtechannelname,j);
	
	/* seek to the start of the frame */
	if (XLALFrSeek(fs,&epoch)) {
	  LogPrintf(LOG_CRITICAL,"%s: unable to seek to start of frame file %s.\n",fn,framefiles->file[i].filename);
	  XLAL_ERROR(fn,XLAL_EINVAL);
	}

	/* read in timeseries from this file - final arg is limit on length of timeseries (0 = no limit) */
	if ((ts = XLALFrReadINT4TimeSeries(fs,channelname,&epoch,duration,0)) == NULL) {
	  LogPrintf(LOG_CRITICAL,"%s: unable to read channel %s from frame file %s.\n",fn,PCUchannel[j],framefiles->file[i].filename);
	  XLAL_ERROR(fn,XLAL_EINVAL);
	}
	LogPrintf(LOG_DEBUG,"%s: reading channel %s from %s.\n",fn,channelname,framefiles->file[i].filename);
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

      /* fill in output structure if there are any operational PCUs */
      if (pcucount) {
	(*pcu)->start[count].gpsSeconds = epoch.gpsSeconds;
	(*pcu)->start[count].gpsNanoSeconds = epoch.gpsNanoSeconds;
	XLALGPSSetREAL8(&end,(REAL8)fs->flist->t0 + duration);
	(*pcu)->end[count].gpsSeconds = end.gpsSeconds;
	(*pcu)->end[count].gpsNanoSeconds = end.gpsNanoSeconds;
	(*pcu)->pcucount[count] = (UINT2)pcucount;
	
	/* increment file count */
	count++;
      }
      printf("pcucount = %d\n",pcucount);
      /* close the frame file */
      XLALFrClose(fs);
            
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
  if (((*pcu)->end = (LIGOTimeGPS *)LALRealloc((*pcu)->end,(*pcu)->length*sizeof(LIGOTimeGPS))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to resize memory for pcu end data.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*pcu)->pcucount = (UINT2 *)LALRealloc((*pcu)->pcucount,(*pcu)->length*sizeof(UINT2))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to resize memory for pcu counts data.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}

/** Finds a subset of frame files within a given time interval 
 */
int XLALFindFramesInInterval(FrameFileList **subframefiles,   /**< [out] a list of filenames containing data within the interval */
			     FrameFileList *framefiles,       /**< [in] the framefile list */
			     LIGOTimeGPS intstart,            /**< [in] the interval start time */
			     LIGOTimeGPS intend               /**< [in] the interval end time */
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
  if ((*subframefiles) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output FileList structure != NULL.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (XLALGPSCmp(&intstart,&intend) >= 0) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, the interval start time is >= interval end time.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
  
  /* allocate memory for output filelist structure */
  if (((*subframefiles) = (FrameFileList *)LALCalloc(1,sizeof(FrameFileList))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for filelist structure.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* allocate memory for filelist metadata */
  if (((*subframefiles)->file = (FrameFile *)LALCalloc(framefiles->length,sizeof(FrameFile))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FrameFile structures.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* loop over the input frame file names */
  for (i=0;i<(INT4)framefiles->length;i++) {
    
    /* check if it is NOT an FS46 file using the filename */
    if (strstr(framefiles->file[i].filename,"_FS46_") == NULL) {
     
      LIGOTimeGPS *framestart = &(framefiles->file[i].epoch);
      LIGOTimeGPS *frameend = &(framefiles->file[i].end);;     

      /* check if file overlaps with interval */
      if ( ! (( XLALGPSCmp(framestart,&intend) > 0 ) || ( XLALGPSCmp(frameend,&intstart) < 0 )) ) {
	
	/* compute actual overlap */
	REAL8 x = MYMIN(intstart.gpsSeconds,framestart->gpsSeconds);
	REAL8 y = MYMAX(intend.gpsSeconds,frameend->gpsSeconds);
	REAL8 overlap = y-x;
	
	/* record overlapping frames */
	strncpy((*subframefiles)->file[count].filename,framefiles->file[i].filename,STRINGLENGTH);
	memcpy(&((*subframefiles)->file[count].epoch),&(framefiles->file[i].epoch),sizeof(LIGOTimeGPS));
 	memcpy(&((*subframefiles)->file[count].end),&(framefiles->file[i].end),sizeof(LIGOTimeGPS));
	(*subframefiles)->file[count].duration = framefiles->file[i].duration;
	(*subframefiles)->file[count].lld = framefiles->file[i].lld;
	(*subframefiles)->file[count].energy[0] = framefiles->file[i].energy[0];
	(*subframefiles)->file[count].energy[1] = framefiles->file[i].energy[1];
	(*subframefiles)->file[count].dt = framefiles->file[i].dt;

	/* copy history information */
	if (((*subframefiles)->file[count].history = (CHAR *)XLALCalloc(strlen(framefiles->file[i].history)+1,sizeof(CHAR))) == NULL ) {
	  LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for history fields.\n",fn);
	  XLAL_ERROR(fn,XLAL_ENOMEM);
	}
	strncpy((*subframefiles)->file[count].history,framefiles->file[i].history,strlen(framefiles->file[i].history));

	LogPrintf(LOG_DEBUG,"%s : Int [%d->%d] : overlapping frame [%d->%d] = %.0f overlap.\n",fn,intstart.gpsSeconds,intend.gpsSeconds,framestart->gpsSeconds,frameend->gpsSeconds,overlap);
	count++;
      }

    }

  }

  /* fill in additional information */
  (*subframefiles)->length = count;
  strncpy((*subframefiles)->dir,framefiles->dir,STRINGLENGTH);
  
  /* resize the output vectors */
  if ((*subframefiles)->length) {
    if (((*subframefiles)->file = (FrameFile *)LALRealloc((*subframefiles)->file,(*subframefiles)->length*sizeof(FrameFile))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: failed to resize memory for subframefile list.\n",fn);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
  }
  else {
    XLALFree((*subframefiles)->file);
    (*subframefiles)->file = NULL;
  }
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}

/** Finds a subset of frame files within a given time interval 
 */
int XLALCreateCombinationPlan(FrameCombinationPlanVector *plans,     /**< [out] a plan of how to combine the frames */
			      FrameFileList *framefiles               /**< [in] the framefile list */
			      )
{
  
  const CHAR *fn = __func__;      /* store function name for log output */ 
  LIGOTimeGPS *lateststart = &(framefiles->file[0].epoch);
  LIGOTimeGPS *earliestend = &(framefiles->file[0].end);
  LIGOTimeGPS *earlieststart = &(framefiles->file[0].epoch);
  LIGOTimeGPS *latestend = &(framefiles->file[0].end);
  INT4 i,s;                         /* counter */
  LIGOTimeGPSVector *epochlist;
  INT4 k = 0;
  REAL8 minspan = 0.0;
  REAL8 maxspan = 0.0;
  INT4 N = 100;
  INT4 pcount = 0;

  /* loop over each file */
  for (i=0;i<(INT4)framefiles->length;i++) {
    
    if (XLALGPSCmp(&(framefiles->file[i].epoch),earlieststart) < 0) earlieststart = &(framefiles->file[i].epoch);
    if (XLALGPSCmp(&(framefiles->file[i].epoch),lateststart) > 0) lateststart = &(framefiles->file[i].epoch);
    if (XLALGPSCmp(&(framefiles->file[i].end),earliestend) < 0) earliestend = &(framefiles->file[i].end);
    if (XLALGPSCmp(&(framefiles->file[i].end),latestend) > 0) latestend = &(framefiles->file[i].end);

  }
  minspan = XLALGPSDiff(earliestend,lateststart);
  maxspan = XLALGPSDiff(latestend,earlieststart);
  for (i=0;i<(INT4)framefiles->length;i++) fprintf(stdout,"%s\n",framefiles->file[i].filename);
  
  /* allocate nmemory for the timestamps list */
  epochlist = XLALCreateTimestampVector(2*framefiles->length);

  /* loop over each file */
  for (i=0;i<(INT4)framefiles->length;i++) {
    
    /* draw (in ascii) the data stretches */
    {
      INT4 sidx = floor(0.5 + N*(XLALGPSGetREAL8(&(framefiles->file[i].epoch)) - XLALGPSGetREAL8(earlieststart))/maxspan);
      INT4 eidx = floor(0.5 + N*(XLALGPSGetREAL8(&(framefiles->file[i].end)) - XLALGPSGetREAL8(earlieststart))/maxspan);
      INT4 j;

      for (j=0;j<sidx;j++) fprintf(stdout," ");
      fprintf(stdout,"|");
      for (j=sidx+1;j<eidx-1;j++) fprintf(stdout,"-");
      fprintf(stdout,"| dt = %6.6f energy = %d-%d lld = %d\n",framefiles->file[i].dt,framefiles->file[i].energy[0],framefiles->file[i].energy[1],framefiles->file[i].lld);
    }
    
    /* add start and end to a list of times that will then be sorted */
    epochlist->data[k].gpsSeconds = framefiles->file[i].epoch.gpsSeconds;
    epochlist->data[k].gpsNanoSeconds = framefiles->file[i].epoch.gpsNanoSeconds;
    epochlist->data[k+1].gpsSeconds = framefiles->file[i].end.gpsSeconds;
    epochlist->data[k+1].gpsNanoSeconds = framefiles->file[i].end.gpsNanoSeconds;
    k = k + 2;

  }
       
  /* sort the epochs */
  qsort(epochlist->data, epochlist->length, sizeof(LIGOTimeGPS), compareGPS);
  for (i=0;i<(INT4)epochlist->length;i++) LogPrintf(LOG_DEBUG,"%s : epoch boundaries %d %d\n",fn,epochlist->data[i].gpsSeconds,epochlist->data[i].gpsNanoSeconds);

  /* allocate memory for the current plan under the assumption that there will be epochlist->length-1 seperate plans */
  plans->length = epochlist->length-1;
  if ((plans->data = (FrameCombinationPlan *)LALCalloc(plans->length,sizeof(FrameCombinationPlan))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FrameCombinationPlan.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  
  /* loop over each epoch span and determine which files belong within it */
  for (i=0;i<(INT4)epochlist->length-1;i++) {

    FrameFileList tempfilelist;
    INT4Vector *badidx = NULL;
    INT4 nfiles = 0;
    INT4 count = 0;

    LogPrintf(LOG_DEBUG,"%s : working on span #%d : %d -> %d\n",fn,i,epochlist->data[i].gpsSeconds,epochlist->data[i+1].gpsSeconds);
    
    /* count number of coincident files */
    for (k=0;k<(INT4)framefiles->length;k++) {
      if ( !((XLALGPSCmp(&(framefiles->file[k].epoch),&(epochlist->data[i+1])) >= 0) || (XLALGPSCmp(&(framefiles->file[k].end),&(epochlist->data[i])) <= 0)) ) {
	LogPrintf(LOG_DEBUG,"%s : coincident file = %s\n",fn,framefiles->file[k].filename);
	nfiles++;
      }
    }
    
    /* allocate memory for the temporary filelist */
    tempfilelist.length = nfiles;
    snprintf(tempfilelist.dir,STRINGLENGTH,"%s",framefiles->dir);
    if ((tempfilelist.file = (FrameFile *)LALCalloc(tempfilelist.length,sizeof(FrameFile))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FrameFile structures.\n",fn);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
    
    /* allocate memory for the list of bad indices and initialise */
    badidx = XLALCreateINT4Vector(nfiles);
    for (k=0;k<nfiles;k++) badidx->data[k] = 0;

    /* find the coincident files again and record the filenames */
    for (k=0;k<(INT4)framefiles->length;k++) {
      if ( !((XLALGPSCmp(&(framefiles->file[k].epoch),&(epochlist->data[i+1])) >= 0) || (XLALGPSCmp(&(framefiles->file[k].end),&(epochlist->data[i])) <= 0)) ) {
	      
	/* fill in temporary filelist */
	strncpy(tempfilelist.file[count].filename,framefiles->file[k].filename,STRINGLENGTH);
	memcpy(&(tempfilelist.file[count].epoch),&(framefiles->file[k].epoch),sizeof(LIGOTimeGPS));
	memcpy(&(tempfilelist.file[count].end),&(framefiles->file[k].end),sizeof(LIGOTimeGPS));
	tempfilelist.file[count].dt = framefiles->file[k].dt;
	tempfilelist.file[count].duration = framefiles->file[k].duration;
	tempfilelist.file[count].energy[0] = framefiles->file[k].energy[0];
	tempfilelist.file[count].energy[1] = framefiles->file[k].energy[1];
	tempfilelist.file[count].lld = framefiles->file[k].lld;
	
	/* copy history information */
	if ((tempfilelist.file[count].history = (CHAR *)XLALCalloc(strlen(framefiles->file[k].history)+1,sizeof(CHAR))) == NULL ) {
	  LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for history fields.\n",fn);
	  XLAL_ERROR(fn,XLAL_ENOMEM);
	}
	strncpy(tempfilelist.file[count].history,framefiles->file[k].history,strlen(framefiles->file[k].history));

	count++;

      }
    }
    LogPrintf(LOG_DEBUG,"%s : recorded the coincident filenames\n",fn);

    /* now check for single lld files */
    if ( (nfiles == 1) && (tempfilelist.file[0].lld == 1) ) badidx->data[0] = 1;
    LogPrintf(LOG_DEBUG,"%s : checked for single lld files\n",fn);

    /* now check for energy consistency */
    if (tempfilelist.file) {
      
      for (k=0;k<(INT4)tempfilelist.length;k++) {
	for (s=0;s<(INT4)k;s++) {
	  /* check for energy overlap if not ldd */
	  if ( ( tempfilelist.file[k].lld == 0 ) && ( tempfilelist.file[s].lld == 0 ) ) {
	    if ( ( ( tempfilelist.file[k].energy[0] < tempfilelist.file[s].energy[0] ) && ( tempfilelist.file[k].energy[1] > tempfilelist.file[s].energy[0] ) )
		 || ( ( tempfilelist.file[k].energy[1] > tempfilelist.file[s].energy[1] ) && ( tempfilelist.file[k].energy[0] < tempfilelist.file[s].energy[1] ) ) ) {
	      
	      /* determine which file to keep based on largest energy range */
	      INT4 range1 = tempfilelist.file[k].energy[1] - tempfilelist.file[k].energy[0];
	      INT4 range2 = tempfilelist.file[s].energy[1] - tempfilelist.file[s].energy[0];
	      INT4 highenergyidx =  tempfilelist.file[k].energy[1] > tempfilelist.file[s].energy[1] ? k : s;
	      if ( range1 > range2 ) badidx->data[s] = 1;
	      else if ( range2 > range1 ) badidx->data[k] = 1;
	      else badidx->data[highenergyidx] = 1;
	      
	    }
	  }
	} 
      }
    }
    LogPrintf(LOG_DEBUG,"%s : checked for energy consistency\n",fn);

    /* now fill in the plan by picking the non-bad indices */
    {
      INT4 sum = badidx->length;
      INT4 newcount = 0;
      REAL8 largestdt = 0.0;
      for (k=0;k<(INT4)badidx->length;k++) sum -= badidx->data[k];
    
      /* allocate memory for the current plan if thee was any good data */
      if (sum) {

	if ((plans->data[pcount].filelist.file = (FrameFile *)LALCalloc(sum,sizeof(FrameFile))) == NULL) {
	  LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FrameCombinationPlan.\n",fn);
	  XLAL_ERROR(fn,XLAL_ENOMEM);
	}
	
	/* loop over coincident files and pick out the good ones */
	for (k=0;k<(INT4)tempfilelist.length;k++) {
	  
	  if (badidx->data[k] == 0) {
	    
	    /* find largest sampling time */
	    if (tempfilelist.file[k].dt > largestdt) largestdt = tempfilelist.file[k].dt;

	    strncpy(plans->data[pcount].filelist.file[newcount].filename,tempfilelist.file[k].filename,STRINGLENGTH);
	    memcpy(&(plans->data[pcount].filelist.file[newcount].epoch),&(tempfilelist.file[k].epoch),sizeof(LIGOTimeGPS));
	    memcpy(&(plans->data[pcount].filelist.file[newcount].end),&(tempfilelist.file[k].end),sizeof(LIGOTimeGPS));
	    plans->data[pcount].filelist.file[newcount].dt = tempfilelist.file[k].dt;
	    plans->data[pcount].filelist.file[newcount].duration = tempfilelist.file[k].duration;
	    plans->data[pcount].filelist.file[newcount].energy[0] = tempfilelist.file[k].energy[0];
	    plans->data[pcount].filelist.file[newcount].energy[1] = tempfilelist.file[k].energy[1];
	    plans->data[pcount].filelist.file[newcount].lld = tempfilelist.file[k].lld;

	    /* copy history information */
	    if ((plans->data[pcount].filelist.file[newcount].history = (CHAR *)XLALCalloc(strlen(tempfilelist.file[k].history)+1,sizeof(CHAR))) == NULL ) {
	      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for history fields.\n",fn);
	      XLAL_ERROR(fn,XLAL_ENOMEM);
	    }
	    strncpy(plans->data[pcount].filelist.file[newcount].history,tempfilelist.file[k].history,strlen(tempfilelist.file[k].history));
	    
	    newcount++;
	    
	  }
	  
	}

	/* fill in extra plan params */
	plans->data[pcount].dt = largestdt;
	memcpy(&(plans->data[pcount].epoch),&(epochlist->data[i]),sizeof(LIGOTimeGPS));
	memcpy(&(plans->data[pcount].end),&(epochlist->data[i+1]),sizeof(LIGOTimeGPS));
	plans->data[pcount].duration = XLALGPSDiff(&(epochlist->data[i+1]),&(epochlist->data[i]));
	plans->data[pcount].filelist.length = newcount;
	snprintf(plans->data[pcount].filelist.dir,STRINGLENGTH,"%s",tempfilelist.dir);
	pcount++;
      }
    
    }
    LogPrintf(LOG_DEBUG,"%s : recorded plan.\n",fn);

    /* free the temp filelist */
    for (k=0;k<(INT4)tempfilelist.length;k++) XLALFree(tempfilelist.file[k].history);
    XLALFree(tempfilelist.file);
    XLALDestroyINT4Vector(badidx);

  }

  /* resize plan */
  plans->length = pcount;
  if ((plans->data = (FrameCombinationPlan *)LALRealloc(plans->data,plans->length*sizeof(FrameCombinationPlan))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FrameCombinationPlan.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* free mem */
  XLALDestroyTimestampVector(epochlist);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}

/** wrapper for XLALGPSCmp for use with qsort
 */
static int compareGPS(const void *p1, const void *p2)
{
  /* The actual arguments to this function are "pointers to
     pointers to char", but strcmp(3) arguments are "pointers
     to char", hence the following cast plus dereference */
  
  return XLALGPSCmp(p1,p2);
  
}

/** this function combines the files listed in the combination plan into a single timeseries 
 */
int XLALCombinationPlanToINT4TimeSeries(INT4TimeSeries **ts,           /**< [out] the timeseries containing the combined data */
					FrameCombinationPlan *plan     /**< [in] the plan describing which files to combine */
					)
{

  const CHAR *fn = __func__;      /* store function name for log output */ 
  INT8 N;
  INT4 i,j,k;

  /* determine how many samples are required for the timeseries */
  N = floor(0.5 + plan->duration/plan->dt);
  LogPrintf(LOG_DEBUG,"%s : combined timeseries requires %ld samples.\n",fn,N);
  
  /* create the output timeseries */
  if ( ((*ts) = XLALCreateINT4TimeSeries("X1:COMBINED_PHOTONCOUNTS",&(plan->epoch),0,plan->dt,&lalDimensionlessUnit,N)) == NULL) {
    LogPrintf(LOG_CRITICAL, "%s : XLALCreateINT4TimeSeries() failed to allocate an %d length timeseries with error = %d.\n",fn,N,xlalErrno);
    XLAL_ERROR(fn,XLAL_ENOMEM);  
  }
  LogPrintf(LOG_DEBUG,"%s : allocated memory for temporary timeseries\n",fn);

  /* initialise the output timeseries */
  memset((*ts)->data->data,0,(*ts)->data->length*sizeof(INT4));

  /* loop over each frame file in the plan */
  for (i=0;i<(INT4)plan->filelist.length;i++) {

    INT4TimeSeries *tempts = NULL;
    FrStream *fs = NULL;
    CHAR channelname[STRINGLENGTH];
    INT4 lldfactor = 1;
 
    /* define channel name */
    snprintf(channelname,STRINGLENGTH,"%s0",xtechannelname);
    
    /* open the frame file */
    if ((fs = XLALFrOpen(plan->filelist.dir,plan->filelist.file[i].filename)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to open frame file %s.\n",fn,plan->filelist.file[i].filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    LogPrintf(LOG_DEBUG,"%s: opened frame file %s.\n",fn,plan->filelist.file[i].filename);
    
    /* seek to the start of the frame */
    if (XLALFrSeek(fs,&(plan->epoch))) {
      LogPrintf(LOG_CRITICAL,"%s: unable to seek to start of frame file %s.\n",fn,plan->filelist.file[i].filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }

    /* read in timeseries from this file - final arg is limit on length of timeseries (0 = no limit) */
    if ((tempts = XLALFrReadINT4TimeSeries(fs,channelname,&(plan->epoch),plan->duration,0)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to read channel %s from frame file %s.\n",fn,plan->filelist.file[i].filename);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
    LogPrintf(LOG_DEBUG,"%s: reading channel %s from %s.\n",fn,channelname,plan->filelist.file[i].filename);
    
    /* if the data is lld data then we multiply by 2 before adding */
    if (plan->filelist.file[i].lld) lldfactor = 2;

    /* compute the rebinning factor */
    {  
      REAL8 realbinratio = (REAL8)((*ts)->deltaT/tempts->deltaT);
      INT4 intbinratio = (INT4)((*ts)->deltaT/tempts->deltaT);
      
      if ( fmod(realbinratio,1.0) > 0.0 ) {
	LogPrintf(LOG_CRITICAL,"%s: sampling rate is not consistent with an integer factor in file %s.\n",fn,plan->filelist.file[i].filename);
	XLAL_ERROR(fn,XLAL_EINVAL);
      }
      
      /* add to output timeseries if the timeseries are consistent */
      if ( XLALGPSCmp(&((*ts)->epoch),&(tempts->epoch)) == 0 ) {
	
	INT4 idx = 0;
	
	/* loop over each output bin */
	for (j=0;j<(INT4)(*ts)->data->length;j++) {
	 	  
	  /* perform rebinning on the fly */
	  INT4 temp = 0;
	  for (k=0;k<intbinratio;k++) {
	    temp += tempts->data->data[idx];
	    idx++;
	  }

	  /* add to output timeseries */
	  (*ts)->data->data[j] += lldfactor*temp;
	}
      
      }
      
    }
      
    /* free memory */
    XLALDestroyINT4TimeSeries(tempts);

    /* close the frame file */
    XLALFrClose(fs);
 
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}

/** this function combines the files listed in the combination plan into a single timeseries 
 */
int XLALINT4TimeSeriesToFrame(CHAR *outputdir,               /**< [in] name of output directory */
			      INT4TimeSeries *ts,            /**< [in] timeseries to output */
			      FrameCombinationPlan *plan,    /**< [in] the plan used to generate the timeseries */
			      CHAR *clargs                   /**< [in] the command line args */
			      )
{
  
  const CHAR *fn = __func__;      /* store function name for log output */ 
  struct FrameH *outFrame = NULL;        /* frame data structure */
  REAL8 T;
  INT4 i;
  CHAR filecomment[STRINGLENGTH];
  CHAR outputfile[STRINGLENGTH];

  /* define observation time */
  T = ts->deltaT*ts->data->length;

  /* generate a frame data structure - last threee inputs are [project, run, frnum, detectorFlags] */
  if ((outFrame = XLALFrameNew(&(ts->epoch),T,"XTE_PCA",1,0,0)) == NULL) {
    LogPrintf(LOG_CRITICAL, "%s : XLALFrameNew() failed with error = %d.\n",fn,xlalErrno);
    XLAL_ERROR(fn,XLAL_EFAILED);
  }
  LogPrintf(LOG_DEBUG,"%s : set-up frame structure\n",fn);

  /* add timeseries to frame structure */
  if (XLALFrameAddINT4TimeSeriesProcData(outFrame,ts)) {
    LogPrintf(LOG_CRITICAL, "%s : XLALFrameAddINT4TimeSeries() failed with error = %d.\n",fn,xlalErrno);
    XLAL_ERROR(fn,XLAL_EFAILED);
  }
  LogPrintf(LOG_DEBUG,"%s : added timeseries to frame structure\n",fn);

  /* Here's where we add extra information into the frame */
  /* we include the current command line args and the original FITS file headers from each contributing file */ 
  /* we also add the git version info */
  {
    CHAR *versionstring = NULL;              /* pointer to a string containing the git version information */ 
    versionstring = XLALGetVersionString(1);
    FrHistoryAdd(outFrame,versionstring);
    FrHistoryAdd(outFrame,clargs);
    for (i=0;i<(INT4)plan->filelist.length;i++) {
      printf("outputting history #%d\n%s\n",i,plan->filelist.file[i].history);
      FrHistoryAdd(outFrame,plan->filelist.file[i].history);
    }
    XLALFree(versionstring);
  }
  
  /* construct file name - we use the LIGO format <DETECTOR>-<COMMENT>-<GPSSTART>-<DURATION>.gwf */
  /* the comment field we sub-format into <INSTRUMENT>_<FRAME>_<SOURCE>_<OBSID_APID> */
  /* first we need to extract parts of the original filenames */
  {
    CHAR originalfile[STRINGLENGTH];
    CHAR *c1,*c2;
    INT4 j;
    INT4 n;
    snprintf(originalfile,STRINGLENGTH,"%s",plan->filelist.file[0].filename);
    c1 = strstr(originalfile,"-");
    c2 = c1;
    for (j=0;j<4;j++) {
      CHAR *temp = strstr(c2+1,"_");
      c2 = temp;	
    }
    n = strlen(c1) - strlen(c2) - 1;
    snprintf(filecomment,n,"%s",c1+1);
  }

  /* generate output filename */
  snprintf(outputfile,STRINGLENGTH,"%s/X1-%s-%d-%d.gwf",outputdir,filecomment,ts->epoch.gpsSeconds,(INT4)T);
  LogPrintf(LOG_DEBUG,"%s : output file = %s\n",fn,outputfile);

  /* write frame structure to file (opens, writes, and closes file) - last argument is compression level */
  if (XLALFrameWrite(outFrame,outputfile,1)) {
    LogPrintf(LOG_CRITICAL, "%s : XLALFrameWrite() failed with error = %d.\n",fn,xlalErrno);
    XLAL_ERROR(fn,XLAL_EFAILED);
  }
 
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}

/** this function reads in the frame history as a string
 */
int XLALReadFrameHistory(CHAR **history_string,     /**< [out] the history field read in as a string */ 
			 FrFile *file              /**< [in] frame file pointer */
			 )
{

  const CHAR *fn = __func__;      /* store function name for log output */ 
  FrameH *frame = NULL;
  INT4 stringlen = 1;
  FrHistory *localhist;

  /* check input */
  if ((*history_string) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input history string is not null.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (file == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input frame file pointer is null.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* read the frame */
  frame = FrameRead(file);

  /* initialise the string to start with */
  (*history_string) = (CHAR *)XLALCalloc(stringlen,sizeof(CHAR));
  
  localhist = frame->history;
  while (localhist) {
    
    /* get length of history string */
    INT4 n = strlen(localhist->comment);
    stringlen += n + 1;
    
    /* extend the length of the output to include the current string */
    if ( ( (*history_string) = (CHAR *)XLALRealloc((*history_string),stringlen*sizeof(CHAR))) == NULL ) {
      LogPrintf(LOG_CRITICAL,"%s : failed to re-allocate memory for history string.\n",fn);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
    
    /* append the current history string to the output */
    strncat((*history_string),localhist->comment,n);
  
    /* point to the next history record */
    localhist = localhist->next;
  }
 
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;

}
