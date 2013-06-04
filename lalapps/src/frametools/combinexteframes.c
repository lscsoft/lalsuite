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
#include <glob.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/LALFrameIO.h>
#include <lal/SFTutils.h>
#include <lal/LALFrStream.h>
#include <lalappsfrutils.h>
#include <lalapps.h>

/* remove me */
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALFrameIO.h>
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/LALDatatypes.h>

/***********************************************************************************************/
/* some global constants */
#define MAX_DT 0.000976562            /* the maximum sample time we will use (equiv max freq 512 Hz) */
#define MIN_DT 0.000244140625         /* the minimum sample time we will use (equiv max freq 2048 Hz) */
#define APIDLENGTH 5                  /* the length of APID HEX strings */
#define STRINGLENGTH 256              /* the length of general string */
#define APIDLENGTH 5                  /* the length of an APID string */
#define LONGSTRINGLENGTH 1024         /* the length of general string */
#define NAPID 6                       /* the number of valid APIDs we can currently read */
#define MINFRAMELENGTH 10             /* the minimum duration of output frame file in seconds */
#define NPCU 5                        /* the number of PCUs on XTE */
#define PCUCOUNTTHRESH 5              /* threshold on the counts per PCU that determine operational status */
#define ARRAY 0                       /* data type codes */
#define EVENT 1                       /* data type codes */
#define PCU_AREA 0.13                 /* the collecting area of a single PCU in square metres */

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

typedef struct { 
  CHAR *header_string;        /**< the ASCII header string */
  INT4 length;                /**< the length of the string */
} Header;

typedef struct { 
  Header *data;               /**<  avector of ASCII headers */
  INT4 length;                /**< the length of the vector */
} HeaderVector;

typedef struct { 
  CHAR filename[STRINGLENGTH];      /**< the name of the file */
  CHAR channelname[STRINGLENGTH];   /**< name of the channel from which the data was read */
  REAL8 dt;                         /**< the time step size in seconds */
  LIGOTimeGPS epoch;                /**< file start time */
  LIGOTimeGPS end;                  /**< file end time */
  REAL8 duration;                   /**< file duration */
  INT4 lld;                         /**< the lld flag */
  INT4 type;                        /**< event or array data */
  INT4 energy[2];                   /**< the energy channel range (0-255) */
  INT4 detconfig[5];                /**< contains detector config flags */
  CHAR OBS_ID[STRINGLENGTH];        /**< the OBS_ID of the interval */      
} FrameChannel;

/** A structure that stores information about a collection of Frame files
 */
typedef struct { 
  UINT4 length;	                    /**< the number of channels */
  CHAR dir[STRINGLENGTH];           /**< the directory containing the files */
  INT4 npcus;                       /**< the number of operational PCUs */ 
  FrameChannel *channel;            /**< a pointer to FrameChannel structures */
} FrameChannelList;

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
  FrameChannelList channellist;     /**< a list of frame channels */
  CHAR OBS_ID[STRINGLENGTH];        /**< the OBS_ID of this observation */
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
  CHAR *pattern;                    /**< the input frame file name pattern */
  CHAR *outputdir;                  /**< name of output directory */
  REAL8 deltat;                     /**< the desired sampling time */
  INT4 goodxenon;                   /**< flag for using goodxenon data */
  BOOLEAN version;	            /**< output version-info */
} UserInput_t;

/***********************************************************************************************/
/* global variables */
extern int vrbflg;	 	/**< defined in lalapps.c */

/* fill in the HEX APID names */
const char *APID[6] = {"FS37","FS3b","FS3f","FS4f","XENO","FS46"};

const char xtechannelname[16] = "X1";
/* char *glob_pattern;    */         /* used to store global frame file name pattern for scandir */

/* the PCU channel names for FS46 fits files */
const char *PCUchannel[5] = {"XeCntPcu0","XeCntPcu1","XeCntPcu2","XeCntPcu3","XeCntPcu4"};


/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
void ReadUserVars(LALStatus *status,int argc,char *argv[],UserInput_t *uvar, CHAR *clargs);
int XLALReadFrameDir(FrameChannelList **framechannels, CHAR *inputdir, CHAR *pattern,INT4 goodxenon);
int XLALReadGoodPCUInterval(GoodPCUIntervals **pcu,FrameChannelList *framechannels);
int XLALFindFramesInInterval(FrameChannelList **subframechannels, FrameChannelList *framechannels,LIGOTimeGPS intstart,LIGOTimeGPS intend,UINT2 npcus);
int XLALCreateCombinationPlan(FrameCombinationPlanVector *plans,FrameChannelList *framechannels,LIGOTimeGPS *intstart,LIGOTimeGPS *intend);
static int compareGPS(const void *p1, const void *p2);
int XLALCombinationPlanToREAL4TimeSeries(REAL4TimeSeries **ts,HeaderVector *history,FrameCombinationPlan *plan,UINT2 npcus);
int XLALREAL4TimeSeriesToFrame(CHAR *outdir,REAL4TimeSeries *ts,FrameCombinationPlan *plan,CHAR *clargs,HeaderVector *history);
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
int main( int argc, char *argv[] )
{
  LALStatus status = blank_status;              /* empty LAL status structure */
  UserInput_t uvar = empty_UserInput;           /* user input variables */
  FrameChannelList *framechannels = NULL;       /* list of input frame channels */
  GoodPCUIntervals *pcu = NULL;                 /* the operational pcu count */ 
  CHAR clargs[LONGSTRINGLENGTH];                /* store the command line args */ 
  UINT4 i;                                       /* counter */

  vrbflg = 1;	                        /* verbose error-messages */

  /* turn off default GSL error handler */
  gsl_set_error_handler_off();

  /* register and read all user-variables */
  LAL_CALL (ReadUserVars(&status,argc,argv,&uvar,clargs), &status);
  LogPrintf(LOG_DEBUG,"%s : read in uservars\n",__func__);
 
  /**********************************************************************************/
  /* READ FILES */
  /**********************************************************************************/

  /* get a list of frame file names */
  if (XLALReadFrameDir(&framechannels,uvar.inputdir,uvar.pattern,uvar.goodxenon)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadFrameDir() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }

  /**********************************************************************************/
  /* FIND TIME INTERVALS CONTAINING PCU DATA */
  /**********************************************************************************/

  /* get good pcu time intervals based on FS46 files photon counts */
  if (XLALReadGoodPCUInterval(&pcu,framechannels)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadGoodPCUInterval() failed with error = %d\n",__func__,xlalErrno);
    return 1;
  }
  for (i=0;i<pcu->length;i++) LogPrintf(LOG_DEBUG,"%s : pcu intervals -> %d %d (%d) #PCU = %d\n",__func__,
					      pcu->start[i].gpsSeconds,pcu->end[i].gpsSeconds,
					      pcu->end[i].gpsSeconds-pcu->start[i].gpsSeconds,
					      pcu->pcucount[i]);
 
  /**********************************************************************************/
  /* FOR EACH INTERVAL FIND DATA AND COMBINE APPROPRIATELY */
  /**********************************************************************************/

  /* for each time interval for which we have PCU data we now find any coincident data stretches */
  for (i=0;i<pcu->length;i++) {

    FrameChannelList *subframechannels = NULL;            /* list of frame channel names within given interval */
    FrameCombinationPlanVector plans;                     /* defines how to combine frames */
    INT4 j;                                               /* counter */

    /* find any frames coincident with this data stretch */
    if (XLALFindFramesInInterval(&subframechannels,framechannels,pcu->start[i],pcu->end[i],pcu->pcucount[i])) {
      LogPrintf(LOG_CRITICAL,"%s : XLALFindFramesInInterval() failed with error = %d\n",__func__,xlalErrno);
      return 1;
    }
  
    /* if any coincident frames were found */
    if (subframechannels->length) {
      
      /* of these frames find out which ones to use and how to combine them */
      if (XLALCreateCombinationPlan(&plans,subframechannels,&(pcu->start[i]),&(pcu->end[i]))) {
	LogPrintf(LOG_CRITICAL,"%s : XLALCreateCombinationPlan() failed with error = %d\n",__func__,xlalErrno);
	return 1;
      }
    
      /**********************************************************************************/
      /* FOR EACH SEGMENT OF DATA WITHIN THE INTERVAL */
      /**********************************************************************************/
      
      /* take each plan and combine the corresponding files appropriately into a timeseries */
      for (j=0;j<plans.length;j++) {
	
	HeaderVector header;
	REAL4TimeSeries *ts = NULL;
	INT4 k;
	LogPrintf(LOG_DEBUG,"%s : generating and ouputting timeseries %d/%d\n",__func__,j+1,plans.length);

	if (XLALCombinationPlanToREAL4TimeSeries(&ts,&header,&(plans.data[j]),pcu->pcucount[i])) {
	  LogPrintf(LOG_CRITICAL,"%s : XLALCombinationPlanToINT4TimeSeries() failed with error = %d\n",__func__,xlalErrno);
	  return 1;
	}

	/**********************************************************************************/
	/* NOW GENERATE AN OUTPUT FRAME */
	/**********************************************************************************/

	if (XLALREAL4TimeSeriesToFrame(uvar.outputdir,ts,&(plans.data[j]),clargs,&header)) {
	  LogPrintf(LOG_CRITICAL,"%s : XLALINT4TimeSeriesToFrame() failed with error = %d\n",__func__,xlalErrno);
	  return 1;
	}
	
	/* store the output frame stats */

	/* free the timeseries */
 	XLALDestroyREAL4TimeSeries(ts);
	for (k=0;k<header.length;k++) {
	  XLALFree(header.data[k].header_string);
	}
	XLALFree(header.data);
	  
      } /* end the loop over plans */
    
      /* free plan memory */
      if (plans.data) {
	for (j=0;j<plans.length;j++) {
	  XLALFree(plans.data[j].channellist.channel);
	}
	XLALFree(plans.data);
      }
      
    } /* end the if statement on whether there were any channels within the PCU interval */

    /* free mem */
    XLALFree(subframechannels->channel);
    XLALFree(subframechannels);

  } 

  /**********************************************************************************/
  /* CLEAN UP */
  /**********************************************************************************/

  /* free filelists */
  XLALFree(framechannels->channel);
  XLALFree(framechannels);
  
  /* free PCU count structure */
  XLALFree(pcu->start);
  XLALFree(pcu->end);
  XLALFree(pcu->pcucount);
  XLALFree(pcu);

  /* Free config-Variables and userInput stuff */
  LAL_CALL (LALDestroyUserVars (&status), &status);

  /* did we forget anything ? */
  LALCheckMemoryLeaks();
  LogPrintf(LOG_DEBUG,"%s : successfully checked memory leaks.\n",__func__);

  LogPrintf(LOG_DEBUG,"%s : successfully completed.\n",__func__);
  return 0;
  
}

/** Read in input user arguments
 */
void ReadUserVars(LALStatus *status,int argc,char *argv[],UserInput_t *uvar,CHAR *clargs)
{
  
  CHAR *version_string;
  INT4 i;

  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /* initialise user variables */
  uvar->inputdir = NULL; 
  uvar->outputdir = NULL; 
  uvar->pattern = NULL;
  uvar->deltat = MIN_DT;
  uvar->goodxenon = 0;

  /* ---------- register all user-variables ---------- */
  LALregBOOLUserStruct(status, 	help, 		'h', UVAR_HELP,     "Print this message");
  LALregSTRINGUserStruct(status,inputdir, 	'i', UVAR_REQUIRED, "The input frame directory"); 
  LALregSTRINGUserStruct(status,pattern, 	'p', UVAR_OPTIONAL, "The frame file name pattern"); 
  LALregSTRINGUserStruct(status,outputdir, 	'o', UVAR_REQUIRED, "The output frame file directory"); 
  LALregREALUserStruct(status,  deltat,         't', UVAR_OPTIONAL, "The output sampling time (in seconds)");
  LALregINTUserStruct(status,  goodxenon,       'x', UVAR_OPTIONAL, "Set this flag to include good xenon data");
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
int XLALReadFrameDir(FrameChannelList **framechannels,    /**< [out] a structure containing a list of all input frame channels */
		     CHAR *inputdir,                      /**< [in] the input frame directory */
		     CHAR *pattern,                       /**< [in] the input frame file name pattern */
		     INT4 goodxenon                       /**< [in] flag for including goodxenon data */  
		     )
{
  INT4 nfiles;                    /* the number of returned files from scandir */
  INT4 i;                         /* counter */
  INT4 totalchannels = 0;         /* the total number of channels */
  INT4 count = 0;                 /* counter */
  CHAR *temp,*c;                  /* temporary char pointers */
  CHAR *channelinfo;              /* string containing channel names */
  LALFrStream *fs = NULL;            /* frame stream pointer */
  glob_t pglob;
  CHAR glob_pattern[STRINGLENGTH];
  CHAR apid[APIDLENGTH];

  /* check input arguments */
  if ((*framechannels) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input frame file list structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  if (inputdir == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : Invalid input, input directory string == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
         
  /* allocate memory for filelist structure */
  if (((*framechannels) = (FrameChannelList *)XLALCalloc(1,sizeof(FrameChannelList))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for channellist structure.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* if we have a filename patttern then set the global variable */
  snprintf(apid,APIDLENGTH,"%s",pattern);
  if (pattern != NULL) snprintf(glob_pattern,STRINGLENGTH,"%s/*%s*.gwf",inputdir,apid); 
  else snprintf(glob_pattern,STRINGLENGTH,"%s/*.gwf",inputdir); 
  LogPrintf(LOG_DEBUG,"%s : searching for file pattern %s\n",__func__,glob_pattern);
  
  /* get the mode1 FITS file (apID 55) */
  if (glob(glob_pattern,0,NULL,&pglob)) {
    LogPrintf(LOG_CRITICAL,"%s : glob() failed to return a filelist.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  nfiles = pglob.gl_pathc;
  if (nfiles>0) {
    for (i=0;i<nfiles;i++) LogPrintf(LOG_DEBUG,"%s : found file %s\n",__func__,pglob.gl_pathv[i]);
  }
  else if (nfiles == -1) {
    LogPrintf(LOG_CRITICAL,"%s : scandir() failed.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  else if (nfiles == 0) {
    LogPrintf(LOG_CRITICAL,"%s : could not find any frame files in director %s.\n",__func__,inputdir);
    XLAL_ERROR(XLAL_EINVAL);
  }
  
  /* open each file and count how many channels we have in total */
  for (i=0;i<nfiles;i++) {

    /* check we do not allow goodxenon and the file has XENO in the filename then we ignore the file */
    if ( ! ( (!goodxenon) && (strstr(pglob.gl_pathv[i],"XENO") != NULL) ) ) {

      /* open the frame file */
      if ((fs = XLALFrStreamOpen(inputdir,pglob.gl_pathv[i])) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s : unable to open FS46 frame file %s.\n",__func__,pglob.gl_pathv[i]);
	XLAL_ERROR(XLAL_EINVAL);
      }
      LogPrintf(LOG_DEBUG,"%s : opened frame file %s.\n",__func__,pglob.gl_pathv[i]);
      
      /* count the number of channels in the frame file */
      /* this returns a string containing the channel names */
      channelinfo = FrFileIGetChannelList(fs->file,(INT4)fs->flist->t0);
      
      /* find all instances of the keyword PROC */
      c = channelinfo;
      while ((temp = strstr(c,"PROC")) != NULL) {
	totalchannels++;
	c = temp + 1;
      }
      free(channelinfo);
      
      /* close the frame file */
      XLALFrStreamClose(fs);
  
    }

  }
  
  if (totalchannels == 0) {
    LogPrintf(LOG_CRITICAL,"%s : could not find any PROC channels in dir %s.\n",__func__,inputdir);
    XLAL_ERROR(XLAL_EINVAL);
  }
  LogPrintf(LOG_DEBUG,"%s : counted %d channels in total.\n",__func__,totalchannels);

  /* allocate memory for channellist metadata */
  if (((*framechannels)->channel = (FrameChannel *)XLALCalloc(totalchannels,sizeof(FrameChannel))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FrameChannel structures.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* open each file (again) and extract info */
  for (i=0;i<nfiles;i++) {
   
    /* check we do not allow goodxenon and the file has XENO in the filename then we ignore the file */
    if ( ! ( (!goodxenon) && (strstr(pglob.gl_pathv[i],"XENO") != NULL) ) ) {
      
      CHAR *temp_obsid;
      
      /* open the frame file */
      if ((fs = XLALFrStreamOpen(inputdir,pglob.gl_pathv[i])) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s : unable to open FS46 frame file %s.\n",__func__,pglob.gl_pathv[i]);
	XLAL_ERROR(XLAL_EINVAL);
      }
      LogPrintf(LOG_DEBUG,"%s : opened frame file %s (file %d/%d).\n",__func__,pglob.gl_pathv[i],i+1,nfiles);
      
      /* get history information from this file - mainly to extract the OBS_ID */
      {
	CHAR * header_string = NULL;
	CHAR *c1,*c2,*c3;
	UINT4 length;
	
	/* read in history from this file */
	if (XLALReadFrameHistory(&header_string,fs->file)) {
	  LogPrintf(LOG_CRITICAL,"%s : XLALReadFrameHistory() unable to read history from frame file %s.\n",__func__,pglob.gl_pathv[i]);
	  XLAL_ERROR(XLAL_EINVAL);
	}
	LogPrintf(LOG_DEBUG,"%s : read history field from file %s.\n",__func__,pglob.gl_pathv[i]);
	
	/* extract OBS_ID */
	c1 = strstr(header_string,"OBS_ID");
	c2 = strstr(c1,"'");
	c3 = strstr(c2+1,"'");
	length = strlen(c2) - strlen(c3);
	temp_obsid = (CHAR *)XLALCalloc(length+1,sizeof(CHAR));
	snprintf(temp_obsid,length,"%s",c2+1);
	XLALFree(header_string);
	
      }
      
      /* count the number of channels in the frame file */
      /* this returns a string containing the channel names */
      channelinfo = FrFileIGetChannelList(fs->file,(INT4)fs->flist->t0);
      
      /* get proc channel name(s) from channel info string */
      {
	char *channelname;
	char *start,*end;
	int length;
	INT4TimeSeries ts;
	
	/* find all instances of the keyword PROC */
	c = channelinfo;
	while ((temp = strstr(c,"PROC")) != NULL) {
	  
	  /* extract the first channel name from the string */
	  start = temp + 5;  
	  end = strstr(start,"\t");
	  length = strlen(start) - strlen(end) + 1;
	  channelname = (CHAR *)XLALCalloc(length+1,sizeof(CHAR));
	  snprintf(channelname,length,"%s",start);
	  snprintf(ts.name,LALNameLength,"%s",channelname);
	  
	  LogPrintf(LOG_DEBUG,"%s : read channel name as %s\n",__func__,channelname);
	  
	  /* read in timeseries metadata from this channel - specifically to get deltaT */
	  if (XLALFrStreamGetINT4TimeSeriesMetadata(&ts,fs)) {
	    LogPrintf(LOG_CRITICAL,"%s : unable to read channel %s from frame file %s.\n",__func__,ts.name,(*framechannels)->channel[i].filename);
	    XLAL_ERROR(XLAL_EINVAL);
	  }
	  LogPrintf(LOG_DEBUG,"%s : read deltaT as %6.12f\n",__func__,ts.deltaT);
	  
	  snprintf((*framechannels)->channel[count].channelname,length,"%s",ts.name);	
	  (*framechannels)->channel[count].dt = ts.deltaT;
	  LogPrintf(LOG_DEBUG,"%s : read channel name as %s (channel %d/%d)\n",__func__,(*framechannels)->channel[count].channelname,count+1,totalchannels);
	  
	  /* fill in the channel name and extract start and duration of frame from stream */
	  XLALGPSSetREAL8(&((*framechannels)->channel[count].epoch),(REAL8)fs->flist->t0);
	  XLALGPSSetREAL8(&((*framechannels)->channel[count].end),(REAL8)(fs->flist->t0 + fs->flist->dt));
	  (*framechannels)->channel[count].duration = (REAL8)fs->flist->dt;
	  
	  /* extract detconfig, energy, LLD and data type information from the channelname */
	  /* the format is X1:<MODE>-<COLNAME>-<LLD>-<DETCONFIG>-<MINENERGY>_<MAXENERGY> */
	  /* the data type is indicated by either "Event" or "Cnt" being present in the channel name */
	  {
	    CHAR *c1 = ts.name;             /* points to start of string */
	    CHAR *c2 = strstr(c1,"-");      /* points to first instance of "-" after mode */
	    CHAR *c3 = strstr(c2+1,"-");    /* points to first instance of "-" after colname */
	    CHAR *c4 = strstr(c3+1,"-");    /* points to first instance of "-" after LLD */
	    CHAR *c5 = strstr(c4+1,"-");    /* points to first instance of "-" after detconfig */
	    CHAR *subs,*e1,*e2;
	    INT4 sublen,elen;
	    
	    /* check for type */
	    if (strstr(ts.name,"Cnt") != NULL) (*framechannels)->channel[count].type = ARRAY;
	    else if (strstr(ts.name,"Event") != NULL) (*framechannels)->channel[count].type = EVENT;
	    else {
	      LogPrintf(LOG_CRITICAL,"%s : unable to read data type from channel %s in file %s.\n",__func__,ts.name,(*framechannels)->channel[i].filename);
	      XLAL_ERROR(XLAL_EINVAL);
	    }
	    
	    /* extract LLD info */
	    sublen = 2;
	    subs = (CHAR *)XLALCalloc(sublen,sizeof(CHAR));
	    snprintf(subs,sublen,"%s",c3+1);	 
	    (*framechannels)->channel[count].lld = atoi(subs);
	    XLALFree(subs);
	    
	    /* extract detconfig info */
	    sublen = strlen(c4) - strlen(c5);   
	    subs = (CHAR *)XLALCalloc(sublen,sizeof(CHAR));
	    snprintf(subs,sublen,"%s",c4+1);
	    
	    /* check if there is are instances of 0,1,2,3,4,5 in the string */
	    if (strstr(subs,"0") != NULL) (*framechannels)->channel[count].detconfig[0]= 1;
	    if (strstr(subs,"1") != NULL) (*framechannels)->channel[count].detconfig[1]= 1;
	    if (strstr(subs,"2") != NULL) (*framechannels)->channel[count].detconfig[2]= 1;
	    if (strstr(subs,"3") != NULL) (*framechannels)->channel[count].detconfig[3]= 1;
	    if (strstr(subs,"4") != NULL) (*framechannels)->channel[count].detconfig[4]= 1;      
	    XLALFree(subs);
	    
	    /* extract energy info */
	    sublen = strlen(c5);
	    subs = (CHAR *)XLALCalloc(sublen,sizeof(CHAR));
	    snprintf(subs,sublen,"%s",c5+1);
	    c1 = strstr(subs,"_");
	    elen = strlen(subs) - strlen(c1) + 1;
	    e1 = (CHAR *)XLALCalloc(elen,sizeof(CHAR));
	    snprintf(e1,elen,"%s",subs);
	    (*framechannels)->channel[count].energy[0] = atoi(e1);
	    elen = strlen(c1);	
	    e2 = (CHAR *)XLALCalloc(elen,sizeof(CHAR));
	    snprintf(e2,elen,"%s",c1+1);	
	    (*framechannels)->channel[count].energy[1] = atoi(e2);
	    XLALFree(subs);
	    XLALFree(e1);
	    XLALFree(e2);
	    
	  } /* end of channel name info extraction */
	  
	  /* copy filename and history to output structure */
	  snprintf((*framechannels)->channel[count].filename,STRINGLENGTH,"%s",pglob.gl_pathv[i]);
	  snprintf((*framechannels)->channel[count].OBS_ID,LALNameLength,"%s",temp_obsid);

	  LogPrintf(LOG_DEBUG,"%s : channel number %d.\n",__func__,count);
	  LogPrintf(LOG_DEBUG,"%s : frame filename = %s.\n",__func__,(*framechannels)->channel[count].filename);
	  LogPrintf(LOG_DEBUG,"%s : frame channelname = %s.\n",__func__,(*framechannels)->channel[count].channelname);
	  LogPrintf(LOG_DEBUG,"%s : frame start time = %d.\n",__func__,(*framechannels)->channel[count].epoch.gpsSeconds);
	  LogPrintf(LOG_DEBUG,"%s : frame end time = %d.\n",__func__,(*framechannels)->channel[count].end.gpsSeconds);
	  LogPrintf(LOG_DEBUG,"%s : frame duration = %f.\n",__func__,(*framechannels)->channel[count].duration);
	  LogPrintf(LOG_DEBUG,"%s : frame lld flag = %d.\n",__func__,(*framechannels)->channel[count].lld);
	  LogPrintf(LOG_DEBUG,"%s : frame data type = %d.\n",__func__,(*framechannels)->channel[count].type);
	  LogPrintf(LOG_DEBUG,"%s : frame min energy = %d.\n",__func__,(*framechannels)->channel[count].energy[0]);
	  LogPrintf(LOG_DEBUG,"%s : frame max energy = %d.\n",__func__,(*framechannels)->channel[count].energy[1]);
	  LogPrintf(LOG_DEBUG,"%s : frame detconfig = %d %d %d %d %d.\n",__func__,(*framechannels)->channel[count].detconfig[0],
		    (*framechannels)->channel[count].detconfig[1],(*framechannels)->channel[count].detconfig[2],
		    (*framechannels)->channel[count].detconfig[3],(*framechannels)->channel[count].detconfig[4]);
	  LogPrintf(LOG_DEBUG,"%s : frame dt = %6.12f.\n",__func__,(*framechannels)->channel[count].dt);
	  count++;
	  
	  /* free mem */
	  XLALFree(channelname);
	  c = temp + 1;
	  
	}
	
	/* free mem */
	free(channelinfo);
	
      }
      
      /* close the frame file and free history */
      XLALFrStreamClose(fs);
      XLALFree(temp_obsid);
    
    }

  }

  /* add inputdir and length to output structure */
  strncpy((*framechannels)->dir,inputdir,STRINGLENGTH);
  (*framechannels)->length = totalchannels;

  /* free original filelist */
  globfree(&pglob);
 
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}


/** Read in pcu counts from FS46 files and generate a vector of PCU numbers 
 */
int XLALReadGoodPCUInterval(GoodPCUIntervals **pcu,              /**< [out] the PCU interval information */
			    FrameChannelList *framechannels      /**< [in] the framefile list */
			    )
{
  UINT4 i,j;
  UINT4 count = 0;
  CHAR **tempFS46 = NULL;         /* used to store unique FS46 filenames */ 
  UINT4 nfiles = 0;                /* the number of unique FS46 files */
  INT4 idx = 0;

  /* check input arguments */
  if (framechannels == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input frame file channel structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ((*pcu) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output pcu structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
 
  /* allocate memory for the output structure */
  if (((*pcu) = (GoodPCUIntervals *)XLALCalloc(1,sizeof(GoodPCUIntervals))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to allocate memory for GoodPCUIntervals structure.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*pcu)->start = (LIGOTimeGPS *)XLALCalloc(framechannels->length,sizeof(LIGOTimeGPS))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to allocate memory for pcu start data.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*pcu)->end = (LIGOTimeGPS *)XLALCalloc(framechannels->length,sizeof(LIGOTimeGPS))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to allocate memory for pcu end data.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*pcu)->pcucount = (UINT2 *)XLALCalloc(framechannels->length,sizeof(UINT2))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to allocate memory for pcu counts data.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  (*pcu)->length = framechannels->length;
 
  /* allocate memory for temporary FS46 filenames */
  if ((tempFS46 = (CHAR **)XLALCalloc(framechannels->length,sizeof(CHAR *))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to allocate memory for temporary FS46 filename pointers.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  
  /* loop over each frame channel and find unique FS46 filenames */
  for (i=0;i<framechannels->length;i++) {
   
    /* check if it is an FS46 file */
    if (strstr(framechannels->channel[i].filename,"_FS46") != NULL) {
      
      INT4 prod = 1;
      for (j=0;j<count;j++) prod *= strcmp(framechannels->channel[i].filename,tempFS46[j]);
      
      /* copy the filename to temp storage if unique */
      if (prod) {
	if ((tempFS46[count] = (CHAR *)XLALCalloc(strlen(framechannels->channel[i].filename)+1,sizeof(CHAR))) == NULL ) {
	  LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for temporary FS46 filename.\n",__func__);
	  XLAL_ERROR(XLAL_ENOMEM);
	}
	strncpy(tempFS46[count],framechannels->channel[i].filename,strlen(framechannels->channel[i].filename));
	count++;
      }
    }
  }
  nfiles = count;
  LogPrintf(LOG_DEBUG,"%s: found %d unique FS46 files.\n",__func__,nfiles);
  
 /* loop over each unique FS46 filename */
  for (j=0;j<nfiles;j++) {

    UINT2 pcucount = 0;
    INT4 FS46channels = 0;
    LALFrStream *fs = NULL;
    LIGOTimeGPS epoch;
    LIGOTimeGPS end;
    REAL8 duration;
    
    /* loop over each channel */
    for (i=0;i<framechannels->length;i++) {
      
      /* check if the filenames match - if so, open the channel and count the PCUs */
      if (strcmp(framechannels->channel[i].filename,tempFS46[j]) == 0) {	
	
      	INT4TimeSeries *ts = NULL;

	/* open the frame file */
	if ((fs = XLALFrStreamOpen(framechannels->dir,framechannels->channel[i].filename)) == NULL) {
	  LogPrintf(LOG_CRITICAL,"%s: unable to open FS46 frame file %s.\n",__func__,framechannels->channel[i].filename);
	  XLAL_ERROR(XLAL_EINVAL);
	}
	LogPrintf(LOG_DEBUG,"%s: opened frame file %s.\n",__func__,framechannels->channel[i].filename);
	
	/* extract start and duration of frame from stream */
	XLALGPSSetREAL8(&epoch,(REAL8)fs->flist->t0);
	duration = (REAL8)fs->flist->dt;
	XLALGPSSetREAL8(&end,(REAL8)fs->flist->t0 + duration);
	LogPrintf(LOG_DEBUG,"%s: frame start time = %d.\n",__func__,epoch.gpsSeconds);
	LogPrintf(LOG_DEBUG,"%s: frame duration = %f.\n",__func__,duration);
	
	/* seek to the start of the frame */
	if (XLALFrStreamSeek(fs,&epoch)) {
	  LogPrintf(LOG_CRITICAL,"%s: unable to seek to start of frame file %s.\n",__func__,framechannels->channel[i].filename);
	  XLAL_ERROR(XLAL_EINVAL);
	}

	/* read in timeseries from this file - final arg is limit on length of timeseries (0 = no limit) */
	if ((ts = XLALFrStreamReadINT4TimeSeries(fs,framechannels->channel[i].channelname,&epoch,duration,0)) == NULL) {
	  LogPrintf(LOG_CRITICAL,"%s: unable to read channel %s from frame file %s.\n",__func__,framechannels->channel[i].channelname,framechannels->channel[i].filename);
	  XLAL_ERROR(XLAL_EINVAL);
	}
	LogPrintf(LOG_DEBUG,"%s: reading channel %s.\n",__func__,framechannels->channel[i].channelname);
	LogPrintf(LOG_DEBUG,"%s: read N = %d, dt = %6.12f, T = %f\n",__func__,ts->data->length,ts->deltaT,ts->data->length*ts->deltaT);
	
	/* determine whether each PCU was working or not */
	/* from ED MORGANS email ----------------------------------------------------------------------------------------------*/
	/* All 5 PCUs are all still functional to some extent. Any given observation may have  between 1 and 5 PCU turned on. */
	/* The average number of PCUs has been declining more or less linearly over the last 13 years. */
	/* The way I  determine the number of PCUs is to look apid 70. The first 5 columns contain the count rate per PCU per .125 s. */
	/* So I rebin to a longer time (e.g. 128s) count the number of non zero columns (of the 5). Actually even off detectors occasionally */
	/* get a count, count the number of columns with cts > a few.  */
	{
	  REAL8 mean = 0;
	  UINT4 k;
	  for (k=0;k<ts->data->length;k++) mean += ts->data->data[k];
	  mean /= (REAL8)ts->data->length;
	  LogPrintf(LOG_DEBUG,"%s: computed mean photon count as %f per bin\n",__func__,mean);
	  if ( (mean > PCUCOUNTTHRESH) && (duration > MINFRAMELENGTH) ) pcucount++;
	}
	
	/* free the PCU timeseries data */
	XLALDestroyINT4TimeSeries(ts);
	
	/* close the frame file */
	XLALFrStreamClose(fs);

	/* increment the number of channels read from this file */
	FS46channels++;

      } /* end of if statement on whether the file is an FS46 file */

    } /* end of loop over every channel */
    
    /* fill in output structure if there are any operational PCUs and we've read 5 channels from this file */
    if (pcucount && (FS46channels == NPCU)) {
      (*pcu)->start[idx].gpsSeconds = epoch.gpsSeconds;
      (*pcu)->start[idx].gpsNanoSeconds = epoch.gpsNanoSeconds;
      (*pcu)->end[idx].gpsSeconds = end.gpsSeconds;
      (*pcu)->end[idx].gpsNanoSeconds = end.gpsNanoSeconds;
      (*pcu)->pcucount[idx] = (UINT2)pcucount;
                  
      /* increment file count */
      idx++;
    }
    LogPrintf(LOG_DEBUG,"%s : computed PCU count for this interval as %d.\n",__func__,pcucount);
  
  }  /* end of loop over unique FS46 file names */
  
  /* check if any FS46 files were found */
  if (idx == 0) {
    LogPrintf(LOG_CRITICAL,"%s: could not find any FS46 (PCU) files.  Exiting.\n",__func__);
    exit(0);
  }

  /* resize the output vectors */
  (*pcu)->length = idx;
  if (((*pcu)->start = (LIGOTimeGPS *)XLALRealloc((*pcu)->start,(*pcu)->length*sizeof(LIGOTimeGPS))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to resize memory for pcu start data.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*pcu)->end = (LIGOTimeGPS *)XLALRealloc((*pcu)->end,(*pcu)->length*sizeof(LIGOTimeGPS))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to resize memory for pcu end data.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  if (((*pcu)->pcucount = (UINT2 *)XLALRealloc((*pcu)->pcucount,(*pcu)->length*sizeof(UINT2))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: failed to resize memory for pcu counts data.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  
  /* free mem */
  for (i=0;i<framechannels->length;i++) XLALFree(tempFS46[i]);
  XLALFree(tempFS46);
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/** Finds a subset of frame files within a given time interval 
 */
int XLALFindFramesInInterval(FrameChannelList **subframechannels,   /**< [out] a list of channel names containing data within the interval */
			     FrameChannelList *framechannels,       /**< [in] the frame channel list */
			     LIGOTimeGPS intstart,                  /**< [in] the interval start time */
			     LIGOTimeGPS intend,                    /**< [in] the interval end time */
			     UINT2 npcus                            /**< [in] the number of operational PCUs */ 
			     )
{
  INT4 i;
  INT4 count = 0;

  /* check input arguments */
  if (framechannels == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input frame channel list structure == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if ((*subframechannels) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output FileList structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (XLALGPSCmp(&intstart,&intend) >= 0) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, the interval start time %d >= %d interval end time.\n",__func__,intstart.gpsSeconds,intend.gpsSeconds);
    XLAL_ERROR(XLAL_EINVAL);
  }
  
  /* allocate memory for output filelist structure */
  if (((*subframechannels) = (FrameChannelList *)XLALCalloc(1,sizeof(FrameChannelList))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for filelist structure.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* allocate memory for filelist metadata */
  if (((*subframechannels)->channel = (FrameChannel *)XLALCalloc(framechannels->length,sizeof(FrameChannel))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FrameChannel structures.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }

  /* loop over the input frame file names */
  for (i=0;i<(INT4)framechannels->length;i++) {
    
    LIGOTimeGPS *framestart = &(framechannels->channel[i].epoch);
    LIGOTimeGPS *frameend = &(framechannels->channel[i].end);
    
    /* check if it is NOT an FS46 file using the filename */
    if (strstr(framechannels->channel[i].filename,"_FS46-") == NULL) {
      	
      /* check if file overlaps with interval */
      if ( ! (( XLALGPSCmp(framestart,&intend) > 0 ) || ( XLALGPSCmp(frameend,&intstart) < 0 )) ) {

	/* check for timing consistency */	
	if ( (framechannels->channel[i].dt >= MIN_DT) && (framechannels->channel[i].dt <= MAX_DT) ) {
	
	  /* compute actual overlap */
	  REAL8 x = MYMAX(intstart.gpsSeconds,framestart->gpsSeconds);
	  REAL8 y = MYMIN(intend.gpsSeconds,frameend->gpsSeconds);
	  REAL8 overlap = y-x;
	  
	  /* record overlapping frames and channels */
	  memcpy(&((*subframechannels)->channel[count].epoch),&(framechannels->channel[i].epoch),sizeof(LIGOTimeGPS));
	  memcpy(&((*subframechannels)->channel[count].end),&(framechannels->channel[i].end),sizeof(LIGOTimeGPS));
	  (*subframechannels)->channel[count].duration = framechannels->channel[i].duration;
	  (*subframechannels)->channel[count].lld = framechannels->channel[i].lld;
	  (*subframechannels)->channel[count].type = framechannels->channel[i].type;
	  (*subframechannels)->channel[count].energy[0] = framechannels->channel[i].energy[0];
	  (*subframechannels)->channel[count].energy[1] = framechannels->channel[i].energy[1];
	  memcpy(&((*subframechannels)->channel[count].detconfig),&(framechannels->channel[i].detconfig),NPCU*sizeof(INT4));
	  (*subframechannels)->channel[count].dt = framechannels->channel[i].dt;
	  strncpy((*subframechannels)->channel[count].OBS_ID,framechannels->channel[i].OBS_ID,STRINGLENGTH);
	  LogPrintf(LOG_DEBUG,"%s : Int [%d->%d] : overlapping frame [%d->%d] = %.0f overlap.\n",__func__,intstart.gpsSeconds,intend.gpsSeconds,framestart->gpsSeconds,frameend->gpsSeconds,overlap);
	  
	  /* copy filename and channelname information */
	  strncpy((*subframechannels)->channel[count].filename,framechannels->channel[i].filename,STRINGLENGTH);
	  strncpy((*subframechannels)->channel[count].channelname,framechannels->channel[i].channelname,STRINGLENGTH);
	  
	  /* increment the count of overlapping channels */
	  count++;
	  
	}
	else {
	  LogPrintf(LOG_DEBUG,"%s : channel %s failed timing consitency check (dt = %f).\n",__func__,framechannels->channel[i].channelname,framechannels->channel[i].dt);
	}
	
      }
      /* else { */
/* 	LogPrintf(LOG_DEBUG,"%s : frame range %d -> %d doesn't overlap with interval %d -> %d.\n",__func__, */
/* 		  framechannels->channel[i].epoch.gpsSeconds,framechannels->channel[i].end.gpsSeconds, */
/* 		  intstart.gpsSeconds,intend.gpsSeconds); */
/*       } */

    }
  
  }
  
  /* fill in additional information */
  (*subframechannels)->length = count;
  (*subframechannels)->npcus = npcus;  
  strncpy((*subframechannels)->dir,framechannels->dir,STRINGLENGTH);
  LogPrintf(LOG_DEBUG,"%s : found %d overlapping channels.\n",__func__,(*subframechannels)->length);

  /* resize the output vectors */
  if ((*subframechannels)->length) {
    if (((*subframechannels)->channel = (FrameChannel *)XLALRealloc((*subframechannels)->channel,(*subframechannels)->length*sizeof(FrameChannel))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: failed to resize memory for subframefile list.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
  }
  else {
    XLALFree((*subframechannels)->channel);
    (*subframechannels)->channel = NULL;
  }
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/** Finds a subset of frame files within a given time interval 
 */
int XLALCreateCombinationPlan(FrameCombinationPlanVector *plans,             /**< [out] a plan of how to combine the frames */
			      FrameChannelList *framechannels,               /**< [in] the framefile list */
			      LIGOTimeGPS *intstart,                          /**< [in] the interval start time */
			      LIGOTimeGPS *intend                             /**< [in] the interval end time */
			      )
{
  LIGOTimeGPS *earlieststart = &(framechannels->channel[0].epoch);
  LIGOTimeGPS *latestend = &(framechannels->channel[0].end);
  INT4 i,j,k,s,r;                         /* counters */
  LIGOTimeGPSVector *epochlist;
  INT4 q = 0;
  REAL8 maxspan = 0.0;
  REAL8 intspan = 0.0;
  INT4 N = 75;
  INT4 pcount = 0;

  /* MOST OF THE FIRST PART IS JUST FOR PLOTTING IN ASCII TO THE SCREEN */

  /* loop over each channel to find the earliest and latest GPS time for this segment */
  for (i=0;i<(INT4)framechannels->length;i++) {
    
    if (XLALGPSCmp(&(framechannels->channel[i].epoch),earlieststart) < 0) earlieststart = &(framechannels->channel[i].epoch);
    if (XLALGPSCmp(&(framechannels->channel[i].end),latestend) > 0) latestend = &(framechannels->channel[i].end);

  }
  /* check against the actual interval start and end */
  if (XLALGPSCmp(intstart,earlieststart) < 0) earlieststart = intstart;
  if (XLALGPSCmp(intend,latestend) > 0) latestend = intend;

  maxspan = XLALGPSDiff(latestend,earlieststart);
  intspan = XLALGPSDiff(intend,intstart);
  for (i=0;i<(INT4)framechannels->length;i++) LogPrintf(LOG_DEBUG,"%s : for this interval we have the channel %s\n",__func__,framechannels->channel[i].channelname);
  
  /* allocate memory for the timestamps list */
  epochlist = XLALCreateTimestampVector(2*framechannels->length);

  /* draw the interval in ascii */
  {
    INT4 sidx = floor(0.5 + N*(XLALGPSGetREAL8(intstart) - XLALGPSGetREAL8(earlieststart))/maxspan);
    INT4 eidx = floor(0.5 + N*(XLALGPSGetREAL8(intend) - XLALGPSGetREAL8(earlieststart))/maxspan);
    for (j=0;j<sidx;j++) fprintf(stdout," ");
    fprintf(stdout,"|");
    for (j=sidx+1;j<eidx-1;j++) fprintf(stdout,"#");
    fprintf(stdout,"| INTERVAL (%.0f)\n",intspan);
  }

  /* loop over each channel */
  for (i=0;i<(INT4)framechannels->length;i++) {
    
    /* draw (in ascii) the data stretches */
    {
      INT4 sidx = floor(0.5 + N*(XLALGPSGetREAL8(&(framechannels->channel[i].epoch)) - XLALGPSGetREAL8(earlieststart))/maxspan);
      INT4 eidx = floor(0.5 + N*(XLALGPSGetREAL8(&(framechannels->channel[i].end)) - XLALGPSGetREAL8(earlieststart))/maxspan);
      
      for (j=0;j<sidx;j++) fprintf(stdout," ");
      fprintf(stdout,"|");
      for (j=sidx+1;j<eidx-1;j++) fprintf(stdout,"-");
      fprintf(stdout,"| dt = %6.6f detconfig = %d%d%d%d%d energy = %d-%d lld = %d type = %d\n",
	      framechannels->channel[i].dt,framechannels->channel[i].detconfig[0],framechannels->channel[i].detconfig[1],
	      framechannels->channel[i].detconfig[2],framechannels->channel[i].detconfig[3],framechannels->channel[i].detconfig[4],
	      framechannels->channel[i].energy[0],framechannels->channel[i].energy[1],framechannels->channel[i].lld,framechannels->channel[i].type);
      
    }
    
    /* add start and end to a list of times that will then be sorted */
    {
      LIGOTimeGPS *tempepoch = &(framechannels->channel[i].epoch);
      LIGOTimeGPS *tempend = &(framechannels->channel[i].end);
      if (XLALGPSCmp(intstart,&(framechannels->channel[i].epoch)) > 0) tempepoch = intstart;
      epochlist->data[q].gpsSeconds = tempepoch->gpsSeconds;
      epochlist->data[q].gpsNanoSeconds = tempepoch->gpsNanoSeconds;
      if (XLALGPSCmp(intend,&(framechannels->channel[i].end)) < 0) tempend = intend;
      epochlist->data[q+1].gpsSeconds = tempend->gpsSeconds;
      epochlist->data[q+1].gpsNanoSeconds = tempend->gpsNanoSeconds;
    }
    q = q + 2;

  }
       
  /* sort the epochs */
  qsort(epochlist->data, epochlist->length, sizeof(LIGOTimeGPS), compareGPS);
  
  /* allocate memory for the current plan under the assumption that there will be epochlist->length-1 seperate plans */
  plans->length = epochlist->length-1;
  if ((plans->data = (FrameCombinationPlan *)XLALCalloc(plans->length,sizeof(FrameCombinationPlan))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FrameCombinationPlan.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  
  /* loop over each epoch span and determine which files belong within it */
  for (i=0;i<(INT4)epochlist->length-1;i++) {
    
    FrameChannelList tempchannellist;
    INT4Vector *badidx = NULL;
    INT4 nchannels = 0;
    INT4 count = 0;
    REAL8 span = 0;                /* span of an interval */

    LogPrintf(LOG_DEBUG,"%s : working on span #%d : %d -> %d (%d)\n",__func__,i,epochlist->data[i].gpsSeconds,epochlist->data[i+1].gpsSeconds,epochlist->data[i+1].gpsSeconds-epochlist->data[i].gpsSeconds);
    
    /* if there is any data spanned greater than the minimum frame length */
    span = XLALGPSDiff(&(epochlist->data[i+1]),&(epochlist->data[i]));
    if (span >= MINFRAMELENGTH) {
    
      /* count number of coincident channels */
      for (k=0;k<(INT4)framechannels->length;k++) {
	if ( !((XLALGPSCmp(&(framechannels->channel[k].epoch),&(epochlist->data[i+1])) >= 0) || (XLALGPSCmp(&(framechannels->channel[k].end),&(epochlist->data[i])) <= 0)) ) {
	  LogPrintf(LOG_DEBUG,"%s : coincident channel %d = %s\n",__func__,k,framechannels->channel[k].channelname);
	  nchannels++;
	}
      }
      LogPrintf(LOG_DEBUG,"%s : %d coincident channels for this subinterval\n",__func__,nchannels);
      
      /* allocate memory for the temporary filelist */
      tempchannellist.length = nchannels;
      tempchannellist.npcus = framechannels->npcus;
      snprintf(tempchannellist.dir,STRINGLENGTH,"%s",framechannels->dir);
      if ((tempchannellist.channel = (FrameChannel *)XLALCalloc(tempchannellist.length,sizeof(FrameChannel))) == NULL) {
	LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for temp FrameChannel structures.\n",__func__);
	XLAL_ERROR(XLAL_ENOMEM);
      }
      LogPrintf(LOG_DEBUG,"%s : allocated memory for the temp channel list\n",__func__);
      
      /* allocate memory for the list of bad indices and initialise */
      if ((badidx = XLALCreateINT4Vector(tempchannellist.length))==NULL) {
	LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for badidx vector.\n",__func__);
	XLAL_ERROR(XLAL_ENOMEM);
      }
      LogPrintf(LOG_DEBUG,"%s : allocated memory for the badidx vector\n",__func__);
      
      /* initialise badidx */
      for (k=0;k<(INT4)badidx->length;k++) badidx->data[k] = 0;
      
      /* find the coincident files again and record the filenames and channel details */
      for (k=0;k<(INT4)framechannels->length;k++) {
	if ( !((XLALGPSCmp(&(framechannels->channel[k].epoch),&(epochlist->data[i+1])) >= 0) || (XLALGPSCmp(&(framechannels->channel[k].end),&(epochlist->data[i])) <= 0)) ) {
	  
	  /* fill in temporary filelist */
	  strncpy(tempchannellist.channel[count].channelname,framechannels->channel[k].channelname,STRINGLENGTH);
	  strncpy(tempchannellist.channel[count].filename,framechannels->channel[k].filename,STRINGLENGTH);
	  memcpy(&(tempchannellist.channel[count].epoch),&(framechannels->channel[k].epoch),sizeof(LIGOTimeGPS));
	  memcpy(&(tempchannellist.channel[count].end),&(framechannels->channel[k].end),sizeof(LIGOTimeGPS));
	  tempchannellist.channel[count].dt = framechannels->channel[k].dt;
	  tempchannellist.channel[count].duration = framechannels->channel[k].duration;
	  tempchannellist.channel[count].energy[0] = framechannels->channel[k].energy[0];
	  tempchannellist.channel[count].energy[1] = framechannels->channel[k].energy[1];
	  snprintf(tempchannellist.channel[count].OBS_ID,STRINGLENGTH,"%s",framechannels->channel[k].OBS_ID);
	  memcpy(&(tempchannellist.channel[count].detconfig),&(framechannels->channel[k].detconfig),NPCU*sizeof(INT4));
	  tempchannellist.channel[count].lld = framechannels->channel[k].lld;
	  tempchannellist.channel[count].type = framechannels->channel[k].type;
	  
	  count++;
	  
	}
      }
    
      LogPrintf(LOG_DEBUG,"%s : recorded the coincident filenames in temp channel list\n",__func__);
    
      /* now check for single lld files */
      if ( (nchannels == 1) && (tempchannellist.channel[0].lld == 1) ) {
	badidx->data[0] = 1;
	LogPrintf(LOG_DEBUG,"%s : found single lld channel\n",__func__);
      }
      LogPrintf(LOG_DEBUG,"%s : checked for single lld channels\n",__func__);

      /* check OBS ID consistency */
      if (tempchannellist.length) { 
	
	/* check each pair of overlapping channels */
	for (k=0;k<(INT4)tempchannellist.length;k++) {
	  for (s=0;s<(INT4)k;s++) {
	    printf("comparing files %s and %s OBS_IDs = %s %s\n",tempchannellist.channel[k].filename,tempchannellist.channel[s].filename,tempchannellist.channel[k].OBS_ID,tempchannellist.channel[s].OBS_ID);
	    if (strcmp(tempchannellist.channel[k].OBS_ID,tempchannellist.channel[s].OBS_ID)) {
	      LogPrintf(LOG_CRITICAL,"%s : interval contains different OBS IDs (%s and %s) !!  Exiting.\n",__func__,tempchannellist.channel[k].OBS_ID,tempchannellist.channel[s].OBS_ID);
	      XLAL_ERROR(XLAL_EINVAL);
	    }
	  }
	}
      }

      /* now check for energy consistency */
      if (tempchannellist.length) { 

	/* check each pair of overlapping channels */
	for (k=0;k<(INT4)tempchannellist.length;k++) {
	  for (s=0;s<(INT4)k;s++) {
	  
	    INT4 pcusum = 0;
	    
	    /* check for energy overlap if not ldd */
	    if ( ( tempchannellist.channel[k].lld == 0 ) && ( tempchannellist.channel[s].lld == 0 ) ) {
	     
	      /* check if we're comparing event data with array data - apparently we can combine in this case */
	      if ( tempchannellist.channel[k].type == tempchannellist.channel[s].type ) {
	
		/* check for energy overlap */
		if ( ! ( ( tempchannellist.channel[k].energy[0] > tempchannellist.channel[s].energy[1] ) || 
			 ( tempchannellist.channel[k].energy[1] < tempchannellist.channel[s].energy[0] ) ) ) {
		
		  /* check for PCU overlap - we may have complimentary PCUs operational i.e 01010 & 10101 */
		  for (r=0;r<NPCU;r++) pcusum += tempchannellist.channel[k].detconfig[r]*tempchannellist.channel[s].detconfig[r];
		  LogPrintf(LOG_DEBUG,"%s : PCU overlap product = %d\n",__func__,pcusum);

		  /* if the PCUs overlap then we have to choose which channel to use based on the energy range */
		  if (pcusum) {

		    /* determine which file to keep based on largest energy range */
		    INT4 range1 = tempchannellist.channel[k].energy[1] - tempchannellist.channel[k].energy[0];
		    INT4 range2 = tempchannellist.channel[s].energy[1] - tempchannellist.channel[s].energy[0];
		    INT4 highenergyidx = tempchannellist.channel[k].energy[1] > tempchannellist.channel[s].energy[1] ? k : s;
		    
		    if ( range1 > range2 ) {
		      badidx->data[s] = 1;
		      LogPrintf(LOG_DEBUG,"%s : selected %d as badidx -> energy overlap [%d - %d],[%d - %d]\n",
				__func__,s,tempchannellist.channel[k].energy[0],tempchannellist.channel[k].energy[1],
				tempchannellist.channel[s].energy[0],tempchannellist.channel[s].energy[1]);
		    }
		    else if ( range2 > range1 ) {
		      badidx->data[k] = 1;
		      LogPrintf(LOG_DEBUG,"%s : selected %d as badidx -> energy overlap [%d - %d],[%d - %d]\n",
				__func__,k,tempchannellist.channel[k].energy[0],tempchannellist.channel[k].energy[1],
				tempchannellist.channel[s].energy[0],tempchannellist.channel[s].energy[1]);
		    }
		    else {
		      badidx->data[highenergyidx] = 1;
		      LogPrintf(LOG_DEBUG,"%s : selected %d as badidx -> energy overlap [%d - %d],[%d - %d]\n",
				__func__,highenergyidx,tempchannellist.channel[k].energy[0],tempchannellist.channel[k].energy[1],
				tempchannellist.channel[s].energy[0],tempchannellist.channel[s].energy[1]);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      LogPrintf(LOG_DEBUG,"%s : checked for energy consistency\n",__func__);
      
      /* now fill in the plan by picking the non-bad indices */
      {
	INT4 sum = badidx->length;
	INT4 newcount = 0;
	REAL8 largestdt = 0.0;
	
	/* compute the total number of usable channels left */
	for (k=0;k<(INT4)badidx->length;k++) sum -= badidx->data[k];
	
	/* allocate memory for the current plan if there was any good data */
	if (sum) {
	  
	  if ((plans->data[pcount].channellist.channel = (FrameChannel *)XLALCalloc(sum,sizeof(FrameChannel))) == NULL) {
	    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FrameCombinationPlan.\n",__func__);
	    XLAL_ERROR(XLAL_ENOMEM);
	  }
	  
	  /* loop over coincident files and pick out the good ones */
	  for (k=0;k<(INT4)tempchannellist.length;k++) {
	    
	    if (badidx->data[k] == 0) {
	      
	      /* find largest sampling time */
	      if (tempchannellist.channel[k].dt > largestdt) largestdt = tempchannellist.channel[k].dt;
	      
	      strncpy(plans->data[pcount].channellist.channel[newcount].channelname,tempchannellist.channel[k].channelname,STRINGLENGTH);
	      strncpy(plans->data[pcount].channellist.channel[newcount].filename,tempchannellist.channel[k].filename,STRINGLENGTH);
	      memcpy(&(plans->data[pcount].channellist.channel[newcount].epoch),&(tempchannellist.channel[k].epoch),sizeof(LIGOTimeGPS));
	      memcpy(&(plans->data[pcount].channellist.channel[newcount].end),&(tempchannellist.channel[k].end),sizeof(LIGOTimeGPS));
	      plans->data[pcount].channellist.channel[newcount].dt = tempchannellist.channel[k].dt;
	      plans->data[pcount].channellist.channel[newcount].duration = tempchannellist.channel[k].duration;
	      plans->data[pcount].channellist.channel[newcount].energy[0] = tempchannellist.channel[k].energy[0];
	      plans->data[pcount].channellist.channel[newcount].energy[1] = tempchannellist.channel[k].energy[1];
	      plans->data[pcount].channellist.channel[newcount].lld = tempchannellist.channel[k].lld;
	      plans->data[pcount].channellist.channel[newcount].type = tempchannellist.channel[k].type;
	      memcpy(&(plans->data[pcount].channellist.channel[newcount].detconfig),&(tempchannellist.channel[k].detconfig),NPCU*sizeof(INT4));
	      
	      LogPrintf(LOG_DEBUG,"%s : plan[%d] -> filename = %s\n",__func__,pcount,plans->data[pcount].channellist.channel[newcount].filename);
	      LogPrintf(LOG_DEBUG,"%s : plan[%d] -> channelname = %s\n",__func__,pcount,plans->data[pcount].channellist.channel[newcount].channelname);
	      LogPrintf(LOG_DEBUG,"%s : plan[%d] -> epoch = %d\n",__func__,pcount,plans->data[pcount].channellist.channel[newcount].epoch.gpsSeconds);
	      LogPrintf(LOG_DEBUG,"%s : plan[%d] -> end = %d\n",__func__,pcount,plans->data[pcount].channellist.channel[newcount].end.gpsSeconds);
	      LogPrintf(LOG_DEBUG,"%s : plan[%d] -> duration = %f\n",__func__,pcount,plans->data[pcount].channellist.channel[newcount].duration);
	      LogPrintf(LOG_DEBUG,"%s : plan[%d] -> dt = %f\n",__func__,pcount,plans->data[pcount].channellist.channel[newcount].dt);
	      LogPrintf(LOG_DEBUG,"%s : plan[%d] -> energy = %d - %d\n",__func__,pcount,plans->data[pcount].channellist.channel[newcount].energy[0],plans->data[pcount].channellist.channel[newcount].energy[1]);
	      LogPrintf(LOG_DEBUG,"%s : plan[%d] -> lld = %d\n",__func__,pcount,plans->data[pcount].channellist.channel[newcount].lld);
	      LogPrintf(LOG_DEBUG,"%s : plan[%d] -> data type = %d\n",__func__,pcount,plans->data[pcount].channellist.channel[newcount].type);
	      newcount++;
	      
	    }
	    
	  }
	  
	  /* fill in extra plan params */
	  plans->data[pcount].dt = largestdt;
	  memcpy(&(plans->data[pcount].epoch),&(epochlist->data[i]),sizeof(LIGOTimeGPS));
	  memcpy(&(plans->data[pcount].end),&(epochlist->data[i+1]),sizeof(LIGOTimeGPS));
	  plans->data[pcount].duration = XLALGPSDiff(&(epochlist->data[i+1]),&(epochlist->data[i]));
	  plans->data[pcount].channellist.length = newcount;
	  plans->data[pcount].channellist.npcus = tempchannellist.npcus;
	  snprintf(plans->data[pcount].OBS_ID,STRINGLENGTH,"%s",tempchannellist.channel[0].OBS_ID);
	  snprintf(plans->data[pcount].channellist.dir,STRINGLENGTH,"%s",tempchannellist.dir);
	  LogPrintf(LOG_DEBUG,"%s : plan[%d] -> common dt = %f.\n",__func__,pcount,plans->data[pcount].dt);
	  LogPrintf(LOG_DEBUG,"%s : plan[%d] -> common epoch = %d\n",__func__,pcount,plans->data[pcount].epoch.gpsSeconds);
	  LogPrintf(LOG_DEBUG,"%s : plan[%d] -> common end = %d\n",__func__,pcount,plans->data[pcount].end.gpsSeconds);
	  LogPrintf(LOG_DEBUG,"%s : plan[%d] -> common duration = %f\n",__func__,pcount,plans->data[pcount].duration);
	  LogPrintf(LOG_DEBUG,"%s : plan[%d] -> common npcus = %d\n",__func__,pcount,plans->data[pcount].channellist.npcus);
	  LogPrintf(LOG_DEBUG,"%s : plan[%d] -> common OBS_ID = %s\n",__func__,pcount,plans->data[pcount].OBS_ID);
	  LogPrintf(LOG_DEBUG,"%s : plan[%d] -> number of channels = %d\n",__func__,pcount,plans->data[pcount].channellist.length);
	  pcount++;
	}
	
      }
      LogPrintf(LOG_DEBUG,"%s : recorded plan.\n",__func__);
      
      /* free the temp channel list */
      XLALFree(tempchannellist.channel);
      XLALDestroyINT4Vector(badidx);
      
    }
    else if ( (span < MINFRAMELENGTH) && (span >= 0) ){
      LogPrintf(LOG_DEBUG,"%s : epoch boundaries span (%f) no data > MINFRAMELENGTH (%d).\n",__func__,span,MINFRAMELENGTH);
    }
    else {
      LogPrintf(LOG_CRITICAL,"%s : epoch boundaries span negative duration (%d - %d) (%f).\n",__func__,epochlist->data[i].gpsSeconds,epochlist->data[i+1].gpsSeconds,span);
      XLAL_ERROR(XLAL_ENOMEM);
    }
    
  }
  
  /* resize plan */
  plans->length = pcount;
  if (plans->length > 0) {
    if ((plans->data = (FrameCombinationPlan *)XLALRealloc(plans->data,plans->length*sizeof(FrameCombinationPlan))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FrameCombinationPlan.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }
  }
  
  /* free mem */
  XLALDestroyTimestampVector(epochlist);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
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

/** this function combines the files listed in the combination plan into a single REAL4 timeseries 
 */
int XLALCombinationPlanToREAL4TimeSeries(REAL4TimeSeries **ts,           /**< [out] the timeseries containing the combined data */
					 HeaderVector *header,          /**< [out] the combined history fields of all files */
					 FrameCombinationPlan *plan,    /**< [in] the plan describing which files to combine */
					 UINT2 npcus                    /**< [in] the number of operational PCUs in this segment */
					 )
{
  INT8 N;
  INT4 i,j,k;
 
  /* check input arguments */
  if ((*ts) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output INT4TimeSeries structure != NULL.\n",__func__);
    XLAL_ERROR(XLAL_EFAULT);
  } 
  if (header == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output history string == NULL.\n",__func__);
    XLAL_ERROR(XLAL_EFAULT);
  } 
  if (plan == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input frame combination plan structure pointer = NULL.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }  
  LogPrintf(LOG_DEBUG,"%s : checked input\n",__func__);
  
  /* determine how many samples are required for the timeseries */
  LogPrintf(LOG_DEBUG,"%s : current plan has duration = %f sec and dt = %f sec.\n",__func__,plan->duration,plan->dt);
  N = floor(0.5 + plan->duration/plan->dt);
  LogPrintf(LOG_DEBUG,"%s : combined timeseries requires %ld samples.\n",__func__,N);
  
  /* create the output timeseries */
  if ( ((*ts) = XLALCreateREAL4TimeSeries("X1:COMBINED_PHOTONFLUX",&(plan->epoch),0,plan->dt,&lalDimensionlessUnit,N)) == NULL) {
    LogPrintf(LOG_CRITICAL, "%s : XLALCreateREAL4TimeSeries() failed to allocate an %d length timeseries with error = %d.\n",__func__,N,xlalErrno);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  LogPrintf(LOG_DEBUG,"%s : allocated memory for temporary timeseries\n",__func__);
  
  /* initialise the output timeseries */
  memset((*ts)->data->data,0,(*ts)->data->length*sizeof(REAL4));
  
  /* allocate mem for the ascii headers */
  header->length = (INT4)plan->channellist.length;
  if ((header->data = (Header *)XLALCalloc(header->length,sizeof(Header))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: unable to allocate memory for the ascii header structures\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  LogPrintf(LOG_DEBUG,"%s: allocated memory for the ascii header structure.\n",__func__);

  /* loop over each channel in the plan */
  for (i=0;i<(INT4)plan->channellist.length;i++) {

    INT4TimeSeries *tempts = NULL;
    LALFrStream *fs = NULL;
    CHAR channelname[STRINGLENGTH];
    INT4 lldfactor = 1;
    INT8 testN = (INT8)(plan->duration/plan->channellist.channel[i].dt);
    REAL4 fluxfactor = 1.0/((REAL4)npcus*(REAL4)PCU_AREA*(REAL4)plan->dt);

    /* define channel name */
    snprintf(channelname,STRINGLENGTH,"%s",plan->channellist.channel[i].channelname);
    LogPrintf(LOG_DEBUG,"%s: channel %d/%d is %s.\n",__func__,i+1,plan->channellist.length,channelname);

    /* open the frame file */
    if ((fs = XLALFrStreamOpen(plan->channellist.dir,plan->channellist.channel[i].filename)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to open frame file %s.\n",__func__,plan->channellist.channel[i].filename);
      XLAL_ERROR(XLAL_EINVAL);
    }
    LogPrintf(LOG_DEBUG,"%s: opened frame file %s.\n",__func__,plan->channellist.channel[i].filename);
    
    /* get history information fom this file */
    if (XLALReadFrameHistory(&(header->data[i].header_string),fs->file)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALReadFrameHistory() unable to read history from frame file %s.\n",__func__,plan->channellist.channel[i].filename);
      XLAL_ERROR(XLAL_EINVAL);
    }
    LogPrintf(LOG_DEBUG,"%s : read history field from file %s.\n",__func__,plan->channellist.channel[i].filename);

    /* seek to the start of the frame */
    if (XLALFrStreamSeek(fs,&(plan->epoch))) {
      LogPrintf(LOG_CRITICAL,"%s: unable to seek to start of frame file %s.\n",__func__,plan->channellist.channel[i].filename);
      XLAL_ERROR(XLAL_EINVAL);
    }

    /* read in timeseries from this file - final arg is limit on length of timeseries (0 = no limit) */
    if ((tempts = XLALFrStreamReadINT4TimeSeries(fs,channelname,&(plan->epoch),plan->duration,testN)) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s: unable to read channel %s from frame file %s.\n",__func__,plan->channellist.channel[i].filename);
      XLAL_ERROR(XLAL_EINVAL);
    }
    LogPrintf(LOG_DEBUG,"%s: reading channel %s\n",__func__,channelname);
    LogPrintf(LOG_DEBUG,"%s: ",__func__);
    for (k=0;k<10;k++) printf("%d ",tempts->data->data[k]);
    printf("\n");
    
    /* close the frame file */
    XLALFrStreamClose(fs);
    LogPrintf(LOG_DEBUG,"%s : closed input frame file.\n",__func__);

    /* if the data is lld data then we multiply by 2 before adding */
    if (plan->channellist.channel[i].lld) lldfactor = 2;

    /* compute the rebinning factor */
    {
      REAL8 realbinratio = (REAL8)((*ts)->deltaT/tempts->deltaT);
      INT4 intbinratio = (INT4)((*ts)->deltaT/tempts->deltaT);
      LogPrintf(LOG_DEBUG,"%s : rebinning factor = %d\n",__func__,intbinratio);

      if ( fmod(realbinratio,1.0) > 0.0 ) {
	LogPrintf(LOG_CRITICAL,"%s: sampling rate is not consistent with an integer factor in file %s.\n",__func__,plan->channellist.channel[i].filename);
	XLAL_ERROR(XLAL_EINVAL);
      }
      
      /* add to output timeseries if the timeseries are consistent */
      if ( XLALGPSCmp(&((*ts)->epoch),&(tempts->epoch)) == 0 ) {
	
	LogPrintf(LOG_DEBUG,"%s : timeseries are consistent so adding to main timeseries\n",__func__);
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
	  (*ts)->data->data[j] += fluxfactor*lldfactor*(REAL4)temp;
	}
      
      }
      else {
	LogPrintf(LOG_DEBUG,"%s : timeseries NOT consistent so not adding to main timeseries\n",__func__);
      }
      
    }
      
    /* free memory */
    XLALDestroyINT4TimeSeries(tempts);
   
  } /* end loop over channels */
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/** this function combines the files listed in the combination plan into a single timeseries 
 */
int XLALREAL4TimeSeriesToFrame(CHAR *outputdir,               /**< [in] name of output directory */
			       REAL4TimeSeries *ts,           /**< [in] timeseries to output */
			       FrameCombinationPlan *plan,    /**< [in] the plan used to generate the timeseries */
			       CHAR *clargs,                  /**< [in] the command line args */
			       HeaderVector *header           /**< [in] the combined history information */
			       )
{
  struct FrameH *outFrame = NULL;        /* frame data structure */
  REAL8 T;
  CHAR filecomment[STRINGLENGTH];
  CHAR outputfile[STRINGLENGTH];
  INT4 i;
  INT8 sum = 0;

  /* define observation time */
  T = ts->deltaT*ts->data->length;
  LogPrintf(LOG_DEBUG,"%s : defined output observation time as %f\n",__func__,T);
 
  /* output initial output samples */
  {
    INT4 k;
    LogPrintf(LOG_DEBUG,"%s: ",__func__);
    for (k=0;k<10;k++) printf("%f ",ts->data->data[k]);
    printf(" ... ");
    for (k=(INT4)ts->data->length-11;k<(INT4)ts->data->length;k++) printf("%f ",ts->data->data[k]);
    printf("\n");
  }

  /* check for an all zero timeseries */
  for (i=0;i<(INT4)ts->data->length;i++) sum += (INT8)ts->data->data[i];

  /* if not full of zeros then proceed */
  if (sum) {
    
    /* generate a frame data structure - last three inputs are [project, run, frnum, detectorFlags] */
    if ((outFrame = XLALFrameNew(&(ts->epoch),T,"XTE_PCA",1,0,0)) == NULL) {
      LogPrintf(LOG_CRITICAL, "%s : XLALFrameNew() failed with error = %d.\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_EFAILED);
    }
    LogPrintf(LOG_DEBUG,"%s : set-up frame structure\n",__func__);
    
    /* add timeseries to frame structure */
    if (XLALFrameAddREAL4TimeSeriesProcData(outFrame,ts)) {
      LogPrintf(LOG_CRITICAL, "%s : XLALFrameAddINT4TimeSeries() failed with error = %d.\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_EFAILED);
    }
    LogPrintf(LOG_DEBUG,"%s : added timeseries to frame structure\n",__func__);
    
    /* Here's where we add extra information into the frame */
    /* we include the current command line args and the original FITS file headers from each contributing file */
    /* we also add the git version info */
    /* we also add the keyord pairs for NPCUS, MINENERGY, MAXENERGY etc.. */
    {
      CHAR *versionstring = NULL;               /* pointer to a string containing the git version information */
      CHAR npcus_string[STRINGLENGTH];
      CHAR OBS_ID_string[STRINGLENGTH];
      CHAR deltat_string[STRINGLENGTH];
      CHAR tobs_string[STRINGLENGTH];
      CHAR nchannels_string[STRINGLENGTH];
      versionstring = XLALGetVersionString(1); 
      FrHistoryAdd(outFrame,versionstring); 
      FrHistoryAdd(outFrame,clargs); 
      for (i=0;i<header->length;i++) FrHistoryAdd(outFrame,header->data[i].header_string);
      XLALFree(versionstring);
      snprintf(OBS_ID_string,STRINGLENGTH,"OBS_ID = %s",plan->OBS_ID);
      snprintf(npcus_string,STRINGLENGTH,"NPCUS = %d",plan->channellist.npcus);
      snprintf(deltat_string,STRINGLENGTH,"DELTAT = %6.12f",ts->deltaT);
      snprintf(tobs_string,STRINGLENGTH,"TOBS = %6.12f",T);
      snprintf(nchannels_string,STRINGLENGTH,"NCHANNELS = %d",plan->channellist.length);

      /* loop over each channel and add some metadata to the history field */
      for (i=(INT4)plan->channellist.length-1;i>=0;i--) {
	CHAR minenergy_string[STRINGLENGTH];
	CHAR maxenergy_string[STRINGLENGTH];
	CHAR dt_string[STRINGLENGTH];
	CHAR filename_string[STRINGLENGTH];
	CHAR channelname_string[STRINGLENGTH];
	CHAR lld_string[STRINGLENGTH];
	CHAR type_string[STRINGLENGTH];
	CHAR detconfig_string[STRINGLENGTH];
	snprintf(detconfig_string,STRINGLENGTH,"DETCONFIG_%d = %d%d%d%d%d",i,plan->channellist.channel[i].detconfig[0],
		 plan->channellist.channel[i].detconfig[1],plan->channellist.channel[i].detconfig[2],
		 plan->channellist.channel[i].detconfig[3],plan->channellist.channel[i].detconfig[4]);
	snprintf(type_string,STRINGLENGTH,"TYPE_%d = %d",i,plan->channellist.channel[i].type);
	snprintf(lld_string,STRINGLENGTH,"LLD_%d = %d",i,plan->channellist.channel[i].lld);
	snprintf(minenergy_string,STRINGLENGTH,"MAXENERGY_%d = %d",i,plan->channellist.channel[i].energy[1]);
	snprintf(maxenergy_string,STRINGLENGTH,"MINENERGY_%d = %d",i,plan->channellist.channel[i].energy[0]);
	snprintf(dt_string,STRINGLENGTH,"DT_%d = %6.12f",i,plan->channellist.channel[i].dt);
	snprintf(channelname_string,STRINGLENGTH,"CHANNELNAME_%d = %s",i,plan->channellist.channel[i].channelname);
	snprintf(filename_string,STRINGLENGTH,"FILENAME_%d = %s",i,plan->channellist.channel[i].filename);
	FrHistoryAdd(outFrame,detconfig_string);
	FrHistoryAdd(outFrame,type_string);
	FrHistoryAdd(outFrame,lld_string);
	FrHistoryAdd(outFrame,maxenergy_string);
	FrHistoryAdd(outFrame,minenergy_string);
	FrHistoryAdd(outFrame,dt_string);
	FrHistoryAdd(outFrame,channelname_string);
	FrHistoryAdd(outFrame,filename_string);
      }
      FrHistoryAdd(outFrame,nchannels_string);
      FrHistoryAdd(outFrame,tobs_string); 
      FrHistoryAdd(outFrame,deltat_string);
      FrHistoryAdd(outFrame,npcus_string); 
      FrHistoryAdd(outFrame,OBS_ID_string); 
    }
      
    /* construct file name - we use the LIGO format <DETECTOR>-<COMMENT>-<GPSSTART>-<DURATION>.gwf */
    /* the comment field we sub-format into <INSTRUMENT>_<FRAME>_<SOURCE>_<OBSID_APID> */
    /* first we need to extract parts of the original filenames */
    {
      CHAR originalfile[STRINGLENGTH];
      CHAR *c1,*c2,*c3;
      INT4 j;
      INT4 n;
      snprintf(originalfile,STRINGLENGTH,"%s",plan->channellist.channel[0].filename);
      if (!(c1 = strrchr(originalfile,'/'))) c1 = originalfile;
      c2 = c1;
      for (j=0;j<3;j++) {
	CHAR *temp = strstr(c1+1,"-");
	c2 = temp;
      }
      c3 = c2;
      for (j=0;j<4;j++) {
	CHAR *temp = strstr(c3+1,"_");
	c3 = temp;
      }
      n = strlen(c2) - strlen(c3);
      snprintf(filecomment,n,"%s",c1+1);
    }
    
    /* generate output filename */
    snprintf(outputfile,STRINGLENGTH,"%s/X1-%s-%d-%d.gwf",outputdir,filecomment,ts->epoch.gpsSeconds,(INT4)T);
    LogPrintf(LOG_DEBUG,"%s : output file = %s\n",__func__,outputfile);
    
    /* write frame structure to file (opens, writes, and closes file) - last argument is compression level */
    if (XLALFrameWrite(outFrame,outputfile,0)) {
      LogPrintf(LOG_CRITICAL, "%s : XLALFrameWrite() failed with error = %d.\n",__func__,xlalErrno);
      XLAL_ERROR(XLAL_EFAILED);
    }
    LogPrintf(LOG_DEBUG,"%s : written frame to output file\n",__func__);

    /* free the frame structure */
    /* there doesn't seem to be an XLAL function for doing this */
    FrameFree(outFrame);

  }
  else {
    LogPrintf(LOG_NORMAL,"%s : timeseries contains only zero counts.  Not generating a frame file.\n",__func__);
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;
  
}

/** this function reads in the frame history as a string
 */
int XLALReadFrameHistory(CHAR **history_string,     /**< [out] the history field read in as a string */ 
			 FrFile *file               /**< [in] frame file pointer */
			 )
{
  FrameH *frame = NULL;
  INT4 stringlen = 1;
  FrHistory *localhist;

  /* check input */
  if ((*history_string) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input history string is not null.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
  }
  if (file == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input frame file pointer is null.\n",__func__);
    XLAL_ERROR(XLAL_EINVAL);
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
      LogPrintf(LOG_CRITICAL,"%s : failed to re-allocate memory for history string.\n",__func__);
      XLAL_ERROR(XLAL_ENOMEM);
    }

    /* append the current history string to the output */
    strncat((*history_string),localhist->comment,n);

    /* point to the next history record */
    localhist = localhist->next;
  }

  /* extend the length of the output to include a new line character */
  if ( ( (*history_string) = (CHAR *)XLALRealloc((*history_string),(stringlen+1)*sizeof(CHAR))) == NULL ) {
    LogPrintf(LOG_CRITICAL,"%s : failed to re-allocate memory for history string.\n",__func__);
    XLAL_ERROR(XLAL_ENOMEM);
  }
  strncat((*history_string),"\n",1);

  LogPrintf(LOG_DEBUG,"%s : length of history string = %d characters .\n",__func__,stringlen);

  /* free the frame */
  FrameFree(frame);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",__func__);
  return XLAL_SUCCESS;

}
