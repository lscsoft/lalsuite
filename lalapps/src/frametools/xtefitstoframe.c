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
 * This code is designed to convert an (R)XTE FITS data file into a frame file. 
 *
 * It reads in either individually time-tagged photons or pre-binned photon counts.  
 * If the input file has been pre-processed such that it also contains barycentric 
 * time information then the user can specify that the output timeseries be 
 * barycentered.  The final timeseries is output in frame file format.
 *
 */

/***********************************************************************************************/
/* includes */
#include <fitsio.h>
#include <math.h>
#include <gsl/gsl_interp.h>        /* needed for the gsl interpolation */
#include <gsl/gsl_spline.h>        /* needed for the gsl interpolation */
#include <gsl/gsl_rng.h>
#include <lal/TimeSeries.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>
#include <lal/UserInput.h>
#include <lal/LogPrintf.h>
#include <lal/LALFrameIO.h>
#include <lal/FrameStream.h>
#include <lalappsfrutils.h>
#include <lalapps.h>

/***********************************************************************************************/
/* some global constants */
#define MAX_DT 0.000976562            /* the maximum sample time we will use (equiv max freq 512 Hz) */
#define MIN_DT 0.000244140625         /* the minimum sample time we will use (equiv max freq 2048 Hz) */
#define GPSMINUSTAI 441417609         /* IMPORTANT, DO NOT TOUCH = difference between GPS and TAI in seconds on January 1, 1994, at 00:00:28 */
#define TIMESTAMPDELTAT 1             /* the spacing used when reducing the number of event time stamps (seconds) */
#define MAXTIMESTAMPDELTAT 256        /* the maximum allowed spacing timestamp spacing for barycentering (seconds) */
#define APIDLENGTH 5                  /* the length of APID HEX strings */
#define STRINGLENGTH 256              /* the length of general string */
#define LONGSTRINGLENGTH 1024         /* the length of general string */
#define NAPID 6                       /* the number of valid APIDs we can currently read */
#define NUSELESSDATAMODE 5            /* the number of useless data modes we know of */
#define MINFRAMELENGTH 1              /* the minimum duration of output frame file in seconds */

/***********************************************************************************************/
/* error codes */
#define XTEFITSTOFRAMEC_ENULL 		1
#define XTEFITSTOFRAMEC_ESYS     	2
#define XTEFITSTOFRAMEC_EINPUT   	3
#define XTEFITSTOFRAMEC_EMEM 		4
#define XTEFITSTOFRAMEC_ENONULL     	5
#define XTEFITSTOFRAMEC_EXLAL   	6
#define XTEFITSTOFRAMEC_MSGENULL 	"Arguments contained an unexpected null pointer"
#define XTEFITSTOFRAMEC_MSGESYS		"System call failed (probably file IO)"
#define XTEFITSTOFRAMEC_MSGEINPUT   	"Invalid input"
#define XTEFITSTOFRAMEC_MSGEMEM   	"Out of memory. Bad."
#define XTEFITSTOFRAMEC_MSGENONULL      "Output pointer is non-NULL"
#define XTEFITSTOFRAMEC_MSGEXLAL	"XLALFunction-call failed"

/***********************************************************************************************/
/* structure for storing the content read from a FITS file */

/** A structure for storing vectors of detector and barycentric frame timestamps 
 *
 * A pre-barycentered FITS file contains an original set of detector frame 
 * timestamps (not always uniformly spaced) and a corresponding set of barycentric
 * timestamps.    
 *
 */
typedef struct {
  REAL8 *dettime;                /**< the detector timestamps */
  REAL8 *barytime;               /**< the barycentered timestamps */
  INT8 length;                   /**< the number of timestamps */ 
  REAL8 dtrow;                  /**< the time steps for the timestamps */
} BarycentricData;

/** The good time interval data read from a FITS file  
 *
 */
typedef struct {
  REAL8 *start;                /**< the start times of good data */
  REAL8 *end;                  /**< the end times of good data */
  INT8 length;                 /**< the number of gti segments */
} GTIData;

/** A structure containing all of the relavent information extracted from the header of 
 * an (R)XTE FITS PCA file.
 *
 */
typedef struct {
  char file[STRINGLENGTH];         /**< the full path of the fits file */
  char filename[STRINGLENGTH];     /**< the name of the fits file */
  int event;                       /**< is the data events (1) or array (0) */
  char objectname[STRINGLENGTH];   /**< the name of the source object */
  char obsid[STRINGLENGTH];        /**< the observation ID */
  char mode[STRINGLENGTH];         /**< stores the data mode */
  char apid[APIDLENGTH];           /**< stores the APID hex string */
  int LLD;                         /**< flag for whether it's LLD2 data */
  long int nbytesperrow;           /**< the number of bytes per row for event data */
  double gpsoffset;                /**< the number to be added to TAI times to get to GPS */
  double ra;                       /**< the source right ascension */
  double dec;                      /**< the source declination */
  int nXeCntcol;                   /**< the number of columns with Cnt or Event data */
  int *XeCntcolidx;                /**< the XeCnt column index */
  char **colname;                  /**< stores the Cnt and Event column names */ 
  long int ncols;                  /**< the total number of columns in the data table */
  long int nrows;                  /**< the number of rows ( = number of timestamps) */
  double deltat;                   /**< the deltat for this time series */
  char *headerdump;                /**< the entire header dumped in ascii */
} FITSHeader;

/** A structure to store TDDES2 DDL data.  
 * 
 * This is a string found in the FITS header that contains information
 * regarding the energy range used in the observation and the timing
 * parameters.
 *
 */
typedef struct {
  INT4 minenergy;                  /**< minimum energy channel (0-255) */
  INT4 maxenergy;                  /**< maximum energy channel (0-255) */
  REAL8 deltat;                    /**< the sampling time */  
  REAL8 offset;                    /**< the time offset */
  INT4 nsamples;                   /**< the number of samples */
} XTETDDESParams;

/** A structure containing all of the relavent information extracted from an (R)XTE FITS PCA file. 
 */
typedef struct {
  UINT4 *data;                     /**< vector of data */
  CHAR *undefined;                 /**< data quality flag */
  INT8 length;                     /**< length of the vector */
  INT8 rowlength;                  /**< the number of data per timestamp */
  INT4 nchannels;                  /**< the number of channels */
  INT8 channelsize;                /**< the length of a channel */
  REAL8 deltat;                    /**< the sampling time */
  INT4 energy[2];                  /**< the energy channel range (0-255) */
} XTEUINT4Vector;

/** A structure containing all of the relavent information extracted from an (R)XTE FITS PCA file. 
 */
typedef struct {
  CHAR *data;                      /**< vector of data */
  CHAR *undefined;                 /**< data quality flag */
  INT8 length;                     /**< length of the vector */
  INT8 nevents;                    /**< the actual number of events */
  INT8 rowlength;                  /**< the number of data per timestamp */
  INT4 nchannels;                  /**< the number of channels */
  INT8 channelsize;                /**< the length of a channel */
  REAL8 deltat;                    /**< the sampling time */
  INT4 energy[2];                  /**< the energy channel range (0-255) */
} XTECHARVector;

/** A structure containing all of the relavent information extracted from an (R)XTE FITS PCA file 
 */
typedef struct {
  FITSHeader *header;              /**< FITS file header information */
  XTEUINT4Vector **array;          /**< array of vectors containing array data */
  XTECHARVector **event;           /**< array of vectors containing event data */
  GTIData *gti;                    /**< good time interval information */
  BarycentricData *stamps;         /**< barycentric timestamps information */   
} FITSData;

/** A structure for storing a timeseries of unsigned integers for XTE data.
 * 
 * This structure contains a distilled version of the information found
 * in the FITS data file.  It stores a continuous timeseries as well as 
 * a data quality vector.
 *
 */
typedef struct {
  CHAR colname[STRINGLENGTH];      /**< name of the column from which the data was read */
  UINT4 *data;                     /**< the data (timeseries of binned photon counts) */
  CHAR *undefined;                 /**< a data quality flag (0=good, 1=bad) */
  INT8 length;                     /**< length of the timeseries */
  REAL8 deltat;                    /**< the time step size in seconds */
  REAL8 tstart;                    /**< the GPS start time */
  REAL8 T;                         /**< the time span in seconds */
  INT4 energy[2];                  /**< the energy channel range (0-255) */
} XTEUINT4TimeSeries;

/** A structure for storing an array of integer timeseries for XTE data.
 *
 * It is designed to be able to store a number of timeseries with identical start, duration and 
 * sampling times.
 */
typedef struct {
  CHAR objectname[STRINGLENGTH];   /**< the source name */
  CHAR obsid[STRINGLENGTH];        /**< the XTE observation ID <proposal>-<target>-<viewing>-<seq no><type> */
  CHAR apid[APIDLENGTH];           /**< the APID */
  CHAR mode[STRINGLENGTH];         /**< the operation mode as a string */
  INT4 bary;                       /**< barycentered flag */
  INT4 lld;                        /**< lld flag */
  XTEUINT4TimeSeries **ts;         /**< pointer to single timeseries vectors */
  INT4 length;                     /**< number of timeseries */
  CHAR *headerdump;                /**< an ascii dump of the original FITS header */
  CHAR *comment;                   /**< a comment field (used to store original clargs) */ 
} XTEUINT4TimeSeriesArray;

/** A structure that stores user input variables 
 */
typedef struct { 
  BOOLEAN help;		            /**< trigger output of help string */
  CHAR *inputfile;                  /**< directory containing input FITS files */
  CHAR *outputdir;                  /**< name of output directory */
  REAL8 deltat;                     /**< the desired sampling time */
  BOOLEAN bary;                     /**< flag for performing barycentering */
  BOOLEAN version;	            /**< output version-info */
} UserInput_t;

/***********************************************************************************************/
/* global variables */
extern int vrbflg;	 	/**< defined in lalapps.c */

/* keywords in FITS file header */
char string_OBJECT[] = "OBJECT";
char string_RA_NOM[] = "RA_NOM";
char string_DEC_NOM[] = "DEC_NOM";
char string_OBS_ID[] = "OBS_ID";
char string_TFIELDS[] = "TFIELDS";
char string_NAXIS1[] = "NAXIS1";
char string_TIMEZERO[] = "TIMEZERO";
char string_DELTAT[] = "DELTAT";
char string_DATAMODE[] = "DATAMODE";
char string_2LLD[] = "2LLD";
char string_HDUCLAS1[] = "HDUCLAS1";
char string_TDDES2[] = "TDDES2";
char string_TIMEDEL[] = "TIMEDEL";
char string_TTYPE1[] = "TTYPE1";

const char *APID[6] = {"FS37","FS3b","FS3f","FS4f","XENO","FS46"};    /* fill in the HEX APID names */

const char xtechannelname[16] = "XTEcounts";

/* a list of useless data modes (specifically for the Sco X-1 analysis) */
const char *USELESSDATAMODE[5] = {"D_1US_0_249_1024_64S_F","D_1US_0_249_128_1S_F","D_1US_0_249_128_1S_2LLD_F","CB","GoodXenon"};

RCSID( "$Id$");		/* FIXME: use git-ID instead to set 'rcsid' */

/***********************************************************************************************/
/* define functions */
int main(int argc,char *argv[]);
void ReadUserVars(LALStatus *status,int argc,char *argv[],UserInput_t *uvar, CHAR *clargs);

/* FITS reading functions */ 
int XLALReadFITSFile(FITSData **data,char *filename);
int XLALReadFITSHeader(FITSHeader *fitsfileheader,fitsfile *fptr);
int XLALReadFITSGTI(GTIData **fitsfilegti, fitsfile *fptr);
int XLALReadFITSArrayData(XTEUINT4Vector **array,fitsfile *fptr,int col);
int XLALReadFITSEventData(XTECHARVector **event,fitsfile *fptr,int col);
int XLALReadFITSTimeStamps(BarycentricData **fitsfilebary,fitsfile *fptr);
int XLALConvertTDDES(XTETDDESParams **params,char *tddes);
int removechar(CHAR *p,CHAR ch);

/* FITS conversion functions */
int XLALArrayDataToXTEUINT4TimeSeries(XTEUINT4TimeSeries **ts,XTEUINT4Vector *event,BarycentricData *stamps,REAL8 dt);
int XLALEventDataToXTEUINT4TimeSeries(XTEUINT4TimeSeries **ts,XTECHARVector *array,BarycentricData *stamps,REAL8 dt);
int XLALEventDataToXTEUINT4TimeSeriesArray(XTEUINT4TimeSeriesArray **ts,FITSData *fits,REAL8 dt);
int XLALArrayDataToXTEUINT4TimeSeriesArray(XTEUINT4TimeSeriesArray **ts,FITSData *fits,REAL8 dt);
int XLALApplyGTIToXTEUINT4TimeSeries(XTEUINT4TimeSeries **ts,GTIData *gti);
int XLALApplyGTIToXTEUINT4TimeSeriesArray(XTEUINT4TimeSeriesArray **ts,GTIData *gti);
int XLALReduceBarycentricData(BarycentricData **stamps);
int XLALBarycenterXTEUINT4TimeSeries(XTEUINT4TimeSeries **ts,BarycentricData *stamps);
int XLALBarycenterXTEUINT4TimeSeriesArray(XTEUINT4TimeSeriesArray **ts,BarycentricData *stamps);
int XLALXTEUINT4TimeSeriesArrayToFrames(XTEUINT4TimeSeriesArray *ts,CHAR *outputdir);
int XLALXTEUINT4TimeSeriesArrayToGTI(GTIData **gti,XTEUINT4TimeSeriesArray *ts);

/* memory handling functions */
int XLALCreateXTEUINT4TimeSeries(XTEUINT4TimeSeries **x,INT8 N);
int XLALCreateXTEUINT4TimeSeriesArray(XTEUINT4TimeSeriesArray **x,INT8 N);
int XLALCreateBarycentricData(BarycentricData **stamps,INT8 N);
int XLALCreateGTIData(GTIData **gti,INT8 N);
int XLALFreeGTIData(GTIData *gti);
int XLALFreeFITSData(FITSData *x);
int XLALFreeXTEUINT4TimeSeries(XTEUINT4TimeSeries *ts);
int XLALFreeXTEUINT4TimeSeriesArray(XTEUINT4TimeSeriesArray *ts);
int XLALReallocXTEUINT4TimeSeries(XTEUINT4TimeSeries **ts,INT8 N);
int XLALFreeBarycentricData(BarycentricData *stamps);

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
  FITSData *fitsdata = NULL;                    /* a FITS data structure */
  XTEUINT4TimeSeriesArray *ts = NULL;           /* a timeseries array structure */ 
  CHAR clargs[LONGSTRINGLENGTH];                /* store the command line args */

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
  /* READ THE DATA */
  /**********************************************************************************/

  /* read in data from the input FITS file */
  if (XLALReadFITSFile(&fitsdata,uvar.inputfile)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadFitsFile() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }   
  LogPrintf(LOG_DEBUG,"%s : Read in fits data from file %s\n",fn,fitsdata->header->filename);  
  
  /* if barycentering was requested check if barycentric data was found */
  if ((uvar.bary == 1) && (fitsdata->stamps->barytime == NULL)) {
    LogPrintf(LOG_CRITICAL,"%s : User requested barycentering but no barycentered timestamps found in file %s.\n",fn,xlalErrno);
    return 1;
  }   

  /**********************************************************************************/
  /* CONVERT THE DATA */
  /**********************************************************************************/

  /* convert fits data to timeseries data */
  if (fitsdata->header->event) {

    if (XLALEventDataToXTEUINT4TimeSeriesArray(&ts,fitsdata,uvar.deltat)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALEventDataToXTEUINT4TimeSeries() failed with error = %d\n",fn,xlalErrno);
      return 1;
    }
    
  }
  else {
    
    if (XLALArrayDataToXTEUINT4TimeSeriesArray(&ts,fitsdata,uvar.deltat)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALArrayDataToXTEUINT4TimeSeries() failed with error = %d\n",fn,xlalErrno);
      return 1;
    }
  }
  LogPrintf(LOG_DEBUG,"%s : converted FITS data structure to a timeseries\n",fn);
    
  /* apply GTI table */
  if (XLALApplyGTIToXTEUINT4TimeSeriesArray(&ts,fitsdata->gti)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALApplyGTITable() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : applied GTI information to data from file %s\n",fn,fitsdata->header->filename);

  if (uvar.bary) {
    
    /* if barycentering then extract timestamps */
    if (XLALReduceBarycentricData(&(fitsdata->stamps))) {
      LogPrintf(LOG_CRITICAL,"%s : XLALReduceBarycentricData() failed with error = %d\n",fn,xlalErrno);
      return 1;
    }
    LogPrintf(LOG_DEBUG,"%s : extracted barycentered timestamps from file %s\n",fn,fitsdata->header->filename);
 
    /* perform barycentering */
    if (XLALBarycenterXTEUINT4TimeSeriesArray(&ts,fitsdata->stamps)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALBarycenterXTEUINT4TimeSeries() failed with error = %d\n",fn,xlalErrno);
      return 1;
    }
    LogPrintf(LOG_DEBUG,"%s : extracted barycentered timestamps from file %s\n",fn,fitsdata->header->filename);

  }
 
  /**********************************************************************************/
  /* OUTPUT THE DATA */
  /**********************************************************************************/

  /* first add command line args to the comment field in the timeseries array */
  if ((ts->comment = (CHAR *)XLALCalloc(LONGSTRINGLENGTH+1,sizeof(CHAR))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for comment field.\n",fn,xlalErrno);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  snprintf(ts->comment,LONGSTRINGLENGTH+1,"\n%s",clargs);

  /* output data */
  if (XLALXTEUINT4TimeSeriesArrayToFrames(ts,uvar.outputdir)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALXTEUINT4TimeSeriesToFrames() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : output data to frame file(s)\n",fn);

  /**********************************************************************************/
  /* CLEAN UP */
  /**********************************************************************************/

  /* free memory */
  if (XLALFreeFITSData(fitsdata)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALFreeFITSData() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }
  if (XLALFreeXTEUINT4TimeSeriesArray(ts)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALFreeXTEUINT4TimeSeries() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }

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
  uvar->inputfile = NULL; 
  uvar->outputdir = NULL; 
  uvar->bary = 0;
  uvar->deltat = MIN_DT;

  /* ---------- register all user-variables ---------- */
  LALregBOOLUserStruct(status, 	help, 		'h', UVAR_HELP,     "Print this message");
  LALregSTRINGUserStruct(status,inputfile, 	'i', UVAR_REQUIRED, "The input FITS file name"); 
  LALregSTRINGUserStruct(status,outputdir, 	'o', UVAR_REQUIRED, "The output frame file directory name"); 
  LALregREALUserStruct(status,  deltat,         't', UVAR_OPTIONAL, "The output sampling time (in seconds)");
  LALregBOOLUserStruct(status, 	bary,   	'b', UVAR_OPTIONAL, "Output barycentered data");
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



/** This function reads in data from a FITS file and returns a FITSdata data structure.
 *
 */
int XLALReadFITSFile(FITSData **fitsfiledata,        /**< [out] FITS file null data structure */
		     CHAR *filepath                  /**< [in] full path to input FITS file */ 
		     )
{
  
  const CHAR *fn = __func__;   /* store function name for log output */
  fitsfile *fptr;              /* pointer to the FITS file, defined in fitsio.h */
  INT4 status = 0;             /* fitsio status flag initialised */
  INT4 i;                      /* counter */
  FITSHeader *header;          /* temporary pointer */

  /* check input arguments */
  if ((*fitsfiledata) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output FITSdata structure != NULL.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }  

  /* allocate memory for output FITS data structure and header */
  if ( ( (*fitsfiledata) = (FITSData *)LALCalloc(1,sizeof(FITSData))) == NULL ) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FITSdata structure.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if ( ( (*fitsfiledata)->header = (FITSHeader *)LALCalloc(1,sizeof(FITSHeader))) == NULL ) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FITSdata structure.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* initialise data pointers */
  (*fitsfiledata)->array = NULL;
  (*fitsfiledata)->event = NULL;
  (*fitsfiledata)->gti = NULL;
  (*fitsfiledata)->stamps = NULL;
  header = (*fitsfiledata)->header;

  /* open the fits file */
  if (fits_open_file(&fptr,filepath,READONLY,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : failed to open FITS file %s for reading.\n",fn,filepath);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* add full file path to the header information */
  strncpy(header->file,filepath,STRINGLENGTH);
  
  /* read the header information from the first extension */
  if (XLALReadFITSHeader(header,fptr)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadFitsHeader() failed to read header information from FITS file %s with error = %d.\n",fn,filepath,xlalErrno);
    XLAL_ERROR(fn,XLAL_EFAULT);
  }
  
  /* read GTI information from the second extension */
  if (XLALReadFITSGTI(&((*fitsfiledata)->gti),fptr)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadFITSGTI() failed to read GTI information from FITS file %s with error = %d.\n",fn,filepath,xlalErrno);
    XLAL_ERROR(fn,XLAL_EFAULT);
  }
  
  /* read the timestamps information */
  if (XLALReadFITSTimeStamps(&((*fitsfiledata)->stamps),fptr)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALReadFITSTimeStamps() failed to read time stamps information from FITS file %s with error = %d.\n",fn,filepath,xlalErrno);
    XLAL_ERROR(fn,XLAL_EFAULT);
  }
  
  /* read in the data */
  if (!(*fitsfiledata)->header->event) {
    
    /* allocate memory for array data pointers */
    if (((*fitsfiledata)->array = (XTEUINT4Vector **)LALCalloc((*fitsfiledata)->header->nXeCntcol,sizeof(XTEUINT4Vector*))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for array data pointers.\n",fn);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
    
    /* loop over XeCnt columns and read in each column */
    for (i=0;i<header->nXeCntcol;i++) {
      
      if (XLALReadFITSArrayData(&((*fitsfiledata)->array[i]),fptr,header->XeCntcolidx[i])) {
	LogPrintf(LOG_CRITICAL,"%s : XLALReadFITSBinnedData() failed to read binned data from FITS file %s with error = %d.\n",fn,filepath,xlalErrno);
	XLAL_ERROR(fn,XLAL_EFAULT);
      }
      
    }
    
  }
  else {
    
    /* allocate memory for event data pointers */
    if (((*fitsfiledata)->event = (XTECHARVector **)LALCalloc(header->nXeCntcol,sizeof(XTECHARVector*))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for event data pointers.\n",fn);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
    
    /* loop over XeCnt columns and read in each column */
    for (i=0;i<header->nXeCntcol;i++) {
      
      if (XLALReadFITSEventData(&((*fitsfiledata)->event[i]),fptr,header->XeCntcolidx[i])) {
	LogPrintf(LOG_CRITICAL,"%s : XLALReadFITSEventData() failed to read event data from FITS file %s with error = %d.\n",fn,filepath,xlalErrno);
	XLAL_ERROR(fn,XLAL_EFAULT);
      }
      
    }
    
  }
  
  /* close the file */
  if (fits_close_file(fptr,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : failed to close FITS file %s.\n",fn,filepath);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;

}

/** Read in GTI (good time interval) table from FITS file.  
 *
 * The second extension HDU2 (3 from the beginning) lists the standard good time intervals, 
 * i.e. the start and stop times of all the potentially useable data in the file. 
 *            
 */
int XLALReadFITSGTI(GTIData **gti,              /**< [out] Good time interval data structure */ 
		    fitsfile *fptr              /**< [in] fitsfile file pointer */
		    )
{
  
  const CHAR *fn = __func__;   /* store function name for log output */
  INT4 i;                      /* counter */
  INT4 hdutype;                /* returned value of the header type */
  REAL8 doublenull = 0.0;      /* dummy variable ? */
  INT4 anynull;                /* dummy ? */
  INT4 status = 0;             /* fitsio status flag initialised */
  INT8 ngtirows;               /* the number of rows in the GTI table */
  REAL8 tzero;                 /* the clock correction */
  REAL8 gpsoffset;             /* the correction needed to convert time values to GPS */
  char comment[80];            /* a temporary comment string */

   /* check input arguments */
  if ((*gti) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output GTIData structure != NULL.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }  
  if (fptr == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input FITS file pointer = NULL.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }  

 /* move to the gti header - the third extension */ 
  if (fits_movabs_hdu(fptr,3,&hdutype,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_movabs_hdu() failed to move to extension header 3.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  
  /* get TIMEZERO for this file - this can be read from any of the extension headers */
  /* Timestamps plus the value of TIMEZERO provide the time in TAI as seconds since January 1, 1994, at 00:00:28 (TAI). */
  if (fits_read_key(fptr,TDOUBLE,string_TIMEZERO,&tzero,comment,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_TIMEZERO);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : read tzero as %6.12f\n",fn,tzero);

  /* convert tzero to GPS - this is simply the TAI time plus the GPS-TAI constant */
  gpsoffset = tzero + (double)GPSMINUSTAI;  
  LogPrintf(LOG_DEBUG,"%s : computed GPS offset as %6.12f\n",fn,gpsoffset);
  
  /* get number of columns in the gti table - there should be 2 */
  {
    int ngticols = 0;
    if (fits_get_num_cols(fptr,&ngticols,&status)) {
      fits_report_error(stderr,status);
      LogPrintf(LOG_CRITICAL,"%s : fits_get_num_cols() failed to return the number of columns.\n",fn);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
    if (ngticols!=2) {
      LogPrintf(LOG_CRITICAL,"%s : The number of GTI columns != 2.\n",fn);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
  }
  
  /* get the number of rows in the gti table - this is the number of good segments in the file */
  if (fits_get_num_rows(fptr,&ngtirows,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_get_num_rows() failed to return the number of rows.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  if (ngtirows<1) {
    LogPrintf(LOG_CRITICAL,"%s : found %ld rows in GTI table, should be >0.\n",fn,ngtirows,fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : found %ld rows in GTI table.\n",fn,ngtirows);
  
  /* allocate memory for GTI table */
  if (XLALCreateGTIData(gti,ngtirows)) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for GTI start time data with error = %d.\n",fn,xlalErrno);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /*  read the gti start and end values  */  
  {
    if (fits_read_col(fptr,TDOUBLE,1,1,1,(*gti)->length,&doublenull,(*gti)->start,&anynull,&status)) {
      fits_report_error(stderr,status);
      LogPrintf(LOG_CRITICAL,"%s : fits_read_col() failed to read GTI table start times.\n",fn);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
    if (fits_read_col(fptr,TDOUBLE,2,1,1,(*gti)->length,&doublenull,(*gti)->end,&anynull,&status)) {
      fits_report_error(stderr,status);
      LogPrintf(LOG_CRITICAL,"%s : fits_read_col() failed to read GTI table end times.\n",fn);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
  }

  /* IMPORTANT : convert all times to GPS using the offset value */
  for (i=0;i<(*gti)->length;i++) {
    (*gti)->start[i] += gpsoffset;
    (*gti)->end[i] += gpsoffset;
    LogPrintf(LOG_DEBUG,"%s : found GTI range %6.12f -> %6.12f.\n",fn,(*gti)->start[i],(*gti)->end[i]);
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;

}

/** Reads in the FITS file first extension header information.
 *
 * The first extension Header HDU1 (2 from the beginning) contains keywords which
 * provide a complete and detailed description of the contents of the first
 * extension. For convenience, it also contains some of the same information as
 * the primary header.
 *
 */
int XLALReadFITSHeader(FITSHeader *header,        /**< [out] The FITS file header structure */
		       fitsfile *fptr             /**< [in] fitsfile file pointer */
		       )
{
  
  const char *fn = __func__;        /* store function name for log output */
  int hdutype;                      /* the return value when moving to headers */
  char comment[80];                 /* a temporary comment string */
  char type[256];                   /* stores the type (array or event) */
  int status = 0;                   /* fitsio status flag initialised */
  int i;                            /* counter */ 

  /* first we extract the filename from the full file path */
  {
    char *c;
    int len = 0;              
    if ((c = strrchr(header->file,'/'))) c++;
    else c = header->file;
    len = strlen(c) + 1;
    strncpy(header->filename,c,len);  
  }
  LogPrintf(LOG_DEBUG,"%s : extracted filename as %s\n",fn,header->filename);

  /* now we extract the APID hex string from the filename */
  snprintf(header->apid,APIDLENGTH,header->filename);
  LogPrintf(LOG_DEBUG,"%s : extracted APID as %s\n",fn,header->apid);

  /* check that APID is one we can currently deal with */
  /* if not then we exit immediately with error code 0 */
  {
    INT4 j = 0;
    for (i=0;i<NAPID;i++) if (!strncmp(header->apid,APID[i],APIDLENGTH-1)) j++;
    if (j==0) {
      LogPrintf(LOG_NORMAL,"%s : Sorry, currently unable to read data with APID = %s.  Exiting.\n",fn,header->apid);
      exit(0);
    }
  }

  /* move to the first extension header */
  if (fits_movabs_hdu(fptr,2,&hdutype,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_movabs_hdu() failed to move to the first extension header.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }

  /* get mode string defining the detector configuration e.g. B_250us_2A_0_17_Q */
  if (fits_read_key(fptr,TSTRING,string_DATAMODE,&(header->mode),comment,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_DATAMODE);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : read data mode as %s\n",fn,header->mode);

  /* check that the data mode is one we can currently deal with */
  /* if not then we exit immediately with error code 0 */
  {
    INT4 j = 0;
    for (i=0;i<NUSELESSDATAMODE;i++) if (!strcasecmp(header->mode,USELESSDATAMODE[i])) j++;
    if (j>0) {
      LogPrintf(LOG_NORMAL,"%s : Sorry, currently unable to read data data mode = %s.  Exiting.\n",fn,header->mode);
      exit(0);
    }
  }

  /* check object name and store in output structure */
  {
    CHAR tempobject[STRINGLENGTH];
    if (fits_read_key(fptr,TSTRING,string_OBJECT,&tempobject,comment,&status)) {
      fits_report_error(stderr,status);
      LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_OBJECT);
      XLAL_ERROR(fn,XLAL_EFAULT);
    }
    /* need to get rid of all non-allowed characters e.g. "-" and "_" */
    removechar(tempobject,'_');
    removechar(tempobject,'-');
    strncpy(header->objectname,tempobject,STRINGLENGTH);
   }
  LogPrintf(LOG_DEBUG,"%s : checked object name is %s.\n",fn,header->objectname);

  /* read in sky position */
  if (fits_read_key(fptr,TDOUBLE,string_RA_NOM,&(header->ra),comment,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_RA_NOM);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  if (fits_read_key(fptr,TDOUBLE,string_DEC_NOM,&(header->dec),comment,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_DEC_NOM);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : read ra = %6.12f dec = %6.12f.\n",fn,header->ra,header->dec);

  /* extract obsid */
  {
    CHAR tempobsid[STRINGLENGTH];
    if (fits_read_key(fptr,TSTRING,string_OBS_ID,&tempobsid,comment,&status)) {
      fits_report_error(stderr,status);
      LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_OBS_ID);
      XLAL_ERROR(fn,XLAL_EFAULT);
    }
    /* need to get rid of all non-allowed characters e.g. "-" and "_" */
    removechar(tempobsid,'_');
    removechar(tempobsid,'-');
    strncpy(header->obsid,tempobsid,STRINGLENGTH);
  }
  LogPrintf(LOG_DEBUG,"%s : read obsid as %s\n",fn,header->obsid);

  /* get the number of rows in the data table */
  if (fits_get_num_rows(fptr,&(header->nrows),&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_get_num_rows() failed to read in number of rows.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : found %ld rows in first extension.\n",fn,header->nrows);

  /* get number of columns in the data table */
  if (fits_read_key(fptr,TINT,string_TFIELDS,&(header->ncols),comment,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_TFIELDS);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : found %ld columns in first extension.\n",fn,header->ncols);

  /* get the type SA or SE */
  /* this tells us whether it is binned (science-array) or event (science-event) data */
  /* if not either mode then we exit immediately with no error */
  if (fits_read_key(fptr,TSTRING,string_HDUCLAS1,&type,comment,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_HDUCLAS1);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  if (strcmp(type,"ARRAY")==0) header->event = 0;
  else if (strcmp(type,"EVENTS")==0) header->event = 1;
  else {
    LogPrintf(LOG_NORMAL,"%s : data type \"%s\" not recognised.  Exiting. \n",fn,type);
    exit(0);
  }
  LogPrintf(LOG_DEBUG,"%s : data type is %s\n",fn,type);

  /* allocate memory for column data */
  if ((header->XeCntcolidx = (INT4 *)LALCalloc(header->ncols,sizeof(INT4))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for column data.\n",fn,xlalErrno);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if ((header->colname = (CHAR **)LALCalloc(header->ncols,sizeof(CHAR *))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for column names.\n",fn,xlalErrno);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  for (i=0;i<header->ncols;i++) {
    if ((header->colname[i] = (CHAR *)LALCalloc(STRINGLENGTH,sizeof(CHAR))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for column names.\n",fn,xlalErrno);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
  }
  LogPrintf(LOG_DEBUG,"%s : allocated memory for column data\n",fn);

  /* record all column numbers and names that have XeCnt data */
  {
    char colnamestring[STRINGLENGTH];
    char keyword[STRINGLENGTH];
    INT4 idx = 0;
    for (i=0;i<header->ncols;i++) {
      snprintf(keyword,STRINGLENGTH,"TTYPE%d",i+1);
      if (fits_read_key(fptr,TSTRING,keyword,&colnamestring,comment,&status)) {
	fits_report_error(stderr,status);
	LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,keyword);
	 XLAL_ERROR(fn,XLAL_EFAULT);
      }
      if ( ((strstr(colnamestring,"Cnt")!=NULL) && (strstr(colnamestring,"MeanCnt")==NULL)) || ((strncmp(colnamestring,"Event",5)==0))) {
	header->XeCntcolidx[idx] = i+1;
	snprintf(header->colname[idx],STRINGLENGTH,colnamestring); 
	idx ++;
	LogPrintf(LOG_DEBUG,"%s : found col %d contains Cnt or Event data.\n",fn,i+1);
      }
    }
    header->nXeCntcol = idx;
    if (!header->nXeCntcol) {
      LogPrintf(LOG_NORMAL,"%s : failed to find XeCnt or Event columns in file.\n",fn,xlalErrno);
      exit(0);
    }
  }
  
  /* resize memory for column data */
  if ((header->XeCntcolidx = (INT4 *)LALRealloc(header->XeCntcolidx,header->nXeCntcol*sizeof(INT4))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to re-allocate memory for column data.\n",fn,xlalErrno);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  for (i=header->nXeCntcol;i<header->ncols;i++) XLALFree(header->colname[i]);
  if ((header->colname = (CHAR **)LALRealloc(header->colname,header->nXeCntcol*sizeof(CHAR *))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to re-allocate memory for column names.\n",fn,xlalErrno);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  LogPrintf(LOG_DEBUG,"%s : re-allocated memory for column data\n",fn);
  
  /* get number of bytes per row in the data table */
  if (fits_read_key(fptr,TLONG,string_NAXIS1,&(header->nbytesperrow),comment,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_NAXIS1);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : read the number of bytes per row as %ld\n",fn,header->nbytesperrow);

  /* if the data mode is equal to 2LLD then we record this in the LLD (Lower Level Discriminator) flag */
  if (strstr(header->mode,string_2LLD)!=NULL) header->LLD = 1;
  else header->LLD = 0;
  LogPrintf(LOG_DEBUG,"%s : set LLD flag = %d\n",fn,header->LLD);

  /* get sampling time - this is the actual sampling time interval in seconds */
  if (fits_read_key(fptr,TDOUBLE,string_TIMEDEL,&(header->deltat),comment,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_TIMEDEL);
    XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : read TIMEMDEL keyword as %6.12e\n",fn,header->deltat);

  /* dump entire header into a long string */
  /* we do extra stuff to add newline characters at the end of each 80 character long stretch */
  /* it makes it look nicer when you do a frame dump */
  {
    INT4 nkeys;
    CHAR *dummy = NULL;
    if (fits_hdr2str(fptr,0,NULL,0,&dummy,&nkeys,&status)) {
      fits_report_error(stderr,status);
      LogPrintf(LOG_CRITICAL,"%s : fits_hdr2str() failed to read header.\n",fn);
      XLAL_ERROR(fn,XLAL_EFAULT);
    }
    
    /* allocate memory for new char vector */
    if ((header->headerdump = (CHAR *)XLALCalloc(1+nkeys*81,sizeof(CHAR))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to re-allocate memory for column names.\n",fn,xlalErrno);
      XLAL_ERROR(fn,XLAL_ENOMEM);
    }
    
    /* add newline characters after every 80 CHARS */
    for (i=0;i<nkeys;i++) {
      
      INT4 idx = i*80;
      CHAR tempstr[81];
      CHAR kwpair[80];
      snprintf(kwpair,80,"%s",dummy+idx);
      snprintf(tempstr,81,"\n%s",kwpair);
      strncat(header->headerdump,tempstr,81);
     
    }
    free(dummy);
  }
  fprintf(stdout,"%s\n",header->headerdump);
 
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;

}
/** Removes a character from a string 
 *
 */
int removechar(CHAR *p,    /**< [in/out] string to edit */ 
	       CHAR ch     /**< [in] character to remove */
	       )
{
  const char *fn = __func__;        /* store function name for log output */
  CHAR *temp = p;
  
  while (*temp) {
    
    if ((*temp) == (CHAR)ch) {
      
      while (*temp) {
	(*temp) = *(temp+1);	
	if (*temp) temp++;	
      }	
      temp = p;
      
    }    
    temp++;
   
  }
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;

} 

/** Reads in FITS file first extension header timestamps information.
 *
 * It also reads in barycentric timestamps if they are present.
 *
 */
int XLALReadFITSTimeStamps(BarycentricData **stamps,     /**< [out] the detector and barycenter timestamps */
			   fitsfile *fptr              /**< [in] a FITS file pointer */
			   ) 
{
  
  const char *fn = __func__;        /* store function name for log output */
  long int j;
  char comment[80];
  double doublenull = 0.0; 
  int anynull;
  int hdutype;
  int status = 0;                   /* fitsio status flag initialised */
  REAL8 tzero;
  REAL8 gpsoffset;
  INT4 ncols;
  INT4 detidx = 0;
  INT4 baryidx = 0;
  
  /* check input arguments */
  if ((*stamps) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output BarycentricData structure != NULL.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  } 
  if (fptr == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input FITS file pointer = NULL.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }  
  
  /* allocate memory for timestamps structure */
  if (((*stamps) = (BarycentricData *)LALCalloc(1,sizeof(BarycentricData))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for timestamps structure.\n",fn);
       XLAL_ERROR(fn,XLAL_ENOMEM);
    }

  /* move to the first extension header */
  if (fits_movabs_hdu(fptr,2,&hdutype,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_movabs_hdu() failed to move to the first extension header.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  
  /* get DELTAT for this file */
  /* DELTAT gives the the time between each row in the file - these are of order second in length */
  if (fits_read_key(fptr,TDOUBLE,string_DELTAT,&((*stamps)->dtrow),comment,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_DELTAT);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : read DELTAT keyword as  %6.12f\n",fn,(*stamps)->dtrow);

  /* get TIMEZERO for this file */
  /* Timestamps plus the value of TIMEZERO provide the time in TAI as seconds since January 1, 1994, at 00:00:28 (TAI). */
  if (fits_read_key(fptr,TDOUBLE,string_TIMEZERO,&tzero,comment,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_TIMEZERO);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : read tzero as %6.12f\n",fn,tzero);

  /* convert tzero to GPS - this is simply the TAI time plus the GPS-TAI constant */
  gpsoffset = tzero + (double)GPSMINUSTAI;  
  LogPrintf(LOG_DEBUG,"%s : converted tzero to gps -> %6.12f\n",fn,gpsoffset);

  /* get the number of rows in the data table */
  if (fits_get_num_rows(fptr,&((*stamps)->length),&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_get_num_rows() failed to read in number of rows.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : found %ld rows in first extension.\n",fn,(*stamps)->length);
  
  /* get number of columns in the data table */
  if (fits_read_key(fptr,TINT,string_TFIELDS,&ncols,comment,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_TFIELDS);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : found %ld fields in first extension.\n",fn,ncols);

  /* check if there are columns called time and barytime */
  {
    char timestring[STRINGLENGTH];
    char keyword[STRINGLENGTH];
    INT4 i;
    for (i=0;i<ncols;i++) {
      snprintf(keyword,STRINGLENGTH,"TTYPE%d",i+1);
      if (fits_read_key(fptr,TSTRING,keyword,&timestring,comment,&status)) {
	fits_report_error(stderr,status);
	LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,keyword);
	 XLAL_ERROR(fn,XLAL_EFAULT);
      }
      if (strncmp(timestring,"Time",STRINGLENGTH)==0) detidx = i+1;
      if (strncmp(timestring,"BARYTIME",8)==0) baryidx = i+1;
    }
    
  }

  /* if we have a time column then allocate memory and read in the timestamps */
  if (detidx) {

    /* allocate memory for detector frame timestamps */
    if (((*stamps)->dettime = (double *)LALCalloc((*stamps)->length,sizeof(double))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for detector frame timestamps.\n",fn);
       XLAL_ERROR(fn,XLAL_ENOMEM);
    }
 
    /*  read the time detector frame timestamps from column 1 */  
    if (fits_read_col(fptr,TDOUBLE,detidx,1,1,(*stamps)->length,&doublenull,(*stamps)->dettime,&anynull,&status)) {
      fits_report_error(stderr,status);
      LogPrintf(LOG_CRITICAL,"%s : fits_read_col() failed to read in detector frame timestamps.\n",fn);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
    
    /* IMPORTANT : convert to GPS using the offset value */
    for (j=0;j<(*stamps)->length;j++) (*stamps)->dettime[j] += gpsoffset;
    LogPrintf(LOG_DEBUG,"%s : read first detector frame timestamp as %6.12f\n",fn,(*stamps)->dettime[0]);
  
  }
  else {
    LogPrintf(LOG_CRITICAL,"%s : failed to find \"Time\" column.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }

  /* if the file contains barycentered times */
  if (baryidx) {
    
    /* allocate memory for barycentric frame timestamps */
    if (((*stamps)->barytime = (double *)LALCalloc((*stamps)->length,sizeof(double))) == NULL) {
      LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for barycentric frame timestamps.\n",fn);
       XLAL_ERROR(fn,XLAL_ENOMEM);
    }
    
    /*  read the barycentered time values from final column */  
    if (fits_read_col(fptr,TDOUBLE,baryidx,1,1,(*stamps)->length,&doublenull,(*stamps)->barytime,&anynull,&status)) {
      fits_report_error(stderr,status);
      LogPrintf(LOG_CRITICAL,"%s : fits_read_col() failed to read in barycentric frame timestamps.\n",fn);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
    
    /* IMPORTANT : convert to GPS using the offset value */
    for (j=0;j<(*stamps)->length;j++) (*stamps)->barytime[j] += gpsoffset;
    LogPrintf(LOG_DEBUG,"%s : read first barycentered time value as %6.12f\n",fn,(*stamps)->barytime[0]);
  
  }
  else (*stamps)->barytime = NULL;

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;

}

/** Reads in array mode data from a FITS file.
 *
 * Science array mode data is binned data and so contains photon counts.
 *
 */
int XLALReadFITSArrayData(XTEUINT4Vector **array,      /**< [out] the output data vector */
			  fitsfile *fptr,             /**< [in] a fitsfile file pointer */
			  INT4 col                     /**< [in] the column to be read */ 
			  )
{
  
  const char *fn = __func__;        /* store function name for log output */
  INT4 status = 0;                  /* fitsio status flag initialised */
  char comment[80];                 /* a temporary comment string */
  INT4 anynull;
  INT4 hdutype;
  INT4 ncols = 0;
  INT8 nrows = 0;

  /* check input arguments */
  if ((*array) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output XTEUINT4Vector structure != NULL.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  } 
  if (fptr == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input FITS file pointer = NULL.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }  
  if (col < 1) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, column number must be > 0.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }  

  /* move to the first extension header */
  if (fits_movabs_hdu(fptr,2,&hdutype,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_movabs_hdu() failed to move to the first extension header.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  
  /* get number of columns in the data table */
  if (fits_get_num_cols(fptr,&ncols,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_get_num_cols() failed to return the number of columns.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  
  /* get the number of rows in the data table */
  if (fits_get_num_rows(fptr,&nrows,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_get_num_rows() failed to return the number of rows.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  
  /* allocate memory for array data */
  if (((*array) = (XTEUINT4Vector *)LALCalloc(1,sizeof(XTEUINT4Vector))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for array.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* read the number of and size of the dimensions of the requested column */
  if (fits_read_tdim(fptr,col,1,&((*array)->nchannels),&((*array)->channelsize),&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_tdim() failed to read in the number and size of dimensions in col %d.\n",fn,col);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
 
  /* compute the number of samples per row */
  LogPrintf(LOG_DEBUG,"%s : read channelsize as %ld for col %d\n",fn,(*array)->channelsize,col);
  LogPrintf(LOG_DEBUG,"%s : read nchannels as %d for col %d\n",fn,(*array)->nchannels,col);
  if ((*array)->nchannels > 1) {
    printf("FOUND CHANNELS = %d.  What do we do Chris ?\n",(*array)->nchannels);
    exit(1);
  }
    

  /* get the TDDES<col> keyword - describes the data in the second column */
  /* this tells us the sampling time and the energy channels used */
  /* A typical string of DDL (Data Description Language) is a concatenation of tokens */
  /* (denoted by single letters) with assigned values (enclosed in square brackets). */
  {
    char *tddes_string;
    XTETDDESParams *tddes = NULL;
    CHAR keyword[STRINGLENGTH];
    snprintf(keyword,STRINGLENGTH,"TDDES%d",col);
    if (fits_read_key_longstr(fptr,keyword,&tddes_string,comment,&status)) {
      fits_report_error(stderr,status);
      LogPrintf(LOG_CRITICAL,"%s : fits_read_key_longstr() failed to read in long string keyword %s.\n",fn,keyword);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
    LogPrintf(LOG_DEBUG,"%s : read TDDES as %s\n",fn,tddes_string);
  
    /* extract energy and sampling time info from the TDDES2 DDL string */
    if (XLALConvertTDDES(&tddes,tddes_string)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALConvertTDDES() failed to convert TDDES string %s with error = %s.\n",fn,tddes_string,xlalErrno);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
    
    LogPrintf(LOG_DEBUG,"%s : energy range extracted as %d - %d\n",fn,tddes->minenergy,tddes->maxenergy);
    if (tddes->deltat) {
      (*array)->deltat = tddes->deltat;
      LogPrintf(LOG_DEBUG,"%s : read time resolution as %6.12e\n",fn,tddes->deltat);
    }
    else {
      LogPrintf(LOG_DEBUG,"%s : could not extract time resolution from TDDES2 DDL string\n",fn);
      
      /* get sampling time from header field TIMEDEL if cannot get it from tddes string */
      if (fits_read_key(fptr,TDOUBLE,string_TIMEDEL,&((*array)->deltat),comment,&status)) {
	fits_report_error(stderr,status);
	LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_TIMEDEL);
	 XLAL_ERROR(fn,XLAL_EFAULT);
      }
      LogPrintf(LOG_DEBUG,"%s : using TIMEDEL keyword value as array dt (%6.12e)\n",fn,(*array)->deltat);
    
    }
    
    /* extract row length from tddes string - if not available then use header information */
    if (tddes->nsamples) {
      (*array)->rowlength = tddes->nsamples;
      LogPrintf(LOG_DEBUG,"%s : read rowlength from tddes as %d\n",fn,(*array)->rowlength);
    }
    else {
      LogPrintf(LOG_DEBUG,"%s : could not extract row length from TDDES DDL string\n",fn);
      (*array)->rowlength = (*array)->channelsize*(*array)->nchannels;
    }
   
    /* store energy information */
    (*array)->energy[0] = tddes->minenergy;
    (*array)->energy[1] = tddes->maxenergy;

    /* free string mem */
    free(tddes_string);
    XLALFree(tddes);

  }
  
  /* define number of elements to read in - this is the total number of expected data values */
  (*array)->length = (INT8)((*array)->rowlength*nrows);
  LogPrintf(LOG_DEBUG,"%s : computed total number of expected data samples as %ld\n",fn,(*array)->length);
  
  /* check if requested column value is valid */
  if (col>ncols) {
    LogPrintf(LOG_CRITICAL,"%s : tried to read more columns than present in file.\n",fn);
    XLAL_ERROR(fn,XLAL_EFAULT);
  }
  
  /* allocate mem for data and data undefined flag */
  if (((*array)->data = (UINT4 *)LALCalloc((*array)->length,sizeof(UINT4))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for array data.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*array)->undefined = (CHAR *)LALCalloc((*array)->length,sizeof(CHAR))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for data quality flag.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  LogPrintf(LOG_DEBUG,"%s : allocated memory for data (approx %.1f MB)\n",fn,(*array)->length*2e-6);
  
  /*  read the complete column data set */  
  if (fits_read_colnull(fptr,TINT,col,1,1,(*array)->length,(*array)->data,(*array)->undefined,&anynull,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_colnull() failed to read col %d in science array data.\n",fn,col);
    XLAL_ERROR(fn,XLAL_EFAULT);
  }

  /* output debugging information */
  {
    int i;
    char temp1[STRINGLENGTH],temp2[STRINGLENGTH];
    int M = (*array)->length > 10 ? 10 : (*array)->length;
    sprintf(temp1,"");
    sprintf(temp2,"");
    for (i=0;i<M;i++) {
      sprintf(temp2,"%s%d(%d),",temp1,(*array)->data[i],(*array)->undefined[i]);
      strcpy(temp1,temp2);
    }
    LogPrintf(LOG_DEBUG,"%s : read array data as : %s ...\n",fn,temp2);
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}

/** Reads in event mode data from a FITS file.
 *
 * The information regarding each event is stored in groups of nbytes.  The first char in each 
 * group defined whether the event is real.  We read in all information (real and non-real events). 
 *
 */
int XLALReadFITSEventData(XTECHARVector **event,      /**< [out] The FITSdata structure */
  			  fitsfile *fptr,             /**< [in] a fitsfile file pointer */
			  INT4 col                     /**< [in] the column to be read */ 
			  )
{
 
  const CHAR *fn = __func__;        /* store function name for log output */
  INT4 hdutype;
  INT4 status = 0;                   /* fitsio status flag initialised */
  INT4 ncols = 0;
  INT8 nrows = 0;
  char comment[80];                 /* a temporary comment string */

  /* check input arguments */
  if ((*event) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, output XTECHARVector structure != NULL.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  } 
  if (fptr == NULL) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, input FITS file pointer = NULL.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }  
  if (col < 1) {
    LogPrintf(LOG_CRITICAL,"%s: Invalid input, column number must be > 0.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }  

  /* move to the first extension header */
  if (fits_movabs_hdu(fptr,2,&hdutype,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_movabs_hdu() failed to move to the first extension header.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }

  /* get number of columns in the data table */
  if (fits_get_num_cols(fptr,&ncols,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_get_num_cols() failed to return the number of columns.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  
  /* get the number of rows in the data table */
  if (fits_get_num_rows(fptr,&nrows,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_get_num_rows() failed to return the number of rows.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }

  /* check if requested column value is valid */
  if (col>ncols) {
    LogPrintf(LOG_CRITICAL,"%s : tried to read more columns than present in file.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }

  /* allocate mem for data and data undefined flag */
  if (((*event) = (XTECHARVector *)LALCalloc(1,sizeof(XTECHARVector))) == NULL ) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for event structure.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  LogPrintf(LOG_DEBUG,"%s : allocated memory for event structure.\n",fn);
  
  /* read the number of and size of the dimensions of the requested column */
  if (fits_read_tdim(fptr,col,1,&((*event)->nchannels),&((*event)->channelsize),&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_tdim() failed to read in the number and size of dimensions in col %d.\n",fn,col);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }

  LogPrintf(LOG_DEBUG,"%s : read channelsize as %ld for col %d\n",fn,(*event)->channelsize,col);
  LogPrintf(LOG_DEBUG,"%s : read nchannels as %d for col %d\n",fn,(*event)->nchannels,col);
  
  /* get the TDDES<col> keyword - describes the data in the second column */
  /* this tells us the sampling time and the energy channels used */
  /* A typical string of DDL (Data Description Language) is a concatenation of tokens */
  /* (denoted by single letters) with assigned values (enclosed in square brackets). */
  {
    char *tddes_string;
    XTETDDESParams *tddes = NULL;
    CHAR keyword[STRINGLENGTH];
    snprintf(keyword,STRINGLENGTH,"TDDES%d",col);
    if (fits_read_key_longstr(fptr,keyword,&tddes_string,comment,&status)) {
      fits_report_error(stderr,status);
      LogPrintf(LOG_CRITICAL,"%s : fits_read_key_longstr() failed to read in long string keyword %s.\n",fn,keyword);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
    LogPrintf(LOG_DEBUG,"%s : read TDDES2 as %s\n",fn,tddes_string);
  
    /* extract energy and sampling time info from the TDDES2 DDL string */
    if (XLALConvertTDDES(&tddes,tddes_string)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALConvertTDDES() failed to convert TDDES string %s with error = %s.\n",fn,tddes_string,xlalErrno);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
    
    LogPrintf(LOG_DEBUG,"%s : energy range extracted as %d - %d\n",fn,tddes->minenergy,tddes->maxenergy);
    if (tddes->deltat) {
      (*event)->deltat = tddes->deltat;
      LogPrintf(LOG_DEBUG,"%s : read time resolution as %6.12e\n",fn,tddes->deltat);
    }
    else {
      LogPrintf(LOG_DEBUG,"%s : could not extract time resolution from TDDES2 DDL string\n",fn);
      
      /* get sampling time from header field TIMEDEL if cannot get it from tddes string */
      if (fits_read_key(fptr,TDOUBLE,string_TIMEDEL,&((*event)->deltat),comment,&status)) {
	fits_report_error(stderr,status);
	LogPrintf(LOG_CRITICAL,"%s : fits_read_key() failed to read in keyword %s.\n",fn,string_TIMEDEL);
	 XLAL_ERROR(fn,XLAL_EFAULT);
      }
      LogPrintf(LOG_DEBUG,"%s : using TIMEDEL keyword value as event dt (%6.12e)\n",fn,(*event)->deltat);
    
    }

    /* extract row length from tddes string - if not available then use header information */
    if (tddes->nsamples) {
      (*event)->rowlength = tddes->nsamples;
      LogPrintf(LOG_DEBUG,"%s : read rowlength from tddes as %d\n",fn,(*event)->rowlength);
    }
    else {
      LogPrintf(LOG_DEBUG,"%s : could not extract row length from TDDES DDL string\n",fn);
      (*event)->rowlength = (*event)->channelsize*(*event)->nchannels;
    }
    
    /* store energy information */
    (*event)->energy[0] = tddes->minenergy;
    (*event)->energy[1] = tddes->maxenergy;

    /* free string mem */
    free(tddes_string);
    XLALFree(tddes);
    
  }
  
  /* define number of elements to read in (each of size char) - this is the total number of expected data values */
  (*event)->nevents = nrows;
  (*event)->length = (INT8)((*event)->rowlength*nrows*(*event)->nchannels);
  LogPrintf(LOG_DEBUG,"%s : preparing to read in %ld events\n",fn,(*event)->nevents);
  
  /* allocate mem for data and data undefined flag */
  if (((*event)->data = (CHAR *)LALCalloc((*event)->length,sizeof(CHAR))) == NULL ) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for event data.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  LogPrintf(LOG_DEBUG,"%s : allocated memory for nbytes*sizeof(char) = %ld\n",fn,(*event)->nevents);
  
  /* read the complete data set - we must remember that any event with the most significant bit = 0 is not a real event !! */
  /* IMPORTANT : we read in rowlength chars for each event and the first char in each row defines whether the event is real */ 
  if (fits_read_col_bit(fptr,col,1,1,(*event)->length,(*event)->data,&status)) {
    fits_report_error(stderr,status);
    LogPrintf(LOG_CRITICAL,"%s : fits_read_col_bit() failed to read in event data.\n",fn);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : read in %ld events\n",fn,(*event)->nevents);
  
  /* output debugging information */
  {
    int i;
    char temp1[STRINGLENGTH],temp2[STRINGLENGTH];
    int M = (*event)->nevents > 10 ? 10 : (*event)->nevents;
    sprintf(temp1,"");
    sprintf(temp2,"");
    for (i=0;i<M;i++) {
      sprintf(temp2,"%s%d,",temp1,(*event)->data[(*event)->rowlength*i]);
      strcpy(temp1,temp2);
    }
    LogPrintf(LOG_DEBUG,"%s : read array data as : %s ...\n",fn,temp2);
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
 
}



/** Extracts energy range and sampling parameters from tddes2 DDL string.
 *
 */
int XLALConvertTDDES(XTETDDESParams **params,     /**< [out] a null TDDES parameter structure */
		     char *tddes                 /**< [in] a TDDES DDL string */
		     )
{
  
  const CHAR *fn = __func__;         /* store function name for log output */
  
  char out1[256], out2[256], new[256];
  char out3[256], out4[256], out5[256];
  char *tempstr1,*tempstr2,*tempstr3;
  int N1, N2, N3;
 
  /* allocate memory for the output params */
  if (((*params) = (XTETDDESParams *)XLALCalloc(1,sizeof(XTETDDESParams))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for an XTETDDESParams structure.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  
  /* initialise energy ranges */
  (*params)->minenergy = 0;
  (*params)->maxenergy = 0;
  (*params)->deltat = 0.0;
  (*params)->nsamples = 0;
  (*params)->offset = 0.0;
  
  /* find start of energy section of string - exit if not found */
  if ((tempstr1 = strstr(tddes,"C["))==NULL)  return XLAL_SUCCESS;

  /* find location of delimiter following first energy */
  N1 = strcspn(tempstr1,"~");
  
  /* copy first energy to new string */
  strncpy (out1,tempstr1+2,N1-2);

  /* find location of end of energy string */
  N2 = strcspn(tempstr1,"]");

  /* copy whole energy string */
  strncpy(new,tempstr1,N2);
  
  /* find location of last delimiter */
  tempstr2 = strrchr(new,'~');
  sprintf(out2,"%s\n",tempstr2+1);

  /* convert strings to numerical output */
  (*params)->minenergy = atoi(out1);
  (*params)->maxenergy = atoi(out2);

  /* find start of time section of string - exit if not found */
  if ((tempstr3 = strstr(tddes,"T["))==NULL)  return XLAL_SUCCESS;

  /* find the number of characters following "T" to the delimiter following the offset time */
  N1 = strcspn(tempstr3,";"); 

  /* find the number of characters following "T" to the delimiter following dt */
  N2 = strcspn(tempstr3+N1+1,";") + N1 + 1;

  /* find the number of characters following "T" to the delimiter following the number of samples per accumulation */
  N3 = strcspn(tempstr3+N2+1,"]") + N2 + 1;

  /* copy time offset to new string */
  strncpy (out3,tempstr3+2,N1-2);  
  
  /* copy sampling time to new string */
  strncpy (out4,tempstr3+1+N1,N2-N1-1);
 
  /* copy sample number to new string */
  strncpy (out5,tempstr3+N2+1,N3-N2-1); 

  /* convert strings to numerical output */
  (*params)->offset = atof(out3); 
  (*params)->deltat = atof(out4);
  (*params)->nsamples = atoi(out5); 

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;

}

/** Converts a FITS event data file to a binned timeseries.
 *
 * This function acts as a wrapper for XLALEventDataToXTEUINT4TimeSeries which deals 
 * with single timeseries.  We enforce that each timeseries has the same start, end, 
 * and sampling time for consistency.
 * 
 */
int XLALEventDataToXTEUINT4TimeSeriesArray(XTEUINT4TimeSeriesArray **ts,   /**< [out] a NULL timeseries */
					   FITSData *fits,                 /**< in] a FITSHeader structure */
					   REAL8 dt                        /**< [in] the desired sampling time */
					   )
{

  const CHAR *fn = __func__;         /* store function name for log output */
  INT8 i;                            /* counter */
  REAL8 newdt = dt;                  /* the modified sampling time */

  /* check input and output pointers */
  if ((*ts) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries structure has non-null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (fits == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input FITSData structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (fits->header == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input FITSHeader structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* check timeseries consistency */
  for (i=0;i<fits->header->nXeCntcol;i++) {
    if ((fits->event[i]->deltat != fits->event[0]->deltat) || (fits->event[i]->length != fits->event[0]->length)) {
      LogPrintf(LOG_CRITICAL,"%s : inconsistent parameters between column timeseries.\n",fn,xlalErrno);
       XLAL_ERROR(fn,XLAL_EINVAL);
    }
  }

  /* check that requested sampling rate is at least as large as the original rate */
  {
    INT4 n = (INT4)ceil(fits->event[0]->deltat/dt);
    newdt = n*dt;
    LogPrintf(LOG_DEBUG,"%s : requested sampling time %6.12f -> %6.12f sec.\n",fn,dt,newdt);
  }
  
  /* allocate mem for output timeseries array */
  if (((*ts) = (XTEUINT4TimeSeriesArray *)LALCalloc(1,sizeof(XTEUINT4TimeSeriesArray))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for XTEUINT4TimeSeriesArray.\n",fn,xlalErrno);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*ts)->ts = (XTEUINT4TimeSeries **)LALCalloc(fits->header->nXeCntcol,sizeof(XTEUINT4TimeSeries *))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for XTEUINT4TimeSeries pointers.\n",fn,xlalErrno);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* loop over each relevant FITS data column */
  for (i=0;i<fits->header->nXeCntcol;i++) {
   
    if (XLALEventDataToXTEUINT4TimeSeries(&((*ts)->ts[i]),fits->event[i],fits->stamps,newdt)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALEventDataToXTEUINT4TimeSeries() failed with error = %d\n",fn,xlalErrno);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
    LogPrintf(LOG_DEBUG,"%s : converted column (%d/%d) of FITS event data to a timeseries.\n",fn,i+1,fits->header->nXeCntcol);
  
    /* add column name info */
    snprintf((*ts)->ts[i]->colname,STRINGLENGTH,fits->header->colname[i]);

  }
  
  /* update header info */
  strncpy((*ts)->objectname,fits->header->objectname,STRINGLENGTH);
  strncpy((*ts)->obsid,fits->header->obsid,STRINGLENGTH);
  strncpy((*ts)->apid,fits->header->apid,APIDLENGTH);
  (*ts)->bary = 0;
  (*ts)->lld = fits->header->LLD;
  (*ts)->length = fits->header->nXeCntcol;

  /* allocate mem for the header dump and copy the original fits header */
  if (((*ts)->headerdump = (CHAR *)XLALCalloc(1,sizeof(fits->header->headerdump))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FITS header dump copy.\n",fn,xlalErrno);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  memcpy((*ts)->headerdump,fits->header->headerdump,sizeof(fits->header->headerdump));
 
  

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
  
}

/** Converts a FITS event data file to a binned timeseries.
 *
 * Using the detector timestamps data we convert the FITS data into a continuous binned timeseries. 
 *
 */
int XLALEventDataToXTEUINT4TimeSeries(XTEUINT4TimeSeries **ts,      /**< [out] a NULL timeseries */
				      XTECHARVector *event,        /**< [in] a FITSdata structure */
				      BarycentricData *stamps,     /**< [in] timestamps data */
				      REAL8 dt                     /**< [in] the sampling time */
				      )
{

  const CHAR *fn = __func__;         /* store function name for log output */
  INT8 i;                            /* counter */

  /* check input and output pointers */
  if ((*ts) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries structure has non-null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (event == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input XTECHARVector structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (stamps == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input BarycentricData structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* compute number of samples in the continuous timeseries */
  /* define the time span in the detector frame from start of first bin to end of last bin */
  {
    REAL8 tstart = stamps->dettime[0];
    REAL8 tempspan = stamps->dettime[stamps->length-1] - tstart + stamps->dtrow;
    INT8 N = (INT8)(tempspan/dt);
    
    /* allocate mem for output timeseries */
    if (XLALCreateXTEUINT4TimeSeries(ts,N)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALCreateXTEUINT4TimeSeries() failed with error = %d.\n",fn,xlalErrno);
       XLAL_ERROR(fn,XLAL_ENOMEM);
    }
   
    /* fill in the timeseries array params */
    (*ts)->tstart = tstart;
    (*ts)->T = dt*N;
    (*ts)->deltat = dt; 
    (*ts)->energy[0] = event->energy[0];
    (*ts)->energy[1] = event->energy[1];
    LogPrintf(LOG_DEBUG,"%s : new binned timeseries parameters :\n",fn);
    LogPrintf(LOG_DEBUG,"%s : tstart = %6.12f\n",fn,(*ts)->tstart);
    LogPrintf(LOG_DEBUG,"%s : T = %f\n",fn,(*ts)->T);
    LogPrintf(LOG_DEBUG,"%s : dt = %e\n",fn,(*ts)->deltat);
    LogPrintf(LOG_DEBUG,"%s : length = %ld\n",fn,(*ts)->length);
    LogPrintf(LOG_DEBUG,"%s : energy = %d - %d\n",fn,(*ts)->energy[0],(*ts)->energy[1]);

  }

  /* initialise undefined as zero - for event data all times are good except for those defined by the GTI table */
  /* time marker events are not real events and are not counted anyway so we treat those times as good */
  for (i=0;i<(*ts)->length;i++) (*ts)->undefined[i] = 0;
  
  /* loop over each event and bin it if its not a time marker */
  for (i=0;i<event->nevents;i++) {

    /* set counter tracking 1st bit for each event */
    long int k = i*event->rowlength;

    /* define bin index corresponding to this event time with reference to the timestamp of the first event */
    long int idx = (long int)floor((stamps->dettime[i]-stamps->dettime[0])/(*ts)->deltat);
   
    /* only add to bin if event was NOT a time marker */
    (*ts)->data[idx] += (unsigned char)event->data[k];

  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;

}

/** Converts a FITS array data file to a binned timeseries.
 *
 * This function acts as a wrapper for XLALArrayDataToXTEUINT4TimeSeries which deals 
 * with single timeseries.  We enforce that each timeseries has the same start, end, 
 * and sampling time for consistency.
 * 
 */
int XLALArrayDataToXTEUINT4TimeSeriesArray(XTEUINT4TimeSeriesArray **ts,   /**< [out] a NULL timeseries */
					   FITSData *fits,                 /**< in] a FITSHeader structure */
					   REAL8 dt                        /**< [in] the desired sampling time */
					   )
{

  const CHAR *fn = __func__;         /* store function name for log output */
  INT8 i;                            /* counter */
  REAL8 newdt = dt;                  /* modified sampling time */

  /* check input and output pointers */
  if ((*ts) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries structure has non-null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (fits == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input FITSData structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (fits->header == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input FITSHeader structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  
  /* check timeseries consistency */
  for (i=0;i<fits->header->nXeCntcol;i++) {
    if ((fits->array[i]->deltat != fits->array[0]->deltat) || (fits->array[i]->length != fits->array[0]->length)) {
      LogPrintf(LOG_CRITICAL,"%s : inconsistent parameters between column timeseries.\n",fn,xlalErrno);
       XLAL_ERROR(fn,XLAL_EINVAL);
    }
  }

  /* check that requested sampling rate is at least as large as the original rate */
  {
    INT4 n = (INT4)ceil(fits->array[0]->deltat/dt);
    newdt = n*dt;
    LogPrintf(LOG_NORMAL,"%s : requested sampling time %6.12f -> %6.12f sec.\n",fn,dt,newdt);
  }
  
  /* allocate mem for output timeseries array */
  if (((*ts) = (XTEUINT4TimeSeriesArray *)LALCalloc(1,sizeof(XTEUINT4TimeSeriesArray))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for XTEUINT4TimeSeriesArray.\n",fn,xlalErrno);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*ts)->ts = (XTEUINT4TimeSeries **)LALCalloc(fits->header->nXeCntcol,sizeof(XTEUINT4TimeSeries *))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for XTEUINT4TimeSeries pointers.\n",fn,xlalErrno);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  
  /* loop over each relevant FITS data column */
  for (i=0;i<fits->header->nXeCntcol;i++) {
     
    if (XLALArrayDataToXTEUINT4TimeSeries(&((*ts)->ts[i]),fits->array[i],fits->stamps,newdt)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALEventDataToXTEUINT4TimeSeries() failed with error = %d\n",fn,xlalErrno);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
    LogPrintf(LOG_DEBUG,"%s : input FITSdata structure has null pointer.\n",fn);
  
    /* add column name info */
    snprintf((*ts)->ts[i]->colname,STRINGLENGTH,fits->header->colname[i]);

  }
  
  /* update header info */
  strncpy((*ts)->objectname,fits->header->objectname,STRINGLENGTH);
  strncpy((*ts)->obsid,fits->header->obsid,STRINGLENGTH);
  strncpy((*ts)->apid,fits->header->apid,APIDLENGTH);
  (*ts)->bary = 0;
  (*ts)->lld = fits->header->LLD;
  (*ts)->length = fits->header->nXeCntcol;
  
  /* allocate mem for the header dump */
  if (((*ts)->headerdump = (CHAR *)XLALCalloc(1+strlen(fits->header->headerdump),sizeof(CHAR))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for FITS header dump copy.\n",fn,xlalErrno);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  strncpy((*ts)->headerdump,fits->header->headerdump,strlen(fits->header->headerdump)*sizeof(CHAR));
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
  
}

/** Converts a FITS file data structure containing science array (binned) data into a timeseries.
 *
 * Using the detector timestamps data we convert the FITS data into a continuous binned timeseries. 
 *
 */
int XLALArrayDataToXTEUINT4TimeSeries(XTEUINT4TimeSeries **ts,      /**< [out] a null timeseries */
				      XTEUINT4Vector *array,        /**< [in] a FITSdata structure */
				      BarycentricData *stamps,      /**< [in] timestamps data */
				      REAL8 dt                      /**< [in] the sampling time */
				      )
{

  const CHAR *fn = __func__;        /* store function name for log output */
  INT4 i,j;                         /* counters */
  INT2 *temp_undefined = NULL;      /* temporary data quality vector */

  /* check input and output pointers */
  if ((*ts) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries structure has non-null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (array == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input FITSdata structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (stamps == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input FITSHeader structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* compute number of samples in the continuous timeseries */
  /* define the time span in the detector frame from start of first bin to end of last bin */
  {
    REAL8 tstart = stamps->dettime[0];
    REAL8 tempspan = stamps->dettime[stamps->length-1] - tstart + stamps->dtrow;
    INT8 N = (INT8)(tempspan/dt);
   
    /* allocate mem for output timeseries */
    if (XLALCreateXTEUINT4TimeSeries(ts,N)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALCreateXTEUINT4TimeSeries() failed with error = %d.\n",fn,xlalErrno);
       XLAL_ERROR(fn,XLAL_ENOMEM);
    }

    /* fill in params */
    (*ts)->tstart = tstart;
    (*ts)->T = dt*N;
    (*ts)->deltat = dt;
    (*ts)->energy[0] = array->energy[0];
    (*ts)->energy[1] = array->energy[1];
    LogPrintf(LOG_DEBUG,"%s : new binned timeseries parameters :\n",fn);
    LogPrintf(LOG_DEBUG,"%s : tstart = %6.12f\n",fn,(*ts)->tstart);
    LogPrintf(LOG_DEBUG,"%s : T = %f\n",fn,(*ts)->T);
    LogPrintf(LOG_DEBUG,"%s : dt = %e\n",fn,(*ts)->deltat);
    LogPrintf(LOG_DEBUG,"%s : length = %ld\n",fn,(*ts)->length);
    LogPrintf(LOG_DEBUG,"%s : energy = %d - %d\n",fn,(*ts)->energy[0],(*ts)->energy[1]);
  }

  /* allocate memory for a temporary data quality vector */
  if ((temp_undefined = (short int *)LALCalloc((*ts)->length,sizeof(short int))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for temporary data quality vector.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* initialise as zero and all undefined - therefore all gaps will be automatically undefined */
  for (i=0;i<(*ts)->length;i++) {
    (*ts)->data[i] = 0;
    (*ts)->undefined[i] = 1;
    temp_undefined[i] = 0;
  }

  /* loop over each timestamp and fill in timeseries */
  /* remember that times are in reference to bin edges */
  for (i=0;i<stamps->length;i++) {
   
    /* loop over data within timestamp span */
    for (j=0;j<array->rowlength;j++) {
      
      /* the index in the fits data in the j'th element in the i'th row */
      long int idxin = i*array->rowlength + j;
          
      /* the index in the output timeseries of the fits data in the j'th element in the i'th row */
      long int idxout = floor((stamps->dettime[i] - stamps->dettime[0] + j*array->deltat)/(*ts)->deltat);
  
      /* add corresponding data to correct time output bin */
      (*ts)->data[idxout] += array->data[idxin];
      temp_undefined[idxout] += array->undefined[idxin];
        
    }
    
  }

  /* make sure the undefined flags are either 0 or 1 - if greater than 1 then set to 1 */
  for (i=0;i<(*ts)->length;i++) {
    if (temp_undefined[i] > 1) (*ts)->undefined[i] = 1;
    else (*ts)->undefined[i] = 0;
  }

  /* free memory */
  XLALFree(temp_undefined);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;

}

/** Allocates memory for a XTEUINT4TimeSeries.
 *
 */
int XLALCreateXTEUINT4TimeSeries(XTEUINT4TimeSeries **x,   /**< [out] a null timeseries */
				 INT8 N                    /**< [in] the desired number of elements in the timeseries */
				 )
{
  
  const char *fn = __func__;        /* store function name for log output */
  
  /* check input */
  if (N<1) {
    LogPrintf(LOG_CRITICAL,"%s : tried to allocate an XTEUINT4TimeSeries with non-positive size.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if ((*x) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries structure does not have null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  if (((*x) = (XTEUINT4TimeSeries *)LALCalloc(1,sizeof(XTEUINT4TimeSeries))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for an XTEUINT4TimeSeries structure.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  (*x)->length = N;
  if (((*x)->data = (UINT4 *)LALCalloc(N,sizeof(UINT4))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for an XTEUINT4TimeSeries data vector.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*x)->undefined = (char *)LALCalloc(N,sizeof(char))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for an XTEUINT4TimeSeries data quality vector.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
  
}

/** Frees a FITS data structure.
 *
 */
int XLALFreeFITSData(FITSData *x    /**< [in/out] a FITS data structure */
		     )
{
  
  const char *fn = __func__;        /* store function name for log output */
  INT4 i;

  /* check input */
  if (x == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input FITSdata structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  if (x->event != NULL) {
    for (i=0;i<x->header->nXeCntcol;i++) {
      if (x->event[i] != NULL) {
	if (x->event[i]->data != NULL) XLALFree(x->event[i]->data);
	if (x->event[i]->undefined != NULL) XLALFree(x->event[i]->undefined);
	XLALFree(x->event[i]);
      }
    }
    XLALFree(x->event);
  }
  if (x->array != NULL) {
    for (i=0;i<x->header->nXeCntcol;i++) {
      if (x->array[i] != NULL) {
	if (x->array[i]->data != NULL) XLALFree(x->array[i]->data);
	if (x->array[i]->undefined != NULL) XLALFree(x->array[i]->undefined);
	XLALFree(x->array[i]);
      }
    }
    XLALFree(x->array);
  }
  if (x->header != NULL) {
    XLALFree(x->header->XeCntcolidx);
    for (i=0;i<x->header->nXeCntcol;i++) XLALFree(x->header->colname[i]);
    XLALFree(x->header->colname);
    XLALFree(x->header->headerdump);
    XLALFree(x->header);
  }
  if (x->gti != NULL) {
    if (XLALFreeGTIData(x->gti)) {
      LogPrintf(LOG_CRITICAL,"%s : unable to free GTI data with error = %d\n",fn,xlalErrno);
       XLAL_ERROR(fn,XLAL_EINVAL);
    }
  }
  if (x->stamps != NULL) {
    if (XLALFreeBarycentricData(x->stamps)) {
      LogPrintf(LOG_CRITICAL,"%s : unable to free barycentric data with error = %d\n",fn,xlalErrno);
       XLAL_ERROR(fn,XLAL_EINVAL);
    }
    
  }
  XLALFree(x);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
  
}

/** Frees an XTEUINT4TimeSeriesArray data structure.
 *
 */
int XLALFreeXTEUINT4TimeSeriesArray(XTEUINT4TimeSeriesArray *ts   /**< [in/out] a timeseries */
				    )
{

  const char *fn = __func__;        /* store function name for log output */
  INT4 i;                           /* counter */

  /* check input */
  if (ts == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries array has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  
  if (ts->ts != NULL) {
    for (i=0;i<ts->length;i++) {    
      if (ts->ts[i] != NULL) {

	if (XLALFreeXTEUINT4TimeSeries(ts->ts[i])) {
	  LogPrintf(LOG_CRITICAL,"%s : unable to free array data timeseries with error = %d\n",fn,xlalErrno);
	   XLAL_ERROR(fn,XLAL_EFAULT);
	}
      
      }
    }
    XLALFree(ts->ts);
  }
  XLALFree(ts->headerdump);
  XLALFree(ts->comment);
  XLALFree(ts);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;

}

/** Frees an XTEUINT4TimeSeries data structure.
 *
 */
int XLALFreeXTEUINT4TimeSeries(XTEUINT4TimeSeries *ts   /**< [in/out] a timeseries */
				  )
{

  const char *fn = __func__;        /* store function name for log output */
  
  /* check input */
  if (ts == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries vector has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  if (ts->data != NULL) XLALFree(ts->data);
  if (ts->undefined != NULL) XLALFree(ts->undefined);
  XLALFree(ts);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;

}

/** Sets timeseries array samples NOT in the GTI table as undefined and zeros the corresponding data.
 *
 * This is a wrapper for XLALApplyGTIToXTEUINT4TimeSeries.
 *
 */
int XLALApplyGTIToXTEUINT4TimeSeriesArray(XTEUINT4TimeSeriesArray **ts,    /**< [in/out] the timeseries */
					  GTIData *gti                     /**< [in] a GTIdata structure containing GTI information */
					  )
{
  
  const CHAR *fn = __func__;        /* store function name for log output */
  INT4 i;

   /* check input and output pointers */
  if ((*ts) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries structure has a null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (gti == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input GTIData structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  
  /* loop over each relavent column */
  for (i=0;i<(*ts)->length;i++) {
    
    if (XLALApplyGTIToXTEUINT4TimeSeries(&((*ts)->ts[i]),gti)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALApplyGTIToXTEUINT4TimeSeries() failed with error = %d\n",fn,xlalErrno);
       XLAL_ERROR(fn,XLAL_EFAULT);
    }
    LogPrintf(LOG_DEBUG,"%s : applied GTI table to data from column (%d/%d).\n",fn,i+1,(*ts)->length);
  }
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;

}

/** Sets timeseries samples NOT in the GTI table as undefined and zeros the corresponding data.
 *
 */
int XLALApplyGTIToXTEUINT4TimeSeries(XTEUINT4TimeSeries **ts,    /**< [in/out] the timeseries */
				     GTIData *gti                /**< [in] a GTIData structure containing GTI information */
				     )
{
  
  const char *fn = __func__;        /* store function name for log output */
  int i,k;
  long int newstartindex,newendindex,newN;
  GTIData *bti = NULL;

  /* check input */
  if ((*ts) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries vector has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (gti == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input GTIData structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* allocate memory for bad time intervals BTIs - more useful than good time intervals for vetoing data */
  if (XLALCreateGTIData(&bti,gti->length+1)) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for GTI start time data with error = %d.\n",fn,xlalErrno);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* compute first BTI from the timeseries start to the first GTI */
  bti->start[0] = (*ts)->tstart;
  bti->end[0] = gti->start[0];
  if (bti->end[0]<bti->start[0]) bti->start[0] = bti->end[0];  /* from rebinning the new start can be after the first gti start point */
  
  /* compute BTI (bad time interval) - inbetween gti times */
  for (i=1;i<gti->length;i++) {
    bti->start[i] = gti->end[i-1];
    bti->end[i] = gti->start[i];
  }

  /* compute last BTI from the last GTI to the timeseries end */
  bti->start[bti->length-1] = gti->end[gti->length-1];
  bti->end[gti->length] = (*ts)->tstart + (*ts)->T;
  if (bti->start[bti->length-1]>bti->end[bti->length-1]) bti->start[bti->length-1] = bti->end[bti->length-1];    /* from rebinning the new end can be before the last gti end point */
  for (i=0;i<bti->length;i++) LogPrintf(LOG_DEBUG,"%s : computed BTI : %f -> %f\n",fn,bti->start[i],bti->end[i]);
  
  /* loop over each BTI */
  for (i=0;i<bti->length;i++) {

    /* define start and end indices for this BTI - make sure that indices are within valid bounds */
    long int sindex = (bti->start[i] - (*ts)->tstart)/(*ts)->deltat;
    long int eindex = (bti->end[i] - (*ts)->tstart)/(*ts)->deltat;
    if (sindex<0) sindex = 0;
    if (eindex>(*ts)->length) eindex = (*ts)->length;
    
    /* loop over the data between GTIs and set as undefined and zero the data */
    for (k=sindex;k<eindex;k++) {
      (*ts)->undefined[k]=1;
      (*ts)->data[k] = 0;
    }

  }
  LogPrintf(LOG_DEBUG,"%s : zeroed and undefined non-GTI data in timeseries.\n",fn);
  
  /* cut out start and end BTI data */
  newstartindex = (bti->end[0]-(*ts)->tstart)/(*ts)->deltat;
  newendindex = (bti->start[gti->length]-(*ts)->tstart)/(*ts)->deltat;
  newN = newendindex - newstartindex;

  /* if the new start index is moved then slide all of the data */
  if (newstartindex>0) {
    long int count = 0;
    for (i=newstartindex;i<newendindex;i++) {
      (*ts)->data[count] = (*ts)->data[i];
      (*ts)->undefined[count] = (*ts)->undefined[i];
      count++;
    }
  }
 
  /* resize timeseries vector */
  if (XLALReallocXTEUINT4TimeSeries(ts,newN)) {
    LogPrintf(LOG_CRITICAL,"%s : failed to resize memory for XTEUINT4TimeSeries.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  
  /* update timeseries params */
  (*ts)->length = newN;
  (*ts)->T = newN*(*ts)->deltat;
  (*ts)->tstart += newstartindex*(*ts)->deltat;
  LogPrintf(LOG_DEBUG,"%s : changed start time and duration to %f and %f.\n",fn,(*ts)->tstart,(*ts)->T);

  /* free the bad time interval data */
  if (XLALFreeGTIData(bti)) {
    LogPrintf(LOG_CRITICAL,"%s : failed to free memory for GTIData with error = %d.\n",fn,xlalErrno);
     XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : freed the bad time interval data.\n",fn);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;

}

/** Resizes an XTEUINT4TimeSeries data structure.
 *
 */  
int XLALReallocXTEUINT4TimeSeries(XTEUINT4TimeSeries **ts,    /**< [in/out] the timeseries */
				  INT8 N                         /**< [in] the new length */
				  )
{
  const char *fn = __func__;        /* store function name for log output */
  
  /* check input */
  if (N<1) {
    LogPrintf(LOG_CRITICAL,"%s : tried to resize an XTEUINT4TimeSeries with non-positive size.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  
  if (((*ts)->data = (UINT4 *)LALRealloc((*ts)->data,N*sizeof(UINT4))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for an XTEUINT4TimeSeries data vector.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*ts)->undefined = (CHAR *)LALRealloc((*ts)->undefined,N*sizeof(CHAR))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for an XTEUINT4TimeSeries undefined vector.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  (*ts)->length = N;
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
  
}

/** Reduces the number of barycentric timestamps from a FITS data structure.
 *
 * We do this because we don't need fast sampling in time to perform accurate barycentering.  For event data 
 * there will be a timestamp for every event which means a lot of timestamps.  Here we reduce the number to
 * one timestamp per TIMESTAMPDELTAT seconds.
 * 
 */
int XLALReduceBarycentricData(BarycentricData **stamps   /**< [in/out] barycentered null timestamps vector */      
			      )
{
  
  const CHAR *fn = __func__;                 /* store function name for log output */
  INT8 i;                                     /* counter */
  INT4 m = 1;                                 /* another count */
  REAL8 thresh = 0.0;

  /* check input */
  if ((*stamps) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timestamps vector has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  
  /* initialise the first timestamp threshold */
  thresh = (*stamps)->dettime[0] + TIMESTAMPDELTAT;

  /* loop over each timestamp and record a timestamp every time the timestamp step threshold is passed */
  /* start at index 1 so that we keep the first timestamp */
  for (i=1;i<(*stamps)->length;i++) {
  
    /* if we've crossed the threshold then record timestamps */
    if ((*stamps)->dettime[i]>thresh) {

      INT4 n = ceil(((*stamps)->dettime[i] - (*stamps)->dettime[0])/TIMESTAMPDELTAT);
      thresh = (*stamps)->dettime[0] + n*TIMESTAMPDELTAT;
      (*stamps)->dettime[m] = (*stamps)->dettime[i];
      (*stamps)->barytime[m] = (*stamps)->barytime[i];
      m++;
    }
    
  }
  
  /* if last timestamp is different to original last timestamp then record last timestamp and new size of vector */
  if ((*stamps)->dettime[m-1]<(*stamps)->dettime[(*stamps)->length-1]) {
    (*stamps)->dettime[m] = (*stamps)->dettime[(*stamps)->length-1];
    (*stamps)->barytime[m] = (*stamps)->barytime[(*stamps)->length-1];
    m++;
  }
  (*stamps)->length = m;
  LogPrintf(LOG_DEBUG,"%s : reduced timestamps vector size to %d.\n",fn,(*stamps)->length);
  
  /* resize stamps vectors */
  if (((*stamps)->dettime = (double *)LALRealloc((*stamps)->dettime,(*stamps)->length*sizeof(double))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to resize memory for detector frame timestamps.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*stamps)->barytime = (double *)LALRealloc((*stamps)->barytime,(*stamps)->length*sizeof(double))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to resize memory for barycentric frame timestamps.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* output debugging information */
  {
    int M = (*stamps)->length > 5 ? 5 : (*stamps)->length;
    LogPrintf(LOG_DEBUG,"%s : extracted timestamps as :\n",fn);
    for (i=0;i<M;i++) LogPrintf(LOG_DEBUG,"%s : %6.12f (%6.12f)\n",
				fn,(*stamps)->dettime[i],(*stamps)->barytime[i]);
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
  
}

/** Allocates memory for a GTI structure.
 *
 */
int XLALCreateGTIData(GTIData **gti,      /**< [out] a null timeseries */
		      INT8 N              /**< [in] the desired number of elements in the timeseries */
		      )
{
  
  const char *fn = __func__;        /* store function name for log output */
  
  /* check input */
  if (N<1) {
    LogPrintf(LOG_CRITICAL,"%s : tried to allocate an barycenteredtimestamps with non-positive size.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if ((*gti) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timestamps structure does not have null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  if (((*gti) = (GTIData *)LALCalloc(1,sizeof(GTIData))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for GTI structure.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*gti)->start = (double *)LALCalloc(N,sizeof(double))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for GTI start vector.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*gti)->end = (double *)LALCalloc(N,sizeof(double))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for GTI end vector.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  (*gti)->length = N;

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
  
}
/** Frees memory for a GTI structure.
 *
 */
int XLALFreeGTIData(GTIData *gti   /**< [out] a null timeseries */
		    )
{
  
  const char *fn = __func__;        /* store function name for log output */
  
  /* check input */
  if (gti == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input GTI structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  if (gti->start != NULL) XLALFree(gti->start);
  if (gti->end != NULL) XLALFree(gti->end);
  XLALFree(gti);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
  
}

/** Allocates memory for a barycentric timestamps structure.
 *
 */
int XLALCreateBarycentricData(BarycentricData **stamps,   /**< [out] a null timeseries */
			      INT8 N                      /**< [in] the desired number of elements in the timeseries */
			      )
{
  
  const char *fn = __func__;        /* store function name for log output */
  
  /* check input */
  if (N<1) {
    LogPrintf(LOG_CRITICAL,"%s : tried to allocate an barycenteredtimestamps with non-positive size.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if ((*stamps) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timestamps structure does not have null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  if (((*stamps) = (BarycentricData *)LALCalloc(1,sizeof(BarycentricData))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for timestamps structure.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*stamps)->dettime = (double *)LALCalloc(N,sizeof(double))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for detector frame timestamps.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if (((*stamps)->barytime = (double *)LALCalloc(N,sizeof(double))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for barycentric frame timestamps.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  (*stamps)->length = N;
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
  
}

/** Allocates memory for a barycentric timestamps structure.
 *
 */
int XLALFreeBarycentricData(BarycentricData *stamps   /**< [out] a null timeseries */
			    )
{
  
  const char *fn = __func__;        /* store function name for log output */
  
  /* check input */
  if (stamps == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timestamps structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  if (stamps->dettime != NULL) XLALFree(stamps->dettime);
  if (stamps->dettime != NULL) XLALFree(stamps->barytime);
  XLALFree(stamps);
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
  
}
/** This function outputs an XTEUINT4TimeSeriesArray to a frame file or files.
 * 
 */
int XLALXTEUINT4TimeSeriesArrayToGTI(GTIData **gti,XTEUINT4TimeSeriesArray *ts) 
{

  static const char *fn = __func__;        /* store function name for log output */
  INT4 ngti = 0;                           /* the number of gti entries */
  CHAR *temp_undefined = NULL;             /* temporary pointer to master undefined vector */
  XTEUINT4TimeSeries *tempts = NULL;       /* temporary pointer to first timeseries */
  INT4 i;                                  /* counter */
  INT8 j;                                  /* counter */

  /* check input */
  if (ts == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if ((*gti) != NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input GTIData structure has non-null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* point temporary timeseries pointer to first timeseries */
  tempts = ts->ts[0];

  /* allocate memory for a master list of undefined samples */
  if ((temp_undefined = (CHAR *)LALCalloc(tempts->length,sizeof(CHAR))) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for master undefined vector.\n",fn);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* initialise master list and then loop over each timeseries and make a master list of undefined samples */
  for (j=0;j<tempts->length;j++) temp_undefined[j] = 0;
  for (i=0;i<ts->length;i++) {
    for (j=0;j<ts->ts[i]->length;j++) temp_undefined[j] += ts->ts[i]->undefined[j];
  }

  /* compute the number of GTIs */  
  {
    BOOLEAN flag = FALSE;
    for (j=0;j<tempts->length;j++) {
      if ((!flag) && (tempts->undefined[j] == 0)) flag = TRUE;
      if (flag && tempts->undefined[j]) {
	flag = FALSE;
	ngti++;
      }
    }
    if (flag) ngti++;
  }
  LogPrintf(LOG_DEBUG,"%s : final GTI table has %d entries.\n",fn,ngti);

  /* allocate memory for a GTI table */
  if (XLALCreateGTIData(gti,ngti)) {
    LogPrintf(LOG_CRITICAL,"%s : failed to allocate memory for GTI with error = %d.\n",fn,xlalErrno);
     XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* go back through the data quality and fill in the GTI table */
  {
    BOOLEAN flag = FALSE;
    i = 0;
    for (j=0;j<tempts->length;j++) {
      if ((!flag) && (tempts->undefined[j] == 0)) {
	flag = TRUE;
	(*gti)->start[i] = tempts->tstart + (REAL8)j*tempts->deltat; 
      }
      if (flag && tempts->undefined[j]) {
	flag = FALSE;
	(*gti)->end[i] = tempts->tstart + (REAL8)j*tempts->deltat; 
	i++;
      }
    }
    if (flag) (*gti)->end[i] = tempts->tstart + tempts->T; 
  }
  
  /* debugging */
  LogPrintf(LOG_DEBUG,"%s : master GTI times :\n",fn);
  for (i=0;i<(*gti)->length;i++)  LogPrintf(LOG_DEBUG,"%s : %6.12f -> %6.12f\n",fn,(*gti)->start[i],(*gti)->end[i]);
  
  /* free memory */
  XLALFree(temp_undefined);

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;

}

/** This function outputs an XTEUINT4TimeSeriesArray to a frame file or files.
 * 
 */
int XLALXTEUINT4TimeSeriesArrayToFrames(XTEUINT4TimeSeriesArray *ts,      /**< [in] the input xte timeseries array */
					CHAR *outputdir                  /**< [in] the output directory name */
					)
{

  /* FIXME - need to put comments into the frame */
 
  static const char *fn = __func__;        /* store function name for log output */  
  INT4 i;                                  /* counter */
  INT4 k;                                  /* counter */
  XTEUINT4TimeSeries *tempts = NULL;       /* pointer to first timeseries */
  GTIData *gti = NULL;                          /* final GTI information */

  /* check input */
  if (ts == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries structure has null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* point temporary pointer to first timeseries */
  tempts = ts->ts[0];

  /* make sure all timeseries have identical start, duration and sampling times */
  for (i=0;i<ts->length;i++) {
    if ((ts->ts[i]->length != tempts->length) || 
	(ts->ts[i]->tstart != tempts->tstart) || 
	(ts->ts[i]->T != tempts->T)) {
      LogPrintf(LOG_CRITICAL,"%s : input timeseries have different parameters.\n",fn);
       XLAL_ERROR(fn,XLAL_EINVAL);
    }
  }

  /* make a new GTI table based on data quality from all timeseries */
  if (XLALXTEUINT4TimeSeriesArrayToGTI(&gti,ts)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALXTEUINT4TimeSeriesToTGI() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }
  LogPrintf(LOG_DEBUG,"%s : generated a master GTI table\n",fn);

  /* loop over each GTI entry such that we output a single frame for each one */
  for (k=0;k<gti->length;k++) {

    char outputfile[STRINGLENGTH];           /* stores name of output frame file */
    struct FrameH *outFrame   = NULL;        /* frame data structure */
    LIGOTimeGPS epoch;                       /* stores the timeseries epoch */
    INT8 sidx;                               /* the integer GPS second starting index of the data */
    INT8 N;                                  /* the new number of data samples */

    /* enforce integer seconds duration */
    INT4 start = (INT4)ceil(gti->start[k]);  /* the new end time of the data */
    INT4 end = (INT4)floor(gti->end[k]);     /* the new end time of the data */
    INT4 T = end - start;                    /* the new observation time */
    LogPrintf(LOG_DEBUG,"%s : segment %d -> %d (%d sec)\n",fn,start,end,T);
    
    /* if the segment is long enough to warrant making a frame */
    if (T>=MINFRAMELENGTH) {

      /* convert start epoch to a GPS structure - we enforce integer GPS start times for simplicity */
      epoch.gpsSeconds = start;
      epoch.gpsNanoSeconds = 0;
      LogPrintf(LOG_DEBUG,"%s : modifying start time to integer epoch %d\n",fn,epoch.gpsSeconds);
      
      /* compute new number of samples and starting index */
      sidx = (INT8)floor(0.5 + (epoch.gpsSeconds - tempts->tstart)/tempts->deltat);
      N = (INT8)floor(T/tempts->deltat);
      LogPrintf(LOG_DEBUG,"%s : outputting %ld samples to frame\n",fn,N);
      
      /* construct file name - we use the LIGO format <DETECTOR>-<COMMENT>-<GPSSTART>-<DURATION>.gwf */
      /* the comment field we sub-format into <INSTRUMENT>_<FRAME>_<SOURCE>_<OBSID_APID>_<LLD>_<MINENERGY>_<MAXENERGY> */
      if (ts->bary) snprintf(outputfile,STRINGLENGTH,"%s/X-PCA_SSB_%s_%s_%s_%d_%d_%d-%d-%d.gwf",
			     outputdir,ts->objectname,ts->obsid,ts->apid,ts->lld,
			     ts->ts[0]->energy[0],ts->ts[0]->energy[1],epoch.gpsSeconds,T);
      else snprintf(outputfile,STRINGLENGTH,"%s/X-PCA_DET_%s_%s_%s_%d_%d_%d-%d-%d.gwf",
			     outputdir,ts->objectname,ts->obsid,ts->apid,ts->lld,
			     ts->ts[0]->energy[0],ts->ts[0]->energy[1],epoch.gpsSeconds,T);
      LogPrintf(LOG_DEBUG,"%s : output file = %s\n",fn,outputfile);
      
      /* generate a frame data structure - last threee inputs are [project, run, frnum, detectorFlags] */
      if ((outFrame = XLALFrameNew(&epoch,(REAL8)T,"XTE_PCA",1,0,0)) == NULL) {
	LogPrintf(LOG_CRITICAL, "%s : XLALFrameNew() failed with error = %d.\n",fn,xlalErrno);
	XLAL_ERROR(fn,XLAL_EFAILED);
      }
      LogPrintf(LOG_DEBUG,"%s : set-up frame structure\n",fn);

      /* loop over original FITS file columns */
      for (i=0;i<ts->length;i++) {

	INT4TimeSeries *output = NULL;           /* temporary output timeseries */
	UINT4 j;                                 /* counter */
	CHAR channelname[STRINGLENGTH];          /* string used to name each channel in the frame */

	/* define current channel name */
	snprintf(channelname,STRINGLENGTH,"%s%d",xtechannelname,i);

	/* create empty timeseries - this is INT4 not UINT4 because there is no frame writing function for UINT4 */
	if ((output = XLALCreateINT4TimeSeries(channelname,&epoch,0,tempts->deltat,&lalDimensionlessUnit,N)) == NULL) {
	  LogPrintf(LOG_CRITICAL, "%s : XLALCreateINT4TimeSeries() failed to allocate an %d length timeseries with error = %d.\n",fn,N,xlalErrno);
	  XLAL_ERROR(fn,XLAL_ENOMEM);  
	}
	LogPrintf(LOG_DEBUG,"%s : allocated memory for temporary timeseries\n",fn);
	
	/* fill in timeseries starting from first integer second */
	for (j=0;j<output->data->length;j++) output->data->data[j] = (INT4)ts->ts[i]->data[j+sidx];
	LogPrintf(LOG_DEBUG,"%s : populated temporary timeseries\n",fn);
	  
	/* add timeseries to frame structure */
	if (XLALFrameAddINT4TimeSeriesProcData(outFrame,output)) {
	  LogPrintf(LOG_CRITICAL, "%s : XLALFrameAddINT4TimeSeries() failed with error = %d.\n",fn,xlalErrno);
	  XLAL_ERROR(fn,XLAL_EFAILED);
	}
	LogPrintf(LOG_DEBUG,"%s : added timeseries from col %s for epoch %d to frame structure\n",fn,ts->ts[i]->colname,epoch.gpsSeconds);

	/* Here's where we add extra information into the frame */ 
	{
	  /* CHAR energycomment[STRINGLENGTH];   */      /* string used to store the energy channel range for output to the frame */
	  CHAR channelcomment[STRINGLENGTH];       /* string used to store the channel/column name for output to the frame */
	 
	  FrHistoryAdd( outFrame,ts->headerdump);
	  FrHistoryAdd( outFrame,ts->comment);
	 /*  snprintf(energycomment,STRINGLENGTH,"\nMINENERGY = %d\nMAXENERGY = %d",ts->ts[i]->energy[0],ts->ts[i]->energy[1]);  */
/* 	  FrHistoryAdd( outFrame,energycomment); */
	  snprintf(channelcomment,STRINGLENGTH,"\n%s = %s",channelname,ts->ts[i]->colname); 
	  FrHistoryAdd( outFrame,channelcomment);
	 	 
	}

	/* free timeseries */
	XLALDestroyINT4TimeSeries(output);
	
      }

      /* write frame structure to file (opens, writes, and closes file) - last argument is compression level */
      if (XLALFrameWrite(outFrame,outputfile,1)) {
	LogPrintf(LOG_CRITICAL, "%s : XLALFrameWrite() failed with error = %d.\n",fn,xlalErrno);
	XLAL_ERROR(fn,XLAL_EFAILED);
      }

    }
    else {
      LogPrintf(LOG_DEBUG, "%s : segment %f -> %f not long enough so no frame generated.\n",fn,gti->start[k],gti->end[k]);
    }
  
  }

  if (XLALFreeGTIData(gti)) {
    LogPrintf(LOG_CRITICAL,"%s : XLALFreeGTIData() failed with error = %d\n",fn,xlalErrno);
    return 1;
  }

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
   return XLAL_SUCCESS;
 
}

/** This function barycenters an XTEUINT4TimeSeriesArray using a set of barycentric timestamps.
 * 
 * This is a wrapper for XLALBarycenterXTEUINT4TimeSeries.
 * 
 */
int XLALBarycenterXTEUINT4TimeSeriesArray(XTEUINT4TimeSeriesArray **ts,       /**< [in/out] timeseries */
					  BarycentricData *stamps             /**< [in] vector of pairs of detector/barycentered timestamps */
					  )
{
  
  static const char *fn = __func__;           /* store function name for log output */
  INT8 i;                            /* counter */
  
  /* check input and output pointers */
  if ((*ts) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries structure has null pointer.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (stamps == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input barycentric timestamps structure has non-null pointer.\n",fn);
     XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* loop over each relevant column */
  for (i=0;i<(*ts)->length;i++) {
    
    if (XLALBarycenterXTEUINT4TimeSeries(&((*ts)->ts[i]),stamps)) {
      LogPrintf(LOG_CRITICAL,"%s : XLALBarycenterXTEUINT4TimeSeries() failed with error = %d\n",fn,xlalErrno);
      XLAL_ERROR(fn,XLAL_EFAULT);
    }
    LogPrintf(LOG_DEBUG,"%s : barycentered data from column (%d/%d).\n",fn,i+1,(*ts)->length);
  
  }

  (*ts)->bary = 1;

  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;

}

/** This function barycenters an XTEUINT4TimeSeries using a set of barycentric timestamps.
 * 
 * The input detector frame timeseries is replaced with the output barycentered timeseries.
 *
 */
int XLALBarycenterXTEUINT4TimeSeries(XTEUINT4TimeSeries **ts,       /**< [in/out] timeseries */
				     BarycentricData *stamps       /**< [in] vector of pairs of detector/barycentered timestamps */
				     )
{
  
  static const char *fn = __func__;           /* store function name for log output */
  long int i;                                 /* counter */
  XTEUINT4TimeSeries *tempts = NULL;       /* temporary storage for barycentered timeseries */
  double start_det,end_det;                   /* start and end times of detector frame timeseries */
  double start_bary,end_bary;                 /* start and end times of barycentric frame timeseries */
  gsl_interp_accel *acc = NULL;               /* structure for accelerating gsl interpolation */
  gsl_spline *spline = NULL;                  /* structure for gsl spline interpolation */
  
  /* check input and output pointers */
  if ((*ts) == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input timeseries structure has null pointer.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
  if (stamps == NULL) {
    LogPrintf(LOG_CRITICAL,"%s : input barycentric timestamps structure has non-null pointer.\n",fn);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }

  /* check timestamp gaps - they can't be too far apart */
  for (i=0;i<stamps->length-1;i++) {
    if (stamps->dettime[i+1]-stamps->dettime[i]>MAXTIMESTAMPDELTAT) {
      LogPrintf(LOG_CRITICAL, "%s : timestamp spacing is too large for accurate barycentering.\n",fn);
      XLAL_ERROR(fn,XLAL_EINVAL);
    }
  }

  /* gsl memory allocation for interpolation of the timestamps */
  if ((acc = gsl_interp_accel_alloc()) == NULL) {
    LogPrintf(LOG_CRITICAL, "%s : gsl_interp_accel_alloc() failed to allocate memory for acceleration.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  if ((spline = gsl_spline_alloc(gsl_interp_cspline,stamps->length)) == NULL) {
    LogPrintf(LOG_CRITICAL, "%s : gsl_spline_alloc() failed to allocate memory for acceleration.\n",fn);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }
  
  /* select latest start time and earliest end time that are within the timestamps AND data limits */
  start_det = (*ts)->tstart > stamps->dettime[0] ? (*ts)->tstart : stamps->dettime[0];
  end_det = ((*ts)->tstart + (*ts)->T) < stamps->dettime[stamps->length-1] ? ((*ts)->tstart + (*ts)->T) : stamps->dettime[stamps->length-1];
  
  /* check that times are still valid */
  if (start_det>=end_det) {
    LogPrintf(LOG_CRITICAL, "%s : GPS start (%f) > GPS end time (%f), unable to generate a barycentric time series.\n",fn,start_det,end_det);
    XLAL_ERROR(fn,XLAL_EINVAL);
  }
   
  /* set up gsl interpolation for det -> bary to get start and end times of barycentric timeseries */
  if (gsl_spline_init(spline,stamps->dettime,stamps->barytime,stamps->length)) {
    LogPrintf(LOG_CRITICAL, "%s : gsl_spline_init() failed for det -> bary interpolation.\n",fn);
    XLAL_ERROR(fn,XLAL_EFAULT);
  }

  /* compute barycentric time of first and last sample */
  start_bary = gsl_spline_eval(spline,start_det,acc);
  end_bary = gsl_spline_eval(spline,end_det,acc);
  LogPrintf(LOG_DEBUG,"%s : barycentric GPS start %f\n",fn,start_bary);
  LogPrintf(LOG_DEBUG,"%s : barycentric GPS end = %f\n",fn,end_bary);
  LogPrintf(LOG_DEBUG,"%s : barycentric duration = %f\n",fn,end_bary-start_bary);

  {
    /* define number of barycentric timesamples */
    long int N = (long int)((end_bary-start_bary)/(*ts)->deltat);
    
    /* allocate temp memory first */
    XLALCreateXTEUINT4TimeSeries(&tempts,N);
    LogPrintf(LOG_DEBUG,"%s : created a temporary short int timeseries with %ld points\n",fn,N);
  }

  /* set up gsl interpolation for bary -> det for actual resampling */
  if (gsl_spline_init(spline,stamps->barytime,stamps->dettime,stamps->length)) {
    LogPrintf(LOG_CRITICAL, "%s : gsl_spline_init() failed for bary -> det interpolation.\n",fn);
    XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : initialised the spline interpolation\n",fn);

  /* loop over each evenly spaced time sample in barycentric frame and interpolate to generate barycentric timeseries */
  for (i=0;i<tempts->length;i++) {
    
    REAL8 bt = start_bary + (REAL8)i*(*ts)->deltat;
    REAL8 delta_dettime = gsl_spline_eval(spline,bt,acc) - (*ts)->tstart;
    INT8 idx = floor(0.5 + delta_dettime/(*ts)->deltat);
    tempts->data[i] = (*ts)->data[idx];
    tempts->undefined[i] = (*ts)->undefined[idx];

  }
  LogPrintf(LOG_DEBUG,"%s : performed barycentering using nearest bin interpolation\n",fn);

  /* reallocate input time series and fill it with barycentric timeseries and params */
  if (XLALReallocXTEUINT4TimeSeries(ts,tempts->length)) {
    LogPrintf(LOG_CRITICAL, "%s : XLALReallocXTEUINT4TimeSeries() failed with error = %d.\n",fn,xlalErrno);
    XLAL_ERROR(fn,XLAL_ENOMEM);
  }

  /* fill it with barycentric timeseries and params */
  for (i=0;i<tempts->length;i++) {
    (*ts)->data[i] = tempts->data[i];
    (*ts)->undefined[i] = tempts->undefined[i];
  }
  (*ts)->length = tempts->length;
  (*ts)->T = (*ts)->length*(*ts)->deltat;
  (*ts)->tstart = start_bary;
  LogPrintf(LOG_DEBUG,"%s : resized original timeseries and replaced it with barycentered timeseries\n",fn);

  /* free temp timeseries */
  if (XLALFreeXTEUINT4TimeSeries(tempts)) {
    LogPrintf(LOG_CRITICAL, "%s : XLALFreeXTEUINT4TimeSeries() failed with error = %d.\n",fn,xlalErrno);
    XLAL_ERROR(fn,XLAL_EFAULT);
  }
  LogPrintf(LOG_DEBUG,"%s : freed temporary memory\n",fn);

  /* free interpolation structures */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  LogPrintf(LOG_DEBUG,"%s : freed interpolation structures\n",fn);
  
  /* debugging */
  {
    int M = (*ts)->length > 5 ? 5 : (*ts)->length;
    LogPrintf(LOG_DEBUG,"%s : barycentered timeseries :\n",fn);
    for (i=0;i<M;i++) LogPrintf(LOG_DEBUG,"%s : %6.12f %d (%d)\n",fn,start_bary+(double)i*(*ts)->deltat,(*ts)->data[i],(*ts)->undefined[i]);
  }
  
  LogPrintf(LOG_DEBUG,"%s : leaving.\n",fn);
  return XLAL_SUCCESS;
  
}
 
