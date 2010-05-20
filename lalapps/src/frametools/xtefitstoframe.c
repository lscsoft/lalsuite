/*
*  Copyright (C) 2007 Badri Krishnan, Chris Messenger, Jolien Creighton, Reinhard Prix
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

/*********************************************************************************/
/*               XTE FITS file to frame file format conversion code              */
/*                                                                               */
/*		                   C. Messenger                                  */
/*                                                                               */
/*                 Albert Einstein Institute - started May 2010                  */
/*********************************************************************************/
/*                                                                               */
/* This code is designed to read in an XTE FITS data file specifically from the  */
/* pca (proportional counting array) onboard the (R)XTE sattelite.  It reads in  */
/* either individually time-tagged photons or pre-binned photon counts.  If the  */
/* input file has been pre-processed such that it also contains barycentric time */
/* information then the user can specify that the output timeseries be           */
/* barycentered.  The final timeseries is output in frame format.                */
/*                                                                               */
/*********************************************************************************/

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
#include <lalappsfrutils.h>

/***********************************************************************************************/
/* some global constants */
#define MAX_DT 0.000976562            /* the maximum sample time we will use (equiv max freq 512 Hz) */
#define MIN_DT 0.000244140625         /* the minimum sample time we will use (equiv max freq 2048 Hz) */
#define GPSMINUSTAI 441417609         /* IMPORTANT, DO NOT TOUCH = difference between GPS and TAI in seconds on January 1, 1994, at 00:00:28 */
#define TIMESTAMP_DELTAT 1            /* the spacing used when reducing the number of event time stamps */

/***********************************************************************************************/
/* error codes */
#define XTEFITSTOFRAMEC_ENULL 		1
#define XTEFITSTOFRAMEC_ESYS     	2
#define XTEFITSTOFRAMEC_EINPUT   	3
#define XTEFITSTOFRAMEC_MSGENULL 	"Arguments contained an unexpected null pointer"
#define XTEFITSTOFRAMEC_MSGESYS		"System call failed (probably file IO)"
#define XTEFITSTOFRAMEC_MSGEINPUT   	"Invalid input"

/***********************************************************************************************/
/* structure for storing the content read from a FITS file */
typedef struct {
  char file[256];                  /* the full path of the fits file */
  char filename[256];              /* the name of the fits file */
  int event;                       /* is the data events (1) or array (0) */
  int bary;                        /* does the data contain barycentering information */
  short int *data;                 /* binned data */
  char *events;                    /* event data */
  char *undefined;                 /* a data quality flag */
  double *gtistart;                /* the start times of good data */
  double *gtiend;                  /* the end times of good data */
  double *dettime;                 /* a vector of times */
  double *barytime;                /* a vector of barycentric times */
  long int ngtirows;               /* the number of gti segments */
  long int rowlength;              /* the number of data per timestamp */
  long int nbytesperrow;           /* the number of bytes per row for event data */
  long int N;                      /* the total number of data elements */
  long int Nevents;                /* the total number of events */
  long int ncols;                  /* the number of columns in the data table */
  long int nrows;                  /* the number of rows = number of timestamps */
  double deltat;                   /* the deltat for this time series */
  double eventdeltat;              /* the event sampling time */
  double dttimestamp;              /* the time steps for the timestamps */
  double T;                        /* the total length of data */
  double tstart;                   /* the start time of the data */
  double gpsoffset;                /* the number to be added to TAI times to get to GPS */
  double ra;                       /* the source right ascension */
  double dec;                      /* the source declination */
  char objectname[256];            /* the name of the source object */
  char obsid[256];                 /* the observation ID */
  char mode[256];                  /* stores the data mode */
  int LLD;                         /* flag for whether it's LLD2 data */
  int nchannels;                   /* the number of channels */
  long int channelsize;            /* the length of a channel */
  int minenergy;                   /* the min energy in the RXTE 255 scale */
  int maxenergy;                   /* the max energy in the RXTE 255 scale */
} FITSdata;

/* structure to store TDDES2 DDL data */
typedef struct {
  int minenergy;                   /* minimum energy channel */
  int maxenergy;                   /* maximum energy channel */
  double deltat;                   /* the sampling time */  
  double offset;                   /* the tim offset */
  int nsamples;                    /* the number of samples */
} tddesparams;

/* structure for storing a timeseries of short integers */
typedef struct {
  short int *data;                /* the data */
  char *undefined;                /* a quality flag */
  double deltat;                  /* the time step size */
  long int N;                     /* the number of data points */ 
  double tstart;                  /* the start time */
  double T;                       /* the time span */
} shortinttimeseries;

/* structure for storing vectors of barycentered timestamps */
typedef struct {
  double *dettime;                /* the detector timestamps */
  double *barytime;               /* the barycentered timestamps */
  long int N;                     /* the number of timestamps */ 
} barycenteredtimestamps;

/* user input variables */
typedef struct { 
  char inputfile[256];             /* directory containing input FITS files */
  char outputdir[256];             /* name of output directory */
  double deltat;                   /* the desired sampling time */
  int bary;                        /* flag for performing barycentering */
} UserVariables_t;

/***********************************************************************************************/
/* global variables */
int debug = 1;                     /* initialise the debug flag */
char clargs[512];                  /* the command line args stored as a string */

/* keywords in FITS file header */
char string_OBJECT[] = "OBJECT";
char string_RA_NOM[] = "RA_NOM";
char string_DEC_NOM[] = "DEC_NOM";
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

/***********************************************************************************************/
/* define functions */
void InitUserVars(int argc, char *argv[],UserVariables_t *uvar);
void converttddes(char *tddes, tddesparams *energytimeparams);
void readfitsfile(char *filename,FITSdata *data,int col);
void printFITSerror(int FITSstatus);
void eventdata_to_shortinttimeseries(FITSdata *fits,shortinttimeseries *ts,double dt);
void arraydata_to_shortinttimeseries(FITSdata *fits,shortinttimeseries *ts,double dt);
void readFITSheader(fitsfile *fptr, FITSdata *fitsfiledata, int status);
void readFITSgti(fitsfile *fptr, FITSdata *fitsfiledata, int status);
void readFITSbinneddata(fitsfile *fptr,FITSdata *fitsfiledata,int col,int status);
void readFITSeventdata(fitsfile *fptr,FITSdata *fitsfiledata,int col,int status);
void readFITStimestamps(fitsfile *fptr,FITSdata *fitsfiledata,int status);
void createshortinttimeseries(shortinttimeseries *x,long int N);
void freeFITSdata(FITSdata *x);
void freeshortinttimeseries(shortinttimeseries *ts);
void applygtitable(shortinttimeseries *ts,FITSdata *fits);
void reallocshortinttimeseries(shortinttimeseries *ts,int N);
void extractbarytimestamps(FITSdata *fits,barycenteredtimestamps *stamps);
void barycentershortinttimeseries(shortinttimeseries *ts,barycenteredtimestamps *stamps);
void shortinttimeseriestoframes(shortinttimeseries *ts,char *outputdir);

/*************************************************************************************************/
/*                                                                                               */
/* start the main code */
/*                                                                                               */
/*************************************************************************************************/
int main( int argc, char *argv[] )  {

  UserVariables_t uvar;                 /* user input variables */
  FITSdata fitsdata;                    /* a FITS data structure */
  shortinttimeseries ts;                /* a timeseries structure */ 
  barycenteredtimestamps barystamps;    /* vector of barycentered timestamps */

  /* initiate the user variables */
  InitUserVars(argc,argv,&uvar);
  if (debug) fprintf(stdout,"STATUS : read in uservars\n");
      
  /* read in data from the input FITS file */
  if (debug) fprintf(stdout,"STATUS : about to read in fits data from file %s\n",uvar.inputfile);
  readfitsfile(uvar.inputfile,&fitsdata,2);
  if (debug) fprintf(stdout,"STATUS : read in fits data from file %s\n",uvar.inputfile);  
  
  /* convert fits data to timeseries data */
  if (debug) fprintf(stdout,"STATUS : about to convert FITS data structure to a timeseries\n");
  if (fitsdata.event) eventdata_to_shortinttimeseries(&fitsdata,&ts,uvar.deltat);
  else arraydata_to_shortinttimeseries(&fitsdata,&ts,uvar.deltat);
  if (debug) fprintf(stdout,"STATUS : converted FITS data structure to a timeseries\n");  
    
  /* apply GTI and undefined data flag */
  if (debug) fprintf(stdout,"STATUS : about to apply GTI information to data from file %s\n",uvar.inputfile);
  applygtitable(&ts,&fitsdata);
  if (debug) fprintf(stdout,"STATUS : applied GTI information to data from file %s\n",uvar.inputfile);  

  if (uvar.bary) {
    
    /* if barycentering then extract timestamps */
    if (debug) fprintf(stdout,"STATUS : about to extract barycentered timestamps from file %s\n",uvar.inputfile);
    extractbarytimestamps(&fitsdata,&barystamps);
    if (debug) fprintf(stdout,"STATUS : extracted barycentered timestamps from file %s\n",uvar.inputfile);  
 
    /* perform barycentering */
    if (debug) fprintf(stdout,"STATUS : about to extract barycentered timestamps from file %s\n",uvar.inputfile);
    barycentershortinttimeseries(&ts,&barystamps);
    if (debug) fprintf(stdout,"STATUS : extracted barycentered timestamps from file %s\n",uvar.inputfile); 

  }

  /* output data */
  if (debug) fprintf(stdout,"STATUS : about to output data to frame file(s)\n");
  shortinttimeseriestoframes(&ts,uvar.outputdir);
  if (debug) fprintf(stdout,"STATUS : output data to frame file(s)\n"); 

  /* free memory */
  freeFITSdata(&fitsdata);
  freeshortinttimeseries(&ts);
  
  return 0;
  
}

/*************************************************************************************************/
/*                                                                                               */
/* read input user arguments */
/*                                                                                               */
/*************************************************************************************************/
void InitUserVars(int argc, char *argv[],UserVariables_t *uvar)
{

  int i;

  /* initialise user variables */
  uvar->inputfile[0] = '\0'; 
  uvar->outputdir[0] = '\0'; 
  uvar->bary = 0;
  uvar->deltat = MIN_DT;

  /* parse command line arguments */
  i = 1;
  while (i<argc) {
    
    if ((strcmp(argv[i],"-i")==0)||(strcmp(argv[i],"--inputfile")==0)) {
      sprintf(uvar->inputfile,"%s",argv[i+1]);
    }
    if ((strcmp(argv[i],"-o")==0)||(strcmp(argv[i],"--outputdir")==0)) {
      sprintf(uvar->outputdir,"%s",argv[i+1]);
    }
    else if ((strcmp(argv[i],"-b")==0)||(strcmp(argv[i],"--bary")==0)) {
       uvar->bary = 1;
       i--;
    }
    else if ((strcmp(argv[i],"-t")==0)||(strcmp(argv[i],"--deltat")==0)) {
      uvar->deltat = atof(argv[i+1]);
    }
    else if ((strcmp(argv[i],"-d")==0)||(strcmp(argv[i],"--debug")==0)) {
      debug = atol(argv[i+1]);
    }
    else if ((strcmp(argv[i],"-h")==0)||(strcmp(argv[i],"--help")==0)) {
              
      fprintf(stdout,"\nUsage: Convertfits [options], where options are:\n\n");
      fprintf(stdout," -h, --help\t\tboolean\tPrint this message\n");
      fprintf(stdout," -i, --inputfile\tstring\tThe input FITS file name [Default = %s]\n",uvar->inputfile);
      fprintf(stdout," -o, --outputdir\tstring\tThe output frame file directory [Default = %s]\n",uvar->outputdir);
      fprintf(stdout," -t, --deltat\t\tfloat\tThe sampling rate [Default = %e]\n",uvar->deltat);
      fprintf(stdout," -b, --bary\t\tboolean\tFlag for outputting barycentered data [Default = %d]\n",uvar->bary);
      fprintf(stdout," -d, --debug\t\tint\tFlag for outputting debug information [Default = %d]\n",debug);
      fprintf(stdout,"\n");
      exit(0);
    }
    i = i+2;
  }
  
  /* put clargs into string and store in general params */
  strcpy(clargs,"");
  for (i=0;i<argc;i++) {
    strcat(clargs,argv[i]);
    strcat(clargs," ");
  }

  return;

}
/*************************************************************************************************
 * 
 * this function outputs a shortinttimeseries to a frame file or files
 * 
 *************************************************************************************************/
void shortinttimeseriestoframes(shortinttimeseries *ts,char *outputdir)
{

  struct FrFile *frOutFile  = NULL;        /* frame file structure */
  struct FrameH *outFrame   = NULL;        /* frame data structure */  
  UINT2TimeSeries *output = NULL;           /* timeseries structure */
  char outputfile[256];                    /* stores name of output frame file */
  LIGOTimeGPS epoch;                       /* stores the timeseries epoch */
  INT4 i;                                  /* counter */

  /* convert start epoch to a GPS structure */
  epoch.gpsSeconds = (INT4)floor(ts->tstart);
  epoch.gpsNanoSeconds = (INT4)1e9*(ts->tstart - epoch.gpsSeconds);
  if (debug) fprintf(stdout,"STATUS : outputting file with epoch %d %d\n",epoch.gpsSeconds,epoch.gpsNanoSeconds); 

  /* constrct file name */
  sprintf(outputfile,"%s/XTE-SCOX1-%.0f-%.0f.gwf",outputdir,ts->tstart,ts->T);
  if (debug) fprintf(stdout,"STATUS : outputting to file name %s\n",outputfile);

  /* create empty timeseries */
  output = XLALCreateUINT2TimeSeries("XTE_PCA_SCOX1",&epoch,0,ts->deltat,&lalDimensionlessUnit,ts->N);
  if (debug) fprintf(stdout,"STATUS : allocated memory for temporary timeseries\n");
		    
  /* fill in timeseries */
  for (i=0;i<ts->N;i++) output->data->data[i] = ts->data[i];
  if (debug) fprintf(stdout,"STATUS : populated temporary timeseries\n");

  /* write out to frame data structure */
  outFrame = fr_add_proc_UINT2TimeSeries(outFrame,output,"ct","XTE_PCA"); 
  if (debug) fprintf(stdout,"STATUS : added timeseries to frame structure\n");

  /* open up a file for outputting to - second arg is compression level */
  frOutFile = FrFileONew(outputfile,0);

  /* write frame to file */
  FrameWrite(outFrame,frOutFile);
  
  /* close frame file */
  FrFileOEnd(frOutFile);
 
  /* free timeseries */
  XLALDestroyUINT2TimeSeries(output);
 
  return;
 
}
/*************************************************************************************************
 * 
 * this function barycenters a short int timeseries using a set of barycentric timestamps
 * 
 *************************************************************************************************/
void barycentershortinttimeseries(shortinttimeseries *ts,barycenteredtimestamps *stamps)
{
  
  long int i;
  long int N;
  shortinttimeseries tempts;
  double Ts,Te;
  double start,end;

  /* interpolation for the timestamps */
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,stamps->N);
  
  /* select latest start time and earliest end time that are within the timestamps and data limits */
  Ts = ts->tstart > stamps->dettime[0] ? ts->tstart : stamps->dettime[0];
  Te = (ts->tstart + ts->T) > stamps->dettime[stamps->N-1] ? (ts->tstart + ts->T) : stamps->dettime[stamps->N-1];

  /* check that times are still valid */
  if (Ts>=Te) {
    fprintf(stdout,"ERROR : GPS start is after GPS end time, not generating a barycentric time series.  Exiting.\n");
    exit(1);
  }

  /* set up gsl interpolation for det -> bary to get start and end times of barycentric timeseries */
  gsl_spline_init(spline,stamps->dettime,stamps->barytime,stamps->N);

  /* compute barycentric time of first and last sample */ 
  start = gsl_spline_eval(spline,Ts,acc);
  end = gsl_spline_eval(spline,Te,acc);

  /* redefine start and end to be integer GPS seconds */
  start = ceil(start);
  end = floor(end);
  if (start==end) {
    fprintf(stdout,"ERROR : integer GPS start equals integer GPS end time, not generating a barycentric time series.  Exiting.\n");
    exit(1);
  }
  if (debug) fprintf(stdout,"STATUS : barycentric GPS start and end times : %f %f -> span = %f\n",start,end,end-start);

  /* define number of barycentric timesamples */
  N = (long int)((end-start)/ts->deltat);
  
  /* allocate temp memory first */
  createshortinttimeseries(&tempts,N);
  if (debug) fprintf(stdout,"STATUS : created a temporary short int timeseries with %ld points\n",N);

  /* set up gsl interpolation for bary -> det for actual resampling */
  gsl_spline_init(spline,stamps->barytime,stamps->dettime,stamps->N);
  if (debug) fprintf(stdout,"STATUS : initialised the spline interpolation\n");

  /* loop over each evenly spaced time sample in barycentric frame */
  for (i=0;i<tempts.N;i++) {
    
    double bt = start + (double)i*ts->deltat;
    double delta_dettime = gsl_spline_eval(spline,bt,acc) - ts->tstart;
    long int idx = floor(0.5 + delta_dettime/ts->deltat);
  /*   printf("barytime = %6.12f delta_dettime = %6.12f index = %ld N = %ld\n",bt,delta_dettime,index,N); */
    tempts.data[i] = ts->data[idx];
    tempts.undefined[i] = ts->undefined[idx];

  }
  if (debug) fprintf(stdout,"STATUS : performed barycentering using nearest bin interpolation\n");

  /* reallocate input time series and fill it with barycentric timeseries and params */
  reallocshortinttimeseries(ts,tempts.N);
  for (i=0;i<tempts.N;i++) {
    ts->data[i] = tempts.data[i];
    ts->undefined[i] = tempts.undefined[i];
  }
  ts->N = tempts.N;
  ts->T = ts->N*ts->deltat;
  ts->tstart = start;
  if (debug) fprintf(stdout,"STATUS : resized original timeseries and replaced it with barycentered timeseries\n");

  /* free temp timeseries */
  freeshortinttimeseries(&tempts);
  if (debug) fprintf(stdout,"STATUS : freed temporary memory\n");

  /* free interpolation structures */
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  if (debug) fprintf(stdout,"STATUS : freed interpolation structures\n");

  {
    for (i=0;i<10;i++) {
      if (debug) fprintf(stdout,"STATUS : barycentered time series %6.12f %d (%d)\n",start+(double)i*ts->deltat,ts->data[i],ts->undefined[i]);
    }
    if (debug) fprintf(stdout,"STATUS : ...\n");
    for (i=ts->N-10;i<ts->N;i++) {
      if (debug) fprintf(stdout,"STATUS : barycentered time series %6.12f %d (%d)\n",start+(double)i*ts->deltat,ts->data[i],ts->undefined[i]);
    } 
  }

  return;

}
/*************************************************************************************************/
/*                                                                                               */
/* This function reads in data from a FITS file and returns a FITSdata data structure. */
/* If the dataflag is not set then it returns only header information otherwise it also */
/* returns the data from the specified column. */
/*                                                                                               */
/*************************************************************************************************/
void readfitsfile(char *filename,FITSdata *fitsfiledata,int col)
{

  fitsfile *fptr;                                    /* pointer to the FITS file, defined in fitsio.h */
  int status;

  /* initialise status */
  status = 0;

  /* initialise pointers */
  fitsfiledata->data = NULL;
  fitsfiledata->undefined = NULL;
  fitsfiledata->events = NULL;
  fitsfiledata->gtistart = NULL;
  fitsfiledata->gtiend = NULL;
  fitsfiledata->dettime = NULL;
  fitsfiledata->barytime = NULL;

  /* open the fits file */
  if (fits_open_file(&fptr,filename,READONLY,&status)) printFITSerror(status);

  /* add file name to the header information */
  sprintf(fitsfiledata->file,"%s",filename);
  
  /* read the header information */
  readFITSheader(fptr,fitsfiledata,status);
  
  /* read GTI information */
  readFITSgti(fptr,fitsfiledata,status);
    
  /* read the timestamps information */
  readFITStimestamps(fptr,fitsfiledata,status);
  
  /* read in the data */
  if (!fitsfiledata->event) readFITSbinneddata(fptr,fitsfiledata,col,status);
  else readFITSeventdata(fptr,fitsfiledata,col,status);
  
  /* close the file */
  if (fits_close_file(fptr,&status)) printFITSerror( status );
 
  return;

}
/*************************************************************************************************
 *
 * check GTI for good time window data.  The second extension HDU2 (3 from the beginning)
 * lists the standard good time intervals, i.e. the start and
 * stop times of all the potentially useable data in the file 
 *                                                                                               
 *************************************************************************************************/
void readFITSgti(fitsfile *fptr, FITSdata *fitsfiledata, int status)
{

  int i;                /* counter */
  int hdutype;
  double doublenull = 0.0; 
  int anynull;

  /* move to the gti header */ 
  if (fits_movabs_hdu(fptr,3,&hdutype,&status)) printFITSerror(status);

  /* get number of columns in the gti table - should be 2, exit if not. */
  {
    int ngticols = 0;
    if (fits_get_num_cols(fptr,&ngticols,&status)) printFITSerror(status);
    if (ngticols!=2) {
      fprintf(stderr,"ERROR : found %d columns in GTI table, should be 2.  Exiting.\n",ngticols);
      exit(1);
    }
  }
  
  /* get the number of rows in the gti table - this is the number of good segments in the file */
  if (fits_get_num_rows(fptr,&(fitsfiledata->ngtirows),&status)) printFITSerror(status);
  if (fitsfiledata->ngtirows<1) {
    fprintf(stderr,"ERROR : found %ld rows in GTI table, should be >0.  Exiting.\n", fitsfiledata->ngtirows);
    exit(1);
  }
  if (debug) fprintf(stdout,"STATUS : found %ld rows in GTI table.\n",fitsfiledata->ngtirows);
  
  /* allocate memory for the gti info */
  if ((fitsfiledata->gtistart = (double *)calloc(fitsfiledata->ngtirows,sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR : unable to allocate memory for GTI start time data.  Exiting.\n");
    exit(1);
  }
  if ((fitsfiledata->gtiend = (double *)calloc(fitsfiledata->ngtirows,sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR : unable to allocate memory for GTI end time data.  Exiting.\n");
    exit(1);
  }
  
  /*  read the gti start and end values  */  
  {
    if (fits_read_col(fptr,TDOUBLE,1,1,1,fitsfiledata->ngtirows,&doublenull,fitsfiledata->gtistart,&anynull,&status)) printFITSerror(status);
    if (fits_read_col(fptr,TDOUBLE,2,1,1,fitsfiledata->ngtirows,&doublenull,fitsfiledata->gtiend,&anynull,&status)) printFITSerror(status);
    for (i=0;i<fitsfiledata->ngtirows;i++) {
      if (debug) fprintf(stdout,"STATUS : found GTI range %6.12f -> %6.12f.\n",fitsfiledata->gtistart[i],fitsfiledata->gtiend[i]);
    }
  }

  /* IMPORTANT : convert to GPS using the offset value */
  for (i=0;i<fitsfiledata->ngtirows;i++) {
    fitsfiledata->gtistart[i] += fitsfiledata->gpsoffset;
    fitsfiledata->gtiend[i] += fitsfiledata->gpsoffset;
  }

  return;

}
/*************************************************************************************************
 *
 * read in first extension header info
 * The first extension Header HDU1 (2 from the beginning) contains keywords which
 * provide a complete and detailed description of the contents of the first
 * extension. For convenience, it also contains some of the same information as
 * the primary header.
 *
 ************************************************************************************************/
void readFITSheader(fitsfile *fptr, FITSdata *fitsfiledata, int status) 
{

  int hdutype;
  tddesparams energytimeparams;
  char comment[80];
  double tzero;
  char type[256];

  /* move to the first extension header */
  if (fits_movabs_hdu(fptr,2,&hdutype,&status)) printFITSerror(status);
   
  /* check object name and store in output structure */
  if (fits_read_key(fptr,TSTRING,string_OBJECT,&(fitsfiledata->objectname),comment,&status)) printFITSerror(status);
  if (debug) fprintf(stdout,"STATUS : checked object name is %s.\n",fitsfiledata->objectname);

  /* read in sky position */
  if (fits_read_key(fptr,TDOUBLE,string_RA_NOM,&(fitsfiledata->ra),comment,&status)) printFITSerror(status);
  if (fits_read_key(fptr,TDOUBLE,string_DEC_NOM,&(fitsfiledata->dec),comment,&status)) printFITSerror(status);
  if (debug) fprintf(stdout,"STATUS : read ra = %6.12f dec = %6.12f.\n",fitsfiledata->ra,fitsfiledata->dec);

  /* get the number of rows in the data table */
  if (fits_get_num_rows(fptr,&(fitsfiledata->nrows),&status)) printFITSerror(status);
  if (debug) fprintf(stdout,"STATUS : found %ld rows in first extension.\n",fitsfiledata->nrows);

  /* get number of columns in the data table */
  if (fits_read_key(fptr,TINT,string_TFIELDS,&(fitsfiledata->ncols),comment,&status)) printFITSerror(status);
  if (debug) fprintf(stdout,"STATUS : found %ld fields in first extension.\n",fitsfiledata->ncols);

  /* read the number of and size of the dimensions of table column 2 (the data column) */
  if (fits_read_tdim(fptr,2,1,&(fitsfiledata->nchannels),&(fitsfiledata->channelsize),&status)) printFITSerror(status);
  fitsfiledata->rowlength = fitsfiledata->channelsize*fitsfiledata->nchannels;
  if (debug) fprintf(stdout,"STATUS : read channelsize as %ld\n",fitsfiledata->channelsize);
  if (debug) fprintf(stdout,"STATUS : read nchannels as %d\n",fitsfiledata->nchannels);
  if (debug) fprintf(stdout,"STATUS : read rowlength as %ld\n",fitsfiledata->rowlength);
  
  /* get number of bytes per row in the data table */
  if (fits_read_key(fptr,TLONG,string_NAXIS1,&(fitsfiledata->nbytesperrow),comment,&status)) printFITSerror(status);
  if (debug) fprintf(stdout,"STATUS : read the number of bytes per row as %ld\n",fitsfiledata->nbytesperrow);

  /* get TIMEZERO for this file */
  /* Timestamps plus the value of TIMEZERO provide the time in TAI as seconds since January 1, 1994, at 00:00:28 (TAI). */
  if (fits_read_key(fptr,TDOUBLE,string_TIMEZERO,&tzero,comment,&status)) printFITSerror(status);
  if (debug) fprintf(stdout,"STATUS : read tzero as %6.12f\n",tzero);

  /* convert tzero to GPS - this is simply the TAI time plus the GPS-TAI constant */
  fitsfiledata->gpsoffset = tzero + (double)GPSMINUSTAI;  
  if (debug) fprintf(stdout,"STATUS : converted tzero to gps -> %6.12f\n",fitsfiledata->gpsoffset);
  
  /* get DELTAT for this file */
  /* DELTAT gives the the time between each row in the file - these are typically 1 or 2 seconds in length */
  if (fits_read_key(fptr,TDOUBLE,string_DELTAT,&(fitsfiledata->dttimestamp),comment,&status)) printFITSerror(status); 
  if (debug) fprintf(stdout,"STATUS : read dttimestamp as  %6.12f\n",fitsfiledata->dttimestamp);

  /* get mode string defining the detector configuration e.g. B_250us_2A_0_17_Q */
  if (fits_read_key(fptr,TSTRING,string_DATAMODE,&(fitsfiledata->mode),comment,&status)) printFITSerror(status); 
  if (debug) fprintf(stdout,"STATUS : read data mode as %s\n",fitsfiledata->mode);

  /* if the data mode is equal to 2LLD then we record this in the LLD (Lower Level Discriminator) flag */
  if (strstr(fitsfiledata->mode,string_2LLD)!=NULL) fitsfiledata->LLD = 1;
  else fitsfiledata->LLD = 0;
  if (debug) fprintf(stdout,"STATUS : set LLD flag = %d\n",fitsfiledata->LLD);
  
  /* get the type SA or SE */
  /* this tells us whether it is binned (science-array) or event (science-event) data */
  if (fits_read_key(fptr,TSTRING,string_HDUCLAS1,&type,comment,&status)) printFITSerror(status);
  if (strcmp(type,"ARRAY")==0) fitsfiledata->event = 0;
  else if (strcmp(type,"EVENTS")==0) fitsfiledata->event = 1;
  else {
    fprintf(stderr,"ERROR : data type \"%s\" not recognised.  Exiting\n",type);
    exit(1);
  }
  if (debug) fprintf(stdout,"STATUS : data type is %s\n",type);

  /* get the TDDES2 keyword - describes the data in the second column */
  /* this tells us the sampling time and the energy channels used */
  /* A typical string of DDL (Data Description Language) is a concatenation of tokens */
  /* (denoted by single letters) with assigned values (enclosed in square brackets). */
  {
    char *tddes2;
    if (fits_read_key_longstr(fptr,string_TDDES2,&tddes2,comment,&status)) printFITSerror(status);
    if (debug) fprintf(stdout,"STATUS : read TDDES2 as %s\n",tddes2);
  
    /* extract energy and sampling time info from the TDDES2 DDL string */
    converttddes(tddes2,&energytimeparams);
    fitsfiledata->minenergy = energytimeparams.minenergy;
    fitsfiledata->maxenergy = energytimeparams.maxenergy;
    if (debug) {
      fprintf(stdout,"STATUS : energy range extracted as %d - %d\n",fitsfiledata->minenergy,fitsfiledata->maxenergy);
      if (energytimeparams.deltat) fprintf(stdout,"STATUS : read time resolution as %6.12e\n",energytimeparams.deltat);
      else fprintf(stdout,"STATUS : could not extract time resolution from TDDES2 DDL string\n");
    }
    
    /* free string mem */
    free(tddes2);
   
  }

  /* check that sampling time is consistent */
  {
    double deltatcheck = 0.0;
    double timedel = 0.0;
    if (fits_read_key(fptr,TDOUBLE,string_TIMEDEL,&timedel,comment,&status)) printFITSerror(status);
    deltatcheck = timedel*fitsfiledata->nchannels/(fitsfiledata->rowlength);
    if (energytimeparams.deltat) {
      if (fabs(deltatcheck-energytimeparams.deltat)>1e-16) {
	fprintf(stderr,"ERROR : sampling time not consistent %6.12e != %6.12e\n",deltatcheck,fitsfiledata->deltat);
	exit(1);
      }
    }
    fitsfiledata->deltat = deltatcheck;
  }

}
/*************************************************************************************************
 *
 * read in first extension header info
 * The first extension Header HDU1 (2 from the beginning) contains keywords which
 * provide a complete and detailed description of the contents of the first
 * extension. For convenience, it also contains some of the same information as
 * the primary header.
 *
 ************************************************************************************************/
void readFITStimestamps(fitsfile *fptr, FITSdata *fitsfiledata, int status) 
{

  long int j;
  char comment[80];
  double doublenull = 0.0; 
  int anynull;
  int hdutype;
  
  /* move to the first extension header */
  if (fits_movabs_hdu(fptr,2,&hdutype,&status)) printFITSerror(status);
  
  /* make sure the first column is the time column otherwise exit */
  {
    char timestring[256];
    if (fits_read_key(fptr,TSTRING,string_TTYPE1,&timestring,comment,&status)) printFITSerror(status); 
    if (strncmp(timestring,"Time",4)!=0) {
      fprintf(stderr,"ERROR : first column label is \"%s\" NOT \"Time\".  Exiting.\n", timestring);
      exit(1);
    }
  }
  
  /* allocate mem for storing timestamps - there is a timestamp for each row NOT for each data value */
  if ((fitsfiledata->dettime = (double *)calloc(fitsfiledata->nrows,sizeof(double))) == NULL ) {
    fprintf(stderr,"ERROR : unable to allocate memory for detector timestamps.  Exiting.\n");
    exit(1);
  }

  /*  read the time detector frame timestamps from column 1 */  
  {
     int i;
     if (fits_read_col(fptr,TDOUBLE,1,1,1,fitsfiledata->nrows,&doublenull,fitsfiledata->dettime,&anynull,&status)) printFITSerror(status);  
     for (i=0;i<10;i++) {
       if (debug) fprintf(stdout,"STATUS : read time values as %6.12f\n",fitsfiledata->dettime[i]);
     }
  }
  
  /* make sure the final column is the time column otherwise exit */
  {
    char timestring[256];
    char keyword[8];
    snprintf(keyword,7,"TTYPE%ld",fitsfiledata->ncols);
    if (fits_read_key(fptr,TSTRING,keyword,&timestring,comment,&status)) printFITSerror(status); 
    if (strncmp(timestring,"BARYTIME",8)!=0) {
      fitsfiledata->bary = 0;
      fprintf(stdout,"WARNING : final column label is \"%s\" NOT \"BARYTIME\".\n", timestring);
    }
    else fitsfiledata->bary = 1;
  }
  
  /* IMPORTANT : convert to GPS using the offset value */
  for (j=0;j<fitsfiledata->nrows;j++) fitsfiledata->dettime[j] += fitsfiledata->gpsoffset;
  
  /* if the file contains barycentered times */
  if (fitsfiledata->bary) {
    
    /* allocate mem for storing timestamps - there is a timestamp for each row NOT for each data value */
    if ((fitsfiledata->barytime = (double *)calloc(fitsfiledata->nrows,sizeof(double))) == NULL) {
      fprintf(stderr,"ERROR : unable to allocate memory for barycentric timestamps.  Exiting.\n");
      exit(1);
    }
  
    /*  read the barycentered time values from final column */  
    {
      int i;
      if (fits_read_col(fptr,TDOUBLE,fitsfiledata->ncols,1,1,fitsfiledata->nrows,&doublenull,fitsfiledata->barytime,&anynull,&status)) printFITSerror(status);
      for (i=0;i<10;i++) {
	if (debug) fprintf(stdout,"STATUS : read barycentered time values as %6.12f\n",fitsfiledata->barytime[i]);
      }
    }

    /* IMPORTANT : convert to GPS using the offset value */
    for (j=0;j<fitsfiledata->nrows;j++) fitsfiledata->barytime[j] += fitsfiledata->gpsoffset;

  }

  /* fill start time and duration based on detector times */
  fitsfiledata->tstart = fitsfiledata->dettime[0];     /* start time */
  fitsfiledata->T = fitsfiledata->dettime[fitsfiledata->nrows-1] - fitsfiledata->dettime[0] + fitsfiledata->dttimestamp;    /* time span */
  if (debug) fprintf(stdout,"STATUS : determined start time as %6.12f and duration as %6.12f\n",fitsfiledata->tstart,fitsfiledata->T);

  return;

}
/*************************************************************************************************
 *
 * read in event data.
 *
 ************************************************************************************************/
void readFITSeventdata(fitsfile *fptr, FITSdata *fitsfiledata, int col, int status)
{
  
  char comment[80];
  long int nbytes;
  int hdutype;

  /* move to the first extension header */
  if (fits_movabs_hdu(fptr,2,&hdutype,&status)) printFITSerror(status);

  /* get sampling time - this is the actual sampling time interval in seconds */
  if (fits_read_key(fptr,TDOUBLE,string_TIMEDEL,&(fitsfiledata->deltat),comment,&status)) printFITSerror(status);
  if (debug) fprintf(stdout,"STATUS : read event data dt as %6.12e\n",fitsfiledata->deltat);
  
  /* define number of elements to read in (each of size char) - this is the total number of expected data values */
  fitsfiledata->Nevents = fitsfiledata->nrows;
  nbytes = (long)(fitsfiledata->rowlength*fitsfiledata->nrows*fitsfiledata->nchannels);
  if (debug) fprintf(stdout,"STATUS : preparing to read in %ld events\n",fitsfiledata->Nevents);
  
  /* allocate mem for data and data undefined flag */
  if ((fitsfiledata->events = (char *)calloc(nbytes,sizeof(char))) == NULL ) {
    fprintf(stderr,"ERROR : unable to allocate memory for event data.  Exiting.\n");
    exit(1);
  }
  if (debug) fprintf(stdout,"STATUS : allocated memory for nbytes*sizeof(char) = %ld\n",fitsfiledata->Nevents);
  
  /* read the complete data set - we must remember that any event with the most significant bit = 0 is not a real event !! */
  /* IMPORTANT : we read in rowlength chars for each event and the first char in each row defines whether the event is real */ 
  if (fits_read_col_bit(fptr,col,1,1,nbytes,fitsfiledata->events,&status)) printFITSerror(status);
  if (debug) fprintf(stdout,"STATUS : read in %ld events\n",fitsfiledata->Nevents);
  {
    int i;
    for (i=0;i<10;i++) {
      if (debug) fprintf(stdout,"STATUS : read data events as %6.12f (%6.12f) %d\n",fitsfiledata->dettime[i],fitsfiledata->barytime[i],fitsfiledata->events[fitsfiledata->rowlength*i]);
    }
    if (debug) fprintf(stdout,"STATUS : ...\n");
    for (i=fitsfiledata->Nevents-10;i<fitsfiledata->Nevents;i++) {
      if (debug) fprintf(stdout,"STATUS : read data events as %6.12f (%6.12f) %d\n",fitsfiledata->dettime[i],fitsfiledata->barytime[i],fitsfiledata->events[fitsfiledata->rowlength*i]);
    }
    
  }

  return;
 
}
/*************************************************************************************************
 *
 * read in binned data.
 *
 ************************************************************************************************/
void readFITSbinneddata(fitsfile *fptr, FITSdata *fitsfiledata, int col, int status)
{

  int anynull;
  int hdutype;

  /* move to the first extension header */
  if (fits_movabs_hdu(fptr,2,&hdutype,&status)) printFITSerror(status);
  
  /* define number of elements to read in - this is the total number of expected data values */
  fitsfiledata->N = (long)(fitsfiledata->rowlength*fitsfiledata->nrows);
  if (debug) fprintf(stdout,"STATUS : computed total number of expected data samples as %ld\n",fitsfiledata->N);
  
  /* check if requested column value is valid */
  if (col>fitsfiledata->ncols) {
    fprintf(stderr,"ERROR : tried to read more columns than present in file,  Exiting.\n");
    exit(1);
  }
  
  /* allocate mem for data and data undefined flag */
  if ((fitsfiledata->data = (short int *)calloc(fitsfiledata->N,sizeof(short int))) == NULL) {
    fprintf(stderr,"ERROR : unable to allocate memory for SCIENCE ARRAY data.  Exiting.\n");
    exit(1);
  }
  if ((fitsfiledata->undefined = (char *)calloc(fitsfiledata->N,sizeof(char))) == NULL) {
    fprintf(stderr,"ERROR : unable to allocate memory for SCIENCE ARRAY data quality.  Exiting.\n");
    exit(1);
  }
  if (debug) fprintf(stdout,"STATUS : allocated memory for data (approx %.0f MB)\n",floor(0.5 + fitsfiledata->N*2/1e6));
  
  /*  read the complete column data set */  
  if (fits_read_colnull(fptr,TSHORT,col,1,1,fitsfiledata->N,fitsfiledata->data,fitsfiledata->undefined,&anynull,&status)) printFITSerror( status );
  if (debug) fprintf(stdout,"STATUS : read in data from file %s\n",fitsfiledata->file);
  
  {
    int i;
    for (i=0;i<5;i++) {
      if (debug) fprintf(stdout,"STATUS : read data and undefined values as %d %d\n",fitsfiledata->data[i],fitsfiledata->undefined[i]);
    }
    if (debug) fprintf(stdout,"STATUS : ...\n");
    for (i=fitsfiledata->N-5;i<fitsfiledata->N;i++) {
      if (debug) fprintf(stdout,"STATUS : read data and undefined values as %d %d\n",fitsfiledata->data[i],fitsfiledata->undefined[i]);
    }
    
  }
  
  return;
 
}
/*************************************************************************************************/
/*                                                                                               */
/* extracts energy range and sampling parameters from tddes2 DDL string */
/*                                                                                               */
/*************************************************************************************************/
void converttddes(char *tddes, tddesparams *params)
{

  char out1[256], out2[256], new[256];
  char out3[256], out4[256], out5[256];
  char *tempstr1,*tempstr2,*tempstr3;
  int N1, N2, N3;
 
  /* initialise energy ranges */
  params->minenergy = 0;
  params->maxenergy = 0;
  params->deltat = 0.0;
  params->nsamples = 0;
  params->offset = 0.0;

  /* find start of energy section of string - exit if not found */
  if ((tempstr1 = strstr(tddes,"C["))==NULL) return;

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
  params->minenergy = atoi(out1);
  params->maxenergy = atoi(out2);

  /* find start of time section of string - exit if not found */
  if ((tempstr3 = strstr(tddes,"T["))==NULL) return;

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
  params->offset = atof(out3); 
  params->deltat = atof(out4);
  params->nsamples = atoi(out5); 

  return;

}
/*************************************************************************************************/
/*                                                                                               */
/* Print out cfitsio error messages and exit program */
/*                                                                                               */
/*************************************************************************************************/
void printFITSerror(int FITSstatus)
{
  
  if (FITSstatus)
    {
      fits_report_error(stderr, FITSstatus); /* print error report */
      
      exit( FITSstatus );    /* terminate the program, returning error status */
    }
  return;
}
/*************************************************************************************************/
/*                                                                                               */
/* this function converts a FITS event data file to binned for the given sampling rate.          */
/* It also reduces the events timestamps vector into a more coarsely sampled vector.             */
/*                                                                                               */
/*************************************************************************************************/
void eventdata_to_shortinttimeseries(FITSdata *fits,shortinttimeseries *ts,double dt)
{

  long int i;                            /* counter */
 
  /* check that requested sampling rate is at least as large as the original rate */
  if (dt<fits->eventdeltat) {
    fprintf(stderr,"ERROR : requested sampling time for event->binned conversion < original sampling time.  Exiting.\n");
    exit(1);
  }

  /* define complete span of time (remember event time is still in non GPS units) */
  /* define the time span in the detector frame from start of first bin to end of last bin */
  {
    double tempspan;
    ts->tstart = fits->dettime[0];                                  
    tempspan = fits->dettime[fits->nrows-1] - fits->tstart + fits->deltat;   
    ts->N = (long int)(tempspan/dt);
    ts->T = dt*ts->N;
    ts->deltat = dt;
  }
  if (debug) printf("STATUS : new binned timeseries parameters -> tstart = %6.12f T = %f dt = %e N = %ld\n",ts->tstart,ts->T,ts->deltat,ts->N);
  
  /* allocate memory for binned data */
  if ((ts->data = (short int *)calloc(ts->N,sizeof(short int))) == NULL) {
    fprintf(stderr,"ERROR : unable to allocate memory for binned SCIENCE EVENT data.  Exiting.\n");
    exit(1);
  }
  if ((ts->undefined = (char *)calloc(ts->N,sizeof(char))) == NULL) {
    fprintf(stderr,"ERROR : unable to allocate memory for binned SCIENCE EVENT data quality.  Exiting.\n");
    exit(1);
  }

  /* initialise undefined as zero - for event data all times are good except for those defined by the GTI table */
  /* time marker events are not real events and are not counted anyway so we treat those times as good */
  for (i=0;i<ts->N;i++) ts->undefined[i] = 0; 

  /* loop over each event and bin it if its not a time marker - ignore last count so we make things end on an integer second */
  /* also select coarsely spaced timestamps for later */
  for (i=0;i<fits->nrows-1;i++) {

    /* set counter tracking 1st bit for each event */
    long int k = i*fits->rowlength;
  
    /* define bin index corresponding to this event time with reference to the timestamp of the first event */
    long int idx = (long int)((fits->dettime[i]-fits->tstart)/dt);
   
    /* only add to bin if event was NOT a time marker */
    ts->data[idx] += (unsigned char)fits->events[k];

  }

  return;

}
/*************************************************************************************************
 * 
 * converts a fits file data structure containing science array (binned) data into a timeseries
 *
 *************************************************************************************************/
void arraydata_to_shortinttimeseries(FITSdata *fits, shortinttimeseries *ts,double dt)
{

  int i,j;                               /* counters */
  short int *temp_undefined = NULL;      /* temporary data quality vector */

  /* check that requested sampling rate is at least as large as the original rate */
  if (dt<fits->deltat) {
    fprintf(stderr,"ERROR : requested sampling time for binned conversion < original sampling time.  Exiting.\n");
    exit(1);
  }

  /* fill in output params */
  ts->tstart = fits->tstart;
  ts->T = fits->T;
  ts->deltat = fits->deltat;

  /* define number of elements in the timeseries - may be longer than the FITS data due to gaps */
  /* note that we use the user requested bin size here */
  ts->N = floor(0.5 + ts->T/dt);
  ts->deltat = ts->T/(double)ts->N;
  if (debug) fprintf(stdout,"STATUS : Defined new timeseries params : [Ts = %f T = %f dt = %6.12f N = %ld]\n",ts->tstart,ts->T,ts->deltat,ts->N);
  
  /* allocate mem for output timeseries */
  createshortinttimeseries(ts,ts->N);
  if (debug) fprintf(stdout,"STATUS : allocated memory for a timeseries\n");
  
  /* allocate memory for a temporary data quality vector */
  if ((temp_undefined = (short int *)calloc(ts->N,sizeof(short int))) == NULL) {
    fprintf(stderr,"ERROR : unable to allocate memory for temporary data quality data.  Exiting.\n");
    exit(1);
  }

  /* initialise as zero and all undefined - therefore all gaps will be automatically undefined */
  for (i=0;i<ts->N;i++) {
    ts->data[i] = 0;  
    ts->undefined[i] = 1; 
    temp_undefined[i] = 0;
  }

  /* loop over each timestamp and fill in timeseries */
  for (i=0;i<fits->nrows;i++) {
     
    /* loop over data within timestamp span */
    for (j=0;j<fits->rowlength;j++) {
      
      /* the index in the fits data in the j'th element in the i'th row */ 
      long int idxin = i*fits->rowlength + j;
      
      /* the index in the output timeseries of the fits data in the j'th element in the i'th row */ 
      long int idxout = floor(0.5 + (fits->dettime[i] - ts->tstart + j*fits->deltat)/ts->deltat);    

      /* add corresponding data to correct time output bin */
      ts->data[idxout] += fits->data[idxin];
      temp_undefined[idxout] += fits->undefined[idxin];
      
    }
    
  }

  /* make sure the undefined flags are either 0 or 1 - if greater than 1 then set to 1 */
  for (i=0;i<ts->N;i++) if (temp_undefined[i] > 1) ts->undefined[i] = 1;

  /* free memory */
  free(temp_undefined);

  return;

}
/*************************************************************************************************/
/*                                                                                               */
/* allocate memory for a shortinttimeseries                                                      */
/*                                                                                               */
/*************************************************************************************************/
void createshortinttimeseries(shortinttimeseries *x,long int N)
{

  x->N = N;
  x->data = (short int *)calloc(N,sizeof(short int));
  x->undefined = (char *)calloc(N,sizeof(char));
  
  return;
  
}
/************************************************************************************************* 
 *
 * free fits data structure
 *
 *************************************************************************************************/
void freeFITSdata(FITSdata *x)
{
 
  if (x->event == 1) free(x->events);
  if (x->data != NULL) free(x->data);
  if (x->undefined != NULL) free(x->undefined);
  if (x->dettime != NULL) free(x->dettime);
  if (x->barytime != NULL) free(x->barytime);
  if (x->gtistart != NULL) free(x->gtistart);
  if (x->gtiend != NULL) free(x->gtiend);
 
  return;
  
}
/************************************************************************************************* 
 *
 * free short int timeseries data structure
 *
 *************************************************************************************************/
void freeshortinttimeseries(shortinttimeseries *ts) 
{
  
  free(ts->data);
  free(ts->undefined);
 
  return;

}
/*************************************************************************************************
 *
 * here we set samples not in the GTI table as undefined and zero the corresponding data
 *
 ************************************************************************************************/
void applygtitable(shortinttimeseries *ts,FITSdata *fits)
{

  int i,k;
  long int newstartindex,newendindex,newN;
  double *BTIs = NULL;
  double *BTIe = NULL;

  /* allocate memory for bad time intervals BTIs - more useful than good time intervals for vetoing data */
  if ((BTIs = (double *)calloc(fits->ngtirows+1,sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR : failed to allocate memory for BTI start vector.  Exiting.\n");
    exit(1);
  }
  if ((BTIe = (double *)calloc(fits->ngtirows+1,sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR : failed to allocate memory for BTI end vector.  Exiting.\n");
    exit(1);
  }
  
  /* compute BTI (bad time interval) - including start end end gaps */
  BTIs[0] = ts->tstart;
  BTIe[0] = fits->gtistart[0];
  if (BTIe[0]<BTIs[0]) BTIs[0]=BTIe[0];  /* from rebinning the new start can be after the first gti start point */
  for (i=1;i<fits->ngtirows;i++) {
    BTIs[i] = fits->gtiend[i-1];
    BTIe[i] = fits->gtistart[i];
  }
  BTIs[fits->ngtirows] = fits->gtiend[fits->ngtirows-1];
  BTIe[fits->ngtirows] = ts->tstart + ts->T;
  if (BTIs[fits->ngtirows]>BTIe[fits->ngtirows]) BTIs[fits->ngtirows]=BTIe[fits->ngtirows];    /* from rebinning the new end can be before the last gti end point */
  for (i=0;i<fits->ngtirows+1;i++) fprintf(stdout,"STATUS : computed BTI : %f -> %f\n",BTIs[i],BTIe[i]);
  
  /* loop over each BTI */
  for (i=0;i<fits->ngtirows+1;i++) {

    /* define start and end indices for this BTI - make sure that indices are within valid bounds */
    long int sindex = (BTIs[i] - ts->tstart)/ts->deltat;
    long int eindex = (BTIe[i] - ts->tstart)/ts->deltat;
    if (sindex<0) sindex = 0;
    if (eindex>ts->N) eindex = ts->N;

    /* loop over the data between GTIs and set as undefined and zero the data */
    for (k=sindex;k<eindex;k++) {
      ts->undefined[k]=1;
      ts->data[k] = 0;
    }

  }
  if (debug) fprintf(stdout,"STATUS : zeroed and undefined non-GTI data in timeseries.\n");
  
  /* cut out start and end BTI data */
  newstartindex = (BTIe[0]-ts->tstart)/ts->deltat;
  newendindex = (BTIs[fits->ngtirows]-ts->tstart)/ts->deltat;
  newN = newendindex - newstartindex;

  /* if the new start index is moved then slide all of the data */
  if (newstartindex>0) {
    long int count = 0;
    for (i=newstartindex;i<newendindex;i++) {
      ts->data[count] = ts->data[i];
      ts->undefined[count] = ts->undefined[i];
      count++;
    }
  }
   
  /* resize timeseries vector */
  reallocshortinttimeseries(ts,newN);
  
  /* update timeseries params */
  ts->N = newN;
  ts->T = newN*ts->deltat;
  ts->tstart += newstartindex*ts->deltat;
  if (debug) fprintf(stdout,"STATUS : using the GTI we've changed start time and duration to %f and %f.\n",ts->tstart,ts->T);

  /* free the bad time interval data */
  free(BTIs);
  free(BTIe);
  
  return;

}
/************************************************************************************************* 
 *
 * resize a short int timeseries data structure
 *
 *************************************************************************************************/
void reallocshortinttimeseries(shortinttimeseries *ts,int N)
{

  ts->data = (short int *)realloc(ts->data,N*sizeof(short int));
  ts->undefined = (char *)realloc(ts->undefined,N*sizeof(char));
  ts->N = N;
  
  return;
  
}
/************************************************************************************************* 
 *
 * extract the timestamps from a FITS data structure
 *
 *************************************************************************************************/
void extractbarytimestamps(FITSdata *fits,barycenteredtimestamps *stamps)
{

  int i;       /* counter */
  int m = 1;   /* another count */

  /* pre-allocate memory for output timestamps */
  if ((stamps->dettime = (double *)calloc(fits->nrows,sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR : failed to allocate memory for detector frame timestamps.  Exiting.\n");
    exit(1);
  }
  if ((stamps->barytime = (double *)calloc(fits->nrows,sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR : failed to allocate memory for barycentered frame timestamps.  Exiting.\n");
    exit(1);
  }
  
  /* record first timestamp - we keep first and last points for maximum coverage */
  stamps->dettime[0] = fits->dettime[0];
  stamps->barytime[0] = fits->barytime[0];

  /* loop over each timestamp and record a timestamp every TIMESTAMP_DELTAT sec */
  for (i=1;i<fits->nrows-1;i++) {

    /* define next threshold timestamp value */
    double temp = fits->dettime[0] + m*TIMESTAMP_DELTAT;
    
    /* if we've crossed the threshold then record timestamps */
    if (fits->dettime[i]>temp) {
      stamps->dettime[m] = fits->dettime[i];
      stamps->barytime[m] = fits->barytime[i];
      m++;
    }

  }

  /* record last timestamp and new size of timestamps vector */
  stamps->dettime[m] = fits->dettime[fits->nrows-1];
  stamps->barytime[m] = fits->barytime[fits->nrows-1];
  stamps->N = m+1;

  /* resize stamps vectors */
  if ((stamps->dettime = (double *)realloc(stamps->dettime,stamps->N*sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR : failed to reallocate memory for detector frame timestamps.  Exiting.\n");
    exit(1);
  }
  if ((stamps->barytime = (double *)realloc(stamps->barytime,stamps->N*sizeof(double))) == NULL) {
    fprintf(stderr,"ERROR : failed to reallocate memory for barycentered frame timestamps.  Exiting.\n");
    exit(1);
  }

  if (debug) {
    for (i=0;i<5;i++) {
      if (debug) fprintf(stdout,"STATUS : extracted timestamps as %6.12f (%6.12f)\n",stamps->dettime[i],stamps->barytime[i]);
    }
    fprintf(stdout,"STATUS : ...\n");
    for (i=stamps->N-5;i<stamps->N;i++) {
      if (debug) fprintf(stdout,"STATUS : extracted timestamps as %6.12f (%6.12f)\n",stamps->dettime[i],stamps->barytime[i]);
    }
  }

  return;
  
}
