/*
*  Copyright (C) 2007 Gregory Mendell
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
/*                                                                               */
/* File: MakeSFTs.c                                                              */
/* Purpose: generate SFTs                                                        */
/* Origin: first written by Xavier Siemens, UWM - May 2005,                      */
/*         based on make_sfts.c by Bruce Allen                                   */
/* $Id: MakeSFTs.c,v 1.16 2007/06/27 01:46:48 gmendell Exp $                     */
/* Time domain cleaning procedure routines: Paola Leaci                          */
/*********************************************************************************/

/* REVISIONS: */
/* 11/02/05 gam; To save memory, change lalDebugLevel to 0; note when this is set to 3 that 3x the memory is used! */
/* 11/02/05 gam; To save memory, do not save the window function in memory; just recompute this for each SFT. */
/* 11/02/05 gam; To save memory, do not hold dataSingle.data and vtilde in memory at the same time. */
/* 11/03/05 gam; Add TRACKMEMUSE preprocessor flag and function for tracking memory usage; copied from make_sfts.c */
/* 11/19/05 gam; Reorganize code so that function appear in order called */
/* 11/19/05 gam; Add command line option to use single rather than double precision; will use LALDButterworthREAL4TimeSeries in single precision case. */
/* 11/19/05 gam; Rename vtilde fftDataDouble; add in fftDatasingle; rename fplan fftPlanDouble; add in fftPlanSingle */
/* 11/29/05 gam; Add PRINTEXAMPLEDATA to print the first NUMTOPRINT, middle NUMTOPRINT, and last NUMTOPRINT input/ouput data at various stages */
/* 12/27/05 gam; Add option make-gps-dirs, -D <num>, to make directory based on this many GPS digits. */
/* 12/28/05 gam; Add option misc-desc, -X <string> giving misc. part of the SFT description field in the filename */
/* 12/28/05 gam; Make FMIN = 48.0 and DF = 2000.0 global variables, add start-freq -F and band -B options to enter these */
/* 12/28/05 gam; Add in window-type options; 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window */
/* 12/28/05 gam; Add option --overlap-fraction -P (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 0.0) */
/* 12/28/05 gam; Add sft-version, -v option to select output SFT version (1 = default is version 1 SFTs; 2 = version 2 SFTs */
/* 12/28/05 gam; Add comment-field, -c option, for comment for version 2 SFTs */
/* 01/05/06 gam; Add in version 2 normalization; add function to print example version 2 SFT data; add memory checking */
/* 01/09/06 gam; Add make-tmp-file, -Z option; write SFT to .*.tmp file, then move to final file name. */
/* 01/10/07 gam; Add -u --frame-struct-type option; specified the input data type in the frames (default ADC_REAL4) */
/* 01/14/07 gam; Add -i --ifo option to specify the ifo independent of the channel name which can begin with H0, L0, or G0. */
/* 06/26/07 gam; Write all command line arguments to commentField of version 2 SFTs, based on /lalapps/src/calibration/ComputeStrainDriver.c */
/* 06/26/07 gam; Use finite to check that data does not contains a non-FINITE (+/- Inf, NaN) values, based on sftlib/SFTvalidate.c */

#include <config.h>
#if !defined HAVE_LIBGSL || !defined HAVE_LIBLALFRAME
#include <stdio.h>
int main(void) {fputs("disabled, no gsl or no lal frame library support.\n", stderr);return 1;}
#else

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <glob.h>
#include <errno.h>
#include <getopt.h>
#include <stdarg.h>

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/AVFactories.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/Window.h>
#include <lal/Calibration.h>
#include <lal/LALConstants.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/SFTfileIO.h>

#include <XLALPSSInterface.h>

extern char *optarg;
extern int optind, opterr, optopt;

#if ! __USE_ISOC99
extern int finite(double);
extern long lroundf(float);
#endif

/* track memory usage under linux */
#define TRACKMEMUSE 0

/* print the first NUMTOPRINT, middle NUMTOPRINT, and last NUMTOPRINT input/ouput data at various stages */
#define PRINTEXAMPLEDATA 0
#define NUMTOPRINT       2

/* 06/26/07 gam; use finite to check that data does not contains a non-FINITE (+/- Inf, NaN) values */
#define CHECKFORINFINITEANDNANS 1

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

/* This repeatedly tries to re-open a file.  This is useful on a
   cluster because automount may fail or time out.  This allows a
   number of tries with a bit of rest in between.
*/
FILE* tryopen(char *name, char *mode){
  int count=0;
  FILE *fp;
  
  while (!(fp=fopen(name, mode))){
    fprintf(stderr,"Unable to open file %s in mode %s.  Will retry...\n", name, mode);
    if (count++<10)
      sleep(10);
    else
      exit(3);
  }
  return fp;
}

/***************************************************************************/

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL8 HPf;              /* High pass filtering frequency */
  INT4 T;                 /* SFT duration */
  char *stringT;          /* 12/27/05 gam; string with SFT duration */
  INT4 GPSStart;
  INT4 GPSEnd;
  INT4 makeGPSDirs;        /* 12/27/05 gam; add option to make directories based on gps time */
  INT4 sftVersion;         /* 12/28/05 gam; output SFT version */
  char *commentField;      /* 12/28/05 gam; string comment for version 2 SFTs */
  BOOLEAN htdata;          /* flag that indicates we're doing h(t) data */
  BOOLEAN makeTmpFile;     /* 01/09/06 gam */
  char *FrCacheFile;       /* Frame cache file */
  char *ChannelName;
  char *IFO;               /* 01/14/07 gam */
  char *SFTpath;           /* path to SFT file location */
  char *miscDesc;          /* 12/28/05 gam; string giving misc. part of the SFT description field in the filename */
  INT4 TDcleaningProc;	   /*1=YES and 0=NO*/
  REAL8 fc;                /*Cut frequency for the bilateral highpass filter. It has to be used only if TDcleaningProc is YES.*/
  REAL8 cr;                /*Critical ratio threshold. It has to be used only if TDcleaningProc is YES.*/
  INT4 windowOption;       /* 12/28/05 gam; window options; 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window */
  REAL8 overlapFraction;   /* 12/28/05 gam; overlap fraction (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 1.0). */
  BOOLEAN useSingle;       /* 11/19/05 gam; use single rather than double precision */
  char *frameStructType;   /* 01/10/07 gam */
} CommandLineArgs;

struct headertag {
  REAL8 endian;
  INT4  gps_sec;
  INT4  gps_nsec;
  REAL8 tbase;
  INT4  firstfreqindex;
  INT4  nsamples;
} header;

typedef struct BilHP{
  REAL8 fcut;
  REAL8 crt;
} BilHP;

typedef struct tagParamOfEvent{
  REAL4 tau;       /*Memory time of the autoregressive average*/   
  REAL4 deadtime;  /*Dead-time in seconds*/
  REAL4 factor;    /*e.g. 10: when re-evaluate mean and std*/
  INT4 iflev;      /*0=no events; 1=begin; 2=event; 3=end of the event*/
  REAL4 *xamed;    /*AUTOREGRESSIVE MEAN*/
  REAL4 *xastd;    /*AUTOREGRESSIVE STANDARD DEVIATION*/
  REAL4 *med;     /*AUTOREGRESSIVE MEAN*/
  REAL4 *std; 
  REAL4 xw;
  REAL4 qw;
  REAL8 w_norm;   /*Normalization factor for mean and std*/
  REAL4 edge;     /*Seconds around (before and after) the event have to be excluded*/
  int total;      /* Samples, due to all events, have been cleaned*/
  int *begin;     /*Sample numbers of the beginning of each event*/
  int *duration;  /*Event duration in number of samples*/
  int *imax;      /*Index at which each event has the maximum*/
  REAL4 *crmax;   /*Critical ratio of the maximum of the event*/
  REAL4 *ener;    /*Event energy (sum of the squared amplitudes)*/
  int number;     /*Number of events*/
  INT2 ar;        /*Algorithm for the autoregressive evaluation; it has to be = 1 */
 } ParamOfEvent;



/***************************************************************************/

/* GLOBAL VARIABLES */
REAL8 FMIN = 48.0; /* default start frequency; 0.0056 */
REAL8 DF = 2000.0; /* 2000.0 default band; 16383.9944 */

static LALStatus status;
INT4 lalDebugLevel=0;        /* 11/02/05 gam; changed from 3 to 0. */
FrCache *framecache;         /* frame reading variables */
FrStream *framestream=NULL;

REAL8TimeSeries dataDouble;
REAL4TimeSeries dataSingle;

REAL8TimeSeries dataDoubleFirstHP;
REAL4TimeSeries dataSingleFirstHP;
REAL8TimeSeries dataDoubleHP;

REAL8TimeSeries databil2;
REAL4TimeSeries databil2Single;

REAL4TimeSeries dataSingleHP;
REAL8TimeSeries dataDoubleClean;
REAL4TimeSeries dataSingleClean;

INT2TimeSeries dataINT2;
INT4TimeSeries dataINT4;
INT8TimeSeries dataINT8;
INT4 SegmentDuration;
LIGOTimeGPS gpsepoch;

/* REAL8Vector *window; */ /* 11/02/05 gam */

REAL8FFTPlan *fftPlanDouble;           /* fft plan and data container, double precision case */
COMPLEX16Vector *fftDataDouble = NULL;

REAL4FFTPlan *fftPlanSingle;           /* 11/19/05 gam; fft plan and data container, single precision case */
COMPLEX8Vector *fftDataSingle = NULL;

CHAR allargs[16384]; /* 06/26/07 gam; copy all command line args into commentField, based on /lalapps/src/calibration/ComputeStrainDriver.c */

int testpao = 0;

/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);

/* Allocates space for data time series */
int AllocateData(struct CommandLineArgsTag CLA);

/* Reads data */
int ReadData(struct CommandLineArgsTag CLA);

/* High passes data */
int HighPass(struct CommandLineArgsTag CLA);

/* Windows data */
int WindowData(struct CommandLineArgsTag CLA);
int WindowDataTukey2(struct CommandLineArgsTag CLA);
int WindowDataHann(struct CommandLineArgsTag CLA);

/*-----------------------------------------------------------*/
/* Useful functions for the time domain cleaning procedure */
/*-----------------------------------------------------------*/

/* Time Domain Cleaning Procedure */
int TDCleaning(struct CommandLineArgsTag CLA);
int TDCleaning_R4(struct CommandLineArgsTag CLA);
int TDCleaning_NoBil_R4(struct CommandLineArgsTag CLA);

/* Time Domain Cleaning with PSS functions */
int PSSTDCleaningDouble(struct CommandLineArgsTag CLA);
int PSSTDCleaningSingle(struct CommandLineArgsTag CLA);

/*Initialization of the event paramaters*/
int EventParamInit(long len, ParamOfEvent *even_param);

/*BILATERAL HIGHPASS FILTER: the SINGLE-PRECISION data are highpassed in the time domain.*/
int BilHighpass_dataSingle( REAL4TimeSeries *seriesHP, REAL4TimeSeries *firsthighpass, REAL4TimeSeries *series, ParamOfEvent *even_param, BilHP *myparams );

/*BILATERAL HIGHPASS FILTER: the DOUBLE-PRECISION data are highpassed in the time domain.*/
int BilHighpass_dataDouble( REAL8TimeSeries *seriesHP, REAL8TimeSeries *firsthighpass, REAL8TimeSeries *series, ParamOfEvent *even_param, BilHP *myparams );

/*COMPUTATION OF MEAN AND STANDARD DEVIATION for data in SINGLE-PRECISION: the calculations are made by means of an autoregressive procedure*/
int MedStd_dataSingle( REAL4TimeSeries *series, ParamOfEvent *even_param );

/*COMPUTATION OF MEAN AND STANDARD DEVIATION for data in DOUBLE-PRECISION: the calculations are made by means of an autoregressive procedure*/
int MedStd_dataDouble( REAL8TimeSeries *seriesHPbil2, REAL8TimeSeries *series, ParamOfEvent *even_param );

/*EVENT-HUNTING: this function looks for events with a critical ratio above a 'fixed' threshold. DATA IN SINGLE-PRECISION*/
int EventSearch_dataSingle( REAL4TimeSeries *series, ParamOfEvent *even_param, BilHP *myparams );

/*EVENT-HUNTING: this function looks for events with a critical ratio above a 'fixed' threshold. DATA IN DOUBLE-PRECISION*/
int EventSearch_dataDouble( REAL8TimeSeries *series, ParamOfEvent *even_param, BilHP *myparams );

/*ELIMINATION OF HIGH FREQUENCY SPIKES: this function identifies and removes large events from the highpassed time series. DATA IN SINGLE-PRECISION*/
int EventRemoval_dataSingle(REAL4TimeSeries *seriesCL, REAL4TimeSeries *series, REAL4TimeSeries *seriesHP, ParamOfEvent *even_param, BilHP *myparams);

/*ELIMINATION OF HIGH FREQUENCY SPIKES: this function identifies and removes large events from the highpassed time series. DATA IN DOUBLE-PRECISION*/
/*int EventRemoval_dataDouble(REAL8TimeSeries *seriesCL, REAL8TimeSeries *series, REAL8TimeSeries *seriesHP, ParamOfEvent *even_param, BilHP *myparams);*/
int EventRemoval_dataDouble(REAL8TimeSeries *seriesCL, REAL8TimeSeries *series, REAL8TimeSeries *seriesHP, REAL8TimeSeries *seriesHPbil2, ParamOfEvent *even_param, BilHP *myparams);
/*-----------------------------------------------------------*/


/* create an SFT */
int CreateSFT(struct CommandLineArgsTag CLA);

/* writes out an SFT */
int WriteSFT(struct CommandLineArgsTag CLA);
int WriteVersion2SFT(struct CommandLineArgsTag CLA);

/* Frees the memory */
int FreeMem(struct CommandLineArgsTag CLA);


/*****************************************paola********/
/* 12/27/05 gam; function to create sftDescField for output directory or filename. */
void getSFTDescField(CHAR *sftDescField, CHAR *numSFTs, CHAR *ifo, CHAR *stringT, CHAR *typeMisc) {
       strcpy(sftDescField,numSFTs);
       strcat(sftDescField, "_");
       strcat(sftDescField,ifo);
       strcat(sftDescField, "_");
       strcat(sftDescField,stringT);
       strcat(sftDescField, "SFT");
       if (typeMisc != NULL) {
          strcat(sftDescField, "_");
          strcat(sftDescField, typeMisc);
       }
}

/* 12/27/05 gam; concat to the sftPath the directory name based on GPS time; make this directory if it does not already exist */
void mkSFTDir(CHAR *sftPath, CHAR *site, CHAR *numSFTs, CHAR *ifo, CHAR *stringT, CHAR *typeMisc,CHAR *gpstime, INT4 numGPSdigits) {
     CHAR sftDescField[256];
     CHAR mkdirCommand[256];
     strcat(sftPath,"/");
     strcat(sftPath,site);
     strcat(sftPath,"-");
     getSFTDescField(sftDescField, numSFTs, ifo, stringT, typeMisc);
     strcat(sftPath,sftDescField);
     strcat(sftPath,"-");
     strncat(sftPath,gpstime,numGPSdigits);
     sprintf(mkdirCommand,"mkdir -p %s",sftPath);
     system(mkdirCommand);
}

/* 12/27/05 gam; make SFT file name according to LIGO T040164-01 specification */
void mkSFTFilename(CHAR *sftFilename, CHAR *site, CHAR *numSFTs, CHAR *ifo, CHAR *stringT, CHAR *typeMisc,CHAR *gpstime) {
     CHAR sftDescField[256];
     strcpy(sftFilename,site);
     strcat(sftFilename,"-");
     getSFTDescField(sftDescField, numSFTs, ifo, stringT, typeMisc);
     strcat(sftFilename,sftDescField);
     strcat(sftFilename,"-");
     strcat(sftFilename,gpstime);
     strcat(sftFilename,"-");
     strcat(sftFilename,stringT);
     strcat(sftFilename,".sft");
}

/* 01/09/06 gam; move filename1 to filename2 */
void mvFilenames(CHAR *filename1, CHAR *filename2) {
     CHAR mvFilenamesCommand[512];
     sprintf(mvFilenamesCommand,"mv %s %s",filename1,filename2);
     system(mvFilenamesCommand);
}

#if TRACKMEMUSE
void printmemuse() {
   pid_t mypid=getpid();
   char commandline[256];
   fflush(NULL);
   sprintf(commandline,"cat /proc/%d/status | /bin/grep Vm | /usr/bin/fmt -140 -u", (int)mypid);
   system(commandline);
   fflush(NULL);
 }
#endif

#if PRINTEXAMPLEDATA
void printExampleDataSingle() {
      INT4 i,dataLength;
      dataLength = dataSingle.data->length;
      fprintf(stdout,"dataSingle.deltaT, 1.0/dataSingle.deltaT = %23.16e, %23.16e\n",dataSingle.deltaT,1.0/dataSingle.deltaT);
      fprintf(stdout,"dataSingle.epoch.gpsSeconds,dataSingle.epoch.gpsNanoSeconds = %i, %i\n",dataSingle.epoch.gpsSeconds,dataSingle.epoch.gpsNanoSeconds);
      for (i=0;i<NUMTOPRINT;i++) {
        fprintf(stdout,"%i %23.16e\n",i,dataSingle.data->data[i]);
      }
      for (i=dataLength/2;i<dataLength/2+NUMTOPRINT;i++) {
        fprintf(stdout,"%i %23.16e\n",i,dataSingle.data->data[i]);
      }
      for (i=dataLength-NUMTOPRINT;i<dataLength;i++) {
        fprintf(stdout,"%i %23.16e\n",i,dataSingle.data->data[i]);
      }
      fflush(stdout);
}    

void printExampleDataDouble() {
      INT4 i,dataLength;
      dataLength = dataDouble.data->length;
      fprintf(stdout,"dataDouble.deltaT, 1.0/dataDouble.deltaT = %23.16e, %23.16e\n",dataDouble.deltaT,1.0/dataDouble.deltaT);
      fprintf(stdout,"dataDouble.epoch.gpsSeconds,dataDouble.epoch.gpsNanoSeconds = %i, %i\n",dataDouble.epoch.gpsSeconds,dataDouble.epoch.gpsNanoSeconds);
      for (i=0;i<NUMTOPRINT;i++) {
        fprintf(stdout,"%i %23.16e\n",i,dataDouble.data->data[i]);
      }
      for (i=dataLength/2;i<dataLength/2+NUMTOPRINT;i++) {
        fprintf(stdout,"%i %23.16e\n",i,dataDouble.data->data[i]);
      }
      for (i=dataLength-NUMTOPRINT;i<dataLength;i++) {
        fprintf(stdout,"%i %23.16e\n",i,dataDouble.data->data[i]);
      }
      fflush(stdout);
}

void printExampleFFTData(struct CommandLineArgsTag CLA)
{
  int firstbin=(INT4)(FMIN*CLA.T+0.5);
  int k;
  INT4 nsamples=(INT4)(DF*CLA.T+0.5);

  if(CLA.useSingle) {
    for (k=0; k<NUMTOPRINT; k++) {
      fprintf(stdout,"%i %23.16e %23.16e\n",k,fftDataSingle->data[k+firstbin].re,fftDataSingle->data[k+firstbin].im);
    }
    for (k=nsamples/2; k<nsamples/2+NUMTOPRINT; k++) {
      fprintf(stdout,"%i %23.16e %23.16e\n",k,fftDataSingle->data[k+firstbin].re,fftDataSingle->data[k+firstbin].im);
    }
    for (k=nsamples-NUMTOPRINT; k<nsamples; k++) {
      fprintf(stdout,"%i %23.16e %23.16e\n",k,fftDataSingle->data[k+firstbin].re,fftDataSingle->data[k+firstbin].im);
    }    
  } else {
    for (k=0; k<NUMTOPRINT; k++) {
      fprintf(stdout,"%i %23.16e %23.16e\n",k,fftDataDouble->data[k+firstbin].re,fftDataDouble->data[k+firstbin].im);
    }
    for (k=nsamples/2; k<nsamples/2+NUMTOPRINT; k++) {
      fprintf(stdout,"%i %23.16e %23.16e\n",k,fftDataDouble->data[k+firstbin].re,fftDataDouble->data[k+firstbin].im);
    }
    for (k=nsamples-NUMTOPRINT; k<nsamples; k++) {
      fprintf(stdout,"%i %23.16e %23.16e\n",k,fftDataDouble->data[k+firstbin].re,fftDataDouble->data[k+firstbin].im);
    }    
  }
  fflush(stdout);
}

void printExampleSFTDataGoingToFile(struct CommandLineArgsTag CLA)
{
  int firstbin=(INT4)(FMIN*CLA.T+0.5);
  int k;
  INT4 nsamples=(INT4)(DF*CLA.T+0.5);
  REAL4 rpw,ipw;

  if(CLA.useSingle) {
    for (k=0; k<NUMTOPRINT; k++) {
      rpw=((REAL4)(((REAL8)DF)/(0.5*(REAL8)(1/dataSingle.deltaT))))
	* fftDataSingle->data[k+firstbin].re;
      ipw=((REAL4)(((REAL8)DF)/(0.5*(REAL8)(1/dataSingle.deltaT))))
	* fftDataSingle->data[k+firstbin].im;
      fprintf(stdout,"%i %23.16e %23.16e\n",k,rpw,ipw);
    }
    for (k=nsamples/2; k<nsamples/2+NUMTOPRINT; k++) {
      rpw=((REAL4)(((REAL8)DF)/(0.5*(REAL8)(1/dataSingle.deltaT))))
	* fftDataSingle->data[k+firstbin].re;
      ipw=((REAL4)(((REAL8)DF)/(0.5*(REAL8)(1/dataSingle.deltaT))))
	* fftDataSingle->data[k+firstbin].im;
      fprintf(stdout,"%i %23.16e %23.16e\n",k,rpw,ipw);
    }
    for (k=nsamples-NUMTOPRINT; k<nsamples; k++) {
      rpw=((REAL4)(((REAL8)DF)/(0.5*(REAL8)(1/dataSingle.deltaT))))
	* fftDataSingle->data[k+firstbin].re;
      ipw=((REAL4)(((REAL8)DF)/(0.5*(REAL8)(1/dataSingle.deltaT))))
	* fftDataSingle->data[k+firstbin].im;
      fprintf(stdout,"%i %23.16e %23.16e\n",k,rpw,ipw);
    }    
  } else {
    for (k=0; k<NUMTOPRINT; k++) {
      rpw=(((REAL8)DF)/(0.5*(REAL8)(1/dataDouble.deltaT))) 
	* fftDataDouble->data[k+firstbin].re;
      ipw=(((REAL8)DF)/(0.5*(REAL8)(1/dataDouble.deltaT))) 
	* fftDataDouble->data[k+firstbin].im;
      fprintf(stdout,"%i %23.16e %23.16e\n",k,rpw,ipw);
    }
    for (k=nsamples/2; k<nsamples/2+NUMTOPRINT; k++) {
      rpw=(((REAL8)DF)/(0.5*(REAL8)(1/dataDouble.deltaT))) 
	* fftDataDouble->data[k+firstbin].re;
      ipw=(((REAL8)DF)/(0.5*(REAL8)(1/dataDouble.deltaT))) 
	* fftDataDouble->data[k+firstbin].im;
      fprintf(stdout,"%i %23.16e %23.16e\n",k,rpw,ipw);
    }
    for (k=nsamples-NUMTOPRINT; k<nsamples; k++) {
      rpw=(((REAL8)DF)/(0.5*(REAL8)(1/dataDouble.deltaT))) 
	* fftDataDouble->data[k+firstbin].re;
      ipw=(((REAL8)DF)/(0.5*(REAL8)(1/dataDouble.deltaT))) 
	* fftDataDouble->data[k+firstbin].im;
      fprintf(stdout,"%i %23.16e %23.16e\n",k,rpw,ipw);
    }    
  }
  fflush(stdout);
}

void printExampleVersion2SFTDataGoingToFile(struct CommandLineArgsTag CLA, SFTtype *oneSFT)
{
  int k;
  INT4 nsamples=(INT4)(DF*CLA.T+0.5);

  for (k=0; k<NUMTOPRINT; k++) {  
    fprintf(stdout,"%i %23.16e %23.16e\n",k,oneSFT->data->data[k].re,oneSFT->data->data[k].im);
  }
  for (k=nsamples/2; k<nsamples/2+NUMTOPRINT; k++) {
    fprintf(stdout,"%i %23.16e %23.16e\n",k,oneSFT->data->data[k].re,oneSFT->data->data[k].im);
  }
  for (k=nsamples-NUMTOPRINT; k<nsamples; k++) {
    fprintf(stdout,"%i %23.16e %23.16e\n",k,oneSFT->data->data[k].re,oneSFT->data->data[k].im);
  }    
}
#endif

int PrintREAL4ArrayToFile(char*name,REAL4*array,UINT4 length) {
  UINT4 i;
  FILE*fp=fopen(name,"w");
  if(!fp) {
    fprintf(stderr,"Could not open file '%s' for writing\n",name);
    return -1;
  } else {
    for(i=0;i<length;i++)
      fprintf(fp,"%23.16e\n",array[i]);
    fclose(fp);
  }
  return 0;
}


/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
  /* int j; */ /* 12/28/05 gam */

  #if TRACKMEMUSE
    printf("Memory use at startup is:\n"); printmemuse();
  #endif

  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;
  SegmentDuration = CommandLineArgs.GPSEnd - CommandLineArgs.GPSStart ;

  /* create Frame cache, open frame stream and delete frame cache */
  LALFrCacheImport(&status,&framecache,CommandLineArgs.FrCacheFile);
  TESTSTATUS( &status );
  LALFrCacheOpen(&status,&framestream,framecache);
  TESTSTATUS( &status );
  LALDestroyFrCache(&status,&framecache);
  TESTSTATUS( &status );

  #if TRACKMEMUSE
    printf("Memory use after reading command line arguments and reading frame cache:\n"); printmemuse();
  #endif

  if( SegmentDuration < CommandLineArgs.T)
    {
      fprintf(stderr,"Cannot fit an SFT of duration %d between %d and %d\n",
	      CommandLineArgs.T, CommandLineArgs.GPSStart, CommandLineArgs.GPSEnd );
      return 0;;
    }

  gpsepoch.gpsSeconds = CommandLineArgs.GPSStart;
  gpsepoch.gpsNanoSeconds = 0;

  /* Allocates space for data */
  if (AllocateData(CommandLineArgs)) return 2;

  /* for(j=0; j < SegmentDuration/CommandLineArgs.T; j++) */ /* 12/28/05 gam */
  while(gpsepoch.gpsSeconds + CommandLineArgs.T <= CommandLineArgs.GPSEnd)
    {

      /* gpsepoch.gpsSeconds = CommandLineArgs.GPSStart + j*CommandLineArgs.T;
      gpsepoch.gpsNanoSeconds = 0; */ /* 12/28/05 gam; now updated at the bottom of while loop */

      /* Reads T seconds of data */
      if (ReadData(CommandLineArgs)) return 3;

      /* High-pass data with Butterworth filter */
      if(HighPass(CommandLineArgs)) return 4;

      /* Window data; 12/28/05 gam; add options */
      if (CommandLineArgs.windowOption==1) {
         if(WindowData(CommandLineArgs)) return 5; /* CommandLineArgs.windowOption==1 is the default */
      } else if (CommandLineArgs.windowOption==2) {
         if(WindowDataTukey2(CommandLineArgs)) return 5;
      } else if (CommandLineArgs.windowOption==3) {
         if(WindowDataHann(CommandLineArgs)) return 5;
      } else {
        /* Continue with no windowing; parsing of command line args makes sure options are one of the above or 0 for now windowing. */
      }

      /* Time Domain cleaning procedure */
      if (CommandLineArgs.TDcleaningProc == 1) {
         if(TDCleaning(CommandLineArgs))
	   return 6;   
      } else if (CommandLineArgs.TDcleaningProc == 2) {
         if(PSSTDCleaningDouble(CommandLineArgs))
	   return 6;   
      } else if (CommandLineArgs.TDcleaningProc == 3) {
         if(TDCleaning_R4(CommandLineArgs))
	   return 6;   
      } else if (CommandLineArgs.TDcleaningProc == 4) {
         if(TDCleaning_NoBil_R4(CommandLineArgs))
	   return 6;   
      }
  
      /* create an SFT */
      if(CreateSFT(CommandLineArgs)) return 7;

      /* write out sft */
      if (CommandLineArgs.sftVersion==2) {
        if(WriteVersion2SFT(CommandLineArgs)) return 8;
      } else {
        if(WriteSFT(CommandLineArgs)) return 8;  /* default is to output version 1 SFTs */
      }

      gpsepoch.gpsSeconds = gpsepoch.gpsSeconds + (INT4)((1.0 - CommandLineArgs.overlapFraction)*((REAL8)CommandLineArgs.T));
      gpsepoch.gpsNanoSeconds = 0;
    }

  if(FreeMem(CommandLineArgs)) return 9;

  #if TRACKMEMUSE
        printf("Memory use after freeing all allocated memory:\n"); printmemuse();
  #endif

  return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/

/*** FUNCTIONS ***/

/*******************************************************************************/
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA)
{
  INT4 errflg=0;
  INT4 i;              /* 06/26/07 gam */
  struct option long_options[] = {
    {"high-pass-freq",       required_argument, NULL,          'f'},
    {"sft-duration",         required_argument, NULL,          't'},
    {"sft-write-path",       required_argument, NULL,          'p'},
    {"frame-cache",          required_argument, NULL,          'C'},
    {"channel-name",         required_argument, NULL,          'N'},
    {"gps-start-time",       required_argument, NULL,          's'},
    {"gps-end-time",         required_argument, NULL,          'e'},
    {"sft-version",          optional_argument, NULL,          'v'},
    {"comment-field",        optional_argument, NULL,          'c'},
    {"start-freq",           optional_argument, NULL,          'F'},
    {"band",                 optional_argument, NULL,          'B'},
    {"make-gps-dirs",        optional_argument, NULL,          'D'},
    {"make-tmp-file",        optional_argument, NULL,          'Z'},
    {"misc-desc",            optional_argument, NULL,          'X'},
    {"frame-struct-type",    optional_argument, NULL,          'u'},
    {"ifo",                  optional_argument, NULL,          'i'},
    {"TD-cleaning-proc",     optional_argument, NULL,          'a'},
    {"cut-frequency",        optional_argument, NULL,          'b'},
    {"crit-ratio-thrs",      optional_argument, NULL,          'r'},
    {"window-type",          optional_argument, NULL,          'w'},
    {"overlap-fraction",     optional_argument, NULL,          'P'},
    {"ht-data",              no_argument,       NULL,          'H'},
    {"use-single",           no_argument,       NULL,          'S'},
    {"help",                 no_argument,       NULL,          'h'},
    {0, 0, 0, 0}
  };
  char args[] = "hHZSf:t:C:N:i:s:e:v:c:F:B:D:a:b:r:X:u:w:P:p:";

  /* Initialize default values */
  CLA->HPf=-1.0;
  CLA->T=0.0;
  CLA->stringT=NULL;  /* 12/27/05 gam */
  CLA->FrCacheFile=NULL;
  CLA->GPSStart=0;
  CLA->GPSEnd=0;
  CLA->makeGPSDirs=0; /* 12/27/05 gam; add option to make directories based on gps time */
  CLA->sftVersion=1;  /* 12/28/05 gam; output SFT version; default is version 1 SFTs */
  CLA->ChannelName=NULL;
  CLA->IFO=NULL;        /* 01/14/07 gam */
  CLA->SFTpath=NULL;
  CLA->miscDesc=NULL;   /* 12/28/05 gam; misc. part of the SFT description field in the filename (also used if makeGPSDirs > 0) */
  CLA->commentField=NULL; /* 12/28/05 gam; comment for version 2 SFT header. */
  CLA->TDcleaningProc=0;  /*Time Domain cleaning procedure. It is zero by default.*/
  CLA->fc=0;              /*Cut frequency for the bilateral highpass filter. It has to be used only if TDcleaningProc is YES.*/
  CLA->cr=0;             /*Critical ratio threshold. It has to be used only if TDcleaningProc is YES.*/
  CLA->windowOption=1;  /* 12/28/05 gam; window options; 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window */
  CLA->overlapFraction=0.0; /* 12/28/05 gam; overlap fraction (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 0.0). */
  CLA->htdata = 0;
  CLA->makeTmpFile = 0; /* 01/09/06 gam */  
  CLA->useSingle = 0; /* 11/19/05 gam; default is to use double precision, not single. */
  CLA->frameStructType=NULL; /* 01/10/07 gam */

  strcat(allargs, "Command line args: "); /* 06/26/07 gam; copy all command line args into commentField */
  for(i = 0; i < argc; i++)
  {
      strcat(allargs,argv[i]);
      strcat(allargs, " ");
  }
  CLA->commentField=allargs;
      
  /* Scan through list of command line arguments */
  while ( 1 )
  {
    int option_index = 0; /* getopt_long stores long option here */
    int c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch ( c )
    {
    case 'H':
      /* high pass frequency */
      CLA->htdata=1;
      break;
    case 'Z':
      /* make tmp file */
      CLA->makeTmpFile=1;
      break;
    case 'S':
      /* use single precision */
      CLA->useSingle=1;
      break;
    case 'f':
      /* high pass frequency */
      CLA->HPf=atof(optarg);
      break;
    case 't':
      /* SFT time */
      CLA->stringT = optarg;  /* 12/27/05 gam; keep pointer to string that gives the SFT duration */
      CLA->T=atoi(optarg);
      break;
    case 'C':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      break;
    case 's':
      /* GPS start */
      CLA->GPSStart=atof(optarg);
      break;
    case 'e':
      /* GPS end */
      CLA->GPSEnd=atof(optarg);
      break;
    case 'F':
      /* 12/28/05 gam; start frequency */
      FMIN=(REAL8)atof(optarg);
      break;
    case 'B':
      /* 12/28/05 gam; band */
      DF=(REAL8)atof(optarg);
      break;
    case 'D':
      /* 12/27/05 gam; make directories based on GPS time */
      CLA->makeGPSDirs=atof(optarg);
      break;
    case 'v':
      /* 12/28/05 gam; output SFT version; default is version 1 SFTs */
      CLA->sftVersion=atoi(optarg);
      break;
    case 'c':
      /* 12/28/05 gam; comment for version 2 SFTs */
      /* CLA->commentField=optarg; */ /* 06/26/07 gam */
      strcat(CLA->commentField, " Additional comment: "); /* 06/26/07 gam; copy all command line args into commentField */      
      strcat(CLA->commentField,optarg);
      break;
    case 'X':
      /* 12/28/05 gam; misc. part of the SFT description field in the filename (also used if makeGPSDirs > 0) */
      CLA->miscDesc=optarg;
      break;
    case 'u':
      CLA->frameStructType=optarg; /* 01/10/07 gam */
      break;
    case 'a':
      /* Time domain cleaning procedure; the default is 0.*/
      CLA->TDcleaningProc=atoi(optarg);
      break;
    case 'b':
      /*Cut frequency for the bilateral highpass filter. It has to be used only if TDcleaningProc is YES.*/
      CLA->fc=atof(optarg);
      break;
   case 'r':
      /*Critical ratio threshold. It has to be used only if TDcleaningProc is YES.*/
      CLA->cr=atof(optarg);
      break;
    case 'w':
      /* 12/28/05 gam; window options; 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window */
      CLA->windowOption=atoi(optarg);
      break;
    case 'P':
      /* 12/28/05 gam; overlap fraction (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 1.0). */
      CLA->overlapFraction=(REAL8)atof(optarg);
      break;
    case 'N':
      CLA->ChannelName=optarg;       
      break;
    case 'i':
      CLA->IFO=optarg; /* 01/14/07 gam */
      break;
    case 'p':
      CLA->SFTpath=optarg;       
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\thigh-pass-freq (-f)\tFLOAT\t High pass filtering frequency in Hz.\n");
      fprintf(stdout,"\tsft-duration (-t)\tFLOAT\t SFT duration in seconds.\n");
      fprintf(stdout,"\tsft-write-path (-p)\tFLOAT\t Location of output SFTs.\n");
      fprintf(stdout,"\tframe-cache (-C)\tSTRING\t Path to frame cache file (including the filename).\n");
      fprintf(stdout,"\tgps-start-time (-s)\tINT\t GPS start time of segment.\n");
      fprintf(stdout,"\tgps-end-time (-e)\tINT\t GPS end time of segment.\n");
      fprintf(stdout,"\tchannel-name (-N)\tSTRING\t Name of channel to read within a frame.\n");
      fprintf(stdout,"\tifo (-i)\t\tSTRING\t (optional) Name of IFO, i.e., H1, H2, L1, or G1; use if channel name begins with H0, L0, or G0; default: use first two characters from channel name.\n");
      fprintf(stdout,"\tsft-version (-v)\tINT\t (optional) SFT version (1 = default = output version 1 SFTs; 2 = output version 2 SFTs.\n");
      fprintf(stdout,"\tcomment-field (-c)\tSTRING\t (optional) Comment for version 2 SFT header.\n");
      fprintf(stdout,"\tstart-freq (-F) \tFLOAT\t (optional) Start frequency of the SFTs (default is 48 Hz).\n");
      fprintf(stdout,"\tband (-B)       \tFLOAT\t (optional) Frequency band of the SFTs (default is 2000 Hz).\n");
      fprintf(stdout,"\tmake-gps-dirs (-D)\tINT\t (optional) Make directories for output SFTs based on this many digits of the GPS time.\n");
      fprintf(stdout,"\tmake-tmp-file (-Z)\tINT\t (optional) Write SFT to .*.tmp file, then move to final filename.\n");
      fprintf(stdout,"\tmisc-desc (-X)   \tSTRING\t (optional) Misc. part of the SFT description field in the filename (also used if make-gps-dirs, -D option, is > 0)\n");
      fprintf(stdout,"\tTD-cleaning-proc (-a)\tINT\t (optional) Time Domain cleaning procedure; 1=YES and 0=NO; the default is 0.\n");
      fprintf(stdout,"\tcut-frequency (-b)\tFLOAT\t (optional) Cut frequency in Hz for the bilateral highpass filter. It has to be used only if TDcleaningProc=1.\n");
      fprintf(stdout,"\tcrit-ratio-thrs (-r)\tFLOAT\t (optional) Critical ratio threshold. It has to be used only if TDcleaningProc=1.\n");
      fprintf(stdout,"\twindow-type (-w)\tINT\t (optional) 0 = no window, 1 = default = Matlab style Tukey window; 2 = make_sfts.c Tukey window; 3 = Hann window\n");
      fprintf(stdout,"\toverlap-fraction (-P)\tFLOAT\t (optional) Overlap fraction (for use with windows; e.g., use -P 0.5 with -w 3 Hann windows; default is 0.0).\n");
      fprintf(stdout,"\tht-data (-H)\t\tFLAG\t (optional) Input data is h(t) data (input is PROC_REAL8 data ).\n");
      fprintf(stdout,"\tuse-single (-S)\t\tFLAG\t (optional) Use single precision for window, plan, and fft; double precision filtering is always done.\n");
      fprintf(stdout,"\tframe-struct-type (-u)\tSTRING\t (optional) String specifying the input frame structure and data type. Must begin with ADC_ or PROC_ followed by REAL4, REAL8, INT2, INT4, or INT8; default: ADC_REAL4; -H is the same as PROC_REAL8.\n");
      fprintf(stdout,"\thelp (-h)\t\tFLAG\t This message.\n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }
  }

  if(CLA->HPf < 0 )
    {
      fprintf(stderr,"No high pass filtering frequency specified.\n"); 
      fprintf(stderr,"If you don't want to high pass filter set the frequency to 0.\n");
      fprintf(stderr,"Try %s -h \n",argv[0]);
      return 1;
    }      
  if(CLA->T == 0.0)
    {
      fprintf(stderr,"No SFT duration specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      
  if(CLA->GPSStart == 0)
    {
      fprintf(stderr,"No GPS start time specified.\n");
       fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->GPSEnd == 0)
    {
      fprintf(stderr,"No GPS end time specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(FMIN < 0.0)
    {
      fprintf(stderr,"Illegal start-freq option given.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(DF < 0.0)
    {
      fprintf(stderr,"Illegal band option given.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->makeGPSDirs < 0)
    {
      fprintf(stderr,"Illegal make-gps-dirs option given.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }

  if( (CLA->TDcleaningProc < 0) || (CLA->TDcleaningProc > 4) )
    {
      fprintf(stderr,"Illegal specification of TDcleaningProc.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }

 if( (CLA->TDcleaningProc == 0) && ((CLA->fc < 0) || (CLA->fc > 0) || (CLA->cr < 0) || (CLA->cr > 0)) ) 
   {
      fprintf(stderr,"You did NOT specify that you need the cleaning procedure.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
   }

  if( (CLA->sftVersion < 1) || (CLA->sftVersion > 2) )
    {
      fprintf(stderr,"Illegal sft-version given.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if( (CLA->windowOption < 0) || (CLA->windowOption > 3) )
    {
      fprintf(stderr,"Illegal window-type given.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if( (CLA->overlapFraction < 0.0) || (CLA->overlapFraction >= 1.0) )
    {
      fprintf(stderr,"Illegal overlap-fraction given.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->FrCacheFile == NULL)
    {
      fprintf(stderr,"No frame cache file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      
  if(CLA->ChannelName == NULL)
    {
      fprintf(stderr,"No data channel name specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      
  if(CLA->SFTpath == NULL)
    {
      fprintf(stderr,"No output path specified for SFTs.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }      

  return errflg;
}
/*******************************************************************************/

/*******************************************************************************/
int AllocateData(struct CommandLineArgsTag CLA)
{
  static FrChanIn chanin;

  chanin.name  = CLA.ChannelName;

  /* These calls just return deltaT for the channel */
  if(CLA.htdata)
    {
      chanin.type  = ProcDataChannel;
      /* Get channel time step size by calling LALFrGetREAL8TimeSeries */
      LALFrSeek(&status,&gpsepoch,framestream);
      TESTSTATUS( &status );
      LALFrGetREAL8TimeSeries(&status,&dataDouble,&chanin,framestream);
      TESTSTATUS( &status );
      dataSingle.deltaT=dataDouble.deltaT;      
    }
  else if (CLA.frameStructType != NULL) 
    {
      /* 01/10/07 gam */  
      if ( strstr(CLA.frameStructType,"ADC_") ) {
         chanin.type  = ADCDataChannel;
      } else if ( strstr(CLA.frameStructType,"PROC_") ) {
         chanin.type  = ProcDataChannel;
      } else {
        return 1; 
      }
      /* Get channel time step size by calling LALFrGet... functions: */
      LALFrSeek(&status,&gpsepoch,framestream);
      TESTSTATUS( &status );      
      if ( strstr(CLA.frameStructType,"REAL8") ) {
         LALFrGetREAL8TimeSeries(&status,&dataDouble,&chanin,framestream);
         TESTSTATUS( &status );
         dataSingle.deltaT=dataDouble.deltaT;
      } else if ( strstr(CLA.frameStructType,"REAL4") ) { 
         LALFrGetREAL4TimeSeries(&status,&dataSingle,&chanin,framestream);
         TESTSTATUS( &status );
         dataDouble.deltaT=dataSingle.deltaT;
      } else if ( strstr(CLA.frameStructType,"INT2") ) {
         LALFrGetINT2TimeSeries(&status,&dataINT2,&chanin,framestream);
         TESTSTATUS( &status );
         dataDouble.deltaT = dataINT2.deltaT;
         dataSingle.deltaT=dataDouble.deltaT;
      } else if ( strstr(CLA.frameStructType,"INT4") ) {
         LALFrGetINT4TimeSeries(&status,&dataINT4,&chanin,framestream);
         TESTSTATUS( &status );
         dataDouble.deltaT = dataINT4.deltaT;
         dataSingle.deltaT=dataDouble.deltaT;
      } else if ( strstr(CLA.frameStructType,"INT8") ) {      
         LALFrGetINT8TimeSeries(&status,&dataINT8,&chanin,framestream);
         TESTSTATUS( &status );
         dataDouble.deltaT = dataINT8.deltaT;
         dataSingle.deltaT=dataDouble.deltaT;
      } else {
        return 1;
      }      
    }
  else
    {
      chanin.type  = ADCDataChannel;
      /* Get channel time step size by calling LALFrGetREAL4TimeSeries */
      LALFrSeek(&status,&gpsepoch,framestream);
      TESTSTATUS( &status );
      LALFrGetREAL4TimeSeries(&status,&dataSingle,&chanin,framestream);
      TESTSTATUS( &status );
      dataDouble.deltaT=dataSingle.deltaT;
    }
  /*paola*/
  /*  11/19/05 gam; will keep dataDouble or dataSingle in memory at all times, depending on whether using single or double precision */

  if ((CLA.useSingle) || (CLA.TDcleaningProc >= 3)) {  
      LALCreateVector(&status,&dataSingle.data,(UINT4)(CLA.T/dataSingle.deltaT +0.5));
      TESTSTATUS( &status );

      if (CLA.TDcleaningProc > 0) {
  
	LALCreateVector(&status,&dataSingleFirstHP.data,(UINT4)((CLA.T/dataSingle.deltaT +0.5)+1));
	TESTSTATUS( &status );

	LALCreateVector(&status,&dataSingleHP.data,(UINT4)(CLA.T/dataSingle.deltaT +0.5));
	TESTSTATUS( &status );
   
	LALCreateVector(&status,&databil2Single.data,(UINT4)(CLA.T/dataDouble.deltaT +0.5 + 20.0/dataDouble.deltaT));
	TESTSTATUS( &status );

	LALCreateVector(&status,&dataSingleClean.data,(UINT4)(CLA.T/dataSingle.deltaT +0.5));
	TESTSTATUS( &status );
      }     


      #if TRACKMEMUSE
         printf("Memory use after creating dataSingle and before calling LALCreateForwardRealFFTPlan.\n"); printmemuse();
      #endif

      LALCreateForwardRealFFTPlan( &status, &fftPlanSingle, dataSingle.data->length, 0 );
      TESTSTATUS( &status );

      #if TRACKMEMUSE
         printf("Memory use after creating dataSingle and after calling LALCreateForwardRealFFTPlan.\n"); printmemuse();
      #endif

  }

  if (!CLA.useSingle) {  
   
      LALDCreateVector(&status,&dataDouble.data,(UINT4)(CLA.T/dataDouble.deltaT +0.5));
      TESTSTATUS( &status );

      if ((CLA.TDcleaningProc > 0) && (CLA.TDcleaningProc < 3)) {  
  
	LALDCreateVector(&status,&dataDoubleFirstHP.data,(UINT4)((CLA.T/dataDouble.deltaT +0.5)+1));
	TESTSTATUS( &status );

	LALDCreateVector(&status,&dataDoubleHP.data,(UINT4)(CLA.T/dataDouble.deltaT +0.5));
	TESTSTATUS( &status );
   
 
	LALDCreateVector(&status,&databil2.data,(UINT4)(CLA.T/dataDouble.deltaT +0.5 + 20.0/dataDouble.deltaT));
	TESTSTATUS( &status );


	LALDCreateVector(&status,&dataDoubleClean.data,(UINT4)(CLA.T/dataDouble.deltaT +0.5));
	TESTSTATUS( &status );
      }

      #if TRACKMEMUSE
         printf("Memory use after creating dataDouble and before calling LALCreateForwardREAL8FFTPlan.\n"); printmemuse();
      #endif
      
      LALCreateForwardREAL8FFTPlan( &status, &fftPlanDouble, dataDouble.data->length, 0 );
      TESTSTATUS( &status );

      #if TRACKMEMUSE
            printf("Memory use after creating dataDouble and after calling LALCreateForwardREAL8FFTPlan.\n"); printmemuse();
      #endif
    
  }

  return 0;
}
/*******************************************************************************/

/*******************************************************************************/
int ReadData(struct CommandLineArgsTag CLA)
{
  static FrChanIn chanin;
  INT4 k;
  INT4 i;/*mydatatest*/
  INT4 itest=0;
  chanin.name  = CLA.ChannelName;
  LALFrSeek(&status,&gpsepoch,framestream);
  TESTSTATUS( &status );

  FILE *outputtest;
 
  /* 11/19/05 gam */
  if(CLA.useSingle) {
    
    if(CLA.htdata)
    {

      #if TRACKMEMUSE
        printf("Memory use before creating dataDouble and calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
      #endif

      LALDCreateVector(&status,&dataDouble.data,(UINT4)(CLA.T/dataDouble.deltaT +0.5));
      TESTSTATUS( &status );

      chanin.type  = ProcDataChannel;
      LALFrGetREAL8TimeSeries(&status,&dataDouble,&chanin,framestream);
      TESTSTATUS( &status );

      #if TRACKMEMUSE
        printf("Memory use after creating dataDouble and calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
      #endif

      /*copy the data into the single precision timeseries */
      for (k = 0; k < (int)dataSingle.data->length; k++) {
          dataSingle.data->data[k] = dataDouble.data->data[k];
      }

      LALDDestroyVector(&status,&dataDouble.data);
      TESTSTATUS( &status );

      #if TRACKMEMUSE
        printf("Memory use after destroying dataDouble.\n"); printmemuse();
      #endif

    }
    else if (CLA.frameStructType != NULL)
    {
      /* 01/10/07 gam */  
      if ( strstr(CLA.frameStructType,"ADC_") ) {
         chanin.type  = ADCDataChannel;
      } else if ( strstr(CLA.frameStructType,"PROC_") ) {
         chanin.type  = ProcDataChannel;
      } else {
        return 1; 
      }
      if ( strstr(CLA.frameStructType,"REAL8") ) {
         #if TRACKMEMUSE
           printf("Memory use before creating dataDouble and calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
         #endif

         LALDCreateVector(&status,&dataDouble.data,(UINT4)(CLA.T/dataDouble.deltaT +0.5));
         TESTSTATUS( &status );

         LALFrGetREAL8TimeSeries(&status,&dataDouble,&chanin,framestream);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after creating dataDouble and calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
         #endif

         /*copy the data into the single precision timeseries */
         for (k = 0; k < (int)dataSingle.data->length; k++) {
             dataSingle.data->data[k] = dataDouble.data->data[k];
         }

         LALDDestroyVector(&status,&dataDouble.data);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after destroying dataDouble.\n"); printmemuse();
         #endif
      } else if ( strstr(CLA.frameStructType,"REAL4") ) {      
         #if TRACKMEMUSE
           printf("Memory use before calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
         #endif

         LALFrGetREAL4TimeSeries(&status,&dataSingle,&chanin,framestream);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
         #endif
      } else if ( strstr(CLA.frameStructType,"INT2") ) {
         #if TRACKMEMUSE
           printf("Memory use before creating dataINT2 and calling LALFrGetINT2TimeSeries.\n"); printmemuse();
         #endif

         LALI2CreateVector(&status,&dataINT2.data,(UINT4)(CLA.T/dataINT2.deltaT +0.5));
         TESTSTATUS( &status );

         LALFrGetINT2TimeSeries(&status,&dataINT2,&chanin,framestream);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after creating dataINT2 and calling LALFrGetINT2TimeSeries.\n"); printmemuse();
         #endif

         /*copy the data into the single precision timeseries */
         for (k = 0; k < (int)dataDouble.data->length; k++) {
             dataSingle.data->data[k] = (REAL4)(dataINT2.data->data[k]);
         }     

         LALI2DestroyVector(&status,&dataINT2.data);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after destroying dataINT2.\n"); printmemuse();
         #endif
      } else if ( strstr(CLA.frameStructType,"INT4") ) {
         #if TRACKMEMUSE
           printf("Memory use before creating dataINT4 and calling LALFrGetINT4TimeSeries.\n"); printmemuse();
         #endif

         LALI4CreateVector(&status,&dataINT4.data,(UINT4)(CLA.T/dataINT4.deltaT +0.5));
         TESTSTATUS( &status );

         LALFrGetINT4TimeSeries(&status,&dataINT4,&chanin,framestream);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after creating dataINT4 and calling LALFrGetINT4TimeSeries.\n"); printmemuse();
         #endif

         /*copy the data into the single precision timeseries */
         for (k = 0; k < (int)dataDouble.data->length; k++) {
             dataSingle.data->data[k] = (REAL4)(dataINT4.data->data[k]);
         }     

         LALI4DestroyVector(&status,&dataINT4.data);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after destroying dataINT4.\n"); printmemuse();
         #endif
      } else if ( strstr(CLA.frameStructType,"INT8") ) {
         #if TRACKMEMUSE
           printf("Memory use before creating dataINT8 and calling LALFrGetINT8TimeSeries.\n"); printmemuse();
         #endif

         LALI8CreateVector(&status,&dataINT8.data,(UINT4)(CLA.T/dataINT8.deltaT +0.5));
         TESTSTATUS( &status );

         LALFrGetINT8TimeSeries(&status,&dataINT8,&chanin,framestream);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after creating dataINT8 and calling LALFrGetINT8TimeSeries.\n"); printmemuse();
         #endif

         /*copy the data into the single precision timeseries */
         for (k = 0; k < (int)dataDouble.data->length; k++) {
             dataSingle.data->data[k] = (REAL4)(dataINT8.data->data[k]);
         }     

         LALI8DestroyVector(&status,&dataINT8.data);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after destroying dataINT8.\n"); printmemuse();
         #endif
      } else {
        return 1;
      }      
    }
    else
    {

      #if TRACKMEMUSE
        printf("Memory use before calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
      #endif

      chanin.type  = ADCDataChannel;
      LALFrGetREAL4TimeSeries(&status,&dataSingle,&chanin,framestream);
      TESTSTATUS( &status );

      #if TRACKMEMUSE
        printf("Memory use after calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
      #endif

    }

    #if PRINTEXAMPLEDATA
        printf("\nExample dataSingle values after reading data from frames in ReadData:\n"); printExampleDataSingle();
    #endif

  } else {

    if(CLA.htdata)
    {
    
      #if TRACKMEMUSE
        printf("Memory use before calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
      #endif
    
      chanin.type  = ProcDataChannel;
      LALFrGetREAL8TimeSeries(&status,&dataDouble,&chanin,framestream);
      TESTSTATUS( &status );


    if(itest==1){
             outputtest=fopen("datadoubleOR.dat","w");
             for(i=0; i< (int)dataDouble.data->length; i++)fprintf(outputtest,"%23.16e\n",dataDouble.data->data[i]);
             fclose(outputtest);
       }
/*mydatatest!HERE!*/


      #if TRACKMEMUSE
        printf("Memory use after calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
      #endif

    }
    else if (CLA.frameStructType != NULL) 
    {
      /* 01/10/07 gam */  
      if ( strstr(CLA.frameStructType,"ADC_") ) {
         chanin.type  = ADCDataChannel;
      } else if ( strstr(CLA.frameStructType,"PROC_") ) {
         chanin.type  = ProcDataChannel;
      } else {
        return 1; 
      }
      if ( strstr(CLA.frameStructType,"REAL8") ) {
         #if TRACKMEMUSE
           printf("Memory use before calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
         #endif

         LALFrGetREAL8TimeSeries(&status,&dataDouble,&chanin,framestream);
         TESTSTATUS( &status );


         #if TRACKMEMUSE
           printf("Memory use after calling LALFrGetREAL8TimeSeries.\n"); printmemuse();
         #endif
      } else if ( strstr(CLA.frameStructType,"REAL4") ) {      
         #if TRACKMEMUSE
           printf("Memory use before creating dataSingle and calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
         #endif

         LALCreateVector(&status,&dataSingle.data,(UINT4)(CLA.T/dataSingle.deltaT +0.5));
         TESTSTATUS( &status );

         LALFrGetREAL4TimeSeries(&status,&dataSingle,&chanin,framestream);
         TESTSTATUS( &status );
 
         #if TRACKMEMUSE
           printf("Memory use after creating dataSingle and calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
         #endif

         /*copy the data into the double precision timeseries */
         for (k = 0; k < (int)dataDouble.data->length; k++) {
             dataDouble.data->data[k] = dataSingle.data->data[k];
         }     

         LALDestroyVector(&status,&dataSingle.data);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after destroying dataSingle.\n"); printmemuse();
         #endif      
      } else if ( strstr(CLA.frameStructType,"INT2") ) {
         #if TRACKMEMUSE
           printf("Memory use before creating dataINT2 and calling LALFrGetINT2TimeSeries.\n"); printmemuse();
         #endif

         LALI2CreateVector(&status,&dataINT2.data,(UINT4)(CLA.T/dataINT2.deltaT +0.5));
         TESTSTATUS( &status );

         LALFrGetINT2TimeSeries(&status,&dataINT2,&chanin,framestream);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after creating dataINT2 and calling LALFrGetINT2TimeSeries.\n"); printmemuse();
         #endif

         /*copy the data into the double precision timeseries */
         for (k = 0; k < (int)dataDouble.data->length; k++) {
             dataDouble.data->data[k] = (REAL8)(dataINT2.data->data[k]);
         }     

         LALI2DestroyVector(&status,&dataINT2.data);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after destroying dataINT2.\n"); printmemuse();
         #endif
      } else if ( strstr(CLA.frameStructType,"INT4") ) {
         #if TRACKMEMUSE
           printf("Memory use before creating dataINT4 and calling LALFrGetINT4TimeSeries.\n"); printmemuse();
         #endif

         LALI4CreateVector(&status,&dataINT4.data,(UINT4)(CLA.T/dataINT4.deltaT +0.5));
         TESTSTATUS( &status );

         LALFrGetINT4TimeSeries(&status,&dataINT4,&chanin,framestream);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after creating dataINT4 and calling LALFrGetINT4TimeSeries.\n"); printmemuse();
         #endif

         /*copy the data into the double precision timeseries */
         for (k = 0; k < (int)dataDouble.data->length; k++) {
             dataDouble.data->data[k] = (REAL8)(dataINT4.data->data[k]);
         }     

         LALI4DestroyVector(&status,&dataINT4.data);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after destroying dataINT4.\n"); printmemuse();
         #endif
      } else if ( strstr(CLA.frameStructType,"INT8") ) {
         #if TRACKMEMUSE
           printf("Memory use before creating dataINT8 and calling LALFrGetINT8TimeSeries.\n"); printmemuse();
         #endif

         LALI8CreateVector(&status,&dataINT8.data,(UINT4)(CLA.T/dataINT8.deltaT +0.5));
         TESTSTATUS( &status );

         LALFrGetINT8TimeSeries(&status,&dataINT8,&chanin,framestream);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after creating dataINT8 and calling LALFrGetINT8TimeSeries.\n"); printmemuse();
         #endif

         /*copy the data into the double precision timeseries */
         for (k = 0; k < (int)dataDouble.data->length; k++) {
             dataDouble.data->data[k] = (REAL8)(dataINT8.data->data[k]);
         }     

         LALI8DestroyVector(&status,&dataINT8.data);
         TESTSTATUS( &status );

         #if TRACKMEMUSE
           printf("Memory use after destroying dataINT8.\n"); printmemuse();
         #endif
      } else {
        return 1;
      }      
    }
    else
    {
      #if TRACKMEMUSE
        printf("Memory use before creating dataSingle and calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
      #endif

      /* 11/02/05 gam; add next two lines */
      LALCreateVector(&status,&dataSingle.data,(UINT4)(CLA.T/dataSingle.deltaT +0.5));
      TESTSTATUS( &status );
      
      chanin.type  = ADCDataChannel;
      LALFrGetREAL4TimeSeries(&status,&dataSingle,&chanin,framestream);
      TESTSTATUS( &status );
 
      #if TRACKMEMUSE
        printf("Memory use after creating dataSingle and calling LALFrGetREAL4TimeSeries.\n"); printmemuse();
      #endif

      /*copy the data into the double precision timeseries */
      for (k = 0; k < (int)dataDouble.data->length; k++) {
          dataDouble.data->data[k] = dataSingle.data->data[k];
      }     

      LALDestroyVector(&status,&dataSingle.data);
      TESTSTATUS( &status );

      #if TRACKMEMUSE
        printf("Memory use after destroying dataSingle.\n"); printmemuse();
      #endif

    }

    #if PRINTEXAMPLEDATA
        printf("\nExample dataDouble values after reading data from frames in ReadData:\n"); printExampleDataDouble();
    #endif

  } /* END if(CLA.useSingle) else */

  return 0;
}
/*******************************************************************************/

/*******************************************************************************/
int HighPass(struct CommandLineArgsTag CLA)
{
  INT4 i;
  FILE *out;

  PassBandParamStruc filterpar;

  filterpar.name  = "Butterworth High Pass";
  filterpar.nMax  = 10;
  filterpar.f2    = CLA.HPf;
  filterpar.a2    = 0.5;
  filterpar.f1    = -1.0;
  filterpar.a1    = -1.0;

  if (CLA.HPf > 0.0 )
    {
      if(CLA.useSingle) {
        #if TRACKMEMUSE
          printf("Memory use before calling LALDButterworthREAL4TimeSeries:\n"); printmemuse();
        #endif

        /* High pass the time series */
        LALDButterworthREAL4TimeSeries(&status,&dataSingle,&filterpar);
        TESTSTATUS( &status );

        #if TRACKMEMUSE
          printf("Memory use after calling LALDButterworthREAL4TimeSeries:\n"); printmemuse();
        #endif

        #if PRINTEXAMPLEDATA
            printf("\nExample dataSingle values after filtering data in HighPass:\n"); printExampleDataSingle();
        #endif
      } else {
        #if TRACKMEMUSE
          printf("Memory use before calling LALButterworthREAL8TimeSeries:\n"); printmemuse();
        #endif

        /* High pass the time series */
        LALButterworthREAL8TimeSeries(&status,&dataDouble,&filterpar);
        TESTSTATUS( &status );
   
        
	if(testpao){
	  out=fopen("dataDligoHP.dat","w");
	  for(i=0; i< (int)dataDouble.data->length; i++)fprintf(out,"%23.16e\n",dataDouble.data->data[i]);
	  fclose(out);
	}



        #if TRACKMEMUSE
          printf("Memory use after calling LALButterworthREAL8TimeSeries:\n"); printmemuse();
        #endif

        #if PRINTEXAMPLEDATA
            printf("\nExample dataDouble values after filtering data in HighPass:\n"); printExampleDataDouble();
        #endif
      }
    }

  return 0;
}


/*******************************************************************************/

/*******************************************************************************/
int WindowData(struct CommandLineArgsTag CLA)
{
  REAL8 r = 0.001;
  INT4 k, N, kl, kh;
  FILE *OUTpao;
  
  /* This implementation of a Tukey window is describes 
     in the Matlab tukey window documentation */

  /* 11/19/05 gam */
  if(CLA.useSingle) {
    N=dataSingle.data->length;
    kl=r/2*(N-1)+1;
    kh=N-r/2*(N-1)+1;
    for(k = 1; k < kl; k++) 
    {
      dataSingle.data->data[k-1]=dataSingle.data->data[k-1]*0.5*( (REAL4)( 1.0 + cos(LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) ) );
    }
    for(k = kh; k <= N; k++) 
    {
      dataSingle.data->data[k-1]=dataSingle.data->data[k-1]*0.5*( (REAL4)( 1.0 + cos(LAL_TWOPI/r - LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) ) );
    }

    #if PRINTEXAMPLEDATA
        printf("\nExample dataSingle values after windowing data in WindowData:\n"); printExampleDataSingle();
    #endif
  } else {
    N=dataDouble.data->length;
    kl=r/2*(N-1)+1;
    kh=N-r/2*(N-1)+1;
    for(k = 1; k < kl; k++) 
    {
      dataDouble.data->data[k-1]=dataDouble.data->data[k-1]*0.5*( 1.0 + cos(LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) );
    }
    for(k = kh; k <= N; k++) 
    {
      dataDouble.data->data[k-1]=dataDouble.data->data[k-1]*0.5*( 1.0 + cos(LAL_TWOPI/r - LAL_TWOPI/r*(k-1)/(N-1) - LAL_PI) );
    }
    #if PRINTEXAMPLEDATA
        printf("\nExample dataDouble values after windowing data in WindowData:\n"); printExampleDataDouble();
    #endif  

   if(testpao){
     OUTpao=fopen("windowT.dat","w");
     for(k =0; k < N; k++)fprintf(OUTpao,"%23.16e\n",dataDouble.data->data[k]);
     fclose(OUTpao);
   }
  
  }

  return 0;
}
/*******************************************************************************/

/*******************************************************************************/
/* Same as window function given in lalapps/src/pulsar/make_sfts.c */
int WindowDataTukey2(struct CommandLineArgsTag CLA)
{  
  /* Define the parameters to make the window */
  INT4 WINSTART = 4096;
  INT4 WINEND = 8192;
  INT4 WINLEN = (WINEND-WINSTART);
  INT4 i;
  REAL8 win;
    
  if(CLA.useSingle) {
        /* window data.  Off portion */
        for (i=0; i<WINSTART; i++) {
          dataSingle.data->data[i] = 0.0;
          dataSingle.data->data[dataSingle.data->length - 1 - i] = 0.0;
        }
        /* window data, smooth turn-on portion */
        for (i=WINSTART; i<WINEND; i++) {
          win=((sin((i - WINSTART)*LAL_PI/(WINLEN)-LAL_PI_2)+1.0)/2.0);
          dataSingle.data->data[i] *= win;
          dataSingle.data->data[dataSingle.data->length - 1 - i]  *= win;
        }
        #if PRINTEXAMPLEDATA
           printf("\nExample dataSingle values after windowing data in WindowDataTukey2:\n"); printExampleDataSingle();
        #endif
  } else {
        /* window data.  Off portion */
        for (i=0; i<WINSTART; i++) {
          dataDouble.data->data[i] = 0.0;
          dataDouble.data->data[dataDouble.data->length - 1 - i] = 0.0;
        }
        /* window data, smooth turn-on portion */
        for (i=WINSTART; i<WINEND; i++) {
          win=((sin((i - WINSTART)*LAL_PI/(WINLEN)-LAL_PI_2)+1.0)/2.0);
          dataDouble.data->data[i] *= win;
          dataDouble.data->data[dataDouble.data->length - 1 - i]  *= win;
        }
        #if PRINTEXAMPLEDATA
            printf("\nExample dataDouble values after windowing data in WindowDataTukey2:\n"); printExampleDataDouble();
        #endif  
  }

  return 0;
}
/*******************************************************************************/

/*******************************************************************************/
/* Hann window based on Matlab, but with C indexing: w[k] = 0.5*( 1 - cos(2*pi*k/(N-1)) ) k = 0, 1, 2,...N-1 */
int WindowDataHann(struct CommandLineArgsTag CLA)
{
  INT4 k;
  REAL8 win,N,Nm1;
  REAL8 real8TwoPi = 2.0*((REAL8)(LAL_PI));

  if(CLA.useSingle) {
        N = ((REAL8)dataSingle.data->length);
        Nm1 = N - 1;
        for (k=0; k<N; k++) {
          win=0.5*( 1.0 - cos(real8TwoPi*((REAL8)(k))/Nm1) );
          dataSingle.data->data[k] *= win;
        }
        #if PRINTEXAMPLEDATA
           printf("\nExample dataSingle values after windowing data in WindowDataHann:\n"); printExampleDataSingle();
        #endif
  } else {
        N = ((REAL8)dataDouble.data->length);
        Nm1 = N - 1;
        for (k=0; k<N; k++) {
          win=0.5*( 1.0 - cos(real8TwoPi*((REAL8)(k))/Nm1) );
          dataDouble.data->data[k] *= win;
        }  
        #if PRINTEXAMPLEDATA
            printf("\nExample dataDouble values after windowing data in WindowDataHann:\n"); printExampleDataDouble();
        #endif  
  }

  return 0;
}


/*******************************************************************************/
/*----------FUNCTIONS for the TIME DOMAIN CLEANING PROCEDURE---------------*/
/*******************************************************************************/
/*********************************************/
/*Initialization of the event paramaters*/
/*********************************************/
int EventParamInit(long len, ParamOfEvent *even_param)
{ 
  int i;
  long lent; 
  lent=len;
  /*EVEN_PARAM *even_param; Parameters to check for events and clean the data*/
  /*even_param=(EVEN_PARAM *)calloc(1,sizeof(EVEN_PARAM)); */
  even_param->tau=20.0;       /*memory time of the autoregressive average; */   
  if(lent <= 131072) even_param->tau=10.0; 
  even_param->deadtime=0.1;  /*Dead time in seconds*/
  even_param->factor=1.0;    /*e.g. 10: when re-evaluate mean and std*/
  even_param->iflev=0;       /*0=no events; 1=begin; 2=event; 3=end of the event*/
  even_param->med=(REAL4*)malloc(lent*sizeof(REAL4));       /*Autoregressive mean*/
  even_param->std=(REAL4*)malloc(lent*sizeof(REAL4));   
  /*even_param->xamed=(REAL4*)malloc(lent*sizeof(REAL4));       Autoregressive mean*/
 /* even_param->xastd=(REAL4*)malloc(lent*sizeof(REAL4));      Autoregressive std*/
  even_param->xw=0.;
  even_param->qw=0.;
  even_param->w_norm=1.;
  even_param->edge=0.00061035;  /*Right value: 0.00061035; 0.15 s Seconds around (before and after) the event have to be excluded"*/
  even_param->total=0;    /*Samples, due to all events, have been cleaned*/
  even_param->number=0;   /*Number of events*/
  even_param->ar=1;       /*Algorithm for the autoregressive evaluation; it has to be = 1*/
  even_param->begin=malloc((size_t) len*sizeof(int)); /*The dimension is the number of events, less than len, but here len is used, just for simplicity*/
  even_param->duration=malloc((size_t) len*sizeof(long)); 
  even_param->imax=malloc((size_t) len*sizeof(long));
  even_param->crmax=malloc((size_t) len*sizeof(REAL4));
  even_param->ener=malloc((size_t) len*sizeof(REAL4));
  for(i=0;i<lent;i++){
   /* even_param->xamed[i]=0.;
    even_param->xastd[i]=0.;*/
    even_param->med[i]=0.;
    even_param->std[i]=0.;
    even_param->begin[i]=-1;
    even_param->duration[i]=0;
    even_param->imax[i]=0;
    even_param->crmax[i]=0;
    even_param->ener[i]=0;
  }
  return 0;
}

int EventParamInit_NoBil(long len, ParamOfEvent *even_param)
{ 
  int i;
  long lent; 
  lent=len;
  /*EVEN_PARAM *even_param; Parameters to check for events and clean the data*/
  /*even_param=(EVEN_PARAM *)calloc(1,sizeof(EVEN_PARAM)); */
  even_param->tau=20.0;       /*memory time of the autoregressive average; */   
  if(lent <= 131072) even_param->tau=10.0; 
  even_param->deadtime=0.1;  /*Dead time in seconds*/
  even_param->factor=1.0;    /*e.g. 10: when re-evaluate mean and std*/
  even_param->iflev=0;       /*0=no events; 1=begin; 2=event; 3=end of the event*/
  even_param->med=(REAL4*)malloc(lent*sizeof(REAL4));       /*Autoregressive mean*/
  even_param->std=(REAL4*)malloc(lent*sizeof(REAL4));   
  even_param->xamed=(REAL4*)malloc(lent*sizeof(REAL4));      /*Autoregressive mean*/
  even_param->xastd=(REAL4*)malloc(lent*sizeof(REAL4));    /*Autoregressive std*/
  even_param->xw=0.;
  even_param->qw=0.;
  even_param->w_norm=1.;
  even_param->edge=0.00061035;  /*Right value: 0.00061035; 0.15 s Seconds around (before and after) the event have to be excluded"*/
  even_param->total=0;    /*Samples, due to all events, have been cleaned*/
  even_param->number=0;   /*Number of events*/
  even_param->ar=1;       /*Algorithm for the autoregressive evaluation; it has to be = 1*/
  even_param->begin=malloc((size_t) len*sizeof(int)); /*The dimension is the number of events, less than len, but here len is used, just for simplicity*/
  even_param->duration=malloc((size_t) len*sizeof(long)); 
  even_param->imax=malloc((size_t) len*sizeof(long));
  even_param->crmax=malloc((size_t) len*sizeof(REAL4));
  even_param->ener=malloc((size_t) len*sizeof(REAL4));
  for(i=0;i<lent;i++){
    even_param->xamed[i]=0.;
    even_param->xastd[i]=0.;
    even_param->med[i]=0.;
    even_param->std[i]=0.;
    even_param->begin[i]=-1;
    even_param->duration[i]=0;
    even_param->imax[i]=0;
    even_param->crmax[i]=0;
    even_param->ener[i]=0;
  }
  return 0;
}


/*********************************************/
/*Initialization of the event paramaters for the 'fake' highpass series*/
/*********************************************/
int Evenbil2(ParamOfEvent *even_param, REAL8TimeSeries *series )
{ 
  int i;
  long newlen; 
  newlen=ceil(even_param->tau/series->deltaT)+series->data->length;
  even_param->xamed=(REAL4*)malloc(newlen*sizeof(REAL4));      
  even_param->xastd=(REAL4*)malloc(newlen*sizeof(REAL4));    
  for(i=0;i<(int)(newlen);i++){
    even_param->xamed[i]=0.;
    even_param->xastd[i]=0.;
  }
  
  return 0;
}

int Evenbil2_R4(ParamOfEvent *even_param, REAL4TimeSeries *series )
{ 
  int i;
  long newlen; 
  newlen=ceil(even_param->tau/series->deltaT)+series->data->length;
  even_param->xamed=(REAL4*)malloc(newlen*sizeof(REAL4));      
  even_param->xastd=(REAL4*)malloc(newlen*sizeof(REAL4));    
  for(i=0;i<(int)(newlen);i++){
    even_param->xamed[i]=0.;
    even_param->xastd[i]=0.;
  }
  
  return 0;
}

/******************************************************************************************************************************************/
/*BILATERAL HIGHPASS FILTER: the SINGLE-PRECISION data are highpassed in the time domain. It's useful to look for high frequency spikes!*/
/******************************************************************************************************************************************/
int BilHighpass_dataSingle( REAL4TimeSeries *seriesHP, REAL4TimeSeries *firsthighpass, REAL4TimeSeries *series, ParamOfEvent *even_param, BilHP *myparams )
{

        INT4 i,k,len;
	REAL4 w,w1,b0,b1,a1;
        int itest;
	       
        FILE *OUTtest;
        itest=0;
 
        seriesHP->data->length=series->data->length;   
	seriesHP->deltaT=series->deltaT; 
	w=exp(-2*LAL_PI*(myparams->fcut)*(series->deltaT));
	w1=(1+w)/2;
	b0=w1;
	b1=-w1;
	a1=-w;
        firsthighpass->data->data[0]=b0*series->data->data[0]; 
	for(i=1;i<(INT4)series->data->length;i++){
	   firsthighpass->data->data[i]=b0*series->data->data[i]+b1*series->data->data[i-1]-a1*firsthighpass->data->data[i-1];        
	}
        len=(series->data->length)-1;
	seriesHP->data->data[len]=b0*firsthighpass->data->data[len]; 
        for(i=len-1;i>=0;i--){
	  seriesHP->data->data[i]=b0*firsthighpass->data->data[i]+b1*firsthighpass->data->data[i+1]-a1*seriesHP->data->data[i+1];
	}
      
        if(itest==1){
         OUTtest=fopen("testhp_singleprec.dat","w");
         for(k=0;k<(INT4)series->data->length;k++)
	   fprintf(OUTtest,"%i %23.16e %23.16e\n",k,seriesHP->data->data[k],series->data->data[k]);
         fclose(OUTtest);

       }

	return i;

}

/********************************************************************************************************************************************/
/*BILATERAL HIGHPASS FILTER: the DOUBLE-PRECISION data are highpassed in the time domain. It's useful to look for high frequency spikes!*/
/********************************************************************************************************************************************/
/*int BilHighpass_dataDouble( REAL8TimeSeries *seriesHP, REAL8TimeSeries *firsthighpass, REAL8TimeSeries *series, BilHP *myparams )*/
int BilHighpass_dataDouble( REAL8TimeSeries *seriesHP, REAL8TimeSeries *firsthighpass, REAL8TimeSeries *series, ParamOfEvent *even_param, BilHP *myparams )
{

        INT4 i,k,len;
	REAL4 w,w1,b0,b1,a1;
        int itest;
        
        FILE *OUTtest;
        itest=0;
	 
	/* FIXME:
	  with timeseries allocated outside the function this should be checked
	  and an error should be returned if they are unequal, but a length should
	  never be modified !
	*/
        seriesHP->data->length=series->data->length;
	/* seriesHP->data->length=series->data->length;*/
	seriesHP->deltaT=series->deltaT; 
	if(0)
	  w=exp(-2*LAL_PI*(myparams->fcut)*(series->deltaT));
	else
	  w=exp(-2*PSS_PI*(myparams->fcut)*(series->deltaT));
	fprintf(stderr,"BilHighpass_dataDouble(): f:%23.16e, n:%23.16e, w:%23.16e\n",myparams->fcut,series->deltaT,w);
	w1=(1+w)/2;
	b0=w1;
	b1=-w1;
	a1=-w;

        printf("SerieLength = %i\n",series->data->length);
        firsthighpass->data->data[0]=b0*series->data->data[0]; 
	for(i=1;i<(INT4)series->data->length;i++){
	   firsthighpass->data->data[i]=b0*series->data->data[i]+b1*series->data->data[i-1]-a1*firsthighpass->data->data[i-1];        
	}
        len=(series->data->length)-1;
	seriesHP->data->data[len]=b0*firsthighpass->data->data[len]; 
        for(i=len-1;i>=0;i--){
	  seriesHP->data->data[i]=b0*firsthighpass->data->data[i]+b1*firsthighpass->data->data[i+1]-a1*seriesHP->data->data[i+1];
	}

            
	/* firsthighpass->data->data[ts]=b0*series->data->data[ts]; 
	for(i=ts;i<(series->data->length-ts);i++){
	   firsthighpass->data->data[i]=b0*series->data->data[i]+b1*series->data->data[i-1]-a1*firsthighpass->data->data[i-1];        
	}
        len=(series->data->length-ts)-1;
	seriesHP->data->data[len]=b0*firsthighpass->data->data[len]; 
        for(i=len-1;i>=ts;i--){
	  seriesHP->data->data[i]=b0*firsthighpass->data->data[i]+b1*firsthighpass->data->data[i+1]-a1*seriesHP->data->data[i+1];
	}*/

      
       if(itest==1){
         OUTtest=fopen("testhp.dat","w");
         for(k=0;k<(INT4)series->data->length;k++)fprintf(OUTtest,"%23.16e %23.16e\n",seriesHP->data->data[k],series->data->data[k]);
         fclose(OUTtest);

       }
       return i;

}

int BilHighpass_dataDouble_R4( REAL4TimeSeries *seriesHP, REAL4TimeSeries *firsthighpass, REAL4TimeSeries *series, BilHP *myparams )
{
        INT4 i,last;
	REAL4 w,w1,b0,b1,a1;
        
        if(seriesHP->data->length != series->data->length) {
	  fprintf(stderr,"ERROR: BilHighpass_dataDouble_R4(): length of timeseries differs: %u != %u\n",
		  seriesHP->data->length, series->data->length);
	  return -1;
	}
        if(firsthighpass->data->length <= series->data->length) {
	  fprintf(stderr,"ERROR: BilHighpass_dataDouble_R4(): length of firsthighpass timeseries too short: %u <= %u\n",
		  firsthighpass->data->length, series->data->length);
	  return -1;
	}

	seriesHP->deltaT = series->deltaT; 

	w  = exp(-2.0f * 3.1415926f * (myparams->fcut) * (series->deltaT));
	w1 = (1.0f + w) / 2.0f;
	b0 = w1;
	b1 = -w1;
	a1 = -w;

	fprintf(stderr,"BilHighpass_dataDouble(): f:%23.16e, n:%23.16e, w:%23.16e\n",myparams->fcut,series->deltaT,w);

        firsthighpass->data->data[0] = b0 * series->data->data[0]; 
	for(i = 1; i < (INT4)series->data->length; i++) {
	  firsthighpass->data->data[i]
	    = b0 * series->data->data[i]
	    + b1 * series->data->data[i-1]
	    - a1 * firsthighpass->data->data[i-1];        
	}

        last = series->data->length - 1;
	seriesHP->data->data[last] = b0 * firsthighpass->data->data[last]; 
        for(i = last - 1; i >= 0; i--) {
	  seriesHP->data->data[i]
	    = b0 * firsthighpass->data->data[i]
	    + b1 * firsthighpass->data->data[i+1]
	    - a1 * seriesHP->data->data[i+1];
	}

	return 0;
}


/********************************************************************************************************************************************/
/*BILATERAL PAOLA HIGHPASS FILTER: the DOUBLE-PRECISION data are highpassed in the time domain. It's useful to look for high frequency spikes!*/
/********************************************************************************************************************************************/
int Bil2( REAL8TimeSeries *seriesHPbil2, REAL8TimeSeries *seriesHP, ParamOfEvent *even_param )
{

  INT4 i,k;
  int itest;
            
  INT8 ts;
  ts=ceil(even_param->tau/seriesHP->deltaT);

  FILE *OUTtest;
  itest=0;

  seriesHPbil2->data->length=(seriesHP->data->length+ts);
  seriesHPbil2->deltaT=seriesHP->deltaT; 

  for(i=0;i<(int)ts;i++){
    seriesHPbil2->data->data[i]=seriesHP->data->data[i];
  } 
  for(i=ts;i<(int)(seriesHP->data->length+ts);i++){
    /* for(j=0;i<(int)(seriesHP->data->length);j++)seriesHPbil2->data->data[j]=seriesHP->data->data[j];*/
    seriesHPbil2->data->data[i]=seriesHP->data->data[i-ts];
  } 
          
  printf("SerieLengthBil2 = %i\n",seriesHPbil2->data->length);

  if(itest==1){
    OUTtest=fopen("testhp.dat","w");
    for(k=0;k<(seriesHP->data->length+ts);k++)fprintf(OUTtest,"%23.16e \n",seriesHPbil2->data->data[k]);
    fclose(OUTtest);
  }

  return i;
}


int Bil2_R4( REAL4TimeSeries *seriesHPbil2, REAL4TimeSeries *seriesHP, ParamOfEvent *even_param )
{
  UINT4 i, ts;

  ts = ceil(even_param->tau / seriesHP->deltaT);

  seriesHPbil2->data->length=(seriesHP->data->length+ts);
  seriesHPbil2->deltaT=seriesHP->deltaT; 

  for(i = 0; i < ts; i++){
    seriesHPbil2->data->data[i] = seriesHP->data->data[i];
  } 
  for(i = ts; i < seriesHP->data->length + ts; i++){
    seriesHPbil2->data->data[i] = seriesHP->data->data[i-ts];
  } 

  return i;
}


/**********************************************************************************************************************************************/
/*COMPUTATION OF MEAN AND STANDARD DEVIATION for data in SINGLE-PRECISION: the calculations are made by means of an autoregressive procedure*/
/**********************************************************************************************************************************************/
int MedStd_dataSingle( REAL4TimeSeries *series, ParamOfEvent *even_param )
{
  INT4  i,i_eval;
  REAL8 w_even;
  REAL4 ad,qd; /*abs of the datum and its square*/
  REAL4 s=0.;
  REAL4 ss=0.;  /*mean and average, without normalization, every i_eval samples*/
  REAL8 itaust;
  REAL8 norm=1.0;
  REAL8 asd_appo;

  if(series->data->length*series->deltaT < even_param->tau)
    even_param->tau = (REAL4) (1.0*series->data->length*series->deltaT);  /*change the memory time, if too long*/
  itaust = (REAL8) 1.0/(series->deltaT/even_param->tau);
  w_even = exp(-1.0/itaust);
  i_eval = (int) ((even_param->tau/series->deltaT)/even_param->factor); /*even_param->factor indicates that the computation of mean and std has to be re-evaluated after i_eval samples */
  if(i_eval<1)i_eval = 1;
  even_param->ar = 1;  /*The algorithm for the autoregressive evaluation is set.*/

  for(i=0; i<(int)series->data->length;i++){
    ad=fabs(series->data->data[i]); 
    qd=ad*ad;
    
    if(even_param->ar==1)even_param->xw=(1.0-w_even)*ad+w_even*even_param->xw;
    if(even_param->ar==1)even_param->qw=(1.0-w_even)*qd+w_even*even_param->qw;
    if(even_param->ar==1)even_param->w_norm=(1.0-w_even)+w_even*even_param->w_norm; 
    if(i%i_eval==0 && i!=0){ 
      s=even_param->xw;
      ss=even_param->qw;      
    } 
    if(even_param->w_norm !=0) {
      if(even_param->ar==1)
	norm=1.0/even_param->w_norm;
    }

    even_param->xamed[i]=s*norm;
    asd_appo=fabs(ss*norm-s*s*norm*norm);
    if (asd_appo !=0)
      even_param->xastd[i]=(REAL4) sqrt(asd_appo);
    else
      even_param->xastd[i]=0.0;
  }
  
  for(i=0; i<i_eval;i++){
     even_param->xamed[i]=even_param->xamed[i_eval];
     even_param->xastd[i]=even_param->xastd[i_eval];
   }
   
  /*Fill the empty values, with the previous ones (that is, the last evaluated)*/
    for(i=i_eval+1; i< (INT4) series->data->length;i++){
      if(i%i_eval!=0){
	even_param->xamed[i]=even_param->xamed[i-1];
	even_param->xastd[i]=even_param->xastd[i-1];
      }
    }

  return i;
}

int MedStd_dataSingle_R4( REAL4TimeSeries *series, ParamOfEvent *even_param )
{
  INT4  i,i_eval;
  REAL4 w_even;
  REAL4 ad,qd; /*abs of the datum and its square*/
  REAL4 s=0.;
  REAL4 ss=0.;  /*mean and average, without normalization, every i_eval samples*/
  REAL4 itaust;
  REAL4 norm=1.0;
  REAL4 asd_appo;

  if(series->data->length*series->deltaT < even_param->tau)
    even_param->tau = (REAL4) (1.0*series->data->length*series->deltaT);  /*change the memory time, if too long*/
  itaust = (REAL4) 1.0/(series->deltaT/even_param->tau);
  w_even = exp(-1.0/itaust);
  i_eval = (int) ((even_param->tau/series->deltaT)/even_param->factor); /*even_param->factor indicates that the computation of mean and std has to be re-evaluated after i_eval samples */
  if(i_eval<1)i_eval = 1;
  even_param->ar = 1;  /*The algorithm for the autoregressive evaluation is set.*/

  for(i=0; i<(int)series->data->length;i++){
    ad=fabs(series->data->data[i]); 
    qd=ad*ad;
    
    if(even_param->ar==1)even_param->xw=(1.0-w_even)*ad+w_even*even_param->xw;
    if(even_param->ar==1)even_param->qw=(1.0-w_even)*qd+w_even*even_param->qw;
    if(even_param->ar==1)even_param->w_norm=(1.0-w_even)+w_even*even_param->w_norm; 
    if(i%i_eval==0 && i!=0){ 
      s=even_param->xw;
      ss=even_param->qw;      
    } 
    if(even_param->w_norm !=0) {
      if(even_param->ar==1)
	norm=1.0/even_param->w_norm;
    }

    even_param->xamed[i]=s*norm;
    asd_appo=fabs(ss*norm-s*s*norm*norm);
    if (asd_appo !=0)
      even_param->xastd[i]=(REAL4) sqrt(asd_appo);
    else
      even_param->xastd[i]=0.0;
  }
  
  for(i=0; i<i_eval;i++){
     even_param->xamed[i]=even_param->xamed[i_eval];
     even_param->xastd[i]=even_param->xastd[i_eval];
   }
   
  /*Fill the empty values, with the previous ones (that is, the last evaluated)*/
    for(i=i_eval+1; i< (INT4) series->data->length;i++){
      if(i%i_eval!=0){
	even_param->xamed[i]=even_param->xamed[i-1];
	even_param->xastd[i]=even_param->xastd[i-1];
      }
    }

  return i;
}

/**********************************************************************************************************************************************/
/*COMPUTATION OF MEAN AND STANDARD DEVIATION for data in DOUBLE-PRECISION: the calculations are made by means of an autoregressive procedure*/
/**********************************************************************************************************************************************/
int MedStd_dataDouble( REAL8TimeSeries *seriesHPbil2, REAL8TimeSeries *series, ParamOfEvent *even_param )
{

  INT8 ts;
  int  i,i_eval;
  REAL8 w_even;
  /*REAL8 mean_e,se,sum_e,std_e,diff,mean,su,std,sum;*/
  REAL8 ad,qd; /*abs of the datum and its square*/
  REAL4 s=0.;
  REAL4 ss=0.;  /*mean and average, without normalization, every i_eval samples*/
  REAL8 itaust;
  REAL8 norm=1.0;
  REAL8 asd_appo;
  
  ts = ceil(even_param->tau/series->deltaT);
  seriesHPbil2->data->length = (series->data->length+ts);
  seriesHPbil2->deltaT = series->deltaT;

  /* if(series->data->length*series->deltaT < even_param->tau)even_param->tau=(REAL4) (1.0*series->data->length*series->deltaT); */ /*change the memory time, if too long*/
  itaust=(REAL8) 1.0/(series->deltaT/even_param->tau);
  w_even=exp(-1.0/itaust);
  i_eval= (int) ((even_param->tau/series->deltaT)/even_param->factor);/*even_param->factor indicates that the computation of mean and std has to be re-evaluated every i_eval samples */
  if(i_eval<1)i_eval=1;
  even_param->ar=1;  /*The algorithm for the autoregressive evaluation is set.*/

  for(i=0; i<(int)(series->data->length+ts);i++){
    ad=seriesHPbil2->data->data[i];
    qd=ad*ad;
    if(even_param->ar==1)even_param->xw=(1.0-w_even)*ad+w_even*even_param->xw;
    if(even_param->ar==1)even_param->qw=(1.0-w_even)*qd+w_even*even_param->qw;
    if(even_param->ar==1)even_param->w_norm=(1.0-w_even)+w_even*even_param->w_norm; 
    if(i%i_eval==0 && i!=0){ 
      s=even_param->xw;
      ss=even_param->qw;      
    } 
      if(even_param->w_norm !=0) {
      if(even_param->ar==1)norm=1.0/even_param->w_norm;
    }
      even_param->xamed[i]=s*norm;
      asd_appo=fabs(ss*norm-s*s*norm*norm);
      if(asd_appo !=0)even_param->xastd[i]=(REAL8) sqrt(asd_appo);
      else
	even_param->xastd[i]=0.0;
  }
  
  /* for(i=ts; i<i_eval;i++){*/
  for(i=0; i<i_eval;i++){
    even_param->xamed[i]=even_param->xamed[i_eval];
    even_param->xastd[i]=even_param->xastd[i_eval];
   }
   
  /*Fill the empty values, with the previous ones (that is, the last evaluated)*/
  for(i=i_eval+1; i<(int) (series->data->length+ts);i++){
    /*for(i=i_eval+1; i< (int) series->data->length;i++){*/
    if(i%i_eval!=0){
      even_param->xamed[i]=even_param->xamed[i-1];
      even_param->xastd[i]=even_param->xastd[i-1];
    }
  }
   
  /* OUTtest2=fopen("testsn2.dat","w");*/
  for(i=ts; i<(int)(series->data->length+ts);i++){
    even_param->med[i-ts]=even_param->xamed[i];
    even_param->std[i-ts]=even_param->xastd[i];
  } 

  return i;

}


int MedStd_R4( REAL4TimeSeries *seriesHPbil2, REAL4TimeSeries *series, ParamOfEvent *even_param )
{

  INT8 ts;
  int  i,i_eval;
  REAL4 w_even;
  /*REAL4 mean_e,se,sum_e,std_e,diff,mean,su,std,sum;*/
  REAL4 ad,qd; /*abs of the datum and its square*/
  REAL4 s=0.;
  REAL4 ss=0.;  /*mean and average, without normalization, every i_eval samples*/
  REAL4 itaust;
  REAL4 norm=1.0;
  REAL4 asd_appo;
  
  ts = ceil(even_param->tau/series->deltaT);
  seriesHPbil2->data->length = (series->data->length+ts);
  seriesHPbil2->deltaT = series->deltaT;

  /* if(series->data->length*series->deltaT < even_param->tau)even_param->tau=(REAL4) (1.0*series->data->length*series->deltaT); */ /*change the memory time, if too long*/
  itaust=(REAL4) 1.0/(series->deltaT/even_param->tau);
  w_even=exp(-1.0/itaust);
  i_eval= (int) ((even_param->tau/series->deltaT)/even_param->factor);/*even_param->factor indicates that the computation of mean and std has to be re-evaluated every i_eval samples */
  if(i_eval<1)i_eval=1;
  even_param->ar=1;  /*The algorithm for the autoregressive evaluation is set.*/

  for(i=0; i<(int)(series->data->length+ts);i++){
    ad=seriesHPbil2->data->data[i];
    qd=ad*ad;
    if(even_param->ar==1)even_param->xw=(1.0-w_even)*ad+w_even*even_param->xw;
    if(even_param->ar==1)even_param->qw=(1.0-w_even)*qd+w_even*even_param->qw;
    if(even_param->ar==1)even_param->w_norm=(1.0-w_even)+w_even*even_param->w_norm; 
    if(i%i_eval==0 && i!=0){ 
      s=even_param->xw;
      ss=even_param->qw;      
    } 
      if(even_param->w_norm !=0) {
      if(even_param->ar==1)norm=1.0/even_param->w_norm;
    }
      even_param->xamed[i]=s*norm;
      asd_appo=fabs(ss*norm-s*s*norm*norm);
      if(asd_appo !=0)even_param->xastd[i]=(REAL4) sqrt(asd_appo);
      else
	even_param->xastd[i]=0.0;
  }
  
  /* for(i=ts; i<i_eval;i++){*/
  for(i=0; i<i_eval;i++){
    even_param->xamed[i]=even_param->xamed[i_eval];
    even_param->xastd[i]=even_param->xastd[i_eval];
   }
   
  /*Fill the empty values, with the previous ones (that is, the last evaluated)*/
  for(i=i_eval+1; i<(int) (series->data->length+ts);i++){
    /*for(i=i_eval+1; i< (int) series->data->length;i++){*/
    if(i%i_eval!=0){
      even_param->xamed[i]=even_param->xamed[i-1];
      even_param->xastd[i]=even_param->xastd[i-1];
    }
  }
   
  /* OUTtest2=fopen("testsn2.dat","w");*/
  for(i=ts; i<(int)(series->data->length+ts);i++){
    even_param->med[i-ts]=even_param->xamed[i];
    even_param->std[i-ts]=even_param->xastd[i];
  } 

  return i;

}



/**************************************************************************************************************************/
/*EVENT-HUNTING: this function looks for events with a critical ratio above a 'fixed' threshold. It uses the mean and the 
standard deviation computed in an autoregressive way and it registers the event characteristics. DATA IN SINGLE-PRECISION*/
/**************************************************************************************************************************/
int EventSearch_dataSingle( REAL4TimeSeries *series, ParamOfEvent *even_param, BilHP *myparams )
{
  INT4 i,j;
  REAL4 ad,rd;
  REAL4 xmax = 0;
  INT8 time_below; /*Number of samples below*/
  INT8 samples_deadtime;
  REAL4 appo;
  appo=even_param->deadtime/series->deltaT;
  samples_deadtime=lroundf(appo);

  for(i=0;i<(int) series->data->length;i++){
    even_param->begin[i]=-1;
    even_param->duration[i]=0;
    even_param->imax[i]=0; 
    even_param->crmax[i]=0;
    even_param->ener[i]=0;
  }
  even_param->iflev=0;
  time_below=0.;
  j=-1;
  
  for(i=0; i<(int) series->data->length;i++){
#ifdef PSS_STYLE_RD
       ad=fabs(series->data->data[i]);
       rd=(ad-even_param->xamed[i])/(even_param->xastd[i]+1e-25);
#else
       ad=(series->data->data[i]);
       rd=fabs((ad-even_param->xamed[i])/(even_param->xastd[i]+1e-25));
#endif
         if((even_param->iflev==0)&&(rd>=myparams->crt) && (ad !=0)){  
	   j+=1;
	   even_param->iflev=1;
	   even_param->begin[j]=i;
	   even_param->duration[j]=1;
	   even_param->ener[j]+=ad*ad;
	   even_param->imax[j]=i;
	   even_param->crmax[j]=rd;
	   xmax=ad;
	   time_below=0;
         }
	 if((even_param->iflev==2)&&(rd>=myparams->crt)){
	   time_below=0;
	   even_param->duration[j]+=1;
	   even_param->ener[j]+=ad*ad;
	   if(ad>=xmax)even_param->imax[j]=i;
	   if(ad>=xmax)even_param->crmax[j]=rd;
	   if(ad>=xmax)xmax=ad;
	 }
	 if((even_param->iflev==2)&&(rd<myparams->crt)){
	   time_below+=1.;
	   even_param->duration[j]+=1;
	   even_param->ener[j]+=ad*ad;
	 }
	 if(time_below==samples_deadtime){
	   even_param->iflev=3;
	   time_below=0.;
	 } 
	 if(time_below>samples_deadtime){
	   even_param->duration[j]-=1;
	   even_param->iflev=3;
	   time_below=0.;
	 }
	 if(even_param->iflev==1)even_param->iflev=2;
	 if(even_param->iflev==3){
	   even_param->duration[j]-=(INT4) samples_deadtime; 
	    if(even_param->duration[j]<=0) even_param->duration[j]=1;    
	    even_param->iflev=0;
	 }
  }
  even_param->number=j+1;  /*Total number of events found (j starts from 0)*/
  return i;
}

/**************************************************************************************************************************/
/*EVENT-HUNTING: this function looks for events with a critical ratio above a 'fixed' threshold. It uses the mean and the 
standard deviation computed in an autoregressive way and it registers the event characteristics. DATA IN DOUBLE-PRECISION*/
/**************************************************************************************************************************/
int EventSearch_dataDouble( REAL8TimeSeries *series, ParamOfEvent *even_param, BilHP *myparams )
{

  /* REAL8 r = 60.0;
     INT8 ts, secstep;*/

  INT4 i,j,i_start;
  REAL8 ad,rd,rdend,rdst;
  REAL8 xmax = 0;
  long int time_below; /*Number of samples below*/
  long int samples_deadtime;
  REAL8 appo; /*REAL8 appo;*/
  appo=even_param->deadtime/series->deltaT;
  samples_deadtime=floor(appo+0.5); /*lround(appo);*/

  printf("deadtime, cr_th, fcut,samples_deadtime, deltaT, SerieLength = %f %f %f %ld %23.16e %i\n",even_param->deadtime,myparams->crt,myparams->fcut,samples_deadtime,series->deltaT,series->data->length);

  int itest=0;
  FILE *energy;
  FILE *ALLEVENT; 
  FILE *EVENDATone;
  FILE *EVENDATtwo;
  FILE *EVENDATthree;
  FILE *EVENDATfour;

  /*float numfl=1;
    double numdou=0.000215;
    float trial;
    double trial2;
    long int sd,sdf;
    trial=numfl/numdou;
    trial2=numfl/numdou;
    sdf=lroundf(trial);
    sd=floor(trial2+0.5); 
    printf("sdf sd: %ld %ld\n",sdf,sd);*/

  energy=fopen("energy.dat","w");
  for(i=0;i<(int) series->data->length;i++){
    even_param->begin[i]=-1;
    even_param->duration[i]=0;
    even_param->imax[i]=0; 
    even_param->crmax[i]=0;
    even_param->ener[i]=0;
  }
  even_param->iflev=0;
  time_below=0.;

  j=-1;

  if(itest) { 
    ALLEVENT=fopen("AllEvents.dat","w");
    EVENDATone=fopen("ARmean.dat","w");
    EVENDATtwo=fopen("ARstd.dat","w");
    EVENDATthree=fopen("NOAbsHp.dat","w");
    EVENDATfour=fopen("CR.dat","w");
  }
  
  for(i=0; i<(int) series->data->length;i++){
    ad=series->data->data[i];
    rd=fabs((ad-even_param->med[i])/(even_param->std[i]+1e-25));

    if(itest) { 
      fprintf(EVENDATone,"%23.16e \n",even_param->med[i]);
      fprintf(EVENDATtwo,"%23.16e \n",even_param->std[i]);
      fprintf(EVENDATthree,"%23.16e \n",series->data->data[i]);
      fprintf(EVENDATfour,"%23.16e \n",rd);
    }

    if((even_param->iflev==0) && (rd>=myparams->crt) && (ad !=0)){  
      j+=1;
      even_param->iflev=1;
      even_param->begin[j]=i;
      even_param->duration[j]=1;
      even_param->ener[j]+=ad*ad;
      even_param->imax[j]=i;
      even_param->crmax[j]=rd;
      xmax=ad;
      fprintf(energy,"%23.16e %23.16e %23.16e %23.16e %23.16e \n",even_param->ener[j],ad,even_param->crmax[j],even_param->med[j],even_param->std[j]);
      time_below=0;
    }

    if((even_param->iflev==2)&&(rd>=myparams->crt)){
      time_below=0;
      even_param->duration[j]+=1;
      even_param->ener[j]+=ad*ad;
      if(ad>=xmax)even_param->imax[j]=i;
      if(ad>=xmax)even_param->crmax[j]=rd;
      if(ad>=xmax)xmax=ad;
    }

    if((even_param->iflev==2)&&(rd<myparams->crt)){
      time_below+=1.;
      even_param->duration[j]+=1;
      even_param->ener[j]+=ad*ad;
    }

    if(time_below==samples_deadtime){
      even_param->iflev=3;
      time_below=0.;
    } 

    if(time_below>samples_deadtime){
      even_param->duration[j]-=1;
      even_param->iflev=3;
      /*puts("PAY ATTENTION: time_below is greater than samples_deadtime !!");*/
      time_below=0.;
    }

    if(even_param->iflev==1)
      even_param->iflev=2;
    if(even_param->iflev==3){
      even_param->duration[j]-=(int) samples_deadtime; 
      if(even_param->duration[j]<=0) even_param->duration[j]=1;    
      /*if(j==0)*/ 
      /*  printf("End of an event (j,i,dur,ad,med,std,rd) %i %i %i %23.16e %23.16e %23.16e %23.16e \n",j,i,even_param->duration[j],xmax,even_param->xamed[i],even_param->xastd[i],rd);*/
      /*EventNumber, Beg, Dur(samples), End, Dur(seconds), xmax, rd, mean, std;*/
        
      i_start=i-samples_deadtime-even_param->duration[j]+1;
      rdend=(xmax-even_param->med[i])/(even_param->std[i]+1e-25);
      rdst=(fabs(series->data->data[i_start])-even_param->med[i_start])/(even_param->std[i_start]+1e-25);

      even_param->iflev=0;

      if(itest)
	fprintf(ALLEVENT,"%i %i %i %i %23.16e %23.16e %23.16e %23.16e %23.16e %23.16e \n",
		j,even_param->begin[j],
		even_param->duration[j],
		even_param->imax[j],
		even_param->crmax[j],
		even_param->med[j],
		even_param->std[j],
		even_param->ener[j],
		ad*ad,xmax);
    }
  }

  even_param->number=j+1;  /*Total number of events found (j starts from 0)*/

  if(itest) { 
    fclose(EVENDATone);  
    fclose(EVENDATtwo); 
    fclose(EVENDATthree); 
    fclose(EVENDATfour); 
    fclose(ALLEVENT);
    fclose(energy);
  }

  return i;
}

/*******************************************************************************************************************************************/
/*ELIMINATION OF HIGH FREQUENCY SPIKES: this function identifies and removes large events from the highpassed time series. DATA IN SINGLE-PRECISION*/
/*******************************************************************************************************************************************/
int EventRemoval_dataSingle(REAL4TimeSeries *seriesCL, REAL4TimeSeries *series, REAL4TimeSeries *seriesHP, ParamOfEvent *even_param, BilHP *myparams)
{
  INT4 i,k;
  INT4 lwind;
  INT4 index1,index2,index3,index4;
  INT4 imax;
  
  i=MedStd_dataSingle_R4(seriesHP,even_param);
  i=EventSearch_dataSingle(seriesHP,even_param,myparams);

  /*copy in seriesCL the data of the original time series*/
  seriesCL->data->length=series->data->length; 
  seriesCL->deltaT=series->deltaT;
  for(i=0;i<(INT4)series->data->length;i++)
    seriesCL->data->data[i]=series->data->data[i];

  lwind=ceil(((REAL8) even_param->edge/series->deltaT));

  for(i=0; i<even_param->number; i++){
    if(even_param->begin[i]!=-1){
      imax=even_param->imax[i];

      index1=even_param->begin[i]-lwind;
      index2=even_param->begin[i]+even_param->duration[i]+lwind;
      index3=even_param->begin[i];
      index4=even_param->begin[i]+even_param->duration[i];
      if(index1 <= 0)
	index1=0;
      if(index2 >= (INT4)series->data->length)
	index2=series->data->length;

      /*In the event*/
      for(k=index3; k<index4; k++)
	seriesCL->data->data[k]-=seriesHP->data->data[k];
      /*At the edges:*/
      for(k=index1; k<index3; k++)seriesCL->data->data[k]-=seriesHP->data->data[k]*(k-1.0*index1)/lwind; 
      for(k=index4; k<index2; k++)seriesCL->data->data[k]-=seriesHP->data->data[k]*(1.0*index2-k)/lwind;
    } 
  }

  return i;
}

/*******************************************************************************************************************************************/
/*ELIMINATION OF HIGH FREQUENCY SPIKES: this function identifies and removes large events from the highpassed time series. DATA IN DOUBLE-PRECISION*/
/*******************************************************************************************************************************************/
int EventRemoval_dataDouble(REAL8TimeSeries *seriesCL, REAL8TimeSeries *series, REAL8TimeSeries *seriesHP, REAL8TimeSeries *seriesHPbil2, ParamOfEvent *even_param, BilHP *myparams)
/*int EventRemoval_dataDouble(REAL8TimeSeries *seriesCL, REAL8TimeSeries *series, REAL8TimeSeries *seriesHP, ParamOfEvent *even_param, BilHP *myparams)*/
{
  INT4 i,kk;
  INT4 lwind;
  INT4 index1,index2,index3,index4;
  INT4 imax;
  REAL8 su,su2,su3,cl;
  REAL8 st1,st2,st3,hp,wr,diffs;
  
  int itest = 0;
  FILE *outputtest;
  FILE *EnergyInfo;
  FILE *EventInfo;
  FILE *TM;
                
  cl = 0; 
  hp = 0;
  wr = 0;

  /* i = MedStd_dataDouble(seriesHP,even_param);*/
  i = MedStd_dataDouble(seriesHPbil2,seriesHP,even_param); 
  i = EventSearch_dataDouble(seriesHP,even_param,myparams);

  /*copy in seriesCL the data of  of the original time series*/
  seriesCL->data->length = series->data->length; 
  seriesCL->deltaT = series->deltaT;
  for(i = 0; i < (INT4)series->data->length; i++)
    seriesCL->data->data[i] = series->data->data[i];

  lwind = ceil(((REAL8) even_param->edge/series->deltaT));
  EventInfo = fopen("even_info.dat","w");
  EnergyInfo = fopen("energyInfo.dat","w");
  TM = fopen("mixT.dat","w");

  for(i = 0; i<even_param->number; i++){
    if(even_param->begin[i] != -1){
      imax=even_param->imax[i];

      index1 = even_param->begin[i]-lwind;
      index2 = even_param->begin[i]+even_param->duration[i]+lwind;
      index3 = even_param->begin[i];
      index4 = even_param->begin[i]+even_param->duration[i];
      if(index1 <= 0)
	index1 = 0;
      if(index2 >= (INT4)series->data->length)
	index2 = series->data->length;

      /*In the event*/
      for(kk = index3; kk<index4; kk++)
	seriesCL->data->data[kk] -= seriesHP->data->data[kk];
      /*At the edges:*/
      for(kk = index1; kk<index3; kk++)
	seriesCL->data->data[kk] -= seriesHP->data->data[kk] * (kk - 1.0 * index1) / lwind; 
      for(kk = index4; kk<index2; kk++)
	seriesCL->data->data[kk] -= seriesHP->data->data[kk] * (1.0 * index2 - kk) / lwind;
           
      fprintf(EventInfo,"%23.16e %23.16e %23.16e %i %23.16e %23.16e\n",
	      even_param->duration[i]*series->deltaT,
	      series->data->data[imax],
	      seriesHP->data->data[imax],
	      even_param->duration[i],
	      even_param->crmax[i],
	      even_param->ener[i]);

      su = 0;
      for(kk = index3; kk<index4; kk++){
	su += (seriesHP->data->data[kk])*(seriesHP->data->data[kk]);
      }
      st1 = 0;
      for(kk = index3; kk<index4; kk++)
	st1 += 2 * (seriesHP->data->data[kk]) * (series->data->data[kk]);
      
      su2 = 0;
      for(kk = index1; kk<index3; kk++)
	su2 += (seriesHP->data->data[kk] * (kk-1.0*index1) / lwind) * (seriesHP->data->data[kk] * (kk-1.0*index1) / lwind);
      st2 = 0;
      for(kk = index1; kk<index3; kk++)
	st2 +=  2 * (seriesHP->data->data[kk] * (kk - 1.0 * index1) / lwind) * (series->data->data[kk]); 
          
      su3 = 0;
      for(kk = index4; kk<index2; kk++)
	su3 += (seriesHP->data->data[kk] * (1.0 * index2 - kk) / lwind) * (seriesHP->data->data[kk] * (1.0*index2-kk) / lwind);
      st3 = 0;
      for(kk = index4; kk<index2; kk++)
	st3 += 2 * (seriesHP->data->data[kk] * (1.0 * index2-kk) / lwind) * (series->data->data[kk]);

      /* printf("%i %23.16e %23.16e %23.16e\n",kk,su,su2,su3); */
      fprintf(EnergyInfo,"%i %23.16e %23.16e %23.16e \n",i,su,su2,su3); 
      /* printf("%i %23.16e %23.16e %23.16e\n",kk,st1,st2,st3);*/
      fprintf(TM,"%i %23.16e %23.16e %23.16e\n",i,st1,st2,st3); 
        
    }
	 
  }
  fclose(EventInfo);
  fclose(EnergyInfo);
  fclose(TM);

  if(itest == 1){
    outputtest = fopen("CleanedTS.dat","w");
    for(i = 0; i < (INT4)series->data->length; i++)
      fprintf(outputtest,"%23.16e\n",seriesCL->data->data[i]);
    fclose(outputtest);      
  }

  for(i = 0;i < (INT4)series->data->length; i++){
    hp += seriesHP->data->data[i] * seriesHP->data->data[i];
    cl += seriesCL->data->data[i] * seriesCL->data->data[i];
    wr += series->data->data[i] * series->data->data[i];
  }

  diffs = wr-cl;
  printf("diffs %23.16e\n", diffs);
  printf("sumHP2 sumcl2 sumw2 %23.16e %23.16e %23.16e\n",hp,cl,wr);

  return i;
}

int EventRemoval_R4(REAL4TimeSeries *seriesCL, REAL4TimeSeries *series, REAL4TimeSeries *seriesHP, REAL4TimeSeries *seriesHPbil2, ParamOfEvent *even_param, BilHP *myparams)
{
  INT4 i,kk;
  INT4 lwind;
  INT4 index1,index2,index3,index4;
  INT4 imax;
  REAL4 su,su2,su3,cl;
  REAL4 st1,st2,st3,hp,wr,diffs;
  
  FILE *EnergyInfo;
  FILE *EventInfo;
  FILE *TM;
                
  cl = 0; 
  hp = 0;
  wr = 0;

  i = MedStd_R4(seriesHPbil2,seriesHP,even_param); 
  i = EventSearch_dataSingle(seriesHP,even_param,myparams);

  /*copy in seriesCL the data of  of the original time series*/
  seriesCL->data->length = series->data->length; 
  seriesCL->deltaT = series->deltaT;
  for(i = 0; i < (INT4)series->data->length; i++)
    seriesCL->data->data[i] = series->data->data[i];

  lwind = ceil(((REAL8) even_param->edge/series->deltaT));
  EventInfo = fopen("even_info.dat","w");
  EnergyInfo = fopen("energyInfo.dat","w");
  TM = fopen("mixT.dat","w");

  for(i = 0; i<even_param->number; i++){
    if(even_param->begin[i] != -1){
      imax=even_param->imax[i];

      index1 = even_param->begin[i]-lwind;
      index2 = even_param->begin[i]+even_param->duration[i]+lwind;
      index3 = even_param->begin[i];
      index4 = even_param->begin[i]+even_param->duration[i];
      if(index1 <= 0)
	index1 = 0;
      if(index2 >= (INT4)series->data->length)
	index2 = series->data->length;

      /*In the event*/
      for(kk = index3; kk<index4; kk++)
	seriesCL->data->data[kk] -= seriesHP->data->data[kk];
      /*At the edges:*/
      for(kk = index1; kk<index3; kk++)
	seriesCL->data->data[kk] -= seriesHP->data->data[kk] * (kk - 1.0 * index1) / lwind; 
      for(kk = index4; kk<index2; kk++)
	seriesCL->data->data[kk] -= seriesHP->data->data[kk] * (1.0 * index2 - kk) / lwind;
           
      fprintf(EventInfo,"%23.16e %23.16e %23.16e %i %23.16e %23.16e\n",
	      even_param->duration[i]*series->deltaT,
	      series->data->data[imax],
	      seriesHP->data->data[imax],
	      even_param->duration[i],
	      even_param->crmax[i],
	      even_param->ener[i]);

      su = 0;
      for(kk = index3; kk<index4; kk++){
	su += (seriesHP->data->data[kk])*(seriesHP->data->data[kk]);
      }
      st1 = 0;
      for(kk = index3; kk<index4; kk++)
	st1 += 2 * (seriesHP->data->data[kk]) * (series->data->data[kk]);
      
      su2 = 0;
      for(kk = index1; kk<index3; kk++)
	su2 += (seriesHP->data->data[kk] * (kk-1.0*index1) / lwind) * (seriesHP->data->data[kk] * (kk-1.0*index1) / lwind);
      st2 = 0;
      for(kk = index1; kk<index3; kk++)
	st2 +=  2 * (seriesHP->data->data[kk] * (kk - 1.0 * index1) / lwind) * (series->data->data[kk]); 
          
      su3 = 0;
      for(kk = index4; kk<index2; kk++)
	su3 += (seriesHP->data->data[kk] * (1.0 * index2 - kk) / lwind) * (seriesHP->data->data[kk] * (1.0*index2-kk) / lwind);
      st3 = 0;
      for(kk = index4; kk<index2; kk++)
	st3 += 2 * (seriesHP->data->data[kk] * (1.0 * index2-kk) / lwind) * (series->data->data[kk]);

      /* printf("%i %23.16e %23.16e %23.16e\n",kk,su,su2,su3); */
      fprintf(EnergyInfo,"%i %23.16e %23.16e %23.16e \n",i,su,su2,su3); 
      /* printf("%i %23.16e %23.16e %23.16e\n",kk,st1,st2,st3);*/
      fprintf(TM,"%i %23.16e %23.16e %23.16e\n",i,st1,st2,st3); 
        
    }
	 
  }
  fclose(EventInfo);
  fclose(EnergyInfo);
  fclose(TM);

  for(i = 0;i < (INT4)series->data->length; i++){
    hp += seriesHP->data->data[i] * seriesHP->data->data[i];
    cl += seriesCL->data->data[i] * seriesCL->data->data[i];
    wr += series->data->data[i] * series->data->data[i];
  }

  diffs = wr-cl;
  printf("diffs %23.16e\n", diffs);
  printf("sumHP2 sumcl2 sumw2 %23.16e %23.16e %23.16e\n",hp,cl,wr);

  return i;
}

/*******************************************************************************/

/**
  returns
  -1 if out of memory
  -2 if some computation failes
  -3 if input parameters are invalid
   0 otherwise (all went well)
*/
int PSSTDCleaningREAL8(REAL8TimeSeries *LALTS, REAL4 highpassFrequency) {
  UINT4 samples;                /**< number of samples in the timeseries */
  PSSTimeseries *originalTS;    /**< the timeseries converted to a PSS timeseries */
  PSSTimeseries *highpassTS;    /**< originalTS after high pass filtering */
  PSSTimeseries *cleanedTS;     /**< originalTS after cleaning */
  PSSEventParams *eventParams;  /**< keeps track of the "events" */
  PSSHeaderParams headerParams; /**< dummy, we don't actually use this to write SFTs */
  int retval = 0;               /**< return value of the function */

  fprintf(stderr,"[DEBUG] PSSTDCleaningREAL8 called\n");

  /* input sanity checks */
  if( !(LALTS) || !(LALTS->data) )
    return -3;


  /* number of samples in the original timeseries */
  samples = LALTS->data->length;

  /* reset errno before calling XLAL functions */
  xlalErrno = 0;

  /* open a log file, currently necessary for PSS */
  XLALPSSOpenLog("-");

  /* creation / memory allocation */
  /* there can't be more events than there are samples,
     so we prepare for as many events as we have samples */
  if( (eventParams = XLALCreatePSSEventParams(samples)) == NULL) {
    fprintf(stderr,"XLALCreatePSSEventParams call failed %s,%d\n",__FILE__,__LINE__);
    retval = -1;
    goto PSSTDCleaningREAL8FreeNothing;
  }
  if( (originalTS = XLALCreatePSSTimeseries(samples)) == NULL) {
    fprintf(stderr,"XLALCreatePSSTimeseries call failed %s,%d\n",__FILE__,__LINE__);
    retval = -1;
    goto PSSTDCleaningREAL8FreeEventParams;
  }
  if( (highpassTS = XLALCreatePSSTimeseries(samples)) == NULL) {
    fprintf(stderr,"XLALCreatePSSTimeseries call failed %s,%d\n",__FILE__,__LINE__);
    retval = -1;
    goto PSSTDCleaningREAL8FreeOriginalTS;
  }
  if( (cleanedTS = XLALCreatePSSTimeseries(samples)) == NULL) {
    fprintf(stderr,"XLALCreatePSSTimeseries call failed %s,%d\n",__FILE__,__LINE__);
    retval = -1;
    goto PSSTDCleaningREAL8FreeHighpassTS;
  }

  if (xlalErrno)
    fprintf(stderr,"PSSTDCleaningREAL8 (after alloc): unhandled XLAL Error %s,%d\n",__FILE__,__LINE__);

  /* the actual cleaning */
  if( XLALConvertREAL8TimeseriesToPSSTimeseries(originalTS, LALTS) == NULL) {
    fprintf(stderr,"XLALConvertREAL8TimeseriesToPSSTimeseries call failed %s,%d\n",__FILE__,__LINE__);
    retval = -2;
    goto PSSTDCleaningREAL8FreeAll;
  }
  if (xlalErrno)
    fprintf(stderr,"PSSTDCleaningREAL8 (after convert): unhandled XLAL Error %s,%d\n",__FILE__,__LINE__);

  XLALPrintREAL8TimeSeriesToFile(LALTS,"LALts.dat",50,-1);
  XLALPrintPSSTimeseriesToFile(originalTS,"originalTS.dat",50);
  if (xlalErrno)
    fprintf(stderr,"PSSTDCleaningREAL8 (after convert): unhandled XLAL Error %s,%d\n",__FILE__,__LINE__);

  if( XLALPSSHighpassData(highpassTS, originalTS, &headerParams, highpassFrequency) == NULL) {
    fprintf(stderr,"XLALPSSHighpassData call failed %s,%d\n",__FILE__,__LINE__);
    retval = -2;
    goto PSSTDCleaningREAL8FreeAll;
  }

  XLALPrintPSSTimeseriesToFile(highpassTS,"highpassTS.dat",50);

  if( XLALIdentifyPSSCleaningEvents(eventParams, highpassTS, &headerParams) == NULL) {
    fprintf(stderr,"XLALIdentifyPSSCleaningEvents call failed %s,%d\n",__FILE__,__LINE__);
    retval = -2;
    goto PSSTDCleaningREAL8FreeAll;
  }
  if( XLALSubstractPSSCleaningEvents(cleanedTS, originalTS, highpassTS, eventParams, &headerParams) == NULL) {
    fprintf(stderr,"XLALSubstractPSSCleaningEvents call failed %s,%d\n",__FILE__,__LINE__);
    retval = -2;
    goto PSSTDCleaningREAL8FreeAll;
  }

  XLALPrintPSSTimeseriesToFile(cleanedTS,"cleanedTS.dat",0);

  if (xlalErrno)
    fprintf(stderr,"PSSTDCleaningREAL8 (before convert): unhandled XLAL Error %s,%d\n",__FILE__,__LINE__);

  if( XLALConvertPSSTimeseriesToREAL8Timeseries(LALTS, cleanedTS) == NULL) {
    fprintf(stderr,"XLALConvertPSSTimeseriesToREAL8Timeseries call failed %s,%d\n",__FILE__,__LINE__);
    retval = -2;
    goto PSSTDCleaningREAL8FreeAll;
  }

  if (xlalErrno)
    fprintf(stderr,"PSSTDCleaningREAL8 (after PSS): unhandled XLAL Error %s,%d\n",__FILE__,__LINE__);

  /* debug: write out autoregressive mean and std */
  PrintREAL4ArrayToFile( "PSS_ARmed.dat", eventParams->xamed, dataDouble.data->length );
  PrintREAL4ArrayToFile( "PSS_ARstd.dat", eventParams->xastd, dataDouble.data->length );

  /* cleanup & return */
 PSSTDCleaningREAL8FreeAll:
  XLALDestroyPSSTimeseries(cleanedTS);
  fprintf(stderr, "[DEBUG] XLALDestroyPSSTimeseries done.\n");
 PSSTDCleaningREAL8FreeHighpassTS:
  XLALDestroyPSSTimeseries(highpassTS);
  fprintf(stderr, "[DEBUG] XLALDestroyPSSTimeseries done.\n");
 PSSTDCleaningREAL8FreeOriginalTS:
  XLALDestroyPSSTimeseries(originalTS);
  fprintf(stderr, "[DEBUG] XLALDestroyPSSTimeseries done.\n");
 PSSTDCleaningREAL8FreeEventParams:
  XLALDestroyPSSEventParams(eventParams);
  fprintf(stderr, "[DEBUG] XLALDestroyPSSEventParams done.\n");
 PSSTDCleaningREAL8FreeNothing:
  XLALPSSCloseLog();
  fprintf(stderr, "[DEBUG] XLALPSSCloseLog done.\n");

  if (xlalErrno)
    fprintf(stderr,"PSSTDCleaningREAL8 (after free()): unhandled XLAL Error %s,%d\n",__FILE__,__LINE__);

  if (retval)
    fprintf(stderr,"PSSTDCleaningREAL8 nonzero retval %d\n", retval);

  return retval;
}


int PSSTDCleaningREAL4(REAL4TimeSeries *LALTS, REAL4 highpassFrequency) {
  UINT4 samples;                /**< number of samples in the timeseries */
  PSSTimeseries *originalTS;    /**< the timeseries converted to a PSS timeseries */
  PSSTimeseries *highpassTS;    /**< originalTS after high pass filtering */
  PSSTimeseries *cleanedTS;     /**< originalTS after cleaning */
  PSSEventParams *eventParams;  /**< keeps track of the "events" */
  PSSHeaderParams headerParams; /**< dummy, we don't actually use this to write SFTs */
  int retval = 0;               /**< return value of the function */

  /* input sanity checks */
  if( !(LALTS) || !(LALTS->data) )
    return -3;

  /* number of samples in the original timeseries */
  samples = LALTS->data->length;

  /* open a log file, currently necessary for PSS */
  XLALPSSOpenLog("-");

  /* creation / memory allocation */
  /* there can't be more events than there are samples,
     so we prepare for as many events as we have samples */
  if( (eventParams = XLALCreatePSSEventParams(samples)) == NULL) {
    retval = -1;
    goto PSSTDCleaningREAL4FreeNothing;
  }
  if( (originalTS = XLALCreatePSSTimeseries(samples)) == NULL) {
    retval = -1;
    goto PSSTDCleaningREAL4FreeEventParams;
  }
  if( (highpassTS = XLALCreatePSSTimeseries(samples)) == NULL) {
    retval = -1;
    goto PSSTDCleaningREAL4FreeOriginalTS;
  }
  if( (cleanedTS = XLALCreatePSSTimeseries(samples)) == NULL) {
    retval = -1;
    goto PSSTDCleaningREAL4FreeHighpassTS;
  }


  /* the actual cleaning */
  if( XLALConvertREAL4TimeseriesToPSSTimeseries(originalTS, LALTS) == NULL) {
    retval = -2;
    goto PSSTDCleaningREAL4FreeAll;
  }
  if( XLALPSSHighpassData(highpassTS, originalTS, &headerParams, highpassFrequency) == NULL) {
    retval = -2;
    goto PSSTDCleaningREAL4FreeAll;
  }
  if( XLALIdentifyPSSCleaningEvents(eventParams, highpassTS, &headerParams) == NULL) {
    retval = -2;
    goto PSSTDCleaningREAL4FreeAll;
  }
  if( XLALSubstractPSSCleaningEvents(cleanedTS, originalTS, highpassTS, eventParams, &headerParams) == NULL) {
    retval = -2;
    goto PSSTDCleaningREAL4FreeAll;
  }
  if( XLALConvertPSSTimeseriesToREAL4Timeseries(LALTS, cleanedTS) == NULL) {
    retval = -2;
    goto PSSTDCleaningREAL4FreeAll;
  }


  /* cleanup & return */
 PSSTDCleaningREAL4FreeAll:
  XLALDestroyPSSTimeseries(cleanedTS);
 PSSTDCleaningREAL4FreeHighpassTS:
  XLALDestroyPSSTimeseries(highpassTS);
 PSSTDCleaningREAL4FreeOriginalTS:
  XLALDestroyPSSTimeseries(originalTS);
 PSSTDCleaningREAL4FreeEventParams:
  XLALDestroyPSSEventParams(eventParams);
 PSSTDCleaningREAL4FreeNothing:
  XLALPSSCloseLog();
  return retval;
}

int PSSTDCleaningSingle(struct CommandLineArgsTag CLA) {
  return(PSSTDCleaningREAL4(&dataSingle, CLA.fc));
}

int PSSTDCleaningDouble(struct CommandLineArgsTag CLA) {
  return(PSSTDCleaningREAL8(&dataDouble, CLA.fc));
}

int TDCleaning(struct CommandLineArgsTag CLA)
{
  /* FIXME: Memory leak even_param */
  INT4 j, k;
  ParamOfEvent *even_param;
  even_param=(ParamOfEvent *)calloc(1,sizeof(ParamOfEvent)); 

  BilHP bilparam;

  bilparam.fcut    = CLA.fc;
  bilparam.crt     = CLA.cr;
  
  EventParamInit(dataDouble.data->length, even_param);
  j=BilHighpass_dataDouble( &dataDoubleHP, &dataDoubleFirstHP, &dataDouble, even_param, &bilparam );

  XLALPrintREAL8TimeSeriesToFile(&dataDouble,"dataDouble.dat",50,-1);
  XLALPrintREAL8TimeSeriesToFile(&dataDoubleFirstHP,"dataDoubleFirstHP.dat",50,-1);
  XLALPrintREAL8TimeSeriesToFile(&dataDoubleHP,"dataDoubleHP.dat",50,-1);

  Evenbil2(even_param, &dataDouble );
  j=Bil2( &databil2, &dataDoubleHP, even_param);
  j=EventRemoval_dataDouble( &dataDoubleClean, &dataDouble, &dataDoubleHP, &databil2, even_param, &bilparam);
  dataDouble.data->length=dataDoubleClean.data->length;   
  dataDouble.deltaT=dataDoubleClean.deltaT;

  for (k = 0; k < (int)dataDoubleClean.data->length; k++) {
    dataDouble.data->data[k] = dataDoubleClean.data->data[k];
  }

  return 0;
}


int TDCleaning_R4(struct CommandLineArgsTag CLA)
{       
  INT4 k;
  ParamOfEvent even_param;
  BilHP bilparam;

  bilparam.fcut    = CLA.fc;
  bilparam.crt     = CLA.cr;
  
  EventParamInit(dataDouble.data->length, &even_param);

  /* debug: write out original dataDouble TS */
  XLALPrintREAL8TimeSeriesToFile(&dataDouble,"dataDouble.dat",50,-1);

  /* debug: convert dataDouble into dataSingle TS */
  XLALConvertREAL8TimeSeriesToREAL4TimeSeries(&dataSingle, &dataDouble);

  /* debug: write out dataSingle to be sure the conversion works */
  XLALPrintREAL4TimeSeriesToFile(&dataSingle,"dataSingle.dat",50);

  /* highpass the data using single precision */
  BilHighpass_dataSingle( &dataSingleHP, &dataSingleFirstHP, &dataSingle, &even_param, &bilparam );

  /* debug: write out the highpass TS */
  XLALPrintREAL4TimeSeriesToFile(&dataSingleFirstHP,"dataSingleFirstHP.dat",50);
  XLALPrintREAL4TimeSeriesToFile(&dataSingleHP,"dataSingleHP.dat",50);

  /* REAL4TimeSeries version of Evenbil2() */
  Evenbil2_R4(&even_param, &dataSingle );

  /* REAL4TimeSeries version of Bil2() */
  Bil2_R4( &databil2Single, &dataSingleHP, &even_param);
  
  EventRemoval_R4( &dataSingleClean, &dataSingle, &dataSingleHP, &databil2Single, &even_param, &bilparam);

  /* write out the cleaned data */
  XLALPrintREAL4TimeSeriesToFile(&dataSingleClean,"dataSingleClean.dat",0);

  /* copy the cleaned data back to dataDouble */
  for (k = 0; k < (int)dataSingleClean.data->length; k++) {
    dataDouble.data->data[k] = dataSingleClean.data->data[k];
  }

  return 0;
}


int TDCleaning_NoBil_R4(struct CommandLineArgsTag CLA)
{       
  INT4 k;
  ParamOfEvent even_param;
  BilHP bilparam;

  bilparam.fcut    = CLA.fc;
  bilparam.crt     = CLA.cr;
  
  EventParamInit_NoBil(dataDouble.data->length, &even_param);

  /* debug: write out original dataDouble TS */
  XLALPrintREAL8TimeSeriesToFile(&dataDouble,"NoBilDataDouble.dat",50,-1);

  /* debug: convert dataDouble into dataSingle TS */
  XLALConvertREAL8TimeSeriesToREAL4TimeSeries(&dataSingle, &dataDouble);

  /* debug: write out dataSingle to be sure the conversion works */
  XLALPrintREAL4TimeSeriesToFile(&dataSingle,"NoBilDataSingle.dat",50);

  /* highpass the data using single precision */
  BilHighpass_dataSingle( &dataSingleHP, &dataSingleFirstHP, &dataSingle, &even_param, &bilparam );

  /* debug: write out the highpass TS */
  XLALPrintREAL4TimeSeriesToFile(&dataSingleFirstHP,"NoBilDataSingleFirstHP.dat",50);
  XLALPrintREAL4TimeSeriesToFile(&dataSingleHP,"NoBilDataSingleHP.dat",50);

  EventRemoval_dataSingle( &dataSingleClean, &dataSingle, &dataSingleHP, &even_param, &bilparam);

  /* debug: write out autoregressive mean and std */
  PrintREAL4ArrayToFile( "NoBil_ARmed.dat", even_param.xamed, dataDouble.data->length );
  PrintREAL4ArrayToFile( "NoBil_ARstd.dat", even_param.xastd, dataDouble.data->length );

  /* write out the cleaned data */
  XLALPrintREAL4TimeSeriesToFile(&dataSingleClean,"NoBilDataSingleClean.dat",0);

  /* copy the cleaned data back to dataDouble */
  for (k = 0; k < (int)dataSingleClean.data->length; k++) {
    dataDouble.data->data[k] = dataSingleClean.data->data[k];
  }

  return 0;
}


/*******************************************************************************/

/*******************************************************************************/
int CreateSFT(struct CommandLineArgsTag CLA)
{
  /*INT4 i;*/
  /* 11/19/05 gam */
  if(CLA.useSingle) {
   
      #if TRACKMEMUSE
        printf("Memory use before creating output vector fftDataSingle and calling LALForwardRealFFT:\n"); printmemuse();
      #endif

      /* 11/02/05 gam; allocate container for SFT data */
      LALCCreateVector( &status, &fftDataSingle, dataSingle.data->length / 2 + 1 );
      TESTSTATUS( &status );  

      /* compute sft */
      LALForwardRealFFT(&status, fftDataSingle, dataSingle.data, fftPlanSingle );
      TESTSTATUS( &status );      

      #if TRACKMEMUSE
        printf("Memory use after creating output vector fftDataSingle and calling LALForwardRealFFT:\n"); printmemuse();
      #endif

      #if PRINTEXAMPLEDATA
        printf("\nExample real and imaginary value of fftDataSingle from CreateSFT:\n"); printExampleFFTData(CLA);
      #endif      
  
  } else {

      #if TRACKMEMUSE
        printf("Memory use before creating output vector fftDataDouble and calling XLALREAL8ForwardFFT:\n"); printmemuse();
      #endif

      /* 11/02/05 gam; allocate container for SFT data */
      LALZCreateVector( &status, &fftDataDouble, dataDouble.data->length / 2 + 1 );
      TESTSTATUS( &status );  

      /* compute sft */
      XLALREAL8ForwardFFT( fftDataDouble, dataDouble.data, fftPlanDouble );

      #if TRACKMEMUSE
        printf("Memory use after creating output vector fftDataDouble and calling XLALREAL8ForwardFFT:\n"); printmemuse();
      #endif

      #if PRINTEXAMPLEDATA
        printf("\nExample real and imaginary value of fftDataDouble from CreateSFT:\n"); printExampleFFTData(CLA);
      #endif      
 
  }
      return 0;
}
/*******************************************************************************/

/*******************************************************************************/
int WriteSFT(struct CommandLineArgsTag CLA)
{
  char sftname[256];
  char sftnameFinal[256]; /* 01/09/06 gam */
  char numSFTs[2]; /* 12/27/05 gam */
  char site[2];    /* 12/27/05 gam */
  char ifo[3];     /* 12/27/05 gam; allow 3rd charactor for null termination */
  int firstbin=(INT4)(FMIN*CLA.T+0.5), k;
  FILE *fpsft;
  int errorcode1=0,errorcode2=0;
  char gpstime[11]; /* 12/27/05 gam; allow for 10 digit GPS times and null termination */
  REAL4 rpw,ipw;

  /* 12/27/05 gam; set up the number of SFTs, site, and ifo as null terminated strings */
  numSFTs[0] = '1';
  numSFTs[1] = '\0'; /* null terminate */
  strncpy( site, CLA.ChannelName, 1 );
  site[1] = '\0'; /* null terminate */
  if (CLA.IFO != NULL) {
    strncpy( ifo, CLA.IFO, 2 );
  } else {
    strncpy( ifo, CLA.ChannelName, 2 );
  }
  ifo[2] = '\0'; /* null terminate */
  sprintf(gpstime,"%09d",gpsepoch.gpsSeconds);

  strcpy( sftname, CLA.SFTpath );
  /* 12/27/05 gam; add option to make directories based on gps time */
  if (CLA.makeGPSDirs > 0) {
     /* 12/27/05 gam; concat to the sftname the directory name based on GPS time; make this directory if it does not already exist */
     mkSFTDir(sftname, site, numSFTs, ifo, CLA.stringT, CLA.miscDesc, gpstime, CLA.makeGPSDirs);
  }

  /* 01/09/06 gam; sftname will be temporary; will move to sftnameFinal. */
  if(CLA.makeTmpFile) {
    /* set up sftnameFinal with usual SFT name */
    strcpy(sftnameFinal,sftname);
    strcat(sftnameFinal,"/SFT_");
    strncat(sftnameFinal,ifo,2);
    strcat(sftnameFinal,".");
    strcat(sftnameFinal,gpstime);
    /* sftname begins with . and ends in .tmp */
    strcat(sftname,"/.SFT_");
    strncat(sftname,ifo,2);
    strcat(sftname,".");
    strcat(sftname,gpstime);
    strcat(sftname,".tmp");
  } else {     
    strcat(sftname,"/SFT_");
    strncat(sftname,ifo,2);
    strcat(sftname,".");
    strcat(sftname,gpstime);
  }  

  /* open SFT file for writing */
  fpsft=tryopen(sftname,"w");

  /* write header */
  header.endian=1.0;
  header.gps_sec=gpsepoch.gpsSeconds;
  header.gps_nsec=gpsepoch.gpsNanoSeconds;
  header.tbase=CLA.T;
  header.firstfreqindex=firstbin;
  header.nsamples=(INT4)(DF*CLA.T+0.5);
  if (1!=fwrite((void*)&header,sizeof(header),1,fpsft)){
    fprintf(stderr,"Error in writing header into file %s!\n",sftname);
    return 7;
  }

  /* 11/19/05 gam */
  if(CLA.useSingle) {
    /* Write SFT */
    for (k=0; k<header.nsamples; k++)
    {
      rpw=((REAL4)(((REAL8)DF)/(0.5*(REAL8)(1/dataSingle.deltaT))))
	* fftDataSingle->data[k+firstbin].re;
      ipw=((REAL4)(((REAL8)DF)/(0.5*(REAL8)(1/dataSingle.deltaT))))
	* fftDataSingle->data[k+firstbin].im;
      /* 06/26/07 gam; use finite to check that data does not contains a non-FINITE (+/- Inf, NaN) values */
      #if CHECKFORINFINITEANDNANS
        if (!finite(rpw) || !finite(ipw)) {
          fprintf(stderr, "Infinite or NaN data at freq bin %d.\n", k);
          return 7;
        }
      #endif
      errorcode1=fwrite((void*)&rpw, sizeof(REAL4),1,fpsft);
      errorcode2=fwrite((void*)&ipw, sizeof(REAL4),1,fpsft);
    }
    #if PRINTEXAMPLEDATA
        printf("\nExample real and imaginary SFT values going to file from fftDataSingle in WriteSFT:\n"); printExampleSFTDataGoingToFile(CLA);
    #endif     
    LALCDestroyVector( &status, &fftDataSingle );
    TESTSTATUS( &status );
    #if TRACKMEMUSE
      printf("Memory use after destroying fftDataSingle in WriteSFT:\n"); printmemuse();
    #endif
  } else {    
    /* Write SFT */
    for (k=0; k<header.nsamples; k++)
    {
      rpw=(((REAL8)DF)/(0.5*(REAL8)(1/dataDouble.deltaT))) 
	* fftDataDouble->data[k+firstbin].re;
      ipw=(((REAL8)DF)/(0.5*(REAL8)(1/dataDouble.deltaT))) 
	* fftDataDouble->data[k+firstbin].im;
      /* 06/26/07 gam; use finite to check that data does not contains a non-FINITE (+/- Inf, NaN) values */
      #if CHECKFORINFINITEANDNANS
        if (!finite(rpw) || !finite(ipw)) {
          fprintf(stderr, "Infinite or NaN data at freq bin %d.\n", k);
          return 7;
        }
      #endif
      errorcode1=fwrite((void*)&rpw, sizeof(REAL4),1,fpsft);
      errorcode2=fwrite((void*)&ipw, sizeof(REAL4),1,fpsft);
    }
    #if PRINTEXAMPLEDATA
        printf("\nExample real and imaginary SFT values going to file from fftDataDouble in WriteSFT:\n"); printExampleSFTDataGoingToFile(CLA);
    #endif
    LALZDestroyVector( &status, &fftDataDouble );
    TESTSTATUS( &status );
    #if TRACKMEMUSE
      printf("Memory use after destroying fftDataDouble in WriteSFT:\n"); printmemuse();
    #endif
  }

  /* Check that there were no errors while writing SFTS */
  if (errorcode1-1 || errorcode2-1)
    {
      fprintf(stderr,"Error in writing data into SFT file %s!\n",sftname);
      return 8;
    }
  
  fclose(fpsft);

  /* 01/09/06 gam; sftname is temporary; move to sftnameFinal. */
  if(CLA.makeTmpFile) {  
    mvFilenames(sftname,sftnameFinal);
  }

  return 0;
}
/*******************************************************************************/

/*******************************************************************************/
/* 12/28/05 gam; write out version 2 SFT */
int WriteVersion2SFT(struct CommandLineArgsTag CLA) 
{
  char sftname[256];
  char sftFilename[256];
  char sftnameFinal[256]; /* 01/09/06 gam */
  char numSFTs[2]; /* 12/27/05 gam */
  char site[2];    /* 12/27/05 gam */
  char ifo[3];     /* 12/27/05 gam; allow 3rd charactor for null termination */
  int firstbin=(INT4)(FMIN*CLA.T+0.5), k;
  char gpstime[11]; /* 12/27/05 gam; allow for 10 digit GPS times and null termination */
  SFTtype *oneSFT = NULL;
  INT4 nBins = (INT4)(DF*CLA.T+0.5);
  REAL4 singleDeltaT = 0.0; /* 01/05/06 gam */
  REAL8 doubleDeltaT = 0.0; /* 01/05/06 gam */

  REAL8 sureim;
   /* 12/27/05 gam; set up the number of SFTs, site, and ifo as null terminated strings */
  numSFTs[0] = '1';
  numSFTs[1] = '\0'; /* null terminate */
  strncpy( site, CLA.ChannelName, 1 );
  site[1] = '\0'; /* null terminate */
  if (CLA.IFO != NULL) {
    strncpy( ifo, CLA.IFO, 2 );
  } else {  
    strncpy( ifo, CLA.ChannelName, 2 );
  }
  ifo[2] = '\0'; /* null terminate */
  sprintf(gpstime,"%09d",gpsepoch.gpsSeconds);

  strcpy( sftname, CLA.SFTpath );
  /* 12/27/05 gam; add option to make directories based on gps time */
  if (CLA.makeGPSDirs > 0) {
     /* 12/27/05 gam; concat to the sftname the directory name based on GPS time; make this directory if it does not already exist */
     mkSFTDir(sftname, site, numSFTs, ifo, CLA.stringT, CLA.miscDesc, gpstime, CLA.makeGPSDirs);
  }

  strcat(sftname,"/");
  mkSFTFilename(sftFilename, site, numSFTs, ifo, CLA.stringT, CLA.miscDesc, gpstime);
  /* 01/09/06 gam; sftname will be temporary; will move to sftnameFinal. */
  if(CLA.makeTmpFile) {
    /* set up sftnameFinal with usual SFT name */
    strcpy(sftnameFinal,sftname);
    strcat(sftnameFinal,sftFilename);
    /* sftname begins with . and ends in .tmp */
    strcat(sftname,".");
    strcat(sftname,sftFilename);
    strcat(sftname,".tmp");
  } else {
    strcat(sftname,sftFilename);
  }  

  /* make container to store the SFT data */
  LALCreateSFTtype (&status, &oneSFT, ((UINT4)nBins));
  TESTSTATUS( &status );
  #if TRACKMEMUSE
      printf("Memory use after creating oneSFT and calling LALCreateSFTtype in WriteVersion2SFT:\n"); printmemuse();
  #endif
  
  /* copy the data to oneSFT */
  strcpy(oneSFT->name,ifo);
  oneSFT->epoch.gpsSeconds=gpsepoch.gpsSeconds;
  oneSFT->epoch.gpsNanoSeconds=gpsepoch.gpsNanoSeconds;
  oneSFT->f0 = FMIN;
  oneSFT->deltaF = 1.0/((REAL8)CLA.T);
  oneSFT->data->length=nBins;
  
  if(CLA.useSingle) {
    singleDeltaT = ((REAL4)dataSingle.deltaT); /* 01/05/06 gam; and normalize SFTs using this below */
    for (k=0; k<nBins; k++)
    {
      oneSFT->data->data[k].re = singleDeltaT*fftDataSingle->data[k+firstbin].re;
      oneSFT->data->data[k].im = singleDeltaT*fftDataSingle->data[k+firstbin].im;
      /* 06/26/07 gam; use finite to check that data does not contains a non-FINITE (+/- Inf, NaN) values */
      #if CHECKFORINFINITEANDNANS
        if (!finite(oneSFT->data->data[k].re) || !finite(oneSFT->data->data[k].im)) {
          fprintf(stderr, "Infinite or NaN data at freq bin %d.\n", k);
          return 7;
        }
      #endif      
    }
    #if PRINTEXAMPLEDATA
        printf("\nExample real and imaginary SFT values going to file from fftDataSingle in WriteVersion2SFT:\n"); printExampleVersion2SFTDataGoingToFile(CLA,oneSFT);
    #endif
    LALCDestroyVector( &status, &fftDataSingle );
    TESTSTATUS( &status );
    #if TRACKMEMUSE
      printf("Memory use after destroying fftDataSingle in WriteVersion2SFT:\n"); printmemuse();
    #endif
  } else {
    doubleDeltaT = ((REAL8)dataDouble.deltaT); /* 01/05/06 gam; and normalize SFTs using this below */
  
    sureim=0;
    for (k=0; k<nBins; k++)
    {
      oneSFT->data->data[k].re = doubleDeltaT*fftDataDouble->data[k+firstbin].re;
      oneSFT->data->data[k].im = doubleDeltaT*fftDataDouble->data[k+firstbin].im;
      /* 06/26/07 gam; use finite to check that data does not contains a non-FINITE (+/- Inf, NaN) values */
      sureim+=(oneSFT->data->data[k].re)*(oneSFT->data->data[k].re)+(oneSFT->data->data[k].im)*(oneSFT->data->data[k].im); 
      
      #if CHECKFORINFINITEANDNANS
        if (!finite(oneSFT->data->data[k].re) || !finite(oneSFT->data->data[k].im)) {
          fprintf(stderr, "Infinite or NaN data at freq bin %d.\n", k);
          return 7;
        }
      #endif
    }
    printf("re2+im2Sum %23.26e \n",sureim); 
   
    #if PRINTEXAMPLEDATA
        printf("\nExample real and imaginary SFT values going to file from fftDataDouble in WriteVersion2SFT:\n"); printExampleVersion2SFTDataGoingToFile(CLA,oneSFT);
    #endif
    LALZDestroyVector( &status, &fftDataDouble );
    TESTSTATUS( &status );
    #if TRACKMEMUSE
      printf("Memory use after destroying fftDataDouble in WriteVersion2SFT:\n"); printmemuse();
    #endif
  }  

  /* write the SFT */
  LALWriteSFT2file(&status, oneSFT, sftname, CLA.commentField);
  TESTSTATUS( &status );

  /* 01/09/06 gam; sftname is temporary; move to sftnameFinal. */
  if(CLA.makeTmpFile) {  
    mvFilenames(sftname,sftnameFinal);
  }

  LALDestroySFTtype (&status,&oneSFT);
  TESTSTATUS( &status );
  #if TRACKMEMUSE
      printf("Memory use after destroying oneSFT and calling LALDestroySFTtype in WriteVersion2SFT:\n"); printmemuse();
  #endif

  return 0;
}
/*******************************************************************************/

/*******************************************************************************/
int FreeMem(struct CommandLineArgsTag CLA)
{

  LALFrClose(&status,&framestream);
  TESTSTATUS( &status );

  /* 11/19/05 gam */
  if(CLA.useSingle || (CLA.TDcleaningProc >= 3)) {

    LALDestroyVector(&status,&dataSingle.data);
    TESTSTATUS( &status );

    if (CommandLineArgs.TDcleaningProc > 0) {  
    
      LALDestroyVector(&status,&dataSingleFirstHP.data);
      TESTSTATUS( &status );

      LALDestroyVector(&status,&dataSingleHP.data);
      TESTSTATUS( &status ); 

      LALDestroyVector(&status,&dataSingleClean.data);
      TESTSTATUS( &status );
    }


    LALDestroyRealFFTPlan( &status, &fftPlanSingle );
    TESTSTATUS( &status );
     
  }

  if(!CLA.useSingle) {

    LALDDestroyVector(&status,&dataDouble.data);
    TESTSTATUS( &status );

    if ((CLA.TDcleaningProc > 0) && (CLA.TDcleaningProc < 3)) {  
    
      LALDDestroyVector(&status,&dataDoubleFirstHP.data);
      TESTSTATUS( &status );

      LALDDestroyVector(&status,&dataDoubleHP.data);
      TESTSTATUS( &status ); 

      LALDDestroyVector(&status,&databil2.data);
      TESTSTATUS( &status ); 

      LALDDestroyVector(&status,&dataDoubleClean.data);
      TESTSTATUS( &status );
    }
    
    LALDestroyREAL8FFTPlan( &status, &fftPlanDouble );
    TESTSTATUS( &status );
    
  }

  LALCheckMemoryLeaks();
 
  return 0;
}
/*******************************************************************************/

#endif /* #if !defined HAVE_LIBGSL || !defined HAVE_LIBLALFRAME */
