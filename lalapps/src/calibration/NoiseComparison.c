/*
*  Copyright (C) 2007 Xavier Siemens
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
/* Noise comparison code between h(t) and frequency domain calibrated DARM_ERR   */
/*                                                                               */
/*                                  X. Siemens                                   */
/*                                                                               */
/*                           Caltech -- December 2006                            */
/*********************************************************************************/

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
#include <lal/ConfigFile.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/LALFrameL.h>

#include <series.h>



/* #define DEBUG 1 */  /* Uncomment this line to enable debugging fprintfs and mem usage */

double hypot(double, double);

extern char *optarg;
extern int optind, opterr, optopt;

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)

COMPLEX8 tmpa, tmpb, tmpc;
REAL4 tmpx, tmpy;

#define MAXLINESRS   60000     /* Maximum # of lines in a Response file */
#define MAXFACTORS   100000       /* Maximum # of factors to be computed */
#define MAXFREQUENCIES 200       /* Maximum number of frequemcies for which to do the comparison */
#define DECIMATE 4 /*Factor by which to decimate the data when the -D (decimate) flag is used */


#define PROGRAM_NAME "NoiseComparison"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

/***************************************************************************/

/* STRUCTURES */

struct CommandLineArgsTag {
  char *freqfile;          /* File with frequencies */
  REAL8 fcal;              /* Calibration line frequency */
  REAL8 b;                 /* band */
  REAL8 t;                 /* Time interval to compute the FFT */
  char *hoftFrCacheFile;   /* Frame cache file for h(t) */
  char *derrFrCacheFile;   /* Frame cache file for DARM_ERR */
  char *derr_chan;         /* DARM_ERR channel name */
  char *exc_chan;          /* DARM_ERR channel name */
  char *darm_chan;         /* DARM_ERR channel name */
  char *asq_chan;          /* DARM_ERR channel name */
  char *hoft_chan;         /* h(t) channel name */
  char *noisefile;         /* output file for the noise */
  char *noisefilebin;      /* output file for the individual bins */
  char *OLGFile;           /* open loop gain file */
  char *SensingFile;       /* sensing function file */
  INT4 GPSStart;           /* Start and end GPS times of the segment to be compared*/
  INT4 GPSEnd;
  REAL8 G0Re;              /* Real part of open loop gain at cal line freq.*/
  REAL8 G0Im;              /* Imaginary part of open loop gain at cal line freq. */
  REAL8 D0Re;              /* Real part of digital filter at cal line freq.*/
  REAL8 D0Im;              /* Imaginary part of digital filter at cal line freq. */
  REAL8 W0Re;              /* Real part of whitening filter at cal line freq.*/
  REAL8 W0Im;              /* Imaginary part of whitening filter at cal line freq. */
  REAL8 gamma_fudgefactor; /* fudge factor to divide gammas by */
  INT4 nofactors;
  INT4 outputphase;
  INT4 decimate;
} CommandLineArgs;

typedef struct ResponseFunctionTag
{
  REAL4 Frequency[MAXLINESRS];
  REAL4 Magnitude[MAXLINESRS];
  REAL4 Phase[MAXLINESRS];
  REAL4 re[MAXLINESRS];
  REAL4 im[MAXLINESRS];
} Response;

/***************************************************************************/

/* GLOBAL VARIABLES */

static LALStatus status;
INT4 lalDebugLevel=0;

FrCache *hoftframecache=NULL;                                           /* frame reading variables */
FrStream *hoftframestream=NULL;

FrCache *derrframecache=NULL;                                           /* frame reading variables */
FrStream *derrframestream=NULL;

LIGOTimeGPS gpsepoch;

INT4 duration;

static REAL4TimeSeries derr;
static REAL8TimeSeries hoft;
static FrChanIn chanin_derr;
static FrChanIn chanin_hoft;
REAL4Window *derrwin=NULL;

FILE *fpout=NULL;
FILE *fpoutbin=NULL;
COMPLEX16Vector *ffthtData = NULL;
COMPLEX8Vector *fftderrData = NULL;
REAL8FFTPlan *fftPlanDouble=NULL;
REAL4FFTPlan *fftPlan=NULL;
FrPos derrpos;
FrPos hoftpos;

Response OLG0, OLG[MAXFREQUENCIES], Sensing0, Sensing[MAXFREQUENCIES];

REAL8 gamma_fac[MAXFACTORS];

REAL8 frequencies[MAXFREQUENCIES]; /* frequency array */
INT4 Nfrequencies;  /* number of frequencies */

#if 0
static void printmemuse(void)
{
   pid_t mypid=getpid();
   char commandline[256];
   fflush(NULL);
   sprintf(commandline,"cat /proc/%d/status | /bin/grep Vm | /usr/bin/fmt -140 -u", (int)mypid);
   system(commandline);
   fflush(NULL);
}
#endif

/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
/* */
int Initialise(struct CommandLineArgsTag CLA);
/* */
int GetFactors(struct CommandLineArgsTag CLA);
/* */
int ReadCalFiles(struct CommandLineArgsTag CLA);
/* */
int ComputeNoise(struct CommandLineArgsTag CLA, int i);
/* */
int Finalise(void);

/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
  int i, N;

   #ifdef DEBUG
     fprintf(stdout,"Made it to: -4\n");
     printmemuse();
   #endif

  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;

   #ifdef DEBUG
     fprintf(stdout,"Made it to: -3\n");
     printmemuse();
   #endif

  if (ReadCalFiles(CommandLineArgs)) return 2;

  duration = CommandLineArgs.GPSEnd-CommandLineArgs.GPSStart;

  N = duration / CommandLineArgs.t ;

  if (N > MAXFACTORS-1 )
    {
      fprintf(stderr,"Too many subsegments in segment, exiting.\n");
      return 10;
    }

  gpsepoch.gpsSeconds=CommandLineArgs.GPSStart;
  gpsepoch.gpsNanoSeconds=0;

  if(GetFactors(CommandLineArgs)) return 4;

  if(Initialise(CommandLineArgs)) return 3;

  for(i=0;i<N;i++)
    {
      gpsepoch.gpsSeconds=CommandLineArgs.GPSStart+i*CommandLineArgs.t;
      gpsepoch.gpsNanoSeconds=0;

      if(ComputeNoise(CommandLineArgs,i)) return 5;
    }

  if(Finalise()) return 6;

  return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/


/*  FUNCTIONS */
/*******************************************************************************/

int Initialise(struct CommandLineArgsTag CLA)
{
  /* create Frame cache, open frame stream and delete frame cache */
  LALFrCacheImport(&status,&derrframecache,CommandLineArgs.derrFrCacheFile);
  TESTSTATUS( &status );
  LALFrCacheOpen(&status,&derrframestream,derrframecache);
  TESTSTATUS( &status );
  LALDestroyFrCache(&status,&derrframecache);
  TESTSTATUS( &status );

  LALFrCacheImport(&status,&hoftframecache,CommandLineArgs.hoftFrCacheFile);
  TESTSTATUS( &status );
  LALFrCacheOpen(&status,&hoftframestream,hoftframecache);
  TESTSTATUS( &status );
  LALDestroyFrCache(&status,&hoftframecache);
  TESTSTATUS( &status );

  chanin_derr.type  = ADCDataChannel;
  chanin_hoft.type = ProcDataChannel;

  chanin_hoft.name  = CLA.hoft_chan;
  chanin_derr.name = CLA.derr_chan;

  /* Get channel time step size by calling LALFrGetREAL4TimeSeries */
  /* DARM_ERR */
  LALFrSeek(&status,&gpsepoch,derrframestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&derrpos,derrframestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&derr,&chanin_derr,derrframestream);
  TESTSTATUS( &status );

  /* h(t) */
  LALFrSeek(&status,&gpsepoch,hoftframestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&hoftpos,hoftframestream);
  TESTSTATUS( &status );
  LALFrGetREAL8TimeSeries(&status,&hoft,&chanin_hoft,hoftframestream);
  TESTSTATUS( &status );

  /* Allocate space for data vectors */
  LALCreateVector(&status,&derr.data,(UINT4)(CLA.t/(derr.deltaT) +0.5));
  TESTSTATUS( &status );

  LALDCreateVector(&status,&hoft.data,(UINT4)(CLA.t/(hoft.deltaT) +0.5));
  TESTSTATUS( &status );


  if (CLA.decimate == 1)
    {
     #ifdef DEBUG
      fprintf(stdout,"(UINT4) CLA.t/(derr.deltaT) + 0.5 is %i\n",(UINT4)(CLA.t/(derr.deltaT*DECIMATE) +0.5));
       printmemuse();
     #endif

    /* make window  */
    derrwin = XLALCreateHannREAL4Window((INT4)(CLA.t/(derr.deltaT*DECIMATE) +0.5));
    }
  else
    {
     #ifdef DEBUG
      fprintf(stdout,"(UINT4) CLA.t/(derr.deltaT) + 0.5 is %i\n",(UINT4)(CLA.t/(derr.deltaT) +0.5));
       printmemuse();
     #endif

    /* make window  */
    #ifdef DEBUG
      fprintf(stdout,"(INT4) CLA.t/(derr.deltaT) + 0.5 is %i\n",(INT4)(CLA.t/(derr.deltaT) +0.5));
       printmemuse();
     #endif

    derrwin = XLALCreateHannREAL4Window((INT4)(CLA.t/(derr.deltaT) +0.5));

     #ifdef DEBUG
      fprintf(stdout,"Made it to -1b\n");
       printmemuse();
     #endif

    }

  /* Open output file */
  fpout=fopen(CLA.noisefile,"w");
  if (fpout==NULL)
    {
      fprintf(stderr,"Could not open %s!\n",CLA.noisefile);
      return 1;
    }
  /*   setvbuf( fpout, NULL, _IONBF, 0 );  */
   #ifdef DEBUG
     fprintf(stdout,"Made it to: 0\n");
     printmemuse();
   #endif
   /* Open individual bin output file */
   fpoutbin=fopen(CLA.noisefilebin,"w");
   if (fpoutbin==NULL)
     {
       fprintf(stderr,"Could not open %s!\n",CLA.noisefilebin);
       return 1;
     }
   #ifdef DEBUG
     fprintf(stdout,"Made it to: 0b \n");
     printmemuse();
   #endif

  if (CLA.decimate == 1)
    {
    LALCreateForwardREAL8FFTPlan( &status, &fftPlanDouble, hoft.data->length/DECIMATE, 0 );
    TESTSTATUS( &status );

     #ifdef DEBUG
       fprintf(stdout,"Made it to: 1\n");
       printmemuse();
     #endif
    LALCreateForwardREAL4FFTPlan( &status, &fftPlan, derr.data->length/DECIMATE, 0 );
    TESTSTATUS( &status );
     #ifdef DEBUG
       fprintf(stdout,"Made it to: 2\n");
       printmemuse();
     #endif
    LALZCreateVector( &status, &ffthtData, (hoft.data->length/DECIMATE) / 2 + 1 );
    TESTSTATUS( &status );

     #ifdef DEBUG
       fprintf(stdout,"Made it to: 3\n");
       printmemuse();
     #endif
    LALCCreateVector( &status, &fftderrData, (hoft.data->length/DECIMATE) / 2 + 1 );
    TESTSTATUS( &status );
     #ifdef DEBUG
       fprintf(stdout,"Made it to: 4\n");
       printmemuse();
     #endif
     }
  else
    {
    LALCreateForwardREAL8FFTPlan( &status, &fftPlanDouble, hoft.data->length, 0 );
    TESTSTATUS( &status );

     #ifdef DEBUG
       fprintf(stdout,"Made it to: 1\n");
       printmemuse();
     #endif
    LALCreateForwardREAL4FFTPlan( &status, &fftPlan, derr.data->length, 0 );
    TESTSTATUS( &status );
     #ifdef DEBUG
       fprintf(stdout,"Made it to: 2\n");
       printmemuse();
     #endif
    LALZCreateVector( &status, &ffthtData, (hoft.data->length) / 2 + 1 );
    TESTSTATUS( &status );

     #ifdef DEBUG
       fprintf(stdout,"Made it to: 3\n");
       printmemuse();
     #endif
    LALCCreateVector( &status, &fftderrData, (hoft.data->length) / 2 + 1 );
    TESTSTATUS( &status );
     #ifdef DEBUG
       fprintf(stdout,"Made it to: 4\n");
       printmemuse();
     #endif
     }
  return 0;
}

/*******************************************************************************/

int GetFactors(struct CommandLineArgsTag CLA)
{

FrPos pos1;

static REAL4TimeSeries darm;
static REAL4TimeSeries exc;

static FrChanIn chanin_darm;
static FrChanIn chanin_exc;

CalFactors factors;
UpdateFactorsParams params;

REAL4Window *excwin=NULL,*darmwin=NULL;  /* window */

INT4 k,m;
LIGOTimeGPS localgpsepoch=gpsepoch; /* Local variable epoch used to calculate the calibration factors */
long double gtime=(long double)(localgpsepoch.gpsSeconds+(long double)localgpsepoch.gpsNanoSeconds*1E-9);

FrCache *framecache=NULL;                                           /* frame reading variables */
FrStream *framestream=NULL;

  /* create Frame cache, open frame stream and delete frame cache */
  LALFrCacheImport(&status,&framecache,CommandLineArgs.derrFrCacheFile);
  TESTSTATUS( &status );
  LALFrCacheOpen(&status,&framestream,framecache);
  TESTSTATUS( &status );
  LALDestroyFrCache(&status,&framecache);
  TESTSTATUS( &status );

  chanin_darm.type = ADCDataChannel;
  chanin_exc.type  = ADCDataChannel;

  chanin_darm.name = CLA.darm_chan;
  chanin_exc.name  = CLA.exc_chan;

  /* Get channel time step size by calling LALFrGetREAL4TimeSeries */
  LALFrSeek(&status,&localgpsepoch,framestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );

  LALFrGetREAL4TimeSeries(&status,&darm,&chanin_darm,framestream);
  TESTSTATUS( &status );
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&exc,&chanin_exc,framestream);
  TESTSTATUS( &status );

  /* Allocate space for data vectors */
  LALCreateVector(&status,&darm.data,(UINT4)(CLA.t/darm.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&exc.data,(UINT4)(CLA.t/exc.deltaT +0.5));
  TESTSTATUS( &status );

  /* windows for time domain channels */
  /* darm */
  darmwin = XLALCreateHannREAL4Window((INT4)(CLA.t/darm.deltaT +0.5));

  /* exc */
  excwin = XLALCreateHannREAL4Window((INT4)(CLA.t/exc.deltaT +0.5));

  for(m=0;m < (INT4)(duration/CLA.t);m++)
    {
      /* Fill data vectors with data */
      LALFrSeek(&status,&localgpsepoch,framestream);
      TESTSTATUS( &status );
      LALFrGetPos(&status,&pos1,framestream);
      TESTSTATUS( &status );

      LALFrSetPos(&status,&pos1,framestream);
      TESTSTATUS( &status );
      LALFrGetREAL4TimeSeries(&status,&darm,&chanin_darm,framestream);
      TESTSTATUS( &status );

      LALFrSetPos(&status,&pos1,framestream);
      TESTSTATUS( &status );
      LALFrGetREAL4TimeSeries(&status,&exc,&chanin_exc,framestream);
      TESTSTATUS( &status );

      /* Window the data */
      for(k=0;k<(INT4)(CLA.t/darm.deltaT +0.5);k++)
	{
	  darm.data->data[k] *= 2.0*darmwin->data->data[k];
	}
      for(k=0;k<(INT4)(CLA.t/exc.deltaT +0.5);k++)
	{
	  exc.data->data[k] *= 2.0*excwin->data->data[k];
	}

      /* set params to call LALComputeCalibrationFactors */
      params.darmCtrl = &darm;
      params.asQ = &darm;
      params.exc = &exc;
      params.lineFrequency = CLA.fcal;
      params.openloop = crect( CLA.G0Re, CLA.G0Im );
      params.digital = crect( CLA.D0Re, CLA.D0Im );
      params.whitener = crect( CLA.W0Re, CLA.W0Im );

      LALComputeCalibrationFactors(&status,&factors,&params);
      TESTSTATUS( &status );

      /* put factors into series */
      factors.alphabeta = crect( creal(factors.alphabeta) / CLA.gamma_fudgefactor, cimag(factors.alphabeta) );
      gamma_fac[m]   = creal(factors.alphabeta);

      if(CLA.nofactors)
	{
	  gamma_fac[m]   = 1.0;
	}

      fprintf(stdout,"%18.9Lf %e %e %e %e %e %e %e %e %e %e %e %e %e \n",gtime,
	      creal(factors.alpha),cimag(factors.alpha),
	      creal(factors.beta),cimag(factors.beta),
	      creal(factors.alphabeta),gamma_fac[m],cimag(factors.alphabeta),
	      creal(factors.asq)*2/CLA.t,cimag(factors.asq)*2/CLA.t,
	      creal(factors.darm)*2/CLA.t,cimag(factors.darm)*2/CLA.t,
	      creal(factors.exc)*2/CLA.t,cimag(factors.exc)*2/CLA.t);

      gtime += CLA.t;
      localgpsepoch.gpsSeconds = (INT4)gtime;
      localgpsepoch.gpsNanoSeconds = (INT4)((gtime-(INT4)gtime)*1E+09);
    }

  LALDestroyVector(&status,&darm.data);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&exc.data);
  TESTSTATUS( &status );

  XLALDestroyREAL4Window(darmwin);
  XLALDestroyREAL4Window(excwin);

  LALFrClose(&status,&framestream);
  TESTSTATUS( &status );

  return 0;
}

/*******************************************************************************/

int ComputeNoise(struct CommandLineArgsTag CLA, int n)
{

  INT4 k,j;
  REAL8 mean_Sh_derr;
  REAL8 mean_Sh_hoft;
  COMPLEX8 R,H,C;

  /* Fill data vectors with data */
  /* DARM_ERR */
  LALFrSeek(&status,&gpsepoch,derrframestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&derrpos,derrframestream);
  TESTSTATUS( &status );

  LALFrSetPos(&status,&derrpos,derrframestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&derr,&chanin_derr,derrframestream);
  TESTSTATUS( &status );

  /* h(t) */
  LALFrSeek(&status,&gpsepoch,hoftframestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&hoftpos,hoftframestream);
  TESTSTATUS( &status );

  LALFrSetPos(&status,&hoftpos,hoftframestream);
  TESTSTATUS( &status );
  LALFrGetREAL8TimeSeries(&status,&hoft,&chanin_hoft,hoftframestream);
  TESTSTATUS( &status );

  if ( derr.epoch.gpsSeconds !=  gpsepoch.gpsSeconds || hoft.epoch.gpsSeconds !=  gpsepoch.gpsSeconds)
    {
      fprintf(stderr,"GPS start time of darm err data (%d) or h(t) data (%d) does not agree with requested start time (%d). Exiting.\n",
	      derr.epoch.gpsSeconds, hoft.epoch.gpsSeconds, gpsepoch.gpsSeconds);
      return 1;
    }

    #ifdef DEBUG
     fprintf(stdout,"Made it to: 5\n");
     printmemuse();
    #endif

  if (CLA.decimate == 1)
    {
      XLALResampleREAL4TimeSeries( &derr, derr.deltaT*DECIMATE);
      TESTSTATUS( &status );

     #ifdef DEBUG
       fprintf(stdout,"Made it to: 6\n");
       printmemuse();
     #endif

      XLALResampleREAL8TimeSeries( &hoft, hoft.deltaT*DECIMATE);
      TESTSTATUS( &status );
    }

   #ifdef DEBUG
     fprintf(stdout,"Made it to: 7\n");
     printmemuse();
   #endif

  /* Window the data */
  for(k=0;k<(INT4)(CLA.t/derr.deltaT +0.5);k++)
    {
      derr.data->data[k] *= 2.0*derrwin->data->data[k];
    }
  for(k=0;k<(INT4)(CLA.t/hoft.deltaT +0.5);k++)
    {
      hoft.data->data[k] *= 2.0*derrwin->data->data[k];
    }

  #ifdef DEBUG
  for(k=0;k<(INT4)(CLA.t/hoft.deltaT +0.5);k++)
    {
      fprintf(stdout,"%e\n",hoft.data->data[k]);
    }
  #endif

  /* FFT the data */
  XLALREAL8ForwardFFT( ffthtData, hoft.data, fftPlanDouble );

   #ifdef DEBUG
     fprintf(stdout,"Made it to: 8\n");
     printmemuse();
   #endif
  XLALREAL4ForwardFFT( fftderrData, derr.data, fftPlan );

   #ifdef DEBUG
     fprintf(stdout,"Made it to: 9\n");
     printmemuse();
   #endif

  fprintf(fpout, "%d ", gpsepoch.gpsSeconds);
  fprintf(fpoutbin, "%d \n", gpsepoch.gpsSeconds);

  for (j=0; j < Nfrequencies; j++)
    {
      /* Compute the noise */
      mean_Sh_derr = 0.0;
      mean_Sh_hoft = 0.0;
      {
	int firstbin=(INT4)(frequencies[j]*CLA.t+0.5);
	int nsamples=(INT4)(CLA.b*CLA.t+0.5);

	for (k=0; k<nsamples; k++)
	  {
	    REAL4 dr= crealf(fftderrData->data[k+firstbin]);
	    REAL4 di= cimagf(fftderrData->data[k+firstbin]);
	    REAL8 dpsq;
	    REAL4 caldr,caldi;

	    REAL4 hr= creal(ffthtData->data[k+firstbin]);
	    REAL4 hi= cimag(ffthtData->data[k+firstbin]);
	    REAL8 hpsq;

	    /* Compute Response */
	    C = crectf( gamma_fac[n]*Sensing[j].re[k], gamma_fac[n]*Sensing[j].im[k] );

	    H = crectf( 1.0+gamma_fac[n]*OLG[j].re[k], gamma_fac[n]*OLG[j].im[k] );
	    R = H / C;

	    /* Apply response to DARM_ERR */
	    caldr=dr*crealf(R)-di*cimagf(R);
	    caldi=di*crealf(R)+dr*cimagf(R);

	    dpsq=pow(caldr,2.0)+pow(caldi,2);
	    hpsq=pow(hr,2.0)+pow(hi,2);

	    mean_Sh_derr += dpsq/nsamples;
	    mean_Sh_hoft += hpsq/nsamples;

	    if(k==0)
	      {
 	        fprintf(fpoutbin,"%e %e %e %e \n",hr*sqrt(2.0*hoft.deltaT/(REAL4)hoft.data->length),hi*sqrt(2.0*hoft.deltaT/(REAL4)hoft.data->length),caldr*sqrt(2.0*derr.deltaT/(REAL4)derr.data->length),caldi*sqrt(2.0*derr.deltaT/(REAL4)derr.data->length));
              }
	  }
	mean_Sh_derr *= 2.0*derr.deltaT/(REAL4)derr.data->length;
	mean_Sh_hoft *= 2.0*hoft.deltaT/(REAL4)hoft.data->length;


      }
      fprintf(fpout, "%e %e ",sqrt(mean_Sh_hoft), sqrt(mean_Sh_derr));

    }

  fprintf(fpout, "\n");

  if (CLA.decimate == 1)
    {
    /* resize the time series back to original length */
    XLALResizeREAL8TimeSeries(&hoft, 0,hoft.data->length*DECIMATE);
    hoft.deltaT= hoft.deltaT/DECIMATE;

    XLALResizeREAL4TimeSeries(&derr, 0,derr.data->length*DECIMATE);
    derr.deltaT= derr.deltaT/DECIMATE;

     #ifdef DEBUG
       fprintf(stdout,"Made it to: 10\n");
       printmemuse();
     #endif
     }

  return 0;
}


/*******************************************************************************/

int ReadCalFiles(struct CommandLineArgsTag CLA)
{
  char line[256];

  INT4 i,j,k;
  FILE *fpR;
  int nsamples=(INT4)(CLA.b*CLA.t+0.5);
  REAL4 f,df;

 /* ------ Open and read Frequency File ------ */
 i=0;
 fpR=fopen(CLA.freqfile,"r");
 if (fpR==NULL)
   {
     fprintf(stderr,"Weird... %s doesn't exist!\n",CLA.freqfile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpR))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i >= MAXFREQUENCIES)
       {
	 fprintf(stderr,"Too many lines in file %s! Exiting... \n", CLA.freqfile);
	 return 1;
       }
     sscanf(line,"%le",&frequencies[i]);
     i++;
   }
 Nfrequencies=i;
 /* -- close OLG file -- */
 fclose(fpR);


 /* ------ Open and read OLG File ------ */
 i=0;
 fpR=fopen(CLA.OLGFile,"r");
 if (fpR==NULL)
   {
     fprintf(stderr,"Weird... %s doesn't exist!\n",CLA.OLGFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpR))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i >= MAXLINESRS)
       {
	 fprintf(stderr,"Too many lines in file %s! Exiting... \n", CLA.OLGFile);
	 return 1;
       }
     sscanf(line,"%e %e %e",&OLG0.Frequency[i],&OLG0.Magnitude[i],&OLG0.Phase[i]);
     OLG0.re[i]=OLG0.Magnitude[i]*cos(OLG0.Phase[i]);
     OLG0.im[i]=OLG0.Magnitude[i]*sin(OLG0.Phase[i]);
     i++;
   }
 /* -- close OLG file -- */
 fclose(fpR);


 /* for each frequency in the ferquency file: */
 for (k=0; k < Nfrequencies; k++)
   {
     /* compute a response appropriate for the frequency spacing of our data */
     i=0;
     for (j=0; j < nsamples;j++)
       {
	 REAL4 a,b;

	 f = frequencies[k]  + j/CLA.t;

	 while (i < MAXLINESRS-1 && f > OLG0.Frequency[i]) i++;

	 if(i == MAXLINESRS-1 && f > OLG0.Frequency[i])
	   {
	     fprintf(stderr,"No calibration info for frequency %f!\n",f);
	     return 1;
	   }
	 /* checks order */

	 /* Since now Rraw.Frequency[i-1] < f =< Rraw.Frequency[i] ||*/
	 /* can compute a weighted average of the raw frequencies at i-1 and i */

	 /* check both bounds!!!!! */
	 if(f < OLG0.Frequency[i-1] || f > OLG0.Frequency[i])
	   {
	     fprintf(stderr,"Frequency %f in SFT does not lie between two lines in OLG file!\n",f);
	     return 1;
	   }

	 /* If the frequencies are closely spaced this may create dangerous floating point errors */
	 df=OLG0.Frequency[i]-OLG0.Frequency[i-1];

	 a=(f-OLG0.Frequency[i-1])/df;
	 if (a>1.0) a=1.0;
	 b=1.0-a;

	 OLG[k].Frequency[j]=f;
	 OLG[k].Magnitude[j]=a*OLG0.Magnitude[i]+b*OLG0.Magnitude[i-1];
	 OLG[k].Phase[j]=a*OLG0.Phase[i]+b*OLG0.Phase[i-1];

	 OLG[k].re[j]=a*OLG0.re[i]+b*OLG0.re[i-1];
	 OLG[k].im[j]=a*OLG0.im[i]+b*OLG0.im[i-1];
       }
   }

 /* ------ Open and read Sensing File ------ */
 i=0;
 fpR=fopen(CLA.SensingFile,"r");
 if (fpR==NULL)
   {
     fprintf(stderr,"Weird... %s doesn't exist!\n",CLA.SensingFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpR))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i >= MAXLINESRS)
       {
	 fprintf(stderr,"Too many lines in file %s! Exiting... \n", CLA.SensingFile);
	 return 1;
       }
     sscanf(line,"%e %e %e",&Sensing0.Frequency[i],&Sensing0.Magnitude[i],&Sensing0.Phase[i]);
     Sensing0.re[i]=Sensing0.Magnitude[i]*cos(Sensing0.Phase[i]);
     Sensing0.im[i]=Sensing0.Magnitude[i]*sin(Sensing0.Phase[i]);
     i++;
   }
 /* -- close Sensing file -- */
 fclose(fpR);


 /* for each frequency in the ferquency file: */
 for (k=0; k < Nfrequencies; k++)
   {
     /* compute a response appropriate for the frequency spacing of our data */
     i=0;
     for (j=0; j < nsamples;j++)
       {
	 REAL4 a,b;

	 f = frequencies[k]  + j/CLA.t;

	 while (i < MAXLINESRS-1 && f > Sensing0.Frequency[i]) i++;

	 if(i == MAXLINESRS-1 && f > Sensing0.Frequency[i])
	   {
	     fprintf(stderr,"No calibration info for frequency %f!\n",f);
	     return 1;
	   }
	 /* checks order */

	 /* Since now Rraw.Frequency[i-1] < f =< Rraw.Frequency[i] ||*/
	 /* can compute a weighted average of the raw frequencies at i-1 and i */

	 /* check both bounds!!!!! */
	 if(f < Sensing0.Frequency[i-1] || f > Sensing0.Frequency[i])
	   {
	     fprintf(stderr,"Frequency %f in SFT does not lie between two lines in Sensing file!\n",f);
	     return 1;
	   }

	 /* If the frequencies are closely spaced this may create dangerous floating point errors */
	 df=Sensing0.Frequency[i]-Sensing0.Frequency[i-1];

	 a=(f-Sensing0.Frequency[i-1])/df;
	 if (a>1.0) a=1.0;
	 b=1.0-a;

	 Sensing[k].Frequency[j]=f;
	 Sensing[k].Magnitude[j]=a*Sensing0.Magnitude[i]+b*Sensing0.Magnitude[i-1];
	 Sensing[k].Phase[j]=a*Sensing0.Phase[i]+b*Sensing0.Phase[i-1];

	 Sensing[k].re[j]=a*Sensing0.re[i]+b*Sensing0.re[i-1];
	 Sensing[k].im[j]=a*Sensing0.im[i]+b*Sensing0.im[i-1];
       }
   }

  return 0;
}

/*******************************************************************************/

int Finalise(void)
{

   #ifdef DEBUG
   fprintf(stdout,"Made it to: 12\n");
   printmemuse();
   #endif

  LALDDestroyVector(&status,&hoft.data);
  TESTSTATUS( &status );

   #ifdef DEBUG
   fprintf(stdout,"Made it to: 11\n");
   printmemuse();
   #endif

  LALDestroyVector(&status,&derr.data);
  TESTSTATUS( &status );

   #ifdef DEBUG
   fprintf(stdout,"Made it to: 13\n");
   printmemuse();
   #endif

  XLALDestroyREAL4Window(derrwin);


   #ifdef DEBUG
   fprintf(stdout,"Made it to: 14\n");
   printmemuse();
   #endif

  LALZDestroyVector( &status, &ffthtData);

   #ifdef DEBUG
   fprintf(stdout,"Made it to: 15\n");
   printmemuse();
   #endif

  LALCDestroyVector( &status, &fftderrData);

   #ifdef DEBUG
   fprintf(stdout,"Made it to: 16\n");
   printmemuse();
   #endif

  LALDestroyREAL8FFTPlan( &status, &fftPlanDouble );

   #ifdef DEBUG
   fprintf(stdout,"Made it to: 17\n");
   printmemuse();
   #endif

  LALDestroyREAL4FFTPlan( &status, &fftPlan );

   #ifdef DEBUG
   fprintf(stdout,"Made it to: 18\n");
   printmemuse();
   #endif

  LALFrClose(&status,&derrframestream);
  TESTSTATUS( &status );

   #ifdef DEBUG
   fprintf(stdout,"Made it to: 19\n");
   printmemuse();
   #endif

  LALFrClose(&status,&hoftframestream);
  TESTSTATUS( &status );

   #ifdef DEBUG
   fprintf(stdout,"Made it to: 20\n");
   printmemuse();
   #endif

  fclose(fpout);
  fclose(fpoutbin);

  LALCheckMemoryLeaks();

  return 0;
}

/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA)
{

  INT4 errflg=0;
  struct option long_options[] = {
    {"freq-file",            required_argument, NULL,           'a'},
    {"time",                 required_argument, NULL,           'b'},
    {"band",                 required_argument, NULL,           'c'},
    {"hoft-cache",           required_argument, NULL,           'd'},
    {"derr-cache",           required_argument, NULL,           'e'},
    {"derr-channel",         required_argument, NULL,           'f'},
    {"darm-channel",         required_argument, NULL,           'g'},
    {"exc-channel",          required_argument, NULL,           'i'},
    {"asq-channel",          required_argument, NULL,           'j'},
    {"hoft-channel",         required_argument, NULL,           'k'},
    {"gps-start-time",       required_argument, NULL,           'l'},
    {"gps-end-time",         required_argument, NULL,           'm'},
    {"output-file",          required_argument, NULL,           'n'},
    {"output-bin",           required_argument, NULL,           'B'},
    {"olg-file",             required_argument, NULL,           'o'},
    {"sensing-file",         required_argument, NULL,           'p'},
    {"olg-re",               required_argument, NULL,           'q'},
    {"olg-im",               required_argument, NULL,           'r'},
    {"servo-re",             required_argument, NULL,           's'},
    {"servo-im",             required_argument, NULL,           't'},
    {"whitener-re",          required_argument, NULL,           'u'},
    {"whitener-im",          required_argument, NULL,           'v'},
    {"fcal",                 required_argument, NULL,           'w'},
    {"version",              no_argument, NULL,                 'x'},
    {"gamma-fudge-factor",   required_argument, NULL,           'y'},
    {"no-factors",           no_argument, NULL,                 'z'},
    {"output-phase",         no_argument, NULL,                 'Q'},
    {"decimate",             required_argument, NULL,           'D'},
    {"help",                 no_argument, NULL,                 'h'},
    {0, 0, 0, 0}
  };
  char args[] = "ha:b:c:d:e:g:i:j:k:l:m:n:B:o:p:q:r:s:t:u:v:w:xy:zQD";

  /* Initialize default values */
  CLA->freqfile=NULL;
  CLA->fcal=0.0;
  CLA->b=0.0;
  CLA->t=0.0;
  CLA->hoftFrCacheFile=NULL;
  CLA->derrFrCacheFile=NULL;
  CLA->derr_chan=NULL;
  CLA->asq_chan=NULL;
  CLA->darm_chan=NULL;
  CLA->exc_chan=NULL;
  CLA->hoft_chan=NULL;
  CLA->noisefile=NULL;
  CLA->noisefilebin=NULL;
  CLA->OLGFile=NULL;
  CLA->SensingFile=NULL;
  CLA->GPSStart = 0;
  CLA->GPSEnd = 0;
  CLA->G0Re=0.0;
  CLA->G0Im=0.0;
  CLA->D0Re=0.0;
  CLA->D0Im=0.0;
  CLA->W0Re=0.0;
  CLA->W0Im=0.0;
  CLA->gamma_fudgefactor=1.0;
  CLA->outputphase=0;
  CLA->nofactors=0;
  CLA->decimate=0;

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
    case 'x':
      fprintf(stdout,"Program: %s\n",PROGRAM_NAME);
      fprintf(stdout,"Revision: %s\n", CVS_REVISION );
      fprintf(stdout,"Source: %s\n", CVS_SOURCE);
      fprintf(stdout,"Date: %s\n", CVS_DATE);
      exit(0);
      break;
    case 'a':
      CLA->freqfile=optarg;
      break;
    case 'b':
      CLA->t=atof(optarg);
      break;
    case 'c':
      CLA->b=atof(optarg);
      break;
    case 'd':
      CLA->hoftFrCacheFile=optarg;
      break;
    case 'e':
      CLA->derrFrCacheFile=optarg;
      break;
    case 'f':
      /* name of darm channel */
      CLA->derr_chan=optarg;
      break;
    case 'g':
      /* name of darm channel */
      CLA->darm_chan=optarg;
      break;
    case 'i':
      /* name of darm channel */
      CLA->exc_chan=optarg;
      break;
    case 'j':
      /* name of darm channel */
      CLA->asq_chan=optarg;
      break;
    case 'k':
      /* name of darm err channel */
      CLA->hoft_chan=optarg;
      break;
    case 'l':
      /* GPS start */
      CLA->GPSStart=atoi(optarg);
      break;
    case 'm':
      /* GPS end */
      CLA->GPSEnd=atof(optarg);
      break;
    case 'n':
      /* real part of OLG */
      CLA->noisefile=optarg;
      break;
    case 'B':
      /* real part of OLG */
      CLA->noisefilebin=optarg;
      break;
    case 'o':
      /* real part of OLG */
      CLA->OLGFile=optarg;
      break;
    case 'p':
      /* real part of OLG */
      CLA->SensingFile=optarg;
      break;
    case 'q':
      /* real part of OLG */
      CLA->G0Re=atof(optarg);
      break;
    case 'r':
      /* imaginary part of OLG */
      CLA->G0Im=atof(optarg);
      break;
    case 's':
      /*  real part of servo*/
      CLA->D0Re=atof(optarg);
      break;
    case 't':
      /*  imaginary part of servo */
      CLA->D0Im=atof(optarg);
      break;
    case 'u':
      /*  real part of servo*/
      CLA->W0Re=atof(optarg);
      break;
    case 'v':
      /*  imaginary part of servo */
      CLA->W0Im=atof(optarg);
      break;
    case 'w':
      /*  imaginary part of servo */
      CLA->fcal=atof(optarg);
      break;
    case 'y':
      CLA->gamma_fudgefactor=atof(optarg);
      break;
    case 'z':
      CLA->nofactors=1;
      break;
    case 'Q':
      CLA->outputphase=1;
      break;
    case 'D':
      CLA->decimate=1;
      break;
    case 'h':
      /* print usage/help message */
      fprintf(stdout,"Arguments are:\n");
      fprintf(stdout,"\t freq-file (-a)\tFLOAT\t Name of file with frequencies.\n");
      fprintf(stdout,"\t time (-b)\tFLOAT\t Integration time (s).\n");
      fprintf(stdout,"\t band (-c)\tFLOAT\t Frequency band (Hz).\n");
      fprintf(stdout,"\t hoft-cache (-d)\tFLOAT\t h(t) cache file name.\n");
      fprintf(stdout,"\t derr-cache (-e)\tFLOAT\t DARM_ERR cache file name.\n");
      fprintf(stdout,"\t derr-channel (-f)\tFLOAT\t DARM_ERR channel name.\n");
      fprintf(stdout,"\t darm-channel (-g)\tFLOAT\t DARM_CTRL channel name.\n");
      fprintf(stdout,"\t exc-channel (-i)\tFLOAT\t EXC channel name.\n");
      fprintf(stdout,"\t asq-channel (-j)\tFLOAT\t ASQ channel name.\n");
      fprintf(stdout,"\t hoft-channel (-k)\tFLOAT\t h(t) channel name.\n");
      fprintf(stdout,"\t gps-start-time (-l)\tINT\t GPS start time.\n");
      fprintf(stdout,"\t gps-end-time (-m)\tINT\t GPS end time.\n");
      fprintf(stdout,"\t output-file (-n)\tSTRING\t Name of output file.\n");
      fprintf(stdout,"\t output-bin (-B)\tSTRING\t Name of individual bin output file.\n");
      fprintf(stdout,"\t olg-file (-o)\tSTRING\t Name of open loop gain file.\n");
      fprintf(stdout,"\t sensing-file (-p)\tSTRING\t Name of sensing function file.\n");
      fprintf(stdout,"\tolg-re (-q)\tFLOAT\t Real part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\tolg-im (-r)\tFLOAT\t Imaginary part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\tservo-re (-s)\tFLOAT\t Real part of the digital filter at the calibration line frequency.\n");
      fprintf(stdout,"\tservo-im (-t)\tFLOAT\t Imaginary part of digital filter at the calibration line frequency.\n");
      fprintf(stdout,"\twhitener-re (-u)\tFLOAT\t Real part of the whitening filter at the calibration line frequency.\n");
      fprintf(stdout,"\twhitener-im (-v)\tFLOAT\t Imaginary part of whitening filter at the calibration line frequency.\n");
      fprintf(stdout,"\tfcal (-w)\tFLOAT\t Calibration line frequency.\n");
      fprintf(stdout,"\tgamma-fudge-factor (-y)\tFLAG\t Fudge factor used to adjust factor values. Gamma is divided by this value.\n");
      fprintf(stdout,"\tno-factors (-z)\tFLAG\t Set factors to 1.\n");
      fprintf(stdout,"\toutput-phase (-Q)\tFLAG\t Outputs real and imaginary parts of each frequency bin (careful! will make big files!).\n");
      fprintf(stdout,"\tdecimate (-D)\tFLAG\t Decimates the data by a factor of 4.\n");
      fprintf(stdout,"\thelp (-h)\tFLAG\t This message\n");
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

  if(CLA->freqfile == NULL)
    {
      fprintf(stderr,"No frequency file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }

  if(CLA->b == 0.0)
    {
      fprintf(stderr,"No frequency band specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }

  if(CLA->t == 0.0)
    {
      fprintf(stderr,"No integration time specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->GPSStart == 0.0)
    {
      fprintf(stderr,"No GPS start time specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->GPSEnd == 0.0)
    {
      fprintf(stderr,"No GPS end time specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->hoftFrCacheFile==NULL)
    {
      fprintf(stderr,"No h(t) frame cache file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->derrFrCacheFile==NULL)
    {
      fprintf(stderr,"No DARM_ERR frame cache file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->derr_chan==NULL)
    {
      fprintf(stderr,"No DARM_ERR channel name specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->exc_chan==NULL)
    {
      fprintf(stderr,"No EXC channel name specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->asq_chan==NULL)
    {
      fprintf(stderr,"No ASQ channel name specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->darm_chan==NULL)
    {
      fprintf(stderr,"No DARM_CTRL channel name specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->hoft_chan==NULL)
    {
      fprintf(stderr,"No h(t) channel name specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->noisefile==NULL)
    {
      fprintf(stderr,"No output file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->noisefilebin==NULL)
    {
      fprintf(stderr,"No individual bin output file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }

  if(CLA->OLGFile==NULL)
    {
      fprintf(stderr,"No response file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->SensingFile==NULL)
    {
      fprintf(stderr,"No response file specified.\n");
      fprintf(stderr,"Try %s -h \n", argv[0]);
      return 1;
    }
  if(CLA->G0Re == 0.0 )
    {
      fprintf(stderr,"No real part of open loop gain specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }
  if(CLA->G0Im == 0.0 )
    {
      fprintf(stderr,"No imaginary part of open loop gain specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }
  if(CLA->D0Re == 0.0 )
    {
      fprintf(stderr,"No real part of digital filter specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }
  if(CLA->D0Im == 0.0 )
    {
      fprintf(stderr,"No imaginary part of digital filter specified.\n");
      fprintf(stderr,"Try ./ComputeStrainDriver -h \n");
      return 1;
    }


  /* Some sanity checks */
  if (CLA->b < 1/CLA->t)
    {
      fprintf(stderr,"Frequency band requested (%e) is smaller than frequency resolution[1/integration time] (%e)\n", CLA->b, 1.0/CLA->t);
      return 1;
    }
  if (CLA->t < (CLA->GPSStart-CLA->GPSEnd) )
    {
      fprintf(stderr,"Integration time (%e) larger than duration (%d) \n", CLA->t, (CLA->GPSStart-CLA->GPSEnd));
      return 1;
    }

  return errflg;
}

/*******************************************************************************/
#endif
