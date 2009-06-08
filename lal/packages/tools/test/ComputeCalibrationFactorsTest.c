/*
*  Copyright (C) 2007 Bernd Machenschalk, Duncan Brown, Jolien Creighton, Xavier Siemens
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

/**** <lalVerbatim file="ComputeCalibrationFactorsTestCV">
 * Author: X. Siemens
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Program \texttt{ComputeCalibrationFactorsTest.c}}
 *
 * Tests the computation of calibration factors.
 *
 * \subsubsection*{Usage}
 *
 * \begin{verbatim}
 * ComputeCalibrationFactorsTest
 * \end{verbatim}
 *
 * \subsubsection*{Description}
 *
 **** </lalLaTeX> */
#include <lal/LALConfig.h>
#ifndef LAL_FRAME_ENABLED
int main( void ) { return 77; }
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
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/AVFactories.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <math.h>
#include <lal/PrintFTSeries.h>
#include <lal/Window.h>
#include <lal/Calibration.h>

NRCSID (COMPUTECALIBRATIONFACTORSTESTC,"$Id: packages/tools/test/ComputeCalibrationFactorsTest.c $");

#define MAXLINERS 76800                   /*HARDWIRED !!*/
#define MAXLINESEGS 10000                 /*HARDWIRED !!*/

struct CommandLineArgsTag {
  REAL8 f;                 /* Frequency of the calibration line */
  REAL8 k;                 /* Value of output matrix to x arm */
  char *RFile;             /* Text file with the response funtion */
  char *CFile;             /* Text file with the sensing function */
  char *AFile;             /* Text file with the actuation function */
  char *FrCacheFile;       /* Frame cache file */
  char *SegmentsFile;      /* Text file with the segments */
  char *exc_chan;          /* excitation channel name */
  char *darm_chan;         /* darm channel name */
  char *asq_chan;          /* asq channel name */
} CommandLineArgs;

typedef
struct SegmentListTag {
  INT4 tgps;          /* GPS Seconds of start of segment group */
  INT4 nseg;          /* number of segments starting at tgps */
  REAL4 seglength;    /* length of segment in seconds */
} SegmentList;

static LALStatus status;
INT4 lalDebugLevel=3,numsegs;

FrCache *framecache;
FrStream *framestream=NULL;

static REAL4TimeSeries darm;
static REAL4TimeSeries asq;
static REAL4TimeSeries exc;

static FrChanIn chanin_darm;
static FrChanIn chanin_asq;
static FrChanIn chanin_exc;

static COMPLEX16 Rf0,Cf0,Af0;

static CalFactors factors;
static UpdateFactorsParams params;

LIGOTimeGPS epoch;

SegmentList SL[MAXLINESEGS];

LALWindowParams winparams;

REAL4Vector *asqwin=NULL,*excwin=NULL,*darmwin=NULL;

/* Function Prototypes */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);
int ReadFiles(struct CommandLineArgsTag CLA);
int GetChannelNames(struct CommandLineArgsTag CLA);
void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f);
void TimeToFloat(REAL8 *f, LIGOTimeGPS *tgps);
int FreeMem(void);

/**************** MAIN PROGRAM ***************/

int main(int argc,char *argv[])
{
  FrPos pos;
  int i,j,k;

  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;
  if (ReadFiles(CommandLineArgs)) return 3;
  if (GetChannelNames(CommandLineArgs)) return 4;

  for(i=0;i<numsegs;i++)
    {
      for(j=0;j<SL[i].nseg;j++)
	{
	  REAL4 magexc, magdarm, magasq;

	  REAL8 t=SL[i].tgps+j*SL[i].seglength;
	  FloatToTime(&epoch, &t);

	  LALFrSeek(&status,&epoch,framestream);
	  LALFrGetPos(&status,&pos,framestream);

	  LALFrGetREAL4TimeSeries(&status,&exc,&chanin_exc,framestream);
	  LALFrSetPos(&status,&pos,framestream);

	  LALFrGetREAL4TimeSeries(&status,&darm,&chanin_darm,framestream);
	  LALFrSetPos(&status,&pos,framestream);

	  LALFrGetREAL4TimeSeries(&status,&asq,&chanin_asq,framestream);
	  LALFrSetPos(&status,&pos,framestream);

 	  /*Windowing*/
	  for(k=0;k<(INT4)(SL[i].seglength/asq.deltaT +0.5);k++)
	    {
	      asq.data->data[k] *= 2.0*asqwin->data[k];
	    }
	  for(k=0;k<(INT4)(SL[i].seglength/darm.deltaT +0.5);k++)
	    {
	      darm.data->data[k] *= 2.0*darmwin->data[k];
	    }
	  for(k=0;k<(INT4)(SL[i].seglength/exc.deltaT +0.5);k++)
	    {
	      exc.data->data[k] *= 2.0*excwin->data[k];
	    }

	  /* set params to call LALComputeCalibrationFactors */
	  params.darmCtrl = &darm;
	  params.asQ = &asq;
	  params.exc = &exc;
	  params.lineFrequency = CommandLineArgs.f;
	  params.outputMatrix = CommandLineArgs.k;
	  params.actuationFactor.re = Af0.re/4000.0; /*ACHTUNG: HARDWIRED !!*/
	  params.actuationFactor.im = Af0.im/4000.0; /*ACHTUNG: HARDWIRED !!*/
	  params.responseFactor = Rf0;
	  params.sensingFactor = Cf0;

	  LALComputeCalibrationFactors(&status,&factors,&params);

	  magexc=sqrt(factors.exc.re*factors.exc.re+factors.exc.im*factors.exc.im)*2/SL[i].seglength;
	  magasq=sqrt(factors.asq.re*factors.asq.re+factors.asq.im*factors.asq.im)*2/SL[i].seglength;
	  magdarm=sqrt(factors.darm.re*factors.darm.re+factors.darm.im*factors.darm.im)*2/SL[i].seglength;

	  fprintf(stdout,"%20.15f %20.15f %20.15f %20.15f %20.15f %20.15f %20.15f %20.15f\n",t,
		  factors.alpha.re,factors.alpha.im,factors.alphabeta.re,
		  factors.alphabeta.im,magexc,magdarm,magasq);

	  fflush(stdout);

	}
    }

  if(FreeMem()) return 4;

  return 0;
}

/*******************************************************************************/
 void FloatToTime(LIGOTimeGPS *tgps, REAL8 *f)
{
  REAL8 temp0, temp2, temp3;
  REAL8 temp1, temp4;

  temp0 = floor(*f);     /* this is tgps.S */
  temp1 = (*f) * 1.e10;
  temp2 = fmod(temp1, 1.e10);
  temp3 = fmod(temp1, 1.e2);
  temp4 = (temp2-temp3) * 0.1;

  tgps->gpsSeconds = (INT4)temp0;
  tgps->gpsNanoSeconds = (INT4)temp4;
}


/*******************************************************************************/

int GetChannelNames(struct CommandLineArgsTag CLA)
{

  FrPos pos1;

  chanin_asq.type=ADCDataChannel;
  chanin_darm.type=ADCDataChannel;
  chanin_exc.type=ADCDataChannel;

  chanin_asq.name = CLA.asq_chan;
  chanin_darm.name= CLA.darm_chan;
  chanin_exc.name = CLA.exc_chan;

  /* Get channel time step size by calling LALFrGetREAL4TimeSeries */

  epoch.gpsSeconds=SL[0].tgps;
  epoch.gpsNanoSeconds=0;
  LALFrSeek(&status,&epoch,framestream);
  LALFrGetPos(&status,&pos1,framestream);
  LALFrGetREAL4TimeSeries(&status,&asq,&chanin_asq,framestream);
  LALFrSetPos(&status,&pos1,framestream);
  LALFrGetREAL4TimeSeries(&status,&darm,&chanin_darm,framestream);
  LALFrSetPos(&status,&pos1,framestream);
  LALFrGetREAL4TimeSeries(&status,&exc,&chanin_exc,framestream);
  LALFrRewind(&status,framestream);

  /* Determine from the sample rate how many data points to allocate for each segment */

  LALCreateVector(&status,&asq.data,(INT4)(SL[0].seglength/asq.deltaT +0.5));
  LALCreateVector(&status,&darm.data,(INT4)(SL[0].seglength/darm.deltaT +0.5));
  LALCreateVector(&status,&exc.data,(INT4)(SL[0].seglength/exc.deltaT +0.5));

  /* Create Window vectors */
  LALCreateVector(&status,&asqwin,(INT4)(SL[0].seglength/asq.deltaT +0.5));
  LALCreateVector(&status,&darmwin,(INT4)(SL[0].seglength/darm.deltaT +0.5));
  LALCreateVector(&status,&excwin,(INT4)(SL[0].seglength/exc.deltaT +0.5));

  winparams.type=Hann;

  /* windows for time domain channels */
  /* asq */
  winparams.length=(INT4)(SL[0].seglength/asq.deltaT +0.5);
  LALWindow(&status,asqwin,&winparams);

  /* darm */
  winparams.length=(INT4)(SL[0].seglength/darm.deltaT +0.5);
  LALWindow(&status,darmwin,&winparams);

  /*exc*/
  winparams.length=(INT4)(SL[0].seglength/exc.deltaT +0.5);
  LALWindow(&status,excwin,&winparams);


  fprintf(stdout,"# GPStime       alpha real           alpha im           alpha*beta real      alpha*beta im        exc. amplitude       darm amplitude      AS_Q amplitude\n");
  fprintf(stdout,"# ----------------------------------------------------------------------------------------------------------------------------------------------------------\n");

  fflush(stdout);


  return 0;
}

/*******************************************************************************/

int ReadFiles(struct CommandLineArgsTag CLA)
{
  char line[256];
  INT4 i;
  FILE *fpS,*fpR,*fpA,*fpSeg;
  REAL8 Cmag,Cphase,Rmag,Rphase,Amag,Aphase,freq,x,y;
  static COMPLEX16FrequencySeries R0;
  static COMPLEX16FrequencySeries C0;
  static COMPLEX16FrequencySeries A0;

  /* Allocate space for response and sensing functions; just enough to read first 1200Hz */
  LALZCreateVector( &status, &R0.data, MAXLINERS);
  LALZCreateVector( &status, &C0.data, MAXLINERS);
  LALZCreateVector( &status, &A0.data, MAXLINERS);

  /* Fill in R0, C0 data */
  R0.f0=0.0;
  R0.deltaF=1.0/64.0;                   /*ACHTUNG: HARDWIRED !!*/
  C0.f0=0.0;
  C0.deltaF=1.0/64.0;                   /*ACHTUNG: HARDWIRED !!*/
  A0.f0=0.0;
  A0.deltaF=1.0/64.0;                   /*ACHTUNG: HARDWIRED !!*/

 /* This is kinda messy... Unfortunately there's no good way of doing this */
 /* ------ Open and read Sensing file ------ */
 i=0;
 fpS=fopen(CLA.CFile,"r");
 if (fpS==NULL)
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.CFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpS))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i > MAXLINERS-1)
       {
	 /* done reading file */
	 break;
       }
     sscanf(line,"%le %le %le",&freq,&Cmag,&Cphase);
     C0.data->data[i].re=Cmag*cos(Cphase);
     C0.data->data[i].im=Cmag*sin(Cphase);
     i++;
   }
 fclose(fpS);
 /* -- close Sensing file -- */

 /* ------ Open and read Response file ------ */
 i=0;
 fpR=fopen(CLA.RFile,"r");
 if (fpR==NULL)
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.RFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpR))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i > MAXLINERS-1)
       {
	 /* done reading file */
	 break;
       }
     sscanf(line,"%le %le %le",&freq,&Rmag,&Rphase);
     R0.data->data[i].re=Rmag*cos(Rphase);
     R0.data->data[i].im=Rmag*sin(Rphase);
     i++;
   }
 fclose(fpR);
 /* -- close Sensing file -- */

 /* ------ Open and read Response file ------ */
 i=0;
 fpA=fopen(CLA.AFile,"r");
 if (fpA==NULL)
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.AFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpA))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i > MAXLINERS-1)
       {
	 /* done reading file */
	 break;
       }
     sscanf(line,"%le %le %le",&freq,&Amag,&Aphase);
     A0.data->data[i].re=Amag*cos(Aphase);
     A0.data->data[i].im=Amag*sin(Aphase);
     i++;
   }
 fclose(fpA);
 /* -- close Sensing file -- */

 /* ------ Open and read Segment file ------ */
 i=0;
 fpSeg=fopen(CLA.SegmentsFile,"r");
 if (fpSeg==NULL)
   {
     fprintf(stderr,"That's weird... %s doesn't exist!\n",CLA.SegmentsFile);
     return 1;
   }
 while(fgets(line,sizeof(line),fpSeg))
   {
     if(*line == '#') continue;
     if(*line == '%') continue;
     if (i > MAXLINESEGS-1)
       {
	 fprintf(stderr,"Too many lines in file %s! Exiting... \n", CLA.SegmentsFile);
	 return 1;
       }
     sscanf(line,"%d %d %f",&SL[i].nseg,&SL[i].tgps,&SL[i].seglength);
     i++;
   }
 numsegs=i;
 fclose(fpSeg);
 /* -- close Sensing file -- */

  /* compute C0 and R0 at correct frequency */
  /* use linear interpolation */

  x = modf( CLA.f / R0.deltaF, &y );
  i = floor( y );
  Rf0.re  = ( 1 - x ) * R0.data->data[i].re;
  Rf0.re += x * R0.data->data[i].re;
  Rf0.im  = ( 1 - x ) * R0.data->data[i].im;
  Rf0.im += x * R0.data->data[i].im;
  x = modf( CLA.f / C0.deltaF, &y );
  i = floor( y );
  Cf0.re  = ( 1 - x ) * C0.data->data[i].re;
  Cf0.re += x * C0.data->data[i].re;
  Cf0.im  = ( 1 - x ) * C0.data->data[i].im;
  Cf0.im += x * C0.data->data[i].im;
  x = modf( CLA.f / A0.deltaF, &y );
  i = floor( y );
  Af0.re  = ( 1 - x ) * A0.data->data[i].re;
  Af0.re += x * A0.data->data[i].re;
  Af0.im  = ( 1 - x ) * A0.data->data[i].im;
  Af0.im += x * A0.data->data[i].im;

 /* create Frame cache */
 LALFrCacheImport(&status,&framecache,CLA.FrCacheFile);
 LALFrCacheOpen(&status,&framestream,framecache);

 LALDestroyFrCache(&status,&framecache);

 LALZDestroyVector(&status,&R0.data);
 LALZDestroyVector(&status,&C0.data);
 LALZDestroyVector(&status,&A0.data);

 return 0;
}

/*******************************************************************************/

extern char *optarg;
   int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA)
{
  INT4 c, errflg = 0;
  optarg = NULL;

  /* Initialize default values */
  CLA->f=0.0;
  CLA->k=0.0;
  CLA->RFile="";
  CLA->CFile="";
  CLA->AFile="";
  CLA->FrCacheFile="";
  CLA->SegmentsFile="";
  CLA->exc_chan="";
  CLA->darm_chan="";
  CLA->asq_chan="";

  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt(argc, argv,"hf:F:r:S:c:A:E:D:a:k:"))!=-1))
    switch (c) {
    case 'f':
      /* frequency bandwidth */
      CLA->f=atof(optarg);
      break;
    case 'k':
      /* frequency bandwidth */
      CLA->k=atof(optarg);
      break;
    case 'F':
      /* starting observation time */
      CLA->FrCacheFile=optarg;
      break;
    case 'r':
      /* starting observation time */
      CLA->RFile=optarg;
      break;
    case 'c':
      /* starting observation time */
      CLA->CFile=optarg;
      break;
    case 'a':
      /* starting observation time */
      CLA->AFile=optarg;
      break;
    case 'S':
      /* starting observation time */
      CLA->SegmentsFile=optarg;
      break;
    case 'E':
      /* frequency bandwidth */
      CLA->exc_chan=optarg;
      break;
    case 'A':
      /* frequency bandwidth */
      CLA->asq_chan=optarg;
      break;
    case 'D':
      /* frequency bandwidth */
      CLA->darm_chan=optarg;
      break;
   case 'h':
      /* print usage/help message */
      fprintf(stdout,"All arguments are required. They are:\n");
      fprintf(stdout,"\t-f\tFLOAT\t Calibration line frequency in Hz.\n");
      fprintf(stdout,"\t-k\tFLOAT\t Darm to etmx output matrix value.\n");
      fprintf(stdout,"\t-F\tSTRING\t Name of frame cache file.\n");
      fprintf(stdout,"\t-r\tSTRING\t Name of response function file.\n");
      fprintf(stdout,"\t-c\tSTRING\t Name of sensing function file.\n");
      fprintf(stdout,"\t-a\tSTRING\t Name of actuation function file.\n");
      fprintf(stdout,"\t-S\tSTRING\t Name of segment list file.\n");
      fprintf(stdout,"\t-A\tSTRING\t AS_Q channel name (eg, L1:LSC-AS_Q).\n");
      fprintf(stdout,"\t-E\tSTRING\t Excitation channel name (eg, L1:LSC-ETMX_EXC_DAQ)\n");
      fprintf(stdout,"\t-D\tSTRING\t Darm channel name (eg, L1:LSC-DARM_CTRL)\n");
      fprintf(stdout,"Caution: There are a few hardwired quantities in this code. Look for 'HARDWIRED!!' in the code.\n");
      exit(0);
      break;
    default:
      /* unrecognized option */
      errflg++;
      fprintf(stderr,"Unrecognized option argument %c\n",c);
      exit(1);
      break;
    }
  if(CLA->f == 0)
    {
      fprintf(stderr,"No calibration line frequency specified.\n");
      fprintf(stderr,"Try ./CCFDriver -h \n");
      exit( 77 );
    }
  if(CLA->k == 0)
    {
      fprintf(stderr,"No value of the output matrix to x arm specified.\n");
      fprintf(stderr,"Try ./CCFDriver -h \n");
      exit( 77 );
    }
  if(CLA->CFile == "")
    {
      fprintf(stderr,"No sensing function file specified.\n");
      fprintf(stderr,"Try ./CCFDriver -h \n");
      exit( 77 );
    }
  if(CLA->RFile == "")
    {
      fprintf(stderr,"No response function file specified.\n");
      fprintf(stderr,"Try ./CCFDriver -h \n");
      exit( 77 );
    }
  if(CLA->AFile == "")
    {
      fprintf(stderr,"No actuation function file specified.\n");
      fprintf(stderr,"Try ./CCFDriver -h \n");
      exit( 77 );
    }
  if(CLA->FrCacheFile=="")
    {
      fprintf(stderr,"No frame cache file specified.\n");
      fprintf(stderr,"Try ./CCFDriver -h \n");
      exit( 77 );
    }
  if(CLA->SegmentsFile=="")
    {
      fprintf(stderr,"No segments file specified.\n");
      fprintf(stderr,"Try ./CCFDriver -h \n");
      exit( 77 );
    }
   if(CLA->exc_chan == "")
    {
      fprintf(stderr,"No excitation channel specified.\n");
      fprintf(stderr,"Try ./CCFDriver -h \n");
      exit( 77 );
    }
   if(CLA->darm_chan == "")
    {
      fprintf(stderr,"No darm channel specified.\n");
      fprintf(stderr,"Try ./CCFDriver -h \n");
      exit( 77 );
    }
   if(CLA->asq_chan == "")
    {
      fprintf(stderr,"No asq channel specified.\n");
      fprintf(stderr,"Try ./CCFDriver -h \n");
      exit( 77 );
    }



  return errflg;
}

/*******************************************************************************/

int FreeMem(void)
{

  LALFrClose(&status,&framestream);

  LALDestroyVector(&status,&darm.data);
  LALDestroyVector(&status,&exc.data);
  LALDestroyVector(&status,&asq.data);

  LALDestroyVector(&status,&asqwin);
  LALDestroyVector(&status,&darmwin);
  LALDestroyVector(&status,&excwin);

  LALCheckMemoryLeaks();

  return 0;
}
#endif

