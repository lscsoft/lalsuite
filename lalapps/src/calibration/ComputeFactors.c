/*
*  Copyright (C) 2007 Jolien Creighton, Xavier Siemens
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
/*                         Calibration factor generation code                    */
/*                                                                               */
/*                                  X. Siemens                                   */
/*                                                                               */
/*                 Albert Einstein Institute/UWM - October 2003                  */
/*********************************************************************************/

#include <config.h>
#if !defined HAVE_LIBLALFRAME
#include <stdio.h>
int main(void) {fputs("disabled, no lal frame library support.\n", stderr);return 1;}
#elif !defined HAVE_FRAMEL_H
#include <stdio.h>
int main(void) {fputs("disabled, no frame library support.\n", stderr);return 1;}
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

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/AVFactories.h>
#include <lal/LALCache.h>
#include <lal/FrameStream.h>
#include <lal/Window.h>
#include <lal/Calibration.h>
#include <lal/LALConstants.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALFrameL.h>

#include <series.h>

extern char *optarg;
extern int optind, opterr, optopt;

#define MAXLINESEGS 10000                 /* Maximum number of science segments */

#define A_CHANNEL "CAL-CAV_FAC"  /* name of alpha channel in frames */
#define AB_CHANNEL "CAL-OLOOP_FAC"  /* name of alphabeta channel in frames */

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)


/***************************************************************************/

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL8 f;                 /* Frequency of the calibration line */
  REAL8 t;                 /* Time interval to calculate factors */
  REAL8 G0Re;              /* Real part of open loop gain at cal line freq.*/
  REAL8 G0Im;              /* Imaginary part of open loop gain at cal line freq. */
  REAL8 D0Re;              /* Real part of digital filter at cal line freq.*/
  REAL8 D0Im;              /* Imaginary part of digital filter at cal line freq. */
  REAL8 W0Re;              /* Real part of whitening filter at cal line freq.*/
  REAL8 W0Im;              /* Imaginary part of whitening filter at cal line freq. */
  char *FrCacheFile;       /* Frame cache file */
  char *SegmentsFile;      /* Text file with the segments */
  char *exc_chan;          /* excitation channel name */
  char *darm_chan;         /* darm channel name */
  char *asq_chan;          /* asq channel name */
  char *alphafile;         /* file to store factors */
  int   outputframes;      /* flag to indicate that frame file output is requested */
  char *version;           /* version of this calibration: needed for frame file output name */
} CommandLineArgs;

typedef
struct SegmentListTag {
  INT4 gpsstart;          /* GPS Seconds of start of segment group */
  INT4 gpsend;            /* number of segments starting at tgps */
  INT4 seglength;         /* length of segment in seconds */
} SegmentList;


typedef
struct GlobalVariablesTag {
COMPLEX16 Rf0,Cf0,Af0;              /* Response, sensing and actuation function values at the frequency of the calibration line */
SegmentList SL[MAXLINESEGS];        /* Structure containing science segement info */
INT4 numsegs;                       /* number of science segments */
LIGOTimeGPS gpsepoch;               /* Global variables epoch and duration used to calculate the calibration factors */
INT4 duration;                      /* they are set every science segment to GPS start time and duration of segment */
char ifo[3];                        /* interferometer name: needed for frame file channels and output name */
} GlobalVariables;

/***************************************************************************/

/* GLOBAL VARIABLES */

static LALStatus status;
LALCache *framecache;                                           /* frame reading variables */
FrStream *framestream=NULL;

GlobalVariables GV;   /* A bunch of stuff is stored in here; mainly to protect it from accidents */

/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);

/* Reads time segment, cache, response,sensing and actuation files */
int ReadFiles(struct CommandLineArgsTag CLA);

/* Produces a table of calibration factors for a science segment;
the factors are calculated every interval CLA.t */
int GetFactors(struct CommandLineArgsTag CLA);

/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
int i;


  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;
  if (ReadFiles(CommandLineArgs)) return 3;

  /* get ifo name from one of the channel names (choose as_q) */
  strncpy( GV.ifo, CommandLineArgs.asq_chan, 2 );
  GV.ifo[2] = 0;

  /* create Frame cache, open frame stream and delete frame cache */
  framecache = XLALCacheImport(CommandLineArgs.FrCacheFile);
  LALFrCacheOpen(&status,&framestream,framecache);
  TESTSTATUS( &status );
  XLALDestroyCache(framecache);

  for(i=0;i<GV.numsegs;i++)
    {
      GV.gpsepoch.gpsSeconds=GV.SL[i].gpsstart; /* Set global variable epoch */
      GV.gpsepoch.gpsNanoSeconds=0;
      GV.duration=GV.SL[i].seglength;

      if(GetFactors(CommandLineArgs)) return 5;
    }


  LALFrClose(&status,&framestream);
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();

  return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/


/*  FUNCTIONS */
/*******************************************************************************/

int GetFactors(struct CommandLineArgsTag CLA)
{
  struct series a;
  struct series ab;
  char a_name[64];
  char ab_name[64];

FrPos pos1;

static REAL4TimeSeries darm;
static REAL4TimeSeries asq;
static REAL4TimeSeries exc;

static FrChanIn chanin_darm;
static FrChanIn chanin_asq;
static FrChanIn chanin_exc;

CalFactors factors;
UpdateFactorsParams params;

REAL4Window *asqwin=NULL,*excwin=NULL,*darmwin=NULL;  /* windows */

INT4 k,m;
LIGOTimeGPS localgpsepoch=GV.gpsepoch; /* Local variable epoch used to calculate the calibration factors */
LIGOTimeGPS tend; /* end time */
INT8 tendns; /* end time in nanoseconds */
long double gtime=(long double)(localgpsepoch.gpsSeconds+(long double)localgpsepoch.gpsNanoSeconds*1E-9);

FILE *fpAlpha=NULL;

  chanin_asq.type  = ADCDataChannel;
  chanin_darm.type = ADCDataChannel;
  chanin_exc.type  = ADCDataChannel;

  chanin_asq.name  = CLA.asq_chan;
  chanin_darm.name = CLA.darm_chan;
  chanin_exc.name  = CLA.exc_chan;

  /* Get channel time step size by calling LALFrGetREAL4TimeSeries */
  LALFrSeek(&status,&localgpsepoch,framestream);
  TESTSTATUS( &status );
  LALFrGetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&asq,&chanin_asq,framestream);
  TESTSTATUS( &status );
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&darm,&chanin_darm,framestream);
  TESTSTATUS( &status );
  LALFrSetPos(&status,&pos1,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&exc,&chanin_exc,framestream);
  TESTSTATUS( &status );


  /* Allocate space for data vectors */
  LALCreateVector(&status,&asq.data,(UINT4)(CLA.t/asq.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&darm.data,(UINT4)(CLA.t/darm.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&exc.data,(UINT4)(CLA.t/exc.deltaT +0.5));
  TESTSTATUS( &status );


  /* windows for time domain channels */
  /* asq */
  asqwin = XLALCreateHannREAL4Window((UINT4)(CLA.t/asq.deltaT +0.5));
  if( !asqwin ) return -1;

  /* darm */
  darmwin = XLALCreateHannREAL4Window((UINT4)(CLA.t/darm.deltaT +0.5));
  if( !darmwin ) return -1;

  /* exc */
  excwin = XLALCreateHannREAL4Window((UINT4)(CLA.t/exc.deltaT +0.5));
  if( !excwin ) return -1;

  /* setup series for frame file output */
  snprintf( a_name, sizeof( a_name ), "%s:" A_CHANNEL, GV.ifo );
  snprintf( ab_name, sizeof( ab_name ), "%s:" AB_CHANNEL, GV.ifo );
  tendns  = (INT8)(1000000000) * (INT8)(GV.gpsepoch.gpsSeconds) + (INT8)(GV.gpsepoch.gpsNanoSeconds);
  tendns += (INT8)(1e9 * GV.duration);
  tend.gpsSeconds = tendns / (INT8)(1000000000);
  tend.gpsNanoSeconds = tendns % (INT8)(1000000000);
  a.name  = a_name;
  ab.name = ab_name;
  a.unit = "none";
  ab.unit = "none";
  a.dom   = ab.dom  = Time;
  a.type  = ab.type = FR_VECT_8C;
  a.tbeg  = ab.tbeg = GV.gpsepoch;
  a.tend  = ab.tend = tend;
  a.step  = ab.step = CLA.t;
  a.size  = ab.size = (INT4)(GV.duration/CLA.t);
  a.data  = calloc( 2 * a.size, sizeof( *a.data ) );
  ab.data = calloc( 2 * ab.size, sizeof( *ab.data ) );
  if ( ! a.data || ! ab.data )
  {
    fprintf(stderr,"Memory allocation error!\n");
    return 1;
  }

  /* Open user input factors file */
  fpAlpha=fopen(CLA.alphafile,"w");
  if (fpAlpha==NULL)
    {
      fprintf(stderr,"Could not open %s!\n",CLA.alphafile);
      return 1;
    }
  setvbuf( fpAlpha, NULL, _IONBF, 0 );  /* This stops buffering -- From Duncan Brown*/

  fprintf(fpAlpha,"%s Re(Go) = %e, Im(Go)= %e \n", "%",CLA.G0Re, CLA.G0Im);
  fprintf(fpAlpha,"%s Re(Do) = %e, Im(Do)= %e \n", "%",CLA.D0Re, CLA.D0Im);
  fprintf(fpAlpha,"%s Re(Wo) = %e, Im(Wo)= %e \n", "%",CLA.W0Re, CLA.W0Im);
  fprintf(fpAlpha,"%s  GPS Time          Re(a)    Im(a)     Re(b)    Im(b)    Re(ab)    Im(ab)   Re(AS_Q) Im(AS_Q) Re(DARM) Im(DARM) Re(EXC)  Im(EXC)\n", "%");

  for(m=0;m < (INT4)(GV.duration/CLA.t);m++)
    {
      /* Fill data vectors with data */
      LALFrSeek(&status,&localgpsepoch,framestream);
      TESTSTATUS( &status );
      LALFrGetPos(&status,&pos1,framestream);
      TESTSTATUS( &status );

      LALFrSetPos(&status,&pos1,framestream);
      TESTSTATUS( &status );
      LALFrGetREAL4TimeSeries(&status,&asq,&chanin_asq,framestream);
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
      for(k=0;k<(INT4)(CLA.t/asq.deltaT +0.5);k++)
	{
	  asq.data->data[k] *= 2.0*asqwin->data->data[k];
	}
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
      params.asQ = &asq;
      params.exc = &exc;
      params.lineFrequency = CLA.f;
      params.openloop.real_FIXME =  CLA.G0Re;
      params.openloop.imag_FIXME =  CLA.G0Im;
      params.digital.real_FIXME = CLA.D0Re;
      params.digital.imag_FIXME = CLA.D0Im;
      params.whitener.real_FIXME = CLA.W0Re;
      params.whitener.imag_FIXME = CLA.W0Im;

      LALComputeCalibrationFactors(&status,&factors,&params);
      TESTSTATUS( &status );

      fprintf(fpAlpha,"%18.9Lf %f %f %f %f %f %f %f %f %f %f %f %f \n",gtime,
	      creal(factors.alpha),cimag(factors.alpha),
	      creal(factors.beta),cimag(factors.beta),
	      creal(factors.alphabeta),cimag(factors.alphabeta),
	      creal(factors.asq)*2/CLA.t,cimag(factors.asq)*2/CLA.t,
	      creal(factors.darm)*2/CLA.t,cimag(factors.darm)*2/CLA.t,
	      creal(factors.exc)*2/CLA.t,cimag(factors.exc)*2/CLA.t);

      /* put factors into series for frame output */
      a.data[2*m]    = creal(factors.alpha);
      a.data[2*m+1]  = cimag(factors.alpha);
      ab.data[2*m]   = creal(factors.alphabeta);
      ab.data[2*m+1] = cimag(factors.alphabeta);

      gtime += CLA.t;
      localgpsepoch.gpsSeconds = (INT4)gtime;
      localgpsepoch.gpsNanoSeconds = (INT4)((gtime-(INT4)gtime)*1E+09);
    }

  /* output frame data if desired */
  if ( CLA.outputframes )
  {
    struct FrFile *frfile = NULL;
    struct FrameH *frame = NULL;
    char frfilename[256];
    sprintf( frfilename, "%c-CAL_FAC_%s_%s-%d-%d.gwf", GV.ifo[0], CLA.version, GV.ifo, GV.gpsepoch.gpsSeconds, (int)ceil(GV.duration) );
    frfile = FrFileONew( frfilename, 0 );
    frame = fr_add_proc_data( frame, &a );
    frame = fr_add_proc_data( frame, &ab );
    FrameWrite( frame, frfile );
    FrFileOEnd( frfile );
  }

  /* Clean up */
  free( a.data );
  free( ab.data );
  LALDestroyVector(&status,&darm.data);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&exc.data);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&asq.data);
  TESTSTATUS( &status );

  XLALDestroyREAL4Window(asqwin);
  XLALDestroyREAL4Window(darmwin);
  XLALDestroyREAL4Window(excwin);

  fclose(fpAlpha);

  return 0;
}

/*******************************************************************************/

int ReadFiles(struct CommandLineArgsTag CLA)
{
  char line[256];
  INT4 i;
  FILE *fpSeg;

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
     sscanf(line,"%d %d %d",&GV.SL[i].gpsstart,&GV.SL[i].gpsend,&GV.SL[i].seglength);
     i++;
   }
 GV.numsegs=i;
 fclose(fpSeg);
 /* -- close Sensing file -- */

 return 0;
}

/*******************************************************************************/

int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA)
{
  INT4 c, errflg = 0;
  optarg = NULL;

  /* Initialize default values */
  CLA->f=0.0;
  CLA->t=0.0;
  CLA->G0Re=0.0;
  CLA->G0Im=0.0;
  CLA->D0Re=0.0;
  CLA->D0Im=0.0;
  CLA->W0Re=0.0;
  CLA->W0Im=0.0;
  CLA->FrCacheFile=NULL;
  CLA->SegmentsFile=NULL;
  CLA->alphafile=NULL;
  CLA->exc_chan=NULL;
  CLA->darm_chan=NULL;
  CLA->asq_chan=NULL;
  CLA->alphafile=NULL;
  CLA->outputframes=0;
  CLA->version=NULL;

  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt(argc, argv,"hof:F:S:A:E:D:b:t:i:j:k:l:m:n:v:"))!=-1))
    switch (c) {
    case 'f':
      /* calibration line frequency */
      CLA->f=atof(optarg);
      break;
    case 't':
      /* calibration line frequency */
      CLA->t=atof(optarg);
      break;
    case 'i':
      /* calibration line frequency */
      CLA->G0Re=atof(optarg);
      break;
    case 'j':
      /* calibration line frequency */
      CLA->G0Im=atof(optarg);
      break;
    case 'k':
      /* calibration line frequency */
      CLA->D0Re=atof(optarg);
      break;
    case 'l':
      /* calibration line frequency */
      CLA->D0Im=atof(optarg);
      break;
    case 'm':
      /* calibration line frequency */
      CLA->W0Re=atof(optarg);
      break;
    case 'n':
      /* calibration line frequency */
      CLA->W0Im=atof(optarg);
      break;
    case 'F':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      break;
    case 'S':
      /* name of segments file */
      CLA->SegmentsFile=optarg;
      break;
    case 'E':
      /* name of excitation channel */
      CLA->exc_chan=optarg;
      break;
    case 'A':
      /* name of as_q channel */
      CLA->asq_chan=optarg;
      break;
    case 'D':
      /* name of darm channel */
      CLA->darm_chan=optarg;
      break;
    case 'b':
      /* name of darm channel */
      CLA->alphafile=optarg;
      break;
    case 'o':
      /* output frame files */
      CLA->outputframes=1;
      break;
    case 'v':
      /* calibration version string */
      CLA->version=optarg;
      break;
   case 'h':
      /* print usage/help message */
      fprintf(stdout,"Flags and optional arguments are:\n");
      fprintf(stdout,"\t-o\t     \t Output frame files.\n");
      fprintf(stdout,"\t-v\tSTRING\t Calibration version for frame files.  Required with -o flag.\n");
      fprintf(stdout,"All of the following arguments are required. They are:\n");
      fprintf(stdout,"\t-f\tFLOAT\t Calibration line frequency in Hz.\n");
      fprintf(stdout,"\t-t\tFLOAT\t Time interval to calculate factors in seconds (>0.0005).\n");
      fprintf(stdout,"\t-i\tFLOAT\t Real part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\t-j\tFLOAT\t Imaginary part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\t-k\tFLOAT\t Real part of the digital filter at the calibration line frequency.\n");
      fprintf(stdout,"\t-l\tFLOAT\t Imaginary part of digital filter at the calibration line frequency.\n");
      fprintf(stdout,"\t-m\tFLOAT\t Real part of the whitening filter at the calibration line frequency.\n");
      fprintf(stdout,"\t-n\tFLOAT\t Imaginary part of whitening filter at the calibration line frequency.\n");
      fprintf(stdout,"\t-F\tSTRING\t Name of frame cache file.\n");
      fprintf(stdout,"\t-S\tSTRING\t Name of segment list file.\n");
      fprintf(stdout,"\t-A\tSTRING\t AS_Q channel name (eg, L1:LSC-AS_Q).\n");
      fprintf(stdout,"\t-E\tSTRING\t Excitation channel name (eg, L1:LSC-ETMX_EXC_DAQ)\n");
      fprintf(stdout,"\t-D\tSTRING\t Darm channel name (eg, L1:LSC-DARM_CTRL)\n");
      fprintf(stdout,"\t-b\tSTRING\t Output file for calibration factors.\n");
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
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
  if(CLA->t == 0)
    {
      fprintf(stderr,"No time interval to calculate factors.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
  if(CLA->t < 0.0005)
    {
      fprintf(stderr,"Time interval to calculate factors too small (<0.0005).\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
  if(CLA->G0Re == 0.0 )
    {
      fprintf(stderr,"No real part of open loop gain specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
  if(CLA->G0Im == 0.0 )
    {
      fprintf(stderr,"No imaginary part of open loop gain specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
  if(CLA->D0Re == 0.0 )
    {
      fprintf(stderr,"No real part of digital filter specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
  if(CLA->D0Im == 0.0 )
    {
      fprintf(stderr,"No imaginary part of digital filter specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
  if(CLA->alphafile == NULL)
    {
      fprintf(stderr,"No output factors file specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
  if(CLA->FrCacheFile == NULL)
    {
      fprintf(stderr,"No frame cache file specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
  if(CLA->SegmentsFile == NULL)
    {
      fprintf(stderr,"No segments file specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
   if(CLA->exc_chan == NULL)
    {
      fprintf(stderr,"No excitation channel specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
   if(CLA->darm_chan == NULL)
    {
      fprintf(stderr,"No darm channel specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
   if(CLA->asq_chan == NULL)
    {
      fprintf(stderr,"No asq channel specified.\n");
      fprintf(stderr,"Try ./ComputeFactors -h \n");
      return 1;
    }
   if (CLA->outputframes)
   {
     if (!CLA->version)
     {
       fprintf(stderr,"Must specify calibration version string.\n");
       fprintf(stderr,"Try ./ComputeFactors -h \n");
       return 1;
     }
   }

  return errflg;
}

/*******************************************************************************/
#endif
