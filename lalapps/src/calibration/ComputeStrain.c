/*********************************************************************************/
/*                         Calibrated h(t) generation code                       */
/*                                                                               */
/*                            B. Allen and X. Siemens                            */
/*                                                                               */
/*                 Albert Einstein Institute/UWM - October 2003                  */
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

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include "filters-H1-S3.h"                      /* files that contains the filter coefficients */
#define CHANNEL "H1:Calibrated-Strain"
#define DATADIR "/scratch4/xavi/hoft/S3/H1/H"
#define FRAMETYPE "H1_RDS_C02_LX"

extern char *optarg;
extern int optind, opterr, optopt;

#define MAXLINERS 76800                   /* Max lines read in Freq Response files */
#define MAXLINESEGS 10000                 /* Maximum number of science segments */
#define To 1                             /* length of time used in table of alphas and betas (seconds)*/
#define T  16                             /* length of time calibrated per iteration (seconds) */
#define SR 16384                          /* Sampling rate of data (Hz) */
#define USR 16                            /* Upsampling factor */
#define MAXALPHAS 100000                    /* Maximum number of calibration factors in a science segment */

#define SLUDGEFACTOR 32764                 /* Number of samples in AS_Q to keep as extra sludge for smooth filtering */

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 100; } else ((void)0)



/***************************************************************************/

/* STRUCTURES */
struct CommandLineArgsTag {
  REAL8 f;                 /* Frequency of the calibration line */
  REAL8 G0Re;              /* Real part of open loop gain at cal line freq.*/
  REAL8 G0Im;              /* Imaginary part of open loop gain at cal line freq. */
  REAL8 D0Re;              /* Real part of digital filter at cal line freq.*/
  REAL8 D0Im;              /* Imaginary part of digital filter at cal line freq. */
  char *FrCacheFile;       /* Frame cache file */
  char *SegmentsFile;      /* Text file with the segments */
  char *exc_chan;          /* excitation channel name */    
  char *darm_chan;         /* darm channel name */ 
  char *asq_chan;          /* asq channel name */
  char *alphafile;          /* asq channel name */
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
REAL8 ta_interp[MAXALPHAS],alpha[MAXALPHAS],beta[MAXALPHAS];   /* Arrays that contain the calibration factors computed in To time intervals */
LIGOTimeGPS gpsepoch;               /* Global variables epoch and duration used to calculate the calibration factors */
INT4 duration;                      /* they are set every science segment to GPS start time and duration of segment */
MyIIRFilter Cinv,G[NGfilt],AA,AX[NAXfilt],AY[NAYfilt];    /* Inverse sensing, servo, analog actuation, digital x actuation  digital y actuation */
MyIIRFilter ADW, AADW;  /* dewhitening and anti dewhiytwening filters in actuation */
REAL8TimeSeries AS_Q,hR,hC,uphR,hCx,hCy,h;                /* AS_Q, Residual, control and total strains; up... time series are upsampled */
REAL8TimeSeries hipasssludge,lopasssludge;                /* Sludge for high and low pass filtering of data */
REAL8TimeSeries hipasssludgeR,hipasssludgeC;              /* Sludge for high filtering of contral and residual signals */
} GlobalVariables;

/***************************************************************************/

/* GLOBAL VARIABLES */

static LALStatus status;
INT4 lalDebugLevel=3;
FrCache *framecache;                                           /* frame reading variables */
FrStream *framestream=NULL;

GlobalVariables GV;   /* A bunch of stuff is stored in here; mainly to protect it from accidents */

/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);

/* Reads time segment, cache, response,sensing and actuation files */
int ReadFiles(struct CommandLineArgsTag CLA);             

/* Create the filters from the filters.h file */
int MakeFilters(void);

/* Create the filters from the filters.h file */
int DeleteFilterHistories(void);

/* Produces a table of calibration factors for a science segment; 
the factors are calculated every interval To*/
int GetFactors(struct CommandLineArgsTag CLA);   

/* Allocates space for data time series AS_Q and DARM_CTRL */
int AllocateData(struct CommandLineArgsTag CLA, int DT);

/* Reads T seconds of AS_Q and DARM_CTRL data */
int ReadDTsecondsofData(struct CommandLineArgsTag CLA, int DT);

/* Multiplies each sample of AS_Q by 1/alpha */
int hROverAlpha(int i);                                        

/* Multiplies each sample of AS_Q by beta */
int hCTimesBeta(int i);

/* Upsamples AS_QR storing result in upAS_Q by putting SR-1 zeros between samples in AS_Q */
int UpsamplehR(int up_factor);   

/* Low passes upsampled AS_Q doing all the sludge accouting */
int LowPasshR(int j);

/* My own filtering function that sets the history correctly */
int FilterSeries(MyIIRFilter *F, REAL8TimeSeries *TSeries, int sludge);

/* Downsamples calibrated AS_Q storing result in hr by taking 1 sample out of USR in upsampled series */
int DownsamplehR(void);   

/* High passes AS_Q data doing all the sludge accouting */
int HighPass1(int j);

/* High passes residual and control signals doing all the sludge accouting */
int HighPass2(int j);

/* Frees the memory */
int DeAllocateData(void);                                        

/* Frees the memory */
int FreeMem(void);                                        


/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
int i,j,DT,p; 

  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;
  if (ReadFiles(CommandLineArgs)) return 3;
  if (MakeFilters()) return 4;

  /* create Frame cache, open frame stream and delete frame cache */
  LALFrCacheImport(&status,&framecache,CommandLineArgs.FrCacheFile);
  TESTSTATUS( &status );
  LALFrCacheOpen(&status,&framestream,framecache);
  TESTSTATUS( &status );
  LALDestroyFrCache(&status,&framecache);
  TESTSTATUS( &status );

  /* Create sludge vectors for BW filtering */
  GV.hipasssludge.deltaT=1.0/SR;
  LALDCreateVector(&status,&GV.hipasssludge.data,SLUDGEFACTOR);
  TESTSTATUS( &status );
  GV.hipasssludgeR.deltaT=1.0/SR;
  LALDCreateVector(&status,&GV.hipasssludgeR.data,SLUDGEFACTOR);
  TESTSTATUS( &status );
  GV.hipasssludgeC.deltaT=1.0/SR;
  LALDCreateVector(&status,&GV.hipasssludgeC.data,SLUDGEFACTOR);
  TESTSTATUS( &status );
  GV.lopasssludge.deltaT=1.0/SR/USR;
  LALDCreateVector(&status,&GV.lopasssludge.data,USR*SLUDGEFACTOR);
  TESTSTATUS( &status );

  for(i=0;i<GV.numsegs;i++)
    {
      if(GV.SL[i].seglength < 3*To)
	{
	  fprintf(stderr,"Throwing out segment %d, it is shorter than %d seconds\n",GV.SL[i].gpsstart,3*To);
	  continue;
	}

      GV.gpsepoch.gpsSeconds=GV.SL[i].gpsstart; /* Set global variable epoch */
      GV.gpsepoch.gpsNanoSeconds=0;
      GV.duration=GV.SL[i].seglength;

      if(GetFactors(CommandLineArgs)) return 5;

      for(j=0; j <= GV.duration/T; j++) 
	{

	  /* Normally we read in and filter T seconds of data */
	  DT=T;
	  /* But the very last bit may be a different length (< T seconds) */
	  if(j == GV.duration/T) DT=GV.duration-(GV.duration/T)*T;
	  /* If by chance the duration is a multiple of T seconds then the 
	     last chunk is zero seconds long and we should exit */
	  if( DT == 0) break;

	  /*  fprintf(stdout,"%d %d  %d  %d\n",j,GV.gpsepoch.gpsSeconds,GV.gpsepoch.gpsSeconds+DT,DT); */

	  /* Allocates space for data */
	  if (AllocateData(CommandLineArgs,DT)) return 4;

	  /* Reads T seconds of AS_Q data */
	  if (ReadDTsecondsofData(CommandLineArgs,DT)) return 4;

	  /* hipass AS_Q with Butterworth filter */
	  if(HighPass1(j)) return 4;

	  /* Copy AS_Q into residual and control signal time series */  
	  for (p=0; p<(int)GV.hR.data->length; p++) {
	    GV.hR.data->data[p]=GV.AS_Q.data->data[p];
	  }
	  for (p=0; p<(int)GV.hC.data->length; p++) {
	    GV.hC.data->data[p]=GV.AS_Q.data->data[p];
	  }

	  /******** COMPUTE RESIDUAL STRAIN **********/
	  
	  /* multiply hR by 1/alpha(t) */
	  if (hROverAlpha(i)) return 4;

	  /* Upsample hR by a factor of USR defined above */
	  if (UpsamplehR(USR)) return 4;

	  /* smooth with low pass Butterworth filter at 7000Hz */
	  if(LowPasshR(j)) return 4;

          /* Filter through inverse of sensing function */
	  if (FilterSeries(&GV.Cinv,&GV.uphR, SLUDGEFACTOR*USR)) return 5;

	  /* Downsample by a factor of 16 */
	  if(DownsamplehR()) return 3;


	  /******** COMPUTE CONTROL STRAIN **********/

	  /* multiply hC by beta(t) */
	  if (hCTimesBeta(i)) return 4; 

          /* Filter hC through servo to get DARM_CTRL */
	  for(p=NGfilt-1;p>=0;p--){
	    if (FilterSeries(&GV.G[p],&GV.hC, SLUDGEFACTOR)) return 5;
	  }
	  /* Adjust to account for servo gain */ 
	  for (p=0; p<(int)GV.hC.data->length;p++) {
	    GV.hC.data->data[p]= ServoGain*GV.hC.data->data[p];
 	  }

	  /* Copy data into x and y time series for parallel filtering */
	  for (p=0; p<(int)GV.hCx.data->length; p++) {
	    GV.hCx.data->data[p]=GV.hC.data->data[p];
	  }
	  for (p=0; p<(int)GV.hCy.data->length; p++) {
	    GV.hCy.data->data[p]=GV.hC.data->data[p];
	  }

	  /* Filter x-arm */
	  for(p=NAXfilt-1;p>=0;p--){
	    if (FilterSeries(&GV.AX[p],&GV.hCx,SLUDGEFACTOR)) return 5;
	  }
	  /* Adjust to account for digital gain on x-arm*/ 
	  for (p=0; p<(int)GV.hC.data->length;p++) {
	    GV.hCx.data->data[p] *= AXGain;
 	  }

	  /* Filter y-arm */
	  for(p=NAYfilt-1;p>=0;p--){
	    if (FilterSeries(&GV.AY[p],&GV.hCy,SLUDGEFACTOR)) return 5;
	  }
	  /* Adjust to account for digital gain on y-arm*/ 
	  for (p=0; p<(int)GV.hC.data->length;p++) {
	    GV.hCy.data->data[p] *= AYGain;
 	  }

	  /* add x-arm and y-arm together */
	  for (p=0; p<(int)GV.hC.data->length; p++) {
	    GV.hC.data->data[p]=(GV.hCx.data->data[p]+GV.hCy.data->data[p])/2;
	  }

 	  /* filter through analog part of actuation */ 
	  if (FilterSeries(&GV.AA,&GV.hC,SLUDGEFACTOR)) return 5;

 	  /* filter through Dewhitener in actuation */ 
	  if (FilterSeries(&GV.ADW,&GV.hC,SLUDGEFACTOR)) return 5;

 	  /* filter through Anti-Dewhitener in actuation */ 
	  if (FilterSeries(&GV.AADW,&GV.hC,SLUDGEFACTOR)) return 5;


	  /******** COMPUTE NET CALIBRATED STRAIN **********/

	  /* hipass control and residual strains with Butterworth filter before adding them together */
	  if(HighPass2(j)) return 4;

	  /* add control and residual signals together */
	  for (p=0; p<(int)GV.h.data->length; p++) {
	    GV.h.data->data[p]= GV.hR.data->data[p]+ GV.hC.data->data[p];
	  }

	  /* WRITE A FRAME */
	  strncpy( GV.h.name, CHANNEL, sizeof( GV.h.name ) );
	  GV.h.epoch.gpsSeconds=GV.gpsepoch.gpsSeconds;
	  
	  {
	    char filename[256];
/* 	    char nodenumber[16]; */

	    strcpy(filename,DATADIR);
/* 	    sprintf(nodenumber,"s%03d/H",(GV.gpsepoch.gpsSeconds-751654515)/((757699249-751654515)/296)+1); */
/* 	    strcat(filename,nodenumber);       */

	    {
	      FrOutPar opar = { filename, FRAMETYPE, ProcDataChannel, 1, 0, 2 };
 
	      LALFrWriteREAL8TimeSeries( &status, &GV.h, &opar );
	      TESTSTATUS( &status );

	      fprintf(stdout,"%s-%s-%d-%d.gwf\n",filename,FRAMETYPE,GV.h.epoch.gpsSeconds,DT);
	    }

	  }


	  
	  GV.gpsepoch.gpsSeconds += DT;

	  if (DeAllocateData()) return 4;

	}
      /* Delete filter histories */
      if (DeleteFilterHistories()) return 4;
    }
  
  if(FreeMem()) return 8;

  return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/


/*  FUNCTIONS */

/*******************************************************************************/

int DeAllocateData(void)
{

  LALDDestroyVector(&status,&GV.AS_Q.data);
  TESTSTATUS( &status );

  LALDDestroyVector(&status,&GV.hC.data);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&GV.hCx.data);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&GV.hCy.data);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&GV.hR.data);
  TESTSTATUS( &status );

  LALDDestroyVector(&status,&GV.h.data);
  TESTSTATUS( &status );

  LALDDestroyVector(&status,&GV.uphR.data);
  TESTSTATUS( &status );

  return 0;
}

/*******************************************************************************/

int DeleteFilterHistories(void)
{
  int l,n;

  /* Inverse of sensing function */
  for(l=0;l<GV.Cinv.yOrder-1;l++) GV.Cinv.yhist[l]=0.0;
  for(l=0;l<GV.Cinv.xOrder-1;l++) GV.Cinv.xhist[l]=0.0;

  /* Servo */
  for(n=0;n<NGfilt;n++){
    for(l=0;l<GV.G[n].yOrder-1;l++) GV.G[n].yhist[l]=0.0;
    for(l=0;l<GV.G[n].xOrder-1;l++) GV.G[n].xhist[l]=0.0;
  }

  /* actuation: analog part */
  for(l=0;l<GV.AA.yOrder-1;l++) GV.AA.yhist[l]=0.0;
  for(l=0;l<GV.AA.xOrder-1;l++) GV.AA.xhist[l]=0.0;

  /* actuation: x-arm digital filters */
  for(n=0;n<NAXfilt;n++){
    for(l=0;l<GV.AX[n].yOrder-1;l++) GV.AX[n].yhist[l]=0.0;
    for(l=0;l<GV.AX[n].xOrder-1;l++) GV.AX[n].xhist[l]=0.0;
  }

  /* actuation: y-arm digital filters */
  for(n=0;n<NAYfilt;n++){
    for(l=0;l<GV.AY[n].yOrder-1;l++) GV.AY[n].yhist[l]=0.0;
    for(l=0;l<GV.AY[n].xOrder-1;l++) GV.AY[n].xhist[l]=0.0;
  }

  return 0;
}

/*******************************************************************************/

int LowPasshR(int j)
{
  int n;
  PassBandParamStruc filterpar;

  filterpar.name  = "Butterworth Low Pass";
  filterpar.nMax  = 12;
  filterpar.f2    = -1.0;
  filterpar.a2    = -1.0;
  filterpar.f1    = 6000.0;
  filterpar.a1    = 0.5;

  LALButterworthREAL8TimeSeries(&status,&GV.uphR,&filterpar);
  TESTSTATUS( &status );

  /* If this chunk of size T is not the first within a science segment copy half of slugde to start of hR */
  if(j>0){  
    for (n=0; n<USR*SLUDGEFACTOR/2;n++) { 
      GV.uphR.data->data[n]=GV.lopasssludge.data->data[n]; 
    }  
  }  
  /* copy sludge at the end of AS_Q to hipasssludge */
  for (n=0; n<USR*SLUDGEFACTOR;n++) { 
    GV.lopasssludge.data->data[n]=GV.uphR.data->data[n+GV.uphR.data->length-USR*SLUDGEFACTOR];
  }

  return 0;
}

/*******************************************************************************/

int HighPass1(int j)
{
  int n;
  PassBandParamStruc filterpar;

  filterpar.name  = "Butterworth High Pass";
  filterpar.nMax  = 10;
  filterpar.f2    = 40.0;
  filterpar.a2    = 0.5;
  filterpar.f1    = -1.0;
  filterpar.a1    = -1.0;

  /* High pass AS_Q signal */
  LALButterworthREAL8TimeSeries(&status,&GV.AS_Q,&filterpar);
  TESTSTATUS( &status );  

  /* If this chunk of size T is not the first within a science segment copy half of slugde to start of AS_Q */
  if(j>0){  
    for (n=0; n<SLUDGEFACTOR/2;n++) { 
      GV.AS_Q.data->data[n]=GV.hipasssludge.data->data[n]; 
    }  
  }  
  /* copy sludge at the end of residual signal to hipasssludge */
  for (n=0; n<SLUDGEFACTOR;n++) { 
    GV.hipasssludge.data->data[n]=GV.AS_Q.data->data[n+GV.AS_Q.data->length-SLUDGEFACTOR];
  }

  return 0;
}

/*******************************************************************************/

int HighPass2(int j)
{
  int n;
  PassBandParamStruc filterpar;

  filterpar.name  = "Butterworth High Pass";
  filterpar.nMax  = 10;
  filterpar.f2    = 40.0;
  filterpar.a2    = 0.5;
  filterpar.f1    = -1.0;
  filterpar.a1    = -1.0;

  /* High pass residual signal */
  LALButterworthREAL8TimeSeries(&status,&GV.hR,&filterpar);
  TESTSTATUS( &status );  

  /* If this chunk of size T is not the first within a science segment copy half of slugde to start of AS_Q */
  if(j>0){  
    for (n=0; n<SLUDGEFACTOR/2;n++) { 
      GV.hR.data->data[n]=GV.hipasssludgeR.data->data[n]; 
    }  
  }  
  /* copy sludge at the end of residual signal to hipasssludge */
  for (n=0; n<SLUDGEFACTOR;n++) { 
    GV.hipasssludgeR.data->data[n]=GV.hR.data->data[n+GV.hR.data->length-SLUDGEFACTOR];
  }

  /* High pass control signal */
  LALButterworthREAL8TimeSeries(&status,&GV.hC,&filterpar);
  TESTSTATUS( &status );  

  /* If this chunk of size T is not the first within a science segment copy half of slugde to start of AS_Q */
  if(j>0){  
    for (n=0; n<SLUDGEFACTOR/2;n++) { 
      GV.hC.data->data[n]=GV.hipasssludgeC.data->data[n]; 
    }  
  }  
  /* copy sludge at the end of residual signal to hipasssludge */
  for (n=0; n<SLUDGEFACTOR;n++) { 
    GV.hipasssludgeC.data->data[n]=GV.hC.data->data[n+GV.hC.data->length-SLUDGEFACTOR];
  }

  return 0;
}

/*******************************************************************************/

int DownsamplehR(void)
{
  int n;

  for (n=0; n<(int)GV.hR.data->length; n++) {
    GV.hR.data->data[n]=GV.uphR.data->data[n*USR];
  }

  return 0;
}

/*******************************************************************************/

int FilterSeries(MyIIRFilter *F, REAL8TimeSeries *TSeries, int sludge)
{
  int n,r;
  REAL8 yn,xn,xsum,ysum;
  MyIIRFilter H;  


  /* It's very important here to only filter up to the end of T: Do not filter sludge yet! 
   History gets messed up */

  for (n=0; n<(int)TSeries->data->length-sludge;n++) {

    xsum=0.0;
    ysum=0.0;

    xn=TSeries->data->data[n];
    
    for(r=0;r<F->xOrder-1;r++){
      xsum += F->xhist[r]*F->b[r+1];
    }
    xsum=xsum+xn*F->b[0];
    
    for(r=0;r<F->yOrder-1;r++){
      ysum -= F->yhist[r]*F->a[r+1];
    }
    
    yn=xsum+ysum;

    TSeries->data->data[n]=yn;

    for(r=F->xOrder-2;r>0;r--){
      F->xhist[r]=F->xhist[r-1];
    }
    for(r=F->yOrder-2;r>0;r--){
      F->yhist[r]=F->yhist[r-1];
    }

    F->yhist[0]=yn;
    F->xhist[0]=xn;
  }

  /* Filter sludge here; we filter it with a separate filter H so history of F does not get messed up */
  H=*F;
  for (n=TSeries->data->length-sludge; n<(int)TSeries->data->length;n++) 
    {
      xsum=0.0;
      ysum=0.0;

      xn=TSeries->data->data[n];
    
      for(r=0;r<H.xOrder-1;r++){
	xsum += H.xhist[r]*H.b[r+1];
      }
      xsum=xsum+xn*H.b[0];
      
      for(r=0;r<H.yOrder-1;r++){
	ysum -= H.yhist[r]*H.a[r+1];
      }
    
      yn=xsum+ysum;

      TSeries->data->data[n]=yn;

      for(r=H.xOrder-2;r>0;r--){
	H.xhist[r]=H.xhist[r-1];
      }
      for(r=H.yOrder-2;r>0;r--){
	H.yhist[r]=H.yhist[r-1];
      }

      H.yhist[0]=yn;
      H.xhist[0]=xn;
  }

  return 0;
}

/*******************************************************************************/

int UpsamplehR(int up_factor)
{
  int n;

  /* Set all values to 0 */
  for (n=0; n<(int)GV.uphR.data->length; n++) {
    GV.uphR.data->data[n] = 0.0;
  }

  /* Set one in every USR to the value of hR x USR */
  for (n=0; n<(int)GV.hR.data->length; n++) {
    GV.uphR.data->data[n * up_factor] = up_factor * GV.hR.data->data[n];
  }

  return 0;
}

/*******************************************************************************/

int hROverAlpha(int i)
{
  int n;
  double time,InterpolatedAlpha;

  gsl_interp_accel *acc_alpha = gsl_interp_accel_alloc();      /* GSL spline interpolation stuff */
  gsl_spline *spline_alpha = gsl_spline_alloc(gsl_interp_cspline,(UINT4)GV.duration/To);
  gsl_spline_init(spline_alpha,GV.ta_interp,GV.alpha,(UINT4)GV.duration/To);

  time=GV.gpsepoch.gpsSeconds-GV.SL[i].gpsstart;    /* time variable */

  for (n = 0; n < (int)GV.hR.data->length; n++) {

    InterpolatedAlpha=gsl_spline_eval(spline_alpha,time,acc_alpha);
    
    GV.hR.data->data[n] /= InterpolatedAlpha;
    time=time+GV.hR.deltaT;

  }

  /* clean up GSL spline interpolation stuff */
  gsl_spline_free(spline_alpha);
  gsl_interp_accel_free(acc_alpha);

  return 0;
}

/*******************************************************************************/

int hCTimesBeta(int i)
{
  int n;
  REAL8 time,InterpolatedBeta;

  gsl_interp_accel *acc_beta = gsl_interp_accel_alloc();      /* GSL spline interpolation stuff */
  gsl_spline *spline_beta = gsl_spline_alloc(gsl_interp_cspline,(UINT4)GV.duration/To);
  gsl_spline_init(spline_beta,GV.ta_interp,GV.beta,(UINT4)GV.duration/To);

  time=GV.gpsepoch.gpsSeconds-GV.SL[i].gpsstart;  /* time variable shifted by alphawings */

  for (n = 0; n < (int)GV.hC.data->length; n++) {

    InterpolatedBeta=gsl_spline_eval(spline_beta,time,acc_beta);
    
    GV.hC.data->data[n] *= InterpolatedBeta;
    time=time+GV.hC.deltaT;
  }

  /* clean up GSL spline interpolation stuff */
  gsl_spline_free(spline_beta);
  gsl_interp_accel_free(acc_beta);

  return 0;
}

/*******************************************************************************/

int ReadDTsecondsofData(struct CommandLineArgsTag CLA, int DT)
{
  static FrChanIn chanin_asq;
  static REAL4TimeSeries lAS_Q;        /* local REAL4 timeseries */
  int p;

  /* Allocate space for data vectors */
  LALCreateVector(&status,&lAS_Q.data,(UINT4)(DT/GV.AS_Q.deltaT +0.5)+SLUDGEFACTOR);
  TESTSTATUS( &status );

  chanin_asq.type  = ADCDataChannel;
  chanin_asq.name  = CLA.asq_chan;

  LALFrSeek(&status,&GV.gpsepoch,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&lAS_Q,&chanin_asq,framestream);
  TESTSTATUS( &status );

  /* Store AS_Q as double */  
  for (p=0; p<(int)GV.AS_Q.data->length; p++) {
    GV.AS_Q.data->data[p]=lAS_Q.data->data[p];
  }

  LALDestroyVector(&status,&lAS_Q.data);
  TESTSTATUS( &status );

  return 0;
}

/*******************************************************************************/

int AllocateData(struct CommandLineArgsTag CLA, int DT)
{
  static FrChanIn chanin_asq;
  static REAL4TimeSeries lAS_Q;       /* local REAL4 timeseries */

  lAS_Q.data=NULL;

  chanin_asq.type  = ADCDataChannel;
  chanin_asq.name  = CLA.asq_chan;

  /* Get channel time step size by calling LALFrGetREAL4TimeSeries */
  LALFrSeek(&status,&GV.gpsepoch,framestream);
  TESTSTATUS( &status );
  LALFrGetREAL4TimeSeries(&status,&lAS_Q,&chanin_asq,framestream);
  TESTSTATUS( &status );

  GV.AS_Q.deltaT=lAS_Q.deltaT;
  LALDCreateVector(&status,&GV.AS_Q.data,(UINT4)(DT/GV.AS_Q.deltaT +0.5)+SLUDGEFACTOR);
  TESTSTATUS( &status );

  GV.hR.deltaT=lAS_Q.deltaT;
  GV.hC.deltaT=lAS_Q.deltaT;
  LALDCreateVector(&status,&GV.hR.data,(UINT4)(DT/GV.hR.deltaT +0.5)+SLUDGEFACTOR);
  TESTSTATUS( &status );
  LALDCreateVector(&status,&GV.hC.data,(UINT4)(DT/GV.hC.deltaT +0.5)+SLUDGEFACTOR);
  TESTSTATUS( &status );

  GV.hCx.deltaT=lAS_Q.deltaT;
  LALDCreateVector(&status,&GV.hCx.data,(UINT4)(DT/GV.hCx.deltaT +0.5)+SLUDGEFACTOR);
  TESTSTATUS( &status );
  GV.hCy.deltaT=lAS_Q.deltaT;
  LALDCreateVector(&status,&GV.hCy.data,(UINT4)(DT/GV.hCy.deltaT +0.5)+SLUDGEFACTOR);
  TESTSTATUS( &status );

  GV.h.deltaT=lAS_Q.deltaT;
  LALDCreateVector(&status,&GV.h.data,(UINT4)(DT/GV.h.deltaT +0.5));
  TESTSTATUS( &status );

  GV.uphR.deltaT=GV.hR.deltaT/USR;
  LALDCreateVector(&status,&GV.uphR.data,(UINT4)(USR*DT/GV.hR.deltaT +0.5)+USR*SLUDGEFACTOR);
  TESTSTATUS( &status );

  return 0;
}

/*******************************************************************************/

int MakeFilters(void)
{
  int l,n;
 
  GV.Cinv.yOrder=CinvRecursOrder;
  GV.Cinv.xOrder=CinvDirectOrder;

  /* Fill coefficient vectors with coefficients from filters.h */
  for(l=0;l<GV.Cinv.xOrder;l++) GV.Cinv.b[l]=CinvDirectCoefs[l];
  for(l=0;l<GV.Cinv.yOrder;l++) GV.Cinv.a[l]=CinvRecursCoefs[l];
  for(l=0;l<GV.Cinv.yOrder-1;l++) GV.Cinv.yhist[l]=0.0;
  for(l=0;l<GV.Cinv.xOrder-1;l++) GV.Cinv.xhist[l]=0.0;

  for(n=0;n<NGfilt;n++){
    GV.G[n].yOrder=G_Dord;
    GV.G[n].xOrder=G_Rord;
    
    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<GV.G[n].xOrder;l++) GV.G[n].b[l]=G_D[n][l];
    for(l=0;l<GV.G[n].yOrder;l++) GV.G[n].a[l]=G_R[n][l];
    for(l=0;l<GV.G[n].yOrder-1;l++) GV.G[n].yhist[l]=0.0;
    for(l=0;l<GV.G[n].xOrder-1;l++) GV.G[n].xhist[l]=0.0;
  }

  GV.AA.yOrder= A_0_Rord;
  GV.AA.xOrder= A_0_Dord;

  /* Fill coefficient vectors with coefficients from filters.h */
  for(l=0;l<GV.AA.xOrder;l++) GV.AA.b[l]=A_0_D[l];
  for(l=0;l<GV.AA.yOrder;l++) GV.AA.a[l]=A_0_R[l];
  for(l=0;l<GV.AA.yOrder-1;l++) GV.AA.yhist[l]=0.0;
  for(l=0;l<GV.AA.xOrder-1;l++) GV.AA.xhist[l]=0.0;

  GV.ADW.yOrder= A_DW_Ord;
  GV.ADW.xOrder= A_DW_Ord;

  /* Fill coefficient vectors with coefficients from filters.h */
  for(l=0;l<GV.ADW.xOrder;l++) GV.ADW.b[l]=A_DW_D[l];
  for(l=0;l<GV.ADW.yOrder;l++) GV.ADW.a[l]=A_DW_R[l];
  for(l=0;l<GV.ADW.yOrder-1;l++) GV.ADW.yhist[l]=0.0;
  for(l=0;l<GV.ADW.xOrder-1;l++) GV.ADW.xhist[l]=0.0;

  GV.AADW.yOrder= A_ADW_Ord;
  GV.AADW.xOrder= A_ADW_Ord;
  
  /* Fill coefficient vectors with coefficients from filters.h */
  for(l=0;l<GV.AADW.xOrder;l++) GV.AADW.b[l]=A_ADW_D[l];
  for(l=0;l<GV.AADW.yOrder;l++) GV.AADW.a[l]=A_ADW_R[l];
  for(l=0;l<GV.AADW.yOrder-1;l++) GV.AADW.yhist[l]=0.0;
  for(l=0;l<GV.AADW.xOrder-1;l++) GV.AADW.xhist[l]=0.0;
    
  for(n=0;n<NAXfilt;n++){
    GV.AX[n].yOrder=A_digital_Rord;
    GV.AX[n].xOrder=A_digital_Dord;
    
    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<GV.AX[n].xOrder;l++) GV.AX[n].b[l]=AX_D[n][l];
    for(l=0;l<GV.AX[n].yOrder;l++) GV.AX[n].a[l]=AX_R[n][l];
    for(l=0;l<GV.AX[n].yOrder-1;l++) GV.AX[n].yhist[l]=0.0;
    for(l=0;l<GV.AX[n].xOrder-1;l++) GV.AX[n].xhist[l]=0.0;
  }

  for(n=0;n<NAYfilt;n++){
    GV.AY[n].yOrder=A_digital_Rord;
    GV.AY[n].xOrder=A_digital_Dord;
    
    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<GV.AY[n].xOrder;l++) GV.AY[n].b[l]=AY_D[n][l];
    for(l=0;l<GV.AY[n].yOrder;l++) GV.AY[n].a[l]=AY_R[n][l];
    for(l=0;l<GV.AY[n].yOrder-1;l++) GV.AY[n].yhist[l]=0.0;
    for(l=0;l<GV.AY[n].xOrder-1;l++) GV.AY[n].xhist[l]=0.0;
  }

  return 0;
}

/*******************************************************************************/

int GetFactors(struct CommandLineArgsTag CLA)
{

FrPos pos1;

static REAL4TimeSeries darm;
static REAL4TimeSeries asq;
static REAL4TimeSeries exc;

static FrChanIn chanin_darm;
static FrChanIn chanin_asq;
static FrChanIn chanin_exc;

CalFactors factors;
UpdateFactorsParams params;

REAL4Vector *asqwin=NULL,*excwin=NULL,*darmwin=NULL;  /* windows */

LALWindowParams winparams;

INT4 k,m,outflag=0;
LIGOTimeGPS localgpsepoch=GV.gpsepoch; /* Local variable epoch used to calculate the calibration factors */
 
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
  LALCreateVector(&status,&asq.data,(UINT4)(To/asq.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&darm.data,(UINT4)(To/darm.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&exc.data,(UINT4)(To/exc.deltaT +0.5));
  TESTSTATUS( &status );

  /* Create Window vectors */
  LALCreateVector(&status,&asqwin,(UINT4)(To/asq.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&darmwin,(UINT4)(To/darm.deltaT +0.5));
  TESTSTATUS( &status );
  LALCreateVector(&status,&excwin,(UINT4)(To/exc.deltaT +0.5));
  TESTSTATUS( &status );

  winparams.type=Hann;
   
  /* windows for time domain channels */
  /* asq */
  winparams.length=(INT4)(To/asq.deltaT +0.5);
  LALWindow(&status,asqwin,&winparams);
  TESTSTATUS( &status );
  
  /* darm */
  winparams.length=(INT4)(To/darm.deltaT +0.5);
  LALWindow(&status,darmwin,&winparams);
  TESTSTATUS( &status );

  /* exc */
  winparams.length=(INT4)(To/exc.deltaT +0.5);
  LALWindow(&status,excwin,&winparams);
  TESTSTATUS( &status );

  /* If we have a user input factors file then open it */
  if(CLA.alphafile != NULL)
    {
      fpAlpha=fopen(CLA.alphafile,"w");
      if (fpAlpha==NULL) 
	{
	  fprintf(stderr,"Could not open %s!\n",CLA.alphafile);
	  return 1;
	}
      setvbuf( fpAlpha, NULL, _IONBF, 0 );  /* This stops buffering -- From Duncan Brown*/
    }

  for(m=0;m < GV.duration/To;m++)
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
      for(k=0;k<(INT4)(To/asq.deltaT +0.5);k++)
	{
	  asq.data->data[k] *= 2.0*asqwin->data[k];
	}
      for(k=0;k<(INT4)(To/darm.deltaT +0.5);k++)
	{
	  darm.data->data[k] *= 2.0*darmwin->data[k];
	}
      for(k=0;k<(INT4)(To/exc.deltaT +0.5);k++)
	{
	  exc.data->data[k] *= 2.0*excwin->data[k];
	}

      /* set params to call LALComputeCalibrationFactors */
      params.darmCtrl = &darm;
      params.asQ = &asq;
      params.exc = &exc;
      params.lineFrequency = CLA.f;
      params.openloop.re =  CLA.G0Re;
      params.openloop.im =  CLA.G0Im;
      params.digital.re = CLA.D0Re;
      params.digital.im = CLA.D0Im;

      LALComputeCalibrationFactors(&status,&factors,&params);
      TESTSTATUS( &status );

      GV.alpha[m]= factors.alpha.re;
      
      if( GV.alpha[m] < 0.3 ) 
	{
	 fprintf(stderr,"There were invalid values of alpha at %d!\n",localgpsepoch.gpsSeconds);
	 GV.alpha[m]=1.0;
	 if (m>0) GV.alpha[m]=GV.alpha[m-1];
	}

      if( GV.alpha[m] > 2.0 ) 
	{
	 fprintf(stderr,"There were invalid values of alpha at %d!\n",localgpsepoch.gpsSeconds);
	 GV.alpha[m]=1.0;
	 if (m>0) GV.alpha[m]=GV.alpha[m-1];
	}
      
      GV.beta[m]= factors.beta.re;

      GV.ta_interp[m]=m*To;

      if(CLA.alphafile != NULL)
	{
	  fprintf(fpAlpha,"%d %d %e %e %e %e %e %e %e %e %e %e %e %e\n",m,localgpsepoch.gpsSeconds,factors.alpha.re,factors.alpha.im,
		  factors.beta.re,factors.beta.im,
		  factors.alphabeta.re,factors.alphabeta.im,
		  factors.asq.re*2/To,factors.asq.im*2/To,
		  factors.darm.re*2/To,factors.darm.im*2/To,
		  factors.exc.re*2/To,factors.exc.im*2/To);
	}

      localgpsepoch.gpsSeconds = localgpsepoch.gpsSeconds+To;
    }

  /* Clean up */
  LALDestroyVector(&status,&darm.data);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&exc.data);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&asq.data);
  TESTSTATUS( &status );

  LALDestroyVector(&status,&asqwin);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&darmwin);
  TESTSTATUS( &status );
  LALDestroyVector(&status,&excwin);
  TESTSTATUS( &status );

  if(CLA.alphafile != NULL)
    {
      fclose(fpAlpha);
    }


  if(outflag == 1) fprintf(stderr,"There were invalid values of alpha at %d! \n", GV.gpsepoch.gpsSeconds);

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
  CLA->G0Re=0.0;
  CLA->G0Im=0.0;
  CLA->D0Re=0.0;
  CLA->D0Im=0.0;
  CLA->FrCacheFile=NULL;
  CLA->SegmentsFile=NULL;
  CLA->exc_chan=NULL;
  CLA->darm_chan=NULL;
  CLA->asq_chan=NULL;
  CLA->alphafile=NULL;

  /* Scan through list of command line arguments */
  while (!errflg && ((c = getopt(argc, argv,"hf:F:S:A:E:D:b:i:j:k:l:"))!=-1))
    switch (c) {
    case 'f':
      /* calibration line frequency */
      CLA->f=atof(optarg);
      break;
    case 'F':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      break;
    case 'S':
      /* name of segments file */
      CLA->SegmentsFile=optarg;
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
   case 'h':
      /* print usage/help message */
      fprintf(stdout,"All arguments are required except -b. They are:\n");
      fprintf(stdout,"\t-f\tFLOAT\t Calibration line frequency in Hz.\n");
      fprintf(stdout,"\t-i\tFLOAT\t Real part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\t-j\tFLOAT\t Imaginary part of the open loop gain at the calibration line frequency.\n");
      fprintf(stdout,"\t-k\tFLOAT\t Real part of the digital filter at the calibration line frequency.\n");
      fprintf(stdout,"\t-l\tFLOAT\t Imaginary part of digital filter at the calibration line frequency.\n");
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
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
  if(CLA->FrCacheFile == NULL)
    {
      fprintf(stderr,"No frame cache file specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
  if(CLA->SegmentsFile == NULL)
    {
      fprintf(stderr,"No segments file specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
   if(CLA->exc_chan == NULL)
    {
      fprintf(stderr,"No excitation channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
   if(CLA->darm_chan == NULL)
    {
      fprintf(stderr,"No darm channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
   if(CLA->asq_chan == NULL)
    {
      fprintf(stderr,"No asq channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
  if(CLA->G0Re == 0.0 )
    {
      fprintf(stderr,"No real part of open loop gain specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
  if(CLA->G0Im == 0.0 )
    {
      fprintf(stderr,"No imaginary part of open loop gain specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
  if(CLA->D0Re == 0.0 )
    {
      fprintf(stderr,"No real part of digital filter specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
  if(CLA->D0Im == 0.0 )
    {
      fprintf(stderr,"No imaginary part of digital filter specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      

  return errflg;
}

/*******************************************************************************/

int FreeMem(void)
{

  LALFrClose(&status,&framestream);
  TESTSTATUS( &status );

  LALDDestroyVector(&status,&GV.hipasssludge.data);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&GV.hipasssludgeR.data);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&GV.hipasssludgeC.data);
  TESTSTATUS( &status );
  LALDDestroyVector(&status,&GV.lopasssludge.data);
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
 
  return 0;
}

/*******************************************************************************/
#endif
