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
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/FileIO.h>
#include <lal/AVFactories.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <math.h>
#include <lal/Window.h>
#include <lal/Calibration.h>
#include <lal/LALConstants.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <math.h>
#include <lal/BandPassTimeSeries.h>

#include "filters-H1-S2.h"                 /* files that contains the filter coefficients */

#define MAXLINERS 76800                   /* Max lines read in Freq Response files */
#define MAXLINESEGS 10000                 /* Maximum number of science segments */
#define To 60                             /* length of time used in table of alphas and betas (seconds)*/
#define T  1                              /* length of time calibrated per iteration (seconds) */
#define SR 16384                          /* Sampling rate of data (Hz) */
#define USR 16                            /* Upsampling factor */
#define MAXALPHAS 4000                    /* Maximum number of calibration factors in a science segment */
#define ALPHAWINGS 4                      /* Number of extra values of alpha at beginning and end of alpha array */

#define SLUDGEFACTOR 2048                 /* Number of samples in AS_Q to keep as extra sludge for smooth filtering */
 
/***************************************************************************/
/* Complex division routine -- used to calculate beta from alpha and alpha*beta */
static COMPLEX16 *cdiv( COMPLEX16 *pc, COMPLEX16 *pa, COMPLEX16 *pb )
{
  COMPLEX16 a = *pa;
  COMPLEX16 b = *pb;
  COMPLEX16 c;
  REAL8 rat;
  REAL8 den;
  if ( fabs( b.re ) > fabs( b.im ) )
  {
    rat = b.im / b.re;
    den = b.re + rat * b.im;
    c.re = ( a.re + rat * a.im ) / den;
    c.im = ( a.im - rat * a.re ) / den;
  }
  else
  {
    rat = b.re / b.im;
    den = b.im + rat * b.re;
    c.re = ( a.re * rat + a.im ) / den;
    c.im = ( a.im * rat - a.re ) / den;
  }
  *pc = c;
  return pc;
}

/***************************************************************************/

/* STRUCTURES */
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
  INT4 gpsstart;          /* GPS Seconds of start of segment group */
  INT4 gpsend;            /* number of segments starting at tgps */
  INT4 seglength;         /* length of segment in seconds */       
} SegmentList;

typedef
struct MyIIRFilter {
  INT4 yOrder;
  INT4 xOrder;
  REAL8 a[20];
  REAL8 b[20];
  REAL8 yhist[20];
  REAL8 xhist[20];
} MyIIRFilter;


/***************************************************************************/

/* GLOBAL VARIABLES */
static LALStatus status;
INT4 lalDebugLevel=3;

FrCache *framecache;                                       /* frame reading variables */
FrStream *framestream=NULL;

COMPLEX16 Rf0,Cf0,Af0;                                     /* Response, sensing and actuation function values at the 
							      frequency of the calibration line */

SegmentList SL[MAXLINESEGS];                               /* Structure containing science segement info */
INT4 numsegs;                                              /* number of science segments */

REAL8 alpha[MAXALPHAS],beta[MAXALPHAS];                    /* Arrays that contain the calibration factors computed in To time intervals */

LIGOTimeGPS gpsepoch;                                      /* Global variables epoch and duration used to calculate the calibration factors */
INT4 duration;                                             /* they are set every science segment to GPS start time and duration of segment */

MyIIRFilter Cinv,G[NGfilt],AA[2],AX[NAXfilt],AY[NAYfilt];     /* Inverse sensing, servo, analog part of actuation, digital x-arm actuation digital y-arm actuation */

REAL8TimeSeries AS_Q;                                      /* Global time series with AS_Q as a double */

REAL8TimeSeries hR,hC,uphR,hCx,hCy,h;                      /* Residual, control and total strains; up... time series are upsampled */

REAL8TimeSeries hipasssludge,lopasssludge;                 /* Sludge for high and low pass filtering of data */

/***************************************************************************/

/* FUNCTION PROTOTYPES */
/* Reads the command line */
int ReadCommandLine(int argc,char *argv[],struct CommandLineArgsTag *CLA);

/* Reads time segment, cache, response,sensing and actuation files */
int ReadFiles(struct CommandLineArgsTag CLA);             

/* Create the filters from the filters.h file */
int MakeFilters(struct CommandLineArgsTag CLA);

/* Allocates space for data time series AS_Q and DARM_CTRL */
int AllocateData(struct CommandLineArgsTag CLA);

/* Produces a table of calibration factors for a science segment; 
the factors are calculated every interval To*/
int GetFactors(struct CommandLineArgsTag CLA);   

/* Reads T seconds of AS_Q and DARM_CTRL data */
int ReadTsecondsofData(struct CommandLineArgsTag CLA);

/* High passes data doing all the sludge accouting */
int HighPassAS_Q(int j);

/* Multiplies each sample of AS_Q by 1/alpha */
int hROverAlpha(int i);                                        

/* Multiplies each sample of AS_Q by beta */
int hCTimesBeta(int i);

/* Upsamples AS_QR storing result in upAS_Q by putting SR-1 zeros between samples in AS_Q */
int UpsamplehR();   

/* Low passes upsampled AS_Q doing all the sludge accouting */
int LowPasshR(int j);

/* My own filtering function that sets the history correctly */
int FilterSeries(MyIIRFilter *F, REAL8TimeSeries *TSeries);

/* Downsamples calibarted AS_Q storing result in hr by takine 1 sample out of USR  zeros between samples in AS_Q */
int DownsamplehR(void);   

/* Frees the memory */
int FreeMem(void);                                        


/************************************* MAIN PROGRAM *************************************/

int main(int argc,char *argv[])
{
  int i,j,p; 

  if (ReadCommandLine(argc,argv,&CommandLineArgs)) return 1;
  if (ReadFiles(CommandLineArgs)) return 3;
  if (MakeFilters(CommandLineArgs)) return 4;

  gpsepoch.gpsSeconds=SL[0].gpsstart; /* Set global variable epoch */
  gpsepoch.gpsNanoSeconds=0;
  if (AllocateData(CommandLineArgs)) return 4;

  for(i=0;i<numsegs;i++)
    {
      if(SL[i].seglength < To)
	{
	  fprintf(stdout,"Throwing out segment #%d, it is shorter than %d seconds\n",i,To);
	  break;
	}

      gpsepoch.gpsSeconds=SL[i].gpsstart; /* Set global variable epoch */
      gpsepoch.gpsNanoSeconds=0;
      duration=SL[i].seglength;

      if(GetFactors(CommandLineArgs)) return 5;
      for(j=0; j<30/*duration/T*/; j++) 
	{

	  /* Reads T seconds of AS_Q data */
	  if (ReadTsecondsofData(CommandLineArgs)) return 4;

	  /* hipass with Butterworth filter at 40Hz */
	  if(HighPassAS_Q(j)) return 4;

	  /* Copy AS_Q into residual and control signal time series */  
	  for (p=0; p<hR.data->length; p++) {
	    hR.data->data[p]=AS_Q.data->data[p];
	  }
	  for (p=0; p<hC.data->length; p++) {
	    hC.data->data[p]=AS_Q.data->data[p];
	  }


	  /******** COMPUTE RESIDUAL STRAIN **********/

	  /* multiply hR by 1/alpha(t) */
	  if (hROverAlpha(i)) return 4;

	  /* Upsample hR by a factor of USR defined above */
	  if (UpsamplehR()) return 4;

	  /* smooth with low pass Butterworth filter at 7000Hz */
	  if(LowPasshR(j)) return 4;

          /* Filter through inverse of sensing function */
	  if (FilterSeries(&Cinv,&hR)) return 5;

	  /* Downsample by a factor of 16 */
	  if(DownsamplehR()) return 3;


	  /******** COMPUTE CONTROL STRAIN **********/

	  /* multiply hC by beta(t) */
	  if (hCTimesBeta(i)) return 4; 

          /* Filter hC through servo to get DARM_CTRL */
	  for(p=NGfilt-1;p>=0;p--){
	  if (FilterSeries(&G[p],&hC)) return 5;
	  }

	  /* Adjust to account for servo gain */ 
	  for (p=0; p<hC.data->length;p++) {
	    hC.data->data[p]= ServoGain*hC.data->data[p];
 	  }

	  /* Copy data into x and y time series for parallel filtering */
	  for (p=0; p<hCx.data->length; p++) {
	    hCx.data->data[p]=hC.data->data[p];
	  }
	  for (p=0; p<hCy.data->length; p++) {
	    hCy.data->data[p]=hC.data->data[p];
	  }

	  /* Filter x-arm */
	  for(p=NAXfilt-1;p>=0;p--){
	  if (FilterSeries(&AX[p],&hCx)) return 5;
	  }

	  /* Filter y-arm */
	  for(p=NAYfilt-1;p>=0;p--){
	  if (FilterSeries(&AY[p],&hCy)) return 5;
	  }

	  /* add them together */
	  for (p=0; p<hC.data->length; p++) {
	    hC.data->data[p]=(hCx.data->data[p]+hCy.data->data[p])/2;
	  }

	  for (p=0; p<hC.data->length; p++) {
	    hC.data->data[p]=0.0;
	  }
	  if(j==0) hC.data->data[0]=1.0;
	  
 	  /* filter through analog part of actuation */ 

	  if (FilterSeries(&AA[0],&hC)) return 5;
	  if (FilterSeries(&AA[1],&hC)) return 5;
	  if (FilterSeries(&AA[2],&hC)) return 5;

	  for (p=0; p<hC.data->length-SLUDGEFACTOR;p++) {
	     fprintf(stdout,"%e  %e\n",j*T+p*hC.deltaT,hC.data->data[p]);
 	  } 

	  /******** COMPUTE NET CALIBRATED STRAIN **********/

	  for (p=0; p<h.data->length; p++) {
	    h.data->data[p]= hR.data->data[p]+ hC.data->data[p];
	  }
	  
/* 	  for (p=0; p<h.data->length-SLUDGEFACTOR;p++) { */
/* 	     fprintf(stdout,"%e  %e\n",j*T+p*h.deltaT,h.data->data[p]); */
/*  	  }  */


	  gpsepoch.gpsSeconds += T;

	}

      //Delete filter histories

    }
  
  if(FreeMem()) return 8;

  return 0;
}

/************************************* MAIN PROGRAM ENDS *************************************/


/*  FUNCTIONS */

/*******************************************************************************/

int LowPasshR(int j)
{
  int n;
  PassBandParamStruc filterpar;

  filterpar.name  = "Butterworth Low Pass at 7 kHz";
  filterpar.nMax  = 5;
  filterpar.f2    = -1.0;
  filterpar.a2    = -1.0;
  filterpar.f1    = 7000.0;
  filterpar.a1    = 0.5;

  LALButterworthREAL8TimeSeries(&status,&uphR,&filterpar);
  
  /* If this chunk of size T is not the first within a science segment copy half of slugde to start of hR */
  if(j>0){  
    for (n=0; n<SLUDGEFACTOR/2;n++) { 
      uphR.data->data[n]=lopasssludge.data->data[n]; 
    }  
  }  
  /* copy sludge at the end of AS_Q to hipasssludge */
  for (n=0; n<SLUDGEFACTOR;n++) { 
    lopasssludge.data->data[n]=uphR.data->data[n+uphR.data->length-SLUDGEFACTOR];
  }

  return 0;
}

/*******************************************************************************/

int HighPassAS_Q(int j)
{
  int n;
  PassBandParamStruc filterpar;

  filterpar.name  = "Butterworth High Pass at 40 Hz";
  filterpar.nMax  = 5;
  filterpar.f2    = 40.0;
  filterpar.a2    = 0.5;
  filterpar.f1    = -1.0;
  filterpar.a1    = -1.0;

  LALButterworthREAL8TimeSeries(&status,&AS_Q,&filterpar);
  
  /* If this chunk of size T is not the first within a science segment copy half of slugde to start of AS_Q */
  if(j>0){  
    for (n=0; n<SLUDGEFACTOR/2;n++) { 
      AS_Q.data->data[n]=hipasssludge.data->data[n]; 
    }  
  }  
  /* copy sludge at the end of AS_Q to hipasssludge */
  for (n=0; n<SLUDGEFACTOR;n++) { 
    hipasssludge.data->data[n]=AS_Q.data->data[n+AS_Q.data->length-SLUDGEFACTOR];
  }

  return 0;
}

/*******************************************************************************/

int DownsamplehR(void)
{
  int n;

  for (n=0; n<hR.data->length-SLUDGEFACTOR; n++) {
    hR.data->data[n]=uphR.data->data[n*USR]/SensingGain;
  }

  return 0;
}
/*******************************************************************************/

int FilterSeries(MyIIRFilter *F, REAL8TimeSeries *TSeries)
{
  int n,r;
  REAL8 yn,xn,xsum,ysum;


  /* It's very important here to only filter up to the end of T: Do not filter sludge! 
   History gets messed up */

  for (n=0; n<TSeries->data->length-SLUDGEFACTOR;n++) {

    xsum=0.0;
    ysum=0.0;

    xn=TSeries->data->data[n];
    
    for(r=0;r<F->xOrder-1;r++){
      xsum += F->xhist[r]*F->b[r+1];
    }
    xsum=xsum+xn*F->b[0];
    
    for(r=0;r<F->yOrder-1;r++){
      ysum += F->yhist[r]*F->a[r+1];
    }
    
    yn=xsum+ysum;

    TSeries->data->data[n]=yn;

    for(r=1;r<F->xOrder-1;r++){
      F->xhist[r]=F->xhist[r-1];
    }
    for(r=1;r<F->yOrder-1;r++){
      F->yhist[r]=F->yhist[r-1];
    }

    F->yhist[0]=yn;
    F->xhist[0]=xn;
  }
  return 0;
}

/*******************************************************************************/

int UpsamplehR(void)
{
  int n;

  /* Set all values to 0 */
  for (n=0; n<uphR.data->length; n++) {
    uphR.data->data[n] = 0.0;
  }

  /* Set one in every USR to the value of hR x USR */
  for (n=0; n<hR.data->length-(SLUDGEFACTOR-SLUDGEFACTOR/USR); n++) {
    uphR.data->data[n*USR] = USR * hR.data->data[n];
  }

  return 0;
}

/*******************************************************************************/

int hROverAlpha(int i)
{
  int n,k;
  REAL8 time,InterpolatedAlpha,d1,d2,dirichlet1,dirichlet2;
  int r,r1,r2;                                /* r is the index corresponding to the alpha
						 previous or at the point at time time */

  time=gpsepoch.gpsSeconds-SL[i].gpsstart+ALPHAWINGS*To;  /* time variable shifted by alphawings */
  r=time/To;
  for (n = 0; n < hR.data->length; n++) {
    InterpolatedAlpha=0.0;
    
    d1=time/To-r;
    d2=time/To-(r+1) ;
    if (d1==0.0) {
      InterpolatedAlpha=alpha[r]  ;
    } else if (d2==0.0) {
      InterpolatedAlpha= alpha[r+1] ;
    } else {
      REAL8 norm=0;
      for(k=0; k < ALPHAWINGS+1; k++) {
	r1=r-k;  
	r2=(r+1)+k;
	d1=time/To-r1;
	d2=time/To-r2;
	
	dirichlet1=sin(LAL_PI*d1)/(LAL_PI*d1);
	dirichlet2=sin(LAL_PI*d2)/(LAL_PI*d2);
      
	InterpolatedAlpha +=  alpha[r1]*  dirichlet1+  alpha[r2]* dirichlet2; 
	norm+=dirichlet1+dirichlet2;
      }
      InterpolatedAlpha /= norm;
    }
    
    hR.data->data[n] /= InterpolatedAlpha;
    //fprintf(stdout,"%f %f\n",time,InterpolatedAlpha);
    time=time+hR.deltaT;
  }
  return 0;
}
/*******************************************************************************/

int hCTimesBeta(int i)
{
  int n,k;
  REAL8 time,InterpolatedBeta,d1,d2,dirichlet1,dirichlet2;
  int r,r1,r2;                                /* r is the index corresponding to the alpha
						 previous or at the point at time time */

  time=gpsepoch.gpsSeconds-SL[i].gpsstart+ALPHAWINGS*To;  /* time variable shifted by alphawings */
  r=time/To;
  for (n = 0; n < hC.data->length; n++) {
    InterpolatedBeta=0.0;
    
    d1=time/To-r;
    d2=time/To-(r+1) ;
    if (d1==0.0) {
      InterpolatedBeta=beta[r]  ;
    } else if (d2==0.0) {
      InterpolatedBeta= beta[r+1] ;
    } else {
      REAL8 norm=0;
      for(k=0; k < ALPHAWINGS+1; k++) {
	r1=r-k;  
	r2=(r+1)+k;
	d1=time/To-r1;
	d2=time/To-r2;
	
	dirichlet1=sin(LAL_PI*d1)/(LAL_PI*d1);
	dirichlet2=sin(LAL_PI*d2)/(LAL_PI*d2);
      
	InterpolatedBeta +=  beta[r1]*  dirichlet1+  beta[r2]* dirichlet2; 
	norm+=dirichlet1+dirichlet2;
      }
      InterpolatedBeta /= norm;
    }
    
    hC.data->data[n] *= InterpolatedBeta;
    //    fprintf(stdout,"%f %f\n",time,InterpolatedBeta);
    time=time+hC.deltaT;
  }
  return 0;
}

/*******************************************************************************/

int ReadTsecondsofData(struct CommandLineArgsTag CLA)
{
  static FrChanIn chanin_asq;
  static FrChanIn chanin_darm;
  static REAL4TimeSeries lAS_Q;        /* local REAL4 timeseries */
  int p;

  /* Allocate space for data vectors */
  LALCreateVector(&status,&lAS_Q.data,(INT4)(T/AS_Q.deltaT +0.5)+SLUDGEFACTOR);

  chanin_asq.type  = ADCDataChannel;
  chanin_asq.name  = CLA.asq_chan;

  chanin_darm.type  = ADCDataChannel;
  chanin_darm.name  = CLA.darm_chan;

  //   fprintf(stdout,"%d\n",gpsepoch.gpsSeconds);

  LALFrSeek(&status,&gpsepoch,framestream);
  LALFrGetREAL4TimeSeries(&status,&lAS_Q,&chanin_asq,framestream);


  /* Store AS_QR and AS_QC as doubles */  
  for (p=0; p<AS_Q.data->length; p++) {
    AS_Q.data->data[p]=lAS_Q.data->data[p];
  }

  LALDestroyVector(&status,&lAS_Q.data);

  return 0;
}

/*******************************************************************************/

int AllocateData(struct CommandLineArgsTag CLA)
{
  static FrChanIn chanin_asq;
  static REAL4TimeSeries lAS_Q;       /* local REAL4 timeseries */

  lAS_Q.data=NULL;

  chanin_asq.type  = ADCDataChannel;
  chanin_asq.name  = CLA.asq_chan;

  /* Get channel time step size by calling LALFrGetREAL4TimeSeries */
  LALFrSeek(&status,&gpsepoch,framestream);
  LALFrGetREAL4TimeSeries(&status,&lAS_Q,&chanin_asq,framestream);

  AS_Q.deltaT=lAS_Q.deltaT;
  LALDCreateVector(&status,&AS_Q.data,(INT4)(T/AS_Q.deltaT +0.5)+SLUDGEFACTOR);

  hR.deltaT=lAS_Q.deltaT;
  hC.deltaT=lAS_Q.deltaT;
  LALDCreateVector(&status,&hR.data,(INT4)(T/hR.deltaT +0.5)+SLUDGEFACTOR);
  LALDCreateVector(&status,&hC.data,(INT4)(T/hC.deltaT +0.5)+SLUDGEFACTOR);

  hCx.deltaT=lAS_Q.deltaT;
  LALDCreateVector(&status,&hCx.data,(INT4)(T/hCx.deltaT +0.5)+SLUDGEFACTOR);
  hCy.deltaT=lAS_Q.deltaT;
  LALDCreateVector(&status,&hCy.data,(INT4)(T/hCy.deltaT +0.5)+SLUDGEFACTOR);

  h.deltaT=lAS_Q.deltaT;
  LALDCreateVector(&status,&h.data,(INT4)(T/h.deltaT +0.5)+SLUDGEFACTOR);

  uphR.deltaT=hR.deltaT/USR;
  LALDCreateVector(&status,&uphR.data,(INT4)(USR*T/hR.deltaT +0.5)+SLUDGEFACTOR);



  hipasssludge.deltaT=AS_Q.deltaT;
  LALDCreateVector(&status,&hipasssludge.data,SLUDGEFACTOR);

  lopasssludge.deltaT=AS_Q.deltaT;
  LALDCreateVector(&status,&lopasssludge.data,SLUDGEFACTOR);

  return 0;
}

/*******************************************************************************/

int MakeFilters(struct CommandLineArgsTag CLA)
{
  int l,n;
 
  Cinv.yOrder=CinvRecursOrder;
  Cinv.xOrder=CinvDirectOrder;

  /* Fill coefficient vectors with coefficients from filters.h */
  for(l=0;l<Cinv.xOrder;l++) Cinv.b[l]=CinvDirectCoefs[l];
  for(l=0;l<Cinv.yOrder;l++) Cinv.a[l]=CinvRecursCoefs[l];
  for(l=0;l<Cinv.yOrder-1;l++) Cinv.yhist[l]=0.0;
  for(l=0;l<Cinv.xOrder-1;l++) Cinv.xhist[l]=0.0;

  for(n=0;n<NGfilt;n++){
    G[n].yOrder=G_Dord;
    G[n].xOrder=G_Rord;
    
    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<G[n].xOrder;l++) G[n].b[l]=G_D[n][l];
    for(l=0;l<G[n].yOrder;l++) G[n].a[l]=G_R[n][l];
    for(l=0;l<G[n].yOrder-1;l++) G[n].yhist[l]=0.0;
    for(l=0;l<G[n].xOrder-1;l++) G[n].xhist[l]=0.0;
  }


  for(n=0;n<NAAfilt;n++){
    AA[n].yOrder= A_0_Rord;
    AA[n].xOrder= A_0_Dord;

    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<AA[n].xOrder;l++) AA[n].b[l]=A_0_D[n][l];
    for(l=0;l<AA[n].yOrder;l++) AA[n].a[l]=A_0_R[n][l];
    for(l=0;l<AA[n].yOrder-1;l++) AA[n].yhist[l]=0.0;
    for(l=0;l<AA[n].xOrder-1;l++) AA[n].xhist[l]=0.0;
  }

  for(n=0;n<NAXfilt;n++){
    AX[n].yOrder=A_digital_Rord;
    AX[n].xOrder=A_digital_Dord;
    
    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<AX[n].xOrder;l++) AX[n].b[l]=AX_D[n][l];
    for(l=0;l<AX[n].yOrder;l++) AX[n].a[l]=AX_R[n][l];
    for(l=0;l<AX[n].yOrder-1;l++) AX[n].yhist[l]=0.0;
    for(l=0;l<AX[n].xOrder-1;l++) AX[n].xhist[l]=0.0;
  }

  for(n=0;n<NAYfilt;n++){
    AY[n].yOrder=A_digital_Rord;
    AY[n].xOrder=A_digital_Dord;
    
    /* Fill coefficient vectors with coefficients from filters.h */
    for(l=0;l<AY[n].xOrder;l++) AY[n].b[l]=AY_D[n][l];
    for(l=0;l<AY[n].yOrder;l++) AY[n].a[l]=AY_R[n][l];
    for(l=0;l<AY[n].yOrder-1;l++) AY[n].yhist[l]=0.0;
    for(l=0;l<AY[n].xOrder-1;l++) AY[n].xhist[l]=0.0;
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

COMPLEX16 be;
INT4 k,m;
LIGOTimeGPS localgpsepoch=gpsepoch; /* Local variable epoch used to calculate the calibration factors */
REAL8 a[MAXALPHAS],b[MAXALPHAS];


  chanin_asq.type  = ADCDataChannel;
  chanin_darm.type = ADCDataChannel;
  chanin_exc.type  = ADCDataChannel;

  chanin_asq.name  = CLA.asq_chan;
  chanin_darm.name = CLA.darm_chan;
  chanin_exc.name  = CLA.exc_chan; 

  /* Get channel time step size by calling LALFrGetREAL4TimeSeries */

  LALFrSeek(&status,&localgpsepoch,framestream);
  LALFrGetPos(&status,&pos1,framestream);
  LALFrGetREAL4TimeSeries(&status,&asq,&chanin_asq,framestream);
  LALFrSetPos(&status,&pos1,framestream);
  LALFrGetREAL4TimeSeries(&status,&darm,&chanin_darm,framestream);
  LALFrSetPos(&status,&pos1,framestream);
  LALFrGetREAL4TimeSeries(&status,&exc,&chanin_exc,framestream);

  /* Allocate space for data vectors */
  LALCreateVector(&status,&asq.data,(INT4)(To/asq.deltaT +0.5));
  LALCreateVector(&status,&darm.data,(INT4)(To/darm.deltaT +0.5));
  LALCreateVector(&status,&exc.data,(INT4)(To/exc.deltaT +0.5));
  
  /* Create Window vectors */
  LALCreateVector(&status,&asqwin,(INT4)(To/asq.deltaT +0.5));
  LALCreateVector(&status,&darmwin,(INT4)(To/darm.deltaT +0.5));
  LALCreateVector(&status,&excwin,(INT4)(To/exc.deltaT +0.5));
   
  winparams.type=Hann;
   
  /* windows for time domain channels */
  /* asq */
  winparams.length=(INT4)(To/asq.deltaT +0.5);
  LALWindow(&status,asqwin,&winparams);
  
  /* darm */
  winparams.length=(INT4)(To/darm.deltaT +0.5);
  LALWindow(&status,darmwin,&winparams);

  /* exc */
  winparams.length=(INT4)(To/exc.deltaT +0.5);
  LALWindow(&status,excwin,&winparams);


  for(m=0;m<duration/To;m++)
    {
      /* Fill data vectors with data */

      LALFrSeek(&status,&localgpsepoch,framestream);
      LALFrGetPos(&status,&pos1,framestream);

      LALFrSetPos(&status,&pos1,framestream);
      LALFrGetREAL4TimeSeries(&status,&asq,&chanin_asq,framestream);

      LALFrSetPos(&status,&pos1,framestream);
      LALFrGetREAL4TimeSeries(&status,&darm,&chanin_darm,framestream);

      LALFrSetPos(&status,&pos1,framestream);
      LALFrGetREAL4TimeSeries(&status,&exc,&chanin_exc,framestream);

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
      params.outputMatrix = CLA.k;
      params.actuationFactor.re = Af0.re/ 3995.064 ;
      params.actuationFactor.im = Af0.im/ 3995.064 ;
      params.responseFactor = Rf0;
      params.sensingFactor = Cf0;
	  
      LALComputeCalibrationFactors(&status,&factors,&params);
      
      a[m]= factors.alpha.re;
      cdiv(&be,&factors.alphabeta,&factors.alpha);
      b[m] = be.re;

      //    fprintf(stdout,"%d %f  %f\n",m*To, a[m],b[m]);

      localgpsepoch.gpsSeconds = localgpsepoch.gpsSeconds+To;

    }

  /* Generate array of alphas with wings (to be used later for interpolations) */
  for(m=0;m<ALPHAWINGS;m++)
    {
      alpha[m]=a[0];
    }
  for(m=0;m<duration/To;m++)
    {
      alpha[m+ALPHAWINGS]=a[m];
    }
  for(m=0;m<ALPHAWINGS+2;m++)
    {
      alpha[m+ALPHAWINGS+duration/To]=a[duration/To-1];
    }

  /* Generate array of betas with wings (to be used later for interpolations) */
  for(m=0;m<ALPHAWINGS;m++)
    {
      beta[m]=b[0];
    }
  for(m=0;m<duration/To;m++)
    {
      beta[m+ALPHAWINGS]=b[m];
    }
  for(m=0;m<ALPHAWINGS+2;m++)
    {
      beta[m+ALPHAWINGS+duration/To]=b[duration/To-1];
    }
 
  /* Clean up */
  LALDestroyVector(&status,&darm.data);
  LALDestroyVector(&status,&exc.data);
  LALDestroyVector(&status,&asq.data);

  LALDestroyVector(&status,&asqwin);
  LALDestroyVector(&status,&darmwin);
  LALDestroyVector(&status,&excwin);

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
  R0.deltaF=1.0/8.0;                   /*ACHTUNG: HARDWIRED !!*/
  C0.f0=0.0;
  C0.deltaF=1.0/8.0;                   /*ACHTUNG: HARDWIRED !!*/
  A0.f0=0.0;
  A0.deltaF=1.0/8.0;                   /*ACHTUNG: HARDWIRED !!*/

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
     sscanf(line,"%d %d %d",&SL[i].gpsstart,&SL[i].gpsend,&SL[i].seglength);
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

 /* create Frame cache and open frame stream */
 LALFrCacheImport(&status,&framecache,CLA.FrCacheFile);
 LALFrCacheOpen(&status,&framestream,framecache);

 LALDestroyFrCache(&status,&framecache);

 LALZDestroyVector(&status,&R0.data);
 LALZDestroyVector(&status,&C0.data);
 LALZDestroyVector(&status,&A0.data);

 return 0;
}

/*******************************************************************************/

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
      /* calibration line frequency */
      CLA->f=atof(optarg);
      break;
    case 'k':
      /* darm to etmx output matrix value */
      CLA->k=atof(optarg);
      break;
    case 'F':
      /* name of frame cache file */
      CLA->FrCacheFile=optarg;
      break;
    case 'r':
      /* name of response file */
      CLA->RFile=optarg;
      break;
    case 'c':
      /* name of sensing file */
      CLA->CFile=optarg;
      break;
    case 'a':
      /* name of actuation file */
      CLA->AFile=optarg;
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
   case 'h':
      /* print usage/help message */
      fprintf(stdout,"All arguments are required. They are:\n");
      fprintf(stdout,"\t-f\tFLOAT\t Calibration line frequency in Hz.\n");
      fprintf(stdout,"\t-k\tFLOAT\t Output matrix darm to etmx.\n");
      fprintf(stdout,"\t-F\tSTRING\t Name of frame cache file.\n");
      fprintf(stdout,"\t-r\tSTRING\t Name of response function file.\n");
      fprintf(stdout,"\t-c\tSTRING\t Name of sensing function file.\n");
      fprintf(stdout,"\t-a\tSTRING\t Name of actuation function file.\n");
      fprintf(stdout,"\t-S\tSTRING\t Name of segment list file.\n");
      fprintf(stdout,"\t-A\tSTRING\t AS_Q channel name (eg, L1:LSC-AS_Q).\n");
      fprintf(stdout,"\t-E\tSTRING\t Excitation channel name (eg, L1:LSC-ETMX_EXC_DAQ)\n");
      fprintf(stdout,"\t-D\tSTRING\t Darm channel name (eg, L1:LSC-DARM_CTRL)\n");
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
  if(CLA->k == 0)
    {
      fprintf(stderr,"No value of the output matrix to x arm specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
  if(CLA->CFile == "")
    {
      fprintf(stderr,"No sensing function file specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
  if(CLA->RFile == "")
    {
      fprintf(stderr,"No response function file specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
  if(CLA->AFile == "")
    {
      fprintf(stderr,"No actuation function file specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
  if(CLA->FrCacheFile=="")
    {
      fprintf(stderr,"No frame cache file specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
  if(CLA->SegmentsFile=="")
    {
      fprintf(stderr,"No segments file specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
   if(CLA->exc_chan == "")
    {
      fprintf(stderr,"No excitation channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
   if(CLA->darm_chan == "")
    {
      fprintf(stderr,"No darm channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      
   if(CLA->asq_chan == "")
    {
      fprintf(stderr,"No asq channel specified.\n");
      fprintf(stderr,"Try ./ComputeStrain -h \n");
      return 1;
    }      

  return errflg;
}

/*******************************************************************************/

int FreeMem(void)
{

  LALFrClose(&status,&framestream);

  LALDDestroyVector(&status,&AS_Q.data);

  LALDDestroyVector(&status,&hC.data);
  LALDDestroyVector(&status,&hCx.data);
  LALDDestroyVector(&status,&hCy.data);
  LALDDestroyVector(&status,&hR.data);

  LALDDestroyVector(&status,&h.data);


  LALDDestroyVector(&status,&uphR.data);

  LALDDestroyVector(&status,&hipasssludge.data);
  LALDDestroyVector(&status,&lopasssludge.data);

  LALCheckMemoryLeaks();
 
  return 0;
}

/*******************************************************************************/

#endif
