/* <lalVerbatim file="SimulatePopcornTestCV">
Author: Tania Regimbau $Id$
</lalVerbatim> */

/*<lalLaTeX>
\subsection{Program \texttt{SimulatePopcornTest.c}}
\label{inject:ss:SimulatePopcornTest.c}

A program to test \texttt{LALSimPopcornTimeSeries()}

</lalLaTeX> */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConfig.h>
#include <lal/LALConstants.h>
#include <lal/LALStatusMacros.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/PrintFTSeries.h>
#include <lal/ReadFTSeries.h>
#include <lal/StreamInput.h>
#include <lal/Units.h>
#include <lal/PrintVector.h>
#include <lal/DetectorSite.h>
#include <lal/Random.h>
#include <lal/FileIO.h>
#include "SimulatePopcorn.h"

NRCSID (SIMULATEPOPCORNTESTC, "$Id$");



/**   define parameters here   **/

/*detectors */

#define LALPOPCORN_START 714265182
#define LALPOPCORN_RESP0 "H2.dat"
#define LALPOPCORN_RESP1 "L1.dat"
#define LALPOPCORN_NDATASET 2
#define LALPOPCORN_LENGTH 4
#define LALPOPCORN_SRATE 1024
#define LALPOPCORN_FREF 100.
#define LALPOPCORN_SEED 3
#define LALPOPCORN_WFDURATION 0.001
#define LALPOPCORN_LAMBDA 0.01
#define LALPOPCORN_SITE0 0
#define LALPOPCORN_SITE1 1
#define POPCORN_FILENAMEOUT0 "H2_"
#define POPCORN_FILENAMEOUT1 "L1_"
#define Pi LAL_PI

/* Commenting out the following definition turns off testing normalizations */
#define TEST
/* Commenting out the following definition turns off most printing */
#define PRINT
/*#define ASCII */
#define ILWD

/* <lalErrTable file="SimulatePopcornTestCE"> */
#define SIMULATEPOPCORNTESTC_ENOM 0
#define SIMULATEPOPCORNTESTC_EARG 1
#define SIMULATEPOPCORNTESTC_ECHK 2
#define SIMULATEPOPCORNTESTC_EFLS 3
#define SIMULATEPOPCORNTESTC_EUSE 4
#define SIMULATEPOPCORNTESTC_MSGENOM "Nominal exit"
#define SIMULATEPOPCORNTESTC_MSGEARG "Error parsing command-line arguments"
#define SIMULATEPOPCORNTESTC_MSGECHK "Error checking failed to catch bad data"
#define SIMULATEPOPCORNTESTC_MSGEFLS "Incorrect answer for valid data"
#define SIMULATEPOPCORNTESTC_MSGEUSE "Bad user-entered data"
/* </lalErrTable> */

/*sine gaussian waveform; */
static void wformfunc (REAL4 *result,REAL4 x)
 {
   REAL4 fo, to, duration;
   fo = 500.;
   duration = 0.001;
   to=duration/6.;
   *result =  exp(-(x*x)/(2.*to*to))*sin(2.*Pi*fo*x);
   return;
 }

int lalDebugLevel = LALNDEBUG;
int main (int argc, char* argv[])
{
  static LALStatus status;
  
  /* input and output structure declaration */
  SimPopcornInputStruc               PopcornInput;
  SimPopcornOutputStruc              PopcornOutput;
  SimPopcornParamsStruc              PopcornParams;
  
  /* output parameters declaration */
  REAL4TimeSeries                    Popcorn0;
  REAL4TimeSeries                    Popcorn1;

  REAL8                              totnorm2zero;
  REAL8                              totnorm2one;
  
  /* vector to store response functions of a pair of detectors  */
  COMPLEX8Vector *response[2]={NULL,NULL};
  COMPLEX8FrequencySeries response0;
  COMPLEX8FrequencySeries response1;
  LIGOTimeGPS t;
  
  /*others*/ 
  UINT4 i, j, n, nfreq, start;
  REAL4 deltat=(1./LALPOPCORN_SRATE);
  char fname[50];
  FILE *pfresp0,*pfresp1,*pfzero,*pfone;
 
  /* input parameters */
   PopcornInput.inputwform = &wformfunc;
   PopcornInput.inputduration = LALPOPCORN_WFDURATION;
   PopcornInput.inputlambda =  LALPOPCORN_LAMBDA; 
   PopcornInput.inputNdataset =  LALPOPCORN_NDATASET;
   PopcornInput.inputsite0 =  LALPOPCORN_SITE0;
   PopcornInput.inputsite1 =  LALPOPCORN_SITE1;
   
   PopcornParams.paramsstarttime =  LALPOPCORN_START;
   PopcornParams.paramslength =  LALPOPCORN_LENGTH;
   PopcornParams.paramssrate =  LALPOPCORN_SRATE;
   PopcornParams.paramsfref =  LALPOPCORN_FREF;
   PopcornParams.paramsseed =  LALPOPCORN_SEED;

  n = LALPOPCORN_LENGTH*LALPOPCORN_SRATE;
  nfreq = (n/2) +1;
 
  start = LALPOPCORN_START;
  for (i=0;i<2;i++)
    {
      LALCCreateVector(&status, &response[i],nfreq); 
       if (status.statusCode != 0) 
	{
	  fprintf(stderr,"[%i]: %s [CreateInVectors:%s]\n",status.statusCode,
		  status.statusDescription, SIMULATEPOPCORNTESTC_MSGEARG);
	  return SIMULATEPOPCORNTESTC_EARG;
	}
    }
 
  /*t.gpsSeconds = LALPOPCORN_START; */
  t.gpsSeconds = start;
  t.gpsNanoSeconds = 0;

  response0.deltaF = (1./LALPOPCORN_LENGTH); 
  response0.epoch = t;
  response0.f0 = 0.;
  response0.data = response[0];
  
  LALCReadFrequencySeries(&status, &response0, LALPOPCORN_RESP0);
  if (status.statusCode != 0) 
    {
      fprintf(stderr,"[%i]: %s [ReadFTSeriesTest1:%s]\n",status.statusCode,
	      status.statusDescription, SIMULATEPOPCORNTESTC_MSGEUSE);
      return SIMULATEPOPCORNTESTC_EUSE;
    }

  response1.deltaF = 1./LALPOPCORN_LENGTH;
  response1.epoch = t;
  response1.f0 = 0.;
  response1.data = response[1];

  LALCReadFrequencySeries(&status, &response1, LALPOPCORN_RESP1);
  if (status.statusCode != 0) 
    {
      fprintf(stderr,"[%i]: %s [ReadFTSeriesTest1:%s]\n",status.statusCode,
	      status.statusDescription, SIMULATEPOPCORNTESTC_MSGEUSE);
      return SIMULATEPOPCORNTESTC_EUSE;
    }
  
  PopcornInput.wfilter0 = &response0;
  PopcornInput.wfilter1 = &response1;
     
  

/* TEST INVALID DATA HERE */ 
#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
    {
     /*test behavior for null pointer to real time series for output */
     LALSimPopcornTimeSeries(&status,NULL,&PopcornInput,&PopcornParams);
     printf("  PASS: null pointer to output series results in error: \n"
	     "\"%s\"\n", SIMULATEPOPCORNH_MSGENULLP);
    }
  
#endif /* LAL_NDEBUG */


/* TEST VALID DATA HERE */
  Popcorn0.data = NULL;
  LALSCreateVector(&status, &(Popcorn0.data), n);
  Popcorn1.data = NULL;
  LALSCreateVector(&status, &(Popcorn1.data), n);
   
  PopcornOutput.SimPopcorn0 = &Popcorn0;
  PopcornOutput.SimPopcorn1 = &Popcorn1; 

  LALSimPopcornTimeSeries(&status,&PopcornOutput,&PopcornInput,&PopcornParams);
  
#ifdef TEST
  /* Mean square */
  totnorm2zero=0.0;
  for (j=0;j<n;j++)
    totnorm2zero+=((Popcorn0.data->data[j])*(Popcorn0.data->data[j]));
  totnorm2zero/=n;
  printf("Mean square of whitened output no. 1 is: %e\n",totnorm2zero);

  totnorm2one=0.0;
  for (j=0;j<n;j++)
    totnorm2one+=((Popcorn1.data->data[j])*(Popcorn1.data->data[j]));
  totnorm2one/=n;
  printf("Mean square of whitened output no. 2 is: %e\n",totnorm2one);
#endif    

#ifndef PRINT  
#ifndef ILWD

 LALSnprintf(fname,50,POPCORN_FILENAMEOUT0"%d.ilwd",start);
 pfzero=LALFopen(fname,"w");
 LALSnprintf(fname,50,POPCORN_FILENAMEOUT1"%d.ilwd",start);
 pfone=LALFopen(fname,"w");
 fprintf(pfzero,"<?ilwd?>\n");fprintf(pfone,"<?ilwd?>\n");
 
 fprintf(pfzero,"<ilwd comment='"POPCORN_FILENAMEOUT0"%d'",start); 
 fprintf(pfzero," name='"POPCORN_FILENAMEOUT0"%d' size='1'>\n",start);
 fprintf(pfzero,
  " <real_8 dims='1048576' name='"POPCORN_FILENAMEOUT0"%d'>",start);
 
 fprintf(pfone,"<ilwd comment='"POPCORN_FILENAMEOUT1"%d'",start); 
 fprintf(pfone," name='"POPCORN_FILENAMEOUT1"%d' size='1'>\n",start);
 fprintf(pfone,
  " <real_8 dims='1048576' name='"POPCORN_FILENAMEOUT1"%d'>",start);
 
 for(j=0;j<n;j++)
  {
   fprintf(pfzero,"% e",Popcorn0.data->data[i]);
   fprintf(pfone,"% e",Popcorn1.data->data[i]);
  }
 fprintf(pfzero,"</real_8>");fprintf(pfone,"</real_8>");
 fprintf(pfzero," </ilwd>");fprintf(pfone," </ilwd>");
#endif

#ifndef ASCII
  LALSnprintf(fname,50,POPCORN_FILENAMEOUT0"%d.dat",start);
  pfzero=LALFopen(fname,"w");
  LALSnprintf(fname,50,POPCORN_FILENAMEOUT1"%d.dat",start);
  pfone=LALFopen(fname,"w");

   for(j=0;j<n;j++)
   {
    fprintf(pfzero,"%f\t%e\n",(i*deltat),Popcorn0.data->data[i]);
    fprintf(pfone,"%f\t%e\n",(i*deltat),Popcorn1.data->data[i]);
    }
#endif 

 LALFclose(pfzero);LALFclose(pfone);
#endif

 /* clean up valid data */
 LALSDestroyVector(&status, &(Popcorn0.data));
 LALSDestroyVector(&status, &(Popcorn1.data));

   
 } 
  





