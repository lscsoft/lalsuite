/*
*  Copyright (C) 2007 Tania Regimbau
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

/*
 * popcorn.c 
 *
 * Tania Regimbau <regimbau@obs-nice.fr>  
 *
 *
 * $Id$
 */


/*
 gcc -Wall -o popcorn popcorn.c 
 -I/opt/lscsoft/non-lsc/include  -I${LAL_PREFIX}/ include 
 -L/usr/local/lib -L${LAL_PREFIX}/lib  -L/opt/lscsoft/libframe/lib  -L/opt/lscsoft/non-lsc/lib 
 -lgsl -lgslcblas -lm -llal -lz -llalframe -lFrame
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multimin.h>
#include <FrameL.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConfig.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/LALStatusMacros.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/ReadFTSeries.h>
#include <lal/StreamInput.h>
#include <lal/PrintVector.h>
#include <lal/VectorOps.h>
#include <lal/FileIO.h> 
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/ResampleTimeSeries.h>

NRCSID (POPCORNC, "$Id$");
RCSID ("$Id$");

/* cvs info */
#define CVS_ID "$Id$"
#define CVS_REVISION "$Revision$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "popcorn"

/* variables for getopt options parsing */
extern char *optarg;
extern int optind;

/* flag for getopt_long */
static int verbose_flag = 0;
static int gaussian_flag = 0;
static int ascii_flag = 0;
static int frame_flag = 0;
static int resample_flag = 0;
static int test_flag = 0;

UINT4 Npt =1000000;
UINT8 startTime = 700000000;
UINT8 stopTime = 700000150;
CHAR frameCache1 [200]= "H1.cache";
CHAR frameCache2[200] = "H2.cache";
CHAR channel1[LALNameLength]= "H1:STRAIN";
CHAR channel2[LALNameLength]= "H2:STRAIN";
UINT4 sampleRate = 16384;
UINT4 resampleRate = 1024;

REAL8 mu = 0.5;
REAL8 sigma = 1.;
REAL8 sigma1 = 1.;
REAL8 sigma2 = 1.;
REAL8 mu0 = 0.1;
REAL8 sigma0 = 0.1;
REAL8 v1, v2;
REAL8Vector *h,*s1,*s2,*e12;

double lambda (gsl_vector *x,void *params);
void parseOptions(INT4 argc, CHAR *argv[]);
void displayUsage(INT4 exitcode);

INT4 main (INT4 argc, CHAR *argv[])
 { 
   /** declaration **/
  /* random number generation */ 
  const gsl_rng_type * Trgn;
  gsl_rng * rrgn;

  /* minization */
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
  gsl_multimin_fminimizer *s;
  gsl_multimin_function minex_func;
  double size;
  size_t iter, k;
  
  int i, j, n;
  /* signal */
  double var, varn, sigman, varmean, snr;
  double pmu[10];
  double var1, var2, sigmaref;
  double muest, sigmaest, varest, varmeanest, sigma1est, sigma2est;
  double value;
  
  /* time series */
  UINT4 resampleFactor, Npt0;
  UINT4 duration, seglength, nsegment=1;
  REAL8TimeSeries n1, n2;
  REAL4TimeSeries n1Temp, n2Temp;

  /* frame variables */
  FrCache *frCache1,*frCache2;
  FrStream *frStream1,*frStream2;
  FrChanIn frChanIn1, frChanIn2;
  LIGOTimeGPS gpsStartTime;
  ResampleTSParams resampleParams;
  
  FILE *pf1,*pf2,*pf3;
  
  int status;
  static LALStatus lalstatus;
  
  /* parse command line options */
  parseOptions(argc, argv);   
  
  /** setup **/
  gsl_rng_env_setup();
  Trgn = gsl_rng_default;
  rrgn = gsl_rng_alloc (Trgn); 
  
  if(resample_flag){
   /* set resample parameters */
   resampleFactor=(UINT4)(sampleRate/resampleRate);
   Npt0=Npt*resampleFactor + 2*sampleRate;
   resampleParams.deltaT = 1./(REAL8)resampleRate;
   resampleParams.filterType = defaultButterworth;
  }
  
   /* set duration and number of jobs */
   duration = stopTime-startTime;
   if(resample_flag){
    seglength=(int)Npt/resampleRate;
    nsegment= (duration-2)/seglength;
   }
   else{
    seglength=(int)Npt/sampleRate;
    nsegment= duration/seglength;
   }
   
  if(verbose_flag){
   fprintf(stdout, "will analyze %d segments of length %d s...\n",nsegment, seglength);}
   
   /* poisson coefficient */
  for(i=0;i<10;i++)
   pmu[i]=gsl_ran_poisson_pdf(i,mu);
	
  /* signal-to-noise ratio rho = varmean*sqrt(N)/(sigma1*sigma2)*/
  var=sigma*sigma;
  varmean=0.; 
  for(i=1;i<10;i++)
   varmean=varmean+pmu[i]*var*i;
  snr=varmean*sqrt(Npt);

  /* open data files */
  if(ascii_flag){
   pf1=fopen("data/n1.dat","r");
   pf2=fopen("data/n2.dat","r");
  }
  
  gpsStartTime.gpsSeconds = startTime;
  gpsStartTime.gpsNanoSeconds = 0.;
  
   if(frame_flag){
   /* set channels */
   frChanIn1.name = channel1;
   frChanIn2.name = channel2;
   frChanIn1.type = SimDataChannel;
   frChanIn2.type = SimDataChannel;
 
   /* open frame cache */
   if(verbose_flag){
	fprintf(stdout,"Opening first frame cache...\n");}
   frCache1=NULL;frStream1=NULL;
   LALFrCacheImport(&lalstatus, &frCache1, frameCache1);
   LALFrCacheOpen( &lalstatus, &frStream1, frCache1);

   if(verbose_flag){
	fprintf(stdout, "Opening second frame cache...\n");}
   frCache2=NULL;frStream2=NULL;
   LALFrCacheImport( &lalstatus, &frCache2, frameCache2);
   LALFrCacheOpen( &lalstatus, &frStream2, frCache2);
  }
  
  /* open output file */
  pf3=fopen("popcorn.txt","a");
  
  fprintf(pf3,"number of segments=%d length of segments=%d s data point per segment=%d\n",nsegment,seglength,Npt);
  if(resample_flag){
   fprintf(pf3,"sampling rate=%d\n",resampleRate);}
  else{
   fprintf(pf3,"sampling rate=%d\n",sampleRate);}
  fprintf(pf3,"mu=%f sigma=%f\n",mu,sigma);
  
  /*** loop over segments ***/
  
  for(j=0;j<nsegment;j++){

   /* allocate memory for time series */	 
   h=NULL;n1.data=NULL; n2.data=NULL; s1=NULL; s2=NULL;e12=NULL;
  
   LALDCreateVector( &lalstatus, &h, Npt);
   memset(h->data, 0,h->length * sizeof(*h->data));
   LALDCreateVector( &lalstatus, &s1, Npt);
   memset(s1->data, 0,s1->length * sizeof(*s1->data));
   LALDCreateVector( &lalstatus, &s2, Npt);
   memset(s2->data, 0,s2->length * sizeof(*s2->data));
   LALDCreateVector( &lalstatus, &e12, Npt);
   memset(e12->data, 0,e12->length * sizeof(*e12->data));

   if(resample_flag){
    LALDCreateVector( &lalstatus, &n1.data, Npt0);
    memset(n1.data->data, 0,n1.data->length * sizeof(*n1.data->data));
    LALDCreateVector( &lalstatus, &n2.data, Npt0);
    memset(n2.data->data, 0,n2.data->length * sizeof(*n2.data->data));
   }
   else{
    LALDCreateVector( &lalstatus, &n1.data, Npt);
    memset(n1.data->data, 0,n1.data->length * sizeof(*n1.data->data));
    LALDCreateVector( &lalstatus, &n2.data, Npt);
    memset(n2.data->data, 0,n2.data->length * sizeof(*n2.data->data));
   }
  
  
   /** noise **/
   gpsStartTime.gpsSeconds = startTime+j*seglength;
   if(gaussian_flag){
    if(verbose_flag){
     fprintf(stdout, "generate gaussian noise with variance sigma1=%f and sigma2=%f...\n",sigma1,sigma2);}
   /* generate gaussian noise */ 
   for(i=0;i<Npt;i++){
    n1.data->data[i] = gsl_ran_gaussian (rrgn,sigma1);
    n2.data->data[i] = gsl_ran_gaussian (rrgn,sigma2);
   }
   sigmaref=sqrt(sigma1*sigma2);
  }
  else{
   /* read noise from file */
  
   /* from ascii */
   if(ascii_flag){
    if(verbose_flag){
	 fprintf(stdout, "read data from ascii files...\n");}
	for (i = 0; i < Npt; i++) {
	 fscanf(pf1,"%lf",&value);
	 n1.data->data[i]=value;
	}
    for (i = 0; i < Npt; i++) {
     fscanf(pf2,"%lf",&value);
	 n2.data->data[i]=value;
	}
   }
   else{
   
    /* read from frames */
	
    /* read channels */
    if(verbose_flag){
     fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn1.name);}
    LALFrSeek( &lalstatus, &gpsStartTime, frStream1);
    LALFrGetREAL8TimeSeries(&lalstatus, &n1, &frChanIn1, frStream1);
	
    if(verbose_flag){		     
     fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn2.name);}
    LALFrSeek(&lalstatus, &gpsStartTime, frStream2);
    LALFrGetREAL8TimeSeries(&lalstatus, &n2,&frChanIn2, frStream2);
	
	 if(test_flag){
   for(i=0;i<10;i++){
    fprintf(stdout,"%e %e\n",n1.data->data[i],n2.data->data[i]);}}
    }
	
	if(resample_flag){
	 
	 n1Temp.data=NULL; n2Temp.data=NULL;
	 LALCreateVector( &lalstatus, &n1Temp.data, Npt0);
	 memset(n1Temp.data->data, 0,n1Temp.data->length * sizeof(*n1Temp.data->data));
	 LALCreateVector( &lalstatus, &n2Temp.data, Npt0);
	 memset(n2Temp.data->data, 0,n2Temp.data->length * sizeof(*n2Temp.data->data));
	 
	 /* cast to REAL4  */ 
	 if(verbose_flag){
	  fprintf(stdout, "cast to single precision...\n");} 
	 for (i=0;i<Npt0;i++){           
	  n1Temp.data->data[i]= (REAL4)n1.data->data[i];         
	  n2Temp.data->data[i]= (REAL4)n2.data->data[i];
	  }
     
     /* resample */
     if(verbose_flag){
      fprintf(stdout, "Resampling first data stream to %d Hz...\n",resampleRate);}
     LALResampleREAL4TimeSeries(&lalstatus,&n1Temp,&resampleParams);	
     if(verbose_flag){
      fprintf(stdout, "Resampling second data stream to %d Hz...\n",resampleRate);}
     LALResampleREAL4TimeSeries(&lalstatus,&n2Temp,&resampleParams);
	 
	 /* cast to REAL8 and throw away first and last seconds */
	 if(verbose_flag){
      fprintf(stdout, "cast to double precision and throw away first and last seconds...\n");}
     for (i=0;i<Npt;i++){           
	  n1.data->data[i]= (REAL8)n1Temp.data->data[i+resampleRate];
	  n2.data->data[i]= (REAL8)n2Temp.data->data[i+resampleRate];
	 }
     LALDestroyVector(&lalstatus,&n1Temp.data);
     LALDestroyVector(&lalstatus,&n2Temp.data);
	 
	}
	 
   /* normalize so that maximal sigma1*sigma2 = 1 **/ 
   if(verbose_flag){
    fprintf(stdout, "calculate noise variances...\n");}
   var1=0.;var2=0.;	   
   for (i = 0; i < Npt; i++) {
    var1 = var1 + n1.data->data[i]*n1.data->data[i];
    var2 = var2 + n2.data->data[i]*n2.data->data[i];
   }
   sigma1=sqrt(var1/Npt);
   sigma2=sqrt(var2/Npt);
   sigmaref=sqrt(sigma1*sigma2);
   if(verbose_flag){
    fprintf(stdout, "normalize noise...\n");}	 
    
   sigma1=sigma1/sigmaref; sigma2=sigma2/sigmaref;	
   for (i = 0; i < Npt; i++) {
    n1.data->data[i] = n1.data->data[i]/sigmaref;
    n2.data->data[i] = n2.data->data[i]/sigmaref;
   }
  }
  	
  /** generate gw signal **/

  if(verbose_flag){
   fprintf(stdout, "generate gw signal with mu=%f and sigma=%f...\n",mu,sigma);}
   	
  v1=0.;v2=0.;	
  for (i = 0; i < Npt; i++) {
   n = gsl_ran_poisson (rrgn,mu);           
   if(n>0){
	varn = var*(double)n;                                          
	sigman = sqrt(varn);        
	h->data[i] = gsl_ran_gaussian (rrgn,sigman);
   }
   s1->data[i] = n1.data->data[i]+h->data[i];                                 
   s2->data[i] = n2.data->data[i]+h->data[i];
   v1=v1+s1->data[i]*s1->data[i];
   v2=v2+s2->data[i]*s2->data[i];                               
  }
  v1=v1/Npt; v2=v2/Npt;
	
  /** maximize likelihood function for parameter estimation **/ 
  /* Nelder-Mead Simplex algorithm */
  
  if(verbose_flag){
   fprintf(stdout, "estimate parameters...\n");}
  
  gsl_vector *ss=gsl_vector_alloc (4), *x=gsl_vector_alloc (4);
  
  gsl_vector_set_all (ss, 0.1);
  gsl_vector_set(x,0,mu0);
  gsl_vector_set(x,1,sigma0);
  gsl_vector_set(x,2,1.);
  gsl_vector_set(x,3,1.);

  minex_func.f = &lambda;
  minex_func.n = 4;
  minex_func.params = NULL;
  
  s=NULL; 
  s = gsl_multimin_fminimizer_alloc (T, 4);
  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
  
  iter=0;
  do{
   iter++;
   status = gsl_multimin_fminimizer_iterate(s);
   if (status){ 
	break;}
   size = gsl_multimin_fminimizer_size (s);
   status = gsl_multimin_test_size (size, 1.e-4);
   if(test_flag){
	if (status == GSL_SUCCESS){
	 fprintf (stdout,"converged to minimum at\n");}
	fprintf (stdout,"%5d ", (int)iter);
	for (k = 0; k < 4; k++){
	 fprintf (stdout,"%10.3e ", gsl_vector_get (s->x, k));}
	fprintf (stdout,"f = %e\n", s->fval);
   }
  }
 while (status == GSL_CONTINUE && iter < 300);

 muest=gsl_vector_get (s->x, 0);
 sigmaest=gsl_vector_get (s->x, 1);
 sigma1est=gsl_vector_get (s->x, 2);
 sigma2est=gsl_vector_get (s->x, 3);
 
 gsl_multimin_fminimizer_free (s);
 gsl_vector_free(x);
 gsl_vector_free(ss);
 
 LALDDestroyVector(&lalstatus,&h);
 LALDDestroyVector(&lalstatus,&n1.data);
 LALDDestroyVector(&lalstatus,&n2.data);
 LALDDestroyVector(&lalstatus,&s1);
 LALDDestroyVector(&lalstatus,&s2);
 LALDDestroyVector(&lalstatus,&e12);
  
 /* poisson coefficient */
 for(i=0;i<10;i++){
  pmu[i]=gsl_ran_poisson_pdf(i,muest);}
 varest=sigmaest*sigmaest;
 varmeanest=0.; 
 for(i=1;i<10;i++){
  varmeanest=varmeanest+pmu[i]*varest*i;}
  
 
 if(verbose_flag){
  fprintf(stdout,"T=%d:",gpsStartTime.gpsSeconds);
  fprintf(stdout,"mu=%f sigma=%f mean(var)=%f sigma1=%f sigma2=%f\n",mu,sigma,varmean,sigma1,sigma2);
  fprintf(stdout,"muest=%f sigmaest=%f mean(varmest)=%f sigma1est=%f sigma2est=%f\n",muest,sigmaest,varmeanest,sigma1est,sigma2est);
  fprintf(stdout,"snr=%f\n",snr);
  fprintf(stdout,"\n");}
  
 fprintf(pf3,"T=%d:",gpsStartTime.gpsSeconds);
 fprintf(pf3,"mean(var)=%f sigma1=%f sigma2=%f\n",varmean,sigma1,sigma2);
 fprintf(pf3,"muest=%f sigmaest=%f varmeanest=%f sigma1est=%f sigma2est=%f\n",muest,sigmaest,varmeanest,sigma1est,sigma2est);
 fprintf(pf3,"snr=%f\n",snr);

 }  
  
 /** cleanup **/
 if(verbose_flag){
  fprintf(stdout,"cleanup and exit\n");}
  
 
 gsl_rng_free (rrgn);
 
 fclose(pf1); fclose(pf2); fclose(pf3);
 
 /* close frame cache */
 LALFrClose(&lalstatus, &frStream1);
 LALFrClose(&lalstatus, &frStream2);
 
 return 0;
}
 
 
 double lambda (gsl_vector *x,void *params)
{

  int i,j,np;
  double y;
  double dsum,psum;
  double pmu,psigma,psigma1,psigma2;
  double pv,pv1,pv2,psig12,pv12,pvv1,pvv2,v1opv1,v2opv2;
  double pc[10],a[10],b[10];  

  pmu=gsl_vector_get(x,0);
  psigma=gsl_vector_get(x,1);
  psigma1=gsl_vector_get(x,2);
  psigma2=gsl_vector_get(x,3);

  if((pmu>0.)&&(pmu<1.)&&(psigma>0.)&&(psigma1>0.)&&(psigma2>0.)){

   pv=psigma*psigma;
   pv1=psigma1*psigma1;
   pv2=psigma2*psigma2;
   psig12=psigma1*psigma2;
   pv12=psig12*psig12;
   pvv1=pv*pv1;
   pvv2=pv*pv2;
   v1opv1=v1/pv1;
   v2opv2=v2/pv2;

   /* poisson coefficient */
   for(i=0;i<10;i++)
     pc[i]=gsl_ran_poisson_pdf(i,pmu);
    
   i=1;
   while(pc[i]>10./(double)Npt){
    a[i]=1./sqrt(pv12+i*(pvv1+pvv2))*pc[i];
    b[i]=2.*(1./pv1+1./pv2+1./((double)i*pv));
    i++;
    }
   np=i;dsum=0.;
   for(i=0;i<Npt;i++){
     e12->data[i]=(s1->data[i]/pv1+s2->data[i]/pv2);
     e12->data[i]= e12->data[i]*e12->data[i];
     psum=0.;   
     for(j=1;j<np;j++){
       psum = psum + a[j]*exp(e12->data[i]/b[j]);
     }
     psum = pc[0]+psig12*psum;
     dsum=dsum+log(psum);
   }
   y=-(0.5*(log(v1opv1*v2opv2)-(v1opv1+v2opv2))+1.+(1./Npt)*dsum);
  }
  else y=10000.;
  return y; 

}

	
 /* parse command line options */
void parseOptions(INT4 argc, CHAR *argv[])
 {
  int c = -1;

  while(1)
   {
    static struct option long_options[] =
     {
	  /* options that set a flag */
      {"verbose", no_argument, &verbose_flag, 1},
	  {"gaussian", no_argument, &gaussian_flag, 1},
      {"ascii", no_argument, &ascii_flag, 1},
	  {"frame", no_argument, &frame_flag, 1},
	  {"resample", no_argument, &resample_flag, 1},
	  {"test", no_argument, &test_flag, 1},
      /* options that don't set a flag */
      {"help", no_argument, 0, 'h'},
      {"gps-start-time", required_argument, 0, 't'},
	  {"gps-stop-time", required_argument, 0, 'T'},
      {"number-points", required_argument, 0, 'N'},
	  {"sampleRate", required_argument, 0, 'r'},
      {"resampleRate", required_argument, 0, 'R'},
      {"mu", required_argument, 0, 'm'},
	  {"mu0", required_argument, 0, 'M'},
      {"sigma", required_argument, 0, 's'},
	  {"sigma0", required_argument, 0, 'S'},
	  {"sigma1", required_argument, 0, 'g'},
	  {"sigma2", required_argument, 0, 'G'},
      {"channel1", required_argument, 0, 'c'},
      {"channel2", required_argument, 0, 'C'},
      {"frame-cache1", required_argument, 0, 'd'},
      {"frame-cache2", required_argument, 0, 'D'},
      {"version", no_argument, 0, 'v'},
      {0, 0, 0, 0}
     };

    /* getopt_long stores the option here */
    int option_index = 0;

    c = getopt_long(argc, argv, 
                  "ht:T:N:r:R:m:M:s:S:g:G:c:C:d:D:v",
 		   long_options, &option_index);

    if (c == -1)
     {
      /* end of options, break loop */
      break;
     }

    switch(c)
     {
      case 0:
             /* If this option set a flag, do nothing else now. */
             if (long_options[option_index].flag != 0)
              break;
             printf ("option %s", long_options[option_index].name);
             if (optarg)
              printf (" with arg %s", optarg);
             printf ("\n");
             break;

      case 'h':
               /* HELP!!! */
               displayUsage(0);
               break;

      case 't':
			   /* start time */
	           startTime = atoi(optarg);
	           break;
			   
	  case 'T':
			   /* stop time */
	           stopTime = atoi(optarg);
	           break;

      case 'N':
			   /* number of points in time serie */
	           Npt = atoi(optarg);
	           break;
			   
	  case 'r':
			   /* sample rate */
	           sampleRate = atof(optarg);
	           break;
			   
	  case 'R':
			   resampleRate= atof(optarg);
	           break;
			   
	  case 'm':
	           mu = atof(optarg);
	           break;
			   
	  case 'M':
	           mu0 = atof(optarg);
	           break;
			   
	  case 's':
	           sigma = atof(optarg);
	           break;
			   
	  case 'S':
	           sigma0 = atof(optarg);
	           break;

			   
	  case 'g':
	           sigma1 = atof(optarg);
	           break;
			   
	  case 'G':
	           sigma2 = atof(optarg);
	           break;
     
      case 'c':
	           /* ifo for first stream */
	           strncpy(channel1, optarg, LALNameLength);
               break;

      case 'C':
	           /* ifo for first stream */
	           strncpy(channel2, optarg, LALNameLength);
               break;

      case 'd':
               /* data cache one */
               strncpy(frameCache1, optarg, 200);
               break;

      case 'D':
               /* data cache two */
               strncpy(frameCache2, optarg, 200);
               break;

      case 'v':
	     /* display version info and exit */
	     fprintf(stdout, "non gaussian stochastic background pipeline\n" CVS_ID "\n");
	     exit(0);
	     break;

     default:
     displayUsage(1);
     }
    }

   if (optind < argc)
    {
     displayUsage(1);
    }

  return;
}

/* display program usage */
void displayUsage(INT4 exitcode)
 {
  fprintf(stderr, "Usage: pipeline [options]\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, " -h                    print this message\n");
  fprintf(stderr, " -v                    display version\n");
  fprintf(stderr, " --verbose             verbose mode\n");
  fprintf(stderr, " --test                test mode\n");
  fprintf(stderr, " --frame               read data from frames\n");
  fprintf(stderr, " --ascii               read data from ascii files\n");
  fprintf(stderr, " --gaussian            simulate gaussian noise\n");
  fprintf(stderr, " --resample            resample data\n");
  fprintf(stderr, " -t                    GPS start time\n");
  fprintf(stderr, " -T                    GPS stop time\n");
  fprintf(stderr, " -N                    number of points per data segment\n");
  fprintf(stderr, " -r                    sampleRate\n");
  fprintf(stderr, " -R                    resampleRate\n");
  fprintf(stderr, " -m                    mu of injected gw signal\n");
  fprintf(stderr, " -M                    initial mu for parameter estimation\n");
  fprintf(stderr, " -s                    sigma of injected gw signal\n");
  fprintf(stderr, " -S                    initial sigma for parameter estimation\n");
  fprintf(stderr, " -g                    sigma1 of gaussian noise \n");
  fprintf(stderr, " -G                    sigma2 of gaussian noise \n");
  fprintf(stderr, " -c                    channel for first stream\n");
  fprintf(stderr, " -C                    channel for second stream\n");
  fprintf(stderr, " -d                    cache file for first stream\n");
  fprintf(stderr, " -D                    cache file for second stream\n");
  exit(exitcode);
}