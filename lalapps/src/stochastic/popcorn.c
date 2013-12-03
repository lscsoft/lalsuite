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
#include <gsl/gsl_statistics.h>
#include <lal/LALFrameL.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConfig.h>
#include <lal/Units.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStatusMacros.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/ReadFTSeries.h>
#include <lal/StreamInput.h>
#include <lal/PrintVector.h>
#include <lal/VectorOps.h>
#include <lal/FileIO.h> 
#include <lal/LALCache.h>
#include <lal/LALFrStream.h>
#include <lal/ResampleTimeSeries.h>


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
static int normalize_flag = 0;
static int test_flag = 0;
static int montecarlo_flag = 0;
static int condor_flag = 0;
static int post_flag = 0;
static int stat_flag = 0;
static int contour_flag = 0;

UINT4 job=1;
UINT4 Npt =1000000;

UINT8 startTime = 700000000;
UINT8 stopTime = 700000150;
UINT4 duration = 10000;
CHAR frameCache1 [200]= "H1.cache";
CHAR frameCache2[200] = "H2.cache";
CHAR channel1[LALNameLength]= "H1:STRAIN";
CHAR channel2[LALNameLength]= "H2:STRAIN";

/*
UINT8 startTime = 816071371 ;
UINT8 stopTime = 816073004;
UINT4 duration = 10000;
CHAR frameCache1[200]= "S5H1.cache";
CHAR frameCache2[200] = "S5L1.cache";
CHAR channel1[LALNameLength]= "H1:LSC-STRAIN";
CHAR channel2[LALNameLength]= "L1:LSC-STRAIN";
*/

UINT4 sampleRate = 16384;
UINT4 resampleRate = 1024;

UINT4 stat=2;
UINT4 mcstat=2;
REAL8 mu = 0.5;
REAL8 ksi = -1.;
REAL8 sigma = 0.5;
REAL8 sigma1_ref = 1.;
REAL8 sigma2_ref = 1.;
REAL8 mu0 = 0.1;
REAL8 ksi0=-1.;
REAL8 sigma0 = 0.1;

UINT4 nmax = 1;
REAL8 v1, v2, v12;
REAL8Vector *h,*s1,*s2,*e12,*hgw;

double lambda0 (gsl_vector *x,void *params);
double lambda1 (gsl_vector *x,void *params);
double lambda2 (gsl_vector *x,void *params);
void parseOptions(INT4 argc, CHAR *argv[]);
void displayUsage(INT4 exitcode);

INT4 main (INT4 argc, CHAR *argv[])
 { 
   /** declaration **/
  /* random number generation */ 
  const gsl_rng_type * Trgn;
  gsl_rng * rrgn;

  /* minization */
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s;
  gsl_multimin_function minex_func;
  double size;
  size_t iter, k;
  
  int i,j,n,l,m,Ne;
  /* signal */
  double var,varn,sigman,varmean,snr;
  double meangw,vargw,skewgw,kurtgw;
  double sigma1,sigma2,var1,var2,sigmaref;
  double ksi_ml,sigma_ml,CP[1000][1000],cp0,pmax;
  double muest, ksiest, sigmaest, varest, varmeanest, sigma1est, sigma2est;
  double muestMean, ksiestMean, sigmaestMean, varmeanestMean, sigma1estMean, sigma2estMean;
  double value;
  
  /* time series */
  UINT4 resampleFactor, Npt0;
  UINT4 seglength, nsegment=1;
  REAL8TimeSeries n1, n2;
  REAL4TimeSeries n1Temp, n2Temp;

  /* frame variables */
  LALCache *frCache1,*frCache2;
  LALFrStream *frStream1,*frStream2;
  FrChanIn frChanIn1, frChanIn2;
  LIGOTimeGPS gpsStartTime;
  ResampleTSParams resampleParams;
  
  FILE *pf1,*pf2,*pf3,*pf4,*pf5;
  /* output file name */
  CHAR fileName[LALNameLength];

  int status;
  static LALStatus lalstatus;
  lalstatus.statusPtr = NULL;
  
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
   if(condor_flag){
    startTime=startTime+(job-1)*duration;
    if(stopTime>startTime+duration)
     stopTime=startTime+duration;
   }
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
   
  /* probability of non null signal */
  if (mcstat==0)
   ksi=1; 
  else if (mcstat==2) 
   ksi=1.-gsl_ran_poisson_pdf(0,mu);
  
  /* signal-to-noise ratio rho = varmean*sqrt(N)/(sigma1*sigma2)*/
  var=sigma*sigma;
  snr=var*sqrt(Npt)/(sigma1_ref*sigma2_ref);
  if (mcstat==1){
   snr=ksi*snr;
   varmean=ksi*var;}
  else if (mcstat==2){
   snr=mu*snr;
   varmean=mu*var;}

  /* open data files */
  if(ascii_flag){
   if(verbose_flag){
	fprintf(stdout,"Opening ascii files...\n");}
   pf1=fopen("data/s1.dat","r");
   pf2=fopen("data/s2.dat","r");
  }
  
  gpsStartTime.gpsSeconds = startTime;
  gpsStartTime.gpsNanoSeconds = 0.;
  
   if(frame_flag){
   /* set channels */
   frChanIn1.name = channel1;
   frChanIn2.name = channel2;
   frChanIn1.type = ProcDataChannel;
   frChanIn2.type = ProcDataChannel;
 
   /* open frame cache */
   if(verbose_flag){
	fprintf(stdout,"Opening first frame cache...\n");}
   frCache1=NULL;frStream1=NULL;
   frCache1 = XLALCacheImport( frameCache1 );
   LALFrCacheOpen( &lalstatus, &frStream1, frCache1);

   if(verbose_flag){
	fprintf(stdout, "Opening second frame cache...\n");}
   frCache2=NULL;frStream2=NULL;
   frCache2 = XLALCacheImport( frameCache2 );
   LALFrCacheOpen( &lalstatus, &frStream2, frCache2);
  }
  
  /* name output file */
  snprintf(fileName, LALNameLength,"output/stat_%d.dat",job);
  pf3 = fopen(fileName,"w");
  
  
  /*** loop over segments ***/
  muestMean=0.;ksiestMean=0.;sigmaestMean=0.;varmeanestMean=0.;
  sigma1estMean=0.;sigma2estMean=0.;
  for(j=0;j<nsegment;j++){
   sigma1=sigma1_ref; sigma2=sigma2_ref;
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
    n2.data->data[i] = gsl_ran_gaussian (rrgn,sigma2);   }
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
   }//end if(ascii_flag)
   else{
    /* read from frames */
	
    /* read channels */
    if(verbose_flag){
     fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn1.name);}
    LALFrSeek( &lalstatus, &gpsStartTime, frStream1);
	//LALFrGetREAL4TimeSeries(&lalstatus, &n1Temp, &frChanIn1, frStream1);
    LALFrGetREAL8TimeSeries(&lalstatus, &n1, &frChanIn1, frStream1);
	
    if(verbose_flag){		     
     fprintf(stdout, "Reading in channel \"%s\"...\n", frChanIn2.name);}
    LALFrSeek(&lalstatus, &gpsStartTime, frStream2);
    LALFrGetREAL8TimeSeries(&lalstatus, &n2,&frChanIn2, frStream2);
    }//end else read from frames
	
    //for(i=0;i<10;i++){
 	 //fprintf(stdout,"%e %e\n",n1.data->data[i],n2.data->data[i]);}
	
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
	}//end if(resample_flag)
	 
	 
   }//end else read from files
  
  
  /* normalize so that maximal sigma1*sigma2 = 1 **/ 
   if((stat_flag)||(normalize_flag)){
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
     fprintf(stdout, "sigma1=%e sigma2=%e sigmaref=%e\n",sigma1,sigma2,sigmaref);}
    	
    if(normalize_flag){
     if(verbose_flag)
      fprintf(stdout, "normalize noise...\n");	 
     sigma1=sigma1/sigmaref; sigma2=sigma2/sigmaref;	
     for (i = 0; i < Npt; i++) {
      n1.data->data[i] = n1.data->data[i]/sigmaref;
      n2.data->data[i] = n2.data->data[i]/sigmaref;
     }
    }
   }
  
  if(verbose_flag){
   if(montecarlo_flag){
    if(mcstat==1)
	 fprintf(stdout, "generate gw signal with ksi=%f and sigma=%f\n",ksi,sigma);
	else if(mcstat==2)
     fprintf(stdout, "generate gw signal with mu=%f and sigma=%f\n",mu,sigma);
    else
	 fprintf(stdout, "generate gw signal with sigma=%f\n",sigma);}
	 	
	fprintf(stdout, "calculate variances v1 and v2...\n");}
	
  v1=0.;v2=0.;v12=0.;Ne=0;
  /* calculate variance */		
  for (i = 0; i < Npt; i++) {
   if(montecarlo_flag){
   /* generate gw signal */
	n=1;
    if(mcstat==1){
	 if(gsl_ran_flat (rrgn,0.,1.)>ksi)
	  n=0;}
    else if(mcstat==2)
     n = gsl_ran_poisson (rrgn,mu);           
    if(n>0){
	 Ne++;
     varn = var*(double)n;                                          
     sigman = sqrt(varn);
	 h->data[i] = gsl_ran_gaussian (rrgn,sigman);
	 //fprintf(stdout, "%e\n",h->data[i]);
    }
	s1->data[i] = n1.data->data[i]+h->data[i];                                 
    s2->data[i] = n2.data->data[i]+h->data[i];
   }
   else{
    s1->data[i] = n1.data->data[i];                                 
    s2->data[i] = n2.data->data[i];
   }	

   v1=v1+s1->data[i]*s1->data[i];
   v2=v2+s2->data[i]*s2->data[i];   
   if(stat==0){v12=v12+s1->data[i]*s2->data[i];}                            
  }
  
   
   /*
   // output data for test with Matlab
   pf5=fopen("data.dat","w");
   for(i=0;i<Npt;i++)
    fprintf(pf5,"%e %e\n",s1->data[i],s2->data[i]);
   */	
  
  /*statistics of the GW signal*/
  
   if(stat_flag){
    fprintf(stdout,"statistics of the GW signal\n");
    fprintf(stdout,"number of data points containing a GW signal: %d or a ratio of %e\n",Ne, (double)Ne/(double)Npt);
	fprintf(stdout,"calculate moments of the distribution of the amplitude...\n");
	hgw=NULL;
    LALDCreateVector( &lalstatus, &hgw, Ne);
    memset(hgw->data, 0,hgw->length * sizeof(*hgw->data));
	l=0;
	for(i=0;i<Npt;i++) {
	  if(fabs(h->data[i])>0.){
	   hgw->data[l]=h->data[i];
	   l++;
	 }}
	
	meangw= gsl_stats_mean(hgw->data, 1, Ne);
    vargw = gsl_stats_variance(hgw->data, 1, Ne);
	skewgw= gsl_stats_skew(hgw->data, 1, Ne);
    kurtgw = gsl_stats_kurtosis(hgw->data, 1, Ne);
	fprintf(stdout,"mean=%e variance=%e skewness=%e kurtosis=%e\n",meangw,vargw,skewgw,kurtgw);
	}
	
  v1=v1/Npt;v2=v2/Npt;v12=v12/Npt;
  
  /*     
  if(test_flag){
   for(i=0;i<Npt;i++){
	fprintf(stdout,"%e %e\n",n1.data->data[i],n2.data->data[i]);}}
  */
  
  /** maximize likelihood function for parameter estimation **/ 
  /* Nelder-Mead Simplex algorithm */
  
  if(verbose_flag){
   fprintf(stdout, "estimate parameters...\n");
   fprintf(stdout,"T=%d:",gpsStartTime.gpsSeconds);
   fprintf(stdout,"snr=%f\n",snr);
   if(mcstat==1)
    fprintf(stdout,"ksi=%e sigma=%e varmean=%e sigma1=%e sigma2=%e\n",
	                ksi,sigma,varmean,sigma1,sigma2);
   else if(mcstat==2)
    fprintf(stdout,"mu=%e ksi=%e sigma=%e varmean=%e sigma1=%e sigma2=%e\n",
	                mu,ksi,sigma,varmean,sigma1,sigma2);
    else 
    fprintf(stdout,"sigma=%e varmean=%e sigma1=%e sigma2=%e\n",
	                sigma,varmean,sigma1,sigma2);}	
	
  fprintf(pf3,"%d ",gpsStartTime.gpsSeconds);
  fprintf(pf3,"%f %f ",sigma1,sigma2); 
  
  if(stat==2){
   if(mcstat==2) 
    {ksi0=1.-gsl_ran_poisson_pdf(0,mu0);}
   pmax=1.;
   for(i=0;i<=nmax;i++)
    {pmax=pmax-gsl_ran_poisson_pdf(i,mu);}
   fprintf(stdout,"nmax=%d (%e percent left out)\n",nmax,pmax*100.);
   gsl_vector *ss=gsl_vector_alloc (4), *x=gsl_vector_alloc (4);
   gsl_vector_set_all (ss, 0.1);
   gsl_vector_set(x,0,ksi0);
   gsl_vector_set(x,1,sigma0);
   gsl_vector_set(x,2,1.);
   gsl_vector_set(x,3,1.);
   minex_func.f = &lambda2;
   minex_func.n = 4;
   minex_func.params = NULL;
   s=NULL; 
   s = gsl_multimin_fminimizer_alloc (T, 4);
   gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
   iter=0;
   do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status){break;}
    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, 1.e-6);
	
    if(test_flag){
 	 if (status == GSL_SUCCESS){
	  fprintf (stdout,"converged to minimum at\n");}
	 fprintf (stdout,"%5d ", (int)iter);
	 for (k = 0; k < 4; k++){
	  fprintf (stdout,"%10.3e ", gsl_vector_get (s->x, k));}
	 fprintf (stdout,"f = %e\n", s->fval);
    }
   }
   while (status == GSL_CONTINUE && iter < 1000);
   ksiest=gsl_vector_get (s->x, 0);
   muest=-log(1.-ksiest);
   sigmaest=gsl_vector_get (s->x, 1);
   sigma1est=gsl_vector_get (s->x, 2);
   sigma2est=gsl_vector_get (s->x, 3);
   varest=sigmaest*sigmaest;
   varmeanest=muest*varest;
   muestMean=muestMean+muest;
   ksiestMean=ksiestMean+ksiest;
   if(verbose_flag)
    fprintf(stdout,"muest=%e ksiest=%e sigmaest=%e varest=%e varmeanest=%e sigma1est=%e sigma2est=%e\n",
	                muest,ksiest,sigmaest,varest,varmeanest,sigma1est,sigma2est);
   fprintf(pf3,"%f %f %f %f %f %f\n",muest,ksiest,sigmaest,varmeanest,sigma1est,sigma2est);
   
   if(stat_flag){   
	gsl_vector_set(x,0,ksiest);
	gsl_vector_set(x,1,sigmaest);
	gsl_vector_set(x,2,sigma1est);
	gsl_vector_set(x,3,sigma2est);
	cp0=lambda2(x,NULL);	
	printf("likelihood function at estimated parameters: %e\n",-(double)Npt*cp0);
    gsl_vector_set(x,0,ksi);
	gsl_vector_set(x,1,sigma);
	gsl_vector_set(x,2,sigma1);
	gsl_vector_set(x,3,sigma2);
	cp0=lambda2(x,NULL);
	printf("likelihood function at injected parameters: %e\n",-(double)Npt*cp0);
   }
   
   //test for comparison with fig 11 of Drasc's paper
   if(contour_flag){
    pf4=fopen("contour.dat","w");
    for(m=1;m<100;m++){
	 for(l=1;l<100;l++){
	  ksi_ml=(double)m*0.01;
	  sigma_ml=sqrt((double)l*0.01);
	  gsl_vector_set(x,0,ksi_ml);
      gsl_vector_set(x,1,sigma_ml);
      gsl_vector_set(x,2,1.);
      gsl_vector_set(x,3,1.);
	  CP[m][l]=-(double)Npt*lambda2(x,NULL);
	  }}
	for(m=1;m<100;m++){
	 for(l=1;l<100;l++){ 
	  fprintf(pf4,"%e\t",CP[m][l]);}
	  fprintf(pf4,"\n");}
    fclose(pf4);
   }
   gsl_multimin_fminimizer_free (s);
   gsl_vector_free(x); 
   gsl_vector_free(ss); 
  }
  else if(stat==1){
   if(mcstat==2)
    ksi0=1.-gsl_ran_poisson_pdf(0,mu0);
	
   gsl_vector *ss=gsl_vector_alloc (4), *x=gsl_vector_alloc (4);
   gsl_vector_set_all (ss, 0.1);
   gsl_vector_set(x,0,ksi0);
   gsl_vector_set(x,1,sigma0);
   gsl_vector_set(x,2,1.);
   gsl_vector_set(x,3,1.);
   minex_func.f = &lambda1;
   minex_func.n = 4;
   minex_func.params = NULL;
   s=NULL; 
   s = gsl_multimin_fminimizer_alloc (T, 4);
   gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
   iter=0;
   do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status){break;}
    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, 1.e-6);
    if(test_flag){
 	 if (status == GSL_SUCCESS){
	  fprintf (stdout,"converged to minimum at\n");}
	 fprintf (stdout,"%5d ", (int)iter);
	 for (k = 0; k < 4; k++){
	  fprintf (stdout,"%10.3e ", gsl_vector_get (s->x, k));}
	 fprintf (stdout,"f = %e\n", s->fval);
    }
   }
   while (status == GSL_CONTINUE && iter < 1000);
   ksiest=gsl_vector_get (s->x, 0);
   sigmaest=gsl_vector_get (s->x, 1);
   sigma1est=gsl_vector_get (s->x, 2);
   sigma2est=gsl_vector_get (s->x, 3);
   varest=sigmaest*sigmaest;
   varmeanest=ksiest*varest;
   ksiestMean=ksiestMean+ksiest;
   if(verbose_flag)
    fprintf(stdout,"ksiest=%e sigmaest=%e varest=%e varmeanest=%e sigma1est=%e sigma2est=%e\n",
	                ksiest,sigmaest,varest,varmeanest,sigma1est,sigma2est);
   fprintf(pf3,"%e %e %e %e %e\n",ksiest,sigmaest,varmeanest,sigma1est,sigma2est);
   
   if(stat_flag){ 
    gsl_vector_set(x,0,ksiest);
	gsl_vector_set(x,1,sigmaest);
	gsl_vector_set(x,2,sigma1est);
	gsl_vector_set(x,3,sigma2est);
	cp0=lambda1(x,NULL);
    fprintf(stdout,"likelihood function at estimated parameters: %e\n",-(double)Npt*cp0);
    gsl_vector_set(x,0,ksi);
	gsl_vector_set(x,1,sigma);
	gsl_vector_set(x,2,sigma1);
	gsl_vector_set(x,3,sigma2);
	cp0=lambda1(x,NULL);
	printf("likelihood function at injected parameters: %e\n",-(double)Npt*cp0);
   }
   
   //test for comparison with fig 11 of Drasc's paper
   if(contour_flag){
    pf4=fopen("contour.dat","w");
    for(m=1;m<500;m++){
	 for(l=1;l<500;l++){
	  ksi_ml=(double)m*0.001;
	  sigma_ml=sqrt((double)l*0.001);
	  gsl_vector_set(x,0,ksi_ml);
      gsl_vector_set(x,1,sigma_ml);
      gsl_vector_set(x,2,1.);
      gsl_vector_set(x,3,1.);
	  CP[m][l]=-(double)Npt*lambda1(x,NULL);
	  }}
	for(m=1;m<500;m++)
	 for(l=1;l<500;l++){
	  if(CP[m][l]>50.)
	  fprintf(pf4,"{%f,%f,%f},",(double)m*0.001,(double)l*0.001,CP[m][l]);}
	  
    fclose(pf4);
   }
  gsl_multimin_fminimizer_free (s);
  gsl_vector_free(x); 
  gsl_vector_free(ss); 	
  }
  
  else if(stat==0){
   gsl_vector *ss=gsl_vector_alloc (3), *x=gsl_vector_alloc (3);
   gsl_vector_set_all (ss, 0.1);
   gsl_vector_set(x,0,sigma0);
   gsl_vector_set(x,1,1.);
   gsl_vector_set(x,2,1.);

   minex_func.f = &lambda0;
   minex_func.n = 3;
   minex_func.params = NULL;
  
   s=NULL; 
   s = gsl_multimin_fminimizer_alloc (T, 3);
   gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

   iter=0;
   do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if (status){break;}
    size = gsl_multimin_fminimizer_size (s);
    status = gsl_multimin_test_size (size, 1.e-4);
    if(test_flag){
 	 if (status == GSL_SUCCESS){
	  fprintf (stdout,"converged to minimum at\n");}
	 fprintf (stdout,"%5d ", (int)iter);
	 for (k = 0; k < 3; k++){
	  fprintf (stdout,"%10.3e ", gsl_vector_get (s->x, k));}
	 fprintf (stdout,"f = %e\n", s->fval);
    }
   }
   while (status == GSL_CONTINUE && iter < 500);
   sigmaest=gsl_vector_get (s->x, 0);
   sigma1est=gsl_vector_get (s->x, 1);
   sigma2est=gsl_vector_get (s->x, 2);
   varmeanest=sigmaest*sigmaest;
   if(verbose_flag)
    fprintf(stdout,"sigmaest=%e varmeanest=%e sigma1est=%e sigma2est=%e\n",
	                sigmaest,varmeanest,sigma1est,sigma2est);
   fprintf(pf3,"%e %e %e %e\n",sigmaest,varmeanest,sigma1est,sigma2est); 
   
   if(stat_flag){
	gsl_vector_set(x,0,sigmaest);
	gsl_vector_set(x,1,sigma1est);
	gsl_vector_set(x,2,sigma2est);
	cp0=lambda0(x,NULL);
	printf("likelihood function at estimated parameters: %e\n",-(double)Npt*cp0);
	gsl_vector_set(x,0,sigma);
	gsl_vector_set(x,1,sigma1);
	gsl_vector_set(x,2,sigma2);
	cp0=lambda0(x,NULL);
	printf("likelihood function at injected parameters: %e\n",-(double)Npt*cp0);
   }
  gsl_multimin_fminimizer_free (s);
  gsl_vector_free(x); 
  gsl_vector_free(ss); 	
 
  }
 
  
 LALDDestroyVector(&lalstatus,&h);
 LALDDestroyVector(&lalstatus,&n1.data);
 LALDDestroyVector(&lalstatus,&n2.data);
 LALDDestroyVector(&lalstatus,&s1);
 LALDDestroyVector(&lalstatus,&s2);
 LALDDestroyVector(&lalstatus,&e12);
 if(stat_flag)
  LALDDestroyVector(&lalstatus,&hgw); 
 
 sigmaestMean=sigmaestMean+sigmaest;
 varmeanestMean=varmeanestMean+varmeanest;
 sigma1estMean=sigma1estMean+sigma1est;
 sigma2estMean=sigma2estMean+sigma2est;
 }  
 
 /** postprocessing **/
 if(post_flag){
  if(verbose_flag){
   fprintf(stdout,"combine the results of the %d segments\n",nsegment);

   if(stat==1)
    fprintf(stdout,"ksiest=%e sigmaest=%e varmeanest=%e sigma1est=%e sigma2est=%e\n",ksiestMean/nsegment,sigmaestMean/nsegment,varmeanestMean/nsegment,sigma1estMean/nsegment,sigma2estMean/nsegment);
   else if(stat==2)
    fprintf(stdout,"muest=%e sigmaest=%e varmeanest=%e sigma1est=%e sigma2est=%e\n",muestMean/nsegment,sigmaestMean/nsegment,varmeanestMean/nsegment,sigma1estMean/nsegment,sigma2estMean/nsegment);
   else	
	fprintf(stdout,"sigmaest=%e varmeanest=%e sigma1est=%e sigma2est=%e\n",sigmaestMean/nsegment,varmeanestMean/nsegment,sigma1estMean/nsegment,sigma2estMean/nsegment);
	}
  }
  

 /** cleanup **/
 if(verbose_flag){
  fprintf(stdout,"cleanup and exit\n");}

 gsl_rng_free (rrgn);
 if(ascii_flag){ 
   fclose(pf1); fclose(pf2);}fclose(pf3);
 
 if(frame_flag){
 /* close frame cache */
 LALFrClose(&lalstatus, &frStream1);
 LALFrClose(&lalstatus, &frStream2);
 }
 
 return 0;
}

double lambda0 (gsl_vector *x,void *params)
{
  double y;
  double psigma,psigma1,psigma2;
  double pv,pv1,pv2,psig12,pv12,pvv1,pvv2;
  psigma=gsl_vector_get(x,0);
  psigma1=gsl_vector_get(x,1);
  psigma2=gsl_vector_get(x,2);

  if((psigma>0.)&&(psigma1>0.)&&(psigma2>0.)){
   pv=psigma*psigma;
   pv1=psigma1*psigma1;
   pv2=psigma2*psigma2;
   psig12=psigma1*psigma2;
   pv12=psig12*psig12;
   pvv1=pv*pv1;
   pvv2=pv*pv2;

   y=-(0.5*(log(v1*v2)-log(pv12+pvv1+pvv2))+(v1/(pv1*pv1)+v2/(pv2*pv2)+2.*v12/pv12)/(2.*(1./pv1+1./pv2+1./pv))-v1/(2.*pv1)-v2/(2.*pv2)+1.);
  }
  else y=10000.;
  return y; 

}

double lambda1 (gsl_vector *x,void *params)
{
  int i;
  double y;
  double dsum,psum;
  double pksi,psigma,psigma1,psigma2;
  double pv,pv1,pv2,psig12,pv12,pvv1,pvv2,v1opv1,v2opv2;
  double a,b;  

  pksi=gsl_vector_get(x,0);
  psigma=gsl_vector_get(x,1);
  psigma1=gsl_vector_get(x,2);
  psigma2=gsl_vector_get(x,3);

  if((pksi>0.)&&(pksi<1.)&&(psigma>0.)&&(psigma1>0.)&&(psigma2>0.)){
   pv=psigma*psigma;
   pv1=psigma1*psigma1;
   pv2=psigma2*psigma2;
   psig12=psigma1*psigma2;
   pv12=psig12*psig12;
   pvv1=pv*pv1;
   pvv2=pv*pv2;
   v1opv1=v1/pv1;
   v2opv2=v2/pv2;

   a=1./sqrt(pv12+pvv1+pvv2);
   b=2.*(1./pv1+1./pv2+1./pv);
   
   dsum=0.;
   for(i=0;i<Npt;i++){
     e12->data[i]=(s1->data[i]/pv1+s2->data[i]/pv2);
     e12->data[i]= e12->data[i]*e12->data[i];
     psum = (1.-pksi)+pksi*psig12*a*exp(e12->data[i]/b);
     dsum=dsum+log(psum);
   }
   y=-(0.5*(log(v1opv1*v2opv2)-(v1opv1+v2opv2))+1.+(1./Npt)*dsum);
  }
  else y=10000.;
  return y; 
}

double lambda2 (gsl_vector *x,void *params)
{

  int i,j;
  double y;
  double dsum,psum,cfd;
  double pksi,pmu,psigma,psigma1,psigma2;
  double pv,pv1,pv2,psig12,pv12,pvv1,pvv2,v1opv1,v2opv2;
  double pc[10],a[10],b[10]; 
  
  pksi=gsl_vector_get(x,0);
  psigma=gsl_vector_get(x,1);
  psigma1=gsl_vector_get(x,2);
  psigma2=gsl_vector_get(x,3);
 
  if((pksi>0.)&&(pksi<1.)&&(psigma>0.)&&(psigma1>0.)&&(psigma2>0.)){

   pv=psigma*psigma;
   pv1=psigma1*psigma1;
   pv2=psigma2*psigma2;
   psig12=psigma1*psigma2;
   pv12=psig12*psig12;
   pvv1=pv*pv1;
   pvv2=pv*pv2;
   v1opv1=v1/pv1;
   v2opv2=v2/pv2;   
   
   pmu=-log(1.-pksi);
   pc[0]=1.-pksi;
   cfd=pksi;
   psum=0.;

   for(j=1;j<nmax;j++){
    pc[j]=pc[j-1]*pmu/(double)j;
	cfd=cfd-pc[j];
	a[j]=1./sqrt(pv12+(double)j*(pvv1+pvv2))*pc[j];
	b[j]=2.*(1./pv1+1./pv2+1./((double)j*pv));
   }
   a[nmax]=1./sqrt(pv12+(double)nmax*(pvv1+pvv2))*cfd;
   b[nmax]=2.*(1./pv1+1./pv2+1./(nmax*pv));
   
   dsum=0.;
   for(i=0;i<Npt;i++){
     e12->data[i]=(s1->data[i]/pv1+s2->data[i]/pv2);
     e12->data[i]= e12->data[i]*e12->data[i];
	 psum=0.;
	 for(j=1;j<=nmax;j++){
      psum=psum+a[j]*exp(e12->data[i]/b[j]);
	  }
	 psum=(1.-pksi)+psig12*psum;
	 dsum=dsum+log(psum);
   }
   y=-(0.5*(log(v1opv1*v2opv2)-(v1opv1+v2opv2))+1.+(1./(double)Npt)*dsum);
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
	  {"normalize", no_argument, &normalize_flag, 1},
	  {"test", no_argument, &test_flag, 1},
	  {"montecarlo", no_argument, &montecarlo_flag, 1},
	  {"condor", no_argument, &condor_flag, 1},
	  {"post", no_argument, &post_flag, 1},
	  {"stat", no_argument, &stat_flag, 1},
	  {"contour", no_argument, &contour_flag, 1},
      /* options that don't set a flag */
      {"help", no_argument, 0, 'h'},
	  {"job-number", required_argument, 0, 'n'},
      {"gps-start-time", required_argument, 0, 't'},
	  {"gps-stop-time", required_argument, 0, 'T'},
	  {"duration", required_argument, 0, 'l'},
      {"number-points", required_argument, 0, 'N'},
	  {"sampleRate", required_argument, 0, 'r'},
      {"resampleRate", required_argument, 0, 'R'},
	  {"stat", required_argument, 0, 'a'},
	  {"mcstat", required_argument, 0, 'A'},
	  {"ksi", required_argument, 0, 'k'},
	  {"ksi0", required_argument, 0, 'K'},
      {"mu", required_argument, 0, 'm'},
	  {"mu0", required_argument, 0, 'M'},
      {"sigma", required_argument, 0, 's'},
	  {"sigma0", required_argument, 0, 'S'},
	  {"sigma1_ref", required_argument, 0, 'g'},
	  {"sigma2_ref", required_argument, 0, 'G'},
	  {"nmax", required_argument, 0, 'p'},
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
                  "hn:t:T:l:N:r:R:a:A:k:K:m:M:s:S:g:G:p:c:C:d:D:v:",
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
			   
      case 'n':
			   /* jub number */
	           job = atoi(optarg);
	           break;
			   
      case 't':
			   /* start time */
	           startTime = atoi(optarg);
	           break;
			   
	  case 'T':
			   /* stop time */
	           stopTime = atoi(optarg);
	           break;
			   
	  case 'l':
			   /* duration if condor flag */
	           duration = atoi(optarg);
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
			   
	  case 'a':
			   /* statistic for analysis */
	           stat = atoi(optarg);
	           break;
			   
	  case 'A':
			   /* statistic for MC */
	           mcstat = atoi(optarg);
	           break;
			   
	  case 'k':
	           ksi = atof(optarg);
	           break;
			   
	  case 'K':
	           ksi0 = atof(optarg);
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
	           sigma1_ref = atof(optarg);
	           break;
			   
	  case 'G':
	           sigma2_ref = atof(optarg);
	           break;
			   
	  case 'p':
	           nmax = atoi(optarg);
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
  fprintf(stderr, " --normalize           normalize noise data (not needed for gaussian noise)\n");
  fprintf(stderr, " --montecarlo          inject gw signal\n");
  fprintf(stderr, " --condor              run on cluster\n");
  fprintf(stderr, " --post                crude post processing\n");
  fprintf(stderr, " --contour             generate matrix for contour plot (when -a 1)\n");
  fprintf(stderr, " --stat                display mean and variance of the amplitude distribution of the GW signal\n");
  fprintf(stderr, " -n                    job number\n");
  fprintf(stderr, " -t                    GPS start time\n");
  fprintf(stderr, " -T                    GPS stop time\n");
  fprintf(stderr, " -l                    length of data set on each node\n");
  fprintf(stderr, " -N                    number of points per data segment\n");
  fprintf(stderr, " -r                    sampleRate\n");
  fprintf(stderr, " -R                    resampleRate\n");
  fprintf(stderr, " -k                    ksi of injected gw signal if option -a 1 -A 1\n");
  fprintf(stderr, " -K                    initial ksi for parameter estimation if -a 1 -A 1\n");
  fprintf(stderr, " -m                    mu of injected gw signal if -A 2\n");
  fprintf(stderr, " -M                    initial mu for parameter estimation if -a 2 -A 2\n");
  fprintf(stderr, " -s                    sigma of injected gw signal\n");
  fprintf(stderr, " -S                    initial sigma for parameter estimation\n");
  fprintf(stderr, " -g                    sigma1 of gaussian noise \n");
  fprintf(stderr, " -G                    sigma2 of gaussian noise \n");
  fprintf(stderr, " -a                    statistic for analysis: 0 for CC, 1 for ML (no overlap), 2 for ML (overlap)\n");
  fprintf(stderr, " -A                    statistic for MC simulations: 2 by default, 1 to test Drasco statistic, 0 to test CC analysis\n");
  fprintf(stderr, " -p                    maximal number of events for the Poisson distribution for the calculation of the likelihood\n");
  fprintf(stderr, " -c                    channel for first stream\n");
  fprintf(stderr, " -C                    channel for second stream\n");
  fprintf(stderr, " -d                    cache file for first stream\n");
  fprintf(stderr, " -D                    cache file for second stream\n");
  exit(exitcode);
}
