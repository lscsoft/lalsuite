/* TBary.c is a prototype program written to perform the calculation of the F-statistic via the use of the algorithm described in the JKS paper on pulsar data analysis. This program uses LAL, LALApps, GSL and FFTW. It also requires files which provide data on the location of the earth and sun during a desired interval of time.

Authors - Pinkesh Patel, Rejean Dupuis & Xavier Siemens
References - Jaranowski,Krolak and Schutz, Phys Rev D. Numerical Recipes in C, FFTW and FFTW3, GNU Scientific Library. 
*/

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <lal/DetResponse.h>
#include <time.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/LALInitBarycenter.h>
#include <lal/LALComputeAM.h>
#include <lal/LALMalloc.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <fftw3.h>
#include <time.h>
#include <string.h>

using namespace std;                      

/****** Ephemeris files, hard coded for now ******/
#define EARTHDATA "/archive/home/ppatel/opt/lscsoft/lal/share/lal/earth00-04.dat"    
#define SUNDATA "/archive/home/ppatel/opt/lscsoft/lal/share/lal/sun00-04.dat"

/****** T0 is the reference time in the Solar System Barycenter ******/     
#define T0  700000000;

/****** A structure specifying frequency derivatives, num -> number of derivatives, fk contains the values ******/
struct freq
{
  int num;
  double *fk;
};

/****** Given a GPS time, returns a double precision value of the time ******/
inline double gettime(LIGOTimeGPS t)
{
  return t.gpsSeconds+t.gpsNanoSeconds*1e-9;
}

/****** Given a double precision time, returns its equivalent GPS time ******/
inline LIGOTimeGPS gettime(double t)
{
  LIGOTimeGPS tgps;
  tgps.gpsSeconds = (INT4)floor(t);
  tgps.gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));
  return(tgps);
}

/****** Returns the magnitude square of a complex number ******/
inline double magsquare(fftw_complex f)
{
  return(f[0]*f[0]+f[1]*f[1]);
}

/****** Returns the inner product of two vectors x and y as defined by JKS in the form of <x||y> ******/
double innerprod(int n_steps,double dt,double *x,double *y)
{
  int i;
  double sum = 0;                                               
  double t = n_steps*dt;                                       
  for(i=0;i<n_steps-1;i++)
    {
      sum += dt/2*(x[i]*y[i]+x[i+1]*y[i+1]);                     
    }
  return(sum*2/t);                              
}

/***** Resamples using linear interpolation (detailed comments only in the spline interpolation function) ******/
void interp_lin(LIGOTimeGPS start,double dt,long int inpsteps,long int n_steps,LIGOTimeGPS *tb,double *real,double *imag,fftw_complex *out)
{
  int i;
  double j;
  int count = 1;
  double *t,t0 = gettime(start);
  t = (double *)malloc(sizeof(double)*inpsteps);
  for(j=0;j<inpsteps;j++)
    t[(int)(j)] = t0 + j*dt;
  gsl_interp_accel *accl = gsl_interp_accel_alloc();
  gsl_interp *lininter = gsl_interp_alloc(gsl_interp_linear,inpsteps);
  gsl_interp_init(lininter,t,real,inpsteps);
  double x,y;
  for(i=0;i<n_steps;i++)
    {
      x = gettime(tb[i]);
      //if(x>t[inpsteps-1])
      //cerr<<"extrapolating";
      out[i][0] = gsl_interp_eval(lininter,t,real,x,accl);
    }
  gsl_interp_free(lininter);
  gsl_interp_accel_free(accl);
  accl = gsl_interp_accel_alloc();
  lininter = gsl_interp_alloc(gsl_interp_linear,inpsteps);
  gsl_interp_init(lininter,t,imag,inpsteps);
  for(i=0;i<n_steps;i++)
    {
      x = gettime(tb[i]);
      //if(x>t[inpsteps-1])
      //cerr<<"extrapolating";
      out[i][1] = gsl_interp_eval(lininter,t,imag,x,accl);
    }
  gsl_interp_free(lininter);
  gsl_interp_accel_free(accl);
}

/***** Resamples using Akima interpolation (detailed comments only in the spline interpolation function) ******/
void interp_akima(LIGOTimeGPS start,double dt,long int inpsteps,long int n_steps,LIGOTimeGPS *tb,double *real,double *imag,fftw_complex *out)
{
  int i;
  double j;
  int count = 1;
  double *t,t0 = gettime(start);
  t = (double *)malloc(sizeof(double)*inpsteps);
  for(j=0;j<inpsteps;j++)
    t[(int)(j)] = t0 + j*dt;
  gsl_interp_accel *accl = gsl_interp_accel_alloc();
  gsl_interp *akimainter = gsl_interp_alloc(gsl_interp_akima,inpsteps);
  gsl_interp_init(akimainter,t,real,inpsteps);
  double x,y;
  for(i=0;i<n_steps;i++)
    {
      x = gettime(tb[i]);
      //if(x>t[inpsteps-1])
      //cerr<<"extrapolating";
      y = gsl_interp_eval(akimainter,t,real,x,accl);
      out[i][0] = y;
    }
  gsl_interp_free(akimainter);
  gsl_interp_accel_free(accl);
  accl = gsl_interp_accel_alloc();
  akimainter = gsl_interp_alloc(gsl_interp_akima,inpsteps);
  gsl_interp_init(akimainter,t,imag,inpsteps);
  for(i=0;i<n_steps;i++)
    {
      x = gettime(tb[i]);
      //if(x>t[inpsteps-1])
      //cerr<<"extrapolating";
      y = gsl_interp_eval(akimainter,t,imag,x,accl);
      out[i][1] = y;
    }
  gsl_interp_free(akimainter);
  gsl_interp_accel_free(accl);
}


/****** Resamples the input complex time series and outputs the resampled time series ******/
/******
        Inputs are (1) Start time of the input time series
                   (2) The sampling time dt
                   (3) The number of data points in the input time series -> inpsteps
                   (4) The number of data points required in the output time series -> n_steps < inpsteps 
		   (5) An array of times specifying the detector times which correspond to the t0+i*dt times (iterate over i)
		   (6) An array containing the real part of the input time series
		   (7) An array containing the imaginary part of the input time series
		   (8) The output in the from of one complex number 
		   
******/

void interp_spline(LIGOTimeGPS start,double dt,long int inpsteps,long int n_steps,LIGOTimeGPS *tb,double *real,double *imag,fftw_complex *out)
{
  /****** Book-keeping variables ******/
  int i;
  double j;
  int count = 1;

  /****** t  -> array of times corresponding to the times at which the input time series was sampled
	  t0 -> double precision of starting time ******/
  double *t,t0 = gettime(start);
  t = (double *)malloc(sizeof(double)*inpsteps);
  /****** t = t0 + dt * steps ******/
  for(j=0;j<inpsteps;j++)
    t[(int)(j)] = t0 + j*dt;
  
  /****** Interpolation declarations with x = t and y = real time series. ******/
  gsl_interp_accel *accl = gsl_interp_accel_alloc();
  gsl_spline *splineinter = gsl_spline_alloc(gsl_interp_cspline,inpsteps);
  gsl_spline_init(splineinter,t,real,inpsteps);

  /****** Now evaluate the interpolation function with the times noted in the barycentered time series. These times correspond to a linear sampling the barycentered frame starting at the barycentered time corresponding to t0 ******/
  double x,y;
  for(i=0;i<n_steps;i++)
    {
      x = gettime(tb[i]);
      //if(x>t[inpsteps-1])
      //cerr<<"extrapolating";
      y = gsl_spline_eval(splineinter,x,accl);
      out[i][0] = y;
    }

  /****** Now do the same thing as above for the imaginary part ******/
  gsl_spline_free(splineinter);
  gsl_interp_accel_free(accl);
  accl = gsl_interp_accel_alloc();
  splineinter = gsl_spline_alloc(gsl_interp_cspline,inpsteps);
  gsl_spline_init(splineinter,t,imag,inpsteps);
  for(i=0;i<n_steps;i++)
    {
      x = gettime(tb[i]);
      //if(x>t[inpsteps-1])
      //cerr<<"extrapolating";
      y = gsl_spline_eval(splineinter,x,accl);
      out[i][1] = y;
    }
  gsl_spline_free(splineinter);
  gsl_interp_accel_free(accl);

}

/* Random function, for debugging, will not be used in the main code */
void interp_spline_reverse(LIGOTimeGPS start,double dt,long int inpsteps,long int n_steps,LIGOTimeGPS *tb,double *real,double *imag,fftw_complex *out)
{
  int count = 1;
  double *t,t0 = gettime(start);
  t = (double *)malloc(sizeof(double)*inpsteps);
  for(int i=0;i<inpsteps;i++)
    t[i] = gettime(tb[i]);
  gsl_interp_accel *accl = gsl_interp_accel_alloc();
  gsl_spline *splineinter = gsl_spline_alloc(gsl_interp_cspline,inpsteps);
  gsl_spline_init(splineinter,t,real,inpsteps);
  double x,y;
  for(int i=0;i<n_steps;i++)
    {
      x = t0 + (double)(i)*dt;
      //cerr<<x<<" "<<t[i]<<endl;
      //if(x>t[inpsteps-1])
      //cerr<<"extrapolating";
      y = gsl_spline_eval(splineinter,x,accl);
      out[i][0] = y;
    }
  gsl_spline_free(splineinter);
  gsl_interp_accel_free(accl);
  accl = gsl_interp_accel_alloc();
  splineinter = gsl_spline_alloc(gsl_interp_cspline,inpsteps);
  gsl_spline_init(splineinter,t,imag,inpsteps);
  for(int i=0;i<n_steps;i++)
    {
      x = t0 + (double)(i)*dt;
      //if(x>t[inpsteps-1])
      //cerr<<"extrapolating";
      y = gsl_spline_eval(splineinter,x,accl);
      out[i][1] = y;
    }
}


/****** Returns the sinc function -> sin(pi*x)/x/pi ******/
double sinc(double t)
{
  if(t == 0)
    return 1.0;
  else
    return sin(M_PI*t)/t/M_PI;
}


/****** This function returns the interpolated valuegiven the two arrays x and t and other parameters like number of terms and starting time and dt ******/
void retband(double t0, double dt, double* t,double* x, double* y,int n,int size, int terms)
{
  int i,j;
  double f_s = 1/dt;
  double Value = 0;
  int lterms = terms;
  int rterms = terms;
  for(i=0;i<size;i++)
    {
      Value = 0;
      int index;
      index = (int)((t[i]-t0)/dt);
      if(index < terms)
	lterms = index;
      if(index > (n-terms))
	rterms = (n-index);
      for(j=0;j<lterms;j++)
	Value += sinc(f_s*(t[i]-t0-(index-j)*dt))*x[index-j];
      for(j=1;j<rterms;j++)
	Value += sinc(f_s*(t[i]-t0-(index+j)*dt))*x[index+j];
      y[i] = Value;
    }
}

/****** Performs the same function as the spline interpolation, but does this using band-limited interpolation ******/
void interp_bandl(LIGOTimeGPS start,double dt,long int inpsteps,long int n_steps,LIGOTimeGPS *tb,double *real,double *imag,fftw_complex *out)
{
  int terms = 100,i = 0;
  double t0 = gettime(start),*y,*t;
  t = (double*)malloc(sizeof(double)*n_steps);
  
  for(i=0;i<n_steps;i++)
    t[i] = gettime(tb[i]);
  
  y = (double *)malloc(sizeof(double)*n_steps);
  retband(t0,dt,t,real,y,inpsteps,n_steps,terms);
  for(i=0;i<n_steps;i++)
    out[i][0] = y[i];
  retband(t0,dt,t,imag,y,inpsteps,n_steps,terms);
  for(i=0;i<n_steps;i++)
    out[i][1] = y[i];
}


/****** This function calculates the values of a(t_b) and b(t_b), given the value of psi and the times at which the value is required ******/
/******
        The inputs to this function are
	(1) tdout, the times in the detector frame corresponding to a linear sampling the barycentric frame.
	(2) n_steps is the number of points in tdout.
	(3) a and b are the arrays which are output.
	(4) LALstatus is for debugging purposes.
	(5) LALDetAndSource encodes all the information about the detector and source location.
	(6) psi is required to calculate the a's and b's from F_+ and F_x.

******/

void getAandB(LIGOTimeGPS *tdout,long int n_steps,double* a,double* b, LALStatus status,LALDetAndSource both,double psi)
{
  /****** Book-keeping variables ******/
  LIGOTimeGPS tgps;
  double sin2psi = sin(2*psi);     
  double cos2psi = cos(2*psi);
  
  /****** atemp and btemp stores the vales of a's and b's at some coarse intervals, say 3 minutes here and then interpolates as this is faster ******/
  double *atemp,*btemp,*T;
  LALGPSandAcc pGPSandAcc;
  LALDetAMResponse response;
  pGPSandAcc.gps.gpsNanoSeconds = 0;

  /****** t0 is the start time and its the first time in the array of barycentered detector times ******/
  double t0 = gettime(tdout[0]);
  
  /****** tstop is the time to stop procuring the a's and b' ******/
  double tstop = gettime(tdout[n_steps-1]);
  
  /****** The number of seconds used to calculate atemp and btemp before interpolation ******/
  double coarse = 180.0;
  
  /****** The number of steps in the coarse calculation ******/
  long int steps = long((tstop-t0)/coarse)+1;
  
  /****** Allocating memory ******/
  atemp = (double *)malloc(sizeof(double)*steps);
  btemp = (double *)malloc(sizeof(double)*steps);
  T = (double *)malloc(sizeof(double)*steps);

  int i = 0;
  double t;
  
  for(t=t0;t<=tstop;t+=coarse)                               
    {
      /****** Convert to LIGOGPS so that LAL functions can use it ******/
      tgps = gettime(t);
      pGPSandAcc.gps = tgps;

      /****** Compute Response ******/
      LALComputeDetAMResponse(&status, &response, &both, &pGPSandAcc);
      if(status.statusCode) 
	{
	  printf("Unexpectedly got error code %d and message %s\n",
		 status.statusCode, status.statusDescription);
	}
      
      /****** Calculate the a's and b's from the plus and cross responses ******/
      atemp[i] = sin2psi*response.plus+cos2psi*response.cross;
      btemp[i] = cos2psi*response.plus-sin2psi*response.cross;
      T[i] = t;
      i++;
    }

  /****** GSL declarations to interpolate a from atemp and T ******/
  gsl_interp_accel *accl = gsl_interp_accel_alloc();
  gsl_spline *splineinter = gsl_spline_alloc(gsl_interp_cspline,steps);
  gsl_spline_init(splineinter,T,atemp,steps);
  double x,y;
  for(i=0;i<n_steps;i++)
    {
      /****** Calculate for a at each time in tdout ******/
      x = gettime(tdout[i]);
      y = gsl_spline_eval(splineinter,x,accl);
      a[i] = y;
    }
  gsl_spline_free(splineinter);
  gsl_interp_accel_free(accl);

  /****** Do the same as above for b ******/
  accl = gsl_interp_accel_alloc();
  splineinter = gsl_spline_alloc(gsl_interp_cspline,steps);
  gsl_spline_init(splineinter,T,btemp,steps);
  for(i=0;i<n_steps;i++)
    {
      x = gettime(tdout[i]);
      y = gsl_spline_eval(splineinter,x,accl);
      b[i] = y;
    }

  gsl_spline_free(splineinter);
  gsl_interp_accel_free(accl);
  free(atemp);
  free(btemp);
  free(T);

}

/****** This is the main function of barcentering, since it calculates the times at the detector, which correspond to a linear sampling at the Solar System Barycenter ******/

 /****** This function recieves the following arguments
	 (1) tdout, an array which will be filled with times which correspond to linear sampling in the detector frame at the rate of dt.
	 (2) Phi_m is the difference in between the times at the SSB and the detector. (Needed for heterodyne correction)
	 (3) start is the starting time in the detector frame.
	 (4) n_steps is the number of data points required in tdout.
	 (5) status is for debugging purposes
	 (6) baryinput conatins all the information of the detector including the alpha, delta of the source and detector location etc.
	 (7) edata contains information of the Earths position and velocity.
 ******/
void getTbary(LIGOTimeGPS*tdout,double *Phi_m,LIGOTimeGPS start,long int n_steps,double dt,LALStatus status,BarycenterInput baryinput,EphemerisData edata)
{
  /***** Start and Stop times ******/
  double tstart,tstop; 
  
  /****** Book-keeping variables ******/
  LIGOTimeGPS tgps;
  double *tb,*tin,tout;          
  EarthState earth;                   
  EmissionTime emit;                              
  int k=0,i;                                       
  double t;

  /****** Start and stop times ******/
  tstart = gettime(start);
  tstop = tstart + n_steps*dt;
  
  /****** The barycentric times t_b(t_d) are calculated only at coarse intervals and then interpolation is used to fill in the gaps. This coarse parameter determines the number of seconds in between such calculations ******/ 
  double coarse = 180.0;

  /****** Number of coarse steps to be taken ******/
  long int steps = long((tstop-tstart)/coarse)+1;

  /****** Allocating memory ******/
  tb = (double *)malloc(sizeof(double)*steps);
  tin = (double *)malloc(sizeof(double)*steps);


  for (t=tstart;t<tstop;t+=coarse)
    {
      /****** Convert to LIGOGPS and set as baryinput ******/
      tgps = gettime(t); 
      baryinput.tgps = tgps;

      /****** Calculate the earth location ******/
      LALBarycenterEarth(&status,&earth,&baryinput.tgps,&edata);         
      if(status.statusCode) 
       {
         printf("Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
       }

      /****** emit contains the t_b and Phi_m, so calculate it ******/
      LALBarycenter(&status,&emit,&baryinput,&earth);         
      if(status.statusCode) 
       {
         printf("Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
       }
      
      /****** Store the times of the barycentering tb -> t_b and tin -> t_d ******/
      tb[k] = gettime(emit.te);
      tin[k++] = t;
    }
  
  
  /******  Now that we have t_b(t_d) and t_d(t_b), we can now use interpolation to calculate the values of t_d which correspond to a linear sampling t_d. So here we declare some GSL internal variables *****/
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,steps);
  gsl_spline_init(spline,tb,tin,steps);

  /****** The start time in the barycentric frame is the time corresponding to the barycentric time of the detector start time ******/
  double x0 = tb[0];
  double x,y;
  
  for(i=0;i<n_steps;i++)
    {
      /****** Here x -> barycentric time and it proceeds in linear steps of dt ******/
      x = x0 + double(i)*dt;
      /****** y is the corresponding detector time ******/
      y = gsl_spline_eval(spline,x,acc);
      tdout[i] = gettime(y);
      /****** This is difference in the two times ******/
      Phi_m[i] = y - x;
    }

  free(tb);
  free(tin);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

/****** Calculates the factorial of an integer ******/
long int factorial(long int x)
{
  if(x==0 || x==1)
    return(1);
  long int prod = 1,i;
  for(i=2;i<=x;i++)
    prod = prod*i;
  return(prod);
}

/****** Returns the maximum of an array x with steps n_steps ******/
double max(double* x,int n_steps)
{
  double M = 0;
  int i;
  for(i=0;i<n_steps;i++)
    {
      if(x[i]>M)
	M = x[i];
    }
  return M;
}

/****** Returns the maximum of an array x with steps n_steps and also the index of the maximum value ******/
double max2(double*x,int n_steps,int &index)
{
  double M = 0;
  int i;
  index = 0;
  for(i=0;i<n_steps;i++)
    {
      if(x[i]>M)
	{
	  M = x[i];
	  index = i;
	}
    }
  return M;
}


/****** This will eventually become a function in CFS_v2, which will take in all sorts of parameters and an array full of times ******/
int main(int argc, char **argv)
{
  int i,j,k;

  /****** LAL Variables ******/
  LALStatus status;
  static LALDetector detector;
  static EphemerisData edata;
  static LALSource source;
  static LALDetAndSource both;
  EarthState earth;
  BarycenterInput baryinput;
  EmissionTime emit;

  /****** RA, DEC, Psi are as they sound :) and F_check is the frequency whose F_stat you want returned, fdot is the first spindown ******/
  REAL8 RA = 0;
  REAL8 DEC = 0; 
  REAL8 psi = 0.5;
  REAL8 F_check = 0;
  REAL8 fdot = 0;
  
  /****** The dt number of steps in the time you want the F_stat you want calculated ******/
  double dt = 1.0/4.0;                               
  long int n_steps = 200000;             
 
  /****** inpdt == dt, the dt of the input timeseries and the number of points in it ******/
  double inpdt = 1.0/4.0;              
  long int inpsteps = 220000;

  /****** FILE I/O variables ******/
  int fileopen = 0;
  ifstream file;
  ifstream freqfile;
  ofstream ofile;
  
  /****** Checking for arguments from command line ******/
  for(i=1;i<argc;i++)             
    {
      if(argv[i][0] == '-')     
	switch(argv[i][1])           
	  {
	  case 'a': 
	    RA = atof(argv[i+1]);
	    i++;
	    break;

	  case 'd':
	    DEC = atof(argv[i+1]);
	    i++;
	    break;
	    
	  case 'p':
	    psi = atof(argv[i+1]);
	    i++;
	    break;

	  case 't':
	    dt = atof(argv[i+1]);
	    i++;
	    break;
	    
	  case 'n':
	    inpsteps = atoi(argv[i+1]);
	    n_steps = inpsteps - (long)(1000*dt);
	    i++;
	    break;

	  case 'l':
	    fdot = atof(argv[i+1]);
	    i++;
	    break;
	    
	  case 'f':
	    F_check = atof(argv[i+1]);
	    i++;
	    break;
	    
	  case 'i':
	    file.open(argv[i+1]);
	    i++;
	    fileopen = 1;
	    break;

	  }
    }

  /****** Set status to NULL to begin with ******/
  status.statusPtr=NULL;                                    
  
  /****** Using a cached detector and so setting the ephemeris files and the detector ******/
  detector = lalCachedDetectors[LALDetectorIndexLHODIFF]; 
  edata.ephiles.earthEphemeris = EARTHDATA;               
  edata.ephiles.sunEphemeris = SUNDATA;                   
  
  /****** Initializing earth data set edata ******/
  LALInitBarycenter( &status, &edata);          
 
  if(status.statusCode) 
    {
      printf("Unexpectedly got error code %d and message %s\n",
	     status.statusCode, status.statusDescription);
      return 0;
    }
	
  /******If file not opened by command line arguments then defaults to "T" ******/
  if(!fileopen)
    file.open("T");

  /****** in1,in2 are arrays which will be FFTed at the end, y is the resampled time series ******/
  fftw_complex *in1,*in2,*y;
  
  /****** Read in real and imaginary parts in x_re and x_im respectively, Two arrays are required as gsl interpolation will not deal with complex numbers ******/
  double *x_re,*x_im;
  
  /****** Allocate memory to x_re and x_im ******/
  x_re = (double *)malloc(sizeof(double)*inpsteps);
  x_im = (double *)malloc(sizeof(double)*inpsteps);
  
  /****** Read in data ******/
  fori=0;i<inpsteps;i++)
    {
      /****** temp is to ignore the first variable, which is time stamp ******/
      double temp;
      file>>temp;
      file>>x_re[i];
      file>>x_im[i];
    }

  /****** A for loop to loop over parameters if we want, its just there as a marker for now, and its not a loop for now ******/
  for(k=0;k<1;k++)
    {
      /****** Declare a Frequency structure F, this can potentially be changed to loop over different frequencies ******/
      freq F;
      /****** Only single fdot specified here ******/
      F.num = 1;
      /****** Now store it ******/
      F.fk = (double*)malloc(sizeof(double));
      F.fk[0] = fdot;
      
      /****** Allocate memory to in1,in2 and y ******/
      in1 = (fftw_complex *)malloc(sizeof(fftw_complex)*n_steps);
      in2 = (fftw_complex *)malloc(sizeof(fftw_complex)*n_steps);
      y = (fftw_complex *)malloc(sizeof(fftw_complex)*n_steps);

      /****** Set the time of first measurement ******/
      LIGOTimeGPS t0;                                        
      t0.gpsSeconds = 700010000;                            
      t0.gpsNanoSeconds = 0;

      /****** Set Source details ******/
      source.equatorialCoords.latitude = DEC;                
      source.equatorialCoords.longitude = RA;                
      source.equatorialCoords.system=COORDINATESYSTEM_EQUATORIAL;
      source.orientation = psi;
      
      /****** set the det and source details ******/
      both.pSource = &source;                    
      both.pDetector = &detector;
      
      /****** Set baryinput (for LAL) ******/
      baryinput.site.location[0] = detector.location[0]/LAL_C_SI; 
      baryinput.site.location[1] = detector.location[1]/LAL_C_SI;
      baryinput.site.location[2] = detector.location[2]/LAL_C_SI;    
      baryinput.alpha = RA;
      baryinput.delta = DEC; 
      
      /****** tdout will store the detectors times corresponding to a linear sampling in SSB starting at t_b(t0) ******/
      LIGOTimeGPS *tdout; 
      /****** Phi_m stores the difference between SSB and detector times ******/
      double *Phi_m;
      
      /****** Set some book-keeping times ******/
      double t=gettime(t0),tb;                               
      double t0d=t;
      
      /****** Allocate memory to these variables ******/
      tdout = (LIGOTimeGPS *)malloc(sizeof(LIGOTimeGPS)*n_steps);
      Phi_m = (double *)malloc(sizeof(double)*n_steps);

      /****** Calculate tdout ******/
      getTbary(tdout,Phi_m,t0,n_steps,dt,status,baryinput,edata);

      //for(int i=0;i<n_steps;i++)
      // tdout[i] = gettime(t+dt*double(i));
     
      /****** Declare and allocate memory to detector responses a and b ******/
      double *a,*b;                                        
      a = (double *)malloc(sizeof(double)*inpsteps);
      b = (double *)malloc(sizeof(double)*inpsteps);
      /****** Now calculate them ******/
      getAandB(tdout,n_steps,a,b,status,both,psi);

      /****** Resample the time series and calculate y ******/
      interp_spline(t0,dt,inpsteps,n_steps,tdout,x_re,x_im,y);

      /****** JUNK (IGNORE) ******/
      /*for(double j=0;j<n_steps;j++)
	{
	int i = j;
	//cout<<j*dt<<" "<<x_re[i]<<" "<<x_im[i]<<endl;
	x_re[i] = y[i][0];
	x_im[i] = y[i][1];
	}
	for(int i = n_steps; i<inpsteps;i++)
	{
	x_re[i] = 0;
	x_im[i] = 0;
	}
	
	interp_spline_reverse(t0,dt,inpsteps,n_steps,tdout,x_re,x_im,y);
	for(int i=0;i<n_steps;i++)
	{
	cout<<(double)(i)*dt<<" "<<x_re[i]<<" "<<y[i][1]<<endl;
	}*/
      //for(int i=0;i<n_steps;i++)
      // {
      //   y[i][0] = x_re[i];
      //   y[i][1] = x_im[i];
      // }
      
      
      
      
      /****** Book-keeping stuff for heterodyne correction and fdot correction ******/
      double y1_re,y1_im,y2_re,y2_im;
      double phi,sinphi,cosphi;

      /****** The heterodyne frequency ******/
      double Het = 106.0;	  
      
      /****** This loop corrects for the heterodyne shift, such that the heterodyning is not done purely in the SSB frame ******/
      for(i=0;i<n_steps;i++)
	{
	  /****** retemp and imtemp store the real and imaginary parts temporarily, shift = Het*Phi*2*Pi are multiplied ******/
	  double retemp,imtemp,shift = fmod(Het*Phi_m[i],1.0)*2*LAL_PI;
	  retemp = y[i][0];
	  imtemp = y[i][1];
	  
	  /****** y = y*exp(i*shift) ******/
	  y[i][0] = retemp*cos(shift) - imtemp*sin(shift);
	  y[i][1] = retemp*sin(shift) + imtemp*cos(shift);

	}
      /****** This loop adds in the frequency spin downs and also multiplies with a's and b's ******/
      for(i=0;i<n_steps;i++)
	{
	  /****** Calculate DT as the time elapsed ******/
	  double DT = gettime(tdout[i])+Phi_m[i]-T0;
	  /****** Multiply with a and b ******/
	  y1_re = a[i]*y[i][0];
	  y2_re = b[i]*y[i][0];
	  y1_im = a[i]*y[i][1];
	  y2_im = b[i]*y[i][1];
	  
	  /****** phi is phase shift due to spindown ******/
	  phi = 0;
	  for(j=1;j<=1;j++)
	    phi += 2*M_PI*fmod(F.fk[j-1]*(pow(DT,(j+1))/factorial(j+1)),1.0);
	  
	  /****** Calculate sins and coss to avoid doing them twice ******/
	  sinphi = sin(phi);
	  cosphi = cos(phi);

	  /****** in1 = in1 *exp(i*phi) and same for in2 ******/
	  in1[i][0] = y1_re*cosphi+y1_im*sinphi;
	  in2[i][0] = y2_re*cosphi+y2_im*sinphi;
	  in1[i][1] = y1_im*cosphi-y1_re*sinphi;
	  in2[i][1] = y2_im*cosphi-y2_re*sinphi;
	}
      
      /****** Following are no longer needed ******/
      free(tdout);
      free(Phi_m);
      free(y);
      
      /****** out1 and out2 are the results of the FFT ******/
      fftw_complex *out1,*out2;
      /****** Scalc is the estimate of the power spectral density ******/
      double *Scalc;
      /****** FFTW plans to calculate out1,out2 ******/
      fftw_plan p1,p2;
      /****** Stores the Fstatistic value ******/
      double *Fstat;
      /****** Inner products as defined in JKS of a and b and combinations of them ******/
      double A,B,C,D;
      A = innerprod(n_steps,dt,a,a);
      B = innerprod(n_steps,dt,b,b);
      C = innerprod(n_steps,dt,a,b);
      D = A*B-C*C;
      /****** a and b no longer required ******/
      free(a);
      free(b);
      
      /****** Allocate memories ******/
      Fstat = (double *) fftw_malloc(sizeof(double)*(n_steps));
      out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(n_steps));
      out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(n_steps));
      Scalc = (double *) fftw_malloc(sizeof(double)*n_steps);

      /****** Declare the plans for fftw ******/
      p1 = fftw_plan_dft_1d(n_steps,in1,out1,FFTW_FORWARD,FFTW_ESTIMATE);
      p2 = fftw_plan_dft_1d(n_steps,in2,out2,FFTW_FORWARD,FFTW_ESTIMATE);
      
      /****** FFT!!!! ******/
      fftw_execute(p1);
      fftw_execute(p2);
      
      /****** Dont need the data any more!!!! ******/
      free(in1);
      free(in2);
      /****** Set Noise Spectral Density ******/
      for(i=0;i<n_steps;i++)
	Scalc[i] = dt;
      /****** Calculate the F-Statistic ******/
      for(i=0;i<(n_steps);i++)
	{
	  Fstat[i] =2*4*dt*dt*(B*magsquare(out1[i])+A*magsquare(out2[i])-2*C*(out1[i][0]*out2[i][0]+out1[i][1]*out2[i][1]))/D/n_steps/dt/Scalc[i];
	}

      

      /****** Output the F-Statistic, there are two parts to it, because of the way FFTW stores the positive and negative parts ******/
      for(double i=0;i<n_steps/2.0;i++)
      if(fabs(F_check-((double)i/(double)(n_steps)/dt+Het))< 1.0/n_steps/dt/2.0)
      printf("%g %g %6.10g %g %g\n",RA,F_check,(double)i/(double)(n_steps)/dt+Het,Fstat[(int)(i)],Fstat[(int)(i)]);
      
      for(double i=n_steps-1;i>n_steps/2.0;i--)
      if(fabs(F_check-((double)(i-n_steps)/(double)(n_steps)/dt+Het)) < 1.0/n_steps/dt/2.0) 
      printf("%g %g %6.10g %g %g\n",RA,F_check,(double)(i-n_steps)/(double)(n_steps)/dt+Het,Fstat[(int)(i)],Fstat[(int)(i)]);
      
      
      


      /*-------------------------------------------------*/
      
      double MAX;
      int index;
      MAX = max2(Fstat,n_steps,index);
      if(index<n_steps/2)
	printf("%6.10g %g %6.10g \n",Het,MAX,(double)index/(double)(n_steps)/dt+Het);
      else
	printf("%6.10g %g %6.10g \n",Het,MAX,(double)(index-n_steps)/(double)(n_steps)/dt+Het);
      
      int h_steps=100;
      double H[h_steps];
      double Hsum = 0;
      double Hmax = 25,Hmin = 0;
      double bsize = (Hmax-Hmin)/double(h_steps);
      double D_KS = 0;
      for(int i=0;i<h_steps;i++)
	{ 
	  for(int j=0;j<(n_steps/2+1);j++)
	    if(Fstat[j] >= double(i)*bsize && Fstat[j] < double(i+1)*bsize)
	      Hsum += 1;
	  H[i] = Hsum/(n_steps/2+1);
	  //cout<<(double(i)+0.5)*bsize<<" "<<H[i]<<endl;
	  double err = fabs(H[i] - (1-(exp(-double(i+1)*bsize/2))*(1+double(i+1)*bsize/2)));
	  if ( err > D_KS)
	    D_KS = err;
	}
      //cerr<<RA<<" "<<DEC<<" "<<D_KS<<endl;
      fftw_destroy_plan(p1);
      fftw_destroy_plan(p2);
      fftw_free(Fstat);
      fftw_free(out1);
      fftw_free(out2);
      fftw_free(S);
      fftw_free(Scalc);
      free(F.fk);
      file.close();
    }
  free(x_re);
  free(x_im);
  LALCheckMemoryLeaks();
}
