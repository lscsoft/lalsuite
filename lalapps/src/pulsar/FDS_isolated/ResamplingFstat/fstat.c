/*
*  Copyright (C) 2007 Pinkesh Patel
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

/* TBary.c is a prototype program written to perform the calculation of the F-statistic via the use of the algorithm described in the JKS paper on pulsar data analysis. This program uses LAL, LALApps, GSL and FFTW. It also requires files which provide data on the location of the earth and sun during a desired interval of time.

Authors - Pinkesh Patel,Xavier Siemens,Rejean Dupuis, Alan Weinstein
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

using namespace std;                       //Standard C++ line, in order to avoid annoying std::

#define EARTHDATA "earth00-04.dat"         //The data file containing Earth's location and velocity for the years 2000-2004.
#define SUNDATA "sun00-04.dat"             //The data file containing Sun's location and velocity for the years 2000-2004.

#define T0  700000000 //751680000 //714153733

// gettime(LIGOTimeGPS t) returns the GPS time in double format.

struct freq
{
  int num;
  double *fk;
};

inline double gettime(LIGOTimeGPS t)
{
  return t.gpsSeconds+t.gpsNanoSeconds*1e-9;
}

// gettime(double t) returns the corresponding time in GPS format.

inline LIGOTimeGPS gettime(double t)
{
  LIGOTimeGPS tgps;
  tgps.gpsSeconds = (INT4)floor(t);
  tgps.gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));
  return(tgps);
}

//mag(fftw_complex f) returns the magnitude squared of the complex number f.

double magsquare(fftw_complex f)
{
  return(f[0]*f[0]+f[1]*f[1]);
}




/*innerprod(int n_steps, double dt, double *x, double *y) return the inner product of x(t),y(t) using the trapezoidal method of integration
  The inner product between two variables is defined as 2/T Integrate[x(t)*y(t),{t,-T,T}] (written in mathematica format).*/

double innerprod(int n_steps,double dt,double *x,double *y)
{
  double sum = 0;                                                // Integral
  double t = n_steps*dt;                                         // Total time = steps * delta of time.
  for(int i=0;i<n_steps-1;i++)
    {
      sum += dt/2*(x[i]*y[i]+x[i+1]*y[i+1]);                     // Integral calculated using trapezoidal rule.
    }
  return(sum*2/t);                                               // Multiply with 2/(total time) to get the correct factors.
}


/* Stroboscopic Sampling function, Takes in start time, the delta time (1/freq) of both input data and output data, number of steps of each and the data streams themselves.*/
/* Basic idea of the function is that it will pick the point nearest to the correct time for output stream, from the input stream. The error is a function of Input freq and Output freq and roughly scales as per (Input Freq)/(Output freq) */
/* We can also use interpolation for this, I feel that interpolation will be better, since it will be very difficult and cumbersome to load lots of data */
  	
void strob(LIGOTimeGPS start,double xdt,double ydt,LIGOTimeGPS *tb,double *x,double *y,long int xsteps,long int ysteps)
{
  for(int i=0;i<ysteps;i++)
    {
      double dt = gettime(tb[i])-gettime(start);                 // Difference between points calculated in time. 
      int j = int(dt/xdt);                                       // int(dt/xdt), gives us the index of the point closest to our time.
      y[i] = x[j];                                               // assign the next value in the output stream. 
      double DT = dt - ydt*double(i);
      //cout<<i<<" "<<j<<" "<<DT<<endl;
    }
}

void interp(LIGOTimeGPS start,double dt,long int inpsteps,long int n_steps,LIGOTimeGPS *tb,double *real,double *imag,fftw_complex *out)
{
  int count = 1;
  double *t,t0 = gettime(start);
  t = (double *)malloc(sizeof(double)*inpsteps);
  for(double i=0;i<inpsteps;i++)
    t[int(i)] = t0 + i*dt;
  gsl_interp_accel *accl = gsl_interp_accel_alloc();
  gsl_spline *splineinter = gsl_spline_alloc(gsl_interp_cspline,inpsteps);
  gsl_spline_init(splineinter,t,real,inpsteps);
  double x,y;
  for(int i=0;i<n_steps;i++)
    {
      x = gettime(tb[i]);
      if(x>t[inpsteps-1])
	cerr<<"extrapolating";
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
      x = gettime(tb[i]);
      if(x>t[inpsteps-1])
	cerr<<"extrapolating";
      y = gsl_spline_eval(splineinter,x,accl);
      out[i][1] = y;
    }
}

void getAandB(LIGOTimeGPS *tdout,long int n_steps,double* a,double* b, LALStatus status,LALDetAndSource both,double psi)
{
  LIGOTimeGPS tgps;
  double sin2psi = sin(2*psi);     
  double cos2psi = cos(2*psi);      
  double *atemp,*btemp,*T;
  LALGPSandAcc pGPSandAcc;
  LALDetAMResponse response;
  pGPSandAcc.gps.gpsNanoSeconds = 0;
  double t0 = gettime(tdout[0]);
  double tstop = gettime(tdout[n_steps-1]);
  double coarse = 1800.0;
  long int steps = long((tstop-t0)/coarse)+1;
  atemp = (double *)malloc(sizeof(double)*steps);
  btemp = (double *)malloc(sizeof(double)*steps);
  T = (double *)malloc(sizeof(double)*steps);
  int i = 0;
  for(double t=t0;t<=tstop;t+=coarse)                               
    {
      tgps = gettime(t);
      pGPSandAcc.gps = tgps;
      LALComputeDetAMResponse(&status, &response, &both, &pGPSandAcc);
      if(status.statusCode) 
	{
	  printf("Unexpectedly got error code %d and message %s\n",
		 status.statusCode, status.statusDescription);
	}
      atemp[i] = sin2psi*response.plus+cos2psi*response.cross;
      btemp[i] = cos2psi*response.plus-sin2psi*response.cross;
      T[i] = t;
      i++;
    }
  
  gsl_interp_accel *accl = gsl_interp_accel_alloc();
  gsl_interp *lininter = gsl_interp_alloc(gsl_interp_linear,steps);
  gsl_interp_init(lininter,T,atemp,steps);
  double x,y;
  for(int i=0;i<n_steps;i++)
    {
      x = gettime(tdout[i]);
      y = gsl_interp_eval(lininter,T,atemp,x,accl);
      a[i] = y;
    }
  gsl_interp_free(lininter);
  gsl_interp_accel_free(accl);
  accl = gsl_interp_accel_alloc();
  lininter = gsl_interp_alloc(gsl_interp_linear,steps);
  gsl_interp_init(lininter,T,atemp,steps);
  for(int i=0;i<n_steps;i++)
    {
      x = gettime(tdout[i]);
      y = gsl_interp_eval(lininter,T,btemp,x,accl);
      b[i] = y;
    }
  //for(int i=0;i<n_steps;i++)
  //cout<<t0+double(i)*dt<<" "<<a[i]<<" "<<b[i]<<endl;
  gsl_interp_free(lininter);
  gsl_interp_accel_free(accl);
  free(atemp);
  free(btemp);
  free(T);

}

double getTB0(double time,LALStatus status,BarycenterInput baryinput,EphemerisData edata)
{
  LIGOTimeGPS tgps = gettime(time);
  EarthState earth;
  EmissionTime emit;
  baryinput.tgps = tgps;
  LALBarycenterEarth(&status,&earth,&baryinput.tgps,&edata);           // Calculates the earth variable.
  if(status.statusCode) 
    {
      printf("Unexpectedly got error code %d and message %s\n",
	     status.statusCode, status.statusDescription);
    }
  LALBarycenter(&status,&emit,&baryinput,&earth);           // Calculates the variable emit, which contains an object called te, or t_b.
  if(status.statusCode) 
    {
      printf("Unexpectedly got error code %d and message %s\n",
	     status.statusCode, status.statusDescription);
    }
  return(gettime(emit.te));
}

/* getTbary takes in a start time, delta t and number of steps in the detector frame and converts it to barycentric time by using LALBarycenter. After which, it interpolates t_d vs t_b and returns values of t_d, which are linear in t_b.   */
  
void getTbary(LIGOTimeGPS*tdout,double *Phi_m,LIGOTimeGPS start,long int n_steps,double dt,LALStatus status,BarycenterInput baryinput,EphemerisData edata)
{
  double tstart,tstop; 
  LIGOTimeGPS tgps;                         
  double *tb,*tin,tout;                            // book-keeping variables; tb,tin are used in the interpolation.
  EarthState earth;                                // Required for LALBarycenter.
  EmissionTime emit;                               // Required for LALBarycenter.
  int k = 0;                                       // Another book-keeping variable
  
  tstart = gettime(start);
  tstop = tstart + n_steps*dt;
  double coarse = 1800.0;
  //cerr<<"start "<<tstart<<" stop "<<tstop<<endl;
  long int steps = long((tstop-tstart)/coarse)+1;
  tb = (double *)malloc(sizeof(double)*steps);
  tin = (double *)malloc(sizeof(double)*steps);
  for (double t=tstart;t<tstop;t+=coarse)
    {
      tgps = gettime(t); 
      baryinput.tgps = tgps;
      //cerr<<t<<endl;
      LALBarycenterEarth(&status,&earth,&baryinput.tgps,&edata);           // Calculates the earth variable.
      if(status.statusCode) 
       {
         printf("Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
       }
      LALBarycenter(&status,&emit,&baryinput,&earth);           // Calculates the variable emit, which contains an object called te, or t_b.
      if(status.statusCode) 
       {
         printf("Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
       }
      tb[k] = gettime(emit.te);
      tin[k++] = t;
      //cerr<<k<<endl;
    }
  
  /* Interpolation routine from gsl */

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline,steps);
  gsl_spline_init(spline,tb,tin,steps);
  double x0 = tb[0];
  double x,y;
  for(int i=0;i<n_steps;i++)
    {
      x = x0 + double(i)*dt;
      y = gsl_spline_eval(spline,x,acc);
      tdout[i] = gettime(y);
      Phi_m[i] = x - y;
      //cerr<<i<<endl;
      //cout<<i<<" "<<x<<" "<<tb[i]<<" "<<y<<" "<<y-x<<endl;
    }
  //cerr<<"done interpolating"<<endl;
  free(tb);
  free(tin);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

long int factorial(long int x)
{
  if(x==0 || x==1)
    return(1);
  long int prod = 1;
  for(long int i = 2;i<=x;i++)
    prod = prod*i;
  return(prod);
}


double max(double* x,int n_steps)
{
  double M = 0;
  for(int i=0;i<n_steps;i++)
    {
      if(x[i]>M)
	M = x[i];
    }
  return M;
}


int main(int argc, char **argv)
{
  int count = 1; //Counting number of templates run through

  cerr<<"5 "<<factorial(5)<<endl;
  int DetIndex = 0;
  LALStatus status;
  static LALDetector detector;
  static EphemerisData edata;
  static LALSource source;
  static LALDetAndSource both;
  EarthState earth;
  BarycenterInput baryinput;
  EmissionTime emit;

  
  REAL8 lat = 0;                                        // Latitude
  REAL8 lon = 0;                                        // Longitude
  REAL8 RA = 0;//4.88671*M_PI/180.0;                              // Right Ascension
  REAL8 DEC = 0;//-0.217584*M_PI/180.0;                             // Declination  
  REAL8 psi = 0;//LAL_PI/2.0;//0;//-0.647939;                                  // Psi - Source Orientation
  
  double dt = 1.0/9.0;                                // dt = time between consecutive data points in used data
  long int n_steps = 16000*9.0;                                  // Number of data points in used data
 
  double inpdt = 1.0/9.0;                       // The time between data points in input data
  long int inpsteps = 162001*9.0/10.0;            // Number of data points in input data
  ifstream file;
  ifstream freqfile;
  ofstream ofile;
  ofstream maxofile("max");
  ofstream binofile("bin");
  
  for(int i=1;i<argc;i++)                                // A loop over command line data  
    {
      if(argv[i][0] == '-')                              // '-' is required before input is accepted
	switch(argv[i][1])                               // a character after that will put in a range of options
	  {
	  case 'n': 
	    n_steps = atoi(argv[i+1]);
	    i++;
	    break;

	  case 'd':
	    dt = atof(argv[i+1]);
	    i++;
	    break;
	    
	  case 'o':
	    file.open(argv[i+1]);
	    i++;
	    break;
	  }
    }

  status.statusPtr=NULL;                                  // status is used to monitor errors while calling LAL routines.  
  
  detector = lalCachedDetectors[LALDetectorIndexLHODIFF]; // lalCachedDetectors[] is a precalculated detector.
  edata.ephiles.earthEphemeris = EARTHDATA;               // setting up ephemeris earth data.
  edata.ephiles.sunEphemeris = SUNDATA;                   // setting up ephemeris sun data.
  edata.leap = 13;                                        // Number of leap seconds adjusted.
  LALInitBarycenter( &status, &edata);                    // Initializing ephemeris data.
 
  if(status.statusCode) 
    {
      printf("Unexpectedly got error code %d and message %s\n",
	     status.statusCode, status.statusDescription);
      return 0;
    }
	
  
  freqfile.open("freq1");
  freq Fcorrect;
  double *tempfk;
  int tempint=0;
  tempfk = (double *)malloc(sizeof(double)*15);
  while(freqfile)
    {
      freqfile>>tempfk[tempint++];
    }
  Fcorrect.num = tempint;
  Fcorrect.fk = (double *)malloc(sizeof(double)*tempint);
  for(int i=0;i<tempint;i++)
    {
      Fcorrect.fk[i] = tempfk[i];
    }
  free(tempfk);

  ofile.open("output");
  ofstream ofileX("X"),ofileY("Y"),ofileZ("Z");
  int doonce = 1;
  
  file.open("plot");
  //file.open("data1",ios::binary);


  fftw_complex *in1,*in2,*y;
  double *x_re,*x_im;
  x_re = (double *)malloc(sizeof(double)*inpsteps);
  x_im = (double *)malloc(sizeof(double)*inpsteps);
  
  for(int i=0;i<inpsteps;i++)
    {
      double temp;
      file>>temp;
      file>>x_re[i];
      //x_im[i] = 0;
      file>>x_im[i];
      //file>>temp;
      //x_im[i] = 0;
    } 
  
  /* for(int i=0;i<n_steps;i++)
    {
      double retemp,imtemp,shift = 0,DT = double(i)*dt+700000000-T0;
      retemp = x_re[i];
      imtemp = x_im[i];
      
      x_re[i] = retemp*cos(2*LAL_PI*fmod(shift*DT,1.0)) - imtemp*sin(2*LAL_PI*fmod(shift*DT,1.0));
      x_im[i] = retemp*sin(2*LAL_PI*fmod(shift*DT,1.0)) + imtemp*cos(2*LAL_PI*fmod(shift*DT,1.0));
      //cout<<x_re[i]<<" "<<retemp<<" "<<x_im[i]<<" "<<imtemp<<endl;
    }
  */
  string si = "data", so = "output",s;
  for(int i=0;i<1;i++)
    for(int j=0;j<1;j++)
      for(int k=1;k<2;k++)
	{
	  freq F;
	  F.num = 1;
	  F.fk = (double*)malloc(sizeof(double));
	  double m = 0.9/50;
	  double c = 1/10.0;
	  F.fk[0] = 0;//-2.54e-8 ;//* (m*double(k)+c);
	  ofileX<<F.fk[0]<<endl;
	  
	  in1 = (fftw_complex *)malloc(sizeof(fftw_complex)*n_steps);
	  in2 = (fftw_complex *)malloc(sizeof(fftw_complex)*n_steps);
	  
	  y = (fftw_complex *)malloc(sizeof(fftw_complex)*n_steps);
	  char stemp[50];
	  sprintf(stemp,"%d",k);
	  s = si + stemp;
	  
	  
	  s = so + stemp;
	  

	 
	  
	  
	  
	  //inpsteps = INPSTEPS/k;
	  //n_steps = N_STEPS/k;
	  
	  


	  cerr<<s<<" "<<inpsteps<<" "<<n_steps<<endl;
	  
	  REAL8 RA = 0;//4.88671*M_PI/180.0;                           
	  DEC = 0;//-0.217584*M_PI/180.0; 
	  LIGOTimeGPS t0;                                        
	  t0.gpsSeconds = 700000000;                            
	  t0.gpsNanoSeconds = 0;
	  source.equatorialCoords.latitude = DEC;                // Setting up Source settings.
	  source.equatorialCoords.longitude = RA;                
	  source.equatorialCoords.system=COORDINATESYSTEM_EQUATORIAL;
	  source.orientation = psi;
	    
	  both.pSource = &source;                                // Det and Source variable both is initialized.
	  both.pDetector = &detector;
	  
	  baryinput.site.location[0] = detector.location[0]/LAL_C_SI; // Initializing baryinput.
	  baryinput.site.location[1] = detector.location[1]/LAL_C_SI;
	  baryinput.site.location[2] = detector.location[2]/LAL_C_SI;    
	  
	  baryinput.alpha = RA;
	  baryinput.delta = DEC; 
	  
	  LIGOTimeGPS *tdout; 
	  double *Phi_m;
	  double t=gettime(t0),tb;                               
	  double t0d=t;
	    
	  time_t start,end;
	  double dif;
	  tdout = (LIGOTimeGPS *)malloc(sizeof(LIGOTimeGPS)*n_steps);
	  Phi_m = (double *)malloc(sizeof(double)*n_steps);
	  getTbary(tdout,Phi_m,t0,n_steps,dt,status,baryinput,edata);
	  
	  //for(int i=0;i<n_steps;i++)
	  // tdout[i] = gettime(t+dt*double(i));
	  
	  cerr<<count++<<" Here\n"<<endl;
	  //for(int i=0;i<inpsteps;i++)
	  // cout<<x_re[i]<<endl;
	  
	  //strob(t0,inpdt,dt,tdout,x,y,inpsteps,n_steps);
	  double *a,*b;                                        
	  a = (double *)malloc(sizeof(double)*inpsteps);
	  b = (double *)malloc(sizeof(double)*inpsteps);

	  getAandB(tdout,n_steps,a,b,status,both,psi);
	  //for(int i =0;i<n_steps;i++)
	  // {
	  //  printf("%15.15g %g %g\n",gettime(tdout[i])-Phi_m[i],a[i],b[i]);
	  // }
	  interp(t0,dt,inpsteps,n_steps,tdout,x_re,x_im,y);

	  
	 
	  

	  double y1_re,y1_im,y2_re,y2_im;
	  double phi,sinphi,cosphi;
	  //double TB0 = getTB0(T0,status,baryinput,edata);
	  //cerr<<T0<<" "<<TB0<<endl;
	  for(int i=0;i<n_steps;i++)
	    {
	      double retemp,imtemp,shift = -5.0 ,DT = Phi_m[i];
	      retemp = y[i][0];
	      imtemp = y[i][1];
	      
	      y[i][0] = retemp*cos(2*LAL_PI*fmod(shift*DT,1.0)) - imtemp*sin(2*LAL_PI*fmod(shift*DT,1.0));
	      y[i][1] = retemp*sin(2*LAL_PI*fmod(shift*DT,1.0)) + imtemp*cos(2*LAL_PI*fmod(shift*DT,1.0));
	      //cout<<x_re[i]<<" "<<retemp<<" "<<x_im[i]<<" "<<imtemp<<endl;
	    }
	  
	  for(int i=0;i<n_steps;i++)
	    {
	      double DT = gettime(tdout[i])+Phi_m[i]-T0;
	      y1_re = a[i]*y[i][0];
	      y2_re = b[i]*y[i][0];
	      y1_im = a[i]*y[i][1];
	      y2_im = b[i]*y[i][1];
	      phi = 0;
	      for(int j=1;j<=1;j++)
	      {
	        phi += 2*M_PI*fmod(F.fk[j-1]*(pow(DT,(j+1))/factorial(j+1)),1.0);
		//if(phi != 0)
		//cerr<<"\n UMMM Problem "<<phi<<endl;
		//cout<<phi<<" "<<gettime(tdout[i])-T0<<endl;
	      }
	      sinphi = sin(phi);
	      cosphi = cos(phi);
	      in1[i][0] = y1_re*cosphi+y1_im*sinphi;
	      in2[i][0] = y2_re*cosphi+y2_im*sinphi;
	      in1[i][1] = y1_im*cosphi-y1_re*sinphi;
	      in2[i][1] = y2_im*cosphi-y2_re*sinphi;
	    }
	  
	  free(tdout);
	  free(Phi_m);
	  free(y);
	  
	  fftw_complex *out1,*out2,*S;
	  double *Scalc;
	  double Sfactor = n_steps/10;
	  fftw_plan p1,p2,p3;
	  double *Fstat;
	  double A,B,C,D;
	  double M1 = 0, M2 = 0, M3 = 0;
	  double *Noise;
	  
	  A = innerprod(n_steps,dt,a,a);
	  B = innerprod(n_steps,dt,b,b);
	  C = innerprod(n_steps,dt,a,b);
	  D = A*B-C*C;
	  free(a);
	  free(b);
	  
	  Fstat = (double *) fftw_malloc(sizeof(double)*(n_steps));
	  out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(n_steps));
	  out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(n_steps));
	  S = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*long(n_steps/Sfactor));
	  
	  p1 = fftw_plan_dft_1d(n_steps,in1,out1,FFTW_FORWARD,FFTW_ESTIMATE);
	  p2 = fftw_plan_dft_1d(n_steps,in2,out2,FFTW_FORWARD,FFTW_ESTIMATE);
	  
	  fftw_execute(p1);
	  fftw_execute(p2);
	  
	  free(in1);
	  free(in2);
	  
	  Scalc = (double *) fftw_malloc(sizeof(double)*long(n_steps/Sfactor));
	  
	  for(int i=0;i<n_steps/Sfactor;i++)
	    Scalc[i] = 2.0*dt;	    
	  int check = 1;
	  
	  
	  /*for(int i=0;i<Sfactor;i++)
	    {
	    Noise = (y + long(i*n_steps/Sfactor));
	    p3 = fftw_plan_dft_r2c_1d(long(n_steps/Sfactor),Noise,S,FFTW_ESTIMATE);
	    fftw_execute(p3);
	    for(int j=0;j<n_steps/Sfactor/2+1;j++)
	    Scalc[j] += magsquare(S[j])*(2.0/(n_steps*dt));
	    }*/
	  
	  for(int i=0;i<(n_steps);i++)
	    {
	      Fstat[i] =2*4*dt*dt*(B*magsquare(out1[i])+A*magsquare(out2[i])-2*C*(out1[i][0]*out2[i][0]+out1[i][1]*out2[i][1]))/D/n_steps/dt/Scalc[long(i/Sfactor)];
	      //ofile<<F.fk[0]<<" "<<double(i)/n_steps/dt<<" "<<Fstat[i]<<endl;
	      //ofile<<Fstat[i]<<endl;
	      /*if(k==1)
		{
		  ofileY<<double(i)/n_steps/dt<<endl;
		}
		ofileZ<<Fstat[i]<<" ";*/
	      if(double(i)/n_steps/dt == 0.02)
		binofile<<Fstat[i]<<endl;
	    }
	  //ofileZ<<endl;
	  for(int i=0;i<n_steps*2.0/8.0;i++)
	    ofile<<(double)i/(double)n_steps*9.0+5.0<<" "<<Fstat[i]<<endl;
	  for(int i=n_steps-1;i>n_steps*6.0/8.0;i--)
	    ofile<<5.0+(double)(i-n_steps)/(double)n_steps*9.0<<" "<<Fstat[i]<<endl;
	  
	  maxofile<<max(Fstat,n_steps)<<endl;
	  
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
	      double err = fabs(H[i] - (1-(exp(-double(i+1)*bsize/2))*(1+double(i+1)*bsize/2)));
	      if ( err > D_KS)
		D_KS = err;
	    }
	  cerr<<RA<<" "<<DEC<<" "<<D_KS<<endl;
	  fftw_destroy_plan(p1);
	  fftw_destroy_plan(p2);
	  fftw_free(Fstat);
	  fftw_free(out1);
	  fftw_free(out2);
	  fftw_free(S);
	  fftw_free(Scalc);
	  free(F.fk);
	  file.close();
	  //ofile.close();
	}
  free(x_re);
  free(x_im);
  LALCheckMemoryLeaks();
  //fftw_destroy_plan(p1);
  //fftw_destroy_plan(p2);
  //dif = end-start;
  //printf ("You have taken %.2lf seconds \n", dif );
}
