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

#define EARTHDATA "/archive/home/ppatel/opt/lscsoft/lal/share/lal/earth00-04.dat"    
#define SUNDATA "/archive/home/ppatel/opt/lscsoft/lal/share/lal/sun00-04.dat"     
#define T0  700000000;//751680000 

struct freq
{
  int num;
  double *fk;
};

/*int LALSnprintf( char *str, size_t size, const char *fmt, ... )
{
  int n;
  va_list ap;
  va_start( ap, fmt );
  n = vsnprintf( str, size, fmt, ap );
  va_end( ap );
  return n;
}
*/
inline double gettime(LIGOTimeGPS t)
{
  return t.gpsSeconds+t.gpsNanoSeconds*1e-9;
}

inline LIGOTimeGPS gettime(double t)
{
  LIGOTimeGPS tgps;
  tgps.gpsSeconds = (INT4)floor(t);
  tgps.gpsNanoSeconds = (INT4)floor((fmod(t,1.0)*1.e9));
  return(tgps);
}

double magsquare(fftw_complex f)
{
  return(f[0]*f[0]+f[1]*f[1]);
}

double innerprod(int n_steps,double dt,double *x,double *y)
{
  double sum = 0;                                               
  double t = n_steps*dt;                                       
  for(int i=0;i<n_steps-1;i++)
    {
      sum += dt/2*(x[i]*y[i]+x[i+1]*y[i+1]);                     
    }
  return(sum*2/t);                              
}


/********************************************************************/

void interp_lin(LIGOTimeGPS start,double dt,long int inpsteps,long int n_steps,LIGOTimeGPS *tb,double *real,double *imag,fftw_complex *out)
{
  int count = 1;
  double *t,t0 = gettime(start);
  t = (double *)malloc(sizeof(double)*inpsteps);
  for(double i=0;i<inpsteps;i++)
    t[int(i)] = t0 + i*dt;
  gsl_interp_accel *accl = gsl_interp_accel_alloc();
  gsl_interp *lininter = gsl_interp_alloc(gsl_interp_linear,inpsteps);
  gsl_interp_init(lininter,t,real,inpsteps);
  double x,y;
  for(int i=0;i<n_steps;i++)
    {
      x = gettime(tb[i]);
      if(x>t[inpsteps-1])
	cerr<<"extrapolating";
      out[i][0] = gsl_interp_eval(lininter,t,real,x,accl);
    }
  gsl_interp_free(lininter);
  gsl_interp_accel_free(accl);
  accl = gsl_interp_accel_alloc();
  lininter = gsl_interp_alloc(gsl_interp_linear,inpsteps);
  gsl_interp_init(lininter,t,imag,inpsteps);
  for(int i=0;i<n_steps;i++)
    {
      x = gettime(tb[i]);
      if(x>t[inpsteps-1])
	cerr<<"extrapolating";
      out[i][1] = gsl_interp_eval(lininter,t,imag,x,accl);
    }
  gsl_interp_free(lininter);
  gsl_interp_accel_free(accl);
}

/**************************************************************************/


void interp_akima(LIGOTimeGPS start,double dt,long int inpsteps,long int n_steps,LIGOTimeGPS *tb,double *real,double *imag,fftw_complex *out)
{
  int count = 1;
  double *t,t0 = gettime(start);
  t = (double *)malloc(sizeof(double)*inpsteps);
  for(double i=0;i<inpsteps;i++)
    t[int(i)] = t0 + i*dt;
  gsl_interp_accel *accl = gsl_interp_accel_alloc();
  gsl_interp *akimainter = gsl_interp_alloc(gsl_interp_akima,inpsteps);
  gsl_interp_init(akimainter,t,real,inpsteps);
  double x,y;
  for(int i=0;i<n_steps;i++)
    {
      x = gettime(tb[i]);
      if(x>t[inpsteps-1])
	cerr<<"extrapolating";
      y = gsl_interp_eval(akimainter,t,real,x,accl);
      out[i][0] = y;
    }
  gsl_interp_free(akimainter);
  gsl_interp_accel_free(accl);
  accl = gsl_interp_accel_alloc();
  akimainter = gsl_interp_alloc(gsl_interp_akima,inpsteps);
  gsl_interp_init(akimainter,t,imag,inpsteps);
  for(int i=0;i<n_steps;i++)
    {
      x = gettime(tb[i]);
      if(x>t[inpsteps-1])
	cerr<<"extrapolating";
      y = gsl_interp_eval(akimainter,t,imag,x,accl);
      out[i][1] = y;
    }
}


void interp_spline(LIGOTimeGPS start,double dt,long int inpsteps,long int n_steps,LIGOTimeGPS *tb,double *real,double *imag,fftw_complex *out)
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

double sinc(double t)
{
  if(t == 0)
    return 1.0;
  else
    return sin(M_PI*t)/t/M_PI;
}

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

void interp_bandl(LIGOTimeGPS start,double dt,long int inpsteps,long int n_steps,LIGOTimeGPS *tb,double *real,double *imag,fftw_complex *out)
{
  int terms = 15,i = 0;
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
  double coarse = 180.0;
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
  gsl_spline *splineinter = gsl_spline_alloc(gsl_interp_cspline,steps);
  gsl_spline_init(splineinter,T,atemp,steps);
  double x,y;
  for(int i=0;i<n_steps;i++)
    {
      x = gettime(tdout[i]);
      y = gsl_spline_eval(splineinter,x,accl);
      a[i] = y;
    }
  gsl_spline_free(splineinter);
  gsl_interp_accel_free(accl);
  accl = gsl_interp_accel_alloc();
  splineinter = gsl_spline_alloc(gsl_interp_cspline,steps);
  gsl_spline_init(splineinter,T,btemp,steps);
  for(int i=0;i<n_steps;i++)
    {
      x = gettime(tdout[i]);
      y = gsl_spline_eval(splineinter,x,accl);
      b[i] = y;
    }
  //for(int i=0;i<n_steps;i++)
  //cout<<t0+double(i)*dt<<" "<<a[i]<<" "<<b[i]<<endl;
  gsl_spline_free(splineinter);
  gsl_interp_accel_free(accl);
  free(atemp);
  free(btemp);
  free(T);

}

  
void getTbary(LIGOTimeGPS*tdout,double *Phi_m,LIGOTimeGPS start,long int n_steps,double dt,LALStatus status,BarycenterInput baryinput,EphemerisData edata)
{
  double tstart,tstop; 
  LIGOTimeGPS tgps;                         
  double *tb,*tin,tout;          
  EarthState earth;                   
  EmissionTime emit;                              
  int k = 0;                                       
  
  tstart = gettime(start);
  tstop = tstart + n_steps*dt;
  double coarse = 180.0;
  //cerr<<"start "<<tstart<<" stop "<<tstop<<endl;
  long int steps = long((tstop-tstart)/coarse)+1;
  tb = (double *)malloc(sizeof(double)*steps);
  tin = (double *)malloc(sizeof(double)*steps);
  for (double t=tstart;t<tstop;t+=coarse)
    {
      tgps = gettime(t); 
      baryinput.tgps = tgps;
      //cerr<<t<<endl;
      LALBarycenterEarth(&status,&earth,&baryinput.tgps,&edata);         
      if(status.statusCode) 
       {
         printf("Unexpectedly got error code %d and message %s\n",
	 status.statusCode, status.statusDescription);
       }
      LALBarycenter(&status,&emit,&baryinput,&earth);         
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
  //cerr<<"here \n";
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

double max2(double*x,int n_steps,int &index)
{
  double M = 0;
  index = 0;
  for(int i=0;i<n_steps;i++)
    {
      if(x[i]>M)
	{
	  M = x[i];
	  index = i;
	}
    }
  return M;
}

int main(int argc, char **argv)
{
  LALStatus status;
  static LALDetector detector;
  static EphemerisData edata;
  static LALSource source;
  static LALDetAndSource both;
  EarthState earth;
  BarycenterInput baryinput;
  EmissionTime emit;

  
  REAL8 RA = 0;
  REAL8 DEC = 0; 
  REAL8 psi = 0;
  REAL8 F_check = 0;
  REAL8 fdot = 0;
  
  double dt = 1.0/2.0;                               
  long int n_steps = (86400*10-500)*2.0;             
 
  double inpdt = 1.0/2.0;              
  long int inpsteps = 2000000;
  int fileopen = 0;
  ifstream file;
  ifstream freqfile;
  ofstream ofile;
  
  for(int i=1;i<argc;i++)             
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

  //cerr<<F_check<<" "<<RA<<" "<<DEC<<" "<<inpsteps<<" "<<dt<<" "<<n_steps<<endl;

  status.statusPtr=NULL;                                    
  
  detector = lalCachedDetectors[LALDetectorIndexLHODIFF]; 
  edata.ephiles.earthEphemeris = EARTHDATA;               
  edata.ephiles.sunEphemeris = SUNDATA;                   
  edata.leap = 13;                                        
  LALInitBarycenter( &status, &edata);          
 
  if(status.statusCode) 
    {
      printf("Unexpectedly got error code %d and message %s\n",
	     status.statusCode, status.statusDescription);
      return 0;
    }
	
  

  ofile.open("output");
  int doonce = 1;
  
  if(!fileopen)
    file.open("T");
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
      //cout<<i<<" "<<x_re[i]<<" "<<x_im[i]<<endl;
    }


  string si = "data", so = "output",s;
  for(int i=0;i<1;i++)
    for(int j=0;j<1;j++)
      for(int k=0;k<1;k++)
	{
	  //fprintf(stderr,"%d\n",k);
	  freq F;
	  F.num = 1;
	  F.fk = (double*)malloc(sizeof(double));
	  F.fk[0] = fdot;
	  
	  in1 = (fftw_complex *)malloc(sizeof(fftw_complex)*n_steps);
	  in2 = (fftw_complex *)malloc(sizeof(fftw_complex)*n_steps);
	  
	  y = (fftw_complex *)malloc(sizeof(fftw_complex)*n_steps);
	  //char stemp[50];
	  //sprintf(stemp,"%d",k);
	  //s = so + stemp;
	  //cerr<<s<<" "<<inpsteps<<" "<<n_steps<<endl;
	  
	  //RA = .1*(i);//4.88671;                           
	  //DEC = .01*(j);//0.217584; 
	  LIGOTimeGPS t0;                                        
	  t0.gpsSeconds = 700010000;                            
	  t0.gpsNanoSeconds = 0;
	  source.equatorialCoords.latitude = DEC;                
	  source.equatorialCoords.longitude = RA;                
	  source.equatorialCoords.system=COORDINATESYSTEM_EQUATORIAL;
	  source.orientation = psi;
	    
	  both.pSource = &source;                    
	  both.pDetector = &detector;
	  
	  baryinput.site.location[0] = detector.location[0]/LAL_C_SI; 
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
	  //cerr<<" here \n";
	  getTbary(tdout,Phi_m,t0,n_steps,dt,status,baryinput,edata);
	  //cerr<<" got here too \n";
	  //for(int i=0;i<n_steps;i++)
	  // tdout[i] = gettime(t+dt*double(i));
	  
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
	  interp_spline(t0,dt,inpsteps,n_steps,tdout,x_re,x_im,y);

	  
	 
	  

	  double y1_re,y1_im,y2_re,y2_im;
	  double phi,sinphi,cosphi;
	  double Het = 106.0 + (double)k/1.0;	  

	  for(int i=0;i<n_steps;i++)
	    {
	      double retemp,imtemp,shift = Het ,DT = -Phi_m[i];
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
	    Scalc[i] = 2*dt;//3.22e-44;	    
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
	    }
	  //ofileZ<<endl;
	  //cerr<<F_check<<endl;
	  for(double i=0;i<n_steps/2.0;i++)
	    //cerr<<(double)i/(double)n_steps/dt<<" "<<(double)(i+1)/(double)n_steps/dt<<endl;
	    if(fabs(F_check-((double)i/(double)(n_steps)/dt+Het))< 1.0/n_steps/dt/2.0)
	      printf("%g %g %6.10g %g %g\n",RA,F_check,(double)i/(double)(n_steps)/dt+Het,Fstat[(int)(i)],Fstat[(int)(i)]);
	    
	      for(double i=n_steps-1;i>n_steps/2.0;i--)
		if(fabs(F_check-((double)(i-n_steps)/(double)(n_steps)/dt+Het)) < 1.0/n_steps/dt/2.0) 
		printf("%g %g %6.10g %g %g\n",RA,F_check,(double)(i-n_steps)/(double)(n_steps)/dt+Het,Fstat[(int)(i)],Fstat[(int)(i)]);
	  
	  
	  /*-------------------------------------------------*/

	  double MAX;
	  int index;
	  //MAX = max2(Fstat,n_steps,index);
	  //if(index<n_steps/2)
	  //  printf("%6.10g %g %6.10g \n",Het,MAX,(double)index/(double)(n_steps)/dt+105.0);
	  // else
	  //  printf("%6.10g %g %6.10g \n",Het,MAX,(double)(index-n_steps)/(double)(n_steps)/dt+105.0);








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
