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

#include <lal/SFTfileIO.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <unistd.h>
#include<fftw3.h>
#include<zlib.h>
#include<gsl/gsl_blas.h>
//#include<sys/types.h>
//#include<sys/stat.h>
//#include<glob.h>



/*---------- empty initializers ---------- */
LALStatus empty_status;
SFTConstraints empty_constraints;
/*---------- Global variables ----------*/

INT4 lalDebugLevel = 1;
SFTVector *sft_vect = NULL;
REAL8 *sinVal,*cosVal;

typedef enum 
  {
    BARTLETT = 1, 
    HANN,
    RECTANGULAR,
    WELCH
  }WindowType;
    

double GPS2REAL(LIGOTimeGPS time)
{
  return(time.gpsSeconds+1e-9*time.gpsNanoSeconds);
}

typedef struct 
{
  int length;
  int cont[10000];
  double gap[10000];
}conti;

void window(WindowType W,REAL8* left,REAL8* right,int N)
{
  int i;
  REAL8 y,z;
  for(i=0;i<N;i++)
    {
      y = -1.0 + (REAL8)i/(REAL8)(N-1);
      z = (REAL8)i/(REAL8)(N-1);
      switch(W)
	{
	case BARTLETT:
	  left[i] = 1.0 + y;
	  right[i] = 1.0 - z;
	  break;
	case HANN:
	  left[i] = cos(M_PI*y/2.0)*cos(M_PI*y/2.0);
	  right[i] = cos(M_PI*z/2.0)*cos(M_PI*z/2.0);
	  break;
	case WELCH:
	  left[i] = 1.0 - y*y;
	  right[i] = 1.0 - z*z;
	  break;
	default: 
	  fprintf(stderr,"Sorry, Wrong Window choice\n");
	}
    }
}

void ApplyWindow(fftw_complex* X,int N,int NW,WindowType W)
{
  REAL8 *left,*right;
  left = (REAL8*)malloc(sizeof(REAL8)*NW);
  right = (REAL8*)malloc(sizeof(REAL8)*NW);
  int i;
  window(W,left,right,NW);

  for(i=0;i<NW;i++)
    X[i][0] *= left[i];
  for(i=N-NW;i<N;i++)
    X[i][0] *= right[NW-N+i];
  for(i=0;i<NW;i++)
    X[i][1] *= left[i];
  for(i=N-NW;i<N;i++)
    X[i][1] *= right[NW-N+i];
  free(left);
  free(right);
}
    
int CSFTs(fftw_complex *L,REAL8 Fmin,REAL8 Fmax,int number,int startindex,REAL8 sTsft,LIGOTimeGPS ref, int Dterms)
{ 
  REAL8 t0 = GPS2REAL(ref);
  REAL8 lTsft = number*sTsft;
  INT4 alpha,m;                 /* loop indices */
  REAL8	xTemp;	                /* temp variable for phase model */
  REAL8	deltaF;	                /* width of SFT band */
  INT4	k, k1;	                /* defining the sum over which is calculated */
  REAL8	x;		        /* local variable for holding x */
  REAL8	realXP, imagXP; 	/* temp variables used in computation of */
  REAL8	realP, imagP;	        /* real and imaginary parts of P, see CVS */
  INT4	nDeltaF;	        /* number of frequency bins per SFT band */
  INT4	sftIndex;	        /* more temp variables */
  REAL8	y;		        /* local variable for holding y */
  REAL8 realQ, imagQ;

  INT4 index;

  COMPLEX16 llSFT;
  COMPLEX8 lSFT;

  REAL8 f;
  int errorcode;

  REAL8 if0 = floor(Fmin*sTsft);
  REAL8 if1 = ceil((Fmax-Fmin+Fmin)*sTsft);

  REAL8 ifmin = floor(Fmin*sTsft)-Dterms;
  REAL8 ifmax = ceil((Fmax-Fmin+Fmin)*sTsft)+Dterms;//ceil(Fmax*sTsft)+Dterms;
  
  

  /* variable redefinitions for code readability */
  deltaF=sft_vect->data->deltaF;
  nDeltaF=sft_vect->data->data->length;

  

  // Loop over frequencies to be demodulated
  for(m = 0 ; m <= number*(if1-if0)  ; m++ )
  {
    llSFT.re =0.0;
    llSFT.im =0.0;

    f=if0*deltaF+m*deltaF/number;

    // Loop over SFTs that contribute to F-stat for a given frequency
    for(alpha=0;alpha<number;alpha++)
      {
	//fprintf(stderr,"start %d\n",m);
	REAL8 tsin, tcos, tempFreq;
	COMPLEX8 *Xalpha = sft_vect->data[alpha+startindex].data->data;
	xTemp = (REAL8)if0+(REAL8)m/(REAL8)number;
	realXP = 0.0;
	imagXP = 0.0;
	      	/* find correct index into LUT -- pick closest point */
	tempFreq = xTemp-(INT4)xTemp;
	index=(INT4)(tempFreq*64+0.5); //just like res above
	      
	{
	  REAL8 d=LAL_TWOPI*(tempFreq-(REAL8)index/64.0);//just like res above
	  REAL8 d2=0.5*d*d;
	  REAL8 ts=sinVal[index];
	  REAL8 tc=cosVal[index];
		
	  tsin=ts+d*tc-d2*ts;
	  tcos=tc-d*ts-d2*tc-1.0;
	}

        tempFreq=LAL_TWOPI*(tempFreq+Dterms-1);
        k1=(INT4)xTemp-Dterms+1;
        /* Loop over terms in dirichlet Kernel */
        for(k=0;k<2*Dterms;k++)
	  {
	    COMPLEX8 Xalpha_k;
	    x = tempFreq-LAL_TWOPI*(REAL8)k;
	    realP = tsin/x;
	    imagP = tcos/x;

	    /* If x is small we need correct x->0 limit of Dirichlet kernel */
	    if(fabs(x) < 0.000001) 
	      {
		realP = 1.0;
		imagP = 0.0;
	      }	 
 
	    sftIndex=k1+k-ifmin+1;

	   
	    /* these four lines compute P*xtilde */
	    Xalpha_k = Xalpha[sftIndex];
	    //fprintf(stderr,"%d\n",sftIndex);
	    realXP += Xalpha_k.re*realP;
	    realXP -= Xalpha_k.im*imagP;
	    imagXP += Xalpha_k.re*imagP;
	    imagXP += Xalpha_k.im*realP;
	  }
	//double time = GPS2REAL(sft_vect->data[alpha+startindex].epoch) - t0;
	//y = -LAL_TWOPI*(time/sTsft)*(if0+(REAL8)m/(REAL8)number);
	//printf("Alpha = %g \n",time/sTsft);
	y = -LAL_TWOPI*alpha*(if0+(REAL8)m/(REAL8)number);
	//fprintf(stderr,"Time %g , time/sTsft %g\n",time,time/sTsft);

	realQ = cos(y);
	imagQ = sin(y);

	/* implementation of amplitude demodulation */
	{
	  REAL8 realQXP = realXP*realQ-imagXP*imagQ;
	  REAL8 imagQXP = realXP*imagQ+imagXP*realQ;
	  llSFT.re += realQXP;
	  llSFT.im += imagQXP;
	}
      }      


    L[m][0] = llSFT.re; 
    L[m][1] = llSFT.im; 
    
    
  }

  return 0;

}

void init()
{
  int k = 0;
  int res=64;
  sinVal=(REAL8 *)malloc((res+1)*sizeof(REAL8));
  cosVal=(REAL8 *)malloc((res+1)*sizeof(REAL8)); 
  for (k=0; k<=res; k++)
    {
      sinVal[k]=sin((LAL_TWOPI*k)/res);
      cosVal[k]=cos((LAL_TWOPI*k)/res);
    }
}

void reshuffle(fftw_complex *X, int N)
{
  fftw_complex *Temp;
  Temp = (fftw_complex*) malloc(sizeof(fftw_complex)*N);
  int i,k = 0;
  for(i=N/2.0;i<N;i++)
    {
      Temp[k][0] = X[i][0];
      Temp[k++][1] = X[i][1];
    }
  for(i=0;i<N/2.0;i++)
    {
      Temp[k][0] = X[i][0];
      Temp[k++][1] = X[i][1];
    }
  for(i=0;i<N;i++)
    {
      X[i][0] = Temp[i][0];
      X[i][1] = Temp[i][1];
    }
  fftw_free(Temp);
}
     

int main(int argc, char **argv)
{
  int Dterms = 32;
  INT4 i,j,k=0;
  FILE *fp,*fp2;
  REAL8 InputST = 500000.0;
  LALStatus status = empty_status;

  WindowType W = HANN;
  SFTCatalog *catalog = NULL;
  UINT4 numBins, nSFT;
  double f_min,f_max; 
  SFTConstraints constraints=empty_constraints;
  LIGOTimeGPS startTime, endTime; 
  double avg;
  char outfile[128],outfile2[128];
  if (argc!=8)
  {
   fprintf(stderr, "startGPS endGPS Detector F_min f_max where InputST?\n"); 
   return(0);
  }  
  
  startTime.gpsSeconds = atoi(argv[1]);
  startTime.gpsNanoSeconds = 0;
  constraints.startTime = &startTime; 
  
  endTime.gpsSeconds = atoi(argv[2]);
  endTime.gpsNanoSeconds = 0;
  constraints.endTime = &endTime;
  constraints.detector = argv[3];
  
  InputST = atof(argv[7]);
  
  f_min = atof(argv[4]);
  f_max = atof(argv[5]);
  
  f_min -= (REAL8)Dterms/InputST;
  f_max += (REAL8)Dterms/InputST;

  fprintf(stderr,"About to Find Files\n");
  
  LALSFTdataFind ( &status, &catalog, argv[6], &constraints );
 
  
  if(status.statusCode) 
    {
      printf("Unexpectedly got error code %d and message %s\n",
             status.statusCode, status.statusDescription);
      return 0;
    }             

  fprintf(stderr,"Catalog made %d\n",catalog->length);

  LALLoadSFTs ( &status, &sft_vect, catalog, f_min, f_max);
  
  if(status.statusCode) 
    {
      printf("Unexpectedly got error code %d and message %s\n",
             status.statusCode, status.statusDescription);
      return 0;
    }             

  fprintf(stderr,"SFTs loaded in\n");

  f_min += (REAL8)Dterms/InputST;
  f_max -= (REAL8)Dterms/InputST;

  double contstarttime;
  double diff;
  int stillcont = 1;
  int first = 1;
  int numcont = 1;
  conti C;
  
  C.length = 0;

  init();
  
  fprintf(stderr,"about to start loop fmax = %g fmin = %g \n",f_max,f_min);
  for(i=0;i<sft_vect->length;i++)
    {
      fprintf(stderr,"GPS Time %lf\n",GPS2REAL(sft_vect->data[i].epoch));
      if(first)
	{
	  contstarttime = GPS2REAL(sft_vect->data[i].epoch);
	  first = 0;
	  numcont = 1;
	}
      else
	{
	   diff = GPS2REAL(sft_vect->data[i].epoch)-GPS2REAL(sft_vect->data[i-1].epoch);
	  //fprintf(stderr,"Diff = %lf\n Cont = %d\n",diff,numcont);
	  if(diff == InputST)
	    numcont++;
	      
	  else
	    {
	      stillcont = 0;
	      first = 1;
	      i--;
	      C.gap[k] = diff;
	      C.cont[k++] = numcont;
	      
	      //fprintf(stderr,"Diff = %lf\n Cont = %d\n",diff,numcont);
	    }
	}
    }
  C.gap[k] = 0;
  C.cont[k++] = numcont;
  C.length = k;

  REAL8 sTsft = 1.0/sft_vect->data->deltaF;
  REAL8 deltaF = sft_vect->data->deltaF;
  fftw_complex *L;
  fftw_complex *TS,*TempT;
  REAL8 tdiff = (GPS2REAL(endTime)-GPS2REAL(startTime))*(f_max-f_min)+1;
  int Tlength = ceil(tdiff);
  

  REAL8 START = GPS2REAL(startTime);
  double dt = 1/(f_max-f_min);
  
  TS = (fftw_complex*)malloc(sizeof(fftw_complex)*Tlength);
 
  for(i=0;i<Tlength;i++)
    {
      TS[i][0] = 0;
      TS[i][1] = 0;
    }

  fprintf(stderr," %lf %d \n",sTsft,Tlength);
  
  
	  
  //fprintf(stderr," Name %s GPSSeconds %lf ",sft_vect->data[i].name,sft_vect->data[i].epoch.gpsSeconds+1e-9*sft_vect->data[i].epoch.gpsNanoSeconds);
  int startindex = 0;
  int sftindex = 0;
  int Tlapsed = 0;
  int err = 0;

  for(i=0;i<C.length;i++)
    {
      fprintf(stderr,"\n %d %g\n",C.cont[i],C.gap[i]);

      int N = (sTsft*C.cont[i]*(f_max-f_min))+1;
      //if(fmod(N,2.0)==0)
      //N++;
      L = (fftw_complex*)malloc(sizeof(fftw_complex)*N);
      TempT = (fftw_complex*)malloc(sizeof(fftw_complex)*N);
      

      //fprintf(stderr,"\n %lf %lf \n",(int)sTsft*C.cont[i]*(f_max-f_min)+1-N,N);

      if(C.cont[i] >  1)
	 err = CSFTs(L,f_min,f_max,C.cont[i],startindex,sTsft,startTime,Dterms);
      else
	{
	  for(j=0;j<N;j++)
	    {
	      //printf("j = %d , Other = %d \n",j,j+Dterms+1);
	      L[j][0] = sft_vect->data[i].data->data[j+Dterms].re;
	      L[j][1] = sft_vect->data[i].data->data[j+Dterms].im;
	    }
	  
	}
      
      for(j=0;j<N;j++)
	{
	  double F = (double)j/(double)(N-1)*(f_max-f_min)+f_min;
	  //printf(" %g %g %g %g \n",F,L[j][0],L[j][1],sqrt(L[j][0]*L[j][0]+L[j][1]*L[j][1]));
	}
      
      fftw_plan p = fftw_plan_dft_1d(N,L,TempT,FFTW_BACKWARD,FFTW_ESTIMATE);
      ApplyWindow(L,N,Dterms,W);
      
      reshuffle(L,N);
      fftw_execute(p);
      startindex+=C.cont[i];
      
      fprintf(stderr,"%d\n",N);

      ApplyWindow(TempT,N,10,W);
      //if(i==0)
      for(j=0;j<N;j++)
	{
	  TS[sftindex][0] = TempT[j][0]*deltaF/(REAL8)C.cont[i];
	  TS[sftindex++][1] = TempT[j][1]*deltaF/(REAL8)C.cont[i];
	  double F = (double)j/(double)(N-1)*(f_max-f_min)+f_min;
	  //printf(" %g %g %g %g \n",F,L[j][0],L[j][1],sqrt(L[j][0]*L[j][0]+L[j][1]*L[j][1]));
	}
      
      if(C.gap[i] != 0)
	{
	  sftindex += (C.gap[i]-InputST)/dt;
	  fprintf(stderr,"SFT Index Jump = %g",(C.gap[i]));
	}
      
//printf(" Real %g , Imag %g \n",sft_vect->data[0].data->data->re,sft_vect->data[0].data->data->im);
      fftw_free(L);
      fftw_free(TempT);
      fftw_destroy_plan(p);
    }

  LALDestroySFTCatalog( &status, &catalog);
  
  numBins = sft_vect->data->data->length;
  nSFT = sft_vect->length;
  

  //heterodyne(TS,Tlength,-f_min,1.0/(f_max-f_min));
  for(i=0;i<Tlength;i++)
   printf("%10.15g %g %g\n",START+dt*i,TS[i][0],TS[i][1]);    
    //printf("%.18g %g %g\n",(REAL8)i/(f_max-f_min)/2.0+START,TS[i][0],TS[i][1]);
  
  fprintf(stderr, "nSFT = %d\t numBins = %d\t df = %f\n", nSFT, numBins,sft_vect->data->deltaF);
  
  free(TS);

  
  
 LALDestroySFTVector (&status, &sft_vect );

 return(0);

}
