/**
 * \file
 * \ingroup pulsarApps
 * \brief
 * This code takes in SFTs as input and creates a heterodyned, downsampled time series out of them.
 * It will fill in the gaps in SFTs with zeros and combine consecutive SFTs are combined coherently.
 */



#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/SFTfileIO.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include<fftw3.h>
#include<zlib.h>
#include<gsl/gsl_blas.h>



/*---------- empty initializers ---------- */
LALStatus empty_status;
SFTConstraints empty_constraints;

/*---------- Global variables ----------*/

INT4 lalDebugLevel = 1;
SFTVector *sft_vect = NULL;
REAL8 *sinVal,*cosVal;

/*------ An alias for Windowtype --------*/
typedef enum 
  {
    BARTLETT = 1, 
    HANN,
    RECTANGULAR,
    WELCH
  }WindowType;
 
/*------ Converts LIGO GPS times to real times ------*/   

double GPS2REAL(LIGOTimeGPS time)
{
  return(time.gpsSeconds+1e-9*time.gpsNanoSeconds);
}

/*------ A book-keeping structure which stores information about the continuity of the SFTs input along with the gaps in between each contiguous block ------*/
          
typedef struct 
{
  /*------- Number of contiguous blocks -------*/
  int length;
  /*------- Number of SFTs in each continous blocks -------*/
  int cont[10000];
  /*------- Gap between each block in seconds -------*/
  double gap[10000];
}conti;

/*------ Window gives you the values of the leading and falling edge of each window type. They rise from 0 to 1 and fall from 1 to 0 ------*/
/******* WindowType W -> Type of Window ********/
/******* left and right -> arrays corresponding to the leading and falling edge *******/
/******* N -> Number of terms in left and right *******/
void window(WindowType W,REAL8* left,REAL8* right,int N)
{
  int i;
  /****** y -> left and z -> right ******/
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

/****** ApplyWindow applies a window to a complex array X, N -> number of terms in X and NW -> Number of terms to keep in the window, W -> Window Type ******/

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


/****** CSFTs combines SFTs coherently to give a longer time baseline SFT, written by Xavier Siemens, modified by Pinkesh Patel ******/
    
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
  REAL8 ifmax = ceil((Fmax-Fmin+Fmin)*sTsft)+Dterms;/*ceil(Fmax*sTsft)+Dterms;*/
  
  

  /* variable redefinitions for code readability */
  deltaF=sft_vect->data->deltaF;
  nDeltaF=sft_vect->data->data->length;

  

  /* Loop over frequencies to be demodulated */
  for(m = 0 ; m <= number*(if1-if0)  ; m++ )
  {
    llSFT.real_FIXME =0.0;
    llSFT.im =0.0;

    f=if0*deltaF+m*deltaF/number;

    /* Loop over SFTs that contribute to F-stat for a given frequency */
    for(alpha=0;alpha<number;alpha++)
      {
	/* fprintf(stderr,"start %d\n",m); */
	REAL8 tsin, tcos, tempFreq;
	COMPLEX8 *Xalpha = sft_vect->data[alpha+startindex].data->data;
	xTemp = (REAL8)if0+(REAL8)m/(REAL8)number;
	realXP = 0.0;
	imagXP = 0.0;
	      	/* find correct index into LUT -- pick closest point */
	tempFreq = xTemp-(INT4)xTemp;
	index=(INT4)(tempFreq*64+0.5); /*just like res above */
	      
	{
	  REAL8 d=LAL_TWOPI*(tempFreq-(REAL8)index/64.0);/*just like res above */
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
	    /*fprintf(stderr,"%d\n",sftIndex); */
	    realXP += crealf(Xalpha_k)*realP;
	    realXP -= cimagf(Xalpha_k)*imagP;
	    imagXP += crealf(Xalpha_k)*imagP;
	    imagXP += cimagf(Xalpha_k)*realP;
	  }
	/* double time = GPS2REAL(sft_vect->data[alpha+startindex].epoch) - t0; */
	/* y = -LAL_TWOPI*(time/sTsft)*(if0+(REAL8)m/(REAL8)number); */
	/*printf("Alpha = %g \n",time/sTsft);*/
	y = -LAL_TWOPI*alpha*(if0+(REAL8)m/(REAL8)number);
	/*fprintf(stderr,"Time %g , time/sTsft %g\n",time,time/sTsft);*/

	realQ = cos(y);
	imagQ = sin(y);

	/* implementation of amplitude demodulation */
	{
	  REAL8 realQXP = realXP*realQ-imagXP*imagQ;
	  REAL8 imagQXP = realXP*imagQ+imagXP*realQ;
	  llSFT.real_FIXME += realQXP;
	  llSFT.im += imagQXP;
	}
      }      


    L[m][0] = creal(llSFT); 
    L[m][1] = llSFT.im; 
    
    
  }

  return 0;

}

/****** Look-up table for sin and cos ******/
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

/****** FFTW3 requires a format such that an array has the positive frequencies first followed by a ascending negative frequencies. This function inverts an array which has negative and positive frequencies lined up in the logical manner to this fftw3 manner ******/
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
  /****** Number of Direchlet Terms used in combining the SFTs ******/
  int Dterms = 32;

  /****** Book-keeping + loop variables ******/
  INT4 i,j,k=0;
  
  /*FILE *fp,*fp2;*/

  /****** The time baseline of input SFTs default to 1800.0 seconds ******/
  REAL8 InputST = 1800.0;

  /****** LALStatus for debugging ******/
  LALStatus status = empty_status;

  /****** Type of Window used (its not used for the whole data, but only at the edges, so its split in two parts and attached to ramp function) ******/
  WindowType W = HANN;

  /****** Catalog of SFTs found at the given location ******/
  SFTCatalog *catalog = NULL;
  
  /****** Number of bins in each SFT and the number of total SFTs ******/
  UINT4 numBins, nSFT;
 
  /****** Minimum and Maximum frequencies of the band of interest ******/
  double f_min,f_max; 
  
  /****** Internal LAL variable ******/
  SFTConstraints constraints=empty_constraints;
  
  /****** Start and End times of the time series to be outputted ******/ 
  LIGOTimeGPS startTime, endTime; 

  /****** Hard Coded (BAD!!) argument variables ******/
  if (argc!=8)
  {
   fprintf(stderr, "startGPS endGPS Detector F_min f_max where InputST?\n"); 
   return(0);
  }  
  
  /****** Assigning arguments and creating SFT constraints ******/
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
  
  /****** Padding up the f_min and f_max to take care of Dterms ******/
  f_min -= (REAL8)Dterms/InputST;
  f_max += (REAL8)Dterms/InputST;

  fprintf(stderr,"About to Find Files\n");
  
  /****** Finding SFTs ******/
  LALSFTdataFind ( &status, &catalog, argv[6], &constraints );
 
  
  if(status.statusCode) 
    {
      printf("Unexpectedly got error code %d and message %s\n",
             status.statusCode, status.statusDescription);
      return 0;
    }             

  fprintf(stderr,"Catalog made %d\n",catalog->length);

  /****** Loading SFTs into sft_vect ******/
  LALLoadSFTs ( &status, &sft_vect, catalog, f_min, f_max);
  
  if(status.statusCode) 
    {
      printf("Unexpectedly got error code %d and message %s\n",
             status.statusCode, status.statusDescription);
      return 0;
    }             

  fprintf(stderr,"SFTs loaded in\n");

  /****** Readjusting f_min and f_max after loading in the SFTs ******/
  f_min += (REAL8)Dterms/InputST;
  f_max -= (REAL8)Dterms/InputST;

  /****** Book-keeping variables to calculate the contiguity structure for the loaded SFTs ******/
  double contstarttime;
  double diff;
  int first = 1;
  int numcont = 1;
  
  /****** The contiguity structure ******/
  conti C;
  
  C.length = 0;

  init();
  
  fprintf(stderr,"about to start loop fmax = %g fmin = %g \n",f_max,f_min);
  for(i=0;i<sft_vect->length;i++)
    {
      fprintf(stderr,"GPS Time %lf\n",GPS2REAL(sft_vect->data[i].epoch));
      /****** If the SFT in question is the first one of a contiguous batch ******/
      if(first)
	{
	  /****** Record the starttime of this continuous bunch ******/
	  contstarttime = GPS2REAL(sft_vect->data[i].epoch);
	  /****** No longer first! ******/
	  first = 0;
	  /****** Number of SFTs in this continuous bunch is now 1 ******/
	  numcont = 1;
	}
      else
	{
	  /****** Calculate the difference between the times of this SFT from the one previous ******/
	  diff = GPS2REAL(sft_vect->data[i].epoch)-GPS2REAL(sft_vect->data[i-1].epoch);

	  /****** If the difference in time = time base line of SFT, then they are contiguous, therefore increase count by 1 ******/
	  if(diff == InputST)
	    numcont++;
	  
	  /****** If not, then this is a new block and so we restart with this SFT as being the first SFT of the next block ******/
	  else
	    {
	      /****** New first SFT ******/
	      first = 1;

	      /****** We must restart the cycle and redo this SFT ******/
	      i--;

	      /****** Record the gap between these two blocks ******/
	      C.gap[k] = diff;

	      /****** Record the number of SFTs contiguous in this block ******/
	      C.cont[k++] = numcont;
	    }
	}
    }

  /****** Record the information for the last SFT ******/
  C.gap[k] = 0;
  C.cont[k++] = numcont;
  C.length = k;

  /****** deltaF and 1/deltaF recorded directly from the SFTs ******/
  REAL8 sTsft = 1.0/sft_vect->data->deltaF;
  REAL8 deltaF = sft_vect->data->deltaF;

  /****** L is the combined SFT for each contiguous block******/
  fftw_complex *L;
  
  /****** TS is the time series which will be output and TempT is the time series generated by each block which is then patched up along with the gaps to create TS ******/
  fftw_complex *TS,*TempT;
  
  /****** tdiff and Tlength are the number of points in TS ******/
  REAL8 tdiff = (GPS2REAL(endTime)-GPS2REAL(startTime))*(f_max-f_min)+1;
  int Tlength = ceil(tdiff);
  
  /****** Start time of the time series ******/
  REAL8 START = GPS2REAL(startTime);

  /****** Since the time series in inherectly complex, dt = 1/df and not dt = 1/(2.0 * df) [ for real time series, this would be the case ] ******/
  double dt = 1/(f_max-f_min);
  
  /****** Assign Memory ******/
  TS = (fftw_complex*)malloc(sizeof(fftw_complex)*Tlength);
 
  /****** Set initially to zeros, this takes care of all gaps and all one has to do now, is to calculate the time series for each continuous chunk, window it and put it in the appropraite place in TS ******/
  for(i=0;i<Tlength;i++)
    {
      TS[i][0] = 0;
      TS[i][1] = 0;
    }

  fprintf(stderr," %lf %d \n",sTsft,Tlength);
  
  
	  
  /*fprintf(stderr," Name %s GPSSeconds %lf ",sft_vect->data[i].name,sft_vect->data[i].epoch.gpsSeconds+1e-9*sft_vect->data[i].epoch.gpsNanoSeconds);*/

  /****** Book-keeping variables for the loop over contiguous blocks ******/
  int startindex = 0;
  int sftindex = 0;
  int Tlapsed = 0;
  int err = 0;

  /****** Loop over contiguous blocks in data ******/
  for(i=0;i<C.length;i++)
    {
      fprintf(stderr,"\n %d %g\n",C.cont[i],C.gap[i]);

      /****** Number of data points in this continuous block ******/
      int N = (sTsft*C.cont[i]*(f_max-f_min))+1;

      /****** Assign temporary memory to the SFT and Temporary time series ******/
      L = (fftw_complex*)malloc(sizeof(fftw_complex)*N);
      TempT = (fftw_complex*)malloc(sizeof(fftw_complex)*N);
      
      /****** Combine SFTs only if the number of SFTs > 1 ******/
      if(C.cont[i] >  1)
	 err = CSFTs(L,f_min,f_max,C.cont[i],startindex,sTsft,startTime,Dterms);

      /****** Else just assign the values to L ******/
      else
	{
	  for(j=0;j<N;j++)
	    {
	      L[j][0] = crealf(sft_vect->data[i].data->data[j+Dterms]);
	      L[j][1] = cimagf(sft_vect->data[i].data->data[j+Dterms]);
	    }
	  
	}
      
      /*for(j=0;j<N;j++)*/
      /*	{*/
      /* double F = (double)j/(double)(N-1)*(f_max-f_min)+f_min;*/
      /*printf(" %g %g %g %g \n",F,L[j][0],L[j][1],sqrt(L[j][0]*L[j][0]+L[j][1]*L[j][1]));*/
      /* }*/

      /****** Create an FFTW3 Plan (Backward, since its an inverse fourier transform) ******/
      fftw_plan p = fftw_plan_dft_1d(N,L,TempT,FFTW_BACKWARD,FFTW_ESTIMATE);

      /****** Apply a window to the SFT , before inverse transforming ******/
      ApplyWindow(L,N,Dterms,W);
      
      /****** Reshuffle to get into the FFTW3 format ******/
      reshuffle(L,N);

      /****** Run the IFFT ******/
      fftw_execute(p);

      /****** Update the startindex, which keeps track of how many SFTs have already been processed ******/
      startindex+=C.cont[i];
      
      fprintf(stderr,"%d\n",N);

      /****** Now apply a window to the data in the time domain, so that it fits into the gaps smoothly ******/
      ApplyWindow(TempT,N,10,W);

      /****** Fit in the TempT in the right place in TS, the right place is kept in track by the variable sftindex ******/
      for(j=0;j<N;j++)
	{
	  TS[sftindex][0] = TempT[j][0]*deltaF/(REAL8)C.cont[i];
	  TS[sftindex++][1] = TempT[j][1]*deltaF/(REAL8)C.cont[i];
	  /*double F = (double)j/(double)(N-1)*(f_max-f_min)+f_min;*/
	  /*printf(" %g %g %g %g \n",F,L[j][0],L[j][1],sqrt(L[j][0]*L[j][0]+L[j][1]*L[j][1]));*/
	}
      
      /****** Update the sftindex to keep track of where in TS we are inserting the next TempT ******/
      if(C.gap[i] != 0)
	{
	  sftindex += (C.gap[i]-InputST)/dt;
	  fprintf(stderr,"SFT Index Jump = %g",(C.gap[i]));
	}
      
      /*printf(" Real %g , Imag %g \n",sft_vect->data[0].data->data->re,sft_vect->data[0].data->data->im);*/
      fftw_free(L);
      fftw_free(TempT);
      fftw_destroy_plan(p);
    }

  LALDestroySFTCatalog( &status, &catalog);
  
  numBins = sft_vect->data->data->length;
  nSFT = sft_vect->length;
  

  for(i=0;i<Tlength;i++)
   printf("%10.15g %g %g\n",START+dt*i,TS[i][0],TS[i][1]);    
  
  fprintf(stderr, "nSFT = %d\t numBins = %d\t df = %f\n", nSFT, numBins,sft_vect->data->deltaF);
  
  free(TS);

  
  
 LALDestroySFTVector (&status, &sft_vect );

 return(0);

}
