/*  <lalVerbatim file="GEO600FrameChannelsCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{GEO600FrameChannels.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
%% \input{GEO600FrameChannelsCP}
%% \idx{GEO600FrameChannels()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GEO600FrameChannelsCV}}

</lalLaTeX>  */

/*
File to contain a routine to get a single channel of Data continuously 
from a directory of frame files 
*/
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define PI 3.141592654
#include <FrameL.h>
#include "GEO600FrameChannels.h"

typedef struct tagComplex
{
    float r;
    float i;
}Complex;

/* Static Function prototypes */
static float mygasdev(long *);
static float ran2(long *idum);
static int SetFileNames(int);
static Complex mult(Complex *, Complex *);
static Complex divi(Complex *, Complex *);
static int GetCalibrationFunc(int, const float *, int, float *, float);

/* Static internal Variables */

static char *baseDirectory=NULL;
static long idum = -784387882;
static char currFile[512];
static const int numFramesInFile = 60;
static const int gpsFid=662342400;
static const int hourFid=0;
static const int dayFid=1;
static const int yearFid=2001;

int
SetBaseDirectory(char *baseDir)
{
    int len;
    struct stat buf;
    int ret;

    ret = stat(baseDir,&buf);
    if(ret<0){
        fprintf(stderr, "base Directory does not exist\n%s\n",baseDir);
        fflush(stdout);
        return -1;
    }
    len = strlen(baseDir);
    baseDirectory = malloc(len+1);
    strcpy(baseDirectory,baseDir);
    return 0;
}

static int
SetFileNames(int gpsTime)
{
    int gpsMin,ret;
    struct stat buf;
    int year=yearFid, day=dayFid, hour=hourFid, gps=gpsFid;

    if(baseDirectory == NULL){
        fprintf(stderr, "BaseDirectory Not Set\n");
        fflush(stdout);
        return -1;
    }
    if(gpsTime < (gpsFid + 60)){
        fprintf(stderr, "Data not available prior to Jan 1st 2001\n");
        fprintf(stderr, "Incorrect GPSTime\n");
        fflush(stdout);
        return(-1);
    }
    gpsMin = gpsTime - gpsTime%60;
    while(gps!=gpsMin){
        gps+=60;
        if(gps%3600 == 0){
            hour++;
            if(hour==24){
                hour=0;
                day++;
                if(day==366){
                    day = 1;
                    year++;
                }
            }
            
        }
    }
    sprintf(currFile,"%s/%d/day%03d/hour%02d/GEO_%d.fr",baseDirectory,year,day,hour,gps);
    ret = stat(currFile,&buf);
    if(ret<0){
        fprintf(stderr, "Not Able to Open File %s \n",currFile);
        fflush(stdout);
        return -1;
    }
    return 0;
}



int
GetSampleRate(int gpsTime, char *chan, float *sampleRate)
{
    FrFile *frFile;
    FrameH *frameH;
    FrAdcData *adcStr;
    
    if(SetFileNames(gpsTime)<0){
        fprintf(stderr, "Recursive Error:GetSampleRate\n");
        fflush(stdout);
        return -1;
    }
    if((frFile=FrFileINew(currFile))==NULL){
        fprintf(stderr, "Cannot Open Frame File %s\n",currFile);
        fflush(stdout);
        return -1;
    }
    if((frameH = FrameReadTAdc(frFile,(double)gpsTime,chan))==NULL){
        fprintf(stderr, "Cannot Read Frame for channel %s for gps\n",chan);
        fprintf(stderr, " time %d in file \n",gpsTime);
        fflush(stdout);
        return -1;
    }
    if((adcStr = FrAdcDataFind(frameH, chan))==NULL){
        fprintf(stderr, "Cant get pointer to adcStruct in Frame for channel %s\n",chan);
        fflush(stdout);
        return(-1);
    }
    *sampleRate = adcStr->sampleRate;
    FrameFree(frameH);
    FrFileIEnd(frFile);
    return 0;
}

   


int 
GetChannelData(int gpsTime, int numSecs, char *chan, float *data, int size)
{
    int frRdInFile;
    int gpsT = gpsTime;
    int secRead=0;
    int iSRate;
    FrFile *frFile=NULL;
    FrameH *frameH=NULL;
    FrAdcData *adcStr=NULL;
    float slope, bias;
    int i,ind=0;
    
    frRdInFile = gpsTime % 60;
    if(SetFileNames(gpsTime)<0){
        fprintf(stderr, "Recursive Error:GetChannelData\n");
        fflush(stdout);
        return -1;
    }
    if((frFile=FrFileINew(currFile))==NULL){
        fprintf(stderr, "Cannot Open Frame File %s\n",currFile);
        fflush(stdout);
        return -1;
    }
    while(secRead < numSecs){
        if((frameH = FrameReadTAdc(frFile,(double)gpsT,chan))==NULL){
            fprintf(stderr, "Cannot Read Frame for channel %s for gps\n",chan);
            fprintf(stderr, " time %d in file \n",gpsTime);
            fflush(stdout);
            return -1;
        }
        if((adcStr = FrAdcDataFind(frameH, chan))==NULL){
            fprintf(stderr, "Cant get pointer to adcStruct in Frame for channel %s\n",chan);
            fflush(stdout);
            return(-1);
        }
        iSRate = adcStr->sampleRate;
        if(iSRate*numSecs > size){
            fprintf(stderr, "Size of array is insufficient: number of Secs = %d\n",numSecs);
            fprintf(stderr, "SampleRate = %d and size = %d \n",iSRate, size);
            fflush(stdout);
            return -1;
        }

        slope = adcStr->slope;
        bias = adcStr->bias;

	switch (adcStr->data->type)
	{
		case FR_VECT_C:
			for(i=0;i<iSRate;i++)
			{
				data[ind++] = adcStr->data->dataU[i] * slope - bias;
			}
			break;
		case FR_VECT_2S:
			for(i=0;i<iSRate;i++)
			{
				data[ind++] = adcStr->data->dataS[i] * slope - bias;
			}
			break;
		case FR_VECT_4R:
			for(i=0;i<iSRate;i++)
			{
				data[ind++] = adcStr->data->dataF[i] * slope - bias;
			}
			break;
		case FR_VECT_4S:
			for(i=0;i<iSRate;i++)
			{
				data[ind++] = adcStr->data->dataL[i] * slope - bias;
			}
			break;
		default:
			break;
	}

        FrameFree(frameH);
        adcStr=NULL;
        frameH=NULL;
        secRead++;
        gpsT++;
        if((gpsT%60 == 0) && (secRead < numSecs)){
            FrFileIEnd(frFile);
            frFile = NULL;
            if(SetFileNames(gpsT)<0){
                fprintf(stderr, "Recursive Error:GetChannelData\n");
                fflush(stdout);
                return -1;
            }
            if((frFile=FrFileINew(currFile))==NULL){
                fprintf(stderr, "Cannot Open Frame File %s\n",currFile);
                fflush(stdout);
                return -1;
            }  
        }
    }
    FrFileIEnd(frFile);
    return 0;
}



int 
GetChannelDoubleData(int gpsTime, int numSecs, char *chan, double *data, int size)
{
    int frRdInFile;
    int gpsT = gpsTime;
    int secRead=0;
    int iSRate;
    FrFile *frFile=NULL;
    FrameH *frameH=NULL;
    FrAdcData *adcStr=NULL;
    float slope, bias;
    int i,ind=0;
    
    frRdInFile = gpsTime % 60;
    if(SetFileNames(gpsTime)<0){
        fprintf(stderr, "Recursive Error:GetChannelData\n");
        fflush(stdout);
        return -1;
    }
    if((frFile=FrFileINew(currFile))==NULL){
        fprintf(stderr, "Cannot Open Frame File %s\n",currFile);
        fflush(stdout);
        return -1;
    }
    while(secRead < numSecs){
        if((frameH = FrameReadTAdc(frFile,(double)gpsT,chan))==NULL){
            fprintf(stderr, "Cannot Read Frame for channel %s for gps\n",chan);
            fprintf(stderr, " time %d in file \n",gpsTime);
            fflush(stdout);
            return -1;
        }
        if((adcStr = FrAdcDataFind(frameH, chan))==NULL){
            fprintf(stderr, "Cant get pointer to adcStruct in Frame for channel %s\n",chan);
            fflush(stdout);
            return(-1);
        }
        iSRate = adcStr->sampleRate;
        if(iSRate*numSecs > size){
            fprintf(stderr, "Size of array is insufficient: number of Secs = %d\n",numSecs);
            fprintf(stderr, "SampleRate = %d and size = %d \n",iSRate, size);
            fflush(stdout);
            return -1;
        }

        slope = adcStr->slope;
        bias = adcStr->bias;

	switch (adcStr->data->type)
	{
		case FR_VECT_C:
			for(i=0;i<iSRate;i++)
			{
				data[ind++] = adcStr->data->dataU[i] * slope - bias;
			}
			break;
		case FR_VECT_2S:
			for(i=0;i<iSRate;i++)
			{
				data[ind++] = adcStr->data->dataS[i] * slope - bias;
			}
			break;
		case FR_VECT_4R:
			for(i=0;i<iSRate;i++)
			{
				data[ind++] = adcStr->data->dataF[i] * slope - bias;
			}
			break;
		case FR_VECT_4S:
			for(i=0;i<iSRate;i++)
			{
				data[ind++] = adcStr->data->dataL[i] * slope - bias;
			}
		case FR_VECT_8R:
			for(i=0;i<iSRate;i++)
			{
				data[ind++] = adcStr->data->dataD[i] * slope - bias;
			}
			break;
		default:
			break;
	}

        FrameFree(frameH);
        adcStr=NULL;
        frameH=NULL;
        secRead++;
        gpsT++;
        if((gpsT%60 == 0) && (secRead < numSecs)){
            FrFileIEnd(frFile);
            frFile = NULL;
            if(SetFileNames(gpsT)<0){
                fprintf(stderr, "Recursive Error:GetChannelData\n");
                fflush(stdout);
                return -1;
            }
            if((frFile=FrFileINew(currFile))==NULL){
                fprintf(stderr, "Cannot Open Frame File %s\n",currFile);
                fflush(stdout);
                return -1;
            }  
        }
    }
    FrFileIEnd(frFile);
    return 0;
}






int
GetCalibFunction(int gpsTime, char *chan, float *data, float deltaF, int size)
{
    int gpsMin;
    FrFile *frFile;
    FrameH *frameH;
    FrStatData *frStat;
    
    if(SetFileNames(gpsTime)<0){
        fprintf(stderr, "Recursive Error:GetCalibFunc\n");
        fflush(stdout);
        return -1;
    }
    gpsMin = gpsTime - gpsTime%60;
    if((frFile=FrFileINew(currFile))==NULL){
        fprintf(stderr, "Cannot Open Frame File %s\n",currFile);
        fflush(stdout);
        return -1;
    }
    if((frameH = FrameRead(frFile))==NULL){
        fprintf(stderr, "Cannot Read Frame for channel %s for gps\n",chan);
        fprintf(stderr, " time %d in file \n",gpsMin);
        fflush(stdout);
        return -1;
    }
    if((frStat = FrStatDataFind(frameH->detectProc,chan,0))==NULL){
        fprintf(stderr, "Cannot find Calibration Info for channel %s\n",chan);
        fprintf(stderr, "gpsMin = %d, currFile = %s\n",gpsMin,currFile);
        fflush(stdout);
        return -1;
    }
    if(GetCalibrationFunc(frStat->data->nx[0],
                          (const float *)frStat->data->dataF,size,
                          data,deltaF)< 0){
        fprintf(stderr, "Recursive Error: GetCalibFunc\n");
        fflush(stdout);
        return -1;
    }
    FrameFree(frameH);
    FrFileIEnd(frFile);    
    return 0;
}

static int
GetCalibrationFunc(int n, const float *f, int size, float *calib, float deltaF)
{
    
    int nSys,nWhit,nTot;
    float *gain;
    float *offset;
    float *fStart, *fEnd;
    int *nz,*np;
    Complex **zero,**pole;
    int i, j, k, ind;
    
    if(size<4){
        fprintf(stderr, "the length of the array passed for the output array \n");
        fprintf(stderr, "is too small (should be > 4), size = %d\n",size);
        fflush(stdout);
        return -1;
    }
    nSys=f[0];nWhit=f[1];
    ind=2;
    nTot = nSys+nWhit;
    gain = malloc(nTot*sizeof(float));
    offset= malloc(nTot*sizeof(float));
    nz = malloc(nTot*sizeof(int));
    np= malloc(nTot*sizeof(int));
    fStart = malloc(nTot*sizeof(float));
    fEnd = malloc(nTot*sizeof(float));
    zero = malloc(sizeof(Complex *) * nTot);
    pole = malloc(sizeof(Complex *) * nTot);
    for(i=0;i<nTot;i++){
        gain[i] = f[ind++];
        if(gain[i]==0.0) gain[i] = 1.0;
        offset[i] = f[ind++];
        np[i] = f[ind++];
        nz[i] = f[ind++];
        fStart[i] = f[ind++];
        fEnd[i] = f[ind++];
        pole[i] = malloc(np[i] * sizeof(Complex));
        for(j=0;j<np[i];j++){
            pole[i][j].r=f[ind++];
            pole[i][j].i=f[ind++];
        }
        zero[i] = malloc(nz[i] * sizeof(Complex));
        for(j=0;j<nz[i];j++){
            zero[i][j].r=f[ind++];
            zero[i][j].i=f[ind++];
        }       
    }
    if(ind !=n ){
        fprintf(stderr, "Error: Calibration Array information seems to be incorrect\n");
        fflush(stdout);
        return(-1);
    }
    for(i=0;i<size;i++)
        calib[i] = 0.0;
    {
        float freq;
        Complex prod, term,tmp;
        for(i=0;i<size/2;i++){
            freq = i*deltaF;
            for(j=0;j<nTot;j++){
                if((freq>=fStart[j])&&(freq<=fEnd[j])){
                    prod.r = gain[j];
                    prod.i = 0.0;
                    for(k=0;k<nz[j];k++){
                        term.r = -1.0 * zero[j][k].r;
                        term.i = 2.0 * PI * freq - zero[j][k].i;
                        prod = mult(&prod, &term);
                        
                    }
                    for(k=0;k<np[j];k++){
                        term.r = -1.0 * pole[j][k].r;
                        term.i = 2.*PI*freq - pole[j][k].i;
                        prod = divi(&prod, &term);
                    }
                    if((calib[2*i] == 0.0) && (calib[2*i+1] == 0.0)){
                        calib[ 2 * i] = prod.r;
                        calib[ 2 * i + 1] = prod.i;
                    }
                    else{
                        tmp.r = calib[2*i];
                        tmp.i = calib[2*i+1];
                        prod = mult ( &prod, &tmp);
                    }
                }
            }
        }
    }

    for(i=0;i<nTot;i++){
        if(nz[i]!=0)
            free(zero[i]);
        if(np[i]!=0)
            free(pole[i]);
    }
    free(gain);
    free(offset);
    free(fStart);
    free(fEnd);
    free(nz);
    free(np);
    return 0;   
}


static Complex mult (Complex *a, Complex *b)
{
    Complex tmp;

    tmp.r = a->r * b->r - a->i * b->i;
    tmp.i = a->i * b->r + a->r * b->i;
    return tmp;
}

static Complex divi (Complex *a, Complex *b)
{
    Complex tmp1,tmp2;

    tmp1.r = b->r/ ( b->r*b->r + b->i * b->i);
    tmp1.i = (-1.0 * b->i) / ( b->r*b->r + b->i * b->i);

    tmp2 = mult(a, &tmp1);
    return tmp2;
}

int 
GetWhiteNoiseData(int n, float *ch)
{
  static int entry=1;
  char *ranSeedString=NULL;
  int i;

  if(entry){
    entry=0;
    if((ranSeedString=getenv("IDUM_SEED"))!=NULL){
      if(atoi(ranSeedString)>=0){
	fprintf(stderr,"The seed should be a large negative number;Change env variable IDUM_SEED\n");
	exit(-1);
      }
      idum = atoi(ranSeedString);
    }
    else {
      fprintf(stderr,"\n**************************\n");
      fprintf(stderr,"Evironment Variable IDUM_SEED should be set\n"); 
      fprintf(stderr,"to a large negative number to\n ");
      fprintf(stderr,"change the hard coded value of Random Seed\n");
      fprintf(stderr,"\n**************************\n");
    }
  }
  for(i=0;i<n;i++)
    ch[i] = mygasdev(&idum);
  return n;
}

static 
float mygasdev(long *idum)
{
   
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  
  if  (iset == 0) {
    do {
      v1=2.0*ran2(idum)-1.0;
      v2=2.0*ran2(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
        iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX



