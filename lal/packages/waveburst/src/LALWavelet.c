/************************************ <lalVerbatim file="LALWaveletCV">
Author: Klimenko, Sergey and Yakushin, Igor
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{LALWavelet.c}}

%This module defines functions used by Waveburst DSO to do standard wavelet
%operations and transforms, percentile transform, coincidence, and clustering.

\subsubsection*{Prototypes}
\input{LALWaveletCP}

\subsubsection*{Description}

%\textbf{LALWavelet} is designed to ...

%\subsubsection*{Algorithm}

%Wavelet, ...

\subsubsection*{Uses}

% List any external functions called by this function.
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}



%\input{LALWaveletCTODO}

\vfill{\footnotesize\input{LALWaveletCV}}

******************************************************* </lalLaTeX> */ 

/******* INCLUDE STANDARD LIBRARY HEADERS; ************/
/* note LALStdLib.h already includes stdio.h and stdarg.h */

/******* INCLUDE ANY LDAS LIBRARY HEADERS ************/

/******* INCLUDE ANY LAL HEADERS ************/
#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALWavelet.h>
#include <lal/Random.h>
#include <lal/AVFactories.h>
#include <math.h>

/******* DEFINE RCS ID STRING ************/
NRCSID( LALWAVELETC, "$Id$" );

/******* DEFINE LOCAL CONSTANTS AND MACROS ************/

#define maxError1 0.005

/******* DECLARE LOCAL (static) FUNCTIONS ************/

#include "wavelet_static.h"


/******* DEFINE GLOBAL FUNCTIONS ************/

/********* <lalVerbatim file="LALWaveletCP"> 

********* </lalVerbatim> *****/

/********* <lalVerbatim file="LALWaveletCTODO"> 
%To compute correlation between clusters.
********* </lalVerbatim> *****/



/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALGetLayerWavelet(LALStatus *status,
		   OutputLayerWavelet **output,
		   InputLayerWavelet *input)
/******** </lalVerbatim> ********/
{
  INITSTATUS( status, "LALGetLayerWavelet", LALWAVELETC);
  ATTATCHSTATUSPTR (status);

  if(output==NULL || input==NULL || input->wavelet==NULL 
     || input->wavelet->data==NULL || input->wavelet->data->data==NULL){
        ABORT( status, LALWAVELETH_ENULLP, LALWAVELETH_MSGENULLP );
  }

  *output=(OutputLayerWavelet*)LALCalloc(1,sizeof(OutputLayerWavelet));
  (*output)->status=_getLayer(&((*output)->layer), input->index, input->wavelet);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALGetMaxLayerWavelet(LALStatus *status,
		      OutputGetMaxLayerWavelet **output,
		      InputGetMaxLayerWavelet *input)
/******** </lalVerbatim> ********/
{
  INITSTATUS( status, "LALGetMaxLayerWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  if(input==NULL || input->wavelet==NULL){
        ABORT( status, LALWAVELETH_ENULLP, LALWAVELETH_MSGENULLP );
  }
  *output=(OutputGetMaxLayerWavelet *)LALCalloc(1,sizeof(OutputGetMaxLayerWavelet));

  (*output)->maxLayer=_getMaxLayer(input->wavelet);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALPercentileWavelet( LALStatus *status,
		      OutputPercentileWavelet **output,
		      InputPercentileWavelet  *input)
/******** </lalVerbatim> ********/
{
  InputTimeIntervalWavelet in;
  OutputTimeIntervalWavelet out;

  INITSTATUS( status, "LALPercentileWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  if(output==NULL || (*output)==NULL || input==NULL){ 
        ABORT( status, LALWAVELETH_ENULLP, LALWAVELETH_MSGENULLP );
  }


  _createClusterWavelet(&(*output)->out);
  (*output)->out->nsubintervals=input->nsubintervals;

  if(input->interpolate==1)
    {
      _interpolate(input);
    }

  _calibrate(input->in, input->R, input->C, input->alpha, input->gamma, 
	     (*output)->out, input->offsetSec);

  _assignWavelet(&(*output)->out->afterCalibration,input->in);

  _whiteAlone(input->in, &(*output)->out->medians, 
	      &(*output)->out->norm50, (*output)->out->nsubintervals, input->offsetSec);

  if(input->wavefilter==1)
    {
      _waveFilter(&input->in, (*output)->out, input->offsetSec, input->extradeep, 
		  input->wf_LPFilterLength,input->wf_HPFilterLength);
    }

  in.w=input->in;
  in.offsetSec=input->offsetSec;
  in.offsetNan=0;
  in.durationSec=(int)(input->in->data->data->length*input->in->data->deltaT + 0.5) - 
    2*input->offsetSec;
  in.durationNan=0;
  in.type=TIME_INTERVAL_WAVELET;
  
  LALGetTimeIntervalWavelet(status->statusPtr,
			    &out, &in);
  CHECKSTATUSPTR(status);

  _assignWavelet(&((*output)->out->wavelet), out.w);

  _freeWavelet(&input->in);
  _assignWavelet(&input->in,out.w);

  (*output)->out->nonZeroFractionAfterPercentile=
    _percentile((*output)->out->wavelet, input->nonZeroFraction, 
		FALSE, 
		&(*output)->out->norm10L, &(*output)->out->norm10R);


  if(fabs((*output)->out->nonZeroFractionAfterPercentile - 
	  input->nonZeroFraction) > maxError1)
    {
        ABORT( status, LALWAVELETH_EDIFF, LALWAVELETH_MSGEDIFF );
    } 

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALPixelSwapWavelet(LALStatus *status,
		    OutputPixelSwapWavelet **output,
		    InputPixelSwapWavelet *input)
/******** </lalVerbatim> ********/
{
  INT4 i, j, M, nS;
  REAL4TimeSeries *a, *b, *ao, *bo;

  INITSTATUS( status, "LALPixelSwapWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _assignClusterWavelet(&(*output)->out, input->in);

  M = _getMaxLayer(input->in->wavelet)+1;
  nS=input->in->wavelet->data->data->length/M;

  for(i=0; i<M; i++){
    _getLayer(&a,i,input->in->wavelet);
    _assignREAL4TimeSeries(&b,a);
    _getLayer(&ao,i,input->in->original);

    _assignREAL4TimeSeries(&bo,ao);

    for(j=0; j<nS; j++){
      b->data->data[j]=a->data->data[(j+nS/2)%nS];
      bo->data->data[j]=ao->data->data[(j+nS/2)%nS];
    }
    _putLayer(b, i, (*output)->out->wavelet);
    _putLayer(bo,i, (*output)->out->original);
  }

  (*output)->out->pixelSwapApplied=TRUE;
  (*output)->out->clusterType=SWAPPED_CL;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALPixelMixerWavelet(LALStatus *status,
		     OutputPixelMixerWavelet **output,
		     InputPixelMixerWavelet *input)
/******** </lalVerbatim> ********/
{
  INT4 i, j, k, nS, M;
  RandomParams *rparams=input->rparams;
  REAL4 x;
  REAL4TimeSeries *a, *b, *ao, *bo;

  INITSTATUS( status, "LALPixelMixerWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _assignClusterWavelet(&(*output)->out, input->in);


/*   bzero((*output)->out->wavelet->data->data->data, */
/* 	sizeof(REAL4)*input->in->wavelet->data->data->length); */
/*   bzero((*output)->out->original->data->data->data, */
/* 	sizeof(REAL4)*input->in->original->data->data->length); */
  

  M = _getMaxLayer(input->in->wavelet)+1;
  nS=input->in->wavelet->data->data->length/M;

/*   nz=_countNonZeroes(input->in->wavelet->data); */
/*   printf("PixelMixer before: %d\n",nz);fflush(stdout); */
/*   printf("M=%d nS=%d\n",M,nS); */

  for(i=0; i<M; i++){
    _getLayer(&a,i,input->in->wavelet);
    _assignREAL4TimeSeries(&b,a);
    _getLayer(&ao,i,input->in->original);
    _assignREAL4TimeSeries(&bo,ao);

    for(j=0;j<nS;j++)
      {
	b->data->data[j]=bo->data->data[j]=0.0;
      }

    for(j=0; j<nS; j++){
      if(a->data->data[j]!=0.0)
	{
	  do
	    {
	      LALUniformDeviate(status->statusPtr,&x,rparams);
	      k=(INT4)(x*(nS-1));
	    }
	  while(b->data->data[k]!=0.0);
	  
	  b->data->data[k]=a->data->data[j];
	  bo->data->data[k]=ao->data->data[j];
	}
    }
    _putLayer(b, i, (*output)->out->wavelet);
    _putLayer(bo,i, (*output)->out->original);

    _freeREAL4TimeSeries(&a);
    _freeREAL4TimeSeries(&b);
    _freeREAL4TimeSeries(&bo);
    _freeREAL4TimeSeries(&ao);
  }

/*   nz=_countNonZeroes((*output)->out->wavelet->data); */
/*   printf("PixelMixer after: %d\n",nz);fflush(stdout); */

  (*output)->out->pixelMixerApplied=TRUE;
  (*output)->out->clusterType=MIXED_CL;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/*
void
LALFractionWavelet( LALStatus *status,
		    OutputFractionWavelet **output,
		    InputFractionWavelet  *input){
  Wavelet *pw=NULL;
  INT4 dim;
  INT4 nF;
  INT4 nL;
  INT4 i,j,k;
  REAL4 *p;
  REAL4 x,dR,d;
  long r,R;
  Slice S;
  INT4 nS,kS,nmax, ndim;
  INT4 nr, n4;
  INT4 nZero = 0;
  INT4 M;
  INT4 nLa, ncut;
  REAL4TimeSeries *a;
  REAL4TimeSeries *b;

  double test;

  INITSTATUS( status, "LALFractionWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  if(output==NULL || (*output)==NULL || input==NULL){ 
        ABORT( status, LALWAVELETH_ENULLP, LALWAVELETH_MSGENULLP );
  }

  dim=input->dim;
  nF=input->nF;
  nL=input->nL;

  if(input->seed != 0)
    {
      srandom(input->seed);
    }

  _createClusterWavelet(&(*output)->out);
  _assignWavelet(&((*output)->out->wavelet),input->in);

  pw=(*output)->out->wavelet;

  if(dim<1) dim = 1;

  M = _getMaxLayer(input->in)+1;

  nF *= dim;
  nL *= dim;

  if(nF==0) nL = 0;
  else if(nF>0 && abs(nL)<nF) nL = nF; 
  else if(nF<0 && abs(nL)<abs(nF) && nL) 
    nL = nL>0 ? -nF : nF;

  nLa  = abs(nL);
  ncut = abs(nF);


  if(nL){                         
    for(i=0; i<M; i++){
      _getLayer(&a,i,pw); 
      _assignREAL4TimeSeries(&b,a);

      nS = a->data->length;

      nr = (nS*3)/4;
      dR = nr/((REAL4)2147483647.0);
      n4 = nS/4;

      _getSliceF(i,pw,&S);

      kS = S.step;
      p = pw->data->data->data+S.start;

      for(j=0; j<nS; j++){
	if(a->data->data[j]<0.) a->data->data[j]*=-1.;
	p[j*kS] = 0;
      } 
	 
      for(j=0; j<nS; j++){
	x = a->data->data[j];

	if(x==0.) { nZero++; continue; }

	k = j-nS/8;
	nmax = nLa;
	ndim = dim;

	do { 
	  R = random();

	  d = dR;
	  r = (long)(R*d);
	  R-= (long)(r/d);

	  if(r-k > 0) r += n4;
	  if(x<=a->data->data[r]) ndim--;

	  if(!(ndim && --nmax)) break;

	  if((d*=nr) > 0.5) continue;
	  r = (long)(R*d);
	  R-= (long)(r/d);
	  if(r-k > 0) r += n4;
	  if(x<=a->data->data[r]) ndim--;
	  if(!(ndim && --nmax)) break;

	  if((d*=nr) > 0.5) continue;
	  r = (long)(R*d);


	  R-= (long)(r/d);

	  if(r-k > 0) r += n4;

	  if(x<=a->data->data[r]) ndim--;

	}  while(ndim && --nmax);


	nmax = nLa - nmax + 1;
	if(nmax<ncut) { nZero++; continue; }

	if(nLa > abs(nF)){
	  x = b->data->data[j]>0 ? nmax-dim : dim-nmax;
	  x /= (REAL4)(dim);     
	}
	else x = b->data->data[j];   

	if(nF>0){                
	  p[j*kS] = x;
	  
	  continue;
	}

	else if(nL>0){           
	  r = j;                 
	  if(nF < 0)  r += nS/2; 
	  if(r >= nS) r -= nS;
	  p[r*kS] = x;

	  continue;
	}

	else {                  
	  d=nS/((REAL4)2147483647.0);
	  do{ if((r=(INT4)(d*random())) == nS) r--; }
	  while(p[r*kS] != 0);
	  p[r*kS] = x; 

	}
      }
      _freeREAL4TimeSeries(&a);
      _freeREAL4TimeSeries(&b);
    }
  }

  else if(nF){                
    M = pw->data->data->length;
    x = ((REAL4)(dim))/((REAL4)(abs(nF)));
    for(i=0; i<M; i++)
      if(drand48() > x) { pw->data->data->data[i] = 0; nZero++; }

   }

   else{                        
      M = pw->data->data->length;
      for(i=0; i<M; i++)
	if(pw->data->data->data[i]==0) nZero++;

   }
  
  (*output)->out->nonZeroFractionAfterPercentile = 
    ((REAL4)(pw->data->data->length-nZero))/((REAL4)(pw->data->data->length));

  DETATCHSTATUSPTR(status);
  RETURN(status);
}
*/


/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALCoincidenceWavelet(LALStatus *status,
		      OutputCoincidenceWavelet **output,
		      InputCoincidenceWavelet *input)
/******** </lalVerbatim> ********/
{
  int i,j,k,l,m,n;
  int status1, status2;
  int swt,ewt,swf,ewf;
  int length;
  int events=0;
  int tw = input->timeWindowPixels;
  int fw = input->freqWindowPixels;
  REAL4TimeSeries *one=NULL;
  REAL4TimeSeries *two=NULL;
  REAL4 **one2D=NULL;
  REAL4 **two2D=NULL;
  REAL4 aone, atwo, energy;
  BOOLEAN CROSS = input->coincidenceLevel==CROSS_CO || input->coincidenceLevel==STRICT_CROSS_CO;
  BOOLEAN STRICT= input->coincidenceLevel==STRICT_CROSS_CO || input->coincidenceLevel==STRICT_BOX_CO;
  int pixels;

  INITSTATUS( status, "LALCoincidenceWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  if(input->one->wavelet->type!=input->two->wavelet->type)
    {
      ABORT( status, LALWAVELETH_ETYPEMISMATCH, LALWAVELETH_MSGETYPEMISMATCH);
    }

  _assignClusterWavelet(&((*output)->one),input->one);
  _assignClusterWavelet(&((*output)->two),input->two);
  
  length=input->one->wavelet->data->data->length;
  
  if(_getMaxLayer(input->one->wavelet) != _getMaxLayer(input->two->wavelet)) { /* dummy */ }
  m = _getMaxLayer(input->one->wavelet);
  n = (int)(length/(m+1.)+0.5)-1;

  _getSpectrogram(input->one->wavelet, &one2D);
  _getSpectrogram(input->two->wavelet, &two2D);
  
  /* cleanup single pixel clusters in both channels*/

  
  for(k=0; k<=m; k++) 
    {

      status1=_getLayer(&one, k, input->one->wavelet);
      status2=_getLayer(&two, k, input->two->wavelet);
    
      if(status1!=status2 || status1==-1 || status2==-1) 
	{
	  fprintf(stderr,"LALCoincidence: bad getLayer status: status1=%d, status2=%d\n",status1, status2);
	  fflush(stderr);
	}
    
      if((int)one->data->length != (int)n+1 || (int)two->data->length != (int)n+1 ) 
	{
	  fprintf(stderr,"LALCoincidence:n1!=n2, n1=%d, n2=%d\n", one->data->length, two->data->length);
	  fflush(stderr);
	}
    
      for(i=0; i<=n; i++) 
	{
	  aone = one2D[k][i];
	  atwo = two2D[k][i];
	  if(aone==0. && atwo==0.) continue;
      
	  swt = i>0 ? i-1 : 0;
	  swf = k>0 ? k-1 : 0;
	  ewt = i<n ? i+1 : n;
	  ewf = k<m ? k+1 : m;
	  
	  /* clean up first channel*/
	  events = 0;
	  if(aone != 0 && fabs(aone)<input->minAmp4SinglePixels) 
	    {
	      for(j=swt; j<=ewt; j++)
		for(l=swf; l<=ewf; l++)
		  if(!(j==i && l==k) && one2D[l][j]!=0.) { events++; break; }
	      if(!events) one->data->data[i] = 0.0;
	    }      

	  /* clean up second channel*/
	  events = 0;
	  if(atwo != 0 && fabs(atwo)<input->minAmp4SinglePixels)
	    {
	      for(j=swt; j<=ewt; j++)
		for(l=swf; l<=ewf; l++)
		  if(!(j==i && l==k) && two2D[l][j]!=0.) { events++; break; } 
	      if(!events) two->data->data[i] = 0.0;
	    }      
      
	}

      _putLayer(one, k, (*output)->one->wavelet);
      _putLayer(two, k, (*output)->two->wavelet);	
    
      _freeREAL4TimeSeries(&one);
      _freeREAL4TimeSeries(&two);
    }


  _freeSpectrogram(input->one->wavelet, &one2D);
  _freeSpectrogram(input->two->wavelet, &two2D);
  

  /* do coincidence*/

  if(input->coincidenceLevel >= 0) 
    {

      /* prepare spectrograms*/
      _getSpectrogram((*output)->one->wavelet, &one2D);
      _getSpectrogram((*output)->two->wavelet, &two2D);
    
      for(k=0; k<=m; k++) 
	{
	  _getLayer(&one, k, (*output)->one->wavelet);
	  _getLayer(&two, k, (*output)->two->wavelet);

	  for(i=0; i<=n; i++) 
	    {
	
	      aone = one2D[k][i];
	      atwo = two2D[k][i];
	      if(aone==0. && atwo==0.) continue;            /* nothing to do*/
	      if(aone!=0. && atwo!=0. && STRICT) continue;  /* strict coincidence*/
	
	      swt = i-tw<0 ? 0 : i-tw;
	      swf = k-fw<0 ? 0 : k-fw;
	      ewt = i+tw>n ? n : i+tw;
	      ewf = k+fw>m ? m : k+fw;
	      
	      energy = 0.;
	      pixels=0;
	      if(one->data->data[i] != 0.) 
		{
		  for(j=swt; j<=ewt; j++) 
		    for(l=swf; l<=ewf; l++) 
		      {
			if(CROSS && !(j==i || l==k)) continue;
			energy += fabs(two2D[l][j]);
			pixels++;
		      }
		  /*
		    energy=energy*sqrt(3)+1.69*pixels;
		    if(energy < input->minAmp4ClusterExtension*
		    input->minAmp4ClusterExtension) one->data->data[i] = 0.;
		  */
		  if(energy < input->minAmp4ClusterExtension) one->data->data[i] = 0.;
		}
	
	      energy = 0.;	
	      pixels = 0;
	      if(two->data->data[i] != 0.) 
		{
		  for(j=swt; j<=ewt; j++) 
		    for(l=swf; l<=ewf; l++)
		      { 
			if(CROSS && !(j==i || l==k)) continue;
			energy += fabs(one2D[l][j]);
			pixels++;
		      }
		  /*
		    energy=energy*sqrt(3)+1.69*pixels;
		    if(energy < input->minAmp4ClusterExtension*
		    input->minAmp4ClusterExtension) two->data->data[i] = 0.;
		  */
                  if(energy < input->minAmp4ClusterExtension) two->data->data[i] = 0.;
		}

	    }

	  _putLayer(one, k, (*output)->one->wavelet);
	  _putLayer(two, k, (*output)->two->wavelet);	
	  
	  _freeREAL4TimeSeries(&one);
	  _freeREAL4TimeSeries(&two);
	}

      _freeSpectrogram((*output)->one->wavelet, &one2D);
      _freeSpectrogram((*output)->two->wavelet, &two2D);
      
    }

  /* cleanup single pixel clusters and set the same pattern in both channels*/

  _getSpectrogram((*output)->one->wavelet, &one2D);
  _getSpectrogram((*output)->two->wavelet, &two2D);
  
  for(k=0; k<=m; k++) 
    {

      _getLayer(&one, k, (*output)->one->wavelet);
      _getLayer(&two, k, (*output)->two->wavelet);
    
      for(i=0;i<=n;i++) 
	{
      
	  aone = one2D[k][i];
	  atwo = two2D[k][i];
	  if(aone==0 && atwo==0) continue;
	  
	  swt = i>0 ? i-1 : 0;
	  swf = k>0 ? k-1 : 0;
	  ewt = i<n ? i+1 : n;
	  ewf = k<m ? k+1 : m;
	  
	  /* clean up first channel*/
	  events = 0;
	  if(aone!=0 && fabs(aone)<input->minAmp4SinglePixels) 
	    {
	      for(j=swt;j<=ewt;j++)
		for(l=swf;l<=ewf;l++)
		  if(!(j==i && l==k) && one2D[l][j]!=0.) { events++; break; }
	      if(!events) one->data->data[i] = 0.0;
	    }
      
	  /* clean up second channel*/
	  events = 0;
	  if(atwo!=0 && fabs(atwo) < input->minAmp4SinglePixels) 
	    {
	      for(j=swt; j<=ewt; j++)
		for(l=swf; l<=ewf; l++)
		  if(!(j==i && l==k) && two2D[l][j]!=0.) { events++; break; }
	      if(!events) two->data->data[i] = 0.0;
	    }
      
	  /* set a partner pixel in second channel if zero */
	  if(one->data->data[i]!=0 && two->data->data[i]==0) two->data->data[i] = 2*length;
      
	  /* set a partner pixel in first channel if zero*/
	  if(two->data->data[i]!=0 && one->data->data[i]==0) one->data->data[i] = 2*length;
	  
	}
      _putLayer(one, k, (*output)->one->wavelet);
      _putLayer(two, k, (*output)->two->wavelet);	
      
      _freeREAL4TimeSeries(&one);
      _freeREAL4TimeSeries(&two);
    }
  
  _freeSpectrogram((*output)->one->wavelet,&one2D);
  _freeSpectrogram((*output)->two->wavelet,&two2D);
  
  (*output)->one->nonZeroFractionAfterCoincidence=
    ((REAL4)_countNonZeroes((*output)->one->wavelet->data))/(*output)->one->wavelet->data->data->length;
  (*output)->two->nonZeroFractionAfterCoincidence=
    ((REAL4)_countNonZeroes((*output)->two->wavelet->data))/(*output)->two->wavelet->data->data->length;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}







/*


void
LALCoincidenceWavelet(LALStatus *status,
		      OutputCoincidenceWavelet **output,
		      InputCoincidenceWavelet *input)

{
  int maxLayer1,maxLayer2, maxLayer, k,i,j,l;
  REAL4TimeSeries *one=NULL;
  REAL4TimeSeries *two=NULL;
  int status1, status2, n1,n2,n=0;
  int swt,ewt,swf,ewf;
  int events=0, eventsOne=0, eventsTwo=0;
  int length;
  REAL4 **one2D;
  REAL4 **two2D;
  int neighbors;

  INITSTATUS( status, "LALCoincidenceWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  if(input->one->wavelet->type!=input->two->wavelet->type){
    ABORT( status, LALWAVELETH_ETYPEMISMATCH, LALWAVELETH_MSGETYPEMISMATCH);
  }

  _assignClusterWavelet(&((*output)->one),input->one);
  _assignClusterWavelet(&((*output)->two),input->two);

  length=input->one->wavelet->data->data->length;

  if(input->coincidenceLevel>=0)
    {

      maxLayer1=_getMaxLayer(input->one->wavelet);
      maxLayer2=_getMaxLayer(input->two->wavelet);
      
      if(maxLayer1!=maxLayer2)
	{
	
	}
      maxLayer=maxLayer1;

      _getSpectrogram(input->one->wavelet, &one2D);
      _getSpectrogram(input->two->wavelet, &two2D);


      for(k=0;k<=maxLayer1;k++)
	{
	  status1=_getLayer(&one, k, input->one->wavelet);
	  status2=_getLayer(&two, k, input->two->wavelet);
	
	  if(status1!=status2 || status1==-1 || status2==-1)
	    {
	    
	    }
	
	  n1=one->data->length;
	  n2=two->data->length;
	
	  if(n1!=n2)
	    {
	      fprintf(stderr,"LALCoincidence:n1!=n2, n1=%d, n2=%d\n",n1,n2);
	      fflush(stderr);
	    }
	  n=n1;


	  for(i=0;i<n;i++)
	    {
	      events=0;
	      if(one2D[k][i]==0) continue;
	      if(two2D[k][i]!=0)
		{
		  eventsOne++;
		  continue;
		}
	      if(fabs(one2D[k][i])<input->minAmp4ClusterExtension)
		{
		  one->data->data[i]=0;
		  continue;
		}
	      
	      swt=i-input->timeWindowPixels;
	      swf=k-input->freqWindowPixels;
	      if(swt<0) swt=0;
	      if(swf<0) swf=0;
	      ewt=i+input->timeWindowPixels;
	      ewf=k+input->freqWindowPixels;
	      if(ewt>=n) ewt=n-1;
	      if(ewf>maxLayer) ewf=maxLayer; 
	      for(j=swt;j<=ewt;j++)
		{
		  for(l=swf;l<=ewf;l++)
		    {
		      if(fabs(two2D[l][j]) > 0)
			{
			  if(input->coincidenceLevel==CROSS_CO && (j==i || l==k))
			    {
			      events++;
			      break;
			    }
			  else if(input->coincidenceLevel==RECTANGLE_CO)
			    {
			      events++;
			      break;
			    }			 
			}
		    }		  
		}

	      if(events==0) one->data->data[i]=0;
	      else
		{
		  eventsOne++;
		  if(two->data->data[i]==0.0)
		    {
		      two->data->data[i]=2*length;
		    }
		}
	    }

	  
	  for(i=0;i<n;i++)
	    {

	      events=0;
	      if(two2D[k][i]==0) continue;
              if(one2D[k][i]!=0)
                {
                  eventsTwo++;
                  continue;
                }
	      if(fabs(two2D[k][i])<input->minAmp4ClusterExtension)
                {
		  two->data->data[i]=0;
                  continue;
                }

	      swt=i-input->timeWindowPixels;
	      swf=k-input->freqWindowPixels;
	      if(swt<0) swt=0;
	      if(swf<0) swf=0;
	      ewt=i+input->timeWindowPixels;
	      ewf=k+input->freqWindowPixels;
	      if(ewt>=n) ewt=n-1;
	      if(ewf>maxLayer) ewf=maxLayer;  
	      for(j=swt;j<=ewt;j++)
		{
		  for(l=swf;l<=ewf;l++)
		    {
		      if(fabs(one2D[l][j]) > 0)
			{
			  if(input->coincidenceLevel==CROSS_CO && (j==i || l==k))
			    {
			      events++;
			      break;
			    }
			  else if(input->coincidenceLevel==RECTANGLE_CO)
			    {
			      events++;
			      break;
			    }
			}
		    }			 
		}
	      
	      if(events==0) two->data->data[i]=0;
	      else
		{
		  eventsTwo++;
		  if(one->data->data[i]==0.0)
		    {
		      one->data->data[i]=2*length;
		    }
		}
	    }
	  _putLayer(one, k, (*output)->one->wavelet);
	  _putLayer(two, k, (*output)->two->wavelet);	

	  _freeREAL4TimeSeries(&one);
	  _freeREAL4TimeSeries(&two);
	}

      _freeSpectrogram(input->one->wavelet, &one2D);
      _freeSpectrogram(input->two->wavelet, &two2D);


      _getSpectrogram((*output)->one->wavelet, &one2D);
      _getSpectrogram((*output)->two->wavelet, &two2D);


      for(k=0;k<=maxLayer;k++)
	{
	  status1=_getLayer(&one,k,(*output)->one->wavelet);
	  status2=_getLayer(&two,k,(*output)->two->wavelet);


	  for(i=0;i<n;i++)
	    {
	      if(one2D[k][i]==0.0 && two2D[k][i]==0.0) continue;
	      swt=i-1;
	      swf=k-1;
	      if(swt<0) swt=0;
	      if(swf<0) swf=0;
	      ewt=i+1;
	      ewf=k+1;
	      if(ewt>=n) ewt=n-1;
	      if(ewf>maxLayer) ewf=maxLayer;
	      neighbors=0;


	      for(j=swt;j<=ewt;j++)
		{
		  for(l=swf;l<=ewf;l++)
		    {
		      if(!(j==i && l==k) && one2D[l][j]) neighbors++;
		    }
		}


	      if(!neighbors)
		{
		  if(fabs(one2D[k][i]) < input->minAmp4SinglePixels ||
		     one2D[k][i]>=2*length-1 ||
		     fabs(two2D[k][i]) < input->minAmp4SinglePixels ||
		     two2D[k][i]>=2*length-1)
		    {
		      one->data->data[i]=0.0;
		      two->data->data[i]=0.0;
		    }
		}



	    }
          _putLayer(one, k, (*output)->one->wavelet);
          _putLayer(two, k, (*output)->two->wavelet);



	  _freeREAL4TimeSeries(&one);
          _freeREAL4TimeSeries(&two);


	}



      _freeSpectrogram((*output)->one->wavelet,&one2D);
      _freeSpectrogram((*output)->two->wavelet,&two2D);



    }


  (*output)->one->nonZeroFractionAfterCoincidence=
    ((REAL4)_countNonZeroes((*output)->one->wavelet->data))/(*output)->one->wavelet->data->data->length;
  (*output)->two->nonZeroFractionAfterCoincidence=
    ((REAL4)_countNonZeroes((*output)->two->wavelet->data))/(*output)->two->wavelet->data->data->length;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

*/
/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALClusterWavelet(LALStatus *status,
		  OutputClusterWavelet **output,
		  InputClusterWavelet *input)
/******** </lalVerbatim> ********/
{

  int ncluster;

  INITSTATUS( status, "LALClusterWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  if(output==NULL || (*output)==NULL || input==NULL || input->w==NULL){
        ABORT( status, LALWAVELETH_ENULLP, LALWAVELETH_MSGENULLP );
  }

  _assignClusterWavelet(&((*output)->w),input->w);

  (*output)->w->nonZeroFractionAfterSetMask = 
    _setMask((*output)->w, input->minClusterSize, input->aura);

  ncluster=_clusterMain((*output)->w);

  _clusterProperties((*output)->w);

  (*output)->w->nonZeroFractionAfterClustering=
    ((REAL4)_countNonZeroes((*output)->w->wavelet->data))/(*output)->w->wavelet->data->data->length;

  DETATCHSTATUSPTR(status);
  RETURN(status);

}



/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALReuseClusterWavelet(LALStatus *status,
		       OutputClusterWavelet **output,
		       InputReuseClusterWavelet *input)
     /******** </lalVerbatim> ********/
{
  INITSTATUS( status, "LALReuseClusterWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  if(output==NULL || (*output)==NULL || input==NULL || input->w==NULL){
    ABORT( status, LALWAVELETH_ENULLP, LALWAVELETH_MSGENULLP );
  }

  _assignClusterWavelet(&((*output)->w),input->w);

  (*output)->w->nonZeroFractionAfterSetMask =
    _duplicateClusterStructure(*output, input);

  _clusterProperties((*output)->w);

  (*output)->w->nonZeroFractionAfterClustering=
    ((REAL4)_countNonZeroes((*output)->w->wavelet->data))/(*output)->w->wavelet->data->data->length;

  DETATCHSTATUSPTR(status);
  RETURN(status);

}


/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALAllocateWavelet(LALStatus *status,
		   Wavelet **wavelet)
/******** </lalVerbatim> ********/
{
  BOOLEAN result;

  INITSTATUS( status, "LALAllocateWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  result=_allocateWavelet(wavelet);

  if(!result)
    {

    }

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALFreeWavelet(LALStatus *status,
	       Wavelet **wavelet)
/******** </lalVerbatim> ********/
{
  INITSTATUS( status, "LALFreeWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _freeWavelet(wavelet);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALFreeREAL4TimeSeries(LALStatus *status,
		       REAL4TimeSeries **t)
/******** </lalVerbatim> ********/
{

  INITSTATUS( status, "LALFreeREAL4TimeSeries", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _freeREAL4TimeSeries(t);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALFreeREAL8TimeSeries(LALStatus *status,
                       REAL8TimeSeries **t)
     /******** </lalVerbatim> ********/
{

  INITSTATUS( status, "LALFreeREAL8TimeSeries", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _freeREAL8TimeSeries(t);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


void
LALFreeREAL4FrequencySeries(LALStatus *status,
                            REAL4FrequencySeries **f)
{

  INITSTATUS( status, "LALFreeREAL4FrequencySeries", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _freeREAL4FrequencySeries(f);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

void
LALFreeCOMPLEX8FrequencySeries(LALStatus *status,
			       COMPLEX8FrequencySeries **f)
{

  INITSTATUS( status, "LALFreeCOMPLEX8FrequencySeries", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _freeCOMPLEX8FrequencySeries(f);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}


/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALFreeClusterWavelet(LALStatus *status,
		      ClusterWavelet **w)
/******** </lalVerbatim> ********/
{
  INITSTATUS( status, "LALFreeClusterWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _freeClusterWavelet(w);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALFreeOutPercentile(LALStatus *status,
		     OutputPercentileWavelet **p)
/******** </lalVerbatim> ********/
{
  INITSTATUS( status, "LALFreeOutPercentile", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _freeOutPercentile(p);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALFreeOutPixelSwap(LALStatus *status,
                    OutputPixelSwapWavelet **ps)
{
  INITSTATUS( status, "LALFreeOutPixelSwap", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _freeOutPixelSwap(ps);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}


/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALFreeOutCoincidence(LALStatus *status,
		      OutputCoincidenceWavelet **co)
/******** </lalVerbatim> ********/
{
  INITSTATUS( status, "LALFreeOutCoincidence", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _freeOutCoincidence(co);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALFreeOutCluster(LALStatus *status,
		  OutputClusterWavelet **cl)
/******** </lalVerbatim> ********/
{
  INITSTATUS( status, "LALFreeOutCluster", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _freeOutCluster(cl);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALSetAmplitudesWavelet(LALStatus *status,
			ClusterWavelet *w)
/******** </lalVerbatim> ********/
{
  INITSTATUS( status, "LALSetAmplitudesWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _setAmplitudes(w);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/******** <lalVerbatim file="LALWaveletCP"> ********/
void
LALAssignREAL4TimeSeries(LALStatus *status,
			 REAL4TimeSeries **left,
			 REAL4TimeSeries *right)
/******** </lalVerbatim> ********/
{
  INITSTATUS( status, "LALAssignREAL4TimeSeries", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _assignREAL4TimeSeries(left,right);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}

void
LALAssignREAL4FrequencySeries(LALStatus *status,
			      REAL4FrequencySeries **left,
			      REAL4FrequencySeries *right)
{
  INITSTATUS( status, "LALAssignREAL4FrequencySeries", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _assignREAL4FrequencySeries(left,right);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


void
LALAssignCOMPLEX8FrequencySeries(LALStatus *status,
                              COMPLEX8FrequencySeries **left,
                              COMPLEX8FrequencySeries *right)
{
  INITSTATUS( status, "LALAssignCOMPLEX8FrequencySeries", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _assignCOMPLEX8FrequencySeries(left,right);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}





void
LALAssignWavelet(LALStatus *status,
                 Wavelet **left,
		 Wavelet*right)
{
  INITSTATUS( status, "LALSetAmplitudesWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _assignWavelet(left,right);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALAssignClusterWavelet(LALStatus *status,
			ClusterWavelet**left,
			ClusterWavelet*right)
{
  INITSTATUS( status, "LALSetAmplitudesWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _assignClusterWavelet(left,right);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}



void LALForwardWavelet(LALStatus *status,
		       InputForwardWavelet *input)
{
  INITSTATUS( status, "LALForwardWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  if(input->w->PForward==NULL && input->w->pLForward==NULL) _setFilter(input->w);
  _forward(input->w, input->level, input->layer);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void LALInverseWavelet(LALStatus *status,
		       InputInverseWavelet *input)
{
  INITSTATUS( status, "LALInverseWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  if(input->w->PForward==NULL && input->w->pLForward==NULL) _setFilter(input->w);
  _inverse(input->w, input->level, input->layer);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void LALt2wWavelet(LALStatus *status,
		   Inputt2wWavelet *input,
		   Outputt2wWavelet **output)
{
  
  int maxLevel;
  int levs;
  int levf;
  int layf;
  int level;
  int layer;

  INITSTATUS( status, "LALt2wWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _assignWavelet(&((*output)->w),input->w);

  if((*output)->w->PForward==NULL && 
    (*output)->w->pLForward==NULL) _setFilter((*output)->w);

  maxLevel = _getMaxLevel((*output)->w,(*output)->w->data->data->length);
  
  levs = (*output)->w->level;
  levf = (*output)->w->level+input->ldeep;
  if((input->ldeep == -1) || (levf > maxLevel)) levf = maxLevel;
  
  for(level=levs; level<levf; level++)
    {
      layf = ((*output)->w->treeType==1) ? 1<<level : 1;
      
      for(layer=0; layer<layf; layer++)
	_forward((*output)->w,level,layer);
      
      (*output)->w->level=level+1;
      
    }
  
  (*output)->w->level=levf;
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void LALw2tWavelet(LALStatus *status,
		   Inputw2tWavelet *input,
		   Outputw2tWavelet **output)
{
  int levs;
  int levf;
  int layf;
  int level;
  int layer;

  INITSTATUS( status, "LALw2tWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _assignWavelet(&((*output)->w),input->w);

  if((*output)->w->PForward==NULL &&
     (*output)->w->pLForward==NULL) _setFilter((*output)->w);

  levs = (*output)->w->level;
  levf = (*output)->w->level - input->ldeep;
  if((input->ldeep == -1) || (levf < 0)) levf = 0;
  
  for(level=levs-1; level>=levf; level--)
    {
      layf = ((*output)->w->treeType==1) ? 1<<level : 1;
      
      for(layer=0; layer<layf; layer++)
	{
	  /*
	  printf("Before: level=%d levs=%d ldeep=%d levf=%d layf=%d layer=%d\n",
		 level,levs,input->ldeep,levf,layf,layer);fflush(stdout);
	  */
	  _inverse((*output)->w,level,layer);
	  /*
	  printf("After:\n");fflush(stdout);
	  */
	}
      
      
      (*output)->w->level=level;
      
    }
  (*output)->w->level=levf;
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}


void
LALGetTimeIntervalWavelet(LALStatus *status,
                          OutputTimeIntervalWavelet *output,
                          InputTimeIntervalWavelet *input)
{
  UINT4 M,i,levels;
  LALTimeInterval interval;

  INITSTATUS( status, "LALGetTimeIntervalWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  M=_getMaxLayer(input->w);  
  levels=input->w->data->data->length/(M+1);

  if(input->w->data->data->length%(M+1)!=0)
    {
      ABORT(status, LALWAVELETH_ENONZEROREMAINDER, LALWAVELETH_MSGENONZEROREMAINDER);
    }

  if(input->type==TIME_INTERVAL_WAVELET)
    {
      input->offsetSteps=
	(int)((input->offsetSec+
	       input->offsetNan*pow(10,-9))/input->w->data->deltaT + 0.5);
      input->durationSteps=
	(int)((input->durationSec + 
	       input->durationNan*pow(10,-9))/input->w->data->deltaT + 0.5);
    }
  if(input->type==STEP_INTERVAL_WAVELET)
    {
      input->offsetSec=(int)(input->offsetSteps*input->w->data->deltaT + 0.5);
      input->offsetNan=
	(int)((input->offsetSteps*input->w->data->deltaT - 
	       input->offsetSec)*pow(10,9) + 0.5);
    }
 
  if( (input->offsetSteps + input->durationSteps)  > 
      input->w->data->data->length )
    {
      char mess[LALNameLength];
      sprintf(mess,"offsetSteps=%d durationSteps=%d length=%d\n",
              input->offsetSteps, input->durationSteps,
	      input->w->data->data->length);
      /*      ABORT(status, LALWAVELETH_EOUTOFBOUNDS, LALWAVELETH_MSGEOUTOFBOUNDS);*/
      ABORT(status, LALWAVELETH_EOUTOFBOUNDS, mess);
    }

  interval.seconds=input->offsetSec;
  interval.nanoSeconds=input->offsetNan;

  _assignWavelet(&output->w,input->w);
  if(output->w->data->data->data!=NULL) 
    {
      LALFree(output->w->data->data->data);
      output->w->data->data->data=NULL;
    }
  output->w->data->data->length=input->durationSteps;
  output->w->data->data->data=(REAL4*)LALCalloc(input->durationSteps,sizeof(REAL4));

  for(i=0;i<input->durationSteps;i++)
    {
      output->w->data->data->data[i]=input->w->data->data->data[i+input->offsetSteps];
    }

  LALIncrementGPS(status->statusPtr, &output->w->data->epoch, &input->w->data->epoch, &interval);
  CHECKSTATUSPTR(status);

  /*
  printf("GetTimeInterval: deltaT=%f offsetSec=%d offsetNan=%d offsetSteps=%d durationSec=%d durationNan=%d durationSteps=%d\n",
	 input->w->data->deltaT, input->offsetSec, input->offsetNan, input->offsetSteps, input->durationSec, input->durationNan,
	 input->durationSteps);
  printf("GetTimeInterval: startTimeSec=%d startTimeNanSec=%d\n",output->w->data->epoch.gpsSeconds, output->w->data->epoch.gpsNanoSeconds);
  fflush(stdout);
  */
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void LALAddTSToWavelet(LALStatus *status,
                       LALAddTSToWaveletIO *inout)
{

  UINT4 M,i;
  UINT4 offsetSteps, durationSteps;
  REAL8 offset;
  REAL8 start1, start2;
  REAL8 waveletstep;

  INITSTATUS( status, "LALAddTSToWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  M=_getMaxLayer(inout->w)+1;

  if(inout->w->data->data->length%M!=0)
    {
      ABORT(status, LALWAVELETH_ENONZEROREMAINDER, LALWAVELETH_MSGENONZEROREMAINDER);
    }

  LALGPStoFloat(status->statusPtr, &start1, &inout->w->data->epoch);
  CHECKSTATUSPTR(status);

  LALGPStoFloat(status->statusPtr, &start2, &inout->ts->epoch);
  CHECKSTATUSPTR(status);

  waveletstep=inout->ts->deltaT*M;
  offset=(start2-start1);
  
  offsetSteps=(int)(offset/waveletstep)*M;
  durationSteps=inout->ts->data->length;

  if(offsetSteps+durationSteps>=inout->w->data->data->length)
    {
      fprintf(stderr, 
	      "AddTStoWavelet: beyond array bounds offsetSteps=%d durationSteps=%d length=%d\n",
	      offsetSteps,durationSteps,inout->w->data->data->length);
      exit(1);
    }

  /*
  printf("%d %d %d %f %f %f %f\n",
	 offsetSteps, durationSteps, M, waveletstep, start1, start2, offset);
  fflush(stdout);
  */

  for(i=0;i<durationSteps;i++)
    {
      inout->w->data->data->data[i+offsetSteps]+=inout->ts->data->data[i]*inout->strain;
    }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALUnitCopy(LALStatus *status, LALUnit *source, LALUnit *destination)
{
  INITSTATUS( status, "LALUnitCopy", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _unitCopy(source,destination);

  DETATCHSTATUSPTR(status);
  RETURN(status);  
} 

void
LALFreeOutPixelMixer(LALStatus*status,
                     OutputPixelMixerWavelet **pm)
{
  INITSTATUS( status, "LALFreeOutPixelMixer", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _freeOutPixelMixer(pm);

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

