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
/*  #include "../include/wavelet_test_static.h" */

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

  *output=(OutputLayerWavelet*)LALMalloc(sizeof(OutputLayerWavelet));
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

  if(output==NULL || (*output)==NULL || input==NULL || input->wavelet==NULL){
        ABORT( status, LALWAVELETH_ENULLP, LALWAVELETH_MSGENULLP );
  }


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
  int i;

  INITSTATUS( status, "LALPercentileWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  if(output==NULL || (*output)==NULL || input==NULL){ 
        ABORT( status, LALWAVELETH_ENULLP, LALWAVELETH_MSGENULLP );
  }

  _createClusterWavelet(&(*output)->out);
  _assignWavelet(&((*output)->out->wavelet),input->in);
  (*output)->out->nonZeroFractionAfterPercentile=
    _percentile((*output)->out->wavelet, input->nonZeroFraction, 
		FALSE, &(*output)->out->medians, &(*output)->out->norm50);


  if(fabs((*output)->out->nonZeroFractionAfterPercentile - input->nonZeroFraction) > maxError1)
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
  REAL4TimeSeries *a, *b;

  INITSTATUS( status, "LALPixelSwapWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _assignClusterWavelet(&(*output)->out, input->in);

  M = _getMaxLayer(input->in->wavelet)+1;
  nS=input->in->wavelet->data->data->length/M;

  for(i=0; i<M; i++){
    _getLayer(&a,i,input->in->wavelet);
    _assignREAL4TimeSeries(&b,a);

    for(j=0; j<nS; j++){
      b->data->data[j]=a->data->data[nS-1-j];
    }
    _putLayer(b, i, (*output)->out->wavelet);
  }

  (*output)->out->pixelSwapApplied=TRUE;

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
  INT4 i, j, nS;
  RandomParams *rparams;
  REAL4 x;

  INITSTATUS( status, "LALPixelMixerWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  nS=input->in->wavelet->data->data->length;

  LALCreateRandomParams(status->statusPtr,&rparams,input->seed);

  _assignClusterWavelet(&(*output)->out, input->in);

  for(i=0; i<nS; i++)
    {
      (*output)->out->wavelet->data->data->data[i] = 0.0;
    }

  for(i=0; i<nS; i++)
    {
      if(input->in->wavelet->data->data->data[i] != 0.0)
	{
	  do
	    {
	      LALUniformDeviate(status->statusPtr,&x,rparams);
	      j=(INT4)(x*(nS-1));
	    }
	  while((*output)->out->wavelet->data->data->data[j] != 0.0);
	  (*output)->out->wavelet->data->data->data[j] = 
	    input->in->wavelet->data->data->data[i];
	}
    }

  (*output)->out->pixelMixerApplied=TRUE;
  LALDestroyRandomParams(status->statusPtr,&rparams);

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
  int maxLayer1,maxLayer2,k,i,j;
  REAL4TimeSeries *one=NULL;
  REAL4TimeSeries *two=NULL;
  int status1, status2, n1,n2;
  int timeWindowSteps;
  int sw,ew;
  int events=0, eventsOne=0, eventsTwo=0;

  INITSTATUS( status, "LALCoincidenceWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  if(input->one->wavelet->type!=input->two->wavelet->type){
    ABORT( status, LALWAVELETH_ETYPEMISMATCH, LALWAVELETH_MSGETYPEMISMATCH);
  }

  _assignClusterWavelet(&((*output)->one),input->one);
  _assignClusterWavelet(&((*output)->two),input->two);

  maxLayer1=_getMaxLayer(input->one->wavelet);
  maxLayer2=_getMaxLayer(input->two->wavelet);

  if(maxLayer1!=maxLayer2){

  }

  for(k=0;k<=maxLayer1;k++){
    status1=_getLayer(&one, k, input->one->wavelet);
    status2=_getLayer(&two, k, input->two->wavelet);

    if(status1!=status2 || status1==-1 || status2==-1){

    }

    n1=one->data->length;
    n2=two->data->length;

    if(n1!=n2){

    }

    timeWindowSteps=_nanoSeconds2steps(one, input->timeWindowNanoSec);

    for(i=0;i<n1;i++){
      events=0;
      if(one->data->data[i]==0) continue;
      sw=i-timeWindowSteps;
      if(sw<0) sw=0;
      ew=i+timeWindowSteps;
      if(ew>n2) ew=n2;
      for(j=sw;j<=ew;j++){
	if(fabs(two->data->data[j]) > 0){
	  events++;
	  break;
	}
      }
      if(events==0) one->data->data[i]=0;
      else eventsOne++;
    }
    _putLayer(one, k, (*output)->one->wavelet);

    status1=_getLayer(&one, k, input->one->wavelet);
    for(i=0;i<n2;i++){
      events=0;
      if(two->data->data[i]==0) continue;
      sw=i-timeWindowSteps;
      if(sw<0) sw=0;
      ew=i+timeWindowSteps;
      if(ew>n1) ew=n1;
      for(j=sw;j<=ew;j++){
	if(fabs(one->data->data[j]) > 0){
	  events++;
	  break;
	}
      }
      if(events==0) two->data->data[i]=0;
      else eventsTwo++;
    }
    _putLayer(two, k, (*output)->two->wavelet);

  }

  (*output)->one->nonZeroFractionAfterCoincidence =
    ((REAL4)eventsOne)/(*output)->one->wavelet->data->data->length;
  (*output)->two->nonZeroFractionAfterCoincidence =
    ((REAL4)eventsTwo)/(*output)->two->wavelet->data->data->length;

  DETATCHSTATUSPTR(status);
  RETURN(status);
}

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
    _setMask((*output)->w, input->minClusterSize, input->aura, input->original);
  ncluster=_clusterMain((*output)->w);
  _clusterProperties((*output)->w);

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
  INITSTATUS( status, "LALSetAmplitudesWavelet", LALWAVELETC );
  ATTATCHSTATUSPTR (status);

  _assignREAL4TimeSeries(left,right);

  DETATCHSTATUSPTR(status);
  RETURN(status);

}
