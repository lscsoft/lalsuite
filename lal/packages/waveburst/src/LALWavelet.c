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
  int maxLayer1,maxLayer2, maxLayer, k,i,j,l;
  REAL4TimeSeries *one=NULL;
  REAL4TimeSeries *two=NULL;
  int status1, status2, n1,n2,n;
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

/*       printf("one\n");fflush(stdout); */

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

/* 	  printf("two: k=%d\n",k);fflush(stdout); */

	  for(i=0;i<n;i++)
	    {
/* 	      printf("three1: i=%d\n",i);fflush(stdout); */
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
/* 	      printf("three2: i=%d\n",i);fflush(stdout); */
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

      /*      printf("four\n");fflush(stdout); */


      _getSpectrogram((*output)->one->wavelet, &one2D);
      _getSpectrogram((*output)->two->wavelet, &two2D);

      /*      printf("five\n");fflush(stdout);*/

      for(k=0;k<=maxLayer;k++)
	{
	  status1=_getLayer(&one,k,(*output)->one->wavelet);
	  status2=_getLayer(&two,k,(*output)->two->wavelet);

	  /*	  printf("five A\n");fflush(stdout);*/

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

	      /*	      printf("six\n");fflush(stdout);*/

	      for(j=swt;j<=ewt;j++)
		{
		  for(l=swf;l<=ewf;l++)
		    {
		      if(!(j==i && l==k) && one2D[l][j]) neighbors++;
		    }
		}

	      /*	      printf("seven\n");fflush(stdout);*/

	      if(!neighbors)
		{
		  if((fabs(one2D[k][i]) < input->minAmp4SinglePixels || one2D[k][i]>=2*length-1) && 
		     (fabs(two2D[k][i]) < input->minAmp4SinglePixels || two2D[k][i]>=2*length-1))
		    {
		      one->data->data[i]=0.0;
		      two->data->data[i]=0.0;
		    }
		}

	    }
          _putLayer(one, k, (*output)->one->wavelet);
          _putLayer(two, k, (*output)->two->wavelet);


	  /*	  printf("eight\n");fflush(stdout);*/

	  _freeREAL4TimeSeries(&one);
          _freeREAL4TimeSeries(&two);

	  /*	  printf("nine\n");fflush(stdout);*/
	}

      /*      printf("ten\n");fflush(stdout);*/

      _freeSpectrogram((*output)->one->wavelet,&one2D);
      _freeSpectrogram((*output)->two->wavelet,&two2D);


      /*      printf("eleven\n");fflush(stdout);*/
    }


  (*output)->one->nonZeroFractionAfterCoincidence=
    ((REAL4)_countNonZeroes((*output)->one->wavelet->data))/(*output)->one->wavelet->data->data->length;
  (*output)->two->nonZeroFractionAfterCoincidence=
    ((REAL4)_countNonZeroes((*output)->two->wavelet->data))/(*output)->two->wavelet->data->data->length;

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

  levs = (*output)->w->level;
  levf = (*output)->w->level - input->ldeep;
  if((input->ldeep == -1) || (levf < 0)) levf = 0;
  
  for(level=levs-1; level>=levf; level--)
    {
      layf = ((*output)->w->treeType==1) ? 1<<level : 1;
      
      for(layer=0; layer<layf; layer++)
	_inverse((*output)->w,level,layer);
      
      (*output)->w->level=level;
      
    }
  (*output)->w->level=levf;
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}


