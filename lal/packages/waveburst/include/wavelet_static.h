#ifndef WAVELET_STATIC_H
#define WAVELET_STATIC_H

#include <math.h>
#include <strings.h>
#include "lal/LALWavelet.h"

static int _f2l(int level, int indx);
static int _l2f(int level, int layer);
static void _getSliceL(int level, int layer, Wavelet *wavelet, Slice *slice);
static void _getSliceF(int indx, Wavelet *wavelet, Slice *slice);
static int _getOffset(int level,int layer);
static int _getMaxLevel(Wavelet *wavelet, UINT4 n);
static int _getMaxLayer(Wavelet *wavelet);
static UINT4 _limit(Slice s);
static int _getLayer(REAL4TimeSeries **layerOut, int indx, Wavelet *wavelet);
static void _putLayer(REAL4TimeSeries *layerData, int layer, Wavelet *wavelet);
static void _assignREAL4TimeSeriesMetadata(REAL4TimeSeries **left, REAL4TimeSeries *right);
static void _assignREAL4TimeSeries(REAL4TimeSeries **left, REAL4TimeSeries *right);
static void _assignWavelet(Wavelet **left,Wavelet *right);
static void  _assignClusterWavelet(ClusterWavelet **left, ClusterWavelet *right);

static void _createClusterWavelet(ClusterWavelet **w);
static double _setMask(ClusterWavelet *wavelet, int nc, BOOLEAN aura);
static void _pMaskPushBack(ClusterWavelet *wavelet, PixelWavelet *pix);
static void _copyPixel(PixelWavelet *to, PixelWavelet *from);
static void _assignPixel(PixelWavelet **to, PixelWavelet *from);
static BOOLEAN _create2DintArray(int ***a, int n1, int n2);
static void _destroy2DintArray(int ***a, int n1);
static void _print2DintArray(int ***a, FILE *out, int n1, int n2);
static int _clusterMain(ClusterWavelet *wavelet);
static int _clusterR(ClusterWavelet *wavelet, int k);

static BOOLEAN _allocateWavelet(Wavelet **wavelet);

static double _percentile(Wavelet *wavelet, double nonZeroFraction, 
			  BOOLEAN keepAmplitudes, REAL4 **median, REAL4 **norm50); 
static int _compare(const void *a, const void *b);
static void _clusterProperties(ClusterWavelet *w);

void _doubleToSecNan(double t, UINT4 *sec, UINT4 *nan);
double _secNanToDouble(UINT4 sec, UINT4 nan);
static int _nanoSeconds2steps(REAL4TimeSeries *ts, int nanoseconds);

static void _freeREAL4TimeSeries(REAL4TimeSeries **t);
static void _freeWavelet(Wavelet **w);
static void _freeClusterWavelet(ClusterWavelet **w);
static void _freeOutPercentile(OutputPercentileWavelet **p);
static void _freeOutCoincidence(OutputCoincidenceWavelet **co);
static void _freeOutCluster(OutputClusterWavelet **cl);

static void _setAmplitudes(ClusterWavelet *w);





static int _f2l(int level, int indx)
{
  int n = indx;
  int j,i;
  for(i=level-1; i>=1; i--) {
    j = ((1<<i) & (n));
    if(j) n = ((1<<i)-1) ^ (n);
  }
  return n;
}

static int _l2f(int level, int layer)
{
  int n = layer;
  int i;
  int j;
  for(i=1; i<level; i++) {
    j = (1<<i) & (n);
    if(j) n = ((1<<i)-1) ^ (n);
  }
  return n;
}

static void _getSliceL(int level, int layer, Wavelet *wavelet, Slice *slice)
{
  UINT4 waveletSize=wavelet->data->data->length;
  
  slice->start=_getOffset(level,layer);
  slice->step=(1<<level);
  slice->size=(waveletSize>>level);
  if(slice->start+(slice->size-1)*slice->step+1 > waveletSize)
    {
      fprintf(stderr,"Inconsistent slice: start=%d, step=%d, size=%d, waveletSize=%d\n",
	      slice->start, slice->step, slice->size, waveletSize);
      slice->size=-1;
    }
}

static void _getSliceF(int indx, Wavelet *wavelet, Slice *slice)
{

  int level=wavelet->level;
  int layer=abs(indx);
  
  int maxLayer=(wavelet->treeType==BINARY) ? ((1<<level) - 1) : level;

  if(layer>maxLayer){
    fprintf(stderr,"layer=%d > maxLayer=%d\n",layer,maxLayer);
    layer=maxLayer;
    fprintf(stderr,"Setting layer=maxLayer\n");
  }

  if(wavelet->treeType==BINARY){
    if(indx>0) layer=_f2l(level,layer);
  }
  else{
      if(layer) {               
	 level -= layer-1; 
	 layer = 1;
      }
      else{                     
	 layer = 0;
      }
  }
  _getSliceL(level,layer,wavelet,slice);
}

static int _getOffset(int level,int layer)
{
  int n=0,i;
  for(i=0; i<level; i++)
    if((layer>>i)&1) n += 1<<(level-1-i);
  return n;
}

static int _getMaxLevel(Wavelet *wavelet, UINT4 n)
{
  int maxLevel = 0;
  for(; (n>=2*wavelet->HPFilterLength) && (n>=2*wavelet->LPFilterLength) && !(n&1); n/=2) maxLevel++;
  return maxLevel;
} 

static int _getMaxLayer(Wavelet *wavelet)
{ 
  return ((wavelet->treeType==BINARY) ? (1<<(int)wavelet->level)-1 : (int)wavelet->level);
}

static UINT4 _limit(Slice s)
{ 
  return s.start + (s.size-1)*s.step + 1; 
}

static int _getLayer(REAL4TimeSeries **layerOut, int indx, Wavelet *wavelet)
{
  char name[LALNameLength];

  Slice s;
  UINT4 i;
  UINT4 waveletSize=wavelet->data->data->length;
  int maxLayer=_getMaxLayer(wavelet);

  sprintf(name,"layer indx=%d of wavelet ",indx);
  strcat(name,wavelet->data->name);

  if(indx>maxLayer) indx=maxLayer;
  _getSliceF(indx,wavelet,&s);

  if(_limit(s)<=waveletSize){
    _assignREAL4TimeSeriesMetadata(layerOut, wavelet->data);
    (*layerOut)->deltaT *= s.step;
    strcpy((*layerOut)->name,name);
    (*layerOut)->data->length=s.size;
    (*layerOut)->data->data=(REAL4*)LALMalloc(s.size*sizeof(REAL4));
    if(*layerOut==NULL)
      {
	fprintf(stderr,"Cannot allocate memory 13\n"); fflush(stderr);
	exit(1);
      }
    if((*layerOut)->data->data==NULL) return -1;
    for(i=0;i<s.size;i++){
      (*layerOut)->data->data[i]=wavelet->data->data->data[s.start+s.step*i];
    }
    return indx;
  }
  else{
    return -1;
  }
}

static void _putLayer(REAL4TimeSeries *layerData, int layer, Wavelet *wavelet)
{
  UINT4 i;
  Slice s;
  _getSliceF(layer,wavelet,&s);
  if(s.size < layerData->data->length || 
     _limit(s) > wavelet->data->data->length) {
    fprintf(stderr,"Invalid layer size\n");
  }
  else{
    for(i=0;i<s.size;i++){
      wavelet->data->data->data[s.start+i*s.step]=layerData->data->data[i];
    }
  }
} 

/* Allocates memory for left and _assignes all the metadata from right. 
   (*left)->name and (*left)->data->data are not _assigned.
   (*left)->data->data is not allocated.
*/

static void _assignREAL4TimeSeriesMetadata(REAL4TimeSeries **left, REAL4TimeSeries *right)
{
  *left=(REAL4TimeSeries*)LALMalloc(sizeof(REAL4TimeSeries));
  if(*left==NULL)
    {
      fprintf(stderr,"Cannot allocate memory 12\n");
      exit(1);
    }
  (*left)->data=(REAL4Sequence*)LALMalloc(sizeof(REAL4Sequence));
  if((*left)->data==NULL)
    {
      fprintf(stderr,"Cannot allocate memory 11\n");
      exit(1);
    }
  (*left)->epoch=right->epoch;
  (*left)->deltaT=right->deltaT;
  (*left)->f0=right->f0;
  (*left)->sampleUnits=right->sampleUnits;
  (*left)->data->length=right->data->length;
  (*left)->data->data=NULL;
  strcpy((*left)->name,"Assign a name manually");
}

/* REAL4TimeSeries copy constructor */
static void _assignREAL4TimeSeries(REAL4TimeSeries **left, REAL4TimeSeries *right)
{
  UINT4 i;
  _assignREAL4TimeSeriesMetadata(left,right);
  strcpy((*left)->name,right->name);
  (*left)->data->data=(REAL4*)LALCalloc(right->data->length,sizeof(REAL4));
  if((*left)->data->data==NULL)
    {
      fprintf(stderr,"Cannot allocate memory 10\n");
      exit(1);
    }
  for(i=0;i<right->data->length;i++){
    (*left)->data->data[i]=right->data->data[i];
  }
}

/* wavelet copy constructor */
static void _assignWavelet(Wavelet **left,Wavelet *right)
{

  *left=(Wavelet*)LALMalloc(sizeof(Wavelet));
  if(*left==NULL)
    {
      fprintf(stderr,"Cannot allocate memory 9\n");
      exit(1);
    }
  (*left)->type=right->type;
  (*left)->border=right->border;
  (*left)->treeType=right->treeType;
  (*left)->level=right->level;
  (*left)->HPFilterLength=right->HPFilterLength;
  (*left)->LPFilterLength=right->LPFilterLength;

  _assignREAL4TimeSeries(&((*left)->data),right->data);
}

static void  _assignClusterWavelet(ClusterWavelet **left, ClusterWavelet *right)
{
  UINT4 i,j;
  int M;
  _createClusterWavelet(left);
  if(right->wavelet!=NULL) _assignWavelet(&(*left)->wavelet,right->wavelet);
  if(right->original!=NULL) _assignWavelet(&(*left)->original,right->original);
  (*left)->pMaskCount=right->pMaskCount;
  (*left)->clusterCount=right->clusterCount;
  (*left)->clusterType=right->clusterType;
  (*left)->simulationType=right->simulationType;
  (*left)->nonZeroFractionAfterPercentile=right->nonZeroFractionAfterPercentile;
  (*left)->nonZeroFractionAfterCoincidence=right->nonZeroFractionAfterCoincidence;
  (*left)->nonZeroFractionAfterSetMask=right->nonZeroFractionAfterSetMask;
  (*left)->nonZeroFractionAfterClustering=right->nonZeroFractionAfterClustering;
  (*left)->nonZeroFractionAfterCuts=right->nonZeroFractionAfterCuts;
  (*left)->nonZeroFractionAfterVetoes=right->nonZeroFractionAfterVetoes;
  (*left)->pixelSwapApplied=right->pixelSwapApplied;
  (*left)->pixelMixerApplied=right->pixelMixerApplied;

  M = _getMaxLayer(right->wavelet)+1;
  if(right->medians!=NULL)
    {
      (*left)->medians=(REAL4*)LALCalloc(M,sizeof(REAL4));
      for(i=0;i<M;i++)
	{
	  (*left)->medians[i]=right->medians[i];
	}
    }
  if(right->norm50!=NULL)
    {
      (*left)->norm50=(REAL4*)LALCalloc(M,sizeof(REAL4));
      for(i=0;i<M;i++)
	{
	  (*left)->norm50[i]=right->norm50[i];
	}
    }
  if(right->pMask!=NULL)
    {
      (*left)->pMask=(PixelWavelet**)LALCalloc(right->pMaskCount,sizeof(PixelWavelet*));
      for(i=0;i<right->pMaskCount;i++)
	{
	  if(right->pMask[i]!=NULL)
	    {
	      _assignPixel(&(*left)->pMask[i],right->pMask[i]);
	    }
	}
    }
  if(right->sCuts!=NULL)
    {
      (*left)->sCuts=(UINT4*)LALCalloc(right->clusterCount,sizeof(UINT4));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->sCuts[i]=right->sCuts[i];
	}
    }

  if(right->volumes!=NULL)
    {
      (*left)->volumes=(UINT4*)LALCalloc(right->clusterCount,sizeof(UINT4));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->volumes[i]=right->volumes[i];
	}
    }

  if(right->cList!=NULL)
    {
      (*left)->cList=(UINT4**)LALCalloc(right->clusterCount,sizeof(UINT4*));
      for(i=0;i<right->clusterCount;i++)
	{
	  if(right->cList[i]!=NULL)
	    {
	      (*left)->cList[i]=(UINT4*)LALCalloc(right->volumes[i],sizeof(UINT4));
	      for(j=0;j<right->volumes[i];j++)
		{
		  (*left)->cList[i][j]=right->cList[i][j];
		}
	    }
	}
    }

  if(right->coreSize!=NULL)
    {
      (*left)->coreSize=(UINT4*)LALCalloc(right->clusterCount,sizeof(UINT4));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->coreSize[i]=right->coreSize[i];
	}
    }

  if(right->correlation!=NULL)
    {
      (*left)->correlation=(REAL4*)LALCalloc(right->clusterCount,sizeof(REAL4));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->correlation[i]=right->correlation[i];
	}
    }

   if(right->likelihood!=NULL)
    {
      (*left)->likelihood=(REAL4*)LALCalloc(right->clusterCount,sizeof(REAL4));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->likelihood[i]=right->likelihood[i];
	}
    }
 
   if(right->power!=NULL)
    {
      (*left)->power=(REAL4*)LALCalloc(right->clusterCount,sizeof(REAL4));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->power[i]=right->power[i];
	}
    }
 
   if(right->maxAmplitude!=NULL)
     {
      (*left)->maxAmplitude=(REAL4*)LALCalloc(right->clusterCount,sizeof(REAL4));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->maxAmplitude[i]=right->maxAmplitude[i];
	}
     }

   if(right->relativeStartTime!=NULL)
     {
      (*left)->relativeStartTime=(REAL8*)LALCalloc(right->clusterCount,sizeof(REAL8));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->relativeStartTime[i]=right->relativeStartTime[i];
	}
     }

   if(right->relativeStopTime!=NULL)
     {
      (*left)->relativeStopTime=(REAL8*)LALCalloc(right->clusterCount,sizeof(REAL8));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->relativeStopTime[i]=right->relativeStopTime[i];
	}
     }
  
   if(right->duration!=NULL)
    {
      (*left)->duration=(REAL8*)LALCalloc(right->clusterCount,sizeof(REAL8));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->duration[i]=right->duration[i];
	}
    }

   if(right->absoluteStartTime!=NULL)
     {
      (*left)->absoluteStartTime=
	(LIGOTimeGPS*)LALCalloc(right->clusterCount,sizeof(LIGOTimeGPS));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->absoluteStartTime[i]=right->absoluteStartTime[i];
	}
     }

   if(right->absoluteStopTime!=NULL)
     {
      (*left)->absoluteStopTime=
	(LIGOTimeGPS*)LALCalloc(right->clusterCount,sizeof(LIGOTimeGPS));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->absoluteStopTime[i]=right->absoluteStopTime[i];
	}
     }

   if(right->startFrequency!=NULL)
    {
      (*left)->startFrequency=
	(REAL4*)LALCalloc(right->clusterCount,sizeof(REAL4));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->startFrequency[i]=right->startFrequency[i];
	}
    }

   if(right->stopFrequency!=NULL)
    {
      (*left)->stopFrequency=
	(REAL4*)LALCalloc(right->clusterCount,sizeof(REAL4));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->stopFrequency[i]=right->stopFrequency[i];
	}
    }

   if(right->bandwidth!=NULL)
    {
      (*left)->bandwidth=(REAL4*)LALCalloc(right->clusterCount,sizeof(REAL4));
      for(i=0;i<right->clusterCount;i++)
	{
	  (*left)->bandwidth[i]=right->bandwidth[i];
	}
    }


}

static int _nanoSeconds2steps(REAL4TimeSeries *ts, int nanoseconds)
{
  return nanoseconds*pow(10,-9)/ts->deltaT+0.5;
}

static double _setMask(ClusterWavelet *w, int nc, BOOLEAN aura)
{
  register int i;
  register int j;
  int x, y, L, k;
  register int* q = NULL;
  register int* p = NULL; 
  int* pp; 
  int* pm;
  PixelWavelet pix;
  REAL4TimeSeries *a;
  int ni;
  int nj;
  int n;
  int m;
  int nPixel;
  int **FT;
  int **XY;
  BOOLEAN status;
  int maxPixels;
  int *f;
  int *t;
  register int* pN;
  int nM;
  Slice S;
  Wavelet *original=w->original;
  
  if(w->wavelet->treeType!=BINARY) return 1.;

  ni = _getMaxLayer(w->wavelet)+1;
  nj = w->wavelet->data->data->length/ni;
  n = ni-1;
  m = nj-1;
  nPixel = 0;
  
  status=_create2DintArray(&FT,ni,nj);

  if(!status){
    printf("Memory allocation problem\n");
    exit(1);
  }
  status=_create2DintArray(&XY,ni,nj);
  if(!status){
    printf("Memory allocation problem\n");
    exit(1);
  }

  for(i=0; i<ni; i++){
    p  = FT[i]; 
    q  = XY[i];
    _getLayer(&a, i, w->wavelet);

    for(j=0; j<nj; j++){
      p[j] = (a->data->data[j]!=0) ? 1 : 0;
      if(p[j]) nPixel++;
      q[j] = 0;
    }
    _freeREAL4TimeSeries(&a);
  }


  if(!nc || ni<3 || nPixel<2) 
    return (double)(nPixel)/((double)(w->wavelet->data->data->length));


  if(FT[0][0]) XY[0][0] = FT[0][1]   + FT[1][0]   + FT[1][1];
  if(FT[0][m]) XY[0][m] = FT[1][m]   + FT[1][m-1] + FT[0][m-1];
  if(FT[n][0]) XY[n][0] = FT[n-1][0] + FT[n-1][1] + FT[n][1];
  if(FT[n][m]) XY[n][m] = FT[n][m-1] + FT[n-1][m] + FT[n-1][m-1];
  

  for(j=1; j<m; j++){
    p  = FT[0]; 
    q  = FT[n];

    if(p[j]){
      pp = FT[1];
      XY[0][j] = p[j-1]+p[j+1] + pp[j-1]+pp[j]+pp[j+1];
    }
    if(q[j]){
      pm = FT[n-1];
      XY[n][j] = q[j-1]+q[j+1] + pm[j-1]+pm[j]+pm[j+1];
    }
  }

  for(i=1; i<n; i++){
    pm = FT[i-1]; p  = FT[i]; pp = FT[i+1]; q = XY[i];

    if(p[0]) 
      q[0] = p[1] + pm[0]+pm[1] + pp[0]+pp[1];

    if(p[m]) 
      q[m] = p[m-1] + pm[m]+pm[m-1] + pp[m]+pp[m-1];

    for(j=1; j<m; j++){
      if(p[j])
	q[j] = pm[j-1]+pm[j]+pm[j+1] + pp[j-1]+pp[j]+pp[j+1] + p[j-1]+p[j+1];
    }
  }

  /**************************************/
  /* remove clusters with 2,3 pixels  */
  /**************************************/

  if(nc>1){


    if(XY[0][0]){ 
      x = XY[0][1]   + XY[1][0]   + XY[1][1];
      if(x==1 || x==4) XY[0][0]=XY[0][1]=XY[1][0]=XY[1][1]=0;
    }
    if(XY[0][m]){ 
      x = XY[1][m]   + XY[1][m-1] + XY[0][m-1];
      if(x==1 || x==4) XY[0][m]=XY[1][m]=XY[1][m-1]=XY[0][m-1]=0;
    }
    if(XY[n][0]){ 
      x = XY[n-1][0] + XY[n-1][1] + XY[n][1];
      if(x==1 || x==4) XY[n][0]=XY[n-1][0]=XY[n-1][1]=XY[n][1]=0;
    }
    if(XY[n][m]){
      x = XY[n-1][m]   + XY[n][m-1] + XY[n-1][m-1];
      if(x==1 || x==4) XY[n][m]=XY[n-1][m]=XY[n][m-1]=XY[n-1][m-1]=0;
    }


    for(j=1; j<m; j++){
      p  = XY[0]; 
      q  = XY[n];

      if(p[j]==1 || p[j]==2){
	if(p[j-1]+p[j+1] < 4){
	  pp = XY[1];
	  L = p[j-1]+p[j+1] + pp[j];
	  x = pp[j-1] + pp[j+1] + L;
	       
	  if(x==1 || (p[j]==2 && nc>2 && (x==2 || (x==4 && L==4))))
	    p[j]=p[j-1]=p[j+1]=pp[j-1]=pp[j]=pp[j+1]=0;
	}
      }

      if(q[j]==1 || q[j]==2){
	if(q[j-1]+q[j+1] < 4){
	  pm = XY[n-1];
	  L = q[j-1]+q[j+1] + pm[j];
	  x = pm[j-1] + pm[j+1] + L;

	  if(x==1 || (q[j]==2 && nc>2 && (x==2 || (x==4 && L==4))))
	    q[j]=q[j-1]=q[j+1]=pm[j-1]=pm[j]=pm[j+1]=0;
	}
      }
    }


    for(i=1; i<n; i++){
      pm = XY[i-1];
      p  = XY[i];
      pp = XY[i+1];
	 
	 
      if(p[0]==1 || p[0]==2){ 
	if(pm[0]+pp[0] < 4){
	  L = p[1] + pm[0] + pp[0];
	  x = pm[1] + pp[1] + L;
	       
	  if(x==1 || (p[0]==2 && nc>2 && (x==2 || (x==4 && L==4))))
	    p[0]=pm[0]=pp[0]=pm[1]=pp[1]=p[1]=0;;
	}
      }

      if(p[m]==1 || p[m]==2){ 
	if(pm[m]+pp[m] < 4){
	  L = p[m-1] + pm[m] + pp[m];
	  x = pm[m-1] + pp[m-1] + L;
	       
	  if(x==1 || (p[m]==2 && nc>2 && (x==2 || (x==4 && L==4))))
	    p[m]=pm[m]=pp[m]=pm[m-1]=pp[m-1]=p[m-1]=0;
	}
      }
	 
      for(j=1; j<m; j++){
	y = p[j];
	if(y == 1 || y == 2){
	  if(pm[j]+pp[j] >3) continue;
	  if(p[j-1]+p[j+1] >3) continue;

	  L  = pm[j]+pp[j] + p[j-1]+p[j+1];
	  x  = pm[j-1]+pm[j+1] + pp[j-1]+pp[j+1] + L;
	       
	  if(x==1 || (y==2 && nc>2 && (x==2 || (x==4 && L==4))))
	    p[j]=p[j-1]=p[j+1]=pm[j-1]=pm[j]=pm[j+1]=pp[j-1]=pp[j]=pp[j+1]=0;
	}
      }
    }
  }

  L = 0;
  maxPixels=(int)(w->wavelet->data->data->length*w->nonZeroFractionAfterCoincidence);

/*    printf("nonzerofraction after coincidence: %f\n",w->nonZeroFractionAfterCoincidence); fflush(stdout); */

  w->pMask=(PixelWavelet**)LALCalloc(maxPixels,sizeof(PixelWavelet*));
  if(w->pMask==NULL)
    {
      fprintf(stderr,"Cannot allocate memory 1\n");
      exit(1);
    }
  w->pMaskCount=0;

  w->sCuts=NULL;
  w->cList=NULL;

  pix.amplitude=-1.0;
  pix.amplitudeOriginal=-1.0;
  pix.clusterID = 0;     
  pix.neighborsCount=0;

  f = &pix.frequency;
  t = &pix.time;

  for(i=0; i<ni; i++){
    p  = FT[i]; q  = XY[i];
    for(j=0; j<nj; j++)
      p[j] = q[j];
  }

  for(i=0; i<ni; i++){
    q  = XY[i];
    p  = FT[i];
    for(j=0; j<nj; j++){
      if(q[j]) {
	*t = j; *f = i;              
	pix.core = TRUE;             
	_pMaskPushBack(w, &pix);    
	p[j] = ++L;                 


	if(aura){                 
	  pix.core = FALSE;       

	  if(i>0 && j>0) { 
	    *t=j-1; *f=i-1; 
	    if(!FT[*f][*t]) {
	      _pMaskPushBack(w, &pix);  		    
	      FT[*f][*t] = ++L;
	    }
	  }
	  if(i>0) {
	    *t=j;   *f=i-1; 
	    if(!FT[*f][*t]) {
	      _pMaskPushBack(w, &pix);   		    
	      FT[*f][*t] = ++L;
	    }
	  }
	  if(i>0 && j<m) { 
	    *t=j+1; *f=i-1; 
	    if(!FT[*f][*t]) {
	      _pMaskPushBack(w, &pix);    		    
	      FT[*f][*t] = ++L;
	    }
	  }
	  if(j>0)        { 
	    *t=j-1; *f=i;   
	    if(!FT[*f][*t]) {
	      _pMaskPushBack(w, &pix);    		    
	      FT[*f][*t] = ++L;
	    }
	  }
	  if(j<m)        { 
	    *t=j+1; *f=i;   
	    if(!FT[*f][*t]) {
	      _pMaskPushBack(w, &pix);   		    
	      FT[*f][*t] = ++L;
	    }
	  }
	  if(i<n && j>0) { 
	    *t=j-1; *f=i+1; 
	    if(!FT[*f][*t]) {
	      _pMaskPushBack(w, &pix);    		    
	      FT[*f][*t] = ++L;
	    }
	  }
	  if(i<n)        { 
	    *t=j;   *f=i+1; 
	    if(!FT[*f][*t]) {
	      _pMaskPushBack(w, &pix);    		    
	      FT[*f][*t] = ++L;
	    }
	  }
	  if(i<n && j<m) { 
	    *t=j+1; *f=i+1; 
	    if(!FT[*f][*t]) {
	      _pMaskPushBack(w, &pix);   		    
	      FT[*f][*t] = ++L;
	    }
	  }
	}             
      }
    }
  }

  nM = w->pMaskCount;


  for(k=0; k<nM; k++){

    i = w->pMask[k]->frequency;
    j = w->pMask[k]->time;

    _getSliceF(i,w->wavelet,&S); 
    w->pMask[k]->amplitude=-1.0;
    w->pMask[k]->amplitudeOriginal=-1.0;
    w->pMask[k]->amplitude=w->wavelet->data->data->data[S.start+S.step*j];
    w->pMask[k]->amplitudeOriginal=original->data->data->data[S.start+S.step*j];

    w->pMask[k]->neighborsCount=0;

    pN = &(w->pMask[k]->neighbors[0]);
    L = 0;

    if(i==0 || i==n){                        
      if(i==0){ p = FT[0]; q = FT[1];}
      if(i==n){ p = FT[n]; q = FT[n-1];}

      if(j==0){                             
	if(p[1]) pN[L++] = p[1];
	if(q[1]) pN[L++] = q[1];
	if(q[0]) pN[L++] = q[0];
      }
      else if(j==m){                       
	if(p[m-1]) pN[L++] = p[m-1];
	if(q[m-1]) pN[L++] = q[m-1];
	if(q[m]<0)   pN[L++] = q[m];
      }
      else{                                
	if(p[j-1]) pN[L++] = p[j-1];
	if(p[j+1]) pN[L++] = p[j+1];
	if(q[j-1]) pN[L++] = q[j-1];
	if(q[j])   pN[L++] = q[j];
	if(q[j+1]) pN[L++] = q[j+1];
      }
    }

    else{          
      pp = FT[i+1];
      p  = FT[i];
      pm = FT[i-1];

      if(j==0){                              
	if(pm[0]) pN[L++] = pm[0];
	if(pp[0]) pN[L++] = pp[0];
	if( p[1]) pN[L++] =  p[1];
	if(pm[1]) pN[L++] = pm[1];
	if(pp[1]) pN[L++] = pp[1];
      }
      else if(j==m){
	if(pm[m])   pN[L++] = pm[m];        
	if(pp[m])   pN[L++] = pp[m];
	if( p[m-1]) pN[L++] =  p[m-1];
	if(pm[m-1]) pN[L++] = pm[m-1];
	if(pp[m-1]) pN[L++] = pp[m-1];
      }
      else{
	if(pm[j-1]) pN[L++] = pm[j-1];
	if(pm[j])   pN[L++] = pm[j];
	if(pm[j+1]) pN[L++] = pm[j+1];
	if( p[j-1]) pN[L++] =  p[j-1];
	if( p[j+1]) pN[L++] =  p[j+1];
	if(pp[j-1]) pN[L++] = pp[j-1];
	if(pp[j])   pN[L++] = pp[j];
	if(pp[j+1]) pN[L++] = pp[j+1];
      }
    }

    w->pMask[k]->neighborsCount=L;

    if(!aura){
      x = XY[i][j];
      
      if(x != L || L>8){ 
	printf("pMask size error: L=%d, x=%d, k=%d, i=%d, j=%d\n",L,x,k,i,j);
      }
    }
  }
/*    printf("-----------------------"); */
/*    printf("FT\n"); */
/*    _print2DintArray(&FT, stdout,ni,nj); */
/*    printf("-----------------------"); */
/*    printf("-----------------------"); */
/*    printf("XY\n"); */
/*    _print2DintArray(&XY, stdout,ni,nj); */
/*    printf("-----------------------"); */

  _destroy2DintArray(&FT,ni);
  _destroy2DintArray(&XY,ni);

/*    printf("pMaskCount=%d maxPixels=%d\n",w->pMaskCount,maxPixels); fflush(stdout); */

  w->pMask=(PixelWavelet**)LALRealloc(w->pMask,w->pMaskCount*sizeof(PixelWavelet*));

  return (double)(w->pMaskCount+1)/(double)(w->wavelet->data->data->length);
}


static void _pMaskPushBack(ClusterWavelet *w, PixelWavelet *pix)
{
  PixelWavelet *newPix=(PixelWavelet*)LALMalloc(sizeof(PixelWavelet));
  if(newPix==NULL)
    {
      fprintf(stderr,"Cannot allocate memory 2\n");
      exit(1);
    }
  _copyPixel(newPix,pix);
  w->pMask[++w->pMaskCount-1]=newPix;
}

static void _copyPixel(PixelWavelet *to, PixelWavelet *from)
{
  int i;
  to->time=from->time;
  to->frequency=from->frequency;
  to->clusterID=from->clusterID;
  to->core=from->core;
  to->amplitude=from->amplitude;
  to->amplitudeOriginal=from->amplitudeOriginal;
  to->neighborsCount=from->neighborsCount;
  for(i=0;i<8;i++){
    to->neighbors[i]=from->neighbors[i];
  }
}

static void _assignPixel(PixelWavelet **to, PixelWavelet *from)
{
  if(*to!=NULL) LALFree(*to);
  *to=(PixelWavelet*)LALCalloc(1,sizeof(PixelWavelet));
  _copyPixel(*to,from);
}

static BOOLEAN _create2DintArray(int ***a, int n1, int n2)
{
  int i;
  BOOLEAN result=TRUE;

  *a=(int**)LALMalloc(n1*sizeof(int*));


  if(*a==NULL){ 
    fprintf(stderr,"Cannot allocate memory 3\n");
    exit(1);
  }

  for(i=0;i<n1;i++){
    (*a)[i]=(int*)LALMalloc(n2*sizeof(int));
    if((*a)[i]==NULL)
      {
	fprintf(stderr,"Cannot allocate memory 4\n");
	exit(1);
      }
  }

  return result;
}

static void _print2DintArray(int ***a, FILE *out, int n1, int n2)
{
  int i,j;
  for(i=0;i<n1;i++)
    {
      for(j=0;j<n2;j++)
	{
	  fprintf(out,"%d ",(*a)[i][j]);
	}
      fprintf(out,"\n");
    }
}



static void _destroy2DintArray(int ***a, int n1)
{
  int i;
  for(i=0;i<n1;i++){
    if((*a)[i]!=NULL) LALFree((*a)[i]);
    (*a)[i]=NULL;
  }
  if((*a)!=NULL) LALFree((*a));
  (*a)=NULL;
}

static int _clusterMain(ClusterWavelet *w)
{
  int volume;
  size_t i,m;
  size_t ncluster = 0;
  size_t n = w->pMaskCount;
  UINT4 k;


  w->cList=(UINT4**)LALCalloc(n,sizeof(UINT4*));
  if(w->cList==NULL)
    {
      fprintf(stderr,"Cannot allocate memory 5\n");
      exit(1);
    }
  w->sCuts=(UINT4*)LALCalloc(n,sizeof(UINT4));
  if(w->sCuts==NULL)
    {
      fprintf(stderr,"Cannot allocate memory 6\n");
      exit(1);
    }
  w->volumes=(UINT4*)LALCalloc(n,sizeof(UINT4));
  if(w->volumes==NULL)
    {
      fprintf(stderr,"Cannot allocate memory 7\n");
      exit(1);
    }

  m=0;
  for(i=0; i<n; i++){
    if(w!=NULL && w->pMask!=NULL && w->pMask[i]!=NULL && !w->pMask[i]->clusterID){
      w->pMask[i]->clusterID = ++ncluster;
/*        printf("HERE!!!\n");fflush(stdout); */
      volume = _clusterR(w,i);
/*        printf("volume=%d\n",volume);fflush(stdout); */
      w->cList[m]=(UINT4*)LALCalloc(volume,sizeof(UINT4));
      if(w->cList[m]==NULL)
	{
	  fprintf(stderr,"Cannot allocate memory 8\n");
	  exit(1);
	}
      w->volumes[m]=volume;
      w->sCuts[m]=0;
      m++;
/*        printf("Here: volume=%d i=%d n=%d m=%d\n",volume,i,n,m); fflush(stdout); */
    }
  }

/*    printf("Ever get here?\n");fflush(stdout); */
  
  w->clusterCount=ncluster;

  for(k=0;k<ncluster;k++){
    m=0;
    for(i=0;i<n;i++){
      if(w->pMask[i]->clusterID==k+1)
	{
	  w->cList[k][m++]=i;
	}
    }

    if(w->volumes[k]!=m){
      fprintf(stderr,"Size mismatch: m=%d, volume=%d\n",m,w->volumes[k]);
    }
  }

  return ncluster;
}

static int _clusterR(ClusterWavelet *w, int k)
{
/*    size_t i,j; */
/*    size_t n = w->pMask[k]->neighborsCount; */

  UINT4 i,j;
  int volume = 1;
  int ncluster;
  UINT4 n;
  UINT4 *p;


/*    if(w!=NULL && w->pMask!=NULL && w->pMask[k]!=NULL) */
/*      { */
      ncluster = w->pMask[k]->clusterID;
      n = w->pMask[k]->neighborsCount;
      p = n ? &(w->pMask[k]->neighbors[0]) : NULL;
/*        printf("In if: ncluster=%d, neighborsCount=%d k=%d\n",ncluster,n,k);fflush(stdout); */
/*      } */
/*    else  */
/*      return volume; */

/*    printf("Are we here?\n");fflush(stdout); */

  for(i=0; i<n; i++){
    j = p[i]-1;
/*      printf("j=%d pMaskCount=%d\n",j,w->pMaskCount);fflush(stdout); */
    if(w!=NULL && w->pMask!=NULL && w->pMask[j]!=NULL && !w->pMask[j]->clusterID){
      w->pMask[j]->clusterID = ncluster;
/*        printf("ncluster=%d i=%d j=%d\n",ncluster,i,j);fflush(stdout); */
/*        printf("Crash here begin \n");fflush(stdout); */
      volume += _clusterR(w,j);
/*        printf("Crash here end\n"); fflush(stdout); */
/*        printf("inside clusterR volume=%d\n", volume);fflush(stdout); */
    }
  }
  return volume;
}

static void _freeREAL4TimeSeries(REAL4TimeSeries **t)
{
  LALFree((*t)->data->data);
  LALFree((*t)->data);
  LALFree((*t));
}

static void _freeWavelet(Wavelet **w)
{
  _freeREAL4TimeSeries(&(*w)->data);
  LALFree((*w));
}


/*  static void _freeWavelet_TS_PM(Wavelet **w) */
/*  { */
/*    UINT4 i; */
/*    for(i=0;i<(*w)->pMaskCount;i++) */
/*      { */
/*        LALFree((*w)->pMask[i]); */
/*      } */
/*    LALFree((*w)->pMask); */
/*    _freeWavelet(w); */
/*  } */

/*  static void _freeWavelet_ALL(Wavelet **w) */
/*  { */
/*    UINT4 i; */

/*    for(i=0;i<(*w)->clusterCount;i++) */
/*      { */
/*        LALFree((*w)->cList[i]); */
/*      } */
/*    LALFree((*w)->cList); */
/*    LALFree((*w)->volumes); */
/*    LALFree((*w)->sCuts); */
/*    _freeWavelet_TS_PM(w); */
/*  } */

static BOOLEAN _allocateWavelet(Wavelet **wavelet)
{
  *wavelet=(Wavelet*)LALMalloc(sizeof(Wavelet));
  if(*wavelet==NULL) return FALSE;
  (*wavelet)->data=(REAL4TimeSeries*)LALMalloc(sizeof(REAL4TimeSeries));
  if((*wavelet)->data==NULL) return FALSE;
  (*wavelet)->data->data=(REAL4Sequence*)LALMalloc(sizeof(REAL4Sequence));
  if((*wavelet)->data->data==NULL) return FALSE;
  return TRUE;
}

static void _freeOutPercentile(OutputPercentileWavelet **p)
{
  _freeClusterWavelet(&(*p)->out);
  LALFree(*p);
  *p=NULL;
}

static void _freeOutCoincidence(OutputCoincidenceWavelet **co)
{
  _freeClusterWavelet(&(*co)->one);
  _freeClusterWavelet(&(*co)->two);
  LALFree(*co);
  *co=NULL;
}

static void _freeOutCluster(OutputClusterWavelet **cl)
{
  _freeClusterWavelet(&(*cl)->w);
  LALFree(*cl);
  *cl=NULL;
}

static double _percentile(Wavelet *wavelet, double nonZeroFraction, 
			  BOOLEAN keepAmplitudes, REAL4 **median, REAL4 **norm50)
{

  INT4 i, j, M, nS, boundary1, boundary2;
  REAL4TimeSeries *a;
  REAL4 **aPtr;

  M = _getMaxLayer(wavelet)+1;
  nS=wavelet->data->data->length/M;
  aPtr=(REAL4**)LALMalloc((sizeof(REAL4*)*nS));

  *median=(REAL4*)LALCalloc(M,sizeof(REAL4));
  *norm50=(REAL4*)LALCalloc(M,sizeof(REAL4));

  for(i=0; i<M; i++){
    _getLayer(&a,i,wavelet);

    for(j=0;j<nS;j++){
      aPtr[j]=a->data->data+j;
    }

    qsort(aPtr,nS,sizeof(REAL4*),_compare);

    (*norm50)[i]=
      -((*aPtr[nS/4]) + (*aPtr[nS/4+1]))/2+
      ((*aPtr[(3*nS)/4]) + (*aPtr[(3*nS)/4-1]))/2;
    (*norm50)[i]/=2;

/*     printf("In percentile: i=%d norm=%f\n",i,(*norm50)[i]);fflush(stdout); */

    if(nS%2==0){
      (*median)[i]=(*aPtr[nS/2-1] + *aPtr[nS/2])/2;
      boundary1=nonZeroFraction*nS/2;
      for(j=0;j<boundary1;j++){
	if(!keepAmplitudes) *aPtr[j]=-((REAL4)nS)/2.0/(j+1.0);
      }
      boundary2=nS-boundary1;
      for(j=boundary1;j<boundary2;j++){
	*aPtr[j]=0.0;
      }
      for(j=boundary2;j<nS;j++){
	if(!keepAmplitudes) *aPtr[j]=((REAL4)nS)/2.0/(nS-j);
      }
    }
    else{
      (*median)[i]=*aPtr[nS/2];
      /* so far no difference below */
      boundary1=(int)(nonZeroFraction*nS/2);
      for(j=0;j<boundary1;j++){
	if(!keepAmplitudes) *aPtr[j]=-((REAL4)nS)/2.0/(j+1);
      }
      boundary2=nS/2+boundary1;
      for(j=boundary1;j<=boundary2;j++){
	*aPtr[j]=0.0;
      }
      for(j=boundary2+1;j<nS;j++){
	if(!keepAmplitudes) *aPtr[j]=((REAL4)nS)/2.0/(nS-j);
      }
    }
    _putLayer(a, i, wavelet);
    _freeREAL4TimeSeries(&a);
  }

  _getLayer(&a,0,wavelet);
  bzero(a->data->data, a->data->length*sizeof(REAL4));
  _putLayer(a,0,wavelet);
  LALFree(aPtr);

  return nonZeroFraction;
}

static int _compare(const void *a, const void *b)
{
  REAL4 **A=(REAL4 **)a;
  REAL4 **B=(REAL4 **)b;
  if(**A < **B) return -1;
  if(**A == **B) return 0;
  if(**A > **B) return 1;
  return -2;
}

static void _clusterProperties(ClusterWavelet *w)
{
  UINT4 i, j, f, ff, t;
  double a,b;
  double delta_t, delta_f;
  double x;
  int N;
  Slice s;

  delta_t=w->wavelet->data->deltaT*(1<<w->wavelet->level);
  delta_f=1/w->wavelet->data->deltaT/(1<<w->wavelet->level)/2.;
  N=w->wavelet->data->data->length/(1<<w->wavelet->level);
  w->delta_t=delta_t;
  w->delta_f=delta_f;

  if(w->coreSize!=NULL) LALFree(w->coreSize);
  if(w->correlation!=NULL) LALFree(w->correlation);
  if(w->likelihood!=NULL) LALFree(w->likelihood);
  if(w->power!=NULL) LALFree(w->power);
  if(w->maxAmplitude!=NULL) LALFree(w->maxAmplitude);
  if(w->relativeStartTime!=NULL) LALFree(w->relativeStartTime);
  if(w->relativeStopTime!=NULL) LALFree(w->relativeStopTime);
  if(w->duration!=NULL)  LALFree(w->duration);
  if(w->absoluteStartTime!=NULL) LALFree(w->absoluteStartTime);
  if(w->absoluteStopTime!=NULL) LALFree(w->absoluteStopTime);
  if(w->startFrequency!=NULL) LALFree(w->startFrequency);
  if(w->stopFrequency!=NULL) LALFree(w->stopFrequency);
  if(w->bandwidth!=NULL) LALFree(w->bandwidth);
  if(w->blobs!=NULL)
    {
      for(i=0;i<w->clusterCount;i++)
	{
	  if(w->blobs[i].pBlob!=NULL) 
	    {
	      LALFree(w->blobs[i].pBlob);
	      w->blobs[i].pBlob=NULL;
	    }
	  if(w->blobs[i].oBlob!=NULL) 
	    {
	      LALFree(w->blobs[i].oBlob);
	      w->blobs[i].oBlob=NULL;
	    }
	}
      LALFree(w->blobs);
      w->blobs=NULL;
    }

  w->coreSize=(UINT4*)LALCalloc(w->clusterCount,sizeof(UINT4));
  w->correlation=(REAL4*)LALCalloc(w->clusterCount,sizeof(REAL4));
  w->likelihood=(REAL4*)LALCalloc(w->clusterCount,sizeof(REAL4));
  w->power=(REAL4*)LALCalloc(w->clusterCount,sizeof(REAL4));
  w->maxAmplitude=(REAL4*)LALCalloc(w->clusterCount,sizeof(REAL4));
  w->relativeStartTime=(REAL8*)LALCalloc(w->clusterCount,sizeof(REAL8));
  w->relativeStopTime=(REAL8*)LALCalloc(w->clusterCount,sizeof(REAL8));
  w->absoluteStartTime=(LIGOTimeGPS*)LALCalloc(w->clusterCount,sizeof(LIGOTimeGPS));
  w->absoluteStopTime=(LIGOTimeGPS*)LALCalloc(w->clusterCount,sizeof(LIGOTimeGPS));
  w->duration=(REAL8*)LALCalloc(w->clusterCount,sizeof(REAL8));
  w->startFrequency=(REAL4*)LALCalloc(w->clusterCount,sizeof(REAL4));
  w->stopFrequency=(REAL4*)LALCalloc(w->clusterCount,sizeof(REAL4));
  w->bandwidth=(REAL4*)LALCalloc(w->clusterCount,sizeof(REAL4));
  w->blobs=(ClusterBlobWavelet*)LALCalloc(w->clusterCount,sizeof(ClusterBlobWavelet));

  for(i=0;i<w->clusterCount;i++)
    {
      if(w->sCuts[i]) continue;

      w->coreSize[i]=0;
      w->correlation[i]=0.0;
      w->likelihood[i]=0.0;
      w->power[i]=0.0;
      w->maxAmplitude[i]=0.0;
      w->relativeStartTime[i]=N*delta_t;
      w->relativeStopTime[i]=0.0;
      w->startFrequency[i]=((REAL4)w->wavelet->data->data->length)/((REAL4)2.0);
      w->stopFrequency[i]=0.0;

      for(j=0;j<w->volumes[i];j++)
	{
	  if(w->pMask[w->cList[i][j]]->core) 
	    {
	      w->coreSize[i]++;
	      a=w->pMask[w->cList[i][j]]->amplitude;
	      b=w->pMask[w->cList[i][j]]->amplitudeOriginal;
	      b-=w->medians[w->pMask[w->cList[i][j]]->frequency];
	      b/=w->norm50[w->pMask[w->cList[i][j]]->frequency];
	      if(a>0) w->correlation[i]++;
	      if(a<0) w->correlation[i]--;
	      w->likelihood[i]+=log(fabs(a));
	      w->power[i]+=b*b;
	      if(w->maxAmplitude[i]<fabs(a)) w->maxAmplitude[i]=fabs(a);

	      f=w->pMask[w->cList[i][j]]->frequency;
	      _getSliceF(f,w->wavelet,&s);
	      ff=s.size/N;

	      t=w->pMask[w->cList[i][j]]->time/((double)ff);
	      x=delta_t*t;
	      if(x<w->relativeStartTime[i]) 
		{
		  w->relativeStartTime[i]=x;
		  w->blobs[i].start_time_indx=t;
		}
	      if(x>w->relativeStopTime[i])
		{
		  w->relativeStopTime[i]=x;
		  w->blobs[i].stop_time_indx=t;
		}

	      x=delta_f*f*ff;
	      if(x<w->startFrequency[i]) 
		{
		  w->startFrequency[i]=x;
		  w->blobs[i].start_freq_indx=f*ff;
		}
	      if(x>w->stopFrequency[i]) 
		{
		  w->stopFrequency[i]=x;
		  w->blobs[i].stop_freq_indx=f*ff;
		}
	    }
	}

      w->blobs[i].time_width = w->blobs[i].stop_time_indx - 
	w->blobs[i].start_time_indx + 1;

      w->blobs[i].freq_width = w->blobs[i].stop_freq_indx -
	w->blobs[i].start_freq_indx + 1;

      w->blobs[i].pBlob = 
	(REAL4*)LALCalloc(w->blobs[i].time_width*w->blobs[i].freq_width,
			  sizeof(REAL4));
      w->blobs[i].oBlob =
	(REAL4*)LALCalloc(w->blobs[i].time_width*w->blobs[i].freq_width,
			  sizeof(REAL4));

      for(j=0;j<w->volumes[i];j++)
	{
	  if(w->pMask[w->cList[i][j]]->core) 
	    {
	      w->blobs[i].pBlob[(w->pMask[w->cList[i][j]]->frequency - w->blobs[i].start_freq_indx) * 
			       w->blobs[i].time_width + (w->pMask[w->cList[i][j]]->time - 
							 w->blobs[i].start_time_indx )] = 
		w->pMask[w->cList[i][j]]->amplitude;
	      w->blobs[i].oBlob[(w->pMask[w->cList[i][j]]->frequency - w->blobs[i].start_freq_indx) * 
			       w->blobs[i].time_width + (w->pMask[w->cList[i][j]]->time - 
							 w->blobs[i].start_time_indx )] = 
		w->pMask[w->cList[i][j]]->amplitudeOriginal;
	    }
	}

      /*      w->likelihood[i]/=w->coreSize[i];*/
      w->correlation[i]/=w->coreSize[i];
      w->power[i]/=w->coreSize[i];

      x=_secNanToDouble(w->wavelet->data->epoch.gpsSeconds, 
			w->wavelet->data->epoch.gpsNanoSeconds) + 
	w->relativeStartTime[i];
      _doubleToSecNan(x, &w->absoluteStartTime[i].gpsSeconds, 
		      &w->absoluteStartTime[i].gpsNanoSeconds);

      x=_secNanToDouble(w->wavelet->data->epoch.gpsSeconds, 
			w->wavelet->data->epoch.gpsNanoSeconds) + 
	w->relativeStopTime[i];
      _doubleToSecNan(x, &w->absoluteStopTime[i].gpsSeconds, 
		      &w->absoluteStopTime[i].gpsNanoSeconds);

      w->duration[i]=w->relativeStopTime[i] - w->relativeStartTime[i];
      w->bandwidth[i]=w->stopFrequency[i] - w->startFrequency[i] + delta_f;
    }
}

void _doubleToSecNan(double t, UINT4 *sec, UINT4 *nan)
{
  *sec=(UINT4)t;
  *nan=(UINT4)((t-*sec)*pow(10,9));
}

double _secNanToDouble(UINT4 sec, UINT4 nan)
{
  return ((double)sec + (double)nan/pow(10,9));
}

static void _createClusterWavelet(ClusterWavelet **w)
{
  *w=(ClusterWavelet*)LALMalloc(sizeof(ClusterWavelet));
  (*w)->wavelet=NULL;
  (*w)->original=NULL;
  (*w)->medians=NULL;
  (*w)->norm50=NULL;
  (*w)->pMaskCount=0;
  (*w)->clusterCount=0;
  (*w)->clusterCountFinal=0;
  (*w)->clusterType=ORIGINAL_CL;
  (*w)->simulationType=0;
  (*w)->pMask=NULL;
  (*w)->sCuts=NULL;  
  (*w)->cList=NULL;
  (*w)->volumes=NULL;
  (*w)->coreSize=NULL;
  (*w)->correlation=NULL;
  (*w)->likelihood=NULL;
  (*w)->power=NULL;
  (*w)->maxAmplitude=NULL;
  (*w)->relativeStartTime=NULL;
  (*w)->relativeStopTime=NULL;
  (*w)->duration=NULL;
  (*w)->absoluteStartTime=NULL;
  (*w)->absoluteStopTime=NULL;
  (*w)->startFrequency=NULL;
  (*w)->stopFrequency=NULL;
  (*w)->bandwidth=NULL;
  (*w)->nonZeroFractionAfterPercentile=-1.0;
  (*w)->nonZeroFractionAfterCoincidence=-1.0;
  (*w)->nonZeroFractionAfterSetMask=-1.0;
  (*w)->nonZeroFractionAfterClustering=-1.0;
  (*w)->nonZeroFractionAfterCuts=-1.0;
  (*w)->nonZeroFractionAfterVetoes=-1.0;
  (*w)->pixelSwapApplied=FALSE;
  (*w)->pixelMixerApplied=FALSE;
  (*w)->clusterType=0;
  (*w)->blobs=NULL;
  (*w)->delta_t=-1.0;
  (*w)->delta_f=-1.0;
}


static void _freeClusterWavelet(ClusterWavelet **w)
{
  UINT4 i;
  if((*w)->wavelet!=NULL)
    {
      _freeWavelet(&(*w)->wavelet);
      (*w)->wavelet=NULL;
    }
  if((*w)->original!=NULL)
    {
      _freeWavelet(&(*w)->original);
      (*w)->original=NULL;
    }
  if((*w)->medians!=NULL)
    {
      LALFree((*w)->medians);
      (*w)->medians=NULL;
    }
  if((*w)->norm50!=NULL)
    {
      LALFree((*w)->norm50);
      (*w)->norm50=NULL;
    }
  if((*w)->pMask!=NULL)
    {
      for(i=0;i<(*w)->pMaskCount;i++)
	{
	  if((*w)->pMask[i]!=NULL)
	    {
	      LALFree((*w)->pMask[i]);
	      (*w)->pMask[i]=NULL;
	    }
	}
      LALFree((*w)->pMask);
    }
  if((*w)->sCuts!=NULL)
    {
      LALFree((*w)->sCuts);
      (*w)->sCuts=NULL;
    }
  if((*w)->cList!=NULL)
    {
      for(i=0;i<(*w)->clusterCount;i++)
	{
	  if((*w)->cList[i]!=NULL)
	    {
	      LALFree((*w)->cList[i]);
	      (*w)->cList[i]=NULL;
	    }
	}
      LALFree((*w)->cList);
    }
  if((*w)->volumes!=NULL)
    {
      LALFree((*w)->volumes);
      (*w)->volumes=NULL;
    }
  if((*w)->coreSize!=NULL)
    {
      LALFree((*w)->coreSize);
      (*w)->coreSize=NULL;
    }
  if((*w)->correlation!=NULL)
    {
      LALFree((*w)->correlation);
      (*w)->correlation=NULL;
    }
  if((*w)->likelihood!=NULL)
    {
      LALFree((*w)->likelihood);
      (*w)->likelihood=NULL;
    }
  if((*w)->power!=NULL)
    {
      LALFree((*w)->power);
      (*w)->power=NULL;
    }
  if((*w)->maxAmplitude!=NULL)
    {
      LALFree((*w)->maxAmplitude);
      (*w)->maxAmplitude=NULL;
    }
  if((*w)->relativeStartTime!=NULL)
    {
      LALFree((*w)->relativeStartTime);
      (*w)->relativeStartTime=NULL;
    }
  if((*w)->relativeStopTime!=NULL)
    {
      LALFree((*w)->relativeStopTime);
      (*w)->relativeStopTime=NULL;
    }
  if((*w)->duration!=NULL)
    {
      LALFree((*w)->duration);
      (*w)->duration=NULL;
    }
  if((*w)->absoluteStartTime!=NULL)
    {
      LALFree((*w)->absoluteStartTime);
      (*w)->absoluteStartTime=NULL;
    }
  if((*w)->absoluteStopTime!=NULL)
    {
      LALFree((*w)->absoluteStopTime);
      (*w)->absoluteStopTime=NULL;
    }
  if((*w)->startFrequency!=NULL)
    {
      LALFree((*w)->startFrequency);
      (*w)->startFrequency=NULL;
    }
  if((*w)->stopFrequency!=NULL)
    {
      LALFree((*w)->stopFrequency);
      (*w)->stopFrequency=NULL;
    }
  if((*w)->bandwidth!=NULL)
    {
      LALFree((*w)->bandwidth);
      (*w)->bandwidth=NULL;
    }
  if((*w)->blobs!=NULL)
    {
      for(i=0;i<(*w)->clusterCount;i++)
	{
	  if((*w)->blobs[i].pBlob!=NULL) 
	    {
	      LALFree((*w)->blobs[i].pBlob);
	      (*w)->blobs[i].pBlob=NULL;
	    }
	  if((*w)->blobs[i].oBlob!=NULL) 
	    {
	      LALFree((*w)->blobs[i].oBlob);
	      (*w)->blobs[i].oBlob=NULL;
	    }
	}
      LALFree((*w)->blobs);
      (*w)->blobs=NULL;
    }
  (*w)->pMaskCount=0;
  (*w)->clusterCount=0;
  LALFree(*w);
}


static void _setAmplitudes(ClusterWavelet *w)
{
  UINT4 i,k,f,t;
  Slice s;


  bzero(w->wavelet->data->data->data,w->wavelet->data->data->length*sizeof(REAL4));

  for(i=0;i<w->pMaskCount;i++)
    {
      f=w->pMask[i]->frequency;
      t=w->pMask[i]->time;
      _getSliceF(f,w->wavelet,&s);
      k=s.start+s.step*t;
      w->wavelet->data->data->data[k]=w->pMask[i]->amplitude;
    }
}


#endif
