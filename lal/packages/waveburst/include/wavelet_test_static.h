#ifndef TEST_WAVELET_STATIC_H
#define TEST_WAVELET_STATIC_H

#include "lal/LALWavelet.h"

#define maxStrLen 100
#define  maxError 0.000000001

typedef struct
tagTest
{
  char waveletFileName[maxStrLen];
  char shouldFileName[maxStrLen];
  char paramsFileName[maxStrLen];
  BOOLEAN result;
} Test;

static void _readTestRecord(FILE *in, Test *t);
/*time series and its metadata*/
static void _readWavelet_TS(Wavelet **wavelet, FILE *in);
/*same as readWavelet + pMask*/
static void _readWavelet_TS_PM(ClusterWavelet **w, FILE *in);
/* everyting */
static void _readWavelet_ALL(ClusterWavelet **w, FILE *in);
static void _readREAL4TimeSeries(REAL4TimeSeries **t, char* fileName);
static int _compareREAL4TimeSeries(REAL4TimeSeries *is, REAL4TimeSeries *should, BOOLEAN debug);
static int _compareWavelet_MD(Wavelet *is, Wavelet *should, BOOLEAN debug);
static int _compareWavelets_MD_TS(Wavelet *is, Wavelet *should, BOOLEAN debug);
static int _compareWavelets_MD_TS_PM(ClusterWavelet *is, ClusterWavelet *should, BOOLEAN debug);
static int _compareWavelets_ALL(ClusterWavelet *is, ClusterWavelet *should, BOOLEAN debug);
static void _writeWavelet(Wavelet *w, FILE *out);
static void _writeREAL4TimeSeries(REAL4TimeSeries *t, FILE *out);
static void _writePixel(PixelWavelet *pix, FILE *out);


static void readTestRecord(FILE *in, Test *t)
{
  fscanf(in, "%s %s %s", t->waveletFileName, t->paramsFileName, t->shouldFileName);
}

static void readWavelet_TS(Wavelet **wavelet, FILE *in)
{
  int i;
  char tmp[maxStrLen];
  char *ptr;

  (*wavelet)=(Wavelet*)LALMalloc(sizeof(Wavelet));
  (*wavelet)->data=(REAL4TimeSeries*)LALMalloc(sizeof(REAL4TimeSeries));
  (*wavelet)->data->data=(REAL4Sequence*)LALMalloc(sizeof(REAL4Sequence));

  fscanf(in,"%d",&(*wavelet)->type);
  fscanf(in,"%d",&(*wavelet)->border);
  fscanf(in,"%d",&(*wavelet)->treeType);
  fscanf(in,"%d",&(*wavelet)->level);
  fscanf(in,"%d",&(*wavelet)->HPFilterLength);
  fscanf(in,"%d",&(*wavelet)->LPFilterLength);
  fscanf(in,"%s",(*wavelet)->data->name);
  fscanf(in,"%d",&(*wavelet)->data->epoch.gpsSeconds);
  fscanf(in,"%d",&(*wavelet)->data->epoch.gpsNanoSeconds);
  
  fscanf(in,"%s",tmp);
  (*wavelet)->data->deltaT=strtod(tmp,&ptr);

  fscanf(in,"%s",tmp);
  (*wavelet)->data->f0=strtod(tmp,&ptr);

  fscanf(in,"%d",&(*wavelet)->data->data->length);
  
  (*wavelet)->data->data->data=(REAL4*)LALCalloc((*wavelet)->data->data->length,sizeof(REAL4));

  for(i=0;i<(*wavelet)->data->data->length;i++)
    {
      fscanf(in,"%s",tmp);
      (*wavelet)->data->data->data[i]=strtod(tmp,&ptr);
    }
}


static void readWavelet_TS_PM(ClusterWavelet **w, FILE *in)
{
  UINT4 i,j;
  char tmp[maxStrLen];
  char *ptr;

  readWavelet_TS(&(*w)->wavelet,in);

  fscanf(in,"%d",&(*w)->pMaskCount);

  (*w)->pMask=(PixelWavelet**)LALMalloc(sizeof(PixelWavelet*) * (*w)->pMaskCount);

  for(i=0;i<(*w)->pMaskCount;i++)
  {
    (*w)->pMask[i]=(PixelWavelet*)LALMalloc(sizeof(PixelWavelet));
    fscanf(in,"%d",&(*w)->pMask[i]->time);
    fscanf(in,"%d",&(*w)->pMask[i]->frequency);
    fscanf(in,"%d",&(*w)->pMask[i]->clusterID);
    fscanf(in,"%d",&(*w)->pMask[i]->core);
    fscanf(in,"%d",&(*w)->pMask[i]->neighborsCount);
    fscanf(in,"%s",tmp);
    (*w)->pMask[i]->amplitude=strtod(tmp,&ptr);

    for(j=0;j<(*w)->pMask[i]->neighborsCount;j++)
      {
	fscanf(in,"%d",&(*w)->pMask[i]->neighbors[j]);
      }
  }
}


static void readWavelet_ALL(ClusterWavelet **w, FILE *in)
{
  UINT4 i,j;
  char tmp[maxStrLen];
  char *ptr;

  readWavelet_TS_PM(w,in);

  fscanf(in,"%d",&(*w)->clusterCount);

  (*w)->sCuts=(INT4*)LALMalloc(sizeof(INT4)*(*w)->clusterCount);

  for(i=0;i<(*w)->clusterCount;i++)
    {
      fscanf(in,"%d",&(*w)->sCuts[i]);
    }

  (*w)->volumes=(UINT4*)LALMalloc(sizeof(UINT4) * (*w)->clusterCount);

  for(i=0;i<(*w)->clusterCount;i++)
    {
      fscanf(in,"%d",&(*w)->volumes[i]);
    }

  (*w)->cList=(UINT4**)LALMalloc(sizeof(UINT4*) * (*w)->clusterCount);

  for(i=0;i<(*w)->clusterCount;i++)
    {
      (*w)->cList[i]=(UINT4*)LALMalloc(sizeof(UINT4) * (*w)->volumes[i]);

      for(j=0;j<(*w)->volumes[i];j++)
	{
	  fscanf(in,"%d",&(*w)->cList[i][j]);
	}
    }  
}

static void readREAL4TimeSeries(REAL4TimeSeries **t, char* fileName)
{
  FILE *in;
  int i;
  char tmp[maxStrLen];
  char *ptr;

  in=LALOpenDataFile(fileName);
  if(in==NULL)
    {
      printf("Cannot open %s\n",fileName);
      exit(1);
    }

  (*t)=(REAL4TimeSeries*)LALMalloc(sizeof(REAL4TimeSeries));

  (*t)->data=(REAL4Sequence*)LALMalloc(sizeof(REAL4Sequence));

  fscanf(in,"%s",(*t)->name);
  fscanf(in,"%d",&(*t)->epoch.gpsSeconds);
  fscanf(in,"%d",&(*t)->epoch.gpsNanoSeconds);

  fscanf(in,"%s",tmp);
  (*t)->deltaT=strtod(tmp,&ptr);

  fscanf(in,"%s",tmp);
  (*t)->f0=strtod(tmp,&ptr);

  fscanf(in,"%d",&(*t)->data->length);
  
  (*t)->data->data=(REAL4*)LALMalloc((*t)->data->length*sizeof(REAL4));
  
  for(i=0;i<(*t)->data->length;i++)
    {

      fscanf(in,"%s",tmp);
      (*t)->data->data[i]=strtod(tmp,&ptr);
    }

  LALFclose(in);
}


static int compareREAL4TimeSeries(REAL4TimeSeries *is, REAL4TimeSeries *should, BOOLEAN debug)
{
  int mismatch=0;
  int i;
  
  if(is->epoch.gpsSeconds!=should->epoch.gpsSeconds || debug) 
    {
      mismatch++;
      printf("Seconds: is=%d should=%d\n",is->epoch.gpsSeconds,should->epoch.gpsSeconds);
    }
  if(is->epoch.gpsNanoSeconds!=should->epoch.gpsNanoSeconds || debug) 
    {
      mismatch++;
      printf("Nano seconds: is=%d should=%d\n",is->epoch.gpsNanoSeconds,should->epoch.gpsNanoSeconds);
    }
  if(fabs(is->deltaT - should->deltaT) > maxError || debug) 
    { 
      mismatch++;
      printf("deltaT: is=%f should=%f\n", is->deltaT, should->deltaT);
    }
  if(fabs(is->f0 - should->f0) > maxError || debug)
    {
      mismatch++;
      printf("f0: is=%f should=%f\n", is->f0, should->f0);
    }
  if(is->data->length != should->data->length || debug)
    {
      mismatch++;
      printf("data->length: is=%d should=%d\n", is->data->length, should->data->length);
    }
  fflush(stdout);
  for(i=0;i<should->data->length;i++)
    {
      if(fabs(is->data->data[i] - should->data->data[i]) > maxError || debug)
	{
	  mismatch++;
	  printf("data[%d]: is=%f should=%f\n", i,is->data->data[i], should->data->data[i]);	  
	}
    }
  fflush(stdout);
  return mismatch;
}


/* compare wavelet metadata (MD) only */
static int compareWavelet_MD(Wavelet *is, Wavelet *should, BOOLEAN debug)
{
  int mismatch=0;

  if(is->type!=should->type || debug) 
    {
      mismatch++;
      printf("Wavelet type: is=%d should=%d\n",is->type,should->type);
    }
  if(is->border!=should->border || debug) 
    {
      mismatch++;
      printf("Wavelet border: is=%d should=%d\n",is->border,should->border);
    }
  if(is->treeType!=should->treeType || debug) 
    {
      mismatch++;
      printf("Wavelet treeType: is=%d should=%d\n",is->treeType,should->treeType);
    }
  if(is->level!=should->level || debug) 
    {
      mismatch++;
      printf("Wavelet level: is=%d should=%d\n",is->level,should->level);
    }
  if(is->HPFilterLength!=should->HPFilterLength || debug) 
    {
      mismatch++;
      printf("Wavelet HPFilterLength: is=%d should=%d\n",is->HPFilterLength,should->HPFilterLength);
    }
  if(is->LPFilterLength!=should->LPFilterLength || debug) 
    {
      mismatch++;
      printf("Wavelet LPFilterLength: is=%d should=%d\n",is->LPFilterLength,should->LPFilterLength);
    }

  return mismatch;
}

/* compare time series (TS) and wavelet metadata (MD) only (no pMask, cCusts, cList);*/
static int compareWavelets_MD_TS(Wavelet *is, Wavelet *should, BOOLEAN debug)
{
  int mismatch=0;
  mismatch += compareWavelet_MD(is,should,debug);
  mismatch += compareREAL4TimeSeries(is->data,should->data,debug);
  return mismatch;
}

static int compareWavelets_MD_TS_PM(ClusterWavelet *is, ClusterWavelet *should, BOOLEAN debug)
{
  int mismatch=0;
  UINT4 i,j;
  mismatch += compareWavelets_MD_TS(is->wavelet,should->wavelet,debug);


  if(is->pMaskCount!=should->pMaskCount || debug)
    {
      mismatch++;
      printf("Wavelet pMaskCount: is=%ld should=%ld\n",is->pMaskCount,should->pMaskCount);      
    }

  for(i=0;i<is->pMaskCount;i++)
    {
      if(is->pMask[i]->time!=should->pMask[i]->time || debug)
	{
	  mismatch++;
	  printf("Wavelet pMask[%ld]->time: is=%ld should=%ld\n",i,
		 is->pMask[i]->time,should->pMask[i]->time);
	}
      if(is->pMask[i]->frequency!=should->pMask[i]->frequency || debug)
	{
	  mismatch++;
	  printf("Wavelet pMask[%ld]->frequency: is=%ld should=%ld\n",i,
		 is->pMask[i]->frequency,should->pMask[i]->frequency);
	}
      if(is->pMask[i]->clusterID!=should->pMask[i]->clusterID || debug)
	{
	  mismatch++;
	  printf("Wavelet pMask[%ld]->clusterID: is=%ld should=%ld\n",i,
		 is->pMask[i]->clusterID,should->pMask[i]->clusterID);
	}
      if(is->pMask[i]->core!=should->pMask[i]->core || debug)
	{
	  mismatch++;
	  printf("Wavelet pMask[%ld]->core: is=%d should=%d\n",i,
		 is->pMask[i]->core,should->pMask[i]->core);
	}
      if(is->pMask[i]->neighborsCount!=should->pMask[i]->neighborsCount || debug)
	{
	  mismatch++;
	  printf("Wavelet pMask[%ld]->neighborsCount: is=%ld should=%ld\n",i,
		 is->pMask[i]->neighborsCount,should->pMask[i]->neighborsCount);
	}
      if(fabs(is->pMask[i]->amplitude - 
	      should->pMask[i]->amplitude) > maxError || debug)
	{
	  mismatch++;
	  printf("Wavelet pMask[%ld]->amplitude: is=%f should=%f\n",i,
		 is->pMask[i]->amplitude, should->pMask[i]->amplitude);
	}

      for(j=0;j<is->pMask[i]->neighborsCount;j++)
	{
	  if(is->pMask[i]->neighbors[j] != should->pMask[i]->neighbors[j] || debug)
	    {
	      mismatch++;
	      printf("Wavelet pMask[%ld]->neighbors[%ld]: is=%ld should=%ld\n",i,j,
		     is->pMask[i]->neighbors[j],should->pMask[i]->neighbors[j]);
	    }
	}
    }

  return mismatch;
}

static int compareWavelets_ALL(ClusterWavelet *is, ClusterWavelet *should, BOOLEAN debug)
{
  int mismatch=0;
  UINT4 i,j;

  mismatch += compareWavelets_MD_TS_PM(is,should,debug);
  if(is->clusterCount!=should->clusterCount || debug)
    {
      mismatch++;
      printf("Wavelet clusterCount: is=%ld should=%ld\n",
	     is->clusterCount,should->clusterCount);
    }

  for(i=0;i<is->clusterCount;i++)
    {
      if(is->sCuts[i]!=should->sCuts[i] || debug)
	{
	  mismatch++;
	  printf("Wavelet sCuts[%ld]: is=%d should=%d\n",
		 i,is->sCuts[i],should->sCuts[i]);
	}
      if(is->volumes[i]!=should->volumes[i] || debug)
	{
	  mismatch++;
	  printf("Wavelet volumes[%ld]: is=%ld should=%ld\n",
		 i,is->volumes[i],should->volumes[i]);
	}
      for(j=0;j<is->volumes[i];j++)
	{
	  if(is->cList[i][j]!=should->cList[i][j] || debug)
	    {
	      mismatch++;
	      printf("Wavelet cList[%ld][%ld]: is=%ld should=%ld\n",
		     i,j,is->cList[i][j],should->cList[i][j]);
	    }
	}
    }

  return mismatch;
}

static void writeWavelet(Wavelet *w, FILE *out)
{
  int i;
  fprintf(out,"%d\n",w->type);
  fprintf(out,"%d\n",w->border);
  fprintf(out,"%d\n",w->treeType);
  fprintf(out,"%d\n",w->level);
  fprintf(out,"%d\n",w->HPFilterLength);
  fprintf(out,"%d\n",w->LPFilterLength);
  fprintf(out,"%s\n",w->data->name);
  fprintf(out,"%d\n",w->data->epoch.gpsSeconds);
  fprintf(out,"%d\n",w->data->epoch.gpsNanoSeconds);
  fprintf(out,"%f\n",w->data->deltaT);
  fprintf(out,"%f\n",w->data->f0);
  fprintf(out,"%d\n",w->data->data->length);
  
  for(i=0;i<w->data->data->length;i++)
    {
      fprintf(out,"%f ",w->data->data->data[i]);
    }
  fprintf(out,"\n\n");
  fflush(stdout);
}

static void writeREAL4TimeSeries(REAL4TimeSeries *t, FILE *out)
{
  int i;
  fprintf(out,"%s\n",t->name);
  fprintf(out,"%d\n",t->epoch.gpsSeconds);
  fprintf(out,"%d\n",t->epoch.gpsNanoSeconds);
  fprintf(out,"%f\n",t->deltaT);
  fprintf(out,"%f\n",t->f0);
  fprintf(out,"%d\n",t->data->length);
  
  for(i=0;i<t->data->length;i++)
    {
      fprintf(out,"%f ",&t->data->data[i]);
    }
  fprintf(out,"\n\n");
}

static void printSlice(Slice *s)
{
  printf("Slice: start=%ld, size=%ld, step=%ld\n",s->start,s->size,s->step);
}

static void _writePixel(PixelWavelet *pix, FILE *out)
{
  fprintf(out,"t=%d f=%d cid=%d core=%d a=%g nc=%d\n",
	  pix->time, pix->frequency, pix->clusterID,
	  pix->core, pix->amplitude, pix->neighborsCount);
  fflush(out);
}



#endif
