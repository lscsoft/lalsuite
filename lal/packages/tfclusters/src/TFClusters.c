/*-----------------------------------------------------------------------
 *
 * File Name: TFClusters.c
 *
 * Author: Julien Sylvestre
 *
 * Revision: $Id$ 
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * TFClusters
 *
 * SYNOPSIS
 * Function to analyze a spectrogram by clustering analysis in order to
 * detect transients. The output is a list of events, significant according
 * to thresholds on the clusters size in the spectrogram, and the total power
 * in a cluster.
 *
 * DIAGNOSTICS
 *
 * CALLS
 *
 * NOTES
 * 
 *-----------------------------------------------------------------------
 */


/******** <lalVerbatim file="TFClustersCV"> ********
Author: Sylvestre, J
$Id$
********* </lalVerbatim> ********/


#include "lal/LALRCSID.h"

NRCSID (TFCLUSTERSC, "$Id$");


#include <lal/AVFactories.h>
#include <lal/TFClusters.h>
#include <math.h>

/***********************MAIN FUNCTIONS*********************************/


/******** <lalLaTeX file="TFClustersC"> ********
\noindent
Compute the spectrogram from a time series.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void 
LALComputeSpectrogram (
		       LALStatus *status, 
		       Spectrogram *out, 
		       TFPlaneParams *tspec, 
		       REAL4TimeSeries *tseries
		       )
}
\idx{LALComputeSpectrogram()}

\subsubsection*{Description}
Computes the spectrogram \texttt{*out} for the time series \texttt{*tseries}, using the parameters defined in \texttt{*tspec}. 
FFTs can overlap if \texttt{deltaT * timeBins} is larger than the time series duration; if they do overlap, a Welch window is applied.
The power is the norm square of the (normalized) discrete Fourier transform.

\subsubsection*{Uses}
\begin{verbatim}
LALCreateForwardRealFFTPlan()
LALCCreateVector()
LALForwardRealFFT()
LALCDestroyVector()
LALDestroyRealFFTPlan()
\end{verbatim}

\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> ********/

void 
LALComputeSpectrogram (
		       LALStatus *status, 
		       Spectrogram *out, 
		       TFPlaneParams *tspec, 
		       REAL4TimeSeries *tseries
		       )
{
  UINT4 sid, olap, NN, minF;
  /* INT4 i,j; */
  INT4 j, k;
  REAL8 *power;
  REAL4 *wwin = NULL, *tdat = NULL;

  /* LALWindowParams winParams; */
  TFPlaneParams *params;

  RealFFTPlan *pfwd = NULL;
  REAL4Vector Pvec;
  COMPLEX8Vector *Hvec = NULL;
  REAL8 norm;

  INITSTATUS (status, "LALComputeSpectrogram", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);

  
  /* Check return structure: out should have NULL pointer */
  ASSERT ( out, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP);
  ASSERT ( out->power == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
  ASSERT ( out->params == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);

  /* Check plane params */
  ASSERT ( tspec, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP);
  ASSERT ( tspec->timeBins > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  ASSERT ( tspec->freqBins > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  ASSERT ( tspec->deltaT > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  ASSERT ( tspec->flow >= 0.0, status, TFCLUSTERSH_EPOS, TFCLUSTERSH_MSGEPOS);


  /* Check time series */
  ASSERT ( tseries, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP);
  ASSERT ( tseries->data, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP);

  ASSERT ( (REAL8)(tseries->data->length) * tseries->deltaT <=
	   (REAL8)(tspec->timeBins) * tspec->deltaT,
	   status, TFCLUSTERSH_EINCOMP, TFCLUSTERSH_MSGEINCOMP);

  /* check compatibility */
  ASSERT ( (REAL8)tspec->freqBins <= 0.5 * tspec->deltaT / tseries->deltaT,
	   status, TFCLUSTERSH_EINCOMP, TFCLUSTERSH_MSGEINCOMP );

  /* determine overlap */
  NN = (UINT4)floor(tspec->deltaT / tseries->deltaT);

  olap = (NN * tspec->timeBins - tseries->data->length) / (tspec->timeBins - 1);

  ASSERT( tseries->data->length == olap + tspec->timeBins * (NN - olap),
	  status, TFCLUSTERSH_EINCOMP, TFCLUSTERSH_MSGEINCOMP );

  minF = (UINT4)floor(tspec->flow * tspec->deltaT);

  ASSERT( tspec->freqBins + minF == NN/2,
	  status, TFCLUSTERSH_EINCOMP, TFCLUSTERSH_MSGEINCOMP );

  /* copy tspec */
  params = (TFPlaneParams *)LALMalloc(sizeof(TFPlaneParams));
  if(!params)
    {ABORT(status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

  params->timeBins = tspec->timeBins;
  params->freqBins = tspec->freqBins;
  params->flow = tspec->flow; 
  params->deltaT = tspec->deltaT;

  out->params = params;


  /* memory for output */
  power = (REAL8 *)LALMalloc(params->timeBins * params->freqBins * sizeof(REAL8));
  if(!power)
    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

  out->power = power;


  /* Set Fourier plans */
  LALCreateForwardRealFFTPlan(status->statusPtr, &pfwd, NN, 0 );
  CHECKSTATUSPTR (status);

  Pvec.length = NN;
  
  LALCCreateVector( status->statusPtr, &Hvec, NN/2 + 1);
  CHECKSTATUSPTR (status);

  norm = (REAL8)NN;

  /* set window */
  if(olap) { 
    REAL4 nn2 = 0.5*(REAL4)NN;
    wwin = (REAL4 *)LALMalloc(NN * sizeof(REAL4));
    ASSERT(wwin, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );

    for(k=0; (int)k<(int)NN; k++) {
      wwin[k] = 1.0 - pow(((REAL4)k-nn2)/nn2,2.0);
    }

    tdat =  (REAL4 *)LALMalloc(NN * sizeof(REAL4));
    ASSERT(tdat, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );
  }


  /* Get TF representation */
  for(sid = 0; sid < (UINT4)params->timeBins; sid++) {
    if(olap) {
      memcpy(tdat, tseries->data->data + sid * (NN - olap), NN * sizeof(REAL4));
      for(k=0; (int)k<(int)NN; k++) {
	tdat[k] *= wwin[k];
      }
      Pvec.data = tdat;
    } else {
      Pvec.data = tseries->data->data + sid * NN;
    }

    LALForwardRealFFT(status->statusPtr, Hvec, &Pvec, pfwd);

    for(j=minF; (UINT4)j< minF + params->freqBins; j++)
      power[sid*params->freqBins + j - minF] = ((REAL8)Hvec->data[j].re * (REAL8)Hvec->data[j].re + (REAL8)Hvec->data[j].im * (REAL8)Hvec->data[j].im) / norm;
  }

  /*
  {
    UINT4 m;
    REAL8 me = 0.0;
    for(m=0; m<params->freqBins * params->timeBins; ++m) me += power[m];
    me /= (REAL8)m;
    printf("%g\n",me);
  }
  */

  if(olap) {
    LALFree(tdat);
    LALFree(wwin);
  }

  LALCDestroyVector( status->statusPtr, &Hvec );
  CHECKSTATUSPTR (status);

  LALDestroyRealFFTPlan( status->statusPtr, &pfwd );
  CHECKSTATUSPTR (status);

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}




/******** <lalLaTeX file="TFClustersC"> ********
\noindent
Compute the spectrogram from a time series.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void 
LALComputeXSpectrogram (
		       LALStatus *status, 
		       Spectrogram *out, 
		       TFPlaneParams *tspec, 
		       REAL4VectorTimeSeries *tseries
		       )
}
\idx{LALComputeXSpectrogram()}

\subsubsection*{Description}
Computes the cross-spectrogram \texttt{*out} for the time series \texttt{*tseries}, using the parameters defined in \texttt{*tspec}. 
FFTs can overlap if \texttt{deltaT * timeBins} is larger than the time series duration; if they do overlap, a Welch window is applied.
The power is the norm square of the (normalized) discrete Fourier transform.

\subsubsection*{Uses}
\begin{verbatim}
LALCreateForwardRealFFTPlan()
LALCCreateVector()
LALForwardRealFFT()
LALCDestroyVector()
LALDestroyRealFFTPlan()
\end{verbatim}

\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> ********/

void 
LALComputeXSpectrogram (
		       LALStatus *status, 
		       Spectrogram *out, 
		       TFPlaneParams *tspec, 
		       REAL4TimeVectorSeries *tseries
		       )
{
  UINT4 sid, olap, NN, minF;
  INT4 j, k;
  REAL8 *power;
  REAL4 *wwin = NULL, *tdat = NULL;

  /* LALWindowParams winParams; */
  TFPlaneParams *params;

  RealFFTPlan *pfwd = NULL;
  REAL4Vector Pvec;
  COMPLEX8Vector *Hvec = NULL;
  COMPLEX8Vector *Ivec = NULL;
  REAL8 norm;

  INITSTATUS (status, "LALComputeSpectrogram", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);

  
  /* Check return structure: out should have NULL pointer */
  ASSERT ( out, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP);
  ASSERT ( out->power == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
  ASSERT ( out->params == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);

  /* Check plane params */
  ASSERT ( tspec, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP);
  ASSERT ( tspec->timeBins > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  ASSERT ( tspec->freqBins > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  ASSERT ( tspec->deltaT > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  ASSERT ( tspec->flow >= 0.0, status, TFCLUSTERSH_EPOS, TFCLUSTERSH_MSGEPOS);


  /* Check time series */
  ASSERT ( tseries, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP);
  ASSERT ( tseries->data, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP);

  ASSERT ( (REAL8)(tseries->data->vectorLength) * tseries->deltaT <=
	   (REAL8)(tspec->timeBins) * tspec->deltaT,
	   status, TFCLUSTERSH_EINCOMP, TFCLUSTERSH_MSGEINCOMP);
  ASSERT ( tseries->data->length == 2, status, TFCLUSTERSH_EINCOMP, TFCLUSTERSH_MSGEINCOMP);
  ASSERT ( tseries->data->data, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP);

  /* check compatibility */
  ASSERT ( (REAL8)tspec->freqBins <= 0.5 * tspec->deltaT / tseries->deltaT,
	   status, TFCLUSTERSH_EINCOMP, TFCLUSTERSH_MSGEINCOMP );

  /* determine overlap */
  NN = (UINT4)floor(tspec->deltaT / tseries->deltaT);

  olap = (NN * tspec->timeBins - tseries->data->vectorLength) / (tspec->timeBins - 1);

  ASSERT( tseries->data->vectorLength == olap + tspec->timeBins * (NN - olap),
	  status, TFCLUSTERSH_EINCOMP, TFCLUSTERSH_MSGEINCOMP );

  minF = (UINT4)floor(tspec->flow * tspec->deltaT);

  ASSERT( tspec->freqBins + minF == NN/2,
	  status, TFCLUSTERSH_EINCOMP, TFCLUSTERSH_MSGEINCOMP );

  /* copy tspec */
  params = (TFPlaneParams *)LALMalloc(sizeof(TFPlaneParams));
  if(!params)
    {ABORT(status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

  params->timeBins = tspec->timeBins;
  params->freqBins = tspec->freqBins;
  params->flow = tspec->flow; 
  params->deltaT = tspec->deltaT;

  out->params = params;


  /* memory for output */
  power = (REAL8 *)LALMalloc(params->timeBins * params->freqBins * sizeof(REAL8));
  if(!power)
    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

  out->power = power;


  /* Set Fourier plans */
  LALCreateForwardRealFFTPlan(status->statusPtr, &pfwd, NN, 0 );
  CHECKSTATUSPTR (status);

  Pvec.length = NN;
  
  LALCCreateVector( status->statusPtr, &Hvec, NN/2 + 1);
  CHECKSTATUSPTR (status);

  LALCCreateVector( status->statusPtr, &Ivec, NN/2 + 1);
  CHECKSTATUSPTR (status);

  norm = (REAL8)NN;

  /* set window */
  if(olap) { 
    REAL4 nn2 = 0.5*(REAL4)NN;
    wwin = (REAL4 *)LALMalloc(NN * sizeof(REAL4));
    ASSERT(wwin, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );

    for(k=0; (int)k<(int)NN; k++) {
      wwin[k] = 1.0 - pow(((REAL4)k-nn2)/nn2,2.0);
    }

    tdat =  (REAL4 *)LALMalloc(NN * sizeof(REAL4));
    ASSERT(tdat, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );
  }


  /* Get TF representation */
  for(sid = 0; sid < (UINT4)params->timeBins; sid++) {
    if(olap) {
      memcpy(tdat, tseries->data->data + sid * (NN - olap), NN * sizeof(REAL4));
      for(k=0; (int)k<(int)NN; k++) {
	tdat[k] *= wwin[k];
      }
      Pvec.data = tdat;
    } else {
      Pvec.data = tseries->data->data + sid * NN;
    }

    LALForwardRealFFT(status->statusPtr, Hvec, &Pvec, pfwd);


    
    if(olap) {
      memcpy(tdat, tseries->data->data + tseries->data->vectorLength + sid * (NN - olap), NN * sizeof(REAL4));
      for(k=0; (int)k<(int)NN; k++) {
	tdat[k] *= wwin[k];
      }
      Pvec.data = tdat;
    } else {
      Pvec.data = tseries->data->data + sid * NN;
    }

    LALForwardRealFFT(status->statusPtr, Ivec, &Pvec, pfwd);



    for(j=minF; (UINT4)j< minF + params->freqBins; j++)
      power[sid*params->freqBins + j - minF] = sqrt( pow((REAL8)Hvec->data[j].re * (REAL8)Ivec->data[j].re + (REAL8)Hvec->data[j].im * (REAL8)Ivec->data[j].im, 2.0) + pow((REAL8)Hvec->data[j].im * (REAL8)Ivec->data[j].re - (REAL8)Hvec->data[j].re * (REAL8)Ivec->data[j].im, 2.0) ) / norm;
  }

  if(olap) {
    LALFree(tdat);
    LALFree(wwin);
  }

  LALCDestroyVector( status->statusPtr, &Hvec );
  CHECKSTATUSPTR (status);

  LALCDestroyVector( status->statusPtr, &Ivec );
  CHECKSTATUSPTR (status);

  LALDestroyRealFFTPlan( status->statusPtr, &pfwd );
  CHECKSTATUSPTR (status);

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/******** <lalLaTeX file="TFClustersC"> ********
\newpage
\noindent
Apply the first two levels of thresholds: (i) cut on power of individual pixels in spectrogram and (ii) cut on cluster sizes in thresholded spectrogram.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void 
LALGetClusters (
		LALStatus *status, 
		CList *clist, 
		Spectrogram *tpower, 
		CListDir *dir
		)
}
\idx{LALGetClusters()}

\subsubsection*{Description}
First, this function transforms \texttt{*tpower} into a binary map, by applying the frequency dependent thresholds \texttt{dir->rho} on the power in the spectrogram. Only frequencies up to \texttt{dir->maxf} are retained. A recursive function is then called to identify the clusters on a `nearest neighbours' basis (i.e., pixels touching by one `edge'). Only clusters with power strictly between \texttt{dir->minf} and \texttt{dir->maxf} are kept (if \texttt{dir->maxf} is negative, only clusters with at least some power between \texttt{dir->minf} and \texttt{dir->maxf} are kept). Clusters larger or equal to \texttt{dir->sigma} are sent to \texttt{*clist}. The remaining clusters are grouped in pairs. Whenever a pair pass the distance thresholds defined by \texttt{dir->s1}, \texttt{dir->s2} and \texttt{dir->d}, the two clusters are fused and are added as a single cluster to \texttt{*clist}. 

\subsubsection*{Uses}
\begin{verbatim}
LALInitCList()
LALFreeCList()
\end{verbatim}

\subsubsection*{Notes}
\begin{itemize}
\item \texttt{*clist} must be initialized by a proper call to \texttt{LALInitCList()} before calling this function.
\item Calling this function destroys \texttt{*tpower}.
\item \texttt{dir->rho[0]} corresponds to \texttt{minf}, not DC. 
\end{itemize}

\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> ********/
static void GetNearestNeighb(CList *, Spectrogram *, REAL8 *, UINT4, UINT4);

static UINT4 ClustDistance(UINT4 s1, UINT4 *t1, UINT4 *f1, UINT4 s2, UINT4 *t2, UINT4 *f2);

static UINT4 min(UINT4 a, UINT4 b);
  
static UINT4 max(UINT4 a, UINT4 b);


void 
LALGetClusters (
		LALStatus *status, 
		CList *clist, 
		Spectrogram *tpower, 
		CListDir *dir
		)
{
  UINT4 i,j, myj0, mins, ndir, dpos, s1, s2, nfriends, nbigs, nt, jmax=0, jmin;
  UINT4 *friends1, *friends2, *bigs, *where, *ts;
  UINT4 **bbox, **tf, **tt;
  REAL8 *rho = dir->rho;
  REAL8 **tP;
  CList tlist;
  CHAR strictFCut = 1;
  REAL8 minf, maxf;

  INITSTATUS (status, "LALGetClusters", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);

  
  /* check I/O structure */
  ASSERT ( dir, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( tpower, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( clist, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  
  ASSERT ( tpower->power, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( tpower->params, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );

  ASSERT ( dir->rho, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );

/*
  ASSERT ( dir->maxf > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  ASSERT ( dir->minf >= 0, status, TFCLUSTERSH_EPOS, TFCLUSTERSH_MSGEPOS);
  ASSERT ( dir->minf <= dir->maxf, status, TFCLUSTERSH_EIARG, TFCLUSTERSH_MSGEIARG );
  */

  if(dir->sigma > 1) {
    ASSERT ( dir->s1, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( dir->s2, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( dir->d, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  }

  ASSERT ( clist->params, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( clist->nclusters == 0, status, TFCLUSTERSH_ENZERO, TFCLUSTERSH_MSGENZERO);
  ASSERT ( clist->sizes == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
  ASSERT ( clist->t == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
  ASSERT ( clist->f == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
  ASSERT ( clist->P == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);


  
  ndir = dir->sigma * (dir->sigma - 1) / 2;

  for(mins=-1, i=0;i<ndir;i++)
    if(dir->d[i] > 0 && dir->s1[i] < mins) mins=dir->s1[i]; 
      


  LALInitCList(status->statusPtr, &tlist, clist->params);
  CHECKSTATUSPTR (status);

  tlist.nclusters = 0;
  tlist.sizes = NULL;
  tlist.t = tlist.f = NULL;
  tlist.P = NULL;

  

  /* Get the raw clusters */
  for(i=0;i<(UINT4)tpower->params->timeBins;i++) /* loop over all points */
    for(j=0;j<(UINT4)tpower->params->freqBins;j++)
      if(tpower->power[i*tpower->params->freqBins + j] > rho[j])
	/* note: this makes the bursts roughly time ordered */
	{
	  tlist.nclusters++;

	  tlist.sizes = (UINT4 *)LALRealloc(tlist.sizes, tlist.nclusters*sizeof(UINT4));
	  ASSERT ( tlist.sizes, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );
	  
	  tlist.t = (UINT4 **)LALRealloc(tlist.t, tlist.nclusters*sizeof(UINT4 *));
	  ASSERT ( tlist.t, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );

	  tlist.f = (UINT4 **)LALRealloc(tlist.f, tlist.nclusters*sizeof(UINT4 *));
	  ASSERT ( tlist.f, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );

	  tlist.P = (REAL8 **)LALRealloc(tlist.P, tlist.nclusters*sizeof(REAL8 *));
	  ASSERT ( tlist.P, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );


	  tlist.sizes[tlist.nclusters - 1] = 0;
	  tlist.t[tlist.nclusters - 1] = NULL;
	  tlist.f[tlist.nclusters - 1] = NULL;
	  tlist.P[tlist.nclusters - 1] = NULL;

	  GetNearestNeighb(&tlist, tpower, rho, i, j);
	}        


  
  /* clean to minf and maxf */
  if(dir->maxf < 0.0) {
    strictFCut = 0;
    maxf = -1.0*dir->maxf;
    minf = -1.0*dir->minf;
  } else {
    strictFCut = 1;
    maxf = dir->maxf;
    minf = dir->minf;
  }

  if(tpower->params->flow + (double)tpower->params->freqBins / tpower->params->deltaT >= maxf &&
     maxf >= tpower->params->flow) {
    REAL8 tmp = (maxf - tpower->params->flow) * tpower->params->deltaT;
    if(tmp - floor(tmp) >= 0.5) jmax = (UINT4)ceil(tmp);
    else jmax = (UINT4)floor(tmp);
  }
  else
    jmax=-1U;


  if(tpower->params->flow + (double)tpower->params->freqBins / tpower->params->deltaT >= minf &&
     minf >= tpower->params->flow) {
    REAL8 tmp = (minf - tpower->params->flow) * tpower->params->deltaT;
    if(tmp - floor(tmp) >= 0.5) jmin = (UINT4)ceil(tmp);
    else jmin = (UINT4)floor(tmp);
  }
  else
    jmin=0;

  /*
  printf("%u\t\t%u\n",jmin,jmax);
  */

  if(jmin>0 || jmax!=-1U)
    {
      ts = NULL;
      tf=tt=NULL;
      tP=NULL;
      nt=0;

      for(i=0;i<tlist.nclusters;i++)
	{
	  BOOLEAN keepIt = strictFCut;

	  for(j=0;j<tlist.sizes[i];j++) {
	    if(strictFCut && (
			      tlist.f[i][j] > jmax ||
			      tlist.f[i][j] < jmin)) {
	      keepIt = 0;
	      break;
	    }

	    if(!strictFCut && (
			       tlist.f[i][j] <= jmax &&
			       tlist.f[i][j] >= jmin)) {
	      keepIt = 1;
	      break;
	    }
	  }

	 if(keepIt) /* keep it */
	  {
	    nt++;
	    ts=(UINT4 *)LALRealloc(ts, nt*sizeof(UINT4));
	    ASSERT ( ts, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );

	    tf=(UINT4 **)LALRealloc(tf, nt*sizeof(UINT4 *));
	    ASSERT ( tf, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );

	    tt=(UINT4 **)LALRealloc(tt, nt*sizeof(UINT4 *));
	    ASSERT ( tt, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );

	    tP=(REAL8 **)LALRealloc(tP, nt*sizeof(REAL8 *));
	    ASSERT ( tP, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );


	    ts[nt-1] = tlist.sizes[i];
	    tf[nt-1] = tlist.f[i];
	    tt[nt-1] = tlist.t[i];
	    tP[nt-1] = tlist.P[i];
	  }
	else  /* trash it */
	  {
	    LALFree(tlist.f[i]);
	    LALFree(tlist.t[i]);
	    LALFree(tlist.P[i]);
	  }
	}
	
      if(tlist.nclusters > 0)
	{LALFree(tlist.sizes);
	LALFree(tlist.f);
	LALFree(tlist.t);
	LALFree(tlist.P);}

      tlist.nclusters=nt;
      tlist.sizes=ts;
      tlist.f = tf;
      tlist.t = tt;
      tlist.P = tP;	  
    }


  
  /* clean to mins */
  if(mins > 1 && mins != -1U)
    {
      ts = NULL;
      tf=tt=NULL;
      tP=NULL;
      nt=0;

      for(i=0;i<tlist.nclusters;i++)
	if(tlist.sizes[i] >= mins)
	  {
	    nt++;
	    ts=(UINT4 *)LALRealloc(ts, nt*sizeof(UINT4));
	    ASSERT ( ts, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );
	    tf=(UINT4 **)LALRealloc(tf, nt*sizeof(UINT4 *));
	    ASSERT ( tf, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );
	    tt=(UINT4 **)LALRealloc(tt, nt*sizeof(UINT4 *));
	    ASSERT ( tt, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );
	    tP=(REAL8 **)LALRealloc(tP, nt*sizeof(REAL8 *));
	    ASSERT ( tP, status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );

	    ts[nt-1] = tlist.sizes[i];
	    tf[nt-1] = tlist.f[i];
	    tt[nt-1] = tlist.t[i];
	    tP[nt-1] = tlist.P[i];
	  }
	else
	  {
	    LALFree(tlist.f[i]);
	    LALFree(tlist.t[i]);
	    LALFree(tlist.P[i]);
	  }

      if(tlist.nclusters > 0)
	{LALFree(tlist.sizes);
	LALFree(tlist.f);
	LALFree(tlist.t);
	LALFree(tlist.P);}

      tlist.nclusters=nt;
      tlist.sizes=ts;
      tlist.f = tf;
      tlist.t = tt;
      tlist.P = tP;
    }



  /* clean using distances */
  nfriends = nbigs = 0;
  friends1 = friends2 = bigs = NULL;
  bbox = NULL;

  if(tlist.nclusters > 0) {
  bbox=(UINT4 **)LALMalloc(tlist.nclusters * sizeof(UINT4 *));
  if(!bbox)
    {ABORT( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

  for(i=0;i<tlist.nclusters;i++)
    {
      bbox[i] = (UINT4 *)LALMalloc(4*sizeof(UINT4));
      if(!bbox[i])
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

      bbox[i][0] = bbox[i][1] = tlist.t[i][0];
      bbox[i][2] = bbox[i][3] = tlist.f[i][0];

      for(j=1;j<tlist.sizes[i];j++)
	{
	  if(tlist.t[i][j] < bbox[i][0]) bbox[i][0] = tlist.t[i][j];
	  if(tlist.t[i][j] > bbox[i][1]) bbox[i][1] = tlist.t[i][j];
	  if(tlist.f[i][j] < bbox[i][2]) bbox[i][2] = tlist.f[i][j];
	  if(tlist.f[i][j] > bbox[i][3]) bbox[i][3] = tlist.f[i][j];
	}
    }

  for(i=0;i<tlist.nclusters;i++)
    {
      if(tlist.sizes[i] >= dir->sigma)
	{
	  nbigs++;

	  bigs = (UINT4 *)LALRealloc(bigs, nbigs * sizeof(UINT4));
	  if(!bigs)
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

	  bigs[nbigs-1] = i;
	}
      else
	for(j=i+1;j<tlist.nclusters;j++)
	  {
	    if((int)bbox[j][0] - (int)bbox[i][1] > (int)dir->mdist) break; /*don't care about clusters that start more than ? pixels after the end of the actual one */ 


	  if(tlist.sizes[j] < dir->sigma)
	    {
	      s1 = min(tlist.sizes[i], tlist.sizes[j]);
	      s2 = max(tlist.sizes[i], tlist.sizes[j]);
	    
	      dpos = s2 - s1 + (2*dir->sigma - s1) * (s1 - 1) / 2;

	      if(dir->d[dpos] > 0)
		{
		  if(min(abs((int)bbox[i][0] - (int)bbox[j][1]), abs((int)bbox[i][1] - (int)bbox[j][0])) <= dir->d[dpos] &&
		     min(abs((int)bbox[i][2] - (int)bbox[j][3]), abs((int)bbox[i][3] - (int)bbox[j][2])) <= dir->d[dpos]) 
		    {
		      if( ClustDistance(tlist.sizes[i], tlist.t[i], tlist.f[i], tlist.sizes[j], tlist.t[j], tlist.f[j]) <= dir->d[dpos])
			{
			  nfriends++;

			  friends1 = (UINT4 *)LALRealloc(friends1, nfriends*sizeof(UINT4));
			  if(!friends1)
			    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

			  friends2 = (UINT4 *)LALRealloc(friends2, nfriends*sizeof(UINT4));
			  if(!friends2)
			    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
			  friends1[nfriends-1] = i;
			  friends2[nfriends-1] = j;
			}	
		    }
		}
	    }
	  }

    }
  }




  /* generate output list */
  clist->nclusters = 0;
  clist->sizes = NULL;
  clist->t = clist->f = NULL;
  clist->P = NULL;

  for(i=0;i<nbigs;i++)
    {
      clist->nclusters++;

      clist->sizes = (UINT4 *)LALRealloc(clist->sizes, clist->nclusters*sizeof(UINT4));
      if(!(clist->sizes))
	{ABORT (  status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

      clist->sizes[clist->nclusters-1] = tlist.sizes[bigs[i]];

      clist->f = (UINT4 **)LALRealloc(clist->f, clist->nclusters*sizeof(UINT4 *));
      if(!(clist->f))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

      clist->t = (UINT4 **)LALRealloc(clist->t, clist->nclusters*sizeof(UINT4 *));
      if(!(clist->t))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

      clist->P = (REAL8 **)LALRealloc(clist->P, clist->nclusters*sizeof(REAL8 *));
      if(!(clist->P))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

      clist->f[clist->nclusters-1] = (UINT4 *)LALMalloc(clist->sizes[clist->nclusters-1] * sizeof(UINT4));
      if(!(clist->f[clist->nclusters-1]))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

      clist->t[clist->nclusters-1] = (UINT4 *)LALMalloc(clist->sizes[clist->nclusters-1] * sizeof(UINT4));
      if(!(clist->t[clist->nclusters-1]))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

      clist->P[clist->nclusters-1] = (REAL8 *)LALMalloc(clist->sizes[clist->nclusters-1] * sizeof(REAL8));
      if(!(clist->P[clist->nclusters-1]))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}


      for(j=0;j<clist->sizes[clist->nclusters-1];j++)
	{clist->f[clist->nclusters-1][j] = tlist.f[bigs[i]][j];
	clist->t[clist->nclusters-1][j] = tlist.t[bigs[i]][j];
	clist->P[clist->nclusters-1][j] = tlist.P[bigs[i]][j];}
    }



  
  /* friends */

  if(tlist.nclusters > 0) {
  where = (UINT4 *)LALMalloc(tlist.nclusters * sizeof(UINT4));
  if(!where)
    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

  for(i=0;i<tlist.nclusters;i++)
    where[i] = -1;

  for(i=0;i<nfriends;i++)
    {
      if(where[friends1[i]] == -1U &&
	 where[friends2[i]] == -1U) /* create new burst */
	{
	  where[friends1[i]] = clist->nclusters;
	  where[friends2[i]] = clist->nclusters;

	  clist->nclusters++;
	  
	  clist->sizes = (UINT4 *)LALRealloc(clist->sizes, clist->nclusters*sizeof(UINT4));
	  if(!(clist->sizes))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

	  clist->sizes[clist->nclusters-1] = tlist.sizes[friends1[i]] + tlist.sizes[friends2[i]];

	  clist->f = (UINT4 **)LALRealloc(clist->f, clist->nclusters*sizeof(UINT4 *));
	  if(!(clist->f))
	    {ABORT (  status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

	  clist->t = (UINT4 **)LALRealloc(clist->t, clist->nclusters*sizeof(UINT4 *));
	  if(!(clist->t))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

	  clist->P = (REAL8 **)LALRealloc(clist->P, clist->nclusters*sizeof(REAL8 *));
	  if(!(clist->P))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}


	  clist->f[clist->nclusters-1] = (UINT4 *)LALMalloc(clist->sizes[clist->nclusters-1] * sizeof(UINT4));
	  if(!(clist->f[clist->nclusters-1]))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

	  clist->t[clist->nclusters-1] = (UINT4 *)LALMalloc(clist->sizes[clist->nclusters-1] * sizeof(UINT4));
	  if(!(clist->t[clist->nclusters-1]))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

	  clist->P[clist->nclusters-1] = (REAL8 *)LALMalloc(clist->sizes[clist->nclusters-1] * sizeof(REAL8));
	  if(!(clist->P[clist->nclusters-1]))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

	  for(j=0;j<tlist.sizes[friends1[i]];j++)
	    {clist->f[clist->nclusters-1][j] = tlist.f[friends1[i]][j];
	    clist->t[clist->nclusters-1][j] = tlist.t[friends1[i]][j];
	    clist->P[clist->nclusters-1][j] = tlist.P[friends1[i]][j];}

	  for(;j<tlist.sizes[friends1[i]]+tlist.sizes[friends2[i]];j++)
	    {clist->f[clist->nclusters-1][j] = tlist.f[friends2[i]][j - tlist.sizes[friends1[i]]];
	    clist->t[clist->nclusters-1][j] = tlist.t[friends2[i]][j - tlist.sizes[friends1[i]]];
	    clist->P[clist->nclusters-1][j] = tlist.P[friends2[i]][j - tlist.sizes[friends1[i]]];}
	}

      if(where[friends1[i]] != -1U &&
	 where[friends2[i]] == -1U) /* append second */
	{
	  where[friends2[i]] = where[friends1[i]];

	  myj0 = clist->sizes[where[friends1[i]]];

	  clist->sizes[where[friends1[i]]] += tlist.sizes[friends2[i]];
	  
	  clist->f[where[friends1[i]]] = (UINT4 *)LALRealloc(clist->f[where[friends1[i]]], clist->sizes[where[friends1[i]]] * sizeof(UINT4));	  
	  if(!(clist->f[where[friends1[i]]]))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

	  clist->t[where[friends1[i]]] = (UINT4 *)LALRealloc(clist->t[where[friends1[i]]], clist->sizes[where[friends1[i]]] * sizeof(UINT4));
	  if(!(clist->t[where[friends1[i]]]))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

	  clist->P[where[friends1[i]]] = (REAL8 *)LALRealloc(clist->P[where[friends1[i]]], clist->sizes[where[friends1[i]]] * sizeof(REAL8));
	  if(!(clist->P[where[friends1[i]]]))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}


	  for(j=0;j<tlist.sizes[friends2[i]];j++)
	    {clist->f[where[friends1[i]]][j+myj0] = tlist.f[friends2[i]][j];
	    clist->t[where[friends1[i]]][j+myj0] = tlist.t[friends2[i]][j];
	    clist->P[where[friends1[i]]][j+myj0] = tlist.P[friends2[i]][j];}
	}      

      if(where[friends1[i]] == -1U &&
	 where[friends2[i]] != -1U) /* append first */
	{
	  where[friends1[i]] = where[friends2[i]];

	  myj0 = clist->sizes[where[friends2[i]]];

	  clist->sizes[where[friends2[i]]] += tlist.sizes[friends1[i]];
	  
	  clist->f[where[friends2[i]]] = (UINT4 *)LALRealloc(clist->f[where[friends2[i]]], clist->sizes[where[friends2[i]]] * sizeof(UINT4));	
	  if(!(clist->f[where[friends2[i]]]))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
  
	  clist->t[where[friends2[i]]] = (UINT4 *)LALRealloc(clist->t[where[friends2[i]]], clist->sizes[where[friends2[i]]] * sizeof(UINT4));
	  if(!(clist->t[where[friends2[i]]]))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

	  clist->P[where[friends2[i]]] = (REAL8 *)LALRealloc(clist->P[where[friends2[i]]], clist->sizes[where[friends2[i]]] * sizeof(REAL8));
	  if(!(clist->P[where[friends2[i]]]))
	    { ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}


	  for(j=0;j<tlist.sizes[friends1[i]];j++)
	    {clist->f[where[friends2[i]]][j+myj0] = tlist.f[friends1[i]][j];
	    clist->t[where[friends2[i]]][j+myj0] = tlist.t[friends1[i]][j];
	    clist->P[where[friends2[i]]][j+myj0] = tlist.P[friends1[i]][j];}
	}      
    }

  if(where != NULL)
    LALFree(where);
  if(friends1 != NULL)
    LALFree(friends1);
  if(friends2 != NULL)
    LALFree(friends2);
  if(bigs != NULL)
    LALFree(bigs);
  }

  for(i=0;i<tlist.nclusters;i++)
    LALFree(bbox[i]);
  if(tlist.nclusters > 0)
    LALFree(bbox);


  LALFreeCList(status->statusPtr, &tlist);
  CHECKSTATUSPTR (status);


  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}





/******** <lalLaTeX file="TFClustersC"> ********
\newpage
\noindent
Apply the final cut by thresholding on the total power in the clusters.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void 
LALClustersPowerThreshold (
			   LALStatus *status, 
			   CList *out,
			   CList *in, 
			   CListDir *dir
			   )
}
\idx{LALClustersPowerThreshold()}

\subsubsection*{Description}
This function loops over all clusters in \texttt{*in}; for each cluster it computes its total power by summing over the pixels of the cluster, and computes the probability for Gaussian noise to produce a cluster with this total power at this stage of the analysis. This probability is compared to \texttt{dir->alpha}; if smaller, the cluster from \texttt{*in} is appended to \texttt{*out}. Therefore, \texttt{dir->alpha} is the fraction of clusters that had survive the first cuts that will pass this one, assuming Gaussian noise as input of the algorithm. When \texttt{dir->alpha} < 0, only clusters which have at least one pixel with power larger or equal to -\texttt{dir->alpha} times the first power threshold will survive.

\subsubsection*{Notes}
\begin{itemize}
\item \texttt{*out} must be initialized by a proper call to \texttt{LALInitCList()} before calling this function.
\end{itemize}

\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> ********/
static void incgam(LALStatus *, REAL4, REAL4, REAL4 *);

void 
LALClustersPowerThreshold (
			   LALStatus *status, 
			   CList *out,
			   CList *in, 
			   CListDir *dir
			   )
{
  BOOLEAN winner;
  UINT4 i,j;
  REAL4 po, P0;
  REAL4 prob /* , norm */ ;

  INITSTATUS (status, "LALClustersPowerThreshold", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);


  /* check I/O structure */
  ASSERT ( out, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( in, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( dir, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );

  ASSERT ( out->nclusters == 0, status, TFCLUSTERSH_ENZERO, TFCLUSTERSH_MSGENZERO);
  ASSERT ( out->sizes == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
  ASSERT ( out->t == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
  ASSERT ( out->f == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
  ASSERT ( out->P == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);

  if(dir->alpha >= 1) {
    LALCopyCList(status->statusPtr, out, in);
    CHECKSTATUSPTR (status);
    DETATCHSTATUSPTR (status);
    RETURN (status); 
  }

  out->nclusters = 0;
  out->t = out->f = NULL;
  out->sizes = NULL;
  out->P = NULL;


  ASSERT ( dir->rho, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  /*
  ASSERT ( dir->maxf > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  */
  if(dir->sigma > 1) {
    ASSERT ( dir->s1, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( dir->s2, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( dir->d, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  }
  ASSERT ( dir->alpha <= 1, status, TFCLUSTERSH_E01, TFCLUSTERSH_MSGE01 );

  ASSERT ( out->params, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( out->params->timeBins == in->params->timeBins &&
	   out->params->freqBins == in->params->freqBins &&
	   out->params->deltaT == in->params->deltaT &&
	   out->params->flow == in->params->flow,
	   status, TFCLUSTERSH_EIARG, TFCLUSTERSH_MSGEIARG );

  /*  ASSERT ( in->nclusters >= 0, status, TFCLUSTERSH_EPOS, TFCLUSTERSH_MSGEPOS); */
  if(!(in->nclusters)) {
    DETATCHSTATUSPTR (status);
    RETURN (status);
  }

  ASSERT ( in->sizes, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( in->t, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( in->f, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( in->P, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( in->params, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );

  for(i=0; i<in->nclusters; i++) { /* loop over input clusters */
  
    winner = 0;

    if(dir->alpha > 0) {
      /* run the power threshold test */
      for(po = 0.0, P0=0.0, j = 0; j<in->sizes[i]; j++) {
	po += in->P[i][j];
	P0 += dir->rho[in->f[i][j]];
      }
     
      po -= P0;

      incgam(status->statusPtr, (float)in->sizes[i], po, &prob);
      CHECKSTATUSPTR (status);

      if(prob < dir->alpha) { 
	winner = 1;
      }
    } else {
      /* run the maximum pixel power test */
      for(j = 0; j<in->sizes[i]; j++) {
	po = -dir->alpha * dir->rho[in->f[i][j]];
	if(in->P[i][j] > po) {
	  winner = 1;
	  break;
	}
      }
    }

    if(winner) { /* we have a winner */

      (out->nclusters)++;

      /* copy input cluster to ouput list */
      out->sizes = (UINT4 *)LALRealloc(out->sizes, out->nclusters * sizeof(UINT4));
      if(!(out->sizes))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC); }

      out->sizes[out->nclusters - 1] = in->sizes[i];

      out->t = (UINT4**)LALRealloc(out->t, out->nclusters * sizeof(UINT4*));
      if(!(out->t))
	{ABORT ( status,TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC);}
      
      out->f = (UINT4**)LALRealloc(out->f, out->nclusters * sizeof(UINT4*));
      if(!(out->f)) {
	ABORT ( status,TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC);
      }

      out->P = (REAL8**)LALRealloc(out->P, out->nclusters * sizeof(REAL8*));
      if(!(out->P)) {
	ABORT ( status,TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC);
      }

      out->t[out->nclusters-1] = (UINT4 *)LALMalloc(in->sizes[i] * sizeof(UINT4));
      if(!(out->t[out->nclusters-1])) {
	ABORT ( status,TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC);
      }

      out->f[out->nclusters-1] = (UINT4 *)LALMalloc(in->sizes[i] * sizeof(UINT4));
      if(!(out->f[out->nclusters-1])) {
	ABORT ( status,TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC);
      }

      out->P[out->nclusters-1] = (REAL8 *)LALMalloc(in->sizes[i] * sizeof(REAL8));
      if(!(out->P[out->nclusters-1])) {
	ABORT ( status,TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC);
      }

      for(j=0; j<in->sizes[i]; j++) {
	out->t[out->nclusters-1][j] = in->t[i][j];
	out->f[out->nclusters-1][j] = in->f[i][j];
	out->P[out->nclusters-1][j] = in->P[i][j];
      }
    }
  }
  

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/****************************UTILITY FUNCTIONS*****************************/
/******** <lalLaTeX file="TFClustersC"> ********
\newpage
\noindent
Merge two cluster lists.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void
LALMergeClusterLists (
		      LALStatus *status, 
		      CList *out,
		      CList *A, 
		      CList *B
		      )
}
\idx{LALMergeClusterLists()}

\subsubsection*{Description}
Merge \texttt{*A} and \texttt{*B} into cluster list \texttt{*out}. The merging is done so that any two clusters that overlapp or that have black pixels that are nearest neighbors will be replaced by the union of the two clusters in \texttt{*out}. The clusters that don't satisfy these two conditions are just copied into \texttt{*out}.

\subsubsection*{Uses}
\begin{verbatim}
LALCopyCList()
LALInitCList()
LALFreeCList()
\end{verbatim}

\subsubsection*{Notes}
\begin{itemize}
\item \texttt{*out} must be initialized by a proper call to \texttt{LALInitCList()} before calling this function.
\end{itemize}
\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> ********/
void
LALMergeClusterLists (
		      LALStatus *status, 
		      CList *out_,
		      CList *A, 
		      CList *B
		      )
{
  UINT4 i,j,k,l;
  UINT4 **bboxA, **bboxB;
  BOOLEAN *cA, *cB;
  BOOLEAN kStop, safeOutput=0;
  CList *out = NULL;

  INITSTATUS (status, "LALMergeClusterLists", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);

  ASSERT ( A, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( B, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( out_, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );

  if(out_ != A && out_ != B) {
    ASSERT ( out_->nclusters == 0, status, TFCLUSTERSH_ENZERO, TFCLUSTERSH_MSGENZERO);
    ASSERT ( out_->sizes == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
    ASSERT ( out_->t == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
    ASSERT ( out_->f == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
    ASSERT ( out_->P == NULL, status, TFCLUSTERSH_ENNULLP, TFCLUSTERSH_MSGENNULLP);
  }


  if(A->nclusters > 0) {
    ASSERT ( A->sizes, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( A->t, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( A->f, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( A->P, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( A->params, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  }
  else {
    if(B != out_) {
      LALCopyCList(status->statusPtr, out_, B);
      CHECKSTATUSPTR (status);
    }
    
    DETATCHSTATUSPTR (status);
    RETURN (status); 
  }

  if(B->nclusters > 0) {
    ASSERT ( B->sizes, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( B->t, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( B->f, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( B->P, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    ASSERT ( B->params, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  }
  else {
    
    if(A!=out_) {
      LALCopyCList(status->statusPtr, out_, A);
      CHECKSTATUSPTR (status);
    }
    
    DETATCHSTATUSPTR (status);
    RETURN (status); 
  }


  /* protect output */
  if(A == out_ || B == out_) {
    out = (CList *)LALMalloc(sizeof(CList));
    if(!out) {
      ABORT(status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );
    }

    LALInitCList(status->statusPtr, out, out_->params);
    CHECKSTATUSPTR (status);

    safeOutput = 1;
  }
  else {out = out_;}


  bboxA = (UINT4 **)LALMalloc(A->nclusters * sizeof(UINT4 *));
  bboxB = (UINT4 **)LALMalloc(B->nclusters * sizeof(UINT4 *));


  for(i=0;i<A->nclusters;i++)
    {
      bboxA[i] = (UINT4 *)LALMalloc(4*sizeof(UINT4));
      if(!bboxA[i])
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

      bboxA[i][0] = bboxA[i][1] = A->t[i][0];
      bboxA[i][2] = bboxA[i][3] = A->f[i][0];

      for(j=1;j<A->sizes[i];j++)
	{
	  if(A->t[i][j] < bboxA[i][0]) bboxA[i][0] = A->t[i][j];
	  if(A->t[i][j] > bboxA[i][1]) bboxA[i][1] = A->t[i][j];
	  if(A->f[i][j] < bboxA[i][2]) bboxA[i][2] = A->f[i][j];
	  if(A->f[i][j] > bboxA[i][3]) bboxA[i][3] = A->f[i][j];
	}
    }

  for(i=0;i<B->nclusters;i++)
    {
      bboxB[i] = (UINT4 *)LALMalloc(4*sizeof(UINT4));
      if(!bboxB[i])
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

      bboxB[i][0] = bboxB[i][1] = B->t[i][0];
      bboxB[i][2] = bboxB[i][3] = B->f[i][0];

      for(j=1;j<B->sizes[i];j++)
	{
	  if(B->t[i][j] < bboxB[i][0]) bboxB[i][0] = B->t[i][j];
	  if(B->t[i][j] > bboxB[i][1]) bboxB[i][1] = B->t[i][j];
	  if(B->f[i][j] < bboxB[i][2]) bboxB[i][2] = B->f[i][j];
	  if(B->f[i][j] > bboxB[i][3]) bboxB[i][3] = B->f[i][j];
	}
    }

  if(!(cA = (BOOLEAN *)LALMalloc(A->nclusters * sizeof(BOOLEAN))))
    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

  if(!(cB = (BOOLEAN *)LALMalloc(B->nclusters * sizeof(BOOLEAN))))
    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
  
  for(i=0; i<A->nclusters; i++) cA[i] = 1;
  for(i=0; i<B->nclusters; i++) cB[i] = 1;

  for(i=0; i<A->nclusters; i++) {
    if(cA[i]) {

      BOOLEAN firstOne = 1;

    for(j=0; j<B->nclusters; j++) { 
      if(cB[j]) {
	if(bboxA[i][0] <= bboxB[j][1] + 1 && 
	   bboxB[j][0] <= bboxA[i][1] + 1 &&
	   bboxA[i][2] <= bboxB[j][3] + 1 &&
	   bboxB[j][2] <= bboxA[i][3] + 1) {

	  for(k=0, kStop=0; k<A->sizes[i] && !kStop; k++)
	    for(l=0; l<B->sizes[j]; l++) {
	      if(abs(A->t[i][k] - B->t[j][l]) +
		 abs(A->f[i][k] - B->f[j][l]) <= 1) {
		kStop = 1;
		break;
	      }
	    }

	  if(kStop) { /* the two are really touching; merge */
	    if(firstOne) {

	    firstOne = 0;
	    cA[i] = cB[j] = 0;

	    out->nclusters++;	  
	    out->sizes = (UINT4 *)LALRealloc(out->sizes, out->nclusters*sizeof(UINT4));
	    if(!(out->sizes))
	      {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	    out->f = (UINT4 **)LALRealloc(out->f, out->nclusters*sizeof(UINT4 *));
	    if(!(out->f))
	      {ABORT (  status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	    out->t = (UINT4 **)LALRealloc(out->t, out->nclusters*sizeof(UINT4 *));
	    if(!(out->t))
	      {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	    out->P = (REAL8 **)LALRealloc(out->P, out->nclusters*sizeof(REAL8 *));
	    if(!(out->P))
	      {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

	    out->sizes[out->nclusters-1] = A->sizes[i];

	    out->f[out->nclusters-1] = (UINT4 *)LALMalloc(out->sizes[out->nclusters-1] * sizeof(UINT4));
	    if(!(out->f[out->nclusters-1]))
	      {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	    out->t[out->nclusters-1] = (UINT4 *)LALMalloc(out->sizes[out->nclusters-1] * sizeof(UINT4));
	    if(!(out->t[out->nclusters-1]))
	      {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	    out->P[out->nclusters-1] = (REAL8 *)LALMalloc(out->sizes[out->nclusters-1] * sizeof(REAL8));
	    if(!(out->P[out->nclusters-1]))
	      {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	    for(k=0;k<A->sizes[i];k++)
	      {out->f[out->nclusters-1][k] = A->f[i][k];
	      out->t[out->nclusters-1][k] = A->t[i][k];
	      out->P[out->nclusters-1][k] = A->P[i][k];}

	    
	    for(k=0;k<B->sizes[j]; k++) {
	      for(l=0; l<A->sizes[i]; l++) {
		if(A->t[i][l] == B->t[j][k] &&
		   A->f[i][l] == B->f[j][k]) {break;}
	      }  
	      if(l == A->sizes[i])
		{

		  (out->sizes[out->nclusters-1])++;

		  out->f[out->nclusters-1] = (UINT4 *)LALRealloc(out->f[out->nclusters-1], out->sizes[out->nclusters-1] * sizeof(UINT4));
		  if(!(out->f[out->nclusters-1]))
		    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
		  out->t[out->nclusters-1] = (UINT4 *)LALRealloc(out->t[out->nclusters-1], out->sizes[out->nclusters-1] * sizeof(UINT4));
		  if(!(out->t[out->nclusters-1]))
		    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
		  out->P[out->nclusters-1] = (REAL8 *)LALRealloc(out->P[out->nclusters-1], out->sizes[out->nclusters-1] * sizeof(REAL8));
		  if(!(out->P[out->nclusters-1]))
		    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

		  out->f[out->nclusters-1][out->sizes[out->nclusters-1]-1] = B->f[j][k];
		  out->t[out->nclusters-1][out->sizes[out->nclusters-1]-1] = B->t[j][k];
		 out->P[out->nclusters-1][out->sizes[out->nclusters-1]-1] = B->P[j][k];
		}
	    }
	    }
	    else /* not firstOne */ 
	      {
		UINT4 ktmp = out->sizes[out->nclusters-1];
		for(l=0; l<B->sizes[j]; l++) {
		  for(k=0;k<ktmp; k++) {

		    if(B->t[j][l] == out->t[out->nclusters-1][k] &&
		       B->f[j][l] == out->f[out->nclusters-1][k]) {break;}
		  }  

		  if(k == ktmp)
		    {

		      (out->sizes[out->nclusters-1])++;

		      out->f[out->nclusters-1] = (UINT4 *)LALRealloc(out->f[out->nclusters-1], out->sizes[out->nclusters-1] * sizeof(UINT4));
		      if(!(out->f[out->nclusters-1]))
			{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
		      out->t[out->nclusters-1] = (UINT4 *)LALRealloc(out->t[out->nclusters-1], out->sizes[out->nclusters-1] * sizeof(UINT4));
		      if(!(out->t[out->nclusters-1]))
			{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
		      out->P[out->nclusters-1] = (REAL8 *)LALRealloc(out->P[out->nclusters-1], out->sizes[out->nclusters-1] * sizeof(REAL8));
		      if(!(out->P[out->nclusters-1]))
			{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
		      
		      out->f[out->nclusters-1][out->sizes[out->nclusters-1]-1] = B->f[j][l];
		      out->t[out->nclusters-1][out->sizes[out->nclusters-1]-1] = B->t[j][l];
		      out->P[out->nclusters-1][out->sizes[out->nclusters-1]-1] = B->P[j][l];
		    }
		}

	      }
	  }
	}
      
	if(cA[i]) { /* doesn't touch to anybody */
	
	  cA[i] = 0;

	  out->nclusters++;	  
	  out->sizes = (UINT4 *)LALRealloc(out->sizes, out->nclusters*sizeof(UINT4));
	  if(!(out->sizes))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	  out->sizes[out->nclusters-1] = A->sizes[i];
	  out->f = (UINT4 **)LALRealloc(out->f, out->nclusters*sizeof(UINT4 *));
	  if(!(out->f))
	    {ABORT (  status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	  out->t = (UINT4 **)LALRealloc(out->t, out->nclusters*sizeof(UINT4 *));
	  if(!(out->t))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	  out->P = (REAL8 **)LALRealloc(out->P, out->nclusters*sizeof(REAL8 *));
	  if(!(out->P))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	  out->f[out->nclusters-1] = (UINT4 *)LALMalloc(out->sizes[out->nclusters-1] * sizeof(UINT4));
	  if(!(out->f[out->nclusters-1]))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	  out->t[out->nclusters-1] = (UINT4 *)LALMalloc(out->sizes[out->nclusters-1] * sizeof(UINT4));
	  if(!(out->t[out->nclusters-1]))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	  out->P[out->nclusters-1] = (REAL8 *)LALMalloc(out->sizes[out->nclusters-1] * sizeof(REAL8));
	  if(!(out->P[out->nclusters-1]))
	    {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
	  for(k=0;k<A->sizes[i];k++)
	    {out->f[out->nclusters-1][k] = A->f[i][k];
	    out->t[out->nclusters-1][k] = A->t[i][k];
	    out->P[out->nclusters-1][k] = A->P[i][k];}
	}
      }
    }
  }
  }
	

  for(i=0;i<B->nclusters; i++) { 
    if(cB[i]) {

      cB[i] = 0;

      out->nclusters++;	  
      out->sizes = (UINT4 *)LALRealloc(out->sizes, out->nclusters*sizeof(UINT4));
      if(!(out->sizes))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
      out->sizes[out->nclusters-1] = B->sizes[i];
      out->f = (UINT4 **)LALRealloc(out->f, out->nclusters*sizeof(UINT4 *));
      if(!(out->f))
	{ABORT (  status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
      out->t = (UINT4 **)LALRealloc(out->t, out->nclusters*sizeof(UINT4 *));
      if(!(out->t))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
      out->P = (REAL8 **)LALRealloc(out->P, out->nclusters*sizeof(REAL8 *));
      if(!(out->P))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
      out->f[out->nclusters-1] = (UINT4 *)LALMalloc(out->sizes[out->nclusters-1] * sizeof(UINT4));
      if(!(out->f[out->nclusters-1]))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
      out->t[out->nclusters-1] = (UINT4 *)LALMalloc(out->sizes[out->nclusters-1] * sizeof(UINT4));
      if(!(out->t[out->nclusters-1]))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
      out->P[out->nclusters-1] = (REAL8 *)LALMalloc(out->sizes[out->nclusters-1] * sizeof(REAL8));
      if(!(out->P[out->nclusters-1]))
	{ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
      for(k=0;k<B->sizes[i];k++)
	{out->f[out->nclusters-1][k] = B->f[i][k];
	out->t[out->nclusters-1][k] = B->t[i][k];
	out->P[out->nclusters-1][k] = B->P[i][k];}
    }
  }

  LALFree(cA);
  LALFree(cB);

  for(i=0;i<A->nclusters; i++) LALFree(bboxA[i]);
  for(i=0;i<B->nclusters; i++) LALFree(bboxB[i]);
  
  LALFree(bboxA);
  LALFree(bboxB);

  if(safeOutput) {
    LALCopyCList(status->statusPtr, out_, out);
    CHECKSTATUSPTR (status);

    LALFreeCList(status->statusPtr, out);
    CHECKSTATUSPTR (status);

    LALFree(out);
  }

  DETATCHSTATUSPTR (status);
  RETURN (status); 
}


/******** <lalLaTeX file="TFClustersC"> ********
\newpage
\noindent
Make a copy of a cluster list.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void
LALCopyCList (
	      LALStatus *status, 
	      CList *dest,
	      CList *src
	      );
}
\idx{LALCopyCList()}

\subsubsection*{Description}
Make a copy of \texttt{*src} onto \texttt{*dest}.

\subsubsection*{Uses}
\begin{verbatim}
LALFreeCList()
\end{verbatim}

\subsubsection*{Notes}
\begin{itemize}
\item \texttt{*src}, if not empty, will be overwritten.
\item \texttt{*dest} must be initialized by a proper call to \texttt{LALInitCList()} before calling this function.
\end{itemize}
\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> ********/

void
LALCopyCList (
	      LALStatus *status, 
	      CList *dest,
	      CList *src
	      )
{
  UINT4 i;

  INITSTATUS (status, "LALCopyCList", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);

  ASSERT ( dest, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( src, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );

  LALFreeCList(status->statusPtr, dest);
  CHECKSTATUSPTR (status);
  
  dest->nclusters = src->nclusters;

  if(src->nclusters > 0) {

    dest->sizes = (UINT4*)LALMalloc(src->nclusters * sizeof(UINT4));
    if(!dest->sizes) {
      ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );
    }
    memcpy(dest->sizes, src->sizes, src->nclusters * sizeof(UINT4));

  
    dest->t = (UINT4 **)LALMalloc(src->nclusters * sizeof(UINT4 *));
    if(!dest->t) {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
    dest->f = (UINT4 **)LALMalloc(src->nclusters * sizeof(UINT4 *));
    if(!dest->f) {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}
    dest->P = (REAL8 **)LALMalloc(src->nclusters * sizeof(REAL8 *));
    if(!dest->P) {ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );}

    for(i=0; i<src->nclusters; i++) {
      dest->t[i] = (UINT4 *)LALMalloc(src->sizes[i] * sizeof(UINT4));
      dest->f[i] = (UINT4 *)LALMalloc(src->sizes[i] * sizeof(UINT4));
      dest->P[i] = (REAL8 *)LALMalloc(src->sizes[i] * sizeof(REAL8));
      
      memcpy(dest->t[i], src->t[i], src->sizes[i] * sizeof(UINT4));
      memcpy(dest->f[i], src->f[i], src->sizes[i] * sizeof(UINT4));
      memcpy(dest->P[i], src->P[i], src->sizes[i] * sizeof(REAL8));
    }
  }

  dest->params = (TFPlaneParams *)LALMalloc(sizeof(TFPlaneParams));
  memcpy(dest->params, src->params, sizeof(TFPlaneParams));

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/******** <lalLaTeX file="TFClustersC"> ********
\newpage
\noindent
Initialize a spectrogram with default values.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void 
LALPlainSpectrogram(
		    LALStatus *status,
		    TFPlaneParams *tspec,
		    REAL4TimeSeries *tseries,
		    REAL8 T
		    )
}
\idx{LALPlainSpectrogram()}

\subsubsection*{Description}
Initialize the spectrogram \texttt{*tspec} so that it has a time resolution \texttt{T} and frequency resolution 1/\texttt{T}, with frequency ranging from 1/\texttt{T} to the Nyquist frequency of the time series \texttt{*tseries}. Also set the length of \texttt{*tspec} so it matches \texttt{*tseries}.

\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> ********/
void 
LALPlainSpectrogram(
		    LALStatus *status,
		    TFPlaneParams *tspec,
		    REAL4TimeSeries *tseries,
		    REAL8 T
		    )
{
  INITSTATUS (status, "LALPlainSpectrogram", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);


  ASSERT ( tseries, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( tspec, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( T > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS );


  tspec->timeBins = (INT4)floor((double)tseries->data->length * tseries->deltaT / T);
  tspec->freqBins = (INT4)ceil(T/(2.0*tseries->deltaT)) - 1; /* do not include Nyquist */
  tspec->deltaT = T;
  tspec->flow = 1.0/T; /* do not include DC */

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/******** <lalLaTeX file="TFClustersC"> ********
\newpage
\noindent
Initialize a spectrogram with default values.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void 
LALPlainSpectrogramWin(
		    LALStatus *status,
		    TFPlaneParams *tspec,
		    REAL4TimeSeries *tseries,
		    REAL8 T
		    )
}
\idx{LALPlainSpectrogramWin()}

\subsubsection*{Description}
Initialize the spectrogram \texttt{*tspec} so that it has a time resolution \texttt{T}/2 and frequency resolution 1/\texttt{T}, with frequency ranging from 1/\texttt{T} to the Nyquist frequency of the time series \texttt{*tseries}. 

\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> ********/
void 
LALPlainSpectrogramWin(
		    LALStatus *status,
		    TFPlaneParams *tspec,
		    REAL4TimeSeries *tseries,
		    REAL8 T
		    )
{
  INITSTATUS (status, "LALPlainSpectrogram", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);


  ASSERT ( tseries, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( tspec, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( T > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS );


  tspec->timeBins = 2*(INT4)floor((double)tseries->data->length * tseries->deltaT / T) - 1;
  tspec->freqBins = (INT4)ceil(T/(2.0*tseries->deltaT)) - 1; /* do not include Nyquist */
  tspec->deltaT = T;
  tspec->flow = 1.0/T; /* do not include DC */

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/******** <lalLaTeX file="TFClustersC"> ********
\newpage
\noindent
Initialize a cluster list structure.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void 
LALInitCList (
	      LALStatus *status,
	      CList *clist,
	      TFPlaneParams *tspec
	      )
}
\idx{LALInitCList()}

\subsubsection*{Description}
Initialize \texttt{*clist} and set its parameters to \texttt{*tspec}.

\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> ********/
void 
LALInitCList (
	      LALStatus *status,
	      CList *clist,
	      TFPlaneParams *tspec
	      )
{
  INITSTATUS (status, "LALInitCList", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);

  ASSERT ( clist, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( tspec, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( tspec->timeBins > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  ASSERT ( tspec->freqBins > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  ASSERT ( tspec->deltaT > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  ASSERT ( tspec->flow >= 0.0, status, TFCLUSTERSH_EPOS, TFCLUSTERSH_MSGEPOS);


  clist->nclusters = 0;
  clist->t = clist->f = NULL;
  clist->sizes = NULL;
  clist->P = NULL;
  clist->params = NULL;

  /* copy tspec */
  clist->params = (TFPlaneParams *)LALMalloc(sizeof(TFPlaneParams));
  if(!(clist->params)) {
    ABORT ( status, TFCLUSTERSH_EMALLOC, TFCLUSTERSH_MSGEMALLOC );
  }

  clist->params->timeBins = tspec->timeBins;
  clist->params->freqBins = tspec->freqBins;
  clist->params->flow = tspec->flow; 
  clist->params->deltaT = tspec->deltaT;


  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}




/******** <lalLaTeX file="TFClustersC"> ********
\newpage
\noindent
Initialize a threshold structure.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void 
LALFillCListDir (
		 LALStatus *status,
		 CListDir *cldir,
		 REAL8 rho
		 )
}
\idx{LALFillCListDir()}

\subsubsection*{Description}
Initialize \texttt{*cldir}. This means allocating memory for \texttt{cldir->s1}, \texttt{cldir->s2}, \texttt{cldir->d} according to the value of \texttt{cldir->sigma}, and for the threshold vector \texttt{cldir->rho}; all the values of \texttt{cldir->rho} are initialized to \texttt{rho}.

\subsubsection*{Notes}
\begin{itemize}
\item Before calling this function, \texttt{cldir->sigma} and \texttt{cldir->freqBins} must be set to there desired value.
\end{itemize}

\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> ********/
void 
LALFillCListDir (
		 LALStatus *status,
		 CListDir *cldir,
		 REAL8 rho
		 )
{
  UINT4 i,j,k;


  INITSTATUS (status, "LALFillCListDir", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);

  ASSERT ( cldir, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
  ASSERT ( cldir->freqBins > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  /*
  ASSERT ( cldir->maxf > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  */
  ASSERT ( cldir->sigma > 0, status, TFCLUSTERSH_ESTRICTPOS, TFCLUSTERSH_MSGESTRICTPOS);
  ASSERT ( rho >= 0, status, TFCLUSTERSH_EPOS, TFCLUSTERSH_MSGEPOS);

  if(cldir->sigma > 1) {

    cldir->s1 = (UINT4 *)LALMalloc(cldir->sigma *(cldir->sigma-1) * sizeof(UINT4)/2);
    if(!(cldir->s1)) {
      ABORT ( status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );
    }

    cldir->s2 = (UINT4 *)LALMalloc(cldir->sigma *(cldir->sigma-1) * sizeof(UINT4)/2);
    if(!(cldir->s2)) {
      ABORT ( status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP ); 
    }

    cldir->d = (UINT4 *)LALMalloc(cldir->sigma *(cldir->sigma-1) * sizeof(UINT4)/2);
    if(!(cldir->d)) {
      ABORT ( status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP ); 
    }


    for(k=0, i=1; i < cldir->sigma; i++)
      for(j=i;j<cldir->sigma;j++)
	{cldir->s1[k] = i;
	cldir->s2[k] = j;
	k++;}
  }
  else {
    cldir->s1 = cldir->s2 = NULL;
    cldir->d = NULL;
  }

  cldir->rho = (REAL8 *)LALMalloc(cldir->freqBins * sizeof(REAL8));
  if(!(cldir->rho)) {
    ABORT ( status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP ); 
  }

  for(i=0; i<cldir->freqBins; i++)
    cldir->rho[i] = rho;

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/******** <lalLaTeX file="TFClustersC"> ********
\newpage
\noindent
Functions to destroy the different structures.
\subsubsection*{Prototype}
\vspace{0.1in}
\texttt{
void 
LALFreeCList(
	     LALStatus *status, 
	     CList *clist
	     )
}
\idx{LALFreeCList()}

\vspace{0.1in}
\texttt{
\noindent
void 
LALFreeSpecgram(
		LALStatus *status, 
		Spectrogram *power
		)
}
\idx{LALFreeSpecgram()}

\vspace{0.1in}
\texttt{
\noindent
void
LALFreeCListDir (
		 LALStatus *status,
		 CListDir *cdir
		 )
}
\idx{LALFreeCListDir()}

\subsubsection*{Description}
Release allocated memory.

\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> ********/

void 
LALFreeCList(
	     LALStatus *status, 
	     CList *clist
	     )
{
  UINT4 i;

  INITSTATUS (status, "LALFreeCList", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);

  ASSERT ( clist, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );

  if(clist->nclusters > 0) {

    for(i=0; i<clist->nclusters; i++) {
      LALFree(clist->t[i]);
      LALFree(clist->f[i]);
      LALFree(clist->P[i]);
    }

    LALFree(clist->sizes);
    LALFree(clist->t);
    LALFree(clist->f);
    LALFree(clist->P);
  }
    
  LALFree(clist->params);

  clist->nclusters = 0;
  clist->t = clist->f = NULL;
  clist->P = NULL;
  clist->params = NULL;

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



void 
LALFreeSpecgram(
		LALStatus *status, 
		Spectrogram *power
		)
{
  INITSTATUS (status, "LALFreeSpecgram", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);

  ASSERT ( power, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );

  LALFree(power->power);
  LALFree(power->params);

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



void
LALFreeCListDir (
		 LALStatus *status,
		 CListDir *cdir
		 )
{
  INITSTATUS (status, "LALFreeCListDir", TFCLUSTERSC);
  ATTATCHSTATUSPTR (status);

  ASSERT ( cdir, status, TFCLUSTERSH_ENULLP, TFCLUSTERSH_MSGENULLP );

  if(cdir->sigma > 1) {
    LALFree(cdir->s1);
    LALFree(cdir->s2);
    LALFree(cdir->d);
  }

  LALFree(cdir->rho);

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}



/**************************INTERNAL FUNCTIONS***************************/
static void GetNearestNeighb(CList *clist, Spectrogram *tpower, REAL8 *rho, UINT4 i, UINT4 j) 
{
  /* first, add point to clist */
  (clist->sizes[clist->nclusters - 1])++;

  clist->t[clist->nclusters - 1] = (UINT4 *)LALRealloc( clist->t[clist->nclusters - 1], clist->sizes[clist->nclusters - 1] * sizeof(UINT4));
  clist->f[clist->nclusters - 1] = (UINT4 *)LALRealloc( clist->f[clist->nclusters - 1], clist->sizes[clist->nclusters - 1] * sizeof(UINT4));
  clist->P[clist->nclusters - 1] = (REAL8 *)LALRealloc( clist->P[clist->nclusters - 1], clist->sizes[clist->nclusters - 1] * sizeof(REAL8));


  clist->t[clist->nclusters - 1][clist->sizes[clist->nclusters - 1] - 1] = i;
  clist->f[clist->nclusters - 1][clist->sizes[clist->nclusters - 1] - 1] = j;
  clist->P[clist->nclusters - 1][clist->sizes[clist->nclusters - 1] - 1] = tpower->power[i*tpower->params->freqBins + j];
  tpower->power[i*tpower->params->freqBins + j] = 0;

  /* next, check for neighbors */

  if(i >= 1)
    if(tpower->power[(i-1)*tpower->params->freqBins + j] > rho[j])
      GetNearestNeighb(clist, tpower, rho, i-1, j);

  if(i+1 < (UINT4)tpower->params->timeBins)
    if(tpower->power[(i+1)*tpower->params->freqBins + j] > rho[j])
      GetNearestNeighb(clist, tpower, rho, i+1, j);

  if(j >= 1)
    if(tpower->power[i*tpower->params->freqBins + j - 1] > rho[j-1])
      GetNearestNeighb(clist, tpower, rho, i, j - 1);

  if(j+1 < (UINT4)tpower->params->freqBins)
    if(tpower->power[i*tpower->params->freqBins + j + 1] > rho[j+1])
      GetNearestNeighb(clist, tpower, rho, i, j + 1);
}




static UINT4 ClustDistance(UINT4 s1, UINT4 *t1, UINT4 *f1, UINT4 s2, UINT4 *t2, UINT4 *f2)
{
  UINT4 d = -1U, i, j, dtmp;

  for(i=0;i<s1;i++)
    for(j=0;j<s2;j++)
      {dtmp = abs((INT4)t1[i] - (INT4)t2[j]) + abs((INT4)f1[i] - (INT4)f2[j]);
      if(dtmp < d)
	d = dtmp;}

  return d;
}



static UINT4 min(UINT4 a, UINT4 b)
{
  if(a<=b) return a;
  return b;
}


static UINT4 max(UINT4 a, UINT4 b)
{
  if(a>b) return a;
  return b;
}


#define ITMAX 1000

/* the following returns
1 - int_0^x exp(-t) t^(a-1) dt / Gamma(a)
*/
static void incgam(LALStatus *status, REAL4 a, REAL4 x, REAL4 *retu)
{
  INT4 n;
  REAL4 sum,del,elln,an,a1,a2,a3;
  REAL8 tmp,ser;
  static REAL8 cof[6]={76.18009172947146,
		       -86.50532032941677,
		       24.01409824083091,
		       -1.231739572450155,
		       0.1208650973866179e-2,
		       -0.5395239384953e-5};

  *retu = 0;

  if(x<0.0 || a <= 0.0) {
    ABORT ( status, TFCLUSTERSH_EIARG, TFCLUSTERSH_MSGEIARG);
  }

  tmp = a + 5.5 - (a + 0.5) * log(a + 5.5);
  ser=1.000000000190015;
  for (n=1;n<=6;n++) ser += cof[n-1]/(a + n);
  elln = -tmp+log(2.5066282746310005*ser/a);

  elln = exp(-x+a*log(x)-elln);


  if (x - 1.0 < a) {
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      del *= x/(a+n);
      sum += del;
      if (fabs(del) < 3e-7 * fabs(sum)) {
	elln *= sum;
	break;
      }
    }
   
    if(n>ITMAX) 
      {ABORT ( status, TFCLUSTERSH_EMAXITE, TFCLUSTERSH_MSGEMAXITE);}

    *retu = 1.0 - elln;
    return;
  } 
  else {
    a1=x+1.0-a;
    a2=1.0/1e-30;
    sum=a3=1.0/a1;

    for (n=1;n<=ITMAX;n++) {
      an = -n*(n-a);
      a1+=2.0;

      a3=an*a3+a1;
      if (fabs(a3) < 1e-30) a3=1e-30;
      
      a2=a1+an/a2;
      if (fabs(a2) < 1e-30) a2=1e-30;
      
      a3=1.0/a3;
      del=a2*a3;
      sum *= del;
      if (fabs(del-1.0) < 3e-7) break;
    }
    
    if(n>ITMAX) {
      ABORT ( status, TFCLUSTERSH_EMAXITE, TFCLUSTERSH_MSGEMAXITE);
    }


    elln *= sum;
    *retu = elln;
    return;
  }
}

#undef ITMAX


/******** <lalLaTeX file="TFClustersC"> ********
\newpage
\subsubsection{Operating Instructions}
\texttt{
  \#include <TFClusters.h> \\
  \#include <Random.h> \\
\\
  static LALStatus status;\\
  REAL4TimeSeries tseries;\\
  CListDir dir;\\
  CList clist, list;\\
  TFPlaneParams tspec;\\
  Spectrogram spower;\\
\\
  RandomParams *params;\\
  REAL4Vector *vect;\\
\\
  REAL8 T, P;\\
  UINT4 i, j, N;\\
  INT4 seed = 0;\\
\\
\\
  {\rm\bf\it first generate a time series of white Gaussian noise of unit variance}\\
  N = 1024 * 128;\\
  LALCreateVector(\&status, \&vect, N);\\
  LALCreateRandomParams(\&status, \&params, seed);\\
  LALNormalDeviates(\&status, vect, params);\\
\\
  tseries.epoch.gpsSeconds = 0;\\
  tseries.epoch.gpsNanoSeconds = 0;\\
  tseries.deltaT = 1.0/1024.0;\\
  tseries.f0 = 0.0;\\
  tseries.data = vect;\\
\\
\\
  {\rm\bf\it Next compute a spectrogram for the time series }\\
  T = 1.0; {\rm\bf\it this is the resolution in seconds of the spectrogram }\\
\\
  LALPlainSpectrogram(\&status, \&tspec, \&tseries, T); {\rm\bf\it this creates spectrogram parameters at the 'Heisenberg limit' from DC+1/T to the Nyquist frequency }\\
  LALComputeSpectrogram(\&status, \&spower, \&tspec, \&tseries);\\
\\
\\
  {\rm\bf\it Set thresholds }\\
  dir.freqBins = tspec.freqBins; {\rm\bf\it number of frequency bins in spectrogram }\\
  dir.sigma = 5; {\rm\bf\it threshold on cluster size }\\
  dir.minf = 0; {\rm\bf\it min frequency to consider (Hz) } \\
  dir.maxf = 512; {\rm\bf\it max frequency to consider (Hz) }\\
\\
  LALFillCListDir(\&status, \&dir, -log(0.1)); {\rm\bf\it allocate memory and set the threshold on power so that 1 every 10 pixel is black }\\
\\
\\
  {\rm\bf\it set thresholds on distance for different size couples }\\
  dir.d[0] = 0; {\rm\bf\it 1,1 }\\
  dir.d[1] = 0; {\rm\bf\it ... }\\
  dir.d[2] = 0;\\
  dir.d[3] = 0; {\rm\bf\it 1,4 }\\
  dir.d[4] = 0; {\rm\bf\it 2,2 }\\
  dir.d[5] = 0; \\
  dir.d[6] = 2; {\rm\bf\it 2,4 }\\
  dir.d[7] = 3; {\rm\bf\it 3,3 }\\
  dir.d[8] = 4; {\rm\bf\it 3,4 }\\
  dir.d[9] = 4; {\rm\bf\it 4,4 }\\
\\
  dir.mdist = 4; {\rm\bf\it no need to worry about things that are more than 4 units away from each other }\\
\\
  {\rm\bf\it run cluster threshold algorithm }\\
  LALInitCList(\&status, \&clist, \&tspec); {\rm\bf\it initialize clist }\\
\\
  LALGetClusters(\&status, \&clist, \&spower, \&dir); {\rm\bf\it generate list of clusters }\\
\\
  LALFreeSpecgram(\&status, \&spower); {\rm\bf\it spectrogram no longer useful }\\
\\
\\
  {\rm\bf\it run threshold on cluster total power }\\
  dir.alpha = 0.25; {\rm\bf\it only 1/4 of all clusters from white noise will make it }\\
\\
  LALInitCList(\&status, \&list, \&tspec); {\rm\bf\it initialize list }\\
  \\
  LALClustersPowerThreshold(\&status, \&list, \&clist, \&dir); {\rm\bf\it generate new list }\\
\\
\\
  {\rm\bf\it clean up a bit }\\
  LALFreeCList(\&status, \&clist);\\
  LALFreeCListDir(\&status, \&dir);\\
\\
\\
  {\rm\bf\it display results to stdout }\\
  printf("Id/t/tSize/t/tPower/n");\\
  for(i=0; i<list.nclusters; i++) {\\
    for(P=0, j=0; j<list.sizes[i]; j++) P += list.P[i][j];\\
    printf("\%i/t/t\%i/t/t\%g/n",i,list.sizes[i],P);\\
  }\\
  \\
  {\rm\bf\it clean up }\\
  LALFreeCList(\&status, \&list);\\
}

\vfill{\footnotesize\input{TFClustersCV}}
********* </lalLaTeX> *******/
