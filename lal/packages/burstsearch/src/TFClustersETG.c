/*-----------------------------------------------------------------------*
 *
 * File Name: TFClustersETG.c
 *
 * Author: Julien Sylvestre
 *
 * Revision: $Id$ 
 *
 *-----------------------------------------------------------------------*/



/******** <lalVerbatim file="TFClustersETGCV"> ********
Author: Sylvestre, J
$Id$
********* </lalVerbatim> ********/

#include <config.h>
#include <lal/StdBurstSearch.h>
#include <lal/Sort.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/TFClusters.h>
#include <lal/AVFactories.h>
#include <lal/TFCThresholds.h>
#include <lal/LALRCSID.h>
#include <string.h>

NRCSID (TFCLUSTERSETGC, "$Id$");

#ifndef WORDS_BIGENDIAN
  static void endian_swap(char * pdata, int dsize, int nelements);
#endif

/******** <lalLaTeX file="TFClustersETGC"> ********
\noindent
Implement the TFCLUSTERS event trigger generator.
\subsubsection*{Prototype}
********* </lalLaTeX> ********/
/* <lalVerbatim> */
void
LALTFClustersETG(
		 LALStatus *status, 
		 EventIDColumn *output, 
		 REAL4TimeVectorSeries *input, 
		 BurstParameter *params
	       ) {
/* </lalVerbatim> */
/******** <lalLaTeX file="TFClustersETGC"> ********
\subsubsection*{Description}
Description of the parameters: 
\begin{center}
\begin{tabular}{l|l|l}
parameter index & type & description \\ \hline 
1 & string & channel name \\
2 & REAL4 & black pixel probability \\
3 & INT4 & 1 (0) for (no) windowing \\
4 & INT4 & thresholds method: \\
  &      & 0: rank method \\
  &      & 1: white noise \\
  &      & 2: Rice fit \\
5 & INT4 & 1 (0) to (not) save tf data \\
6 & REAL4 & time resolution (s) \\
7 & REAL4 & min frequency (Hz) \\
8 & REAL4 & max frequency (Hz) \\
9 & REAL4 & alpha \\
10 & INT4 & sigma \\
11 & INT4 & delta(1,1) \\
... & ... & ... \\
10 + sigma*(sigma-1)/2 & INT4 & delta(sigma-1,sigma-1) 
\end{tabular}
\end{center}

\subsubsection*{Uses}
\begin{verbatim}
...a bunch of stuff.
\end{verbatim}
********* </lalLaTeX> ********/


  INT4 thr_method;
  UINT4 i, j, delL, minF, ip;
  REAL4 p,T,fs;
  REAL4 mint, maxt, minf, maxf, trez;
  REAL8Vector *tmpv;
  TFPlaneParams tspec;
  TFCSpectrogram spower;
  CList clist, list;
  CListDir dir;
  CHAR ifo[LIGOMETA_IFO_MAX];
  CHAR channel[LIGOMETA_CHANNEL_MAX];
  CHAR  search[LIGOMETA_SEARCH_MAX] = "STD_TFCLUSTERS";
  BOOLEAN blob=0, win=0;

  INITSTATUS (status, "LALBurstOutput", TFCLUSTERSETGC);
  ATTATCHSTATUSPTR (status);

  /* trivial input check */
  ASSERT ( output, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  ASSERT ( input, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  ASSERT ( params, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);

  /* parse parameters */
#define SetParameter(par,type) if(!(params)) {ABORT(status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI);} \
  if(!(params->type)) {ABORT(status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI);} \
  par = *(params->type); \
  params = params->next

#define SetStringParameter(par) if(!(params)) {ABORT(status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI);} \
  if(!(params->char_)) {ABORT(status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI);} \
  ASSERT ( par, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP); \
  strcpy(par,params->char_); \
  params = params->next

  params = params->next;
  
  SetStringParameter(channel);
  strncpy(ifo,channel,2);

  SetParameter(p,real4_);
  SetParameter(win,int4_);
  SetParameter(thr_method,int4_);
  SetParameter(blob,int4_);

  SetParameter(T,real4_);

  SetParameter(dir.minf,real4_);
  SetParameter(dir.maxf,real4_);
  SetParameter(dir.alpha,real4_);
  SetParameter(dir.sigma,int4_);

  fs = 1.0 / input->deltaT;
  dir.freqBins = 1 + (UINT4)(fs*T/2.0);

  LALFillCListDir(status->statusPtr, &dir, -log(p)); /* note: allocates memory for dir.rho from DC to Nyquist, and set to -log(p) */
  CHECKSTATUSPTR (status);

  delL = dir.sigma * (dir.sigma - 1) / 2;
  for(i=0;i<delL;i++) {
    SetParameter(dir.d[i],int4_);
  }

#undef SetParameter
#undef SetStringParameter

  dir.mdist = dir.d[0];
  for(i=1; i<delL; i++) {
    if(dir.d[i] > dir.mdist) {dir.mdist = dir.d[i];}
  }

  /***********TFCLUSTERS code************/
  /* construct spectrogram */
  {
    REAL4TimeSeries tmp;
    REAL4Vector tv;
    tmp.epoch = input->epoch;
    tmp.deltaT = input->deltaT;
    tmp.f0 = input->f0;
    tmp.data = &tv;
    tv.length = input->data->vectorLength;
    tv.data = input->data->data;
    
    if(win) {
      LALPlainTFCSpectrogramWin(status->statusPtr, &tspec, &tmp, T);
      CHECKSTATUSPTR (status);
      trez = T / 2.0;
    } else {
      LALPlainTFCSpectrogram(status->statusPtr, &tspec, &tmp, T);
      CHECKSTATUSPTR (status);
      trez = T;
    }

    spower.power = NULL;
    spower.params = NULL;

    if(input->data->length == 1) {
      LALComputeTFCSpectrogram(status->statusPtr, &spower, &tspec, &tmp);
      CHECKSTATUSPTR (status);
    } else {
      LALComputeXTFCSpectrogram(status->statusPtr, &spower, &tspec, input);
      CHECKSTATUSPTR (status);
    }
  }

  /* Estimate thresholds */
  switch(thr_method) {
  case 0:
    minF = (UINT4)floor(tspec.flow / fs);
    tmpv = NULL;
    LALDCreateVector(status->statusPtr, &tmpv, spower.params->timeBins);
    CHECKSTATUSPTR (status);
    ip = (UINT4)floor((1.0-p) * (spower.params->timeBins - 1));
    for(i=0; (int)i<(int)spower.params->freqBins; i++) {

      for(j=0; (int)j<(int)spower.params->timeBins; j++) {
	tmpv->data[j] = spower.power[j*spower.params->freqBins + i];
      }

      LALDHeapSort(status->statusPtr, tmpv);
      CHECKSTATUSPTR (status);
      dir.rho[i] = tmpv->data[ip];
    }
    LALDDestroyVector(status->statusPtr, &tmpv);
    CHECKSTATUSPTR (status);
    break;
  case 1:
    /* white noise of unit variance; nothing to do */
    break;
  case 2:
    /* rice fit */
    {
      UINT4 NP;
      REAL8 P, P2, R;
      RiceThresholdParams rtp;

      REAL8 Nsigma = 9.5; /* hardcoded !! */

      rtp.bpp = p;
      rtp.eGoal = 1e-3 * p; 
      rtp.nFreq = spower.params->freqBins;

      rtp.P0 = (REAL8 *)LALCalloc(rtp.nFreq, sizeof(REAL8));
      if(!(rtp.P0)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);} 

      rtp.Q = (REAL8 *)LALCalloc(rtp.nFreq, sizeof(REAL8));
      if(!(rtp.Q)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);} 

      for(i=0; (int)i<(int)spower.params->freqBins; i++) {

	P = P2 = 0.0;

	for(j=0; (int)j<(int)spower.params->timeBins; j++) {
	  P += spower.power[j*spower.params->freqBins + i];
	  P2 += spower.power[j*spower.params->freqBins + i] * spower.power[j*spower.params->freqBins + i];
	}

	P /= (REAL8)spower.params->timeBins;
	P2 /= (REAL8)spower.params->timeBins;

	R = 2.0*P*P - P2;
	    
	if(R>0.0) {
	  R = sqrt(R);
	} else {
	  R = 0.0;
	}

	rtp.Q[i] = R;
	rtp.P0[i] = P - R;

	/* rerun, computing mean & variance without last NSigma */
	P = P2 = 0.0;
	NP = 0;

	for(j=0; (int)j<(int)spower.params->timeBins; j++) {
	  if(spower.power[j*spower.params->freqBins + i] < rtp.Q[j] + Nsigma * rtp.P0[j]) {
	    NP++;
	    P += spower.power[j*spower.params->freqBins + i];
	    P2 += spower.power[j*spower.params->freqBins + i] * spower.power[j*spower.params->freqBins + i];
	  }
	}

	P /= (REAL8)NP;
	P2 /= (REAL8)NP;

	R = 2.0*P*P - P2;
	
	if(R>0.0) {
	  R = sqrt(R);
	} else {
	  R = 0.0;
	}

	rtp.Q[i] = R;
	rtp.P0[i] = P - R;
	
      }

      {
	UINT4 k;
	REAL4 *rho;
	
	rho = (REAL4 *)LALMalloc(spower.params->freqBins * sizeof(REAL4));
	if(!(rho)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}

	LALTFCRiceThreshold(status->statusPtr, rho, &rtp);
	CHECKSTATUSPTR (status);
      
	for(k=0;(int)k<(int)spower.params->freqBins;k++) {
	  dir.rho[k] = rho[k];
	}
      
	LALFree(rho);
      }

      LALFree(rtp.P0);
      LALFree(rtp.Q);
    }

    break;

  default:
    ABORT( status, STDBURSTSEARCHH_EUINP, STDBURSTSEARCHH_MSGEUINP);
  }

  /* run tfclusters */
  LALInitCList(status->statusPtr, &clist, &tspec); 
  CHECKSTATUSPTR (status);

  LALGetClusters(status->statusPtr, &clist, &spower, &dir);
  CHECKSTATUSPTR (status);

  LALFreeSpecgram(status->statusPtr, &spower);
  CHECKSTATUSPTR (status);

  LALInitCList(status->statusPtr, &list, &tspec); 
  CHECKSTATUSPTR (status);

  if(clist.nclusters > 0) { /* this should be fixed in lal */
    LALClustersPowerThreshold(status->statusPtr, &list, &clist, &dir);
    CHECKSTATUSPTR (status);
  }

  LALFreeCList(status->statusPtr, &clist);
  CHECKSTATUSPTR (status);

  /* generate output with start_time, duration, central_freq & bandwidth */
  bzero(output, sizeof(EventIDColumn));

  for(i=0; i<list.nclusters; i++) {

    SnglBurstTable *boutput;
    
    output->next = (EventIDColumn *)LALCalloc(1,sizeof(EventIDColumn));
    if(!(output->next)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}
    
    output = output->next;

    boutput = output->snglBurstTable = (SnglBurstTable *)LALCalloc(1,sizeof(SnglBurstTable));
    if(!(boutput)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}

    strncpy(boutput->search,search,LIGOMETA_SEARCH_MAX);
    strncpy(boutput->channel, channel, LIGOMETA_CHANNEL_MAX);
    strncpy(boutput->ifo, ifo, LIGOMETA_IFO_MAX);
    boutput->ifo[2] = 0;

    mint = maxt = list.t[i][0];
    minf = maxf = list.f[i][0];

    for(j=1; j<list.sizes[i]; j++) {
      if(list.t[i][j] < mint) mint = list.t[i][j];
      if(list.t[i][j] > maxt) maxt = list.t[i][j];
      if(list.f[i][j] < minf) minf = list.f[i][j];
      if(list.f[i][j] > maxf) maxf = list.f[i][j];
    }

    boutput->start_time.gpsSeconds = input->epoch.gpsSeconds + (INT4)floor((REAL4)mint * trez);
    boutput->start_time.gpsNanoSeconds = (INT4)floor(1E9*((REAL4)mint * trez - floor((REAL4)mint * trez))) + input->epoch.gpsNanoSeconds;
    if(boutput->start_time.gpsNanoSeconds >= 1000000000) {
      (boutput->start_time.gpsSeconds)++;
      (boutput->start_time.gpsNanoSeconds)-=1000000000;
    }

    boutput->duration = (REAL4)(maxt-mint+1) * trez;

    boutput->central_freq = list.params->flow + 0.5*(REAL4)(minf + maxf)/T;
    boutput->bandwidth = (REAL4)(maxf - minf + 1)/T;

    if(blob) {

      SnglTransdataTable *trans_data;

      trans_data = output->snglTransdataTable = (SnglTransdataTable *)LALCalloc(1,sizeof(SnglTransdataTable));
      if(!(output->snglTransdataTable)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}

      strncpy(trans_data->ifo, boutput->ifo, LIGOMETA_IFO_MAX);
      trans_data->ifo[2] = 0;
      strncpy(trans_data->name, boutput->channel, LIGOMETA_TRANSDATA_NAME_MAX);
      
      trans_data->dimensions = 2;

      trans_data->x_bins = list.sizes[i];
      trans_data->x_start = (REAL8)(boutput->start_time.gpsSeconds) + 1E-9*(REAL8)(boutput->start_time.gpsNanoSeconds);
      trans_data->x_end = trans_data->x_start + boutput->duration;      
      strncpy(trans_data->x_units, "GPS seconds", LIGOMETA_TRANSDATA_UNITS_MAX);

      trans_data->y_bins = list.sizes[i];
      trans_data->y_start = boutput->central_freq - 0.5*boutput->bandwidth;
      trans_data->y_end = boutput->central_freq + 0.5*boutput->bandwidth;
      strncpy(trans_data->y_units, "Hz", LIGOMETA_TRANSDATA_UNITS_MAX);

      strncpy(trans_data->data_type,"tf map", LIGOMETA_TRANSDATA_DATA_MAX);
      strncpy(trans_data->data_units,"arbitrary", LIGOMETA_TRANSDATA_DATA_MAX);
      
      trans_data->transdata_length = sizeof(UINT4) + list.sizes[i] * (2*sizeof(UINT4) + sizeof(REAL8));

      /* save cluster shape in trans_data */
      trans_data->trans_data = (UCHAR *)LALMalloc(sizeof(UINT4) + list.sizes[i] * (2*sizeof(UINT4) + sizeof(REAL8)));
      if(!(trans_data->trans_data)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);} 

      memcpy(trans_data->trans_data, &(list.sizes[i]), sizeof(UINT4));

#ifndef WORDS_BIGENDIAN
      endian_swap((char *)(trans_data->trans_data), sizeof(UINT4), 1);
#endif

      for(j=0; j<list.sizes[i]; j++) {
	memcpy(trans_data->trans_data + sizeof(UINT4)+(sizeof(REAL8)+2*sizeof(UINT4))*j, &(list.t[i][j]), sizeof(UINT4));
	memcpy(trans_data->trans_data + 2*sizeof(UINT4)+(sizeof(REAL8)+2*sizeof(UINT4))*j, &(list.f[i][j]), sizeof(UINT4));
	memcpy(trans_data->trans_data + 3*sizeof(UINT4)+(sizeof(REAL8)+2*sizeof(UINT4))*j, &(list.P[i][j]), sizeof(REAL8));

#ifndef WORDS_BIGENDIAN
	endian_swap((char *)(trans_data->trans_data + sizeof(UINT4)+(sizeof(REAL8)+2*sizeof(UINT4))*j), sizeof(UINT4), 2);
	endian_swap((char *)(trans_data->trans_data + 3*sizeof(UINT4)+(sizeof(REAL8)+2*sizeof(UINT4))*j), sizeof(REAL8), 1);
#endif

      }
      
    }

  }

  /* clean up */
  LALFreeCList(status->statusPtr, &list);
  CHECKSTATUSPTR (status);

  LALFreeCListDir(status->statusPtr, &dir);
  CHECKSTATUSPTR (status);

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


#ifndef WORDS_BIGENDIAN
static void endian_swap(char * pdata, int dsize, int nelements)

{

        int i,j,indx;
        char tempbyte;

        if (dsize <= 1) return;

        for (i=0; i<nelements; i++)
        {
                indx = dsize;
                for (j=0; j<dsize/2; j++)
                {
                        tempbyte = pdata[j];
                        indx = indx - 1;
                        pdata[j] = pdata[indx];
                        pdata[indx] = tempbyte;
                }

                pdata = pdata + dsize;
        }

        return;

}
#endif
