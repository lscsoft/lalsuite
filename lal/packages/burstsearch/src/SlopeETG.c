/*-----------------------------------------------------------------------*
 *
 * File Name: SlopeETG.c
 *
 * Author: Julien Sylvestre
 *
 * Revision: $Id$ 
 *
 *-----------------------------------------------------------------------*/



/******** <lalVerbatim file="SlopeETGCV"> ********
Author: Sylvestre, J
$Id$
********* </lalVerbatim> ********/

#include <lal/StdBurstSearch.h>
#include <lal/Sort.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/TFClusters.h>
#include <lal/AVFactories.h>

#include "lal/LALRCSID.h"

#include <strings.h>

NRCSID (SLOPEETGC, "$Id$");

/******** <lalLaTeX file="SlopeETGC"> ********
\noindent
Implement the SLOPE event trigger generator.
\subsubsection*{Prototype}
********* </lalLaTeX> ********/
/* <lalVerbatim> */
void
LALSlopeETG(
	    LALStatus *status, 
	    EventIDColumn *output, 
	    REAL4TimeVectorSeries *input, 
	    BurstParameter *params
	    ) {
/* </lalVerbatim> */
/******** <lalLaTeX file="SlopeETGC"> ********
\subsubsection*{Description}
Description of the parameters: 
\begin{center}
\begin{tabular}{l|l|l}
parameter index & type & description \\ \hline 
1 & string & channel name \\
2 & REAL4 & threshold (if negative, normalized by -rms) \\
3 & REAL4Vector & coefficients \\
4 & INT4 & points in cluster \\
5 & REAL4 & min frequency to report (Hz) \\
6 & REAL4 & max frequency to report (Hz) 
\end{tabular}
\end{center}

\subsubsection*{Uses}
\begin{verbatim}
...a bunch of stuff.
\end{verbatim}
********* </lalLaTeX> ********/
  REAL4 thr, minf, maxf;
  REAL4Vector *coef;
  INT4 nclus;
  UINT4 i;
  REAL8 sum2;

  REAL4Vector buffer;
  REAL4Vector data1, *data = &data1;

  CHAR ifo[LIGOMETA_IFO_MAX];
  CHAR channel[LIGOMETA_CHANNEL_MAX];

  SnglBurstTable *boutput;

  INITSTATUS (status, "LALSlopeETG", SLOPEETGC);
  ATTATCHSTATUSPTR (status);

  /* trivial input check */
  ASSERT ( output, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  ASSERT ( input, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  ASSERT ( params, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);

  ASSERT ( input->data->length == 1, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);

  data->length = input->data->vectorLength;
  data->data = input->data->data;

  /* parse parameters */
#define SetParameter(par,type) if(!params) {ABORT ( status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI);} \
  if(!(params->type)) {ABORT ( status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI); } \
  par = *(params->type); \
  params = params->next
  
#define SetStringParameter(par) if(!params) {ABORT ( status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI);} \
  if(!(params->char_)) {ABORT(status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI);} \
  ASSERT ( par, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP); \
  strcpy(par,params->char_); \
  params = params->next

#define SetVectorParameter(par) if(!params) {ABORT ( status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI);} \
  if(!(params->real4vector_)) {ABORT(status, STDBURSTSEARCHH_ENULLPI, STDBURSTSEARCHH_MSGENULLPI);} \
  par = params->real4vector_; \
  params = params->next
  
  params = params->next;
  
  SetStringParameter(channel);
  strncpy(ifo,channel,2);

  SetParameter(thr,real4_);
  SetVectorParameter(coef);
  SetParameter(nclus,int4_);
  SetParameter(minf,real4_);
  SetParameter(maxf,real4_);

#undef SetParameter
#undef SetStringParameter
#undef SetVectorParameter


  /***********SLOPE code************/
  /* save a copy of the input */
  buffer.length = data->length;
  buffer.data = (REAL4 *)LALMalloc(buffer.length * sizeof(REAL4));
  if(!(buffer.data)) {
    ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM); 
  }
  memcpy(buffer.data, data->data, buffer.length * sizeof(REAL4));
  
  /* filter */
  {
    REAL4IIRFilter f;
    REAL4Vector zerov;
    REAL4 r40;

    f.deltaT = -1.0; /* take from data */
    f.recursCoef = &zerov;
    zerov.length = 0;
    zerov.data = &r40;
    r40 = 0.0;
    
    f.directCoef = coef;
    
    f.history = NULL;
    LALCreateVector(status->statusPtr, &(f.history), 1+coef->length);
    CHECKSTATUSPTR (status);
    bzero(f.history->data, f.history->length * sizeof(REAL4));

    LALIIRFilterREAL4Vector(status->statusPtr, &buffer, &f);
    CHECKSTATUSPTR (status);

    LALDestroyVector(status->statusPtr, &(f.history));
    CHECKSTATUSPTR (status);
  }

  /* adaptive threshold */
  if(thr<0.0) {
    for(i=0, sum2=0.0; i<buffer.length; i++) {
      sum2 += pow(buffer.data[i],2.0);
    }
    sum2 /= (REAL8)buffer.length; 
    sum2 = sqrt(sum2); /* RMS */
    thr *= -sum2;
  }

  /* search filtered data and generate output with start_time, duration, central_freq & bandwidth */
  bzero(output, sizeof(SnglBurstTable));

  for(i=0; i<buffer.length; i++) {

    if(buffer.data[i] >= thr) {

      UINT4 i0 = i;
      UINT4 i1 = 0;

      i++;

      while(i<buffer.length) {
	if(buffer.data[i] < thr) {
	  i1 = i-1;
	  break;
	}
	if(i >= i0 + nclus) {
	  i1 = i;
	  i--; /* to avoid missing one sample */
	  break;
	}
	
	i++;
      }

      /* generate output */
      output->next = (EventIDColumn *)LALCalloc(1,sizeof(EventIDColumn));
      if(!(output->next)) {
	ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM); 
      }

      output = output->next;

      boutput = output->snglBurstTable = (SnglBurstTable *)LALCalloc(1,sizeof(SnglBurstTable));
      if(!(boutput)) {
	ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM); 
      }

      strncpy(boutput->channel, channel, LIGOMETA_CHANNEL_MAX);
      strncpy(boutput->ifo, ifo, LIGOMETA_IFO_MAX);
      boutput->ifo[2] = 0;

      boutput->start_time.gpsSeconds = input->epoch.gpsSeconds + (INT4)floor((REAL4)i0 * input->deltaT);
      boutput->start_time.gpsNanoSeconds = (INT4)floor(1E9*((REAL4)i0 * input->deltaT - floor((REAL4)i0 * input->deltaT))) + input->epoch.gpsNanoSeconds;
      if(boutput->start_time.gpsNanoSeconds >= 1000000000) {
	(boutput->start_time.gpsSeconds)++;
	(boutput->start_time.gpsNanoSeconds)-=1000000000;
      }

      boutput->duration = (REAL4)(i1-i0+1) * input->deltaT;

      boutput->central_freq = 0.5*(minf + maxf);
      boutput->bandwidth = maxf - minf;

    } /* end if above threshold */
  } /* end loop over data */

  /* clean up */
  LALFree(buffer.data);

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);

}
