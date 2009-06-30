/*
*  Copyright (C) 2007 Jolien Creighton, Julien Sylvestre
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*-----------------------------------------------------------------------*
 *
 * File Name: StdBurstSearch.c
 *
 * Author: Julien Sylvestre
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------*/



/******** <lalVerbatim file="StdBurstSearchCV"> ********
Author: Sylvestre, J
$Id$
********* </lalVerbatim> ********/

#include <lal/StdBurstSearch.h>
#include <lal/FindRoot.h>
#include <lal/CoarseGrainFrequencySeries.h>
#include "lal/LALRCSID.h"
#include <lal/LPC.h>
#include <lal/AVFactories.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/Sort.h>
#include <math.h>
#include <strings.h>

NRCSID (STDBURSTSEARCHC, "$Id$");

#define HALF_LOG_2PI 0.918938533204673

static REAL8 LALi0(REAL8 x);
/*static REAL8 LALi1(REAL8 x); */
static REAL8 LALi0e(REAL8 x);
static REAL8 LALi1e(REAL8 x);

/******** <lalLaTeX file="StdBurstSearchC"> ********
\noindent
Compute standardized burst information.
\subsubsection*{Prototype}
********* </lalLaTeX> ********/
/* <lalVerbatim> */
void
LALBurstOutput(
	       LALStatus *status,
	       EventIDColumn *output,
	       EventIDColumn *Input,
	       BurstOutputParameters *params
	       ) {
/* </lalVerbatim> */
/******** <lalLaTeX file="StdBurstSearchC"> ********
\subsubsection*{Description}
Description of the function...

\subsubsection*{Uses}
\begin{verbatim}
...
\end{verbatim}
********* </lalLaTeX> ********/

  /* static variables for good speed; params == NULL resets them */
  static REAL4 *e1 = NULL,   /* lpc error, 1st half */
    *e2 = NULL;              /* lpc error, 2nd half */
  static UINT4 e1_end = 0;
  static UINT4 e2_start = 0; /* where 2nd half starts */
  static UINT4 e1_filter_delay = 0,
    e2_filter_delay = 0;
  static BurstOutputSpecStat *bsstat = NULL; /* linked list of spectrogram stats */
  static COMPLEX8Vector *Hvec = NULL; /* buffer for fft */
  static REAL4Vector *Dvec = NULL;    /* buffer for fft */
  static UINT4 Hlen = 0,
    Dlen = 0;

  BurstOutputDataSegment *data;

  INITSTATUS (status, "LALBurstOutput", STDBURSTSEARCHC);
  ATTATCHSTATUSPTR (status);

  if(!params) {
    /* reset */

    BurstOutputSpecStat *bptr;

    while(bsstat) {

      if(bsstat->P0) {
	LALFree(bsstat->P0);
      }

      if(bsstat->Q) {
	LALFree(bsstat->Q);
      }

      if(bsstat->wwin) {
	LALFree(bsstat->wwin);
      }

      if(bsstat->pfwd) {
	LALDestroyRealFFTPlan(status->statusPtr, &(bsstat->pfwd));
	CHECKSTATUSPTR (status);
      }

      if(bsstat->resp) {
	LALCDestroyVector(status->statusPtr, &(bsstat->resp->data));
	CHECKSTATUSPTR (status);
	LALFree(bsstat->resp);
      }

      bptr = bsstat;
      bsstat = bsstat->next;
      LALFree(bptr);

    }

    if(e1) {
      LALFree(e1);
    }

    if(e2) {
      LALFree(e2 + e2_start);
    }

    e1 = e2 = NULL;

    Dvec->length = Dlen;
    LALDestroyVector(status->statusPtr, &Dvec);
    CHECKSTATUSPTR (status);

    Hvec->length = Hlen;
    LALCDestroyVector( status->statusPtr, &Hvec);
    CHECKSTATUSPTR (status);

    DETATCHSTATUSPTR (status);
    RETURN (status);
  }

  /* trivial input check */
  ASSERT ( output, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  ASSERT ( Input, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  ASSERT ( params, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  ASSERT ( params->data, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  ASSERT ( params->data->data, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);
  ASSERT ( params->data->data->data, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);

  ASSERT( params->data->data->data->length == 1, status, STDBURSTSEARCHH_ENULLP, STDBURSTSEARCHH_MSGENULLP);


  data = params->data;



  if(!Dvec) {
    Dlen = data->data->data->vectorLength;
    LALCreateVector(status->statusPtr, &Dvec, data->data->data->vectorLength);
    CHECKSTATUSPTR (status);
  }

  if(!Hvec) {
    Hlen = data->data->data->vectorLength/2 + 1;
    LALCCreateVector( status->statusPtr, &Hvec, data->data->data->vectorLength/2 + 1);
    CHECKSTATUSPTR (status);
  }

  /* loop over input */
  Input = Input->next; /* first element in list is empty */

  while(Input) {

    LIGOTimeGPS   start_time;
    REAL4         duration = 0;
    REAL4         central_freq;
    REAL4         bandwidth;
    REAL4         amplitude;
    REAL4         snr;
    REAL4         confidence;
    REAL4 totalPower;

    SnglBurstTable *boutput, *input;

    input = Input->snglBurstTable;

    /* estimate start time, etc. */
    switch(params->method) {

    case 0:
      /* plain copy from ETG */
      start_time = input->start_time;
      duration = input->duration;
      central_freq = input->central_freq;
      bandwidth = input->bandwidth;
      amplitude = 0.0;
      snr = 0.0;
      confidence = 0.0;

      break;

    case 1:
      /* the real thing */
/******** <lalLaTeX file="StdBurstSearchC"> ********
\subsection*{Algorithm}
\begin{itemize}
********* </lalLaTeX> ********/

      {
	BurstOutputSpecStat *bptr, *blg;
	UINT4 start_index;
	UINT4 mfi;
	REAL8 mlik;
	REAL8 llikf, llikt;
	INT4 llo, lhi;
	REAL4 *r4ptr;
	UINT4 i, k, l, nTime;
	INT4 didel;
	UINT4 filter_delay;
	REAL8Vector *FPower;

	/*****************************************************************/
	/**                      default parameters                     **/
	/*****************************************************************/

/******** <lalLaTeX file="StdBurstSearchC"> ********
\item Hardcoded parameters:
********* </lalLaTeX> ********/
/* <lalVerbatim> */
	UINT4 wforder = 16; /* order of FIR whitening filter */
	UINT4 eolap = (UINT4)(0.5/data->data->deltaT); /* e1 goes from 0 to N/2; e2 goes from N/2-eolap to N-1 */
	REAL8 duration_increase = 2.0; /* factor increasing duration of a burst in psd calculations, to account for windowing */

	REAL4 bandSafe = 1000.0 * data->data->deltaT; /* extra data for bandpass, in seconds */
	REAL4 bandDf = 6e-3*0.5/data->data->deltaT;  /* width of transition band at 500 Hz; scaled linearly in frequency */
	REAL4 bandAttenuation = 0.01; /* attenuation outside of band */

	REAL4 bandwidthPF = 0.5; /* fraction of burst power in bandwidth */

	REAL4 durationPF = 0.5; /* fraction of burst power in duration */
/* </lalVerbatim> */
	/*****************************************************************/
	/**                          frequency                          **/
	/*****************************************************************/

/******** <lalLaTeX file="StdBurstSearchC"> ********
\item The code keeps a bank of statistics for the power in frequency bands, with one spectrogram for every burst duration (within an accuracy of 10 $\mu$s). The spectrograms are constructed starting at the beginning of the analysis segment, and include the ETG trigger.
********* </lalLaTeX> ********/

	/* scan to see if we have the stats for that duration */
	blg = NULL;
	bptr = bsstat;
	while(bptr) {
	  if(fabs(bptr->duration - input->duration) < 1E-5) {
	    break;
	  }
	  blg = bptr;
	  bptr = bptr->next;
	}

	if(!bptr) { 	  /* we don't have the stats; create them */

	  RealFFTPlan *pfwd = NULL;
	  REAL4 *wwin;
	  REAL4 nn2;
	  FrequencySamplingParams fsp;

	  bptr = (BurstOutputSpecStat *)LALCalloc(1, sizeof(BurstOutputSpecStat));
	  if(!(bptr)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}

	  if(blg) {
	    blg->next = bptr;
	  } else {
	    bsstat = bptr;
	  }

/******** <lalLaTeX file="StdBurstSearchC"> ********
To create a new entry for a given burst duration, a segment size of (burst duration * duration\_increase) is defined.
********* </lalLaTeX> ********/

	  bptr->duration = input->duration;
	  bptr->nTime = (UINT4)floor(input->duration * duration_increase / data->data->deltaT);
	  if(bptr->nTime < 2) {
	    bptr->nTime = 2;
	  }
	  didel = (INT4)floor((REAL4)bptr->nTime - input->duration / data->data->deltaT);
	  if(didel < 0) {
	    didel = 0;
	  }

	  bptr->nFreq = bptr->nTime / 2;

	  bptr->P0 = (REAL8 *)LALCalloc(bptr->nFreq, sizeof(REAL8));
	  if(!(bptr->P0)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}

	  bptr->Q = (REAL8 *)LALCalloc(bptr->nFreq, sizeof(REAL8));
	  if(!(bptr->Q)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}

	  LALCreateForwardRealFFTPlan(status->statusPtr, &pfwd, bptr->nTime, 0 );
	  CHECKSTATUSPTR (status);
	  bptr->pfwd = pfwd;

/******** <lalLaTeX file="StdBurstSearchC"> ********
A parabolic window function is used.
********* </lalLaTeX> ********/

	  bptr->wwin = wwin = (REAL4 *)LALMalloc(bptr->nTime * sizeof(REAL4));
	  if(!(wwin)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}

	  nn2 = 0.5*(REAL4)bptr->nTime;
	  bptr->norm = 0.0;
	  for(k=0; k<bptr->nTime; k++) {
	    wwin[k] = 1.0 - pow(((REAL4)k-nn2)/nn2,2.0);
	    bptr->norm += wwin[k] * wwin[k];
	  }

	  /* for psd normalization */
	  bptr->norm = pow(bptr->norm, 2.0);


/******** <lalLaTeX file="StdBurstSearchC"> ********
The response function is resampled to the frequency resolution defined by the segment duration.
********* </lalLaTeX> ********/

	  bptr->resp = (COMPLEX8FrequencySeries *)LALCalloc(1, sizeof(COMPLEX8FrequencySeries));
	  if(!(bptr->resp)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}

	  LALCCreateVector(status->statusPtr, &(bptr->resp->data), bptr->nFreq);
	  CHECKSTATUSPTR (status);
	  fsp.length = bptr->nFreq;
	  fsp.f0 = 0.0;
	  fsp.deltaF = 0.5 / (data->data->deltaT * (REAL8)bptr->nFreq);
	  if(data->resp) {
	    LALCCoarseGrainFrequencySeries(status->statusPtr, bptr->resp, data->resp, &fsp);
	    CHECKSTATUSPTR (status);
	  }

	  Dvec->length = bptr->nTime;
	  Hvec->length = bptr->nTime / 2 + 1;

/******** <lalLaTeX file="StdBurstSearchC"> ********
The whole analysis segment is divided in 50\%-overlapping subsegments which are windowed, FFTed, squared and normalized, multiplied by the response function.
********* </lalLaTeX> ********/

	  for(k=0;k<2*data->data->data->vectorLength / bptr->nTime - 1; k++) {

	    memcpy(Dvec->data, data->data->data->data + k*bptr->nTime/2, Dvec->length * sizeof(REAL4));

	    /* apply window */
	    for(l=0;l<bptr->nTime;l++) {
	      Dvec->data[l] *= wwin[l];
	    }

	    LALForwardRealFFT(status->statusPtr, Hvec, Dvec, pfwd);
	    CHECKSTATUSPTR (status);

	    for(l=0;l<bptr->nFreq;l++) {
	      REAL8 P = (Hvec->data[l].re*Hvec->data[l].re + Hvec->data[l].im*Hvec->data[l].im) / bptr->norm;

	      if(data->resp) {
		P *= pow(bptr->resp->data->data[l].re,2.0) + pow(bptr->resp->data->data[l].im,2.0);
	      }

	      bptr->P0[l] += P;
	      bptr->Q[l] += P*P;
	    }

	  }

/******** <lalLaTeX file="StdBurstSearchC"> ********
The sum and sum-squared of the power in each frequency band are saved and used to compute the P0 and Q parameters of the Rice distribution.
********* </lalLaTeX> ********/

	  for(l=0;l<bptr->nFreq;l++) {

	    REAL8 Pm, P2, R;

	    Pm = bptr->P0[l] / (REAL8)k;
	    P2 = bptr->Q[l] / (REAL8)k;

	    R = 2.0*Pm*Pm - P2;

	    if(R>0.0) {
	      R = sqrt(R);
	    } else {
	      R = 0.0;
	    }

	    bptr->Q[l] = R;
	    bptr->P0[l] = Pm - R;

	  }

	} else {
	  didel = (INT4)floor((REAL4)bptr->nTime - input->duration / data->data->deltaT);
	  if(didel < 0) {
	    didel = 0;
	  }
	}


/******** <lalLaTeX file="StdBurstSearchC"> ********
\item The frequency representation of the segment containing the burst is calculated exactly as above.
********* </lalLaTeX> ********/
	/* now compute fft of burst */
	start_index = (UINT4)floor((REAL8)(input->start_time.gpsSeconds - data->data->epoch.gpsSeconds)/data->data->deltaT + 1E-9*((REAL8)input->start_time.gpsNanoSeconds - (REAL8)data->data->epoch.gpsNanoSeconds)/data->data->deltaT);

	if((int)start_index >= (int)didel) {
	  start_index -= didel; /* correct for duration_increase */
	} else {
	  start_index = 0;
	}

	memcpy(Dvec->data, data->data->data->data + start_index, bptr->nTime * sizeof(REAL4));
	Dvec->length = bptr->nTime;

	Hvec->length = bptr->nTime / 2 + 1;

	for(l=0;l<bptr->nTime;l++) {
	  Dvec->data[l] *= bptr->wwin[l];
	}

	LALForwardRealFFT(status->statusPtr, Hvec, Dvec, bptr->pfwd);
	CHECKSTATUSPTR (status);


/******** <lalLaTeX file="StdBurstSearchC"> ********
\item For all frequencies in the band identified by the ETG, a maximum likelihood estimation (using the Rice fit) of the signal power is produced.
********* </lalLaTeX> ********/
	/* get max likelihood */
	mfi = 0;
	mlik = -1.0;
	llikf = 0.0;
	snr = 0.0;
	totalPower = 0.0;

	llo = (INT4)floor((input->central_freq - 0.5 * input->bandwidth) * data->data->deltaT * (REAL8)bptr->nTime);
	if(llo < 0) {
	  llo = 0;
	}

	lhi = (INT4)ceil((input->central_freq + 0.5 * input->bandwidth) * data->data->deltaT * (REAL8)bptr->nTime);
	if((int)lhi > (int)bptr->nFreq) {
	  lhi = bptr->nFreq;
	}

	if(lhi < llo) {
	  lhi = llo;
	}

	FPower = NULL;
	if(lhi>llo) {
	  LALDCreateVector(status->statusPtr, &(FPower),lhi-llo);
	  CHECKSTATUSPTR (status);
	}

	for(l=llo;(int)l<(int)lhi;l++) {
	  RiceLikelihoodParams rp;
	  REAL8 Pmax, q0;
	  DFindRootIn dfri;

	  rp.P = (Hvec->data[l].re*Hvec->data[l].re + Hvec->data[l].im*Hvec->data[l].im) / bptr->norm;

	  if(data->resp) {
	    rp.P *= pow(bptr->resp->data->data[l].re,2.0) + pow(bptr->resp->data->data[l].im,2.0);
	  }

	  rp.P0 = bptr->P0[l];
	  rp.Q = bptr->Q[l];

	  LALRiceLikelihood(status->statusPtr, &q0, 0.0, &rp);
	  if(status->statusPtr->statusCode) {
	    char msg[1024];
	    float bw = 1.0 / (data->data->deltaT * (REAL4)bptr->nTime);
	    float fr = (REAL4)l / (data->data->deltaT * (REAL4)bptr->nTime);
	    sprintf(msg,"Error in LALRiceLikelihood at frequency %g Hz, BW = %g Hz; %g %g %g %g %g %g",fr,bw,bptr->norm,Hvec->data[l].re,Hvec->data[l].im,rp.P0,rp.Q,rp.P);
	    ABORT(status, 111, msg);
	  }
	  /*
	  CHECKSTATUSPTR (status);
	  */

	  if(q0 > 0.0) {
	    dfri.xacc = 1e-4 * rp.P;
	    dfri.xmin = 0.0;
	    dfri.xmax = rp.P;
	    dfri.function = LALRiceLikelihood;

	    LALDBisectionFindRoot(status->statusPtr, &Pmax, &dfri, &rp);
	    CHECKSTATUSPTR (status);
	  } else {
	    Pmax = 0.0;
	  }

	  if(Pmax > mlik) {
	    mlik = Pmax;
	    mfi = l;
	  }

	  /*
	  fprintf(stderr,"%g\t%g\t%g\t%g\t%g\t%g\n",(REAL4)l / (data->data->deltaT * (REAL4)bptr->nTime),q0,rp.P0,rp.Q,rp.P,Pmax);
	  */

	  snr += Pmax / (rp.P0 + rp.Q);
	  totalPower += Pmax;

	  {
	    REAL8 ar = 2.0*sqrt(rp.P*rp.Q)/rp.P0;

	    if(ar > 30.0) {
	      llikf += -log(rp.P0) - (rp.P + rp.Q)/rp.P0 + 27.3847;
	    } else {
	      llikf += -log(rp.P0) - (rp.P + rp.Q)/rp.P0 + log(LALi0(ar));
	    }

	  }

	  FPower->data[l-llo] = Pmax;

	}


/******** <lalLaTeX file="StdBurstSearchC"> ********
\item The snr is the sqrt of the sum (over in-band frequencies) of (estimated signal power / background noise (i.e., P0) + line power (i.e. Q)).
********* </lalLaTeX> ********/
	/* get linear snr */
	snr = sqrt(snr);

/******** <lalLaTeX file="StdBurstSearchC"> ********
\item central\_freq is set to the in-band frequency with the largest signal power.
********* </lalLaTeX> ********/
	/* set "central_freq" to maximum of power */
	if(mfi>0) {
	  central_freq = (REAL4)mfi / (data->data->deltaT * (REAL4)bptr->nTime);
	} else {
	  central_freq = input->central_freq;
	}

/******** <lalLaTeX file="StdBurstSearchC"> ********
\item amplitude is set to sqrt(total in-band power).
********* </lalLaTeX> *******/

	amplitude = sqrt(totalPower * (data->data->deltaT * bptr->nTime));

/*set amplitude to sqrt(maximum of power); get strain per rtHz
  amplitude = sqrt(mlik * 2.0 * data->data->deltaT); */

	/*****************************************************************/
	/**                           bandwidth                         **/
	/*****************************************************************/
/******** <lalLaTeX file="StdBurstSearchC"> ********
\item The bandwidth is the smallest band containing the peak frequency such that 50% of the estimated purst power is in that band.
********* </lalLaTeX> ********/

	if(FPower) {
	  REAL8 tPower = 0.0;
	  REAL8 thr;
	  INT4 pind, sdur, tdel;
	  INT4 msdur, mpeak = mfi - llo;

	  msdur = FPower->length - mpeak;
	  if(mpeak < msdur) {
	    msdur = mpeak;
	  }

	  for(sdur=0;(int)sdur<(int)FPower->length;sdur++) {
	    tPower += FPower->data[sdur];
	  }
	  thr = tPower * (1.0 - bandwidthPF);

	  for(sdur = 1; sdur < msdur; sdur++) {
	    tPower = 0.0;
	    for(pind = mpeak - sdur + 1; pind < mpeak; pind++) {
	      tPower += FPower->data[pind];
	    }

	    if(tPower >= thr) {
	      break;
	    }

	    for(tdel = 1; tdel < sdur; tdel++) {
	      tPower -= FPower->data[mpeak - sdur + tdel];
	      tPower += FPower->data[mpeak + tdel];

	      if(tPower >= thr) {
		break;
	      }
	    }

	  }

	  bandwidth = (REAL4)sdur / (data->data->deltaT * (REAL4)bptr->nTime);
	} else {

	  bandwidth = 0.0;

	}

	if(FPower) {
	  LALDDestroyVector(status->statusPtr, &FPower);
	  CHECKSTATUSPTR (status);
	}

	/*****************************************************************/
	/**                             time                            **/
	/*****************************************************************/
/******** <lalLaTeX file="StdBurstSearchC"> ********
\item In the time domain, the full segment is splitted in to overlapping subsegments, and a whitening filter is trained on every subsegment. The impulse response of each filter is measured.
********* </lalLaTeX> ********/

	start_index = (UINT4)floor((REAL8)(input->start_time.gpsSeconds - data->data->epoch.gpsSeconds)/data->data->deltaT + 1E-9*((REAL8)input->start_time.gpsNanoSeconds - (REAL8)data->data->epoch.gpsNanoSeconds)/data->data->deltaT);
	nTime = (UINT4)floor(input->duration / data->data->deltaT);
	if(nTime == 0) {
	  nTime = 1;
	}

	if(!e1) {
	  /* need to train the whitening filter */
	  REAL4Vector a1, a2; /* coefficients of whitening filter for 1st and 2nd halves */
	  REAL4Vector e1v, e2v; /*vectors of errors */
	  REAL4Vector zerov;
	  REAL4 *e1t, *e2t;
	  REAL4 r40;
	  REAL4IIRFilter f; /* whitening filter */
	  CHAR dummyName[2] = {1,0};
	  REAL4Vector *delta = NULL;

	  /* compute errors */
	  f.name = (CHAR *)dummyName;
	  f.deltaT = -1.0; /* take from data */
	  f.recursCoef = &zerov;
	  zerov.length = 0;
	  zerov.data = &r40;
	  r40 = 0.0;

	  /* first half */
	  e1_end = e1v.length = data->data->data->vectorLength / 2;
	  e1 = e1v.data = (REAL4 *)LALCalloc(e1v.length, sizeof(REAL4));
	  if(!e1) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}
	  memcpy(e1, data->data->data->data, e1v.length * sizeof(REAL4));

	  e1t = (REAL4 *)LALCalloc(e1v.length, sizeof(REAL4));
	  if(!e1t) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}
	  memcpy(e1t, data->data->data->data, e1v.length * sizeof(REAL4));

	  a1.data = NULL;
	  a1.length = wforder;

	  LALLPC(status->statusPtr, &a1, &e1v, wforder);
	  CHECKSTATUSPTR (status);

	  f.directCoef = &a1;

	  f.history = NULL;
	  LALCreateVector(status->statusPtr, &(f.history), 1+a1.length);
	  CHECKSTATUSPTR (status);
	  bzero(f.history->data, f.history->length * sizeof(REAL4));

	  a1.data[0] = 0.0;
	  for(i=1;i<a1.length;i++) {
	    a1.data[i] = -a1.data[i];
	  }

	  LALIIRFilterREAL4Vector(status->statusPtr, &e1v, &f);
	  CHECKSTATUSPTR (status);


	  /* estimate filter delay */
	  LALCreateVector(status->statusPtr, &delta, wforder);
	  CHECKSTATUSPTR (status);

	  bzero(delta->data+1, (delta->length-1)*sizeof(REAL4));
	  delta->data[0] = 1.0;

	  bzero(f.history->data, f.history->length * sizeof(REAL4));

	  LALIIRFilterREAL4Vector(status->statusPtr, delta, &f);
	  CHECKSTATUSPTR (status);

	  {
	    REAL4 tr4;
	    UINT4 im = 0;
	    tr4 = fabs(delta->data[0]);
	    for(i=1;i<delta->length;i++) {
	      if(fabs(delta->data[i]) > tr4) {
		tr4 = fabs(delta->data[i]);
		im = i;
	      }
	    }
	    e1_filter_delay = im;
	  }


	  for(i=0; i<e1v.length; i++) {
	    e1[i] = e1t[i] - e1[i];
	  }

	  /* second half */
	  e2_start = e1v.length - eolap;

	  e2v.length = data->data->data->vectorLength - e1v.length + eolap;
	  e2 = e2v.data = (REAL4 *)LALCalloc(e2v.length, sizeof(REAL4));
	  if(!e2) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}
	  memcpy(e2, data->data->data->data + e2_start, e2v.length * sizeof(REAL4));

	  e2t = (REAL4 *)LALCalloc(e2v.length, sizeof(REAL4));
	  if(!e2t) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}
	  memcpy(e2t, data->data->data->data + e2_start, e2v.length * sizeof(REAL4));

	  a2.data = NULL;
	  a2.length = wforder;

	  LALLPC(status->statusPtr, &a2, &e2v, wforder);
	  CHECKSTATUSPTR (status);

	  f.directCoef = &a2;

	  bzero(f.history->data, f.history->length * sizeof(REAL4));

	  a2.data[0] = 0.0;
	  for(i=1;i<a2.length;i++) {
	    a2.data[i] = -a2.data[i];
	  }

	  LALIIRFilterREAL4Vector(status->statusPtr, &e2v, &f);
	  CHECKSTATUSPTR (status);

	  /* estimate filter delay */
	  bzero(delta->data+1, (delta->length-1)*sizeof(REAL4));
	  delta->data[0] = 1.0;

	  bzero(f.history->data, f.history->length * sizeof(REAL4));

	  LALIIRFilterREAL4Vector(status->statusPtr, delta, &f);
	  CHECKSTATUSPTR (status);

	  {
	    REAL4 tr4;
	    UINT4 im = 0;
	    tr4 = fabs(delta->data[0]);
	    for(i=1;i<delta->length;i++) {
	      if(fabs(delta->data[i]) > tr4) {
		tr4 = fabs(delta->data[i]);
		im = i;
	      }
	    }
	    e2_filter_delay = im;
	  }


	  LALDestroyVector(status->statusPtr, &delta);
	  CHECKSTATUSPTR (status);

	  for(i=0; i<e2v.length; i++) {
	    e2[i] = e2t[i] - e2[i];
	  }

	  e2 = e2 - e2_start;

	  if(f.history) {
	    LALDestroyVector(status->statusPtr, &(f.history));
	    CHECKSTATUSPTR (status);
	    f.history = NULL;
	  }

	  /* clean up */
	  LALFree(a1.data);
	  LALFree(a2.data);

	  LALFree(e1t);
	  LALFree(e2t);

	}

/******** <lalLaTeX file="StdBurstSearchC"> ********
\item For each trigger, the data are then whitened and bandpassed to the band identified by the ETG.
********* </lalLaTeX> ********/

	/* look at whitened data near the burst */
	llo = start_index;
	lhi = start_index + nTime;

	if(start_index >= e2_start) {
	  filter_delay = e2_filter_delay;
	  r4ptr = e2;

	  if((int)lhi > (int)data->data->data->vectorLength) {
	    lhi = data->data->data->vectorLength;
	  }
	} else {
	  r4ptr = e1;
	  filter_delay = e1_filter_delay;

	  if((int)lhi > (int)e1_end) {
	    lhi = e1_end;
	  }
	}

	{
	  REAL4TimeSeries banddata;
	  REAL4Vector *bdat = NULL;
	  INT4 fi,li;
	  PassBandParamStruc hipas, lopas;
	  REAL4 norm = input->bandwidth * 2.0 * data->data->deltaT;

	  banddata.deltaT = data->data->deltaT;
	  fi = llo - (INT4)ceil(bandSafe / banddata.deltaT);
	  if((int)llo >= (int)e2_start && (int)fi < (int)e2_start) {
	    fi = e2_start;
	  }
	  if(fi < 0) {
	    fi = 0;
	  }

	  li = lhi + (INT4)ceil(bandSafe / banddata.deltaT);
	  if((int)llo < (int)e2_start && (int)li > (int)e1_end) {
	    li = e1_end;
	  }
	  if((int)li > (int)data->data->data->vectorLength) {
	    li = data->data->data->vectorLength;
	  }

	  LALCreateVector(status->statusPtr,&bdat,li-fi+1);
	  CHECKSTATUSPTR (status);

	  memcpy(bdat->data, r4ptr + fi, (li-fi+1)*sizeof(REAL4));
	  banddata.data = bdat;

	  /* bandpass filter */
	  hipas.nMax = 0;
	  hipas.f2 = input->central_freq - 0.5 * input->bandwidth;
	  hipas.a2 = 0.95;
	  hipas.f1 = hipas.f2 - bandDf * input->central_freq / 500.0;
	  hipas.a1 = bandAttenuation;

	  if(hipas.f1 > 0.0) {
	    LALDButterworthREAL4TimeSeries(status->statusPtr,&banddata,&hipas);
	    CHECKSTATUSPTR (status);
	  }

	  lopas.nMax = 0;
	  lopas.f1 = input->central_freq + 0.5 * input->bandwidth;
	  lopas.a1 = 0.95;
	  lopas.f2 = lopas.f1 + bandDf * input->central_freq / 500.0;
	  lopas.a2 = bandAttenuation;
	  if(lopas.f2 < 0.5/banddata.deltaT) {
	    LALDButterworthREAL4TimeSeries(status->statusPtr,&banddata,&lopas);
	    CHECKSTATUSPTR (status);
	  }

	  r4ptr = bdat->data - fi;

	  for(l=llo+1; (int)l<(int)lhi; l++) {
	    r4ptr[l] /= norm;
	  }

	  mfi = llo;
	  mlik = r4ptr[llo];

	  if(fabs(r4ptr[llo]) < 1E15) {
	    llikt = -HALF_LOG_2PI -0.5*pow(r4ptr[llo],2.0);
	  } else {
	    llikt = -HALF_LOG_2PI -0.5e30;
	  }

	  for(l=llo+1; (int)l<(int)lhi; l++) {

	    if(fabs(r4ptr[l]) > mlik) {
	      mlik = fabs(r4ptr[l]);
	      mfi = l;
	    }

	    if(fabs(r4ptr[l]) < 1E15) {
	      llikt += -HALF_LOG_2PI -0.5*pow(r4ptr[l],2.0);
	    } else {
	      llikt += -HALF_LOG_2PI -0.5e30;
	    }

	  }

	  /* correct estimated time for filter delay */
	  mfi -= filter_delay;

/******** <lalLaTeX file="StdBurstSearchC"> ********
\item The time of maximum excursion in the resulting time series is corrected for the filter delay, and is save in start\_time.
********* </lalLaTeX> ********/

	  /* put max time into "start_time" */
	  start_time.gpsSeconds = data->data->epoch.gpsSeconds + (INT4)floor((REAL8)mfi * data->data->deltaT);
	  start_time.gpsNanoSeconds = (INT4)floor(1E9*((REAL8)mfi * data->data->deltaT - floor((REAL8)mfi * data->data->deltaT))) + data->data->epoch.gpsNanoSeconds;
	  if(start_time.gpsNanoSeconds >= 1000000000) {
	    (start_time.gpsSeconds)++;
	    (start_time.gpsNanoSeconds)-=1000000000;
	  }


	/*****************************************************************/
	/**                          duration                           **/
	/*****************************************************************/
/******** <lalLaTeX file="StdBurstSearchC"> ********
\item The duration is computed exactly like the bandwidth, except with the whitened and bandpassed data.
********* </lalLaTeX> ********/

	  FPower = NULL;
	  if(lhi > llo+1) {
	    LALDCreateVector(status->statusPtr, &(FPower),lhi-llo-1);
	    CHECKSTATUSPTR (status);
	  }

	  for(l=llo+1; (int)l<(int)lhi; l++) {
	    FPower->data[l-llo-1] = r4ptr[l] * r4ptr[l];
	  }

	  if(FPower) {
	    REAL8 tPower = 0.0;
	    REAL8 thr;
	    INT4 pind, sdur, tdel;
	    INT4 msdur, mpeak = mfi - llo - 1;

	    msdur = FPower->length - mpeak;
	    if(mpeak < msdur) {
	      msdur = mpeak;
	    }

	    for(sdur=0;(int)sdur<(int)FPower->length;sdur++) {
	      tPower += FPower->data[sdur];
	    }
	    thr = tPower * (1.0 - durationPF);

	    for(sdur = 1; sdur < msdur; sdur++) {
	      tPower = 0.0;
	      for(pind = mpeak - sdur + 1; pind < mpeak; pind++) {
		tPower += FPower->data[pind];
	      }

	      if(tPower >= thr) {
		break;
	      }

	    for(tdel = 1; tdel < sdur; tdel++) {
	      tPower -= FPower->data[mpeak - sdur + tdel];
	      tPower += FPower->data[mpeak + tdel];

	      if(tPower >= thr) {
		break;
	      }
	    }

	  }

	  duration = (REAL4)sdur * data->data->deltaT;

	  }

	  if(FPower) {
	    LALDDestroyVector(status->statusPtr, &FPower);
	    CHECKSTATUSPTR (status);
	  }

	  LALDestroyVector(status->statusPtr,&bdat);
	  CHECKSTATUSPTR (status);
	}


	/*****************************************************************/
	/**                         confidence                          **/
	/*****************************************************************/
/******** <lalLaTeX file="StdBurstSearchC"> ********
\item The confidence is the log10 of the confidence calculated in the frequency domain.
********* </lalLaTeX> ********/
/*
	if(llikf < llikt) {
	  confidence = llikf;
	} else {
	  confidence = llikt;
	}
*/
	confidence = llikf;

	confidence /= log(10.0);

      }

      break;

    default:

      ABORT ( status, STDBURSTSEARCHH_EIN, STDBURSTSEARCHH_MSGEIN);
    }

    /* allocate output for event*/
    output->next = (EventIDColumn *)LALCalloc(1,sizeof(EventIDColumn));
    if(!(output->next)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}
    output = output->next;

    boutput = output->snglBurstTable = (SnglBurstTable *)LALCalloc(1,sizeof(SnglBurstTable));
    if(!(boutput)) {ABORT(status, STDBURSTSEARCHH_EMEM, STDBURSTSEARCHH_MSGEMEM);}

    /* copy results in output event */
    if(!(params->skip & STDBURSTSEARCHSKIP_STARTTIME)) {
      boutput->start_time = start_time;
    }
    if(!(params->skip & STDBURSTSEARCHSKIP_DURATION)) {
      boutput->duration = duration;
    }
    if(!(params->skip & STDBURSTSEARCHSKIP_CENTRALFREQ)) {
      boutput->central_freq = central_freq;
    }
    if(!(params->skip & STDBURSTSEARCHSKIP_BANDWIDTH)) {
      boutput->bandwidth = bandwidth;
    }
    if(!(params->skip & STDBURSTSEARCHSKIP_AMPLITUDE)) {
      boutput->amplitude = amplitude;
    }
    if(!(params->skip & STDBURSTSEARCHSKIP_SNR)) {
      boutput->snr = snr;
    }
    if(!(params->skip & STDBURSTSEARCHSKIP_CONFIDENCE)) {
      boutput->confidence = confidence;
    }

    /*
    printf("%u:%u\t%g\t%g\t%g\t%g\n",start_time.gpsSeconds,start_time.gpsNanoSeconds, central_freq, confidence, snr, amplitude);
    */

    strncpy(boutput->ifo,input->ifo,LIGOMETA_IFO_MAX);
    strncpy(boutput->search,input->search,LIGOMETA_SEARCH_MAX);
    strncpy(boutput->channel,input->channel,LIGOMETA_CHANNEL_MAX);

    if(Input->snglTransdataTable) {
      output->snglTransdataTable = Input->snglTransdataTable;
    }

    /* advance to next event */
    Input = Input->next;

  }

/******** <lalLaTeX file="StdBurstSearchC"> ********
\end{itemize}
********* </lalLaTeX> ********/

  /* Normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}

void
LALRiceLikelihood(
		  LALStatus *status,
		  REAL8 *llik,
		  REAL8 Qs,
		  void *Params
		  ) {
  REAL8 x, Q;

  RiceLikelihoodParams *params = (RiceLikelihoodParams *)Params;

  INITSTATUS (status, "LALRiceLikelihood", STDBURSTSEARCHC);
  ATTATCHSTATUSPTR (status);

  Q = params->Q + Qs;

  x = 2.0*sqrt(params->P*Q)/params->P0;

  if(Q>0.0) {
    /*
    *llik = exp(-(params->P+Q)/params->P0) * (-LALi0(x) + params->P*LALi1(x)/sqrt(params->P*Q)) / (params->P0 * params->P0);
    */
    *llik = exp(-pow(sqrt(params->P)-sqrt(Q),2.0)/params->P0)*(-LALi0e(x) + sqrt(params->P/Q)*LALi1e(x)) / (params->P0 * params->P0);
  } else {
    *llik = exp(-params->P/params->P0) * (params->P-params->P0) / pow(params->P0,3.0);
  }

  /* JC: isnan is not allowed.
  if(isnan(*llik)) {
    CHAR ebuf[2048];
    sprintf(ebuf,"Q=%g P=%g P0=%g x=%g *llik=%g io=%g i1=%g\n",Q,params->P,params->P0,x,*llik,LALi0e(x),LALi1e(x));
    ABORT(status,100,ebuf);
  }
  */

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


/* Chebyshev coefficients for exp(-x) LALI0(x)
 * in the interval [0,8].
 *
 * lim(x->0){ exp(-x) LALI0(x) } = 1.
 */

static REAL8 A[] =
{
-4.41534164647933937950E-18,
 3.33079451882223809783E-17,
-2.43127984654795469359E-16,
 1.71539128555513303061E-15,
-1.16853328779934516808E-14,
 7.67618549860493561688E-14,
-4.85644678311192946090E-13,
 2.95505266312963983461E-12,
-1.72682629144155570723E-11,
 9.67580903537323691224E-11,
-5.18979560163526290666E-10,
 2.65982372468238665035E-9,
-1.30002500998624804212E-8,
 6.04699502254191894932E-8,
-2.67079385394061173391E-7,
 1.11738753912010371815E-6,
-4.41673835845875056359E-6,
 1.64484480707288970893E-5,
-5.75419501008210370398E-5,
 1.88502885095841655729E-4,
-5.76375574538582365885E-4,
 1.63947561694133579842E-3,
-4.32430999505057594430E-3,
 1.05464603945949983183E-2,
-2.37374148058994688156E-2,
 4.93052842396707084878E-2,
-9.49010970480476444210E-2,
 1.71620901522208775349E-1,
-3.04682672343198398683E-1,
 6.76795274409476084995E-1
};


/* Chebyshev coefficients for exp(-x) sqrt(x) LALI0(x)
 * in the inverted interval [8,infinity].
 *
 * lim(x->inf){ exp(-x) sqrt(x) LALI0(x) } = 1/sqrt(2pi).
 */


static REAL8 B[] =
{
-7.23318048787475395456E-18,
-4.83050448594418207126E-18,
 4.46562142029675999901E-17,
 3.46122286769746109310E-17,
-2.82762398051658348494E-16,
-3.42548561967721913462E-16,
 1.77256013305652638360E-15,
 3.81168066935262242075E-15,
-9.55484669882830764870E-15,
-4.15056934728722208663E-14,
 1.54008621752140982691E-14,
 3.85277838274214270114E-13,
 7.18012445138366623367E-13,
-1.79417853150680611778E-12,
-1.32158118404477131188E-11,
-3.14991652796324136454E-11,
 1.18891471078464383424E-11,
 4.94060238822496958910E-10,
 3.39623202570838634515E-9,
 2.26666899049817806459E-8,
 2.04891858946906374183E-7,
 2.89137052083475648297E-6,
 6.88975834691682398426E-5,
 3.36911647825569408990E-3,
 8.04490411014108831608E-1
};

static REAL8 chbevl(REAL8 x, REAL8 array[], INT4 n);

static REAL8 LALi0(REAL8 x)
{
REAL8 y;

if( x < 0 )
        x = -x;
if( x <= 8.0 )
        {
        y = (x/2.0) - 2.0;
        return( exp(x) * chbevl( y, A, 30 ) );
        }

return(  exp(x) * chbevl( 32.0/x - 2.0, B, 25 ) / sqrt(x) );

}




static REAL8 LALi0e(REAL8 x)
{
REAL8 y;

if( x < 0 )
        x = -x;
if( x <= 8.0 )
        {
        y = (x/2.0) - 2.0;
        return( chbevl( y, A, 30 ) );
        }

return(  chbevl( 32.0/x - 2.0, B, 25 ) / sqrt(x) );

}


static REAL8 chbevl(REAL8 x, REAL8 array[], INT4 n)
{
REAL8 b0, b1, b2, *p;
INT4 i;

p = array;
b0 = *p++;
b1 = 0.0;
i = n - 1;

do
        {
        b2 = b1;
        b1 = b0;
        b0 = x * b1  -  b2  + *p++;
        }
while( --i );

return( 0.5*(b0-b2) );
}


static REAL8 AA[] =
{
 2.77791411276104639959E-18,
-2.11142121435816608115E-17,
 1.55363195773620046921E-16,
-1.10559694773538630805E-15,
 7.60068429473540693410E-15,
-5.04218550472791168711E-14,
 3.22379336594557470981E-13,
-1.98397439776494371520E-12,
 1.17361862988909016308E-11,
-6.66348972350202774223E-11,
 3.62559028155211703701E-10,
-1.88724975172282928790E-9,
 9.38153738649577178388E-9,
-4.44505912879632808065E-8,
 2.00329475355213526229E-7,
-8.56872026469545474066E-7,
 3.47025130813767847674E-6,
-1.32731636560394358279E-5,
 4.78156510755005422638E-5,
-1.61760815825896745588E-4,
 5.12285956168575772895E-4,
-1.51357245063125314899E-3,
 4.15642294431288815669E-3,
-1.05640848946261981558E-2,
 2.47264490306265168283E-2,
-5.29459812080949914269E-2,
 1.02643658689847095384E-1,
-1.76416518357834055153E-1,
 2.52587186443633654823E-1
};

static REAL8 BB[] =
{
 7.51729631084210481353E-18,
 4.41434832307170791151E-18,
-4.65030536848935832153E-17,
-3.20952592199342395980E-17,
 2.96262899764595013876E-16,
 3.30820231092092828324E-16,
-1.88035477551078244854E-15,
-3.81440307243700780478E-15,
 1.04202769841288027642E-14,
 4.27244001671195135429E-14,
-2.10154184277266431302E-14,
-4.08355111109219731823E-13,
-7.19855177624590851209E-13,
 2.03562854414708950722E-12,
 1.41258074366137813316E-11,
 3.25260358301548823856E-11,
-1.89749581235054123450E-11,
-5.58974346219658380687E-10,
-3.83538038596423702205E-9,
-2.63146884688951950684E-8,
-2.51223623787020892529E-7,
-3.88256480887769039346E-6,
-1.10588938762623716291E-4,
-9.76109749136146840777E-3,
 7.78576235018280120474E-1
};

/*
static REAL8 LALi1(REAL8 x)
{
REAL8 y, z;

z = fabs(x);
if( z <= 8.0 )
	{
	y = (z/2.0) - 2.0;
	z = chbevl( y, AA, 29 ) * z * exp(z);
	}
else
	{
	z = exp(z) * chbevl( 32.0/z - 2.0, BB, 25 ) / sqrt(z);
	}
if( x < 0.0 )
	z = -z;
return( z );
}
*/

static REAL8 LALi1e( REAL8 x )
{
REAL8 y, z;

z = fabs(x);
if( z <= 8.0 )
	{
	y = (z/2.0) - 2.0;
	z = chbevl( y, AA, 29 ) * z;
	}
else
	{
	z = chbevl( 32.0/z - 2.0, BB, 25 ) / sqrt(z);
	}
if( x < 0.0 )
	z = -z;
return( z );
}
