/*  <lalVerbatim file="GEO600InspiralSearchFuncsCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{GEO600InspiralSearchFuncs.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{ReadInspiralFindEventsInCP}
  \idx{ReadInspiralFindEventsIn}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GEO600InspiralSearchFuncsCV}}

</lalLaTeX>  */

#include <config.h>
#if defined HAVE_GSL_GSL_FFT_REAL_H && defined HAVE_MYSQL_H

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/AVFactories.h>
#include <time.h>
#include <mysql.h>
#include <sys/types.h>
#include <unistd.h>
#include <gsl/gsl_fft_real.h>
 
#include "GEO600InspiralSearch.h"
#include "GEO600FrameChannels.h"
 
NRCSID (INSPIRALSEARCHC, "$Id$");

/*  <lalVerbatim file="AppendLALProcessTableCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{AppendLALProcessTable.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{AppendLALProcessTableCP}
\idx{AppendLALProcessTable()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{AppenLALProcessTableCV}}

</lalLaTeX>  */

 
/*  <lalVerbatim file="AppendLALProcessTableCP"> */
void
AppendLALProcessTable (LALDBConnect conn) 
{  /*  </lalVerbatim>  */


	MYSQL mysql;
	char string[128];
	time_t tt;
	char *tm;
	
	time(&tt);
	tm = ctime(&tt);

	sprintf(string, "update process set end_time = '%s' where process_id = %d", 
		         tm, conn.insert_id);

        mysql_init(&mysql);
        if(!(mysql_real_connect(&mysql,conn.host,conn.user,conn.password,
             conn.database,0,NULL,0)))
             fprintf(stderr, "error connecting, code: %s\n",mysql_error(&mysql));

	if (mysql_query(&mysql,string))
            fprintf(stderr, "error writing to database, code: %s\n",mysql_error(&mysql));


	mysql_close(&mysql);
	return;

}

/*  <lalVerbatim file="ReadInspiralFindEventsInCP"> */
void 
ReadInspiralFindEventsIn
   (
   LALStatus            *status, 
   InspiralFindEventsIn *eventsin
   )
{  /*  </lalVerbatim>  */

   INITSTATUS (status, "ReadInspiralFindEventsIn", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);
   ASSERT (eventsin,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

   /* 
    * eventsin->Threshold 
    * VERY, VERY, IMPORTANT: Don't choose the threshold too small
    * as otherwise we will end up with too many events producing
    * memory allocation failure; ideally should be >= 5
    */
   eventsin->Threshold = 10.0;
   eventsin->ClusterThreshold = 5.0;
   eventsin->displayPSD = 0;
   eventsin->displayData = 0;
   eventsin->displayTemplates = 0;
   eventsin->displayCorrelation = 0;
   eventsin->displayCorrelationStats = 0;

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/*  <lalVerbatim file="ReadInspiralTemplateBankParamsCP"> */
void 
ReadInspiralTemplateBankParams
   (
   LALStatus            *status, 
   InspiralCoarseBankIn *coarseIn
   )
{  /*  </lalVerbatim>  */

   INITSTATUS (status, "ReadInspiralTemplateBankParams", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);
   ASSERT (coarseIn,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   coarseIn->mMin =   4.0;
   coarseIn->mMax =   6.0;
   coarseIn->MMax =  coarseIn->mMax*2.;
   coarseIn->massRange =  MinComponentMassMaxTotalMass;
   coarseIn->massRange =  MinMaxComponentMass;
   coarseIn->mmCoarse = 0.97;
   coarseIn->mmFine = 0.99;
   coarseIn->fLower = 55.;
   coarseIn->fUpper = 1700.;
   coarseIn->iflso = 0;
   coarseIn->tSampling = 16384.;
   coarseIn->order = twoPN;
   coarseIn->approximant = TaylorT3;
   coarseIn->space = Tau0Tau3;
   coarseIn->etamin = coarseIn->mMin * ( coarseIn->MMax - coarseIn->mMin) /
      pow(coarseIn->MMax,2.);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/*  <lalVerbatim file="ReadInspiralSearchParamsCP"> */
void 
ReadInspiralSearchParams
   (
   LALStatus            *status, 
   InspiralSearchParams *searchParams
   )
{  /*  </lalVerbatim>  */

   float samplingRate;
   INT4 isuccess;

   INITSTATUS (status, "ReadInspiralSearchParams", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);
   ASSERT(searchParams,status,LALNOISEMODELSH_ENULL,LALNOISEMODELSH_MSGENULL);

		   
   /*
   searchParams->channelName = "SEI_NBC_SEIS-X";
   searchParams->channelName = "LSC_MID_EP-P";
   searchParams->findEventsType = LoudestEvent;
   searchParams->findEventsType = AllTemplates;
   */
   searchParams->findEventsType = AllEvents;
   searchParams->startGPSTime = 715434000;
   searchParams->endGPSTime =   715434200;
   searchParams->channelName = "G1:DER_H_HP-EP";
   searchParams->currentGPSTime = searchParams->startGPSTime;

   searchParams->numDataSetsForPSD = 8;         /* should be >=1 */
   searchParams->paddingFactor = 4;             /* should be >=4 */

   searchParams->fLo =  55.;                    /* lower frequency cutoff (in Hz) for bandpass */
   searchParams->fHi = 1700.;                   /* upper frequency cutoff (in Hz) for bandpass */

   /*
    * The parameters below are used only by Monte-Carlo codes
   */
   searchParams->NumTrials = 10;
   searchParams->PSDType = syntheticPSD;
   searchParams->PSDType = realPSD;

   /*
    * Set the base directory where data is situated 
    */
   if (searchParams->PSDType == realPSD)
   {
	   isuccess = SetBaseDirectory("/data/frames/");
	   if (isuccess) 
	   {
		   fprintf(stderr, "Cannot open base directory to read data\n");
		   exit(0);
	   }

	   /*
	    * Get the sample rate for the specified Channel 
	    */
	   isuccess=GetSampleRate(searchParams->currentGPSTime, searchParams->channelName, &samplingRate);
	   searchParams->samplingRate = (UINT4) samplingRate;
   
	   if (isuccess) 
	   {
		   fprintf(stderr, "Cannot find channel %s in frame file\n", searchParams->channelName);
		   exit(0);
	   }
   }

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/*  <lalVerbatim file="SetUpProcessTableInfoCP"> */
void
SetUpProcessTableInfo (LALProcessTable *procTable) 
{  /*  </lalVerbatim>  */


	time_t tt;
	char *tm;
	
	time(&tt);
	tm = ctime(&tt);

	procTable->creator_db = 1;
	procTable->program = "GEO600InspiralSearch";
	procTable->version = "0.1";
	procTable->cvs_repository = "lal"; 
	procTable->cvs_entry_time = 1;
	procTable->comment = "Playground 2-lockPSD";
	procTable->is_online = 0;
	procTable->node = "localhost";
	procTable->username = "sathya";
	procTable->unix_procid = getpid(); 
	procTable->start_time = tm;
	procTable->end_time = "null";
	procTable->jobid = 10; 
	procTable->node = "Cardiff";
	procTable->process_id = "null";
	procTable->param_set = 0;
	procTable->ifos = "G1";

	return;
       
	
}

/*  <lalVerbatim file="ReadRandomInspiralSignalParamsCP"> */
void 
ReadRandomInspiralSignalParams
   (
   LALStatus              *status, 
   RandomInspiralSignalIn *randIn,
   InspiralCoarseBankIn   *coarseIn
   )
{  /*  </lalVerbatim>  */

   InspiralBankParams bankPars;

   REAL8 mass1Sq, mass2Sq, spin1Frac, spin2Frac, spin1Theta, spin1Phi, spin2Theta, spin2Phi;

   INITSTATUS (status, "ReadRandomInspiralSignalParams", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);

   ASSERT (coarseIn,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (randIn,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

   LALInspiralSetSearchLimits(status->statusPtr, &bankPars, *coarseIn);
   randIn->type = 0;
   randIn->SignalAmp = 10.;
   randIn->NoiseAmp = 1.;
   randIn->useed = 274829;
   randIn->t0Min = bankPars.x0Min;
   randIn->t0Max = bankPars.x0Max;
   randIn->tnMin = bankPars.x1Min;
   randIn->tnMax = bankPars.x1Max;
   randIn->param.startTime = 0.;
   randIn->param.startPhase = 0.;
   randIn->param.nStartPad = 0;
   randIn->param.nEndPad = 0;
   randIn->param.signalAmplitude = 1.;
   randIn->param.ieta = 1;
   randIn->param.Theta = 0.;
   randIn->param.OmegaS = 0.;
   randIn->mMin = coarseIn->mMin;
   randIn->mMax = coarseIn->mMax;
   randIn->MMax = coarseIn->MMax;
   randIn->param.fLower = coarseIn->fLower;
   randIn->param.fCutoff = coarseIn->fUpper;
   randIn->param.tSampling = coarseIn->tSampling;
   randIn->param.order = coarseIn->order;
   randIn->param.approximant = EOB;
   randIn->param.massChoice = m1Andm2;
   randIn->param.massChoice = t03;

   randIn->param.sourceTheta = LAL_PI/6.L;
   randIn->param.sourcePhi = LAL_PI/6.L;
   randIn->param.distance = 1.e8 * LAL_PC_SI/LAL_C_SI;
   spin1Frac = 0.9;
   spin2Frac = 0.9;

   mass1Sq = pow(randIn->param.mass1*LAL_MTSUN_SI,2.L);
   mass2Sq = pow(randIn->param.mass2*LAL_MTSUN_SI,2.L);
   randIn->param.orbitTheta0 = 0.L;
   randIn->param.orbitPhi0 = 0.L;

   spin1Theta = LAL_PI_2/3.L;
   spin1Phi = -LAL_TWOPI/4.L; 
   spin2Theta = LAL_PI_2/5.L; 
   spin2Phi = -LAL_TWOPI/6.L;

   randIn->param.spin1[0] =  mass1Sq * spin1Frac * sin(spin1Theta) * cos(spin1Phi);
   randIn->param.spin1[1] =  mass1Sq * spin1Frac * sin(spin1Theta) * sin(spin1Phi);
   randIn->param.spin1[2] =  mass1Sq * spin1Frac * cos(spin1Theta);
   randIn->param.spin2[0] =  mass2Sq * spin2Frac * sin(spin2Theta) * cos(spin2Phi);
   randIn->param.spin2[1] =  mass2Sq * spin2Frac * sin(spin2Theta) * sin(spin2Phi);
   randIn->param.spin2[2] =  mass2Sq * spin2Frac * cos(spin2Theta);
   if (randIn->param.massChoice != m1Andm2) 
   {
      if (coarseIn->space==Tau0Tau2) 
         randIn->param.massChoice = t02;
      else
         randIn->param.massChoice = t03;
   }
   randIn->etaMin = randIn->mMin*(randIn->MMax - randIn->mMin)/pow(randIn->MMax,2.);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/*  <lalVerbatim file="RunningAverageDoublePSDCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{RunningAverageDoublePSD.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{RunningAverageDoublePSDCP}
\idx{RunningAverageDoublePSD()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{RunningAverageDoublePSDCV}}

</lalLaTeX>  */

/*  <lalVerbatim file="RunningAverageDoublePSDCP"> */
void 
RunningAverageDoublePSD
   (
   LALStatus           *status, 
   REAL8Vector         *psd, 
   INT4                *ret, 
   RunningAveragePSDIn *psdIn
   )
{  /*  </lalVerbatim>  */

   int currentTime, nSecs, dataSetLength;
   char* chan;
   REAL8Vector *dataset=NULL, *buffer1=NULL;
   REAL4Vector *buffer2=NULL;
   REAL4Vector *buffer3=NULL;
   int nSetsby2, n, j, i, nby2;
   REAL8 dt, x, nSets, df, norm;

   INITSTATUS (status, "RunningAveragePSD", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);

   ASSERT(psd,status,LALNOISEMODELSH_ENULL,LALNOISEMODELSH_MSGENULL);
   ASSERT(psdIn,status,LALNOISEMODELSH_ENULL,LALNOISEMODELSH_MSGENULL);

   dt = 1./psdIn->searchParams.samplingRate;
   n = dataSetLength = psdIn->searchParams.dataSetLength;
   nSecs = psdIn->searchParams.dataSetDuration;
   chan = psdIn->searchParams.channelName;
   nSetsby2 = psdIn->searchParams.numDataSetsForPSD/2;
   nby2 = dataSetLength/2;

   ASSERT(psd->length == (UINT4) (nby2+1), status, GEO600INSPIRALSEARCHH_ECHOICE, GEO600INSPIRALSEARCHH_MSGECHOICE);
   ASSERT(psdIn->norm != 0, status, GEO600INSPIRALSEARCHH_ECHOICE, GEO600INSPIRALSEARCHH_MSGECHOICE);

   LALDCreateVector(status->statusPtr, &dataset, dataSetLength);
   LALDCreateVector(status->statusPtr, &buffer1, dataSetLength);
   LALCreateVector(status->statusPtr, &buffer2, dataSetLength);
   LALCreateVector(status->statusPtr, &buffer3, dataSetLength);

   nSets = 2*nSetsby2 + 1;
   df = 1./(double)nSecs;

   /*
   fprintf(stderr, "________________________________\n");
   fprintf(stderr, "DataSet CurrentTime for PSD=%d\n", psdIn->searchParams.currentGPSTime);
   fprintf(stderr, "________________________________\n");
   */
   /* Overlapping won't work properly if nSecs is less than 2 because if nSecs<2 we will
    * have to deal with fractions of seconds; so stop code if nSecs < 2.
    */
   if (nSecs < 2) exit(1);

   for (j=-nSetsby2; j<=nSetsby2; j++)
   {

      if (j != 0)
      {

	      int k;
	      k = j % nSetsby2;
	      if (j<0) currentTime = psdIn->searchParams.currentGPSTime - (int)((1. - (float)k/2.)*nSecs);
	      else currentTime = psdIn->searchParams.currentGPSTime + (int)((1. + (float)k/2.)*nSecs);
	      fprintf(stderr, "%d %d %d %e\n", j, k, currentTime, psdIn->norm);
	      *ret=GetChannelDoubleData(currentTime, nSecs, chan, buffer1->data, dataSetLength);
	      if (*ret) 
	      {
		      goto endAnalysis;
	      }
	      for (i=0; i<dataSetLength; i++) 
	      {
		      buffer1->data[i] *= psdIn->window->data[i];
	      }
      
	      /* double precision FFT using gsl */ 
	      gsl_fft_real_radix2_transform(buffer1->data, 1, buffer1->length);
	      /*
	      for (i=0; i<buffer1->length; i++) printf("%e\n", buffer1->data[i]);
	      for (i=0; i<buffer1->length; i++) buffer2->data[i] = (float) (buffer1->data[i]);
	      LALREAL4VectorFFT(status->statusPtr, buffer3, buffer2, psdIn->revp);
	      printf("&\n");
	      for (i=0; i<buffer3->length; i++) printf("%e\n", buffer3->data[i]/dataSetLength);
	      printf("&\n");
	      exit(0);
	      */
   
   /* Now need to convert the double precision result of this
    * fft to single precision for the rest of the code to 
    * work with */

   for (i=0; i<(INT4) buffer1->length; ++i)
   {
	   dataset->data[i] = (REAL4)buffer1->data[i];
   }

   /*
	      LALREAL4VectorFFT(status->statusPtr, dataset, buffer1, psdIn->fwdp);
  */
      
	      if (j==-nSetsby2)
	      { 
		      for (i=1; i<nby2; i++) 
		      {
			      x = dataset->data[i]*dataset->data[i] 
				      + dataset->data[n-i]*dataset->data[n-i];
			      psd->data[i] = x;
		      }
   
		      psd->data[0] = dataset->data[0] * dataset->data[0];
		      psd->data[nby2] = dataset->data[nby2] * dataset->data[nby2];
      
	      } else {
         
		      for (i=1; i<nby2; i++) 
		      {
			      x = dataset->data[i]*dataset->data[i] 
				      + dataset->data[n-i]*dataset->data[n-i];
			      psd->data[i] += x;
		      }
		      psd->data[0] += dataset->data[0] * dataset->data[0];
		      psd->data[nby2] += dataset->data[nby2] * dataset->data[nby2];
	      }
      }
   }
		      
   /* Calibration is not reliable at frequencies below and above 
    * fmin, fmax, respectively; therefore set the PSD to 'zero' at
    * those frequencies so that correlations, etc., will not have 
    * contribution from those regions. (The correlation programme
    * LALInspiralWaveCorrelate, ignores those bins of PSD that are zero.)
    */
   norm = 2.*dt/( (double) (psdIn->norm*(nSets-1)) );

   for (i=0; i<= nby2; i++) 
   {
	   REAL8 f;

	   f = i*df;

	   if (f < psdIn->searchParams.fLo || f > psdIn->searchParams.fHi) 
	   {
		   psd->data[i] = 0.;
	   }
	   else
	   {
		   psd->data[i] *= norm;
	   }
	   /* Checking the support of the response in the time domain */
	   /*
	   INT4 k;
	   if (i==0) k = 0; else k = n-i;
	   if (psd->data[i]>0) {
		   buffer2->data[i] = (float) (cos(2.*LAL_PI*f*nSecs/2.)/sqrt(psd->data[i]));
		   buffer2->data[k] = (float) (sin(2.*LAL_PI*f*nSecs/2.)/sqrt(psd->data[i]));
	   }
	   else 
	   {
		   buffer2->data[i] = buffer2->data[k] = 0.;
	   }
	   */
   }
	      
   /*
   for (i=1; i<psd->length; i++) printf("%e %e\n", df*i, psd->data[i]);
   */
   /*
   LALREAL4VectorFFT(status->statusPtr, buffer3, buffer2, psdIn->revp);
   for (i=0; i<n; i++) printf("%e %e\n", i+1., buffer3->data[i]);
   */

   endAnalysis:

   LALDDestroyVector(status->statusPtr, &dataset);
   LALDDestroyVector(status->statusPtr, &buffer1);
   LALDestroyVector(status->statusPtr, &buffer2);
   LALDestroyVector(status->statusPtr, &buffer3);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/*  <lalVerbatim file="SetUpDBConnectionInfoCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{SetUpDBConnectionInfo.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SetUpDBConnectionInfoCP}
\idx{SetUpDBConnectionInfo()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{SetUpDBConnectionInfoCV}}

</lalLaTeX>  */

/*  <lalVerbatim file="SetUpDBConnectionInfo"> */
void
SetUpDBConnectionInfo (LALDBConnect *conn) 
{  /*  </lalVerbatim>  */


	conn->host = "localhost";
	conn->database = "inspirals";
	conn->user = "sathya";
	conn->password = "sathya"; 
       
	return;
	
}


/*  <lalVerbatim file="WriteInspiralEventsCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{WriteInspiralEvents.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{WriteInspiralEventsCP}
\idx{WriteInspiralEvents()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{WriteInspiralEventsCV}}

</lalLaTeX>  */

/*  <lalVerbatim file="WriteInspiralEventsCP"> */
void
WriteInspiralEvents
   (
   LALStatus          *status, 
   InspiralEventsList *eventList,
   INT4               nEvents,
   FILE               *evtfile
   )
{  /*  </lalVerbatim>  */

   INT4 i;

   INITSTATUS (status, "WriteInspiralEvents", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);
   ASSERT (evtfile,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

   for (i=0; i<nEvents; i++)
   {
      fprintf(evtfile, "%10d %10d %10d %10d %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10d\n", 
	eventList[i].endTime, 
	eventList[i].endTimeNS, 
	eventList[i].impulseTime, 
	eventList[i].impulseTimeNS, 
	eventList[i].amplitude, 
	eventList[i].effDistance,
        eventList[i].phase,   
        eventList[i].param.mass1,   
        eventList[i].param.mass2,   
        eventList[i].param.chirpMass,
        eventList[i].param.eta,   
	eventList[i].param.t0, 
	eventList[i].param.t2, 
	eventList[i].param.t3,
        eventList[i].param.t4,   
        eventList[i].param.t5,   
        eventList[i].param.tC,   
	eventList[i].snr, 
	eventList[i].chisq,
	eventList[i].sigmasq, 
        eventList[i].chisqDOF); 
   }
   fflush(evtfile);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
/*  <lalVerbatim file="WriteInspiralEventsDBCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{WriteInspiralEventsDB.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{WriteInspiralEventsDBCP}
\idx{WriteInspiralEventsDB()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{WriteInspiralEventsDBCV}}

</lalLaTeX>  */

/*  <lalVerbatim file="WriteInspiralEventsDBCP"> */
void
WriteInspiralEventsDB
   (
   LALStatus          *status, 
   InspiralEventsList *eventList,
   INT4               nEvents,
   LALDBConnect       conn
   )
{  /*  </lalVerbatim>  */

   INT4 i;
   MYSQL mysql;
   char string[512];

   INITSTATUS (status, "WriteInspiralEvents", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);


   mysql_init(&mysql);
   if(!(mysql_real_connect(&mysql,conn.host,conn.user,conn.password,
        conn.database,0,NULL,0)))
        fprintf(stderr, "error connecting, code: %s\n",mysql_error(&mysql));

   for (i=0; i<nEvents; i++)
   {
   sprintf(string, "insert into snglInspiral values (%d,%d,'%s','%s','%s',%d,%d,%d,%d,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%d,%e,'%s')",
   1,
   conn.insert_id,
   "filter_id?",
   "GEO",
   "TemplatedSearch",
   eventList[i].endTime,
   eventList[i].endTimeNS,
   eventList[i].impulseTime,
   eventList[i].impulseTimeNS,
   eventList[i].amplitude,   
   eventList[i].effDistance,   
   eventList[i].phase,   
   eventList[i].param.mass1,   
   eventList[i].param.mass2,   
   eventList[i].param.chirpMass,
   eventList[i].param.eta,   
   eventList[i].param.t0,   
   eventList[i].param.t2,   
   eventList[i].param.t3,   
   eventList[i].param.t4,   
   eventList[i].param.t5,   
   eventList[i].param.tC,   
   eventList[i].snr,   
   eventList[i].chisq,   
   eventList[i].chisqDOF,   
   eventList[i].sigmasq,   
   "null");
 
   if (mysql_query(&mysql,string))
       fprintf(stderr, "error writing to database, code: %s\n",mysql_error(&mysql));

   }

   mysql_close(&mysql);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

/*  <lalVerbatim file="WriteInspiralFindEventsInCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{WriteInspiralFindEventsIn.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{WriteInspiralFindEventsInCP}
\idx{WriteInspiralFindEventsIn()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{WriteInspiralFindEventsInCV}}

</lalLaTeX>  */

/*  <lalVerbatim file="WriteInspiralFindEventsInCP"> */
void 
WriteInspiralFindEventsIn
   (
   LALStatus            *status, 
   FILE                 *logfile,
   InspiralFindEventsIn *eventsin
   )
{  /*  </lalVerbatim>  */

   INITSTATUS (status, "WriteInspiralFindEventsIn", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);
   ASSERT (logfile,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

   fprintf(logfile, "____________________________________________________\n\n"); 
   fprintf(logfile, "        Threshold = %e\n", eventsin->Threshold);
   fprintf(logfile, "____________________________________________________\n\n"); 
   fflush(logfile);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}


/*  <lalVerbatim file="WriteInspiralSearchParamsCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{WriteInspiralSearchParams.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{WriteInspiralSearchParamsCP}
\idx{WriteInspiralSearchParams()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{WriteInspiralSearchParamsCV}}

</lalLaTeX>  */

/*  <lalVerbatim file="WriteInspiralSearchParamsCP"> */
void 
WriteInspiralSearchParams
   (
   LALStatus            *status, 
   FILE                 *logfile,
   InspiralSearchParams *searchParams
   )
{  /*  </lalVerbatim>  */

   INITSTATUS (status, "WriteInspiralSearchParams", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);
   ASSERT(logfile, status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

   fprintf(logfile, "____________________________________________________\n\n"); 
   fprintf (logfile,"               InspiralSearchParams                 \n");
   fprintf(logfile, "____________________________________________________\n\n"); 

   /*
   fprintf (logfile, "         psdModel = %d\n", searchParams->psdModel);
   fprintf (logfile, "    simulatedData = %d\n", searchParams->simulatedData);
   */
   if (searchParams->channelName!=NULL)
   {
   fprintf (logfile, "      channelName = %s\n", searchParams->channelName);
   }
   fprintf (logfile, "     startGPSTime = %d\n", searchParams->startGPSTime);
   fprintf (logfile, "       endGPSTime = %d\n", searchParams->endGPSTime);
   fprintf (logfile, "   currentGPSTime = %d\n", searchParams->currentGPSTime);
   fprintf (logfile, "numDataSetsForPSD = %d\n", searchParams->numDataSetsForPSD);
   fprintf (logfile, "    paddingFactor = %d\n", searchParams->paddingFactor);
   fprintf (logfile, "              fLo = %f\n", searchParams->fLo);
   fprintf (logfile, "              fHi = %f\n", searchParams->fHi);
   fprintf (logfile, "        NumTrials = %d\n", searchParams->NumTrials);
   fprintf (logfile, "          PSDType = %d\n", searchParams->PSDType);

   fflush(logfile);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
/*  <lalVerbatim file="WriteInspiralTemplateBankParamsCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{WriteInspiralTemplateBankParams.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{WriteInspiralTemplateBankParamsCP}
\idx{WriteInspiralTemplateBankParams()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{WriteInspiralTemplateBankParamsCV}}

</lalLaTeX>  */

/*  <lalVerbatim file="WriteInspiralTemplateBankParamsCP"> */
void 
WriteInspiralTemplateBankParams
   (
   LALStatus            *status, 
   FILE                 *logfile,
   InspiralCoarseBankIn *coarseIn
   )
{  /*  </lalVerbatim>  */


   INITSTATUS (status, "WriteInspiralTemplateBankParams", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);
   ASSERT (logfile,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

   /*
   char psdchoice[64];
   sprintf(psdchoice, "No known PSD");

   if (coarseIn->NoisePsd == LALGEOPsd)
   {
	   sprintf(psdchoice, "LALGEOPsd");
   }
   if (coarseIn->NoisePsd == LALLIGOIPsd)
   {
	   sprintf(psdchoice, "LALLIGOIPsd");
   }
   if (coarseIn->NoisePsd == LALTAMAPsd)
   {
	   sprintf(psdchoice, "LALTAMAPsd");
   }
   if (coarseIn->NoisePsd == LALVIRGOPsd)
   {
	   sprintf(psdchoice, "LALVIRGOPsd");
   }
   */

   fprintf(logfile, "____________________________________________________\n\n"); 
   fprintf(logfile, "               InspiralTemplateBankParams           \n"); 
   fprintf(logfile, "____________________________________________________\n\n"); 

   /*
   fprintf (logfile, "NoisePsdModel = %s\n", psdchoice);
   */
   fprintf (logfile, "         mMin = %e\n", coarseIn->mMin);
   fprintf (logfile, "         mMax = %e\n", coarseIn->mMax);
   fprintf (logfile, "         MMax = %e\n", coarseIn->MMax);
   fprintf (logfile, "     mmCoarse = %e\n", coarseIn->mmCoarse);
   fprintf (logfile, "       mmFine = %e\n", coarseIn->mmFine);
   fprintf (logfile, "       fLower = %e\n", coarseIn->fLower);
   fprintf (logfile, "       fUpper = %e\n", coarseIn->fUpper);
   fprintf (logfile, "        iflso = %d\n", coarseIn->iflso);
   fprintf (logfile, "    tSampling = %e\n", coarseIn->tSampling);
   fprintf (logfile, "        order = %d\n", coarseIn->order);
   fprintf (logfile, "  approximant = %d\n", coarseIn->approximant);
   fprintf (logfile, "        space = %d\n", coarseIn->space);
   fprintf (logfile, "       etamin = %e\n", coarseIn->etamin);
   if (coarseIn->massRange==MinComponentMassMaxTotalMass)
   {
   fprintf (logfile, "    massRange = MinComponentMassMaxTotalMass\n");
   }
   else
   {
   fprintf (logfile, "    massRange = MinMaxComponentMass\n");
   }
   fflush(logfile);

   DETATCHSTATUSPTR(status);
   RETURN(status);
}
/*  <lalVerbatim file="WriteInspiralTemplateListCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{WriteInspiralTemplateList.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{WriteInspiralTemplateListCP}
\idx{WriteInspiralTemplateList()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{WriteInspiralTemplateListCV}}

</lalLaTeX>  */

/*  <lalVerbatim file="WriteInspiralTemplateListCP"> */
void 
WriteInspiralTemplateList 
   (
   LALStatus *status, 
   FILE *logfile, 
   INT4 nlist, 
   InspiralTemplateList *list
   )
{  /*  </lalVerbatim>  */
   INT4 i;

   INITSTATUS (status, "WriteInspiralTemplateList", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);
   ASSERT (logfile,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (list,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (nlist>0,  status, LALNOISEMODELSH_ESIZE, LALNOISEMODELSH_MSGESIZE);

   fprintf(logfile, "___________________________________________________________________\n\n"); 

   fprintf(logfile, "  Number of templates in the bank=%d\n\n", nlist);

   fprintf(logfile, "   No.    mass1       mass2      chirptime0  chirptime2  chirptime3\n");
   fprintf(logfile, "___________________________________________________________________\n"); 
	
   for (i=0; i<nlist; i++)
   {
	   fprintf(logfile, "%5d %e %e %e %e %e\n", i, list[i].params.mass1, list[i].params.mass2,
			   list[i].params.t0, list[i].params.t2, list[i].params.t3);
   }
   fprintf(logfile, "___________________________________________________________________\n"); 
   fflush(logfile);

   DETATCHSTATUSPTR(status);
   RETURN(status);

}
/*  <lalVerbatim file="WriteRandomInspiralSignalParamsCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{WriteRandomInspiralSignalParams.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{WriteRandomInspiralSignalParamsCP}
\idx{WriteRandomInspiralSignalParams()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{WriteRandomInspiralSignalParamsCV}}

</lalLaTeX>  */

/*  <lalVerbatim file="WriteRandomInspiralSignalParamsCP"> */
void 
WriteRandomInspiralSignalParams
   (
   LALStatus              *status, 
   FILE                   *logfile,
   RandomInspiralSignalIn *randIn
   )
{  /*  </lalVerbatim>  */

   INITSTATUS (status, "WriteRandomInspiralSignalParams", INSPIRALSEARCHC);
   ATTATCHSTATUSPTR(status);

   ASSERT (randIn,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);
   ASSERT (logfile,  status, LALNOISEMODELSH_ENULL, LALNOISEMODELSH_MSGENULL);

   fprintf(logfile, "____________________________________________________\n\n"); 
   fprintf (logfile,"                 RandomInspiralSignalParams         \n");
   fprintf(logfile, "____________________________________________________\n\n"); 

   fprintf (logfile, "       simulation type = %d\n", randIn->type);
   fprintf (logfile, "             SignalAmp = %e\n", randIn->SignalAmp);
   fprintf (logfile, "              NoiseAmp = %e\n", randIn->NoiseAmp);
   fprintf (logfile, "                 useed = %d\n", randIn->useed);
   fprintf (logfile, "                  mMin = %e\n", randIn->mMin);
   fprintf (logfile, "                  MMax = %e\n", randIn->MMax);
   fprintf (logfile, "                etaMin = %e\n", randIn->etaMin);

   fprintf (logfile, "       param.startTime = %e\n", randIn->param.startTime);
   fprintf (logfile, "      param.startPhase = %e\n", randIn->param.startPhase);
   fprintf (logfile, "       param.nStartPad = %d\n", randIn->param.nStartPad);
   fprintf (logfile, "         param.nEndPad = %d\n", randIn->param.nEndPad);
   fprintf (logfile, " param.signalAmplitude = %e\n", randIn->param.signalAmplitude);
   fprintf (logfile, "            param.ieta = %d\n", randIn->param.ieta);
   fprintf (logfile, "          param.fLower = %e\n", randIn->param.fLower);
   fprintf (logfile, "         param.fCutoff = %e\n", randIn->param.fCutoff);
   fprintf (logfile, "       param.tSampling = %e\n", randIn->param.tSampling);
   fprintf (logfile, "           param.order = %d\n", randIn->param.order);
   fprintf (logfile, "     param.approximant = %d\n", randIn->param.approximant);
   switch(randIn->param.massChoice)
   {
	   case m1Andm2:
		   fprintf (logfile, "      param.massChoice = m1Andm2\n");
		   break;
	   case totalMassAndEta:
		   fprintf (logfile, "      param.massChoice = totalMassAndEta\n");
		   break;
	   case t02:
		   fprintf (logfile, "      param.massChoice = t02\n");
		   break;
	   case t03:
		   fprintf (logfile, "      param.massChoice = t03\n");
		   break;
	   case t01:
	   case t04:
	   case totalMassAndMu:
		   fprintf (logfile, "      Undefined param.massChoice\n");
		   break;
   } 

   fflush(logfile);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
/*  <lalVerbatim file="WriteToLALProcessTableCV">
Author: Sathyaprakash, B.S.
</lalVerbatim>  */

/*  <lalLaTeX>

\subsection{Module \texttt{WriteToLALProcessTable.c}}

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{WriteToLALProcessTableCP}
\idx{WriteToLALProcessTable()}

\subsubsection*{Description}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{WriteToLALProcessTableCV}}

</lalLaTeX>  */

/*  <lalVerbatim file="WriteToLALProcessTable"> */
void
WriteToLALProcessTable (LALProcessTable procTable, LALDBConnect *conn) 
{  /*  </lalVerbatim>  */


	MYSQL mysql;
	char string[4096];

	sprintf(string, "insert into process values (%d,'%s','%s','%s',%d,'%s',%d,'%s','%s',%d,'%s','%s',%d,'%s','%s',%d,'%s')",
	procTable.creator_db,
	procTable.program,
	procTable.version,
	procTable.cvs_repository, 
	procTable.cvs_entry_time,
	procTable.comment,
	procTable.is_online,
	procTable.node,
	procTable.username,
	procTable.unix_procid, 
	procTable.start_time,
	procTable.end_time, 
	procTable.jobid, 
	procTable.domain,
	procTable.process_id,
	procTable.param_set,
	procTable.ifos);
       
        mysql_init(&mysql);
        if(!(mysql_real_connect(&mysql,conn->host,conn->user,conn->password,
		conn->database,0,NULL,0)))
            fprintf(stderr, "error writing to database, code: %s\n",mysql_error(&mysql));


	if (mysql_query(&mysql,string))
            fprintf(stderr, "error writing to database, code: %s\n",mysql_error(&mysql));


        conn->insert_id = mysql_insert_id(&mysql);

	mysql_close(&mysql);
	return;
}
#endif
