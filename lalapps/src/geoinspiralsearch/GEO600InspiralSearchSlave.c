/******************************** <lalVerbatim file="GEO600InspiralSearchSlaveCV">
Author: Sathyaprakash, B. S.
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{GEO600InspiralSearchSlave.c}}
\label{ss:GEO600InspiralSearchSlave.c}

LAL code for inspiral searches.

\subsubsection*{Usage}
\begin{verbatim}
GEO600InspiralSearchSlave
   (
   LALStatus *status
   )
\end{verbatim}

\subsubsection*{Description}

This is the slave code to filter the GEO600 data 
through a template bank.
\subsubsection*{Exit codes}
\input{GEO600InspiralSearchSlaveCE}

\subsubsection*{Uses}
This code directly uses the following functions:
\begin{verbatim}
INITSTATUS
ATTATCHSTATUSPTR
CHECKSTATUSPTR
DETATCHSTATUSPTR
RETURN

LALMPIRecvUINT2
LALMPISendUINT4
LALMPIRecvUINT4
LALMPISendINT4 
LALMPIRecvINT4
LALMPIRecvREAL8
LALMPIRecvREAL8Vector 
LALMPIRecvREAL4Vector 

MPI_Send
MPI_Recv
MPI_Comm_rank

LALDCreateVector 
LALSCreateVector 
LALCreateForwardRealFFTPlan
LALCreateReverseRealFFTPlan

LALREAL4VectorFFT

LALInspiralFindEvents 
LALInspiralFindEventsCluster
LALInspiralFindLoudestEvent

LALMalloc
LALRealloc
LALFree
LALDestroyRealFFTPlan
LALDestroyVector
LALSDestroyVector
LALDDestroyVector
LALDestroyRealFFTPlan
LALDestroyRealFFTPlan

\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GEO600InspiralSearchSlaveCV}}
******************************************************* </lalLaTeX> */

/***************************** <lalErrTable file="GEO600InspiralSearchSlaveCE"> */
/***************************** </lalErrTable> */

#include <config.h>
#if defined HAVE_GSL_GSL_FFT_REAL_H && defined HAVE_LIBGSL \
  && defined HAVE_MYSQL_H && defined HAVE_LIBMYSQLCLIENT \
  && defined HAVE_FRAMEL_H && defined HAVE_LIBFRAME \
  && defined HAVE_MPI_H

#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/Comm.h>
#include <mpi.h>
#include <GEO600InspiralSearch.h>
#include <gsl/gsl_fft_real.h>

#define WORKTAG 1
 
/*
 * Slave code
 */

NRCSID (INSPIRALSEARCHSLAVEC, "$Id$");

/* <lalVerbatim file="GEO600InspiralSearchSlaveCP"> */
void 
GEO600InspiralSearchSlave 
   (
   LALStatus *status
   )
{ /* </lalVerbatim> */
   UINT4 dataSetDuration;
   INT4 i;
   INT4 j;
   INT4 nlist;
   INT4 eventsListLength;
   INT4 tmpltListLength;
   INT4 sourceM=0;
   INT4 last=0;
   INT4 slaveNum;
   INT4 initMe;
   MPI_Status mpistatus;

   InspiralTemplateList *list;         /* parameter structure of templates */
   LALWindowParams      winParams;     /* Parameters of the window */

   InspiralFindEventsIn eventsin;      /* parameter structure for finding events */
   RunningAveragePSDIn  pwrSpecIn;     /* parameter struct used in computing the PSD */
   InspiralSearchParams searchParams;  /* parameters of the search */

   REAL8Vector *buffer1=NULL;          /* Vectors to hold data set */
   REAL8Vector *pwrSpec=NULL;
   REAL4Vector *dataset=NULL;
   RealFFTPlan *fwdp=NULL;
   RealFFTPlan *revp=NULL;
   REAL4Vector *winVect=NULL;          /* Vector to hold window used in computing psd */


   UINT4 dataSetLength;
   UINT4 pwrSpecLength;
   UINT4 searchParamsLength;
   UINT4 currentGPSTime;
   UINT4 samplingRate;
   int sizeint;


   INITSTATUS (status, "GEO600InspiralSearchSlave", INSPIRALSEARCHSLAVEC);
   ATTATCHSTATUSPTR(status);
/* 
 * Find out who I am
 */
   MPI_Comm_rank(MPI_COMM_WORLD, &slaveNum);

/*
 * is master ready to initialise me?
 */
   LALMPIRecvINT4 (status->statusPtr, &initMe, sourceM, MPI_COMM_WORLD); 
   CHECKSTATUSPTR(status);

   if (initMe) 
   {
      fprintf(stderr, "slave %d terminating ...\n", slaveNum);
      DETATCHSTATUSPTR(status);
      RETURN(status);
   }

/* 
 * read the parameters of search and threshold and whether to display info
 */
              sizeint = sizeof(INT4);
   searchParamsLength = sizeof(InspiralSearchParams);

   /*
    * Read the parameters of search 
    */

   ReadInspiralSearchParams(status->statusPtr, &searchParams);
   CHECKSTATUSPTR(status);
   ReadInspiralFindEventsIn (status->statusPtr, &eventsin);
   CHECKSTATUSPTR(status);

/* 
 * receive the sampling rate of the channel
 */
   LALMPIRecvUINT4(status->statusPtr, &samplingRate, sourceM, MPI_COMM_WORLD);
   CHECKSTATUSPTR(status);
/*
 * receive the number of templates in the bank
 */

   LALMPIRecvINT4 (status->statusPtr, &nlist, sourceM, MPI_COMM_WORLD); 
   CHECKSTATUSPTR(status);

/*
 * receive the length of dataset 
 */
   LALMPIRecvINT4 (status->statusPtr, &dataSetLength, sourceM, MPI_COMM_WORLD); 
   CHECKSTATUSPTR(status);
   LALMPIRecvINT4 (status->statusPtr, &pwrSpecLength, sourceM, MPI_COMM_WORLD); 
   CHECKSTATUSPTR(status);

/*
 * allocate memory for the template bank parameters
 */

   tmpltListLength = sizeof(InspiralTemplateList);
   list = (InspiralTemplateList *) LALMalloc(tmpltListLength * nlist);

   fprintf(stderr, "slave %d starting %d %d %d\n", slaveNum, dataSetLength, pwrSpecLength, nlist);
/*
 * receive template-bank parameters one-by-one
 */

   for (i=0; i<nlist; i++)
   {
      
	   MPI_Recv(&list[i], tmpltListLength, MPI_BYTE, 
			   sourceM, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus);
   }

/*
 * create vectors and estimate FFTW plans 
 */

   LALDCreateVector (status->statusPtr, &pwrSpec, pwrSpecLength);
   CHECKSTATUSPTR(status);
   LALSCreateVector (status->statusPtr, &dataset, dataSetLength);
   CHECKSTATUSPTR(status);
   LALCreateVector(status->statusPtr, &winVect, dataSetLength);
   CHECKSTATUSPTR(status);
   LALCreateForwardRealFFTPlan(status->statusPtr, &fwdp, dataSetLength, 0);
   CHECKSTATUSPTR(status);
   LALCreateReverseRealFFTPlan(status->statusPtr, &revp, dataSetLength, 0);
   CHECKSTATUSPTR(status);
   LALDCreateVector(status->statusPtr, &buffer1, dataSetLength);
   CHECKSTATUSPTR(status);
   /*
    * Create window to be used for PSD computation 
    */

   winParams.length = dataSetLength;
   winParams.type = Hamming;
   LALWindow (status->statusPtr, winVect, &winParams);
   CHECKSTATUSPTR(status);


/*
 * cast the structure needed by LALInspiralFindEvents 
 */

   eventsin.nBegin = dataSetLength / searchParams.paddingFactor;
   eventsin.nEnd = dataSetLength / searchParams.paddingFactor;
   eventsin.fwdp = fwdp;
   eventsin.revp = revp;
   eventsin.psd = *pwrSpec;
   eventsin.signal = *dataset;
   eventsListLength = sizeof(InspiralEventsList);
   tmpltListLength = sizeof(InspiralTemplateList);
      
   pwrSpecIn.norm = winParams.sumofsquares; 
   pwrSpecIn.fwdp = fwdp; 
   pwrSpecIn.window = winVect; 
   dataSetDuration = dataSetLength / samplingRate;
   searchParams.samplingRate = (UINT4) samplingRate;
   searchParams.dataSetLength = dataSetLength;
   searchParams.dataSetDuration = dataSetDuration;
   searchParams.df = 1. / (REAL8) dataSetDuration;

   do
   {
      static InspiralEventsList currentEvent, loudestEvent;
      InspiralEventsList *totEventsList=NULL;
      INT4 nEvents=0, first=1, totEvents=0, inputPSD=0, inputData=0;
/*
 * (9) and (16) receive dataset and psd 
 */
      LALMPIRecvUINT4(status->statusPtr, &currentGPSTime, sourceM, MPI_COMM_WORLD);
      CHECKSTATUSPTR(status);
	   
      searchParams.currentGPSTime = currentGPSTime;
      pwrSpecIn.searchParams = searchParams; 
   
      RunningAverageDoublePSD(status->statusPtr, pwrSpec, &inputPSD, &pwrSpecIn);
      /* LALNoiseSpectralDensity (status->statusPtr, pwrSpec, noisemodel, df); */
      CHECKSTATUSPTR(status);
   
      if (inputPSD) 
      {
	      fprintf(stderr, "problem computing running average PSD of data\n");
	      exit(0);
      }
      /* 
       * get/generate dataset
       */
   
      inputData = GetChannelDoubleData(currentGPSTime, dataSetDuration, 
		      searchParams.channelName, buffer1->data, dataSetLength); 
      if (inputData) 
      {
	      fprintf(stderr, "problem finding specified data\n");
	      exit(0);
      } 

   
      /*
       * and compute its FFT 
       * Need to use gsl to do this in double precision 
       */
   
      gsl_fft_real_radix2_transform(buffer1->data, 1, buffer1->length);
   
      /* Now need to convert the double precision result of this
       * fft to single precision for the rest of the code to 
       * work with 
       */

      for (i=0; i<(INT4)buffer1->length; i++)
      {
	      dataset->data[i] = (REAL4)buffer1->data[i];
      }
   
      /*
       * Old Fourier transform used:
       LALREAL4VectorFFT(status->statusPtr, dataset, buffer1, fwdp);
       CHECKSTATUSPTR(status);
       */

      eventsin.currentGPSTime = currentGPSTime;
      switch (searchParams.findEventsType)
      {
	      case LoudestEvent:

		      /* 
		       * for zeroeth template in the list ...
		       */
      
		      eventsin.param = loudestEvent.param = list[0].params;
		      LALInspiralFindLoudestEvent(status->statusPtr, &nEvents, &loudestEvent, &eventsin);

		      /* 
		       * for each template in the list ...
		       */
      
		      for (j=1; j<nlist; j++) 
		      { 
			      eventsin.param = list[j].params;
			      /* 
			       * filter and find the loudest event
			       */
			      LALInspiralFindLoudestEvent(status->statusPtr, &nEvents, &currentEvent, &eventsin);
			      totEvents += nEvents;
			      /*
			      fprintf (stderr, "slaveNum=%d nEvents=%d max=%e\n", slaveNum, nEvents, currentEvent.snr);
			      */
			      if (currentEvent.snr > loudestEvent.snr)
			      {
				      loudestEvent = currentEvent;
			      }
		      }
      
		      /* send master number of events */ 
		      fprintf (stderr, "slaveNum=%d nEvents=%d max=%e chisq=%e\n", slaveNum, totEvents, loudestEvent.snr, loudestEvent.chisq);
      
		      nEvents=1;
		      MPI_Send(&nEvents, sizeint, MPI_BYTE, sourceM, WORKTAG, MPI_COMM_WORLD);
		      MPI_Send(&loudestEvent, eventsListLength, MPI_BYTE, sourceM, WORKTAG, MPI_COMM_WORLD);
		      break;

	      case AllTemplates:
				      
		      totEventsList = (InspiralEventsList *) LALMalloc(eventsListLength * nlist);
		      for (j=0; j<nlist; j++)
		      {
			      eventsin.param = list[j].params;
			      /* 
			       * filter and find the loudest event in template j
			       */
			      LALInspiralFindLoudestEvent(status->statusPtr, &nEvents, &totEventsList[j], &eventsin);
		      }
      
		      /* send master number of events */ 

		      MPI_Send(&nlist, sizeint, MPI_BYTE, sourceM, WORKTAG, MPI_COMM_WORLD);

		      /* 
		       * send master param struct of each event
		       */

		      for (i=0; i<nlist; i++) 
		      { 
			      /*
			      fprintf(stderr, "eventID=%d, snr=%e\n", i, totEventsList[i].max);
			      */
			      MPI_Send(&totEventsList[i], eventsListLength, MPI_BYTE, sourceM, WORKTAG, MPI_COMM_WORLD);
		      }
		      LALFree(totEventsList);
		      break;
	      case AllEvents:
		            
		      for (j=0; j<nlist; j++)
		      {
			      InspiralEventsList *eventsList=NULL;
			      INT4 NEvents = 0;
			      eventsin.param = list[j].params;
										          
			      /*
			       * filter and find number of events crossing eventsin.Threshold
			       */
			      /*
			      LALInspiralFindEvents(status->statusPtr, &NEvents, &eventsList, &eventsin);
			      */
			      LALInspiralFindEventsCluster(status->statusPtr, &NEvents, &eventsList, &eventsin);

			      if (NEvents)
			      {
				      if (first)
				      {
					      totEventsList = (InspiralEventsList *) 
						      LALMalloc(eventsListLength * NEvents);
					      first = 0;
				      } 
				      else 
				      {
					      totEventsList = (InspiralEventsList *) 
						      LALRealloc(totEventsList, eventsListLength *(NEvents+totEvents));
				      }
				      for (i=0; i<NEvents; i++)
				      {
					      totEventsList[totEvents+i] = eventsList[i];
				      }
				      LALFree(eventsList);
				      eventsList=NULL;
			      }
			      totEvents+=NEvents;
		      }
      
		      /* send master number of events */ 

		      nEvents = totEvents;
		      fprintf(stderr, "slaveNum=%d, NEvents=%d\n", slaveNum, totEvents);
		      MPI_Send(&nEvents, sizeint, MPI_BYTE, sourceM, WORKTAG, MPI_COMM_WORLD);

		      /* 
		       * send master param struct of each event
		       */

		      for (i=0; i<nEvents; i++) 
		      { 
			      /*
			      fprintf(stderr, "eventID=%d, snr=%e\n", i, totEventsList[i].max);
			      */
			      MPI_Send(&totEventsList[i], eventsListLength, MPI_BYTE, sourceM, WORKTAG, MPI_COMM_WORLD);
		      }
		      if (totEvents)
		      {
			      LALFree(totEventsList);
		      }
		      break;
      }

      LALMPISendUINT4(status->statusPtr, &currentGPSTime, sourceM, MPI_COMM_WORLD);
      CHECKSTATUSPTR(status);

      /* 
       * are we receiving anymore data?
       */

      LALMPIRecvINT4 (status->statusPtr, &last, sourceM, MPI_COMM_WORLD); 

   } while (!last);

   LALDestroyVector(status->statusPtr, &winVect);
   CHECKSTATUSPTR(status);
   LALDDestroyVector(status->statusPtr, &buffer1);
   CHECKSTATUSPTR(status);
   LALSDestroyVector(status->statusPtr, &dataset);
   CHECKSTATUSPTR(status);
   LALDDestroyVector(status->statusPtr, &pwrSpec);
   CHECKSTATUSPTR(status);
   LALDestroyRealFFTPlan(status->statusPtr,&fwdp);   
   CHECKSTATUSPTR(status);
   LALDestroyRealFFTPlan(status->statusPtr,&revp);
   LALFree(list);
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
#endif
