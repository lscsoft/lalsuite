/******************************** <lalVerbatim file="GEO600InspiralSearchDoubleMasterDBCV">
Author: Sathyaprakash, B. S.
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{GEO600InspiralSearchDoubleMasterDB.c}}
\label{ss:GEO600InspiralSearchDoubleMasterDB.c}

LAL code for inspiral searches.

\subsubsection*{Usage}
\begin{verbatim}
GEO600InspiralSearchDoubleMasterDB
   (
   LALStatus *status
   );
\end{verbatim}

\subsubsection*{Description}

This is the master code that filters the GEO600 
data through a template bank.

\subsubsection*{Exit codes}
\input{GEO600InspiralSearchDoubleMasterDBCE}

\subsubsection*{Uses}
This code directly uses the following functions:

\begin{verbatim}
INITSTATUS
ATTATCHSTATUSPTR
CHECKSTATUSPTR
DETATCHSTATUSPTR
RETURN

LALMPISendUINT2
LALMPISendUINT4
LALMPIRecvUINT4
LALMPISendINT4 
LALMPIRecvINT4
LALMPISendREAL8
LALMPISendREAL8Vector 
LALMPISendREAL4Vector 

MPI_Send
MPI_Recv
MPI_Comm_rank

LALDCreateVector 
LALSCreateVector 
LALCreateForwardRealFFTPlan
LALCreateReverseRealFFTPlan

LALREAL4VectorFFT

LALInspiralFindEvents 

LALMalloc
LALRealloc
LALFree
LALDestroyRealFFTPlan
LALDestroyVector
LALSDestroyVector
LALDDestroyVector
LALDestroyRealFFTPlan
LALDestroyRealFFTPlan

LALInspiralCreateCoarseBank
LALWindow 

LongestTemplateInBank
RunningAverageDoublePSD

GetSampleRate
ReadInspiralFindEventsIn

WriteInspiralFindEventsIn 
WriteInspiralSearchParams
WriteInspiralTemplateList 
WriteInspiralTemplateBankParams
WriteInspiralEvents
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GEO600InspiralSearchDoubleMasterDBCV}}
******************************************************* </lalLaTeX> */

/***************************** <lalErrTable file="GEO600InspiralSearchDoubleMasterDBCE"> */
/***************************** </lalErrTable> */

#include <config.h>
#if defined HAVE_GSL_GSL_FFT_REAL_H && defined HAVE_LIBGSL \
 && defined HAVE_MYSQL_H && defined HAVE_LIBMYSQLCLIENT \
 && defined HAVE_FRAMEL_H && defined HAVE_LIBFRAME \
 && defined HAVE_MPI_H

#include <lal/LALInspiralBank.h>
#include <lal/LALNoiseModels.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/Comm.h>
#include <mpi.h>
#include <sys/stat.h>
#include <GEO600InspiralSearch.h>
#include <gsl/gsl_fft_real.h>

#define maxNumEvents 100000
#define WORKTAG 2
 
/*
 * Master code
 */

NRCSID (INSPIRALSEARCHMASTERC, "$Id$");

/* <lalVerbatim file="GEO600InspiralSearchDoubleMasterDBCP"> */
void 
GEO600InspiralSearchDoubleMasterDB
   (
   LALStatus *status
   )
{ /* </lalVerbatim> */

/* 
 * Database structures
 */
   
   LALDBConnect conn;
   LALProcessTable procTable;

   MPI_Status mpistatus;

   UINT4 nBegin, nEnd;
   UINT4 slaveNum;                     /* slave ID */
   UINT4 nActiveTasks=1;               /* in the beginning there is only master no slaves */
   UINT4 nTasks;                       /* number of tasks running = number of slaves-1 */
   UINT4 eventListLength;              /* length, in bytes, of struct InspiralEventsList */
   UINT4 tmpltListLength;              /* length, in bytes, of struct InspiralTemplateList */
   UINT4 currentGPSTime;
   UINT4 slaveGPSTime;

   INT4 nlist;                         /* number of templates in coarse bank */
   INT4 j, k, nby2;
   INT4 last=0;
/* 
 * ingeter variables to keep track of success/failure of certain operations 
 */

   INT4 inputPSD=0;
   INT4 q, i;
/* 
 * Various parameter structures
 */

   InspiralSearchParams   searchParams;/* parameters of the search */
   InspiralCoarseBankIn   coarseIn;    /* parameter stucture to find coarse bank templates */
   InspiralTemplateList   *list=NULL;  /* list of templates */
   InspiralFindEventsIn   eventsin;    /* parameter structure for finding events */
   LALWindowParams        winParams;   /* Parameters of the window */
   RunningAveragePSDIn    pwrSpecIn;   /* parameter struct used in computing the PSD */
   InspiralEventsList     *eventList=NULL;  /* param struct to hold events */
/* 
 * Various lengths of arrays
 */
   UINT4 pwrSpecLength;                /* Lengths of template, data and S_h(f) */
   UINT4 dataSetLength;
   UINT4 templateLength;
   UINT4 dataSetDuration;
   UINT4 templateDuration;             /* Duration of data and template  */
   UINT4 totPSDDuration;               /* Duration for which data sets are needed to compute PSD */
   UINT4 PSDStartTime;                 /* Starting time of the data set to compute PSD at currentGPSTime */
   UINT4 overlapDuration;              /* Duration for which consecutive data sets must overlap */

   REAL8Vector *pwrSpec=NULL;          /* Vector to hold S_h(f) */
   REAL8Vector *buffer1=NULL;          /* Vectors to hold data set */
   REAL4Vector *dataset=NULL;
   REAL4Vector *winVect=NULL;          /* Vector to hold window used in computing psd */

                                       /* Vectors to compute forward and inverse FFT */
   RealFFTPlan *fwdp=NULL;
   RealFFTPlan *revp=NULL;             
   /*
   */

   REAL8 df;                           /* Resolution in time and frequency spaces */
   REAL8 dt;

   static float samplingRate;          /* Sampling rate of the channel */

   FILE         *logfile;              /* File in which all log will be written */

   static char  EvtFile[512];
   static char  LogFile[512];

   int sizeint;

   void *noisemodel=LALLIGOIPsd;
	
   char *channelName = "G1:DER_H_QUALITY";
   REAL8Vector *quality=NULL;


   INITSTATUS (status, "GEO600InspiralSearchMaster", INSPIRALSEARCHMASTERC);
   ATTATCHSTATUSPTR(status);

   SetUpDBConnectionInfo(&conn);
   SetUpProcessTableInfo(&procTable);
   WriteToLALProcessTable(procTable, &conn);

   MPI_Comm_size(MPI_COMM_WORLD, &nTasks);

   /*
    * length in bytes of important structures
    */
   tmpltListLength = sizeof(InspiralTemplateList);
   eventListLength = sizeof(InspiralEventsList);
           sizeint = sizeof(INT4);


   /*
    * allocate memory for the template bank parameters
    */
         
   eventList = (InspiralEventsList *) LALMalloc(eventListLength * maxNumEvents);

   if (eventList==NULL) 
   {
	   fprintf(stderr, "memory allocation failure for eventList in master\n");
	   goto failAnalysis;
   }

   /*
    * Read the parameters of search 
    */

   ReadInspiralSearchParams(status->statusPtr, &searchParams);
   CHECKSTATUSPTR(status);
   ReadInspiralFindEventsIn (status->statusPtr, &eventsin);
   CHECKSTATUSPTR(status);

   /*
    * Open files to read and write 
    */

   sprintf(EvtFile,"GEO600InspiralSearch%d.evt",searchParams.currentGPSTime);
   sprintf(LogFile,"GEO600InspiralSearch%d.log",searchParams.currentGPSTime);

   logfile = fopen(LogFile,"w");

   if (logfile == NULL)
   {
      fprintf(stderr, "Could not open output file(s)\n%s and/or\n%s\n\n", EvtFile, LogFile);
      goto failAnalysis;
   }


   samplingRate = (float) searchParams.samplingRate;

   /*
    * Read the parameters of the template bank
    */

   ReadInspiralTemplateBankParams(status->statusPtr, &coarseIn);
   CHECKSTATUSPTR(status);
   coarseIn.tSampling = (REAL8) samplingRate;

   /*
    * Get the length of the longest template 
    */

   LALInspiralLongestTemplateInBank(status->statusPtr, &dataSetLength, &coarseIn);
   CHECKSTATUSPTR(status);

   /*
    * Set array lengths 
    */

   templateLength = dataSetLength;
   templateDuration = templateLength / samplingRate;
   dataSetLength *= searchParams.paddingFactor;

   if (dataSetLength < (UINT4) samplingRate)  
   {
       /* We want at least 1 s of data */
       dataSetLength = (UINT4) samplingRate; 
   }

   dataSetDuration = dataSetLength / samplingRate;
   pwrSpecLength = dataSetLength/2 + 1;

   searchParams.dataSetLength = dataSetLength;
   searchParams.dataSetDuration = dataSetDuration;
   nBegin = dataSetLength / searchParams.paddingFactor;
   nEnd = dataSetLength / searchParams.paddingFactor;
   /* 
    * In filtered data we will discard chunks both at the beginning 
    * and end; hence an additional factor of 2 in overlapDuration
    */
   overlapDuration = 2 * dataSetDuration/searchParams.paddingFactor;
   /* The duration for which we need to check quality is the 
    * full span of the data segment from which the PSD is computed
    * including the current data segment, which is not included in
    * the PSD calculation
    */ 
   totPSDDuration = (searchParams.numDataSetsForPSD/2+1) * dataSetDuration + dataSetDuration;

   nby2 = dataSetLength / 2;
   dt = 1. / (REAL8) samplingRate;
   df = 1. / (REAL8) dataSetDuration;
   searchParams.df = df;
   /*
    * Estimate FFTW plans 
    */

   LALCreateForwardRealFFTPlan(status->statusPtr, &fwdp, dataSetLength, 0);
   CHECKSTATUSPTR(status);
   LALCreateReverseRealFFTPlan(status->statusPtr, &revp, dataSetLength, 0);
   CHECKSTATUSPTR(status);
   /*
   */

   /*
    * Create dataSet, buffer, window and psd vectors 
    */


   LALCreateVector(status->statusPtr, &dataset, dataSetLength);
   CHECKSTATUSPTR(status);
   LALCreateVector(status->statusPtr, &winVect, dataSetLength);
   CHECKSTATUSPTR(status);
   LALDCreateVector(status->statusPtr, &buffer1, dataSetLength);
   CHECKSTATUSPTR(status);
   LALDCreateVector(status->statusPtr, &pwrSpec, pwrSpecLength);
   CHECKSTATUSPTR(status);
   LALDCreateVector(status->statusPtr, &quality, totPSDDuration);
   CHECKSTATUSPTR(status);

   /*
    * Create window to be used for PSD computation 
    */

   winParams.length = dataSetLength;
   winParams.type = Hann;
   LALWindow (status->statusPtr, winVect, &winParams);
   CHECKSTATUSPTR(status);


   /*
    * Get the template bank 
    */

   memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries));
   LALDCreateVector (status->statusPtr, &(coarseIn.shf.data), pwrSpecLength);
   CHECKSTATUSPTR(status);
   coarseIn.shf.f0 = 0.;
   coarseIn.shf.deltaF = coarseIn.tSampling / (REAL8) coarseIn.shf.data->length;
   pwrSpecIn.fwdp = fwdp; 
   pwrSpecIn.revp = revp; 
   /*
   */
   pwrSpecIn.window = winVect; 
   pwrSpecIn.searchParams = searchParams; 
   pwrSpecIn.norm = winParams.sumofsquares; 
   RunningAverageDoublePSD(status->statusPtr, coarseIn.shf.data, &inputPSD, &pwrSpecIn);
   CHECKSTATUSPTR(status);
   /*
   LALNoiseSpectralDensity (status->statusPtr, coarseIn.shf.data, noisemodel, df);
   */

   if (eventsin.displayPSD)
   {
	   for (k=0; k<(INT4)coarseIn.shf.data->length; k++) 
	   {
		   fprintf(stdout, "%e %e\n", df*k, pow(coarseIn.shf.data->data[k],0.5));
	   }
           fprintf(stdout, "&\n");
   }

   LALInspiralCreateCoarseBank(status->statusPtr, &list, &nlist, coarseIn);
   CHECKSTATUSPTR(status);

   /*
    * Write some useful information onto logfile
    */

   WriteInspiralSearchParams(status->statusPtr, logfile, &searchParams);
   CHECKSTATUSPTR(status);
   WriteInspiralFindEventsIn (status->statusPtr, logfile, &eventsin);
   CHECKSTATUSPTR(status);
   WriteInspiralTemplateBankParams(status->statusPtr, logfile, &coarseIn);
   CHECKSTATUSPTR(status);
   fprintf(logfile, "____________________________________________________\n\n"); 

   fprintf(logfile, "              Template bank parameters              \n"); 
   fprintf(logfile, "____________________________________________________\n"); 
   fprintf(logfile, "        templateLength = %d\n", templateLength); 
   fprintf(logfile, "      templateDuration = %d\n", templateDuration); 
   fprintf(logfile, "         dataSetLength = %d\n", dataSetLength); 
   fprintf(logfile, "       dataSetDuration = %d\n", dataSetDuration); 
   fprintf(logfile, "       totPSDDuration = %d\n", totPSDDuration); 
   fprintf(logfile, "____________________________________________________\n"); 
   WriteInspiralTemplateList (status->statusPtr, logfile, nlist, list);
   CHECKSTATUSPTR(status);

   /*
    * if number of templates is zero we might as well terminate
    * analysis
    */
   if (nlist==0) 
   {
      fprintf(logfile, "____________________________________________________\n"); 

      fprintf(logfile, " parameter space too small: no templates to analyse\n\n");
      fprintf(logfile, "____________________________________________________\n"); 
      goto terminateAnalysis;
   }

   /* 
    * if the currentGPSTime is greater than endGPSTime we stop 
    */

   if (searchParams.currentGPSTime > searchParams.endGPSTime) 
   {
           fprintf(logfile, "endGPSTime=%d < startGPSTime=%d\n",searchParams.endGPSTime, searchParams.currentGPSTime);
	   last = 1;
	   goto endAnalysis;
   }

   currentGPSTime = searchParams.currentGPSTime;

   /* 
    * get a locked streatch
    */

   q=0;
   while (!q && currentGPSTime <= searchParams.endGPSTime)
   {
	   PSDStartTime = currentGPSTime - dataSetDuration/2*(searchParams.numDataSetsForPSD/2+1);
	   GetChannelDoubleData(PSDStartTime, totPSDDuration, channelName, quality->data, totPSDDuration); 
	   q = 1; for (i=0; i<(INT4)totPSDDuration; i++) if ( (q=quality->data[i]) == 0 ) break;
	   fprintf (logfile, "quality=%d at GPS=%d PSD start Time=%d\n", q, currentGPSTime, PSDStartTime);
	   if (q) break;
	   searchParams.currentGPSTime += (dataSetDuration-overlapDuration);
	   currentGPSTime = searchParams.currentGPSTime;

   }
	   
   if (currentGPSTime > searchParams.endGPSTime) 
   {
           fprintf(logfile, "endGPSTime=%d<startGPSTime=%d\n",searchParams.endGPSTime, 
			   currentGPSTime);
	   last = 1;
	   goto endAnalysis;
   }

	   
   for (slaveNum = 1; slaveNum < nTasks; ++slaveNum) 
   { 
	   if (last)
	   {
		   goto endAnalysis;
	   }

	   nActiveTasks++;

	   /* 
	    * send slaves the following initialization info 
	    * parameters of the search 
	    */ 

           LALMPISendINT4(status->statusPtr, &last, slaveNum, MPI_COMM_WORLD);

           /*
            * the sampling rate of the channel 
            */

           LALMPISendUINT4(status->statusPtr, &(searchParams.samplingRate), 
			   slaveNum, MPI_COMM_WORLD);
           CHECKSTATUSPTR(status);

           /*
            * number of templates in the bank
            */

           LALMPISendINT4 (status->statusPtr, &nlist, slaveNum, MPI_COMM_WORLD); 
           CHECKSTATUSPTR(status);

           /*
            * send the length of the dataset and psd,
            */

           LALMPISendINT4 (status->statusPtr, &dataSetLength, slaveNum, MPI_COMM_WORLD); 
           CHECKSTATUSPTR(status);
           LALMPISendINT4 (status->statusPtr, &pwrSpecLength, slaveNum, MPI_COMM_WORLD); 
           CHECKSTATUSPTR(status);

           /*
            * and the template bank parameters one-by-one,
            */

           for (j=0; j<nlist; j++) 
	   { 
		   MPI_Send(&list[j], tmpltListLength, MPI_BYTE, slaveNum, 
				   WORKTAG, MPI_COMM_WORLD); 
	   }

           /*
            * and current GPS time.
            */

	   fprintf(stderr, "Sending slave %d data from %d GPS sec\n",slaveNum,currentGPSTime);
	   fprintf(logfile, "Sending slave %d data from %d GPS sec\n",slaveNum,currentGPSTime);

           LALMPISendUINT4(status->statusPtr, &currentGPSTime, slaveNum, MPI_COMM_WORLD);
           CHECKSTATUSPTR(status);

	   /* 
	    * increment clock
	    */

           searchParams.currentGPSTime += (dataSetDuration-overlapDuration);
            if (searchParams.currentGPSTime > searchParams.endGPSTime) 
            {
                    fprintf(logfile, "Analysis stopped after segment GPS time %d \n", currentGPSTime);
	            last = 1;
	            goto endAnalysis;
            }
           currentGPSTime = searchParams.currentGPSTime;

	   /* 
	    * read a new dataset and compute its FFT 
	    * Need to use gsl to do this in double precision
	    */ 

	   q=0;
	   while (!q && currentGPSTime <= searchParams.endGPSTime)
	   {
		   PSDStartTime = currentGPSTime - dataSetDuration/2*(searchParams.numDataSetsForPSD/2+1);
		   GetChannelDoubleData(PSDStartTime, totPSDDuration, channelName, quality->data, totPSDDuration); 
		   q = 1; for (i=0; i<(INT4)totPSDDuration; i++) if ( (q=quality->data[i]) == 0 ) break;
		   fprintf (logfile, "quality=%d at GPS=%d PSD start Time=%d\n", q, currentGPSTime, PSDStartTime);
		   if (q) break;
		   searchParams.currentGPSTime += (dataSetDuration-overlapDuration);
		   currentGPSTime = searchParams.currentGPSTime;
	   }

	   if (currentGPSTime > searchParams.endGPSTime) last = 1;
   }

   /* 
    * here follows the continuous mode
    */

   while (!last)
   { 
	   INT4 nEvents=0;                      /* Number of events */ 

	   /* 
	    * close the log and events files so that they have updated info at
	    * any point in the run - a crude check-pointing; open again to write. 
	    */

	   fflush(logfile);
   
	   /*
	    * receive number of events from any slave ready to send
	    */
      
	   MPI_Recv(&nEvents, sizeint, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus);
           slaveNum = mpistatus.MPI_SOURCE;
      
	   /*
	    * receive param struct of each event
	    */ 

	   for (j=0; j<nEvents; j++)
	   { 
		   MPI_Recv(&eventList[j], eventListLength, MPI_BYTE, slaveNum, 
				   MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus); 
	   }

	   LALMPIRecvUINT4(status->statusPtr, &slaveGPSTime, slaveNum, MPI_COMM_WORLD);
	   CHECKSTATUSPTR(status);
	
	   /*
	    * tell the same slave if a new dataset should be processed
	    */
      
	   LALMPISendINT4 (status->statusPtr, &last, slaveNum, MPI_COMM_WORLD);
	   CHECKSTATUSPTR(status);

	   /*
	    * send the gps second of the next dataset to the same slave
	    */
      
	   LALMPISendUINT4(status->statusPtr, &currentGPSTime, slaveNum, MPI_COMM_WORLD);
	   CHECKSTATUSPTR(status);

           fprintf(stderr, "from slave %d, in data GPSTime=%d, number of events=%d\n", 
           		slaveNum, slaveGPSTime, nEvents); 
           fprintf(logfile, "from slave %d, in data GPSTime=%d, number of events=%d\n", 
           		slaveNum, slaveGPSTime, nEvents); 

	   fprintf(stderr, "Sending slave %d data from %d GPS sec\n",slaveNum, currentGPSTime);

	   fprintf(logfile, "Sending slave %d data from %d GPS sec\n",slaveNum, currentGPSTime);

	   WriteInspiralEventsDB(status->statusPtr, eventList, nEvents, conn);
	   CHECKSTATUSPTR(status);
	   
	   /*
	    * increment clock; stop if we have no more datasets to analyse
	    */
      
	   searchParams.currentGPSTime += (dataSetDuration-overlapDuration);

	   if (searchParams.currentGPSTime > searchParams.endGPSTime) 
	   {
                   fprintf(logfile, "Analysis stopped after segment GPS time %d \n", currentGPSTime);
		   last = 1;
	   }
	   currentGPSTime = searchParams.currentGPSTime;
      
	   q=0;
	   while (!q && currentGPSTime <= searchParams.endGPSTime)
	   {
		   PSDStartTime = currentGPSTime - dataSetDuration/2*(searchParams.numDataSetsForPSD/2+1);
		   GetChannelDoubleData(PSDStartTime, totPSDDuration, channelName, quality->data, totPSDDuration); 
		   q = 1; for (i=0; i<(INT4)totPSDDuration; i++) if ( (q=quality->data[i]) == 0 ) break;
		   fprintf (logfile, "quality=%d at GPS=%d PSD start Time=%d\n", q, currentGPSTime, PSDStartTime);
		   if (q) break;
		   searchParams.currentGPSTime += (dataSetDuration-overlapDuration);
		   currentGPSTime = searchParams.currentGPSTime;
	   }

	   if (currentGPSTime > searchParams.endGPSTime) last = 1;
   }

   endAnalysis:
   fprintf(stderr, "#nActiveTasks=%d\n", nActiveTasks); 
   fprintf(logfile, "#nActiveTasks=%d\n", nActiveTasks); 

   /* 
    * close the log and events files so that they have updated info at
    * any point in the run - a crude check-pointing; open again to write. 
    */

    fflush(logfile);
   
   /*
    * Receive events from all the processes 
    */

   
   for (k=1; k<(INT4)nActiveTasks; k++) 
   {
   
	   INT4 nEvents=0;
   
	   /*
	    * receive number of events from slave
	    */
      
	   MPI_Recv(&nEvents, sizeint, MPI_BYTE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus);
           slaveNum = mpistatus.MPI_SOURCE;

	   /*
	    * receive param struct of each event
	    */ 
		   
	   for (j=0; j<nEvents; j++) 
	   {
		   MPI_Recv(&eventList[j], eventListLength, MPI_BYTE, 
				   slaveNum, MPI_ANY_TAG, MPI_COMM_WORLD, &mpistatus); 
	   }

	   fprintf(stderr, "Recd %d events from slave %d\n",nEvents, slaveNum);
	   fprintf(logfile, "Recd %d events from slave %d\n",nEvents, slaveNum);
	   LALMPIRecvUINT4(status->statusPtr, &slaveGPSTime, slaveNum, MPI_COMM_WORLD);
	   CHECKSTATUSPTR(status);

	   /*
           fprintf(stderr, "#from slave %d, in data GPSTime=%d, number of events=%d\n", 
           		slaveNum, slaveGPSTime, nEvents); 

	   */

	   WriteInspiralEventsDB(status->statusPtr, eventList, nEvents, conn);
	   CHECKSTATUSPTR(status);
   }

   terminateAnalysis:
   AppendLALProcessTable(conn);
   LALDDestroyVector(status->statusPtr, &quality);
   CHECKSTATUSPTR(status);
   LALDestroyVector(status->statusPtr, &winVect);
   CHECKSTATUSPTR(status);
   LALDestroyVector(status->statusPtr, &dataset);
   CHECKSTATUSPTR(status);
   LALDDestroyVector(status->statusPtr, &pwrSpec);
   CHECKSTATUSPTR(status);
   LALDDestroyVector(status->statusPtr, &buffer1);
   CHECKSTATUSPTR(status);
   LALDestroyRealFFTPlan(status->statusPtr,&fwdp);   
   CHECKSTATUSPTR(status);
   LALDestroyRealFFTPlan(status->statusPtr,&revp);
   CHECKSTATUSPTR(status);
   /*
   */
   LALDDestroyVector (status->statusPtr, &(coarseIn.shf.data));
   CHECKSTATUSPTR(status);

   if (eventList !=NULL ) LALFree(eventList);
   if (list !=NULL ) LALFree(list);
   /*
    * finally close files that may have been open
    */

   fclose(logfile);


   failAnalysis:

   /*
    * tell all slaves no dataset will be sent, 
    */

   last = 1;
   for (slaveNum=1; slaveNum<nTasks; slaveNum++)
   { 
	   LALMPISendINT4 (status->statusPtr, &last, slaveNum, MPI_COMM_WORLD);
	   CHECKSTATUSPTR(status);
	   fprintf(stderr, "Asked slave %d to shut\n",slaveNum);
	   fprintf(logfile, "Asked slave %d to shut\n",slaveNum);

   }
   DETATCHSTATUSPTR(status);
   RETURN(status);
}
#endif
