/******************************** <lalVerbatim file="GEO600InspiralSearch">
Author: Sathyaprakash, B. S.
< $Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{GEO600InspiralSearch.c}}
\label{ss:GEO600InspiralSearch.c}

LAL code for inspiral searches.

\subsubsection*{Usage}
\begin{verbatim}
GEO600InspiralSearch
\end{verbatim}

\subsubsection*{Description}

This code filters the GEO600 data through a template bank.
\subsubsection*{Exit codes}
\input{GEO600InspiralSearchCE}

\subsubsection*{Uses}
This code directly uses the following functions (see those functions
to find out what they call in turn):
\begin{verbatim}
LALInspiralWaveLength
LALInspiralCreateCoarseBank
LALNoiseSpectralDensity 
RunningAverageDoublePSD
LALCreateForwardRealFFTPlan
LALCreateReverseRealFFTPlan
LALREAL4VectorFFT
LALDestroyRealFFTPlan
LALDestroyVector
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GEO600InspiralSearchCV}}
******************************************************* </lalLaTeX> */

/***************************** <lalErrTable file="GEO600InspiralSearchCE"> */
/***************************** </lalErrTable> */

#include <config.h>
#include <stdio.h>

#if !defined HAVE_GSL_GSL_FFT_REAL_H || !defined HAVE_LIBGSL \
 || !defined HAVE_MYSQL_H || !defined HAVE_LIBMYSQLCLIENT \
 || !defined HAVE_FRAMEL_H || !defined HAVE_LIBFRAME

int main( void ) { return fputs( "disabled\n", stderr ), 77; }

#else

#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>
#include <lal/GEO600InspiralSearch.h>
#include <lal/FrameRead.h>
#include <lal/AVFactories.h>
#include <gsl/gsl_fft_real.h>
 
NRCSID (INSPIRALSEARCHC, "$Id$");

INT4 lalDebugLevel=1;
     
int 
main () 
{
   InspiralSearchParams searchParams; /* parameters of the search */
   InspiralCoarseBankIn coarseIn; /* parameter stucture to find coarse bank templates */
   InspiralTemplateList *list=NULL; /* list of templates */
   InspiralFindEventsIn findeventsin; /* parameter structure for finding events */
   UINT4 templateLength, dataSetLength, psdLength; /* Lengths of template, data and S_h(f) */
   UINT4 dataSetDuration, templateDuration; /* Duration of data and template  */
   UINT4 overlapDuration; /* Duration for which consecutive data sets must overlap */
   INT4 currentGPSTime; /* GPS time at the beginning of the data set currently being analysed */
   INT4 nlist; /* Number of templates in coarse bank */
   INT4 nEvents; /* Number of events found in the current data set */
   REAL8Vector *psd=NULL; /* Vector to hold S_h(f) */
   REAL4Vector *dataset=NULL; /* Vectors to hold data set */
   REAL8Vector *buff=NULL; /* Vectors to hold data set */
   REAL4Vector *hannWindow=NULL; /* Vector to hold window used in computing psd */
   RealFFTPlan *fwdp=NULL,*revp=NULL;
   FILE *GEO600InspiralSearch;
   char* channelName;
   LALWindowParams winParams;
   RunningAveragePSDIn psdIn;
   INT4 i, nby2, j, ret, last;
   static float sampleRate;
   REAL8 dt, df;
   static LALStatus status; 

   if ((GEO600InspiralSearch = fopen("GEO600InspiralSearch.out", "w"))==NULL)
     {
      fprintf(stderr, "Cannot open file to write output\n");
      exit(0);
     }
/*------------------------------------------------------------------*/
/*-------- Read the parameters of search ---------------------------*/
/*------------------------------------------------------------------*/
   ReadInspiralSearchParams(&status, &searchParams);
   sampleRate = (float) searchParams.samplingRate;
/*------------------------------------------------------------------*/
/*-------- Read the parameters of the template bank ----------------*/
/*------------------------------------------------------------------*/
   ReadInspiralTemplateBankParams(&status, &coarseIn);
   coarseIn.tSampling = sampleRate;
   dt = 1./sampleRate;
/*------------------------------------------------------------------*/
/*-------- Get the length of the longest template ------------------*/
/*------------------------------------------------------------------*/
   LALInspiralLongestTemplateInBank(&status, &templateLength, &coarseIn);
/*------------------------------------------------------------------*/
/*-------- Set array lengths ---------------------------------------*/
/*------------------------------------------------------------------*/
   if (templateLength < (UINT4) sampleRate) 
     {
       templateLength = (UINT4) sampleRate;
     }
   templateDuration = templateLength / sampleRate;
   dataSetDuration = searchParams.paddingFactor * templateDuration; 
   dataSetLength = dataSetDuration * sampleRate;
   overlapDuration = 2 * dataSetDuration/searchParams.paddingFactor;
   nby2 = dataSetLength/2;
   psdLength = dataSetLength/2+1;
   searchParams.dataSetLength = dataSetLength;
   searchParams.dataSetDuration = dataSetDuration;
/*------------------------------------------------------------------*/
/*-------- Write some useful information on screen -----------------*/
/*------------------------------------------------------------------*/
   fprintf(stderr, "SamplingRate=%e,\n",sampleRate);
   fprintf(stderr, "TemplateLength=%d, dataSetLength=%d\n",
                    templateLength,dataSetLength); 
   fprintf(stderr, "TemplateDuration=%d, dataSetDuration=%d\n",
                    templateDuration,dataSetDuration); 
/*------------------------------------------------------------------*/
/*--------- Estimate FFTW plans ------------------------------------*/
/*------------------------------------------------------------------*/
   LALCreateForwardRealFFTPlan(&status, &fwdp, dataSetLength, 0);
   LALCreateReverseRealFFTPlan(&status, &revp, dataSetLength, 0);
/*------------------------------------------------------------------*/
/*----------Create dataSet, buffer, window and psd vectors ---------*/
/*------------------------------------------------------------------*/
   LALCreateVector(&status, &dataset, dataSetLength);
   LALCreateVector(&status, &hannWindow, dataSetLength);
   LALDCreateVector(&status, &buff, dataSetLength);
   LALDCreateVector(&status, &psd, psdLength);
/*------------------------------------------------------------------*/
/*-------- Create the window to be used for PSD computation --------*/
/*------------------------------------------------------------------*/
   winParams.length = dataSetLength;
   winParams.type = Hann;
   LALWindow (&status, hannWindow, &winParams);
   psdIn.norm = winParams.sumofsquares; 
   psdIn.fwdp = fwdp; 
   psdIn.revp = revp; 
   psdIn.window = hannWindow; 
   psdIn.searchParams = searchParams; 
   psdIn.norm = winParams.sumofsquares; 
/*
   for (j=0; j<dataSetLength; j++) printf("%e %e\n", j*dt, hannWindow->data[j]);
   agraphout(0.0, (float) (j*dt), last);
   exit(0);
*/
   memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries));
   LALDCreateVector (&status, &(coarseIn.shf.data), psdLength);
   coarseIn.shf.f0 = 0.;
   coarseIn.shf.deltaF = coarseIn.tSampling / (REAL8) coarseIn.shf.data->length;
/*------------------------------------------------------------------*/
/*-------- Compute the power spec-----------------------------------*/
/*------------------------------------------------------------------*/
   RunningAverageDoublePSD(&status, coarseIn.shf.data, &ret, &psdIn);
/*------------------------------------------------------------------*/
/*-------- Get the template bank -----------------------------------*/
/*------------------------------------------------------------------*/
   LALInspiralCreateCoarseBank(&status, &list, &nlist, coarseIn);
   /*
   fprintf(stderr, "Number of Templates=%d\n", nlist);
   scanf("%d", &nlist);
   printf("Reading %d templates\n", nlist);
   list = (InspiralTemplateList*) malloc(sizeof(InspiralTemplateList)*(nlist));
   for (i=0; i<nlist; i++) 
   {
	   REAL8 t1, t2, t3;
	   scanf("%d %le %le %le %le %le", &j, &(list[i].params.mass1), &(list[i].params.mass2), &t1, &t2, &t3);
	   list[i].params.massChoice = m1Andm2;
	   list[i].params.ieta = 1;
	   list[i].params.signalAmplitude = 1.;
	   list[i].params.tSampling = coarseIn.tSampling;
	   list[i].params.fLower = coarseIn.fLower;
	   list[i].params.fCutoff = coarseIn.fUpper;
	   list[i].params.order = coarseIn.order;
	   list[i].params.approximant = coarseIn.approximant;
	   list[i].params.nStartPad = 0;
	   list[i].params.nEndPad = 0;
	   list[i].params.Theta = 0.;
	   list[i].params.OmegaS = 0.;
	   list[i].params.startTime = 0.;
	   list[i].params.startPhase = 0.;
	   list[i].params.spin1[0] = list[i].params.spin1[1] = list[i].params.spin1[2] = 0.;
	   list[i].params.spin2[0] = list[i].params.spin2[1] = list[i].params.spin2[2] = 0.;
	   list[i].params.inclination = 0.;
	   list[i].params.eccentricity = 0.;
	   list[i].params.rInitial = 0.;
	   list[i].params.vInitial = 0.;
	   list[i].params.rFinal = 0.;
	   list[i].params.vFinal = 0.;
	   list[i].params.rLightRing = 0.;
	   LALInspiralParameterCalc(&status, &(list[i].params));
   }
   */
   for (i=0; i<nlist; i++) 
	   fprintf(stderr, "%e %e %e\n", list[i].params.mass1, list[i].params.mass2, list[i].params.totalMass);

   /*
   for (i=1; i<psdLength; i++) 
   {
	   scanf("%d %le\n", &j, &psd->data[i]);
	   coarseIn.shf.data->data[i] = psd->data[i];
   }
   LALInspiralCreateCoarseBank(&status, &list, &nlist, coarseIn);
   for (i=0; i<nlist; i++) 
	   printf("%e %e %e\n", list[i].params.mass1, list[i].params.mass2, list[i].params.totalMass);
   */

   if (nlist==0) 
   {
	   fprintf(stderr, "No templates to analyse\n");
	   exit(0);
   }
/*------------------------------------------------------------------*/
/*----------Cast the structure needed by LALInspiralFindEvents -----*/
/*------------------------------------------------------------------*/
   findeventsin.fwdp = fwdp;
   findeventsin.revp = revp;
   findeventsin.signal = *dataset;
   findeventsin.psd = *psd;
   findeventsin.nBegin = dataSetLength / searchParams.paddingFactor;
   findeventsin.nEnd = dataSetLength / searchParams.paddingFactor;
   findeventsin.Threshold = 8.;
   findeventsin.ClusterThreshold = 3.;
   findeventsin.displayData = 0;
   findeventsin.displayPSD = 0;
   findeventsin.displayTemplates = 0;
   findeventsin.displayCorrelation = 0;
   findeventsin.displayCorrelationStats = 1;
/*------------------------------------------------------------------*/
   df = coarseIn.tSampling/(REAL8)dataSetLength;
   psdIn.fwdp = fwdp;
   psdIn.window = hannWindow;
/*
   RunningAverageDoublePSD(&status, psd, &psdIn);
   for (i=1; i<nby2; i++) 
      printf("%e %e\n", i*df, psd->data[i]);
   exit(0);
*/
   last = 0;
   channelName = searchParams.channelName;
   while (searchParams.currentGPSTime <= searchParams.endGPSTime)
   {
      currentGPSTime = searchParams.currentGPSTime;
      fprintf(stderr, "currentGPSTime=%d\n",currentGPSTime); 
/*------------------------------------------------------------------*/
/*----------Compute noise PSD --------------------------------------*/
/*------------------------------------------------------------------*/
      psdIn.searchParams = searchParams;
      RunningAverageDoublePSD(&status, psd, &ret, &psdIn);
      findeventsin.psd = *psd;
/*
      LALNoiseSpectralDensity(&status, psd, coarseIn.NoisePsd, df);
*/
      if (findeventsin.displayPSD)
      {
         for (j=1; j<(INT4)psd->length-1; j++) 
            printf("%e %e\n", j*df, psd->data[j]);
         printf("&\n");
	 exit(0);
      }
/*------------------------------------------------------------------*/
/*----------Get next data set and increment time -------------------*/
/*------------------------------------------------------------------*/
	      
      ret=GetChannelDoubleData(currentGPSTime, dataSetDuration, channelName, buff->data, dataSetLength);
      if (ret) exit(0);
/*
      list[0].params.nStartPad = dataSetLength/2;
      LALInspiralWave (&status, buff, &(list[0].params));
      list[0].params.nStartPad = 0;
*/
      if (findeventsin.displayData)
      {
         for (j=0; j<(INT4)dataSetLength; j++) printf("%e %e\n", j*dt, buff->data[j]);
      }
/*------------------------------------------------------------------*/
/*----------Compute FFT of data and convert to sngl precession------*/
/*------------------------------------------------------------------*/
      gsl_fft_real_radix2_transform(buff->data, 1, buff->length);
      for (j=0; j<(INT4)dataSetLength; j++) 
      {
	      dataset->data[j] = (REAL4) buff->data[j];
	      /*
	      printf("%e %e\n", j*dt, dataset->data[j]);
	      */
      }
      findeventsin.currentGPSTime = currentGPSTime;
      findeventsin.signal = *dataset;
/*------------------------------------------------------------------*/
/*----------Filter data through template bank ----------------------*/
/*------------------------------------------------------------------*/
      for (j=0; j<nlist; j++) 
      { 
          InspiralEventsList *eventList=NULL;
          findeventsin.param = list[j].params;
	  /*
	  LALInspiralFindLoudestEvent(&status, &nEvents, &eventList, &findeventsin);
	  */
	  LALInspiralFindEventsCluster(&status, &nEvents, &eventList, &findeventsin);
          fprintf(stderr, "GPS=%d, Template=%d, numEvents=%d\n", currentGPSTime, j, nEvents); 
          if (nEvents) 
	  {
		  WriteInspiralEvents(&status, eventList, nEvents, GEO600InspiralSearch);
		  LALFree(eventList);
		  eventList=NULL;
	  }
      }
      searchParams.currentGPSTime += (dataSetDuration-overlapDuration);
   }
   fclose(GEO600InspiralSearch);
   LALDestroyVector(&status, &hannWindow);
   LALDestroyVector(&status, &dataset);
   LALDDestroyVector(&status, &psd);
   LALDDestroyVector(&status, &buff);
   LALDestroyRealFFTPlan(&status,&fwdp);   
   LALDestroyRealFFTPlan(&status,&revp);
   LALDDestroyVector(&status, &(coarseIn.shf.data));
   if (list != NULL) LALFree(list);
   LALCheckMemoryLeaks();
   return 0;
}
#endif
