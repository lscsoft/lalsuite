/******************************** <lalVerbatim file="BankEfficiencyCV">
Author: Cokelaer, T. and Sathyaprakash, B. S.
< $Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{BankEfficiency.c}}
\label{ss:BankEfficiency.c}

Test code for the inspiral bank modules.

\subsubsection*{Usage}
\begin{verbatim}
BankEfficiency [options]

-n : number of trials
-seed : seed for random generation 
-mm : minimal match 
-fl : lower frequency cutoff
-mMin : minimal mass of component stars 
-mMax : maximal mass of component stars 
-alpha : BCV amplitude correction parameter 
-approximant : Post-Newtonian model such as TaylorT1, PadeT1, EOB, BCV ...
\end{verbatim}

\subsubsection*{Description}

This test code gives an example of how one might generate inspiral
waveforms and use them to compute the overlap of a random signal
(with or witnout simulated noise) of a certain strength. The parameter
{\texttt randIn.type=0} generates only signal
{\texttt randIn.type=1} generates only noise
{\texttt randIn.type=2} generates {\texttt randIn.SignalAmp * signal(t)
+ randIn.NoiseAmp * noise(t)}
Note that one must calculate the length of the waveform and allocate memory for it
\emph{before} calling
\texttt{InspiralWave}. The length of the waveform can be calculated by calling the function
\texttt{InspiralWaveLength} beforehand, as shown.

There are only two functions which one can call to generate waveforms. These are
\texttt{InspiralWave},
which will return a \emph{single} waveform, and \texttt{InspiralWaveTemplates}, which
returns a \emph{pair}
of waveforms which have phases which differ by $\pi/2$.

\subsubsection*{Exit codes}
\input{BankEfficiencyCE}

\subsubsection*{Uses}
This code directly uses the following functions (see those functions
to find out what they call in turn):
\begin{verbatim}
LALInspiralWaveLength
LALInspiralCreateCoarseBank
LALRandomInspiralSignal
LALInspiralParameterCalc
LALNoiseSpectralDensity 
LALCreateForwardRealFFTPlan
LALCreateReverseRealFFTPlan
LALForwardRealFFT
LALReverseRealFFT
LALDestroyRealFFTPlan
LALInspiralWaveOverlap
LALInspiralParameterCalc
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{BankEfficiencyCV}}
******************************************************* </lalLaTeX> */

/***************************** <lalErrTable file="BankEfficiencyCE"> */
/***************************** </lalErrTable> */

#include <stdio.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
 
void  LALInspiralCreateFlatBank(LALStatus *status, REAL4VectorSequence *list, InspiralBankParams *bankParams);
void printf_timeseries (INT4 n, REAL4 *signal, REAL8 delta, REAL8 t0);

INT4 lalDebugLevel=0;
     
int 
main (  int argc, char **argv ) 
{
   REAL4VectorSequence *list=NULL;
   static LALStatus status;
   static InspiralBankParams bankParams;
	   
   INT4 ntrials, i, j, nlist, jmax, approx;
   REAL8 df, omax;
/*
   REAL8 dt0, dt1, g00, g11, h00, h11, h01, ang, match;
*/
   REAL4Vector signal, correlation;
   UINT4 numFcutTemplates = 5, quietFlag = 0;
   void *noisemodel = LALLIGOIPsd;
   static RandomInspiralSignalIn randIn;
   static InspiralWaveOverlapIn overlapin;
   static InspiralWaveOverlapOut overlapout, overlapoutmax;
   static CreateVectorSequenceIn in; 
   REAL8FrequencySeries shf;
   static InspiralMetric metric;
   RealFFTPlan *fwdp=NULL,*revp=NULL;
   static InspiralTemplateList *tmpltList=NULL;
   
   
   quietFlag = 0;	
   numFcutTemplates  = 5;
   randIn.useed = 128092;
   ntrials=2;
   randIn.mMin = 5.0;
   randIn.mMax = 20.0;
   randIn.param.fLower = 60;
   randIn.param.alpha = 0;  
   randIn.param.fendBCV = 300;
   randIn.type = 0;
   randIn.SignalAmp = 10.0;
   randIn.param.approximant = PadeT1;
   bankParams.minimalMatch = 0.80;
   approx = 0;
   randIn.param.order = 4;
   
   bankParams.x0Min = 0.0;
   bankParams.x0Max = 2.5e5;
   bankParams.x1Min = -2.2e3;
   bankParams.x1Max = 0.0;


   i=1;
   while(i <argc)
   {
	   if (strcmp(argv[i],"-fl")==0)
		   randIn.param.fLower = atof(argv[++i]);
	   else if (strcmp(argv[i],"-mMin")==0)
		   randIn.mMin = atof(argv[++i]); 
	   else if (strcmp(argv[i],"-numFcut")==0)
		   numFcutTemplates = atof(argv[++i]); 
	   else if (strcmp(argv[i],"-simType")==0)
		   randIn.type = atoi(argv[++i]);
	   else if (strcmp(argv[i],"-quiet")==0)
		   quietFlag = 1;
	   else if (strcmp(argv[i],"-mMax")==0)
		   randIn.mMax = atoi(argv[++i]);
	   else if (strcmp(argv[i],"-sigAmp")==0)
		   randIn.SignalAmp = atof(argv[++i]);
	   else if (strcmp(argv[i],"-alpha")==0)
		   randIn.param.alpha = atof(argv[++i]); 
	   else if (strcmp(argv[i],"-order")==0)
		   randIn.param.order = atoi(argv[++i]); 
	   else if (strcmp(argv[i],"-n")==0)
		   ntrials = atoi(argv[++i]);       
	   else if (strcmp(argv[i],"-seed")==0)
		   randIn.useed = atoi(argv[++i])+1;
	   else if (strcmp(argv[i],"-x0Max")==0)
		   bankParams.x0Max =atof(argv[++i]);
	   else if (strcmp(argv[i],"-x1Min")==0)
		   bankParams.x1Min =atof(argv[++i]);
	   else if (strcmp(argv[i],"-mm")==0)
		   bankParams.minimalMatch =atof(argv[++i]);
	   else if (strcmp(argv[i], "-approximant")==0)
	   {
		   if (strcmp(argv[++i],"TaylorT1")==0)
			   approx = 0;
		   else if (strcmp(argv[i],"TaylorT2")==0)
			   approx = 1;
		   else if (strcmp(argv[i],"TaylorT3")==0)
			   approx = 2;
		   else if (strcmp(argv[i],"TaylorF1")==0)
			   approx = 3;
		   else if (strcmp(argv[i],"TaylorF2")==0)
			   approx = 4;
		   else if (strcmp(argv[i],"PadeT1")==0)
			   approx = 5;
		   else if (strcmp(argv[i],"PadeF1")==0)
			   approx = 6;
		   else if (strcmp(argv[i],"EOB")==0)
			   approx = 7;
		   else if (strcmp(argv[i],"BCV")==0)
			   approx = 8;
		   else if (strcmp(argv[i],"SpinTaylorT3")==0)
			   approx = 9;
		   randIn.param.approximant = approx;	
	   }	
	   else 
	   {
		   fprintf(stderr,"\nUSAGE: %s [options]\n", argv[0]);
		   fprintf(stderr,"The options are (with default values in brackets)\n");
		   fprintf(stderr,"-approximant : Post-Newtonian model such as TaylorT1, PadeT1, EOB, BCV ...  (PadeT1)\n");
		   fprintf(stderr,"       -mMin : minimal mass of component stars    (%7.2f) SolarMass\n", randIn.mMin);
		   fprintf(stderr,"       -mMax : maximal mass of component stars    (%7.2f) SolarMass\n", randIn.mMax);
		   fprintf(stderr,"     -sigAmp : amplitude of the signal            (%7.2f)\n", randIn.SignalAmp);
		   fprintf(stderr,"      -alpha : BCV amplitude correction parameter (%7.2f)\n",randIn.param.alpha);
		   fprintf(stderr,"         -fl : lower frequency cutoff             (%7.2f) Hz\n", randIn.param.fLower);
		   fprintf(stderr,"      -order : order of PN model                  (%7.2d)\n", randIn.param.order);
		   fprintf(stderr,"    -numFcut : number of layers in Fcut dimension (%7.2d)\n", numFcutTemplates);
		   fprintf(stderr,"         -mm : minimal match for template bank    (%7.3f)\n", bankParams.minimalMatch);
		   fprintf(stderr,"       -seed : seed for random generation         (%d)\n", randIn.useed);
		   fprintf(stderr,"          -n : number of trials                   (%d)\n", ntrials);
		   fprintf(stderr,"    -simType : type of simulation, 0, 1 or 2      (%d)\n\n", randIn.type);
		   return 1;	

	   }
	   i++;       
   }

   /*
   fprintf(stderr, "This test code does three things:\n");
   fprintf(stderr, "(a) Creates a filter bank at a certain minimal match\n");
   fprintf(stderr, "(b) Generates signals with random parmeters\n");
   fprintf(stderr, "(c) Filters each of these signals with the a \n");
   fprintf(stderr, "subset of templates close to the random signal \n");
   fprintf(stderr, "and reports the best SNR achived\n");
   fprintf(stderr, "Results of the run are written on to stdout\n");
   */


/*---------------------------------------------------------------------------*/
   overlapin.nBegin = 0;
   overlapin.nEnd = 0;
/*---------------------------------------------------------------------------*/
/* User can choose allowed values of the various parameters below this line  */
/*---------------------------------------------------------------------------*/
   randIn.NoiseAmp = 1.0;

   randIn.param.startTime=0.0; 
   randIn.param.startPhase=0.88189; 
   randIn.param.nStartPad=2000;
   randIn.psi0Min = bankParams.x0Min;
   randIn.psi0Max = bankParams.x0Max;
   randIn.psi3Max = bankParams.x1Max;
   randIn.psi3Min = bankParams.x1Min;

   randIn.param.fCutoff = 1000.;
   randIn.param.tSampling = 2048;
   randIn.param.signalAmplitude = 1.0;
   randIn.param.nEndPad = 2000;
   randIn.param.massChoice = m1Andm2;
   randIn.MMax = 2.*randIn.mMax;
   randIn.param.mass1 = randIn.mMin;
   randIn.param.mass2 = randIn.mMin;

   signal.length = 0.;
   randIn.param.approximant = EOB;
   LALInspiralWaveLength (&status, &signal.length, randIn.param);
   randIn.param.approximant = approx;
   if (!quietFlag)
	fprintf(stdout, "signal length = %d\n", signal.length);
/* REPORTSTATUS(&status); */
   correlation.length = signal.length;
   randIn.psd.length = signal.length/2 + 1;
   signal.data = (REAL4*) LALMalloc(sizeof(REAL4)*signal.length);
   correlation.data = (REAL4*) LALMalloc(sizeof(REAL4)*correlation.length);
   memset( &(shf), 0, sizeof(REAL8FrequencySeries) );
   shf.f0 = 0;
   LALDCreateVector( &status, &(shf.data), randIn.psd.length );
   shf.deltaF = randIn.param.tSampling / signal.length;
   randIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8)*randIn.psd.length);
   df = randIn.param.tSampling/(float) signal.length;
   LALNoiseSpectralDensity (&status, shf.data, noisemodel, df);
	   
   /*---------------------------------------------------------------------------*/
   /* Prepare parameters needed to create a verctor sequence to store BCV templates */
   /*---------------------------------------------------------------------------*/
   in.length = 1;
   in.vectorLength = 2;
   LALSCreateVectorSequence(&status, &list, &in);
   list->vectorLength = 2;
   LALInspiralComputeMetricBCV(&status, &metric, &shf, &randIn.param);
   bankParams.metric = &metric;
   LALInspiralCreateFlatBank(&status, list, &bankParams);
   nlist = list->length;

   tmpltList = (InspiralTemplateList *) LALMalloc (sizeof (InspiralTemplateList) * nlist);

   /* Print out the template parameters */
   for (i=0; i<nlist; i++)
   {
   	/*Retain only those templates that have meaningful chirptimes*/
   	tmpltList[i].params.psi0 = (REAL8) list->data[2*i];
   	tmpltList[i].params.psi3 = (REAL8) list->data[2*i+1];
   	tmpltList[i].params.fLower = randIn.param.fLower;
   	tmpltList[i].params.nStartPad = randIn.param.nStartPad;
   	tmpltList[i].params.startPhase = randIn.param.startPhase;
	tmpltList[i].params.nEndPad = randIn.param.nEndPad;
	tmpltList[i].params.tSampling = randIn.param.tSampling;
	tmpltList[i].params.fendBCV = randIn.param.fendBCV;
	tmpltList[i].params.massChoice = psi0Andpsi3;
	tmpltList[i].params.approximant = BCV;
   }

	LALInspiralBCVFcutBank( &status, &tmpltList, &nlist, numFcutTemplates) ;

   if (nlist==0) exit(0);
   if (!quietFlag){
   for (i=0; i<nlist; i++)
   {
	   fprintf(stdout, "%e %e %e %e\n", 
			   tmpltList[i].params.psi0, 
			   tmpltList[i].params.psi3, 
			   tmpltList[i].params.totalMass, 
			   tmpltList[i].params.fendBCV);
   }
   }
/* REPORTSTATUS(&status); */

   randIn.psd = *shf.data;
   overlapin.psd = randIn.psd;
/*--------------------------
   Estimate the plans 
--------------------------*/
   LALCreateForwardRealFFTPlan(&status, &fwdp, signal.length, 0);
/* REPORTSTATUS(&status); */
   LALCreateReverseRealFFTPlan(&status, &revp, signal.length, 0);
/* REPORTSTATUS(&status); */
   overlapin.fwdp = randIn.fwdp = fwdp;
   overlapin.revp = revp;
   
   if (!quietFlag){
   	fprintf(stdout, "#---------------------------------------------------------------\n");
	fprintf(stdout, "#Number of Coarse Bank Templates=%d\n",nlist);
   	fprintf(stdout, "#Signal Length=%d, Number of sims=%d\n", signal.length, ntrials);
   	fprintf(stdout, "#psi0Min=%e, psi0Max=%e, psi3Min=%e, psi3Max=%e\n",
		   randIn.psi0Min,randIn.psi0Max,randIn.psi3Min,randIn.psi3Max);
	fprintf(stdout,"---------------------------------------------------------------\n");
   	fprintf(stdout, "#psi0Tmplt   psi0Sig    psi3Tmplt    psi3Sig    fFinalTmplt   fFinalSig  m1 m2 Overlap/SNR\n");
   	fprintf(stdout,"---------------------------------------------------------------\n");
	}
   
   while (ntrials--) 
   {
      randIn.param.approximant = approx;
      if (approx==BCV) 
	      randIn.param.massChoice = psi0Andpsi3;
      else
	      randIn.param.massChoice = m1Andm2;

      LALRandomInspiralSignal(&status, &signal, &randIn);

      overlapin.signal = signal;
      jmax = 0;
      omax = 0.0;
	      
      for (j=0; j<nlist; j++) 
      {
     	      overlapin.param = tmpltList[j].params;
	      if (overlapin.param.approximant==BCV) 
		      overlapin.param.fCutoff = tmpltList[j].params.fendBCV;
	      else
		      overlapin.param.fCutoff = randIn.param.fCutoff;

	      LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
              tmpltList[j].params.fFinal = overlapin.param.fFinal;
	      
	      if (omax < overlapout.max) 
	      {
		      omax = overlapout.max;
		      overlapoutmax = overlapout;
		      jmax = j;
	      }
      }

      LALInspiralParameterCalc(&status, &tmpltList[jmax].params);
      fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e\n", 
		      tmpltList[jmax].params.psi0, 
		      randIn.param.psi0, 
                      tmpltList[jmax].params.psi3,
		      randIn.param.psi3, 
		      tmpltList[jmax].params.fFinal, 
		      randIn.param.fFinal,
		      tmpltList[jmax].params.totalMass,
		      randIn.param.totalMass,
		      randIn.param.mass1,
		      randIn.param.mass2,
		      omax);
      fflush(stdout);
   }
   
   /* destroy the plans, correlation and signal */

   /*
      LALDDestroyVector( &status, &(shf.data) );
      if (signal.data != NULL) LALFree(signal.data);
      if (correlation.data != NULL) LALFree(correlation.data);
      if (randIn.psd.data != NULL) LALFree(randIn.psd.data);
      if (list!= NULL) LALFree(list);
      LALDestroyRealFFTPlan(&status,&fwdp);   
      LALDestroyRealFFTPlan(&status,&revp);
      LALCheckMemoryLeaks();
    */
   
   return(0);

}

void printf_timeseries (INT4 n, REAL4 *signal, double delta, double t0) 
{
  int i=0;

  do 
     printf ("%e %e\n", i*delta+t0, *(signal+i));
  while (n-++i); 
  printf("&\n");
}
