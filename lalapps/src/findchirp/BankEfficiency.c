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
The options are :
USAGE: lalapps_BankEfficiency [options]
The options are (with default values in brackets)\n");
      -quiet : if this flag is present, the output is restricted to the minimum
          -n : number of trials
         -fl : lower frequency cutoff
       -mMin : minimal mass of component stars
       -mMax : maximal mass of component stars
      -x0Max : Max value of psi0 (%7.2f)
      -x1Min : Min value of psi3 (%7.2f)
      -alpha : BCV amplitude correction parameter
    -numFcut : number of layers in Fcut dimension
         -mm : minimal match for template bank
    -simType : type of simulation, 0, 1 or 2
-approximant : PN model for Monte-Carlo (TaylorT1, ...)
     -signal : same as -approximant
   -template : PN model for template bank (TaylorT1, ...)
      -order : order of PN model
       -seed : seed for random generation
     -sigAmp : amplitude of the signal

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
 
void printf_timeseries (INT4 n, REAL4 *signal, REAL8 delta, REAL8 t0);

INT4 lalDebugLevel=0;
     
int 
main (  int argc, char **argv ) 
{
   static LALStatus status;
	   
   INT4 ntrials, i, j, nlist, jmax, approx=0, template=8;
   REAL8 df, omax;
/*
   REAL8 dt0, dt1, g00, g11, h00, h11, h01, ang, match;
*/
   REAL4Vector signal, correlation;
   UINT4 quietFlag = 0;
   void *noisemodel = LALLIGOIPsd;
   static RandomInspiralSignalIn randIn;
   static InspiralWaveOverlapIn overlapin;
   static InspiralWaveOverlapOut overlapout, overlapoutmax;
   RealFFTPlan *fwdp=NULL,*revp=NULL;
   InspiralTemplateList *list=NULL;
   static InspiralCoarseBankIn  coarseIn; 
   
   quietFlag = 0;	
   ntrials=2;
   approx = 0;
  
   coarseIn.fLower = 40.L;
   coarseIn.fUpper = 1000.L;
   coarseIn.tSampling = 2048.L;
   coarseIn.order = twoPN;
   coarseIn.space = Tau0Tau3;
   coarseIn.mmCoarse = 0.95;
   coarseIn.mmFine = 0.97;
   coarseIn.iflso = 0.0L;
   coarseIn.mMin = 3.0;
   coarseIn.mMax = 10.0;
   coarseIn.MMax = coarseIn.mMax * 2.;
   coarseIn.massRange = MinMaxComponentMass; 
   /* coarseIn.massRange = MinComponentMassMaxTotalMass;*/
   /* minimum value of eta */
   coarseIn.etamin = coarseIn.mMin * ( coarseIn.MMax - coarseIn.mMin) / pow(coarseIn.MMax,2.);
   coarseIn.psi0Min = 1.e0;
   coarseIn.psi0Max = 3e5;
   coarseIn.psi3Min = -3.5e3;
   coarseIn.psi3Max = 0;
   coarseIn.alpha = 0.L;
   coarseIn.numFcutTemplates = 5;
   coarseIn.approximant = template;

   randIn.useed = 128092;
   randIn.type = 0;
   randIn.SignalAmp = 10.0;
   randIn.param.approximant = PadeT1;
   randIn.param.order = twoPN;
   /*
   randIn.param.fendBCV = 300;
   */
   randIn.param.ieta=1.0; 
   randIn.param.startTime=0.0; 
   randIn.param.startPhase=0.88189; 
   randIn.param.nStartPad=1000;
   randIn.param.signalAmplitude = 1.0;
   randIn.param.nEndPad = 0;

   i=1;
   while(i <argc)
   {
	   if (strcmp(argv[i],"-fl")==0)
		   coarseIn.fLower = atof(argv[++i]);
	   else if (strcmp(argv[i],"-mMin")==0)
		   coarseIn.mMin = atof(argv[++i]); 
	   else if (strcmp(argv[i],"-mMax")==0)
		   coarseIn.mMax = atoi(argv[++i]);
	   else if (strcmp(argv[i],"-x0Max")==0)
		   coarseIn.psi0Max =atof(argv[++i]);
	   else if (strcmp(argv[i],"-x1Min")==0)
		   coarseIn.psi3Min =atof(argv[++i]);
	   else if (strcmp(argv[i],"-mm")==0)
		   coarseIn.mmCoarse =atof(argv[++i]);
	   else if (strcmp(argv[i],"-numFcut")==0)
		   coarseIn.numFcutTemplates = atof(argv[++i]); 
	   else if (strcmp(argv[i],"-simType")==0)
		   randIn.type = atoi(argv[++i]);
	   else if (strcmp(argv[i],"-sigAmp")==0)
		   randIn.SignalAmp = atof(argv[++i]);
	   else if (strcmp(argv[i],"-alpha")==0)
		   coarseIn.alpha = atof(argv[++i]); 
	   else if (strcmp(argv[i],"-order")==0)
		   randIn.param.order = atoi(argv[++i]); 
	   else if (strcmp(argv[i],"-n")==0)
		   ntrials = atoi(argv[++i]);       
	   else if (strcmp(argv[i],"-seed")==0)
		   randIn.useed = atoi(argv[++i])+1;
	   else if (strcmp(argv[i],"-quiet")==0)
		   quietFlag = 1;
	   else if (strcmp(argv[i], "-template")==0)
	   {
		   if (strcmp(argv[++i],"TaylorT1")==0)
			   template = 0;
		   else if (strcmp(argv[i],"TaylorT2")==0)
			   template = 1;
		   else if (strcmp(argv[i],"TaylorT3")==0)
			   template = 2;
		   else if (strcmp(argv[i],"TaylorF1")==0)
			   template = 3;
		   else if (strcmp(argv[i],"TaylorF2")==0)
			   template = 4;
		   else if (strcmp(argv[i],"PadeT1")==0)
			   template = 5;
		   else if (strcmp(argv[i],"PadeF1")==0)
			   template = 6;
		   else if (strcmp(argv[i],"EOB")==0)
			   template = 7;
		   else if (strcmp(argv[i],"BCV")==0)
			   template = 8;
		   else if (strcmp(argv[i],"SpinTaylorT3")==0)
			   template = 9;
		   coarseIn.approximant = template;	
	   }	
	   else if (strcmp(argv[i], "-approximant")==0 || (strcmp(argv[i], "-signal")==0))
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
		   fprintf(stderr,"All options should be followed by a number except -quiet\n");
		   fprintf(stderr,"      -quiet : if this flag is present, the output is restricted to the minimum\n");
		   fprintf(stderr,"          -n : number of trials                   (%d)\n",       ntrials);
		   fprintf(stderr,"       -seed : seed for random generation         (%d)\n",       randIn.useed);
		   fprintf(stderr,"    -simType : type of simulation, 0, 1 or 2      (%d)\n\n",     randIn.type);
		   fprintf(stderr,"         -fl : lower frequency cutoff             (%7.2f) Hz\n", coarseIn.fLower);
		   fprintf(stderr,"       -mMin : minimal mass of component stars    (%7.2f) Mo\n", coarseIn.mMin);
		   fprintf(stderr,"       -mMax : maximal mass of component stars    (%7.2f) Mo\n", coarseIn.mMax);
		   fprintf(stderr,"      -x0Max : Max value of psi0 (%7.2f)\n",                     coarseIn.psi0Max);
		   fprintf(stderr,"      -x1Min : Min value of psi3 (%7.2f)\n",                     coarseIn.psi3Min);
		   fprintf(stderr,"      -alpha : BCV amplitude correction parameter (%7.2f)\n",    coarseIn.alpha);
		   fprintf(stderr,"    -numFcut : number of layers in Fcut dimension (%7.2d)\n",    coarseIn.numFcutTemplates);
		   fprintf(stderr,"         -mm : minimal match for template bank    (%7.3f)\n",    coarseIn.mmCoarse);
		   fprintf(stderr,"   -template : PN model for template bank, e.g. BCV (%d)\n\n",   coarseIn.approximant);
		   fprintf(stderr,"-approximant : PN model for Monte-Carlo, e.g. EOB (%d)\n",       randIn.param.approximant);
		   fprintf(stderr,"     -signal : same as -approximant\n");
		   fprintf(stderr,"      -order : order of PN model                  (%7.2d)\n",    randIn.param.order);
		   fprintf(stderr,"     -sigAmp : amplitude of the signal            (%7.2f)\n",    randIn.SignalAmp);
		   return 1;	

	   }
	   i++;       
   }


/*---------------------------------------------------------------------------*/
   overlapin.nBegin = 0;
   overlapin.nEnd = 0;
/*---------------------------------------------------------------------------*/
/* User can choose allowed values of the various parameters below this line  */
/*---------------------------------------------------------------------------*/
   randIn.NoiseAmp = 1.0;

   randIn.mMin = coarseIn.mMin;
   randIn.mMax = coarseIn.mMax;
   randIn.MMax = 2.*randIn.mMax;
   randIn.param.fLower = coarseIn.fLower;
   randIn.param.alpha = coarseIn.alpha;
   randIn.psi0Min = coarseIn.psi0Min;
   randIn.psi0Max = coarseIn.psi0Max;
   randIn.psi3Min = coarseIn.psi3Min;
   randIn.psi3Max = coarseIn.psi3Max;
   randIn.param.fCutoff = coarseIn.fUpper;
   randIn.param.tSampling = coarseIn.tSampling;
   randIn.param.massChoice = m1Andm2;
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
   memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
   coarseIn.shf.f0 = 0;
   LALDCreateVector( &status, &(coarseIn.shf.data), randIn.psd.length );
   coarseIn.shf.deltaF = randIn.param.tSampling / signal.length;
   randIn.psd.data = (REAL8*) LALMalloc(sizeof(REAL8)*randIn.psd.length);
   df = randIn.param.tSampling/(float) signal.length;
   LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, df);
   LALInspiralCreateCoarseBank(&status,&list, &nlist, coarseIn);
   if (nlist==0) exit(0);
	   
   if (!quietFlag)
   {
   for (i=0; i<nlist; i++)
   {
	   if (coarseIn.approximant == BCV)
	   {
	   
		   fprintf(stdout, "%e %e %e %e\n", list[i].params.psi0, list[i].params.psi3, list[i].params.totalMass, list[i].params.fendBCV);
	   }
	   else
	   {
		   fprintf(stdout, "%e %e %e %e\n", list[i].params.mass1, list[i].params.mass2, list[i].params.t0, list[i].params.t3);
	   }
   }
   }
/* REPORTSTATUS(&status); */

   randIn.psd = *coarseIn.shf.data;
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
   
    if (!quietFlag) 
   {
   	fprintf(stdout, "#---------------------------------------------------------------\n");
	fprintf(stdout, "#Number of Coarse Bank Templates=%d\n",nlist);
   	fprintf(stdout, "#Signal Length=%d, Number of sims=%d\n", signal.length, ntrials);
	fprintf(stdout,"---------------------------------------------------------------\n");
   	if (coarseIn.approximant==BCV)
	{
		fprintf(stdout, "#psi0Min=%e, psi0Max=%e, psi3Min=%e, psi3Max=%e, signal=%d, template=%d\n", coarseIn.psi0Min,coarseIn.psi0Max,coarseIn.psi3Min,coarseIn.psi3Max,approx,template);
		fprintf(stdout, "#psi0Tmplt   psi0Sig    psi3Tmplt    psi3Sig    fFinalTmplt   fFinalSig  m1 m2 Overlap/SNR\n");
	}
	else
	{
		fprintf(stdout, "#mMin=%e, mMax=%e, signal=%d, template=%d\n", coarseIn.mMin,coarseIn.mMax,approx,template);
		fprintf(stdout, "#m1Tmplt   m1Sig    m2Tmplt    m2Sig    fFinalTmplt   fFinalSig  Overlap/SNR\n");
	}
   	fprintf(stdout,"---------------------------------------------------------------\n");
   }
   
   while (ntrials--) 
   {
      randIn.param.approximant = approx;
      randIn.param.fCutoff = coarseIn.fUpper;;
      if (approx==BCV) 
	      randIn.param.massChoice = psi0Andpsi3;
      else
	      randIn.param.massChoice = m1Andm2;

      for (i=0; i<signal.length; i++) signal.data[i] = 0.;
      LALRandomInspiralSignal(&status, &signal, &randIn);

      overlapin.signal = signal;
      jmax = 0;
      omax = 0.0;
	      
      for (j=0; j<nlist; j++) 
      {
     	      overlapin.param = list[j].params;
	      if (overlapin.param.approximant==BCV) 
		      overlapin.param.fCutoff = list[j].params.fendBCV;
	      else
		      overlapin.param.fCutoff = randIn.param.fCutoff;

	      for (i=0; i<signal.length; i++) correlation.data[i] = 0.;
	      LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
              list[j].params.fFinal = overlapin.param.fFinal;
	      
	      if (omax < overlapout.max) 
	      {
		      omax = overlapout.max;
		      overlapoutmax = overlapout;
		      jmax = j;
	      }
      }

      LALInspiralParameterCalc(&status, &list[jmax].params);
      if (coarseIn.approximant==BCV)
	      fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e\n", 
		      list[jmax].params.psi0, 
		      randIn.param.psi0, 
                      list[jmax].params.psi3,
		      randIn.param.psi3, 
		      list[jmax].params.fFinal, 
		      randIn.param.fFinal,
		      list[jmax].params.totalMass,
		      randIn.param.totalMass,
		      randIn.param.mass1,
		      randIn.param.mass2,
		      omax);
      else
	      fprintf(stdout, "%16.12f %16.12f %16.12f %16.12f %e %e %e %e %e %e %e\n", 
		      list[jmax].params.mass1, 
		      randIn.param.mass1, 
                      list[jmax].params.mass2,
		      randIn.param.mass2, 
                      list[jmax].params.totalMass,
		      randIn.param.totalMass, 
                      list[jmax].params.eta,
		      randIn.param.eta, 
		      list[jmax].params.fFinal, 
		      randIn.param.fFinal,
		      omax);
      fflush(stdout);
   }
   
   /* destroy the plans, correlation and signal */

   
      LALDDestroyVector( &status, &(coarseIn.shf.data) );
      if (signal.data != NULL) LALFree(signal.data);
      if (correlation.data != NULL) LALFree(correlation.data);
      if (randIn.psd.data != NULL) LALFree(randIn.psd.data);
      LALDestroyRealFFTPlan(&status,&fwdp);   
      LALDestroyRealFFTPlan(&status,&revp);
      LALCheckMemoryLeaks();
    
   
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
