/******************************** <lalVerbatim file="BankEfficiencyCV">
Author: Cokelaer, T. and Sathyaprakash, B. S.
< $Id$
$Tag:$
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
      -psi0Max : Max value of psi0 (%7.2f)
      -psi1Min : Min value of psi3 (%7.2f)
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


#define BANKEFFICIENCY_FLOWER       	  40.
#define BANKEFFICIENCY_TSAMPLING    	2048.
#define BANKEFFICIENCY_FUPPER       	1000.
#define BANKEFFICIENCY_ORDER       	twoPN
#define BANKEFFICIENCY_MMCOARSE     	0.8
#define BANKEFFICIENCY_MMFINE       	0.97
#define BANKEFFICIENCY_MMIN            	5.
#define BANKEFFICIENCY_MMAX           	20.
#define BANKEFFICIENCY_ALPHA           	0.
#define BANKEFFICIENCY_FENDBCV       	300.
#define BANKEFFICIENCY_SPACE    	Tau0Tau3
#define BANKEFFICIENCY_IFLSO           	0.
#define BANKEFFICIENCY_PSI0MIN        	10.
#define BANKEFFICIENCY_PSI0MAX    	250000.
#define BANKEFFICIENCY_PSI3MIN     	-2200.
#define BANKEFFICIENCY_PSI3MAX       	-10.
#define BANKEFFICIENCY_NFCUT           	5
#define BANKEFFICIENCY_SIGNAL  		TaylorT1
#define BANKEFFICIENCY_BANK          	BCV
#define BANKEFFICIENCY_TYPE            	0
#define BANKEFFICIENCY_SIGNALAMP      	10.
#define BANKEFFICIENCY_IETA            	1.
#define BANKEFFICIENCY_STARTTIME       	0.
#define BANKEFFICIENCY_STARTPHASE    	0.9
#define BANKEFFICIENCY_NSTARTPHASE  	1000
#define BANKEFFICIENCY_SIGNALAMPLITUDE 	1.
#define BANKEFFICIENCY_NENDPAD         	0
#define BANKEFFICIENCY_NOISEAMP        	1.
#define BANKEFFICIENCY_NTRIALS         	2
#define BANKEFFICIENCY_QUIETFLAG       	0


typedef struct	
{
INT4 signal; /*random signal*/
INT4 template;/*template and bank*/
INT4 bank;
INT4 ntrials;
INT4 quietFlag;
} OtherParamIn;


void 	printf_timeseries (INT4 n, REAL4 *signal, REAL8 delta, REAL8 t0);
void 	Init2DummyCoarse(InspiralCoarseBankIn *coarseIn);
void 	Init2DummyRandom(RandomInspiralSignalIn *randIn);
void 	Init2DummyOtherParam(OtherParamIn *otherIn);
void 	ParseParameters(int 			*argc, 
			char 			**argv,
		     	InspiralCoarseBankIn 	*coarseIn,
		     	RandomInspiralSignalIn 	*randIn,
		     	OtherParamIn 		*otherIn);		     

void 	CheckParams(InspiralCoarseBankIn 	coarseIn,
		    RandomInspiralSignalIn 	randIn,
		    OtherParamIn 		otherIn);

void 	SetDefault(InspiralCoarseBankIn 	*coarseIn, 
		   RandomInspiralSignalIn 	*randIn,
		   OtherParamIn 		*otherIn);

void 	Help(InspiralCoarseBankIn 	coarseIn,
	     RandomInspiralSignalIn 	randIn,
	     OtherParamIn 		otherIn);

INT4 lalDebugLevel=1;
     

int 
main (  int argc, char **argv ) 
{
   /* --- Variables ---*/
   static LALStatus status;
   long	i, j;
   INT4 jmax, nlist;
   REAL8 df, omax;
   REAL4Vector  signal, correlation;
   void   			*noisemodel = LALLIGOIPsd;
   RandomInspiralSignalIn	 randIn;			/* random signal waveform to inject*/
   InspiralWaveOverlapIn 	 overlapin;			/* structure for Overlap*/
   InspiralWaveOverlapOut        overlapout, 
			         overlapoutmax;
   OtherParamIn 		 otherIn;			/* personal structure to parse some independant parameters*/
   RealFFTPlan  	        *fwdp=NULL,
		     	        *revp=NULL;
   InspiralTemplateList         *list=NULL;
   InspiralCoarseBankIn          coarseIn; 			/* strcture for the bank of templates*/


   
   /* --- Initialisation of parameters and variables --- */
   Init2DummyCoarse(&coarseIn);					/* the bank*/
   Init2DummyRandom(&randIn);					/* the waveform to inject*/
   Init2DummyOtherParam(&otherIn);				/* optional parameters to parse*/
   ParseParameters(&argc,argv, &coarseIn, &randIn, &otherIn);	/* Read Parameters from user */
   SetDefault(&coarseIn, &randIn, &otherIn);			/* Set default values to variables which have not been given by the user*/
   CheckParams(coarseIn, randIn, otherIn);			/* Check validity of some variables. */
   
   /* --- Some others variables to initialize --- */
   randIn.param.massChoice = m1Andm2;  				/* Only to compute the length of "signal"*/ 
   signal.length 	= 0.;
   randIn.param.approximant = EOB;				/* Only to compute the length of "signal"*/
   LALInspiralWaveLength (&status, &signal.length, randIn.param);
   randIn.param.approximant = otherIn.signal;			/* Retrieve the actual approximant of the waveform to inject*/
   if (!otherIn.quietFlag)
       fprintf(stdout, "signal length = %d\n", signal.length);  /* Optional output*/


   /* --- Allocate memory --- */
   correlation.length 	= signal.length;
   randIn.psd.length 	= signal.length/2 + 1;
   signal.data 		= (REAL4*) LALMalloc(sizeof(REAL4)*signal.length);
   correlation.data 	= (REAL4*) LALMalloc(sizeof(REAL4)*correlation.length);
   memset( &(coarseIn.shf), 0, sizeof(REAL8FrequencySeries) );
   coarseIn.shf.f0 	= 0;
   LALDCreateVector( &status, &(coarseIn.shf.data), randIn.psd.length );
   coarseIn.shf.deltaF 	= randIn.param.tSampling / signal.length;
   randIn.psd.data 	= (REAL8*) LALMalloc(sizeof(REAL8)*randIn.psd.length);

   /* --- Compute Noise Spectral Density --- */
   df = randIn.param.tSampling/(float) signal.length;
   LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, df);
   for (i=0; i<randIn.psd.length; i++)
     randIn.psd.data[i] = coarseIn.shf.data->data[i];

   /* --- And the bank of templates --- */
   LALInspiralCreateCoarseBank(&status, &list, &nlist, coarseIn);
   if (nlist==0) exit(0);					/* If no points in the bank then noting to do here */


   /* --- Let's print some optional comments : the position of templates in the bank ---*/
   if (!otherIn.quietFlag)
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
       fprintf(stdout, "#---------------------------------------------------------------\n");
       fprintf(stdout, "#Number of Coarse Bank Templates=%d\n",nlist);
       fprintf(stdout, "#Signal Length=%d, Number of sims=%d\n", signal.length, otherIn.ntrials);
       fprintf(stdout,"---------------------------------------------------------------\n");
       if (coarseIn.approximant==BCV)
	 {
	       fprintf(stdout, "#psi0Min=%e, psi0Max=%e, psi3Min=%e, psi3Max=%e, signal=%d, template=%d\n", 
		       coarseIn.psi0Min,coarseIn.psi0Max,coarseIn.psi3Min,coarseIn.psi3Max,otherIn.signal,randIn.param.approximant);
	       fprintf(stdout, "#psi0Tmplt   psi0Sig    psi3Tmplt    psi3Sig    fFinalTmplt   fFinalSig  m1 m2 Overlap/SNR\n");
	 }
       else
	 {
	   
	   fprintf(stdout, "#mMin=%e, mMax=%e, signal=%d, template=%d\n",
		       coarseIn.mMin,coarseIn.mMax,otherIn.signal,randIn.param.approximant);
	   fprintf(stdout, "#m1Tmplt   m1Sig    m2Tmplt    m2Sig    fFinalTmplt   fFinalSig  Overlap/SNR\n");
	 }
       fprintf(stdout,"---------------------------------------------------------------\n");
     }

   /* --- Estimate the plans --- */
   LALCreateForwardRealFFTPlan(&status, &fwdp, signal.length, 0);
   LALCreateReverseRealFFTPlan(&status, &revp, signal.length, 0);
   
   
   /* --- The overlap structure --- */
   overlapin.nBegin 	= 0;
   overlapin.nEnd 	= 0;
   overlapin.psd 	= randIn.psd;
   overlapin.fwdp 	= randIn.fwdp = fwdp;
   overlapin.revp 	= revp;
   
   
   
 
   /* --- The main loop --- */ 
   while (otherIn.ntrials--) 
     {
       randIn.param.approximant 	= otherIn.signal;  /* The waveform to inject */
       randIn.param.fCutoff 	= coarseIn.fUpper; /* its cutoff frequency */
       /* What kind of parameters do we use ?*/
       if (otherIn.signal == BCV) 
	 randIn.param.massChoice = psi0Andpsi3;
       else
	 randIn.param.massChoice = m1Andm2;
       /* Let's compute the random parameters of the waveform to inject*/
       for (i=0; i<signal.length; i++) signal.data[i] = 0.;
       LALRandomInspiralSignal(&status, &signal, &randIn);
       
       overlapin.signal = signal;
       jmax = 0;
       omax = 0.0;
       /* Process through the bank */
       for (j=0; j<nlist; j++) 
	 {
	   overlapin.param 	= list[j].params;
	   /**/
	   if (overlapin.param.approximant==BCV) 
	     overlapin.param.fCutoff = list[j].params.fendBCV ;
	   else
	     overlapin.param.fCutoff = randIn.param.fCutoff;
	   
	   for (i=0; i<signal.length; i++) correlation.data[i] = 0.;
	   LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin); /* The computation of the overlap */
	   list[j].params.fFinal = overlapin.param.fFinal;
	   
	   if (omax < overlapout.max) /*storing overlap here*/
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
		 list[jmax].params.fFinal, 
		 randIn.param.fFinal,
		 list[jmax].params.totalMass,
		 randIn.param.totalMass, 
		 list[jmax].params.eta,
		 randIn.param.eta, 
		 omax);
       fflush(stdout);
     }
   
   /* --- destroy the plans, correlation and signal --- */
   
   LALFree(randIn.psd.data);
   LALDDestroyVector( &status, &(coarseIn.shf.data) );
   LALFree(signal.data);
   LALFree(correlation.data);
   LALFree(list);
   LALDestroyRealFFTPlan(&status,&fwdp);   
   LALDestroyRealFFTPlan(&status,&revp);
   LALCheckMemoryLeaks();    
   
   return(0);

}



/* --- Print function --- */
void printf_timeseries (INT4 n, REAL4 *signal, double delta, double t0) 
{
  int i=0;

  do 
     printf ("%e %e\n", i*delta+t0, *(signal+i));
  while (n-++i); 
  printf("&\n");
}


/* --- Affect dummy values to a structure ---*/
void Init2DummyCoarse(InspiralCoarseBankIn *coarseIn)
{
  coarseIn->fLower    = -1;
  coarseIn->fUpper    = -1;
  coarseIn->tSampling = -1;
  coarseIn->order     = -1;
  coarseIn->space     = -1;
  coarseIn->mmCoarse  = -1;
  coarseIn->mmFine    = -1;
  coarseIn->iflso     = -1;
  coarseIn->mMin      = -1;
  coarseIn->mMax      = -1;
  coarseIn->MMax      = -1;
  coarseIn->massRange = -1;
  coarseIn->etamin    = -1;
  coarseIn->psi0Min   = +1;
  coarseIn->psi0Max   = -1;
  coarseIn->psi3Min   = +1;
  coarseIn->psi3Max   = -1;
  coarseIn->alpha     = -1;
  /*unsigned variable ... */
  coarseIn->numFcutTemplates = 0;    /* can't be negative */
  coarseIn->approximant      = 100;  /* can't be negative */
  coarseIn->order            = 100;  /* can't be negative */
}

/* --- Affect dummy values to a structure ---*/
void Init2DummyRandom(RandomInspiralSignalIn *randIn)
{
  randIn->useed             = 128092;
  randIn->type              = -1;
  randIn->SignalAmp         = -1;
  randIn->param.approximant = 100;   /* can't be negative*/
  randIn->param.order       = 100;   /* can't be negative */
  randIn->param.ieta        = -1; 
  randIn->param.startTime   = -1; 
  randIn->param.startPhase  = -1; 
  randIn->param.nStartPad   = -1;
  randIn->param.signalAmplitude = -1;
  randIn->param.nEndPad = -1;
  randIn->mMin          = -1;
  randIn->mMax          = -1;
  randIn->NoiseAmp      = -1;
  randIn->etaMin        = -1;

}

/* --- Affect dummy values to a structure ---*/
void Init2DummyOtherParam(OtherParamIn *otherIn)
{
  otherIn->template   = -1;
  otherIn->bank       = -1;
  otherIn->signal     = -1;
  otherIn->quietFlag  = -1;
  otherIn->ntrials    = -1;
}



/* --- Function to parse user parameters --- */
void 
ParseParameters(int *argc, 
		char **argv,
		InspiralCoarseBankIn *coarseIn,
		RandomInspiralSignalIn *randIn,
		OtherParamIn    *otherIn
		)
{
  INT4 i        = 1;

  while(i < *argc)
    {
      if (strcmp(argv[i],"-fl")==0)
	{
	  coarseIn->fLower = atof(argv[++i]);
	  randIn->param.fLower = coarseIn->fLower;
	}

      else if (strcmp(argv[i],"-tSampling")==0)
	{
	  coarseIn->tSampling  = atof(argv[++i]);
	  randIn->param.tSampling = coarseIn->tSampling;
	}
      else if (strcmp(argv[i],"-mMin")==0)
	{
	  coarseIn->mMin = atof(argv[++i]); 
	  randIn->mMin = coarseIn->mMin;
	  randIn->param.mass1 = randIn->mMin;
	  randIn->param.mass2 = randIn->mMin;
	}
      else if (strcmp(argv[i],"-mMax")==0)
	{
	  coarseIn->mMax = atoi(argv[++i]);
	  coarseIn->MMax = 2. * coarseIn->mMax;
	  randIn->mMax   = coarseIn->mMax;
	  randIn->MMax   = 2.*randIn->mMax;	 
	}
      else if (strcmp(argv[i],"-psi0Max")==0)
	{
	  coarseIn->psi0Max =atof(argv[++i]);
	  randIn->psi0Max = coarseIn->psi0Max;
	}
      else if (strcmp(argv[i],"-psi3Min")==0)
	{
	  coarseIn->psi3Min =atof(argv[++i]);
	  randIn->psi3Min = coarseIn->psi3Min;
	}
      else if (strcmp(argv[i],"-psi0Min")==0)
	{
	  coarseIn->psi0Min =atof(argv[++i]);
	  randIn->psi0Min = coarseIn->psi0Min;
	}
      else if (strcmp(argv[i],"-psi3Max")==0)
	{
	  coarseIn->psi3Max =atof(argv[++i]);
	  randIn->psi3Max = coarseIn->psi3Max;
	}      
      else if (strcmp(argv[i],"-mm")==0)
	{
	  coarseIn->mmCoarse =atof(argv[++i]);

	}
      else if (strcmp(argv[i],"-numFcut")==0)
	coarseIn->numFcutTemplates = atof(argv[++i]); 
      else if (strcmp(argv[i],"-simType")==0)
	randIn->type = atoi(argv[++i]);
      else if (strcmp(argv[i],"-sigAmp")==0)
	randIn->SignalAmp = atof(argv[++i]);
      else if (strcmp(argv[i],"-alpha")==0)
	{
	  coarseIn->alpha = atof(argv[++i]); 
	  randIn->param.alpha = coarseIn->alpha;
	}
      else if (strcmp(argv[i],"-order")==0)
	randIn->param.order = atoi(argv[++i]); 
      else if (strcmp(argv[i],"-n")==0)
	otherIn->ntrials = atoi(argv[++i]);       
      else if (strcmp(argv[i],"-seed")==0)
	randIn->useed = atoi(argv[++i])+1;
      else if (strcmp(argv[i],"-quiet")==0)
	otherIn->quietFlag = 1;
      else if (strcmp(argv[i], "-template")==0)
	{
	  if (strcmp(argv[++i],"TaylorT1")==0)
	    otherIn->template = 0;
	  else if (strcmp(argv[i],"TaylorT2")==0)
	    otherIn->template = 1;
	  else if (strcmp(argv[i],"TaylorT3")==0)
	    otherIn->template = 2;
	  else if (strcmp(argv[i],"TaylorF1")==0)
	    otherIn->template = 3;
	  else if (strcmp(argv[i],"TaylorF2")==0)
	    otherIn->template = 4;
	  else if (strcmp(argv[i],"PadeT1")==0)
	    otherIn->template = 5;
	  else if (strcmp(argv[i],"PadeF1")==0)
	    otherIn->template = 6;
	  else if (strcmp(argv[i],"EOB")==0)
	    otherIn->template = 7;
	  else if (strcmp(argv[i],"BCV")==0)
	    otherIn->template = 8;
	  else if (strcmp(argv[i],"SpinTaylorT3")==0)
	    otherIn->template = 9;
	  coarseIn->approximant = otherIn->template;	

	}	
      else if ( (strcmp(argv[i], "-signal")==0))
	{
	  if (strcmp(argv[++i],"TaylorT1")==0)
	    otherIn->signal = 0;
	  else if (strcmp(argv[i],"TaylorT2")==0)
	    otherIn->signal = 1;
	  else if (strcmp(argv[i],"TaylorT3")==0)
	    otherIn->signal = 2;
	  else if (strcmp(argv[i],"TaylorF1")==0)
	    otherIn->signal = 3;
	  else if (strcmp(argv[i],"TaylorF2")==0)
	    otherIn->signal = 4;
	  else if (strcmp(argv[i],"PadeT1")==0)
	    otherIn->signal = 5;
	  else if (strcmp(argv[i],"PadeF1")==0)
	    otherIn->signal = 6;
	  else if (strcmp(argv[i],"EOB")==0)
	    otherIn->signal = 7;
	  else if (strcmp(argv[i],"BCV")==0)
	    otherIn->signal = 8;
	  else if (strcmp(argv[i],"SpinTaylorT3")==0)
	    otherIn->signal = 9;
	  randIn->param.approximant = otherIn->signal;	
	  
	}	
      else 
	{
	  Help(*coarseIn, *randIn, *otherIn);
	  exit(0);
	}
      
      i++;       
    } 
}





/* --- Function to check validity of parsed parameters (TODO)--- */
void CheckParams(InspiralCoarseBankIn coarseIn,
		 RandomInspiralSignalIn randIn,
		 OtherParamIn otherIn)
{

  if (coarseIn.psi0Min > coarseIn.psi0Max
      || coarseIn.psi3Min >coarseIn.psi3Max
      || otherIn.ntrials <=0 )   
    {
      Help(coarseIn, randIn, otherIn);
      exit(0);
    }
}


/* --- Set default values to parameters --- */
void SetDefault(InspiralCoarseBankIn *coarseIn,
		RandomInspiralSignalIn *randIn, 
		OtherParamIn *otherIn)
{

  if (coarseIn->fLower == -1)
    {
      coarseIn->fLower     = BANKEFFICIENCY_FLOWER; 
      randIn->param.fLower = coarseIn->fLower;
    }
  if (coarseIn->fUpper == -1)
    {
      if (coarseIn->tSampling == -1)
	{
	  coarseIn->tSampling     = BANKEFFICIENCY_TSAMPLING;
	  randIn->param.tSampling = coarseIn->tSampling;
	}       
      coarseIn->fUpper       =  coarseIn->tSampling/2. - 1;
      coarseIn->fUpper       =  BANKEFFICIENCY_FUPPER;
      randIn->param.fCutoff  =  coarseIn->fUpper;
    }
  else if (coarseIn->tSampling == -1)
    {
      coarseIn->tSampling     = BANKEFFICIENCY_TSAMPLING;
      randIn->param.tSampling = coarseIn->tSampling;
    } 
  /*ajout order 1 et 2*/
  if (coarseIn->order == 100 )
    {
      coarseIn->order       = BANKEFFICIENCY_ORDER;
      randIn->param.order  = coarseIn->order; 
    }
  if (coarseIn->space == -1)
    {
      coarseIn->space      = BANKEFFICIENCY_SPACE ;
    }
  if(coarseIn->mmCoarse == -1)
    {
      coarseIn->mmCoarse   = BANKEFFICIENCY_MMCOARSE ;
    }
  if(coarseIn->mmFine == -1)
    {
      coarseIn->mmFine     = BANKEFFICIENCY_MMFINE ;
    }
  if(coarseIn->iflso == -1)
    {
      coarseIn->iflso     = BANKEFFICIENCY_IFLSO ;
    }
  if(coarseIn->mMin == -1)
    {
      coarseIn->mMin     = BANKEFFICIENCY_MMIN ;
      randIn->mMin       = coarseIn->mMin;
      /*In order to allocate memory, we need the longest template*/
      randIn->param.mass1 = coarseIn->mMin;
      randIn->param.mass2 = coarseIn->mMin;
    }
  if (coarseIn->mMax == -1)
    {
      coarseIn->mMax     = BANKEFFICIENCY_MMAX ;  
      coarseIn->MMax     = coarseIn->mMax * 2.;
      randIn->mMax       = coarseIn->mMax;
      randIn->MMax       = coarseIn->mMax * 2;

    }
  coarseIn->massRange = MinMaxComponentMass; 
  coarseIn->etamin = coarseIn->mMin 
    * ( coarseIn->MMax - coarseIn->mMin) / pow(coarseIn->MMax,2.);
  randIn->etaMin = coarseIn->etamin;

  if (coarseIn->psi0Min == 1)
    {
      coarseIn->psi0Min   = BANKEFFICIENCY_PSI0MIN;
      randIn->psi0Min     = coarseIn->psi0Min;
    }
  if (coarseIn->psi0Max == -1)
    {
      coarseIn->psi0Max   = BANKEFFICIENCY_PSI0MAX;
      randIn->psi0Max     = coarseIn->psi0Max;
    }
  if (coarseIn->psi3Min == 1)
    {
      coarseIn->psi3Min   = BANKEFFICIENCY_PSI3MIN;
      randIn->psi3Min     = coarseIn->psi3Min;
    }
  if (coarseIn->psi3Max == -1)
    {
      coarseIn->psi3Max   = BANKEFFICIENCY_PSI3MAX;
      randIn->psi3Max     = coarseIn->psi3Max;
    }
  if (coarseIn->numFcutTemplates == 0)
    {
      coarseIn->numFcutTemplates = BANKEFFICIENCY_NFCUT;
    }
  if(coarseIn->approximant == 100)
    {
      coarseIn->approximant = BANKEFFICIENCY_BANK;
      otherIn->template     = coarseIn->approximant; 
    }
  if(randIn->param.approximant  == 100)
    {
      randIn->param.approximant = BANKEFFICIENCY_SIGNAL;
      otherIn->signal = randIn->param.approximant;
    }
  if(coarseIn->alpha  == -1)
    {
      coarseIn->alpha = BANKEFFICIENCY_ALPHA;
      randIn->param.alpha = coarseIn->alpha;
    }
  if(randIn->param.fendBCV  == -1)
    {
      randIn->param.fendBCV = BANKEFFICIENCY_FENDBCV;
    }


  /* OtherIn structure */
  if (otherIn->quietFlag == -1)
    {
      otherIn->quietFlag = BANKEFFICIENCY_QUIETFLAG;
    }
  if (otherIn->ntrials == -1)
    {
      otherIn->ntrials = BANKEFFICIENCY_NTRIALS;
    }

  /*no parsing for those parameters*/

  randIn->type                  = BANKEFFICIENCY_TYPE;
  randIn->SignalAmp             = BANKEFFICIENCY_SIGNALAMP;
  randIn->param.ieta            = BANKEFFICIENCY_IETA; 
  randIn->param.startTime       = BANKEFFICIENCY_STARTTIME; 
  randIn->param.startPhase      = BANKEFFICIENCY_STARTPHASE; 
  randIn->param.nStartPad       = BANKEFFICIENCY_NSTARTPHASE;
  randIn->param.signalAmplitude = BANKEFFICIENCY_SIGNALAMPLITUDE;
  randIn->param.nEndPad         = BANKEFFICIENCY_NENDPAD;
  randIn->NoiseAmp              = BANKEFFICIENCY_NOISEAMP;
}




/* --- Documenation (TODO) --- */
void Help(InspiralCoarseBankIn   coarseIn,
		 RandomInspiralSignalIn randIn,
		 OtherParamIn           otherIn)
{
  INT4 temp;

  fprintf(stderr,"\nUSAGE:  [options]\n");
  fprintf(stderr,"The options are (with default values in brackets)\n");
  fprintf(stderr,"All options should be followed by a number except -quiet\n");
  fprintf(stderr,"      -quiet : if this flag is present, the output is restricted to the minimum\n");
  fprintf(stderr,"          -n : number of trials                   (%d)\n",       BANKEFFICIENCY_NTRIALS);
  /*  fprintf(stderr,"       -seed : seed for random generation         (%d)\n",       randIn.useed);*/
  fprintf(stderr,"    -simType : type of simulation, 0, 1 or 2      (%d)\n\n",     BANKEFFICIENCY_TYPE);
  fprintf(stderr,"         -fl : lower frequency cutoff             (%7.2f) Hz\n", BANKEFFICIENCY_FLOWER);
  fprintf(stderr,"       -mMin : minimal mass of component stars    (%7.2f) Mo\n", BANKEFFICIENCY_MMIN);
  fprintf(stderr,"       -mMax : maximal mass of component stars    (%7.2f) Mo\n", BANKEFFICIENCY_MMAX);
  fprintf(stderr,"      -x0Max : Max value of psi0 (%7.2f)\n",                     BANKEFFICIENCY_PSI0MAX);
  fprintf(stderr,"      -x1Min : Min value of psi3 (%7.2f)\n",                     BANKEFFICIENCY_PSI3MIN);
  fprintf(stderr,"      -alpha : BCV amplitude correction parameter (%7.2f)\n",    BANKEFFICIENCY_ALPHA);
  fprintf(stderr,"    -numFcut : number of layers in Fcut dimension (%7.2d)\n",    BANKEFFICIENCY_NFCUT);
  fprintf(stderr,"         -mm : minimal match for template bank    (%7.3f)\n",    BANKEFFICIENCY_MMCOARSE);
  fprintf(stderr,"   -template : PN model for template bank, e.g. BCV (%d)\n\n",   BANKEFFICIENCY_BANK);
  fprintf(stderr,"-approximant : PN model for Monte-Carlo, e.g. EOB (%d)\n",       BANKEFFICIENCY_SIGNAL);
  fprintf(stderr,"     -signal : same as -approximant\n");
  fprintf(stderr,"      -order : order of PN model                  (%7.2d)\n",    BANKEFFICIENCY_ORDER);;
  fprintf(stderr,"     -sigAmp : amplitude of the signal            (%7.2f)\n",    BANKEFFICIENCY_SIGNALAMP);

  /*to avoid boring warning*/
  temp = coarseIn.approximant;
  temp =  randIn.useed;
  temp =  otherIn.signal;
}
