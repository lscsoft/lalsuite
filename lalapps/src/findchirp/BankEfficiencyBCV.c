
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
#include <lal/LALNoiseModels.h>
#include <lal/LALInspiralBank.h>
#include <lal/RealFFT.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>


NRCSID( BANKEFFICIENCYC, "$Id$");


#define BANKEFFICIENCY_ENORM  0
#define BANKEFFICIENCY_ESUB   1
#define BANKEFFICIENCY_EARG   2
#define BANKEFFICIENCY_EVAL   3
#define BANKEFFICIENCY_EFILE  4
#define BANKEFFICIENCY_EINPUT 5
#define BANKEFFICIENCY_EMEM   6

#define BANKEFFICIENCY_MSGENORM  "Normal exit"
#define BANKEFFICIENCY_MSGESUB   "Subroutine failed"
#define BANKEFFICIENCY_MSGEARG   "Error parsing arguments"
#define BANKEFFICIENCY_MSGEVAL   "Input argument out of valid range"
#define BANKEFFICIENCY_MSGEFILE  "Could not open file"
#define BANKEFFICIENCY_MSGEINPUT "Error reading file"
#define BANKEFFICIENCY_MSGEMEM   "Out of memory"


#define BANKEFFICIENCY_FLOWER       	  40.
#define BANKEFFICIENCY_TSAMPLING    	2048.
#define BANKEFFICIENCY_FUPPER       	1000.
#define BANKEFFICIENCY_ORDER_SIGNAL     twoPN
#define BANKEFFICIENCY_ORDER_TEMPLATE   twoPN
#define BANKEFFICIENCY_MMCOARSE     	0.8
#define BANKEFFICIENCY_MMFINE       	0.97
#define BANKEFFICIENCY_MMIN            	5.
#define BANKEFFICIENCY_MMAX           	20.
#define BANKEFFICIENCY_ALPHA           	0.
#define BANKEFFICIENCY_FENDBCV       	300.
#define BANKEFFICIENCY_SPACE    	Psi0Psi3
#define BANKEFFICIENCY_IFLSO           	0.
#define BANKEFFICIENCY_PSI0MIN        	10.
#define BANKEFFICIENCY_PSI0MAX    	250000.
#define BANKEFFICIENCY_PSI3MIN     	-2200.
#define BANKEFFICIENCY_PSI3MAX       	-10.
#define BANKEFFICIENCY_NFCUT           	5
#define BANKEFFICIENCY_SIGNAL  		TaylorT1
#define BANKEFFICIENCY_BANK          	BCV
#define BANKEFFICIENCY_TEMPLATE        	BCV
#define BANKEFFICIENCY_TYPE            	0
#define BANKEFFICIENCY_SIGNALAMP      	10.
#define BANKEFFICIENCY_IETA            	1.
#define BANKEFFICIENCY_STARTTIME       	0.
#define BANKEFFICIENCY_STARTPHASE    	0.
#define BANKEFFICIENCY_NSTARTPHASE  	1000
#define BANKEFFICIENCY_SIGNALAMPLITUDE 	10.
#define BANKEFFICIENCY_NENDPAD         	0
#define BANKEFFICIENCY_NOISEAMP        	1.
#define BANKEFFICIENCY_NTRIALS         	2
#define BANKEFFICIENCY_SEED        122888
/* Other Parameters */
#define BANKEFFICIENCY_QUIETFLAG       	        0
#define BANKEFFICIENCY_AMBIGUITYFUNCTION      	0
#define BANKEFFICIENCY_NORMALISATION      	0
#define BANKEFFICIENCY_FMAXIMIZATION      	0
#define BANKEFFICIENCY_PRINTOVERLAP             0
#define BANKEFFICIENCY_PRINTFILTER              0
#define BANKEFFICIENCY_CHECK                    0

typedef struct{
	double a11,a21,a22;
}BCVMaximizationMatrix;
 


typedef enum {
  InQuadrature,
  AlphaMaximization
} OverlapMethodIn;

typedef struct	
{
  INT4 signal; /*random signal*/
  INT4 template;/*template and bank*/
  INT4 bank;
  INT4 ntrials;
  INT4 quietFlag;
  INT4 ambiguityFunction;
  INT4 normalisation;
  INT4 FMaximization;
  INT4 PrintOverlap;
  INT4 PrintFilter;
  OverlapMethodIn overlapMethod;
  INT4 check;
  double m1,m2, psi0,psi3;
  char *inputPSD;
} OtherParamIn;

void LALCreateMomentVector(LALStatus             *status,
			   REAL4Vector           *a11,
			   REAL4Vector           *a21,
			   REAL4Vector           *a22,
			   REAL8FrequencySeries  *psd,
			   InspiralTemplate      *params);





void 
LALCreateVectorFreqPower(
				 REAL4Vector *vector,
				 InspiralTemplate params,
				 int a,
				 int b);


void LALMoments(LALStatus         *status,
		REAL8             *moment,
		InspiralMomentsIn pars);

void LALGenerateWaveform(LALStatus *status, REAL4Vector *signal, RandomInspiralSignalIn *params);
void LALGetOrthoNormalFilter(REAL4Vector *filter2, REAL4Vector filter1);

void
LALWaveOverlapBCV(
		LALStatus *status,
   		REAL4Vector             *correlation,
		InspiralWaveOverlapOut  *overlapout,
		InspiralWaveOverlapIn   *overlapin,
		REAL4Vector              VectorPhase, 
		REAL4Vector              VectorAmplitude1,
		REAL4Vector              VectorAmplitude2,
		BCVMaximizationMatrix    matrix,
		OtherParamIn otherIn
		);




void
LALCreateBCVTemplate(
		REAL4Vector *signal,
		REAL4Vector VectorAmplitude,
		REAL4Vector VectorPhase,
		InspiralTemplate params);


void
LALMatrix(
      		LALStatus               *status,
      		BCVMaximizationMatrix   *output,
      		REAL8FrequencySeries    *psd,
      		InspiralTemplate        params );



void 
LALCreateBCVPhase(
		REAL4Vector *VectorPhase,
		REAL4Vector *VectorPowerFm5_3,
		REAL4Vector *VectorPowerFm2_3,
		InspiralTemplate params);



void 
LALBCVAmplitude1(
		 REAL4Vector *VectorAmplitude,
		 REAL4Vector *VectorPowerFm7_6,
		 InspiralTemplate params,
		 REAL8 cte);
void
LALBCVAmplitude2(
		 REAL4Vector *VectorAmplitude,
		 REAL4Vector VectorPowerFm7_6,
		 REAL4Vector VectorPowerFm1_2,
		 InspiralTemplate params,
		 REAL8 c1, REAL8 c2);

void
printf_timeseries (
		   INT4 n,
		   REAL4 *signal,
		   REAL8 delta,
		   REAL8 t0);

void
Init2DummyCoarse(
		InspiralCoarseBankIn *coarseIn);

void
Init2DummyRandom(
		RandomInspiralSignalIn *randIn);

void
Init2DummyOtherParam(
		OtherParamIn *otherIn);

void
ParseParameters(
		int 			*argc, 
		char 			**argv,
		InspiralCoarseBankIn 	*coarseIn,
		RandomInspiralSignalIn 	*randIn,
		OtherParamIn 		*otherIn);		     

void 	
CheckParams(
		InspiralCoarseBankIn 	coarseIn,
		RandomInspiralSignalIn 	randIn,
		OtherParamIn 		otherIn);

void 
SetDefault(
		InspiralCoarseBankIn 	*coarseIn, 
		RandomInspiralSignalIn 	*randIn,
		OtherParamIn 		*otherIn);

void 	
Help(
		InspiralCoarseBankIn 	coarseIn,
    		RandomInspiralSignalIn 	randIn,
   		OtherParamIn 		otherIn);


int
main (int argc, char **argv ) 
{
lalDebugLevel=1;

   /* --- Variables ---*/
  double f,frequency, frequencyMax=0, overlapMaxAlpha=0;
  INT4	i, j=0,n, k;
   INT4 jmax, nlist;
   REAL8 df, omax;
   REAL4Vector  signal, correlation;
   void   			*noisemodel = LALLIGOIPsd;
   /*   void   			*noisemodel = LALVIRGOPsd;*/
   RandomInspiralSignalIn	 randIn;			/* random signal waveform to inject*/
   InspiralWaveOverlapIn 	 overlapin;			/* structure for Overlap*/
   InspiralWaveOverlapOut        overlapout, 
			         overlapoutmax;
   OtherParamIn 		 otherIn;			/* personal structure to parse some independant parameters*/
   RealFFTPlan  	        *fwdp=NULL,
		     	        *revp=NULL;
   InspiralTemplateList         *list=NULL;
   InspiralCoarseBankIn          coarseIn; 			/* strcture for the bank of templates*/

   REAL4Vector                   VectorAmplitude1, VectorAmplitude2, VectorPhase;
   REAL4Vector                   VectorPowerFm5_3,  VectorPowerFm2_3;
   REAL4Vector                   VectorPowerFm7_6,  VectorPowerFm1_2;

   REAL4Vector                   VectorA11, VectorA21, VectorA22;

   static LALStatus status;
   BCVMaximizationMatrix matrix, bestMatrix;   

   
  
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
   
   
   for (i=0; i< (int)randIn.psd.length; i++)
     {
       randIn.psd.data[i] = coarseIn.shf.data->data[i];
     }
   


   /* --- And the bank of templates --- */
    LALInspiralCreateCoarseBank(&status, &list, &nlist, coarseIn);

   if (nlist==0) exit(0);	
				/* If no points in the bank then noting to do here */
   if (!otherIn.quietFlag)
     fprintf(stderr,"Bank Generated\n");

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
		       coarseIn.psi0Min,coarseIn.psi0Max,coarseIn.psi3Min,coarseIn.psi3Max,otherIn.signal,otherIn.template);
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
   
   /* If ambiguity function flag equal to 1 we don't want to compute 
      tons of trials since th number of ouput could be too important*/
  

   if (otherIn.ambiguityFunction  == 1)
     {
       otherIn.ntrials=1;
     }

   list[1].params.startPhase =0;   
  

   
   /* Then we can  compute the frequency vector*/
   VectorPhase.length = VectorAmplitude1.length = VectorAmplitude2.length= signal.length/2;
   VectorPhase.data = (REAL4*) LALMalloc(sizeof(REAL4) * VectorPhase.length);

   VectorAmplitude1.data = (REAL4*) LALMalloc(sizeof(REAL4) * VectorAmplitude1.length);
   VectorAmplitude2.data = (REAL4*) LALMalloc(sizeof(REAL4) * VectorAmplitude2.length);
   

   VectorA11.length = VectorA22.length = VectorA21.length = signal.length / 2 ;

   VectorA11.data = (REAL4*) LALMalloc(sizeof(REAL4) * VectorA11.length);
   VectorA21.data = (REAL4*) LALMalloc(sizeof(REAL4) * VectorA21.length);
   VectorA22.data = (REAL4*) LALMalloc(sizeof(REAL4) * VectorA22.length);
   


   VectorPowerFm2_3.length  = VectorPowerFm5_3.length = signal.length/2;       
   VectorPowerFm7_6.length  = VectorPowerFm1_2.length = signal.length/2;       
   VectorPowerFm2_3.data = (REAL4*) LALMalloc(sizeof(REAL4) * VectorPowerFm2_3.length);
   VectorPowerFm5_3.data = (REAL4*) LALMalloc(sizeof(REAL4) * VectorPowerFm5_3.length);
   VectorPowerFm7_6.data = (REAL4*) LALMalloc(sizeof(REAL4) * VectorPowerFm7_6.length);
   VectorPowerFm1_2.data = (REAL4*) LALMalloc(sizeof(REAL4) * VectorPowerFm1_2.length);

   n =VectorAmplitude1.length;
   df = randIn.param.tSampling /(REAL8)n/2.;      

   LALCreateVectorFreqPower(&VectorPowerFm5_3, randIn.param, -5, 3);
   LALCreateVectorFreqPower(&VectorPowerFm2_3, randIn.param, -2, 3);
   LALCreateVectorFreqPower(&VectorPowerFm7_6, randIn.param, -7, 6);
   LALCreateVectorFreqPower(&VectorPowerFm1_2, randIn.param, -1, 2);
   
   LALCreateMomentVector(&status, &VectorA11, &VectorA21, &VectorA22, &coarseIn.shf, &(list[j].params));
   

   /* --- The main loop --- */ 
   while (otherIn.ntrials--) 
     {
       
       /* Generate Random waveform here*/
       randIn.param.approximant    = otherIn.signal;  /* The waveform to inject */
       randIn.param.fCutoff 	   = coarseIn.fUpper; /* its cutoff frequency */
       /* What kind of parameters do we use ?*/
       if (otherIn.signal == BCV) 
	 randIn.param.massChoice = psi0Andpsi3;
       else
	 randIn.param.massChoice = m1Andm2;
       /* Let's compute the random parameters of the waveform to inject*/
       for (i=0; i<(INT4)signal.length; i++) signal.data[i] = 0.;       

       /* we might force to compute a non random waveform giving the two
	  input parameters m1 and m2 */
       if (otherIn.m1==-1 && otherIn.m2 ==-1)
	 LALRandomInspiralSignal(&status, &signal, &randIn);
       else
	 {
	   randIn.param.mass1 = otherIn.m1;
	   randIn.param.mass2 = otherIn.m2;
	   LALGenerateWaveform(&status,&signal, &randIn);
	 }
       



       overlapin.signal = signal;
       jmax = 0;
       omax = 0.0;

       if(otherIn.check == 1)
	 {
	   nlist =1;	   
	 }

       /* we might force to use only one template giving the input psi0 and psi3. 
	  Then, we don't need to use the bank*/
       if (otherIn.psi0 != -1 && otherIn.psi3 != -1)
	 {
	   list[0].params.psi0 = otherIn.psi0;
	   list[0].params.psi3 = otherIn.psi3;
	   nlist  = 1;
	 }

       /* Process hrough the bank */      
       for (j=0; j<nlist; j++) 
	 {	   
	   overlapin.param 	= list[j].params;	   
	   
	   switch(otherIn.overlapMethod)
	     {
	     case AlphaMaximization:	   	   	   
	       df    = randIn.param.tSampling / (REAL8)(VectorAmplitude1.length)/2.;
	       
	       if (otherIn.check==1)
		 {
		   list[j].params.psi0    = randIn.param.psi0;
		   list[j].params.psi3    = randIn.param.psi3;
		   list[j].params.fendBCV = randIn.param.fFinal;		   
		   list[j].params.fFinal  = randIn.param.fFinal;		   
		   list[j].params.fCutoff = randIn.param.fFinal;		   
		 }
	       else
		 {
		   list[j].params.fFinal    =	 
		     list[j].params.fCutoff =
		     list[j].params.fendBCV ; 	   	   
		 }
	       
	       overlapin.param.fFinal  = list[j].params.fendBCV;
	       overlapin.param.fendBCV = list[j].params.fendBCV;
	       overlapin.param.fCutoff = list[j].params.fendBCV;
	       
	       
	       
	       
	       k = floor(list[j].params.fendBCV / df);
	       matrix.a21 = VectorA21.data[k];
	       matrix.a11 = VectorA11.data[k];
	       matrix.a22 = VectorA22.data[k];
	       
	       /*
		 
	       
	       LALMatrix(&status,
	       &matrix,
	       &coarseIn.shf,
	       list[j].params);
	       */
	       
	       for (i=0; i< n; i++)
		 {
		   f = i*df;
		   if (f < randIn.param.fLower  
		       || f >list[j].params.fendBCV)	 
		     {
		       VectorAmplitude1.data[i] = 0.;
		       VectorAmplitude2.data[i] = 0.;
		     }
		   else
		     {
		       VectorAmplitude2.data[i] = matrix.a21 * VectorPowerFm7_6.data[i]
			 + matrix.a22 * VectorPowerFm1_2.data[i];
		       VectorAmplitude1.data[i] = matrix.a11 * VectorPowerFm7_6.data[i];
		     }
		 }   
	       
	       LALCreateBCVPhase(&VectorPhase,
				 &VectorPowerFm5_3, 
				 &VectorPowerFm2_3,
				 list[j].params);	  	   
	       
	       overlapin.param.fFinal = list[j].params.fendBCV;
	       overlapin.param.fendBCV = list[j].params.fendBCV;
	       overlapin.param.fCutoff = list[j].params.fendBCV;
	       
	       
	       LALWaveOverlapBCV(&status,
				 &correlation,
				 &overlapout,
				 &overlapin,
				 VectorPhase,
				 VectorAmplitude1,
				 VectorAmplitude2,
				 matrix,
				 otherIn);		 	       
	       break;
	     case InQuadrature:
	       overlapin.param.fFinal = randIn.param.fCutoff;
	       overlapin.param.fendBCV =randIn.param.fCutoff; 
	       overlapin.param.fCutoff = randIn.param.fCutoff;
	       
	       for (i=0; i<(INT4)signal.length; i++) correlation.data[i] = 0.;	   	   	  
	       
	       
	       LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
	       overlapout.alpha = -1;	       
	       break;	     
	     }
	 
	 
 
       
       if (otherIn.ambiguityFunction  == 1 && overlapout.max!=0)
	 printf("%lf %lf %lf %lf %lf\n", list[j].params.psi0, list[j].params.psi3, overlapout.max,randIn.param.fendBCV, list[j].params.fendBCV);
       
       if (omax < overlapout.max) /*storing overlap here*/
	 {
	   omax = overlapout.max;
	   overlapoutmax = overlapout;
	   jmax = j;
	   bestMatrix = matrix; 
	 }
     }

   
   LALInspiralParameterCalc(&status, &list[jmax].params);
   if (coarseIn.approximant==BCV)
     {
	   fprintf(stdout, "%e %e %e %e %e  %e %e %e %e %e %e %e %e %e %e %e %e\n", 
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
		   omax,
		   overlapout.phase,
		   overlapout.alpha,
		   overlapout.alpha*pow(list[jmax].params.fFinal,2./3.), bestMatrix.a11, bestMatrix.a21, bestMatrix.a22);


	   if (otherIn.FMaximization)
	     {
	       for (frequency = 300. ; frequency < 600; frequency+=5)
		 {	   
		   /* we need fcutoff in LALMatrix*/
		   list[jmax].params.fCutoff = frequency;
		   list[jmax].params.fendBCV = frequency;
		   list[jmax].params.fFinal = frequency;
		   LALMatrix(&status,
			     &matrix,
			     &coarseIn.shf,
				       list[jmax].params);

		   for (i=0; i< n; i++)
		     {
		       f = i*df;
		       if (f < randIn.param.fLower  || randIn.psd.data[i]==0 || f >frequency)	 
			 {
			   VectorAmplitude1.data[i] = 0.;
			   VectorAmplitude2.data[i] = 0.;
			 }
		       else
			 {
			   VectorAmplitude2.data[i] = matrix.a21 * pow(f,-7./6.) 
			     + matrix.a22 * pow(f,-1./2.);
			   VectorAmplitude1.data[i] = matrix.a11 * pow(f,-7./6.);		   		   		   
			 }
		     }
		   
		   LALCreateBCVPhase(&VectorPhase,
				     &VectorPowerFm5_3, 
				     &VectorPowerFm2_3,
				     list[jmax].params);	  	   
		   
		   
		   
		   
		   overlapin.param.fFinal = list[jmax].params.fendBCV;
		   overlapin.param.fendBCV = list[jmax].params.fendBCV;
		   overlapin.param.fCutoff = list[jmax].params.fendBCV;
		   
		   otherIn.PrintOverlap = 0;
		   otherIn.PrintFilter = 0;
		   LALWaveOverlapBCV(&status,
				     &correlation,
				     &overlapout,
				     &overlapin,
				     VectorPhase,
				     VectorAmplitude1,
				     VectorAmplitude2,
				     matrix,
				     otherIn);
		   
	       fprintf(stdout, "##%e  %e %e %e %e %e %e %e %e %e %e %e %e %e\n", 
		       list[jmax].params.psi0, 
		       randIn.param.psi0, 
		       list[jmax].params.psi3,
		       randIn.param.psi3, 
		       frequency,
		       randIn.param.fFinal,
		       list[jmax].params.totalMass,
		       randIn.param.totalMass,
		       randIn.param.mass1,
		       randIn.param.mass2,
		       overlapout.max, 
		       overlapout.phase,
		       overlapout.alpha,
		       overlapout.alpha*pow(list[jmax].params.fFinal,2./3.));
		       

		   
		   if (omax < overlapout.max) /*storing overlap here*/
		     {
		       omax = overlapout.max;
		       overlapoutmax = overlapout;
		       overlapMaxAlpha = overlapout.alpha;
		       frequencyMax = frequency;
		     }
		 }
	       fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e\n", 
		       list[jmax].params.psi0, 
		       randIn.param.psi0, 
		       list[jmax].params.psi3,
		       randIn.param.psi3, 
		       frequencyMax,
		       randIn.param.fFinal,
		       list[jmax].params.totalMass,
		       randIn.param.totalMass,
		       randIn.param.mass1,
		       randIn.param.mass2,
		       omax, overlapMaxAlpha * pow(frequencyMax,2./3.));
	     }

	 }
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

   LALFree(VectorPowerFm5_3.data);
   LALFree(VectorPowerFm2_3.data);
   LALFree(VectorPowerFm1_2.data);
   LALFree(VectorPowerFm7_6.data);
   LALFree(VectorAmplitude1.data);
   LALFree(VectorAmplitude2.data);
   LALFree(VectorPhase.data);

      LALFree(VectorA11.data);
   LALFree(VectorA21.data);
   LALFree(VectorA22.data);
   
   
   LALFree(randIn.psd.data);
   LALDDestroyVector( &status, &(coarseIn.shf.data) );

   LALFree(signal.data);
   LALFree(correlation.data);
   LALFree(list);

   LALDestroyRealFFTPlan(&status,&fwdp);
   LALDestroyRealFFTPlan(&status,&revp);
   LALCheckMemoryLeaks();    


return 0;
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
  coarseIn->fLower    = BANKEFFICIENCY_FLOWER;
  coarseIn->fUpper    = -1;
  coarseIn->tSampling = -1;
  coarseIn->space     = BANKEFFICIENCY_SPACE;
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
  coarseIn->approximant      = BANKEFFICIENCY_TEMPLATE;  
  coarseIn->order            = BANKEFFICIENCY_ORDER_TEMPLATE;  
}

/* --- Affect dummy values to a structure ---*/
void Init2DummyRandom(RandomInspiralSignalIn *randIn)
{
  randIn->useed             = -1;
  randIn->type              = -1;
  randIn->SignalAmp         = -1;
  randIn->param.approximant = BANKEFFICIENCY_SIGNAL;   /* can't be negative*/
  randIn->param.order       = BANKEFFICIENCY_ORDER_SIGNAL;   /* can't be negative */
  randIn->param.ieta        = -1; 
  randIn->param.startTime   = -1; 
  randIn->param.startPhase  = -1; 
  randIn->param.nStartPad   = -1;
  randIn->param.signalAmplitude = -1;
  randIn->param.nEndPad = -1;
  randIn->param.fLower  = BANKEFFICIENCY_FLOWER;
  randIn->mMin          = -1;
  randIn->mMax          = -1;
  randIn->NoiseAmp      = -1;
  randIn->etaMin        = -1;

}

/* --- Affect dummy values to a structure ---*/
void Init2DummyOtherParam(OtherParamIn *otherIn)
{
  otherIn->template      = BANKEFFICIENCY_TEMPLATE;
  otherIn->bank          = -1;
  otherIn->signal        = BANKEFFICIENCY_SIGNAL;
  otherIn->quietFlag     = -1;
  otherIn->ntrials       = -1;
  otherIn->ambiguityFunction = -1;
  otherIn->m1            = -1;
  otherIn->m2            = -1;
  otherIn->psi0          = -1;
  otherIn->psi3          = -1;
  otherIn->normalisation = -1;  
  otherIn->FMaximization = -1;
  otherIn->PrintOverlap  = -1;
  otherIn->PrintFilter  = -1;
  otherIn->overlapMethod  = AlphaMaximization;
  otherIn->check          = -1;
  otherIn->inputPSD       = NULL;
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
      else if (strcmp(argv[i],"-sampling")==0)
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
      else if (strcmp(argv[i], "-m1")==0)
	{
	  otherIn->m1 = atof(argv[++i]);
	}
      else if (strcmp(argv[i], "-m2")==0)
	{
	  otherIn->m2 = atof(argv[++i]);	  
	}
      else if (strcmp(argv[i], "-psi0")==0)
	{
	  otherIn->psi0 = atof(argv[++i]);	  
	}
      else if (strcmp(argv[i], "-psi3")==0)
	{
	  otherIn->psi3 = atof(argv[++i]);	  
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
	{
	  coarseIn->numFcutTemplates = atoi(argv[++i]);
	}
      else if (strcmp(argv[i],"-simType")==0)
	{
	  randIn->type = atoi(argv[++i]);
	}
      else if (strcmp(argv[i],"-sigAmp")==0)
	{
	  randIn->SignalAmp = atof(argv[++i]);
	}
      else if (strcmp(argv[i],"-alpha")==0)
	{
	  coarseIn->alpha = atof(argv[++i]); 
	  randIn->param.alpha = coarseIn->alpha;
	}
      else if (strcmp(argv[i],"--TemplateOrder")==0)
        {
          coarseIn->order = atoi(argv[++i]); 
	}
      else if (strcmp(argv[i],"--SignalOrder")==0)	
	{
	  randIn->param.order = atoi(argv[++i]); 
	}
      else if (strcmp(argv[i],"-n")==0)
	{
	  otherIn->ntrials = atoi(argv[++i]);       
	}
      else if (strcmp(argv[i],"-seed")==0)
	{
	  randIn->useed = atoi(argv[++i])+1;
	}
      else if (strcmp(argv[i],"-quiet")==0)       
	  otherIn->quietFlag = 1;
      else if (strcmp(argv[i],"--InQuadrature")==0)       
	  otherIn->overlapMethod = InQuadrature;
      else if (strcmp(argv[i],"--AmbiguityFunction")==0)       
	  otherIn->ambiguityFunction = 1;
      else if (strcmp(argv[i],"--normalisation")==0)       
	  otherIn->ambiguityFunction = 1;
      else if (strcmp(argv[i],"--FMaximization")==0)       
	  otherIn->FMaximization = 1;
      else if (strcmp(argv[i],"--PrintOverlap")==0)       
	  otherIn->PrintOverlap = 1;
      else if (strcmp(argv[i],"--PrintFilter")==0)       
	  otherIn->PrintFilter = 1;
      else if (strcmp(argv[i],"--check")==0)       
	  otherIn->check = 1;
      else if (strcmp(argv[i],"--InputPSD")==0)       
	  otherIn->inputPSD = argv[++i];


      else if (strcmp(argv[i], "--Template")==0)
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
	  if ( coarseIn->approximant == BCV ) 
		  coarseIn->space = Psi0Psi3; 
	  else 
 		  coarseIn->space = Tau0Tau3;
	}	
      else if ( (strcmp(argv[i], "--Signal")==0))
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

  if (coarseIn.approximant != BCV && otherIn.overlapMethod != InQuadrature)
    {
      fprintf(stderr,"If the template are not BCV, the overlap method must be InQuadrature (use the \"--InQuadrature\" option)\n");
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
  if(coarseIn->alpha  == -1)
    {
      coarseIn->alpha = BANKEFFICIENCY_ALPHA;
      randIn->param.alpha = coarseIn->alpha;
    }
  if(randIn->param.fendBCV  == -1)
    {
      randIn->param.fendBCV = BANKEFFICIENCY_FENDBCV;
    }
  if (randIn->type == -1)
    {
      randIn->type = BANKEFFICIENCY_TYPE;
    }
  if (randIn->SignalAmp == -1)
    {
      randIn->SignalAmp = BANKEFFICIENCY_SIGNALAMP;
    }
  if (randIn->useed == -1)
    {
      randIn->useed = BANKEFFICIENCY_SEED;
    }

  /* OtherIn structure */
  if (otherIn->quietFlag == -1)
    {
      otherIn->quietFlag = BANKEFFICIENCY_QUIETFLAG;
    }
  if (otherIn->ambiguityFunction == -1)
    {
      otherIn->ambiguityFunction = BANKEFFICIENCY_AMBIGUITYFUNCTION;
    }
  if (otherIn->ntrials == -1)
    {
      otherIn->ntrials = BANKEFFICIENCY_NTRIALS;
    }
  if (otherIn->normalisation == -1)
    {
      otherIn->normalisation = BANKEFFICIENCY_NORMALISATION;
    }
  if (otherIn->FMaximization == -1)
    {
      otherIn->FMaximization = BANKEFFICIENCY_FMAXIMIZATION;
    }
  if (otherIn->PrintOverlap == -1)
    {
      otherIn->PrintOverlap = BANKEFFICIENCY_PRINTOVERLAP;
    }
  if (otherIn->PrintFilter == -1)
    {
      otherIn->PrintFilter = BANKEFFICIENCY_PRINTFILTER;
    }
  if (otherIn->check == -1)
    {
      otherIn->check = BANKEFFICIENCY_CHECK;
    }




  
  /*no parsing for those parameters*/
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
  fprintf(stderr,"     -simType : type of simulation, 0, 1 or 2      (%d)\n\n",     BANKEFFICIENCY_TYPE);
  fprintf(stderr,"           -fl : lower frequency cutoff             (%7.2f) Hz\n", BANKEFFICIENCY_FLOWER);
  fprintf(stderr,"         -mMin : minimal mass of component stars    (%7.2f) Mo\n", BANKEFFICIENCY_MMIN);
  fprintf(stderr,"         -mMax : maximal mass of component stars    (%7.2f) Mo\n", BANKEFFICIENCY_MMAX);
  fprintf(stderr,"        -psi0Max : Max value of psi0 (%7.2f)\n",                     BANKEFFICIENCY_PSI0MAX);
  fprintf(stderr,"        -psi3Min : Min value of psi3 (%7.2f)\n",                     BANKEFFICIENCY_PSI3MIN);
  fprintf(stderr,"        -alpha : BCV amplitude correction parameter (%7.2f)\n",    BANKEFFICIENCY_ALPHA);
  fprintf(stderr,"      -numFcut : number of layers in Fcut dimension (%7.2d)\n",    BANKEFFICIENCY_NFCUT);
  fprintf(stderr,"           -mm : minimal match for template bank    (%7.3f)\n",    BANKEFFICIENCY_MMCOARSE);
  fprintf(stderr,"     --Template : PN model for template bank, e.g. BCV (%d)\n\n",   BANKEFFICIENCY_TEMPLATE);
  fprintf(stderr,"       --Signal : model to inject in Monte-Carlo, e.g. EOB (%d)\n",       BANKEFFICIENCY_SIGNAL);
  fprintf(stderr,"  --SignalOrder : order of PN model                  (%7.2d)\n",    BANKEFFICIENCY_ORDER_SIGNAL);
  fprintf(stderr,"--TemplateOrder : order of signal to injec           (%7.2d)\n",    BANKEFFICIENCY_ORDER_TEMPLATE);
  fprintf(stderr,"     -sigAmp : amplitude of the signal            (%7.2f)\n",    BANKEFFICIENCY_SIGNALAMP);

  fprintf(stderr,"     -quiet  : extra information on screen             (%d)\n",    BANKEFFICIENCY_QUIETFLAG);
  fprintf(stderr,"     -amb    : print overlap in every point of the mesh(%d)\n",    BANKEFFICIENCY_AMBIGUITYFUNCTION);

  
  fprintf(stderr,"Verbose Options \n");
  fprintf(stderr,"     --FMaximization     : Maximize over Ending frequency(%d)\n",    BANKEFFICIENCY_FMAXIMIZATION);
  fprintf(stderr,"     --PrintOverlap      : Print Overlap given by the best template (%d)\n",    BANKEFFICIENCY_PRINTOVERLAP);
  fprintf(stderr,"     --PrintFilter       : Print filters giving the best overlap template (%d)\n",    BANKEFFICIENCY_PRINTFILTER);
  fprintf(stderr,"     --AmbiguityFunction : print overlap corresponding to the meshgrid(%d)\n",    BANKEFFICIENCY_AMBIGUITYFUNCTION);
  fprintf(stderr,"     --InQuadrature : Overlap is compute without alpha maximization\n");
  fprintf(stderr,"     --check : For a given random waveform, we compute overlap with same parameter (should get overlap=1)\n");
  fprintf(stderr,"  --InputPSD: to give the name of real PSD\n");
  /*to avoid boring warning*/
  temp = coarseIn.approximant;
  temp =  randIn.useed;
  temp =  otherIn.signal;
}




void 
LALBCVAmplitude1(
		 REAL4Vector *VectorAmplitude,
		 REAL4Vector *VectorPowerFm7_6,
		 InspiralTemplate params,
		 REAL8 cte)
{
  int i, n = VectorAmplitude->length;
  REAL8 f, df    = params.tSampling / (REAL8)n;
  
  
  
  /* AFAIRE sampling pour valeur de la requence*/
  for (i=0; i< n; i++)
    {
      f = i*df;
      if (f < params.fLower || f > params.fendBCV)
	{
	  VectorAmplitude->data[i] = 0.;
	}
      else
	{
	  VectorAmplitude->data[i] = cte * VectorPowerFm7_6->data[i];            
	}
    } 
}

void 
LALBCVAmplitude2(
			  REAL4Vector *VectorAmplitude,
			  REAL4Vector VectorPowerFm7_6,
			  REAL4Vector VectorPowerFm1_2,
			  InspiralTemplate params,
			  REAL8 c1, REAL8 c2)
{
  int i, n = VectorAmplitude->length;
  REAL8 f, df    = params.tSampling/(REAL8)n;
  
  
  /* AFAIRE sampling pour valeur de la requence*/
  for (i=0; i< n; i++)
    {
      f = i*df;
      if (f < params.fLower || f > params.fendBCV)
	{
	  VectorAmplitude->data[i] = 0.;
	}
      else
	{
	  VectorAmplitude->data[i] = c1 * VectorPowerFm7_6.data[i] + c2 * VectorPowerFm1_2.data[i];
	} 
    }
}



void 
LALCreateBCVPhase(
			  REAL4Vector *VectorPhase,
			  REAL4Vector *VectorPowerFm5_3,
			  REAL4Vector *VectorPowerFm2_3,
			  InspiralTemplate params)
{

  int    i,n   = VectorPhase->length;
  REAL8 psi0  = params.psi0;
  REAL8 psi3  = params.psi3;
  
  REAL8 shift = -LAL_TWOPI * ((REAL8)params.nStartPad/params.tSampling);
  shift = 0;
  REAL8 f,df    = params.tSampling/(REAL8)n/2;


 
  for (i=0; i< n; i++)
    {
      f = i*df;
      /*fcutoff ou fendbcv ?*/
      if (f < params.fLower || f > params.fCutoff)
	{
	  VectorPhase->data[i] = 0.;
	}
      else
	{
	  VectorPhase->data[i] = (shift*f 
				  + psi0*VectorPowerFm5_3->data[i]
				  + psi3*VectorPowerFm2_3->data[i]
				  );
	  
	  
	}            
    } 
}


void LALCreateMomentVector(LALStatus             *status,
			   REAL4Vector           *a11,
			   REAL4Vector           *a21,
			   REAL4Vector           *a22,
			   REAL8FrequencySeries  *psd,
			   InspiralTemplate      *params)
{
  REAL8 
    m7=0,
    m5=0,
    m3=0,
    fMax,fMin, f;
  INT4 kMin,kMax,k;

  InspiralMomentsIn in;
  
  
  INITSTATUS (status, "LALInspiralGetBCVMaximizationMatrix", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  
  in.shf = psd;  
  in.xmin = params->fLower;
  in.xmax = params->tSampling/2.;
  in.norm =    ( psd->deltaF);

  /* the minimum and maximum frequency where we have four points */
  fMin = psd->f0 + psd->deltaF;
  fMax = psd->f0 + ((REAL8) psd->data->length - 2 ) * psd->deltaF;

  kMin = (UINT8) floor( (in.xmin - psd->f0) / psd->deltaF );
  kMax = (UINT8) floor( (in.xmax - psd->f0) / psd->deltaF );
  
  /* Skip frequency less than fLower.*/
  for (k=0; k<kMin ; k++)
    {
      a11->data[k] = 0;
      a21->data[k] = 0;
      a22->data[k] = 0;
    }
  for ( k = kMin; k < kMax; ++k )
    {
      f = psd->f0 + (REAL8) k * psd->deltaF;

      if (psd->data->data[k])
	{
	  m7 +=  pow( f, -(7./3.) ) / psd->data->data[k] * psd->deltaF / in.norm;
	  m5 +=  pow( f, -(5./3) )  / psd->data->data[k] * psd->deltaF / in.norm;
	  m3 +=  pow( f, -(1.) )    / psd->data->data[k] * psd->deltaF / in.norm;


	  a11->data[k] = 1./sqrt(m7);  
	  a22->data[k] = 1. / sqrt(m3 - m5*m5/m7);
	  a21->data[k] = -m5 / m7 * a22->data[k];
	}
    }


  DETATCHSTATUSPTR(status);
  RETURN (status);
}

  

  


void
LALMatrix(
	  LALStatus               *status,
	  BCVMaximizationMatrix   *output,
	  REAL8FrequencySeries    *psd,
	  InspiralTemplate        params 
	  )
     /*  </lalVerbatim>  */
{
  REAL8 m3, m5, m7;
  InspiralMomentsIn in;
  
  
  INITSTATUS (status, "LALInspiralGetBCVMaximizationMatrix", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  
  /* 
     First compute the moments needed for the maximization 
     We need only I7, I5 and I3 
  */

  in.shf = psd;  
  in.xmin = params.fLower;
  in.xmax = params.fCutoff;


  in.norm =    ( psd->deltaF);
  

  in.ndx = 7. / 3. ;
  LALMoments(status->statusPtr, &m7, in); 
  CHECKSTATUSPTR(status);
  


  in.ndx = 5. / 3.;
  LALMoments(status->statusPtr, &m5, in); 
  CHECKSTATUSPTR(status);


    in.ndx = 3. / 3.;
  LALMoments(status->statusPtr, &m3, in); 
  CHECKSTATUSPTR(status);
 
  /* Then we compute the values of the matrix */



  
  output->a11 = 1./ sqrt(m7);
  output->a22 = 1. / sqrt(m3 - m5*m5/m7);
  output->a21 = -m5 / m7 * output->a22;


  
  DETATCHSTATUSPTR(status);
  RETURN(status);
  
}

















void
LALCreateBCVTemplate(
		     REAL4Vector *signal,
		     REAL4Vector VectorAmplitude,
		     REAL4Vector VectorPhase,
		     InspiralTemplate params)
{
  int i, n = signal->length;
  REAL8 f, df = params.tSampling / (REAL8)(VectorPhase.length)/2.;
  
  signal->data[0] = 0.0;
  signal->data[n/2] = 0.0;
  
  for(i=1;i<n/2;i++)
    {
      f = i*df;
      if (f < params.fLower || f > params.fendBCV)
	{
	  /* 
	   * All frequency components below params->fLower and above fn are set to zero  
	   */
	  signal->data[i]   = 0.;
	  signal->data[n-i] = 0.;
	}
      else
	{
	  signal->data[i]   = (REAL4) ( VectorAmplitude.data[i] * cos(VectorPhase.data[i]));	  
	  signal->data[n-i] = (REAL4) (-VectorAmplitude.data[i] * sin(VectorPhase.data[i]));
	}
    }
  
}


void LALGetOrthoNormalFilter(REAL4Vector *filter2, REAL4Vector filter1)
{
	int i,n,nby2;

	n = filter1.length;
	nby2 = n/2;
	filter2->data[0] = filter1.data[0];
	filter2->data[nby2] = filter1.data[nby2];
	for (i=1; i<nby2; i++)
	{
		filter2->data[i] = -filter1.data[n-i];
		filter2->data[n-i] = filter1.data[i];
	}
}

  

/*  <lalVerbatim file="LALInspiralWaveOverlapBCVCP"> */
void
LALWaveOverlapBCV(
			     LALStatus *status,
			     REAL4Vector             *correlation,
			     InspiralWaveOverlapOut  *overlapout,
			     InspiralWaveOverlapIn   *overlapin,
			     REAL4Vector              VectorPhase, 
			     REAL4Vector              VectorAmplitude1,
			     REAL4Vector              VectorAmplitude2,
			     BCVMaximizationMatrix    matrix,
			     OtherParamIn             otherIn 
			     )
     /*  </lalVerbatim>  */
{
  
  REAL4Vector x1,x2,x3,x4, filter1, filter2, filterextra;
  InspiralWaveCorrelateIn corrin;
  UINT4 i, nBegin, nEnd, n = correlation->length; 
  REAL8 phi, phase, rho;
  REAL8 x1_2,x2_2,x3_2,x4_2, V0,V1, V2, a,b,c;
   REAL8 ax1,ax2,ax3,ax4;

  x1.length = x2.length = x3.length = x4.length = n;
  filter1.length = filter2.length = filterextra.length = n;
  
  x1.data = (REAL4*) LALMalloc(sizeof(REAL4) * x1.length);
  x2.data = (REAL4*) LALMalloc(sizeof(REAL4) * x2.length);
  x3.data = (REAL4*) LALMalloc(sizeof(REAL4) * x3.length);
  x4.data = (REAL4*) LALMalloc(sizeof(REAL4) * x4.length);
  
  filter1.data = (REAL4*) LALMalloc(sizeof(REAL4) * filter1.length);
  filter2.data = (REAL4*) LALMalloc(sizeof(REAL4) * filter2.length);
  filterextra.data = (REAL4*) LALMalloc(sizeof(REAL4) * filterextra.length);
  
  INITSTATUS (status, "LALInspiralWaveOverlapForBCV", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  
  
  overlapin->param.nStartPad = 0;
  overlapin->param.startPhase = 0;

  /* to compute the scalar products xi,= <x,hk> we need TF(x)==========================*/
  corrin.fCutoff      = overlapin->param.fFinal;
  corrin.samplingRate = overlapin->param.tSampling;
  corrin.df           = overlapin->param.tSampling / (REAL8) overlapin->signal.length;
  corrin.psd          = overlapin->psd;
  corrin.signal1      = overlapin->signal;
  corrin.revp         = overlapin->revp;
  
  /*Filter h1 and h3 */
  


  LALCreateBCVTemplate(&filter1, VectorAmplitude1, VectorPhase,overlapin->param);  

  if (otherIn.PrintFilter == 1)
    {
      for (i=0; i< filter1.length; i++)
	printf("%d %e\n", i,filter1.data[i]);
      printf("&\n");	       
    }


  corrin.signal2      = filter1;
  LALInspiralWaveCorrelate(status->statusPtr, &x1, corrin);
  CHECKSTATUSPTR(status);
  LALGetOrthoNormalFilter( &filter2,filter1);
  if (otherIn.PrintFilter == 1)
    {
      for (i=0; i< filter1.length; i++)
	printf("%d %e\n", i,filter2.data[i]);
      printf("&\n");	       
    }

  corrin.signal2      = filter2;
  LALInspiralWaveCorrelate(status->statusPtr, &x3, corrin);
  CHECKSTATUSPTR(status);
  
  LALCreateBCVTemplate(&filter1, VectorAmplitude2, VectorPhase, overlapin->param);  
  if (otherIn.PrintFilter == 1)
    {
      for (i=0; i< filter1.length; i++)
	printf("%d %e\n", i,filter1.data[i]);
      printf("&\n");	       
    }

  corrin.signal2      = filter1;
  LALInspiralWaveCorrelate(status->statusPtr, &x2, corrin);
  CHECKSTATUSPTR(status);
  LALGetOrthoNormalFilter( &filter2,filter1);
  if (otherIn.PrintFilter == 1)
    {
      for (i=0; i< filter1.length; i++)
	printf("%d %e\n", i,filter2.data[i]);
      printf("&\n");	       
    }

  corrin.signal2      = filter2;
  LALInspiralWaveCorrelate(status->statusPtr, &x4, corrin);
  
  
  /*nbegin et nEnd ici*/
  
  nBegin            = overlapin->nBegin;
  nEnd              = filter1.length - overlapin->nEnd;
  overlapout->max   = 0.;
  overlapout->bin   = nBegin;
  overlapout->phase = 0;

  if (otherIn.PrintOverlap == 1)
    {
      for (i=0; i< x1.length; i++)
	printf("%d %e\n", i,fabs(x1.data[i]));
      printf("&\n");	       
      
      for (i=1; i< filter1.length; i++)
	printf("%d %e\n", i,fabs(x2.data[i]));
      printf("&\n");	       
      
      for (i=0; i< filter1.length; i++)
	printf("%d %e\n", i,fabs(x3.data[i]));
      printf("&\n");	       
      
      for (i=0; i< filter1.length-1; i++)
	printf("%d %e\n", i,fabs(x4.data[i]));
      printf("&\n");	       
    }



  for (i=nBegin; i<nEnd; i++) 
    {
      ax1 = x1.data[i];
      ax2 = x2.data[i] ;
      ax3 = x3.data[i];
      ax4 = x4.data[i] ;

      x1_2 = x1.data[i] * x1.data[i];
      x2_2 = x2.data[i] * x2.data[i];
      x3_2 = x3.data[i] * x3.data[i];
      x4_2 = x4.data[i] * x4.data[i];      
      
      a = x1_2 + x3_2;
      b = x2_2 + x4_2;
      c = a + b;

      V0 = x1_2 + x2_2 + x3_2 + x4_2;
      V1 = x1_2 + x3_2 - x2_2 - x4_2;
      V2 = 2*ax1*ax2 + 2*ax3*ax4;

      /* 1/2 * sqrt(V0 + sqrt(V1^2 +V2^2)) */
      rho = 0.5 * sqrt (c + 
			sqrt(
			     a*a + b*b -2*a*b 
			     + 4. * (x1_2*x2_2 
				     + 2*x1.data[i]*x2.data[i]*x3.data[i]*x4.data[i]
				     + x3_2*x4_2)
			     )
			);
      rho =  sqrt(V0 +sqrt(V1*V1+V2*V2))/sqrt(2);


      if (otherIn.PrintOverlap == 1)
	{         
	  printf("%d %e %e %e %e\n", i, rho,V0,V1,V2);
	}
      
      if (rho > overlapout->max)
	{
	  overlapout->max = rho;
	  overlapout->bin = i;
	  overlapout->phase =  .5*atan2(V2,V1) ;
	  overlapout->phase =  .5*atan2(V2,V1) ;
	}          
    }


  while (overlapout->phase > LAL_PI || overlapout->phase > 0 ) overlapout->phase -= LAL_PI;
  while (overlapout->phase < -LAL_PI || overlapout->phase < 0 ) overlapout->phase += LAL_PI;


  /* Check whether phase is in the appropriate range */


  if (otherIn.PrintOverlap == 1)
    printf("&\n");
  
  phase = overlapout->phase;


    phi   = overlapin->param.startPhase;
    for (i=0; i<(UINT4)correlation->length; i++) 
      {
	correlation->data[i] = x1.data[i] * cos(phase) *cos(phi)
	  + x2.data[i] * sin(phase) * cos(phi)
	  + x3.data[i] * cos(phase) * sin(phi)
	  + x4.data[i] * sin(phase) * sin(phi);	
      }
    

    overlapout->alpha = -(matrix.a22 * tan(phase)) / (matrix.a11 + matrix.a21* tan(phase));
    
    
    /*    printf("#%e %e %e\n",overlapout->alpha, overlapout->max, phase);    */

    LALFree(filter1.data);
    LALFree(filter2.data);
    LALFree(filterextra.data);

    LALFree(x1.data);
    LALFree(x2.data);
    LALFree(x3.data);
    LALFree(x4.data);

    DETATCHSTATUSPTR(status);
    RETURN(status);
  
  
}



void LALMoments(LALStatus         *status,
                        REAL8             *moment,
                        InspiralMomentsIn pars)
{ /* </lalVerbatim> */

  REAL8 f;
  REAL8 momentTmp;
  REAL8 fMin;
  REAL8 fMax;
  REAL8 deltaF;
  UINT4 k;
  UINT4 kMin;
  UINT4 kMax;

  INITSTATUS (status, "LALInspiralMoments", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);

  ASSERT (pars.shf, status, 
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT (pars.shf->data, status, 
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);
  ASSERT (pars.shf->data->data, status, 
      LALINSPIRALBANKH_ENULL, LALINSPIRALBANKH_MSGENULL);

  /* make sure that the minimum and maximum of the integral are within */
  /* the frequency series                                              */
  fMax = pars.shf->f0 + (REAL8) pars.shf->data->length * 
    pars.shf->deltaF;
  if ( pars.xmin < pars.shf->f0 || pars.xmax > fMax )
  {
    ABORT( status, 999, "limits outside range of frequency series" );
  }

  /* the minimum and maximum frequency where we have four points */
  deltaF = pars.shf->deltaF;
  fMin = pars.shf->f0 + deltaF;
  fMax = pars.shf->f0 + ((REAL8) pars.shf->data->length - 2 ) * deltaF;

  if ( pars.xmin <= fMin )
  {
    kMin = 1;
  }
  else
  {
    kMin = (UINT8) floor( (pars.xmin - pars.shf->f0) / deltaF );
  }

  if ( pars.xmax >= fMax )
  {
    kMax = pars.shf->data->length - 1;
  }
  else
  {
    kMax = (UINT8) floor( (pars.xmax - pars.shf->f0) / deltaF );
  }

  momentTmp = 0.;
   *moment =0; 

  for ( k = kMin; k <= kMax; ++k )
  {
    momentTmp = 0.;
    f = pars.shf->f0 + (REAL8) k * deltaF;
    if (pars.shf->data->data[k]) momentTmp = pow( f, -(pars.ndx) ) / pars.shf->data->data[k];

     (*moment) += momentTmp;
  }

  *moment *= deltaF;

  /* now divide the moment by the specified norm */
  *moment /= pars.norm;

  DETATCHSTATUSPTR(status);
  RETURN (status);
}





void 
LALCreateVectorFreqPower(
				 REAL4Vector *vector,
				 InspiralTemplate params,
				 int a,
				 int b)
{
  int i , n = vector->length;
  REAL8 power = (REAL8)a / (REAL8)b ;
  REAL8 f, df = params.tSampling/(REAL8)n/2.;;


  vector->data[0] = 0.;
  /* inclure sampling pour la valuer de f*/  
  for (i=1; i<n; i++)
    {
      f = i *df; 
      vector->data[i] = pow(f, power);
    }
}

void LALGenerateWaveform(LALStatus              *status,
			 REAL4Vector            *signal,
			 RandomInspiralSignalIn *randIn)


{
  REAL8  norm;
  REAL4Vector  buff;
  InspiralWaveNormaliseIn normin;
  
  INITSTATUS (status, "LALGenerateWaveform", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  
  buff.length = signal->length;
  if (!(buff.data = (REAL4*) LALMalloc(sizeof(REAL4)*buff.length))) {
    ABORT (status, LALNOISEMODELSH_EMEM, LALNOISEMODELSH_MSGEMEM);
  }
  
  /* set up the structure for normalising the signal */
  normin.psd = &(randIn->psd);
  normin.df = randIn->param.tSampling / (REAL8) signal->length;
  normin.fCutoff = randIn->param.fCutoff;
  normin.samplingRate = randIn->param.tSampling;


  LALInspiralParameterCalc(status->statusPtr, &(randIn->param));
  
  if (randIn->param.approximant == TaylorF1 ||
      randIn->param.approximant == TaylorF2 ||
      randIn->param.approximant == PadeF1)
    {
      LALInspiralWave(status->statusPtr, signal, &randIn->param);
      CHECKSTATUSPTR(status);
    }
  else
    {
      randIn->param.fFinal=0;
      LALInspiralWave(status->statusPtr, &buff, &randIn->param);
      CHECKSTATUSPTR(status);
      LALREAL4VectorFFT(status->statusPtr, signal, &buff, randIn->fwdp);
      CHECKSTATUSPTR(status);
    }
  LALInspiralWaveNormaliseLSO(status->statusPtr, signal, &norm, &normin);
  CHECKSTATUSPTR(status);
  
  
  if (buff.data != NULL) LALFree(buff.data);
  DETATCHSTATUSPTR(status);
   RETURN(status);
}








/*

void 
LALReadSpectralDensity 
(
 LALStatus    *status, 
 REAL8Vector  *psd,  
 char         *name,
 REAL8        df
   ) 
{ 
  INT4 n,i;  
  REAL4 freq;
  FILE *inputFile; 
  
  INITSTATUS(status, "LALNoiseSpectralDensity", BANKEFFICIENCYC);
  ATTATCHSTATUSPTR(status);
  
  n = psd->length;
  inputFile = fopen( name, "r");  
  
  for (i=0; i<n; i++) 
    {
      fscanf(inputFile, "%e %le", &freq, &psd->data[i]);	  
    } 
  
  fclose(inputFile);
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
*/
