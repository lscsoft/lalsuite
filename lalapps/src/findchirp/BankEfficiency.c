 
/******************************** <lalVerbatim file="BankEfficiencyCV">
Author: Cokelaer, T. and Sathyaprakash, B. S.
< $Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{BankEfficiency.c}}
\label{ss:BankEfficiency.c}

\subsubsection*{Description}


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

/* Some Error messages */
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

/* Some constantes to fill InspiralTemplate Structure, Bank structure and internal structure */
/* Bank structure first: */
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
#define BANKEFFICIENCY_LOWGM                    3
#define BANKEFFICIENCY_HIGHGM                    6
/* Other Parameters */
#define BANKEFFICIENCY_QUIETFLAG       	        0
#define BANKEFFICIENCY_AMBIGUITYFUNCTION      	0
#define BANKEFFICIENCY_NORMALISATION      	0
#define BANKEFFICIENCY_FMAXIMIZATION      	0
#define BANKEFFICIENCY_PRINTOVERLAP             0
#define BANKEFFICIENCY_PRINTFILTER              0
#define BANKEFFICIENCY_CHECK                    0




typedef enum{
  LIGOI,
    LIGOII,
    GEO,
    TAMA,
    VIRGO
}DetectorName;


/* Structure to store the value of a 2 by 2 matrix. 
 * This matrix is used in the alpha maximization 
 * process for BCV templates  */
typedef struct{
  double a11,a21,a22;
}BCVMaximizationMatrix;


/* Choice of the overlap method. InQuadrature maximize over 
 * the phase parameter; AlphaMaximization maximize over both 
 * phase and alpha parameter (BCV templates) */
typedef enum {
  InQuadrature,
  AlphaMaximization
} OverlapMethodIn;


/* Internal parameters for the Bank Efficiency code:
 * signal : name of the random signal to inject
 * template: name of the template to use
 * bank: type of bank to use
 * ntrials: number of signal to inject for Monte Carlo Simulation
 * quietFlag: some extra output 
 * ambiguityFunction: print the AF (ntrial is then forced to be equal to 1
 * normalization: useless ? 
 * FMaximization: Once the bank give a detection, with psi0 and psi3 optimal, 
 * 		  we maximize over the last frequency
 * PrintOverlap:
 * PrintFilter:
 * overlapMethod: InQuadrature(classic overlap)  or AlphaMaximization
 * check: force the unique template to have the same psi0, psi3 and FFinal as the injected signal
 * m1: mass1 to inject
 * m2: mass2 to inject
 * psi0: psi0 to inject
 * psi3: psi3 to inject
 * inputPSD: name of an input file for the psd.*/
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
  DetectorName NoiseModel;

} OtherParamIn;


/* As to be cleaned
 * Function to store the optimal  parameters of the overlap 
 * lmax : number of the layer
 * fMax : frequency cutoff of the template
 * jmax: number of the best template. use if one want to use the FMaximization option */
void
KeepHighestValues(InspiralWaveOverlapOut in , 
		  INT4 j,INT4 l , double frequency,
		  InspiralWaveOverlapOut *out, 
		  INT4 *jmax, INT4 *lmax, double *fMax);


/* function to create the filters in the BCV overlap.
 * Input: 	- Vectors F^(-5/3) and so on
 * 		- Matrix of the moments for the amplitude
 * 		- indice at which the template starts
 * 		- psi0 and psi3 for the phase
 * Output	- The two orthogonal filters */
void 
LALCreateFilters(REAL4Vector *Filter1,
		 REAL4Vector *Filter2,
		 REAL4Vector VectorPowerFm5_3,
		 REAL4Vector VectorPowerFm2_3,
		 REAL4Vector VectorPowerFm7_6,
		 REAL4Vector VectorPowerFm1_2,
		 BCVMaximizationMatrix    matrix,
		 INT4 kMin,
		 double psi0,
		 double psi3
		 );


/* Function to compute a orthogonal vector */
void
LALGetOrthoNormalFilter2(REAL4Vector *filter);


/* Functon to compute the overlap between an injected signal 
 * and a set of templates based on BCV templates. 
 * It returns the overlap, the phase and alpha parameter. */
void
LALWaveOverlapBCV(LALStatus *status,
		  REAL4Vector             *correlation,
		  InspiralWaveOverlapOut  *overlapout,
		  InspiralWaveOverlapIn   *overlapin,
		  REAL4Vector             *Filter1, 
		  REAL4Vector             *Filter2,
		  BCVMaximizationMatrix    matrix,
		  OtherParamIn             otherIn 
		  );

/* Function to store the moments needed by the BCV overlap process 
 * a11, a22, a21 are three vectors which stored the three components 
 * of the matrix needed in the BCV maximization process. */
void
LALCreateMomentVector(LALStatus             *status,
		      REAL4Vector           *a11,
		      REAL4Vector           *a21,
		      REAL4Vector           *a22,
		      REAL8FrequencySeries  *psd,
		      InspiralTemplate      *params);

/* Function to create Vectors of the form y(f)= f^(-a/b). 
 * where f is the frequency. */
void 
LALCreateVectorFreqPower( REAL4Vector *vector,
			  InspiralTemplate params,
			  int a,
			  int b);


/* function to generate waveform (has to be tested). */
void
LALGenerateWaveform(LALStatus *status,
		    REAL4Vector *signal,
		    RandomInspiralSignalIn *params);

/* print a time series */
void
printf_timeseries (INT4 n,
		   REAL4 *signal,
		   REAL8 delta,
		   REAL8 t0);


/* Init the CoarseBank structure to dummy values */
void
Init2DummyCoarse(
		 InspiralCoarseBankIn *coarseIn);


/* Init the Random strucutre to dummy values */
void
Init2DummyRandom(
		 RandomInspiralSignalIn *randIn);

/* Init Other Parameters to dummy values */
void
Init2DummyOtherParam(
		     OtherParamIn *otherIn);

/* Parsing function */
void
ParseParameters(int 			*argc, 
		char 			**argv,
		InspiralCoarseBankIn 	*coarseIn,
		RandomInspiralSignalIn 	*randIn,
		OtherParamIn 		*otherIn);		     

/* function to check validity of parameters*/
void 	
CheckParams(InspiralCoarseBankIn 	coarseIn,
	    RandomInspiralSignalIn 	randIn,
	    OtherParamIn 		otherIn);


/* Default values*/
void 
SetDefault(InspiralCoarseBankIn 	*coarseIn, 
	   RandomInspiralSignalIn 	*randIn,
	   OtherParamIn 		*otherIn);

/* Help Function */
void 	
Help(     InspiralCoarseBankIn 	coarseIn,
	  RandomInspiralSignalIn 	randIn,
	  OtherParamIn 		otherIn);

/* Print Function for final results */
void
PrintResults(
	    InspiralTemplate       bank,
	    InspiralTemplate       injected,
	    InspiralWaveOverlapOut overlapout, double a, int b);

/* Print the bank */
void
PrintBank(
	  InspiralCoarseBankIn coarse, 
	  InspiralTemplateList **list,
	  UINT4 nlist,
	  UINT4 signalLength);


extern int lalDebugLevel=1;


/* ===================================== THE MAIN PROGRAM ===================================================== */
int
main (int argc, char **argv ) 
{
  /* --- Variables ---*/
  INT4
    nn,
    nby2,
    l,
    lmax,
    nbank=0,
    i,
    j=0,
    n,
    k,
    jmax,
    nlist,
    kMin;
  
  REAL8
    df,
    fendBCV,
    fMax;
  
  REAL4Vector 
    signal,
    correlation,
    VectorAmplitude1,
    VectorAmplitude2,
    VectorPhase,
    VectorPowerFm5_3,
    VectorPowerFm2_3,
    VectorPowerFm7_6,
    VectorPowerFm1_2,
    VectorA11,
    VectorA21,
    VectorA22,
    Filter1, 
    Filter2;
  
  
  void   			*noisemodel;
  
  RandomInspiralSignalIn
    randIn;			/* random signal waveform to inject*/
  InspiralWaveOverlapIn 	
    overlapin;			/* structure for Overlap*/
  
  InspiralWaveOverlapOut 
    overlapout, 
    overlapoutmax;
  
  OtherParamIn 		 
    otherIn;			/* personal structure to parse some independant parameters*/
  
  RealFFTPlan 
    *fwdp=NULL,
    *revp=NULL;
  
  InspiralTemplateList  
    *list=NULL;
   
  InspiralCoarseBankIn    
    coarseIn; 			/* strcture for the bank of templates*/

  static LALStatus status;
  
  BCVMaximizationMatrix
    matrix;

   
  
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
  randIn.param.OmegaS  = 0.;
  randIn.param.Theta   = 0.;
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
  switch (otherIn.NoiseModel)
    {
    case LIGOI:  noisemodel = LALLIGOIPsd ; break;
    case VIRGO: noisemodel = LALVIRGOPsd;  break;
    case GEO:   noisemodel = LALGEOPsd;    break;
    case TAMA:  noisemodel = LALTAMAPsd;   break;
    default:    noisemodel = LALLIGOIPsd; 
    }
  LALNoiseSpectralDensity (&status, coarseIn.shf.data, noisemodel, df);
  
  
  for (i=0; i< (int)randIn.psd.length; i++)
    {
      randIn.psd.data[i] = coarseIn.shf.data->data[i];
    }
  
  
  
  /* --- And the bank of templates --- */
  l =  coarseIn.numFcutTemplates;
  /*   coarseIn.numFcutTemplates = 1;*/
  LALInspiralCreateCoarseBank(&status, &list, &nlist, coarseIn);
  coarseIn.numFcutTemplates = l;
  if (nlist==0) exit(0);	
  



  
  /* --- Let's print some optional comments : the position of templates in the bank ---*/
  if (!otherIn.quietFlag) 
    {
      /*PrintBank(coarseIn, &list, nlist, signal.length);*/
      for (i=0; i<nlist; i++)
	{
	  if (coarseIn.approximant == BCV)
	    {      		   
	      fprintf(stdout, "%e %e %e %e\n", (list[i]).params.psi0, (list[i]).params.psi3, (list[i]).params.totalMass, (list[i]).params.fFinal);
	    } 
	  else
	    {
	      fprintf(stdout, "%e %e %e %e\n", (list[i]).params.mass1, (list[i]).params.mass2, (list[i]).params.t0, (list[i]).params.t3);
	    }
	}	
    
      fprintf(stdout, "#---------------------------------------------------------------\n");
      fprintf(stdout, "#Number of Coarse Bank Templates=%d\n",nlist);
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
  
  
  if (otherIn.ambiguityFunction  == 1)      otherIn.ntrials=1;



   list[1].params.startPhase =0;   
  

   
   /* Then we can  compute the frequency vector*/
   VectorPhase.length      = signal.length/2;
   VectorAmplitude1.length = signal.length/2;
   VectorAmplitude2.length = signal.length/2;
   VectorA11.length        = signal.length/2 ;
   VectorA22.length        = signal.length/2;
   VectorA21.length        = signal.length/2 ;
   VectorPowerFm2_3.length = signal.length/2 ;
   VectorPowerFm5_3.length = signal.length/2;       
   VectorPowerFm7_6.length = signal.length/2 ;
   VectorPowerFm1_2.length = signal.length/2;       
   
   Filter1.length          = correlation.length;
   Filter2.length          = correlation.length;
   
   VectorPhase.data        = (REAL4*) LALMalloc(sizeof(REAL4) * VectorPhase.length);   
   VectorAmplitude1.data   = (REAL4*) LALMalloc(sizeof(REAL4) * VectorAmplitude1.length);
   VectorAmplitude2.data   = (REAL4*) LALMalloc(sizeof(REAL4) * VectorAmplitude2.length);
   VectorA11.data          = (REAL4*) LALMalloc(sizeof(REAL4) * VectorA11.length);
   VectorA21.data          = (REAL4*) LALMalloc(sizeof(REAL4) * VectorA21.length);
   VectorA22.data          = (REAL4*) LALMalloc(sizeof(REAL4) * VectorA22.length);
   VectorPowerFm2_3.data   = (REAL4*) LALMalloc(sizeof(REAL4) * VectorPowerFm2_3.length);
   VectorPowerFm5_3.data   = (REAL4*) LALMalloc(sizeof(REAL4) * VectorPowerFm5_3.length);
   VectorPowerFm7_6.data   = (REAL4*) LALMalloc(sizeof(REAL4) * VectorPowerFm7_6.length);
   VectorPowerFm1_2.data   = (REAL4*) LALMalloc(sizeof(REAL4) * VectorPowerFm1_2.length);

   
   Filter1.data            = (REAL4*) LALMalloc(sizeof(REAL4) * Filter1.length);
   Filter2.data            = (REAL4*) LALMalloc(sizeof(REAL4) * Filter2.length);

   for (i=0; i<(INT4)Filter1.length; i++) Filter1.data[i] = 0.;	   	   	  
   for (i=0; i<(INT4)Filter2.length; i++) Filter2.data[i] = 0.;	   	   	  

   n = VectorAmplitude1.length;
   df = randIn.param.tSampling /(REAL8)n/2.;      



   /* Create useful vector once for all */
   LALCreateVectorFreqPower(&VectorPowerFm5_3, randIn.param, -5, 3);
   LALCreateVectorFreqPower(&VectorPowerFm2_3, randIn.param, -2, 3);
   LALCreateVectorFreqPower(&VectorPowerFm7_6, randIn.param, -7, 6);
   LALCreateVectorFreqPower(&VectorPowerFm1_2, randIn.param, -1, 2);
   /* Create Matrix for BCV Maximization once for all c*/
   LALCreateMomentVector(&status, &VectorA11, &VectorA21, &VectorA22, &coarseIn.shf, &(list[j].params));
  
  
   /* If check option is on, we compute the best overlap (we should get one if template and approximant are the same) */
   if (otherIn.check==1)
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
       jmax           = 0;
       overlapout.max = 0.0;
       
       nby2 = Filter1.length/2;
       nn = Filter1.length;
       kMin  = floor(randIn.param.fLower /df );
       df    = randIn.param.tSampling / (REAL8)(VectorAmplitude1.length)/2.;       
       
       overlapin.param 	= randIn.param;
       list[1].params       = randIn.param;
       overlapin.param.fFinal  = list[1].params.fFinal;
       overlapin.param.fCutoff = list[1].params.fFinal;
       
       k = floor(overlapin.param.fFinal / df);
       matrix.a21 = VectorA21.data[k];
       matrix.a11 = VectorA11.data[k];
       matrix.a22 = VectorA22.data[k];
       
       LALCreateFilters(&Filter1,
			&Filter2,
			VectorPowerFm5_3,
			VectorPowerFm2_3,
			VectorPowerFm7_6,
			VectorPowerFm1_2,
			matrix,
			kMin,					       
			list[1].params.psi0,
			list[1].params.psi3
			);
       
       
       /* The overlap given two filters and the input signal*/
       LALWaveOverlapBCV(&status, 
			 &correlation,
			 &overlapout,
			 &overlapin,
			 &Filter1,
			 &Filter2,
			 matrix,
			 otherIn);
       PrintResults(list[1].params, randIn.param,overlapout, list[1].params.fFinal, 1);
       j = nlist +1;
     }
   

   else /* the real code is here*/
     
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
	 jmax           = 0;
	 overlapoutmax.max = 0.0;
	 
	 
	 /* we might force to use only one template giving the input psi0 and psi3. 
	    Then, we don't need to use the bank*/
	 if (otherIn.psi0 != -1 && otherIn.psi3 != -1)
	   {
	     list[0].params.psi0 = otherIn.psi0;
	     list[0].params.psi3 = otherIn.psi3;
	     nlist  = 1;
	   }
	 
	 
	 nby2 = Filter1.length/2;
	 nn = Filter1.length;
	 kMin  = floor(randIn.param.fLower /df );
	 df    = randIn.param.tSampling / (REAL8)(VectorAmplitude1.length)/2.;
	 
	 
	 nbank = 0;
	 /* otherwise, we use the bank*/                         
	 for (j=0; j<nlist; j++) 
	   {	   	     	     	     
	     overlapin.param 	= list[j].params;	   
	     
	     
	     switch(otherIn.overlapMethod)
	       {
	       case AlphaMaximization:	   	   	   
		 
		 /*  fend ? */
		 
		 /*	         for (l = 0; l< (INT4)coarseIn.numFcutTemplates; l++)
				 {
				 if (coarseIn.numFcutTemplates == 1) 
				 frac = 1;
				 else
				 frac = (1.L - 1.L/pow(2.L, 1.5L)) / (coarseIn.numFcutTemplates-1.L);
				 
				 
				 fendBCV  = list[j].params.fFinal * 
				 (1.L - (REAL8) l * frac);
				 
				 
				 
				 if (fendBCV  > randIn.param.fLower &&
				 fendBCV < list[j].params.tSampling / 2.0)
				 {			   
		 */    nbank++;
	       fendBCV = list[j].params.fFinal;
	       overlapin.param.fFinal   =  fendBCV;
	       overlapin.param.fCutoff  =  fendBCV;
	       
	       
	       /* Extract some values for template creation*/
	       k = floor(fendBCV / df);
	       matrix.a21 = VectorA21.data[k];
	       matrix.a11 = VectorA11.data[k];
	       matrix.a22 = VectorA22.data[k];
	       
	       /* Template creation*/
	       LALCreateFilters(&Filter1,
				&Filter2,
				VectorPowerFm5_3,
				VectorPowerFm2_3,
				VectorPowerFm7_6,
				VectorPowerFm1_2,
				matrix,
				kMin,					       
				list[j].params.psi0,
				list[j].params.psi3
				);
	       
	       /* The overlap given two filters and the input signal*/
	       
	       
	       LALWaveOverlapBCV(&status, 
				 &correlation,
				 &overlapout,
				 &overlapin,
				 &Filter1,
				 &Filter2,
				 matrix,
				 otherIn);
	       
	       KeepHighestValues(overlapout, j, l, fendBCV, 
				 &overlapoutmax, &jmax, &lmax, &fMax);
	       
	       /*		       		     }
						     }
	       */    	       
	       break;
	       case InQuadrature:
		 if (coarseIn.approximant ==BCV)
		   {
		     /* for (l = 0; l< (INT4)coarseIn.numFcutTemplates; l++)
			{
			if (coarseIn.numFcutTemplates == 1) 
			frac = 1;
			else
			frac = (1.L - 1.L/pow(2.L, 1.5L)) / (coarseIn.numFcutTemplates-1.L);
			
			fendBCV = list[j].params.fFinal * (1.L - (REAL8) l * frac);
			
			if (fendBCV > randIn.param.fLower &&
			fendBCV < list[j].params.tSampling / 2.0)
			{			   
		     */
		     fendBCV  =  list[j].params.fFinal;   
		     overlapin.param.fFinal = fendBCV;
		     for (i=0; i<(INT4)signal.length; i++) correlation.data[i] = 0.;	   	   	  
		     
		     LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
		     
		     KeepHighestValues(overlapout, j, l, fendBCV, 
				       &overlapoutmax, &jmax, &lmax, &fMax);
		     
		     
		     
		   }
		 else
		   {
		     /* SHOULD be replace by flso of the template in the bank*/
		     fendBCV   = 1./LAL_PI/pow(6, 1.5)/(list[j].params.mass1+list[j].params.mass2)/LAL_MTSUN_SI;
		     
		     overlapin.param.fFinal  = fendBCV;
		     overlapin.param.fCutoff = fendBCV;
		     
		     for (i=0; i<(INT4)signal.length; i++) correlation.data[i] = 0.;	   	   	  
		     
		     
		     LALInspiralWaveOverlap(&status,&correlation,&overlapout,&overlapin);
		     KeepHighestValues(overlapout, j, l, fendBCV, &overlapoutmax,  &jmax, &lmax, &fMax);
		     
		   }
		 
		 
		 overlapout.alpha = -1;	       
		 break;	     
	       }/* switch end */
	     
	     if (otherIn.ambiguityFunction  == 1 && overlapout.max!=0)
	       printf("%lf %lf %lf %lf %d %lf\n", 
		      list[j].params.psi0,
		      list[j].params.psi3, 
		      overlapout.max,
		      fendBCV, list[j].nLayer, randIn.param.fFinal);
	     
	   }/*end of the bank process*/
	 
	 
	 /*Now print some results or other stuffs*/
	 
	 
	 
	 LALInspiralParameterCalc(&status, &list[jmax].params);
	 PrintResults(list[jmax].params, randIn.param, overlapoutmax, fMax, list[jmax].nLayer);
	 /* we might want to maximize the optimal psi0-psi3 over a frequency range*/
	 if (coarseIn.approximant==BCV &&  otherIn.FMaximization)
	   {
	     for (fendBCV = list[jmax].params.fLower ;
		  fendBCV < 800;
		  fendBCV+=10)
	       {	   
		 /* we need fcutoff in LALMatrix*/
		 
		 
		 k = floor(fendBCV / df);
		 matrix.a21 = VectorA21.data[k];
		 matrix.a11 = VectorA11.data[k];
		 matrix.a22 = VectorA22.data[k];
		 
		 /* Template creation*/
		 LALCreateFilters(&Filter1,
				  &Filter2,
				  VectorPowerFm5_3,
				  VectorPowerFm2_3,
				  VectorPowerFm7_6,
				  VectorPowerFm1_2,
				  matrix,
				  kMin,					       
				  list[jmax].params.psi0,
				  list[jmax].params.psi3
				  );
		 
		 
		 
		 overlapin.param.fFinal = fendBCV;
		 overlapin.param.fCutoff = fendBCV;
		 
		 otherIn.PrintOverlap = 0;
		 otherIn.PrintFilter = 0;
		 
		 LALWaveOverlapBCV(&status,
				   &correlation,
				   &overlapout,
				   &overlapin,
				   &Filter1,
				   &Filter2,
				   matrix,
				   otherIn
				   );
		 
		 PrintResults(list[jmax].params, randIn.param,overlapout, fendBCV, 1);
		 /*j is just a dummy value to not erase jmax; idem for l*/
		 KeepHighestValues(overlapout, j, l, fendBCV, &overlapoutmax, &j, &l, &fMax);
		 
	       }
	     PrintResults(list[jmax].params, randIn.param, overlapoutmax, fendBCV, 1);
	     
	   }
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
   
   LALFree(Filter2.data);
   LALFree(Filter1.data);
   
   
   
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
  coarseIn->LowGM     = BANKEFFICIENCY_LOWGM;
  coarseIn->HighGM    = BANKEFFICIENCY_HIGHGM;
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
      if      ( strcmp(argv[i],"--LowGM")   == 0 ){	  coarseIn->LowGM      = atof(argv[++i]);  }
      else if ( strcmp(argv[i],"--HighGM")  == 0 ){        coarseIn->HighGM     = atof(argv[++i]);  }     
      else if ( strcmp(argv[i],"-fl")       == 0 ){       
	coarseIn->fLower     = atof(argv[++i]);  
	randIn->param.fLower = coarseIn->fLower;
      }	
      else if ( strcmp(argv[i],"-sampling") ==0) {	  
	coarseIn->tSampling  = atof(argv[++i]);  
	randIn->param.fCutoff = coarseIn->tSampling/2-1;
      }      
      else if ( strcmp(argv[i],"-mMin") == 0 )
	{
	  coarseIn->mMin = atof(argv[++i]); 
	  randIn->mMin = coarseIn->mMin;
	  randIn->param.mass1 = randIn->mMin;
	  randIn->param.mass2 = randIn->mMin;
	}
      else if ( strcmp(argv[i], "-m1") == 0 )  {	  otherIn->m1 = atof(argv[++i]);}	
      else if ( strcmp(argv[i], "-m2") == 0 )  {	  otherIn->m2 = atof(argv[++i]); }       
      else if ( strcmp(argv[i], "-psi0") == 0 ){	  otherIn->psi0 = atof(argv[++i]); }       
      else if ( strcmp(argv[i], "-psi3") == 0 ){          otherIn->psi3 = atof(argv[++i]); }	
      else if ( strcmp(argv[i],"-mMax") == 0 ) {
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
      else if ( strcmp(argv[i],"--mm")==0)	     {	  coarseIn->mmCoarse =atof(argv[++i]);	}
      else if ( strcmp(argv[i],"--numFcut")==0)      {	  coarseIn->numFcutTemplates = atoi(argv[++i]);	}
      else if ( strcmp(argv[i],"-simType")==0)	     {	  randIn->type = atoi(argv[++i]);	}
      else if ( strcmp(argv[i],"-sigAmp")==0)	     {	  randIn->SignalAmp = atof(argv[++i]);	}
      else if ( strcmp(argv[i],"--alpha")==0)	     {	  coarseIn->alpha = atof(argv[++i]); 	  randIn->param.alpha = coarseIn->alpha;	}
      else if ( strcmp(argv[i],"--TemplateOrder")==0){     coarseIn->order = atoi(argv[++i]); 	}
      else if ( strcmp(argv[i],"--SignalOrder")==0)  {	  randIn->param.order = atoi(argv[++i]); 	}
      else if ( strcmp(argv[i],"--ntrial")==0)	     {	  otherIn->ntrials = atoi(argv[++i]);       	}
      else if ( strcmp(argv[i],"--n")==0)	     {	  otherIn->ntrials = atoi(argv[++i]);       	}
      else if ( strcmp(argv[i],"--seed")==0)	     {	  randIn->useed = atoi(argv[++i])+1;	}
      else if ( strcmp(argv[i],"--quiet")==0)        {	  otherIn->quietFlag = 1;}
      else if ( strcmp(argv[i],"--noiseModel")==0) 
	{
	  i++;
	  if (strcmp(argv[i], "LIGOI") == 0)       otherIn->NoiseModel = LIGOI;
	  else if (strcmp(argv[i], "LIGOII") == 0) otherIn->NoiseModel = LIGOII;
	  else if (strcmp(argv[i], "VIRGO") == 0)  otherIn->NoiseModel = VIRGO;
	  else if (strcmp(argv[i], "TAMA") == 0)   otherIn->NoiseModel = TAMA;
	  else if (strcmp(argv[i], "GEO") == 0)    otherIn->NoiseModel = GEO;
      }

      else if ( strcmp(argv[i],"--InQuadrature")==0)    { otherIn->overlapMethod = InQuadrature;}
      else if ( strcmp(argv[i],"--AmbiguityFunction")==0){otherIn->ambiguityFunction = 1;}
      else if ( strcmp(argv[i],"--normalisation")==0) {   otherIn->ambiguityFunction = 1;}
      else if ( strcmp(argv[i],"--FMaximization")==0) {	  otherIn->FMaximization = 1;}
      else if ( strcmp(argv[i],"--PrintOverlap")==0)  {	  otherIn->PrintOverlap = 1;}
      else if ( strcmp(argv[i],"--PrintFilter")==0) {	  otherIn->PrintFilter = 1;}
      else if ( strcmp(argv[i],"--check")==0) {	          otherIn->check = 1;}
      else if ( strcmp(argv[i],"--InputPSD")==0) {	  otherIn->inputPSD = argv[++i];}
      else if ( strcmp(argv[i],"--Template")==0)
	{
	  if ( strcmp(argv[++i],"TaylorT1")==0)	            otherIn->template = 0;
	  else if ( strcmp(argv[i],"TaylorT2")==0)	    otherIn->template = 1;
	  else if ( strcmp(argv[i],"TaylorT3")==0)	    otherIn->template = 2;
	  else if ( strcmp(argv[i],"TaylorF1")==0)	    otherIn->template = 3;
	  else if ( strcmp(argv[i],"TaylorF2")==0)	    otherIn->template = 4;
	  else if ( strcmp(argv[i],"PadeT1")==0)	    otherIn->template = 5;
	  else if ( strcmp(argv[i],"PadeF1")==0)	    otherIn->template = 6;
	  else if ( strcmp(argv[i],"EOB")==0)	            otherIn->template = 7;
	  else if ( strcmp(argv[i],"BCV")==0) 	            otherIn->template = 8;
	  else if ( strcmp(argv[i],"SpinTaylorT3")==0)	    otherIn->template = 9;

	  coarseIn->approximant = otherIn->template;
	  if ( coarseIn->approximant == BCV ) 	
	    coarseIn->space = Psi0Psi3; 
	  else 
	    coarseIn->space = Tau0Tau3;
	}	
      else if ( (strcmp(argv[i], "--Signal")==0))
	{
	  if (strcmp(argv[++i],"TaylorT1")==0)	    otherIn->signal = 0;
	  else if (strcmp(argv[i],"TaylorT2")==0)	    otherIn->signal = 1;
	  else if (strcmp(argv[i],"TaylorT3")==0)	    otherIn->signal = 2;
	  else if (strcmp(argv[i],"TaylorF1")==0)	    otherIn->signal = 3;
	  else if (strcmp(argv[i],"TaylorF2")==0)	    otherIn->signal = 4;
	  else if (strcmp(argv[i],"PadeT1")==0)	    otherIn->signal = 5;
	  else if (strcmp(argv[i],"PadeF1")==0)	    otherIn->signal = 6;
	  else if (strcmp(argv[i],"EOB")==0)	    otherIn->signal = 7;
	  else if (strcmp(argv[i],"BCV")==0)	            otherIn->signal = 8;
	  else if (strcmp(argv[i],"SpinTaylorT3")==0)	    otherIn->signal = 9;
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
  fprintf(stderr,"        --ntrial : number of trials                   (%d)\n",       BANKEFFICIENCY_NTRIALS);
  fprintf(stderr,"         --seed  : seed for random generation         (%d)\n",           BANKEFFICIENCY_SEED);
  fprintf(stderr,"       -simType  : type of simulation, 0, 1 or 2      (%d)\n\n",     BANKEFFICIENCY_TYPE);
  fprintf(stderr,"           --fl  : lower frequency cutoff             (%7.2f) Hz\n", BANKEFFICIENCY_FLOWER);
  fprintf(stderr,"          -mMin  : minimal mass of component stars    (%7.2f) Mo\n", BANKEFFICIENCY_MMIN);
  fprintf(stderr,"          -mMax  : maximal mass of component stars    (%7.2f) Mo\n", BANKEFFICIENCY_MMAX);
  fprintf(stderr,"       -psi0Max  : Max value of psi0 (%7.2f)\n",                     BANKEFFICIENCY_PSI0MAX);
  fprintf(stderr,"       -psi3Min  : Min value of psi3 (%7.2f)\n",                     BANKEFFICIENCY_PSI3MIN);
  fprintf(stderr,"         -alpha  : BCV amplitude correction parameter (%7.2f)\n",    BANKEFFICIENCY_ALPHA);
  fprintf(stderr,"      --numFcut  : number of layers in Fcut dimension (%7.2d)\n",    BANKEFFICIENCY_NFCUT);
  fprintf(stderr,"           --mm  : minimal match for template bank    (%7.3f)\n",    BANKEFFICIENCY_MMCOARSE);
  fprintf(stderr,"     --Template  : PN model for template bank, e.g. BCV (%d)\n\n",   BANKEFFICIENCY_TEMPLATE);
  fprintf(stderr,"       --Signal  : model to inject in Monte-Carlo, e.g. EOB            (%d)\n",       BANKEFFICIENCY_SIGNAL);
  fprintf(stderr,"  --SignalOrder  : order of PN model                                   (%7.2d)\n",    BANKEFFICIENCY_ORDER_SIGNAL);
  fprintf(stderr,"--TemplateOrder  : order of signal to injec                            (%7.2d)\n",    BANKEFFICIENCY_ORDER_TEMPLATE);
  fprintf(stderr,"        -sigAmp  : amplitude of the signal                             (%7.2f)\n",    BANKEFFICIENCY_SIGNALAMP);
  fprintf(stderr,"       --HighGM  : Higher distance at which the coalescnce is stopped  (%7.2d)\n",    BANKEFFICIENCY_HIGHGM);
  fprintf(stderr,"        --LowGM  : Lower distance at which the coalescence is stopped  (%7.2d)\n",    BANKEFFICIENCY_LOWGM);



  
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
KeepHighestValues(InspiralWaveOverlapOut in , 
		  INT4 j, INT4 l, double frequency,
		  InspiralWaveOverlapOut *out, 
		  INT4 *jmax, INT4 *lmax, double *fMax)
{
  if (in.max > out->max)
    { 
      out->max   = in.max;
      out->phase = in.phase;
      out->alpha  = in.alpha;
      out->bin   = in.bin;
      *jmax       = j;
      *lmax = l;
      *fMax = frequency;
    }
}

void
PrintResults(
	    InspiralTemplate       bank,
	    InspiralTemplate       injected,
	    InspiralWaveOverlapOut overlapout,
	    double fendBCV,
	    int layer)
{
  
  fprintf(stdout, "%e %e %e %e   ", bank.psi0, injected.psi0, bank.psi3, injected.psi3);
  fprintf(stdout, "%e %e    %e %e   ", fendBCV, injected.fFinal, bank.totalMass, injected.totalMass);
  fprintf(stdout, "%e %e    %e %e ", injected.mass1, injected.mass2, overlapout.max, overlapout.phase);
  fprintf(stdout, "%e %e   %d\n", overlapout.alpha, overlapout.alpha*pow(fendBCV,2./3.), layer);
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
  in.norm =    1./4.;

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

  

  


void LALGetOrthoNormalFilter2(REAL4Vector *filter)
{
  int i,n,nby2;
  REAL4 temp;
  
  n = filter->length;
  nby2 = n/2;


  for (i=1; i<nby2; i++)
    {
      temp = filter->data[i];
      filter->data[i] = - filter->data[n-i];
      filter->data[n-i] = temp;
    }
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



void 
LALCreateFilters(REAL4Vector *Filter1,
		 REAL4Vector *Filter2,
		 REAL4Vector VectorPowerFm5_3,
		 REAL4Vector VectorPowerFm2_3,
		 REAL4Vector VectorPowerFm7_6,
		 REAL4Vector VectorPowerFm1_2,
		 BCVMaximizationMatrix    matrix,
		 INT4 kMin,
		 double psi0,
		 double psi3
		 )
{
  int i;
  double amplitude, cos_phase, sin_phase;
  int n = Filter1->length;
  int nby2 ;
  nby2 = n/2;
  
  /* Create the templates */	       
  for (i = 0; i< kMin-1; i++)
    {
      Filter1->data[i] = 0;
      Filter2->data[i] = 0;
    }
  for (i=kMin; i< nby2; i++)
    {
      
      amplitude  =   psi0 * VectorPowerFm5_3.data[i] +psi3 * VectorPowerFm2_3.data[i];
      
      cos_phase = cos( amplitude);
      sin_phase = sin(amplitude);
            
      amplitude =  matrix.a11 * VectorPowerFm7_6.data[i];
      
      Filter1->data[i]   = amplitude *  cos_phase;      
      Filter1->data[n-i] = -amplitude * sin_phase;
      
      amplitude =  matrix.a21 * VectorPowerFm7_6.data[i] + 
	matrix.a22 * VectorPowerFm1_2.data[i];
      
      Filter2->data[i]   = amplitude *  cos_phase;      
      Filter2->data[n-i] = -amplitude * sin_phase;
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









void
LALWaveOverlapBCV(
			     LALStatus *status,
			     REAL4Vector             *correlation,
			     InspiralWaveOverlapOut  *overlapout,
			     InspiralWaveOverlapIn   *overlapin,
			     REAL4Vector *Filter1,
			     REAL4Vector *Filter2,

			     BCVMaximizationMatrix    matrix,
			     OtherParamIn             otherIn 
		   			     )
     /*  </lalVerbatim>  */
{
  REAL4 maxrho = 0, maxV1=0, maxV2=0;
  INT4 maxbin = 0;
  REAL4Vector x1,x2,x3,x4;
  InspiralWaveCorrelateIn corrin;
  UINT4 i, nBegin, nEnd, n = correlation->length; 
  REAL4 phi, phase;
  REAL4 x1_2,x2_2,x3_2,x4_2, V0,V1, V2, rho=0.;

  x1.length = x2.length = x3.length = x4.length = n;
  
  x1.data = (REAL4*) LALMalloc(sizeof(REAL4) * x1.length);
  x2.data = (REAL4*) LALMalloc(sizeof(REAL4) * x2.length);
  x3.data = (REAL4*) LALMalloc(sizeof(REAL4) * x3.length);
  x4.data = (REAL4*) LALMalloc(sizeof(REAL4) * x4.length);


  
  INITSTATUS (status, "LALInspiralWaveOverlapBCV", BANKEFFICIENCYC);
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
  


    
  
  

  
  corrin.signal2      = *Filter1;
  LALInspiralWaveCorrelate(status->statusPtr, &x1, corrin);
  CHECKSTATUSPTR(status);
  LALGetOrthoNormalFilter2( Filter1);
  corrin.signal2      = *Filter1;
  LALInspiralWaveCorrelate(status->statusPtr, &x3, corrin);
  CHECKSTATUSPTR(status);
  


  corrin.signal2      = *Filter2;
  LALInspiralWaveCorrelate(status->statusPtr, &x2, corrin);
  CHECKSTATUSPTR(status);
  LALGetOrthoNormalFilter2( Filter2);
  corrin.signal2      = *Filter2;
  LALInspiralWaveCorrelate(status->statusPtr, &x4, corrin);
  
  
  /*nbegin et nEnd ici*/
  
  nBegin            = overlapin->nBegin;
  nEnd              = Filter1->length - overlapin->nEnd;
  overlapout->max   = 0.;
  overlapout->bin   = nBegin;
  overlapout->phase = 0;

  if (otherIn.PrintOverlap == 1)
    {
      for (i=0; i< x1.length; i++)
	printf("%d %e\n", i,fabs(x1.data[i]));
      printf("&\n");	       
      
      for (i=1; i< Filter1->length; i++)
	printf("%d %e\n", i,fabs(x2.data[i]));
      printf("&\n");	       
      
      for (i=0; i< Filter1->length; i++)
	printf("%d %e\n", i,fabs(x3.data[i]));
      printf("&\n");	       
      
      for (i=0; i< Filter1->length-1; i++)
	printf("%d %e\n", i,fabs(x4.data[i]));
      printf("&\n");	       
    }




  for (i=nBegin; i<nEnd; i++) 
    {

      x1_2 = x1.data[i] * x1.data[i];
      x2_2 = x2.data[i] * x2.data[i];
      x3_2 = x3.data[i] * x3.data[i];
      x4_2 = x4.data[i] * x4.data[i];      
      
      
      V0 = x1_2 + x2_2 + x3_2 + x4_2;
      V1 = x1_2 + x3_2 - x2_2 - x4_2;
      V2 = 2*x1.data[i]*x2.data[i] + 2*x3.data[i]*x4.data[i];

      rho = V0 + sqrt(V1*V1+V2*V2);

      if ( rho > maxrho)
	{
	  maxV2 = V2;
	  maxV1 = V1;
	  maxrho = V0 + sqrt(V1*V1+V2*V2);
	  maxbin = i;
	}          
    }

  
  overlapout->bin = maxbin;
  overlapout->phase =  .5*atan2(maxV2,maxV1) ;
  overlapout->max = sqrt(maxrho /2.) ; 

  while (overlapout->phase > LAL_PI || overlapout->phase > 0 ) overlapout->phase -= LAL_PI;
  while (overlapout->phase < -LAL_PI || overlapout->phase < 0 ) overlapout->phase += LAL_PI;


  /* Check whether phase is in the appropriate range */


  if (otherIn.PrintOverlap == 1)
    printf("&\n");
  
  phase = overlapout->phase;
  phi   = overlapin->param.startPhase;
  
  
  /*   for (i=0; i<(UINT4)correlation->length; i++) 
       {
       cphase = cos(phase);
       cphi = cos(phi);
       sphase = sqrt(1-cphase*cphase);
       sphi = sqrt(1-cphi*cphi);
       
       correlation->data[i] = x1.data[i] * cphase * cphi
       + x2.data[i] * sphase * cphi
       + x3.data[i] * cphase * sphi
       + x4.data[i] * sphase * sphi;	
       }
  */
  
  overlapout->alpha = -(matrix.a22 * tan(phase)) / (matrix.a11 + matrix.a21* tan(phase));
  
  LALFree(x1.data);
  LALFree(x2.data);
  LALFree(x3.data);
  LALFree(x4.data);
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
  
  
}



void
PrintBank(
	  InspiralCoarseBankIn coarseIn,
	  InspiralTemplateList **list,
	  UINT4 nlist, 
	  UINT4 signalLength)
{
  UINT4 i;
  
  for (i=0; i<nlist; i++)
    {
      if (coarseIn.approximant == BCV)
	{      		   
	  fprintf(stdout, "%e %e %e %e\n", (*list[i]).params.psi0, (*list[i]).params.psi3, (*list[i]).params.totalMass, (*list[i]).params.fFinal);
	} 
      else
	{
	  fprintf(stdout, "%e %e %e %e\n", (*list[i]).params.mass1, (*list[i]).params.mass2, (*list[i]).params.t0, (*list[i]).params.t3);
	}
    }	   
  fprintf(stdout, "#---------------------------------------------------------------\n");
  fprintf(stdout, "#Number of Coarse Bank Templates=%d\n",nlist);
  fprintf(stdout, "#Signal Length=%d\n", signalLength);
  fprintf(stdout,"---------------------------------------------------------------\n");
  if (coarseIn.approximant==BCV)
    {
      fprintf(stdout, "#psi0Min=%e, psi0Max=%e, psi3Min=%e, psi3Max=%e\n", 
	      coarseIn.psi0Min,coarseIn.psi0Max,coarseIn.psi3Min,coarseIn.psi3Max);
      /*	      coarseIn.psi0Min,coarseIn.psi0Max,coarseIn.psi3Min,coarseIn.psi3Max,otherIn.signal,otherIn.template);*/
      fprintf(stdout, "#psi0Tmplt   psi0Sig    psi3Tmplt    psi3Sig    fFinalTmplt   fFinalSig  m1 m2 Overlap phase alpha alpha*f^2/3\n");
    }
  else
    {
      
      fprintf(stdout, "#mMin=%e, mMax=%e\n",
	      coarseIn.mMin,coarseIn.mMax);
      fprintf(stdout, "#psi0Tmplt   psi0Sig    psi3Tmplt    psi3Sig    fFinalTmplt   fFinalSig  m1 m2 Overlap phase alpha alpha*f^2/3\n");
    }
  fprintf(stdout,"---------------------------------------------------------------\n");
}


