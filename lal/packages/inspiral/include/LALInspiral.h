/* <lalVerbatim file="LALInspiralHV">

Author: Churches, D. K and B. S. Sathyaprakash.
$Id$

</lalVerbatim>  */

/*  <lalLaTeX>

\section{Header \texttt{LALInspiral.h}}
\label{s:LALInspiral.h}

Header file for the template generation codes.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/LALInspiral.h>
\end{verbatim}

\noindent This header file covers routines that are used in template generation.

</lalLaTeX> */

#ifndef _LALINSPIRAL_H
#define _LALINSPIRAL_H

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <lal/LALStdlib.h>
# include <lal/LALConstants.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( LALINSPIRALH, "$Id$" );

# define oneby3    0.333333333333333333333333333
# define twoby3    0.666666666666666666666666667
# define fourby3   1.333333333333333333333333333
# define fiveby3   1.666666666666666666666666667
# define sevenby3  2.333333333333333333333333333
# define eightby3  2.666666666666666666666666667
# define tenby3    3.333333333333333333333333333
# define elevenby3 3.666666666666666666666666666
# define threeby8  0.375
# define fiveby8   0.625
# define threeby4  0.75
# define sevenby8  0.875

/*  <lalLaTeX>

\subsection*{Error codes}

</lalLaTeX>  */

/* <lalErrTable> */

#define LALINSPIRALH_ENULL 1
#define LALINSPIRALH_EMEM 2
#define LALINSPIRALH_ECHOICE 3
#define LALINSPIRALH_EDIV0 4
#define LALINSPIRALH_ESIZE 8
#define LALINSPIRALH_MSGENULL "Arguments contained an unexpected null pointer"
#define LALINSPIRALH_MSGEMEM "Memory allocation error"
#define LALINSPIRALH_MSGECHOICE "Invalid choice for an input parameter"
#define LALINSPIRALH_MSGEDIV0 "Division by zero"
#define LALINSPIRALH_MSGESIZE "Invalid input range"


/* </lalErrTable> */


/* <lalLaTeX>

\section*{Structures}
\input{LALInspiralHS}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagEtaTau04In
{
   REAL8 t4;
   REAL8 A4;
   REAL8 B4;
   REAL8 C4;
} EtaTau04In;   
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{EtaTau04In}}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagEtaTau02In
{
   REAL8 t2;
   REAL8 A2;
   REAL8 B2;
} EtaTau02In;   
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{EtaTau02In}}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralHS">  */
typedef enum {
  m1Andm2,
  totalMassAndEta,
  totalMassAndMu,
  t01,
  t02,
  t03,
  t04,
 } InputMasses;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InputMasses}}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralHS">  */
typedef enum {
  TimeDomain,
  FrequencyDomain
 } Domain;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{Domain}}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralHS">  */
typedef enum {
  newtonian,
  oneHalfPN,
  onePN,
  onePointFivePN,
  twoPN,
  twoPointFivePN
 } Order;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{Order}}
</lalLaTeX>  */



/* <lalVerbatim file="LALInspiralHS">  */
typedef enum {
  taylor,
  pade
 } Approximant;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{Approximant}}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef enum {
  one,
  two,
  three,
  eob,
  best
 } Method;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{Method}}
</lalLaTeX>  */



/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralTemplate
{
  INT4 ieta;
  INT4 level;
  INT4 number;
  INT4 nStartPad;
  INT4 nEndPad;
  REAL8 mass1; 
  REAL8 mass2;
  REAL8 spin1[3];
  REAL8 spin2[3];
  REAL8 inclination;
  REAL8 eccentricity;
  REAL8 totalMass; 
  REAL8 chirpMass; 
  REAL8 t0; 
  REAL8 t2; 
  REAL8 t3; 
  REAL8 t4; 
  REAL8 tC; 
  REAL8 mu; 
  REAL8 eta;
  REAL8 fLower;
  REAL8 fCutoff;
  REAL8 tSampling;
  REAL8 startPhase;
  REAL8 startTime;
  REAL8 signalAmplitude;
  REAL8 rInitial;
  REAL8 vInitial;
  REAL8 rFinal;
  REAL8 vFinal;
  REAL8 rLightRing;
  InputMasses massChoice;
  Method method;
  Order order;
  Domain domain;
  Approximant approximant;
  struct tagInspiralTemplate *next;
  struct tagInspiralTemplate *fine;
 }
 InspiralTemplate;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InspiralTemplate}}
</lalLaTeX>  */



/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralwavePhaseInput
{
  REAL8 etaby2;
  REAL8 td;
  REAL8 a2;
  REAL8 a3;
  REAL8 a4;
}
InspiralwavePhaseInput;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InspiralwavePhaseInput}}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralwavePhaseOutput
{
  REAL8 phase;
}
InspiralwavePhaseOutput;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InspiralwavePhaseOutput}}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralwaveFrequencyInput
{
  REAL8 td;
  REAL8 eightPiM;
  REAL8 a2;
  REAL8 a3;
  REAL8 a4;
}
InspiralwaveFrequencyInput;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InspiralwaveFrequencyInput}}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralwaveFrequencyOutput
{
  REAL8 frequency;
}
InspiralwaveFrequencyOutput;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InspiralwaveFrequencyOutput}}
</lalLaTeX>  */



/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralParamsInput
{
  REAL8 m1; 
  REAL8 m2; 
  REAL8 totalMass; 
  REAL8 mu; 
  REAL8 eta;
  REAL8 fLower;
  InputMasses massChoice;
  Order order;
 }
InspiralParamsInput;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InspiralParamsInput}}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralParamsOutput
{
  REAL8 m1; 
  REAL8 m2; 
  REAL8 totalMass; 
  REAL8 mu; 
  REAL8 eta;
  REAL8 chirpMass;
  REAL8 tau0;
  REAL8 tau2;
  REAL8 tau3;
  REAL8 tau4;
  REAL8 tauC;
 }
InspiralParamsOutput;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InspiralParamsOutput}}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralPhasesInput
{
  REAL8 p0;
  REAL8 p2;
  REAL8 p3;
  REAL8 p4;
  REAL8 pc;
  REAL8 f;
 }
InspiralPhasesInput;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InspiralPhasesInput}}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralToffInput
{
  REAL8 t0;
  REAL8 t2;
  REAL8 t3;
  REAL8 t4;
  REAL8 tc;
  REAL8 t;
 }
InspiralToffInput;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InspiralToffInput}}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct {
   int n;
   int ieta;

   double eTaN;
   double eTa0;
   double eTa1;
   double eTa2;
   double eTa3;
   double eTa4;
   double eTa5;
   double ePaN;
   double ePa0;
   double ePa1;
   double ePa2;
   double ePa3;
   double ePa4;
   double ePa5;

   double ETaN;
   double ETa0;
   double ETa1;
   double ETa2;
   double ETa3;
   double ETa4;
   double ETa5;

   double dETaN;
   double dETa0;
   double dETa1;
   double dETa2;
   double dETa3;
   double dETa4;
   double dETa5;

   double fTaN;
   double fTa0;
   double fTa1;
   double fTa2;
   double fTa3;
   double fTa4;
   double fTa5;
   double fPaN;
   double fPa0;
   double fPa1;
   double fPa2;
   double fPa3;
   double fPa4;
   double fPa5;

   double FTaN;
   double FTa0;
   double FTa1;
   double FTa2;
   double FTa3;
   double FTa4;
   double FTa5;

   double samplingrate;
   double samplinginterval;

   double eta;
   double totalmass;
   double m1;
   double m2;

   double f0;
   double fn;
   double t0;
   double tn;
   double v0;
   double vn;

   double vf;

   double vlso;
   double vlsoT0;
   double vlsoT2;
   double vlsoT4;
   double vlsoP0;
   double vlsoP2;
   double vlsoP4;
   double vpoleP4;
} expnCoeffs;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{expnCoeffs}}
</lalLaTeX>  */


typedef double EnergyFunction(double, expnCoeffs *);
typedef double FluxFunction(double, expnCoeffs *);
typedef void (TestFunction) (REAL8Vector *,REAL8Vector *,void *);


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagexpnFunc
{
  EnergyFunction *dEnergy;
  FluxFunction *flux;
}
expnFunc;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{expnFunc}}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagTofVIntegrandIn
{
  EnergyFunction *dEnergy;
  FluxFunction *flux;
  expnCoeffs *coeffs;
}
TofVIntegrandIn;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{TofVIntegrandIn}}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralDerivativesIn
{
  REAL8 totalmass;
  EnergyFunction *dEnergy;
  FluxFunction *flux;
  expnCoeffs *coeffs;
}
InspiralDerivativesIn;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InspiralDerivativesIn}}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagrk4In
{
  TestFunction *function;
  REAL8 x;
  REAL8Vector *y;
  REAL8Vector *dydx;
  REAL8Vector *yt;
  REAL8Vector *dym;
  REAL8Vector *dyt;
  REAL8 h;
  INT4 n;
}
rk4In;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{rk4In}}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagTimeDomainIn
{
  REAL8 v0;
  REAL8 phi0;
  REAL8 vlso;
  REAL8 t0;
  REAL8 totalmass;
  EnergyFunction *dEnergy;
  FluxFunction *flux;
  expnCoeffs *coeffs;
}
TimeDomainIn;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{TimeDomainIn}}
</lalLaTeX>  */



/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagTofVIn
{
  REAL8 t;
  REAL8 v0;
  REAL8 t0;
  REAL8 vlso;
  REAL8 totalmass;
  EnergyFunction *dEnergy;
  FluxFunction *flux;
  expnCoeffs *coeffs;
}
TofVIn;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{TofVIn}}
</lalLaTeX>  */




/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralPhaseIn
{
  REAL8 v0;
  REAL8 phi0;
  EnergyFunction *dEnergy;
  FluxFunction *flux;
  expnCoeffs *coeffs;
}
InspiralPhaseIn;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{InspiralPhaseIn}}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagPhiofVIntegrandIn
{
  EnergyFunction *dEnergy;
  FluxFunction *flux;
  expnCoeffs *coeffs;
}
PhiofVIntegrandIn;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{PhiofVIntegrandIn}}
</lalLaTeX>  */

/*  <lalLaTeX>
\vfill{\footnotesize\input{LALInspiralHV}}
</lalLaTeX>  */



/* Function prototypes */

/*  <lalLaTeX>
\newpage\input{LALChooseModelC}
</lalLaTeX>  */

void LALChooseModel(LALStatus *,
		    expnFunc *,
		    expnCoeffs *,
		    InspiralTemplate *);

/*  <lalLaTeX>
\newpage\input{LALEtaTau02C}
</lalLaTeX>  */

void LALEtaTau02(LALStatus *status,
                 REAL8 *x,
                 REAL8 eta, 
                 void  *in);

/* <lalLaTeX>
\newpage\input{LALEtaTau04C}
</lalLaTeX>  */

void LALEtaTau04(LALStatus *status,
                 REAL8 *x,
                 REAL8 eta, 
                 void  *in);

/*  <lalLaTeX>
\newpage\input{LALInspiralDerivativesC}
</lalLaTeX>  */

void LALInspiralDerivatives (REAL8Vector *,
			     REAL8Vector *,
			     void *);

/*  <lalLaTeX>
\newpage\input{LALInspiralParameterCalcC}
</lalLaTeX>  */

void LALInspiralParameterCalc (LALStatus *,
			       InspiralTemplate *);

/*  <lalLaTeX>
\newpage\input{LALInspiralPhaseC}
</lalLaTeX>  */

void LALInspiralPhase (LALStatus *,
	               REAL8 *,
	               REAL8,
	               void *);

/*  <lalLaTeX>
\newpage\input{LALInspiralSetupC}
</lalLaTeX>  */

void LALInspiralSetup (LALStatus *,
		       expnCoeffs *,
		       InspiralTemplate *);

/*  <lalLaTeX>
\newpage\input{LALInspiralVelocityC}
</lalLaTeX>  */

void LALInspiralVelocity (LALStatus *,
	                  REAL8 *,
	                  TofVIn *);

/*  <lalLaTeX>
\newpage\input{LALInspiralWaveC}
</lalLaTeX>  */

void LALInspiralWave(LALStatus *,
         	     REAL4Vector *,
		     InspiralTemplate *);

/*  <lalLaTeX>
\newpage\input{LALInspiralWaveLengthC}
</lalLaTeX>  */

void LALInspiralWaveLength (LALStatus *,
			    INT4 *,
			    InspiralTemplate);

/*  <lalLaTeX>
\newpage\input{LALInspiralWaveTemplatesC}
</lalLaTeX>  */

void LALInspiralWaveTemplates(LALStatus *,
		              REAL4Vector *,
		              REAL4Vector *,
		              InspiralTemplate *);

/*  <lalLaTeX>
\newpage\input{LALPhiofVIntegrandC}
</lalLaTeX>  */

void LALPhiofVIntegrand (LALStatus *,
		         REAL8 *,
		         REAL8,
		         void *);


/*  <lalLaTeX>
\newpage\input{LALTimeDomain2C}
</lalLaTeX>  */

void LALTimeDomain2(LALStatus *,
		    REAL4Vector *,
		    InspiralTemplate *);


/*  <lalLaTeX>
\newpage\input{LALTimeDomain2TemplatesC}
</lalLaTeX>  */

void LALTimeDomain2Templates(LALStatus *,
		             REAL4Vector *,
		             REAL4Vector *,
		 	     InspiralTemplate *);


/*  <lalLaTeX>
\newpage\input{LALTofVC}
</lalLaTeX>  */

void LALTofV (LALStatus *,
	      REAL8 *,
	      REAL8,
	      void *);


/*  <lalLaTeX>
\newpage\input{LALTofVIntegrandC}
</lalLaTeX>  */

void LALTofVIntegrand (LALStatus *,
		       REAL8 *,
		       REAL8,
		       void *);

void LALPrintTimeseries(int n,
                         double *signal,
                         double dt,
                         double t0);

/*  <lalLaTeX>
\newpage\input{LALRungeKutta4C}
</lalLaTeX>  */

void LALRungeKutta4(LALStatus *,
	            REAL8Vector *,
	            rk4In *,
	            void *);

/*  <lalLaTeX>
\newpage\input{LALTappRpnTdomFreqC}
</lalLaTeX>  */

void LALTappRpnTdomFreq (LALStatus *status,
                         REAL4Vector *signal, 
                         InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALTappRpnTdomFreqTemplatesC}
</lalLaTeX>  */

void LALTappRpnTdomFreqTemplates (LALStatus *status,
                                  REAL4Vector *signal1, 
                                  REAL4Vector *signal2, 
                                  InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALTappRpnTdomFreqTofFC}
</lalLaTeX>  */

void LALTappRpnTdomFreqTofF0PN (LALStatus *,
                                REAL8 *,
			        REAL8,
                                void *);

void LALTappRpnTdomFreqTofF1PN (LALStatus *,
                                REAL8 *,
			        REAL8,
                                void *);

void LALTappRpnTdomFreqTofF2PN (LALStatus *,
                                REAL8 *,
			        REAL8,
                                void *);

void LALTappRpnTdomFreqTofF3PN (LALStatus *,
                                REAL8 *,
			        REAL8,
                                void *);

void LALTappRpnTdomFreqTofF4PN (LALStatus *,
                                REAL8 *,
			        REAL8,
                                void *);

/*  <lalLaTeX>
\newpage\input{LALTappRpnTdomFreqPhaseC}
</lalLaTeX>  */

void LALTappRpnTdomFreqPhase0PN (LALStatus *,
                                REAL8 *, 
                                InspiralPhasesInput *);

void LALTappRpnTdomFreqPhase1PN (LALStatus *,
                                 REAL8 *, 
                                 InspiralPhasesInput *);

void LALTappRpnTdomFreqPhase2PN (LALStatus *,
                                 REAL8 *, 
                                 InspiralPhasesInput *);

void LALTappRpnTdomFreqPhase3PN (LALStatus *,
                                 REAL8 *, 
                                 InspiralPhasesInput *);

void LALTappRpnTdomFreqPhase4PN (LALStatus *,
                                 REAL8 *, 
                                 InspiralPhasesInput *);

/*  <lalLaTeX>
\newpage\input{LALTappRpnTdomTimeC}
</lalLaTeX>  */

void LALTappRpnTdomTime (LALStatus *status,
                         REAL4Vector *signal, 
                         InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALTappRpnTdomTimeTemplatesC}
</lalLaTeX>  */

void LALTappRpnTdomTimeTemplates (LALStatus *status,
                                  REAL4Vector *signal1, 
                                  REAL4Vector *signal2, 
                                  InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALTappRpnTdomTimeFrequencyC}
</lalLaTeX>  */

void LALTappRpnTdomTimeFrequency0PN (LALStatus *,
                                     InspiralwaveFrequencyOutput *,
			             InspiralwaveFrequencyInput *);

void LALTappRpnTdomTimeFrequency1PN (LALStatus *,
                                     InspiralwaveFrequencyOutput *,
			             InspiralwaveFrequencyInput *);

void LALTappRpnTdomTimeFrequency2PN (LALStatus *,
                                     InspiralwaveFrequencyOutput *,
			             InspiralwaveFrequencyInput *);

void LALTappRpnTdomTimeFrequency3PN (LALStatus *,
                                     InspiralwaveFrequencyOutput *,
			             InspiralwaveFrequencyInput *);

void LALTappRpnTdomTimeFrequency4PN (LALStatus *,
                                     InspiralwaveFrequencyOutput *,
			             InspiralwaveFrequencyInput *);

/*  <lalLaTeX>
\newpage\input{LALTappRpnTdomTimePhaseC}
</lalLaTeX>  */

void LALTappRpnTdomTimePhase0PN (LALStatus *,
                                 InspiralwavePhaseOutput *,
			         InspiralwavePhaseInput *);

void LALTappRpnTdomTimePhase1PN (LALStatus *,
                                 InspiralwavePhaseOutput *,
			         InspiralwavePhaseInput *);

void LALTappRpnTdomTimePhase2PN (LALStatus *,
                                 InspiralwavePhaseOutput *,
			         InspiralwavePhaseInput *);

void LALTappRpnTdomTimePhase3PN (LALStatus *,
                                 InspiralwavePhaseOutput *,
			         InspiralwavePhaseInput *);

void LALTappRpnTdomTimePhase4PN (LALStatus *,
                                 InspiralwavePhaseOutput *,
			         InspiralwavePhaseInput *);


/*  <lalLaTeX>
\newpage\input{LALInspiralTestC}
</lalLaTeX>  */

#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRAL_H */
