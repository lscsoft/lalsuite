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
#define LALINSPIRALH_EDIV0 4
#define LALINSPIRALH_ESIZE 8
#define LALINSPIRALH_ECHOICE 16
#define LALINSPIRALH_MSGENULL "Arguments contained an unexpected null pointer"
#define LALINSPIRALH_MSGEMEM "Memory allocation error"
#define LALINSPIRALH_MSGEDIV0 "Division by zero"
#define LALINSPIRALH_MSGESIZE "Invalid input range"
#define LALINSPIRALH_MSGECHOICE "Invalid choice for an input parameter"


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
\idx[Type]{EtaTau04In}
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
\idx[Type]{EtaTau02In}
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
\idx[Type]{InputMasses}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralHS">  */
typedef enum {
  newtonian,
  oneHalfPN,
  onePN,
  onePointFivePN,
  twoPN,
  twoPointFivePN,
  threePN,
  threePointFivePN
 } Order;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{Order}
</lalLaTeX>  */



/* <lalVerbatim file="LALInspiralHS">  */
typedef enum {
  TaylorT1,
  TaylorT2,
  TaylorT3,
  TaylorF1,
  TaylorF2,
  PadeT1,
  PadeF1,
  EOB,
  DJS,
  INSPA,
  IRSPA
 } Approximant;
/* </lalVerbatim>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagInspiralTemplate
{
  INT4 ieta;
  INT4 level;
  INT4Vector *segmentIdVec;
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
  REAL8 t5; 
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
  REAL8 OmegaS, Theta;
  InputMasses massChoice;
  Order order;
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
tagInspiralToffInput
{
  REAL8 tN;
  REAL8 t2;
  REAL8 t3;
  REAL8 t4;
  REAL8 t5;
  REAL8 t6;
  REAL8 t7;
  REAL8 tl6;
  REAL8 piM;
  REAL8 tc;
  REAL8 t;
 }
InspiralToffInput;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{InspiralToffInput}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct {
   int ieta;
   double eTaN, eTa1, eTa2, eTa3;
   double ePaN, ePa1, ePa2, ePa3;
   double ETaN, ETa1, ETa2, ETa3;
   double dETaN, dETa1, dETa2, dETa3;

   double fTaN, fTa1, fTa2, fTa3, fTa4, fTa5, fTa6, fTa7;
   double fPaN, fPa1, fPa2, fPa3, fPa4, fPa5, fPa6, fPa7;
   double FTaN, FTa1, FTa2, FTa3, FTa4, FTa5, FTa6, FTa7, FTl6;

   double tvaN, tva2, tva3, tva4, tva5, tva6, tva7, tvl6; 
   double pvaN, pva2, pva3, pva4, pva5, pva6, pva7, pvl6; 
   double ptaN, pta2, pta3, pta4, pta5, pta6, pta7, ptl6; 
   double ftaN, fta2, fta3, fta4, fta5, fta6, fta7; 

   double samplingrate, samplinginterval;
   double eta, totalmass, m1, m2;
   double lambda, theta, EulerC, omegaS;
   double f0, fn, t0, tn, v0, vn, vf, vlso, flso, phiC;
   double vlsoT0, vlsoT2, vlsoT4, vlsoT6;
   double vlsoP0, vlsoP2, vlsoP4, vlsoP6;
   double vpoleP4, vpoleP6;
} expnCoeffs;

/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{expnCoeffs}
</lalLaTeX>  */

typedef double EnergyFunction(
   double v, 
   expnCoeffs *ak);

typedef double FluxFunction(
   double v, 
   expnCoeffs *ak);

typedef void (TestFunction)(
   REAL8Vector *vector1,
   REAL8Vector *vector2,
   void *params);

typedef void (InspiralPhasing2)(
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

typedef void (InspiralPhasing3)(
   LALStatus *status,
   REAL8 *f,
   REAL8 td,
   expnCoeffs *ak);

typedef void (InspiralFrequency3)(
   LALStatus *status,
   REAL8 *f,
   REAL8 td,
   expnCoeffs *ak);

typedef void (InspiralTiming2) (
   LALStatus *status,
   REAL8 *toff,
   REAL8 f,
   void *params);


/* <lalVerbatim file="LALInspiralHS">  */
typedef struct
tagexpnFunc
{
  EnergyFunction *dEnergy;
  FluxFunction *flux;
  InspiralTiming2 *timing2;
  InspiralPhasing2 *phasing2;
  InspiralPhasing3 *phasing3;
  InspiralFrequency3 *frequency3;
}
expnFunc;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{expnFunc}
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
\idx[Type]{TofVIntegrandIn}
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
\idx[Type]{InspiralDerivativesIn}
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
\idx[Type]{TofVIn}
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
\idx[Type]{InspiralPhaseIn}
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
\idx[Type]{PhiofVIntegrandIn}
</lalLaTeX>  */

/*  <lalLaTeX>
\vfill{\footnotesize\input{LALInspiralHV}}
</lalLaTeX>  */

/* Function prototypes */

/*  <lalLaTeX>
\newpage\input{LALInspiralChooseModelC}
</lalLaTeX>  */

void LALInspiralChooseModel(
   LALStatus *status,
   expnFunc *func,
   expnCoeffs *ak,
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALEtaTau02C}
</lalLaTeX>  */

void LALEtaTau02(
   LALStatus *status,
   REAL8 *x,
   REAL8 eta, 
   void  *in);

/* <lalLaTeX>
\newpage\input{LALEtaTau04C}
</lalLaTeX>  */

void LALEtaTau04(
   LALStatus *status,
   REAL8 *x,
   REAL8 eta, 
   void  *in);

/*  <lalLaTeX>
\newpage\input{LALInspiralDerivativesC}
</lalLaTeX>  */

void LALInspiralDerivatives (
   REAL8Vector *vec1,
   REAL8Vector *vec2,
   void *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralParameterCalcC}
</lalLaTeX>  */

void LALInspiralParameterCalc (
   LALStatus *status,
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralPhasing1C}
</lalLaTeX>  */

void LALInspiralPhasing1 (
   LALStatus *status,
   REAL8 *phase,
   REAL8 v,
   void *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralSetupC}
</lalLaTeX>  */

void LALInspiralSetup (
   LALStatus *status,
   expnCoeffs *ak,
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralVelocityC}
</lalLaTeX>  */

void LALInspiralVelocity (
   LALStatus *status,
   REAL8 *v,
   TofVIn *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralWaveC}
</lalLaTeX>  */

void LALInspiralWave(
   LALStatus *status,
   REAL4Vector *signal,
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralWaveLengthC}
</lalLaTeX>  */

void LALInspiralWaveLength (
   LALStatus *status,
   INT4 *n,
   InspiralTemplate params);

/*  <lalLaTeX>
\newpage\input{LALInspiralWaveTemplatesC}
</lalLaTeX>  */

void LALInspiralWaveTemplates(
   LALStatus *status,
   REAL4Vector *filter1,
   REAL4Vector *filter2,
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralPhiofVIntegrandC}
</lalLaTeX>  */

void LALInspiralPhiofVIntegrand (
   LALStatus *status,
   REAL8 *,
   REAL8,
   void *);


/*  <lalLaTeX>
\newpage\input{LALInspiralWave1C}
</lalLaTeX>  */

void LALInspiralWave1(
   LALStatus *status,
   REAL4Vector *signal,
   InspiralTemplate *params);


/*  <lalLaTeX>
\newpage\input{LALInspiralWave1TemplatesC}
</lalLaTeX>  */

void LALInspiralWave1Templates(
   LALStatus *status,
   REAL4Vector *signal1,
   REAL4Vector *signal2,
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralTofVC}
</lalLaTeX>  */

void LALInspiralTofV (
   LALStatus *,
   REAL8 *,
   REAL8,
   void *);


/*  <lalLaTeX>
\newpage\input{LALInspiralTofVIntegrandC}
</lalLaTeX>  */

void LALInspiralTofVIntegrand (
   LALStatus *status,
   REAL8 *,
   REAL8,
   void *);

void LALPrintTimeseries(
   int n,
   double *signal,
   double dt,
   double t0);

/*  <lalLaTeX>
\newpage\input{LALRungeKutta4C}
</lalLaTeX>  */

void LALRungeKutta4(
  LALStatus *,
  REAL8Vector *,
  rk4In *,
  void *);

/*  <lalLaTeX>
\newpage\input{LALInspiralWave2C}
</lalLaTeX>  */

void LALInspiralWave2(
  LALStatus *status,
  REAL4Vector *signal,
  InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralWave2TemplatesC}
</lalLaTeX>  */

void LALInspiralWave2Templates (
  LALStatus *status,
  REAL4Vector *signal1, 
  REAL4Vector *signal2, 
  InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralTiming2C}
</lalLaTeX>  */

void LALInspiralTiming2_0PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_1PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_2PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_3PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_4PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);


void LALInspiralTiming2_5PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_6PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_7PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralPhasing3C}
</lalLaTeX>  */

void LALInspiralPhasing2_0PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_1PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_2PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_3PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_4PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_5PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_6PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_7PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

/*  <lalLaTeX>
\newpage\input{LALInspiralWave3C}
</lalLaTeX>  */

void LALInspiralWave3 (
   LALStatus *status,
   REAL4Vector *signal, 
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralWave3TemplatesC}
</lalLaTeX>  */

void LALInspiralWave3Templates (
   LALStatus *status,
   REAL4Vector *signal1, 
   REAL4Vector *signal2, 
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralFrequency3C}
</lalLaTeX>  */

void LALInspiralFrequency3_0PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_1PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_2PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_3PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_4PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_5PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_6PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_7PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

/*  <lalLaTeX>
\newpage\input{LALInspiralPhasing3C}
</lalLaTeX>  */

void LALInspiralPhasing3_0PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_1PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_2PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_3PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_4PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_5PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_6PN (
   LALStatus *,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_7PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALRungeKutta4(
  LALStatus *,
  REAL8Vector *,
  rk4In *,
  void *);

void LALEOBWaveform(
	LALStatus *status,
	REAL4Vector *signal,
	InspiralTemplate *params);

void LALHCapDerivatives(
	REAL8Vector *values,
	REAL8Vector *dvalues,
	void *funcParams);

void LALprInit(
	REAL8 *pr,
	REAL8 r,
	InspiralDerivativesIn *ak);

void LALpphiInit(
	REAL8 *phase,
	REAL8 r,
	REAL8 eta);

void LALlightRingRadius(
	LALStatus *status,
	REAL8 *x,
	REAL8 r,
	void *params);

void LALrOfOmega(
	LALStatus *status,
	REAL8 *x,
	REAL8 r,
	void *params);



/*  <lalLaTeX>
\newpage\input{LALInspiralWave2C}
</lalLaTeX>  */

void LALInspiralWave2(
  LALStatus *status,
  REAL4Vector *signal,
  InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralWave2TemplatesC}
</lalLaTeX>  */

void LALInspiralWave2Templates (
  LALStatus *status,
  REAL4Vector *signal1, 
  REAL4Vector *signal2, 
  InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralTiming2C}
</lalLaTeX>  */

void LALInspiralTiming2_0PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_1PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_2PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_3PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_4PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);


void LALInspiralTiming2_5PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_6PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

void LALInspiralTiming2_7PN (
   LALStatus *,
   REAL8 *toff,
   REAL8 f,
   void *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralPhasing3C}
</lalLaTeX>  */

void LALInspiralPhasing2_0PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_1PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_2PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_3PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_4PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_5PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_6PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

void LALInspiralPhasing2_7PN (
   LALStatus *status,
   REAL8 *phase, 
   REAL8 v,
   expnCoeffs *ak);

/*  <lalLaTeX>
\newpage\input{LALInspiralWave3C}
</lalLaTeX>  */

void LALInspiralWave3 (
   LALStatus *status,
   REAL4Vector *signal, 
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralWave3TemplatesC}
</lalLaTeX>  */

void LALInspiralWave3Templates (
   LALStatus *status,
   REAL4Vector *signal1, 
   REAL4Vector *signal2, 
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralFrequency3C}
</lalLaTeX>  */

void LALInspiralFrequency3_0PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_1PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_2PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_3PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_4PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_5PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_6PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralFrequency3_7PN (
   LALStatus *status,
   REAL8 *frequency,
   REAL8 td,
   expnCoeffs *ak);

/*  <lalLaTeX>
\newpage\input{LALInspiralPhasing3C}
</lalLaTeX>  */

void LALInspiralPhasing3_0PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_1PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_2PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_3PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_4PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_5PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_6PN (
   LALStatus *,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);

void LALInspiralPhasing3_7PN (
   LALStatus *status,
   REAL8 *phase,
   REAL8 td,
   expnCoeffs *ak);


/*  <lalLaTeX>
\newpage\input{LALInspiralTestC}
</lalLaTeX>  */

#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRAL_H */
