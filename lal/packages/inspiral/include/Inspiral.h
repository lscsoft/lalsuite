/* 
    Date: Februray 15, 2000
    Author: B.S. Sathyaprakash, Cardiff University
    Purpose: Definition structures in T- and P-approximant codes.
*/
/* Physical and mathematical constants */

#ifndef _INSPIRAL_H
#define _INSPIRAL_H

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <lal/LALStdlib.h>
# include <lal/LALConstants.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( INSPIRALH, "$Id$" );

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


#define INSPIRALWAVE_ENULL 1
#define INSPIRALSETUP_ENULL 1
#define CHOOSEMODEL_ENULL 1
#define CHOOSEMODEL_ECHOICE 2
#define TIMEDOMAIN2_ENULL 1
#define INSPIRALVELOCITY_ENULL 1
#define TOFV_ENULL 1
#define TOFVINTEGRAND_ENULL 1
#define INSPIRALPHASE_ENULL 1
#define PHIOFVINTEGRAND_ENULL 1
#define INSPIRALDERIVATIVES_ENULL 1
#define RUNGEKUTTA4_ENULL 1
#define TAPPRPNTDOMTIME_ENULL 1
#define TAPPRPNTDOMTIME_ESIZE 2
#define TAPPRPNTDOMTIMEFREQUENCY_ENULL 1
#define TAPPRPNTDOMTIMEFREQUENCY_ESIZE 2
#define TAPPRPNTDOMTIMEPHASE_ENULL 1
#define TAPPRPNTDOMTIMEPHASE_ESIZE 2
#define INSPIRALPARAMETERCALCALC_ENULL 1
#define TAPPRPNTDOMFREQ_ENULL 1
#define TAPPRPNTDOMFREQ_ESIZE 2
#define TAPPRPNTDOMFREQPHASE_ENULL 1
#define TAPPRPNTDOMFREQPHASE_ESIZE 2
#define TAPPRPNTDOMFREQTOFF_ENULL 1
#define TAPPRPNTDOMFREQTOFF_ESIZE 2



#define INSPIRALWAVE_MSGENULL "Null pointer"
#define INSPIRALSETUP_MSGENULL "Null pointer"
#define CHOOSEMODEL_MSGENULL "Null pointer"
#define CHOOSEMODEL_MSGECHOICE "Invalid choice of input integer made"
#define TIMEDOMAIN2_MSGENULL "Null pointer"
#define INSPIRALVELOCITY_MSGENULL "Null pointer"
#define TOFV_MSGENULL "Null pointer"
#define TOFVINTEGRAND_MSGENULL "Null pointer"
#define INSPIRALPHASE_MSGENULL "Null pointer"
#define PHIOFVINTEGRAND_MSGENULL "Null pointer"
#define INSPIRALDERIVATIVES_MSGENULL "Null pointer"
#define RUNGEKUTTA4_MSGENULL "Null pointer"
#define TAPPRPNTDOMTIME_MSGENULL "Null pointer"
#define TAPPRPNTDOMTIME_MSGESIZE "Invalid input range"
#define TAPPRPNTDOMTIMEFREQUENCY_MSGENULL "Null pointer"
#define TAPPRPNTDOMTIMEFREQUENCY_MSGESIZE "Invalid input range"
#define TAPPRPNTDOMTIMEPHASE_MSGENULL "Null pointer"
#define TAPPRPNTDOMTIMEPHASE_MSGESIZE "Invalid input range"
#define INSPIRALPARAMETERCALCALC_MSGENULL "Null pointer"
#define TAPPRPNTDOMFREQ_MSGENULL "Null pointer"
#define TAPPRPNTDOMFREQ_MSGESIZE "Invalid input range"
#define TAPPRPNTDOMFREQPHASE_MSGENULL "Null pointer"
#define TAPPRPNTDOMFREQPHASE_MSGESIZE "Invalid input range"
#define TAPPRPNTDOMFREQTOFF_MSGENULL "Null pointer"
#define TAPPRPNTDOMFREQTOFF_MSGESIZE "Invalid input range"


typedef enum {
  m1Andm2,
  totalMassAndEta,
  totalMassAndMu
 } InputMasses;

typedef enum {
  time,
  frequency
 } Domain;

typedef enum {
  newtonian,
  oneHalfPN,
  onePN,
  onePointFivePN,
  twoPN,
  twoPointFivePN
 } Order;



typedef enum {
  taylor,
  pade
 } Approximant;


typedef enum {
  one,
  two,
  three,
  best
 } Method;



typedef struct
tagInspiralTemplate
{
  INT4 number;
  REAL8 mass1; 
  REAL8 mass2;
  REAL8 spin1[3];
  REAL8 spin2[3];
  REAL8 inclination;
  REAL8 eccentricity;
  REAL8 totalMass; 
  REAL8 mu; 
  REAL8 eta;
  REAL8 fLower;
  REAL8 fCutoff;
  REAL8 tSampling;
  REAL8 startPhase;
  REAL8 startTime;
  REAL8 signalAmplitude;
  INT4 nStartPad;
  INT4 ieta;
  INT4 nEndPad;
  InputMasses massChoice;
  Method method;
  Order order;
  Domain domain;
  Approximant approximant;
  struct tagInspiralTemplate *next;
  struct tagInspiralTemplate *fine;
 }
 InspiralTemplate;



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


typedef struct
tagInspiralwavePhaseOutput
{
  REAL8 phase;
}
InspiralwavePhaseOutput;


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


typedef struct
tagInspiralwaveFrequencyOutput
{
  REAL8 frequency;
}
InspiralwaveFrequencyOutput;



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


typedef double EnergyFunction(double, expnCoeffs *);
typedef double FluxFunction(double, expnCoeffs *);
typedef void (TestFunction) (LALStatus *,REAL8Vector *,REAL8Vector *,void *);


typedef struct
tagexpnFunc
{
  EnergyFunction *dEnergy;
  FluxFunction *flux;
}
expnFunc;


typedef struct
tagTofVIntegrandIn
{
  EnergyFunction *dEnergy;
  FluxFunction *flux;
  expnCoeffs *coeffs;
}
TofVIntegrandIn;

typedef struct
tagInspiralDerivativesIn
{
  REAL8 totalmass;
  EnergyFunction *dEnergy;
  FluxFunction *flux;
  expnCoeffs *coeffs;
}
InspiralDerivativesIn;

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


typedef struct
tagPhiofVIntegrandIn
{
  EnergyFunction *dEnergy;
  FluxFunction *flux;
  expnCoeffs *coeffs;
}
PhiofVIntegrandIn;








void LALInspiralSetup (LALStatus *,
		    expnCoeffs *,
		    InspiralTemplate *);

void LALChooseModel(LALStatus *,
		 expnFunc *,
		 expnCoeffs *,
		 InspiralTemplate *);

/*
void  printf_timeseries(int n, double *data, double dt, double t0);
*/

void LALTofV (LALStatus *,
	   REAL8 *,
	   REAL8,
	   void *);

void LALInspiralVelocity (LALStatus *,
	               REAL8 *,
	               TofVIn *);

void LALInspiralDerivatives (LALStatus *,
			  REAL8Vector *,
			  REAL8Vector *,
			  void *);

void LALRungeKutta4(LALStatus *,
	         REAL8Vector *,
	         rk4In *,
	         void *);

void LALTimeDomain2(LALStatus *,
		 REAL8Vector *,
		 InspiralTemplate *);



void LALInspiralPhase (LALStatus *,
	            REAL8 *,
	            REAL8,
	            void *);

void LALPhiofVIntegrand (LALStatus *,
		    REAL8 *,
		    REAL8,
		    void *);

void LALTofVIntegrand (LALStatus *,
		    REAL8 *,
		    REAL8,
		    void *);


void LALInspiralWave(LALStatus *,
		  REAL8Vector *,
		  InspiralTemplate *);

void timedomain(LALStatus *,
		REAL8Vector *,
		InspiralTemplate *,
	        TimeDomainIn *);


void TappRpnTdomTime0PN (LALStatus *,
                         REAL8Vector *, 
		         InspiralTemplate *);
void TappRpnTdomTime1PN (LALStatus *,
                         REAL8Vector *, 
		         InspiralTemplate *);
void TappRpnTdomTime2PN (LALStatus *,
                         REAL8Vector *, 
		         InspiralTemplate *);
void TappRpnTdomTime3PN (LALStatus *,
                         REAL8Vector *, 
		         InspiralTemplate *);
void TappRpnTdomTime4PN (LALStatus *,
                         REAL8Vector *, 
		         InspiralTemplate *);

void LALInspiralParameterCalc (LALStatus *,
			     InspiralParamsOutput *,
			     InspiralParamsInput *);

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

void TappRpnTdomFreq0PN (LALStatus *, 
                         REAL8Vector *, 
                         InspiralTemplate *);

void TappRpnTdomFreq1PN (LALStatus *, 
                         REAL8Vector *, 
                         InspiralTemplate *);

void TappRpnTdomFreq2PN (LALStatus *, 
                         REAL8Vector *, 
                         InspiralTemplate *);

void TappRpnTdomFreq3PN (LALStatus *, 
                         REAL8Vector *, 
                         InspiralTemplate *);

void TappRpnTdomFreq4PN (LALStatus *, 
                         REAL8Vector *, 
                         InspiralTemplate *);

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

void LALTappRpnTdomTime (LALStatus *status,
                      REAL8Vector *signal, 
                      InspiralTemplate *params);
void LALTappRpnTdomFreq (LALStatus *status,
                      REAL8Vector *signal, 
                      InspiralTemplate *params);

#ifdef  __cplusplus
}
#endif

#endif /* _INSPIRAL_H */
