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

\begin{enumerate} 

\item \texttt{EtaTau02In, EtaTau04In:}
These are the input structures needed in solving for the mass
ratio $\eta$ given the chirptimes $\tau_0$ and $\tau_2,$ or
to solve for $\eta$ given the chirptimes $\tau_2$ and $\tau_4.$

\input{LALEtaTau0Tau2InH}

Here, \texttt{t2}~$ = \tau_2,$ \texttt{A2}~$ = A_2 ({\tau_0}/{A_0})^{3/5},$  and 
\texttt{B2}~$=B_2$, 
where $A_0 = 5/[256 (\pi f_{s} )^{8/3}],$ $A_2 = 3715 / [64512 (\pi f_s)^2],$
$B_2 = 4620/3715.$  

Similarly, \texttt{t4}~$ = \tau_4,$ \texttt{A4}~$ = A_4 ({\tau_0}/{A_0})^{1/5},$  
\texttt{B4}~$=B_4$ and \texttt{C4}~$=C_4,$ where 
where $A_0 = 5/[256 (\pi f_{s} )^{8/3}],$ 
$A_4 = 5 \times 3058673/ [128 \times 1016064  (\pi f_s)^{4/3}],$
$B_4 = 5429 \times 1016064 /(1008 \times 3058673),$ and $C_4 = 617 \times
1016064/(144 \times 3058673).$

\item \texttt{InputMasses:}
This structure is one of the members of the \texttt{InspiralTemplate} 
structure. 

\input{LALInputMassesH}

A user can specify the parameters of a binary using any two of the
following combination of masses
\begin{itemize}
\item \texttt{m1Andm2:} component masses
\item \texttt{totalMassAndEta:} total mass and symmetric mass ratio
\item \texttt{totalMassAndMu:} total mass and reduced mass
\item \texttt{t01:} unused; shouldn't be used.
\item \texttt{t02:} chirp times $\tau_0$ and $\tau_2$
\item \texttt{t03:} chirp times $\tau_0$ and $\tau_3$, and 
\item \texttt{t04:} chirp times $\tau_0$ and $\tau_4$
\end{itemize}


\item \texttt{InputMasses:}
Enum that tells which post-Newtonian order is being used.
\input{LALInspiralOrderH}
\begin{itemize}
\item \texttt{newtonian:} Newtonain order, flux and enrgy both to the lowest order.
\item \texttt{oneHalfPN:} same as before
\item \texttt{onePN:} Both energy and flux to order $O(v^2)$ beyond the Newtonian order.
\item \texttt{onePointFivePN:} Energy to order $O(v^2)$ and flux to order $O(v^3)$
\item \texttt{twoPN:} Both energy and flux to order $O(v^4)$
\item \texttt{twoPointFivePN:} Energy to order $O(v^4)$ and flux to order $O(v^5)$
\item \texttt{threePN:} Both energy and flux to order $O(v^6)$
\item \texttt{threePointFivePN:} Energy to order $O(v^6)$ and flux to order $O(v^7)$
\end{itemize}
In all cases, the gravitational wave phase (also frequency and time)
as an expansion of the gauge invariant parameter $v$ is given up to 
the order specified by flux.  Note that there are certain undetermined 
parameters at \texttt{threePN} and \texttt{threePointFivePN.} The waveform 
generation codes use a specific 
value of those parameters while generating the wave.


\item \texttt{Approximant:}
\input{LALInspiralApproximantH}
\begin{itemize}
\item \texttt{TaylorT1:} Time domain Taylor approximant in which 
	the energy and flux are both kept as Taylor expansions
	and a first order ordinary differential equation is solved
	for the GW phase as a function of $t.$
\item \texttt{TaylorT2:} Time domain Taylor approximant in which 
	the phase evolution $\varphi(t)$ is obtained by iteratively 
	solving post-Newtonian expansions $\varphi(v)$ and $t(v).$
\item \texttt{TaylorT3:} Time domain Taylor approximant in which 
phase is explicitly given as a function of time.
\item \texttt{TaylorF1:} The standard stationary phase approximation.
\item \texttt{TaylorF2:} The stationary phase approximation that
correctly represents, in the Fourier domain, the waveform given 
by \texttt{TaylorT1} approximant (see Damour, Iyer, Sathyaprakash,
		Phys. Rev. D . 63, 44023 (2001) for details).
\item \texttt{PadeT1:} Time-domain P-approximant.
\item \texttt{PadeF1:} Frequency-domain P-approximant (not yet implemented).
\item \texttt{EOB:} Effective one-body waveform 
\item \texttt{DJS:} Effective one-body waveform to 3.5 PN order 
\item \texttt{INSPA:} Improved stationary phase approximation (not implemented yet)
\item \texttt{IRSPA:} Improved relativistic stationary phase approximation (not implemented yet)
\end{itemize}


\item \texttt{InspiralTemplate:}
The inspiral waveform parameter structure containing information about the
	waveform to be generated.
\input{LALInspiralTemplateH}

\begin{itemize}
  \item \texttt { ieta:} parameter that tells whether the symmetric mass ratio $\eta$
	  should be set to zero in the PN expansions of GW flux and binding energy.
	  If \texttt{ieta=0} $\eta$ will be set to zero, otherwise the appropriate
	  value of $\eta$ from the given parameters will be used.
	  
  \item \texttt { level:} (introduced by Duncan Brown?)
  \item \texttt { *segmentIdVec:} (introduced by Duncan Brown?)
  \item \texttt { number:} (introduced by Duncan Brown?)
  \item \texttt { nStartPad:} Number of leading elements to be set to zero (input).
  \item \texttt { nEndPad:} Number of trailing bins to be set to zero, the 
  resulting waveform will have at least this many bins zero at the end, probably
  more since we always deal with an integer power of 2 array (input). 
  \item \texttt { mass1:}  Mass of the primary in solar mass (input/output).
  \item \texttt { mass2:}  Mass of the secondary in solar mass 
  (\texttt{mass1} need not be larger than \texttt{mass2} (input/output).
  \item \texttt { spin1[3]:} Spin vector of the primary (currently not in use)
  \item \texttt { spin2[3]:} Spin vector of the secondary (currently not in use)
  \item \texttt { inclination:} Inclination of the orbit  (currently not in use)
  \item \texttt { eccentricity:} initial eccentricity of the orbit  (currently not in use)
  \item \texttt { totalMass:} total mass of the binary $m=m_1+m_2$ in solar mass (input/output).
  \item \texttt { eta:} symmetric mass ratio $\eta=m_1m_2/m^2.$ (input/output). 
  \item \texttt { chirpMass:} chirp mass of the binary $=\eta^{3/5} m$ in solar mass (output).
  \item \texttt { t0:} Newtonain chirp time in seconds (input/output).
  \item \texttt { t2:} first post-Newtonian chirp time in seconds (input/output).
  \item \texttt { t3:} 1.5 post-Newtonian chirp time in seconds (input/output).
  \item \texttt { t4:} second post-Newtonian chirp time in seconds (output).
  \item \texttt { t5:} 2.5 post-Newtonian chirp time in seconds (output).
  \item \texttt { tC:} total chirp time seconds (output).
  \item \texttt { mu:} reduced mass (in solar mass) (input/output)
  \item \texttt { fLower:} lower frequency cutoff of the detector in Hz (input)
  \item \texttt { fCutoff:} upper frequency cutoff in Hz to be used in generating the waveform.
  If the last stable orbit frequency is smaller than the upper cutoff it will be used
  in terminating the waveform instead of fCutoff (input).
  \item \texttt { tSampling:} Sampling rate in Hz (input)
  \item \texttt { startPhase:} starting phase of the waveform in radians (input)
  \item \texttt { startTime:} starting time of the waveform (in sec); if different from zero, the
  waveform will start with an instantaneous frequency different from fLower and reach 
  fLower at time (approximately) zero (input, not used in Stationary phase approximation)
  \item \texttt { signalAmplitude:} dimensionless amplitude of the signal (input, currently unused.)
  \item \texttt { rInitial:} initial radial separation of the two, in units of total mass
  bodies (used only in EOB waveforms) (output)
  \item \texttt { vInitial:} initial velocity parameter, in units of the speed of light (output)
  \item \texttt { rFinal:} final 'separation' between the bodies, in units of total mass (output)
  \item \texttt { vFinal:} final velocity parameter, in units of the speed of light (output)
  \item \texttt { fFinal:} final frequency reached, in units of Hz (output)
  \item \texttt { rLightRing:} radial coordinate at the light ring, in units of total mass (output)
  \item \texttt { OmegaS:} The 3PN (unknown) parameter; calculated to be equal to zero
  by Damour, Jaranowski and Schaffer (input).
  \item \texttt { Theta:} The 3PN unknown flux parameter; likely to be around unity;
  most waveform generation routines take theta to be zero. Robustness of the EOB waveform
	  has been demonstrated for $-2 < $ \texttt{Theta} $< 2.$ (input)
  \item \texttt { massChoice:} The pair of (mass) parameters given (see structure
		  defining this member for more details) (input).
  \item \texttt { order:} Post-Newtonain order to be used in generating the wave (input).
  \item \texttt { approximant:} Post-Newtonain approximant to be used in generating the wave (input).
  \item \texttt { tagInspiralTemplate *next:} Linked list to the next coarse bank template 
  (currently not filled by inspiral or bank codes)
  \item \texttt { tagInspiralTemplate *fine:} Linked list to the next fine bank template
  (currently not filled by inspiral or bank codes)
\end{itemize}

\item \texttt{InspiralACSTParams:}
\input{LALInspiralACSTParamsH} 
This is a structure needed to generate solve the differential equation
	giving the evolution of the orbital angular momentum and the
	spin angular momenta in the case of spinning black hole binaries.
\begin{itemize}
  \item	\texttt {v:} parameter of 'integration': v=sqrt(M/r) 
  \item {magS1:} The constant magnitudes of the primary.
  \item {magS2:} The constant magnitudes of the secondary.
  \item {NCap[3]:} Source direction (unit vector) in detector coordinate system.
  \item {M:} Total mass of the binary (in seconds).
  \item {fourM1Plus:} = $(4 m_1+3 m_2)/(2 m_1 M^3)$ (all masses expressed in seconds).
  \item {fourM2Plus:} = $(4 m_2+3 m_1)/(2 m_2 M^3)$ (all masses expressed in seconds).
  \item {oneBy2Mcube:} = $1/(2 M^3)$
  \item {threeBy2Mcube:}  = $3/(2 M^3)$
  \item {thirtytwoBy5etc:}=  $(32/5) \eta^2 M$
\end{itemize}

\item \texttt{InspiralToffInput:}
\input{LALInspiralToffInputH} 
This is a structure needed by the inner workings of the inspiral wave generation code.


\item \texttt{expnCoeffs:}
This structure contains various post-Newtonian and P-approximant expansion 
	coefficients; the meanings of the coefficients is indicated as comments
	before each list.
\input{LALExpnCoeffsH} 

\item {Energy, flux, phase, time and frequency functions:} The following
	functions are generic function definitions that will be used in 
	template generation. The function \texttt{LALInspiralChooseModel,}
        which is called by wave generation interface code, points these 
	functions to the appropriate specific functions depending on the
	choices made by the user.
\input{LALEnergyAndFluxFunctionsH}

\item \texttt{expnFunc:} Structure to hold all the pointer to generic
       	functions defined above.
\input{LALexpnFuncH}

\item {Misc structures:} 
The rest of the structures below define other structures needed by wave
generation codes. (Will produce a documentation some day so that future code
maintainers can understand what precisely has been done.)

\end{enumerate}
</lalLaTeX>  */

/* <lalVerbatim file="LALEtaTau0Tau2InH">  */
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

/* <lalVerbatim file="LALEtaTau0Tau2InH">  */
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


/* <lalVerbatim file="LALInputMassesH">  */
typedef enum {
  m1Andm2,
  totalMassAndEta,
  totalMassAndMu,
  t01,
  t02,
  t03,
  t04
 } InputMasses;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{InputMasses}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralOrderH">  */
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



/* <lalVerbatim file="LALInspiralApproximantH">  */
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


/* <lalVerbatim file="LALInspiralTemplateH">  */
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
  REAL8 sourceTheta;
  REAL8 sourcePhi;
  REAL8 orbitTheta0;
  REAL8 orbitPhi0;
  REAL8 distance; /* in seconds */
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
  REAL8 fFinal;
  REAL8 rLightRing;
  REAL8 OmegaS;
  REAL8 Theta;
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

/* <lalVerbatim file="LALInspiralToffInputH">  */
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

/* <lalVerbatim file="LALInspiralACSTParamsH">  */
typedef struct
tagInspiralACSTParams
{
	REAL8 v;               
	REAL8 magS1;           
	REAL8 magS2;
	REAL8 NCap[3];         
	REAL8 M;               
	REAL8 fourM1Plus;      
	REAL8 fourM2Plus;      
	REAL8 oneBy2Mcube;     
	REAL8 threeBy2Mcube;   
	REAL8 thirtytwoBy5etc;  
} InspiralACSTParams;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{InspiralACSTParamsH}
</lalLaTeX>  */


/* <lalVerbatim file="LALExpnCoeffsH">  */
typedef struct {
   int ieta;
   /* coefficients in the Taylor expansion of new energy function*/
   double eTaN, eTa1, eTa2, eTa3;
   /* coefficients in the Pade expression of new energy function*/
   double ePaN, ePa1, ePa2, ePa3;
   /* coefficients in the Taylor expansion of usual energy function*/
   double ETaN, ETa1, ETa2, ETa3;
   /* coefficients in the Taylor expansion of the derivative of the 
    usual energy function*/
   double dETaN, dETa1, dETa2, dETa3;

   /* Taylor expansion coefficients of energy flux*/
   double FTaN, FTa1, FTa2, FTa3, FTa4, FTa5, FTa6, FTa7, FTl6;
   /* Taylor expansion coefficients of factored flux*/
   double fTaN, fTa1, fTa2, fTa3, fTa4, fTa5, fTa6, fTa7;
   /* Coefficients of the corresponding P-approximant*/
   double fPaN, fPa1, fPa2, fPa3, fPa4, fPa5, fPa6, fPa7;

   /* Taylor expansion coefficents in t(v)*/
   double tvaN, tva2, tva3, tva4, tva5, tva6, tva7, tvl6; 
   /* Taylor expansion coefficents in phi(v)*/
   double pvaN, pva2, pva3, pva4, pva5, pva6, pva7, pvl6; 
   /* Taylor expansion coefficents in phi(t)*/
   double ptaN, pta2, pta3, pta4, pta5, pta6, pta7, ptl6; 
   /* Taylor expansion coefficents in f(t)*/
   double ftaN, fta2, fta3, fta4, fta5, fta6, fta7, ftl6; 
   /* Taylor expansion coefficents in psi(f) in the Fourier phase*/
   double pfaN, pfa2, pfa3, pfa4, pfa5, pfa6, pfa7, pfl5; 

   /* sampling rate and interval*/
   double samplingrate, samplinginterval;
   /* symmetric mass ratio, total mass, component masses*/
   double eta, totalmass, m1, m2;
   /* unknown 3PN parameters, euler constant*/
   double lambda, theta, EulerC, omegaS;

   /* initial and final values of frequency, time, velocity; lso
    values of velocity and frequency; final phase.*/
   double f0, fn, t0, tn, v0, vn, vf, vlso, flso, phiC;

   /* last stable orbit and pole defined by various Taylor and P-approximants*/
   double vlsoT0, vlsoT2, vlsoT4, vlsoT6;
   double vlsoP0, vlsoP2, vlsoP4, vlsoP6;
   double vpoleP4, vpoleP6;
} expnCoeffs;

/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{expnCoeffs}
</lalLaTeX>  */

/* <lalVerbatim file="LALEnergyAndFluxFunctionsH"> */

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

/* </lalVerbatim> */

/* <lalVerbatim file="LALexpnFuncH">  */
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
   UINT4 *n,
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
\newpage\input{LALInspiralSpinModulatedWaveC}
</lalLaTeX>  */

void 
LALInspiralSpinModulatedWave(
		LALStatus        *status, 
		REAL4Vector      *signal, 
		InspiralTemplate *in);
/*  </lalLaTeX> */

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
/*  <lalLaTeX>
\newpage\input{LALStationaryPhaseApprox1C}
</lalLaTeX>  */
void 
LALInspiralStationaryPhaseApprox1 (
   LALStatus *status,
   REAL4Vector *signal,
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALEOBWaveformC}
</lalLaTeX>  */

void LALEOBWaveform(
	LALStatus *status,
	REAL4Vector *signal,
	InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALEOBWaveformTemplatesC}
</lalLaTeX>  */

void LALEOBWaveformTemplates(
	LALStatus *status,
	REAL4Vector *signal1,
	REAL4Vector *signal2,
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
\newpage\input{LALStationaryPhaseApprox2C}
</lalLaTeX>  */
void 
LALInspiralStationaryPhaseApprox2 (
   LALStatus *status,
   REAL4Vector *signal,
   InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralTestC}
</lalLaTeX>  */

#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRAL_H */
