/*
*  Copyright (C) 2007 Stas Babak, David Churches, Drew Keppel, Duncan Brown, Jolien Creighton, David McKechan, Patrick Brady, Peter Shawhan, Reinhard Prix, B.S. Sathyaprakash, Anand Sengupta, Craig Robinson , Sean Seader, Thomas Cokelaer
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/* <lalVerbatim file="LALInspiralHV">

Author: Churches, D. K ,  B. S. Sathyaprakash,  T. Cokelaer.
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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <lal/LALGSL.h>

# include <lal/LALStdlib.h>
# include <lal/LALConstants.h>
# include <lal/SimulateCoherentGW.h>
# include <lal/GeneratePPNInspiral.h>
# include <lal/LIGOMetadataTables.h>

 
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
#define  ninty4by3etc  18.687902694437592603 /* (94/3 -41/31*pi*pi) */

/*  <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>  */

/* <lalErrTable> */
#define LALINSPIRALH_ENULL           1
#define LALINSPIRALH_EMEM            2
#define LALINSPIRALH_EDIV0           3
#define LALINSPIRALH_ESIZE           4
#define LALINSPIRALH_ECHOICE         5
#define LALINSPIRALH_EORDER          6 
#define LALINSPIRALH_EAPPROXIMANT    7 
#define LALINSPIRALH_EPSI0           8 
#define LALINSPIRALH_EPSI3           9 
#define LALINSPIRALH_EALPHA         10
#define LALINSPIRALH_EFCUTOFF       11
#define LALINSPIRALH_ENOWAVEFORM    12
#define LALINSPIRALH_ESTOPPED       13
#define LALINSPIRALH_EROOTINIT      14
#define LALINSPIRALH_EFLOWER        15  
#define LALINSPIRALH_EVECTOR        16
#define LALINSPIRALH_EFLOWERINJ     17
#define LALINSPIRALH_EORDERMISSING  18
#define LALINSPIRALH_EBPERR         19


#define LALINSPIRALH_MSGENULL         "Arguments contained an unexpected null pointer"
#define LALINSPIRALH_MSGEMEM          "Memory allocation error"
#define LALINSPIRALH_MSGEDIV0         "Division by zero"
#define LALINSPIRALH_MSGESIZE         "Invalid input range"
#define LALINSPIRALH_MSGECHOICE       "Invalid choice for an input parameter"
#define LALINSPIRALH_MSGEORDER        "unknown order specified"
#define LALINSPIRALH_MSGEAPPROXIMANT  "Invalid model"
#define LALINSPIRALH_MSGEPSI0         "psi0 must be > 0"
#define LALINSPIRALH_MSGEPSI3         "psi3 must be < 0"
#define LALINSPIRALH_MSGEALPHA        "alpha must be defined positive"
#define LALINSPIRALH_MSGEFCUTOFF      "fcutoff must be defined and > 0"
#define LALINSPIRALH_MSGENOWAVEFORM   "No Waveform generated"
#define LALINSPIRALH_MSGESTOPPED      "Waveform generation stopped"
#define LALINSPIRALH_MSGEROOTINIT     "Can't find good bracket for BisectionFindRoot"
#define LALINSPIRALH_MSGEFLOWER       "fLower too low in comparison to flso"
#define LALINSPIRALH_MSGEVECTOR       "Attempting to write beyond the end of vector"
#define LALINSPIRALH_MSGEFLOWERINJ    "flower for the injection must be greater than zero"
#define LALINSPIRALH_MSGEORDERMISSING "The PN order requested is not implemented for this approximant"
#define LALINSPIRALH_MSGEBPERR        "Error in band passing signal"

/** ---------------------------------------------------------------------  </lalErrTable> */



/* <lalLaTeX>

\section*{Structures}

\begin{enumerate} 

\item \texttt{Order:}
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


\item \texttt{Approximant:} Enum that specifies the PN approximant to
be used in computing the waveform.
\input{LALInspiralApproximantH}
\begin{itemize}
\item \texttt{TaylorT1:} Time domain Taylor approximant in which 
	the energy and flux are both kept as Taylor expansions
	and a first order ordinary differential equation is solved
	for the GW phase as a function of $t.$ Outputs a time-domain wave.
\item \texttt{TaylorT2:} Time domain Taylor approximant in which 
	the phase evolution $\varphi(t)$ is obtained by iteratively 
	solving post-Newtonian expansions $\varphi(v)$ and $t(v).$ Outputs a time-domain wave.
\item \texttt{TaylorT3:} Time domain Taylor approximant in which 
phase is explicitly given as a function of time. Outputs a time-domain wave.
\item \texttt{TaylorF1:} The stationary phase approximation that
correctly represents, in the Fourier domain, the waveform given 
by \texttt{TaylorT1} approximant (see Ref. \cite{dis2} for details). Outputs a frequency-domain wave.
\item \texttt{TaylorF2:} The standard stationary phase approximation. Outputs a frequency-domain wave.
\item \texttt{PadeT1:} Time-domain P-approximant. Outputs a time-domain wave.
\item \texttt{PadeF1:} Frequency-domain P-approximant (not yet implemented).
\item \texttt{EOB:} Effective one-body waveform  Outputs a time-domain wave.
\item \texttt{BCV:} Detection template family of Buonanno, Chen and 
                    Vallisneri \cite{BCV03}. Outputs a frequency-domain wave.
\item \texttt{BCVSpin:} Detection template family of Buonanno, Chen and 
                    Vallisneri including  spin effects\cite{BCV03b}. Outputs a frequency-domain wave.
\item \texttt{SpinTaylorT3} Spinning case T3 models
\item \texttt{SpinTaylor} Spinning case PN models (should replace SpinTaylorT3 in the future)
\item \texttt{FindChirpSP} The stationary phase templates implemented by FindChirpSPTemplate in the findchirp package (equivalent to TaylorF2 at twoPN order).
\item \texttt{GeneratePPN} The time domain templates generated by LALGeneratePPNInspiral() in the inject package (equivalent to TaylorT3 at twoPN order).
\item \texttt{FrameFile} The waveform contains arbitrary data read from a frame file.


\end{itemize}
\input{LALInputMassesH}
\texttt{InputMasses:}
This structure is one of the members of the \texttt{InspiralTemplate} 
structure. 


A user can specify the parameters of a binary using any of the
following combination of {\it masses:}
\begin{itemize}
\item \texttt{m1Andm2:} component masses
\item \texttt{totalMassAndEta:} total mass and symmetric mass ratio
\item \texttt{totalMassUAndEta:} total mass and eta but uniform distribution in totalMass
\item \texttt{totalMassAndMu:} total mass and reduced mass
\item \texttt{t01:} unused; shouldn't be used.
\item \texttt{t02:} chirptimes $\tau_0$ and $\tau_2$
\item \texttt{t03:} chirptimes $\tau_0$ and $\tau_3$, and 
\item \texttt{t04:} chirptimes $\tau_0$ and $\tau_4$
\item \texttt{psi0Andpsi3:} BCV parameters $\psi_0$ and $\psi_3$
\end{itemize}
The LALRandomInspiralSignal uses that structure as an input. Since the injected
waveform are not necessarely wanted to be random, we also provide the following 
options
\begin{itemize}
\item \texttt{bhns:} One of the mass is a Neutron star and the other a black 
hole. (m1 $\in$ [minMass-3] and m2 $\in$ [3-maxMass]).
\item \texttt{fixedMasses:} The two masses are given by the input parameter structure.
\item \texttt{fixedPsi:} The two psi values are given by the input parameter structure.
\item \texttt{fixedTau:} The two tau values are given by the input parameter structure.
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
	  
  \item \texttt { level:} Flag used in heirarical serached to indicate if this is a coarse or a fine template
  \item \texttt { *segmentIdVec:} Vector of segment that have been filtered against this template needed for the LDAS implementation of the inspiral search.
  \item \texttt { number:} Unique ID number for this template needed for the LDAS implementation of the inspiral search.
  \item \texttt { minMatch:} The minimal match specified by the user when the bank that contains this template was created.
  \item \texttt { nStartPad:} Number of leading elements in the signal generation to be set to zero (input). If template is requested, that value must be set to zero. In the injection routines related to inject package, that nStartPad is set to zero. However, for injection performed using the inspiral package, that value can be set to non zero. 
  \item \texttt { nEndPad:} Number of trailing bins to be set to zero, the 
  resulting waveform will have at least this many bins zero at the end, probably
  more since we always deal with an integer power of 2 array (input). 
  \item \texttt { mass1:}  Mass of the primary in solar mass (input/output).
  \item \texttt { mass2:}  Mass of the secondary in solar mass 
  (\texttt{mass1} need not be larger than \texttt{mass2} (input/output).
  \item \texttt { spin1[3]:} Spin vector of the primary (currently not in use)
  \item \texttt { spin2[3]:} Spin vector of the secondary (currently not in use)
  \item \texttt { sourceTheta:} Co-latitute in the direction to the source.
  \item \texttt { sourcePhi:} Azimuth angle in the direction to the source.
  \item \texttt { orbitTheta0:} Initial co-latitute of the orbit.
  \item \texttt { orbitPhi0:} Initial azimuth angle of the orbit.
  \item \texttt { inclination:} Inclination of the orbit  (currently not in use)
  \item \texttt { distance:} Distance to the binary in seconds
  \item \texttt { psi0:} BCV parameter $\psi_0.$
  \item \texttt { psi3:} BCV parameter $\psi_3.$
  \item \texttt { alpha:} BCV amplitude correction factor $\alpha f_{\rm cut}^{2/3}$
  \item \texttt { eccentricity:} initial eccentricity of the orbit  (currently not in use)
  \item \texttt { totalMass:} total mass of the binary $m=m_1+m_2$ in solar mass (input/output).
  \item \texttt { eta:} symmetric mass ratio $\eta=m_1m_2/m^2.$ (input/output). 
  \item \texttt { chirpMass:} chirp mass of the binary $=\eta^{3/5} m$ in solar mass (output).
  \item \texttt { t0:} Newtonain chirp time in seconds (input/output).
  \item \texttt { t2:} first post-Newtonian chirp time in seconds (input/output).
  \item \texttt { t3:} 1.5 post-Newtonian chirp time in seconds (input/output).
  \item \texttt { t4:} second post-Newtonian chirp time in seconds (output).
  \item \texttt { t5:} 2.5 post-Newtonian chirp time in seconds (output).
  \item \texttt { t6:} third post-Newtonian chirp time in seconds (output).
  \item \texttt { t7:} 3.5 post-Newtonian chirp time in seconds (output).
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
This structure is needed to solve the differential equation
	giving the evolution of the orbital angular momentum and the
	spin angular momenta in the case of spinning black hole binaries.
\input{LALInspiralACSTParamsH} 
\begin{itemize}
  \item	\texttt {v:} parameter of 'integration': v=sqrt(M/r) 
  \item {\tt magS1:} The constant spin magnitude of the primary.
  \item {\tt magS2:} The constant spin magnitude of the secondary.
  \item {\tt NCap[3]:} Source direction (unit vector) in detector coordinate system.
  \item {\tt spin1[3]:} Spin of the larger body.
  \item {\tt M:} Total mass of the binary (in seconds).
  \item {\tt fourM1Plus:} = $(4 m_1+3 m_2)/(2 m_1 M^3)$ (all masses expressed in seconds).
  \item {\tt fourM2Plus:} = $(4 m_2+3 m_1)/(2 m_2 M^3)$ (all masses expressed in seconds).
  \item {\tt oneBy2Mcube:} = $1/(2 M^3)$
  \item {\tt threeBy2Mcube:}  = $3/(2 M^3)$
  \item {\tt thirtytwoBy5etc:}=  $(32/5) \eta^2 M$
\end{itemize}

\item \texttt{EtaTau02In, EtaTau04In:}
These are the input structures needed to solve for the mass
ratio $\eta$ given the chirptimes $(\tau_0,\, \tau_2)$ or
$(\tau_0, \, \tau_4).$

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

\item \texttt{InspiralToffInput:}
This is a structure needed by the inner workings of the inspiral wave generation code.
\input{LALInspiralToffInputH} 


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

\item \texttt{expnFunc:} Structure to hold the pointers to the generic
       	functions defined above.
\input{LALexpnFuncH}

\item \texttt{TofVIn {\rm and} TofVIntegrandIn:} Structures needed to 
	compute the time elapsed
	from/to the starting epoch of the waveform when the velocity 
	parameter was $v_0,$ to/from the current epoch when velocity 
	parameter is $v.$
\input{LALInspiralTofVH}

\item {\tt InspiralPhaseIn {\rm and} PhiofVIntegrandIn:} Structures used 
	to compute the phase of the signal from the `beginning', when the
	veolcity parameter is $v_0,$ to a time when the velocity parameter 
	has evolved to a user input value $v.$
\input{LALInspiralPhaseH}

\item {\tt InspiralDerivativesIn:} Structure used as an input to compute
	the derivatives needed in solving the phasing formula when the
	{\tt approximant} is {\tt TaylorT1, TaylorP1} or {\tt EOB.}
\input{LALInspiralDerivativesH}

\item {\tt rk4GSLIntegrator:} Structure containing steps and controls
for the GSL Runge-Kutta solver
\input{LALInspiralGSLRungeKuttaH}

\item {\tt rk4In:} Structure used as an input to Runge-Kutta solver.
\input{LALInspiralRungeKuttaH}

\end{enumerate}
--------------------------------------------------------------------- </lalLaTeX>  */

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



/* <lalVerbatim file="LALInspiralOrderH">  */
typedef enum {
   newtonian,
   oneHalfPN,
   onePN,
   onePointFivePN,
   twoPN,
   twoPointFivePN,
   threePN,
   threePointFivePN,
   pseudoFourPN
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
   BCV,
   BCVSpin,
   SpinTaylorT3,
   SpinTaylor,
   FindChirpSP,
   FindChirpPTF,
   GeneratePPN,
   BCVC,
   FrameFile,
   AmpCorPPN,
   NumRel,
   Eccentricity,
   EOBNR,
   IMRPhenomA,
   TaylorEt,
   TaylorT4,
   TaylorN,
   NumApproximants
 } Approximant;
/* </lalVerbatim>  */


/* <lalVerbatim file="LALInputMassesH">  */
typedef enum {
   m1Andm2,
   totalMassAndEta,
   totalMassUAndEta, 
   totalMassAndMu,
   t01,
   t02,
   t03,
   t04,
   psi0Andpsi3,
   bhns, 
   fixedMasses, 
   fixedPsi, 
   fixedTau,
   massesAndSpin,
   minmaxTotalMass,
   spinOnly
 } InputMasses;
/* </lalVerbatim>  */



/* <lalLaTeX>
\idx[Type]{InputMasses}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralTemplateH">  */
typedef struct
tagInspiralTemplate
{ 
/*  Parameters needed to generate Taylor/Pade waveforms */
  Approximant approximant;
  Order order;
  Order ampOrder;
  REAL8 mass1; 
  REAL8 mass2;
  REAL8 fCutoff;
  REAL8 fLower;
  REAL8 tSampling;
  REAL8 distance; 
  REAL8 signalAmplitude;
  REAL8 startPhase;
  REAL8 startTime;
  INT4  ieta;
 
/* Additional parameters for EOB waveforms */

  REAL8 Theta;
  REAL8 Zeta2;

/* Parameters for BCV1 template */

  REAL8 alpha;
  REAL8 psi0;
  REAL8 psi3;

/* Additional parameters for BCV2 template */

  REAL8 beta;
  REAL8 alpha1;
  REAL8 alpha2;
  REAL8 alpha3;
  REAL8 alpha4;
  REAL8 alpha5;
  REAL8 alpha6;

/* Parameters for spinning BH waveform */

  REAL8 inclination;
  REAL8 orbitTheta0;
  REAL8 orbitPhi0;
  REAL8 spin1[3];
  REAL8 spin2[3];
  REAL8 sourceTheta;
  REAL8 sourcePhi;
  REAL8 polarisationAngle;

/* Spin parameters for the PTF template */
  REAL8 chi; /* dimensionless spin of black hole (i.e. mass1) */
  REAL8 kappa; /* cosine of angle between spin of mass1 and orb ang mom */

/* Parameters which are currently might be used */

  REAL8 eccentricity;

/* Paramters which are computed using LALInspiralParameterCalc */

  REAL8 chirpMass; 
  REAL8 eta;
  REAL8 totalMass; 
  REAL8 fFinal;
  REAL8 t0; 
  REAL8 t2; 
  REAL8 t3; 
  REAL8 t4; 
  REAL8 t5; 
  REAL8 t6; 
  REAL8 t7; 
  REAL8 tC; 
 
/* Note that tc and fFinal are computed during waveform generation!!! */
 
  REAL4 minMatch;
  REAL8 mu; 
  INT4  level;
  INT4  number;
  INT4  nStartPad;
  INT4  nEndPad;
  REAL8 OmegaS;
  REAL8 vFinal;
/*  REAL8 vInitial;
  REAL8 rFinal;
  REAL8 rInitial;
  REAL8 rLightRing;*/
  InputMasses massChoice;
  INT4Vector *segmentIdVec; 
  LIGOTimeGPS end_time;
  EventIDColumn *event_id;
  CHAR ifo[LIGOMETA_IFO_MAX];

  /* Gamma[] is a vector that stores the upper triangular part of the metric in
   * the space of parameters. For time domain searches, Gamma[0,...,5] stores
   * the following information :
   *    Gamma[0] -> (tc,tc) metric component
   *    Gamma[1] -> (tc,t0) metric component
   *    Gamma[2] -> (tc,t3) metric component
   *    Gamma[3] -> (t0,t0) metric component
   *    Gamma[4] -> (t0,t3) metric component
   *    Gamma[5] -> (t3,t3) metric component
   * For spinBCV searches, (in 4 dimensions) Gamma[0,...,9] would be required.
   */
  REAL4  Gamma[10];
  
  struct tagInspiralTemplate *next;
  struct tagInspiralTemplate *fine; 
} InspiralTemplate;
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
 } InspiralToffInput;
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
   REAL8 spin1[3];         
   REAL8 M;               
   REAL8 fourM1Plus;      
   REAL8 fourM2Plus;      
   REAL8 oneBy2Mcube;     
   REAL8 threeBy2Mcube;   
   REAL8 thirtytwoBy5etc;  
}  InspiralACSTParams;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{InspiralACSTParamsH}
</lalLaTeX>  */


/* <lalVerbatim file="LALExpnCoeffsH">  */
typedef struct 
tagexpnCoeffs {
   int ieta;
   /* coefficients in the Taylor expansion of new energy function*/
   REAL8 eTaN, eTa1, eTa2, eTa3;
   /* coefficients in the Pade expression of new energy function*/
   REAL8 ePaN, ePa1, ePa2, ePa3;
   /* coefficients in the Taylor expansion of usual energy function*/
   REAL8 ETaN, ETa1, ETa2, ETa3;
   /* coefficients in the Taylor expansion of the derivative of the 
    usual energy function*/
   REAL8 dETaN, dETa1, dETa2, dETa3;

   /* Taylor expansion coefficients of energy flux*/
   REAL8 FTaN, FTa1, FTa2, FTa3, FTa4, FTa5, FTa6, FTa7, FTa8, FTl6, FTl8;
   /* Taylor expansion coefficients of factored flux*/
   REAL8 fTaN, fTa1, fTa2, fTa3, fTa4, fTa5, fTa6, fTa7, fTa8;
   /* Coefficients of the corresponding P-approximant*/
   REAL8 fPaN, fPa1, fPa2, fPa3, fPa4, fPa5, fPa6, fPa7, fPa8;

   /* Taylor expansion coefficents in t(v)*/
   REAL8 tvaN, tva2, tva3, tva4, tva5, tva6, tva7, tvl6; 
   /* Taylor expansion coefficents in phi(v)*/
   REAL8 pvaN, pva2, pva3, pva4, pva5, pva6, pva7, pvl6; 
   /* Taylor expansion coefficents in phi(t)*/
   REAL8 ptaN, pta2, pta3, pta4, pta5, pta6, pta7, ptl6; 
   /* Taylor expansion coefficents in f(t)*/
   REAL8 ftaN, fta2, fta3, fta4, fta5, fta6, fta7, ftl6; 
   /* Taylor expansion coefficents in psi(f) in the Fourier phase*/
   REAL8 pfaN, pfa2, pfa3, pfa4, pfa5, pfa6, pfa7, pfl5, pfl6; 
   /* Taylor expansion for the spinning case */
   REAL8 ST[9], thetahat ;


   /* sampling rate and interval*/
   REAL8 samplingrate, samplinginterval;
   /* symmetric mass ratio, total mass, component masses*/
   REAL8 eta, totalmass, m1, m2;
   /* unknown 3PN parameters, euler constant*/
   REAL8 lambda, theta, EulerC, omegaS, zeta2;

   /* initial and final values of frequency, time, velocity; lso
    values of velocity and frequency; final phase.*/
   REAL8 f0, fn, t0, tn, v0, vn, vf, vlso, flso, phiC;

   /* last stable orbit and pole defined by various Taylor and P-approximants*/
   REAL8 vlsoT0, vlsoT2, vlsoT4, vlsoT6;
   REAL8 vlsoP0, vlsoP2, vlsoP4, vlsoP6;
   REAL8 vpoleP4, vpoleP6;
}  expnCoeffs;

/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{expnCoeffs}
</lalLaTeX>  */

/* <lalVerbatim file="LALEnergyAndFluxFunctionsH"> */

typedef REAL8 EnergyFunction(
   REAL8 v, 
   expnCoeffs *ak);

typedef REAL8 FluxFunction(
   REAL8 v, 
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
} expnFunc;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{expnFunc}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralTofVH">  */
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
} TofVIn;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{TofVIn}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralTofVH">  */
typedef struct
tagTofVIntegrandIn
{
   EnergyFunction *dEnergy;
   FluxFunction *flux;
   expnCoeffs *coeffs;
} TofVIntegrandIn;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{TofVIntegrandIn}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralDerivativesH">  */
typedef struct
tagInspiralDerivativesIn
{
   REAL8 totalmass;
   EnergyFunction *dEnergy;
   FluxFunction *flux;
   expnCoeffs *coeffs;
} InspiralDerivativesIn;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{InspiralDerivativesIn}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralRungeKuttaH">  */
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
} rk4In;
/* </lalVerbatim>  */

/* <lalLaTeX>
\index{\texttt{rk4In}}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralGSLRungeKuttaH"> */
typedef struct
tagrk4GSLIntegrator
{
   const gsl_odeiv_step_type *type;
   gsl_odeiv_step *step;
   gsl_odeiv_control *control;
   gsl_odeiv_evolve *evolve;
   REAL8 *y;
   rk4In *input;
} rk4GSLIntegrator;
/* </lalVerbatim> */

/* <lalLaTeX>
\index{\texttt{rk4GSLIntegrator}}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralPhaseH">  */
typedef struct
tagInspiralPhaseIn
{
   REAL8 v0;
   REAL8 phi0;
   EnergyFunction *dEnergy;
   FluxFunction *flux;
   expnCoeffs *coeffs;
} InspiralPhaseIn;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{InspiralPhaseIn}
</lalLaTeX>  */


/* <lalVerbatim file="LALInspiralPhaseH">  */
typedef struct
tagPhiofVIntegrandIn
{
   EnergyFunction *dEnergy;
   FluxFunction *flux;
   expnCoeffs *coeffs;
}  PhiofVIntegrandIn;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{PhiofVIntegrandIn}
</lalLaTeX>  */



/* <lalVerbatim file="LALInspiralInitH">  */
typedef struct
tagInspiralInit
{
  UINT4      nbins;
  expnCoeffs ak;
  expnFunc   func;

}  InspiralInit;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{InspiralInit}
</lalLaTeX>  */

/* <lalVerbatim file="LALInspiralApplyTaperH">  */
typedef enum
{
  INSPIRAL_TAPER_NONE,
  INSPIRAL_TAPER_START,
  INSPIRAL_TAPER_END,
  INSPIRAL_TAPER_STARTEND,
  INSPIRAL_TAPER_NUM_OPTS
}  InspiralApplyTaper;
/* </lalVerbatim>  */

/* <lalLaTeX>
\idx[Type]{InspiralInit}
</lalLaTeX>  */ 



/*  <lalLaTeX>
\vfill{\footnotesize\input{LALInspiralHV}}
</lalLaTeX>  */




/* Function prototypes */


/* --- HERE ARE SOME USEFUL PROTOTYPE FOR LENGTH, PARAMETER CALCULATION... --- */
/*  <lalLaTeX>
\newpage\input{LALInspiralParameterCalcC}
</lalLaTeX>  */

void LALInspiralParameterCalc (
     LALStatus *status,
     InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralAmplitudeC}
</lalLaTeX>  */

void LALInspiralRestrictedAmplitude(
     LALStatus *status,
     InspiralTemplate  *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralWaveLengthC}
</lalLaTeX>  */

void LALInspiralWaveLength (
     LALStatus *status,
     UINT4 *n,
     InspiralTemplate params);


/*  <lalLaTeX>
\newpage\input{LALInspiralChooseModelC}
</lalLaTeX>  */

void LALInspiralChooseModel(
     LALStatus *status,
     expnFunc *func,
     expnCoeffs *ak,
     InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralSetupC}
</lalLaTeX>  */

void LALInspiralSetup (
     LALStatus *status,
     expnCoeffs *ak,
     InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralInitC}
</lalLaTeX>  */
void 
LALInspiralInit(
	LALStatus        *status, 
	InspiralTemplate *params, 
	InspiralInit     *paramsInit);

/*  <lalLaTeX>
\newpage\input{LALInspiralWaveTaperC}
</lalLaTeX>  */

void LALInspiralWaveTaper(
     LALStatus    *status,
     REAL4Vector  *signal,
     UINT4       bookends
     );

int XLALInspiralWaveTaper(
                   REAL4Vector         *signal,
                   InspiralApplyTaper  bookends);

/* --- HERE ARE THE WAVEFORMS/MODELS PROTOTYPES --- */
/*  <lalLaTeX>
\newpage\input{LALInspiralWaveC}
</lalLaTeX>  */

void LALInspiralAmplitudeCorrectedWave(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALInspiralAmplitudeCorrectedWaveTemplates(
     LALStatus *status,
     REAL4Vector *filter1,
     REAL4Vector *filter2,
     InspiralTemplate *params);

void 
LALInspiralAmplitudeCorrectedWaveForInjection(
   LALStatus        *status,
   CoherentGW       *waveform,
   InspiralTemplate *params,
   PPNParamStruc  *ppnParams);

/*  <lalLaTeX>
\newpage\input{LALInspiralWaveC}
</lalLaTeX>  */

void LALInspiralWave(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALInspiralWaveTemplates(
     LALStatus *status,
     REAL4Vector *filter1,
     REAL4Vector *filter2,
     InspiralTemplate *params);

void 
LALInspiralWaveForInjection(
   LALStatus        *status,
   CoherentGW       *waveform,
   InspiralTemplate *params,
   PPNParamStruc  *ppnParams);

/*  <lalLaTeX>
\newpage\input{LALInspiralWave1C}
</lalLaTeX>  */

void LALInspiralWave1(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALInspiralWave1Templates(
     LALStatus *status,
     REAL4Vector *signalvec1,
     REAL4Vector *signalvec2,
     InspiralTemplate *params);

void LALInspiralWave1ForInjection(
     LALStatus        *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams			     
     );

void LALInspiralEccentricity(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALInspiralEccentricityTemplates(
     LALStatus *status,
     REAL4Vector *signalvec1,
     REAL4Vector *signalvec2,
     InspiralTemplate *params);

void LALInspiralEccentricityForInjection(
     LALStatus        *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams			     
     );


/*  <lalLaTeX>
\newpage\input{LALInspiralWave2C}
</lalLaTeX>  */

void LALInspiralWave2(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALInspiralWave2Templates (
     LALStatus *status,
     REAL4Vector *signalvec1, 
     REAL4Vector *signalvec2, 
     InspiralTemplate *params);

void LALInspiralWave2ForInjection(
     LALStatus        *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams			     
     );

/*  <lalLaTeX>
\newpage\input{LALInspiralWave3C}
</lalLaTeX>  */

void LALInspiralWave3 (
     LALStatus *status,
     REAL4Vector *signalvec, 
     InspiralTemplate *params);

void LALInspiralWave3Templates (
     LALStatus *status,
     REAL4Vector *signalvec1, 
     REAL4Vector *signalvec2, 
     InspiralTemplate *params);

void LALInspiralWave3ForInjection(
     LALStatus        *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams);

/*  <lalLaTeX>
\newpage\input{LALInspiralStationaryPhaseApprox1C}
</lalLaTeX>  */
void LALInspiralStationaryPhaseApprox1 (
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralStationaryPhaseApprox2C}
</lalLaTeX>  */
void LALInspiralStationaryPhaseApprox2 (
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALEOBWaveformC}
</lalLaTeX>  */

void LALEOBWaveform(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);

void LALEOBWaveformTemplates(
     LALStatus *status,
     REAL4Vector *signalvec1,
     REAL4Vector *signalvec2,
     InspiralTemplate *params);

void LALEOBWaveformForInjection(
     LALStatus *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams);

/*  <lalLaTeX>
\newpage\input{LALBCVWaveformC}
</lalLaTeX>  */

void LALBCVWaveform(
     LALStatus *status,
     REAL4Vector *signalvec, 
     InspiralTemplate *params);

void LALTaylorEtWaveform(
     LALStatus *status,
     REAL4Vector *signalvec, 
     InspiralTemplate *params);

void LALTaylorT4Waveform(
     LALStatus *status, 
     REAL4Vector *signalvec, 
     InspiralTemplate *params);

void LALBCVSpinWaveform(
     LALStatus *status,
     REAL4Vector *signalvec, 
     InspiralTemplate *params);

/* <lalLaTeX>
\newpage\input{LALInspiralTaylorNWaveformC}
</lalLaTeX> */

void LALTaylorNWaveform (
     LALStatus        *status,
     REAL4Vector      *signal,
     InspiralTemplate *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralSpinningBHBinaryC}
</lalLaTeX>  */

void LALInspiralSpinModulatedWave(
     LALStatus        *status, 
     REAL4Vector      *signalvec, 
     InspiralTemplate *in);


void LALInspiralSpinModulatedWaveForInjection(
     LALStatus *status,
     CoherentGW *waveform,
     InspiralTemplate *params,
     PPNParamStruc  *ppnParams
     );


/*  <lalLaTeX>
 \newpage\input{LALSTPNWaveformC}
</lalLaTeX>  */
void 
LALSTPNWaveformForInjection (
			    LALStatus        *status,
			    CoherentGW       *waveform,
			    InspiralTemplate *params,
			    PPNParamStruc  *ppnParams);
void 
LALSTPNWaveformEngine (
                LALStatus        *status,
                REAL4Vector      *signal1,
                REAL4Vector      *signal2,
                REAL4Vector      *a,
                REAL4Vector      *ff,
                REAL8Vector      *phi,
		REAL4Vector      *shift,
                UINT4            *countback,
                InspiralTemplate *params,
                InspiralInit     *paramsInit
                );
void 
LALSTPNWaveformTemplates (
   LALStatus        *status,
   REAL4Vector      *signal1,
   REAL4Vector      *signal2,
   InspiralTemplate *params
   ) ;

void LALSTPNWaveform(
     LALStatus *status,
     REAL4Vector *signalvec,
     InspiralTemplate *params);


/* Phenomenological waveform generation functions */

void LALBBHPhenWaveFreqDom ( LALStatus        *status,
			     REAL4Vector      *signal,
			     InspiralTemplate *params);

void LALBBHPhenWaveFreqDomTemplates( LALStatus        *status,
				     REAL4Vector      *signal1,
				     REAL4Vector      *signal2,
				     InspiralTemplate *params);

void LALBBHPhenWaveTimeDom ( LALStatus        *status,
			     REAL4Vector      *signal,
			     InspiralTemplate *template);

void LALBBHPhenWaveTimeDomTemplates( LALStatus        *status,
				     REAL4Vector      *signal1,
				     REAL4Vector      *signal2,
				     InspiralTemplate *params);

void LALBBHPhenTimeDomEngine( LALStatus        *status,
			      REAL4Vector      *signal1,
			      REAL4Vector      *signal2,
			      REAL4Vector      *h,
			      REAL4Vector      *a,
			      REAL4Vector      *f,
			      REAL8Vector      *phiOut,
			      InspiralTemplate *params);

void LALBBHPhenWaveTimeDomForInjection (LALStatus        *status,
					CoherentGW       *waveform,
					InspiralTemplate *params,
					PPNParamStruc    *ppnParams);


/* --- OTHER PROTOTYPES --- */

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
\newpage\input{LALInspiralVelocityC}
</lalLaTeX>  */

void LALInspiralVelocity (
     LALStatus *status,
     REAL8 *v,
     TofVIn *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralPhasing1C}
</lalLaTeX>  */

void LALInspiralPhasing1 (
     LALStatus *status,
     REAL8 *phase,
     REAL8 v,
     void *params);

/*  <lalLaTeX>
\newpage\input{LALInspiralPhiofVIntegrandC}
</lalLaTeX>  */

void LALInspiralPhiofVIntegrand (
     LALStatus *status,
     REAL8 *,
     REAL8,
     void *);


/*  <lalLaTeX>
\newpage\input{LALInspiralPhasing2C}
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

REAL4 LALInspiralHPlusPolarization( 
     REAL8 phase, 
     REAL8 v, 
     InspiralTemplate *params );

REAL4 LALInspiralHCrossPolarization( 
     REAL8 phase, 
     REAL8 v, 
     InspiralTemplate *params );

/*  <lalLaTeX>
\newpage\input{LALRungeKutta4C}
</lalLaTeX>  */

rk4GSLIntegrator * XLALRungeKutta4Init(
		   INT4 n,
                   rk4In *input);

void LALRungeKutta4(
     LALStatus *,
     REAL8Vector *,
     rk4GSLIntegrator *,
     void *);

void XLALRungeKutta4Free(
     rk4GSLIntegrator *integrator);

/* --- PARSING PROTOTYPE FOR INSPIRALTEMPLATE STRCUTURE --- */

/*  <lalLaTeX>
\newpage\input{LALInspiralParseParametersC}
</lalLaTeX>  */
void
LALInspiralITStructureParseParameters(
	LALStatus *status,
	UINT4 argc,
	CHAR **argv,
	InspiralTemplate *params);

void 
LALInspiralITStructureSetDefault(
	LALStatus *status, 
	InspiralTemplate *params);

void
LALInspiralITStructurePrint(
	LALStatus *status, 
	InspiralTemplate  params);

void
LALInspiralITStructureHelp(void);

/* --- TEST PROTOTYPES --- */

/*  <lalLaTeX>
\newpage\input{GenerateInspiralWaveformC}
</lalLaTeX>  */

INT4 XLALInspiralRingdownWave (
	REAL4Vector			*rdwave1,
	REAL4Vector			*rdwave2,
	InspiralTemplate		*params,
	REAL4VectorSequence		*inspwave1,
	REAL4VectorSequence		*inspwave2,
	COMPLEX8Vector			*modefreqs,
	UINT4				nmodes
	);
	
INT4 XLALGenerateWaveDerivatives (
	REAL4Vector		*dwave,
	REAL4Vector		*ddwave,
	REAL4Vector		*wave,
	InspiralTemplate	*params
	);
	
INT4 XLALGenerateQNMFreq(
	COMPLEX8Vector		*modefreqs,
	InspiralTemplate	*params,
	UINT4			l,
	UINT4			m,
	UINT4			nmodes
	);
	
INT4 XLALFinalMassSpin(
	REAL8			*finalMass,
	REAL8			*finalSpin,
	InspiralTemplate	*params
	);

INT4 XLALInspiralAttachRingdownWave (
        REAL4Vector 	 *Omega,
        REAL4Vector 	 *signal1,
        REAL4Vector  	 *signal2,
        InspiralTemplate *params);

int XLALInspiralGetApproximantString( CHAR        *output,
                                      UINT4       length,
                                      Approximant approx,
                                      Order       order
                                    );

int XLALBandPassInspiralTemplate(
        REAL4Sequence  *sequence,
        REAL4          fLow,
        REAL4          fHigh,
        REAL4          fSampling
        );


#ifdef  __cplusplus
}
#endif

#endif /* _LALINSPIRAL_H */
