/*********************** <lalVerbatim file="SimulateSBHV">
Author: Sukanta Bose
$Id$ 
*********************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{SimulateSB.h}}
\label{stochastic:s:SimulateSB.h}

Provides prototype and error code information for the modules needed
to simulate a whitened stochastic background signal in a pair of 
interferometric detectors, given the appropriate representations of the 
detector transfer function in each detector. 

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/SimulateSB.h>
\end{verbatim}

\noindent 

\subsection*{Error conditions}
\input{SimulateSBHE}

\subsection*{Structures}

*********************************************************** </lalLaTeX> */

#ifndef _SIMULATESB_H
#define _SIMULATESB_H


#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/Units.h>
#include <lal/RealFFT.h>

#ifdef  __cplusplus
extern "C" {
#endif

  NRCSID( SIMULATESBH, "$Id$" );
  
/***************** <lalErrTable file="SimulateSBHE"> */

#define SIMULATESBH_ENULLP          1
#define SIMULATESBH_ENONPOSLEN      2
#define SIMULATESBH_ENONPOSDELTAF   3
#define SIMULATESBH_ENONPOSDELTAT   4
#define SIMULATESBH_ENEGFMIN        5
#define SIMULATESBH_EMMTIME         6
#define SIMULATESBH_EMMHETERO       7
#define SIMULATESBH_EMMFMIN         8
#define SIMULATESBH_EMMDELTAF       9
#define SIMULATESBH_EMMLEN         10
#define SIMULATESBH_EOORFREF       11
#define SIMULATESBH_ENONPOSOMEGA   12
#define SIMULATESBH_EALOC          13
#define SIMULATESBH_ENONZEROHETERO 14
#define SIMULATESBH_EWRONGUNITS    15
#define SIMULATESBH_ECOMPTIME      16
#define SIMULATESBH_ENOTYETHETERO 255

#define SIMULATESBH_MSGENULLP         "Null pointer"
#define SIMULATESBH_MSGENONPOSLEN     "Negative or zero length for data member of time series" 
#define SIMULATESBH_MSGENONPOSDELTAF  "Negative or zero frequency spacing" 
#define SIMULATESBH_MSGENONPOSDELTAT  "Negative or zero time spacing" 
#define SIMULATESBH_MSGENEGFMIN       "Negative start frequency" 
#define SIMULATESBH_MSGEMMTIME        "Mismatch in epochs"
#define SIMULATESBH_MSGEMMHETERO      "Mismatch in heterodyning frequencies"
#define SIMULATESBH_MSGEMMFMIN        "Mismatch in start frequencies"
#define SIMULATESBH_MSGEMMDELTAF      "Mismatch in frequency spacings"
#define SIMULATESBH_MSGEMMLEN         "Mismatch in sequence lengths"
#define SIMULATESBH_MSGEOORFREF       "Out of range reference frequency"
#define SIMULATESBH_MSGENONPOSOMEGA   "Negative stochastic background strength"
#define SIMULATESBH_MSGEALOC         "Memory allocation error"
#define SIMULATESBH_MSGENONZEROHETERO "Non-zero heterodyning frequency specified for real time series"
#define SIMULATESBH_MSGEWRONGUNITS    "Inconsistent input units"
#define SIMULATESBH_MSGECOMPTIME      "Time domain data complex instead of real"
#define SIMULATESBH_MSGENOTYETHETERO  "Non-zero heterodyning frequency not yet implemented"

/************************************ </lalErrTable> */

  /*************************************************************
   *                                                           *
   *       Structures and prototypes associated with           *
   *             SimulateSB.c                  *
   *                                                           *
   *************************************************************/

/********************************************************** <lalLaTeX>

\subsubsection*{Structures associated with 
  \texttt{SimulateSB.c}
  (Sec.~\ref{stochastic:ss:SimulateSB.c})}

\subsubsection*{\texttt{struct SimulateSBOutput}}
\idx[Type]{SimulateSBOutput}

\noindent Contains the output data produced by 
\texttt{LALSimulateSB()}. It contains 
the simulated whitened stochastic background signal component of the output of 
an interferometric detector.
The fields are:

\begin{description}
\item[\texttt{REAL4TimeSeries *whitenedSimulatedSB1}]
The simulated whitened stochastic background signal component of the output of 
the first interferometric detector.
\item[\texttt{REAL4TimeSeries *whitenedSimulatedSB2}]
The simulated whitened stochastic background signal component of the output of 
the second interferometric detector.
\end{description}

*********************************************************** </lalLaTeX> */

  typedef struct tagSimulateSBOutput {
    REAL4TimeSeries    *whitenedSimulatedSB1;
    REAL4TimeSeries    *whitenedSimulatedSB2;
  } SimulateSBOutput;
  
  
  /*********************************************************** <lalLaTeX> 
							       
\subsubsection*{\texttt{struct SimulateSBInput}}
\idx[Type]{SimulateSBInput}
							       
\noindent Contains the input data needed by 
\texttt{LALSimulateSB()}
to calculate the whitened stochastic background signal in the output of 
an interferometric detector.
The fields are:

\begin{description}
\item[\texttt{REAL4FrequencySeries *overlapReductionFunction}]
The overlap reduction function $\gamma(f)$ describing the pair of detector
sites.
\item[\texttt{REAL4FrequencySeries *omegaGW}] The spectrum 
$\Omega_{\scriptstyle{\rm GW}}(f)$ of the stochastic gravitational-wave
background.
\item[\texttt{LALDetectorPair detectors}]
The site location information of the pair of detectors involved in the 
stochastic background search.
\item[\texttt{COMPLEX8FrequencySeries *whiteningFilter1}]
The frequency-domain response function $\tilde{R}_1(f)$ for the first detector.
\item[\texttt{COMPLEX8FrequencySeries *whiteningFilter2}]
The frequency-domain response function $\tilde{R}_2(f)$ for the second detector.
\end{description}

*********************************************************** </lalLaTeX> */

  typedef struct tagSimulateSBInput {
    REAL4FrequencySeries     *overlapReductionFunction;
    REAL4FrequencySeries     *omegaGW;
    const LALDetectorPair    *detectors; 
    COMPLEX8FrequencySeries  *whiteningFilter1;
    COMPLEX8FrequencySeries  *whiteningFilter2;
  } SimulateSBInput;
  
  
/*********************************************************** <lalLaTeX> 


\subsubsection*{\texttt{struct SimulateSBParams}}
\idx[Type]{SimulateSBParams}

\noindent Contains the parameters used by
\texttt{LALSimulateSB()} to determine the whitened stochastic background 
signal in the output of an interferometric detector.
The fields are:

\begin{description}
\item[\texttt{UINT4 length}]
The number of points in the output time series.

\item[\texttt{REAL8 deltaT}]
The temporal spacing of the output time series.

\item[\texttt{REAL8 f0}]
The start frequency of the frequency series 
$\Omega_{\scriptstyle{\rm GW}}(f)$.

\item[\texttt{REAL4 alpha}]
The exponent in the frequency power law of $\Omega_{\scriptstyle{\rm GW}}(f) 
= \Omega_{\scriptstyle{\rm R}} (f/f_{\scriptstyle{\rm R}})^\alpha$.

\item[\texttt{REAL4 fRef}]
The reference frequency in the frequency series 
$\Omega_{\scriptstyle{\rm GW}}(f)$.

\item[\texttt{REAL4 omegaRef}]
The value of $\Omega_{\scriptstyle{\rm GW}}(f)$ at the reference frequency
$f_{\scriptstyle{\rm R}}$.

\item[\texttt{COMPLEX8Vector *ccounts[2]}]
A pair of frequency series representing a whitened 
simulated stochastic background in a pair of detectors.

\item[\texttt{COMPLEX8Vector *ccountsTmp[2]}]
A pair of frequency series required to temporarily store data before 
computing \texttt{*ccounts[2]}.

\item[\texttt{REALFFTPlan *invPlan}]
The FFT plan needed to obtain time series representation of data in
\texttt{*ccounts[2]}.

\item[\texttt{REAL4Vector *gaussdevs}]
Vector for storing real random numbers with a normal distribution of
zero mean and unit variance.

\end{description}
*********************************************************** </lalLaTeX> */

  typedef struct tagSimulateSBParams {
    UINT4     length;   /* time length of output vector data samples */
    REAL8     deltaT;   /* time spacing */
    REAL8     f0;       /* start frequency */
    REAL4     alpha;    /* exponent in power law: omegaGW(f) = f^alpha */
    REAL4     fRef;     /* reference normalization frequency */
    REAL4     omegaRef; /* reference omega coefficient for normalization */
  }
  SimulateSBParams;

/********************************************************** <lalLaTeX>

\subsubsection*{\texttt{struct SimulateSBInitParams}}
\idx[Type]{SimulateSBInitParams}

\noindent Contains the parameters used by
\texttt{LALSimulateSB()} to create its output structure.
The fields are:

\begin{description}
\item[\texttt{UINT4 length}]
The number of points in the output time series.

\end{description}
*********************************************************** </lalLaTeX> */

  typedef struct tagSimulateSBInitParams {
    UINT4     length;   /* time length of output vector data samples */
  }
  SimulateSBInitParams;
  
  void
  LALCreateSimulateSBOutput (
			     LALStatus                     *status,
			     SimulateSBOutput             **output,
			     SimulateSBInitParams          *params
			     );
  
  void
  LALDestroySimulateSBOutput (
			      LALStatus                      *status,
			      SimulateSBOutput              **output
			      );
  
  void
  LALSimulateSBInit (
		     LALStatus                     *status,
		     SimulateSBParams             **output,
		     SimulateSBInitParams          *params
		     );
  
  void
  LALSimulateSBFinalize (
			 LALStatus                     *status,
			 SimulateSBParams             **output
			 );
  
  void
  LALCreateSimulateSBInput (
			    LALStatus                     *status,
			    SimulateSBInput              **output,
			    SimulateSBInitParams          *params
			    );
  
  void
  LALDestroySimulateSBInput (
			     LALStatus                      *status,
			     SimulateSBInput               **output
			     );
  
  void
  LALSimulateSB( LALStatus                  *status,
		 SimulateSBOutput          **output,
		 SimulateSBInput            *input,
		 SimulateSBParams           *params );
  
#ifdef  __cplusplus
}
#endif

#endif /* _SIMULATESB_H */
  
/********************************************************** <lalLaTeX>
							      
\vfill{\footnotesize\input{SimulateSBHV}}
							      
\newpage\input{SimulateSBC}
%\newpage\input{SimulateSBTestC}
*********************************************************** </lalLaTeX> */
