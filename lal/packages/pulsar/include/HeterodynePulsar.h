/********************************* <lalVerbatim file="HeterodynePulsarHV">
Author: Dupuis, R. J.
$Id$
********************************** </lalVerbatim> */

/********************************* <lalLaTeX>

\section{Header \texttt{HeterodynePulsar.h}}

Provides routines to heterodyne, average, and resample the data as required for time domain known pulsar search.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/HeterodynePulsar.h>
\end{verbatim}
The gravitational wave signal from a non-precessing pulsar at twice its frequency can be modeled as

\begin{equation}
h(t) = F_{+}(t;\psi)h_{0}(1 + \cos^{2}\iota)\cos 2 \Psi(t) + 2 F_{\times} h_{0}\cos \iota \sin 2 \Psi(t)
\end{equation}

where $F_{+,\times}$ are the amplitude responses of the detectors, $\psi$ is the polarization angle,
$\iota$ describes the inclination of the pulsar with respect to the line of sight, and $\Psi(t) =
\phi_{0} + \phi(t)$ describes the phase of the pulsar.

The phase $\Psi(t)$ of the pulsar is calculated as 
\begin{equation}
\Psi(t) = \phi_{0} +2\pi \left(f_{0}(T - T_{0}) + 
\frac{1}{2}\dot{f_{0}} (T - T_{0})^{2} + 
\frac{1}{6}\ddot{f_{0}}(T - T_{0})^{3}\right)
\end{equation}
where 
\begin{equation}
T = t + \delta t= t + \frac{\vec{r} \cdot \vec{n}}{c}  + \Delta_{E\odot}
\end{equation}
where {\emph T} is the time in a frame inertial with respect to the pulsar and $\phi_{0}$ is the phase of at time $T_{0}$.
 The time difference $\delta t$ due to the motion of the earth in the solar system is calculated using
 \texttt{LALBarycenter()}.

The function \texttt{LALCoarseHeterodyne()} heterodynes, averages, and resamples the data at a fixed frequency
near the signal. 

Let the calibrated data from the inteferometer be $d(t) = h(t) + n(t)$ where $n(t)$ is the noise.  The first step is to mix the
time series with $e^{-2\pi if_{h}}$ where $f_{h}$ is a fixed frequency near the signal.

\begin{equation}
V_{h}(t) = d(t)e^{-2\pi if_{h}}
\end{equation}

The function \texttt{LALFineHeterodyneToPulsar()} applies a second heterodyne to the data which removes the spindown
and the Doppler shifts.

After applying \texttt{LALCoarseHeterodyne()} and \texttt{LALFineHeterodyneToPulsar()} to the data we have

\begin{equation}
d(t)' = n(t)' + F_{+}(t;\psi)h_{0}(1 + \cos^{2}\iota) e^{i 2\phi_{0}} + 2 F_{\times} h_{0}\cos \iota e^{i 2\phi_{0}}.
\end{equation}


More documentation soon.

\begin{verbatim}
void
LALCoarseHeterodyne(        LALStatus                   *status,
                            CoarseHeterodyneOutput      *output,
                            CoarseHeterodyneInput       *input,
                            CoarseHeterodyneParams      *params );
		     
void
LALFineHeterodyneToPulsar(  LALStatus                   *status,
                            FineHeterodyneOutput        *output,
                            FineHeterodyneInput         *input,
                            FineHeterodyneParams        *params );	
\end{verbatim}

\subsection*{Error conditions}
\input{HeterodynePulsarHE}

\subsection*{Types}

\subsubsection*{Structure \texttt{CoarseHeterodyneInput}}
\idx[Type]{CoarseHeterodyneInput}

\noindent This structure stores the original calibrated gw data.
 
\begin{description}
\item[\texttt{REAL4TimeSeries V}] calibrated strain data from detector

\item[\texttt{REAL4 f0}] heterodyning base frequency

\end{description}

\subsubsection*{Structure \texttt{CoarseHeterodyneOutput}}
\idx[Type]{CoarseHeterodyneOutput}

\noindent This structure stores the output of the heterodyned data.

\begin{description}
\item[\texttt{COMPLEX8TimeSeries Vh}] heterodyned data

\item[\texttt{COMPLEX16 varh}] variance of Vh

\item[\texttt{REAL4 phase}] phase of the reference signal f0 at t0

\item[\texttt{COMPLEX16 avg}] average of Vh

\item[\texttt{COMPLEX16 kurt}] kurtosis of Vh

\item[\texttt{COMPLEX16 skew}] skewness of Vh

\item[\texttt{COMPLEX16 covar}] first term of covariance of Vh
\end{description}

\subsubsection*{Structure \texttt{CoarseHeterodyneParams}}
\idx[Type]{CoarseHeterodyneParams}

\noindent This structure stores parameters for the coarse heterodyne.

\begin{description}
\item[\texttt{UINT4 boxM}] first decimation factor (and order of boxcar)

\item[\texttt{REAL4IIRFilter *iirFilter1Re}] first IIR filter to be applied to real part of complex heterodyned data

\item[\texttt{REAL4IIRFilter *iirFilter1Im}] first IIR filter to be applied to imaginary part of complex heterodyned data

\item[\texttt{UINT4 iirM}] second decimation factor

\item[\texttt{REAL4IIRFilter *iirFilter2Re}] second IIR filter to be applied to real part of complex heterodyned data

\item[\texttt{REAL4IIRFilter *iirFilter2Im}] second IIR filter to be applied to imaginary part of complex heterodyned data

\item[\texttt{UINT4 stats}] set to 1 to calculate only Vh and variance; 2 for Vh, var, kurt, skew, covar; 0 else
\end{description}

\subsubsection*{Structure \texttt{FineHeterodyneInput}}
\idx[Type]{FineHeterodyneInput}
\noindent This structure stores the input for the fine heterodyne.

\begin{description}
\item[\texttt{COMPLEX8TimeSeries Vh}]    heterodyned, averaged and resampled data 
\item[\texttt{COMPLEX8TimeSeries varh}]   variance of corresponding Vh
\item[\texttt{REAL4 f0}]  frequency of the signal 
\item[\texttt{REAL4 f1}]  first time derivative of frequency 
\item[\texttt{REAL4 f2}] second time derivative of frequency 
\item[\texttt{REAL8 fEpochGPS}] epoch of the frequency
\item[\texttt{SkyPosition source}] location of pulsar in sky - equatorial coordinate system 
\item[\texttt{REAL4 pmRA}] proper motion RA (radians / year)
\item[\texttt{REAL4 pmDEC}] proper motion DEC (radians / year)
\item[\texttt{REAL8 posEpochGPS}] epoch of RA and DEC
\item[\texttt{UINT4 model}] 0 for isolated pulsar, 1 for binary
\item[\texttt{REAL8 e}] eccentricity of orbit
\item[\texttt{REAL8 w}]
\item[\texttt{REAL8 T0}]
\item[\texttt{REAL8 Pb}]
\item[\texttt{REAL8 x}]
\item[\texttt{REAL8 lg}]
\end{description}

\subsubsection*{Structure \texttt{FineHeterodyneOutput}}
\idx[Type]{FineHeterodyneOutput}
\noindent This structure stores the output of the fine heterodyne.

\begin{description}
\item[\texttt{COMPLEX8TimeSeries B}] bin value 
\item[\texttt{COMPLEX8TimeSeries var}]  variance 
\item[\texttt{REAL4 phase}] phase
\end{description}



\subsubsection*{Structure \texttt{FineHeterodyneParams}}
\idx[Type]{FineHeterodyneParams}

\noindent This structure stores the params of the fine heterodyne.

\begin{description}
\item[\texttt{EphemerisData *edat}]
\item[\texttt{LALDetector detector}]
\item[\texttt{REAL4IIRFilter *iirFilterRe}] IIR filter to be applied to real part of complex heterodyned data 
\item[\texttt{REAL4IIRFilter *iirFilterIm}] IIR filter to be applied to imaginary part of complex heterodyned data 
\item[\texttt{UINT4 M}]  decimation factor
\item[\texttt{UINT2 iirFlag}]  1 to apply iir filter, 0 for no iir filter
\end{description}

\vfill{\footnotesize\input{HeterodynePulsarHV}}
\newpage\input{HeterodynePulsarC}
\newpage\input{HeterodynePulsarTestC}

********************************** </lalLaTeX> */

#ifndef _HETERODYNEPULSAR_H
#define _HETERODYNEPULSAR_H

#include <lal/LALStdlib.h>
/******* INCLUDE ANY OTHER LAL HEADERS needed for header (NOT module) ****/

#include <lal/IIRFilter.h>
#include <lal/ZPGFilter.h>
#include <lal/LALBarycenter.h>
#include <lal/SkyCoordinates.h> 
#include <lal/AVFactories.h>
#include <lal/BinaryPulsarTiming.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (HETERODYNEPULSARH, "$Id$");

/******************************** <lalErrTable file="HeterodynePulsarHE"> */
#define HETERODYNEPULSARH_ENULLINPUT 1
#define HETERODYNEPULSARH_ENULLOUTPUT 2
#define HETERODYNEPULSARH_ENULLPARAMS 3
#define HETERODYNEPULSARH_ERFACTOR 4
#define HETERODYNEPULSARH_EINVALIDF0 5
#define HETERODYNEPULSARH_ELENGTH 6
#define HETERODYNEPULSARH_EBINARY 7

#define HETERODYNEPULSARH_MSGENULLINPUT "Input was Null"
#define HETERODYNEPULSARH_MSGENULLOUTPUT "Output was Null"
#define HETERODYNEPULSARH_MSGENULLPARAMS "Params was Null"
#define HETERODYNEPULSARH_MSGERFACTOR "The decimation factor supplied was invalid"
#define HETERODYNEPULSARH_MSGEINVALIDF0 "Invalid input f0"
#define HETERODYNEPULSARH_MSGELENGTH "Input vectors were not the same length"
#define HETERODYNEPULSARH_MSGEBINARY "Binary model not yet implemented"

/************************************ </lalErrTable> */

/****** DEFINE NEW STRUCTURES AND TYPES ************/
typedef struct
tagCoarseHeterodyneInput
{
  REAL4TimeSeries 	V;	        /* calibrated strain data from detector */
  REAL4 		f0;	        /* heterodyning base frequency*/
} CoarseHeterodyneInput;

typedef struct
tagCoarseHeterodyneOutput
{
  COMPLEX8TimeSeries 	Vh;             /* heterodyned, averaged and resampled data */
  REAL4       		phase;          /* phase of the reference signal f0 at t0 */  
  COMPLEX16         	varh;           /* variance */
  COMPLEX16 	        avg;  		/* average */
  COMPLEX16	        kurt;		/* kurtosis */
  COMPLEX16             skew;		/* skewness */
  COMPLEX16             covar;          /* covariance */
  
} CoarseHeterodyneOutput;

typedef struct
tagCoarseHeterodyneParams
{  
  UINT4  		boxM;           /* 1st decimation factor */
  REAL4IIRFilter        *iirFilter1Re;   /* IIR filter to be applied to real part of complex heterodyned data */
  REAL4IIRFilter        *iirFilter1Im;   /* IIR filter to be applied to imaginary part of complex heterodyned data */
  UINT4			iirM; 		/* 2nd decimation factor */
  REAL4IIRFilter 	*iirFilter2Re;   /* IIR filter to be applied to real part of complex heterodyned data */
  REAL4IIRFilter 	*iirFilter2Im;   /* IIR filter to be applied to imaginary part of complex heterodyned data */
  UINT4			stats;          /* set 1 for var; to 2 to calculate var, average, kurtosis, and skewness; 0
  else */
} CoarseHeterodyneParams;


typedef struct
tagFineHeterodyneInput
{
  COMPLEX8TimeSeries    Vh;             /* heterodyned, averaged and resampled data */
  COMPLEX8TimeSeries    varh;           /* variance of the rFactor points that were averaged */
  REAL4			f0;             /* frequency of the signal */
  REAL4			f1;             /* first time derivative of frequency */
  REAL4			f2;             /* second time derivative of frequency */
  REAL8 		fEpochGPS;      /* epoch of the frequency at SSB */
  SkyPosition           source;         /* location of pulsar in sky - equatorial coordinate system */
  REAL4 		pmRA;           /* proper motion RA radians / year */
  REAL4			pmDEC;		/* proper motion DEC radians / year*/
  REAL8 		posEpochGPS;    /* epoch of RA and DEC */
  INT4			model;		/* 0 for isolated, 1 for keplerian , 2 for ... */
  REAL8			e;
  REAL8 		w;
  REAL8 		T0;
  REAL4			Pb;
  REAL8			x;
  REAL8 		lg;
} FineHeterodyneInput;

typedef struct
tagFineHeterodyneOutput
{
  COMPLEX16TimeSeries 	B;		/* bin value */
  COMPLEX16TimeSeries 	var;  		/* variance */
  REAL4			phase;		/* phase */
} FineHeterodyneOutput;

typedef struct
tagFineHeterodyneParams
{  
  EphemerisData *edat;
  LALDetector detector;
  REAL4IIRFilter 	*iirFilterRe;    /* IIR filter to be applied to real part of complex heterodyned data */
  REAL4IIRFilter 	*iirFilterIm;    /* IIR filter to be applied to imaginary part of complex heterodyned data */
  UINT4  		M;         	/* decimation factor */
  UINT2 		iirFlag;	/* if 1 iir, if 0 no iir */
} FineHeterodyneParams;

/****** INCLUDE EXTERNAL GLOBAL VARIABLES ************/

void
LALCoarseHeterodyne( LALStatus                   *status,
		     CoarseHeterodyneOutput      *output,
		     CoarseHeterodyneInput       *input,
		     CoarseHeterodyneParams      *params );

void
LALFineHeterodyneToPulsar ( LALStatus                   *status,
		  	    FineHeterodyneOutput        *output,
		   	    FineHeterodyneInput         *input,
		   	    FineHeterodyneParams        *params );		     

#ifdef  __cplusplus
}
#endif

#endif /* _HETERODYNEPULSAR_H */
