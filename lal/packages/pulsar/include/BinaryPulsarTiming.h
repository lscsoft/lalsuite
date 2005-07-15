/********************************* <lalVerbatim file="BinaryPulsarTimingHV">
Author: Dupuis, R. J.
$Id$
********************************** </lalVerbatim> */

/********************************* <lalLaTeX>

\section{Header \texttt{BinaryPulsarTiming.h}}

Provides routines to calculated time delay in binaries...

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/BinaryPulsarTiming.h>
\end{verbatim}
The gravitational wave signal from

The function \texttt{LALBinaryPulsarDeltaT()} ...

\begin{verbatim}
void
LALBinaryPulsarDeltaT(	LALStatus              	*status,
		  	REAL8			deltaT,
		   	BinaryPulsarParameters 	*input,
		   	INT2	   		model )
				
\end{verbatim}

\subsection*{Error conditions}
\input{BinaryPulsarTimingHE}

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
\item[\texttt{REAL4TimeSeries Vh}] Heterodyned data.

\item[\texttt{REAL4 varh}] Variance corresponding to each resampled data point.

\item[\texttt{REAL4 phase}] phase of the reference signal f0 at t0
\end{description}

\subsubsection*{Structure \texttt{CoarseHeterodyneParams}}
\idx[Type]{CoarseHeterodyneParams}

\noindent This structure stores parameters for the coarse heterodyne.

\begin{description}
\item[\texttt{REAL4IIRFilter *iirFilterRe}] IIR filter to be applied to real part of complex heterodyned data

\item[\texttt{REAL4IIRFilter *iirFilterIm}] IIR filter to be applied to imaginary part of complex heterodyned data

\item[\texttt{UINT4 rFactor}] reduction factor - number of points to average before resampling
\end{description}

\subsubsection*{Structure \texttt{FineHeterodyneInput}}
\idx[Type]{FineHeterodyneInput}
\noindent This structure stores the input for the fine heterodyne.

\begin{description}
\item[\texttt{COMPLEX8TimeSeries Vh}]    heterodyned, averaged and resampled data 
\item[\texttt{COMPLEX8TimeSeries varh}]   variance of the rFactor points that were averaged 
\item[\texttt{REAL4 f0}]  frequency of the signal 
\item[\texttt{REAL4 f1}]  first time derivative of frequency 
\item[\texttt{REAL4 f2}] second time derivative of frequency 
\item[\texttt{LIGOTimeGPS tEpochF}] epoch of the frequency
\item[\texttt{SkyPosition source}] location of pulsar in sky - equatorial coordinate system 
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
\item[\texttt{UINT4 rFactor}]  reduction factor - number of points to average before resampling 
\end{description}

\vfill{\footnotesize\input{BinaryPulsarTimingHV}}
\newpage\input{BinaryPulsarTimingC}
\newpage\input{BinaryPulsarTimingTestC}

********************************** </lalLaTeX> */

#ifndef _BinaryPulsarTiming_H
#define _BinaryPulsarTiming_H

#include <lal/LALStdlib.h>
/******* INCLUDE ANY OTHER LAL HEADERS needed for header (NOT module) ****/

#include <lal/IIRFilter.h>
#include <lal/ZPGFilter.h>
#include <lal/LALBarycenter.h>
#include <lal/SkyCoordinates.h> 
#include <lal/AVFactories.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (BinaryPulsarTimingH, "$Id$");

/******************************** <lalErrTable file="BinaryPulsarTimingHE"> */
#define BINARYPULSARTIMINGH_ENULLINPUT 1
#define BINARYPULSARTIMINGH_ENULLOUTPUT 2
#define BINARYPULSARTIMINGH_ENULLPARAMS 3
#define BINARYPULSARTIMINGH_ERFACTOR 4
#define BINARYPULSARTIMINGH_EINVALIDF0 5
#define BINARYPULSARTIMINGH_ELENGTH 6

#define BINARYPULSARTIMINGH_MSGENULLINPUT "Input was Null"
#define BINARYPULSARTIMINGH_MSGENULLOUTPUT "Output was Null"
#define BINARYPULSARTIMINGH_MSGENULLPARAMS "Params was Null"
#define BINARYPULSARTIMINGH_MSGERFACTOR "The decimation factor supplied was invalid"
#define BINARYPULSARTIMINGH_MSGEINVALIDF0 "Invalid input f0"
#define BINARYPULSARTIMINGH_MSGELENGTH "Input vectors were not the same length"

/************************************ </lalErrTable> */

/****** DEFINE NEW STRUCTURES AND TYPES ************/
typedef struct
tagBinaryPulsarParameters
{
  REAL8 	Pb;
  REAL8		x;
  REAL8		e;
  REAL8		w;
  REAL8		T0;
  REAL8 	lg; /*longitude at T0 (usually 0)*/
  INT4	 	model;
} BinaryPulsarParameters;

typedef struct
tagBinaryPulsarTiming
{
  REAL8 	deltaT;

} BinaryPulsarTiming;
/****** INCLUDE EXTERNAL GLOBAL VARIABLES ************/
void
LALBinaryPulsarDeltaT(	LALStatus              	*status,
		  	BinaryPulsarTiming	*output,
			REAL8 			tgps,
		   	BinaryPulsarParameters 	*params);
    

#ifdef  __cplusplus
}
#endif

#endif /* _BinaryPulsarTiming_H */
