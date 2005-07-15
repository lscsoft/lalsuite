

/*
<lalVerbatim file="SimulatePopcornHV">
Author: Tania Regimbau 
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{SimulatePopcorn.h}}
\label{s:SimulatePopcorn.h}

Provides prototype for simulating whitened time-domain signals in a pair 
of detectors that arises from low duty cycle astrophysical backgrounds.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/SimulatePopcorn.h>
\end{verbatim}

</lalLaTeX>
*/


#ifndef _SIMULATEPOPCOR_H
#define _SIMULATEPOPCOR_H


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/StochasticCrossCorrelation.h> 
#include <lal/AVFactories.h>
#include <lal/RealFFT.h>
#include <lal/ComplexFFT.h>
#include <lal/Units.h>
#include <lal/Random.h>
#include <lal/DetectorSite.h>


#ifdef __cplusplus
extern "C" {
#endif

NRCSID (SIMULATEPOPCORNH, "$Id$");
/*
<lalLaTeX>
\subsection*{Error conditions}
\input{SimulatePopcornHErrTab}
</lalLaTeX>
*/
/*
<lalErrTable file="SimulatePopcornHErrTab"> 
*/

#define SIMULATEPOPCORNH_ENULLP          1
#define SIMULATEPOPCORNH_ENONNULLFMIN    2
#define SIMULATEPOPCORNH_EMMDELTA        3
#define SIMULATEPOPCORNH_EMMLEN          4
#define SIMULATEPOPCORNH_EBV             5


#define SIMULATEPOPCORNH_MSGENULLP         "Null pointer"
#define SIMULATEPOPCORNH_MSGENONNULLFMIN   "Non zero start frequency" 
#define SIMULATEPOPCORNH_MSGEMMDELTA       "Mismatch in sequence spacings"
#define SIMULATEPOPCORNH_MSGEMMLEN         "Mismatch in sequence lengths"
#define SIMULATEPOPCORNH_MSGEBV            "Bad input or parameter"

/*</lalErrTable> */
/*<lalLaTeX>
\subsection*{Structures}
These constants define the cosmological model
\begin{verbatim}*/
#define SIMULATEPOPCORN_ho 0.7
#define SIMULATEPOPCORN_OMEGAMATTER 0.3
#define SIMULATEPOPCORN_OMEGAVACUUM 0.7
/*\end{verbatim}*/
/*</lalLaTeX>*/
/*<lalLaTeX>
\subsection*{Structures}
These are function pointers to functions that model burst waveforms.

\begin{verbatim}*/
typedef void (REAL4LALWform) (REAL4 *output, REAL4 input);
/*\end{verbatim}*/
/*</lalLaTeX>*/

/*<lalLaTeX>
The following structure contains the input of the simulation.
\begin{verbatim}*/
typedef struct tagSimPopcornInputStruc {
  REAL4LALWform   *inputwform; /*waveform of a single burst*/
  REAL4   inputduration; /*mean duration of a single burst*/
  REAL4   inputlambda; /*mean tims interval between successive bursts*/
  UINT4   inputNdataset; /*number of detector sites 1 for H1/H2, 2 for H/L*/
  INT2   inputsite0; /*first detector code*/
  INT2   inputsite1; /*second detector code*/  
  COMPLEX8FrequencySeries   *wfilter0; /*response of the first detector*/
  COMPLEX8FrequencySeries   *wfilter1; /*response of the second detector*/
  } SimPopcornInputStruc;
  /*\end{verbatim}*/
/*</lalLaTeX>*/
/*<lalLaTeX>

The following structure contains the parameters of the simulation.
\begin{verbatim}*/
typedef struct tagSimPopcornParamsStruc {
  UINT4   paramsstarttime; /*starting time*/
  UINT4   paramslength; /*length of the time serie in s*/
  UINT4   paramssrate; /*sampling rate of the time serie in Hz*/
  UINT4   paramsseed; /*random generator seed*/
  REAL8   paramsfref; /*reference frequency if normalization, -1 otherwise*/
  } SimPopcornParamsStruc;
/*\end{verbatim}*/
/*</lalLaTeX>*/
/*<lalLaTeX>

The following structure contains the simulated pair time series and $\Omega$ spectrum
\begin{verbatim}*/
typedef struct tagSimPopcornOutputStruc {
  REAL4TimeSeries   *SimPopcorn0;  
  REAL4TimeSeries   *SimPopcorn1;
  REAL4FrequencySeries   *omega0;
  REAL4FrequencySeries   *omega1;
  } SimPopcornOutputStruc;
/*\end{verbatim}*/
/*</lalLaTeX>*/

void
LALSimPopcornTimeSeries (LALStatus *status,SimPopcornOutputStruc *output,
    SimPopcornInputStruc *input,SimPopcornParamsStruc *params);
#ifdef  __cplusplus
}
#endif
#endif
