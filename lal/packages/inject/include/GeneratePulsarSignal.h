/************************************ <lalVerbatim file="GeneratePulsarSignalHV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{GeneratePulsarSignal.h}}
\label{s:GeneratePulsarSignal.h}

Provides structures and prototypes for the pulsar signal-generation
routines in GeneratePulsarSignal.c.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/GeneratePulsarSignal.h>
\end{verbatim}

\noindent 

******************************************************* </lalLaTeX> */

#ifndef _GENERATEPULSARSIGNAL_H  /* Double-include protection. */
#define _GENERATEPULSARSIGNAL_H

#include <lal/LALDatatypes.h>
#include <lal/DetectorSite.h>
#include <lal/GenerateSpinOrbitCW.h>
#include <lal/Date.h>
#include <lal/LALBarycenter.h>
#include <lal/PulsarDataTypes.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( GENERATEPULSARSIGNALH, "$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
***************************************************** <lalErrTable> */
#define GENERATEPULSARSIGNALH_ENULL 		1
#define GENERATEPULSARSIGNALH_ENONULL		2
#define GENERATEPULSARSIGNALH_EMEM		3
#define GENERATEPULSARSIGNALH_ESAMPLING		4
#define GENERATEPULSARSIGNALH_ESSBCONVERT	5
#define GENERATEPULSARSIGNALH_ESYS		6
#define GENERATEPULSARSIGNALH_ETIMEBOUND	7
#define GENERATEPULSARSIGNALH_ENUMSFTS		8
#define GENERATEPULSARSIGNALH_EINCONSBAND	9
#define GENERATEPULSARSIGNALH_ENOISEDELTAF	10
#define GENERATEPULSARSIGNALH_ENOISEBAND	11
#define GENERATEPULSARSIGNALH_ENOISEBINS	12


#define GENERATEPULSARSIGNALH_MSGENULL 		"Arguments contained an unexpected null pointer"
#define GENERATEPULSARSIGNALH_MSGENONULL	"Output pointer is not NULL"
#define GENERATEPULSARSIGNALH_MSGEMEM		"Out of memory"
#define GENERATEPULSARSIGNALH_MSGESAMPLING	"Waveform sampling interval too large."
#define GENERATEPULSARSIGNALH_MSGESSBCONVERT	"SSB->GPS iterative conversion failed"
#define GENERATEPULSARSIGNALH_MSGESYS		"System error, probably while File I/O"
#define GENERATEPULSARSIGNALH_MSGETIMEBOUND	"Timestamp outside of allowed time-interval"
#define GENERATEPULSARSIGNALH_MSGENUMSFTS	"Inconsistent number of SFTs in timestamps and noise-SFTs"
#define GENERATEPULSARSIGNALH_MSGEINCONSBAND	"Inconsistent values of sampling-rate and Tsft"
#define GENERATEPULSARSIGNALH_MSGENOISEDELTAF	"Frequency resolution of noise-SFTs inconsistent with signal"
#define GENERATEPULSARSIGNALH_MSGENOISEBAND	"Frequency band of noise-SFTs inconsistent with signal"
#define GENERATEPULSARSIGNALH_MSGENOISEBINS	"Frequency bins of noise-SFTs inconsistent with signal"

/*************************************************** </lalErrTable> */

/*************************************************** 
<lalLaTeX>

\subsection*{Types}
\idx[Type]{PulsarSignalParams}
\idx[Type]{SFTParams}

\subsubsection*{Structure \texttt{PulsarSignalParams}}

Parameter-structure for \verb+LALGeneratePulsarSignal()+. Defines the
source-parameters (of the pulsar and its binary orbit, if any), the
detector-location and the time-series to be produced. 

</lalLaTeX> */
/* <lalVerbatim> */
typedef struct {
  /* source-parameters */
  PulsarSourceParams pulsar;	/* defining the actual pulsar-source */
  BinaryOrbitParams *orbit;	/* and its binary orbit (NULL if isolated pulsar) */
  
  /* characterize the detector */
  COMPLEX8FrequencySeries *transfer; /* detector transfer function (NULL if not used) */	
  LALDetector *site;		/* detector location and orientation */  
  EphemerisData *ephemerides;	/* Earth and Sun ephemerides */
  
  /* characterize the output time-series */
  LIGOTimeGPS startTimeGPS;     /* start time of output time series */
  UINT4 duration;           	/* length of time series in seconds */
  REAL8 samplingRate;		/* sampling rate of time-series (= 2 * frequency-Band) */
  REAL8 fHeterodyne;		/* heterodyning frequency for output time-series */
} PulsarSignalParams;
/* </lalVerbatim> */
/*<lalLaTeX>

\subsubsection*{Structure \texttt{SFTParams}}

Parameters defining the SFTs to be returned from \verb+LALSignalToSFTs()+.

</lalLaTeX> */
/* <lalVerbatim> */
typedef struct {
  REAL8 Tsft;			 /* length of each SFT in seconds */
  LIGOTimeGPSVector *timestamps; /* timestamps to produce SFTs for (can be NULL) */
  SFTVector *noiseSFTs;		 /* noise SFTs to be added (can be NULL) */
} SFTParams;
/* </lalVerbatim> */


/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{GeneratePulsarSignalHV}}
\newpage\input{GeneratePulsarSignalC}
******************************************************* </lalLaTeX> */


/* Function prototypes */
void LALGeneratePulsarSignal (LALStatus *stat, REAL4TimeSeries **signal, const PulsarSignalParams *params);
void LALSignalToSFTs (LALStatus *stat, SFTVector **outputSFTs, const REAL4TimeSeries *signal, const SFTParams *params);

void LALConvertGPS2SSB (LALStatus* stat, LIGOTimeGPS *SSBout, LIGOTimeGPS GPSin, const PulsarSignalParams *params);
void LALConvertSSB2GPS (LALStatus *stat, LIGOTimeGPS *GPSout, LIGOTimeGPS GPSin, const PulsarSignalParams *params);

/* interanl debug- und testing functions */
void dump_SFT (LALStatus *stat, const SFTtype *sft, const CHAR *fname);
void write_SFT (LALStatus *stat, const SFTtype *sft, const CHAR *fname);
void LALwriteSFTtoXMGR (LALStatus *stat, const SFTtype *sft, const CHAR *fname);
void PrintR4TimeSeries (LALStatus *stat, const REAL4TimeSeries *series, const CHAR *fname);
void PrintGWSignal (LALStatus *stat, const CoherentGW *signal, const CHAR *fname);
void compare_SFTs (const SFTtype *sft1, const SFTtype *sft2);


/********************************************************** <lalLaTeX>
%% \newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
