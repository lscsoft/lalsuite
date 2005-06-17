/************************************ <lalVerbatim file="GeneratePulsarSignalHV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/* NOTES: */
/* 07/14/04 gam; add functions LALFastGeneratePulsarSFTs and LALComputeSkyAndZeroPsiAMResponse */
/* 10/08/04 gam; fix indexing into trig lookup tables (LUTs) by having table go from -2*pi to 2*pi */

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
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/GenerateSpinOrbitCW.h>
#include <lal/Date.h>
#include <lal/LALBarycenter.h>
#include <lal/PulsarDataTypes.h>
#include <lal/ComputeSky.h>
#include <lal/ComputeSkyBinary.h>

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
#define GENERATEPULSARSIGNALH_EBADCOORDS	13
#define GENERATEPULSARSIGNALH_ELUTS		14

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
#define GENERATEPULSARSIGNALH_MSGEBADCOORDS	"Current code requires sky position in equatorial coordinates"
#define GENERATEPULSARSIGNALH_MSGELUTS		"Lookup tables (LUTs) for trig functions must be defined on domain -2pi to 2pi inclusive"
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
/** input parameters to GeneratePulsarSignal(), defining the source and the time-series */
typedef struct {
  /* source-parameters */
  PulsarSourceParams pulsar;	/**< the actual pulsar-source */
  BinaryOrbitParams *orbit;	/**< and its binary orbit (NULL if isolated pulsar) */
  
  /* characterize the detector */
  COMPLEX8FrequencySeries *transfer; /**< detector transfer function (NULL if not used) */	
  LALDetector *site;		/**< detector location and orientation */  
  EphemerisData *ephemerides;	/**< Earth and Sun ephemerides */
  
  /* characterize the output time-series */
  LIGOTimeGPS startTimeGPS;     /**< start time of output time series */
  UINT4 duration;           	/**< length of time series in seconds */
  REAL8 samplingRate;		/**< sampling rate of time-series (= 2 * frequency-Band) */
  REAL8 fHeterodyne;		/**< heterodyning frequency for output time-series */
} PulsarSignalParams;
/* </lalVerbatim> */
/*<lalLaTeX>

\subsubsection*{Structure \texttt{SFTParams}}

Parameters defining the SFTs to be returned from \verb+LALSignalToSFTs()+.

</lalLaTeX> */
/* <lalVerbatim> */
/** input parameters to LALSignalToSFTs() */
typedef struct {
  REAL8 Tsft;			 /**< length of each SFT in seconds */
  LIGOTimeGPSVector *timestamps; /**< timestamps to produce SFTs for (can be NULL) */
  SFTVector *noiseSFTs;		 /**< noise SFTs to be added (can be NULL) */
} SFTParams;
/* </lalVerbatim> */

/*<lalLaTeX>
\subsubsection*{Structure \texttt{SFTandSignalParams}}

Parameters defining the pulsar signal and SFTs used by \verb+LALFastGeneratePulsarSFTs()+.  Lookup tables (LUTs) are
used for trig functions if \verb+resTrig+ $> 0$; the user must then initialize \verb+trigArg+, \verb+sinVal+, and
\verb+cosVal+ on the domain -2pi to 2pi inclusive.  See \verb+GeneratePulsarSignalTest.c+ for an example.

</lalLaTeX> */
/* <lalVerbatim> */
typedef struct {
   PulsarSignalParams *pSigParams; 
   SFTParams *pSFTParams;
   INT4  nSamples;  /* nsample from noise SFT header; 2x this equals effective number of time samples  */
   INT4  resTrig;   /* length sinVal, cosVal; domain: -2pi to 2pi; resolution = 4pi/resTrig */
   REAL8 *trigArg;  /* array of arguments to hold lookup table (LUT) values for doing trig calls */
   REAL8 *sinVal;   /* sinVal holds lookup table (LUT) values for doing trig sin calls */
   REAL8 *cosVal;   /* cosVal holds lookup table (LUT) values for doing trig cos calls */
} SFTandSignalParams;
/* </lalVerbatim> */

/*<lalLaTeX>
\subsubsection*{Structure \texttt{SkyConstAndZeroPsiAMResponse}}

Sky Constants and beam pattern response functions used by \verb+LALFastGeneratePulsarSFTs()+.
These are output from \verb+LALComputeSkyAndZeroPsiAMResponse()+.

</lalLaTeX> */
/* <lalVerbatim> */
typedef struct {
      REAL8  *skyConst;      /* vector of A and B sky constants */
      REAL4  *fPlusZeroPsi;  /* vector of Fplus values for psi = 0 at midpoint of each SFT */
      REAL4  *fCrossZeroPsi; /* vector of Fcross values for psi = 0 at midpoint of each SFT */
} SkyConstAndZeroPsiAMResponse;
/* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{GeneratePulsarSignalHV}}

\newpage\input{GeneratePulsarSignalC}

\newpage\input{SimulatePulsarSignalC}
******************************************************* </lalLaTeX> */


/* Function prototypes */
void LALGeneratePulsarSignal (LALStatus *, REAL4TimeSeries **signal, const PulsarSignalParams *params);
void LALSimulatePulsarSignal (LALStatus *, REAL8TimeSeries **timeSeries, const PulsarSignalParams *params);

void LALSignalToSFTs (LALStatus *, SFTVector **outputSFTs, const REAL4TimeSeries *signal, const SFTParams *params);

void LALComputeSkyAndZeroPsiAMResponse (LALStatus *, SkyConstAndZeroPsiAMResponse *output, const SFTandSignalParams *params);
void LALFastGeneratePulsarSFTs (LALStatus *, SFTVector **outputSFTs, const SkyConstAndZeroPsiAMResponse *input, const SFTandSignalParams *params);

void LALConvertGPS2SSB (LALStatus* , LIGOTimeGPS *SSBout, LIGOTimeGPS GPSin, const PulsarSignalParams *params);
void LALConvertSSB2GPS (LALStatus *, LIGOTimeGPS *GPSout, LIGOTimeGPS GPSin, const PulsarSignalParams *params);

void LALMakeTimestamps (LALStatus *, LIGOTimeGPSVector **timestamps, const LIGOTimeGPS tStart, REAL8 duration, REAL8 Tsft);

/********************************************************** <lalLaTeX>
%% \newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
