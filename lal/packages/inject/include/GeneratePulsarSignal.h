/************************************ <lalVerbatim file="GeneratePulsarSignalHV">
Author: Prix, Reinhard
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{GeneratePulsarSignal.h}}
\label{s:GeneratePulsarSignal.h}

Header file for GeneratePulsarSignal.

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

  
/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{GeneratePulsarSignalHV}}
\newpage\input{GeneratePulsarSignalC}
******************************************************* </lalLaTeX> */

/* New structures and types */

typedef struct {
  LIGOTimeGPS TRefSSB;	/* reference time for pulsar parameters (in SSB time!) */
  			/* if ZERO, startTimeGPS is used instead */
  SkyPosition position;	/* source location (in radians!) */
  REAL4 psi;            /* polarization angle (radians) at TRef */
  REAL4 aPlus, aCross;    /* polarization amplitudes at TRef */
  REAL8 phi0;             /* initial phase (radians) at TRef */
  REAL8 f0;               /* initial frequency (Hz) at TRef */
  REAL8Vector *spindown;  /* frequency spindowns at TRef (NOT f0-normalized!) */
} PulsarSourceParams;

typedef struct {
  LIGOTimeGPS orbitEpoch; /* time of a periapsis passage */
  REAL8 omega;            /* argument of periapsis (radians) */
  REAL8 rPeriNorm;        /* projected, normalized periapsis (s) */
  REAL8 oneMinusEcc;      /* 1 - orbital eccentricity */
  REAL8 angularSpeed;     /* angular speed at periapsis (Hz) */
} BinaryOrbitParams;

typedef struct {
  /* source-params */
  PulsarSourceParams pulsar;	  /* defining the actual pulsar-source */
  BinaryOrbitParams *orbit;	  /* and its binary orbit if applicable (NULL if not) */
  
  /* characterize the detector */
  COMPLEX8FrequencySeries *transferFunction;    /* frequency transfer function (NULL if not used) */	
  LALDetector *site;				/* detector location and orientation */  
  EphemerisData *ephemerides;			/* Earth and Sun ephemerides */
  
  /* characterize the output time-series */
  LIGOTimeGPS startTimeGPS;     /* start time of output time series */
  UINT4 duration;           	/* length of time series in seconds */
  REAL8 samplingRate;		/* sampling rate of time-series (= 2 * frequency-Band) */
  REAL8 fHeterodyne;		/* heterodyning frequency for output time-series */
} PulsarSignalParams;

/* we need a type for a vector of timestamps */
/* FIXME: move type into LAL */
typedef struct {
  UINT4 	length;
  LIGOTimeGPS 	*data;
} LIGOTimeGPSVector;


/* just for reference: this is the FrequencySeries type: 
 * { 
 *   CHAR              name[LALNameLength]; 
 *   LIGOTimeGPS       epoch; 
 *   REAL8             f0;	 
 *   REAL8             deltaF; 
 *   LALUnit           sampleUnits
 *   COMPLEX8Sequence *data; 
 * } 
 * COMPLEX8FrequencySeries; 
 */

/* we need a "Vector"-type of frequency-series  */
/* FIXME: move type into LAL */
typedef struct {
  UINT4 			length;		/* number of frequency-series */
  COMPLEX8FrequencySeries 	*data;		/* array of frequency-series */
} COMPLEX8FrequencySeriesVector;

/* for the lazy */
typedef COMPLEX8FrequencySeries 	SFTtype;	
typedef COMPLEX8FrequencySeriesVector 	SFTVector;

/* this is the current SFT-header in the file-format for storing SFTs
  struct headertag {
    REAL8 endian;
    INT4  gps_sec;
    INT4  gps_nsec;
    REAL8 tbase;
    INT4  firstfreqindex;
    INT4  nsamples;
  } header;
*/

/* parameter-struct for LALSignalToSFTs() */
typedef struct {
  REAL8 Tsft;				/* length of an SFT in seconds */
  LIGOTimeGPSVector *timestamps;	/* timestamps to use for SFT's (can be NULL) */
  SFTVector *noiseSFTs;			/* noise SFTs to be added to the output (can be NULL) */
} SFTParams;


/* Function prototypes */
void LALGeneratePulsarSignal (LALStatus *stat, REAL4TimeSeries **signal, const PulsarSignalParams *params);
void LALSignalToSFTs (LALStatus *stat, SFTVector **outputSFTs, const REAL4TimeSeries *signal, const SFTParams *params);

void LALNormalizeSkyPosition (LALStatus *stat, SkyPosition *posOut, const SkyPosition *posIn);
void dump_SFT (LALStatus *stat, const SFTtype *sft, const CHAR *fname);
void write_SFT (LALStatus *stat, const SFTtype *sft, const CHAR *fname);
void LALwriteSFTtoXMGR (LALStatus *stat, const SFTtype *sft, const CHAR *fname);
void PrintR4TimeSeries (LALStatus *stat, const REAL4TimeSeries *series, const CHAR *fname);
void PrintGWSignal (LALStatus *stat, const CoherentGW *signal, const CHAR *fname);
void ConvertGPS2SSB (LALStatus* stat, LIGOTimeGPS *SSBout, LIGOTimeGPS GPSin, const PulsarSignalParams *params);
void ConvertSSB2GPS (LALStatus *stat, LIGOTimeGPS *GPSout, LIGOTimeGPS GPSin, const PulsarSignalParams *params);
void compare_SFTs (const SFTtype *sft1, const SFTtype *sft2);
void LALCreateSFT (LALStatus *stat, SFTtype **outputSFT, UINT4 length);

void LALCreateSFTVector (LALStatus *stat, SFTVector **output, UINT4 numSFTs, UINT4 SFTlen);
void LALDestroySFTVector (LALStatus *stat, SFTVector **vect);
void LALDestroyTimestampVector (LALStatus *stat, LIGOTimeGPSVector **vect);

/********************************************************** <lalLaTeX>
%% \newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
