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
#include "GeneratePulsarSignal.h"
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
/* #ifdef  __cplusplus    */
/* extern "C" { */
/* #endif */

NRCSID( PULSARSIGNALH, "$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
***************************************************** <lalErrTable> */
#define PULSARSIGNALH_ENULL 		1
#define PULSARSIGNALH_ENONULL		2
#define PULSARSIGNALH_EMEM		3
#define PULSARSIGNALH_ESAMPLING		4
#define PULSARSIGNALH_ESSBCONVERT	5
#define PULSARSIGNALH_ESYS		6
#define PULSARSIGNALH_ETIMEBOUND	7
#define PULSARSIGNALH_ENUMSFTS		8
#define PULSARSIGNALH_EINCONSBAND	9

#define PULSARSIGNALH_MSGENULL 		"Arguments contained an unexpected null pointer"
#define PULSARSIGNALH_MSGENONULL	"Output pointer is not NULL"
#define PULSARSIGNALH_MSGEMEM		"Out of memory"
#define PULSARSIGNALH_MSGESAMPLING	"Waveform sampling interval too large."
#define PULSARSIGNALH_MSGESSBCONVERT	"SSB->GPS iterative conversion failed"
#define PULSARSIGNALH_MSGESYS		"System error, probably while File I/O"
#define PULSARSIGNALH_MSGETIMEBOUND	"Timestamp outside of allowed time-interval"
#define PULSARSIGNALH_MSGENUMSFTS	"Inconsistent number of SFTs in timestamps and noise-SFTs"
#define PULSARSIGNALH_MSGEINCONSBAND	"Inconsistent values of sampling-rate and Tsft: number of samples is not integer!"

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
  REAL8Vector *f;         /* f0-normalized Taylor parameters at TRef */
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
  COMPLEX8FrequencySeries *transferFunction;    /* frequency transfer function */	
  						/* --> FIXME do we need this? */
  LALDetector *site;        		   	/* detector location and orientation */  
  EphemerisData *ephemerides;  			/* Earth and Sun ephemerides */
  
  /* characterize the output time-series */
  LIGOTimeGPS startTimeGPS;     /* start time of output time series */
  UINT4 duration;           	/* length of time series in s*/
  REAL8 samplingRate;		/* sampling rate of time-series (= 2 * frequency-Band) */
  REAL8 fHeterodyne;		/* heterodyning frequency for output time-series */
} PulsarSignalParams;

/* we need a type for a vector of timestamps */
typedef struct {
  UINT4 length;
  LIGOTimeGPS *data;
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
typedef COMPLEX8FrequencySeries SFTtype;	/* for the lazy */

/* we need a "Vector"-type of SFTs  */
typedef struct {
  UINT4 	numSFTs;	/* number of SFTs */
  SFTtype 	*SFTlist;	/* array of SFTs */
} SFTVector;


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


typedef struct {
  REAL8 Tsft;				/* length of an SFT in seconds */
  LIGOTimeGPSVector *timestamps;	/* timestamps to use for SFT's (can be NULL) */
  SFTVector *noiseSFTs;			/* noise SFTs to be added to the output (can be NULL) */
} SFTParams;



/* Function prototypes */
void LALGeneratePulsarSignal (LALStatus *stat, REAL4TimeSeries *signal, const PulsarSignalParams *params);
void LALSignalToSFTs (LALStatus *stat, SFTVector **outputSFTs, const REAL4TimeSeries *signal, const SFTParams *params);

void write_SFT (LALStatus *stat, const SFTtype *sft, const CHAR *fname);
void LALPrintR4TimeSeries (LALStatus *stat, const REAL4TimeSeries *series, const CHAR *fname);
void PrintGWSignal (LALStatus *stat, const CoherentGW *signal, const CHAR *fname);
void ConvertGPS2SSB (LALStatus* stat, LIGOTimeGPS *SSBout, LIGOTimeGPS GPSin, const PulsarSignalParams *params);
void ConvertSSB2GPS (LALStatus *stat, LIGOTimeGPS *GPSout, LIGOTimeGPS GPSin, const PulsarSignalParams *params);
void compare_SFTs (const SFTtype *sft1, const SFTtype *sft2);
void LALCreateSFT (LALStatus *stat, SFTtype **outputSFT, UINT4 length);

void LALCreateSFTVector (LALStatus *stat, SFTVector **output, UINT4 numSFTs, UINT4 SFTlen);
void LALDestroySFTVector (LALStatus *stat, SFTVector **vect);
void LALDestroyTimestampVector (LALStatus *stat, LIGOTimeGPSVector **vect);

/********************************************************** <lalLaTeX>
\newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

/* #ifdef  __cplusplus */
/* } */
/* #endif   */  
/* C++ protection. */

#endif  /* Double-include protection. */
