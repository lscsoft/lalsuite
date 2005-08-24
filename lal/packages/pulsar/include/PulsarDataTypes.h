/************************************ <lalVerbatim file="PulsarDataTypesHV">
Author: Prix, Reinhard, Sintes, A.M.
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\section{Header \texttt{PulsarDataTypes.h}}
\label{s:PulsarDataTypes.h}

Provides general structures and data-types used in pulsar-searches.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/PulsarDataTypes.h>
\end{verbatim}

\noindent 

******************************************************* </lalLaTeX> */

#ifndef _PULSARDATATYPES_H  /* Double-include protection. */
#define _PULSARDATATYPES_H

#include <lal/LALDatatypes.h>
#include <lal/DetectorSite.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>

/* C++ protection. */
#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( PULSARDATATYPESH, "$Id$");

/*************************************************** 
<lalLaTeX>

\subsection*{Types}
\idx[Type]{PulsarSourceParams}
\idx[Type]{BinaryOrbitParams}
\idx[Type]{SFTtype}
\idx[Type]{SFTVector}
\idx[Type]{COMPLEX8FrequencySeriesVector}
\idx[Type]{LIGOTimeGPSVector}
</lalLaTeX> */

/* put those back from LALDatatypes.h, where Jolien didn't want them..*/
/* we also need a bunch of FFTs */
/* <lalVerbatim> */
/** a vector of 'SFTs', which are of type COMPLEX8FrequencySeries */
typedef struct {
  UINT4 			length;		/**< number of SFTs */
  COMPLEX8FrequencySeries 	*data;		/**< array of SFTs */
} COMPLEX8FrequencySeriesVector;

/** a vector of periodograms */
typedef struct {
  UINT4                  length;
  REAL8FrequencySeries   *data;
} REAL8FrequencySeriesVector;

/** a vector of 'timestamps' of type LIGOTimeGPS */
typedef struct {
  UINT4 	length;	/**< number of timestamps */
  LIGOTimeGPS 	*data;	/**< array of timestamps */
} LIGOTimeGPSVector;
/* </lalVerbatim> */
/*<lalLaTeX>

\subsubsection*{Structure \texttt{PulsarSourceParams}}

Defines the astrophysical parameters of the pulsar. 

</lalLaTeX> */
/* <lalVerbatim> */
/** type defining the parameters of a pulsar-source of Gravitational waves */
typedef struct {
  LIGOTimeGPS tRef;	/**< reference time of pulsar parameters (in SSB!) */
  SkyPosition position;	/**< source location (in radians) */
  REAL4 psi;            /**< polarization angle (radians) at tRef */
  REAL4 aPlus; 		/**< plus-polarization amplitude at tRef */
  REAL4 aCross;  	/**< cross-polarization amplitude at tRef */
  REAL8 phi0;           /**< initial phase (radians) at tRef */
  REAL8 f0;             /**< WAVE-frequency(!) at tRef (in Hz) */
  REAL8Vector *spindown;/**< wave-frequency spindowns at tRef (NOT f0-normalized!) */
} PulsarSourceParams;
/* </lalVerbatim> */
/*<lalLaTeX>

\subsubsection*{Structure \texttt{BinaryOrbitParams}}

Defines the astrophysical parameters of the binary orbit of the pulsar.

</lalLaTeX> */
/* <lalVerbatim> */
/** type defining the orbital parameters of a binary pulsar */
typedef struct {
  LIGOTimeGPS orbitEpoch; /**< time of periapsis passage (in SSB) */
  REAL8 omega;            /**< argument of periapsis (radians) */
  REAL8 rPeriNorm;        /**< projected, normalized periapsis (s) */
  REAL8 oneMinusEcc;      /**< 1 - orbital eccentricity */
  REAL8 angularSpeed;     /**< angular speed at periapsis (Hz) */
} BinaryOrbitParams;
/* </lalVerbatim> */
/*<lalLaTeX>

\subsubsection*{Structures \texttt{SFTtype} and \texttt{SFTVector}}

These are trivial typedefs used here for simplicity.

</lalLaTeX> */
/* <lalVerbatim> */
typedef COMPLEX8FrequencySeries 	SFTtype;	
typedef COMPLEX8FrequencySeriesVector 	SFTVector;
/* </lalVerbatim> */


/* ----------------------------------------------------------------------
 *  some prototypes for general functions handling these data-types 
 *----------------------------------------------------------------------*/
void LALCreateSFTtype (LALStatus *status, SFTtype **sft, UINT4 SFTlen);	
void LALCreateSFTVector (LALStatus *status, SFTVector **sftvect, UINT4 numSFTs, UINT4 SFTlen);
void LALDestroySFTtype (LALStatus *status, SFTtype **sft);
void LALDestroySFTVector (LALStatus *status, SFTVector **sftvect);
void LALCopySFTtype (LALStatus *status, SFTtype *dest, const SFTtype *src);

void LALCreateTimestampVector (LALStatus *status, LIGOTimeGPSVector **vect, UINT4 len);
void LALDestroyTimestampVector (LALStatus *status, LIGOTimeGPSVector **vect);

  
/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{PulsarDataTypesHV}}
%%\newpage\input{PulsarDataTypesC}
******************************************************* </lalLaTeX> */


/********************************************************** <lalLaTeX>
%% \newpage\input{LALSampleTestC}
******************************************************* </lalLaTeX> */

#ifdef  __cplusplus
}
#endif  
/* C++ protection. */

#endif  /* Double-include protection. */
