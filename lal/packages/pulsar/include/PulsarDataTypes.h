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
typedef struct {
  UINT4 			length;	
  COMPLEX8FrequencySeries 	*data;	
} COMPLEX8FrequencySeriesVector;

/* and a bunch of timestamps */
typedef struct {
  UINT4 	length;
  LIGOTimeGPS 	*data;
} LIGOTimeGPSVector;
/* </lalVerbatim> */
/*<lalLaTeX>

\subsubsection*{Structure \texttt{PulsarSourceParams}}

Defines the astrophysical parameters of the pulsar. 

</lalLaTeX> */
/* <lalVerbatim> */
typedef struct {
  LIGOTimeGPS tRef;	/* reference GPS time of pulsar parameters */
  SkyPosition position;	/* source location (in radians) */
  REAL4 psi;            /* polarization angle (radians) at TRef */
  REAL4 aPlus, aCross;  /* polarization amplitudes at TRef */
  REAL8 phi0;           /* initial phase (radians) at TRef */
  REAL8 f0;             /* initial frequency (Hz) at TRef */
  REAL8Vector *spindown;/* frequency spindowns at TRef (NOT f0-normalized!) */
} PulsarSourceParams;
/* </lalVerbatim> */
/*<lalLaTeX>

\subsubsection*{Structure \texttt{BinaryOrbitParams}}

Defines the astrophysical parameters of the binary orbit of the pulsar.

</lalLaTeX> */
/* <lalVerbatim> */
typedef struct {
  LIGOTimeGPS orbitEpoch; /* time of periapsis passage */
  REAL8 omega;            /* argument of periapsis (radians) */
  REAL8 rPeriNorm;        /* projected, normalized periapsis (s) */
  REAL8 oneMinusEcc;      /* 1 - orbital eccentricity */
  REAL8 angularSpeed;     /* angular speed at periapsis (Hz) */
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
void LALCreateSFTtype (LALStatus *stat, SFTtype **sft, UINT4 SFTlen);	
void LALCreateSFTVector (LALStatus *stat, SFTVector **sftvect, UINT4 numSFTs, UINT4 SFTlen);
void LALDestroySFTtype (LALStatus *stat, SFTtype **sft);
void LALDestroySFTVector (LALStatus *stat, SFTVector **sftvect);
void LALCopySFTtype (LALStatus *stat, SFTtype *dest, const SFTtype *src);

void LALCreateTimestampVector (LALStatus *stat, LIGOTimeGPSVector **vect, UINT4 len);
void LALDestroyTimestampVector (LALStatus *stat, LIGOTimeGPSVector **vect);

  
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
