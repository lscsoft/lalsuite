/************************************ <lalVerbatim file="ComputeSkyHV">
Author:  Landry, M., and Mendell, G.
$Id$
************************************* </lalVerbatim> */

/* REVISIONS: */
/* 04/26/04 gam; Change LALStackSlide to StackSlide and LALSTACKSLIDE to STACKSLIDE for initial entry to LALapps. */
/* 06/05/04 gam; Add gpsStartTimeSec and gpsStartTimeNan to StackSlideSkyParams; set these to epoch that gives T0 at SSB. */
/* 12/03/04 gam; Clean up indentation; remove extraneous or obsolete comments. */
/* 12/03/04 gam; Add parameter: BOOLEAN divideSUMsByNumSTKs. */

/* <lalLaTeX> 
\section{Header \texttt{StackSlide.h}}
\label{s:StackSlide.h}
Computes frequency model, slide stacks accordingly, and sum them.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/StackSlide.h>
\end{verbatim}
\noindent  This is a short summary...
</lalLaTeX> */

#ifndef _STACKSLIDE_H
#define _STACKSLIDE_H
 
#include <lal/LALStdlib.h>
#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALBarycenter.h>
 
#ifdef __cplusplus
extern "C" {
#endif
  
NRCSID (STACKSLIDEH, "$Id$");
  
/* Author-defined error codes */
/* <lalLaTeX>  
     \subsection*{Error conditions}
     \vspace{0.1in}
     \input{StackSlideHErrorTable}
     
</lalLaTeX> */
  
/* <lalErrTable file="ComputeSkyHErrorTable"> */
#define STACKSLIDEH_ENULL 1
#define STACKSLIDEH_ENNUL 2
#define STACKSLIDEH_ENEGA 4
#define STACKSLIDEH_MSGENULL "Null Pointer"
#define STACKSLIDEH_MSGENNUL "Non-Null Pointer"
#define STACKSLIDEH_MSGENEGA "Bad Negative Value"
#define STACKSLIDECOMPUTESKYH_ENULL 6
#define STACKSLIDECOMPUTESKYH_ENNUL 8
#define STACKSLIDECOMPUTESKYH_ENEGA 10
#define STACKSLIDECOMPUTESKYH_MSGENULL "Null Pointer in StackSlideComputeSky"
#define STACKSLIDECOMPUTESKYH_MSGENNUL "Non-Null Pointer in StackSlideComputeSky"
#define STACKSLIDECOMPUTESKYH_MSGENEGA "Bad Negative Value in StackSlideComputeSky"
/* </lalErrTable>  */

/* <lalLaTeX>
\subsection*{Structures}

\begin{verbatim}
struct StackSlideParams
\end{verbatim}
\index{\texttt{StackSlideParams}}
\noindent This structure contains the parameters for the \verb@StackSlide()@ routine.  The parameters are:
\begin{description}
\end{description}
</lalLaTeX> */

typedef struct
tagStackSlideParams
{
	REAL8 **skyPosData;  
	REAL8 **freqDerivData;  
	INT4 numSkyPosTotal;
	INT4 numFreqDerivTotal;
	REAL8 f0STK;
	REAL8 f0SUM;
	REAL8 tSTK;
	REAL8 tSUM;
	INT4  nBinsPerSUM;
	INT4  numSTKs;
	REAL8 dfSUM;
	UINT4 gpsStartTimeSec;
	UINT4 gpsStartTimeNan;
	LIGOTimeGPS *timeStamps;
	INT4 numSpinDown;
	EphemerisData *edat;
	BarycenterInput *baryinput;
	BOOLEAN divideSUMsByNumSTKs;
}
StackSlideParams;

typedef struct
tagStackSlideSkyParams
{
	INT8            spinDwnOrder;   /* max spindown parameter order */
	INT8            mObsSFT;        /* number of SFTs */
	REAL8           tSFT;           /* timescale of SFT */
	LIGOTimeGPS     *tGPS;          /* GPS time of 1st data sample of each SFT */
	UINT4 gpsStartTimeSec;          /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
	UINT4 gpsStartTimeNan;          /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
	REAL8           *skyPos;        /* array of sky positions */
	BarycenterInput *baryinput;
	EmissionTime *emit;
	EarthState *earth;
	EphemerisData *edat;
}
StackSlideSkyParams;

typedef struct
tagTdotsAndDeltaTs
{
	REAL8  *vecTDots;    /* 1-d array of (dT/dt)'s for frequency calculation */
	REAL8  **vecDeltaTs; /* 2-d array of (T-T_0)'s for frequency calculation */
}
TdotsAndDeltaTs;

void StackSlide(	LALStatus *status, 
			REAL4FrequencySeries **SUMData, 
			REAL4FrequencySeries **STKData, 
			StackSlideParams *params);

void StackSlideComputeSky (LALStatus 	*status, 
			TdotsAndDeltaTs 	*pTdotsAndDeltaTs,
			StackSlideSkyParams 	*params);

#ifdef __cplusplus
}
#endif

#endif /* _STACKSLIDE_H */
