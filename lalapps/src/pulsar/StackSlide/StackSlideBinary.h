/************************************ <lalVerbatim file="ComputeSkyHV">
Author:  Landry, M., and Mendell, G.
$Id$
************************************* </lalVerbatim> */

/* REVISIONS: */
/* 04/26/04 gam; Change LALStackSlide to StackSlide and LALSTACKSLIDE to STACKSLIDE for initial entry to LALapps. */
/* 06/05/04 gam; Add gpsStartTimeSec and gpsStartTimeNan to StackSlideSkyParams; set these to epoch that gives T0 at SSB. */

/* <lalLaTeX> 
\section{Header \texttt{ComputeSky.h}}
\label{s:ComputeSky.h}
Computes phase coefficients necessary for a correct demodulation.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/ComputeSky.h>
\end{verbatim}

\noindent  This is a short summary of the analytical calculations which form the basis for the code in this routine.

Recall that a demodulated Fourier Transform (DeFT) is given by 
\begin{equation}
\hat{x}_b({\vec{\lambda}})=
\sum_{\alpha =0}^{M-1}\sum_{k=0}^{N-1}\tilde{x}_{\alpha k}\left[\frac{1}{N}\sum_{j=0}^{N-1}e^{-2\pi i(\Phi_{\alpha jb}(\vec{\lambda})-\frac{jk}{N})}\right]
\label{e4}
\end{equation}
The index $b$ defines the DeFT frequency bin, the index $\alpha$ loops through
the SFTs that build the DeFT, $k$ runs on all the SFT frequency bins, and $j$
is a time index that runs on each SFT.  As shown in section
\ref{s:LALDemod.h}, the next step in the development of the demodulation
technique involves Taylor expanding the phase model about the temporal
midpoint of each short segment of data, while retaining only first order
terms.  The Taylor expansion of $\Phi (t)$ about the temporal midpoint
$t_{\alpha,1/2}$ is
\begin{equation}
\Phi_{\alpha}(t) = \Phi(t_{\alpha,1/2})+\left[t-t_{\alpha,1/2}\right]\frac{d\Phi}{dt}(t_{\alpha,1/2})\label{taylor2} \\
\end{equation}
For each value of $\alpha$, this expression consist of either constant or linear terms in time.  With the particular time discretization chosen in this code, $t=t_{0}+(N\alpha+j)\ T_{obs}/NM$, we have
\begin{equation}
\label{time}
\left[t-t_{\alpha,1/2}\right]=\frac{\ T_{obs}}{M}\left(\frac{j}{N}-\frac{1}{2}\right)=\mathcal{T}_{s}\left(\frac{j}{N}-\frac{1}{2}\right),
\end{equation}
where $\mathcal{T}_{s}$ is the short time baseline of the $M$ short FTs.  On
the other hand, the phase can also be expressed as a function of SSB time $T$
(i.e. the time at the solar system barycenter).  We will assume the source to
be at rest in this reference frame.  Now, if one adopts the notation $\Delta
T_{\alpha}\equiv\left[T(t_{\alpha,1/2})-
T(t_{0})\right]$ and $\dot{T}_{\alpha}\equiv
dT/dt(t_{\alpha,1/2})$
the phase terms in the above equation are (neglecting constants)
\begin{eqnarray}
\Phi(t_{\alpha,1/2})                     & = & f_{0}\Delta T_{\alpha}+\frac{1}{2}f_{1}\Delta T_{\alpha}^{2}
+\frac{1}{3}f_{2}\Delta T_{\alpha}^{3}+\frac{1}{4}f_{3}\Delta T_{\alpha}^{4}+\frac{1}{5}f_{4}\Delta T_{\alpha}^{5}
+\frac{1}{6}f_{5}\Delta T_{\alpha}^{6} \nonumber\label{phi} \\ 
                                         &   & \\
\frac{d\Phi}{dt}(t_{\alpha,1/2})         & = & \dot{T}_{\alpha}\left(f_{0}+ f_{1}\Delta T_{\alpha}
+f_{2}\Delta T_{\alpha}^{2}+f_{3}\Delta T_{\alpha}^{3}
+f_{4}\Delta T_{\alpha}^{4}+f_{5}\Delta T_{\alpha}^{5}\right). \label{dphi}
\end{eqnarray}
These constants, for each value of $\alpha$, require $\dot{T}_{\alpha}$ and
$\Delta T_{\alpha}$, which are calculated by a suitable timing routine.  For
this demodulation package, this timing routine is provided by \verb@tdb()@.
Thus, for a given sky position, the timing routine will be called once for
each short time chunk, each call returning a specific  $\dot{T}_{\alpha}$ and
$\Delta T_{\alpha}$.  By substituting Eq.s \ref{time}, \ref{phi} and
\ref{dphi} in Eq. \ref{taylor2} and grouping together the terms in $j$ (linear
in $t$) in order to save computations, we have
\begin{equation}
\Phi_{\alpha}(t)=\sum_{s=0}^{n_{spin}}f_{s}A_{s\alpha}+\frac{j}{N}\sum_{s=0}^{n_{spin}}f_{s}B_{s\alpha},
\label{phasecalc}
\end{equation}
where $n_{spin}$ is the maximum order of spindown parameter.  Rather than
store the values of $\dot{T}_{\alpha}$ and $\Delta T_{\alpha}$ for each value
of $\alpha$, it is more efficient to calculate the constants $A_{s\alpha}$ and
$B_{s\alpha}$ only once, and then use these values for every spindown
parameter set used when searching in a given sky position.  Analytical
formulae for these constants are easily derived:
\begin{equation}
A_{s \alpha}=\frac{1}{s+1}\Delta T_{\alpha}^{s+1}-\frac{1}{2}\mathcal{T}_{SFT}\dot{T}_{\alpha}\Delta T_{\alpha}^{s}
\end{equation}
\begin{equation}
B_{s \alpha}=\mathcal{T}_{SFT}\dot{T}_{\alpha}\Delta T_{\alpha}^{s}
\end{equation}

</lalLaTeX> */

#ifndef _STACKSLIDE_H
#define _STACKSLIDE_H

#include <stdio.h>
#include <stdlib.h>
#include <lal/FileIO.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <math.h>
#include <string.h> 
#include <lal/LALConstants.h>
#include <lal/StreamInput.h>
#include <lal/SeqFactories.h>
#include <lal/LALBarycenter.h>
#include <lal/LALInitBarycenter.h>
#include <lal/Random.h>
#include <getopt.h>
#include <lal/Date.h>
#include <lal/Units.h>
#include <lal/LALDatatypes.h>
#include <lal/FindRoot.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/VectorOps.h>


#ifdef __cplusplus
extern "C" {
#endif
  
  NRCSID (STACKSLIDEH, "$Id$");
  
   
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
struct CSParams
\end{verbatim}
\index{\texttt{CSParams}}

\noindent This structure contains the parameters for the \verb@ComputeSky()@ routine.  The parameters are:

\begin{description}
\item[\texttt{INT8 spinDwnOrder}] The maximal number of spindown parameters per spindown parameter set.
\item[\texttt{INT8 mObsSFT}] The number of SFTs in the observation time.
\item[\texttt{REAL8 tSFT}] The timescale of one SFT.
\item[\texttt{LIGOTimeGPS *tGPS}] An array containing the GPS times of the first datum from each SFT.
\item[\texttt{REAL8 *skyPos}] The array containing the sky patch coordinates.
\item[\texttt{CHAR *sw}] A switch which turns modulation on/off. 
\item[\texttt{void (*funcName)(REAL8 , REAL8 , REAL8 , REAL8 *, REAL8 *, const CHAR *sw)}] A function pointer, to make the use of different timing routines easy.
\end{description}

</lalLaTeX> */

/* <lalErrTable file="ComputeSkyBinaryHErrorTable"> */
#define COMPUTESKYBINARYH_ENULL 1
#define COMPUTESKYBINARYH_ENNUL 2
#define COMPUTESKYBINARYH_ERANG 3
#define COMPUTESKYBINARYH_ENEGA 4
#define COMPUTESKYBINARYH_MSGENULL "Null Pointer"
#define COMPUTESKYBINARYH_MSGENNUL "Non-Null Pointer"
#define COMPUTESKYBINARYH_MSGERANG "Input parameter out of range"
#define COMPUTESKYBINARYH_MSGENEGA "Bad Negative Value"
/* </lalErrTable>  */

  
#define ACC 1e-9

typedef struct
tagStackSlideParams /* substituted tagStackSlideBinaryParams*/
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
	INT2 binaryFlag;
	REAL8 dfSUM;
	UINT4 gpsStartTimeSec;
	UINT4 gpsStartTimeNan;
	INT4 numSpinDown;
	EphemerisData *edat;
	LIGOTimeGPS *timeStamps;
	/* INT4 patchNumber; */
	BarycenterInput *baryinput;
        REAL8 		SemiMajorAxis;  /* orbital radius of binary (in sec) */
        REAL8           OrbitalPeriod;         /* Period of binary (in sec) */
        REAL8           OrbitalEccentricity;   /* Orbital eccentricy */
        REAL8           ArgPeriapse;    /* Argument of Periapse */
        UINT4 TperiapseSSBSec;
	UINT4 TperiapseSSBNanoSec;
	REAL8 deltaSMA;
	REAL8 SMAcentral;
        INT4 iFreqDeriv;
}
StackSlideParams; 


typedef struct
tagStackSlideSkyParams
{
	INT8		spinDwnOrder;	/* max spindown parameter order */
	INT8		mObsSFT;	/* number of coherent timescales */
	REAL8		tSFT;		/* timescale of SFT */
	LIGOTimeGPS	*tGPS;		/* GPS time of 1st data sample of each SFT */
	UINT4 gpsStartTimeSec;          /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
	UINT4 gpsStartTimeNan;          /* 06/05/04 gam; set these to epoch that gives T0 at SSB. */
	REAL8 		*skyPos; 	/* array of sky positions */
        REAL8 		SemiMajorAxis;  /* orbital radius of binary (in sec) */
        REAL8           OrbitalPeriod;         /* Period of binary (in sec) */
        REAL8           OrbitalEccentricity;   /* Orbital eccentricy */
        REAL8           ArgPeriapse;    /* Argument of Periapse */
        LIGOTimeGPS     TperiapseSSB;   /* Instance of periapse passage measured in the SSB frame */
	UINT4 gpsTperiapseSSBSec;
	BarycenterInput *baryinput;	
	EmissionTime *emit;
	EarthState *earth;
	EphemerisData *edat;
	/*REAL8 alphaS-X1;
	REAL8 deltaS-X1;*/
	REAL8 dInv;
       }
StackSlideSkyParams;

typedef struct
tagTdotsAndDeltaTs
{
	REAL8		*vecTDots;	/* 1-d array of (dT/dt)'s for frequency calculation */
	REAL8		**vecDeltaTs;	/* 2-d array of (T-T_0)'s for frequency calculation */
}
TdotsAndDeltaTs;


void StackSlideBinary(     LALStatus *status, 
		       StackSlideParams *stksldParams,
		       REAL4FrequencySeries **STKData,
		       REAL4FrequencySeries **SUMData); /*make sure about the argument for the function*/

void StackSlide(	LALStatus *status, 
			REAL4FrequencySeries **SUMData, 
			REAL4FrequencySeries **STKData,
			TdotsAndDeltaTs *pTdotsAndDeltaTs,
			StackSlideParams *params);


	
void StackSlideComputeSkyBinary (LALStatus 	*status, 
			TdotsAndDeltaTs 	*pTdotsAndDeltaTs, 
			INT8 		iSkyCoh, 
			StackSlideSkyParams 	*params);
		

void StackSlideComputeSky (LALStatus 	*status, 
			TdotsAndDeltaTs 	*pTdotsAndDeltaTs, 
			/*INT8 		iSkyCoh, */
			StackSlideSkyParams 	*params);




#ifdef __cplusplus
}
#endif

#endif /* _STACKSLIDE_H */
