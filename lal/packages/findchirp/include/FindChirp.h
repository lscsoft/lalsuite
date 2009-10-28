/*
*  Copyright (C) 2007 Sukanta Bose, Chad Hanna, Darren Woods, Diego Fazi, Drew Keppel, Duncan Brown, Eirini Messaritaki, Gareth Jones, Jolien Creighton, Patrick Brady, Anand Sengupta, Stephen Fairhurst, Craig Robinson , Sean Seader, Thomas Cokelaer
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirp.h
 *
 * Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpHV">
Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
$Id$
</lalVerbatim>

<lalLaTeX>
\section{Header \texttt{FindChirp.h}}
\label{s:FindChirp.h}

\noindent This header provides core prototypes, structures and functions to
filter interferometer data for binary inspiral chirps.

\subsubsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirp.h>
\end{verbatim}

\noindent Each function in findchirp falls into one of four classes:
\begin{enumerate}
\item Generate management functions which are independent of the type of
filtering implemented. The prototypes for these functions are provided by
this header file.

\item Functions for filtering data for time domain and frequency domain
templates with an unknown amplitude and phase. These are the functions
that implement matched filtering for time domain templates (TaylorT1,
TaylorT2, TaylorT2, PadeT1, EOB and GeneratePPN) and matched filtering
for post-Newtonian frequency domain templates (FindChirpSP). The main
filter function \texttt{FindChirpFilterSegment()} is prototyped in
this header file. The template generation and data conditioning functions
are prototyped in \texttt{FindChirpSP.h} and \texttt{FindChirpTD.h}.
Full documentation of the filtering algorithm used can be found in the
documentation of the module \texttt{FindChirpFilter.c}.

\item Functions to filter interferometer data for using the frequency
domain non-spinning black hole detection template family known as BCV.
These functions are protoyped by the header \texttt{FindChirpBCV.h}
which contains documentation of the algorithms used.

\item Functions to filter interferometer data for using the frequency domain
spinning black hole detection template family known as BCVSpin.  These
functions are protoyped by the header \texttt{FindChirpBCVSpin.h} which
contains documentation of the algorithms used.
\end{enumerate}

\noindent The goal of all the filtering functions is to determine if the
(calibrated) output of the interferometer $s(t)$ contains a gravitational wave
$h(t)$ in the presence of the detector noise $n(t)$. When the interferometer
is operating properly
\begin{equation}
s(t) = \left\{ \begin{array}{ll}
n(t) + h(t) & \textrm{signal present},\\
n(t) & \textrm{signal absent}.
\end{array}\right.
\end{equation}
The detection of signals of known form in noise is a classic problem of signal
processing\cite{wainstein:1962} and can be answered by the construction of a
\emph{detection statistic} and a test to see if the statistic is above some
pre-assigned threshold. The construction of the various detection
statistics used for each the three types of search are described in the modules
that implement the search.

</lalLaTeX>
#endif

#ifndef _FINDCHIRPH_H
#define _FINDCHIRPH_H

#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/FindChirpDatatypes.h>
#include <lal/FindChirpChisq.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPH, "$Id$");

#if 0
<lalLaTeX>
\subsection*{Error codes}
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPH_ENULL 1
#define FINDCHIRPH_ENNUL 2
#define FINDCHIRPH_EALOC 3
#define FINDCHIRPH_ENUMZ 5
#define FINDCHIRPH_ESEGZ 6
#define FINDCHIRPH_ECHIZ 7
#define FINDCHIRPH_EDTZO 8
#define FINDCHIRPH_ETRNC 10
#define FINDCHIRPH_EFLOW 11
#define FINDCHIRPH_EFREE 12
#define FINDCHIRPH_ERHOT 15
#define FINDCHIRPH_ECHIT 16
#define FINDCHIRPH_ECRUP 17
#define FINDCHIRPH_ESMSM 18
#define FINDCHIRPH_EHETR 19
#define FINDCHIRPH_EDFDT 20
#define FINDCHIRPH_EAPRX 21
#define FINDCHIRPH_EUAPX 22
#define FINDCHIRPH_ECHTZ 23
#define FINDCHIRPH_EMASS 24
#define FINDCHIRPH_EWVFM 25
#define FINDCHIRPH_EBCVC 25
#define FINDCHIRPH_EMAPX 26
#define FINDCHIRPH_EPTFW 27
#define FINDCHIRPH_EIGEN 28
#define FINDCHIRPH_EIMRW 29

#define FINDCHIRPH_MSGENULL "Null pointer"
#define FINDCHIRPH_MSGENNUL "Non-null pointer"
#define FINDCHIRPH_MSGEALOC "Memory allocation error"
#define FINDCHIRPH_MSGENUMZ "Invalid number of points in segment"
#define FINDCHIRPH_MSGESEGZ "Invalid number of segments"
#define FINDCHIRPH_MSGECHIZ "Invalid number of chi squared bins"
#define FINDCHIRPH_MSGEDTZO "deltaT is zero or negative"
#define FINDCHIRPH_MSGETRNC "Duration of inverse spectrum in time domain is negative"
#define FINDCHIRPH_MSGEFLOW "Inverse spectrum low frequency cutoff is negative"
#define FINDCHIRPH_MSGEFREE "Error freeing memory"
#define FINDCHIRPH_MSGERHOT "Rhosq threshold is negative"
#define FINDCHIRPH_MSGECHIT "Chisq threshold is negative"
#define FINDCHIRPH_MSGECRUP "Chirp length or invSpecTrunc too long for length of data segment"
#define FINDCHIRPH_MSGESMSM "Size mismatch between vectors"
#define FINDCHIRPH_MSGEHETR "Attempting to simulate heterodyned GW"
#define FINDCHIRPH_MSGEDFDT "Waveform sampling interval is too large"
#define FINDCHIRPH_MSGEAPRX "Incorrect waveform approximant"
#define FINDCHIRPH_MSGEUAPX "Unknown waveform approximant"
#define FINDCHIRPH_MSGECHTZ "Length of chirp is zero or negative"
#define FINDCHIRPH_MSGEMASS "Invalid mass parameters for template generation"
#define FINDCHIRPH_MSGEWVFM "Unknown injection waveform"
#define FINDCHIRPH_MSGEBCVC "BCVC code: thetav not in [-pi, pi]."
#define FINDCHIRPH_MSGEMAPX "Mismatch in waveform approximant"
#define FINDCHIRPH_MSGEPTFW "Error generating PTF waveform"
#define FINDCHIRPH_MSGEIGEN "Error computing eigenvalues"
#define FINDCHIRPH_MSGEIMRW "Error computing IMR waveform"

/* </lalErrTable> */


#if 0
<lalLaTeX>
\subsection*{Types}

\subsubsection*{Structure \texttt{FindChirpInitParams}}
\idx[Type]{FindChirpInitParams}

\noindent This structure provides the essential information for the
filter initialisation and memory allocation functions used by findchirp.

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpInitParams
{
  UINT4                         numSegments;
  UINT4                         numPoints;
  UINT4                         ovrlap;
  UINT4                         numChisqBins;
  BOOLEAN                       createRhosqVec;
  BOOLEAN                       createCVec;
  Approximant                   approximant;
  LALPNOrder                    order;
}
FindChirpInitParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>

\begin{description}
\item[\texttt{UINT4 numSegments}] The number of data segments in the input
\texttt{DataSegmentVector} and a the \texttt{FindChirpSegmentVector}.

\item[\texttt{UINT4 numPoints}] The number of discrete data points $N$ in each
data segment.

\item[\texttt{UINT4 ovrlap}] The number of sample points by which each
data segment overlaps.

\item[\texttt{UINT4 numChisqBins}] The number of bins $p$ used to contruct the
$\chi^2$ veto.

\item[\texttt{BOOLEAN createRhosqVec}] Flag that controls whether or not the
filter function should store the output of the matched filter, $\rho^2(t)$, as
well as the events. Memory is allocated for this vector if the flag is set to
1.

\item[\texttt{BOOLEAN createCVec}] Flag that controls whether or not the
filter function should store the complex filter output $x(t) + i y(t)$ needed
by the coherent inspiral code.  Memory is allocated for this vector if the
flag is set to 1.

\item[\texttt{Approximant approximant}] Initialize the findchirp routines
to fiter with templates of type \texttt{approximant}. Valid approximants are
TaylorT1, TaylorT2, TaylorT3, PadeT1, EOB, FindChirpSP, BCV and BCVSpin.

\end{description}

\subsubsection*{Structure \texttt{FindChirpDataParams}}
\idx[Type]{FindChirpDataParams}

\noindent This structure contains the parameters needed to call the data
conditioning functions \texttt{FindChirpSPData()}, \texttt{FindChirpTDData()},
\texttt{FindChirpBCVData()} or \texttt{FindChirpBCVSpinData()}. It should be
initialized by \texttt{FindChirpDataInit()} and destroyed by
\texttt{FindChirpDataFinalize()}.

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpDataParams
{
  REAL4Vector                  *ampVec;
  REAL4Vector                  *ampVecBCV;
  REAL8Vector                  *ampVecBCVSpin1;
  REAL8Vector                  *ampVecBCVSpin2;
  RealFFTPlan                  *fwdPlan;
  RealFFTPlan                  *invPlan;
  REAL4Vector                  *wVec;
  COMPLEX8Vector               *wtildeVec;
  REAL4Vector                  *tmpltPowerVec;
  REAL4Vector                  *tmpltPowerVecBCV;
  REAL4                         fLow;
  REAL4                         dynRange;
  UINT4                         invSpecTrunc;
  Approximant                   approximant;
}
FindChirpDataParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>

\begin{description}
\item[\texttt{REAL4Vector *ampVec}] A vector containing the frequency domain
quantity $(k/N)^{-7/6}$, where $k$ is the frequency series index and $N$ is the
number of points in a data segment. NB: for time domain templates, this is set
to unity by the function \texttt{FindChirpTDData()}.

\item[\texttt{REAL4Vector *ampVecBCV}] A vector containing the frequency domain
quantity $(k/N)^{-1/2}$, where $k$ is the frequency series index and $N$ is the
number of points in a data segment.

\item[\texttt{REAL4Vector *ampVecBCVSpin1}] Undocumented spinning BCV
amplitude vector.

\item[\texttt{REAL4Vector *ampVecBCVSpin2}] Undocumented spinning BCV
amplitude vector.

\item[\texttt{REAL4Vector *fwdPlan}] An FFTW plan used to transform the
time domain interferometer data $v(t_j)$ into its DFT $\tilde{v}_k$.

\item[\texttt{REAL4Vector *invPlan}] An FFTW plan used to transform the
dimensionless frequency domain interferometer strain $\tilde{w}_k$ into
the quantity $N w(t_j)$ to allow time domain trunction of the inverse
power spectrum.

\item[\texttt{REAL4Vector *wVec}] A vector used as workspace when truncating
the imverse power spectrum in the time domain.

\item[\texttt{COMPLEX8Vector *wtildeVec}] A vector which on exit from
the data conditioning function contains the inverse of the strain one sided
power spectral density, after trunction in the time domain, \emph{for the last
data segment conditioned.} Typically all the data segments are conditioned
using the same power spectrum, so this quantity is identical for all data
segments. It contains:
\begin{equation}
\tilde{w}_k = {1}/{\ospsd}.
\end{equation}

\item[\texttt{REAL4Vector *tmpltPowerVec}] A vector which on exit from
\texttt{FindChirpSPData()} or from \texttt{FindChirpBCVData()}
contains the quantity
\begin{equation}
\mathtt{tmpltPower[k]} = \frac{f^{-7/3}}{\ospsd}
\end{equation}

\item[\texttt{REAL4Vector *tmpltPowerVecBCV}] A vector which on exit from
\texttt{FindChirpBCVData()}
contains the quantity
\begin{equation}
\mathtt{tmpltPowerBCV[k]} = \frac{f^{-1}}{\ospsd}
\end{equation}

\item[\texttt{REAL4 fLow}] The frequency domain low frequency cutoff
$f_\mathrm{low}$. All frequency domain data is set to zero below this
frequency.

\item[\texttt{REAL4 dynRange}] A dynamic range factor $d$ which cancels from
the filter output.  This allows quantities to be stored in the range of
\texttt{REAL4} rather than \texttt{REAL8}. This must be set to the same value
as \texttt{dynRange} in the \texttt{FindChirpTmpltParams}. For LIGO data a
value of $d = 2^{69}$ is appropriate.

\item[\texttt{UINT4 invSpecTrunc}] The length to which to truncate the inverse
power spectral density of the data in the time domain. If set to zero, no
truncation is performed.

\item[\texttt{Approximant approximant}] Condition the data for templates of
type \texttt{approximant}. Valid approximants are TaylorT1, TaylorT2,
TaylorT3, PadeT1, EOB, FindChirpSP, BCV and BCVSpin.

\end{description}

\subsubsection*{Structure \texttt{FindChirpTmpltParams}}
\idx[Type]{FindChirpTmpltParams}

\noindent This structure contains the parameters for generation of templates
by the various template generation functions provided in finchirp.

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpTmpltParams
{
  REAL8                         deltaT;
  REAL4                         fLow;
  REAL4                         dynRange;
  REAL4Vector                  *xfacVec;
  REAL4VectorSequence          *ACTDVecs;
  REAL4VectorSequence          *PTFQ;
  REAL4Vector                  *PTFphi;
  REAL4Vector                  *PTFomega_2_3;
  REAL4VectorSequence          *PTFe1;
  REAL4VectorSequence          *PTFe2;
  RealFFTPlan                  *fwdPlan;
  Approximant                   approximant;
  LALPNOrder                    order;
  INT4                          reverseChirpBank;
  INT4                          bandPassTmplt;
  InspiralApplyTaper            taperTmplt;
}
FindChirpTmpltParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>

\begin{description}
\item[\texttt{REAL8 deltaT}] The sampling interval $\Delta t$ of the input
data channel.

\item[\texttt{REAL4 fLow}] The frequency domain low frequency cutoff
$f_\mathrm{low}$. All frequency domain data is zero below this frequency.

\item[\texttt{REAL4 dynRange}] A dynamic range factor $d$ which cancels from
the filter output.  This allows quantities to be stored in the range of
\texttt{REAL4} rather than \texttt{REAL8}. This must be set to the same value
as \texttt{dynRange} in the \texttt{FindChirpDataParams}. For LIGO data a
value of $d = 2^{69}$ is appropriate.

\item[\texttt{REAL4Vector *xfacVec}] For frequency domain templates, this is a
vector of length $N/2+1$ which contains the quantity $k^{-1/3}$. For time
domain templates, this is a workspace vector of length $N$ which contains the
time domain template generated by the inspiral package, shifted so that the
end of the template is at the end of the vector. This vector is Fourier
transformed to obtain the quantity findchirp template $\tilde{T}_k$.

\item[\texttt{REAL4Vector *fwdPlan}] For time domain templates, an FFTW plan
used to transform the time domain data stored in \texttt{xfacVec} into its DFT
which is stored in the findchirp template.

\item[\texttt{Approximant approximant}] Generate templates of type
\texttt{approximant}. Valid approximants are TaylorT1, TaylorT2, TaylorT3,
PadeT1, EOB, FindChirpSP, BCV and BCVSpin. For time domain templates the
post-Newtonian order is always two. For stationary phase templates, the
post-Newtonian order is specified by \texttt{order}.

\item[\texttt{LALPNOrder order}] Specifies the post-Newtonian order of the
templates. Valid pN orders are twoPN, twoPointFivePN, threePN,
threePointFivePN, pseudoFourPN. pseudoFourPN is not the true four PN
correction, but may increase the fitting factor between stationary phase and
numerical relativity waveforms.

\item[\texttt{INT4 reverseChirpBank}] Switches a FindChirpSP template bank
to be a reverse chirp template bank if true.
\end{description}

\subsubsection*{Structure \texttt{FindChirpClustering}}
\idx[Type]{FindChirpClustering}

\noindent This structure contains the possible methods by which
to maximize over a chirp in a data segment.
</lalLaTeX>
#endif

/* <lalVerbatim> */
typedef enum {
   FindChirpClustering_none,
   FindChirpClustering_tmplt,
   FindChirpClustering_window,
   FindChirpClustering_tmpltwindow
}
FindChirpClustering;
/* </lalVerbatim> */

#if 0
<lalLaTeX>

\begin{description}
\item[\texttt{FindChirpClustering\_none}]
The decision to do no clustering
of events.

\item[\texttt{FindChirpClustering\_tmplt}]
Cluster over the length of the data segment.

\item[\texttt{FindChirpClustering\_window}]
Cluster over a given number of seconds
given by the argument to the flag \texttt{--cluster-window}
(required to be less than the length of the data segment).

\end{description}

\subsubsection*{Structure \texttt{FindChirpFilterOutputVetoParams}}
\idx[Type]{FindChirpFilterOutputVetoParams}

\noindent This structure provides the parameters for the filter
output veto.

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpFilterOutputVetoParams
{
  REAL4          rsqvetoWindow;
  REAL4          rsqvetoThresh;
  REAL4          rsqvetoTimeThresh;
  REAL4          rsqvetoMaxSNR;
  REAL4          rsqvetoCoeff;
  REAL4          rsqvetoPow;
}
FindChirpFilterOutputVetoParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>

\begin{description}
\item[\texttt{REAL4 rsqvetoWindow}] Width of the $r^2$ veto window in units of
seconds.

\item[\texttt{REAL4 rsqvetoThresh}] Threshold of the $r^2$ veto test analogous to the
$r^2$ threshold employed in the bns and macho inspiral searches.
\end{description}
\subsubsection*{Structure \texttt{FindChirpFilterParams}}
\idx[Type]{FindChirpFilterParams}

\noindent This structure provides the parameters used by the
\texttt{FindChirpFilterSegment()} function.

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpFilterParams
{
  REAL8                         deltaT;
  REAL4                         clusterWindow;
  REAL4                         rhosqThresh;
  REAL4                         chisqThresh;
  REAL4                         chisqDelta;
  UINT4                         maximiseOverChirp;
  UINT4                         ignoreIndex;
  FindChirpClustering           clusterMethod;
  Approximant                   approximant;
  LALPNOrder                    order;
  COMPLEX8Vector               *qVec;
  COMPLEX8Vector               *qVecBCV;
  COMPLEX8Vector               *qVecBCVSpin1;
  COMPLEX8Vector               *qVecBCVSpin2;
  COMPLEX8Vector               *qtildeVec;
  COMPLEX8Vector               *qtildeVecBCV;
  COMPLEX8Vector               *qtildeVecBCVSpin1;
  COMPLEX8Vector               *qtildeVecBCVSpin2;
  COMPLEX8Vector              **qVecACTD;
  COMPLEX8Vector              **qtildeVecACTD;
  COMPLEX8VectorSequence       *PTFqVec;
  COMPLEX8Vector               *PTFsnrVec;
  REAL4Array                   *PTFA;
  REAL4Array                   *PTFMatrix;
  ComplexFFTPlan               *invPlan;
  REAL4TimeSeries              *rhosqVec;
  COMPLEX8TimeSeries           *cVec;
  REAL4Vector                  *chisqVec;
  FindChirpChisqParams         *chisqParams;
  FindChirpChisqInput          *chisqInput;
  FindChirpChisqInput          *chisqInputBCV;
  FindChirpFilterOutputVetoParams      *filterOutputVetoParams;
}
FindChirpFilterParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>

\begin{description}
\item[\texttt{REAL8 deltaT}] The sampling interval $\Delta t$.

\item[\texttt{REAL4 rhosqThresh}] The signal-to-noise ratio squared threshold
$\rho^2_\ast$. If the matched filter output exceeds this value, that is
$\rho^2(t_j) > \rho^2_\ast$, the event processing algorithm is entered and
triggers may be generated (subject to addition vetoes such as the $\chi^2$
veto). The value of $\rho^2_\ast0$ must be greater than or equal to zero.

\item[\texttt{REAL4 chisqThresh}] The $\chi^2$ veto threshold on. This
threshold is described in details in the documentation for the $\chi^2$
veto.

\item[\texttt{REAL4 norm}] On exit this contains the normalisation constant
that relates the quantity $|q_j|^2$ with the signal to noise squared,
$\rho^2(t_j)$ by
\begin{equation}
\rho^2(t_j) = \textrm{norm} \times \left|q_j\right|^2.
\end{equation}

\item[\texttt{UINT4 maximiseOverChirp}] If not zero, use the
maximise over chirp length algorithm to decide which time $t_j$ should
have an inspiral trigger generated. Otherwise record all points that pass the
$\rho^2$ and $\chi^2$ threshold as triggers (this may generate may triggers).

\item[\texttt{Approximant approximant}] Filter the data using templates of
type \texttt{approximant}. Valid approximants are TaylorT1, TaylorT2,
TaylorT3, PadeT1, EOB, FindChirpSP, BCV and BCVSpin. The value of
\texttt{approximant} here must match that in the findchirp data segment and
findchirp template used as input.

\item[\texttt{COMPLEX8Vector *qVec}] Pointer to vector of length $N$ allocated
by \texttt{FindChirpFilterInit()} to store the quantity $q_j$. The pointer
must not be NULL on entry, but the vetor may contain garbage which will be
overwritten with the value of $q_j$ for the segment filtered on exit.

\item[\texttt{COMPLEX8Vector *qVecBCV}] Pointer to the additional vector
required for the BCV templates, allocated by \texttt{FindChirpFilterInit()}.

\item[\texttt{COMPLEX8Vector *qVecBCVSpin1}] Pointer to the additional vector
required for filtering spinning BCV templates, allocated by
\texttt{FindChirpFilterInit()}.

\item[\texttt{COMPLEX8Vector *qVecBCVSpin2}] Pointer to the additional vector
required for filtering spinning BCV templates, allocated by
\texttt{FindChirpFilterInit()}.

\item[\texttt{COMPLEX8Vector *qtildeVec}] Pointer to vector of length $N$
allocated by \texttt{FindChirpFilterInit()} to store the quantity
$\tilde{q}_k$, given by
\begin{equation}
\tilde{q}_k = \left\{
\begin{array}{ll}
\tilde{F}_k \tilde{T}_k^\ast & \quad 0 < k < \frac{N}{2} \\,
0 & \quad \textrm{otherwise}.
\end{array}
\right.
\end{equation}
The pointer must not be NULL on entry, but the vetor may contain garbage which
will be overwritten with the value of $\tilde{q}_k$ for the segment filtered
on exit.

\item[\texttt{COMPLEX8Vector *qtildeVecBCV}] Pointer to the additional
vector required for filtering BCV templates, allocated by
\texttt{FindChirpFilterInit()}.

\item[\texttt{COMPLEX8Vector *qtildeVecBCVSpin1}] Pointer to the additional
vector required for filtering spinning BCV templates, allocated by
\texttt{FindChirpFilterInit()}.

\item[\texttt{COMPLEX8Vector *qtildeVecBCVSpin2}] Pointer to the additional
vector required for filtering spinning BCV templates, allocated by
\texttt{FindChirpFilterInit()}.

\item[\texttt{ComplexFFTPlan *invPlan}] Pointer to FFTW plan created by
\texttt{FindChirpFilterInit()} to transform the quantity $\tilde{q}_k$ to
${q}_j$ usimg the inverse DFT. Must not be NULL.

\item[\texttt{REAL4TimeSeries *rhosqVec}] Pointer to a time series which
contains a vector of length $N$. If this is not NULL, the filter output $\rho^2(t_j)$
is stored in the vector.

\item[\texttt{COMPLEX8Vector *rhosqVec}] Pointer to a time series which
contains a vector of length $N$. If this is not NULL, the complex filter
output $\rho(t_j) = x(t_j) + iy(t_j)$ is stored in the vector. This
quantity can be used by the coherent filtering code.

\item[\texttt{REAL4Vector *chisqVec}] Workspace vector of length $N$ used to
compute and store $\chi^2(t_j)$. Must not be NULL if \texttt{numChisqBins} is
greater than zero. Contains $\chi^2(t_j)$ on exit.

\item[\texttt{FindChirpChisqParams *chisqParams}] Pointer to parameter
structure for the $\chi^2$ veto. Must not be NULL if \texttt{numChisqBins} is
greater than zero.

\item[\texttt{FindChirpChisqInput *chisqInput}] Pointer to input data
structure for the $\chi^2$ veto.  Must not be NULL if \texttt{numChisqBins} is
greater than zero.

\item[\texttt{FindChirpChisqInput *chisqInputBCV}] Pointer to input data
structure for the BCV $\chi^2$ veto. Must not be NULL if the approximant is
BCV and \texttt{numChisqBins} is greater than zero.

\item[\texttt{FindChirpFilterOutputVetoParams *filterOutputVetoParams}]
Pointer to the parameter structure for the additional signal based veto
function.

\end{description}
</lalLaTeX>
#endif


/*
 *
 * typedefs of input structures used by functions in findchirp
 *
 */


#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpFilterInput}}
\idx[Type]{FindChirpFilterInput}

\noindent This structure groups the input data required for the
\texttt{FindChirpFilterSegment()} function into a single structure.

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpFilterInput
{
  FindChirpTemplate            *fcTmplt;
  FindChirpSegment             *segment;
}
FindChirpFilterInput;
/* </lalVerbatim> */
#if 0
<lalLaTeX>

\begin{description}
\item[\texttt{FindChirpTemplate *fcTmplt}] Pointer to the input template
in a form that can be used by \texttt{FindChirpFilterSegment()}

\item[\texttt{FindChirpSegment *segment}] Pointer to the input data segment
in a form that can be used by \texttt{FindChirpFilterSegment()}
\end{description}

\subsubsection*{Structure \texttt{FindChirpSubBankData}}
\idx[Type]{FindChirpSubBankData}

\noindent This structure contains data needed for the bank veto.
</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpBankVetoData
{
  UINT4                   length;
  COMPLEX8Vector        **qtildeVecArray;
  COMPLEX8Vector        **qVecArray;
  FindChirpFilterInput  **fcInputArray;
  REAL4Vector            *ccMat;
  REAL4Vector		 *normMat;
  REAL4Vector		 *spec;
  COMPLEX8Vector         *resp;
}
FindChirpBankVetoData;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\begin{description}

\item[\texttt{struct tagFindChirpSubBankData *next}] The next structure in
the linked list.
\end{description}

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpBankSimParams
{
  Approximant           approx;
  LALPNOrder            order;
  REAL4                 minMass;
  REAL4                 maxMass;
  RandomParams         *randParams;
  INT4                  maxMatch;
  CHAR                 *frameName;
  CHAR                 *frameChan;
  REAL4			f_lower;
}
FindChirpBankSimParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>

\begin{description}
\item[\texttt{Approximant approx}] Waveform pproximant to use for injection.

\item[\texttt{Approximant order}] Waveform order to use for injection.

\item[\texttt{REAL4 minMass}] Minimum mass of injected signals.

\item[\texttt{REAL4 maxMass}] Maximum mass of injected signals.
\end{description}


\vfill{\footnotesize\input{FindChirpHV}}
</lalLaTeX>
#endif


/*
 *
 * function prototypes for memory management functions
 *
 */


#if 0
<lalLaTeX>
\newpage\input{FindChirpLinkedListC}
</lalLaTeX>
#endif

void
LALFindChirpCreateTmpltNode (
    LALStatus                  *status,
    InspiralTemplate           *thistmplt,
    InspiralTemplateNode      **tmpltNode
    );

void
LALFindChirpDestroyTmpltNode (
    LALStatus                  *status,
    InspiralTemplateNode      **tmpltNode
    );

#if 0
<lalLaTeX>
\newpage\input{FindChirpMemoryC}
</lalLaTeX>
#endif

void
LALInitializeDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **dataSegVec,
    REAL4TimeSeries            *chan,
    REAL4FrequencySeries       *spec,
    COMPLEX8FrequencySeries    *resp,
    FindChirpInitParams        *params
    );

void
LALFinalizeDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector
    );

void
LALCreateDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector,
    FindChirpInitParams        *params
    );

void
LALDestroyDataSegmentVector (
    LALStatus                  *status,
    DataSegmentVector         **vector
    );

void
LALCreateFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector,
    FindChirpInitParams        *params
    );

void
LALDestroyFindChirpSegmentVector (
    LALStatus                  *status,
    FindChirpSegmentVector    **vector
    );


/*
 *
 * function prototypes for initialization, finalization of data functions
 *
 */


#if 0
<lalLaTeX>
\newpage\input{FindChirpDataC}
</lalLaTeX>
#endif


void
LALFindChirpDataInit (
    LALStatus                  *status,
    FindChirpDataParams       **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpDataFinalize (
    LALStatus                  *status,
    FindChirpDataParams       **output
    );


/*
 *
 * function prototypes for initialization, finalization of template functions
 *
 */


#if 0
<lalLaTeX>
\newpage\input{FindChirpTemplateC}
</lalLaTeX>
#endif

void
LALFindChirpTemplateInit (
    LALStatus                  *status,
    FindChirpTmpltParams      **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpTemplateFinalize (
    LALStatus                  *status,
    FindChirpTmpltParams      **output
    );



/*
 *
 * function prototypes for initialization, finalization and filter functions
 *
 */


#if 0
<lalLaTeX>
\newpage\input{FindChirpFilterC}
</lalLaTeX>
#endif

void
LALFindChirpFilterInit (
    LALStatus                  *status,
    FindChirpFilterParams     **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpFilterFinalize (
    LALStatus                  *status,
    FindChirpFilterParams     **output
    );

void
LALCreateFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output,
    FindChirpInitParams        *params
    );

void
LALDestroyFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output
    );

void
LALFindChirpCreateCoherentInput(
     LALStatus                  *status,
     COMPLEX8TimeSeries         **coherentInputData,
     COMPLEX8TimeSeries         *input,
     SnglInspiralTable          *templt,
     REAL4                      coherentSegmentLength,
     INT4                       corruptedDataLength
     );

void
LALFindChirpFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    );

void
LALFindChirpStoreEvent (
    LALStatus                  *status,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    SnglInspiralTable          *thisEvent,
    COMPLEX8                   *q,
    UINT4                       kmax,
    REAL4                       norm,
    UINT4                       eventStartIdx,
    UINT4                       numChisqBins,
    CHAR                       *searchName
    );

void
LALFindChirpClusterEvents (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    FindChirpBankVetoData      *bankVetoData,
    UINT4                       subBankIndex
    );

#if 0
<lalLaTeX>
\newpage\input{FindChirpFilterOutputVetoC}
</lalLaTeX>
#endif

void
LALFindChirpFilterOutputVeto(
    LALStatus                          *status,
    SnglInspiralTable                 **eventList,
    FindChirpFilterInput               *input,
    FindChirpFilterParams              *fcParams
    );

#if 0
<lalLaTeX>
\newpage\input{FindChirpSimulationC}
</lalLaTeX>
#endif

void
LALFindChirpInjectSignals (
    LALStatus                  *status,
    REAL4TimeSeries            *chan,
    SimInspiralTable           *events,
    COMPLEX8FrequencySeries    *resp
    );

void
LALFindChirpInjectIMR (
    LALStatus                     *status,
    REAL4TimeSeries               *chan,
    SimInspiralTable              *events,
    SimRingdownTable              *ringdownevents,
    COMPLEX8FrequencySeries       *resp,
    INT4                           injectSignalType
    );

INT4
XLALFindChirpTagTemplateAndSegment (
    DataSegmentVector       *dataSegVec,
    InspiralTemplate        *tmpltHead,
    SnglInspiralTable       **events,
    CHAR                    *ifo,
    REAL4                   tdFast,
    UINT4                   *analyseThisTmplt
    );


INT4
XLALFindChirpSetAnalyzeSegment (
    DataSegmentVector          *dataSegVec,
    SimInspiralTable           *injections
    );

INT4
XLALFindChirpSetFollowUpSegment (
    DataSegmentVector          *dataSegVec,
    SnglInspiralTable          **events
    );

void
LALFindChirpSetAnalyseTemplate (
    LALStatus                    *status,
    UINT4                        *analyseThisTmplt,
    REAL4                        mmFast,
    REAL8                        deltaF,
    INT4                         sampleRate,
    FindChirpDataParams          *fcDataParams,
    int                          numTmplts,
    InspiralTemplate             *tmpltHead,
    int                          numInjections,
    SimInspiralTable             *injections
    );

UINT4
XLALCmprSgmntTmpltFlags (
    UINT4 numInjections,
    UINT4 TmpltFlag,
    UINT4 SgmntFlag
    );

UINT4
XLALFindChirpBankSimInitialize (
    REAL4FrequencySeries       *spec,
    COMPLEX8FrequencySeries    *resp,
    REAL8                       fLow
    );

SimInspiralTable *
XLALFindChirpBankSimInjectSignal (
    DataSegmentVector          *dataSegVec,
    COMPLEX8FrequencySeries    *resp,
    SimInspiralTable           *injParams,
    FindChirpBankSimParams     *simParams
    );

REAL4
XLALFindChirpBankSimSignalNorm(
    FindChirpDataParams         *fcDataParams,
    FindChirpSegmentVector      *fcSegVec,
    UINT4                        cut
    );

SimInstParamsTable *
XLALFindChirpBankSimMaxMatch (
    SnglInspiralTable         **bestTmplt,
    REAL4                       matchNorm
    );

SimInstParamsTable *
XLALFindChirpBankSimComputeMatch (
    SnglInspiralTable   *inputTmplt,
    REAL4                matchNorm
    );

FindChirpSubBank*
XLALFindChirpCreateSubBanks(
    UINT4                      *maxSubBankSize,
    UINT4                       subBankSize,
    UINT4                       bankSize,
    InspiralTemplate           *bankHead
    );

void
XLALBankVetoCCMat (
    FindChirpBankVetoData 	*bankVetoData,
    FindChirpSubBank            *vetoBank,
    FindChirpDataParams         *params,
    REAL4 			dynRange,
    REAL4 			fLow,
    REAL4 			deltaF,
    REAL4			deltaT
    );

REAL4
XLALComputeBankVeto( FindChirpBankVetoData *bankVetoData,
                     UINT4 i,
                     UINT4 snrIX,
                     UINT4 *dof);

REAL4
XLALComputeFullChisq(
    FindChirpBankVetoData      *bankVetoData,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params,
    COMPLEX8                   *q,
    UINT4                       i,
    UINT4                       snrIX,
    UINT4                      *dof,
    REAL4                       norm
);

InspiralTemplate *
XLALFindChirpSortTemplates(
  InspiralTemplate *bankHead,
  UINT4 num,
  UINT4 subbanksize
);

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPH_H */
