/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpSP.h
 *
 * Author: Brown, D. A., BCV-Modifications by Messaritaki E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpSPHV">
Author: Brown, D. A., BCV-Modifications by Messaritaki E.
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{FindChirpSP.h}}
\label{s:FindChirpSP.h}

Provides structures and functions to condition interferometer data
and generate binary inspiral chirps using the stationary phase approximation.
Recent addition deals with the conditioning of the data using the BCV
templates.

\subsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpSP.h>
\end{verbatim}

\input{FindChirpSPHDoc}
</lalLaTeX>
#endif


#ifndef _FINDCHIRPSPH_H
#define _FINDCHIRPSPH_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPSPH, "$Id$");

#if 0
<lalLaTeX> 
\subsection*{Error codes} 
</lalLaTeX>
#endif
/* <lalErrTable> */
#define FINDCHIRPSPH_ENULL 1
#define FINDCHIRPSPH_ENNUL 2
#define FINDCHIRPSPH_EALOC 3
#define FINDCHIRPSPH_ENUMZ 4
#define FINDCHIRPSPH_ESEGZ 5
#define FINDCHIRPSPH_EMISM 6
#define FINDCHIRPSPH_EDELT 7
#define FINDCHIRPSPH_EFLOW 8
#define FINDCHIRPSPH_EDYNR 9
#define FINDCHIRPSPH_EISTN 10
#define FINDCHIRPSPH_EDIVZ 11
#define FINDCHIRPSPH_MSGENULL "Null pointer"
#define FINDCHIRPSPH_MSGENNUL "Non-null pointer"
#define FINDCHIRPSPH_MSGEALOC "Memory allocation error"
#define FINDCHIRPSPH_MSGENUMZ "Invalid number of segments"
#define FINDCHIRPSPH_MSGESEGZ "Invalid number of points in segments"
#define FINDCHIRPSPH_MSGEMISM "Mismatch between number of points in segments"
#define FINDCHIRPSPH_MSGEDELT "deltaT is zero or negative"
#define FINDCHIRPSPH_MSGEFLOW "Low frequency cutoff is negative"
#define FINDCHIRPSPH_MSGEDYNR "Dynamic range scaling is zero or negative"
#define FINDCHIRPSPH_MSGEISTN "Truncation of inverse power spectrum is negative"
#define FINDCHIRPSPH_MSGEDIVZ "Attempting to divide by zero"
/* </lalErrTable> */


#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif

/* --- the parameter structure for the data conditioning function -------- */
/* <lalVerbatim file="FindChirpSPHFindChirpSPDataParams"> */
typedef struct
tagFindChirpSPDataParams
{
  REAL4Vector                  *ampVec;
  REAL4Vector                  *ampVecBCV;
  RealFFTPlan                  *fwdPlan;
  RealFFTPlan                  *invPlan;
  REAL4Vector                  *wVec;
  COMPLEX8Vector               *wtildeVec;
  REAL4Vector                  *tmpltPowerVec;
  REAL4Vector                  *tmpltPowerVecBCV;
  REAL4                         deltaT;
  REAL4                         fLow;
  REAL4                         dynRange;
  UINT4                         invSpecTrunc;
}
FindChirpSPDataParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpSPDataParams}}
\idx[Type]{FindChirpSPDataParams}

\input{FindChirpSPHFindChirpSPDataParams}

\noindent This structure contains the parameters needed to call the
\texttt{FindChirpSPData()} function. It should be initialized by
\texttt{FindChirpSPDataInit()} and destroyed by
\texttt{FindChirpSPDataFinalize()}. The fields are:

\begin{description}
\item[\texttt{REAL4Vector *ampVec}] A vector containing the frequency domain
quantity $(k/N)^{-7/6}$, where $k$ is the frequency series index and $N$ is the
number of points in a data segment.

\item[\texttt{REAL4Vector *ampVecBCV}] A vector containing the frequency domain
quantity $(k/N)^{-1/2}$, where $k$ is the frequency series index and $N$ is the
number of points in a data segment.

\item[\texttt{REAL4Vector *fwdPlan}] An FFTW plan used to transform the
time domain interferometer data $v(t_j)$ into its DFT $\tilde{v}_k$.

\item[\texttt{REAL4Vector *fwdPlan}] An FFTW plan used to transform the
dimensionless frequency domain interferometer strain $\tilde{w}_k$ into 
the quantity $N w(t_j)$ to allow time domain trunction of the inverse 
power spectrum.

\item[\texttt{REAL4Vector *vVec}] {\color{red} FIXME} A vector to contain
the time domain interferometer output $v(t_j)$. This is obsolete since LIGO
gives us $v(t_j)$ as floats. The 40m prototype gave integers which needed to
be cast to floats.

\item[\texttt{REAL4Vector *wVec}] A vector used as workspace when truncating
the imverse power spectrum in the time domain.

\item[\texttt{COMPLEX8Vector *wtildeVec}] A vector which on exit from
\texttt{FindChirpSPData()} contains the inverse of the strain one sided power
spectral density, after trunction in the time domain, that is
$ \tilde{w}_k = {1}/{\ospsd}$.

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


\item[\texttt{REAL4 deltaT}] {\color{red} FIXME} The sampling interval 
$\Delta t$. Should be a \texttt{REAL8} or derived from the input time series
\texttt{chan}.

\item[\texttt{REAL4 fLow}] The low frequency cutoff for the algorithm. All
data is zero below this frequency.

\item[\texttt{REAL4 dynRange}] A dynamic range factor which cancells from
the filter output (if set correctly in \texttt{FindChirpSPTmplt()} as well).
This allows quantities to be stored in the range of \texttt{REAL4} rather
than \texttt{REAL8}.

\item[\texttt{UINT4 invSpecTrunc}] The length to which to truncate the inverse
power spectral density of the data in the time domain. If set to zero, no
truncation is performed.
\end{description}
</lalLaTeX>
#endif

/* --- vector of DataSegment, as defined the framedata package ----------- */
/* <lalVerbatim file="FindChirpSPHFindChirpSPTmpltParams"> */
typedef struct
tagFindChirpSPTmpltParams
{
  REAL4                         deltaT;
  REAL4                         fLow;
  REAL4                         dynRange;
  REAL4Vector                  *xfacVec;
}
FindChirpSPTmpltParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpSPTmpltParams}}
\idx[Type]{FindChirpSPTmpltParams}

\input{FindChirpSPHFindChirpSPTmpltParams}

\noindent This structure contains the parameters for generation of
stationary phase templates by the function \texttt{FindChirpSPTmplt()}
It should be initialized by \texttt{FindChirpSPTmpltInit()} and destroyed by
\texttt{FindChirpSPTmpltFinalize()}. The fields are:

\begin{description}
\item[\texttt{REAL4 *deltaT}] {\color{red} FIXME} The sampling interval 
$\Delta t$. Should be a \texttt{REAL8}.

\item[\texttt{REAL4 fLow}] The low frequency cutoff for the algorithm. All
data is zero below this frequency.

\item[\texttt{REAL4 dynRange}] A dynamic range factor which cancells from
the filter output (if set correctly in \texttt{FindChirpSPData()} as well).
This allows quantities to be stored in the range of \texttt{REAL4} rather
than \texttt{REAL8}.

\item[\texttt{REAL4Vector *xfacVec}] A vector containing the frequency 
domain quantity $k^{-7/6}$.
\end{description}
</lalLaTeX>
#endif

#if 0
<lalLaTeX>
\vfill{\footnotesize\input{FindChirpSPHV}}
</lalLaTeX> 
#endif

#if 0
<lalLaTeX>
\newpage\input{FindChirpSPDataC}
</lalLaTeX>
#endif

void
LALFindChirpSPDataInit (
    LALStatus                  *status,
    FindChirpSPDataParams     **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpSPData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpSPDataParams      *params
    );

void
LALFindChirpSPDataFinalize (
    LALStatus                  *status,
    FindChirpSPDataParams     **output
    );

void
LALFindChirpBCVData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpSPDataParams      *params
    );


#if 0
<lalLaTeX>
\newpage\input{FindChirpSPTemplateC}
</lalLaTeX>
#endif

void
LALFindChirpSPTemplateInit (
    LALStatus                  *status,
    FindChirpSPTmpltParams    **output,
    FindChirpInitParams        *params
    );

void
LALFindChirpSPTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpSPTmpltParams     *params
    );

void
LALFindChirpSPTemplateFinalize (
    LALStatus                  *status,
    FindChirpSPTmpltParams    **output
    );

void
LALFindChirpBCVTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpSPTmpltParams     *params
    );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPSPH_H */
