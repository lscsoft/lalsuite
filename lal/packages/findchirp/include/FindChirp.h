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

\noindent Provides core protypes, structures and functions to filter
interferometer data for binary inspiral chirps.  

\subsection*{Binary Neutron Stars}

The important definitions
are as follows:

\begin{enumerate}
\item The equation for a stationary phase 2-pN inspiral chirp
\begin{equation}
\tilde{h_c}(f_k) = \left(\frac{5\mu}{96M_\odot}\right)^{\frac{1}{2}}
                   \left(\frac{M}{\pi^2M_\odot}\right)^{\frac{1}{3}}
                   f_k^{-\frac{7}{6}} T_\odot^{-\frac{1}{6}}
                   \exp\left(i\Psi(f_k;M,\eta)\right)
\end{equation}
where
\begin{equation}
\Psi(f_k;M,\eta) = C_0 x_k \left(C_{2.0} + x_k \left(C_{1.5} 
             + x_k \left(C_{1.0} + x_k^2 \right)\right)\right),
\end{equation}
\begin{equation}
x_k = \left(\pi M f_k\right)^\frac{1}{3},
\end{equation}
and the post-Newtonian coefficents are
\begin{eqnarray}
C_{0} &=& \frac{3}{128\eta} \\
C_{1.0} &=& \frac{3\,715}{756} + \frac{55\eta}{9} \\
C_{1.5} &=& -16\pi \\
C_{2.0} &=& \frac{15\,293\,365}{508\,032} + \frac{27\,145\eta}{504} +
            \frac{3\,085\eta^2}{72}
\end{eqnarray}

\item Signal-to-noise ratio $\rho$. We actually compute $\rho^2$:
given by
\begin{equation}
\rho^2(t_j) = \frac{16}{\sigma^2} \left( \frac{\Delta t}{N} \right)^2
              \left(\frac{2T_\odot c}{1\mathrm{Mpc}}\right)^2 d^2 A^2(M,\eta)
  \left|\sum_{k=0}^{N/2} e^{2\pi ijk/N}
       \frac{d R\tilde{v}_k k^{-\frac{7}{6}} e^{-i\Psi(f_k;M,\eta)}}
         {|d R|^2 S_v(|f_k|)}
  \right|^2
\end{equation}

\item The matched filter normalization $\sigma^2$:
\begin{equation}
\sigma^2 = 4 \left(\frac{\Delta t}{N}\right)
           \left(\frac{2T_\odot c}{1\mathrm{Mpc}}\right)^2 d^2 A^2(M,\eta)
           \sum_{k=0}^{N/2} \frac{k^{-\frac{7}{3}}}{|dR|^2S_v(|f_k|)}
\end{equation}

\item The effective distance to a source $D_\mathrm{eff}$:
\begin{equation}
D_\mathrm{eff} = \frac{\sigma^2}{\rho^2}
\end{equation}

\item The template dependent normalization
\begin{equation}
\mathtt{tmpltNorm} = \left(\frac{2T_\odot c}{1\mathrm{Mpc}}\right)^2
                     d^2 A^2(M,\eta)
\end{equation}

\item The segment dependent normalization
\begin{equation}
\mathtt{segNorm} = \sum_{k=0}^{N/2} \frac{k^{-\frac{7}{3}}}{|dR|^2 S_v(|f_k|)}
\end{equation}

\item The un-normalized matched filter output
\begin{equation}
q_j = \sum_{k=0}^{N/2} e^{2\pi ijk/N} 
      \frac{ dR \tilde{v}_k k^{-\frac{7}{6}} e^{-i\Psi(f_k;M,\eta)}}
           {|dR|^2 S_v(|f_k|)}
\end{equation}
\end{enumerate}

Then the quantities that we compute in the code are just:
\begin{equation}
\rho^2(t_j) = \frac{16}{\sigma^2}\left(\frac{\Delta t}{N}\right)^2
\cdot\mathtt{tmpltNorm}\cdot|q_j|^2,
\end{equation}
\begin{equation}
\sigma^2 = 4\left(\frac{\Delta t}{N}\right) \cdot
           \mathtt{tmpltNorm}\cdot\mathtt{segNorm}
\end{equation}
and
\begin{equation}
D_\mathrm{eff}^2 = \mathtt{tmpltNorm}\cdot \mathtt{segNorm}^2\cdot |q_j|^{-2}
\end{equation}

\subsubsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirp.h>
\end{verbatim}

\input{FindChirpHDoc}

\input{FindChirpBCVHDoc}
</lalLaTeX>
#endif

#ifndef _FINDCHIRPH_H
#define _FINDCHIRPH_H

#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/DataBuffer.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LALInspiral.h>
#include <lal/LALInspiralBank.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/FindChirpDatatypes.h>
#include <lal/FindChirpChisq.h>
#include <lal/FindChirpFilterOutputVeto.h>

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
/* </lalErrTable> */


#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif

/* --- parameter structure for all init funtions ------------------------- */
/* <lalVerbatim file="FindChirpHFindChirpInitParams"> */
typedef struct
tagFindChirpInitParams
{
  UINT4                         numSegments;
  UINT4                         numPoints;
  UINT4                         numChisqBins;
  UINT4                         ovrlap;
  BOOLEAN                       createRhosqVec;
  BOOLEAN                       createCVec;
  Approximant                   approximant;
}
FindChirpInitParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpInitParams}}
\idx[Type]{FindChirpInitParams}

\input{FindChirpHFindChirpInitParams}

\noindent This structure provides the essential information for the
filter initialisation and memory allocation functions used in the
\texttt{FindChirp} package.

\begin{description}
\item[\texttt{UINT4 numSegments}] The number of data segments to allocate
storage for.

\item[\texttt{UINT4 numPoints}] The number of discrete data points in each
data segment. 

\item[\texttt{UINT4 numChisqBins}] The number of bins used to contruct the
$\chi^2$ veto.

\item[\texttt{BOOLEAN createRhosqVec}] Debugging flag that controls whether
or not the function \texttt{FindChirpFilterSegment()} should store the output
of the filter, $\rho^2(t)$, as well as the events. Memory is only allocated
for this vector if the flag is set to 1.
\end{description}
</lalLaTeX>
#endif


/* --- the parameter structure for the data conditioning function -------- */
/* <lalVerbatim file="FindChirpHFindChirpDataParams"> */
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
\subsubsection*{Structure \texttt{FindChirpDataParams}}
\idx[Type]{FindChirpDataParams}

\input{FindChirpHFindChirpDataParams}

\noindent This structure contains the parameters needed to call the
\texttt{FindChirpSPData()} function. It should be initialized by
\texttt{FindChirpDataInit()} and destroyed by
\texttt{FindChirpDataFinalize()}. The fields are:

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
/* <lalVerbatim file="FindChirpHFindChirpTmpltParams"> */
typedef struct
tagFindChirpTmpltParams
{
  REAL8                         deltaT;
  REAL4                         fLow;
  REAL4                         dynRange;
  REAL4Vector                  *xfacVec;
  Approximant                   approximant;
  RealFFTPlan                  *fwdPlan;
}
FindChirpTmpltParams;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpTmpltParams}}
\idx[Type]{FindChirpTmpltParams}

\input{FindChirpHFindChirpTmpltParams}

\noindent This structure contains the parameters for generation of stationary
phase templates by the function \texttt{FindChirpSPTmplt()},
\texttt{FindChirpBCVTmplt()} or \texttt{FindChirpBCVSpinTmplt()}
It should be initialized by \texttt{FindChirpTmpltInit()} and destroyed by
\texttt{FindChirpTmpltFinalize()}. The fields are:

\begin{description}
\item[\texttt{REAL8 deltaT}] The sampling interval $\Delta t$.

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



/* --- parameter structure for the filtering function -------------------- */
/* <lalVerbatim file="FindChirpHFindChirpFilterParams"> */
typedef struct
tagFindChirpFilterParams
{
  REAL8                         deltaT;
  REAL4                         rhosqThresh;
  REAL4                         chisqThresh;
  REAL4                         norm;
  UINT4                         maximiseOverChirp;
  Approximant                   approximant;
  COMPLEX8Vector               *qVec;
  COMPLEX8Vector               *qVecBCV;
  COMPLEX8Vector               *qVecBCVSpin1;
  COMPLEX8Vector               *qVecBCVSpin2;
  COMPLEX8Vector               *qtildeVec;
  COMPLEX8Vector               *qtildeVecBCV;
  COMPLEX8Vector               *qtildeVecBCVSpin1;
  COMPLEX8Vector               *qtildeVecBCVSpin2;
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
\subsubsection*{Structure \texttt{FindChirpFilterParams}}
\idx[Type]{FindChirpFilterParams}

\input{FindChirpHFindChirpFilterParams}

\noindent This structure provides the parameters used by the
\texttt{FindChirpFilterSegment()} function.

\begin{description}
\item[\texttt{REAL8 deltaT}] The timestep for the sampled data. Must be
set on entry.

\item[\texttt{REAL4 rhosqThresh}] The value to threshold signal to noise
ratio square, $\rho^2$, on. If the signal to noise exceeds this value, then a
candidate event is generated. A $\chi^2$ veto is performed, if requested,
otherwise an \texttt{InspiralEvent} is generated. Must be $\ge 0$ on entry.

\item[\texttt{REAL4 chisqThresh}] The value to threshold the $\chi^2$ veto on.
If the $chi^2$ veto is below this threshold the candidate event an
\texttt{InspiralEvent} is generated. Must be $\ge 0$ on entry.

\item[\texttt{REAL4 norm}] On exit this contains the normalisation constant
that relates the quantity $q_j$ with the signal to noise squared, 
$\rho^2(t_j)$ by 
\begin{equation}
\rho^2(t_j) = \textrm{norm} \times \left|q_j\right|^2.
\end{equation}

\item[\texttt{UINT4 maximiseOverChirp}] If not zero, use algorithm that
maximised over chirp lengths. Otherwise record all points that pass
the $\rho^2$ threshold as events.

\item[\texttt{COMPLEX8Vector *qVec}] Pointer to vector allocated by 
\texttt{FindChirpFilterInit()} to store the quantity $q_j$. Set to the
value of $q_j$ on exit. Must not be NULL.

\item[\texttt{COMPLEX8Vector *qVecBCV}] Pointer to the additional vector
required for the BCV templates, allocated by
\texttt{FindChirpFilterInit()} to store the quantity $q_j$. Set to the
value of $q_j$ on exit. Must not be NULL.

\item[\texttt{COMPLEX8Vector *qtildeVec}] Pointer to vector allocated by 
\texttt{FindChirpFilterInit()} to store the quantity $\tilde{q}_k$. Set to the
value of $\tilde{q}_k$ on exit. Must not be NULL

\item[\texttt{COMPLEX8Vector *qtildeVecBCV}] Pointer to the additional
vector required for the BCV templates, allocated by 
\texttt{FindChirpFilterInit()} to store the quantity $\tilde{q}_k$. Set to the
value of $\tilde{q}_k$ on exit. Must not be NULL

\item[\texttt{ComplexFFTPlan *invPlan}] Pointer to FFTW plan created by 
\texttt{FindChirpFilterInit()} to transform the quantity $\tilde{q}_k$ to
${q}_j$ usimg the inverse DFT. Must not be NULL.

\item[\texttt{REAL4Vector *rhosqVec}] Pointer to a vector that is set to
$\rho^2(t_j)$ on exit. If NULL $\rho^2(t_j)$ is not stored.

\item[\texttt{REAL4Vector *chisqVec}] Workspace vector used to compute and
store $\chi^2(t_j)$. Must not be NULL if \texttt{numChisqBins} is greater than
zero. Contains $\chi^2(t_j)$ on exit.

\item[\texttt{FindChirpChisqParams *chisqParams}] Pointer to parameter
structure for \texttt{FindChirpChisqVeto()} function. Must not be NULL if
\texttt{numChisqBins} is greater than zero.

\item[\texttt{FindChirpChisqInput *chisqInput}] Pointer to input data
structure for the \texttt{FindChirpChisqVeto()} function
and the \texttt{FindCHirpBCVChisqVeto()}. Must not be NULL if
\texttt{numChisqBins} is greater than zero.

\item[\texttt{FindChirpChisqInput *chisqInputBCV}] Pointer to input data
structure for the \texttt{FindChirpBCVChisqVeto()} function. Must not be NULL 
if \texttt{numChisqBins} is greater than zero.

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
\subsubsection*{Filter function input structures}
</lalLaTeX>
#endif


/* --- input to the filtering functions --------------------------------- */
/* <lalVerbatim file="FindChirpHFindChirpFilterInput"> */
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
\subsubsection*{Structure \texttt{FindChirpFilterInput}}
\idx[Type]{FindChirpSegmentVector}

\input{FindChirpHFindChirpFilterInput}

\noindent This structure groups the input data required for the 
\texttt{FindChirpFilterSegment()} function into a single structure.

\begin{description}
\item[\texttt{InspiralTemplate *tmplt}] Pointer the structure that contains
the parameters of the template chirp.

\item[\texttt{FindChirpTemplate *fcTmplt}] Pointer to the input template
in a form that can be used by \texttt{FindChirpFilterSegment()}

\item[\texttt{FindChirpSegment *segment}] Pointer to the input data segment
in a form that can be used by \texttt{FindChirpFilterSegment()}
\end{description}
</lalLaTeX>
#endif

#if 0
<lalLaTeX>
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

#if 0
<lalLaTeX>
\newpage\input{FindChirpLinkedListC}
</lalLaTeX>
#endif

void
LALFindChirpCreateTmpltNode (
    LALStatus                  *status,
    InspiralTemplate           *tmplt,
    InspiralTemplateNode      **tmpltNode
    );

void
LALFindChirpDestroyTmpltNode ( 
    LALStatus                  *status,
    InspiralTemplateNode      **tmpltNode
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
LALFindChirpFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    );

void
LALFindChirpBCVFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
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

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPH_H */
