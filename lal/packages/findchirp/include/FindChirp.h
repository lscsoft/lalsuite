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
#define FINDCHIRPH_MSGERHOT "Rhosq threshold is zero or negative"
#define FINDCHIRPH_MSGECHIT "Chisq threshold is zero or negative"
#define FINDCHIRPH_MSGECRUP "Chirp length or invSpecTrunc too long for length of data segment"
#define FINDCHIRPH_MSGESMSM "Size mismatch between vectors"
#define FINDCHIRPH_MSGEHETR "Attempting to simulate heterodyned GW"
#define FINDCHIRPH_MSGEDFDT "Waveform sampling interval is too large"
#define FINDCHIRPH_MSGEAPRX "Incorrect waveform approximant"
/* </lalErrTable> */


#if 0
<lalLaTeX>
\subsection*{Types}
</lalLaTeX>
#endif


/*
 *
 * typedefs of structures used by the findchirp functions
 *
 */


#if 0
<lalLaTeX>
\subsubsection*{Input and output data structures}
</lalLaTeX>
#endif

/* --- structure for storing the parameters needed to do an injection ---- */
typedef struct
tagFindChirpStandardCandle
{
  CHAR                          ifo[2];
  InspiralTemplate              tmplt;
  REAL4                         rhosq;
  REAL4                         sigmasq;
  REAL4                         effDistance;
}
FindChirpStandardCandle;

/* --- vector of DataSegment, as defined the framedata package ----------- */
/* <lalVerbatim file="FindChirpHDataSegmentVector"> */
typedef struct
tagDataSegmentVector
{
  UINT4                         length;
  DataSegment                  *data;
}
DataSegmentVector;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{DataSegmentVector}}
\idx[Type]{DataSegmentVector}

\input{FindChirpHDataSegmentVector}

\noindent This structure provides a LAL like vector structure for the
\texttt{DataSegment} structure defined in the package \texttt{framedata}

\begin{description}
\item[\texttt{UINT4 length}] Number of \texttt{DataSegment} structres in the
vector

\item[\texttt{DataSegment *data}] Pointer to the data.
\end{description}
</lalLaTeX>
#endif

/* --- processed data segment used by FindChirp filter routine ----------- */
/* <lalVerbatim file="FindChirpHFindChirpSegment"> */
typedef struct
tagFindChirpSegment
{
  COMPLEX8FrequencySeries      *data;
  COMPLEX8FrequencySeries      *dataBCV;
  UINT4Vector                  *chisqBinVec;
  UINT4Vector 		       *chisqBinVecBCV;
  REAL8                         deltaT;
  REAL4                         segNorm;
  REAL4                         a1;     
  REAL4                         b1;
  REAL4                         b2;     
  REAL4                         fLow;
  UINT4                         invSpecTrunc;
  UINT4                         number;
  INT4                          level;
  Approximant                   approximant;
}
FindChirpSegment;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpSegment}}
\idx[Type]{FindChirpSegment}

\input{FindChirpHFindChirpSegment}

\noindent This structure contains the conditioned input data and its
parameters for the \texttt{FindChirpFilter()} function.

\begin{description}
\item[\texttt{COMPLEX8FrequencySeries *data}] The conditioned input data, 
used for the stationary phase chirps and the BCV templates.
The conditioning performed is as described in the documentation for the
module \texttt{FindChirpSPData.c}

\item[\texttt{COMPLEX8FrequencySeries *dataBCV}] The conditioned input data,
used only for the BCV templates.
The conditioning performed is as described in the documentation for the
module \texttt{FindChirpBCVData.c}

\item[\texttt{UINT4Vector *chisqBinVec}] A vector containing the indices of
the boundaries of the bins of equal power for the $\chi^2$ veto created by 
\texttt{FindChirpSPData()} or \texttt{FindChirpBCVData()}.

\item[\texttt{UINT4Vector *chisqBinVecBCV}] A vector containing the indices of
the boundaries of the bins of equal power for the second contribution to the 
$\chi^2$ statistic for the BCV templates, created by 
\texttt{FindChirpBCVData()}

\item[\texttt{REAL8 deltaT}] The time step $\Delta$ of the time series 
input data.

\item[\texttt{REAL4 segNorm}] The template independent part of the 
normalisation constant $\sigma$.

\item[\texttt{REAL4 a1}] BCV-template normalization parameter.

\item[\texttt{REAL4 b1}] BCV-template normalization parameter.

\item[\texttt{REAL4 b2}] BCV-template normalization parameter.

\item[\texttt{UINT4 invSpecTrunc}] The number of points to which the inverse 
power spectrum \ospsd is truncated to in the time domain in order to smooth
out high $Q$ features in the power spectrum.

\item[\texttt{UINT4 number}] A unique identification number for the 
\texttt{FindChirpDataSegment}. This will generally correspond to the number in
the \texttt{DataSegment} from which the conditioned data was computed.
\end{description}
</lalLaTeX>
#endif

/* --- vector of FindChirpSegment defined above -------------------------- */
/* <lalVerbatim file="FindChirpHFindChirpSegmentVector"> */
typedef struct
tagFindChirpSegmentVector
{
  UINT4                         length;
  FindChirpSegment             *data;
}
FindChirpSegmentVector;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpSegmentVector}}
\idx[Type]{FindChirpSegmentVector}

\input{FindChirpHFindChirpSegmentVector}

\noindent This structure provides a LAL like vector structure for the
\texttt{FindChirpSegment} structure defined above.

\begin{description}
\item[\texttt{UINT4 length}] Number of \texttt{FindChirpSegment} structres in
the vector

\item[\texttt{DataSegment *data}] Pointer to the data.
\end{description}
</lalLaTeX>
#endif

/* --- structure to contain an inspiral template ------------------------- */
/* <lalVerbatim file="FindChirpHFindChirpTemplate"> */
typedef struct
tagFindChirpTemplate
{
  COMPLEX8Vector               *data;
  REAL4                         tmpltNorm;
  Approximant                   approximant;
}
FindChirpTemplate;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\subsubsection*{Structure \texttt{FindChirpTemplate}}
\idx[Type]{FindChirpTemplate}


\noindent This structure provides contains the frequency domain representation
of the cosine phase inspiral template $\tilde{h_c}(f)$.

\begin{description}
\item[\texttt{COMPLEX8Vector *data}] A vector containing $\tilde{h_c}(f)$. Note
that in the future, this will be changed to a \texttt{COMPLEX8FrequencySeries}.

\item[\texttt{REAL4 tmpltNorm}] The template dependent part of the 
normalisation constant $\sigma$.
\end{description}
</lalLaTeX>
#endif


/*
 *
 * typedefs of parameter structures used by functions in findchirp
 *
 */


#if 0
<lalLaTeX>
\subsubsection*{Initalisation and parameter structures}
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


/* --- parameter structure for the filtering function -------------------- */
/* <lalVerbatim file="FindChirpHFindChirpFilterParams"> */
typedef struct
tagFindChirpFilterParams
{
  REAL4                         deltaT;
  REAL4                         rhosqThresh;
  REAL4                         chisqThresh;
  REAL4                         norm;
  UINT4                         maximiseOverChirp;
  BOOLEAN                       computeNegFreq;
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
  REAL4Vector                  *chisqVec;
  FindChirpChisqParams         *chisqParams;
  FindChirpChisqInput          *chisqInput;
  FindChirpChisqInput          *chisqInputBCV;
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
\item[\texttt{REAL4 deltaT}] The timestep for the sampled data. Must be
set on entry.  {FIXME: This should be a \texttt{REAL8}}

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

\item[\texttt{BOOLEAN computeNegFreq}] Currently unused. Must be set to
$0$ on entry.

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
  InspiralTemplate             *tmplt;
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

void
LALRandomPPNParamStruc (
    LALStatus                  *status,
    PPNParamStruc              *PPNparams,
    InspiralCoarseBankIn       *massParams,
    RandomParams               *randomParams
    );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPH_H */
