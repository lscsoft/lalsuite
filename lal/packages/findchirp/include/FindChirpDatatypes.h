/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpDatatypes.h
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpDatatypesHV">
Author: Brown, D. A.
$Id$
</lalVerbatim> 

<lalLaTeX>
\section{Header \texttt{FindChirpDatatypes.h}}
\label{s:FindChirpDatatypes.h}

\noindent Provides core protypes for the core datatypes using in
findchirp.

\subsubsection*{Synopsis}

\begin{verbatim}
#include <lal/FindChirpDatatypes.h>
\end{verbatim}

</lalLaTeX>
#endif

#ifndef _FINDCHIRPDATATYPESH_H
#define _FINDCHIRPDATATYPESH_H

#include <lal/LALDatatypes.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif


NRCSID (FINDCHIRPDATATYPESH, "$Id$");

#if 0
<lalLaTeX> 
\subsection*{Error codes} 

\noindent None.

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
\subsubsection*{Structure \texttt{FindChirpStandardCandle}}
\idx[Type]{FindChirpStandardCandle}

\noindent Struture used to contain the standard candle data for an inspiral 
signal. The standard candle is typically the distance to which the 
interferometer is sensitive to an optimally oriented inspiral with the masses
stored in the \texttt{tmplt} structure at a signal-to-noise ratio given by
\texttt{rhosq}.

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpStandardCandle
{
  CHAR                          ifo[3];
  InspiralTemplate              tmplt;
  REAL4                         rhosq;
  REAL4                         sigmasq;
  REAL4                         effDistance;
}
FindChirpStandardCandle;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\begin{description}
\item[\texttt{CHAR ifo[3]}] NULL terminated ifo name.

\item[\texttt{InspiralTemplate tmplt}] Inspiral template used to compute the
standard candle.

\item[\texttt{REAL4 rhosq}] The desired signal-to-noise ratio squared $\rho^2$
of the candle in the data.

\item[\texttt{REAL4 sigmasq}] The variance of the matched filter $\sigma^2$ 
for the data used to calculate the standard candle.

\item[\texttt{REAL4 effDistance}] The distance at which an optimally oriented
inspiral with the masses given by \texttt{tmplt} would give the
signal-to-noise ratio squared \texttt{rhosq}.  
\end{description} 

\subsubsection*{Structure \texttt{DataSegmentVector}}
\idx[Type]{DataSegmentVector}

\noindent Struture used to contain an array of \texttt{DataSegments}.
\texttt{DataSegments} are defined in the header \texttt{DataBuffer.h} of the
framedata package and contains an \texttt{INT4 number} used to identify the data
segment as well as pointers to a data channel (\texttt{REAL4TimeSeries
*chan}), a power spectral estimate (\texttt{REAL4FrequencySeries *spec})
and a response function  (\texttt{COMPLEX8FrequencySeries *resp}). A vector
of \texttt{DataSegmentVector} therefore contains all the input data from 
the interferometer that findchirp needs.

</lalLaTeX>
#endif
/* <lalVerbatim> */
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
\begin{description}
\item[\texttt{UINT4 length}] Number of \texttt{DataSegment} structures in the
vector.

\item[\texttt{DataSegment *data}] Pointer to an array of \texttt{DataSegment}
structures.
\end{description}

\subsubsection*{Structure \texttt{InspiralTemplateNode}}
\idx[Type]{InspiralTemplateNode}

\noindent This structure provides a method of constucting doubly linked
lists of \texttt{InspiralTemplate} structures.
</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagInspiralTemplateNode
{
  struct tagInspiralTemplateNode       *next;
  struct tagInspiralTemplateNode       *prev;
  InspiralTemplate                     *tmpltPtr;
}
InspiralTemplateNode;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\begin{description}
\item[\texttt{struct tagInspiralTemplateNode *next}] The next structure in
the linked list.

\item[\texttt{struct tagInspiralTemplateNode *prev}] The previous structure in
the linked list.

\item[\texttt{InspiralTemplate *tmpltPtr}] A pointer to an 
\texttt{InspiralTemplate} structure containing the template parameters.
\end{description}

\subsubsection*{Structure \texttt{FindChirpSegment}}
\idx[Type]{FindChirpSegment}

\noindent This structure contains the conditioned input data and its
parameters and is one of the required inputs to the \texttt{FindChirpFilter()}
function.

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpSegment
{
  COMPLEX8FrequencySeries      *data;
  COMPLEX8FrequencySeries      *dataBCV;
  UINT4Vector                  *chisqBinVec;
  UINT4Vector                  *chisqBinVecBCV;
  REAL8                         deltaT;
  REAL4Vector                  *segNorm;
  REAL4Vector                  *a1;     
  REAL4Vector                  *b1;
  REAL4Vector                  *b2;     
  REAL4Vector                  *tmpltPowerVec;
  REAL4Vector                  *tmpltPowerVecBCV;
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

\begin{description}
\item[\texttt{COMPLEX8FrequencySeries *data}] The conditioned data used as
part of the matched filter correllation. The exact content of this structure 
is determined by which data conditioning routine is called (stationary phase,
time domain, BCV or spinning BCV). The data in this structure is denoted
$\tilde{F}_k$. It typically contains 
\begin{equation}
\tilde{F}_k = \frac{d\tilde{v}_k f(k)}
{d^2|R|^2S_v\left(\left|f_k\right|\right)}
\end{equation}
where for stationary phase and BCV templates
\begin{equation}
f(k) = \left(\frac{k}{N}\right)^{-\frac{7}{6}}
\end{equation}
and for time domain templates $f(k) \equiv 1$.

\item[\texttt{COMPLEX8FrequencySeries *dataBCV}] Conditioned input data used
only for the BCV templates.  The conditioning performed is as described in the
documentation for the module \texttt{FindChirpBCVData.c}

\item[\texttt{UINT4Vector *chisqBinVec}] A vector containing the indices of
the boundaries of the bins of equal power for the $\chi^2$ veto.

\item[\texttt{UINT4Vector *chisqBinVecBCV}] A vector containing the indices of
the boundaries of the bins of equal power for the second contribution to the 
$\chi^2$ statistic for the BCV templates.

\item[\texttt{REAL8 deltaT}] The time step $\Delta t$ of the input
data channel used to create the \texttt{FindChirpSegment}.

\item[\texttt{REAL4Vector *segNorm}] The quantity segment dependent
normalization quantity $\mathcal{S(k)}$. This for stationary phase 
templates this is 
\begin{equation}
\mathcal{S}(k) = \sum_{k'=1}^{k} 
\frac{\left(\frac{k'}{N}\right)^{-\frac{7}{3}}}{d^2|R|^2S_v\left(\left|f_k'\right|\right)}
\end{equation}
where $1 \le k \le N/2$. For time domain templates, this quantity is recomputed for each template, with $|h(f)|^2$ replacing $(k'/N)^{-7/3}$.

\item[\texttt{REAL4 a1}] BCV-template normalization parameter.

\item[\texttt{REAL4 b1}] BCV-template normalization parameter.

\item[\texttt{REAL4 b2}] BCV-template normalization parameter.

\item[\texttt{REAL4Vector *tmpltPowerVec}] The weighted power in the
template. This is
\begin{equation}
\frac{\left(\frac{k}{N}\right)^{-\frac{7}{3}}}{d^2|R|^2S_v\left(\left|f_k\right|\right)}
\end{equation}
for stationary phase and BCV templates. For time domain templates, this
quantity is recomputed for each template, with $|h(f)|^2$ replacing
$(k'/N)^{-7/3}$. This quantity is used in the computation of the $\chi^2$ bin
boundaries.

\item[\texttt{REAL4Vector *tmpltPowerVecBCV}] Additional weighted template
power for BCV templates.

\item[\texttt{REAL4 flow}] The (frequency domain) low frequency cutoff for the
matched filter, $f_\mathrm{low}$.

\item[\texttt{UINT4 invSpecTrunc}] The number of points to which the inverse 
power spectrum \ospsd is truncated to in the time domain in order to truncate
the impulse response time of the matched filter.

\item[\texttt{UINT4 number}] A unique identification number for the 
\texttt{FindChirpDataSegment}. This will generally correspond to the number in
the \texttt{DataSegment} from which the conditioned data was computed.

\item[\texttt{INT4 level}] A search level, used by the heirarchical search
engine to determine if a particular data segment should be filtered against
a particular template.
\end{description}

\subsubsection*{Structure \texttt{FindChirpSegmentVector}}
\idx[Type]{FindChirpSegmentVector}

\noindent A vector of \texttt{FindChirpSegment} structures, defined above.

</lalLaTeX>
#endif
/* <lalVerbatim> */
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

\begin{description}
\item[\texttt{UINT4 length}] Number of \texttt{FindChirpSegment} structures in
the vector

\item[\texttt{DataSegment *data}] Pointer to an array of
\texttt{FindChirpSegment} structures.  
\end{description}

\subsubsection*{Structure \texttt{FindChirpTemplate}}
\idx[Type]{FindChirpTemplate}

\noindent This structure contains a frequency domain template used as input
to the \texttt{FindChirpFilter()} routine. This may either be a template
generated in the frequency domain or the Fourier transform of template
generated in the time domain.

</lalLaTeX>
#endif
/* <lalVerbatim> */
typedef struct
tagFindChirpTemplate
{
  InspiralTemplate              tmplt;
  COMPLEX8Vector               *data;
  REAL4                         tmpltNorm;
  REAL8                         momentI;
  REAL8                         momentJ;
  REAL8                         momentK;
  REAL8                         rootMomentI;
  REAL8                         numFactor;  
  REAL8                         numFactor1;
  REAL8                         numFactor2;
  REAL8                         numFactor3;  
  REAL8Vector                  *A1BCVSpin;
  REAL8Vector                  *A2BCVSpin;
  REAL8Vector                  *A3BCVSpin;
}
FindChirpTemplate;
/* </lalVerbatim> */
#if 0
<lalLaTeX>
\begin{description}
\item[\texttt{InspiralTemplate tmplt}] contains the parameters of the 
\texttt{FindChirpTemplate}, such as the binary masses and the waveform
approximant.

\item[\texttt{COMPLEX8Vector *data}] Frequency domain data containing the
\emph{findchirp template} $\tilde{F}_k$. For a stationary phase template
\begin{equation}
\tilde{T}_k = \exp\left[i\Psi(f_k;M,\eta)\right] \Theta\left(k-k_\mathrm{isco}\right).
\end{equation}
For a time domain template $\tilde{T}_k$ contains the Fourier transform of 
a waveform generated by the inspiral package.

\item[\texttt{REAL4 tmpltNorm}] The template dependent normalisation constant
$\mathcal{T}$.

\item[\texttt{REAL8 momentI}] Undocumented BCV normalization constant.

\item[\texttt{REAL8 momentJ}] Undocumented BCV normalization constant.

\item[\texttt{REAL8 momentK}] Undocumented BCV normalization constant.

\item[\texttt{REAL8 rootMomentI}] Undocumented BCV normalization constant.

\item[\texttt{REAL8 numFactor}] Undocumented BCV normalization constant.

\item[\texttt{REAL8 numFactor1}] Undocumented BCV normalization constant.

\item[\texttt{REAL8 numFactor2}] Undocumented BCV normalization constant.

\item[\texttt{REAL8 numFactor3}] Undocumented BCV normalization constant.

\item[\texttt{REAL8 A1BCVSpin}] Undocumented spinning BCV template data.

\item[\texttt{REAL8 A2BCVSpin}] Undocumented spinning BCV template data.

\item[\texttt{REAL8 A3BCVSpin}] Undocumented spinning BCV template data.
\end{description}
</lalLaTeX>
#endif

#if 0
<lalLaTeX>
\vfill{\footnotesize\input{FindChirpDatatypesHV}}
</lalLaTeX> 
#endif

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _FINDCHIRPDATATYPESH_H */
