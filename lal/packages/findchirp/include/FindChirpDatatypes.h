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

\noindent Provides core protypes for the core datatypes using in FindChirp.

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
</lalLaTeX>
#endif

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

/* --- structure for managing a list of inspiral templates --------------- */
/* <lalVerbatim file="FindChirpHInspiralTemplateNode"> */
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
\subsubsection*{Structure \texttt{InspiralTemplateNode}}
\idx[Type]{InspiralTemplateNode}

\input{FindChirpHInspiralTemplateNode}

\noindent This structure provides a method of constucting doubly linked
lists of \texttt{InspiralTemplate} structures. The fields are:

\begin{description}
\item[\texttt{struct tagInspiralTemplateNode *next}] The next structure in
the linked list.

\item[\texttt{struct tagInspiralTemplateNode *prev}] The previous structure in
the linked list.

\item[\texttt{InspiralTemplate *tmpltPtr}] A pointer to an \texttt{InspiralTemplate} structure.
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
