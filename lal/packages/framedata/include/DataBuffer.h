#if 0 /* autodoc block */

<lalVerbatim file="DataBufferHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{DataBuffer.h}}
\label{s:DataBuffer.h}

Root finding routines.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/DataBuffer.h>
\end{verbatim}

\noindent Gets segments of IFO\_DMRO data along with its (averaged over several
segments) power spectrum and response function.

</lalLaTeX>

#endif /* autodoc block */

#ifndef _DATABUFFER_H
#define _DATABUFFER_H

#include <lal/LALDatatypes.h>
#include <lal/FrameData.h>
#include <lal/SpecBuffer.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (DATABUFFERH, "$Id$");

#if 0 /* autodoc block */

<lalLaTeX>
\subsection*{Error conditions}
\input{DataBufferHErrTab}
</lalLaTeX>

<lalErrTable file="DataBufferHErrTab">

#endif /* autodoc block */

#define DATABUFFERH_ENULL 1
#define DATABUFFERH_ENNUL 2
#define DATABUFFERH_ESIZE 4
#define DATABUFFERH_ESZMM 8

#define DATABUFFERH_MSGENULL "Null pointer"
#define DATABUFFERH_MSGENNUL "Non-null pointer"
#define DATABUFFERH_MSGESIZE "Invalid input size"
#define DATABUFFERH_MSGESZMM "Size mismatch"

#if 0 /* autodoc block */

</lalErrTable>

<lalLaTeX>

\subsection*{Structures}

\begin{verbatim}
typedef struct
tagDataSegment
{
  REAL4TimeSeries         *chan;
  REAL4FrequencySeries    *spec;
  COMPLEX8FrequencySeries *resp;
  INT4                     number;
  UINT4                    analyzeSegment;
}
DataSegment;
\end{verbatim}

%\begin{verbatim}
%typedef struct
%tagDataSegment
%{
%  INT2TimeSeries          *data;
%  REAL4FrequencySeries    *spec;
%  COMPLEX8FrequencySeries *resp;
%  INT4                     endOfData;
%  INT4                     newLock;
%  INT4                     newCal;
%  INT4                     number;
%  INT4                     level;
%}
%DataSegment;
%\end{verbatim}

The data structure returned by the data aquisition routine.  The fields are:
\begin{description}
\item[\texttt{data}] The time-series data.
\item[\texttt{spec}] The corresponding estimated (averaged) power spectrum.
\item[\texttt{resp}] The current response function, built from the swept-sine
    calibration information.
\item[\texttt{endOfData}] Boolean that is non-zero if there is no more data.
\item[\texttt{newLock}] Boolean that is non-zero if a new locked segment has
    been entered.
\item[\texttt{newCal}] Boolean that is non-zero if a new calibration
    information has been read.
\item[\texttt{number}] The segment number.
\item[\texttt{level}] A field that Duncan wanted.
\end{description}

\begin{verbatim}
typedef struct
tagDataBufferPar
{
  INT4         numSpec;
  INT4         numPoints;
  WindowType   windowType;
  RealFFTPlan *plan;
  CHAR        *framePath;
}
DataBufferPar;
\end{verbatim}

This is the parameter structure used when creating a data buffer.  The fields
are:
\begin{description}
\item[\texttt{numSpec}] The number of spectra to average to get an average
    spectrum.
\item[\texttt{numPoints}] The number of points of data in each segment.
\item[\texttt{windowType}] The type of window to use when creating the
    spectra (see the \texttt{window} package).
\item[\texttt{plan}] The FFT plan to use when creating the spectra (see the
    \texttt{fft} package).
\item[\texttt{framePath}] The path of the directory containing the frame files.
\end{description}

\begin{verbatim}
typedef struct
tagDataBuffer
{
  INT4 endOfData;
  INT4 newLock;
  INT4 newCal;
  /* ... private data ... */
}
\end{verbatim}
This is the buffer of data segments.  The ``public'' fields are:
\begin{description}
\item[\texttt{endOfData}] Boolean that is non-zero if there is no more data.
\item[\texttt{newLock}] Boolean that is non-zero if a new locked segment has
    been entered.
\item[\texttt{newCal}] Boolean that is non-zero if a new calibration
    information has been read.
\end{description}

\begin{verbatim}
typedef struct tagDataBlock DataBlock;
\end{verbatim}

An internal structure.

</lalLaTeX>

#endif /* autodoc block */

typedef struct
tagDataBlock
{
  INT4 number;
  INT4 continuous;
  INT4 anomalous;
  INT2TimeSeries *framedata;
}
DataBlock;

typedef struct
tagDataBuffer
{
  /*
   * public data
   */
  INT4 endOfData;
  INT4 newLock;
  INT4 newCal;
  /*
   * private data
   */
  FrameData          *frameData;
  SpectrumBuffer     *specBuffer;
  INT2VectorSequence *dataBuffer;
  DataBlock          *blockArray;
  INT4                blockArraySize;
  INT4                blocksFilled;
  INT4                nextData;
  INT4                nextBlock;
  INT4                deltaPoints;
  INT4                numSent;
  INT4                first;
}
DataBuffer;

typedef struct
tagDataBufferPar
{
  INT4         numSpec;
  INT4         numPoints;
  WindowType   windowType;
  RealFFTPlan *plan;
  CHAR        *framePath;
}
DataBufferPar;

typedef struct
tagDataSegment
{
  REAL4TimeSeries         *chan;
  REAL4FrequencySeries    *spec;
  COMPLEX8FrequencySeries *resp;
  INT4                     number;
  UINT4                    analyzeSegment;
}
DataSegment;

#if 0
/* this is the old data segment structure */
typedef struct
tagDataSegment
{
  INT2TimeSeries          *data;
  REAL4TimeSeries         *real4Data;
  REAL4FrequencySeries    *spec;
  COMPLEX8FrequencySeries *resp;
  INT4                     endOfData;
  INT4                     newLock;
  INT4                     newCal;
  INT4                     number;
  INT4                     level;
}
DataSegment;
#endif

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{DataBufferC}
</lalLaTeX>

#endif /* autodoc block */

void
LALCreateDataBuffer (
    LALStatus         *status,
    DataBuffer    **buffer,
    DataBufferPar  *params
    );

void
LALDestroyDataBuffer (
    LALStatus      *status,
    DataBuffer **buffer
    );

void
LALGetData (
    LALStatus      *status,
    DataSegment *output,
    INT4         advance,
    DataBuffer  *buffer
    );

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{DataBufferTestC}
</lalLaTeX>

#endif /* autodoc block */

#ifdef  __cplusplus
}
#endif

#endif /* _DATABUFFER_H */
