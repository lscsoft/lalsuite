#if 0 /* autodoc block */

<lalVerbatim file="SpecBufferHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{SpecBuffer.h}}
\label{s:SpecBuffer.h}

Root finding routines.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/SpecBuffer.h>
\end{verbatim}

\noindent
Generates and stores a buffer of spectra and returns an average spectrum.

</lalLaTeX>

#endif /* autodoc block */


#ifndef _SPECBUFFER_H
#define _SPECBUFFER_H

#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (SPECBUFFERH, "$Id$");

#if 0 /* autodoc block */

<lalLaTeX>
\subsection*{Error conditions}
\input{SpecBufferHErrTab}
</lalLaTeX>

<lalErrTable file="SpecBufferHErrTab">

#endif /* autodoc block */

#define SPECBUFFERH_ENULL 1
#define SPECBUFFERH_ENNUL 2
#define SPECBUFFERH_ESIZE 4
#define SPECBUFFERH_ESZMM 8
#define SPECBUFFERH_EZERO 16
#define SPECBUFFERH_ENONE 32

#define SPECBUFFERH_MSGENULL "Null pointer"
#define SPECBUFFERH_MSGENNUL "Non-null pointer"
#define SPECBUFFERH_MSGESIZE "Invalid input size"
#define SPECBUFFERH_MSGESZMM "Size mismatch"
#define SPECBUFFERH_MSGEZERO "Zero divide"
#define SPECBUFFERH_MSGENONE "No stored spectra"

#if 0 /* autodoc block */

</lalErrTable>

<lalLaTeX>

\subsection*{Structures}

\begin{verbatim}
typedef struct
tagSpectrumBufferPar
{
  INT4         numSpec;
  INT4         numPoints;
  WindowType   windowType;
  RealFFTPlan *plan;
}
SpectrumBufferPar;
\end{verbatim}

The parameter structure for creating the spectrum buffer.  The fields are:
\begin{description}
\item[\texttt{numSpec}] The number of spectra to hold in the buffer.
\item[\texttt{numPoints}] The number of points of data used to generate each
    spectrum in the buffer (note: this is not the number of points in the
    spectrum).
\item[\texttt{windowType}] The type of window used in generating the spectrum
    (see the package \texttt{window}).
\item[\texttt{plan}] The FFT plan to use in generating the spectrum
    (see the package \texttt{fft}).
\end{description}

\begin{verbatim}
typedef struct tagSpectrumBuffer SpectrumBuffer;
\end{verbatim}

The structure that holds the spectrum buffer.  The user should not modify
the fields of this structure.

\begin{verbatim}
typedef struct
tagComputeSpectrumPar
{
  WindowType   windowType;
  REAL4Vector *window;
  REAL4        wss;
  RealFFTPlan *plan;
}
ComputeSpectrumPar;
\end{verbatim}

The parameter structure for computing a windowed power spectrum from a segment
of data.  (This is for use in a lower-level function.)  The fields are:
\begin{description}
\item[\texttt{windowType}] The type of window to be used (see the package
    \texttt{window}).
\item[\texttt{window}] The actual window data.
\item[\texttt{wss}] The sum-of-squares of the window data (needed for
    normalization).
\item[\texttt{plan}] The FFT plan to use in constructing the power spectrum
    (see the package \texttt{fft}).
\end{description}

</lalLaTeX>

#endif /* autodoc block */

typedef struct
tagComputeSpectrumPar
{
  WindowType   windowType;
  REAL4Vector *window;
  REAL4        wss;
  RealFFTPlan *plan;
}
ComputeSpectrumPar;

typedef struct
tagSpectrumBuffer
{
  INT4                  numSpec;
  INT4                  specFilled;
  REAL4FrequencySeries *specBuffer;
  ComputeSpectrumPar   *specParams;
}
SpectrumBuffer;

typedef struct
tagSpectrumBufferPar
{
  INT4         numSpec;
  INT4         numPoints;
  WindowType   windowType;
  RealFFTPlan *plan;
}
SpectrumBufferPar;

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{SpecBufferC}
</lalLaTeX>

#endif /* autodoc block */

void
LALComputeSpectrum (
    LALStatus            *status,
    REAL4FrequencySeries *spectrum,
    INT2TimeSeries       *timeSeries,
    ComputeSpectrumPar   *parameters
    );

void
LALCreateSpectrumBuffer (
    LALStatus          *status,
    SpectrumBuffer    **buffer,
    SpectrumBufferPar  *params
    );

void
LALDestroySpectrumBuffer (
    LALStatus       *status,
    SpectrumBuffer **buffer
    );

void
LALAddSpectrum (
    LALStatus      *status,
    SpectrumBuffer *specBuffer,
    INT2TimeSeries *timeSeries
    );

void
LALAverageSpectrum (
    LALStatus            *status,
    REAL4FrequencySeries *spectrum,
    SpectrumBuffer       *buffer
    );

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{SpecBufferTestC}
</lalLaTeX>

#endif /* autodoc block */

#ifdef  __cplusplus
}
#endif

#endif /* _SPECBUFFER_H */
