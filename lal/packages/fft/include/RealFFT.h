#if 0 /* autodoc block */

<lalVerbatim file="RealFFTHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{RealFFT.h}}
\label{s:RealFFT.h}

Generates random numbers.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/RealFFT.h>
\end{verbatim}

\noindent Perform real-to-complex and complex-to-real fast Fourier transforms of
vectors, and sequences of vectors using the package FFTW~\cite{fj:1998}.

</lalLaTeX>

#endif /* autodoc block */


#ifndef _REALFFT_H
#define _REALFFT_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (REALFFTH, "$Id$");

#if 0 /* autodoc block */

<lalLaTeX>
\subsection*{Error conditions}
\input{RealFFTHErrTab}
</lalLaTeX>

<lalErrTable file="RealFFTHErrTab">

#endif /* autodoc block */

#define REALFFTH_ENULL 1
#define REALFFTH_ENNUL 2
#define REALFFTH_ESIZE 4
#define REALFFTH_ESZMM 8
#define REALFFTH_ESLEN 16
#define REALFFTH_ESAME 32
#define REALFFTH_ESIGN 64
#define REALFFTH_EDATA 128

#define REALFFTH_MSGENULL "Null pointer"
#define REALFFTH_MSGENNUL "Non-null pointer"
#define REALFFTH_MSGESIZE "Invalid input size"
#define REALFFTH_MSGESZMM "Size mismatch"
#define REALFFTH_MSGESLEN "Invalid/mismatched sequence lengths"
#define REALFFTH_MSGESAME "Input/Output data vectors are the same"
#define REALFFTH_MSGESIGN "Incorrect plan sign"
#define REALFFTH_MSGEDATA "Bad input data: DC/Nyquist should be real"

#if 0 /* autodoc block */

</lalErrTable>

<lalLaTeX>

\subsection*{Structures}

\begin{verbatim}
typedef struct tagRealFFTPlan RealFFTPlan;
\end{verbatim}

This structure contains the parameters necessary for performing an FFT of a
given size and direction.  The contents should not be manually adjusted.

</lalLaTeX>

#endif /* autodoc block */

typedef struct
tagRealFFTPlan
{
  INT4   sign;
  UINT4  size;
  void  *plan;
}
RealFFTPlan;

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{RealFFTC}
</lalLaTeX>

#endif /* autodoc block */

void
LALEstimateFwdRealFFTPlan (
    LALStatus       *stat,
    RealFFTPlan **plan,
    UINT4         size
    );

void
LALEstimateInvRealFFTPlan (
    LALStatus       *stat,
    RealFFTPlan **plan,
    UINT4         size
    );

void
LALMeasureFwdRealFFTPlan (
    LALStatus       *stat,
    RealFFTPlan **plan,
    UINT4         size
    );

void
LALMeasureInvRealFFTPlan (
    LALStatus       *stat,
    RealFFTPlan **plan,
    UINT4         size
    );

void
LALDestroyRealFFTPlan (
    LALStatus       *stat,
    RealFFTPlan **plan
    );


void
LALREAL4VectorFFT (
    LALStatus      *stat,
    REAL4Vector *vout,
    REAL4Vector *vinp,
    RealFFTPlan *plan
    );

void
LALREAL4VectorSequenceFFT (
    LALStatus              *stat,
    REAL4VectorSequence *vout,
    REAL4VectorSequence *vinp,
    RealFFTPlan         *plan
    );


void
LALFwdRealFFT (
    LALStatus         *stat,
    COMPLEX8Vector *vout,
    REAL4Vector    *vinp,
    RealFFTPlan    *plan
    );

void
LALInvRealFFT (
    LALStatus         *stat,
    REAL4Vector    *vout,
    COMPLEX8Vector *vinp,
    RealFFTPlan    *plan
    );

void
LALRealPowerSpectrum (
    LALStatus      *stat,
    REAL4Vector *vout,
    REAL4Vector *vinp,
    RealFFTPlan *plan
    );


void
LALFwdRealSequenceFFT (
    LALStatus                 *stat,
    COMPLEX8VectorSequence *vout,
    REAL4VectorSequence    *vinp,
    RealFFTPlan            *plan
    );

void
LALInvRealSequenceFFT (
    LALStatus                 *stat,
    REAL4VectorSequence    *vout,
    COMPLEX8VectorSequence *vinp,
    RealFFTPlan            *plan
    );

void
LALRealSequencePowerSpectrum (
    LALStatus              *stat,
    REAL4VectorSequence *vout,
    REAL4VectorSequence *vinp,
    RealFFTPlan         *plan
    );

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{RealFFTTestC}
\newpage\input{RealPowerSpectrumTestC}
</lalLaTeX>

#endif /* autodoc block */

#ifdef  __cplusplus
}
#endif

#endif /* _REALFFT_H */
