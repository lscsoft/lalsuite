#if 0 /* autodoc block */

<lalVerbatim file="ComplexFFTHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{ComplexFFT.h}}
\label{s:ComplexFFT.h}

Performs complex-to-complex FFTs.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/ComplexFFT.h>
\end{verbatim}

Perform complex-to-complex fast Fourier transforms of vectors using the
package FFTW~\cite{fj:1998}.

</lalLaTeX>

#endif /* autodoc block */

#ifndef _COMPLEXFFT_H
#define _COMPLEXFFT_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( COMPLEXFFTH, "$Id$" );

#if 0 /* autodoc block */

<lalLaTeX>
\subsection*{Error conditions}
\input{ComplexFFTHErrTab}
</lalLaTeX>

<lalErrTable file="ComplexFFTHErrTab">

#endif /* autodoc block */

#define COMPLEXFFTH_ENULL 1
#define COMPLEXFFTH_ENNUL 2
#define COMPLEXFFTH_ESIZE 4
#define COMPLEXFFTH_ESZMM 8
#define COMPLEXFFTH_ESLEN 16
#define COMPLEXFFTH_ESAME 32

#define COMPLEXFFTH_MSGENULL "Null pointer"
#define COMPLEXFFTH_MSGENNUL "Non-null pointer"
#define COMPLEXFFTH_MSGESIZE "Invalid input size"
#define COMPLEXFFTH_MSGESZMM "Size mismatch"
#define COMPLEXFFTH_MSGESLEN "Invalid/mismatched sequence lengths"
#define COMPLEXFFTH_MSGESAME "Input/Output data vectors are the same"

#if 0 /* autodoc block */

</lalErrTable>

<lalLaTeX>

\subsection*{Structures}

\begin{verbatim}
typedef struct tagComplexFFTPlan ComplexFFTPlan;
\end{verbatim}

This structure contains the parameters necessary for performing an FFT of a
given size and direction.  The contents should not be manually adjusted.

</lalLaTeX>

#endif /* autodoc block */

typedef struct
tagComplexFFTPlan
{
  INT4   sign;
  UINT4  size;
  void  *plan;
}
ComplexFFTPlan;

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{ComplexFFTC}
</lalLaTeX>

#endif /* autodoc block */

void
LALEstimateFwdComplexFFTPlan (
    LALStatus          *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    );

void
LALEstimateInvComplexFFTPlan (
    LALStatus          *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    );

void
LALMeasureFwdComplexFFTPlan (
    LALStatus          *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    );

void
LALMeasureInvComplexFFTPlan (
    LALStatus          *stat,
    ComplexFFTPlan **plan,
    UINT4            size
    );

void
LALDestroyComplexFFTPlan (
    LALStatus          *stat,
    ComplexFFTPlan **plan
    );


void
LALCOMPLEX8VectorFFT (
    LALStatus         *stat,
    COMPLEX8Vector *vout,
    COMPLEX8Vector *vinp,
    ComplexFFTPlan *plan
    );

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{ComplexFFTTestC}
</lalLaTeX>

#endif /* autodoc block */

#ifdef  __cplusplus
}
#endif

#endif /* _COMPLEXFFT_H */
