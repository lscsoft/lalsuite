/**** <lalVerbatim file="ComplexFFTHV">
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * 
 * \section{Header \texttt{ComplexFFT.h}}
 * \label{s:ComplexFFT.h}
 * 
 * Performs complex-to-complex FFTs.
 * 
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/ComplexFFT.h>
 * \end{verbatim}
 * 
 * Perform complex-to-complex fast Fourier transforms of vectors using the
 * package FFTW~\cite{fj:1998}.
 * 
 **** </lalLaTeX> */

#ifndef _COMPLEXFFT_H
#define _COMPLEXFFT_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( COMPLEXFFTH, "$Id$" );

/**** <lalLaTeX>
 * \subsection*{Error conditions}
 **** </lalLaTeX> */
/**** <lalErrTable> */

#define COMPLEXFFTH_ENULL 1
#define COMPLEXFFTH_ENNUL 2
#define COMPLEXFFTH_ESIZE 4
#define COMPLEXFFTH_ESZMM 8
#define COMPLEXFFTH_ESLEN 16
#define COMPLEXFFTH_ESAME 32
#define COMPLEXFFTH_EALOC 64
#define COMPLEXFFTH_EFFTW 128
#define COMPLEXFFTH_ESNGL 256
#define COMPLEXFFTH_EINTL 512
#define COMPLEXFFTH_ESIGN 1024

#define COMPLEXFFTH_MSGENULL "Null pointer"
#define COMPLEXFFTH_MSGENNUL "Non-null pointer"
#define COMPLEXFFTH_MSGESIZE "Invalid input size"
#define COMPLEXFFTH_MSGESZMM "Size mismatch"
#define COMPLEXFFTH_MSGESLEN "Invalid/mismatched sequence lengths"
#define COMPLEXFFTH_MSGESAME "Input/Output data vectors are the same"
#define COMPLEXFFTH_MSGEFFTW "Error in FFTW"
#define COMPLEXFFTH_MSGEALOC "Memory allocation failed"
#define COMPLEXFFTH_MSGESNGL "FFTW library is not single-precision"
#define COMPLEXFFTH_MSGEINTL "Error in Intel FFT library"
#define COMPLEXFFTH_MSGESIGN "Unknown sign of transform in plan"

/**** </lalErrTable> */
/**** <lalLaTeX>
 * 
 * \subsection*{Structures}
 * 
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct tagComplexFFTPlan ComplexFFTPlan;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * 
 * This structure contains the parameters necessary for performing an FFT of a
 * given size and direction.  The contents should not be manually adjusted.
 * 
 * \newpage\input{ComplexFFTC}
 * \newpage\input{ComplexFFTTestC}
 **** </lalLaTeX> */

/* #define KEEP_OLD_COMPLEX_FFT */

void
LALCreateForwardComplexFFTPlan(
    LALStatus       *status,
    ComplexFFTPlan **plan,
    UINT4            size,
    INT4             measure
    );

void
LALCreateReverseComplexFFTPlan(
    LALStatus       *status,
    ComplexFFTPlan **plan,
    UINT4            size,
    INT4             measure
    );

void
LALDestroyComplexFFTPlan (
    LALStatus       *status,
    ComplexFFTPlan **plan
    );

void
LALCOMPLEX8VectorFFT (
    LALStatus      *status,
    COMPLEX8Vector *output,
    COMPLEX8Vector *input,
    ComplexFFTPlan *plan
    );

#ifdef KEEP_OLD_FFT
#define KEEP_OLD_COMPLEX_FFT
#endif

#ifdef KEEP_OLD_COMPLEX_FFT

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

#endif /* KEEP_OLD_COMPLEX_FFT */

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _COMPLEXFFT_H */
