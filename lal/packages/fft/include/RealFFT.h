/**** <lalVerbatim file="RealFFTHV">
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{RealFFT.h}}
 * \label{s:RealFFT.h}
 * 
 * Performs real-to-complex and complex-to-real FFTs.
 * 
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/RealFFT.h>
 * \end{verbatim}
 * 
 * \noindent Perform real-to-complex and complex-to-real fast Fourier
 * transforms of vectors, and sequences of vectors using the package
 * FFTW~\cite{fj:1998}.
 * 
 **** </lalLaTeX> */

#ifndef _REALFFT_H
#define _REALFFT_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( REALFFTH, "$Id$" );


/**** <lalLaTeX>
 * \subsection*{Error conditions}
 **** </lalLaTeX> */
/**** <lalErrTable> */

#define REALFFTH_ENULL 1
#define REALFFTH_ENNUL 2
#define REALFFTH_ESIZE 4
#define REALFFTH_ESZMM 8
#define REALFFTH_ESLEN 16
#define REALFFTH_ESAME 32
#define REALFFTH_ESIGN 64
#define REALFFTH_EDATA 128
#define REALFFTH_EALOC 256
#define REALFFTH_EFFTW 512

#define REALFFTH_MSGENULL "Null pointer"
#define REALFFTH_MSGENNUL "Non-null pointer"
#define REALFFTH_MSGESIZE "Invalid input size"
#define REALFFTH_MSGESZMM "Size mismatch"
#define REALFFTH_MSGESLEN "Invalid/mismatched sequence lengths"
#define REALFFTH_MSGESAME "Input/Output data vectors are the same"
#define REALFFTH_MSGESIGN "Incorrect plan sign"
#define REALFFTH_MSGEDATA "Bad input data: DC/Nyquist should be real"
#define REALFFTH_MSGEALOC "Memory allocation failed"
#define REALFFTH_MSGEFFTW "Error in FFTW"

/**** </lalErrTable> */
/**** <lalLaTeX>
 * 
 * \subsection*{Structures}
 * 
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct tagRealFFTPlan RealFFTPlan;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * 
 * This structure contains the parameters necessary for performing an FFT of a
 * given size and direction.  The contents should not be manually adjusted.
 * 
 * \newpage\input{RealFFTC}
 * \newpage\input{RealFFTTestC}
 **** </lalLaTeX> */

/* #define KEEP_OLD_REAL_FFT */

void
LALCreateForwardRealFFTPlan(
    LALStatus    *status,
    RealFFTPlan **plan,
    UINT4         size,
    INT4          measure
    );

void
LALCreateReverseRealFFTPlan(
    LALStatus    *status,
    RealFFTPlan **plan,
    UINT4         size,
    INT4          measure
    );

void
LALDestroyRealFFTPlan(
    LALStatus    *status,
    RealFFTPlan **plan
    );

void
LALForwardRealFFT(
    LALStatus      *status,
    COMPLEX8Vector *output,
    REAL4Vector    *input,
    RealFFTPlan    *plan
    );

void
LALReverseRealFFT(
    LALStatus      *status,
    REAL4Vector    *output,
    COMPLEX8Vector *input,
    RealFFTPlan    *plan
    );

void
LALRealPowerSpectrum (
    LALStatus   *status,
    REAL4Vector *spec,
    REAL4Vector *data,
    RealFFTPlan *plan
    );

void
LALREAL4VectorFFT(
    LALStatus   *status,
    REAL4Vector *output,
    REAL4Vector *input,
    RealFFTPlan *plan
    );

/** OLD ROUTINES **/
#ifdef KEEP_OLD_REAL_FFT

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

#endif /*  KEEP_OLD_REAL_FFT */

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _REALFFT_H */
