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
#define REALFFTH_ESNGL 1024
#define REALFFTH_EINTL 2048

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
#define REALFFTH_MSGESNGL "FFTW library is not single-precision"
#define REALFFTH_MSGEINTL "Error in Intel FFT library"


/**** </lalErrTable> */
/**** <lalLaTeX>
 * 
 * \subsection*{Structures}
 * 
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct tagREAL4FFTPlan REAL4FFTPlan;
typedef struct tagREAL8FFTPlan REAL8FFTPlan;
#define tagRealFFTPlan tagREAL4FFTPlan
#define RealFFTPlan REAL4FFTPlan
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * 
 * This structure contains the parameters necessary for performing an FFT of a
 * given size and direction.  The contents should not be manually adjusted.
 * 
 * \newpage\input{RealFFTC}
 * \newpage\input{RealFFTTestC}
 **** </lalLaTeX> */

/*
 *
 * XLAL REAL4 functions
 *
 */

REAL4FFTPlan * XLALCreateREAL4FFTPlan( UINT4 size, int fwdflg, int measurelvl );
REAL4FFTPlan * XLALCreateForwardREAL4FFTPlan( UINT4 size, int measurelvl );
REAL4FFTPlan * XLALCreateReverseREAL4FFTPlan( UINT4 size, int measurelvl );
void XLALDestroyREAL4FFTPlan( REAL4FFTPlan *plan );

int XLALREAL4ForwardFFT( COMPLEX8Vector *output, REAL4Vector *input,
    REAL4FFTPlan *plan );
int XLALREAL4ReverseFFT( REAL4Vector *output, COMPLEX8Vector *input,
    REAL4FFTPlan *plan );
int XLALREAL4VectorFFT( REAL4Vector *output, REAL4Vector *input,
    REAL4FFTPlan *plan );
int XLALREAL4PowerSpectrum( REAL4Vector *spec, REAL4Vector *data,
    REAL4FFTPlan *plan );

/*
 *
 * XLAL REAL8 functions
 *
 */

REAL8FFTPlan * XLALCreateREAL8FFTPlan( UINT4 size, int fwdflg, int measurelvl );
REAL8FFTPlan * XLALCreateForwardREAL8FFTPlan( UINT4 size, int measurelvl );
REAL8FFTPlan * XLALCreateReverseREAL8FFTPlan( UINT4 size, int measurelvl );
void XLALDestroyREAL8FFTPlan( REAL8FFTPlan *plan );

int XLALREAL8ForwardFFT( COMPLEX16Vector *output, REAL8Vector *input,
    REAL8FFTPlan *plan );
int XLALREAL8ReverseFFT( REAL8Vector *output, COMPLEX16Vector *input,
    REAL8FFTPlan *plan );
int XLALREAL8VectorFFT( REAL8Vector *output, REAL8Vector *input,
    REAL8FFTPlan *plan );
int XLALREAL8PowerSpectrum( REAL8Vector *spec, REAL8Vector *data,
    REAL8FFTPlan *plan );

/*
 *
 * LAL REAL4 functions
 *
 */

void
LALCreateForwardREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan,
    UINT4         size,
    INT4          measure
    );
#define LALCreateForwardRealFFTPlan LALCreateForwardREAL4FFTPlan

void
LALCreateReverseREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan,
    UINT4         size,
    INT4          measure
    );
#define LALCreateReverseRealFFTPlan LALCreateReverseREAL4FFTPlan

void
LALDestroyREAL4FFTPlan(
    LALStatus    *status,
    REAL4FFTPlan **plan
    );
#define LALDestroyRealFFTPlan LALDestroyREAL4FFTPlan

void
LALForwardREAL4FFT(
    LALStatus      *status,
    COMPLEX8Vector *output,
    REAL4Vector    *input,
    REAL4FFTPlan    *plan
    );
#define LALForwardRealFFT LALForwardREAL4FFT

void
LALReverseREAL4FFT(
    LALStatus      *status,
    REAL4Vector    *output,
    COMPLEX8Vector *input,
    REAL4FFTPlan    *plan
    );
#define LALReverseRealFFT LALReverseREAL4FFT

void
LALREAL4PowerSpectrum (
    LALStatus   *status,
    REAL4Vector *spec,
    REAL4Vector *data,
    REAL4FFTPlan *plan
    );
#define LALRealPowerSpectrum LALREAL4PowerSpectrum

void
LALREAL4VectorFFT(
    LALStatus   *status,
    REAL4Vector *output,
    REAL4Vector *input,
    REAL4FFTPlan *plan
    );

/*
 *
 * LAL REAL8 functions
 *
 */

void
LALCreateForwardREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan,
    UINT4         size,
    INT4          measure
    );

void
LALCreateReverseREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan,
    UINT4         size,
    INT4          measure
    );

void
LALDestroyREAL8FFTPlan(
    LALStatus    *status,
    REAL8FFTPlan **plan
    );

void
LALForwardREAL8FFT(
    LALStatus      *status,
    COMPLEX16Vector *output,
    REAL8Vector    *input,
    REAL8FFTPlan    *plan
    );

void
LALReverseREAL8FFT(
    LALStatus      *status,
    REAL8Vector    *output,
    COMPLEX16Vector *input,
    REAL8FFTPlan    *plan
    );

void
LALREAL8PowerSpectrum (
    LALStatus   *status,
    REAL8Vector *spec,
    REAL8Vector *data,
    REAL8FFTPlan *plan
    );

void
LALREAL8VectorFFT(
    LALStatus   *status,
    REAL8Vector *output,
    REAL8Vector *input,
    REAL8FFTPlan *plan
    );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _REALFFT_H */
