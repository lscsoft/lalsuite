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
 *
 * \newpage\subsection{XLAL Functions}
 *
 * \subsubsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/RealFFT.h>
 *
 * REAL4FFTPlan * XLALCreateREAL4FFTPlan( UINT4 size, int fwdflg, int measurelvl );
 * REAL4FFTPlan * XLALCreateForwardREAL4FFTPlan( UINT4 size, int measurelvl );
 * REAL4FFTPlan * XLALCreateReverseREAL4FFTPlan( UINT4 size, int measurelvl );
 * void XLALDestroyREAL4FFTPlan( REAL4FFTPlan *plan );
 * 
 * int XLALREAL4ForwardFFT( COMPLEX8Vector *output, REAL4Vector *input,
 *     REAL4FFTPlan *plan );
 * int XLALREAL4ReverseFFT( REAL4Vector *output, COMPLEX8Vector *input,
 *     REAL4FFTPlan *plan );
 * int XLALREAL4VectorFFT( REAL4Vector *output, REAL4Vector *input,
 *     REAL4FFTPlan *plan );
 * int XLALREAL4PowerSpectrum( REAL4Vector *spec, REAL4Vector *data,
 *     REAL4FFTPlan *plan );
 *
 * REAL8FFTPlan * XLALCreateREAL8FFTPlan( UINT4 size, int fwdflg, int measurelvl );
 * REAL8FFTPlan * XLALCreateForwardREAL8FFTPlan( UINT4 size, int measurelvl );
 * REAL8FFTPlan * XLALCreateReverseREAL8FFTPlan( UINT4 size, int measurelvl );
 * void XLALDestroyREAL8FFTPlan( REAL8FFTPlan *plan );
 * 
 * int XLALREAL8ForwardFFT( COMPLEX16Vector *output, REAL8Vector *input,
 *     REAL8FFTPlan *plan );
 * int XLALREAL8ReverseFFT( REAL8Vector *output, COMPLEX16Vector *input,
 *     REAL8FFTPlan *plan );
 * int XLALREAL8VectorFFT( REAL8Vector *output, REAL8Vector *input,
 *     REAL8FFTPlan *plan );
 * int XLALREAL8PowerSpectrum( REAL8Vector *spec, REAL8Vector *data,
 *     REAL8FFTPlan *plan );
 * \end{verbatim}
 * \idx{XLALCreateREAL4FFTPlan}
 * \idx{XLALCreateForwardREAL4FFTPlan}
 * \idx{XLALCreateReverseREAL4FFTPlan}
 * \idx{XLALDestroyREAL4FFTPlan}
 * \idx{XLALREAL4ForwardFFT}
 * \idx{XLALREAL4ReverseFFT}
 * \idx{XLALREAL4VectorFFT}
 * \idx{XLALREAL4PowerSpectrum}
 * \idx{XLALCreateREAL8FFTPlan}
 * \idx{XLALCreateForwardREAL8FFTPlan}
 * \idx{XLALCreateReverseREAL8FFTPlan}
 * \idx{XLALDestroyREAL8FFTPlan}
 * \idx{XLALREAL8ForwardFFT}
 * \idx{XLALREAL8ReverseFFT}
 * \idx{XLALREAL8VectorFFT}
 * \idx{XLALREAL8PowerSpectrum}
 * 
 * \subsubsection*{Description}
 *
 * The \verb+REAL4+ routines are described below.  These use single-precision
 * FFTs, i.e., they convert \verb+REAL4Vector+s into \verb+COMPLEX8Vector+s
 * and vice-versa.  The \verb+REAL8+ versions of the routines are the same
 * but they are double-precision versions, i.e., they convert
 * \verb+REAL8Vector+s into \verb+COMPLEX16Vector+s.
 *
 * The routine \verb+XLALCreateREAL4FFTPlan+ creates a \verb+REAL4FFTPlan+
 * structure to perform FFTs of vectors of length \verb+size+.  If 
 * \verb+fwdflg+ is non-zero then the plan is created to perform forward
 * (real-to-complex) FFTs with a negative exponential sign.  Otherwise
 * the plan is created to perform reverse (complex-to-real) FFTs with a
 * positive exponential sign.  The value of \verb+measurelvl+ determines
 * how much optimization of the plan FFTW will do with the most optimization
 * taking the most amount of time.  Reasonable values for \verb+measurelvl+
 * would be 0 for the fasted plan creation (FFTW does not measure the speed
 * of any transform with this level but rather estimates which plan will
 * be the fastet) or 1 to measure a few likely plans to determine the fastest.
 *
 * \verb+XLALCreateForwardREAL4FFTPlan+ is equivalent to
 * \verb+XLALCreateREAL4FFTPlan+ with \verb+fwdflg+ set to 1.
 * \verb+XLALCreateReverseREAL4FFTPlan+ is equivalent to
 * \verb+XLALCreateREAL4FFTPlan+ with \verb+fwdflg+ set to 0.
 *
 * \verb+XLALDestroyREAL4FFTPlan+ is used to destroy the plan, freeing all
 * memory that was allocated in the structure as well as the structure
 * itself.  It can be used on either forward or reverse plans.
 *
 * \verb+XLALREAL4ForwardFFT+ and 
 * \verb+XLALREAL4ReverseFFT+ perform forward (real to complex) and
 * reverse (complex to real) transforms respectively.  The plan supplied
 * to these routines must be correctly generated for the direction of the
 * transform.  I.e., \verb+XLALREAL4ForwardFFT+ cannot be supplied with
 * a plan generated by \verb+XLALCreateReverseREAL4FFTPlan+.
 *
 * \verb+XLALREAL4VectorFFT+ is a low-level routine that transforms
 * a real vector to a half-complex real vector (with a forward plan) or
 * a half-complex real vector to a real vector (with a reverse plan).
 * If you're not sure what this means, don't use this routine.
 * The input and output vectors (and their data) must be distinct pointers.
 *
 * \verb+XLALREAL4PowerSpectrum+ computes a real power spectrum of the
 * input real vector and a forward FFT plan.
 *
 * \subsubsection*{Return Values}
 *
 * Upon success, 
 * \verb+XLALCreateREAL4FFTPlan+,
 * \verb+XLALCreateForwardREAL4FFTPlan+, and 
 * \verb+XLALCreateReverseREAL4FFTPlan+ return a pointer to a newly-allocated
 * FFT plan.  Upon failure, they return a \verb+NULL+ pointer and set 
 * \verb+xlalErrno+ to one of the following values:
 * \verb+XLAL_EBADLEN+ if \verb+size+ is not greater than zero,
 * \verb+XLAL_ENOMEM+ if a memory allocation failed, or
 * \verb+XLAL_EFAILED+ if the FFTW plan creation routine failed.
 *
 * \verb+XLALDestroyREAL4FFTPlan+ does not return any value but, upon
 * failure, it will set \verb+xlalErrno+ to one of the following values:
 * \verb+XLAL_EFAULT+ if the routine is provided a \verb+NULL+ pointer, or
 * \verb+XLAL_EINVAL+ if the contents of the plan are invalid (e.g., if the
 * routine is provided a plan that had been previously destroyed).
 *
 * \verb+XLALREAL4ForwardFFT+,
 * \verb+XLALREAL4ReverseFFT+,
 * \verb+XLALREAL4VectorFFT+, and
 * \verb+XLALREAL4PowerSpectrum+ return the value 0 upon succes; upon
 * failure they return \verb+XLAL_FAILURE+ and set \verb+xlalErrno+ to
 * one of the following values:
 * \verb+XLAL_EFAULT+ if one of the input pointers is \verb+NULL+,
 * \verb+XLAL_EINVAL+ if the input, output, or plan structures appears
 * invalid or if the routine is passed a plan for the wrong transform
 * directions or if the input and output data pointers are not distinct
 * for \verb+XLALREAL4VectorFFT+,
 * \verb+XLAL_EBADLEN+ if the input and output vectors and the plan have
 * incompatible lengths,
 * \verb+XLAL_ENOMEM+ if a memory allocation of temporary internal memory
 * fails.
 *
 * As before, the \verb+REAL8+ versions of these routines behave the
 * same way but for double-precision transforms.
 *
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
