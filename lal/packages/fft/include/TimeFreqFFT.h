/**** <lalVerbatim file="TimeFreqFFTHV">
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{TimeFreqFFT.h}}
 * \label{s:TimeFreqFFT.h}
 * 
 * Performs real-to-complex and complex-to-real FFTs.
 * 
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/TimeFreqFFT.h>
 * \end{verbatim}
 * 
 * \noindent Perform time-to-frequency and frequency-to-time fast Fourier
 * transforms.
 * 
 **** </lalLaTeX> */

#ifndef _TIMEFREQFFT_H
#define _TIMEFREQFFT_H

#include <lal/LALDatatypes.h>
#include <lal/ComplexFFT.h>
#include <lal/RealFFT.h>
#include <lal/Window.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( TIMEFREQFFTH, "$Id$" );

/**** <lalLaTeX>
 * \subsection*{Error conditions}
 **** </lalLaTeX> */
/**** <lalErrTable> */

#define TIMEFREQFFTH_ENULL 1
#define TIMEFREQFFTH_ESIZE 2
#define TIMEFREQFFTH_ERATE 4
#define TIMEFREQFFTH_ESIGN 4
#define TIMEFREQFFTH_EALLOC 16
#define TIMEFREQFFTH_EPOSARG 32
#define TIMEFREQFFTH_EMALLOC 64
#define TIMEFREQFFTH_EINCOMP 128

#define TIMEFREQFFTH_MSGENULL "Null pointer"
#define TIMEFREQFFTH_MSGESIZE "Invalid size"
#define TIMEFREQFFTH_MSGERATE "Invalid rate"
#define TIMEFREQFFTH_MSGESIGN "Incorrect plan sign"
#define TIMEFREQFFTH_MSGEALLOC "Pointer has already been allocated, should be null"
#define TIMEFREQFFTH_MSGEPOSARG "Argument must be positive"
#define TIMEFREQFFTH_MSGEMALLOC "Malloc failure"
#define TIMEFREQFFTH_MSGEINCOMP "Incompatible arguments"

/**** </lalErrTable> */

typedef struct tagRealDFTParams
{
  WindowType               windowType;
  REAL4Vector              *window;
  REAL4                    sumofsquares;
  RealFFTPlan              *plan;
}
RealDFTParams;

typedef enum
{
  useMean,
  useMedian,
  useUnity
} AvgSpecMethod;


/**** <lalLaTeX>
 * \newpage\input{TimeFreqFFTC}
 * \newpage\input{TimeFreqFFTTestC}
 **** </lalLaTeX> */

void
LALCreateRealDFTParams ( 
        LALStatus                         *status, 
        RealDFTParams                  **dftParams, 
        LALWindowParams                   *params,
        INT2                           sign
        );

void
LALDestroyRealDFTParams (
        LALStatus                        *status, 
        RealDFTParams                 **dftParams
        );

void
LALTimeFreqRealFFT(
    LALStatus               *status,
    COMPLEX8FrequencySeries *freq,
    REAL4TimeSeries         *time,
    RealFFTPlan             *plan
    );

void
LALFreqTimeRealFFT(
    LALStatus               *status,
    REAL4TimeSeries         *time,
    COMPLEX8FrequencySeries *freq,
    RealFFTPlan             *plan
    );

void
LALRealAverageSpectrum(
    LALStatus            *status,
    REAL4FrequencySeries *fSeries,
    REAL4TimeSeries      *tSeries,
    RealDFTParams        *params,
    AvgSpecMethod              method
    );

void
LALTimeFreqComplexFFT(
    LALStatus               *status,
    COMPLEX8FrequencySeries *freq,
    COMPLEX8TimeSeries      *time,
    ComplexFFTPlan          *plan
    );

void
LALFreqTimeComplexFFT(
    LALStatus               *status,
    COMPLEX8TimeSeries      *time,
    COMPLEX8FrequencySeries *freq,
    ComplexFFTPlan          *plan
    );

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _TIMEFREQFFT_H */
