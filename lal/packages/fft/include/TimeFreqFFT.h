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

#define TIMEFREQFFTH_MSGENULL "Null pointer"
#define TIMEFREQFFTH_MSGESIZE "Invalid size"
#define TIMEFREQFFTH_MSGERATE "Invalid rate"
#define TIMEFREQFFTH_MSGESIGN "Incorrect plan sign"

/**** </lalErrTable> */
/**** <lalLaTeX>
 * \newpage\input{TimeFreqFFTC}
 * \newpage\input{TimeFreqFFTTestC}
 **** </lalLaTeX> */

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
