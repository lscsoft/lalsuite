/************************************ <lalVerbatim file="ZeroPadAndFFTCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{ZeroPadAndFFT.c}}
\label{stochastic:ss:ZeroPadAndFFT.c}

Routines for zero-padding and Fourier transforming a time series.

\subsubsection*{Prototypes}
\idx{LALSZeroPadAndFFT()}
\idx{LALCZeroPadAndFFT()}
\input{ZeroPadAndFFTCP}

\subsubsection*{Description}

As described in \ref{stochastic:ss:StochasticCrossCorrelation.c}, data
streams to be cross-correlated need to be zero-padded to the same
length as the optimal filter via
% 
\begin{equation}
\bar{h}[k]=\ 
\left\{ \begin{array}{cl} 
w[k] h[k]  &    k = 0, \ldots, N-1 \\ 
0     &    k = N \ldots, M-1 
\end{array} 
\right. 
\end{equation} 
%
(where $w[k]$ is a windowing function)
before being Fourier transformed via
\begin{equation}
\widetilde{h}[\ell] := \sum_{\ell=0}^{M-1} 
\delta t\,h[k]\,e^{-i2\pi k\ell/M}
\ .
\end{equation}

\texttt{LALSZeroPadAndFFT()} performs this operaton on a
\texttt{REAL4TimeSeries} of length $N$, zero-padding it to length
$M$ and Fourier-transforming it into a
\texttt{COMPLEX8FrequencySeries} of length $[M/2]+1$.

\texttt{LALCZeroPadAndFFT()} performs this operaton on a
\texttt{COMPLEX8TimeSeries} of length $N$, zero-padding it to length
$M$ and Fourier-transforming it into a
\texttt{COMPLEX8FrequencySeries} of length $M$.

\subsubsection*{Algorithm}

\texttt{LALSZeroPadAndFFT()} constructs the sequence $\bar{h}[k]$, and
then applies a real-to-complex time-to-frequency discrete Fourier
transform from the \texttt{fft} package.

\texttt{LALCZeroPadAndFFT()} constructs the sequence $\bar{h}[k]$, and
then applies a complex-to-complex time-to-frequency discrete Fourier
transform from the \texttt{fft} package.

\subsubsection*{Uses}
\noindent
{\tt LALSZeroPadAndFFT()\/} calls:

\begin{verbatim}
LALSCreateVector()
LALSDestroyVector()
LALTimeFreqRealFFT()
memset()
strncpy()
\end{verbatim}

{\tt LALSZeroPadAndFFT()\/} calls:

\begin{verbatim}
LALCCreateVector()
LALCDestroyVector()
LALTimeFreqComplexFFT()
memset()
strncpy()
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}

\item The Fourier transform is defined to be the discrete
approximation of a continuous Fourier transorm, which makes it $\delta
t$ times the discrete Fourier transform.

\item The Fourier transform of a series of $M$ points is calculated
  with the FFTW \cite{stochastic:Frigo:1998} (via the interfaces in
  the \texttt{fft} package), which is efficient for products of small
  primes, so $M$ should be chosen to have this property.  The minimum
  value, $2N-1$, is odd and can thus be at best a power of 3.
  Additionally, if $2N-1$ is a convenient number, $N$ will likely not
  be, which is one reason it might be convenient to work with $M=2N$
  instead.

\item \texttt{LALCZeroPadAndFFT()} inherits its behavior from
  \texttt{LALTimeFreqComplexFFT()}, which currently does not use the
  initial phase of the reference oscillator.  The calling routine must
  therefore remove the effects of this phase explicitly in order to
  obtain the band-limited FFT of the unheterodyned data.

\item The output units are determined from the input units, but under
  normal circumstances in the context of a stochastic background
  search, we will have
  \begin{eqnarray}
    {} [h(t)] &=& \textrm{count}\\
    {} [\widetilde{\bar{h}}(f)] &:=& [h(t)] \,\textrm{Hz}^{-1}
    = \textrm{count}\,\textrm{Hz}^{-1}
  \end{eqnarray}


\end{itemize}

\vfill{\footnotesize\input{ZeroPadAndFFTCV}}

******************************************************* </lalLaTeX> */
/**************************************** <lalLaTeX file="ZeroPadAndFFTCB">
\bibitem{stochastic:Frigo:1998}
  M. Frigo and S. G. Johnson,
  \textit{FFTW User's Manual},
  (Massachusetts Institute of Technology, Cambridge, USA, 1998).
  URL: \href{http://www.fftw.org/doc/}{\texttt{http://www.fftw.org/doc/}}
******************************************************* </lalLaTeX> */ 

#include <lal/LALStdlib.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/AVFactories.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <string.h>

NRCSID(ZEROPADANDFFTC,
       "$Id$");

/* <lalVerbatim file="ZeroPadAndFFTCP"> */
void
LALSZeroPadAndFFT(
    LALStatus                *status,
    COMPLEX8FrequencySeries  *output,
    const REAL4TimeSeries    *input,
    SZeroPadAndFFTParameters *parameters)
/* </lalVerbatim> */
{
  UINT4 length, fullLength;
  REAL8 deltaT;
  REAL4TimeSeries  hBar;
  REAL4 *sPtr, *sStopPtr, *hBarPtr, *windowPtr;

  /* initialize status structure */
  INITSTATUS(status, "LALSZeroPadAndFFT", ZEROPADANDFFTC);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING --------------------------------------------------- */

  /* check that pointer to real timer series for input is non-null */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of real time series for input is 
   * non-null */
  ASSERT(input->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of data member of real time series for input is 
   * not equal to zero */
  length = input->data->length;
  ASSERT(length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that pointer to data-data member of real time series for input 
   * is non-null */
  ASSERT(input->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to parameter structure is non-null */
  ASSERT(parameters != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to FFT plan parameter is non-null */
  ASSERT(parameters->fftPlan != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* These checks are only relevant if a window function is specified */
  if (parameters->window != NULL)
  {
    /* check that pointer to data member of window function is non-null */
    ASSERT(parameters->window->data != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* check that window function is same length as input time series */
    if (parameters->window->length != length)
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
          STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
    }
  }

  /* check that zero-padded output is not shorter than input */
  fullLength = parameters->length;
  if (fullLength < length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to complex frequency series for output is non-null */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of complex frequency series for
   * output is non-null */
  ASSERT(output->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of complex frequency series for output
     is consistent with length of zero-padded input */
  if (fullLength/2 + 1 != output->data->length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to data-data member of complex frequency series 
   * for output is non-null */
  ASSERT(output->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
  
  /* check that heterodyning frequency of real frequency series for
   * the input is equal to zero */
  if (input->f0 != 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENONZEROHETERO, \
           STOCHASTICCROSSCORRELATIONH_MSGENONZEROHETERO);
  }

  /* check that frequency spacing is positive */
  deltaT = input->deltaT;
  ASSERT(deltaT > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAT, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAT);

  /* EVERYTHING OKAY HERE! -------------------------------------------- */

  /* replicate input for zero-padding */
  hBar = *input;
  hBar.data = NULL;

  /* allocate memory for zero-padded vector */
  TRY(LALSCreateVector(status->statusPtr, &(hBar.data), fullLength), status);

  if (parameters->window == NULL)
  {
    sPtr = memcpy(hBar.data->data, input->data->data, length * sizeof(REAL4));
    if (sPtr != hBar.data->data)
    {
      TRY(LALSDestroyVector(status->statusPtr, &(hBar.data)), status);
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EMEMORY, \
          STOCHASTICCROSSCORRELATIONH_MSGEMEMORY );
    }
  }
  else 
  {
    /* window data */
    sStopPtr = input->data->data + input->data->length;
    for (sPtr = input->data->data, hBarPtr = hBar.data->data, \
        windowPtr = parameters->window->data ; sPtr < sStopPtr ; \
        ++sPtr, ++hBarPtr, ++windowPtr)
    {
      *(hBarPtr) = *(sPtr) * *(windowPtr);
    }
  }

  /* zero pad */
  sPtr = memset(hBar.data->data + length, 0, \
      (fullLength - length) * sizeof(REAL4));

  if (sPtr != hBar.data->data + length)
  {
    TRY(LALSDestroyVector(status->statusPtr, &(hBar.data)), status);
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMEMORY, \
        STOCHASTICCROSSCORRELATIONH_MSGEMEMORY);
  }
  
  /* take DFT */
  LALTimeFreqRealFFT(status->statusPtr, output, &hBar, parameters->fftPlan);

  /* Can't use TRY because we have memory allocated */
  BEGINFAIL(status)
    TRY(LALSDestroyVector(status->statusPtr, &(hBar.data)), status);
  ENDFAIL(status); 

  /* fill output parameters */
  strncpy(output->name, "Fourier Transform of Zero-Padded Time Series", \
      LALNameLength);

  /* clean up */
  TRY(LALSDestroyVector(status->statusPtr, &(hBar.data)), status);

  /* normal exit*/
  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* SZeroPadAndFFT() */

/* <lalVerbatim file="ZeroPadAndFFTCP"> */
void
LALCZeroPadAndFFT(
    LALStatus                *status,
    COMPLEX8FrequencySeries  *output, 
    const COMPLEX8TimeSeries *input, 
    CZeroPadAndFFTParameters *parameters)
/* </lalVerbatim> */
{
  UINT4 length, fullLength;
  REAL8 deltaT;
  COMPLEX8TimeSeries  hBar;
  COMPLEX8 *cPtr, *cStopPtr, *hBarPtr;
  REAL4 *windowPtr;

  /* initialize status structure */
  INITSTATUS(status, "LALCZeroPadAndFFT", ZEROPADANDFFTC);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING --------------------------------------------------- */

  /* check that pointer to real time series for input is non-null */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of real time series for input is 
   * non-null */
  ASSERT(input->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of data member of real time series for input is 
   * not equal to zero */
  length = input->data->length;
  ASSERT(length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that pointer to data-data member of real time series for input 
   * is non-null */
  ASSERT(input->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to parameter structure is non-null */
  ASSERT(parameters != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to FFT plan parameter is non-null */
  ASSERT(parameters->fftPlan != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* These checks are only relevant if a window function is specified */
  if (parameters->window != NULL)
  {
    /* check that pointer to data member of window function is non-null */
    ASSERT(parameters->window->data != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    
    /* check that window function is same length as input time series */
    if (parameters->window->length != length)
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
          STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
    }
  }    

  /* check that zero-padded output is not shorter than input */
  fullLength = parameters->length;
  if (fullLength < length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to complex frequency series for output is non-null */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of complex frequency series for 
   * output is non-null */
  ASSERT(output->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of complex frequency series for output
   * is consistent with length of zero-padded input */
  if (fullLength != output->data->length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to data-data member of complex frequency series 
   * for output is non-null */
  ASSERT(output->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that frequency spacing is positive */
  deltaT = input->deltaT;
  ASSERT(deltaT > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAT, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAT);

  /* EVERYTHING OKAY HERE! -------------------------------------------- */

  /* replicate input for zero-padding */
  hBar = *input;
  hBar.data = NULL;

  /* allocate memory for zero-padded vector */
  TRY(LALCCreateVector(status->statusPtr, &(hBar.data), fullLength), status);

  if (parameters->window == NULL) 
  {
    cPtr = memcpy(hBar.data->data, input->data->data, \
        length * sizeof(COMPLEX8));
    if (cPtr != hBar.data->data)
    {
      TRY(LALCDestroyVector(status->statusPtr, &(hBar.data)), status);
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EMEMORY, \
          STOCHASTICCROSSCORRELATIONH_MSGEMEMORY );
    }
  }
  else 
  {
    /* window data */
    cStopPtr = input->data->data + input->data->length;
    for (cPtr = input->data->data, hBarPtr = hBar.data->data, \
        windowPtr = parameters->window->data ; cPtr < cStopPtr ; \
        ++cPtr, ++hBarPtr, ++windowPtr )
    {
      hBarPtr->re = cPtr->re * *(windowPtr);
      hBarPtr->im = cPtr->im * *(windowPtr);
    }
  }

  /* zero pad */
  cPtr = memset(hBar.data->data + length, 0, \
      (fullLength - length) * sizeof(COMPLEX8));

  if (cPtr != hBar.data->data + length)
  {
    TRY(LALCDestroyVector(status->statusPtr, &(hBar.data)), status);
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMEMORY, \
        STOCHASTICCROSSCORRELATIONH_MSGEMEMORY );
  }

  /* take DFT */
  LALTimeFreqComplexFFT(status->statusPtr, output, &hBar, parameters->fftPlan);

  /* Can't use TRY because we have memory allocated */
  BEGINFAIL(status) 
    TRY(LALCDestroyVector(status->statusPtr, &(hBar.data)), status);
  ENDFAIL(status); 

  /* fill output parameters */
  strncpy(output->name, "Fourier Transform of Zero-Padded Time Series", \
      LALNameLength);

  /* clean up */
  TRY(LALCDestroyVector(status->statusPtr, &(hBar.data)), status);

  /* normal exit*/
  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* CZeroPadAndFFT() */
