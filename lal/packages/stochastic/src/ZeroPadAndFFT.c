/************************************ <lalVerbatim file="ZeroPadAndFFTCV">
Author: UTB Relativity Group; contact J. T. Whelan
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{ZeroPadAndFFT.c}}
\label{stochastic:ss:ZeroPadAndFFT.c}

{\bf {\Large WARNING} The functionality of this module has been expanded
and modified, so the documentation is not yet complete and/or correct.}

Routines for zero-padding and Fourier transforming a time series.

\subsubsection*{Prototypes}
\input{ZeroPadAndFFTCP}
\index{\texttt{LALSZeroPadAndFFT()}}
\index{\texttt{LALCZeroPadAndFFT()}}

\subsubsection*{Description}

As described in \ref{stochastic:ss:StochasticCrossCorrelation.c}, data
streams to be cross-correlated need to be zero-padded to the same
length as the optimal filter via
% 
\begin{equation}
\bar{h}[k]=\ 
\left\{ \begin{array}{cl} 
h[k]  &    k = 0, 1, \cdots, N-1 \\ 
0     &    k = -1, -2, \cdots, -(N-1) 
\end{array} 
\right. 
\end{equation} 
%
before being Fourier transformed via
\begin{equation}
\widetilde{h}[\ell] := \sum_{\ell=-(N-1)}^{N-1} 
\delta t\,h[k]\,e^{-i2\pi k\ell/(2N-1)}
\ .
\end{equation}

\texttt{LALSZeroPadAndFFT()} performs this operaton on a
\texttt{REAL4TimeSeries} of length $N$, zero-padding it to length
$2N-1$ and Fourier-transforming it into a
\texttt{COMPLEX8FrequencySeries} of length $N$.

\texttt{LALCZeroPadAndFFT()} performs this operaton on a
\texttt{COMPLEX8TimeSeries} of length $N$, zero-padding it to length
$2N-1$ and Fourier-transforming it into a
\texttt{COMPLEX8FrequencySeries} of length $2N-1$.

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
LALUnitMultiply()
memcpy()
memset()
strncpy()
\end{verbatim}
% LALUnitPair
% lalSecondUnit

{\tt LALSZeroPadAndFFT()\/} calls:

\begin{verbatim}
LALCCreateVector()
LALCDestroyVector()
LALTimeFreqComplexFFT()
LALUnitMultiply()
memcpy()
memset()
strncpy()
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}

\item The Fourier transform is defined to be the discrete
approximation of a continuous Fourier transorm, which makes it $\delta
t$ times the discrete Fourier transform.

\item The Fourier transform of a series of $2N-1$ points is calculated
with the FFTW \cite{stochastic:Frigo:1998} (via the interfaces in the
\texttt{fft} package), which is efficient for products of small
primes, so $N$ should be chosen so that $2N-1$ is a power of 3 or
other suitable number.

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
LALSZeroPadAndFFT(LALStatus                *status, 
                  COMPLEX8FrequencySeries  *output, 
                  const REAL4TimeSeries    *input, 
                  RealFFTPlan              *fftPlan)
/* </lalVerbatim> */
{

  UINT4          length, fullLength; /* fullLength = 2 * length - 1 */

  REAL8         deltaT, deltaF; /* deltaT * deltaF = 1/fullLength */

  REAL4TimeSeries  hBar;

  /* initialize status structure */
  INITSTATUS(status, "LALSZeroPadAndFFT", ZEROPADANDFFTC);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING --------------------------------------------------- */

  /* check that pointer to real timer series for input is non-null */
  ASSERT(input != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that pointer to data member of real time series for input is 
     non-null */
  ASSERT(input->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that length of data member of real time series for input is 
     not equal to zero */
  length = input->data->length;
  ASSERT(length != 0, status,
         STOCHASTICCROSSCORRELATIONH_EZEROLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that pointer to data-data member of real time series for input 
     is non-null */
  ASSERT(input->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that pointer to FFT plan variable is non-null */
  ASSERT(fftPlan != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that length of FFT plan is correct */
  fullLength = 2 * length - 1;
  /* OMITTED -- JC
  if (fullLength != fftPlan->size) {
    ABORT(  status,
         STOCHASTICCROSSCORRELATIONH_EMMLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  */

  /* check that pointer to complex frequency series for output is non-null */
  ASSERT(output != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that pointer to data member of complex frequency series for 
     output is non-null */
  ASSERT(output->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that lengths of data member of real time series for input 
     and data member of complex frequency series for output are equal */
  if(length != output->data->length)
  {
    ABORT(  status,
         STOCHASTICCROSSCORRELATIONH_EMMLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to data-data member of complex frequency series 
     for output is non-null */
  ASSERT(output->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);
  
  /* check that heterodyning frequency of real frequency series for
     the input is equal to zero */
  if (input->f0 != 0)
  {
     ABORT(status,STOCHASTICCROSSCORRELATIONH_ENONZEROHETERO,
           STOCHASTICCROSSCORRELATIONH_MSGENONZEROHETERO); 
  }

  /* check that frequency spacing is positive */
  deltaT = input->deltaT;
  ASSERT(deltaT > 0, status, 
         STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAT,
         STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAT);


  /* EVERYTHING OKAY HERE! -------------------------------------------- */

  /* replicate input for zero-padding */
  hBar = *input;
  hBar.data = NULL;

  /* allocate memory for zero-padded vector */
  TRY(LALSCreateVector(status->statusPtr, &(hBar.data), fullLength), status);

  /* copy data */
  memcpy(hBar.data->data, input->data->data, length*sizeof(REAL4));

  /* zero pad */
  memset(hBar.data->data + length, 0, (fullLength - length)*sizeof(REAL4));
    
  /* take DFT */
  LALTimeFreqRealFFT(status->statusPtr, output, &hBar, fftPlan);
  /* Can't use TRY because we have memory allocated */
  BEGINFAIL( status ) 
    TRY(LALSDestroyVector(status->statusPtr, &(hBar.data)), status);
  ENDFAIL( status ); 

  /* fill output parameters */
  strncpy(output->name,"Fourier Transform of Zero-Padded Time Series",
          LALNameLength);

  /* clean up */
  TRY(LALSDestroyVector(status->statusPtr, &(hBar.data)), status);

  /* normal exit*/
  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* SZeroPadAndFFT() */

/* <lalVerbatim file="ZeroPadAndFFTCP"> */
void 
LALCZeroPadAndFFT(LALStatus                *status, 
                  COMPLEX8FrequencySeries  *output, 
                  const COMPLEX8TimeSeries *input, 
                  ComplexFFTPlan           *fftPlan)
/* </lalVerbatim> */
{

  UINT4          length, fullLength; /* fullLength = 2 * length - 1 */

  REAL8         deltaT, deltaF; /* deltaT * deltaF = 1/fullLength */

  COMPLEX8TimeSeries  hBar;

  /* initialize status structure */
  INITSTATUS(status, "LALCZeroPadAndFFT", ZEROPADANDFFTC);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING --------------------------------------------------- */

  /* check that pointer to real timer series for input is non-null */
  ASSERT(input != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that pointer to data member of real time series for input is 
     non-null */
  ASSERT(input->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that length of data member of real time series for input is 
     not equal to zero */
  length = input->data->length;
  ASSERT(length != 0, status,
         STOCHASTICCROSSCORRELATIONH_EZEROLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that pointer to data-data member of real time series for input 
     is non-null */
  ASSERT(input->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that pointer to FFT plan variable is non-null */
  ASSERT(fftPlan != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that length of FFT plan is correct */
  fullLength = 2 * length - 1;
  /* OMITTED -- JC
  if (fullLength != fftPlan->size) {
    ABORT(  status,
         STOCHASTICCROSSCORRELATIONH_EMMLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  */

  /* check that pointer to complex frequency series for output is non-null */
  ASSERT(output != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that pointer to data member of complex frequency series for 
     output is non-null */
  ASSERT(output->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that lengths of data member of zero-padded time series
     and data member of complex frequency series for output are equal */
  if(fullLength != output->data->length)
  {
    ABORT(  status,
         STOCHASTICCROSSCORRELATIONH_EMMLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to data-data member of complex frequency series 
     for output is non-null */
  ASSERT(output->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* check that frequency spacing is positive */
  deltaT = input->deltaT;
  ASSERT(deltaT > 0, status, 
         STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAT,
         STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAT);


  /* EVERYTHING OKAY HERE! -------------------------------------------- */

  /* replicate input for zero-padding */
  hBar = *input;
  hBar.data = NULL;

  /* allocate memory for zero-padded vector */
  TRY(LALCCreateVector(status->statusPtr, &(hBar.data), fullLength), status);

  /* copy data */
  memcpy(hBar.data->data, input->data->data, length*sizeof(COMPLEX8));

  /* zero pad */
  memset(hBar.data->data + length, 0, (fullLength - length)*sizeof(COMPLEX8));
    
  /* take DFT */
  LALTimeFreqComplexFFT(status->statusPtr, output, &hBar, fftPlan);
  /* Can't use TRY because we have memory allocated */
  BEGINFAIL( status ) 
    TRY(LALCDestroyVector(status->statusPtr, &(hBar.data)), status);
  ENDFAIL( status ); 

  /* fill output parameters */
  strncpy(output->name,"Fourier Transform of Zero-Padded Time Series",
          LALNameLength);

  /* clean up */
  TRY(LALCDestroyVector(status->statusPtr, &(hBar.data)), status);

  /* normal exit*/
  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* CZeroPadAndFFT() */
