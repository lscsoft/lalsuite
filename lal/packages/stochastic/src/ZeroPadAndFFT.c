/************************************ <lalVerbatim file="ZeroPadAndFFTCV">
Author: UTB Relativity Group; contact J. T. Whelan
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{ZeroPadAndFFT.c}}
\label{stochastic:ss:ZeroPadAndFFT.c}

Routines for zero-padding and Fourier transforming a time series.

\subsubsection*{Prototypes}
\input{ZeroPadAndFFTCP}
\index{\texttt{LALSZeroPadAndFFT()}}
% \index{\texttt{LALDZeroPadAndFFT()}}

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
\delta t\,h[k]\,e^{-2\pi k\ell/(2N-1)}
\ .
\end{equation}


\texttt{LALSZeroPadAndFFT()} performs this operaton on a
\texttt{REAL4TimeSeries} of length $N$, zero-padding it to length
$2N-1$ and Fourier-transforming it into a
\texttt{COMPLEX8FrequencySeries} of length $N$.

\subsubsection*{Algorithm}

\texttt{LALSZeroPadAndFFT()} constructs the sequence $\delta
t\,\bar{h}[k]$, and then applies a real-to-complex discrete Fourier
transform from the \texttt{fft} package.

\subsubsection*{Uses}
\noindent
% {\tt LALSZeroPadAndFFT()\/} calls:

\begin{verbatim}
LALSCreateVector()
LALSDestroyVector()
LALFwdRealFFT()
LALUnitMultiply()
strncpy()
\end{verbatim}
% LALUnitPair
% lalSecondUnit

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
#include <lal/RealFFT.h>
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

  REAL4Vector  *hBarData = NULL;
  
  REAL4        *sPtrHBar = NULL;
  REAL4        *sPtrH = NULL;
  REAL4        *sStopPtr = NULL;

  LALUnitPair   unitPair;

  /* initialize status structure */
  INITSTATUS(status, "ZeroPadAndFFT.c", ZEROPADANDFFTC);
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
  if (fullLength != fftPlan->size) {
    ABORT(  status,
	 STOCHASTICCROSSCORRELATIONH_EMMLEN,
	 STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

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

  /* unit manipulation */
  unitPair.unitOne =lalSecondUnit;
  unitPair.unitTwo = input->sampleUnits;
  TRY(LALUnitMultiply(status->statusPtr,&(output->sampleUnits),
		      &unitPair),status);

  /* fill output parameters */
  strncpy(output->name,"Fourier Transform of Zero-Padded Time Series",
	  LALNameLength);
  
  output->epoch = input->epoch;
  output->f0 = input->f0;
  output->deltaF = 1.0 / ( (REAL8) fullLength * input->deltaT );

  /* allocate memory for zero-padded vector */
  TRY(LALSCreateVector(status->statusPtr, &(hBarData), fullLength), status);

  /* initialize pointers for copying */
  sStopPtr = input->data->data + length;

  /* copy data and normalize by deltaT */
  for (sPtrH = input->data->data, sPtrHBar = hBarData->data;
       sPtrH < sStopPtr; 
       ++sPtrH, ++sPtrHBar) 
  {
    *sPtrHBar = *sPtrH * deltaT;
  }

  /* initialize pointers for zero-padding */
  sStopPtr = hBarData->data + fullLength;

  /* zero pad */
  for (sPtrHBar = hBarData->data + length; sPtrHBar < sStopPtr; ++sPtrHBar) 
  {
    *sPtrHBar = 0.0;
  }
    
  /* take DFT */
  LALFwdRealFFT(status->statusPtr, output->data, hBarData, fftPlan);
  /* Can't use TRY because we have memory allocated */
  if (status->statusPtr->statusCode) 
  {
    TRY(LALSDestroyVector(status->statusPtr, &hBarData), status);
    ABORT( status, -1, "Recursive error" );
  }

  /* clean up */
  TRY(LALSDestroyVector(status->statusPtr, &hBarData), status);

  /* normal exit*/
  DETATCHSTATUSPTR(status);
  RETURN(status);


}
