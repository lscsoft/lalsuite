/************************** <lalVerbatim file="StochasticCrossCorrelationCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu (original by S. Drasco)
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{StochasticCrossCorrelation.c}}
\label{stochastic:ss:StochasticCrossCorrelation.c}

Calculates the value of the standard optimally-filtered 
cross-correlation statistic for stochastic background searches.

\subsubsection*{Prototypes}
\idx{LALStochasticCrossCorrelationStatistic()}
\idx{LALStochasticHeterodynedCrossCorrelationStatistic()}
\idx{LALStochasticCrossCorrelationSpectrum()}
\input{StochasticCrossCorrelationCP}

\subsubsection*{Description}

\subsubsection*{\texttt{LALStochasticCrossCorrelationStatistic()}}

The default version of the function, for handling non-heterodyned
data, calculates the value of the standard optimally-filtered
cross-correlation statistic
% 
\begin{eqnarray} 
Y
&:=&\int_{t_0}^{t_0+T} dt_1\int_{t_0}^{t_0+T} dt_2\,
w_1(t_1)\, h_1(t_1)\, Q(t_1-t_2)\, w_2(t_2)\, h_2(t_2) \nonumber \\
&\approx& \sum_{j=0}^{N-1}\delta t\sum_{k=0}^{N-1}\delta t\,
w_1[j]\, h_1[j]\, Q[j-k]\, w_2[j]\, h_2[k] \nonumber \\
&=& \sum_{\ell=0}^{M-1} \delta f\,
\widetilde{\bar{h}}_{1}[\ell]^* \,\widetilde{Q}[\ell]\,
\widetilde{\bar{h}}_{2}[\ell],
\label{stochastic:e:ymax}
\end{eqnarray}
%
where the sampling period is $\delta t=T/N$, the frequency spacing is
$\delta f = [M\delta t]^{-1}$, the tilde indicates a discrete
Fourier transform normalized to approximate the continuous Fourier
transform:
\begin{equation}
\widetilde{Q}[\ell] := \sum_{k=0}^{N-1} \delta t\,
Q[k]\, e^{-i2\pi k\ell/M}
\end{equation}
the asterisk indicates complex conjugation, and the overbar indicates
windowing and zero-padding:
% 
\begin{equation}
\bar{h}[k]=\ 
\left\{ \begin{array}{cl} 
w[k]\,h[k]  &    k = 0, \ldots, N-1 \\ 
0     &    k = N, \ldots, M-1
\end{array} 
\right. 
\end{equation} 
% 
which is needed because the range of indices for $h[k]$ and $Q[k]$ do
not match.  $M$ should be at least $2N-1$, but may be chosen to be,
\textit{e.g.}, $2M$ for convenience.
%
The inputs to \texttt{LALStochasticCrossCorrelationStatistic()} are
the (windowed) zero-padded, FFTed data streams
$\widetilde{\bar{h}}_{1}[\ell]$ and $\widetilde{\bar{h}}_{2}[\ell]$,
along with the optimal filter $\widetilde{Q}[\ell]$.  Since the
underlying time series are real, the input series only need to include
the values for $\ell=0,\ldots,P-1$ (where
$P=\left[\frac{M+1}{2}\right]$ is the number of independent elements
in the frequency series) with the elements corresponding to negative
frequencies determined by complex conjugation.  This allows $Y$ to be
computed as
\begin{equation}
Y=\ 
\delta f\ 
\left( 
\widetilde{\bar{h}}_{1}[0]\ 
\widetilde{Q}[0]\ 
\widetilde{\bar{h}}_{2}[0]+\ 
2\sum_{\ell=1}^{P-1}\ 
{\mathrm{Re}} \left\{
\widetilde{\bar{h}}_{1}[\ell]^* \ 
\widetilde{Q}[\ell]\ 
\widetilde{\bar{h}}_{2}[\ell] 
\right\} 
\right)\ .
\label{stochastic:e:shortcut}
\end{equation} 

The routine \texttt{LALStochasticCrossCorrelationStatistic()} is
designed for analyzing non-heterodyned data, so if the input FFTed
datasets have a positive start frequency, and thus represent a range
of frequencies $f_0\le f< f_0 + (P-1)\delta f$, it is assumed that
they were produced by discarding frequencies below $f_0$ from a longer
frequency series, which was still the Fourier transform of a real time
series.  In this case the cross-correlation statistic is calculated
as\footnote{Note that the $P$th frequency bin is not treated
  specially, as would be expected for the Nyquist frequency.  This is
  the appropriate behavior if $M$ is an odd number (so that there is
  no Nyquist bin) or if, as a result of coarse-graining, the Nyquist
  bin has been removed from $\widetilde{Q}$.  At any rate, if there's
  a significant contribution to the cross-correlation statistic from
  the Nyquist frequency, something is wrong.}
\begin{eqnarray}
Y&=&\ 
\delta f\ 
2\sum_{\ell=0}^{P-1}\ 
{\mathrm{Re}} \left\{
\widetilde{\bar{h}}_{1}[\ell]^* \ 
\widetilde{Q}[\ell]\ 
\widetilde{\bar{h}}_{2}[\ell] 
\right\} \nonumber
\\
&\approx&
\int_{-f_0-P\delta f}^{-f_0} df\ 
\widetilde{h}_1(f)^*\ \widetilde{Q}(f)\ \widetilde{h}_2(f)
+ \int_{f_0}^{f_0+P\delta f} df\ 
\widetilde{h}_1(f)^*\ \widetilde{Q}(f)\ \widetilde{h}_2(f)
\label{stochastic:e:bandlimited}
\end{eqnarray}

The frequency sampling parameters (start frequency, frequency spacing,
and number of points) must be the same for both data streams, but if
the optimal filter is more coarsely sampled (for instance, if it
varies in frequency too slowly to warrant the finer resolution), the
data streams will be multiplied in the frequency domain and their
product coarse-grained
(cf.~Sec.~\ref{stochastic:ss:CoarseGrainFrequencySeries.c}) to the
optimal filter resolution before calculating
(\ref{stochastic:e:bandlimited}).

If the \texttt{epochsMatch} boolean variable is set to a true value,
the function will confirm that the start times for both time series
agree.  It can be set to false to allow for cross-correlation of
time-shifted data as a control case.

\subsubsection*{\texttt{LALStochasticHeterodynedCrossCorrelationStatistic()}}

In the case of heterodyned data, one wishes to calculate 
% 
\begin{eqnarray} 
Y
&:=&\int_{t_0}^{t_0+T} dt_1\int_{t_0}^{t_0+T} dt_2\,
w_1(t_1)\, h_1(t_1)^*\, Q(t_1-t_2)\, w_2(t_2)\, h_2(t_2) \nonumber \\
&\approx& \sum_{j=0}^{N-1}\delta t\sum_{k=0}^{N-1}\delta t\,
w_1[k]\, h_1[j]^*\, Q[j-k]\, w_2[k]\, h_2[k] \nonumber \\
&=& \sum_{\ell=0}^{M-1} \delta f\,
\widetilde{\bar{h}}_{1}[\ell]^* \,\widetilde{Q}[\ell]\,
\widetilde{\bar{h}}_{2}[\ell],
\label{stochastic:e:ymaxhet}
\end{eqnarray}
%
In this case, the Fourier transforms of the zero-padded data streams
have $M$ independent elements, which must all be included in the sum,
which is calculated as
\begin{equation}
Y=\ 
\sum_{\ell=0}^{M-1}\ 
\widetilde{\bar{h}}_{1}[\ell]^* \ 
\widetilde{Q}[\ell]\ 
\widetilde{\bar{h}}_{2}[\ell] 
\ .
\label{stochastic:e:heterodyned}
\end{equation} 
While the mean value of the cross-correlation statistic for
heterodyned data should be real (assuming both series were heterodyned
with the same phase), the value for an individual stretch of data will
be complex, so the output is returned as \texttt{COMPLEX8WithUnits}.

\subsubsection*{\texttt{LALStochasticCrossCorrelationSpectrum()}}

For diagnostic purposes, this function calculates the integrand of
(\ref{stochastic:e:ymax}) or (\ref{stochastic:e:ymaxhet}), i.e.
\begin{equation}
  \label{stochastic:e:ccspec}
Y(f)=
\widetilde{\bar{h}}_{1}(f)^* \ 
\widetilde{Q}(f)\ 
\widetilde{\bar{h}}_{2}(f)=
Y[\ell]=
\widetilde{\bar{h}}_{1}[\ell]^* \ 
\widetilde{Q}[\ell]\ 
\widetilde{\bar{h}}_{2}[\ell]  
\end{equation}
and returns it as a frequency series.
 
\subsubsection*{Algorithm}

The function \texttt{LALStochasticCrossCorrelationSpectrum()}
calculates the integrand (\ref{stochastic:e:ccspec}) as follows: First
it calculates $\widetilde{\bar{h}}_{1}[\ell]^* \, 
\widetilde{\bar{h}}_{2}[\ell]$ with
\texttt{LALCCVectorMultiplyConjugate()} and matches the resolution of
the result to that of \verb+input->optimalFilter+ with
\texttt{LALCCoarseGrainFrequencySeries()}.  Then it uses
\texttt{LALCCVectorMultiply()} to calculate
(\ref{stochastic:e:ccspec}) from the input $\widetilde{Q}[\ell]$ and
the coarse-grained $\widetilde{\bar{h}}_{1}[\ell]^* \ 
\widetilde{\bar{h}}_{2}[\ell]$

The functions \texttt{LALStochasticCrossCorrelationStatistic()} and\\
\texttt{LALStochasticHeterodynedCrossCorrelationStatistic()} 
call \\
\texttt{LALStochasticCrossCorrelationSpectrum()} and then
integrate over all frequencies to  calculate
(\ref{stochastic:e:shortcut}) or (\ref{stochastic:e:heterodyned}),
respectively.

\subsubsection*{Uses}

\texttt{LALStochasticCrossCorrelationSpectrum()} calls
\begin{verbatim}
LALCCreateVector()
LALCCVectorMultiplyConjugate()
LALCDestroyVector()
LALCCoarseGrainFrequencySeries()
LALUnitMultiply()
\end{verbatim}

\texttt{LALStochasticCrossCorrelationStatistic()} and 
\texttt{LALStochasticHeterodynedCrossCorrelationStatistic()}
call
\begin{verbatim}
LALCCreateVector()
LALStochasticCrossCorrelationSpectrum()
LALCDestroyVector()
LALUnitMultiply()
\end{verbatim}

\subsubsection*{Notes}
\begin{itemize}
\item When $f_0=0$, $\widetilde{\bar{h}}_{1}[0]$, $\widetilde{Q}[0]$,
  and $\widetilde{\bar{h}}_{2}[0]$ are assumed to be real, but this is
  not checked.
\item The optimal filter $\widetilde{Q}(f)$ is represented by a
  complex frequency series because it will in general be applied to
  whitened data include the different complex whitening filters for
  the two streams.  
  (cf.\ Sec.~\ref{stochastic:ss:StochasticOptimalFilter.c}.)
\item The coarse-graining technique produces the same
  cross-correlation statistic as fine-graining the optimal filter by
  assuming it is zero outside the coarse-grained frequency range and
  constant across each coarse-grained frequency bin.
\item The output units are constructed by combining the input units,
  but under normal circumstances the units will be as follows:
  \begin{eqnarray}
    {} [\widetilde{Q}] &=& \textrm{count}^{-2} \\
    {} [\widetilde{\bar{h}}_{1,2}] &=& \textrm{count}\,\textrm{Hz}^{-1} \\
    {} [Y(f)] &:=& [\widetilde{\bar{h}}_1] 
    [\widetilde{Q}]  [\widetilde{\bar{h}}_2]
    = \textrm{s}^2 \\
    {} [Y] &:=& [Y(f)]\, \textrm{Hz}^{-1} = \textrm{s}
  \end{eqnarray}
\end{itemize}

\vfill{\footnotesize\input{StochasticCrossCorrelationCV}}

******************************************************* </lalLaTeX> */ 
/**************************** <lalLaTeX file="StochasticCrossCorrelationCB">

% \bibitem{stochastic:}

******************************************************* </lalLaTeX> */ 
#include <lal/LALStdlib.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/CoarseGrainFrequencySeries.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>

NRCSID(STOCHASTICCROSSCORRELATIONC, 
"$Id$");

void
LALStochasticCrossCorrelationStatisticCal(
    LALStatus                                *status,
    REAL4WithUnits                           *output,
    const StochasticCrossCorrelationCalInput *input,
    BOOLEAN                                  epochsMatch)
{
  COMPLEX8FrequencySeries ccSpec;
  LALUnitPair unitPair;

  COMPLEX8 *cPtr;
  COMPLEX8 *cStopPtr;

  /* initialize status structure */
  INITSTATUS(status, "LALStochasticCrossCorrelationStatisticCal", \
    STOCHASTICCROSSCORRELATIONC);
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /* output structure */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* input structure */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* optimal filter */
  ASSERT(input->optimalFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* first response function */
  ASSERT(input->responseFunctionOne != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* second response function */
  ASSERT(input->responseFunctionTwo != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* first data stream */
  ASSERT(input->hBarTildeOne != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* second data stream */
  ASSERT(input->hBarTildeTwo != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for optimal filter */
  ASSERT(input->optimalFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member first response function */
  ASSERT(input->responseFunctionOne->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member second response function */
  ASSERT(input->responseFunctionTwo->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for first data stream */
  ASSERT(input->hBarTildeOne->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for second data stream */
  ASSERT(input->hBarTildeTwo->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for optimal filter */
  ASSERT(input->optimalFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member first response function */
  ASSERT(input->responseFunctionOne->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member second response function */
  ASSERT(input->responseFunctionTwo->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for first data stream */
  ASSERT(input->hBarTildeOne->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for second data stream */
  ASSERT(input->hBarTildeTwo->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check for legality of values */

  /* length must be positive */
  ASSERT(input->hBarTildeOne->data->length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  ASSERT(input->optimalFilter->data->length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* start frequency must not be negative */
  if (input->hBarTildeOne->f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }
  if (input->optimalFilter->f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* frequency spacing must be positive */
  ASSERT(input->hBarTildeOne->deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  ASSERT(input->optimalFilter->deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check for mismatches */

  if (input->hBarTildeOne->data->length != input->hBarTildeTwo->data->length) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  if (input->responseFunctionOne->data->length != \
      input->optimalFilter->data->length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  if (input->responseFunctionTwo->data->length != \
      input->optimalFilter->data->length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  

  /* start frequency */
  if (input->hBarTildeOne->f0 != input->hBarTildeTwo->f0) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  if (input->responseFunctionOne->f0 != input->optimalFilter->f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  if (input->responseFunctionTwo->f0 != input->optimalFilter->f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->hBarTildeOne->deltaF != input->hBarTildeTwo->deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  if (input->responseFunctionOne->deltaF != input->optimalFilter->deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  if (input->responseFunctionTwo->deltaF != input->optimalFilter->deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* epoch (start time) */
  if (epochsMatch && \
      ((input->hBarTildeOne->epoch.gpsSeconds != \
        input->hBarTildeTwo->epoch.gpsSeconds) || \
       (input->hBarTildeOne->epoch.gpsNanoSeconds != \
        input->hBarTildeTwo->epoch.gpsNanoSeconds)))
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMTIME, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMTIME);
  }

  ccSpec.data = NULL;

  TRY(LALCCreateVector(status->statusPtr, &(ccSpec.data), \
        input->optimalFilter->data->length), status);

  LALStochasticCrossCorrelationSpectrumCal(status->statusPtr, \
      &ccSpec, input, epochsMatch);

  BEGINFAIL(status)
    TRY(LALCDestroyVector(status->statusPtr, &(ccSpec.data)), status);
  ENDFAIL(status);
 
  if (ccSpec.f0 == 0)
  {
    /* DC contribution */
    output->value = ccSpec.data->data->re / 2.0;

    /* We might want to check that the imaginary parts of the DC
     * components of all the series vanish */

    /* initialize pointers */
    cPtr = ccSpec.data->data + 1;
  } /* if f0 == 0 */
  else
  {
    output->value = 0.0;

    /* initialize pointers */
    cPtr = ccSpec.data->data;
  }
  cStopPtr = ccSpec.data->data + ccSpec.data->length;
  
  /* contributions from positive and (negative) frequency components */
  for (; cPtr < cStopPtr; ++cPtr) 
  {
    output->value += cPtr->re;
  }

  /* normalize */
  output->value *= ccSpec.deltaF * 2.0;

  TRY(LALCDestroyVector(status->statusPtr, &(ccSpec.data)), status);
  
  /* Set output units */
  unitPair.unitOne = &ccSpec.sampleUnits;
  unitPair.unitTwo = &lalHertzUnit;
  TRY(LALUnitMultiply(status->statusPtr, &(output->units), &unitPair), \
      status);
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALStochasticCrossCorrelationStatisticCal() */



void
LALStochasticHeterodynedCrossCorrelationStatisticCal(
    LALStatus                                *status,
    COMPLEX8WithUnits                        *output,
    const StochasticCrossCorrelationCalInput *input,
    BOOLEAN                                  epochsMatch)
{
  COMPLEX8FrequencySeries ccSpec;
  LALUnitPair unitPair;

  COMPLEX8 *cPtr;
  COMPLEX8 *cStopPtr;

  /* initialize status structure */
  INITSTATUS(status, "LALStochasticHeterodynedCrossCorrelationStatisticCal", \
      STOCHASTICCROSSCORRELATIONC);
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /* output structure */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* input structure */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* optimal filter */
  ASSERT(input->optimalFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* first response function */
  ASSERT(input->responseFunctionOne != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* second response function */
  ASSERT(input->responseFunctionTwo != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* first data stream */
  ASSERT(input->hBarTildeOne != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* second data stream */
  ASSERT(input->hBarTildeTwo != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for optimal filter */
  ASSERT(input->optimalFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR,
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member first response function */
  ASSERT(input->responseFunctionOne->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member second response function */
  ASSERT(input->responseFunctionTwo->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for first data stream */
  ASSERT(input->hBarTildeOne->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for second data stream */
  ASSERT(input->hBarTildeTwo->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for optimal filter */
  ASSERT(input->optimalFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member first response function */
  ASSERT(input->responseFunctionOne->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member second response function */
  ASSERT(input->responseFunctionTwo->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for first data stream */
  ASSERT(input->hBarTildeOne->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for second data stream */
  ASSERT(input->hBarTildeTwo->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check for legality of values */

  /* length must be positive */
  ASSERT(input->hBarTildeOne->data->length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  ASSERT(input->optimalFilter->data->length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* start frequency must not be negative */
  if (input->hBarTildeOne->f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }
  if (input->optimalFilter->f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN );
  }

  /* frequency spacing must be positive */
  ASSERT(input->hBarTildeOne->deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  ASSERT(input->optimalFilter->deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check for mismatches */

  if (input->hBarTildeOne->data->length != input->hBarTildeTwo->data->length) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  if (input->responseFunctionOne->data->length != \
      input->optimalFilter->data->length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  if (input->responseFunctionTwo->data->length != \
      input->optimalFilter->data->length) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }  

  /* start frequency */
  if (input->hBarTildeOne->f0 != input->hBarTildeTwo->f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  if (input->responseFunctionOne->f0 != input->optimalFilter->f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  if (input->responseFunctionTwo->f0 != input->optimalFilter->f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->hBarTildeOne->deltaF != input->hBarTildeTwo->deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  if (input->responseFunctionOne->deltaF != input->optimalFilter->deltaF) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  if (input->responseFunctionTwo->deltaF != input->optimalFilter->deltaF) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* epoch (start time) */
  if (epochsMatch && \
      ((input->hBarTildeOne->epoch.gpsSeconds != \
        input->hBarTildeTwo->epoch.gpsSeconds) || \
       (input->hBarTildeOne->epoch.gpsNanoSeconds != \
        input->hBarTildeTwo->epoch.gpsNanoSeconds)))
  {
     ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMTIME, \
         STOCHASTICCROSSCORRELATIONH_MSGEMMTIME );
  }

  ccSpec.data = NULL;

  TRY(LALCCreateVector(status->statusPtr, &(ccSpec.data), \
        input->optimalFilter->data->length), status);

  LALStochasticCrossCorrelationSpectrumCal(status->statusPtr, \
      &ccSpec, input, epochsMatch);

  BEGINFAIL(status)
    TRY(LALCDestroyVector(status->statusPtr, &(ccSpec.data)), status);
  ENDFAIL(status); 

  output->value.re = output->value.im = 0.0;

  /* initialize pointers */
  cPtr = ccSpec.data->data;

  cStopPtr = ccSpec.data->data + ccSpec.data->length;
  
  /* contributions from positive and (negative) frequency components */
  for (; cPtr < cStopPtr; ++cPtr) 
  {
    output->value.re += cPtr->re;
    output->value.im += cPtr->im;
  }
  
  /* normalize */
  output->value.re *= ccSpec.deltaF;
  output->value.im *= ccSpec.deltaF;

  TRY(LALCDestroyVector(status->statusPtr, &(ccSpec.data)), status);
  
  /* Set output units */  
  unitPair.unitOne = &ccSpec.sampleUnits;
  unitPair.unitTwo = &lalHertzUnit;
  TRY(LALUnitMultiply(status->statusPtr, &(output->units), &unitPair), \
      status);
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALStochasticHeterodynedCrossCorrelationStatisticCal() */


void
LALStochasticCrossCorrelationSpectrumCal(
    LALStatus                                *status,
    COMPLEX8FrequencySeries                  *output,
    const StochasticCrossCorrelationCalInput *input,
    BOOLEAN                                  epochsMatch)
{
  LALUnitPair unitPair;
  LALUnit h1H2Units, hc1Hc2Units, r1R2Units, invR1R2Units;

  COMPLEX8FrequencySeries h1StarH2, h1StarH2Coarse, hc1StarHc2Coarse;
  COMPLEX8FrequencySeries r1StarR2;

  FrequencySamplingParams freqParams;
  REAL8 streamF0;
  REAL8 streamDF;
  UINT4 streamLength;
  
  INT4 i;

  /* initialize status structure */
  INITSTATUS(status, "LALStochasticCrossCorrelationSpectrumCal", \
      STOCHASTICCROSSCORRELATIONC);
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /* output series */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for output series */
  ASSERT(output->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for output series */
  ASSERT(output->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* input structure */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR,
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* optimal filter */
  ASSERT(input->optimalFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* first response function */
  ASSERT(input->responseFunctionOne != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* second response function */
  ASSERT(input->responseFunctionTwo != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* first data stream */
  ASSERT(input->hBarTildeOne != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* second data stream */
  ASSERT(input->hBarTildeTwo != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for optimal filter */
  ASSERT(input->optimalFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for first response function */
  ASSERT(input->responseFunctionOne->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for second response function */
  ASSERT(input->responseFunctionTwo->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for first data stream */
  ASSERT(input->hBarTildeOne->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for second data stream */
  ASSERT(input->hBarTildeTwo->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for optimal filter */
  ASSERT(input->optimalFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for first response function */
  ASSERT(input->responseFunctionOne->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for second response function */
  ASSERT(input->responseFunctionTwo->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for first data stream */
  ASSERT(input->hBarTildeOne->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for second data stream */
  ASSERT(input->hBarTildeTwo->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* Check for duplicate pointers */

  /* output series = first response function */
  ASSERT(output != input->responseFunctionOne, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = second response function */
  ASSERT(output != input->responseFunctionTwo, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = first data stream */
  ASSERT(output != input->hBarTildeOne, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = second data stream */
  ASSERT(output != input->hBarTildeTwo, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = first response function */
  ASSERT(output->data != input->responseFunctionOne->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = second response function */
  ASSERT(output->data != input->responseFunctionTwo->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = first data stream */
  ASSERT(output->data != input->hBarTildeOne->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = second data stream */
  ASSERT(output->data != input->hBarTildeTwo->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = first response function */
  ASSERT(output->data->data != input->responseFunctionOne->data->data, \
      status, STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = second response function */
  ASSERT(output->data->data != input->responseFunctionTwo->data->data, \
      status, STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = first data stream */
  ASSERT(output->data->data != input->hBarTildeOne->data->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = second data stream */
  ASSERT(output->data->data != input->hBarTildeTwo->data->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* extract parameters */
  freqParams.length = input->optimalFilter->data->length;
  freqParams.f0 = input->optimalFilter->f0;
  freqParams.deltaF = input->optimalFilter->deltaF;

  /* extract parameters */
  streamLength = input->hBarTildeOne->data->length;
  streamF0 = input->hBarTildeOne->f0;
  streamDF = input->hBarTildeOne->deltaF;

  /* check for legality of values */

  /* length must be positive */
  ASSERT(streamLength != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  ASSERT(freqParams.length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* start frequency must not be negative */
  if (streamF0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }
  if (freqParams.f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* frequency spacing must be positive */
  ASSERT(streamDF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  ASSERT(freqParams.deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check for mismatches */
  
  /* length */
  if (output->data->length != freqParams.length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  if (input->hBarTildeTwo->data->length != streamLength)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  if (input->responseFunctionOne->data->length != freqParams.length) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  if (input->responseFunctionTwo->data->length != freqParams.length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* start frequency */
  if (input->hBarTildeTwo->f0 != streamF0) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  if (input->responseFunctionOne->f0 != freqParams.f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  if (input->responseFunctionTwo->f0 != freqParams.f0) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->hBarTildeTwo->deltaF != streamDF) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  if (input->responseFunctionOne->deltaF != freqParams.deltaF) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  if (input->responseFunctionTwo->deltaF != freqParams.deltaF) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* epoch (start time) */
  if (epochsMatch && \
      ((input->hBarTildeOne->epoch.gpsSeconds != \
        input->hBarTildeTwo->epoch.gpsSeconds) || \
       (input->hBarTildeOne->epoch.gpsNanoSeconds != \
        input->hBarTildeTwo->epoch.gpsNanoSeconds)))
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMTIME, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMTIME );
  }

  h1StarH2 = *(input->hBarTildeOne);
  r1StarR2 = *(input->responseFunctionOne);
  h1StarH2.data = NULL;
  r1StarR2.data = NULL;

  TRY(LALCCreateVector(status->statusPtr, &(h1StarH2.data), streamLength), \
      status);

  TRY(LALCCreateVector(status->statusPtr, &(r1StarR2.data), \
        freqParams.length), status);

  LALCCVectorMultiplyConjugate(status->statusPtr, h1StarH2.data, \
      input->hBarTildeTwo->data, input->hBarTildeOne->data);

  BEGINFAIL(status) 
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2.data)), status);
  ENDFAIL(status);

  LALCCVectorMultiplyConjugate(status->statusPtr, r1StarR2.data, \
      input->responseFunctionTwo->data, input->responseFunctionOne->data);

  BEGINFAIL(status) 
    TRY(LALCDestroyVector(status->statusPtr, &(r1StarR2.data)), status);
  ENDFAIL(status); 

  h1StarH2Coarse.data = NULL;

  LALCCreateVector(status->statusPtr, &(h1StarH2Coarse.data), \
      freqParams.length);

  BEGINFAIL(status) 
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2.data)), status);
  ENDFAIL(status); 

  LALCCoarseGrainFrequencySeries(status->statusPtr, &h1StarH2Coarse, \
      &h1StarH2, &freqParams);

  BEGINFAIL(status) {
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2.data)), status);
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2Coarse.data)), status);
  } ENDFAIL(status); 

  LALCDestroyVector(status->statusPtr, &(h1StarH2.data));

  BEGINFAIL(status) {
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2.data)), status);
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2Coarse.data)), status);
  } ENDFAIL(status); 


  hc1StarHc2Coarse.data = NULL;

  LALCCreateVector(status->statusPtr, &(hc1StarHc2Coarse.data), \
      freqParams.length);

  BEGINFAIL(status) {
    TRY(LALCDestroyVector(status->statusPtr, &(hc1StarHc2Coarse.data)), \
        status);
  } ENDFAIL(status); 

  LALCCVectorDivide(status->statusPtr, hc1StarHc2Coarse.data, \
      h1StarH2Coarse.data, r1StarR2.data);

  BEGINFAIL(status) 
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2Coarse.data)), status);
  ENDFAIL(status); 

  TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2Coarse.data)), status);

  LALSCVectorMultiply(status->statusPtr, output->data, \
      input->optimalFilter->data, hc1StarHc2Coarse.data);

  BEGINFAIL(status) {
    TRY(LALCDestroyVector(status->statusPtr, &(hc1StarHc2Coarse.data)), \
        status);
  } ENDFAIL( status ); 

  TRY(LALCDestroyVector(status->statusPtr, &(hc1StarHc2Coarse.data)), status);

  /* Fill fields of output frequency series */
  output->deltaF = input->optimalFilter->deltaF;
  output->epoch = input->hBarTildeOne->epoch;
  output->f0 = input->optimalFilter->f0;
  strncpy(output->name, "Integrand of cross-correlation statistic", \
      LALNameLength );
  
  /* Set output units */  
  unitPair.unitOne = &(input->hBarTildeOne->sampleUnits);
  unitPair.unitTwo = &(input->hBarTildeTwo->sampleUnits);
  TRY(LALUnitMultiply(status->statusPtr, &h1H2Units, &unitPair), status);
  unitPair.unitOne = &(input->responseFunctionOne->sampleUnits);
  unitPair.unitTwo = &(input->responseFunctionTwo->sampleUnits);
  TRY(LALUnitMultiply(status->statusPtr, &r1R2Units, &unitPair), status);

  invR1R2Units.powerOfTen = -r1R2Units.powerOfTen;
  for (i = 0; i < LALNumUnits; ++i)
  {
    invR1R2Units.unitNumerator[i] = -r1R2Units.unitNumerator[i];
    invR1R2Units.unitDenominatorMinusOne[i] = r1R2Units.unitDenominatorMinusOne[i];
   }
  unitPair.unitOne = &h1H2Units;
  unitPair.unitTwo = &invR1R2Units;
  TRY(LALUnitMultiply(status->statusPtr, &hc1Hc2Units, &unitPair), status);
  unitPair.unitOne = &hc1Hc2Units;
  unitPair.unitTwo = &(input->optimalFilter->sampleUnits);
  TRY(LALUnitMultiply(status->statusPtr, &(output->sampleUnits), &unitPair), \
      status);

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALStochasticCrossCorrelationSpectrumCal() */

/* <lalVerbatim file="StochasticCrossCorrelationCP"> */
void
LALStochasticCrossCorrelationStatistic(
    LALStatus                             *status,
    REAL4WithUnits                        *output,
    const StochasticCrossCorrelationInput *input,
    BOOLEAN                               epochsMatch)
/* </lalVerbatim> */
{
  COMPLEX8FrequencySeries ccSpec;
  LALUnitPair unitPair;

  COMPLEX8 *cPtr;
  COMPLEX8 *cStopPtr;

  /* initialize status structure */
  INITSTATUS(status, "LALStochasticCrossCorrelationStatistic", \
      STOCHASTICCROSSCORRELATIONC);
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /* output structure */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* input structure */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* optimal filter */
  ASSERT(input->optimalFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* first data stream */
  ASSERT(input->hBarTildeOne != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* second data stream */
  ASSERT(input->hBarTildeTwo != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for optimal filter */
  ASSERT(input->optimalFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for first data stream */
  ASSERT(input->hBarTildeOne->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for second data stream */
  ASSERT(input->hBarTildeTwo->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for optimal filter */
  ASSERT(input->optimalFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for first data stream */
  ASSERT(input->hBarTildeOne->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for second data stream */
  ASSERT(input->hBarTildeTwo->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check for legality of values */

  /* length must be positive */
  ASSERT(input->hBarTildeOne->data->length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  ASSERT(input->optimalFilter->data->length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* start frequency must not be negative */
  if (input->hBarTildeOne->f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN );
  }
  if (input->optimalFilter->f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN );
  }

  /* frequency spacing must be positive */
  ASSERT(input->hBarTildeOne->deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  ASSERT(input->optimalFilter->deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check for mismatches */

  if (input->hBarTildeOne->data->length != input->hBarTildeTwo->data->length) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* start frequency */
  if (input->hBarTildeOne->f0 != input->hBarTildeTwo->f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->hBarTildeOne->deltaF != input->hBarTildeTwo->deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* epoch (start time) */
  if (epochsMatch && \
      ((input->hBarTildeOne->epoch.gpsSeconds != \
        input->hBarTildeTwo->epoch.gpsSeconds) || \
       (input->hBarTildeOne->epoch.gpsNanoSeconds != \
        input->hBarTildeTwo->epoch.gpsNanoSeconds)))
  {
     ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMTIME, \
         STOCHASTICCROSSCORRELATIONH_MSGEMMTIME );
  }

  ccSpec.data = NULL;

  TRY(LALCCreateVector(status->statusPtr, &(ccSpec.data), \
        input->optimalFilter->data->length), status);

  LALStochasticCrossCorrelationSpectrum(status->statusPtr, &ccSpec, \
      input, epochsMatch);

  BEGINFAIL(status) 
    TRY(LALCDestroyVector(status->statusPtr, &(ccSpec.data)), status);
  ENDFAIL(status); 
 
  if (ccSpec.f0 == 0) 
  {
    /* DC contribution */
    output->value = ccSpec.data->data->re / 2.0;

    /* We might want to check that the imaginary parts of the DC
       components of all the series vanish */

    /* initialize pointers */
    cPtr = ccSpec.data->data + 1;
  } /* if f0 == 0 */
  else 
  {
    output->value = 0.0;

    /* initialize pointers */
    cPtr = ccSpec.data->data;
  }    
  cStopPtr = ccSpec.data->data + ccSpec.data->length;

  /* contributions from positive and (negative) frequency components */
  for (; cPtr < cStopPtr; ++cPtr) 
  {
    output->value += cPtr->re;
  }
  
  /* normalize */
  output->value *= ccSpec.deltaF * 2.0;

  TRY(LALCDestroyVector(status->statusPtr, &(ccSpec.data)), status);
  
  /* Set output units */  
  unitPair.unitOne = &ccSpec.sampleUnits;
  unitPair.unitTwo = &lalHertzUnit;
  TRY(LALUnitMultiply(status->statusPtr, &(output->units), &unitPair), \
      status);
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALStochasticCrossCorrelationStatistic() */


/* <lalVerbatim file="StochasticCrossCorrelationCP"> */
void
LALStochasticHeterodynedCrossCorrelationStatistic(
    LALStatus                             *status,
    COMPLEX8WithUnits                     *output,
    const StochasticCrossCorrelationInput *input,
    BOOLEAN                               epochsMatch)
/* </lalVerbatim> */
{
  COMPLEX8FrequencySeries ccSpec;
  LALUnitPair unitPair;

  COMPLEX8 *cPtr;
  COMPLEX8 *cStopPtr;

  /* initialize status structure */
  INITSTATUS(status, "LALStochasticHeterodynedCrossCorrelationStatistic", \
      STOCHASTICCROSSCORRELATIONC);
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /* output structure */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* input structure */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* optimal filter */
  ASSERT(input->optimalFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* first data stream */
  ASSERT(input->hBarTildeOne != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* second data stream */
  ASSERT(input->hBarTildeTwo != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for optimal filter */
  ASSERT(input->optimalFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for first data stream */
  ASSERT(input->hBarTildeOne->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for second data stream */
  ASSERT(input->hBarTildeTwo->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for optimal filter */
  ASSERT(input->optimalFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for first data stream */
  ASSERT(input->hBarTildeOne->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for second data stream */
  ASSERT(input->hBarTildeTwo->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check for legality of values */

  /* length must be positive */
  ASSERT(input->hBarTildeOne->data->length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  ASSERT(input->optimalFilter->data->length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* start frequency must be positive */
  if (input->hBarTildeOne->f0 <= 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }
  if (input->optimalFilter->f0 <= 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN );
  }

  /* frequency spacing must be positive */
  ASSERT(input->hBarTildeOne->deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);
  
  ASSERT(input->optimalFilter->f0 > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check for mismatches */

  if (input->hBarTildeOne->data->length != input->hBarTildeTwo->data->length) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* start frequency */
  if (input->hBarTildeOne->f0 != input->hBarTildeTwo->f0) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->hBarTildeOne->deltaF != input->hBarTildeTwo->deltaF) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* epoch (start time) */
  if (epochsMatch && \
      ((input->hBarTildeOne->epoch.gpsSeconds != \
        input->hBarTildeTwo->epoch.gpsSeconds) || \
       (input->hBarTildeOne->epoch.gpsNanoSeconds != \
        input->hBarTildeTwo->epoch.gpsNanoSeconds)))
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMTIME, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMTIME );
  }

  ccSpec.data = NULL;

  TRY(LALCCreateVector(status->statusPtr, &(ccSpec.data), \
        input->optimalFilter->data->length), status);

  LALStochasticCrossCorrelationSpectrum(status->statusPtr, &ccSpec, \
      input, epochsMatch);

  BEGINFAIL(status)
    TRY(LALCDestroyVector(status->statusPtr, &(ccSpec.data)), status);
  ENDFAIL(status); 

  output->value.re = output->value.im = 0.0;

  /* initialize pointers */
  cPtr = ccSpec.data->data;

  cStopPtr = ccSpec.data->data + ccSpec.data->length;
  /* contributions from positive and (negative) frequency components */
  for (; cPtr < cStopPtr; ++cPtr) 
  {
    output->value.re += cPtr->re;
    output->value.im += cPtr->im;
  }
  
  /* normalize */
  output->value.re *= ccSpec.deltaF;
  output->value.im *= ccSpec.deltaF;

  TRY(LALCDestroyVector(status->statusPtr, &(ccSpec.data)), status);
  
  /* Set output units */  
  unitPair.unitOne = &ccSpec.sampleUnits;
  unitPair.unitTwo = &lalHertzUnit;
  TRY(LALUnitMultiply(status->statusPtr, &(output->units), &unitPair), status);
  
  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALStochasticHeterodynedCrossCorrelationStatistic() */


/* <lalVerbatim file="StochasticCrossCorrelationCP"> */
void
LALStochasticCrossCorrelationSpectrum(
    LALStatus                             *status,
    COMPLEX8FrequencySeries               *output,
    const StochasticCrossCorrelationInput *input,
    BOOLEAN                               epochsMatch)
/* </lalVerbatim> */
{
  LALUnitPair unitPair;
  LALUnit h1H2Units;

  COMPLEX8FrequencySeries h1StarH2, h1StarH2Coarse;

  FrequencySamplingParams freqParams;
  REAL8 streamF0;
  REAL8 streamDF;
  UINT4 streamLength;

  /* initialize status structure */
  INITSTATUS(status, "LALStochasticCrossCorrelationSpectrum", \
      STOCHASTICCROSSCORRELATIONC);
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /* output series */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for output series */
  ASSERT(output->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for output series */
  ASSERT(output->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* input structure */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* optimal filter */
  ASSERT(input->optimalFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* first data stream */
  ASSERT(input->hBarTildeOne != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* second data stream */
  ASSERT(input->hBarTildeTwo != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for optimal filter */
  ASSERT(input->optimalFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for first data stream */
  ASSERT(input->hBarTildeOne->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member for second data stream */
  ASSERT(input->hBarTildeTwo->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for optimal filter */
  ASSERT(input->optimalFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for first data stream */
  ASSERT(input->hBarTildeOne->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member for second data stream */
  ASSERT(input->hBarTildeTwo->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* Check for duplicate pointers */

  /* output series = optimal filter */
  ASSERT(output != input->optimalFilter, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = first data stream */
  ASSERT(output != input->hBarTildeOne, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = second data stream */
  ASSERT(output != input->hBarTildeTwo, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = optimal filter */
  ASSERT(output->data != input->optimalFilter->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = first data stream */
  ASSERT(output->data != input->hBarTildeOne->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = second data stream */
  ASSERT(output->data != input->hBarTildeTwo->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = optimal filter */
  ASSERT(output->data->data != input->optimalFilter->data->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = first data stream */
  ASSERT(output->data->data != input->hBarTildeOne->data->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* output series = second data stream */
  ASSERT(output->data->data != input->hBarTildeTwo->data->data, status, \
      STOCHASTICCROSSCORRELATIONH_ESAMEPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGESAMEPTR);

  /* extract parameters */
  freqParams.length = input->optimalFilter->data->length;
  freqParams.f0 = input->optimalFilter->f0;
  freqParams.deltaF = input->optimalFilter->deltaF;

  /* extract parameters */
  streamLength = input->hBarTildeOne->data->length;
  streamF0 = input->hBarTildeOne->f0;
  streamDF = input->hBarTildeOne->deltaF;

  /* check for legality of values */

  /* length must be positive */
  ASSERT(streamLength != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  ASSERT(freqParams.length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* start frequency must not be negative */
  if (streamF0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }
  if (freqParams.f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* frequency spacing must be positive */
  ASSERT(streamDF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  ASSERT(freqParams.deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check for mismatches */
  
  /* length */
  if (output->data->length != freqParams.length) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  if (input->hBarTildeTwo->data->length != streamLength) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* start frequency */
  if (input->hBarTildeTwo->f0 != streamF0) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->hBarTildeTwo->deltaF != streamDF) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* epoch (start time) */
  if (epochsMatch && \
      ((input->hBarTildeOne->epoch.gpsSeconds != \
        input->hBarTildeTwo->epoch.gpsSeconds) || \
       (input->hBarTildeOne->epoch.gpsNanoSeconds != \
        input->hBarTildeTwo->epoch.gpsNanoSeconds)))
  {
     ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMTIME, \
         STOCHASTICCROSSCORRELATIONH_MSGEMMTIME );
  }

  h1StarH2 = *(input->hBarTildeOne);

  h1StarH2.data = NULL;

  TRY(LALCCreateVector(status->statusPtr, &(h1StarH2.data), streamLength), \
      status);

  LALCCVectorMultiplyConjugate(status->statusPtr, h1StarH2.data, \
      input->hBarTildeTwo->data, input->hBarTildeOne->data);

  BEGINFAIL(status) 
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2.data)), status);
  ENDFAIL(status); 

  h1StarH2Coarse.data = NULL;

  LALCCreateVector(status->statusPtr, &(h1StarH2Coarse.data), \
      freqParams.length);

  BEGINFAIL(status) 
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2.data)), status);
  ENDFAIL(status); 

  LALCCoarseGrainFrequencySeries(status->statusPtr, &h1StarH2Coarse, \
      &h1StarH2, &freqParams);

  BEGINFAIL(status) {
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2.data)), status);
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2Coarse.data)), status);
  } ENDFAIL(status); 

  LALCDestroyVector(status->statusPtr, &(h1StarH2.data));

  BEGINFAIL(status) {
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2.data)), status);
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2Coarse.data)), status);
  } ENDFAIL(status); 

  LALCCVectorMultiply(status->statusPtr, output->data, \
      h1StarH2Coarse.data, input->optimalFilter->data);

  BEGINFAIL(status) 
    TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2Coarse.data)), status);
  ENDFAIL(status); 

  TRY(LALCDestroyVector(status->statusPtr, &(h1StarH2Coarse.data)), status);

  /* Fill fields of output frequency series */

  output->deltaF = input->optimalFilter->deltaF;
  output->epoch = input->hBarTildeOne->epoch;
  output->f0 = input->optimalFilter->f0;
  strncpy(output->name, "Integrand of cross-correlation statistic", \
      LALNameLength);
  
  /* Set output units */  
  unitPair.unitOne = &(input->hBarTildeOne->sampleUnits);
  unitPair.unitTwo = &(input->hBarTildeTwo->sampleUnits);
  TRY(LALUnitMultiply(status->statusPtr, &h1H2Units, &unitPair), status);
  unitPair.unitOne = &h1H2Units;
  unitPair.unitTwo = &(input->optimalFilter->sampleUnits);
  TRY(LALUnitMultiply(status->statusPtr, &(output->sampleUnits), \
        &unitPair), status );

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALStochasticCrossCorrelationSpectrum() */
