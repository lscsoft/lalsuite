/********************* <lalVerbatim file="StochasticInverseNoiseCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{StochasticInverseNoise.c}}
\label{stochastic:ss:StochasticInverseNoise.c}

Calculates the values of the calibrated and half-calibrated inverse
noise power spectra from the uncalibrated noise power spectrum and the
frequency-domain instrument response function

\subsubsection*{Prototypes}
\idx{LALStochasticInverseNoise()}
\input{StochasticInverseNoiseCP}

\subsubsection*{Description}

As described in Sec.~\ref{stochastic:ss:StochasticOptimalFilter.c},
the most convenient combinations of the noise $P(f)$ (defined by $\langle h(f)h(f')^*\rangle=\delta(f-f')P(f)$) and
instrument response 
$\widetilde{R}(f)=h(f)/h(f)$ to use in
constructing an optimal filter are the inverse half-calibrated power
spectral density
\begin{equation}
  \label{stochastic:e:halfCalibratedPSD}
  \frac{1}{P^{\scriptstyle{\rm HC}}(f)}=\frac{1}{\widetilde{R}(f)
  \,P^{\scriptstyle{\rm C}}(f)}
  =\frac{\widetilde{R}(f)^*}{P(f)}
\end{equation}
and the inverse calibrated PSD
\begin{equation}
  \label{stochastic:e:calibratedPSD}
  \frac{1}{P^{\scriptstyle{\rm C}}(f)}
  =\frac{|\widetilde{R}(f)|^2}{P(f)}
\end{equation}
The function \texttt{LALStochasticInverseNoise()} takes in a
\texttt{REAL4FrequencySeries} describing the uncalibrated PSD
$P(f)$ along with a
\texttt{COMPLEX8FrequencySeries} describing the frequency-domain
response $\widetilde{R}(f)$, and outputs a
\texttt{REAL4FrequencySeries} describing the calibrated inverse PSD
$1/P^{\scriptstyle{\rm C}}(f)$
along with a \texttt{COMPLEX8FrequencySeries} describing the
half-calibrated inverse PSD $1/P^{\scriptstyle{\rm HC}}(f)$.

\subsubsection*{Algorithm}

The output series are filled according to a straightforward
implemementation of
(\ref{stochastic:e:halfCalibratedPSD}-\ref{stochastic:e:calibratedPSD}).
The DC components, if included in the series, are set to zero.

\subsubsection*{Uses}
\begin{verbatim}
LALUnitRaise()
LALUnitMultiply()
strncpy()
\end{verbatim}

\subsubsection*{Notes}
\begin{itemize}
\item Note that although $P^{\scriptstyle{\rm C}}(f)$
 and $P(f)$
  are real, $P^{\scriptstyle{\rm HC}}(f)$ is \emph{complex}.
\item The output units are constructed by combining the input units,
  but under normal circumstances the units will be as follows:
  \begin{eqnarray}
    {} [P] &=& \textrm{count}^{2}\, \textrm{Hz}^{-1}\\
    {} [\widetilde{R}] &=& 10^{18}\,\textrm{strain}^{-1}\,\textrm{count} \\
    {} [1/P^{\scriptstyle{\rm C}}] 
    &:=& [\widetilde{R}]^2 [P]
    = 10^{36}\,\textrm{Hz}\,\textrm{strain}^{-2} \\
    {} [1/P^{\scriptstyle{\rm HC}}] 
    &:=&  [\widetilde{R}] [P]
    = 10^{18}\,\textrm{Hz}\,\textrm{strain}^{-1}\,\textrm{count}^{-1}
  \end{eqnarray}
\end{itemize}

\vfill{\footnotesize\input{StochasticInverseNoiseCV}}

******************************************************* </lalLaTeX> */ 
/********************** <lalLaTeX file="StochasticInverseNoiseCB">

% \bibitem{stochastic:}

******************************************************* </lalLaTeX> */ 

#include <lal/LALStdlib.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/FrequencySeries.h>
#include <string.h>

#define invNoise output->calibratedInverseNoisePSD
#define hwInvNoise output->halfCalibratedInverseNoisePSD
#define wNoise input->unCalibratedNoisePSD
#define wFilter input->responseFunction



NRCSID(STOCHASTICINVERSENOISEC, "$Id$");

void
LALStochasticInverseNoiseCal(
    LALStatus                         *status,
    StochasticInverseNoiseCalOutput   *output,
    const StochasticInverseNoiseInput *input )
{
  REAL8 deltaF;
  REAL8 f0;
  UINT4 length;

  REAL4 *sPtrPW, *sPtrIP, *sStopPtr;
  COMPLEX8 *cPtrR, *cPtrIPHC;
  
  RAT4 power;
  LALUnitPair unitPair;
  LALUnit wInvNoiseUnits;

  COMPLEX8FrequencySeries *hcInvNoise;

  /* initialize status structure */
  INITSTATUS(status, "LALStochasticInverseNoiseCal", STOCHASTICINVERSENOISEC);
  ATTATCHSTATUSPTR(status);

  /*****************************************************************
   *                                                               *
   *                    Test validity of inputs                    *
   *                                                               *
   *****************************************************************/

  /* check that pointer to input structure is not null */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to output structure is not null */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
 
  /* check that pointers to members of input structure are not null */
  ASSERT(wNoise != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointers to members of output structure are not null */
  ASSERT(invNoise != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
  
  
  /* check that pointers to data members of series are not null */
  ASSERT(wNoise->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(invNoise->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
  
  /* check that pointers to data-data members of series are not null */
  ASSERT(wNoise->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(invNoise->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length is not zero */
  length = wNoise->data->length;
  ASSERT(length > 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that lengths of all series match */
  if (wFilter->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (invNoise->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that frequency spacing is positive */
  deltaF = wNoise->deltaF;
  ASSERT(deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check that frequency spacings of input series match */
  if (wFilter->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* set frequency spacing of output series */
  invNoise->deltaF = deltaF;

  /* check that initial frequency is non-negative */
  f0 = wNoise->f0;
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* check that initial frequency of input series match */
  if (wFilter->f0 != f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* set initial frequency of output series */
  invNoise->f0 = f0;

  /* set epochs */
  invNoise->epoch = wNoise->epoch;

  /*---------------Valid data here---------------*/

  strncpy(invNoise->name, "Calibrated invserse noise PSD", LALNameLength);

  /* allocate memory for half calibrated inverse noise */
	TRY(LALCreateCOMPLEX8FrequencySeries(status->statusPtr, &hcInvNoise, \
        "half-calibrated invserse noise PSD", wNoise->epoch, f0, deltaF, \
        lalDimensionlessUnit, length), status);

  /* unit structure manipulation */

  /* Find units of uncalibrated inverse power spectrum */
  power.numerator = -1;
  power.denominatorMinusOne = 0;
  TRY(LALUnitRaise(status->statusPtr, &wInvNoiseUnits, \
        &(wNoise->sampleUnits), &power), status);

  /* multiply by response function units to get half-calibrated inv noise
   * units */
  unitPair.unitOne = &(wFilter->sampleUnits);
  unitPair.unitTwo = &wInvNoiseUnits;
  TRY(LALUnitMultiply(status->statusPtr, &(hcInvNoise->sampleUnits), \
        &unitPair), status);

  /* multiply by response function units to get calibrated inv noise units */
  unitPair.unitTwo = &(hcInvNoise->sampleUnits);
  TRY(LALUnitMultiply(status->statusPtr, &(invNoise->sampleUnits), \
        &unitPair), status);

	
  sStopPtr = wNoise->data->data + length;

  if (f0 == 0)
  {
    /* set DC channel to zero */
    hcInvNoise->data->data[0].re = 0;
    hcInvNoise->data->data[0].im = 0;
    invNoise->data->data[0] = 0;

    /* initialize pointers */
    sPtrPW = wNoise->data->data + 1;
    cPtrR = wFilter->data->data + 1;
    sPtrIP = invNoise->data->data + 1;
    cPtrIPHC = hcInvNoise->data->data + 1;
  } /* if (f0 == 0) */
  else
  {
    /* initialize pointers */
    sPtrPW = wNoise->data->data;
    cPtrR = wFilter->data->data;
    sPtrIP = invNoise->data->data;
    cPtrIPHC = hcInvNoise->data->data;
  }

  for (; sPtrPW < sStopPtr ; ++sPtrPW, ++cPtrR, ++sPtrIP, ++cPtrIPHC)
  {
    *sPtrIP = (cPtrR->re*cPtrR->re + cPtrR->im*cPtrR->im) / *sPtrPW;
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALStochasticInverseNoiseCal() */

/* <lalVerbatim file="StochasticInverseNoiseCP"> */
void
LALStochasticInverseNoise(
    LALStatus                         *status,
    StochasticInverseNoiseOutput      *output,
    const StochasticInverseNoiseInput *input )
/* </lalVerbatim> */
{
  REAL8 deltaF;
  REAL8 f0;
  UINT4 length;

  REAL4 *sPtrPW, *sPtrIP, *sStopPtr;
  COMPLEX8 *cPtrR, *cPtrIPHC;

  RAT4 power;
  LALUnitPair unitPair;
  LALUnit wInvNoiseUnits;

  /* initialize status structure */
  INITSTATUS(status, "LALStochasticInverseNoise", STOCHASTICINVERSENOISEC);
  ATTATCHSTATUSPTR(status);

  /*****************************************************************
   *                                                               *
   *                    Test validity of inputs                    *
   *                                                               *
   *****************************************************************/

  /* check that pointer to input structure is not null */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to output structure is not null */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
 
  /* check that pointers to members of input structure are not null */
  ASSERT(wNoise != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointers to members of output structure are not null */
  ASSERT(invNoise != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(hwInvNoise != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
  
  /* check that pointers to data members of series are not null */
  ASSERT(wNoise->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(invNoise->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(hwInvNoise->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointers to data-data members of series are not null */
  ASSERT(wNoise->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(wFilter->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(invNoise->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  ASSERT(hwInvNoise->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length is not zero */
  length = wNoise->data->length;
  ASSERT(length > 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that lengths of all series match */
  if (wFilter->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (invNoise->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (hwInvNoise->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that frequency spacing is positive */
  deltaF = wNoise->deltaF;
  ASSERT(deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check that frequency spacings of input series match */
  if (wFilter->deltaF != deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* set frequency spacing of output series */
  invNoise->deltaF = deltaF;
  hwInvNoise->deltaF = deltaF;

  /* check that initial frequency is non-negative */
  f0 = wNoise->f0;
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* check that initial frequency of input series match */
  if (wFilter->f0 != f0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* set initial frequency of output series */
  invNoise->f0 = f0;
  hwInvNoise->f0 = f0;

  /* set epochs */
  invNoise->epoch = hwInvNoise->epoch = wNoise->epoch;

  /*---------------Valid data here---------------*/
  
  strncpy(invNoise->name, "Calibrated invserse noise PSD", LALNameLength);
  strncpy(hwInvNoise->name, "half-calibrated invserse noise PSD", \
      LALNameLength);

  /* unit structure manipulation */

  /* Find units of uncalibrated inverse power spectrum */
  power.numerator = -1;
  power.denominatorMinusOne = 0;
  TRY(LALUnitRaise(status->statusPtr, &wInvNoiseUnits, \
        &(wNoise->sampleUnits), &power), status);

  /* multiply by response function units to get half-calibrated inv noise
   * units */
  unitPair.unitOne = &(wFilter->sampleUnits);
  unitPair.unitTwo = &wInvNoiseUnits;
  TRY(LALUnitMultiply(status->statusPtr, &(hwInvNoise->sampleUnits), \
        &unitPair), status);

  /* multiply by response function units to get calibrated inv noise units */
  unitPair.unitTwo = &(hwInvNoise->sampleUnits);
  TRY(LALUnitMultiply(status->statusPtr, &(invNoise->sampleUnits), \
        &unitPair), status);

  sStopPtr = wNoise->data->data + length;

  if (f0 == 0)
  {
    /* set DC channel to zero */
    hwInvNoise->data->data[0].re = 0;
    hwInvNoise->data->data[0].im = 0;
    invNoise->data->data[0] = 0;

    /* initialize pointers */
    sPtrPW = wNoise->data->data + 1;
    cPtrR = wFilter->data->data + 1;
    sPtrIP = invNoise->data->data + 1;
    cPtrIPHC = hwInvNoise->data->data + 1;
  } /* if (f0 == 0) */
  else
  {
    /* initialize pointers */
    sPtrPW = wNoise->data->data;
    cPtrR = wFilter->data->data;
    sPtrIP = invNoise->data->data;
    cPtrIPHC = hwInvNoise->data->data;
  }

  for (; sPtrPW < sStopPtr ; ++sPtrPW, ++cPtrR, ++sPtrIP, ++cPtrIPHC)
  {
    *sPtrIP = (cPtrR->re*cPtrR->re + cPtrR->im*cPtrR->im) / *sPtrPW;
    cPtrIPHC->re = cPtrR->re / *sPtrPW;

    /* minus sign because of complex conjugate */
    cPtrIPHC->im = -cPtrR->im / *sPtrPW; 
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
} /* LALStochasticInverseNoise() */
