/********************* <lalVerbatim file="StochasticInverseNoiseCV">
Author: UTB Relativity Group; contact J. T. Whelan
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{StochasticInverseNoise.c}}
\label{stochastic:ss:StochasticInverseNoise.c}

Calculates the values of the unwhitened and half-whitened inverse
noise power spectra from the whitened noise power spectrum and the
whitening filter (frequency-domain instrument response).

\subsubsection*{Prototypes}
\idx{LALStochasticInverseNoise()}
\input{StochasticInverseNoiseCP}

\subsubsection*{Description}

As described in Sec.~\ref{stochastic:ss:StochasticOptimalFilter.c},
the most convenient combinations of the noise $P^{\scriptstyle{\rm
    W}}(f)$ (defined by $\langle h^{\scriptstyle{\rm
    W}}(f)h^{\scriptstyle{\rm
    W}}(f')^*\rangle=\delta(f-f')P^{\scriptstyle{\rm W}}(f)$) and
instrument response (known in the frequency domain as the whitening
filter) $\widetilde{R}(f)=h^{\scriptstyle{\rm W}}(f)/h(f)$ to use in
constructing an optimal filter are the inverse half-whitened power
spectral density
\begin{equation}
  \label{stochastic:e:halfWhitenedPSD}
  \frac{1}{P^{\scriptstyle{\rm HW}}(f)}=\frac{1}{\widetilde{R}(f)\,P(f)}
  =\frac{\widetilde{R}(f)^*}{P^{\scriptstyle{\rm W}}(f)}
\end{equation}
and the inverse unwhitened PSD
\begin{equation}
  \label{stochastic:e:unWhitenedPSD}
  \frac{1}{P(f)}
  =\frac{|\widetilde{R}(f)|^2}{P^{\scriptstyle{\rm W}}(f)}
\end{equation}
The function \texttt{LALStochasticInverseNoise()} takes in a
\texttt{REAL4FrequencySeries} describing the whitened PSD
$P^{\scriptstyle{\rm W}}(f)$ along with a
\texttt{COMPLEX8FrequencySeries} describing the frequency-domain
response $\widetilde{R}(f)$, and outputs a
\texttt{REAL4FrequencySeries} describing the unwhitened inverse PSD
$1/P(f)$ along with a \texttt{COMPLEX8FrequencySeries} describing the
half-whitened inverse PSD $1/P^{\scriptstyle{\rm HW}}(f)$.

\subsubsection*{Algorithm}

The output series are filled according to a straightforward
implemementation of
(\ref{stochastic:e:halfWhitenedPSD}-\ref{stochastic:e:unWhitenedPSD}).
The DC components, if included in the series, are set to zero.

\subsubsection*{Uses}
\begin{verbatim}
LALUnitRaise()
LALUnitMultiply()
strncpy()
\end{verbatim}

\subsubsection*{Notes}
\begin{itemize}
\item Note that although $P(f)$ and $P^{\scriptstyle{\rm W}}(f)$
  are real, $P^{\scriptstyle{\rm HW}}(f)$ is \emph{complex}.
\item The output units are constructed by combining the input units,
  but under normal circumstances the units will be as follows:
  \begin{eqnarray}
    {} [P^{\scriptstyle{\rm W}}] &=& \textrm{count}^{2}\, \textrm{Hz}^{-1}\\
    {} [\widetilde{R}] &=& 10^{18}\,\textrm{strain}^{-1}\,\textrm{count} \\
    {} [1/P] &:=& [\widetilde{R}]^2 [P^{\scriptstyle{\rm W}}]
    = 10^{36}\,\textrm{Hz}\,\textrm{strain}^{-2} \\
    {} [1/P^{\scriptstyle{\rm HW}}] 
    &:=&  [\widetilde{R}] [P^{\scriptstyle{\rm W}}]
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
#include <string.h>

NRCSID(STOCHASTICINVERSENOISEC, "$Id$");

/* <lalVerbatim file="StochasticInverseNoiseCP"> */
void
LALStochasticInverseNoise( LALStatus                          *status,
                           StochasticInverseNoiseOutput       *output,
                           const StochasticInverseNoiseInput  *input )
/* </lalVerbatim> */
{
  REAL8          deltaF;
  REAL8          f0;
  UINT4          length;

  REAL4FrequencySeries          *invNoise;
  COMPLEX8FrequencySeries       *hwInvNoise;
  REAL4FrequencySeries          *wNoise;
  COMPLEX8FrequencySeries       *wFilter;

  REAL4          *sPtrPW, *sPtrIP, *sStopPtr;
  COMPLEX8       *cPtrR, *cPtrIPHW;

  RAT4        power;
  LALUnitPair unitPair;

  /* initialize status structure */
  INITSTATUS( status, "LALStochasticInverseNoise", STOCHASTICINVERSENOISEC );
  ATTATCHSTATUSPTR (status);

  /*****************************************************************
   *                                                               *
   *                    Test validity of inputs                    *
   *                                                               *
   *****************************************************************/



  /* check that pointer to input structure is not null */
  ASSERT( input != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );

  /* check that pointer to output structure is not null */
  ASSERT( output != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );
 
  /* check that pointers to members of input structure are not null */
  wNoise      = input->whitenedNoisePSD;
  ASSERT( wNoise != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );
  wFilter     = input->whiteningFilter;
  ASSERT( wFilter != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );

  /* check that pointers to members of output structure are not null */
  invNoise    = output->unWhitenedInverseNoisePSD;
  ASSERT( invNoise != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );
  hwInvNoise  = output->halfWhitenedInverseNoisePSD;
  ASSERT( hwInvNoise != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );
  
  /* check that pointers to data members of series are not null */
  ASSERT( wNoise->data != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, 
          STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );
  ASSERT( wFilter->data != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, 
          STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );
  ASSERT( invNoise->data != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, 
          STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );
  ASSERT( hwInvNoise->data != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, 
          STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );

  /* check that pointers to data-data members of series are not null */
  ASSERT( wNoise->data->data != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, 
          STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );
  ASSERT( wFilter->data->data != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, 
          STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );  
  ASSERT( invNoise->data->data != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, 
          STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );
  ASSERT( hwInvNoise->data->data != NULL, status, STOCHASTICCROSSCORRELATIONH_ENULLPTR, 
          STOCHASTICCROSSCORRELATIONH_MSGENULLPTR );

  /* check that length is not zero */
  length = wNoise->data->length;
  ASSERT(length > 0, status, STOCHASTICCROSSCORRELATIONH_EZEROLEN, STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that lengths of all series match */
  if ( wFilter->data->length != length ) {
    ABORT( status, STOCHASTICCROSSCORRELATIONH_EMMLEN, STOCHASTICCROSSCORRELATIONH_MSGEMMLEN );
  }
  if ( invNoise->data->length != length ) {
    ABORT( status, STOCHASTICCROSSCORRELATIONH_EMMLEN, STOCHASTICCROSSCORRELATIONH_MSGEMMLEN );
  }
  if ( hwInvNoise->data->length != length ) {
    ABORT( status, STOCHASTICCROSSCORRELATIONH_EMMLEN, STOCHASTICCROSSCORRELATIONH_MSGEMMLEN );
  }

  /* check that frequency spacing is positive */
  deltaF = wNoise->deltaF;
  ASSERT(deltaF > 0, status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, 
         STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check that frequency spacings of input series match */
  if ( wFilter->deltaF != deltaF ) {
    ABORT( status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF );
  }
  /* set frequency spacing of output series */
  invNoise->deltaF = deltaF;
  hwInvNoise->deltaF = deltaF;

  /* check that initial frequency is non-negative */
  f0 = wNoise->f0;
   if ( f0 < 0 ) {
     ABORT( status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
            STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN );
   }

  /* check that initial frequency of input series match */
  if ( wFilter->f0 != f0 ) {
    ABORT( status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN );
  }

  /* set initial frequency of output series */
  invNoise->f0 = f0;
  hwInvNoise->f0 = f0;

  /* set epochs */
  invNoise->epoch = hwInvNoise->epoch = wNoise->epoch;

  /*---------------Valid data here---------------*/
  
  strncpy(invNoise->name,"Unwhitened invserse noise PSD",LALNameLength);
  strncpy(hwInvNoise->name,"Half-whitened invserse noise PSD",LALNameLength);

  /* unit structure manipulation */
  power.numerator = 2;
  power.denominatorMinusOne = 0;
  TRY(LALUnitRaise(status->statusPtr, &(unitPair.unitOne),
                   &(wFilter->sampleUnits),
                   &power),status);
  power.numerator = -1;
  power.denominatorMinusOne = 0;
  TRY(LALUnitRaise(status->statusPtr, &(unitPair.unitTwo),
                   &(wNoise->sampleUnits), 
                   &power),status);
  TRY(LALUnitMultiply(status->statusPtr,&(invNoise->sampleUnits),
                      &unitPair),status);
  unitPair.unitOne = wFilter->sampleUnits;
  TRY(LALUnitMultiply(status->statusPtr, &(hwInvNoise->sampleUnits), 
                      &unitPair),status);

  sStopPtr = wNoise->data->data + length;

  if (f0 == 0)
  {
      /* set DC channel to zero */
      hwInvNoise->data->data[0].re = hwInvNoise->data->data[0].im =
      invNoise->data->data[0] = 0;

      /* initialize pointers */
      sPtrPW = wNoise->data->data + 1;
      cPtrR = wFilter->data->data + 1;
      sPtrIP = invNoise->data->data + 1;
      cPtrIPHW = hwInvNoise->data->data + 1;
  } /* if (f0 == 0) */
  else
  {
    /* initialize pointers */
    sPtrPW = wNoise->data->data;
    cPtrR = wFilter->data->data;
    sPtrIP = invNoise->data->data;
    cPtrIPHW = hwInvNoise->data->data;
  }


  for ( sPtrPW, cPtrR, sPtrIP, cPtrIPHW;
        sPtrPW < sStopPtr ;
        ++sPtrPW, ++cPtrR, ++sPtrIP, ++cPtrIPHW )
  {
    *sPtrIP = ( cPtrR->re*cPtrR->re + cPtrR->im*cPtrR->im ) / *sPtrPW;
    cPtrIPHW->re = cPtrR->re / *sPtrPW;
    /* minus sign because of complex conjugate */
    cPtrIPHW->im = - cPtrR->im / *sPtrPW; 
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);
};




