/************************** <lalVerbatim file="StochasticCrossCorrelationVarianceCV">
Author: UTB Relativity Group; contact whelan@oates.utb.edu
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{StochasticCrossCorrelationVariance.c}}
\label{stochastic:ss:StochasticCrossCorrelationVariance.c}

Calculates the expected variance per unit time of the standard
optimally-filtered cross-correlation statistic for stochastic
background searches.

\subsubsection*{Prototypes}
\idx{LALStochasticCrossCorrelationVariance()}
\input{StochasticCrossCorrelationVarianceCP}

\subsubsection*{Description}
 
\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{StochasticCrossCorrelationVarianceCV}}

******************************************************* </lalLaTeX> */ 
/**************************** <lalLaTeX file="StochasticCrossCorrelationVarianceCB">

% \bibitem{stochastic:}

******************************************************* </lalLaTeX> */ 
#include <lal/LALStdlib.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/CoarseGrainFrequencySeries.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>

NRCSID(STOCHASTICCROSSCORRELATIONVARIANCEC, 
"$Id$");

/* <lalVerbatim file="StochasticCrossCorrelationVarianceCP"> */
void 
LALStochasticCrossCorrelationVariance(
            LALStatus                                      *status,
            REAL4WithUnits                                 *output,
            const StochasticCrossCorrelationVarianceInput  *input,
            BOOLEAN                                         heterodyned)
/* </lalVerbatim> */
{
  LALUnitPair   unitPair1, unitPair2;

  REAL4FrequencySeries      p1P2, p1P2Coarse;
  FrequencySamplingParams   freqParams;

  REAL8         psdF0;
  REAL8         psdDF;
  UINT4         psdLength;

  REAL4      *sPtr;
  COMPLEX8   *cPtr, *cStopPtr;

  /* initialize status structure */
  INITSTATUS( status, "LALStochasticCrossCorrelationVariance",
              STOCHASTICCROSSCORRELATIONVARIANCEC );
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /*    output structure */
  ASSERT(output != NULL, status,
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*    input structure */
  ASSERT(input != NULL, status,
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*       optimal filter */
  ASSERT(input->optimalFilter != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*       first PSD */
  ASSERT(input->noisePSDOne != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*       second PSD */
  ASSERT(input->noisePSDTwo != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*          data member for optimal filter */
  ASSERT(input->optimalFilter->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*          data member for first PSD */
  ASSERT(input->noisePSDOne->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*          data member for second PSD */
  ASSERT(input->noisePSDTwo->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*             data-data member for optimal filter */
  ASSERT(input->optimalFilter->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*             data-data member for first PSD */
  ASSERT(input->noisePSDOne->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*             data-data member for second PSD */
  ASSERT(input->noisePSDTwo->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* extract parameters */
  freqParams.length = input->optimalFilter->data->length;
  freqParams.f0     = input->optimalFilter->f0;
  freqParams.deltaF = input->optimalFilter->deltaF;

  /* extract parameters */
  psdLength = input->noisePSDOne->data->length;
  psdF0     = input->noisePSDOne->f0;
  psdDF = input->noisePSDOne->deltaF;

  /* check for legality of values */

  /*    length must be positive */
  ASSERT(psdLength != 0, status,
         STOCHASTICCROSSCORRELATIONH_EZEROLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  ASSERT(freqParams.length != 0, status,
         STOCHASTICCROSSCORRELATIONH_EZEROLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /*    start frequency must not be negative */
  if (psdF0 < 0)
  {
    ABORT( status,
         STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
         STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN );
  }
  if (freqParams.f0 < 0)
  {
    ABORT( status,
         STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
         STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN );
  }

  /*    frequency spacing must be positive */
  ASSERT(psdDF > 0, status, 
         STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
         STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  ASSERT(freqParams.deltaF > 0, status, 
         STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
         STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check for mismatches */
  
  /* length */
  if (input->noisePSDTwo->data->length != psdLength) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* start frequency */
  if (input->noisePSDTwo->f0 != psdF0) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMFMIN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->noisePSDTwo->deltaF != psdDF) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  p1P2 = *(input->noisePSDOne);

  p1P2.data = NULL;

  TRY(LALSCreateVector(status->statusPtr, &(p1P2.data), psdLength),
      status);

  LALSSVectorMultiply(status->statusPtr, p1P2.data,
		      input->noisePSDTwo->data,
		      input->noisePSDOne->data);

  BEGINFAIL( status ) 
    TRY(LALSDestroyVector(status->statusPtr, &(p1P2.data)), status);
  ENDFAIL( status ); 

  p1P2Coarse.data = NULL;

  LALSCreateVector(status->statusPtr, &(p1P2Coarse.data),
                   freqParams.length);

  BEGINFAIL( status ) 
    TRY(LALSDestroyVector(status->statusPtr, &(p1P2.data)), status);
  ENDFAIL( status ); 

  LALSCoarseGrainFrequencySeries(status->statusPtr, &p1P2Coarse,
                                 &p1P2, &freqParams);

  BEGINFAIL( status ) {
    TRY(LALSDestroyVector(status->statusPtr, &(p1P2.data)), status);
    TRY(LALSDestroyVector(status->statusPtr, &(p1P2Coarse.data)), status);
  } ENDFAIL( status ); 

  if (input->optimalFilter->f0 == 0) 
  {
    /* DC contribution */
    output->value = p1P2Coarse.data->data[0]
      * input->optimalFilter->data->data->re 
      * input->optimalFilter->data->data->re
      / 2.0;

    /* We might want to check that the imaginary parts of the DC
       components of all the series vanish */

    /* initialize pointers */
    sPtr = p1P2Coarse.data->data + 1;
    cPtr = input->optimalFilter->data->data + 1;
  } /* if f0 == 0 */
  else 
  {
    output->value = 0.0;

    /* initialize pointers */
    sPtr = p1P2Coarse.data->data;
    cPtr = input->optimalFilter->data->data;
  }    
  cStopPtr = input->optimalFilter->data->data 
    + input->optimalFilter->data->length;
  /* contributions from positive and (negative) frequency components */
  for ( ; cPtr < cStopPtr; ++cPtr, ++sPtr ) 
  {
    output->value += *sPtr * (cPtr->re * cPtr->re + cPtr->im * cPtr->im);
  }
  
  /* normalize */
  output->value *= p1P2.deltaF * 2.0;

  LALSDestroyVector(status->statusPtr, &(p1P2.data));

  BEGINFAIL( status ) {
    TRY(LALSDestroyVector(status->statusPtr, &(p1P2.data)), status);
    TRY(LALSDestroyVector(status->statusPtr, &(p1P2Coarse.data)), status);
  } ENDFAIL( status ); 

  TRY(LALSDestroyVector(status->statusPtr, &(p1P2Coarse.data)), status);

  /* Set output units */  
  unitPair1.unitOne = input->noisePSDOne->sampleUnits;
  unitPair1.unitTwo = input->noisePSDTwo->sampleUnits;
  TRY( LALUnitMultiply(status->statusPtr, &(unitPair2.unitOne), &unitPair1),
       status );
  unitPair1.unitOne = input->optimalFilter->sampleUnits;
  unitPair1.unitTwo = lalHertzUnit;
  TRY( LALUnitMultiply(status->statusPtr, &(unitPair2.unitTwo), &unitPair1),
       status );
  TRY( LALUnitMultiply(status->statusPtr, &(output->units), &unitPair2),
       status );

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALStochasticCrossCorrelationVariance() */
