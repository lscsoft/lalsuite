/************************** <lalVerbatim file="StochasticCrossCorrelationCV">
Author: UTB Relativity Group; contact J. T. Whelan (original by S. Drasco)
$Id$  
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{StochasticCrossCorrelation.c}}
\label{stochastic:ss:StochasticCrossCorrelation.c}

Calculates the value of the standard optimally-filtered 
cross-correlation statistic for stochastic background searches.

\subsubsection*{Prototypes}
\input{StochasticCrossCorrelationCP}
\index{\texttt{LALStochasticCrossCorrelationStatistic()}}

\subsubsection*{Description}
\texttt{LALStochasticCrossCorrelationStatistic()} calculates the value
of the standard optimally-filtered cross-correlation statistic
% 
\begin{eqnarray} 
Y
&:=&\int_{t_0}^{t_0+T} dt_1\int_{t_0}^{t_0+T} dt_2\,
h_1(t_1)\, Q(t_1-t_2)\, h_2(t_2) \nonumber \\
&\approx& \sum_{j=0}^{N-1}\delta t\sum_{k=0}^{N-1}\delta t\,
h_1[j]\, Q[j-k]\, h_2[k] \nonumber \\
&=& \sum_{\ell=-(N-1)}^{N-1} \delta f\,
\widetilde{\bar{h}}_{1}[\ell]^* \,\widetilde{Q}[\ell]\,
\widetilde{\bar{h}}_{2}[\ell],
\label{stochastic:e:ymax}
\end{eqnarray}
%
where the sampling period is $\delta t=T/N$, the frequency spacing is
$\delta f = [(2N-1)\delta t]^{-1}$, the tilde indicates a discrete
Fourier transform normalized to approximate the continuous Fourier
transform:
\begin{equation}
\widetilde{Q}[\ell] := \sum_{k=-(N-1)}^{N-1} \delta t\,
Q[k]\, e^{i2\pi k\ell/(2N-1)}
\end{equation}
the asterisk indicates complex conjugation, and the overbar indicates
zero-padding:
% 
\begin{equation}
\bar{h}[k]=\ 
\left\{ \begin{array}{cl} 
h[k]  &    k = 0, 1, \ldots, N-1 \\ 
0     &    k = -1, -2, \ldots, -(N-1) 
\end{array} 
\right. 
\end{equation} 
% 
which is needed because the range of indices for $h[k]$ and $Q[k]$ do
not match.
%
The inputs to \texttt{LALStochasticCrossCorrelationStatistic()} are
the zero-padded, FFTed data streams $\widetilde{\bar{h}}_{1}[\ell]$
and $\widetilde{\bar{h}}_{2}[\ell]$, along with the optimal filter
$\widetilde{Q}[\ell]$.  Since the underlying time series are real, the
input series only need to include the values for $\ell=0,\ldots,N-1$,
with the elements corresponding to negative indices determined by
complex conjugation.  This allows $Y$ to be computed as
\begin{equation}
Y=\ 
\delta f\ 
\left( 
\widetilde{\bar{h}}_{1}[0]\ 
\widetilde{Q}[0]\ 
\widetilde{\bar{h}}_{2}[0]+\ 
2\sum_{\ell=1}^{N-1}\ 
{\mathrm{Re}} \left\{
\widetilde{\bar{h}}_{1}[\ell]^* \ 
\widetilde{Q}[\ell]\ 
\widetilde{\bar{h}}_{2}[\ell] 
\right\} 
\right)\ .
\label{stochastic:e:shortcut}
\end{equation} 

If the input data streams represent a range of positive frequencies
$f_0\le f \le f_0+(N-1)\delta f$, $f_0>0$ (and thus were produced by
Fourier transforming heterodyned data) we calculate the
cross-correlation statistic
\begin{eqnarray}
Y&=&\ 
\delta f\ 
2\sum_{\ell=1}^{N-1}\ 
{\mathrm{Re}} \left\{
\widetilde{\bar{h}}_{1}[\ell]^* \ 
\widetilde{Q}[\ell]\ 
\widetilde{\bar{h}}_{2}[\ell] 
\right\} \nonumber
\\
&\approx&
\int_{-f_0-N\delta f}^{-f_0} df\ 
\widetilde{h}_1(f)^*\ \widetilde{Q}(f)\ \widetilde{h}_2(f)
+ \int_{f_0}^{f_0+N\delta f} df\ 
\widetilde{h}_1(f)^*\ \widetilde{Q}(f)\ \widetilde{h}_2(f)
\label{stochastic:e:heterodyned}
\end{eqnarray}
 
\subsubsection*{Algorithm}
The value  of $Y$ is calculated via  a straightforward implementation
of (\ref{stochastic:e:shortcut}) or (\ref{stochastic:e:heterodyned}).

\subsubsection*{Uses}
\begin{verbatim}
LALUnitMultiply()
\end{verbatim}
% LALUnitPair
% lalHertzUnit

\subsubsection*{Notes}
\begin{itemize}
\item When $f_0=0$, $\widetilde{\bar{h}}_{1}[0]$, $\widetilde{Q}[0]$,
  and $\widetilde{\bar{h}}_{2}[0]$ are assumed to be real, but this is
  not checked.
\end{itemize}

\vfill{\footnotesize\input{StochasticCrossCorrelationCV}}

******************************************************* </lalLaTeX> */ 
/**************************** <lalLaTeX file="StochasticCrossCorrelationCB">

% \bibitem{stochastic:}

******************************************************* </lalLaTeX> */ 
#include <lal/LALStdlib.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/Units.h>

NRCSID(STOCHASTICCROSSCORRELATIONC, 
       "$Id$");

/* <lalVerbatim file="StochasticCrossCorrelationCP"> */
void
LALStochasticCrossCorrelationStatistic(LALStatus                              *status,
                                       REAL4WithUnits                         *output,
                                       const StochasticCrossCorrelationInput  *input)
/* </lalVerbatim> */
{
  COMPLEX8     *cPtrFilter;
  COMPLEX8     *cPtrH1;
  COMPLEX8     *cPtrH2;
  COMPLEX8     *cStopPtr;

  UINT4         length;
  REAL8         f0;
  REAL8         deltaF;
  LALUnitPair   unitPair1, unitPair2;

  /* initialize status structure */
  INITSTATUS( status, "LALStochasticCrossCorrelation",
              STOCHASTICCROSSCORRELATIONC );
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /*    output structure */
  ASSERT(output != NULL, status,
	 STOCHASTICCROSSCORRELATIONH_ENULLP,
	 STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /*    input structure */
  ASSERT(input != NULL, status,
	 STOCHASTICCROSSCORRELATIONH_ENULLP,
	 STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /*       optimal filter */
  ASSERT(input->optimalFilter != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
	 STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /*       first data stream */
  ASSERT(input->hBarTildeOne != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
	 STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /*       second data stream */
  ASSERT(input->hBarTildeTwo != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
	 STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /*          data member for optimal filter */
  ASSERT(input->optimalFilter->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
	 STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /*          data member for first data stream */
  ASSERT(input->hBarTildeOne->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
	 STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /*          data member for second data stream */
  ASSERT(input->hBarTildeTwo->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
	 STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /*             data-data member for first data stream */
  ASSERT(input->optimalFilter->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
	 STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /*             data-data member for first data stream */
  ASSERT(input->hBarTildeOne->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
	 STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /*             data-data member for second data stream */
  ASSERT(input->hBarTildeTwo->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
	 STOCHASTICCROSSCORRELATIONH_MSGENULLP);



  /* extract parameters (from optimal filter) */
  length = input->optimalFilter->data->length;
  f0     = input->optimalFilter->f0;
  deltaF = input->optimalFilter->deltaF;

  /* check for legality of values */

  /*    length must be positive */
  ASSERT(length != 0, status,
	 STOCHASTICCROSSCORRELATIONH_EZEROLEN,
	 STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /*    start frequency must not be negative */
  if (f0 < 0)
  {
    ABORT( status,
	 STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
	 STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN );
  }

  /*    frequency spacing must be positive */
  ASSERT(deltaF > 0, status, 
         STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
	 STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check for mismatches */
  
  /* length */
  if (input->hBarTildeOne->data->length != length) 
  {
    ABORT(status,
	 STOCHASTICCROSSCORRELATIONH_EMMLEN,
	 STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->hBarTildeTwo->data->length != length) 
  {
    ABORT(status,
	 STOCHASTICCROSSCORRELATIONH_EMMLEN,
	 STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* start frequency */
  if (input->hBarTildeOne->f0 != f0) 
  {
    ABORT(status,
	 STOCHASTICCROSSCORRELATIONH_EMMFMIN,
	 STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }
  if (input->hBarTildeTwo->f0 != f0) 
  {
    ABORT(status,
	 STOCHASTICCROSSCORRELATIONH_EMMFMIN,
	 STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->hBarTildeOne->deltaF != deltaF) 
  {
    ABORT(status,
	 STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
	 STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  if (input->hBarTildeTwo->deltaF != deltaF) 
  {
    ABORT(status,
	 STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
	 STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }

  /* epoch (start time) */
  if ( (input->hBarTildeOne->epoch.gpsSeconds != 
        input->hBarTildeTwo->epoch.gpsSeconds) 
       ||
       (input->hBarTildeOne->epoch.gpsNanoSeconds != 
        input->hBarTildeTwo->epoch.gpsNanoSeconds) )
  {
     ABORT( status,
	 STOCHASTICCROSSCORRELATIONH_EMMTIME,
	 STOCHASTICCROSSCORRELATIONH_MSGEMMTIME );
  }
 
  if (f0 == 0) 
  {
    /* DC contribution */
    output->value = ( input->hBarTildeOne->data->data->re 
                      * input->optimalFilter->data->data->re
                      * input->hBarTildeTwo->data->data->re
                      / 2.0 );

    /* We might want to check that the imaginary parts of the DC
       components of all the series vanish */

    /* initialize pointers */
    cPtrFilter = input->optimalFilter->data->data + 1;
    cPtrH1     = input->hBarTildeOne->data->data + 1;
    cPtrH2     = input->hBarTildeTwo->data->data + 1;
  } /* if f0 == 0 */
  else 
  {
    output->value = 0.0;

    /* initialize pointers */
    cPtrFilter = input->optimalFilter->data->data;
    cPtrH1     = input->hBarTildeOne->data->data;
    cPtrH2     = input->hBarTildeTwo->data->data;
  }    
  cStopPtr    = input->optimalFilter->data->data + length;
  /* contributions from positive and (negative) frequency components */
  for ( ; cPtrFilter < cStopPtr; ++cPtrFilter, ++cPtrH1, ++cPtrH2) 
  {
    output->value += (  cPtrH1->re * (  cPtrFilter->re * cPtrH2->re 
				       - cPtrFilter->im * cPtrH2->im )
                        + cPtrH1->im * (  cPtrFilter->re * cPtrH2->im
					 + cPtrFilter->im * cPtrH2->re ) );
  }
  
  /* normalize */
  output->value *= deltaF * 2.0 ;
  
  /* Set output units */  
  unitPair1.unitOne = input->hBarTildeOne->sampleUnits;
  unitPair1.unitTwo = input->hBarTildeTwo->sampleUnits;
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

} /* LALStochasticCrossCorrelationStatistic() */

