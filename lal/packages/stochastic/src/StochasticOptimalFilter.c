/*************************** <lalVerbatim file="StochasticOptimalFilterCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{StochasticOptimalFilter.c}}
\label{stochastic:ss:StochasticOptimalFilter.c}

Calculates the values of the optimal filter function for the
standard cross-correlation statistic.

\subsubsection*{Prototypes}
\idx{LALStochasticOptimalFilter()}
\input{StochasticOptimalFilterCP}

\subsubsection*{Description}

As described in
\cite{stochastic:Allen:1997,stochastic:Allen:1999,stochastic:Finn:2001},
the optimal filter $\widetilde{Q}(f)$ which maximizes the ratio of the
mean $\mu=\langle Y\rangle$ to the standard deviation
$\sigma=\sqrt{\langle (Y-\mu)^2\rangle}$ of the cross-correlation
statistic (\ref{stochastic:e:ymax}) is
%
\begin{equation}
\widetilde{Q}(f)=\lambda\,
\frac{\gamma(f)\,\Omega_{\scriptstyle{\rm GW}}(f)}
{|f|^3\,P_1(f)\,P_2(f)}
\end{equation}
%
where $\lambda$ is a normalization constant, $\gamma(f)$ is the
overlap reduction function (\textit{cf}
Sec.~\ref{stochastic:ss:OverlapReductionFunction.c}) for the two
detectors, $\Omega_{\scriptstyle{\rm GW}}(f)$ is the stochastic
gravitational wave background strength (\textit{cf}
Sec.~\ref{stochastic:ss:OverlapReductionFunction.c}), and $P_i(f)$ is
the power spectral density ($\langle
h_i(f)h_i(f')^*\rangle=\delta(f-f')P_i(f)$)
% (Fourier transform of the auto-correlation function
% $C_{ii}(\tau):=\langle h_i(t+\tau)h_i(t)\rangle$) 
for the $i$th detector.

However, in practice, the data stream coming out of the $i$th detector
is not the strain $h_i(t)=h_{ab}(t,\vec{x}_i)d^{ab}$, but that
convolved with an instrumental response function $R_i(\tau)$ to
produce a ``whitened''\footnote{So named because the instrumental
  response is often tuned to give $h^{\scriptstyle{\rm W}}$ a lower
  dynamic range, for which reason $\widetilde{R}(f)$ is referred to as
  a ``whitening filter''.} data stream
\begin{equation}
h_i^{\scriptstyle{\rm W}}(t) = \int_0^{\infty} d\tau\, R_i(\tau)\,
h_i(t-\tau)
\end{equation}
which has the simpler frequency-domain representation
\begin{equation}
\widetilde{h}_i^{\scriptstyle{\rm W}}(f) 
=  \widetilde{R}_i(f)\, \widetilde{h}_i(f)
\end{equation}
If we want to calculate the cross-correlation statistic $Y$ using the
whitened detector output, the expression is
\begin{equation}
Y
= \int_{-\infty}^{\infty} df\,
\left(
  \frac{\widetilde{\bar{h}}{}^{\scriptstyle{\rm W}}_{1}(f)}
  {\widetilde{R}_{1}(f)}
\right)^*
\,
\widetilde{Q}(f)\,
\left(
  \frac{\widetilde{\bar{h}}{}^{\scriptstyle{\rm W}}_{2}(f)}
  {\widetilde{R}_{2}(f)}
\right)
= \int_{-\infty}^{\infty} df\,
\widetilde{\bar{h}}{}^{\scriptstyle{\rm W}}_{1}(f) ^*
\,
\widetilde{Q}^{\scriptstyle{\rm W}}(f)\,
\widetilde{\bar{h}}{}^{\scriptstyle{\rm W}}_{2}(f)
\end{equation}
where the ``whitened optimal filter'' is
\begin{eqnarray}
\label{stochastic:e:QW}
\widetilde{Q}^{\scriptstyle{\rm W}}(f)
&=&\frac{\widetilde{Q}(f)}{\widetilde{R}_1(f)^*\widetilde{R}_2(f)}
=\lambda\,\left(\frac{1}{\widetilde{R}_1(f)^*P_1(f)}\right)
\frac{\gamma(f)\,\Omega_{\scriptstyle{\rm GW}}(f)}
{|f|^3}\left(\frac{1}{\widetilde{R}_2(f)P_2(f)}\right)
\nonumber
\\
&=&\lambda\,
\frac{\gamma(f)\,\Omega_{\scriptstyle{\rm GW}}(f)}
{|f|^3\,P^{\scriptstyle{\rm HW}}_1(f)^*\,P^{\scriptstyle{\rm HW}}_2(f)}
\ ,
\end{eqnarray}
where $P^{\scriptstyle{\rm HW}}_i(f)=\widetilde{R}_i(f)\,P_i(f)$ is the
``half-whitened'' PSD.  (The whitened PSD is 
$P^{\scriptstyle{\rm W}}_i(f)=|\widetilde{R}_i(f)|^2\,P_i(f)$.)

\texttt{LALStochasticOptimalFilter()} generates a complex frequency
series containing the whitened optimal filter
$\widetilde{Q}^{\scriptstyle{\rm W}}(f)$, taking as inputs real
frequency series representing the overlap reduction function
$\gamma(f)$ and the stochastic gravitational wave background spectrum
${h_{100}}^2\Omega_{\scriptstyle{\rm GW}}(f)$, as well as complex
frequency series representing the half-whitened (inverse) PSDs
$\{1/P^{\scriptstyle{\rm HW}}_i(f)|i=1,2\}$, and as a real parameter
the normalization constant $\lambda$.

\subsubsection*{Algorithm}

The routine \texttt{LALStochasticOptimalFilter()} fills its output
series is filled with the values corresponding to the definition
(\ref{stochastic:e:QW}).

\subsubsection*{Uses}

\begin{verbatim}
LALUnitMultiply()
LALUnitRaise()
LALUnitCompare()
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}
\item If $f_0=0$, the DC element $Q(0)$ is set to zero, regardless of
  the values of the inputs, because the $f^3$ term would make it
  diverge otherwise, and because any conceivable realistic noise
  spectrum will end up blowing up at zero frequency fast enough to
  kill the optimal filter.
\item The implementation of the optimal filter function given here
  assumes a large observation time continuum-limit approximation.  In
  this limit, the Dirichlet kernels (which appear in an exact
  expression for the standard cross-correlation statistic, when
  evaluated in discrete time \cite{stochastic:Finn:2001}; see also
  Sec.~\ref{stochastic:ss:Dirichlet.c}) may be replaced by Dirac delta
  functions.
\item Although $Q(f)$ is real by construction, the whitened optimal
  filter $\widetilde{Q}^{\scriptstyle{\rm W}}(f)$ will in general be
  complex because the whitening filters $\widetilde{R}_i(f)$ for the
  two sites will be different.
\item The expected units for the inputs and output of this function
  are as follows (although the actual output units will be constructed
  from the input units):
  \begin{eqnarray}
    {} [\lambda] &=& 10^{-36}\,\textrm{s}^{-1}\\
    {} [\gamma] &=& \textrm{strain}^{2} \\
    {} [\Omega_{\scriptstyle{\rm GW}}] &=& 1 \\
    {} [1/P^{\scriptstyle{\rm HW}}_{1,2}] 
    &=& 10^{18}\,\textrm{Hz}\,\textrm{strain}^{-1}\,\textrm{count}^{-1} \\
    {} [\widetilde{Q}^{\scriptstyle{\rm W}}] &:=&
    [\lambda] [\gamma][\Omega_{\scriptstyle{\rm GW}}]
    \left[\frac{1}{P^{\scriptstyle{\rm HW}}_1}\right]
    \left[\frac{1}{P^{\scriptstyle{\rm HW}}_2}\right]
    \,\textrm{s}^3
    =
    \textrm{count}^{-2}
  \end{eqnarray}
\end{itemize}

\vfill{\footnotesize\input{StochasticOptimalFilterCV}}

******************************************************* </lalLaTeX> */ 

/**************************** <lalLaTeX file="StochasticOptimalFilterCB">
\bibitem{stochastic:Allen:1997}
  B.~Allen
  ``The stochastic gravity-wave background: sources and detection''
  in \textit{Proceedings of the Les Houches School on Astrophysical Sources of 
  Gravitational Waves}, 
  eds. J.~A.~Marck and J.~P.~Lasota, Cambridge, 373 (1997);
  \href{http://www.arXiv.org/abs/gr-qc/9604033}{gr-qc/9604033}
\bibitem{stochastic:Allen:1999}
  B.~Allen and J.~D.~Romano, ``Detecting a stochastic background of
  gravitational radiation: Signal processing strategies and
  sensitivities''
  Phys.\ Rev.\ D {\bf 59}, 102001 (1999);
  \href{http://www.arXiv.org/abs/gr-qc/9710117}{gr-qc/9710117}
\bibitem{stochastic:Finn:2001}
  L.~S.~Finn and J.~D.~Romano, ``Detecting stochastic gravitational waves:
  Performance of maximum-likelihood and cross-correlation statistics'',
  unpublished.
******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <math.h>
#include <lal/StochasticCrossCorrelation.h>
#include <lal/Units.h>

NRCSID (STOCHASTICOPTIMALFILTERC,
"$Id$");

/* <lalVerbatim file="StochasticOptimalFilterCP"> */
void
LALStochasticOptimalFilter(
            LALStatus                                  *status,
            COMPLEX8FrequencySeries                    *optimalFilter,
            const StochasticOptimalFilterInput         *input,
            const REAL4WithUnits                       *lambda)
/* </lalVerbatim> */
{

  REAL4     gamma;
  REAL4     omega;
  COMPLEX8  p1HWInv;
  COMPLEX8  p2HWInv;

  COMPLEX8  *cPtrOptimalFilter;

  REAL8      f;
  REAL8      f0;
  REAL8      f3;
  REAL8      deltaF;

  /* normalization factor */
  UINT4      i;
  REAL8      realFactor;

  UINT4       length;
  
  RAT4        power;
  LALUnitPair unitPair;
  LALUnit     tmpUnit1, tmpUnit2, checkUnit;

  /* initialize status pointer */
  INITSTATUS(status, "LALStochasticOptimalFilter", STOCHASTICOPTIMALFILTERC);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING ----------------------------------------------------- */

  /***** check for null pointers *****/
  /* input structure */
  ASSERT(input != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* output structure */
  ASSERT(optimalFilter != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
 
  /* overlap member of input */
  ASSERT(input->overlapReductionFunction != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* omega member of input */
  ASSERT(input->omegaGW != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* half-whitened inverse noise 1 of input */
  ASSERT(input->halfWhitenedInverseNoisePSD1 != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* half-whitened inverse noise 2 of input */
  ASSERT(input->halfWhitenedInverseNoisePSD2 != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of overlap */
  ASSERT(input->overlapReductionFunction->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of omega */
  ASSERT(input->omegaGW->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of half-whitened inverse noise 1 */
  ASSERT(input->halfWhitenedInverseNoisePSD1->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of half-whitened inverse noise 2 */
  ASSERT(input->halfWhitenedInverseNoisePSD2->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of output */
  ASSERT(optimalFilter->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of overlap */
  ASSERT(input->overlapReductionFunction->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of omega */
  ASSERT(input->omegaGW->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of half whitened inverse noise 1 */
  ASSERT(input->halfWhitenedInverseNoisePSD1->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of half whitened inverse noise 2 */
  ASSERT(input->halfWhitenedInverseNoisePSD2->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of output structure */
  ASSERT(optimalFilter->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLPTR,
         STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /*** done with null pointers ***/

  /* extract parameters from overlap */
  length = input->overlapReductionFunction->data->length;
  f0     = input->overlapReductionFunction->f0;
  deltaF = input->overlapReductionFunction->deltaF;

  /**** check for legality ****/
  /* length must be positive */
  ASSERT(length != 0, status,
         STOCHASTICCROSSCORRELATIONH_EZEROLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* start frequency must not be negative */
  if (f0 < 0)
  {
    ABORT( status,
         STOCHASTICCROSSCORRELATIONH_ENEGFMIN,
         STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN );
  }

  /* frequency spacing must be positive */
  ASSERT(deltaF > 0, status, 
         STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF,
         STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /** check for mismatches **/
  /* length */
  if (input->omegaGW->data->length != length) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->halfWhitenedInverseNoisePSD1->data->length != length) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->halfWhitenedInverseNoisePSD2->data->length != length) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (optimalFilter->data->length != length) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* initial frequency */
  if (input->omegaGW->f0 != f0) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMFMIN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }
  if (input->halfWhitenedInverseNoisePSD1->f0 != f0) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMFMIN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }
  if (input->halfWhitenedInverseNoisePSD2->f0 != f0) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMFMIN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->omegaGW->deltaF != deltaF) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  if (input->halfWhitenedInverseNoisePSD1->deltaF != deltaF) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  if (input->halfWhitenedInverseNoisePSD2->deltaF != deltaF) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  /* EVERYHTING OKAY HERE! ---------------------------------------------- */

  /* assign parameters to optimalFilter */
  optimalFilter->f0                   = f0;
  optimalFilter->deltaF               = deltaF;
  optimalFilter->epoch.gpsSeconds     = 0;
  optimalFilter->epoch.gpsNanoSeconds = 0;
  strncpy( optimalFilter->name, "Optimal filter for stochastic search",
           LALNameLength );

  /* All the powers we use are integers, so we can do this once here */
  power.denominatorMinusOne = 0;

  /******* Set tmpUnit1 to dims of Omega/H0^2 ******/

  /* First, set it to dims of H0 */

  tmpUnit1 = lalHertzUnit;

  /* Account for scaled units of Hubble constant */
  tmpUnit1.powerOfTen -= 18;

  /* Now set tmpUnit2 to dims of H0^-2 */
  power.numerator    = -2;
  TRY( LALUnitRaise(status->statusPtr, &tmpUnit2, &tmpUnit1, &power)
       , status );

  unitPair.unitOne = &(input->omegaGW->sampleUnits);
  unitPair.unitTwo = &tmpUnit2;
  
  TRY( LALUnitMultiply(status->statusPtr, &tmpUnit1, &unitPair), status );
  /* Now tmpUnit1 has units of Omega/H0^2 */

  /* Now we need to set the Optimal Filter Units equal to the units of */
  /* lambda*gamma*Omega*f^-3*P1HW^-1*P2HW^-1) */

  unitPair.unitOne = &(input->halfWhitenedInverseNoisePSD1->sampleUnits);
  unitPair.unitTwo = &(input->halfWhitenedInverseNoisePSD2->sampleUnits);
  TRY( LALUnitMultiply(status->statusPtr, &tmpUnit1, &unitPair)
       , status );
  /* tmpUnit1 now holds the units of P1HW^-1*P2HW^-1 */

  power.numerator = -3;
  TRY( LALUnitRaise(status->statusPtr, &tmpUnit2, &lalHertzUnit, &power)
       , status );
  /* tmpUnit2 now holds the units of f^-3 */

  unitPair.unitOne = &tmpUnit1;
  unitPair.unitTwo = &tmpUnit2;
  TRY( LALUnitMultiply(status->statusPtr, &checkUnit, &unitPair), status );
  /* checkUnit now holds the units of f^-3*P1HW^-1*P2HW^-1) */

  unitPair.unitOne = &checkUnit;
  unitPair.unitTwo = &(input->omegaGW->sampleUnits);
  TRY( LALUnitMultiply(status->statusPtr, &tmpUnit1, &unitPair), status );
  /* tmpUnit1 now holds units of Omega*f^-3*P1HW^-1*P2HW^-1) */

  unitPair.unitOne = &tmpUnit1;
  unitPair.unitTwo = &(input->overlapReductionFunction->sampleUnits);
  TRY( LALUnitMultiply(status->statusPtr, &tmpUnit2, &unitPair), status );
  /* tmpUnit2 now holds units of gamma*Omega*f^-3*P1HW^-1*P2HW^-1) */

  unitPair.unitOne = &(lambda->units);
  unitPair.unitTwo = &tmpUnit2;
  TRY( LALUnitMultiply(status->statusPtr, &(optimalFilter->sampleUnits),
                       &unitPair)
       , status );

  /* Done with unit manipulation */

  optimalFilter->data->data[0].re = 0;
  optimalFilter->data->data[0].im = 0;
  
  /* calculate optimal filter values */  
  for (i = ( f0 == 0 ? 1 : 0 ) ; i < length; ++i)
  {
    f = f0 + deltaF * (REAL8) i;
    
    f3 = f * f * f;
    
    omega = input->omegaGW->data->data[i];
    gamma = input->overlapReductionFunction->data->data[i];
    p1HWInv = input->halfWhitenedInverseNoisePSD1->data->data[i];
    p2HWInv = input->halfWhitenedInverseNoisePSD2->data->data[i];
    
    cPtrOptimalFilter = &(optimalFilter->data->data[i]);
    
    realFactor = gamma * omega * lambda->value / f3;
    
    cPtrOptimalFilter->re = realFactor * 
      ( (p1HWInv.re) * (p2HWInv.re) + (p1HWInv.im) * (p2HWInv.im) );
    
    cPtrOptimalFilter->im = realFactor * 
      ( (p1HWInv.re) * (p2HWInv.im) - (p1HWInv.im) * (p2HWInv.re) );
  }
  
  DETATCHSTATUSPTR(status);
  RETURN(status);

}
