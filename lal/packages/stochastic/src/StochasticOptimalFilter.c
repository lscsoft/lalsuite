/*************************** <lalVerbatim file="StochasticOptimalFilterCV">
Author: UTB Relativity Group; contact J. T. Whelan
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{StochasticOptimalFilter.c}}
\label{stochastic:ss:StochasticOptimalFilter.c}

Calculates the values of the optimal filter function for the
standard cross-correlation statistic.

\subsubsection*{Prototypes}
\input{StochasticOptimalFilterCP}
\index{\texttt{LALStochasticOptimalFilter()}}

\subsubsection*{Description}

As described in
\cite{stochastic:Allen:1997,stochastic:Allen:1999,stochastic:Finn:2001},
the optimal filter $\widetilde{Q}(f)$ which maximizes the ratio of the
mean $\mu=\langle Y\rangle$ to the variance $\sigma=\sqrt{\langle
  (Y-\mu)^2\rangle}$ of the cross-correlation statistic
(\ref{stochastic:e:ymax}) is
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

The normalization constant $\lambda$ is chosen so that the expected mean
value of the cross-correlation statistic is \cite{stochastic:Allen:1999}
\begin{equation}
\mu = \frac{3 {H_0}^2}{20\pi^2}\, T \int_{-\infty}^{\infty} df\, |f|^{-3}\,
\gamma(f)\,\Omega_{\scriptstyle{\rm GW}}(f) 
\widetilde{Q}(f) = \Omega_{\scriptstyle{\rm R}} T
\label{stochastic:e:mu}
\end{equation}
where $T$ is the integration time
(\textit{cf.}~(\ref{stochastic:e:ymax})), and
$\Omega_{\scriptstyle{\rm R}} =\Omega_{\scriptstyle{\rm
    GW}}(f_{\scriptstyle{\rm R}})$ is the overall strength of the
stochastic background (see
Sec.~\ref{stochastic:ss:OverlapReductionFunction.c}).  This sets the
value at
\begin{equation}
\label{stochastic:e:lambda}
\lambda = \frac{20\pi^2\, \Omega_{\scriptstyle{\rm R}}}{3\,{H_0}^2}
\left(
      \int_{-\infty}^\infty \frac{df}{f^6} 
      \frac{[\gamma(f)\,\Omega_{\scriptstyle{\rm GW}}(f)]^2}{P_1(f)P_2(f)}
\right)^{-1}
\end{equation}
Note that since the mean value of the cross-correlation statistic
should be the same whether it's calculated from the whitened or % '
unwhitened data, the normalization constant $\lambda$ involves only
the unwhitened PSDs $\{P_i(f)|i=1,2\}$.

\texttt{LALStochasticOptimalFilter()} generates a complex frequency
series containing the whitened optimal filter
$\widetilde{Q}^{\scriptstyle{\rm W}}(f)$, taking as inputs real
frequency series representing the overlap reduction function
$\gamma(f)$, the stochastic gravitational wave background spectrum
${h_{100}}^2\Omega_{\scriptstyle{\rm GW}}(f)$, and the unwhitened
(inverse) PSDs $\{1/P_i(f)|i=1,2\}$ of the two detectors, as well as
complex frequency series representing the half-whitened (inverse) PSDs
$\{1/P^{\scriptstyle{\rm HW}}_i(f)|i=1,2\}$.

\subsubsection*{Algorithm}

The routine \texttt{LALStochasticOptimalFilter()} first uses
(\ref{stochastic:e:lambda}) to find the normalization constant
$\lambda$ (the amplitude ${h_{100}}^2\Omega_{\scriptstyle{\rm R}}$ is
found from the specified reference frequency $f_{\scriptstyle{\rm R}}$
and the input series representing ${h_{100}}^2\Omega_{\scriptstyle{\rm
    GW}}(f)$ by logarithmic interpolation), then the output series is
filled with the values corresponding to the definition
(\ref{stochastic:e:QW}).

\subsubsection*{Uses}

\begin{verbatim}
LALUnitMultiply()
LALUnitRaise()
LALUnitCompare()
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}
\item This routine does not yet support non-zero heterodyning frequencies.
\item The implementation of the optimal filter function given here
  assumes a large observation time continuum-limit approximation.  In
  this limit, the Dirichlet kernels (which appear in an exact
  expression for the standard cross-correlation statistic, when
  evaluated in discrete time \cite{stochastic:Finn:2001}; see also
  Sec.~\ref{stochastic:ss:Dirichlet.c}) may be replaced by Dirac delta
  functions.
\item The units of the input series are checked for consistency; since
  \cite{stochastic:Allen:1999}
  \begin{equation}
    \langle\widetilde{h}_1(f)^*\widetilde{h}_2(f')\rangle
    = \frac{3H_0^2}{20\pi^2}\delta(f-f')
    \gamma(|f|)\Omega_{\scriptstyle{\rm GW}}(|f|)
  \end{equation}
  and
  \begin{equation}
    \langle\widetilde{h}_i(f)^*\widetilde{h}_i(f')\rangle
    = \delta(f-f')P_i(f)
    \,
  \end{equation}
  we demand that, up to powers of ten,
  \begin{equation}
    ([\gamma][\Omega_{\scriptstyle{\rm GW}}])^2
    =([\widetilde{h}_1][\widetilde{h}_2][T]^{-2})^2
    =[P_1][P_2][T]^{-2}
  \end{equation}
\item This routine, like all those in the \texttt{stochastic} package,
  uses with single precision arithmetic, consistent with the accuracy
  of the continuous approximation for the intended number of data
  points.  However, the limited dynamic range can pose a problem,
  especially in this routine, which uses numbers like $H_0^2$ which
  are far from unity when expressed in the standard SI units.  To
  avoid this problem, the application of (\ref{stochastic:e:lambda})
  takes the units $H_0$ of to be $10^{-18}\,\textrm{s}^{-1}$ rather
  than $\textrm{s}^{-1}$.  In these units $h_{100}H_0$ has a numerical
  value of $3.2407792903$, and inclusion of $(h_{100}H_0)^2$ in an
  expression doesn't risk single-precision overflow.  When %'
  constructing the unit structure of its output,
  \texttt{LALStochasticOptimalFilter()} uses the power-of-ten feature
  of the \texttt{LALUnit} structure to account for the units of $H_0$.
\item Although $Q(f)$ is real by construction, the whitened optimal
  filter $\widetilde{Q}^{\scriptstyle{\rm W}}(f)$ will in general be
  complex because the whitening filters $\widetilde{R}_i(f)$ for the
  two sites will be different.
\item The expected units for the inputs and output of this function
  are as follows (although the actual output units will be constructed
  from the input units):
  \begin{equation}
    ([\gamma][\Omega_{\scriptstyle{\rm GW}}])^2
    =([\widetilde{h}_1][\widetilde{h}_2][T]^{-2})^2
    =[P_1][P_2][T]^{-2}
  \end{equation}
  \begin{eqnarray}
    {} [\gamma] &=& \textrm{strain}^{2} \\
    {} [\Omega_{\scriptstyle{\rm GW}}] &=& 1 \\
    {} [1/P_{1,2}] &=& 10^{36}\,\textrm{Hz}\,\textrm{strain}^{-2} \\
    {} [1/P^{\scriptstyle{\rm HW}}_{1,2}] 
    &=& 10^{18}\,\textrm{Hz}\,\textrm{strain}^{-1}\,\textrm{count}^{-1} \\
    {} [\widetilde{Q}^{\scriptstyle{\rm W}}] &:=&
    10^{36} \left[\frac{1/P^{\scriptstyle{\rm HW}}_1}{1/P_1}\right]
    \left[\frac{1/P^{\scriptstyle{\rm HW}}_2}{1/P_2}\right]
    [\gamma]^{-1}
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
LALStochasticOptimalFilter( LALStatus                          *status,
                            COMPLEX8FrequencySeries            *optimalFilter,
                            const StochasticOptimalFilterInput *input,
                            const REAL8                        fRef )
/* </lalVerbatim> */
{

  REAL4     gamma;
  REAL4     omega;
  REAL4     omegaTimesGamma;
  REAL4     p1Inv;
  REAL4     p2Inv;
  COMPLEX8  p1HWInv;
  COMPLEX8  p2HWInv;

  COMPLEX8  *cPtrOptimalFilter;

  REAL8      f;
  REAL8      f0;
  REAL8      f3;
  REAL8      deltaF;

  /* normalization factor */
  UINT4      i;
  UINT4       xRef;
  REAL4      omega1;
  REAL4      omega2;
  REAL8      freq1;
  REAL8      freq2;
  REAL4      exponent;
  REAL4      omegaRef;   
  REAL8      realFactor;
  REAL8      constant;     /* constants used to calculate lambda */
  REAL8      lambdaInv;
  REAL8      f6; 

  INT4       length;
  
  RAT4        power;
  LALUnitPair unitPair1, unitPair2, unitPair3;
  LALUnit     unit;
  BOOLEAN     result;

  /* initialize status pointer */
  INITSTATUS(status, "LALStochasticOptimalFilter", STOCHASTICOPTIMALFILTERC);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING ----------------------------------------------------- */

  /***** check for null pointers *****/
  /* input structure */
  ASSERT(input != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* output structure */
  ASSERT(optimalFilter != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);
 
  /* overlap member of input */
  ASSERT(input->overlapReductionFunction != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* omega member of input */
  ASSERT(input->omegaGW != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* inverse noise 1 of input */
  ASSERT(input->unWhitenedInverseNoisePSD1 != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* inverse noise 2 of input */
  ASSERT(input->unWhitenedInverseNoisePSD2 != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);
  /* half-whitened inverse noise 1 of input */
  ASSERT(input->halfWhitenedInverseNoisePSD1 != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* half-whitened inverse noise 2 of input */
  ASSERT(input->halfWhitenedInverseNoisePSD2 != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data member of overlap */
  ASSERT(input->overlapReductionFunction->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data member of omega */
  ASSERT(input->omegaGW->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data member of inverse noise 1 */
  ASSERT(input->unWhitenedInverseNoisePSD1->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data member of inverse noise 2 */
  ASSERT(input->unWhitenedInverseNoisePSD2->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data member of half-whitened inverse noise 1 */
  ASSERT(input->halfWhitenedInverseNoisePSD1->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data member of half-whitened inverse noise 2 */
  ASSERT(input->halfWhitenedInverseNoisePSD2->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data member of output */
  ASSERT(optimalFilter->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data-data member of overlap */
  ASSERT(input->overlapReductionFunction->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data-data member of omega */
  ASSERT(input->omegaGW->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data-data member of inverse noise 1 */
  ASSERT(input->unWhitenedInverseNoisePSD1->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data-data member of inverse noise 2 */
  ASSERT(input->unWhitenedInverseNoisePSD2->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data-data member of half whitened inverse noise 1 */
  ASSERT(input->halfWhitenedInverseNoisePSD1->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data-data member of half whitened inverse noise 2 */
  ASSERT(input->halfWhitenedInverseNoisePSD2->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

  /* data-data member of output structure */
  ASSERT(optimalFilter->data->data != NULL, status, 
         STOCHASTICCROSSCORRELATIONH_ENULLP,
         STOCHASTICCROSSCORRELATIONH_MSGENULLP);

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
  if (input->unWhitenedInverseNoisePSD1->data->length != length) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMLEN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->unWhitenedInverseNoisePSD2->data->length != length) 
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
  if (input->unWhitenedInverseNoisePSD1->f0 != f0) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMFMIN,
         STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }
  if (input->unWhitenedInverseNoisePSD2->f0 != f0) 
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
  if (input->unWhitenedInverseNoisePSD1->deltaF != deltaF) 
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EMMDELTAF,
         STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  if (input->unWhitenedInverseNoisePSD2->deltaF != deltaF) 
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
  /** check for reference frequency lower and upper limits **/ 
  if ( fRef < deltaF )
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EOORFREF,
         STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  }
  if ( fRef > ((length-1)*deltaF) )
  {
    ABORT(status,
         STOCHASTICCROSSCORRELATIONH_EOORFREF,
         STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  }
  /* EVERYHTING OKAY HERE! ---------------------------------------------- */

  /* assign parameters to optimalFilter */
  optimalFilter->f0                   = f0;
  optimalFilter->deltaF               = deltaF;
  optimalFilter->epoch.gpsSeconds     = 0;
  optimalFilter->epoch.gpsNanoSeconds = 0;
  strncpy( optimalFilter->name, "Optimal filter for stochastic search",
           LALNameLength );

  /* check that units of gamma, Omega, P1 and P2 are consistent */
  /* we must have (gamma*Omega)^2 with the same units as f^2*P1*P2 */
  /* up to a power of ten */
  unitPair1.unitOne = input->overlapReductionFunction->sampleUnits;
  unitPair1.unitTwo = input->omegaGW->sampleUnits;
  TRY( LALUnitMultiply(status->statusPtr, &unit, &unitPair1)
      , status );
  power.numerator = -2;
  power.denominatorMinusOne = 0;
  TRY( LALUnitRaise(status->statusPtr, &(unitPair3.unitOne), &unit,
                    &power)
       , status );
  power.numerator = 2;
  TRY( LALUnitRaise(status->statusPtr, &(unitPair2.unitOne), &lalHertzUnit,
                    &power)
       , status );
  /* unitPair2.unitOne holds Hz^2 */
  unitPair1.unitOne = input->unWhitenedInverseNoisePSD1->sampleUnits;
  unitPair1.unitTwo = input->unWhitenedInverseNoisePSD2->sampleUnits;
  TRY( LALUnitMultiply(status->statusPtr, &unit, &unitPair1)
      , status );
  power.numerator   = -1;
  TRY( LALUnitRaise(status->statusPtr, &(unitPair2.unitTwo), &unit,
                    &power)
       , status );
  TRY( LALUnitMultiply(status->statusPtr, &(unitPair3.unitTwo), &unitPair2)
      , status );
  TRY( LALUnitMultiply(status->statusPtr, &unit, &unitPair3), status );

  for (i=0; i<LALNumUnits; ++i)
  {
    if ( unit.unitNumerator[i] || unit.unitDenominatorMinusOne[i] )
      {
        ABORT(status,
              STOCHASTICCROSSCORRELATIONH_EWRONGUNITS,
              STOCHASTICCROSSCORRELATIONH_MSGEWRONGUNITS);    
      }
  }

  /* assign correct units to normalized optimal filter */
  power.numerator = -1;
  power.denominatorMinusOne = 0;

  unitPair1.unitOne = input->unWhitenedInverseNoisePSD1->sampleUnits;
  TRY( LALUnitRaise(status->statusPtr, &(unitPair1.unitTwo),
                    &(input->halfWhitenedInverseNoisePSD1->sampleUnits),
                    &power)
       , status );
  TRY( LALUnitMultiply(status->statusPtr, &(unitPair2.unitOne), &unitPair1)
      , status );
  /* unitPair2.unitOne now holds the units of whitening filter #1 */
  unitPair1.unitOne = input->unWhitenedInverseNoisePSD1->sampleUnits;
  TRY( LALUnitRaise(status->statusPtr, &(unitPair1.unitTwo),
                    &(input->halfWhitenedInverseNoisePSD1->sampleUnits),
                    &power)
       , status );
  TRY( LALUnitMultiply(status->statusPtr, &(unitPair2.unitTwo), &unitPair1)
      , status );
  /* unitPair2.unitTwo now holds the units of whitening filter #2 */
  unitPair3.unitOne = input->overlapReductionFunction->sampleUnits;
  TRY( LALUnitMultiply(status->statusPtr, &(unitPair3.unitTwo), &unitPair2)
      , status );
  TRY( LALUnitMultiply(status->statusPtr, &unit, &unitPair3), status );
  TRY( LALUnitRaise(status->statusPtr, &(optimalFilter->sampleUnits),
                    &unit, &power)
       , status );

  /* Now change power of ten to account for scaled Hubble constant */

  optimalFilter->sampleUnits.powerOfTen += 36;
  
  if (f0 == 0)
  {
    /* calculate lambda */
    /* find omegaRef */
    xRef = (UINT4)((fRef)/deltaF);
    if ( ((fRef)/deltaF) == (REAL8)(xRef) )
    {
      omegaRef = input->omegaGW->data->data[xRef];
    }
    else
    {
      omega1   = input->omegaGW->data->data[xRef];
      omega2   = input->omegaGW->data->data[xRef+1];
      freq1    = xRef * deltaF; 
      freq2    = (xRef + 1) * deltaF;
      if (omega1 <= 0)
      {
        ABORT( status,
         STOCHASTICCROSSCORRELATIONH_ENONPOSOMEGA,
         STOCHASTICCROSSCORRELATIONH_MSGENONPOSOMEGA );
      }
      else
      {
        exponent = ((log((fRef)/freq1))/(log(freq2/freq1)));
      }
      omegaRef = omega1*(pow((omega2/omega1),exponent));
    }
    
    /* calculate inverse lambda value */
    lambdaInv = 0.0; 
    for (i = 1; i < length; ++i)
    {
      f = f0 + deltaF * (REAL8) i;

      f3 = f * f * f;
      f6 = f3 * f3;

      omegaTimesGamma = input->omegaGW->data->data[i]
        * input->overlapReductionFunction->data->data[i];
      p1Inv = input->unWhitenedInverseNoisePSD1->data->data[i];
      p2Inv = input->unWhitenedInverseNoisePSD2->data->data[i];

      lambdaInv += 2.0 * omegaTimesGamma * omegaTimesGamma * p1Inv * p2Inv 
        / f6;
    }

    lambdaInv /= omegaRef 
      * ( (20.0L * LAL_PI * LAL_PI) 
          / ( 3.0L * (LAL_H0FAC_SI*1e+18) * (LAL_H0FAC_SI*1e+18) )
          );

    optimalFilter->data->data[0].re = 0;
    optimalFilter->data->data[0].im = 0;
      
    /* calculate optimal filter values */  
    for (i = 1; i < length; ++i)
    {
      f = f0 + deltaF * (REAL8) i;

      f3 = f * f * f;

      omega = input->omegaGW->data->data[i];
      gamma = input->overlapReductionFunction->data->data[i];
      p1HWInv = input->halfWhitenedInverseNoisePSD1->data->data[i];
      p2HWInv = input->halfWhitenedInverseNoisePSD2->data->data[i];
      
      cPtrOptimalFilter = &(optimalFilter->data->data[i]);

      realFactor = (gamma * omega)/(f3 * lambdaInv);

      cPtrOptimalFilter->re = realFactor * 
        ( (p1HWInv.re) * (p2HWInv.re) + (p1HWInv.im) * (p2HWInv.im) );

      cPtrOptimalFilter->im = realFactor * 
        ( (p1HWInv.re) * (p2HWInv.im) - (p1HWInv.im) * (p2HWInv.re) );
    }

  } /* if f0 == 0 */
  else
  {
    /*****This abort should be replaced with the correct *****/
       /***** non-zero heterodyne frequency procedure******/
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENOTYETHETERO,
          STOCHASTICCROSSCORRELATIONH_MSGENOTYETHETERO); 
  }

  DETATCHSTATUSPTR(status);
  RETURN(status);

}
