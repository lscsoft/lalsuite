/*************************** <lalVerbatim file="StochasticOptimalFilterNormalizationCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{StochasticOptimalFilterNormalization.c}}
\label{stochastic:ss:StochasticOptimalFilterNormalization.c}

Calculates the normalization factor for the optimal filter and the
expected variance per unit time of the standard cross-correlation
statistic.

\subsubsection*{Prototypes}
\idx{LALStochasticOptimalFilterNormalization()}
\input{StochasticOptimalFilterNormalizationCP}

\subsubsection*{Description}

As described in Section~\ref{stochastic:ss:StochasticOptimalFilter.c},
the optimal filter for stochastic searches is defined as
%
\begin{equation}\label{stochastic:e:Q}
\widetilde{Q}{}^{\scriptstyle{\rm C}}(f)=\lambda\,
\frac{\gamma(f)\,\Omega_{\scriptstyle{\rm GW}}(f)}
{|f|^3\,P^{\scriptstyle{\rm C}}_1(f)\,P^{\scriptstyle{\rm C}}_2(f)}
\end{equation}
%
The normalization constant $\lambda$ is chosen so that the expected mean
value of the cross-correlation statistic is \cite{stochastic:Allen:1999}
\begin{equation}
\mu = \frac{3 {H_0}^2}{20\pi^2}\, T \,\overline{w_1w_2}
\int_{-\infty}^{\infty} df\, |f|^{-3}\,
\gamma(f)\,\Omega_{\scriptstyle{\rm GW}}(f) 
\widetilde{Q}{}^{\scriptstyle{\rm C}}(f) = \Omega_{\scriptstyle{\rm R}} T
\label{stochastic:e:mu}
\end{equation}
where $T$ is the integration time
(\textit{cf.}~(\ref{stochastic:e:ymax})), $w_1$ and $w_2$ are the functions
used to window the data, and 
$\Omega_{\scriptstyle{\rm R}} =\Omega_{\scriptstyle{\rm
    GW}}(f_{\scriptstyle{\rm R}})$ is the overall strength of the
stochastic background (see
Sec.~\ref{stochastic:ss:OverlapReductionFunction.c}).  This sets the
value at
\begin{equation}
\label{stochastic:e:lambda}
\lambda = \frac{20\pi^2\, \Omega_{\scriptstyle{\rm R}}}
               {3\,{H_0}^2 \overline{w_1w_2}}
\left(
      \int_{-\infty}^\infty \frac{df}{f^6} 
      \frac{[\gamma(f)\,\Omega_{\scriptstyle{\rm GW}}(f)]^2}{P^{\scriptstyle{\rm C}}_1(f)P^{\scriptstyle{\rm C}}_2(f)}
\right)^{-1}
\end{equation}

The same integral used to calculate $\lambda$ also allows one to
calculate the expected variance per unit integration time of the
cross-correlation statistic, since
\begin{eqnarray}
  \frac{\sigma^2}{T} 
  &=& \frac{\overline{(w_1w_2)^2}}{4T}\int_{-\infty}^{\infty} df
  \, P^{\scriptstyle{\rm C}}_1(f)\, P^{\scriptstyle{\rm C}}_2(f)\,
  \left(
    \widetilde{Q}{}^{\scriptstyle{\rm C}}(f)
  \right)^2
  = \frac{\lambda}{4T} \overline{(w_1w_2)^2} 
\int_{-\infty}^{\infty} \frac{df}{|f|^3}\,\gamma(f)\,
  \Omega_{\scriptstyle{\rm GW}}(f)\,\widetilde{Q}{}^{\scriptstyle{\rm C}}(f) 
  \nonumber
  \\
\label{stochastic:e:variance}
  &=& \frac{5\pi^2}{3 {H_0}^2}
  \,\frac{\overline{(w_1w_2)^2}}{\overline{w_1w_2}}
\,\Omega_{\scriptstyle{\rm R}} \,\lambda
\end{eqnarray}
where we have used (\ref{stochastic:e:Q}) to replace one of the two
factors of $\widetilde{Q}{}^{\scriptstyle{\rm C}}(f)$ and (\ref{stochastic:e:mu}) to replace
the integral.

\texttt{LALStochasticOptimalFilterNormalization()} uses
(\ref{stochastic:e:lambda}) to calculate the normalization constant
$\lambda$ and (\ref{stochastic:e:variance}) to calculate the expected
variance per unit time $\sigma^2/T$ of the cross-correlation
statistic.

\subsubsection*{Algorithm}

The routine \texttt{LALStochasticOptimalFilterNormalization()} first uses
(\ref{stochastic:e:lambda}) to find the normalization constant
$\lambda$ (the amplitude ${h_{100}}^2\Omega_{\scriptstyle{\rm R}}$ is
found by logarithmic interpolation using the reference frequency
$f_{\scriptstyle{\rm R}}$ specified in the parameter structure and the
input series representing ${h_{100}}^2\Omega_{\scriptstyle{\rm GW}}(f)$).

The precise behavior of the normalization depends on the boolean
parameter \verb+parameters->heterodyned+, which indicates whether the
filter is to be used on heterodyned data or not.  In the case of
heterodyned data, the integral is approximated by the sum
\begin{eqnarray}
\lambda &\approx& \frac{20\pi^2\, \Omega_{\scriptstyle{\rm R}}}
                     {3\,{H_0}^2 \overline{w_1w_2}}
\left(
  \delta f\sum_{k=0}^{N-1}
  (f_0 + k\,\delta f)^{-6}
  \frac{(\gamma[k]\,\Omega_{\scriptstyle{\rm GW}}[k])^2}{P^{\scriptstyle{\rm C}}_1[k]P^{\scriptstyle{\rm C}}_2[k]}
\right)^{-1}
\nonumber
\\
&\approx&
\frac{20\pi^2\, \Omega_{\scriptstyle{\rm R}}}{3\,{H_0}^2}
\left(
      \int_{f_0}^{f_0+N\delta f} \frac{df}{f^6} 
      \frac{[\gamma(f)\,\Omega_{\scriptstyle{\rm GW}}(f)]^2}{P^{\scriptstyle{\rm C}}_1(f)P^{\scriptstyle{\rm C}}_2(f)}
\right)^{-1}
\end{eqnarray}
(Leaving out frequencies outside the band is equivalent to assuming
one or both of the noise PSDs $P^{\scriptstyle{\rm C}}_{1,2}(f)$ blows up outside that
range.)

In the case of non-heterodyned data with $f_0=0$, we calculate
\begin{equation}
\lambda \approx \frac{20\pi^2\, \Omega_{\scriptstyle{\rm R}}}
                     {3\,{H_0}^2 \overline{w_1w_2}}
\left(
  \delta f\, 2\ {\mathrm{Re}}  \sum_{k=0 \scriptstyle{\rm or } 1}^{N-1}
  (k\,\delta f)^{-6}
  \frac{(\gamma[k]\,\Omega_{\scriptstyle{\rm GW}}[k])^2}{P^{\scriptstyle{\rm C}}_1[k]P^{\scriptstyle{\rm C}}_2[k]}
\right)^{-1}
\end{equation}
which includes negative frequencies as well.  The difference
between the two is because the cross-correlation statistic appearing
in the definition (\ref{stochastic:e:mu}) is the one calculated by
\texttt{StochasticHeterodynedCrossCorrelationStatistic()} in the case
of heterodyned and \texttt{StochasticCrossCorrelationStatistic()} in
the case of non-heterodyned data.

\subsubsection*{Uses}

\begin{verbatim}
LALUnitMultiply()
LALUnitRaise()
LALUnitCompare()
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}
\item The reference frequency $f_{\scriptstyle{\rm R}}$ must lie
  safely enough in the frequency range of the inputs to allow the
  value of ${h_{100}}^2\Omega_{\scriptstyle{\rm R}}$ to be determined
  by interpolation.
\item The implementation of the optimal filter function given here
  assumes a large observation time continuum-limit approximation.  In
  this limit, the Dirichlet kernels (which appear in an exact
  expression for the standard cross-correlation statistic, when
  evaluated in discrete time \cite{stochastic:Finn:2001}; see also
  the documentation for the module Dirichlet.c in the utilities package)
  may be replaced by Dirac delta functions.
\item The units of the input series are checked for consistency; since
  \cite{stochastic:Allen:1999}
  \begin{equation}
    \langle\widetilde{h}{}^{\scriptstyle{\rm C}}_1(f)^*
    \widetilde{h}{}^{\scriptstyle{\rm C}}_2(f')\rangle
    = \frac{3H_0^2}{20\pi^2}\delta(f-f')
    \gamma(|f|)\Omega_{\scriptstyle{\rm GW}}(|f|)
  \end{equation}
  and
  \begin{equation}
    \langle\widetilde{h}_i(f)^*\widetilde{h}_i(f')\rangle
    = \delta(f-f')P^{\scriptstyle{\rm C}}_i(f)
    \,
  \end{equation}
  we demand that, up to powers of ten,
  \begin{equation}
    ([\gamma][\Omega_{\scriptstyle{\rm GW}}])^2
    =([\widetilde{h}_1][\widetilde{h}_2][T]^{-2})^2
    =[P^{\scriptstyle{\rm C}}_1][P^{\scriptstyle{\rm C}}_2][T]^{-2}
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
  constructing the unit structures of its outputs,
  \texttt{LALStochasticOptimalFilterNormalization()} uses the
  power-of-ten feature of the \texttt{LALUnit} structure to account
  for the units of $H_0$.
\item The expected units for the inputs and outputs of this function
  are as follows (although the actual output units will be constructed
  from the input units):
  \begin{equation}
    ([\gamma][\Omega_{\scriptstyle{\rm GW}}])^2
    =([\widetilde{h}_1][\widetilde{h}_2][T]^{-2})^2
    =[P^{\scriptstyle{\rm C}}_1][P^{\scriptstyle{\rm C}}_2][T]^{-2}
  \end{equation}
  \begin{eqnarray}
    {} [\gamma] &=& \textrm{strain}^{2} \\
    {} [\Omega_{\scriptstyle{\rm GW}}] &=& 1 \\
    {} [1/P^{\scriptstyle{\rm C}}_{1,2}] &=& 10^{36}\,\textrm{Hz}\,\textrm{strain}^{-2} \\
    {} [\lambda] &=&
    10^{36}\, [P^{\scriptstyle{\rm C}}_1]\,[P^{\scriptstyle{\rm C}}_2]
    \,[\gamma]^{-2}\,[\Omega_{\scriptstyle{\rm GW}}]^{-1}\,\textrm{s}^{-3}
    =  10^{-36}\,\textrm{s}^{-1}\\
    {} [\sigma^2/T] &=&
    10^{36}\, [\lambda]\, [\Omega_{\scriptstyle{\rm GW}}] \,\textrm{s}^2
    = \textrm{s}
  \end{eqnarray}
\end{itemize}

\vfill{\footnotesize\input{StochasticOptimalFilterNormalizationCV}}

******************************************************* </lalLaTeX> */ 

/**************************** <lalLaTeX file="StochasticOptimalFilterNormalizationCB">
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

/* <lalVerbatim file="StochasticOptimalFilterNormalizationCP"> */
void
LALStochasticOptimalFilterNormalization(
    LALStatus                                            *status,
    StochasticOptimalFilterNormalizationOutput           *output,
    const StochasticOptimalFilterNormalizationInput      *input,
    const StochasticOptimalFilterNormalizationParameters *parameters)
/* </lalVerbatim> */
{
  REAL4 omegaTimesGamma;
  REAL4 p1Inv;
  REAL4 p2Inv;

  REAL8 f;
  REAL8 f0;
  REAL8 f3;
  REAL8 deltaF;

  /* normalization factor */
  UINT4 i;
  UINT4 xRef;
  REAL4 omega1;
  REAL4 omega2;
  REAL8 freq1;
  REAL8 freq2;
  REAL4 exponent;
  REAL4 omegaRef;   
  REAL8 lambdaInv;
  REAL8 f6; 

  /* windowing */
  REAL4 w2, meanW2, meanW4;
  REAL4 *sPtrW1, *sPtrW2, *sStopPtr;
  UINT4 winLength = 0;

  UINT4 length;

  RAT4 power;
  LALUnitPair unitPair;
  LALUnit tmpUnit1, tmpUnit2, checkUnit;
  REAL4WithUnits *lamPtr;

  /* initialize status pointer */
  INITSTATUS(status, "LALStochasticOptimalFilterNormalization", \
      STOCHASTICOPTIMALFILTERC);
  ATTATCHSTATUSPTR(status);

  /* ERROR CHECKING ----------------------------------------------------- */

  /***** check for null pointers *****/
  /* input structure */
  ASSERT(input != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* output structure */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* parameter structure */
  ASSERT(parameters != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
 
  /* overlap member of input */
  ASSERT(input->overlapReductionFunction != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* omega member of input */
  ASSERT(input->omegaGW != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* inverse noise 1 of input */
  ASSERT(input->inverseNoisePSD1 != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* inverse noise 2 of input */
  ASSERT(input->inverseNoisePSD2 != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of overlap */
  ASSERT(input->overlapReductionFunction->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of omega */
  ASSERT(input->omegaGW->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of inverse noise 1 */
  ASSERT(input->inverseNoisePSD1->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data member of inverse noise 2 */
  ASSERT(input->inverseNoisePSD2->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of overlap */
  ASSERT(input->overlapReductionFunction->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of omega */
  ASSERT(input->omegaGW->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of inverse noise 1 */
  ASSERT(input->inverseNoisePSD1->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* data-data member of inverse noise 2 */
  ASSERT(input->inverseNoisePSD2->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* Checks that only apply if windowing is specified */

  if (parameters->window1 || parameters->window2)
  {
    /* window 1 parameter */
    ASSERT(parameters->window1 != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    
    /* window 2 parameter */
    ASSERT(parameters->window2 != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    
    /* data member of window 1 */
    ASSERT(parameters->window1->data != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

    /* data member of window 2 */
    ASSERT(parameters->window2->data != NULL, status, \
        STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
        STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);
    
    winLength = parameters->window1->length;

    ASSERT(winLength != 0, status, \
        STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);
    
    /* Check that windows are the same length */
    if (parameters->window2->length != winLength)
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
          STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
    }
  }

  /*** done with null pointers ***/

  /* extract parameters from overlap */
  length = input->overlapReductionFunction->data->length;
  f0 = input->overlapReductionFunction->f0;
  deltaF = input->overlapReductionFunction->deltaF;

  /**** check for legality ****/
  /* length must be positive */
  ASSERT(length != 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* start frequency must not be negative */
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* frequency spacing must be positive */
  ASSERT(deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /** check for mismatches **/
  /* length */
  if (input->omegaGW->data->length != length) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->inverseNoisePSD1->data->length != length) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }
  if (input->inverseNoisePSD2->data->length != length) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* initial frequency */
  if (input->omegaGW->f0 != f0) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }
  if (input->inverseNoisePSD1->f0 != f0) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }
  if (input->inverseNoisePSD2->f0 != f0) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMFMIN);
  }

  /* frequency spacing */
  if (input->omegaGW->deltaF != deltaF) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  if (input->inverseNoisePSD1->deltaF != deltaF) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  if (input->inverseNoisePSD2->deltaF != deltaF) 
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMDELTAF, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMDELTAF);
  }
  /** check for reference frequency lower and upper limits **/ 
  if (parameters->fRef < f0 + deltaF)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EOORFREF, \
        STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  }
  if (parameters->fRef > f0 + ((REAL8)(length-1)*deltaF))
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EOORFREF, \
        STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  }

  /* EVERYHTING OKAY HERE! ---------------------------------------------- */

  /* All the powers we use are integers, so we can do this once here */
  power.denominatorMinusOne = 0;

  /* check that units of gamma, Omega, P1 and P2 are consistent
   * we must have (gamma*Omega)^2 with the same units as f^2*P1*P2
   * up to a power of ten. We check this by constructing
   * checkUnit = f^2*P1*P2*(gamma*Omega)^-2 */
  unitPair.unitOne = &(input->inverseNoisePSD1->sampleUnits);
  unitPair.unitTwo = &(input->inverseNoisePSD2->sampleUnits);
  TRY(LALUnitMultiply(status->statusPtr, &tmpUnit1, &unitPair), status);

  /* tmpUnit1 holds units of 1/(P1*P2) */

  power.numerator = -1;
  TRY(LALUnitRaise(status->statusPtr, &tmpUnit2, &tmpUnit1, &power), status);

  /* tmpUnit2 holds units of P1*P2 */

  unitPair.unitOne = &(input->overlapReductionFunction->sampleUnits);
  unitPair.unitTwo = &(input->omegaGW->sampleUnits);
  TRY(LALUnitMultiply(status->statusPtr, &tmpUnit1, &unitPair), status);

  /* tmpUnit1 holds units of Omega*Gamma */

  unitPair.unitOne = &tmpUnit1;
  unitPair.unitTwo = &lalSecondUnit;
  TRY(LALUnitMultiply(status->statusPtr, &checkUnit, &unitPair), status);

  /* checkUnit holds units of f^-1*Omega*Gamma */

  power.numerator = -2;
  TRY(LALUnitRaise(status->statusPtr, &tmpUnit1, &checkUnit, &power), status);

  /* tmpUnit1 holds units of f^2*(Omega*Gamma)^-2 */

  /** unitPair.unitOne = &tmpUnit1; **/
  /** Commented out because it's redundant with the last assignment **/
  
  unitPair.unitTwo = &tmpUnit2;

  TRY(LALUnitMultiply(status->statusPtr, &checkUnit, &unitPair), status);

  /* checkUnit holds units of f^2*P1*P2(Omega*Gamma)^-2 */

  /*** Check that checkUnit is dimensionless up to a power of ten ***/

  for (i=0; i<LALNumUnits; ++i)
  {
    if (checkUnit.unitNumerator[i] || checkUnit.unitDenominatorMinusOne[i])
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_EWRONGUNITS, \
          STOCHASTICCROSSCORRELATIONH_MSGEWRONGUNITS);
    }
  }

  /******* Set tmpUnit1 to dims of Omega/H0^2 ******/

  /* First, set it to dims of H0 */

  tmpUnit1 = lalHertzUnit;

  /* Account for scaled units of Hubble constant */
  tmpUnit1.powerOfTen -= 18;

  /* Now set tmpUnit2 to dims of H0^-2 */
  power.numerator = -2;
  TRY(LALUnitRaise(status->statusPtr, &tmpUnit2, &tmpUnit1, &power), status);
  unitPair.unitOne = &(input->omegaGW->sampleUnits);
  unitPair.unitTwo = &tmpUnit2;
  TRY(LALUnitMultiply(status->statusPtr, &tmpUnit1, &unitPair), status);

  /* Now tmpUnit1 has units of Omega/H0^2 */

  /******* assign correct units to normalization constant ********/
  /* These are Omega/H0^2*f^5*P1*P2*(gamma*Omega)^-2
   * which is the same as tmpUnit1*f^3*checkUnit */
  unitPair.unitOne = &tmpUnit1;
  unitPair.unitTwo = &checkUnit;
  TRY(LALUnitMultiply(status->statusPtr, &tmpUnit2, &unitPair), status);

  /* tmpUnit2 has units of Omega/H0^2*f^2*P1*P2*(gamma*Omega)^-2 */

  /* Now that we've used checkUnit, we don't need it any more and */
  /* can use it for temp storage (of f^3) */
  power.numerator = 3;
  TRY(LALUnitRaise(status->statusPtr, &checkUnit, &lalHertzUnit, &power), \
      status);

  /* In case the normalization output was NULL, we need to allocate it
   * since we still use it as an intermediate; we do so here to
   * minimize the BEGINFAIL/ENDFAIL pairs needed */

  if (output->normalization != NULL) 
  {
    lamPtr = output->normalization;
  }
  else 
  {
    lamPtr = (REAL4WithUnits*)LALMalloc(sizeof(REAL4WithUnits));
  }

  unitPair.unitOne = &tmpUnit2;
  
  /* unitPair.unitTwo = &checkUnit; */
  /** Commented out because it's redundant with the last assignment **/

  LALUnitMultiply(status->statusPtr, &(lamPtr->units), &unitPair);

  BEGINFAIL(status) 
    if (output->normalization == NULL) LALFree(lamPtr);
  ENDFAIL(status);

  if (output->variance != NULL) 
  {
    /******* assign correct units to variance per time of CC stat ********/
    unitPair.unitOne = &(lamPtr->units);
    unitPair.unitTwo = &tmpUnit1;
    
    LALUnitMultiply(status->statusPtr, &(output->variance->units), &unitPair);
    
    BEGINFAIL(status) 
      if (output->normalization == NULL) LALFree(lamPtr);
    ENDFAIL(status);
  }

  /* Calculate properties of windows */
  if (parameters->window1 == NULL)
  {
    meanW2 = meanW4 = 1.0;
  }
  else
  {
    meanW2 = meanW4 = 0.0;
    for (sPtrW1 = parameters->window1->data, \
        sPtrW2 = parameters->window2->data, sStopPtr = sPtrW1 + winLength; \
        sPtrW1 < sStopPtr; ++sPtrW1, ++sPtrW2)
    {
      w2 = *sPtrW1 * *sPtrW2;
      meanW2 += w2;
      meanW4 += w2 * w2;
    }
    meanW2 /= winLength;
    meanW4 /= winLength;

    if (meanW2 <= 0)
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_ENONPOSWIN, \
            STOCHASTICCROSSCORRELATIONH_MSGENONPOSWIN);
    }
  }

  /************** calculate lambda ***********************/
  /* find omegaRef */
  xRef = (UINT4)((parameters->fRef - f0)/deltaF);
  if (((parameters->fRef - f0)/deltaF) == (REAL8)(xRef))
  {
    omegaRef = input->omegaGW->data->data[xRef];
  }
  else
  {
    omega1 = input->omegaGW->data->data[xRef];
    omega2 = input->omegaGW->data->data[xRef+1];
    freq1 = f0 + xRef * deltaF; 
    freq2 = f0 + (xRef + 1) * deltaF;
    if (omega1 <= 0)
    {
      ABORT(status, STOCHASTICCROSSCORRELATIONH_ENONPOSOMEGA, \
          STOCHASTICCROSSCORRELATIONH_MSGENONPOSOMEGA);
    }
    else
    {
      exponent = ((log((parameters->fRef)/freq1))/(log(freq2/freq1)));
    }
    omegaRef = omega1*(pow((omega2/omega1),exponent));
  }

  /* calculate inverse lambda value */
  lambdaInv = 0.0; 

  for (i = (f0 == 0 ? 1 : 0) ; i < length; ++i)
  {
    f = f0 + deltaF * (REAL8) i;

    f3 = f * f * f;
    f6 = f3 * f3;
    
    omegaTimesGamma = input->omegaGW->data->data[i] * \
                      input->overlapReductionFunction->data->data[i];
    p1Inv = input->inverseNoisePSD1->data->data[i];
    p2Inv = input->inverseNoisePSD2->data->data[i];
    lambdaInv += (omegaTimesGamma * omegaTimesGamma * p1Inv * p2Inv) / f6;
  }

  lambdaInv /= (omegaRef / deltaF) * ((20.0L * LAL_PI * LAL_PI) / \
      (3.0L * (LAL_H0FAC_SI*1e+18) * (LAL_H0FAC_SI*1e+18) * meanW2));

  if (!parameters->heterodyned)
    lambdaInv *= 2.0;
  
  lamPtr->value = 1 / lambdaInv;

  if (output->variance != NULL) 
  {
    output->variance->value = (meanW4/meanW2) * ((5.0L * LAL_PI * LAL_PI) / \
        (3.0L * (LAL_H0FAC_SI*1e+18) * (LAL_H0FAC_SI*1e+18))) * \
                              omegaRef / lambdaInv;
  }
  if (output->normalization == NULL)
    LALFree(lamPtr);
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
