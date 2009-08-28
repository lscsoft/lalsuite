/*
*  Copyright (C) 2007 Jolien Creighton, Robert Adam Mercer, John Whelan
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/******************* <lalVerbatim file="OverlapReductionFunctionCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{OverlapReductionFunction.c}}
\label{stochastic:ss:OverlapReductionFunction.c}

Calculates the values of the overlap reduction function for a pair
of gravitational wave detectors.

\subsubsection*{Prototypes}
\idx{LALOverlapReductionFunction()}
\input{OverlapReductionFunctionCP}

\subsubsection*{Description}

\texttt{LALOverlapReductionFunction()} calculates the values of the overlap reduction:
%
\begin{equation}
\gamma(f):={5\over 8\pi}\sum_A\int_{S^2}d\hat\Omega\
e^{i2\pi f\hat\Omega\cdot\Delta \vec x/c}\
F_1^A(\hat\Omega)F_2^A(\hat\Omega)\ ,
\label{stochastic:e:gamma(f)}
\end{equation}
%
where $\hat \Omega$ is a unit vector specifying a direction on the
two-sphere, $\Delta\vec x:=\vec x_1-\vec x_2$ is the separation vector
between the two detectors, and
%
\begin{equation}
F_i^A(\hat\Omega):=e_{ab}^A(\hat\Omega)\ d_i^{ab}
% :=
% e_{ab}^A(\hat\Omega)\ {1\over 2}\left(\hat X_i^a\hat X_i^b-
% \hat Y_i^a\hat Y_i^b\right)\
\label{stochastic:e:F_i^A}
\end{equation}
%
is the response of the $i$th detector $(i=1,2)$ to the $A=+,\times$
polarization.
Here $d_i^{ab}$ is the response tensor for the $i$th detector, which
relates the ``strain'' $h$ measured by the detector to the metric
perturbation $h_{ab}$ due to gravitational waves by
%
\begin{equation}
h = d_i^{ab} h_{ab}
\end{equation}
%
The Cartesian components of $d_{ab}$ are constant in an earth-fixed
rotating co\"{o}rdinate system.  $\{e_{ab}^A(\hat\Omega)|A=+,\times\}$
are the spin-2 polarization tensors for the ``plus'' and ``cross''
polarization states, normalized so that $e_{ab}^A e^{Bab}=2\delta^{AB}$.
With this definition,
%
\begin{equation}
\gamma(f)=d_{1ab}d_2^{cd}{5\over 4\pi}\int_{S^2}d\hat\Omega\
e^{i2\pi f\hat\Omega\cdot\Delta \vec x/c}\
P^{ab}_{cd}(\hat\Omega)
\end{equation}
%
where $P^{ab}_{cd}(\hat\Omega)$ is a projection operator onto the space
of symmetric, traceless second-rank tensors orthogonal to $\hat\Omega$.

The overlap reduction function for a pair of identical detectors is a
maximum when they are co\"{\i}ncident and co\"{a}ligned; it
decreases when the detectors are shifted apart (so there is a phase
shift between the signals in the two detectors), or rotated out of
co\"{a}lignment (so the detectors are sensitive to different
polarizations).  The overlap reduction function arises naturally when
calculating the cross-correlated signal due to an isotropic and
unpolarized stochastic gravitational-wave background.

Given a choice of two detector sites, a frequency spacing $\delta f$,
a start frequency $f_0$, and the number of desired values $N$, {\tt
LALOverlapReductionFunction()\/} calculates the values of $\gamma(f)$ at the discrete
frequencies $f_i=f_0 + i\Delta f$, $i=0,1,\cdots, N-1$.

\subsubsection*{Algorithm}

As shown in Appendix B of \cite{stochastic:Flanagan:1993} and Sec.~III.B of
\cite{stochastic:Allen:1999}, the overlap reduction function can be written in
closed form in terms of the traceless response tensor
$D_{ab}=d_{ab}-\delta_{ab} d^c_c/3$
as sum of three spherical Bessel functions:
%
\begin{equation}
\gamma(f)=\rho_1(\alpha)\ D_1^{ab}D_{2ab}
+\rho_2(\alpha)\ D_1^{ab}D_{2a}{}^c s_b s_c
+\rho_3(\alpha)\ D_1^{ab}D_2^{cd}s_a s_b s_c s_d\ ,
\label{stochastic:e:closed1}
\end{equation}
%
where
%
\begin{equation}
\left[
\begin{array}{c}
\rho_1\\
\rho_2\\
\rho_3
\end{array}
\right]
=
{1\over 2\alpha^2}
\left[
\begin{array}{rrr}
 10\alpha^2 & -20\alpha   & 10\\
-20\alpha^2 &  80\alpha   & -100\\
  5\alpha^2 & -50\alpha   & 175
\end{array}
\right]\
\left[
\begin{array}{c}
j_0\\
j_1\\
j_2
\end{array}
\right]\ ,
\label{stochastic:e:closed2}
\end{equation}
%
$j_0$, $j_1$, and $j_2$ are the standard spherical Bessel functions:
%
\begin{eqnarray}
j_0(\alpha)&=&{\sin\alpha\over\alpha} ,\nonumber\\
j_1(\alpha)&=&{\sin\alpha\over\alpha^2}-{\cos\alpha\over\alpha}\ ,
\label{stochastic:e:closed3}\\
j_2(\alpha)&=&3\ {\sin\alpha\over\alpha^3}-3\ {\cos\alpha\over\alpha^2}
-{\sin\alpha\over\alpha}\ ,\nonumber
\end{eqnarray}
%
$\vec s$ is a unit vector pointing in the direction of
$\Delta \vec x:=\vec x_1-\vec x_2$, and $\alpha:=2\pi f|\Delta\vec x|/c$.

{\tt LALOverlapReductionFunction()\/} calculates the values of $\gamma(f)$
as follows:

\begin{enumerate}

\item Gets the locations and response tensors for the two detectors
  from the \texttt{LALDetector} structures in the input.

\item Constructs the traceless parts $D_{iab}$ of the two detector
  response tensors and finds the distance $|\Delta\vec x|$ and
  direction $s^a$ between the sites.

\item Calculates the frequency-independent co\"{e}fficients
  $D_1^{ab}D_{2ab}$, $D_1^{ab}D_{2a}{}^c s_b s_c$, and
  $D_1^{ab}D_2^{cd}s_a s_b s_c s_d$ that appear in
  Eq.~(\ref{stochastic:e:closed1}).

\item Calculates $\gamma(f)$ at each discrete frequency
  $f_i:=f_0+i\Delta f$, $i=0,1,\cdots N-1$, using the power series
  expansion
  \begin{eqnarray}
    j_0(\alpha)&=& 1 - \frac{\alpha^2}{6} + \frac{\alpha^4}{120}
    + \mathcal{O}(\alpha^6) \\
    \frac{j_1(\alpha)}{\alpha}
    &=& \frac{1}{3} - \frac{\alpha^2}{30} + \frac{\alpha^4}{840}
    + \mathcal{O}(\alpha^6) \\
    \frac{j_2(\alpha)}{\alpha^2}
    &=& \frac{1}{15} - \frac{\alpha^2}{210} + \frac{\alpha^4}{7560}
    + \mathcal{O}(\alpha^6)
  \end{eqnarray}
  for the spherical Bessel functions $j_0(\alpha_i)$,
  $j_a(\alpha_i)$, $j_2(\alpha_i)$ when $\alpha_i=2\pi f_i
  |\Delta\vec x|/c<0.01$.

\end{enumerate}

\subsubsection*{Uses}
\begin{verbatim}
LALUnitRaise()
sin()
cos()
sqrt()
strncpy()
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}

\item The $\gamma(f)$ here is related to the unnormalized $\Gamma(f)$
  defined by Maggiore
  \cite{stochastic:Maggiore:2000a,stochastic:Maggiore:2000b} by
  $\gamma(f) = \frac{5}{2}\Gamma(f)$.  This normalization, which
  agrees with the literature
  \cite{stochastic:Flanagan:1993,stochastic:Allen:1997,stochastic:Allen:1999}
  on interferometers, is chosen so that $\gamma(f)\equiv 1$ for a pair
  of co\"{\i}ncident, co\"{a}ligned interferometers with perpendicular
  arms.  It means that, for combinations other than a pair of
  interferometers, our $\gamma(f)$ is \emph{not} equal to the
  generalization of $\gamma(f)$ defined by Maggiore, whose
  relationship to $\Gamma(f)$ depends on the type of detector.
  Defining $\gamma(f)$ as we do allows us to use the formul{\ae} from,
  e.g., \cite{stochastic:Allen:1999}, irrespective of the detector
  type in question.

\item While $\gamma(f)$ is usually considered to be dimensionless,
  this routine attaches to it units of strain$^2$.  This is because it
  contains two powers of the response tensor $d^{ab}$, which converts
  the dimensionless metric perturbation $h_{ab}$ to $h=h_{ab}d^{ab}$,
  which has units of strain.

\end{itemize}

\vfill{\footnotesize\input{OverlapReductionFunctionCV}}

******************************************************* </lalLaTeX> */
/**************************************** <lalLaTeX file="OverlapReductionFunctionCB">
\bibitem{stochastic:Flanagan:1993}
  \'{E}.~\'{E}.~Flanagan, ``The Sensitivity of Ligo to a Stochastic
  Background, and its Dependence on the Detector Orientations''
  Phys.\ Rev.\ D.\ {\bf 48}, 2389 (1993);
  \href{http://www.arXiv.org/abs/astro-ph/9305029}{astro-ph/9305029}
% \bibitem{stochastic:Allen:1997}
%   B.~Allen
%   ``The stochastic gravity-wave background: sources and detection''
%   in \textit{Proceedings of the Les Houches School on Astrophysical Sources of
%   Gravitational Waves},
%   eds. J.~A.~Marck and J.~P.~Lasota, Cambridge, 373 (1997);
%   \href{http://www.arXiv.org/abs/gr-qc/9604033}{gr-qc/9604033}
% \bibitem{stochastic:Allen:1999}
%   B.~Allen and J.D.~Romano, ``Detecting a stochastic background of
%   gravitational radiation: Signal processing strategies and
%   sensitivities''
%   Phys.\ Rev.\ D {\bf 59}, 102001 (1999);
%   \href{http://www.arXiv.org/abs/gr-qc/9710117}{gr-qc/9710117}
\bibitem{stochastic:Maggiore:2000a}
  M.~Maggiore and A.~Nicholis, ``Detection strategies for scalar
  gravitational waves with interferometers and resonant spheres'',
  Phys.\ Rev.\ D \textbf{62}, 024004 (2000);
  \href{http://www.arXiv.org/abs/gr-qc/9907055}{gr-qc/9907055}
\bibitem{stochastic:Maggiore:2000b}
  M.~Maggiore, ``Gravitational Wave Experiments and Early Universe Cosmology'',
  Phys.\ Rept.\ \textbf{331}, 283-367 (2000);
  \href{http://www.arXiv.org/abs/gr-qc/9909001}{gr-qc/9909001}
******************************************************* </lalLaTeX> */
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <math.h>
#include <string.h>
#include <lal/StochasticCrossCorrelation.h>

NRCSID(OVERLAPREDUCTIONFUNCTIONC, "$Id$");

static void evaluateBessels(REAL4 rho[3], REAL4 alpha);
static REAL4 cartesianInnerProduct(REAL4 a[3], REAL4 b[3]);

/* <lalVerbatim file="OverlapReductionFunctionCP"> */
void
LALOverlapReductionFunction(
    LALStatus                                *status,
    REAL4FrequencySeries                     *output,
    const LALDetectorPair                    *detectors,
    const OverlapReductionFunctionParameters *parameters)
/* </lalVerbatim> */
{
  UINT4 length;
  REAL8 deltaF;
  REAL8 f0;
  UINT4 i, j;
  REAL4 trace1, trace2;
  REAL4 d1DotS[3], d2DotS[3];
  REAL4 s[3];
  REAL4 distance;
  REAL4 c1, c2, c3;
  REAL4 alpha, alpha0, deltaAlpha;
  REAL4 rho[3];
  REAL4 d1[3][3], d2[3][3];
  RAT4 power;

  /* initialize status structure */
  INITSTATUS(status, "LALOverlapReductionFunction", OVERLAPREDUCTIONFUNCTIONC);
  ATTATCHSTATUSPTR(status);

  /* check that pointer to parameters is not null */
  ASSERT(parameters!=NULL, status, \
			STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that specified length of output vector is > 0 */
  length = parameters->length;
  ASSERT(length > 0, status, \
			STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that frequency spacing is > 0 */
  deltaF = parameters->deltaF;
  ASSERT(deltaF > 0, status, STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check that minimum frequency is >= 0 */
  f0 = parameters->f0;
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
				STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* check that pointer to output frequency series is not null */
  ASSERT(output != NULL, status, \
			STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
			STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of output frequency series is
	 * not null */
  ASSERT(output->data != NULL, status, \
			STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
			STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of the data member of output frequency series
	 * agrees with length specified in input parameters */
  if(output->data->length != length)
	{
		ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
				STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to data-data member of output vector is not null */
  ASSERT(output->data->data != NULL, status, \
			STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
			STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to input structure is not null */
  ASSERT(detectors != NULL, status, \
			STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
			STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* everything okay here -------------------------------------------- */

  /* the overlap reduction function has units of strain^2 */
  power.numerator = 2;
  power.denominatorMinusOne = 0;
  TRY(LALUnitRaise(status->statusPtr, &(output->sampleUnits), \
				&lalStrainUnit, &power), status);

  /* set parameters for output */
  strncpy(output->name, "Overlap reduction function", LALNameLength);
  output->epoch.gpsSeconds = 0.0;
  output->epoch.gpsNanoSeconds = 0.0;
  output->deltaF = parameters->deltaF;
  output->f0 = parameters->f0;

  /* calculate separation vector between sites */
  for (i = 0; i < 3; i++)
	{
		s[i] = (REAL4)(detectors->detectorOne.location[i] - \
				detectors->detectorTwo.location[i]);
		d1[i][i] = detectors->detectorOne.response[i][i];
    d2[i][i] = detectors->detectorTwo.response[i][i];

		for (j = i; j<3; j++) {
      d1[i][j] = d1[j][i] = detectors->detectorOne.response[i][j];
      d2[i][j] = d2[j][i] = detectors->detectorTwo.response[i][j];

			/* check for non symmetric response tensor */
      ASSERT(d1[j][i] == detectors->detectorOne.response[j][i], status, \
					STOCHASTICCROSSCORRELATIONH_ENONSYMDIJ, \
					STOCHASTICCROSSCORRELATIONH_MSGENONSYMDIJ);
      ASSERT(d2[j][i] == detectors->detectorTwo.response[j][i], status, \
					STOCHASTICCROSSCORRELATIONH_ENONSYMDIJ, \
					STOCHASTICCROSSCORRELATIONH_MSGENONSYMDIJ);
		}
	}

	/* calculate distance between sites, in meters */
  distance = sqrt(cartesianInnerProduct(s,s));

  /* calculate unit separation vector */
  if (distance != 0)
	{
		for (i = 0; i < 3; i++)
		{
      s[i] /= distance;
    }
  }

  trace1 = d1[0][0] + d1[1][1] + d1[2][2];
  if (trace1)
	{
    trace1 /= 3.0;
    for (i = 0; i < 3; i++)
		{
      d1[i][i] -= trace1;
    }
  }
  trace2 = d2[0][0] + d2[1][1] + d2[2][2];
  if (trace2)
	{
    trace2 /= 3.0;
    for (i = 0; i < 3; i++)
		{
      d2[i][i] -= trace2;
    }
  }

  /* calculate coefficients c1, c2, c3 for overlap reduction funtion */

  /* c1 = d1 : d2 */
  c1 = 0;
  for (i = 0; i < 3; i++)
	{
    for (j = 0; j < 3; j++)
		{
      c1 += d1[i][j] * d2[i][j];
    }
  }

  for (i = 0; i < 3; i++)
	{
    d1DotS[i] = cartesianInnerProduct(d1[i], s);
    d2DotS[i] = cartesianInnerProduct(d2[i], s);
  }

  /* c2 = s . d1 . d2 . s */
  c2 = cartesianInnerProduct(d1DotS, d2DotS);

  /* c3 = (s . d1 . s)(s . d2 . s) */
  c3 = cartesianInnerProduct(s, d1DotS) * cartesianInnerProduct(s, d2DotS);

  distance *= (2 * LAL_PI / LAL_C_SI);
  deltaAlpha = deltaF * distance;
  alpha0 = f0 * distance;

  if (f0 == 0)
  {
    for (i = 0; i < length; ++i)
    {
      alpha = deltaAlpha * (REAL4)i;
      evaluateBessels(rho, alpha);
      output->data->data[i] = c1 * rho[0] + c2 * rho[1] + c3 * rho[2];
    }
  }
  else
  {
    for (i = 0; i < length; ++i)
    {
      alpha = alpha0 + deltaAlpha * (REAL4)i;
      evaluateBessels(rho, alpha);
      output->data->data[i] = c1 * rho[0] + c2 * rho[1] + c3 * rho[2];
    }
  }

  /* normal exit */
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

static void evaluateBessels(REAL4 rho[3], REAL4 alpha)
{
  REAL8 alpha2, alpha4;
  REAL8 s, c;
  REAL8 b0, b1, b2;

  alpha2 = alpha * alpha;

  /* if the argument is close to zero, use power series */
  if (alpha < 0.01)
	{
		alpha4 = alpha2 * alpha2;
    b0 = 1.0 - alpha2/6.0 + alpha4/120.0;
    b1 = 1.0/3.0 - alpha2/30.0 + alpha4/840.0;
    b2 = 1.0/15.0 - alpha2/210.0 + alpha4/7560.0;
  }
  else
	{
		s = sin(alpha);
    c = cos(alpha);

    /* define spherical bessel functions j0, j1, j2 */

    b0 = s/alpha; /* = j0 */
    b1 = (b0 - c)/alpha;
    b2 = ((3.0 * b1) - s) / alpha;
    b1 /= alpha; /* = j1/alpha */
    b2 /= alpha2; /* = j2/alpha2 */
  }

  rho[0] = 5.0*b0 - 10.0*b1 + 5.0*b2;
  rho[1] = -10.0*b0 + 40.0*b1 - 50.0*b2;
  rho[2] = 2.5*b0 - 25.0*b1 + 87.5*b2;

  return;
}

static REAL4 cartesianInnerProduct(REAL4 a[3], REAL4 b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
