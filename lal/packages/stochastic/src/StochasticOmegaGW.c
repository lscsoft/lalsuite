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

/************************************ <lalVerbatim file="StochasticOmegaGWCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{StochasticOmegaGW.c}}
\label{stochastic:ss:StochasticOmegaGW.c}

Generates a frequency series containing a simple power law spectrum.

\subsubsection*{Prototypes}
\idx{LALStochasticOmegaGW()}
\input{StochasticOmegaGWCP}

\subsubsection*{Description}
The strength of a stochastic gravitational wave background is
defined as
\begin{equation}
\Omega_{\scriptstyle{\rm GW}}(f)
:=\frac{f}{\rho_{\scriptstyle{\rm crit}}}
\frac{d\rho_{\scriptstyle{\rm GW}}}{df}
\ ,
\end{equation}
where
\begin{equation}
\rho_{\scriptstyle{\rm crit}} = \frac{3 H_0^2 c^2}{8\pi G}
\end{equation}
is the critical density needed to close the universe.  Since the value
of $\rho_{\scriptstyle{\rm crit}}$ depends on the observed value of the
Hubble constant $H_0$, it is traditional to remove this experimental
uncertainty from the definition by working with
${h_{100}}^2\Omega_{\scriptstyle{\rm GW}}$, where $h_{100}$ is the Hubble
constant divided by
$100\,\textrm{km}\,\textrm{s}^{-1}\,\textrm{Mpc}^{-1}$.

\texttt{LALStochasticOmegaGW()} generates a simple power law spectrum
\begin{equation}
{h_{100}}^2\Omega_{\scriptstyle{\rm GW}}(f)
=\Omega_{\scriptstyle{\rm R}}
\left(
  \frac{f}{f_{\scriptstyle{\rm R}}}
\right)^\alpha
\end{equation}
The parameter \texttt{parameters.omegaRef} specifies the amplitude
${h_{100}}^2\Omega_{\scriptstyle{\rm R}}$ of the spectrum.  This is simply
defined as the value of ${h_{100}}^2\Omega_{\scriptstyle{\rm GW}}(f)$ at the
reference frequency $f_{\scriptstyle{\rm R}}$ which is specified in
\texttt{parameters.omegaRef}.

\subsubsection*{Algorithm}

\texttt{LALStochasticOmegaGW()} treats the constant spectrum $\alpha=0$ as a
special case, and simply sets every element of the output series to
$\Omega_{\scriptstyle{\rm R}}$.

If $\alpha\ne 0$, the DC value \texttt{output->data->data[0]} is set
to 0 or \texttt{LAL\_REAL4\_MAX}, depending on whether $\alpha$ is
positive or negative, respectively.

The output units are set to be dimensionless.

\subsubsection*{Uses}
\begin{verbatim}
LAL_REAL4_MAX
\end{verbatim}

\subsubsection*{Notes}

\begin{itemize}
\item This routine will eventually be generalized to
include ``broken'' power law spectra
\begin{equation}
{h_{100}}^2\Omega_{\scriptstyle{\rm GW}}
= \left\{
\begin{array}{cc}
{h_{100}}^2\Omega_1 f^{\alpha_1} & f\le f_c\\
{h_{100}}^2\Omega_2 f^{\alpha_2} & f\ge f_c
\end{array}
\right.
\end{equation}
\end{itemize}


\vfill{\footnotesize\input{StochasticOmegaGWCV}}

******************************************************* </lalLaTeX> */
/**************************************** <lalLaTeX file="StochasticOmegaGWCB">

% \bibitem{stochastic:}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>

#include <lal/LALConstants.h>
#include <math.h>
#include <lal/StochasticCrossCorrelation.h>

NRCSID (STOCHASTICOMEGAGWC, "$Id$");

/* <lalVerbatim file="StochasticOmegaGWCP"> */
void
LALStochasticOmegaGW(
    LALStatus                         *status,
    REAL4FrequencySeries              *output,
    const StochasticOmegaGWParameters *parameters)
/* </lalVerbatim> */
{
  REAL4* sPtr;
  REAL4* sStopPtr;

  REAL4 alpha;
  REAL8 f0;
  REAL8 deltaF;
  REAL4 x, deltaX, x0;   /* x = f / fRef */
  REAL8 fRef;
  REAL4 omegaRef;
  UINT4 i;

  UINT4 length;

  INITSTATUS(status, "LALStochasticOmegaGW", STOCHASTICOMEGAGWC);

  /* ERROR CHECKING -------------------------------------------------- */

  /* check that pointer to input parameters is non-null */
  ASSERT(parameters != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that frequency spacing is greater than zero */
  deltaF = parameters->deltaF;
  ASSERT(deltaF > 0, status, \
      STOCHASTICCROSSCORRELATIONH_ENONPOSDELTAF, \
      STOCHASTICCROSSCORRELATIONH_MSGENONPOSDELTAF);

  /* check that length parameter is greater than zero */
  length = parameters->length;
  ASSERT(length > 0, status, \
      STOCHASTICCROSSCORRELATIONH_EZEROLEN, \
      STOCHASTICCROSSCORRELATIONH_MSGEZEROLEN);

  /* check that heterodyning doesn't include negative physical frequencies */

  f0 = parameters->f0;
  if (f0 < 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENEGFMIN, \
        STOCHASTICCROSSCORRELATIONH_MSGENEGFMIN);
  }

  /* check that pointer to real frequency series for output is non-null */
  ASSERT(output != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that pointer to data member of real frequency series for
   * output is non-null */
  ASSERT(output->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that length of data member of real frequency series for output
   * equals length specified in input parameters */
  if (output->data->length != length)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EMMLEN, \
        STOCHASTICCROSSCORRELATIONH_MSGEMMLEN);
  }

  /* check that pointer to data-data member of real frequency series for
   * output is non-null */
  ASSERT(output->data->data != NULL, status, \
      STOCHASTICCROSSCORRELATIONH_ENULLPTR, \
      STOCHASTICCROSSCORRELATIONH_MSGENULLPTR);

  /* check that the fRef value is positive */
  fRef = parameters->fRef;
  if (fRef <= 0.0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_EOORFREF, \
        STOCHASTICCROSSCORRELATIONH_MSGEOORFREF);
  }

  /* check that omegaRef is larger than zero */
  omegaRef = parameters->omegaRef;
  if (omegaRef <= 0)
  {
    ABORT(status, STOCHASTICCROSSCORRELATIONH_ENONPOSOMEGA, \
        STOCHASTICCROSSCORRELATIONH_MSGENONPOSOMEGA);
  }

  /* EVERYTHING OKAY HERE! -------------------------------------------- */

  alpha = parameters->alpha;

  /* assign output to be dimensionless */
  output->sampleUnits = lalDimensionlessUnit;

  strncpy(output->name, "Gravitational wave strength OmegaGW", LALNameLength);

  /* assign parameters to frequency series */
  output->epoch.gpsSeconds = 0;
  output->epoch.gpsNanoSeconds = 0;
  output->f0 = f0;
  output->deltaF = deltaF;

  deltaX = deltaF / fRef;
  x0 = f0 / fRef;

  /* assign pointers */
  sStopPtr = output->data->data + length;

  /* calculate output(f) values */
  if (alpha == 0){
    for (sPtr = output->data->data; sPtr < sStopPtr; ++sPtr)
    {
      *sPtr = omegaRef;
    }
  }
  else
  {
    if (f0 == 0)
    {
      output->data->data[0] = (alpha>0 ? 0 : LAL_REAL4_MAX);
      for (i=1 ; i < length ; ++i)
      {
        x = deltaX * (REAL4) i;
        output->data->data[i] = omegaRef * pow(x,alpha);
      }
    }
    else
    {
      for (i=0 ; i < length ; ++i)
      {
        x = x0 + deltaX * (REAL4) i;
        output->data->data[i] = omegaRef * pow(x,alpha);
      }
    }
  }

  RETURN(status);
}
