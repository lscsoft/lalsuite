/*
*  Copyright (C) 2007 John Whelan
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

/************************** <lalVerbatim file="CoarseGrainFrequencySeriesCV">
Author: UTB Relativity Group; contact whelan@phys.utb.edu (original by S. Drasco)
$Id$
************************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Module \texttt{CoarseGrainFrequencySeries.c}}
\label{utilities:ss:CoarseGrainFrequencySeries.c}

``Coarse grains'' a frequency series to produce a series with a lower
frequncy resolution.

\subsubsection*{Prototypes}
\input{CoarseGrainFrequencySeriesCP}
\idx{LALSCoarseGrainFrequencySeries()}
\idx{LALDCoarseGrainFrequencySeries()}
\idx{LALCCoarseGrainFrequencySeries()}
\idx{LALZCoarseGrainFrequencySeries()}

\subsubsection*{Description}

These functions are designed to facilitate approximation of integrals
such as
$$
\int g(f)\,h(f)\,df
$$
when $g(f)$ and $h(f)$ are sampled with different frequency
resolutions.  If the frequency resolution were the same for both
functions, e.g., a frequency spacing of $\delta f$ and a start
frequency of $f_0$, so that the $k$th element corresponded to a
frequency $f_k = f_0 + k\delta f$, the approximation would be defined
as
$$
\int g(f)\,h(f)\,df \approx \delta f \sum_k g_k h_k
$$
whose contribution from the $k$th element is\footnote{It is
  important to make the limits of integration symmetric about $f_k$ to
  maintain reality conditions when dealing with Fourier transforms of
  real quantities.}
$$
\int_{f_k-\delta f/2}^{f_k+\delta f/2} g(f)\,h(f)\,df \approx
\delta f g_k h_k
\ .
$$
The central idea in our definitions of coarse graining will thus be
the correspondence
\begin{equation}
  \label{utilities:e:coarse}
  h_k \approx \frac{1}{\delta f}
  \int_{f_k-\delta f/2}^{f_k+\delta f/2} h(f)\,df
\end{equation}

The purpose of this function is to obtain a frequency series $\{h_k\}$
with start frequency $f_0$ and frequency spacing $\delta f$ from a
finer-grained frequency series $\{h'_\ell\}$ with start frequency
$f'_0$ and frequency spacing $\delta f'$.  Focussing on the $k$th
element of the coarse-grained series, which represents a frequency
range from $f_k-\delta f/2$ to $f_k+\delta f/2$, we consider the
elements of the fine-grained series whose frequency ranges overlap
with this.  (Fig.~\ref{utilities:f:coarse})
\begin{figure}[htbp]
  \begin{center}
    \begin{picture}(200,60)(-50,0)
      \put(-50,50){$\ell^{\scriptstyle{\rm min}}_k-1$}
      \put(-10,50){\vector(1,0){20}}
      \put(0,40){\framebox(20,20){}}
      \put(20,40){\framebox(20,20){$\ell^{\scriptstyle{\rm min}}_k$}}
      \put(40,40){\framebox(30,20){$\cdots$}}
      \put(70,40){\framebox(20,20){$\ell^{\scriptstyle{\rm max}}_k$}}
      \put(90,40){\framebox(20,20){}}
      \put(120,50){\vector(-1,0){20}}
      \put(130,50){$\ell^{\scriptstyle{\rm min}}_k+1$}
      \put(13,20){\framebox(90,20){$k$}}
      \put(13,0){\framebox(20,20){$\lambda^{\scriptstyle{\rm min}}_k$}}
      \put(83,0){\framebox(20,20){$\lambda^{\scriptstyle{\rm max}}_k$}}
    \end{picture}
  \end{center}
  \caption{Coarse graining a frequency series}
  \label{utilities:f:coarse}
\end{figure}
We define $\ell^{\scriptstyle{\rm min}}_k$ and $\ell^{\scriptstyle{\rm
    min}}_k$ to be the indices of the first and last elements of
$h'_\ell$ which overlap \emph{completely} with the frequency range
corresponding to $h_k$.  These are most easily defined in terms of
non-integer indices $\lambda^{\scriptstyle{\rm min}}_k$ and
$\lambda^{\scriptstyle{\rm max}}_k$ which correspond to the locations
of fine-grained elements which would exactly reach the edges of the
coarse-grained element with index $k$.  These are defined by
\begin{eqnarray*}
  f_0 + \left(k-\frac{1}{2}\right) \delta f
  &=& f'_0 + \left(\lambda^{\scriptstyle{\rm min}}_k-\frac{1}{2}\right)
  \delta f' \\
  f_0 + \left(k+\frac{1}{2}\right) \delta f
  &=& f'_0 + \left(\lambda^{\scriptstyle{\rm max}}_k+\frac{1}{2}\right)
  \delta f'
\end{eqnarray*}
or, defining the offset $\Omega=(f_0-f'_0)/\delta f'$ and the coarse
graining ratio $\rho = \delta f / \delta f'$,
\begin{eqnarray*}
  \lambda^{\scriptstyle{\rm min}}_k &=&
  \Omega + \left(k-\frac{1}{2}\right) \rho + \frac{1}{2}\\
  \lambda^{\scriptstyle{\rm max}}_k &=&
  \Omega + \left(k+\frac{1}{2}\right) \rho - \frac{1}{2}
\ .
\end{eqnarray*}
Examination of Fig.~\ref{utilities:f:coarse} shows that
$\ell^{\scriptstyle{\rm min}}_k$ is the smallest integer not less than
$\lambda^{\scriptstyle{\rm min}}_k$ and $\ell^{\scriptstyle{\rm
    min}}_k$ is the largest integer not greater than
$\lambda^{\scriptstyle{\rm min}}_k$.

With these definitions, approximating the integral in
(\ref{utilities:e:coarse}) gives
\begin{equation}\label{utilities:e:coarseapprox}
h_k = \frac{1}{\rho}
\left(
  (\ell^{\scriptstyle{\rm min}}_k - \lambda^{\scriptstyle{\rm min}}_k)
  h'_{\ell^{\scriptscriptstyle{\rm min}}_k-1}
  + \sum_{\ell=\ell^{\scriptscriptstyle{\rm min}}_k}
  ^{\ell^{\scriptscriptstyle{\rm max}}_k}
  h'_\ell
  + (\lambda^{\scriptstyle{\rm max}}_k - \ell^{\scriptstyle{\rm max}}_k)
  h'_{\ell^{\scriptscriptstyle{\rm max}}_k+1}
\right)
\end{equation}

In the special case $f_0=f'_0$, we assume both frequency series
represent the independent parts of larger frequency series
$\{h_k|k=-(N-1)\ldots(N-1)\}$ and $\{h'_\ell|\ell=-(N-1)\ldots(N-1)\}$
which obey $h_{-k}=h_k^*$ and $h'_{-\ell}{}=h'_\ell{}^*$ (e.g.,
fourier transforms of real data).  In that case, the DC element of the
coarse-grained series can be built out of both positive- and implied
negative-frequency elements in the fine-grained series.
\begin{equation}
  h_0 = \frac{1}{\rho}
  \left[
    h'_0
    + 2\ \mathrm{Re}
    \left(
      \sum_{\ell=1}^{\ell^{\scriptscriptstyle{\rm max}}_0}
      h'_\ell
      + (\lambda^{\scriptstyle{\rm max}}_0 - \ell^{\scriptstyle{\rm max}}_0)
      h'_{\ell^{\scriptscriptstyle{\rm max}}_0+1}
    \right)
  \right]
\end{equation}

\subsubsection*{Algorithm}

These routines move through the output series, using
(\ref{utilities:e:coarseapprox}) to add up the contributions from the
bins in the fine-grained series.

\subsubsection*{Uses}
\begin{verbatim}
strncpy()
\end{verbatim}

\subsubsection*{Notes}
\begin{itemize}
\item The coarse graining ratio must obey $\rho\ge 1$ (so the
  coarse-grained frequency spacing must be less than the fine-grained
  one).  Additionally, the bins in the fine-grained frequency series
  must \emph{completely} overlap those in the coarse-grained frequency
  series.  In particular, since the lowest frequency in the first bin
  of the coarse-grained series is $f_{\scriptstyle{\rm
      min}}=f_0-\delta f/2$ and the last is $f_{\scriptstyle{\rm
      max}}=f_0 + (N-1) \delta f +\delta f/2$ (taking into account the
  width of the bins), the conitions are
  \begin{eqnarray*}
    f_0 - \frac{\delta f}{2} &\ge& f'_0 - \frac{\delta f'}{2}\\
    f_0 + \left(N-\frac{1}{2}\right)\,\delta f &\le&
     f'_0 + \left(N'-\frac{1}{2}\right)\,\delta f'
  \end{eqnarray*}
  (The special case $f_0=f'_0=0$ is an
  exception to the condition on the minimum frequency.)
\item The routines return an error if either minimum frequency
  ($f_{\scriptstyle{\rm min}}$ or $f'_{\scriptstyle{\rm min}}$) is
  negative (unless $f_0=0$ or $f'_0=0$, respectively).
\end{itemize}

\vfill{\footnotesize\input{CoarseGrainFrequencySeriesCV}}

******************************************************* </lalLaTeX> */
/**************************** <lalLaTeX file="CoarseGrainFrequencySeriesCB">

% \bibitem{utilities:}

******************************************************* </lalLaTeX> */
#include <lal/LALStdlib.h>
#include <lal/Units.h>
#include <lal/CoarseGrainFrequencySeries.h>
#include <math.h>

NRCSID(COARSEGRAINFREQUENCYSERIESC,
       "$Id$");

/* <lalVerbatim file="CoarseGrainFrequencySeriesCP"> */
void
LALSCoarseGrainFrequencySeries(LALStatus                      *status,
                               REAL4FrequencySeries           *output,
                               const REAL4FrequencySeries     *input,
                               const FrequencySamplingParams  *params)
/* </lalVerbatim> */
{
  UINT4         lengthCoarse, lengthFine;
  REAL8         f0Coarse, f0Fine;
  REAL8         fMinCoarse, fMinFine;
  REAL8         deltaFCoarse, deltaFFine;
  REAL8         offset, resRatio;

  UINT4         k, l;
  UINT4         lMin, lMax;
  REAL8         lamMin, lamMax;
  REAL4         value;

  /* initialize status structure */
  INITSTATUS( status, "LALSCoarseGrainFrequencySeries",
              COARSEGRAINFREQUENCYSERIESC );
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /*    output series */
  ASSERT(output != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of output series */
  ASSERT(output->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of output series */
  ASSERT(output->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    input series */
  ASSERT(input != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of input series */
  ASSERT(input->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of input series */
  ASSERT(input->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    parameter structure */
  ASSERT(params != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /* checks for duplicate pointers: */
  ASSERT(output != input, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data != input->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data->data != input->data->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  /* extract coarse-grained parameters  */
  lengthCoarse = params->length;
  f0Coarse     = params->f0;
  deltaFCoarse = params->deltaF;
  fMinCoarse = ( f0Coarse == 0.0 ? 0.0 : f0Coarse - deltaFCoarse / 2.0 );

  /* extract fine-grained parameters */
  lengthFine   = input->data->length;
  f0Fine       = input->f0;
  deltaFFine   = input->deltaF;
  fMinFine = ( f0Fine == 0.0 ? 0.0 : f0Fine - deltaFFine / 2.0 );

  /* check for legality of values */

  /*    length must be positive */

  ASSERT(lengthCoarse != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  ASSERT(lengthFine != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  /*    start frequency must not be negative */

  if (fMinCoarse < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  if (fMinFine < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  /*    frequency spacing must be positive */

  ASSERT(deltaFCoarse > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  ASSERT(deltaFFine > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  /* check for length mismatch */

  if (output->data->length != lengthCoarse)
  {
    ABORT(status,
         COARSEGRAINFREQUENCYSERIESH_EMMLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEMMLEN);
  }

  /* Calculate coarse-graining parameters */

  offset = ( f0Coarse - f0Fine ) / deltaFFine;

  resRatio = deltaFCoarse / deltaFFine;

  /* Check that coarse-graining makes sense */

  /* make sure coarse-grained series is not finer than fine-grained */
  if ( resRatio < 1.0 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure minimum frequency in coarse-grained series is not
     less than minimum frequency in fine-grained series */
  if ( fMinCoarse < fMinFine )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure maximum frequency in coarse-grained series is not
     more than maximum frequency in fine-grained series */
  if ( offset + resRatio * ( (REAL8) lengthCoarse - 0.5 )
       > lengthFine - 0.5 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  if (f0Coarse == 0.0)
  {
    /* DC component */

    lamMax = (resRatio / 2.0) - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    if ( lamMax != (REAL8) lMax )
    {
      value = ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1];
    }
    else {
      value = 0.0;
    }

    for ( l = 1 ; l <= lMax ; ++l)
    {
      value += input->data->data[l];
    }

    output->data->data[0] = ( input->data->data[0] + 2.0 * value )
      / resRatio;

    k = 1;
  }
  else
  {
    k = 0;
  }

  for ( ; k < lengthCoarse ; ++k )
  {
    lamMin = offset + ( (REAL8) k - 0.5 ) * resRatio + 0.5 ;
    lMin = (UINT4) ceil(lamMin);

    if ( lamMin != (REAL8) lMin ) {
      value = ( (REAL8) lMin - lamMin ) * input->data->data[lMin-1];
    }
    else
    {
      value = 0.0;
    }

    lamMax = offset + ( (REAL8) k + 0.5 ) * resRatio - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    for ( l = lMin ; l <= lMax ; ++l)
    {
      value += input->data->data[l];
    }

    if ( lamMax != (REAL8) lMax ) {
      value += ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1];
    }

    output->data->data[k] = value / resRatio;

  }

  /* Set output properties */
  output->sampleUnits = input->sampleUnits;
  strncpy( output->name, input->name, LALNameLength );
  output->f0 = f0Coarse;
  output->deltaF = deltaFCoarse;
  output->epoch = input->epoch;

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALSCoarseGrainFrequencySeries() */

/* <lalVerbatim file="CoarseGrainFrequencySeriesCP"> */
void
LALDCoarseGrainFrequencySeries(LALStatus                      *status,
                               REAL8FrequencySeries           *output,
                               const REAL8FrequencySeries     *input,
                               const FrequencySamplingParams  *params)
/* </lalVerbatim> */
{
  UINT4         lengthCoarse, lengthFine;
  REAL8         f0Coarse, f0Fine;
  REAL8         fMinCoarse, fMinFine;
  REAL8         deltaFCoarse, deltaFFine;
  REAL8         offset, resRatio;

  UINT4         k, l;
  UINT4         lMin, lMax;
  REAL8         lamMin, lamMax;
  REAL8         value;

  /* initialize status structure */
  INITSTATUS( status, "LALDCoarseGrainFrequencySeries",
              COARSEGRAINFREQUENCYSERIESC );
  ATTATCHSTATUSPTR(status);

  /* checks for null pointers: */

  /*    output series */
  ASSERT(output != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of output series */
  ASSERT(output->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of output series */
  ASSERT(output->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    input series */
  ASSERT(input != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of input series */
  ASSERT(input->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of input series */
  ASSERT(input->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    parameter structure */
  ASSERT(params != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /* checks for duplicate pointers: */
  ASSERT(output != input, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data != input->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data->data != input->data->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  /* extract coarse-grained parameters  */
  lengthCoarse = params->length;
  f0Coarse     = params->f0;
  deltaFCoarse = params->deltaF;
  fMinCoarse = ( f0Coarse == 0.0 ? 0.0 : f0Coarse - deltaFCoarse / 2.0 );

  /* extract fine-grained parameters */
  lengthFine   = input->data->length;
  f0Fine       = input->f0;
  deltaFFine   = input->deltaF;
  fMinFine = ( f0Fine == 0.0 ? 0.0 : f0Fine - deltaFFine / 2.0 );

  /* check for legality of values */

  /*    length must be positive */

  ASSERT(lengthCoarse != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  ASSERT(lengthFine != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  /*    start frequency must not be negative */

  if (fMinCoarse < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  if (fMinFine < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  /*    frequency spacing must be positive */

  ASSERT(deltaFCoarse > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  ASSERT(deltaFFine > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  /* check for length mismatch */

  if (output->data->length != lengthCoarse)
  {
    ABORT(status,
         COARSEGRAINFREQUENCYSERIESH_EMMLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEMMLEN);
  }

  /* Calculate coarse-graining parameters */

  offset = ( f0Coarse - f0Fine ) / deltaFFine;

  resRatio = deltaFCoarse / deltaFFine;

  /* Check that coarse-graining makes sense */

  /* make sure coarse-grained series is not finer than fine-grained */
  if ( resRatio < 1.0 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure minimum frequency in coarse-grained series is not
     less than minimum frequency in fine-grained series */
  if ( fMinCoarse < fMinFine )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure maximum frequency in coarse-grained series is not
     more than maximum frequency in fine-grained series */
  if ( offset + resRatio * ( (REAL8) lengthCoarse - 0.5 )
       > lengthFine - 0.5 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  if (f0Coarse == 0.0)
  {
    /* DC component */

    lamMax = (resRatio / 2.0) - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    if ( lamMax != (REAL8) lMax )
    {
      value = ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1];
    }
    else {
      value = 0.0;
    }

    for ( l = 1 ; l <= lMax ; ++l)
    {
      value += input->data->data[l];
    }

    output->data->data[0] = ( input->data->data[0] + 2.0 * value )
      / resRatio;

    k = 1;
  }
  else
  {
    k = 0;
  }

  for ( ; k < lengthCoarse ; ++k )
  {
    lamMin = offset + ( (REAL8) k - 0.5 ) * resRatio + 0.5 ;
    lMin = (UINT4) ceil(lamMin);

    if ( lamMin != (REAL8) lMin ) {
      value = ( (REAL8) lMin - lamMin ) * input->data->data[lMin-1];
    }
    else
    {
      value = 0.0;
    }

    lamMax = offset + ( (REAL8) k + 0.5 ) * resRatio - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    for ( l = lMin ; l <= lMax ; ++l)
    {
      value += input->data->data[l];
    }

    if ( lamMax != (REAL8) lMax ) {
      value += ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1];
    }

    output->data->data[k] = value / resRatio;

  }

  /* Set output properties */
  output->sampleUnits = input->sampleUnits;
  strncpy( output->name, input->name, LALNameLength );
  output->f0 = f0Coarse;
  output->deltaF = deltaFCoarse;
  output->epoch = input->epoch;

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALDCoarseGrainFrequencySeries() */

/* <lalVerbatim file="CoarseGrainFrequencySeriesCP"> */
void
LALCCoarseGrainFrequencySeries(LALStatus                      *status,
                               COMPLEX8FrequencySeries        *output,
                               const COMPLEX8FrequencySeries  *input,
                               const FrequencySamplingParams  *params)
/* </lalVerbatim> */
{
  UINT4         lengthCoarse, lengthFine;
  REAL8         f0Coarse, f0Fine;
  REAL8         fMinCoarse, fMinFine;
  REAL8         deltaFCoarse, deltaFFine;
  REAL8         offset, resRatio;

  UINT4         k, l;
  UINT4         lMin, lMax;
  REAL8         lamMin, lamMax;
  COMPLEX8      value;

  /* initialize status structure */
  INITSTATUS( status, "LALCCoarseGrainFrequencySeries",
              COARSEGRAINFREQUENCYSERIESC );
  ATTATCHSTATUSPTR(status);

  /* printf("entering function\n"); */

  /* checks for null pointers: */

  /*    output series */
  ASSERT(output != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of output series */
  ASSERT(output->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of output series */
  ASSERT(output->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    input series */
  ASSERT(input != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of input series */
  ASSERT(input->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of input series */
  ASSERT(input->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    parameter structure */
  ASSERT(params != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /* checks for duplicate pointers: */

  ASSERT(output != input, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data != input->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data->data != input->data->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  /* extract coarse-grained parameters  */
  lengthCoarse = params->length;
  f0Coarse     = params->f0;
  deltaFCoarse = params->deltaF;
  fMinCoarse = ( f0Coarse == 0.0 ? 0.0 : f0Coarse - deltaFCoarse / 2.0 );

  /* extract fine-grained parameters */
  lengthFine   = input->data->length;
  f0Fine       = input->f0;
  deltaFFine   = input->deltaF;
  fMinFine = ( f0Fine == 0.0 ? 0.0 : f0Fine - deltaFFine / 2.0 );

  /* check for legality of values */

  /*    length must be positive */

  ASSERT(lengthCoarse != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  ASSERT(lengthFine != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  /*    start frequency must not be negative */

  if (fMinCoarse < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  if (fMinFine < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  /*    frequency spacing must be positive */

  ASSERT(deltaFCoarse > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  ASSERT(deltaFFine > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  /* check for length mismatch */

  if (output->data->length != lengthCoarse)
  {
    ABORT(status,
         COARSEGRAINFREQUENCYSERIESH_EMMLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEMMLEN);
  }

  /* Calculate coarse-graining parameters */

  offset = ( f0Coarse - f0Fine ) / deltaFFine;

  resRatio = deltaFCoarse /deltaFFine;

  /* Check that coarse-graining makes sense */

  /* make sure coarse-grained series is not finer than fine-grained */
  if ( resRatio < 1.0 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure minimum frequency in coarse-grained series is not
     less than minimum frequency in fine-grained series */
  if ( fMinCoarse < fMinFine )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure maximum frequency in coarse-grained series is not
     more than maximum frequency in fine-grained series */
  if ( offset + resRatio * ( (REAL8) lengthCoarse - 0.5 )
       > lengthFine - 0.5 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /*
  printf("survived checks\n");

  printf("res ratio %f/%f = %f\n",deltaFCoarse,deltaFFine,resRatio);
  printf("offset (%f-%f)/%f = %f\n",f0Coarse,f0Fine,deltaFFine,offset);
  */

  if (f0Coarse == 0.0)
  {
    /* DC component */

    lamMax = (resRatio / 2.0) - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    /* printf("%f %d %d\n",lamMax,lMax,lengthFine); */

    if ( lamMax != (REAL8) lMax )
    {
      value.re = ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1].re;
    }
    else {
      value.re = 0.0;
    }

    for ( l = 1 ; l <= lMax ; ++l)
    {
      value.re += input->data->data[l].re;
    }

    output->data->data[0].re = ( input->data->data[0].re + 2.0 * value.re )
      / resRatio;
    output->data->data[0].im = 0.0;

    /* :TODO: ?  check that imaginary parts of DC vanish? */

    k = 1;
  }
  else
  {
    k = 0;
  }

  for ( ; k < lengthCoarse ; ++k )
  {
    lamMin = offset + ( (REAL8) k - 0.5 ) * resRatio + 0.5 ;
    lMin = (UINT4) ceil(lamMin);

    /* printf("%f %d\n",lamMin,lMin); */

    if ( lamMin != (REAL8) lMin ) {
      value.re = ( (REAL8) lMin - lamMin ) * input->data->data[lMin-1].re;
      value.im = ( (REAL8) lMin - lamMin ) * input->data->data[lMin-1].im;
    }
    else
    {
      value.re = value.im = 0.0;
    }

    lamMax = offset + ( (REAL8) k + 0.5 ) * resRatio - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    /*    printf("%f %d %d\n",lamMax,lMax,lengthFine); */

    for ( l = lMin ; l <= lMax ; ++l)
    {
      value.re += input->data->data[l].re;
      value.im += input->data->data[l].im;
    }

    if ( lamMax != (REAL8) lMax ) {
      value.re += ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1].re;
      value.im += ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1].im;
    }

    output->data->data[k].re = value.re / resRatio;
    output->data->data[k].im = value.im / resRatio;

  }

  /* Set output properties */
  output->sampleUnits = input->sampleUnits;
  strncpy( output->name, input->name, LALNameLength );
  output->f0 = f0Coarse;
  output->deltaF = deltaFCoarse;
  output->epoch = input->epoch;

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALCCoarseGrainFrequencySeries() */

/* <lalVerbatim file="CoarseGrainFrequencySeriesCP"> */
void
LALZCoarseGrainFrequencySeries(LALStatus                      *status,
                               COMPLEX16FrequencySeries        *output,
                               const COMPLEX16FrequencySeries  *input,
                               const FrequencySamplingParams  *params)
/* </lalVerbatim> */
{
  UINT4         lengthCoarse, lengthFine;
  REAL8         f0Coarse, f0Fine;
  REAL8         fMinCoarse, fMinFine;
  REAL8         deltaFCoarse, deltaFFine;
  REAL8         offset, resRatio;

  UINT4         k, l;
  UINT4         lMin, lMax;
  REAL8         lamMin, lamMax;
  COMPLEX16      value;

  /* initialize status structure */
  INITSTATUS( status, "LALZCoarseGrainFrequencySeries",
              COARSEGRAINFREQUENCYSERIESC );
  ATTATCHSTATUSPTR(status);

  /* printf("entering function\n"); */

  /* checks for null pointers: */

  /*    output series */
  ASSERT(output != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of output series */
  ASSERT(output->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of output series */
  ASSERT(output->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    input series */
  ASSERT(input != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of input series */
  ASSERT(input->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    data member of data member of input series */
  ASSERT(input->data->data != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /*    parameter structure */
  ASSERT(params != NULL, status,
         COARSEGRAINFREQUENCYSERIESH_ENULLPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR);

  /* checks for duplicate pointers: */

  ASSERT(output != input, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data != input->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  ASSERT(output->data->data != input->data->data, status,
         COARSEGRAINFREQUENCYSERIESH_ESAMEPTR,
         COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR);

  /* extract coarse-grained parameters  */
  lengthCoarse = params->length;
  f0Coarse     = params->f0;
  deltaFCoarse = params->deltaF;
  fMinCoarse = ( f0Coarse == 0.0 ? 0.0 : f0Coarse - deltaFCoarse / 2.0 );

  /* extract fine-grained parameters */
  lengthFine   = input->data->length;
  f0Fine       = input->f0;
  deltaFFine   = input->deltaF;
  fMinFine = ( f0Fine == 0.0 ? 0.0 : f0Fine - deltaFFine / 2.0 );

  /* check for legality of values */

  /*    length must be positive */

  ASSERT(lengthCoarse != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  ASSERT(lengthFine != 0, status,
         COARSEGRAINFREQUENCYSERIESH_EZEROLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN);

  /*    start frequency must not be negative */

  if (fMinCoarse < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  if (fMinFine < 0.0)
  {
    ABORT( status,
         COARSEGRAINFREQUENCYSERIESH_ENEGFMIN,
         COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN );
  }

  /*    frequency spacing must be positive */

  ASSERT(deltaFCoarse > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  ASSERT(deltaFFine > 0.0, status,
         COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF,
         COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF);

  /* check for length mismatch */

  if (output->data->length != lengthCoarse)
  {
    ABORT(status,
         COARSEGRAINFREQUENCYSERIESH_EMMLEN,
         COARSEGRAINFREQUENCYSERIESH_MSGEMMLEN);
  }

  /* Calculate coarse-graining parameters */

  offset = ( f0Coarse - f0Fine ) / deltaFFine;

  resRatio = deltaFCoarse /deltaFFine;

  /* Check that coarse-graining makes sense */

  /* make sure coarse-grained series is not finer than fine-grained */
  if ( resRatio < 1.0 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure minimum frequency in coarse-grained series is not
     less than minimum frequency in fine-grained series */
  if ( fMinCoarse < fMinFine )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /* make sure maximum frequency in coarse-grained series is not
     more than maximum frequency in fine-grained series */
  if ( offset + resRatio * ( (REAL8) lengthCoarse - 0.5 )
       > lengthFine - 0.5 )
  {
    ABORT( status,
           COARSEGRAINFREQUENCYSERIESH_EOORCOARSE,
           COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE );
  }

  /*
  printf("survived checks\n");

  printf("res ratio %f/%f = %f\n",deltaFCoarse,deltaFFine,resRatio);
  printf("offset (%f-%f)/%f = %f\n",f0Coarse,f0Fine,deltaFFine,offset);
  */

  if (f0Coarse == 0.0)
  {
    /* DC component */

    lamMax = (resRatio / 2.0) - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    /* printf("%f %d %d\n",lamMax,lMax,lengthFine); */

    if ( lamMax != (REAL8) lMax )
    {
      value.re = ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1].re;
    }
    else {
      value.re = 0.0;
    }

    for ( l = 1 ; l <= lMax ; ++l)
    {
      value.re += input->data->data[l].re;
    }

    output->data->data[0].re = ( input->data->data[0].re + 2.0 * value.re )
      / resRatio;
    output->data->data[0].im = 0.0;

    /* :TODO: ?  check that imaginary parts of DC vanish? */

    k = 1;
  }
  else
  {
    k = 0;
  }

  for ( ; k < lengthCoarse ; ++k )
  {
    lamMin = offset + ( (REAL8) k - 0.5 ) * resRatio + 0.5 ;
    lMin = (UINT4) ceil(lamMin);

    /* printf("%f %d\n",lamMin,lMin); */

    if ( lamMin != (REAL8) lMin ) {
      value.re = ( (REAL8) lMin - lamMin ) * input->data->data[lMin-1].re;
      value.im = ( (REAL8) lMin - lamMin ) * input->data->data[lMin-1].im;
    }
    else
    {
      value.re = value.im = 0.0;
    }

    lamMax = offset + ( (REAL8) k + 0.5 ) * resRatio - 0.5 ;
    lMax = (UINT4) floor(lamMax);

    /*    printf("%f %d %d\n",lamMax,lMax,lengthFine); */

    for ( l = lMin ; l <= lMax ; ++l)
    {
      value.re += input->data->data[l].re;
      value.im += input->data->data[l].im;
    }

    if ( lamMax != (REAL8) lMax ) {
      value.re += ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1].re;
      value.im += ( lamMax - (REAL8) lMax ) * input->data->data[lMax+1].im;
    }

    output->data->data[k].re = value.re / resRatio;
    output->data->data[k].im = value.im / resRatio;

  }

  /* Set output properties */
  output->sampleUnits = input->sampleUnits;
  strncpy( output->name, input->name, LALNameLength );
  output->f0 = f0Coarse;
  output->deltaF = deltaFCoarse;
  output->epoch = input->epoch;

  DETATCHSTATUSPTR(status);
  RETURN(status);

} /* LALZCoarseGrainFrequencySeries() */
