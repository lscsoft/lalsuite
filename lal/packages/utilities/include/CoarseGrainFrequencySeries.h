/*
*  Copyright (C) 2007 Jolien Creighton, John Whelan
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

#include <lal/LALStdlib.h>

#ifndef _COARSEGRAINFREQUENCYSERIES_H
#define _COARSEGRAINFREQUENCYSERIES_H

#ifdef  __cplusplus
extern "C" {
#endif


/**
   \addtogroup CoarseGrainFrequencySeries_h
   \author UTB Relativity Group; contact whelan@phys.utb.edu (original by S. Drasco)

   \brief Provides prototype, structure and error code information for routines which coarse-grain a frequency series.

   \heading{Synopsis}
   \code
   #include <lal/CoarseGrainFrequencySeries.h>
   \endcode

\heading{Description}

These functions are designed to facilitate approximation of integrals
such as
\f[
\int g(f)\,h(f)\,df
\f]
when \f$g(f)\f$ and \f$h(f)\f$ are sampled with different frequency
resolutions.  If the frequency resolution were the same for both
functions, e.g., a frequency spacing of \f$\delta f\f$ and a start
frequency of \f$f_0\f$, so that the \f$k\f$th element corresponded to a
frequency \f$f_k = f_0 + k\delta f\f$, the approximation would be defined
as
\f[
\int g(f)\,h(f)\,df \approx \delta f \sum_k g_k h_k
\f]
whose contribution from the \f$k\f$th element is [It is
  important to make the limits of integration symmetric about \f$f_k\f$ to
  maintain reality conditions when dealing with Fourier transforms of
  real quantities.]
\f[
\int_{f_k-\delta f/2}^{f_k+\delta f/2} g(f)\,h(f)\,df \approx
\delta f g_k h_k
\ .
\f]
The central idea in our definitions of coarse graining will thus be
the correspondence
\anchor utilities_e_coarse \f{equation}{
  \tag{utilities_e_coarse}
  h_k \approx \frac{1}{\delta f}
  \int_{f_k-\delta f/2}^{f_k+\delta f/2} h(f)\,df
\f}

The purpose of this function is to obtain a frequency series \f$\{h_k\}\f$
with start frequency \f$f_0\f$ and frequency spacing \f$\delta f\f$ from a
finer-grained frequency series \f$\{h'_\ell\}\f$ with start frequency
\f$f'_0\f$ and frequency spacing \f$\delta f'\f$.  Focussing on the \f$k\f$th
element of the coarse-grained series, which represents a frequency
range from \f$f_k-\delta f/2\f$ to \f$f_k+\delta f/2\f$, we consider the
elements of the fine-grained series whose frequency ranges overlap
with this.  (Fig.\figref{utilities_f_coarse}

\floatfig{htbp,utilities_f_coarse}
\image html  utilitiesCoarseGrain.png "Fig. [utilities_f_coarse]: Coarse graining a frequency series"
\image latex utilitiesCoarseGrain.pdf "Coarse graining a frequency series"

We define \f$\ell^{\scriptstyle\textrm{min}}_k\f$ and \f$\ell^{\scriptstyle{\rm  min}}_k\f$
to be the indices of the first and last elements of
\f$h'_\ell\f$ which overlap \e completely with the frequency range
corresponding to \f$h_k\f$.  These are most easily defined in terms of
non-integer indices \f$\lambda^{\scriptstyle\textrm{min}}_k\f$ and
\f$\lambda^{\scriptstyle\textrm{max}}_k\f$ which correspond to the locations
of fine-grained elements which would exactly reach the edges of the
coarse-grained element with index \f$k\f$.  These are defined by
\f{eqnarray*}{
  f_0 + \left(k-\frac{1}{2}\right) \delta f
  &=& f'_0 + \left(\lambda^{\scriptstyle\textrm{min}}_k-\frac{1}{2}\right)
  \delta f' \\
  f_0 + \left(k+\frac{1}{2}\right) \delta f
  &=& f'_0 + \left(\lambda^{\scriptstyle\textrm{max}}_k+\frac{1}{2}\right)
  \delta f'
\f}
or, defining the offset \f$\Omega=(f_0-f'_0)/\delta f'\f$ and the coarse
graining ratio \f$\rho = \delta f / \delta f'\f$,
\f{eqnarray*}{
  \lambda^{\scriptstyle\textrm{min}}_k &=&
  \Omega + \left(k-\frac{1}{2}\right) \rho + \frac{1}{2}\\
  \lambda^{\scriptstyle\textrm{max}}_k &=&
  \Omega + \left(k+\frac{1}{2}\right) \rho - \frac{1}{2}
\ .
\f}
Examination of Fig.\figref{utilities_f_coarse} shows that
\f$\ell^{\scriptstyle\textrm{min}}_k\f$ is the smallest integer not less than
\f$\lambda^{\scriptstyle\textrm{min}}_k\f$ and \f$\ell^{\scriptstyle{\rm
    min}}_k\f$ is the largest integer not greater than
\f$\lambda^{\scriptstyle\textrm{min}}_k\f$.

With these definitions, approximating the integral in
\eqref{utilities_e_coarse} gives
\anchor utilities_e_coarseapprox \f{equation}{\tag{utilities_e_coarseapprox}
h_k = \frac{1}{\rho}
\left(
  (\ell^{\scriptstyle\textrm{min}}_k - \lambda^{\scriptstyle\textrm{min}}_k)
  h'_{\ell^{\scriptscriptstyle\textrm{min}}_k-1}
  + \sum_{\ell=\ell^{\scriptscriptstyle\textrm{min}}_k}
  ^{\ell^{\scriptscriptstyle\textrm{max}}_k}
  h'_\ell
  + (\lambda^{\scriptstyle\textrm{max}}_k - \ell^{\scriptstyle\textrm{max}}_k)
  h'_{\ell^{\scriptscriptstyle\textrm{max}}_k+1}
\right)
\f}

In the special case \f$f_0=f'_0\f$, we assume both frequency series
represent the independent parts of larger frequency series
\f$\{h_k|k=-(N-1)\ldots(N-1)\}\f$ and \f$\{h'_\ell|\ell=-(N-1)\ldots(N-1)\}\f$
which obey \f$h_{-k}=h_k^*\f$ and \f$h'_{-\ell}{}=h'_\ell{}^*\f$ (e.g.,
fourier transforms of real data).  In that case, the DC element of the
coarse-grained series can be built out of both positive- and implied
negative-frequency elements in the fine-grained series.
\f{equation}{
  h_0 = \frac{1}{\rho}
  \left[
    h'_0
    + 2\ \mathrm{Re}
    \left(
      \sum_{\ell=1}^{\ell^{\scriptscriptstyle\textrm{max}}_0}
      h'_\ell
      + (\lambda^{\scriptstyle\textrm{max}}_0 - \ell^{\scriptstyle\textrm{max}}_0)
      h'_{\ell^{\scriptscriptstyle\textrm{max}}_0+1}
    \right)
  \right]
\f}

\heading{Algorithm}

These routines move through the output series, using
\eqref{utilities_e_coarseapprox} to add up the contributions from the
bins in the fine-grained series.

\heading{Notes}
<ul>
<li> The coarse graining ratio must obey \f$\rho\ge 1\f$ (so the
  coarse-grained frequency spacing must be less than the fine-grained
  one).  Additionally, the bins in the fine-grained frequency series
  must \e completely overlap those in the coarse-grained frequency
  series.  In particular, since the lowest frequency in the first bin
  of the coarse-grained series is \f$f_{\scriptstyle{\rm min}}=f_0-\delta f/2\f$
  and the last is \f$f_{\scriptstyle{\rm max}}=f_0 + (N-1) \delta f +\delta f/2\f$
  (taking into account the width of the bins), the conitions are
  \f{eqnarray*}{
    f_0 - \frac{\delta f}{2} &\ge& f'_0 - \frac{\delta f'}{2}\\
    f_0 + \left(N-\frac{1}{2}\right)\,\delta f &\le&
     f'_0 + \left(N'-\frac{1}{2}\right)\,\delta f'
  \f}
  (The special case \f$f_0=f'_0=0\f$ is an
  exception to the condition on the minimum frequency.)</li>
<li> The routines return an error if either minimum frequency
  (\f$f_{\scriptstyle\textrm{min}}\f$ or \f$f'_{\scriptstyle\textrm{min}}\f$) is
  negative (unless \f$f_0=0\f$ or \f$f'_0=0\f$, respectively).</li>
</ul>

*/
/*@{*/

/**\name Error Codes */
/*@{*/
#define COARSEGRAINFREQUENCYSERIESH_ENULLPTR        1		/**< Null pointer */
#define COARSEGRAINFREQUENCYSERIESH_ESAMEPTR        2		/**< Input and Output pointers the same */
#define COARSEGRAINFREQUENCYSERIESH_EZEROLEN        3		/**< Zero length for data member of series */
#define COARSEGRAINFREQUENCYSERIESH_ENONPOSDELTAF   4		/**< Negative or zero frequency spacing */
#define COARSEGRAINFREQUENCYSERIESH_ENEGFMIN        5		/**< Negative start frequency */
#define COARSEGRAINFREQUENCYSERIESH_EMMHETERO       7		/**< Mismatch in heterodyning frequencies */
#define COARSEGRAINFREQUENCYSERIESH_EMMFMIN         8		/**< Mismatch in start frequencies */
#define COARSEGRAINFREQUENCYSERIESH_EMMDELTAF       9		/**< Mismatch in frequency spacings */
#define COARSEGRAINFREQUENCYSERIESH_EMMLEN         10		/**< Mismatch in sequence lengths */
#define COARSEGRAINFREQUENCYSERIESH_EOORCOARSE     16		/**< Coarse-graining paramaters out of range */
/*@}*/

/** \cond DONT_DOXYGEN */
#define COARSEGRAINFREQUENCYSERIESH_MSGENULLPTR    "Null pointer"
#define COARSEGRAINFREQUENCYSERIESH_MSGESAMEPTR    "Input and Output pointers the same"
#define COARSEGRAINFREQUENCYSERIESH_MSGEZEROLEN    "Zero length for data member of series"
#define COARSEGRAINFREQUENCYSERIESH_MSGENONPOSDELTAF "Negative or zero frequency spacing"
#define COARSEGRAINFREQUENCYSERIESH_MSGENEGFMIN "Negative start frequency"
#define COARSEGRAINFREQUENCYSERIESH_MSGEMMHETERO   "Mismatch in heterodyning frequencies"
#define COARSEGRAINFREQUENCYSERIESH_MSGEMMFMIN     "Mismatch in start frequencies"
#define COARSEGRAINFREQUENCYSERIESH_MSGEMMDELTAF   "Mismatch in frequency spacings"
#define COARSEGRAINFREQUENCYSERIESH_MSGEMMLEN      "Mismatch in sequence lengths"
#define COARSEGRAINFREQUENCYSERIESH_MSGEOORCOARSE  "Coarse-graining paramaters out of range"
/** \endcond */


  /* ***********************************************************
   *                                                           *
   *       Structures and prototypes associated with           *
   *             CoarseGrainFrequencySeries.c                  *
   *                                                           *
   *************************************************************/

/** Contains the parameters needed to specify the sampling of a frequency series */
typedef struct
tagFrequencySamplingParams
{
  REAL8         f0;		/**< The start frequency of the frequency series */
  REAL8         deltaF;		/**< The frequency spacing of the frequency series */
  UINT4         length;		/**< The number of points in the frequency series */
}
FrequencySamplingParams;

/** \see See \ref CoarseGrainFrequencySeries_h for documetation */
void
LALSCoarseGrainFrequencySeries(LALStatus                      *status,
			       REAL4FrequencySeries           *output,
			       const REAL4FrequencySeries     *input,
			       const FrequencySamplingParams  *params);

/** \see See \ref CoarseGrainFrequencySeries_h for documetation */
void
LALDCoarseGrainFrequencySeries(LALStatus                      *status,
			       REAL8FrequencySeries           *output,
			       const REAL8FrequencySeries     *input,
			       const FrequencySamplingParams  *params);

/** \see See \ref CoarseGrainFrequencySeries_h for documetation */
void
LALCCoarseGrainFrequencySeries(LALStatus                        *status,
			       COMPLEX8FrequencySeries          *output,
			       const COMPLEX8FrequencySeries    *input,
			       const FrequencySamplingParams    *params);

/** \see See \ref CoarseGrainFrequencySeries_h for documetation */
void
LALZCoarseGrainFrequencySeries(LALStatus                        *status,
			       COMPLEX16FrequencySeries          *output,
			       const COMPLEX16FrequencySeries    *input,
			       const FrequencySamplingParams    *params);

/*@}*/

#ifdef  __cplusplus
}
#endif /* C++ protection */


#endif /* _COARSEGRAINFREQUENCYSERIES_H */
