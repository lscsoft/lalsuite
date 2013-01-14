/*
*  Copyright (C) 2007 Teviet Creighton
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

#ifndef _SIMULATEINSPIRAL_H
#define _SIMULATEINSPIRAL_H

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/SkyCoordinates.h>
#include <lal/LALBarycenter.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/**
   \addtogroup SimulateInspiral_h
   \author Creighton, T. D.

   \brief Provides a routine to inject inspirals into time series data.

   \section synopsis Synopsis
   \code
   #include <lal/SimulateInspiral.h>
   \endcode

The routines in \ref GeneratePPNInspiral_h, \ref SimulateCoherentGW_h,
and \ref Inject_h provide a powerful
mechanism for simulating the instrumental response to a physical
inspiral event, including such considerations as the polarization
response and propagation delay for a particular instrument, which are
necessary if one wants to model coincident detection in a network of
detectors.  In many cases, though, one simply wants to generate a
signal with a given signal-to-noise ratio and a given coalescence
time, and inject it into a time series.  This header provides a
steamlined interface to accomplish this with a minimum of fuss to the
user.

In order to provide this streamlined interface, two calculations have
to be internalized.  First, the waveform must be time-shifted so that
it coalesces at the specified time.  This is straightforward and
requires no explanation.  Second, the waveform must be scaled to have
some specified amplitude.  To do this, we must first state what we
mean by "amplitude".

We define the <em>characteristic detection amplitude</em> \f$A_c\f$ of a
gravitational-wave signal to be its root summed squared contribution
to the sampled detector output.  That is, if the detector output can
be written as \f$o(t_k)=n(t_k)+s(t_k)\f$, where \f$n\f$ is the contribution
due to noise, \f$s\f$ is the contribution due to signal, and \f$t_k=k\Delta
t\f$ are the discrete time samples, then:
\anchor eq_SimulateInspiralH_characteristic_amplitude \f{equation}{
\label{eq_SimulateInspiralH_characteristic_amplitude}
A_c \equiv \sqrt{\sum_{k=-\infty}^\infty |s(t_k)|^2} \;.
\f}
If \f$T(f)\f$ is the detector transfer function (such that a gravitational
have signal \f$\tilde{h}(f)\f$ in the frequency domain produces an output
\f$\tilde{o}(f)=\tilde{n}(f)+T(f)\tilde{h}(f)\f$), the characteristic
detection amplitude has the not-so-obvious relation that:
\anchor eq_SimulateInspiralH_characteristic_gw_amplitude \f{equation}{
\label{eq_SimulateInspiralH_characteristic_gw_amplitude}
A_c^2 = \int_{-\infty}^\infty \frac{df}{\Delta t}
	|T(f)\tilde{h}(f)|^2 \;.
\f}
So why use this quantity to specify the signal amplitude?  First, it
is easy for a simulation/injection routine to calculate.  Second, if
we assume that \f$T(f)\f$ is a true whitening filter such that the output
power spectral density is flat \f$S_o(f)=S_o=\f$constant, then the sampled
noise output is uncorrelated noise with a mean of zero and a variance
of \f$\sigma_n^2=S_o/2\Delta t\f$.  The intrinsic signal-to-noise power is
then given by the simple relation:
\anchor eq_SimulateInspiralH_instrinsic_snr_power \f{eqnarray}{
(h|h) & = & 2\int_{-\infty}^\infty df\,\frac{|\tilde{h}(f)|^2}{S_h(f)}
	\nonumber\\
      & = & 2\int_{-\infty}^\infty df\,\frac{|T(f)\tilde{h}(f)|^2}{S_o}
	\nonumber\\
      & = & \frac{A_c^2}{\sigma_n^2} \;.
	\label{eq_SimulateInspiralH_instrinsic_snr_power}
\f}
Thus to simulate a signal with an intrinsic signal-to-noise amplitude
\f$\sqrt{(h|h)}\f$, simply fill a data vector with uncorrelated noise with
variance \f$\sigma_n^2\f$ and inject a signal with a characteristic
detection amplitude of \f$A_c=\sigma_n\sqrt{(h|h)}\f$.

We refer the reader to the Conventions section for a more detailed derivation of these
specifications.

\section sec_si_conv Conventions

We define here the conventions we use when talking about
signal-to-noise ratios in coloured and white noise.  You may also want
to read the signal processing conventions in Secs.\ \ref ss_conventions and \ref ss_psdconv of
\ref pkg_findchirp, since this section is esentially a summary and extension of those conventions.

\subsection sec_si_SNR Signal-to-noise definitions

We first reiterate the standard definitions (given in the
\ref pkg_findchirp) of the Fourier transform pair:
\anchor eq_SimulateInspiralH_fourier_transforms \f{equation}{
\label{eq_SimulateInspiralH_fourier_transforms}
\tilde{a}(f) = \int_{-\infty}^\infty dt\,a(t)e^{-2\pi ift}
\qquad\rightleftharpoons\qquad
a(t) = \int_{-\infty}^\infty df\,\tilde{a}(f)e^{2\pi ift} \;,
\f}
of the power spectral density \f$S_n(f)\f$ of a stationary random process
(noise) \f$n(t)\f$:
\anchor eq_SimulateInspiralH_noise_psd \f{equation}{
\label{eq_SimulateInspiralH_noise_psd}
\langle\tilde{n}(f)\tilde{n}^*(f')\rangle
	= \frac{1}{2}S_n(f)\delta(f-f')
\f}
(where \f$\langle\ldots\rangle\f$ denotes an enseble average over
instantiations of the random process), and finally of the
noise-weighted inner product \f$(a|b)\f$ of two time series \f$a(t)\f$ and
\f$b(t)\f$:
\anchor eq_SimulateInspiralH_inner_product \f{equation}{
\label{eq_SimulateInspiralH_inner_product}
(a|b) = \int_{-\infty}^\infty df\,\frac{\tilde{a}(f)\tilde{b}(f)^*
	+ \tilde{a}(f)^*\tilde{b}(f)}{S_n(f)} \;.
\f}
In the case where the time series are all real and the noise has zero
mean \f$\langle n(t)\rangle=0\f$, the weighting on the inner product leads
to the property that:
\anchor eq_SimulateInspiralH_inner_product_normalisation \f{equation}{
\label{eq_SimulateInspiralH_inner_product_normalisation}
\langle(n|s)^2\rangle = (s|s) \;.
\f}
We call this quantity the <em>intrinsic signal-to-noise power</em> of
the waveform \f$s(t)\f$, and its square root \f$\sqrt{(s|s)}\f$ the
<em>intrinsic signal-to-noise amplitude</em>.

In the theory of signal processing, if one has a data stream
\f$o(t)=n(t)+s(t)\f$ and can determine a <em>matched filter</em>
\f$a(t)\propto s(t)\f$, one can then define a signal-to-noise estimator
\f$r=(o|a)/\sqrt{(a|a)}\f$ with the property that \f$\langle
r\rangle=\sqrt{(s|s)}\f$ and \f$\sigma_r=\sqrt{\langle r^2\rangle-\langle
r\rangle^2}=1\f$ (these follow from
Eq.\eqref{eq_SimulateInspiralH_inner_product_normalisation} and
\f$\langle n\rangle=0\f$).  This is the justification for calling
\f$\sqrt{(s|s)}\f$ a signal-to-noise ratio.

However, in many cases one can only define a set of orthogonal
waveforms \f$\{a_1,\ldots,a_N:(a_i|a_j)=(a_i|a_i)\delta_{ij}\}\f$ spanning
the space of possible signals.  Specifically, for inspiral signals,
one does not know in advance the phase of the waveform, and so one
must filter the data using two orthogonal waveforms \f$\{a_s,a_c\}\f$ (the
sine and cosine quadratures) that are \f$90^\circ\f$ out of phase.  As
described in the
\ref FindChirp_h, the optimal statistic in this case is:
\anchor eq_SimulateInspiralH_rhosq \f{equation}{
\label{eq_SimulateInspiralH_rhosq}
\rho^2 = \frac{(o|a_s)^2}{(a_s|a_s)} + \frac{(o|a_c)^2}{(a_c|a_c)}
\f}
(where we have implicitly already maximised over any filter
parameters, including time-of-arrival).  This statistic no longer has
unit variance, but instead has the property that:
\anchor eq_SimulateInspiralH_rhosq_expectation \f{equation}{
\label{eq_SimulateInspiralH_rhosq_expectation}
\langle\rho^2\rangle = 2 + (s|s)\;.
\f}
By comparison, for the perfectly-matched filter one has \f$\langle
r^2\rangle=1+(s|s)\f$, which is why one often says that the search over
phase halves the effective signal-to-noise power.  However, this
statement is ambiguous, as the two statistics \f$\rho\f$ and \f$r\f$ have
completely diffferent probability distributions in the presence of
noise, and even in the absence of any signal we have
\f$\langle\rho\rangle>0\f$ (its value depends on the particular noise
probability distribution).

\subsection sec_si_wn Specification to white noise

White noise is stationary zero-mean noise whose power spectral density
\f$S_n(f)\f$ is independent of \f$f\f$.  This property allows inner products
to be expressed equivalently in the time domain as well as in the
frequency domain:
\anchor eq_SimulateInspiralH_inner_product_time \f{equation}{
\label{eq_SimulateInspiralH_inner_product_time}
(a|b) = \int_{-\infty}^\infty dt\,\frac{a(t)b(t)^* + a(t)^*b(t)}{S_n} \;.
\f}
If the white noise process \f$n(t)\f$ is discretely sampled at intervals
\f$\Delta t\f$, we find that different time samples are uncorrelated:
\anchor eq_SimulateInspiralH_noise_correlation \f{equation}{
\label{eq_SimulateInspiralH_noise_correlation}
\langle n(t_j)n(t_{j'})^*\rangle = \sigma_n^2 \delta_{jj'} \;,
\f}
where \f$\sigma_n^2\f$ is the variance in the sampled noise.  This can be
proven using the definitions of the discrete inverse FFT and discrete
power spectrum given in Eqs.\ (.9) and\ (.19) of the \c findchirp
package:
\f{eqnarray}{
\langle n(t_j)n(t_{j'})^*\rangle
& = & \frac{1}{N^2}\sum_{k=0}^{N-1}\sum_{k'=0}^{N-1}
	e^{2\pi i(jk-j'k')/N} \langle\tilde{n}_k\tilde{n}_{k'}\rangle
	\nonumber\\
& = & \frac{1}{N^2}\sum_{k=0}^{N-1}\sum_{k'=0}^{N-1}
	e^{2\pi i(jk-j'k')/N} \frac{N}{2\Delta t}S_n\delta_{kk'}
	\nonumber\\
& = & \frac{S_n}{2N\Delta t}\sum_{k=0}^{N-1} e^{2\pi i(j-j')k/N}
	\nonumber\\
& = & \frac{S_n}{2\Delta t}\delta_{jj'} \;,\nonumber
\f}
whence we have \f$\sigma_n^2=S_n/2\Delta t\f$.  We note that the
divergence as \f$\Delta t\rightarrow0\f$ reflects the fact that the
idealized white power spectral density integrates to an infinite
amount of power as \f$f\rightarrow\infty\f$.  Real noise processes always
have some minimum timescale below which correlations between
subsequent measurements are inevitable.

Reexpressing the inner product in
Eq.\eqref{eq_SimulateInspiralH_inner_product_time} for
discretely-sampled time series, we have:
\anchor eq_SimulateInspiralH_inner_product_sampled \f{equation}{
\label{eq_SimulateInspiralH_inner_product_sampled}
(a|b) = \sum_k \Delta t\,\frac{a(t_k)b(t_k)^* + a(t_k)^*b(t_k)}{S_n}
      = \frac{\mathrm{Re}\left[\sum_k
		a(t_k)b(t_k)^*\right]}{\sigma_n^2} \;.
\f}
The intrinsic signal-to-noise amplitude also takes an intuitive form:
\anchor eq_SimulateInspiralH_snr_sampled \f{equation}{
\label{eq_SimulateInspiralH_snr_sampled}
\sqrt{(s|s)} = \frac{\sqrt{\sum_k |s(t_k)|^2}}{\sigma_n}
\f}
This, then, is the key formula that we will use to determine the
signal-to-noise amplitude of a whitened signal that is to be injected
into white noise.  It applies to any stationary white noise (recall
that "white" implies zero mean), but we will almost always be
concerned with white Gaussian noise, where the differential
probability of a noise sample \f$n(t_k)\f$ lying in an infinitesimal range
\f$(n,n+dn)\f$ is:
\anchor eq_SimulateInspiralH_gaussian_pdf \f{equation}{
\label{eq_SimulateInspiralH_gaussian_pdf}
\frac{dP[n(t_k)\in(n,n+dn)]}{dn} = \frac{1}{\sqrt{2\pi\sigma_n^2}}
	e^{-n^2/2\sigma_n^2} \;.
\f}

*/
/*@{*/

/** \name Error Codes */
/*@{*/
#define SIMULATEINSPIRALH_ENUL 1	/**< Unexpected null pointer in arguments */
#define SIMULATEINSPIRALH_EMEM 2	/**< Memory allocation error */
#define SIMULATEINSPIRALH_EDF  3	/**< Transfer frequency interval is zero */
#define SIMULATEINSPIRALH_EBAD 4	/**< Bad parameters: ac and dEff are negative */
/*@}*/

/** \cond DONT_DOXYGEN */
#define SIMULATEINSPIRALH_MSGENUL "Unexpected null pointer in arguments"
#define SIMULATEINSPIRALH_MSGEMEM "Memory allocation error"
#define SIMULATEINSPIRALH_MSGEDF  "Transfer frequency interval is zero"
#define SIMULATEINSPIRALH_MSGEBAD "Bad parameters: ac and dEff are negative"
/** \endcond */

/** This structure stores the parameters required to simulate a
 * set of inspiral signal in white Gaussian noise.  It can be part of a
 * linked list of inspiral events, to allow for multiple injections.
 */
typedef struct tagSimulateInspiralParamStruc {
  LIGOTimeGPS timeC;      /**< time of coalescence */
  REAL4 phiC;             /**< The wave phase at coalescence, in radians */
  REAL4 mass1, mass2;     /**< The masses of the binary components, in \f$M_\odot\f$ */

  /** \name Amplitude.
   * The characteristic detection amplitude \f$A_c\f$, in ADC counts, and the effective distance
   * in Mpc of an optimally-oriented source that would give that amplitude.
   * A negative number means the quantity is unspecified.  In general only
   * one of these must be specified by the user; the simulation routine
   * will set the other to be consistent with the first
   */
  /*@{*/
  REAL4 signalAmplitude;
  REAL4 effDist;
  /*@}*/
  REAL4 fStart;           /**< The lower cutoff frequency at which
                           * waveform generation will begin, in Hz;  If \f$\leq0\f$, the cutoff
                           * frequency will be taken as the point where the instrument sensitivity
                           * function is \f$\sim10^{-6}\f$ of its optimal value, as determined from the
                           * transfer function
                           */
  struct tagSimulateInspiralParamStruc *next; /**< Pointer to another inspiral event to be injected, or \c NULL if this is the last (or only) injection */
} SimulateInspiralParamStruc;


/* ---------- Function prototypes ---------- */

void
LALSimulateInspiral( LALStatus                  *,
		     REAL4TimeSeries            *output,
		     COMPLEX8FrequencySeries    *transfer,
		     SimulateInspiralParamStruc *params );

/*@}*/ /* end:SimulateInspiral_h */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _SIMULATEINSPIRAL_H */
