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

/************************ <lalVerbatim file="SimulateInspiralHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{SimulateInspiral.h}}
\label{s:SimulateInspiral.h}

Provides a routine to inject inspirals into time series data.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/SimulateInspiral.h>
\end{verbatim}

The routines in \verb@GeneratePPNInspiral.h@,
\verb@SimulateCoherentGW.h@, and \verb@Inject.h@ provide a powerful
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
mean by ``amplitude''.

We define the \emph{characteristic detection amplitude} $A_c$ of a
gravitational-wave signal to be its root summed squared contribution
to the sampled detector output.  That is, if the detector output can
be written as $o(t_k)=n(t_k)+s(t_k)$, where $n$ is the contribution
due to noise, $s$ is the contribution due to signal, and $t_k=k\Delta
t$ are the discrete time samples, then:
\begin{equation}
\label{eq:SimulateInspiralH:characteristic-amplitude}
A_c \equiv \sqrt{\sum_{k=-\infty}^\infty |s(t_k)|^2} \;.
\end{equation}
If $T(f)$ is the detector transfer function (such that a gravitational
have signal $\tilde{h}(f)$ in the frequency domain produces an output
$\tilde{o}(f)=\tilde{n}(f)+T(f)\tilde{h}(f)$), the characteristic
detection amplitude has the not-so-obvious relation that:
\begin{equation}
\label{eq:SimulateInspiralH:characteristic-gw-amplitude}
A_c^2 = \int_{-\infty}^\infty \frac{df}{\Delta t}
	|T(f)\tilde{h}(f)|^2 \;.
\end{equation}
So why use this quantity to specify the signal amplitude?  First, it
is easy for a simulation/injection routine to calculate.  Second, if
we assume that $T(f)$ is a true whitening filter such that the output
power spectral density is flat $S_o(f)=S_o=$constant, then the sampled
noise output is uncorrelated noise with a mean of zero and a variance
of $\sigma_n^2=S_o/2\Delta t$.  The intrinsic signal-to-noise power is
then given by the simple relation:
\begin{eqnarray}
(h|h) & = & 2\int_{-\infty}^\infty df\,\frac{|\tilde{h}(f)|^2}{S_h(f)}
	\nonumber\\
      & = & 2\int_{-\infty}^\infty df\,\frac{|T(f)\tilde{h}(f)|^2}{S_o}
	\nonumber\\
      & = & \frac{A_c^2}{\sigma_n^2} \;.
	\label{eq:SimulateInspiralH:instrinsic-snr-power}
\end{eqnarray}
Thus to simulate a signal with an intrinsic signal-to-noise amplitude
$\sqrt{(h|h)}$, simply fill a data vector with uncorrelated noise with
variance $\sigma_n^2$ and inject a signal with a characteristic
detection amplitude of $A_c=\sigma_n\sqrt{(h|h)}$.

We refer the reader to the Conventions section at the end of this
header documentation for a more detailed derivation of these
specifications.

******************************************************* </lalLaTeX> */

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

NRCSID( SIMULATEINSPIRALH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define SIMULATEINSPIRALH_ENUL 1
#define SIMULATEINSPIRALH_EMEM 2
#define SIMULATEINSPIRALH_EDF  3
#define SIMULATEINSPIRALH_EBAD 4

#define SIMULATEINSPIRALH_MSGENUL "Unexpected null pointer in arguments"
#define SIMULATEINSPIRALH_MSGEMEM "Memory allocation error"
#define SIMULATEINSPIRALH_MSGEDF  "Transfer frequency interval is zero"
#define SIMULATEINSPIRALH_MSGEBAD "Bad parameters: ac and dEff are negative"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{SimulateInspiralParamStruc}}
\idx[Type]{SimulateInspiralParamStruc}

\noindent This structure stores the parameters required to simulate a
set of inspiral signal in white Gaussian noise.  It can be part of a
linked list of inspiral events, to allow for multiple injections.  It
consists of the following fields:

\begin{description}
\item[\texttt{LIGOTimeGPS timeC}] The time of coalescence.

\item[\texttt{REAL4 phiC}] The wave phase at coalescence, in radians.

\item[\texttt{REAL4 mass1, mass2}] The masses of the binary
components, in $M_\odot$.

\item[\texttt{REAL4 signalAmplitude, effDist}] The characteristic
detection amplitude $A_c$, in ADC counts, and the effective distance
in Mpc of an optimally-oriented source that would give that amplitude.
A negative number means the quantity is unspecified.  In general only
one of these must be specified by the user; the simulation routine
will set the other to be consistent with the first.

\item[\texttt{REAL4 fStart}] The lower cutoff frequency at which
waveform generation will begin, in Hz.  If $\leq0$, the cutoff
frequency will be taken as the point where the instrument sensitivity
function is $\sim10^{-6}$ of its optimal value, as determined from the
transfer function.

\item[\texttt{SimulateInspiralParamStruc *next}] Pointer to another
inspiral event to be injected, or \verb@NULL@ if this is the last (or
only) injection.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagSimulateInspiralParamStruc {
  LIGOTimeGPS timeC;      /* time of coalescence */
  REAL4 phiC;             /* phase at coalescence */
  REAL4 mass1, mass2;     /* binary masses (solar masses) */
  REAL4 signalAmplitude;  /* characteristic amplitude (counts) */
  REAL4 effDist;          /* effective distance (Mpc) */
  REAL4 fStart;           /* waveform start frequency (Hz) */
  struct tagSimulateInspiralParamStruc *next; /* next node in list */
} SimulateInspiralParamStruc;

/********************************************************** <lalLaTeX>

\subsection*{Conventions}

We define here the conventions we use when talking about
signal-to-noise ratios in coloured and white noise.  You may also want
to read the signal processing conventions in Secs.~.1.1 and~.1.2 of
the \verb@findchirp@ package, since this section is esentially a
summary and extension of those conventions.

\subsubsection*{Signal-to-noise definitions}

We first reiterate the standard definitions (given in the
\verb@findchirp@ package) of the Fourier transform pair:
\begin{equation}
\label{eq:SimulateInspiralH:fourier-transforms}
\tilde{a}(f) = \int_{-\infty}^\infty dt\,a(t)e^{-2\pi ift}
\qquad\rightleftharpoons\qquad
a(t) = \int_{-\infty}^\infty df\,\tilde{a}(f)e^{2\pi ift} \;,
\end{equation}
of the power spectral density $S_n(f)$ of a stationary random process
(noise) $n(t)$:
\begin{equation}
\label{eq:SimulateInspiralH:noise-psd}
\langle\tilde{n}(f)\tilde{n}^*(f')\rangle
	= \frac{1}{2}S_n(f)\delta(f-f')
\end{equation}
(where $\langle\ldots\rangle$ denotes an enseble average over
instantiations of the random process), and finally of the
noise-weighted inner product $(a|b)$ of two time series $a(t)$ and
$b(t)$:
\begin{equation}
\label{eq:SimulateInspiralH:inner-product}
(a|b) = \int_{-\infty}^\infty df\,\frac{\tilde{a}(f)\tilde{b}(f)^*
	+ \tilde{a}(f)^*\tilde{b}(f)}{S_n(f)} \;.
\end{equation}
In the case where the time series are all real and the noise has zero
mean $\langle n(t)\rangle=0$, the weighting on the inner product leads
to the property that:
\begin{equation}
\label{eq:SimulateInspiralH:inner-product-normalisation}
\langle(n|s)^2\rangle = (s|s) \;.
\end{equation}
We call this quantity the \emph{intrinsic signal-to-noise power} of
the waveform $s(t)$, and its square root $\sqrt{(s|s)}$ the
\emph{intrinsic signal-to-noise amplitude}.

In the theory of signal processing, if one has a data stream
$o(t)=n(t)+s(t)$ and can determine a \emph{matched filter}
$a(t)\propto s(t)$, one can then define a signal-to-noise estimator
$r=(o|a)/\sqrt{(a|a)}$ with the property that $\langle
r\rangle=\sqrt{(s|s)}$ and $\sigma_r=\sqrt{\langle r^2\rangle-\langle
r\rangle^2}=1$ (these follow from
Eq.~(\ref{eq:SimulateInspiralH:inner-product-normalisation}) and
$\langle n\rangle=0$).  This is the justification for calling
$\sqrt{(s|s)}$ a signal-to-noise ratio.

However, in many cases one can only define a set of orthogonal
waveforms $\{a_1,\ldots,a_N:(a_i|a_j)=(a_i|a_i)\delta_{ij}\}$ spanning
the space of possible signals.  Specifically, for inspiral signals,
one does not know in advance the phase of the waveform, and so one
must filter the data using two orthogonal waveforms $\{a_s,a_c\}$ (the
sine and cosine quadratures) that are $90^\circ$ out of phase.  As
described in the
\verb@FindChirp.h@ header, the optimal statistic in this case is:
\begin{equation}
\label{eq:SimulateInspiralH:rhosq}
\rho^2 = \frac{(o|a_s)^2}{(a_s|a_s)} + \frac{(o|a_c)^2}{(a_c|a_c)}
\end{equation}
(where we have implicitly already maximised over any filter
parameters, including time-of-arrival).  This statistic no longer has
unit variance, but instead has the property that:
\begin{equation}
\label{eq:SimulateInspiralH:rhosq-expectation}
\langle\rho^2\rangle = 2 + (s|s)\;.
\end{equation}
By comparison, for the perfectly-matched filter one has $\langle
r^2\rangle=1+(s|s)$, which is why one often says that the search over
phase halves the effective signal-to-noise power.  However, this
statement is ambiguous, as the two statistics $\rho$ and $r$ have
completely diffferent probability distributions in the presence of
noise, and even in the absence of any signal we have
$\langle\rho\rangle>0$ (its value depends on the particular noise
probability distribution).

\subsubsection*{Specification to white noise}

White noise is stationary zero-mean noise whose power spectral density
$S_n(f)$ is independent of $f$.  This property allows inner products
to be expressed equivalently in the time domain as well as in the
frequency domain:
\begin{equation}
\label{eq:SimulateInspiralH:inner-product-time}
(a|b) = \int_{-\infty}^\infty dt\,\frac{a(t)b(t)^* + a(t)^*b(t)}{S_n} \;.
\end{equation}
If the white noise process $n(t)$ is discretely sampled at intervals
$\Delta t$, we find that different time samples are uncorrelated:
\begin{equation}
\label{eq:SimulateInspiralH:noise-correlation}
\langle n(t_j)n(t_{j'})^*\rangle = \sigma_n^2 \delta_{jj'} \;,
\end{equation}
where $\sigma_n^2$ is the variance in the sampled noise.  This can be
proven using the definitions of the discrete inverse FFT and discrete
power spectrum given in Eqs.~(.9) and~(.19) of the \verb@findchirp@
package:
\begin{eqnarray}
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
\end{eqnarray}
whence we have $\sigma_n^2=S_n/2\Delta t$.  We note that the
divergence as $\Delta t\rightarrow0$ reflects the fact that the
idealized white power spectral density integrates to an infinite
amount of power as $f\rightarrow\infty$.  Real noise processes always
have some minimum timescale below which correlations between
subsequent measurements are inevitable.

Reexpressing the inner product in
Eq.~(\ref{eq:SimulateInspiralH:inner-product-time}) for
discretely-sampled time series, we have:
\begin{equation}
\label{eq:SimulateInspiralH:inner-product-sampled}
(a|b) = \sum_k \Delta t\,\frac{a(t_k)b(t_k)^* + a(t_k)^*b(t_k)}{S_n}
      = \frac{\mathrm{Re}\left[\sum_k
		a(t_k)b(t_k)^*\right]}{\sigma_n^2} \;.
\end{equation}
The intrinsic signal-to-noise amplitude also takes an intuitive form:
\begin{equation}
\label{eq:SimulateInspiralH:snr-sampled}
\sqrt{(s|s)} = \frac{\sqrt{\sum_k |s(t_k)|^2}}{\sigma_n}
\end{equation}
This, then, is the key formula that we will use to determine the
signal-to-noise amplitude of a whitened signal that is to be injected
into white noise.  It applies to any stationary white noise (recall
that ``white'' implies zero mean), but we will almost always be
concerned with white Gaussian noise, where the differential
probability of a noise sample $n(t_k)$ lying in an infinitesimal range
$(n,n+dn)$ is:
\begin{equation}
\label{eq:SimulateInspiralH:gaussian-pdf}
\frac{dP[n(t_k)\in(n,n+dn)]}{dn} = \frac{1}{\sqrt{2\pi\sigma_n^2}}
	e^{-n^2/2\sigma_n^2} \;.
\end{equation}

******************************************************* </lalLaTeX> */

/* <lalLaTeX>
\vfill{\footnotesize\input{SimulateInspiralHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{SimulateInspiralC}
</lalLaTeX> */
void
LALSimulateInspiral( LALStatus                  *,
		     REAL4TimeSeries            *output,
		     COMPLEX8FrequencySeries    *transfer,
		     SimulateInspiralParamStruc *params );

/* <lalLaTeX>
%\newpage\input{SimulateInspiralTestC}
</lalLaTeX> */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _SIMULATEINSPIRAL_H */
