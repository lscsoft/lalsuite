/********************************** <lalVerbatim file="StackMetricHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\newcommand{\bm}[1]{\mbox{\boldmath$#1$\unboldmath}}

\section{Header \texttt{StackMetric.h}}
\label{s:StackMetric.h}

Provides routines to compute parameter-space metrics for coherent or
stacked pulsar searches.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/StackMetric.h>
\end{verbatim}

\noindent This header covers routines that determine the metric
coefficients for the mismatch function (ambiguity function) on the
parameter space for a pulsar search.  The assumed search method is
stacking one or more Fourier power spectra, after some suitable
demodulation.

The method for determining the parameter metric is discussed in detail
in Sec.~II of~\cite{Brady_P:2000}; we present the key results here in
brief.  We assume that a model waveform in our search is described by
an overall frequency scale $f_0$, and by some modulation about that
frequency described by ``shape'' parameters
$\vec\lambda=(\lambda^1,\ldots,\lambda^n)$, such that the
parameterized phase of the waveform is $\phi[t;\bm{\lambda}] = 2\pi
f_0\tau[t;\vec\lambda]$.  Here $\bm{\lambda} = (\lambda^0,\vec\lambda)
= (f_0,\lambda^1,\ldots,\lambda^n)$ represents the total set of
parameters we must search over, and $\tau[t;\vec\lambda]$ is a
canonical time coordinate describing the shape of the waveform.

A (local) maximum in detected power $P$ occurs if a signal is filtered
(or demodulated in time and the Fourier spectrum is sampled) using a
phase model that matches the true phase of the signal.  If the
parameters $\bm{\lambda}$ do not match the true parameters of the
signal, then the detected power will be degraded.  The fractional
power loss $\Delta P/P$ thus has a (local) minimum of 0 for matched
parameters and increases for mismatched parameters; it can be thought
of as describing a distance between the two (nearby) parameter sets.
The \emph{metric} of this distance measure is simply the set of
quadratic coefficients of the Taylor expansion of $\Delta P/P$ about
its minimum.

Clearly the power will degrade rapidly with variation in some
parameter $\lambda^\alpha$ if the phase function $\phi$ depends
strongly on that parameter.  It turns out that if the detected power
is computed from a coherent power spectrum of a time interval $\Delta
t$, then the metric components are given simply by the covariances of
the phase derivatives $\partial\phi/\partial\lambda^\alpha$ over the
time interval:
\begin{equation}
g_{\alpha\beta}(\bm\lambda) =
	\left\langle
	\frac{\partial\phi[t;\bm{\lambda}]}{\partial\lambda^\alpha}
	\frac{\partial\phi[t;\bm{\lambda}]}{\partial\lambda^\beta}
	\right\rangle
	-
	\left\langle
	\frac{\partial\phi[t;\bm{\lambda}]}{\partial\lambda^\alpha}
	\right\rangle
	\left\langle
	\frac{\partial\phi[t;\bm{\lambda}]}{\partial\lambda^\beta}
	\right\rangle \; ,
\label{eq:gab-phi}
\end{equation}
where $\langle\ldots\rangle$ denotes a time average over the interval
$\Delta t$, and $\alpha$ and $\beta$ are indecies running from 0 to
$n$.  The partial derivatives are evaluated at the point \bm{\lambda}
in parameter space.  If instead the detected power is computed from
the sum of several power spectra computed from separate time intervals
(of the same length), then the overall metric is the \emph{average} of
the metrics from each time interval.

When power spectra are computed using fast Fourier transforms, the
entire frequency band from DC to Nyquist is computed at once; one then
scans all frequencies for significant peaks.  In this case one is
concerned with how the peak power (maximized over frequency) is
reduced by mismatch in the remaining ``shape'' parameters
$\vec\lambda$.  This is given by the the \emph{projected} metric
$\gamma_{ij}(\vec\lambda)$, where $i$ and $j$ run from 1 to $n$:
\begin{equation}
\gamma_{ij}(\vec\lambda) = \left[g_{ij}-\frac{g_{0i}g_{0j}}{g_{00}}
	\right]_{\lambda^0=f_\mathrm{max}} \; .
\label{eq:gij-gab}
\end{equation}
Here $f_\mathrm{max}$ is the highest-frequency signal expected to be
present, which ensures that, for lower-frequency signals,
$\gamma_{ij}$ will \emph{overestimate} the detection scheme's
sensitivity to the ``shape'' parameters.

******************************************************* </lalLaTeX> */

#ifndef _STACKMETRIC_H
#define _STACKMETRIC_H

#include <lal/LALStdlib.h>
#include <lal/PulsarTimes.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID(STACKMETRICH,"$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define STACKMETRICH_ENUL 1
#define STACKMETRICH_EBAD 2

#define STACKMETRICH_MSGENUL "Null pointer"
#define STACKMETRICH_MSGEBAD "Bad parameter values"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Structures}
\begin{verbatim}
MetricParamStruc
\end{verbatim}
\index{\texttt{MetricParamStruc}}

\noindent This structure stores and passes parameters for computing a
parameter-space metric.  It points to the canonical time function used
to compute the metric and to the parameters required by this function.
In addition, this structure must indicate the timespan over which the
timing differences accumulate, and whether this accumulation is
coherent or divided into stacks which are summed in power.  The fields
are:

\begin{description}
\item[\texttt{void *dtCanon( LALStatus *, REAL8Vector *, REAL8Vector
*, PulsarTimesParamStruc * )}] The function to compute the canonical
time coordinate and its derivatives.

\item[\texttt{PulsarTimesParamStruc *constants}] The constant
parameters used by \verb@*dt()@.

\item[\texttt{REAL8 start}] Start time of search, measured relative to
\verb@constants->epoch@.

\item[\texttt{REAL8 deltaT}] Length of each stack, in s.

\item[\texttt{UINT4 n}] Number of stacks.

\item[\texttt{BOOLEAN errors}] Whether to estimate errors in the
metric components.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagMetricParamStruc{
  void (*dtCanon)( LALStatus *, REAL8Vector *, REAL8Vector *,
		   PulsarTimesParamStruc * );
  /* The function to compute the canonical time coordinate and its
     derivatives. */
  PulsarTimesParamStruc *constants;
  /* The constant parameters used by *dt(). */
  REAL8 start; 
  /* Start time of search, measured relative to constants->epoch. */
  REAL8 deltaT;
  /* Length of each stack, in s. */
  UINT4 n;
  /* Number of stacks. */
  BOOLEAN errors;
  /* Whether to estimate errors in metric. */
} MetricParamStruc;


/* <lalLaTeX>
\vfill{\footnotesize\input{StackMetricHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{CoherentMetricC}
</lalLaTeX> */
void
LALCoherentMetric( LALStatus *stat,
		   REAL8Vector *metric,
		   REAL8Vector *lambda,
		   MetricParamStruc *params );

/* <lalLaTeX>
\newpage\input{StackMetricC}
</lalLaTeX> */
void
LALStackMetric( LALStatus *stat,
		REAL8Vector *metric,
		REAL8Vector *lambda,
		MetricParamStruc *params );

/* <lalLaTeX>
\newpage\input{ProjectMetricC}
</lalLaTeX> */
void
LALProjectMetric( LALStatus *stat, REAL8Vector *metric, BOOLEAN errors );

/* <lalLaTeX>
\newpage\input{StackMetricTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _STACKMETRIC_H */
