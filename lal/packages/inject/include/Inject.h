/*************************************** <lalVerbatim file="InjectHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{Inject.h}}
\label{s:Inject.h}

Provides routines to inject a signal into detector output.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Inject.h>
\end{verbatim}

This header provides simple routines to inject a signal, stored as a
floating-point time series, into an integer time series that
represents the ADC output of a detector channel.

The basic concept at work here is that of \emph{dithering}.  That is,
to add a real signal $x(t_k)$ to the integer output $n(t_k)$ of an
ADC, we cannot simply compute $x+n$ at each time $t_k$ and round to
the nearest integer.  To see this, consider injecting a sinusoid with
an amplitude less than half the ADC resolution, $|x|<0.5$.  Then
adding $x+n$ and rounding will always give back $n$, and the signal
injection will have no effect on the output.

Instead, what we would like to do is to add $x$ to the ADC
\emph{input} before rounding occurs, and then round to an integer.
That is, the ADC input was actually $n+d$, where $d\in[-0.5,0.5)$; we
want to compute the rounded value of $n+d+x$, which may occasionally
be different from $n$.  In principle, a Fourier transforms can detect
a sinusoidal signal with an output much smaller than the ADC
resolution: by integrating enough data, one can eventually detect a
statistically significant phase correlation between the occasional
increments and decrements in $n$.

Of course given the output time series $n(t_k)$ we can only guess at
the input series $n(t_k)+d(t_k)$ that produced it.  The simplest guess
is to assume that each $d(t_k)$ is an independent random variable with
a flat distribution over the range $[-0.5,0.5)$.  This is a reasonable
guess to make if the root mean square variation between succesive
output values $n$ is a few or more ADC counts; i.e. if the
dimensionless power spectral density $\sqrt{fS(f)}$ has a value of a
few or more around the sampling freqeuncy $f$.  This is almost always
true of any detector designed to work at or near its noise limit: the
input to the ADC will first be whitened so that $\sqrt{fS(f)}$ is
nearly flat, and then amplified so that $\sqrt{fS(f)}$ is on the order
of several (or more) ADC counts.

In the routines covered by this header we will take it for granted
that the above is a reasonable approximation, and will not check for
it.  We will further assume that the signal to be injected has already
been subjected to the same whitening and amplification, so that the
units of $x(t_k)$ are normalized ADC counts (although it is still a
real number, not an integer).

******************************************************* </lalLaTeX> */

#ifndef _INJECT_H
#define _INJECT_H

#include <lal/LALStdlib.h>
#include <lal/Random.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( INJECTH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define INJECTH_ENUL 1
#define INJECTH_EBAD 2

#define INJECTH_MSGENUL "Unexpected null pointer in arguments"
#define INJECTH_MSGEBAD "A sampling interval is (effectively) zero"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Structures}
******************************************************* </lalLaTeX> */

/* <lalLaTeX>
\vfill{\footnotesize\input{InjectHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{InjectVectorC}
</lalLaTeX> */
void
LALSI2InjectVector( LALStatus    *stat,
		    INT2Vector   *output,
		    REAL4Vector  *signal,
		    RandomParams *params );

/* <lalLaTeX>
\newpage\input{InjectTimeSeriesC}
</lalLaTeX> */
void
LALSI2InjectTimeSeries( LALStatus       *stat,
			INT2TimeSeries  *output,
			REAL4TimeSeries *signal,
			RandomParams    *params );

/* <lalLaTeX>
\newpage\input{BasicInjectTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _INJECT_H */
