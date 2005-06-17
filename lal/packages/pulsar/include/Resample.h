/************************************* <lalVerbatim file="ResampleHV">
Author: Creighton, T. D.
Revision: $Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{Resample.h}}
\label{s:Resample.h}

Provides routines for resampling time series according to a new
canonical time coordinate.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Resample.h>
\end{verbatim}

\noindent One of the crucial problems in searching for
constant-frequency astrophysical signals is removing the effects of
Doppler modulation due to the Earth's motion.  This is normally
accomplished by constructing a canonical time coordinate $\tau$ of an
inertial frame (i.e.\ the \emph{barycentred time}), and
decimating/resampling the data at fixed intervals in $\tau$.  The
reconstructed $\tau$ depends on the direction to the source relative
to the Earth's motion; in addition, slow intrinsic parameterized
modulations in the source frequency can also be corrected by this
coordinate transformation.

Most of the routines in this module assume that $\tau$ can be
piecewise expanded as a Taylor series in $t$.  That is, one defines a
set of fitting \emph{regions} $T_i=[t_{\mathrm{bound}(i-1)},
t_{\mathrm{bound}(i)}]$, and a set of fitting \emph{points}
$t_{(i)}\in T_i$.  In each region one then writes:
\begin{equation}
\label{eq:tau}
\tau(t) = \sum_{k=0} \frac{1}{k!}c_{k(i)}(t-t_{(i)})^k \; .
\end{equation}
Since one is normally interested in tracking the difference
$\tau(t)-t$, one can also write the expansion as:
\begin{equation}
\label{eq:delta-tau}
\tau(t)-t = \sum_{k=0} a_{k(i)}(t-t_{(i)})^k \; ,
\end{equation}
where
\begin{eqnarray}
a_{0(i)} & = & c_{0(i)}-t_{(i)}           \; , \nonumber\\
a_{1(i)} & = & c_{1(i)}-1                 \; , \nonumber\\
a_{k(i)} & = & c_{k(i)}/k! \; , \; k\geq2 \; . \nonumber
\label{eq:a_c}
\end{eqnarray}
These are the polynomial coefficients normally assumed in the modules
under this header.

The procedure for resampling according to $\tau$ is normally combined
with \emph{decimating} the time series.  That is, one takes a time
series sampled at constant intervals $\Delta t$ in $t$, and samples it
at constant intervals $d\Delta t$ in $\tau$, where the
\emph{decimation factor} $d$ is normally taken to be an integer
$\geq1$.  When $\tau$ and $t$ are drifting out of phase relatively
slowly, this means that most of the time every $d^\mathrm{th}$ sample
in the original time series becomes the next sample in the decimated
time series.  However, when $\tau$ and $t$ drift out of synch by an
amount $\pm\Delta t$, one can force the decimated time series to track
$\tau$ (rather than $t$) by sampling the $d\pm1^\mathrm{th}$ next
datum (rather than the $d^\mathrm{th}$).  If the drift is sufficiently
rapid or $d$ is sufficiently large, one may be forced to choose the
point $d\pm2$, $d\pm3$, etc.; the size of this adjustment is called
the correction \emph{shift}.  The number of (resampled) time intervals
between one correction point and the next is called the correction
\emph{interval}.

Unless otherwise specified, all time variables and parameters in the
functions under this header can be assumed to measure the detector
time coordinate $t$.  Canonical times are specified by giving the
difference $\tau-t$.

\paragraph{Caveat emptor:} The inclusion of this header and its
associated modules into LAL is provisional at this time.  The routines
and the test code appear to work, but a later standalone code,
operating on much larger datasets, appeared to encounter a memory
leak.  I have not yet determined whether this leak was in the
standalone code or in these LAL routines.

******************************************************* </lalLaTeX> */

#ifndef _RESAMPLE_H
#define _RESAMPLE_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID(RESAMPLEH,"$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define RESAMPLEH_ENUL    1
#define RESAMPLEH_EOUT    2
#define RESAMPLEH_EMEM    3
#define RESAMPLEH_EDTPOS  4
#define RESAMPLEH_ELENGTH 5
#define RESAMPLEH_ETIME   6

#define RESAMPLEH_MSGENUL    "Unexpected null pointer in arguments"
#define RESAMPLEH_MSGEOUT    "Output handle points to a non-null pointer"
#define RESAMPLEH_MSGEMEM    "Memory allocation error"
#define RESAMPLEH_MSGELENGTH "Vector lengths in polyco structure don't argree"
#define RESAMPLEH_MSGEDTPOS  "Sampling interval is not positive"
#define RESAMPLEH_MSGETIME   "Requested output time span extends beyond range of validity of input"
/*************************************************** </lalErrTable> */

/********************************************************** <lalLaTeX>
\subsection*{Types}

\subsubsection*{Structure \texttt{ResampleRules}}
\idx[Type]{ResampleRules}

\noindent This structure stores the rules for taking a time series
$t$, sampled at constant intervals $\Delta t$, and resampling it at
constant intervals $d\Delta t$ in the canonical time coordinate $\tau$,
as described above.  The fields in this structure are as follows:

\begin{description}
\item[\texttt{LIGOTimeGPS start}] The initial time for which the rules
	apply.

\item[\texttt{LIGOTimeGPS stop}] The final time for which the rules
	apply.

\item[\texttt{INT4 length}] The number of correction points, i.e.\
	points where the resampling interval is adjusted from $d\Delta
	t$ to $(d\pm n)\Delta t$.

\item[\texttt{INT4 *interval}] An array giving the number of resampled
	time intervals between correction points.

\item[\texttt{INT2 *shift}] An array giving the size of the correction
	shift (i.e.\ the number $n$ above) at each correction point.

\item[\texttt{INT4 decimate}] The decimation factor $d$.

\item[\texttt{REAL8 deltaT}] The sampling interval before decimation,
	in seconds.

\item[\texttt{REAL8 startDiff}] The difference $\tau-t$ at the time
	\verb@start@, in seconds.

\item[\texttt{REAL8 stopDiff}] The difference $\tau-t$ at the time
	\verb@stop@, in seconds.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagResampleRules{
  LIGOTimeGPS start; /* Detector time at start of resample rules. */
  LIGOTimeGPS stop;  /* Last time for which resample rules apply. */
  INT4 length;       /* Size of the following two arrays. */
  INT4 *interval;    /* Number of samples to the next shift point. */
  INT2 *shift;       /* Size of shift (usually +/- 1). */
  INT4 decimate;     /* Decimation factor. */
  REAL8 deltaT;      /* Sampling rate of the unresampled data. */
  REAL8 startDiff;   /* Offset between tau and t at the start. */
  REAL8 stopDiff;    /* Offset between tau and t at the end. */
} ResampleRules;


/********************************************************** <lalLaTeX>
\subsubsection*{Structure \texttt{PolycoStruc}}
\idx[Type]{PolycoStruc}

\noindent This structure stores the parameters of the piecewise
polynomial fit of $\tau-t$ as a function of $t$.  See
Eq.~\ref{eq:delta-tau} for notation.  The fields of this structure
are:

\begin{description}
\item[\texttt{REAL4 ra}] The right ascension angle of the source, in
	\emph{radians} in the range $[0,2\pi)$.

\item[\texttt{REAL4 dec}] The declination angle of the source, in
	\emph{radians} in the range $[-\pi/2,pi/2]$.

\item[\texttt{REAL4Vector *spindown}] A vector
	$\vec\lambda=(\lambda_0,\ldots,\lambda_{n-1})$ of parameters
	describing a slow intrisic frequency drift $f=f(t)$ of the
	source: $\lambda_k=f^{-1}d^{k+1}f/dt^{k+1}$ at the time given
	by \verb@start@ (below).

\item[\texttt{LIGOTimeGPS start}] The initial time over which the
	polynomial fit applies.

\item[\texttt{REAL4Sequence *tBound}] The sequence of times
	$t_{\mathrm{bound}(i)}$ defining the endpoints of the fitting
	regions, given in seconds after the time \verb@start@.  The
	first fitting region $i=0$ runs from \verb@start@ to
	\verb@start@+$t_{\mathrm{bound}(0)}$, the next from there to
	\verb@start@+$t_{\mathrm{bound}(1)}$, and so on.

\item[\texttt{REAL4Sequence *t0}] The sequence of times $t_{(i)}$ in
	each fitting region at which the polynomial fits are computed,
	given in seconds after the time \verb@start@.

\item[\texttt{REAL4VectorSequence *polyco}] A sequence of vectors
	$\vec a_{(i)}=(a_{0(i)},a_{1(i)},\ldots)$ giving the
	coefficients of the polynomial fit at each time $t_{(i)}$.
	Each element $a_{k(i)}$ has units of $\mathrm{s}^{1-k}$.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagPolycoStruc{
  REAL4 ra;  /* Right ascension of source */
  REAL4 dec; /* Declination of source */
  REAL4Vector *spindown; /* Spindown terms: f0^{-1} d^n f/(dt)^n */
  LIGOTimeGPS start;  /* Start (reference) time of the polyco fit */
  REAL4Sequence *tBound; /* End times of each fitting region */
  REAL4Sequence *t0;     /* Fitting times in each region */
  REAL4VectorSequence *polyco; /* Polynomial fitting parameters for
                                  each fitting region */
} PolycoStruc;


/********************************************************** <lalLaTeX>
\subsubsection*{Structure \texttt{ResampleParamStruc}}
\idx[Type]{ResampleParamStruc}

\noindent This structure stores extra parameters required to construct
a \verb@ResampleRules@ object from a \verb@PolycoStruc@ object.  The
fields of this structure are:

\begin{description}
\item[\texttt{LIGOTimeGPS start}] The initial time for which the
	resample rules will apply.

\item[\texttt{LIGOTimeGPS stop}] The final time for which the resample
	rules will apply.

\item[\texttt{REAL8 deltaT}] The sampling interval before decimation,
	in seconds.

\item[\texttt{INT4 decimate}] The decimation factor.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagResampleParamStruc{
  LIGOTimeGPS start;    /* Initial time for which the rules apply. */
  LIGOTimeGPS stop;     /* Final time for which the rules apply. */
  REAL8       deltaT;   /* Base (oversampled) sampling interval. */
  INT4        decimate; /* Decimation factor. */
} ResampleParamStruc;


/* <lalLaTeX>
\vfill{\footnotesize\input{ResampleHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{CreateResampleRulesC}
</lalLaTeX> */
void
LALCreateResampleRules( LALStatus          *status,
			ResampleRules      **rules,
			PolycoStruc        *polyco,
			ResampleParamStruc *params );

/* <lalLaTeX>
\newpage\input{DestroyResampleRulesC}
</lalLaTeX> */
void
LALDestroyResampleRules( LALStatus     *status,
			 ResampleRules **rules );

/* <lalLaTeX>
\newpage\input{ApplyResampleRulesC}
</lalLaTeX> */
void
LALApplyResampleRules( LALStatus       *status,
		       REAL4TimeSeries *output,
		       REAL4TimeSeries *input,
		       ResampleRules   *rules );

/* Patrick: I don't know if you want ApplyResampleRules() to be a
   lower-level routine that operates on vectors, or whether time
   series are okay. */

/* <lalLaTeX>
\newpage\input{PolycoToTimingDifferenceC}
</lalLaTeX> */
void
LALPolycoToTimingDifference( LALStatus       *status,
			     REAL4TimeSeries *difference,
			     PolycoStruc     *polyco );

/* <lalLaTeX>
\newpage\input{RulesToTimingDifferenceC}
</lalLaTeX> */
void
LALRulesToTimingDifference( LALStatus       *status,
			    REAL4TimeSeries *difference,
			    ResampleRules   *rules );

/* <lalLaTeX>
\newpage\input{ResampleTestC}
</lalLaTeX> */

#ifdef __cplusplus
}
#endif

#endif /* _RESAMPLE_H */
