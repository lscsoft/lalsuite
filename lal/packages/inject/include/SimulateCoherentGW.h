/*************************** <lalVerbatim file="SimulateCoherentGWHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\bm}[1]{\mbox{\boldmath$#1$\unboldmath}}
\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\section{Header \texttt{SimulateCoherentGW.h}}
\label{s:SimulateCoherentGW.h}

Provides routines to simulate generic gravitational waveforms
originating from a particular source.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/SimulateCoherentGW.h>
\end{verbatim}

This header covers generic routines and structures to represent and
simulate the effects of a plane gravitational wave propagating from a
distinct point on the sky.

Any plane gravitational wave is specified by a direction
$\bm{\hat{n}}$ to its apparent source (i.e.\ opposite to its direction
of propagation), and by the inistantaneous values $h_+(t)$,
$h_\times(t)$ of its plus and cross polarizations as functions of
(retarded) time $t=t_0+\bm{\hat{n}}\cdot(\bm{x}-\bm{x}_0)$, where
$t_0$ is the time measured at some local reference point $\bm{x}_0$,
and $t$ is the time measured by a synchronized clock at $\bm{x}$.  We
adopt the standard meaning of the instantaneous strain amplitudes
$h_{+,\times}$: in some reference transverse $x$-$y$ coordinate system
oriented such that $\bm{\hat{x}}\times\bm{\hat{y}}=-\bm{\hat{n}}$
points in the direction of propagation, two free observers originally
separated by a displacement $(x,y)$ will experience an additional
tidal displacement $\delta x=(xh_+ + yh_\times)/2$, $\delta
y=(xh_\times - yh_+)/2$.

\paragraph{Quasiperiodic waves:} Most astrophysical sources of
gravitational radiation are described as \emph{quasiperiodic} (or,
less accurately, as ``adiabatic''), in that they can be said to have
an instantaneous frequency, amplitude, and polarization, all of which
vary on timescales much longer than a wave period. Mathematically we
write this as:
$$
h_{+,\times}(t) = A_{+,\times}(t)\cos\phi(t)
	+ B_{+,\times}(t)\sin\phi(t) \; ,
$$
\begin{wrapfigure}{r}{0.47\textwidth}
\vspace{-4ex}
\begin{center}
\resizebox{0.42\textwidth}{!}{\includegraphics{inject_phase_diagram}} \\
\parbox{0.42\textwidth}{\caption{\label{fig:phase-diagram}
Polarization phase diagram for a quasiperiodic gravitational wave.
The phase point $p(t)$ traces out the indicated ellipse in the
$h_+,h_\times$ plane; the parameters $A_1$, $A_2$, and $\Phi$
remain roughly constant over many cycles in $\phi$.}}
\end{center}
\vspace{-2ex}
\end{wrapfigure}
where $\phi(t)=2\pi\int f(t)\,dt$, and the \emph{evolution timescale}
$\tau=\min\{A/\dot{A},B/\dot{B},f/\dot{f}\}$ is much greater than
$h/\dot{h}\sim1/f$.  Obviously it is mathematically impossible for the
physical functions $h_{+,\times}(t)$ to specify uniquely more than two
other functions of time; we rely on the notion of quasiperiodicity to
define ``natural'' choices of instantaneous frequency and amplitude.
The ambiguity in this choice is on the order of the amount that these
quantities change over a cycle.

While the above formula appears to have five degrees of freedom (two
quadrature amplitudes $A$ and $B$ for each polarization, plus a common
phase function $\phi$), there is a degeneracy between the two
quadrature amplitudes and a shift in phase.  One could simply treat
each polarization independently and represent the system with two
amplitude functions $A_{+,\times}$ and two phase functions
$\phi_{+,\times}$, but we would like to preserve the notion that the
phases of the two waveforms derive from a single underlying
instantaneous frequency.  We therefore write the waveforms in terms of
two polarization amplitudes $A_1(t)$ and $A_2(t)$, a single phase
function $\phi(t)$, and a polarization shift $\Phi(t)$:
\begin{eqnarray}
\label{eq:quasiperiodic-hplus}
h_+(t) & = & A_1(t)\cos\Phi(t)\cos\phi(t)
		- A_2(t)\sin\Phi(t)\sin\phi(t) \; , \\
\label{eq:quasiperiodic-hcross}
h_\times(t) & = & A_1(t)\sin\Phi(t)\cos\phi(t)
		+ A_2(t)\cos\Phi(t)\sin\phi(t) \; .
\end{eqnarray}
The physical meaning of these functions is shown in
Fig.~\ref{fig:phase-diagram}.

There is a close relationship between the polarization shift $\Phi$
and the orientation of the $x$-$y$ coordinates used to define our
polarization basis: if we rotate the $x$ and $y$ axes by an angle
$\Delta\psi$, we change $\Phi$ by an amount $-2\Delta\psi$.  (The
factor of 2 comes from the fact that the + and $\times$ modes are
quadrupolar: a + mode rotated $45^\circ$ is a $\times$ mode, and a
mode rotated $90^\circ$ is the opposite of itself.)  We use the
\emph{polarization angle} $\psi$ to define the orientation of the
$x$-axis of the polarization basis relative to an Earth-fixed
reference frame (see the coordinate conventions below).  If $\Phi$ is
constant, one can redefine $\psi$ such that $\Phi=0$; however, when
$\Phi$ changes with time, we would nonetheless like our polarization
basis to remain fixed.  We therefore retain the constant $\psi$ and
the function $\Phi(t)$ as distinct quantities.

The advantage of this quasiperiodic representation of a gravitational
wave is that a physical sampling of the parameters $A_1$, $A_2$,
$\phi$, and $\Phi$ need only be done on timescales $\Delta
t\lessim\tau$, whereas the actual wave functions $h_{+,\times}$ need
to be sampled on timescales $\Delta t\lessim1/f$.

The following coordinate conventions are assumed:
\begin{enumerate}
\item Fig.~7 of~\cite{Will_C:1996} defines standard coordinate
conventions for nonprecessing binaries, and by extension, for any
fixed-axis rotating source: If $\bm{\hat{z}}$ points in the direction
of wave propagation (away from the source), and $\bm{\hat{l}}$ points
in the (constant) direction of the source's angular momentum vector,
then the $x$-$y$ coordinates used to define the + and $\times$
polarizations are given by $\bm{\hat{x}}=|\csc
i|\bm{\hat{z}}\times\bm{\hat{l}}$ and
$\bm{\hat{y}}=\bm{\hat{z}}\times\bm{\hat{x}}$, where
$i=\arccos(\bm{\hat{z}}\cdot\bm{\hat{l}})$ is the inclination angle
between $\bm{\hat{l}}$ and $\bm{\hat{z}}$.  Such a system will
generically have $A_1(t)=A(t)(1+\cos^2i)$, $A_2(t)=2A(t)\cos i$,
$\Phi(t)=0$, and $f(t)>0$ (i.e.\ $\phi(t)$ increasing with time).  For
precessing systems, prescriptions for $\bm{\hat{x}}$ and
$\bm{\hat{y}}$ become ambiguous, but they \emph{must} be fixed; the
relations for $A_1$, $A_2$, and $\Phi$ will no longer be maintained.

\item Appendix~B of~\cite{Anderson_W:2000} defines a convention for
the overal polarization angle $\psi$: Let $\bm{\hat{N}}$ be the
direction of the Earth's north celestial pole, and define the
direction of the \emph{ascending node}
$\bm{\hat{\Omega}}=|\csc\alpha|\bm{\hat{N}}\times\bm{\hat{z}}$, where
$\alpha$ is the right ascension of the source.  Then $\psi$ is the
angle, right-handed about $\bm{\hat{z}}$, from $\bm{\hat{\Omega}}$ to
$\bm{\hat{x}}$.

\item The direction of propagation of the wave is defined by the right
ascension $\alpha$ and declination $\delta$ of the \emph{source}, as
seen from the point of measurement.  See \verb@SkyCoordinates.h@ for a
definition of these quantities.  We expect that these will be
effectively constant for almost any gravitational wave source of
interest.
\end{enumerate}

\paragraph{The polarization response:} The relative strain induced in
the test masses of a detector by a passing gravitational wave depends
not only on the amplitudes $h_{+,\times}$ of the gravitational wave,
but also on the design of the detector and its orientation with
relative to the $x$-$y$ coordinate system used to define the + and
$\times$ polarizations.  For a given detector, the response to each
polarization thus depends on the right ascension $\alpha$, declination
$\delta$, and polarization angle $\psi$ of the source (which define
the orientation of the + and $\times$ polarization axes relative to
the Earth), and on the time $t$ (which determines the orientation of
the detector as the Earth rotates).  The strain $h(t)$ induced in the
detector is thus given by two polarization response functions
$F_{+,\times}(\alpha,\delta,\psi;t)$ by:
$$
h(t) = h_+(t)F_+(\alpha,\delta,\psi;t) +
	h_\times(t)F_\times(\alpha,\delta,\psi;t) \; .
$$
We will not discuss the computation of these functions $F_{+,\times}$,
as these are covered under the header \verb@DetResponse.h@.

\paragraph{The transfer function:} All gravitational wave detectors
incorporate a set of analog and digital filters that convert a
gravitational excitation on the test masses into a measurable output
time series.  The effects of these functions are aggregated into a
complex-valued \emph{transfer function} ${\cal T}(f)$, which gives the
instrumental response (in units of ``counts'' from an
analog$\rightarrow$digital converter) to gravitational waves of unit
amplitued at the frequency $f$.  Specifically, if the strain exerted
on the antenna is given by $h(t)=\mathrm{Re}[{\cal H}e^{2\pi ift}]$
(where the complex amplitude $\cal H$ includes the phase of the wave),
then the ADC output of the instrument is given by:
$$
o(t) = \mathrm{Re}\left[ {\cal T}(f) {\cal H}e^{2\pi ift} \right] \; .
$$
The transfer function has a strong frequency dependence in order to
``whiten'' the highly-coloured instrumental noise, and thus preserve
instrumental sensitivity across a broad band of frequencies.

We note that although the transfer function measures the response of
the instrument to a gravitational wave, the term \emph{response
function} refers to inverse transformation of taking an instrumental
response and computing a gravitational waveform; that is, ${\cal
R}(f)=1/{\cal T}(f)$.  This confusing bit of nomenclature arises from
the fact that most data analysis deals with extracting gravitational
waveforms from the instrumental output, rather than injecting
waveforms into the output.

For quasiperiodic waveforms with a well-defined instantaneous
frequency $f(t)$ and phase $\phi(t)$, we can compute the response of
the instrument entirely in the time domain in the adiabatic limit: if
our instrumental excitation is a linear superposition of waveforms
$h(t)=\mathrm{Re}[{\cal H}(t)e^{i\phi(t)}]$, then the output is a
superposition of waves of the form
$$
o(t) \approx \mathrm{Re}\left[ {\cal T}\{f(t)\}
	{\cal H}(t)e^{i\phi(t)} \right] \; .
$$
This expression is approximate to the extent that ${\cal T}(f)$ varies
over the range $f\pm1/\tau$, where $\tau$ is the evolution timescale
of ${\cal H}(t)$ and $f(t)$.  Since the transfer function and
polarization response (above) are linear operators, we can apply them
in either order.

\paragraph{A note on terminology:} We use the word ``coherent'' in the
name of this header in the loosest possible sense, refering to any
wave with a well-defined direction of propagation, whose wave
amplitudes $h_{+,\times}$ are deterministic functions of retarded
time.  Given a knowledge of these parameters, such a waveform is
amenable to ``coherent'' detection in a network of detectors, through
time-shifted matched filtering.

However, coherence is often used to refer to a more restricted class
of waveforms that are ``effectively monochromatic'' over some
coherence timescale $t_\mathrm{coh}$; i.e.\ in any timespan
$t_\mathrm{coh}$ there is a fixed-frequency sinusoid that is never
more than $90^\circ$ out of phase with the waveform.  This is more
retrictive even than our concept of quasiperiodic waves; for
smoothly-varying waveforms one has $t_\mathrm{coh}\sim\dot{f}^{-1/2}$,
which is much shorter than the evolution timescale $\tau\sim
f/\dot{f}$ (provided $\tau\gg1/f$, as we have assumed).

******************************************************* </lalLaTeX> */

#ifndef _SIMULATECOHERENTGW_H
#define _SIMULATECOHERENTGW_H

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/SkyCoordinates.h>
#include <lal/LALBarycenter.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( SIMULATECOHERENTGWH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define SIMULATECOHERENTGWH_ENUL  1
#define SIMULATECOHERENTGWH_EBAD  2
#define SIMULATECOHERENTGWH_ESIG  3
#define SIMULATECOHERENTGWH_EDIM  4
#define SIMULATECOHERENTGWH_EMEM  5
#define SIMULATECOHERENTGWH_EUNIT 6

#define SIMULATECOHERENTGWH_MSGENUL  "Unexpected null pointer in arguments"
#define SIMULATECOHERENTGWH_MSGEBAD  "A sampling interval is (effectively) zero"
#define SIMULATECOHERENTGWH_MSGESIG  "Input signal must specify amplitude and phase functions"
#define SIMULATECOHERENTGWH_MSGEDIM  "Amplitude must be a 2-dimensional vector"
#define SIMULATECOHERENTGWH_MSGEMEM  "Memory allocation error"
#define SIMULATECOHERENTGWH_MSGEUNIT "Bad input units"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{CoherentGW}}
\idx[Type]{CoherentGW}

\noindent This structure stores a representation of a plane
gravitational wave propagating from a particular point on the sky.
Several alternate representations are permitted to allow a more
natural characterization of quasiperiodic waveforms.  The fields are:

\begin{description}
\item[\texttt{SkyPosition position}] The location of the source in the
sky.  This should be in equatorial celestial coordinates, but routines
may be able to do the conversion.

\item[\texttt{REAL4 psi}] The polarization angle $\psi$, in radians,
as defined in Appendix~B of~\cite{Anderson_W:2000}.

\item[\texttt{REAL4TimeVectorSeries *h}] A time-sampled
two-dimensional vector storing the waveforms $h_+(t)$ and
$h_\times(t)$, in dimensionless strain.

\item[\texttt{REAL4TimeVectorSeries *a}] A time-sampled
two-dimensional vector storing the amplitudes $A_1(t)$ and $A_2(t)$,
in dimensionless strain.

\item[\texttt{REAL4TimeSeries *f}] A time-sampled sequence storing the
instantaneous frequency $f(t)$, in Hz.

\item[\texttt{REAL8TimeSeries *phi}] A time-sampled sequence storing
the phase function $\phi(t)$, in radians.

\item[\texttt{REAL4TimeSeries *shift}] A time-sampled sequence storing
the polarization shift $\Phi(t)$, in radians.
\end{description}

\noindent It is permissible to set only some of the
\verb@REAL4TimeSeries@ or \verb@REAL4TimeVectorSeries@ fields above,
but the waveform is treated as being zero except during those times
when either \verb@h@, or both \verb@a@ and \verb@phi@, are defined.
Where \verb@shift@ is not specified, it is assumed that $\Phi$ is
zero; where \verb@f@ is not specified but \verb@phi@ is, $f(t)$ can be
computed as $\dot{\phi}(t)/2\pi$.  Where \verb@f@ and \verb@phi@
overlap, or where \verb@h@ and any other time series overlap, they
must be defined consistently.

******************************************************* </lalLaTeX> */
/**  This structure stores a representation of a plane
gravitational wave propagating from a particular point on the sky.
Several alternate representations are permitted to allow a more
natural characterization of quasiperiodic waveforms.*/
typedef struct tagCoherentGW {
  SkyPosition position;     /**< sky position of source */
  REAL4 psi;                /**< polarization angle of source */
  REAL4TimeVectorSeries *h; /**< sampled waveforms \f$h_+, h_\times\f$ */
  REAL4TimeVectorSeries *a; /**< amplitudes \f$A_+, A_\times\f$ */
  REAL4TimeSeries *f;       /**< instantaneous frequency */
  REAL8TimeSeries *phi;     /**< phase function */
  REAL4TimeSeries *shift;   /**< polarization shift Phi */
} CoherentGW;

/********************************************************** <lalLaTeX>

\subsubsection*{Structure \texttt{DetectorResponse}}
\idx[Type]{DetectorResponse}

\noindent This structure contains information required to determine
the response of a detector to a gravitational waveform.  The fields
are:

\begin{description}
\item[\texttt{COMPLEX8FrequencySeries *transfer}] The
frequency-dependent transfer function of the interferometer, in ADC
counts per unit strain amplitude at any given frequency.  If absent,
the response will be given in raw strain rather than ADC output.

\item[\texttt{LALDetector *site}] A structure storing site and
polarization information, used to compute the polarization response
and the propagation delay.  If absent, the response will be computed
to the plus mode waveform with no time delay.

\item[\texttt{EphemerisData *ephemerides}] A structure storing the
positions, velocities, and accelerations of the Earth and Sun centres
of mass, used to compute the propagation delay to the solar system
barycentre.  If absent, the propagation delay will be computed to the
Earth centre (rather than a true barycentre).

\item[\texttt{LIGOTimeGPS heterodyneEpoch}] A reference time for
heterodyned detector output time series, where the phase of the mixing
signal is zero.  This parameter is only used when generating detector
output time series with nonzero heterodyne frequency \verb@f0@.
(Note: This should really be a parameter stored in the
\verb@TimeSeries@ structure along with \verb@f0@, but it isnt, so we
have to add it here.)
\end{description}

******************************************************* </lalLaTeX> */
/** This structure contains information required to determine
    the response of a detector to a gravitational waveform. */
typedef struct tagDetectorResponse {
  COMPLEX8FrequencySeries *transfer; 	/**< frequency transfer function */
  LALDetector *site;           		/**< detector location and orientation */
  EphemerisData *ephemerides;  		/**< Earth and Sun ephemerides */
  LIGOTimeGPS heterodyneEpoch; 		/**< reference time for heterodyning */
} DetectorResponse;


/* <lalLaTeX>
\vfill{\footnotesize\input{SimulateCoherentGWHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{SimulateCoherentGWC}
</lalLaTeX> */
void
LALSimulateCoherentGW( LALStatus        *,
		       REAL4TimeSeries  *output,
		       CoherentGW       *input,
		       DetectorResponse *detector );

void
LALSimulateCoherentGW_exp (LALStatus        *stat,
			   REAL4TimeSeries  *output,
			   CoherentGW       *input,
			   DetectorResponse *detector );

/* <lalLaTeX>
%\newpage\input{SimulateCoherentGWTestC}
</lalLaTeX> */

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _SIMULATECOHERENTGW_H */
