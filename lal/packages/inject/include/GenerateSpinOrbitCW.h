/************************** <lalVerbatim file="GenerateSpinOrbitCWHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{GenerateSpinOrbitCW.h}}
\label{s:GenerateSpinOrbitCW.h}

Provides routines to generate continuous waveforms with spindown and
orbital modulation.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/GenerateSpinOrbitCW.h>
\end{verbatim}

This header covers routines to generate continuous quasiperiodic
waveforms with a smoothly-varying intrinsic frequency modulated by
orbital motions around a binary companion.  The intrinsic frequency is
modeled by Taylor series coefficients as in \verb@GenerateTaylorCW.h@,
and the orbital modulation is described by a reduced set of orbital
parameters.  Note that the routines do \emph{not} account for spin
precession, accretion processes, or other complicating factors; they
simply Doppler-modulate a polynomial frequency function.

The frequency and phase of the wave in the source's rest frame are
given by Eqs.~(\ref{eq:taylorcw-freq}) and~(\ref{eq:taylorcw-phi}) of
\verb@GenerateTaylorCW.h@, where $t$ is the proper time in this rest
frame.  The frequency and phase of the wave fronts crossing a
reference point in an inertial frame (e.g.\ the Solar system
barycentre) are simply $f[t(t_r)]$ and $\phi[t(t_r)]$, where
\begin{equation}
\label{eq:spinorbit-tr}
t_r = t + R(t)/c
\end{equation}
is the (retarded) time measured at the inertial reference point a
distance $r$ from the source.

The generation of the waveform thus consists of computing the radial
component $R(t)$ of the orbital motion of the source in the binary
centre-of-mass frame, inverting Eq.~(\ref{eq:spinorbit-tr}) to find
the ``emission time'' $t$ for a given ``detector time'' $t_r$, and
plugging this into the Taylor expansions to generate the instantaneous
frequency and phase.  The received frequency is also multiplied by the
instantaneous Doppler shift $[1+\dot{R}(t)/c]^{-1}$ at the time of
emission.

Since we do not include precession effects, the polarization state of
the wave is constant: we simply specify the polarization amplitudes
$A_+$, $A_\times$ and the polarization phase $\psi$ based on the
(constant) orientation of the source's \emph{rotation}.  The following
discussion defines a set of parameters for the source's orbital
\emph{revolution}, which we regard as completely independent from its
rotation.

\subsubsection*{Orbital motion}

\begin{wrapfigure}{r}{0.52\textwidth}
\vspace{-4ex}
\begin{center}
\resizebox{0.47\textwidth}{!}{\includegraphics{inject_binary}} \\
\parbox{0.47\textwidth}{\caption{\label{fig:binary-orbit} Binary orbit
orientation parameters.}}
\end{center}
\vspace{-2ex}
\end{wrapfigure}
Fig.~\ref{fig:binary-orbit} illustrates the notation conventions
defining a binary orbit.  We define a radial axis $R$ directed
\emph{from} the observer (Earth) \emph{to} the source, as shown.  The
horizontal plane is thus the plane of the sky, and the direction
marked $N$ is the direction along a meridian towards the North
celestial pole.  The tilted plane is the plane of the binary orbit,
and the axis labeled $z$ is the normal to this plane directed such
that the orbit is right-handed about this axis.  The \emph{ascending
node} of the orbit, denoted by
\raisebox{-0.5pt}{\includegraphics{inject_ascend}}, is the direction
defined by $\hat{\mathbf{\mathit{R}}}\times\hat{\mathbf{\mathit{z}}}$.
The binary orbit itself is shown as an off-centred ellipse, with the
barycentre at one of its foci; the wave-emitting source is also shown.

The \emph{inclination angle} $i$ is the angle between the sky and
orbital planes.  The \emph{longitude of the ascending node} $\Omega$
is the angle in the plane of the sky from the North direction to the
ascending node, measured right-handed about
$\hat{\mathbf{\mathit{R}}}$.  The \emph{argument of the periapsis}
$\omega$ is the angle in the orbital plane from the ascending node to
the direction of periapsis (point where the source is closest to the
system barycentre), and the \emph{true anomaly} $\upsilon(t)$ of the
source is the angle from the periapsis to the current location of the
source; both angles are measured right-handed about
$\hat{\mathbf{\mathit{z}}}$ (i.e.\ prograde).  The \emph{periapsis
separation} $r_p$ is the distance from the periapsis to the
barycentre, and we denote the \emph{eccentricity} of the orbital
ellipse as $e$, so that the separation between the source and the
barycentre at any time is $r=r_p(1+e)/(1+e\cos\upsilon)$.

In this convention, $i\in[0,\pi]$ and $\Omega\in[0,2\pi)$.  Another
convention common in astronomy is to restrict $\Omega$ to the range
$[0,\pi)$, refering to whichever node (ascending or descending) lies
in this range.  The argument of the periapsis $\omega$ is then also
measured from this node.  In this case the range of $i$ must be
extended to $(-\pi,\pi]$; it is negative if the reference node is
descending, and positive if it is ascending.  The formulae that follow
are the same in either convention, though, since one can verify that
adding $\pi$ to $\Omega$ and $\omega$ is equivalent to reversing the
sign on $i$.

Some spherical trigonometry gives us $R=r\sin(\omega+\upsilon)\sin i$.
We can differentiate $R$ with respect to $t$, and apply Keplers
second law
$r^2\dot{\upsilon}=r_p^2\dot{\upsilon}_p=\mathrm{constant}$, where
$\dot{\upsilon}_p$ is the angular speed at periapsis, to get:
\begin{eqnarray}
\label{eq:orbit-r}
R & = & R_0 + \frac{(1+e) r_p\sin i}{1+e\cos\upsilon}
	\sin(\omega+\upsilon) \;,\\
\label{eq:orbit-rdot}
\dot{R} & = & \dot{R}_0 + \frac{\dot{\upsilon}_p r_p\sin i}{1+e}
	\left[ \cos(\omega+\upsilon) + e\cos\omega \right] \;.
\end{eqnarray}
Without loss of generality, we will henceforth drop the offsets $R_0$
and (constant) $\dot{R}_0$ from these equations.  This means that we
ignore the overall propagation delay between the $R=R_0$ plane and the
observer, and incorporate any (constant) Doppler shifts due to net
centre-of-mass motions into the values of $f$ and $\dot{\upsilon}_p$.
The resulting times and parameter values are referred to as being in
the \emph{barycentric} frame.  The only time delays and Doppler shifts
that we explicitly treat are those arising from the motion of the
source relative to the $R=R_0$ sky plane passing through the system
barycentre.

All we need now to determine the orbital motion is an equation for
$\upsilon(t)$.  Many basic astronomy textbooks give exact but
transcendental expressions relating $\upsilon$ and $t$ for elliptical
orbits with $0\leq e<1$, and/or series expansions of $\upsilon(t)$ for
$e\ll1$.  However, for a generic binary system we cannot guarantee
that $e\ll1$, and for now we would like to retain the possibility of
modeling open orbits with $e\geq1$.  For now we will simply present
the exact formulae, and discuss the numerical solution methods in the
modules under this header.

Let $t_p$ be the time of a periapsis passage (preferably a recent one
in the case of closed orbits).  We express both $t$ and $\upsilon$ in
terms of an intermediate variable $E$ (called the \emph{eccentric
anomaly} for elliptic orbits, unnamed for open orbits).  The formulae
are:
\begin{equation}
\label{eq:spinorbit-t}
t - t_p = \left\{ \begin{array}{l@{\qquad}c}
	\frac{1}{\dot{\upsilon}_p} \sqrt{\frac{1+e}{(1-e)^3}}
		\left( E - e\sin E \right) & 0 \leq e < 1 \\ & \\
	 \frac{1}{\dot{\upsilon}_p} E
		\left( 1 + \frac{E^2}{12} \right) & e = 1 \\ & \\
	 \frac{1}{\dot{\upsilon}_p} \sqrt{\frac{e+1}{(e-1)^3}}
		\left( e\sinh E - E \right) & e > 1
\end{array} \right.
\end{equation}
\begin{equation}
\label{eq:spinorbit-upsilon}
\begin{array}{c} \tan\left(\frac{\upsilon}{2}\right) \end{array}
= \left\{ \begin{array}{l@{\qquad}c}
	\sqrt{\frac{1+e}{1-e}}\tan\left(\frac{E}{2}\right)
		& 0 \leq e < 1 \\ & \\
	\frac{E}{2} & e = 1 \\ & \\
	\sqrt{\frac{e+1}{e-1}}\tanh\left(\frac{E}{2}\right) & e > 1
\end{array} \right.
\end{equation}

Thus to solve for $\upsilon(t)$ one typically inverts the equation for
$t-t_p$ numerically or by series expansion, finds the corresponding
$E$, and then plugs this into the expression for $\upsilon$.  However,
in our case we would then need to do another numerical inversion to
find the retarded time $t_r$ from Eq.~(\ref{eq:spinorbit-tr}).  A more
efficient approach is thus to take an initial guess for $E$, compute
both $t$, $\upsilon$, and hence $t_r$, and then refine directly on
$E$.

\subsubsection*{Other notation conventions}

Since we may deal with highly eccentric or open orbits, we will
specify these orbits with parameters that are definable for all
classes of orbit.  Thus we specify the size of the orbit with the
periapsis separation $r_p$ rather than the semimajor axis $a$, and the
speed of the orbit with the angular speed at periapsis
$\dot{\upsilon}_p$ rather than with the period $P$.  These parameters
are related by:
\begin{eqnarray}
\label{eq:spinorbit-a}
a & = & \frac{r_p}{1-e} \;,\\
\label{eq:spinorbit-p}
P & = & \frac{2\pi}{\dot{\upsilon}_p} \sqrt{\frac{1+e}{(1-e)^3}} \;.
\end{eqnarray}
Furthermore, for improved numerical precision when dealing with
near-parabolic orbits, we specify the value of $1-e$ rather than the
value of $e$.  We note that $1-e$ has a maximum value of $1$ for a
circular orbit, positive for closed elliptical orbits, zero for
parabolic orbits, and negative (unbounded) for hyperbolic orbits.

******************************************************* </lalLaTeX> */

#ifndef _GENERATESPINORBITCW_H
#define _GENERATESPINORBITCW_H

#include <lal/LALStdlib.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>

#ifdef  __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( GENERATESPINORBITCWH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define GENERATESPINORBITCWH_ENUL 1
#define GENERATESPINORBITCWH_EOUT 2
#define GENERATESPINORBITCWH_EMEM 3
#define GENERATESPINORBITCWH_EECC 4
#define GENERATESPINORBITCWH_EFTL 5
#define GENERATESPINORBITCWH_ESGN 6

#define GENERATESPINORBITCWH_MSGENUL "Unexpected null pointer in arguments"
#define GENERATESPINORBITCWH_MSGEOUT "Output field a, f, phi, or shift already exists"
#define GENERATESPINORBITCWH_MSGEMEM "Out of memory"
#define GENERATESPINORBITCWH_MSGEECC "Eccentricity out of range"
#define GENERATESPINORBITCWH_MSGEFTL "Periapsis motion is faster than light"
#define GENERATESPINORBITCWH_MSGESGN "Sign error: positive parameter expected"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{SpinOrbitCWParamStruc}}
\idx[Type]{SpinOrbitCWParamStruc}

This structure stores the parameters for constructing a gravitational
waveform with both a Taylor-polynomial intrinsic frequency and phase,
and a binary-orbit modulation.  As with the \verb@PPNParamStruc@ type
in \verb@GeneratePPNInspiral.h@, we divide the fields into passed
fields (which are supplied to the final \verb@CoherentGW@ structure
but not used in any calculations), input fields (that are used by the
waveform generator), and output fields (that are set by the waveform
generator).  They are:

\bigskip\noindent\textit{Passed fields:}
\begin{description}
\item[\texttt{SkyPosition position}] The location of the source on the
sky, normally in equatorial coordinates.

\item[\texttt{REAL4 psi}] The polarization angle of the source, in
radians.
\end{description}

\medskip\noindent\textit{Input fields:}
\begin{description}
\item[\texttt{LIGOTimeGPS epoch}] The start time of the output series.

\item[\texttt{LIGOTimeGPS spinEpoch}] A reference time
$t_\mathrm{ref}$ (in the barycentric frame) at which the rotational
properties of the source are specified.

\item[\texttt{LIGOTimeGPS orbitEpoch}] A time $t_\mathrm{peri}$ (in
the barycentric frame) at which the source passes through periapsis.
Note that this is the proper or ``true'' time of passage; the
\emph{observed} periapsis passage occurs at time
$t_\mathrm{peri}+r(t_\mathrm{peri})/c$.

\item[\texttt{REAL8 deltaT}] The requested sampling interval of the
waveform, in s.

\item[\texttt{UINT4 length}] The number of samples in the generated
waveform.

\item[\texttt{REAL4 aPlus, aCross}] The polarization amplitudes $A_+$,
$A_\times$, in dimensionless strain units.

\item[\texttt{REAL8 phi0}] The phase of the wave emitted at time
$t_\mathrm{ref}$, in radians.

\item[\texttt{REAL8 f0}] The frequency of the wave emitted at time
$t_\mathrm{ref}$ (and incorporating any Doppler shift due to
$\dot{R}_0$), in Hz.

\item[\texttt{REAL8Vector *f}] The spin-normalized Taylor parameters
$f_k$, as defined in Eq.~\ref{eq:taylorcw-freq} of
\verb@GenerateTaylorCW.h@.  If \verb@f@=\verb@NULL@, the (proper) spin
of the source is assumed to be constant.

\item[\texttt{REAL8 omega}] The argument of the periapsis, $\omega$,
in radians.

\item[\texttt{REAL8 rPeriNorm}] The projected,
speed-of-light-normalized periapsis separation of the orbit,
$(r_p/c)\sin i$, in s.

\item[\texttt{REAL8 oneMinusEcc}] The value of $1-e$.

\item[\texttt{REAL8 angularSpeed}] The angular speed at periapsis,
$\dot{\upsilon}_p$, in Hz.
\end{description}

\medskip\noindent\textit{Output fields:}
\begin{description}
\item[\texttt{REAL4 dfdt}] The maximum value of $\Delta f\Delta t$
encountered over any timestep $\Delta t$ used in generating the
waveform.
\end{description}

******************************************************* </lalLaTeX> */

/**
 * This structure stores the parameters for constructing a gravitational
 * waveform with both a Taylor-polynomial intrinsic frequency and phase,
 * and a binary-orbit modulation.  As with the PPNParamStruc type
 * in GeneratePPNInspiral.h, we divide the fields into passed
 * fields (which are supplied to the final CoherentGW structure
 * but not used in any calculations), input fields (that are used by the
 * waveform generator), and output fields (that are set by the waveform
 * generator).
 */
typedef struct tagSpinOrbitCWParamStruc {
  /* Passed parameters. */
  SkyPosition position;   /**< location of source on sky */
  REAL4 psi;              /**< polarization angle (radians) */

  /* Input parameters. */
  LIGOTimeGPS epoch;      /**< start time of output time series */
  LIGOTimeGPS spinEpoch;  /**< reference time for rotational parameters */
  LIGOTimeGPS orbitEpoch; /**< time of a periapsis passage */
  REAL8 deltaT;           /**< requested sampling interval (s) */
  UINT4 length;           /**< length of time series */
  REAL4 aPlus, aCross;    /**< polarization amplitudes */
  REAL8 phi0;             /**< initial phase (radians) */
  REAL8 f0;               /**< initial frequency (Hz) */
  REAL8Vector *f;         /**< f0-normalized Taylor parameters */
  REAL8 omega;            /**< argument of periapsis (radians) */
  REAL8 rPeriNorm;        /**< projected, normalized periapsis (s) */
  REAL8 oneMinusEcc;      /**< 1 - orbital eccentricity */
  REAL8 angularSpeed;     /**< angular speed at periapsis (Hz) */

  /* Output parameters. */
  REAL4 dfdt;             /**< [OUT:] maximum value of df*dt over any timestep */
} SpinOrbitCWParamStruc;


/* <lalLaTeX>
\vfill{\footnotesize\input{GenerateSpinOrbitCWHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{GenerateSpinOrbitCWC}
</lalLaTeX> */
void
LALGenerateSpinOrbitCW( LALStatus             *,
			CoherentGW            *output,
			SpinOrbitCWParamStruc *params );

/* <lalLaTeX>
\newpage\input{GenerateEllipticSpinOrbitCWC}
</lalLaTeX> */
void
LALGenerateEllipticSpinOrbitCW( LALStatus             *,
				CoherentGW            *output,
				SpinOrbitCWParamStruc *params );

/* <lalLaTeX>
\newpage\input{GenerateParabolicSpinOrbitCWC}
</lalLaTeX> */
void
LALGenerateParabolicSpinOrbitCW( LALStatus             *,
				 CoherentGW            *output,
				 SpinOrbitCWParamStruc *params );

/* <lalLaTeX>
\newpage\input{GenerateHyperbolicSpinOrbitCWC}
</lalLaTeX> */
void
LALGenerateHyperbolicSpinOrbitCW( LALStatus             *,
				  CoherentGW            *output,
				  SpinOrbitCWParamStruc *params );

/* <lalLaTeX>
%\newpage\input{SimulateSpinOrbitCWTestC}
</lalLaTeX> */

#ifdef  __cplusplus
#pragma {
}
#endif

#endif /* _GENERATESPINORBITCW_H */
