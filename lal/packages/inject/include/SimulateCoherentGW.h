/*************************** <lalVerbatim file="SimulateCoherentGWHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

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

Any plane gravitational wave can be represented by a direction of
propagation (e.g. by the right ascension and declination of its source
on the sky), and by the instantaneous values $h_+(t)$, $h_\times(t)$
of its plus and cross polarizations as functions of time.  However,
representing these functions smoothly as sampled time series can
require very high sampling rates, and many signals of astrophysical
interest can be characterized much more coarsely.  A good generic
representation of a gravitational waveform should allow for other more
natural characterizations.

In many cases a waveform is ``coherent'' or ``adiabatic''; i.e.\ one
can define amplitudes $A_{+,\times}(t)$ and frequencies
$f_{+,\times(t)}$ that vary over timescales much longer than
$h_{+,\times}(t)$, such that
$$
h_{+,\times}(t) = A_{+,\times}(t) \cos\phi_{+,\times}(t) \; ,
$$
where $\phi_{+,\times}(t)=2\pi\int f_{+,\times}(t)\,dt$.  Obviously,
the two amplitudes and frequencies are not uniquely defined from the
two waveforms $h_{+,\times}(t)$; one requires the notion of
adiabaticity to define a ``natural'' choice of instantaneous frequency
and amplitude.  However, many signals of interest are coherent in this
sense.  The structure \verb@CoherentGW@, below, allows waves to be
represented in this highly-undersampled way by specifying
$A_{+,\times}(t)$ and $\phi_{+,\times}(t)$ instead of
$h_{+,\times}(t)$.

Furthermore, most coherent astrophysical waves come from rotating
systems, which produce waveforms that are circularly polarized about
the axis of rotation.  According to the standard polarization
convention given in Fig.~7 of~\cite{Will_C:1996}, where the $z$-axis
points in the direction of propagation, the reference $x$-axis for the
+ and $\times$ polarization tensors is the ascending line of nodes
between the plane transverse to $z$ and the invariant plane of
rotation (i.e.\ it is a ray formed by the the intersection of the
transverse and invariant rotation planes, oriented such that the
rotation of the system carries bodies in the $+z$ direction as they
cross that ray).  In this scheme, $h_+(t)$ and $h_\times(t)$ are
uniformly $90^\circ$ out of phase.  The structure \verb@CoherentGW@
therefore allows one to specify a single phase function $\phi_+(t)$,
on the assumption that $\phi_\times(t)=\phi_+(t)-90^\circ$.

In addition to the coordinate convention defining the + and $\times$
polarization tensors, we impose the following additional conventions.
Sky positions refer to the normal geocentric right ascension (measured
east along the equatorial plane from the vernal equinox) and
declination (measured north from the equatorial plane).  The
polarization angle $\psi$ follows the convention specified in
Appendix~B in~\cite{Anderson_W:2000}, which is as follows: Define the
line of nodes as the intersection between the equatorial plane and the
plane transverse to the wave propagation, and the ascending node as
the ray from the origin along this line such that a point, rotating in
the right-hand sense about the propagation direction $z$, moves from
south to north as it crosses that ray.  Then $\psi$ is the angle
between the ascending node and the $x$-axis defined above, measured in
the right-hand sense about the $z$-axis.  Routines using the
\verb@CoherentGW@ structure (below) with different conventions should
document those conventions.

******************************************************* </lalLaTeX> */

#ifndef _SIMULATECOHERENTGW_H
#define _SIMULATECOHERENTGW_H

#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( SIMULATECOHERENTGWH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define SIMULATECOHERENTGWH_ENUL 1
#define SIMULATECOHERENTGWH_EBAD 2
#define SIMULATECOHERENTGWH_ESIG 3
#define SIMULATECOHERENTGWH_EDIM 4

#define SIMULATECOHERENTGWH_MSGENUL "Unexpected null pointer in arguments"
#define SIMULATECOHERENTGWH_MSGEBAD "A sampling interval is (effectively) zero"
#define SIMULATECOHERENTGWH_MSGESIG "Input signal must specify amplitude and phase functions"
#define SIMULATECOHERENTGWH_MSGEDIM "Amplitude and phase must be 1- or 2-dimensional vectors"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Structures}
\begin{verbatim}
CoherentGW
\end{verbatim}
\index{\texttt{CoherentGW}}

\noindent This structure stores a representation of a plane
gravitational wave propagating from a particular point on the sky.
Several alternate representations are permitted to allow a more
natural characterization of coherent waveforms.  The fields are:

\begin{description}
\item[\texttt{REAL4 ra, dec}] The right ascension and declination of
the point from which the wave is propagating, in radians.

\item[\texttt{REAL4 psi}] The polarization angle $\psi$, in radians,
as defined in Appendix~B of~\cite{Anderson_W:2000}.

\item[\texttt{REAL4TimeVectorSeries *h}] A time-sampled
two-dimensional vector storing the waveforms $h_+(t)$ and
$h_\times(t)$.

\item[\texttt{REAL4TimeVectorSeries *a}] A time-sampled
two-dimensional vector storing the amplitudes $A_+(t)$ and
$A_\times(t)$.  Optionally, the vector can be one-dimensional, in
which case it is assumed that the polarizations have equal amplitude.

\item[\texttt{REAL4TimeVectorSeries *phi}] A time-sampled
two-dimensional vector storing the phase functions $\phi_+(t)$ and
$\phi_\times(t)$, in radians.  Optionally, the vector can be
one-dimensional, in which case it is assumed that waveform is
circularly polarized with $\phi_{\times}=\phi_{+}-\pi/2$.
\end{description}

\noindent It is permitted to set \verb@*h@ while leaving
\verb@a=phi=NULL@, or to set \verb@*a@ and \verb@*phi@ while leaving
\verb@h=NULL@, but if all three are specified then they must be
consistent over any time period of overlap.  Furthermore, \verb@*a@
and \verb@*phi@ are only meaningful over time periods where
\emph{both} are defined.

******************************************************* </lalLaTeX> */

typedef struct tagCoherentGW {
  REAL4 ra, dec;              /* sky position of source */
  REAL4 psi;                  /* polarization angle of source */
  REAL4TimeVectorSeries *h;   /* sampled waveforms h_+, h_x */
  REAL4TimeVectorSeries *a;   /* amplitudes A_+, A_x */
  REAL4TimeVectorSeries *phi; /* phase functions phi_+, phi_x */
} CoherentGW;

/********************************************************** <lalLaTeX>

\begin{verbatim}
DetectorResponse
\end{verbatim}
\index{\texttt{CoherentGW}}

\noindent This structure contains information required to determine
the response of a detector to a gravitational waveform.  The fields
are:

\begin{description}
\item[\texttt{COMPLEX8FrequencySeries *transfer}] The
frequency-dependent transfer function of the interferometer, in ADC
counts per unit strain amplitude at any given frequency.

\item[\texttt{??? *detector}] A structure storing site and
polarization information, used to compute the polarization response at
any time.  This structure has yet to be defined.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagDetectorResponse {
  COMPLEX8FrequencySeries *transfer; /* frequency transfer function */
  LALDetector *site;           /* detector location and orientation */
} DetectorResponse;


/* <lalLaTeX>
\vfill{\footnotesize\input{SimulateCoherentGWHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{SimulateCoherentGWC}
</lalLaTeX> */
void
LALSimulateCoherentGW( LALStatus        *stat,
		       REAL4TimeSeries  *output,
		       CoherentGW       *signal,
		       DetectorResponse *detector );

/* <lalLaTeX>
%\newpage\input{SimulateCoherentGWTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _SIMULATECOHERENTGW_H */
