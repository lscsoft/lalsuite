/*
*  Copyright (C) 2007 Jolien Creighton, Teviet Creighton
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

/***************************** <lalVerbatim file="GenerateTaylorCWHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{GenerateTaylorCW.h}}
\label{s:GenerateTaylorCW.h}

Provides routines to generate Taylor-parameterized continuous
waveforms.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/GenerateTaylorCW.h>
\end{verbatim}

This header covers routines to generate continuous quasiperiodic
waveforms whose frequency varies slowly and smoothly with time.  For
such sources the frequency function is normally described by its
Taylor ``spindown'' (or spin-up) coefficients.  This type of waveform
may be typical of objects such as neutron stars that are gradually
shedding angular momentum, or are accelerating in the gravitational
potential of a star cluster.  The Taylor expansion is likely
\emph{not} suitable for long-term modelling of the frequency of waves
from glitching neutron stars, neutron stars in close binary orbits, or
neutron stars that are accreting or shedding angular momentum in a
stochastic manner.

The frequency and phase of such slowly-varying quasiperiodic sources
are given by their Taylor series:
\begin{eqnarray}
\label{eq:taylorcw-freq}
f(t)    & = & f_0 \left[ 1 + \sum_{k=1}^n f_k(t-t_0)^k \right] \;, \\
\label{eq:taylorcw-phi}
\phi(t) & = & \phi_0 + 2\pi f_0 \left[ (t-t_0) +
		\sum_{k=1}^n \frac{f_k}{k+1}(t-t_0)^{k+1} \right] \;,
\end{eqnarray}
where $f_k$ are the spin-normalized Taylor coefficients.  If the
source's spin is varying over some timescale $\tau$, one typically
expects that $f_k\sim\tau^{-k}$.  Note that in this and later
discussions, $f$ and $\phi$ refer to the frequency and phase of the
gravitational wave, which are typically some constant multiple of
(often twice) the frequency and phase of the rotating source.

The \verb@CoherentGW@ structure allows for a very general
description of waveforms with modulations in the amplitudes or
relative phases of the wave polarizations, as described in
\verb@SimulateCoherentGW.h@.  However, in this simplest model of
quasiperiodic waveforms, we neglect such phenomena as precession that
would produce these effects.  Thus for any given source one can choose
a polarization basis (described by some polarization angle $\psi$) in
which the wave has a constant elliptical polarization of the form:
\begin{eqnarray}
\label{eq:taylorcw-hplus}
h_+(t)      & = & A_+      \cos\phi(t) \;, \\
\label{eq:taylorcw-hcross}
h_\times(t) & = & A_\times \sin\phi(t) \;.
\end{eqnarray}

******************************************************* </lalLaTeX> */

#ifndef _GENERATETAYLORCW_H
#define _GENERATETAYLORCW_H

#include <lal/LALStdlib.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/SkyCoordinates.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

NRCSID( GENERATETAYLORCWH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define GENERATETAYLORCWH_ENUL 1
#define GENERATETAYLORCWH_EOUT 2
#define GENERATETAYLORCWH_EMEM 3

#define GENERATETAYLORCWH_MSGENUL "Unexpected null pointer in arguments"
#define GENERATETAYLORCWH_MSGEOUT "Output field a, f, phi, or shift already exists"
#define GENERATETAYLORCWH_MSGEMEM "Out of memory"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{TaylorCWParamStruc}}
\idx[Type]{TaylorCWParamStruc}

This structure stores the parameters for constructing a gravitational
waveform with a Taylor-polynomial frequency and phase.  As with the
\verb@PPNParamStruc@ type in \verb@GeneratePPNInspiral.h@, we divide
the fields into passed fields (which are supplied to the final
\verb@CoherentGW@ structure but not used in any calculations), input
fields (that are used by the waveform generator), and output fields
(that are set by the waveform generator).  They are:

\bigskip\noindent\textit{Passed fields:}
\begin{description}
\item[\texttt{SkyPosition position}] The location of the source on the
sky, normally in equatorial coordinates.

\item[\texttt{REAL4 psi}] The polarization angle of the source, in
radians.

\item[\texttt{LIGOTimeGPS epoch}] The start time $t_0$ of the output
series.
\end{description}

\medskip\noindent\textit{Input fields:}
\begin{description}
\item[\texttt{REAL8 deltaT}] The requested sampling interval of the
waveform, in s.

\item[\texttt{UINT4 length}] The number of samples in the generated
waveform.

\item[\texttt{REAL4 aPlus, aCross}] The polarization amplitudes $A_+$,
$A_\times$, in dimensionless strain units.

\item[\texttt{REAL8 phi0}] The wave phase at time $t_0$, in radians.

\item[\texttt{REAL8 f0}] The wave frequency at time $t_0$, in Hz.

\item[\texttt{REAL8Vector *f}] The spin-normalized Taylor parameters
$f_k$, as defined in Eq.~\ref{eq:taylorcw-freq}, above.  If
\verb@f@=\verb@NULL@, a monochromatic wave is generated.
\end{description}

\medskip\noindent\textit{Output fields:}
\begin{description}
\item[\texttt{REAL4 dfdt}] The maximum value of $\Delta f\Delta t$
encountered over any timestep $\Delta t$ used in generating the
waveform.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagTaylorCWParamStruc {
  /* Passed parameters. */
  SkyPosition position; /* location of source on sky */
  REAL4 psi;            /* polarization angle (radians) */
  LIGOTimeGPS epoch;    /* start time of output time series */

  /* Input parameters. */
  REAL8 deltaT;         /* requested sampling interval (s) */
  UINT4 length;         /* length of time series */
  REAL4 aPlus, aCross;  /* polarization amplitudes */
  REAL8 phi0;           /* initial phase */
  REAL8 f0;             /* initial frequency */
  REAL8Vector *f;       /* f0-normalized Taylor parameters */

  /* Output parameters. */
  REAL4 dfdt;           /* maximum value of df*dt over any timestep */
} TaylorCWParamStruc;


/* <lalLaTeX>
\vfill{\footnotesize\input{GenerateTaylorCWHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{GenerateTaylorCWC}
</lalLaTeX> */
void
LALGenerateTaylorCW( LALStatus          *,
		     CoherentGW         *output,
		     TaylorCWParamStruc *params );

/* <lalLaTeX>
\newpage\input{SimulateTaylorCWTestC}
</lalLaTeX> */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GENERATETAYLORCW_H */
