/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Reinhard Prix
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

/************************* <lalVerbatim file="CelestialCoordinatesCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{CelestialCoordinates.c}}
\label{ss:CelestialCoordinates.c}

Converts among Galactic, ecliptic, and equatorial coordinates.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{CelestialCoordinatesCP}
\idx{LALGalacticToEquatorial()}
\idx{LALEquatorialToGalactic()}
\idx{LALEclipticToEquatorial()}
\idx{LALEquatorialToEcliptic()}

\subsubsection*{Description}

These functions perform the specified coordinate transformation on the
contents of \verb@*input@ and store the result in \verb@*output@.  The
two pointers may point to the same object, in which case the
conversion is done in place.  The functions will also check
\verb@input->system@ and set \verb@output->system@ as appropriate.

These routines are collected together because they involve fixed,
absolute coordinate systems, so the transformations require no
additional parameters such as the time or site of observation.  We
also note that there are no direct conversions between Galactic and
ecliptic coordinates.  At the risk of additional computational
overhead, it is simple to use the equatorial coordinate system as an
intermediate step.

\subsubsection*{Algorithm}

These routines follow the spherical angle relations on p.~13
of~\cite{Lang_K:1999}.  Note that the actual formulae for Galactic
longitude and right ascension in this reference are wrong; we give
corrected formulae below derived from the sine and cosine equations.
(The Galactic to equatorial transformations can also be found in
Sec.~12.3 of~\cite{GRASP_1.9.8:2000}.)  All positions are assumed to
be in the J2000 epoch.

\paragraph{Galactic coordinates:} The following formulae relate
Galactic latitude $b$ and longitude $l$ to declination $\delta$ and
right ascension $\alpha$:
\begin{eqnarray}
\label{eq:b-galactic}
b & = & \arcsin[\cos\delta\cos\delta_\mathrm{NGP}
		\cos(\alpha-\alpha_\mathrm{NGP}) +
		\sin\delta\sin\delta_\mathrm{NGP}] \;,\\
l & = & \arctan\!2[\sin\delta\cos\delta_\mathrm{NGP} -
		\cos\delta\cos(\alpha-\alpha_\mathrm{NGP})
			\sin\delta_\mathrm{NGP},
		\cos\delta\sin(\alpha-\alpha_\mathrm{NGP})] \nonumber\\
\label{eq:l-galactic}
& & \quad + \; l_\mathrm{ascend} \;,
\end{eqnarray}
where $\arctan\!2(y,x)$ can be thought of as the argument of the
complex number $x+iy$; unlike $\arctan(y/x)$, it ranges over the full
range $[0,2\pi)$ instead of just half of it.  The inverse
transformations are:
\begin{eqnarray}
\label{eq:delta-galactic}
\delta & = & \arcsin[\cos b\cos\delta_\mathrm{NGP}\sin(l-l_\mathrm{ascend}) +
		\sin b\sin\delta_\mathrm{NGP}] \;,\\
\alpha & = & \arctan\!2[\cos b\cos(l-l_\mathrm{ascend}),
		\sin b\cos\delta_\mathrm{NGP} -
		\cos b\sin(l-l_\mathrm{ascend})\sin\delta_\mathrm{NGP}]
		\nonumber\\
\label{eq:alpha-galactic}
& & \quad + \; \alpha_\mathrm{NGP} \;.
\end{eqnarray}
In these equations we have defined the orientation of the Galaxy with
the following parameters (which should eventually be placed in
\verb@LALConstants.h@):
$$
\begin{array}{r@{\quad=\quad}l@{\quad=\quad}l}
\alpha_\mathrm{NGP} & 192.8594813^\circ &
\mbox{the right ascension (epoch J2000) of the north Galactic pole} \\
\delta_\mathrm{NGP} & 27.1282511^\circ &
\mbox{the declination (epoch J2000) of the north Galactic pole} \\
l_\mathrm{ascend} & 33^\circ &
\mbox{the longitude of the ascending node of the Galactic plane}
\end{array}
$$
The ascending node of the Galactic plane is defined as the direction
along the intersection of the Galactic and equatorial planes where
rotation in the positive sense about the Galactic $z$ axis carries a
point from the southern to northern equatorial hemisphere.  That is,
if \textbf{\textit{u}} points in the direction $\delta=90^\circ$
(celestial north), and \textbf{\textit{v}} points in the direction
$b=90^\circ$ (Galactic north), then
\textbf{\textit{u}}$\times$\textbf{\textit{v}} points along the
ascending node.

\paragraph{Ecliptic coordinates:} The following formulae relate
Ecliptic latitude $\beta$ and longitude $\lambda$ to declination
$\delta$ and right ascension $\alpha$:
\begin{eqnarray}
\label{eq:beta-ecliptic}
\beta & = & \arcsin(\sin\delta\cos\epsilon -
		\cos\delta\sin\alpha\sin\epsilon) \;, \\
\label{eq:l-ecliptic}
\lambda & = & \arctan\!2(\cos\delta\sin\alpha\cos\epsilon +
		\sin\delta\sin\epsilon, \cos\delta\cos\alpha) \;.
\end{eqnarray}
The inverse transformations are:
\begin{eqnarray}
\label{eq:delta-ecliptic}
\delta & = & \arcsin(\cos\beta\sin\lambda\sin\epsilon +
		\sin\beta\cos\epsilon) \;, \\
\label{eq:alpha-ecliptic}
\alpha & = & \arctan\!2(\cos\beta\sin\lambda\cos\epsilon -
		\sin\beta\sin\epsilon, \cos\beta\cos\lambda) \;.
\end{eqnarray}
Here $\epsilon$ is the obliquity (inclination) of the ecliptic plane,
which varies over time; at epoch J200 it has a mean value of:
$$
\epsilon = 23.4392911^\circ \; .
$$

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{CelestialCoordinatesCV}}

******************************************************* </lalLaTeX> */

/* <lalVerbatim file="CelestialCoordinatesCP">
void
LALGalacticToEquatorial( LALStatus   *stat,
			 SkyPosition *output,
			 SkyPosition *input )
{ </lalVerbatim> */

/* <lalVerbatim file="CelestialCoordinatesCP">
void
LALEquatorialToGalactic( LALStatus   *stat,
			 SkyPosition *output,
			 SkyPosition *input )
{ </lalVerbatim> */

/* <lalVerbatim file="CelestialCoordinatesCP">
void
LALEclipticToEquatorial( LALStatus   *stat,
			 SkyPosition *output,
			 SkyPosition *input )
{ </lalVerbatim> */

/* <lalVerbatim file="CelestialCoordinatesCP">
void
LALEquatorialToEcliptic( LALStatus   *stat,
			 SkyPosition *output,
			 SkyPosition *input )
{ </lalVerbatim> */

#include <stddef.h> /* just so the contents aren't empty... */

#include <lal/LALRCSID.h>
NRCSID (CELESTIALCOORDINATESLALDOCC,"$Id$");
