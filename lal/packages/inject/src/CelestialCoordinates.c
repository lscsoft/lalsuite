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
\index{\texttt{LALGalacticToEquatorial()}}
\index{\texttt{LALEquatorialToGalactic()}}
\index{\texttt{LALEclipticToEquatorial()}}
\index{\texttt{LALEquatorialToEcliptic()}}

\subsubsection*{Description}

These functions convert celestial coordinates from one spherical
coordinate system to another; the definitions of the coordinate
systems are given in \verb@SkyCoordinates.h@.  For simplicity the
transformation is done in place; \verb@*position@ stores the
coordinates in the first coordinate system at the start of the call,
and the second coordinate system upon return.

These routines are collected together because they involve fixed,
absolute coordinate systems, so the transformations require no
additional parameters such as the time or site of observation.  We
also note that there are no direct conversions between Galactic and
ecliptic coordinates.  At the risk of additional computational
overhead, it is simple to use the equatorial coordinate system as an
intermediate step.

\subsubsection*{Algorithm}

These routines follow the formulae on p.~13 of~\cite{Lang_K:1998},
which we reproduce below; these are also found in Sec.~12.3
of~\cite{GRASP_1.9.8:2000}.  All positions are assumed to be in the
J2000 epoch.

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
		\cos b\sin(l-l_\mathrm{ascend})
		- \sin b\cos\delta_\mathrm{NGP}] \nonumber\\
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

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SkyCoordinates.h>

#define LAL_ALPHAGAL 3.366032942
#define LAL_DELTAGAL 0.473477302
#define LAL_LGAL 0.576

NRCSID( CELESTIALCOORDINATESC, "$Id$" );

/* <lalVerbatim file="CelestialCoordinatesCP"> */
void
LALGalacticToEquatorial( LALStatus   *stat,
			 SkyPosition *position )
{ /* </lalVerbatim> */
  REAL8 sinDGal = sin( LAL_DELTAGAL ); /* sin(delta_NGP) */
  REAL8 cosDGal = cos( LAL_DELTAGAL ); /* cos(delta_NGP) */
  REAL8 l = -LAL_LGAL;          /* will be l-l(ascend) */
  REAL8 sinB, cosB, sinL, cosL; /* sin and cos of b and l */
  REAL8 sinD, sinA, cosA;       /* sin and cos of delta and alpha */

  INITSTATUS( stat, "LALGalacticToEquatorial", CELESTIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( position, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Compute intermediates. */
  l += position->longitude;
  sinB = sin( position->latitude );
  cosB = cos( position->latitude );
  sinL = sin( l );
  cosL = cos( l );

  /* Compute components. */
  sinD = cosB*cosDGal*sinL + sinB*sinDGal;
  sinA = cosB*cosL;
  cosA = cosB*sinL - sinB*cosDGal;

  /* Compute final results. */
  position->latitude = asin( sinD );
  l = atan2( sinA, cosA ) + LAL_ALPHAGAL;

  /* Optional phase correction. */
  if ( l < 0.0 )
    l += LAL_TWOPI;
  position->longitude = l;

  RETURN( stat );
}


/* <lalVerbatim file="CelestialCoordinatesCP"> */
void
LALEquatorialToGalactic( LALStatus   *stat,
			 SkyPosition *position )
{ /* </lalVerbatim> */
  REAL8 sinDGal = sin( LAL_DELTAGAL ); /* sin(delta_NGP) */
  REAL8 cosDGal = cos( LAL_DELTAGAL ); /* cos(delta_NGP) */
  REAL8 a = -LAL_ALPHAGAL;      /* will be alpha-alpha_NGP */
  REAL8 sinD, cosD, sinA, cosA; /* sin and cos of delta and alpha */
  REAL8 sinB, sinL, cosL;       /* sin and cos of b and l */

  INITSTATUS( stat, "LALGalacticToEquatorial", CELESTIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( position, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Compute intermediates. */
  a += position->longitude;
  sinD = sin( position->latitude );
  cosD = cos( position->latitude );
  sinA = sin( a );
  cosA = cos( a );

  /* Compute components. */
  sinB = cosD*cosDGal*cosA + sinD*sinDGal;
  sinL = sinD*cosDGal - cosD*cosA*sinDGal;
  cosL = cosD*sinA;

  /* Compute final results. */
  position->latitude = asin( sinB );
  a = atan2( sinL, cosL ) + LAL_LGAL;

  /* Optional phase correction. */
  if ( a < 0.0 )
    a += LAL_TWOPI;
  position->longitude = a;

  RETURN( stat );
}


/* <lalVerbatim file="CelestialCoordinatesCP"> */
void
LALEclipticToEquatorial( LALStatus   *stat,
			 SkyPosition *position )
{ /* </lalVerbatim> */
  REAL8 sinE = sin( LAL_IEARTH ); /* sin(epsilon) */
  REAL8 cosE = cos( LAL_IEARTH ); /* cos(epsilon) */
  REAL8 sinB, cosB, sinL, cosL;   /* sin and cos of b and l */
  REAL8 sinD, sinA, cosA;         /* sin and cos of delta and alpha */

  INITSTATUS( stat, "LALGalacticToEquatorial", CELESTIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( position, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Compute intermediates. */
  sinB = sin( position->latitude );
  cosB = cos( position->latitude );
  sinL = sin( position->longitude );
  cosL = cos( position->longitude );

  /* Compute components. */
  sinD = cosB*sinL*sinE + sinB*cosE;
  sinA = cosB*sinL*cosE - sinB*sinE;
  cosA = cosB*cosL;

  /* Compute final results. */
  position->latitude = asin( sinD );
  position->longitude = atan2( sinA, cosA );

  /* Optional phase correction. */
  if ( position->longitude < 0.0 )
    position->longitude += LAL_TWOPI;

  RETURN( stat );
}


/* <lalVerbatim file="CelestialCoordinatesCP"> */
void
LALEquatorialToEcliptic( LALStatus   *stat,
			 SkyPosition *position )
{ /* </lalVerbatim> */
  REAL8 sinE = sin( LAL_IEARTH ); /* sin(epsilon) */
  REAL8 cosE = cos( LAL_IEARTH ); /* cos(epsilon) */
  REAL8 sinD, cosD, sinA, cosA;   /* sin and cos of delta and alpha */
  REAL8 sinB, sinL, cosL;         /* sin and cos of b and l */

  INITSTATUS( stat, "LALGalacticToEquatorial", CELESTIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( position, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Compute intermediates. */
  sinD = sin( position->latitude );
  cosD = cos( position->latitude );
  sinA = sin( position->longitude );
  cosA = cos( position->longitude );

  /* Compute components. */
  sinB = sinD*cosE - cosD*sinA*sinE;
  sinL = cosD*sinA*cosE + sinD*sinE;
  cosL = cosD*cosA;

  /* Compute final results. */
  position->latitude = asin( sinB );
  position->longitude = atan2( sinL, cosL );

  /* Optional phase correction. */
  if ( position->longitude < 0.0 )
    position->longitude += LAL_TWOPI;

  RETURN( stat );
}
