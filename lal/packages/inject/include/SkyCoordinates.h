/******************************* <lalVerbatim file="SkyCoordinatesHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{SkyCoordinates.h}}
\label{s:SkyCoordinates.h}

Provides routines to convert among various sky coordinate systems.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/SkyCoordinates.h>
\end{verbatim}

This header covers routines to perform coordinate transformations
among the various spherical coordinate systems used in astronomy.
Most of these routines are discussed in Sec.~5.1
of~\cite{Lang_K:1998}; we reproduce here some of the essential
elements of this discussion.

A general spatial coordinate system is defined by six parameters:
three positions specifying the location of the origin, and three
angles specifying the orientations of the coordinate axes.  In
astronomy it is normally assumed that the centre of the coordinate
system is the centre of the Earth (geocentric coordinates).  Once the
origin has been specified, the orientation is fixed by defining one
direction as the \emph{north pole} or $z$-axis (two degrees of
freedom), and another orthogonal direction as the \emph{reference
meridian} or $x$-axis (one degree of freedom).  A $y$-axis can also be
defined such that $x$, $y$, and $z$ form an orthogonal right-handed
coordinate system; however, in astronomy it is more conventional to
use spherical coordinates, defined as follows:

For any given point, define a plane (called its meridian) containing
both the $z$-axis and the point in question.  The \emph{longitude} is
the angle in the $x$-$y$ plane from the reference meridian to the line
where the object's meridian crosses the $x$-$y$ plane.  The
\emph{latitude} is the angle in the meridian plane between this line
and the direction to the object.  The \emph{distance} is simply
measured in a straight line from the origin to the point.  Longitude
is defined to increase in the right-handed direction about the
$z$-axis (i.e.\ the $y$-axis lies at \emph{positive} $\pi/2$~radians
longitude), and is typically given in the range $[0,2\pi)$~radians.
Latitude is defined to increase towards the $z$-axis (i.e.\ the
$z$-axis lies at \emph{positive} $\pi/2$~radians latitude), and is
typically given in the range $[-\pi/2,\pi/2]$~radians; a point with
latitude $\pm\pi/2$~radians has arbitrary (undefined) longitude.
Distance should always be positive.  This convention is shown in
Fig.~\ref{fig:lat-long}.
\begin{figure}
\centerline{\resizebox{0.8\textwidth}{!}{%
\includegraphics{inject_lat_long}}}
\caption{\label{fig:lat-long} Definition of latitude $\beta$ and
longitude $\lambda$ for a point $P$ in an arbitrary coordinate system.
The point $O$ is the origin of the coordinate system, with $z$ and $x$
axes as shown.  Both $\beta$ and $\lambda$ are drawn in the positive
sense.  The distance $r$ to the point is not used by any of the
routines under this header, since the origin of the coordinate system
is taken to be fixed.}
\end{figure}

In the routines in this module, we do not perform transformations
between coordinate systems having different origins.  By default, all
coordinates are assumed to be \emph{geogentric} (having their origin
at the centre of the Earth), unless otherwise specified.  Other common
specifications for the coordinate origin are \emph{heliocentric}
(origin at the centre of the Sun), \emph{barycentric} (origin at the
centre of mass of the solar system), and \emph{Galactocentric} (origin
at the centre of our Galaxy).  Since we ignore translations in the
coordinate origin, distances remain unchanged, so these routines only
consider transformations in latitude and longitude.  These are
generically stored in the \verb@SkyPosition@ structure, defined below.

The coordinate systems that we consider are defined as follows:

\paragraph{Horizon coordinates:} This is defined assuming an observer
at a particular point on the Earth.  The $z$-axis is defined to be the
direction opposite to the local acceleration due to gravity.  The
$x$-axis is defined to lie in the plane formed by the $z$-axis and the
Earth's rotational axis, and to be directed into the northern
hemisphere; if the $z$-axis and rotation axis are parallel, the
$x$-axis is defined to be the same as for geographic coordinates
(below).  In this coordinate system, the latitude coordinate is called
the \emph{altitude} and the longitude coordinate is called the
\emph{azimuth}.

\paragraph{Geographic coordinates:} The $z$-axis is defined to be
parallel to the Earth's axis, in the direction of the Earth's north
pole.  The $x$-axis is defined to be parallel to the direction
perpendicular from the Earth's rotation axis to a reference point in
Greenwich, UK.  Note that we adopt a longitude convention that is
opposite to that in~\cite{Lang_K:1998}, in that our geographic
longitudes increase \emph{eastward} like the rest of our longitudinal
coordinates.

\paragraph{Equatorial coordinates:} The $z$-axis is defined as for
geographic coordinates, above; the plane orthogonal to this passing
through the Earth's centre is called the \emph{equator}.  The $x$-axis
is defined to be the direction, as viewed from the centre of the
Earth, where the Sun appears to cross the equator moving north in
spring.  In this coordinate system, the latitude coordinate is called
the \emph{declination} $\delta$ and the longitude coordinate is called
the \emph{right ascension} $\alpha$.

\paragraph{Ecliptic coordinates:} The $z$-axis is defined to be the
direction orthogonal to the orbital plane of the Earth about the Sun,
directed such that the Earth orbits in a right-handed sense.  The
$x$-axis is defined as for equatorial coordinates, above; we note that
by definition it lies parallel to the intersection of the equatorial
and orbital planes of the Earth.  The ecliptic latitude is normally
denoted as $\beta$ and the ecliptic longitude as $\lambda$.

\paragraph{Galactic coordinates:} The $z$-axis is defined to be the
direction orthogonal to the plane of our Galaxy and pointing into the
northern hemisphere of the equatorial coordinate system.
(Unfortunately this convention has the result that Galactic rotation
is left-handed about this axis.)  The $x$-axis is defined to be the
direction of the Galactic centre as viewed from the Earth.  The
Galactic latitude coordinate is normally denoted as $b$ and the
Galactic longitude as $l$.

******************************************************* </lalLaTeX> */

#ifndef _SKYCOORDINATES_H
#define _SKYCOORDINATES_H

#include <lal/LALStdlib.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID( SKYCOORDINATESH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define SKYCOORDINATESH_ENUL  1

#define SKYCOORDINATESH_MSGENUL  "Unexpected null pointer in arguments"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Structures}
\begin{verbatim}
SkyPosition
\end{verbatim}
\index{\texttt{SkyPosition}}

\noindent This structure stores the two spherical coordinates of a sky
position; i.e.\ a generic latitude and longitude.  The structure is
not defined specific to a particular coordinate system; it is up to
the defining routines to specify whether the coordinates are, say, an
ecliptic latitude and longitude, or a declination and a right
ascension.  The fields are:

\begin{description}
\item[\texttt{REAL8 longitude}] The longitudinal coordinate, as
defined above.

\item[\texttt{REAL8 latitude}] The latitudinal coordinate, as defined
above.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagSkyPosition {
  REAL8 longitude;
  REAL8 latitude;
} SkyPosition;

/* <lalLaTeX>
\vfill{\footnotesize\input{SkyCoordinatesHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{CelestialCoordinatesC}
</lalLaTeX> */
void
LALGalacticToEquatorial( LALStatus   *stat,
			 SkyPosition *position );

void
LALEquatorialToGalactic( LALStatus   *stat,
			 SkyPosition *position );

void
LALEclipticToEquatorial( LALStatus   *stat,
			 SkyPosition *position );

void
LALEquatorialToEcliptic( LALStatus   *stat,
			 SkyPosition *position );

/* <lalLaTeX>
%\newpage\input{TerrestrialCoordinatesC}
</lalLaTeX> */
void
LALGeographicToEquatorial( LALStatus   *stat,
			   SkyPosition *position,
			   LIGOTimeGPS *time );

void
LALEquatorialToGeographic( LALStatus   *stat,
			   SkyPosition *position,
			   LIGOTimeGPS *time );

void
LALGeographicToHorizon( LALStatus   *stat,
			SkyPosition *position );
  /* plus some site parameter */

void
LALHorizonToGeographic( LALStatus   *stat,
			SkyPosition *position );
  /* plus some site parameter */

/* <lalLaTeX>
%\newpage\input{SkyCoordinatesTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _SKYCOORDINATES_H */
