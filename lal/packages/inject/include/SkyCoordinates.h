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

%\begin{wrapfigure}{r}{0.55\textwidth}
\begin{figure}[bp]
\begin{center}
\resizebox{0.5\textwidth}{!}{\includegraphics{inject_lat_long}} \\
\parbox{0.5\textwidth}{\caption{\label{fig:lat-long} Definition of the
latitude $\phi$ and longitude $\lambda$ of a point $P$, in an
arbitrary coordinate system specified by an origin $O$ and axes $z$
and $x$, as shown.  These routines generally assume a fixed $O$, so
the distance coordinate $r$ is not transformed.  (Caution: The
characters $\phi$ and $\lambda$ do not seem to print in hardcopy.)}}
\end{center}
%\end{wrapfigure}
\end{figure}
A general spatial coordinate system is defined by six parameters:
three positions specifying the location of the origin, and three
angles specifying the orientations of the coordinate axes.  In
astronomy it is normally assumed that the centre of the coordinate
system is the centre of the Earth (geocentric coordinates).  Once the
origin has been specified, the orientation is fixed by defining one
direction as the \emph{pole} or $z$-axis (two degrees of freedom), and
another orthogonal direction as the \emph{reference meridian} or
$x$-axis (one degree of freedom).  A $y$-axis can also be defined such
that $x$, $y$, and $z$ form an orthogonal right-handed coordinate
system; however, in astronomy it is more conventional to use spherical
coordinates, defined as follows:

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

In the routines in this module, we do not perform transformations
between coordinate systems having different origins.  By default, all
coordinates are assumed to be centred on the observer; however, one
may also consider coordinate systems that are \emph{geogentric}
(having their origin at the centre of the Earth), \emph{heliocentric}
(origin at the centre of the Sun), \emph{barycentric} (origin at the
centre of mass of the solar system), and \emph{Galactocentric} (origin
at the centre of our Galaxy).  Since we ignore translations in the
coordinate origin, distances remain unchanged, so these routines only
consider transformations in latitude and longitude.  To put it another
way, these routines transform \emph{directions} in space, not
\emph{locations} in space.  These directions are generically stored in
the \verb@SkyPosition@ structure, defined below.

The coordinate systems that we consider are defined as follows:

\paragraph{Horizon coordinates:} This is defined assuming an observer
at a particular point on the Earth.  The $z$-axis is defined to be the
direction opposite to the local acceleration due to gravity.  The
$x$-axis is defined to lie in the plane formed by the $z$-axis and the
Earth's rotational axis, and to be directed into the northern
hemisphere.  In this coordinate system, the latitude coordinate is
called the \emph{altitude} and the longitude coordinate is the
\emph{negative} of what astronomers call the \emph{azimuth}; this sign
reversal is due to the fact that astronomers define azimuth to
increase clockwise, and our longitudinal coordinates uniformly
increase counterclockwise about the $z$-axis.

\paragraph{Geographic coordinates:} The $z$-axis is defined to be
parallel to the Earth's axis, in the direction of the Earth's north
pole.  The $x$-axis is defined to be parallel to the direction
perpendicular from the Earth's rotation axis to a reference point in
Greenwich, UK.  Note that we adopt a longitude convention that is
consistent with the \textit{Astronomical Almanac}, but opposite to
that in~\cite{Lang_K:1998}, in that our geographic longitudes increase
\emph{eastward} (counterclockwise) like the rest of our longitudinal
coordinates.

Geographic latitude and longitude are often referred to simply as
latitude and longitude, and are represented in~\cite{Lang_K:1998} by
the symbols $\lambda$ and $\phi$, as in Fig.~\ref{fig:lat-long}.
However, we emphasize once again that geographic latitude and
longitude as defined above refer to directions in space, not to
locations on the Earth's surface.  This can lead to some confusion.
The \emph{geodetic} latitude and longitude of a point on the Earth's
surface are the latitude and longitude of its vertical direction,
while the \emph{geocentric} latitude and longitude of the point are
the latitude and longitude of the line from the geometric centre of
the Earth through that point.  These are not necessarily the same, due
to the Earth's ellipticity.

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
\newpage\input{TerrestrialCoordinatesC}
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
			SkyPosition *position,
			SkyPosition *zenith );

void
LALHorizonToGeographic( LALStatus   *stat,
			SkyPosition *position,
			SkyPosition *zenith );

/* <lalLaTeX>
%\newpage\input{SkyCoordinatesTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _SKYCOORDINATES_H */
