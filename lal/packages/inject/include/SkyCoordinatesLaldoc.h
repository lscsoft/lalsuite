/*
*  Copyright (C) 2007 Bernd Machenschalk, Reinhard Prix
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

/******************************* <lalVerbatim file="SkyCoordinatesLaldocHV">
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
of~\cite{Lang_K:1999}; we reproduce here some of the essential
elements of this discussion.

\begin{wrapfigure}{r}{0.55\textwidth}
\vspace{-2ex}
\begin{center}
\resizebox{0.5\textwidth}{!}{\includegraphics{inject_lat_long}} \\
\parbox{0.5\textwidth}{\caption{\label{fig:lat-long} Definition of
latitude $\phi$ and longitude $\lambda$ for an arbitrary coordinate
system.}}
\end{center}
\vspace{-2ex}
\end{wrapfigure}
A general spatial coordinate system is shown in
Fig.~\ref{fig:lat-long}.  It is defined by six parameters: three
positions specifying the location of the origin $O$, two angles
specifying the direction of the \emph{pole} or $z$-axis, and one
further angle specifying a \emph{reference meridian} or $x$-axis
orthogonal to the $z$-axis.  A $y$-axis can also be defined such that
$x$, $y$, and $z$ form an orthogonal right-handed coordinate system;
however, in astronomy it is more conventional to use spherical
coordinates, defined as follows:

For any given point $P$, define a plane (called its meridian)
containing both the $z$-axis and the point in question.  The
\emph{longitude} $\lambda$ is the angle in the $x$-$y$ plane from the
$x$-axis to the line where the object's meridian crosses the $x$-$y$
plane.  The \emph{latitude} $\phi$ is the angle in the meridian plane
between this line and the direction to the object.  The
\emph{distance} $r$ is simply measured in a straight line from the
origin to the point.  Longitude is defined to increase in the
right-handed direction about the $z$-axis (i.e.\ the $y$-axis lies at
\emph{positive} $\pi/2$~radians longitude), and is typically given in
the range $[0,2\pi)$~radians.  Latitude is defined to increase towards
the $z$-axis (i.e.\ the $z$-axis lies at \emph{positive}
$\pi/2$~radians latitude), and is typically given in the range
$[-\pi/2,\pi/2]$~radians; a point with latitude $\pm\pi/2$~radians has
arbitrary (undefined) longitude.  Distance should always be positive.
This convention is shown in Fig.~\ref{fig:lat-long}.

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


\newpage
\begin{wrapfigure}{r}{0.63\textwidth}
\vspace{-2ex}
\begin{center}
\resizebox{0.58\textwidth}{!}{\includegraphics{inject_horizon}} \\
\parbox{0.58\textwidth}{\caption{\label{fig:horizon} Definition of the
horizon coordinate system.}}
\end{center}
\vspace{-2ex}
\end{wrapfigure}
\paragraph{Horizon coordinates:} This is a local coordinate system for
a particular observation point $O$ on the Earth, as shown in
Fig.~\ref{fig:horizon}.  The $z$-axis is defined to be the direction
opposite to the local acceleration due to gravity.  The $x$-axis is
defined to lie in the plane formed by the $z$-axis and the Earth's
rotational axis, and to be directed into the northern hemisphere.  In
this coordinate system, the latitude coordinate is called the
\emph{altitude} and the longitude coordinate is the \emph{negative} of
what astronomers call the \emph{azimuth}; this sign reversal is due to
the fact that astronomers define azimuth to increase clockwise, and
our longitudinal coordinates uniformly increase counterclockwise about
the $z$-axis.

This coordinate system is related to the geographic coordinate system
(below) by the geographic latitude $\phi_z$ and longitude $\lambda_z$
of the observer's $z$-axis direction.

\begin{wrapfigure}{l}{0.43\textwidth}
\vspace{-4ex}
\begin{center}
\resizebox{0.38\textwidth}{!}{\includegraphics{inject_geographic}} \\
\parbox{0.38\textwidth}{\caption{\label{fig:geographic} Definition of
the geographic (Earth-fixed) coordinate system.}}
\end{center}
\end{wrapfigure}
\paragraph{Geographic coordinates:} This is a planetwide Earth-fixed
coordinate system, shown in Fig.~\ref{fig:geographic}.  The $z$-axis
is defined to be parallel to the Earth's axis, in the direction of the
Earth's north pole.  The $x$-axis is defined to be parallel to the
direction perpendicular from the Earth's rotation axis to a reference
point in Greenwich, UK (the \emph{prime meridian}.  Note that we adopt
a longitude convention that is consistent with the
\textit{Astronomical Almanac}, but opposite to that
in~\cite{Lang_K:1999}, in that our geographic longitudes increase
\emph{eastward} (counterclockwise) like the rest of our longitudinal
coordinates.

The terms ``latitude'' and ``longitude'' without qualification
normally refer to geographic latitude and longitude.  However, we
emphasize once again that geographic latitude and longitude as defined
above refer to directions in space, not to locations on the Earth's
surface.  This can lead to some confusion.  The \emph{geodetic}
latitude and longitude of a point on the Earth's surface are the
latitude and longitude of its vertical direction; this is the standard
meaning used by cartographers, and relates directly to the
horizon-based coordinate system above.  However, one can also define a
\emph{geocentric} latitude and longitude for a point on the surface,
which are the latitude and longitude of the direction from the
geometric centre of the Earth through that point.  These angles are
not necessarily the same, due to the Earth's ellipticity, as shown in
Fig.~\ref{fig:geodetic} in \verb@TerrestrialCoordinates.h@.

Geographic coordinates are related to sky-fixed equatorial coordinates
by specifying the counterclockwise angle \emph{to} the prime meridian
\emph{from} the reference meridian $\Upsilon$ of the sky-fixed
coordinates, as defined below.  This angle is called the Greenwich
Mean Sidereal Time (GMST), and is often specified in hours, minutes,
and seconds.

\newpage
\begin{wrapfigure}{r}{0.46\textwidth}
\vspace{-2ex}
\begin{center}
\resizebox{0.41\textwidth}{!}{\includegraphics{inject_ecliptic}} \\
\parbox{0.41\textwidth}{\caption{\label{fig:ecliptic} Definition of
the ecliptic sky-fixed coordinate systems.}}
\end{center}
\vspace{-2ex}
\end{wrapfigure}
\paragraph{Equatorial coordinates:} This is the standard sky-fixed
coordinate system.  The $z$-axis is defined as for geographic
coordinates, above; the plane orthogonal to this passing through the
Earth's centre is called the \emph{equator}.  The $x$-axis is defined
to be the direction, as viewed from the centre of the Earth, where the
Sun appears to cross the equator moving north in spring.  This is
called the \emph{vernal equinox} $\Upsilon$, and is shown in
Fig.~\ref{fig:ecliptic}.  In this coordinate system, the latitude
coordinate is called the \emph{declination} $\delta$ and the longitude
coordinate is called the \emph{right ascension} $\alpha$.

\paragraph{Ecliptic coordinates:} This is another sky-fixed coordinate
system, shown in Fig.~\ref{fig:ecliptic}.  The $z$-axis is defined to
be the direction orthogonal to the orbital plane of the Earth about
the Sun, directed such that the Earth orbits in a right-handed sense.
The $x$-axis is defined as for equatorial coordinates, above; we note
that by definition it lies parallel to the intersection of the
equatorial and orbital planes of the Earth.

The equatorial and ecliptic coordinate systems are related by a single
angle $\epsilon$, called the \emph{obliquity of the ecliptic} (that
is, the inclination of the Earth's rotation axis relative to its
orbital axis).  Ecliptic latitude is normally denoted as $\beta$ and
ecliptic longitude as $\lambda$.

\begin{wrapfigure}{l}{0.505\textwidth}
\vspace{-4ex}
\begin{center}
\resizebox{0.455\textwidth}{!}{\includegraphics{inject_galactic}} \\
\parbox{0.455\textwidth}{\caption{\label{fig:galactic} Definition of
the Galactic coordinate system.}}
\end{center}
\end{wrapfigure}
\paragraph{Galactic coordinates:} This coordinate system is shown in
Fig.~\ref{fig:galactic}.  The $z$-axis is defined to be the direction
orthogonal to the plane of our Galaxy and pointing into the northern
hemisphere of the equatorial coordinate system.  (Unfortunately this
convention has the unintuitive result that the physical rotation of
the Galaxy is left-handed about this axis.)  The $x$-axis is defined
to be the direction of the Galactic centre as viewed from the Earth.
The Galactic latitude coordinate is normally denoted as $b$ and the
Galactic longitude as $l$.

The definition of the Galactic coordinate system is completely
unrelated to any of the the other coordinate systems; thus, the
relationship between Galactic and equatorial coodinates requires one
to specify three arbitrary (but constant) angles.  Two of these are
the right ascension $\alpha_\mathrm{NGP}$ and declination
$\delta_\mathrm{NPG}$ of the North Galactic Pole ($z$-axis) in the
equatorial coordinate system, the third is the Galactic longitude
$l_\mathrm{ascend}$ of the point where the Galactic plane ascends
through the equatorial plane; i.e.\ the $l$ value for the direction
along the intersection of the Galactic and equatorial planes, such
that right-handed rotation about the Galactic $z$-axis moves you from
south to north through the equator.

******************************************************* </lalLaTeX> */

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
/*
#define SKYCOORDINATESH_ENUL  1
#define SKYCOORDINATESH_ESYS  2
#define SKYCOORDINATESH_EZERO 3
#define SKYCOORDINATESH_ESING 4

#define SKYCOORDINATESH_MSGENUL  "Unexpected null pointer in arguments"
#define SKYCOORDINATESH_MSGESYS  "Wrong coordinate system in input"
#define SKYCOORDINATESH_MSGEZERO "Angular coordinates undefined at origin"
#define SKYCOORDINATESH_MSGESING "Point is inside singular ellipsoid"
*/
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Enumeration \texttt{CoordinateSystem}}
\idx[Type]{CoordinateSystem}

This enumerated type is used to identify data as being in one of the
coordinate systems discussed above.  The allowed values are:

\idx[Constant]{COORDINATESYSTEM\_HORIZON}
\idx[Constant]{COORDINATESYSTEM\_GEOGRAPHIC}
\idx[Constant]{COORDINATESYSTEM\_EQUATORIAL}
\idx[Constant]{COORDINATESYSTEM\_ECLIPTIC}
\idx[Constant]{COORDINATESYSTEM\_GALACTIC}
\medskip\noindent
\begin{tabular}{ll}
\verb@COORDINATESYSTEM_HORIZON@ & A horizon coordinate system. \\
\verb@COORDINATESYSTEM_GEOGRAPHIC@ & The Earth-fixed geographic
coordinate system. \\
\verb@COORDINATESYSTEM_EQUATORIAL@ & The sky-fixed equatorial
coordinate system. \\
\verb@COORDINATESYSTEM_ECLIPTIC@ & The ecliptic coordinate system. \\
\verb@COORDINATESYSTEM_GALACTIC@ & The galactic coordinate system.
\end{tabular}
\bigskip

******************************************************* </lalLaTeX> */

/********************************************************** <lalLaTeX>

\subsubsection*{Structure \texttt{SkyPosition}}
\idx[Type]{SkyPosition}

This structure stores the two spherical coordinates of a sky position;
i.e.\ a generic latitude and longitude.  The structure is not defined
specific to a particular coordinate system, but maintains a tag
indicating which coordinate system it is expressed in.  The fields
are:

\begin{description}
\item[\texttt{REAL8 longitude}] The longitudinal coordinate (in
radians), as defined above.

\item[\texttt{REAL8 latitude}] The latitudinal coordinate (in
radians), as defined above.

\item[\texttt{CoordinateSystem system}] The coordinate system in which
the latitude and longitude have been expressed.
\end{description}

******************************************************* </lalLaTeX> */

/********************************************************** <lalLaTeX>

\subsubsection*{Structure \texttt{EarthPosition}}
\idx[Type]{EarthPosition}

This structure stores the location of a point on (or near) the surface
of the Earth in both geodetic and geocentric coordinates, as described
in \verb@TerrestrialCoordinates.c@.  The fields are:

\begin{description}
\item[\texttt{SkyPosition geodetic}] The geographic coordinates of the
upward vertical direction from the point; that is, the point's
\emph{geodetic} latitude and longitude.

\item[\texttt{REAL8 elevation}] The vertical distance of the point above
the reference ellipsoid, in metres.

\item[\texttt{REAL8 x, y, z}] The Earth-fixed geocentric Cartesian
coordinates of the point, in metres.

\item[\texttt{REAL8 radius}] The distance of the point from the
geocentre, in metres.

\item[\texttt{SkyPosition geocentric}] The geographic coordinates of
the direction from the centre of the Earth through the point; that is,
the point's \emph{geocentric} latitude and longitude.
\end{description}

******************************************************* </lalLaTeX> */

/********************************************************** <lalLaTeX>

\subsubsection*{Structure \texttt{ConvertSkyParams}}
\idx[Type]{ConvertSkyParams}

This structure stores parameters for the function
\verb@LALConvertSkyPosition()@.  The fields are:

\begin{description}
\item[\texttt{CoordinateSystem system}] The coordinate system to which
one is transforming.

\item[\texttt{SkyPosition *zenith}] The position of the zenith of the
horizon coordinate system; may be \verb@NULL@ if one is neither
converting to nor from a horizon system.

\item[\texttt{LIGOTimeGPS *gpsTime}] The GPS time for conversions
between Earth-fixed and sky-fixed coordinates; may be \verb@NULL@ if
no such conversion is required (or if one is transforming to or from
horizon coordinates and \verb@*zenith@ is given in the sky-fixed
equatorial system).
\end{description}

******************************************************* </lalLaTeX> */

/* <lalLaTeX>
\vfill{\footnotesize\input{SkyCoordinatesLaldocHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{CelestialCoordinatesLaldocC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{TerrestrialCoordinatesC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{SkyCoordinatesC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{SkyCoordinatesTestC}
</lalLaTeX> */

/* <lalLaTeX>
\newpage\input{GeocentricGeodeticTestC}
</lalLaTeX> */

#include <lal/LALRCSID.h>
NRCSID (SKYCOORDINATESLALDOCH,"$Id$");
