/*********************** <lalVerbatim file="TerrestrialCoordinatesCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{TerrestrialCoordinates.c}}
\label{ss:TerrestrialCoordinates.c}

Converts among equatorial, geographic, and horizon coordinates.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{TerrestrialCoordinatesCP}
\index{\texttt{LALEquatorialToGeographic()}}
\index{\texttt{LALGeographicToEquatorial()}}
\index{\texttt{LALGeographicToHorizon()}}
\index{\texttt{LALHorizonToGeographic()}}

\subsubsection*{Description}

The functions \verb@LALEquatorialToGeographic()@ and
\verb@LALGeographicToEquatorial()@ convert between equatorial and
geographic coordinate systems, reading coordinates in the first system
from \verb@*input@ and storing the new coordinates in \verb@*output@.
The two pointers may point to the same object, in which case the
conversion is done in place.  The functions will also check
\verb@input->system@ and set \verb@output->system@ as appropriate.
Because the geographic coordinate system is not fixed, one must also
specify the time of the transformation in \verb@*gpsTime@.

The function \verb@LALSystemToHorizon()@ transforms coordinates from
an arbitrary coordinate to a horizon coordinate system, reading
coordinates in the first system from \verb@*input@ and storing the
horizon coordinates in \verb@*output@, as above.  The parameter
\verb@*zenith@ specifies the direction of the vertical axis \emph{in
the original coordinate system}; the routine checks to see that
\verb@input->system@ and \verb@zenith->system@ agree.  Normally this
routine is used to convert from \emph{geographic} latitude and
longitude to a horizon system, in which case \verb@*zenith@ simply
stores the geographic (geodetic) coordinates of the observer.

The function \verb@LALHorizonToSystem()@ does the reverse of the
above, transforming coordinates from horizon coordinates to an
arbitrary coordinate system specified by \verb@zenith->system@; the
value of \verb@output->system@ is set to agree with
\verb@zenith->system@.

Although it is conventional to specify an observation location by its
\emph{geodetic} coordinates, some routines may provide or require
\emph{geocentric} coordinates.  The routines
\verb@LALGeocentricToGeodetic()@ and \verb@LALGeodeticToGeocentric()@
perform this computation, reading and writing to the variable
parameter structure \verb@*location@.  The function
\verb@LALGCentricToGDetic()@ reads the fields \verb@location->x@,
\verb@y@, \verb@z@, and computes \verb@location->zenith@ and
\verb@location->altitude@.  The function \verb@LALGDeticToGCentric()@
does the reverse, and also sets the fields \verb@location->position@
and \verb@location->radius@.

\subsubsection*{Algorithm}

These routines follow the formulae in Sec.~5.1 of~\cite{Lang_K:1998},
which we reproduce below.

\paragraph{Geographic coordinates:} Since geographic and equatorial
coordinates share the same $z$-axis, the geographic latitude $\phi$ of
a direction in space is the same as its declination $\delta$, and
longitude $\lambda$ and right ascension $\alpha$ differ only through
the rotation of the Earth:
\begin{equation}
\label{eq:lambda-geographic}
\lambda = \alpha - \left(\frac{2\pi\,\mathrm{radians}}
	{24\,\mathrm{hours}}\right)\times\mathrm{GMST} \; ,
\end{equation}
where GMST is Greenwich mean sidereal time.  The conversion routines
here simply use the functions in the \verb@date@ package to compute
GMST for a given GPS time, and add it to the longitude.  While this is
simple enough, it does involve several function calls, so it is
convenient to collect these into one routine.

\paragraph{Horizon coordinates:} We correct a typographical
error on the second line of Eq.~5.45 of~\cite{Lang_K:1998} (it should
have $\cos A$, not $\sin A$).  We also note that while our latitudinal
coordinate is just the altitude $a$ in this system, our longitudinal
coordinate increases counterclockwise, and thus corresponds to the
\emph{negative} of the azimuth $A$ as defined by~\cite{Lang_K:1998}.
So we have:
\begin{eqnarray}
\label{eq:altitude-horizon}
a & = & \arcsin(\sin\delta\sin\phi + \cos\delta\cos\phi\cos h) \; , \\
\label{eq:azimuth-horizon}
-A & = & \arctan\!2(\cos\delta\sin h, \sin\delta\cos\phi -
		\cos\delta\sin\phi\cos h) \; ,
\end{eqnarray}
where $\delta$ is the declination (geographic latitude) of the
direction being transformed, $\phi$ is the geographic latitude of the
observer's zenith (i.e.\ the observer's \emph{geodetic} latitude), and
$h$ is the \emph{hour angle} of the direction being transformed.  This
is defined as:
\begin{eqnarray}
h & = & \lambda_\mathrm{zenith} - \lambda \nonumber\\
  & = & \mathrm{LMST} - \alpha \nonumber
\end{eqnarray}
where LMST is the local mean sidereal time at the point of
observation.  The inverse transformation is:
\begin{eqnarray}
\label{eq:delta-horizon}
\delta & = & \arcsin(\sin a\sin\phi + \cos a\cos A\cos\phi) \; , \\
\label{eq:h-horizon}
h & = & \arctan\!2[\cos a\sin(-A), \sin a\cos\phi -
		\cos a\cos A\sin\phi] \; .
\end{eqnarray}
As explained in \verb@CelestialCoordinates.c@, the function
$\arctan\!2(y,x)$ returns the argument of the complex number $x+iy$.

\newpage

\begin{wrapfigure}{r}{0.35\textwidth}
\vspace{-2ex}
\begin{center}
\resizebox{0.3\textwidth}{!}{\includegraphics{inject_geodetic}} \\
\parbox{0.3\textwidth}{\caption{\label{fig:geodetic} The difference
between geodetic and geocentric latitude.}}
\end{center}
\vspace{-2ex}
\end{wrapfigure}
\paragraph{Geocentric coordinates:} As shown in
Fig.~\ref{fig:geodetic}, the ellipticity of the Earth means that the
vertical axis of a point on the Earth's surface does not pass through
the geometric centre of the Earth.  This means that the geodetic
latitude of a location (defined as the latitude angle
$\phi_\mathrm{geodetic}$ of that location's zenith direction) is
typically some 10~arcminutes larger than its geocentric latitude
(defined as the latitude angle $\phi_\mathrm{geographic}$ of the
position vector from the geocentre through the location).
Cartographers traditionally refer to locations by their geodetic
coordinates, since these can be determined locally; however,
geocentric coordinates are required if one wants to construct a
uniform Cartesian system for the Earth as a whole.

To transform from geodetic to geocentric coordinates, one first
defines a ``reference ellipsoid'', the best-fit ellipsoid to the
surface of the Earth.  This is specified by the polar and equatorial
radii of the Earth $r_p$ and $r_e$, or equivalently by $r_e$ and a
flattening factor:
$$
f \equiv 1 - \frac{r_p}{r_e} = 0.00335281 \; .
$$
(This constant will eventually migrate into \verb@LALConstants.h@.)
The surface of the ellipsoid is then specified by the equation
$$
r = r_e ( 1 - f\sin^2\phi ) \; ,
$$
where $\phi=\phi_\mathrm{geodetic}$ is the geodetic latitude.  For
points off of the reference ellipsoid, the transformation from
geodetic coordinates $\lambda$, $\phi$ to geocentric Cartesian
coordinates is:
\begin{eqnarray}
x & = & ( r_e C + h ) \cos\phi\cos\lambda \; , \\
y & = & ( r_e C + h ) \cos\phi\sin\lambda \; , \\
z & = & ( r_e S + h ) \sin\phi \; ,
\end{eqnarray}
where
\begin{eqnarray}
C & = & \frac{1}{\sqrt{\cos^2\phi + (1-f)^2\sin^2\phi}} \; , \\
S & = & (1-f)^2 C \; ,
\end{eqnarray}
and $h$ is the perpendicular height of the location above the
reference ellipsoid.  The geocentric spherical coordinates are given
simply by:
\begin{eqnarray}
r & = & \sqrt{ x^2 + y^2 + z^2 } \; , \\
\lambda_\mathrm{geocentric} & = & \lambda \quad = \quad
	\lambda_\mathrm{geodetic} \; , \\
\phi_\mathrm{geocentric} & = & \arcsin( z/r ) \; .
\end{eqnarray}
When computing $r$ we are careful to factor out the largest component
before computing the sum of squares, to avoid floating-point overflow;
however this should be unnecessary for radii near the surface of the
Earth.

The inverse transformation is somewhat trickier, and I'm still working
on it.  I will commit it when I'm done.

%Given geocentric Cartesian coordinates
%$x=r\cos\phi_\mathrm{geodetic}\cos\lambda$,
%$y=r\cos\phi_\mathrm{geodetic}\sin\lambda$, and
%$z=r\sin\phi_\mathrm{geodetic}$, one computes a number of
%intermediates:
%\begin{eqnarray}
%\varpi & = & \sqrt{ x^2 + y^2 } \;,\nonumber\\
%E & = & (1-f)\frac{z}{r} - f(2-f)\frac{r_e}{r} \;,\nonumber\\
%F & = & (1-f)\frac{z}{r} + f(2-f)\frac{r_e}{r} \;,\nonumber\\
%P & = & \frac{4}{3}( EF + 1 ) \;,\nonumber\\
%Q & = & 2

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TerrestrialCoordinatesCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALError.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>

#define LAL_EARTHFLAT (0.00335281)

NRCSID( TERRESTRIALCOORDINATESC, "$Id$" );

/* <lalVerbatim file="TerrestrialCoordinatesCP"> */
void
LALEquatorialToGeographic( LALStatus   *stat,
			   SkyPosition *output,
			   SkyPosition *input,
			   LIGOTimeGPS *gpsTime )
{ /* </lalVerbatim> */
  LIGOTimeUnix unixTime; /* Unix time */
  LALDate date;          /* LALDate time */
  REAL8 gmst;            /* siderial time (radians) */

  INITSTATUS( stat, "LALEquatorialToGeographic", TERRESTRIALCOORDINATESC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures exist. */
  ASSERT( input, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( output, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( gpsTime, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Make sure we're given the right coordinate system. */
  if ( input->system != COORDINATESYSTEM_EQUATORIAL ) {
    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  /* Compute the Greenwich mean sidereal time. */
  TRY( LALGPStoU( stat->statusPtr, &unixTime, gpsTime ), stat );
  TRY( LALUtime( stat->statusPtr, &date, &unixTime ), stat );
  TRY( LALGMST1( stat->statusPtr, &gmst, &date, MST_RAD ), stat );

  /* Add to longitude, and exit. */
  output->system = COORDINATESYSTEM_GEOGRAPHIC;
  output->longitude = input->longitude - gmst;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="TerrestrialCoordinatesCP"> */
void
LALGeographicToEquatorial( LALStatus   *stat,
			   SkyPosition *output,
			   SkyPosition *input,
			   LIGOTimeGPS *gpsTime )
{ /* </lalVerbatim> */
  LIGOTimeUnix unixTime; /* Unix time */
  LALDate date;          /* LALDate time */
  REAL8 gmst;            /* siderial time (radians) */

  INITSTATUS( stat, "LALEquatorialToGeographic", TERRESTRIALCOORDINATESC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures exist. */
  ASSERT( input, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( output, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( gpsTime, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Make sure we're given the right coordinate system. */
  if ( input->system != COORDINATESYSTEM_GEOGRAPHIC ) {
    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  /* Compute the Greenwich mean sidereal time. */
  TRY( LALGPStoU( stat->statusPtr, &unixTime, gpsTime ), stat );
  TRY( LALUtime( stat->statusPtr, &date, &unixTime ), stat );
  TRY( LALGMST1( stat->statusPtr, &gmst, &date, MST_RAD ), stat );

  /* Subtract from longitude, and exit. */
  output->system = COORDINATESYSTEM_EQUATORIAL;
  output->longitude = input->longitude + gmst;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="TerrestrialCoordinatesCP"> */
void
LALSystemToHorizon( LALStatus   *stat,
		    SkyPosition *output,
		    SkyPosition *input,
		    const SkyPosition *zenith )
{ /* </lalVerbatim> */
  REAL8 h, sinH, cosH; /* hour angle, and its sine and cosine */
  REAL8 sinP, cosP;    /* sin and cos of zenith latitude */
  REAL8 sinD, cosD;    /* sin and cos of position latitude (declination) */
  REAL8 sina, sinA, cosA; /* sin and cos of altitude and -azimuth */

  INITSTATUS( stat, "LALSystemToHorizon", TERRESTRIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( input, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( output, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( zenith, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Make sure we have consistent input coordinate systems. */
  if ( input->system != zenith->system ) {
    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  /* Compute intermediates. */
  h = zenith->longitude - input->longitude;
  sinH = sin( h );
  cosH = cos( h );
  sinP = sin( zenith->latitude );
  cosP = cos( zenith->latitude );
  sinD = sin( input->latitude );
  cosD = cos( input->latitude );

  /* Compute components. */
  sina = sinD*sinP + cosD*cosP*cosH;
  sinA = cosD*sinH;
  cosA = sinD*cosP - cosD*sinP*cosH;

  /* Compute final results. */
  output->system = COORDINATESYSTEM_HORIZON;
  output->latitude = asin( sina );
  output->longitude = atan2( sinA, cosA );

  /* Optional phase correction. */
  if ( output->longitude < 0.0 )
    output->longitude += LAL_TWOPI;

  RETURN( stat );
}


/* <lalVerbatim file="TerrestrialCoordinatesCP"> */
void
LALHorizonToSystem( LALStatus   *stat,
		    SkyPosition *output,
		    SkyPosition *input,
		    const SkyPosition *zenith )
{ /* </lalVerbatim> */
  REAL8 sinP, cosP;       /* sin and cos of zenith latitude */
  REAL8 sina, cosa;       /* sin and cos of altitude */
  REAL8 sinA, cosA;       /* sin and cos of -azimuth */
  REAL8 sinD, sinH, cosH; /* sin and cos of declination and hour angle */

  INITSTATUS( stat, "LALHorizonToSystem", TERRESTRIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( input, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( output, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( zenith, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Make sure we're given the right coordinate system. */
  if ( input->system != COORDINATESYSTEM_HORIZON ) {
    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  /* Compute intermediates. */
  sinP = sin( zenith->latitude );
  cosP = cos( zenith->latitude );
  sina = sin( input->latitude );
  cosa = cos( input->latitude );
  sinA = sin( input->longitude );
  cosA = cos( input->longitude );

  /* Compute components. */
  sinD = sina*sinP + cosa*cosA*cosP;
  sinH = cosa*sinA;
  cosH = sina*cosP - cosa*cosA*sinP;

  /* Compute final results. */
  output->system = zenith->system;
  output->latitude = asin( sinD );
  output->longitude = zenith->longitude - atan2( sinH, cosH );

  /* Optional phase correction. */
  if ( output->longitude < 0.0 )
    output->longitude += LAL_TWOPI;

  RETURN( stat );
}


/* <lalVerbatim file="TerrestrialCoordinatesCP"> */
void
LALGeodeticToGeocentric( LALStatus *stat, EarthPosition *location )
{ /* </lalVerbatim> */
  REAL8 c, s; /* position components in and orthogonal to the equator */
  REAL8 cosP, sinP; /* cosine and sine of latitude */
  REAL8 fFac;       /* ( 1 - f )^2 */
  REAL8 x, y;       /* Cartesian coordinates */
  REAL8 maxComp, r; /* max{x,y,z}, and sqrt(x^2+y^2+z^2) */

  INITSTATUS( stat, "LALGeodeticToGeocentric", TERRESTRIALCOORDINATESC );

  /* Make sure parameter structure exists. */
  ASSERT( location, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Make sure the geodetic coordinates are in the right system. */
  if ( location->geodetic.system != COORDINATESYSTEM_GEOGRAPHIC ) {
    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  /* Compute intermediates. */
  fFac = 1.0 - LAL_EARTHFLAT;
  fFac *= fFac;
  cosP = cos( location->geodetic.latitude );
  sinP = sin( location->geodetic.latitude );
  c = sqrt( 1.0 / ( cosP*cosP + fFac*sinP*sinP ) );
  s = fFac*c;
  c = ( LAL_REARTH_SI*c + location->height )*cosP;
  s = ( LAL_REARTH_SI*s + location->height )*sinP;

  /* Compute Cartesian coordinates. */
  location->x = x = c*cos( location->geodetic.longitude );
  location->y = y = c*sin( location->geodetic.longitude );
  location->z = s;

  /* Compute the radius. */
  maxComp = x;
  if ( y > maxComp )
    maxComp = y;
  if ( s > maxComp )
    maxComp = s;
  x /= maxComp;
  y /= maxComp;
  s /= maxComp;
  r = sqrt( x*x + y*y + s*s );

  /* Compute the spherical coordinates, and exit. */
  location->radius = maxComp*r;
  location->geocentric.longitude = location->geodetic.longitude;
  location->geocentric.latitude = asin( s / r );
  location->geocentric.system = COORDINATESYSTEM_GEOGRAPHIC;
  RETURN( stat );
}


/* <lalVerbatim file="TerrestrialCoordinatesCP"> */
void
LALGeocentricToGeodetic( LALStatus *stat, EarthPosition *location )
{ /* </lalVerbatim> */

  INITSTATUS( stat, "LALGeocentricToGeodetic", TERRESTRIALCOORDINATESC );

  /* Make sure parameter structure exists. */
  ASSERT( location, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Make sure the geodetic coordinates are in the right system. */
  if ( location->geocentric.system != COORDINATESYSTEM_GEOGRAPHIC ) {
    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  LALWarning( stat, "This function is just a stub." );

  RETURN( stat );
}
