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

These functions allow one to convert celestial coordinates in the
equatorial coordinate system, to geographic coordinates in an
Earth-fixed coordinate system, to local horizon coordinates, and back
again.  The definitions of the coordinate systems are given in
\verb@SkyCoordinates.h@.  For simplicity all transformations are done
in place; \verb@*position@ stores the coordinates in the first
coordinate system at the start of the call, and the second coordinate
system upon return.  The parameter \verb@*gpsTime@ stores the time of
observation, and \verb@*zenith@ stores the geographic coordinates of
the observer's zenith (i.e.\ the observer's \emph{geodetic} longitude
and latitude).

\subsubsection*{Algorithm}

These routines follow the formulae on p.~15 of~\cite{Lang_K:1998},
which we reproduce below.  All positions are assumed to be in the
J2000 epoch.

\paragraph{Geographic coordinates:} Since geographic and equatorial
coordinates share the same $z$-axis, the geographic latitude $\phi$ of
a direction in space is the same as its declination $\delta$, and
longitude $\lambda$ and right ascension $\alpha$ differ only through
the rotation of the Earth:
\begin{equation}
\label{eq:lambda-geographic}
\lambda = \alpha + \left(\frac{2\pi\,\mathrm{radians}}
	{24\,\mathrm{hours}}\right)\times\mathrm{GMST} \; ,
\end{equation}
where GMST is Greenwich mean sidereal time.  The conversion routines
here simply use the functions in the \verb@date@ package to compute
GMST for a given GPS time, and add it to the longitude.  While this is
simple enough, it does involve several function calls, so it is
convenient to collect these into one routine.

\paragraph{Horizon coordinates:} We correct a typographical
error in~\cite{Lang_K:1998} (the second line of Eq.~5.45 should have a
$\cos A$), to obtain the following equations for the altitude $a$ and
azimuth $A$:
\begin{eqnarray}
\label{eq:altitude-horizon}
a & = & \arcsin(\sin\delta\sin\phi + \cos\delta\cos\phi\cos h) \; , \\
\label{eq:azimuth-horizon}
A & = & \arctan\!2(-\cos\delta\sin h, \sin\delta\cos\phi -
		\cos\delta\sin\phi\cos h) \; ,
\end{eqnarray}
where $\delta$ is the declination (geographic latitude) of the
position being transformed, $\phi$ is the geographic latitude of the
observer's zenith (i.e.\ the observer's \emph{geodetic} latitude), and
$h$ is the \emph{hour angle} of the position being transformed.  This
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
h & = & \arctan\!2(-\cos a\sin A, \sin a\cos\phi -
		\cos a\cos A\sin\phi) \; .
\end{eqnarray}
As explained in \verb@CelestialCoordinates.c@, the function
$\arctan\!2(y,x)$ returns the argument of the complex number $x+iy$.

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TerrestrialCoordinatesCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>

NRCSID( TERRESTRIALCOORDINATESC, "$Id$" );

/* <lalVerbatim file="TerrestrialCoordinatesCP"> */
void
LALEquatorialToGeographic( LALStatus   *stat,
			   SkyPosition *position,
			   LIGOTimeGPS *gpsTime )
{ /* </lalVerbatim> */
  LIGOTimeUnix unixTime; /* Unix time */
  LALDate date;          /* LALDate time */
  REAL8 gmst;            /* siderial time (radians) */

  INITSTATUS( stat, "LALEquatorialToGeographic", TERRESTRIALCOORDINATESC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures exist. */
  ASSERT( position, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( gpsTime, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Compute the Greenwich mean sidereal time. */
  TRY( LALGPStoU( stat->statusPtr, &unixTime, gpsTime ), stat );
  TRY( LALUtime( stat->statusPtr, &date, &unixTime ), stat );
  TRY( LALGMST1( stat->statusPtr, &gmst, &date, MST_RAD ), stat );

  /* Add to longitude, and exit. */
  position->longitude += gmst;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="TerrestrialCoordinatesCP"> */
void
LALGeographicToEquatorial( LALStatus   *stat,
			   SkyPosition *position,
			   LIGOTimeGPS *gpsTime )
{ /* </lalVerbatim> */
  LIGOTimeUnix unixTime; /* Unix time */
  LALDate date;          /* LALDate time */
  REAL8 gmst;            /* siderial time (radians) */

  INITSTATUS( stat, "LALEquatorialToGeographic", TERRESTRIALCOORDINATESC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structures exist. */
  ASSERT( position, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( gpsTime, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Compute the Greenwich mean sidereal time. */
  TRY( LALGPStoU( stat->statusPtr, &unixTime, gpsTime ), stat );
  TRY( LALUtime( stat->statusPtr, &date, &unixTime ), stat );
  TRY( LALGMST1( stat->statusPtr, &gmst, &date, MST_RAD ), stat );

  /* Subtract from longitude, and exit. */
  position->longitude -= gmst;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="TerrestrialCoordinatesCP"> */
void
LALGeographicToHorizon( LALStatus   *stat,
			SkyPosition *position,
			SkyPosition *zenith )
{ /* </lalVerbatim> */
  REAL8 h, sinH, cosH; /* hour angle, and its sine and cosine */
  REAL8 sinP, cosP;    /* sin and cos of zenith latitude */
  REAL8 sinD, cosD;    /* sin and cos of position latitude (declination) */
  REAL8 sina, sinA, cosA; /* sin and cos of altitude and azimuth */

  INITSTATUS( stat, "LALGeographicToHorizon", TERRESTRIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( position, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( zenith, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Compute intermediates. */
  h = zenith->longitude - position->longitude;
  sinH = sin( h );
  cosH = cos( h );
  sinP = sin( zenith->latitude );
  cosP = cos( zenith->latitude );
  sinD = sin( position->latitude );
  cosD = cos( position->latitude );

  /* Compute components. */
  sina = sinD*sinP + cosD*cosP*cosH;
  sinA = -cosD*sinH;
  cosA = sinD*cosP - cosD*sinP*cosH;

  /* Compute final results. */
  position->latitude = asin( sina );
  position->longitude = atan2( sinA, cosA );

  /* Optional phase correction. */
  if ( position->longitude < 0.0 )
    position->longitude += LAL_TWOPI;

  RETURN( stat );
}


/* <lalVerbatim file="TerrestrialCoordinatesCP"> */
void
LALHorizonToGeographic( LALStatus   *stat,
			SkyPosition *position,
			SkyPosition *zenith )
{ /* </lalVerbatim> */
  REAL8 sinP, cosP;       /* sin and cos of zenith latitude */
  REAL8 sina, cosa;       /* sin and cos of altitude */
  REAL8 sinA, cosA;       /* sin and cos of altitude */
  REAL8 sinD, sinH, cosH; /* sin and cos of declination and hour angle */

  INITSTATUS( stat, "LALHorizonToGeographic", TERRESTRIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( position, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( zenith, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Compute intermediates. */
  sinP = sin( zenith->latitude );
  cosP = cos( zenith->latitude );
  sina = sin( position->latitude );
  cosa = cos( position->latitude );
  sinA = sin( position->longitude );
  cosA = cos( position->longitude );

  /* Compute components. */
  sinD = sina*sinP + cosa*cosA*cosP;
  sinH = -cosa*sinA;
  cosH = sina*cosP - cosa*cosA*sinP;

  /* Compute final results. */
  position->latitude = asin( sinD );
  position->longitude = zenith->longitude - atan2( sinH, cosH );

  /* Optional phase correction. */
  if ( position->longitude < 0.0 )
    position->longitude += LAL_TWOPI;

  RETURN( stat );
}
