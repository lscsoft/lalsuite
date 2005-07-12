/** \file
 *  \ingroup SkyCoordinates
 *  \author Creighton, T. D.
 *  \date $Date$
 *  \brief Converts among Galactic, ecliptic, and equatorial coordinates.
 * 
 * $Id$
 * 
These functions perform the specified coordinate transformation on the
contents of \a input and store the result in \a *output.  The
two pointers may point to the same object, in which case the
conversion is done in place.  The functions will also check
<tt>input->system</tt> and set <tt>output->system</tt> as appropriate.

These routines are collected together because they involve fixed,
absolute coordinate systems, so the transformations require no
additional parameters such as the time or site of observation.  We
also note that there are no direct conversions between Galactic and
ecliptic coordinates.  At the risk of additional computational
overhead, it is simple to use the equatorial coordinate system as an
intermediate step.

\par Description

These functions perform the specified coordinate transformation on the
contents of \a input and store the result in \a output.  The
two pointers may point to the same object, in which case the
conversion is done in place.  The functions will also check
<tt>input->system</tt> and set <tt>output->system</tt> as appropriate.

These routines are collected together because they involve fixed,
absolute coordinate systems, so the transformations require no
additional parameters such as the time or site of observation.  We
also note that there are no direct conversions between Galactic and
ecliptic coordinates.  At the risk of additional computational
overhead, it is simple to use the equatorial coordinate system as an
intermediate step.

\par Algorithm

These routines follow the spherical angle relations on p. 13
of \ref Lang_K1999.  Note that the actual formulae for Galactic
longitude and right ascension in this reference are wrong; we give
corrected formulae below derived from the sine and cosine equations.
(The Galactic to equatorial transformations can also be found in
Sec. 12.3 of \ref GRASP2000.  All positions are assumed to
be in the J2000 epoch.


<b>Galactic coordinates:</b> The following formulae relate
Galactic latitude \f$b\f$ and longitude \f$l\f$ to declination \f$\delta\f$ and
right ascension \f$\alpha\f$:
\f{eqnarray}
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
\f}
where \f$\arctan\!2(y,x)\f$ can be thought of as the argument of the
complex number \f$x+iy\f$; unlike \f$\arctan(y/x)\f$, it ranges over the full
range \f$[0,2\pi)\f$ instead of just half of it.  The inverse
transformations are:
\f{eqnarray}
\label{eq:delta-galactic}
\delta & = & \arcsin[\cos b\cos\delta_\mathrm{NGP}\sin(l-l_\mathrm{ascend}) +
		\sin b\sin\delta_\mathrm{NGP}] \;,\\
\alpha & = & \arctan\!2[\cos b\cos(l-l_\mathrm{ascend}),
		\sin b\cos\delta_\mathrm{NGP} -
		\cos b\sin(l-l_\mathrm{ascend})\sin\delta_\mathrm{NGP}]
		\nonumber\\
\label{eq:alpha-galactic}
& & \quad + \; \alpha_\mathrm{NGP} \;.
\f}
In these equations we have defined the orientation of the Galaxy with
the following parameters (which should eventually be placed in <tt>LALConstants.h</tt>:
\f[ 
\begin{array}{r@{\quad=\quad}l@{\quad=\quad}l}
\alpha_\mathrm{NGP} & 192.8594813^\circ &
\mbox{the right ascension (epoch J2000) of the north Galactic pole} \\
\delta_\mathrm{NGP} & 27.1282511^\circ &
\mbox{the declination (epoch J2000) of the north Galactic pole} \\
l_\mathrm{ascend} & 33^\circ &
\mbox{the longitude of the ascending node of the Galactic plane}
\end{array}
\f]
The ascending node of the Galactic plane is defined as the direction
along the intersection of the Galactic and equatorial planes where
rotation in the positive sense about the Galactic \f$z\f$ axis carries a
point from the southern to northern equatorial hemisphere.  That is,
if \f$\mathbf{u}\f$ points in the direction \f$\delta=90^\circ\f$
(celestial north), and \f$\mathbf{v}\f$ points in the direction
\f$b=90^\circ\f$ (Galactic north), then
\f$\mathbf{u} \times \mathbf{v}\f$ points along the ascending node.

<b>Ecliptic coordinates:</b> The following formulae relate
Ecliptic latitude \f$\beta\f$ and longitude \f$\lambda\f$ to declination
\f$\delta\f$ and right ascension \f$\alpha\f$:
\f{eqnarray}
\label{eq:beta-ecliptic}
\beta & = & \arcsin(\sin\delta\cos\epsilon -
		\cos\delta\sin\alpha\sin\epsilon) \;, \\
\label{eq:l-ecliptic}
\lambda & = & \arctan\!2(\cos\delta\sin\alpha\cos\epsilon +
		\sin\delta\sin\epsilon, \cos\delta\cos\alpha) \;.
\f}
The inverse transformations are:
\f{eqnarray}
\label{eq:delta-ecliptic}
\delta & = & \arcsin(\cos\beta\sin\lambda\sin\epsilon +
		\sin\beta\cos\epsilon) \;, \\
\label{eq:alpha-ecliptic}
\alpha & = & \arctan\!2(\cos\beta\sin\lambda\cos\epsilon -
		\sin\beta\sin\epsilon, \cos\beta\cos\lambda) \;.
\f}
Here \f$\epsilon\f$ is the obliquity (inclination) of the ecliptic plane,
which varies over time; at epoch J200 it has a mean value of:
\f[
\epsilon = 23.4392911^\circ \; .
\f]

*/

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/SkyCoordinates.h>

#define LAL_ALPHAGAL (3.366032942)
#define LAL_DELTAGAL (0.473477302)
#define LAL_LGAL     (0.576)

NRCSID( CELESTIALCOORDINATESC, "$Id$" );

void
LALGalacticToEquatorial( LALStatus   *stat,
			 SkyPosition *output,
			 SkyPosition *input )
{
  REAL8 sinDGal = sin( LAL_DELTAGAL ); /* sin(delta_NGP) */
  REAL8 cosDGal = cos( LAL_DELTAGAL ); /* cos(delta_NGP) */
  REAL8 l = -LAL_LGAL;          /* will be l-l(ascend) */
  REAL8 sinB, cosB, sinL, cosL; /* sin and cos of b and l */
  REAL8 sinD, sinA, cosA;       /* sin and cos of delta and alpha */

  INITSTATUS( stat, "LALGalacticToEquatorial", CELESTIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( input, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( output, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Make sure we're given the right coordinate system. */
  if ( input->system != COORDINATESYSTEM_GALACTIC ) {
    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  /* Compute intermediates. */
  l += input->longitude;
  sinB = sin( input->latitude );
  cosB = cos( input->latitude );
  sinL = sin( l );
  cosL = cos( l );

  /* Compute components. */
  sinD = cosB*cosDGal*sinL + sinB*sinDGal;
  sinA = cosB*cosL;
  cosA = sinB*cosDGal - cosB*sinL*sinDGal;

  /* Compute final results. */
  output->system = COORDINATESYSTEM_EQUATORIAL;
  output->latitude = asin( sinD );
  l = atan2( sinA, cosA ) + LAL_ALPHAGAL;

  /* Optional phase correction. */
  if ( l < 0.0 )
    l += LAL_TWOPI;
  output->longitude = l;

  RETURN( stat );
}


void
LALEquatorialToGalactic( LALStatus   *stat,
			 SkyPosition *output,
			 SkyPosition *input )
{
  REAL8 sinDGal = sin( LAL_DELTAGAL ); /* sin(delta_NGP) */
  REAL8 cosDGal = cos( LAL_DELTAGAL ); /* cos(delta_NGP) */
  REAL8 a = -LAL_ALPHAGAL;      /* will be alpha-alpha_NGP */
  REAL8 sinD, cosD, sinA, cosA; /* sin and cos of delta and alpha */
  REAL8 sinB, sinL, cosL;       /* sin and cos of b and l */

  INITSTATUS( stat, "LALGalacticToEquatorial", CELESTIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( input, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( output, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Make sure we're given the right coordinate system. */
  if ( input->system != COORDINATESYSTEM_EQUATORIAL ) {
    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  /* Compute intermediates. */
  a += input->longitude;
  sinD = sin( input->latitude );
  cosD = cos( input->latitude );
  sinA = sin( a );
  cosA = cos( a );

  /* Compute components. */
  sinB = cosD*cosDGal*cosA + sinD*sinDGal;
  sinL = sinD*cosDGal - cosD*cosA*sinDGal;
  cosL = cosD*sinA;

  /* Compute final results. */
  output->system = COORDINATESYSTEM_GALACTIC;
  output->latitude = asin( sinB );
  a = atan2( sinL, cosL ) + LAL_LGAL;

  /* Optional phase correction. */
  if ( a < 0.0 )
    a += LAL_TWOPI;
  output->longitude = a;

  RETURN( stat );
}


void
LALEclipticToEquatorial( LALStatus   *stat,
			 SkyPosition *output,
			 SkyPosition *input )
{
  REAL8 sinE = sin( LAL_IEARTH ); /* sin(epsilon) */
  REAL8 cosE = cos( LAL_IEARTH ); /* cos(epsilon) */
  REAL8 sinB, cosB, sinL, cosL;   /* sin and cos of b and l */
  REAL8 sinD, sinA, cosA;         /* sin and cos of delta and alpha */

  INITSTATUS( stat, "LALGalacticToEquatorial", CELESTIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( input, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( output, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Make sure we're given the right coordinate system. */
  if ( input->system != COORDINATESYSTEM_ECLIPTIC ) {
    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  /* Compute intermediates. */
  sinB = sin( input->latitude );
  cosB = cos( input->latitude );
  sinL = sin( input->longitude );
  cosL = cos( input->longitude );

  /* Compute components. */
  sinD = cosB*sinL*sinE + sinB*cosE;
  sinA = cosB*sinL*cosE - sinB*sinE;
  cosA = cosB*cosL;

  /* Compute final results. */
  output->system = COORDINATESYSTEM_EQUATORIAL;
  output->latitude = asin( sinD );
  output->longitude = atan2( sinA, cosA );

  /* Optional phase correction. */
  if ( output->longitude < 0.0 )
    output->longitude += LAL_TWOPI;

  RETURN( stat );
}


void
LALEquatorialToEcliptic( LALStatus   *stat,
			 SkyPosition *output,
			 SkyPosition *input )
{
  REAL8 sinE = sin( LAL_IEARTH ); /* sin(epsilon) */
  REAL8 cosE = cos( LAL_IEARTH ); /* cos(epsilon) */
  REAL8 sinD, cosD, sinA, cosA;   /* sin and cos of delta and alpha */
  REAL8 sinB, sinL, cosL;         /* sin and cos of b and l */

  INITSTATUS( stat, "LALGalacticToEquatorial", CELESTIALCOORDINATESC );

  /* Make sure parameter structures exist. */
  ASSERT( input, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );
  ASSERT( output, stat, SKYCOORDINATESH_ENUL, SKYCOORDINATESH_MSGENUL );

  /* Make sure we're given the right coordinate system. */
  if ( input->system != COORDINATESYSTEM_EQUATORIAL ) {
    ABORT( stat, SKYCOORDINATESH_ESYS, SKYCOORDINATESH_MSGESYS );
  }

  /* Compute intermediates. */
  sinD = sin( input->latitude );
  cosD = cos( input->latitude );
  sinA = sin( input->longitude );
  cosA = cos( input->longitude );

  /* Compute components. */
  sinB = sinD*cosE - cosD*sinA*sinE;
  sinL = cosD*sinA*cosE + sinD*sinE;
  cosL = cosD*cosA;

  /* Compute final results. */
  output->system = COORDINATESYSTEM_ECLIPTIC;
  output->latitude = asin( sinB );
  output->longitude = atan2( sinL, cosL );

  /* Optional phase correction. */
  if ( output->longitude < 0.0 )
    output->longitude += LAL_TWOPI;

  RETURN( stat );
}
