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

/**
\author Creighton, T. D.
\file
\ingroup pulsarTODO

\heading{Header \ref PulsarCat.h}
\latexonly\label{s_PulsarCat_h}\endlatexonly

Provides structures and routines to store and manipulate pulsar
properties.

\heading{Synopsis}
\code
#include <lal/PulsarCat.h>
\endcode

This header covers structures to store pulsar properties in
a standard format, and routines to manipulate and update these
properties.  The set of properties stored in the catalogue is based on
radio pulsar catalogues, with some additions and subtractions specific
to gravitational wave observations.  The list of properties can be
expanded in future by adding more fields to the structure.

All properties are those that would be measured by an observer at the
solar system barycentre.  For properties that depend on time (e.g.\
time-varying position, periods, etc.), an epoch is specified.  The
properties are then those that would be observed at the specified
instant of time at the solar system barycentre; i.e.\ when the wave
fronts carrying that information pass the solar system barycentre.

\heading{A note on companion orbits:} Several known pulsars exist in
multiple systems, and radio-pulsar catalogues include detailed models
of the companion orbits, as determined from the pulsar timing.  See
\ref GenerateSpinOrbitCW.h in the \c inject package for a
discussion of the parameters defining the orientation of a companion
orbit.

Radio-pulsar observations rarely determine the inclination \f$i\f$ of the
orbit to the sky plane, and thus cannot resolve the longitude of the
ascending node \f$\Omega\f$ and the argument of the periapsis \f$\omega\f$ as
independent parameters.  Instead, they list the longitude of the
periapsis \f$w\f$, which is the angle in the plane of the sky from the
North direction towards the West direction, to the ray from the system
barycentre to the periapsis projected onto the plane of the sky.  If
any three of \f$i\f$, \f$\Omega\f$, \f$\omega\f$, and \f$w\f$ are known, the fourth
can be determined from the relation:
\f[
w - \Omega = \arctan\;2(\sin\omega\cos i,\cos\omega) \;,
\f]
or equivalently:
\f[
\omega = \arctan\;2(\cos[w-\Omega],\sin[w-\Omega]/\cos i) \;.
\f]

In addition to these Keplerian orbital parameters, some radio-pulsar
systems have measured post-Keplerian relativistic orbital parameters.
Some of these are obvious: \f$\dot{w}\f$ and \f$\dot{P}\f$ are the rate of
change in \f$w\f$ (periapsis precession) and the orbital period \f$P\f$ (due
to gravitational radiation reaction).  The catalogue also lists
post-Keplerian parameters "sin" and "r", whose meanings I don't
know.

*/

#ifndef _PULSARCAT_H
#define _PULSARCAT_H

#include <lal/LALStdlib.h>
#include <lal/StringInput.h>
#include <lal/SkyCoordinates.h>
#include <lal/Date.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

NRCSID( PULSARCATH, "$Id$" );

/**
 \name Error Codes */ /*@{*/
#define PULSARCATH_ENUL   1
#define PULSARCATH_EOUT   2
#define PULSARCATH_EMEM   3
#define PULSARCATH_EPARSE 4

#define PULSARCATH_MSGENUL   "Unexpected null pointer in arguments"
#define PULSARCATH_MSGEOUT   "Output handle points to a non-null pointer"
#define PULSARCATH_MSGEMEM   "Memory allocation error"
#define PULSARCATH_MSGEPARSE "Error parsing input file"
/*@}*/


/** Structure \c CompanionNode

This structure stores the orbital parameters of a companion
to a pulsar in a multiple system.  If there is more than one
companion, these structures form a linked list.

<dl>
<dt><tt>LIGOTimeGPS epoch</tt></dt><dd> Epoch of companion periapsis.</dd>

<dt><tt>REAL8 x</tt></dt><dd> Projected orbital semimajor axis \f$(a/c)\sin
i\f$, in seconds.</dd>

<dt><tt>REAL8 period</tt></dt><dd> Orbital period, in seconds, measured at
\c epoch.</dd>

<dt><tt>REAL8 periodDot</tt></dt><dd> First time derivative of orbital
period (dimensionless).</dd>

<dt><tt>REAL8 omega</tt></dt><dd> Longitude of periapsis, in radians,
measured at \c epoch.</dd>

<dt><tt>REAL8 omegaDot</tt></dt><dd> Rate of advance of periapsis, in
radians/s.</dd>

<dt><tt>REAL8 ecc</tt></dt><dd> Orbital eccentricity.</dd>

<dt><tt>REAL8 gamma</tt></dt><dd> Post-Keplerian "gamma" term, in seconds.</dd>

<dt><tt>REAL8 sin</tt></dt><dd> Post-Keplerian "s" term.</dd>

<dt><tt>REAL8 r</tt></dt><dd> Post-Keplerian "r" term.</dd>

<dt><tt>CompanionNode *next</tt></dt><dd> Pointer to next companion's data;
\c NULL if there are no further companions in the system.</dd>
</dl>

*/
typedef struct tagCompanionNode{
  LIGOTimeGPS epoch; /* epoch of periapsis */
  REAL8 x;           /* projected semimajor axis (a/c)sin(i), (s) */
  REAL8 p;           /* orbital period at epoch (s) */
  REAL8 pDot;        /* time derivative of period */
  REAL8 w;           /* longitude of periapsis at epoch (rad) */
  REAL8 wDot;        /* rate of advance of periapsis (rad/s) */
  REAL8 ecc;         /* orbital eccentricity */
  REAL8 gamma;       /* post-Keplerian gamma term (s) */
  REAL8 sin;         /* post-Keplerian s term */
  REAL8 r;           /* post-Keplerian r term */
  struct tagCompanionNode *next; /* pointer to next companion */
} CompanionNode;

/** Structure \c PulsarCatNode

This structure represents a single node in a linked list of
pulsar data, storing data for a single pulsar.  The fields are:

<dl>
<dt><tt>CHAR bname[10]</tt></dt><dd> The B1950 pulsar name (e.g.\
<tt>B0021-72C</tt>), terminated by a <tt>'\0'</tt> character.</dd>

<dt><tt>CHAR jname[12]</tt></dt><dd> The J2000 pulsar name (e.g.\
<tt>J0024-7203U</tt>), terminated by a <tt>'\0'</tt> character.</dd>

<dt><tt>SkyPosition pos</tt></dt><dd> The J2000 pulsar position, in radians.</dd>

<dt><tt>SkyPosition dpos</tt></dt><dd> Uncertainty in \c pos, in
radians.</dd>

<dt><tt>SkyPosition pm</tt></dt><dd> The pulsar proper motion, in radians
per second.</dd>

<dt><tt>SkyPosition dpm</tt></dt><dd> Uncertainty in \c pm, in radians
per second.</dd>

<dt><tt>LIGOTimeGPS posepoch</tt></dt><dd> The epoch of the postion
measurement.</dd>

<dt><tt>REAL8Vector *f</tt></dt><dd> The pulsar spin frequency
<tt>f->data[0]</tt>, and its time derivatives
<tt>f->data[1]</tt>\f$\ldots\f$<tt>f->data[</tt>\f$k\f$<tt>]</tt>\f$\ldots\f$, in units
of Hz\f${}^{k+1}\f$.</dd>

<dt><tt>REAL8Vector *df</tt></dt><dd> The uncertainty in the frequency and
its time derivatives, in the same units.</dd>

<dt><tt>LIGOTimeGPS fepoch</tt></dt><dd> The epoch of the spin and phase
measurements.</dd>

<dt><tt>REAL4 dist</tt></dt><dd> Distance to pulsar, in m.  If negative,
only a lower or upper limit has been established.</dd>

<dt><tt>REAL4 dmin</tt></dt><dd> Lower-limit distance to pulsar, in m.  If
negative, no lower limit has been specified.</dd>

<dt><tt>REAL4 dmax</tt></dt><dd> Upper-limit distance to pulsar, in m.  If
negative, no upper limit has been specified.</dd>

<dt><tt>CHAR lcode</tt></dt><dd> Reliability of distance measurement on low
side, from <tt>'a'</tt> (best) to <tt>'d'</tt> (worst).</dd>

<dt><tt>CHAR ucode</tt></dt><dd> Reliability of distance measurement on high
side, from <tt>'a'</tt> (best) to <tt>'d'</tt> (worst).</dd>

<dt><tt>CompanionNode *companion</tt></dt><dd> Pointer to head of linked
list of orbital parameters for other components of a multiple system;
\c NULL if the pulsar has no known companion.  See below for the
contents of these data nodes.</dd>

<dt><tt>UINT2 typecode</tt></dt><dd> Binary code for additional pulsar
properties.  The typecode is the logical "`or" (i.e.\ the numerical
sum) of the following property codes:
<ul></dd>
<dt>1</dt><dd> Globular cluster association</dd>
<dt>2</dt><dd> Supernova remnant association</dd>
<dt>4</dt><dd> Glitches in period</dd>
<dt>8</dt><dd> Binary or multiple pulsar</dd>
<dt>16</dt><dd> Millisecond pulsar</dd>
<dt>32</dt><dd> Recycled pulsar</dd>
<dt>64</dt><dd> Radio interpulse</dd>
<dt>128</dt><dd> Optical, xray, or gamma-ray pulsed emission
</ul></dd>

<dt><tt>PulsarCatNode *next</tt></dt><dd> Next pulsar in the catalogue's
linked list; \c NULL if this is the last (or only) pulsar in the
list.</dd>
</dl>

*/
typedef struct tagPulsarCatNode {
  CHAR bname[10];   /* B1950 pulsar name */
  CHAR jname[12];   /* J2000 pulsar name */
  SkyPosition pos;  /* J2000 pulsar position */
  SkyPosition dpos; /* J2000 pulsar position uncertainty */
  SkyPosition pm;   /* pulsar proper motion (rad/s) */
  SkyPosition dpm;  /* pulsar proper motion uncertainty (rad/s) */
  LIGOTimeGPS posepoch; /* epoch of the postion measurement */
  REAL8Vector *f;       /* spin frequency and its time derivatives */
  REAL8Vector *df;      /* uncertainty in frequency and derivatives */
  LIGOTimeGPS fepoch;   /* epoch of the spin measurements */
  REAL4 dist; /* distance to pulsar (m) */
  REAL4 dmin; /* lower-limit distance to pulsar (m) */
  REAL4 dmax; /* upper-limit distance to pulsar (m) */
  CHAR lcode; /* reliability of distance on low side ('a' to 'd'@) */
  CHAR ucode; /* reliability of distance on high side ('a' to 'd'@) */
  CompanionNode *companion;      /* linked list of companion data */
  UINT2 typecode;                /* code for additional properties */
  struct tagPulsarCatNode *next; /* next node in list */
} PulsarCatNode;

/** Enumeration \c PulsarCatIndex

This enumerated type is used to give a default ordering to
the fields in the pulsar catalogue.  This is used, for instance, when
reading pulsar catalogue data from a file.  The values are of the form
\c PULSARCATINDEX_\f$\langle\mathrm{label}\rangle\f$, where the
(currently) allowed values of \f$\langle\mathrm{label}\rangle\f$ are:

<table><tr><td>
\c NAME</td><td>pulsar name</td><td></td><td></td></tr>
<tr><td>\c RAJ</td><td>J2000 right ascension</td><td>\c RAJERR</td><td>its uncertainty</td></tr>
<tr><td>\c DECJ</td><td>J2000 declination</td><td>\c DECJERR</td><td>its uncertainty</td></tr>
<tr><td>\c PMRA</td><td>right ascension proper motion</td><td>\c PMRAERR</td><td>its uncertainty</td></tr>
<tr><td>\c PMDEC</td><td>declination proper motion</td><td>\c PMDECERR</td><td>its uncertainty</td></tr>
<tr><td>\c POSEPOCH</td><td>position measurement epoch</td><td></td><td></td></tr>
<tr><td>\c F</td><td>spin frequency</td><td>\c FERR</td><td>its uncertainty</td></tr>
<tr><td>\c F1</td><td>spin frequency derivative</td><td>\c F1ERR</td><td>its uncertainty</td></tr>
<tr><td>\c F1</td><td>spin frequency second derivative</td><td>\c F2ERR</td><td>its uncertainty</td></tr>
<tr><td>\c PEPOCH</td><td>spin measurement epoch</td><td></td><td></td></tr>
<tr><td>\c Dist</td><td>distance</td><td></td><td></td></tr>
<tr><td>\c NUM</td><td>number of enum values</td><td></td><td>
</td></tr></table>

*/
typedef enum {
  PULSARCATINDEX_NAME,
  PULSARCATINDEX_RAJ,   PULSARCATINDEX_RAJERR,
  PULSARCATINDEX_DECJ,  PULSARCATINDEX_DECJERR,
  PULSARCATINDEX_PMRA,  PULSARCATINDEX_PMRAERR,
  PULSARCATINDEX_PMDEC, PULSARCATINDEX_PMDECERR,
  PULSARCATINDEX_POSEPOCH,
  PULSARCATINDEX_F,     PULSARCATINDEX_FERR,
  PULSARCATINDEX_F1,    PULSARCATINDEX_F1ERR,
  PULSARCATINDEX_F2,    PULSARCATINDEX_F2ERR,
  PULSARCATINDEX_PEPOCH,
  PULSARCATINDEX_Dist,
  PULSARCATINDEX_NUM
} PulsarCatIndex;

/* Function prototypes. */




void
LALUpdatePulsarCatNode( LALStatus      *status,
			PulsarCatNode  *node,
			LALPlaceAndGPS *detectorTime,
			EphemerisData  *edat );

void
LALUpdatePulsarCat( LALStatus      *status,
		    PulsarCatNode  *head,
		    LALPlaceAndGPS *detectorTime,
		    EphemerisData  *edat );

void
LALDestroyPulsarCat( LALStatus     *status,
		     PulsarCatNode **head );




void
LALReadPulsarCatHead( LALStatus *status,
		      INT4      indx[PULSARCATINDEX_NUM],
		      TokenList *list );

void
LALReadPulsarCatLine( LALStatus     *status,
		      PulsarCatNode *node,
		      TokenList     *list,
		      INT4          indx[PULSARCATINDEX_NUM] );





#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _PULSARCAT_H */
