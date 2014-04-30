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
 * \author Creighton, T. D.
 * \file
 * \ingroup pulsarTODO
 *
 * \brief Provides structures and routines to store and manipulate pulsar properties.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/PulsarCat.h>
 * \endcode
 *
 * This header covers structures to store pulsar properties in
 * a standard format, and routines to manipulate and update these
 * properties.  The set of properties stored in the catalogue is based on
 * radio pulsar catalogues, with some additions and subtractions specific
 * to gravitational wave observations.  The list of properties can be
 * expanded in future by adding more fields to the structure.
 *
 * All properties are those that would be measured by an observer at the
 * solar system barycentre.  For properties that depend on time
 * (e.g.\ time-varying position, periods, etc.), an epoch is specified.
 * The properties are then those that would be observed at the specified
 * instant of time at the solar system barycentre; i.e.\ when the wave
 * fronts carrying that information pass the solar system barycentre.
 *
 * \par A note on companion orbits:
 * Several known pulsars exist in
 * multiple systems, and radio-pulsar catalogues include detailed models
 * of the companion orbits, as determined from the pulsar timing.  See
 * \ref GenerateSpinOrbitCW.h in the \c inject package for a
 * discussion of the parameters defining the orientation of a companion
 * orbit.
 *
 * Radio-pulsar observations rarely determine the inclination \f$i\f$ of the
 * orbit to the sky plane, and thus cannot resolve the longitude of the
 * ascending node \f$\Omega\f$ and the argument of the periapsis \f$\omega\f$ as
 * independent parameters.  Instead, they list the longitude of the
 * periapsis \f$w\f$, which is the angle in the plane of the sky from the
 * North direction towards the West direction, to the ray from the system
 * barycentre to the periapsis projected onto the plane of the sky.  If
 * any three of \f$i\f$, \f$\Omega\f$, \f$\omega\f$, and \f$w\f$ are known, the fourth
 * can be determined from the relation:
 * \f[
 * w - \Omega = \arctan\;2(\sin\omega\cos i,\cos\omega) \;,
 * \f]
 * or equivalently:
 * \f[
 * \omega = \arctan\;2(\cos[w-\Omega],\sin[w-\Omega]/\cos i) \;.
 * \f]
 *
 * In addition to these Keplerian orbital parameters, some radio-pulsar
 * systems have measured post-Keplerian relativistic orbital parameters.
 * Some of these are obvious: \f$\dot{w}\f$ and \f$\dot{P}\f$ are the rate of
 * change in \f$w\f$ (periapsis precession) and the orbital period \f$P\f$ (due
 * to gravitational radiation reaction).  The catalogue also lists
 * post-Keplerian parameters \"sin\" and \"r\", whose meanings I don't
 * know.
 *
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

/** \name Error Codes */ /*@{*/
#define PULSARCATH_ENUL   1
#define PULSARCATH_EOUT   2
#define PULSARCATH_EMEM   3
#define PULSARCATH_EPARSE 4

#define PULSARCATH_MSGENUL   "Unexpected null pointer in arguments"
#define PULSARCATH_MSGEOUT   "Output handle points to a non-null pointer"
#define PULSARCATH_MSGEMEM   "Memory allocation error"
#define PULSARCATH_MSGEPARSE "Error parsing input file"
/*@}*/


/**
 * This structure stores the orbital parameters of a companion
 * to a pulsar in a multiple system.  If there is more than one
 * companion, these structures form a linked list.
 */
typedef struct tagCompanionNode{
  LIGOTimeGPS epoch; /**< Epoch of companion periapsis */
  REAL8 x;           /**< Projected orbital semimajor axis \f$(a/c)\sin i\f$, in seconds */
  REAL8 p;           /**< Orbital period, in seconds, measured at \c epoch */
  REAL8 pDot;        /**< First time derivative of orbital period (dimensionless) */
  REAL8 w;           /**< Longitude of periapsis, in radians, measured at \c epoch */
  REAL8 wDot;        /**< Rate of advance of periapsis, in radians/s */
  REAL8 ecc;         /**< Orbital eccentricity */
  REAL8 gamma;       /**< Post-Keplerian \"gamma\" term, in seconds */
  REAL8 sin;         /**< Post-Keplerian \"s\" term */
  REAL8 r;           /**< Post-Keplerian \"r\" term */
  struct tagCompanionNode *next; /**< Pointer to next companion's data; \c NULL if there are no further companions in the system. */
} CompanionNode;

/**
 * This structure represents a single node in a linked list of
 * pulsar data, storing data for a single pulsar.
 */
typedef struct tagPulsarCatNode {
  CHAR bname[10];   	/**< The B1950 pulsar name (e.g.\ <tt>B0021-72C</tt>), terminated by a <tt>'\\0'</tt> character */
  CHAR jname[12];   	/**< The J2000 pulsar name (e.g.\ <tt>J0024-7203U</tt>), terminated by a <tt>'\\0'</tt> character */
  SkyPosition pos;  	/**< The J2000 pulsar position, in radians */
  SkyPosition dpos; 	/**< Uncertainty in \c pos, in radians */
  SkyPosition pm;   	/**< The pulsar proper motion, in radians per second */
  SkyPosition dpm;  	/**< Uncertainty in \c pm, in radians per second */
  LIGOTimeGPS posepoch; /**< The epoch of the postion measurement */
  REAL8Vector *f;       /**< The pulsar spin frequency <tt>f-\>data[0]</tt>, and its time derivatives
                         * <tt>f-\>data[1]...f-\>data[k]...</tt>, in units of \f$\mathrm{Hz}^{k+1}\f$
                         */
  REAL8Vector *df;      /**< The uncertainty in the frequency and its time derivatives, in the same units */
  LIGOTimeGPS fepoch;   /**< The epoch of the spin and phase measurements */
  REAL4 dist; 		/**< Distance to pulsar, in m;  If negative, only a lower or upper limit has been established */
  REAL4 dmin; 		/**< Lower-limit distance to pulsar, in m; If negative, no lower limit has been specified */
  REAL4 dmax; 		/**< Upper-limit distance to pulsar, in m;  If negative, no upper limit has been specified */
  CHAR lcode; 		/**< Reliability of distance measurement on low side, from <tt>'a'</tt> (best) to <tt>'d'</tt> (worst) */
  CHAR ucode; 		/**< Reliability of distance measurement on high side, from <tt>'a'</tt> (best) to <tt>'d'</tt> (worst) */
  CompanionNode *companion;/**< Pointer to head of linked list of orbital parameters for other components of a multiple system;
                            * \c NULL if the pulsar has no known companion;  See below for the contents of these data nodes
                            */
  UINT2 typecode;	/**< Binary code for additional pulsar properties. The typecode is the logical \"or\" (i.e.\ the numerical sum)
                         * of the following property codes:
                         * <dl>
                         * <dt>1</dt><dd> Globular cluster association</dd>
                         * <dt>2</dt><dd> Supernova remnant association</dd>
                         * <dt>4</dt><dd> Glitches in period</dd>
                         * <dt>8</dt><dd> Binary or multiple pulsar</dd>
                         * <dt>16</dt><dd> Millisecond pulsar</dd>
                         * <dt>32</dt><dd> Recycled pulsar</dd>
                         * <dt>64</dt><dd> Radio interpulse</dd>
                         * <dt>128</dt><dd> Optical, xray, or gamma-ray pulsed emission</dd>
                         * </dl>
                         */
  struct tagPulsarCatNode *next; /**< Next pulsar in the catalogue's linked list; \c NULL if this is the last (or only) pulsar in the list */
} PulsarCatNode;

/**
 * This enumerated type is used to give a default ordering to
 * the fields in the pulsar catalogue.  This is used, for instance, when
 * reading pulsar catalogue data from a file.  The values are of the form
 * \c PULSARCATINDEX_\f$\langle\mathrm{label}\rangle\f$, where the
 * (currently) allowed values of \f$\langle\mathrm{label}\rangle\f$ are:
 *
 * <table>
 * <tr><th>NAME</th><th>pulsar name</th><th></th><th></th></tr>
 * <tr><td>\c RAJ</td><td>J2000 right ascension</td><td>\c RAJERR</td><td>its uncertainty</td></tr>
 * <tr><td>\c DECJ</td><td>J2000 declination</td><td>\c DECJERR</td><td>its uncertainty</td></tr>
 * <tr><td>\c PMRA</td><td>right ascension proper motion</td><td>\c PMRAERR</td><td>its uncertainty</td></tr>
 * <tr><td>\c PMDEC</td><td>declination proper motion</td><td>\c PMDECERR</td><td>its uncertainty</td></tr>
 * <tr><td>\c POSEPOCH</td><td>position measurement epoch</td><td></td><td></td></tr>
 * <tr><td>\c F</td><td>spin frequency</td><td>\c FERR</td><td>its uncertainty</td></tr>
 * <tr><td>\c F1</td><td>spin frequency derivative</td><td>\c F1ERR</td><td>its uncertainty</td></tr>
 * <tr><td>\c F1</td><td>spin frequency second derivative</td><td>\c F2ERR</td><td>its uncertainty</td></tr>
 * <tr><td>\c PEPOCH</td><td>spin measurement epoch</td><td></td><td></td></tr>
 * <tr><td>\c Dist</td><td>distance</td><td></td><td></td></tr>
 * <tr><td>\c NUM</td><td>number of enum values</td><td></td><td></td></tr>
 * </table>
 *
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
