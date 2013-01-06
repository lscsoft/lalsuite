/*
*  Copyright (C) 2007 David M. Whitbeck, Jolien Creighton, Ian Jones, Reinhard Prix, Teviet Creighton
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


#ifndef _PULSARTIMES_H
#define _PULSARTIMES_H

#include <lal/LALStdlib.h>
#include <lal/LALBarycenter.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
   \defgroup PulsarTimes_h  Header PulsarTimes.h
   \author Creighton, T. D.
   \ingroup pkg_pulsarCommon
   \brief Provides routines to transform among various time coordinates used in a pulsar search.

   \heading{Synopsis}
   \code
   #include <lal/PulsarTimes.h>
   \endcode

This module covers routines that computes time coordinate
transformations, and derivatives of these transformations with respect
to their parameters.
The motivation is to provide a number of useful
transformations for doing pulsar searches.  For instance, one
transformation might take you from the time measured at a point on
Earth to the time measured in an inertial frame, synchronized
according to a signal arriving from a particular direction in space.
Another might transform from the inertial time coordinate to the
proper time coordinate of a pulsar in a binary system.  Another might
gradually stretch the proper time coordinate of a pulsar that is
spinning down, so that it appears to be at a constant frequency.
Other transformations might be compositions of these, allowing you to
go straight from detector time to the canonical rotation-defined time
of the pulsar.

Mathematically, the transformation from one time \f$t\f$ to another \f$\tau\f$
is written as a function \f$\tau(t)\f$.  In general this function will
depend on other parameters, such as the right ascension and
declination of the source on the sky, or the latitude and longitude of
the observer.  Since in pulsar searches one is often concerned with
how the transformation depends on these parameters, it is necessary to
specify which parameters are allowed to vary and which will be treated
as constants.  We write the transformation as:
\f[
\tau(t,\vec\lambda,\{p\}) \; ,
\f]
where \f$\vec\lambda=(\lambda^1,\ldots,\lambda^n)\f$ are the parameters
that we will allow to vary, while \f$\{p\}\f$ are the parameters that will
be treated as constant.  As the notation suggests, the variable
parameters must be representable in a real vector space, while the
constant parameters need not be real numbers; they may be integers,
names, flags, or anything else required by the transformation
function.

The modules under this header will typically provide function pairs of
the form:
\code
void LALTau( LALStatus             *status,
             REAL8                 *tau,
             REAL8Vector           *variables,
             PulsarTimesParamStruc *constants );

void LALDTau( LALStatus             *status,
              REAL8Vector           *dTau,
              REAL8Vector           *variables,
              PulsarTimesParamStruc *constants );
\endcode
The actual function names will be different; these are just examples.
The function <tt>LALTau()</tt> computes the transformation, while
<tt>LALDTau()</tt> computes the transformation \e and its
derivatives with respect to the parameters in <tt>*variables</tt>.  The
arguments are described below:
<dl>
<dt><tt>*status</tt></dt><dd> This is the universal status structure required
by all LAL functions.</dd>

<dt><tt>*tau</tt></dt><dd> This stores the returned value of
\f$\tau(t,\vec\lambda,\{p\})\f$.</dd>

<dt><tt>*variables</tt></dt><dd> This is a length \f$n+1\f$ vector storing the
arguments of \f$\tau\f$ that are considered to be ``variable''; that is,
\f$t\f$ and \f$\vec\lambda\f$.  They are stored as follows:<br>
	<tt>variables->data[0]</tt>\f$=t\f$,<br>
        <tt>variables->data[</tt>\f$i\f$<tt>]</tt>\f$=\lambda^i\f$, where \f$i=1,\ldots,n\f$.<br>
</dd>

<dt><tt>*constants</tt></dt><dd> This stores the constant parameters
\f$\{p\}\f$, in a format described in the Structures section below.</dd>

<dt><tt>*dTau</tt></dt><dd>This is a length \f$n+2\f$ vector storing the value
of \f$\tau(t,\vec\lambda,\{p\})\f$ and its derivatives with respect to its
variable arguments, in the following format:<br>
	<tt>dTau->data[0]</tt>\f$=\tau\f$<br>
        <tt>dTau->data[1]</tt>\f$=\partial\tau/\partial t\f$<br>
        <tt>dTau->data[</tt>\f$i+1\f$<tt>]</tt>\f$=\partial\tau/
	\partial\lambda^i\f$, where \f$i=1,\ldots,n\f$.<br>
</dd>
</dl>

It may seem redundant that both <tt>LALTau()</tt> and <tt>LALDTau()</tt>
compute and return the value of \f$\tau\f$, especially since returning
\f$\tau\f$ in <tt>*dTau</tt> messes up an otherwise elegant indexing scheme.
The reason is that many of the calculations involved in computing the
derivatives of \f$\tau\f$ are also used in calculating \f$\tau\f$ itself, and
it would be inefficient to have to repeat them by calling both
functions.  <tt>LALTau()</tt> is provided simply as a shortcut for those
occasions when you do not need the derivatives.

It is also worth noting that different pulsar searches may involve the
same transformations but vary different sets of parameters.  There are
two approches to dealing with this.  The simplest is usually to use a
function that treats all of the parameters as variable, and then
ignore (or set to zero) the derivatives that aren't used.  The more
computationally efficient approach is to write separate pairs of
routines that place different parameters in <tt>*variables</tt> and
<tt>*constants</tt>, although this may require recompiling the library.

*/
/*@{*/

/** \name Error Codes */
/*@{*/
/** \ingroup PulsarTimes_h */
#define PULSARTIMESH_ENUL 1
#define PULSARTIMESH_EBAD 2

#define PULSARTIMESH_MSGENUL "Null pointer"
#define PULSARTIMESH_MSGEBAD "Bad parameter values"
/*@}*/

/**
 * This structure stores a superset of all constant parameters
 * required by the functions provided by the header PulsarTimes.h.
 * Although the structure is quite large, it is
 * supposed to store only constants, so only one structure need be
 * allocated for each transformation type used in the overall algorithm.
 * It is passed by pointer, so its size entails no overhead.  A basic
 * description of the current fields is given below, but a detailed
 * discussion of those fields must be deferred to the documentation of
 * the individual modules that use those particular fields.
 */
typedef struct tagPulsarTimesParamStruc {
  LIGOTimeGPS epoch; 	/**< A reference detector time; all
                         * other times in the transformation are represented as \c REAL8
                         * numbers, giving the time in seconds since #epoch.
                         */

  REAL8 t0;		/**<  A reference time for a particular transformation,
                         * normally defined such that \f$\tau(\f$\c t0\f$)=0\f$.
                         */

  REAL8 tAutumn; 	/**< Time of the first autumnal equinox following #epoch. */

  REAL8 tMidnight; 	/**< Time of the first sidereal midnight following #epoch. */

  REAL8 latitude; 	/**< Detector north latitude, in radians. */
  REAL8 longitude; 	/**< Detector east longitude (i.e.\ counterclockwise about the north pole), in radians. */

  const EphemerisData *ephemeris; /**< Ephemeris data containing positions, velocities, etc... of Earth and Sun
                                   * for the year under consideration. */

  const LALDetector *site;        /**< The particular detector under consideration. */

  /** \name Transformation composition.
   * The following fields are used by the module TComp.c, which
   * composes two transformations \f$t_1(t)\f$ and \f$t_2(t)\f$ into an overall
   * transformation \f$t_c(t)=t_2(t_1(t))\f$. */
  /*@{ */
  void (*t1)( LALStatus *, REAL8 *, REAL8Vector *, struct tagPulsarTimesParamStruc * ); /**< The first of the pair of transformations to be composed */

  void (*t2)( LALStatus *, REAL8 *, REAL8Vector *, struct tagPulsarTimesParamStruc * ); /**< The second of the pair of transformations to be composed */

  void (*dt1)( LALStatus *, REAL8Vector *, REAL8Vector *, struct tagPulsarTimesParamStruc * ); /**< The time derivative function corresponding to t1() */

  void (*dt2)( LALStatus *, REAL8Vector *, REAL8Vector *, struct tagPulsarTimesParamStruc * ); /**< The time derivative function corresponding to t2() */

  struct tagPulsarTimesParamStruc *constants1;	/**< The constant parameters used by t1() */

  struct tagPulsarTimesParamStruc *constants2;	/**< The constant parameters used by t2() */

  UINT4 nArgs; 		/**< The number of variable parameters \f$\lambda^k\f$ to be sent to the function t1() */
  /*@} */
} PulsarTimesParamStruc;

int XLALGetEarthTimes( const LIGOTimeGPS *tepoch, REAL8 *tMidnight, REAL8 *tAutumn );

/** \cond DONT_DOXYGEN */

/* Function prototypes. */
void
LALGetEarthTimes( LALStatus *, PulsarTimesParamStruc *times );




void
LALTBaryPtolemaic( LALStatus             *,
		   REAL8                 *tBary,
		   REAL8Vector           *variables,
		   PulsarTimesParamStruc *constants );

void
LALDTBaryPtolemaic( LALStatus             *,
		    REAL8Vector           *dtBary,
		    REAL8Vector           *variables,
		    PulsarTimesParamStruc *constants );




void
LALTSpin( LALStatus             *,
	  REAL8                 *tSpin,
	  REAL8Vector           *variables,
	  PulsarTimesParamStruc *constants );

void
LALDTSpin( LALStatus             *,
	   REAL8Vector           *dtSpin,
	   REAL8Vector           *variables,
	   PulsarTimesParamStruc *constants );




void
LALTComp( LALStatus             *,
	  REAL8                 *tComp,
	  REAL8Vector           *variables,
	  PulsarTimesParamStruc *constants );

void
LALDTComp( LALStatus             *,
	   REAL8Vector           *dtComp,
	   REAL8Vector           *variables,
	   PulsarTimesParamStruc *constants );




void
LALDTEphemeris( LALStatus             *,
	        REAL8Vector           *tBary,
	        REAL8Vector           *variables,
	        PulsarTimesParamStruc *constants );

void
LALTEphemeris(LALStatus *,
	      REAL8 *tBary,
	      REAL8Vector *variables,
	      PulsarTimesParamStruc *constants);

/** \endcond */

/*@}*/

#ifdef __cplusplus
}
#endif

#endif
