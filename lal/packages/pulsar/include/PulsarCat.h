/************************************ <lalVerbatim file="PulsarCatHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\section{Header \texttt{PulsarCat.h}}
\label{s:PulsarCat.h}

Provides structures and routines to store and manipulate pulsar
properties.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/PulsarCat.h>
\end{verbatim}

\noindent This header covers structures to store pulsar properties in
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

\paragraph{A note on companion orbits:} Several known pulsars exist in
multiple systems, and radio-pulsar catalogues include detailed models
of the companion orbits, as determined from the pulsar timing.  See
\verb@GenerateSpinOrbitCW.h@ in the \verb@inject@ package for a
discussion of the parameters defining the orientation of a companion
orbit.

Radio-pulsar observations rarely determine the inclination $i$ of the
orbit to the sky plane, and thus cannot resolve the longitude of the
ascending node $\Omega$ and the argument of the periapsis $\omega$ as
independent parameters.  Instead, they list the longitude of the
periapsis $w$, which is the angle in the plane of the sky from the
North direction towards the West direction, to the ray from the system
barycentre to the periapsis projected onto the plane of the sky.  If
any three of $i$, $\Omega$, $\omega$, and $w$ are known, the fourth
can be determined from the relation:
$$
w - \Omega = \arctan\;2(\sin\omega\cos i,\cos\omega) \;,
$$
or equivalently:
$$
\omega = \arctan\;2(\cos[w-\Omega],\sin[w-\Omega]/\cos i) \;.
$$

In addition to these Keplerian orbital parameters, some radio-pulsar
systems have measured post-Keplerian relativistic orbital parameters.
Some of these are obvious: $\dot{w}$ and $\dot{P}$ are the rate of
change in $w$ (periapsis precession) and the orbital period $P$ (due
to gravitational radiation reaction).  The catalogue also lists
post-Keplerian parameters ``sin'' and ``r'', whose meanings I don't
know.

******************************************************* </lalLaTeX> */

#ifndef _PULSARCAT_H
#define _PULSARCAT_H

#include <lal/LALStdlib.h>
#include <lal/SkyCoordinates.h>
#include <lal/Date.h>

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID( PULSARCATH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define PULSARCATH_ENUL    1
#define PULSARCATH_EOUT    2
#define PULSARCATH_EMEM    3

#define PULSARCATH_MSGENUL    "Unexpected null pointer in arguments"
#define PULSARCATH_MSGEOUT    "Output handle points to a non-null pointer"
#define PULSARCATH_MSGEMEM    "Memory allocation error"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{CompanionNode}}
\idx[Type]{CompanionNode}

\noindent This structure stores the orbital parameters of a companion
to a pulsar in a multiple system.  If there is more than one
companion, these structures form a linked list.

\begin{description}
\item[\texttt{LIGOTimeGPS epoch}] Epoch of companion periapsis.

\item[\texttt{REAL8 x}] Projected orbital semimajor axis $(a/c)\sin
i$, in seconds.

\item[\texttt{REAL8 period}] Orbital period, in seconds, measured at
\verb@epoch@.

\item[\texttt{REAL8 periodDot}] First time derivative of orbital
period (dimensionless).

\item[\texttt{REAL8 omega}] Longitude of periapsis, in radians,
measured at \verb@epoch@.

\item[\texttt{REAL8 omegaDot}] Rate of advance of periapsis, in
radians/s.

\item[\texttt{REAL8 ecc}] Orbital eccentricity.

\item[\texttt{REAL8 gamma}] Post-Keplerian ``gamma'' term, in seconds.

\item[\texttt{REAL8 sin}] Post-Keplerian ``s'' term.

\item[\texttt{REAL8 r}] Post-Keplerian ``r'' term.

\item[\texttt{CompanionNode *next}] Pointer to next companion's data;
\verb@NULL@ if there are no further companions in the system.
\end{description}

******************************************************* </lalLaTeX> */

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

/******************************************************** <lalLaTeX>

\subsubsection*{Structure \texttt{PulsarCatNode}}
\idx[Type]{PulsarCatNode}

\noindent This structure represents a single node in a linked list of
pulsar data, storing data for a single pulsar.  The fields are:

\begin{description}
\item[\texttt{CHAR bname[10]}] The B1950 pulsar name (e.g.\
\verb@B0021-72C@), terminated by a \verb@'\0'@ character.

\item[\texttt{CHAR jname[12]}] The J2000 pulsar name (e.g.\
\verb@J0024-7203U@), terminated by a \verb@'\0'@ character.

\item[\texttt{SkyPosition bpos}] The B1950 pulsar position.

\item[\texttt{SkyPosition jpos}] The J2000 pulsar position.

\item[\texttt{REAL4 pmra}] The proper motion in right ascension, in
radians per year.

\item[\texttt{REAL4 pmdec}] The proper motion in declination, in
radians per year.

\item[\texttt{LIGOTimeGPS posepoch}] The epoch of the postion
measurement.

\item[\texttt{REAL8Vector *f}] The pulsar spin frequency
\verb@f->data[0]@, and its time derivatives
\verb@f->data[1]@$\ldots$\verb@f->data[@$k$\verb@]@$\ldots$, in units
of Hz${}^{k+1}$.

\item[\texttt{LIGOTimeGPS fepoch}] The epoch of the spin and phase
measurements.

\item[\texttt{REAL4 dist}] Distance to pulsar, in m.  If negative,
only a lower or upper limit has been established.

\item[\texttt{REAL4 dmin}] Lower-limit distance to pulsar, in m.  If
negative, no lower limit has been specified.

\item[\texttt{REAL4 dmax}] Upper-limit distance to pulsar, in m.  If
negative, no upper limit has been specified.

\item[\texttt{CHAR lcode}] Reliability of distance measurement on low
side, from \verb@'a'@ (best) to \verb@'d'@ (worst).

\item[\texttt{CHAR ucode}] Reliability of distance measurement on high
side, from \verb@'a'@ (best) to \verb@'d'@ (worst).

\item[\texttt{CompanionNode *companion}] Pointer to head of linked
list of orbital parameters for other components of a multiple system;
\verb@NULL@ if the pulsar has no known companion.  See below for the
contents of these data nodes.

\item[\texttt{UINT2 typecode}] Binary code for additional pulsar
properties.  The typecode is the logical ```or'' (i.e.\ the numerical
sum) of the following property codes:
\begin{itemize}
\item[1] Globular cluster association
\item[2] Supernova remnant association
\item[4] Glitches in period
\item[8] Binary or multiple pulsar
\item[16] Millisecond pulsar
\item[32] Recycled pulsar
\item[64] Radio interpulse
\item[128] Optical, xray, or gamma-ray pulsed emission
\end{itemize}

\item[\texttt{PulsarCatNode *next}] Next pulsar in the catalogue's
linked list; \verb@NULL@ if this is the last (or only) pulsar in the
list.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagPulsarCatNode {
  CHAR bname[10];  /* B1950 pulsar name */
  CHAR jname[12];  /* J2000 pulsar name */
  SkyPosition pos; /* J2000 pulsar position */
  REAL4 pmra;      /* proper motion in right ascension (rad/yr) */
  REAL4 pmdec;     /* proper motion in declination (rad/yr) */
  LIGOTimeGPS posepoch; /* epoch of the postion measurement */
  REAL8Vector *f;       /* spin frequency and its time derivatives */
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


/* <lalLaTeX>
\vfill{\footnotesize\input{PulsarCatHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{PulsarCatC}
</lalLaTeX> */
void
LALUpdatePulsarCatNode( LALStatus      *stat,
			PulsarCatNode  *node,
			LALPlaceAndGPS *detectorTime,
			EphemerisData  *edat );

void
LALUpdatePulsarCat( LALStatus      *stat,
		    PulsarCatNode  *head,
		    LALPlaceAndGPS *detectorTime,
		    EphemerisData  *edat );

void
LALDestroyPulsarCat( LALStatus     *stat,
		     PulsarCatNode **head );

/* <lalLaTeX>
\newpage\input{PulsarCatTestC}
</lalLaTeX> */

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _PULSARCAT_H */
