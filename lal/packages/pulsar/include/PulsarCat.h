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
#include <lal/StringInput.h>
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
#define PULSARCATH_ENUL   1
#define PULSARCATH_EOUT   2
#define PULSARCATH_EMEM   3
#define PULSARCATH_EPARSE 4

#define PULSARCATH_MSGENUL   "Unexpected null pointer in arguments"
#define PULSARCATH_MSGEOUT   "Output handle points to a non-null pointer"
#define PULSARCATH_MSGEMEM   "Memory allocation error"
#define PULSARCATH_MSGEPARSE "Error parsing input file"
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

\item[\texttt{SkyPosition pos}] The J2000 pulsar position, in radians.

\item[\texttt{SkyPosition dpos}] Uncertainty in \verb@pos@, in
radians.

\item[\texttt{SkyPosition pm}] The pulsar proper motion, in radians
per second.

\item[\texttt{SkyPosition dpm}] Uncertainty in \verb@pm@, in radians
per second.

\item[\texttt{LIGOTimeGPS posepoch}] The epoch of the postion
measurement.

\item[\texttt{REAL8Vector *f}] The pulsar spin frequency
\verb@f->data[0]@, and its time derivatives
\verb@f->data[1]@$\ldots$\verb@f->data[@$k$\verb@]@$\ldots$, in units
of Hz${}^{k+1}$.

\item[\texttt{REAL8Vector *df}] The uncertainty in the frequency and
its time derivatives, in the same units.

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

/******************************************************** <lalLaTeX>

\subsubsection*{Enumeration \texttt{PulsarCatIndex}}
\idx[Type]{PulsarCatIndex}

\noindent This enumerated type is used to give a default ordering to
the fields in the pulsar catalogue.  This is used, for instance, when
reading pulsar catalogue data from a file.  The values are of the form
\verb@PULSARCATINDEX_@$\langle\mathrm{label}\rangle$, where the
(currently) allowed values of $\langle\mathrm{label}\rangle$ are:

\idx[Constant]{PULSARCATINDEX\_NAME}
\idx[Constant]{PULSARCATINDEX\_RAJ}
\idx[Constant]{PULSARCATINDEX\_RAJERR}
\idx[Constant]{PULSARCATINDEX\_DECJ}
\idx[Constant]{PULSARCATINDEX\_DECJERR}
\idx[Constant]{PULSARCATINDEX\_PMRA}
\idx[Constant]{PULSARCATINDEX\_PMRAERR}
\idx[Constant]{PULSARCATINDEX\_PMDEC}
\idx[Constant]{PULSARCATINDEX\_PMDECERR}
\idx[Constant]{PULSARCATINDEX\_POSEPOCH}
\idx[Constant]{PULSARCATINDEX\_F}
\idx[Constant]{PULSARCATINDEX\_FERR}
\idx[Constant]{PULSARCATINDEX\_F1}
\idx[Constant]{PULSARCATINDEX\_F1ERR}
\idx[Constant]{PULSARCATINDEX\_F2}
\idx[Constant]{PULSARCATINDEX\_F2ERR}
\idx[Constant]{PULSARCATINDEX\_PEPOCH}
\idx[Constant]{PULSARCATINDEX\_Dist}
\idx[Constant]{PULSARCATINDEX\_NUM}
\medskip\noindent
\begin{tabular}{ll@{\qquad}ll}
\verb@NAME@  & pulsar name & & \\
\verb@RAJ@   & J2000 right ascension & \verb@RAJERR@  & its uncertainty \\
\verb@DECJ@  & J2000 declination     & \verb@DECJERR@ & its uncertainty \\
\verb@PMRA@  & right ascension proper motion & \verb@PMRAERR@  & its uncertainty \\
\verb@PMDEC@ & declination proper motion     & \verb@PMDECERR@ & its uncertainty \\
\verb@POSEPOCH@ & position measurement epoch & & \\
\verb@F@  & spin frequency                   & \verb@FERR@  & its uncertainty \\
\verb@F1@ & spin frequency derivative        & \verb@F1ERR@ & its uncertainty \\
\verb@F1@ & spin frequency second derivative & \verb@F2ERR@ & its uncertainty \\
\verb@PEPOCH@ & spin measurement epoch & & \\
\verb@Dist@ & distance & & \\
\verb@NUM@ & number of enum values & &
\end{tabular}


******************************************************* </lalLaTeX> */

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

/* <lalLaTeX>
\vfill{\footnotesize\input{PulsarCatHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{PulsarCatC}
</lalLaTeX> */
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

/* <lalLaTeX>
\newpage\input{PulsarCatInputC}
</lalLaTeX> */
void
LALReadPulsarCatHead( LALStatus *status,
		      INT4      indx[PULSARCATINDEX_NUM],
		      TokenList *list );

void
LALReadPulsarCatLine( LALStatus     *status,
		      PulsarCatNode *node,
		      TokenList     *list,
		      INT4          indx[PULSARCATINDEX_NUM] );

/* <lalLaTeX>
\newpage\input{PulsarCatTestC}
</lalLaTeX> */

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _PULSARCAT_H */
