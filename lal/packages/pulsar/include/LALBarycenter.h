/********************************** <lalVerbatim file="LALBarycenterHV">
Author: Cutler, C.
$Id$
*********************************** </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{LALBarycenter.h}}
\label{s:LALBarycenter.h}

Provides routines for transforming from arrival time
at detector (GPS) to pulse emission time (TDB); i.e.,
for ``barycentering'' the measured astronomical time series.
 

\subsection*{Synopsis}

\begin{verbatim}
#include "LALBarycenter.h"
\end{verbatim}

\noindent This header covers the routine
\verb@LALBarycenter.c@.

</lalLaTeX> */



#ifndef _LALBARYCENTER_H    /* Protect against double-inclusion */
#define _LALBARYCENTER_H

#include <stdio.h>
#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/DetectorSite.h>

#ifdef  __cplusplus
extern "C" {
#endif

NRCSID (LALBARYCENTERH,"$Id$");

/* <lalErrTable file="LALBarycenterHErrorTable"> */
#define LALBARYCENTERH_ENULL  2

#define LALBARYCENTERH_MSGENULL  "Null input to Barycenter routine."
/* </lalErrTable> */

#define EphemTableDimE 2500
#define EphemTableDimS 2500
/*Curt: should really allocate less space for EphemTableDimS, but this
  was for testing purposes */ 
/*Curt: rem .h files should not allocate storage */


/* <lalLaTeX>
\subsection*{Error conditions}
\vspace{0.1in}
\input{LALBarycenterHErrorTable}
</lalLaTeX> */


/* Structures. */

/* <lalLaTeX>
\subsection*{Structures}

\begin{verbatim}
struct EphemerisFilenames
\end{verbatim}
\index{\verb&EphemerisFilenames&}
\noindent This structure contains two pointers
to data files containing arrays of center-of-mass
positions for the Earth and Sun, respectively.
The tables are derived from the JPL ephemeris.


\noindent Files tabulate positions for one calendar year 
(actually, a little more than one year, to deal 
with overlaps).  The first line of each table summarizes
what is in it. Subsequent lines give the time (GPS) and the
Earth's position $(x,y,z)$,
velocity $(v_x, v_y, v_z)$, and acceleration $(a_x, a_y, a_z)$
at that instant.  All units of seconds; e.g. positiions have
units of seconds, and accelerations have units 1/sec.

\begin{description}
\item[\texttt{CHAR *earthEphemeris}] File containing Earth's position. 

\item[\texttt{CHAR *sunEphemeris}] File containing Sun's position. 
\end{description}


\begin{verbatim}
struct EphemerisData
\end{verbatim}
\index{\verb&EphemerisData&}
\noindent This structure contains all information about the
center-of-mass positions of the Earth and Sun, listed at regular
time intervals. The fields are

\begin{description}
\item[\texttt{EphemerisFilenames ephiles}] Stucture giving names of
the two files containing positions of Earth and Sun, resp., at evenly spaced times.
\item[\texttt{INT2  leap}]  The number of leap seconds that have
been inserted into UTC between Jan. 6, 1980 (= start of GPS calendar) 
and Jan. 2 of year covered by this ephemeris file; 
e.g. leap = 13 for year 2000.

\item[\texttt{REAL8 gpsE[EphemTableDimE]}] Instant of time, labelled by 
seconds in GPS.

\item[\texttt{REAL8 position[EphemTableDimE][3]}] x,y,z components
of Earth's position at that instant (with respect to SSB, in the ICRS J2000 
frame).  Based on JPL DE405 ephemerid. Units are seconds.

\item[\texttt{REAL8 velocity[EphemTableDimE][3]}] x,y,z components
of Earth's velocity at that instant. Dimensionless (c=1).

\item[\texttt{REAL8 acceleration[EphemTableDimE][3]}] x,y,z components
of Earth's acceleration at that instant. Units are 1/sec.

\item[\texttt{REAL8 gpsS[EphemTableDimS]}] Instant of time, labelled by 
seconds in GPS.

\item[\texttt{REAL8 sunPos[EphemTableDimS][3]}] x,y,z components
of Sun's position  at that instant. Units are seconds.

\item[\texttt{REAL8 sunVel[EphemTableDimS][3]}] x,y,z components
of Earth's velocity at that instant. Dimensionless (c=1).

\item[\texttt{REAL8 sunAccel[EphemTableDimS][3]}] x,y,z components
of Sun's acceleration at that instant. Units are 1/sec. 
\end{description}
  
\begin{verbatim}
struct EarthState
\end{verbatim}
\index{\verb&EarthState&}
\noindent Basic output structure of LALBarycenterEarth.c.
\begin{description}
\item[\texttt{REAL8  einstein}] the einstein delay equiv TDB - TDT 
\item[\texttt{REAL8 deinstein}] d(einstein)/d(tgps) 

\item[\texttt{REAL8 posNow[3]}] Cartesian coords of Earth's center at tgps, extrapolated from JPL DE405 ephemeris; units= sec.

\item[\texttt{REAL8 velNow[3]}] dimensionless velocity of Earth's center at 
tgps, extrapolated from JPL DE405 ephemeris 

\item[\texttt{REAL8 gastRad}] Greenwich Apparent Sidereal Time, 
in radians, at tgps. It's basically the angle thru which Earth has 
spun at given time. gast is like gmst, but has additional correction 
for short-term nutation. 

\item[\texttt{REAL8 M}] variable describing effect of lunisolar precession, at tgps
\item[\texttt{REAL8 N}] variable describing effect of lunisolar precession, at tgps
\item[\texttt{REAL8 delpsi}] variable describing effect of Earth nutation, at tgps
\item[\texttt{REAL8 deleps}] variable describing effect of Earth nutation, at tgps
\item[\texttt{REAL8 se[3]}] vector that points from Sun to Earth at instant tgps, in DE405 coords; units = sec 
\item[\texttt{REAL8 dse[3]}] d(se[3])/d(tgps). Dimensionless
\item[\texttt{REAL8 rse}] length of vector se[3]; units = sec 
\item[\texttt{REAL8 drse}] d(rse)/d(tgps); dimensionless 
\end{description}

\begin{verbatim}
struct BarycenterInput
\end{verbatim}
\index{\verb&BarycenterInput&}
\noindent Basic input structure to LALBarycenter.c.
\begin{description}
\item[\texttt{LIGOTimeGPS  tgps}] input GPS arrival time. 

\item[\texttt{LALDetector site}]  detector site info.

\item[\texttt{REAL8 alpha}] Source right ascension in ICRS J2000 
coords (radians).

\item[\texttt{REAL8 delta}] Source declination in ICRS J2000 coords (radians).

\item[\texttt{REAL8 dInv}] 1/(distance to source), in 1/sec.
\end{description}

\begin{verbatim}
struct EmissionTime
\end{verbatim}
\index{\verb&EmissionTime&}
\noindent Basic output structure produced by LALBarycenter.c.
\begin{description}
\item[\texttt{  REAL8  deltaT}] $t_e$(TDB) - $t_a$(GPS) (+ constant = ``light-travel-time from source to SSB'') 
                               

\item[\texttt{  REAL8 te}]   pulse emission time $t_e$ in TDB (plus constant =
``light-travel-time from source to SSB''), in format of LIGOTImeGPS structure.

\item[\texttt{  REAL8 tDot}]   d(emission time in TDB)/d(arrival time in GPS)  

\item[\texttt{  REAL8 rDetector[3]}]  Cartesian coords (0=x,1=y,2=z) of 
detector position at $t_a$ (GPS), in ICRS J2000 coords. Units = sec. 

\item[\texttt{  REAL8 vDetector[3]}] Cartesian coords (0=x,1=y,2=z) of 
detector velocity at $t_a$ (GPS), in ICRS J2000 coords. Dimensionless. 
\end{description}
</lalLaTeX> */


typedef struct
tagEphemerisFilenames
{
   CHAR *earthEphemeris; 
   CHAR *sunEphemeris;
}
EphemerisFilenames;


typedef struct
tagEphemerisData
{
  EphemerisFilenames ephiles;
  INT2  leap;
  REAL8 gpsE[EphemTableDimE];
  REAL8 position[EphemTableDimE][3];
  REAL8 velocity[EphemTableDimE][3];
  REAL8 acceleration[EphemTableDimE][3];
  REAL8 gpsS[EphemTableDimS];
  REAL8 sunPos[EphemTableDimS][3];
  REAL8 sunVel[EphemTableDimS][3];
  REAL8 sunAccel[EphemTableDimS][3];
}
EphemerisData;

typedef struct
tagEarthState
{
  REAL8  einstein; /*the einstein delay equiv TDB - TDT */
  REAL8 deinstein; /*d(einstein)/d(tgps) */

  REAL8 posNow[3]; /* Cartesian coords of Earth's center at tgps, 
                       extrapolated from JPL DE405 ephemeris; units= sec */
  REAL8 velNow[3];  /* dimensionless velocity of Earth's center at tgps, 
                       extrapolated from JPL DE405 ephemeris */

  REAL8 gastRad;    /*Greenwich Apparent Sidereal Time, in radians, at tgps;
		      Is basically the angle thru which Earth has spun at 
                      given time. gast is like gmst, but has 
                      additional correction for short-term nutation */
  REAL8 M;  /*variable describing effect of lunisolar precession, at tgps*/
  REAL8 N;  /*variable describing effect of lunisolar precession, at tgps*/
  REAL8 delpsi;  /*variable describing effect of Earth nutation, at tgps*/
  REAL8 deleps;  /*variable describing effect of Earth nutation, at tgps*/

  REAL8 se[3];      /*vector that points from Sun to Earth at instant tgps, 
                      in DE405 coords; units = sec */
  REAL8 dse[3];      /*d(se[3])/d(tgps). Dimensionless*/
  REAL8 rse;         /*length of vector se[3]; units = sec */
  REAL8 drse;        /* d(rse)/d(tgps); dimensionless */
}
EarthState;


typedef struct
tagBarycenterInput
{
  LIGOTimeGPS  tgps; /*input GPS arrival time. I use tgps (lower case)
                       to remind that here the LAL structure is a
                       field in the larger structure BarycenterInput.
                       I use tGPS as an input structure (by itself) to 
                       LALBarycenterEarth */
                       
  LALDetector site;      /*DetectorSite structure*/

  REAL8 alpha;      /* source right ascension in ICRS 
                        J2000 coords (radians)*/
  REAL8 delta;      /* source declination in ICRS J2000 coords (radians)*/
  REAL8 dInv;       /* 1/(distance to source), in 1/sec 
                       This is needed to calculate the parallax for very 
                       nearby sources, but that part of code not yet  
                       implemented. */
}
BarycenterInput;


typedef struct
tagEmissionTime
{
  REAL8 deltaT;   /* $t_e$(TDB) - $t_a$(GPS)  
                    +(light-travel-time from source to SSB) */

  LIGOTimeGPS te;/* pulse emission time (TDB); also sometimes called
		    ``arrival time (TDB) of same wavefront at SSB'' */
  REAL8 tDot;    /* d(emission time in TDB)/d(arrival time in GPS)  */

  REAL8 rDetector[3]; /* Cartesian coords (0=x,1=y,2=z) of detector position
                        at $t_a$ (GPS), in ICRS J2000 coords. Units = sec. */

  REAL8 vDetector[3]; /* Cartesian coords (0=x,1=y,2=z) of detector velocity
                        at $t_a$ (GPS), in ICRS J2000 coords. Dimensionless. */
}
EmissionTime;


/*Curt: I should probably take the 1.0 OUT of tDot--ie., output tDot-1 
Is that right, or will users just immediately add back the one anyway??
*/

/*say output REAL8 T is ``time pulse would arrive at a GPS clock 
way out in empty space, if you renormalized  and zero-ed the latter
to give, on average, the same arrival time as the GPS clock on Earth'' */


/* Function prototypes. */

void LALBarycenterEarth(LALStatus *, EarthState *, LIGOTimeGPS *, EphemerisData *);

void LALBarycenter(LALStatus *, EmissionTime *, BarycenterInput *, EarthState *);


/********************************************************** <lalLaTeX>
\vfill{\footnotesize\input{LALBarycenterHV}}
\newpage\input{LALBarycenterC}
******************************************************* </lalLaTeX> */

/* Test program. */

/* <lalLaTeX>
\newpage\input{LALInitBarycenterH}
\newpage\input{LALBarycenterTestC}
</lalLaTeX> */

#ifdef  __cplusplus
}
#endif      /* Close C++ protection */

#endif      /* Close double-include protection */























