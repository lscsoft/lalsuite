/********************************** <lalVerbatim file="PulsarTimesHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{PulsarTimes.h}}
\label{s:PulsarTimes.h}

Provides routines to transform among various time coordinates used in
a pulsar search.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/PulsarTimes.h>
\end{verbatim}

\noindent This header covers routines that computes time coordinate
transformations, and derivatives of these transformations with respect
to their parameters.  The motivation is to provide a number of useful
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

Mathematically, the transformation from one time $t$ to another $\tau$
is written as a function $\tau(t)$.  In general this function will
depend on other parameters, such as the right ascension and
declination of the source on the sky, or the latitude and longitude of
the observer.  Since in pulsar searches one is often concerned with
how the transformation depends on these parameters, it is necessary to
specify which parameters are allowed to vary and which will be treated
as constants.  We write the transformation as:
$$
\tau(t,\vec\lambda,\{p\}) \; ,
$$
where $\vec\lambda=(\lambda^1,\ldots,\lambda^n)$ are the parameters
that we will allow to vary, while $\{p\}$ are the parameters that will
be treated as constant.  As the notation suggests, the variable
parameters must be representable in a real vector space, while the
constant parameters need not be real numbers; they may be integers,
names, flags, or anything else required by the transformation
function.

The modules under this header will typically provide function pairs of
the form:
\begin{verbatim}
void LALTau( LALStatus             *,
             REAL8                 *tau,
             REAL8Vector           *variables,
             PulsarTimesParamStruc *constants );

void LALDTau( LALStatus             *,
              REAL8Vector           *dTau,
              REAL8Vector           *variables,
              PulsarTimesParamStruc *constants );
\end{verbatim}
The actual function names will be different; these are just examples.
The function \verb@LALTau()@ computes the transformation, while
\verb@LALDTau()@ computes the transformation \emph{and} its
derivatives with respect to the parameters in \verb@*variables@.  The
arguments are described below:
\begin{description}
\item[\texttt{*stat}] This is the universal status structure required
by all LAL functions.

\item[\texttt{*tau}] This stores the returned value of
$\tau(t,\vec\lambda,\{p\})$.

\item[\texttt{*variables}] This is a length $n+1$ vector storing the
arguments of $\tau$ that are considered to be ``variable''; that is,
$t$ and $\vec\lambda$.  They are stored as follows:
	\begin{description}
	\item[\texttt{variables->data[0]}]$=t$,
	\item[\texttt{variables->data[}$i$\texttt{]}]$=\lambda^i$,
	where $i=1,\ldots,n$.
	\end{description}

\item[\texttt{*constants}] This stores the constant parameters
$\{p\}$, in a format described in the Structures section below.

\item[\texttt{*dTau}] This is a length $n+2$ vector storing the value
of $\tau(t,\vec\lambda,\{p\})$ and its derivatives with respect to its
variable arguments, in the following format:
	\begin{description}
	\item[\texttt{dTau->data[0]}]$=\tau$
	\item[\texttt{dTau->data[1]}]$=\partial\tau/\partial t$
	\item[\texttt{dTau->data[}$i+1$\texttt{]}]$=\partial\tau/
	\partial\lambda^i$, where $i=1,\ldots,n$.
	\end{description}
\end{description}

It may seem redundant that both \verb@LALTau()@ and \verb@LALDTau()@
compute and return the value of $\tau$, especially since returning
$\tau$ in \verb@*dTau@ messes up an otherwise elegant indexing scheme.
The reason is that many of the calculations involved in computing the
derivatives of $\tau$ are also used in calculating $\tau$ itself, and
it would be inefficient to have to repeat them by calling both
functions.  \verb@LALTau()@ is provided simply as a shortcut for those
occasions when you do not need the derivatives.

It is also worth noting that different pulsar searches may involve the
same transformations but vary different sets of parameters.  There are
two approches to dealing with this.  The simplest is usually to use a
function that treats all of the parameters as variable, and then
ignore (or set to zero) the derivatives that aren't used.  The more
computationally efficient approach is to write separate pairs of
routines that place different parameters in \verb@*variables@ and
\verb@*constants@, although this may require recompiling the library.

******************************************************* </lalLaTeX> */

#ifndef _PULSARTIMES_H
#define _PULSARTIMES_H

#include <lal/LALStdlib.h>

#include <lal/LALBarycenter.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID(PULSARTIMESH,"$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define PULSARTIMESH_ENUL 1
#define PULSARTIMESH_EBAD 2

#define PULSARTIMESH_MSGENUL "Null pointer"
#define PULSARTIMESH_MSGEBAD "Bad parameter values"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{PulsarTimesParamStruc}}
\idx[Type]{PulsarTimesParamStruc}

\noindent This structure stores a superset of all constant parameters
required by the functions provided by this header
\verb@PulsarTimes.h@.  Although the structure is quite large, it is
supposed to store only constants, so only one structure need be
allocated for each transformation type used in the overall algorithm.
It is passed by pointer, so its size entails no overhead.  A basic
description of the current fields is given below, but a detailed
discussion of those fields must be deferred to the documentation of
the individual modules that use those particular fields.  The fields
are:

\begin{description}
\item[\texttt{LIGOTimeGPS epoch}] A reference detector time.  All
other times in the transformation are represented as \verb@REAL8@
numbers, giving the time in seconds since \verb@epoch@.

\item[\texttt{REAL8 t0}] A reference time for a particular
transformation, normally defined such that $\tau($\verb@t0@$)=0$.

\item[\texttt{REAL8 tAutumn}] Time of the first autumnal equinox
following \verb@epoch@.

\item[\texttt{REAL8 tMidnight}] Time of the first sidereal midnight
following \verb@epoch@.

\item[\texttt{REAL8 latitude}] Detector north latitude, in radians.

\item[\texttt{REAL8 longitude}] Detector east longitude (i.e.\
counterclockwise about the north pole), in radians.

\item[\texttt{EphemerisData ephemeris}] Ephemeris data containing positions, 
velocities, etc... of Earth and Sun for the year under consideration.

\item[\texttt{LALDetector site}] The particular detector under consideration.


\end{description}
The following fields are used by the module \verb@TComp.c@, which
composes two transformations $t_1(t)$ and $t_2(t)$ into an overall
transformation $t_c(t)=t_2(t_1(t))$.
\begin{description}
\item[\texttt{void *t1( LALStatus *, REAL8 *, REAL8Vector *,
PulsarTimesParamStruc * )}] The first of the pair of transformations
to be composed.

\item[\texttt{void *t2( LALStatus *, REAL8 *, REAL8Vector *,
PulsarTimesParamStruc * )}] The second of the pair of transformations
to be composed.

\item[\texttt{void *dt1( LALStatus *, REAL8Vector *, REAL8Vector *,
PulsarTimesParamStruc * )}] The time derivative function corresponding
to \verb@*t1()@.

\item[\texttt{void *dt2( LALStatus *, REAL8Vector *, REAL8Vector *,
PulsarTimesParamStruc * )}] The time derivative function corresponding
to \verb@*t2()@.

\item[\texttt{PulsarTimesParamStruc *constants1}] The constant
parameters used by \verb@*t1()@.

\item[\texttt{PulsarTimesParamStruc *constants2}] The constant
parameters used by \verb@*t2()@.

\item[\texttt{UINT4 nArgs}] The number of variable parameters
$\lambda^k$ to be sent to the function \verb@*t1()@.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagPulsarTimesParamStruc {
  LIGOTimeGPS epoch; /* A reference detector time.  All other times in
			the transformation are represented as REAL8
			numbers, giving the time in seconds since
			epoch. */
  REAL8 t0; /* A reference time for a particular transformation,
	       normally defined such that tau=0. */
  REAL8 tAutumn; /* Time of the first autumnal equinox following
		    epoch. */
  REAL8 tMidnight; /* Time of the first sidereal midnight following
		      epoch */
  REAL8 latitude; /* Detector north latitude, in radians. */
  REAL8 longitude; /* Detector east longitude (i.e. counterclockwise
		      about the north pole), in radians. */
  EphemerisData *ephemeris; /* Ephemeris data containing positions, */
                           /* velocities, etc... of Earth and Sun */
  LALDetector *site;        /* The particular detector under consideration */

  void (*t1)( LALStatus *, REAL8 *, REAL8Vector *,
	      struct tagPulsarTimesParamStruc * );
  /* The first of the pair of transformations to be composed. */
  void (*t2)( LALStatus *, REAL8 *, REAL8Vector *,
	      struct tagPulsarTimesParamStruc * );
  /* The second of the pair of transformations to be composed. */
  void (*dt1)( LALStatus *, REAL8Vector *, REAL8Vector *,
	       struct tagPulsarTimesParamStruc * );
  /* The time derivative function corresponding to *t1(). */
  void (*dt2)( LALStatus *, REAL8Vector *, REAL8Vector *,
	       struct tagPulsarTimesParamStruc * );
  /* The time derivative function corresponding to *t2(). */
  struct tagPulsarTimesParamStruc *constants1;
  /* The constant parameters used by *t1(). */
  struct tagPulsarTimesParamStruc *constants2;
  /* The constant parameters used by *t2(). */
  UINT4 nArgs;
  /* The number of variable parameters to be sent to *t1(). */
} PulsarTimesParamStruc;

/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{GetEarthTimesC}
</lalLaTeX> */
void
LALGetEarthTimes( LALStatus *, PulsarTimesParamStruc *times );

/* <lalLaTeX>
\newpage\input{TBaryPtolemaicC}
</lalLaTeX> */
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

/* <lalLaTeX>
\newpage\input{TSpinC}
</lalLaTeX> */
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

/* <lalLaTeX>
\newpage\input{TCompC}
</lalLaTeX> */
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

/* <lalLaTeX>
\newpage\input{DTEphemerisC}
</lalLaTeX> */
void
LALDTEphemeris( LALStatus             *,
	        REAL8Vector           *tBary,
	        REAL8Vector           *variables,
	        PulsarTimesParamStruc *constants );

void
LALTEphemeris(LALStatus *,
	      REAL8Vector *tBary,
	      REAL8Vector *variables,
	      PulsarTimesParamStruc *constants);

#ifdef __cplusplus
}
#endif

#endif
