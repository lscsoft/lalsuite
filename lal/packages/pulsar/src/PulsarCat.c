/************************************ <lalVerbatim file="PulsarCatCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{PulsarCat.c}}
\label{ss:PulsarCat.c}

Manipulates a catalogue of pulsar data.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{PulsarCatCP}
\idx{LALUpdatePulsarCatNode()}
\idx{LALUpdatePulsarCat()}
\idx{LALDestroyPulsarCat()}

\subsubsection*{Description}

The routine \verb@LALUpdatePulsarCatNode()@ updates all time-varying
properties of the pulsar system to a new epoch specified by the input
\verb@*time@.  The interpretation of this input is specified below,
and can involve the time-varying position of the Earth, as specified
in \verb@*edat@.  Right ascension and declination are udated based on
the specified proper motion, and the pulsar frequency and its
derivatives are updated based on the higher-order derivatives.  For
companion objects, a new periapsis epoch is chosen that is as close as
possible (within half an orbit) of the desired epoch, and all other
time-dependent orbital parameters are updated to this epoch.  All
updates are done ``in place'', to eliminate memory usage and
computation that are, in most cases, unnecessary.

The routine \verb@LALUpdatePulsarCat()@ does the same thing as above,
but to all nodes in the list pointed to by \verb@head@.

The routine \verb@LALDestroyPulsarCat()@ iteratively frees all memory
allocated to the list pointed to by \verb@*head@, and then sets
\verb@*head@ to \verb@NULL@.

\paragraph{Interpretation of \texttt{*time}:} The epoch of a
catalogue update is specified by a \verb@LALPlaceAndGPS@ structure
\verb@*time@, which contains both a GPS time and a detector site.  The
interpretation is as follows: The desired properties are those
properties of the system, as measured by an observer at the solar
system barycentre, that the system had when it emitted the waves that
arrive at the specified detector at the specified GPS time.  By
contrast, the properties listed in the catalogue are those properties
of the system, as measured by an observer at the solar system
barycentre, that the system had when it emitted the waves that arrive
at the solar system barycentre at the GPS time given in the catalogue.
Having specified (in the \verb@*time@ structure) the instant that the
wave fronts reach the detector, \verb@LALUpdatePulsarCatNode()@ first
computes the time when those waves pass the solar system barycentre,
and uses that time as the new catalogue epoch.  Thus, after calling
\verb@LALUpdatePulsarCatNode()@, the GPS times in \verb@*node@ will in
general \emph{not} be the same as the GPS time in \verb@*time@, but
will differ by a light propagation time.

If the \verb@time->p_detector@ field is \verb@NULL@, then
\verb@time->p_gps@ is assumed to be the time when the waves reach the
solar system barycentre, and the complication described above does not
arise.  If the \verb@time->p_gps@ field is \verb@NULL@, then the
routine will return an error.  If \verb@edat@ is \verb@NULL@, then
\verb@time->p@ must also be \verb@NULL@, or an error is returned.

\subsubsection*{Algorithm}

The function \verb@LALUpdatePulsarCatNode()@ first computes the
correct epoch for the pulsar data, taking into account the difference
between the detector time specified in \verb@*time@ and the
barycentric time specified in \verb@*node@: a propagation time delay
id computed using \verb@LALTimeDelayFromEarthCenter()@,
\verb@LALBarycenterEarth()@, and \verb@LALBarycenter()@, with the
pulsar position taken from \verb@*node@ and the Earth ephemeris given
in \verb@*edat@.  This is done in a loop (since updating the epoch can
conceivably change the pulsar location), until the correct epoch is
determined to within $3\mu$s, the percision of \verb@LALBarycenter()@.

Next, the pulsar location is updated using its proper motions, and the
frequency and spindown are updated using the following formula:
$$
\left.\frac{d^k f}{dt^k}\right|_{t=t_2} = \sum_{j=k}^N {j\choose k}
	(t_2-t_1)^{j-k} \left.\frac{d^j f}{dt^j}\right|_{t=t_1} \;.
$$

Uncertainties are propagated using linear propagation methods, under
the conservative presumption that all errors are independent.  Thus
the uncertainty $\sigma_k$ in the frequency component $d^kf/dt^k$ is
given by:
$$
\sigma_k|_{t=t_2} = \sqrt{\sum_{j=k}^N\left[{j\choose k}(t_2-t_1)^{j-k}
	\sigma_j|_{t=t_1}\right]^2 + \left[\sigma_t\sum_{j=k+1}^N
	(j-k){j\choose k}(t_2-t_1)^{j-k-1}
	\left.\frac{d^j f}{dt^j}\right|_{t=t_1}\right]^2} \;,
$$
where $\sigma_t$ is the uncertainty in the time interval, either
3$\mu$s (the precision limit of \verb@LALBarycenter()@) or $10^{-16}$
times the interval (the precision limit of \verb@REAL8@ arithmetic).

An additional uncertainty is assessed based on the assumption that
there may be unmeasured higher-order derivatives of the spin.  The
inverse spin timescale $\tau^{-1}=\max_{k=1}^N\{f^{-1}(d^k
f/dt^k)^{1/k}\}$ is roughly the time that the spin frequency will
change by an amount comparable to its initial value.  Na\"ively, the
next higher frequency derivative will be of order
$d^{N+1}f/dt^{N+1}\sim f(\tau^{-1})^{N+1}$, and will introduce a
further error in $d^k f/dt^k$ equal to: %"
$$
\delta_k \approx f\left[\tau^{-1}(t_2-t_1)\right]^{N+1}
	(t_2-t_1)^{-k} \;.
$$
This uncertainty is added in quadrature to the other errors.

Finally, the properties of any companion orbits are updated using
their first time derivatives, under the implicit assumption that
higher-order derivatives have negligible effects over the interval of
the update.  In the catalogue, all orbital properties are referred to
an epoch of periapsis passage $t_0$.  Keeping only first-order time
derivative corrections to the period, the number of orbits at some
later time $t_0+\Delta t$ is:
$$
n = \frac{\Delta t}{P_0}\left(1-\frac{\dot{P}t}{2P_0}\right) \;,
$$
where $P_0=P(t=t_0)$.  This number is rounded to an integer to get a
periapsis passage near the desired time, and the epoch of this passage
is determined by inverting the formula (and again expanding only to
first order):
$$
t = nP_0\left(1+\mbox{$\frac{1}{2}$}n\dot{P}\right) \;.
$$
Once the new epoch is determined, the period and longitude of
periapsis will be updated using their first derivatives.

\emph{Note:} I am assuming that the orbital period given in the pulsar
catalogue is from periapsis to periapsis.  If it is not, then a
$\dot{w}$ correction will have to be included when computing the
updated epoch.

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel                 LALConvertSkyCoordinates()
LALINT8toGPS()                LALGPStoINT8()
LALBarycenterEarth()          LALBarycenter()
LALDDestroyVector()           LALFree()
\end{verbatim}

\subsubsection*{Notes}

The error propagation has not yet been implemented.

\vfill{\footnotesize\input{PulsarCatCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/LALBarycenter.h>
#include <lal/SkyCoordinates.h>
#include <lal/PulsarCat.h>

NRCSID( PULSARCATC, "$Id$" );


/* First, define a function to compute C(a,b) = (a!)/[(b!)*(a-b)!] */
static UINT4
choose( UINT4 a, UINT4 b );
static UINT4
choose( UINT4 a, UINT4 b )
{
  UINT4 numer = 1;
  UINT4 denom = 1;
  UINT4 index = b + 1;
  while ( --index ) {
    numer *= a - b + index;
    denom *= index;
  }
  return numer/denom;
}


/* <lalVerbatim file="PulsarCatCP"> */
void
LALUpdatePulsarCatNode( LALStatus      *stat,
			PulsarCatNode  *node,
			LALPlaceAndGPS *detectorTime,
			EphemerisData  *edat )
{ /* </lalVerbatim> */
  UINT4 i, j;  /* indecies */
  INT8 t1, t2; /* old and new SSB reference times (ns) */
  REAL8 dt;    /* (new SSB time) - (old SSB time) (s) */
  ConvertSkyParams params; /* sky coordinate parameters */
  CompanionNode *here;     /* paramaters of system companions */

  INITSTATUS( stat, "LALUpdatePulsarCatNode", PULSARCATC );
  ATTATCHSTATUSPTR( stat );

  /* Check that required input parameters exist. */
  ASSERT( node, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );
  ASSERT( detectorTime, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );
  ASSERT( detectorTime->p_gps, stat, PULSARCATH_ENUL,
	  PULSARCATH_MSGENUL );
  if ( !edat ) {
    ASSERT( !(detectorTime->p_detector), stat, PULSARCATH_ENUL,
	    PULSARCATH_MSGENUL );
  }
  if ( node->f ) {
    ASSERT( node->f->data, stat, PULSARCATH_ENUL,
	    PULSARCATH_MSGENUL );
  }

  /* Make sure that the sky positions is in equatorial celestial
     coordinates, since that's what the subroutines require. */
  memset( &params, 0, sizeof(ConvertSkyParams) );
  params.system = COORDINATESYSTEM_EQUATORIAL;
  TRY( LALConvertSkyCoordinates( stat->statusPtr, &(node->pos),
				 &(node->pos), &params ), stat );

  /* Compute the new epoch, if this isn't trivial. */
  TRY( LALGPStoINT8( stat->statusPtr, &t1, &(node->posepoch) ),
       stat );
  TRY( LALGPStoINT8( stat->statusPtr, &t2, detectorTime->p_gps ),
       stat );
  if ( detectorTime->p_detector ) {
    INT8 dt1; /* (new detector time) - (old SSB time) (ns) */
    INT8 dt2; /* (new SSB time) - (new detector time) (ns) */
    EmissionTime emit;     /* output from LALBarycenter() */
    BarycenterInput input; /* input to LALBarycenter() */

    /* Set up barycentring. */
    memset( &input, 0, sizeof(BarycenterInput) );
    input.tgps = *(detectorTime->p_gps);
    input.site = *(detectorTime->p_detector);
    for ( i = 0; i < 3; i++ )
      input.site.location[i] /= LAL_C_SI;
    if ( node->dist > 0.0 )
      input.dInv = 1.0e6*LAL_PC_SI/node->dist;
    else if ( node->dmax > 0.0 )
      input.dInv = 1.0e6*LAL_PC_SI/node->dmax;
    else if ( node->dmin > 0.0 )
      input.dInv = 1.0e6*LAL_PC_SI/node->dmin;
    dt1 = t2 - t1;
    emit.deltaT = 0.0;

    /* Repeatedly apply barycentring and adjust pulsar position until
       we converge to a consistent SSB epoch. */
    do {
      EarthState earth; /* parameters for LALBarycenter() */
      dt2 = (INT8)( 1.0e9*emit.deltaT );
      dt = (1.0e-9)*(REAL8)( dt1 + dt2 );
      input.alpha = node->pos.longitude + node->pmra*dt;
      input.delta = node->pos.latitude + node->pmdec*dt;
      TRY( LALBarycenterEarth( stat->statusPtr, &earth,
			       detectorTime->p_gps, edat ), stat );
      TRY( LALBarycenter( stat->statusPtr, &emit, &input, &earth ),
	   stat );
    } while ( fabs( 1.0e-9*(REAL8)( dt2 ) - emit.deltaT ) > 3.0e-6 );
    t2 += (INT8)( 1.0e9*emit.deltaT );
  }

  /* Adjust pulsar position. */
  dt = (1.0e-9)*(REAL8)( t2 - t1 );
  node->pos.longitude += node->pmra*dt;
  node->pos.latitude += node->pmdec*dt;
  TRY( LALINT8toGPS( stat->statusPtr, &(node->posepoch), &t2 ),
       stat );

  /* Adjust pulsar spin frequency and derivatives. */
  if ( node->f ) {
    TRY( LALGPStoINT8( stat->statusPtr, &t1, &(node->fepoch) ),
	 stat );
    dt = (1.0e-9)*(REAL8)( t2 - t1 );
    for ( i = 0; i < node->f->length; i++ ) {
      REAL8 dtN = 1.0;
      for ( j = i + 1; j < node->f->length; j++ )
	node->f->data[i] += node->f->data[j]*( dtN *= dt )
	  *choose( j, i );
    }
  }
  TRY( LALINT8toGPS( stat->statusPtr, &(node->fepoch), &t2 ),
       stat );

  /* Adjust properties of companion orbits to a nearby epoch: */
  here = node->companion;
  while ( here ) {
    REAL8 n;   /* number of orbital cycles to desired start time */
    INT4 nInt; /* nearest integral number of orbital cycles */

    /* Find the epoch of periapsis passage nearest to the desired
       reference time. */
    TRY( LALGPStoINT8( stat->statusPtr, &t1, &(here->epoch) ), stat );
    dt = (1.0e-9)*(REAL8)( t2 - t1 );
    n = dt/here->p;
    n -= 0.5*(here->pDot)*n*n;
    if ( n > 0.0 )
      nInt = (INT4)( n + 0.5 );
    else
      nInt = (INT4)( n - 0.5 );
    dt = (here->p)*nInt*( 1.0 + 0.5*(here->pDot)*nInt );
    t2 = t1 + (INT8)( 1.0e9*dt );
    TRY( LALINT8toGPS( stat->statusPtr, &(here->epoch), &t2 ), stat );

    /* Update period and longitude of periapsis. */
    here->p += here->pDot*dt;
    here->w += here->wDot*dt;

    /* Proceed to the next companion, if any. */
    here = here->next;
  }

  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="PulsarCatCP"> */
void
LALUpdatePulsarCat( LALStatus      *stat,
		    PulsarCatNode  *head,
		    LALPlaceAndGPS *detectorTime,
		    EphemerisData  *edat )
{ /* </lalVerbatim> */
  INITSTATUS( stat, "LALUpdatePulsarCat", PULSARCATC );
  ATTATCHSTATUSPTR( stat );

  /* Check that head points to a list.  All further pointer checks are
     done in the subroutine. */
  ASSERT( head, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );

  /* Iteratively update each node. */
  while ( head ) {
    TRY( LALUpdatePulsarCatNode( stat->statusPtr, head, detectorTime,
				 edat ), stat );
  }

  /* Done. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}


/* <lalVerbatim file="PulsarCatCP"> */
void
LALDestroyPulsarCat( LALStatus    *stat,
		     PulsarCatNode **head )
{ /* </lalVerbatim> */
  INITSTATUS( stat, "LALDestroyPulsarCat", PULSARCATC );
  ATTATCHSTATUSPTR( stat );

  /* Check that the list exists. */
  ASSERT( head, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );
  ASSERT( *head, stat, PULSARCATH_ENUL, PULSARCATH_MSGENUL );

  /* Free everything iteratively. */
  while ( *head ) {
    PulsarCatNode *here = *head;                /* current node */
    CompanionNode *companion = here->companion; /* companion sublist */
    if ( here->f ) {
      TRY( LALDDestroyVector( stat->statusPtr, &(here->f) ), stat );
    }
    while ( companion ) {
      CompanionNode *here2 = companion; /* current node */
      companion = companion->next;
      LALFree( here2 );
    }
    *head = here->next;
    LALFree( here );
  }

  /* If we got here without sigsegving, we're done. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
