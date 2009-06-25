/*
*  Copyright (C) 2007 Reinhard Prix, Teviet Creighton
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

/****************** <lalVerbatim file="GenerateEllipticSpinOrbitCWCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}
\providecommand{\greatersim}{\stackrel{>}{\scriptstyle\sim}}

\subsection{Module \texttt{GenerateEllipticSpinOrbitCW.c}}
\label{ss:GenerateEllipticSpinOrbitCW.c}

Computes a continuous waveform with frequency drift and Doppler
modulation from an elliptical orbital trajectory.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GenerateEllipticSpinOrbitCWCP}
\idx{LALGenerateEllipticSpinOrbitCW()}

\subsubsection*{Description}

This function computes a quaiperiodic waveform using the spindown and
orbital parameters in \verb@*params@, storing the result in
\verb@*output@.

In the \verb@*params@ structure, the routine uses all the ``input''
fields specified in \verb@GenerateSpinOrbitCW.h@, and sets all of the
``output'' fields.  If \verb@params->f@=\verb@NULL@, no spindown
modulation is performed.  If \verb@params->oneMinusEcc@$\notin(0,1]$
(an open orbit), or if
\verb@params->rPeriNorm@$\times$\verb@params->angularSpeed@$\geq1$
(faster-than-light speed at periapsis), an error is returned.

In the \verb@*output@ structure, the field \verb@output->h@ is
ignored, but all other pointer fields must be set to \verb@NULL@.  The
function will create and allocate space for \verb@output->a@,
\verb@output->f@, and \verb@output->phi@ as necessary.  The
\verb@output->shift@ field will remain set to \verb@NULL@.

\subsubsection*{Algorithm}

For elliptical orbits, we combine Eqs.~(\ref{eq:spinorbit-tr}),
(\ref{eq:spinorbit-t}), and~(\ref{eq:spinorbit-upsilon}) to get $t_r$
directly as a function of the eccentric anomaly $E$:
\begin{eqnarray}
t_r = t_p & + & \left(\frac{r_p \sin i}{c}\right)\sin\omega \nonumber\\
\label{eq:tr-e1}
	& + & \left(\frac{P}{2\pi}\right) \left( E +
		\left[v_p(1-e)\cos\omega - e\right]\sin E
		+ \left[v_p\sqrt{\frac{1-e}{1+e}}\sin\omega\right]
			[\cos E - 1]\right) \;,
\end{eqnarray}
where $v_p=r_p\dot{\upsilon}_p\sin i/c$ is a normalized velocity at
periapsis and $P=2\pi\sqrt{(1+e)/(1-e)^3}/\dot{\upsilon}_p$ is the
period of the orbit.  For simplicity we write this as:
\begin{equation}
\label{eq:tr-e2}
t_r = T_p + \frac{1}{n}\left( E + A\sin E + B[\cos E - 1] \right) \;,
\end{equation}
\begin{wrapfigure}{r}{0.28\textwidth}
\vspace{-4ex}
\begin{center}
\resizebox{0.23\textwidth}{!}{\includegraphics{inject_eanomaly}} \\
%\parbox{0.23\textwidth}{\caption{\label{fig:binary-orbit} Function to
%be inverted to find eccentric anomaly.}}
\end{center}
\vspace{-4ex}
\end{wrapfigure}
where $T_p$ is the \emph{observed} time of periapsis passage and
$n=2\pi/P$ is the mean angular speed around the orbit.  Thus the key
numerical procedure in this routine is to invert the expression
$x=E+A\sin E+B(\cos E - 1)$ to get $E(x)$.  We note that
$E(x+2n\pi)=E(x)+2n\pi$, so we only need to solve this expression in
the interval $[0,2\pi)$, sketched to the right.

We further note that $A^2+B^2<1$, although it approaches 1 when
$e\rightarrow1$, or when $v_p\rightarrow1$ and either $e=0$ or
$\omega=\pi$.  Except in this limit, Newton-Raphson methods will
converge rapidly for any initial guess.  In this limit, though, the
slope $dx/dE$ approaches zero at the point of inflection, and an
initial guess or iteration landing near this point will send the next
iteration off to unacceptably large or small values.  However, by
restricting all initial guesses and iterations to the domain
$E\in[0,2\pi)$, one will always end up on a trajectory branch that
will converge uniformly.  This should converge faster than the more
generically robust technique of bisection.

In this algorithm, we start the computation with an arbitrary initial
guess of $E=0$, and refine it until the we get agreement to within
0.01 parts in part in $N_\mathrm{cyc}$ (where $N_\mathrm{cyc}$ is the
larger of the number of wave cycles in an orbital period, or the
number of wave cycles in the entire waveform being generated), or one
part in $10^{15}$ (an order of magnitude off the best precision
possible with \verb@REAL8@ numbers).  The latter case indicates that
\verb@REAL8@ precision may fail to give accurate phasing, and one
should consider modeling the orbit as a set of Taylor frequency
coefficients \'{a} la \verb@LALGenerateTaylorCW()@.  On subsequent
timesteps, we use the previous timestep as an initial guess, which is
good so long as the timesteps are much smaller than an orbital period.
This sequence of guesses will have to readjust itself once every orbit
(as $E$ jumps from $2\pi$ down to 0), but this is relatively
infrequent; we don't bother trying to smooth this out because the
additional tests would probably slow down the algorithm overall.

Once a value of $E$ is found for a given timestep in the output
series, we compute the system time $t$ via Eq.~(\ref{eq:spinorbit-t}),
and use it to determine the wave phase and (non-Doppler-shifted)
frequency via Eqs.~(\ref{eq:taylorcw-freq})
and~(\ref{eq:taylorcw-phi}).  The Doppler shift on the frequency is
then computed using Eqs.~(\ref{eq:spinorbit-upsilon})
and~(\ref{eq:orbit-rdot}).  We use $\upsilon$ as an intermediate in
the Doppler shift calculations, since expressing $\dot{R}$ directly in
terms of $E$ results in expression of the form $(1-e)/(1-e\cos E)$,
which are difficult to simplify and face precision losses when
$E\sim0$ and $e\rightarrow1$.  By contrast, solving for $\upsilon$ is
numerically stable provided that the system \verb@atan2()@ function is
well-designed.

The routine does not account for variations in special relativistic or
gravitational time dilation due to the elliptical orbit, nor does it
deal with other gravitational effects such as Shapiro delay.  To a
very rough approximation, the amount of phase error induced by
gravitational redshift goes something like $\Delta\phi\sim
fT(v/c)^2\Delta(r_p/r)$, where $f$ is the typical wave frequency, $T$
is either the length of data or the orbital period (whichever is
\emph{smaller}), $v$ is the \emph{true} (unprojected) speed at
periapsis, and $\Delta(r_p/r)$ is the total range swept out by the
quantity $r_p/r$ over the course of the observation.  Other
relativistic effects such as special relativistic time dilation are
comparable in magnitude.  We make a crude estimate of when this is
significant by noting that $v/c\greatersim v_p$ but
$\Delta(r_p/r)\lessim 2e/(1+e)$; we take these approximations as
equalities and require that $\Delta\phi\lessim\pi$, giving:
\begin{equation}
\label{eq:relativistic-orbit}
f_0Tv_p^2\frac{4e}{1+e}\lessim1 \;.
\end{equation}
When this critereon is violated, a warning is generated.  Furthermore,
as noted earlier, when $v_p\geq1$ the routine will return an error, as
faster-than-light speeds can cause the emission and reception times to
be non-monotonic functions of one another.

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                   LALFree()
LALSCreateVectorSequence()    LALSDestroyVectorSequence()
LALSCreateVector()            LALSDestroyVector()
LALDCreateVector()            LALDDestroyVector()
snprintf()                 LALWarning()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GenerateEllipticSpinOrbitCWCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateSpinOrbitCW.h>

NRCSID( GENERATEELLIPTICSPINORBITCWC, "$Id$" );

/* <lalVerbatim file="GenerateEllipticSpinOrbitCWCP"> */
void
LALGenerateEllipticSpinOrbitCW( LALStatus             *stat,
				CoherentGW            *output,
				SpinOrbitCWParamStruc *params )
{ /* </lalVerbatim> */
  UINT4 n, i;              /* number of and index over samples */
  UINT4 nSpin = 0, j;      /* number of and index over spindown terms */
  REAL8 t, dt, tPow;       /* time, interval, and t raised to a power */
  REAL8 phi0, f0, twopif0; /* initial phase, frequency, and 2*pi*f0 */
  REAL8 f, fPrev;          /* current and previous values of frequency */
  REAL4 df = 0.0;          /* maximum difference between f and fPrev */
  REAL8 phi;               /* current value of phase */
  REAL8 p, vDotAvg;        /* orbital period, and 2*pi/(period) */
  REAL8 vp;                /* projected speed at periapsis */
  REAL8 upsilon, argument; /* true anomaly, and argument of periapsis */
  REAL8 eCosOmega;         /* eccentricity * cosine of argument */
  REAL8 tPeriObs;          /* time of observed periapsis */
  REAL8 spinOff;           /* spin epoch - orbit epoch */
  REAL8 x;                 /* observed mean anomaly */
  REAL8 dx, dxMax;         /* current and target errors in x */
  REAL8 a, b;              /* constants in equation for x */
  REAL8 ecc;               /* orbital eccentricity */
  REAL8 oneMinusEcc, onePlusEcc; /* 1 - ecc and 1 + ecc */
  REAL8 e = 0.0;                 /* eccentric anomaly */
  REAL8 sine = 0.0, cose = 0.0;  /* sine of e, and cosine of e minus 1 */
  REAL8 *fSpin = NULL;           /* pointer to Taylor coefficients */
  REAL4 *fData;                  /* pointer to frequency data */
  REAL8 *phiData;                /* pointer to phase data */

  INITSTATUS( stat, "LALGenerateEllipticSpinOrbitCW",
	      GENERATEELLIPTICSPINORBITCWC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter and output structures exist. */
  ASSERT( params, stat, GENERATESPINORBITCWH_ENUL,
	  GENERATESPINORBITCWH_MSGENUL );
  ASSERT( output, stat, GENERATESPINORBITCWH_ENUL,
	  GENERATESPINORBITCWH_MSGENUL );

  /* Make sure output fields don't exist. */
  ASSERT( !( output->a ), stat, GENERATESPINORBITCWH_EOUT,
	  GENERATESPINORBITCWH_MSGEOUT );
  ASSERT( !( output->f ), stat, GENERATESPINORBITCWH_EOUT,
	  GENERATESPINORBITCWH_MSGEOUT );
  ASSERT( !( output->phi ), stat, GENERATESPINORBITCWH_EOUT,
	  GENERATESPINORBITCWH_MSGEOUT );
  ASSERT( !( output->shift ), stat, GENERATESPINORBITCWH_EOUT,
	  GENERATESPINORBITCWH_MSGEOUT );

  /* If Taylor coeficients are specified, make sure they exist. */
  if ( params->f ) {
    ASSERT( params->f->data, stat, GENERATESPINORBITCWH_ENUL,
	    GENERATESPINORBITCWH_MSGENUL );
    nSpin = params->f->length;
    fSpin = params->f->data;
  }

  /* Set up some constants (to avoid repeated calculation or
     dereferencing), and make sure they have acceptable values. */
  oneMinusEcc = params->oneMinusEcc;
  ecc = 1.0 - oneMinusEcc;
  onePlusEcc = 1.0 + ecc;
  if ( ecc < 0.0 || oneMinusEcc <= 0.0 ) {
    ABORT( stat, GENERATESPINORBITCWH_EECC,
	   GENERATESPINORBITCWH_MSGEECC );
  }
  vp = params->rPeriNorm*params->angularSpeed;
  n = params->length;
  dt = params->deltaT;
  f0 = fPrev = params->f0;
  vDotAvg = params->angularSpeed
    *sqrt( oneMinusEcc*oneMinusEcc*oneMinusEcc/onePlusEcc );
  if ( vp >= 1.0 ) {
    ABORT( stat, GENERATESPINORBITCWH_EFTL,
	   GENERATESPINORBITCWH_MSGEFTL );
  }
  if ( vp <= 0.0 || dt <= 0.0 || f0 <= 0.0 || vDotAvg <= 0.0 ||
       n == 0 ) {
    ABORT( stat, GENERATESPINORBITCWH_ESGN,
	   GENERATESPINORBITCWH_MSGESGN );
  }

  /* Set up some other constants. */
  twopif0 = f0*LAL_TWOPI;
  phi0 = params->phi0;
  argument = params->omega;
  p = LAL_TWOPI/vDotAvg;
  a = vp*oneMinusEcc*cos( argument ) + oneMinusEcc - 1.0;
  b = vp*sqrt( oneMinusEcc/( onePlusEcc ) )*sin( argument );
  eCosOmega = ecc*cos( argument );
  if ( n*dt > p )
    dxMax = 0.01/( f0*n*dt );
  else
    dxMax = 0.01/( f0*p );
  if ( dxMax < 1.0e-15 ) {
    dxMax = 1.0e-15;
#ifndef NDEBUG
    LALWarning( stat, "REAL8 arithmetic may not have sufficient"
		" precision for this orbit" );
#endif
  }
#ifndef NDEBUG
  if ( lalDebugLevel & LALWARNING ) {
    REAL8 tau = n*dt;
    if ( tau > p )
      tau = p;
    if ( f0*tau*vp*vp*ecc/onePlusEcc > 0.25 )
      LALWarning( stat, "Orbit may have significant relativistic"
		  " effects that are not included" );
  }
#endif

  /* Compute offset between time series epoch and observed periapsis,
     and betweem true periapsis and spindown reference epoch. */
  tPeriObs = (REAL8)( params->orbitEpoch.gpsSeconds -
		      params->epoch.gpsSeconds );
  tPeriObs += 1.0e-9 * (REAL8)( params->orbitEpoch.gpsNanoSeconds -
				params->epoch.gpsNanoSeconds );
  tPeriObs += params->rPeriNorm*sin( params->omega );
  spinOff = (REAL8)( params->orbitEpoch.gpsSeconds -
		     params->spinEpoch.gpsSeconds );
  spinOff += 1.0e-9 * (REAL8)( params->orbitEpoch.gpsNanoSeconds -
			       params->spinEpoch.gpsNanoSeconds );

  /* Allocate output structures. */
  if ( ( output->a = (REAL4TimeVectorSeries *)
	 LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL ) {
    ABORT( stat, GENERATESPINORBITCWH_EMEM,
	   GENERATESPINORBITCWH_MSGEMEM );
  }
  memset( output->a, 0, sizeof(REAL4TimeVectorSeries) );
  if ( ( output->f = (REAL4TimeSeries *)
	 LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    ABORT( stat, GENERATESPINORBITCWH_EMEM,
	   GENERATESPINORBITCWH_MSGEMEM );
  }
  memset( output->f, 0, sizeof(REAL4TimeSeries) );
  if ( ( output->phi = (REAL8TimeSeries *)
	 LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL ) {
    LALFree( output->a ); output->a = NULL;
    LALFree( output->f ); output->f = NULL;
    ABORT( stat, GENERATESPINORBITCWH_EMEM,
	   GENERATESPINORBITCWH_MSGEMEM );
  }
  memset( output->phi, 0, sizeof(REAL8TimeSeries) );

  /* Set output structure metadata fields. */
  output->position = params->position;
  output->psi = params->psi;
  output->a->epoch = output->f->epoch = output->phi->epoch
    = params->epoch;
  output->a->deltaT = n*params->deltaT;
  output->f->deltaT = output->phi->deltaT = params->deltaT;
  output->a->sampleUnits = lalStrainUnit;
  output->f->sampleUnits = lalHertzUnit;
  output->phi->sampleUnits = lalDimensionlessUnit;
  snprintf( output->a->name, LALNameLength, "CW amplitudes" );
  snprintf( output->f->name, LALNameLength, "CW frequency" );
  snprintf( output->phi->name, LALNameLength, "CW phase" );

  /* Allocate phase and frequency arrays. */
  LALSCreateVector( stat->statusPtr, &( output->f->data ), n );
  BEGINFAIL( stat ) {
    LALFree( output->a );   output->a = NULL;
    LALFree( output->f );   output->f = NULL;
    LALFree( output->phi ); output->phi = NULL;
  } ENDFAIL( stat );
  LALDCreateVector( stat->statusPtr, &( output->phi->data ), n );
  BEGINFAIL( stat ) {
    TRY( LALSDestroyVector( stat->statusPtr, &( output->f->data ) ),
	 stat );
    LALFree( output->a );   output->a = NULL;
    LALFree( output->f );   output->f = NULL;
    LALFree( output->phi ); output->phi = NULL;
  } ENDFAIL( stat );

  /* Allocate and fill amplitude array. */
  {
    CreateVectorSequenceIn in; /* input to create output->a */
    in.length = 2;
    in.vectorLength = 2;
    LALSCreateVectorSequence( stat->statusPtr, &(output->a->data), &in );
    BEGINFAIL( stat ) {
      TRY( LALSDestroyVector( stat->statusPtr, &( output->f->data ) ),
	   stat );
      TRY( LALDDestroyVector( stat->statusPtr, &( output->phi->data ) ),
	   stat );
      LALFree( output->a );   output->a = NULL;
      LALFree( output->f );   output->f = NULL;
      LALFree( output->phi ); output->phi = NULL;
    } ENDFAIL( stat );
    output->a->data->data[0] = output->a->data->data[2] = params->aPlus;
    output->a->data->data[1] = output->a->data->data[3] = params->aCross;
  }

  /* Fill frequency and phase arrays. */
  fData = output->f->data->data;
  phiData = output->phi->data->data;
  for ( i = 0; i < n; i++ ) {
    INT4 nOrb; /* number of orbits since the specified orbit epoch */

    /* First, find x in the range [0,2*pi]. */
    x = vDotAvg*( i*dt - tPeriObs );
    nOrb = (INT4)( x/LAL_TWOPI );
    if ( x < 0.0 )
      nOrb -= 1;
    x -= LAL_TWOPI*nOrb;

    /* Newton-Raphson iteration to find E(x). */
    while ( fabs( dx = e + a*sine + b*cose - x ) > dxMax ) {
      e -= dx/( 1.0 + a*cose + a - b*sine );
      if ( e < 0.0 )
	e = 0.0;
      else if ( e > LAL_TWOPI )
	e = LAL_TWOPI;
      sine = sin( e );
      cose = cos( e ) - 1.0;
    }

    /* Compute source emission time, phase, and frequency. */
    phi = t = tPow =
      ( e + LAL_TWOPI*nOrb - ecc*sine )/vDotAvg + spinOff;
    f = 1.0;
    for ( j = 0; j < nSpin; j++ ) {
      f += fSpin[j]*tPow;
      phi += fSpin[j]*( tPow*=t )/( j + 2.0 );
    }

    /* Appy frequency Doppler shift. */
    upsilon = 2.0*atan2( sqrt( -onePlusEcc*cose ), sqrt( oneMinusEcc*( cose + 2.0 ) ) );
    f *= f0 / ( 1.0 + vp*( cos( argument + upsilon ) + eCosOmega ) /onePlusEcc );
    phi *= twopif0;
    if ( (i > 0) && (fabs( f - fPrev ) > df) )
      df = fabs( f - fPrev );
    *(fData++) = fPrev = f;
    *(phiData++) = phi + phi0;

  } /* for i < n */

  /* Set output field and return. */
  params->dfdt = df*dt;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
