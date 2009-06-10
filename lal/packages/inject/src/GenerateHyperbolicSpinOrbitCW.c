/*
*  Copyright (C) 2007 Teviet Creighton
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

/**************** <lalVerbatim file="GenerateHyperbolicSpinOrbitCWCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{GenerateHyperbolicSpinOrbitCW.c}}
\label{ss:GenerateHyperbolicSpinOrbitCW.c}

Computes a continuous waveform with frequency drift and Doppler
modulation from a hyperbolic orbital trajectory.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GenerateHyperbolicSpinOrbitCWCP}
\idx{LALGenerateHyperbolicSpinOrbitCW()}

\subsubsection*{Description}

This function computes a quaiperiodic waveform using the spindown and
orbital parameters in \verb@*params@, storing the result in
\verb@*output@.

In the \verb@*params@ structure, the routine uses all the ``input''
fields specified in \verb@GenerateSpinOrbitCW.h@, and sets all of the
``output'' fields.  If \verb@params->f@=\verb@NULL@, no spindown
modulation is performed.  If \verb@params->oneMinusEcc@$\not<0$ (a
non-hyperbolic orbit), or if
\verb@params->rPeriNorm@$\times$\verb@params->angularSpeed@$\geq1$
(faster-than-light speed at periapsis), an error is returned.

In the \verb@*output@ structure, the field \verb@output->h@ is
ignored, but all other pointer fields must be set to \verb@NULL@.  The
function will create and allocate space for \verb@output->a@,
\verb@output->f@, and \verb@output->phi@ as necessary.  The
\verb@output->shift@ field will remain set to \verb@NULL@.

\subsubsection*{Algorithm}

For hyperbolic orbits, we combine Eqs.~(\ref{eq:spinorbit-tr}),
(\ref{eq:spinorbit-t}), and~(\ref{eq:spinorbit-upsilon}) to get $t_r$
directly as a function of $E$:
\begin{eqnarray}
t_r = t_p & + & \left(\frac{r_p \sin i}{c}\right)\sin\omega \nonumber\\
\label{eq:tr-e3}
	& + & \frac{1}{n} \left( -E +
		\left[v_p(e-1)\cos\omega + e\right]\sinh E
		- \left[v_p\sqrt{\frac{e-1}{e+1}}\sin\omega\right]
			[\cosh E - 1]\right) \;,
\end{eqnarray}
where $v_p=r_p\dot{\upsilon}_p\sin i/c$ is a normalized velocity at
periapsis and $n=\dot{\upsilon}_p\sqrt{(1-e)^3/(1+e)}$ is a normalized
angular speed for the orbit (the hyperbolic analogue of the mean
angular speed for closed orbits).  For simplicity we write this as:
\begin{equation}
\label{eq:tr-e4}
t_r = T_p + \frac{1}{n}\left( E + A\sinh E + B[\cosh E - 1] \right) \;,
\end{equation}
\begin{wrapfigure}{r}{0.28\textwidth}
\vspace{-5ex}
\begin{center}
\resizebox{0.23\textwidth}{!}{\includegraphics{inject_hanomaly}} \\
%\parbox{0.23\textwidth}{\caption{\label{fig:binary-orbit} Function to
%be inverted to find eccentric anomaly.}}
\end{center}
\vspace{-5ex}
\end{wrapfigure}
where $T_p$ is the \emph{observed} time of periapsis passage.  Thus
the key numerical procedure in this routine is to invert the
expression $x=E+A\sinh E+B(\cosh E - 1)$ to get $E(x)$.  This function
is sketched to the right (solid line), along with an approximation
used for making an initial guess (dotted line), as described later.

We note that $A^2-B^2<1$, although it approaches 1 when
$e\rightarrow1$, or when $v_p\rightarrow1$ and either $e=0$ or
$\omega=\pi$.  Except in this limit, Newton-Raphson methods will
converge rapidly for any initial guess.  In this limit, though, the
slope $dx/dE$ approaches zero at $E=0$, and an initial guess or
iteration landing near this point will send the next iteration off to
unacceptably large or small values.  A hybrid root-finding strategy is
used to deal with this, and with the exponential behaviour of $x$ at
large $E$.

First, we compute $x=x_{\pm1}$ at $E=\pm1$.  If the desired $x$ lies
in this range, we use a straightforward Newton-Raphson root finder,
with the constraint that all guesses of $E$ are restricted to the
domain $[-1,1]$.  This guarantees that the scheme will eventually find
itself on a uniformly-convergent trajectory.

Second, for $E$ outside of this range, $x$ is dominated by the
exponential terms: $x\approx\frac{1}{2}(A+B)\exp(E)$ for $E\gg1$, and
$x\approx-\frac{1}{2}(A-B)\exp(-E)$ for $E\ll-1$.  We therefore do an
\emph{approximate} Newton-Raphson iteration on the function $\ln|x|$,
where the approximation is that we take $d\ln|x|/d|E|\approx1$.  This
involves computing an extra logarithm inside the loop, but gives very
rapid convergence to high precision, since $\ln|x|$ is very nearly
linear in these regions.

At the start of the algorithm, we use an initial guess of
$E=-\ln[-2(x-x_{-1})/(A-B)-\exp(1)]$ for $x<x_{-1}$, $E=x/x_{-1}$ for
$x_{-1}\leq x\leq0$, $E=x/x_{+1}$ for $0\leq x\leq x_{+1}$, or
$E=\ln[2(x-x_{+1})/(A+B)-\exp(1)]$ for $x>x_{+1}$.  We refine this
guess until we get agreement to within 0.01 parts in part in
$N_\mathrm{cyc}$ (where $N_\mathrm{cyc}$ is the larger of the number
of wave cycles in a time $2\pi/n$, or the number of wave cycles in the
entire waveform being generated), or one part in $10^{15}$ (an order
of magnitude off the best precision possible with \verb@REAL8@
numbers).  The latter case indicates that \verb@REAL8@ precision may
fail to give accurate phasing, and one should consider modeling the
orbit as a set of Taylor frequency coefficients \'{a} la
\verb@LALGenerateTaylorCW()@.  On subsequent timesteps, we use the
previous timestep as an initial guess, which is good so long as the
timesteps are much smaller than $1/n$.

Once a value of $E$ is found for a given timestep in the output
series, we compute the system time $t$ via Eq.~(\ref{eq:spinorbit-t}),
and use it to determine the wave phase and (non-Doppler-shifted)
frequency via Eqs.~(\ref{eq:taylorcw-freq})
and~(\ref{eq:taylorcw-phi}).  The Doppler shift on the frequency is
then computed using Eqs.~(\ref{eq:spinorbit-upsilon})
and~(\ref{eq:orbit-rdot}).  We use $\upsilon$ as an intermediate in
the Doppler shift calculations, since expressing $\dot{R}$ directly in
terms of $E$ results in expression of the form $(e-1)/(e\cosh E-1)$,
which are difficult to simplify and face precision losses when
$E\sim0$ and $e\rightarrow1$.  By contrast, solving for $\upsilon$ is
numerically stable provided that the system \verb@atan2()@ function is
well-designed.

This routine does not account for relativistic timing variations, and
issues warnings or errors based on the criterea of
Eq.~(\ref{eq:relativistic-orbit}) in
\verb@GenerateEllipticSpinOrbitCW.c@.

\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                   LALFree()
LALSCreateVectorSequence()    LALSDestroyVectorSequence()
LALSCreateVector()            LALSDestroyVector()
LALDCreateVector()            LALDDestroyVector()
LALSnprintf()                 LALWarning()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GenerateHyperbolicSpinOrbitCWCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateSpinOrbitCW.h>

NRCSID( GENERATEHYPERBOLICSPINORBITCWC, "$Id$" );

/* <lalVerbatim file="GenerateHyperbolicSpinOrbitCWCP"> */
void
LALGenerateHyperbolicSpinOrbitCW( LALStatus             *stat,
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
  REAL8 vDotAvg;           /* nomalized orbital angular speed */
  REAL8 vp;                /* projected speed at periapsis */
  REAL8 upsilon, argument; /* true anomaly, and argument of periapsis */
  REAL8 eCosOmega;         /* eccentricity * cosine of argument */
  REAL8 tPeriObs;          /* time of observed periapsis */
  REAL8 spinOff;           /* spin epoch - orbit epoch */
  REAL8 x;                 /* observed ``mean anomaly'' */
  REAL8 xPlus, xMinus;     /* limits where exponentials dominate */
  REAL8 dx, dxMax;         /* current and target errors in x */
  REAL8 a, b;              /* constants in equation for x */
  REAL8 ecc;               /* orbital eccentricity */
  REAL8 eccMinusOne, eccPlusOne; /* ecc - 1 and ecc + 1 */
  REAL8 e;                       /* hyperbolic anomaly */
  REAL8 sinhe, coshe;            /* sinh of e, and cosh of e minus 1 */
  REAL8 *fSpin = NULL;           /* pointer to Taylor coefficients */
  REAL4 *fData;                  /* pointer to frequency data */
  REAL8 *phiData;                /* pointer to phase data */

  INITSTATUS( stat, "LALGenerateHyperbolicSpinOrbitCW",
	      GENERATEHYPERBOLICSPINORBITCWC );
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
  eccMinusOne = -params->oneMinusEcc;
  ecc = 1.0 + eccMinusOne;
  eccPlusOne = 2.0 + eccMinusOne;
  if ( eccMinusOne <= 0.0 ) {
    ABORT( stat, GENERATESPINORBITCWH_EECC,
	   GENERATESPINORBITCWH_MSGEECC );
  }
  vp = params->rPeriNorm*params->angularSpeed;
  vDotAvg = params->angularSpeed
    *sqrt( eccMinusOne*eccMinusOne*eccMinusOne/eccPlusOne );
  n = params->length;
  dt = params->deltaT;
  f0 = fPrev = params->f0;
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
  a = vp*eccMinusOne*cos( argument ) + ecc;
  b = -vp*sqrt( eccMinusOne/eccPlusOne )*sin( argument );
  eCosOmega = ecc*cos( argument );
  if ( n*dt*vDotAvg > LAL_TWOPI )
    dxMax = 0.01/( f0*n*dt );
  else
    dxMax = 0.01/( f0*LAL_TWOPI/vDotAvg );
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
    if ( tau > LAL_TWOPI/vDotAvg )
      tau = LAL_TWOPI/vDotAvg;
    if ( f0*tau*vp*vp*ecc/eccPlusOne > 0.25 )
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

  /* Determine bounds of hybrid root-finding algorithm, and initial
     guess for e. */
  xMinus = 1.0 + a*sinh( -1.0 ) + b*cosh( -1.0 ) - b;
  xPlus = -1.0 + a*sinh( 1.0 ) + b*cosh( 1.0 ) - b;
  x = -vDotAvg*tPeriObs;
  if ( x < xMinus )
    e = -log( -2.0*( x - xMinus )/( a - b ) - exp( 1.0 ) );
  else if ( x <= 0 )
    e = x/xMinus;
  else if ( x <= xPlus )
    e = x/xPlus;
  else
    e = log( 2.0*( x - xPlus )/( a + b ) - exp( 1.0 ) );
  sinhe = sinh( e );
  coshe = cosh( e ) - 1.0;

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
  LALSnprintf( output->a->name, LALNameLength, "CW amplitudes" );
  LALSnprintf( output->f->name, LALNameLength, "CW frequency" );
  LALSnprintf( output->phi->name, LALNameLength, "CW phase" );

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

    x = vDotAvg*( i*dt - tPeriObs );

    /* Use approximate Newton-Raphson method on ln|x| if |x| > 1. */
    if ( x < xMinus ) {
      x = log( -x );
      while ( fabs( dx = log( e - a*sinhe - b*coshe ) - x ) > dxMax ) {
	e += dx;
	sinhe = sinh( e );
	coshe = cosh( e ) - 1.0;
      }
    }
    else if ( x > xPlus ) {
      x = log( x );
      while ( fabs( dx = log( -e + a*sinhe + b*coshe ) - x ) > dxMax ) {
	e -= dx;
	sinhe = sinh( e );
	coshe = cosh( e ) - 1.0;
      }
    }

    /* Use ordinary Newton-Raphson method on x if |x| <= 1. */
    else {
      while ( fabs( dx = -e + a*sinhe + b*coshe - x ) > dxMax ) {
	e -= dx/( -1.0 + a*coshe + a + b*sinhe );
	if ( e < -1.0 )
	  e = -1.0;
	else if ( e > 1.0 )
	  e = 1.0;
	sinhe = sinh( e );
	coshe = cosh( e ) - 1.0;
      }
    }

    /* Compute source emission time, phase, and frequency. */
    phi = t = tPow =
      ( ecc*sinhe - e )/vDotAvg + spinOff;
    f = 1.0;
    for ( j = 0; j < nSpin; j++ ) {
      f += fSpin[j]*tPow;
      phi += fSpin[j]*( tPow*=t )/( j + 2.0 );
    }

    /* Appy frequency Doppler shift. */
    upsilon = 2.0*atan2( sqrt( eccPlusOne*coshe ),
			 sqrt( eccMinusOne*( coshe + 2.0 ) ) );
    f *= f0 / ( 1.0 + vp*( cos( argument + upsilon ) + eCosOmega )
		/eccPlusOne );
    phi *= twopif0;
    if ( fabs( f - fPrev ) > df )
      df = fabs( f - fPrev );
    *(fData++) = fPrev = f;
    *(phiData++) = phi + phi0;
  }

  /* Set output field and return. */
  params->dfdt = df*dt;
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
