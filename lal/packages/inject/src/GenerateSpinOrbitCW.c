/************************** <lalVerbatim file="GenerateSpinOrbitCWCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\subsection{Module \texttt{GenerateSpinOrbitCW.c}}
\label{ss:GenerateSpinOrbitCW.c}

Computes a spindown- and Doppler-modulated continuous waveform.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GenerateSpinOrbitCWCP}
\idx{LALGenerateSpinOrbitCW()}

\subsubsection*{Description}

This function computes a quaiperiodic waveform using the spindown and
orbital parameters in \verb@*params@, storing the result in
\verb@*output@.

In the \verb@*params@ structure, the routine uses all the ``input''
fields specified in \verb@GenerateSpinOrbitCW.h@, and sets all of the
``output'' fields.  If \verb@params->f@=\verb@NULL@, no spindown
modulation is performed.  If \verb@params->rPeriNorm@=0, no Doppler
modulation is performed.

In the \verb@*output@ structure, the field \verb@output->h@ is
ignored, but all other pointer fields must be set to \verb@NULL@.  The
function will create and allocate space for \verb@output->a@,
\verb@output->f@, and \verb@output->phi@ as necessary.  The
\verb@output->shift@ field will remain set to \verb@NULL@.

\subsubsection*{Algorithm}

This routine calls \verb@LALGenerateCircularSpinOrbitCW()@,
\verb@LALGenerateCircularSpinOrbitCW()@,
\verb@LALGenerateCircularSpinOrbitCW()@, or
\verb@LALGenerateCircularSpinOrbitCW()@, depending on the value of
\verb@params->oneMinusEcc@.  See the other modules under
\verb@GenerateSpinOrbitCW.h@ for descriptions of these routines'
algorithms.

If \verb@params->rPeriNorm@=0, this routine will call
\verb@LALGenerateTaylorCW()@ to generate the waveform.  It creates a
\verb@TaylorCWParamStruc@ from the values in \verb@*params@, adjusting
the values of \verb@params->phi0@, \verb@params->f0@, and
\verb@params->f@ from the reference time \verb@params->spinEpoch@ to
the time \verb@params->epoch@, as follows: Let $\Delta
t=t^{(2)}-t^{(1)}$ be the time between the old epoch $t^{(1)}$ and the
new one $t^{(2)}$.  Then the phase, base frequency, and spindown
parameters for the new epoch are:
\begin{eqnarray}
\phi_0^{(2)} & = & \phi_0^{(1)} + 2\pi f_0^{(1)}t \left( 1 +
	\sum_{k=1}^N \frac{1}{k+1}f_k^{(1)} \Delta t^k \right)
	\nonumber\\
f_0^{(2)} & = & f_0^{(1)} \left( 1 +
	\sum_{k=1}^N f_k^{(1)} \Delta t^k \right)
	\nonumber\\
f_k^{(2)} & = & \frac{f_0^{(1)}}{f_0^{(2)}} \left( f_k^{(1)} +
	\sum_{j=k+1}{N} {j\choose k} f_j^{(1)}\Delta t^{j-k} \right)
	\nonumber
\end{eqnarray}
The phase function $\phi(t)=\phi_0^{(i)}+2\pi
f_0^{(i)}\left[t-t^{(i)}+\sum_{k=1}^N
\frac{f_k^{(i)}}{k+1}\left(t-t^{(i)}\right)^{k+1}\right]$ then has the
same functional dependence on $t$ for either $i=1$ or~2.

\subsubsection*{Uses}
\begin{verbatim}
LALDCreateVector()                      LALDDestroyVector()
LALGenerateCircularSpinOrbitCW()        LALGenerateEllipticSpinOrbitCW()
LALGenerateTaylorCW()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GenerateSpinOrbitCWCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GenerateTaylorCW.h>
#include <lal/GenerateSpinOrbitCW.h>

NRCSID( GENERATESPINORBITCWC, "$Id$" );

/* First, define a function to compute C(a,b) = (a!)/[(b!)*(a-b)!] */
static UINT4
choose( UINT4 a, UINT4 b );
static UINT4
choose( UINT4 a, UINT4 b )
{
  UINT4 numer = 1;
  UINT4 denom = 1;
  UINT4 myindex = b + 1;
  while ( --myindex ) {
    numer *= a - b + myindex;
    denom *= myindex;
  }
  return numer/denom;
}


/* <lalVerbatim file="GenerateSpinOrbitCWCP"> */
void
LALGenerateSpinOrbitCW( LALStatus             *stat,
			CoherentGW            *output,
			SpinOrbitCWParamStruc *params )
{ /* </lalVerbatim> */

  INITSTATUS( stat, "LALGenerateSpinOrbitCW", GENERATESPINORBITCWC );
  ATTATCHSTATUSPTR( stat );

  /* Make sure parameter structure exists (output structure will be
     tested by subroutine). */
  ASSERT( params, stat, GENERATESPINORBITCWH_ENUL,
	  GENERATESPINORBITCWH_MSGENUL );

  /* If there is no orbital motion, use LALGenerateTaylorCW() to
     compute the waveform. */
  if ( params->rPeriNorm == 0.0 ) {
    TaylorCWParamStruc taylorParams; /* subroutine parameters */
    REAL8 t;                         /* time shift */
    memset( &taylorParams, 0, sizeof(TaylorCWParamStruc) );
    taylorParams.position = params->position;
    taylorParams.psi = params->psi;
    taylorParams.epoch = params->epoch;
    taylorParams.deltaT = params->deltaT;
    taylorParams.length = params->length;
    taylorParams.aPlus = params->aPlus;
    taylorParams.aCross = params->aCross;
    taylorParams.phi0 = params->phi0;
    taylorParams.f0 = params->f0;
    t = (REAL8)( params->epoch.gpsSeconds -
		 params->spinEpoch.gpsSeconds );
    t += ( 1.0e-9 )*(REAL8)( params->epoch.gpsNanoSeconds -
			     params->spinEpoch.gpsNanoSeconds );

    /* Adjust epochs. */
    if ( params->f ) {
      UINT4 length = params->f->length; /* number of coefficients */
      UINT4 i, j;       /* indecies over coefficients */
      REAL8 tN = 1.0;   /* t raised to various powers */
      REAL8 fFac = 1.0; /* fractional change in frequency */
      REAL8 tFac = 1.0; /* time integral of fFac */
      REAL8 *f1Data = params->f->data; /* pointer to coeficients */
      REAL8 *f2Data;    /* pointer to corrected coeficients */

      TRY( LALDCreateVector( stat->statusPtr, &(taylorParams.f),
			     length ), stat );
      f2Data = taylorParams.f->data;
      memcpy( f2Data, f1Data, length*sizeof(REAL8) );
      for ( i = 0; i < length; i++ ) {
        REAL8 tM = 1.0; /* t raised to various powers */
        fFac += f1Data[i]*( tN *= t );
        tFac += f1Data[i]*tN/( i + 2.0 );
        for ( j = i + 1; j < length; j++ )
          f2Data[i] += choose( j + 1, i + 1 )*f1Data[j]*( tM *= t );
      }
      taylorParams.phi0 += LAL_TWOPI*taylorParams.f0*t*tFac;
      taylorParams.f0 *= fFac;
      for ( i = 0; i < length; i++ )
        f2Data[i] /= fFac;
    } else
      taylorParams.phi0 += LAL_TWOPI*taylorParams.f0*t;

    /* Generate waveform. */
    LALGenerateTaylorCW( stat->statusPtr, output, &taylorParams );
    BEGINFAIL( stat ) {
      if ( taylorParams.f ) {
	TRY( LALDDestroyVector( stat->statusPtr, &(taylorParams.f) ),
	     stat );
      }
    } ENDFAIL( stat );
    if ( taylorParams.f ) {
      TRY( LALDDestroyVector( stat->statusPtr, &(taylorParams.f) ),
	   stat );
    }
    params->dfdt = taylorParams.dfdt;
  }

  /* If there is orbital motion, call the appropriate subroutine. */

  /* e < 0 is out of range. */
  else if ( params->oneMinusEcc > 1.0 ) {
    ABORT( stat, GENERATESPINORBITCWH_EECC,
	   GENERATESPINORBITCWH_MSGEECC );
  }

  /* 0 <= e < 1 is elliptic. */
  else if ( params->oneMinusEcc > 0.0 ) {
    TRY( LALGenerateEllipticSpinOrbitCW( stat->statusPtr, output,
					 params ), stat );
  }

  /* e = 1 is parabolic. */
  else if ( params->oneMinusEcc == 0.0 ) {
    TRY( LALGenerateParabolicSpinOrbitCW( stat->statusPtr, output,
					  params ), stat );
  }

  /* e > 1 is hyperbolic. */
  else {
    TRY( LALGenerateHyperbolicSpinOrbitCW( stat->statusPtr, output,
					   params ), stat );
  }

  /* That is all. */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
