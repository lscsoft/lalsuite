/**************************************** <lalVerbatim file="PtoleMetricCV">
Author: Owen, B. J.
$Id$
********************************************************** </lalVerbatim> */

/**************************************************************** <lalLaTeX>

\subsection{Module \texttt{LALPtoleMetric}}
\label{ss:PtoleMetric}

Computes metric components for a pulsar search in the ``Ptolemaic''
approximation.

\subsubsection*{Prototypes}
\input{PtoleMetricCP}
\idx{PtoleMetric}

\subsubsection*{Description}

This function computes metric components in a way that yields results very
similar to those of \texttt{CoherentMetric} called with the output of
\texttt{TBaryPtolemaic}. The CPU demand, however, is less, and the metric
components can be expressed analytically, lending themselves to better
understanding of the behavior of the parameter space.

\subsubsection*{Algorithm}

For speed and checking reasons, a minimum of numerical computation is
involved. The metric components can be expressed analytically (though not
tidily) in terms of trig functions.

With a subscript 0 referring to the $f_0$ parameter and subscript $j$
referring to the spindown parameter $f_j$, the metric components (without
orbital motion) are given by
\begin{eqnarray}
g_{00} &=& \pi^2T^2/3,
\\
g_{0\alpha} &=& 2\pi^2 {f_0r \over \omega} \cos\delta \left\{ -\cos\psi_1
-\cos\psi_0 + {2\over\omega T} \left[ \sin\psi_1 - \sin\psi_0 \right]
\right\},
\\
g_{0\delta} &=& -2\pi^2 {f_0r \over \omega} \sin\delta \left\{ \sin\psi_1 +
\sin\psi_0 + {2\over \omega T} \left[ \cos\psi_1 - \cos\psi_0 \right]
\right\},
\\
g_{0j} &=& {2\pi^2f_0T^2 \over (j+2)(j+3)},
\\
g_{\alpha\alpha} &=& 2(\pi f_0r\cos\delta)^2 \left\{1 - {1\over\omega T}
\left[ \sin\psi_1 \cos\psi_1 - \sin\psi_0 \cos\psi_0 \right] -
{2\over(\omega T)^2} \left[ \cos\psi_1 - \cos\psi_0 \right]^2 \right\},
\\
g_{\alpha\delta} &=& -2(\pi f_0r)^2 \sin\delta \cos\delta \left\{
{1\over\omega T} \left[ \sin^2\psi_1 - \sin^2\psi_0 \right] + {2\over\omega
T}^2 \left[ \sin\psi_1 - \sin\psi_0 \right] \left[ \cos\psi_1 - \cos\psi_0
\right] \right\},
\\
g_{\delta\delta} &=& 2(\pi f_0r\sin\delta)^2 \left\{1 + {1\over\omega T}
\left[ \sin\psi_1 \cos\psi_1 - \sin\psi_0 \cos\psi_0 \right] -
{2\over(\omega T)^2} \left[ \sin\psi_1 - \sin\psi_0 \right]^2 \right\},
\\
g_{\alpha j} &=& -{4\pi^2 \over j+1} {f_0^2r\cos\delta \over \omega} \left[
(-)^{j//2} {(j+1)! \over (\omega T)^{j+1}} \left\{
\begin{array}{ll} \cos\psi_0, & j\mbox{ odd}
\\ \sin\psi_0, & j\mbox{ even}
\end{array}
\right\} \right. \nonumber\\
&& + \sum_{k=0}^{j+1} (-)^{(k+1)//2} \left(1\over \omega T\right)^k
{(j+1)!\over (j+1-k)!} \left\{
\begin{array}{ll}
\cos(\omega T+\psi_0), & k\mbox{ even}\\
\sin(\omega T+\psi_0), & k\mbox{ odd}
\end{array}
\right\} \nonumber\\
&& \left. - {1\over j+2} \left[ \cos(\omega T + \psi_0) - \cos(\psi_0)
\right] \right],
\\
g_{\delta j} &=& -{4\pi^2 \over j+1} {f_0^2r\sin\delta \over T^j}
\left[ \int_0^T {dt\over T}\, t^{j+1} \cos(\omega t + \psi_0) \right.
\nonumber\\
&& \left. - \int_0^T {dt\over T}\, t^{j+1} \int_0^T {dt\over T}\,
\cos(\omega t + \psi_0) \right]
\nonumber\\
&=& {4\pi^2 \over j+1} {f_0^2r\sin\delta \over \omega} \left[
(-)^{(j+1)//2} {(j+1)!\over (\omega T)^{j+1}} \left\{
\begin{array}{ll}
\cos\psi_0, & j\mbox{ even} \\
\sin\psi_0, & j\mbox{ odd}
\end{array}
\right\} \right. \nonumber\\
&& - \sum_{k=0}^{j+1} (-)^{k//2} \left(1\over \omega T\right)^k
{(j+1)!\over (j+1-k)!} \left\{
\begin{array}{ll}
\cos(\omega T+\psi_0), & k\mbox{ odd} \\
\sin(\omega T+\psi_0), & k\mbox{ even}
\end{array}
\right\} \nonumber\\
&& \left. + {1\over j+2} \left[ \sin(\omega T +\psi_0) - \sin(\psi_0)
\right] \right],
\\
g_{jk} &=& {(2\pi f_0T)^2 \over (j+2)(k+2)(j+k+3)},
\end{eqnarray}
where $j//2$ indicates integer division. The duration of integration is $T$,
$\omega$ and $r$ are the Earth's rotational angular velocity and equatorial
radius, $\psi_0$ is the phase at the beginning of the observation, and
$\psi_1 = \psi_0 + \omega T$ is the phase at the end.  This expression,
which neglects orbital motion, is sufficient to get a feel for the needs of
the first stage of a hierarchical search, where the coherent integration is
of order one day.

On output, the \texttt{metric->data} is arranged with the same indexing
scheme as in \texttt{CoherentMetric()}. The order of the parameters is
$(f_0, \alpha, \delta, f_1, \ldots)$.

\subsubsection*{Uses}

LALGetEarthTimes()

\subsubsection*{Notes}

The awkward series where $j//2$ and similar things appear can be represented
as incomplete gamma functions of imaginary arguments, but in the code the
trig series is the easiest to implement. I could write them more cleanly in
terms of complex exponentials, but the trig expressions are what is actually
in the code.

Orbital motion is not yet included.

\vfill{\footnotesize\input{PtoleMetricCV}}

************************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/AVFactories.h>
#include <lal/DetectorSite.h>
#include <lal/LALStdlib.h>
#include <lal/PtoleMetric.h>
#include <lal/PulsarTimes.h>

NRCSID( PTOLEMETRICC, "$Id$" );

/* Bounds on acceptable parameters, may be somewhat arbitrary */
#define MIN_DURATION (LAL_DAYSID_SI/LAL_TWOPI) /* Metric acts funny if
                                                * duration too short */
#define MIN_MAXFREQ  1.                        /* Arbitrary */

/* A private factorial function */
static int factrl( int );

/* <lalVerbatim file="PtoleMetricCP"> */
void LALPtoleMetric( LALStatus *status,
                     REAL8Vector *metric,
                     PtoleMetricIn *input )
{ /* </lalVerbatim> */
  INT2 j, k;         /* loop counters */
  REAL4 r, R;        /* radii of rotation and revolution, in seconds */
  REAL4 w, W;        /* angular frequencies of rotation and revolution, in radians per second */
  REAL4 t0, T;       /* initial time and duration of integration, in seconds */
  REAL4 psi0, psi1;  /* initial and final phases of rotation, in radians */
  REAL4 dpsi;        /* elapsed phase of rotation, in radians */
  REAL4 phi0, phi1;  /* initial and final phases of revolution, in radians */
  REAL4 dphi;        /* elapsed phase of revolution, in radians */
  REAL4 s0, c0;      /* sin(psi0) and cos(psi0) */
  REAL4 s1, c1;      /* sin(psi1) and cos(psi1) */
  REAL4 S0, C0;      /* sin(phi0) and cos(phi0) */
  REAL4 S1, C1;      /* sin(phi1) and cos(phi1) */
  REAL4 cosa, sina;  /* cosine and sine of right ascension */
  REAL4 cosb, sinb;  /* cosine and sine of detector latitude */
  REAL4 cosd, sind;  /* cosine and sine of declination */
  REAL4 cosi, sini;  /* cosine and sine of Earth's inclination */
  UINT2 dim;         /* dimension of parameter space */
  REAL4 f0;          /* maximum frequency to search */
  REAL4 Gx0, Gxa, G00, G0a, Gaa;
  PulsarTimesParamStruc tev;

  INITSTATUS( status, "LALPtoleMetric", PTOLEMETRICC );
  ATTATCHSTATUSPTR( status );

  /* Check for valid input structure. */
  ASSERT( input != NULL, status, PTOLEMETRICH_ENULL,
          PTOLEMETRICH_MSGENULL );

  /* Check for valid sky position. */
  ASSERT( input->position.system == COORDINATESYSTEM_EQUATORIAL, status,
          PTOLEMETRICH_EPARM, PTOLEMETRICH_MSGEPARM );
  ASSERT( input->position.longitude >= 0, status, PTOLEMETRICH_EPARM,
          PTOLEMETRICH_MSGEPARM );
  ASSERT( input->position.longitude < LAL_TWOPI, status, PTOLEMETRICH_EPARM,
          PTOLEMETRICH_MSGEPARM );
  ASSERT( abs(input->position.latitude) <= LAL_PI_2, status,
          PTOLEMETRICH_EPARM, PTOLEMETRICH_MSGEPARM );

  /* Check for valid duration. */
  ASSERT( input->duration > MIN_DURATION, status, PTOLEMETRICH_EPARM,
          PTOLEMETRICH_MSGEPARM );

  /* Check for valid maximum frequency. */
  ASSERT( input->maxFreq > MIN_MAXFREQ, status, PTOLEMETRICH_EPARM,
          PTOLEMETRICH_MSGEPARM );

  /* Check for valid detector location. */
  ASSERT( abs(input->site.vertexLatitudeRadians) <= LAL_PI_2, status,
          PTOLEMETRICH_EPARM, PTOLEMETRICH_MSGEPARM );
  ASSERT( abs(input->site.vertexLongitudeRadians) <= LAL_PI, status,
          PTOLEMETRICH_EPARM, PTOLEMETRICH_MSGEPARM );

  /* Check that metric has been provided. */
  if( input->spindown )
    dim = 2+input->spindown->length;
  else
    dim = 2;
  ASSERT( metric != NULL, status, PTOLEMETRICH_ENULL,
          PTOLEMETRICH_MSGENULL );
  ASSERT( metric->data != NULL, status, PTOLEMETRICH_ENULL,
          PTOLEMETRICH_MSGENULL );
  ASSERT( metric->length == (UINT4)(dim+1)*(dim+2)/2, status, PTOLEMETRICH_EDIM,
          PTOLEMETRICH_MSGEDIM );

  /* Get times until next midnight and autumnal equinox. */
  tev.epoch = input->epoch;
  tev.latitude = input->site.vertexLatitudeRadians;
  tev.longitude = input->site.vertexLongitudeRadians;
  TRY( LALGetEarthTimes( status->statusPtr, &tev ), status );

  /* These assignments make later equations more legible. */
  f0 = input->maxFreq;
  r = LAL_REARTH_SI / LAL_C_SI;
  R = LAL_AU_SI / LAL_C_SI;
  w = LAL_TWOPI / LAL_DAYSID_SI;
  W = LAL_TWOPI / LAL_YRSID_SI;
  t0 = input->epoch.gpsSeconds + 1e-9*input->epoch.gpsNanoSeconds;
  T = input->duration;
  psi0 = -w*tev.tMidnight + input->site.vertexLongitudeRadians - input->position.longitude;
  dpsi = w*T;
  psi1 = psi0+dpsi;
  phi0 = -W*tev.tAutumn;
  dphi = W*T;
  phi1 = phi0+dphi;

  /* These are for legibility and cut CPU trig time. */
  cosa = cos(input->position.longitude);
  sina = sin(input->position.longitude);
  cosd = cos(input->position.latitude);
  sind = sin(input->position.latitude);
  cosb = cos(input->site.vertexLatitudeRadians);
  sinb = sin(input->site.vertexLatitudeRadians);
  cosi = cos(LAL_IEARTH);
  sini = sin(LAL_IEARTH);
  c0 = cos(psi0);
  s0 = sin(psi0);
  c1 = cos(psi1);
  s1 = sin(psi1);
  C0 = cos(phi0);
  S0 = sin(phi0);
  C1 = cos(phi1);
  S1 = sin(phi1);

#if 0
  Gx0 = LAL_PI*(T+2*t0);
  Gxa = -LAL_TWOPI*f0*( r/dpsi*cosb*cosd*(c1-c0)
        + R/dphi*sina*cosd*(S1-S0) + R/dphi*cosi*cosa*cosd*(C1-C0) );
  G00 = 4/3*LAL_PI*LAL_PI * (T*T+3*T*t0+3*t0*t0);
  G0a = LAL_TWOPI*LAL_TWOPI*f0*(
          r/w/dpsi*cosb*cosd*(-psi1*c1 + psi0*c0 + s1 - s0)
        - R/W/dphi*sina*cosd*(phi1*S1 - phi0*S0 + C1 - C0)
        + R/W/dphi*cosi*cosa*cosd*(-phi1*C1 + phi0*C0 + S1 - S0) );
  Gaa = 2*LAL_PI*LAL_PI*f0*f0*cosd*cosd*(
          r*r*cosb*cosb*(1 - (s1*c1-s0*c0)/dpsi)
        + R*R*sina*sina*(1 + (S1*C1-S0*C0)/dphi)
        + R*R*cosi*cosi*cosa*cosa*(1 - (S1*C1-S0*C0)/dphi)
        - R*R/dphi*cosi*sina*cosa*(S1*S1 - S0*S0)
        + 2*r*R/dpsi*cosb*sina*(c1*C1 - c0*C0)
        - 2*r*R/dpsi*cosb*cosi*cosa*(c1*S1 - c0*S0) );

//  This line has overflow problems due to large t0.
//  metric->data[0] = G00 - Gx0*Gx0;
  // Why is this one not agreeing with Tev anymore?
  metric->data[0] = LAL_PI*LAL_PI/3*T*T;
  metric->data[2] = Gaa - Gxa*Gxa;
#else
  /* Spindown-spindown metric components, before projection of f0 */
  if( input->spindown )
    for (j=1; j<=dim-2; j++)
      for (k=1; k<=j; k++)
        metric->data[(k+2)+(j+2)*(j+3)/2] = pow(LAL_TWOPI*f0*T,2)/(j+2)/(k+2)/(j+k+3);

  /* angle-angle metric components, before projecting out f0 */
  /* RA-RA : 1+1*2/2 = 2 */
  metric->data[2] = pow(LAL_PI*f0*r*cosb*cosd,2);
  metric->data[2] *= 2 - 2*(s1*c1-s0*c0)/dpsi - pow(2*(c1-c0)/dpsi,2);
  /* RA-dec : 1+2*3/2 = 4 */
  metric->data[4] = -pow(LAL_PI*f0*r*cosb,2)*2*sind*cosd;
  metric->data[4] *= (s1*s1-s0*s0 + 2*(s1-s0)*(c1-c0)/dpsi)/dpsi;
  /* dec-dec : 2+2*3/2 = 5 */
  metric->data[5] = pow(LAL_PI*f0*r*cosb*sind,2);
  metric->data[5] *= 2 + 2*(s1*c1-s0*c0)/dpsi - pow(2*(s1-s0)/dpsi,2);

  /* spindown-angle metric components, before projecting out f0 */
  if( input->spindown )
    for (j=1; j<=(INT4)input->spindown->length; j++) {

      /* Spindown-RA: 1+(j+2)*(j+3)/2 */
      metric->data[1+(j+2)*(j+3)/2] = 0;
      for (k=j+1; k>=0; k--)
        metric->data[1+(j+2)*(j+3)/2] += pow(-1,(k+1)/2)*factrl(j+1)
          /factrl(j+1-k)/pow(dpsi,k)*((k%2)?s1:c1);
      metric->data[1+(j+2)*(j+3)/2] += pow(-1,j/2)/pow(dpsi,j+1)*factrl(j+1)
        *((j%2)?c0:s0);
      metric->data[1+(j+2)*(j+3)/2] -= (c1-c0)/(j+2);
      metric->data[1+(j+2)*(j+3)/2] *= -pow(LAL_TWOPI*f0,2)*r*cosb*cosd/w/(j+1);

      /* Spindown-dec: 2+(j+2)*(j+3)/2 */
      metric->data[2+(j+2)*(j+3)/2] = 0;
      for (k=j+1; k>=0; k--)
        metric->data[2+(j+2)*(j+3)/2] -= pow(-1,k/2)*factrl(j+1)/factrl(j+1-k)
          /pow(dpsi,k)*((k%2)?c1:s1);
      metric->data[2+(j+2)*(j+3)/2] += pow(-1,(j+1)/2)/pow(dpsi,j+1)
        *factrl(j+1)*((j%2)?s0:c0);
      metric->data[2+(j+2)*(j+3)/2] += (s1-s0)/(j+2);
      metric->data[2+(j+2)*(j+3)/2] *= pow(LAL_TWOPI*f0,2)*r*cosb*sind/w/(j+1);
    } /* for( j... ) */

  /* f0-f0 : 0 */
  metric->data[0] = pow(LAL_PI*T,2)/3;
  /* f0-ra : 0+1*2/2 = 1 */
  metric->data[1] = 2*pow(LAL_PI,2)*f0*r*cosb/w*cosd;
  metric->data[1] *= -c1-c0 + 2/dpsi*(s1-s0);
  /* f0-dec : 0+2*3/2 = 3 */
  metric->data[3] = -2*pow(LAL_PI,2)*f0*r*cosb/w*sind;
  metric->data[3] *= s1+s0 + 2/dpsi*(c1-c0);
  /* f0-spindown : 0+(j+2)*(j+3)/2 */
  if( input->spindown )
    for (j=1; j<=dim-2; j++)
      metric->data[(j+2)*(j+3)/2] = 2*pow(LAL_PI,2)*f0*T*T/(j+2)/(j+3);
#endif

  /* Clean up and leave. */
  DETATCHSTATUSPTR( status );
  RETURN( status );
} /* LALPtoleMetric() */


/* This is a dead simple, no error-checking, private factorial function. */

static int factrl( int arg )
{
  int ans = 1;

  if (arg==0) return 1;
  do {
    ans *= arg;
  }
  while(--arg>0);
  return ans;
} /* factrl() */
