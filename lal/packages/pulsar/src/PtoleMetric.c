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
understanding of the behavior of the parameter space. For example, the
density of patches goes as $\sin^22\delta$.

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

\subsubsection*{Notes}

Orbital motion is not yet included.

The awkward series where $j//2$ and similar things appear can be represented
as incomplete gamma functions of imaginary arguments, but in the code the
trig series is the easiest to implement. I could write them more cleanly in
terms of complex exponentials, but the trig expressions are what is actually
in the code.

\vfill{\footnotesize\input{PtoleMetricCV}}

************************************************************* </lalLaTeX> */

#include <math.h>
#include <stdlib.h>
#include <lal/AVFactories.h>
#include <lal/DetectorSite.h>
#include <lal/LALStdlib.h>
#include <lal/PtoleMetric.h>

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
  INT2 j, k;         /* Loop counters */
  REAL4 dayR, yearR; /* Amplitude of daily/yearly modulation, s */
  REAL4 dayW, yearW; /* Angular frequency of daily/yearly modulation, rad/s */
  REAL4 psi0;        /* Phase of Earth's rotation at beginning (epoch) */
  REAL4 psi1;        /* Phase of Earth's rotation at end (epoch+duration) */
  REAL4 dpsi;        /* psi1 - psi0 */
  REAL4 cosd, sind;  /* Cosine and sine of source declination */
  REAL4 s0, c0;      /* sin(psi0) and cos(psi0) */
  REAL4 s1, c1;      /* sin(psi1) and cos(psi1) */
  UINT2 dim;         /* Dimension of parameter space */
  REAL4 lat, lon;    /* latitude and longitude of detector site */

  INITSTATUS( status, "LALPtoleMetric", PTOLEMETRICC );

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
  ASSERT( metric->length == (UINT4)(dim+2)*(dim+3)/2, status, PTOLEMETRICH_EDIM,
          PTOLEMETRICH_MSGEDIM );

  /* This section saves typing and makes the equations more legible. */
  lat = input->site.vertexLatitudeRadians;
  lon = input->site.vertexLongitudeRadians;

  /* Spindown-spindown metric components, before projection */
  if( input->spindown )
    for (j=1; j<=dim-2; j++)
      for (k=1; k<=j; k++)
        metric->data[(k+2)+(j+2)*(j+3)/2] = pow(LAL_TWOPI*input->maxFreq
          *input->duration,2)/(j+2)/(k+2)/(j+k+3);

  /* Is equatorial radius good enough? */
  dayR = LAL_REARTH_SI / LAL_C_SI * cos(lat);
  /* Is a circular orbit good enough? */
  yearR = LAL_AU_SI / LAL_C_SI;
  /* Is 1994 sidereal day good enough? */
  dayW = LAL_TWOPI / LAL_DAYSID_SI;
  /* Is 1994 sidereal year good enough? */
  yearW = LAL_TWOPI / LAL_YRSID_SI;

  /* This assumes that gpsSeconds=0 is when the source is on the meridian */
  psi0 = dayW*input->epoch.gpsSeconds - input->position.longitude + lon;
  dpsi = dayW*input->duration;
  psi1 = psi0+dpsi;

  /* Trig-savers */
  cosd = cos(input->position.latitude);
  sind = sin(input->position.latitude);
  c0 = cos(psi0);
  s0 = sin(psi0);
  c1 = cos(psi0+dpsi);
  s1 = sin(psi0+dpsi);

  /* Angle-angle metric components, before projection */
  /* RA-RA : 1+1*2/2 = 2 */
  metric->data[2] = pow(LAL_PI*input->maxFreq*dayR*cosd,2);
  metric->data[2] *= 2 - 2*(s1*c1-s0*c0)/dpsi - pow(2*(c1-c0)/dpsi,2);
  /* RA-dec : 1+2*3/2 = 4 */
  metric->data[4] = -pow(LAL_PI*input->maxFreq*dayR,2)*2*sind*cosd;
  metric->data[4] *= (s1*s1-s0*s0 + 2*(s1-s0)*(c1-c0)/dpsi)/dpsi;
  /* dec-dec : 2+2*3/2 = 5 */
  metric->data[5] = pow(LAL_PI*input->maxFreq*dayR*sind,2);
  metric->data[5] *= 2 + 2*(s1*c1-s0*c0)/dpsi - pow(2*(s1-s0)/dpsi,2);

  /* Spindown-angle metric components, before projection */
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
      metric->data[1+(j+2)*(j+3)/2] *= -pow(LAL_TWOPI*input->maxFreq,2)*dayR
        *cosd/dayW/(j+1);

      /* Spindown-dec: 2+(j+2)*(j+3)/2 */
      metric->data[2+(j+2)*(j+3)/2] = 0;
      for (k=j+1; k>=0; k--)
        metric->data[2+(j+2)*(j+3)/2] -= pow(-1,k/2)*factrl(j+1)/factrl(j+1-k)
          /pow(dpsi,k)*((k%2)?c1:s1);
      metric->data[2+(j+2)*(j+3)/2] += pow(-1,(j+1)/2)/pow(dpsi,j+1)
        *factrl(j+1)*((j%2)?s0:c0);
      metric->data[2+(j+2)*(j+3)/2] += (s1-s0)/(j+2);
      metric->data[2+(j+2)*(j+3)/2] *= pow(LAL_TWOPI*input->maxFreq,2)*dayR
        *sind/dayW/(j+1);
    } /* for( j... ) */

  /* f0-f0 : 0 */
  metric->data[0] = pow(LAL_PI*input->duration,2)/3;
  /* f0-ra : 0+1*2/2 = 1 */
  metric->data[1] = 2*pow(LAL_PI,2)*input->maxFreq*dayR/dayW*cosd;
  metric->data[1] *= -c1-c0 + 2/dpsi*(s1-s0);
  /* f0-dec : 0+2*3/2 = 3 */
  metric->data[3] = -2*pow(LAL_PI,2)*input->maxFreq*dayR/dayW*sind;
  metric->data[3] *= s1+s0 + 2/dpsi*(c1-c0);
  /* f0-spindown : 0+(j+2)*(j+3)/2 */
  if( input->spindown )
    for (j=1; j<=dim-2; j++)
      metric->data[(j+2)*(j+3)/2] = 2*pow(LAL_PI,2)*input->maxFreq
        *pow(input->duration,2)/(j+2)/(j+3);

  /* All done */
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
