/**************************************** <lalVerbatim file="PtoleMetricCV">
Author: Jones D. I.,   Owen, B. J.
$Id$
********************************************************** </lalVerbatim> */

/**************************************************************** <lalLaTeX>

\subsection{Module \texttt{LALPtoleMetric.c}}
\label{ss:PtoleMetric}

Computes metric components for a pulsar search in the ``Ptolemaic''
approximation.  Both the Earth's spin and orbit are included.

\subsubsection*{Prototypes}
\input{PtoleMetricCP}
\idx{PtoleMetric}

\subsubsection*{Description}

This function computes metric components in a way that yields results
very similar to those of \texttt{CoherentMetric} called with the
output of \texttt{TBaryPtolemaic}. The CPU demand, however, is less,
and the metric components can be expressed analytically, lending
themselves to better understanding of the behavior of the parameter
space.  For short durations (less than about 70000 seconds) a Taylor
series expansion is used to improve the accuracy with which several
terms are computed.

\subsubsection*{Algorithm}

For speed and checking reasons, a minimum of numerical computation is
involved. The metric components can be expressed analytically (though
not tidily) in terms of trig functions and are much too large to be
worth writing down here.  Jones \& Owen will write up a note on the
calculation for submission to the preprint server.

The function \texttt{GetEarthTimes} is used to calculate the spin and
rotational phase of the Earth at the beginning of the observation.

On output, the \texttt{metric->data} is arranged with the same indexing
scheme as in \texttt{CoherentMetric()}. The order of the parameters is
$(f_0, \alpha, \delta)$.

\subsubsection*{Uses}

LALGetEarthTimes()

\subsubsection*{Notes}

The analytic metric components were derived separately by the
two authors and found to agree.  Also, the output of this function has
been compared against that of the function combination (CoherentMetric
+ TDBaryPtolemaic), which is a non-analytic implementation of the
Ptolemaic approximation, and found to agree up to the fourth
significant figure or better.

Spindown is not yet included.

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
INT2 j, k;              /* Loop counters */
REAL8 R_o, R_s;         /* Amplitude of daily/yearly modulation, s */
REAL8 lat, lon;         /* latitude and longitude of detector site */
REAL8 omega_s, omega_o; /* Ang freq of daily/yearly modulation, rad/s */
REAL8 phi_o_i;          /* Phase of Earth's orbit at begining (epoch) */
REAL8 phi_o_f;          /* Phase of Earth's orbit at end (epoch+duration) */
REAL8 phi_s_i;          /* Phase of Earth's spin at beginning (epoch) */
REAL8 phi_s_f;          /* Phase of Earth's spin at end (epoch+duration) */
REAL8 cos_l, sin_l, sin_2l;  /* Trig fns for detector lattitude */
REAL8 cos_i, sin_i;     /* Trig fns for inclination of Earth's spin */
REAL8 cos_a, sin_a, sin_2a;  /* Trig fns for source RA */
REAL8 cos_d, sin_d, sin_2d;  /* Trig fns for source declination */
REAL8 A[22];     /* Array of intermediate quantities */
REAL8 B[10];     /* Array of intermediate quantities */
UINT2 dim;         /* Dimension of parameter space */
PulsarTimesParamStruc zero_phases; /* Needed to calculate phases of spin*/
                                   /* and orbit at t_gps =0             */
REAL8Vector *big_metric; /* 10-dim metric in (phi,f,a,d) for internal use */
REAL8 T;  /* Duration of observation */
REAL8 f; /* Frequency */
/* Some useful short-hand notations involving the orbital phase: */
REAL8 D_p_o;
REAL8 sin_p_o;
REAL8 cos_p_o;
REAL8 D_sin_p_o;
REAL8 D_cos_p_o;
REAL8 D_sin_2p_o;
REAL8 D_cos_2p_o;
/* Some useful short-hand notations involving the spin phase: */
REAL8 D_p_s;
REAL8 sin_p_s;
REAL8 cos_p_s;
REAL8 D_sin_p_s;
REAL8 D_cos_p_s;
REAL8 D_sin_2p_s;
REAL8 D_cos_2p_s;
/* Some useful short-hand notations involving spin and orbital phases: */
REAL8 D_sin_p_o_plus_s;
REAL8 D_p_o_plus_s;
REAL8 D_sin_p_o_minus_s;
REAL8 D_p_o_minus_s;
REAL8 D_cos_p_o_plus_s;
REAL8 D_cos_p_o_minus_s;
/* Quantities related to the short-time Tatlor expansions : */
REAL8 D_p_o_crit;       /* The orbital phase BELOW which series used.  */ 
REAL8 sum1, sum2, sum3;
REAL8 next1, next2, next3;
 INT2 j_max;  /* Number of terms in series */


INITSTATUS( status, "LALPtoleMetric", PTOLEMETRICC );

/* A 4-dim (phi, f, a, d) metric for internal use only:  */
big_metric = NULL;
LALDCreateVector( status, &big_metric, 10 );

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
ASSERT( abs(input->site->frDetector.vertexLatitudeRadians) <= LAL_PI_2, status,
	PTOLEMETRICH_EPARM, PTOLEMETRICH_MSGEPARM );
ASSERT( abs(input->site->frDetector.vertexLongitudeRadians) <= LAL_PI, status,
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

/* Detector location: */
lat = input->site->frDetector.vertexLatitudeRadians;
lon = input->site->frDetector.vertexLongitudeRadians;
cos_l = cos(lat);
sin_l = sin(lat);
sin_2l = sin(2*lat);

/* Inclination of Earth's spin-axis to ecliptic */
sin_i = sin(LAL_IEARTH);
cos_i = cos(LAL_IEARTH);

/* Radii of circular motions in seconds:  */
R_s = LAL_REARTH_SI / LAL_C_SI;
R_o = LAL_AU_SI / LAL_C_SI;

/* To switch off orbital motion uncomment this: */
/* R_o = 0.0; */

/* To switch off spin motion uncomment this: */
/* R_s = 0.0; */

/* Angular velocities: */
omega_s = LAL_TWOPI / LAL_DAYSID_SI;
omega_o = LAL_TWOPI / LAL_YRSID_SI;

/* Duration of observation */
T = input->duration;

/* Frequency: */
f = input->maxFreq;

/* Source RA: */
sin_a  = sin(input->position.longitude);
cos_a  = cos(input->position.longitude);
sin_2a = sin(2*(input->position.longitude));

/* Source declination: */
sin_d  = sin(input->position.latitude);
cos_d  = cos(input->position.latitude);
sin_2d = sin(2*(input->position.latitude));

/* Calculation of phases of spin and orbit at start: */
zero_phases.epoch.gpsSeconds = input->epoch.gpsSeconds;
zero_phases.epoch.gpsNanoSeconds = input->epoch.gpsNanoSeconds;
LALGetEarthTimes( status, &zero_phases);
phi_o_i = -zero_phases.tAutumn/LAL_YRSID_SI*LAL_TWOPI;
phi_s_i = -zero_phases.tMidnight/LAL_DAYSID_SI*LAL_TWOPI + lon;


/* Quantities involving the orbital phase: */
phi_o_f   = phi_o_i + omega_o*T;
D_p_o     = omega_o*T;
sin_p_o   = sin(phi_o_f);
cos_p_o   = cos(phi_o_f);
D_sin_p_o = (sin(phi_o_f) - sin(phi_o_i))/D_p_o;
D_cos_p_o = (cos(phi_o_f) - cos(phi_o_i))/D_p_o;
D_sin_2p_o = (sin(2*phi_o_f) - sin(2*phi_o_i))/2.0/D_p_o;
D_cos_2p_o = (cos(2*phi_o_f) - cos(2*phi_o_i))/2.0/D_p_o;

/* Quantities involving the spin phase: */
phi_s_f    = phi_s_i + omega_s*T;
D_p_s      = omega_s*T;
sin_p_s    = sin(phi_s_f);
cos_p_s    = cos(phi_s_f);
D_sin_p_s  = (sin(phi_s_f) - sin(phi_s_i))/D_p_s;
D_cos_p_s  = (cos(phi_s_f) - cos(phi_s_i))/D_p_s;
D_sin_2p_s = (sin(2*phi_s_f) - sin(2*phi_s_i))/2.0/D_p_s;
D_cos_2p_s = (cos(2*phi_s_f) - cos(2*phi_s_i))/2.0/D_p_s;

/* Some mixed quantities: */
D_sin_p_o_plus_s  
  = (sin(phi_o_f+phi_s_f) - sin(phi_o_i+phi_s_i))/(D_p_o + D_p_s);
D_p_o_plus_s      = D_p_o + D_p_s;
D_sin_p_o_minus_s 
  = (sin(phi_o_f-phi_s_f) - sin(phi_o_i-phi_s_i))/(D_p_o - D_p_s);
D_p_o_minus_s     = D_p_o - D_p_s;
D_cos_p_o_plus_s  
  = (cos(phi_o_f+phi_s_f) - cos(phi_o_i+phi_s_i))/(D_p_o + D_p_s);
D_cos_p_o_minus_s 
  = (cos(phi_o_f-phi_s_f) - cos(phi_o_i-phi_s_i))/(D_p_o - D_p_s);



/***************************************************************/
/* Small D_p_o overwrite:                                      */
/***************************************************************/
 j_max = 7;
 D_p_o_crit = 1.4e-2; /* This corresponds to about 70000 seconds */

 sum1  = next1 = D_p_o/2.0;
 sum2  = next2 = D_p_o/3.0/2.0;
 sum3  = next3 = D_p_o;

 for(j=1; j<=j_max; j++)
   {
     next1 *= -pow(D_p_o,2.0)/(2.0*j+1.0)/(2.0*j+2.0);
     sum1  += next1; 
     next2 *= -pow(D_p_o,2.0)/(2.0*j+2.0)/(2.0*j+3.0);
     sum2  += next2; 
     next3 *= -pow(2.0*D_p_o,2.0)/(2.0*j+1.0)/(2.0*j+2.0);
     sum3  += next3;
   }
 
 if(D_p_o < D_p_o_crit)
   {
     D_sin_p_o = sin(phi_o_f)*sum1 + cos(phi_o_f)*sin(D_p_o)/D_p_o;
     D_cos_p_o = cos(phi_o_f)*sum1 - sin(phi_o_f)*sin(D_p_o)/D_p_o;
     D_sin_2p_o 
       = sin(2.0*phi_o_f)*sum3 + cos(2.0*phi_o_f)*sin(2.0*D_p_o)/2.0/D_p_o;
     D_cos_2p_o 
       = cos(2.0*phi_o_f)*sum3 - sin(2.0*phi_o_f)*sin(2.0*D_p_o)/2.0/D_p_o;
   }
/****************************************************************/


/* The A[i] quantities: */
A[1] = 
R_o*D_sin_p_o + R_s*cos_l*D_sin_p_s;

A[2] = 
R_o*cos_i*D_cos_p_o + R_s*cos_l*D_cos_p_s;

A[3] =
-R_o*sin_i*D_cos_p_o + R_s*sin_l;

A[4] =
R_o*(sin_p_o/D_p_o + D_cos_p_o/D_p_o);

A[5] =
R_s*(sin_p_s/D_p_s + D_cos_p_s/D_p_s);

A[6] =
R_o*(-cos_p_o/D_p_o + D_sin_p_o/D_p_o);

/* Special overwrite for A4 and A6: *********************/
 if(D_p_o < D_p_o_crit)
   {
     A[4] = R_o*(cos(phi_o_f)*sum1/D_p_o + sin(phi_o_f)*sum2);
     A[6] = R_o*(sin(phi_o_f)*sum1/D_p_o - cos(phi_o_f)*sum2);
   }
 /***********************************************************/

A[7] =
R_s*(-cos_p_s/D_p_s + D_sin_p_s/D_p_s);

A[8] =
R_o*R_o*(1 + D_sin_2p_o);

A[9] =
R_o*R_s*(D_sin_p_o_minus_s + D_sin_p_o_plus_s);

A[10] =
R_s*R_s*(1 + D_sin_2p_s);

A[11] =
R_o*R_o*D_cos_2p_o;

A[12] =
R_o*R_s*(-D_cos_p_o_minus_s + D_cos_p_o_plus_s);

A[13] =
R_o*R_s*(D_cos_p_o_minus_s + D_cos_p_o_plus_s);

A[14] =
R_s*R_s*D_cos_2p_s;

A[15] =
R_o*R_o*(1 - D_sin_2p_o);

A[16] =
R_o*R_s*(D_sin_p_o_minus_s - D_sin_p_o_plus_s);

A[17] =
R_s*R_s*(1 - D_sin_2p_s);

A[18] =
R_o*R_s*D_sin_p_o;

A[19] =
R_s*R_s*D_sin_p_s;

A[20] =
R_o*R_s*D_cos_p_o;

A[21] =
R_s*R_s*D_cos_p_s;


/* The B[i] quantities: */
B[1] =
A[4] + A[5]*cos_l;

B[2] =
A[6]*cos_i + A[7]*cos_l;

B[3] =
A[6]*sin_i + R_s*sin_l/2;

B[4] =
A[8] + 2*A[9]*cos_l + A[10]*cos_l*cos_l;

B[5] =
A[11]*cos_i + A[12]*cos_l + A[13]*cos_i*cos_l + A[14]*cos_l*cos_l;

B[6] =
A[15]*cos_i*cos_i +2*A[16]*cos_i*cos_l + A[17]*cos_l*cos_l;

B[7] =
-A[11]*sin_i + 2*A[18]*sin_l - A[13]*sin_i*cos_l + A[19]*sin_2l; 

B[8] =
A[15]*sin_i*cos_i - 2*A[20]*cos_i*sin_l + A[16]*sin_i*cos_l - A[21]*sin_2l;

B[9] =
A[15]*sin_i*sin_i - 4*A[20]*sin_i*sin_l + 2*R_s*R_s*sin_l*sin_l;


/* The 4-dim metric components: */
/* g_pp = */
big_metric->data[0] =
1;

/* g_pf = */
big_metric->data[1] =
LAL_PI*T;

/* g_pa = */
big_metric->data[3] =
-LAL_TWOPI*f*cos_d*(A[1]*sin_a + A[2]*cos_a);

/* g_pd = */
big_metric->data[6] =
LAL_TWOPI*f*(-A[1]*sin_d*cos_a + A[2]*sin_d*sin_a + A[3]*cos_d);

/* g_ff = */
big_metric->data[2] =
pow(LAL_TWOPI*T, 2)/3;

/* g_fa = */
big_metric->data[4] =
pow(LAL_TWOPI,2)*f*cos_d*T*(-B[1]*sin_a + B[2]*cos_a);

/* g_fd = */
big_metric->data[7] =
pow(LAL_TWOPI,2)*f*T*(-B[1]*sin_d*cos_a - B[2]*sin_d*sin_a + B[3]*cos_d);

/* g_aa = */
big_metric->data[5] =
2*pow(LAL_PI*f*cos_d,2) * (B[4]*sin_a*sin_a + B[5]*sin_2a + B[6]*cos_a*cos_a);

/* g_ad =  */
big_metric->data[8] =
2*pow(LAL_PI*f,2)*cos_d*(B[4]*sin_a*cos_a*sin_d - B[5]*sin_a*sin_a*sin_d
-B[7]*sin_a*cos_d + B[5]*cos_a*cos_a*sin_d - B[6]*sin_a*cos_a*sin_d
+B[8]*cos_a*cos_d);

/* g_dd = */
big_metric->data[9] =
2*pow(LAL_PI*f,2)*(B[4]*pow(cos_a*sin_d,2) + B[6]*pow(sin_a*sin_d,2)
+B[9]*pow(cos_d,2) - B[5]*sin_2a*pow(sin_d,2) - B[8]*sin_a*sin_2d
-B[7]*cos_a*sin_2d);


/**********************************************************/
/* Spin-down stuff not consistent with rest of code - don't uncomment! */
/* Spindown-spindown metric components, before projection */
if( input->spindown )
     for (j=1; j<=dim-2; j++)
     for (k=1; k<=j; k++)
     metric->data[(k+2)+(j+2)*(j+3)/2] = pow(LAL_TWOPI*input->maxFreq
     *input->duration,2)/(j+2)/(k+2)/(j+k+3);

/* Spindown-angle metric components, before projection */
if( input->spindown )
     for (j=1; j<=(INT4)input->spindown->length; j++) {

/* Spindown-RA: 1+(j+2)*(j+3)/2 */
metric->data[1+(j+2)*(j+3)/2] = 0;
for (k=j+1; k>=0; k--)
     metric->data[1+(j+2)*(j+3)/2] += pow(-1,(k+1)/2)*factrl(j+1)
     /factrl(j+1-k)/pow(D_p_s,k)*((k%2)?sin_p_s:cos_p_s);
     metric->data[1+(j+2)*(j+3)/2] += pow(-1,j/2)/pow(D_p_s,j+1)*factrl(j+1)
     *((j%2)?cos(phi_s_i):sin(phi_s_i));
     metric->data[1+(j+2)*(j+3)/2] -= (cos_p_s-cos(phi_s_i))/(j+2);
     metric->data[1+(j+2)*(j+3)/2] *= -pow(LAL_TWOPI*input->maxFreq,2)*R_s*cos_l
     *cos_d/omega_s/(j+1);

/* Spindown-dec: 2+(j+2)*(j+3)/2 */
metric->data[2+(j+2)*(j+3)/2] = 0;
for (k=j+1; k>=0; k--)
     metric->data[2+(j+2)*(j+3)/2] -= pow(-1,k/2)*factrl(j+1)/factrl(j+1-k)
     /pow(D_p_s,k)*((k%2)?cos_p_s:sin_p_s);
     metric->data[2+(j+2)*(j+3)/2] += pow(-1,(j+1)/2)/pow(D_p_s,j+1)
     *factrl(j+1)*((j%2)?sin(phi_s_i):cos(phi_s_i));
     metric->data[2+(j+2)*(j+3)/2] += (sin_p_s-sin(phi_s_i))/(j+2);
     metric->data[2+(j+2)*(j+3)/2] *= pow(LAL_TWOPI*input->maxFreq,2)*R_s*cos_l
     *sin_d/omega_s/(j+1);
} /* for( j... ) */

  /* f0-spindown : 0+(j+2)*(j+3)/2 */
if( input->spindown )
     for (j=1; j<=dim-2; j++)
     metric->data[(j+2)*(j+3)/2] = 2*pow(LAL_PI,2)*input->maxFreq
     *pow(input->duration,2)/(j+2)/(j+3);
/*************************************************************/


/* Project down to 3-dim metric: */
metric->data[0] =  big_metric->data[2] 
    - big_metric->data[1]*big_metric->data[1]/big_metric->data[0];
metric->data[1] =  big_metric->data[4] 
    - big_metric->data[1]*big_metric->data[3]/big_metric->data[0];
metric->data[2] =  big_metric->data[5] 
    - big_metric->data[3]*big_metric->data[3]/big_metric->data[0];
metric->data[3] =  big_metric->data[7] 
    - big_metric->data[6]*big_metric->data[1]/big_metric->data[0];
metric->data[4] =  big_metric->data[8] 
    - big_metric->data[6]*big_metric->data[3]/big_metric->data[0];
metric->data[5] =  big_metric->data[9] 
    - big_metric->data[6]*big_metric->data[6]/big_metric->data[0]; 


LALDDestroyVector( status, &big_metric ); 
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
