/*

*  Copyright (C) 2007 David McKechan, Thomas Cokelaer
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

/************************** <lalVerbatim file="GeneratePPNAmpCorInspiralCV">
Author: Creighton, T. D., McKechan David, Van Den Broeck Chris
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\providecommand{\lessim}{\stackrel{<}{\scriptstyle\sim}}

\subsection{Module \texttt{GeneratePPNampCorInspiral.c}}
\label{ss:GeneratePPNAmpCorInspiral.c}

Computes a parametrized post-Newtonian inspiral waveform
with ampltidude corrections.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{GeneratePPNAmpCorInspiralCP}
\idx{LALGeneratePPNampCorInspiral()}

\subsubsection*{Description}

See GeneratePPNInspiral.c

Phase computed to 3.5PN.
Amplitude computed to 2.5PN.

The ampitude corrrected gravitaional wave polarizations $h_+$ and $h_x$
are stored in output.h.

Warning! output.a is used to store the first three harmonics in
alternate values, (i.e. [a1(0),a2(0),a3(0),a1(dt),a2(dt),a3(dt)...]) as
this will be used for filtering with higher harmonic waveforms.

Although $h_{+,\times} are computed, $output.phi is also required for filtering.


\subsubsection*{Algorithm}


\subsubsection*{Uses}
\begin{verbatim}
LALMalloc()                   LALFree()
LALSCreateVectorSequence()    LALSDestroyVectorSequence()
LALSCreateVector()            LALSDestroyVector()
LALDCreateVector()            LALDDestroyVector()
LALSBisectionFindRoot()       snprintf()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GeneratePPNAmpCorInspiralCV}}

******************************************************* </lalLaTeX> */

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/FindRoot.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/GeneratePPNInspiral.h>

NRCSID( GENERATEPPNAMPCORINSPIRALC, "$Id$" );

/* Define some constants used in this module. */
#define MAXORDER 8               /* Maximum number of N and PN terms */
#define AMPMAXORDER 6            /* Maximum PN order in amplitude (plus one) */
#define BUFFSIZE 1024            /* Number of timesteps buffered */
#define ACCURACY (1.0e-8)        /* Accuracy of root finder */
#define TWOTHIRDS (0.6666666667) /* 2/3 */
#define ONEMINUSEPS (0.99999)    /* Something close to 1 */

/*
   A macro to computing the (normalized) frequency.  It appears in
   many places, including in the main loop, and I don't want the
   overhead of a function call.  The following variables are required
   to be defined and set outside of the macro:

   REAL4 c0, c1, c2, c3, c4, c5, c6, c7, p6; PN frequency coefficients
   BOOLEAN b0, b1, b2, b3, b4, b5, b6, b7;   whether to include each PN term

   The following variables must be defined outside the macro, but are
   set inside it:

   REAL4 x2, x3, x4, x5, x6, x7;  the input x raised to power
   2, 3, 4, 5, 6 and 7
*/

#define FREQ( f, x )                                                 \
do {                                                                 \
  x2 = (x)*(x);                                                      \
  x3 = x2*(x);                                                       \
  x4 = x3*(x);                                                       \
  x5 = x4*(x);                                                       \
  x6 = x5*(x);                                                       \
  x7 = x6*(x);                                                       \
  (f) = 0;                                                           \
  if ( b0 )                                                          \
    (f) += c0;                                                       \
  if ( b1 )                                                          \
    (f) += c1*(x);                                                   \
  if ( b2 )                                                          \
    (f) += c2*x2;                                                    \
  if ( b3 )                                                          \
    (f) += c3*x3;                                                    \
  if ( b4 )                                                          \
    (f) += c4*x4;                                                    \
  if ( b5 )                                                          \
    (f) += c5*x5;                                                    \
  if ( b6 )                                                          \
    (f) += (c6 + p6*( -107.0/2240.0*(-8.0)*(log(2.0*x))))*x6;        \
  if ( b7 )                                                          \
    (f) += c7*x7;                                                    \
  (f) *= x3;                                                         \
} while (0)

/* Definition of a data structure used by FreqDiff() below. */
typedef struct tagFreqDiffParamStruc
{
  REAL4 *c;   /* PN coefficients of frequency series */
  BOOLEAN *b; /* whether to include each PN term */
  REAL4 p6;   /* synonym for p[6] */
  REAL4 y0;   /* normalized frequency being sought */
} FreqDiffParamStruc;

/* A function to compute the difference between the current and
   requested normalized frequency, used by the root bisector. */
static void
FreqDiff( LALStatus *stat, REAL4 *y, REAL4 x, void *p )
{
  FreqDiffParamStruc *par;                  /* *p cast to its proper type */
  REAL4 c0, c1, c2, c3, c4, c5, c6, c7, p6; /* PN frequency coefficients */
  BOOLEAN b0, b1, b2, b3, b4, b5, b6, b7;   /* whether each order is nonzero */
  REAL4 x2, x3, x4, x5, x6, x7;             /* x^2, x^3, x^4, x^5, x^6, & x^7 */

  INITSTATUS( stat, "FreqDiff", GENERATEPPNAMPCORINSPIRALC );
  ASSERT( p, stat, 1, "Null pointer" );

  /* Set constants used by FREQ() macro. */
  par = (FreqDiffParamStruc *)( p );
  c0 = par->c[0];
  c1 = par->c[1];
  c2 = par->c[2];
  c3 = par->c[3];
  c4 = par->c[4];
  c5 = par->c[5];
  c6 = par->c[6];
  c7 = par->c[7];
  b0 = par->b[0];
  b1 = par->b[1];
  b2 = par->b[2];
  b3 = par->b[3];
  b4 = par->b[4];
  b5 = par->b[5];
  b6 = par->b[6];
  b7 = par->b[7];
  p6 = par->p6;

  /* Evaluate frequency and compare with reference. */
  FREQ( *y, x );
  *y -= par->y0;
  RETURN( stat );
}

/* Definition of a data buffer list for storing the waveform. */
typedef struct tagPPNInspiralBuffer
{
  REAL4 h[2*BUFFSIZE];               /* polarisation data */
  REAL4 a[3*BUFFSIZE];               /* first three harmonics data */
  REAL8 phi[BUFFSIZE];               /* phase data */
  REAL4 f[BUFFSIZE];                 /* frequency data */
  struct tagPPNInspiralBuffer *next; /* next buffer in list */
} PPNInspiralBuffer;

/* Definition of a macro to free the tail of said list, from a given
   node onward. */
#define FREELIST( node )                                             \
do {                                                                 \
  PPNInspiralBuffer *herePtr = (node);                               \
  while ( herePtr ) {                                                \
    PPNInspiralBuffer *lastPtr = herePtr;                            \
    herePtr = herePtr->next;                                         \
    LALFree( lastPtr );                                              \
  }                                                                  \
} while (0)


/*********************************************************************
 * MAIN FUNCTION                                                     *
 *********************************************************************/

/* <lalVerbatim file="GeneratePPNAmpCorInspiralCP"> */
void
LALGeneratePPNAmpCorInspiral(
                              LALStatus     *stat,
                              CoherentGW    *output,
                              PPNParamStruc *params
                            )
{ /* </lalVerbatim> */

  /* System-derived constants. */
  BOOLEAN b0, b1, b2, b3, b4, b5, b6, b7; /* whether each order is nonzero */
  BOOLEAN b[MAXORDER];                    /* vector of above coefficients */
  REAL4 c0, c1, c2, c3, c4, c5, c6, c7;   /* PN frequency coefficients */
  REAL4 c[MAXORDER];                      /* vector of above coefficients */
  REAL4 d0, d1, d2, d3, d4, d5, d6, d7;   /* PN phase coefficients */
  REAL4 e0, e1, e2, e3, e4, e5, e6, e7;   /* PN dy/dx coefficients */
  REAL4 p[MAXORDER];                      /* PN parameter values in phase */
  REAL4 q[AMPMAXORDER];                   /* PN parameter values in amplitude */
  REAL4 p6;                        /* synonym for p[6] */
  UINT4 ampOrder;                  /* Amplitude Order */
  INT4 harmonics;                  /* Number of harmonics */
  REAL4 mTot, mu;                  /* total mass and reduced mass */
  REAL4 eta, etaInv;               /* mass ratio and its inverse */
  REAL4 phiC;                      /* phase at coalescence */
  REAL4 cosI, cos2I, cos4I, cos6I; /* cosine of system inclination */
  REAL4 sinI, sin2I, sin4I, sin5I; /* sine of system inclination */

  REAL4 fFac;          /* SI normalization for f and t */
  REAL4 f2aFac;        /* factor multiplying f in amplitude function */
  REAL4 fthree, ffour, ffive, fsix, fseven; /* powers in f2a to speed up
                                               waveform construction */
  REAL4 preFac;        /* Overall prefactor in waveforms */
  REAL4 delta;         /* relative mass difference */
  REAL4 sd, scd;       /* sinI*delta, sd*cosI*/

  /* Integration parameters. */
  UINT4 i;               /* index over PN terms */
  UINT4 j;               /* index of leading nonzero PN term */
  UINT4 n, nMax;         /* index over timesteps, and its maximum + 1 */
  UINT4 nNext;           /* index where next buffer starts */
  REAL8 t, t0, dt;       /* dimensionless time, start time, and increment */
  REAL4 tStop = 0.0625;  /* time when orbit reaches minimum radius */
  REAL4 x, xStart, xMax; /* x = t^(-1/8), and its maximum range */
  REAL4 y, yStart, yMax; /* normalized frequency and its range and start time */
  REAL4 yOld, dyMax;     /* previous timestep y, and maximum y - yOld */
  REAL4 *f;              /* pointer to generated frequency data */
  REAL4 x2, x3, x4, x5, x6, x7;  /* x^2, x^3, x^4, x^5, x^6 and x^7 */

  /* Harmonic terms */
  REAL4 a1Pthree, a1Pfive, a1Psix, a1Pseven, a1Pmixsix;
  REAL4 a2Ptwo, a2Pfour, a2Pfive, a2Psix, a2Pseven, a2Pmixseven;
  REAL4 a3Pthree, a3Pfive, a3Psix, a3Pseven, a3Pmixsix;
  REAL4 a4Pfour, a4Psix, a4Pseven, a4Pmixseven;
  REAL4 a5Pfive, a5Pseven;
  REAL4 a6Psix;
  REAL4 a7Pseven;
  REAL4 a1Cthree, a1Cfive, a1Csix, a1Cseven, a1Cmixsix;
  REAL4 a2Ctwo, a2Cfour, a2Cfive, a2Csix, a2Cseven, a2Cmixseven;
  REAL4 a3Cthree, a3Cfive, a3Csix, a3Cseven, a3Cmixsix;
  REAL4 a4Cfour, a4Csix, a4Cseven, a4Cmixseven;
  REAL4 a5Cfive, a5Cseven;
  REAL4 a6Csix;
  REAL4 a7Cseven;

  REAL4 a1, a2, a3, a4, a5, a6, a7; /* generated amplitudes of harmonics */
  REAL4 a1mix, a2mix, a3mix, a4mix; /* generated amplitudes of harmonic
                                       mixed terms */
  REAL4 *h;                /* pointer to generated hplus and hcross */
  REAL4 *a;                /* pointer to generated first three plus harmonics*/
  REAL8 *phi;              /* pointer to generated phase data */

  PPNInspiralBuffer *head, *here; /* pointers to buffered data */

  INITSTATUS( stat, "LALGeneratePPNAmpCorInspiral", GENERATEPPNAMPCORINSPIRALC);
  ATTATCHSTATUSPTR( stat );

  /*******************************************************************
   * CHECK INPUT PARAMETERS                                          *
   *******************************************************************/

  /* Dumb initialization to shut gcc up. */
  head = here = NULL;
  b0 = b1 = b2 = b3 = b4 = b5 = b6 = b7 = 0.0;
  c0 = c1 = c2 = c3 = c4 = c5 = c6 = c7 = 0.0;
  d0 = d1 = d2 = d3 = d4 = d5 = d6 = d7 = 0.0;

  /* Make sure parameter and output structures exist. */
  ASSERT( params, stat, GENERATEPPNINSPIRALH_ENUL,
    GENERATEPPNINSPIRALH_MSGENUL );
  ASSERT( output, stat, GENERATEPPNINSPIRALH_ENUL,
    GENERATEPPNINSPIRALH_MSGENUL );

  /* Make sure output fields don't exist. */
  ASSERT( !( output->h ), stat, GENERATEPPNINSPIRALH_EOUT,
    GENERATEPPNINSPIRALH_MSGEOUT );
  ASSERT( !( output->a ), stat, GENERATEPPNINSPIRALH_EOUT,
    GENERATEPPNINSPIRALH_MSGEOUT );
  ASSERT( !( output->f ), stat, GENERATEPPNINSPIRALH_EOUT,
    GENERATEPPNINSPIRALH_MSGEOUT );
  ASSERT( !( output->phi ), stat, GENERATEPPNINSPIRALH_EOUT,
    GENERATEPPNINSPIRALH_MSGEOUT );
  ASSERT( !( output->shift ), stat, GENERATEPPNINSPIRALH_EOUT,
    GENERATEPPNINSPIRALH_MSGEOUT );

  /* Get PN parameters, if they are specified; otherwise use
     2PN */
  if ( params->ppn ) {
    ASSERT( params->ppn->data, stat, GENERATEPPNINSPIRALH_ENUL,
        GENERATEPPNINSPIRALH_MSGENUL );
    j = params->ppn->length;
    if ( j > MAXORDER )
      j = MAXORDER;
    for ( i = 0; i < j; i++ )
      p[i] = params->ppn->data[i];
    for ( ; i < MAXORDER; i++ )
      p[i] = 0.0;
  }
  else {
    p[0] = 1.0;
    p[1] = 0.0;
    p[2] = 1.0;
    p[3] = 1.0;
    p[4] = 1.0;
    for ( i = 5; i < MAXORDER; i++ )
      p[i] = 0.0;
  }



  /* Set PN parameters for amplitude */
  if ( (params->ampOrder < 0) || (params->ampOrder >= AMPMAXORDER) )
  {
    if (params->ppn->length - 1 >= 5 )
      ampOrder = 5;
    else
      ampOrder = params->ppn->length - 1;
  }
  else
    ampOrder = params->ampOrder;

  q[0] = 1.0;
  for(i = 1; i < AMPMAXORDER; i++)
    q[i] = ( i <= ampOrder? 1.0 : 0.0 );

  /* Set number of harmonics in accordance with params->ampOrder*/
  /* Dominant harmonic is the 2nd */
  harmonics = ampOrder + 2;


  /*******************************************************************
   * COMPUTE SYSTEM PARAMETERS                                       *
   *******************************************************************/

  /* Compute parameters of the system. */
  mTot = params->mTot;
  ASSERT( mTot != 0.0, stat, GENERATEPPNINSPIRALH_EMBAD,
    GENERATEPPNINSPIRALH_MSGEMBAD );
  eta = params->eta;
  ASSERT( eta != 0.0, stat, GENERATEPPNINSPIRALH_EMBAD,
    GENERATEPPNINSPIRALH_MSGEMBAD );
  etaInv = 2.0 / eta;
  mu = eta*mTot;

  sinI = sin( params->inc );
  sin2I = sinI*sinI;
  sin4I = sin2I*sin2I;
  sin5I = sin4I*sinI;

  cosI = cos( params->inc );
  cos2I = cosI*cosI;
  cos4I = cos2I*cos2I;
  cos6I = cos4I*cos2I;

  phiC = params->phi;

  preFac = -2.0*mu*LAL_MRSUN_SI/params->d;
  delta = pow((1-4*eta), 0.5);
  sd = sinI*delta;
  scd = sd*cosI;



  /* *************************************************************
   * COEFFICIENTS IN h_+ & h_x ***********************************
   *************************************************************** */

  /* The variables are annotated by the harmonic, then the polarisation and
     then the frequency power by which they are multiplied. For instance
     'a1Pthree' is the first harmonic, plus polarisation and is multiplied
     by fthree */

  /* First harmonic plus */
  a1Pthree = sd*(5.0 + cos2I)/8.0;

  a1Pfive = - sd*(
                 19.0/64.0 + 5.0/16.0*cos2I - 1.0/192.0*cos4I
                 + eta*(-49.0/96.0 + 1.0/8.0*cos2I + 1.0/96.0*cos4I)
                 );

  a1Psix = + sd*LAL_PI*(5.0 + cos2I)/8.0;

  a1Pseven = - sd*(
                  (1771.0 - 1667.0*cos2I)/5120 + (217*cos4I - cos6I)/9216.0
                  + eta*(
                        681.0/256.0 + (13.0*cos2I - 35*cos4I)/768.0
                        + cos6I/2304.0
                        )
                  + eta*eta*(
                            -(3451.0 + 5.0*cos4I)/9216.0
                            + (673.0*cos2I - cos6I)/3072.0
                            )
                  );

  a1Pmixsix = - sd*(
                   11.0/40.0 + 5.0*log(2.0)/4.0
                   + cos2I*(7.0/40.0 + log(2.0)/4.0)
                   );



  /* Second harmonic plus */
  a2Ptwo = (1.0 + cos2I);

  a2Pfour = - (
              19.0/6.0 + 3.0/2.0*cos2I - 1.0/3.0*cos4I
              + eta*(-19.0/6.0 + 11.0/6.0*cos2I + cos4I)
              );

  a2Pfive = 2.0*LAL_PI*(1.0 + cos2I);

  a2Psix = - (
             11.0/60.0 + 33.0/10.0*cos2I + (29.0*cos4I - 1.0*cos6I)/24.0
             + eta*(
                   353.0/36.0 - 3.0*cos2I - 251.0/72.0*cos4I
                   + 5.0/24.0*cos6I
                   )
             + eta*eta*(-49.0/12.0 + 9.0/2.0*cos2I
             - cos4I*(7.0 + 5.0*cos2I)/24.0)
             );

  a2Pseven = - LAL_PI*(
                      19.0/3.0 + 3.0*cos2I - 2.0/3.0*cos4I
                      + eta*((-16.0 + 14.0*cos2I)/3.0 + 2*cos4I)
                      );

  a2Pmixseven = - (
                  -9.0 + 14.0*cos2I + 7.0*cos4I
                  + eta*(96.0 - 8.0*cos2I - 28.0*cos4I)
                  )/5.0;



  /* Third harmonic plus */
  a3Pthree = -9.0/8.0*sd*(1.0 + cos2I);

  a3Pfive = - sd*(
                 - 657.0/128.0 -45.0/16.0*cos2I + 81.0/128.0*cos4I
                 + eta*(225.0/64.0 - 9.0/8.0*cos2I - 81.0/64.0*cos4I)
                 );

  a3Psix = - sd*LAL_PI*27.0/8.0*(1.0 + cos2I);

  a3Pseven = - sd*(
                  3537.0/1024.0
                  - (22977*cos2I + 15309.0*cos4I - 729.0*cos6I)/5120
                  + eta*(
                        -23829.0 + 5529.0*cos2I
                        + 7749.0*cos4I -729.0*cos6I
                        )/1280.0
                  + eta*eta*(
                            29127.0 - 27267.0*cos2I - 1647.0*cos4I
                            + 2187.0*cos6I
                            )/5120.0
                  );

  a3Pmixsix = - sd*(-189.0/40.0 + 27.0/4.0*log(1.5))*(1 + cos2I);



  /* Fourth harmonic plus */
  a4Pfour = 4.0/3.0*sin2I*(1.0 + cos2I)*(1.0 - 3.0*eta);

  a4Psix = - (
             118.0/15.0 - 16.0/5.0*cos2I - cos4I*(86.0 - 16.0*cos2I)/15.0
             + eta*(
                   -262.0/9.0 + 16.0*cos2I + 166.0/9.0*cos4I
                   - 16.0/3.0*cos6I
                   )
             + eta*eta*(14.0 - 16.0*cos2I + (-10.0*cos4I + 16.0*cos6I)/3.0)
             );

  a4Pseven = + 16.0*LAL_PI/3.0*(1.0 + cos2I)*sin2I*(1.0 - 3.0*eta);

  a4Pmixseven = - sin2I*(1.0 + cos2I)*(
                                      56.0/5.0 - 32.0*log(2.0)/3.0
                                      - eta*(1193.0/30.0 -32.0*log(2.0))
                                      );



  /* Fifth harmonic plus */
  a5Pfive = - sd*(625.0/384.0*sin2I*(1.0 + cos2I)*(1.0 - 2.0*eta));

  a5Pseven = - sd*(
                  (-108125.0 + 40625.0*cos2I + 83125.0*cos4I
                  - 15625.0*cos6I
                  )/9216.0
                  + eta*(8125.0/256.0 - (
                                        40625.0*cos2I + 48125.0*cos4I
                                        - 15625.0*cos6I
                                        )/2304.0
                        )
                  + eta*eta*(
                            (44375.0*cos4I - 119375.0)/9216.0
                            + (40625.0*cos2I - 15625*cos6I)/3072.0)
                            );



  /* Sixth harmonic plus */
  a6Psix = 81.0/40.0*sin4I*(1.0 + cos2I)*(1.0 + 5.0*eta*(eta - 1.0));



  /* Seventh harmonic plus */
  a7Pseven = delta*sin5I*117649.0/46080.0*(1 + cos2I)*(1 + eta*(3.0*eta - 4.0));



  /* First harmonic cross */
  a1Cthree = 3.0/4.0*scd;

  a1Cfive  = - scd*(21.0/32.0 - 5.0/96.0*cos2I + eta*(-23.0 + 5.0*cos2I)/48.0);

  a1Csix   = scd*3.0*LAL_PI/4.0;

  a1Cseven =  - scd*(-913.0/7680.0 + 1891.0/11520.0*cos2I - 7.0/4608.0*cos4I
                 + eta*(1165.0/384.0 - 235.0/576.0*cos2I + 7.0/1152.0*cos4I)
           + eta*eta*(-1301.0/4608.0 + 301.0/2304.0*cos2I - 7.0/1536.0*cos4I));

  a1Cmixsix   = scd*(9.0/20.0 + 3.0*log(2.0)/2.0);



  /* Second harmonic cross */

  a2Ctwo = 2.0*cosI;

  a2Cfour = - cosI*(17.0/3.0 - 4.0/3.0*cos2I + eta*(-13.0/3.0 + 4.0*cos2I));
  a2Cfive = 4*LAL_PI*cosI;

  a2Csix = - cosI*(
                  17.0/15.0 + 113.0/30.0*cos2I - 0.25*cos4I
                  + eta*(143.0/9.0 - 245.0/18.0*cos2I + 5.0/4.0*cos4I)
                  + eta*eta*(-14.0/3.0 + 35.0/6.0*cos2I - 5.0/4.0*cos4I)
                  );

  a2Cseven = - LAL_PI*cosI*(
                           (34.0 - 8.0*cos2I)/3.0
                           - eta*(20.0/3.0 - 8.0*cos2I)
                           );

  a2Cmixseven = - cosI*(2.0 + (-22.0*cos2I + eta*(-154.0 + 94.0*cos2I))/5.0);



  /* Third harmonic cross */
  a3Cthree = - 9.0/4.0*scd;

  a3Cfive = - scd*(
                  -603.0/64.0 + 135.0/64.0*cos2I
                  + eta*(171.0 - 135.0*cos2I)/32.0
                  );

  a3Csix = - scd*27.0/4.0*LAL_PI;

  a3Cseven = - scd*(
                   (12501.0 - 24138.0*cos2I + 1701.0*cos4I)/2560.0
                   + eta*(-19581.0 + 15642.0*cos2I - 1701.0*cos4I)/640.0
                   + eta*eta*(18903.0 - 22806.0*cos2I + 5103.0*cos4I)/2560.0
                   );

  a3Cmixsix = - scd*(189.0/20.0 - 27.0/2.0*log(1.5));



  /* Fourth harmonic cross */
  a4Cfour = cosI*sin2I*8.0/3.0*(1.0 - 3.0*eta);

  a4Csix = - cosI*(
                  44.0/3.0 - 268.0/15.0*cos2I + 16.0/5.0*cos4I
                  + eta*((-476.0 + 620.0*cos2I)/9.0 - 16.0*cos4I)
                  + eta*eta*((68.0 - 116.0*cos2I)/3.0 + 16.0*cos4I)
                  );

  a4Cseven = sin2I*cosI*32.0/3.0*LAL_PI*(1.0 - 3.0*eta);

  a4Cmixseven = - cosI*sin2I*(
                             -112.0/5.0 + 64.0*log(2.0)/3.0
                             + eta*(1193.0/15.0 - 64.0*log(2.0))
                             );



  /* Fifth harmonic cross */
  a5Cfive = - scd*(625.0/192.0*(1.0 - 2.0*eta)*sin2I);

  a5Cseven = - scd*(
                   6875.0/256.0*cos2I
                   - (101875.0 + 21875.0*cos4I)/4608.0
                   + eta*(
                         (66875.0 + 21875.0*cos4I)/1152.0
                         - 44375.0/576.0*cos2I
                         )
                   + eta*eta*(
                             -100625.0/4608.0 + 83125.0/2304.0*cos2I
                             - 21875.0/1536.0*cos4I
                             )
                   );



  /* Sixth harmonic cross */
  a6Csix = cosI*81.0/20.0*sin4I*(1.0 + 5.0*eta*(eta - 1.0));



  /* Seventh harmonic cross */
  a7Cseven = - scd*sin4I*117649.0/23040.0*(1.0 + eta*(3.0*eta - 4.0));



  /* *************************************************************
   * END OF COEFFICIENTS IN h_+ & h_x ****************************
   *************************************************************** */



  /* Compute frequency, phase, and amplitude factors. */
  fFac = 1.0 / ( 4.0*LAL_TWOPI*LAL_MTSUN_SI*mTot );
  dt   = -params->deltaT * eta / ( 5.0*LAL_MTSUN_SI*mTot );
  ASSERT( dt < 0.0, stat, GENERATEPPNINSPIRALH_ETBAD,
    GENERATEPPNINSPIRALH_MSGETBAD );

  f2aFac = LAL_PI*LAL_MTSUN_SI*mTot*fFac;
  ASSERT( params->d != 0.0, stat, GENERATEPPNINSPIRALH_EDBAD,
    GENERATEPPNINSPIRALH_MSGEDBAD );

  /* Compute PN expansion coefficients.
     - Correction to the c5 term berlow in accordance with
       erratum BFIJ, PRD 71 129902
     - c6 does not include log at this stage but is added later
     - p6 = p[6] used in FREQ() macro
  */
  c0 = c[0] =  p[0];
  c1 = c[1] =  p[1];
  c2 = c[2] =  p[2]*( 743.0/2688.0 + eta*11.0/32.0 );
  c3 = c[3] = -p[3]*( 3.0*LAL_PI/10.0 );
  c4 = c[4] =  p[4]*(
                    1855099.0/14450688.0 + eta*56975.0/258048.0 +
                    eta*eta*371.0/2048.0
                    );
  c5 = c[5] =  p[5]*( -7729.0/21504.0 + eta*13.0/256.0 )*LAL_PI;
  c6 = c[6] = -p[6]*(
                    720817631400877.0/288412611379200.0
                    - 107.0*LAL_GAMMA/280.0 - LAL_PI*LAL_PI*53.0/200.0
                    + eta*(
                          - 25302017977.0/4161798144.0
                          + LAL_PI*LAL_PI*451.0/2048.0
                          )
                    + eta*eta*30913.0/1835008.0
                    + eta*eta*eta*235925.0/1769472.0
                    );
  c7 = c[7] = -p[7]*LAL_PI*(
                           377033378.0/867041280.0 + eta*977650.0/2580480.0
                           - eta*eta*283538.0/2580480.0
                           );
  p6 = p[6];

  /* Compute expansion coefficients for series in phi and dy/dx. */
  d0 =  c0;
  d1 =  c1*5.0/4.0;
  d2 =  c2*5.0/3.0;
  d3 =  c3*5.0/2.0;
  d4 =  c4*5.0;
  d5 =  c5*5.0/8.0;
  d6 =  p6*(
           831032450749357.0/57682522275840.0 - LAL_PI*LAL_PI*53.0/40.0
           - 107.0*LAL_GAMMA/56.0
           + eta*(
                 -123292747421.0/4161798144.0 + LAL_PI*LAL_PI*2255.0/2048.0
                 + 385.0/48.0*(-1987.0/3080) -55.0/16.0*(-11831.0/9240.0)
                 )
           + eta*eta*(154565.0/1835008.0 - eta*1179625.0/1769472.0)
           );
  d7 = -c7*5.0/2.0;
  e0 =  c0*3.0;
  e1 =  c1*4.0;
  e2 =  c2*5.0;
  e3 =  c3*6.0;
  e4 =  c4*7.0;
  e5 =  c5*8.0;
  e6 =  c6*9.0;
  e7 =  c7*10.0;

  /* Use Boolean variables to exclude terms that are zero. */
  b0 = b[0] = ( c0 == 0.0 ? 0 : 1 );
  b1 = b[1] = ( c1 == 0.0 ? 0 : 1 );
  b2 = b[2] = ( c2 == 0.0 ? 0 : 1 );
  b3 = b[3] = ( c3 == 0.0 ? 0 : 1 );
  b4 = b[4] = ( c4 == 0.0 ? 0 : 1 );
  b5 = b[5] = ( c5 == 0.0 ? 0 : 1 );
  b6 = b[6] = ( c6 == 0.0 ? 0 : 1 );
  b7 = b[7] = ( c7 == 0.0 ? 0 : 1 );

  /* Find the leading-order frequency term. */
  for ( j = 0; ( j < MAXORDER ) && ( b[j] == 0 ); j++ )
    ;
  if ( j == MAXORDER ) {
    ABORT( stat, GENERATEPPNINSPIRALH_EPBAD,
     GENERATEPPNINSPIRALH_MSGEPBAD );
  }


  /*******************************************************************
   * COMPUTE START TIME                                              *
   *******************************************************************/

  /* First, find the normalized start frequency, and the best guess as
     to the start times from each term.  We require the
     frequency to be increasing. */
  yStart =  2.0/(REAL4)(harmonics)*params->fStartIn / fFac;

  if ( params->fStopIn == 0.0 )
    yMax = 1.0/(LAL_PI*pow(6.0, 1.5)*mTot*LAL_MTSUN_SI) / fFac;
  else if( params->fStopIn < 0.0 )
    yMax = -1.0*params->fStopIn / fFac;
  else {
    ASSERT( fabs( params->fStopIn ) > params->fStartIn, stat,
      GENERATEPPNINSPIRALH_EFBAD, GENERATEPPNINSPIRALH_MSGEFBAD );
    yMax = fabs( params->fStopIn ) / fFac;
  }

  if ( ( c[j]*fFac < 0.0 ) || ( yStart < 0.0 ) || ( yMax < 0.0 ) ) {
     ABORT( stat, GENERATEPPNINSPIRALH_EPBAD,
     GENERATEPPNINSPIRALH_MSGEPBAD );
  }

  xStart = pow( yStart/c[j], 1.0/( j + 3.0 ) );
  xMax = LAL_SQRT2;

  /* The above is exact if the leading-order term is the only one in
   the expansion.  Check to see if there are any other terms. */
  for ( i = j + 1; ( i < MAXORDER ) && ( b[i] == 0 ); i++ )
  ;
  if ( i < MAXORDER )
  {
    /* There are other terms, so we have to use bisection to find the
       start time. */
    REAL4 xLow, xHigh; /* ultimately these will bracket xStart */
    REAL4 yLow, yHigh; /* the normalized frequency at these times */

    if ( xStart > 0.39*xMax )
      xStart = 0.39*xMax;

    /* If we are ignoring PN breakdown, adjust xMax (so that it won't
       interfere with the start time search) and tStop. */
    if ( params->fStopIn < 0.0 )
    {
      xMax = LAL_REAL4_MAX;
      tStop = 0.0;
    }

    /* If our frequency is too high, step backwards and/or forwards
       until we have bracketed the correct frequency. */
    xLow = xHigh = xStart;
    FREQ( yHigh, xStart);
    yLow = yHigh;
    while ( yLow > yStart )
    {
      xHigh = xLow;
      yHigh = yLow;
      xLow *= 0.95;
      FREQ( yLow, xLow );
    }

    while ( yHigh < yStart )
    {
      xLow = xHigh;
      yLow = yHigh;
      xHigh *= 1.05;
      FREQ( yHigh, xHigh );
      /* Check for PN breakdown. */
      if ( ( yHigh < yLow ) || ( xHigh > xMax ) )
      {
        ABORT( stat, GENERATEPPNINSPIRALH_EFBAD,
         GENERATEPPNINSPIRALH_MSGEFBAD );
      }
    }

    /* We may have gotten lucky and nailed the frequency right on.
         Otherwise, find xStart by root bisection. */
    if ( yLow == yStart )
      xStart = xLow;
    else if ( yHigh == yStart )
      xStart = xHigh;
    else
    {
      SFindRootIn in;
      FreqDiffParamStruc par;
      in.xmax = xHigh;
      in.xmin = xLow;
      in.xacc = ACCURACY;
      in.function = FreqDiff;
      par.c = c;
      par.b = b;
      par.p6 = p6;
      par.y0 = yStart;

      TRY( LALSBisectionFindRoot( stat->statusPtr, &(xStart), &in,
          (void *)( &par ) ), stat );
    }

  }

  /* If we are ignoring PN breakdown, adjust xMax and tStop, if they
     haven't been adjusted already. */
  else if ( params->fStopIn < 0.0 )
  {
    xMax = LAL_REAL4_MAX;
    tStop = 0.0;
  }

  /* Compute initial dimensionless time, record actual initial
     frequency (in case it is different), and record dimensional
     time-to-coalescence. */
  t0 = pow( xStart, -8.0 );
  FREQ( yStart, xStart );
  if ( yStart >= yMax )
  {
    ABORT( stat, GENERATEPPNINSPIRALH_EFBAD,
                   GENERATEPPNINSPIRALH_MSGEFBAD );
  }
  params->fStart = yStart*fFac;
  params->tc = t0 * ( 5.0*LAL_MTSUN_SI*mTot ) / eta;



  /*******************************************************************
   * GENERATE WAVEFORM                                               *
   *******************************************************************/

  /* Set up data pointers and storage. */
  here = head = (PPNInspiralBuffer *)
  LALMalloc( sizeof(PPNInspiralBuffer) );
  if ( !here )
  {
    ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
     GENERATEPPNINSPIRALH_MSGEMEM );
  }
  here->next = NULL;
  h = here->h;
  a = here->a;
  f = here->f;

  phi = here->phi;
  nMax = (UINT4)( -1 );
  if ( params->lengthIn > 0 )
    nMax = params->lengthIn;
  nNext = BUFFSIZE;
  if ( nNext > nMax )
    nNext = nMax;


  /* Start integrating!  Inner loop exits each time a new buffer is
     required.  Outer loop has no explicit test; when a termination
     condition is met, we jump directly from the inner loop using a
     goto statement.  All goto statements jump to the terminate: label
     at the end of the outer loop. */
  n = 0;
  t = t0;
  dyMax = 0.0;
  y = yOld = 0.0;
  x = xStart;

  while ( 1 )
  {
    while ( n < nNext )
    {
      REAL4 f2a; /* value inside 2/3 power in amplitude functions */
      REAL4 phase = 0.0; /* wave phase excluding overall constants */
      REAL4 dydx2 = 0.0; /* dy/dx divided by x^2 */

      /* Check if we're still in a valid PN regime. */
      if ( x > xMax )
      {
        params->termCode = GENERATEPPNINSPIRALH_EPNFAIL;
        params->termDescription = GENERATEPPNINSPIRALH_MSGEPNFAIL;
        goto terminate;
      }

     /* Compute the normalized frequency.  This also computes the
        variables x2, x3, x4, x5, x6 and x7 which are used later. */
     FREQ( y, x );

     if ( y > yMax )
     {
       params->termCode = GENERATEPPNINSPIRALH_EFSTOP;
       params->termDescription = GENERATEPPNINSPIRALH_MSGEFSTOP;
       goto terminate;
     }


     /* Check that frequency is still increasing. */
     if ( b0 )
       dydx2 += e0;
     if ( b1 )
       dydx2 += e1*x;
     if ( b2 )
       dydx2 += e2*x2;
     if ( b3 )
       dydx2 += e3*x3;
     if ( b4 )
       dydx2 += e4*x4;
     if ( b5 )
       dydx2 += e5*x5;
     if ( b6 )
       dydx2 += (e6 + (856.0/2240.0*(2.0 + 9.0*log(2.0*x))))*x6;
     if ( b7 )
       dydx2 += e7*x7;

    if ( dydx2 < 0.0 )
    {
       params->termCode = GENERATEPPNINSPIRALH_EFNOTMON;
       params->termDescription = GENERATEPPNINSPIRALH_MSGEFNOTMON;
       goto terminate;
    }

     if ( y - yOld > dyMax )
       dyMax = y - yOld;
     *(f++) = fFac*y;

     /* Compute the phase. */
     if ( b0 )
       phase += d0;
     if ( b1 )
       phase += d1*x;
     if ( b2 )
       phase += d2*x2;
     if ( b3 )
       phase += d3*x3;
     if ( b4 )
       phase += d4*x4;
     if ( b5 )
       phase += d5*log(t)*x5;
     if ( b6 )
       phase += (d6 - 8.0*107.0*log(2.0*x)/448.0)*x6;
     if ( b7 )
       phase += d7*x7;
     /* etaInv absorbs the factor of 2! */
     phase *= t*x3*etaInv;
     *(phi++) = phiC - phase;

     /* Compute hplus and hcross */
     f2a = pow(f2aFac*y, TWOTHIRDS);

     /* powers of frequency */
     fthree = pow(f2a, 1.5);
     ffour  = pow(f2a, 2.0);
     ffive  = pow(f2a, 2.5);
     fsix   = pow(f2a, 3.0);
     fseven = pow(f2a, 3.5);

     /* PLUS */
     a1 =   q[1]*a1Pthree*fthree + q[3]*a1Pfive*ffive + q[4]*a1Psix*fsix
          + q[5]*a1Pseven*fseven;
     a1mix = q[4]*a1Pmixsix*fsix;

     a2 =   q[0]*a2Ptwo*f2a + q[2]*a2Pfour*ffour + q[3]*a2Pfive*ffive
          + q[4]*a2Psix*fsix + q[5]*a2Pseven*fseven;
     a2mix = q[5]*a2Pmixseven*fseven;

     a3 =   q[1]*a3Pthree*fthree + q[3]*a3Pfive*ffive + q[4]*a3Psix*fsix
          + q[5]*a3Pseven*fseven;
     a3mix = q[4]*a3Pmixsix*fsix;

     a4 = q[2]*a4Pfour*ffour + q[4]*a4Psix*fsix + q[5]*a4Pseven*fseven;
     a4mix = q[5]*a4Pmixseven*fseven;

     a5 = q[3]*a5Pfive*ffive + q[5]*a5Pseven*fseven;

     a6 = q[4]*a6Psix*fsix;

     a7 = q[5]*a7Pseven*fseven;

     /* Store first three harmonics for filtering */
     *(a++) = preFac*a1Pthree*fthree;
     *(a++) = preFac*a2Ptwo*f2a;
     *(a++) = preFac*a3Pthree*fthree;

     *(h++) = preFac*(
                     a1*cos(1.0/2.0*(phiC - phase))
                     + a2*cos(2.0/2.0*(phiC - phase))
                     + a3*cos(3.0/2.0*(phiC - phase))
                     + a4*cos(4.0/2.0*(phiC - phase))
                     + a5*cos(5.0/2.0*(phiC - phase))
                     + a6*cos(6.0/2.0*(phiC - phase))
                     + a7*cos(7.0/2.0*(phiC - phase))
                     + a1mix*sin(1.0/2.0*(phiC - phase))
                     + a2mix*sin(2.0/2.0*(phiC - phase))
                     + a3mix*sin(3.0/2.0*(phiC - phase))
                     + a4mix*sin(4.0/2.0*(phiC - phase))
                     );


     /* CROSS */
     a1 =   q[1]*a1Cthree*fthree + q[3]*a1Cfive*ffive + q[4]*a1Csix*fsix
          + q[5]*a1Cseven*fseven;
     a1mix = q[4]*a1Cmixsix*fsix;

     a2 =   q[0]*a2Ctwo*f2a + q[2]*a2Cfour*ffour + q[3]*a2Cfive*ffive
          + q[4]*a2Csix*fsix + q[5]*a2Cseven*fseven;
     a2mix = q[5]*a2Cmixseven*fseven;

     a3 =   q[1]*a3Cthree*fthree + q[3]*a3Cfive*ffive + q[4]*a3Csix*fsix
          + q[5]*a3Cseven*fseven;
     a3mix = q[4]*a3Cmixsix*fsix;

     a4 = q[2]*a4Cfour*ffour + q[4]*a4Csix*fsix + q[5]*a4Cseven*fseven;
     a4mix = q[5]*a4Cmixseven*fseven;

     a5 = q[3]*a5Cfive*ffive + q[5]*a5Cseven*fseven;

     a6 = q[4]*a6Csix*fsix;

     a7 = q[5]*a7Cseven*fseven;

     *(h++) = preFac*(
                     a1*sin(1.0/2.0*(phiC - phase))
                     + a2*sin(2.0/2.0*(phiC - phase))
                     + a3*sin(3.0/2.0*(phiC - phase))
                     + a4*sin(4.0/2.0*(phiC - phase))
                     + a5*sin(5.0/2.0*(phiC - phase))
                     + a6*sin(6.0/2.0*(phiC - phase))
                     + a7*sin(7.0/2.0*(phiC - phase))
                     + a1mix*cos(1.0/2.0*(phiC - phase))
                     + a2mix*cos(2.0/2.0*(phiC - phase))
                     + a3mix*cos(3.0/2.0*(phiC - phase))
                     + a4mix*cos(4.0/2.0*(phiC - phase))
                     );

     n++;
     t = t0 + n*dt;
     yOld = y;
     if ( t <= tStop ) {
       params->termCode = GENERATEPPNINSPIRALH_ERTOOSMALL;
       params->termDescription = GENERATEPPNINSPIRALH_MSGERTOOSMALL;
       goto terminate;
     }
     x = pow( t, -0.125 );

   }

   /* We've either filled the buffer or we've exceeded the maximum
      length.  If the latter, we're done! */
   if ( n >= nMax ) {
     params->termCode = GENERATEPPNINSPIRALH_ELENGTH;
     params->termDescription = GENERATEPPNINSPIRALH_MSGELENGTH;

   }


   /* Otherwise, allocate the next buffer. */
   here->next =
      (PPNInspiralBuffer *)LALMalloc( sizeof(PPNInspiralBuffer) );
   here = here->next;
   if ( !here ) {
     FREELIST( head );
     ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
       GENERATEPPNINSPIRALH_MSGEMEM );
   }

   here->next = NULL;
   h = here->h;
   a = here->a;
   f = here->f;
   phi = here->phi;
   nNext += BUFFSIZE;
   if ( nNext > nMax )
     nNext = nMax;


 }


  /*******************************************************************
   * CLEANUP                                                         *
   *******************************************************************/

  /* The above loop only exits by triggering one of the termination
     conditions, which jumps to the following point for cleanup and
     return. */
 terminate:

  /* First, set remaining output parameter fields. */
  params->dfdt = dyMax*fFac*params->deltaT;
  params->fStop = yOld*fFac;
  params->length = n;

  /* Allocate the output structures. */
  if ( ( output->h = (REAL4TimeVectorSeries *)
       LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL )
  {
    FREELIST( head );
    ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
     GENERATEPPNINSPIRALH_MSGEMEM );
  }
  memset( output->h, 0, sizeof(REAL4TimeVectorSeries) );

  if ( ( output->a = (REAL4TimeVectorSeries *)
       LALMalloc( sizeof(REAL4TimeVectorSeries) ) ) == NULL )
  {
    FREELIST( head );
    LALFree( output->h ); output->h = NULL;
    ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
     GENERATEPPNINSPIRALH_MSGEMEM );
  }
  memset( output->a, 0, sizeof(REAL4TimeVectorSeries) );

  if ( ( output->f = (REAL4TimeSeries *)
       LALMalloc( sizeof(REAL4TimeSeries) ) ) == NULL )
  {
    FREELIST( head );
    LALFree( output->h ); output->h = NULL;
    LALFree( output->a ); output->a = NULL;
    ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
     GENERATEPPNINSPIRALH_MSGEMEM );
  }
  memset( output->f, 0, sizeof(REAL4TimeSeries) );

  if ( ( output->phi = (REAL8TimeSeries *)
       LALMalloc( sizeof(REAL8TimeSeries) ) ) == NULL )
  {
    FREELIST( head );
    LALFree( output->h ); output->h = NULL;
    LALFree( output->a ); output->a = NULL;
    LALFree( output->f ); output->f = NULL;
    ABORT( stat, GENERATEPPNINSPIRALH_EMEM,
     GENERATEPPNINSPIRALH_MSGEMEM );
  }
  memset( output->phi, 0, sizeof(REAL8TimeSeries) );

  /* Allocate the output data fields. */
  {
    CreateVectorSequenceIn in;
    in.length = n;
    in.vectorLength = 2;
    LALSCreateVectorSequence( stat->statusPtr, &( output->h->data ), &in );
    BEGINFAIL( stat )
    {
      FREELIST( head );
      LALFree( output->h );   output->h = NULL;
      LALFree( output->a );   output->a = NULL;
      LALFree( output->f );   output->f = NULL;
      LALFree( output->phi ); output->phi = NULL;
    } ENDFAIL( stat );
    /* Adjust in for waveform.a which stores the first 3 harmonics */
    in.vectorLength = 3;
    LALSCreateVectorSequence( stat->statusPtr, &( output->a->data ), &in );
    BEGINFAIL( stat )
    {
      FREELIST( head );
      TRY( LALSDestroyVectorSequence( stat->statusPtr, &( output->h->data ) ),
           stat );
      LALFree( output->h );   output->h = NULL;
      LALFree( output->a );   output->a = NULL;
      LALFree( output->f );   output->f = NULL;
      LALFree( output->phi ); output->phi = NULL;
    } ENDFAIL( stat );
    LALSCreateVector( stat->statusPtr, &( output->f->data ), n );
    BEGINFAIL( stat )
    {
      TRY( LALSDestroyVectorSequence( stat->statusPtr, &( output->h->data ) ),
           stat );
      TRY( LALSDestroyVectorSequence( stat->statusPtr, &( output->a->data ) ),
           stat );
      FREELIST( head );
      LALFree( output->h );   output->h = NULL;
      LALFree( output->a );   output->a = NULL;
      LALFree( output->f );   output->f = NULL;
      LALFree( output->phi ); output->phi = NULL;
    } ENDFAIL( stat );
    LALDCreateVector( stat->statusPtr, &( output->phi->data ), n );
    BEGINFAIL( stat )
    {
      TRY( LALSDestroyVectorSequence( stat->statusPtr, &( output->h->data ) ),
           stat );
      TRY( LALSDestroyVectorSequence( stat->statusPtr, &( output->a->data ) ),
           stat );
      TRY( LALSDestroyVector( stat->statusPtr, &( output->f->data ) ),
           stat );
      FREELIST( head );
      LALFree( output->h );   output->h = NULL;
      LALFree( output->a );   output->a = NULL;
      LALFree( output->f );   output->f = NULL;
      LALFree( output->phi ); output->phi = NULL;
    } ENDFAIL( stat );
  }


  /* Structures have been successfully allocated; now fill them.  We
     deallocate the list as we go along. */
  output->position = params->position;
  output->psi = params->psi;
  output->h->epoch = output->a->epoch = output->f->epoch = output->phi->epoch
    = params->epoch;
  output->h->deltaT = output->a->deltaT = output->f->deltaT
    = output->phi->deltaT = params->deltaT;
  output->h->sampleUnits = lalStrainUnit;
  output->a->sampleUnits = lalStrainUnit;
  output->f->sampleUnits = lalHertzUnit;
  output->phi->sampleUnits = lalDimensionlessUnit;
  snprintf( output->h->name, LALNameLength,
               "PPN inspiral waveform polarisations" );
  snprintf( output->a->name, LALNameLength,
               "First three harmonics a1_+, a2_+ & a3_+" );
  snprintf( output->f->name, LALNameLength, "PPN inspiral frequency" );
  snprintf( output->phi->name, LALNameLength, "PPN inspiral phase" );
  h = output->h->data->data;
  a = output->a->data->data;
  f = output->f->data->data;
  phi = output->phi->data->data;
  here = head;
  while ( here && ( n > 0 ) )
  {
    PPNInspiralBuffer *last = here;
    UINT4 nCopy = BUFFSIZE;
    if ( nCopy > n )
      nCopy = n;
    memcpy( h, here->h, 2*nCopy*sizeof(REAL4) );
    memcpy( a, here->a, 3*nCopy*sizeof(REAL4) );
    memcpy( f, here->f, nCopy*sizeof(REAL4) );
    memcpy( phi, here->phi, nCopy*sizeof(REAL8) );
    h += 2*nCopy;
    a += 3*nCopy;
    f += nCopy;
    phi += nCopy;
    n -= nCopy;
    here = here->next;
    LALFree( last );
  }

  /* This shouldn't happen, but free any extra buffers in the
     list. */
  FREELIST( here );

  /* Everything's been stored and cleaned up, so there's nothing left
     to do but quit! */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
