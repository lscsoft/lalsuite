#if 0  /* autodoc block */

<lalVerbatim file="RandomCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{Random.c}}
\label{ss:Random.c}

Functions for generating random numbers.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{RandomCP}
\index{\texttt{LALCreateRandomParams()}}
\index{\texttt{LALDestroyRandomParams()}}
\index{\texttt{LALUniformDeviate()}}
\index{\texttt{LALNormalDeviates()}}

\subsubsection*{Description}

The routines \verb+LALCreateRandomParams()+ and \verb+LALDestroyRandomParams()+
create and destroy a parameter structure for the generation of random
variables.  The creation routine requires a random number seed \verb+seed+.
If the seed is zero then a seed is generated using the current time.

The routine \verb+LALUniformDeviate()+ returns a single random deviate
distributed uniformly between zero and unity.  

The routine \verb+LALNormalDeviates()+ fills a vector with normal (Gaussian)
deviates with zero mean and unit variance.

\subsubsection*{Operating Instructions}

\begin{verbatim}
static LALStatus     status;
static RandomParams *params;
static REAL4Vector  *vector;
UINT4 i;
INT4 seed = 0;

LALCreateVector( &status, &vector, 9999 );
LALCreateRandomParams( &status, &params, seed );

/* fill vector with uniform deviates */
for ( i = 0; i < vector->length; ++i )
{
  LALUniformDeviate( &status, vector->data + i, params );
}

/* fill vector with normal deviates */
LALNormalDeviates( &status, vector, params );

LALDestroyRandomParams( &status, &params );
LALDestroyVector( &status, &vector );
\end{verbatim}

\subsubsection*{Algorithm}

This is an implementation of the random number generators \verb+ran1+ and
\verb+gasdev+ described in Numerical Recipes~\cite{ptvf:1992}.

\subsubsection*{Uses}

\subsubsection*{Notes}
\vfill{\footnotesize\input{RandomCV}}

</lalLaTeX>

#endif /* autodoc block */


#include <time.h>
#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/Random.h>

NRCSID (RANDOMC, "$Id$");

static const INT4 a = 16807;
static const INT4 m = 2147483647;
static const INT4 q = 127773;
static const INT4 r = 2836;

static const REAL4 eps = 1.2e-7;

static
INT4
BasicRandom (INT4 i)
{
  INT4 k;
  k = i/q;
  i = a*(i - k*q) - r*k;
  if (i < 0)
    i += m;
  return i;
}

/* <lalVerbatim file="RandomCP"> */
void
LALCreateRandomParams (
    LALStatus     *status,
    RandomParams **params,
    INT4           seed
    )
{ /* </lalVerbatim> */
  INT4 n;

  INITSTATUS (status, "LALCreateRandomParams", RANDOMC);

  ASSERT (params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (!*params, status, RANDOMH_ENNUL, RANDOMH_MSGENNUL);

  *params = (RandomParams *) LALMalloc (sizeof(RandomParams));
  ASSERT (*params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);

  while (seed == 0)
  {
    seed = time (NULL);
  }

  if (seed < 0)
  {
    seed = -seed;
  }

  (*params)->i = seed;

  for (n = 0; n < 8; ++n)
  {
    (*params)->i = BasicRandom ((*params)->i);
  }

  for (n = 0; n < (int)(sizeof((*params)->v)/sizeof(INT4)); ++n)
  {
    (*params)->v[n] = (*params)->i = BasicRandom ((*params)->i);
  }

  (*params)->y = (*params)->v[0];

  RETURN (status);
}

/* <lalVerbatim file="RandomCP"> */
void
LALDestroyRandomParams (
    LALStatus     *status,
    RandomParams **params
    )
{ /* </lalVerbatim> */
  INITSTATUS (status, "LALDestroyRandomParams", RANDOMC);

  ASSERT (params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (*params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);

  LALFree (*params);
  *params = NULL;

  RETURN (status);
}


/* <lalVerbatim file="RandomCP"> */
void
LALUniformDeviate (
    LALStatus    *status,
    REAL4        *deviate,
    RandomParams *params
    )
{ /* </lalVerbatim> */
  INT4 ndiv;
  INT4 n;

  INITSTATUS (status, "LALUniformDeviate", RANDOMC);

  ASSERT (deviate, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);

  ndiv = 1 + (m - 1)/(sizeof(params->v)/sizeof(INT4));
  n    = params->y/ndiv;

  params->y = params->v[n];
  params->v[n] = params->i = BasicRandom (params->i);

  *deviate = params->y/(REAL4)m;
  if (*deviate > 1 - eps)
  {
    *deviate = 1 - eps;
  }

  RETURN (status);
}


/* <lalVerbatim file="RandomCP"> */
void
LALNormalDeviates (
    LALStatus    *status,
    REAL4Vector  *deviates,
    RandomParams *params
    )
{ /* </lalVerbatim> */
  REAL4 *data;
  INT4   half;

  INITSTATUS (status, "LALNormalDeviates", RANDOMC);
  ATTATCHSTATUSPTR (status);

  ASSERT (params, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (deviates, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (deviates->data, status, RANDOMH_ENULL, RANDOMH_MSGENULL);
  ASSERT (deviates->length > 0, status, RANDOMH_ESIZE, RANDOMH_MSGESIZE);

  data = deviates->data;
  half = deviates->length/2;

  while (half-- > 0)
  {
    REAL4 u;
    REAL4 v;
    REAL4 x;
    REAL4 y;
    REAL4 rsq;
    REAL4 fac;

    do {
      LALUniformDeviate (status->statusPtr, &u, params);
      CHECKSTATUSPTR (status);
      LALUniformDeviate (status->statusPtr, &v, params);
      CHECKSTATUSPTR (status);
      x   = 2*u - 1;
      y   = 2*v - 1;
      rsq = x*x + y*y;
    }
    while (rsq >= 1 || rsq == 0);

    fac     = sqrt(-2*log(rsq)/rsq);
    *data++ = fac*x;
    *data++ = fac*y;
  }

  /* do it again if there is an odd amount of data */
  if (deviates->length % 2)
  {
    REAL4 u;
    REAL4 v;
    REAL4 x;
    REAL4 y;
    REAL4 rsq;
    REAL4 fac;

    do {
      LALUniformDeviate (status->statusPtr, &u, params);
      CHECKSTATUSPTR (status);
      LALUniformDeviate (status->statusPtr, &v, params);
      CHECKSTATUSPTR (status);
      x   = 2*u - 1;
      y   = 2*v - 1;
      rsq = x*x + y*y;
    }
    while (rsq >= 1 || rsq == 0);

    fac   = sqrt(-2*log(rsq)/rsq);
    *data = fac*x;
    /* throw away y */
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
