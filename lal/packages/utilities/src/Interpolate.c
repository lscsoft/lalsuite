#if 0  /* autodoc block */

<lalVerbatim file="InterpolateCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{Interpolate.c}}
\label{ss:Interpolate.c}

Functions for generating random numbers.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{InterpolateCP}
\idx{LALSPolynomialInterpolation()}
\idx{LALDPolynomialInterpolation()}

\subsubsection*{Description}

The routine \verb+LALSPolynomialInterpolation()+ computes the interpolated $y$
value \verb+output+ at the $x$ value \verb+target+ by fitting a polynomial of
order \verb+params.n-1+ to the data.  The result \verb+output+ is of type
\verb+SInterpolateOut+, which contains the value \verb+output.y+ as well as
an estimate of the error \verb+output.dy+.  The routine
\verb+LALDPolynomialInterpolation()+ is the same but for double precision.


\subsubsection*{Operating Instructions}

The following program fits a fourth-order polynomial to the five data points
$\{(0,0),(1,1),(2,3),(3,4),(4,3)\}$, and interpolates the value at $x=2.4$.

\begin{verbatim}
#include <lal/LALStdlib.h>
#include <lal/Interpolate.h>

int main ()
{
  enum { ArraySize = 5 };
  static LALStatus status;
  REAL4            x[ArraySize] = {0,1,2,3,4};
  REAL4            y[ArraySize] = {0,1,3,4,3};
  REAL4            target       = 2.4;
  SInterpolatePar  intpar       = {ArraySize, x, y};
  SInterpolateOut  intout;

  LALSPolynomialInterpolation( &status, &intout, target, &intpar );

  return 0;
}
\end{verbatim}

\subsubsection*{Algorithm}

This is an implementation of the Neville algroithm, see \verb+polint+ in
Numerical Recipes~\cite{ptvf:1992}.

\subsubsection*{Uses}

\subsubsection*{Notes}
\vfill{\footnotesize\input{InterpolateCV}}

</lalLaTeX>

#endif /* autodoc block */


#include <math.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/Interpolate.h>

NRCSID (INTERPOLATEC, "$Id$");

/* <lalVerbatim file="InterpolateCP"> */
void
LALSPolynomialInterpolation (
    LALStatus       *status,
    SInterpolateOut *output,
    REAL4            target,
    SInterpolatePar *params
    )
{ /* </lalVerbatim> */
  REAL4 *dn;   /* difference in a step down */
  REAL4 *up;   /* difference in a step up   */
  REAL4  diff;
  UINT4  near;
  UINT4  order;
  UINT4  n;
  UINT4  i;

  INITSTATUS (status, "LALSPolynomialInterpolation", INTERPOLATEC);

  ASSERT (output, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params->x, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params->y, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);

  n = params->n;
  ASSERT (n > 1, status, INTERPOLATEH_ESIZE, INTERPOLATEH_MSGESIZE);

  dn = (REAL4 *) LALMalloc (n*sizeof(REAL4));
  ASSERT (dn, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);

  up = (REAL4 *) LALMalloc (n*sizeof(REAL4));
  ASSERT (up, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);


  /*
   *
   * Initialize dn[] and up[] and find element of domain nearest the target.
   *
   */


  memcpy (dn, params->y, n*sizeof(REAL4));
  memcpy (up, params->y, n*sizeof(REAL4));
  diff = target - params->x[near = 0];
  for (i = 1; diff > 0 && i < n; ++i)
  {
    REAL4 tmp = target - params->x[i];
    diff = (fabs(tmp) < fabs(diff) ? near = i, tmp : diff);
  }
  output->y = params->y[near];


  /*
   *
   * Recompute dn[] and up[] for each polynomial order.
   *
   */


  for (order = 1; order < n; ++order)
  {
    UINT4 imax = n - order;
    for (i = 0; i < imax; ++i)
    {
      REAL4 xdn = params->x[i];
      REAL4 xup = params->x[i + order];
      REAL4 den = xdn - xup;
      REAL4 fac;
      ASSERT (den != 0, status, INTERPOLATEH_EZERO, INTERPOLATEH_MSGEZERO);
      fac   = (dn[i + 1] - up[i])/den;
      dn[i] = fac*(xdn - target);
      up[i] = fac*(xup - target);
    }

    /* go down unless impossible */
    output->y += (near < imax ? dn[near] : up[--near]);
  }

  output->dy = fabs(dn[0]) < fabs(up[0]) ? fabs(dn[0]) : fabs(up[0]);

  LALFree (dn);
  LALFree (up);
  RETURN  (status);
}


/* <lalVerbatim file="InterpolateCP"> */
void
LALDPolynomialInterpolation (
    LALStatus       *status,
    DInterpolateOut *output,
    REAL8            target,
    DInterpolatePar *params
    )
{ /* </lalVerbatim> */
  REAL8 *dn;   /* difference in a step down */
  REAL8 *up;   /* difference in a step up   */
  REAL8  diff;
  UINT4  near;
  UINT4  order;
  UINT4  n;
  UINT4  i;

  INITSTATUS (status, "LALDPolynomialInterpolation", INTERPOLATEC);

  ASSERT (output, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params->x, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);
  ASSERT (params->y, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);

  n = params->n;
  ASSERT (n > 1, status, INTERPOLATEH_ESIZE, INTERPOLATEH_MSGESIZE);

  dn = (REAL8 *) LALMalloc (n*sizeof(REAL8));
  ASSERT (dn, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);

  up = (REAL8 *) LALMalloc (n*sizeof(REAL8));
  ASSERT (up, status, INTERPOLATEH_ENULL, INTERPOLATEH_MSGENULL);


  /*
   *
   * Initialize dn[] and up[] and find element of domain nearest the target.
   *
   */


  memcpy (dn, params->y, n*sizeof(REAL8));
  memcpy (up, params->y, n*sizeof(REAL8));
  diff = target - params->x[near = 0];
  for (i = 1; diff > 0 && i < n; ++i)
  {
    REAL8 tmp = target - params->x[i];
    diff = (fabs(tmp) < fabs(diff) ? near = i, tmp : diff);
  }
  output->y = params->y[near];


  /*
   *
   * Recompute dn[] and up[] for each polynomial order.
   *
   */


  for (order = 1; order < n; ++order)
  {
    UINT4 imax = n - order;
    for (i = 0; i < imax; ++i)
    {
      REAL8 xdn = params->x[i];
      REAL8 xup = params->x[i + order];
      REAL8 den = xdn - xup;
      REAL8 fac;
      ASSERT (den != 0, status, INTERPOLATEH_EZERO, INTERPOLATEH_MSGEZERO);
      fac   = (dn[i + 1] - up[i])/den;
      dn[i] = fac*(xdn - target);
      up[i] = fac*(xup - target);
    }

    /* go down unless impossible */
    output->y += (near < imax ? dn[near] : up[--near]);
  }

  output->dy = fabs(dn[0]) < fabs(up[0]) ? fabs(dn[0]) : fabs(up[0]);

  LALFree (dn);
  LALFree (up);
  RETURN  (status);
}
