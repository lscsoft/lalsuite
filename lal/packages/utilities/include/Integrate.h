#if 0 /* autodoc block */

<lalVerbatim file="IntegrateHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{Integrate.h}}
\label{s:Integrate.h}

Integrates a function.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Integrate.h>
\end{verbatim}

\noindent This header covers the routines for integrating a function.

</lalLaTeX>

#endif /* autodoc block */


#ifndef _INTEGRATE_H
#define _INTEGRATE_H

#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID (INTEGRATEH, "$Id$");

#if 0 /* autodoc block */

<lalLaTeX>
\subsection*{Error conditions}
\input{IntegrateHErrTab}
</lalLaTeX>

<lalErrTable file="IntegrateHErrTab">

#endif /* autodoc block */

#define INTEGRATEH_ENULL 1
#define INTEGRATEH_ETYPE 2
#define INTEGRATEH_EIDOM 4
#define INTEGRATEH_EMXIT 8

#define INTEGRATEH_MSGENULL "Null pointer"
#define INTEGRATEH_MSGETYPE "Unknown integral type"
#define INTEGRATEH_MSGEIDOM "Invalid domain"
#define INTEGRATEH_MSGEMXIT "Maximum iterations exceeded"

#if 0 /* autodoc block */

</lalErrTable>

<lalLaTeX>

\subsection*{Structures}

\begin{verbatim}
typedef enum
{
  ClosedInterval,
  OpenInterval,
  SingularLowerLimit,
  SingularUpperLimit,
  InfiniteDomainPow,
  InfiniteDomainExp
}
IntegralType;

typedef struct
tagSIntegrateIn
{
  void (*function)(LALStatus *s, REAL4 *y, REAL4 x, void *p);
  REAL4         xmax;
  REAL4         xmin;
  IntegralType  type;
}
SIntegrateIn;

typedef struct
tagDIntegrateIn
{
  void (*function)(LALStatus *s, REAL8 *y, REAL8 x, void *p);
  REAL8         xmax;
  REAL8         xmin;
  IntegralType  type;
}
DIntegrateIn;
\end{verbatim}

These are input structures to the integration routines.  The fields are:
\begin{description}
\item[\texttt{function}] The function to integrate.
\item[\texttt{xmax}] The maximum value of the domain of integration.
\item[\texttt{xmax}] The minimum value of the domain of integration.
\item[\texttt{type}] The type of integration.  This is an enumerated type
  which can take on one of the following values:

  \noindent
  \begin{tabular}{|l|l|}
  \hline
  Enumeration constant & Meaning \\
  \hline
  \verb+ClosedInterval+     & Evaluate integral on a closed interval \\
  \verb+OpenInterval+       & Evaluate integral on an open interval \\
  \verb+SingularLowerLimit+ & Integrate an inverse square-root singularity at lower limit \\
  \verb+SingularUpperLimit+ & Integrate an inverse square-root singularity at upper limit \\
  \verb+InfiniteDomainPow+  & Integrate an infinite domain with power-law falloff \\
  \verb+InfiniteDomainExp+  & Integrate an infinite domain with exponential falloff \\
  \hline
  \end{tabular}

  \noindent
  The types of integration are the following: I.\@ \verb+ClosedInterval+
  indicates that the integral should be computed on equal-spaced domain
  intervals including the boundary.  II.\@ \verb+OpenInterval+ indicates that
  the integral should be computed on intervals of the domain not including the
  boundary.  III.\@ \verb+SingularLowerLimit+ indicates that the integral
  should be evaluated on an open interval with a transformation so that a
  inverse-square-root singularity at the lower limit can be integrated.  IV.\@
  \verb+SingularUpperLimit+ is the same as above but for a singularity at the
  upper limit.  V.\@ \verb+InfiniteDomainPow+ indicates that the integral
  should be evaluated over an semi-infinite domain---appropriate when both
  limits have the same sign (though one is very large) and when the integrand
  vanishes faster than $x^{-1}$ at infinity.  VI.\@ \verb+InfiniteDomainExp+
  indicates that the integral should be evaluated over an infinite domain
  starting at \verb+xmin+ and going to infinity (\verb+xmax+ is ignored)---the
  integrand should vanish exponentially for large $x$.
\end{description}

</lalLaTeX>

#endif /* autodoc block */


typedef enum
{
  ClosedInterval,     /* evaluate integral on a closed interval             */
  OpenInterval,       /* evaluate integral on an open interval              */
  SingularLowerLimit, /* integrate an inv sqrt singularity at lower limit   */
  SingularUpperLimit, /* integrate an inv sqrt singularity at upper limit   */
  InfiniteDomainPow,  /* integrate infinite domain with power-law falloff   */
  InfiniteDomainExp   /* integrate infinite domain with exponential falloff */
}
IntegralType;


typedef struct
tagSIntegrateIn
{
  void (*function)(LALStatus *s, REAL4 *y, REAL4 x, void *p);
  REAL4         xmax;
  REAL4         xmin;
  IntegralType  type;
}
SIntegrateIn;


typedef struct
tagDIntegrateIn
{
  void (*function)(LALStatus *s, REAL8 *y, REAL8 x, void *p);
  REAL8         xmax;
  REAL8         xmin;
  IntegralType  type;
}
DIntegrateIn;

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{IntegrateC}
</lalLaTeX>

#endif /* autodoc block */

void
LALSRombergIntegrate (
    LALStatus       *status,
    REAL4        *result,
    SIntegrateIn *input,
    void         *params
    );


void
LALDRombergIntegrate (
    LALStatus       *status,
    REAL8        *result,
    DIntegrateIn *input,
    void         *params
    );

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{IntegrateTestC}
</lalLaTeX>

#endif /* autodoc block */

#ifdef __cplusplus
}
#endif

#endif /* _INTEGRATE_H */
