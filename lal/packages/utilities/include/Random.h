#if 0 /* autodoc block */

<lalVerbatim file="RandomHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{Random.h}}
\label{s:Random.h}

Generates random numbers.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Random.h>
\end{verbatim}

\noindent This header covers the routines for generating random numbers.

</lalLaTeX>

#endif /* autodoc block */


#ifndef _RANDOM_H
#define _RANDOM_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (RANDOMH, "$Id$");

#if 0 /* autodoc block */

<lalLaTeX>
\subsection*{Error conditions}
\input{RandomHErrTab}
</lalLaTeX>

<lalErrTable file="RandomHErrTab">

#endif /* autodoc block */

#define RANDOMH_ENULL 1
#define RANDOMH_ENNUL 2
#define RANDOMH_ESIZE 4

#define RANDOMH_MSGENULL "Null pointer"
#define RANDOMH_MSGENNUL "Non-null pointer"
#define RANDOMH_MSGESIZE "Invalid size"

#if 0 /* autodoc block */

</lalErrTable>

<lalLaTeX>

\subsection*{Structures}

\begin{verbatim}
typedef struct tagRandomParams RandomParams;
\end{verbatim}

This structure contains the parameters necessary for generating the current
sequence of random numbers (based on the initial seed).  The contents should
not be manually adjusted.

</lalLaTeX>

#endif /* autodoc block */

typedef struct
tagRandomParams
{
  INT4 i;
  INT4 y;
  INT4 v[32];
}
RandomParams;

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{RandomC}
</lalLaTeX>

#endif /* autodoc block */

void
LALCreateRandomParams (
    LALStatus        *status,
    RandomParams **params,
    INT4           seed
    );

void
LALDestroyRandomParams (
    LALStatus        *status,
    RandomParams **params
    );

void
LALUniformDeviate (
    LALStatus       *status,
    REAL4        *deviate,
    RandomParams *params
    );

void
LALNormalDeviates (
    LALStatus       *status,
    REAL4Vector  *deviates,
    RandomParams *params
    );

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{RandomTestC}
</lalLaTeX>

#endif /* autodoc block */


#ifdef  __cplusplus
}
#endif

#endif /* _RANDOM_H */
