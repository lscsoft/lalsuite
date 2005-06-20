#if 0 /* autodoc block */

<lalVerbatim file="FindRootHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{FindRoot.h}}
\label{s:FindRoot.h}

Root finding routines.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/FindRoot.h>
\end{verbatim}

\noindent This header covers the routines for root finding.

</lalLaTeX>

#endif /* autodoc block */

#ifndef _FINDROOT_H
#define _FINDROOT_H

#include <lal/LALDatatypes.h>

#ifdef __cplusplus
extern "C" {
#endif


NRCSID (FINDROOTH, "$Id$");

#if 0 /* autodoc block */

<lalLaTeX>
\subsection*{Error conditions}
\input{FindRootHErrTab}
</lalLaTeX>

<lalErrTable file="FindRootHErrTab">

#endif /* autodoc block */

#define FINDROOTH_ENULL 1
#define FINDROOTH_EIDOM 2
#define FINDROOTH_EMXIT 4
#define FINDROOTH_EBRKT 8

#define FINDROOTH_MSGENULL "Null pointer"
#define FINDROOTH_MSGEIDOM "Invalid initial domain"
#define FINDROOTH_MSGEMXIT "Maximum iterations exceeded"
#define FINDROOTH_MSGEBRKT "Root not bracketed"

#if 0 /* autodoc block */

</lalErrTable>

<lalLaTeX>

\subsection*{Structures}

\begin{verbatim}
typedef void (REAL4LALFunction) (LALStatus *s, REAL4 *y, REAL4 x, void *p);
typedef void (REAL8LALFunction) (LALStatus *s, REAL8 *y, REAL8 x, void *p);
\end{verbatim}

These are function pointers to functions that map real numbers to real numbers.

\begin{verbatim}
typedef struct
tagSFindRootIn
{
  void (*function)(LALStatus *s, REAL4 *y, REAL4 x, void *p);
  REAL4 xmax;
  REAL4 xmin;
  REAL4 xacc;
}
SFindRootIn;

typedef struct
tagDFindRootIn
{
  void (*function)(LALStatus *s, REAL8 *y, REAL8 x, void *p);
  REAL8 xmax;
  REAL8 xmin;
  REAL8 xacc;
}
DFindRootIn;
\end{verbatim}

These are the input structures to the root finding routines.  The fields are:

\begin{description}
\item[\texttt{function}] The function to find the root of.
\item[\texttt{xmax}] The maximum value of the domain interval to look for the root.
\item[\texttt{xmax}] The minimum value of the domain interval to look for the root.
\item[\texttt{xacc}] The accuracy desired for the root.
\end{description}

</lalLaTeX>

#endif /* autodoc block */


typedef struct
tagSFindRootIn
{
  void (*function)(LALStatus *s, REAL4 *y, REAL4 x, void *p);
  REAL4 xmax;
  REAL4 xmin;
  REAL4 xacc;
}
SFindRootIn;

typedef struct
tagDFindRootIn
{
  void (*function)(LALStatus *s, REAL8 *y, REAL8 x, void *p);
  REAL8 xmax;
  REAL8 xmin;
  REAL8 xacc;
}
DFindRootIn;

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{FindRootC}
</lalLaTeX>

#endif /* autodoc block */

void
LALSBracketRoot (
    LALStatus      *status,
    SFindRootIn *inout,
    void        *params
    );

int
XLALDBracketRoot(
    REAL8 (*y)(REAL8, void *),
    REAL8 *xmin,
    REAL8 *xmax,
    void *params
);

void
LALDBracketRoot (
    LALStatus      *status,
    DFindRootIn *inout,
    void        *params
    );

void
LALSBisectionFindRoot (
    LALStatus      *status,
    REAL4       *root,
    SFindRootIn *input,
    void        *params
    );

REAL8
XLALDBisectionFindRoot (
    REAL8 (*y)(REAL8, void *),
    REAL8 xmin,
    REAL8 xmax,
    REAL8 xacc,
    void *params
);

void
LALDBisectionFindRoot (
    LALStatus      *status,
    REAL8       *root,
    DFindRootIn *input,
    void        *params
    );

#if 0 /* autodoc block */

<lalLaTeX>
\newpage\input{FindRootTestC}
</lalLaTeX>

#endif /* autodoc block */

#ifdef __cplusplus
}
#endif

#endif /* _FINDROOT_H */
