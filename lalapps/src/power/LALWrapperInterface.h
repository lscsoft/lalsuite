#if 0  /* autodoc block */

<lalVerbatim file="LALWrapperInterfaceHV">
$Id$
</lalVerbatim>

<lalLaTeX>

\section{Header \texttt{LALWrapperInterface.h}}
\label{s:LALWrapperInterface.h}

Provides the generic LAL function interface with \texttt{wrapperAPI}.

\subsection*{Synopsis}
\begin{verbatim}
#include <LALWrapperInterface.h>
\end{verbatim}

\noindent This header covers the generic LAL routines that are ultimately
called by \texttt{wrapperAPI}.  These routines have the same form for all
searches, though their implementation is search specific.

\subsection*{Structures}

\begin{verbatim}
struct LALSearchInput
struct LALSearchOutput
struct LALInitSearchParams
struct LALMPIParams
struct LALApplySearchParams
\end{verbatim}

\subsection*{Prototypes}
\input{LALWrapperInterfaceHP}

\noindent These are the generic LAL routines that form the LAL-Wrapper
interface.

\begin{description}
\item[\texttt{LALInitSearch}] \relax
\item[\texttt{LALConditionData}] \relax
\item[\texttt{LALApplySearch}] \relax
\item[\texttt{LALFinalizeSearch}] \relax
\end{description}

</lalLaTeX>

#endif /* autodoc block */

#ifndef _LALWRAPPERINTERFACE_H
#define _LALWRAPPERINTERFACE_H

#include <wrapperInterfaceDatatypes.h>
#include <lal/LALDatatypes.h>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID( LALWRAPPERINTERFACEH, "$Id$" );


typedef inPut        LALSearchInput;
typedef SearchOutput LALSearchOutput;
typedef InitParams   LALInitSearchParams;
typedef SearchParams LALMPIParams;

typedef struct
tagLALApplySearchParams
{
  LALMPIParams *mpiParams;
  void         *searchParams;
}
LALApplySearchParams;


/* <lalVerbatim file="LALWrapperInterfaceHP"> */

void
LALInitSearch(
    LALStatus             *status,
    void                 **searchParams,
    LALInitSearchParams   *initSearchParams
    );


void
LALConditionData(
    LALStatus             *status,
    LALSearchInput        *inout,
    void                  *searchParams,
    LALMPIParams          *mpiParams
    );


void
LALApplySearch(
    LALStatus             *status,
    LALSearchOutput       *output,
    LALSearchInput        *input,
    LALApplySearchParams  *params
    );


void
LALFinalizeSearch(
    LALStatus             *status,
    void                 **searchParams
    );

/* </lalVerbatim> */

#ifdef __cplusplus
}
#endif

#if 0  /* autodoc block */

<lalLaTeX>
\vfill{\footnotesize\input{LALWrapperInterfaceHV}}
\newpage\input{LALWrapperInterfaceC}
\newpage\input{happyAPIC}
</lalLaTeX>

#endif /* autodoc block */

#endif /* _LALWRAPPERINTERFACE_H */
