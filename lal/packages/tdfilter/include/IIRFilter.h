/************************************ <lalVerbatim file="IIRFilterHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{IIRFilter.h}}
\label{s:IIRFilter.h}

Provides routines to make and apply IIR filters.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/IIRFilter.h>
\end{verbatim}

\noindent This header covers routines that create, destroy, and apply
generic time-domain filters, given by objects of type
\verb@<datatype>IIRFilter@, where \verb@<datatype>@ is either
\verb@REAL4@ or \verb@REAL8@.

An IIR (Infinite Impulse Response) filter is a generalized linear
causal time-domain filter, in which the filter output $y_n=y(t_n)$ at
any sampled time $t_n=t_0+n\Delta t$ is a linear combination of the
input $x$ \emph{and} output $y$ at previous sampled times:
$$
y_n = \sum_{k=0}^M c_k x_{n-k} + \sum_{l=1}^N d_l y_{n-l} \; .
$$
The coefficients $c_k$ are called the direct filter coefficients, and
the coefficients $d_l$ are the recursive filter coefficients.  The
filter order is the larger of $M$ or $N$, and determines how far back
in time the filter must look to determine its next output.  However,
the recursive nature of the filter means that the output can depend on
input arbitrarily far in the past; hence the name ``infinite impulse
response''.  Nonetheless, for a well-designed, stable filter, the
actual filter response to an impulse should diminish rapidly beyond
some characteristic timescale.

Note that nonrecursive FIR (Finite Impulse Response) filters are
considered a subset of IIR filters, having $N=0$.

For practical implementation, it is convenient to express the bilinear
equation above as two linear equations involving an auxiliary sequence
$w$:
$$
w_n = x_n + \sum_{l=1}^N d_l w_{n-l} \; ,
$$
$$
y_n = \sum_{k=0}^M c_k w_{n-k} \; .
$$
The equivalence of this to the first expression is not obvious, but
can be proven by mathematical induction.  The advantage of the
auxiliary variable representation is twofold.  First, when one is
feeding data point by point to the filter, the filter needs only
``remember'' the previous $M$ or $N$ (whichever is larger) values of
$w$, rather than remembering the previous $M$ values of $x$ \emph{and}
the previous $N$ values of $y$.  Second, when filtering a large stored
data vector, the filter response can be computed in place: one first
runs forward through the vector replacing $x$ with $w$, and then
backward replacing $w$ with $y$.

Although the IIR filters in these routines are explicitly real, one
can consider formally their complex response.  A sinusoidal input can
thus be written as $x_n=X\exp(2\pi ifn\Delta t)=Xz^n$, where $X$ is a
complex amplitude and $z=\exp(2\pi if\Delta t)$ is a complex
parametrization of the frequency.  By linearity, the output must also
be sinusoidal: $y_m=Y\exp(2\pi ifm\Delta t)=Yz^m$.  Putting these into
the bilinear equation, one can easily compute the filter's complex
transfer function:
$$
T(z) = \frac{Y}{X} = \frac{\sum_{k=0}^M c_k z^{-k}}
                      {1 - \sum_{l=1}^N d_l z^{-l}}
$$
This can be readily converted to and from the ``zeros, poles, gain''
representation of a filter, which expresses $T(z)$ as a factored
rational function of $z$.

It should also be noted that, in the routines covered by this header,
I have adopted the convention of including a redundant recursive
coefficient $d_0$, in order to make the indexing more intuitive.  For
formal correctness $d_0$ should be set to $-1$, although the filtering
routines never actually use this coefficient.

******************************************************* </lalLaTeX> */

#ifndef _IIRFILTER_H
#define _IIRFILTER_H

#include <lal/LALStdlib.h>
#include <lal/ZPGFilter.h>

#ifdef __cplusplus
extern "C" {
#pragma }
#endif

NRCSID(IIRFILTERH,"$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define IIRFILTERH_ENUL  1
#define IIRFILTERH_EOUT  2
#define IIRFILTERH_EMEM  3
#define IIRFILTERH_EPAIR 4

#define IIRFILTERH_MSGENUL  "Unexpected null pointer in arguments"
#define IIRFILTERH_MSGEOUT  "Output handle points to a non-null pointer"
#define IIRFILTERH_MSGEMEM  "Memory allocation error"
#define IIRFILTERH_MSGEPAIR "Input has unpaired nonreal poles or zeros"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{<datatype>IIRFilter}}
\idx[Type]{REAL4IIRFilter}
\idx[Type]{REAL8IIRFilter}

This structure stores the direct and recursive filter coefficients, as
well as the history of the auxiliary sequence $w$.  \verb@<datatype>@
may be \verb@REAL4@ or \verb@REAL8@.  The length of the history vector
gives the order of the filter.  The fields are:

\begin{description}
\item[\texttt{const CHAR *name}] A user-assigned name.

\item[\texttt{REAL8 deltaT}] The sampling time interval of the filter.
  If $\leq0$, it will be ignored (i.e.\ it will be taken from the data
  stream).

\item[\texttt{<datatype>Vector *directCoef}] The direct filter
  coefficients.

\item[\texttt{<datatype>Vector *recursCoef}] The recursive filter
  coefficients.

\item[\texttt{<datatype>Vector *history}] The previous values of $w$.
\end{description}

******************************************************* </lalLaTeX> */

typedef struct tagREAL4IIRFilter{
  const CHAR *name;              /* User assigned name. */
  REAL8 deltaT;            /* Sampling time interval. */
  REAL4Vector *directCoef; /* The direct filter coefficients. */
  REAL4Vector *recursCoef; /* The recursive filter coefficients. */
  REAL4Vector *history;    /* The previous values of w. */
} REAL4IIRFilter;

typedef struct tagREAL8IIRFilter{
  const CHAR *name;              /* User assigned name. */
  REAL8 deltaT;            /* Sampling time interval. */
  REAL8Vector *directCoef; /* The direct filter coefficients. */
  REAL8Vector *recursCoef; /* The recursive filter coefficients. */
  REAL8Vector *history;    /* The previous values of w. */
} REAL8IIRFilter;


/* <lalLaTeX>
\vfill{\footnotesize\input{IIRFilterHV}}
</lalLaTeX> */


/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{CreateIIRFilterC}
</lalLaTeX> */
void
LALCreateREAL4IIRFilter( LALStatus         *stat,
			 REAL4IIRFilter    **output,
			 COMPLEX8ZPGFilter *input );

void
LALCreateREAL8IIRFilter( LALStatus          *stat,
			 REAL8IIRFilter     **output,
			 COMPLEX16ZPGFilter *input );

/* <lalLaTeX>
\newpage\input{DestroyIIRFilterC}
</lalLaTeX> */
void
LALDestroyREAL4IIRFilter( LALStatus      *stat,
			  REAL4IIRFilter **input );

void
LALDestroyREAL8IIRFilter( LALStatus      *stat,
			  REAL8IIRFilter **input );

/* <lalLaTeX>
\newpage\input{IIRFilterC}
</lalLaTeX> */
void
LALIIRFilterREAL4( LALStatus      *stat,
		   REAL4          *output,
		   REAL4          input,
		   REAL4IIRFilter *filter );

void
LALIIRFilterREAL8( LALStatus      *stat,
		   REAL8          *output,
		   REAL8          input,
		   REAL8IIRFilter *filter );

REAL4
LALSIIRFilter( REAL4 x, REAL4IIRFilter *filter );

REAL8
LALDIIRFilter( REAL8 x, REAL8IIRFilter *filter );

/* <lalLaTeX>
\newpage\input{IIRFilterVectorC}
</lalLaTeX> */
void
LALIIRFilterREAL4Vector( LALStatus      *stat,
			 REAL4Vector    *vector,
			 REAL4IIRFilter *filter );

void
LALIIRFilterREAL8Vector( LALStatus      *stat,
			 REAL8Vector    *vector,
			 REAL8IIRFilter *filter );

void
LALDIIRFilterREAL4Vector( LALStatus      *stat,
			  REAL4Vector    *vector,
			  REAL8IIRFilter *filter );

/* <lalLaTeX>
\newpage\input{IIRFilterVectorRC}
</lalLaTeX> */
void
LALIIRFilterREAL4VectorR( LALStatus      *stat,
			  REAL4Vector    *vector,
			  REAL4IIRFilter *filter );

void
LALIIRFilterREAL8VectorR( LALStatus      *stat,
			  REAL8Vector    *vector,
			  REAL8IIRFilter *filter );

void
LALDIIRFilterREAL4VectorR( LALStatus      *stat,
			   REAL4Vector    *vector,
			   REAL8IIRFilter *filter );

/* <lalLaTeX>
\newpage\input{IIRFilterTestC}
</lalLaTeX> */

#ifdef __cplusplus
#pragma {
}
#endif

#endif /* _IIRFILTER_H */
