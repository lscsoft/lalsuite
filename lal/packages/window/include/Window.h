/************************************************ <lalVerbatim file="WindowHV">
Authors: Allen, B. and Brown, D. A.
$Id$
**************************************************** </lalVerbatim> */
/*************************************************** <lalLaTeX>

\section{Header \texttt{Window.h}}
\label{s:Window.h}
\idx[Constant]{Rectangular}
\idx[Constant]{Hann}
\idx[Constant]{Welch}
\idx[Constant]{Bartlett}
\idx[Constant]{Parzen}
\idx[Constant]{Papoulis}
\idx[Constant]{Kaiser}
\index{Apodize}
\index{Taper}
\index{Power Spectrum}
\index{Bias in power spectrum}
\index{Spectral Estimation}

This header file contains enums that define the different types of windows,
and a parameter block which is used as input to the window-making function.
This allows you to crate and destroy a REAL4Window structure containing a
window (also called a taper, lag window, or apodization function).

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Window.h>
\end{verbatim}

\noindent
The routine prototyped in this header file creates REAL4Window structure
containing a window (also called a taper, lag window, or apodization
function) and the sum of squares of the window elements.  The choices
currently available are: 
\begin{itemize}
\item Rectangular
\item Hann
\item Welch
\item Bartlett
\item Parzen
\item Papoulis
\item Hamming
\item Kaiser
\end{itemize}
It should be straighforward to add additional window functions if
they are desired.

The window functions are defined for $j=0,\cdots,N-1$ by the following
formulae.  Note that $N$ is the length of the vector.  In these formulae, let
$x=2 \pi j/N$, and $y=|2j/N-1|$,
\begin{eqnarray*}
{\rm Rectangular:\ } w_j &=& 1 \\
{\rm Hann:\ } w_j &=& {1 \over 2} ( 1 - \cos  x  ) \\
{\rm Welch:\ } w_j &=& 1 -  y^2 \\
{\rm Bartlett:\ } w_j &=& 1 -  y \\
{\rm Parzen:\ } w_j &=&  1 - 6 y^2 + 6 y^3  {\rm\ if\ } y\le 1/2\\
                    &=&  2 (1-y)^3 {\rm\ if\ } y>1/2\\
{\rm Papoulis:\ } w_j &=& {1 \over \pi} \sin (\pi  y  ) + ( 1 -  y  ) \cos (\pi  y  )\\
{\rm Hamming:\ } w_j &=& 1-0.46 (1 + \cos x ) \\
{\rm Kaiser:\ } w_j &=& I_0\left( \beta\sqrt{1 - (y - 1)^2} \right)/I_0(\beta) \\
\end{eqnarray*}
where in the case of the Kaiser window $I_0(x)$ is the $0$th order, modified
Bessel function of the first kind, and $\beta (>=0)$ is a shape parameter
related to the amplitude of the sidelobes of the Fourier transform of the
window.  When $\beta=0$, a Kaiser window reduces to a rectangular window.

These window functions are shown in Fig.~\ref{f:window} for $N=1024$.

A couple of comments and warnings may be useful.  First, the definitions given
here are taken from {\it Numerical Recipes in C} \cite{numrec} and {\it
Spectral Analysis for Physical Applications} \cite{pw}.  The definitions of
windows are {\it not standard}.  In particular, some authors (e.g. J.G.
Proakis and D.G. Manolakis, {\it Digital Signal Processing} \cite{pm}) and
some standard computer applications (e.g. {\tt Matlab}) use definitions in
which $N$ is replaced by $N-1$ in the definitions of $x$ and $y$, with $j$
covering the range $j=0,\cdots,N-1$.  This has the advantage of making the
window function ``more symmetric'', typically be appending an extra ``0'' to
the end of the array.  It has the disadvantage that it throws out more of your
precious data.  The definition of the Kaiser window comes from ``Discrete-time
Signal Processing'' by Oppenheim and Schafer, p.474.

If you want to get a window function that agrees with either Proakis \&
Manolakis, or with Matlab, just call \verb|LALCreateREAL4Window()| with length
parameter $M-1$.  Then create an array of length $M$, copy the $M-1$ elements
of the window array returned by \verb|LALCreateREAL4Window()|, into it, and
finally copy the {\it first} element of the array returned by
\verb|LALCreateREAL4Window()| into the last element of your new array.

**********************************************************************</lalLaTeX> */

#ifndef _WINDOW_H
#define _WINDOW_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (WINDOWH, "$Id$");

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
*********************************************************** </lalLaTeX>*/
  /* <lalErrTable> */
#define WINDOWH_ENULLPARAM     1
#define WINDOWH_ENULLVECTOR    2
#define WINDOWH_EEALLOCATE     4
#define WINDOWH_EELENGTH       8
#define WINDOWH_ETYPEUNKNOWN  16
#define WINDOWH_ENULLHANDLE   32
#define WINDOWH_EWRONGLENGTH  64
#define WINDOWH_ENULLDATA    128
#define WINDOWH_ENNUL        256
#define WINDOWH_ENULL        512
#define WINDOWH_EBETA        1024

#define WINDOWH_MSGENULLPARAM    "null input parameter structure pointer"
#define WINDOWH_MSGENULLVECTOR   "null output vector pointer"
#define WINDOWH_MSGEEALLOCATE    "unable to allocate vector to store window"
#define WINDOWH_MSGEELENGTH      "length of window is <=0, must be positive"
#define WINDOWH_MSGETYPEUNKNOWN  "window is of unknown type"
#define WINDOWH_MSGENULLHANDLE   "input vector is null"
#define WINDOWH_MSGEWRONGLENGTH  "input vector is the wrong length"
#define WINDOWH_MSGENULLDATA     "data area of input vector is null"
#define WINDOWH_MSGENULL         "null pointer"
#define WINDOWH_MSGENNUL         "non-null pointer"
#define WINDOWH_MSGEBETA         "Invalid Kaiser window shape parameter"
  /*********************************************************** </lalErrTable>*/
/*<lalLaTeX> 

\subsection*{Types}
\subsubsection*{\texttt{enum WindowType}}
\idx[Type]{WindowType}
\idx[Constant]{WINDOWNAMELIST}

This enum defines the different possible types of windows that can be
generated.  Any code should take into account that this list may grow if
someone adds their favorite window to the list. {\bf  WARNING:} additional
window types must be added just before \texttt{NumberWindowTypes} and after
the existing window types.  Note that since an enum by default gives integers
starting at zero and incrementing by one, the enum \texttt{NumberWindowTypes}
will always give the correct number of window types, even if the list is
extended in the future.  The definition of the enum is:
\begin{verbatim}
typedef enum {Rectangular,
              Hann,
              Welch,
              Bartlett,
              Parzen,
              Papoulis,
              Hamming,
              Kaiser,
              NumberWindowTypes} WindowType;
\end{verbatim}
For convenience, the following macro is also defined
\begin{verbatim}
#define WINDOWNAMELIST {"Rectangular","Hann","Welch","Bartlett","Parzen","Papoulis",
                        "Hamming", "Kaiser"}
\end{verbatim}
This string can be used to print out the name of any of the windows (see the
test program for an example of this).  If a new window is added, be sure to
put its name onto the end of the array.


*****************************************************************</lalLaTeX> */
/* Define the types of available windows */
/* WARNING: additional window types must be added just before */
/* NumberWindowTypes, and after the existing window types */

typedef enum {Rectangular,
              Hann,
              Welch,
              Bartlett,
              Parzen,
              Papoulis,
              Hamming,
              Kaiser,
              /* add any new window types just before this comment */
              NumberWindowTypes} WindowType;

/* if you add new windows above, be sure to add a descriptive name below. */
#define WINDOWNAMELIST \
{"Rectangular","Hann","Welch","Bartlett","Parzen","Papoulis","Hamming","Kaiser"}

/*******************************************************************<lalLaTeX>

\subsubsection*{\texttt{Structure LALWindowParams}}
\idx[Type]{LALWindowParams}
This structure stores the parameters used to call the window function.
It is also used to return the sum of the vector squared.
The structure is defined by
\begin{verbatim}
typedef struct tagLALWindowParams {
  INT4        length;       <==> length of window (input)
  WindowType  type;         <==> type of window (input)
  REAL4       beta          <==> parameter used ony for Kaiser window (input)
  REAL8       sumofsquares; <==> sum of window squared  (output)
  CHAR*       windowname;   <==> pointer to a char string with window name (output)
} LALWindowParams;
\end{verbatim}
 The four fields are:
\begin{description}
\item[\texttt{INT4 length}] The length of the window. This is used as input.
\item[\texttt{WindowType type}] The type of the window. This is used as input.
\item[\texttt{REAL4 beta}] The shape parameter $\beta$ of a Kaiser window.
This is used as input.  
\item[\texttt{REAL8 sumofsquares}] The sum of the squares of the window
vector. This is used as output.
\item[\texttt{CHAR* windowname}] A pointer to a character string containing the window name. This is used as output.
\end{description}

****************************************************************</lalLaTeX> */

typedef struct tagLALWindowParams {
  INT4        length;       /* length of window (input) */
  WindowType  type;         /* type of window (input) */
  REAL4       beta;         /* shape parameters for Kaiser Window */
  REAL8       sumofsquares; /* sum of window squared  (output) */
  const CHAR* windowname;   /* pointer to a char string with window name (output) */
} LALWindowParams;

/*******************************************************************<lalLaTeX>

\subsubsection*{\texttt{Structure REAL4Window}}
\idx[Type]{REAL4Window}
This structure stores the window and it's parameters, such as the sum of the
vector squared.  The structure is defined by
\begin{verbatim}
typedef struct 
tagLALWindowParams 
{
  WindowType  type;                 <==> type of window
  REAL4Vector *data;                <==> window vector
  CHAR windowname[LALNameLength];   <==> char array with window name
  REAL4       beta;                 <==> shape parameter of Kaiser window
  REAL8       sumofsquares;         <==> sum of window squared
} LALWindowParams;
\end{verbatim}
 The four fields are:
\begin{description}
\item[\texttt{WindowType type}] The type of the window.
\item[\texttt{REAL4Vector *data}] The window data.
\item[\texttt{REAL4 beta}] The shape parameter $\beta$ of a Kaiser window.
\item[\texttt{REAL8 sumofsquares}] The sum of the squares of the window
vector.
\item[\texttt{CHAR windowname[LALNameLength]}] A pointer to a character string
containing the window name.  
\end{description}

****************************************************************</lalLaTeX> */
typedef struct
tagREAL4Window
{
  WindowType    type;
  REAL4Vector  *data;
  CHAR          windowname[LALNameLength];
  REAL4         beta;
  REAL8         sumofsquares;
}
REAL4Window;
  
void LALWindow (
    LALStatus *,
    REAL4Vector *, 
    LALWindowParams *
    );

void LALCreateREAL4Window (
    LALStatus          *status,
    REAL4Window       **output,
    LALWindowParams    *params
    );

void LALDestroyREAL4Window (
    LALStatus          *status,
    REAL4Window       **output
    );

/**************************************************************************** <lalLaTeX>
\vfill{\footnotesize\input{WindowHV}}
\newpage
\input{WindowC}
\newpage
\input{WindowTestC}
**************************************************************************** </lalLaTeX> */

#ifdef  __cplusplus
}
#endif

#endif /* _WINDOW_H */
