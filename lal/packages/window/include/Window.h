/************************************************ <lalVerbatim file="WindowHV">
Authors: Allen, B., Brown, D. A., and Creighton, T.
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
\idx[Constant]{Creighton}
\index{Apodize}
\index{Taper}
\index{Power Spectrum}
\index{Bias in power spectrum}
\index{Spectral Estimation}

This header file provides routines and structures to create and store
window functions (also called a taper, lag window, or apodization
function).

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Window.h>
\end{verbatim}

\noindent
The routines prototyped in this header file create and destroy
\verb@REAL4Window@ and \verb@REAL8Window@ structures, each containing
a window and ancillary information such as the the sum of squares of
the window elements.  It is conventional to express windows as
functions on the normalized domain $y\in[-1,1]$, where the window is
implicitly zero outside this domain.  The available windows and their
formulae are:
\begin{center}\begin{tabular}{l@{\qquad$w(y)\;=\;$}l}
Rectanglar: & $1$ \\[1ex]
Hann:       & $\cos^2(\pi y/2)$ \\[1ex]
Welch:      & $1-y^2$ \\[1ex]
Bartlett:   & $1-|y|$ \\[1ex]
Parzen:     & $\left\{\begin{array}{c@{\qquad}c}
		1-6y^2(1-|y|) & |y|\leq1/2 \\
		2(1-|y|)^3    & |y|>1/2 \end{array}\right.$ \\[3ex]
Papoulis:   & $\frac{1}{\pi}\sin(\pi |y|)+(1-|y|)\cos(\pi |y|)$ \\[1ex]
Hamming:    & $1-0.46[1-\cos(\pi y)]$ \\[1ex]
Kaiser:     & $I_0\left(\beta\sqrt{1-(1-|y|)^2}\right)
		/I_0\left(\beta\right)$ \\[1ex]
Creighton:  & $\exp\left[-\beta y^2/(1-y^2)\right]$
\end{tabular}\end{center}

\begin{wrapfigure}{r}{0.65\textwidth}
\vspace{-4ex}
\begin{center}
\resizebox{0.6\textwidth}{!}{\includegraphics{window_t}} \\
\parbox{0.6\textwidth}{\caption{\label{f:window-t} Various windows as
functions of the normalized independend variable $y$, choosing
$\beta=6$ for Kaiser and $\beta=2$ for Creighton.}}
\end{center}
\vspace{-2ex}
\end{wrapfigure}
\noindent where $I_0(x)$ is the $0$th order, modified Bessel function
of the first kind, and $\beta(\geq0)$ is a continuous shape parameter
related to the amplitude of the sidelobes of the Fourier transform of
the window.  As $\beta\rightarrow0$, the Kaiser and Creighton windows
approach a rectangular window.

It is straightforward to add new window functions as the need arises.

For discretely-sampled data $w_j$, $j=0,\ldots,N-1$, the
transformation to the normalized domain $y$ is \emph{not
standardized}.  In these routines we adopt the convention that
$y=2j/N-1$, which is the one used by \textit{Numerical Recipes in
C}~\cite{numrec} and \textit{Spectral Analysis for Physical
Applications}~\cite{pw}.  Other authors (e.g.\ J.G.  Proakis and
D.G. Manolakis, \textit{Digital Signal Processing}~\cite{pm}) and some
standard computer applications (e.g.\ \texttt{Matlab}) adopt the
convention that $y=2j/(N-1)-1$.  This has the advantage of making the
window function exactly symmetric about $N/2$, typically by setting
the final array element to zero.  It has the disadvantage that it
throws out one more of your precious data.

If you want to get a window function that agrees with either Proakis \&
Manolakis, or with Matlab, just call \verb|LALCreateREAL4Window()| with length
parameter $M-1$.  Then create an array of length $M$, copy the $M-1$ elements
of the window array returned by \verb|LALCreateREAL4Window()|, into it, and
finally copy the {\it first} element of the array returned by
\verb|LALCreateREAL4Window()| into the last element of your new array.

\begin{wrapfigure}{l}{0.65\textwidth}
\vspace{-2ex}
\begin{center}
\resizebox{0.6\textwidth}{!}{\includegraphics{window_f}} \\
\parbox{0.6\textwidth}{\caption{\label{f:window-f} Frequency behaviour
of various windows as functions of the inverse of the normalized
independend variable $y$, choosing $\beta=6$ for Kaiser and $\beta=2$
for Creighton.  Solid lines demark the central lobe, circles mark the
peaks of the sidelobes.}}
\end{center}
\vspace{-2ex}
\end{wrapfigure}
These window functions are shown as functions of $y$ in
Fig.~\ref{f:window-t}.  Their Fourier transforms are shown as
functions of $1/y$ in Fig.~\ref{f:window-f}.  Since the Fourier
transform of windowed data is the Fourier transform of the data
convolved with the Fourier transform of the window,
Fig.~\ref{f:window-f} is the major guideline for selecting a window.
One can see that windows with a narrow central lobe tend to have
higher sidelobes, and windows which suppress their low-order sidelobes
tend to have more power in the high-order sidelobes.  The choice of
window thus depends on whether one is trying to resolve nearby
spectral features of comparable magnitude (suggesting a rectangular or
a Welch window), to reduce spectral bias and low-order sidelobes (a
Hamming or Kaiser window), or to measure a broad spectrum with a large
dynamical range (a Creighton or a Papoulis window).  A second
consideration in some circumstances is computational cost: a
rectangular window is trivial, and polynomial windows (e.g.\ Welch,
Bartlett, Parzen) are somewhat cheaper than windows requiring a
function call (e.g.\ Hann, Papoulis, Kaiser, Creighton).

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
#define WINDOWH_EBETA       1024

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
#define WINDOWH_MSGEBETA         "Invalid window shape parameter"
  /*********************************************************** </lalErrTable>*/
/*<lalLaTeX> 

\subsection*{Types}
\subsubsection*{Enumeration \texttt{WindowType}}
\idx[Type]{WindowType}
\idx[Constant]{WINDOWNAMELIST}

This enumeration defines the different possible types of windows that can be
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
              Creighton,
              NumberWindowTypes} WindowType;
\end{verbatim}
For convenience, the following macro is also defined
\begin{verbatim}
#define WINDOWNAMELIST {"Rectangular","Hann","Welch","Bartlett","Parzen","Papoulis",
                        "Hamming", "Kaiser", "Creighton"}
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
	      Creighton,
              /* add any new window types just before this comment */
              NumberWindowTypes} WindowType;

/* if you add new windows above, be sure to add a descriptive name below. */
#define WINDOWNAMELIST \
{"Rectangular","Hann","Welch","Bartlett","Parzen","Papoulis","Hamming","Kaiser","Creighton"}

/*******************************************************************<lalLaTeX>

\subsubsection*{Structure \texttt{LALWindowParams}}
\idx[Type]{LALWindowParams}
This structure stores the parameters used to call the window function.
It is also used to return the sum of the vector squared and the window
name.  The fields of the structure are:
\begin{description}
\item[\texttt{INT4 length}] The length of the window. This is used as input.
\item[\texttt{WindowType type}] The type of the window. This is used as input.
\item[\texttt{REAL4 beta}] The shape parameter $\beta$ of certain windows.
This is used as input.  
\item[\texttt{REAL8 sumofsquares}] The sum of the squares of the window
vector. This is used as output.
\item[\texttt{CHAR* windowname}] A pointer to a character string containing the window name. This is used as output.
\end{description}

****************************************************************</lalLaTeX> */

typedef struct tagLALWindowParams {
  INT4        length;       /* length of window (input) */
  WindowType  type;         /* type of window (input) */
  REAL4       beta;         /* shape parameters for certain windows */
  REAL8       sumofsquares; /* sum of window squared  (output) */
  const CHAR* windowname;   /* pointer to a char string with window name (output) */
} LALWindowParams;

/*******************************************************************<lalLaTeX>

\subsubsection*{Structure \texttt{<datatype>Window}}
\idx[Type]{REAL4Window}
This structure stores the window and it's parameters (above), where
\verb@<datatype>@ can be either \verb@REAL4@ or \verb@REAL8@.  The
fields of the structure are:
\begin{description}
\item[\texttt{WindowType type}] The type of the window.
\item[\texttt{<datatype>Vector *data}] The window data.
\item[\texttt{REAL4 beta}] The shape parameter $\beta$ of the window.
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

typedef struct
tagREAL8Window
{
  WindowType    type;
  REAL8Vector  *data;
  CHAR          windowname[LALNameLength];
  REAL4         beta;
  REAL8         sumofsquares;
}
REAL8Window;


REAL4Window *XLALCreateREAL4Window( UINT4 length, WindowType type, REAL4 beta );
REAL8Window *XLALCreateREAL8Window( UINT4 length, WindowType type, REAL4 beta );
REAL4Window *XLALCreateRectangularREAL4Window( UINT4 length );
REAL4Window *XLALCreateHannREAL4Window( UINT4 length );
REAL4Window *XLALCreateWelchREAL4Window( UINT4 length );
REAL4Window *XLALCreateBartlettREAL4Window( UINT4 length );
REAL4Window *XLALCreateParzenREAL4Window( UINT4 length );
REAL4Window *XLALCreatePapoulisREAL4Window( UINT4 length );
REAL4Window *XLALCreateHammingREAL4Window( UINT4 length );
REAL4Window *XLALCreateKaiserREAL4Window( UINT4 length, REAL4 beta );
REAL4Window *XLALCreateCreightonREAL4Window( UINT4 length, REAL4 beta );
REAL8Window *XLALCreateRectangularREAL8Window( UINT4 length );
REAL8Window *XLALCreateHannREAL8Window( UINT4 length );
REAL8Window *XLALCreateWelchREAL8Window( UINT4 length );
REAL8Window *XLALCreateBartlettREAL8Window( UINT4 length );
REAL8Window *XLALCreateParzenREAL8Window( UINT4 length );
REAL8Window *XLALCreatePapoulisREAL8Window( UINT4 length );
REAL8Window *XLALCreateHammingREAL8Window( UINT4 length );
REAL8Window *XLALCreateKaiserREAL8Window( UINT4 length, REAL4 beta );
REAL8Window *XLALCreateCreightonREAL8Window( UINT4 length, REAL4 beta );
void XLALDestroyREAL4Window( REAL4Window *window );
void XLALDestroyREAL8Window( REAL8Window *window );



  
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

void LALCreateREAL8Window (
    LALStatus          *status,
    REAL8Window       **output,
    LALWindowParams    *params
    );

void LALDestroyREAL8Window (
    LALStatus          *status,
    REAL8Window       **output
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
