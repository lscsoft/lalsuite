/************************************ <lalVerbatim file="ZPGFilterHV">
Author: Creighton, T. D.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>

\section{Header \texttt{ZPGFilter.h}}
\label{s:ZPGFilter.h}

Provides routines to manipulate ZPG filters.

\subsection*{Synopsis}
\begin{verbatim}
#include "ZPGFilter.h"
\end{verbatim}

\noindent This header covers routines that create, destroy, and
transform objects of type \verb@<datatype>ZPGFilter@, where
\verb@<datatype>@ is either \verb@COMPLEX8@ or \verb@COMPLEX16@.
Generically, these data types can be used to store any rational
complex function in a factored form.  Normally this function is a
filter response, or ``transfer function'' $T(z)$, expressed in terms
of a complex frequency parameter $z=\exp(2\pi if\Delta t)$, where
$\Delta t$ is the sampling interval.  The rational function is
factored as follows:
$$
T(f) = g\times\frac{\prod_k (z-a_k)}{\prod_l (z-b_l)}
$$
where $g$ is the gain, $a_k$ are the (finite) zeros, and $b_l$ are the
(finite) poles.  It should be noted that rational functions always
have the same number of zeros as poles if one includes the point
$z=\infty$; any excess in the number of finite zeros or poles in the
rational expression simply indicates that there is a corresponding
pole or zero of that order at infinity.  It is also worth pointing out
that the ``gain'' is just the overall prefactor of this rational
function, and is not necessarily equal to the actual gain of the
transfer function at any particular frequency.

Another common complex frequency space is the $w$-space, obtained
from the $z$-space by the bilinear transformation:
$$
w = i\left(\frac{1-z}{1+z}\right) = \tan(\pi f\Delta t) , \quad
z = \frac{1+iw}{1-iw} \; .
$$
Other variables can also be used to represent the complex frequency
plane.  The \verb@<datatype>ZPGFilter@ structure can be used to
represent the transfer function in any of these spaces by transforming
the coordinates of the zeros and poles, and incorporating any residual
factors into the gain.  Care must be taken to include any zeros or
poles that are brought in from infinity by the transformation, and to
remove any zeros or poles which were sent to infinity.  Thus the
number of zeros and poles of the \verb@<datatype>ZPGFilter@ is not
necessarily constant under transformations!  Routines invoking the
\verb@<datatype>ZPGFilter@ data types should document which complex
variable is assumed.

</lalLaTeX> */

#ifndef _ZPGFILTER_H
#define _ZPGFILTER_H

#include "LALStdlib.h"

#ifdef __cplusplus
extern "C" {
#endif

NRCSID(ZPGFILTERH,"$Id$");

/* <lalLaTeX>

\subsection*{Error conditions}
\begin{tabular}{|c|l|l|}
\hline
status & status      & Explanation \\
 code  & description &             \\
\hline
\tt 1  & \tt Null pointer            & Missing a required pointer.           \\
\tt 2  & \tt Output already exists   & Can't allocate to a non-null pointer. \\
\tt 3  & \tt Memory allocation error & Could not allocate memory.            \\
\tt 4  & \tt Bad filter parameters   & Filter creation parameters outside of \\
       &                             & acceptable ranges.                    \\
\hline
\end{tabular}

</lalLaTeX> */

#define ZPGFILTER_ENUL 1
#define ZPGFILTER_EOUT 2
#define ZPGFILTER_EMEM 3
#define ZPGFILTER_EBAD 4

#define ZPGFILTER_MSGENUL "Null pointer"
#define ZPGFILTER_MSGEOUT "Output already exists"
#define ZPGFILTER_MSGEMEM "Memory allocation error"
#define ZPGFILTER_MSGEBAD "Bad filter parameters"

/* <lalLaTeX>
\subsection*{Structures}
</lalLaTeX> */

/* <lalLaTeX>
\vfill{\footnotesize\input{ZPGFilterHV}}
</lalLaTeX> */

/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{CreateZPGFilterC}
</lalLaTeX> */
void LALCreateCOMPLEX8ZPGFilter(LALStatus            *stat,
			     COMPLEX8ZPGFilter **output,
			     INT4              numZeros,
			     INT4              numPoles);

void LALCreateCOMPLEX16ZPGFilter(LALStatus             *stat,
			      COMPLEX16ZPGFilter **output,
			      INT4               numZeros,
			      INT4               numPoles);

/* <lalLaTeX>
\newpage\input{DestroyZPGFilterC}
</lalLaTeX> */
void LALDestroyCOMPLEX8ZPGFilter(LALStatus            *stat,
			      COMPLEX8ZPGFilter **input);

void LALDestroyCOMPLEX16ZPGFilter(LALStatus             *stat,
			       COMPLEX16ZPGFilter **input);

/* <lalLaTeX>
\newpage\input{BilinearTransformC}
</lalLaTeX> */
void LALWToZCOMPLEX8ZPGFilter(LALStatus            *stat,
			   COMPLEX8ZPGFilter *filter);

void LALWToZCOMPLEX16ZPGFilter(LALStatus             *stat,
			    COMPLEX16ZPGFilter *filter);

#ifdef __cplusplus
}
#endif

#endif /* _ZPGFILTER_H */
