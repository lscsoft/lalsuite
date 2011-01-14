/***************************************** <lalVerbatim file="GridCV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\subsection{Module \texttt{Grid.c}}
\label{ss:Grid.c}

Creates or destroys a LAL grid structure.

\subsubsection*{Prototypes}
\vspace{0.1in}
\begin{verbatim}
void
LAL<typecode>CreateGrid( LALStatus      *stat,
                         <datatype>Grid **grid,
                         UINT4Vector    *dimLength,
                         UINT4          dimension )

void
LAL<typecode>DestroyGrid( LALStatus      *stat,
                          <datatype>Grid **grid )
\end{verbatim}

\idx{LALI2CreateGrid()}
\idx{LALI4CreateGrid()}
\idx{LALI8CreateGrid()}
\idx{LALU2CreateGrid()}
\idx{LALU4CreateGrid()}
\idx{LALU8CreateGrid()}
\idx{LALSCreateGrid()}
\idx{LALDCreateGrid()}
\idx{LALCCreateGrid()}
\idx{LALZCreateGrid()}
\idx{LALI2DestroyGrid()}
\idx{LALI4DestroyGrid()}
\idx{LALI8DestroyGrid()}
\idx{LALU2DestroyGrid()}
\idx{LALU4DestroyGrid()}
\idx{LALU8DestroyGrid()}
\idx{LALSDestroyGrid()}
\idx{LALDDestroyGrid()}
\idx{LALCDestroyGrid()}
\idx{LALZDestroyGrid()}

\subsubsection*{Description}

These routines create or destroy a \verb@<datatype>Grid@ structure.
The input vector \verb@dimLength@ stores the lengths of each dimension
of the grid \emph{and} of the array at each grid point: in the
notation defined in \verb@Grid.h@, \verb@dimLength->length@$=M$.  The
parameter \verb@dimension@ gives the dimension $m$ of the physical
grid space; if $M>m$, then the remaining dimensions refer to a tangent
space at each grid point.  When creating a grid, the routines allocate
space for all the internal vectors and arrays, but no data are filled
in, with the exception of the \verb@(*grid)->data->dimLength@ vector
(which will contain exactly the same data as the \verb@dimLength@
input parameter).  When calling the \verb@LAL<typecode>CreateGrid()@
routines, or on returning from the \verb@LAL<typecode>DestroyGrid()@
routines, \verb@grid@ should be a non-\verb@NULL@ handle to a
\verb@NULL@-valued pointer.

For each of these prototype templates there are in fact 10 separate
routines corresponding to all the numerical atomic datatypes
\verb@<datatype>@ referred to by \verb@<typecode>@:
\begin{center}
\begin{tabular}{|c@{\qquad}c|c@{\qquad}c|}
\hline
\tt <typecode> & \tt <datatype> & \tt <typecode> & \tt <datatype> \\
\hline
\tt I2 & \tt  INT2 & \tt U2 & \tt    UINT2  \\
\tt I4 & \tt  INT4 & \tt U4 & \tt    UINT4  \\
\tt I8 & \tt  INT8 & \tt U8 & \tt    UINT8  \\
\tt  S & \tt REAL4 & \tt  C & \tt COMPLEX8  \\
\tt  D & \tt REAL8 & \tt  Z & \tt COMPLEX16 \\
\hline
\end{tabular}
\end{center}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALMalloc()                     LALFree()
LALDCreateVector()              LALDDestroyVector()
LAL<typecode>CreateArray()      LAL<typecode>DestroyArray()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{GridCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Grid.h>

NRCSID( GRIDC, "$Id$" );

#define TYPECODE Z
#define TYPE COMPLEX16
#include "Grid_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE C
#define TYPE COMPLEX8
#include "Grid_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE D
#define TYPE REAL8
#include "Grid_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "Grid_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I2
#define TYPE INT2
#include "Grid_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I4
#define TYPE INT4
#include "Grid_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I8
#define TYPE INT8
#include "Grid_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U2
#define TYPE UINT2
#include "Grid_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U4
#define TYPE UINT4
#include "Grid_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U8
#define TYPE UINT8
#include "Grid_source.c"
#undef TYPECODE
#undef TYPE
