/***************************************** <lalVerbatim file="GridHV">
Author: Creighton, T. D.
$Id$
**************************************************** </lalVerbatim> */

/********************************************************** <lalLaTeX>

\section{Header \texttt{Grid.h}}
\label{s:Grid.h}

Provides a structured datatype for a multidimensional rectilinear
grid.

\subsection*{Synopsis}
\begin{verbatim}
#include <lal/Grid.h>
\end{verbatim}

\noindent This header provides a new structured datatype storing data
on a multidimensional rectilinear grid.  It is in some sense a
generalization of the series datatypes (frequency series, time series,
etc.), representing evenly-sampled data over some physical parameter
space.

******************************************************* </lalLaTeX> */

#ifndef _GRID_H
#define _GRID_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID( GRIDH, "$Id$" );

/********************************************************** <lalLaTeX>
\subsection*{Error conditions}
****************************************** </lalLaTeX><lalErrTable> */
#define GRIDH_ENUL 1
#define GRIDH_EOUT 2
#define GRIDH_EMEM 3

#define GRIDH_MSGENUL "Unexpected null pointer in arguments"
#define GRIDH_MSGEOUT "Output handle points to a non-null pointer"
#define GRIDH_MSGEMEM "Memory allocation error"
/******************************************** </lalErrTable><lalLaTeX>

\subsection*{Types}

\subsubsection*{Structure \texttt{<datatype>Grid}}
\idx[Type]{INT2Grid}
\idx[Type]{INT4Grid}
\idx[Type]{INT8Grid}
\idx[Type]{UINT2Grid}
\idx[Type]{UINT4Grid}
\idx[Type]{UINT8Grid}
\idx[Type]{REAL4Grid}
\idx[Type]{REAL8Grid}
\idx[Type]{COMPLEX8Grid}
\idx[Type]{COMPLEX16Grid}

\noindent This structure is a generalization of the LAL series types,
storing data on an $m$-dimensional rectangular grid on a physical
parameter space.  The values on the grid are of type \verb@<datatype>@
which can be any LAL primitive \emph{numerical} datatype (\verb@INT2@,
\verb@INT4@, \verb@INT8@, \verb@UINT2@, \verb@UINT4@, \verb@UINT8@,
\verb@REAL4@, \verb@REAL8@, \verb@COMPLEX8@, \verb@COMPLEX16@).  The
data are stored in an array of dimension $M\geq m$: if $M=m$, then the
structure stores a single value for each grid point; if $M>m$, then
the structure stores a vector or array of values on the ``tangent
space'' of each grid point.  We refer to $m$ as the \emph{grid}
dimension and $M$ as the \emph{data} dimension.  The fields of the
structure are:

\begin{description}
\item[\texttt{CHAR name[LALNameLength]}] A name identifying the grid
and/or the data being sampled.

\item[\texttt{LALUnit sampleUnits}] The physical units of the
quantities on the grid.

\item[\texttt{LALUnit *dimUnits}] The physical units of the grid axes.
This must be allocated as an array of length $m$.

\item[\texttt{REAL8Vector *offset}] A vector $\mathbf{p}_0$
specifying the location of the grid point indexed by $(0,\ldots,0)$.
Must have dimensionality $m$.

\item[\texttt{REAL8Vector *interval}] The vector $\Delta\mathbf{p}$
specifying the grid spacing in each dimension.  Must have
dimensionality $m$.

\item[\texttt{<datatype>Array *data}] Pointer to an array storing the
data values at the corresponding grid points.  The data dimension
$M=$\verb@data->dimLength->length@ must be greater than or equal to
the grid dimension $m=$\verb@offset->length@=\verb@interval->length@.
An index $\mathbf{i}=(i_0,\ldots,i_{M-1})$, where $i_k$ are integers
from 0 to \verb@data->dimLength->data@${}_k$, specifies a grid point
located at $\mathbf{p}=\mathbf{p}_0+\sum_{k=0}^{n-1}\hat{\mathbf{e}}_k
i_k \Delta p_k$ if $M=m$, or an array element $\mathsf{A}_{i_m\cdots
i_{M-1}}$ at that grid point if $M>m$.  The values in
\verb@data->data@ are the value stored at each grid point (or array
element at each grid point), arranged in the manner discussed in
\verb@LALDatatypes.h@.
\end{description}
******************************************************* </lalLaTeX> */

typedef struct tagINT2Grid {
  CHAR name[LALNameLength];
  LALUnit sampleUnits;
  LALUnit *dimUnits;
  REAL8Vector *offset;
  REAL8Vector *interval;
  INT2Array *data;
} INT2Grid;

typedef struct tagINT4Grid {
  CHAR name[LALNameLength];
  LALUnit sampleUnits;
  LALUnit *dimUnits;
  REAL8Vector *offset;
  REAL8Vector *interval;
  INT4Array *data;
} INT4Grid;

typedef struct tagINT8Grid {
  CHAR name[LALNameLength];
  LALUnit sampleUnits;
  LALUnit *dimUnits;
  REAL8Vector *offset;
  REAL8Vector *interval;
  INT8Array *data;
} INT8Grid;

typedef struct tagUINT2Grid {
  CHAR name[LALNameLength];
  LALUnit sampleUnits;
  LALUnit *dimUnits;
  REAL8Vector *offset;
  REAL8Vector *interval;
  UINT2Array *data;
} UINT2Grid;

typedef struct tagUINT4Grid {
  CHAR name[LALNameLength];
  LALUnit sampleUnits;
  LALUnit *dimUnits;
  REAL8Vector *offset;
  REAL8Vector *interval;
  UINT4Array *data;
} UINT4Grid;

typedef struct tagUINT8Grid {
  CHAR name[LALNameLength];
  LALUnit sampleUnits;
  LALUnit *dimUnits;
  REAL8Vector *offset;
  REAL8Vector *interval;
  UINT8Array *data;
} UINT8Grid;

typedef struct tagREAL4Grid {
  CHAR name[LALNameLength];
  LALUnit sampleUnits;
  LALUnit *dimUnits;
  REAL8Vector *offset;
  REAL8Vector *interval;
  REAL4Array *data;
} REAL4Grid;

typedef struct tagREAL8Grid {
  CHAR name[LALNameLength];
  LALUnit sampleUnits;
  LALUnit *dimUnits;
  REAL8Vector *offset;
  REAL8Vector *interval;
  REAL8Array *data;
} REAL8Grid;

typedef struct tagCOMPLEX8Grid {
  CHAR name[LALNameLength];
  LALUnit sampleUnits;
  LALUnit *dimUnits;
  REAL8Vector *offset;
  REAL8Vector *interval;
  COMPLEX8Array *data;
} COMPLEX8Grid;

typedef struct tagCOMPLEX16Grid {
  CHAR name[LALNameLength];
  LALUnit sampleUnits;
  LALUnit *dimUnits;
  REAL8Vector *offset;
  REAL8Vector *interval;
  COMPLEX16Array *data;
} COMPLEX16Grid;


/* <lalLaTeX>
\vfill{\footnotesize\input{GridHV}}
</lalLaTeX> */

/* Function prototypes. */

/* <lalLaTeX>
\newpage\input{GridC}
</lalLaTeX> */
void
LALI2CreateGrid( LALStatus *status, INT2Grid **grid, UINT4Vector *dimLength, UINT4 dimension );

void
LALI4CreateGrid( LALStatus *status, INT4Grid **grid, UINT4Vector *dimLength, UINT4 dimension );

void
LALI8CreateGrid( LALStatus *status, INT8Grid **grid, UINT4Vector *dimLength, UINT4 dimension );

void
LALU2CreateGrid( LALStatus *status, UINT2Grid **grid, UINT4Vector *dimLength, UINT4 dimension );

void
LALU4CreateGrid( LALStatus *status, UINT4Grid **grid, UINT4Vector *dimLength, UINT4 dimension );

void
LALU8CreateGrid( LALStatus *status, UINT8Grid **grid, UINT4Vector *dimLength, UINT4 dimension );

void
LALSCreateGrid( LALStatus *status, REAL4Grid **grid, UINT4Vector *dimLength, UINT4 dimension );

void
LALDCreateGrid( LALStatus *status, REAL8Grid **grid, UINT4Vector *dimLength, UINT4 dimension );

void
LALCCreateGrid( LALStatus *status, COMPLEX8Grid **grid, UINT4Vector *dimLength, UINT4 dimension );

void
LALZCreateGrid( LALStatus *status, COMPLEX16Grid **grid, UINT4Vector *dimLength, UINT4 dimension );


void
LALI2DestroyGrid( LALStatus *status, INT2Grid **grid );

void
LALI4DestroyGrid( LALStatus *status, INT4Grid **grid );

void
LALI8DestroyGrid( LALStatus *status, INT8Grid **grid );

void
LALU2DestroyGrid( LALStatus *status, UINT2Grid **grid );

void
LALU4DestroyGrid( LALStatus *status, UINT4Grid **grid );

void
LALU8DestroyGrid( LALStatus *status, UINT8Grid **grid );

void
LALSDestroyGrid( LALStatus *status, REAL4Grid **grid );

void
LALDDestroyGrid( LALStatus *status, REAL8Grid **grid );

void
LALCDestroyGrid( LALStatus *status, COMPLEX8Grid **grid );

void
LALZDestroyGrid( LALStatus *status, COMPLEX16Grid **grid );

/* <lalLaTeX>
%\newpage\input{GridTestC}
</lalLaTeX> */

#ifdef __cplusplus
}
#endif

#endif /* _GRID_H */
