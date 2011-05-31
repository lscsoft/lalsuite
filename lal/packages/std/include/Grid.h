/*
*  Copyright (C) 2007 Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**
\author Creighton, T. D.
\file

\brief Provides a structured datatype for a multidimensional rectilinear grid.

\heading{Synopsis}
\code
#include <lal/Grid.h>
\endcode

This header provides a new structured datatype storing data
on a multidimensional rectilinear grid.  It is in some sense a
generalization of the series datatypes (frequency series, time series,
etc.), representing evenly-sampled data over some physical parameter
space.

*/

#ifndef _GRID_H
#define _GRID_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

NRCSID( GRIDH, "$Id$" );

/**
\heading{Error conditions}
 \name Error Codes */ /*@{*/
#define GRIDH_ENUL 1
#define GRIDH_EOUT 2
#define GRIDH_EMEM 3

#define GRIDH_MSGENUL "Unexpected null pointer in arguments"
#define GRIDH_MSGEOUT "Output handle points to a non-null pointer"
#define GRIDH_MSGEMEM "Memory allocation error"
/*@}*//**

\heading{Types}

\heading{Structure <tt>\<datatype\>Grid</tt>}

This structure is a generalization of the LAL series types,
storing data on an \f$m\f$-dimensional rectangular grid on a physical
parameter space.  The values on the grid are of type <tt>\<datatype\></tt>
which can be any LAL primitive \e numerical datatype (\c INT2,
\c INT4, \c INT8, \c UINT2, \c UINT4, \c UINT8,
\c REAL4, \c REAL8, \c COMPLEX8, \c COMPLEX16).  The
data are stored in an array of dimension \f$M\geq m\f$: if \f$M=m\f$, then the
structure stores a single value for each grid point; if \f$M>m\f$, then
the structure stores a vector or array of values on the ``tangent
space'' of each grid point.  We refer to \f$m\f$ as the \e grid
dimension and \f$M\f$ as the \e data dimension.  The fields of the
structure are:

<dl>
<dt><tt>CHAR name[LALNameLength]</tt></dt><dd> A name identifying the grid
and/or the data being sampled.</dd>

<dt><tt>LALUnit sampleUnits</tt></dt><dd> The physical units of the
quantities on the grid.</dd>

<dt><tt>LALUnit *dimUnits</tt></dt><dd> The physical units of the grid axes.
This must be allocated as an array of length \f$m\f$.</dd>

<dt><tt>REAL8Vector *offset</tt></dt><dd> A vector \f$\mathbf{p}_0\f$
specifying the location of the grid point indexed by \f$(0,\ldots,0)\f$.
Must have dimensionality \f$m\f$.</dd>

<dt><tt>REAL8Vector *interval</tt></dt><dd> The vector \f$\Delta\mathbf{p}\f$
specifying the grid spacing in each dimension.  Must have
dimensionality \f$m\f$.</dd>

<dt><tt>\<datatype\>Array *data</tt></dt><dd> Pointer to an array storing the
data values at the corresponding grid points.  The data dimension
\f$M=\f$<tt>data->dimLength->length</tt> must be greater than or equal to
the grid dimension \f$m=\f$<tt>offset->length</tt>=<tt>interval->length</tt>.
An index \f$\mathbf{i}=(i_0,\ldots,i_{M-1})\f$, where \f$i_k\f$ are integers
from 0 to <tt>data->dimLength->data</tt>\f${}_k\f$, specifies a grid point
located at \f$\mathbf{p}=\mathbf{p}_0+\sum_{k=0}^{n-1}\hat{\mathbf{e}}_k
i_k \Delta p_k\f$ if \f$M=m\f$, or an array element \f$\mathsf{A}_{i_m\cdots
i_{M-1}}\f$ at that grid point if \f$M>m\f$.  The values in
<tt>data->data</tt> are the value stored at each grid point (or array
element at each grid point), arranged in the manner discussed in
\ref LALDatatypes.h.</dd>
</dl>
*/

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






/* Function prototypes. */




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

/**
%
*/

#ifdef __cplusplus
}
#endif

#endif /* _GRID_H */
