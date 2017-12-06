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

#ifndef _GRID_H
#define _GRID_H

#include <lal/LALStdlib.h>

#ifdef __cplusplus
extern "C" {
#elif 0
}       /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup Grid_h
 * \author Creighton, T. D.
 *
 * \brief Provides a structured datatype for a multidimensional rectilinear grid.
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/Grid.h>
 * \endcode
 *
 * This header provides a new structured datatype storing data
 * on a multidimensional rectilinear grid.  It is in some sense a
 * generalization of the series datatypes (frequency series, time series,
 * etc.), representing evenly-sampled data over some physical parameter
 * space.
 *
 * ### Structure <tt>\<datatype\>Grid</tt> ###
 *
 * This structure is a generalization of the LAL series types,
 * storing data on an \f$m\f$-dimensional rectangular grid on a physical
 * parameter space.  The values on the grid are of type <tt>\<datatype\></tt>
 * which can be any LAL primitive \e numerical datatype (#INT2,
 * #INT4, #INT8, #UINT2, #UINT4, #UINT8,
 * #REAL4, #REAL8, #COMPLEX8, #COMPLEX16).  The
 * data are stored in an array of dimension \f$M\geq m\f$: if \f$M=m\f$, then the
 * structure stores a single value for each grid point; if \f$M>m\f$, then
 * the structure stores a vector or array of values on the ``tangent
 * space'' of each grid point.  We refer to \f$m\f$ as the \e grid
 * dimension and \f$M\f$ as the \e data dimension.  The fields of the
 * structure are:
 *
 * <dl>
 * <dt><tt>CHAR name[LALNameLength]</tt></dt><dd> A name identifying the grid
 * and/or the data being sampled.</dd>
 *
 * <dt><tt>LALUnit sampleUnits</tt></dt><dd> The physical units of the
 * quantities on the grid.</dd>
 *
 * <dt><tt>LALUnit *dimUnits</tt></dt><dd> The physical units of the grid axes.
 * This must be allocated as an array of length \f$m\f$.</dd>
 *
 * <dt><tt>REAL8Vector *offset</tt></dt><dd> A vector \f$\mathbf{p}_0\f$
 * specifying the location of the grid point indexed by \f$(0,\ldots,0)\f$.
 * Must have dimensionality \f$m\f$.</dd>
 *
 * <dt><tt>REAL8Vector *interval</tt></dt><dd> The vector \f$\Delta\mathbf{p}\f$
 * specifying the grid spacing in each dimension.  Must have
 * dimensionality \f$m\f$.</dd>
 *
 * <dt><tt>\<datatype\>Array *data</tt></dt><dd> Pointer to an array storing the
 * data values at the corresponding grid points.  The data dimension
 * \f$M=\f$<tt>data->dimLength->length</tt> must be greater than or equal to
 * the grid dimension \f$m=\f$<tt>offset->length</tt>=<tt>interval->length</tt>.
 * An index \f$\mathbf{i}=(i_0,\ldots,i_{M-1})\f$, where \f$i_k\f$ are integers
 * from 0 to <tt>data->dimLength->data</tt>\f${}_k\f$, specifies a grid point
 * located at \f$\mathbf{p}=\mathbf{p}_0+\sum_{k=0}^{n-1}\hat{\mathbf{e}}_k
 * i_k \Delta p_k\f$ if \f$M=m\f$, or an array element \f$\mathsf{A}_{i_m\cdots
 * i_{M-1}}\f$ at that grid point if \f$M>m\f$.  The values in
 * <tt>data->data</tt> are the value stored at each grid point (or array
 * element at each grid point), arranged in the manner discussed in
 * \ref LALDatatypes.h.</dd>
 * </dl>
 *
 * ### Prototypes ###
 *
 * \code
 * void
 * LAL<typecode>CreateGrid( LALStatus      *stat,
 * <datatype>Grid **grid,
 * UINT4Vector    *dimLength,
 * UINT4          dimension )
 *
 * void
 * LAL<typecode>DestroyGrid( LALStatus      *stat,
 * <datatype>Grid **grid )
 * \endcode
 *
 * ### Description ###
 *
 * These routines create or destroy a <tt>\<datatype\>Grid</tt> structure.
 * The input vector \c dimLength stores the lengths of each dimension
 * of the grid \e and of the array at each grid point: in the
 * notation defined in \ref Grid.h, <tt>dimLength->length</tt>\f$=M\f$.  The
 * parameter \c dimension gives the dimension \f$m\f$ of the physical
 * grid space; if \f$M>m\f$, then the remaining dimensions refer to a tangent
 * space at each grid point.  When creating a grid, the routines allocate
 * space for all the internal vectors and arrays, but no data are filled
 * in, with the exception of the <tt>(*grid)->data->dimLength</tt> vector
 * (which will contain exactly the same data as the \c dimLength
 * input parameter).  When calling the <tt>LAL\<typecode\>CreateGrid()</tt>
 * routines, or on returning from the <tt>LAL\<typecode\>DestroyGrid()</tt>
 * routines, \c grid should be a non-\c NULL handle to a
 * \c NULL-valued pointer.
 *
 * For each of these prototype templates there are in fact 10 separate
 * routines corresponding to all the numerical atomic datatypes
 * <tt>\<datatype\></tt> referred to by <tt>\<typecode\></tt>:
 *
 * <table>
 * <tr><th>\<typecode\></th><th>\<datatype\></th><th>\<typecode\></th><th>\<datatype\></th></tr>
 * <tr><td>I2</td><td>  INT2</td><td> U2</td><td>    UINT2</td></tr>
 * <tr><td>I4</td><td>  INT4</td><td> U4</td><td>    UINT4</td></tr>
 * <tr><td>I8</td><td>  INT8</td><td> U8</td><td>    UINT8</td></tr>
 * <tr><td> S</td><td> REAL4</td><td>  C</td><td> COMPLEX8</td></tr>
 * <tr><td> D</td><td> REAL8</td><td>  Z</td><td> COMPLEX16</td></tr>
 * </table>
 *
 */
/*@{*/

/** \name Error Codes */ /*@{ */
#define GRIDH_ENUL 1    /**< Unexpected null pointer in arguments */
#define GRIDH_EOUT 2    /**< Output handle points to a non-null pointer */
#define GRIDH_EMEM 3    /**< Memory allocation error */
/*@}*/
/*@}*/

/** \cond DONT_DOXYGEN */
#define GRIDH_MSGENUL "Unexpected null pointer in arguments"
#define GRIDH_MSGEOUT "Output handle points to a non-null pointer"
#define GRIDH_MSGEMEM "Memory allocation error"
/** \endcond */

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
LALI2CreateGrid(LALStatus * status, INT2Grid ** grid,
                UINT4Vector * dimLength, UINT4 dimension);

void
LALI4CreateGrid(LALStatus * status, INT4Grid ** grid,
                UINT4Vector * dimLength, UINT4 dimension);

void
LALI8CreateGrid(LALStatus * status, INT8Grid ** grid,
                UINT4Vector * dimLength, UINT4 dimension);

void
LALU2CreateGrid(LALStatus * status, UINT2Grid ** grid,
                UINT4Vector * dimLength, UINT4 dimension);

void
LALU4CreateGrid(LALStatus * status, UINT4Grid ** grid,
                UINT4Vector * dimLength, UINT4 dimension);

void
LALU8CreateGrid(LALStatus * status, UINT8Grid ** grid,
                UINT4Vector * dimLength, UINT4 dimension);

void
LALSCreateGrid(LALStatus * status, REAL4Grid ** grid,
               UINT4Vector * dimLength, UINT4 dimension);

void
LALDCreateGrid(LALStatus * status, REAL8Grid ** grid,
               UINT4Vector * dimLength, UINT4 dimension);

void
LALCCreateGrid(LALStatus * status, COMPLEX8Grid ** grid,
               UINT4Vector * dimLength, UINT4 dimension);

void
LALZCreateGrid(LALStatus * status, COMPLEX16Grid ** grid,
               UINT4Vector * dimLength, UINT4 dimension);


void LALI2DestroyGrid(LALStatus * status, INT2Grid ** grid);

void LALI4DestroyGrid(LALStatus * status, INT4Grid ** grid);

void LALI8DestroyGrid(LALStatus * status, INT8Grid ** grid);

void LALU2DestroyGrid(LALStatus * status, UINT2Grid ** grid);

void LALU4DestroyGrid(LALStatus * status, UINT4Grid ** grid);

void LALU8DestroyGrid(LALStatus * status, UINT8Grid ** grid);

void LALSDestroyGrid(LALStatus * status, REAL4Grid ** grid);

void LALDDestroyGrid(LALStatus * status, REAL8Grid ** grid);

void LALCDestroyGrid(LALStatus * status, COMPLEX8Grid ** grid);

void LALZDestroyGrid(LALStatus * status, COMPLEX16Grid ** grid);

#if 0
{       /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _GRID_H */
