/**
\author Creighton, T. D.
\file

\heading{Module \ref Grid.c}
\latexonly\label{ss_Grid_c}\endlatexonly

Creates or destroys a LAL grid structure.

\heading{Prototypes}

\code
void
LAL<typecode>CreateGrid( LALStatus      *stat,
                         <datatype>Grid **grid,
                         UINT4Vector    *dimLength,
                         UINT4          dimension )

void
LAL<typecode>DestroyGrid( LALStatus      *stat,
                          <datatype>Grid **grid )
\endcode






















\heading{Description}

These routines create or destroy a <tt><datatype>Grid</tt> structure.
The input vector \c dimLength stores the lengths of each dimension
of the grid \e and of the array at each grid point: in the
notation defined in \ref Grid.h, <tt>dimLength->length</tt>\f$=M\f$.  The
parameter \c dimension gives the dimension \f$m\f$ of the physical
grid space; if \f$M>m\f$, then the remaining dimensions refer to a tangent
space at each grid point.  When creating a grid, the routines allocate
space for all the internal vectors and arrays, but no data are filled
in, with the exception of the <tt>(*grid)->data->dimLength</tt> vector
(which will contain exactly the same data as the \c dimLength
input parameter).  When calling the <tt>LAL<typecode>CreateGrid()</tt>
routines, or on returning from the <tt>LAL<typecode>DestroyGrid()</tt>
routines, \c grid should be a non-\c NULL handle to a
\c NULL-valued pointer.

For each of these prototype templates there are in fact 10 separate
routines corresponding to all the numerical atomic datatypes
<tt><datatype></tt> referred to by <tt><typecode></tt>:

<table><tr><td>

\tt <typecode></td><td>\tt <datatype></td><td>\tt <typecode></td><td>\tt <datatype></td></tr>
<tr><td>
\tt I2</td><td>\tt  INT2</td><td>\tt U2</td><td>\tt    UINT2</td></tr>
<tr><td>\tt I4</td><td>\tt  INT4</td><td>\tt U4</td><td>\tt    UINT4</td></tr>
<tr><td>\tt I8</td><td>\tt  INT8</td><td>\tt U8</td><td>\tt    UINT8</td></tr>
<tr><td>\tt  S</td><td>\tt REAL4</td><td>\tt  C</td><td>\tt COMPLEX8</td></tr>
<tr><td>\tt  D</td><td>\tt REAL8</td><td>\tt  Z</td><td>\tt COMPLEX16</td></tr>
</table>


\heading{Algorithm}

\heading{Uses}
\code
lalDebugLevel
LALMalloc()                     LALFree()
LALDCreateVector()              LALDDestroyVector()
LAL<typecode>CreateArray()      LAL<typecode>DestroyArray()
\endcode

\heading{Notes}



*/

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
