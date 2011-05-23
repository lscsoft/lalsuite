/**
\author Tibbits, M. M.
\file

\heading{Module \ref MatrixMultiply.c}
\latexonly\label{s_MatrixMultiply_c}\endlatexonly

This file is dedicated to reproducing the matlab function " .* ".  This file
has several declarations of the same function taking all forms of available
input.  This being said, I have yet to script the complex actions and their
counterparts.

\heading{Description}

This file is to help make the conversion from Matlab to c much earier.
In this file, we have created all of the versions of .* that we plan on
using.

\heading{Algorithms}

The algorithm is the same as it is in matlab.  The dot in front of an operator
in matlab signifies that if either or both of the operands are vectors, then
the operation will be carried out member by member.  For instance

vector a[25];
vector b[25];
vector c[25];

c = a .* b;

The result of this is:

c[0] =	a[0] *	b[0];
c[1] =	a[1] *	b[1];
.	.	.
.	.	.
.	.	.

etc.

\heading{Uses}

<ul>
<li> \c LALDCreateVector</li>
</ul>

\heading{Notes}

At the current time none of the operations have been specified for neither the
complex datatypes nor the unsigned datatypes.

\heading{Prototypes}

*/

#include <math.h>
#include "Matrix.h"

NRCSID( MATLABMATRIXMULTC, "$Id$");

#define TYPECODE D
#define TYPE REAL8
#define TYPECODE2 D
#define TYPE2 REAL8
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE D
#define TYPE REAL8
#define TYPECODE2 S
#define TYPE2 REAL4
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE D
#define TYPE REAL8
#define TYPECODE2 I2
#define TYPE2 INT2
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE D
#define TYPE REAL8
#define TYPECODE2 I4
#define TYPE2 INT4
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE D
#define TYPE REAL8
#define TYPECODE2 I8
#define TYPE2 INT8
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE S
#define TYPE REAL4
#define TYPECODE2 S
#define TYPE2 REAL4
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE S
#define TYPE REAL4
#define TYPECODE2 I2
#define TYPE2 INT2
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE S
#define TYPE REAL4
#define TYPECODE2 I4
#define TYPE2 INT4
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE S
#define TYPE REAL4
#define TYPECODE2 I8
#define TYPE2 INT8
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I2
#define TYPE INT2
#define TYPECODE2 I2
#define TYPE2 INT2
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I4
#define TYPE INT4
#define TYPECODE2 I2
#define TYPE2 INT2
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I4
#define TYPE INT4
#define TYPECODE2 I4
#define TYPE2 INT4
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I8
#define TYPE INT8
#define TYPECODE2 I2
#define TYPE2 INT2
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I8
#define TYPE INT8
#define TYPECODE2 I4
#define TYPE2 INT4
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2

#define TYPECODE I8
#define TYPE INT8
#define TYPECODE2 I8
#define TYPE2 INT8
#include "MatrixMultiply_source.c"
#undef TYPECODE
#undef TYPE
#undef TYPECODE2
#undef TYPE2






